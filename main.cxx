#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPoints2D.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData2DFS.h>
#include <vtkLine.h>
#include <vtkDataSetMapper.h>
#include <vtkTetra.h>
#include <vtkProperty.h>


#include <string>
#include <iostream>  
#include <vector> 

#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkMetaImageWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPlane.h>
#include <vtkKdTree.h>
#include <vtkPointLocator.h>
#include <vtkPolyLine.h>
#include <vtkPointLocator.h>
#include <vtkMatrix4x4.h>

#include "main.h"
 

using namespace std;


// Define interaction style
class MouseInteractorStyle: public vtkInteractorStyleTrackballCamera
{
public:

	static MouseInteractorStyle* New();
	vtkTypeMacro(MouseInteractorStyle, vtkInteractorStyleTrackballCamera);

	MouseInteractorStyle() 
	{
		SamplePoly = vtkSmartPointer<vtkPolyData>::New();
		isControlPoint = vtkSmartPointer<vtkIdTypeArray>::New();
		ControlPointCoord = vtkSmartPointer<vtkDoubleArray>::New();
		RelaxShape = vtkSmartPointer<vtkPoints>::New();

		LeftButtonDown = false;
	}

	~MouseInteractorStyle()	{}

	bool InitialPointcloud()
	{
		//if (BoundaryPolynormalGenerator == NULL)
		//{
		//	std::cerr << "cannot find BoundaryPolynormalGenerator" << std::endl;
		//	return false;
		//}
		//vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPolynormalGenerator->GetOutput()->GetPointData()->GetArray("Normals"));
		if (BoundaryNormalArray == NULL)
		{
			std::cerr << "cannot find BoundaryNormalArray" << std::endl;
			return false;
		}

		connectionCellArray->Reset();
	//	relationships.clear();
	//	parameters.clear();

		vtkSmartPointer<vtkPoints> SamplePoints = SamplePoly->GetPoints();
		for (int i = 0; i < SamplePoints->GetNumberOfPoints();)
		{
			double coordi[3] = { 0.0, 0.0, 0.0 };
			coordi[0] = vtkMath::Random(-3.0 * CRADIUS, 3.0 * CRADIUS);
			coordi[1] = vtkMath::Random(-3.0 * CRADIUS, 3.0 * CRADIUS);
			coordi[2] = vtkMath::Random(-3.0 * CRADIUS, 3.0 * CRADIUS);

			vtkIdType nearestboundaryPID = boundarypointLocator->FindClosestPoint(coordi);
			double boundarycoord[3];
			BoundaryPoly->GetPoint(nearestboundaryPID, boundarycoord);
			double dir_p2boundary[3];
			vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
			vtkMath::Normalize(dir_p2boundary);
			double boundarynormal[3];
			BoundaryNormalArray->GetTuple(nearestboundaryPID, boundarynormal);
			if (vtkMath::Dot(dir_p2boundary, boundarynormal) < 0)
				continue;
			SamplePoints->SetPoint(i, coordi);
			i ++;
		}

		for (int i = 0; i < R->GetNumberOfTuples(); i ++)
		{
			R->SetValue(i, R_INITIAL);
		}

		return true;
	}

	bool UniformRedistribution(int iterationnumber)
	{
		//if (BoundaryPolynormalGenerator == NULL)
		//{
		//	std::cerr << "cannot find BoundaryPolynormalGenerator" << std::endl;
		//	return false;
		//}
		//vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPolynormalGenerator->GetOutput()->GetPointData()->GetArray("Normals"));
		if (BoundaryNormalArray == NULL)
		{
			std::cerr << "cannot find BoundaryNormalArray" << std::endl;
			return false;
		}
		
		for (int iter = 0; iter < iterationnumber; iter++)
		{
			for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				double coordi[3];
				SamplePoly->GetPoints()->GetPoint(i, coordi);
				double Ri = R->GetValue(i);

				double Fi[3] = { 0.0, 0.0, 0.0 };
				double Fi_c[3] = { 0.0, 0.0, 0.0 };
				double Fi_p[3] = { 0.0, 0.0, 0.0 };
								
				// force between center and this point
				{
					//vtkIdType nearestboundaryPID = boundarypointLocator->FindClosestPoint(coordi);
					vtkSmartPointer<vtkIdList> NeighorboundaryIds = vtkSmartPointer<vtkIdList>::New();
					//pointLocator->FindPointsWithinRadius(4.0, coordi, NeighorpIds);
					boundarypointLocator->FindClosestNPoints(5, coordi, NeighorboundaryIds);

					double boundarycoord[3];
					BoundaryPoly->GetPoint(NeighorboundaryIds->GetId(0), boundarycoord);
					double dir_p2boundary[3];
					vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
					vtkMath::Normalize(dir_p2boundary);
					double boundarynormal[3] = {0.0, 0.0, 0.0};
					for (int idxj = 0; idxj < NeighorboundaryIds->GetNumberOfIds(); idxj++)
					{
						vtkIdType j = NeighorboundaryIds->GetId(idxj);
						double boundarynormalj[3];
						BoundaryNormalArray->GetTuple(j, boundarynormalj);
						vtkMath::Add(boundarynormal, boundarynormalj, boundarynormal);
					}
					vtkMath::Normalize(boundarynormal);

					double F_mode_c = 0.0;
					if (vtkMath::Dot(dir_p2boundary, boundarynormal) > 0)
						F_mode_c = 0.0;
					else
						F_mode_c = -1000.0;

					for (int l = 0; l < 3; l++) Fi_c[l] += F_mode_c * boundarynormal[l];
				}

				// force between points
				double Fpi_abssum = 0.0;
				double Fpi_sumabs = 0.0;

				vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New();
				//pointLocator->FindPointsWithinRadius(4.0, coordi, NeighorpIds);
				pointLocator->FindClosestNPoints(10, coordi, NeighorpIds);

				for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj++)
				{
					vtkIdType j = NeighorpIds->GetId(idxj);
					if (i == j) continue;

					double coordj[3];
					SamplePoly->GetPoint(j, coordj);

					double dir[3];
					vtkMath::Subtract(coordi, coordj, dir);
					double dis = vtkMath::Norm(dir);
					vtkMath::Normalize(dir);

					double Rj = R->GetValue(j);
					double R_total = Ri + Rj;

					if (dis < 0.55 * R_total)
					{
						double temp = 0.5 / R_total * dis + 0.5;
						double Fij_p_mode = (1.0 / (temp*temp*temp*temp*temp*temp) - 1.0 / (temp));
						double Fij_p[3] = { 0.0, 0.0, 0.0 };
						for (int l = 0; l < 3; l++) Fij_p[l] = Fij_p_mode * dir[l];
						for (int l = 0; l < 3; l++) Fi_p[l] += Fij_p[l];
						Fpi_sumabs += vtkMath::Norm(Fij_p);

						//	connections[i].push_back(j);
					}
				}

				Fpi_abssum = vtkMath::Norm(Fi_p);

				// combine Fi_cl and Fi_p
				for (int l = 0; l < 3; l++) Fi[l] = Fi_c[l] + 3.0 * Fi_p[l];

				// move points based on force
				if (vtkMath::Norm(Fi) > 20.0)
				{
					vtkMath::Normalize(Fi);
					for (int l = 0; l < 3; l++) Fi[l] = 20.0 * Fi[l];
				}
				
				if (isControlPoint->GetValue(i) != 2  // if it is not a tumor center point
					&& isControlPoint->GetValue(i) != 3) // if it is not a Jbar 
				{
					double coordi_new[3];
					for (int l = 0; l < 3; l++) coordi_new[l] = coordi[l] + Fi[l] / POINTMASS * TIMESTEP;
					SamplePoly->GetPoints()->SetPoint(i, coordi_new);
				}

				F_abssum->SetValue(i, Fpi_abssum);
				F_sumabs->SetValue(i, Fpi_sumabs);
			}

			// update radius Ri
			//double meanR = CRADIUS * CRADIUS / SamplePoints->GetNumberOfPoints();
			for (int i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				double Ri = R->GetValue(i);
				double delta_Ri_Max = 3.0 * TIMESTEP;

				double delta_Ri = delta_Ri_Max;
				double CondenseForce = F_sumabs->GetValue(i) - F_abssum->GetValue(i) - 5.5 * BorderForce; // 6 means hexagon grid mesh
				//	delta_Ri = (delta_Ri_Max + 1.0 - 0.1 * exp(CondenseForce / Hardness)) * timestep;	
				double temp = 1.0 + 0.003 * CondenseForce;
				if (temp <= 0.05)
					delta_Ri = delta_Ri_Max;
				else
					delta_Ri = -3.0 * log(temp) * TIMESTEP;
				if (delta_Ri > delta_Ri_Max)
					delta_Ri = delta_Ri_Max;

				Ri += delta_Ri;
				Ri = Ri < 0.05 ? 0.05 : Ri;
				Ri = Ri > 10.0 ? 10.0 : Ri;

				R->SetValue(i, Ri);
			}
		}

		return true;
	}

	bool BuildConnections()
	{
		this->Connections.clear();

		// draw the connections
		connectionCellArray->Reset();
		isTumorConnections->Reset();

		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
		{
			double coordi[3];
			SamplePoly->GetPoints()->GetPoint(i, coordi);
			double Ri = R->GetValue(i);

			vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New();
			//pointLocator->FindPointsWithinRadius(4.0, coordi, NeighorpIds);
			pointLocator->FindClosestNPoints(10, coordi, NeighorpIds);

			vector<vtkIdType> Connection_i;

			for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj++)
			{
				vtkIdType j = NeighorpIds->GetId(idxj);
				if (i == j) continue;

				double coordj[3];
				SamplePoly->GetPoints()->GetPoint(j, coordj);

				double dir[3];
				vtkMath::Subtract(coordi, coordj, dir);
				double dis = vtkMath::Norm(dir);

				double Rj = R->GetValue(j);
				double R_total = Ri + Rj;

				if (dis < 0.90 * R_total)
					Connection_i.push_back(j);
			}

			this->Connections.push_back(Connection_i);
		}

		// make sure tumor and Jbar have connection
		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
		{
			if (isControlPoint->GetValue(i) == 2)
			{
				for (vtkIdType j = i + 1; j < SamplePoly->GetPoints()->GetNumberOfPoints(); j++)
				{
					if (isControlPoint->GetValue(j) == 3)
						this->Connections.at(i).push_back(j);
				}
			}
		}

		// refine connections
		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			for (unsigned int idj = 0; idj < this->Connections[i].size(); idj++)
		{
			vtkIdType j = this->Connections[i].at(idj);
			bool flag_iinjfile = false;
			for (unsigned int idk = 0; idk < this->Connections[j].size(); idk++)
			{
				if (i == this->Connections[j].at(idk))
				{
					flag_iinjfile = true;
					break;
				}
			}
			if (flag_iinjfile == false)
				this->Connections[j].push_back(i);
		}

		// build connectionCellArray for display
		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			for (unsigned int idj = 0; idj < this->Connections[i].size(); idj++)
			{
				vtkIdType j = this->Connections[i].at(idj);
				if (j <= i)	continue;
				vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
				line->GetPointIds()->SetId(0, i);
				line->GetPointIds()->SetId(1, j);
				connectionCellArray->InsertNextCell(line);
				//if (isControlPoint->GetValue(i) == 2 && isControlPoint->GetValue(j) == 3
				//	|| isControlPoint->GetValue(i) == 3 && isControlPoint->GetValue(j) == 2)
				//{
				//	isTumorConnections->InsertNextValue(2);
				//	continue;
				//}

				if (isControlPoint->GetValue(i) == 2 || isControlPoint->GetValue(j) == 2)
					isTumorConnections->InsertNextValue(1);
				else
					isTumorConnections->InsertNextValue(0);

			}

		return true;
	}
	
	bool FindSurfaceControlPoints()
	{
		RelatedBoundaryPids.clear();

		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
		{
			if (isControlPoint->GetValue(i) == 2   // if i is a tumor control point
				|| isControlPoint->GetValue(i) == 3) // if i is a J-bar control point
			{				
			//	double tempcoord[3] = { Jbarx, Jbary, Jbarz };
			//	isControlPoint->SetValue(i, 2);
			//	ControlPointCoord->SetTuple(i, tempcoord);
				RelatedBoundaryPids.push_back(-1);
				continue;
			}

			double coordi[3];
			SamplePoly->GetPoint(i, coordi);

			vtkIdType nearestboundaryPID = boundarypointLocator->FindClosestPoint(coordi);
			double boundarycoord[3];
			BoundaryPoly->GetPoint(nearestboundaryPID, boundarycoord);
			double dir_p2boundary[3];
			vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
			double dis_p2boundary = vtkMath::Norm(dir_p2boundary);
			
			if (dis_p2boundary < 2e-1)
			{
				isControlPoint->SetValue(i, 1);
				ControlPointCoord->SetTuple(i, boundarycoord);
				RelatedBoundaryPids.push_back(nearestboundaryPID);
			}
			else
			{
				double tempcoord[3] = { 0.0, 0.0, 0.0 };
				isControlPoint->SetValue(i, 0);
				ControlPointCoord->SetTuple(i, tempcoord);
				RelatedBoundaryPids.push_back(-1);
			}
		}

		return true;
	}

	bool DeformationMotion(int controltag1, int controltag2) // set controltag2 = -1 if not using this parameter
	{
		if (this->Connections.size() == 0)
			return false;
		if (isControlPoint == NULL)
			return false;
		if (RelaxShape == NULL)
			return false;

		//if (BoundaryPolynormalGenerator == NULL)
		//{
		//	std::cerr << "cannot find BoundaryPolynormalGenerator" << std::endl;
		//	return false;
		//}
		//vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPolynormalGenerator->GetOutput()->GetPointData()->GetArray("Normals"));
		if (BoundaryNormalArray == NULL)
		{
			std::cerr << "cannot find BoundaryNormalArray" << std::endl;
			return false;
		}
		
		{
			vtkSmartPointer<vtkPoints> SamplePoints = SamplePoly->GetPoints();
			std::memset(&forces.at(0), 0, forces.size() * sizeof(CForce));

			for (int i = 0; i < RelaxShape->GetNumberOfPoints(); i++)
			{
				double RelaxCenter[3];
				RelaxShape->GetPoint(i, RelaxCenter);
				double CurrentCenter[3];
				SamplePoints->GetPoint(i, CurrentCenter);

				double** RelaxShape_local = new double*[Connections[i].size()];
				double** CurrentShape_local = new double*[Connections[i].size()];
				for (unsigned int j = 0; j < Connections[i].size(); j ++)
				{
					double RelaxCoord[3], CurrentCoord[3];

					RelaxShape->GetPoint(Connections[i][j], RelaxCoord);
					SamplePoints->GetPoint(Connections[i][j], CurrentCoord);
					
					vtkMath::Subtract(RelaxCoord, RelaxCenter, RelaxCoord);
					vtkMath::Subtract(CurrentCoord, CurrentCenter, CurrentCoord);
					
					RelaxShape_local[j] = new double[3];
					CurrentShape_local[j] = new double[3];

					for (int l = 0; l < 3; l ++)
					{
						RelaxShape_local[j][l] = RelaxCoord[l];
						CurrentShape_local[j][l] = CurrentCoord[l];
					}
				}



				// add a syn rotation
				//double** rot = new double*[3];
				//for (int j = 0; j < 3; j++)
				//	rot[j] = new double[3];
				//rot[0][0] = 0.8754; rot[0][1] = 0.2346; rot[0][2] = -0.4226;
				//rot[1][0] = -0.2588; rot[1][1] = 0.9659; rot[1][2] = 0;
				//rot[2][0] = 0.4082; rot[2][1] = 0.1094; rot[2][2] = 0.9063;
				//vtkMath::MultiplyMatrix(RelaxShape_local, rot, Connections[i].size(), 3, 3, 3, CurrentShape_local);
				//
				//for (int j = 0; j < 3; j++)
				//	delete[] rot[j];
				//delete[] rot;

				// calculate rotation
				double** RelaxShape_local_T = new double*[3];
				for (int j = 0; j < 3; j++)
					RelaxShape_local_T[j] = new double[Connections[i].size()];

				for (int j = 0; j < 3; j++)
					for (unsigned int k = 0; k < Connections[i].size(); k++)
					{
						RelaxShape_local_T[j][k] = RelaxShape_local[k][j];
					}

				double** C = new double*[3];
				for (int j = 0; j < 3; j++)
				{
					C[j] = new double[3];
				}

				vtkMath::MultiplyMatrix(RelaxShape_local_T, CurrentShape_local, 3, (unsigned int)Connections[i].size(), (unsigned int)Connections[i].size(), 3, C);
				double A[3][3];
				memcpy(&A[0][0], &C[0][0], 3 * sizeof(double));
				memcpy(&A[1][0], &C[1][0], 3 * sizeof(double));
				memcpy(&A[2][0], &C[2][0], 3 * sizeof(double));

				double U[3][3], w[3], VT[3][3];
				vtkMath::SingularValueDecomposition3x3(A, U, w, VT);

				double UVT[3][3];
				vtkMath::Multiply3x3(U, VT, UVT);
				
				double** rot = new double*[3];
				for (int j = 0; j < 3; j++)
					rot[j] = new double[3];
				memcpy(&rot[0][0], &UVT[0][0], 3 * sizeof(double));
				memcpy(&rot[1][0], &UVT[1][0], 3 * sizeof(double));
				memcpy(&rot[2][0], &UVT[2][0], 3 * sizeof(double));

				if (isControlPoint->GetValue(i) == 2) // if this is tumor center point
				{
					memcpy(&tumorRotate[0][0], &rot[0][0], 3 * sizeof(double));
					memcpy(&tumorRotate[1][0], &rot[1][0], 3 * sizeof(double));
					memcpy(&tumorRotate[2][0], &rot[2][0], 3 * sizeof(double));
				}

				double** RelaxShape_local_rot = new double*[Connections[i].size()];
				for (unsigned int j = 0; j < Connections[i].size(); j++)
				{
					RelaxShape_local_rot[j] = new double[3];
				}
				
				vtkMath::MultiplyMatrix(RelaxShape_local, rot, (unsigned int)Connections[i].size(), 3, 3, 3, RelaxShape_local_rot);
				
				// calculate force sent from i
				for (unsigned int j = 0; j < Connections[i].size(); j++)
				{
					double pid = Connections[i][j];
					double CurrentCoord[3], GoalCoord[3];
					memcpy(CurrentCoord, &CurrentShape_local[j][0], 3 * sizeof(double));
					memcpy(GoalCoord, &RelaxShape_local_rot[j][0], 3 * sizeof(double));

					double dir2goal[3];
					vtkMath::Subtract(GoalCoord, CurrentCoord, dir2goal);
					double dis2goal = vtkMath::Norm(dir2goal);
					vtkMath::Normalize(dir2goal);

					double force_i2pid[3];
					double force_i2pidnorm = 5.0 * dis2goal;
					//force_i2pidnorm = force_i2pidnorm > 10.0 ? 10.0 : force_i2pidnorm;

				//	if (isControlPoint->GetValue(i) == 1) // Surface point will have a much larger force to its connections
				//		force_i2pidnorm = 5.0 * force_i2pidnorm;
				//	if (isControlPoint->GetValue(i) == 2) // Jbar point will have a much larger force to its connections
				//		force_i2pidnorm = 20.0 * force_i2pidnorm;

					if (isControlPoint->GetValue(i) == 2 && isControlPoint->GetValue(pid) == 3
						|| isControlPoint->GetValue(i) == 3 && isControlPoint->GetValue(pid) == 2) // forces between tumor and Jbar
					{
						force_i2pidnorm = 1.0 * force_i2pidnorm;
					}

					for (int l = 0; l < 3; l++)	force_i2pid[l] = force_i2pidnorm * dir2goal[l];

					forces[pid].x = forces[pid].x + force_i2pid[0];
					forces[pid].y = forces[pid].y + force_i2pid[1];
					forces[pid].z = forces[pid].z + force_i2pid[2];
					forces[i].x = forces[i].x - force_i2pid[0];
					forces[i].y = forces[i].y - force_i2pid[1];
					forces[i].z = forces[i].z - force_i2pid[2];
				}

				// delete memory
				for (unsigned int j = 0; j < Connections[i].size(); j++)
				{
					delete[] RelaxShape_local[j];
					delete[] CurrentShape_local[j];
					delete[] RelaxShape_local_rot[j];
				}
				delete[] RelaxShape_local;
				delete[] CurrentShape_local;
				delete[] RelaxShape_local_rot;

				for (int j = 0; j < 3; j++)
					delete[] RelaxShape_local_T[j];
				delete[] RelaxShape_local_T;
				for (int j = 0; j < 3; j++)
					delete[] C[j];
				delete[] C;
				for (int j = 0; j < 3; j++)
					delete[] rot[j];
				delete[] rot;

			}

			//// add boundary force
			//for (int i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			//{
			//	if (isControlPoint->GetValue(i) == 1)
			//		continue;
			//	double coordi[3];
			//	SamplePoly->GetPoint(i, coordi);
			//	vtkIdType nearestboundaryPID = boundarypointLocator->FindClosestPoint(coordi);
			//	double boundarycoord[3];
			//	BoundaryPoly->GetPoint(nearestboundaryPID, boundarycoord);
			//	double dir_p2boundary[3];
			//	vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
			//	double dis_p2boundary = vtkMath::Norm(dir_p2boundary);
			//	vtkMath::Normalize(dir_p2boundary);
			//
			//	double boundarynormal[3];
			//	BoundaryNormalArray->GetTuple(nearestboundaryPID, boundarynormal);
			//
			//	double thisboundaryforce_norm = 0.0;
			//	if (vtkMath::Dot(dir_p2boundary, boundarynormal) > 0)
			//		thisboundaryforce_norm = 0.0;
			//	else
			//		thisboundaryforce_norm = 100.0;
			//
			//	for (int l = 0; l < 3; l++)
			//	{
			//		forces[i][l] = forces[i][l] + thisboundaryforce_norm * dir_p2boundary[l];
			//	}
			//}

			// move points according to forces
			for (int i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				if (isControlPoint->GetValue(i) == controltag1
					|| isControlPoint->GetValue(i) == controltag2)
				{
					double GoalCoord[3];
					ControlPointCoord->GetTuple(i, GoalCoord);
					double CurrentCoord[3];
					SamplePoly->GetPoint(i, CurrentCoord);

					double dir2goal[3];
					vtkMath::Subtract(GoalCoord, CurrentCoord, dir2goal);
					double dis2goal = vtkMath::Norm(dir2goal);
					vtkMath::Normalize(dir2goal);

					double Fi[3];
					double force_i2pidnorm = 20 * dis2goal;
				//	force_i2pidnorm = force_i2pidnorm > 50.0 ? 50.0 : force_i2pidnorm;
					for (int l = 0; l < 3; l++)	Fi[l] = force_i2pidnorm * dir2goal[l];

					double coordnew[3];
					for (int l = 0; l < 3; l++) coordnew[l] = CurrentCoord[l] + Fi[l] / POINTMASS * TIMESTEP;

					SamplePoly->GetPoints()->SetPoint(i, coordnew);

					if (isControlPoint->GetValue(i) == 3) // if this is a J-bar center point
					{
						JbarsphereSource->SetCenter(coordnew);
						JbarsphereSource->Update();
					}
				}
				else
				{
					double Fi[3];
					Fi[0] = forces[i].x;
					Fi[1] = forces[i].y;
					Fi[2] = forces[i].z;

					//if (vtkMath::Norm(Fi) > 50.0)
					//{
					//	vtkMath::Normalize(Fi);
					//	vtkMath::MultiplyScalar(Fi, 50.0);
					//}

					double coordi[3];
					SamplePoly->GetPoint(i, coordi);
					double movei[3];
					for (int l = 0; l < 3; l++) movei[l] = Fi[l] / POINTMASS * TIMESTEP;
					vtkMath::Add(coordi, movei, coordi);
					SamplePoly->GetPoints()->SetPoint(i, coordi);

					if (isControlPoint->GetValue(i) == 2) // if this is a tumor center point
					{
						vtkMath::Transpose3x3(tumorRotate, tumorRotate);
						double R_new[3][3];
						vtkMath::Multiply3x3(tumorRotate, tumorScaledRotate_benchmark, R_new);

						double T_new[3];
						for (int l = 0; l < 3; l++)
						{
							T_new[l] = coordi[l] - vtkMath::Dot(R_new[l], tumorCenter_benchmark);
						}
					
						vtkMatrix4x4* m = tumorTransform->GetMatrix();
						vtkSmartPointer<vtkMatrix4x4> newm = vtkSmartPointer<vtkMatrix4x4>::New();
						for (int l = 0; l < 3; l++)
							for (int k = 0; k < 3; k++)
								newm->SetElement(l, k, R_new[l][k]);

						for (int l = 0; l < 3; l++)	newm->SetElement(l, 3, T_new[l]);
						for (int l = 0; l < 3; l++)	newm->SetElement(3, l, 0.0);
						newm->SetElement(3, 3, 1.0);

						tumorTransform->SetMatrix(newm);

					}
					else if (isControlPoint->GetValue(i) == 3) // if this is a J-bar center point
					{
						JbarsphereSource->SetCenter(coordi);
						JbarsphereSource->Update();
					}
				}				
			}
		}
	
		return true;
	}

	bool CubeRedistributionandBuildConnections()
	{
		const double cubestep = CUBESTEP;

		vtkSmartPointer<vtkPoints> SamplePoints = SamplePoly->GetPoints(); // vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkIdTypeArray> xyzindex = vtkSmartPointer<vtkIdTypeArray>::New();
		xyzindex->SetNumberOfComponents(3);

		int DPts = int(2.0 * CRADIUS / cubestep + 1);
		vtkIdType idx = 0;
		for (vtkIdType xidx = 0; xidx < DPts; xidx ++)
			for (vtkIdType yidx = 0; yidx < DPts; yidx++)
				for (vtkIdType zidx = 0; zidx < DPts; zidx++)
				{
					double x = xidx * cubestep - CRADIUS;
					double y = yidx * cubestep - CRADIUS;
					double z = zidx * cubestep - CRADIUS;
					SamplePoints->SetPoint(idx, x, y, z);

					xyzindex->InsertNextTuple3(xidx, yidx, zidx);

					//if (xidx == 0 || xidx == DPts - 1
					//	|| yidx == 0 || yidx == DPts - 1
					//	|| zidx == 0 || zidx == DPts - 1)
					if (zidx == DPts - 1)
					{
						isControlPoint->SetValue(idx, 1);
						double coordtemp[3] = { x, y, z };
						ControlPointCoord->SetTuple(idx, coordtemp);
						R->SetValue(idx, cubestep);
					}
					else if (xidx == int(0.5 * DPts) && yidx == int(0.5 * DPts) && zidx == 0)
					{
						isControlPoint->SetValue(idx, 1);
						double coordtemp[3] = { x, y, z };
						ControlPointCoord->SetTuple(idx, coordtemp);
						R->SetValue(idx, cubestep);
					}
					else
					{
						isControlPoint->SetValue(idx, 0);
						double coordtemp[3] = { 0.0, 0.0, 0.0 };
						ControlPointCoord->SetTuple(idx, coordtemp);
						R->SetValue(idx, cubestep);
					}

					idx++;
				}

		
		// build connections
		this->Connections.clear();
		connectionCellArray->Reset();
		isTumorConnections->Reset();

		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
		{
		/*	double* xyzidx = xyzindex->GetTuple3(i);
			vtkIdType xi = xyzidx[0];
			vtkIdType yi = xyzidx[1];
			vtkIdType zi = xyzidx[2];

		//	std::cout << xi << ", " << yi << ", " << zi << std::endl;

			vector<vtkIdType> Connection_i;
			// 6 connections
			if (xi - 1 >= 0)
			{
				vtkIdType neighorhoodpid = (xi - 1) * DPts * DPts + (yi) * DPts + zi;
				Connection_i.push_back(neighorhoodpid);
			}
			if (xi + 1 < DPts)
			{
				vtkIdType neighorhoodpid = (xi + 1) * DPts * DPts + (yi) * DPts + zi;
				Connection_i.push_back(neighorhoodpid);
			}
			if (yi - 1 >= 0)
			{
				vtkIdType neighorhoodpid = (xi)* DPts * DPts + (yi - 1) * DPts + zi;
				Connection_i.push_back(neighorhoodpid);
			}
			if (yi + 1 < DPts)
			{
				vtkIdType neighorhoodpid = (xi)* DPts * DPts + (yi + 1)* DPts + zi;
				Connection_i.push_back(neighorhoodpid);
			}
			if (zi - 1 >= 0)
			{
				vtkIdType neighorhoodpid = (xi)* DPts * DPts + (yi) * DPts + zi - 1;
				Connection_i.push_back(neighorhoodpid);
			}
			if (zi + 1 < DPts)
			{
				vtkIdType neighorhoodpid = (xi)* DPts * DPts + (yi)* DPts + zi + 1;
				Connection_i.push_back(neighorhoodpid);
			}
		*/
			
			double coordi[3];
			SamplePoly->GetPoints()->GetPoint(i, coordi);
			double Ri = R->GetValue(i);

			vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New();
			//pointLocator->FindPointsWithinRadius(4.0, coordi, NeighorpIds);
			pointLocator->FindClosestNPoints(20, coordi, NeighorpIds);

			vector<vtkIdType> Connection_i;

			for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj++)
			{
				vtkIdType j = NeighorpIds->GetId(idxj);
				if (i == j) continue;

				double coordj[3];
				SamplePoly->GetPoints()->GetPoint(j, coordj);

				double dir[3];
				vtkMath::Subtract(coordi, coordj, dir);
				double dis = vtkMath::Norm(dir);

				double Rj = R->GetValue(j);
				double R_total = Ri + Rj;

				if (dis < 0.9 * R_total)
					Connection_i.push_back(j);
			}
			
			this->Connections.push_back(Connection_i);
		}

		// build connectionCellArray for display
		for (vtkIdType i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
		{
			double* xyzidx = xyzindex->GetTuple3(i);
			double zi = xyzidx[2];
			//std::cout << i << ", " << xyzidxi[0] << ", " << xyzidxi[1] << ", " << xyzidxi[2] << std::endl;

			for (unsigned int idj = 0; idj < this->Connections[i].size(); idj++)
			{
				vtkIdType j = this->Connections[i].at(idj);
						
				if (j <= i)	continue;

				xyzidx = xyzindex->GetTuple3(j);
				double zj = xyzidx[2];
				if (int(zi) != int(zj))
				{
					continue;
				}

				vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
				line->GetPointIds()->SetId(0, i);
				line->GetPointIds()->SetId(1, j);
				connectionCellArray->InsertNextCell(line);

				if (isControlPoint->GetValue(i) == 1 && isControlPoint->GetValue(j) == 1)
					isTumorConnections->InsertNextValue(1);
				else
					isTumorConnections->InsertNextValue(0);

			}
		}

		return true;
	}

	virtual void OnKeyPress()
	{
		std::string key = this->Interactor->GetKeySym();
	
		if (key == "s") // start to uniformly sampling
		{
		//	if (InitialPointcloud() == false)
		//		return;

			if (UniformRedistribution(10) == true)
			{
				BuildConnections();
				FindSurfaceControlPoints();

				connectionCellArray->Modified();
				connectionPolyData->Modified();
				SamplePoly->Modified();
//				tumorPoly->Modified();
				renderWindow->Render();

				RelaxShape->DeepCopy(SamplePoly->GetPoints());
				//	SavePolyData(SamplePoly, "C:\\work\\smooth_deformation_3D\\testdata\\SamplePoly.vtp");
			}
		}

		else if (key == "g")
		{
			forces.resize(SamplePoly->GetPoints()->GetNumberOfPoints());
			for (int iter = 0; iter < 10; iter++)
			{
				DeformationMotion(1, -1);
			}

			//const double cubestep = CUBESTEP;
			//int DPts = int(2.0 * CRADIUS / cubestep + 1);
			//double TuberCenterCoord[3];
			//SamplePoly->GetPoint(TumorCenterIdx[0] * DPts * DPts + TumorCenterIdx[1] * DPts + TumorCenterIdx[2], TuberCenterCoord);
			//TumorSphereSource->SetCenter(TuberCenterCoord);
			//TumorSphereSource->Update();

			connectionCellArray->Modified();
			connectionPolyData->Modified();
			SamplePoly->Modified();
			tumorTransformFilter->GetOutput()->Modified();
			renderWindow->Render();
		}

		else if (key == "z" || key == "x" || key == "i" || key == "k" || key == "j" || key == "l")
		{
			for (int i = 0; i < SamplePoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				if (isControlPoint->GetValue(i) == 3)
				{
					double coordjbar[3];
					SamplePoly->GetPoint(i, coordjbar);
					const double movestep = 0.05;
					if (key == "z")
						coordjbar[0] -= movestep;
					else if (key == "x")
						coordjbar[0] += movestep;
					else if (key == "i")
						coordjbar[1] += movestep;
					else if (key == "k")
						coordjbar[1] -= movestep;
					else if (key == "j")
						coordjbar[2] -= movestep;
					else if (key == "l")
						coordjbar[2] += movestep;
					SamplePoly->GetPoints()->SetPoint(i, coordjbar);
					ControlPointCoord->SetTuple(i, coordjbar);
				}
			}

			forces.resize(SamplePoly->GetPoints()->GetNumberOfPoints());
			for (int iter = 0; iter < 8; iter++)
			{
				DeformationMotion(1, 3);
			}

			connectionCellArray->Modified();
			connectionPolyData->Modified();
			SamplePoly->Modified();
			//JbarsphereSource->;
			tumorTransformFilter->GetOutput()->Modified();
			renderWindow->Render();
		}

		else if (key == "Left" || key == "Right")
		{
			if (ControlPointCoord == NULL)
				return;

			for (int i = 0; i < BoundaryPoly->GetPoints()->GetNumberOfPoints(); i++)
			{
				double boundarycoordi[3];
				BoundaryPoly->GetPoint(i, boundarycoordi);
				if (key == "Left")
				{
					for (int l = 0; l < 1; l++)
						boundarycoordi[l] = 0.95 * boundarycoordi[l];
				}
				else
				{
					for (int l = 0; l < 1; l++)
						boundarycoordi[l] = 1.05 * boundarycoordi[l];
				}
				BoundaryPoly->GetPoints()->SetPoint(i, boundarycoordi);
			}
			BoundaryPoly->GetPoints()->Modified();
			BoundaryPoly->Modified();
			renderWindow->Render();

			// BoundaryPolynormalGenerator->Update();

			for (int i = 0; i < ControlPointCoord->GetNumberOfTuples(); i++)
			{
				if (isControlPoint->GetValue(i) != 1)
					continue;

				double newcontrolcoord[3];
				BoundaryPoly->GetPoint(RelatedBoundaryPids[i], newcontrolcoord);
				ControlPointCoord->SetTuple(i, newcontrolcoord);
			}
		}
	}

	bool Pick(double picked[3])
	{
		int x = this->Interactor->GetEventPosition()[0];
		int y = this->Interactor->GetEventPosition()[1];

		if (this->CurrentRenderer == NULL) return false;

		this->Interactor->GetPicker()->Pick(x, y, 0, this->CurrentRenderer);
		this->Interactor->GetPicker()->GetPickPosition(picked);
		return true;
	}
	
	virtual void OnLeftButtonDown()
	{
		this->Superclass::OnLeftButtonDown();

		if (!SamplePoly && !BoundaryPoly)
			return;

		if (this->Pick(lastpickpos))
		{
			this->LeftButtonDown = true;
			this->PickedBoundaryPID = boundarypointLocator->FindClosestPointWithinRadius(2.0, lastpickpos, this->pickboundarydistance);
			this->PickedSamplePID = pointLocator->FindClosestPointWithinRadius(2.0, lastpickpos, this->picksampledistance);
		//	std::cout << "this->PickedBoundaryPID = " << this->PickedBoundaryPID << std::endl;
		}
	}

	virtual void OnLeftButtonUp()
	{
		this->Superclass::OnLeftButtonUp();

		this->PickedBoundaryPID = -1;
		this->PickedSamplePID = -1;
		this->LeftButtonDown = false;
	}

	virtual void OnMouseMove()
	{
		if (this->Interactor->GetControlKey() == false
			&& this->Interactor->GetShiftKey() == false)
		{
			this->Superclass::OnMouseMove();
			return;
		}

		if (this->LeftButtonDown == false)
			return;

		if (this->Interactor->GetControlKey() != false)
		{
			if (this->PickedBoundaryPID == -1)
				return;

		//	vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPolynormalGenerator->GetOutput()->GetPointData()->GetArray("Normals"));
			if (BoundaryNormalArray == NULL)
			{
				std::cerr << "cannot find BoundaryNormalArray" << std::endl;
				return;
			}

			if (Pick(pickpos))
			{
				//if (this->pickboundarydistance > 5.0)
				//{
				//	std::cout << "this->pickboundarydistance > 5.0" << std::endl;
				//	return;
				//}
				double move[3];
				vtkMath::Subtract(pickpos, lastpickpos, move);
						
				//double boundarynormal[3];
				//BoundaryNormalArray->GetTuple(this->PickedBoundaryPID, boundarynormal);
				//double moveproj = vtkMath::Dot(move, boundarynormal);
				//for (int l = 0; l < 3; l++) move[l] = moveproj * boundarynormal[l];

				double PickedBoundaryCoord[3];
				BoundaryPoly->GetPoint(this->PickedBoundaryPID, PickedBoundaryCoord);

				vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New();
				boundarypointLocator->FindPointsWithinRadius(10.0, PickedBoundaryCoord, NeighorpIds);

				for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj++)
				{
					vtkIdType j = NeighorpIds->GetId(idxj);
					double coordj[3];
					BoundaryPoly->GetPoint(j, coordj);
					double boundarynormalj[3];
					BoundaryNormalArray->GetTuple(j, boundarynormalj);

					//if (vtkMath::Dot(boundarynormal, boundarynormalj) < 0.0)
					//	continue;

					double dis = sqrt(vtkMath::Distance2BetweenPoints(coordj, PickedBoundaryCoord));
					double w = exp(-0.2 * (dis * dis));
					double wmove[3];
					for (int l = 0; l < 3; l++) wmove[l] = w * move[l];

					//double moveprojj = vtkMath::Dot(wmove, boundarynormalj);
					//for (int l = 0; l < 3; l++) wmove[l] = moveprojj * boundarynormalj[l];

					vtkMath::Add(coordj, wmove, coordj);
					BoundaryPoly->GetPoints()->SetPoint(j, coordj);
				}

				BoundaryPoly->GetPoints()->Modified();
				BoundaryPoly->Modified();
				renderWindow->Render();

				for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj++)
				{
					vtkIdType j = NeighorpIds->GetId(idxj);
					for (int i = 0; i < ControlPointCoord->GetNumberOfTuples(); i++)
					{
						if (RelatedBoundaryPids[i] != j)
							continue;

						double newcontrolcoord[3];
						BoundaryPoly->GetPoint(RelatedBoundaryPids[i], newcontrolcoord);
						ControlPointCoord->SetTuple(i, newcontrolcoord);
					}
				}

				std::swap(lastpickpos, pickpos);

				//forces.resize(SamplePoly->GetPoints()->GetNumberOfPoints());
				//for (int iter = 0; iter < 10; iter++)
				//{
				//	DeformationMotion(1, 2);
				//}
				//connectionCellArray->Modified();
				//connectionPolyData->Modified();
				//SamplePoly->Modified();
				//renderWindow->Render();
			}
		}

		if (this->Interactor->GetShiftKey() != false)
		{
			if (this->PickedSamplePID == -1)
				return;

			if (isControlPoint->GetValue(this->PickedSamplePID) == 0)
				return;

			if (Pick(pickpos))
			{
				if (this->picksampledistance > 1.0)
					return;

				double move[3];
				vtkMath::Subtract(pickpos, lastpickpos, move);
				
				double PickedSampleCoord[3];
				SamplePoly->GetPoint(this->PickedSamplePID, PickedSampleCoord);

				vtkSmartPointer<vtkIdList> NeighorpIds = vtkSmartPointer<vtkIdList>::New();
				pointLocator->FindPointsWithinRadius(20.0, PickedSampleCoord, NeighorpIds);

				for (int idxj = 0; idxj < NeighorpIds->GetNumberOfIds(); idxj++)
				{
					vtkIdType j = NeighorpIds->GetId(idxj);
					if (isControlPoint->GetValue(j) == 0) continue;

					double coordj[3];
					SamplePoly->GetPoint(j, coordj);

					double dis = sqrt(vtkMath::Distance2BetweenPoints(coordj, PickedSampleCoord));
					double w = exp(-0.2 * (dis * dis));
					double wmove[3];
					for (int l = 0; l < 3; l++) wmove[l] = w * move[l];
	
					vtkMath::Add(coordj, wmove, coordj);
					SamplePoly->GetPoints()->SetPoint(j, coordj);

					ControlPointCoord->SetTuple(j, coordj);
				}

				std::swap(lastpickpos, pickpos);

				forces.resize(SamplePoly->GetPoints()->GetNumberOfPoints());
				for (int iter = 0; iter < 10; iter++)
				{
					DeformationMotion(1, 2);
				}
				connectionCellArray->Modified();
				connectionPolyData->Modified();
				SamplePoly->Modified();
				renderWindow->Render();
			}
		}

	}

public:

	// system
	int ClickCount;
	vtkRenderWindow* renderWindow;
	string dataname;
	
	// boundary
	vtkPolyData* BoundaryPoly;
	vtkPointLocator* boundarypointLocator;	
	vtkSmartPointer<vtkPolyDataNormals> BoundaryPolynormalGenerator;
	vtkFloatArray* BoundaryNormalArray;

	// p
	vtkPointLocator* pointLocator;
	vtkSmartPointer<vtkPolyData> SamplePoly;
	vtkDoubleArray* R;	// radius
	vtkDoubleArray* F_sumabs;	// force with respect to other points
	vtkDoubleArray* F_abssum;	// force with respect to other points
	vtkSmartPointer<vtkPoints> RelaxShape;

	vtkSmartPointer<vtkIdTypeArray> isControlPoint;
	vtkSmartPointer<vtkDoubleArray> ControlPointCoord;
	vector<vtkIdType> RelatedBoundaryPids;

	vector< vector<vtkIdType> > Connections;
//	vector< vector<vtkIdType> > relationships; // each one has two points ids.
//	vector< vector<double> > parameters; // each one has some parameters for calculating forces
	vector< CForce > forces;		// forces used in phrase 2
	
	vtkPolyData* connectionPolyData;
	vtkCellArray* connectionCellArray;
	vtkIdTypeArray* isTumorConnections;

	//
	bool LeftButtonDown;

	double pickpos[3];
	double lastpickpos[3];

	double pickboundarydistance;
	vtkIdType PickedBoundaryPID;

	double picksampledistance;
	vtkIdType PickedSamplePID;

	vtkTransformPolyDataFilter* tumorTransformFilter;
	vtkTransform* tumorTransform;
	double tumorRotate[3][3];
	
	double tumorScaledRotate_benchmark[3][3];
	double tumorTranslate_benchmark[3];
	double tumorCenter_benchmark[3];

	vtkSphereSource* JbarsphereSource;
};
vtkStandardNewMacro(MouseInteractorStyle);


int main(int argc, char *argv[])
{	
	vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(CRADIUS);
	sphereSource->SetPhiResolution(60);
	sphereSource->SetThetaResolution(60);
	sphereSource->Update();
	vtkSmartPointer<vtkPolyData> BoundaryPoly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkTransform> Boundarytranslation = vtkSmartPointer<vtkTransform>::New();
	vtkSmartPointer<vtkTransformPolyDataFilter> BoundaryTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();

	{
		vtkSmartPointer< vtkXMLPolyDataReader > reader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
		reader->SetFileName("C:\\work\\Lung_Models\\FMA7333.vtp");
		try
		{
			reader->Update();
		}
		catch (...)
		{
			std::cerr << "Error occurs when reading" << std::endl;
			return 0;
		}

		//	BoundaryPoly->DeepCopy(reader->GetOutput());

		vtkSmartPointer<vtkTriangleFilter> trianglefilter = vtkSmartPointer<vtkTriangleFilter>::New();
		trianglefilter->SetInputData(reader->GetOutput());
		trianglefilter->Update();

		vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanPolyData->SetInputData(trianglefilter->GetOutput());
		cleanPolyData->Update();
		BoundaryPoly->DeepCopy(cleanPolyData->GetOutput());

		double* BoundaryBounds = BoundaryPoly->GetBounds();
		double BoundaryScaleX = BoundaryBounds[1] - BoundaryBounds[0];
		double BoundaryScaleY = BoundaryBounds[3] - BoundaryBounds[2];
		double BoundaryScaleZ = BoundaryBounds[5] - BoundaryBounds[4];
		double BoundaryScale = vtkMath::Max(BoundaryScaleX, BoundaryScaleY);
		BoundaryScale = vtkMath::Max(BoundaryScale, BoundaryScaleZ);
		double* BoundaryCenter = BoundaryPoly->GetCenter();

		// change scale from BoundaryScale to 2*CRADIUS, change center from BoundaryCenter to [0,0,0]
		Boundarytranslation->Scale(2 * CRADIUS / BoundaryScale, 2 * CRADIUS / BoundaryScale, 2 * CRADIUS / BoundaryScale);
		Boundarytranslation->Translate(-BoundaryCenter[0], -BoundaryCenter[1], -BoundaryCenter[2]);
		BoundaryTransformFilter->SetInputData(BoundaryPoly);
		BoundaryTransformFilter->SetTransform(Boundarytranslation);
		BoundaryTransformFilter->Update();

		//	SavePolyData(BoundaryPoly, "C:\\work\\smooth_deformation_3D\\testdata\\lungboundarypoly1.vtp");
		//	smoothvtkpolydata(BoundaryPoly, 50);
		vtkSmartPointer<vtkLoopSubdivisionFilter> loopsubdivisionfilter = vtkSmartPointer<vtkLoopSubdivisionFilter>::New();
		loopsubdivisionfilter->SetInputData(BoundaryTransformFilter->GetOutput());
		loopsubdivisionfilter->SetNumberOfSubdivisions(1);
		loopsubdivisionfilter->Update();
		BoundaryPoly = loopsubdivisionfilter->GetOutput();
	}

	vtkSmartPointer<vtkPolyDataNormals> BoundaryPolynormalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	BoundaryPolynormalGenerator->SetInputData(BoundaryPoly);
	BoundaryPolynormalGenerator->ComputePointNormalsOn();
	BoundaryPolynormalGenerator->ConsistencyOn();
	BoundaryPolynormalGenerator->AutoOrientNormalsOn();
	BoundaryPolynormalGenerator->Update();
	vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPolynormalGenerator->GetOutput()->GetPointData()->GetArray("Normals"));
	//vtkFloatArray* BoundaryNormalArray = vtkFloatArray::SafeDownCast(BoundaryPoly->GetPointData()->GetArray("Normals"));

	BoundaryPoly->GetPointData()->SetNormals(BoundaryNormalArray);
//	SavePolyData(BoundaryPoly, "C:\\work\\smooth_deformation_3D\\testdata\\lungboundarypoly2.vtp");

	if (BoundaryNormalArray == NULL)
	{
		std::cerr << "cannot find BoundaryNormalArray" << std::endl;
		return false;
	}	
		
	vtkSmartPointer<vtkPointLocator> boundarypointLocator = vtkSmartPointer<vtkPointLocator>::New();
	boundarypointLocator->SetDataSet(BoundaryPoly);
	boundarypointLocator->AutomaticOn();
	boundarypointLocator->SetNumberOfPointsPerBucket(5);
	boundarypointLocator->BuildLocator();
	
	// generate initial SamplePoints
	vtkSmartPointer<vtkPolyData> SamplePoly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> SamplePoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> SampleCell = vtkSmartPointer<vtkCellArray>::New();
	SamplePoly->SetPoints(SamplePoints);
	SamplePoly->SetVerts(SampleCell);

	vtkSmartPointer<vtkIdTypeArray> isControlPoint = vtkSmartPointer<vtkIdTypeArray>::New();
	isControlPoint->SetName("isControlPoint");
	isControlPoint->SetNumberOfComponents(1);
	isControlPoint->SetNumberOfValues(SamplePoints->GetNumberOfPoints());
	vtkSmartPointer<vtkDoubleArray> ControlPointCoord = vtkSmartPointer<vtkDoubleArray>::New();
	ControlPointCoord->SetName("ControlPointCoord");
	ControlPointCoord->SetNumberOfComponents(3);
	ControlPointCoord->SetNumberOfTuples(SamplePoints->GetNumberOfPoints());


	for (int i = 0; i < N; )
	{
		double coordi[3] = {0.0, 0.0, 0.0};
		coordi[0] = vtkMath::Random(-CRADIUS, CRADIUS);
		coordi[1] = vtkMath::Random(-CRADIUS, CRADIUS);
		coordi[2] = vtkMath::Random(-CRADIUS, CRADIUS);
		//	if (vtkMath::Norm(coordi) > CRADIUS) continue;
		
		vtkSmartPointer<vtkIdList> NeighorboundaryIds = vtkSmartPointer<vtkIdList>::New();
		boundarypointLocator->FindClosestNPoints(5, coordi, NeighorboundaryIds);
		double boundarycoord[3];
		BoundaryPoly->GetPoint(NeighorboundaryIds->GetId(0), boundarycoord);
		double dir_p2boundary[3];
		vtkMath::Subtract(boundarycoord, coordi, dir_p2boundary);
		vtkMath::Normalize(dir_p2boundary);
		double boundarynormal[3] = { 0.0, 0.0, 0.0 };
		for (int idxj = 0; idxj < NeighorboundaryIds->GetNumberOfIds(); idxj++)
		{
			vtkIdType j = NeighorboundaryIds->GetId(idxj);
			double boundarynormalj[3];
			BoundaryNormalArray->GetTuple(j, boundarynormalj);
			vtkMath::Add(boundarynormal, boundarynormalj, boundarynormal);
		}
		vtkMath::Normalize(boundarynormal);

		if (vtkMath::Dot(dir_p2boundary, boundarynormal) < 0)
			continue;	

		SamplePoints->InsertNextPoint(coordi);
		isControlPoint->InsertNextValue(0);
		double temp[3] = { 0.0, 0.0, 0.0 };
		ControlPointCoord->InsertNextTuple(temp);
		i ++;
	}
	

//	SavePolyData(SamplePoly, "C:\\work\\smooth_deformation_3D\\testdata\\SamplePoly.vtp");
	
	// insert tumor mesh model
	vtkSmartPointer<vtkPolyData> tumorPoly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkTransform> tumorTransform = vtkSmartPointer<vtkTransform>::New();
	vtkSmartPointer<vtkTransformPolyDataFilter> tumorTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	{
		vtkSmartPointer< vtkXMLPolyDataReader > tumorreader = vtkSmartPointer< vtkXMLPolyDataReader >::New();
		tumorreader->SetFileName("C:\\work\\Lung_Models\\Tumor.vtp");
		try
		{
			tumorreader->Update();
		}
		catch (...)
		{
			std::cerr << "Error occurs when reading" << std::endl;
			return 0;
		}
		vtkSmartPointer<vtkTriangleFilter> trianglefilter = vtkSmartPointer<vtkTriangleFilter>::New();
		trianglefilter->SetInputData(tumorreader->GetOutput());
		trianglefilter->Update();

		vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanPolyData->SetInputData(trianglefilter->GetOutput());
		cleanPolyData->Update();
		tumorPoly = cleanPolyData->GetOutput();

		double* tumorBounds = tumorPoly->GetBounds();
		double tumorScaleX = tumorBounds[1] - tumorBounds[0];
		double tumorScaleY = tumorBounds[3] - tumorBounds[2];
		double tumorScaleZ = tumorBounds[5] - tumorBounds[4];
		double tumorScale_Original = vtkMath::Max(tumorScaleX, tumorScaleY);
		tumorScale_Original = vtkMath::Max(tumorScale_Original, tumorScaleZ);
		double* tumorCenter = tumorPoly->GetCenter();

		// change scale from tumorScale totumorgoalscale, change center from tumorCenter to [0,0,0]
		const double tumorgoalscale = 1.0;
		tumorTransform->Scale(tumorgoalscale / tumorScale_Original, tumorgoalscale / tumorScale_Original, tumorgoalscale / tumorScale_Original);
		tumorTransform->Translate(15.0 - tumorCenter[0], 15.0 - tumorCenter[1], 35.0 - tumorCenter[2]);
		
		//vtkMatrix4x4* m = tumorTransform->GetMatrix();
		//std::cout << m->GetElement(0, 0) << ", " << m->GetElement(0, 1) << ", " << m->GetElement(0, 2) << ", " << m->GetElement(0, 3) << std::endl;
		//std::cout << m->GetElement(1, 0) << ", " << m->GetElement(1, 1) << ", " << m->GetElement(1, 2) << ", " << m->GetElement(1, 3) << std::endl;
		//std::cout << m->GetElement(2, 0) << ", " << m->GetElement(2, 1) << ", " << m->GetElement(2, 2) << ", " << m->GetElement(2, 3) << std::endl;
		//std::cout << m->GetElement(3, 0) << ", " << m->GetElement(3, 1) << ", " << m->GetElement(3, 2) << ", " << m->GetElement(3, 3) << std::endl;
		
		tumorTransformFilter->SetInputData(tumorPoly);
		tumorTransformFilter->SetTransform(tumorTransform);
		tumorTransformFilter->Update();
		tumorCenter = tumorTransformFilter->GetOutput()->GetCenter();

	//	SavePolyData(tumorTransformFilter->GetOutput(), "C:\\work\\smooth_deformation_3D\\testdata\\tumorTransformFilter.vtp");
		{
			SamplePoints->InsertNextPoint(tumorCenter);
			isControlPoint->InsertNextValue(2);
			ControlPointCoord->InsertNextTuple(tumorCenter);
		}
	}

	// insert the J-bar 
	vtkSmartPointer<vtkSphereSource> JbarsphereSource = vtkSmartPointer<vtkSphereSource>::New();
	{
		double* tumorCenter = tumorTransformFilter->GetOutput()->GetCenter();
		double coord_jbar[3] = { tumorCenter[0] + 0.5, tumorCenter[1] - 0.5, tumorCenter[2]};
		SamplePoints->InsertNextPoint(coord_jbar);
		isControlPoint->InsertNextValue(3);
		ControlPointCoord->InsertNextTuple(coord_jbar);

		JbarsphereSource->SetCenter(coord_jbar);
		JbarsphereSource->SetRadius(0.2);
		JbarsphereSource->SetPhiResolution(6);
		JbarsphereSource->SetThetaResolution(6);
		JbarsphereSource->Update();
	}

	for (vtkIdType i = 0; i < SamplePoints->GetNumberOfPoints(); i++)
	{
		SampleCell->InsertNextCell(1, &i);
	}

	// add color lookup table
	SamplePoly->GetCellData()->SetScalars(isControlPoint);
	vtkSmartPointer<vtkLookupTable> lookupTable = vtkSmartPointer<vtkLookupTable>::New();
	lookupTable->SetNumberOfTableValues(4);
	lookupTable->SetRange(0.0, 3.0);
	lookupTable->SetScaleToLinear();
	lookupTable->SetTableValue(0, 1.0, 0.0, 0.0, 1); 
	lookupTable->SetTableValue(1, 0.0, 0.0, 1.0, 1); 
	lookupTable->SetTableValue(2, 0.0, 0.0, 0.0, 1);
	lookupTable->SetTableValue(3, 1.0, 0.0, 1.0, 1);
	lookupTable->Build();

	vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
	pointLocator->SetDataSet(SamplePoly);
	pointLocator->AutomaticOn();
	pointLocator->SetNumberOfPointsPerBucket(5);
	pointLocator->BuildLocator();
	
	vtkSmartPointer<vtkDoubleArray> R = vtkSmartPointer<vtkDoubleArray>::New(); // radius
	R->SetName("radius");
	R->SetNumberOfComponents(1);
	vtkSmartPointer<vtkDoubleArray> F_abssum = vtkSmartPointer<vtkDoubleArray>::New(); // radius
	F_abssum->SetName("F_abssum");
	F_abssum->SetNumberOfComponents(1);
	vtkSmartPointer<vtkDoubleArray> F_sumabs = vtkSmartPointer<vtkDoubleArray>::New(); // radius
	F_sumabs->SetName("F_sumabs");
	F_sumabs->SetNumberOfComponents(1);

	for (int i = 0; i < SamplePoints->GetNumberOfPoints(); i ++)
	{
		R->InsertNextValue(R_INITIAL);	// R is initialized with a small number
		F_abssum->InsertNextValue(0.0);
		F_sumabs->InsertNextValue(0.0);
	}

	SamplePoly->GetPointData()->AddArray(R);
	SamplePoly->GetPointData()->AddArray(F_abssum);
	SamplePoly->GetPointData()->AddArray(F_sumabs);

	// 
	vtkSmartPointer<vtkCellArray> connectionCellArray =	vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPolyData> connectionPolyData = vtkSmartPointer<vtkPolyData>::New();
	connectionPolyData->SetPoints(SamplePoints);
	connectionPolyData->SetLines(connectionCellArray);
	vtkSmartPointer<vtkIdTypeArray> isTumorConnections = vtkSmartPointer<vtkIdTypeArray>::New();
	connectionPolyData->GetCellData()->SetScalars(isTumorConnections);

	// add connections color lookup table
	vtkSmartPointer<vtkLookupTable> connectionslookupTable = vtkSmartPointer<vtkLookupTable>::New();
	connectionslookupTable->SetNumberOfTableValues(3);
	connectionslookupTable->SetRange(0.0, 2.0);
	connectionslookupTable->SetScaleToLinear();
	connectionslookupTable->SetTableValue(0, 0.0, 1.0, 0.0, 1);
	connectionslookupTable->SetTableValue(1, 0.0, 0.0, 0.0, 1);
	connectionslookupTable->SetTableValue(2, 1.0, 1.0, 0.0, 1);
	connectionslookupTable->Build();

	// Render window
	vtkSmartPointer<vtkRenderWindow> renderWindow =	vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->SetSize(200* 1,200); //(width, height)
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkRenderer> renderer =	vtkSmartPointer<vtkRenderer>::New();
	renderWindow->AddRenderer(renderer);
	renderer->SetViewport(static_cast<double>(0)/1,0,static_cast<double>(0+1)/1,1);
	vtkSmartPointer<vtkPolyDataMapper> mapper =	vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(SamplePoly); 
	mapper->SetScalarRange(0.0, 3.0);
	mapper->SetLookupTable(lookupTable);
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
	actor->GetProperty()->SetPointSize(5.0);
//	renderer->AddActor(actor);

	// the boundary lung model
	vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper1->SetInputData(BoundaryPoly);
	vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
	actor1->SetMapper(mapper1);
	actor1->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
	actor1->GetProperty()->SetOpacity(0.1);
	actor1->GetProperty()->SetDiffuse(1);
	actor1->GetProperty()->SetSpecular(1);
	renderer->AddActor(actor1);

	// the connections
	vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New(); 
	mapper2->SetInputData(connectionPolyData);
	mapper2->SetScalarRange(0.0, 2.0);
	mapper2->SetLookupTable(connectionslookupTable);
	vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
	actor2->SetMapper(mapper2);
	actor2->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
	actor2->GetProperty()->SetLineWidth(0.3);
//	actor2->GetProperty()->SetOpacity(0.5);
	renderer->AddActor(actor2);


	// tumor model
	vtkSmartPointer<vtkPolyDataMapper> mapper_tumor = vtkSmartPointer<vtkPolyDataMapper>::New(); // the tumor sphere
	mapper_tumor->SetInputConnection(tumorTransformFilter->GetOutputPort());
	vtkSmartPointer<vtkActor> actor_tumor = vtkSmartPointer<vtkActor>::New();
	actor_tumor->SetMapper(mapper_tumor);
	actor_tumor->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
	actor_tumor->GetProperty()->SetOpacity(1.0);
	renderer->AddActor(actor_tumor);

	// J-bar
	vtkSmartPointer<vtkPolyDataMapper> mapper_Jbar = vtkSmartPointer<vtkPolyDataMapper>::New(); // the tumor sphere
	mapper_Jbar->SetInputConnection(JbarsphereSource->GetOutputPort());
	vtkSmartPointer<vtkActor> actor_Jbar = vtkSmartPointer<vtkActor>::New();
	actor_Jbar->SetMapper(mapper_Jbar);
	actor_Jbar->GetProperty()->SetColor(1.0, 1.0, 0.0); //(R,G,B)
	actor_Jbar->GetProperty()->SetOpacity(1.0);
	renderer->AddActor(actor_Jbar);


	renderer->SetBackground(1.0, 1.0, 1.0);
	renderer->SetAutomaticLightCreation(1);

	renderer->ResetCamera();
	renderWindow->Render();
	renderWindow->SetWindowName("Show Smoothed Points");

	vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
	vtkSmartPointer<vtkOrientationMarkerWidget> axeswidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
	axeswidget->SetOutlineColor(0.9300, 0.5700, 0.1300);
	axeswidget->SetOrientationMarker(axes);
	axeswidget->SetInteractor(renderWindowInteractor);
	axeswidget->SetViewport(0.0, 0.0, 0.4, 0.4);
	axeswidget->SetEnabled(1);
	axeswidget->InteractiveOff();

	vtkSmartPointer<vtkPointPicker> pointPicker = vtkSmartPointer<vtkPointPicker>::New();
	pointPicker->SetTolerance(1e-5);
	renderWindowInteractor->SetPicker(pointPicker);

	vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
	style->BoundaryPoly = BoundaryPoly;
	style->boundarypointLocator = boundarypointLocator;
//	style->BoundaryPolynormalGenerator = BoundaryPolynormalGenerator;
	style->BoundaryNormalArray = BoundaryNormalArray;

	style->pointLocator = pointLocator;

	style->SamplePoly = SamplePoly;
	style->R = R;
	style->F_abssum = F_abssum;
	style->F_sumabs = F_sumabs;

	style->connectionPolyData = connectionPolyData;
	style->connectionCellArray = connectionCellArray;
	style->isTumorConnections = isTumorConnections;
	
	style->isControlPoint = isControlPoint;
	style->ControlPointCoord = ControlPointCoord;

	style->tumorTransform = tumorTransform;
	style->tumorTransformFilter = tumorTransformFilter;
	
	style->JbarsphereSource = JbarsphereSource;

	vtkMatrix4x4* m = tumorTransform->GetMatrix();
	for (int l = 0; l < 3; l++)
		for (int k = 0; k < 3; k++)
			style->tumorScaledRotate_benchmark[l][k] = m->GetElement(l, k);

	for (int l = 0; l < 3; l++) style->tumorTranslate_benchmark[l] = m->GetElement(l, 3);
	for (int l = 0; l < 3; l++) style->tumorCenter_benchmark[l] = tumorPoly->GetCenter()[l];

//	std::cout << "tumorCenter_benchmark = " << style->tumorCenter_benchmark[0] << ", " << style->tumorCenter_benchmark[1] << ", " << style->tumorCenter_benchmark[2] << std::endl;

	style->renderWindow = renderWindow;
	style->ClickCount = 0;

	renderWindowInteractor->SetInteractorStyle( style );		
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}

void SaveVTKImage(vtkImageData *image, const char* fileName)
{
	vtkSmartPointer< vtkMetaImageWriter > writer = vtkSmartPointer< vtkMetaImageWriter >::New();
	writer->SetFileName(fileName);
	writer->SetInputData(image);
	try
	{
		writer->Write();
	}
	catch(...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}

void SavePolyData(vtkPolyData *poly, const char* fileName)
{
	if (!poly) return;
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(poly);
	writer->SetFileName(fileName);
	writer->SetDataModeToBinary();
	try
	{
		writer->Write();
	}
	catch (...)
	{
		std::cerr << "Error occurs when writing " << fileName << std::endl;
		return;
	}
}

int smoothvtkpolydata(vtkPolyData* Poly, int iternum, int TYPE)
{
	vtkPoints* Points = Poly->GetPoints();
	vtkCellArray* Strips = NULL;
	if (TYPE == 1)
		Strips = Poly->GetPolys();
	else
		Strips = Poly->GetStrips();

	vtkIdType CellId = 0;
	vtkIdType npts = 0, *pts = NULL;

	int* adjcent = (int*)malloc(Points->GetNumberOfPoints() * 10 * sizeof(int) );
	int* num_adjcent = (int*)malloc(Points->GetNumberOfPoints() * sizeof(int));


/*	for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid ++)
	{
		num_adjcent[pid] = 0;
		for(CellId = 0, Strips->InitTraversal(); Strips->GetNextCell(npts, pts); CellId ++)
		{
			if (pts[0] != pid && pts[1] != pid && pts[2] != pid)
				continue;

			for (vtkIdType j = 0; j < npts; j ++)
			{
				if (pts[j] == pid)
					continue;

				bool find_ptsj_in_adjofpid = false;
				for (int k = 0; k < num_adjcent[pid]; k ++)
				{
					if (adjcent[pid * 10 + k] == pts[j])
					{
						find_ptsj_in_adjofpid = true;
						break;
					}
				}

				if (find_ptsj_in_adjofpid == false)
				{
					adjcent[pid * 10 + num_adjcent[pid]] = pts[j];
					num_adjcent[pid] ++;
				}
			}
		}
	}
*/
	for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid ++)
	{
		num_adjcent[pid] = 0;
	}
	for(CellId = 0, Strips->InitTraversal(); Strips->GetNextCell(npts, pts); CellId ++)
	{
		if (npts != 3)
		{
			std::cout << "not triangle, smooth cannot work!" << std::endl;
			return 0;
		}
		for (int i = 0; i < npts; i ++)
		{
			int p[2];
			int pidx = 0;
			for (int k = 0; k < npts; k ++)
			{
				if (k != i)
				{
					p[pidx] = k;
					pidx ++;
				}
			}
			for (int l = 0; l < 2; l ++)
			{
				bool find_pl_in_adjofptsi = false;
				for (int k = 0; k < num_adjcent[pts[i]]; k ++)
				{
					if (adjcent[pts[i] * 10 + k] == pts[p[l]])
					{
						find_pl_in_adjofptsi = true;
						break;
					}
				}

				if (find_pl_in_adjofptsi == false)
				{
					adjcent[pts[i] * 10 + num_adjcent[pts[i]]] = pts[p[l]];
					num_adjcent[pts[i]] ++;
				}
			}
		}
	}
	

	// the smooth algorithm
	{
		vtkSmartPointer<vtkPoints> Points_orig = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPoints> Points_last = vtkSmartPointer<vtkPoints>::New();
		Points_orig->DeepCopy(Points);

		double* b = (double*)malloc(Points->GetNumberOfPoints() * 3 * sizeof(double));
		double pi[3], qi[3], oi[3];
		const double alpha = 0.1, beta = 0.2;

		for (int iter = 0; iter < iternum; iter ++)
		{		
			Points_last->DeepCopy(Points);
			for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid ++)
			{
				if (num_adjcent[pid] > 0)
				{
					pi[0] = 0;
					pi[1] = 0;
					pi[2] = 0;
					for (int j = 0; j < num_adjcent[pid]; j ++)
					{
						Points_last->GetPoint(adjcent[pid * 10 + j], qi);
						pi[0] += qi[0];
						pi[1] += qi[1];
						pi[2] += qi[2];
					}
					pi[0] = pi[0] / num_adjcent[pid];
					pi[1] = pi[1] / num_adjcent[pid];
					pi[2] = pi[2] / num_adjcent[pid];
					Points->SetPoint(pid, pi);
					Points_orig->GetPoint(pid, oi);
					Points_last->GetPoint(pid, qi);

					b[pid * 3 + 0] = pi[0] - (alpha * oi[0] + (1.0 - alpha) * qi[0]);
					b[pid * 3 + 1] = pi[1] - (alpha * oi[1] + (1.0 - alpha) * qi[1]);
					b[pid * 3 + 2] = pi[2] - (alpha * oi[2] + (1.0 - alpha) * qi[2]);

				}
			}
			for (vtkIdType pid = 0; pid < Points->GetNumberOfPoints(); pid ++)
			{
				if (num_adjcent[pid] > 0)
				{
					double sumbj[3];
					sumbj[0] = 0;
					sumbj[1] = 0;
					sumbj[2] = 0;

					for (int j = 0; j < num_adjcent[pid]; j ++)
					{					
						sumbj[0] += b[adjcent[pid * 10 + j] * 3 + 0];
						sumbj[1] += b[adjcent[pid * 10 + j] * 3 + 1];
						sumbj[2] += b[adjcent[pid * 10 + j] * 3 + 2];
					}
					sumbj[0] = sumbj[0] / num_adjcent[pid];
					sumbj[1] = sumbj[1] / num_adjcent[pid];
					sumbj[2] = sumbj[2] / num_adjcent[pid];

					Points->GetPoint(pid, pi);

					pi[0] = pi[0] - (beta * b[pid * 3 + 0] + (1.0 - beta) * sumbj[0]);
					pi[1] = pi[1] - (beta * b[pid * 3 + 1] + (1.0 - beta) * sumbj[1]);
					pi[2] = pi[2] - (beta * b[pid * 3 + 2] + (1.0 - beta) * sumbj[2]);

					Points->SetPoint(pid, pi);
				}
			}
		}
		Points_orig->Reset();
		Points_last->Reset();
		free(b);
	}

	Poly->SetPoints(Points);
	return 1;
}









