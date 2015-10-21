/*!
 * \file dual_grid_structure.cpp
 * \brief Main classes for defining the dual grid
 * \author F. Palacios, T. Economon
 * \version 4.0.1 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/dual_grid_structure.hpp"

unsigned short CDualGrid::nDim = 0;

CDualGrid::CDualGrid(unsigned short val_nDim) { nDim = val_nDim;}

CDualGrid::~CDualGrid() {}

CPoint::CPoint(unsigned short val_nDim, unsigned long val_globalindex, CConfig *config) : CDualGrid(val_nDim) {
	unsigned short iDim, jDim;
	
	/*--- Element, point and edge structures initialization ---*/
	Elem.clear(); nElem = 0;
	Point.clear(); nPoint = 0;
	Edge.clear();
  
  Volume = NULL;  vertex = NULL;
	coord = NULL; Coord_old = NULL; Coord_sum = NULL;
	Coord_n = NULL; Coord_n1 = NULL;  Coord_p1 = NULL;
	GridVel = NULL; GridVel_Grad = NULL;

	/*--- Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume ---*/
	if (config->GetUnsteady_Simulation() == NO) { Volume = new su2double[1]; Volume[0] = 0.0; }
	else { Volume = new su2double[3]; Volume[0] = 0.0; Volume[1] = 0.0; Volume[2] = 0.0; }
	coord = new su2double[nDim];

	/*--- Indicator if the control volume has been agglomerated ---*/
	Agglomerate = false;
	
  /*--- Flip the normal orientation ---*/
  Flip_Orientation = false;
  
	/*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
	Move = true;

	/*--- Identify boundaries, physical boundaries (not send-receive 
	 condition), detect if an element belong to the domain or it must 
	 be computed with other processor  ---*/
	Boundary = false;
	PhysicalBoundary = false;
	SolidBoundary = false;
	Domain = true;

  /*--- Set the global index in the parallel simulation ---*/
	GlobalIndex = val_globalindex;
  
	/*--- Set the color for mesh partitioning ---*/
	color = 0;

	/*--- For smoothing the numerical grid coordinates ---*/
	if (config->GetSmoothNumGrid()) {
		Coord_old = new su2double[nDim];
		Coord_sum = new su2double[nDim];
	}
	
	/*--- Storage of grid velocities for dynamic meshes ---*/
	if (config->GetGrid_Movement()) {
		GridVel  = new su2double[nDim];
			for (iDim = 0; iDim < nDim; iDim ++) 
		GridVel[iDim] = 0.0;
    
    /*--- Gradient of the grid velocity ---*/
    GridVel_Grad = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      GridVel_Grad[iDim] = new su2double[nDim];
      for (jDim = 0; jDim < nDim; jDim++)
        GridVel_Grad[iDim][jDim] = 0.0;
    }
    
    /*--- Structures for storing old node coordinates for computing grid 
     velocities via finite differencing with dynamically deforming meshes. ---*/
    if (config->GetUnsteady_Simulation() != NO) {
      Coord_p1 = new su2double[nDim];
      Coord_n  = new su2double[nDim];
      Coord_n1 = new su2double[nDim];
    }
	}
  
  /*--- Intialize the value of the curvature ---*/
  Curvature = 0.0;

}

CPoint::CPoint(su2double val_coord_0, su2double val_coord_1, unsigned long val_globalindex, CConfig *config) : CDualGrid(2) {
	unsigned short iDim, jDim;

	/*--- Element, point and edge structures initialization ---*/
	Elem.clear(); nElem = 0;
	Point.clear(); nPoint = 0;
	Edge.clear();
  
  Volume = NULL;  vertex = NULL;
	coord = NULL; Coord_old = NULL; Coord_sum = NULL;
	Coord_n = NULL; Coord_n1 = NULL;  Coord_p1 = NULL;
	GridVel = NULL; GridVel_Grad = NULL;

	/*--- Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume ---*/
	if (config->GetUnsteady_Simulation() == NO) { Volume = new su2double[1]; Volume[0] = 0.0; }
	else { Volume = new su2double[3]; Volume[0] = 0.0; Volume[1] = 0.0; Volume[2] = 0.0; }
	coord = new su2double[nDim]; coord[0] = val_coord_0; coord[1] = val_coord_1;
	
	/*--- Indicator if the control volume has been agglomerated ---*/
	Agglomerate = false;
	
  /*--- Flip the normal orientation ---*/
  Flip_Orientation = false;
  
	/*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
	Move = true;
	
	/*--- Identify boundaries, physical boundaries (not send-receive 
	 condition), detect if an element belong to the domain or it must 
	 be computed with other processor  ---*/
	Boundary = false;
  PhysicalBoundary = false;
  SolidBoundary = false;
	Domain = true;
	
	/*--- Set the color for mesh partitioning ---*/
	color = 0;
	
	/*--- Set the global index in the parallel simulation ---*/
	GlobalIndex = val_globalindex;
	
	/*--- For smoothing the numerical grid coordinates ---*/
	if (config->GetSmoothNumGrid()) {
		Coord_old = new su2double[nDim];
		Coord_sum = new su2double[nDim];
	}
	
	/*--- Storage of grid velocities for dynamic meshes ---*/
	if (config->GetGrid_Movement()) {
		GridVel  = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim ++)
      GridVel[iDim] = 0.0;
    
    /*--- Gradient of the grid velocity ---*/
    GridVel_Grad = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      GridVel_Grad[iDim] = new su2double[nDim];
      for (jDim = 0; jDim < nDim; jDim++)
        GridVel_Grad[iDim][jDim] = 0.0;
    }
    
    /*--- Structures for storing old node coordinates for computing grid
     velocities via finite differencing with dynamically deforming meshes. ---*/
    if (config->GetUnsteady_Simulation() != NO) {
      Coord_p1 = new su2double[nDim];
      Coord_n  = new su2double[nDim];
      Coord_n1 = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim ++) {
        Coord_p1[iDim] = coord[iDim];
        Coord_n[iDim]  = coord[iDim];
        Coord_n1[iDim] = coord[iDim];
      }
    }
	}
  
  /*--- Intialize the value of the curvature ---*/
  Curvature = 0.0;
  
}

CPoint::CPoint(su2double val_coord_0, su2double val_coord_1, su2double val_coord_2, unsigned long val_globalindex, CConfig *config) : CDualGrid(3) {
	unsigned short iDim, jDim;

	/*--- Element, point and edge structures initialization ---*/
	Elem.clear(); nElem = 0;
	Point.clear(); nPoint = 0;
	Edge.clear();
  
	Volume = NULL;  vertex = NULL;
	coord = NULL; Coord_old = NULL; Coord_sum = NULL;
	Coord_n = NULL; Coord_n1 = NULL;  Coord_p1 = NULL;
	GridVel = NULL; GridVel_Grad = NULL;
  
	/*--- Volume (0 -> Vol_nP1, 1-> Vol_n, 2 -> Vol_nM1 ) and coordinates of the control volume ---*/
	if (config->GetUnsteady_Simulation() == NO) { Volume = new su2double[1]; Volume[0] = 0.0; }
	else { Volume = new su2double[3]; Volume[0] = 0.0; Volume[1] = 0.0; Volume[2] = 0.0; }
	coord = new su2double[nDim]; coord[0] = val_coord_0; coord[1] = val_coord_1; coord[2] = val_coord_2;

	/*--- Indicator if the control volume has been agglomerated ---*/
	Agglomerate = false;
	
	/*--- Indicator if the point is going to be moved in a volumetric deformation ---*/
	Move = true;
	
  /*--- Flip the normal orientation ---*/
  Flip_Orientation = false;

	/*--- Identify boundaries, physical boundaries (not send-receive 
	 condition), detect if an element belong to the domain or it must 
	 be computed with other processor  ---*/
	Boundary = false;
  PhysicalBoundary = false;
  SolidBoundary = false;
	Domain = true;
	
	/*--- Set the color for mesh partitioning ---*/
	color = 0;
	
	/*--- Set the global index in the parallel simulation ---*/
	GlobalIndex = val_globalindex;
	
	/*--- For smoothing the numerical grid coordinates ---*/
	if (config->GetSmoothNumGrid()) {
		Coord_old = new su2double[nDim];
		Coord_sum = new su2double[nDim];
	}
	
	/*--- Storage of grid velocities for dynamic meshes ---*/
	if (config->GetGrid_Movement()) {
		GridVel = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim ++)
      GridVel[iDim] = 0.0;
    
    /*--- Gradient of the grid velocity ---*/
    GridVel_Grad = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      GridVel_Grad[iDim] = new su2double[nDim];
      for (jDim = 0; jDim < nDim; jDim++)
        GridVel_Grad[iDim][jDim] = 0.0;
    }
    
    /*--- Structures for storing old node coordinates for computing grid
     velocities via finite differencing with dynamically deforming meshes. ---*/
    if (config->GetUnsteady_Simulation() != NO) {
      Coord_p1 = new su2double[nDim];
      Coord_n  = new su2double[nDim];
      Coord_n1 = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim ++) {
        Coord_p1[iDim] = coord[iDim];
        Coord_n[iDim]  = coord[iDim];
        Coord_n1[iDim] = coord[iDim];
      }
    }
	}
  
  /*--- Intialize the value of the curvature ---*/
  Curvature = 0.0;
  
}

CPoint::~CPoint() {
  
	Elem.~vector();
	Point.~vector();
	Edge.~vector();
  Children_CV.~vector();

	if (Volume != NULL) delete[] Volume;
	if (vertex != NULL) delete[] vertex;
	if (coord != NULL) delete[] coord;
	if (Coord_old != NULL) delete[] Coord_old;
	if (Coord_sum != NULL) delete[] Coord_sum;
	if (Coord_n != NULL) delete[] Coord_n;
	if (Coord_n1 != NULL) delete[] Coord_n1;
	if (Coord_p1 != NULL) delete[] Coord_p1;
	if (GridVel != NULL) delete[] GridVel;
  if (GridVel_Grad != NULL) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      delete [] GridVel_Grad[iDim];
    delete [] GridVel_Grad;
  }
  
}

void CPoint::SetPoint(unsigned long val_point) {
	unsigned short iPoint;
	bool new_point;
	
	/*--- Look for the point in the list ---*/
	new_point = true;
	for (iPoint = 0; iPoint < GetnPoint(); iPoint++)
		if (Point[iPoint] == val_point) {
			new_point = false; 
			break;
		}

	/*--- Store the point structure and dimensionalizate edge structure ---*/
	if (new_point) {
		Point.push_back(val_point);
		Edge.push_back(-1);
		nPoint = Point.size();
	}
}

void CPoint::SetBoundary(unsigned short val_nmarker) {
	unsigned short imarker;
	
	/*--- To be sure that we are not goint to initializate twice the same vertex ---*/
	if (!Boundary) {
		vertex = new long[val_nmarker];
		/*--- The initialization is made with -1 ---*/
		for (imarker = 0; imarker < val_nmarker; imarker++) 
			vertex[imarker] = -1;
	}
	Boundary = true;
}

CEdge::CEdge(unsigned long val_iPoint, unsigned long val_jPoint, unsigned short val_nDim) : CDualGrid(val_nDim) {
	unsigned short iDim;
	
  /*--- Pointers initialization ---*/
  Coord_CG = NULL;
	Normal = NULL;
	Nodes = NULL;
  
	/*--- Allocate center of gravity coordinates, nodes, and face normal ---*/
	Coord_CG = new su2double[nDim];
	Nodes = new unsigned long[2];
	Normal = new su2double [nDim];

	/*--- Initializate the structure ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_CG[iDim] = 0.0;
		Normal[iDim] = 0.0;
	}
  
	Nodes[0] = val_iPoint; 
	Nodes[1] = val_jPoint;

}

CEdge::~CEdge() {
  
	if (Coord_CG != NULL) delete[] Coord_CG;
	if (Normal != NULL) delete[] Normal;
	if (Nodes != NULL) delete[] Nodes;
  
}

void CEdge::SetCG(su2double **val_coord) {
	unsigned short iDim, iNode;
	
	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_CG[iDim] = 0.0;
		for (iNode = 0; iNode < 2;  iNode++)
			Coord_CG[iDim] += val_coord[iNode][iDim]/2.0;
	}
}

su2double CEdge::GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point) {
	unsigned short iDim;
  su2double vec_a[3] = {0.0,0.0,0.0}, vec_b[3] = {0.0,0.0,0.0}, vec_c[3] = {0.0,0.0,0.0}, vec_d[3] = {0.0,0.0,0.0}, Local_Volume;

	for (iDim = 0; iDim < nDim; iDim++) {
		vec_a[iDim] = val_coord_Edge_CG[iDim]-val_coord_Point[iDim];
		vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Point[iDim];
		vec_c[iDim] = val_coord_Elem_CG[iDim]-val_coord_Point[iDim];
	}

	vec_d[0] = vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1];
	vec_d[1] = -(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0]);
	vec_d[2] = vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0];

	Local_Volume = fabs(vec_c[0]*vec_d[0] + vec_c[1]*vec_d[1] + vec_c[2]*vec_d[2])/6.0;
	
	return Local_Volume;
}

su2double CEdge::GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point) {
	unsigned short iDim;
	su2double vec_a[2] = {0.0,0.0}, vec_b[2] = {0.0,0.0}, Local_Volume;

	for (iDim = 0; iDim < nDim; iDim++) {
		vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Point[iDim];
		vec_b[iDim] = val_coord_Edge_CG[iDim]-val_coord_Point[iDim];
	}

	Local_Volume = 0.5*fabs(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0]);
	
	return Local_Volume;
}

void CEdge::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG) {
	unsigned short iDim;
	su2double vec_a[3] = {0.0,0.0,0.0}, vec_b[3] = {0.0,0.0,0.0}, Dim_Normal[3];

	for (iDim = 0; iDim < nDim; iDim++) {
		vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Edge_CG[iDim];
		vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Edge_CG[iDim];
	}

	Dim_Normal[0] = 0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1]);
	Dim_Normal[1] = -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0]);
	Dim_Normal[2] = 0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0]);
	
	Normal[0] += Dim_Normal[0]; 
	Normal[1] += Dim_Normal[1];		
	Normal[2] += Dim_Normal[2];
  
}

void CEdge::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG) {
	su2double Dim_Normal[2];

	Dim_Normal[0] = val_coord_Elem_CG[1]-val_coord_Edge_CG[1];
	Dim_Normal[1] = -(val_coord_Elem_CG[0]-val_coord_Edge_CG[0]);
	
	Normal[0] += Dim_Normal[0]; 
	Normal[1] += Dim_Normal[1];
  
}

CVertex::CVertex(unsigned long val_point, unsigned short val_nDim) : CDualGrid(val_nDim) {
	unsigned short iDim;
	
  /*--- Pointers initialization ---*/
  Nodes = NULL;
	Normal = NULL;
  
	/*--- Allocate node, and face normal ---*/
	Nodes = new unsigned long[1]; 
	Normal = new su2double [nDim];

	/*--- Initializate the structure ---*/
	Nodes[0] = val_point;
	for (iDim = 0; iDim < nDim; iDim ++) Normal[iDim] = 0.0;
	
	/*--- Set to zero the variation of the coordinates ---*/
	VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;

}

CVertex::~CVertex() {
  
	if (Normal != NULL) delete[] Normal;
	if (Nodes != NULL) delete[] Nodes;
  
}

void CVertex::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG) {
  su2double vec_a[3] = {0.0,0.0,0.0}, vec_b[3] = {0.0,0.0,0.0}, Dim_Normal[3] = {0.0,0.0,0.0};
	unsigned short iDim;

	for (iDim = 0; iDim < nDim; iDim++) {
		vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Edge_CG[iDim];
		vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Edge_CG[iDim];
	}

	Dim_Normal[0] = 0.5*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1]);
	Dim_Normal[1] = -0.5*(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0]);
	Dim_Normal[2] = 0.5*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0]);
		
	Normal[0] += Dim_Normal[0]; 
	Normal[1] += Dim_Normal[1];	
	Normal[2] += Dim_Normal[2];
  
}

void CVertex::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG) {
	su2double Dim_Normal[2];

	Dim_Normal[0] = val_coord_Elem_CG[1]-val_coord_Edge_CG[1];
	Dim_Normal[1] = -(val_coord_Elem_CG[0]-val_coord_Edge_CG[0]);

	Normal[0] += Dim_Normal[0]; 
	Normal[1] += Dim_Normal[1];
  
}

void CVertex::AddNormal(su2double *val_face_normal) {

	Normal[0] += val_face_normal[0]; 
	Normal[1] += val_face_normal[1];
	if (nDim == 3) Normal[2] += val_face_normal[2];
}
