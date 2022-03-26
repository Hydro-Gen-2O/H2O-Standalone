/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/

#ifndef DEF_POINT_SET
	#define DEF_POINT_SET

	#include <iostream>
	#include <vector>
	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	
	#include "common_defs.h"
	#include "vector.h"
	#include "fluid.h"

	typedef signed int		xref;
	
	#define MAX_NEIGHBOR		80
	
	#define MAX_PARAM			21

	// Scalar params
	#define PNT_DRAWMODE		0
	#define PNT_DRAWSIZE		1		
	#define POINT_GRAV			2
	#define PLANE_GRAV			3

	// Vector params
	#define POINT_GRAV_POS		5	
	#define PLANE_GRAV_DIR		6

	typedef unsigned short		ushort;

	class GeomBuf {
	public:
		ushort stride;
		char* data;
	};

	class PointSet {
	public:
		PointSet ();

		// Point Sets
		virtual void Draw ( float* view_mat, float rad );		
		virtual void Reset ();
		
		// Parameters			
		void SetParam (int p, float v )		{ m_Param[p] = v; }
		void SetParam (int p, int v )		{ m_Param[p] = (float) v; }
		float GetParam ( int p )			{ return (float) m_Param[p]; }

		// Spatial Subdivision
		void Grid_Setup ( Vector3DF min, Vector3DF max, float sim_scale, float cell_size, float border );		
		void Grid_InsertParticles ();	
		void Grid_FindCells ( Vector3DF p, float radius );

	protected:
		GeomBuf mBuf;
		std::vector<std::unique_ptr<Fluid>> fluidPs;

		// Parameters
		double						m_Param [ MAX_PARAM ];			// see defines above
		Vector3DF					m_Vec[MAX_PARAM];
		
		// Particle System
		double						m_DT;
		double						m_Time;

		// Spatial Grid
		std::vector< int >			m_Grid;
		std::vector< int >			m_GridCnt;
		int							m_GridTotal;			// total # cells
		Vector3DF					m_GridMin;				// volume of grid (may not match domain volume exactly)
		Vector3DF					m_GridMax;
		Vector3DF					m_GridRes;				// resolution in each axis
		Vector3DF					m_GridSize;				// physical size in each axis
		Vector3DF					m_GridDelta;
		float						m_GridCellsize;
		int							m_GridCell[27];

		// Neighbor Table
		unsigned short				m_NC[65536];			// neighbor table (600k)
		unsigned short				m_Neighbor[65536][MAX_NEIGHBOR];	
		float						m_NDist[65536][MAX_NEIGHBOR];

		static int m_pcurr;
	};

#endif
