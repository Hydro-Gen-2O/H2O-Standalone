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

#include "gl_helper.h"

#include "point_set.h"

PointSet::PointSet ()
{	
	m_GridRes.Set ( 0, 0, 0 );
	Reset ();
}

void PointSet::Reset ()
{
	// Reset number of particles
//	ResetBuffer ( 0 );	

	m_Time = 0;
	m_DT = 0.1;
	m_Param[POINT_GRAV] = 100.0;
	m_Param[PLANE_GRAV] = 0.0;
	
	m_Vec[POINT_GRAV_POS].Set(0, 0, 50.0);
	m_Vec[PLANE_GRAV_DIR].Set(0, 0, -9.8);
}


void PointSet::Draw ( float* view_mat, float rad )
{
	glEnable ( GL_NORMALIZE );	

	// point drawmode
//glLoadMatrixf(view_mat);
//dat = mBuf.data;
//glBegin(GL_POINTS);
//	for (auto& p : fluidPs) {
//	glColor3f(RED(p->clr), GRN(p->clr), BLUE(p->clr));
//	glVertex3f(p->pos.x, p->pos.y, p->pos.z);
//	dat += mBuf.stride;
//}
//glEnd();
		glLoadMatrixf(view_mat);
		for (auto& f : fluidPs) {
			glPushMatrix();
			glTranslatef(f->pos.x, f->pos.y, f->pos.z);
			glScalef(0.2, 0.2, 0.2);
			glColor4f(RED(f->clr), GRN(f->clr), BLUE(f->clr), ALPH(f->clr));
			drawSphere();
			glPopMatrix();
		}
}

// Ideal grid cell size (gs) = 2 * smoothing radius = 0.02*2 = 0.04
// Ideal domain size = k*gs/d = k*0.02*2/0.005 = k*8 = {8, 16, 24, 32, 40, 48, ..}
//    (k = number of cells, gs = cell size, d = simulation scale)
void PointSet::Grid_Setup ( Vector3DF min, Vector3DF max, float sim_scale, float cell_size, float border )
{
	float world_cellsize = cell_size / sim_scale;
	m_Grid.clear();
	m_GridCnt.clear();
	m_GridMin = min;	m_GridMin -= border;
	m_GridMax = max;	m_GridMax += border;
	m_GridSize = m_GridMax;
	m_GridSize -= m_GridMin;
	m_GridRes.x = ceil ( m_GridSize.x / world_cellsize );		// Determine grid resolution
	m_GridRes.y = ceil ( m_GridSize.y / world_cellsize );
	m_GridRes.z = ceil ( m_GridSize.z / world_cellsize );
	m_GridSize.x = m_GridRes.x * cell_size / sim_scale;				// Adjust grid size to multiple of cell size
	m_GridSize.y = m_GridRes.y * cell_size / sim_scale;
	m_GridSize.z = m_GridRes.z * cell_size / sim_scale;
	m_GridDelta = m_GridRes;		// delta = translate from world space to cell #
	m_GridDelta /= m_GridSize;
	m_GridTotal = (int)(m_GridSize.x * m_GridSize.y * m_GridSize.z);

	m_Grid.reserve ( m_GridTotal );
	m_GridCnt.reserve ( m_GridTotal );	
	for (int n = 0; n < m_GridTotal; n++) {
		m_Grid.push_back(-1);
		m_GridCnt.push_back(0);
	}
}

void PointSet::Grid_InsertParticles ()
{
	for (int n = 0; n < m_GridTotal; n++) {
		m_Grid.at(n) = -1;
		m_GridCnt.at(n) = 0;
	}
	for (int n = 0; n < fluidPs.size(); ++n) {
		std::unique_ptr<Fluid>& p = fluidPs.at(n);
		int gx = (int)((p->pos.x - m_GridMin.x) * m_GridDelta.x);		// Determine grid cell
		int gy = (int)((p->pos.y - m_GridMin.y) * m_GridDelta.y);
		int gz = (int)((p->pos.z - m_GridMin.z) * m_GridDelta.z);
		int gs = (int)((gz * m_GridRes.y + gy) * m_GridRes.x + gx);
		if (gs >= 0 && gs < m_GridTotal) {
			p->next = m_Grid[gs];
			m_Grid[gs] = n;
			m_GridCnt[gs]++;
		}
	}
}

void PointSet::Grid_FindCells ( Vector3DF p, float radius )
{
	Vector3DI sph_min;
	// Compute sphere range
	sph_min.x = (int)((-radius + p.x - m_GridMin.x) * m_GridDelta.x);
	sph_min.y = (int)((-radius + p.y - m_GridMin.y) * m_GridDelta.y);
	sph_min.z = (int)((-radius + p.z - m_GridMin.z) * m_GridDelta.z);
	if ( sph_min.x < 0 ) sph_min.x = 0;
	if ( sph_min.y < 0 ) sph_min.y = 0;
	if ( sph_min.z < 0 ) sph_min.z = 0;

	m_GridCell[0] = (int)((sph_min.z * m_GridRes.y + sph_min.y) * m_GridRes.x + sph_min.x);
	m_GridCell[1] = m_GridCell[0] + 1;
	m_GridCell[2] = (int)(m_GridCell[0] + m_GridRes.x);
	m_GridCell[3] = m_GridCell[2] + 1;

	if ( sph_min.z+1 < m_GridRes.z ) {
		m_GridCell[4] = (int)(m_GridCell[0] + m_GridRes.y*m_GridRes.x);
		m_GridCell[5] = m_GridCell[4] + 1;
		m_GridCell[6] = (int)(m_GridCell[4] + m_GridRes.x);
		m_GridCell[7] = m_GridCell[6] + 1;
	}
	if ( sph_min.x+1 >= m_GridRes.x ) {
		m_GridCell[1] = -1;		m_GridCell[3] = -1;		
		m_GridCell[5] = -1;		m_GridCell[7] = -1;
	}
	if ( sph_min.y+1 >= m_GridRes.y ) {
		m_GridCell[2] = -1;		m_GridCell[3] = -1;
		m_GridCell[6] = -1;		m_GridCell[7] = -1;
	}
}
