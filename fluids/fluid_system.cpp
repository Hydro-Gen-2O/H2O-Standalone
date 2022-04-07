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
#include <gl/glut.h>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "common_defs.h"
#include "fluid_system.h"

FluidSystem::FluidSystem()
{
}

void FluidSystem::Draw(float* view_mat)
{
	glEnable(GL_NORMALIZE);
	glLoadMatrixf(view_mat);
	for (auto& f : fluidPs) {
		glm::vec3 scaledPos = f->pos;
		scaledPos /= SPH_RADIUS;
		glPushMatrix();
		glTranslatef(scaledPos.x, scaledPos.y, scaledPos.z);
		glScalef(SPH_RADIUS, SPH_RADIUS, SPH_RADIUS);
		glColor4f(RED(f->clr), GRN(f->clr), BLUE(f->clr), ALPH(f->clr));
		drawSphere();
		glPopMatrix();
	}
}

void FluidSystem::SPH_DrawDomain()
{
	glm::vec3 min, max;
	min = m_Vec[SPH_VOLMIN];
	max = m_Vec[SPH_VOLMAX];
	min.z += 0.5;

	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(min.x, min.y, min.z);	glVertex3f(max.x, min.y, min.z);
	glVertex3f(min.x, max.y, min.z);	glVertex3f(max.x, max.y, min.z);
	glVertex3f(min.x, min.y, min.z);	glVertex3f(min.x, max.y, min.z);
	glVertex3f(max.x, min.y, min.z);	glVertex3f(max.x, max.y, min.z);
	glEnd();
}

void FluidSystem::SPH_CreateExample(int n, int nmax)
{
	m_Poly6Kern = 315.0f / (64.0f * 3.141592 * pow(SPH_RADIUS, 9));	// Wpoly6 kernel (denominator part) - 2003 Muller, p.4

	m_Vec[SPH_VOLMIN] = glm::vec3(-4, -4, 0);
	m_Vec[SPH_VOLMAX] = glm::vec3(4, 4, 10);

	m_Vec[SPH_INITMIN] = glm::vec3(-2, -2, 3);
	m_Vec[SPH_INITMAX] = glm::vec3(2, 2, 7);  
	float ss = 1.f / 2;
	//ss = 1;
	for (float x = m_Vec[SPH_INITMIN].x; x < m_Vec[SPH_INITMAX].x; x += ss) {
		for (float y = m_Vec[SPH_INITMIN].y; y < m_Vec[SPH_INITMAX].y; y += ss) {
			for (float z = m_Vec[SPH_INITMIN].z; z < m_Vec[SPH_INITMAX].z; z += ss) {
				//glm::vec3 init = glm::vec3(x, y, z) * SPH_RADIUS * 0.5f;
				glm::vec3 init = glm::vec3(x, y, z) * SPH_RADIUS;
				fluidPs.push_back(std::make_unique<Fluid>(init, COLORA(0.5, 0.5, 0.6, 1)));
			}
		}
	}

	// SETUP GRID
	gridSpaceDiag = m_Vec[SPH_VOLMAX] - m_Vec[SPH_VOLMIN];
	totalGridCells = (int)(gridSpaceDiag.x * gridSpaceDiag.y * gridSpaceDiag.z);
	grid.reserve(totalGridCells);
	for (int _ = 0; _ < totalGridCells; ++_) {
		std::vector<int>& test = std::vector<int>();
		test.reserve(MAX_NEIGHBOR);
		grid.push_back(test);
	}

	for (int _ = 0; _ < fluidPs.size(); ++_) {
		neighbors.push_back(std::vector<int>());
	}
}

glm::vec3 FluidSystem::GetGridPos(const glm::vec3 &pos) {
	// round - replace w/ other things check??
	//glm::vec3 out = glm::vec3(std::lround(pos.x / SPH_RADIUS - m_Vec[SPH_VOLMIN].x),
	//	std::lround(pos.y / SPH_RADIUS - m_Vec[SPH_VOLMIN].y),
	//	std::lround(pos.z / SPH_RADIUS - m_Vec[SPH_VOLMIN].z));
	// ?????
	//out = pos / SPH_RADIUS - m_Vec[SPH_VOLMIN];
	//std::cout << "pos: " << glm::to_string(pos) << std::endl;
	//std::cout << "sca: " << glm::to_string(pos / SPH_RADIUS) << std::endl;
	//std::cout << "out: " << glm::to_string(out) << std::endl;
	//return out;
	return glm::round(pos / SPH_RADIUS) - m_Vec[SPH_VOLMIN];
}

int FluidSystem::GetGridIndex(const glm::vec3 &gridPos) {
	return gridPos.z * gridSpaceDiag.y * gridSpaceDiag.x + gridPos.y * gridSpaceDiag.x + gridPos.x;
}

void FluidSystem::Run()
{
	PredictPositions();
	FindNeighbors();
	for (int _ = 0; _ < FLUID_ITERS; ++_) {
		ComputeDensity();
		ComputeLambda();
		ComputeCorrections();
		ApplyCorrections();
	}
	Advance();
}

void FluidSystem::PredictPositions() {
	for (std::unique_ptr<Fluid>& p : fluidPs) {
		glm::vec3 deltaVel = GRAVITY * m_DT; // a * dt = change in v

		// apply force to velocity (gravity)
		//p->vel += deltaVel;

		// consider limiting velocity
		//if (length vec > m_Param[SPH_LIMIT]) { p->vel = m_Param[SPH_LIMIT]; }
		p->predictPos = p->pos + (p->vel * m_DT);
	}
}

void FluidSystem::FindNeighbors() {
	for (int i = 0; i < totalGridCells; ++i) {
		grid.at(i).clear();
	}

	//equivalinet of insertgrid / update grid finding the postns within the grid
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);

		glm::vec3 gridPos = GetGridPos(p->predictPos);
		int gIndex = GetGridIndex(gridPos);
		if (0 <= gIndex && gIndex < totalGridCells) {
			grid.at(gIndex).push_back(i); // maybe set a limit? (see MAX_NEIGHBOR)
		}
	}

	// equiv. Finding the Neighbors
	for (int i = 0; i < fluidPs.size(); ++i) {
		neighbor_loop:
		neighbors.at(i).clear(); // clear neighbors from prev. while were at it

		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		glm::vec3 gridPos = GetGridPos(p->predictPos);

		int SEARCH_SIZE = 1;
		// 2x2 neighborhood.
		for (int z = -SEARCH_SIZE; z <= SEARCH_SIZE; z++) {
		for (int y = -SEARCH_SIZE; y <= SEARCH_SIZE; y++) {
		for (int x = -SEARCH_SIZE; x <= SEARCH_SIZE; x++) {
			glm::vec3 n = gridPos + glm::vec3(x, y, z);
			if (n.x >= 0 && n.x < gridSpaceDiag.x &&
				n.y >= 0 && n.y < gridSpaceDiag.y &&
				n.z >= 0 && n.z < gridSpaceDiag.z) {
				int gIndex = GetGridIndex(n);
				for (int pIndex : grid.at(gIndex)) { // each 
					std::unique_ptr<Fluid>& pcurr = fluidPs.at(pIndex);
					float lenR = glm::length(p->predictPos - pcurr->predictPos);
					if (lenR <= SPH_RADIUS && lenR != 0) {
						neighbors.at(i).push_back(pIndex);
						if (neighbors.at(i).size() >= MAX_NEIGHBOR) { goto neighbor_loop; }
					}
				}
			}
		}}}
	}
	// debug print lines
	//std::cout << GetGridIndex(GetGridPos(fluidPs.at(0)->predictPos)) << std::endl;
	//for (int i = 0; i < neighbors.at(0).size(); ++i) {
	//	std::cout << neighbors.at(0).at(i) << " ";
	//}
	//std::cout << std::endl;
}

void FluidSystem::ComputeDensity() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->density = 0.0;
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			// Poly Kernel
			float dist = glm::length(p->predictPos - pcurr->predictPos);
			float c = SPH_RADIUS * SPH_RADIUS - dist * dist;
			p->density += c * c * c * m_Poly6Kern;
		}
	}
}

void FluidSystem::SpikyKernel(glm::vec3 &r) {
	float dist = glm::length(r);
	r *= -45.0f / (3.141592 * pow(SPH_RADIUS, 6));			// grad of spiky (no minus)
	r *= (SPH_RADIUS - dist) * (SPH_RADIUS - dist) / dist;
	r /= REST_DENSITY;
}

void FluidSystem::ComputeLambda() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		float sumGradients = 0.0f;
		glm::vec3 pGrad = glm::vec3(0.f);
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);

			// Spiky Kernel - modifies r by ref
			glm::vec3 r = (p->predictPos - pcurr->predictPos);
			SpikyKernel(r);
			// End Spiky Kernel

			sumGradients += glm::length2(r);
			pGrad += r; // -= r; ?? - i think += b/c -45
		}
		sumGradients += glm::length2(pGrad);
		float constraint = p->density / REST_DENSITY - 1.f; // real scale constraint
		p->lambda = -constraint / (sumGradients + 600.f); // maybe + 500 or so
	}
}

void FluidSystem::ComputeCorrections() {
	float delQ = 0.2f * SPH_RADIUS;
	float cDen = SPH_RADIUS * SPH_RADIUS - delQ * delQ;
	float polyDen = cDen * cDen * cDen * m_Poly6Kern;
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		p->deltaPos = glm::vec3(0.f);
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			//---------Calculate SCORR-----
			float dist = glm::length(p->predictPos - pcurr->predictPos);
			float cNum = SPH_RADIUS * SPH_RADIUS - dist * dist;
			float polyNum = cNum * cNum * cNum * m_Poly6Kern;

			float frac = polyNum / polyDen;
			float sCorr = -K_CORR * frac * frac * frac * frac;
			//------------End SCORR calculation-------

			// Spiky Kernel - modifies r by ref
			glm::vec3 r = (p->predictPos - pcurr->predictPos);
			SpikyKernel(r);
			// End Spiky Kernel

			p->deltaPos += r * (p->lambda + pcurr->lambda + sCorr);
		}
	}
}

void FluidSystem::ApplyCorrections() {
	for (std::unique_ptr<Fluid>& p : fluidPs) {
		p->predictPos += p->deltaPos;
		glm::vec3 scaledMin = m_Vec[SPH_VOLMIN] * SPH_RADIUS;
		glm::vec3 scaledMax = m_Vec[SPH_VOLMAX] * SPH_RADIUS;
		// Perform collision detect and response
		if (p->predictPos.y < scaledMin.y) { p->vel.y = 0.0; p->predictPos.y = scaledMin.y + 0.001; }
		if (p->predictPos.y > scaledMax.y) { p->vel.y = 0.0; p->predictPos.y = scaledMax.y - 0.001; }
																								
		if (p->predictPos.z < scaledMin.z) { p->vel.z = 0.0; p->predictPos.z = scaledMin.z + 0.001; }
		if (p->predictPos.z > scaledMax.z) { p->vel.z = 0.0; p->predictPos.z = scaledMax.z - 0.001; }
																								
		if (p->predictPos.x < scaledMin.x) { p->vel.x = 0.0; p->predictPos.x = scaledMin.x + 0.001; }
		if (p->predictPos.x > scaledMax.x) { p->vel.x = 0.0; p->predictPos.x = scaledMax.x - 0.001; }
	}
}

void FluidSystem::Advance() {
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		//	std::cout << "p pos: " << glm::to_string(p->pos) << std::endl;
		//	std::cout << "p pred: " << glm::to_string(p->predictPos) << std::endl;
		//	std::cout << "p delP: " << glm::to_string(p->deltaPos) << std::endl;
		//	std::cout << "p vel: " << glm::to_string(p->vel) << std::endl;
		//	std::cout << "p dens: " << p->density << std::endl;
		//	std::cout << "p lamb: " << p->lambda << std::endl << std::endl;
		p->vel = (p->predictPos - p->pos) / m_DT;
		p->pos = p->predictPos;
	}
	//std::cout << "before: " << glm::to_string(fluidPs.at(0)->vel) << std::endl;
	// VISCOSITY
	for (int i = 0; i < fluidPs.size(); ++i) {
		std::unique_ptr<Fluid>& p = fluidPs.at(i);
		glm::vec3 acc(0.f, 0.f, 0.f);
		for (int j : neighbors.at(i)) { // for each neighbor
			std::unique_ptr<Fluid>& pcurr = fluidPs.at(j);
			float dist = glm::length(p->predictPos - pcurr->predictPos);
			float c = SPH_RADIUS * SPH_RADIUS - dist * dist;
			acc += (pcurr->vel - p->vel) * c * c * c * m_Poly6Kern;
		}
		p->vel += VISC_CONST * acc * m_DT;
	}
	// END VISCOSITY
	//std::cout << "after: " << glm::to_string(fluidPs.at(0)->vel) << std::endl;
}