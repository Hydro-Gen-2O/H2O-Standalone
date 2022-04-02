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

#ifndef DEF_FLUID_SYS
	#define DEF_FLUID_SYS

	#include <iostream>
	#include <vector>
	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	#include "fluid.h"
	
	// Scalar params
	#define SPH_VISC			6

	// Vector params
	#define SPH_VOLMIN			7
	#define SPH_VOLMAX			8
	#define SPH_INITMIN			9
	#define SPH_INITMAX			10
	
	#define MAX_PARAM			21

#define MAX_NEIGHBOR		80

#define MAX_PARAM			21

	class FluidSystem {
	public:
		FluidSystem ();

		void Draw(float* view_mat, float rad);

		virtual void Run ();

		void SPH_CreateExample(int n, int nmax);
		void SPH_DrawDomain();

	private:
		void PredictPositions();
		void FindNeighbors();
		void ComputeDensity();
		void ComputeLambda();
		void ComputeCorrections();
		void ApplyCorrections();
		void Advance();

		glm::vec3 GetGridPos(glm::vec3 pos);
		int GetGridIndex(glm::vec3 gridPos);

		float						m_Param[MAX_PARAM];			// see defines above
		glm::vec3					m_Vec[MAX_PARAM];
		float m_DT;

		std::vector<std::unique_ptr<Fluid>> fluidPs;
		// grid maps indexSpace To vector of fluid there
		std::vector<std::vector<int>> grid;
		std::vector<std::vector<int>> neighbors;

		// the axis between the bounds of the fluid https://en.wikipedia.org/wiki/Space_diagonal
		glm::vec3 gridSpaceDiag;
		int totalGridCells;
		// Smoothed Particle Hydrodynamics
		float m_Poly6Kern, m_SpikyKern; // Kernel functions
	};

#endif
