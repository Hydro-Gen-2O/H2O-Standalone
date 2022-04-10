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

	#include <vector>
	#include "fluid.h"
	
	// Physical constants
	#define GRAVITY glm::vec3(0, 0, -9.8)
	
	// Tunable(ish) parameters
	#define FLUID_ITERS 1
	#define m_DT 0.0083f
	#define SPH_RADIUS 0.1f
	#define REST_DENSITY 6378.f
	#define MAX_NEIGHBOR 50
	
	#define K_CORR 0.00001f
	#define VISC_CONST 0.01f

	// Vector params
	#define SPH_VOLMIN glm::vec3(-2, -2, 0)//glm::vec3(-4, -4, 0)
	#define SPH_VOLMAX glm::vec3(2, 2, 10)
	#define SPH_INITMIN	glm::vec3(-1, -1, 3)
	#define SPH_INITMAX	glm::vec3(1, 1, 8)

	class FluidSystem {
	public:
		FluidSystem ();

		void Draw(float* view_mat);

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

		float PolyKernel(float dist);
		void SpikyKernel(glm::vec3 &r);

		glm::ivec3 GetGridPos(const glm::vec3 &pos);
		// get index in grid space
		int GetGridIndex(const glm::ivec3 &gridPos);

		std::vector<std::unique_ptr<Fluid>> fluidPs;
		// grid maps indexSpace To vector of fluid there
		std::vector<std::vector<int>> grid;
		std::vector<std::vector<int>> neighbors;

		// the axis between the bounds of the fluid https://en.wikipedia.org/wiki/Space_diagonal
		glm::ivec3 gridSpaceDiag;
		int totalGridCells;
	};
#endif
