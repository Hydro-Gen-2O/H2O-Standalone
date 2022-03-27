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

	#include "point_set.h"
	#include "fluid.h"
	
	// Scalar params
	#define SPH_SIMSCALE		5
	#define SPH_VISC			6
	#define SPH_RESTDENSITY		7
	#define SPH_PMASS			8
	#define SPH_PRADIUS			9
	#define SPH_PDIST			10
	#define SPH_SMOOTHRADIUS	11
	#define SPH_INTSTIFF		12
	#define SPH_EXTSTIFF		13
	#define SPH_EXTDAMP			14
	#define SPH_LIMIT			15
	#define BOUND_ZMIN_SLOPE	16
	#define FORCE_XMAX_SIN		17
	#define FORCE_XMIN_SIN		18
	#define MAX_FRAC			19
	#define CLR_MODE			20

	// Vector params
	#define SPH_VOLMIN			7
	#define SPH_VOLMAX			8
	#define SPH_INITMIN			9
	#define SPH_INITMAX			10
	
	#define MAX_PARAM			21
	#define BFLUID				2

	class FluidSystem : public PointSet {
	public:
		FluidSystem ();

		// Basic Particle System
		virtual void Run ();
		virtual void Advance ();
		
		// Smoothed Particle Hydrodynamics
		void SPH_Setup(int n);
		void SPH_CreateExample(int n, int nmax);
		void SPH_DrawDomain ();

		void SPH_ComputePressureGrid ();			// O(kn) - spatial grid
		void SPH_ComputeDensity();
		void SPH_ComputeLambda();
		void SPH_ComputeCorrections();
		void SPH_ComputeForceGridNC ();				// O(cn) - neighbor table

		void SPH_ComputeVorticityAndViscosity();

	private:
		// Smoothed Particle Hydrodynamics
		double m_R2, m_Poly6Kern, m_LapKern, m_SpikyKern; // Kernel functions
	};

#endif
