/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	inelastic_solid_dynamics.h
 * @brief 	Here, we define the algorithm classes for inelastic_solid dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Xiaojing Tang, Chi Zhang and Xiangyu Hu
 * @version	1.0
 *			Try to implement EIGEN libaary for base vector, matrix and 
 *			linear algebra operation.  
 *			-- Chi ZHANG
 */

#pragma once

#include "elastic_dynamics.h"
#include "inelastic_solid.h"

namespace SPH
{
	namespace solid_dynamics
	{
		/**
		 * @class PlasticStressRelaxationFirstHalf
		 * @brief computing stress relaxation process by verlet time stepping
		 * This is the first step
		 */
		class PlasticStressRelaxationFirstHalf
			: public StressRelaxationFirstHalf
		{
		public:
			PlasticStressRelaxationFirstHalf(BaseInnerRelation &inner_relation);
			virtual ~PlasticStressRelaxationFirstHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);

		protected:
			PlasticSolid &plastic_solid_;
		};
	}
}
