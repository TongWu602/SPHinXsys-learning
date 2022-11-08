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
 * @file 	base_data_type.h
 * @brief 	This is the class for Wenland kernel. 
 * @author	Chi ZHang and Xiangyu Hu
 * @version 0.1
 *			NThis kernel has compact support of 2h.
 * 			The smoothing length h can be variable when variable h functions are applied.
 *			Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	1.0
 *			Try to implement EIGEN libaary for base vector, matrix and 
 *			linear algebra operation.  
 *			-- Chi ZHANG
 */

#ifndef KERNEL_WENLAND_C2_H
#define KERNEL_WENLAND_C2_H

#include "base_kernel.h"

namespace SPH
{
	/**
	 * @class KernelWendlandC2
	 * @brief Kernel WendlandC2
	 */
	class KernelWendlandC2 : public Kernel
	{
	public:
		explicit KernelWendlandC2(Real h);

		/** 
		 * Calculates the kernel value for 
		 * the given distance of two particles 
		 */
		virtual Real W_1D(const Real q) const override;
		virtual Real W_2D(const Real q) const override;
		virtual Real W_3D(const Real q) const override;

		virtual Real dW_1D(const Real q) const override;
		virtual Real dW_2D(const Real q) const override;
		virtual Real dW_3D(const Real q) const override;

		virtual Real d2W_1D(const Real q) const override;
		virtual Real d2W_2D(const Real q) const override;
		virtual Real d2W_3D(const Real q) const override;
	};
}
#endif //KERNEL_WENLAND_C2_H