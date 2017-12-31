/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


/// \file
/// \brief Defines TFactor<> and Factor classes which represent factors in probability distributions.

#ifndef __defined_libdai_factor_h
#define __defined_libdai_factor_h

// dfactor.h also includes realfactor.h, because DFactor has asTFactor() function
#include <dai/dfactor.h>
typedef dai::DFactor Factor;

namespace dai {

/// Returns a binary unnormalized single-variable factor \f$ \exp(hx) \f$ where \f$ x = \pm 1 \f$
/** \param x Variable (should be binary)
 *  \param h Field strength
 */
Factor createFactorIsing( const Var &x, Real h );


/// Returns a binary unnormalized pairwise factor \f$ \exp(J x_1 x_2) \f$ where \f$ x_1, x_2 = \pm 1 \f$
/** \param x1 First variable (should be binary)
 *  \param x2 Second variable (should be binary)
 *  \param J Coupling strength
 */
Factor createFactorIsing( const Var &x1, const Var &x2, Real J );


/// Returns a random factor on the variables \a vs with strength \a beta
/** Each entry are set by drawing a normally distributed random with mean
 *  0 and standard-deviation \a beta, and taking its exponent.
 *  \param vs Variables
 *  \param beta Factor strength (inverse temperature)
 */
Factor createFactorExpGauss( const VarSet &vs, Real beta );


/// Returns a pairwise Potts factor \f$ \exp( J \delta_{x_1, x_2} ) \f$
/** \param x1 First variable
 *  \param x2 Second variable (should have the same number of states as \a x1)
 *  \param J  Factor strength
 */
Factor createFactorPotts( const Var &x1, const Var &x2, Real J );


/// Returns a Kronecker delta point mass
/** \param v Variable
 *  \param state The state of \a v that should get value 1
 */
Factor createFactorDelta( const Var &v, size_t state );


/// Returns a Kronecker delta point mass
/** \param vs Set of variables
 *  \param state The state of \a vs that should get value 1
 */
Factor createFactorDelta( const VarSet& vs, size_t state );

}

#endif
