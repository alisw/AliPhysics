/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliHLTMUONCalculations.h"
#include <cmath>

AliHLTFloat32_t AliHLTMUONCalculations::fgZf = -975.0;  // cm

AliHLTFloat32_t AliHLTMUONCalculations::fgQBLScaled
	= 3.0 * 2.99792458e8 / 1e9; // T.m.*c/1e9
	
AliHLTMUONParticleSign AliHLTMUONCalculations::fgSign = kSignUnknown;
AliHLTFloat32_t AliHLTMUONCalculations::fgPx = 0;
AliHLTFloat32_t AliHLTMUONCalculations::fgPy = 0;
AliHLTFloat32_t AliHLTMUONCalculations::fgPz = 0;


bool AliHLTMUONCalculations::ComputeMomentum(
		AliHLTFloat32_t x1,
		AliHLTFloat32_t y1, AliHLTFloat32_t y2,
		AliHLTFloat32_t z1, AliHLTFloat32_t z2
	)
{
	AliHLTFloat64_t z2mz1 = z2 - z1;
	if (z2mz1 == 0 or z1 == 0)
	{
		fgSign = kSignUnknown;
		fgPx = fgPy = fgPz = 0;
		return false;
	}
	AliHLTFloat64_t thetaTimesZf = (y1*z2 - y2*z1) / z2mz1;
	AliHLTFloat64_t xf = x1 * fgZf / z1;
	AliHLTFloat64_t yf = y2 - ((y2-y1) * (z2-fgZf)) / z2mz1;

	if (thetaTimesZf == 0)
	{
		fgSign = kSignUnknown;
		fgPx = fgPy = fgPz = 0;
		return false;
	}
	AliHLTFloat64_t pDivZf = (fgQBLScaled / thetaTimesZf);
	AliHLTFloat64_t p = pDivZf * fgZf;
	pDivZf = fabs(pDivZf);
	
	if (p < 0)
		fgSign = kSignMinus;
	else if (p > 0)
		fgSign = kSignPlus;
	else
		fgSign = kSignUnknown;
	
	fgPx = AliHLTFloat32_t( pDivZf * xf );
	fgPy = AliHLTFloat32_t( pDivZf * yf );
	fgPz = AliHLTFloat32_t( sqrt(p*p - fgPx*fgPx - fgPy*fgPy) );
	// fgPz must be the same sign as fgZf else it could not have been measured.
	if (fgZf < 0) fgPz = -fgPz;

	return true;
}


AliHLTFloat32_t AliHLTMUONCalculations::QBL()
{
	// We have to convert back into Tesla metres.
	return fgQBLScaled * 1e9 / 2.99792458e8;
}


void AliHLTMUONCalculations::QBL(AliHLTFloat32_t value)
{
	// Note: 2.99792458e8/1e9 is the conversion factor for GeV.
	// It is c/1e9, where c is the speed of light.
	fgQBLScaled = value * 2.99792458e8 / 1e9;
}
