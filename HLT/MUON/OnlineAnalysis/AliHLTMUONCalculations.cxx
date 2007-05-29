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


AliHLTFloat32_t AliHLTMUONCalculateSignedPt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2,
		AliHLTFloat32_t& p
	)
{
	register AliHLTFloat32_t qBL = 3.0;
	register AliHLTFloat32_t zf = 975.0;
	AliHLTFloat32_t thetaTimesZf = - (y1*z2 - y2*z1) / (z2-z1);
	AliHLTFloat32_t xf = x1 * zf / z1;
	AliHLTFloat32_t yf = y2 - ((y2-y1) * (z2-zf)) / (z2-z1);

	AliHLTFloat32_t pDivZf = (qBL * 0.3 / thetaTimesZf);
	p = (AliHLTFloat32_t) fabs( pDivZf * zf );
	
	// Note: the 0.3 is a conversion factor to GeV. it is 1e9 / c where c is
	// the speed of light.
	AliHLTFloat32_t pt = pDivZf * sqrt(xf*xf+yf*yf);
	return pt;
};


AliHLTFloat32_t AliHLTMUONCalculateSignedPt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2,
		register AliHLTFloat32_t zf, register AliHLTFloat32_t qBL,
		AliHLTFloat32_t& p
	)
{
	AliHLTFloat32_t thetaTimesZf = - (y1*z2 - y2*z1) / (z2-z1);
	AliHLTFloat32_t xf = x1 * zf / z1;
	AliHLTFloat32_t yf = y2 - ((y2-y1) * (z2-zf)) / (z2-z1);

	AliHLTFloat32_t pDivZf = (qBL * 0.3 / thetaTimesZf);
	p = (AliHLTFloat32_t) fabs( pDivZf * zf );

	// Note: the 0.3 is a conversion factor to GeV. it is 1e9 / c where c is
	// the speed of light.
	AliHLTFloat32_t pt = pDivZf * sqrt(xf*xf+yf*yf);
	return pt;
};


AliHLTFloat32_t AliHLTMUONCalculateSignedPt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2
	)
{
	register AliHLTFloat32_t qBL = 3.0;
	register AliHLTFloat32_t zf = 975.0;
	AliHLTFloat32_t thetaTimesZf = - (y1*z2 - y2*z1) / (z2-z1);
	AliHLTFloat32_t xf = x1 * zf / z1;
	AliHLTFloat32_t yf = y2 - ((y2-y1) * (z2-zf)) / (z2-z1);

	// Note: the 0.3 is a conversion factor to GeV. it is 1e9 / c where c is
	// the speed of light.
	AliHLTFloat32_t pt = (qBL * 0.3 / thetaTimesZf) * sqrt(xf*xf+yf*yf);
	return pt;
};


AliHLTFloat32_t AliHLTMUONCalculateSignedPt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2,
		register AliHLTFloat32_t zf, register AliHLTFloat32_t qBL
	)
{
	AliHLTFloat32_t thetaTimesZf = - (y1*z2 - y2*z1) / (z2-z1);
	AliHLTFloat32_t xf = x1 * zf / z1;
	AliHLTFloat32_t yf = y2 - ((y2-y1) * (z2-zf)) / (z2-z1);

	// Note: the 0.3 is a conversion factor to GeV. it is 1e9 / c where c is
	// the speed of light.
	AliHLTFloat32_t pt = (qBL * 0.3 / thetaTimesZf) * sqrt(xf*xf+yf*yf);
	return pt;
};


AliHLTFloat32_t AliHLTMUONCalculatePt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2
	)
{
	return (AliHLTFloat32_t) fabs(AliHLTMUONCalculateSignedPt(x1, y1, y2, z1, z2));
};


AliHLTFloat32_t AliHLTMUONCalculatePt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2,
		register AliHLTFloat32_t zf, register AliHLTFloat32_t qBL
	)
{
	return (AliHLTFloat32_t) fabs(AliHLTMUONCalculateSignedPt(x1, y1, y2, z1, z2, zf, qBL));
};

