#ifndef ALIHLTMUONCOREPOINT_H
#define ALIHLTMUONCOREPOINT_H
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
// A 2D point structure using floats.
// These are used to store impact points on the trigger chambers and 
// cluster centroids.
// 
////////////////////////////////////////////////////////////////////////////////

#include "AliHLTMUONBasicTypes.h"

class AliHLTMUONCorePoint
{
public:

	AliHLTMUONCorePoint() : fX(0.0), fY(0.0) {}

	AliHLTMUONCorePoint(Float x, Float y) : fX(x), fY(y) {}

	Float X() const { return fX; }
	Float& X() { return fX; }
	void X(Float x) { fX = x; }
	Float Y() const { return fY; }
	Float& Y() { return fY; }
	void Y(Float y) { fY = y; }

private:

	Float fX, fY; // X and Y coordinate.
};

#endif // ALIHLTMUONCOREPOINT_H
