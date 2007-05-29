#ifndef ALIHLTMUONCALCULATIONS_H
#define ALIHLTMUONCALCULATIONS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliHLTMUONDataTypes.h"


/* Computes the Pt (transverse mementum) based on the equations given in the
   ALICE dimuon spectrometer Technical Design Report (TDR-5): trigger section.
   
   Reference: 
     "CERN/LHCC 2000-046
      Addendum 1 to ALICE TDR 5
      15 Dec 2000"
     
     Section 3.1.2 pages 144 and 145.
     
   Input can be in meters, cm or mm. 
   Output is in GeV.
 */
AliHLTFloat32_t AliHLTMUONCalculatePt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2
	);

/* Performs the same calculation as above however alows the zf and qBL
   parameters to be specified.
 */
AliHLTFloat32_t AliHLTMUONCalculatePt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2,
		register AliHLTFloat32_t zf, register AliHLTFloat32_t qBL
	);

/* The same Pt calculation as above however the sign of the result indicates
   the sign of the particle.
 */
AliHLTFloat32_t AliHLTMUONCalculateSignedPt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2
	);

/* Performs the same calculation as above however alows the zf and qBL
   parameters to be specified.
 */
AliHLTFloat32_t AliHLTMUONCalculateSignedPt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2,
		register AliHLTFloat32_t zf, register AliHLTFloat32_t qBL
	);

/* The same Pt calculation as above however the sign of the result indicates
   the sign of the particle. The momentum is also computed and returned.
 */
AliHLTFloat32_t AliHLTMUONCalculateSignedPt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2,
		AliHLTFloat32_t& p
	);

/* Performs the same calculation as above however alows the zf and qBL
   parameters to be specified.
 */
AliHLTFloat32_t AliHLTMUONCalculateSignedPt(
		register AliHLTFloat32_t x1,
		register AliHLTFloat32_t y1, register AliHLTFloat32_t y2,
		register AliHLTFloat32_t z1, register AliHLTFloat32_t z2,
		register AliHLTFloat32_t zf, register AliHLTFloat32_t qBL,
		AliHLTFloat32_t& p
	);

#endif // ALIHLTMUONCALCULATIONS_H
