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

/*
 * Note: this class is not reentrant so thread protection must be explicit in
 * any multi-threaded usage.
 */
class AliHLTMUONCalculations
{
public:

	/* Computes the momentum components based on the equations given in the
	   ALICE dimuon spectrometer Technical Design Report (TDR-5): trigger section.
	   
	   Reference: 
	     "CERN/LHCC 2000-046
	      Addendum 1 to ALICE TDR 5
	      15 Dec 2000"
	     
	     Section 3.1.2 pages 144 and 145.
	     
	   Input can be in meters, cm or mm. 
	   Output is in GeV.
	 */
	static bool ComputeMomentum(
			AliHLTFloat32_t x1,
			AliHLTFloat32_t y1, AliHLTFloat32_t y2,
			AliHLTFloat32_t z1, AliHLTFloat32_t z2
		);
		
	static AliHLTFloat32_t Zf() { return fgZf; }
	static void Zf(AliHLTFloat32_t value) { fgZf = value; }
	static AliHLTFloat32_t QBL();
	static void QBL(AliHLTFloat32_t value);
	
	static AliHLTMUONParticleSign Sign() { return fgSign; }
	static AliHLTFloat32_t Px() { return fgPx; }
	static AliHLTFloat32_t Py() { return fgPy; }
	static AliHLTFloat32_t Pz() { return fgPz; }
		
private:

	// Prevent destroying or creating of this object.
	AliHLTMUONCalculations();
	~AliHLTMUONCalculations();

	static AliHLTFloat32_t fgZf;  // The Z coordinate of the middle of the dipole magnetic field.
	static AliHLTFloat32_t fgQBLScaled; // The integrated field strength times units of charge (T.m.*c/1e9)
	
	static AliHLTMUONParticleSign fgSign;  // The calculated sign.
	static AliHLTFloat32_t fgPx;  // The calculated X momentum (GeV/c).
	static AliHLTFloat32_t fgPy;  // The calculated Y momentum (GeV/c).
	static AliHLTFloat32_t fgPz;  // The calculated Z momentum (GeV/c).
};

#endif // ALIHLTMUONCALCULATIONS_H
