#ifndef ALIHLTMUONCALCULATIONS_H
#define ALIHLTMUONCALCULATIONS_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliHLTMUONDataTypes.h"

extern "C" struct AliHLTMUONTriggerRecordStruct;

/*
 * Note: this class uses static global variables so thread protection must be
 * explicit in any multi-threaded usage. Or the class should be rewritten.
 */
class AliHLTMUONCalculations
{
public:

	/// Calculates the momentum estimate given two track points behind the
	/// dipole magnet and assuming origin is the interaction point.
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
	
	/// Calculates the invariant mass for a pair of particles.
	static AliHLTFloat32_t ComputeMass(
			AliHLTFloat32_t massA,
			AliHLTFloat32_t pxA,
			AliHLTFloat32_t pyA,
			AliHLTFloat32_t pzA,
			AliHLTFloat32_t massB,
			AliHLTFloat32_t pxB,
			AliHLTFloat32_t pyB,
			AliHLTFloat32_t pzB
		);
	
	static bool FitLineToTriggerRecord(const AliHLTMUONTriggerRecordStruct& trigger);
	
	static bool FitLineToTriggerRecord(
			const AliHLTMUONTriggerRecordStruct& trigger,
			const bool hitset[4]
		);
	
	static bool FitLine(
			const AliHLTMUONTriggerRecordStruct& trigger,
			const bool hitset[4]
		);
	
	static AliHLTFloat32_t IdealZ1() { return fgIdealZ1; }
	static void IdealZ1(AliHLTFloat32_t value) { fgIdealZ1 = value; }
	static AliHLTFloat32_t IdealZ2() { return fgIdealZ2; }
	static void IdealZ2(AliHLTFloat32_t value) { fgIdealZ2 = value; }
	
	static AliHLTFloat32_t IdealX1() { return fgIdealX1; }
	static AliHLTFloat32_t IdealY1() { return fgIdealY1; }
	static AliHLTFloat32_t IdealX2() { return fgIdealX2; }
	static AliHLTFloat32_t IdealY2() { return fgIdealY2; }
	
	static bool FitLineToData(
			const AliHLTFloat32_t* x, const AliHLTFloat32_t* y,
			const AliHLTFloat32_t* z, AliHLTUInt32_t n
		);
	
	static bool FitLineToData(
			const AliHLTFloat32_t* x, const AliHLTFloat32_t* z,
			AliHLTUInt32_t n
		);
	
	static AliHLTFloat32_t Mzx() { return fgMzx; }
	static AliHLTFloat32_t Mzy() { return fgMzy; }
	static AliHLTFloat32_t Czx() { return fgCzx; }
	static AliHLTFloat32_t Czy() { return fgCzy; }
	
	static AliHLTFloat32_t ComputeChi2(
			const AliHLTFloat32_t* x, const AliHLTFloat32_t* y,
			const AliHLTFloat32_t* z, AliHLTUInt32_t n
		);
	
	static AliHLTFloat32_t ComputeChi2(
			const AliHLTMUONTriggerRecordStruct& trigger,
			const bool hitset[4]
		);
	
	static AliHLTFloat32_t SigmaX2() { return fgSigmaX2; }
	static void SigmaX2(AliHLTFloat32_t value) { fgSigmaX2 = (value != 0 ? value : 1.); }
	static AliHLTFloat32_t SigmaY2() { return fgSigmaY2; }
	static void SigmaY2(AliHLTFloat32_t value) { fgSigmaY2 = (value != 0 ? value : 1.); }
	
private:

	// Prevent destroying or creating of this object.
	AliHLTMUONCalculations();
	~AliHLTMUONCalculations();

	static AliHLTFloat32_t fgZf;  /// The Z coordinate of the middle of the dipole magnetic field.
	static AliHLTFloat32_t fgQBLScaled; /// The integrated field strength times units of charge (T.m.*c/1e9)
	
	static AliHLTMUONParticleSign fgSign;  /// The calculated sign.
	static AliHLTFloat32_t fgPx;  /// The calculated X momentum (GeV/c).
	static AliHLTFloat32_t fgPy;  /// The calculated Y momentum (GeV/c).
	static AliHLTFloat32_t fgPz;  /// The calculated Z momentum (GeV/c).
	
	static AliHLTFloat32_t fgSigmaX2;  /// The sigma squared value for the variance / uncertainty in X coordinates.
	static AliHLTFloat32_t fgSigmaY2;  /// The sigma squared value for the variance / uncertainty in Y coordinates.
	
	static AliHLTFloat32_t fgMzx;  /// Calculated slope of the line fitted to the ZX plane. (x = fgMzx * z + fgCzx)
	static AliHLTFloat32_t fgMzy;  /// Calculated slope of the line fitted to the ZY plane. (y = fgMzy * z + fgCzy)
	static AliHLTFloat32_t fgCzx;  /// Calculated coefficient of the line fitted to the ZX plane. (x = fgMzx * z + fgCzx)
	static AliHLTFloat32_t fgCzy;  /// Calculated coefficient of the line fitted to the ZY plane. (y = fgMzy * z + fgCzy)
	
	static AliHLTFloat32_t fgIdealX1;  /// Ideal X coordinate of the point on MT1
	static AliHLTFloat32_t fgIdealY1;  /// Ideal Y coordinate of the point on MT1
	static AliHLTFloat32_t fgIdealZ1;  /// Ideal Z coordinate of the point on MT1
	static AliHLTFloat32_t fgIdealX2;  /// Ideal X coordinate of the point on MT2
	static AliHLTFloat32_t fgIdealY2;  /// Ideal Y coordinate of the point on MT2
	static AliHLTFloat32_t fgIdealZ2;  /// Ideal Z coordinate of the point on MT2
};

#endif // ALIHLTMUONCALCULATIONS_H
