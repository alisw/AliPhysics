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

// $Id$

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliHLTMUONCalculations.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include <cmath>

AliHLTFloat32_t AliHLTMUONCalculations::fgZf = -975.0;  // cm

AliHLTFloat32_t AliHLTMUONCalculations::fgQBLScaled
	= 3.0 * 2.99792458e8 / 1e9; // T.m.*c/1e9
	
AliHLTMUONParticleSign AliHLTMUONCalculations::fgSign = kSignUnknown;
AliHLTFloat32_t AliHLTMUONCalculations::fgPx = 0;  // GeV/c
AliHLTFloat32_t AliHLTMUONCalculations::fgPy = 0;  // GeV/c
AliHLTFloat32_t AliHLTMUONCalculations::fgPz = 0;  // GeV/c

AliHLTFloat32_t AliHLTMUONCalculations::fgSigmaX2 = 1.;  // cm^2
AliHLTFloat32_t AliHLTMUONCalculations::fgSigmaY2 = 1.;  // cm^2

AliHLTFloat32_t AliHLTMUONCalculations::fgMzx = 0;
AliHLTFloat32_t AliHLTMUONCalculations::fgMzy = 0;
AliHLTFloat32_t AliHLTMUONCalculations::fgCzx = 0;
AliHLTFloat32_t AliHLTMUONCalculations::fgCzy = 0;

AliHLTFloat32_t AliHLTMUONCalculations::fgIdealX1 = 0;  // cm
AliHLTFloat32_t AliHLTMUONCalculations::fgIdealY1 = 0;  // cm
AliHLTFloat32_t AliHLTMUONCalculations::fgIdealZ1 = -1603.5f;  // cm
AliHLTFloat32_t AliHLTMUONCalculations::fgIdealX2 = 0;  // cm
AliHLTFloat32_t AliHLTMUONCalculations::fgIdealY2 = 0;  // cm
AliHLTFloat32_t AliHLTMUONCalculations::fgIdealZ2 = -1703.5f;  // cm


bool AliHLTMUONCalculations::ComputeMomentum(
		AliHLTFloat32_t x1,
		AliHLTFloat32_t y1, AliHLTFloat32_t y2,
		AliHLTFloat32_t z1, AliHLTFloat32_t z2
	)
{
	/// Computes the momentum components based on the equations given in the
	///   ALICE dimuon spectrometer Technical Design Report (TDR-5): trigger section.
	///
	///   Reference: 
	///     "CERN/LHCC 2000-046
	///      Addendum 1 to ALICE TDR 5
	///      15 Dec 2000"
	///     Section 3.1.2 pages 144 and 145.
	///
	/// Input can be in meters, cm or mm. Output is in GeV/c.
	///
	/// \param x1  X coordinate of hit point 1 on the track.
	/// \param y1  Y coordinate of hit point 1 on the track.
	/// \param z1  Z coordinate of hit point 1 on the track.
	/// \param y2  Y coordinate of hit point 2 on the track.
	/// \param z2  Z coordinate of hit point 2 on the track.
	/// \return true if the momentum could be calculated and false otherwise.
	///    If true is returned then the estimated momentum can be fetched by the
	///    method calls: Px(), Py() and Pz() for the individual components.
	
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
	AliHLTFloat64_t k = p*p - fgPx*fgPx - fgPy*fgPy;
	if (k > 0)
		fgPz = AliHLTFloat32_t( sqrt(k) );
	else
		fgPz = 0;
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


AliHLTFloat32_t AliHLTMUONCalculations::ComputeMass(
		AliHLTFloat32_t massA,
		AliHLTFloat32_t pxA,
		AliHLTFloat32_t pyA,
		AliHLTFloat32_t pzA,
		AliHLTFloat32_t massB,
		AliHLTFloat32_t pxB,
		AliHLTFloat32_t pyB,
		AliHLTFloat32_t pzB
	)
{
	/// Calculates the invariant mass for a pair of particles.
	/// \param massA Mmass in GeV/c of particle A.
	/// \param pxA  X component of the momentum in GeV/c for particle A.
	/// \param pyA  Y component of the momentum in GeV/c for particle A.
	/// \param pzA  Z component of the momentum in GeV/c for particle A.
	/// \param massB  Mass in GeV/c of particle B.
	/// \param pxB  X component of the momentum in GeV/c for particle B.
	/// \param pyB  Y component of the momentum in GeV/c for particle B.
	/// \param pzB  Z component of the momentum in GeV/c for particle B.
	/// \return  The invariant mass in GeV/c^2 or -1 if there was a problem
	///          in the calculation due to bad input parameters.
	
	AliHLTFloat32_t massA2 = massA*massA;
	AliHLTFloat32_t massB2 = massB*massB;
	AliHLTFloat32_t energyA = sqrt(massA2 + pxA*pxA + pyA*pyA + pzA*pzA);
	AliHLTFloat32_t energyB = sqrt(massB2 + pxB*pxB + pyB*pyB + pzB*pzB);
	AliHLTFloat32_t mass2 = massA2 + massB2 + 2. * (energyA*energyB - pxA*pxB - pyA*pyB - pzA*pzB);
	if (mass2 < 0.) return -1.;
	return sqrt(mass2);
}


bool AliHLTMUONCalculations::FitLineToTriggerRecord(
		const AliHLTMUONTriggerRecordStruct& trigger
	)
{
	/// Straight lines are fitted in the ZX and ZY planes using a least
	/// squares fit to the coordinates in the trigger record.
	/// http://mathworld.wolfram.com/LeastSquaresFitting.html
	/// If this method returns true, then the fitted parameters can fetched
	/// using the method calls Mzx(), Mzy(), Czx() and Czy(). The lines are
	/// then given by: x = Mzx() * z + Czx() and y = Mzy() * z + Czy()
	/// The ideal coordinates are also calculated and can be fetched with
	/// the method calls: IdealX1(), IdealY1() and IdealZ1() for point on MT1,
	/// and IdealX2(), IdealY2() and IdealZ2() for point on MT2.
	/// \param trigger  The trigger record structure to which we fit a line.
	/// \return  true if the line could be fitted or false otherwise.
	///     The reason for failure could be either too few hits or the slopes
	///     Mzx() or Mzy() would be infinite, implying a line that is
	///     perpendicular to the z axis.
	
	AliHLTMUONParticleSign sign;
	bool hitset[4];
	AliHLTMUONUtils::UnpackTriggerRecordFlags(trigger.fFlags, sign, hitset);
	DebugTrace("hitset = {" << hitset[0] << ", " << hitset[1] << ", "
		<< hitset[2] << ", " << hitset[3] << "}"
	);
	
	return FitLineToTriggerRecord(trigger, hitset);
}


bool AliHLTMUONCalculations::FitLineToTriggerRecord(
		const AliHLTMUONTriggerRecordStruct& trigger,
		const bool hitset[4]
	)
{
	/// Performs a straight line fit like FitLineToTriggerRecord(trigger)
	/// but requires pree-decoded flags indicating which hits were set.
	/// \param trigger  The trigger record structure to which we fit a line.
	/// \param hitset  Flags indicating which hits were set in the trigger record.
	/// \return  true if the line could be fitted or false otherwise.
	
	bool lineOk = FitLine(trigger, hitset);
	if (lineOk)
	{
		// Calculate ideal points on chambers 11 and 13:
		fgIdealX1 = fgMzx * fgIdealZ1 + fgCzx;
		fgIdealY1 = fgMzy * fgIdealZ1 + fgCzy;
		fgIdealX2 = fgMzx * fgIdealZ2 + fgCzx;
		fgIdealY2 = fgMzy * fgIdealZ2 + fgCzy;
	}
	return lineOk;
}


bool AliHLTMUONCalculations::FitLine(
		const AliHLTMUONTriggerRecordStruct& trigger,
		const bool hitset[4]
	)
{
	/// Performs a straight line fit to the trigger record hits which are indicated
	/// by the hitset flags array.
	/// \param trigger  The trigger record structure to which we fit a line.
	/// \param hitset  Flags indicating which hits to use and were set in the trigger record.
	/// \return  true if the line could be fitted or false otherwise.
	
	AliHLTFloat32_t sumX = 0;
	AliHLTFloat32_t sumY = 0;
	AliHLTFloat32_t sumZ = 0;
	int n = 0;
	for (int i = 0; i < 4; i++)
	{
		if (hitset[i])
		{
			sumX += trigger.fHit[i].fX;
			sumY += trigger.fHit[i].fY;
			sumZ += trigger.fHit[i].fZ;
			n++;
		}
	}
	if (n < 2) return false;
	AliHLTFloat32_t meanX = sumX / AliHLTFloat32_t(n);
	AliHLTFloat32_t meanY = sumY / AliHLTFloat32_t(n);
	AliHLTFloat32_t meanZ = sumZ / AliHLTFloat32_t(n);
	
	AliHLTFloat32_t vSSzz = 0;
	AliHLTFloat32_t vSSzx = 0;
	AliHLTFloat32_t vSSzy = 0;
	for (int i = 0; i < 4; i++)
	{
		if (hitset[i])
		{
			vSSzz += (trigger.fHit[i].fZ - meanZ)*(trigger.fHit[i].fZ - meanZ);
			vSSzx += (trigger.fHit[i].fZ - meanZ)*(trigger.fHit[i].fX - meanX);
			vSSzy += (trigger.fHit[i].fZ - meanZ)*(trigger.fHit[i].fY - meanY);
		}
	}
	
	// Calculate params for lines x = fgMzx * z + fgCzx and y = fgMzy * z + fgCzy.
	if (vSSzz == 0) return false;
	fgMzx = vSSzx / vSSzz;
	fgMzy = vSSzy / vSSzz;
	fgCzx = meanX - fgMzx * meanZ;
	fgCzy = meanY - fgMzy * meanZ;
	
	return true;
}


bool AliHLTMUONCalculations::FitLineToData(
		const AliHLTFloat32_t* x, const AliHLTFloat32_t* y,
		const AliHLTFloat32_t* z, AliHLTUInt32_t n
	)
{
	/// Straight lines are fitted in the ZX and ZY planes using a least
	/// squares fit for the (x, y, z) data points.
	/// http://mathworld.wolfram.com/LeastSquaresFitting.html
	/// If this method returns true, then the fitted parameters can fetched
	/// using the method calls Mzx(), Mzy(), Czx() and Czy(). The lines are
	/// then given by: x = Mzx() * z + Czx() and y = Mzy() * z + Czy()
	/// \param x  This must point to the array of x data values.
	/// \param y  This must point to the array of y data values.
	/// \param z  This must point to the array of z data values.
	/// \param n  Specifies the number of data points in the x, y and z arrays.
	/// \return  true if the line could be fitted or false otherwise.
	///     The reason for failure could be either too few data points or the
	///     slopes Mzx() or Mzy() would be infinite, implying a line that is
	///     perpendicular to the z axis.
	
	if (n < 2) return false;
	
	AliHLTFloat32_t sumX = 0;
	AliHLTFloat32_t sumY = 0;
	AliHLTFloat32_t sumZ = 0;
	for (AliHLTUInt32_t i = 0; i < n; i++)
	{
		sumX += x[i];
		sumY += y[i];
		sumZ += z[i];
	}
	AliHLTFloat32_t meanX = sumX / AliHLTFloat32_t(n);
	AliHLTFloat32_t meanY = sumY / AliHLTFloat32_t(n);
	AliHLTFloat32_t meanZ = sumZ / AliHLTFloat32_t(n);
	
	AliHLTFloat32_t vSSzz = 0;
	AliHLTFloat32_t vSSzx = 0;
	AliHLTFloat32_t vSSzy = 0;
	for (AliHLTUInt32_t i = 0; i < n; i++)
	{
		vSSzz += (z[i] - meanZ)*(z[i] - meanZ);
		vSSzx += (z[i] - meanZ)*(x[i] - meanX);
		vSSzy += (z[i] - meanZ)*(y[i] - meanY);
	}
	
	// Calculate params for lines x = fgMzx * z + fgCzx and y = fgMzy * z + fgCzy.
	if (vSSzz == 0) return false;
	fgMzx = vSSzx / vSSzz;
	fgMzy = vSSzy / vSSzz;
	fgCzx = meanX - fgMzx * meanZ;
	fgCzy = meanY - fgMzy * meanZ;
	
	return true;
}


bool AliHLTMUONCalculations::FitLineToData(
		const AliHLTFloat32_t* x, const AliHLTFloat32_t* z, AliHLTUInt32_t n
	)
{
	/// A straight line is fitted in the X, Z data points using a least squares fit.
	/// http://mathworld.wolfram.com/LeastSquaresFitting.html
	/// If this method returns true, then the fitted parameters can fetched using the
	/// method calls Mzx() and Czx(). The line is then given by: x = Mzx() * z + Czx()
	/// \param x  This must point to the array of x data values.
	/// \param z  This must point to the array of z data values.
	/// \param n  Specifies the number of data points in the x and z arrays.
	/// \return  true if the line could be fitted or false otherwise.
	///     The reason for failure could be either too few data points or the slopes
	///     Mzx() would be infinite, implying a line that is perpendicular to the z axis.
	
	if (n < 2) return false;
	
	AliHLTFloat32_t sumX = 0;
	AliHLTFloat32_t sumZ = 0;
	for (AliHLTUInt32_t i = 0; i < n; i++)
	{
		sumX += x[i];
		sumZ += z[i];
	}
	AliHLTFloat32_t meanX = sumX / AliHLTFloat32_t(n);
	AliHLTFloat32_t meanZ = sumZ / AliHLTFloat32_t(n);
	
	AliHLTFloat32_t vSSzz = 0;
	AliHLTFloat32_t vSSzx = 0;
	for (AliHLTUInt32_t i = 0; i < n; i++)
	{
		vSSzz += (z[i] - meanZ)*(z[i] - meanZ);
		vSSzx += (z[i] - meanZ)*(x[i] - meanX);
	}
	
	// Calculate params for line x = fgMzx * z + fgCzx.
	if (vSSzz == 0) return false;
	fgMzx = vSSzx / vSSzz;
	fgCzx = meanX - fgMzx * meanZ;
	
	return true;
}


AliHLTFloat32_t AliHLTMUONCalculations::AliHLTMUONCalculations::ComputeChi2(
		const AliHLTFloat32_t* x, const AliHLTFloat32_t* y,
		const AliHLTFloat32_t* z, AliHLTUInt32_t n
	)
{
	/// Calculates the chi squared value for the set of data points given
	/// the fitted slope and coefficient parameters previously fitted by
	/// one of FitLine(x, y, z, n) or FitLineToTriggerRecord
	/// The fgSigmaX2 and fgSigmaY2 are used as the variance for the X and
	/// Y coordinates respectively. Note we assume that the covariance terms
	/// are zero.
	/// \param x  This must point to the array of x data values.
	/// \param y  This must point to the array of y data values.
	/// \param z  This must point to the array of z data values.
	/// \param n  Specifies the number of data points in the x, y and z arrays.
	/// \return  The chi squared value.
	
	AliHLTFloat32_t chi2 = 0;
	for (AliHLTUInt32_t i = 0; i < n; i++)
	{
		AliHLTFloat32_t residualX = fgMzx * z[i] + fgCzx - x[i];
		AliHLTFloat32_t residualY = fgMzy * z[i] + fgCzy - y[i];
		chi2 += residualX*residualX/fgSigmaX2 + residualY*residualY/fgSigmaY2;
	}
	return chi2;
}


AliHLTFloat32_t AliHLTMUONCalculations::AliHLTMUONCalculations::ComputeChi2(
		const AliHLTMUONTriggerRecordStruct& trigger,
		const bool hitset[4]
	)
{
	/// Calculates the chi squared value for trigger record using the hits
	/// indicated by the hitset array.
	/// \param trigger  The trigger record structure for which we compute the chi squared value.
	/// \param hitset  Flags indicating which hits to use and were set in the trigger record.
	/// \return  The chi squared value or -1 if it could not be calculated.
	
	if (not FitLine(trigger, hitset)) return -1;
	AliHLTFloat32_t chi2 = 0;
	for (AliHLTUInt32_t i = 0; i < 4; i++)
	{
		if (hitset[i])
		{
			AliHLTFloat32_t residualX = fgMzx * trigger.fHit[i].fZ + fgCzx - trigger.fHit[i].fX;
			AliHLTFloat32_t residualY = fgMzy * trigger.fHit[i].fZ + fgCzy - trigger.fHit[i].fY;
			chi2 += residualX*residualX/fgSigmaX2 + residualY*residualY/fgSigmaY2;
		}
	}
	return chi2;
}
