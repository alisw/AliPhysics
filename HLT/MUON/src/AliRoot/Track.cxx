////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/Track.hpp"
#include <TMath.h>
#include "AliHLTMUONUtils.h"

ClassImp(AliHLTMUONTrack)


AliHLTMUONTrack::AliHLTMUONTrack() :
	TObject(), fTriggerID(-1), fParticleSign(0), fP(0), fPt(0)
{
// Default constructor initialises everything to zero and fTriggerID to -1.

	Init();
}


AliHLTMUONTrack::AliHLTMUONTrack(
		Int_t triggerid, Int_t sign, Float_t momentum, Float_t pt,
		const AliHLTMUONPoint hits[10], const AliHLTMUONRegion regions[10]
	) :
	TObject(), fTriggerID(-1), fParticleSign(0), fP(0), fPt(0)
{
// Create a track object from the given parameters.
// This constructor checks that momentum >= pt and sign is one of the
// following values: -1, 0 or +1. If these conditions are violated then
// the internal data is initialised as in the default constructor and an
// error message is displayed.

	if (sign < -1 || +1 < sign)
	{
		Init();
		Error("AliHLTMUONTrack", "The particle sign was not one of -1, 0 or +1. Got %d", sign);
	}
	else if (momentum < pt)
	{
		Init();
		Error("AliHLTMUONTrack", "The momentum (%f) must be larger or equal to the pt (%f).",
			momentum, pt
		);
	}
	else if (pt < 0.0)
	{
		Init();
		Error("AliHLTMUONTrack", "The pt must be a positive number. Got: %f", pt);
	}
	else
	{
		fTriggerID = triggerid;
		fParticleSign = sign;
		fP = momentum;
		fPt = pt;
		for (Int_t i = 0; i < 10; i++)
		{
			fHit[i] = hits[i];
			fRegionOfInterest[i] = regions[i];
		}
	}
}


void AliHLTMUONTrack::Init()
{
// Internal initialisation routine used by the constructors.

	fTriggerID = -1;
	fParticleSign = 0;
	fP = fPt = 0.0;
}


void AliHLTMUONTrack::ParticleSign(Int_t value)
{
// Set method for the particle sign. The particle sign must be one
// of the following values: -1, 0 or +1
// If the new value is not in this range then an error message is
// displayed and the internal value remain unchanged.

	if (-1 <= value && value <= +1)
		fParticleSign = value;
	else
		Error("ParticleSign",
			"The particle sign must be a value of -1, 0 or +1, but got %d",
			value
		);
}


void AliHLTMUONTrack::P(Float_t value)
{
// The set method for the momentum.
// This method checks that the momentum is always equal or larger than
// the pt. If not then the internal values are left unchanged and an
// error message is displayed. The numbers must also be positive.

	if (value >= fPt)
		fP = value;
	else
		Error("P",
			"Trying to assing momentum (%f) which is smaller than the pt value (%f).",
			value, fPt
		);
}

void AliHLTMUONTrack::Pt(Float_t value)
{
// The set method for the pt.
// This method checks that the momentum is always equal or larger than
// the pt. If not then the internal values are left unchanged and an
// error message is displayed. The numbers must also be positive.

	if (value >= 0.0)
	{
		if (value <= fP)
			fPt = value;
		else
			Error("Pt",
				"Trying to assign pt (%f) which is larger than the momentum value (%f).",
				value, fP
			);
	}
	else
		Error("Pt", "Cannot have a negative value pt. Got: %f", value);
};


AliHLTMUONPoint& AliHLTMUONTrack::Hit(UInt_t chamber)
{
// Returns the hit point for the specified chamber.
// If the chamber number in out of bounds the point on the first
// chamber is returned and an error message displayed.

	if (chamber < 10)
		return fHit[chamber];
	else
	{
		Error("Hit",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
		return fHit[0];
	}
}


const AliHLTMUONPoint& AliHLTMUONTrack::Hit(UInt_t chamber) const
{
// Returns a constant hit object for the specified chamber.
// If the chamber number in out of bounds the point on the first
// chamber is returned and an error message displayed.

	if (chamber < 10)
		return fHit[chamber];
	else
	{
		Error("Hit",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
		return fHit[0];
	};
};


void AliHLTMUONTrack::Hit(UInt_t chamber, const AliHLTMUONPoint& value)
{
// Set method for hits. The chamber must be in the range [0..9]

	if (chamber < 10)
		fHit[chamber] = value;
	else
		Error("Hit",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
}


AliHLTMUONRegion& AliHLTMUONTrack::RegionOfInterest(UInt_t chamber)
{
// Returns the region of interest for the specified chamber.
// If the chamber number in out of bounds the region on the first
// chamber is returned and an error message displayed.

	if (chamber < 10)
		return fRegionOfInterest[chamber];
	else
	{
		Error("RegionOfInterest",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
		return fRegionOfInterest[0];
	}
}


const AliHLTMUONRegion& AliHLTMUONTrack::RegionOfInterest(UInt_t chamber) const
{
// Returns the constatn region of interest object for the specified chamber.
// If the chamber number in out of bounds the region on the first
// chamber is returned and an error message displayed.

	if (chamber < 10)
		return fRegionOfInterest[chamber];
	else
	{
		Error("RegionOfInterest",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
		return fRegionOfInterest[0];
	}
}


void AliHLTMUONTrack::RegionOfInterest(UInt_t chamber, const AliHLTMUONRegion& value)
{
// Set method for regions. The chamber must be in the range [0..9]

	if (chamber < 10)
		fRegionOfInterest[chamber] = value;
	else
		Error("RegionOfInterest",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
}


Bool_t AliHLTMUONTrack::HitsInRegions() const
{
// Checks to see if the all the hits are within their respective regions
// of interest for each chamber. kTRUE is returned if they are and kFALSE
// otherwise.

	for (Int_t i = 0; i < 10; i++)
	{
		if ( ! fRegionOfInterest[i].Contains(fHit[i]) )
			return kFALSE;
	}
	return kTRUE;
}


ostream& operator << (ostream& os, const AliHLTMUONTrack& t)
{
	os << "{trigid: " << t.fTriggerID << ", sign: " << t.fParticleSign
	   << ", p: " << t.fP << ", pt: " << t.fPt << "}";
	return os;
}

