////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/Track.hpp"
#include <TMath.h>
#include "Utils.hpp"

ClassImp(AliHLTMUONTrack)


AliHLTMUONTrack::AliHLTMUONTrack() : TObject()
{
	Init();
}


AliHLTMUONTrack::AliHLTMUONTrack(
		Int_t triggerid, Int_t sign, Float_t momentum, Float_t pt,
		const AliHLTMUONPoint hits[10], const AliHLTMUONRegion regions[10]
	) : TObject()
{
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
	fTriggerID = -1;
	fParticleSign = 0;
	fP = fPt = 0.0;
}


void AliHLTMUONTrack::ParticleSign(Int_t value)
{
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

