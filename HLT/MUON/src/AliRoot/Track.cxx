////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/Track.hpp"
#include "TMath.h"

ClassImp(AliMUONHLT::Track);

namespace AliMUONHLT
{


Track::Track() : TObject()
{
	Init();
};


Track::Track(
		const Int_t triggerid, const Int_t sign, const Float_t momentum, const Float_t pt,
		const Point hits[10], const Region regions[10]
	) : TObject()
{
	if (sign < -1 or +1 < sign)
	{
		Init();
		Error("Track", "The particle sign was not one of -1, 0 or +1. Got %d", sign);
	}
	else if (momentum < pt)
	{
		Init();
		Error("Track", "The momentum (%f) must be larger or equal to the pt (%f).",
			momentum, pt
		);
	}
	else if (pt < 0.0)
	{
		Init();
		Error("Track", "The pt must be a positive number. Got: %f", pt);
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
		};
	};
};


void Track::Init()
{
	fTriggerID = -1;
	fParticleSign = 0;
	fP = fPt = 0.0;
};


void Track::ParticleSign(const Int_t value)
{
	if (-1 <= value and value <= +1)
		fParticleSign = value;
	else
		Error("ParticleSign",
			"The particle sign must be a value of -1, 0 or +1, but got %d",
			value
		);
};


void Track::P(const Float_t value)
{
	if (value >= fPt)
		fP = value;
	else
		Error("P",
			"Trying to assing momentum (%f) which is smaller than the pt value (%f).",
			value, fPt
		);
};

void Track::Pt(const Float_t value)
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


Point& Track::Hit(const UInt_t chamber)
{
	if (0 <= chamber and chamber < 10)
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


const Point& Track::Hit(const UInt_t chamber) const
{
	if (0 <= chamber and chamber < 10)
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


void Track::Hit(const UInt_t chamber, const Point& value)
{
	if (0 <= chamber and chamber < 10)
		fHit[chamber] = value;
	else
		Error("Hit",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
};


Region& Track::RegionOfInterest(const UInt_t chamber)
{
	if (0 <= chamber and chamber < 10)
		return fRegionOfInterest[chamber];
	else
	{
		Error("RegionOfInterest",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
		return fRegionOfInterest[0];
	};
};


const Region& Track::RegionOfInterest(const UInt_t chamber) const
{
	if (0 <= chamber and chamber < 10)
		return fRegionOfInterest[chamber];
	else
	{
		Error("RegionOfInterest",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
		return fRegionOfInterest[0];
	};
};


void Track::RegionOfInterest(const UInt_t chamber, const Region& value)
{
	if (0 <= chamber and chamber < 10)
		fRegionOfInterest[chamber] = value;
	else
		Error("RegionOfInterest",
			"The chamber is out of range. Got: %d, but should be in [0..9].",
			chamber
		);
};


Bool_t Track::HitsInRegions() const
{
	for (Int_t i = 0; i < 10; i++)
	{
		if ( not fRegionOfInterest[i].Contains(fHit[i]) )
			return kFALSE;
	};
	return kTRUE;
};


std::ostream& operator << (std::ostream& os, const Track& t)
{
	os << "{trigid: " << t.fTriggerID << ", sign: " << t.fParticleSign
	   << ", p: " << t.fP << ", pt: " << t.fPt << "}";
	return os;
};


}; // AliMUONHLT
