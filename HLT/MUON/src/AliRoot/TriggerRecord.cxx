////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/TriggerRecord.hpp"
#include "TMath.h"

ClassImp(AliMUONHLT::TriggerRecord)

namespace AliMUONHLT
{

TriggerRecord::TriggerRecord()
{
	Init();
}


TriggerRecord::TriggerRecord(
		const Int_t triggernumber, const Int_t sign, const Float_t pt,
		const Point& station1point, const Point& station2point
	)
{
	if (triggernumber < 0)
	{
		Init();
		Error("TriggerRecord",
			"The trigger number must be a positive number. Got: %d",
			triggernumber
		);
	}
	else if (sign < -1 or +1 < sign)
	{
		Init();
		Error("TriggerRecord",
			"The particle sign must a value of -1, 0 or +1. Got: %d",
			sign
		);
	}
	else if (pt < 0.0)
	{
		Init();
		Error("TriggerRecord",
			"The transverse momentum must be a positive number. Got: %f",
			pt
		);
	}
	else
	{
		fTriggerNumber = triggernumber;
		fParticleSign = sign;
		fPt = pt;
		fSt1Point = station1point;
		fSt2Point = station2point;
	}
}


void TriggerRecord::Init()
{
	fTriggerNumber = -1;
	fParticleSign = 0;
	fPt = 0.0;
}


void TriggerRecord::TriggerNumber(const Int_t value)
{
	if (value >= 0)
		fTriggerNumber = value;
	else
		Error("TriggerNumber",
			"The trigger number must be a positive number. Got: %d",
			value
		);
}


void TriggerRecord::ParticleSign(const Int_t value)
{
	if (-1 <= value and value <= +1)
		fParticleSign = value;
	else
		Error("ParticleSign",
			"The particle sign must a value of -1, 0 or +1. Got: %d",
			value
		);
}


void TriggerRecord::Pt(const Float_t value)
{
	if (value >= 0)
		fPt = value;
	else
		Error("Pt",
			"The transverse momentum must be a positive number. Got: %f",
			value
		);
}


ostream& operator << (ostream& os, const TriggerRecord& r)
{
	os << "{trig#: " << r.fTriggerNumber << ", sign: " << r.fParticleSign
	   << ", pt: " << r.fPt << ", st1: " << r.fSt1Point << ", st2: "
	   << r.fSt2Point << "}";
	return os;
}


} // AliMUONHLT
