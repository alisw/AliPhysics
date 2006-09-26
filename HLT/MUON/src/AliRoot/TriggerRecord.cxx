////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/TriggerRecord.hpp"
#include "TMath.h"

ClassImp(AliHLTMUONTriggerRecord)


AliHLTMUONTriggerRecord::AliHLTMUONTriggerRecord() :
	TObject(), fTriggerNumber(-1), fParticleSign(0), fPt(0), fSt1Point(), fSt2Point()
{
// Default constructor initialises everything to zero and the trigger number to -1.

	Init();
}


AliHLTMUONTriggerRecord::AliHLTMUONTriggerRecord(
		Int_t triggernumber, Int_t sign, Float_t pt,
		const AliHLTMUONPoint& station1point, const AliHLTMUONPoint& station2point
	) :
	TObject(), fTriggerNumber(-1), fParticleSign(0), fPt(0), fSt1Point(), fSt2Point()
{
// Creates a trigger record from the specified parameters.
// Note: the trigger number must be greater or equal to zero. The particle
// sign must also be one of the following values: -1, 0 or +1
// Pt must be a positive number.
// If these conditions are not met then an error message is displayed and
// the object is filled like it is in the default constructor. 

	if (triggernumber < 0)
	{
		Init();
		Error("AliHLTMUONTriggerRecord",
			"The trigger number must be a positive number. Got: %d",
			triggernumber
		);
	}
	else if (sign < -1 || +1 < sign)
	{
		Init();
		Error("AliHLTMUONTriggerRecord",
			"The particle sign must a value of -1, 0 or +1. Got: %d",
			sign
		);
	}
	else if (pt < 0.0)
	{
		Init();
		Error("AliHLTMUONTriggerRecord",
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


void AliHLTMUONTriggerRecord::Init()
{
// Performs internal initialisation for the constructors.

	fTriggerNumber = -1;
	fParticleSign = 0;
	fPt = 0.0;
}


void AliHLTMUONTriggerRecord::TriggerNumber(Int_t value)
{
// Set method for the trigger number. 
// The trigger number must be positive when assigning the trigger number.
// If it is not then an error message is displayed and the internal value
// remains untouched.

	if (value >= 0)
		fTriggerNumber = value;
	else
		Error("TriggerNumber",
			"The trigger number must be a positive number. Got: %d",
			value
		);
}


void AliHLTMUONTriggerRecord::ParticleSign(Int_t value)
{
// Set method for the particle sign.
// The particle sign must be one of the following values: -1, 0 or +1
// If it is not then an error message is displayed and the internal value
// remains untouched.

	if (-1 <= value && value <= +1)
		fParticleSign = value;
	else
		Error("ParticleSign",
			"The particle sign must a value of -1, 0 or +1. Got: %d",
			value
		);
}


void AliHLTMUONTriggerRecord::Pt(Float_t value)
{
// Set method for the particle Pt, as measured by the L0 trigger.
// The pt must be a positive number when assigning the pt.
// If it is not then an error message is displayed and the internal value
// remains untouched.

	if (value >= 0)
		fPt = value;
	else
		Error("Pt",
			"The transverse momentum must be a positive number. Got: %f",
			value
		);
}


ostream& operator << (ostream& os, const AliHLTMUONTriggerRecord& r)
{
	os << "{trig#: " << r.fTriggerNumber << ", sign: " << r.fParticleSign
	   << ", pt: " << r.fPt << ", st1: " << r.fSt1Point << ", st2: "
	   << r.fSt2Point << "}";
	return os;
}

