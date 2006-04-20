////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONTRIGGERRECORD_H
#define ALIHLTMUONTRIGGERRECORD_H

#include <TObject.h>
#include "AliRoot/Point.hpp"
#include "Utils.hpp"


class AliHLTMUONTriggerRecord : public TObject
{
	// ostream operator usefull for text output.
	friend ostream& operator << (ostream& os, const AliHLTMUONTriggerRecord& r);
	
public:

	/* Default constructor initialises everything to zero and the trigger
	   number to -1.
	 */
	AliHLTMUONTriggerRecord();
	
	/* Creates a trigger record from the specified parameters.
	   Note: the trigger number must be greater or equal to zero. The particle
	   sign must also be one of the following values: -1, 0 or +1
	   Pt must be a positive number.
	   If these conditions are not met then an error message is displayed and
	   the object is filled like it is in the default constructor. 
	 */
	AliHLTMUONTriggerRecord(
			Int_t triggernumber, Int_t sign, Float_t pt,
			const AliHLTMUONPoint& station1point, const AliHLTMUONPoint& station2point
		);

	virtual ~AliHLTMUONTriggerRecord() {};
	
	/* Get/Set method for the trigger number. 
	   The trigger number must be positive when assigning the trigger number.
	   If it is not then an error message is displayed and the internal value
	   remains untouched.
	 */
	void TriggerNumber(Int_t value);
	Int_t TriggerNumber() const { return fTriggerNumber; };
	
	/* Get/Set method for the particle sign.
	   The particle sign must be one of the following values: -1, 0 or +1
	   If it is not then an error message is displayed and the internal value
	   remains untouched.
	 */
	void ParticleSign(Int_t value);
	Int_t ParticleSign() const { return fParticleSign; };
	
	/* Get/Set method for the particle Pt, as measured by the L0 trigger.
	   The pt must be a positive number when assigning the pt.
	   If it is not then an error message is displayed and the internal value
	   remains untouched.
	 */
	void Pt(Float_t value);
	Float_t Pt() const { return fPt; };
	
	/* Get/Set methods for the two trigger stations.
	 */
	void Station1Point(const AliHLTMUONPoint& value) { fSt1Point = value; };
	AliHLTMUONPoint& Station1Point() { return fSt1Point; };
	const AliHLTMUONPoint& Station1Point() const { return fSt1Point; };
	void Station2Point(const AliHLTMUONPoint& value) { fSt2Point = value; };
	AliHLTMUONPoint& Station2Point() { return fSt2Point; };
	const AliHLTMUONPoint& Station2Point() const { return fSt2Point; };
	

private:

	// Performs internal initialisation for the constructors.
	void Init();

	Int_t fTriggerNumber;  // The trigger number/index in AliMUONDataInterface.
	Int_t fParticleSign;   // -1 = negative, 0 = unknown, +1 = positive.
	Float_t fPt;           // Transverse momentum of the particle.
	AliHLTMUONPoint fSt1Point;       // Coordinate on trigger station 1.
	AliHLTMUONPoint fSt2Point;       // Coordinate on trigger station 2.

	ClassDef(AliHLTMUONTriggerRecord, 1)  // Trigger information.
};


#endif // ALIHLTMUONTRIGGERRECORD_H
