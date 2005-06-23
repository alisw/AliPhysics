////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_TRIGGER_RECORD_HPP
#define dHLT_ALIROOT_TRIGGER_RECORD_HPP

#include <TObject.h>
#include "AliRoot/Point.hpp"
#include "Utils.hpp"

namespace AliMUONHLT
{


class TriggerRecord : public TObject
{
public:

	/* Default constructor initialises everything to zero and the trigger
	   number to -1.
	 */
	TriggerRecord();
	
	/* Creates a trigger record from the specified parameters.
	   Note: the trigger number must be greater or equal to zero. The particle
	   sign must also be one of the following values: -1, 0 or +1
	   Pt must be a positive number.
	   If these conditions are not met then an error message is displayed and
	   the object is filled like it is in the default constructor. 
	 */
	TriggerRecord(
			const Int_t triggernumber, const Int_t sign, const Float_t pt,
			const Point& station1point, const Point& station2point
		);

	virtual ~TriggerRecord() {};
	
	/* Get/Set method for the trigger number. 
	   The trigger number must be positive when assigning the trigger number.
	   If it is not then an error message is displayed and the internal value
	   remains untouched.
	 */
	void TriggerNumber(const Int_t value);
	Int_t TriggerNumber() const { return fTriggerNumber; };
	
	/* Get/Set method for the particle sign.
	   The particle sign must be one of the following values: -1, 0 or +1
	   If it is not then an error message is displayed and the internal value
	   remains untouched.
	 */
	void ParticleSign(const Int_t value);
	Int_t ParticleSign() const { return fParticleSign; };
	
	/* Get/Set method for the particle Pt, as measured by the L0 trigger.
	   The pt must be a positive number when assigning the pt.
	   If it is not then an error message is displayed and the internal value
	   remains untouched.
	 */
	void Pt(const Float_t value);
	Float_t Pt() const { return fPt; };
	
	/* Get/Set methods for the two trigger stations.
	 */
	void Station1Point(const Point& value) { fSt1Point = value; };
	Point& Station1Point() { return fSt1Point; };
	const Point& Station1Point() const { return fSt1Point; };
	void Station2Point(const Point& value) { fSt2Point = value; };
	Point& Station2Point() { return fSt2Point; };
	const Point& Station2Point() const { return fSt2Point; };

	// ostream operator usefull for text output.
	friend ostream& operator << (ostream& os, const TriggerRecord& r);

private:

	// Performs internal initialisation for the constructors.
	void Init();

	Int_t fTriggerNumber;  // The trigger number/index in AliMUONDataInterface.
	Int_t fParticleSign;   // -1 = negative, 0 = unknown, +1 = positive.
	Float_t fPt;           // Transverse momentum of the particle.
	Point fSt1Point;       // Coordinate on trigger station 1.
	Point fSt2Point;       // Coordinate on trigger station 2.

	ClassDef(TriggerRecord, 1)  // Trigger information.
};


} // AliMUONHLT

#endif // dHLT_ALIROOT_TRIGGER_RECORD_HPP
