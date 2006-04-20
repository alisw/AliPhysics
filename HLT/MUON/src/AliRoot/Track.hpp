////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONTRACK_H
#define ALIHLTMUONTRACK_H

#include <TObject.h>
#include <Riostream.h>

#include "AliRoot/Point.hpp"
#include "AliRoot/Region.hpp"


class AliHLTMUONTrack : public TObject
{
	// ostream operator usefull for text output.
	friend ostream& operator << (ostream& os, const AliHLTMUONTrack& t);
	
public:

	/* Default constructor initialises everything to zero and fTriggerID to -1.
	 */
	AliHLTMUONTrack();
	
	/* Create a track object from the given parameters.
	   This constructor checks that momentum >= pt and sign is one of the
	   following values: -1, 0 or +1. If these conditions are violated then
	   the internal data is initialised as in the default constructor and an
	   error message is displayed.
	 */
	AliHLTMUONTrack(
		Int_t triggerid, Int_t sign, Float_t momentum, Float_t pt,
		const AliHLTMUONPoint hits[10], const AliHLTMUONRegion regions[10]
	);

	virtual ~AliHLTMUONTrack() {}
	
	
	// Get/et methods for the trigger ID.
	void TriggerID(Int_t value) { fTriggerID = value; };
	Int_t TriggerID() const { return fTriggerID; };
	
	/* Get/Set method for the particle sign. The particle sign must be one
	   of the following values: -1, 0 or +1
	   If the new value is not in this range then an error message is
	   displayed and the internal value remain unchanged.
	 */
	void ParticleSign(Int_t value);
	Int_t ParticleSign() const { return fParticleSign; };
	
	/* The get and set methods for the momentum and transverse momentum pt.
	   These methods check that the momentum is always equal or larger than
	   the pt. If not then the internal values are left unchanged and an
	   error message is displayed. The numbers must also be positive.
	 */
	void P(Float_t value);
	Float_t P() const { return fP; };
	void Pt(Float_t value);
	Float_t Pt() const { return fPt; };
	
	/* Returns the hit point for the specified chamber.
	   If the chamber number in out of bounds the point on the first
	   chamber is returned and an error message displayed.
	 */
	AliHLTMUONPoint& Hit(UInt_t chamber);
	const AliHLTMUONPoint& Hit(UInt_t chamber) const;
	
	/* Set method for hits. The chamber must be in the range [0..9]
	 */
	void Hit(UInt_t chamber, const AliHLTMUONPoint& value);
	
	/* Returns the region of interest for the specified chamber.
	   If the chamber number in out of bounds the region on the first
	   chamber is returned and an error message displayed.
	 */
	AliHLTMUONRegion& RegionOfInterest(UInt_t chamber);
	const AliHLTMUONRegion& RegionOfInterest(UInt_t chamber) const;
	
	/* Set method for regions. The chamber must be in the range [0..9]
	 */
	void RegionOfInterest(UInt_t chamber, const AliHLTMUONRegion& value);
	
	/* Checks to see if the all the hits are within their respective regions
	   of interest for each chamber. kTRUE is returned if they are and kFALSE
	   otherwise.
	 */
	Bool_t HitsInRegions() const;

private:

	// Internal initialisation routine used by the constructors.
	void Init();

	Int_t fTriggerID;     // The ID number of the trigger record this track was found for.
	Int_t fParticleSign;  // The sign of the particle.
	Float_t fP;       // momentum of particle.
	Float_t fPt;      // transverse momentum.
	AliHLTMUONPoint fHit[10];   // Fitted track hit coordinates for the 10 tracking chambers.
	AliHLTMUONRegion fRegionOfInterest[10];  // Region of interest per chamber.

	ClassDef(AliHLTMUONTrack, 1)  // Track data.
};


#endif // ALIHLTMUONTRACK_H
