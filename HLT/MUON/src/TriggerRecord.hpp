////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCORETRIGGERRECORD_H
#define ALIHLTMUONCORETRIGGERRECORD_H

#include "BasicTypes.hpp"
#include "Point.hpp"


/* The sign of the particle as given by L0.
 */
enum AliHLTMUONCoreParticleSign
{
    kSignMinus   = -1,
    kUnknownSign = 0,
    kSignPlus    = 1
};


typedef UInt AliHLTMUONCoreTriggerRecordID;


/* Data structure containing information about L0 validated trigger hit.
 */
struct AliHLTMUONCoreTriggerRecord
{
	AliHLTMUONCoreParticleSign fSign;     // The sign of the particle.
	Float fPt;              // Transverse momentum of the particle.
	AliHLTMUONCorePoint fStation1impact;  // Impact point of particle on trigger station 1.
	AliHLTMUONCorePoint fStation2impact;  // Impact point of particle on trigger station 2.


	/* Default constructor.
	   Sets the fSign to UnknownSign, fPt to -1 and the impact points are set to zero.
	 */
	AliHLTMUONCoreTriggerRecord()
		: fStation1impact(), fStation2impact()
	{
		fSign = kUnknownSign;
		fPt = -1.0;
	};

	/* Creates a trigger record with the specifed particle sign, pt and impact points.
	   The impactpoint1 corresponds to trigger station 1 and simmilarly impactpoint2
	   corresponds to station 2.
	 */
	AliHLTMUONCoreTriggerRecord(
			const AliHLTMUONCoreParticleSign sign, const Float pt,
			const AliHLTMUONCorePoint impactpoint1, const AliHLTMUONCorePoint impactpoint2
		)
	{
		fSign = sign;
		fPt = pt;
		fStation1impact = impactpoint1;
		fStation2impact = impactpoint2;
	};
};


#endif // ALIHLTMUONCORETRIGGERRECORD_H
