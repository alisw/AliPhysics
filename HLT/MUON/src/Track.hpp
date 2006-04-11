////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCORETRACK_H
#define ALIHLTMUONCORETRACK_H

#include "Point.hpp"
#include "RegionOfInterest.hpp"
#include "TriggerRecord.hpp"


typedef UInt AliHLTMUONCoreTrackID;


struct AliHLTMUONCoreTrack
{

	AliHLTMUONCoreTriggerRecordID fTriggerid;
	AliHLTMUONCoreParticleSign fSign;
	Float fP;   // momentum.
	Float fPt;  // transverse momentum.
	AliHLTMUONCorePoint fPoint[10];  // Computed track coordinates on the 10 tracking chambers.
	AliHLTMUONCoreROI fRegion[10];   // Regions of interest from which clusters were used to compute this track.

};


#endif // ALIHLTMUONCORETRACK_H
