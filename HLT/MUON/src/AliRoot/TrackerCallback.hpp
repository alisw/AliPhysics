////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONTRACKERCALLBACK_H
#define ALIHLTMUONTRACKERCALLBACK_H

#include "TObject.h"


class AliHLTMUONTrackerCallback : public TObject
{
public:

	virtual void RequestClusters(
			Float_t left, Float_t right, Float_t bottom, Float_t top,
			Int_t chamber, const void* tag = NULL
		) = 0;
	virtual void EndOfClusterRequests() = 0;
	virtual void FoundTrack() = 0;
	virtual void NoTrackFound() = 0;
	
	ClassDef(AliHLTMUONTrackerCallback, 0)  // Abstract tracker callback class.
};


#endif // ALIHLTMUONTRACKERCALLBACK_H
