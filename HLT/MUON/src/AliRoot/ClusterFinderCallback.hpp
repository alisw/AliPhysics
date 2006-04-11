////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCLUSTERFINDERCALLBACK_H
#define ALIHLTMUONCLUSTERFINDERCALLBACK_H

#include "TObject.h"


class AliHLTMUONClusterFinderCallback : public TObject
{
public:

	virtual ~AliHLTMUONClusterFinderCallback() {};

	virtual void FoundClusters(UInt_t numberfound) = 0;
	virtual void NoClustersFound() = 0;
	
	ClassDef(AliHLTMUONClusterFinderCallback, 0)  // Abstract cluster finder callback class.
};


#endif // ALIHLTMUONCLUSTERFINDERCALLBACK_H
