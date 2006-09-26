////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCLUSTERFINDERINTERFACE_H
#define ALIHLTMUONCLUSTERFINDERINTERFACE_H

#include "Rtypes.h"

class AliHLTMUONPoint;
class AliHLTMUONADCStream;
class AliHLTMUONClusterFinderCallback;
class AliHLTMUONDummyClusterFinder;


class AliHLTMUONClusterFinderInterface
{
public:
	AliHLTMUONClusterFinderInterface(AliHLTMUONDummyClusterFinder* clusterfinder)
		: fClusterFinder(clusterfinder)
	{
		fClusterFinder = clusterfinder;
	};
	
	const AliHLTMUONDummyClusterFinder* GetClusterFinder() const
	{
		return fClusterFinder;
	};
	
        void FindClusters(const AliHLTMUONADCStream* stream);
 	UInt_t FillClusterData(AliHLTMUONPoint* clusters, UInt_t arraysize);
 	void Reset();
        void SetCallback(AliHLTMUONClusterFinderCallback* callback);

private:

	
	AliHLTMUONClusterFinderInterface(const AliHLTMUONClusterFinderInterface& /*clusterfinder*/)
		: fClusterFinder(NULL)
	{}

	AliHLTMUONClusterFinderInterface& operator = (const AliHLTMUONClusterFinderInterface& /*clusterfinder*/)
	{
		return *this;
	}


	AliHLTMUONDummyClusterFinder* fClusterFinder;   //! Pointer to interpreted cluster finder class.
};


#endif // ALIHLTMUONCLUSTERFINDERINTERFACE_H
