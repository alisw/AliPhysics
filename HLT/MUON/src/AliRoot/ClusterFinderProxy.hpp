////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCLUSTERFINDERPROXY_H
#define ALIHLTMUONCLUSTERFINDERPROXY_H

#include "Clustering/ClusterFinder.hpp"
#include "AliRoot/ClusterFinderCallback.hpp"
#include "AliRoot/ClusterFinderInterface.hpp"


class AliHLTMUONClusterFinderProxy : public AliHLTMUONCoreClusterFinder, public AliHLTMUONClusterFinderCallback
{
public:

	AliHLTMUONClusterFinderProxy(AliHLTMUONClusterFinderInterface* client);
	virtual ~AliHLTMUONClusterFinderProxy() {};

	// inherited methods from Clustering::ClusterFinder:
	virtual void FindClusters(const AliHLTMUONCoreADCStream* stream);
	virtual UInt FillClusterData(AliHLTMUONCoreClusterPoint* clusters, UInt arraysize);
	virtual void Reset();

	// inherited methods from AliMUONHLT::ClusterFinderCallback:
	virtual void FoundClusters(UInt_t numberfound);
	virtual void NoClustersFound();

private:
	// Do not allow copying.
	AliHLTMUONClusterFinderProxy(const AliHLTMUONClusterFinderProxy& /*object*/)
		: AliHLTMUONCoreClusterFinder(), AliHLTMUONClusterFinderCallback(),
		  fClusterFinder(NULL)
	{}

	AliHLTMUONClusterFinderProxy& operator = (const AliHLTMUONClusterFinderProxy& /*object*/)
	{
		return *this;
	}


	AliHLTMUONClusterFinderInterface* fClusterFinder;  // The clusterfinder we are proxying for.
};


#endif // ALIHLTMUONCLUSTERFINDERPROXY_H
