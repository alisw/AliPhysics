////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_CLUSTER_FINDER_PROXY_HPP
#define dHLT_ALIROOT_CLUSTER_FINDER_PROXY_HPP

#include "Clustering/ClusterFinder.hpp"
#include "AliRoot/ClusterFinderCallback.hpp"
#include "AliRoot/ClusterFinderInterface.hpp"

namespace dHLT
{
namespace AliRoot
{


class ClusterFinderProxy : public Clustering::ClusterFinder, public AliMUONHLT::ClusterFinderCallback
{
public:

	ClusterFinderProxy(AliMUONHLT::ClusterFinderInterface* client);
	virtual ~ClusterFinderProxy() {};

	// inherited methods from Clustering::ClusterFinder:
	virtual void FindClusters(const ADCStream* stream);
	virtual UInt FillClusterData(ClusterPoint* clusters, const UInt arraysize);
	virtual void Reset();

	// inherited methods from AliMUONHLT::ClusterFinderCallback:
	virtual void FoundClusters(const UInt_t numberfound);
	virtual void NoClustersFound();

private:

	AliMUONHLT::ClusterFinderInterface* clusterfinder;  // The clusterfinder we are proxying for.
};


} // AliRoot
} // dHLT

#endif // dHLT_ALIROOT_CLUSTER_FINDER_PROXY_HPP
