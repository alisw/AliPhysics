////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_CLUSTER_FINDER_INTERFACE_HPP
#define dHLT_ALIROOT_CLUSTER_FINDER_INTERFACE_HPP

#include "Rtypes.h"


namespace AliMUONHLT
{

class Point;
class ADCStream;
class ClusterFinder;
class ClusterFinderCallback;


class ClusterFinderInterface
{
public:
	ClusterFinderInterface(ClusterFinder* clusterfinder)
	{
		fClusterFinder = clusterfinder;
	};
	
	const ClusterFinder* GetClusterFinder() const
	{
		return fClusterFinder;
	};
	
	void FindClusters(const ADCStream* stream);
	UInt_t FillClusterData(Point* clusters, const UInt_t arraysize);
	void Reset();
	void SetCallback(ClusterFinderCallback* callback);

private:

	ClusterFinder* fClusterFinder;   //! Pointer to interpreted cluster finder class.
};


} // AliMUONHLT

#endif // dHLT_ALIROOT_CLUSTER_FINDER_INTERFACE_HPP
