////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_CLUSTER_FINDER_CALLBACK_HPP
#define dHLT_ALIROOT_CLUSTER_FINDER_CALLBACK_HPP

#include "TObject.h"


namespace AliMUONHLT
{


class ClusterFinderCallback : public TObject
{
public:

	virtual void FoundClusters(const UInt_t numberfound) = 0;
	virtual void NoClustersFound() = 0;
	
	ClassDef(ClusterFinderCallback, 0);  // Abstract cluster finder callback class.
};


} // AliMUONHLT

#endif // dHLT_ALIROOT_CLUSTER_FINDER_CALLBACK_HPP
