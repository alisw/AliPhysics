////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/ClusterFinderProxy.hpp"
#include "AliRoot/ADCStream.hpp"
#include "AliRoot/convert.hpp"
#include "Utils.hpp"
#include "new.hpp"

namespace dHLT
{
namespace AliRoot
{


ClusterFinderProxy::ClusterFinderProxy(AliMUONHLT::ClusterFinderInterface* client)
	: Clustering::ClusterFinder(), AliMUONHLT::ClusterFinderCallback()
{
	clusterfinder = client;
}


void ClusterFinderProxy::FindClusters(const ADCStream* stream)
{
	// TODO: perform conversion
	//AliMUONHLT::ADCStream stream = Convert(stream);
	AliMUONHLT::ADCStream adc;
	DebugMsg(6, "ClusterFinderProxy::FindClusters: " << adc);
	clusterfinder->FindClusters(&adc);
}


UInt ClusterFinderProxy::FillClusterData(ClusterPoint* clusters, const UInt arraysize)
{
	UInt result;
	AliMUONHLT::Point* points = new AliMUONHLT::Point[arraysize];
	try
	{
		DebugMsg(6, "ClusterFinderProxy::FillClusterData");
		result = clusterfinder->FillClusterData(points, arraysize);
		for (UInt i = 0; i < arraysize; i++)
		{
			clusters[i] = Convert(points[i]);
			DebugMsg(6, "\tpoints[" << i << "] = " << points[i] );
		}
	}
	finally
	(
		delete [] points;
	)
	return result;
}


void ClusterFinderProxy::Reset()
{
	DebugMsg(6, "ClusterFinderProxy::Reset");
	clusterfinder->Reset();
}


void ClusterFinderProxy::FoundClusters(const UInt_t numberfound)
{
	DebugMsg(6, "ClusterFinderProxy::FoundClusters");
	Clustering::ClusterFinder::FoundClusters(numberfound);
}


void ClusterFinderProxy::NoClustersFound()
{
	DebugMsg(6, "ClusterFinderProxy::NoClustersFound");
	Clustering::ClusterFinder::NoClustersFound();
}


} // AliRoot
} // dHLT
