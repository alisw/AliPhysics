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


AliHLTMUONClusterFinderProxy::AliHLTMUONClusterFinderProxy(AliHLTMUONClusterFinderInterface* client)
	: AliHLTMUONCoreClusterFinder(), AliHLTMUONClusterFinderCallback()
{
	clusterfinder = client;
}


void AliHLTMUONClusterFinderProxy::FindClusters(const AliHLTMUONCoreADCStream* /*stream*/)
{
	// TODO: perform conversion
	//ADCStream adc = AliHLTMUONConvert(stream);
	AliHLTMUONADCStream adc;
	DebugMsg(6, "AliHLTMUONClusterFinderProxy::FindClusters: " << adc);
	clusterfinder->FindClusters(&adc);
}


UInt AliHLTMUONClusterFinderProxy::FillClusterData(AliHLTMUONCoreClusterPoint* clusters, UInt arraysize)
{
	UInt result;
	AliHLTMUONPoint* points = new AliHLTMUONPoint[arraysize];
	try
	{
		DebugMsg(6, "AliHLTMUONClusterFinderProxy::FillClusterData");
		result = clusterfinder->FillClusterData(points, arraysize);
		for (UInt i = 0; i < arraysize; i++)
		{
			clusters[i] = AliHLTMUONConvert(points[i]);
			DebugMsg(6, "\tpoints[" << i << "] = " << points[i] );
		}
	}
	finally
	(
		delete [] points;
	)
	return result;
}


void AliHLTMUONClusterFinderProxy::Reset()
{
	DebugMsg(6, "AliHLTMUONClusterFinderProxy::Reset");
	clusterfinder->Reset();
}


void AliHLTMUONClusterFinderProxy::FoundClusters(UInt_t numberfound)
{
	DebugMsg(6, "AliHLTMUONClusterFinderProxy::FoundClusters");
	AliHLTMUONCoreClusterFinder::FoundClusters(numberfound);
}


void AliHLTMUONClusterFinderProxy::NoClustersFound()
{
	DebugMsg(6, "AliHLTMUONClusterFinderProxy::NoClustersFound");
	AliHLTMUONCoreClusterFinder::NoClustersFound();
}
