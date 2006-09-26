////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONDUMMYCLUSTERFINDER_H
#define ALIHLTMUONDUMMYCLUSTERFINDER_H

#ifndef __CINT__
#include "AliRoot/Point.hpp"
#include "AliRoot/ADCStream.hpp"
#endif // __CINT__


class AliHLTMUONDummyClusterFinder
{
public:

	AliHLTMUONDummyClusterFinder() : fInterface(this), fCallback(NULL)
	{
		fCallback = NULL;
	};

	virtual ~AliHLTMUONDummyClusterFinder() {};

	/* This is the starting point of the cluster finding algorithm.
	   Deriving cluster finders should implement all processing in this method
	   to find clusters in the specified ADC stream. When all clusters are found
	   the FoundClusters method should be called to indicate that processing is
	   complete. If no clusters could be found then call NoClustersFound instead.
	 */
	virtual void FindClusters(const AliHLTMUONADCStream* stream) = 0;

	/* After a call to FoundClusters this method will be called to retreive the
	   cluster points. The clusters array should be filled consecutively with 
	   the points that were found. However no more than 'arraysize' number of 
	   points should be written to the clusters array.
	   This method should also return the actual number of cluster points written
	   to the array.
	   If the number of clusters written is less that the number specified in the
	   'numberfound' parameter of the FoundClusters method, then this method will
	   be called again by the framework. Thus on successive calls to this method,
	   the cluster finder must resume writing clusters at the point it stopped at
	   in the previous call to FillClusterData.
	  */
	virtual UInt_t FillClusterData(AliHLTMUONPoint* clusters, UInt_t arraysize) = 0;

	/* This is called when the cluster finder should be reset to an initial state.
	   All extra internal memory allocated during processing should be released.
	 */
	virtual void Reset() = 0;

	/* Sets the ClusterFinderCallback callback interface.
	 */
	inline void SetCallback(AliHLTMUONClusterFinderCallback* callback)
	{
		this->fCallback = callback;
	};

	AliHLTMUONClusterFinderInterface* Interface()
	{
		return &fInterface;
	};

private:

	AliHLTMUONClusterFinderInterface fInterface;
	AliHLTMUONClusterFinderCallback* fCallback;

	// Hide copy constructor and assignment operator
	AliHLTMUONDummyClusterFinder(const AliHLTMUONDummyClusterFinder& /*clusterfinder*/)
		: fInterface(this), fCallback(NULL)
	{};

	AliHLTMUONDummyClusterFinder& operator = (const AliHLTMUONDummyClusterFinder& /*clusterfinder*/)
	{
		return *this;
	};
};


void AliHLTMUONClusterFinderInterface::FindClusters(const AliHLTMUONADCStream* stream)
{
	fClusterFinder->FindClusters(stream);
};

UInt_t AliHLTMUONClusterFinderInterface::FillClusterData(AliHLTMUONPoint* clusters, UInt_t arraysize)
{
	return fClusterFinder->FillClusterData(clusters, arraysize);
};

void AliHLTMUONClusterFinderInterface::Reset()
{
        fClusterFinder->Reset();
};

void AliHLTMUONClusterFinderInterface::SetCallback(AliHLTMUONClusterFinderCallback* callback)
{
	fClusterFinder->SetCallback(callback);
};

void AliHLTMUONMicrodHLT::SetClusterFinder(AliHLTMUONDummyClusterFinder* clusterfinder)
{
	SetClusterFinder(clusterfinder->Interface());
};


#endif // ALIHLTMUONDUMMYCLUSTERFINDER_H
