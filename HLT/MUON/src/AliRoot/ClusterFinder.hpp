////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_CLUSTERING_CLUSTER_FINDER_HPP
#define dHLT_CLUSTERING_CLUSTER_FINDER_HPP

#ifndef __CINT__
#include "BasicTypes.hpp"
#include "Cluster.hpp"
#include "ADCStream.hpp"
#include "Utils.hpp"
#endif // __CINT__

namespace AliMUONHLT
{

class ClusterFinder
{
public:

	ClusterFinder() :fInterface(this)
	{
		fCallback = NULL;
	};

	virtual ~ClusterFinder() {};

	/* This is the starting point of the cluster finding algorithm.
	   Deriving cluster finders should implement all processing in this method
	   to find clusters in the specified ADC stream. When all clusters are found
	   the FoundClusters method should be called to indicate that processing is
	   complete. If no clusters could be found then call NoClustersFound instead.
	 */
	virtual void FindClusters(const ADCStream * stream) = 0;

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
	virtual UInt_t FillClusterData(Point * clusters, UInt_t arraysize) = 0;

	/* This is called when the cluster finder should be reset to an initial state.
	   All extra internal memory allocated during processing should be released.
	 */
	virtual void Reset() = 0;

	/* Sets the ClusterFinderCallback callback interface.
	 */
	inline void SetCallback(ClusterFinderCallback * callback) {
		this->fCallback = callback;
	};

	ClusterFinderInterface* Interface()
	{
		return &fInterface;
	};

private:

	ClusterFinderInterface fInterface;
	ClusterFinderCallback * fCallback;
};


void ClusterFinderInterface::FindClusters(const ADCStream * stream)
{
	fClusterFinder->FindClusters(stream);
};

UInt_t ClusterFinderInterface::FillClusterData(Point * clusters, UInt_t arraysize)
{
	return fClusterFinder->FillClusterData(clusters,arraysize);
};

void ClusterFinderInterface::Reset()
{
        fClusterFinder->Reset();
};

void ClusterFinderInterface::SetCallback(ClusterFinderCallback* callback)
{
	fClusterFinder->SetCallback(callback);
};

void MicrodHLT::SetClusterFinder(ClusterFinder* clusterfinder)
{
	SetClusterFinder(clusterfinder->Interface());
};

}				// AliMUONHLT

#endif				// dHLT_CLUSTERING_CLUSTER_FINDER_HPP
