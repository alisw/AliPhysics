////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCORECLUSTERFINDER_H
#define ALIHLTMUONCORECLUSTERFINDER_H

#include "AliHLTMUONBasicTypes.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONCoreCluster.h"
#include "AliHLTMUONCoreADCStream.h"


class AliHLTMUONCoreClusterFinder;


class AliHLTMUONCoreClusterFinderCallback
{
public:

	virtual ~AliHLTMUONCoreClusterFinderCallback() {};

	/* Called when the cluster finder has found all clusters in the ADC stream.
	   At this point the ADC stream is no longer is use by the cluster finder and
	   the stream can be released.
	   The numberfound parameter indicated how many clusters were actually found.
	*/
	virtual void FoundClusters(AliHLTMUONCoreClusterFinder* clusterfinder, UInt numberfound) = 0;

	/* Called when the cluster finder has finished its job however no clusters were
	   found in the ADC stream. At this point the ADC stream is no longer is use by
	   the cluster finder and the stream can be released.
	 */
	virtual void NoClustersFound(AliHLTMUONCoreClusterFinder* clusterfinder) = 0;

};


class AliHLTMUONCoreClusterFinder
{
public:

	AliHLTMUONCoreClusterFinder() : fCallback(NULL) {};


	AliHLTMUONCoreClusterFinder(const AliHLTMUONCoreClusterFinder& cf)
		: fCallback(cf.fCallback)
	{};

	AliHLTMUONCoreClusterFinder& operator = (const AliHLTMUONCoreClusterFinder& cf)
	{
		fCallback = cf.fCallback;
		return *this;
	};


	virtual ~AliHLTMUONCoreClusterFinder() {};


	/* This is the starting point of the cluster finding algorithm.
	   Deriving cluster finders should implement all processing in this method
	   to find clusters in the specified ADC stream. When all clusters are found
	   the FoundClusters method should be called to indicate that processing is
	   complete. If no clusters could be found then call NoClustersFound instead.
	 */
	virtual void FindClusters(const AliHLTMUONCoreADCStream* stream) = 0;

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
	virtual UInt FillClusterData(AliHLTMUONCoreClusterPoint* clusters, UInt arraysize) = 0;

	/* This is called when the cluster finder should be reset to an initial state.
	   All extra internal memory allocated during processing should be released.
	 */
	virtual void Reset() = 0;

	/* Sets the ClusterFinderCallback callback interface.
	 */
	inline void SetCallback(AliHLTMUONCoreClusterFinderCallback* callback)
	{
		fCallback = callback;
	};

protected:

	/* Called by the cluster finder algorithm when all clusters in the ADC stream
	   were found. At this point the ADC stream should be considered released and
	   the memory block MUST NOT be accessed.
	   The numberfound parameter should indicate how many clusters were found.
	 */
	inline void FoundClusters(UInt numberfound)
	{
		Assert(fCallback != NULL);
		fCallback->FoundClusters(this, numberfound);
	};

	/* When the cluster finding algorithm is finished processing but no clusters
	   could be found then this method should be called. At this point the ADC stream
	   should be considered released and the memory block MUST NOT be accessed.
	 */
	inline void NoClustersFound()
	{
		Assert(fCallback != NULL);
		fCallback->NoClustersFound(this);
	};

private:

	AliHLTMUONCoreClusterFinderCallback* fCallback;
};


#endif  // ALIHLTMUONCORECLUSTERFINDER_H
