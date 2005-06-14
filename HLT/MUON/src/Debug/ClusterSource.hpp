////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DEBUG_CLUSTER_SOURCE_HPP
#define dHLT_DEBUG_CLUSTER_SOURCE_HPP

#include "Clustering/IOInterface.hpp"
#include <vector>

namespace dHLT
{
namespace Debug
{


class ClusterSource :
	public Clustering::IOInterface
{
public:

	ClusterSource();


	virtual ~ClusterSource();


	virtual void AddADCStream(const EventID event, const ADCStream* adcstream);

	virtual void EndOfADCStreams(const EventID event);

	
	void NewEvent(const EventID id);
	
	void NewClusterBlock(const ROI region);
	
	void AddCluster(const ClusterPoint& p);
	
	void Clear();
	
	void Process();
	
	void Dump();
	
	void SetCallback(Clustering::IOCallback* callback)
	{
		framework = callback;
	};
	
	Clustering::IOCallback* GetCallback() const
	{
		return framework;
	};


protected:

	Clustering::IOCallback* framework;

	struct DataBlock
	{
		ROI region;
		UInt count;
		std::vector<ClusterPoint> cluster;
	};
	
	struct EventBlock
	{
		EventID eventid;
		std::vector<DataBlock> block;
	};
	
	std::vector<EventBlock> event;
	
	EventBlock* currentevent;
	DataBlock* currentblock;

};


} // Debug
} // dHLT

#endif // dHLT_DEBUG_CLUSTER_SOURCE_HPP
