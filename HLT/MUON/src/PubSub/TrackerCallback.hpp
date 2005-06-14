////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// Author: Gareth de Vaux
// Email:  devaux@lhc.phy.uct.ac.za | dhlt@lordcow.org
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_PUBSUB_TRACKERCALLBACK_HPP
#define dHLT_PUBSUB_TRACKERCALLBACK_HPP

#include "Tracking/Tracker.hpp"  
#include "Tracking/MansoTracker.hpp" 
#include "Cluster.hpp"

#include <vector>

namespace dHLT
{
namespace PubSub
{


class TrackerCallback: public Tracking::TrackerCallback
{
public:

	TrackerCallback();

	virtual ~TrackerCallback() {}
	
	/* All clusters that fall within the specified boundary box on the specified
	   chamber should be returned the the tracker, by calling the ReturnClusters
	   method of the given tracker. The same tag parameter must be passed on the 
	   ReturnClusters method's parameter list.
	 */
	virtual void RequestClusters(
			Tracking::Tracker* tracker, 
			const Float left, const Float right, const Float bottom, const Float top,
			const ChamberID chamber, const void* tag
		);

	/* When this method is called then one knows no more RequestClusters method
	   calls are expected.
	 */
	virtual void EndOfClusterRequests(Tracking::Tracker* tracker);

	/* This method is called when the tracker has found a track. The FillTrackData
	   method of the given tracker should be called to received the track data.
	   At this point all cluster blocks can be released.
	 */
	virtual void FoundTrack(Tracking::Tracker* tracker);
	
	/* When the tracker is finished with its work but no track was found then
	   this method is called. At this point no more work should be performed by
	   the tracker and all cluster blocks can be released.
	 */
	virtual void NoTrackFound(Tracking::Tracker* tracker);

	/* Store a pointer to cluster data
	*/
	void AddClusterBlock(const ClusterPoint *clusters, const UInt size, const UInt chamber);

	/* gives tracks access to the externally allocated output buffer
	*/
	void SetTracks(Track* output) {tracks = output;}

	UInt TrackCount() const {return trackcount;}

private:

	struct ClusterBlock
	{
		const ClusterPoint* data;
		UInt count;
	};

	std::vector<ClusterBlock> clusterblocks[10];  // cluster blocks for the 10 tracking chambers

	Track* tracks;              // array of tracks
	UInt trackcount;

	UInt cyclic_count;          // TODO
};


} // PubSub
} // dHLT

#endif // dHLT_TRACKING_TRACKERCALLBACK_HPP
