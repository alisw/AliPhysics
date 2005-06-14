////////////////////////////////////////////////////////////////////////////////
//
// Author: Gareth de Vaux
// Email:  devaux@lhc.phy.uct.ac.za | dhlt@lordcow.org
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "PubSub/TrackerCallback.hpp"
#include "Cluster.hpp"


namespace dHLT
{
namespace PubSub
{


TrackerCallback::TrackerCallback()
{
	for (UInt i=0; i<10; i++)
		clusterblocks[i].erase(clusterblocks[i].begin(), clusterblocks[i].end());

	tracks = NULL;
	trackcount = 0;

	cyclic_count = 1; // TEMP
}

void TrackerCallback::RequestClusters(Tracking::Tracker* tracker, const Float /*left*/, const Float /*right*/, const Float /*bottom*/, const Float /*top*/, const ChamberID chamber, const void* tag)
{
	Assert(0 <= chamber and chamber <= 9);

	UInt i;

	for (i=0; i<clusterblocks[chamber].size(); i++)
	{
		ClusterBlock &block = clusterblocks[chamber][i];

		if (block.data != NULL)
			tracker->ReturnClusters(const_cast<void*>(tag), block.data, block.count);
		else
			DebugMsg(1, "null cluster pointer");

	}

	tracker->EndOfClusters(const_cast<void*>(tag)); // put this in the else block?
}

void TrackerCallback::EndOfClusterRequests(Tracking::Tracker* /*tracker*/)
{
	DebugMsg(1, "End of cluster requests");
}

void TrackerCallback::FoundTrack(Tracking::Tracker* tracker)
{
	//tracks[] must get set by SetTracks() before this happens:

	Track& newtrack = tracks[trackcount++];
	tracker->FillTrackData(newtrack);

	newtrack.triggerid = cyclic_count++;
}
	
void TrackerCallback::NoTrackFound(Tracking::Tracker* /*tracker*/)
{
	DebugMsg(1, "No tracks found");
}

void TrackerCallback::AddClusterBlock(const ClusterPoint *clusters, const UInt size, const UInt chamber)
{
	Assert(chamber<10);

	ClusterBlock newrec;
	newrec.data = clusters;
	newrec.count = size / sizeof(ClusterPoint);
	clusterblocks[chamber].push_back(newrec);
}


} // PubSub
} // dHLT
