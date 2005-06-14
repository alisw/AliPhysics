////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_FRAMEWORK_CLUSTER_LOOKUP_TABLE_HPP
#define dHLT_FRAMEWORK_CLUSTER_LOOKUP_TABLE_HPP

/*
#include <vector>
#include <algorithm>
*/

#include "Buffers/RegionTree.hpp"
#include "Buffers/List.hpp"
#include "Buffers/LookupTable.hpp"

namespace dHLT
{
namespace Framework
{


using dHLT::Buffers::List;
using dHLT::Buffers::LookupTable;


class ClusterLookupTable
{
public:
	
	struct ClusterBlockRecord
	{
		ROI region;
		ClusterPoint* clusters;
		UInt count;

		bool operator == (const ClusterBlockRecord& rec) const
		{
			return region == rec.region and clusters == rec.clusters and count == rec.count;
		};
	};

	typedef List<ClusterBlockRecord> ClusterBlockList;


	/* Adds the specified cluster block to the lookup table.
	   True is returned if it was found that the cluster block just added has
	   already been requested.
	 */
	bool AddClusters(const EventID& event, const ROI region, const ClusterPoint* clusters, const UInt count)
	{
		EventRecord* eventrec = FetchEvent(event);

		ChamberID chamber; UChar level; UInt left, bottom;
		RegionOfInterest::Decode(region, chamber, level, left, bottom);

		AddClustersHandler addblock(region, clusters, count);
		eventrec->tree.ForAllOverlapping(chamber, left, bottom, level, addblock);

		return addblock.requested;
	};
	

	/* Adds the a request to the lookup table for a specific region.
	   The foundlist should be initialised to empty before calling this method.
	   It will be filled with all the cluster blocks already in the lookup table
	   that the caller is interested in.
	   True is returned if no more clusters can be expected for the specified region.
	 */
	bool AddRequest(const EventID& event, const ROI region, ClusterBlockList& foundlist)
	{
		EventRecord* eventrec = FetchEvent(event);

		ChamberID chamber; UChar level; UInt left, bottom;
		RegionOfInterest::Decode(region, chamber, level, left, bottom);

		AddRequestHandler addrequest(&foundlist);
		eventrec->tree.ForAllOverlapping(chamber, left, bottom, level, addrequest);

		return addrequest.endofclusters;
	};


	/* Marks the end of clusters for a whole chamber.
	 */
	bool MarkEndOfClusters(const EventID& event, const ChamberID chamber)
	{
		EventRecord* eventrec = FetchEvent(event);

		MarkEndOfClustersHandler markendofclusters;
		eventrec->tree.ForAllOverlapping(chamber, 0, 0, 0, markendofclusters);

		return markendofclusters.requested;
	};


	/* Sets the end of clusters flag for all nodes in the specified region.
	   All traversed nodes are closed so no more clusters can be added to them.
	   Children nodes, whose parent has end of clusters flag set, are also closed.
	   True is returned if any of the traversed nodes was requested.
	 */
	bool MarkEndOfClusters(const EventID& event, const ROI region)
	{
		EventRecord* eventrec = FetchEvent(event);

		ChamberID chamber; UChar level; UInt left, bottom;
		RegionOfInterest::Decode(region, chamber, level, left, bottom);

		MarkEndOfClustersHandler markendofclusters;
		eventrec->tree.ForAllOverlapping(chamber, left, bottom, level, markendofclusters);

		return markendofclusters.requested;
	};
	
	
	/* Removes the specified event from the lookup table. All internal data
	   structures are deleted. This method returns true if the event was actually
	   found to be deleted and false if no physical deletion took place.
	 */
	bool RemoveEvent(const EventID& event)
	{
		EventRecord* eventrec = eventtable.Find(event);
		if (eventrec != NULL)
		{
			eventtable.Remove(event);
			return true;
		}
		else
			return false;
	};


private:

	struct TreeRecord
	{
		bool requested;   // Indicated this clusters were requested for this sub-region.
		bool endofclusters;  // Indicated that the end of clusters signal was received for this sub-region.
		bool closed;  // Indicated that no more clusters should be added to this node.
		ClusterBlockList clusterblocks;

		TreeRecord()
		{
			requested = endofclusters = closed = false;
		};
	};

	typedef Buffers::RegionTree<TreeRecord> Tree;

	struct EventRecord
	{
		EventID event;
		Tree tree;

		EventRecord(const EventID& id)
		{
			event = id;
		};

		EventID Key() const { return event; };
	};

	typedef LookupTable<EventID, EventRecord> EventLookupTable;


	class AddClustersHandler : public Tree::Handler
	{
	public:

		inline AddClustersHandler(const ROI region, const ClusterPoint* clusters, const UInt count)
		{
			data.region = region;
			data.clusters = const_cast<ClusterPoint*>( clusters );
			data.count = count;
			requested = false;
		};

		inline void AddTo(Tree::Node* node)
		{
			if (node->data.closed)
			{
				// Received a previous signal indicating no more clusters are
				// to be expected in this sub region, so this is an error.
				// TODO: throw exception.
			};
			
			// Add the new cluster block to the appropriate node.
			node->data.clusterblocks.Add(data);
		};

		inline void Initialise(Tree::Node* node, Tree::Node* parent)
		{
			// Need to inherit the parents end of clusters flag.
			node->data.endofclusters = parent->data.endofclusters;
			
			// Close the node only if the parent got the end of clusters signal.
			node->data.closed = parent->data.endofclusters;
		};

		inline void FoundOverlapping(Tree::Node* node)
		{
			// Check if anyone has requested this new cluster block in an
			// overlapping region of interest.
			requested |= node->data.requested;
		};

		inline void FoundContained(Tree::Node* node)
		{
			// Check if anyone has requested this new cluster block in a
			// contained (child) region of interest.
			requested |= node->data.requested;
		};

		ClusterBlockRecord data;
		bool requested;
	};


	class AddRequestHandler : public Tree::Handler
	{
	public:

		inline AddRequestHandler(ClusterBlockList* output)
		{
			Assert( output != NULL );
			foundlist = output;

			// since we do an 'and' operation on endofclusters, it must be 
			// initialised to true.
			endofclusters = true;
		};

		inline void AddTo(Tree::Node* node)
		{
			// Check to see if there are no more clusters in the region.
			// Remember that all 4 tree nodes that form the region must have
			// their endofclusters flag set, so we use the 'and' operation.
			endofclusters &= node->data.endofclusters;
			
			// Mark that we are interested in this sub region.
			node->data.requested = true;
		};

		inline void Initialise(Tree::Node* node, Tree::Node* parent)
		{
			// Need to inherit the parents end of clusters flag.
			node->data.endofclusters = parent->data.endofclusters;
			
			// Close the node only if the parent got the end of clusters signal.
			node->data.closed = parent->data.endofclusters;
		};

		inline void FoundOverlapping(Tree::Node* node)
		{
			// For all cluster blocks already in the tree we need to add them
			// to the lis of found blocks.
			for (ClusterBlockList::Iterator rec = node->data.clusterblocks.First(); 
				rec != node->data.clusterblocks.End(); 
				rec++
				)
			{
				foundlist->AddUniquely(*rec);
			};
		};

		inline void FoundContained(Tree::Node* node)
		{
			// For all cluster blocks already in the tree we need to add them
			// to the lis of found blocks.
			for (ClusterBlockList::Iterator rec = node->data.clusterblocks.First(); 
				rec != node->data.clusterblocks.End(); 
				rec++
				)
			{
				foundlist->AddUniquely(*rec);
			};
		};

		ClusterBlockList* foundlist;
		bool endofclusters;
	};
	

	class MarkEndOfClustersHandler : public Tree::Handler
	{
	public:

		inline MarkEndOfClustersHandler()
		{
			requested = false;
		};

		inline void AddTo(Tree::Node* node)
		{
			//node->data.endofclusters = true;  // Performed in Initialise or FoundOverlapping.
			
			// All traversed nodes must be closed since they can no longer receive clusters.
			node->data.closed = true;
		};

		inline void InitialiseFirst(Tree::Node* node)
		{
			// All nodes overlapping the region must have the end of clusters flag set.
			node->data.endofclusters = true;
			
			// All traversed nodes must be closed since they can no longer receive clusters.
			node->data.closed = true;
		};

		inline void Initialise(Tree::Node* node, Tree::Node* parent)
		{
			// All nodes overlapping the region must have the end of clusters flag set.
			node->data.endofclusters = true;
			
			// All traversed nodes must be closed since they can no longer receive clusters.
			node->data.closed = true;
		};

		inline void FoundOverlapping(Tree::Node* node)
		{
			// Check if any requests were made to this node.
			requested |= node->data.requested;
			
			// All nodes overlapping the region must have the end of clusters flag set.
			node->data.endofclusters = true;
			
			// All traversed nodes must be closed since they can no longer receive clusters.
			node->data.closed = true;
		};

		inline void FoundContained(Tree::Node* node)
		{
			// Check if any requests were made to this node.
			requested |= node->data.requested;
			
			// All nodes contained in the region must have the end of clusters flag set.
			node->data.endofclusters = true;
			
			// All traversed nodes must be closed since they can no longer receive clusters.
			node->data.closed = true;
		};

		bool requested;
	};
	
	
	/* Checks to see if the event exists. If it does, the EventRecord is returend 
	   otherwise a new event is created and returned.
	 */
	EventRecord* FetchEvent(const EventID& event)
	{
		EventRecord* eventrec = eventtable.Find(event);
		if (eventrec == NULL)
			eventrec = eventtable.Add(event);
		return eventrec;
	};


	EventLookupTable eventtable;
};


} // Framework
} // dHLT

#endif // dHLT_FRAMEWORK_CLUSTER_LOOKUP_TABLE_HPP
