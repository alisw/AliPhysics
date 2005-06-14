////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_TRACKING_EVENT_HANDLER_HPP
#define dHLT_TRACKING_EVENT_HANDLER_HPP

#include <vector>
#include <algorithm>

#include "Buffers/RegionTree.hpp"
#include "Buffers/List.hpp"
#include "Buffers/CountedList.hpp"
#include "Buffers/LookupTable.hpp"

#include "IOInterface.hpp"
#include "Tracking/Tracker.hpp"
#include "Error.hpp"


/*
ostream& operator << (ostream& os, const dHLT::EventID& id)
{
	os << "[" << id.bunch << ", " << id.timestamp << "]";
	return os;
};
*/


namespace dHLT
{
namespace Tracking
{


using dHLT::Buffers::List;
using dHLT::Buffers::CountedList;
using dHLT::Buffers::LookupTable;


class IOHandler;



/////////////////////////////////////////////////////////////////////////////////////////////////

class RegionRegistrationTable
{
public:

	struct ClusterBlockRecord
	{
		ClusterPoint* clusters;
		UInt count;

		ClusterBlockRecord() {};

		ClusterBlockRecord(ClusterPoint* clusters, UInt count)
		{
			this->clusters = clusters;
			this->count = count;
		};

		friend bool operator == (const ClusterBlockRecord& a, const ClusterBlockRecord& b)
		{
			return a.clusters == b.clusters and a.count == b.count;
		};
	};
	
	typedef List<ClusterBlockRecord> ClusterBlockList;


	struct TrackerRecord
	{
		Tracker* tracker;
		void* tag;

		TrackerRecord(Tracker* newtracker = NULL, void* newtag = NULL)
		{
			tracker = newtracker;
			tag = newtag;
		};

		friend bool operator == (const TrackerRecord& a, const TrackerRecord& b)
		{
			return a.tracker == b.tracker and a.tag == b.tag;
		};

		friend bool operator != (const TrackerRecord& a, const TrackerRecord& b)
		{
			return not (a == b);
		};
	};

	typedef List<TrackerRecord> TrackerRecordList;


	bool AddRequest(
			Tracker* tracker, const void* tag,
			const ChamberID chamber, const UInt left, const UInt bottom, const UChar level,
			ClusterBlockList& foundlist
		)
	{
		AddRequestHandler addrequest(tracker, tag, &foundlist);
		regions.ForAllOverlapping(chamber, left, bottom, level, addrequest);
		return addrequest.endofclusters;
	};
	

	void AddClusters(
			const ClusterPoint* clusters, const UInt count,
			const ChamberID chamber, const UInt left, const UInt bottom, const UChar level,
			TrackerRecordList& foundlist
		)
	{
		AddClustersHandler addblock(clusters, count, &foundlist);
		regions.ForAllOverlapping(chamber, left, bottom, level, addblock);
	};


	void MarkEndOfClusters(
			const ChamberID chamber, const UInt left, const UInt bottom, const UChar level,
			TrackerRecordList& foundlist
		)
	{
		MarkEndOfClustersHandler markendofclusters(&foundlist);
		regions.ForAllOverlapping(chamber, left, bottom, level, markendofclusters);
	};


private:

	struct TreeRecord
	{
		bool endofclusters;
		bool forchildren;     // Set to true if end of clusters flag applies for all children.
		ClusterBlockList clusterblocks;  // The cluster blocks returned so far.
		TrackerRecordList trackers;      // Trackers interested in the region.

		TreeRecord()
		{
			endofclusters = forchildren = false;
		};
	};

	typedef Buffers::RegionTree<TreeRecord> Tree;


	class AddRequestHandler : public Tree::Handler
	{
	public:

		inline AddRequestHandler(Tracker* tracker, const void* tag, ClusterBlockList* output)
		{
			Assert( output != NULL );
			data.tracker = tracker;
			data.tag =  const_cast<void*>( tag );
			foundlist = output;

			// since we do an and operation on endofclusters, it must be initialised to true.
			endofclusters = true;
		};

		inline void AddTo(Tree::Node* node)
		{
			endofclusters = endofclusters & node->data.endofclusters;
			node->data.trackers.Add(data);
		};

		inline void Initialise(Tree::Node* node, Tree::Node* parent)
		{
			node->data.endofclusters = parent->data.endofclusters & parent->data.forchildren;
		};

		inline void FoundOverlapping(Tree::Node* node)
		{
			for (ClusterBlockList::Iterator block = node->data.clusterblocks.First();
				block != node->data.clusterblocks.End(); 
				block++
				)
			{
				foundlist->AddUniquely(*block);
			};
		};

		inline void FoundContained(Tree::Node* node)
		{
			for (ClusterBlockList::Iterator block = node->data.clusterblocks.First();
				block != node->data.clusterblocks.End(); 
				block++
				)
			{
				foundlist->AddUniquely(*block);
			};
		};

		TrackerRecord data;
		ClusterBlockList* foundlist;
		bool endofclusters;
	};


	class AddClustersHandler : public Tree::Handler
	{
	public:

		inline AddClustersHandler(const ClusterPoint* clusters, const UInt count, TrackerRecordList* output)
		{
			Assert( output != NULL );
			data.clusters = const_cast<ClusterPoint*>( clusters );
			data.count = count;
			foundlist = output;
		};

		inline void AddTo(Tree::Node* node)
		{
			node->data.clusterblocks.Add(data);
		};

		inline void Initialise(Tree::Node* node, Tree::Node* parent)
		{
			node->data.endofclusters = parent->data.endofclusters & parent->data.forchildren;
		};

		inline void FoundOverlapping(Tree::Node* node)
		{
			for (TrackerRecordList::Iterator tracker = node->data.trackers.First();
				tracker != node->data.trackers.End(); 
				tracker++
				)
			{
				foundlist->AddUniquely(*tracker);
			};
		};

		inline void FoundContained(Tree::Node* node)
		{
			for (TrackerRecordList::Iterator tracker = node->data.trackers.First();
				tracker != node->data.trackers.End(); 
				tracker++
				)
			{
				foundlist->AddUniquely(*tracker);
			};
		};

		ClusterBlockRecord data;
		TrackerRecordList* foundlist;
	};


	class MarkEndOfClustersHandler : public Tree::Handler
	{
	public:

		inline MarkEndOfClustersHandler(TrackerRecordList* output)
		{
			Assert( output != NULL );
			foundlist = output;
		};

		inline void AddTo(Tree::Node* node)
		{
			node->data.forchildren = true;
		};

		inline void InitialiseFirst(Tree::Node* node)
		{
			node->data.endofclusters = true;
			node->data.forchildren = false;
		};

		inline void Initialise(Tree::Node* node, Tree::Node* parent)
		{
			node->data.endofclusters = true;
			node->data.forchildren = false;
		};

		inline void FoundOverlapping(Tree::Node* node)
		{
			node->data.endofclusters = true;
			node->data.forchildren = false;
		};

		inline void FoundContained(Tree::Node* node)
		{
			node->data.endofclusters = true;
			node->data.forchildren = true;

			for (TrackerRecordList::Iterator tracker = node->data.trackers.First();
				tracker != node->data.trackers.End(); 
				tracker++
				)
			{
				foundlist->AddUniquely(*tracker);
			};
		};

		TrackerRecordList* foundlist;
	};


	Tree regions;
};


/////////////////////////////////////////////////////////////////////////////////////////////


class EventHandler : public TrackerCallback
{
public:

	EventHandler(const EventID event)
	{
		eventid = event;
		iohandler = NULL;
		trackblock = NULL;
		maxtrackcount = currenttrackcount = newtrackcount = 0;
		endoftriggers = false;
		endofrequests_count = 0;
	};
	
	virtual ~EventHandler() {};

	void AddTriggers(const TriggerRecord* triggers, const UInt count);

	void EndOfTriggers();

	void AddClusters(const ROI region, const ClusterPoint* clusters, const UInt count);

	void EndOfClusters(const ROI region);


	// TrackerCallback methods

	virtual void RequestClusters(
			Tracker* tracker,
			const Float left, const Float right, const Float bottom, const Float top,
			const ChamberID chamber, const void* tag
		);

	virtual void EndOfClusterRequests(Tracker* tracker);

	virtual void FoundTrack(Tracker* tracker);
	
	virtual void NoTrackFound(Tracker* tracker);

	void CheckForTotalEndOfRequests();

	void CheckForEndOfTracks();

	EventID Key() const { return eventid; };
	EventID Event() const { return eventid; };

	IOHandler* GetIOHandler() const { return iohandler; };
	void SetIOHandler(IOHandler* handler) { iohandler = handler; };

private:

	typedef CountedList<Tracker*> TrackerList;
	
	/*
	class TrackerList : public List<Tracker*>
	{
	public:
		TrackerList() : List<Tracker*>()
		{
			count = 0;
		}
		
		Tracker** Add()
		{
			count++;
			return List<Tracker*>::Add();
		};

		
		Tracker** AddNew(Tracker* data)
		{
			count++;
			return List<Tracker*>::AddNew(data);
		};

		Tracker** AddUniquely(Tracker* data)
		{
			Tracker** result = Find(data);
			if (result == NULL)
			{
				count++;
				return AddNew(data);
			}
			else
				return result;
		};

		void Remove(const UInt index)
		{
			List<Tracker*>::Remove(index);
			count--;
		};

		bool Remove(Tracker* data)
		{
			bool removed = List<Tracker*>::Remove(data);
			if (removed)
				count--;
			return removed;
		};

		void Remove(Iterator& iter)
		{
			List<Tracker*>::Remove(iter);
			count--;
		};

		void Clear()
		{
			List<Tracker*>::Clear();
			count = 0;
		};
		
		UInt Count() const
		{
			return count;
		};
	
	private:
	
		UInt count;
	};
	*/


	EventID eventid;
	IOHandler* iohandler;

	Track* trackblock;
	UInt maxtrackcount, currenttrackcount, newtrackcount;

	bool endoftriggers;
	UInt endofrequests_count;

	TrackerList trackers;
	RegionRegistrationTable registration_table;
};


} // Tracking
} // dHLT

#endif // dHLT_TRACKING_EVENT_HANDLER_HPP
