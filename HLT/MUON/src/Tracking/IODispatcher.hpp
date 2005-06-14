////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_TRACKING_IO_DISPATCHER_HPP
#define dHLT_TRACKING_IO_DISPATCHER_HPP

#include "Tracking/Tracker.hpp"
#include "Tracking/IOInterface.hpp"
#include "Error.hpp"

#include <vector>
#include <algorithm>


namespace dHLT
{
namespace Tracking
{


template <typename KeyType, typename DataType>
class LookupTable
{
public:

	DataType* Add(const KeyType key)
	{
		for (UInt i = 0; i < list.size(); i++)
		{
			if ( list[i].key == key )
				return &list[i].data;
		};
		list.push_back( Record(key) );
		return &list.end()->data;
	};


	void Remove(const KeyType key)
	{
		for (UInt i = 0; i < list.size(); i++)
		{
			if ( list[i].key == key )
			{
				list.erase( &list[i] );
				return;
			};
		};
	};


	DataType* Find(const KeyType key)
	{
		for (UInt i = 0; i < list.size(); i++)
		{
			if ( list[i].key == key )
				return &list[i].data;
		};
		return NULL;
	};


	void Clear()
	{
		list.erase(list.begin(), list.end());
	};


private:

	struct Record
	{
		KeyType key;
		DataType data;

		Record() {};

		Record(KeyType k)
		{
			key = k;
		};
	};

	std::vector<Record> list;
};


template <typename DataType>
class List
{
private:

	typedef std::vector<DataType> ListType;
	ListType list;

public:

	class Iterator
	{
	public:

		Iterator() {};

		Iterator(const ListType::iterator& iter)
		{
			it = iter;
		};

		DataType& operator * () const
		{
			return *it;
		};

		DataType* operator -> () const
		{
			return it;
		};

		Iterator& operator ++ ()
		{
			it++;
			return *this;
		};

		Iterator operator ++ (int)
		{
			Iterator copy = *this;
			it++;
			return copy;
		};

		Iterator& operator -- ()
		{
			it--;
			return *this;
		};

		Iterator operator -- (int)
		{
			Iterator copy = *this;
			it--;
			return copy;
		};


		friend bool operator < (const Iterator& a, const Iterator& b)
		{
			return a.it < b.it;
		};

		friend bool operator <= (const Iterator& a, const Iterator& b)
		{
			return a.it <= b.it;
		};

		friend bool operator == (const Iterator& a, const Iterator& b)
		{
			return a.it == b.it;
		};

		friend bool operator != (const Iterator& a, const Iterator& b)
		{
			return a.it != b.it;
		};

		friend bool operator >= (const Iterator& a, const Iterator& b)
		{
			return a.it >= b.it;
		};

		friend bool operator > (const Iterator& a, const Iterator& b)
		{
			return a.it > b.it;
		};


		operator DataType () const
		{
			return *it;
		};

		operator DataType& ()
		{
			return *it;
		};


	private:

		ListType::iterator it;
	};


	void Add(const DataType data)
	{
		list.push_back(data);
	};


	void Remove(const DataType data)
	{
		for (UInt i = 0; i < list.size(); i++)
		{
			if ( list[i] == data )
			{
				list.erase( &list[i] );
				return;
			};
		};
	};


	DataType* Find(const DataType data)
	{
		for (UInt i = 0; i < list.size(); i++)
		{
			if ( list[i] == data )
				return &list[i];
		};
		return NULL;
	};


	bool Contains(const DataType data)
	{
		return Find(data) != NULL;
	};


	void Clear()
	{
		list.erase(list.begin(), list.end());
	};


	UInt Size()
	{
		return list.size();
	};


	Iterator First()
	{
		return list.begin();
	};


	Iterator Last()
	{
		return list.end();
	};

};



class ClusterBlockList
{
private:
	
	struct Record
	{
		ClusterPoint* clusters;
		UInt count;

		Record() {};

		Record(ClusterPoint* clusters, UInt count)
		{
			this->clusters = clusters;
			this->count = count;
		};

		friend bool operator == (const Record& a, const Record& b)
		{
			return a.clusters == b.clusters and a.count == b.count;
		};
	};

	List<Record> list;

public:

	typedef List<Record>::Iterator Iterator;

	void Add(ClusterPoint* clusters, UInt count)
	{
		list.Add( Record(clusters, count) );
	};

	bool Contains(ClusterPoint* clusters, UInt count)
	{
		return list.Contains( Record(clusters, count) );
	};

	void Clear()
	{
		list.Clear();
	};

	UInt Size()
	{
		return list.Size();
	};

	Iterator First()
	{
		return list.First();
	};

	Iterator Last()
	{
		return list.Last();
	};

};


class TrackerList : public List<Tracker*> {};


class ROILookupTable
{
public:

	struct Record
	{		
		ClusterBlockList clusterblocks;
		TrackerList trackers;
	};


	Record* Add(ROI region)
	{
		return table.Add(region);
	};

	void Remove(ROI region)
	{
		table.Remove(region);
	};

	Record* Find(ROI region)
	{
		return table.Find(region);
	};

	void Clear()
	{
		table.Clear();
	};

private:

	LookupTable<ROI, Record> table;
};


class EventHandler : public TrackerCallback
{
public:

	virtual void RequestClusters(const ROI regions[14], const UInt count)
	{
	};

	virtual void EndOfClusterRequests()
	{
	};

	virtual void FoundTrack()
	{
	};


	EventID Event() const
	{
		return eventid;
	};


	ROI RequestClusters(const Float left, const Float right, const Float bottom, const Float top, const UChar chamber) const
	{
		Assert( chamber < 10 );
		Assert( left <= right );
		Assert( bottom <= top );

		ROI regions[14];
		UInt regionscount;

		RegionOfInterest roi(left, right, bottom, top, chamber);
		regionscount = roi.EncodeROITree(regions);
		//callback->RequestClusters(regions, regionscount);

		Assert( regionscount > 0 );
		return regions[regionscount-1];
	};



private:

	EventID eventid;

	ROILookupTable roitable;
};


class EventLookupTable : public LookupTable<EventID, EventHandler*> {};



class TrackerStack
{
public:

	Tracker* NewTracker();
	void FreeTracker(Tracker* tracker);
};



class IOHandler : public TrackerCallback, public IOInterface
{
public:

private:

	EventLookupTable eventtable;

};


} // Tracking
} // dHLT

#endif // dHLT_TRACKING_IO_DISPATCHER_HPP
