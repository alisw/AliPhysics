////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Debug/ClusterSource.hpp"
#include "Debug/print.hpp"

#include <iostream>
using std::endl;
using std::cout;

namespace dHLT
{
namespace Debug
{


ClusterSource::ClusterSource()
{
	framework = NULL;
	currentevent = NULL;
	currentblock = NULL;
};


ClusterSource::~ClusterSource()
{
	Clear();
};


void ClusterSource::AddADCStream(const EventID event, const ADCStream* adcstream)
{
	// Do nothing
};


void ClusterSource::EndOfADCStreams(const EventID event)
{
	// Do nothing
};


void ClusterSource::NewEvent(const EventID id)
{
	EventBlock newevent;
	newevent.eventid = id;
	event.push_back(newevent);
	currentevent = &event[event.size() -1];
	currentblock = NULL;
};

void ClusterSource::NewClusterBlock(const ROI region)
{
	Assert( currentevent != NULL );
	DataBlock newblock;
	newblock.region = region;
	currentevent->block.push_back(newblock);
	currentblock = &currentevent->block[currentevent->block.size() -1];
};

void ClusterSource::AddCluster(const ClusterPoint& p)
{
	Assert( currentblock != NULL );
	currentblock->cluster.push_back(p);
};

void ClusterSource::Clear()
{
	event.clear();
	currentevent = NULL;
	currentblock = NULL;
};


void ClusterSource::Process()
{
	Assert( framework != NULL );
	
	for (UInt i = 0; i < event.size(); i++)
	{
		for (UInt j = 0; j < event[i].block.size(); j++)
		{
			UInt count = event[i].block[j].cluster.size();
			UInt size = count * sizeof(ClusterPoint);
			ClusterPoint* clusterblock = framework->AllocateClusterBlock(size);
			for (UInt k = 0; k < count; k++)
				clusterblock[k] = event[i].block[j].cluster[k];

			framework->ReturnClusters(
					event[i].eventid,
					event[i].block[j].region,
					clusterblock,
					count
				);
			// Note: clusterblock is deleted by framework.
		};
		framework->EndOfClusters(event[i].eventid);
	};
};


void ClusterSource::Dump()
{
	cout << "==================== Clusters ====================" << endl;
	for (UInt i = 0; i < event.size(); i++)
	{
		cout << "event = " << event[i].eventid << endl;
		for (UInt j = 0; j < event[i].block.size(); j++)
		{
			cout << "\tregion = " << event[i].block[j].region << endl << "\t\t";
			for (UInt k = 0; k < event[i].block[j].cluster.size(); k++)
			{
				cout << event[i].block[j].cluster[k] << " ";
				if (k % 5 == 4) cout << endl << "\t\t";
			};
			cout << endl;
		};
	};
};


} // Debug
} // dHLT
