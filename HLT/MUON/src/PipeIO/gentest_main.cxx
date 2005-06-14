////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

// This program tests the debug generator and large parts of the framework.


#include <iostream>

using std::endl;
using std::cout;


#include "Tracking/IOInterface.hpp"
#include "Tracking/Tracker.hpp"
#include "Clustering/ClusterFinder.hpp"
#include "Decision/DecisionMaker.hpp"
#include "Decision/IOInterface.hpp"
#include "Utils.hpp"
#include "Debug/ClusterSource.hpp"
#include "System/Thread.hpp"
#include "System/Mutex.hpp"
#include "System/MutexCondition.hpp"
#include "Buffers/Queue.hpp"
#include "Tracking/IOHandler.hpp"
#include "Tracking/MansoTracker.hpp"
#include "Buffers/RegionTree.hpp"

#include "Framework/Global.hpp"
#include "Debug/DataGenerator.hpp"


namespace dHLT
{


class MyTrackerStack : public Tracking::TrackerStack
{
public:

	virtual Tracking::Tracker* NewTracker()
	{
		return new Tracking::MansoTracker();
	};
	
	virtual void FreeTracker(Tracking::Tracker* tracker)
	{
		delete tracker;
	};
};

MyTrackerStack trackerstack;


};


int main()
{
	using namespace dHLT;
	
	Framework::Global framework;

	Debug::TriggerSource triggers;
	Debug::ClusterSource clusters;
	Debug::DataGenerator generator;
	
	EventID id(1, 1);
	
	generator.GenerateData(id, 1, triggers, clusters);
	
	/*
	DebugCode(
	
	triggers.Dump();
	clusters.Dump();
	
	cout << "////////////////////////////////////////////////////////" << endl;
	framework.Dump();
	DDL::TriggerInputCallback* a = framework.AddTriggerInput();
	cout << "////////////////////////////////////////////////////////" << endl;
	framework.Dump();
	DDL::TriggerInputCallback* b = framework.AddTriggerInput();
	cout << "////////////////////////////////////////////////////////" << endl;
	framework.Dump();
	DDL::TriggerInputCallback* c = framework.AddTriggerInput();
	cout << "////////////////////////////////////////////////////////" << endl;
	framework.Dump();
	framework.RemoveTriggerInput( b );
	cout << "////////////////////////////////////////////////////////" << endl;
	framework.Dump();
	framework.RemoveTriggerInput( a );
	cout << "////////////////////////////////////////////////////////" << endl;
	framework.Dump();
	framework.RemoveTriggerInput( c );
	cout << "////////////////////////////////////////////////////////" << endl;
	framework.Dump();
	
	);
	*/
	
	clusters.SetCallback( framework.AddClusterFinder(&clusters) );
	triggers.SetCallback( framework.AddTriggerInput() );
	
	
	Tracking::IOHandler iohandler;
	iohandler.SetTrackerStack(&trackerstack);
	iohandler.SetCallback( framework.AddTracker(&iohandler) );
	
	
	
	clusters.Process();
	triggers.Process();
	
	framework.RemoveClusterFinder(&clusters);
	framework.RemoveTriggerInput( triggers.GetCallback() );
	framework.RemoveTracker( &iohandler );
	
	DebugCode(
		cout << "///////////////////// These lists should be clear! ///////////////////" << endl;
		framework.Dump();
	);

	cout << "done." << endl;
	

	return 0;
};
