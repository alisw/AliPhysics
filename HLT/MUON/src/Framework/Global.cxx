////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Framework/Global.hpp"
#include "Buffers/GarbageCollector.hpp"
#include "Error.hpp"


#ifdef DEBUG
#	include "Debug/print.hpp"
#	include <iostream>
	using std::endl;
	using std::cout;
#endif // DEBUG


namespace
{

class WherePointerEquals
{
public:
	WherePointerEquals(void* pointer = NULL)
	{
		ptr = pointer;
	};

	template <typename DataType>
	bool operator () (const DataType& data)
	{
		return (void*)(&data) == ptr;
	};
	
private:
	void* ptr;
};


class WhereInterfaceEquals
{
public:
	WhereInterfaceEquals(void* interface = NULL)
	{
		interf = interface;
	};

	template <typename DataType>
	bool operator () (const DataType& data)
	{
		return data.interface == interf;
	};
	
private:
	void* interf;
};


template <class CallbackType, class ListType, class IteratorType>
bool RemoveCallbackFromList(CallbackType interface, ListType& list, IteratorType dummy)
{
	IteratorType interfacei = list.Find( WherePointerEquals(interface) );
	if ( interfacei != list.End() )
	{
		list.Remove(interfacei);
		return true;
	}
	else
		return false;
};


template <class InterfaceType, class ListType,  class IteratorType>
bool RemoveInterfaceFromList(InterfaceType interface, ListType& list, IteratorType dummy)
{
	IteratorType interfacei = list.Find( WhereInterfaceEquals(interface) );
	if ( interfacei != list.End() )
	{
		list.Remove(interfacei);
		return true;
	}
	else
		return false;
};

template <class IteratorType, class ListType>
bool RemoveFromList(IteratorType interface, ListType& list)
{
	if ( interface != list.End() )
	{
		list.Remove(interface);
		return true;
	}
	else
		return false;
};


#ifdef DEBUG

std::ostream& operator << (std::ostream& os, const dHLT::Framework::Global::TriggerInput& ti)
{
	os << (void*)(&ti) << "\t" << ti.framework;
	return os;
};

std::ostream& operator << (std::ostream& os, const dHLT::Framework::Global::ADCInput& ai)
{
	os << (void*)(&ai) << "\t" << ai.framework;
	return os;
};

std::ostream& operator << (std::ostream& os, const dHLT::Framework::Global::ClusterFinder& cf)
{
	os << (void*)(&cf) << "\t" << cf.interface << "\t" << cf.framework;
	return os;
};

std::ostream& operator << (std::ostream& os, const dHLT::Framework::Global::Tracker& tr)
{
	os << (void*)(&tr) << "\t" << tr.interface << "\t" << tr.framework;
	return os;
};

std::ostream& operator << (std::ostream& os, const dHLT::Framework::Global::DecisionMaker& dm)
{
	os << (void*)(&dm) << "\t" << dm.interface << "\t" << dm.framework;
	return os;
};

std::ostream& operator << (std::ostream& os, const dHLT::Framework::Global::TrackOutput& to)
{
	os << (void*)(&to) << "\t" << to.interface << "\t" << to.framework;
	return os;
};

std::ostream& operator << (std::ostream& os, const dHLT::Framework::Global::DecisionOutput& dout)
{
	os << (void*)(&dout) << "\t" << dout.interface << "\t" << dout.framework;
	return os;
};

#endif // DEBUG

}; // namespace


namespace dHLT
{
namespace Framework
{


TriggerRecord* Global::TriggerInput::AllocateTriggerBlock(const UInt size)
{
	DebugCode( cout << "Global::TriggerInputCallback::AllocateTriggerBlock(" << size << ")" << endl; );
	return static_cast<TriggerRecord*>( GarbageCollector::Allocate(size) );
};


void Global::TriggerInput::ReturnTriggers(const EventID event, TriggerRecord* triggers, const UInt count)
{
	Assert( framework != NULL );
	framework->ReturnTriggers(this, event, triggers, count);
};


void Global::TriggerInput::EndOfTriggers(const EventID event)
{
	Assert( framework != NULL );
	framework->EndOfTriggers(this, event);
};


////////////////////////////////////////////////////////////////////////////////

ADCStream* Global::ADCInput::AllocateADCStream(const UInt size)
{
	DebugCode( cout << "Global::ADCInput::AllocateADCStream(" << size << ")" << endl; );
	return static_cast<ADCStream*>( GarbageCollector::Allocate(size) );
};


void Global::ADCInput::ReturnADCStream(const EventID event, ADCStream* adcstream)
{
	Assert( framework != NULL );
	framework->ReturnADCStream(this, event, adcstream);
};


void Global::ADCInput::EndOfADCStreams(const EventID event)
{
	Assert( framework != NULL );
	framework->EndOfADCStreams(this, event);
};

////////////////////////////////////////////////////////////////////////////////

ClusterPoint* Global::ClusterFinder::AllocateClusterBlock(const UInt size)
{
	DebugCode( cout << "Global::ClusterFinder::AllocateClusterBlock(" << size << ")" << endl; );
	return static_cast<ClusterPoint*>( GarbageCollector::Allocate(size) );
};


void Global::ClusterFinder::ReturnClusters(const EventID event, const ROI region, ClusterPoint* clusters, const UInt count)
{
	Assert( framework != NULL );
	framework->ReturnClusters(this, event, region, clusters, count);
};


void Global::ClusterFinder::EndOfClusters(const EventID event)
{
	Assert( framework != NULL );
	framework->EndOfClusters(this, event);
};


void Global::ClusterFinder::ReleaseADCStream(const ADCStream* stream)
{
	Assert( framework != NULL );
	framework->ReleaseADCStream(this, stream);
};

////////////////////////////////////////////////////////////////////////////////

void Global::Tracker::RequestClusters(const EventID event, const ROI region)
{
	Assert( framework != NULL );
	framework->RequestClusters(this, event, region);
};


void Global::Tracker::EndOfClusterRequests(const EventID event)
{
	Assert( framework != NULL );
	framework->EndOfClusterRequests(this, event);
};


Track* Global::Tracker::AllocateTrackBlock(const UInt size)
{
	DebugCode( cout << "Global::Tracker::AllocateTrackBlock(" << size << ")" << endl; );
	return static_cast<Track*>( GarbageCollector::Allocate(size) );
};


void Global::Tracker::ReturnTracks(const EventID event, Track* newtracks, const UInt count)
{
	Assert( framework != NULL );
	framework->ReturnTracks(this, event, newtracks, count);
};


void Global::Tracker::EndOfTracks(const EventID event)
{
	Assert( framework != NULL );
	framework->EndOfTracks(this, event);
};


void Global::Tracker::ReleaseTriggers(const TriggerRecord* triggers)
{
	Assert( framework != NULL );
	framework->ReleaseTriggers(this, triggers);
};


void Global::Tracker::ReleaseClusters(const ClusterPoint* clusters)
{
	Assert( framework != NULL );
	framework->ReleaseClusters(this, clusters);
};

////////////////////////////////////////////////////////////////////////////////

DecisionRecord* Global::DecisionMaker::AllocateDecisionBlock(const UInt size)
{
	DebugCode( cout << "Global::DecisionMaker::AllocateDecisionBlock(" << size << ")" << endl; );
	return static_cast<DecisionRecord*>( GarbageCollector::Allocate(size) );
};


void Global::DecisionMaker::ReturnDecision(const EventID event, DecisionRecord* decision)
{
	Assert( framework != NULL );
	framework->ReturnDecision(this, event, decision);
};


void Global::DecisionMaker::EndOfDecisions(const EventID event)
{
	Assert( framework != NULL );
	framework->EndOfDecisions(this, event);
};


void Global::DecisionMaker::ReleaseTracks(const Track* tracks)
{
	Assert( framework != NULL );
	framework->ReleaseTracks(this, tracks);
};

////////////////////////////////////////////////////////////////////////////////

void Global::TrackOutput::ReleaseTracks(const Track* tracks)
{
	Assert( framework != NULL );
	framework->ReleaseTracks(this, tracks);
};

////////////////////////////////////////////////////////////////////////////////

void Global::DecisionOutput::ReleaseDecision(const DecisionRecord* decision)
{
	Assert( framework != NULL );
	framework->ReleaseDecision(this, decision);
};

////////////////////////////////////////////////////////////////////////////////


DDL::TriggerInputCallback* Global::AddTriggerInput()
{
	TriggerInput* newinput = triggerinputlist.Add();
	newinput->framework = this;
	return newinput;
};


bool Global::RemoveTriggerInput(DDL::TriggerInputCallback* interface)
{
	TriggerInput* input = dynamic_cast<TriggerInput*>(interface);
	return RemoveCallbackFromList( input, triggerinputlist, TriggerInputList::Iterator() );
};


DDL::ADCInputCallback* Global::AddADCInput()
{
	ADCInput* newinput = adcinputlist.Add();
	newinput->framework = this;
	return newinput;
};


bool Global::RemoveADCInput(DDL::ADCInputCallback* interface)
{
	ADCInput* input = dynamic_cast<ADCInput*>(interface);
	return RemoveCallbackFromList( input, adcinputlist, ADCInputList::Iterator() );
};


Clustering::IOCallback* Global::AddClusterFinder(Clustering::IOInterface* interface)
{
	// TODO
	ClusterFinder* cf = clusterfinderlist.Add();
	cf->framework = this;
	cf->region = INVALID_ROI;
	cf->interface = interface;
	return cf;
};


Clustering::IOCallback* Global::AddClusterFinder(Clustering::IOInterface* interface, const ROI region)
{
	// TODO
	ClusterFinder* cf = clusterfinderlist.Add();
	cf->framework = this;
	cf->region = region;
	cf->interface = interface;
	return cf;
};


Clustering::IOCallback* Global::AddClusterFinder(Clustering::IOInterface* interface, const ROI* region, const UInt count)
{
	// TODO
	ClusterFinder* cf = clusterfinderlist.Add();
	cf->framework = this;
	cf->region = INVALID_ROI;
	cf->interface = interface;
	return cf;
};


bool Global::RemoveClusterFinder(Clustering::IOInterface* interface)
{
	return RemoveInterfaceFromList( interface, clusterfinderlist, ClusterFinderList::Iterator() );
};


Tracking::IOCallback* Global::AddTracker(Tracking::IOInterface* interface)
{
	Tracker* tr = trackerlist.Add();
	tr->framework = this;
	tr->interface = interface;
	currenttracker = trackerlist.First();  // Reset the round robin.
	return tr;
};


bool Global::RemoveTracker(Tracking::IOInterface* interface)
{
	bool result = RemoveInterfaceFromList( interface, trackerlist, TrackerList::Iterator() );
	currenttracker = trackerlist.First();  // Reset the round robin.
	return result;
};


Decision::IOCallback* Global::AddDecision(Decision::IOInterface* interface)
{
	DecisionMaker* dm = decisionmakerlist.Add();
	dm->framework = this;
	dm->interface = interface;
	return dm;
};


bool Global::RemoveDecision(Decision::IOInterface* interface)
{
	return RemoveInterfaceFromList( interface, decisionmakerlist, DecisionMakerList::Iterator() );
};


DDL::TrackOutputCallback* Global::AddTrackOutput(DDL::TrackOutputInterface* interface)
{
	TrackOutput* to = trackoutputlist.Add();
	to->framework = this;
	to->interface = interface;
	return to;
};


bool Global::RemoveTrackOutput(DDL::TrackOutputInterface* interface)
{
	return RemoveInterfaceFromList( interface, trackoutputlist, TrackOutputList::Iterator() );
};


DDL::DecisionOutputCallback* Global::AddDecisionOutput(DDL::DecisionOutputInterface* interface)
{
	DecisionOutput* dout = decisionoutputlist.Add();
	dout->framework = this;
	dout->interface = interface;
	return dout;
};


bool Global::RemoveDecisionOutput(DDL::DecisionOutputInterface* interface)
{
	return RemoveInterfaceFromList( interface, decisionoutputlist, DecisionOutputList::Iterator() );
};


#ifdef DEBUG

void Global::Dump()
{
	cout << "============== triggerinputlist ================" << endl;
	triggerinputlist.Dump();
	cout << "================ adcinputlist ==================" << endl;
	adcinputlist.Dump();
	cout << "============= clusterfinderlist ================" << endl;
	clusterfinderlist.Dump();
	cout << "================ trackerlist ===================" << endl;
	trackerlist.Dump();
	cout << "============== decisionmakerlist ===============" << endl;
	decisionmakerlist.Dump();
	cout << "=============== trackoutputlist ================" << endl;
	trackoutputlist.Dump();
	cout << "============= decisionoutputlist ===============" << endl;
	decisionoutputlist.Dump();
};

#endif // DEBUG

//======================= DDL::TriggerInputCallback ============================

void Global::ReturnTriggers(const TriggerInput* sender, const EventID event, TriggerRecord* triggers, const UInt count)
{
	DebugCode( cout << "Global::ReturnTriggers(" << sender << ", " << event << ", " << triggers << ", " << count << ")" << endl; );

	// Trigger records are passed on to the trackers in a round robin manner.
	if (currenttracker == trackerlist.End() )
		// Loop to the first tracker if at the end of the list.
		currenttracker = trackerlist.First();
	if (currenttracker != trackerlist.End() )   // Must check, we might have a empty list.
	{
		Tracking::IOInterface* tracker = currenttracker->interface;
		tracker->AddTriggers(event, triggers, count);
	};
	
	// The sender should be done with the trigger block so free it.
	GarbageCollector::Free(triggers);
};


void Global::EndOfTriggers(const TriggerInput* sender, const EventID event)
{
	DebugCode( cout << "Global::EndOfTriggers(" << sender << ", " << event << ")" << endl; );
};


//========================= DDL::ADCInputCallback ==============================

void Global::ReturnADCStream(const ADCInput* sender, const EventID event, ADCStream* adcstream)
{
	DebugCode( cout << "Global::ReturnADCStream(" << sender << ", " << event << ", " << adcstream << ")" << endl; );
	
	// TODO
	
	// The sender should be done with the ADC stream so free it.
	GarbageCollector::Free(adcstream);
};


void Global::EndOfADCStreams(const ADCInput* sender, const EventID event)
{
	// TODO
	
	DebugCode( cout << "Global::EndOfADCStreams(" << sender << ", " << event << ")" << endl; );
};


//======================== Clustering::IOCallback ==============================

void Global::ReturnClusters(const ClusterFinder* sender, const EventID event, const ROI region, ClusterPoint* clusters, const UInt count)
{
	DebugCode(
		cout << "Global::ReturnClusters("
		     << sender << ", "
		     << event << ", " 
		     << region << ", "
		     << clusters << ", "
		     << count << ")" << endl;
	);
	
	bool clusters_were_requested = clustertable.AddClusters(event, region, clusters, count);
	if (clusters_were_requested)
	{
		//tr
	};
	
	// The sender should be done with the cluster block so free it.
	GarbageCollector::Free(clusters);
};


void Global::EndOfClusters(const ClusterFinder* sender, const EventID event)
{
	DebugCode( cout << "Global::EndOfClusters(" << sender << ", " << event << ")" << endl; );
};


void Global::ReleaseADCStream(const ClusterFinder* sender, const ADCStream* stream)
{
	DebugCode( cout << "Global::ReleaseADCStream(" << sender << ", " << stream << ")" << endl; );
	GarbageCollector::Free( const_cast<ADCStream*>(stream) );
};


//========================= Tracking::IOCallback ===============================

void Global::RequestClusters(const Tracker* sender, const EventID event, const ROI region)
{
	DebugCode( cout << "Global::RequestClusters(" << sender << ", " << event << ", " << region << ")" << endl; );
};


void Global::EndOfClusterRequests(const Tracker* sender, const EventID event)
{
	DebugCode( cout << "Global::EndOfClusterRequests(" << sender << ", " << event << ")" << endl; );
};


void Global::ReturnTracks(const Tracker* sender, const EventID event, Track* newtracks, const UInt count)
{
	DebugCode( cout << "Global::ReturnTracks(" << sender << ", " << event << ", " << newtracks << ", " << count << ")" << endl; );
	// The sender should be done with the tracks block so free it.
	GarbageCollector::Free(newtracks);
};


void Global::EndOfTracks(const Tracker* sender, const EventID event)
{
	DebugCode( cout << "Global::EndOfTracks(" << sender << ", " << event << ")" << endl; );
};


void Global::ReleaseTriggers(const Tracker* sender, const TriggerRecord* triggers)
{
	DebugCode( cout << "Global::ReleaseTriggers(" << sender << ", " << triggers << ")" << endl; );
	GarbageCollector::Free( const_cast<TriggerRecord*>(triggers) );
};


void Global::ReleaseClusters(const Tracker* sender, const ClusterPoint* clusters)
{
	DebugCode( cout << "Global::ReleaseClusters(" << sender << ", " << clusters << ")" << endl; );
	GarbageCollector::Free( const_cast<ClusterPoint*>(clusters) );
};


//========================= Decision::IOCallback ===============================

void Global::ReturnDecision(const DecisionMaker* sender, const EventID event, DecisionRecord* decision)
{
	DebugCode( cout << "Global::ReturnDecision(" << sender << ", " << event << ", " << decision << ")" << endl; );
	// The sender should be done with the decision block so free it.
	GarbageCollector::Free(decision);
};


void Global::EndOfDecisions(const DecisionMaker* sender, const EventID event)
{
	DebugCode( cout << "Global::EndOfDecisions(" << sender << ", " << event << ")" << endl; );
};


void Global::ReleaseTracks(const DecisionMaker* sender, const Track* tracks)
{
	DebugCode( cout << "Global::ReleaseTracks(" << sender << ", " << tracks << ")" << endl; );
	GarbageCollector::Free( const_cast<Track*>(tracks) );
};


//======================= DDL::TrackOutputCallback =============================

void Global::ReleaseTracks(const TrackOutput* sender, const Track* tracks)
{
	DebugCode( cout << "Global::ReleaseTracks(" << sender << ", " << tracks << ")" << endl; );
	GarbageCollector::Free( const_cast<Track*>(tracks) );
};


//====================== DDL::DecisionOutputCallback ===========================

void Global::ReleaseDecision(const DecisionOutput* sender, const DecisionRecord* decision)
{
	DebugCode( cout << "Global::ReleaseDecision(" << sender << ", " << decision << ")" << endl; );
	GarbageCollector::Free( const_cast<DecisionRecord*>(decision) );
};


} // Framework
} // dHLT
