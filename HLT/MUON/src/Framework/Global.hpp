////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_FRAMEWORK_GLOBAL_HPP
#define dHLT_FRAMEWORK_GLOBAL_HPP

#include "BasicTypes.hpp"
#include "Buffers/List.hpp"
#include "DDL/TriggerInputInterface.hpp"
#include "DDL/ADCInputInterface.hpp"
#include "Clustering/IOInterface.hpp"
#include "Tracking/IOInterface.hpp"
#include "Decision/IOInterface.hpp"
#include "DDL/TrackOutputInterface.hpp"
#include "DDL/DecisionOutputInterface.hpp"
#include "Framework/ClusterLookupTable.hpp"

namespace dHLT
{
namespace Framework
{


using namespace dHLT::Buffers;


class Global
{
public:

	Global()
	{
	};
	
	virtual ~Global()
	{
		// clusterfinderlist, trackerlist and decisionlist free their memory implicitly.
		// Note the objects to which these lists point to are not deleted because they 
		// are not owned by this object.
	};
	
	DDL::TriggerInputCallback* AddTriggerInput();
	bool RemoveTriggerInput(DDL::TriggerInputCallback* interface);
	
	DDL::ADCInputCallback* AddADCInput();
	bool RemoveADCInput(DDL::ADCInputCallback* interface);
	
	Clustering::IOCallback* AddClusterFinder(Clustering::IOInterface* interface);
	Clustering::IOCallback* AddClusterFinder(Clustering::IOInterface* interface, const ROI region);
	Clustering::IOCallback* AddClusterFinder(Clustering::IOInterface* interface, const ROI* region, const UInt count);
	bool RemoveClusterFinder(Clustering::IOInterface* interface);
	
	Tracking::IOCallback* AddTracker(Tracking::IOInterface* interface);
	bool RemoveTracker(Tracking::IOInterface* interface);
	
	Decision::IOCallback* AddDecision(Decision::IOInterface* interface);
	bool RemoveDecision(Decision::IOInterface* interface);
	
	DDL::TrackOutputCallback* AddTrackOutput(DDL::TrackOutputInterface* interface);
	bool RemoveTrackOutput(DDL::TrackOutputInterface* interface);
	
	DDL::DecisionOutputCallback* AddDecisionOutput(DDL::DecisionOutputInterface* interface);
	bool RemoveDecisionOutput(DDL::DecisionOutputInterface* interface);
	
	
#ifdef DEBUG
	void Dump();

public:
#else // DEBUG
private:
#endif // DEBUG

	class FrameworkCallback
	{
	public:
		FrameworkCallback()
		{
			framework = NULL;
		};
		
		friend std::ostream& operator << (std::ostream& os, const FrameworkCallback& fc)
		{
			os << (void*) fc.framework;
			return os;
		};
		
		Global* framework;
	};

	class TriggerInput : public DDL::TriggerInputCallback, public FrameworkCallback
	{
	public:
		virtual TriggerRecord* AllocateTriggerBlock(const UInt size);
		virtual void ReturnTriggers(const EventID event, TriggerRecord* triggers, const UInt count);
		virtual void EndOfTriggers(const EventID event);
	};

	class ADCInput : public DDL::ADCInputCallback, public FrameworkCallback
	{
	public:
		virtual ADCStream* AllocateADCStream(const UInt size);
		virtual void ReturnADCStream(const EventID event, ADCStream* adcstream);
		virtual void EndOfADCStreams(const EventID event);
	};
	
	class ClusterFinder : public Clustering::IOCallback, public FrameworkCallback
	{
	public:
		virtual ClusterPoint* AllocateClusterBlock(const UInt size);
		virtual void ReturnClusters(const EventID event, const ROI region, ClusterPoint* clusters, const UInt count);
		virtual void EndOfClusters(const EventID event);
		virtual void ReleaseADCStream(const ADCStream* stream);
		
		ClusterFinder() : Clustering::IOCallback(), FrameworkCallback() { region = INVALID_ROI; };
		ROI region;   // The region this cluster finder covers.
		Clustering::IOInterface* interface;
	};
	
	class Tracker : public Tracking::IOCallback, public FrameworkCallback
	{
	public:
		virtual void RequestClusters(const EventID event, const ROI region);
		virtual void EndOfClusterRequests(const EventID event);
		virtual Track* AllocateTrackBlock(const UInt size);
		virtual void ReturnTracks(const EventID event, Track* newtracks, const UInt count);
		virtual void EndOfTracks(const EventID event);
		virtual void ReleaseTriggers(const TriggerRecord* triggers);
		virtual void ReleaseClusters(const ClusterPoint* clusters);
		
		Tracking::IOInterface* interface;
	};
	
	class DecisionMaker : public Decision::IOCallback, public FrameworkCallback
	{
	public:
		virtual DecisionRecord* AllocateDecisionBlock(const UInt size);
		virtual void ReturnDecision(const EventID event, DecisionRecord* decision);
		virtual void EndOfDecisions(const EventID event);
		virtual void ReleaseTracks(const Track* tracks);
		
		Decision::IOInterface* interface;
	};
	
	class TrackOutput : public DDL::TrackOutputCallback, public FrameworkCallback
	{
	public:
		virtual void ReleaseTracks(const Track* tracks);
		
		DDL::TrackOutputInterface* interface;
	};
	
	class DecisionOutput : public DDL::DecisionOutputCallback, public FrameworkCallback
	{
	public:
		virtual void ReleaseDecision(const DecisionRecord* decision);
		
		DDL::DecisionOutputInterface* interface;
	};

	
private:
	
	//==================== DDL::TriggerInputCallback =========================
	void ReturnTriggers(const TriggerInput* sender, const EventID event, TriggerRecord* triggers, const UInt count);
	void EndOfTriggers(const TriggerInput* sender, const EventID event);
	
	//====================== DDL::ADCInputCallback ===========================
	void ReturnADCStream(const ADCInput* sender, const EventID event, ADCStream* adcstream);
	void EndOfADCStreams(const ADCInput* sender, const EventID event);
	
	//===================== Clustering::IOCallback ===========================
	void ReturnClusters(const ClusterFinder* sender, const EventID event, const ROI region, ClusterPoint* clusters, const UInt count);
	void EndOfClusters(const ClusterFinder* sender, const EventID event);
	void ReleaseADCStream(const ClusterFinder* sender, const ADCStream* stream);
	
	//====================== Tracking::IOCallback ============================
	void RequestClusters(const Tracker* sender, const EventID event, const ROI region);
	void EndOfClusterRequests(const Tracker* sender, const EventID event);
	void ReturnTracks(const Tracker* sender, const EventID event, Track* newtracks, const UInt count);
	void EndOfTracks(const Tracker* sender, const EventID event);
	void ReleaseTriggers(const Tracker* sender, const TriggerRecord* triggers);
	void ReleaseClusters(const Tracker* sender, const ClusterPoint* clusters);
	
	//====================== Decision::IOCallback ============================
	void ReturnDecision(const DecisionMaker* sender, const EventID event, DecisionRecord* decision);
	void EndOfDecisions(const DecisionMaker* sender, const EventID event);
	void ReleaseTracks(const DecisionMaker* sender, const Track* tracks);
	
	//==================== DDL::TrackOutputCallback ==========================
	void ReleaseTracks(const TrackOutput* sender, const Track* tracks);
	
	//=================== DDL::DecisionOutputCallback ========================
	void ReleaseDecision(const DecisionOutput* sender, const DecisionRecord* decision);
	
	
	typedef List<TriggerInput> TriggerInputList;
	typedef List<ADCInput> ADCInputList;
	typedef List<ClusterFinder> ClusterFinderList; 
	typedef List<Tracker> TrackerList;
	typedef List<DecisionMaker> DecisionMakerList;
	typedef List<TrackOutput> TrackOutputList;
	typedef List<DecisionOutput> DecisionOutputList;
	
	TriggerInputList triggerinputlist;     // List of trigger record input sources.
	ADCInputList adcinputlist;             // List of ADC stream input sources.
	ClusterFinderList clusterfinderlist;   // List of cluster finders.
	
	TrackerList::Iterator currenttracker;  // Pointer to the current tracker to receive triggers.
	TrackerList trackerlist;               // List of traker objects.
	
	DecisionMakerList decisionmakerlist;   // List of decision making objects.
	
	TrackOutputList trackoutputlist;       // List of track block output sinks.
	DecisionOutputList decisionoutputlist; // List of decision block output sinks.

	ClusterLookupTable clustertable;
};


} // Framework
} // dHLT

#endif // dHLT_FRAMEWORK_GLOBAL_HPP
