////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/MicrodHLT.hpp"

#include "TMethodCall.h"

#include "Version/Version.hpp"
#include "Clustering/ClusterFinder.hpp"
#include "Clustering/CenterOfGravityFinder.hpp"
#include "Tracking/Tracker.hpp"
#include "Tracking/MansoTracker.hpp"
#include "Decision/DecisionMaker.hpp"
#include "AliRoot/convert.hpp"
#include "AliRoot/ClusterFinderProxy.hpp"
#include "AliRoot/TrackerProxy.hpp"

#include "Utils.hpp"
#include "new.hpp"


namespace dHLT
{

using namespace dHLT::AliRoot;


class MicroFramework : public Tracking::TrackerCallback
{
public:

	virtual void RequestClusters(
			Tracking::Tracker* tracker, 
			Float /*left*/, Float /*right*/, Float /*bottom*/, Float /*top*/,
			ChamberID chamber, const void* tag
		)
	{
		Assert( 0 <= chamber and chamber <= 10 );
		DebugMsg(2, "RequestClusters: tag = " << tag);
		register void* ctag = const_cast<void*>(tag);  // We never modify tag so this is OK.
		if (clusters[chamber] != NULL)
			tracker->ReturnClusters(ctag, clusters[chamber], clustercount[chamber]);
		tracker->EndOfClusters(ctag);
	};

	
  virtual void EndOfClusterRequests(Tracking::Tracker* /*tracker*/)
	{
		DebugMsg(2, "EndOfClusterRequests");
		// We can ignore this. Nothing special to do here.
	};

	
	virtual void FoundTrack(Tracking::Tracker* tracker)
	{
		DebugMsg(2, "FoundTrack");
		
		// Fetch the track data from the tracker.
		dHLT::Track newtrack;
		tracker->FillTrackData(newtrack);
		
		// Don't forget to fill the trigger ID number, this is not done by the tracker.
		newtrack.triggerid = currenttriggernumber;

		*trackoutput->AddTrack() = Convert(newtrack);
	};
	
	
	virtual void NoTrackFound(Tracking::Tracker* tracker)
	{
		DebugMsg(2, "NoTrackFound");
		// Again nothing special to do here. The allocated memory is released
		// in the Run method.
	};


	void Run(
			const AliMUONHLT::TriggerSource* triggersource,
			const AliMUONHLT::ClusterSource* clustersource,
			AliMUONHLT::TrackSink* tracksink,
			Tracking::Tracker* tracker
		)
	{
		if (not triggersource->GetFirstEvent()) return;
		while (triggersource->MoreEvents())
		{
			Run(triggersource, clustersource, tracksink, tracker, triggersource->CurrentEvent());
			triggersource->GetNextEvent();
		};
	};


	void Run(
			const AliMUONHLT::TriggerSource* triggersource,
			const AliMUONHLT::ClusterSource* clustersource,
			AliMUONHLT::TrackSink* tracksink,
			Tracking::Tracker* tracker, const Int eventnumber
		)
	{
		// Tell the tracker to make callbacks to this framework object.
		tracker->SetCallback(this);
		trackoutput = tracksink;
		trackoutput->SetNames(triggersource);  // Want the file and folder names to correspond.
		
		CreateClusterBlocks(clustersource, eventnumber);
		try
		{
			trackoutput->AddEvent(eventnumber);
			trackoutput->AddBlock();
			ProcessTriggers(triggersource, tracker);
		}
		finally
		(
			FreeClusterBlocks();
		);
	};


private:

	void CountClusterPoints(const AliMUONHLT::ClusterSource* cs)
	{
		for (Int i = 0; i < 10; i++)
			clustercount[i] = 0;

		cs->GetFirstBlock();
		while (cs->MoreBlocks())
		{
			Int chamber = cs->Chamber();
			if (0 <= chamber and chamber < 10)
			{
				clustercount[chamber] += cs->NumberOfClusters();
			};
			cs->GetNextBlock();
		};
	};
	

	void CreateClusterBlocks(const AliMUONHLT::ClusterSource* cs, Int eventnumber)
	{
		// Must select the proper event before counting or filling the arrays.
		if ( not cs->GetEvent(eventnumber) ) return;
		
		CountClusterPoints(cs);
		
		UInt currentcount[10];
		for (Int i = 0; i < 10; i++)
		{
			// Initialise currentcount.
			currentcount[i] = 0;

			// Allocate arrays.
			if (clustercount[i] > 0)
				clusters[i] = new ClusterPoint[ clustercount[i] ];
			else
				clusters[i] = NULL;
		};

		// Copy all the cluster data into arrays.
		cs->GetFirstBlock();
		while (cs->MoreBlocks())
		{
			Int chamber = cs->Chamber();
			if (0 <= chamber and chamber < 10)
			{
				cs->GetFirstCluster();
				while (cs->MoreClusters())
				{
					ClusterPoint newpoint;
					cs->FetchCluster(newpoint.x, newpoint.y);
					clusters[chamber][currentcount[chamber]++] = newpoint;
					cs->GetNextCluster();
				};
			};
			cs->GetNextBlock();
		};
	};


	void FreeClusterBlocks()
	{
		for (Int i = 0; i < 10; i++)
		{
			if (clusters[i] != NULL)
				delete [] clusters[i];
		};
	};
	
	
	void ProcessTriggers(const AliMUONHLT::TriggerSource* ts, Tracking::Tracker* tracker)
	{
		// The proper event must be selected before calling this method.
	
		ts->GetFirstBlock();
		while (ts->MoreBlocks())
		{
			ts->GetFirstTrigger();
			while (ts->MoreTriggers())
			{
				const AliMUONHLT::TriggerRecord* trigdata = ts->GetTrigger();
				Assert( trigdata != NULL );
				currenttriggernumber = (UInt)trigdata->TriggerNumber();

				TriggerRecord trigger = Convert(*trigdata);

				DebugMsg(2, "Finding track:");
				tracker->FindTrack(trigger);

				DebugMsg(2, "Reset tracker.");
				tracker->Reset();
				
				ts->GetNextTrigger();
			};
			ts->GetNextBlock();
		};
	};
	
	
	UInt clustercount[10];
	ClusterPoint* clusters[10];
	
	AliMUONHLT::TrackSink* trackoutput;  // The current track output object.
	UInt currenttriggernumber;   // The current trigger, trigger number to use.
};


} // dHLT


////////////////////////////////////////////////////////////////////////////////


ClassImp(AliMUONHLT::MicrodHLT)

namespace AliMUONHLT
{


TString Version()
{
	TString str = dHLT::VersionString();
	return str;
};

UInt_t MajorVersion()
{
	return dHLT::MajorVersion();
};

UInt_t MinorVersion()
{
	return dHLT::MinorVersion();
};

UInt_t BuildNumber()
{
	return dHLT::BuildNumber();
};


MicrodHLT::MicrodHLT() : TObject()
{
	fTriggerSource = NULL;
	fClusterSource = NULL;
	fTrackSink = NULL;
	fClusterFinder = NULL;
	fTracker = NULL;
};


MicrodHLT::~MicrodHLT()
{
	// Nothing to do here.
};


void MicrodHLT::SetTriggerSource(const TriggerSource* source)
{
	fTriggerSource = const_cast<TriggerSource*>( source );
};

void MicrodHLT::SetClusterSource(const ClusterSource* source)
{
	fClusterSource = const_cast<ClusterSource*>( source );
};


void MicrodHLT::Run()
{
	DebugMsg(1, "Run");

	//dHLT::Clustering::ClusterFinder* clusterfinder;
	//dHLT::Decision::DecisionMaker* decisionmaker;
	//dHLT::Tracking::Tracker* tracker;
	
	if (fTriggerSource == NULL)
	{
		Error("Run", "The trigger source was not set.");
		return;
	};
	if (fClusterSource == NULL)
	{
		Error("Run", "The cluster source was not set.");
		return;
	};
	if (fTrackSink == NULL)
	{
		Error("Run", "The track output sink was not set.");
		return;
	};
	
	if (fTriggerSource->FileName() != fClusterSource->FileName() or
	    fTriggerSource->FolderName() != fClusterSource->FolderName())
	{
		Warning("Run", "The file and folder names of the trigger source and cluster source do not correspond.");
	};

	dHLT::Clustering::ClusterFinder* clusterfinder = NULL;
	dHLT::Tracking::Tracker* tracker = NULL;
	try
	{
		// Assign the dHLT cluster finder object. If the fClusterFinder field is
		// not set then use the default CenterOfGravityFinder.
		if (fClusterFinder == NULL)
		{
			clusterfinder = new dHLT::Clustering::CenterOfGravityFinder();
		}
		else
		{
			dHLT::AliRoot::ClusterFinderProxy* clusterfinderproxy
				= new dHLT::AliRoot::ClusterFinderProxy(fClusterFinder);
			fClusterFinder->SetCallback(clusterfinderproxy);
			clusterfinder = clusterfinderproxy;
		};

		// Assign the dHLT tracker object. If the fTracker field was not set just
		// use the default MansoTracker implementation.
		if (fTracker == NULL)
		{
			tracker = new dHLT::Tracking::MansoTracker();
		}
		else
		{
			dHLT::AliRoot::TrackerProxy* trackerproxy = new dHLT::AliRoot::TrackerProxy(fTracker);
			fTracker->SetCallback(trackerproxy);
			tracker = trackerproxy;
		};

		dHLT::MicroFramework framework;
		framework.Run(fTriggerSource, fClusterSource, fTrackSink, tracker);
	}
	finally
	(
		if (tracker != NULL) delete tracker;
		if (clusterfinder != NULL) delete clusterfinder;
	);
};


void MicrodHLT::DebugLevel(Int_t value)
{
	DebugCode( dHLT::DebugLevel = value );
};


Int_t MicrodHLT::DebugLevel()
{
#ifdef DEBUG
	return dHLT::DebugLevel;
#else // DEBUG
	return -1;
#endif // DEBUG
};


} // AliMUONHLT
