////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* AliHLTMUONMicrodHLT is a minimalist framework for the dHLT tracking algorithm
   to run in. It integrates all the internal components to execute the algorithm
   properly.
   To run the dHLT with the default Manso algorithm one performs the following
   steps (Assuming we have an initialise AliMUONDataInterface object called
   data):

     // Create a trigger source for the track seeds and populate it with data.
     AliHLTMUONTriggerSource ts(data);
     // We also need the cluster points.
     AliHLTMUONClusterSource cs(data);
     // Need a track sink to store the output data.
     AliHLTMUONTrackSink output;

     // Now we need the framework created and hook up the input and output.
     AliHLTMUONMicrodHLT dhlt;
     dhlt.SetTriggerSource(&ts);
     dhlt.SetClusterSource(&ts);
     dhlt.SetTrackSink(&output);

     // And finally run the algorithm.
     dhlt.Run();
 */

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

#include "AliHLTMUONUtils.h"
#include "AliHLTMUONOutOfMemory.h"


namespace  // AliMicroFramework should be hidden.
{


class AliMicroFramework : public AliHLTMUONCoreTrackerCallback
{
public:

	AliMicroFramework();
	virtual ~AliMicroFramework() {}

	virtual void RequestClusters(
			AliHLTMUONCoreTracker* tracker,
			Float /*left*/, Float /*right*/, Float /*bottom*/, Float /*top*/,
			AliHLTMUONCoreChamberID chamber, const void* tag
		);

	virtual void EndOfClusterRequests(AliHLTMUONCoreTracker* /*tracker*/)
	{
		DebugMsg(2, "EndOfClusterRequests");
		// We can ignore this. Nothing special to do here.
	}

	
	virtual void FoundTrack(AliHLTMUONCoreTracker* tracker);
	
	virtual void NoTrackFound(AliHLTMUONCoreTracker* /*tracker*/)
	{
		DebugMsg(2, "NoTrackFound");
		// Again nothing special to do here. The allocated memory is released
		// in the Run method.
	}


	void Run(
			const AliHLTMUONTriggerSource* triggersource,
			const AliHLTMUONClusterSource* clustersource,
			AliHLTMUONTrackSink* tracksink,
			AliHLTMUONCoreTracker* tracker
		);

	void Run(
			const AliHLTMUONTriggerSource* triggersource,
			const AliHLTMUONClusterSource* clustersource,
			AliHLTMUONTrackSink* tracksink,
			AliHLTMUONCoreTracker* tracker, const Int eventnumber
		);

private:

	// Do not allow copying
	AliMicroFramework(const AliMicroFramework& /*object*/)
		: AliHLTMUONCoreTrackerCallback(), fTrackOutput(NULL), fCurrentTriggerNumber(0)
	{}
	
	AliMicroFramework& operator = (const AliMicroFramework& /*object*/) { return *this; }


	void CountClusterPoints(const AliHLTMUONClusterSource* cs);
	void CreateClusterBlocks(const AliHLTMUONClusterSource* cs, Int eventnumber);
	void FreeClusterBlocks();
	void ProcessTriggers(const AliHLTMUONTriggerSource* ts, AliHLTMUONCoreTracker* tracker);
	
	UInt fClusterCount[10];     // Number of clusters in fClusters.
	AliHLTMUONCoreClusterPoint* fClusters[10];  // Buffers for cluster points received.
	AliHLTMUONTrackSink* fTrackOutput;  // The current track output object.
	UInt fCurrentTriggerNumber;   // The current trigger, trigger number to use.
};

//-----------------------------------------------------------------------------

AliMicroFramework::AliMicroFramework() :
	AliHLTMUONCoreTrackerCallback(), fTrackOutput(NULL), fCurrentTriggerNumber(0)
{
// Default constructor.

	for (Int_t i = 0; i < 10; i++)
	{
		fClusterCount[i] = 0;
		fClusters[i] = NULL;
	}
	fTrackOutput = NULL;
	fCurrentTriggerNumber = 0;
}


void AliMicroFramework::RequestClusters(
		AliHLTMUONCoreTracker* tracker,
		Float /*left*/, Float /*right*/, Float /*bottom*/, Float /*top*/,
		AliHLTMUONCoreChamberID chamber, const void* tag
	)
{
// We return the clusters requested by the tracker.

	Assert( 0 <= chamber && chamber <= 10 );
	DebugMsg(2, "RequestClusters: tag = " << tag);
	register void* ctag = const_cast<void*>(tag);  // We never modify tag so this is OK.
	if (fClusters[chamber] != NULL)
		tracker->ReturnClusters(ctag, fClusters[chamber], fClusterCount[chamber]);
	tracker->EndOfClusters(ctag);
}


void AliMicroFramework::FoundTrack(AliHLTMUONCoreTracker* tracker)
{
// Callback method for the tracker.
// We need to add the found track to the track sink.

	DebugMsg(2, "FoundTrack");
	
	// Fetch the track data from the tracker.
	AliHLTMUONCoreTrack newtrack;
	tracker->FillTrackData(newtrack);
	
	// Don't forget to fill the trigger ID number, this is not done by the tracker.
	newtrack.fTriggerid = fCurrentTriggerNumber;

	*fTrackOutput->AddTrack() = AliHLTMUONConvert(newtrack);
}
	

void AliMicroFramework::Run(
		const AliHLTMUONTriggerSource* triggersource,
		const AliHLTMUONClusterSource* clustersource,
		AliHLTMUONTrackSink* tracksink,
		AliHLTMUONCoreTracker* tracker
	)
{
// Run the micro dHLT chain.

	if ( ! triggersource->GetFirstEvent()) return;
	while (triggersource->MoreEvents())
	{
		Run(triggersource, clustersource, tracksink, tracker, triggersource->CurrentEvent());
		triggersource->GetNextEvent();
	}
}


void AliMicroFramework::Run(
		const AliHLTMUONTriggerSource* triggersource,
		const AliHLTMUONClusterSource* clustersource,
		AliHLTMUONTrackSink* tracksink,
		AliHLTMUONCoreTracker* tracker, const Int eventnumber
	)
{
// Process a particlar event from the trigger and cluster source.

	// Tell the tracker to make callbacks to this framework object.
	tracker->SetCallback(this);
	fTrackOutput = tracksink;
	fTrackOutput->SetNames(triggersource);  // Want the file and folder names to correspond.
	
	CreateClusterBlocks(clustersource, eventnumber);
	try
	{
		fTrackOutput->AddEvent(eventnumber);
		fTrackOutput->AddBlock();
		ProcessTriggers(triggersource, tracker);
	}
	finally
	(
		FreeClusterBlocks();
	);
}


void AliMicroFramework::CountClusterPoints(const AliHLTMUONClusterSource* cs)
{
// Count the number of clusters on each chamber in the cluster source.

	for (Int i = 0; i < 10; i++)
		fClusterCount[i] = 0;

	cs->GetFirstBlock();
	while (cs->MoreBlocks())
	{
		Int chamber = cs->Chamber();
		if (0 <= chamber && chamber < 10)
		{
			fClusterCount[chamber] += cs->NumberOfClusters();
		}
		cs->GetNextBlock();
	}
}
	

void AliMicroFramework::CreateClusterBlocks(
		const AliHLTMUONClusterSource* cs, Int eventnumber
	)
{
// Fill the fClusters buffers from cluster source.

	// Must select the proper event before counting or filling the arrays.
	if ( ! cs->GetEvent(eventnumber) ) return;
	
	CountClusterPoints(cs);
	
	UInt currentcount[10];
	for (Int i = 0; i < 10; i++)
	{
		// Initialise currentcount.
		currentcount[i] = 0;

		// Allocate arrays.
		if (fClusterCount[i] > 0)
			fClusters[i] = new AliHLTMUONCoreClusterPoint[ fClusterCount[i] ];
		else
			fClusters[i] = NULL;
	}

	// Copy all the cluster data into arrays.
	cs->GetFirstBlock();
	while (cs->MoreBlocks())
	{
		Int chamber = cs->Chamber();
		if (0 <= chamber && chamber < 10)
		{
			cs->GetFirstCluster();
			while (cs->MoreClusters())
			{
				AliHLTMUONCoreClusterPoint newpoint;
				cs->FetchCluster(newpoint.X(), newpoint.Y());
				fClusters[chamber][currentcount[chamber]++] = newpoint;
				cs->GetNextCluster();
			}
		}
		cs->GetNextBlock();
	}
}


void AliMicroFramework::FreeClusterBlocks()
{
// Release the fClusters buffers.

	for (Int i = 0; i < 10; i++)
	{
		if (fClusters[i] != NULL)
		{
			delete [] fClusters[i];
			fClusters[i] = NULL;
		}
	}
}
	
	
void AliMicroFramework::ProcessTriggers(
		const AliHLTMUONTriggerSource* ts, AliHLTMUONCoreTracker* tracker
	)
{
// Go through the list of triggers and have tracker try find the track using
// the trigger as a seed.
// Note: The proper event must be selected before calling this method.
	
	ts->GetFirstBlock();
	while (ts->MoreBlocks())
	{
		ts->GetFirstTrigger();
		while (ts->MoreTriggers())
		{
			const AliHLTMUONTriggerRecord* trigdata = ts->GetTrigger();
			Assert( trigdata != NULL );
			fCurrentTriggerNumber = (UInt)trigdata->TriggerNumber();

			AliHLTMUONCoreTriggerRecord trigger = AliHLTMUONConvert(*trigdata);

			DebugMsg(2, "Finding track:");
			tracker->FindTrack(trigger);

			DebugMsg(2, "Reset tracker.");
			tracker->Reset();
			
			ts->GetNextTrigger();
		}
		ts->GetNextBlock();
	}
}


} // end of namespace


////////////////////////////////////////////////////////////////////////////////


ClassImp(AliHLTMUONMicrodHLT)


TString AliHLTMUONVersion()
{
	TString str = dHLT::VersionString();
	return str;
}

UInt_t AliHLTMUONMajorVersion()
{
	return dHLT::MajorVersion();
}

UInt_t AliHLTMUONMinorVersion()
{
	return dHLT::MinorVersion();
}

UInt_t AliHLTMUONBuildNumber()
{
	return dHLT::BuildNumber();
}


AliHLTMUONMicrodHLT::AliHLTMUONMicrodHLT() :
	TObject(),
	fTriggerSource(NULL), fClusterSource(NULL), fTrackSink(NULL),
	fClusterFinder(NULL), fTracker(NULL)
{
// Default constructor

	fTriggerSource = NULL;
	fClusterSource = NULL;
	fTrackSink = NULL;
	fClusterFinder = NULL;
	fTracker = NULL;
}


AliHLTMUONMicrodHLT::~AliHLTMUONMicrodHLT()
{
	// Nothing to do here.
}


void AliHLTMUONMicrodHLT::SetTriggerSource(const AliHLTMUONTriggerSource* source)
{
	fTriggerSource = const_cast<AliHLTMUONTriggerSource*>( source );
}

void AliHLTMUONMicrodHLT::SetClusterSource(const AliHLTMUONClusterSource* source)
{
	fClusterSource = const_cast<AliHLTMUONClusterSource*>( source );
}


void AliHLTMUONMicrodHLT::Run()
{
// The entry point routine for the dHLT algorithm in the micro format.
// To run the dHLT set the input and output objects with the set methods
// provided and then call this Run method.

	DebugMsg(1, "Run");

	//AliHLTMUONCoreClusterFinder* clusterfinder;
	//dHLT::Decision::DecisionMaker* decisionmaker;
	//AliHLTMUONCoreTracker* tracker;
	
	if (fTriggerSource == NULL)
	{
		Error("Run", "The trigger source was not set.");
		return;
	}
	if (fClusterSource == NULL)
	{
		Error("Run", "The cluster source was not set.");
		return;
	}
	if (fTrackSink == NULL)
	{
		Error("Run", "The track output sink was not set.");
		return;
	}
	
	if (fTriggerSource->FileName() != fClusterSource->FileName() ||
	    fTriggerSource->FolderName() != fClusterSource->FolderName())
	{
		Warning("Run", "The file and folder names of the trigger source and cluster source do not correspond.");
	}

	AliHLTMUONCoreClusterFinder* clusterfinder = NULL;
	AliHLTMUONCoreTracker* tracker = NULL;
	try
	{
		// Assign the dHLT cluster finder object. If the fClusterFinder field is
		// not set then use the default CenterOfGravityFinder.
		if (fClusterFinder == NULL)
		{
			clusterfinder = new AliHLTMUONCoreCenterOfGravityFinder();
		}
		else
		{
			AliHLTMUONClusterFinderProxy* clusterfinderproxy
				= new AliHLTMUONClusterFinderProxy(fClusterFinder);
			fClusterFinder->SetCallback(clusterfinderproxy);
			clusterfinder = clusterfinderproxy;
		}

		// Assign the dHLT tracker object. If the fTracker field was not set just
		// use the default MansoTracker implementation.
		if (fTracker == NULL)
		{
			tracker = new AliHLTMUONCoreMansoTracker();
		}
		else
		{
			AliHLTMUONTrackerProxy* trackerproxy = new AliHLTMUONTrackerProxy(fTracker);
			fTracker->SetCallback(trackerproxy);
			tracker = trackerproxy;
		}

		AliMicroFramework framework;
		framework.Run(fTriggerSource, fClusterSource, fTrackSink, tracker);
	}
	finally
	(
		if (tracker != NULL) delete tracker;
		if (clusterfinder != NULL) delete clusterfinder;
	);
}


#ifdef DEBUG
void AliHLTMUONMicrodHLT::DebugLevel(Int_t value)
{
	gAliHLTMUONDebugLevel = value;
}
#else // DEBUG
void AliHLTMUONMicrodHLT::DebugLevel(Int_t /*value*/) {}
#endif // DEBUG


Int_t AliHLTMUONMicrodHLT::DebugLevel()
{
// Returns the debug level for the dHLT module.

#ifdef DEBUG
	return gAliHLTMUONDebugLevel;
#else // DEBUG
	return -1;
#endif // DEBUG
}
