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

#ifndef ALIHLTMUONMICRODHLT_H
#define ALIHLTMUONMICRODHLT_H

#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "AliRoot/TriggerSource.hpp"
#include "AliRoot/ClusterSource.hpp"
#include "AliRoot/TrackSink.hpp"
#include "AliRoot/TrackerCallback.hpp"


class AliHLTMUONDummyClusterFinder;
class AliHLTMUONClusterFinderInterface;
class AliHLTMUONDummyTracker;
class AliHLTMUONTrackerInterface;


/* Routines for getting dHLT version information.
*/
TString AliHLTMUONVersion();
UInt_t AliHLTMUONMajorVersion();
UInt_t AliHLTMUONMinorVersion();
UInt_t AliHLTMUONBuildNumber();


class AliHLTMUONMicrodHLT : public TObject
{
public:

	AliHLTMUONMicrodHLT();
	virtual ~AliHLTMUONMicrodHLT();

	/* Get/Set methods for the trigger data source.
	   Note: The source object must be cleaned up by the caller.
	 */
	void SetTriggerSource(const AliHLTMUONTriggerSource* source);
	const AliHLTMUONTriggerSource* GetTriggerSource() const { return fTriggerSource; };
	
	/* Get/Set methods for the cluster data source.
	   Note: The source object must be cleaned up by the caller.
	 */
	void SetClusterSource(const AliHLTMUONClusterSource* source);
	const AliHLTMUONClusterSource* GetClusterSource() const { return fClusterSource; };
	
	/* Get/Set methods for the track data sink (output target).
	   Note: The output object must be cleaned up by the caller.
	 */
	void SetTrackSink(AliHLTMUONTrackSink* sink)    { fTrackSink = sink; };
	const AliHLTMUONTrackSink* GetTrackSink() const { return fTrackSink; };
	
	/* Get/Set methods for the cluster finder interface.
	   Note: The cluster finder object must be cleaned up by the caller.
	   This field is optional. Use it if a custom cluster finder should be used.
	 */
	void SetClusterFinder(AliHLTMUONClusterFinderInterface* clusterfinder)
	{
		fClusterFinder = clusterfinder;
	}

	const AliHLTMUONClusterFinderInterface* GetClusterFinder() const { return fClusterFinder; };
	void SetClusterFinder(AliHLTMUONDummyClusterFinder* clusterfinder);
	
	/* Get/Set methods for the tracker interface.
	   Note: The tracker object must be cleaned up by the caller.
	   This field is optional. Use it if a custom tracker should be used.
	 */
	void SetTracker(AliHLTMUONTrackerInterface* tracker) { fTracker = tracker; };
	const AliHLTMUONTrackerInterface* GetTracker() const { return fTracker; };
	void SetTracker(AliHLTMUONDummyTracker* tracker);
	
	/* The entry point routine for the dHLT algorithm in the micro format.
	   To run the dHLT set the input and output objects with the set methods
	   provided and then call this Run method.
	 */
	void Run();
	
	// Get and set methods for the dHLT debug level.
	static void DebugLevel(Int_t value);
	static Int_t DebugLevel();
	
private:
	
	// Do not allow copying of this object.
	AliHLTMUONMicrodHLT(const AliHLTMUONMicrodHLT& /*object*/)
		: TObject(), fTriggerSource(NULL), fClusterSource(NULL),
		  fTrackSink(NULL), fClusterFinder(NULL), fTracker(NULL)
	{}

	AliHLTMUONMicrodHLT& operator = (const AliHLTMUONMicrodHLT& /*object*/) { return *this; }

	AliHLTMUONTriggerSource* fTriggerSource;           //! Trigger record input source.
	AliHLTMUONClusterSource* fClusterSource;           //! Cluster point input source.
	AliHLTMUONTrackSink* fTrackSink;                   //! Track output sink.
	AliHLTMUONClusterFinderInterface* fClusterFinder;  //! Interface to a custom cluster finder.
	AliHLTMUONTrackerInterface* fTracker;              //! Interface to a custom tracker.

	ClassDef(AliHLTMUONMicrodHLT, 0);  // A very minimal implementation of the dHLT algorithm.
};


#endif // ALIHLTMUONMICRODHLT_H
