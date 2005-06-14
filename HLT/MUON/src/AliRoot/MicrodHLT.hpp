////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_MICRODHLT_HPP
#define dHLT_ALIROOT_MICRODHLT_HPP

#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "AliRoot/TriggerSource.hpp"
#include "AliRoot/ClusterSource.hpp"
#include "AliRoot/TrackSink.hpp"
#include "AliRoot/TrackerCallback.hpp"

namespace AliMUONHLT
{


class ClusterFinder;
class ClusterFinderInterface;
class Tracker;
class TrackerInterface;


/* Routines for getting dHLT version information.
*/
TString Version();
UInt_t MajorVersion();
UInt_t MinorVersion();
UInt_t BuildNumber();


class MicrodHLT : public TObject
{
public:

	MicrodHLT();
	virtual ~MicrodHLT();

	/* Get/Set methods for the trigger data source.
	   Note: The source object must be cleaned up by the caller.
	 */
	void SetTriggerSource(const TriggerSource* source);
	const TriggerSource* GetTriggerSource() const { return fTriggerSource; };
	
	/* Get/Set methods for the cluster data source.
	   Note: The source object must be cleaned up by the caller.
	 */
	void SetClusterSource(const ClusterSource* source);
	const ClusterSource* GetClusterSource() const { return fClusterSource; };
	
	/* Get/Set methods for the track data sink (output target).
	   Note: The output object must be cleaned up by the caller.
	 */
	void SetTrackSink(TrackSink* sink)    { fTrackSink = sink; };
	const TrackSink* GetTrackSink() const { return fTrackSink; };
	
	/* Get/Set methods for the cluster finder interface.
	   Note: The cluster finder object must be cleaned up by the caller.
	   This field is optional. Use it if a custom cluster finder should be used.
	 */
	void SetClusterFinder(ClusterFinderInterface* clusterfinder) { fClusterFinder = clusterfinder; };
	const ClusterFinderInterface* GetClusterFinder() const { return fClusterFinder; };
	void SetClusterFinder(ClusterFinder* clusterfinder);
	
	/* Get/Set methods for the tracker interface.
	   Note: The tracker object must be cleaned up by the caller.
	   This field is optional. Use it if a custom tracker should be used.
	 */
	void SetTracker(TrackerInterface* tracker) { fTracker = tracker; };
	const TrackerInterface* GetTracker() const { return fTracker; };
	void SetTracker(Tracker* tracker);
	
	/* The entry point routine for the dHLT algorithm in the micro format.
	   To run the dHLT set the input and output objects with the set methods
	   provided and then call this Run method.
	 */
	void Run();
	
	// Get and set methods for the dHLT debug level.
	static void DebugLevel(Int_t value);
	static Int_t DebugLevel();
	
private:

	TriggerSource* fTriggerSource;           //! Trigger record input source.
	ClusterSource* fClusterSource;           //! Cluster point input source.
	TrackSink* fTrackSink;                   //! Track output sink.
	ClusterFinderInterface* fClusterFinder;  //! Interface to a custom cluster finder.
	TrackerInterface* fTracker;              //! Interface to a custom tracker.

	ClassDef(MicrodHLT, 0);  // A very minimal implementation of the dHLT algorithm.
};


} // AliMUONHLT

#endif // dHLT_ALIROOT_MICRODHLT_HPP
