////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_TRACKING_IO_INTERFACE_HPP
#define dHLT_TRACKING_IO_INTERFACE_HPP

#include "BasicTypes.hpp"
#include "EventID.hpp"
#include "TriggerRecord.hpp"
#include "RegionOfInterest.hpp"
#include "Cluster.hpp"
#include "Track.hpp"

namespace dHLT
{
namespace Tracking
{


class IOInterface
{
public:

	// Do not free memory of triggers nor clusters for a given event until the
	// corresponding ReleaseTriggers or ReleaseClusters method is called.

	/* Called by the framework when a new trigger record is to be processed.
	 */
	virtual void AddTriggers(const EventID event, const TriggerRecord* triggers, const UInt count) = 0;

	/* When no more triggers are expected by the framework then this method is called.
	 */
	virtual void EndOfTriggers(const EventID event) = 0;

	/* Called by the framework when the requested cluster blocks are to be processed.
	 */
	virtual void AddClusters(const EventID event, const ROI region, const ClusterPoint* clusters, const UInt count) = 0;

	/* When no more clusters are expected by the framework then this method is called.
	 */
	virtual void EndOfClusters(const EventID event, const ROI region) = 0;

};


class IOCallback
{
public:

	/* Called when clusters from the specified region of interest are required
	   to reconstruct the track.
	 */
	virtual void RequestClusters(const EventID event, const ROI region) = 0;

	/* When no more clusters are required for a given event then this method
	   is called.
	 */
	virtual void EndOfClusterRequests(const EventID event) = 0;

	// Do not free the returned memory until ReturnTracks is called
	// with the same memory pointer.
	virtual Track* AllocateTrackBlock(const UInt size) = 0;

	/* When tracks have been found they are returned to the framework with
	   this method.
	 */
	virtual void ReturnTracks(const EventID event, Track* newtracks, const UInt count) = 0;

	/* When no more tracks can be found for this event then this method is called.
	 */
	virtual void EndOfTracks(const EventID event) = 0;

	// Do not free trigger blocks until this method is called for the given block of memory.
	virtual void ReleaseTriggers(const TriggerRecord* triggers) = 0;

	// Do not free cluster blocks until this method is called for the given block of memory.
	virtual void ReleaseClusters(const ClusterPoint* clusters) = 0;

};


} // Tracking
} // dHLT

#endif // dHLT_TRACKING_IO_INTERFACE_HPP
