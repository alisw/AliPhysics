////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DECISION_IO_INTERFACE_HPP
#define dHLT_DECISION_IO_INTERFACE_HPP

#include "BasicTypes.hpp"
#include "EventID.hpp"
#include "DecisionRecord.hpp"
#include "Track.hpp"

namespace dHLT
{
namespace Decision
{


class IOCallback;


class IOInterface
{
public:

	// Do not free memory of tracks until the corresponding ReleaseTracks
	// method is called.
	
	/* When found tracks are made available to the framework then it should
	   pass them on to a decision object using this method.
	 */
	virtual void AddTracks(const EventID event, const Track* tracks, const UInt count) = 0;

	/* When no more tracks can be found for the given event then this method
	   should be called.
	 */
	virtual void EndOfTracks(const EventID event) = 0;

};


class IOCallback
{
public:

	// Do not free the returned memory block until the ReturnDecision method
	// is called for the memory block.
	virtual DecisionRecord* AllocateDecisionBlock(const UInt size) = 0;

	virtual void ReturnDecision(const EventID event, DecisionRecord* decision) = 0;

	virtual void EndOfDecisions(const EventID event) = 0;

	// Do not free track blocks until this method is called for the given block of memory.
	virtual void ReleaseTracks(const Track* tracks) = 0;

};


} // Decision
} // dHLT

#endif // dHLT_DECISION_IO_INTERFACE_HPP
