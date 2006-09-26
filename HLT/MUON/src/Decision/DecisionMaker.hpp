////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DECISION_DECISION_MAKER_HPP
#define dHLT_DECISION_DECISION_MAKER_HPP

#include "AliHLTMUONBasicTypes.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONCoreTrack.h"
#include "DecisionRecord.hpp"


namespace dHLT
{
namespace Decision
{


class DecisionMaker;


class DecisionMakerCallback
{
public:

	virtual ~DecisionMakerCallback() {};

	/* Called when the decision maker has finished processing.
	   At this point and tracks passed to the decision maker are no longer is use
	   by the decision maker and the tracks can be released.
	   'decisionsize' indicates the number of bytes required to store the decision.
	   If decisionsize == 0 then no memory is allocated and FillDecisionData should
	   never be called.
	 */
	virtual void MadeDecision(DecisionMaker* decisionmaker, const UInt decisionsize) = 0;

};


class DecisionMaker
{
public:

	DecisionMaker() : fCallback(NULL) {};
	

	DecisionMaker(const DecisionMaker& dm) : fCallback(dm.fCallback) {};

	DecisionMaker& operator = (const DecisionMaker& dm)
	{
		fCallback = dm.fCallback;
		return *this;
	}


	virtual ~DecisionMaker() {};

	/* This is the starting point of the decision algorithm. Any initialisation
	   processing should be performed here.
	 */
	virtual void StartDecision() = 0;
	
	/* For every track avaiable for a particular event or sub-event this method
	   is called with the relevant track. Most of the actual decision algorithm 
	   will go into this method on a per track basis.
	 */
	virtual void AddTrack(const AliHLTMUONCoreTrack& track) = 0;
	
	/* When no more tracks are available for the decision then this method is called.
	   Final processing of the decision algorithm should go in to this method.
	   When the decision is made by the algorithm then the MadeDecision should be
	   called.
	 */
	virtual void EndOfTracks() = 0;
	
	/* After a call to MadeDecision this method will be called to retreive the
	   decision. The 'decision' pointer should point to a memory block of at least
	   the number of bytes specified in the call to MadeDecision, with the 
	   decisionsize parameter. The decision algorithm should not write more data
	   than requested.
	 */
	virtual void FillDecisionData(DecisionRecord* decision) = 0;

	/* This is called when the decision maker should be reset to an initial state.
	   All extra internal memory allocated during processing should be released.
	 */
	virtual void Reset() = 0;

	/* Sets the DecisionMakerCallback callback interface.
	 */
	inline void SetCallback(DecisionMakerCallback* callback)
	{
		this->fCallback = callback;
	};

protected:

	/* When the decision algorithm is finished processing and a decision is made,
	   then this method should be called. At this point all the track blocks
	   should be considered released and the track memory MUST NOT be accessed.
	   The number of bytes required to store the decision should be indecated by
	   'decisionsize'.
	 */
	inline void MadeDecision(const UInt decisionsize)
	{
		Assert( fCallback != NULL );
		fCallback->MadeDecision(this, decisionsize);
	};

private:

	DecisionMakerCallback* fCallback;
};


} // Decision
} // dHLT

#endif // dHLT_DECISION_DECISION_MAKER_HPP
