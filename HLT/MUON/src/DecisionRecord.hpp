////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DECISION_RECORD_HPP
#define dHLT_DECISION_RECORD_HPP

#include "BasicTypes.hpp"
#include "Track.hpp"

namespace dHLT
{


enum DecisionType
{
	UnknownDecision = 0,
	EventAccepted,
	EventRejected
};


struct DecisionRecord
{
	DecisionType type;
};


} // dHLT

#endif // dHLT_DECISION_RECORD_HPP
