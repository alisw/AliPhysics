////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_DECISION_OUTPUT_INTERFACE_HPP
#define dHLT_DDL_DECISION_OUTPUT_INTERFACE_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "EventID.hpp"
#include "DecisionRecord.hpp"

namespace dHLT
{
namespace DDL
{


class DecisionOutputInterface
{
public:

	virtual void ReturnDecision(const EventID event, const DecisionRecord* decision) = 0;

	virtual void EndOfDecisions(const EventID event) = 0;

};


class DecisionOutputCallback
{
public:
	virtual void ReleaseDecision(const DecisionRecord* decision) = 0;
};


} // DDL
} // dHLT

#endif // dHLT_DDL_DECISION_OUTPUT_INTERFACE_HPP
