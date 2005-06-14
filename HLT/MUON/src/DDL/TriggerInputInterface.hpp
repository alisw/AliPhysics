////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_TRIGGER_INPUT_INTERFACE_HPP
#define dHLT_DDL_TRIGGER_INPUT_INTERFACE_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "EventID.hpp"
#include "TriggerRecord.hpp"

namespace dHLT
{
namespace DDL
{


class TriggerInputCallback
{
public:

	virtual TriggerRecord* AllocateTriggerBlock(const UInt size) = 0;

	virtual void ReturnTriggers(const EventID event, TriggerRecord* triggers, const UInt count) = 0;

	virtual void EndOfTriggers(const EventID event) = 0;

};


} // DDL
} // dHLT

#endif // dHLT_DDL_CLUSTER_INPUT_INTERFACE_HPP
