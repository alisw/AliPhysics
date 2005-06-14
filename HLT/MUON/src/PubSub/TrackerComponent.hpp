////////////////////////////////////////////////////////////////////////////////
//
// Author: Gareth de Vaux
// Email:  devaux@lhc.phy.uct.ac.za | dhlt@lordcow.org
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_PUBSUB_TRACKERCOMPONENT_HPP
#define dHLT_PUBSUB_TRACKERCOMPONENT_HPP

#include "AliHLTProcessingComponent.hpp"
#include "AliHLTProcessingSubscriber.hpp"
#include "AliHLTControlledComponent.hpp"

#include "PubSub/TrackerCallback.hpp" // need this?
#include "TriggerRecord.hpp"

using namespace dHLT;

namespace dHLT
{
namespace PubSub
{


class TrackerComponent: public AliHLTProcessingSubscriber
{                                                                  
public:

	TrackerComponent(const char* name, bool sendEventDone, AliUInt32_t minBlockSize, int eventCountAlloc);

protected:

	bool ProcessEvent(AliEventID_t eventID, AliUInt32_t blockCnt, AliHLTSubEventDescriptor::BlockData* blocks, const AliHLTSubEventDataDescriptor* sedd, const AliHLTEventTriggerStruct* etsp, AliUInt8_t* outputPtr, AliUInt32_t& size, vector<AliHLTSubEventDescriptor::BlockData>& outputBlocks, AliHLTEventDoneData*& edd);

private:

	void AddTriggerBlock(const UInt triggerIDoffset, const TriggerRecord* triggers, const UInt size, const UInt chamber);

	TrackerCallback callback;
	dHLT::Tracking::MansoTracker tracker;
	
	struct TriggerBlock
	{
		UInt triggerIDoffset;
		const TriggerRecord *data;
		UInt count;
	};

	std::vector<TriggerBlock> triggerblocks;  // array of trigger record blocks
};


} // PubSub
} // dHLT

#endif // dHLT_PUBSUB_TRACKERCOMPONENT_HPP
