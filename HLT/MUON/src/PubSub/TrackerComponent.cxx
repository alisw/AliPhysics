////////////////////////////////////////////////////////////////////////////////
//
// Author: Gareth de Vaux
// Email:  devaux@lhc.phy.uct.ac.za | dhlt@lordcow.org
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "PubSub/TrackerComponent.hpp"
#include "Cluster.hpp"

#define CLUSTERS_DATAID (((AliUInt64_t)'CLUS')<<32 | 'TERS')
#define TRIGGERS_DATAID (((AliUInt64_t)'TRIG')<<32 | 'GERS')

using namespace dHLT;

namespace dHLT
{
namespace PubSub
{


TrackerComponent::TrackerComponent(const char* name, bool sendEventDone, AliUInt32_t minBlockSize, int eventCountAlloc)
	: AliHLTProcessingSubscriber(name, sendEventDone, minBlockSize, eventCountAlloc)
{  
	tracker.SetCallback(&callback);
	triggerblocks.erase(triggerblocks.begin(), triggerblocks.end());
}


bool TrackerComponent::ProcessEvent(AliEventID_t /*eventID*/, AliUInt32_t blockCnt, AliHLTSubEventDescriptor::BlockData* blocks, const AliHLTSubEventDataDescriptor* /*sedd*/, const AliHLTEventTriggerStruct* /*etsp*/, AliUInt8_t* outputPtr, AliUInt32_t& size, vector<AliHLTSubEventDescriptor::BlockData>& outputBlocks, AliHLTEventDoneData*& edd)
{
	edd = NULL; // for now

	// go through the list of input blocks and fill internal pointers to them

	AliUInt32_t i;

	for (i=0; i<blockCnt; i++) 
	{
		switch (blocks[i].fDataType.fID)
		{
		case CLUSTERS_DATAID:

			LOG(AliHLTLog::kDebug, "TrackerComponent::ProcessEvent", "Trace")
				<< "Got CLUSTERS block with specification = " << blocks[i].fDataSpecification
				<< " , of length = " << blocks[i].fSize << ENDLOG;

			if (blocks[i].fDataSpecification < 10)
				callback.AddClusterBlock((ClusterPoint*)blocks[i].fData, blocks[i].fSize, blocks[i].fDataSpecification);
			else
				LOG(AliHLTLog::kError, "TrackerComponent::ProcessEvent", "Unknown Chamber")
					<< "Invalid data specification: " << blocks[i].fDataSpecification
					<< " , expected a chamber value in the range 0..9" << ENDLOG;

			break;

		case TRIGGERS_DATAID:

			LOG(AliHLTLog::kDebug, "TrackerComponent::ProcessEvent", "Trace")
				<< "Got TRIGGERS block with specification = " << blocks[i].fDataSpecification
				<< " , of length = " << blocks[i].fSize << ENDLOG;

			UInt *trig_ptr, triggerIDoffset;

			trig_ptr = (UInt*)blocks[i].fData;			
			triggerIDoffset = (UInt)*trig_ptr++;

			AddTriggerBlock(triggerIDoffset, (TriggerRecord*)trig_ptr, blocks[i].fSize - sizeof(UInt), blocks[i].fDataSpecification);

			break;

		default:

			LOG(AliHLTLog::kError, "TrackerComponent::ProcessEvent", "Unknown Data")
				<< "Received unknown data type: 0x" << AliHLTLog::kHex 
				<< blocks[i].fDataType.fID << " (" << AliHLTLog::kDec
				<< blocks[i].fDataType.fID << ")." << ENDLOG;
		}
	}

	callback.SetTracks((Track*)outputPtr);

	// now we can go through the list of trigger records and for each try find
	// a dHLT track using the MansoTracker

	for (i=0; i<triggerblocks.size(); i++)
	{
		TriggerBlock &tblock = triggerblocks[i];

		for (UInt j=0; j<tblock.count; j++)
		{
			LOG(AliHLTLog::kDebug, "TrackerComponent::ProcessEvent", "Trace")
				<< "tracker.FindTrack called for TriggerRecord# " << j     
				<< " in block# " << i << ENDLOG;   

			tracker.FindTrack(tblock.data[j]);
			tracker.Reset();

			LOG(AliHLTLog::kDebug, "TrackerComponent::ProcessEvent", "Trace")
				<< "End of tracker step, trackcount = " << callback.TrackCount() << ENDLOG;
		}
	}

	AliUInt32_t outputsize = callback.TrackCount() * sizeof(Track);

	if (outputsize <= size)
	{
		AliHLTSubEventDescriptor::BlockData block;          // to be inserted into outputBlocks

		block.fShmKey.fShmType = kInvalidShmType;                 // framework will insert this
		block.fShmKey.fKey.fID =~ (AliUInt64_t)0;
		block.fOffset = 0;                         // this block's at the beginning ouput block
		block.fSize = outputsize;
		block.fProducerNode = (AliHLTNodeID_t)-1;                 // framework will insert this
		block.fDataType.fID = TRACKS_DATAID;
		block.fDataOrigin = blocks[0].fDataOrigin;                    // can set this if needed
		block.fDataSpecification = blocks[0].fDataSpecification;      // can set this if needed
		block.fByteOrder = kAliHLTNativeByteOrder;

		outputBlocks.insert(outputBlocks.end(), block);
		size = outputsize; // since it's only this one block in the output
		return true;
	}

	else
	{
		LOG(AliHLTLog::kError, "TrackerComponent::ProcessEvent", "Insufficient output buffer")
			<< "Insufficient output buffer: " << size << " allocated, "
			<< outputsize << " to be written" << ENDLOG;

		return false;
	}
}

void TrackerComponent::AddTriggerBlock(const UInt triggerIDoffset, const TriggerRecord* triggers, const UInt size, const UInt chamber)
{
	Assert(chamber<10);

	TriggerBlock newrec;
	newrec.triggerIDoffset = triggerIDoffset;
	newrec.data = triggers;
	newrec.count = size / sizeof(TriggerRecord);
	triggerblocks.push_back(newrec);
}


} // PubSub
} // dHLT


using namespace dHLT::PubSub;

int main(int argc, char** argv)
{
	AliHLTControlledDefaultProcessorComponent<TrackerComponent> procComp("Tracker", argc, argv);
	procComp.Run();                                                                                                                   
	return 0;    
}
