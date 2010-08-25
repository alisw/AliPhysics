//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCDATACHECKERCOMPONENT_H
#define ALIHLTTPCDATACHECKERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// \file   AliHLTTPCDataCheckerComponent.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   9 Aug 2010
/// \brief  Declaration of the AliHLTTPCDataCheckerComponent component class.

#include "AliHLTProcessor.h"

class AliTPCRawStreamV3;
class AliRawReaderMemory;

/**
 * \class AliHLTTPCDataCheckerComponent
 * The TPC data checker component is used to validate the raw data entering the HLT
 * from the TPC detector. Basic sanity and data integrity checks are performed.
 * Any problems are logged with explanatory error messages. By default all data blocks
 * are forwarded, but one can also optionally filter on the bad data blocks.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCDataChecker <br>
 * Library: \b libAliHLTTPC.so   <br>
 * Input Data Types:  ::kAliHLTAnyDataType | ::kAliHLTDataOriginTPC <br>
 * Output Data Types: ::kAliHLTAnyDataType | ::kAliHLTDataOriginTPC <br>
 *
 * <h2>Mandatory arguments:</h2>
 * None.
 *
 * <h2>Optional arguments:</h2>
 * \li -filter <i>flag</i> <br>
 *      If specified then all data blocks for which a problem was found are forwarded
 *      and all others are dropped. If the optional <i>flag</i> is specified then it
 *      can be one of:
 *        forwardbad  - only the bad data blocks are forwarded (default option).
 *        forwardgood - only the good data blocks are forwarded.
 * \li -ignoretype <br>
 *      If set then the check of the data type is not performed.
 * \li -ignoreorigin <br>
 *      If set then the check of the origin is not performed.
 * \li -ignorespec <br>
 *      If set then the check of the block specification is not performed.
 * \li -handle-all-events <br>
 *      If set then all events are handled and not just data events.
 *
 * <h2>Configuration:</h2>
 * Can only be configured with the command line arguments.
 *
 * <h2>Default CDB entries:</h2>
 * None.
 *
 * <h2>Performance:</h2>
 * Can run over 3kHz in HLT online system.
 *
 * <h2>Memory consumption:</h2>
 * Negligible.
 *
 * <h2>Output size:</h2>
 * The same as the input data size.
 *
 * \ingroup alihlt_tpc_components
 */
class AliHLTTPCDataCheckerComponent : public AliHLTProcessor
{
public:
	
	AliHLTTPCDataCheckerComponent();
	virtual ~AliHLTTPCDataCheckerComponent();
	
	// Methods inherited from AliHLTComponent:
	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();
	virtual Int_t DoInit(int argc, const char** argv);
	virtual Int_t DoDeinit();
	
protected:
	
	// Method inherited from AliHLTProcessor:
	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks, 
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr, 
			AliHLTUInt32_t& size,
			AliHLTComponentBlockDataList& outputBlocks
		);
	
	using AliHLTProcessor::DoEvent;
	
private:
	
	// Do not allow copying of this class.
	AliHLTTPCDataCheckerComponent(const AliHLTTPCDataCheckerComponent& obj);
	AliHLTTPCDataCheckerComponent& operator = (const AliHLTTPCDataCheckerComponent& obj);
	
	/**
	 * Checks the structure of a TPC raw DDL data block.
	 * \param event  The event this data block was found in.
	 * \param index  The index number of the data block as found in the
	 *               DoEvent::blocks parameter.
	 * \param block  The data block to check. Must correspond to 'index'.
	 * \returns true if the raw data is OK and false otherwise.
	 */
	bool CheckRawDataBlock(AliHLTEventID_t event, AliHLTUInt32_t index, const AliHLTComponentBlockData* block);
	
	AliTPCRawStreamV3* fRawStream;  /// Raw stream to read TPC data.
	AliRawReaderMemory* fRawReader;  /// Raw reader object for the TPC stream.
	bool fForwardBadBlocks;  /// Flag indicating if the bad blocks should be forwarded.
	bool fForwardGoodBlocks;  /// Flag indicating if the good blocks should be forwarded.
	bool fIgnoreType;     /// Indicates if the data block type should not be checked.
	bool fIgnoreOrigin;   /// Indicates if the data block origin should not be checked.
	bool fIgnoreSpec;     /// Indicates if the data block specification bits should not be checked.
	bool fHandleAllEvents; /// Indicates if all event types are processed and not just data events.
	
	ClassDef(AliHLTTPCDataCheckerComponent, 0)  // Data sanity and integrity checker component for TPC data.
};

#endif // ALIHLTTPCDATACHECKERCOMPONENT_H
