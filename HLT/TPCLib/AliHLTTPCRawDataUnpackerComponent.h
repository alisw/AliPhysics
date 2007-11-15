// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCRAWDATAUNPACKERCOMPONENT_H
#define ALIHLTTPCRAWDATAUNPACKERCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

class AliRawReaderMemory;
class AliTPCRawStream;

/**
 * @class AliHLTTPCRawDataUnpackerComponent
 * Unpacker component for TPC RAW data.
 */
class AliHLTTPCRawDataUnpackerComponent : public AliHLTProcessor
    {
    public:
	AliHLTTPCRawDataUnpackerComponent();
	virtual ~AliHLTTPCRawDataUnpackerComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
	AliHLTComponent* Spawn();

    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

	using AliHLTProcessor::DoEvent;
	
    private:
      /** copy constructor prohibited */
      AliHLTTPCRawDataUnpackerComponent(const AliHLTTPCRawDataUnpackerComponent&);
      /** assignment operator prohibited */
      AliHLTTPCRawDataUnpackerComponent& operator=(const AliHLTTPCRawDataUnpackerComponent&);

	// Initialize AliROOT TPC raw stream parsing class
      AliRawReaderMemory *fRawMemoryReader; //! transient
      AliTPCRawStream *fTPCRawStream; //! transient

      ClassDef(AliHLTTPCRawDataUnpackerComponent, 0);

    };
#endif
