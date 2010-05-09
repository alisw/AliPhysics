//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSRCUHISTOGRAMPRODUCERCOMPONENT_H
#define ALIHLTPHOSRCUHISTOGRAMPRODUCERCOMPONENT_H 

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTCaloProcessor.h"
#include "AliHLTPHOSUtilities.h"

//
// Class tp produce calibration data
// too be sendt to the HLT file exhange server
// and to the PHOS HLT monitoring GUI
//


class AliHLTPHOSRcuHistogramProducer;
class AliHLTPHOSRcuCellAccumulatedEnergyDataStruct;
class AliHLTPHOSSharedMemoryInterfacev2;
class AliHLTPHOSChannelDataHeaderStruct;

class AliHLTPHOSRcuHistogramProducerComponent:public AliHLTCaloProcessor
{
 public:
  AliHLTPHOSRcuHistogramProducerComponent();
  virtual ~AliHLTPHOSRcuHistogramProducerComponent();
  virtual int DoInit( int argc, const char** argv );
  virtual int Deinit();
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn();
  virtual const char* GetComponentID();

 protected:
  //  using AliHLTPHOSRcuProcessor::DoEvent;

 private:
  AliHLTPHOSRcuHistogramProducerComponent(const AliHLTPHOSRcuHistogramProducerComponent &);
  AliHLTPHOSRcuHistogramProducerComponent & operator = (const AliHLTPHOSRcuHistogramProducerComponent &);
  int fHistoWriteFrequency;
  AliHLTPHOSRcuHistogramProducer* fRcuHistoProducerPtr;   /**<Pointer to a phos histoproducer object*/
  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct*  fOutPtr; /**<Pointer to outputbuffer to write results from the component into shared memory*/
  AliHLTPHOSSharedMemoryInterfacev2 *fShmPtr; // Interface to read altro channel data from shared memory
};

#endif
