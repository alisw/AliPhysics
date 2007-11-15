#ifndef ALIHLTPHOSRCUHISTOGRAMPRODUCERCOMPONENT_H
#define ALIHLTPHOSRCUHISTOGRAMPRODUCERCOMPONENT_H 

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 


#include "AliHLTPHOSRcuProcessor.h"
//#include "AliHLTPHOSDefinitions.h"
//#include "AliHLTPHOSCommonDefs.h"
//#include "AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h"
//#include "Rtypes.h"



class AliHLTPHOSRcuHistogramProducer;
class AliHLTPHOSRcuCellAccumulatedEnergyDataStruct;

class AliHLTPHOSRcuHistogramProducerComponent:public AliHLTPHOSRcuProcessor
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
  using AliHLTPHOSRcuProcessor::DoEvent;

 private:
  int fHistoWriteFrequency;

  /*
  AliHLTPHOSRcuHistogramProducerComponent(const AliHLTPHOSRcuHistogramProducerComponent & );
  AliHLTPHOSRcuHistogramProducerComponent & operator = (const AliHLTPHOSRcuHistogramProducerComponent &)
   {
      return *this;
   };
  */

  AliHLTPHOSRcuHistogramProducer* fRcuHistoProducerPtr;   /**<Pointer to a phos histoproducer object*/
  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct*  fOutPtr; /**<Pointer to outputbuffer to write results from the component into shared memory*/
  //  static const AliHLTComponentDataType fgkIinputDataTypes[];
  //  static const AliHLTComponentDataType fgkOutputDataType;   /**<List of  datatypes that can be given to this component*/
};

#endif
