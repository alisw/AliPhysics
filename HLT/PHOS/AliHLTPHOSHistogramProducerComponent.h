#ifndef ALIHLTPHOSHISTOGRAMPRODUCERCOMPONENT_H
#define ALIHLTPHOSHISTOGRAMPRODUCERCOMPONENT_H 


/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTPHOSProcessor.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"


class AliHLTPHOSModuleCellAccumulatedEnergyDataStruct;

class AliHLTPHOSHistogramProducerComponent:public AliHLTPHOSProcessor
{
 public:
  AliHLTPHOSHistogramProducerComponent();
  virtual ~AliHLTPHOSHistogramProducerComponent();
  AliHLTPHOSHistogramProducerComponent(const AliHLTPHOSHistogramProducerComponent & );
  AliHLTPHOSHistogramProducerComponent & operator = (const AliHLTPHOSHistogramProducerComponent &)
   {
      return *this;
   };
  virtual int DoInit( int argc = 0, const char** argv = 0);
  virtual int Deinit();
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

  void DumpData(int gain = 0);
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes( vector <AliHLTComponentDataType>& list);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn();
 protected:
  void Reset();
  void ResetDataPtr();

 private:
  Double_t fEnergyAverageValues[N_ZROWS_MOD][N_XCOLUMNS_MOD][N_GAINS]; /**<Accumulated energy divided by the number of hits for each readout channel*/  
  Double_t fAccumulatedValues[N_ZROWS_MOD][N_XCOLUMNS_MOD][N_GAINS];   /**<Accumulated energy for each readout channel of one RCU*/
  //  Double_t fTimingAverageValues[N_ZROWS_MOD][N_XCOLUMNS_MOD][N_GAINS]; 
  AliHLTUInt32_t fHits[N_ZROWS_MOD][N_XCOLUMNS_MOD][N_GAINS];         /**<Total number of hits for each cell of one RCU*/
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];                        /**<Array to temporarily store dat fro a single altro channel*/                        
  AliHLTPHOSModuleCellAccumulatedEnergyDataStruct*  fOutPtr;          /**<Pointer to outputbuffer to write results from the component into shared memory*/
  static const AliHLTComponentDataType fgkInputDataTypes[];           /**<List of  datatypes that can be given to this component*/  
  static const AliHLTComponentDataType fgkOutputDataType;             /**<Output datatype produced by this component*/
};

#endif
