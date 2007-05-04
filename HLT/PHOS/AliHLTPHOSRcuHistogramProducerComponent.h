#ifndef ALIHLTPHOSRCUHISTOGRAMPRODUCERCOMPONENT_H
#define ALIHLTPHOSRCUHISTOGRAMPRODUCERCOMPONENT_H 

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice  */ 

#include "AliHLTProcessor.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"
//#include "AliHLTPHOSModuleCellAccumulatedEnergyDataStruct.h"
#include "AliHLTPHOSRcuCellAccumulatedEnergyDataStruct.h"
#include "Rtypes.h"

class AliHLTPHOSRcuHistogramProducer;

class AliHLTPHOSRcuHistogramProducerComponent:public AliHLTProcessor
{
 public:
  AliHLTPHOSRcuHistogramProducerComponent();
  virtual ~AliHLTPHOSRcuHistogramProducerComponent();
  virtual int DoInit( int argc, const char** argv );
  virtual int Deinit();
  virtual int DoDeinit();
  //  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*,
  //		      AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

  //  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual void GetInputDataTypes( std::vector <AliHLTComponentDataType>&);

  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn();

  virtual const char* GetComponentID();
  int GetEquippmentId();
  void SetRcuX(AliHLTUInt8_t X);
  void SetRcuZ(AliHLTUInt8_t Z);
  void SetModuleID(AliHLTUInt8_t moduleID);
  void SetEquippmentId(int id);


 private:
  AliHLTPHOSRcuHistogramProducerComponent(const AliHLTPHOSRcuHistogramProducerComponent & );
  AliHLTPHOSRcuHistogramProducerComponent & operator = (const AliHLTPHOSRcuHistogramProducerComponent &)
   {
      return *this;
   };
  int fEventCount;
  AliHLTUInt8_t fRcuX;
  AliHLTUInt8_t fRcuZ; 
  AliHLTUInt8_t fModuleID;
  AliHLTPHOSRcuHistogramProducer* fRcuHistoProducerPtr;
  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct*  fOutPtr;

  //  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct*  fInnPtr;
  //  AliHLTPHOSRcuCellAccumulatedEnergyDataStruct&  fInnPtr;

  static const AliHLTComponentDataType inputDataTypes[];
  static const AliHLTComponentDataType outputDataType;

  int fEquippmentID; 
};

#endif
