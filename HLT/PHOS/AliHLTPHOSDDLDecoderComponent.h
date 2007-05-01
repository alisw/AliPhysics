#ifndef ALIHLTPHOSDDLDECODERCOMPONENT_H
#define ALIHLTPHOSDDLDECODERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include "AliHLTProcessor.h"
#include "AliHLTPHOSRawAnalyzer.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSRcuChannelDataStruct.h"

class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSRcuChannelDataStruct;

class AliHLTPHOSDDLDecoderComponent:public AliHLTProcessor

{
 public:


 AliHLTPHOSDDLDecoderComponent();
  virtual ~AliHLTPHOSDDLDecoderComponent();
 AliHLTPHOSDDLDecoderComponent(const AliHLTPHOSDDLDecoderComponent & );
 AliHLTPHOSDDLDecoderComponent & operator = (const AliHLTPHOSDDLDecoderComponent &)
   {
      return *this;
   };

  virtual int DoInit( int argc, const char** argv );
  virtual int Deinit();
  virtual int DoDeinit();
  void SetEquippmentID(AliHLTUInt16_t id);
  AliHLTUInt16_t  GetEquippmentID();
  void SetCoordinates(AliHLTUInt16_t equippmentID);
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn();

  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);

 private:
  static int fgEventCount;
  AliHLTUInt16_t fEquippmentID;
  AliHLTUInt8_t  fRcuX;
  AliHLTUInt8_t  fRcuZ;
  AliHLTUInt8_t  fRcuZOffset;
  AliHLTUInt8_t  fRcuXOffset;
  AliHLTUInt8_t  fModuleID;
  Bool_t fSendChannelData;
  Bool_t fPrintInfo;
  Double_t fTmpChannelData[ALTRO_MAX_SAMPLES];
  int fPrintInfoFrequncy;
  AliCaloRawStream *fPHOSRawStream;
  AliRawReaderMemory *fRawMemoryReader;
  AliHLTPHOSRcuChannelDataStruct*  fOutPtr;
  static const AliHLTComponentDataType fgkInputDataTypes[];
};
#endif

