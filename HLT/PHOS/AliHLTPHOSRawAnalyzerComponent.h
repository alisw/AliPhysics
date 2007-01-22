

#ifndef ALIHLTPHOSRAWANALYZERCOMPONENT_H
#define ALIHLTPHOSRAWANALYZERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include "AliHLTProcessor.h"
#include "AliHLTPHOSRawAnalyzer.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStream.h"
#include "AliHLTPHOSDefinitions.h"
#include "TH2.h"


/*
#include "AliHLTTPCRawDataUnpackerComponent.h"
#include "AliTPCRawStream.h"
#include "AliRawDataHeader.h"
#include "AliRawReaderMemory.h"
#include "AliHLTTPCRawDataFormat.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCTransform.h"
#include <stdlib.h>
#include <errno.h>
*/

class AliHLTPHOSRawAnalyzerComponent: public AliHLTProcessor
{
 public:
  AliHLTPHOSRawAnalyzerComponent();
  ~AliHLTPHOSRawAnalyzerComponent();
  AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & );
  AliHLTPHOSRawAnalyzerComponent & operator = (const AliHLTPHOSRawAnalyzerComponent &)
   {
      return *this;
   };



  virtual int DoInit( int argc, const char** argv );
  virtual int Deinit();
  virtual int DoDeinit();
  void DumpData();

  virtual const char* GetComponentID() = 0;

  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  virtual AliHLTComponent* Spawn() = 0;

  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);

  // private:
 protected:
 AliHLTPHOSRawAnalyzer *analyzerPtr; 
 void Reset();
 void ResetDataPtr();

 private:
 

 int eventCount;


 //	AliRawReaderMemory *fRawMemoryReader;
 //	AliTPCRawStream *fTPCRawStream;
 Double_t fTmpChannelData[1008];
 // Int_t fMaxValues[5][64][56][2];
 Double_t fMaxValues[5][64][56][2];

  //  Int_t fMaxValuesLG[5][64][56][2];
  TH2S *legoPlotPtr;
  //  TH2S *legoPlotLgPtr;
 AliCaloRawStream *fPHOSRawStream;
 AliRawReaderMemory *fRawMemoryReader;
 static const AliHLTComponentDataType inputDataTypes[];
 static const AliHLTComponentDataType outputDataType;

};
#endif
