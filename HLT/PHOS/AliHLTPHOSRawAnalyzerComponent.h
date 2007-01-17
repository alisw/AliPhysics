/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#ifndef ALIHLTPHOSANALYZERCOMPONENT_H
#define ALIHLTPHOSANALYZERCOMPONENT_H

//#include <Rtypes.h>
//#include "TObject.h"
#include "AliHLTProcessor.h"
#include "AliHLTPHOSRawAnalyzer.h"
#include "AliHLTPHOSRawAnalyzerCrude.h"

class AliHLTPHOSRawAnalyzerComponent: public AliHLTProcessor
{
 public:
  AliHLTPHOSRawAnalyzerComponent();
  ~AliHLTPHOSRawAnalyzerComponent();
  AliHLTPHOSRawAnalyzerComponent(const AliHLTPHOSRawAnalyzerComponent & );
  AliHLTPHOSRawAnalyzerComponent & operator = (const AliHLTPHOSRawAnalyzerComponent)
   {
      return *this;
   };

  virtual int Deinit();
  virtual int DoDeinit();
  virtual const char* GetComponentID();
  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(long unsigned int&, double&);
  virtual void GetOutputDataSize(long  int&, double&);
  virtual AliHLTComponent* Spawn();
  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);

  // private:
 protected:
 AliHLTPHOSRawAnalyzer *analyzerPtr; 

 // ClassDef(AliHLTPHOSRawAnalyzerComponent, 2) 
};
#endif
