

#ifndef ALIHLTPHOSRAWANALYZERCOMPONENT_H
#define ALIHLTPHOSRAWANALYZERCOMPONENT_H

/* Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

//#include <Rtypes.h>
//#include "TObject.h"
#include "AliHLTProcessor.h"
#include "AliHLTPHOSRawAnalyzer.h"
//#include "AliHLTPHOSRawAnalyzerCrude.h"

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

  virtual const char* GetComponentID() = 0;

  virtual void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);
  virtual AliHLTComponentDataType GetOutputDataType();
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  //  virtual void GetOutputDataSize(long  int&, double&);

  //  virtual AliHLTComponent* Spawn();
  virtual AliHLTComponent* Spawn() = 0;

  virtual int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*, AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&, std::vector<AliHLTComponentBlockData, std::allocator<AliHLTComponentBlockData> >&);

  // private:
 protected:
 AliHLTPHOSRawAnalyzer *analyzerPtr; 

 // ClassDef(AliHLTPHOSRawAnalyzerComponent, 2) 
 private:
 int eventCount;
  static const AliHLTComponentDataType inputDataTypes[];
  static const AliHLTComponentDataType outputDataType;

};
#endif
