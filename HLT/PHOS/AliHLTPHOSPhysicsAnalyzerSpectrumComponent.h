
#ifndef ALIHLTPHOSPHYSICSANALYZERSPECTRUMCOMPONENT
#define ALIHLTPHOSPHYSICSANALYZERSPECTRUMCOMPONENT

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliHLTProcessor.h"
#include "AliHLTPHOSPhysicsAnalyzerPeakFitter.h"
#include "AliHLTPHOSPhysicsAnalyzerSpectrum.h"
#include "AliHLTPHOSClusterDataStruct.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Rtypes.h"


class AliHLTPHOSPhysicsAnalyzerSpectrumComponent: public AliHLTProcessor
{
 public:

  AliHLTPHOSPhysicsAnalyzerSpectrumComponent();
  ~AliHLTPHOSPhysicsAnalyzerSpectrumComponent();
  AliHLTPHOSPhysicsAnalyzerSpectrumComponent(const AliHLTPHOSPhysicsAnalyzerSpectrumComponent &);
  AliHLTPHOSPhysicsAnalyzerSpectrumComponent & operator = (const AliHLTPHOSPhysicsAnalyzerSpectrumComponent &)
    {
      return *this;
    }
  const char* GetComponentID();
  void GetInputDataTypes(std::vector<AliHLTComponentDataType, std::allocator<AliHLTComponentDataType> >&);

  AliHLTComponentDataType GetOutputDataType();

  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);


  Int_t DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*,
		AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&,
		std::vector<AliHLTComponentBlockData>&);
  

  AliHLTComponent* Spawn();

 protected:

  Int_t DoInit(int argc, const char** argv);
  Int_t Deinit();
  Int_t DoDeinit();
  
 private:
  
  AliHLTPHOSPhysicsAnalyzerSpectrum* fAnalyzerPtr;
  AliHLTPHOSPhysicsAnalyzerPeakFitter* fPeakFitter;
  TH1F* fRootHistPtr;
  AliHLTPHOSClusterDataStruct* fClusterArrayPtr[10000];
  Int_t fWriteInterval;

  static const AliHLTComponentDataType inputDataTypes[];
  static int fEventCount;

};

#endif
