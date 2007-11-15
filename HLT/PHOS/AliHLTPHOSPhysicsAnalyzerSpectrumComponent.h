/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                              */
//Component for create an invariant mass spectrum for pi0's

#ifndef ALIHLTPHOSPHYSICSANALYZERSPECTRUMCOMPONENT_H
#define ALIHLTPHOSPHYSICSANALYZERSPECTRUMCOMPONENT_H

// removed  PTH#include "AliHLTProcessor.h"
#include "AliHLTPHOSProcessor.h" // added by PTH
#include "AliHLTPHOSBase.h"

class TH1F;
class AliHLTPHOSPhysicsAnalyzerSpectrum;
class AliHLTPHOSPhysicsAnalyzerPeakFitter;
class Rtypes;
class AliHLTPHOSDefinitions;
class AliHLTPHOSPhysicsDefinitions;
class TFile;

struct AliHLTPHOSClusterDataStruct;


// PTH class AliHLTPHOSPhysicsAnalyzerSpectrumComponent: public AliHLTPHOSBase, public AliHLTProcessor
class AliHLTPHOSPhysicsAnalyzerSpectrumComponent: public AliHLTPHOSProcessor // added by PTH
{
 public:

  AliHLTPHOSPhysicsAnalyzerSpectrumComponent();
  ~AliHLTPHOSPhysicsAnalyzerSpectrumComponent();

  // PTH  AliHLTPHOSPhysicsAnalyzerSpectrumComponent(const AliHLTPHOSPhysicsAnalyzerSpectrumComponent &);
  //  AliHLTPHOSPhysicsAnalyzerSpectrumComponent & operator = (const AliHLTPHOSPhysicsAnalyzerSpectrumComponent &)
  //   {
  //     return *this;
  //  }
  const char* GetComponentID();
  void GetInputDataTypes(vector<AliHLTComponentDataType>& list);

  AliHLTComponentDataType GetOutputDataType();

  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /* 
  Int_t DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*,
		AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&,
		std::vector<AliHLTComponentBlockData>&);*/
  int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  AliHLTComponent* Spawn();

 protected:

  Int_t DoInit(int argc, const char** argv);
  Int_t Deinit();
  Int_t DoDeinit();
  using AliHLTPHOSProcessor::DoEvent;
  
 private:
  
  AliHLTPHOSPhysicsAnalyzerSpectrum* fAnalyzerPtr;                //! /**<Pointer to spectrum analyzer*
  AliHLTPHOSPhysicsAnalyzerPeakFitter* fPeakFitter;               //! /**<Pointer to peak fitter*/
  TH1F* fRootHistPtr;                                             //! /**<Pointer to histogram*/
  AliHLTPHOSClusterDataStruct* fClusterArrayPtr[10000];           //! /**<Pointer to array of clusters*/
  Int_t fWriteInterval;                                               /**<Interval for writing to disk*/

  static const AliHLTComponentDataType fgkInputDataTypes[];           /**<Data types*/                          
  static UInt_t fgCount;                                              /**<Event count*/
};

#endif
