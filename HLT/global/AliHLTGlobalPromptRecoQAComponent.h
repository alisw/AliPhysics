#ifndef ALIHLTGLOBALPROMPTRECOQACOMPONENT_H
#define ALIHLTGLOBALPROMPTRECOQACOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"
#include <vector>

// forward declarations
class AliESDEvent;
class AliESDfriend;
class TTree;
struct AliHLTTracksData;
class AliTPCclusterMI;
class TH2I;
class TH2F;
class AliHLTTPCHWCFData;

/**
 * @class AliHLTGlobalPromptRecoQAComponent
 * simple global data QA
 *
 */
class AliHLTGlobalPromptRecoQAComponent : public AliHLTProcessor
{
 public:
  /** standard constructor */
  AliHLTGlobalPromptRecoQAComponent();
  /** destructor */
  virtual ~AliHLTGlobalPromptRecoQAComponent();

  // interface methods of base class
  const char* GetComponentID() {return "PromptRecoQA";};
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  AliHLTComponent* Spawn() {return new AliHLTGlobalPromptRecoQAComponent;}

 protected:
  // interface methods of base class
  int DoInit(int argc, const char** argv);
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData,
		       const AliHLTComponentBlockData* blocks, 
		       AliHLTComponentTriggerData& trigData,
		       AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size,
		       AliHLTComponentBlockDataList& outputBlocks );

  

  using AliHLTProcessor::DoEvent;

 private:
  /** copy constructor prohibited */
  AliHLTGlobalPromptRecoQAComponent(const AliHLTGlobalPromptRecoQAComponent&);
  /** assignment operator prohibited */
  AliHLTGlobalPromptRecoQAComponent& operator=(const AliHLTGlobalPromptRecoQAComponent&);

  /**
   * (Re)Configure from the CDB
   * Loads the following objects:
   */
  int Reconfigure(const char* cdbEntry, const char* chainId);

  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */
  int Configure(const char* arguments);
  
  //Root cannot do templates...
  //template <class T> void FillHist(int check, T* hist, S val1, U val2, int& flag);

protected:

  int fVerbosity; //!transient
  AliHLTComponentBenchmark fBenchmark; // benchmark
  
  AliHLTTPCHWCFData* fpHWCFData;
  
  Int_t fSkipEvents;
  Int_t fPrintStats; //print status messages: 0: never, 1: when pushing histograms (respect pushback-period), 2: always
  Int_t fEventsSinceSkip;

  TH2I* fHistSPDclusters_SPDrawSize;
  TH2I* fHistSSDclusters_SSDrawSize;
  TH2I* fHistSDDclusters_SDDrawSize;
  TH2I* fHistITSSAtracks_SPDclusters;
  TH2I* fHistSPDclusters_SSDclusters;
  TH2F* fHistTPCHLTclusters_TPCCompressionRatio;
  TH2F* fHistTPCHLTclusters_TPCFullCompressionRatio;
  TH2F* fHistHLTSize_HLTInOutRatio;
  TH2I* fHistTPCtracks_TPCtracklets;
  TH2I* fHistITStracks_ITSOutTracks;
  TH2I* fHistTPCClusterSize_TPCCompressedSize;
  TH2I* fHistTPCRawSize_TPCCompressedSize;
  TH2I* fHistHLTInSize_HLTOutSize;
  TH2F* fHistZNA_VZEROTrigChargeA;
  TH2F* fHistZNC_VZEROTrigChargeC;
  TH2F* fHistZNT_VZEROTrigChargeT;
  TH2F* fHistVZERO_SPDClusters;
  TH2F* fHistVZERO_ITSSAPTracks;

  ClassDef(AliHLTGlobalPromptRecoQAComponent, 0)
};
#endif
