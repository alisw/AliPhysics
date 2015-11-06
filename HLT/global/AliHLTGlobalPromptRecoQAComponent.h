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


protected:

  int fVerbosity; //!transient
  AliHLTComponentBenchmark fBenchmark; // benchmark
  
  Int_t fSkipEvents;
  Int_t fPrintStats;
  Int_t fEventsSinceSkip;

  TH2I* fHistSPDclusters_SPDrawSize;
  TH2I* fHistSSDclusters_SSDrawSize;
  TH2I* fHistSDDclusters_SDDrawSize;
  TH2I* fHistITSSAtracks_SPDclusters;
  TH2I* fHistSPDclusters_SSDclusters;
  TH2F* fHistTPCHLTclusters_TPCCompressionRatio;
  TH2I* fHistTPCtracks_TPCtracklets;
  TH2I* fHistITStracks_ITSOutTracks;
  TH2I* fHistTPCClusterSize_TPCCompressedSize;

  ClassDef(AliHLTGlobalPromptRecoQAComponent, 0)
};
#endif
