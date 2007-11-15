/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef ALIHLTPHOSBASELINEANALYZERCOMPONENT_H
#define ALIHLTPHOSBASELINEANALYZERCOMPONENT_H

#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSBaselineAnalyzer;
class TTree;

class AliHLTPHOSBaselineAnalyzerComponent : public AliHLTPHOSProcessor
{
public:
  AliHLTPHOSBaselineAnalyzerComponent();
  virtual ~AliHLTPHOSBaselineAnalyzerComponent();
  
  const char* GetComponentID();

  void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

  AliHLTComponentDataType GetOutputDataType();

  void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
	      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
	      std::vector<AliHLTComponentBlockData>& outputBlocks);
  
  AliHLTComponent* Spawn();
  
protected:
  using AliHLTPHOSProcessor::DoEvent;
  int DoInit(int argc, const char** argv);
  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor
  
private:
  
  void CalculateAll();

  AliHLTPHOSBaselineAnalyzer *fBaselineAnalyzerPtr; //comment
  TTree *fTreePtr; //comment
  TClonesArray *fBaselineArrayPtr; //comment
  UInt_t fEvCnt; //comment
  UInt_t fWriteInterval; //comment
  UInt_t fFillInterval; //comment
  char *fFilename; //comment
  char* fDirectory; //comment
  char* fHistPath; //comment
  Int_t fRunNb; //comment
  Bool_t fCalculateAll; //comment

  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type

};
#endif
