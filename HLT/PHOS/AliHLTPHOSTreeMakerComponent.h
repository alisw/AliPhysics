
 /**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/



#ifndef ALIHLTPHOSTREEMAKERCOMPONENT_H
#define ALIHLTPHOSTREEMAKERCOMPONENT_H

#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSTreeMaker;
class TTree;

class AliHLTPHOSTreeMakerComponent : public AliHLTPHOSProcessor
{
 public:
  
  AliHLTPHOSTreeMakerComponent();
  virtual ~AliHLTPHOSTreeMakerComponent();
    
  const char* GetComponentID();
  
  void  GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);
  
  AliHLTComponentDataType GetOutputDataType();
  
  void GetOutputDataSize(unsigned long& constBase, double& inputmultiplier);
/*
  int DoEvent(const AliHLTComponentEventData&,
	      AliHLTComponentTriggerData&);
  */
  
  int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
	      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t&size,
       std::vector<AliHLTComponentBlockData>& outputBlocks);
  
  AliHLTComponent* Spawn();

  void Write(); 
  void ResetTrees();
   

 protected:
  using AliHLTPHOSProcessor::DoEvent;
  int DoInit(int argc, const char** argv);

  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor
  
 private:  
  AliHLTPHOSTreeMaker *fTreeMakerPtr; //comment
  TTree *fDigitTreePtr; //comment
  UInt_t fEventCount; //comment
  UInt_t fWriteInterval; //comment
  UInt_t fRunNb; //comment
  char *fDirectory; //comment
 
  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type
  
};
#endif
  
