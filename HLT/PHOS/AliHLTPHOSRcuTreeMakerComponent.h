
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



#ifndef ALIHLTPHOSRCUTREEMAKERCOMPONENT_H
#define ALIHLTPHOSRCUTREEMAKERCOMPONENT_H

# include "AliHLTPHOSRcuProcessor.h"

class AliHLTPHOSRcuTreeMaker;
class TTree;

class AliHLTPHOSRcuTreeMakerComponent : public AliHLTPHOSRcuProcessor
{
 public:
  
  AliHLTPHOSRcuTreeMakerComponent();
  ~AliHLTPHOSRcuTreeMakerComponent();
    
  const char* GetComponentID();
  
  void  GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);
  
  AliHLTComponentDataType GetOutputDataType();
  
  void GetOutputDataSize(unsigned long& constBase, double& inputmultiplier);
/*
  int DoEvent(const AliHLTComponentEventData&,
	      AliHLTComponentTriggerData&);
  */
  
  int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*,
	      AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&,
       std::vector<AliHLTComponentBlockData>&);
  
  AliHLTComponent* Spawn();

  void Write(); 
  void ResetTrees();
   
 protected:
  int DoInit(int argc, const char** argv);

  virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor
  
 private:  
  AliHLTPHOSRcuTreeMaker *fTreeMakerPtr;
  TTree *fDigitTreePtr;
  UInt_t fEventCount;
  UInt_t fWriteInterval;
  UInt_t fRunNb;
  char *fDirectory;
 
  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type
  
};
#endif
  
