
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

#ifndef ALIHLTPHOSMIPCOUNTERCOMPONENT_H
#define ALIHLTPHOSMIPCOUNTERCOMPONENT_H

#include "AliHLTPHOSProcessor.h"

class AliHLTPHOSMIPCounter;
class TH1F;
class TH1I;
class TH2I;

/**
	@author Øystein Djuvsland <oystein.djuvsland@gmail.com>
*/
class AliHLTPHOSMIPCounterComponent : public AliHLTPHOSProcessor
{
public:
    AliHLTPHOSMIPCounterComponent();

    virtual ~AliHLTPHOSMIPCounterComponent();

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
  
  Int_t fEvtCnt;          //comment
  Int_t fInterval; //comment
  Int_t fMIPCount; //comment
  Char_t *fTRUThreshold; //comment
  Float_t fMIPCountInterval; //comment
  Char_t *fPath; //comment
  AliHLTPHOSMIPCounter *fMIPCounterPtr; //comment
  TH1I *fHistPtr; //comment
  TH1I *fIntervalHistPtr; //comment
  TH1F *fRateHistPtr; //comment
  TH2I *fChannelHistPtr; //comment
  TH1F *fRatioHistPtr; //comment
  
  static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type
};
#endif

    
