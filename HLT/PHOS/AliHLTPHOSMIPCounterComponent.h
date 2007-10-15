
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

    ~AliHLTPHOSMIPCounterComponent();

    const char* GetComponentID();

    void GetInputDataTypes(std::vector<AliHLTComponentDataType>& list);

    AliHLTComponentDataType GetOutputDataType();

    void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

    int DoEvent(const AliHLTComponentEventData&, const AliHLTComponentBlockData*,
		AliHLTComponentTriggerData&, AliHLTUInt8_t*, AliHLTUInt32_t&,
		std::vector<AliHLTComponentBlockData>&);

    AliHLTComponent* Spawn();
  
  protected:
    int DoInit(int argc, const char** argv);
    virtual int Deinit(); ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor

  private:

    Int_t fEventCount;
    Int_t fInterval;
    Int_t fMIPCount;
    Char_t *fTRUThreshold;
    Float_t fMIPCountInterval;
    Char_t *fPath;
    AliHLTPHOSMIPCounter *fMIPCounterPtr;
    TH1I *fHistPtr;
    TH1I *fIntervalHistPtr;
    TH1F *fRateHistPtr;
    TH2I *fChannelHistPtr;
    TH1F *fRatioHistPtr;
    static const AliHLTComponentDataType fgkInputDataTypes[];     //HLT input data type
};
#endif

    
