/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALDigitsMonitor.h"
#include "AliHLTEMCALDigitsMonitorComponent.h"

#include <TH1.h>
#include <TObjArray.h>
#include <TString.h>
#include <TFile.h>

AliHLTEMCALDigitsMonitorComponent::AliHLTEMCALDigitsMonitorComponent():
    AliHLTCaloProcessor(),
	fLocalEventCount(0),
    fHistoResetOnPush(kTRUE),
    fDigitsMonitor(0x0)

{

}

AliHLTEMCALDigitsMonitorComponent::~AliHLTEMCALDigitsMonitorComponent(){

}
    
const char* AliHLTEMCALDigitsMonitorComponent::GetComponentID(){
	//see header file for documentation
    return "EmcalDigitsMonitor";
}

void AliHLTEMCALDigitsMonitorComponent::GetInputDataTypes(std::vector<AliHLTComponentDataType>& list){
	//see header file for documentation
	list.clear();
	list.push_back(AliHLTEMCALDefinitions::fgkDigitDataType);
}

AliHLTComponentDataType AliHLTEMCALDigitsMonitorComponent::GetOutputDataType(){
	//see header file for documentation
	return kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL;
}

void AliHLTEMCALDigitsMonitorComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier){
    //constBase = 250000000;
    constBase = 250000;
    // to be reviewed later
    inputMultiplier = 0;
}

int AliHLTEMCALDigitsMonitorComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/,
			std::vector<AliHLTComponentBlockData>& /*outputBlocks*/){
    if(!blocks) {
        return 0;
    }

    //patch in order to skip calib events
    if(! IsDataEvent()){
        return 0;
    }
    
    AliHLTCaloDigitDataStruct *digits = NULL;
    const AliHLTComponentBlockData* iter = 0;
    Int_t nDigits = 0;
    for(Int_t ndx = 0; ndx < evtData.fBlockCnt; ndx++){
        iter = blocks + ndx;

        if(!CheckInputDataType(iter->fDataType)){
            continue;
        }

        if(iter->fDataType == AliHLTEMCALDefinitions::fgkDigitDataType) {
            digits = reinterpret_cast<AliHLTCaloDigitDataStruct*>(iter->fPtr);
            nDigits = iter->fSize/sizeof(AliHLTCaloDigitDataStruct);
            HLTDebug("Block %d received %d digits\n", ndx, nDigits);
            fDigitsMonitor->ProcessDigits(nDigits, digits);
        }
    }
    fLocalEventCount += 1;

    PushHistograms(fDigitsMonitor->GetListOfHistograms());
    return 0;
}

void AliHLTEMCALDigitsMonitorComponent::PushHistograms(TCollection* list){
  TIter next(list);

  TH1* histo = 0;

  while ((histo = static_cast<TH1*>(next()))) {
    if (histo->GetEntries() > 0) {
      Int_t nbytes = PushBack(histo, kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL , 0);
      if (fHistoResetOnPush && nbytes > 0) {
        histo->Reset();
      }
    }
  }
}

bool AliHLTEMCALDigitsMonitorComponent::CheckInputDataType(const AliHLTComponentDataType &datatype){
  vector <AliHLTComponentDataType> validTypes;
  GetInputDataTypes(validTypes);

  for(UInt_t i = 0; i < validTypes.size(); i++) {
    if (datatype == validTypes.at(i)) {
      return true;
    }
  }

  HLTDebug("Invalid Datatype");
  return false;
}

AliHLTComponent* AliHLTEMCALDigitsMonitorComponent::Spawn() {
    return new AliHLTEMCALDigitsMonitorComponent();
}

int AliHLTEMCALDigitsMonitorComponent::DoInit(int argc, const char** argv){
    fDigitsMonitor = new AliHLTEMCALDigitsMonitor;

}

int AliHLTEMCALDigitsMonitorComponent::Deinit(){
    if(fDigitsMonitor) delete fDigitsMonitor;
} 