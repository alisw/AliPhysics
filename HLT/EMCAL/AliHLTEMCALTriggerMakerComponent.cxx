/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Markus Fasel                                          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloTriggerPatchContainerStruct.h"
#include "AliHLTCaloTriggerPatchDataStruct.h"

#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALTriggerMaker.h"
#include "AliHLTEMCALTriggerMakerComponent.h"

ClassImp(AliHLTEMCALTriggerMakerComponent)

//AliHLTEMCALTriggerMakerComponent gAliHLTEMCALTriggerMakerComponent;

AliHLTEMCALTriggerMakerComponent::AliHLTEMCALTriggerMakerComponent() :
AliHLTCaloProcessor(),
fTriggerMakerPtr(NULL),
fTriggerPatchPtr(NULL),
fGeometry(NULL)
{
}

AliHLTEMCALTriggerMakerComponent::~AliHLTEMCALTriggerMakerComponent() {
  if(fTriggerMakerPtr) delete fTriggerMakerPtr;
  if(fGeometry) delete fGeometry;
}

int AliHLTEMCALTriggerMakerComponent::DoEvent ( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
                AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
                std::vector<AliHLTComponentBlockData>& outputBlocks ){
  if(!blocks) {
    return 0;
  }

  //patch in order to skip calib events
  if(! IsDataEvent()){
    return 0;
  }

  //see header file for documentation
  UInt_t offset           = 0;
  UInt_t mysize           = 0;
  Int_t digitCount        = 0;
  Int_t ret               = 0;

  UInt_t specification = 0;

  fTriggerMakerPtr->ResetADC();

  const AliHLTComponentBlockData* iter = 0;
  AliHLTCaloDigitDataStruct *digit = NULL;

  Int_t nDigits = 0;
  for(Int_t ndx = 0; ndx < evtData.fBlockCnt; ndx++){
    iter = blocks + ndx;

    if(iter->fDataType != AliHLTEMCALDefinitions::fgkDigitDataType) {
      HLTDebug("Invalid data type received");
      continue;
    }


    specification |= iter->fSpecification;
    digit = reinterpret_cast<AliHLTCaloDigitDataStruct*>(iter->fPtr);
    nDigits = iter->fSize/sizeof(AliHLTCaloDigitDataStruct);
    HLTDebug("Block %d received %d digits", ndx, nDigits);
    for(Int_t idigit = 0; idigit < nDigits; idigit++){
      fTriggerMakerPtr->AddDigit(digit);
      digit++;
    }
  }

  fTriggerMakerPtr->SetTriggerPatchDataPtr(reinterpret_cast<AliHLTCaloTriggerPatchDataStruct *>(outputPtr));
  Int_t npatches = fTriggerMakerPtr->FindPatches(size);
  HLTDebug("Found %d patches", npatches);
  mysize = npatches*sizeof(AliHLTCaloTriggerPatchDataStruct);

  if(mysize != 0){
    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset = offset;
    bd.fSize = mysize;
    bd.fDataType = AliHLTEMCALDefinitions::fgkTriggerPatchDataType;
    bd.fSpecification = specification;
    outputBlocks.push_back(bd);
  }
  size = mysize;
  return 0;
}

const char* AliHLTEMCALTriggerMakerComponent::GetComponentID(){
  //See headerfile for documentation
  return "EmcalTriggerMaker";
}

void AliHLTEMCALTriggerMakerComponent::GetInputDataTypes(std::vector<AliHLTComponentDataType>& list){
  list.clear();
  list.push_back(AliHLTEMCALDefinitions::fgkDigitDataType);
}

AliHLTComponentDataType AliHLTEMCALTriggerMakerComponent::GetOutputDataType(){
  return AliHLTEMCALDefinitions::fgkTriggerPatchDataType;
}

void AliHLTEMCALTriggerMakerComponent::GetOutputDataSize ( unsigned long& constBase, double& inputMultiplier ){
  constBase = 1000 *(float)sizeof(AliHLTCaloTriggerPatchDataStruct);
  inputMultiplier = 0; // (float)sizeof(AliHLTCaloTriggerPatchDataStruct)/sizeof(AliHLTCaloDigitDataStruct)+1;
}

AliHLTComponent* AliHLTEMCALTriggerMakerComponent::Spawn(){
  return new AliHLTEMCALTriggerMakerComponent();
}


int AliHLTEMCALTriggerMakerComponent::DoInit ( int argc, const char** argv ){
  InitialiseGeometry();
  fTriggerMakerPtr = new AliHLTEMCALTriggerMaker;
  fTriggerMakerPtr->Initialise(fGeometry);
  return 0;
}

int AliHLTEMCALTriggerMakerComponent::Deinit(){
  if(fTriggerMakerPtr){
    delete fTriggerMakerPtr;
    fTriggerMakerPtr = 0;
  }
  if(fGeometry){
    delete fGeometry;
    fGeometry = 0;
  }
  return 0;
}

void AliHLTEMCALTriggerMakerComponent::InitialiseGeometry(){
  fGeometry = new AliHLTEMCALGeometry(GetRunNo());
}
