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
#include "AliHLTCaloTriggerHeaderStruct.h"
#include "AliHLTCaloTriggerDataStruct.h"
#include "AliHLTCaloTriggerPatchDataStruct.h"

#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALTriggerMaker.h"
#include "AliHLTEMCALTriggerMakerComponent.h"

ClassImp(AliHLTEMCALTriggerMakerComponent)

//AliHLTEMCALTriggerMakerComponent gAliHLTEMCALTriggerMakerComponent;

AliHLTEMCALTriggerMakerComponent::AliHLTEMCALTriggerMakerComponent() :
AliHLTCaloProcessor(),
fTriggerMakerPtrCells(NULL),
fTriggerMakerPtrFastor(NULL),
fGeometry(NULL)
{
}

AliHLTEMCALTriggerMakerComponent::~AliHLTEMCALTriggerMakerComponent() {
  if(fTriggerMakerPtrCells) delete fTriggerMakerPtrCells;
  if(fTriggerMakerPtrFastor) delete fTriggerMakerPtrFastor;
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

  fTriggerMakerPtrCells->ResetADC();
  fTriggerMakerPtrFastor->ResetADC();

  const AliHLTComponentBlockData* iter = 0;
  AliHLTCaloDigitDataStruct *digit = NULL;
  AliHLTCaloTriggerHeaderStruct *triggerhead = NULL;
  AliHLTCaloTriggerDataStruct *dataptr = NULL;

  Int_t nDigits = 0, nDigitsGlob = 0, nfastor = 0;
  for(Int_t ndx = 0; ndx < evtData.fBlockCnt; ndx++){
    iter = blocks + ndx;

    if(!CheckInputDataType(iter->fDataType)){
      continue;
    }

    if(iter->fDataType == AliHLTEMCALDefinitions::fgkDigitDataType) {
      digit = reinterpret_cast<AliHLTCaloDigitDataStruct*>(iter->fPtr);
      nDigits = iter->fSize/sizeof(AliHLTCaloDigitDataStruct);
      HLTDebug("Block %d received %d digits\n", ndx, nDigits);
      for(Int_t idigit = 0; idigit < nDigits; idigit++){
        fTriggerMakerPtrCells->AddDigit(digit);
        digit++;
      }
      nDigitsGlob += nDigits;
    } else if(iter->fDataType == (kAliHLTDataTypeCaloTrigger | kAliHLTDataOriginEMCAL)){
      triggerhead = reinterpret_cast<AliHLTCaloTriggerHeaderStruct *>(iter->fPtr);
      AliHLTCaloTriggerDataStruct *dataptr = reinterpret_cast<AliHLTCaloTriggerDataStruct *>(reinterpret_cast<AliHLTUInt8_t* >(iter->fPtr) + sizeof(AliHLTCaloTriggerHeaderStruct));
      HLTDebug("Block %d received %d fastor triggers", ndx, triggerhead->fNfastor);
      for(Int_t datacount = 0; datacount < triggerhead->fNfastor; datacount++) {
        fTriggerMakerPtrFastor->SetADC(dataptr->fCol, dataptr->fRow, static_cast<Float_t>(dataptr->fL1TimeSum));
        dataptr++;
      }
      nfastor += triggerhead->fNfastor;
    }
  }

  Int_t npatches = 0;
  if(nfastor){
    fTriggerMakerPtrFastor->SetTriggerPatchDataPtr(reinterpret_cast<AliHLTCaloTriggerPatchDataStruct *>(outputPtr));
    npatches = fTriggerMakerPtrFastor->FindPatches();
    HLTDebug("Found %d patches from fastors\n", npatches);
    outputPtr += npatches;
    mysize += sizeof(AliHLTCaloTriggerPatchDataStruct) * npatches;
  }

  if(nDigitsGlob){
    fTriggerMakerPtrCells->SetTriggerPatchDataPtr(reinterpret_cast<AliHLTCaloTriggerPatchDataStruct *>(outputPtr));
    npatches = fTriggerMakerPtrCells->FindPatches();
    HLTDebug("Found %d patches from cells\n", npatches);
    mysize += sizeof(AliHLTCaloTriggerPatchDataStruct) * npatches;
    outputPtr += npatches;
  }

  if(mysize != 0){
    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset = offset;
    bd.fSize = mysize;
    bd.fDataType = AliHLTEMCALDefinitions::fgkTriggerPatchDataType;
    bd.fSpecification = 0;
    outputBlocks.push_back(bd);
  }
  size = mysize;
  return 0;
}

bool AliHLTEMCALTriggerMakerComponent::CheckInputDataType(const AliHLTComponentDataType &datatype) {
  //comment
  vector <AliHLTComponentDataType> validTypes;
  GetInputDataTypes(validTypes);

  for(UInt_t i=0; i < validTypes.size(); i++) {
    if ( datatype  ==  validTypes.at(i) ) {
      return true;
    }
  }

  HLTDebug("Invalid Datatype");
  return false;
}

const char* AliHLTEMCALTriggerMakerComponent::GetComponentID(){
  //See headerfile for documentation
  return "EmcalTriggerMaker";
}

void AliHLTEMCALTriggerMakerComponent::GetInputDataTypes(std::vector<AliHLTComponentDataType>& list){
  list.clear();
  list.push_back(AliHLTEMCALDefinitions::fgkDigitDataType);
  list.push_back(kAliHLTDataTypeCaloTrigger | kAliHLTDataOriginEMCAL);
}

AliHLTComponentDataType AliHLTEMCALTriggerMakerComponent::GetOutputDataType(){
  return AliHLTEMCALDefinitions::fgkTriggerPatchDataType;
}

void AliHLTEMCALTriggerMakerComponent::GetOutputDataSize ( unsigned long& constBase, double& inputMultiplier ){
  constBase = 3000 *(float)sizeof(AliHLTCaloTriggerPatchDataStruct);
  inputMultiplier = 0; // (float)sizeof(AliHLTCaloTriggerPatchDataStruct)/sizeof(AliHLTCaloDigitDataStruct)+1;
}

AliHLTComponent* AliHLTEMCALTriggerMakerComponent::Spawn(){
  return new AliHLTEMCALTriggerMakerComponent();
}


int AliHLTEMCALTriggerMakerComponent::DoInit ( int argc, const char** argv ){
  InitialiseGeometry();
  Int_t onlinethresh[2]; memset(onlinethresh, 0, sizeof(Int_t) *2);
  Float_t offlinethresh[2]; memset(offlinethresh, 0, sizeof(Float_t) *2);
  for(Int_t iarg = 0; iarg < argc; iarg++){
    TString argstring(argv[iarg]);
    if(argstring.Contains("-gammaoffthresh")){
      offlinethresh[0] = TString(argv[iarg+1]).Atof();
    } else if(argstring.Contains("-gammaonthresh")){
      onlinethresh[0] = TString(argv[iarg+1]).Atoi();
    } else if(argstring.Contains("-jetoffthresh")){
      offlinethresh[1] = TString(argv[iarg+1]).Atof();
    } else if(argstring.Contains("-jetonthresh")){
      onlinethresh[1] = TString(argv[iarg+1]).Atoi();
    }
  }
  fTriggerMakerPtrCells = new AliHLTEMCALTriggerMaker;
  fTriggerMakerPtrCells->SetOrigin(AliHLTEMCALTriggerMaker::kOriginCELLS);
  fTriggerMakerPtrCells->SetGammaThreshold(offlinethresh[0]);
  fTriggerMakerPtrCells->SetJetThreshold(offlinethresh[1]);
  fTriggerMakerPtrCells->Initialise(fGeometry);
  fTriggerMakerPtrFastor = new AliHLTEMCALTriggerMaker;
  fTriggerMakerPtrFastor->SetOrigin(AliHLTEMCALTriggerMaker::kOriginRECAL);
  fTriggerMakerPtrCells->SetGammaThreshold(onlinethresh[0]);
  fTriggerMakerPtrCells->SetJetThreshold(onlinethresh[1]);
  fTriggerMakerPtrFastor->Initialise(fGeometry);
  return 0;
}

int AliHLTEMCALTriggerMakerComponent::Deinit(){
  if(fTriggerMakerPtrCells){
    delete fTriggerMakerPtrCells;
    fTriggerMakerPtrCells = 0;
  }
  if(fTriggerMakerPtrFastor){
    delete fTriggerMakerPtrFastor;
    fTriggerMakerPtrFastor = 0;
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
