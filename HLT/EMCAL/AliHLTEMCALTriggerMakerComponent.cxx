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

#include <TStopwatch.h>

#include <bitset>
#include <iostream>

ClassImp(AliHLTEMCALTriggerMakerComponent)

//AliHLTEMCALTriggerMakerComponent gAliHLTEMCALTriggerMakerComponent;

AliHLTEMCALTriggerMakerComponent::AliHLTEMCALTriggerMakerComponent() :
  AliHLTCaloProcessor(),
  fTriggerMakerPtr(NULL),
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

#ifdef __PROFILE__
  TStopwatch profile;
  profile.Start();
#endif

  //see header file for documentation
  UInt_t offset           = 0;
  UInt_t mysize           = 0;
  Int_t digitCount        = 0;
  Int_t ret               = 0;

  UInt_t specification = 0;

  fTriggerMakerPtr->ResetADC();

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
        fTriggerMakerPtr->AddDigit(digit);
        digit++;
      }
      nDigitsGlob += nDigits;
    } else if(iter->fDataType == (kAliHLTDataTypeCaloTrigger | kAliHLTDataOriginEMCAL)){
      triggerhead = reinterpret_cast<AliHLTCaloTriggerHeaderStruct *>(iter->fPtr);
      AliHLTCaloTriggerDataStruct *dataptr = reinterpret_cast<AliHLTCaloTriggerDataStruct *>(reinterpret_cast<AliHLTUInt8_t* >(iter->fPtr) + sizeof(AliHLTCaloTriggerHeaderStruct));
      HLTDebug("Block %d received %d fastor triggers", ndx, triggerhead->fNfastor);
      for(Int_t datacount = 0; datacount < triggerhead->fNfastor; datacount++) {
        fTriggerMakerPtr->SetADC(dataptr->fCol, dataptr->fRow, static_cast<Float_t>(dataptr->fL1TimeSum));
        fTriggerMakerPtr->SetL0Amplitude(dataptr->fCol, dataptr->fRow, static_cast<Float_t>(dataptr->fAmplitude));
        fTriggerMakerPtr->SetBitMask(dataptr->fCol, dataptr->fRow, dataptr->fTriggerBits);
        if (dataptr->fNL0Times > 0) {
          fTriggerMakerPtr->SetL0Time(Int_t(dataptr->fCol), Int_t(dataptr->fRow), dataptr->fL0Times[0]);
        }
        dataptr++;
      }
      nfastor += triggerhead->fNfastor;
    }
  }
  
  Int_t npatches = 0;
  if(nfastor || nDigitsGlob){
    fTriggerMakerPtr->SetTriggerPatchDataPtr(reinterpret_cast<AliHLTCaloTriggerPatchDataStruct *>(outputPtr), size - sizeof(AliHLTCaloTriggerHeaderStruct));
    npatches = fTriggerMakerPtr->FindPatches();
    HLTDebug("Found %d patches from fastors\n", npatches);
    outputPtr += npatches;
    mysize += sizeof(AliHLTCaloTriggerPatchDataStruct) * npatches;
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

#ifdef __PROFILE__
  profile.Stop();
  printf("End of trigger maker component: %f (Wall) / %f (CPU)\n", profile.RealTime(), profile.CpuTime());
#endif
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
  constBase = 10000 *(float)sizeof(AliHLTCaloTriggerPatchDataStruct) + sizeof(AliHLTCaloTriggerHeaderStruct);
  inputMultiplier = 0; // (float)sizeof(AliHLTCaloTriggerPatchDataStruct)/sizeof(AliHLTCaloDigitDataStruct)+1;
}

AliHLTComponent* AliHLTEMCALTriggerMakerComponent::Spawn(){
  return new AliHLTEMCALTriggerMakerComponent();
}


int AliHLTEMCALTriggerMakerComponent::DoInit ( int argc, const char** argv ){
  InitialiseGeometry();
  Float_t jetTh[2*AliHLTEMCALTriggerMaker::kNthresholds] = {0};
  Float_t gammaTh[2*AliHLTEMCALTriggerMaker::kNthresholds] = {0};
  Float_t bkgTh[2] = {0}, l0Th[2] = {0};
  Bool_t isPbPb = GetRunNo() > 244823 && GetRunNo() < 246995; // For the moment quick hack to distinguish PbPb from pp
  Bool_t runBkgAlgo = isPbPb;
  Int_t jetpatchsize = isPbPb ? 8 : 16;
  for(Int_t iarg = 0; iarg < argc; iarg++){
    TString argstring(argv[iarg]);
    if(argstring.Contains("-gammalowonthresh")){
      gammaTh[AliHLTEMCALTriggerMaker::kLowThreshold] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-gammalowoffthresh")){
      gammaTh[AliHLTEMCALTriggerMaker::kLowThreshold + AliHLTEMCALTriggerMaker::kNthresholds] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-gammahighonthresh")){
      gammaTh[AliHLTEMCALTriggerMaker::kHighThreshold] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-gammahighoffthresh")){
      gammaTh[AliHLTEMCALTriggerMaker::kHighThreshold + AliHLTEMCALTriggerMaker::kNthresholds] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-jetlowonthresh")){
      jetTh[AliHLTEMCALTriggerMaker::kLowThreshold] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-jetlowoffthresh")){
      jetTh[AliHLTEMCALTriggerMaker::kLowThreshold + AliHLTEMCALTriggerMaker::kNthresholds] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-jethighonthresh")){
      jetTh[AliHLTEMCALTriggerMaker::kHighThreshold] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-jethighoffthresh")){
      jetTh[AliHLTEMCALTriggerMaker::kHighThreshold + AliHLTEMCALTriggerMaker::kNthresholds] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-runbkgalgo")){
      runBkgAlgo = true;
    } else if(argstring.Contains("-bkgonthresh")){
      bkgTh[0] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-bkgoffthresh")){
      bkgTh[1] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-l0onthresh")){
      l0Th[0] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-l0offthresh")){
      l0Th[1] = TString(argv[++iarg]).Atof();
    } else if(argstring.Contains("-jetpatchsize")){
      jetpatchsize = TString(argv[++iarg]).Atoi();
    }
  }
  if(runBkgAlgo)
    HLTDebug("running jet background algorithm\n");
  else
    HLTDebug("do not run the jet background algorithm\n");
  fTriggerMakerPtr = new AliHLTEMCALTriggerMaker;
  fTriggerMakerPtr->SetJetPatch(jetpatchsize, 4);
  fTriggerMakerPtr->SetRunBkgAlgorithm(runBkgAlgo);
  for(UInt_t ithresh = AliHLTEMCALTriggerMaker::kHighThreshold; ithresh < AliHLTEMCALTriggerMaker::kNthresholds; ithresh++){
    fTriggerMakerPtr->SetGammaThresholds(AliHLTEMCALTriggerMaker::ThresholdType_t(ithresh), gammaTh[ithresh], gammaTh[ithresh + AliHLTEMCALTriggerMaker::kNthresholds]);
    fTriggerMakerPtr->SetJetThresholds(AliHLTEMCALTriggerMaker::ThresholdType_t(ithresh), jetTh[ithresh], jetTh[ithresh + AliHLTEMCALTriggerMaker::kNthresholds]);
  }
  fTriggerMakerPtr->SetLevel0Thresholds(l0Th[0], l0Th[1]);
  fTriggerMakerPtr->SetBkgThresholds(bkgTh[0], bkgTh[1]);
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
