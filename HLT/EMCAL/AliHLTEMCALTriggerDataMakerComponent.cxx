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
#include "AliCaloRawStreamV3.h"
#include "AliEMCALTriggerData.h"
#include "AliEMCALTriggerSTURawStream.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALMapper.h"
#include "AliHLTEMCALTriggerDataMakerComponent.h"
#include "AliHLTEMCALTriggerHeaderStruct.h"
#include "AliHLTEMCALTriggerFastorDataStruct.h"
#include "AliHLTEMCALTriggerRawDigitDataStruct.h"
#include "AliHLTEMCALTriggerRawDigitMaker.h"
#include "AliRawReaderMemory.h"

ClassImp(AliHLTEMCALTriggerDataMakerComponent)

AliHLTEMCALTriggerDataMakerComponent::AliHLTEMCALTriggerDataMakerComponent():
AliHLTCaloProcessor(),
AliHLTCaloConstantsHandler("EMCAL"),
fMapperPtr(NULL),
fCurrentSpec(0),
fRawReaderMemoryPtr(NULL),
fCaloRawStreamPtr(NULL),
fRawDigitMaker(NULL),
fGeometry(NULL),
fTriggerData(NULL)
{
}

AliHLTEMCALTriggerDataMakerComponent::~AliHLTEMCALTriggerDataMakerComponent() {
  if(fMapperPtr) delete fMapperPtr;
  if(fGeometry) delete fGeometry;
  if(fRawDigitMaker) delete fRawDigitMaker;
  if(fRawReaderMemoryPtr) delete fRawReaderMemoryPtr;
  if(fCaloRawStreamPtr) delete fCaloRawStreamPtr;
  if(fTriggerData) delete fTriggerData;
}

int AliHLTEMCALTriggerDataMakerComponent::DoInit(int argc, const char **argv){
  fGeometry = new AliHLTEMCALGeometry(GetRunNo());

  fRawDigitMaker = new AliHLTEMCALTriggerRawDigitMaker;
  fRawDigitMaker->SetGeometry(fGeometry);

  fRawReaderMemoryPtr = new AliRawReaderMemory();
  fCaloRawStreamPtr = new AliCaloRawStreamV3(fRawReaderMemoryPtr, "EMCAL");

  fTriggerData = new AliEMCALTriggerData;

  return 0;
}

int AliHLTEMCALTriggerDataMakerComponent::DoDeinit(){
  if(fMapperPtr) delete fMapperPtr;
  if(fGeometry) delete fGeometry;
  if(fRawDigitMaker) delete fRawDigitMaker;
  if(fRawReaderMemoryPtr) delete fRawReaderMemoryPtr;
  if(fCaloRawStreamPtr) delete fCaloRawStreamPtr;
  if(fTriggerData) delete fTriggerData;
  return 0;
}

const char* AliHLTEMCALTriggerDataMakerComponent::GetComponentID(){
  return "EmcalTriggerRawDigitMaker";
}

void AliHLTEMCALTriggerDataMakerComponent::GetInputDataTypes( std::vector <AliHLTComponentDataType>& list){
  list.clear();
  list.push_back( AliHLTEMCALDefinitions::fgkDDLRawDataType   | kAliHLTDataOriginEMCAL );
}

AliHLTComponentDataType AliHLTEMCALTriggerDataMakerComponent::GetOutputDataType(){
  return AliHLTEMCALDefinitions::fgkTriggerRawDigitDataType;
}

void AliHLTEMCALTriggerDataMakerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier){
  constBase = 0;
  inputMultiplier = 1.5;
}

AliHLTComponent* AliHLTEMCALTriggerDataMakerComponent::Spawn(){
  return new AliHLTEMCALTriggerDataMakerComponent;
}

bool AliHLTEMCALTriggerDataMakerComponent::CheckInputDataType(const AliHLTComponentDataType &datatype) {
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

int AliHLTEMCALTriggerDataMakerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
         AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr,
         AliHLTUInt32_t& size, std::vector<AliHLTComponentBlockData>& outputBlocks ){

  if(!IsDataEvent()) {
    size = 0;
    return 0;
  }

  UInt_t specification      = 0;
  UInt_t totSize            = 0;
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;

  // Get pointers to output buffer
  AliHLTEMCALTriggerHeaderStruct *headerPtr = reinterpret_cast<AliHLTEMCALTriggerHeaderStruct *>(outputPtr);
  AliHLTEMCALTriggerFastorDataStruct *dataIter = reinterpret_cast<AliHLTEMCALTriggerFastorDataStruct *>(outputPtr + sizeof(AliHLTEMCALTriggerHeaderStruct)),
      *nextFastor(NULL);
  totSize += sizeof(AliHLTEMCALTriggerHeaderStruct);

  bool headerInitialized = false;
  Int_t nfastor = 0;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ ) {
    iter = blocks+ndx;
    if(  ! CheckInputDataType(iter->fDataType) ) {
      continue;
    }

    if(iter->fSpecification != fCurrentSpec) {
      fCurrentSpec = iter->fSpecification;
      InitMapping(iter->fSpecification);
    }
    specification |= iter->fSpecification;

    // Initialize raw reader from input data
    fRawReaderMemoryPtr->SetMemory(reinterpret_cast<UChar_t*>( iter->fPtr ), static_cast<ULong_t>(iter->fSize));
    fRawReaderMemoryPtr->SetEquipmentID(fMapperPtr->GetDDLFromSpec(iter->fSpecification) + fCaloConstants->GetDDLOFFSET());
    fRawReaderMemoryPtr->Reset();
    fRawReaderMemoryPtr->NextEvent();

    fTriggerData->Reset();
    fRawDigitMaker->Reset();

    AliEMCALTriggerSTURawStream stustream(fRawReaderMemoryPtr);
    fRawDigitMaker->SetIO(fRawReaderMemoryPtr, *fCaloRawStreamPtr, stustream, fTriggerData);

    Int_t caloFlag(0);
    while (fCaloRawStreamPtr->NextDDL()) {
      while (fCaloRawStreamPtr->NextChannel()) {
        caloFlag = fCaloRawStreamPtr->GetCaloFlag();

        if ( caloFlag != 2 ) continue; // Only FALTRO
        vector<AliCaloBunchInfo> bunchlist;

        while (fCaloRawStreamPtr->NextBunch())
          bunchlist.push_back( AliCaloBunchInfo(fCaloRawStreamPtr->GetStartTimeBin(), fCaloRawStreamPtr->GetBunchLength(), fCaloRawStreamPtr->GetSignals() ) );

        if (bunchlist.size() == 0) continue;

        fRawDigitMaker->Add(bunchlist);
      } // End while over channel
    } // End while over DDL's, of input stream
    fRawDigitMaker->PostProcess();

    if(!headerInitialized){
      // Set Header
      for (int i = 0; i < 2; i++) {
        headerPtr->fL1Threshold[2*i] = fTriggerData->GetL1JetThreshold(i);
        headerPtr->fL1Threshold[2*i+1] = fTriggerData->GetL1GammaThreshold(i);
      }
      headerPtr->fL1FrameMask = fTriggerData->GetL1FrameMask();
      headerInitialized = true;
    }

    // Write out stuff
    const std::vector<AliHLTEMCALTriggerRawDigitDataStruct> &indigits = fRawDigitMaker->GetRawDigits();
    for(std::vector<AliHLTEMCALTriggerRawDigitDataStruct>::const_iterator digiter = indigits.begin(); digiter != indigits.end(); ++digiter){
      Int_t col, row;
      if(fGeometry->GetGeometryPtr()->GetPositionInEMCALFromAbsFastORIndex(digiter->fID, col, row)){
        nextFastor = dataIter + 1;
        ConvertRawDigit(dataIter, &(*digiter), col, row);
        dataIter = nextFastor;
        nfastor++;
        totSize += sizeof(AliHLTEMCALTriggerFastorDataStruct);
      }
    }
  }
  headerPtr->fNfastor = nfastor;

  AliHLTComponentBlockData bdChannelData;
  FillBlockData( bdChannelData );
  bdChannelData.fOffset = 0; //FIXME
  bdChannelData.fSize = totSize;
  bdChannelData.fDataType = GetOutputDataType();
  bdChannelData.fSpecification = specification;
  outputBlocks.push_back(bdChannelData);
  outputPtr += totSize; //Updating position of the output buffer
  size = totSize;
  return 0;
}

void AliHLTEMCALTriggerDataMakerComponent::ConvertRawDigit(AliHLTEMCALTriggerFastorDataStruct *target, const AliHLTEMCALTriggerRawDigitDataStruct *source, Int_t col, Int_t row) {
  target->fCol = col;
  target->fRow = row;
  Int_t amplitude, time;
  GetRawDigitMaximumAmplitude(*source, amplitude, time);
  target->fAmplitude = static_cast<Float_t>(amplitude);
  target->fTime = static_cast<Float_t>(time);
  target->fL1TimeSum = source->fL1TimeSum;
  target->fNL0Times = source->fNL0Times;
  memcpy(target->fL0Times, source->fL0Times, sizeof(UChar_t) * 10);
  target->fTriggerBits = source->fTriggerBits;
}

void AliHLTEMCALTriggerDataMakerComponent::InitMapping( const int specification ) {
  if (!fMapperPtr)  fMapperPtr =  new AliHLTEMCALMapper( specification );

  if(fMapperPtr->GetIsInitializedMapping() == false ) {
    HLTError("%d:%d, ERROR, mapping not initialized ", __FILE__, __LINE__ );
    exit(-2);
  }
}
