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
#include "AliHLTCaloTriggerDataStruct.h"
#include "AliHLTCaloTriggerHeaderStruct.h"
#include "AliHLTCaloTriggerRawDigitDataStruct.h"
#include "AliEMCALTriggerSTURawStream.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALTriggerDataMakerComponent.h"

ClassImp(AliHLTEMCALTriggerDataMakerComponent)

AliHLTEMCALTriggerDataMakerComponent::AliHLTEMCALTriggerDataMakerComponent():
AliHLTCaloProcessor(),
AliHLTCaloConstantsHandler("EMCAL"),
fGeometry(NULL),
fSTUHeader(),
fNRawDigitsTRU(0),
fNRawDigitsSTU(0),
fMaxChannel(0)
{
  for(Short_t iter = 0; iter < kMaxChannels; iter++){
    fRawIndexesTRU[iter] = -1;
    fRawIndexesSTU[iter] = -1;
  }
}

AliHLTEMCALTriggerDataMakerComponent::~AliHLTEMCALTriggerDataMakerComponent() {
  if(fGeometry) delete fGeometry;
}

int AliHLTEMCALTriggerDataMakerComponent::DoInit(int argc, const char **argv){
  fGeometry = new AliHLTEMCALGeometry(GetRunNo());
  return 0;
}

int AliHLTEMCALTriggerDataMakerComponent::DoDeinit(){
  if(fGeometry) delete fGeometry;
  fGeometry = NULL;
  return 0;
}

const char* AliHLTEMCALTriggerDataMakerComponent::GetComponentID(){
  return "EmcalTriggerDataMaker";
}

void AliHLTEMCALTriggerDataMakerComponent::GetInputDataTypes( std::vector <AliHLTComponentDataType>& list){
  list.clear();
  list.push_back( AliHLTEMCALDefinitions::fgkTriggerRawDigitDataType   | kAliHLTDataOriginEMCAL );
  list.push_back( AliHLTEMCALDefinitions::fgkTriggerSTUDataType | kAliHLTDataOriginEMCAL );
}

AliHLTComponentDataType AliHLTEMCALTriggerDataMakerComponent::GetOutputDataType(){
  return kAliHLTDataTypeCaloTrigger | kAliHLTDataOriginEMCAL;
}

void AliHLTEMCALTriggerDataMakerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier){
  constBase = sizeof(AliHLTCaloTriggerHeaderStruct);
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

  UInt_t totSize            = 0;
  const AliHLTComponentBlockData* iter = NULL;

  // Get pointers to output buffer
  AliHLTCaloTriggerHeaderStruct *headerPtr = reinterpret_cast<AliHLTCaloTriggerHeaderStruct *>(outputPtr);
  AliHLTCaloTriggerDataStruct *dataIter = reinterpret_cast<AliHLTCaloTriggerDataStruct *>(outputPtr + sizeof(AliHLTCaloTriggerHeaderStruct));
  totSize += sizeof(AliHLTCaloTriggerHeaderStruct);

  Reset();
  AliHLTCaloTriggerRawDigitDataStruct *dataptr = NULL;
  for(ULong_t ndx = 0; ndx < evtData.fBlockCnt; ndx++){
    iter = blocks + ndx;

    if(!this->CheckInputDataType(iter->fDataType)){
      continue;
    }

    if(iter->fDataType == AliHLTEMCALDefinitions::fgkTriggerRawDigitDataType){
      // Handle TRU data
      Int_t ndigits = iter->fSize / sizeof(AliHLTCaloTriggerRawDigitDataStruct);
      HLTDebug("Data containing %d TRU digits", ndigits);
      dataptr = reinterpret_cast<AliHLTCaloTriggerRawDigitDataStruct *>(iter->fPtr);
      ReadTRUData(ndigits, dataptr);
    } else if(iter->fDataType == AliHLTEMCALDefinitions::fgkTriggerSTUDataType){
      // Handle STU data
      AliHLTEMCALSTUHeaderStruct *stuheader = reinterpret_cast<AliHLTEMCALSTUHeaderStruct *>(iter->fPtr);
      AliHLTInt32_t sizeExpected = sizeof(AliHLTEMCALSTUHeaderStruct) + sizeof(AliHLTCaloTriggerRawDigitDataStruct) * stuheader->fNRawDigits;
      if(iter->fSize != sizeExpected){
        HLTWarning("STU-reader: Size of the input buffer not matching for amount of digits, expected %d, obtained %d", sizeExpected, iter->fSize);
        continue;
      }
      dataptr = reinterpret_cast<AliHLTCaloTriggerRawDigitDataStruct *>(reinterpret_cast<AliHLTUInt8_t*>(iter->fPtr) + sizeof(AliHLTEMCALSTUHeaderStruct));
      HLTDebug("Data containing %d STU digits", stuheader->fNRawDigits);
      ReadSTUData(stuheader, dataptr);
    }
  }

  // Write header
  memcpy(headerPtr->fL1Threshold, fSTUHeader.fL1Threshold, sizeof(Int_t) *4);
  memcpy(headerPtr->fL1V0, fSTUHeader.fL1V0, sizeof(Int_t) * 2);
  headerPtr->fL1FrameMask = fSTUHeader.fL1FrameMask;

  AliHLTUInt32_t availableSize = size - sizeof(AliHLTCaloTriggerHeaderStruct);

  // Write data
  Int_t dataSize = MakeTriggerData(dataIter, availableSize);
  totSize += dataSize;
  headerPtr->fNfastor = dataSize / sizeof(AliHLTCaloTriggerDataStruct);

  AliHLTComponentBlockData bdChannelData;
  FillBlockData( bdChannelData );
  bdChannelData.fOffset = 0; //FIXME
  bdChannelData.fSize = totSize;
  bdChannelData.fDataType = GetOutputDataType();
  outputBlocks.push_back(bdChannelData);
  outputPtr += totSize; //Updating position of the output buffer
  size = totSize;
  return 0;
}


void AliHLTEMCALTriggerDataMakerComponent::ReadSTUData(AliHLTEMCALSTUHeaderStruct *headerptr, AliHLTCaloTriggerRawDigitDataStruct *dataptr){
  fSTUHeader = *headerptr;
  for(UShort_t idig = 0; idig < headerptr->fNRawDigits; idig++){
	if(dataptr->fID > kMaxChannels || dataptr->fID < 0){
		HLTWarning("Invalid TRU index: %d", dataptr->fID);
		dataptr++;
		continue;
	}
    fRawIndexesSTU[dataptr->fID] = fNRawDigitsSTU;
    if(dataptr->fID > fMaxChannel) fMaxChannel = dataptr->fID;
    fSTURawDigitBuffer[fNRawDigitsSTU] = *dataptr;
    dataptr++;
    fNRawDigitsSTU++;
  }
  HLTDebug("Successfully read in %d STU digits", fNRawDigitsSTU);
}

void AliHLTEMCALTriggerDataMakerComponent::ReadTRUData(UShort_t ndigits, AliHLTCaloTriggerRawDigitDataStruct *triggerdata){
  for(UShort_t idig = 0; idig < ndigits; idig++){
    fRawIndexesTRU[triggerdata->fID] = fNRawDigitsTRU;
    if(triggerdata->fID > fMaxChannel) fMaxChannel = triggerdata->fID;
    fTRURawDigitBuffer[fNRawDigitsTRU] = *triggerdata;
    triggerdata++;
    fNRawDigitsTRU++;
  }
  HLTDebug("Successfully read in %d TRU digits", fNRawDigitsTRU);
}

Int_t AliHLTEMCALTriggerDataMakerComponent::MakeTriggerData(AliHLTCaloTriggerDataStruct *outputdata, AliHLTUInt32_t &availableSize) {
  if(availableSize < sizeof(AliHLTCaloTriggerDataStruct)){
	  HLTWarning("Not enough space in buffer to write triggers");
	  return 0;
  }
  Int_t outputsize = 0, col = 0, row = 0, ntriggers = 0 ;
  AliHLTCaloTriggerRawDigitDataStruct tmpdigit;
  for(UShort_t indcounter = 0; indcounter <= fMaxChannel; indcounter++){
    if(availableSize < sizeof(AliHLTCaloTriggerDataStruct)){
      HLTWarning("Buffer exceeded after %d triggers", ntriggers);
      break;
    }
    fGeometry->GetGeometryPtr()->GetPositionInEMCALFromAbsFastORIndex(indcounter, col, row);
    if(fRawIndexesTRU[indcounter] >= 0 && fRawIndexesSTU[indcounter] >=0){
      CombineTRUSTUDigit(tmpdigit, fTRURawDigitBuffer[fRawIndexesTRU[indcounter]], fSTURawDigitBuffer[fRawIndexesSTU[indcounter]]);
      ConvertRawDigit(outputdata, &tmpdigit, col, row);
      outputsize += sizeof(AliHLTCaloTriggerDataStruct);
      availableSize -= sizeof(AliHLTCaloTriggerDataStruct);
      outputdata++;
      ntriggers++;
    } else if(fRawIndexesTRU[indcounter] >= 0){
      ConvertRawDigit(outputdata, &(fTRURawDigitBuffer[fRawIndexesTRU[indcounter]]), col, row);
      outputsize += sizeof(AliHLTCaloTriggerDataStruct);
      availableSize -= sizeof(AliHLTCaloTriggerDataStruct);
      outputdata++;
      ntriggers++;
    } else if(fRawIndexesSTU[indcounter] >= 0){
      ConvertRawDigit(outputdata, &(fSTURawDigitBuffer[fRawIndexesSTU[indcounter]]), col, row);
      outputsize += sizeof(AliHLTCaloTriggerDataStruct);
      availableSize -= sizeof(AliHLTCaloTriggerDataStruct);
      outputdata++;
      ntriggers++;
    }
  }
  return outputsize;
}

void AliHLTEMCALTriggerDataMakerComponent::CombineTRUSTUDigit(
    AliHLTCaloTriggerRawDigitDataStruct &target,
    const AliHLTCaloTriggerRawDigitDataStruct &trudigit,
    const AliHLTCaloTriggerRawDigitDataStruct &studigit){
  AliHLTCaloTriggerRawDigitDataStruct merged;
  target.fID = trudigit.fID;
  target.fNTimeSamples = trudigit.fNTimeSamples;
  memcpy(target.fTimeSamples, trudigit.fTimeSamples, sizeof(target.fTimeSamples));
  target.fNL0Times = trudigit.fNL0Times;
  memcpy(target.fL0Times, trudigit.fL0Times, sizeof(target.fL0Times));
  target.fL1TimeSum = studigit.fL1TimeSum;
  target.fTriggerBits = trudigit.fTriggerBits | studigit.fTriggerBits;
}

void AliHLTEMCALTriggerDataMakerComponent::Reset(){
  for(Short_t iter = 0; iter < kMaxChannels; iter++){
    fRawIndexesTRU[iter] = -1;
    fRawIndexesSTU[iter] = -1;
  }
  fNRawDigitsTRU = 0;
  fNRawDigitsSTU = 0;
  fMaxChannel = 0;
}

void AliHLTEMCALTriggerDataMakerComponent::ConvertRawDigit(AliHLTCaloTriggerDataStruct *target, const AliHLTCaloTriggerRawDigitDataStruct *source, Int_t col, Int_t row) {
  target->fCol = col;
  target->fRow = row;
  Int_t amplitude(0), time(0);
  GetRawDigitMaximumAmplitude(*source, amplitude, time);
  target->fAmplitude = static_cast<Float_t>(amplitude);
  target->fTime = static_cast<Float_t>(time);
  target->fL1TimeSum = source->fL1TimeSum;
  target->fNL0Times = source->fNL0Times;
  memcpy(target->fL0Times, source->fL0Times, sizeof(UChar_t) * 10);
  target->fTriggerBits = source->fTriggerBits;
}
