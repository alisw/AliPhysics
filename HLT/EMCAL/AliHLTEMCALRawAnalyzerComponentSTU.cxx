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
#include "AliHLTCaloTriggerRawDigitDataStruct.h"
#include "AliCaloRawStreamV3.h"
#include "AliEMCALTriggerData.h"
#include "AliEMCALTriggerSTURawStream.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALMapper.h"
#include "AliHLTEMCALRawAnalyzerComponentSTU.h"
#include "AliHLTEMCALSTURawDigitMaker.h"
#include "AliHLTEMCALSTUHeaderStruct.h"
#include "AliRawReaderMemory.h"
#include "AliDAQ.h"

ClassImp(AliHLTEMCALRawAnalyzerComponentSTU)

AliHLTEMCALRawAnalyzerComponentSTU::AliHLTEMCALRawAnalyzerComponentSTU():
AliHLTCaloProcessor(),
AliHLTCaloConstantsHandler("EMCAL"),
fRawReaderMemoryPtr(NULL),
fSTURawDigitMaker(NULL),
fGeometry(NULL)
{
}

AliHLTEMCALRawAnalyzerComponentSTU::~AliHLTEMCALRawAnalyzerComponentSTU() {
  if(fGeometry) delete fGeometry;
  if(fRawReaderMemoryPtr) delete fRawReaderMemoryPtr;
}

int AliHLTEMCALRawAnalyzerComponentSTU::DoInit(int argc, const char **argv){
  fGeometry = new AliHLTEMCALGeometry(GetRunNo());

  fSTURawDigitMaker = new AliHLTEMCALSTURawDigitMaker;
  fSTURawDigitMaker->SetGeometry(fGeometry);

  fRawReaderMemoryPtr = new AliRawReaderMemory();

  return 0;
}

int AliHLTEMCALRawAnalyzerComponentSTU::DoDeinit(){
  if(fGeometry) delete fGeometry;
  fGeometry = NULL;
  if(fRawReaderMemoryPtr) delete fRawReaderMemoryPtr;
  fRawReaderMemoryPtr = NULL;
  return 0;
}

const char* AliHLTEMCALRawAnalyzerComponentSTU::GetComponentID(){
  return "EmcalStuAnalyzer";
}

void AliHLTEMCALRawAnalyzerComponentSTU::GetInputDataTypes( std::vector <AliHLTComponentDataType>& list){
  list.clear();
  list.push_back( AliHLTEMCALDefinitions::fgkDDLRawDataType   | kAliHLTDataOriginEMCAL );
}

AliHLTComponentDataType AliHLTEMCALRawAnalyzerComponentSTU::GetOutputDataType(){
  return AliHLTEMCALDefinitions::fgkTriggerSTUDataType;
}

void AliHLTEMCALRawAnalyzerComponentSTU::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier){
  // Allocate size for 1000 digits
  constBase = sizeof(AliHLTEMCALSTUHeaderStruct) + 1000 * sizeof(AliHLTCaloTriggerRawDigitDataStruct);
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTEMCALRawAnalyzerComponentSTU::Spawn(){
  return new AliHLTEMCALRawAnalyzerComponentSTU;
}

bool AliHLTEMCALRawAnalyzerComponentSTU::CheckInputDataType(const AliHLTComponentDataType &datatype) {
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

int AliHLTEMCALRawAnalyzerComponentSTU::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
         AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr,
         AliHLTUInt32_t& size, std::vector<AliHLTComponentBlockData>& outputBlocks ){

  if(!IsDataEvent()) {
    size = 0;
    return 0;
  }

  UInt_t totSize            = 0;
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;

  // Get pointers to output buffer
  AliHLTEMCALSTUHeaderStruct *headerPtr = reinterpret_cast<AliHLTEMCALSTUHeaderStruct *>(outputPtr);
  AliHLTCaloTriggerRawDigitDataStruct *dataIter = reinterpret_cast<AliHLTCaloTriggerRawDigitDataStruct *>(outputPtr + sizeof(AliHLTEMCALSTUHeaderStruct)),
      *nextFastor(NULL);
  totSize += sizeof(AliHLTEMCALSTUHeaderStruct);

  bool headerInitialized = false;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ ) {
    iter = blocks+ndx;
    if(  ! CheckInputDataType(iter->fDataType) ) {
      continue;
    }
    if(iter->fSpecification < AliDAQ::GetFirstSTUDDL() || iter->fSpecification > AliDAQ::GetLastSTUDDL()) // check for STU DDLs
      continue;

    // Initialize raw reader from input data
    fRawReaderMemoryPtr->SetMemory(reinterpret_cast<UChar_t*>( iter->fPtr ), static_cast<ULong_t>(iter->fSize));
    fRawReaderMemoryPtr->SetEquipmentID(iter->fSpecification + fCaloConstants->GetDDLOFFSET());

    fSTURawDigitMaker->Reset();

    AliEMCALTriggerSTURawStream stustream(fRawReaderMemoryPtr);
    fSTURawDigitMaker->ProcessSTUStream(fRawReaderMemoryPtr, &stustream);
    const AliEMCALTriggerData *triggerData = fSTURawDigitMaker->GetTriggerData();

    if(!headerInitialized){
      // Set Header
      for (int i = 0; i < 2; i++) {
        headerPtr->fL1Threshold[2*i] = triggerData->GetL1JetThreshold(i);
        headerPtr->fL1Threshold[2*i+1] = triggerData->GetL1GammaThreshold(i);
      }
      headerPtr->fL1FrameMask = triggerData->GetL1FrameMask();
      headerInitialized = true;
    }
    totSize += fSTURawDigitMaker->WriteRawDigitsBuffer(dataIter);
    dataIter += fSTURawDigitMaker->GetNumberOfRawDigits();
  }

  headerPtr->fNRawDigits = fSTURawDigitMaker->GetNumberOfRawDigits();
  HLTDebug("Successfully decoded %d digits.", headerPtr->fNRawDigits);

  AliHLTComponentBlockData bdChannelData;
  FillBlockData( bdChannelData );
  bdChannelData.fOffset = 0; //FIXME
  bdChannelData.fSize = totSize;
  bdChannelData.fDataType = GetOutputDataType();
  bdChannelData.fSpecification = 0;
  outputBlocks.push_back(bdChannelData);
  outputPtr += totSize; //Updating position of the output buffer
  size = totSize;
  return 0;
}
