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


#include "AliCaloRawAnalyzer.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliDAQ.h"
#include "AliHLTEMCALRawAnalyzerComponentTRU.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTCaloChannelDataStruct.h"
#include "AliHLTEMCALMapper.h"
#include "AliHLTCaloSanityInspector.h"
#include "AliHLTEMCALTRURawDigitMaker.h"
#include "AliRawReaderMemory.h"
#include "AliCaloRawStreamV3.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTCaloChannelRawDataStruct.h"
#include "AliLog.h"
#include "TStopwatch.h"

#include "AliCaloRawAnalyzerFactory.h"

//#include "AliCaloConstants.h"

//#include "AliCaloRawAnalyzer.h"

//using namespace Algo;

#include <vector>
using namespace std;

ClassImp(AliHLTEMCALRawAnalyzerComponentTRU);

AliHLTEMCALRawAnalyzerComponentTRU::AliHLTEMCALRawAnalyzerComponentTRU( ):AliHLTCaloProcessor(),
    AliHLTCaloConstantsHandler("EMCAL"),
    fTRUhandler(0),
    fDebug(false)
{
  //Constructor

}


AliHLTEMCALRawAnalyzerComponentTRU::~AliHLTEMCALRawAnalyzerComponentTRU()
{
  //destructor
}

void
AliHLTEMCALRawAnalyzerComponentTRU::GetInputDataTypes( vector <AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back( AliHLTEMCALDefinitions::fgkDDLRawDataType   | kAliHLTDataOriginEMCAL );
}


AliHLTComponentDataType
AliHLTEMCALRawAnalyzerComponentTRU::GetOutputDataType()
{
  //comment
  return AliHLTEMCALDefinitions::fgkTriggerRawDigitDataType;
}

AliHLTComponent* AliHLTEMCALRawAnalyzerComponentTRU::Spawn(){
  return new AliHLTEMCALRawAnalyzerComponentTRU;
}

const char*
AliHLTEMCALRawAnalyzerComponentTRU::GetComponentID()
{
  // component id
  return "EmcalTruAnalyzer";
}


int
AliHLTEMCALRawAnalyzerComponentTRU::DoInit( int argc, const char** argv )
{
  //See base class for documentation
  int iResult=0;

  for(int i = 0; i < argc; i++)
  {
    if(!strcmp("-suppressalilogwarnings", argv[i]))
    {
      AliLog::SetGlobalLogLevel(AliLog::kError);  //PHOS sometimes produces bad data -> Fill up the HLT logs...
    }
  }

  fTRUhandler = new AliHLTEMCALTRURawDigitMaker();
  fTRUhandler->Initialize(GetRunNo());

  return iResult;
}


int
AliHLTEMCALRawAnalyzerComponentTRU::DoDeinit()
{
  //comment
  if(fTRUhandler) delete fTRUhandler;

  return 0;
}

void
AliHLTEMCALRawAnalyzerComponentTRU::PrintDebugInfo()
{
  //comment
  static TStopwatch  watch; //CRAP PTH
  static double wlast = -1;
  static double wcurrent = 0;

  if( true == fDebug )
  {
    if( fCaloEventCount %1000 == 0  )
    {
      cout << __FILE__ << __LINE__ << " : Processing event "  << fCaloEventCount << endl;
      wlast =  wcurrent;
      wcurrent =  watch.RealTime();
      cout << __FILE__ << __LINE__ << "The event rate is " <<
          1000/( wcurrent  -  wlast ) << "  Hz" << endl;    watch.Start(kFALSE);
    }
  }
}


void
AliHLTEMCALRawAnalyzerComponentTRU::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  //comment
  constBase = 0;
  inputMultiplier = 1.5;
}


bool
AliHLTEMCALRawAnalyzerComponentTRU::CheckInputDataType(const AliHLTComponentDataType &datatype)
{
  //comment
  vector <AliHLTComponentDataType> validTypes;
  GetInputDataTypes(validTypes);

  for(UInt_t i=0; i < validTypes.size(); i++ )
  {
    if ( datatype  ==  validTypes.at(i) )
    {
      return true;
    }
  }

  HLTDebug("Invalid Datatype");
  return false;
}


int
AliHLTEMCALRawAnalyzerComponentTRU::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
    AliHLTComponentTriggerData& /*trigData*/,
    AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks )
{
  //comment
  if(!IsDataEvent())
  {
    size = 0;
    return 0;
  }
  if( true == fDebug )
  { PrintDebugInfo();
  };

  Int_t blockSize          = -1;
  UInt_t totSize           = 0;
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;
  AliHLTUInt32_t availableSize = size;

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
  {
    iter = blocks+ndx;
    if(  ! CheckInputDataType(iter->fDataType) )
    {
      continue;
    }

    if(iter->fSpecification >= AliDAQ::GetFirstSTUDDL()) continue;

    blockSize = DoIt(iter, outputPtr, availableSize, totSize); // Processing the block
    totSize += blockSize; //Keeping track of the used size
    availableSize -= blockSize;
    AliHLTComponentBlockData bdChannelData;
    FillBlockData( bdChannelData );
    bdChannelData.fOffset = totSize-blockSize; //FIXME
    bdChannelData.fSize = blockSize;
    bdChannelData.fDataType = GetOutputDataType();
    bdChannelData.fSpecification = iter->fSpecification;
    outputBlocks.push_back(bdChannelData);
    outputPtr += blockSize; //Updating position of the output buffer
  }

  fCaloEventCount++;
  size = totSize; //telling the framework how much buffer space we have used.

  return 0;

}//end DoEvent



Int_t
AliHLTEMCALRawAnalyzerComponentTRU::DoIt(const AliHLTComponentBlockData* iter, AliHLTUInt8_t* outputPtr, const AliHLTUInt32_t size, UInt_t& totSize)
{
  //comment
  int tmpsize=  0;
  Int_t crazyness          = 0;
  Int_t nSamples           = 0;
  Short_t channelCount     = 0;

  AliHLTCaloTriggerRawDigitDataStruct *digitDataPtr = reinterpret_cast<AliHLTCaloTriggerRawDigitDataStruct*>(outputPtr);
  AliRawReaderMemory rawReaderMemoryPtr;
  rawReaderMemoryPtr.SetMemory(         reinterpret_cast<UChar_t*>( iter->fPtr ),  static_cast<ULong_t>( iter->fSize )  );
  rawReaderMemoryPtr.SetEquipmentID(    iter->fSpecification + fCaloConstants->GetDDLOFFSET() );

  AliCaloRawStreamV3 altroRawStreamPtr(&rawReaderMemoryPtr, "EMCAL");
  rawReaderMemoryPtr.Reset();
  rawReaderMemoryPtr.NextEvent();


  fTRUhandler->SetRawReader(&altroRawStreamPtr);
  fTRUhandler->Reset();

  if(altroRawStreamPtr.NextDDL()) {
    int cnt = 0;
    while(altroRawStreamPtr.NextChannel()) {
      if( altroRawStreamPtr.GetCaloFlag() == 2) { // These are the level1 triggers
        std::vector <AliCaloBunchInfo> bunchlist;
        while (altroRawStreamPtr.NextBunch())
          bunchlist.push_back( AliCaloBunchInfo(altroRawStreamPtr.GetStartTimeBin(), altroRawStreamPtr.GetBunchLength(), altroRawStreamPtr.GetSignals() ) );
        if (bunchlist.size() == 0) continue;

        fTRUhandler->Add(bunchlist);
      }
    }
  }

  AliHLTUInt32_t availableSize = size;
  HLTDebug("Number of TRU digits: %d",fTRUhandler->GetNumberOfRawDigits());
  tmpsize = fTRUhandler->WriteRawDigitsBuffer(digitDataPtr, availableSize);

  /*
  std::cout << "Found TRU  raw digits: " << std::endl;
  for(Int_t idig = 0; idig < fTRUhandler->GetNumberOfRawDigits(); idig++){
    PrintRawDigit(digitDataPtr[idig]);
  }
  */
  return  tmpsize;
}



