 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliTPCCalibRawBase.cxx */

/////////////////////////////////////////////////////////////////////////////////////////
//                                                                                     //
//          Base class for the calibration algorithms using raw data as input          //
//                                                                                     //
//   Origin: Jens Wiechula   J.Wiechula@gsi.de                                         //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

//Root includes
#include <TDirectory.h>
#include <TFile.h>

//Aliroot includes
#include "AliRawReaderDate.h"
#include "AliRawReader.h"
#include "AliRawEventHeaderBase.h"
#include "AliAltroMapping.h"
#include "AliAltroRawStream.h"
#include "AliTPCROC.h"
#include "AliTPCRawStreamFast.h"
#include "AliTPCRawStream.h"
#include "TTreeStream.h"
#include "event.h"

#include "AliTPCCalibRawBase.h"

ClassImp(AliTPCCalibRawBase)

AliTPCCalibRawBase::AliTPCCalibRawBase() :
  TNamed(),
  fFirstTimeBin(0),
  fLastTimeBin(1000),
  fNevents(0),
  fDebugLevel(0),
  fStreamLevel(0),
  fTimeStamp(0),
  fRunNumber(0),
  fEventType(0),
  fAltroL1Phase(0),
  fAltroL1PhaseTB(0),
  fUseL1Phase(kTRUE),
  fDebugStreamer(0x0),
  fAltroRawStream(0x0),
  fMapping(0x0),
  fROC(AliTPCROC::Instance())
{
    //
    // default ctor
    //

}
//_____________________________________________________________________
AliTPCCalibRawBase::AliTPCCalibRawBase(const AliTPCCalibRawBase &calib) :
  TNamed(calib),
  fFirstTimeBin(calib.fFirstTimeBin),
  fLastTimeBin(calib.fLastTimeBin),
  fNevents(calib.fNevents),
  fDebugLevel(calib.fDebugLevel),
  fStreamLevel(calib.fStreamLevel),
  fTimeStamp(0),
  fRunNumber(0),
  fEventType(0),
  fAltroL1Phase(0),
  fAltroL1PhaseTB(0),
  fUseL1Phase(kTRUE),
  fDebugStreamer(0x0),
  fAltroRawStream(0x0),
  fMapping(0x0),
  fROC(AliTPCROC::Instance())
{
    //
    // copy ctor
    //
  
}
//_____________________________________________________________________
AliTPCCalibRawBase::~AliTPCCalibRawBase()
{
  //
  // dtor
  //
  if (fDebugStreamer) delete fDebugStreamer;
}
//_____________________________________________________________________
  AliTPCCalibRawBase& AliTPCCalibRawBase::operator = (const  AliTPCCalibRawBase &source)
  {
    //
    // assignment operator
    //
    if (&source == this) return *this;
    new (this) AliTPCCalibRawBase(source);
    
    return *this;
  }
//_____________________________________________________________________
Bool_t AliTPCCalibRawBase::ProcessEventFast(AliTPCRawStreamFast *rawStreamFast)
{
  //
  // Event Processing loop - AliTPCRawStreamFast
  //
  ResetEvent();
  Bool_t withInput = kFALSE;
  while ( rawStreamFast->NextDDL() ){
    while ( rawStreamFast->NextChannel() ){
      Int_t isector  = rawStreamFast->GetSector();                       //  current sector
      Int_t iRow     = rawStreamFast->GetRow();                          //  current row
      Int_t iPad     = rawStreamFast->GetPad();                          //  current pad
      
      while ( rawStreamFast->NextBunch() ){
        Int_t startTbin = (Int_t)rawStreamFast->GetStartTimeBin();
        Int_t endTbin = (Int_t)rawStreamFast->GetEndTimeBin();
        for (Int_t iTimeBin = startTbin; iTimeBin < endTbin; iTimeBin++){
          Float_t signal=(Float_t)rawStreamFast->GetSignals()[iTimeBin-startTbin];
          Update(isector,iRow,iPad,iTimeBin+1,signal);
          withInput = kTRUE;
        }
      }
    }
  }
  if (withInput){
    EndEvent();
  }
  return withInput;
}
//_____________________________________________________________________
Bool_t AliTPCCalibRawBase::ProcessEventFast(AliRawReader *rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //
  AliRawEventHeaderBase* eventHeader = (AliRawEventHeaderBase*)rawReader->GetEventHeader();
  if (eventHeader){
    fTimeStamp   = eventHeader->Get("Timestamp");
    fRunNumber = eventHeader->Get("RunNb");
    fEventType = eventHeader->Get("Type");
  }
  AliTPCRawStreamFast *rawStreamFast = new AliTPCRawStreamFast(rawReader, (AliAltroMapping**)fMapping);
  Bool_t res=ProcessEventFast(rawStreamFast);
  delete rawStreamFast;
  return res;
}
//_____________________________________________________________________
Bool_t AliTPCCalibRawBase::ProcessEvent(AliTPCRawStream *rawStream)
{
  //
  // Event Processing loop - AliTPCRawStream
  //

  ResetEvent();

  Bool_t withInput = kFALSE;
  fAltroL1Phase=0;
  fAltroL1PhaseTB=0;
  fAltroRawStream = static_cast<AliAltroRawStream*>(rawStream);
  while (rawStream->Next()) {
    if (fUseL1Phase){
      fAltroL1Phase  = fAltroRawStream->GetL1Phase();
      fAltroL1PhaseTB = (fAltroL1Phase*1e09/100.);
    }
    Int_t isector  = rawStream->GetSector();                       //  current sector
    Int_t iRow     = rawStream->GetRow();                          //  current row
    Int_t iPad     = rawStream->GetPad();                          //  current pad
    Int_t iTimeBin = rawStream->GetTime();                         //  current time bin
    Float_t signal = rawStream->GetSignal();                       //  current ADC signal
    
    Update(isector,iRow,iPad,iTimeBin,signal);
    withInput = kTRUE;
  }
  fAltroRawStream=0x0;
  if (withInput){
    EndEvent();
  }
  return withInput;
}
//_____________________________________________________________________
Bool_t AliTPCCalibRawBase::ProcessEvent(AliRawReader *rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //
  AliRawEventHeaderBase* eventHeader = (AliRawEventHeaderBase*)rawReader->GetEventHeader();
  if (eventHeader){
    fTimeStamp   = eventHeader->Get("Timestamp");
    fRunNumber = eventHeader->Get("RunNb");
    fEventType = eventHeader->Get("Type");
  }

  AliTPCRawStream rawStream(rawReader, (AliAltroMapping**)fMapping);
  rawReader->Select("TPC");
  return ProcessEvent(&rawStream);
}
//_____________________________________________________________________
Bool_t AliTPCCalibRawBase::ProcessEvent(eventHeaderStruct *event)
{
  //
  //  Event processing loop - date event
  //
    AliRawReader *rawReader = new AliRawReaderDate((void*)event);
    Bool_t result=ProcessEvent(rawReader);
    delete rawReader;
    return result;

}
//_____________________________________________________________________
void AliTPCCalibRawBase::DumpToFile(const Char_t *filename, const Char_t *dir, Bool_t append)
{
    //
    //  Write class to file
    //
  
  TString sDir(dir);
  TString option;
  
  if ( append )
    option = "update";
  else
    option = "recreate";
  
  TDirectory *backup = gDirectory;
  TFile f(filename,option.Data());
  f.cd();
  if ( !sDir.IsNull() ){
    f.mkdir(sDir.Data());
    f.cd(sDir);
  }
  this->Write();
  f.Close();
  
  if ( backup ) backup->cd();
}
//_____________________________________________________________________
TTreeSRedirector *AliTPCCalibRawBase::GetDebugStreamer(){
  //
  // Get Debug streamer
  // In case debug streamer not yet initialized and StreamLevel>0 create new one
  //
  if (fStreamLevel==0) return 0;
  if (fDebugStreamer) return fDebugStreamer;
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ","");
  fDebugStreamer = new TTreeSRedirector(dsName.Data());
  return fDebugStreamer;
}
