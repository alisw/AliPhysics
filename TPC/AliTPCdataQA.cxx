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


/* $Id$ */


//Root includes
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TFile.h>
//AliRoot includes
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"
#include "AliTPCRawStream.h"
#include "AliTPCCalROC.h"
#include "AliTPCROC.h"
#include "AliMathBase.h"
#include "TTreeStream.h"
#include "AliTPCRawStreamFast.h"

//date
#include "event.h"
#include "AliTPCCalPad.h"

//header file
#include "AliTPCdataQA.h"



ClassImp(AliTPCdataQA)

AliTPCdataQA::AliTPCdataQA() : /*FOLD00*/
  TObject(),
  fFirstTimeBin(60),
  fLastTimeBin(1000),
  fAdcMin(1),
  fAdcMax(100),
  fOldRCUformat(kTRUE),
  fROC(AliTPCROC::Instance()),
  fMapping(NULL),
  fMaxCharge(0),
  fOverThreshold0(0),
  fOverThreshold5(0),
  fOverThreshold10(0),
  fOverThreshold20(0),
  fOverThreshold30(0),
  fEventCounter(0)
{
  //
  // default constructor
  //
}


//_____________________________________________________________________
AliTPCdataQA::AliTPCdataQA(const AliTPCdataQA &ped) : /*FOLD00*/
  TObject(ped),
  fFirstTimeBin(ped.GetFirstTimeBin()),
  fLastTimeBin(ped.GetLastTimeBin()),
  fAdcMin(ped.GetAdcMin()),
  fAdcMax(ped.GetAdcMax()),
  fOldRCUformat(ped.fOldRCUformat),
  fROC(AliTPCROC::Instance()),
  fMapping(NULL)
{
  //
  // copy constructor
  //
 
}


//_____________________________________________________________________
AliTPCdataQA& AliTPCdataQA::operator = (const  AliTPCdataQA &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTPCdataQA(source);

  return *this;
}


//_____________________________________________________________________
AliTPCdataQA::~AliTPCdataQA() /*FOLD00*/
{
  //
  // destructor
  //

  // do not delete fMapping, because we do not own it.

}




//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEventFast(AliTPCRawStreamFast *rawStreamFast)
{
  //
  // Event Processing loop - AliTPCRawStream
  //
  Bool_t withInput = kFALSE;

  while ( rawStreamFast->NextDDL() ){
      while ( rawStreamFast->NextChannel() ){
	  Int_t isector  = rawStreamFast->GetSector();                       //  current sector
	  Int_t iRow     = rawStreamFast->GetRow();                          //  current row
	  Int_t iPad     = rawStreamFast->GetPad();                          //  current pad
	  Int_t startTbin = (Int_t)rawStreamFast->GetStartTimeBin();
          Int_t endTbin = (Int_t)rawStreamFast->GetEndTimeBin();

	  while ( rawStreamFast->NextBunch() ){
	      for (Int_t iTimeBin = startTbin; iTimeBin < endTbin; iTimeBin++){
		  Float_t signal=(Float_t)rawStreamFast->GetSignals()[iTimeBin-startTbin];
		  Update(isector,iRow,iPad,iTimeBin+1,signal);
		  withInput = kTRUE;
	      }
	  }
      }
  }

  return withInput;
}
//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEventFast(AliRawReader *rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //
  AliTPCRawStreamFast *rawStreamFast = new AliTPCRawStreamFast(rawReader, (AliAltroMapping**)fMapping);
  Bool_t res=ProcessEventFast(rawStreamFast);
  delete rawStreamFast;
  return res;
}

//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEvent(AliTPCRawStream *rawStream)
{
  //
  // Event Processing loop - AliTPCRawStream
  //

  rawStream->SetOldRCUFormat(fOldRCUformat);

  Bool_t withInput = kFALSE;

  while (rawStream->Next()) {

    Int_t iSector  = rawStream->GetSector();      //  current ROC
    Int_t iRow     = rawStream->GetRow();         //  current row
    Int_t iPad     = rawStream->GetPad();         //  current pad
    Int_t iTimeBin = rawStream->GetTime();        //  current time bin
    Float_t signal = rawStream->GetSignal();      //  current ADC signal
    
    Update(iSector,iRow,iPad,iTimeBin,signal);
    withInput = kTRUE;
  }

  return withInput;
}


//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEvent(AliRawReader *rawReader)
{
  //
  //  Event processing loop - AliRawReader
  //

  // if fMapping is NULL the rawstream will crate its own mapping
  AliTPCRawStream rawStream(rawReader, (AliAltroMapping**)fMapping);
  rawReader->Select("TPC");
  return ProcessEvent(&rawStream);
}


//_____________________________________________________________________
Bool_t AliTPCdataQA::ProcessEvent(eventHeaderStruct *event)
{
  //
  //  process date event
  //

  AliRawReader *rawReader = new AliRawReaderDate((void*)event);
  Bool_t result=ProcessEvent(rawReader);
  delete rawReader;
  return result;
}



//_____________________________________________________________________
void AliTPCdataQA::DumpToFile(const Char_t *filename, const Char_t *dir, Bool_t append) /*FOLD00*/
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
Int_t AliTPCdataQA::Update(const Int_t icsector, /*FOLD00*/
				  const Int_t icRow,
				  const Int_t icPad,
				  const Int_t icTimeBin,
				  const Float_t csignal)
{
  //
  // Signal filling method
  //
  if (icTimeBin<fFirstTimeBin) return 0;
  if (icTimeBin>fLastTimeBin) return 0;
  if (!fMaxCharge) fMaxCharge = new AliTPCCalPad("MaxCharge","MaxCharge");
  if (!fOverThreshold0) fOverThreshold0 = new AliTPCCalPad("OverThreshold0","OverThreshold0");
  if (!fOverThreshold5) fOverThreshold5 = new AliTPCCalPad("OverThreshold5","OverThreshold5");
  if (!fOverThreshold10) fOverThreshold10 = new AliTPCCalPad("OverThreshold10","OverThreshold10");
  if (!fOverThreshold20) fOverThreshold20 = new AliTPCCalPad("OverThreshold20","OverThreshold20");
  if (!fOverThreshold30) fOverThreshold30 = new AliTPCCalPad("OverThreshold30","OverThreshold30");
  //
  if (csignal>fMaxCharge->GetCalROC(icsector)->GetValue(icRow, icPad)){
    fMaxCharge->GetCalROC(icsector)->SetValue(icRow, icPad,csignal);
  }
  
  //
  if (csignal>0){
    Int_t count = fOverThreshold0->GetCalROC(icsector)->GetValue(icRow, icPad);
    fOverThreshold0->GetCalROC(icsector)->SetValue(icRow, icPad,count+1);
  };
  //
  if (csignal>5){
    Int_t count = fOverThreshold5->GetCalROC(icsector)->GetValue(icRow, icPad);
    fOverThreshold5->GetCalROC(icsector)->SetValue(icRow, icPad,count+1);
  };
  if (csignal>10){
    Int_t count = fOverThreshold10->GetCalROC(icsector)->GetValue(icRow, icPad);
    fOverThreshold10->GetCalROC(icsector)->SetValue(icRow, icPad,count+1);
  };
  if (csignal>20){
    Int_t count = fOverThreshold20->GetCalROC(icsector)->GetValue(icRow, icPad);
    fOverThreshold20->GetCalROC(icsector)->SetValue(icRow, icPad,count+1);
  };
  if (csignal>30){
    Int_t count = fOverThreshold30->GetCalROC(icsector)->GetValue(icRow, icPad);
    fOverThreshold30->GetCalROC(icsector)->SetValue(icRow, icPad,count+1);
  };

  return 0;
}
