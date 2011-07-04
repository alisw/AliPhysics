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

//---
//  Produces the data needed to calculate the quality assurance. 
//  Alla.Maevskaya@cern.ch
//  
//---

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h> 
#include <TDirectory.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliT0digit.h" 
#include "AliT0hit.h"
#include "AliT0RecPoint.h"
#include "AliT0QADataMakerSim.h"
#include "AliQAChecker.h"
#include "AliT0RawReader.h"

#include <Riostream.h>

ClassImp(AliT0QADataMakerSim)
           
//____________________________________________________________________________ 
  AliT0QADataMakerSim::AliT0QADataMakerSim() : 
  AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kT0), "T0 Quality Assurance Data Maker")

{
  // ctor
 //   fDetectorDir = fOutput->GetDirectory(GetName()) ;  
//   if (!fDetectorDir) 
//     fDetectorDir = fOutput->mkdir(GetName()) ;  
}

//____________________________________________________________________________ 
AliT0QADataMakerSim::AliT0QADataMakerSim(const AliT0QADataMakerSim& qadm) :
  AliQADataMakerSim() 
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliT0QADataMakerSim& AliT0QADataMakerSim::operator = (const AliT0QADataMakerSim& qadm )
{
  // Equal operator.
  this->~AliT0QADataMakerSim();
  new(this) AliT0QADataMakerSim(qadm);
  return *this; 
}
//____________________________________________________________________________
void AliT0QADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kT0, task, list) ;
}

//____________________________________________________________________________
void AliT0QADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

}
 
//____________________________________________________________________________ 
void AliT0QADataMakerSim::InitHits()
{
  // create Hits histograms in Hits subdir
  // create Hits histograms in Hits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TString timename;
  
  TH2F *fhHitsTimeA = new TH2F("hHitsTimeA", "Hits Efficiency;#PMT; Time [ns];", 13, 12, 25, 100,12,15 );
  fhHitsTimeA->SetOption("COLZ");
  Add2HitsList(fhHitsTimeA,0, !expert, image);
  TH2F *fhHitsTimeC = new TH2F("hHitsTimeC", "Hits Efficiency;#PMT; Time [ns];", 13, 0, 13, 100,2,5 );
  fhHitsTimeC->SetOption("COLZ");
  Add2HitsList(fhHitsTimeC,1, !expert, image);
  //
  ClonePerTrigClass(AliQAv1::kHITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliT0QADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH2F * fhDigCFD = new TH2F("fhDigCFD", " CFD digits; #PMT; CFD time [#channel]",25,-0.5,24.5,100,0,1000);
  fhDigCFD->SetOption("COLZ");
  Add2DigitsList( fhDigCFD,0);
  TH2F *fhDigLEDamp = new TH2F("fhDigLEDamp", " LED-CFD digits; #PMT; amplitude  LED-CFD [#channel]",25,-0.5,24.5,100,100,1000);
  fhDigLEDamp->SetOption("COLZ");
  Add2DigitsList( fhDigLEDamp,1, !expert, image);
  TH2F * fhDigQTC = new TH2F("fhDigQTC", " QTC digits; #PMT; amplitude QTC [#channel]",25,-0.5,24.5,200,500,10000);
  fhDigQTC->SetOption("COLZ");
  Add2DigitsList( fhDigQTC,2, !expert, image);
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line 
}

//____________________________________________________________________________

void AliT0QADataMakerSim::MakeHits(TTree *hitTree)
{
  //fills QA histos for Hits
  if (fHitsArray) 
    fHitsArray->Clear() ; 
  else 
    fHitsArray = new TClonesArray("AliT0hit", 1000);
  
  TBranch * branch = hitTree->GetBranch("T0") ;
  if ( ! branch ) {
    AliWarning("T0 branch in Hit Tree not found") ;
  } else {

   if (branch) {
      branch->SetAddress(&fHitsArray);
    }else{
      AliError("Branch T0 hit not found");
      exit(111);
    } 
    Int_t ntracks    = (Int_t) hitTree->GetEntries();
    
    if (ntracks<=0) return;
    // Start loop on tracks in the hits containers

    for (Int_t track=0; track<ntracks;track++) {
      branch->GetEntry(track);
      Int_t nhits = fHitsArray->GetEntriesFast();
      for (Int_t ihit=0;ihit<nhits;ihit++) 
	{
	  AliT0hit  * startHit   = (AliT0hit*) fHitsArray->UncheckedAt(ihit);
	  if (!startHit) {
 	    AliError("The unchecked hit doesn't exist");
	    continue;
	  }
	  Int_t pmt=startHit->Pmt();
	  Float_t time = 0.001 * startHit->Time();
	  if(pmt<13) FillHitsData(1,pmt,time) ;
	  if(pmt>12) FillHitsData(0,pmt,time) ;
	}
    }
  }
  //
  IncEvCountCycleHits();
  IncEvCountTotalHits();
  //
}

//____________________________________________________________________________
void AliT0QADataMakerSim::MakeDigits( TTree *digitsTree)
{
  //fills QA histos for Digits
 
  TArrayI *digCFD = new TArrayI(24);
  TArrayI *digLED = new TArrayI(24);
  TArrayI *digQT0 = new TArrayI(24);
  TArrayI *digQT1 = new TArrayI(24);
  Int_t refpoint=0;

  TBranch *brDigits=digitsTree->GetBranch("T0");
  AliT0digit *fDigits = new AliT0digit() ;
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    AliError(Form("EXEC Branch T0 digits not found"));
     return;
  }

  digitsTree->GetEvent(0);
  digitsTree->GetEntry(0);
  brDigits->GetEntry(0);
  fDigits->GetTimeCFD(*digCFD);
  fDigits->GetTimeLED(*digLED);
  fDigits->GetQT0(*digQT0);
  fDigits->GetQT1(*digQT1);
  refpoint = fDigits->RefPoint();

   for (Int_t i=0; i<24; i++)
    {
      if (digCFD->At(i)>0) {
	Int_t cfd=digCFD->At(i)- refpoint;
	FillDigitsData(0, i,cfd);
	FillDigitsData(1, i,(digLED->At(i) - digCFD->At(i)));
	FillDigitsData(2, i, (digQT1->At(i) - digQT0->At(i)));

      }
    }  
      
  delete digCFD;
  delete digLED;
  delete digQT0;
  delete digQT1;
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
}
