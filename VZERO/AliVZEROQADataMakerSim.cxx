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


/* $Id: AliVZEROQADataMakerSim.cxx 23123 2007-12-18 09:08:18Z hristov $ */

//---
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  Author : BC
//---

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TDirectory.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliVZEROdigit.h" 
#include "AliVZEROhit.h"
#include "AliVZEROQADataMakerSim.h"
#include "AliQAChecker.h"

ClassImp(AliVZEROQADataMakerSim)
           
//____________________________________________________________________________ 
  AliVZEROQADataMakerSim::AliVZEROQADataMakerSim() : 
  AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kVZERO), "VZERO Quality Assurance Data Maker")

{
  // constructor

  
}

//____________________________________________________________________________ 
AliVZEROQADataMakerSim::AliVZEROQADataMakerSim(const AliVZEROQADataMakerSim& qadm) :
  AliQADataMakerSim() 
{
  //copy constructor 

  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliVZEROQADataMakerSim& AliVZEROQADataMakerSim::operator = (const AliVZEROQADataMakerSim& qadm )
{
  // Assign operator.
  this->~AliVZEROQADataMakerSim();
  new(this) AliVZEROQADataMakerSim(qadm);
  return *this;
}
//____________________________________________________________________________
void AliVZEROQADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  ResetEventTrigClasses();
  AliQAChecker::Instance()->Run(AliQAv1::kVZERO, task, list) ;
}

 
//____________________________________________________________________________ 
void AliVZEROQADataMakerSim::InitHits()
{
 
  // create Hits histograms in Hits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I * h0 = new TH1I("hHitMultiplicity", "Hit multiplicity distribution in VZERO;# of Hits;Entries", 300, 0, 299) ; 
  h0->Sumw2() ;
  Add2HitsList(h0, 0, !expert, image) ;  
  
  TH1I * h1 = new TH1I("hHitCellNumber", "Hit cell distribution in VZERO;# of Hits;Entries", 80, 0, 79) ; 
  h1->Sumw2() ;
  Add2HitsList(h1, 1, !expert, image) ;  
  //
  ClonePerTrigClass(AliQAv1::kHITS); // this should be the last line    
}

//____________________________________________________________________________ 
void AliVZEROQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I *fhDigTDC[64]; 
  TH1I *fhDigADC[64]; 

  // create Digits histograms in Digits subdir
  TH1I * h0 = new TH1I("hDigitMultiplicity", "Digits multiplicity distribution in VZERO;# of Digits;Entries", 100, 0, 99) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0, !expert, image) ;
  
  for (Int_t i=0; i<64; i++)
    {
       fhDigTDC[i] = new TH1I(Form("hDigitTDC%d", i),Form("Digit TDC in cell %d; TDC value;Entries",i),300,0.,149.);
       fhDigADC[i]= new TH1I(Form("hDigitADC%d", i),Form("Digit ADC in cell %d;ADC value;Entries",i),1024,0.,1023.);
       
       Add2DigitsList(fhDigTDC[i],i+1, !expert, image);
       Add2DigitsList(fhDigADC[i],i+1+64, !expert, image);  
     }  
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}


//____________________________________________________________________________
void AliVZEROQADataMakerSim::MakeHits()
{
	//make QA data from Hits

  Int_t nhits = fHitsArray->GetEntriesFast();
  FillHitsData(0,nhits) ;    // fills Hit multiplicity
  for (Int_t ihit=0;ihit<nhits;ihit++) 
    {
	   AliVZEROhit  * VZEROHit   = (AliVZEROhit*) fHitsArray->UncheckedAt(ihit);
	   if (!VZEROHit) {
 	      AliError("The unchecked hit doesn't exist");
	      break;
	   }
	   FillHitsData(1,VZEROHit->Cell());
	}
}


//____________________________________________________________________________

void AliVZEROQADataMakerSim::MakeHits(TTree *hitTree)
{
  //fills QA histos for Hits
 if (fHitsArray)
   fHitsArray->Clear() ; 
  else 
    fHitsArray = new TClonesArray("AliVZEROhit", 1000);
  
  TBranch * branch = hitTree->GetBranch("VZERO") ;
  if ( ! branch ) {
    AliWarning("VZERO branch in Hit Tree not found") ;
  } else {

   if (branch) {
      branch->SetAddress(&fHitsArray);
    }else{
      AliError("Branch VZERO hit not found");
      exit(111);
    } 
    // Check id histograms already created for this Event Specie
    if ( ! GetHitsData(0) )
      InitHits() ;
    
    Int_t ntracks    = (Int_t) hitTree->GetEntries();
    
    if (ntracks<=0) return;
    // Start loop on tracks in the hits containers
    for (Int_t track=0; track<ntracks;track++) {
      branch->GetEntry(track);
      Int_t nhits = fHitsArray->GetEntriesFast();
      FillHitsData(0,nhits) ;    // fills Hit multiplicity
      for (Int_t ihit=0;ihit<nhits;ihit++) 
	{
	  AliVZEROhit  * VZEROHit   = (AliVZEROhit*) fHitsArray->UncheckedAt(ihit);
	  if (!VZEROHit) {
 	    AliError("The unchecked hit doesn't exist");
	    break;
	  }
	  FillHitsData(1,VZEROHit->Cell());	 
	}
    }
  }
  //
  IncEvCountCycleHits();
  IncEvCountTotalHits();
  //
}


//____________________________________________________________________________
void AliVZEROQADataMakerSim::MakeDigits()
{
  // makes data from Digits

  FillDigitsData(0,fDigitsArray->GetEntriesFast()) ; 
  TIter next(fDigitsArray) ; 
    AliVZEROdigit *VZERODigit ; 
    while ( (VZERODigit = dynamic_cast<AliVZEROdigit *>(next())) ) {
         Int_t   PMNumber  = VZERODigit->PMNumber();         
         FillDigitsData(PMNumber +1, VZERODigit->Time()) ;    // in 100 of picoseconds
	 FillDigitsData(PMNumber +1+64, VZERODigit->ADC()) ;
    }  
}


//____________________________________________________________________________
void AliVZEROQADataMakerSim::MakeDigits(TTree *digitTree)
{
    // makes data from Digit Tree
	
  if (fDigitsArray)
    fDigitsArray->Clear() ; 
  else 
    fDigitsArray = new TClonesArray("AliVZEROdigit", 1000) ; 

    TBranch * branch = digitTree->GetBranch("VZERODigit") ;
    if ( ! branch ) {
         AliWarning("VZERO branch in Digit Tree not found") ; 
    } else {
         branch->SetAddress(&fDigitsArray) ;
         branch->GetEntry(0) ; 
         MakeDigits() ; 
    }  
    //
    IncEvCountCycleDigits();
    IncEvCountTotalDigits();
    //    
}


//____________________________________________________________________________
void AliVZEROQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

}
