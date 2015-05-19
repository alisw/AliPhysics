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


/* $Id: AliADQADataMakerSim.cxx 23123 2007-12-18 09:08:18Z hristov $ */

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
#include "AliADdigit.h"
#include "AliADSDigit.h" 
#include "AliADhit.h"
#include "AliADQADataMakerSim.h"
#include "AliQAChecker.h"

ClassImp(AliADQADataMakerSim)
           
//____________________________________________________________________________ 
  AliADQADataMakerSim::AliADQADataMakerSim() : 
  AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kAD), "AD Quality Assurance Data Maker")

{
  // constructor

  
}

//____________________________________________________________________________ 
AliADQADataMakerSim::AliADQADataMakerSim(const AliADQADataMakerSim& qadm) :
  AliQADataMakerSim() 
{
  //copy constructor 

  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliADQADataMakerSim& AliADQADataMakerSim::operator = (const AliADQADataMakerSim& qadm )
{
  // Assign operator.
  this->~AliADQADataMakerSim();
  new(this) AliADQADataMakerSim(qadm);
  return *this;
}
//____________________________________________________________________________
void AliADQADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  ResetEventTrigClasses();
 // AliQAChecker::Instance()->Run(AliQAv1::kAD, task, list) ;
}

 
//____________________________________________________________________________ 
void AliADQADataMakerSim::InitHits()
{
 
  // create Hits histograms in Hits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I * h0 = new TH1I("hHitMultiplicity", "Hit multiplicity distribution in AD;# of Hits;Entries", 300, 0, 299) ; 
  h0->Sumw2() ;
  Add2HitsList(h0, 0, !expert, image) ;  
  
  TH1I * h1 = new TH1I("hHitCellNumber", "Hit cell distribution in AD;Cell;# of Hits", 16, 0, 16) ; 
  h1->Sumw2() ;
  Add2HitsList(h1, 1, !expert, image) ; 
  
  TH1I * h2 = new TH1I("hHitNPhotons", "Number of photons per hit in AD;# of Photons;Entries", 100000, 0, 100000) ; 
  h2->Sumw2() ;
  Add2HitsList(h2, 2, !expert, image) ;
  
  TH2I * h3 = new TH2I("hCellNPhotons", "Number of photons per cell in AD;Cell;# of Photons", 16, 0, 16, 100000, 0, 100000) ; 
  h2->Sumw2() ;
  Add2HitsList(h3, 3, !expert, image) ;
  
  TH2D * h4 = new TH2D("hCellTof", "Time of flight per cell in AD;Cell;Time of Flight [ns]", 16, 0, 16, 6000,40,100) ; 
  h2->Sumw2() ;
  Add2HitsList(h4, 4, !expert, image) ; 
   
  //
  ClonePerTrigClass(AliQAv1::kHITS); // this should be the last line    
}


//____________________________________________________________________________ 
void AliADQADataMakerSim::InitSDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I *fhSDigCharge[16]; 

  // create SDigits histograms in SDigits subdir
  TH1I * h0 = new TH1I("hSDigitMultiplicity", "SDigits multiplicity distribution in AD;# of Digits;Entries", 17,-0.5,16.5) ; 
  h0->Sumw2() ;
  Add2SDigitsList(h0, 0, !expert, image) ;
  
  TH2D * h1 = new TH2D("hSDigitChargePerPM", "SDigits total amplified charge per PM in AD;PM number;Charge [C]", 16, 0, 16, 1000,1e-13, 1e-9); 
  h1->Sumw2() ;
  Add2SDigitsList(h1, 1, !expert, image) ;
  
  //
  ClonePerTrigClass(AliQAv1::kSDIGITS); // this should be the last line
}



//____________________________________________________________________________ 
void AliADQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  // create Digits histograms in Digits subdir
  TH1I * h0 = new TH1I("hDigitMultiplicity", "Digits multiplicity distribution in AD;# of Digits;Entries", 17,-0.5,16.5) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0, !expert, image) ;
     
  TH2D * h1 = new TH2D("hDigitLeadingTimePerPM", "Leading time distribution per PM in AD;PM number;Leading Time [ns]",16,0,16, 3062, 0.976562, 300); 
  h1->Sumw2() ;
  Add2DigitsList(h1, 1, !expert, image) ; 
  
  TH2D * h2 = new TH2D("hDigitTimeWidthPerPM", "Time width distribution per PM in AD;PM number;Time width [ns]",16,0,16, 153, 2.343750, 121.875000); 
  h2->Sumw2() ;
  Add2DigitsList(h2, 2, !expert, image) ;
  
  TH2I * h3 = new TH2I("hDigitChargePerClockPerPM", "Charge array per PM in AD;PM number; Clock",16,0,16,21, -10.5, 10.5);
  h3->Sumw2();
  Add2DigitsList(h3, 3, !expert, image) ;
  
  TH1I * h4 = new TH1I("hDigitBBflagsAD","Number of BB flags in AD; # of BB flags; Entries",17,-0.5,16.5);
  h4->Sumw2();
  Add2DigitsList(h4, 4, !expert, image) ;
  
  TH1I * h5 = new TH1I("hDigitBBflagsADA","Number of BB flags in ADA; # of BB flags; Entries",9,-0.5,8.5);
  h5->Sumw2();
  Add2DigitsList(h5, 5, !expert, image) ;
  
  TH1I * h6 = new TH1I("hDigitBBflagsADC","Number of BB flags in ADC; # of BB flags; Entries",9,-0.5,8.5);
  h6->Sumw2();
  Add2DigitsList(h6, 6, !expert, image) ;
  
  TH2D * h7 = new TH2D("hDigitTotalChargePerPM", "Total Charge per PM in AD;PM number; Charge [ADC counts]",16,0,16,10000,0,10000);
  h7->Sumw2();
  Add2DigitsList(h7, 7, !expert, image) ;
  
  TH2I * h8 = new TH2I("hDigitMaxChargeClockPerPM", "Clock with maximum charge per PM in AD;PM number; Clock ",16,0,16,21, -10.5, 10.5);
  h8->Sumw2();
  Add2DigitsList(h8, 8, !expert, image) ;
   
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}


//____________________________________________________________________________
void AliADQADataMakerSim::MakeHits()
{
	//make QA data from Hits

  Int_t nhits = fHitsArray->GetEntriesFast();
  FillHitsData(0,nhits) ;    // fills Hit multiplicity
  for (Int_t ihit=0;ihit<nhits;ihit++) 
    {
	   AliADhit  * ADHit   = (AliADhit*) fHitsArray->UncheckedAt(ihit);
	   if (!ADHit) {
 	      AliError("The unchecked hit doesn't exist");
	      break;
	   }
	   FillHitsData(1,ADHit->GetCell());
	   FillHitsData(2,ADHit->GetNphot());
	   FillHitsData(3,ADHit->GetCell(),ADHit->GetNphot());
	   FillHitsData(4,ADHit->GetCell(),ADHit->GetTof());
	}
}


//____________________________________________________________________________

void AliADQADataMakerSim::MakeHits(TTree *hitTree)
{
  //fills QA histos for Hits
 if (fHitsArray)
   fHitsArray->Clear() ; 
  else 
    fHitsArray = new TClonesArray("AliADhit", 1000);
  
  TBranch * branch = hitTree->GetBranch("AD") ;
  if ( ! branch ) {
    AliWarning("AD branch in Hit Tree not found") ;
  } else {

   if (branch) {
      branch->SetAddress(&fHitsArray);
    }else{
      AliError("Branch AD hit not found");
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
	  AliADhit  * ADHit   = (AliADhit*) fHitsArray->UncheckedAt(ihit);
	  if (!ADHit) {
 	    AliError("The unchecked hit doesn't exist");
	    break;
	  }
	  FillHitsData(1,ADHit->GetCell());
	  FillHitsData(2,ADHit->GetNphot());
	  FillHitsData(3,ADHit->GetCell(),ADHit->GetNphot());
	  FillHitsData(4,ADHit->GetCell(),ADHit->GetTof());	 
	}
    }
  }
  //
  IncEvCountCycleHits();
  IncEvCountTotalHits();
  //
}



//____________________________________________________________________________
void AliADQADataMakerSim::MakeSDigits(TTree *sdigitTree)
{
    // makes data from Digit Tree
	
  if (fSDigitsArray)
    fSDigitsArray->Clear() ; 
  else 
    fSDigitsArray = new TClonesArray("AliADSDigit", 1000) ; 

    TBranch * branch = sdigitTree->GetBranch("ADSDigit") ;
    if ( ! branch ) {
         AliWarning("AD branch in SDigit Tree not found") ; 
    } else {
         branch->SetAddress(&fSDigitsArray) ;
         branch->GetEntry(0) ; 
         MakeSDigits() ; 
    }  
    //
    IncEvCountCycleDigits();
    IncEvCountTotalDigits();
    //    
}

//____________________________________________________________________________
void AliADQADataMakerSim::MakeSDigits()
{
  // makes data from SDigits

  FillSDigitsData(0,fSDigitsArray->GetEntriesFast()) ; 
  TIter next(fSDigitsArray) ; 
    AliADSDigit *ADSDigit ; 
    while ( (ADSDigit = dynamic_cast<AliADSDigit *>(next())) ) {
         Int_t   PMNumber  = ADSDigit->PMNumber();
	 Int_t   Nbins = ADSDigit->GetNBins();
	 Double_t totCharge = 0;
	
	 for(Int_t i = 0; i<Nbins; i++)totCharge += ADSDigit->GetCharges()[i];       
         FillSDigitsData(1, PMNumber, totCharge) ;
    }  
}

//____________________________________________________________________________
void AliADQADataMakerSim::MakeDigits()
{
  // makes data from Digits

  FillDigitsData(0,fDigitsArray->GetEntriesFast()) ; 
  TIter next(fDigitsArray) ; 
    AliADdigit *ADDigit ; 
    Int_t nBBflagsADA = 0;
    Int_t nBBflagsADC = 0;
    
    while ( (ADDigit = dynamic_cast<AliADdigit *>(next())) ) {
         Int_t totCharge = 0;
         Int_t   PMNumber  = ADDigit->PMNumber();

	 if(PMNumber<8 && ADDigit->GetBBflag()) nBBflagsADC++;
	 if(PMNumber>7 && ADDigit->GetBBflag()) nBBflagsADA++;
	 
	 Short_t adc[21];
	 for(Int_t iClock=0; iClock<21; iClock++) { 
	 adc[iClock]= ADDigit->ChargeADC(iClock);
	 FillDigitsData(3, PMNumber,(float)iClock-10,(float)adc[iClock]);
	 totCharge += adc[iClock];
	 }
	    
         FillDigitsData(1,PMNumber,ADDigit->Time()); 
	 FillDigitsData(2,PMNumber,ADDigit->Width());
	 FillDigitsData(7,PMNumber,totCharge);
	 FillDigitsData(8,PMNumber,TMath::LocMax(21,adc)-10); 
	 
    }
    FillDigitsData(4,nBBflagsADA+nBBflagsADC);
    FillDigitsData(5,nBBflagsADA);
    FillDigitsData(6,nBBflagsADC);  
}

//____________________________________________________________________________
void AliADQADataMakerSim::MakeDigits(TTree *digitTree)
{
    // makes data from Digit Tree
	
  if (fDigitsArray)
    fDigitsArray->Clear() ; 
  else 
    fDigitsArray = new TClonesArray("AliADdigit", 1000) ; 

    TBranch * branch = digitTree->GetBranch("ADDigit") ;
    if ( ! branch ) {
         AliWarning("AD branch in Digit Tree not found") ; 
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
void AliADQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle

}
