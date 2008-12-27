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

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH1I.h> 
#include <TH2F.h> 
#include <TTree.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDCaloCluster.h"
#include "AliLog.h"
#include "AliPHOSDigit.h"
#include "AliPHOSHit.h"
#include "AliPHOSQADataMakerSim.h"
#include "AliQAChecker.h"

ClassImp(AliPHOSQADataMakerSim)
           
//____________________________________________________________________________ 
AliPHOSQADataMakerSim::AliPHOSQADataMakerSim() : 
  AliQADataMakerSim(AliQA::GetDetName(AliQA::kPHOS), "PHOS Quality Assurance Data Maker"),
  fHits(0x0)
{
  // ctor
  fHits = new TClonesArray("AliPHOSHit", 1000);
}

//____________________________________________________________________________ 
AliPHOSQADataMakerSim::AliPHOSQADataMakerSim(const AliPHOSQADataMakerSim& qadm) :
  AliQADataMakerSim(),
  fHits(0x0)
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
  fHits = new TClonesArray("AliPHOSHit", 1000);
}

//__________________________________________________________________
AliPHOSQADataMakerSim& AliPHOSQADataMakerSim::operator = (const AliPHOSQADataMakerSim& qadm )
{
  // Assign operator.
  this->~AliPHOSQADataMakerSim();
  new(this) AliPHOSQADataMakerSim(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliPHOSQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kPHOS, task, list) ;  
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerSim::InitHits()
{
  // create Hits histograms in Hits subdir
  Bool_t expert   = kTRUE ; 
  TH1F * h0 = new TH1F("hPhosHits",    "Hits energy distribution in PHOS",       100, 0., 100.) ; 
  h0->Sumw2() ;
  Add2HitsList(h0, kHits, !expert) ;
  TH1I * h1 = new TH1I("hPhosHitsMul", "Hits multiplicity distribution in PHOS", 500, 0., 10000) ; 
  h1->Sumw2() ;
  Add2HitsList(h1, kHitsMul, !expert) ;
  
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  Bool_t expert   = kTRUE ; 
  TH1I * h0 = new TH1I("hPhosDigits",    "Digits amplitude distribution in PHOS",    500, 0, 1000) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, kDigits, !expert) ;
  TH1I * h1 = new TH1I("hPhosDigitsMul", "Digits multiplicity distribution in PHOS", 2000, 0, 10000) ; 
  h1->Sumw2() ;
  Add2DigitsList(h1, kDigitsMul, !expert) ;
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerSim::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
  Bool_t expert   = kTRUE ; 
  TH1F * h0 = new TH1F("hPhosSDigits",    "SDigits energy distribution in PHOS",       500, 0., 1000.) ; 
  h0->Sumw2() ;
  Add2SDigitsList(h0, kSDigits, !expert) ;
  TH1I * h1 = new TH1I("hPhosSDigitsMul", "SDigits multiplicity distribution in PHOS", 500, 0,  1000) ; 
  h1->Sumw2() ;
  Add2SDigitsList(h1, kSDigitsMul, !expert) ;
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeHits()
{
  //make QA data from Hits
  
  TIter next(fHits) ; 
  AliPHOSHit * hit ; 
  while ( (hit = dynamic_cast<AliPHOSHit *>(next())) ) {
    GetHitsData(kHits)->Fill( hit->GetEnergy()) ;
  }
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeHits(TTree * hitTree)
{
  // make QA data from Hit Tree
  
  TBranch * branch = hitTree->GetBranch("PHOS") ;
  if ( ! branch ) {
    AliWarning("PHOS branch in Hit Tree not found") ; 
  } else {
    Int_t nHits = 0;
    branch->SetAddress(&fHits) ;
    for (Int_t ientry = 0 ; ientry < branch->GetEntries() ; ientry++) {
      branch->GetEntry(ientry) ;
      nHits += fHits->GetEntriesFast();
      MakeHits() ; 
      fHits->Clear();
    } 	
    GetHitsData(1)->Fill(nHits) ;
  }
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeDigits(TClonesArray * digits)
{
  // makes data from Digits

    GetDigitsData(1)->Fill(digits->GetEntriesFast()) ; 
    TIter next(digits) ; 
    AliPHOSDigit * digit ; 
    while ( (digit = dynamic_cast<AliPHOSDigit *>(next())) ) {
      GetDigitsData(kDigits)->Fill( digit->GetEnergy()) ;
    }  
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeDigits(TTree * digitTree)
{
	// makes data from Digit Tree
	TClonesArray * digits = new TClonesArray("AliPHOSDigit", 1000) ; 

	TBranch * branch = digitTree->GetBranch("PHOS") ;
	if ( ! branch ) {
		AliWarning("PHOS branch in Digit Tree not found") ; 
	} else {
		branch->SetAddress(&digits) ;
		branch->GetEntry(0) ; 
		MakeDigits(digits) ; 
	}
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeSDigits(TClonesArray * sdigits)
{
  // makes data from SDigits
  
	GetSDigitsData(1)->Fill(sdigits->GetEntriesFast()) ; 
    TIter next(sdigits) ; 
    AliPHOSDigit * sdigit ; 
    while ( (sdigit = dynamic_cast<AliPHOSDigit *>(next())) ) {
      GetSDigitsData(kSDigits)->Fill( sdigit->GetEnergy()) ;
    } 
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeSDigits(TTree * sdigitTree)
{
	// makes data from SDigit Tree
	TClonesArray * sdigits = new TClonesArray("AliPHOSDigit", 1000) ; 

	TBranch * branch = sdigitTree->GetBranch("PHOS") ;
	if ( ! branch ) {
		AliWarning("PHOS branch in SDigit Tree not found") ; 
	} else {
		branch->SetAddress(&sdigits) ;
		branch->GetEntry(0) ;
		MakeSDigits(sdigits) ; 
	}
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
