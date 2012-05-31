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
  AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kPHOS), "PHOS Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliPHOSQADataMakerSim::AliPHOSQADataMakerSim(const AliPHOSQADataMakerSim& qadm) :
  AliQADataMakerSim()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
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
void AliPHOSQADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kPHOS, task, list) ;  
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerSim::InitHits()
{
  // create Hits histograms in Hits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  TH1F * h0 = new TH1F("hPhosHits",    "Hits energy distribution in PHOS;Energy [MeV];Counts",       100, 0., 100.) ; 
  h0->Sumw2() ;
  Add2HitsList(h0, kHits, !expert, image) ;
  TH1I * h1 = new TH1I("hPhosHitsMul", "Hits multiplicity distribution in PHOS;# of Hits;Entries", 500, 0., 10000) ; 
  h1->Sumw2() ;
  Add2HitsList(h1, kHitsMul, !expert, image) ;
  //
  ClonePerTrigClass(AliQAv1::kHITS); // this should be the last line  
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  TH1I * h0 = new TH1I("hPhosDigits",    "Digits amplitude distribution in PHOS;Amplitude [ADC counts];Counts",    500, 0, 1000) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, kDigits, !expert, image) ;
  TH1I * h1 = new TH1I("hPhosDigitsMul", "Digits multiplicity distribution in PHOS;# of Digits;Entries", 2000, 0, 10000) ; 
  h1->Sumw2() ;
  Add2DigitsList(h1, kDigitsMul, !expert, image) ;
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerSim::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  TH1F * h0 = new TH1F("hPhosSDigits",    "SDigits energy distribution in PHOS; Energy [MeV];Counts",       500, 0., 1000.) ; 
  h0->Sumw2() ;
  Add2SDigitsList(h0, kSDigits, !expert, image) ;
  TH1I * h1 = new TH1I("hPhosSDigitsMul", "SDigits multiplicity distribution in PHOS;# of SDigits;Entries", 500, 0,  1000) ; 
  h1->Sumw2() ;
  Add2SDigitsList(h1, kSDigitsMul, !expert, image) ;
  //
  ClonePerTrigClass(AliQAv1::kSDIGITS); // this should be the last line
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeHits()
{
  //make QA data from Hits
  
  TIter next(fHitsArray) ; 
  AliPHOSHit * hit ; 
  while ( (hit = dynamic_cast<AliPHOSHit *>(next())) ) {
    FillHitsData(kHits, hit->GetEnergy()) ;
  }
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeHits(TTree * hitTree)
{
  // make QA data from Hit Tree
  
  if (fHitsArray)
    fHitsArray->Clear() ; 
  else
    fHitsArray = new TClonesArray("AliPHOSHit", 1000);

  TBranch * branch = hitTree->GetBranch("PHOS") ;
  if ( ! branch ) { AliWarning("PHOS branch in Hit Tree not found"); return;}
  //
  Int_t nHits = 0;
  branch->SetAddress(&fHitsArray) ;
  for (Int_t ientry = 0 ; ientry < branch->GetEntries() ; ientry++) {
    branch->GetEntry(ientry) ;
    nHits += fHitsArray->GetEntriesFast();
    MakeHits() ; 
    fHitsArray->Clear();
  } 	
  FillHitsData(1,nHits) ;
  //  
  IncEvCountCycleHits();
  IncEvCountTotalHits();
  //
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeDigits()
{
  // makes data from Digits
 
  FillDigitsData(1,fDigitsArray->GetEntriesFast()) ; 
  TIter next(fDigitsArray) ; 
  AliPHOSDigit * digit ; 
  while ( (digit = dynamic_cast<AliPHOSDigit *>(next())) ) {
    FillDigitsData(kDigits, digit->GetEnergy()) ;
  }  
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeDigits(TTree * digitTree)
{
  // makes data from Digit Tree
  if (fDigitsArray) 
    fDigitsArray->Clear() ; 
  else
    fDigitsArray = new TClonesArray("AliPHOSDigit", 1000) ; 
  
  TBranch * branch = digitTree->GetBranch("PHOS") ;
  if ( ! branch ) {AliWarning("PHOS branch in Digit Tree not found"); return;}
  branch->SetAddress(&fDigitsArray) ;
  branch->GetEntry(0) ; 
  MakeDigits() ; 
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeSDigits()
{
  // makes data from SDigits

  FillSDigitsData(1,fSDigitsArray->GetEntriesFast()) ; 
  TIter next(fSDigitsArray) ; 
  AliPHOSDigit * sdigit ; 
  while ( (sdigit = dynamic_cast<AliPHOSDigit *>(next())) ) {
    FillSDigitsData(kSDigits, sdigit->GetEnergy()) ;
  } 
}

//____________________________________________________________________________
void AliPHOSQADataMakerSim::MakeSDigits(TTree * sdigitTree)
{
	// makes data from SDigit Tree
  if (fSDigitsArray) 
    fSDigitsArray->Clear() ; 
  else
    fSDigitsArray = new TClonesArray("AliPHOSDigit", 1000) ; 
  
  TBranch * branch = sdigitTree->GetBranch("PHOS") ;
  if ( ! branch ) {AliWarning("PHOS branch in SDigit Tree not found"); return;}
  branch->SetAddress(&fSDigitsArray) ;
  branch->GetEntry(0) ;
  MakeSDigits() ; 
  //
  IncEvCountCycleSDigits();
  IncEvCountTotalSDigits();
  //
}

//____________________________________________________________________________ 
void AliPHOSQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
