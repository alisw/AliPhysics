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

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.

  Based on PHOS code written by
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
#include "AliEMCALDigit.h"
#include "AliEMCALHit.h"
#include "AliEMCALQADataMakerSim.h"
#include "AliQAChecker.h"
#include "AliEMCALSDigitizer.h"

ClassImp(AliEMCALQADataMakerSim)
           
//____________________________________________________________________________ 
  AliEMCALQADataMakerSim::AliEMCALQADataMakerSim() : 
  AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kEMCAL), "EMCAL Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliEMCALQADataMakerSim::AliEMCALQADataMakerSim(const AliEMCALQADataMakerSim& qadm) :
  AliQADataMakerSim()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliEMCALQADataMakerSim& AliEMCALQADataMakerSim::operator = (const AliEMCALQADataMakerSim& qadm )
{
  // Assign operator.
  this->~AliEMCALQADataMakerSim();
  new(this) AliEMCALQADataMakerSim(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliEMCALQADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kEMCAL, task, list) ;  
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerSim::InitHits()
{
  // create Hits histograms in Hits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F * h0 = new TH1F("hEmcalHits",    "Hits energy distribution in EMCAL",       200, 0., 2.) ; //GeV
  h0->Sumw2() ;
  Add2HitsList(h0, 0, !expert, image) ;
  TH1I * h1  = new TH1I("hEmcalHitsMul", "Hits multiplicity distribution in EMCAL", 1000, 0, 10000) ; 
  h1->Sumw2() ;
  Add2HitsList(h1, 1, !expert, image) ;
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1I * h0 = new TH1I("hEmcalDigits",    "Digits amplitude distribution in EMCAL",    500, 0, 500) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0, !expert, image) ;
  TH1I * h1 = new TH1I("hEmcalDigitsMul", "Digits multiplicity distribution in EMCAL", 200, 0, 2000) ; 
  h1->Sumw2() ;
  Add2DigitsList(h1, 1, !expert, image) ;
}

//____________________________________________________________________________ 
void AliEMCALQADataMakerSim::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F * h0 = new TH1F("hEmcalSDigits",    "SDigits energy distribution in EMCAL",       200, 0., 20.) ; 
  h0->Sumw2() ;
  Add2SDigitsList(h0, 0, !expert, image) ;
  TH1I * h1 = new TH1I("hEmcalSDigitsMul", "SDigits multiplicity distribution in EMCAL", 500, 0,  5000) ; 
  h1->Sumw2() ;
  Add2SDigitsList(h1, 1, !expert, image) ;
}

//____________________________________________________________________________
void AliEMCALQADataMakerSim::MakeHits(TClonesArray * hits)
{
  //make QA data from Hits

  GetHitsData(1)->Fill(hits->GetEntriesFast()) ; 
  TIter next(hits) ; 
  AliEMCALHit * hit ; 
  while ( (hit = dynamic_cast<AliEMCALHit *>(next())) ) {
    GetHitsData(0)->Fill( hit->GetEnergy()) ;
  }

}

//____________________________________________________________________________
void AliEMCALQADataMakerSim::MakeHits(TTree * hitTree)
{
  // make QA data from Hit Tree
  
  TClonesArray * hits = new TClonesArray("AliEMCALHit", 1000);
  
  TBranch * branch = hitTree->GetBranch("EMCAL") ;
  if ( ! branch ) {
    AliWarning("EMCAL branch in Hit Tree not found") ; 
  } else {
    TClonesArray * tmp =  new TClonesArray("AliEMCALHit", 1000) ;
    branch->SetAddress(&tmp) ;
    Int_t index = 0 ;  
    for (Int_t ientry = 0 ; ientry < branch->GetEntries() ; ientry++) {
      branch->GetEntry(ientry) ; 
      for (Int_t ihit = 0 ; ihit < tmp->GetEntries() ; ihit++) {
	AliEMCALHit * hit = dynamic_cast<AliEMCALHit *> (tmp->At(ihit)) ; 
	new((*hits)[index]) AliEMCALHit(*hit) ; 
	index++ ;
      } 
    } 	
    tmp->Delete() ; 
    delete tmp ; 
    MakeHits(hits) ; 
  }
}

//____________________________________________________________________________
void AliEMCALQADataMakerSim::MakeDigits(TClonesArray * digits)
{
  // makes data from Digits

  GetDigitsData(1)->Fill(digits->GetEntriesFast()) ; 
  TIter next(digits) ; 
  AliEMCALDigit * digit ; 
  while ( (digit = dynamic_cast<AliEMCALDigit *>(next())) ) {
    GetDigitsData(0)->Fill( digit->GetAmp()) ;
  }  

}

//____________________________________________________________________________
void AliEMCALQADataMakerSim::MakeDigits(TTree * digitTree)
{
  // makes data from Digit Tree
  TClonesArray * digits = new TClonesArray("AliEMCALDigit", 1000) ; 
  
  TBranch * branch = digitTree->GetBranch("EMCAL") ;
  if ( ! branch ) {
    AliWarning("EMCAL branch in Digit Tree not found") ; 
  } else {
    branch->SetAddress(&digits) ;
    branch->GetEntry(0) ; 
    MakeDigits(digits) ; 
  }

}

//____________________________________________________________________________
void AliEMCALQADataMakerSim::MakeSDigits(TClonesArray * sdigits)
{
  // makes data from SDigits
  //Need a copy of the SDigitizer to calibrate the sdigit amplitude to
  //energy in GeV
  AliEMCALSDigitizer* sDigitizer = new AliEMCALSDigitizer();

  GetSDigitsData(1)->Fill(sdigits->GetEntriesFast()) ; 
  TIter next(sdigits) ; 
  AliEMCALDigit * sdigit ; 
  while ( (sdigit = dynamic_cast<AliEMCALDigit *>(next())) ) {
    GetSDigitsData(0)->Fill( sDigitizer->Calibrate(sdigit->GetAmp())) ;
  } 

  delete sDigitizer;
}

//____________________________________________________________________________
void AliEMCALQADataMakerSim::MakeSDigits(TTree * sdigitTree)
{
  // makes data from SDigit Tree
  TClonesArray * sdigits = new TClonesArray("AliEMCALDigit", 1000) ; 
  
  TBranch * branch = sdigitTree->GetBranch("EMCAL") ;
  if ( ! branch ) {
    AliWarning("EMCAL branch in SDigit Tree not found") ; 
  } else {
    branch->SetAddress(&sdigits) ;
    branch->GetEntry(0) ;
    MakeSDigits(sdigits) ; 
  }

}

//____________________________________________________________________________ 
void AliEMCALQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
