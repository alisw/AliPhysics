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
  B.K. Nandi
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

#include "AliLog.h"
#include "AliPMDhit.h"
#include "AliPMDsdigit.h"
#include "AliPMDdigit.h"
#include "AliPMDQADataMakerSim.h"
#include "AliQAChecker.h"

ClassImp(AliPMDQADataMakerSim)
           
//____________________________________________________________________________ 
AliPMDQADataMakerSim::AliPMDQADataMakerSim() : 
    AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kPMD), "PMD Quality Assurance Data Maker")
{
    // ctor
}

//____________________________________________________________________________ 
AliPMDQADataMakerSim::AliPMDQADataMakerSim(const AliPMDQADataMakerSim& qadm) :
    AliQADataMakerSim()
{
    //copy ctor 
    SetName((const char*)qadm.GetName()) ; 
    SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliPMDQADataMakerSim& AliPMDQADataMakerSim::operator = (const AliPMDQADataMakerSim& qadm )
{
    // Assign operator.
    this->~AliPMDQADataMakerSim();
    new(this) AliPMDQADataMakerSim(qadm);
    return *this;
}
 
//____________________________________________________________________________ 
void AliPMDQADataMakerSim::InitHits()
{
    // create Hits histograms in Hits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *h0 = new TH1F("hPreHitsEdep","Hits energy distribution PRE(PMD);Energy [keV];Counts", 500, 0., 500.); 
  h0->Sumw2() ;
  Add2HitsList(h0, 0, !expert, image) ;
  
  TH1F *h1 = new TH1F("hCpvHitsEdep","Hits energy distribution CPV(PMD);Energy [keV];Counts", 500, 0., 500.); 
  h1->Sumw2() ;
  Add2HitsList(h1, 1, !expert, image) ;
  
  TH1I *h2 = new TH1I("hPreHitsMult","Hits multiplicity distribution in PRE(PMD);# of Hits;Entries", 500, 0, 3000) ; 
  h2->Sumw2() ;
  Add2HitsList(h2, 2, !expert, image) ;
  
  TH1I *h3 = new TH1I("hCpvHitsMult","Hits multiplicity distribution in PRE(PMD);# of Hits;Entries", 500, 0, 3000) ; 
  h2->Sumw2() ;
  Add2HitsList(h3, 3, !expert, image) ;
  //
  ClonePerTrigClass(AliQAv1::kHITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliPMDQADataMakerSim::InitSDigits()
{
    // create SDigits histograms in SDigits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *h0 = new TH1F("hPreSDigitsEdep","SDigits energy distribution in(keV) PRE(PMD);Energy [keV];Counts", 500, 0., 500.);
  h0->Sumw2();
  Add2SDigitsList(h0, 0, !expert, image);
  
  TH1F *h1 = new TH1F("hCpvSDigitsEdep","SDigits energy distribution in (keV)CPV(PMD);Energy [keV];Counts", 500, 0., 500.);
  h1->Sumw2();
  Add2SDigitsList(h1, 1, !expert, image);
  
  TH1I *h2 = new TH1I("hPreSDigitsMult","SDigits multiplicity distribution in PRE(PMD);# of SDigits;Entries", 500, 0., 1000.);
  h2->Sumw2();
  Add2SDigitsList(h2, 2, !expert, image);
  
  TH1I *h3 = new TH1I("hCpvSDigitsMult","SDigits multiplicity distribution in CPV(PMD);# of SDigits;Entries", 500, 0., 1000.);
  h3->Sumw2();
  Add2SDigitsList(h3, 3, !expert, image);
  //
  ClonePerTrigClass(AliQAv1::kSDIGITS); // this should be the last line  
}

//____________________________________________________________________________
void AliPMDQADataMakerSim::InitDigits()
{
    // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *h0 = new TH1F("hPreDigitsEdep","Digits energy distribution in PRE(PMD);Amplitude [ADC counts];Counts", 100, 0., 2000.);
  h0->Sumw2();
  Add2DigitsList(h0, 0, !expert, image);
  
  TH1F *h1 = new TH1F("hCpvDigitsEdep","Digits energy distribution in CPV(PMD);Amplitude [ADC counts];Counts", 100, 0., 2000.); 
  h1->Sumw2();
  Add2DigitsList(h1, 1, !expert, image);

  TH1I *h2 = new TH1I("hPreDigitsMult","Digits multiplicity distribution in PRE(PMD);# of Digits;Entries", 500, 0, 1000) ; 
  h2->Sumw2();
  Add2DigitsList(h2, 2, !expert, image);
  
  TH1I *h3 = new TH1I("hCpvDigitsMult","Digits multiplicity distribution in CPV(PMD);# of Digits;Entries", 500, 0, 1000);
  h3->Sumw2();
  Add2DigitsList(h3, 3, !expert, image);
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line  
}

//____________________________________________________________________________ 
void AliPMDQADataMakerSim::MakeHits()
{
    //make QA data from Hits

  Int_t premul = 0, cpvmul = 0;
  Float_t edepkev = 0.;
  TIter next(fHitsArray); 
  AliPMDhit * hit; 
    
  while ( (hit = dynamic_cast<AliPMDhit *>(next())) )
    {
      if (hit->Z() > 361.5)
	{
	  edepkev = (hit->GetEnergy())/1000.;
	  FillHitsData(0,edepkev);
	  premul++;
	}
      else if (hit->Z() < 361.5)
	{
	  edepkev = (hit->GetEnergy())/1000.;
	  FillHitsData(1,edepkev);
	  cpvmul++;
	}
    }
  
  if(premul <= 0)
    {
      FillHitsData(2,-1.); 
    }
  else
    {
      FillHitsData(2,premul); 
    }
  
  if(cpvmul <= 0)
    {
      FillHitsData(3,-1.); 
    }
  else
    {
      FillHitsData(3,cpvmul); 
    }
  
}

//____________________________________________________________________________
void AliPMDQADataMakerSim::MakeHits(TTree * hitTree)
{
  // make QA data from Hit Tree
  
  TBranch * branch = hitTree->GetBranch("PMD") ;
  if ( ! branch )
    {
      AliWarning("PMD branch in Hit Tree not found") ;
      return;
    }

  if (fHitsArray) 
    fHitsArray->Clear() ; 
  else 
    fHitsArray = new TClonesArray("AliPMDhit", 1000);

  branch->SetAddress(&fHitsArray);

  for (Int_t ientry = 0 ; ientry < branch->GetEntries() ; ientry++) {
    branch->GetEntry(ientry) ; 
    MakeHits();
    fHitsArray->Clear() ; 
  } 	
  //
  IncEvCountCycleHits();
  IncEvCountTotalHits();
  //
}
//____________________________________________________________________________
void AliPMDQADataMakerSim::MakeSDigits()
{
    // makes data from SDigits

  Int_t cpvmul = 0, premul = 0;
  Float_t edepkev = 0.;
  
  TIter next(fSDigitsArray) ; 
  AliPMDsdigit * sdigit ; 
  while ( (sdigit = dynamic_cast<AliPMDsdigit *>(next())) )
    {
      if(sdigit->GetDetector() == 0)
	{
	  edepkev = (sdigit->GetCellEdep())/1000.;
	  FillSDigitsData(0,edepkev);
	  premul++;
	}
      if(sdigit->GetDetector() == 1)
	{
	  edepkev = (sdigit->GetCellEdep())/1000.;
	  FillSDigitsData(1,edepkev);
	  cpvmul++;
	}
	
    } 
  if (premul > 0) FillSDigitsData(2,premul);
  if (cpvmul > 0) FillSDigitsData(3,cpvmul);
  
}

//____________________________________________________________________________
void AliPMDQADataMakerSim::MakeSDigits(TTree * sdigitTree)
{
    // makes data from SDigit Tree
  
  if (fSDigitsArray) 
    fSDigitsArray->Clear() ; 
  else 
    fSDigitsArray = new TClonesArray("AliPMDsdigit", 1000) ; 
    
  TBranch * branch = sdigitTree->GetBranch("PMDSDigit") ;
  if ( ! branch )
    {
      AliWarning("PMD branch in SDigit Tree not found") ; 
    }
  else
    {
      branch->SetAddress(&fSDigitsArray) ;
      branch->GetEntry(0) ;
      MakeSDigits() ; 
    }
  //
  IncEvCountCycleSDigits();
  IncEvCountTotalSDigits();
  //
}

//____________________________________________________________________________
void AliPMDQADataMakerSim::MakeDigits()
{
    // makes data from Digits

  Int_t cpvmul = 0, premul = 0;
  
  TIter next(fDigitsArray) ; 
  AliPMDdigit * digit ; 
  while ( (digit = dynamic_cast<AliPMDdigit *>(next())) )
    {
      if(digit->GetDetector() == 0)
	{
	  FillDigitsData(0, digit->GetADC()) ;
	  premul++;
	}
      if(digit->GetDetector() == 1)
	{
	  FillDigitsData(1, digit->GetADC());
	  cpvmul++;
	}
    }  
  
  if (premul > 0) FillDigitsData(2,premul);
  if (cpvmul > 0) FillDigitsData(3,cpvmul);
  
  
}

//____________________________________________________________________________
void AliPMDQADataMakerSim::MakeDigits(TTree * digitTree)
{
    // makes data from Digit Tree

  if (fDigitsArray) 
    fDigitsArray->Clear() ; 
  else 
    fDigitsArray = new TClonesArray("AliPMDdigit", 1000) ; 
  
  TBranch * branch = digitTree->GetBranch("PMDDigit") ;

  if ( ! branch )
    {
      AliWarning("PMD branch in Digit Tree not found") ; 
    }
  else
    {
      branch->SetAddress(&fDigitsArray) ;
      for (Int_t ient = 0; ient < branch->GetEntries(); ient++)
	{
	  branch->GetEntry(ient) ; 
	  MakeDigits() ; 
	  fDigitsArray->Clear() ; 
	  
	}
      
    }
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
}


//____________________________________________________________________________ 

void AliPMDQADataMakerSim::StartOfDetectorCycle()
{
    //Detector specific actions at start of cycle
    
}
//____________________________________________________________________________ 

void AliPMDQADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
    //Detector specific actions at end of cycle
    // do the QA checking
    AliQAChecker::Instance()->Run(AliQAv1::kPMD, task, list) ;  
}
