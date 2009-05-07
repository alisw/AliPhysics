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

    TH1F *h0 = new TH1F("hPreHitsEdep","Hits energy distribution in (keV)PRE(PMD)", 500, 0., 500.); 
    h0->Sumw2() ;
    Add2HitsList(h0, 0) ;

    TH1F *h1 = new TH1F("hCpvHitsEdep","Hits energy distribution in (keV)CPV(PMD)", 500, 0., 500.); 
    h1->Sumw2() ;
    Add2HitsList(h1, 1) ;

    TH1I *h2 = new TH1I("hPreHitsMult","Hits multiplicity distribution in PRE(PMD)", 500, 0, 3000) ; 
    h2->Sumw2() ;
    Add2HitsList(h2, 2) ;

    TH1I *h3 = new TH1I("hCpvHitsMult","Hits multiplicity distribution in PRE(PMD)", 500, 0, 3000) ; 
    h2->Sumw2() ;
    Add2HitsList(h3, 3) ;
}

//____________________________________________________________________________ 
void AliPMDQADataMakerSim::InitSDigits()
{
    // create SDigits histograms in SDigits subdir

    TH1F *h0 = new TH1F("hPreSDigitsEdep","SDigits energy distribution in(keV) PRE(PMD)", 500, 0., 500.);
    h0->Sumw2();
    Add2SDigitsList(h0, 0);

    TH1F *h1 = new TH1F("hCpvSDigitsEdep","SDigits energy distribution in (keV)CPV(PMD)", 500, 0., 500.);
    h1->Sumw2();
    Add2SDigitsList(h1, 1);

    TH1I *h2 = new TH1I("hPreSDigitsMult","SDigits multiplicity distribution in PRE(PMD)", 500, 0., 1000.);
    h2->Sumw2();
    Add2SDigitsList(h2, 2);

    TH1I *h3 = new TH1I("hCpvSDigitsMult","SDigits multiplicity distribution in CPV(PMD)", 500, 0., 1000.);
    h3->Sumw2();
    Add2SDigitsList(h3, 3);

}

//____________________________________________________________________________
void AliPMDQADataMakerSim::InitDigits()
{
    // create Digits histograms in Digits subdir

    TH1F *h0 = new TH1F("hPreDigitsEdep","Digits energy distribution in PRE(PMD)", 100, 0., 2000.);
    h0->Sumw2();
    Add2DigitsList(h0, 0);

    TH1F *h1 = new TH1F("hCpvDigitsEdep","Digits energy distribution in CPV(PMD)", 100, 0., 2000.); 
    h1->Sumw2();
    Add2DigitsList(h1, 1);

    TH1I *h2 = new TH1I("hPreDigitsMult","Digits multiplicity distribution in PRE(PMD)", 500, 0, 1000) ; 
    h2->Sumw2();
    Add2DigitsList(h2, 2);

    TH1I *h3 = new TH1I("hCpvDigitsMult","Digits multiplicity distribution in CPV(PMD)", 500, 0, 1000);
    h3->Sumw2();
    Add2DigitsList(h3, 3);

}

//____________________________________________________________________________ 
void AliPMDQADataMakerSim::MakeHits(TClonesArray *hits)
{
    //make QA data from Hits

    Int_t premul = 0, cpvmul = 0;
    Float_t edepkev = 0.;
    TIter next(hits); 
    AliPMDhit * hit; 
    
    while ( (hit = dynamic_cast<AliPMDhit *>(next())) )
      {
	if (hit->Z() > 361.5)
	  {
	    edepkev = (hit->GetEnergy())/1000.;
	    GetHitsData(0)->Fill(edepkev);
	    premul++;
	  }
	else if (hit->Z() < 361.5)
	  {
	    edepkev = (hit->GetEnergy())/1000.;
	    GetHitsData(1)->Fill(edepkev);
	    cpvmul++;
	  }
    }

    if(premul <= 0)
    {
	GetHitsData(2)->Fill(-1.); 
    }
    else
    {
	GetHitsData(2)->Fill(premul); 
    }

    if(cpvmul <= 0)
    {
	GetHitsData(3)->Fill(-1.); 
    }
    else
    {
	GetHitsData(3)->Fill(cpvmul); 
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

    static TClonesArray statichits("AliPMDhit", 1000);
    statichits.Clear();
    TClonesArray *hits = &statichits;
    static TClonesArray staticdummy("AliPMDhit", 1000);
    staticdummy.Clear();
    TClonesArray *dummy = &staticdummy;
    branch->SetAddress(&dummy);
    Int_t index = 0 ;  

    for (Int_t ientry = 0 ; ientry < branch->GetEntries() ; ientry++) {
	branch->GetEntry(ientry) ; 
	for (Int_t ihit = 0 ; ihit < dummy->GetEntries() ; ihit++) {
	    AliPMDhit * hit = dynamic_cast<AliPMDhit *> (dummy->At(ihit)) ; 
	    new((*hits)[index]) AliPMDhit(*hit) ; 

	    index++ ;
	} 
    } 	

    MakeHits(hits);

}
//____________________________________________________________________________
void AliPMDQADataMakerSim::MakeSDigits(TClonesArray * sdigits)
{
    // makes data from SDigits

    Int_t cpvmul = 0, premul = 0;
    Float_t edepkev = 0.;

    TIter next(sdigits) ; 
    AliPMDsdigit * sdigit ; 
    while ( (sdigit = dynamic_cast<AliPMDsdigit *>(next())) )
    {
	if(sdigit->GetDetector() == 0)
	{
	  edepkev = (sdigit->GetCellEdep())/1000.;
	  GetSDigitsData(0)->Fill(edepkev);
	  premul++;
	}
	if(sdigit->GetDetector() == 1)
	{
	  edepkev = (sdigit->GetCellEdep())/1000.;
	  GetSDigitsData(1)->Fill(edepkev);
	  cpvmul++;
	}
	
    } 
    if (premul > 0) GetSDigitsData(2)->Fill(premul);
    if (cpvmul > 0) GetSDigitsData(3)->Fill(cpvmul);
    
}

//____________________________________________________________________________
void AliPMDQADataMakerSim::MakeSDigits(TTree * sdigitTree)
{
    // makes data from SDigit Tree

    TClonesArray * sdigits = new TClonesArray("AliPMDsdigit", 1000) ; 
    
    TBranch * branch = sdigitTree->GetBranch("PMDSDigit") ;
    branch->SetAddress(&sdigits) ;

    if ( ! branch )
    {
	AliWarning("PMD branch in SDigit Tree not found") ; 
    }
    else
    {
	for (Int_t ient = 0; ient < branch->GetEntries(); ient++)
	{
	    branch->GetEntry(ient) ;
	    MakeSDigits(sdigits) ; 
	}
    }
}

//____________________________________________________________________________
void AliPMDQADataMakerSim::MakeDigits(TClonesArray * digits)
{
    // makes data from Digits
    
    Int_t cpvmul = 0, premul = 0;

    TIter next(digits) ; 
    AliPMDdigit * digit ; 
    while ( (digit = dynamic_cast<AliPMDdigit *>(next())) )
    {
	if(digit->GetDetector() == 0)
	{
	    GetDigitsData(0)->Fill( digit->GetADC()) ;
	    premul++;
	}
	if(digit->GetDetector() == 1)
	{
	    GetDigitsData(1)->Fill( digit->GetADC());
	    cpvmul++;
	}
    }  

    if (premul > 0) GetDigitsData(2)->Fill(premul);
    if (cpvmul > 0) GetDigitsData(3)->Fill(cpvmul);


}

//____________________________________________________________________________
void AliPMDQADataMakerSim::MakeDigits(TTree * digitTree)
{
    // makes data from Digit Tree

    TClonesArray * digits = new TClonesArray("AliPMDdigit", 1000) ; 
    
    TBranch * branch = digitTree->GetBranch("PMDDigit") ;
    branch->SetAddress(&digits) ;

    if ( ! branch )
    {
	AliWarning("PMD branch in Digit Tree not found") ; 
    }
    else
    {
	for (Int_t ient = 0; ient < branch->GetEntries(); ient++)
	{
	    branch->GetEntry(ient) ; 
	    MakeDigits(digits) ; 
	}
	
    }
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
