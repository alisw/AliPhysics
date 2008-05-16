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

// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliEMCALDigit.h"
#include "AliEMCALHit.h"
#include "AliEMCALQADataMaker.h"
#include "AliQAChecker.h"
#include "AliEMCALRecPoint.h" 
#include "AliEMCALRawUtils.h"
#include "AliEMCALReconstructor.h"
#include "AliEMCALRecParam.h"

ClassImp(AliEMCALQADataMaker)
           
//____________________________________________________________________________ 
  AliEMCALQADataMaker::AliEMCALQADataMaker() : 
  AliQADataMaker(AliQA::GetDetName(AliQA::kEMCAL), "EMCAL Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliEMCALQADataMaker::AliEMCALQADataMaker(const AliEMCALQADataMaker& qadm) :
  AliQADataMaker()
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliEMCALQADataMaker& AliEMCALQADataMaker::operator = (const AliEMCALQADataMaker& qadm )
{
  // Equal operator.
  this->~AliEMCALQADataMaker();
  new(this) AliEMCALQADataMaker(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliEMCALQADataMaker::EndOfDetectorCycle(AliQA::TASKINDEX task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kEMCAL, task, list) ;  
}

//____________________________________________________________________________ 
void AliEMCALQADataMaker::InitESDs()
{
  //create ESDs histograms in ESDs subdir
	
  TH1F * h0 = new TH1F("hEmcalESDs",    "ESDs energy distribution in EMCAL",       100, 0., 100.) ;  
  h0->Sumw2() ; 
  Add2ESDsList(h0, 0) ;
  TH1I * h1  = new TH1I("hEmcalESDsMul", "ESDs multiplicity distribution in EMCAL", 100, 0., 100) ; 
  h1->Sumw2() ;
  Add2ESDsList(h1, 1) ;
}

//____________________________________________________________________________ 
void AliEMCALQADataMaker::InitHits()
{
  // create Hits histograms in Hits subdir
  TH1F * h0 = new TH1F("hEmcalHits",    "Hits energy distribution in EMCAL",       100, 0., 100.) ; 
  h0->Sumw2() ;
  Add2HitsList(h0, 0) ;
  TH1I * h1  = new TH1I("hEmcalHitsMul", "Hits multiplicity distribution in EMCAL", 500, 0., 10000) ; 
  h1->Sumw2() ;
  Add2HitsList(h1, 1) ;
}

//____________________________________________________________________________ 
void AliEMCALQADataMaker::InitDigits()
{
  // create Digits histograms in Digits subdir
  TH1I * h0 = new TH1I("hEmcalDigits",    "Digits amplitude distribution in EMCAL",    500, 0, 5000) ; 
  h0->Sumw2() ;
  Add2DigitsList(h0, 0) ;
  TH1I * h1 = new TH1I("hEmcalDigitsMul", "Digits multiplicity distribution in EMCAL", 500, 0, 1000) ; 
  h1->Sumw2() ;
  Add2DigitsList(h1, 1) ;
}

//____________________________________________________________________________ 
void AliEMCALQADataMaker::InitRecPoints()
{
  // create Reconstructed Points histograms in RecPoints subdir
  TH1F * h0 = new TH1F("hEmcalRecPoints",    "RecPoints energy distribution in EMCAL",       100, 0., 100.) ; 
  h0->Sumw2() ;
  Add2RecPointsList(h0, 0) ;
  TH1I * h1 = new TH1I("hEmcalRecPointsMul", "RecPoints multiplicity distribution in EMCAL", 100, 0,  100) ; 
  h1->Sumw2() ;
  Add2RecPointsList(h1, 1) ;

}

//____________________________________________________________________________ 
void AliEMCALQADataMaker::InitRaws()
{
  AliInfo(Form("Raw QA infor for EMCAL not yet implemented"));
}

//____________________________________________________________________________ 
void AliEMCALQADataMaker::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
  TH1F * h0 = new TH1F("hEmcalSDigits",    "SDigits energy distribution in EMCAL",       100, 0., 100.) ; 
  h0->Sumw2() ;
  Add2SDigitsList(h0, 0) ;
  TH1I * h1 = new TH1I("hEmcalSDigitsMul", "SDigits multiplicity distribution in EMCAL", 500, 0,  10000) ; 
  h1->Sumw2() ;
  Add2SDigitsList(h1, 1) ;
}

//____________________________________________________________________________
void AliEMCALQADataMaker::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs

  Int_t count = 0 ; 
  for ( Int_t index = 0; index < esd->GetNumberOfCaloClusters() ; index++ ) {
	AliESDCaloCluster * clu = esd->GetCaloCluster(index) ;
	if ( clu->IsEMCAL() ) {
		GetESDsData(0)->Fill(clu->E()) ;
		count++ ;
	} 
  }
  GetESDsData(1)->Fill(count) ;
}

//____________________________________________________________________________
void AliEMCALQADataMaker::MakeHits(TClonesArray * hits)
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
void AliEMCALQADataMaker::MakeHits(TTree * hitTree)
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
void AliEMCALQADataMaker::MakeDigits(TClonesArray * digits)
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
void AliEMCALQADataMaker::MakeDigits(TTree * digitTree)
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
void AliEMCALQADataMaker::MakeRaws(AliRawReader* rawReader)
{
  //Raw QA info not yet implemented for EMCAL

}

//____________________________________________________________________________
void AliEMCALQADataMaker::MakeRecPoints(TTree * clustersTree)
{
  // makes data from RecPoints
  TBranch *emcbranch = clustersTree->GetBranch("EMCALECARP");
  if (!emcbranch) { 
    AliError("can't get the branch with the EMCAL clusters !");
    return;
  }
  TObjArray * emcrecpoints = new TObjArray(100) ;
  emcbranch->SetAddress(&emcrecpoints);
  emcbranch->GetEntry(0);
  
  GetRecPointsData(1)->Fill(emcrecpoints->GetEntriesFast()) ; 
  TIter next(emcrecpoints) ; 
  AliEMCALRecPoint * rp ; 
  while ( (rp = dynamic_cast<AliEMCALRecPoint *>(next())) ) {
    GetRecPointsData(0)->Fill( rp->GetEnergy()) ;
  }
  emcrecpoints->Delete();
  delete emcrecpoints;

}

//____________________________________________________________________________
void AliEMCALQADataMaker::MakeSDigits(TClonesArray * sdigits)
{
  // makes data from SDigits
  
  GetSDigitsData(1)->Fill(sdigits->GetEntriesFast()) ; 
  TIter next(sdigits) ; 
  AliEMCALDigit * sdigit ; 
  while ( (sdigit = dynamic_cast<AliEMCALDigit *>(next())) ) {
    GetSDigitsData(0)->Fill( sdigit->GetAmp()) ;
  } 
}

//____________________________________________________________________________
void AliEMCALQADataMaker::MakeSDigits(TTree * sdigitTree)
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
void AliEMCALQADataMaker::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}
