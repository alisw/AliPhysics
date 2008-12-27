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

// --- Standard library ---
#include <Riostream.h>
// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h>
#include <TBranch.h>
#include <TTree.h>
// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAChecker.h"
#include "AliZDCQADataMakerSim.h"
#include "AliZDCHit.h"
#include "AliZDCDigit.h"
#include "AliZDCRawStream.h"

ClassImp(AliZDCQADataMakerSim)
           
//____________________________________________________________________________ 
  AliZDCQADataMakerSim::AliZDCQADataMakerSim() : 
      AliQADataMakerSim(AliQA::GetDetName(AliQA::kZDC), "ZDC Quality Assurance Data Maker"),
      fHits(0),
      fDigit(0)
{
  // ctor
}

//____________________________________________________________________________ 
AliZDCQADataMakerSim::AliZDCQADataMakerSim(const AliZDCQADataMakerSim& qadm) :
    AliQADataMakerSim(), 
    fHits(0),
    fDigit(0) 
{
  //copy ctor 
  SetName((const char*)qadm.GetName()); 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliZDCQADataMakerSim& AliZDCQADataMakerSim::operator = (const AliZDCQADataMakerSim& qadm )
{
  // Equal operator.
  this->~AliZDCQADataMakerSim();
  new(this) AliZDCQADataMakerSim(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliZDCQADataMakerSim::InitHits()
{
  // create Hits histograms in Hits subdir
  //
  TH2F * hZNCh  = new TH2F("hZNCh", "Hits centroid in ZNC", 100, -5.,5.,100,-5.,5.);
  TH2F * hZNAh  = new TH2F("hZNAh", "Hits centroid in ZNA", 100, -5.,5.,100,-5.,5.);
  // NB -> For the moment no check is performesd on ZP centroids
//  TH2F * hZPCh  = new TH2F("hZPCh", "Hits centroid in ZPC", 100,-12.,12.,100,-12.,12.); 
//  TH2F * hZPAh  = new TH2F("hZPAh", "Hits centroid in ZPA", 100,-12.,12.,100,-12.,12.); 
  Add2HitsList(hZNCh, 0);
  Add2HitsList(hZNAh, 1);
//  Add2HitsList(hZPCh, 2);
//  Add2HitsList(hZPAh, 3);
}


//____________________________________________________________________________ 
void AliZDCQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  //
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hDigZNCTot = new TH1F("hDigZNCTot", "Signal in ZNC", 100, 0., 6000.);
  TH1F * hDigZNATot = new TH1F("hDigZNATot", "Signal in ZNA", 100, 0., 6000.);
  TH1F * hDigZPCTot = new TH1F("hDigZPCTot", "Signal in ZPC", 100, 0., 6000.);
  TH1F * hDigZPATot = new TH1F("hDigZPATot", "Signal in ZPA", 100, 0., 6000.);
  Add2DigitsList(hDigZNCTot, 0);
  Add2DigitsList(hDigZNATot, 1);
  Add2DigitsList(hDigZPCTot, 2);
  Add2DigitsList(hDigZPATot, 3);
  //
  TH1F * hDigSumQZNC = new TH1F("hDigSumQZNC", "Signal in 4 ZNC PMQ",100, 0., 4000.);
  TH1F * hDigSumQZNA = new TH1F("hDigSumQZNA", "Signal in 4 ZNA PMQ",100, 0., 4000.);
  TH1F * hDigSumQZPC = new TH1F("hDigSumQZPC", "Signal in 4 ZPC PMQ",100, 0., 4000.);
  TH1F * hDigSumQZPA = new TH1F("hDigSumQZPA", "Signal in 4 ZPA PMQ",100, 0., 4000.);
  Add2DigitsList(hDigSumQZNC, 4, kTRUE);
  Add2DigitsList(hDigSumQZNA, 5, kTRUE);
  Add2DigitsList(hDigSumQZPC, 6, kTRUE);
  Add2DigitsList(hDigSumQZPA, 7, kTRUE);
  //
  TH1F * hDigPMCZNC = new TH1F("hDigPMCZNC", "Signal in ZNC PMC",100, 0., 4000.);
  TH1F * hDigPMCZNA = new TH1F("hDigPMCZNA", "Signal in ZNA PMC",100, 0., 4000.);
  TH1F * hDigPMCZPC = new TH1F("hDigPMCZPC", "Signal in ZPC PMC",100, 0., 4000.);
  TH1F * hDigPMCZPA = new TH1F("hDigPMCZPA", "Signal in ZPA PMC",100, 0., 4000.);
  Add2DigitsList(hDigPMCZNC, 8, kTRUE);
  Add2DigitsList(hDigPMCZNA, 9, kTRUE);
  Add2DigitsList(hDigPMCZPC, 10, kTRUE);
  Add2DigitsList(hDigPMCZPA, 11, kTRUE);
  // 
  // ------------------- LOW GAIN CHAIN ---------------------------
/*  TH1F * hDigZNCTotlg = new TH1F("hDigZNCTotlg", "Digit lg signal in ZNC", 100, 0., 6000.);
  TH1F * hDigZNATotlg = new TH1F("hDigZNATotlg", "Digit lg signal in ZNA", 100, 0., 6000.);
  TH1F * hDigZPCTotlg = new TH1F("hDigZPCTotlg", "Digit lg signal in ZPC", 100, 0., 6000.);
  TH1F * hDigZPATotlg = new TH1F("hDigZPATotlg", "Digit lg signal in ZPA", 100, 0., 6000.);
  Add2DigitsList(hDigZNCTotlg, 12);
  Add2DigitsList(hDigZNATotlg, 13);
  Add2DigitsList(hDigZPCTotlg, 14);
  Add2DigitsList(hDigZPATotlg, 15);
  //
  TH1F * hDigSumQZNClg = new TH1F("hDigSumQZNClg", "Signal in 4 ZNC PMQlg",100, 0., 4000.);
  TH1F * hDigSumQZNAlg = new TH1F("hDigSumQZNAlg", "Signal in 4 ZNA PMQlg",100, 0., 4000.);
  TH1F * hDigSumQZPClg = new TH1F("hDigSumQZPClg", "Signal in 4 ZPC PMQlg",100, 0., 4000.);
  TH1F * hDigSumQZPAlg = new TH1F("hDigSumQZPAlg", "Signal in 4 ZPA PMQlg",100, 0., 4000.);
  Add2DigitsList(hDigSumQZNClg, 16, kTRUE);
  Add2DigitsList(hDigSumQZNAlg, 17, kTRUE);
  Add2DigitsList(hDigSumQZPClg, 18, kTRUE);
  Add2DigitsList(hDigSumQZPAlg, 19, kTRUE);
  //
  TH1F * hDigPMCZNClg = new TH1F("hDigPMCZNClg", "Signal in ZNC PMClg",100, 0., 4000.);
  TH1F * hDigPMCZNAlg = new TH1F("hDigPMCZNAlg", "Signal in ZNA PMClg",100, 0., 4000.);
  TH1F * hDigPMCZPClg = new TH1F("hDigPMCZPClg", "Signal in ZPC PMClg",100, 0., 4000.);
  TH1F * hDigPMCZPAlg = new TH1F("hDigPMCZPAlg", "Signal in ZPA PMClg",100, 0., 4000.);
  Add2DigitsList(hDigPMCZNClg, 20, kTRUE);
  Add2DigitsList(hDigPMCZNAlg, 21, kTRUE);
  Add2DigitsList(hDigPMCZPClg, 22, kTRUE);
  Add2DigitsList(hDigPMCZPAlg, 23, kTRUE);
*/
}

//____________________________________________________________________________
void AliZDCQADataMakerSim::MakeHits(TClonesArray * /*data*/)
{
  //filling QA histos for Hits
  //
    TIter next(fHits); 
    AliZDCHit * hit; 
    while((hit = dynamic_cast<AliZDCHit *>(next()))){
      if(hit->GetVolume(0)==1) GetHitsData(0)->Fill(hit->GetXImpact(),hit->GetYImpact());
      else if(hit->GetVolume(0)==4) GetHitsData(1)->Fill(hit->GetXImpact(), hit->GetYImpact());
//      else if(hit->GetVolume(0)==2) GetHitsData(1)->Fill(hit->GetXImpact(), hit->GetYImpact());
//      else if(hit->GetVolume(0)==5) GetHitsData(3)->Fill(hit->GetXImpact(), hit->GetYImpact());
    }
}

//___________________________________________________________________________
void AliZDCQADataMakerSim::MakeHits(TTree * hitTree)
{
  // make QA data from Hit Tree
  //
  if(!hitTree){
    AliError("Hit Tree not found!"); 
    return;
  }
  //

  TBranch * branch = hitTree->GetBranch("ZDC") ;

  if(!branch){
    AliError("ZDC branch in Hit Tree not found!"); 
    return;
  } 
  else{
    char** add = (char**) (branch->GetAddress());
    if(add){
        fHits = (TClonesArray*)(*add);
    } 
    else{
        if(!fHits) fHits = new TClonesArray("AliZDCHit", 1000);
        branch->SetAddress(&fHits);
    }
    Int_t ntracks = (Int_t) hitTree->GetEntries();
    //printf("\n\t *** no.track %d\n",ntracks);
    if (ntracks<=0) return;
    //
    for(Int_t itrack=0; itrack<ntracks; itrack++){
        
        branch->GetEntry(itrack);
        //
        //printf("\t *** track %d",itrack);
        //hits->Print("");
        //printf("\n");
        //
        MakeHits(); 
        fHits->Clear();
    }	
  }
}

//___________________________________________________________________________
void AliZDCQADataMakerSim::MakeDigits(TTree *digitTree )
{
  // makes data from Digit Tree
  TBranch * branch = digitTree->GetBranch("ZDC");
  if(!branch){
    AliError("ZDC branch in Digit Tree not found"); 
    return;
  } 
  char** add = (char**) (branch->GetAddress());
  if(add){
      fDigit = (AliZDCDigit*)(*add);
  } 
  else{
      if(!fDigit) fDigit = new AliZDCDigit();
      branch->SetAddress(&fDigit);
  }
  
  Int_t ndig = digitTree->GetEntries();
   
  Float_t adcSum_ZNC=0., adcSum_ZNA=0., adcSum_ZPC=0., adcSum_ZPA=0.;
  Float_t adcSumQ_ZNC=0., adcSumQ_ZNA=0., adcSumQ_ZPC=0., adcSumQ_ZPA=0.;
  //Float_t adcSum_ZNC_lg=0., adcSum_ZNA_lg=0., adcSum_ZPC_lg=0., adcSum_ZPA_lg=0.;
  //Float_t adcSumQ_ZNC_lg=0., adcSumQ_ZNA_lg=0., adcSumQ_ZPC_lg=0., adcSumQ_ZPA_lg=0.;
  //
  for(Int_t i = 0; i < ndig; i++){
      digitTree->GetEntry(i);
      if(fDigit->GetSector(0)==1){
	  adcSum_ZNC += fDigit->GetADCValue(0);
	  //adcSum_ZNC_lg += fDigit->GetADCValue(1);
	  //
	  if(fDigit->GetSector(1)!=0){
	      adcSumQ_ZNC += fDigit->GetADCValue(0);
	      //adcSumQ_ZNC_lg+= fDigit->GetADCValue(1);
	  }
	  else{
	      GetDigitsData(8)->Fill(fDigit->GetADCValue(0));
	      //GetDigitsData(20)->Fill(fDigit->GetADCValue(1));
	  }
      }
      else if(fDigit->GetSector(0)==2){
	  adcSum_ZPC += fDigit->GetADCValue(0);
	  //adcSum_ZPC_lg += fDigit->GetADCValue(1);
	  //
	  if(fDigit->GetSector(1)!=0){
	      adcSumQ_ZPC += fDigit->GetADCValue(0);
	      //adcSumQ_ZPC_lg+= fDigit->GetADCValue(1);
	  }
	  else{
	      GetDigitsData(10)->Fill(fDigit->GetADCValue(0));
	      //GetDigitsData(22)->Fill(fDigit->GetADCValue(1));
	  }
      }
      else if(fDigit->GetSector(0)==4){
	  adcSum_ZNA += fDigit->GetADCValue(0);
	  //adcSum_ZNA_lg += fDigit->GetADCValue(1);
	  //
	  if(fDigit->GetSector(1)!=0){
	      adcSumQ_ZNA += fDigit->GetADCValue(0);
	      //adcSumQ_ZNA_lg+= fDigit->GetADCValue(1);
	  }
	  else{
	      GetDigitsData(9)->Fill(fDigit->GetADCValue(0));
	      //GetDigitsData(21)->Fill(fDigit->GetADCValue(1));
	  }
      }
      else if(fDigit->GetSector(0)==5){
	  adcSum_ZPA += fDigit->GetADCValue(0);
	  //adcSum_ZPA_lg += fDigit->GetADCValue(1);
	  //
	  if(fDigit->GetSector(1)!=0){
	      adcSumQ_ZPA += fDigit->GetADCValue(0);
	      //adcSumQ_ZPA_lg+= fDigit->GetADCValue(1);
	  }
	  else{
	      GetDigitsData(11)->Fill(fDigit->GetADCValue(0));
	      //GetDigitsData(23)->Fill(fDigit->GetADCValue(1));
	  }
      }
  }
  //
  GetDigitsData(0)->Fill(adcSum_ZNC);
  GetDigitsData(1)->Fill(adcSum_ZNA);
  GetDigitsData(2)->Fill(adcSum_ZPC);
  GetDigitsData(3)->Fill(adcSum_ZPA);
  //
  GetDigitsData(4)->Fill(adcSumQ_ZNC);
  GetDigitsData(5)->Fill(adcSumQ_ZNA);
  GetDigitsData(6)->Fill(adcSumQ_ZPC);
  GetDigitsData(7)->Fill(adcSumQ_ZPA);
  //
  /*GetDigitsData(12)->Fill(adcSum_ZNC_lg);
  GetDigitsData(13)->Fill(adcSum_ZNA_lg);
  GetDigitsData(14)->Fill(adcSum_ZPC_lg);
  GetDigitsData(15)->Fill(adcSum_ZPA_lg);
  //
  GetDigitsData(16)->Fill(adcSumQ_ZNC_lg);
  GetDigitsData(17)->Fill(adcSumQ_ZNA_lg);
  GetDigitsData(18)->Fill(adcSumQ_ZPC_lg);
  GetDigitsData(19)->Fill(adcSumQ_ZPA_lg);*/
}

//____________________________________________________________________________
void AliZDCQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

//____________________________________________________________________________ 
void AliZDCQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray ** list)
{
  // Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kZDC, task, list);  
}
