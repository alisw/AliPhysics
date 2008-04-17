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
  AliQADataMakerSim(AliQA::GetDetName(AliQA::kZDC), "ZDC Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliZDCQADataMakerSim::AliZDCQADataMakerSim(const AliZDCQADataMakerSim& qadm) :
  AliQADataMakerSim() 
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
  TH2F * hZPCh  = new TH2F("hZPCh", "Hits centroid in ZPC", 100,-12.,12.,100,-12.,12.); 
  TH2F * hZPAh  = new TH2F("hZPAh", "Hits centroid in ZPA", 100,-12.,12.,100,-12.,12.); 
  Add2HitsList(hZNCh, 0);
  Add2HitsList(hZPCh, 1);
  Add2HitsList(hZNAh, 2);
  Add2HitsList(hZPAh, 3);
}


//____________________________________________________________________________ 
void AliZDCQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  //
  TH1F * hDigZNCTot = new TH1F("hDigZNCTot", "Digit signal in ZNC", 100, 0., 6000.);
  TH1F * hDigZNATot = new TH1F("hDigZNATot", "Digit signal in ZNA", 100, 0., 6000.);
  TH1F * hDigZPCTot = new TH1F("hDigZPCTot", "Digit signal in ZPC", 100, 0., 10000.);
  TH1F * hDigZPATot = new TH1F("hDigZPATot", "Digit signal in ZPA", 100, 0., 10000.);
  Add2DigitsList(hDigZNCTot, 0);
  Add2DigitsList(hDigZPCTot, 1);
  Add2DigitsList(hDigZNATot, 2);
  Add2DigitsList(hDigZPATot, 3);
  //
  TH1F * hDigSumQZNC = new TH1F("hDigSumQZNC", "Signal in 4 ZNC PMQ[i]",100, 0., 4000.);
  TH1F * hDigSumQZPC = new TH1F("hDigSumQZPC", "Signal in 4 ZPC PMQ[i]",100, 0., 4000.);
  TH1F * hDigSumQZNA = new TH1F("hDigSumQZNA", "Signal in 4 ZNA PMQ[i]",100, 0., 4000.);
  TH1F * hDigSumQZPA = new TH1F("hDigSumQZPA", "Signal in 4 ZPA PMQ[i]",100, 0., 4000.);
  Add2DigitsList(hDigSumQZNC, 4);
  Add2DigitsList(hDigSumQZPC, 5);
  Add2DigitsList(hDigSumQZNA, 6);
  Add2DigitsList(hDigSumQZPA, 7);
  //
  TH1F * hDigPMCZNC = new TH1F("hDigPMCZNC", "Signal in 4 ZNC PMQ[i]",100, 0., 4000.);
  TH1F * hDigPMCZPC = new TH1F("hDigPMCZPC", "Signal in 4 ZPC PMQ[i]",100, 0., 4000.);
  TH1F * hDigPMCZNA = new TH1F("hDigPMCZNA", "Signal in 4 ZNA PMQ[i]",100, 0., 4000.);
  TH1F * hDigPMCZPA = new TH1F("hDigPMCZPA", "Signal in 4 ZPA PMQ[i]",100, 0., 4000.);
  Add2DigitsList(hDigPMCZNC, 8);
  Add2DigitsList(hDigPMCZPC, 9);
  Add2DigitsList(hDigPMCZNA, 10);
  Add2DigitsList(hDigPMCZPA, 11);
  // 
  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hDigZNCTotlg = new TH1F("hDigZNCTotlg", "Digit lg signal in ZNC", 100, 0., 6000.);
  TH1F * hDigZNATotlg = new TH1F("hDigZNATotlg", "Digit lg signal in ZNA", 100, 0., 6000.);
  TH1F * hDigZPCTotlg = new TH1F("hDigZPCTotlg", "Digit lg signal in ZPC", 100, 0., 10000.);
  TH1F * hDigZPATotlg = new TH1F("hDigZPATotlg", "Digit lg signal in ZPA", 100, 0., 10000.);
  Add2DigitsList(hDigZNCTotlg, 12);
  Add2DigitsList(hDigZPCTotlg, 13);
  Add2DigitsList(hDigZNATotlg, 14);
  Add2DigitsList(hDigZPATotlg, 15);
  //
  TH1F * hDigSumQZNClg = new TH1F("hDigSumQZNClg", "Signal in 4 ZNC PMQlg[i]",100, 0., 4000.);
  TH1F * hDigSumQZPClg = new TH1F("hDigSumQZPClg", "Signal in 4 ZPC PMQlg[i]",100, 0., 4000.);
  TH1F * hDigSumQZNAlg = new TH1F("hDigSumQZNAlg", "Signal in 4 ZNA PMQlg[i]",100, 0., 4000.);
  TH1F * hDigSumQZPAlg = new TH1F("hDigSumQZPAlg", "Signal in 4 ZPA PMQlg[i]",100, 0., 4000.);
  Add2DigitsList(hDigSumQZNClg, 16);
  Add2DigitsList(hDigSumQZPClg, 17);
  Add2DigitsList(hDigSumQZNAlg, 18);
  Add2DigitsList(hDigSumQZPAlg, 19);
  //
  TH1F * hDigPMCZNClg = new TH1F("hDigPMCZNClg", "Signal in 4 ZNC PMQlg[i]",100, 0., 4000.);
  TH1F * hDigPMCZPClg = new TH1F("hDigPMCZPClg", "Signal in 4 ZPC PMQlg[i]",100, 0., 4000.);
  TH1F * hDigPMCZNAlg = new TH1F("hDigPMCZNAlg", "Signal in 4 ZNA PMQlg[i]",100, 0., 4000.);
  TH1F * hDigPMCZPAlg = new TH1F("hDigPMCZPAlg", "Signal in 4 ZPA PMQlg[i]",100, 0., 4000.);
  Add2DigitsList(hDigPMCZNClg, 20);
  Add2DigitsList(hDigPMCZPClg, 21);
  Add2DigitsList(hDigPMCZNAlg, 22);
  Add2DigitsList(hDigPMCZPAlg, 23);
}

//____________________________________________________________________________
void AliZDCQADataMakerSim::MakeHits(TClonesArray * data)
{
  //filling QA histos for Hits
  //
  TClonesArray * hits = dynamic_cast<TClonesArray *>(data); 
  if(!hits){
    AliError("Wrong type of hits container"); 
  } 
  else {
    TIter next(hits); 
    AliZDCHit * hit; 
    while((hit = dynamic_cast<AliZDCHit *>(next()))){
      if(hit->GetVolume(0)==1) GetHitsData(0)->Fill(hit->GetXImpact(),hit->GetYImpact());
      else if(hit->GetVolume(0)==2) GetHitsData(1)->Fill(hit->GetXImpact(), hit->GetYImpact());
      else if(hit->GetVolume(0)==4) GetHitsData(2)->Fill(hit->GetXImpact(), hit->GetYImpact());
      else if(hit->GetVolume(0)==5) GetHitsData(3)->Fill(hit->GetXImpact(), hit->GetYImpact());
    }
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
  TClonesArray * hits = new TClonesArray("AliZDCHit", 1000);
  TBranch * branch = hitTree->GetBranch("ZDC") ;
  if(!branch){
    AliError("ZDC branch in Hit Tree not found!"); 
    return;
  }
  //
  TClonesArray * tmp = new TClonesArray("AliZDCHit", 1000);
  branch->SetAddress(&tmp) ;
  Int_t index = 0 ;  
  for(Int_t ientry = 0; ientry<branch->GetEntries(); ientry++) {
      branch->GetEntry(ientry) ; 
      for(Int_t ihit =0 ; ihit<tmp->GetEntries(); ihit++){
  	 AliZDCHit * hit = dynamic_cast<AliZDCHit *> (tmp->At(ihit)); 
  	 new((*hits)[index]) AliZDCHit(*hit); 
  	 index++;
      } 
  }	  
  tmp->Delete(); 
  delete tmp; 
  MakeHits(hits); 
}

//____________________________________________________________________________
void AliZDCQADataMakerSim::MakeDigits(TClonesArray * digits)
{
  // makes data from Digits
  //
  TIter next(digits) ; 
  AliZDCDigit * digit ; 
  //
  Float_t ADCSum_ZNC=0., ADCSum_ZNA=0., ADCSum_ZPC=0., ADCSum_ZPA=0.;
  Float_t ADCSumQ_ZNC=0., ADCSumQ_ZNA=0., ADCSumQ_ZPC=0., ADCSumQ_ZPA=0.;
  Float_t ADCSum_ZNC_lg=0., ADCSum_ZNA_lg=0., ADCSum_ZPC_lg=0., ADCSum_ZPA_lg=0.;
  Float_t ADCSumQ_ZNC_lg=0., ADCSumQ_ZNA_lg=0., ADCSumQ_ZPC_lg=0., ADCSumQ_ZPA_lg=0.;
  //
  while((digit = dynamic_cast<AliZDCDigit *>(next()))){
    if(digit->GetSector(0)==1){
      ADCSum_ZNC += digit->GetADCValue(0);
      ADCSum_ZNC_lg += digit->GetADCValue(1);
      //
      if(digit->GetSector(1)!=0){
        ADCSumQ_ZNC += digit->GetADCValue(0);
        ADCSumQ_ZNC_lg+= digit->GetADCValue(1);
      }
      else{
        GetDigitsData(8)->Fill(digit->GetADCValue(0));
        GetDigitsData(20)->Fill(digit->GetADCValue(1));
      }
    }
    else if(digit->GetSector(0)==2){
      ADCSum_ZPC += digit->GetADCValue(0);
      ADCSum_ZPC_lg += digit->GetADCValue(1);
      //
      if(digit->GetSector(1)!=0){
        ADCSumQ_ZPC += digit->GetADCValue(0);
        ADCSumQ_ZPC_lg+= digit->GetADCValue(1);
      }
      else{
        GetDigitsData(9)->Fill(digit->GetADCValue(0));
        GetDigitsData(21)->Fill(digit->GetADCValue(1));
      }
    }
    else if(digit->GetSector(0)==4){
      ADCSum_ZNA += digit->GetADCValue(0);
      ADCSum_ZNA_lg += digit->GetADCValue(1);
      //
      if(digit->GetSector(1)!=0){
        ADCSumQ_ZNA += digit->GetADCValue(0);
        ADCSumQ_ZNA_lg+= digit->GetADCValue(1);
      }
      else{
        GetDigitsData(10)->Fill(digit->GetADCValue(0));
        GetDigitsData(22)->Fill(digit->GetADCValue(1));
      }
    }
    else if(digit->GetSector(0)==5){
      ADCSum_ZPA += digit->GetADCValue(0);
      ADCSum_ZPA_lg += digit->GetADCValue(1);
      //
      if(digit->GetSector(1)!=0){
        ADCSumQ_ZPA += digit->GetADCValue(0);
        ADCSumQ_ZPA_lg+= digit->GetADCValue(1);
      }
      else{
        GetDigitsData(11)->Fill(digit->GetADCValue(0));
        GetDigitsData(23)->Fill(digit->GetADCValue(1));
      }
    }
  }
  //
  GetDigitsData(0)->Fill(ADCSum_ZNC);
  GetDigitsData(1)->Fill(ADCSum_ZPC);
  GetDigitsData(2)->Fill(ADCSum_ZNA);
  GetDigitsData(3)->Fill(ADCSum_ZPA);
  //
  GetDigitsData(4)->Fill(ADCSumQ_ZNC);
  GetDigitsData(5)->Fill(ADCSumQ_ZPC);
  GetDigitsData(6)->Fill(ADCSumQ_ZNA);
  GetDigitsData(7)->Fill(ADCSumQ_ZPA);
  //
  GetDigitsData(12)->Fill(ADCSum_ZNC_lg);
  GetDigitsData(13)->Fill(ADCSum_ZPC_lg);
  GetDigitsData(14)->Fill(ADCSum_ZNA_lg);
  GetDigitsData(15)->Fill(ADCSum_ZPA_lg);
  //
  GetDigitsData(16)->Fill(ADCSumQ_ZNC_lg);
  GetDigitsData(17)->Fill(ADCSumQ_ZPC_lg);
  GetDigitsData(18)->Fill(ADCSumQ_ZNA_lg);
  GetDigitsData(19)->Fill(ADCSumQ_ZPA_lg);
   
}
//___________________________________________________________________________
void AliZDCQADataMakerSim::MakeDigits(TTree *digitTree )
{
   // makes data from Digit Tree
   TClonesArray * digits = new TClonesArray("AliZDCDigit", 1000); 
   //
   TBranch * branch = digitTree->GetBranch("ZDC");
   if(!branch){
      AliError("ZDC branch in Digit Tree not found"); 
      return;
   } 
   branch->SetAddress(&digits) ;
   branch->GetEntry(0) ; 
   MakeDigits(digits) ; 
}

//____________________________________________________________________________
void AliZDCQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

//____________________________________________________________________________ 
void AliZDCQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kZDC, task, list);  
}
