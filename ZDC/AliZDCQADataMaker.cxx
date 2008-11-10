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
#include "AliZDCQADataMaker.h"
#include "AliZDCHit.h"
#include "AliZDCDigit.h"
#include "AliZDCRawStream.h"
#include "AliESDZDC.h"
#include "AliESDEvent.h"

ClassImp(AliZDCQADataMaker)
           
//____________________________________________________________________________ 
  AliZDCQADataMaker::AliZDCQADataMaker() : 
      AliQADataMaker(AliQA::GetDetName(AliQA::kZDC), "ZDC Quality Assurance Data Maker"),
      fHits("AliZDCHit", 1000),
      fDigits("AliZDCDigit", 1000)

{
  // ctor
}

//____________________________________________________________________________ 
AliZDCQADataMaker::AliZDCQADataMaker(const AliZDCQADataMaker& qadm) :
  AliQADataMaker(), 
    fHits("AliZDCHit", 1000),
    fDigits("AliZDCDigit", 1000) 
 
{
  //copy ctor 
  SetName((const char*)qadm.GetName()); 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliZDCQADataMaker& AliZDCQADataMaker::operator = (const AliZDCQADataMaker& qadm )
{
  // Equal operator.
  this->~AliZDCQADataMaker();
  new(this) AliZDCQADataMaker(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliZDCQADataMaker::InitHits()
{
  // create Hits histograms in Hits subdir
  //
  TH2F * hZNCh  = new TH2F("hZNCh", "Hits centroid in ZNC", 100, -5.,5.,100,-5.,5.);
  TH2F * hZNAh  = new TH2F("hZNAh", "Hits centroid in ZNA", 100, -5.,5.,100,-5.,5.);
  TH2F * hZPCh  = new TH2F("hZPCh", "Hits centroid in ZPC", 100,-12.,12.,100,-12.,12.) 
  TH2F * hZPAh  = new TH2F("hZPAh", "Hits centroid in ZPA", 100,-12.,12.,100,-12.,12.) 
  Add2HitsList(hZNCh, 0);
  Add2HitsList(hZPCh, 1);
  Add2HitsList(hZNAh, 2);
  Add2HitsList(hZPAh, 3);
}

//____________________________________________________________________________ 
void AliZDCQADataMaker::InitDigits()
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

void AliZDCQADataMakerRec::InitRaws()
{
  // create Digits histograms in Digits subdir
  //
  TH1F * hRawZNCTot = new TH1F("hRawZNCTot", "Raw signal in ZNC", 100, 0., 6000.);
  TH1F * hRawZNATot = new TH1F("hRawZNATot", "Raw signal in ZNA", 100, 0., 6000.);
  TH1F * hRawZPCTot = new TH1F("hRawZPCTot", "Raw signal in ZPC", 100, 0., 10000.);
  TH1F * hRawZPATot = new TH1F("hRawZPATot", "Raw signal in ZPA", 100, 0., 10000.);
  Add2RawsList(hRawZNCTot, 0);
  Add2RawsList(hRawZPCTot, 1);
  Add2RawsList(hRawZNATot, 2);
  Add2RawsList(hRawZPATot, 3);
  //
  TH1F * hRawSumQZNC = new TH1F("hRawSumQZNC", "Raw summed 4 ZNC quadrants",100, 0., 4000.);
  TH1F * hRawSumQZPC = new TH1F("hRawSumQZPC", "Raw summed 4 ZPC quadrants",100, 0., 4000.);
  TH1F * hRawSumQZNA = new TH1F("hRawSumQZNA", "Raw summed 4 ZNA quadrants",100, 0., 4000.);
  TH1F * hRawSumQZPA = new TH1F("hRawSumQZPA", "Raw summed 4 ZPA quadrants",100, 0., 4000.);
  Add2RawsList(hRawSumQZNC, 4);
  Add2RawsList(hRawSumQZPC, 5);
  Add2RawsList(hRawSumQZNA, 6);
  Add2RawsList(hRawSumQZPA, 7);
  //
  TH1F * hRawPMCZNC = new TH1F("hRawPMCZNC", "Raw common ZNC PMT",100, 0., 4000.);
  TH1F * hRawPMCZPC = new TH1F("hRawPMCZPC", "Raw common ZPC PMT",100, 0., 4000.);
  TH1F * hRawPMCZNA = new TH1F("hRawPMCZNA", "Raw common ZNA PMT",100, 0., 4000.);
  TH1F * hRawPMCZPA = new TH1F("hRawPMCZPA", "Raw common ZPA PMT",100, 0., 4000.);
  Add2RawsList(hRawPMCZNC, 8);
  Add2RawsList(hRawPMCZPC, 9);
  Add2RawsList(hRawPMCZNA, 10);
  Add2RawsList(hRawPMCZPA, 11);
  // 
  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hRawZNCTotlg = new TH1F("hRawZNCTotlg", "Rawit lg signal in ZNC", 100, 0., 6000.);
  TH1F * hRawZNATotlg = new TH1F("hRawZNATotlg", "Rawit lg signal in ZNA", 100, 0., 6000.);
  TH1F * hRawZPCTotlg = new TH1F("hRawZPCTotlg", "Rawit lg signal in ZPC", 100, 0., 10000.);
  TH1F * hRawZPATotlg = new TH1F("hRawZPATotlg", "Rawit lg signal in ZPA", 100, 0., 10000.);
  Add2RawsList(hRawZNCTotlg, 12);
  Add2RawsList(hRawZPCTotlg, 13);
  Add2RawsList(hRawZNATotlg, 14);
  Add2RawsList(hRawZPATotlg, 15);
  //
  TH1F * hRawSumQZNClg = new TH1F("hRawSumQZNClg", "Raw summed 4 lg ZNC quadrants",100, 0., 4000.);
  TH1F * hRawSumQZPClg = new TH1F("hRawSumQZPClg", "Raw summed 4 lg ZPC quadrants",100, 0., 4000.);
  TH1F * hRawSumQZNAlg = new TH1F("hRawSumQZNAlg", "Raw summed 4 lg ZNA quadrants",100, 0., 4000.);
  TH1F * hRawSumQZPAlg = new TH1F("hRawSumQZPAlg", "Raw summed 4 lg ZPA quadrants",100, 0., 4000.);
  Add2RawsList(hRawSumQZNClg, 16);
  Add2RawsList(hRawSumQZPClg, 17);
  Add2RawsList(hRawSumQZNAlg, 18);
  Add2RawsList(hRawSumQZPAlg, 19);
  //
  TH1F * hRawPMCZNClg = new TH1F("hRawPMCZNClg", "Raw common lg ZNC PMT",100, 0., 4000.);
  TH1F * hRawPMCZPClg = new TH1F("hRawPMCZPClg", "Raw common lg ZPC PMT",100, 0., 4000.);
  TH1F * hRawPMCZNAlg = new TH1F("hRawPMCZNAlg", "Raw common lg ZNA PMT",100, 0., 4000.);
  TH1F * hRawPMCZPAlg = new TH1F("hRawPMCZPAlg", "Raw common lg ZPA PMT",100, 0., 4000.);
  Add2RawsList(hRawPMCZNClg, 20);
  Add2RawsList(hRawPMCZPClg, 21);
  Add2RawsList(hRawPMCZNAlg, 22);
  Add2RawsList(hRawPMCZPAlg, 23);
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::InitESDs()
{
  //Booking ESDs histograms
  //
  TH2F * hZNC  = new TH2F("hZNC", "Centroid in ZNC", 100, -5.,5.,100,-5.,5.);
  TH2F * hZNA  = new TH2F("hZNA", "Centroid in ZNA", 100, -5.,5.,100,-5.,5.);
  Add2DigitsList(hZNC, 0);
  Add2DigitsList(hZNA, 1);
  //
  TH1F * hESDZNCTot = new TH1F("hESDZNCTot", "ESD signal in ZNC", 100, 0., 6000.);
  TH1F * hESDZPCTot = new TH1F("hESDZPCTot", "ESD signal in ZPC", 100, 0., 10000.);
  TH1F * hESDZNATot = new TH1F("hESDZNATot", "ESD signal in ZNA", 100, 0., 6000.);
  TH1F * hESDZPATot = new TH1F("hESDZPATot", "ESD signal in ZPA", 100, 0., 10000.);
  Add2ESDsList(hESDZNCTot, 2);
  Add2ESDsList(hESDZPCTot, 3);
  Add2ESDsList(hESDZNATot, 4);
  Add2ESDsList(hESDZPATot, 5);
  //
  TH1F * hESDSumQZNC = new TH1F("hESDSumQZNC", "Sum of 4 ZNC sectors",100, 0., 4000.);
  TH1F * hESDSumQZPC = new TH1F("hESDSumQZPC", "Sum of 4 ZPC sectors",100, 0., 4000.);
  TH1F * hESDSumQZNA = new TH1F("hESDSumQZNA", "Sum of 4 ZNA sectors",100, 0., 4000.);
  TH1F * hESDSumQZPA = new TH1F("hESDSumQZPA", "Sum of 4 ZPA sectors",100, 0., 4000.);
  Add2ESDsList(hESDSumQZNC, 6);
  Add2ESDsList(hESDSumQZPC, 7);
  Add2ESDsList(hESDSumQZNA, 8);
  Add2ESDsList(hESDSumQZPA, 9);
  //
  TH1F * hESDPMCZNC = new TH1F("hESDPMCZNC", "Signal in common ZNC PMT",100, 0., 4000.);
  TH1F * hESDPMCZPC = new TH1F("hESDPMCZPC", "Signal in common ZPC PMT",100, 0., 4000.);
  TH1F * hESDPMCZNA = new TH1F("hESDPMCZNA", "Signal in common ZNA PMT",100, 0., 4000.);
  TH1F * hESDPMCZPA = new TH1F("hESDPMCZPA", "Signal in common ZPA PMT",100, 0., 4000.);
  Add2ESDsList(hESDPMCZNC, 10);
  Add2ESDsList(hESDPMCZPC, 11);
  Add2ESDsList(hESDPMCZNA, 12);
  Add2ESDsList(hESDPMCZPA, 13);
  // 
  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hESDZNCTotlg = new TH1F("hESDZNCTotlg", "ESD lg signal in ZNC", 100, 0., 6000.);
  TH1F * hESDZNATotlg = new TH1F("hESDZNATotlg", "ESD lg signal in ZNA", 100, 0., 6000.);
  TH1F * hESDZPCTotlg = new TH1F("hESDZPCTotlg", "ESD lg signal in ZPC", 100, 0., 10000.);
  TH1F * hESDZPATotlg = new TH1F("hESDZPATotlg", "ESD lg signal in ZPA", 100, 0., 10000.);
  Add2ESDsList(hESDZNCTotlg, 14);
  Add2ESDsList(hESDZPCTotlg, 15);
  Add2ESDsList(hESDZNATotlg, 16);
  Add2ESDsList(hESDZPATotlg, 17);
  //
  TH1F * hESDSumQZNClg = new TH1F("hESDSumQZNClg", "Sum of 4 lg ZNC sectors",100, 0., 4000.);
  TH1F * hESDSumQZPClg = new TH1F("hESDSumQZPClg", "Sum of 4 lg ZPC sectors",100, 0., 4000.);
  TH1F * hESDSumQZNAlg = new TH1F("hESDSumQZNAlg", "Sum of 4 lg ZNA sectors",100, 0., 4000.);
  TH1F * hESDSumQZPAlg = new TH1F("hESDSumQZPAlg", "Sum of 4 lg ZPA sectors",100, 0., 4000.);
  Add2ESDsList(hESDSumQZNClg, 18);
  Add2ESDsList(hESDSumQZPClg, 19);
  Add2ESDsList(hESDSumQZNAlg, 20);
  Add2ESDsList(hESDSumQZPAlg, 21);
  //
  TH1F * hESDPMCZNClg = new TH1F("hESDPMCZNClg", "Signal in common ZNC lg PMT",100, 0., 4000.);
  TH1F * hESDPMCZPClg = new TH1F("hESDPMCZPClg", "Signal in common ZPC lg PMT",100, 0., 4000.);
  TH1F * hESDPMCZNAlg = new TH1F("hESDPMCZNAlg", "Signal in common ZNA lg PMT",100, 0., 4000.);
  TH1F * hESDPMCZPAlg = new TH1F("hESDPMCZPAlg", "Signal in common ZPA lg PMT",100, 0., 4000.);
  Add2ESDsList(hESDPMCZNClg, 22);
  Add2ESDsList(hESDPMCZPClg, 23);
  Add2ESDsList(hESDPMCZNAlg, 24);
  Add2ESDsList(hESDPMCZPAlg, 25);
}
  

//____________________________________________________________________________
void AliZDCQADataMakerSim::MakeHits(TClonesArray */*data*/)
{
  //filling QA histos for Hits
  //
  TIter next(&fHits); 
  AliZDCHit * hit; 
  while((hit = dynamic_cast<AliZDCHit *>(next()))){
    if(hit->GetVolume(0)==1) GetHitsData(0)->Fill(hit->GetXImpact(),hit->GetYImpact());
    else if(hit->GetVolume(0)==2) GetHitsData(1)->Fill(hit->GetXImpact(), hit->GetYImpact());
    else if(hit->GetVolume(0)==4) GetHitsData(2)->Fill(hit->GetXImpact(), hit->GetYImpact());
    else if(hit->GetVolume(0)==5) GetHitsData(3)->Fill(hit->GetXImpact(), hit->GetYImpact());
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
    Int_t ntracks = (Int_t) hitTree->GetEntries();
    //printf("\n\t *** no.track %d\n",ntracks);
    if (ntracks<=0) return;
    //
    for(Int_t itrack=0; itrack<ntracks; itrack++){
      fHits = new TClonesArray("AliZDCHit", 1000);
      //
      branch->SetAddress(&fHits) ;
      branch->GetEntry(itrack);
      //
      //printf("\t *** track %d",itrack);
      //hits->Print("");
      //printf("\n");
      //
      MakeHits(fHits); 
      fHits->Clear();
    }	  
  }
}

//____________________________________________________________________________
void AliZDCQADataMakerSim::MakeDigits(TClonesArray * /*digits*/)
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
   TBranch * branch = digitTree->GetBranch("ZDC");
   if(!branch){
      AliError("ZDC branch in Digit Tree not found"); 
      return;
   } 
   branch->SetAddress(&fDigits);
   branch->GetEntry(0); 
   MakeDigits(fDigits); 
   fDigits->Clear();
}

//____________________________________________________________________________

void AliZDCQADataMakerRec::MakeRaws(AliRawReader *rawReader)
{
  // Filling Raws QA histos
  //
  Float_t Sum_ZNC=0., Sum_ZNA=0., Sum_ZPC=0., Sum_ZPA=0.;
  Float_t SumQ_ZNC=0., SumQ_ZNA=0., SumQ_ZPC=0., SumQ_ZPA=0.;
  Float_t Sum_ZNC_lg=0., Sum_ZNA_lg=0., Sum_ZPC_lg=0., Sum_ZPA_lg=0.;
  Float_t SumQ_ZNC_lg=0., SumQ_ZNA_lg=0., SumQ_ZPC_lg=0., SumQ_ZPA_lg=0.;
  //
  AliZDCRawStream stream(rawReader);
  while(stream.Next()){
    if(stream.IsADCDataWord() && 
     (stream.GetADCModule()==0 || stream.GetADCModule()==1)){
       if(stream.GetSector(0)==1){
         if(stream.GetADCGain()==0){
	   Sum_ZNC += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) SumQ_ZNC += stream.GetADCValue();
	   else GetRawsData(8)->Fill(stream.GetADCValue());
	 }
	 else{
	   Sum_ZNC_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) SumQ_ZNC_lg += stream.GetADCValue();
	   else GetRawsData(20)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==2){
         if(stream.GetADCGain()==0){
	   Sum_ZPC += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) SumQ_ZPC += stream.GetADCValue();
	   else GetRawsData(9)->Fill(stream.GetADCValue());
	 }
	 else{
	   Sum_ZPC_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) SumQ_ZPC_lg += stream.GetADCValue();
	   else GetRawsData(21)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==4){
         if(stream.GetADCGain()==0){
	   Sum_ZNA += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) SumQ_ZNA += stream.GetADCValue();
	   else GetRawsData(10)->Fill(stream.GetADCValue());
	 }
	 else{
	   Sum_ZNA_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) SumQ_ZNA_lg += stream.GetADCValue();
	   else GetRawsData(22)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==5){
         if(stream.GetADCGain()==0){
	   Sum_ZPA += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) SumQ_ZPA += stream.GetADCValue();
	   else GetRawsData(11)->Fill(stream.GetADCValue());
	 }
	 else{
	   Sum_ZPA_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) SumQ_ZPA_lg += stream.GetADCValue();
	   else GetRawsData(23)->Fill(stream.GetADCValue());
	 }
       }
    }
  }
  //
  GetRawsData(0)->Fill(Sum_ZNC);
  GetRawsData(1)->Fill(Sum_ZPC);
  GetRawsData(2)->Fill(Sum_ZNA);
  GetRawsData(3)->Fill(Sum_ZPA);
  //
  GetRawsData(4)->Fill(SumQ_ZNC);
  GetRawsData(5)->Fill(SumQ_ZPC);
  GetRawsData(6)->Fill(SumQ_ZNA);
  GetRawsData(7)->Fill(SumQ_ZPA);
  //
  GetRawsData(12)->Fill(Sum_ZNC_lg);
  GetRawsData(13)->Fill(Sum_ZPC_lg);
  GetRawsData(14)->Fill(Sum_ZNA_lg);
  GetRawsData(15)->Fill(Sum_ZPA_lg);
  //
  GetRawsData(16)->Fill(SumQ_ZNC_lg);
  GetRawsData(17)->Fill(SumQ_ZPC_lg);
  GetRawsData(18)->Fill(SumQ_ZNA_lg);
  GetRawsData(19)->Fill(SumQ_ZPA_lg);
  //
  stream.Delete();
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
  //
  AliESDZDC * zdcESD =  esd->GetESDZDC();
  //
  const Float_t * Centr_ZNC, * Centr_ZNA;
  Int_t NSpecnC = (Int_t) (esd->GetZDCN1Energy()/2.7);
  Int_t NSpecnA = (Int_t) (esd->GetZDCN2Energy()/2.7);
  Centr_ZNC = zdcESD->GetZNCCentroid();
  Centr_ZNA = zdcESD->GetZNACentroid();
  GetESDsData(0)->Fill(Centr_ZNC[0], Centr_ZNC[1]);
  GetESDsData(1)->Fill(Centr_ZNA[0], Centr_ZNA[1]);
  //
  GetESDsData(2)->Fill(esd->GetZDCN1Energy());
  GetESDsData(3)->Fill(esd->GetZDCP1Energy());
  GetESDsData(4)->Fill(esd->GetZDCN2Energy());
  GetESDsData(5)->Fill(esd->GetZDCP2Energy());
  //
  Double_t ZNCSumQ=0., ZPCSumQ=0., ZNASumQ=0., ZPASumQ=0.;
  Double_t ZNCSumQ_lg=0., ZPCSumQ_lg=0., ZNASumQ_lg=0., ZPASumQ_lg=0.;
  //
  const Double_t *ZNCTow, *ZPCTow, *ZNATow, *ZPATow;
  const Double_t *ZNCTow_lg, *ZPCTow_lg, *ZNATow_lg, *ZPATow_lg;
  //
  ZNCTow = zdcESD->GetZN1TowerEnergy();
  ZPCTow = zdcESD->GetZP1TowerEnergy();
  ZNATow = zdcESD->GetZN2TowerEnergy();
  ZPATow = zdcESD->GetZP2TowerEnergy();
  //
  ZNCTow_lg = zdcESD->GetZN1TowerEnergyLR();
  ZPCTow_lg = zdcESD->GetZP1TowerEnergyLR();
  ZNATow_lg = zdcESD->GetZN2TowerEnergyLR();
  ZPATow_lg = zdcESD->GetZP2TowerEnergyLR();
  //
  for(Int_t i=0; i<5; i++){
     if(i==0){
       GetESDsData(10)->Fill(ZNCTow[i]);
       GetESDsData(11)->Fill(ZPCTow[i]);
       GetESDsData(12)->Fill(ZNATow[i]);
       GetESDsData(13)->Fill(ZPATow[i]);
       //
       GetESDsData(22)->Fill(ZNCTow_lg[i]);
       GetESDsData(23)->Fill(ZPCTow_lg[i]);
       GetESDsData(24)->Fill(ZNATow_lg[i]);
       GetESDsData(25)->Fill(ZPATow_lg[i]);
     }
     else{
       ZNCSumQ += ZNCTow[i];
       ZPCSumQ += ZPCTow[i];
       ZNASumQ += ZNATow[i];
       ZPASumQ += ZPATow[i];
       //
       ZNCSumQ_lg += ZNCTow_lg[i];
       ZPCSumQ_lg += ZPCTow_lg[i];
       ZNASumQ_lg += ZNATow_lg[i];
       ZPASumQ_lg += ZPATow_lg[i];
     }
  }
  GetESDsData(6)->Fill(ZNCSumQ);
  GetESDsData(7)->Fill(ZPCSumQ);
  GetESDsData(8)->Fill(ZNASumQ);
  GetESDsData(9)->Fill(ZPASumQ);
  //
  GetESDsData(18)->Fill(ZNCSumQ_lg);
  GetESDsData(19)->Fill(ZPCSumQ_lg);
  GetESDsData(20)->Fill(ZNASumQ_lg);
  GetESDsData(21)->Fill(ZPASumQ_lg);
}

//____________________________________________________________________________
void AliZDCQADataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

//____________________________________________________________________________ 
void AliZDCQADataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQA::kZDC, task, list) ;  
