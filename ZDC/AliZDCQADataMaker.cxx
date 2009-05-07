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
      AliQADataMaker(AliQAv1::GetDetName(AliQAv1::kZDC), "ZDC Quality Assurance Data Maker"),
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
  TH1F * hDigsumQZNC = new TH1F("hDigsumQZNC", "Signal in 4 ZNC PMQ[i]",100, 0., 4000.);
  TH1F * hDigsumQZPC = new TH1F("hDigsumQZPC", "Signal in 4 ZPC PMQ[i]",100, 0., 4000.);
  TH1F * hDigsumQZNA = new TH1F("hDigsumQZNA", "Signal in 4 ZNA PMQ[i]",100, 0., 4000.);
  TH1F * hDigsumQZPA = new TH1F("hDigsumQZPA", "Signal in 4 ZPA PMQ[i]",100, 0., 4000.);
  Add2DigitsList(hDigsumQZNC, 4);
  Add2DigitsList(hDigsumQZPC, 5);
  Add2DigitsList(hDigsumQZNA, 6);
  Add2DigitsList(hDigsumQZPA, 7);
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
  TH1F * hDigsumQZNClg = new TH1F("hDigsumQZNClg", "Signal in 4 ZNC PMQlg[i]",100, 0., 4000.);
  TH1F * hDigsumQZPClg = new TH1F("hDigsumQZPClg", "Signal in 4 ZPC PMQlg[i]",100, 0., 4000.);
  TH1F * hDigsumQZNAlg = new TH1F("hDigsumQZNAlg", "Signal in 4 ZNA PMQlg[i]",100, 0., 4000.);
  TH1F * hDigsumQZPAlg = new TH1F("hDigsumQZPAlg", "Signal in 4 ZPA PMQlg[i]",100, 0., 4000.);
  Add2DigitsList(hDigsumQZNClg, 16);
  Add2DigitsList(hDigsumQZPClg, 17);
  Add2DigitsList(hDigsumQZNAlg, 18);
  Add2DigitsList(hDigsumQZPAlg, 19);
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

void AliZDCQADataMaker::InitRaws()
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
  TH1F * hRawsumQZNC = new TH1F("hRawsumQZNC", "Raw summed 4 ZNC quadrants",100, 0., 4000.);
  TH1F * hRawsumQZPC = new TH1F("hRawsumQZPC", "Raw summed 4 ZPC quadrants",100, 0., 4000.);
  TH1F * hRawsumQZNA = new TH1F("hRawsumQZNA", "Raw summed 4 ZNA quadrants",100, 0., 4000.);
  TH1F * hRawsumQZPA = new TH1F("hRawsumQZPA", "Raw summed 4 ZPA quadrants",100, 0., 4000.);
  Add2RawsList(hRawsumQZNC, 4);
  Add2RawsList(hRawsumQZPC, 5);
  Add2RawsList(hRawsumQZNA, 6);
  Add2RawsList(hRawsumQZPA, 7);
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
  TH1F * hRawsumQZNClg = new TH1F("hRawsumQZNClg", "Raw summed 4 lg ZNC quadrants",100, 0., 4000.);
  TH1F * hRawsumQZPClg = new TH1F("hRawsumQZPClg", "Raw summed 4 lg ZPC quadrants",100, 0., 4000.);
  TH1F * hRawsumQZNAlg = new TH1F("hRawsumQZNAlg", "Raw summed 4 lg ZNA quadrants",100, 0., 4000.);
  TH1F * hRawsumQZPAlg = new TH1F("hRawsumQZPAlg", "Raw summed 4 lg ZPA quadrants",100, 0., 4000.);
  Add2RawsList(hRawsumQZNClg, 16);
  Add2RawsList(hRawsumQZPClg, 17);
  Add2RawsList(hRawsumQZNAlg, 18);
  Add2RawsList(hRawsumQZPAlg, 19);
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
void AliZDCQADataMaker::InitESDs()
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
  TH1F * hESDsumQZNC = new TH1F("hESDsumQZNC", "sum of 4 ZNC sectors",100, 0., 4000.);
  TH1F * hESDsumQZPC = new TH1F("hESDsumQZPC", "sum of 4 ZPC sectors",100, 0., 4000.);
  TH1F * hESDsumQZNA = new TH1F("hESDsumQZNA", "sum of 4 ZNA sectors",100, 0., 4000.);
  TH1F * hESDsumQZPA = new TH1F("hESDsumQZPA", "sum of 4 ZPA sectors",100, 0., 4000.);
  Add2ESDsList(hESDsumQZNC, 6);
  Add2ESDsList(hESDsumQZPC, 7);
  Add2ESDsList(hESDsumQZNA, 8);
  Add2ESDsList(hESDsumQZPA, 9);
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
  TH1F * hESDsumQZNClg = new TH1F("hESDsumQZNClg", "sum of 4 lg ZNC sectors",100, 0., 4000.);
  TH1F * hESDsumQZPClg = new TH1F("hESDsumQZPClg", "sum of 4 lg ZPC sectors",100, 0., 4000.);
  TH1F * hESDsumQZNAlg = new TH1F("hESDsumQZNAlg", "sum of 4 lg ZNA sectors",100, 0., 4000.);
  TH1F * hESDsumQZPAlg = new TH1F("hESDsumQZPAlg", "sum of 4 lg ZPA sectors",100, 0., 4000.);
  Add2ESDsList(hESDsumQZNClg, 18);
  Add2ESDsList(hESDsumQZPClg, 19);
  Add2ESDsList(hESDsumQZNAlg, 20);
  Add2ESDsList(hESDsumQZPAlg, 21);
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
void AliZDCQADataMaker::MakeHits(TClonesArray */*data*/)
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
void AliZDCQADataMaker::MakeHits(TTree * const hitTree)
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
void AliZDCQADataMaker::MakeDigits(TClonesArray * /*digits*/)
{
  // makes data from Digits
  //
  TIter next(digits) ; 
  AliZDCDigit * digit ; 
  //
  Float_t sumADC_ZNC=0., sumADC_ZNA=0., sumADC_ZPC=0., sumADC_ZPA=0.;
  Float_t sumADCQ_ZNC=0., sumADCQ_ZNA=0., sumADCQ_ZPC=0., sumADCQ_ZPA=0.;
  Float_t sumADC_ZNC_lg=0., sumADC_ZNA_lg=0., sumADC_ZPC_lg=0., sumADC_ZPA_lg=0.;
  Float_t sumADCQ_ZNC_lg=0., sumADCQ_ZNA_lg=0., sumADCQ_ZPC_lg=0., sumADCQ_ZPA_lg=0.;
  //
  while((digit = dynamic_cast<AliZDCDigit *>(next()))){
    if(digit->GetSector(0)==1){
      sumADC_ZNC += digit->GetADCValue(0);
      sumADC_ZNC_lg += digit->GetADCValue(1);
      //
      if(digit->GetSector(1)!=0){
        sumADCQ_ZNC += digit->GetADCValue(0);
        sumADCQ_ZNC_lg+= digit->GetADCValue(1);
      }
      else{
        GetDigitsData(8)->Fill(digit->GetADCValue(0));
        GetDigitsData(20)->Fill(digit->GetADCValue(1));
      }
    }
    else if(digit->GetSector(0)==2){
      sumADC_ZPC += digit->GetADCValue(0);
      sumADC_ZPC_lg += digit->GetADCValue(1);
      //
      if(digit->GetSector(1)!=0){
        sumADCQ_ZPC += digit->GetADCValue(0);
        sumADCQ_ZPC_lg+= digit->GetADCValue(1);
      }
      else{
        GetDigitsData(9)->Fill(digit->GetADCValue(0));
        GetDigitsData(21)->Fill(digit->GetADCValue(1));
      }
    }
    else if(digit->GetSector(0)==4){
      sumADC_ZNA += digit->GetADCValue(0);
      sumADC_ZNA_lg += digit->GetADCValue(1);
      //
      if(digit->GetSector(1)!=0){
        sumADCQ_ZNA += digit->GetADCValue(0);
        sumADCQ_ZNA_lg+= digit->GetADCValue(1);
      }
      else{
        GetDigitsData(10)->Fill(digit->GetADCValue(0));
        GetDigitsData(22)->Fill(digit->GetADCValue(1));
      }
    }
    else if(digit->GetSector(0)==5){
      sumADC_ZPA += digit->GetADCValue(0);
      sumADC_ZPA_lg += digit->GetADCValue(1);
      //
      if(digit->GetSector(1)!=0){
        sumADCQ_ZPA += digit->GetADCValue(0);
        sumADCQ_ZPA_lg+= digit->GetADCValue(1);
      }
      else{
        GetDigitsData(11)->Fill(digit->GetADCValue(0));
        GetDigitsData(23)->Fill(digit->GetADCValue(1));
      }
    }
  }
  //
  GetDigitsData(0)->Fill(sumADC_ZNC);
  GetDigitsData(1)->Fill(sumADC_ZPC);
  GetDigitsData(2)->Fill(sumADC_ZNA);
  GetDigitsData(3)->Fill(sumADC_ZPA);
  //
  GetDigitsData(4)->Fill(sumADCQ_ZNC);
  GetDigitsData(5)->Fill(sumADCQ_ZPC);
  GetDigitsData(6)->Fill(sumADCQ_ZNA);
  GetDigitsData(7)->Fill(sumADCQ_ZPA);
  //
  GetDigitsData(12)->Fill(sumADC_ZNC_lg);
  GetDigitsData(13)->Fill(sumADC_ZPC_lg);
  GetDigitsData(14)->Fill(sumADC_ZNA_lg);
  GetDigitsData(15)->Fill(sumADC_ZPA_lg);
  //
  GetDigitsData(16)->Fill(sumADCQ_ZNC_lg);
  GetDigitsData(17)->Fill(sumADCQ_ZPC_lg);
  GetDigitsData(18)->Fill(sumADCQ_ZNA_lg);
  GetDigitsData(19)->Fill(sumADCQ_ZPA_lg);
   
}
//___________________________________________________________________________
void AliZDCQADataMaker::MakeDigits(TTree * const digitTree )
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

void AliZDCQADataMaker::MakeRaws(AliRawReader * const rawReader)
{
  // Filling Raws QA histos
  //
  Float_t sum_ZNC=0., sum_ZNA=0., sum_ZPC=0., sum_ZPA=0.;
  Float_t sumQ_ZNC=0., sumQ_ZNA=0., sumQ_ZPC=0., sumQ_ZPA=0.;
  Float_t sum_ZNC_lg=0., sum_ZNA_lg=0., sum_ZPC_lg=0., sum_ZPA_lg=0.;
  Float_t sumQ_ZNC_lg=0., sumQ_ZNA_lg=0., sumQ_ZPC_lg=0., sumQ_ZPA_lg=0.;
  //
  AliZDCRawStream stream(rawReader);
  while(stream.Next()){
    if(stream.IsADCDataWord() && 
     (stream.GetADCModule()==0 || stream.GetADCModule()==1)){
       if(stream.GetSector(0)==1){
         if(stream.GetADCGain()==0){
	   sum_ZNC += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZNC += stream.GetADCValue();
	   else GetRawsData(8)->Fill(stream.GetADCValue());
	 }
	 else{
	   sum_ZNC_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZNC_lg += stream.GetADCValue();
	   else GetRawsData(20)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==2){
         if(stream.GetADCGain()==0){
	   sum_ZPC += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZPC += stream.GetADCValue();
	   else GetRawsData(9)->Fill(stream.GetADCValue());
	 }
	 else{
	   sum_ZPC_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZPC_lg += stream.GetADCValue();
	   else GetRawsData(21)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==4){
         if(stream.GetADCGain()==0){
	   sum_ZNA += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZNA += stream.GetADCValue();
	   else GetRawsData(10)->Fill(stream.GetADCValue());
	 }
	 else{
	   sum_ZNA_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZNA_lg += stream.GetADCValue();
	   else GetRawsData(22)->Fill(stream.GetADCValue());
	 }
       }
       else if(stream.GetSector(0)==5){
         if(stream.GetADCGain()==0){
	   sum_ZPA += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZPA += stream.GetADCValue();
	   else GetRawsData(11)->Fill(stream.GetADCValue());
	 }
	 else{
	   sum_ZPA_lg += stream.GetADCValue();
	   if(stream.GetSector(1)!=0) sumQ_ZPA_lg += stream.GetADCValue();
	   else GetRawsData(23)->Fill(stream.GetADCValue());
	 }
       }
    }
  }
  //
  GetRawsData(0)->Fill(sum_ZNC);
  GetRawsData(1)->Fill(sum_ZPC);
  GetRawsData(2)->Fill(sum_ZNA);
  GetRawsData(3)->Fill(sum_ZPA);
  //
  GetRawsData(4)->Fill(sumQ_ZNC);
  GetRawsData(5)->Fill(sumQ_ZPC);
  GetRawsData(6)->Fill(sumQ_ZNA);
  GetRawsData(7)->Fill(sumQ_ZPA);
  //
  GetRawsData(12)->Fill(sum_ZNC_lg);
  GetRawsData(13)->Fill(sum_ZPC_lg);
  GetRawsData(14)->Fill(sum_ZNA_lg);
  GetRawsData(15)->Fill(sum_ZPA_lg);
  //
  GetRawsData(16)->Fill(sumQ_ZNC_lg);
  GetRawsData(17)->Fill(sumQ_ZPC_lg);
  GetRawsData(18)->Fill(sumQ_ZNA_lg);
  GetRawsData(19)->Fill(sumQ_ZPA_lg);
  //
  stream.Delete();
}

//____________________________________________________________________________
void AliZDCQADataMaker::MakeESDs(AliESDEvent * esd)
{
  // make QA data from ESDs
  //
  AliESDZDC * zdcESD =  esd->GetESDZDC();
  //
  const Float_t * centrZNC, * centrZNA;
  Int_t nSpecnC = (Int_t) (esd->GetZDCN1Energy()/2.7);
  Int_t nSpecnA = (Int_t) (esd->GetZDCN2Energy()/2.7);
  centrZNC = zdcESD->GetZNCCentroid();
  centrZNA = zdcESD->GetZNACentroid();
  GetESDsData(0)->Fill(centrZNC[0], centrZNC[1]);
  GetESDsData(1)->Fill(centrZNA[0], centrZNA[1]);
  //
  GetESDsData(2)->Fill(esd->GetZDCN1Energy());
  GetESDsData(3)->Fill(esd->GetZDCP1Energy());
  GetESDsData(4)->Fill(esd->GetZDCN2Energy());
  GetESDsData(5)->Fill(esd->GetZDCP2Energy());
  //
  Double_t sumZNCQ=0., sumZPCQ=0., sumZNAQ=0., sumZPAQ=0.;
  Double_t sumZNCQlg=0., sumZPCQlg=0., sumZNAQlg=0., sumZPAQlg=0.;
  //
  const Double_t *towZNC, *towZPC, *towZNA, *towZPA;
  const Double_t *towZNClg, *towZPClg, *towZNAlg, *towZPAlg;
  //
  towZNC = zdcESD->GetZN1TowerEnergy();
  towZPC = zdcESD->GetZP1TowerEnergy();
  towZNA = zdcESD->GetZN2TowerEnergy();
  towZPA = zdcESD->GetZP2TowerEnergy();
  //
  towZNClg = zdcESD->GetZN1TowerEnergyLR();
  towZPClg = zdcESD->GetZP1TowerEnergyLR();
  towZNAlg = zdcESD->GetZN2TowerEnergyLR();
  towZPAlg = zdcESD->GetZP2TowerEnergyLR();
  //
  for(Int_t i=0; i<5; i++){
     if(i==0){
       GetESDsData(10)->Fill(towZNC[i]);
       GetESDsData(11)->Fill(towZPC[i]);
       GetESDsData(12)->Fill(towZNA[i]);
       GetESDsData(13)->Fill(towZPA[i]);
       //
       GetESDsData(22)->Fill(towZNClg[i]);
       GetESDsData(23)->Fill(towZPClg[i]);
       GetESDsData(24)->Fill(towZNAlg[i]);
       GetESDsData(25)->Fill(towZPAlg[i]);
     }
     else{
       sumZNCQ += towZNC[i];
       sumZPCQ += towZPC[i];
       sumZNAQ += towZNA[i];
       sumZPAQ += towZPA[i];
       //
       sumZNCQlg += towZNClg[i];
       sumZPCQlg += towZPClg[i];
       sumZNAQlg += towZNAlg[i];
       sumZPAQlg += towZPAlg[i];
     }
  }
  GetESDsData(6)->Fill(sumZNCQ);
  GetESDsData(7)->Fill(sumZPCQ);
  GetESDsData(8)->Fill(sumZNAQ);
  GetESDsData(9)->Fill(sumZPAQ);
  //
  GetESDsData(18)->Fill(sumZNCQlg);
  GetESDsData(19)->Fill(sumZPCQlg);
  GetESDsData(20)->Fill(sumZNAQlg);
  GetESDsData(21)->Fill(sumZPAQlg);
}

//____________________________________________________________________________
void AliZDCQADataMaker::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

//____________________________________________________________________________ 
void AliZDCQADataMaker::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray * list)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kZDC, task, list) ;  
