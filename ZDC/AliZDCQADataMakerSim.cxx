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

ClassImp(AliZDCQADataMakerSim)
           
//____________________________________________________________________________ 
  AliZDCQADataMakerSim::AliZDCQADataMakerSim() : 
      AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kZDC), "ZDC Quality Assurance Data Maker")
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
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1F * hHitsZNCTot = new TH1F("hHitsZNCTot", "Signal in ZNC; N_{phe}", 100, 0., 6000.);
  TH1F * hHitsZNATot = new TH1F("hHitsZNATot", "Signal in ZNA; N_{phe}", 100, 0., 6000.);
  TH1F * hHitsZPCTot = new TH1F("hHitsZPCTot", "Signal in ZPC; N_{phe}", 100, 0., 6000.);
  TH1F * hHitsZPATot = new TH1F("hHitsZPATot", "Signal in ZPA; N_{phe}", 100, 0., 6000.);
  Add2HitsList(hHitsZNCTot, 0, !expert, image);
  Add2HitsList(hHitsZNATot, 1, !expert, image);
  Add2HitsList(hHitsZPCTot, 2, !expert, image);
  Add2HitsList(hHitsZPATot, 3, !expert, image);
  //
  TH1F * hHitsSumQZNC = new TH1F("hHitsSumQZNC", "Signal in 4 ZNC PMQ; N_{phe}",100, 0., 4000.);
  TH1F * hHitsSumQZNA = new TH1F("hHitsSumQZNA", "Signal in 4 ZNA PMQ; N_{phe}",100, 0., 4000.);
  TH1F * hHitsSumQZPC = new TH1F("hHitsSumQZPC", "Signal in 4 ZPC PMQ; N_{phe}",100, 0., 4000.);
  TH1F * hHitsSumQZPA = new TH1F("hHitsSumQZPA", "Signal in 4 ZPA PMQ; N_{phe}",100, 0., 4000.);
  Add2HitsList(hHitsSumQZNC, 4, expert, !image);
  Add2HitsList(hHitsSumQZNA, 5, expert, !image);
  Add2HitsList(hHitsSumQZPC, 6, expert, !image);
  Add2HitsList(hHitsSumQZPA, 7, expert, !image);
  //
  TH1F * hHitsPMCZNC = new TH1F("hHitsPMCZNC", "Signal in ZNC PMC; N_{phe}",100, 0., 4000.);
  TH1F * hHitsPMCZNA = new TH1F("hHitsPMCZNA", "Signal in ZNA PMC; N_{phe}",100, 0., 4000.);
  TH1F * hHitsPMCZPC = new TH1F("hHitsPMCZPC", "Signal in ZPC PMC; N_{phe}",100, 0., 4000.);
  TH1F * hHitsPMCZPA = new TH1F("hHitsPMCZPA", "Signal in ZPA PMC; N_{phe}",100, 0., 4000.);
  Add2HitsList(hHitsPMCZNC, 8, expert, !image);
  Add2HitsList(hHitsPMCZNA, 9, expert, !image);
  Add2HitsList(hHitsPMCZPC, 10, expert, !image);
  Add2HitsList(hHitsPMCZPA, 11, expert, !image);
  
  TH2F * hHitsZNCh  = new TH2F("hHitsZNCh", "Hits centroid in ZNC;Centroid position [cm];Counts", 100, -5.,5.,100,-5.,5.);
  TH2F * hHitsZNAh  = new TH2F("hHitsZNAh", "Hits centroid in ZNA;Centroid position [cm];Counts", 100, -5.,5.,100,-5.,5.);
  Add2HitsList(hHitsZNCh, 12, !expert, image);
  Add2HitsList(hHitsZNAh, 13, !expert, image);
  // NB -> For the moment no check is performesd on ZP centroids
  TH2F * hHitsZPCh  = new TH2F("hHitsZPCh", "Hits centroid in ZPC;Centroid position [cm];Counts", 100,-12.,12.,100,-12.,12.); 
  TH2F * hHitsZPAh  = new TH2F("hHitsZPAh", "Hits centroid in ZPA;Centroid position [cm];Counts", 100,-12.,12.,100,-12.,12.); 
  Add2HitsList(hHitsZPCh, 14, !expert, image);
  Add2HitsList(hHitsZPAh, 15, !expert, image);
}


//____________________________________________________________________________ 
void AliZDCQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  //
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  // ------------------- HIGH GAIN CHAIN ---------------------------
  TH1F * hDigZNCTot = new TH1F("hDigZNCTot", "Signal in ZNC;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  TH1F * hDigZNATot = new TH1F("hDigZNATot", "Signal in ZNA;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  TH1F * hDigZPCTot = new TH1F("hDigZPCTot", "Signal in ZPC;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  TH1F * hDigZPATot = new TH1F("hDigZPATot", "Signal in ZPA;Amplitude [ADC counts];Counts", 100, 0., 6000.);
  Add2DigitsList(hDigZNCTot, 0, !expert, image);
  Add2DigitsList(hDigZNATot, 1, !expert, image);
  Add2DigitsList(hDigZPCTot, 2, !expert, image);
  Add2DigitsList(hDigZPATot, 3, !expert, image);
  //
  TH1F * hDigSumQZNC = new TH1F("hDigSumQZNC", "Signal in 4 ZNC PMQ;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigSumQZNA = new TH1F("hDigSumQZNA", "Signal in 4 ZNA PMQ;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigSumQZPC = new TH1F("hDigSumQZPC", "Signal in 4 ZPC PMQ;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigSumQZPA = new TH1F("hDigSumQZPA", "Signal in 4 ZPA PMQ;Amplitude [ADC counts];Counts",100, 0., 4000.);
  Add2DigitsList(hDigSumQZNC, 4, expert, !image);
  Add2DigitsList(hDigSumQZNA, 5, expert, !image);
  Add2DigitsList(hDigSumQZPC, 6, expert, !image);
  Add2DigitsList(hDigSumQZPA, 7, expert, !image);
  //
  TH1F * hDigPMCZNC = new TH1F("hDigPMCZNC", "Signal in ZNC PMC;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigPMCZNA = new TH1F("hDigPMCZNA", "Signal in ZNA PMC;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigPMCZPC = new TH1F("hDigPMCZPC", "Signal in ZPC PMC;Amplitude [ADC counts];Counts",100, 0., 4000.);
  TH1F * hDigPMCZPA = new TH1F("hDigPMCZPA", "Signal in ZPA PMC;Amplitude [ADC counts];Counts",100, 0., 4000.);
  Add2DigitsList(hDigPMCZNC, 8, expert, !image);
  Add2DigitsList(hDigPMCZNA, 9, expert, !image);
  Add2DigitsList(hDigPMCZPC, 10, expert, !image);
  Add2DigitsList(hDigPMCZPA, 11, expert, !image);
  // 
  // ------------------- LOW GAIN CHAIN ---------------------------
  TH1F * hDigZNCTotlg = new TH1F("hDigZNCTotlg", "Digit lg signal in ZNC", 100, 0., 6000.);
  TH1F * hDigZNATotlg = new TH1F("hDigZNATotlg", "Digit lg signal in ZNA", 100, 0., 6000.);
  TH1F * hDigZPCTotlg = new TH1F("hDigZPCTotlg", "Digit lg signal in ZPC", 100, 0., 6000.);
  TH1F * hDigZPATotlg = new TH1F("hDigZPATotlg", "Digit lg signal in ZPA", 100, 0., 6000.);
  Add2DigitsList(hDigZNCTotlg, 12, expert, !image);
  Add2DigitsList(hDigZNATotlg, 13, expert, !image);
  Add2DigitsList(hDigZPCTotlg, 14, expert, !image);
  Add2DigitsList(hDigZPATotlg, 15, expert, !image);
  //
  TH1F * hDigSumQZNClg = new TH1F("hDigSumQZNClg", "Signal in 4 ZNC PMQlg",100, 0., 4000.);
  TH1F * hDigSumQZNAlg = new TH1F("hDigSumQZNAlg", "Signal in 4 ZNA PMQlg",100, 0., 4000.);
  TH1F * hDigSumQZPClg = new TH1F("hDigSumQZPClg", "Signal in 4 ZPC PMQlg",100, 0., 4000.);
  TH1F * hDigSumQZPAlg = new TH1F("hDigSumQZPAlg", "Signal in 4 ZPA PMQlg",100, 0., 4000.);
  Add2DigitsList(hDigSumQZNClg, 16, expert, !image);
  Add2DigitsList(hDigSumQZNAlg, 17, expert, !image);
  Add2DigitsList(hDigSumQZPClg, 18, expert, !image);
  Add2DigitsList(hDigSumQZPAlg, 19, expert, !image);
  //
  TH1F * hDigPMCZNClg = new TH1F("hDigPMCZNClg", "Signal in ZNC PMClg",100, 0., 4000.);
  TH1F * hDigPMCZNAlg = new TH1F("hDigPMCZNAlg", "Signal in ZNA PMClg",100, 0., 4000.);
  TH1F * hDigPMCZPClg = new TH1F("hDigPMCZPClg", "Signal in ZPC PMClg",100, 0., 4000.);
  TH1F * hDigPMCZPAlg = new TH1F("hDigPMCZPAlg", "Signal in ZPA PMClg",100, 0., 4000.);
  Add2DigitsList(hDigPMCZNClg, 20, expert, !image);
  Add2DigitsList(hDigPMCZNAlg, 21, expert, !image);
  Add2DigitsList(hDigPMCZPClg, 22, expert, !image);
  Add2DigitsList(hDigPMCZPAlg, 23, expert, !image);

}

//____________________________________________________________________________
void AliZDCQADataMakerSim::MakeHits()
{
  //filling QA histos for Hits

  // Check id histograms already created for this Event Specie
  if( !GetHitsData(0) ) InitHits();
  
  TIter next(fHitsArray); 
  AliZDCHit * hit; 
  Float_t adcSum_ZNC=0., adcSum_ZNA=0., adcSum_ZPC=0., adcSum_ZPA=0.;
  Float_t adcSumQ_ZNC=0., adcSumQ_ZNA=0., adcSumQ_ZPC=0., adcSumQ_ZPA=0.;
  while((hit = dynamic_cast<AliZDCHit *>(next()))){
    if(hit->GetVolume(0)==1){
       adcSumQ_ZNC += hit->GetLightPMQ();
       adcSum_ZNC  += hit->GetLightPMC() + hit->GetLightPMQ();
       //
       GetHitsData(8)->Fill(hit->GetLightPMQ());
       //
       GetHitsData(12)->Fill(hit->GetXImpact(),hit->GetYImpact());        
    }
    else if(hit->GetVolume(0)==4){
       adcSumQ_ZNA += hit->GetLightPMQ();
       adcSum_ZNA  += hit->GetLightPMC() + hit->GetLightPMQ();
       //
       GetHitsData(9)->Fill(hit->GetLightPMQ());
       //
       GetHitsData(13)->Fill(hit->GetXImpact(), hit->GetYImpact());
    }
    else if(hit->GetVolume(0)==2){
       adcSumQ_ZNC += hit->GetLightPMQ();
       adcSum_ZNC  += hit->GetLightPMC() + hit->GetLightPMQ();
       //
       GetHitsData(10)->Fill(hit->GetLightPMQ());
       //
       GetHitsData(14)->Fill(hit->GetXImpact(), hit->GetYImpact());
    }
    else if(hit->GetVolume(0)==5){
       adcSumQ_ZNC += hit->GetLightPMQ();
       adcSum_ZNC  += hit->GetLightPMC() + hit->GetLightPMQ();
       //
       GetHitsData(11)->Fill(hit->GetLightPMQ());
       //
       GetHitsData(15)->Fill(hit->GetXImpact(), hit->GetYImpact());
    }
    //
    GetHitsData(0)->Fill(adcSum_ZNC);
    GetHitsData(1)->Fill(adcSum_ZNA);
    GetHitsData(2)->Fill(adcSum_ZPC);
    GetHitsData(3)->Fill(adcSum_ZPA);
    //
    GetHitsData(4)->Fill(adcSumQ_ZNC);
    GetHitsData(5)->Fill(adcSumQ_ZNA);
    GetHitsData(6)->Fill(adcSumQ_ZPC);
    GetHitsData(7)->Fill(adcSumQ_ZPA);
  }
}

//___________________________________________________________________________
void AliZDCQADataMakerSim::MakeHits(TTree * hitTree)
{
  // make QA data from Hit Tree
  
  if(fHitsArray) fHitsArray->Clear() ; 
  else fHitsArray = new TClonesArray("AliZDCHit", 1000);

  TBranch * branch = hitTree->GetBranch("ZDC") ;

  if(!branch){
    AliError("ZDC branch in Hit Tree not found!"); 
    return;
  } 
  else{
    Int_t nHits = 0;
    branch->SetAddress(&fHitsArray) ;
    for (Int_t ientry = 0 ; ientry < branch->GetEntries() ; ientry++) {
      branch->GetEntry(ientry) ;
      nHits += fHitsArray->GetEntriesFast();
      MakeHits() ; 
      fHitsArray->Clear();
    } 	
  }
}

//___________________________________________________________________________
void AliZDCQADataMakerSim::MakeDigits()
{
  // makes data from Digits
  if( !GetDigitsData(0) ) InitDigits();
  
  Int_t nentries = fDigitsArray->GetEntriesFast();
  if(nentries==0) printf(" AliZDCQADataMakerSim: NO entries in digit array\n\n");
  
  TIter next(fDigitsArray); 
  AliZDCDigit * digit;
     
  Float_t adcSum_ZNC=0., adcSum_ZNA=0., adcSum_ZPC=0., adcSum_ZPA=0.;
  Float_t adcSumQ_ZNC=0., adcSumQ_ZNA=0., adcSumQ_ZPC=0., adcSumQ_ZPA=0.;
  Float_t adcSum_ZNC_lg=0., adcSum_ZNA_lg=0., adcSum_ZPC_lg=0., adcSum_ZPA_lg=0.;
  Float_t adcSumQ_ZNC_lg=0., adcSumQ_ZNA_lg=0., adcSumQ_ZPC_lg=0., adcSumQ_ZPA_lg=0.;
  
  while ( (digit = dynamic_cast<AliZDCDigit *>(next())) ) {
      if(digit->GetSector(0)==1){
	  adcSum_ZNC += digit->GetADCValue(0);
	  adcSum_ZNC_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZNC += digit->GetADCValue(0);
	      adcSumQ_ZNC_lg+= digit->GetADCValue(1);
	  }
	  else{
	      GetDigitsData(8)->Fill(digit->GetADCValue(0));
	      GetDigitsData(20)->Fill(digit->GetADCValue(1));
	  }
      }
      else if(digit->GetSector(0)==2){
	  adcSum_ZPC += digit->GetADCValue(0);
	  adcSum_ZPC_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZPC += digit->GetADCValue(0);
	      adcSumQ_ZPC_lg+= digit->GetADCValue(1);
	  }
	  else{
	      GetDigitsData(10)->Fill(digit->GetADCValue(0));
	      GetDigitsData(22)->Fill(digit->GetADCValue(1));
	  }
      }
      else if(digit->GetSector(0)==4){
	  adcSum_ZNA += digit->GetADCValue(0);
	  adcSum_ZNA_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZNA += digit->GetADCValue(0);
	      adcSumQ_ZNA_lg+= digit->GetADCValue(1);
	  }
	  else{
	      GetDigitsData(9)->Fill(digit->GetADCValue(0));
	      GetDigitsData(21)->Fill(digit->GetADCValue(1));
	  }
      }
      else if(digit->GetSector(0)==5){
	  adcSum_ZPA += digit->GetADCValue(0);
	  adcSum_ZPA_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZPA += digit->GetADCValue(0);
	      adcSumQ_ZPA_lg+= digit->GetADCValue(1);
	  }
	  else{
	      GetDigitsData(11)->Fill(digit->GetADCValue(0));
	      GetDigitsData(23)->Fill(digit->GetADCValue(1));
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
  GetDigitsData(12)->Fill(adcSum_ZNC_lg);
  GetDigitsData(13)->Fill(adcSum_ZNA_lg);
  GetDigitsData(14)->Fill(adcSum_ZPC_lg);
  GetDigitsData(15)->Fill(adcSum_ZPA_lg);
  //
  GetDigitsData(16)->Fill(adcSumQ_ZNC_lg);
  GetDigitsData(17)->Fill(adcSumQ_ZNA_lg);
  GetDigitsData(18)->Fill(adcSumQ_ZPC_lg);
  GetDigitsData(19)->Fill(adcSumQ_ZPA_lg);

}

//___________________________________________________________________________
void AliZDCQADataMakerSim::MakeDigits(TTree *digitTree)
{
  // makes data from Digit Tree
  if(fDigitsArray) fDigitsArray->Clear() ; 
  else fDigitsArray = new TClonesArray("AliZDCDigit", 1000) ; 
   
  TBranch * branch = digitTree->GetBranch("ZDC");
  if(!branch){
    AliError("ZDC branch in Digit Tree not found"); 
    return;
  } 
  
  branch->SetAddress(&fDigitsArray);
  branch->GetEntry(0) ; 
  MakeDigits() ; 
  fDigitsArray->Clear();
}

//____________________________________________________________________________
void AliZDCQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

//____________________________________________________________________________ 
void AliZDCQADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray ** list)
{
  // Detector specific actions at end of cycle
  // do the QA checking
  AliQAChecker::Instance()->Run(AliQAv1::kZDC, task, list);  
}
