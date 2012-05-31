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
  
  ClonePerTrigClass(AliQAv1::kHITS); // this should be the last line
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
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
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
       FillHitsData(8,hit->GetLightPMC());
    }
    else if(hit->GetVolume(0)==4){
       adcSumQ_ZNA += hit->GetLightPMQ();
       adcSum_ZNA  += hit->GetLightPMC() + hit->GetLightPMQ();
       //
       FillHitsData(9,hit->GetLightPMC());
    }
    else if(hit->GetVolume(0)==2){
       adcSumQ_ZNC += hit->GetLightPMQ();
       adcSum_ZNC  += hit->GetLightPMC() + hit->GetLightPMQ();
       //
       FillHitsData(10,hit->GetLightPMC());
    }
    else if(hit->GetVolume(0)==5){
       adcSumQ_ZNC += hit->GetLightPMQ();
       adcSum_ZNC  += hit->GetLightPMC() + hit->GetLightPMQ();
       //
       FillHitsData(11,hit->GetLightPMC());
    }
    //
    FillHitsData(0,adcSum_ZNC);
    FillHitsData(1,adcSum_ZNA);
    FillHitsData(2,adcSum_ZPC);
    FillHitsData(3,adcSum_ZPA);
    //
    FillHitsData(4,adcSumQ_ZNC);
    FillHitsData(5,adcSumQ_ZNA);
    FillHitsData(6,adcSumQ_ZPC);
    FillHitsData(7,adcSumQ_ZPA);
  }
}

//___________________________________________________________________________
void AliZDCQADataMakerSim::MakeHits(TTree * hitTree)
{
  // make QA data from Hit Tree
  if(!hitTree){
    AliError("Can't get ZDC hit tree!!");
    return; 
  }	

  TBranch * branch = hitTree->GetBranch("ZDC") ;

  if(!branch){
    AliError("ZDC branch in Hit Tree not found!"); 
    return;
  } 
  
  if(fHitsArray) fHitsArray->Clear() ; 
  else fHitsArray = new TClonesArray("AliZDCHit", 1000);
 
  branch->SetAddress(&fHitsArray) ;
  for (Int_t ientry = 0 ; ientry < branch->GetEntries() ; ientry++) {
    branch->GetEntry(ientry) ;
    MakeHits() ; 
    fHitsArray->Clear() ; 
  }   
  //
  IncEvCountCycleHits();
  IncEvCountTotalHits();
  //
}

//___________________________________________________________________________
void AliZDCQADataMakerSim::MakeDigits(TTree *digitTree)
{
  // makes data from Digit Tree
  if( !GetDigitsData(0) ) InitDigits();

  if(!digitTree){
    AliError("Can't get ZDC digit tree!!");
    return; 
  }	
   
  TBranch * branch = digitTree->GetBranch("ZDC");
  if(!branch){
    AliError("ZDC branch in digit tree not found"); 
    return;
  } 
    
  AliZDCDigit *digit = 0x0;
  branch->SetAddress(&digit);
     
  Float_t adcSum_ZNC=0., adcSum_ZNA=0., adcSum_ZPC=0., adcSum_ZPA=0.;
  Float_t adcSumQ_ZNC=0., adcSumQ_ZNA=0., adcSumQ_ZPC=0., adcSumQ_ZPA=0.;
  //  Float_t adcSum_ZNC_lg=0., adcSum_ZNA_lg=0., adcSum_ZPC_lg=0., adcSum_ZPA_lg=0.;
  //  Float_t adcSumQ_ZNC_lg=0., adcSumQ_ZNA_lg=0., adcSumQ_ZPC_lg=0., adcSumQ_ZPA_lg=0.;
  
  Int_t ndig = digitTree->GetEntries();
  for(Int_t i=0; i<ndig; i++){
      branch->GetEntry(i);
      
      if(digit->GetSector(0)==1 && digit->GetSector(1)!=5){
	  adcSum_ZNC += digit->GetADCValue(0);
	  //adcSum_ZNC_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZNC += digit->GetADCValue(0);
	      //adcSumQ_ZNC_lg+= digit->GetADCValue(1);
	  }
	  else{
	      FillDigitsData(8,digit->GetADCValue(0));
	      //FillDigitsData(20,digit->GetADCValue(1));
	  }
      }
      else if(digit->GetSector(0)==2){
	  adcSum_ZPC += digit->GetADCValue(0);
	  //adcSum_ZPC_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZPC += digit->GetADCValue(0);
	      //adcSumQ_ZPC_lg+= digit->GetADCValue(1);
	  }
	  else{
	      FillDigitsData(10,digit->GetADCValue(0));
	      //FillDigitsData(22,digit->GetADCValue(1));
	  }
      }
      else if(digit->GetSector(0)==4 && digit->GetSector(1)!=5){
	  adcSum_ZNA += digit->GetADCValue(0);
	  //adcSum_ZNA_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZNA += digit->GetADCValue(0);
	      //adcSumQ_ZNA_lg+= digit->GetADCValue(1);
	  }
	  else{
	      FillDigitsData(9,digit->GetADCValue(0));
	      //FillDigitsData(21,digit->GetADCValue(1));
	  }
      }
      else if(digit->GetSector(0)==5){
	  adcSum_ZPA += digit->GetADCValue(0);
	  //adcSum_ZPA_lg += digit->GetADCValue(1);
	  //
	  if(digit->GetSector(1)!=0){
	      adcSumQ_ZPA += digit->GetADCValue(0);
	      //adcSumQ_ZPA_lg+= digit->GetADCValue(1);
	  }
	  else{
	      FillDigitsData(11,digit->GetADCValue(0));
	      //FillDigitsData(23,digit->GetADCValue(1));
	  }
      }
  }
  //
  FillDigitsData(0,adcSum_ZNC);
  FillDigitsData(1,adcSum_ZNA);
  FillDigitsData(2,adcSum_ZPC);
  FillDigitsData(3,adcSum_ZPA);
  //
  FillDigitsData(4,adcSumQ_ZNC);
  FillDigitsData(5,adcSumQ_ZNA);
  FillDigitsData(6,adcSumQ_ZPC);
  FillDigitsData(7,adcSumQ_ZPA);
  //
  IncEvCountCycleDigits();
  IncEvCountTotalDigits();
  //
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
  ResetEventTrigClasses();
  AliQAChecker::Instance()->Run(AliQAv1::kZDC, task, list);  
}
