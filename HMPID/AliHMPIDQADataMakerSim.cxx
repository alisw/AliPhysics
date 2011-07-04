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


/* $Id$ */

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h>
#include <Riostream.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliQAChecker.h"
#include "AliLog.h"
#include "AliHMPIDDigit.h"
#include "AliHMPIDHit.h"
#include "AliHMPIDCluster.h"
#include "AliHMPIDQADataMakerSim.h"
#include "AliHMPIDParam.h"
#include "AliHMPIDRawStream.h"
#include "AliLog.h"

//.
// HMPID AliHMPIDQADataMakerSim base class
// for QA of simulation
// here also errors are calculated
//.

ClassImp(AliHMPIDQADataMakerSim)
           
//____________________________________________________________________________ 
  AliHMPIDQADataMakerSim::AliHMPIDQADataMakerSim() : 
  AliQADataMakerSim(AliQAv1::GetDetName(AliQAv1::kHMPID), "HMPID Quality Assurance Data Maker"), fChannel(0)
{
  // ctor
}

//____________________________________________________________________________ 
AliHMPIDQADataMakerSim::AliHMPIDQADataMakerSim(const AliHMPIDQADataMakerSim& qadm) :
  AliQADataMakerSim(),fChannel(0)
{
  //copy ctor 
  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliHMPIDQADataMakerSim& AliHMPIDQADataMakerSim::operator = (const AliHMPIDQADataMakerSim& qadm )
{
  // Equal operator.
  this->~AliHMPIDQADataMakerSim();
  new(this) AliHMPIDQADataMakerSim(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliHMPIDQADataMakerSim::InitHits()
{
  // create Hits histograms in Hits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *hHitQdc=new TH1F("HitQdc","HMPID Hit Qdc all chamber;QDC;Entries",500,0,4000);
  Add2HitsList(hHitQdc,0, !expert, image);
  TH2F *hHitMap[7];
  for(Int_t iCh=0;iCh<7;iCh++) {
    hHitMap[iCh]=new TH2F(Form("HMPID HitMap%i",iCh),Form("Ch%i;x_{Hit};y_{Hit};Entries",iCh),162,-1,161,146,-1,145);   
    Add2HitsList(hHitMap[iCh],iCh+1,expert,!image);
  }
  //
  ClonePerTrigClass(AliQAv1::kHITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliHMPIDQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 
  
  TH1F *hDigChEvt = new TH1F("hDigChEvt","Chamber occupancy per event;Occupanc [%];Entries",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);
  TH1F *hDigPcEvt = new TH1F("hDigPcEvt","PC occupancy",156,-1,77);
  TH2F *hDigMap[7];
  TH1F *hDigQ[42];
  for(Int_t iCh =0; iCh < 7; iCh++){
    hDigMap[iCh] = new TH2F(Form("MapCh%i",iCh),Form("Digit Map in Chamber %i;Digit #;Entries",iCh),159,0,159,143,0,143);
    for(Int_t iPc =0; iPc < 6; iPc++ ){
      hDigQ[iCh*6+iPc] = new TH1F(Form("QCh%iPc%i        ",iCh,iPc),Form("Charge of digits (ADC) in Chamber %i and PC %i;Charge;Entries",iCh,iPc),4100,0,4100);
    }
  }
  
  Add2DigitsList(hDigChEvt,0, !expert, image);
  Add2DigitsList(hDigPcEvt,1,expert, !image);
  for(Int_t iMap=0; iMap < 7; iMap++) Add2DigitsList(hDigMap[iMap],2+iMap,expert, !image);
  for(Int_t iH =0; iH < 42 ; iH++) Add2DigitsList(hDigQ[iH]    ,9+iH,expert,!image);
  //
  ClonePerTrigClass(AliQAv1::kDIGITS); // this should be the last line
}

//____________________________________________________________________________ 
void AliHMPIDQADataMakerSim::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
  const Bool_t expert   = kTRUE ; 
  const Bool_t image    = kTRUE ; 

  TH1F   *hSDigits     = new TH1F("hHmpidSDigits",    "SDigits Q  distribution in HMPID;QDC;Entries",  500, 0., 5000.) ; 
  Add2SDigitsList(hSDigits,0, !expert, image);
  //
  ClonePerTrigClass(AliQAv1::kSDIGITS); // this should be the last line}
}

//____________________________________________________________________________ 

void AliHMPIDQADataMakerSim::MakeHits()
{
 //
 //filling QA histos for Hits
 //
  
  TIter next(fHitsArray); 
  AliHMPIDHit * hit ; 
  while ( (hit = dynamic_cast<AliHMPIDHit *>(next())) ) {
    if(hit->Pid()<500000) FillHitsData(0,hit->Q()) ;
    if(hit->Pid()<500000) FillHitsData(hit->Ch()+1,hit->LorsX(),hit->LorsY());
  }
} 

//___________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeHits(TTree * data)
{
//
//Opening of the Hit TTree 
//
 if (fHitsArray) 
   fHitsArray->Clear() ; 
  else 
    fHitsArray=new TClonesArray("AliHMPIDHit");  
  data->SetBranchAddress("HMPID",&fHitsArray);
  for(Int_t iEnt=0;iEnt<data->GetEntriesFast();iEnt++){//entries loop
    data->GetEntry(iEnt);
    MakeHits();
  }//entries loop
}
//___________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeDigits()
{
  //
  //filling QA histos for Digits
  //
   
  Int_t i = fChannel ; 
  FillDigitsData(0,i,fDigitsArray->GetEntriesFast()/(48.*80.*6.));
  TIter next(fDigitsArray); 
  AliHMPIDDigit * digit; 
  while ( (digit = dynamic_cast<AliHMPIDDigit *>(next())) ) {
    FillDigitsData(1,10.*i+digit->Pc(),1./(48.*80.));
    FillDigitsData(2+i,digit->PadChX(),digit->PadChY());
    FillDigitsData(9+i*6+digit->Pc(),digit->Q());
  }  
}  
//___________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeDigits(TTree * data)
{
  //
  //Opening the Digit Tree
  //
  
  if(fDigitsArray) 
    fDigitsArray->Clear() ; 
  else
    fDigitsArray=new TClonesArray("AliHMPIDDigit");
  
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    fChannel = iCh ; 
    data->SetBranchAddress(Form("HMPID%i",iCh),&fDigitsArray);
    data->GetEntry(0);
    MakeDigits();
    fDigitsArray->Clear() ; 
  }
}

//____________________________________________________________________________

void AliHMPIDQADataMakerSim::MakeSDigits()
{
 //
 //filling QA histos for SDigits
 //
 
  TIter next(fSDigitsArray) ; 
  AliHMPIDDigit * sdigit ; 
  while ( (sdigit = dynamic_cast<AliHMPIDDigit *>(next())) ) {
    FillSDigitsData(0,sdigit->Q());
  } 
}
//___________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeSDigits(TTree * data)
{
 //
 // Opening the SDigit Tree
 //
 if (fSDigitsArray)
   fSDigitsArray->Clear() ; 
  else 
    fSDigitsArray = new TClonesArray("AliHMPIDDigit", 1000) ;

  TBranch * branch = data->GetBranch("HMPID") ;
  if ( ! branch ) {
    AliError("HMPID SDigit Tree not found") ;
    return;
  }
  branch->SetAddress(&fSDigitsArray) ;
  branch->GetEntry(0) ;
  MakeSDigits() ;
}
//____________________________________________________________________________
void AliHMPIDQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

void AliHMPIDQADataMakerSim::EndOfDetectorCycle(AliQAv1::TASKINDEX_t task, TObjArray **obj)
{
  //Detector specific actions at end of cycle
  // do the QA checking
  ResetEventTrigClasses(); // reset triggers list to select all histos
  AliQAChecker::Instance()->Run(AliQAv1::kHMPID, task, obj) ;  
}

