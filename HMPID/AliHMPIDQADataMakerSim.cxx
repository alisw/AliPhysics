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
#include <TProfile.h>
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
ClassImp(AliHMPIDQADataMakerSim)
           
//____________________________________________________________________________ 
  AliHMPIDQADataMakerSim::AliHMPIDQADataMakerSim() : 
  AliQADataMakerSim(AliQA::GetDetName(AliQA::kHMPID), "HMPID Quality Assurance Data Maker")
{
  // ctor
}

//____________________________________________________________________________ 
AliHMPIDQADataMakerSim::AliHMPIDQADataMakerSim(const AliHMPIDQADataMakerSim& qadm) :
  AliQADataMakerSim() 
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
     TH1F *hHitQdc=new TH1F("HitQdc","HMPID Hit Qdc all chamber;QDC",500,0,4000);
     Add2HitsList(hHitQdc,0);
     TH2F *hHitMap[7];
     for(Int_t iCh=0;iCh<7;iCh++) {
     hHitMap[iCh]=new TH2F(Form("HMPID HitMap%i",iCh),Form("Ch%i;x_{Hit};y_{Hit}",iCh),162,-1,161,146,-1,145);   
    Add2HitsList(hHitMap[iCh],iCh+1);
    }

}

//____________________________________________________________________________ 
void AliHMPIDQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
      TH1F *hDigPcEvt = new TH1F("hDigPcEvt","PC occupancy",156,-1,77);
      TH1F *hDigQ     = new TH1F("Q        ","Charge of digits (ADC)     ",3000,0,3000);
      TH1F *hDigChEvt = new TH1F("hDigChEvt","Chamber occupancy per event",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);

      TProfile *tDigHighQ = new TProfile("tDigHighQ","Highest charge in chamber  ",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);
      TProfile *tDigChEvt = new TProfile("tDigChEvt","Chamber occupancy per event (profile)",AliHMPIDParam::kMaxCh+1,AliHMPIDParam::kMinCh,AliHMPIDParam::kMaxCh+1);

Add2DigitsList(hDigPcEvt,0);
Add2DigitsList(hDigQ    ,1);
Add2DigitsList(hDigChEvt,2);
Add2DigitsList(tDigHighQ,3);
Add2DigitsList(tDigChEvt,4);
}

//____________________________________________________________________________ 
void AliHMPIDQADataMakerSim::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
   TH1F   *hSDigits     = new TH1F("hHmpidSDigits",    "SDigits Q  distribution in HMPID",  500, 0., 5000.) ; 

Add2SDigitsList(hSDigits,0);
}

//____________________________________________________________________________ 

void AliHMPIDQADataMakerSim::MakeHits(TClonesArray * data)
{
 //
 //filling QA histos for Hits
 //
  TClonesArray * hits = dynamic_cast<TClonesArray *>(data) ; 
  if (!hits){
    AliError("Wrong type of hits container") ; 
  } else {
    TIter next(hits); 
    AliHMPIDHit * hit ; 
    while ( (hit = dynamic_cast<AliHMPIDHit *>(next())) ) {
      if(hit->Pid()<500000) GetHitsData(0)->Fill(hit->Q()) ;
      if(hit->Pid()<500000) GetHitsData(hit->Ch()+1)->Fill(hit->LorsX(),hit->LorsY());
    }
  } 

}
//___________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeHits(TTree * data)
{
//
//Opening of the Hit TTree 
//
 TClonesArray *pHits=new TClonesArray("AliHMPIDHit");  data->SetBranchAddress("HMPID",&pHits);
  for(Int_t iEnt=0;iEnt<data->GetEntriesFast();iEnt++){//entries loop
    data->GetEntry(iEnt);
    MakeHits(pHits);
  }//entries loop
}
//____________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeDigits(TClonesArray * data)
{
 //
 //filling QA histos for Digits
 //
  TObjArray *chamber = dynamic_cast<TObjArray*>(data);
  if ( !chamber) {
    AliError("Wrong type of digits container") ; 
  } else {
    for(Int_t i =0; i< chamber->GetEntries(); i++)
      {
	TClonesArray * digits = dynamic_cast<TClonesArray*>(chamber->At(i)); 
	GetDigitsData(2)->Fill(i,digits->GetEntriesFast()/(48.*80.*6.));
        GetDigitsData(4)->Fill(i,digits->GetEntriesFast()/(48.*80.*6.));
        Double_t highQ=0;
	TIter next(digits); 
	AliHMPIDDigit * digit; 
	while ( (digit = dynamic_cast<AliHMPIDDigit *>(next())) ) {
	  GetDigitsData(0)->Fill(10.*i+digit->Pc(),1./(48.*80.));
	  GetDigitsData(1)->Fill(digit->Q());
          if(digit->Q()>highQ) highQ = digit->Q();
	}  
      GetDigitsData(3)->Fill(i,highQ);
 
      }
  }
}
//___________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeDigits(TTree * data)
{
//
//Opening the Digit Tree
//
 TObjArray *pObjDig=new TObjArray(AliHMPIDParam::kMaxCh+1);
  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    TClonesArray *pCA=new TClonesArray("AliHMPIDDigit");
    pObjDig->AddAt(pCA,iCh);
  }

  pObjDig->SetOwner(kTRUE);

  for(Int_t iCh=AliHMPIDParam::kMinCh;iCh<=AliHMPIDParam::kMaxCh;iCh++){
    data->SetBranchAddress(Form("HMPID%i",iCh),&(*pObjDig)[iCh]);
  }
  data->GetEntry(0);

   MakeDigits((TClonesArray *)pObjDig);
}
//____________________________________________________________________________

void AliHMPIDQADataMakerSim::MakeSDigits(TClonesArray * data)
{
 //
 //filling QA histos for SDigits
 //
  TClonesArray * sdigits = dynamic_cast<TClonesArray *>(data) ; 
  if (!sdigits) {
    AliError("Wrong type of sdigits container") ; 
  } else {
    TIter next(sdigits) ; 
    AliHMPIDDigit * sdigit ; 
    while ( (sdigit = dynamic_cast<AliHMPIDDigit *>(next())) ) {
	    GetSDigitsData(0)->Fill(sdigit->Q());
    } 
  }
}
//___________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeSDigits(TTree * data)
{
 //
 // Opening the SDigit Tree
 //
 TClonesArray * sdigits = new TClonesArray("AliHMPIDDigit", 1000) ;

  TBranch * branch = data->GetBranch("HMPID") ;
  if ( ! branch ) {
    AliError("HMPID SDigit Tree not found") ;
    return;
  }
  branch->SetAddress(&sdigits) ;
  branch->GetEntry(0) ;
  MakeSDigits(sdigits) ;
}
//____________________________________________________________________________
void AliHMPIDQADataMakerSim::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  
}

void AliHMPIDQADataMakerSim::EndOfDetectorCycle(AliQA::TASKINDEX, TObjArray *)
{
  //Detector specific actions at end of cycle
  // do the QA checking
//  AliQAChecker::Instance()->Run(AliQA::kHMPID, task, obj) ;  
}

