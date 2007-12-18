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

//---
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  A. Mastroserio
//---

// --- ROOT system ---
#include <TClonesArray.h>
#include <TFile.h> 
#include <TH1F.h> 
#include <TH2F.h>
#include <TH1I.h> 
#include <TDirectory.h>
#include <Riostream.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliHMPIDDigit.h"
#include "AliHMPIDHit.h"
#include "AliHMPIDQADataMakerSim.h"

ClassImp(AliHMPIDQADataMakerSim)
           
//____________________________________________________________________________ 
  AliHMPIDQADataMakerSim::AliHMPIDQADataMakerSim() : 
  AliQADataMakerSim(AliQA::GetDetName(AliQA::kHMPID), "HMPID Quality Assurance Data Maker"),
  fhHitQdc(0x0), 
  fhSDigits(0x0),
  fhDigPcEvt(0x0),
  fhDigChEvt(0x0),
  fhDigQ(0x0)
 // fhCluEvt(0x0),
//  fhCluChi2(0x0),
//  fhCluQ(0x0),
//  fhCluFlg(0x0), 
//  fhCluSize(0x0),  
//  fhMipCluSize(0x0),
//  fhCkovP(0x0),
//  fhSigP(0x0),
//  fhMipXY(0x0),
//  fhDifXY(0x0)
{
  // ctor
  for(Int_t i=0; i<7; i++) fhHitMap[i]=0x0;
//  for(Int_t j=0; j<5; j++) fhPid[j]=0x0;
//   fDetectorDir = fOutput->GetDirectory(GetName()) ;  
//   if (!fDetectorDir) 
//     fDetectorDir = fOutput->mkdir(GetName()) ;  
}

//____________________________________________________________________________ 
AliHMPIDQADataMakerSim::AliHMPIDQADataMakerSim(const AliHMPIDQADataMakerSim& qadm) :
  AliQADataMakerSim(), 
  fhHitQdc(qadm.fhHitQdc), 
  fhSDigits(qadm.fhSDigits),
  fhDigPcEvt(qadm.fhDigPcEvt),
  fhDigChEvt(qadm.fhDigChEvt),
  fhDigQ(qadm.fhDigQ)
//  fhCluEvt(qadm.fhCluEvt),
//  fhCluChi2(qadm.fhCluChi2),
//  fhCluQ(qadm.fhCluQ),
//  fhCluFlg(qadm.fhCluFlg),
//  fhCluSize(qadm.fhCluSize),
//  fhMipCluSize(qadm.fhMipCluSize),
//  fhCkovP(qadm.fhCkovP),
//  fhSigP(qadm.fhSigP),
//  fhMipXY(qadm.fhMipXY),
//  fhDifXY(qadm.fhDifXY)
{
  //copy ctor 
  for(Int_t i=0; i<7; i++) fhHitMap[i]=qadm.fhHitMap[i];
//  for(Int_t j=0; j<5; j++) fhPid[j]=qadm.fhPid[j];

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
     fhHitQdc=new TH1F("HitQdc","HMPID Hit Qdc all chamber;QDC",500,0,4000);
     for(Int_t iCh=0;iCh<7;iCh++) 
		 fhHitMap[iCh]=new TH2F(Form("HMPID HitMap%i",iCh),Form("Ch%i;x_{Hit};y_{Hit}",iCh),162,-1,161,146,-1,145);   
}

//____________________________________________________________________________ 
void AliHMPIDQADataMakerSim::InitDigits()
{
  // create Digits histograms in Digits subdir
      fhDigPcEvt=new TH1F("hDigPcEvt","PC occupancy",156,-1,77);
      fhDigChEvt=new TH1F("hDigChEvt","Chamber occupancy",32,-1,7);
      fhDigQ  =new TH1F("Q        "," digit charge               ",3000,0,3000);
}

//____________________________________________________________________________ 
void AliHMPIDQADataMakerSim::InitSDigits()
{
  // create SDigits histograms in SDigits subdir
      fhSDigits     = new TH1F("hHmpidSDigits",    "SDigits Q  distribution in HMPID",  500, 0., 5000.) ; 
}

//____________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeHits(TObject * data)
{
  //fills QA histos for Hits
  TClonesArray * hits = dynamic_cast<TClonesArray *>(data) ; 
  if (!hits){
    AliError("Wrong type of hits container") ; 
  } else {
    TIter next(hits); 
    AliHMPIDHit * hit ; 
    while ( (hit = dynamic_cast<AliHMPIDHit *>(next())) ) {
      if(hit->Pid()<500000) fhHitQdc->Fill(hit->Q()) ;
      if(hit->Pid()<500000) fhHitMap[hit->Ch()]->Fill(hit->LorsX(),hit->LorsY());
    }
  } 
}

//____________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeDigits( TObject * data)
{
  //fills QA histos for Digits
  TObjArray *chambers = dynamic_cast<TObjArray*>(data);
  if ( !chambers) {
    AliError("Wrong type of digits container") ; 
  } else {
    for(Int_t i =0; i< chambers->GetEntries(); i++)
      {
	TClonesArray * digits = dynamic_cast<TClonesArray*>(chambers->At(i)); 
	fhDigChEvt->Fill(i,digits->GetEntriesFast()/(48.*80.*6.));
	TIter next(digits); 
	AliHMPIDDigit * digit; 
	while ( (digit = dynamic_cast<AliHMPIDDigit *>(next())) ) {
	  fhDigPcEvt->Fill(10.*i+digit->Pc(),1./(48.*80.));
	  fhDigQ->Fill(digit->Q());
	}  
      }
  }
}

//____________________________________________________________________________
void AliHMPIDQADataMakerSim::MakeSDigits( TObject * data)
{
  //fills QA histos for SDigits
  TClonesArray * sdigits = dynamic_cast<TClonesArray *>(data) ; 
  if (!sdigits) {
    AliError("Wrong type of sdigits container") ; 
  } else {
    AliHMPIDDigit *ref = (AliHMPIDDigit *)sdigits->At(0);
    Float_t zero = ref->GetTrack(0); 
    TIter next(sdigits) ; 
    AliHMPIDDigit * sdigit ; 
    while ( (sdigit = dynamic_cast<AliHMPIDDigit *>(next())) ) {
      fhSDigits->Fill(sdigit->Q()) ;
      if(zero == sdigit->GetTrack(0)) continue;
      else zero = sdigit->GetTrack(0);
    } 
  }
}
