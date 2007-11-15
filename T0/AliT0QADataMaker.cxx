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
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliT0digit.h"
#include "AliT0hit.h"
#include "AliT0RecPoint.h"
#include "AliT0QADataMaker.h"

ClassImp(AliT0QADataMaker)
           
//____________________________________________________________________________ 
  AliT0QADataMaker::AliT0QADataMaker() : 
  AliQADataMaker(AliQA::GetDetName(AliQA::kT0), "T0 Quality Assurance Data Maker"),
  //  fhHitsTime(0x0),
  // fhHitsAmp(0x0),
  fhHitsEff(0x0),
  //  fhDigCFD(0x0),
  //  fhDigLEDamp(0x0),
  //  fhDigQTC(0x0),
  fhDigMean(0x0),
  fhDigEff(0x0),
  //  fhRecCFD(0x0),
  //  fhRecLEDamp(0x0),
  //  fhRecQTC(0x0),
  fhRecMean(0x0),
  fhRecEff(0x0),
  fhESDMean(0x0),
  fhESDVertex(0x0)

{
  // ctor
  for(Int_t i=0; i<24; i++) {
    fhHitsTime[i]=0x0;
   fhDigCFD[i]=0x0;
    fhDigLEDamp[i]=0x0;
    fhRecCFD[i]=0x0;
    fhRecLEDamp[i]=0x0;
    fhRecQTC[i]=0x0;
  }
 //   fDetectorDir = fOutput->GetDirectory(GetName()) ;  
//   if (!fDetectorDir) 
//     fDetectorDir = fOutput->mkdir(GetName()) ;  
}

//____________________________________________________________________________ 
AliT0QADataMaker::AliT0QADataMaker(const AliT0QADataMaker& qadm) :
  AliQADataMaker(), 
 //  fhHitsTime(0x0),
  fhHitsEff(0x0),
  //  fhDigCFD(0x0),
  //  fhDigLEDamp(0x0),
  //  fhDigQTC(0x0),
  fhDigMean(0x0),
  fhDigEff(0x0),
  //  fhRecCFD(0x0),
  //  fhRecLEDamp(0x0),
  //  fhRecQTC(0x0),
  fhRecMean(0x0),
  fhRecEff(0x0),
  fhESDMean(0x0),
  fhESDVertex(0x0)
{
  //copy ctor 
  for(Int_t i=0; i<24; i++) {
    fhHitsTime[i]=0x0;
    fhDigCFD[i]=0x0;
    fhDigLEDamp[i]=0x0;
    fhRecCFD[i]=0x0;
    fhRecLEDamp[i]=0x0;
    fhRecQTC[i]=0x0;
  }

  SetName((const char*)qadm.GetName()) ; 
  SetTitle((const char*)qadm.GetTitle()); 
}

//__________________________________________________________________
AliT0QADataMaker& AliT0QADataMaker::operator = (const AliT0QADataMaker& qadm )
{
  // Equal operator.
  this->~AliT0QADataMaker();
  new(this) AliT0QADataMaker(qadm);
  return *this;
}
 
//____________________________________________________________________________ 
void AliT0QADataMaker::InitHits()
{
  // create Hits histograms in Hits subdir
  TString timename;
  timename ="hHitTime";
  for (Int_t i=0; i<24; i++)
    {
      timename += i;
       if(i<12)  fhHitsTime[i] = new TH1F(timename.Data(),timename.Data(),100,2,3);
     else  
       fhHitsTime[i] = new TH1F(timename.Data(),timename.Data(),100,12,13);
    }
  fhHitsEff = new TH1F("hHitsEff", "Hits Efficiency", 25,-0.5,24.5);
}

//____________________________________________________________________________ 
void AliT0QADataMaker::InitDigits()
{
  // create Digits histograms in Digits subdir
  TString timename, ampname, qtcname;
  timename ="hDigCFD";
  ampname = "hDigLED";
  qtcname = "hDigQTC";
  for (Int_t i=0; i<24; i++)
    {
      timename += i;
      ampname += i;
      qtcname += i;
      fhDigCFD[i] = new TH1F(timename.Data(), timename.Data(),100,100,1000);
      fhDigLEDamp[i] = new TH1F(ampname.Data(), ampname.Data(),100,100,1000);
      fhDigQTC[i] = new TH1F(qtcname.Data(), qtcname.Data(),100,100,1000);
    }
  fhDigEff = new TH1F("hDigEff","digits efficiency", 25,-0.5,24.5);
  fhDigMean = new TH1F("hDigMean","online mean signal", 100,500,600);

}

//____________________________________________________________________________ 

void AliT0QADataMaker::InitRecPoints()
{
  // create cluster histograms in RecPoint subdir
 
  // create Digits histograms in Digits subdir
  TString timename,ampname, qtcname;
  timename ="hRecCFD";
  ampname = "hRecLED";
  qtcname = "hRecQTC";
  for (Int_t i=0; i<24; i++)
    {
      timename += i;
      ampname += i;
      qtcname += i;
      fhRecCFD[i] = new TH1F(timename.Data(), timename.Data(),100,100,1000);
      fhRecLEDamp[i] = new TH1F(ampname.Data(), ampname.Data(),100,100,1000);
      fhRecQTC[i] = new TH1F(qtcname.Data(), qtcname.Data(),100,100,1000);
    }
  fhRecEff = new TH1F("hRecEff","Efficiency rec.points",25,-0.5,24.5);
  fhRecMean = new TH1F("hRecMean"," reconstructed mean signal",100,500,600);

}
//____________________________________________________________________________
void AliT0QADataMaker::InitESDs()
{
  //create ESDs histograms in ESDs subdir
  fhESDMean = new TH1F("hESDmean"," ESD mean",100,500,600);
  fhESDVertex = new TH1F("hESDvertex","EAD vertex",100,-50,50);


}
//____________________________________________________________________________
void AliT0QADataMaker::MakeHits(TObject * data)
{
  //fills QA histos for Hits
  TClonesArray * hits = dynamic_cast<TClonesArray *>(data) ; 
  if (!hits){
    AliError("Wrong type of hits container") ; 
  } else {
    TIter next(hits); 
    AliT0hit * hit ; 
    while ( (hit = dynamic_cast<AliT0hit *>(next())) ) {
      Int_t pmt=hit->Pmt();
      fhHitsTime[pmt]->Fill(hit->Time()) ;
       fhHitsEff->Fill(pmt);
   }
  } 
}

//____________________________________________________________________________
void AliT0QADataMaker::MakeDigits( TObject * data)
{
  //fills QA histos for Digits
  TArrayI *digCFD = new TArrayI(24);
  TArrayI *digLED = new TArrayI(24);
  TArrayI *digQT0 = new TArrayI(24);
  TArrayI *digQT1 = new TArrayI(24);
  Int_t refpoint=0;

  TClonesArray * digits = dynamic_cast<TClonesArray *>(data) ; 
  
  if ( !digits) {
    AliError("Wrong type of digits container") ; 
  } else {
    TIter next(digits) ; 
    AliT0digit * digit ; 
    while ( (digit = dynamic_cast<AliT0digit *>(next())) ) {
      digit->GetTimeCFD(*digCFD);
      digit->GetTimeLED(*digLED);
      digit->GetQT0(*digQT0);
      digit->GetQT1(*digQT1);
      refpoint =  digit->RefPoint();
      for (Int_t i=0; i<24; i++)
	{
	  if (digCFD->At(i)>0) {
	    Int_t cfd=digCFD->At(i)- refpoint;
	    fhDigCFD[i] -> Fill(cfd);
	    fhDigLEDamp[i] -> Fill(digLED->At(i) - digCFD->At(i));
	    fhDigQTC [i]-> Fill(digQT1->At(i) - digQT0->At(i));
	    fhDigEff->Fill(i);
	  }
	}  
    }
  }
  delete digCFD;
  delete digLED;
  delete digQT0;
  delete digQT1;

}


//____________________________________________________________________________
void AliT0QADataMaker::MakeRecPoints(TTree * clustersTree)
{
  //fills QA histos for clusters

   AliT0RecPoint* frecpoints= new AliT0RecPoint ();
    if (!frecpoints) {
    AliError("Reconstruct Fill ESD >> no recpoints found");
    return;
  }
  TBranch *brRec =clustersTree ->GetBranch("T0");
  if (brRec) {
    brRec->SetAddress(&frecpoints);
  }else{
    cerr<<"EXEC Branch T0 rec not found"<<endl;
    // exit(111);
    return;
  } 
    
  brRec->GetEntry(0);
  
  for ( Int_t i=0; i<24; i++) {
    fhRecCFD[i] -> Fill(frecpoints -> GetTime(i)); 
    fhRecQTC[i] -> Fill(frecpoints -> GetAmp(i));
    fhRecLEDamp[i] -> Fill( frecpoints->AmpLED(i));
    if(frecpoints -> GetTime(i) > 0) fhRecEff->Fill(i);
  }
  fhRecMean->Fill(frecpoints->GetMeanTime());
  
}

//____________________________________________________________________________
void AliT0QADataMaker::MakeESDs(AliESDEvent * esd)
{
  //fills QA histos for ESD

  fhESDMean -> Fill(esd->GetT0());
  fhESDVertex -> Fill(esd->GetT0zVertex());

}

