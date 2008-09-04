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

#ifndef AliHMPIDAnalysisTASK_CXX
#define AliHMPIDAnalysisTASK_CXX

#include "AliHMPIDAnalysisTask.h"
#include "TCanvas.h"
#include "AliStack.h"
#include "TParticle.h"
#include "Riostream.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH1F.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "TChain.h"
#include "AliESDtrack.h"
#include "AliLog.h"

ClassImp(AliHMPIDAnalysisTask)

//__________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask() :
  fHistList(0x0),
  fHistEventsProcessed(0x0),
  fNevts(0),
  fTrigNevts(0),
  fTrigger(0)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fHistList(0x0),
  fHistEventsProcessed(0x0),
  fNevts(0),
  fTrigNevts(0),
  fTrigger(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliHMPIDAnalysisTask","Calling Constructor");

  /*
    DefineInput(0) and DefineOutput(0)
    are taken care of by AliAnalysisTaskSE constructor
  */
  DefineOutput(1,TH1I::Class());
  DefineOutput(2,TList::Class());
}

//___________________________________________________________________________
AliHMPIDAnalysisTask& AliHMPIDAnalysisTask::operator=(const AliHMPIDAnalysisTask& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fHistList = c.fHistList ;
    fHistEventsProcessed = c.fHistEventsProcessed;
    fNevts = c.fNevts;
    fTrigNevts = c.fTrigNevts;
    fTrigger = c.fTrigger;
  }
  return *this;
}

//___________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask(const AliHMPIDAnalysisTask& c) :
  AliAnalysisTaskSE(c),
  fHistList(c.fHistList),
  fHistEventsProcessed(c.fHistEventsProcessed),
  fNevts(c.fNevts),
  fTrigNevts(c.fTrigNevts),
  fTrigger(c.fTrigger)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliHMPIDAnalysisTask::~AliHMPIDAnalysisTask() {
  //
  //destructor
  //
  Info("~AliHMPIDAnalysisask","Calling Destructor");
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fHistList) {fHistList->Clear(); delete fHistList;}
}

//_________________________________________________
void AliHMPIDAnalysisTask::UserExec(Option_t *)
{
  //
  // Main loop function, executed on Event basis
  //
  Info("UserExec","") ;

  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if (!fESD) {
    Error("UserExec","NO ESD FOUND!");
    return;
  }
  fNevts++;
 
  if(fESD->GetTriggerMask()&fTrigger == fTrigger) fTrigNevts++;


  for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {

    AliESDtrack* track = fESD->GetTrack(iTrack);
    if(!track) continue;
    Double_t rin[3], rout[3];
    track->GetInnerXYZ(rin);
    track->GetOuterXYZ(rout);
  
     ((TH2F *)fHistList->At(22))->Fill(rin[0],rin[1]);
     ((TH2F *)fHistList->At(23))->Fill(rout[0],rout[1]);


     TH1F *h = ((TH1F*)(fHistList->At(0)));
     if(track->GetHMPIDsignal() == -20) h->Fill(0);
     else if(track->GetHMPIDsignal() == -9) h->Fill(1);
     else if(track->GetHMPIDsignal() == -5) h->Fill(2);
     else if(track->GetHMPIDsignal() == -11) h->Fill(3);
     else if(track->GetHMPIDsignal() == -22) h->Fill(4);
     else h->Fill(4);
    
     if(track->GetHMPIDcluIdx() == 99099999) continue;
     if(track->GetHMPIDcluIdx() < 0) continue;

    Int_t ch = track->GetHMPIDcluIdx()/1000000; 

    Float_t x, y; Int_t q, nph; 
    track->GetHMPIDmip(x,y,q,nph);

    if(x==0 && y==0) continue;
    Float_t xpc, ypc, th, ph;
    track->GetHMPIDtrk(xpc,ypc,th,ph); //special setting in local cosmic rec!!!
                                       //do not use with standard rec

    Double_t dist = TMath::Sqrt( (xpc-x)*(xpc-x) + (ypc - y)*(ypc - y));

     if(ch >=0 && ch < 7) {
     if(dist < 3) ((TH2F *)fHistList->At(ch+1))->Fill(x,y);
     ((TH1F *)fHistList->At(ch+7+1))->Fill(q);
     ((TH2F *)fHistList->At(ch+14+1))->Fill(dist,q);
     }
  
    }//track loop

  fHistEventsProcessed->Fill(0);
  
  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHistEventsProcessed) ;
  PostData(2,fHistList) ;
}


//___________________________________________________________________________
void AliHMPIDAnalysisTask::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  Info("Terminate","");
  AliAnalysisTaskSE::Terminate();

}


//___________________________________________________________________________
void AliHMPIDAnalysisTask::UserCreateOutputObjects() {
  //
  //HERE ONE CAN CREATE OUTPUT OBJECTS
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  Info("CreateOutputObjects","CreateOutputObjects of task %s", GetName());

  //slot #1
  OpenFile(1);
  fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;

   OpenFile(2);
   fHistList = new TList();
   
   //0
   TH1F *trkH = new TH1F("trkH","signal flags in HMPID",6,0,6);
   TString summary[6] =  {"NotPerformed","MipDistCut", "MipQdcCut", "NoPhotAccept", "kNoRad", "other"};
   for(Int_t ibin = 0; ibin < 6; ibin++) trkH->GetXaxis()->SetBinLabel(ibin+1,Form("%i  %s",ibin+1,summary[ibin].Data()));
   fHistList->Add(trkH);

  
   TH2F *mod[7], *dq[7];
   TH1F *q[7];

   //1-7
   for(Int_t i=0; i< 7 ; i++) {// to preserve the histo sorting
   mod[i] = new TH2F(Form("mod%i",i),Form("MIP position in chamber %i",i),180,0,180,180,0,180);
   mod[i]->SetMarkerStyle(8);
   mod[i]->SetMarkerColor(2);
   fHistList->Add(mod[i]);
   }
   //8-14
   for(Int_t i=0; i< 7 ; i++) {//to reserve the histo sorting
   q[i] = new TH1F(Form("q%i",i),Form("MIP charge distribution in chamber %i",i),5000,0,5000);
   q[i]->SetMarkerColor(2);
   fHistList->Add(q[i]);
   }
   //15-21
   for(Int_t i=0; i< 7 ; i++) {//to reserve the histo sorting
   dq[i] = new TH2F(Form("dq%i",i),Form("#Delta(mip-track) vs Q_{mip} in chamber %i",i),1000,0,100,5000,0,5000),
   dq[i]->SetMarkerStyle(6);
   fHistList->Add(dq[i]);
   }
   //22
   TH2F *inner = new TH2F("inner","inner track XY",800,-400,400,800,-400,400);
   inner->SetMarkerStyle(6);
   inner->SetXTitle("X cm");
   inner->SetYTitle("Y cm");
   fHistList->Add(inner); 
   //23
   TH2F *outer = new TH2F("outer","outer track XY",800,-400,400,800,-400,400);
   outer->SetMarkerStyle(6);
   outer->SetXTitle("X cm");
   outer->SetYTitle("Y cm");
   fHistList->Add(outer);

}

#endif
