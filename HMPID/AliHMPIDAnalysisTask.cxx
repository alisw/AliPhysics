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


#include "TH1.h"
#include "TH2.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliHMPIDAnalysisTask.h"

ClassImp(AliHMPIDAnalysisTask)

//__________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask() :
  fESD(0x0),fHmpHistList(0x0),
  fNevts(0),
  fTrigNevts(0),
  fTrigger(0),
  fHmpInner(0x0),fHmpPesdPhmp(0x0),fHmpCkovPesd(0x0),fHmpCkovPhmp(0x0),
  fHmpMipTrkDist(0x0),fHmpMipTrkDistX(0x0),fHmpMipTrkDistY(0x0),fHmpMipCharge3cm(0x0),fHmpMipCharge1cm(0x0),fHmpNumPhots(0x0),
  fHmpTrkFlags(0x0)
{
  //
  //Default ctor
  //
  
}
//___________________________________________________________________________
AliHMPIDAnalysisTask::AliHMPIDAnalysisTask(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fESD(0),fHmpHistList(0x0), fNevts(0),
  fTrigNevts(0),
  fTrigger(0),
  fHmpInner(0x0),fHmpPesdPhmp(0x0),fHmpCkovPesd(0x0),fHmpCkovPhmp(0x0),
  fHmpMipTrkDist(0x0),fHmpMipTrkDistX(0x0),fHmpMipTrkDistY(0x0),fHmpMipCharge3cm(0x0),fHmpMipCharge1cm(0x0),fHmpNumPhots(0x0),
  fHmpTrkFlags(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  
  DefineOutput(0,TList::Class());
 }
//___________________________________________________________________________
AliHMPIDAnalysisTask::~AliHMPIDAnalysisTask() {
  //
  //destructor
  //
  Info("~AliHMPIDAnalysisask","Calling Destructor");
  if (fHmpHistList) {fHmpHistList->Clear(); delete fHmpHistList;}
}
//___________________________________________________________________________
void AliHMPIDAnalysisTask::ConnectInputData(Option_t *) 
{
  // Connect ESD here
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
      AliDebug(2,Form("ERROR: Could not get ESDInputHandler")); 
    } else
      fESD = esdH->GetEvent();
 }
//_________________________________________________
void AliHMPIDAnalysisTask::Exec(Option_t *)
{
  //
  // Main loop function, executed on Event basis
  //
  for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {

    AliESDtrack* track = fESD->GetTrack(iTrack);
    if(!track) continue;
    Double_t rin[3], rout[3];
    track->GetInnerXYZ(rin);
    track->GetOuterXYZ(rout);
        /*
     ((TH2F *)fHmpHistList->At(22))->Fill(rin[0],rin[1]);
     ((TH2F *)fHmpHistList->At(23))->Fill(rout[0],rout[1]);
*/

     if(track->GetHMPIDsignal() == -20) fHmpTrkFlags->Fill(0);
     else if(track->GetHMPIDsignal() == -9) fHmpTrkFlags->Fill(1);
     else if(track->GetHMPIDsignal() == -5) fHmpTrkFlags->Fill(2);
     else if(track->GetHMPIDsignal() == -11) fHmpTrkFlags->Fill(3);
     else if(track->GetHMPIDsignal() == -22) fHmpTrkFlags->Fill(4);
     else fHmpTrkFlags->Fill(4);
  
    if(track->GetHMPIDsignal()== -20) continue;
    if(track->GetHMPIDcluIdx() < 0) continue;
   
    Int_t ch = track->GetHMPIDcluIdx()/1000000; 
    Float_t x, y; Int_t q, nph; 
    Float_t xpc, ypc, th, ph;
    
    track->GetHMPIDmip(x,y,q,nph);
    track->GetHMPIDtrk(xpc,ypc,th,ph);
     
    if(x==0 && y==0 && xpc == 0 && ypc == 0) continue;
    //Printf("%s  %s Good track is found",(char*)__FILE__,__LINE__);
     
    Double_t dist = TMath::Sqrt( (xpc-x)*(xpc-x) + (ypc - y)*(ypc - y));
    fHmpMipTrkDist->Fill(dist);
    fHmpMipTrkDistX->Fill(xpc-x);
    fHmpMipTrkDistY->Fill(ypc - y);
    if(dist<=3.0) fHmpMipCharge3cm->Fill(q);
    
    if(track->GetHMPIDsignal() > 0 )
    {
      Printf("EvtNumInFile: %d HMPID ThetaC: %lf x: %lf xpc: %lf y: %lf ypx: %lf Q: %d nPh: %d IncTheta; %lf IncPhi: %lf",fESD->GetEventNumberInFile(),track->GetHMPIDsignal(),x,xpc,y,ypc,q,nph,th,ph);
      Double_t pHmp[3] = {0},pmod = 0;if(track->GetOuterHmpPxPyPz(pHmp))  pmod = TMath::Sqrt(pHmp[0]*pHmp[0]+pHmp[1]*pHmp[1]+pHmp[2]*pHmp[2]); 
      fHmpPesdPhmp->Fill(track->P(),pmod);  
      if(dist<=1.0)  fHmpMipCharge1cm->Fill(q);
      fHmpNumPhots->Fill(nph);
      fHmpCkovPesd->Fill(track->P(),track->GetHMPIDsignal());
      fHmpCkovPesd->Fill(pmod,track->GetHMPIDsignal());
    }//there is signal
    
    }//track loop
  
  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(0,fHmpHistList) ;
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
void AliHMPIDAnalysisTask::CreateOutputObjects() {
  //
  //HERE ONE CAN CREATE OUTPUT OBJECTS
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //

  //slot #1
   OpenFile(0);
   fHmpHistList = new TList();
   fHmpInner =new TH2F("fHmpInner","HMPID: Inner track XY;X (cm);Y(cm)",800,-400,400,800,-400,400);
   fHmpHistList->Add(fHmpInner);
   
   fHmpPesdPhmp = new TH2F("fHmpPesdPhmp","HMPID: ESD p - running p;HMP p (GeV/c);ESD p (Gev/c)",100,0,10,100,0,10);
   fHmpHistList->Add(fHmpPesdPhmp);
   
   fHmpCkovPesd = new TH2F("fHmpCkovPesd","HMPID: ThetaCherenkov vs P;p_esd (GeV/c);#Theta_C;Entries",100,0,10,110,0,1.1);
   fHmpHistList->Add(fHmpCkovPesd);
   
   fHmpCkovPhmp = new TH2F("fHmpCkovPhmp","HMPID: ThetaCherenkov vs P;p_hmp (GeV/c);#Theta_C;Entries",100,0,10,110,0,1.1);
   fHmpHistList->Add(fHmpCkovPhmp);  
   
   fHmpMipTrkDist = new TH1F("fHmpMipTrkDist","HMPID MIP-Track distance;distance (cm);Entries",800,-20,20);
   fHmpHistList->Add(fHmpMipTrkDist);
   fHmpMipTrkDistX = new TH1F("fHmpMipTrkDistX","HMPID MIP-Track distance in local X;distance (cm);Entries",800,-20,20);
   fHmpHistList->Add(fHmpMipTrkDistX);
   fHmpMipTrkDistY = new TH1F("fHmpMipTrkDistY","HMPID MIP-Track distance in local Y;distance (cm);Entries",800,-20,20);
   fHmpHistList->Add(fHmpMipTrkDistY);
   
   fHmpMipCharge3cm = new TH1F("fHmpMipCharge3cm","HMPID MIP Charge;MIP Charge (ADC);Entries",5001,-0.5,5000.5);
   fHmpHistList->Add(fHmpMipCharge3cm);
   
   fHmpMipCharge1cm = new TH1F("fHmpMipCharge1cm","HMPID MIP Charge;MIP Charge (ADC);Entries",5001,-0.5,5000.5);
   fHmpHistList->Add(fHmpMipCharge1cm);
   
   fHmpNumPhots = new TH1F("fHmpNumPhots","HMPID Number of photon clusters on ring;#photon clus.;Entries",51,-0.5,50.5);
   fHmpHistList->Add(fHmpNumPhots);
   
   fHmpTrkFlags = new TH1F("fHmpTrkFlags","HMPID track flags",6,0,6);
   TString summary[6] =  {"NotPerformed","MipDistCut", "MipQdcCut", "NoPhotAccept", "kNoRad", "other"};
   for(Int_t ibin = 0; ibin < 6; ibin++) fHmpTrkFlags->GetXaxis()->SetBinLabel(ibin+1,Form("%i  %s",ibin+1,summary[ibin].Data()));
   fHmpHistList->Add(fHmpTrkFlags);
   /*
   //0
   TH1F *trkH = new TH1F("trkH","signal flags in HMPID",6,0,6);
   TString summary[6] =  {"NotPerformed","MipDistCut", "MipQdcCut", "NoPhotAccept", "kNoRad", "other"};
   for(Int_t ibin = 0; ibin < 6; ibin++) trkH->GetXaxis()->SetBinLabel(ibin+1,Form("%i  %s",ibin+1,summary[ibin].Data()));
   fHmpHistList->Add(trkH);

  
   TH2F *mod[7], *dq[7];
   TH1F *q[7];

   //1-7
   for(Int_t i=0; i< 7 ; i++) {// to preserve the histo sorting
   mod[i] = new TH2F(Form("mod%i",i),Form("MIP position in chamber %i",i),180,0,180,180,0,180);
   mod[i]->SetMarkerStyle(8);
   mod[i]->SetMarkerColor(2);
   fHmpHistList->Add(mod[i]);
   }
   //8-14
   for(Int_t i=0; i< 7 ; i++) {//to reserve the histo sorting
   q[i] = new TH1F(Form("q%i",i),Form("MIP charge distribution in chamber %i",i),5000,0,5000);
   q[i]->SetMarkerColor(2);
   fHmpHistList->Add(q[i]);
   }
   //15-21
   for(Int_t i=0; i< 7 ; i++) {//to reserve the histo sorting
   dq[i] = new TH2F(Form("dq%i",i),Form("#Delta(mip-track) vs Q_{mip} in chamber %i",i),1000,0,100,5000,0,5000),
   dq[i]->SetMarkerStyle(6);
   fHmpHistList->Add(dq[i]);
   }
   //22
   TH2F *inner = new TH2F("inner","inner track XY",800,-400,400,800,-400,400);
   inner->SetMarkerStyle(6);
   inner->SetXTitle("X cm");
   inner->SetYTitle("Y cm");
   fHmpHistList->Add(inner); 
   //23
   TH2F *outer = new TH2F("outer","outer track XY",800,-400,400,800,-400,400);
   outer->SetMarkerStyle(6);
   outer->SetXTitle("X cm");
   outer->SetYTitle("Y cm");
   fHmpHistList->Add(outer);
*/
}

#endif
