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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber                                                  //
//  Comparison macro for ESD                                                 //
//  responsible: 
//  marian.ivanov@cern.ch                                                    //
//
//

/*
marian.ivanov@cern.ch
Usage:


.L $ALICE_ROOT/STEER/AliGenInfo.C+
//be sure you created genTracks file before
.L $ALICE_ROOT/STEER/AliESDComparisonMI.C+

//
ESDCmpTr *t2 = new ESDCmpTr("genTracks.root","cmpESDTracks.root","galice.root",-1,1,0);
t2->Exec();

//
//some cuts definition
TCut cprim("cprim","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<0.00005&&abs(MC.fVDist[2])<0.0005")
//TCut cprim("cprim","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<0.5&&abs(MC.fVDist[2])<0.5")
TCut citsin("citsin","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<3.9");

TCut csec("csec","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)>0.5");


TCut crec("crec","fReconstructed==1");
TCut cteta1("cteta1","abs(MC.fTrackRef.Theta()/3.1415-0.5)<0.25");
TCut cpos1("cpos1","abs(MC.fParticle.fVz/sqrt(MC.fParticle.fVx*MC.fParticle.fVx+MC.fParticle.fVy*MC.fParticle.fVy))<1");
TCut csens("csens","abs(sqrt(fVDist[0]**2+fVDist[1]**2)-170)<50");
TCut cmuon("cmuon","abs(MC.fParticle.fPdgCode==-13)");
TCut cchi2("cchi2","fESDTrack.fITSchi2MIP[0]<7.&&fESDTrack.fITSchi2MIP[1]<5.&&fESDTrack.fITSchi2MIP[2]<7.&&fESDTrack.fITSchi2MIP[3]<7.5&&fESDTrack.fITSchi2MIP[4]<6.")

AliESDComparisonDraw comp;  
comp.SetIO(); 

//
//example
comp.fTree->SetAlias("radius","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)");

TH1F his("his","his",100,0,20);
TH1F hpools("hpools","hpools",100,-7,7);
TH1F hfake("hfake","hfake",1000,0,150);
TProfile profp0("profp0","profp0",20,-0.4,0.9)

comp.DrawLogXY("fTPCinP0[3]","fTPCDelta[4]/fTPCinP1[3]","fReconstructed==1&&abs(fPdg)==211"+cprim,"1",4,0.2,1.5,-0.06,0.06)
comp.fRes->Draw();
comp.fMean->Draw();  
comp.Eff("fTPCinP0[3]","fRowsWithDigits>120&&abs(fPdg)==211"+cteta1+cpos1+cprim,"fTPCOn",20,0.2,1.5)
comp.fRes->Draw();

comp.fTree->Draw("fESDTrack.fITSsignal/fESDTrack.fTPCsignal","fITSOn&&fTPCOn&&fESDTrack.fITSFakeRatio==0") 

TH1F his("his","his",100,0,20);
TH1F hpools("hpools","hpools",100,-7,7);

TH2F * hdedx0 = new TH2F("dEdx0","dEdx0",100, 0,2,200,0,250); hdedx0->SetMarkerColor(2);
TH2F * hdedx1 = new TH2F("dEdx1","dEdx1",100, 0,2,200,0,250); hdedx1->SetMarkerColor(3);
TH2F * hdedx2 = new TH2F("dEdx2","dEdx2",100, 0,2,200,0,250); hdedx2->SetMarkerColor(4);
TH2F * hdedx3 = new TH2F("dEdx3","dEdx3",100, 0,2,200,0,250); hdedx3->SetMarkerColor(5);
comp.fTree->Draw("fESDTrack.fITSsignal:MC.fParticle.P()>>dEdx0","fITSOn&&MC.fParticle.P()<2&&abs(fPdg)==211&&fITStrack.fN==6"+cprim) 
comp.fTree->Draw("fESDTrack.fITSsignal:MC.fParticle.P()>>dEdx1","fITSOn&&MC.fParticle.P()<2&&abs(fPdg)==2212&&fITStrack.fN==6"+cprim) 
comp.fTree->Draw("fESDTrack.fITSsignal:MC.fParticle.P()>>dEdx2","fITSOn&&MC.fParticle.P()<2&&abs(fPdg)==321&&fITStrack.fN==6"+cprim) 
comp.fTree->Draw("fESDTrack.fITSsignal:MC.fParticle.P()>>dEdx3","fITSOn&&MC.fParticle.P()<2&&abs(fPdg)==11&&fITStrack.fN==6"+cprim) 
hdedx0->Draw(); hdedx1->Draw("same"); hdedx2->Draw("same"); hdedx3->Draw("same");

comp.DrawXY("fITSinP0[3]","fITSPools[4]","fReconstructed==1&&fPdg==-211&&fITSOn"+cprim,"1",4,0.2,1.0,-8,8)

*/


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string.h>
//ROOT includes
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TString.h"
#include "TBenchmark.h"
#include "TStopwatch.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TTimer.h"
#include "TVector3.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TText.h"
#include "Getline.h"
#include "TStyle.h"

//ALIROOT includes
#include "AliRun.h"
#include "AliStack.h"
#include "AliESDtrack.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliTPC.h"
#include "AliTPCLoader.h"
#include "AliDetector.h"
#include "AliTrackReference.h"
#include "AliRun.h"
#include "AliTPCParamSR.h"
#include "AliTracker.h"
#include "AliComplexCluster.h"
#include "AliMagF.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliITStrackV2.h"
#include "AliTRDtrack.h"
#endif
#include "AliGenInfo.h"
#include "AliESDComparisonMI.h"



//
//
void AliESDRecInfo::Update(AliMCInfo* info,AliTPCParam * par, Bool_t reconstructed)
{
  //
  //
  //calculates derived variables
  //  
  //
  AliTrackReference * ref = &(info->fTrackRef);
  fTPCinR0[0] = info->fTrackRef.X();	
  fTPCinR0[1] = info->fTrackRef.Y();	
  fTPCinR0[2] = info->fTrackRef.Z();
  fTPCinR0[3] = TMath::Sqrt(fTPCinR0[0]*fTPCinR0[0]+fTPCinR0[1]*fTPCinR0[1]);
  fTPCinR0[4] = TMath::ATan2(fTPCinR0[1],fTPCinR0[0]);
  //
  fTPCinP0[0] = ref->Px();
  fTPCinP0[1] = ref->Py();
  fTPCinP0[2] = ref->Pz();
  fTPCinP0[3] = ref->Pt();
  fTPCinP0[4] = ref->P();
  fDeltaP     = (ref->P()-info->fParticle.P())/info->fParticle.P();
  //
  //
  if (fTPCinP0[3]>0.0000001){
    //
    fTPCAngle0[0] = TMath::ATan2(fTPCinP0[1],fTPCinP0[0]);
    fTPCAngle0[1] = TMath::ATan(fTPCinP0[2]/fTPCinP0[3]);
  }
  //
  //
  fITSinP0[0]=info->fParticle.Px();
  fITSinP0[1]=info->fParticle.Py();
  fITSinP0[2]=info->fParticle.Pz();
  fITSinP0[3]=info->fParticle.Pt();    
  //
  fITSinR0[0]=info->fParticle.Vx();
  fITSinR0[1]=info->fParticle.Vy();
  fITSinR0[2]=info->fParticle.Vz();
  fITSinR0[3] = TMath::Sqrt(fITSinR0[0]*fITSinR0[0]+fITSinR0[1]*fITSinR0[1]);
  fITSinR0[4] = TMath::ATan2(fITSinR0[1],fITSinR0[0]);
  //
  //
  if (fITSinP0[3]>0.0000001){
    fITSAngle0[0] = TMath::ATan2(fITSinP0[1],fITSinP0[0]);
    fITSAngle0[1] = TMath::ATan(fITSinP0[2]/fITSinP0[3]);
  }
  //
  if (reconstructed==kFALSE){
    fReconstructed = kFALSE;
    fTPCOn = kFALSE;
    fITSOn = kFALSE;
    fTRDOn = kFALSE;
    return;
  }
  fReconstructed = kTRUE;
  fTPCOn = ((fESDTrack.GetStatus()&AliESDtrack::kTPCrefit)>0) ? kTRUE : kFALSE;
  fITSOn = ((fESDTrack.GetStatus()&AliESDtrack::kITSrefit)>0) ? kTRUE : kFALSE;
  fTRDOn = ((fESDTrack.GetStatus()&AliESDtrack::kTRDrefit)>0) ? kTRUE : kFALSE;
  
  if (fTPCOn){
    //TPC
    fESDTrack.GetInnerXYZ(fTPCinR1);
    fTPCinR1[3] = TMath::Sqrt(fTPCinR1[0]*fTPCinR1[0]+fTPCinR1[1]*fTPCinR1[1]);
    fTPCinR1[4] = TMath::ATan2(fTPCinR1[1],fTPCinR1[0]);	
    fESDTrack.GetInnerPxPyPz(fTPCinP1);
    fTPCinP1[3] = TMath::Sqrt(fTPCinP1[0]*fTPCinP1[0]+fTPCinP1[1]*fTPCinP1[1]);
    fTPCinP1[4] = TMath::Sqrt(fTPCinP1[3]*fTPCinP1[3]+fTPCinP1[2]*fTPCinP1[2]);
    //
    //
    if (fTPCinP1[3]>0.0000001){
      fTPCAngle1[0] = TMath::ATan2(fTPCinP1[1],fTPCinP1[0]);
      fTPCAngle1[1] = TMath::ATan(fTPCinP1[2]/fTPCinP1[3]);  
    }    
    Double_t cov[15], param[5],x;
    fESDTrack.GetInnerExternalCovariance(cov);
    fESDTrack.GetInnerExternalParameters(x,param);    
    //
    fTPCDelta[0] = (fTPCinR0[4]-fTPCinR1[4])*fTPCinR1[3];  //delta rfi
    fTPCPools[0] = fTPCDelta[0]/TMath::Sqrt(cov[0]);
    fTPCDelta[1] = (fTPCinR0[2]-fTPCinR1[2]);              //delta z
    fTPCPools[1] = fTPCDelta[1]/TMath::Sqrt(cov[2]);
    fTPCDelta[2] = (fTPCAngle0[0]-fTPCAngle1[0]);
    fTPCPools[2] = fTPCDelta[2]/TMath::Sqrt(cov[5]);
    fTPCDelta[3] = (TMath::Tan(fTPCAngle0[1])-TMath::Tan(fTPCAngle1[1]));
    fTPCPools[3] = fTPCDelta[3]/TMath::Sqrt(cov[9]);
    fTPCDelta[4] = (fTPCinP0[3]-fTPCinP1[3]);
    Double_t sign = (param[4]>0)? 1.:-1; 
    fTPCPools[4] = sign*(1./fTPCinP0[3]-1./fTPCinP1[3])/TMath::Sqrt(cov[14]);    
  }
  if (fITSOn){
    // ITS 
    Double_t param[5],x;
    //fESDTrack.GetExternalParameters(x,param);   
    fESDTrack.GetConstrainedExternalParameters(x,param);   
    Double_t cov[15];
    fESDTrack.GetExternalCovariance(cov);
    fESDTrack.GetConstrainedExternalCovariance(cov);
    if (TMath::Abs(param[4])<0.0000000001) return;

    fESDTrack.GetXYZ(fITSinR1);
    fESDTrack.GetPxPyPz(fITSinP1);
    fITSinP1[3] = TMath::Sqrt(fITSinP1[0]*fITSinP1[0]+fITSinP1[1]*fITSinP1[1]);
    //
    fITSinR1[3] = TMath::Sqrt(fITSinR1[0]*fITSinR1[0]+fITSinR1[1]*fITSinR1[1]);
    fITSinR1[4] = TMath::ATan2(fITSinR1[1],fITSinR1[0]);
    //
    //
    if (fITSinP1[3]>0.0000001){
      fITSAngle1[0] = TMath::ATan2(fITSinP1[1],fITSinP1[0]);
      fITSAngle1[1] = TMath::ATan(fITSinP1[2]/fITSinP1[3]);  
    }
    //
    //
    fITSDelta[0] = (fITSinR0[4]-fITSinR1[4])*fITSinR1[3];  //delta rfi
    fITSPools[0] = fITSDelta[0]/TMath::Sqrt(cov[0]);
    fITSDelta[1] = (fITSinR0[2]-fITSinR1[2]);              //delta z
    fITSPools[1] = fITSDelta[1]/TMath::Sqrt(cov[2]);
    fITSDelta[2] = (fITSAngle0[0]-fITSAngle1[0]);
    fITSPools[2] = fITSDelta[2]/TMath::Sqrt(cov[5]);
    fITSDelta[3] = (TMath::Tan(fITSAngle0[1])-TMath::Tan(fITSAngle1[1]));
    fITSPools[3] = fITSDelta[3]/TMath::Sqrt(cov[9]);
    fITSDelta[4] = (fITSinP0[3]-fITSinP1[3]);    
    Double_t sign = (param[4]>0) ? 1:-1; 
    fITSPools[4] = sign*(1./fITSinP0[3]-1./fITSinP1[3])/TMath::Sqrt(cov[14]);    
  }
  
}


void  AliESDRecV0Info::Update(Float_t vertex[3])
{
  Float_t v[3] = {fXr[0]-vertex[0],fXr[1]-vertex[1],fXr[2]-vertex[2]};
  Float_t p[3] = {fPdr[0]+fPm[0], fPdr[1]+fPm[1],fPdr[2]+fPm[2]};

  Float_t vnorm2 = v[0]*v[0]+v[1]*v[1];
  Float_t vnorm3 = TMath::Sqrt(v[2]*v[2]+vnorm2);
  vnorm2 = TMath::Sqrt(vnorm2);
  Float_t pnorm2 = p[0]*p[0]+p[1]*p[1];
  Float_t pnorm3 = TMath::Sqrt(p[2]*p[2]+pnorm2);
  pnorm2 = TMath::Sqrt(pnorm2);
  
  fPointAngleFi = (v[0]*p[0]+v[1]*p[1])/(vnorm2*pnorm2);
  fPointAngleTh = (v[2]*p[2]+vnorm2*pnorm2)/(vnorm3*pnorm3);  
  fPointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3);

}


////////////////////////////////////////////////////////////////////////
ESDCmpTr::ESDCmpTr()
{
  Reset();
}

////////////////////////////////////////////////////////////////////////
ESDCmpTr::ESDCmpTr(const char* fnGenTracks,
		   const char* fnCmp,
		   const char* fnGalice, Int_t direction,
		   Int_t nEvents, Int_t firstEvent)
{
  Reset();
  //  fFnGenTracks = fnGenTracks;
  //  fFnCmp = fnCmp;
  sprintf(fFnGenTracks,"%s",fnGenTracks);
  sprintf(fFnCmp,"%s",fnCmp);

  fFirstEventNr = firstEvent;
  fEventNr = firstEvent;
  fNEvents = nEvents;
  fDirection = direction;
  //
  fLoader = AliRunLoader::Open(fnGalice);
  if (gAlice){
    //delete gAlice->GetRunLoader();
    delete gAlice;
    gAlice = 0x0;
  }
  if (fLoader->LoadgAlice()){
    cerr<<"Error occured while l"<<endl;
  }
  Int_t nall = fLoader->GetNumberOfEvents();
  if (nEvents==0) {
    nEvents =nall;
    fNEvents=nall;
    fFirstEventNr=0;
  }    

  if (nall<=0){
    cerr<<"no events available"<<endl;
    fEventNr = 0;
    return;
  }
  if (firstEvent+nEvents>nall) {
    fEventNr = nall-firstEvent;
    cerr<<"restricted number of events availaible"<<endl;
  }
  AliMagF * magf = gAlice->Field();
  AliTracker::SetFieldMap(magf);

}


////////////////////////////////////////////////////////////////////////
ESDCmpTr::~ESDCmpTr()
{
  if (fLoader) {
    delete fLoader;
  }
}

//////////////////////////////////////////////////////////////
Int_t ESDCmpTr::SetIO()
{
  //
  // 
  CreateTreeCmp();
  if (!fTreeCmp) return 1;
  fParamTPC = GetTPCParam();
  //
  if (!ConnectGenTree()) {
    cerr<<"Cannot connect tree with generated tracks"<<endl;
    return 1;
  }
  return 0;
}

//////////////////////////////////////////////////////////////

Int_t ESDCmpTr::SetIO(Int_t eventNr)
{
  //
  // 
  // SET INPUT
  //
  TFile f("AliESDs.root");
  //
 
  TTree* tree = (TTree*) f.Get("esdTree");
  if (!tree) { 
    Char_t ename[100]; 
    sprintf(ename,"%d",eventNr);
    fEvent = (AliESD*)f.Get(ename);
    if (!fEvent){
      sprintf(ename,"ESD%d",eventNr);
      fEvent = (AliESD*)f.Get(ename);
    }
  }
  else{
    tree->SetBranchAddress("ESD", &fEvent);
    tree->GetEntry(eventNr);
  }


  /*
  Char_t ename[100]; 
  sprintf(ename,"%d",eventNr);
  fEvent = (AliESD*)f.Get(ename);
  if (!fEvent){
    sprintf(ename,"ESD%d",eventNr);
    fEvent = (AliESD*)f.Get(ename);
  }
  
  TTree* tree = (TTree*) f.Get("esdTree");
  if (!tree) {
    Error("CheckESD", "no ESD tree found");
    return kFALSE;
  }
  tree->SetBranchAddress("ESD", &fEvent);
  tree->GetEntry(eventNr);
  */

  if (!fEvent) return 1;

  return 0;
}



////////////////////////////////////////////////////////////////////////
void ESDCmpTr::Reset()
{
  fEventNr = 0;
  fNEvents = 0;
  fTreeCmp = 0;
  //  fFnCmp = "cmpTracks.root";
  fFileGenTracks = 0;
  fDebug = 0;
  //
  fParamTPC = 0;
  fEvent =0;
}

////////////////////////////////////////////////////////////////////////
Int_t ESDCmpTr::Exec(Int_t nEvents, Int_t firstEventNr)
{
  fNEvents = nEvents;
  fFirstEventNr = firstEventNr;
  return Exec();
}

////////////////////////////////////////////////////////////////////////
Int_t ESDCmpTr::Exec()
{
  TStopwatch timer;
  timer.Start();

  if (SetIO()==1) 
    return 1;
   
  fNextTreeGenEntryToRead = 0;
  cerr<<"fFirstEventNr, fNEvents: "<<fFirstEventNr<<" "<<fNEvents<<endl;
  for (Int_t eventNr = fFirstEventNr; eventNr < fFirstEventNr+fNEvents;
       eventNr++) {
    fEventNr = eventNr;
    SetIO(fEventNr);
    fNParticles = gAlice->GetEvent(fEventNr);    

    fIndexRecTracks = new Int_t[fNParticles*4];  //write at maximum 4 tracks corresponding to particle
    fFakeRecTracks = new Int_t[fNParticles];
    fMultiRecTracks = new Int_t[fNParticles];
    for (Int_t i = 0; i<fNParticles; i++) {
      fIndexRecTracks[4*i] = -1;
      fIndexRecTracks[4*i+1] = -1;
      fIndexRecTracks[4*i+2] = -1;
      fIndexRecTracks[4*i+3] = -1;

      fFakeRecTracks[i] = 0;
      fMultiRecTracks[i] = 0;
    }
  
    cout<<"Start to process event "<<fEventNr<<endl;
    cout<<"\tfNParticles = "<<fNParticles<<endl;
    if (fDebug>2) cout<<"\tStart loop over TreeT"<<endl;
    if (TreeTLoop()>0) return 1;

    if (fDebug>2) cout<<"\tStart loop over tree genTracks"<<endl;
    if (TreeGenLoop(eventNr)>0) return 1;
    if (fDebug>2) cout<<"\tEnd loop over tree genTracks"<<endl;

    delete [] fIndexRecTracks;
    delete [] fFakeRecTracks;
    delete [] fMultiRecTracks;
  }

  CloseOutputFile();

  cerr<<"Exec finished"<<endl;
  timer.Stop();
  timer.Print();
  return 0;

}
////////////////////////////////////////////////////////////////////////
Bool_t ESDCmpTr::ConnectGenTree()
{
//
// connect all branches from the genTracksTree
// use the same variables as for the new cmp tree, it may work
//
  fFileGenTracks = TFile::Open(fFnGenTracks,"READ");
  if (!fFileGenTracks) {
    cerr<<"Error in ConnectGenTree: cannot open file "<<fFnGenTracks<<endl;
    return kFALSE;
  }
  fTreeGenTracks = (TTree*)fFileGenTracks->Get("genTracksTree");
  if (!fTreeGenTracks) {
    cerr<<"Error in ConnectGenTree: cannot find genTracksTree in the file "
	<<fFnGenTracks<<endl;
    return kFALSE;
  }
  //
  fMCInfo = new AliMCInfo;
  fTreeGenTracks->SetBranchAddress("MC",&fMCInfo);
  //
  if (fDebug > 1) {
    cout<<"Number of gen. tracks with TR: "<<fTreeGenTracks->GetEntries()<<endl;
  }
  return kTRUE;
}


////////////////////////////////////////////////////////////////////////
void ESDCmpTr::CreateTreeCmp() 
{
  fFileCmp = TFile::Open(fFnCmp,"RECREATE");
  if (!fFileCmp) {
    cerr<<"Error in CreateTreeCmp: cannot open file "<<fFnCmp<<endl;
    return;
  }


  fTreeCmp    = new TTree("ESDcmpTracks","ESDcmpTracks");
  //
  fMCInfo = new AliMCInfo;
  fRecInfo = new AliESDRecInfo;
  //
  AliESDtrack * esdTrack = new AliESDtrack;
  AliITStrackV2 * itsTrack = new AliITStrackV2;  
  //
  fTreeCmp->Branch("MC","AliMCInfo",&fMCInfo);
  fTreeCmp->Branch("RC","AliESDRecInfo",&fRecInfo);
  fTreeCmp->Branch("fESDTrack","AliESDtrack",&esdTrack);
  //  fTreeCmp->Branch("ITS","AliITStrackV2",&itsTrack);
  delete esdTrack;
  //
  fTreeCmp->AutoSave(); 

}
////////////////////////////////////////////////////////////////////////
void ESDCmpTr::CloseOutputFile()  
{
  if (!fFileCmp) {
    cerr<<"File "<<fFnCmp<<" not found as an open file."<<endl;
    return;
  }
  fFileCmp->cd();
  fTreeCmp->Write();    
  delete fTreeCmp;
  
  fFileCmp->Close();
  delete fFileCmp;
  return;
}
////////////////////////////////////////////////////////////////////////

TVector3 ESDCmpTr::TR2Local(AliTrackReference *trackRef,
			    AliTPCParam *paramTPC) {

  Float_t x[3] = { trackRef->X(),trackRef->Y(),trackRef->Z()};
  Int_t index[4];
  paramTPC->Transform0to1(x,index);
  paramTPC->Transform1to2(x,index);
  return TVector3(x);
}
////////////////////////////////////////////////////////////////////////

Int_t ESDCmpTr::TreeTLoop()
{
  //
  // loop over all ESD reconstructed tracks and store info in memory
  //
  TStopwatch  timer;
  timer.Start();
  //  
  Int_t nEntries = (Int_t)fEvent->GetNumberOfTracks();  
  //
  fTracks      = new TObjArray(nEntries);
  //
  //load tracks to the memory
  for (Int_t i=0; i<nEntries;i++){
    AliESDtrack * track =fEvent->GetTrack(i);
    fTracks->AddAt(track,i);
  }
  //
  //
  AliESDtrack * track=0;
  for (Int_t iEntry=0; iEntry<nEntries;iEntry++){
    track = (AliESDtrack*)fTracks->UncheckedAt(iEntry);
    //
    Int_t label = track->GetLabel();
    Int_t absLabel = abs(label);
    if (absLabel < fNParticles) {
      //      fIndexRecTracks[absLabel] =  iEntry;
      if (label < 0) fFakeRecTracks[absLabel]++;      
      if (fMultiRecTracks[absLabel]>0){
	if (fMultiRecTracks[absLabel]<4)
	  fIndexRecTracks[absLabel*4+fMultiRecTracks[absLabel]] =  iEntry; 	
      }
      else      
	fIndexRecTracks[absLabel*4] =  iEntry;
      fMultiRecTracks[absLabel]++;
    }
  }  
  printf("Time spended in TreeTLoop\n");
  timer.Print();
  
  if (fDebug > 2) cerr<<"end of TreeTLoop"<<endl;  
  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t ESDCmpTr::TreeGenLoop(Int_t eventNr)
{
//
// loop over all entries for a given event, find corresponding 
// rec. track and store in the fTreeCmp
//
  TStopwatch timer;
  timer.Start();
  Int_t entry = fNextTreeGenEntryToRead;
  Double_t nParticlesTR = fTreeGenTracks->GetEntriesFast();
  cerr<<"fNParticles, nParticlesTR, fNextTreeGenEntryToRead: "<<fNParticles<<" "
      <<nParticlesTR<<" "<<fNextTreeGenEntryToRead<<endl;
  TBranch * branch = fTreeCmp->GetBranch("RC");
  branch->SetAddress(&fRecInfo); // set all pointers

  while (entry < nParticlesTR) {
    fTreeGenTracks->GetEntry(entry);
    entry++;
    if (eventNr < fMCInfo->fEventNr) continue;
    if (eventNr > fMCInfo->fEventNr) continue;;
    //
    fNextTreeGenEntryToRead = entry-1;
    if (fDebug > 2 && fMCInfo->fLabel < 10) {
      cerr<<"Fill track with a label "<<fMCInfo->fLabel<<endl;
    }
    //
    fRecInfo->Reset();
    AliESDtrack * track=0;
    fRecInfo->fReconstructed =0;
    TVector3 local = TR2Local(&(fMCInfo->fTrackRef),fParamTPC);
    local.GetXYZ(fRecInfo->fTRLocalCoord);	
    //
    if (fIndexRecTracks[fMCInfo->fLabel*4] >= 0) {
      track= (AliESDtrack*)fTracks->UncheckedAt(fIndexRecTracks[fMCInfo->fLabel*4]);
      //
      //
      // find nearest track if multifound
      if (fIndexRecTracks[fMCInfo->fLabel*4]+1){
	Double_t xyz[3];
	track->GetInnerXYZ(xyz);
	Float_t dz = TMath::Abs(local.Z()-xyz[2]);
	for (Int_t i=1;i<4;i++){
	  if (fIndexRecTracks[fMCInfo->fLabel*4+i]>=0){
	    AliESDtrack * track2 = (AliESDtrack*)fTracks->UncheckedAt(fIndexRecTracks[fMCInfo->fLabel*4+i]);
	    track2->GetInnerXYZ(xyz);
	    Float_t dz2 = TMath::Abs(local.Z()-xyz[2]);
	    if  (TMath::Abs(dz2)<dz)
	      track = track2;		   
	  }
	}
      }	
      //
      fRecInfo->fESDTrack =*track; 
      if (track->GetITStrack())
	fRecInfo->fITStrack = *((AliITStrackV2*)track->GetITStrack());
      if (track->GetTRDtrack()){
	fRecInfo->fTRDtrack = *((AliTRDtrack*)track->GetTRDtrack());
      }
      else{
	fRecInfo->fTRDtrack.SetdEdx(-1);
      }
      fRecInfo->fReconstructed = 1;
      fRecInfo->fFake     = fFakeRecTracks[fMCInfo->fLabel];
      fRecInfo->fMultiple = fMultiRecTracks[fMCInfo->fLabel];
      //
      fRecInfo->Update(fMCInfo,fParamTPC,kTRUE);          
    }
    else{
      fRecInfo->Update(fMCInfo,fParamTPC,kFALSE);
    }

    fTreeCmp->Fill();
  }
  fTreeCmp->AutoSave();
  fTracks->Delete();
  printf("Time spended in TreeGenLoop\n");
  timer.Print();
  if (fDebug > 2) cerr<<"end of TreeGenLoop"<<endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void AliESDComparisonDraw::SetIO(const char *fname)
{
  //
   TFile* file = TFile::Open(fname);
   //
   fTree = (TTree*) file->Get("ESDcmpTracks");
   if (!fTree) {
    printf("no track comparison tree found\n");
    file->Close();
    delete file;
  }
}


