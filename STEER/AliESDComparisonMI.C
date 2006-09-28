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
ESDCmpTr *t2 = new ESDCmpTr("genTracks.root","cmpESDTracks.root","galice.root",-1,0,0);
t2->Exec();

//
//some cuts definition
TCut cprim("cprim","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<0.01&&abs(MC.fVDist[2])<0.01")
//TCut cprim("cprim","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<0.5&&abs(MC.fVDist[2])<0.5")
//TCut citsin("citsin","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<3.9");
TCut citsin("citsin","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)<5");
TCut csec("csec","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)>0.5");


TCut crec("crec","fReconstructed==1");
TCut cteta1("cteta1","abs(MC.fParticle.Theta()/3.1415-0.5)<0.25");
TCut cteta05("cteta05","abs(MC.fParticle.Theta()/3.1415-0.5)<0.1");

TCut cpos1("cpos1","abs(MC.fParticle.fVz/sqrt(MC.fParticle.fVx*MC.fParticle.fVx+MC.fParticle.fVy*MC.fParticle.fVy))<1");
TCut csens("csens","abs(sqrt(fVDist[0]**2+fVDist[1]**2)-170)<50");
TCut cmuon("cmuon","abs(MC.fParticle.fPdgCode==-13)");
TCut cchi2("cchi2","fESDTrack.fITSchi2MIP[0]<7.&&fESDTrack.fITSchi2MIP[1]<5.&&fESDTrack.fITSchi2MIP[2]<7.&&fESDTrack.fITSchi2MIP[3]<7.5&&fESDTrack.fITSchi2MIP[4]<6.")

AliESDComparisonDraw comp;  
comp.SetIO(); 
TFile f("genHits.root");
TTree * treel = (TTree*)f.Get("HitLines");
if (treel) comp->fTree->AddFriend(treel,"L");

//
//example
comp.fTree->SetAlias("radius","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)");
comp.fTree->SetAlias("direction","MC.fParticle.fVx*MC.fParticle.fPx+MC.fParticle.fVy*MC.fParticle.fPy");
comp.fTree->SetAlias("decaydir","MC.fTRdecay.fX*MC.fTRdecay.fPx+MC.fTRdecay.fY*MC.fTRdecay.fPy");
comp.fTree->SetAlias("theta","MC.fTrackRef.Theta()");
comp.fTree->SetAlias("primdca","sqrt(RC.fITStrack.fD[0]**2+RC.fITStrack.fD[1]**2)");
comp.fTree->SetAlias("trdchi2","fTRDtrack.fChi2/fTRDtrack.fN");
comp.fTree->SetAlias("trdchi2","fTRDtrack.fChi2/fTRDtrack.fN");


TH1F his("his","his",100,0,20);
TH1F hpools("hpools","hpools",100,-7,7);
TH1F hfake("hfake","hfake",1000,0,150);
TProfile profp0("profp0","profp0",20,-0.4,0.9)

comp.DrawXY("fTPCinP0[3]","fTPCDelta[4]/fTPCinP1[3]","fReconstructed==1"+cprim,"1",4,0.2,1.5,-0.06,0.06)
comp.fRes->Draw();
comp.fMean->Draw();  

comp.DrawXY("fITSinP0[3]","fITSDelta[4]/fITSinP1[3]","fReconstructed==1&&fITSOn"+cprim,"1",4,0.2,1.5,-0.06,0.06)
comp.fRes->Draw();

comp.Eff("fTPCinP0[3]","fRowsWithDigits>120"+cteta1+cpos1+cprim,"fTPCOn",20,0.2,1.5)
comp.fRes->Draw();

comp.Eff("fTPCinP0[3]","fRowsWithDigits>120"+cteta1+cpos1+cprim,"fTPCOn&&fITSOn&&fESDTrack.fITSFakeRatio<0.1",10,0.2,1.5)
comp.fRes->Draw();
comp.Eff("fTPCinP0[3]","fRowsWithDigits>120"+cteta1+cpos1+cprim,"fTPCOn&&fITSOn&&fESDTrack.fITSFakeRatio>0.1",10,0.2,1.5)
comp.fRes->Draw();

comp.fTree->Draw("fESDTrack.fITSsignal/fESDTrack.fTPCsignal","fITSOn&&fTPCOn&&fESDTrack.fITSFakeRatio==0") 

TH1F his("his","his",100,0,20);
TH1F hpools("hpools","hpools",100,-7,7);

TH2F * hdedx0 = new TH2F("dEdx0","dEdx0",100, 0,2,200,0,550); hdedx0->SetMarkerColor(1);
TH2F * hdedx1 = new TH2F("dEdx1","dEdx1",100, 0,2,200,0,550); hdedx1->SetMarkerColor(4);
TH2F * hdedx2 = new TH2F("dEdx2","dEdx2",100, 0,2,200,0,550); hdedx2->SetMarkerColor(3);
TH2F * hdedx3 = new TH2F("dEdx3","dEdx3",100, 0,2,200,0,550); hdedx3->SetMarkerColor(2);

comp.fTree->Draw("fESDTrack.fITSsignal:MC.fParticle.P()>>dEdx0","fITSOn&&abs(fPdg)==211&&fITStrack.fN==6"+cprim) 
comp.fTree->Draw("fESDTrack.fITSsignal:MC.fParticle.P()>>dEdx1","fITSOn&&abs(fPdg)==2212&&fITStrack.fN==6"+cprim) 
comp.fTree->Draw("fESDTrack.fITSsignal:MC.fParticle.P()>>dEdx2","fITSOn&&abs(fPdg)==321&&fITStrack.fN==6"+cprim) 
comp.fTree->Draw("fESDTrack.fITSsignal:MC.fParticle.P()>>dEdx3","fITSOn&&abs(fPdg)==11&&fITStrack.fN==6"+cprim) 


comp.fTree->Draw("fESDTrack.fTRDsignal:MC.fParticle.P()>>dEdx0","fTRDOn&&abs(fPdg)==211&&fTRDtrack.fN>40&&fStatus[2]>1") 
comp.fTree->Draw("fESDTrack.fTRDsignal:MC.fParticle.P()>>dEdx1","fTRDOn&&abs(fPdg)==2212&&fTRDtrack.fN>40&&fStatus[2]>1") 
comp.fTree->Draw("fESDTrack.fTRDsignal:MC.fParticle.P()>>dEdx2","fTRDOn&&abs(fPdg)==321&&fTRDtrack.fN>40&&fStatus[2]>1") 
comp.fTree->Draw("fESDTrack.fTRDsignal:MC.fParticle.P()>>dEdx3","fTRDOn&&abs(fPdg)==11&&fTRDtrack.fN>40&&fStatus[2]>1") 

comp.fTree->Draw("fESDTrack.fTPCsignal:fTPCinP0[4]>>dEdx0","fTPCOn&&abs(fPdg)==211&&fESDTrack.fTPCncls>180&&fESDTrack.fTPCsignal>10"+cteta1); 
comp.fTree->Draw("fESDTrack.fTPCsignal:fTPCinP0[4]>>dEdx1","fTPCOn&&abs(fPdg)==2212&&fESDTrack.fTPCncls>180&&fESDTrack.fTPCsignal>10"+cteta1); 
comp.fTree->Draw("fESDTrack.fTPCsignal:fTPCinP0[4]>>dEdx2","fTPCOn&&abs(fPdg)==321&&fESDTrack.fTPCncls>180&&fESDTrack.fTPCsignal>10"+cteta1); 
comp.fTree->Draw("fESDTrack.fTPCsignal:fTPCinP0[4]>>dEdx3","fTPCOn&&abs(fPdg)==11&&fESDTrack.fTPCncls>180&&fESDTrack.fTPCsignal>10"+cteta1); 

hdedx3->SetXTitle("P(GeV/c)");
hdedx3->SetYTitle("dEdx(unit)");
hdedx3->Draw(); hdedx1->Draw("same"); hdedx2->Draw("same"); hdedx0->Draw("same");

comp.DrawXY("fITSinP0[3]","fITSPools[4]","fReconstructed==1&&fPdg==-211&&fITSOn"+cprim,"1",4,0.2,1.0,-8,8)

TProfile prof("prof","prof",10,0.5,5);




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
#include "AliITStrackMI.h"
#include "AliTRDtrack.h"
#include "AliHelix.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliESDkink.h"
#include "AliESDv0.h"

#endif
#include "AliGenInfo.h"
#include "AliESDComparisonMI.h"





void MakeAliases(AliESDComparisonDraw&comp)
{
  //
  // aliases definition
  //
  comp.fTree->SetAlias("radius","TMath::Sqrt(MC.fVDist[0]**2+MC.fVDist[1]**2)");
  comp.fTree->SetAlias("direction","MC.fParticle.fVx*MC.fParticle.fPx+MC.fParticle.fVy*MC.fParticle.fPy");
  comp.fTree->SetAlias("decaydir","MC.fTRdecay.fX*MC.fTRdecay.fPx+MC.fTRdecay.fY*MC.fTRdecay.fPy");
  comp.fTree->SetAlias("theta","MC.fTrackRef.Theta()");
  comp.fTree->SetAlias("primdca","sqrt(RC.fITStrack.fD[0]**2+RC.fITStrack.fD[1]**2)");
  comp.fTree->SetAlias("trdchi2","fTRDtrack.fChi2/fTRDtrack.fN");
  comp.fTree->SetAlias("trdchi2","fTRDtrack.fChi2/fTRDtrack.fN");
  
  comp.fTree->SetAlias("trddedx","(RC.fESDTrack.fTRDsignals[0]+RC.fESDTrack.fTRDsignals[1]+RC.fESDTrack.fTRDsignals[2]+RC.fESDTrack.fTRDsignals[3]+RC.fESDTrack.fTRDsignals[4]+RC.fESDTrack.fTRDsignals[5])/6.");
  
  comp.fTree->SetAlias("dtofmc2","fESDTrack.fTrackTime[2]-(10^12*MC.fTOFReferences[0].fTime)");
  comp.fTree->SetAlias("dtofrc2","(fESDTrack.fTrackTime[2]-fESDTrack.fTOFsignal)");

  comp.fTree->SetAlias("psum","fESDTrack.fTOFr[4]+fESDTrack.fTOFr[3]+fESDTrack.fTOFr[2]+fESDTrack.fTOFr[1]+fESDTrack.fTOFr[0]");
  comp.fTree->SetAlias("P0","fESDTrack.fTOFr[0]/psum");
  comp.fTree->SetAlias("P1","fESDTrack.fTOFr[1]/psum");
  comp.fTree->SetAlias("P2","fESDTrack.fTOFr[2]/psum");
  comp.fTree->SetAlias("P3","fESDTrack.fTOFr[3]/psum");
  comp.fTree->SetAlias("P4","fESDTrack.fTOFr[4]/psum");
  comp.fTree->SetAlias("MaxP","max(max(max(P0,P1),max(P2,P3)),P4)");
}




void  AliESDRecInfo::UpdatePoints(AliESDtrack*track)
{
  //
  //
  Int_t iclusters[200];
  Float_t density[160];
  for (Int_t i=0;i<160;i++) density[i]=-1.;
  fTPCPoints[0]= 160;
  fTPCPoints[1] = -1;
  //
  if (fTPCPoints[0]<fTPCPoints[1]) return;
  //  Int_t nclusters=track->GetTPCclusters(iclusters);

  Int_t ngood=0;
  Int_t undeff=0;
  Int_t nall =0;
  Int_t range=20;
  for (Int_t i=0;i<160;i++){
    Int_t last = i-range;
    if (nall<range) nall++;
    if (last>=0){
      if (iclusters[last]>0&& (iclusters[last]&0x8000)==0) ngood--;
      if (iclusters[last]==-1) undeff--;
    }
    if (iclusters[i]>0&& (iclusters[i]&0x8000)==0)   ngood++;
    if (iclusters[i]==-1) undeff++;
    if (nall==range &&undeff<range/2) density[i-range/2] = Float_t(ngood)/Float_t(nall-undeff);
  }
  Float_t maxdens=0;
  Int_t indexmax =0;
  for (Int_t i=0;i<160;i++){
    if (density[i]<0) continue;
    if (density[i]>maxdens){
      maxdens=density[i];
      indexmax=i;
    }
  }
  //
  //max dens point
  fTPCPoints[3] = maxdens;
  fTPCPoints[1] = indexmax;
  //
  // last point
  for (Int_t i=indexmax;i<160;i++){
    if (density[i]<0) continue;
    if (density[i]<maxdens/2.) {
      break;
    }
    fTPCPoints[2]=i;
  }
  //
  // first point
  for (Int_t i=indexmax;i>0;i--){
    if (density[i]<0) continue;
    if (density[i]<maxdens/2.) {
      break;
    }
    fTPCPoints[0]=i;
  }
  //
  // Density at the last 30 padrows
  //
  // 
  nall  = 0;
  ngood = 0;
  for (Int_t i=159;i>0;i--){
    if (iclusters[i]==-1) continue; //dead zone
    nall++;
    if (iclusters[i]>0)   ngood++;
    if (nall>20) break;
  }
  fTPCPoints[4] = Float_t(ngood)/Float_t(nall);
  //
  if ((track->GetStatus()&AliESDtrack::kITSrefit)>0) fTPCPoints[0]=-1;
  //
  //
  // check TRDPoints
  /*
  nclusters=track->GetTRDclusters(iclusters);
  for (Int_t i=nclusters;i>0;i--){
    
  }
  */


}

//
//
void AliESDRecInfo::Update(AliMCInfo* info,AliTPCParam * /*par*/, Bool_t reconstructed, AliESD *event)
{
  //
  //
  //calculates derived variables
  //  
  //
  UpdatePoints(&fESDTrack);
  fBestTOFmatch=1000;
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
  for (Int_t i=0;i<4;i++) fStatus[i] =0;
  fReconstructed = kFALSE;
  fTPCOn = kFALSE;
  fITSOn = kFALSE;
  fTRDOn = kFALSE;  
  if (reconstructed==kFALSE) return;

  fLabels[0] = info->fLabel;
  fLabels[1] = info->fPrimPart;
  fReconstructed = kTRUE;
  fTPCOn = ((fESDTrack.GetStatus()&AliESDtrack::kTPCrefit)>0) ? kTRUE : kFALSE;
  fITSOn = ((fESDTrack.GetStatus()&AliESDtrack::kITSrefit)>0) ? kTRUE : kFALSE;
  fTRDOn = ((fESDTrack.GetStatus()&AliESDtrack::kTRDrefit)>0) ? kTRUE : kFALSE;
  //
  //  
  if ((fESDTrack.GetStatus()&AliESDtrack::kTPCrefit)>0){
    fStatus[1] =3;
  }
  else{
    if ((fESDTrack.GetStatus()&AliESDtrack::kTPCout)>0){
      fStatus[1] =2;
    }
    else{
      if ((fESDTrack.GetStatus()&AliESDtrack::kTPCin)>0)
	fStatus[1]=1;
    }      
  }
  //
  if ((fESDTrack.GetStatus()&AliESDtrack::kITSout)>0){
    fStatus[0] =2;
  }
  else{
    if ((fESDTrack.GetStatus()&AliESDtrack::kITSrefit)>0){
      fStatus[0] =1;
    }
    else{
      fStatus[0]=0;
    }      
  }

  //
  //
  if ((fESDTrack.GetStatus()&AliESDtrack::kTRDrefit)>0){
    fStatus[2] =2;
  }
  else{
    if ((fESDTrack.GetStatus()&AliESDtrack::kTRDout)>0){
      fStatus[2] =1;
    }
  }
  if ((fESDTrack.GetStatus()&AliESDtrack::kTRDStop)>0){
    fStatus[2] =10;
  }

  //
  //TOF 
  // 
  if (((fESDTrack.GetStatus()&AliESDtrack::kTOFout)>0)){
    //
    // best tof match
    Double_t times[5];
    fESDTrack.GetIntegratedTimes(times);    
    for (Int_t i=0;i<5;i++){
      if ( TMath::Abs(fESDTrack.GetTOFsignal()-times[i]) <TMath::Abs(fBestTOFmatch) ){
	fBestTOFmatch = fESDTrack.GetTOFsignal()-times[i];
      }
    }
    Int_t toflabel[3];
    fESDTrack.GetTOFLabel(toflabel);
    Bool_t toffake=kTRUE;
    Bool_t tofdaughter=kFALSE;
    for (Int_t i=0;i<3;i++){
      if (toflabel[i]<0) continue;      
      if (toflabel[i]== TMath::Abs(fESDTrack.GetLabel()))  toffake=kFALSE;	
      if (toflabel[i]==info->fParticle.GetDaughter(0) || (toflabel[i]==info->fParticle.GetDaughter(1))) tofdaughter=kTRUE;  // decay product of original particle
      fStatus[3]=1;
    }
    if (toffake) fStatus[3] =3;       //total fake
    if (tofdaughter) fStatus[3]=2;    //fake because of decay
  }else{
    fStatus[3]=0;
  }


  if (fStatus[1]>0 &&info->fNTPCRef>0&&TMath::Abs(fTPCinP0[3])>0.0001){
    //TPC
    fESDTrack.GetInnerXYZ(fTPCinR1);
    fTPCinR1[3] = TMath::Sqrt(fTPCinR1[0]*fTPCinR1[0]+fTPCinR1[1]*fTPCinR1[1]);
    fTPCinR1[4] = TMath::ATan2(fTPCinR1[1],fTPCinR1[0]);	
    fESDTrack.GetInnerPxPyPz(fTPCinP1);
    fTPCinP1[3] = TMath::Sqrt(fTPCinP1[0]*fTPCinP1[0]+fTPCinP1[1]*fTPCinP1[1]);
    fTPCinP1[4] = TMath::Sqrt(fTPCinP1[3]*fTPCinP1[3]+fTPCinP1[2]*fTPCinP1[2]);
    //
    //
    if (fTPCinP1[3]>0.000000000000001){
      fTPCAngle1[0] = TMath::ATan2(fTPCinP1[1],fTPCinP1[0]);
      fTPCAngle1[1] = TMath::ATan(fTPCinP1[2]/fTPCinP1[3]);  
    }    
    Double_t cov[15], param[5],x, alpha;
    fESDTrack.GetInnerExternalCovariance(cov);
    fESDTrack.GetInnerExternalParameters(alpha, x,param);
    if (x<50) return ;
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
    fSign =sign;
    fTPCPools[4] = sign*(1./fTPCinP0[3]-1./fTPCinP1[3])/TMath::Sqrt(TMath::Abs(cov[14]));
  }
  if (fITSOn){
    // ITS 
    Double_t param[5],x;
    fESDTrack.GetExternalParameters(x,param);   
    //    fESDTrack.GetConstrainedExternalParameters(x,param);   
    Double_t cov[15];
    fESDTrack.GetExternalCovariance(cov);
    //fESDTrack.GetConstrainedExternalCovariance(cov);
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
    fSign = sign;
    fITSPools[4] = sign*(1./fITSinP0[3]-1./fITSinP1[3])/TMath::Sqrt(cov[14]);    
  }
  
}


void  AliESDRecV0Info::Update(Float_t vertex[3])
{ 

  if ( (fT1.fStatus[1]>0)&& (fT2.fStatus[1]>0)){
    Float_t distance1,distance2;
    Double_t xx[3],pp[3];
    //
    Double_t xd[3],pd[3],signd;
    Double_t xm[3],pm[3],signm;
    //
    //
    if (fT1.fITSOn&&fT2.fITSOn){
      for (Int_t i=0;i<3;i++){
	xd[i] = fT2.fITSinR1[i];
	pd[i] = fT2.fITSinP1[i];
	xm[i] = fT1.fITSinR1[i];
	pm[i] = fT1.fITSinP1[i];
      }
    }
    else{
      
      for (Int_t i=0;i<3;i++){
	xd[i] = fT2.fTPCinR1[i];
	pd[i] = fT2.fTPCinP1[i];
	xm[i] = fT1.fTPCinR1[i];
	pm[i] = fT1.fTPCinP1[i];
      }
    }
    //
    //
    signd =  fT2.fSign<0 ? -1:1;
    signm =  fT1.fSign<0 ? -1:1;

    AliHelix dhelix1(xd,pd,signd);
    dhelix1.GetMomentum(0,pp,0);
    dhelix1.Evaluate(0,xx);      
    // 
    //  Double_t x2[3],p2[3];
    //            
    AliHelix mhelix(xm,pm,signm);    
    //
    //find intersection linear
    //
    Double_t phase[2][2],radius[2];
    Int_t  points = dhelix1.GetRPHIintersections(mhelix, phase, radius,200);
    Double_t delta1=10000,delta2=10000;  

    if (points==1){
      fRs[0] = TMath::Sqrt(radius[0]);
      fRs[1] = TMath::Sqrt(radius[0]);
    }
    if (points==2){
      fRs[0] =TMath::Min(TMath::Sqrt(radius[0]),TMath::Sqrt(radius[1]));
      fRs[1] =TMath::Max(TMath::Sqrt(radius[0]),TMath::Sqrt(radius[1]));
    }
    
    if (points>0){
      dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
      dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
      dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    }
    if (points==2){    
      dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
      dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
      dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    }
    if (points==1){
      fRs[0] = TMath::Sqrt(radius[0]);
      fRs[1] = TMath::Sqrt(radius[0]);
      fDistMinR = delta1;
    }
    if (points==2){
      if (radius[0]<radius[1]){
	fRs[0] = TMath::Sqrt(radius[0]);
	fRs[1] = TMath::Sqrt(radius[1]);
	fDistMinR = delta1;
      }
      else{
	fRs[0] = TMath::Sqrt(radius[1]);
	fRs[1] = TMath::Sqrt(radius[0]);
	fDistMinR = delta2;
      }
    }
    //
    //
    distance1 = TMath::Min(delta1,delta2);
    //
    //find intersection parabolic
    //
    points = dhelix1.GetRPHIintersections(mhelix, phase, radius);
    delta1=10000,delta2=10000;  
    
    if (points>0){
      dhelix1.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    }
    if (points==2){    
      dhelix1.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    }
    
    distance2 = TMath::Min(delta1,delta2);
    if (distance2>100) fDist2 =100;
    return;
    if (delta1<delta2){
      //get V0 info
      dhelix1.Evaluate(phase[0][0],fXr);
      dhelix1.GetMomentum(phase[0][0],fPdr);
      mhelix.GetMomentum(phase[0][1],fPm);
      dhelix1.GetAngle(phase[0][0],mhelix,phase[0][1],fAngle);
      fRr = TMath::Sqrt(radius[0]);
    }
    else{
      dhelix1.Evaluate(phase[1][0],fXr);
      dhelix1.GetMomentum(phase[1][0], fPdr);
      mhelix.GetMomentum(phase[1][1], fPm);
      dhelix1.GetAngle(phase[1][0],mhelix,phase[1][1],fAngle);
      fRr = TMath::Sqrt(radius[1]);
    }
    fDist1 = TMath::Sqrt(distance1);
    fDist2 = TMath::Sqrt(distance2);      
    
    if (fDist2<10.5){
      Double_t x,alpha,param[5],cov[15];
      //
      fT1.fESDTrack.GetInnerExternalParameters(alpha,x,param);
      fT1.fESDTrack.GetInnerExternalCovariance(cov);
      AliExternalTrackParam paramm(x,alpha,param,cov);
      //
      fT2.fESDTrack.GetInnerExternalParameters(alpha,x,param);
      fT2.fESDTrack.GetInnerExternalCovariance(cov);
      AliExternalTrackParam paramd(x,alpha,param,cov);
    }    
    //            
    //   
    
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
}

////
void  AliESDRecKinkInfo::Update()
{

  if ( (fT1.fTPCOn)&& (fT2.fTPCOn)){
    //
    // IF BOTH RECONSTRUCTED
    Float_t distance1,distance2;
    Double_t xx[3],pp[3];
    //
    Double_t xd[3],pd[3],signd;
    Double_t xm[3],pm[3],signm;
    for (Int_t i=0;i<3;i++){
      xd[i] = fT2.fTPCinR1[i];
      pd[i] = fT2.fTPCinP1[i];
      xm[i] = fT1.fTPCinR1[i];
      pm[i] = fT1.fTPCinP1[i];
    }
    signd =  fT2.fSign<0 ? -1:1;
    signm =  fT1.fSign<0 ? -1:1;

    AliHelix dhelix1(xd,pd,signd);
    dhelix1.GetMomentum(0,pp,0);
    dhelix1.Evaluate(0,xx);      
    // 
    //  Double_t x2[3],p2[3];
    //            
    AliHelix mhelix(xm,pm,signm);    
    //
    //find intersection linear
    //
    Double_t phase[2][2],radius[2];
    Int_t  points = dhelix1.GetRPHIintersections(mhelix, phase, radius,200);
    Double_t delta1=10000,delta2=10000;  

    if (points==1){
      fMinR = TMath::Sqrt(radius[0]);
    }
    if (points==2){
      fMinR =TMath::Min(TMath::Sqrt(radius[0]),TMath::Sqrt(radius[1]));
    }
    
    if (points>0){
      dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
      dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
      dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    }
    if (points==2){    
      dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
      dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
      dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    }
    if (points==1){
      fMinR = TMath::Sqrt(radius[0]);
      fDistMinR = delta1;
    }
    if (points==2){
      if (radius[0]<radius[1]){
	fMinR = TMath::Sqrt(radius[0]);
	fDistMinR = delta1;
      }
      else{
	fMinR = TMath::Sqrt(radius[1]);
	fDistMinR = delta2;
      }
    }
    //
    //
    distance1 = TMath::Min(delta1,delta2);
    //
    //find intersection parabolic
    //
    points = dhelix1.GetRPHIintersections(mhelix, phase, radius);
    delta1=10000,delta2=10000;  
    
    if (points>0){
      dhelix1.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    }
    if (points==2){    
      dhelix1.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    }
    
    distance2 = TMath::Min(delta1,delta2);
    if (delta1<delta2){
      //get V0 info
      dhelix1.Evaluate(phase[0][0],fXr);
      dhelix1.GetMomentum(phase[0][0],fPdr);
      mhelix.GetMomentum(phase[0][1],fPm);
      dhelix1.GetAngle(phase[0][0],mhelix,phase[0][1],fAngle);
      fRr = TMath::Sqrt(radius[0]);
    }
    else{
      dhelix1.Evaluate(phase[1][0],fXr);
      dhelix1.GetMomentum(phase[1][0], fPdr);
      mhelix.GetMomentum(phase[1][1], fPm);
      dhelix1.GetAngle(phase[1][0],mhelix,phase[1][1],fAngle);
      fRr = TMath::Sqrt(radius[1]);
    }
    fDist1 = TMath::Sqrt(distance1);
    fDist2 = TMath::Sqrt(distance2);      
    
    if (fDist2<10.5){
      Double_t x,alpha,param[5],cov[15];
      //
      fT1.fESDTrack.GetInnerExternalParameters(alpha,x,param);
      fT1.fESDTrack.GetInnerExternalCovariance(cov);
      AliExternalTrackParam paramm(x,alpha,param,cov);
      //
      fT2.fESDTrack.GetInnerExternalParameters(alpha,x,param);
      fT2.fESDTrack.GetInnerExternalCovariance(cov);
      AliExternalTrackParam paramd(x,alpha,param,cov);
      /*
      AliESDkink kink;
      kink.Update(&paramm,&paramd);
      //      kink.Dump();
      Double_t diff  = kink.fRr-fRr;
      Double_t diff2 = kink.fDist2-fDist2;
      printf("Diff\t%f\t%f\n",diff,diff2);
      */
    }
    
    //            
    //
  }

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
  AliTracker::SetFieldMap(magf,0);

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
  fTreeCmpKinks =0;
  fTreeCmpV0 =0;
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
  fNextKinkToRead = 0;
  fNextV0ToRead   =0;
  cerr<<"fFirstEventNr, fNEvents: "<<fFirstEventNr<<" "<<fNEvents<<endl;
  for (Int_t eventNr = fFirstEventNr; eventNr < fFirstEventNr+fNEvents;
       eventNr++) {
    fEventNr = eventNr;
    SetIO(fEventNr);
    fNParticles = gAlice->GetEvent(fEventNr);    

    fIndexRecTracks = new Short_t[fNParticles*20];  //write at maximum 4 tracks corresponding to particle
    fIndexRecKinks  = new Short_t[fNParticles*20];  //write at maximum 20 tracks corresponding to particle
    fIndexRecV0  = new Short_t[fNParticles*20];  //write at maximum 20 tracks corresponding to particle

    fFakeRecTracks = new Short_t[fNParticles];
    fMultiRecTracks = new Short_t[fNParticles];
    fMultiRecKinks = new Short_t[fNParticles];
    fMultiRecV0 = new Short_t[fNParticles];

    for (Int_t i = 0; i<fNParticles; i++) {
      for (Int_t j=0;j<20;j++){
	fIndexRecTracks[20*i+j] = -1;
	fIndexRecKinks[20*i+j]  = -1;
	fIndexRecV0[20*i+j]  = -1;
      }
      fFakeRecTracks[i] = 0;
      fMultiRecTracks[i] = 0;
      fMultiRecKinks[i] = 0;
      fMultiRecV0[i] = 0;      
    }
  
    cout<<"Start to process event "<<fEventNr<<endl;
    cout<<"\tfNParticles = "<<fNParticles<<endl;
    if (fDebug>2) cout<<"\tStart loop over TreeT"<<endl;
    if (TreeTLoop()>0) return 1;

    if (fDebug>2) cout<<"\tStart loop over tree genTracks"<<endl;
    if (TreeGenLoop(eventNr)>0) return 1;
    BuildKinkInfo0(eventNr);
    BuildV0Info(eventNr);
    fRecArray->Delete();

    if (fDebug>2) cout<<"\tEnd loop over tree genTracks"<<endl;

    delete [] fIndexRecTracks;
    delete [] fIndexRecKinks;
    delete [] fIndexRecV0;
    delete [] fFakeRecTracks;
    delete [] fMultiRecTracks;
    delete [] fMultiRecKinks;
    delete [] fMultiRecV0;
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
  //
  fTreeGenKinks = (TTree*)fFileGenTracks->Get("genKinksTree");
  if (!fTreeGenKinks) {
    cerr<<"Error in ConnectGenTree: cannot find genTracksTree in the file "
	<<fFnGenTracks<<endl;
    //return kFALSE;
  }
  else{
    fGenKinkInfo = new AliGenKinkInfo;
    fTreeGenKinks->SetBranchAddress("MC",&fGenKinkInfo);
  }

  fTreeGenV0 = (TTree*)fFileGenTracks->Get("genV0Tree");
  if (!fTreeGenV0) {
    cerr<<"Error in ConnectGenTree: cannot find genTracksTree in the file "
	<<fFnGenTracks<<endl;
    //return kFALSE;
  }
  else{
    fGenV0Info = new AliGenV0Info;
    fTreeGenV0->SetBranchAddress("MC",&fGenV0Info);
  }
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
  //
  //
  fTreeCmp    = new TTree("ESDcmpTracks","ESDcmpTracks");
  fMCInfo = new AliMCInfo;
  fRecInfo = new AliESDRecInfo;
  AliESDtrack * esdTrack = new AliESDtrack;
  //  AliITStrackMI * itsTrack = new AliITStrackMI;  
  fTreeCmp->Branch("MC","AliMCInfo",&fMCInfo,256000);
  fTreeCmp->Branch("RC","AliESDRecInfo",&fRecInfo,256000);
  //  fTreeCmp->Branch("fESDTrack","AliESDtrack",&esdTrack);
  //  fTreeCmp->Branch("ITS","AliITStrackMI",&itsTrack);
  delete esdTrack;
  //
  //
  fTreeCmpKinks    = new TTree("ESDcmpKinks","ESDcmpKinks"); 
  fGenKinkInfo     = new AliGenKinkInfo;
  fRecKinkInfo     = new AliESDRecKinkInfo;
  fTreeCmpKinks->Branch("MC.","AliGenKinkInfo",&fGenKinkInfo,256000);
  fTreeCmpKinks->Branch("RC.","AliESDRecKinkInfo",&fRecKinkInfo,256000);
  //
  //
  fTreeCmpV0       = new TTree("ESDcmpV0","ESDcmpV0"); 
  fGenV0Info     = new AliGenV0Info;
  fRecV0Info     = new AliESDRecV0Info;
  fTreeCmpV0->Branch("MC.","AliGenV0Info",   &fGenV0Info,256000);
  fTreeCmpV0->Branch("RC.","AliESDRecV0Info",&fRecV0Info,256000);
  //
  fTreeCmp->AutoSave(); 
  fTreeCmpKinks->AutoSave(); 
  fTreeCmpV0->AutoSave(); 
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
  paramTPC->Transform1to2Ideal(x,index);
  return TVector3(x);
}
////////////////////////////////////////////////////////////////////////

Int_t ESDCmpTr::TreeTLoop()
{
  //
  // loop over all ESD reconstructed tracks and store info in memory
  //
  // + loop over all reconstructed kinks
  TStopwatch  timer;
  timer.Start();
  //  
  Int_t nEntries = (Int_t)fEvent->GetNumberOfTracks();  
  Int_t nKinks = (Int_t) fEvent->GetNumberOfKinks();
  Int_t nV0MIs = (Int_t) fEvent->GetNumberOfV0s();
  fSignedKinks = new Short_t[nKinks];
  fSignedV0    = new Short_t[nV0MIs];
  //
  // load kinks to the memory
  for (Int_t i=0; i<nKinks;i++){
    AliESDkink * kink =fEvent->GetKink(i);
    fSignedKinks[i]=0;
    if (kink->fStatus<0) continue;
  }
  //
  for (Int_t i=0; i<nV0MIs;i++){
    AliESDv0 * v0MI =fEvent->GetV0(i);
    fSignedV0[i]=0;
    if (v0MI->fStatus<0) continue;
  }
  
  //
  //
  AliESDtrack * track=0;
  for (Int_t iEntry=0; iEntry<nEntries;iEntry++){
    //track = (AliESDtrack*)fTracks->UncheckedAt(iEntry);
    track = (AliESDtrack*)fEvent->GetTrack(iEntry);
    //
    Int_t label = track->GetLabel();
    Int_t absLabel = abs(label);
    if (absLabel < fNParticles) {
      //      fIndexRecTracks[absLabel] =  iEntry;
      if (label < 0) fFakeRecTracks[absLabel]++;      
      if (fMultiRecTracks[absLabel]>0){
	if (fMultiRecTracks[absLabel]<20)
	  fIndexRecTracks[absLabel*20+fMultiRecTracks[absLabel]] =  iEntry; 	
      }
      else      
	fIndexRecTracks[absLabel*20] =  iEntry;
      fMultiRecTracks[absLabel]++;
    }
  }
  // sort reconstructed kinks  
  //
  AliESDkink * kink=0;
  for (Int_t iEntry=0; iEntry<nKinks;iEntry++){
    kink = (AliESDkink*)fEvent->GetKink(iEntry);
    if (!kink) continue;
    //
    Int_t label0 = TMath::Abs(kink->fLab[0]);
    Int_t label1 = TMath::Abs(kink->fLab[1]);
    Int_t absLabel = TMath::Min(label0,label1);
    if (absLabel < fNParticles) {
      if (fMultiRecKinks[absLabel]>0){
	if (fMultiRecKinks[absLabel]<20)
	  fIndexRecKinks[absLabel*20+fMultiRecKinks[absLabel]] =  iEntry; 	
      }
      else      
	fIndexRecKinks[absLabel*20] =  iEntry;
      fMultiRecKinks[absLabel]++;
    }
  }  
  // --sort reconstructed V0
  //
  AliESDv0 * v0MI=0;
  for (Int_t iEntry=0; iEntry<nV0MIs;iEntry++){
    v0MI = fEvent->GetV0(iEntry);
    if (!v0MI) continue;
    //
    //    Int_t label0 = TMath::Abs(v0MI->fLab[0]);
    //Int_t label1 = TMath::Abs(v0MI->fLab[1]);
    //
    for (Int_t i=0;i<2;i++){
      // Int_t absLabel = TMath::Min(label0,label1);
      Int_t absLabel =  TMath::Abs(v0MI->fLab[i]);
      if (absLabel < fNParticles) {
	if (fMultiRecV0[absLabel]>0){
	  if (fMultiRecV0[absLabel]<20)
	    fIndexRecV0[absLabel*20+fMultiRecV0[absLabel]] =  iEntry; 	
	}
	else      
	  fIndexRecV0[absLabel*20] =  iEntry;
	fMultiRecV0[absLabel]++;
      }
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
  fRecArray = new TObjArray(fNParticles);
  AliESDtrack dummytrack;  //

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
    //    if (fMCInfo->fNTPCRef<1) continue; // not TPCref
    //
    fRecInfo->Reset();
    AliESDtrack * track=0;
    fRecInfo->fReconstructed =0;
    TVector3 local = TR2Local(&(fMCInfo->fTrackRef),fParamTPC);
    local.GetXYZ(fRecInfo->fTRLocalCoord);	
    //
    if (fIndexRecTracks[fMCInfo->fLabel*20] >= 0) {
      //track= (AliESDtrack*)fTracks->UncheckedAt(fIndexRecTracks[fMCInfo->fLabel*4]);
      track= (AliESDtrack*)fEvent->GetTrack(fIndexRecTracks[fMCInfo->fLabel*20]);
      //
      //
      // find nearest track if multifound
      //Int_t sign = Int_t(track->GetSign()*fMCInfo->fCharge);
      //
      Int_t status = 0;
      if  ((track->GetStatus()&AliESDtrack::kITSrefit)>0) status++;
      if  ((track->GetStatus()&AliESDtrack::kTPCrefit)>0) status++;
      if  ((track->GetStatus()&AliESDtrack::kTRDrefit)>0) status++;

      //
      if (fIndexRecTracks[fMCInfo->fLabel*20+1]>0){
	//
	Double_t p[3];
	track->GetInnerPxPyPz(p);
	Float_t maxp = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
	//
	for (Int_t i=1;i<20;i++){
	  if (fIndexRecTracks[fMCInfo->fLabel*20+i]>=0){
	    AliESDtrack * track2 = (AliESDtrack*)fEvent->GetTrack(fIndexRecTracks[fMCInfo->fLabel*20+i]);
	    if (!track2) continue;
	    //Int_t sign2 = track->GetSign()*fMCInfo->fCharge; //	    
	    //if (sign2<0) continue;
	    track2->GetInnerPxPyPz(p);
	    Float_t mom = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
	    /*
	    if (sign<0){
	      sign = sign2;
	      track = track2;
	      maxp = mom;
	      continue;
	    }
	    */
	    //
	    Int_t status2 = 0;
	    if  ((track2->GetStatus()&AliESDtrack::kITSrefit)>0) status2++;
	    if  ((track2->GetStatus()&AliESDtrack::kTPCrefit)>0) status2++;
	    if  ((track2->GetStatus()&AliESDtrack::kTRDrefit)>0) status2++;
	    if (status2<status) continue;
	    //
	    if (mom<maxp) continue;
	    maxp = mom;
	    track = track2;
	    //
	  }
	}
      }	
      //
      if (track) {
	new (&(fRecInfo->fESDTrack)) AliESDtrack(*track);
      }else{
	fRecInfo->fESDTrack = dummytrack;
      }
      
      if (track->GetITStrack())
	fRecInfo->fITStrack = *((AliITStrackMI*)track->GetITStrack());
      else{
	fRecInfo->fITStrack = *track;
      }
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
      fRecInfo->Update(fMCInfo,fParamTPC,kTRUE, fEvent);          
    }
    else{
      fRecInfo->fESDTrack = dummytrack;
      fRecInfo->Update(fMCInfo,fParamTPC,kFALSE, fEvent);
    }
    fRecArray->AddAt(new AliESDRecInfo(*fRecInfo),fMCInfo->fLabel);
    fTreeCmp->Fill();
  }
  fTreeCmp->AutoSave();
  //fTracks->Delete();
  printf("Time spended in TreeGenLoop\n");
  timer.Print();
  if (fDebug > 2) cerr<<"end of TreeGenLoop"<<endl;

  return 0;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
Int_t ESDCmpTr::BuildKinkInfo0(Int_t eventNr)
{
//
// loop over all entries for a given event, find corresponding 
// rec. track and store in the fTreeCmp
//
  TStopwatch timer;
  timer.Start();
  Int_t entry = fNextKinkToRead;
  Double_t nParticlesTR = fTreeGenKinks->GetEntriesFast();
  cerr<<"fNParticles, nParticlesTR, fNextKinkToRead: "<<fNParticles<<" "
      <<nParticlesTR<<" "<<fNextKinkToRead<<endl;
  //
  TBranch * branch = fTreeCmpKinks->GetBranch("RC.");
  branch->SetAddress(&fRecKinkInfo); // set all pointers
  
  //
  while (entry < nParticlesTR) {
    fTreeGenKinks->GetEntry(entry);
    entry++;
    if (eventNr < fGenKinkInfo->fMCm.fEventNr) continue;
    if (eventNr > fGenKinkInfo->fMCm.fEventNr) continue;;
    //
    fNextKinkToRead = entry-1;
    //
    //
    AliESDRecInfo*  fRecInfo1 = (AliESDRecInfo*)fRecArray->At(fGenKinkInfo->fMCm.fLabel);
    AliESDRecInfo*  fRecInfo2 = (AliESDRecInfo*)fRecArray->At(fGenKinkInfo->fMCd.fLabel);
    fRecKinkInfo->fT1 = (*fRecInfo1);
    fRecKinkInfo->fT2 = (*fRecInfo2);
    fRecKinkInfo->fStatus =0;
    if (fRecInfo1 && fRecInfo1->fTPCOn) fRecKinkInfo->fStatus+=1;
    if (fRecInfo2 && fRecInfo2->fTPCOn) fRecKinkInfo->fStatus+=2;
    if (fRecKinkInfo->fStatus==3&&fRecInfo1->fSign!=fRecInfo2->fSign) fRecKinkInfo->fStatus*=-1;
    
    if (fRecKinkInfo->fStatus==3){
      fRecKinkInfo->Update();    
    }
    Int_t label =  TMath::Min(fGenKinkInfo->fMCm.fLabel,fGenKinkInfo->fMCd.fLabel);
    Int_t label2 = TMath::Max(fGenKinkInfo->fMCm.fLabel,fGenKinkInfo->fMCd.fLabel);
    
    AliESDkink *kink=0;
    fRecKinkInfo->fRecStatus   =0;
    fRecKinkInfo->fMultiple    = fMultiRecKinks[label];
    fRecKinkInfo->fKinkMultiple=0;
    //
    if (fMultiRecKinks[label]>0){

      //      for (Int_t j=0;j<TMath::Min(fMultiRecKinks[label],100);j++){
      for (Int_t j=TMath::Min(fMultiRecKinks[label],Short_t(20))-1;j>=0;j--){
	Int_t index = fIndexRecKinks[label*20+j];
	//AliESDkink *kink2  = (AliESDkink*)fKinks->At(index);
	AliESDkink *kink2  = (AliESDkink*)fEvent->GetKink(index);
	if (TMath::Abs(kink2->fLab[0])==label &&TMath::Abs(kink2->fLab[1])==label2) {
	  fRecKinkInfo->fKinkMultiple++;
	  fSignedKinks[index]=1;
	  Int_t c0=0;
	  if (kink){
	    //	    if (kink->fTRDOn) c0++;
	    //if (kink->fITSOn) c0++;
	    if (kink->GetStatus(2)>0) c0++;
	    if (kink->GetStatus(0)>0) c0++;
	  }
	  Int_t c2=0;
	  //	  if (kink2->fTRDOn) c2++;
	  //if (kink2->fITSOn) c2++;
	  if (kink2->GetStatus(2)>0) c2++;
	  if (kink2->GetStatus(0)>0) c2++;

	  if (c2<c0) continue;
	  kink =kink2;
	}
	if (TMath::Abs(kink2->fLab[1])==label &&TMath::Abs(kink2->fLab[0])==label2) {
	  fRecKinkInfo->fKinkMultiple++;
	  fSignedKinks[index]=1;
	  Int_t c0=0;
	  if (kink){
	    //if (kink->fTRDOn) c0++;
	    //if (kink->fITSOn) c0++;
	    if (kink->GetStatus(2)>0) c0++;
	    if (kink->GetStatus(0)>0) c0++;

	  }
	  Int_t c2=0;
	  //	  if (kink2->fTRDOn) c2++;
	  //if (kink2->fITSOn) c2++;
	  if (kink2->GetStatus(2)>0) c2++;
	  if (kink2->GetStatus(0)>0) c2++;

	  if (c2<c0) continue;
	  kink =kink2;
	}
      }
    }
    if (kink){
      fRecKinkInfo->fKink = *kink;
      fRecKinkInfo->fRecStatus=1;
    }
    fTreeCmpKinks->Fill();
  }
  //  Int_t nkinks = fKinks->GetEntriesFast();
  Int_t nkinks = fEvent->GetNumberOfKinks();
  for (Int_t i=0;i<nkinks;i++){
    if (fSignedKinks[i]==0){
      //      AliESDkink *kink  = (AliESDkink*)fKinks->At(i);
      AliESDkink *kink  = (AliESDkink*)fEvent->GetKink(i);
      if (!kink) continue;
      //
      fRecKinkInfo->fKink = *kink;
      fRecKinkInfo->fRecStatus =-2;
      //
      AliESDRecInfo*  fRecInfo1 = (AliESDRecInfo*)fRecArray->At(TMath::Abs(kink->fLab[0]));
      AliESDRecInfo*  fRecInfo2 = (AliESDRecInfo*)fRecArray->At(TMath::Abs(kink->fLab[1]));
      if (fRecInfo1 && fRecInfo2){
	fRecKinkInfo->fT1 = (*fRecInfo1);
	fRecKinkInfo->fT2 = (*fRecInfo2);
	fRecKinkInfo->fRecStatus =-1;
      }
      fTreeCmpKinks->Fill();
    }
  }


  fTreeCmpKinks->AutoSave();
  printf("Time spended in BuilKinkInfo Loop\n");
  timer.Print();
  if (fDebug > 2) cerr<<"end of BuildKinkInfo Loop"<<endl;
  return 0;
}



void   ESDCmpTr::MakePoints(AliESDtrack * track, AliPointsMI &points)
{
  //
  // make points in global coordinate frame
  //
  return;
  /*
  UInt_t itscl[10];
  Int_t nits = track->GetITSclusters(itscl);
  Int_t tpccl[1000];
  Int_t ntpc = track->GetTPcclusters(tpccl);
  UInt_t trdcl[1000];
  Int_t ntrd = track->GetTRDclusters(trdcl);
  //
  AliLoader *itsloader = fLoader->GetLoader("ITSLoader");
  AliLoader *tpcloader = fLoader->GetLoader("TPCLoader");
  AliLoader *trdloader = fLoader->GetLoader("TRDLoader");
  //
  AliITStrackerMI itstracker(); 
  */
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



Int_t ESDCmpTr::BuildV0Info(Int_t eventNr)
{
//
// loop over all entries for a given event, find corresponding 
// rec. track and store in the fTreeCmp
//
  TStopwatch timer;
  timer.Start();
  Int_t entry = fNextV0ToRead;
  Double_t nParticlesTR = fTreeGenV0->GetEntriesFast();
  cerr<<"fNParticles, nParticlesTR, fNextV0ToRead: "<<fNParticles<<" "
      <<nParticlesTR<<" "<<fNextV0ToRead<<endl;
  //
  TBranch * branch = fTreeCmpV0->GetBranch("RC.");
  branch->SetAddress(&fRecV0Info); // set all pointers
  const AliESDVertex * esdvertex = fEvent->GetVertex();
  Float_t vertex[3]= {esdvertex->GetXv(), esdvertex->GetYv(),esdvertex->GetZv()};
  
  //
  while (entry < nParticlesTR) {
    fTreeGenV0->GetEntry(entry);
    entry++;
    if (eventNr < fGenV0Info->fMCm.fEventNr) continue;
    if (eventNr > fGenV0Info->fMCm.fEventNr) continue;;
    //
    fNextV0ToRead = entry-1;
    //
    //
    AliESDRecInfo*  fRecInfo1 = (AliESDRecInfo*)fRecArray->At(fGenV0Info->fMCm.fLabel);
    AliESDRecInfo*  fRecInfo2 = (AliESDRecInfo*)fRecArray->At(fGenV0Info->fMCd.fLabel);
    if (fGenV0Info->fMCm.fCharge*fGenV0Info->fMCd.fCharge>0) continue;  // interactions
    if (!fRecInfo1 || !fRecInfo2) continue;
    fRecV0Info->fT1 = (*fRecInfo1);
    fRecV0Info->fT2 = (*fRecInfo2);
    fRecV0Info->fV0Status =0;
    if (fRecInfo1 && fRecInfo1->fStatus[1]>0) fRecV0Info->fV0Status+=1;
    if (fRecInfo2 && fRecInfo2->fStatus[1]>0) fRecV0Info->fV0Status+=2;

    if (fRecV0Info->fV0Status==3&&fRecInfo1->fSign==fRecInfo2->fSign) fRecV0Info->fV0Status*=-1;


    if (abs(fRecV0Info->fV0Status)==3){
      fRecV0Info->Update(vertex);
      {
	//
	// TPC V0 Info
	Double_t x,alpha, param[5],cov[15];
	if ( fRecV0Info->fT1.fESDTrack.GetInnerParam() && fRecV0Info->fT2.fESDTrack.GetInnerParam()){
	  fRecV0Info->fT1.fESDTrack.GetInnerExternalParameters(alpha,x,param);
	  fRecV0Info->fT1.fESDTrack.GetInnerExternalCovariance(cov);
	  AliExternalTrackParam paramP(x,alpha,param,cov);
	  //
	  fRecV0Info->fT2.fESDTrack.GetInnerExternalParameters(alpha,x,param);
	  fRecV0Info->fT2.fESDTrack.GetInnerExternalCovariance(cov);
	  AliExternalTrackParam paramM(x,alpha,param,cov);
	  //
	  fRecV0Info->fV0tpc.SetM(paramM);
	  fRecV0Info->fV0tpc.SetP(paramP);
	  Double_t pid1[5],pid2[5];
	  fRecV0Info->fT1.fESDTrack.GetESDpid(pid1);
	  fRecV0Info->fT1.fESDTrack.GetESDpid(pid2);
	  //
	  fRecV0Info->fV0tpc.UpdatePID(pid1,pid2);
	  fRecV0Info->fV0tpc.Update(vertex);
	
	  //
	  //
	  fRecV0Info->fT1.fESDTrack.GetExternalParameters(x,param);
	  fRecV0Info->fT1.fESDTrack.GetExternalCovariance(cov);
	  alpha = fRecV0Info->fT1.fESDTrack.GetAlpha();
	  new (&paramP) AliExternalTrackParam(x,alpha,param,cov);
	  //
	  fRecV0Info->fT2.fESDTrack.GetExternalParameters(x,param);
	  fRecV0Info->fT2.fESDTrack.GetExternalCovariance(cov);
	  alpha = fRecV0Info->fT2.fESDTrack.GetAlpha();
	  new (&paramM) AliExternalTrackParam(x,alpha,param,cov);
	  //
	  fRecV0Info->fV0its.SetM(paramM);
	  fRecV0Info->fV0its.SetP(paramP);
	  fRecV0Info->fV0its.UpdatePID(pid1,pid2);
	  fRecV0Info->fV0its.Update(vertex);
	}
      }
      if (TMath::Abs(fGenV0Info->fMCm.fPdg)==11 &&TMath::Abs(fGenV0Info->fMCd.fPdg)==11){
	if (fRecV0Info->fDist2>10){
	  fRecV0Info->Update(vertex);
	}
	if (fRecV0Info->fDist2>10){
	  fRecV0Info->Update(vertex);
	}
      }
    }   
    //
    // take the V0 from reconstruction
 
    Int_t label =  TMath::Min(fGenV0Info->fMCm.fLabel,fGenV0Info->fMCd.fLabel);
    Int_t label2 = TMath::Max(fGenV0Info->fMCm.fLabel,fGenV0Info->fMCd.fLabel);    
    AliESDv0 *v0MI=0;
    fRecV0Info->fRecStatus   =0;
    fRecV0Info->fMultiple    = fMultiRecV0[label];
    fRecV0Info->fV0Multiple=0;
    //
    if (fMultiRecV0[label]>0 || fMultiRecV0[label2]>0){

      //      for (Int_t j=0;j<TMath::Min(fMultiRecV0s[label],100);j++){
      for (Int_t j=TMath::Min(fMultiRecV0[label],Short_t(20))-1;j>=0;j--){
	Int_t index = fIndexRecV0[label*20+j];
	if (index<0) continue;
	AliESDv0 *v0MI2  = fEvent->GetV0(index);
	if (TMath::Abs(v0MI2->fLab[0])==label &&TMath::Abs(v0MI2->fLab[1])==label2) {
	  v0MI =v0MI2;
	  fRecV0Info->fV0Multiple++;
	  fSignedV0[index]=1;
	}
	if (TMath::Abs(v0MI2->fLab[1])==label &&TMath::Abs(v0MI2->fLab[0])==label2) {
	  v0MI =v0MI2;
	  fRecV0Info->fV0Multiple++;
	  fSignedV0[index]=1;
	}
      }
    }
    if (v0MI){
      fRecV0Info->fV0rec = *v0MI;
      fRecV0Info->fRecStatus=1;
    }

    fTreeCmpV0->Fill();
  }
  //
  // write fake v0s

  Int_t nV0MIs = fEvent->GetNumberOfV0s();
  for (Int_t i=0;i<nV0MIs;i++){
    if (fSignedV0[i]==0){
      AliESDv0 *v0MI  = fEvent->GetV0(i);
      if (!v0MI) continue;
      //
      fRecV0Info->fV0rec = *v0MI;
      fRecV0Info->fV0Status  =-10;
      fRecV0Info->fRecStatus =-2;
      //
      AliESDRecInfo*  fRecInfo1 = (AliESDRecInfo*)fRecArray->At(TMath::Abs(v0MI->fLab[0]));
      AliESDRecInfo*  fRecInfo2 = (AliESDRecInfo*)fRecArray->At(TMath::Abs(v0MI->fLab[1]));
      if (fRecInfo1 && fRecInfo2){
	fRecV0Info->fT1 = (*fRecInfo1);
	fRecV0Info->fT2 = (*fRecInfo2);
	fRecV0Info->fRecStatus =-1;
      }
      fRecV0Info->Update(vertex);
      fTreeCmpV0->Fill();
    }
  }



  fTreeCmpV0->AutoSave();
  printf("Time spended in BuilV0Info Loop\n");
  timer.Print();
  if (fDebug > 2) cerr<<"end of BuildV0Info Loop"<<endl;
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


