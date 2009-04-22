/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//      Implementation of the on-the-fly v0 finder for ITS tracker MI
//          Origin: Marian Ivanov, CERN, Marian.Ivanov@cern.ch 
//          Extraction from AliITStrackerMI: Andrea Dainese
//          Current support and development: 
//-------------------------------------------------------------------------

#include <TMatrixD.h>
#include <TTree.h>
#include <TString.h>
#include <TRandom.h>
#include <TTreeStream.h>


#include "AliLog.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliV0.h"
#include "AliHelix.h"
#include "AliITSRecPoint.h"
#include "AliITSReconstructor.h"
#include "AliITStrackerMI.h"
#include "AliITSV0Finder.h"

ClassImp(AliITSV0Finder)

AliITSV0Finder::AliITSV0Finder() {
  //Default constructor
}
//------------------------------------------------------------------------
void AliITSV0Finder::UpdateTPCV0(const AliESDEvent *event,
				   AliITStrackerMI *tracker)
{
  //
  //try to update, or reject TPC  V0s
  //
  TObjArray *trackHypothesys = tracker->GetTrackHypothesys();

  Int_t nv0s = event->GetNumberOfV0s();
  Int_t nitstracks = trackHypothesys->GetEntriesFast();

  for (Int_t i=0;i<nv0s;i++){
    AliESDv0 * vertex = event->GetV0(i);
    Int_t ip = vertex->GetIndex(0);
    Int_t im = vertex->GetIndex(1);
    //
    TObjArray * arrayp = (ip<nitstracks) ? (TObjArray*)trackHypothesys->At(ip):0;
    TObjArray * arraym = (im<nitstracks) ? (TObjArray*)trackHypothesys->At(im):0;
    AliITStrackMI * trackp = (arrayp!=0) ? (AliITStrackMI*)arrayp->At(0):0;
    AliITStrackMI * trackm = (arraym!=0) ? (AliITStrackMI*)arraym->At(0):0;
    //
    //
    if (trackp){
      if (trackp->GetNumberOfClusters()+trackp->GetNDeadZone()>5.5){
	if (trackp->GetConstrain()&&trackp->GetChi2MIP(0)<3) vertex->SetStatus(-100);
	if (!trackp->GetConstrain()&&trackp->GetChi2MIP(0)<2) vertex->SetStatus(-100); 
      }
    }

    if (trackm){
      if (trackm->GetNumberOfClusters()+trackm->GetNDeadZone()>5.5){
	if (trackm->GetConstrain()&&trackm->GetChi2MIP(0)<3) vertex->SetStatus(-100);
	if (!trackm->GetConstrain()&&trackm->GetChi2MIP(0)<2) vertex->SetStatus(-100); 
      }
    }
    if (vertex->GetStatus()==-100) continue;
    //
    Double_t xrp[3]; vertex->GetXYZ(xrp[0],xrp[1],xrp[2]);  //I.B.
    Int_t clayer = tracker->GetNearestLayer(xrp);                    //I.B.
    vertex->SetNBefore(clayer);        //
    vertex->SetChi2Before(9*clayer);   //
    vertex->SetNAfter(6-clayer);       //
    vertex->SetChi2After(0);           //
    //
    if (clayer >1 ){ // calculate chi2 before vertex
      Float_t chi2p = 0, chi2m=0;  
      //
      if (trackp){
	for (Int_t ilayer=0;ilayer<clayer;ilayer++){
	  if (trackp->GetClIndex(ilayer)>=0){
	    chi2p+=trackp->GetDy(ilayer)*trackp->GetDy(ilayer)/(trackp->GetSigmaY(ilayer)*trackp->GetSigmaY(ilayer))+
	      trackp->GetDz(ilayer)*trackp->GetDz(ilayer)/(trackp->GetSigmaZ(ilayer)*trackp->GetSigmaZ(ilayer));
	  }
	  else{
	    chi2p+=9;
	  }
	}
      }else{
	chi2p = 9*clayer;
      }
      //
      if (trackm){
	for (Int_t ilayer=0;ilayer<clayer;ilayer++){
	  if (trackm->GetClIndex(ilayer)>=0){
	    chi2m+=trackm->GetDy(ilayer)*trackm->GetDy(ilayer)/(trackm->GetSigmaY(ilayer)*trackm->GetSigmaY(ilayer))+
	      trackm->GetDz(ilayer)*trackm->GetDz(ilayer)/(trackm->GetSigmaZ(ilayer)*trackm->GetSigmaZ(ilayer));
	  }
	  else{
	    chi2m+=9;
	  }
	}
      }else{
	chi2m = 9*clayer;
      }
      vertex->SetChi2Before(TMath::Min(chi2p,chi2m));
      if (TMath::Min(chi2p,chi2m)/Float_t(clayer)<4) vertex->SetStatus(-10);  // track exist before vertex
    }
    
    if (clayer < 5 ){ // calculate chi2 after vertex
      Float_t chi2p = 0, chi2m=0;  
      //
      if (trackp&&TMath::Abs(trackp->GetTgl())<1.){
	for (Int_t ilayer=clayer;ilayer<6;ilayer++){
	  if (trackp->GetClIndex(ilayer)>=0){
	    chi2p+=trackp->GetDy(ilayer)*trackp->GetDy(ilayer)/(trackp->GetSigmaY(ilayer)*trackp->GetSigmaY(ilayer))+
	      trackp->GetDz(ilayer)*trackp->GetDz(ilayer)/(trackp->GetSigmaZ(ilayer)*trackp->GetSigmaZ(ilayer));
	  }
	  else{
	    chi2p+=9;
	  }
	}
      }else{
	chi2p = 0;
      }
      //
      if (trackm&&TMath::Abs(trackm->GetTgl())<1.){
	for (Int_t ilayer=clayer;ilayer<6;ilayer++){
	  if (trackm->GetClIndex(ilayer)>=0){
	    chi2m+=trackm->GetDy(ilayer)*trackm->GetDy(ilayer)/(trackm->GetSigmaY(ilayer)*trackm->GetSigmaY(ilayer))+
	      trackm->GetDz(ilayer)*trackm->GetDz(ilayer)/(trackm->GetSigmaZ(ilayer)*trackm->GetSigmaZ(ilayer));
	  }
	  else{
	    chi2m+=9;
	  }
	}
      }else{
	chi2m = 0;
      }
      vertex->SetChi2After(TMath::Max(chi2p,chi2m));
      if (TMath::Max(chi2m,chi2p)/Float_t(6-clayer)>9) vertex->SetStatus(-20);  // track not found in ITS
    }
  }
  //
}
//------------------------------------------------------------------------
void AliITSV0Finder::FindV02(AliESDEvent *event,
			       AliITStrackerMI *tracker) 
{
  //
  // V0 finder
  //
  //  Cuts on DCA -  R dependent
  //          max distance DCA between 2 tracks cut 
  //          maxDist = TMath::Min(kMaxDist,kMaxDist0+pvertex->GetRr()*kMaxDist);
  //
  const Float_t kMaxDist0      = 0.1;    
  const Float_t kMaxDist1      = 0.1;     
  const Float_t kMaxDist       = 1;
  const Float_t kMinPointAngle  = 0.85;
  const Float_t kMinPointAngle2 = 0.99;
  const Float_t kMinR           = 0.5;
  const Float_t kMaxR           = 220;
  //const Float_t kCausality0Cut   = 0.19;
  //const Float_t kLikelihood01Cut = 0.25;
  //const Float_t kPointAngleCut   = 0.9996;
  const Float_t kCausality0Cut   = 0.19;
  const Float_t kLikelihood01Cut = 0.45;
  const Float_t kLikelihood1Cut  = 0.5;
  const Float_t kCombinedCut     = 0.55;

  TObjArray *trackHypothesys = tracker->GetTrackHypothesys();
  TObjArray *bestHypothesys  = tracker->GetBestHypothesys();
  TObjArray *originalTracks  = tracker->GetOriginal();
  //
  //
  TTreeSRedirector &cstream = *(tracker->GetDebugStreamer());
  Int_t ntracks    = event->GetNumberOfTracks(); 
  Int_t nitstracks = trackHypothesys->GetEntriesFast();
  originalTracks->Expand(ntracks);
  trackHypothesys->Expand(ntracks);
  bestHypothesys->Expand(ntracks);
  //
  AliHelix * helixes   = new AliHelix[ntracks+2];
  TObjArray trackarray(ntracks+2);     //array with tracks - with vertex constrain
  TObjArray trackarrayc(ntracks+2);    //array of "best    tracks" - without vertex constrain
  TObjArray trackarrayl(ntracks+2);    //array of "longest tracks" - without vertex constrain
  Bool_t * forbidden   = new Bool_t [ntracks+2];
  Int_t   *itsmap      = new Int_t  [ntracks+2];
  Float_t *dist        = new Float_t[ntracks+2];
  Float_t *normdist0   = new Float_t[ntracks+2];
  Float_t *normdist1   = new Float_t[ntracks+2];
  Float_t *normdist    = new Float_t[ntracks+2];
  Float_t *norm        = new Float_t[ntracks+2];
  Float_t *maxr        = new Float_t[ntracks+2];
  Float_t *minr        = new Float_t[ntracks+2];
  Float_t *minPointAngle= new Float_t[ntracks+2];
  //
  AliV0 *pvertex      = new AliV0;
  AliITStrackMI * dummy= new AliITStrackMI;
  dummy->SetLabel(0);
  AliITStrackMI  trackat0;    //temporary track for DCA calculation
  //
  Float_t primvertex[3]={tracker->GetX(),tracker->GetY(),tracker->GetZ()};
  //
  // make ITS -  ESD map
  //
  for (Int_t itrack=0;itrack<ntracks+2;itrack++) {
    itsmap[itrack]        = -1;
    forbidden[itrack]     = kFALSE;
    maxr[itrack]          = kMaxR;
    minr[itrack]          = kMinR;
    minPointAngle[itrack] = kMinPointAngle;
  }
  for (Int_t itrack=0;itrack<nitstracks;itrack++){
    AliITStrackMI * original =   (AliITStrackMI*)(originalTracks->At(itrack));
    Int_t           esdindex =   original->GetESDtrack()->GetID();
    itsmap[esdindex]         =   itrack;
  }
  //
  // create ITS tracks from ESD tracks if not done before
  //
  for (Int_t itrack=0;itrack<ntracks;itrack++){
    if (itsmap[itrack]>=0) continue;
    AliITStrackMI * tpctrack = new AliITStrackMI(*(event->GetTrack(itrack)));
    //tpctrack->fD[0] = tpctrack->GetD(GetX(),GetY());
    //tpctrack->fD[1] = tpctrack->GetZat(GetX())-GetZ(); 
    tpctrack->GetDZ(tracker->GetX(),tracker->GetY(),tracker->GetZ(),tpctrack->GetDP());   //I.B.
    if (tpctrack->GetD(0)<20 && tpctrack->GetD(1)<20){
      // tracks which can reach inner part of ITS
      // propagate track to outer its volume - with correction for material
      AliITStrackerMI::CorrectForTPCtoITSDeadZoneMaterial(tpctrack);  
    }
    itsmap[itrack] = nitstracks;
    originalTracks->AddAt(tpctrack,nitstracks);
    nitstracks++;
  }
  //
  // fill temporary arrays
  //
  for (Int_t itrack=0;itrack<ntracks;itrack++){
    AliESDtrack *  esdtrack = event->GetTrack(itrack);
    Int_t          itsindex = itsmap[itrack];
    AliITStrackMI *original = (AliITStrackMI*)originalTracks->At(itsmap[itrack]);
    if (!original) continue;
    AliITStrackMI *bestConst  = 0;
    AliITStrackMI *bestLong   = 0;
    AliITStrackMI *best       = 0;    
    //
    //
    TObjArray * array    = (TObjArray*)  trackHypothesys->At(itsindex);
    Int_t       hentries = (array==0) ?  0 : array->GetEntriesFast();
    // Get best track with vertex constrain
    for (Int_t ih=0;ih<hentries;ih++){
      AliITStrackMI * trackh = (AliITStrackMI*)array->At(ih);
      if (!trackh->GetConstrain()) continue;
      if (!bestConst) bestConst = trackh;
      if (trackh->GetNumberOfClusters()>5.0){
	bestConst  = trackh;                         // full track -  with minimal chi2
	break;
      }
      if (trackh->GetNumberOfClusters()+trackh->GetNDeadZone()<=bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone())  continue;      
      bestConst = trackh;
      break;
    }
    // Get best long track without vertex constrain and best track without vertex constrain
    for (Int_t ih=0;ih<hentries;ih++){
      AliITStrackMI * trackh = (AliITStrackMI*)array->At(ih);
      if (trackh->GetConstrain()) continue;
      if (!best)     best     = trackh;
      if (!bestLong) bestLong = trackh;
      if (trackh->GetNumberOfClusters()>5.0){
	bestLong  = trackh;                         // full track -  with minimal chi2
	break;
      }
      if (trackh->GetNumberOfClusters()+trackh->GetNDeadZone()<=bestLong->GetNumberOfClusters()+bestLong->GetNDeadZone())  continue;      
      bestLong = trackh;	
    }
    if (!best) {
      best     = original;
      bestLong = original;
    }
    //I.B. trackat0 = *bestLong;
    new (&trackat0) AliITStrackMI(*bestLong);
    Double_t xx,yy,zz,alpha; 
    if (!bestLong->GetGlobalXYZat(bestLong->GetX(),xx,yy,zz)) continue;
    alpha = TMath::ATan2(yy,xx);    
    if (!trackat0.Propagate(alpha,0)) continue;      
    // calculate normalized distances to the vertex 
    //
    Float_t ptfac  = (1.+100.*TMath::Abs(trackat0.GetC()));
    if ( bestLong->GetNumberOfClusters()>3 ){      
      dist[itsindex]      = trackat0.GetY();
      norm[itsindex]      = ptfac*TMath::Sqrt(trackat0.GetSigmaY2());
      normdist0[itsindex] = TMath::Abs(trackat0.GetY()/norm[itsindex]);
      normdist1[itsindex] = TMath::Abs((trackat0.GetZ()-primvertex[2])/(ptfac*TMath::Sqrt(trackat0.GetSigmaZ2())));
      normdist[itsindex]  = TMath::Sqrt(normdist0[itsindex]*normdist0[itsindex]+normdist1[itsindex]*normdist1[itsindex]);
      if (!bestConst){
	if (bestLong->GetNumberOfClusters()+bestLong->GetNDeadZone()<6) normdist[itsindex]*=2.;
	if (bestLong->GetNumberOfClusters()+bestLong->GetNDeadZone()<5) normdist[itsindex]*=2.;
	if (bestLong->GetNumberOfClusters()+bestLong->GetNDeadZone()<4) normdist[itsindex]*=2.;
      }else{
	if (bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()<6) normdist[itsindex]*=1.5;
	if (bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()<5) normdist[itsindex]*=1.5;
      }
    }
    else{      
      if (bestConst&&bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()>4.5){
	dist[itsindex] = bestConst->GetD(0);
	norm[itsindex] = bestConst->GetDnorm(0);
	normdist0[itsindex] = TMath::Abs(bestConst->GetD(0)/norm[itsindex]);
	normdist1[itsindex] = TMath::Abs(bestConst->GetD(0)/norm[itsindex]);
	normdist[itsindex]  = TMath::Sqrt(normdist0[itsindex]*normdist0[itsindex]+normdist1[itsindex]*normdist1[itsindex]);
      }else{
	dist[itsindex]      = trackat0.GetY();
	norm[itsindex]      = ptfac*TMath::Sqrt(trackat0.GetSigmaY2());
	normdist0[itsindex] = TMath::Abs(trackat0.GetY()/norm[itsindex]);
	normdist1[itsindex] = TMath::Abs((trackat0.GetZ()-primvertex[2])/(ptfac*TMath::Sqrt(trackat0.GetSigmaZ2())));
	normdist[itsindex]  = TMath::Sqrt(normdist0[itsindex]*normdist0[itsindex]+normdist1[itsindex]*normdist1[itsindex]);
	if (TMath::Abs(trackat0.GetTgl())>1.05){
	  if (normdist[itsindex]<3) forbidden[itsindex]=kTRUE;
	  if (normdist[itsindex]>3) {
	    minr[itsindex] = TMath::Max(Float_t(40.),minr[itsindex]);
	  }
	}
      }
    }
    //
    //-----------------------------------------------------------
    //Forbid primary track candidates - 
    //
    //treetr->SetAlias("forbidden0","Tr0.fN<4&&Tr1.fN+Tr1.fNDeadZone>4.5");
    //treetr->SetAlias("forbidden1","ND<3&&Tr1.fN+Tr1.fNDeadZone>5.5");
    //treetr->SetAlias("forbidden2","ND<2&&Tr1.fClIndex[0]>0&&Tr1.fClIndex[0]>0");
    //treetr->SetAlias("forbidden3","ND<1&&Tr1.fClIndex[0]>0");
    //treetr->SetAlias("forbidden4","ND<4&&Tr1.fNormChi2[0]<2");
    //treetr->SetAlias("forbidden5","ND<5&&Tr1.fNormChi2[0]<1");
    //-----------------------------------------------------------
    if (bestConst){
      if (bestLong->GetNumberOfClusters()<4       && bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()>4.5)               forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<3 && bestConst->GetNumberOfClusters()+bestConst->GetNDeadZone()>5.5)               forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<2 && bestConst->GetClIndex(0)>=0 && bestConst->GetClIndex(1)>=0 ) forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<1 && bestConst->GetClIndex(0)>=0)                              forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<4 && bestConst->GetNormChi2(0)<2)                             forbidden[itsindex]=kTRUE;
      if (normdist[itsindex]<5 && bestConst->GetNormChi2(0)<1)                             forbidden[itsindex]=kTRUE;      
      if (bestConst->GetNormChi2(0)<2.5) {
	minPointAngle[itsindex]= 0.9999;
	maxr[itsindex]         = 10;
      }
    }
    //
    //forbid daughter kink candidates
    //
    if (esdtrack->GetKinkIndex(0)>0) forbidden[itsindex] = kTRUE;
    Bool_t isElectron = kTRUE;
    Bool_t isProton   = kTRUE;
    Double_t pid[5];
    esdtrack->GetESDpid(pid);
    for (Int_t i=1;i<5;i++){
      if (pid[0]<pid[i]) isElectron= kFALSE;
      if (pid[4]<pid[i]) isProton= kFALSE;
    }
    if (isElectron){
      forbidden[itsindex]=kFALSE;	
      normdist[itsindex]*=-1;
    }
    if (isProton){
      if (normdist[itsindex]>2) forbidden[itsindex]=kFALSE;	
      normdist[itsindex]*=-1;
    }

    //
    // Causality cuts in TPC volume
    //
    if (esdtrack->GetTPCdensity(0,10) >0.6)  maxr[itsindex] = TMath::Min(Float_t(110),maxr[itsindex]);
    if (esdtrack->GetTPCdensity(10,30)>0.6)  maxr[itsindex] = TMath::Min(Float_t(120),maxr[itsindex]);
    if (esdtrack->GetTPCdensity(20,40)>0.6)  maxr[itsindex] = TMath::Min(Float_t(130),maxr[itsindex]);
    if (esdtrack->GetTPCdensity(30,50)>0.6)  maxr[itsindex] = TMath::Min(Float_t(140),maxr[itsindex]);
    //
    if (esdtrack->GetTPCdensity(0,60)<0.4&&bestLong->GetNumberOfClusters()<3) minr[itsindex]=100;    
    //
    //
    if (kFALSE){
      cstream<<"Track"<<
	"Tr0.="<<best<<
	"Tr1.="<<((bestConst)? bestConst:dummy)<<
	"Tr2.="<<bestLong<<
	"Tr3.="<<&trackat0<<
	"Esd.="<<esdtrack<<
	"Dist="<<dist[itsindex]<<
	"ND0="<<normdist0[itsindex]<<
	"ND1="<<normdist1[itsindex]<<
	"ND="<<normdist[itsindex]<<
	"Pz="<<primvertex[2]<<
	"Forbid="<<forbidden[itsindex]<<
	"\n";
      //
    }
    trackarray.AddAt(best,itsindex);
    trackarrayc.AddAt(bestConst,itsindex);
    trackarrayl.AddAt(bestLong,itsindex);
    new (&helixes[itsindex]) AliHelix(*best);
  }
  //
  //
  //
  // first iterration of V0 finder
  //
  for (Int_t iesd0=0;iesd0<ntracks;iesd0++){
    Int_t itrack0 = itsmap[iesd0];
    if (forbidden[itrack0]) continue;
    AliITStrackMI * btrack0 = (AliITStrackMI*)trackarray.At(itrack0);
    if (!btrack0) continue;    
    if (btrack0->GetSign()>0 && !AliITSReconstructor::GetRecoParam()->GetStoreLikeSignV0s()) continue;
    AliITStrackMI *trackc0 = (AliITStrackMI*)trackarrayc.At(itrack0);
    //
    for (Int_t iesd1=iesd0+1;iesd1<ntracks;iesd1++){
      Int_t itrack1 = itsmap[iesd1];
      if (forbidden[itrack1]) continue;

      AliITStrackMI * btrack1 = (AliITStrackMI*)trackarray.At(itrack1); 
      if (!btrack1) continue;
      if (btrack1->GetSign()<0 && !AliITSReconstructor::GetRecoParam()->GetStoreLikeSignV0s()) continue;
      Bool_t isGold = kFALSE;
      if (TMath::Abs(TMath::Abs(btrack0->GetLabel())-TMath::Abs(btrack1->GetLabel()))==1){
	isGold = kTRUE;
      }
      AliITStrackMI *trackc1 = (AliITStrackMI*)trackarrayc.At(itrack1);
      AliHelix &h1 = helixes[itrack0];
      AliHelix &h2 = helixes[itrack1];
      //
      // find linear distance
      Double_t rmin =0;
      //
      //
      //
      Double_t phase[2][2],radius[2];
      Int_t  points = h1.GetRPHIintersections(h2, phase, radius);
      if    (points==0)  continue;
      Double_t delta[2]={1000000,1000000};        
      rmin = radius[0];
      h1.ParabolicDCA(h2,phase[0][0],phase[0][1],radius[0],delta[0]);
      if (points==2){    
	if (radius[1]<rmin) rmin = radius[1];
	h1.ParabolicDCA(h2,phase[1][0],phase[1][1],radius[1],delta[1]);
      }
      rmin = TMath::Sqrt(rmin);
      Double_t distance = 0;
      Double_t radiusC  = 0;
      Int_t    iphase   = 0;
      if (points==1 || delta[0]<delta[1]){
	distance = TMath::Sqrt(delta[0]);
	radiusC  = TMath::Sqrt(radius[0]);
      }else{
	distance = TMath::Sqrt(delta[1]);
	radiusC  = TMath::Sqrt(radius[1]);
	iphase=1;
      }
      if (radiusC<TMath::Max(minr[itrack0],minr[itrack1]))    continue;
      if (radiusC>TMath::Min(maxr[itrack0],maxr[itrack1]))     continue; 
      Float_t maxDist  = TMath::Min(kMaxDist,Float_t(kMaxDist0+radiusC*kMaxDist1));      
      if (distance>maxDist) continue;
      Float_t pointAngle = h1.GetPointAngle(h2,phase[iphase],primvertex);
      if (pointAngle<TMath::Max(minPointAngle[itrack0],minPointAngle[itrack1])) continue;
      //
      //
      //      Double_t distance = TestV0(h1,h2,pvertex,rmin);
      //
      //       if (distance>maxDist)           continue;
      //       if (pvertex->GetRr()<kMinR)     continue;
      //       if (pvertex->GetRr()>kMaxR)     continue;
      AliITStrackMI * track0=btrack0;
      AliITStrackMI * track1=btrack1;
      //      if (pvertex->GetRr()<3.5){  
      if (radiusC<3.5){  
	//use longest tracks inside the pipe
	track0 = (AliITStrackMI*)trackarrayl.At(itrack0);
	track1 = (AliITStrackMI*)trackarrayl.At(itrack1);
      }      
      //
      //
      pvertex->SetParamN(*track0);
      pvertex->SetParamP(*track1);
      pvertex->Update(primvertex);
      pvertex->SetClusters(track0->ClIndex(),track1->ClIndex());  // register clusters

      if (pvertex->GetRr()<kMinR) continue;
      if (pvertex->GetRr()>kMaxR) continue;
      if (pvertex->GetV0CosineOfPointingAngle()<kMinPointAngle) continue;
//Bo:      if (pvertex->GetDist2()>maxDist) continue;
      if (pvertex->GetDcaV0Daughters()>maxDist) continue;
//Bo:        pvertex->SetLab(0,track0->GetLabel());
//Bo:        pvertex->SetLab(1,track1->GetLabel());
      pvertex->SetIndex(0,track0->GetESDtrack()->GetID());
      pvertex->SetIndex(1,track1->GetESDtrack()->GetID());
      //      
      AliITStrackMI * htrackc0 = trackc0 ? trackc0:dummy;      
      AliITStrackMI * htrackc1 = trackc1 ? trackc1:dummy;

      //
      //
      TObjArray * array0b     = (TObjArray*)bestHypothesys->At(itrack0);
      if (!array0b&&pvertex->GetRr()<40 && TMath::Abs(track0->GetTgl())<1.1) {
	tracker->SetCurrentEsdTrack(itrack0);
	tracker->FollowProlongationTree((AliITStrackMI*)originalTracks->At(itrack0),itrack0, kFALSE);
      }
      TObjArray * array1b    = (TObjArray*)bestHypothesys->At(itrack1);
      if (!array1b&&pvertex->GetRr()<40 && TMath::Abs(track1->GetTgl())<1.1) { 
	tracker->SetCurrentEsdTrack(itrack1);
	tracker->FollowProlongationTree((AliITStrackMI*)originalTracks->At(itrack1),itrack1, kFALSE);
      }
      //
      AliITStrackMI * track0b = (AliITStrackMI*)originalTracks->At(itrack0);
      AliITStrackMI * track1b = (AliITStrackMI*)originalTracks->At(itrack1);
      AliITStrackMI * track0l = (AliITStrackMI*)originalTracks->At(itrack0);
      AliITStrackMI * track1l = (AliITStrackMI*)originalTracks->At(itrack1);
      
      Float_t minchi2before0=16;
      Float_t minchi2before1=16;
      Float_t minchi2after0 =16;
      Float_t minchi2after1 =16;
      Double_t xrp[3]; pvertex->GetXYZ(xrp[0],xrp[1],xrp[2]);  //I.B.
      Int_t maxLayer = tracker->GetNearestLayer(xrp);                   //I.B.
      
      if (array0b) for (Int_t i=0;i<5;i++){
	// best track after vertex
	AliITStrackMI * btrack = (AliITStrackMI*)array0b->At(i);
	if (!btrack) continue;
	if (btrack->GetNumberOfClusters()>track0l->GetNumberOfClusters()) track0l = btrack;     
	//	if (btrack->fX<pvertex->GetRr()-2.-0.5/(0.1+pvertex->GetAnglep()[2])) {
	if (btrack->GetX()<pvertex->GetRr()-2.) {
	  if ( (maxLayer>i+2|| (i==0)) && btrack->GetNumberOfClusters()==(6-i)&&i<3){
	    Float_t sumchi2= 0;
	    Float_t sumn   = 0;
	    if (maxLayer<3){   // take prim vertex as additional measurement
	      if (normdist[itrack0]>0 && htrackc0){
		sumchi2 += TMath::Min((3.-maxLayer)*normdist[itrack0]*normdist[itrack0],16.);
	      }else{
		sumchi2 +=  TMath::Min((3.-maxLayer)*(3*normdist[itrack0]*normdist[itrack0]+3.),16.);
	      }
	      sumn    +=  3-maxLayer;
	    }
	    for (Int_t ilayer=i;ilayer<maxLayer;ilayer++){
	      sumn+=1.;	      
	      if (btrack->GetClIndex(ilayer)<0){
		sumchi2+=25;
		continue;
	      }else{
		Int_t c=( btrack->GetClIndex(ilayer) & 0x0fffffff);
		for (Int_t itrack=0;itrack<4;itrack++){
		  AliITStrackerMI::AliITSlayer &layer=tracker->GetLayer(ilayer);
		  if (layer.GetClusterTracks(itrack,c)>=0 && layer.GetClusterTracks(itrack,c)!=itrack0){
		    sumchi2+=18.;  //shared cluster
		    break;
		  }
		}
		sumchi2+=btrack->GetDy(ilayer)*btrack->GetDy(ilayer)/(btrack->GetSigmaY(ilayer)*btrack->GetSigmaY(ilayer));
		sumchi2+=btrack->GetDz(ilayer)*btrack->GetDz(ilayer)/(btrack->GetSigmaZ(ilayer)*btrack->GetSigmaZ(ilayer));	       
	      }
	    }
	    sumchi2/=sumn;
	    if (sumchi2<minchi2before0) minchi2before0=sumchi2; 
	  }
	  continue;   //safety space - Geo manager will give exact layer
	}
	track0b       = btrack;
	minchi2after0 = btrack->GetNormChi2(i);
	break;
      }
      if (array1b) for (Int_t i=0;i<5;i++){
	// best track after vertex
	AliITStrackMI * btrack = (AliITStrackMI*)array1b->At(i);
	if (!btrack) continue;
	if (btrack->GetNumberOfClusters()>track1l->GetNumberOfClusters()) track1l = btrack;     
	//	if (btrack->fX<pvertex->GetRr()-2-0.5/(0.1+pvertex->GetAnglep()[2])){
	if (btrack->GetX()<pvertex->GetRr()-2){
	  if ((maxLayer>i+2 || (i==0))&&btrack->GetNumberOfClusters()==(6-i)&&(i<3)){
	    Float_t sumchi2= 0;
	    Float_t sumn   = 0;
	    if (maxLayer<3){   // take prim vertex as additional measurement
	      if (normdist[itrack1]>0 && htrackc1){
		sumchi2 +=  TMath::Min((3.-maxLayer)*normdist[itrack1]*normdist[itrack1],16.);
	      }else{
		sumchi2 += TMath::Min((3.-maxLayer)*(3*normdist[itrack1]*normdist[itrack1]+3.),16.);
	      }
	      sumn    +=  3-maxLayer;
	    }
	    for (Int_t ilayer=i;ilayer<maxLayer;ilayer++){
	      sumn+=1.;
	      if (btrack->GetClIndex(ilayer)<0){
		sumchi2+=30;
		continue;
	      }else{
		Int_t c=( btrack->GetClIndex(ilayer) & 0x0fffffff);
		for (Int_t itrack=0;itrack<4;itrack++){
		  AliITStrackerMI::AliITSlayer &layer=tracker->GetLayer(ilayer);
		  if (layer.GetClusterTracks(itrack,c)>=0 && layer.GetClusterTracks(itrack,c)!=itrack1){
		    sumchi2+=18.;  //shared cluster
		    break;
		  }
		}
		sumchi2+=btrack->GetDy(ilayer)*btrack->GetDy(ilayer)/(btrack->GetSigmaY(ilayer)*btrack->GetSigmaY(ilayer));
		sumchi2+=btrack->GetDz(ilayer)*btrack->GetDz(ilayer)/(btrack->GetSigmaZ(ilayer)*btrack->GetSigmaZ(ilayer));	       
	      }
	    }
	    sumchi2/=sumn;
	    if (sumchi2<minchi2before1) minchi2before1=sumchi2; 
	  }
	  continue;   //safety space - Geo manager will give exact layer	   
	}
	track1b = btrack;
	minchi2after1 = btrack->GetNormChi2(i);
	break;
      }
      //
      // position resolution - used for DCA cut
      Float_t sigmad = track0b->GetSigmaY2()+track0b->GetSigmaZ2()+track1b->GetSigmaY2()+track1b->GetSigmaZ2()+
	(track0b->GetX()-pvertex->GetRr())*(track0b->GetX()-pvertex->GetRr())*(track0b->GetSigmaSnp2()+track0b->GetSigmaTgl2())+
	(track1b->GetX()-pvertex->GetRr())*(track1b->GetX()-pvertex->GetRr())*(track1b->GetSigmaSnp2()+track1b->GetSigmaTgl2());
      sigmad =TMath::Sqrt(sigmad)+0.04;
      if (pvertex->GetRr()>50){
	Double_t cov0[15],cov1[15];
	track0b->GetESDtrack()->GetInnerExternalCovariance(cov0);
	track1b->GetESDtrack()->GetInnerExternalCovariance(cov1);
	sigmad = cov0[0]+cov0[2]+cov1[0]+cov1[2]+
	  (80.-pvertex->GetRr())*(80.-pvertex->GetRr())*(cov0[5]+cov0[9])+
	  (80.-pvertex->GetRr())*(80.-pvertex->GetRr())*(cov1[5]+cov1[9]);
	sigmad =TMath::Sqrt(sigmad)+0.05;
      }
      //       
      AliV0 vertex2;
      vertex2.SetParamN(*track0b);
      vertex2.SetParamP(*track1b);
      vertex2.Update(primvertex);
      //Bo:      if (vertex2.GetDist2()<=pvertex->GetDist2()&&(vertex2.GetV0CosineOfPointingAngle()>=pvertex->GetV0CosineOfPointingAngle())){
      if (vertex2.GetDcaV0Daughters()<=pvertex->GetDcaV0Daughters()&&(vertex2.GetV0CosineOfPointingAngle()>=pvertex->GetV0CosineOfPointingAngle())){
	pvertex->SetParamN(*track0b);
	pvertex->SetParamP(*track1b);
	pvertex->Update(primvertex);
	pvertex->SetClusters(track0b->ClIndex(),track1b->ClIndex());  // register clusters
	pvertex->SetIndex(0,track0->GetESDtrack()->GetID());
	pvertex->SetIndex(1,track1->GetESDtrack()->GetID());
      }
      pvertex->SetDistSigma(sigmad);
      //Bo:      pvertex->SetDistNorm(pvertex->GetDist2()/sigmad);       
      pvertex->SetNormDCAPrim(normdist[itrack0],normdist[itrack1]);
      //
      // define likelihhod and causalities
      //
      Float_t pa0=1, pa1=1, pb0=0.26, pb1=0.26;      
      if (maxLayer<1){
	Float_t fnorm0 = normdist[itrack0];
	if (fnorm0<0) fnorm0*=-3;
	Float_t fnorm1 = normdist[itrack1];
	if (fnorm1<0) fnorm1*=-3;
 	if ((pvertex->GetAnglep()[2]>0.1) || ( (pvertex->GetRr()<10.5)&& pvertex->GetAnglep()[2]>0.05 ) || (pvertex->GetRr()<3)){
 	  pb0    =  TMath::Exp(-TMath::Min(fnorm0,Float_t(16.))/12.);
 	  pb1    =  TMath::Exp(-TMath::Min(fnorm1,Float_t(16.))/12.);
 	}
	pvertex->SetChi2Before(normdist[itrack0]);
	pvertex->SetChi2After(normdist[itrack1]);       
	pvertex->SetNAfter(0);
	pvertex->SetNBefore(0);
      }else{
	pvertex->SetChi2Before(minchi2before0);
	pvertex->SetChi2After(minchi2before1);
	 if (pvertex->GetAnglep()[2]>0.1 || ( pvertex->GetRr()<10.5 && pvertex->GetAnglep()[2]>0.05) || pvertex->GetRr()<3){
	   pb0    =  TMath::Exp(-TMath::Min(minchi2before0,Float_t(16))/12.);
	   pb1    =  TMath::Exp(-TMath::Min(minchi2before1,Float_t(16))/12.);
	 }
	 pvertex->SetNAfter(maxLayer);
	 pvertex->SetNBefore(maxLayer);      
      }
      if (pvertex->GetRr()<90){
	pa0  *= TMath::Min(track0->GetESDtrack()->GetTPCdensity(0,60),Double_t(1.));
	pa1  *= TMath::Min(track1->GetESDtrack()->GetTPCdensity(0,60),Double_t(1.));
      }
      if (pvertex->GetRr()<20){
	pa0  *= (0.2+TMath::Exp(-TMath::Min(minchi2after0,Float_t(16))/8.))/1.2;
	pa1  *= (0.2+TMath::Exp(-TMath::Min(minchi2after1,Float_t(16))/8.))/1.2;
      }
      //
      pvertex->SetCausality(pb0,pb1,pa0,pa1);
      //
      //  Likelihood calculations  - apply cuts
      //         
      Bool_t v0OK = kTRUE;
      Float_t p12 = pvertex->GetParamP()->GetParameter()[4]*pvertex->GetParamP()->GetParameter()[4];
      p12        += pvertex->GetParamN()->GetParameter()[4]*pvertex->GetParamN()->GetParameter()[4];
      p12         = TMath::Sqrt(p12);                             // "mean" momenta
      Float_t    sigmap0   = 0.0001+0.001/(0.1+pvertex->GetRr()); 
      Float_t    sigmap    = 0.5*sigmap0*(0.6+0.4*p12);           // "resolution: of point angle - as a function of radius and momenta

      Float_t causalityA  = (1.0-pvertex->GetCausalityP()[0])*(1.0-pvertex->GetCausalityP()[1]);
      Float_t causalityB  = TMath::Sqrt(TMath::Min(pvertex->GetCausalityP()[2],Double_t(0.7))*
					TMath::Min(pvertex->GetCausalityP()[3],Double_t(0.7)));
      //
      //Bo:      Float_t likelihood0 = (TMath::Exp(-pvertex->GetDistNorm())+0.1) *(pvertex->GetDist2()<0.5)*(pvertex->GetDistNorm()<5);
      Float_t lDistNorm = pvertex->GetDcaV0Daughters()/pvertex->GetDistSigma();
      Float_t likelihood0 = (TMath::Exp(-lDistNorm)+0.1) *(pvertex->GetDcaV0Daughters()<0.5)*(lDistNorm<5);

      Float_t likelihood1 = TMath::Exp(-(1.0001-pvertex->GetV0CosineOfPointingAngle())/sigmap)+
	0.4*TMath::Exp(-(1.0001-pvertex->GetV0CosineOfPointingAngle())/(4.*sigmap))+
	0.4*TMath::Exp(-(1.0001-pvertex->GetV0CosineOfPointingAngle())/(8.*sigmap))+
	0.1*TMath::Exp(-(1.0001-pvertex->GetV0CosineOfPointingAngle())/0.01);
      //
      if (causalityA<kCausality0Cut)                                          v0OK = kFALSE;
      if (TMath::Sqrt(likelihood0*likelihood1)<kLikelihood01Cut)              v0OK = kFALSE;
      if (likelihood1<kLikelihood1Cut)                                        v0OK = kFALSE;
      if (TMath::Power(likelihood0*likelihood1*causalityB,0.33)<kCombinedCut) v0OK = kFALSE;
      
      //
      //
      if (kFALSE){	
	Bool_t gold = TMath::Abs(TMath::Abs(track0->GetLabel())-TMath::Abs(track1->GetLabel()))==1;
	cstream<<"It0"<<
	  "Tr0.="<<track0<<                       //best without constrain
	  "Tr1.="<<track1<<                       //best without constrain  
	  "Tr0B.="<<track0b<<                     //best without constrain  after vertex
	  "Tr1B.="<<track1b<<                     //best without constrain  after vertex 
	  "Tr0C.="<<htrackc0<<                    //best with constrain     if exist
	  "Tr1C.="<<htrackc1<<                    //best with constrain     if exist
	  "Tr0L.="<<track0l<<                     //longest best           
	  "Tr1L.="<<track1l<<                     //longest best
	  "Esd0.="<<track0->GetESDtrack()<<           // esd track0 params
	  "Esd1.="<<track1->GetESDtrack()<<           // esd track1 params
	  "V0.="<<pvertex<<                       //vertex properties
	  "V0b.="<<&vertex2<<                       //vertex properties at "best" track
	  "ND0="<<normdist[itrack0]<<             //normalize distance for track0
	  "ND1="<<normdist[itrack1]<<             //normalize distance for track1
	  "Gold="<<gold<<                         //
	  //	  "RejectBase="<<rejectBase<<             //rejection in First itteration
	  "OK="<<v0OK<<
	  "rmin="<<rmin<<
	  "sigmad="<<sigmad<<
	  "\n";
      }      
      //if (rejectBase) continue;
      //
      pvertex->SetStatus(0);
      //      if (rejectBase) {
      //	pvertex->SetStatus(-100);
      //}
      if (pvertex->GetV0CosineOfPointingAngle()>kMinPointAngle2) {
	//Bo:	  pvertex->SetESDindexes(track0->GetESDtrack()->GetID(),track1->GetESDtrack()->GetID());
	pvertex->SetIndex(0,track0->GetESDtrack()->GetID());//Bo: consistency 0 for neg
	pvertex->SetIndex(1,track1->GetESDtrack()->GetID());//Bo: consistency 1 for pos
	if (v0OK){
	  //	  AliV0vertex vertexjuri(*track0,*track1);
	  //	  vertexjuri.SetESDindexes(track0->fESDtrack->GetID(),track1->fESDtrack->GetID());
	  //	  event->AddV0(&vertexjuri);
	  pvertex->SetStatus(100);
	}
        pvertex->SetOnFlyStatus(kTRUE);
        pvertex->ChangeMassHypothesis(kK0Short);
	event->AddV0(pvertex);
      }
    }
  }
  //
  //
  // delete temporary arrays
  //  
  delete[] forbidden;
  delete[] minPointAngle;
  delete[] maxr;
  delete[] minr;
  delete[] norm;
  delete[] normdist;
  delete[] normdist1;
  delete[] normdist0;
  delete[] dist;
  delete[] itsmap;
  delete[] helixes;
  delete   pvertex;
}
//------------------------------------------------------------------------
void AliITSV0Finder::RefitV02(const AliESDEvent *event,
				AliITStrackerMI *tracker) 
{
  //
  //try to refit  V0s in the third path of the reconstruction
  //
  TTreeSRedirector &cstream = *(tracker->GetDebugStreamer());
  //
  Int_t  nv0s = event->GetNumberOfV0s();
  Float_t primvertex[3]={tracker->GetX(),tracker->GetY(),tracker->GetZ()};
  AliV0 v0temp;
  for (Int_t iv0 = 0; iv0<nv0s;iv0++){
    AliV0 * v0mi = (AliV0*)event->GetV0(iv0);
    if (!v0mi) continue;
    Int_t     itrack0   = v0mi->GetIndex(0);
    Int_t     itrack1   = v0mi->GetIndex(1);
    AliESDtrack *esd0   = event->GetTrack(itrack0);
    AliESDtrack *esd1   = event->GetTrack(itrack1);
    if (!esd0||!esd1) continue;
    AliITStrackMI tpc0(*esd0);  
    AliITStrackMI tpc1(*esd1);
    Double_t x,y,z; v0mi->GetXYZ(x,y,z); //I.B. 
    Double_t alpha =TMath::ATan2(y,x);   //I.B.
    if (v0mi->GetRr()>85){
      if (tpc0.Propagate(alpha,v0mi->GetRr())&&tpc1.Propagate(alpha,v0mi->GetRr())){
	v0temp.SetParamN(tpc0);
	v0temp.SetParamP(tpc1);
	v0temp.Update(primvertex);
	if (kFALSE) cstream<<"Refit"<<
	  "V0.="<<v0mi<<
	  "V0refit.="<<&v0temp<<
	  "Tr0.="<<&tpc0<<
	  "Tr1.="<<&tpc1<<
	  "\n";
	//Bo:	if (v0temp.GetDist2()<v0mi->GetDist2() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	if (v0temp.GetDcaV0Daughters()<v0mi->GetDcaV0Daughters() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	  v0mi->SetParamN(tpc0);
	  v0mi->SetParamP(tpc1);
	  v0mi->Update(primvertex);
	}
      }
      continue;
    }
    if (v0mi->GetRr()>35){
      AliITStrackerMI::CorrectForTPCtoITSDeadZoneMaterial(&tpc0);
      AliITStrackerMI::CorrectForTPCtoITSDeadZoneMaterial(&tpc1);
      if (tpc0.Propagate(alpha,v0mi->GetRr())&&tpc1.Propagate(alpha,v0mi->GetRr())){
	v0temp.SetParamN(tpc0);
	v0temp.SetParamP(tpc1);
	v0temp.Update(primvertex);
	if (kFALSE) cstream<<"Refit"<<
	  "V0.="<<v0mi<<
	  "V0refit.="<<&v0temp<<
	  "Tr0.="<<&tpc0<<
	  "Tr1.="<<&tpc1<<
	  "\n";
	//Bo:	if (v0temp.GetDist2()<v0mi->GetDist2() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	if (v0temp.GetDcaV0Daughters()<v0mi->GetDcaV0Daughters() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	  v0mi->SetParamN(tpc0);
	  v0mi->SetParamP(tpc1);
	  v0mi->Update(primvertex);
	}	
      }
      continue;
    }
    AliITStrackerMI::CorrectForTPCtoITSDeadZoneMaterial(&tpc0);
    AliITStrackerMI::CorrectForTPCtoITSDeadZoneMaterial(&tpc1);    
    //    if (tpc0.Propagate(alpha,v0mi->GetRr())&&tpc1.Propagate(alpha,v0mi->GetRr())){
    if (tracker->RefitAt(v0mi->GetRr(),&tpc0, v0mi->GetClusters(0)) && 
	tracker->RefitAt(v0mi->GetRr(),&tpc1, v0mi->GetClusters(1))){
      v0temp.SetParamN(tpc0);
      v0temp.SetParamP(tpc1);
      v0temp.Update(primvertex);
      if (kFALSE) cstream<<"Refit"<<
	"V0.="<<v0mi<<
	"V0refit.="<<&v0temp<<
	"Tr0.="<<&tpc0<<
	"Tr1.="<<&tpc1<<
	"\n";
      //Bo:      if (v0temp.GetDist2()<v0mi->GetDist2() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
      if (v0temp.GetDcaV0Daughters()<v0mi->GetDcaV0Daughters() || v0temp.GetV0CosineOfPointingAngle()>v0mi->GetV0CosineOfPointingAngle()){
	v0mi->SetParamN(tpc0);
	v0mi->SetParamP(tpc1);
	v0mi->Update(primvertex);
      }	
    }    
  }
}
//------------------------------------------------------------------------
