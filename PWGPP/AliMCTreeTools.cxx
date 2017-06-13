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

/*
  Class to enable easy correlation of the reconstructed information and MC(true)    
  WARNING : this code is for fast visualization and prototyping NOT FOR PRODUCTION SOFTWARE
          : Goals was flexibility not speed
  Example usage:
    1.) GetParticle properties at vertex (using ESd tree) 
        - e.g to compare reconstructed track momenta and MC momenta  ( mentioned in PWGPP-205)     
        esdTree->Draw("Tracks[].fIp.Pt()/AliMCTreeTools::GetValueAt(Entry$,abs(Tracks[].fLabel),0,0,0,8,1,0)","(Tracks[].fTPCncls)>80&&Tracks[].fITSncls>3","")	
    2.) Find TPC space point porperty using closest MC information  ( mentioned in ATO-168) 
        - used for loop finder as a ground cluster/tracklet true
	-- fast as used  togeter with MakeCacheTree(TTree * tree, TString varList, TString outFile, TString outTree, TCut selection);          
   author marian.ivanov@cern.ch
*/


/*
  .L $AliPhysics_SRC/PWGPP/AliMCTreeTools.cxx+
  AliMCTreeTools::SetWDir(gSystem->pwd());
  TFile * f = TFile::Open("AliESDs.root");
  esdTree->Draw("Tracks[].fIp.Pt()/AliMCTreeTools::GetValueAt(Entry$,abs(Tracks[].fLabel),0,0,0,8,1,0)","(Tracks[].fTPCncls)>80&&Tracks[].fITSncls>3","")

*/


#include "TStopwatch.h"
#include "TTree.h" 
#include "TChain.h"
#include "TVectorF.h"
#include "AliStack.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TTreeStream.h"
#include "AliRunLoader.h"
#include "AliTrackReference.h"
#include "AliExternalTrackParam.h"
#include "AliHelix.h"
#include "TCut.h"
#include "AliTreePlayer.h"
#include "THn.h"
#include "TF3.h"
#include "TStatToolkit.h"
#include <stdarg.h>
#include "AliNDLocalRegression.h"
#include "AliMCTreeTools.h"

ClassImp(AliMCTreeTools)


const Double_t kMaxRadius=300; 
TClonesArray* trackRefs = 0;
map<Int_t, TClonesArray*> mapTR;
map<Int_t, AliTrackReference*> mapFirstTr;
map<Int_t, TClonesArray*> mapHelix;
//

TString wdir="";
AliStack *stack=NULL;
AliRunLoader *rl=NULL;
TTree *treeTR=0;
TTree * treeCl=0;
TTree * treeMC=0;
Double_t bz=-5;


void AliMCTreeTools::SetWDir(TString dir){
  wdir=dir;
}

//
void AliMCTreeTools::ClearCache(){
  // clear cache needed to speedup acces to MC information
  mapTR.clear();
  mapHelix.clear(); 
  mapFirstTr.clear();
}



void  AliMCTreeTools::InitStack(Int_t iEvent){
  if (rl==NULL){
    rl = AliRunLoader::Open(Form("%s/galice.root", wdir.Data())); 
  }
  if (treeTR && rl->GetEventNumber()==iEvent) return; 
  rl->GetEvent(iEvent);
  rl->LoadTrackRefs();
  rl->LoadKinematics();
  //
  treeTR=rl->TreeTR();
  stack=rl->Stack();
  AliMCTreeTools::ClearCache();
}


Double_t  AliMCTreeTools::GetValueAt(Int_t iEvent, Int_t itrack, Double_t x, Double_t y, Double_t z, Int_t returnValue, Int_t interpolationType, Int_t verbose){
  //
  // GetValueAt
  // return values 
  //    0 - gx
  //    1 - gy
  //    2 - gz
  //    3 - r
  //    4 - phi
  //    5 - deltaPhi (pos-dir)
  //    6 - minimal distance to track ref.
  //    7 - weighted interpolated  P
  //    8 - weighted interpolated  Pt
  //    9  - P at first tpc ref
  //    10 - dEdx at local P
  //    11 - pdg code
  //    12 - P at vertex
  //    13 - P at vertex

  Int_t index= FindNearestReference(iEvent,itrack,x,y,z,0,verbose);
  if (returnValue==11 || returnValue==12 ||  returnValue==13){
    if (stack==NULL) return -1;
    TParticle *particle = stack->Particle(itrack);
    Int_t pdgCode=(particle!=NULL) ? particle->GetPdgCode():-1;
    if (returnValue==11) return pdgCode;
    if (returnValue==12) return (particle!=NULL) ? particle->P():-1;
    if (returnValue==13) return (particle!=NULL) ? particle->Pt():-1;
  }
  TClonesArray *helixArray = mapHelix[itrack];
  TClonesArray *trArray = mapTR[itrack];
  if (helixArray==NULL) return 0;
  if (trArray==NULL) return 0;
  if (helixArray->GetEntriesFast()<2) return 0;
  AliHelix *helix0=(AliHelix*)helixArray->At(index);
  AliHelix *helix1=(AliHelix*)helixArray->At(index+1);
  AliTrackReference *ref0= (AliTrackReference *)trArray->At(index);
  AliTrackReference *ref1= (AliTrackReference *)trArray->At(index+1);
  if (ref0->GetTrack()!=itrack){
    static Int_t wCounter=0;
    wCounter++;
    return -1;
  }
  Double_t d0=TMath::Sqrt((ref0->X()-x)*(ref0->X()-x)+(ref0->Y()-y)*(ref0->Y()-y)+(ref0->Z()-z)*(ref0->Z()-z));
  Double_t d1=TMath::Sqrt((ref1->X()-x)*(ref1->X()-x)+(ref1->Y()-y)*(ref1->Y()-y)+(ref1->Z()-z)*(ref1->Z()-z));
  Double_t w0=d1/(d1+d0);
  Double_t w1=d0/(d1+d0);
  //
  if (helix0==NULL || helix1==NULL) return 0;
  if (returnValue==6) {
    Double_t value=TMath::Sqrt(TMath::Min(d0,d1));
    return value;
  }
  if (returnValue==7) {
    Double_t value =w0*ref0->P()+w1*ref1->P();
    return value;
  }
  if (returnValue==8) {
    Double_t value =w0*ref0->Pt()+w1*ref1->Pt();
    return value;
  }
  //
  if (returnValue==9){
    AliTrackReference * ref =  mapFirstTr[itrack];
    return (ref!=NULL) ? ref->P():0;
  }
   if (returnValue==10){
     TParticle *particle = stack->Particle(itrack);
     if (particle==NULL) return 0;
     TParticlePDG *mcParticle = TDatabasePDG::Instance()->GetParticle(particle->GetPdgCode());
     if (mcParticle==NULL) return 0;
     Double_t mass = mcParticle->Mass();
     Double_t wP =w0*ref0->P()+w1*ref1->P();
     Double_t dEdx= AliExternalTrackParam::BetheBlochAleph(wP/mass);
     return dEdx;
  }
  //

  Double_t xyz[3],  dxyz[3], ddxyz[3]; 
  if (interpolationType==0){
    helix0->Evaluate(helix0->GetPhase(x,y),xyz,dxyz,ddxyz);
    Double_t zLoop=TMath::TwoPi()/helix0->GetHelix(4);
    if (TMath::Abs(zLoop)>0.5){
      Int_t nLoops=TMath::Nint((xyz[2]-z)/zLoop);
      xyz[2]-=nLoops*zLoop;
    }
  }
  if (interpolationType==1){
    helix1->Evaluate(helix1->GetPhase(x,y),xyz,dxyz,ddxyz);
    Double_t zLoop=TMath::TwoPi()/helix1->GetHelix(4);
    if (TMath::Abs(zLoop)>0.5){
      Int_t nLoops=TMath::Nint((xyz[2]-z)/zLoop);
      xyz[2]-=nLoops*zLoop;
    }
  }
  if (interpolationType==2){
    helix0->Evaluate(helix0->GetPhase(x,y),xyz,dxyz,ddxyz);
    Double_t zLoop=TMath::TwoPi()/helix0->GetHelix(4); 
    if (TMath::Abs(zLoop)>0.5){
      Int_t nLoops=TMath::Nint((xyz[2]-z)/zLoop);
      xyz[2]-=nLoops*zLoop;
    }
    Double_t xyz1[3],  dxyz1[3], ddxyz1[3]; 
    helix1->Evaluate(helix1->GetPhase(x,y),xyz1,dxyz1,ddxyz1);
    zLoop=TMath::TwoPi()/helix1->GetHelix(4); 
    if (TMath::Abs(zLoop)>0.5){
      Int_t nLoops=TMath::Nint((xyz1[2]-z)/zLoop);
      xyz1[2]-=nLoops*zLoop;
    }
    for (Int_t i=0; i<3;i++){
      xyz[i]*=w0;
      xyz[i]+=w1*xyz1[i];
    }
  }
  if (returnValue<3) return xyz[returnValue];
  if (returnValue==3){
    return TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  }
  if (returnValue==4){
    return TMath::ATan2(xyz[1],xyz[0]);
  }
  if (returnValue==5){
    return TMath::ATan2(xyz[1],xyz[0])-TMath::ATan2(dxyz[1],dxyz[0]);
  }

}


Double_t AliMCTreeTools::FindNearestReference(Int_t iEvent, Int_t itrack, Double_t x, Double_t y, Double_t z, Int_t returnValue, Int_t verbose){
  //  
  TDatabasePDG *pdg =  TDatabasePDG::Instance();
  if (itrack<0) return 0;
  InitStack(iEvent); 
  TClonesArray *trefs=mapTR[itrack];
  if (trefs==NULL){ // cache Trackrefs if not done before
    treeTR->SetBranchAddress("TrackReferences", &trackRefs); 
    treeTR->GetEntry(stack->TreeKEntry(itrack));
    mapTR[itrack]=(TClonesArray*)trackRefs->Clone();
    trefs= mapTR[itrack];
    TClonesArray * helixArray = new TClonesArray("AliHelix",trackRefs->GetEntries());
    helixArray->ExpandCreateFast(trackRefs->GetEntries());
    TParticle *particle = stack->Particle(itrack);
    TParticlePDG *mcparticle = pdg->GetParticle(particle->GetPdgCode());
    if (mcparticle==NULL) return 0;
    //
    Float_t conversion = -1000/0.299792458/bz; // AliTracker::GetBz();
    if (mcparticle->Charge()!=0){
      for (Int_t itr=0; itr<trackRefs->GetEntriesFast(); itr++){
	AliTrackReference *ref=  (AliTrackReference *)trackRefs->At(itr);
	if (ref->GetTrack()!=itrack) continue; 
	if ( ref->DetectorId()==AliTrackReference::kTPC  && mapFirstTr[itrack]==NULL){
	  mapFirstTr[itrack]=(AliTrackReference*)ref->Clone();
	}
	Double_t xyz[3]={ref->X(),ref->Y(),ref->Z()};
	Double_t pxyz[3]={ref->Px(),ref->Py(),ref->Pz()};
	if (ref->P()==0) pxyz[0]+=0.00000000001; // create dummy track reference in case of 0 moment (track disappeared)
	new ((*helixArray)[itr]) AliHelix(xyz,pxyz,mcparticle->Charge()/3.,conversion);
      }
      mapHelix[itrack]=helixArray;
    }
  }

  Int_t nTrackRefs = trefs->GetEntriesFast();
  Int_t nTPCRef=0;
  AliTrackReference *refNearest=0;
  TVectorF fdist(nTrackRefs);
  for (Int_t itrR = 0; itrR < nTrackRefs; ++itrR) {
    AliTrackReference* ref = static_cast<AliTrackReference*>(trefs->UncheckedAt(itrR));	
    Double_t lDist=(ref->X()-x)*(ref->X()-x)+(ref->Y()-y)*(ref->Y()-y)+(ref->Z()-z)*(ref->Z()-z);
    fdist[itrR]=lDist;
  }
  Double_t dist=250*250;
  Int_t index0=0,index1=0;
  for (Int_t itrR = 1; itrR < nTrackRefs-1; ++itrR){
    if (fdist[itrR]<dist){      
      dist=fdist[itrR];
      if (fdist[itrR-1]<fdist[itrR+1]){
	index0=itrR-1;
	index1=itrR;	
      }else{
	index0=itrR;
	index1=itrR+1;	
      }
    }
    refNearest=static_cast<AliTrackReference*>(trefs->UncheckedAt(index0));
    if (verbose ) {
      trefs->UncheckedAt(index0)->Print();
      trefs->UncheckedAt(index1)->Print();
    }
  }
  if (returnValue==0) return  index0;
  if (returnValue==1) return  TMath::Sqrt(dist);
}

