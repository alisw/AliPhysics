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
//  Comparison macro for TPC                                                 //
//  responsible: 
//  marian.ivanov@cern.ch                                                    //
//
//
/*
//


.L $ALICE_ROOT/STEER/AliGenInfo.C+
.L $ALICE_ROOT/STEER/AliTPCComparisonMI.C+


TPCCmpTr *t2 = new TPCCmpTr("genTracks.root","cmpTracks.root","galice.root",-1,1,0);
t2->Exec();


TCut cprim("cprim","MC.fVDist[3]<1");
TCut csec("cprim","MC.fVDist[3]>1");
TCut crec("crec","fReconstructed==1");
TCut cteta1("cteta1","abs(MC.fTrackRef.Theta()/3.1415-0.5)<0.25");
TCut cpos1("cpos1","abs(MC.fParticle.fVz/sqrt(MC.fParticle.fVx*MC.fParticle.fVx+MC.fParticle.fVy*MC.fParticle.fVy))<1");
TCut csens("csens","abs(sqrt(fVDist[0]**2+fVDist[1]**2)-170)<50");
TCut cmuon("cmuon","abs(MC.fParticle.fPdgCode==-13)");
AliTPCComparisonDraw comp; 
comp.SetIO();


//example
comp.DrawXY("fTPCinP0[3]","fTPCDelta[4]/fTPCinP1[3]","fReconstructed==1&&fPdg==-211"+cprim,"1",4,0.2,1.5,-0.06,0.06)
comp.fRes->Draw();
comp.fMean->Draw();  
comp.Eff("fTPCinP0[3]","fRowsWithDigits>120"+cteta1+cpos1,"fReconstructed>0",10,0,1.5)
comp.fRes->Draw();
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
#include "AliTPCtrack.h"
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
#endif
#include "AliGenInfo.h"
#include "AliTPCComparisonMI.h"

//
//
  


void AliTPCRecInfo::Update(AliMCInfo* info,AliTPCParam * param, Bool_t reconstructed, Int_t direction)
{
  //
  //
  //calculates derived variables
  //
  
  if (reconstructed==kFALSE){
    fReconstructed = kFALSE;
    return;
  }
  fReconstructed = kTRUE;
  // Find the nearest track reference
  //
  Double_t radius = TMath::Sqrt(fTrack.GetX()*fTrack.GetX()+fTrack.GetY()*fTrack.GetY());
  TClonesArray * tpcreferences = info->fTPCReferences;
  Int_t nentries = tpcreferences->GetEntriesFast();
  AliTrackReference * ref = 0;
  Double_t dr = 1000000;
  for (Int_t i=0;i<nentries;i++){
    AliTrackReference * refi = (AliTrackReference*)tpcreferences->UncheckedAt(i);
    if (!refi) continue;
    //if  ( (direction<0) && (refi->R()<radius-10)) continue;
    //if  ( (direction>0) && (refi->R()>radius+10)) continue;    
    if  (TMath::Abs(refi->R()-radius)<dr){
      dr = TMath::Abs(refi->R()-radius);
      ref = refi;
    }
  }  
  if (ref==0) {
    fReconstructed = kFALSE;
    return;
  }
  //
  //
  fdEdx = fTrack.GetdEdx();
  //
  // propagate track to nearest reference in direction
  TVector3 local = TPCCmpTr::TR2Local(ref,param);
  local.GetXYZ(fTRLocalCoord);
  //
  // rotate if neccessary
  Float_t ymax = local.X()*TMath::Tan(0.5*param->GetInnerAngle());  
  Float_t y = fTrack.GetYat(local.X());
  if (y > ymax) {
    fTrack.Rotate(param->GetInnerAngle());
  } else if (y <-ymax) {
    fTrack.Rotate(-param->GetInnerAngle());
  }   
  Double_t xtrack = fTrack.GetX();
  Double_t delta = (local.X()-xtrack)*0.1;
  for (Int_t i=0;i<9;i++){
    fTrack.PropagateTo(xtrack+delta*float(i));
  }
  fTrack.PropagateTo(local.X());

  Double_t par[5], cov[15], localX;
  fTrack.GetExternalParameters(localX,par);
  fTrack.GetExternalCovariance(cov);
  //
  //
  fRecPhi=TMath::ASin(par[2]) + fTrack.GetAlpha();
  Double_t phi2 =TMath::ATan2(ref->Py(),ref->Px());
  if (phi2<0) phi2+=2*TMath::Pi();
  fGenPhi =phi2;
  //
  if (fRecPhi<0) fRecPhi+=2*TMath::Pi();
  if (fRecPhi>=2*TMath::Pi()) fRecPhi-=2*TMath::Pi();
  fLambda = TMath::ATan(par[3]);
  fRecPt_1 = TMath::Abs(par[4]);  
  //TPC
  //
  //
  Double_t phi=TMath::ATan2(par[0],localX) + fTrack.GetAlpha();
  Double_t r=TMath::Sqrt(localX*localX + par[0]*par[0]);
  fTPCinR1[0]= r*TMath::Cos(phi); 
  fTPCinR1[1]= r*TMath::Sin(phi); 
  fTPCinR1[2]= par[1];
  fTPCinR1[3] = TMath::Sqrt(fTPCinR1[0]*fTPCinR1[0]+fTPCinR1[1]*fTPCinR1[1]);
  fTPCinR1[4] = TMath::ATan2(fTPCinR1[1],fTPCinR1[0]);	
  //
  fTPCinP1[0] = fTrack.Px();
  fTPCinP1[1] = fTrack.Py();
  fTPCinP1[2] = fTrack.Pz();
  fTPCinP1[3] = TMath::Sqrt(fTPCinP1[0]*fTPCinP1[0]+fTPCinP1[1]*fTPCinP1[1]);
  //
  fTPCinR0[0] = ref->X();	
  fTPCinR0[1] = ref->Y();	
  fTPCinR0[2] = ref->Z();
  fTPCinR0[3] = TMath::Sqrt(fTPCinR0[0]*fTPCinR0[0]+fTPCinR0[1]*fTPCinR0[1]);
  fTPCinR0[4] = TMath::ATan2(fTPCinR0[1],fTPCinR0[0]);
  //
  fTPCinP0[0]=ref->Px();
  fTPCinP0[1]=ref->Py();
  fTPCinP0[2]=ref->Pz();
  fTPCinP0[3]=ref->Pt();
  //
  //
  if (fTPCinP1[3]>0.0000001){
    fTPCAngle1[0] = TMath::ATan2(fTPCinP1[1],fTPCinP1[0]);
    fTPCAngle1[1] = TMath::ATan(fTPCinP1[2]/fTPCinP1[3]);  
    //
    fTPCAngle0[0] = TMath::ATan2(fTPCinP0[1],fTPCinP0[0]);
    fTPCAngle0[1] = TMath::ATan(fTPCinP0[2]/fTPCinP0[3]);
  }    
  //
  fTPCDelta[0] = (fTPCinR0[4]-fTPCinR1[4])*fTPCinR1[3];  //delta rfi
  fTPCPools[0] = fTPCDelta[0]/TMath::Sqrt(cov[0]);
  fTPCDelta[1] = (fTPCinR0[2]-fTPCinR1[2]);              //delta z
  fTPCPools[1] = fTPCDelta[1]/TMath::Sqrt(cov[2]);
  //
  fTPCDelta[2] = (fTPCAngle0[0]-fTPCAngle1[0]);
  fTPCPools[2] = fTPCDelta[2]/TMath::Sqrt(cov[5]);
  //
  fTPCDelta[3] = TMath::Tan(fTPCAngle0[1])-TMath::Tan(fTPCAngle1[1]);
  fTPCPools[3] = fTPCDelta[3]/TMath::Sqrt(cov[9]);
  fTPCDelta[4] = (fTPCinP0[3]-fTPCinP1[3]);
  Double_t sign = (fTrack.GetC()>0) ? 1:-1;
  fTPCPools[4] = sign*(1./fTPCinP0[3]-1./fTPCinP1[3])/TMath::Sqrt(cov[14]);      
  
}




////////////////////////////////////////////////////////////////////////
TPCCmpTr::TPCCmpTr()
{
  Reset();
}

////////////////////////////////////////////////////////////////////////
TPCCmpTr::TPCCmpTr(const char* fnGenTracks,
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
TPCCmpTr::~TPCCmpTr()
{
  //if (fLoader) {
  //  delete fLoader;
  //}
}

//////////////////////////////////////////////////////////////
Int_t TPCCmpTr::SetIO()
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

Int_t TPCCmpTr::SetIO(Int_t eventNr)
{
  //
  // 
  // SET INPUT
  //gAlice->GetEvent(eventNr);
  //  fLoader->SetEventNumber(eventNr);  
  //fLoader = AliRunLoader::Open("galice.root");
  fLoader->GetEvent(eventNr);  
  //
  AliTPCLoader * tpcl = (AliTPCLoader*)fLoader->GetLoader("TPCLoader");
  tpcl->LoadTracks();
  fTreeRecTracks = tpcl->TreeT();
  
  //
  return 0;
}



////////////////////////////////////////////////////////////////////////
void TPCCmpTr::Reset()
{
  fEventNr = 0;
  fNEvents = 0;
  fTreeCmp = 0;
  //
  fFileGenTracks = 0;
  fDebug = 0;
  //
  fParamTPC = 0;
  fTreeRecTracks = 0;
  fTreePoints =0;

  fTPCTrack = 0; 
  fTracks   = 0;
  fTrackPoints =0;
}

////////////////////////////////////////////////////////////////////////
Int_t TPCCmpTr::Exec(Int_t nEvents, Int_t firstEventNr)
{
  fNEvents = nEvents;
  fFirstEventNr = firstEventNr;
  return Exec();
}

////////////////////////////////////////////////////////////////////////
Int_t TPCCmpTr::Exec()
{
  TStopwatch timer;
  timer.Start();

  if (SetIO()==1) 
    return 1;
   
  cerr<<"fFirstEventNr, fNEvents: "<<fFirstEventNr<<" "<<fNEvents<<endl;
  for (Int_t eventNr = fFirstEventNr; eventNr < fFirstEventNr+fNEvents;
       eventNr++) {
    fEventNr=eventNr;
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
    if (TreeTLoop(eventNr)>0) return 1;

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
Bool_t TPCCmpTr::ConnectGenTree()
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


  if (fDebug > 1) {
    cout<<"Number of gen. tracks with TR: "<<fTreeGenTracks->GetEntries()<<endl;
  }
  return kTRUE;
}


////////////////////////////////////////////////////////////////////////
void TPCCmpTr::CreateTreeCmp() 
{
  fFileCmp = TFile::Open(fFnCmp,"RECREATE");
  if (!fFileCmp) {
    cerr<<"Error in CreateTreeCmp: cannot open file "<<fFnCmp<<endl;
    return;
  }


  fTreeCmp    = new TTree("TPCcmpTracks","TPCcmpTracks");
  //
  fMCInfo = new AliMCInfo;
  fRecInfo = new AliTPCRecInfo;
  //
  fTPCTrack = new AliTPCtrack;
   //
  fTreeCmp->Branch("MC","AliMCInfo",&fMCInfo);
  fTreeCmp->Branch("RC","AliTPCRecInfo",&fRecInfo);
  fTreeCmp->Branch("fTPCTrack","AliTPCtrack",&fTPCTrack);
  //
  fTreeCmp->AutoSave(); 

}
////////////////////////////////////////////////////////////////////////
void TPCCmpTr::CloseOutputFile()  
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

TVector3 TPCCmpTr::TR2Local(AliTrackReference *trackRef,
			    AliTPCParam *paramTPC) {

  Float_t x[3] = { trackRef->X(),trackRef->Y(),trackRef->Z()};
  Int_t index[4];
  paramTPC->Transform0to1(x,index);
  paramTPC->Transform1to2(x,index);
  return TVector3(x);
}
////////////////////////////////////////////////////////////////////////

Int_t TPCCmpTr::TreeTLoop(Int_t eventNr)
{
  //
  // loop over all TPC reconstructed tracks and store info in memory
  //
  TStopwatch  timer;
  timer.Start();
  
  if (!fTreeRecTracks) {
    cerr<<"Can't get a tree with TPC rec. tracks  "<<endl;
    return 1;
  }
  //fTreePoints=(TTree*)fFileRecTracks->Get("trackDebug");
  
  Int_t nEntries = (Int_t) fTreeRecTracks->GetEntries();
  if (fDebug > 2) cout<<"Event, rec. tracks: "<<eventNr<<" "
		      <<nEntries<<endl;
  TBranch * br= fTreeRecTracks->GetBranch("tracks");
  br->SetAddress(&fTPCTrack);
  TBranch *brp = 0;
  if (fTreePoints) brp = fTreePoints->GetBranch("debug");

  if (fTracks){
    fTracks->Delete();    
    delete fTracks;
  }
  if (fTrackPoints){
    fTrackPoints->Delete();
    delete fTrackPoints;
    fTrackPoints =0;
  }
  fTracks      = new TObjArray(nEntries);
  if (brp){
    fTrackPoints = new TObjArray(nEntries);
  }
  else fTrackPoints = 0;

  //
  //load tracks to the memory
  for (Int_t i=0; i<nEntries;i++){
    AliTPCtrack * track = new AliTPCtrack;
    br->SetAddress(&track);
    br->GetEntry(i);
    fTracks->AddAt(track,i);
  }
  //
  //load track points to the memory
  if (brp) for (Int_t i=0; i<nEntries;i++){
    TClonesArray * arr = new TClonesArray("AliTPCTrackPoint2");
    brp->SetAddress(&arr);
    brp->GetEntry(i);
    if (arr!=0)
      for (Int_t j=0;j<arr->GetEntriesFast();j++){
	AliTPCTrackPoint2 * point = (AliTPCTrackPoint2*)arr->UncheckedAt(j);
	if (point && point->fID>=0){
	  fTrackPoints->AddAt(arr,point->fID);
	  break;
	}
      }    
  }
  //

  //
  for (Int_t iEntry=0; iEntry<nEntries;iEntry++){
    //br->GetEntry(iEntry);
    fTPCTrack = (AliTPCtrack*)fTracks->UncheckedAt(iEntry);
    //
    Int_t label = fTPCTrack->GetLabel();
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
Int_t TPCCmpTr::TreeGenLoop(Int_t eventNr)
{
//
// loop over all entries for a given event, find corresponding 
// rec. track and store in the fTreeCmp
//
  TStopwatch timer;
  timer.Start();
  Int_t entry = 0;
  Double_t nParticlesTR = fTreeGenTracks->GetEntriesFast();
  TBranch * branch = fTreeCmp->GetBranch("RC");
  branch->SetAddress(&fRecInfo); // set all pointers

  while (entry < nParticlesTR) {
    fTreeGenTracks->GetEntry(entry);
    entry++;
    if (eventNr < fMCInfo->fEventNr) continue;
    if (eventNr > fMCInfo->fEventNr) continue;
    //
    if (fDebug > 2 && fMCInfo->fLabel < 10) {
      cerr<<"Fill track with a label "<<fMCInfo->fLabel<<endl;
    }
    fRecInfo->Reset();
    //
    if (fIndexRecTracks[fMCInfo->fLabel*4] >= 0) {
      fTPCTrack= (AliTPCtrack*)fTracks->UncheckedAt(fIndexRecTracks[fMCInfo->fLabel*4]);
      
      //      if (nBytes > 0) {
      if (fTPCTrack) {
	//
	//
	//if multifound find track with biggest pt - returning tracks - loos energy
	//
	if (fIndexRecTracks[fMCInfo->fLabel*4]+1){
	  Double_t pt = TMath::Abs(1/fTPCTrack->Get1Pt());
	  for (Int_t i=1;i<4;i++){
	    if (fIndexRecTracks[fMCInfo->fLabel*4+i]>=0){
	      AliTPCtrack * track = (AliTPCtrack*)fTracks->UncheckedAt(fIndexRecTracks[fMCInfo->fLabel*4+i]);
	      Double_t pt2 = TMath::Abs(1/track->Get1Pt());	      
	      if  (pt2>pt)
		fTPCTrack = track;		   
	    }
	  }
	}
	//
	fRecInfo->fTP=0;
	fRecInfo->fTrack =*fTPCTrack; 
	fRecInfo->fReconstructed = 1;
	fRecInfo->fFake = fFakeRecTracks[fMCInfo->fLabel];
	fRecInfo->fMultiple = fMultiRecTracks[fMCInfo->fLabel];
	fRecInfo->Update(fMCInfo, fParamTPC, kTRUE, 1);
      }      
    } 
    else{
      fRecInfo->Update(fMCInfo, fParamTPC, kFALSE, 1);      
    }    
    fTreeCmp->Fill();
  }
  fTreeCmp->AutoSave();
  fTracks->Delete();
  fTPCTrack =0;
  if (fTrackPoints){
    fTrackPoints->Delete();
    delete fTrackPoints;
    fTrackPoints =0;
  } 
  printf("Time spended in TreeGenLoop\n");
  timer.Print();
  if (fDebug > 2) cerr<<"end of TreeGenLoop"<<endl;  
  return 0;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void AliTPCComparisonDraw::SetIO(const char *fname)
{
  TFile *file = TFile::Open(fname);
  if (!file){
    printf("Comparison file couldn't be generated\n"); 
    return;
   }
  //
  fTree = (TTree*) file->Get("TPCcmpTracks");
  if (!fTree) {
    printf("no track comparison tree found\n");
    file->Close();
    delete file;
  }

}

void AliTPCComparisonDraw::Eff()
{
  
}

void AliTPCComparisonDraw::ResPt()
{
}
