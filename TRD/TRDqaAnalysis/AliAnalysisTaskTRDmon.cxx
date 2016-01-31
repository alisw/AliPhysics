//**************************************************************************
//* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//*                                                                        *
//* Author: The ALICE Off-line Project.                                    *
//* Contributors are mentioned in the code where appropriate.              *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************/
//////////////////////////////////////////////////
// Analysis task used for TRD monitoring        //
// Input:   TChain of ESDs                      //
// Output:  TObjArray of histograms             //
//////////////////////////////////////////////////

#include <TObjArray.h>
#include <TObject.h>
#include <TParticle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliTRDcluster.h"
#include "AliTrackPointArray.h"

#include <cstdio>

#ifndef ALIANALYSISTASKTRDMON_H
#include "AliAnalysisTaskTRDmon.h"
#endif

ClassImp(AliAnalysisTaskTRDmon)

AliAnalysisTaskTRDmon::AliAnalysisTaskTRDmon(const Char_t *name):
  AliAnalysisTask(name, ""),
  fESD(0x0),
  fESDfriend(0x0),
  fEventTriggerName("CINT1B-ABCE-NOPF-ALL"),
  fIsCollisionEvent(kTRUE),
  fOutStorage(0x0),
  fHzvert1(0x0),
  fHzvert2(0x0),
  fHntracks(0x0),
  fHntracks2(0x0),
  fHntracks3(0x0),
  fHdca(0x0),
  fHdcaz(0x0),
  fHpt(0x0),
  fHpt2(0x0),
  fHpt3(0x0),
  fHpt3n(0x0),
  fHpt4(0x0),
  fHpt4n(0x0),
  fHtheta(0x0),
  fHphi(0x0),
  fHtpccl(0x0),
  fHtpccl2(0x0),
  fHdedxp(0x0),
  fHetaphi(0x0),
  fHetancl(0x0),
  fHphincl(0x0),
  fHtrdtr(0x0),
  fHtrdtr2(0x0),
  fHph(0x0),
  fHph2d(0x0),
  fHncltrkl(0x0),
  fHntrkl(0x0),
  fHntrklVsP(0x0),
  fHsm(0x0),
  fHetantr(0x0),
  fHphintr(0x0),
  fHcltime(0x0),
  fHcldiff(0x0),
  fHxyA(0x0),
  fHxyC(0x0),
  fHPropagXY(0x0),
  fHPropagRZ(0x0),
  fHnFriendTracks(0x0),
  fHnCalibObjects(0x0),
  fHTpcEtaPhiLocaln(0x0),
  fHTrdEtaPhiLocaln(0x0),
  fHTpcEtaPhiLocalp(0x0),
  fHTrdEtaPhiLocalp(0x0),
  fHTpcEtaPhiSagitan(0x0),
  fHTrdEtaPhiSagitan(0x0),
  fHTpcEtaPhiSagitap(0x0),
  fHTrdEtaPhiSagitap(0x0),
  fHTpcEtaPhiDeltan(0x0),
  fHTrdEtaPhiDeltan(0x0),
  fHTpcEtaPhiDeltap(0x0),
  fHTrdEtaPhiDeltap(0x0),
  fHTpcEtaPhin(0x0),
  fHTrdEtaPhin(0x0),
  fHTpcEtaPhip(0x0),
  fHTrdEtaPhip(0x0),
  fHTpcRef3Dpos(0x0),
  fHTpcRef3Dneg(0x0),
  fHTrdRef3Dpos(0x0),
  fHTrdRef3Dneg(0x0),
  fHTrdEtaPhiLocalNTracklets(0x0),
  fHTrdEtaPhiSagitaNTracklets(0x0),
  fHTrdEtaPhiDeltaNTracklets(0x0),
  fHTrdEtaPhiNTracklets(0x0)
{
  //
  // constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(0, TObjArray::Class());
}

AliAnalysisTaskTRDmon::AliAnalysisTaskTRDmon(const AliAnalysisTaskTRDmon &task) :
  AliAnalysisTask(task),
  fESD(0x0),
  fESDfriend(0x0),
  fEventTriggerName(task.GetTriggerName()),
  fIsCollisionEvent(task.GetIsCollisionEvent()),
  fOutStorage(0x0),
  fHzvert1(0x0),
  fHzvert2(0x0),
  fHntracks(0x0),
  fHntracks2(0x0),
  fHntracks3(0x0),
  fHdca(0x0),
  fHdcaz(0x0),
  fHpt(0x0),
  fHpt2(0x0),
  fHpt3(0x0),
  fHpt3n(0x0),
  fHpt4(0x0),
  fHpt4n(0x0),
  fHtheta(0x0),
  fHphi(0x0),
  fHtpccl(0x0),
  fHtpccl2(0x0),
  fHdedxp(0x0),
  fHetaphi(0x0),
  fHetancl(0x0),
  fHphincl(0x0),
  fHtrdtr(0x0),
  fHtrdtr2(0x0),
  fHph(0x0),
  fHph2d(0x0),
  fHncltrkl(0x0),
  fHntrkl(0x0),
  fHntrklVsP(0x0),
  fHsm(0x0),
  fHetantr(0x0),
  fHphintr(0x0),
  fHcltime(0x0),
  fHcldiff(0x0),
  fHxyA(0x0),
  fHxyC(0x0),
  fHPropagXY(0x0),
  fHPropagRZ(0x0),
  fHnFriendTracks(0x0),
  fHnCalibObjects(0x0),
  fHTpcEtaPhiLocaln(0x0),
  fHTrdEtaPhiLocaln(0x0),
  fHTpcEtaPhiLocalp(0x0),
  fHTrdEtaPhiLocalp(0x0),
  fHTpcEtaPhiSagitan(0x0),
  fHTrdEtaPhiSagitan(0x0),
  fHTpcEtaPhiSagitap(0x0),
  fHTrdEtaPhiSagitap(0x0),
  fHTpcEtaPhiDeltan(0x0),
  fHTrdEtaPhiDeltan(0x0),
  fHTpcEtaPhiDeltap(0x0),
  fHTrdEtaPhiDeltap(0x0),
  fHTpcEtaPhin(0x0),
  fHTrdEtaPhin(0x0),
  fHTpcEtaPhip(0x0),
  fHTrdEtaPhip(0x0),
  fHTpcRef3Dpos(0x0),
  fHTpcRef3Dneg(0x0),
  fHTrdRef3Dpos(0x0),
  fHTrdRef3Dneg(0x0),
  fHTrdEtaPhiLocalNTracklets(0x0),
  fHTrdEtaPhiSagitaNTracklets(0x0),
  fHTrdEtaPhiDeltaNTracklets(0x0),
  fHTrdEtaPhiNTracklets(0x0)
{
  //
  // copy constructor
  //
}

AliAnalysisTaskTRDmon& AliAnalysisTaskTRDmon::operator=(const AliAnalysisTaskTRDmon &task) 
{
  //
  // Assignment operator
  //
  if(&task == this) return *this;
  AliAnalysisTask::operator=(task);
  fEventTriggerName = task.GetTriggerName();
  fIsCollisionEvent = task.GetIsCollisionEvent();

  return *this;
}

void AliAnalysisTaskTRDmon::ConnectInputData(Option_t *){
  //
  // Connect input data. Overloads AliAnalysisTask::ConnectInputData()
  //
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0)); 
  if (!tree) Printf("ERROR: Could not read chain from input slot 0"); 
  else tree->SetBranchStatus("ESDfriend*", kTRUE); 
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!esdH){
    AliError("ERROR: Could not get ESDInputHandler");
  } else {
    fESD = (AliESDEvent*)esdH->GetEvent();
    fESDfriend = dynamic_cast<AliESDfriend*> (fESD->FindListObject("AliESDfriend"));
    
  }
  
}

void AliAnalysisTaskTRDmon::CreateOutputObjects(){
  //
  // Creates the output object list, initializes histograms, adds histograms to the output list
  //

  fOutStorage = new TObjArray(40);
  fHzvert1 = new TH1F("fHzvert1", "Zvertex", 100, -25., 25.);
  fHzvert2 = new TH1F("fHzvert2", "Zvertex TPC", 100, -25., 25.);
  fHdca = new TH1F("fHdca", "DCA", 100, -100., 100.);
  fHdcaz = new TH1F("fHdcaz", "DCA", 100, -100., 100.);
  fHntracks = new TH1F("fHntracks", "Ntracks", 100, 0., 100.);
  fHntracks2 = new TH1F("fHntracks2", "Ntracks dca.lt.1", 100, 0., 100.);
  fHntracks3 = new TH1F("fHntracks3", "Ntracks dca.lt.1, NclTPC", 100, 0., 100.);
  Float_t binLimits[33] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			   1.0, 1.1, 1.2, 1.3, 1.4, 
			   1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
			   3.4, 3.8, 4.2, 4.6, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  fHpt = new TH1F("fHpt", "Pt", 32, binLimits);
  fHpt2 = new TH1F("fHpt2", "Pt dca.lt.1", 32, binLimits);
  fHpt2->Sumw2();
  fHpt3 = new TH1F("fHpt3", "Pt dca.lt.1, NclTPC +", 32, binLimits);
  fHpt3->Sumw2();
  fHpt3n = new TH1F("fHpt3n", "Pt dca.lt.1, NclTPC -", 32, binLimits);
  fHpt3n->Sumw2();
  fHpt4 = new TH1F("fHpt4", "Pt dca.lt.1, NclTPC & TRD +", 32, binLimits);
  fHpt4->Sumw2();
  fHpt4n = new TH1F("fHpt4n", "Pt dca.lt.1, NclTPC & TRD -", 32, binLimits);
  fHpt4n->Sumw2();
  fHtheta = new TH1F("fHtheta", "theta", 220,.5,2.7);
  fHphi = new TH1F("fHphi", "phi", 157,0,6.28);
  fHtpccl = new TH1F("fHtpccl", "nClTPC", 160, 0., 160.);
  fHtpccl2 = new TH1F("fHtpccl2", "nClTPC p.gt.1", 160, 0., 160.);
  
  //Float_t binx[21]={.1,.13,.16,.22,.3,.4,.55,.7,1,1.3,1.6,2.2,3,4,5.5,7,10,13,16,22,30};
  
  fHdedxp = new TH2F("fHdedxp", "dE/dx vs. p TPC", 100, 0.1,10.1, 120, 0,600.);
  // fHdedxp->SetBit(TH1::kCanRebin); 
  fHetaphi = new TH2F("fHetaphi", "eta-phi TPC", 50, -1, 1, 157, 0, 6.28);
  fHetancl = new TH2F("fHetancl", "eta-Ncl TPC",50, -1, 1, 160, 0, 160.);
  fHphincl = new TH2F("fHphincl", "phi-Ncl TPC",157, 0, 6.28, 160, 0, 160.);
  fHtrdtr = new TH1F("fHtrdtr", "nTRDlayers", 7, -.5, 6.5);
  fHtrdtr2 = new TH1F("fHtrdtr2", "nTRDlayers p.gt.1", 7, -.5, 6.5);
  fHph = new TProfile("fHph", "ph", 30,0,30);
  fHph2d = new TH2F("fHph2d", "PH2d TRD", 30, 0, 30, 200, 0, 1000.);
  fHncltrkl = new TH1F("fHncltrkl", "Nclusters/tracklet", 50,0,50.);
  fHntrkl = new TH1F("fHntrkl", "Ntracklets", 7,-.5,6.5);
  Double_t binsP[19] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0,
			2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0};
  fHntrklVsP = new TH2F("fHntrklVsP", "Ntracklets vs momentum", 18, binsP, 7, -0.5, 6.5);
  fHsm = new TH1F("fHsm", "SM dep", 18,-.5,17.5);
  fHetantr = new TH2F("fHetantr", "PH2d TRD", 50, -1, 1, 7, -.5, 6.5);
  fHphintr = new TH2F("fHphintr", "PH2d TRD", 157, 0,6.28, 7, -.5, 6.5);
  fHcltime = new TH1F("fHcltime", "Cluster-time", 30,0,30);
  fHcldiff = new TH1F("fHcldiff", "Cluster time diff", 30,0,30);
  fHxyA = new TH2F("fHxyA","x-y side A",800,-400,400,800,-400,400);
  fHxyC = new TH2F("fHxyC","x-y side C",800,-400,400,800,-400,400);
  fHPropagXY = new TH2F("fHPropagXY","x-y propagated point to TRD",1000,-500,500,1000,-500,500);
  fHPropagRZ = new TH2F("fHPropagRZ","r-z propagated point to TRD",700, -350., 350.,500,0.,500.);
  fHnFriendTracks = new TH1F("fHnFriendTracks", "NFriendTracks", 100, 0., 100.);
  fHnCalibObjects = new TH1F("fHnCalibObjects", "NCalibObjects", 100, 0., 100.);
  fHTpcEtaPhiLocaln = new TH2F("fHTpcEtaPhiLocaln", "(#eta - #varphi_{local}) distribution of neg TPC tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTrdEtaPhiLocaln = new TH2F("fHTrdEtaPhiLocaln", "(#eta - #varphi_{local}) distribution of neg TRD tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTpcEtaPhiLocalp = new TH2F("fHTpcEtaPhiLocalp", "(#eta - #varphi_{local}) distribution of pos TPC tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTrdEtaPhiLocalp = new TH2F("fHTrdEtaPhiLocalp", "(#eta - #varphi_{local}) distribution of pos TRD tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTpcEtaPhiSagitan = new TH2F("fHTpcEtaPhiSagitan", "(#eta - #varphi_{det}) distribution of neg TPC tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTrdEtaPhiSagitan = new TH2F("fHTrdEtaPhiSagitan", "(#eta - #varphi_{det}) distribution of neg TRD tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTpcEtaPhiSagitap = new TH2F("fHTpcEtaPhiSagitap", "(#eta - #varphi_{det}) distribution of pos TPC tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTrdEtaPhiSagitap = new TH2F("fHTrdEtaPhiSagitap", "(#eta - #varphi_{det}) distribution of pos TRD tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTpcEtaPhiDeltan = new TH2F("fHTpcEtaPhiDeltan", "(#eta - #Delta#varphi) distribution of neg TPC tracks", 
			  50, -1.0, 1.0, 50, -0.4*TMath::Pi(), 0.4*TMath::Pi());
  fHTrdEtaPhiDeltan = new TH2F("fHTrdEtaPhiDeltan", "(#eta - #Delta#varphi) distribution of neg TRD tracks", 
			  50, -1.0, 1.0, 50, -0.4*TMath::Pi(), 0.4*TMath::Pi());
  fHTpcEtaPhiDeltap = new TH2F("fHTpcEtaPhiDeltap", "(#eta - #Delta#varphi) distribution of pos TPC tracks", 
			  50, -1.0, 1.0, 50, -0.4*TMath::Pi(), 0.4*TMath::Pi());
  fHTrdEtaPhiDeltap = new TH2F("fHTrdEtaPhiDeltap", "(#eta - #Delta#varphi) distribution of pos TRD tracks", 
			  50, -1.0, 1.0, 50, -0.4*TMath::Pi(), 0.4*TMath::Pi());
  fHTpcEtaPhin = new TH2F("fHTpcEtaPhin", "(#eta - #varphi) distribution of neg TPC tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTrdEtaPhin = new TH2F("fHTrdEtaPhin", "(#eta - #varphi) distribution of neg TRD tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTpcEtaPhip = new TH2F("fHTpcEtaPhip", "(#eta - #varphi) distribution of pos TPC tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTrdEtaPhip = new TH2F("fHTrdEtaPhip", "(#eta - #varphi) distribution of pos TRD tracks", 
			  50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  Float_t etaBinLimits[51];
  for(Int_t i=0; i<51; i++) etaBinLimits[i] = -1.0 + i*2.0/50.;
  Float_t phiBinLimits[151];
  for(Int_t i=0; i<151; i++) phiBinLimits[i] = -1.1*TMath::Pi() + i*2.2*TMath::Pi()/150.;
  fHTpcRef3Dpos = new TH3F("fHTpcRef3Dpos", "TPC positive ref. tracks", 
			   50, etaBinLimits, 150, phiBinLimits, 32, binLimits);
  fHTpcRef3Dneg = new TH3F("fHTpcRef3Dneg", "TPC negative ref. tracks",
			   50, etaBinLimits, 150, phiBinLimits, 32, binLimits);
  fHTrdRef3Dpos = new TH3F("fHTrdRef3Dpos", "TRD positive tracks",
			   50, etaBinLimits, 150, phiBinLimits, 32, binLimits);
  fHTrdRef3Dneg = new TH3F("fHTrdRef3Dneg", "TRD negative tracks",
			   50, etaBinLimits, 150, phiBinLimits, 32, binLimits);
  fHTrdEtaPhiLocalNTracklets = new TProfile2D("fHTrdEtaPhiLocalNTracklets", "(#eta - #varphi_{local}) distribution of <ntracklets>", 
					      50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTrdEtaPhiSagitaNTracklets = new TProfile2D("fHTrdEtaPhiSagitaNTracklets", "(#eta - #varphi_{det}) distribution of <ntracklets>", 
					      50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  fHTrdEtaPhiDeltaNTracklets = new TProfile2D("fHTrdEtaPhiDeltaNTracklets", "(#eta - #Delta#varphi) distribution of <ntracklets>", 
					      50, -1.0, 1.0, 50, -0.4*TMath::Pi(), 0.4*TMath::Pi());
  fHTrdEtaPhiNTracklets = new TProfile2D("fHTrdEtaPhiNTracklets", "(#eta - #varphi) distribution of <ntracklets>", 
					 50, -1.0, 1.0, 150, -1.1*TMath::Pi(), 1.1*TMath::Pi());
  
  fOutStorage->Add(fHzvert1);
  fOutStorage->Add(fHzvert2);
  fOutStorage->Add(fHntracks);
  fOutStorage->Add(fHntracks2);
  fOutStorage->Add(fHntracks3);
  fOutStorage->Add(fHdca);
  fOutStorage->Add(fHdcaz);
  fOutStorage->Add(fHpt);
  fOutStorage->Add(fHpt2);
  fOutStorage->Add(fHpt3);
  fOutStorage->Add(fHpt3n);
  fOutStorage->Add(fHpt4);
  fOutStorage->Add(fHpt4n);
  fOutStorage->Add(fHtheta);
  fOutStorage->Add(fHphi);
  fOutStorage->Add(fHtpccl);
  fOutStorage->Add(fHtpccl2);
  fOutStorage->Add(fHdedxp);
  fOutStorage->Add(fHetaphi);
  fOutStorage->Add(fHetancl);
  fOutStorage->Add(fHphincl);
  fOutStorage->Add(fHtrdtr);
  fOutStorage->Add(fHtrdtr2);
  fOutStorage->Add(fHph);
  fOutStorage->Add(fHph2d);
  fOutStorage->Add(fHncltrkl);
  fOutStorage->Add(fHntrkl);
  fOutStorage->Add(fHntrklVsP);
  fOutStorage->Add(fHetantr);
  fOutStorage->Add(fHphintr);
  fOutStorage->Add(fHcltime);
  fOutStorage->Add(fHcldiff);
  fOutStorage->Add(fHxyA);
  fOutStorage->Add(fHxyC);
  fOutStorage->Add(fHPropagXY);
  fOutStorage->Add(fHPropagRZ);
  fOutStorage->Add(fHnFriendTracks);
  fOutStorage->Add(fHnCalibObjects);
  fOutStorage->Add(fHTpcEtaPhiLocaln);
  fOutStorage->Add(fHTpcEtaPhiLocalp);
  fOutStorage->Add(fHTrdEtaPhiLocaln);
  fOutStorage->Add(fHTrdEtaPhiLocalp);
  fOutStorage->Add(fHTpcEtaPhiSagitan);
  fOutStorage->Add(fHTpcEtaPhiSagitap);
  fOutStorage->Add(fHTrdEtaPhiSagitan);
  fOutStorage->Add(fHTrdEtaPhiSagitap);
  fOutStorage->Add(fHTpcEtaPhiDeltan);
  fOutStorage->Add(fHTpcEtaPhiDeltap);
  fOutStorage->Add(fHTrdEtaPhiDeltan);
  fOutStorage->Add(fHTrdEtaPhiDeltap);
  fOutStorage->Add(fHTpcEtaPhin);
  fOutStorage->Add(fHTpcEtaPhip);
  fOutStorage->Add(fHTrdEtaPhin);
  fOutStorage->Add(fHTrdEtaPhip);
  fOutStorage->Add(fHTpcRef3Dpos);
  fOutStorage->Add(fHTpcRef3Dneg);
  fOutStorage->Add(fHTrdRef3Dpos);
  fOutStorage->Add(fHTrdRef3Dneg);
  fOutStorage->Add(fHTrdEtaPhiLocalNTracklets);
  fOutStorage->Add(fHTrdEtaPhiSagitaNTracklets);
  fOutStorage->Add(fHTrdEtaPhiDeltaNTracklets);
  fOutStorage->Add(fHTrdEtaPhiNTracklets);
}

void AliAnalysisTaskTRDmon::Exec(Option_t *){
  //
  // Implementation of the Exec() function
  //
  if(!fESD){
	  AliError("ESD Event missing");
	  return;
  }
  //  const Float_t kCutDCAxy = 40.;      // cm
  const Float_t kCutDCAxy = 3; //40.;      // cm
  const Float_t kCutDCAz = 15.;      // cm
  const Int_t kNclTPCmin = 100;      // Nclusters TPC
  const Float_t kPmin = 0.2;      // min. pt
  
  // event trigger cut
  if(fEventTriggerName.Data()[0]!='\0' && !fESD->IsTriggerClassFired(fEventTriggerName.Data())) {
    PostData(0, fOutStorage);
    return; 
  }

  fESDfriend = dynamic_cast<AliESDfriend*> (fESD->FindListObject("AliESDfriend"));
  if (!fESDfriend) return;
  AliESDtrack *fTrack = 0x0;
  AliESDfriendTrack *fFriendTrack = 0x0;
  AliTRDtrackV1 *fTRDtrack=0x0; //TRD "calibration" object

  //const AliESDVertex *EvVertex = fESD->GetPrimaryVertexTPC();
  const AliESDVertex *evVertex = fESD->GetPrimaryVertex();
  Float_t zvert1 = evVertex->GetZ();
  Float_t nvtxContr = evVertex->GetNContributors();

  // if the required trigger is a collision trigger then apply event vertex cut
  if(fIsCollisionEvent && (TMath::Abs(zvert1)>15. || TMath::Abs(zvert1)<1.0e-10 || nvtxContr<=1)) {
    PostData(0, fOutStorage);
    return;
  }
  
  Float_t zvert2 = fESD->GetPrimaryVertex()->GetZ();
  
  fHzvert1->Fill(zvert1);
  fHzvert2->Fill(zvert2);
 
  Int_t ntracks = fESD->GetNumberOfTracks();
  fHntracks->Fill((Float_t)ntracks);
  
  Int_t ntracksDca=0;
  Int_t ntracksNcl=0;
  Int_t nFriendTracks = 0;
  Int_t nCalibObjects = 0;
  
  for(Int_t itrk = ntracks; itrk--;){ 
    fTrack = fESD->GetTrack(itrk);
    if(!fTrack) continue;
    //printf("%d, Label %d\n",itrk,fTrack->GetLabel()); // ?
    if(!(fTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
    if(fTrack->GetKinkIndex(0) > 0) continue; //rejects kinks
    
    Double_t pt = fTrack->Pt();
    Double_t p = fTrack->GetP();
    
    if ( pt < kPmin ) continue;
    
    Float_t theta=fTrack->Theta();
    Float_t phi=fTrack->Phi();
    Float_t eta=-TMath::Log(TMath::Tan(theta/2.));
    
    if ( TMath::Abs(eta)>0.9 || fTrack->GetTPCNcls()<10) continue;
    
    Int_t nTPCclus = fTrack->GetTPCNcls();
    fHetancl->Fill(eta,(Float_t)nTPCclus);	    
    fHphincl->Fill(phi,(Float_t)nTPCclus);	    
    fHpt->Fill(pt);
    
    Float_t imPar[2], cov[3];
    fTrack->GetImpactParameters(imPar, cov);
    //Float_t dca = TMath::Sqrt(imPar[0]*imPar[0] + imPar[1]*imPar[1]);
    Float_t dca =imPar[0]; //transverse
    Float_t dcaz =imPar[1]; //z
    fHdca->Fill(dca);
    fHdcaz->Fill(dcaz);
    
    if (fIsCollisionEvent && TMath::Abs(dca) > kCutDCAxy ) continue;	
    if (fIsCollisionEvent && TMath::Abs(dcaz) > kCutDCAz ) continue;	
    
    fHpt2->Fill(p);
    ntracksDca++;
    
    Double_t dedxTPC=fTrack->GetTPCsignal();
    
    fHtpccl->Fill((Float_t)nTPCclus); //?
    if ( pt > 1 ) fHtpccl2->Fill((Float_t)nTPCclus);
    
    // track points for x-y histograms
    fFriendTrack = (AliESDfriendTrack*)fTrack->GetFriendTrack(); //this works
    Bool_t hasFriends = (fFriendTrack ? kTRUE : kFALSE);
    if(hasFriends) nFriendTracks++;
    const AliTrackPointArray *array = (hasFriends ? fFriendTrack->GetTrackPointArray() : 0);
    if (array) {
      Int_t npoints = array->GetNPoints();
      for (int ipoint=0; ipoint<npoints; ipoint++) {
	AliTrackPoint tp;
	array->GetPoint(tp, ipoint);
	if (tp.GetZ()<0) fHxyC->Fill(tp.GetX(),tp.GetY());
	if (tp.GetZ()>0) fHxyA->Fill(tp.GetX(),tp.GetY());
      }
    }
    
    if ( nTPCclus < kNclTPCmin ) continue; 
    
    //    Bool_t propagGood = fTrack->PropagateTo(298., fESD->GetMagneticField());
    Double_t localCoord[3] = {0., 0., 0.};
    Double_t localMom[3]   = {0., 0., 0.};
    Bool_t localCoordGood = fTrack->GetXYZAt(298., fESD->GetMagneticField(), localCoord);
    if(localCoordGood) {
      fHPropagXY->Fill(localCoord[0], localCoord[1]);
      fHPropagRZ->Fill(localCoord[2], TMath::Sqrt(localCoord[0]*localCoord[0]+localCoord[1]*localCoord[1]));
    }
    Bool_t localMomGood   = fTrack->GetPxPyPzAt(298., fESD->GetMagneticField(), localMom);
    Double_t localPhi = (localMomGood ? TMath::ATan2(localMom[1], localMom[0]) : 0.0);
    Double_t localSagitaPhi = (localCoordGood ? TMath::ATan2(localCoord[1], localCoord[0]) : 0.0);    
    //    if(localPhi<0.) localPhi = TMath::Abs(localPhi)+TMath::Pi();   
    //    if(localSagitaPhi<0.) localSagitaPhi = TMath::Abs(localSagitaPhi)+TMath::Pi();     

    if (fTrack->Charge()>0) {
      fHpt3->Fill(p);
      if(localCoordGood && localMomGood) {
	fHTpcEtaPhiLocalp->Fill(eta,localPhi);
	fHTpcEtaPhiSagitap->Fill(eta,localSagitaPhi);
	fHTpcEtaPhiDeltap->Fill(eta,localPhi-localSagitaPhi);
	fHTpcRef3Dpos->Fill(eta, localSagitaPhi, p);
      }
      fHTpcEtaPhip->Fill(eta,phi-TMath::Pi());
    }
    if (fTrack->Charge()<0) {
      fHpt3n->Fill(p);
      if(localCoordGood && localMomGood) {
	fHTpcEtaPhiLocaln->Fill(eta,localPhi);
	fHTpcEtaPhiSagitan->Fill(eta,localSagitaPhi);
	fHTpcEtaPhiDeltan->Fill(eta,localPhi-localSagitaPhi);
	fHTpcRef3Dneg->Fill(eta, localSagitaPhi, p);
      }
      fHTpcEtaPhin->Fill(eta,phi-TMath::Pi());
    }
    ntracksNcl++;
    fHetaphi->Fill(eta,phi);
    fHdedxp->Fill(p,dedxTPC);
    Int_t nTRDtrkl = fTrack->GetTRDntracklets();  
    
    if ( nTRDtrkl < 1 ) continue;
    if (fTrack->Charge()>0) {
      fHpt4->Fill(p);
      if(localCoordGood && localMomGood) {
	fHTrdEtaPhiLocalp->Fill(eta,localPhi);
	fHTrdEtaPhiSagitap->Fill(eta,localSagitaPhi);
	fHTrdEtaPhiDeltap->Fill(eta,localPhi-localSagitaPhi);
	fHTrdRef3Dpos->Fill(eta, localSagitaPhi, p);
      }
      fHTrdEtaPhip->Fill(eta,phi-TMath::Pi());
    }
    if (fTrack->Charge()<0) {
      fHpt4n->Fill(p);
      if(localCoordGood && localMomGood) {
	fHTrdEtaPhiLocaln->Fill(eta,localPhi);
	fHTrdEtaPhiSagitan->Fill(eta,localSagitaPhi);
	fHTrdEtaPhiDeltan->Fill(eta,localPhi-localSagitaPhi);
	fHTrdRef3Dneg->Fill(eta, localSagitaPhi, p);
      }
      fHTrdEtaPhin->Fill(eta,phi-TMath::Pi());
    }
    fHetantr->Fill(eta,(Float_t)nTRDtrkl);	    
    fHphintr->Fill(phi,(Float_t)nTRDtrkl);	    
    if(localCoordGood && localMomGood) {
      fHTrdEtaPhiLocalNTracklets->Fill(eta, localPhi, (Float_t)nTRDtrkl);
      fHTrdEtaPhiSagitaNTracklets->Fill(eta, localSagitaPhi, (Float_t)nTRDtrkl);
      fHTrdEtaPhiDeltaNTracklets->Fill(eta, localPhi-localSagitaPhi, (Float_t)nTRDtrkl);
    }
    fHTrdEtaPhiNTracklets->Fill(eta, phi-TMath::Pi(), (Float_t)nTRDtrkl);
    fHtrdtr->Fill((Float_t)nTRDtrkl);
    if ( pt > 1 ) fHtrdtr2->Fill((Float_t)nTRDtrkl);
    fHtheta->Fill(theta);
    fHphi->Fill(phi);
    
    // start looking into TRD ...
    TObject *o = 0x0;
    Int_t ical = 0; 
    fTRDtrack = 0x0;
    while(hasFriends && (o = fFriendTrack->GetCalibObject(ical++))){
      //while((o = (const_cast<AliESDfriendTrack*>(fFriendTrack))->GetCalibObject(ical++))) { //...this is fishy...
      //printf("Object type: %s\n", o->IsA()->GetName());
      if(strcmp(o->IsA()->GetName(),"AliTRDtrackV1") != 0){
	continue;
      }
      fTRDtrack = dynamic_cast<AliTRDtrackV1 *>(o);
      break;
    }
    
    Bool_t hasTRDfriends = (fTRDtrack ? kTRUE : kFALSE);
    //    if(!fTRDtrack) continue;
    if(hasTRDfriends) nCalibObjects++;
    Int_t nTracklets = 0;
    Int_t tbinPrev=0;
    
    //    Int_t ntls = fTRDtrack->GetNumberOfTracklets();
    Int_t ntls = fTrack->GetTRDntracklets();
    //Int_t ncls = track->GetNumberOfClusters(); //per track
    if(hasTRDfriends) {
      for(Int_t itl = 0; itl < 6; itl++){ // TRD layers
	AliTRDseedV1 * tracklet = fTRDtrack->GetTracklet(itl);
	if(!tracklet || !tracklet->IsOK()) continue;
	nTracklets++;
	Int_t nclsTracklet = 0;
	for(Int_t itb = 0; itb < AliTRDseedV1::kNtb; itb++){ // timebins
	  AliTRDcluster *c = tracklet->GetClusters(itb);
	  AliTRDcluster *c1 = tracklet->GetClusters(itb + AliTRDseedV1::kNtb);  // second pad row
	  //AliTRDcluster *c = tracklet->GetdQdl(itb);
	  if(!(c || c1)) continue;
	  AliTRDcluster *cptr = c ? c : c1;
	  Int_t tcal = cptr->GetLocalTimeBin();
	  fHcldiff->Fill((Float_t)(tcal-tbinPrev)); 
	  tbinPrev=tcal;
	  Int_t idet = cptr->GetDetector();
	  //If (itrk<100 && itb==0) cout << "...Detector  "<<idet<< endl;
	  Int_t sm = idet/30;
	  fHsm->Fill((Float_t)sm); 
	  idet -= sm * 30;
	  float sig = 0;
	  if(c) sig += TMath::Abs(c->GetQ()); //
	  if(c1) sig += TMath::Abs(c1->GetQ());
	  
	  fHph->Fill((Float_t)tcal, sig);  
	  fHph2d->Fill((Float_t)tcal, sig);  
	  fHcltime->Fill((Float_t)tcal); 
	  nclsTracklet++;
	}
	fHncltrkl->Fill(nclsTracklet);
      } //tracklets
    } // end if(hasTRDfriends)
    //fHntrkl->Fill((Float_t)nTracklets);
    fHntrkl->Fill((Float_t)ntls);
    fHntrklVsP->Fill(p,(Float_t)ntls);
  } //tracks
  
  fHntracks2->Fill((Float_t)ntracksDca);
  fHntracks3->Fill((Float_t)ntracksNcl);
  fHnFriendTracks->Fill(nFriendTracks);
  fHnCalibObjects->Fill(nCalibObjects);
  
  PostData(0, fOutStorage);
}
void AliAnalysisTaskTRDmon::Terminate(Option_t *)
{    
  //
  // Do nothing
  //
}

