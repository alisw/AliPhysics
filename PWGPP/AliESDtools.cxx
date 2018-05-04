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
    author marian.ivanov@cern.ch
*/


/*
  TFile * f = TFile::Open("/lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/ATO-436/data/LHC15o_pass1_esds/alice/data/2015/LHC15o/000246272/pass1/15000246272021.7701/AliESDs.root");
  TTree * tree = (TTree*)f->Get("esdTree");
  // tree = AliXRDPROOFtoolkit::MakeChainRandom("esd.list","esdTree",0,10)
  .L $AliPhysics_SRC/PWGPP/AliESDtools.cxx+
  AliESDtools tools;
  tools.Init(tree);
   tree->GetEntry(3);
  tree->Scan("Tracks@.GetEntries():SPDVertex.fNContributors:Entry$","SPDVertex.fNContributors>100&&Tracks@.GetEntries()/PrimaryVertex.fNContributors>10")
tools->CacheTPCEventInformation()

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
#include "AliESDEvent.h"
#include "AliESDtools.h"

ClassImp(AliESDtools)
AliESDtools fgktools;

AliESDtools::AliESDtools():
  fVerbose(0),
  fESDtree(NULL),
  fEvent(NULL),
  fHisTPCVertexA(nullptr),
  fHisTPCVertexC(nullptr),
  fHisTPCVertexACut(nullptr),
  fHisTPCVertexCCut(nullptr),
  fHisTPCVertex(nullptr),
  fHistPhiTPCcounterA(nullptr),         // helper histogram phi counteres
  fHistPhiTPCcounterC(nullptr),         // helper histogram for TIdentity tree
  fHistPhiTPCcounterAITS(nullptr),      // helper histogram for TIdentity tree
  fHistPhiTPCcounterCITS(nullptr),      // helper histogram for TIdentity tree
  fHistPhiITScounterA(nullptr),         // helper histogram for TIdentity tree
  fHistPhiITScounterC(nullptr),         // helper histogram for TIdentity tree
  fCacheTrackCounters(nullptr),         // track counter
  fCacheTrackTPCCountersZ(nullptr),         // track counter
  fCacheTrackdEdxRatio(nullptr),        // dedx info counter
  fCacheTrackNcl(nullptr),              // ncl counter
  fCacheTrackChi2(nullptr),             // chi2 counter
  fCacheTrackMatchEff(nullptr),         // matchEff counter
  fLumiGraph(nullptr)                  // graph for the interaction rate info for a run
{

}

void AliESDtools::Init(TTree *tree) {
  AliESDtools & tools = *this;
  if (tools.fESDtree) delete tools.fESDtree;
  if (!tools.fEvent) tools.fEvent = new AliESDEvent();
  tools.fESDtree = tree;
  tools.fEvent->ReadFromTree(tree);
  if (fHisTPCVertexA == nullptr) {
    tools.fHisTPCVertexA = new TH1F("hisTPCZA", "hisTPCZA", 1000, -250, 250);
    tools.fHisTPCVertexC = new TH1F("hisTPCZC", "hisTPCZC", 1000, -250, 250);
    tools.fHisTPCVertex = new TH1F("hisTPCZ", "hisTPCZ", 1000, -250, 250);
    tools.fHisTPCVertexACut = new TH1F("hisTPCZACut", "hisTPCZACut", 1000, -250, 250);
    tools.fHisTPCVertexCCut = new TH1F("hisTPCZCCut", "hisTPCZCCut", 1000, -250, 250);
    tools.fHisTPCVertex->SetLineColor(1);
    tools.fHisTPCVertexA->SetLineColor(2);
    tools.fHisTPCVertexC->SetLineColor(4);
    tools.fHisTPCVertexACut->SetLineColor(3);
    tools.fHisTPCVertexCCut->SetLineColor(6);
    tools.fCacheTrackCounters = new TVectorF(20);
    tools.fCacheTrackTPCCountersZ = new TVectorF(8);
    tools.fCacheTrackdEdxRatio = new TVectorF(20);
    tools.fCacheTrackNcl = new TVectorF(20);
    tools.fCacheTrackChi2 = new TVectorF(20);
    tools.fCacheTrackMatchEff = new TVectorF(20);
    // **************************** Event histograms **************************
    fHistPhiTPCcounterA = new TH1F("hPhiTPCcounterC", "control histogram to count tracks on the A side in phi ", 36, 0., 18.);
    fHistPhiTPCcounterC = new TH1F("hPhiTPCcounterA", "control histogram to count tracks on the C side in phi ", 36, 0., 18.);
    fHistPhiTPCcounterAITS = new TH1F("hPhiTPCcounterAITS", "control histogram to count tracks on the A side in phi ", 36, 0., 18.);
    fHistPhiTPCcounterCITS = new TH1F("hPhiTPCcounterCITS", "control histogram to count tracks on the C side in phi ", 36, 0., 18.);
    fHistPhiITScounterA = new TH1F("hPhiITScounterA", "control histogram to count tracks on the A side in phi ", 36, 0., 18.);
    fHistPhiITScounterC = new TH1F("hPhiITScounterC", "control histogram to count tracks on the C side in phi ", 36, 0., 18.);


  }
}

/// cache TPC event information
/// \return
Int_t AliESDtools::CacheTPCEventInformation(){
  AliESDtools &tools=*this;
  const Int_t kNCRCut=80;
  const Double_t kDCACut=5;
  const Float_t knTrackletCut=1.5;
  // FILL DCA histograms
  tools.fHisTPCVertexA->Reset();
  tools.fHisTPCVertexC->Reset();
  tools.fHisTPCVertexACut->Reset();
  tools.fHisTPCVertexCCut->Reset();
  tools.fHisTPCVertex->Reset();
  Int_t nTracks=tools.fEvent->GetNumberOfTracks();
  Int_t selected=0;
  for (Int_t iTrack=0; iTrack<nTracks; iTrack++){
    AliESDtrack * pTrack = tools.fEvent->GetTrack(iTrack);
    Float_t dcaxy,dcaz;
    if (pTrack== nullptr) continue;
    if (pTrack->IsOn(AliVTrack::kTPCin)==0) continue;
    if (pTrack->GetTPCClusterInfo(3,1)<kNCRCut) continue;
    pTrack->GetImpactParameters(dcaxy,dcaz);
    if (TMath::Abs(dcaxy)>kDCACut) continue;
    selected++;
    if ((pTrack->GetNumberOfTRDClusters()/20.+pTrack->GetNumberOfITSClusters())>knTrackletCut){
      tools.fHisTPCVertex->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()>0) tools.fHisTPCVertexACut->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()<0) tools.fHisTPCVertexCCut->Fill(pTrack->GetTPCInnerParam()->GetZ());
    }else{
      if (pTrack->GetTgl()>0) tools.fHisTPCVertexA->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()<0) tools.fHisTPCVertexC->Fill(pTrack->GetTPCInnerParam()->GetZ());

    }
  }
  TPCVertexFit(tools.fHisTPCVertex);
  TPCVertexFit(tools.fHisTPCVertexA);
  TPCVertexFit(tools.fHisTPCVertexC);
  TPCVertexFit(tools.fHisTPCVertexACut);
  TPCVertexFit(tools.fHisTPCVertexCCut);
  if (fVerbose&0x10) printf("%d\n",selected);
  //
  return selected;
}

void AliESDtools::TPCVertexFit(TH1F *hisVertex){
  //0.5 cm 1 bin
  // hisVertex=tools.fHisTPCVertexACut;
  TAxis * axis = hisVertex->GetXaxis();
  for (Int_t iBin=5; iBin<axis->GetNbins()-5; iBin++){
    Double_t median10=TMath::Median(10,&(hisVertex->GetArray()[iBin-5]));
    Double_t rms10=TMath::RMS(10,&(hisVertex->GetArray()[iBin-5]));
    Double_t val0=TMath::Mean(3,&(hisVertex->GetArray()[iBin-2]));
    Double_t val1=TMath::Mean(3,&(hisVertex->GetArray()[iBin-1]));
    Double_t val2=TMath::Mean(3,&(hisVertex->GetArray()[iBin+0]));
    if (val1>=val0 && val1>=val2 && val1>3+1.5*median10){
      Double_t xbin=axis->GetBinCenter(iBin);
      printf("Ibin %d\t%f\t%f\t%f\t%f\n", iBin,xbin,val1,median10,rms10);
      hisVertex->Fit("gaus","qnrsame+","qnr",xbin-2,xbin+2);
    }
  }
}



//________________________________________________________________________
void AliESDtools::CalculateEventVariables()
{

  //AliVEvent *event=InputEvent();
  fCacheTrackCounters->Zero();   // track counter
  fCacheTrackCounters->Zero();   // track counter
  fCacheTrackdEdxRatio->Zero();; // dedx info counter
  fCacheTrackNcl->Zero();;       // ncl counter
  fCacheTrackChi2->Zero();;      // chi2 counter
  fCacheTrackMatchEff->Zero();;  // matchEff counter
  //
  //
  if (fHistPhiTPCcounterA) fHistPhiTPCcounterA->Reset();
  if (fHistPhiTPCcounterC) fHistPhiTPCcounterC->Reset();
  if (fHistPhiTPCcounterAITS) fHistPhiTPCcounterAITS->Reset();
  if (fHistPhiTPCcounterCITS) fHistPhiTPCcounterCITS->Reset();
  if (fHistPhiITScounterA) fHistPhiITScounterA->Reset();
  if (fHistPhiITScounterC) fHistPhiITScounterC->Reset();
  //
  //
  const Int_t kNclTPCcut=60;
  const Int_t kTglCut=1.5;
  const Int_t kDCACut=5;  // 5 cm primary cut
  const Int_t kMindEdxClustersRegion=15;
  const Int_t kPtCut=0.100;
  //
  //
  AliTPCdEdxInfo tpcdEdxInfo;
  for (Int_t itrack=0;itrack<fEvent->GetNumberOfTracks();++itrack)
  {   // Track loop

    AliESDtrack *track = fEvent->GetTrack(itrack);
    // AliESDPid *pid = track->GetDetPid();
    Double_t eta=-100., phiTPC=0.,sectorNumber=0.;
    Double_t tgl     = track->Pz()/track->Pt();
    Double_t phi  = track->Phi()-TMath::Pi();
    phi = track->GetParameterAtRadius(85,5,7);
    Double_t sectorNumbertmp = (9*phi/TMath::Pi()+18*(phi<0));
    if (track == NULL) continue;
    eta = track->Eta();
    if (TMath::Abs(eta)>0.9) continue;
    Bool_t isOnITS = track->IsOn(AliESDtrack::kITSrefit);
    Bool_t isOnTRD = track->IsOn(AliESDtrack::kTRDrefit);
    Bool_t isOnTPC = track->IsOn(AliESDtrack::kTPCrefit);
    //
    //
    if (track->GetInnerParam()) {
      phiTPC = track->GetInnerParam()->GetParameterAtRadius(85,5,7);
      sectorNumber = (9*phiTPC/TMath::Pi()+18*(phiTPC<0));
    }
    //
    // Count only ITS tracks
    if ( isOnITS && !isOnTPC ) {
      if (TMath::Abs(phi)>1e-10){
        if (tgl>0) fHistPhiITScounterA->Fill(sectorNumbertmp);
        if (tgl<0) fHistPhiITScounterC->Fill(sectorNumbertmp);
      }
    }
    //
    if (!track->GetInnerParam()) continue;
    if (track->IsOn(AliVTrack::kTPCout)==kFALSE)  continue;
    (*fCacheTrackCounters)[4]++;      // all TPC track with out flag
    // TPC track counters with DCAZ
    for (Int_t izCut=1; izCut<4; izCut++){
      Float_t impactParam[2];
      track->GetImpactParameters(impactParam[0],impactParam[1]);
      if (TMath::Abs(impactParam[0])>kDCACut) continue;
      if (TMath::Abs(track->GetInnerParam()->GetParameter()[1])<10.*(izCut+1.)) (*fCacheTrackTPCCountersZ)[izCut]++;
      if (TMath::Abs(impactParam[1])<10.*(izCut+1.)) (*fCacheTrackTPCCountersZ)[izCut+4]++;
    }
    //
    //
    Float_t dcaRPhi, dcaZ;
    track->GetImpactParameters(dcaRPhi, dcaZ);
    Int_t nclTPC    = track->GetTPCncls();
    Int_t nclITS    = track->GetITSNcls();
    Int_t nclTRD    = track->GetTRDncls();
    Int_t nclTOF    = track->IsOn(AliVTrack::kTOFout);
    if (nclTRD<1) nclTRD=-1;
    if (nclITS<1) nclITS=-1;
    //if (fNcl<1) fNcl=-1;
    Double_t chi2TPC = TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/nclTPC));
    Double_t chi2ITS = TMath::Sqrt(TMath::Abs(track->GetITSchi2()));
    Double_t chi2TRD = TMath::Sqrt(TMath::Abs(track->GetTRDchi2()));
    Double_t ptot0   = track->GetP();
    Double_t qP      = track->Charge()/track->P();
    //
    //
    if (nclTPC<kNclTPCcut) continue;
    if (TMath::Abs(tgl)>kTglCut) continue;
    if (track->Pt()<kPtCut) continue;
    if (TMath::Abs(dcaRPhi)>kDCACut || TMath::Abs(dcaZ)>kDCACut) continue;
    (*fCacheTrackCounters)[5]++;
    if (TMath::Abs(phiTPC)>1e-10){
      if (tgl>0) fHistPhiTPCcounterA->Fill(sectorNumber);
      if (tgl<0) fHistPhiTPCcounterC->Fill(sectorNumber);
      if(isOnITS){
        if (tgl>0) fHistPhiTPCcounterAITS->Fill(sectorNumber);
        if (tgl<0) fHistPhiTPCcounterCITS->Fill(sectorNumber);
      }
    }
    //
    //
    Bool_t pileUpCut=  ( (nclITS>2) || (nclTRD>40));
    if (pileUpCut==kFALSE) continue;
    if (TMath::Min(chi2TPC,100.)<0) continue;
    (*fCacheTrackCounters)[1]++;

    //
    Bool_t itsOK=track->IsOn(AliVTrack::kITSout) && nclITS>2  && chi2ITS>0;
    Bool_t trdOK=track->IsOn(AliVTrack::kTRDout) && nclTRD>35 && chi2TRD>0;
    Bool_t tofOK=track->IsOn(AliVTrack::kTOFout);
    //
    (*fCacheTrackNcl)[4]+=track->GetTPCncls(0, 63);
    (*fCacheTrackNcl)[5]+=track->GetTPCncls(64, 127);
    (*fCacheTrackNcl)[6]+=track->GetTPCncls(128, 159);
    (*fCacheTrackNcl)[1] += nclTPC;
    (*fCacheTrackChi2)[1]+= (chi2TPC>0) ? TMath::Sqrt(chi2TPC):2;   // sometimes negative chi2?

    if (itsOK && track->GetTPCdEdxInfo(tpcdEdxInfo)){

      Bool_t isOK=(tpcdEdxInfo.GetNumberOfCrossedRows(0)>kMindEdxClustersRegion);
      isOK&=(tpcdEdxInfo.GetNumberOfCrossedRows(1)>kMindEdxClustersRegion);
      isOK&=(tpcdEdxInfo.GetNumberOfCrossedRows(2)>kMindEdxClustersRegion);
      isOK&=((tpcdEdxInfo.GetSignalMax(0)>0) && (tpcdEdxInfo.GetSignalMax(1)>0) && (tpcdEdxInfo.GetSignalMax(2)>0));
      isOK&=((tpcdEdxInfo.GetSignalTot(0)>0) && (tpcdEdxInfo.GetSignalTot(1)>0) && (tpcdEdxInfo.GetSignalTot(2)>0));
      isOK&=(itsOK||trdOK);      // stronger pile-up cut requiring ITS or TRD

      if (isOK) {
        (*fCacheTrackCounters)[6]+=1;         // Counter with accepted TPC dEdx info
        (*fCacheTrackdEdxRatio)[0]+=TMath::Log(tpcdEdxInfo.GetSignalMax(3));
        (*fCacheTrackdEdxRatio)[1]+=TMath::Log(tpcdEdxInfo.GetSignalTot(3));
        (*fCacheTrackdEdxRatio)[2]+=TMath::Log(tpcdEdxInfo.GetSignalMax(0)/tpcdEdxInfo.GetSignalTot(0));
        (*fCacheTrackdEdxRatio)[3]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalTot(1));
        (*fCacheTrackdEdxRatio)[4]+=TMath::Log(tpcdEdxInfo.GetSignalMax(2)/tpcdEdxInfo.GetSignalTot(2));
        (*fCacheTrackdEdxRatio)[5]+=TMath::Log(tpcdEdxInfo.GetSignalMax(3)/tpcdEdxInfo.GetSignalTot(3));
        (*fCacheTrackdEdxRatio)[6]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalMax(0));
        (*fCacheTrackdEdxRatio)[7]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalMax(2));
        (*fCacheTrackdEdxRatio)[8]+=TMath::Log(tpcdEdxInfo.GetSignalTot(1)/tpcdEdxInfo.GetSignalTot(0));
        (*fCacheTrackdEdxRatio)[9]+=TMath::Log(tpcdEdxInfo.GetSignalTot(1)/tpcdEdxInfo.GetSignalTot(2));
      }
    }

    if (itsOK) {  // ITS
      (*fCacheTrackCounters)[0]++;
      (*fCacheTrackNcl)[0] += nclITS;
      (*fCacheTrackChi2)[0] += TMath::Min(TMath::Sqrt(chi2ITS),10.); // cutoff chi2 10
      (*fCacheTrackMatchEff)[2]+=trdOK;
      (*fCacheTrackMatchEff)[3]+=tofOK;
      (*fCacheTrackChi2)[4]+= (chi2TPC>0) ? TMath::Sqrt(chi2TPC):2; // TPC chi2 in case prolongation to ITS
      // long tracks properties
      if (nclITS>4){
        (*fCacheTrackCounters)[7]++;
        (*fCacheTrackNcl)[7] += nclITS;
        (*fCacheTrackChi2)[7]+=TMath::Min(TMath::Sqrt(chi2ITS),10.);
      }
    }
    if (trdOK) {// TRD    ///TODO - why chi2TRD could be smaller than 0?
      (*fCacheTrackCounters)[2]++;
      (*fCacheTrackNcl)[2] += nclTRD;
      (*fCacheTrackChi2)[2] += TMath::Sqrt(chi2TRD);
      (*fCacheTrackMatchEff)[0]+=itsOK;
      (*fCacheTrackChi2)[5]+= (chi2TPC>0) ? TMath::Sqrt(chi2TPC):2; // TPC chi2 in case prolongation to TRD
      if (nclTRD>80){
        (*fCacheTrackCounters)[8]++;
        (*fCacheTrackNcl)[8] += nclTRD;
        (*fCacheTrackChi2)[8]+=TMath::Min(TMath::Sqrt(chi2TRD),10.);
      }
    }
    if (tofOK) {  // TOF
      (*fCacheTrackCounters)[3]++;
      (*fCacheTrackNcl)[3] += 1;   // dummy for the moment
      (*fCacheTrackChi2)[3]+= 1;   //
    }
  } // end of track LOOP

  for (Int_t i=0; i<9; i++) if ((*fCacheTrackCounters)[i]>0) (*fCacheTrackNcl)[i]/=(*fCacheTrackCounters)[i];
  for (Int_t i=0; i<4; i++) if ((*fCacheTrackCounters)[i]>0) (*fCacheTrackChi2)[i]/=(*fCacheTrackCounters)[i];

  for (Int_t i=4; i<7; i++)  if ((*fCacheTrackCounters)[1]>0) (*fCacheTrackNcl)[i]/=(*fCacheTrackCounters)[1];
  if ((*fCacheTrackCounters)[6]>0){
    for (Int_t i=0; i<10; i++)   (*fCacheTrackdEdxRatio)[i]/=(*fCacheTrackCounters)[6];
  }
  //
  // conditional matching efficiency and chi2
  if ((*fCacheTrackCounters)[0]>0){
    (*fCacheTrackMatchEff)[2]/=(*fCacheTrackCounters)[0];  // TRD if ITS
    (*fCacheTrackMatchEff)[3]/=(*fCacheTrackCounters)[0];  // TOF if ITS
    (*fCacheTrackChi2)[4]/=(*fCacheTrackCounters)[0];
  }
  if ((*fCacheTrackCounters)[2]>0) {
    (*fCacheTrackMatchEff)[0]/=(*fCacheTrackCounters)[2];
    (*fCacheTrackChi2)[5]/=(*fCacheTrackCounters)[2];
  } //ITS if TRD
  (*fCacheTrackCounters)[9]=fEvent->GetNumberOfTracks();  // original number of ESDtracks

}