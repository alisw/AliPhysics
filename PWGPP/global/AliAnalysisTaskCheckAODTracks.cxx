#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODEvent.h"
#include "AliStack.h"
#include "AliPID.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAODTracklets.h"
#include "AliMultSelection.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include <AliAODMCParticle.h>
#include "AliAnalysisUtils.h"
#include <TSystem.h>
#include <TTree.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TChain.h>
#include "AliPIDResponse.h"
#include "AliAnalysisTaskCheckAODTracks.h"


/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Implementation of class AliAnalysisTaskCheckAODTracks
// AliAnalysisTaskSE to extract QA and performance histos for tracks
// 
//
// Authors: 
//          F. Prino, prino@to.infn.it
//          
//*************************************************************************

ClassImp(AliAnalysisTaskCheckAODTracks)
//______________________________________________________________________________
AliAnalysisTaskCheckAODTracks::AliAnalysisTaskCheckAODTracks() : 
  AliAnalysisTaskSE("QAofAODtracks"), 
  fOutput{nullptr},
  fHistNEvents{nullptr},
  fHistNTracks{nullptr},
  fHistNTracksVsTPCclusters{nullptr},
  fHistNTracksVsITSclusters{nullptr},
  fHistITSclustersVsTPCclusters{nullptr},
  fHistNTracksFB4VsTPCclusters{nullptr},
  fHistNTracksFB4VsITSclusters{nullptr},
  fHistFilterBits{nullptr},
  fHistITSnClusTPCsel{nullptr},
  fHistITSnClusITSsa{nullptr},
  fHistITSnClusITSPureSA{nullptr},
  fHistITSCluInLayTPCsel{nullptr},
  fHistITSCluInLayITSsa{nullptr},
  fHistITSCluInLayITSPureSA{nullptr},
  fHistNtracksFb4VsV0befEvSel{nullptr},
  fHistNtracksFb5VsV0befEvSel{nullptr},
  fHistNtracksFb4VsV0aftEvSel{nullptr},
  fHistNtracksFb5VsV0aftEvSel{nullptr},
  fHistNtracksFb0{nullptr},
  fHistNtracksFb1{nullptr},
  fHistNtracksFb4{nullptr},
  fHistNtracksFb5{nullptr},
  fHistNtracksFb6{nullptr},
  fHistNtracksFb7{nullptr},
  fHistNtracksFb8{nullptr},
  fHistEtaPhiPtTPCsel{nullptr},
  fHistEtaPhiPtTPCselITSref{nullptr},
  fHistEtaPhiPtTPCselSPDany{nullptr},
  fHistEtaPhiPtPosChargeTPCsel{nullptr},
  fHistEtaPhiPtPosChargeTPCselITSref{nullptr},
  fHistEtaPhiPtPosChargeTPCselSPDany{nullptr},
  fHistEtaPhiPtNegChargeTPCsel{nullptr},
  fHistEtaPhiPtNegChargeTPCselITSref{nullptr},
  fHistEtaPhiPtNegChargeTPCselSPDany{nullptr},
  fHistEtaPhiPtTPCselTOFbc{nullptr},
  fHistEtaPhiPtTPCselITSrefTOFbc{nullptr},
  fHistEtaPhiPtTPCselSPDanyTOFbc{nullptr},
  fHistEtaPhiPtTPCselSPDanyTOFpid{nullptr},
  fHistTPCchi2PerClusPhiPtTPCsel{nullptr},
  fHistTPCchi2PerClusPhiPtTPCselITSref{nullptr},
  fHistTPCchi2PerClusPhiPtTPCselSPDany{nullptr},
  fHistSig1ptCovMatPhiPtTPCsel{nullptr},
  fHistSig1ptCovMatPhiPtTPCselITSref{nullptr},
  fHistSig1ptCovMatPhiPtTPCselSPDany{nullptr},
  fHistImpParXYPtMulPionTPCselSPDany{nullptr},
  fHistImpParXYPtMulKaonTPCselSPDany{nullptr},
  fHistImpParXYPtMulProtonTPCselSPDany{nullptr},
  fHistPtResidVsPtTPCselAll{nullptr},
  fHistPtResidVsPtTPCselITSrefAll{nullptr},
  fHistOneOverPtResidVsPtTPCselAll{nullptr},
  fHistOneOverPtResidVsPtTPCselITSrefAll{nullptr},
  fHistPtResidVsPtTPCselPrim{nullptr},
  fHistPtResidVsPtTPCselITSrefPrim{nullptr},
  fHistOneOverPtResidVsPtTPCselPrim{nullptr},
  fHistOneOverPtResidVsPtTPCselITSrefPrim{nullptr},
  fHistPtResidVsPtTPCselSecDec{nullptr},
  fHistPtResidVsPtTPCselITSrefSecDec{nullptr},
  fHistOneOverPtResidVsPtTPCselSecDec{nullptr},
  fHistOneOverPtResidVsPtTPCselITSrefSecDec{nullptr},
  fHistPtResidVsPtTPCselSecMat{nullptr},
  fHistPtResidVsPtTPCselITSrefSecMat{nullptr},
  fHistOneOverPtResidVsPtTPCselSecMat{nullptr},
  fHistOneOverPtResidVsPtTPCselITSrefSecMat{nullptr},
  fHistOneOverPtResidVsPtTPCsel{nullptr},
  fHistOneOverPtResidVsPtTPCselITSref{nullptr},
  fHistPzResidVsPtTPCselAll{nullptr},
  fHistPzResidVsPtTPCselITSrefAll{nullptr},
  fHistPzResidVsEtaTPCselAll{nullptr},
  fHistPzResidVsEtaTPCselITSrefAll{nullptr},
  fHistEtaPhiPtTPCselITSrefGood{nullptr},
  fHistEtaPhiPtTPCselITSrefFake{nullptr}, 
  fHistImpParXYPtMulTPCselSPDanyGood{nullptr},
  fHistImpParXYPtMulTPCselSPDanyFake{nullptr},
  fHistImpParXYPtMulTPCselSPDanyPrim{nullptr},
  fHistImpParXYPtMulTPCselSPDanySecDec{nullptr},
  fHistImpParXYPtMulTPCselSPDanySecMat{nullptr},
  fHistInvMassK0s{nullptr},
  fHistInvMassLambda{nullptr},
  fHistInvMassAntiLambda{nullptr},
  fFillTree(kFALSE),
  fTrackTree{nullptr},
  fTreeVarFloat{nullptr},
  fTreeVarInt{nullptr},
  fTrCutsTPC{nullptr},
  fMinNumOfTPCPIDclu(0),
  fUsePhysSel(kTRUE),
  fUsePileupCut(kTRUE),
  fTriggerMask(AliVEvent::kAnyINT),
  fSelectOnCentrality(kFALSE),
  fMinCentrality(-1.),
  fMaxCentrality(110.),
  fCentrEstimator("V0M"),
  fNEtaBins(10),
  fNPhiBins(144),
  fNPtBins(100),
  fMinPt(0.),
  fMaxPt(25.),
  fMaxMult(500.),
  fRequireITSforV0dau(0),
  fReadMC{kFALSE},
  fUseMCId{kFALSE},
  fUseGenPt{kFALSE}
{
  //

  fTrCutsTPC = new AliESDtrackCuts("esdtrackCutsTPC");
  fTrCutsTPC->SetMinNCrossedRowsTPC(50);
  fTrCutsTPC->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fTrCutsTPC->SetEtaRange(-0.8,0.8);
  // fTrCutsTPC->SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0); 
  // fTrCutsTPC->SetCutOutDistortedRegionsTPC(kTRUE);
  fTrCutsTPC->SetMaxChi2PerClusterTPC(4);
  fTrCutsTPC->SetAcceptKinkDaughters(kFALSE);
  fTrCutsTPC->SetRequireTPCRefit(kTRUE);
  fTrCutsTPC->SetDCAToVertex2D(kFALSE);
  fTrCutsTPC->SetRequireSigmaToVertex(kFALSE);
  //  fTrCutsTPC->SetMaxChi2TPCConstrainedGlobal(36);
  //  fTrCutsTPC->SetMaxChi2PerClusterITS(36);
  fTrCutsTPC->SetMaxDCAToVertexXY(2.);
  fTrCutsTPC->SetMaxDCAToVertexZ(3.);
  for(Int_t jb=0; jb<kNumOfFilterBits; jb++){
    fHistImpParXYPtMulFiltBit[jb]=0x0;
    fHistEtaPhiPtFiltBit[jb]=0x0;
    fHistITScluPtFiltBit[jb]=0x0;
    fHistSPDcluPtFiltBit[jb]=0x0;
    fHistTPCcluPtFiltBit[jb]=0x0;
    fHistTPCcrrowsPtFiltBit[jb]=0x0;
    fHistTPCCrowOverFindPtFiltBit[jb]=0x0;
    fHistTPCChi2clusPtFiltBit[jb]=0x0;
    fHistChi2TPCConstrVsGlobPtFiltBit[jb]=0x0;
    fHistSig1ptCovMatPtFiltBit[jb]=0x0;
  }
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


//___________________________________________________________________________
AliAnalysisTaskCheckAODTracks::~AliAnalysisTaskCheckAODTracks(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistNTracks;
    delete fHistNTracksVsTPCclusters;
    delete fHistNTracksVsITSclusters;
    delete fHistITSclustersVsTPCclusters;
    delete fHistNTracksFB4VsTPCclusters;
    delete fHistNTracksFB4VsITSclusters;
    delete fHistFilterBits;
    delete fHistITSnClusTPCsel;
    delete fHistITSnClusITSsa;
    delete fHistITSnClusITSPureSA;
    delete fHistITSCluInLayTPCsel;
    delete fHistITSCluInLayITSsa;
    delete fHistITSCluInLayITSPureSA;
    delete fHistNtracksFb4VsV0befEvSel;
    delete fHistNtracksFb5VsV0befEvSel;
    delete fHistNtracksFb4VsV0aftEvSel;
    delete fHistNtracksFb5VsV0aftEvSel;
    delete fHistNtracksFb0;
    delete fHistNtracksFb1;
    delete fHistNtracksFb4;
    delete fHistNtracksFb5;
    delete fHistNtracksFb6;
    delete fHistNtracksFb7;
    delete fHistNtracksFb8;
    delete fHistEtaPhiPtTPCsel;
    delete fHistEtaPhiPtTPCselITSref;
    delete fHistEtaPhiPtTPCselSPDany;
    delete fHistEtaPhiPtPosChargeTPCsel;
    delete fHistEtaPhiPtPosChargeTPCselITSref;
    delete fHistEtaPhiPtPosChargeTPCselSPDany;
    delete fHistEtaPhiPtNegChargeTPCsel;
    delete fHistEtaPhiPtNegChargeTPCselITSref;
    delete fHistEtaPhiPtNegChargeTPCselSPDany;
    delete fHistEtaPhiPtTPCselTOFbc;
    delete fHistEtaPhiPtTPCselITSrefTOFbc;
    delete fHistEtaPhiPtTPCselSPDanyTOFbc;
    delete fHistEtaPhiPtTPCselSPDanyTOFpid;
    delete fHistTPCchi2PerClusPhiPtTPCsel;
    delete fHistTPCchi2PerClusPhiPtTPCselITSref;
    delete fHistTPCchi2PerClusPhiPtTPCselSPDany;
    delete fHistSig1ptCovMatPhiPtTPCsel;
    delete fHistSig1ptCovMatPhiPtTPCselITSref;
    delete fHistSig1ptCovMatPhiPtTPCselSPDany;
    delete fHistImpParXYPtMulPionTPCselSPDany;
    delete fHistImpParXYPtMulKaonTPCselSPDany;
    delete fHistImpParXYPtMulProtonTPCselSPDany;
    delete fHistPtResidVsPtTPCselPrim;
    delete fHistPtResidVsPtTPCselITSrefPrim;
    delete fHistOneOverPtResidVsPtTPCselPrim;
    delete fHistOneOverPtResidVsPtTPCselITSrefPrim;
    delete fHistPtResidVsPtTPCselSecDec;
    delete fHistPtResidVsPtTPCselITSrefSecDec;
    delete fHistOneOverPtResidVsPtTPCselSecDec;
    delete fHistOneOverPtResidVsPtTPCselITSrefSecDec;
    delete fHistPtResidVsPtTPCselSecMat;
    delete fHistPtResidVsPtTPCselITSrefSecMat;
    delete fHistOneOverPtResidVsPtTPCselSecMat;
    delete fHistOneOverPtResidVsPtTPCselITSrefSecMat;
    for (int iS = 0; iS < AliPID::kSPECIESC;++iS) {
      delete fHistPtResidVsPtTPCsel[iS];
      delete fHistPtResidVsPtTPCselITSref[iS];
      delete fHistOneOverPtResidVsPtTPCsel[iS];
      delete fHistOneOverPtResidVsPtTPCselITSref[iS];
    }
    delete fHistPzResidVsPtTPCselAll;
    delete fHistPzResidVsPtTPCselITSrefAll;
    delete fHistPzResidVsEtaTPCselAll;
    delete fHistPzResidVsEtaTPCselITSrefAll;
    delete fHistEtaPhiPtTPCselITSrefGood;
    delete fHistEtaPhiPtTPCselITSrefFake;
    delete fHistImpParXYPtMulTPCselSPDanyGood;
    delete fHistImpParXYPtMulTPCselSPDanyFake;
    delete fHistImpParXYPtMulTPCselSPDanyPrim;
    delete fHistImpParXYPtMulTPCselSPDanySecDec;
    delete fHistImpParXYPtMulTPCselSPDanySecMat;
    delete fHistInvMassK0s;
    delete fHistInvMassLambda;
    delete fHistInvMassAntiLambda;

    for(Int_t jb=0; jb<kNumOfFilterBits; jb++){
      delete fHistImpParXYPtMulFiltBit[jb];
      delete fHistEtaPhiPtFiltBit[jb];
      delete fHistITScluPtFiltBit[jb];
      delete fHistSPDcluPtFiltBit[jb];
      delete fHistTPCcluPtFiltBit[jb];
      delete fHistTPCcrrowsPtFiltBit[jb];
      delete fHistTPCCrowOverFindPtFiltBit[jb];      
      delete fHistTPCChi2clusPtFiltBit[jb];
      delete fHistChi2TPCConstrVsGlobPtFiltBit[jb];
      delete fHistSig1ptCovMatPtFiltBit[jb];
    }
    delete fTrackTree;
  }
  delete fOutput;
  delete fTrCutsTPC;
  delete [] fTreeVarFloat;
  delete [] fTreeVarInt;
}
 
//___________________________________________________________________________
void AliAnalysisTaskCheckAODTracks::UserCreateOutputObjects() {
  // create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fTrackTree = new TTree("trackTree", "Tree for analysis");
  TString floatVarName[kNumOfFloatVar];
  fTreeVarFloat = new Float_t[kNumOfFloatVar];
  floatVarName[0]="xvert";
  floatVarName[1]="yvert";
  floatVarName[2]="zvert";
  floatVarName[3]="px";
  floatVarName[4]="py";
  floatVarName[5]="pz";
  floatVarName[6]="pt";
  floatVarName[7]="p";
  floatVarName[8]="eta";
  floatVarName[9]="phi";
  floatVarName[10]="d0xy";
  floatVarName[11]="d0z";
  floatVarName[12]="chi2clustpc";
  floatVarName[13]="chi2tpcconstglob";
  floatVarName[14]="croverfind";
  floatVarName[15]="dedxTPC";
  floatVarName[16]="nsigel";
  floatVarName[17]="nsigpi";
  floatVarName[18]="nsigk";
  floatVarName[19]="nsigp";
  floatVarName[20]="pxgen";
  floatVarName[21]="pygen";
  floatVarName[22]="pzgen";
  floatVarName[23]="ptgen";
  floatVarName[24]="pgen";
  floatVarName[25]="etagen";
  floatVarName[26]="phigen";
  Int_t usedVar=kNumOfFloatVar-7;
  if(fReadMC) usedVar=kNumOfFloatVar;
  for(Int_t ivar=0; ivar<usedVar; ivar++){
    fTrackTree->Branch(floatVarName[ivar].Data(),&fTreeVarFloat[ivar],Form("%s/F",floatVarName[ivar].Data()));
  }
  

  TString intVarName[kNumOfIntVar];
  fTreeVarInt = new Int_t[kNumOfIntVar];
  intVarName[0]="ntrack";
  intVarName[1]="ntracklets";
  intVarName[2]="trackcharge";
  intVarName[3]="ITSrefit";
  intVarName[4]="ITSclumap";
  intVarName[5]="nTPCclu";
  intVarName[6]="TOFbc";
  intVarName[7]="filbits";
  intVarName[8]="particlecharge";
  intVarName[9]="label";
  intVarName[10]="truePID";
  intVarName[11]="isPhysPrim";
  usedVar=kNumOfIntVar-4;
  if(fReadMC) usedVar=kNumOfIntVar;
  for(Int_t ivar=0; ivar<usedVar; ivar++){
    fTrackTree->Branch(intVarName[ivar].Data(),&fTreeVarInt[ivar],Form("%s/I",intVarName[ivar].Data()));
  }

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",15,-0.5,14.5);
  //fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"PhysSel");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Centrality");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Good vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Pass zSPD-zTrk vert sel");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"|zvert|<10");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Pileup cut");
  fOutput->Add(fHistNEvents);

  fHistNTracks = new TH1F("hNTracks", "Number of tracks in AOD events ; N_{tracks}",(Int_t)(fMaxMult+1.00001),-0.5,fMaxMult+0.5);
  fOutput->Add(fHistNTracks);
  fHistNTracksVsTPCclusters = new TH2F("hNTracksVsTPCclusters"," ; N_{TPCclusters} ; N_{tracks}",100,0.,300.*fMaxMult,100,0.,fMaxMult);
  fHistNTracksVsITSclusters = new TH2F("hNTracksVsITSclusters"," ; N_{ITSclusters} ; N_{tracks}",100,0.,10*fMaxMult,100,0.,fMaxMult);
  fHistITSclustersVsTPCclusters = new TH2F("hITSclustersVsTPCclusters"," ;  N_{ITSclusters} ; N_{TPCclusters}",100,0.,300.*fMaxMult,100,0.,10*fMaxMult);
  fHistNTracksFB4VsTPCclusters = new TH2F("hNTracksFB4VsTPCclusters"," ; N_{TPCclusters} ; N_{tracks,FB4}",100,0.,300.*fMaxMult,100,0.,fMaxMult);
  fHistNTracksFB4VsITSclusters = new TH2F("hNTracksFB4VsITSclusters"," ; N_{ITSclusters} ; N_{tracks,FB4}",100,0.,10*fMaxMult,100,0.,fMaxMult);
  fOutput->Add(fHistNTracksVsTPCclusters);
  fOutput->Add(fHistNTracksVsITSclusters);
  fOutput->Add(fHistITSclustersVsTPCclusters);
  fOutput->Add(fHistNTracksFB4VsTPCclusters);
  fOutput->Add(fHistNTracksFB4VsITSclusters);

  fHistFilterBits = new TH2D("hFilterBits", " ; Filter Bit ; Id ; N_{tracks}",kNumOfFilterBits,-0.5,kNumOfFilterBits-0.5,2,-1,1);
  fHistFilterBits->GetYaxis()->SetBinLabel(1,"Neg. ID");
  fHistFilterBits->GetYaxis()->SetBinLabel(2,"Pos. ID");
  fOutput->Add(fHistFilterBits);

  fHistITSnClusTPCsel=new TH1F("hITSnClusTPCsel","N ITS clusters (TPC+ITS)",7,-0.5,6.5);
  fHistITSnClusITSsa=new TH1F("hITSnClusITSsa","N ITS clusters (ITSsa)",7,-0.5,6.5);
  fHistITSnClusITSPureSA=new TH1F("hITSnClusITSPureSA","N ITS clusters (ITSpuresa)",7,-0.5,6.5);
  fOutput->Add(fHistITSnClusTPCsel);
  fOutput->Add(fHistITSnClusITSsa);
  fOutput->Add(fHistITSnClusITSPureSA);
  fHistITSCluInLayTPCsel = new TH1F("hITSCluInLayTPCsel", "N tracks with point on Layer (TPC+ITS); Layer; N tracks",6, -0.5, 5.5);
  fHistITSCluInLayITSsa  = new TH1F("hITSCluInLayITSsa", "N tracks with point on Layer (ITSsa); Layer; N tracks",6, -0.5, 5.5);
  fHistITSCluInLayITSPureSA  = new TH1F("hITSCluInLayITSPureSA", "N tracks with point on Layer (ITSpuresa); Layer; N tracks",6, -0.5, 5.5);
  fOutput->Add(fHistITSCluInLayTPCsel);
  fOutput->Add(fHistITSCluInLayITSsa);
  fOutput->Add(fHistITSCluInLayITSPureSA);

  fHistNtracksFb4VsV0befEvSel=new TH2F("hNtracksFb4VsV0befEvSel"," ; V0 signal ; N_{tracks,FilBit4}",250,0.,15000.,250,0.,5000.);
  fHistNtracksFb5VsV0befEvSel=new TH2F("hNtracksFb5VsV0befEvSel"," ; V0 signal ; N_{tracks,FilBit5}",250,0.,15000.,250,0.,5000.);
  fHistNtracksFb4VsV0aftEvSel=new TH2F("hNtracksFb4VsV0aftEvSel"," ; V0 signal ; N_{tracks,FilBit4}",250,0.,15000.,250,0.,5000.);
  fHistNtracksFb5VsV0aftEvSel=new TH2F("hNtracksFb5VsV0aftEvSel"," ; V0 signal ; N_{tracks,FilBit5}",250,0.,15000.,250,0.,5000.);
  fOutput->Add(fHistNtracksFb4VsV0befEvSel);
  fOutput->Add(fHistNtracksFb5VsV0befEvSel);
  fOutput->Add(fHistNtracksFb4VsV0aftEvSel);
  fOutput->Add(fHistNtracksFb5VsV0aftEvSel);

  Int_t nMultBins=(Int_t)(fMaxMult+1.00001);
  if(fMaxMult>1000.) nMultBins=500;
  fHistNtracksFb0=new TH1F("hNtracksFb0"," ; N_{tracks,FilBit0}",nMultBins,-0.5,fMaxMult+0.5);
  fHistNtracksFb1=new TH1F("hNtracksFb1"," ; N_{tracks,FilBit1}",nMultBins,-0.5,fMaxMult+0.5);
  fHistNtracksFb4=new TH1F("hNtracksFb4"," ; N_{tracks,FilBit4}",nMultBins,-0.5,fMaxMult+0.5);
  fHistNtracksFb5=new TH1F("hNtracksFb5"," ; N_{tracks,FilBit5}",nMultBins,-0.5,fMaxMult+0.5);
  fHistNtracksFb6=new TH1F("hNtracksFb6"," ; N_{tracks,FilBit6}",nMultBins,-0.5,fMaxMult+0.5);
  fHistNtracksFb7=new TH1F("hNtracksFb7"," ; N_{tracks,FilBit7}",nMultBins,-0.5,fMaxMult+0.5);
  fHistNtracksFb8=new TH1F("hNtracksFb8"," ; N_{tracks,FilBit8}",nMultBins,-0.5,fMaxMult+0.5);
  fOutput->Add(fHistNtracksFb0);
  fOutput->Add(fHistNtracksFb1);
  fOutput->Add(fHistNtracksFb4);
  fOutput->Add(fHistNtracksFb5);
  fOutput->Add(fHistNtracksFb6);
  fOutput->Add(fHistNtracksFb7);
  fOutput->Add(fHistNtracksFb8);

  fHistEtaPhiPtTPCsel = new TH3F("hEtaPhiPtTPCsel"," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtTPCselITSref = new TH3F("hEtaPhiPtTPCselITSref"," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtTPCselSPDany = new TH3F("hEtaPhiPtTPCselSPDany"," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fOutput->Add(fHistEtaPhiPtTPCsel);
  fOutput->Add(fHistEtaPhiPtTPCselITSref);
  fOutput->Add(fHistEtaPhiPtTPCselSPDany);

  fHistEtaPhiPtPosChargeTPCsel = new TH3F("hEtaPhiPtPosChargeTPCsel"," Positive charged tracks ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtPosChargeTPCselITSref = new TH3F("hEtaPhiPtPosChargeTPCselITSref"," Positive charged tracks ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtPosChargeTPCselSPDany = new TH3F("hEtaPhiPtPosChargeTPCselSPDany"," Positive charged tracks ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fOutput->Add(fHistEtaPhiPtPosChargeTPCsel);
  fOutput->Add(fHistEtaPhiPtPosChargeTPCselITSref);
  fOutput->Add(fHistEtaPhiPtPosChargeTPCselSPDany);

  fHistEtaPhiPtNegChargeTPCsel = new TH3F("hEtaPhiPtNegChargeTPCsel"," Negative charged tracks ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtNegChargeTPCselITSref = new TH3F("hEtaPhiPtNegChargeTPCselITSref"," Negative charged tracks ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtNegChargeTPCselSPDany = new TH3F("hEtaPhiPtNegChargeTPCselSPDany"," Negative charged tracks ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fOutput->Add(fHistEtaPhiPtNegChargeTPCsel);
  fOutput->Add(fHistEtaPhiPtNegChargeTPCselITSref);
  fOutput->Add(fHistEtaPhiPtNegChargeTPCselSPDany);

  fHistEtaPhiPtTPCselTOFbc = new TH3F("hEtaPhiPtTPCselTOFbc"," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtTPCselITSrefTOFbc = new TH3F("hEtaPhiPtTPCselITSrefTOFbc"," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtTPCselSPDanyTOFbc = new TH3F("hEtaPhiPtTPCselSPDanyTOFbc"," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtTPCselSPDanyTOFpid = new TH3F("hEtaPhiPtTPCselSPDanyTOFpid"," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fOutput->Add(fHistEtaPhiPtTPCselTOFbc);
  fOutput->Add(fHistEtaPhiPtTPCselITSrefTOFbc);
  fOutput->Add(fHistEtaPhiPtTPCselSPDanyTOFbc);
  fOutput->Add(fHistEtaPhiPtTPCselSPDanyTOFpid);

  fHistTPCchi2PerClusPhiPtTPCsel = new TH3F("hTPCchi2PerClusPhiPtTPCsel"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistTPCchi2PerClusPhiPtTPCselITSref = new TH3F("hTPCchi2PerClusPhiPtTPCselITSref"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistTPCchi2PerClusPhiPtTPCselSPDany = new TH3F("hTPCchi2PerClusPhiPtTPCselSPDany"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fOutput->Add(fHistTPCchi2PerClusPhiPtTPCsel);
  fOutput->Add(fHistTPCchi2PerClusPhiPtTPCselITSref);
  fOutput->Add(fHistTPCchi2PerClusPhiPtTPCselSPDany);
  
  fHistSig1ptCovMatPhiPtTPCsel = new TH3F("hSig1ptCovMatPhiPtTPCsel"," ; p_{T}*#sigma(1/p_{T}); p_{T} (GeV/c) ; #varphi",100, 0, 0.3, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistSig1ptCovMatPhiPtTPCselITSref = new TH3F("hSig1ptCovMatPhiPtTPCselITSref"," ; p_{T}*#sigma(1/p_{T}); p_{T} (GeV/c) ; #varphi",100, 0, 0.3, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistSig1ptCovMatPhiPtTPCselSPDany = new TH3F("hSig1ptCovMatPhiPtTPCselSPDany"," ; p_{T}*#sigma(1/p_{T}); p_{T} (GeV/c) ; #varphi",100, 0, 0.3, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fOutput->Add(fHistSig1ptCovMatPhiPtTPCsel);
  fOutput->Add(fHistSig1ptCovMatPhiPtTPCselITSref);
  fOutput->Add(fHistSig1ptCovMatPhiPtTPCselSPDany);

  const Int_t nPtBins4ip=90;
  Double_t ptBins4ip[nPtBins4ip+1];
  for(Int_t jjj=0; jjj<=50; jjj++) ptBins4ip[jjj]=0.020*jjj;
  for(Int_t jjj=51; jjj<=70; jjj++) ptBins4ip[jjj]=ptBins4ip[50]+0.050*(jjj-50);
  for(Int_t jjj=71; jjj<=80; jjj++) ptBins4ip[jjj]=ptBins4ip[70]+0.1*(jjj-70);
  for(Int_t jjj=81; jjj<=84; jjj++) ptBins4ip[jjj]=ptBins4ip[80]+0.5*(jjj-80);
  for(Int_t jjj=85; jjj<=90; jjj++) ptBins4ip[jjj]=ptBins4ip[84]+1.*(jjj-84);
  const Int_t nMultBins4ip=12;
  Double_t multBins4ip[nMultBins4ip+1]={0.,20.,50.,100.,500.,1000.,
					1500.,2000.,3000.,4000.,5000.,7500.,10000.};
  const Int_t nIPBins4ip=400;
  Double_t ipBins4ip[nIPBins4ip+1];
  for(Int_t jjj=0; jjj<=nIPBins4ip; jjj++) ipBins4ip[jjj]=-1500.+(3000./(Double_t)nIPBins4ip)*(Double_t)jjj;

  fHistImpParXYPtMulPionTPCselSPDany = new TH3F("hImpParXYPtMulPionTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fHistImpParXYPtMulKaonTPCselSPDany = new TH3F("hImpParXYPtMulKaonTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fHistImpParXYPtMulProtonTPCselSPDany = new TH3F("hImpParXYPtMulProtonTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fOutput->Add(fHistImpParXYPtMulPionTPCselSPDany);
  fOutput->Add(fHistImpParXYPtMulKaonTPCselSPDany);
  fOutput->Add(fHistImpParXYPtMulProtonTPCselSPDany);
  for(Int_t jb=0; jb<kNumOfFilterBits; jb++){
    fHistEtaPhiPtFiltBit[jb] =  new TH3F(Form("hEtaPhiPtFiltBit%d",jb)," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistImpParXYPtMulFiltBit[jb] = new TH3F(Form("hImpParXYPtMulPionFiltBit%d",jb)," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
    fHistITScluPtFiltBit[jb] = new TH2F(Form("hITScluPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; n ITS clusters",50,0.,10.,7,-0.5,6.5);
    fHistSPDcluPtFiltBit[jb] = new TH2F(Form("hSPDcluPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; SPD clusters",50,0.,10.,4,-0.5,3.5);
    fHistSPDcluPtFiltBit[jb]->GetYaxis()->SetBinLabel(1,"kNone");
    fHistSPDcluPtFiltBit[jb]->GetYaxis()->SetBinLabel(2,"kOnlyFirst");
    fHistSPDcluPtFiltBit[jb]->GetYaxis()->SetBinLabel(3,"kOnlySecond");
    fHistSPDcluPtFiltBit[jb]->GetYaxis()->SetBinLabel(4,"kBoth");
    fHistTPCcluPtFiltBit[jb] = new TH2F(Form("hTPCcluPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; n TPC clusters",50,0.,10.,161,-0.5,160.5);  
    fHistTPCcrrowsPtFiltBit[jb] = new TH2F(Form("hTPCcrrowsPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; n TPC Crossed Rows",50,0.,10.,161,-0.5,160.5);  
    fHistTPCCrowOverFindPtFiltBit[jb] = new TH2F(Form("hTPCCrowOverFindPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; nTPCCrossedRows / nTPCFindableClusters",50,0.,10.,100,0.,2.);    
    fHistTPCChi2clusPtFiltBit[jb] = new TH2F(Form("hTPCChi2clusPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; #chi^{2}/nTPCclusters",50,0.,10.,160,0.,8.);
    fHistChi2TPCConstrVsGlobPtFiltBit[jb] = new TH2F(Form("hChi2TPCConstrVsGlobPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; golden #chi^{2}",50,0.,10.,160,0.,40.);
    fHistSig1ptCovMatPtFiltBit[jb] = new TH2F(Form("hSig1ptCovMatPtFiltBit%d",jb)," ;  p_{T} (GeV/c) ; p_{T}*#sigma(1/p_{T})",50,0.,10.,100, 0, 0.3);
    fOutput->Add(fHistEtaPhiPtFiltBit[jb]);
    fOutput->Add(fHistImpParXYPtMulFiltBit[jb]);
    fOutput->Add(fHistITScluPtFiltBit[jb]);
    fOutput->Add(fHistSPDcluPtFiltBit[jb]);
    fOutput->Add(fHistTPCcluPtFiltBit[jb]);
    fOutput->Add(fHistTPCcrrowsPtFiltBit[jb]);
    fOutput->Add(fHistTPCChi2clusPtFiltBit[jb]);
    fOutput->Add(fHistChi2TPCConstrVsGlobPtFiltBit[jb]);
    fOutput->Add(fHistTPCCrowOverFindPtFiltBit[jb]);
    fOutput->Add(fHistSig1ptCovMatPtFiltBit[jb]);
  }

  TString xTit="p_{T,rec} (GeV/c)";
  if(fUseGenPt) xTit="p_{T,gen} (GeV/c)";
  fHistPtResidVsPtTPCselAll = new TH2F("hPtResidVsPtTPCselAll",Form("All particles ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefAll = new TH2F("hPtResidVsPtTPCselITSrefAll",Form("All particles ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselAll = new TH2F("hOneOverPtResidVsPtTPCselAll",Form("All particles ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselITSrefAll = new TH2F("hOneOverPtResidVsPtTPCselITSrefAll",Form("All particles ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistPtResidVsPtTPCselPrim = new TH2F("hPtResidVsPtTPCselPrim",Form("Primaries ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefPrim = new TH2F("hPtResidVsPtTPCselITSrefPrim",Form("Primaries ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselPrim = new TH2F("hOneOverPtResidVsPtTPCselPrim",Form("Primaries ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselITSrefPrim = new TH2F("hOneOverPtResidVsPtTPCselITSrefPrim",Form("Primaries ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistPtResidVsPtTPCselSecDec = new TH2F("hPtResidVsPtTPCselSecDec",Form("Secondaries from decays ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefSecDec = new TH2F("hPtResidVsPtTPCselITSrefSecDec",Form("Secondaries from decays ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselSecDec = new TH2F("hOneOverPtResidVsPtTPCselSecDec",Form("Secondaries from decays ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselITSrefSecDec = new TH2F("hOneOverPtResidVsPtTPCselITSrefSecDec",Form("Secondaries from decays ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistPtResidVsPtTPCselSecMat = new TH2F("hPtResidVsPtTPCselSecMat",Form("Secondaries from material ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefSecMat = new TH2F("hPtResidVsPtTPCselITSrefSecMat",Form("Secondaries from material ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselSecMat = new TH2F("hOneOverPtResidVsPtTPCselSecMat",Form("Secondaries from material ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselITSrefSecMat = new TH2F("hOneOverPtResidVsPtTPCselITSrefSecMat",Form("Secondaries from material ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",xTit.Data()),fNPtBins,fMinPt,fMaxPt,200,-0.5,0.5);
  fOutput->Add(fHistPtResidVsPtTPCselAll);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefAll);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselAll);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselITSrefAll);
  fOutput->Add(fHistPtResidVsPtTPCselPrim);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefPrim);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselPrim);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselITSrefPrim);
  fOutput->Add(fHistPtResidVsPtTPCselSecDec);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefSecDec);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselSecDec);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselITSrefSecDec);
  fOutput->Add(fHistPtResidVsPtTPCselSecMat);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefSecMat);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselSecMat);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselITSrefSecMat);

  for (int iS = 0; iS < AliPID::kSPECIESC; ++iS) {
    fHistPtResidVsPtTPCsel[iS] = new TH2F(Form("hPtResidVsPtTPCsel%s",AliPID::ParticleShortName(iS)), Form("%s ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",AliPID::ParticleLatexName(iS),xTit.Data()),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistPtResidVsPtTPCselITSref[iS] = new TH2F(Form("hPtResidVsPtTPCselITSref%s",AliPID::ParticleShortName(iS)),Form("%s ; %s ; p_{T,rec}-p_{T,gen} (GeV/c)",AliPID::ParticleLatexName(iS),xTit.Data()),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistOneOverPtResidVsPtTPCsel[iS] = new TH2F(Form("hOneOverPtResidVsPtTPCsel%s",AliPID::ParticleShortName(iS)), Form("%s ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",AliPID::ParticleLatexName(iS),xTit.Data()),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistOneOverPtResidVsPtTPCselITSref[iS] = new TH2F(Form("hOneOverPtResidVsPtTPCselITSref%s",AliPID::ParticleShortName(iS)),Form("%s ; %s ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",AliPID::ParticleLatexName(iS),xTit.Data()),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fOutput->Add(fHistPtResidVsPtTPCsel[iS]);
    fOutput->Add(fHistPtResidVsPtTPCselITSref[iS]);
    fOutput->Add(fHistOneOverPtResidVsPtTPCsel[iS]);
    fOutput->Add(fHistOneOverPtResidVsPtTPCselITSref[iS]);
  }
  fHistPzResidVsPtTPCselAll = new TH2F("hPzResidVsPtTPCselAll",Form("; %s ; p_{z,rec}-p_{z,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,150,-0.3,0.3);
  fHistPzResidVsPtTPCselITSrefAll = new TH2F("hPzResidVsPtTPCselITSrefAll",Form("; %s ; p_{z,rec}-p_{z,gen} (GeV/c)",xTit.Data()),fNPtBins,fMinPt,fMaxPt,150,-0.3,0.3);
  fHistPzResidVsEtaTPCselAll = new TH2F("hPzResidVsEtaTPCselAll","; #eta ; p_{z,rec}-p_{z,gen} (GeV/c)",fNEtaBins,-1.,1.,150,-0.3,0.3);
  fHistPzResidVsEtaTPCselITSrefAll = new TH2F("hPzResidVsEtaTPCselITSrefAll","; #eta ; p_{z,rec}-p_{z,gen} (GeV/c)",fNEtaBins,-1.,1.,150,-0.3,0.3);
  fOutput->Add(fHistPzResidVsPtTPCselAll);
  fOutput->Add(fHistPzResidVsPtTPCselITSrefAll);
  fOutput->Add(fHistPzResidVsEtaTPCselAll);
  fOutput->Add(fHistPzResidVsEtaTPCselITSrefAll);
  fHistEtaPhiPtTPCselITSrefGood = new TH3F("hEtaPhiPtTPCselITSrefGood"," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtTPCselITSrefFake = new TH3F("hEtaPhiPtTPCselITSrefFake"," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fOutput->Add(fHistEtaPhiPtTPCselITSrefGood);
  fOutput->Add(fHistEtaPhiPtTPCselITSrefFake);

  fHistImpParXYPtMulTPCselSPDanyGood = new TH3F("hImpParXYPtMulTPCselSPDanyGood",Form(" ; %s ; d_{0}^{xy} (#mum) ; N_{CL1}",xTit.Data()),nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fHistImpParXYPtMulTPCselSPDanyFake = new TH3F("hImpParXYPtMulTPCselSPDanyFake",Form(" ; %s ; d_{0}^{xy} (#mum) ; N_{CL1}",xTit.Data()),nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyGood);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyFake);

  fHistImpParXYPtMulTPCselSPDanyPrim = new TH3F("hImpParXYPtMulTPCselSPDanyPrim",Form(" ; %s ; d_{0}^{xy} (#mum) ; N_{CL1}",xTit.Data()),nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fHistImpParXYPtMulTPCselSPDanySecDec = new TH3F("hImpParXYPtMulTPCselSPDanySecDec",Form(" ; %s ; d_{0}^{xy} (#mum) ; N_{CL1}",xTit.Data()),nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fHistImpParXYPtMulTPCselSPDanySecMat = new TH3F("hImpParXYPtMulTPCselSPDanySecMat",Form(" ; %s ; d_{0}^{xy} (#mum) ; N_{CL1}",xTit.Data()),nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyPrim);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanySecDec);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanySecMat);


  fHistInvMassK0s = new TH3F("hInvMassK0s"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(K0s) ; R (cm)",200,0.4,0.6,50,0.,10.,50,0.,50.);
  fHistInvMassLambda = new TH3F("hInvMassLambda"," ;Inv.Mass (GeV/c^{2}) ; p_{T}(#Lambda) ; R (cm)",200,1.0,1.2,50,0.,10.,50,0.,50.);
  fHistInvMassAntiLambda = new TH3F("hInvMassAntiLambda"," ;Inv.Mass (GeV/c^{2}) ; p_{T}(#bar{#Lambda}) ; R (cm)",200,1.0,1.2,50,0.,10.,50,0.,50.);
  fOutput->Add(fHistInvMassK0s);
  fOutput->Add(fHistInvMassLambda);
  fOutput->Add(fHistInvMassAntiLambda);

  PostData(1,fOutput);
  PostData(2,fTrackTree);

}
//______________________________________________________________________________
void AliAnalysisTaskCheckAODTracks::UserExec(Option_t *)
{
  //

  AliAODEvent *aod = (AliAODEvent*) (InputEvent());
  if(!aod) {
    printf("AliAnalysisTaskCheckAODTracks::UserExec(): bad AOD\n");
    return;
  } 
  Double_t magField  = aod->GetMagneticField();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
  
  TClonesArray *arrayMC=0;

  if(fReadMC){
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      Printf("ERROR: MC particles branch not found!\n");
      return;
    }
  }

  
  fHistNEvents->Fill(0);
  if(fUsePhysSel){
    Bool_t isPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
    if(!isPhysSel) return;
  }
  fHistNEvents->Fill(1);

  if(fSelectOnCentrality){
    AliMultSelection* mulSel = (AliMultSelection*)aod->FindListObject("MultSelection");
    if(mulSel){
      Double_t centr=mulSel->GetMultiplicityPercentile(fCentrEstimator.Data());
      if(centr<fMinCentrality || centr>fMaxCentrality) return;
    }
  }
  fHistNEvents->Fill(2);

  Int_t ntracks = aod->GetNumberOfTracks();
  Int_t ntracklets = 0;
  AliAODTracklets *mult=aod->GetTracklets();
  if(mult) ntracklets=mult->GetNumberOfTracklets();
  Int_t ncl1 = aod->GetNumberOfITSClusters(1);
  Int_t totITSclusters=0;
  for(Int_t il=0; il<6; il++) totITSclusters+=aod->GetNumberOfITSClusters(il);
  Int_t totTPCclusters=aod->GetNumberOfTPCClusters();
  Int_t ntracksFB0=0;
  Int_t ntracksFB1=0;
  Int_t ntracksFB4=0;
  Int_t ntracksFB5=0;
  Int_t ntracksFB6=0;
  Int_t ntracksFB7=0;
  Int_t ntracksFB8=0;
  for (Int_t iTrack=0; iTrack < ntracks; iTrack++) {
    AliAODTrack * track = (AliAODTrack*)aod->GetTrack(iTrack);
    if (!track) continue;
    if(track->TestFilterBit(1)) ntracksFB0++;
    if(track->TestFilterBit(1<<1)) ntracksFB1++;
    if(track->TestFilterBit(1<<4)) ntracksFB4++;
    if(track->TestFilterBit(1<<5)) ntracksFB5++;
    if(track->TestFilterBit(1<<6)) ntracksFB6++;
    if(track->TestFilterBit(1<<7)) ntracksFB7++;
    if(track->TestFilterBit(1<<8)) ntracksFB8++;
  }
  Double_t vZEROampl=0;
  for(Int_t i=0;i<64;i++) vZEROampl+=aod->GetVZEROData()->GetMultiplicity(i);
  fHistNtracksFb4VsV0befEvSel->Fill(vZEROampl,ntracksFB4);
  fHistNtracksFb5VsV0befEvSel->Fill(vZEROampl,ntracksFB5);

  const AliVVertex* vtTrc = aod->GetPrimaryVertex();
  const AliVVertex* vtSPD = aod->GetPrimaryVertexSPD();
  TString titTrc=vtTrc->GetTitle();
  if(titTrc.IsNull() || titTrc=="vertexer: 3D" || titTrc=="vertexer: Z") return;
  if (vtSPD->GetNContributors()<1) return;
  fHistNEvents->Fill(3);

  double covTrc[6],covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = vtTrc->GetZ()-vtSPD->GetZ();
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
  if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return; // bad vertexing
  fHistNEvents->Fill(4);

  Float_t xvert=vtTrc->GetX();
  Float_t yvert=vtTrc->GetY();
  Float_t zvert=vtTrc->GetZ();
  if(TMath::Abs(zvert)>10) return;
  fHistNEvents->Fill(5);

  if(fUsePileupCut){
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(5);
    utils.SetMaxPlpChi2MV(5.);
    utils.SetMinWDistMV(15.);
    utils.SetCheckPlpFromDifferentBCMV(kTRUE);
    Bool_t isPUMV = utils.IsPileUpMV(aod);
    if(isPUMV) return;
    fHistNEvents->Fill(6);
  }

  fHistNtracksFb4VsV0aftEvSel->Fill(vZEROampl,ntracksFB4);
  fHistNtracksFb5VsV0aftEvSel->Fill(vZEROampl,ntracksFB5);

  Bool_t iskINT7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  if(iskINT7){
    fHistNtracksFb0->Fill(ntracksFB0);
    fHistNtracksFb1->Fill(ntracksFB1);
    fHistNtracksFb4->Fill(ntracksFB4);
    fHistNtracksFb5->Fill(ntracksFB5);
    fHistNtracksFb6->Fill(ntracksFB6);
    fHistNtracksFb7->Fill(ntracksFB7);
    fHistNtracksFb8->Fill(ntracksFB8);
  }
  
  Double_t pos[3],cov[6];
  vtTrc->GetXYZ(pos);
  vtTrc->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);

  fHistNTracks->Fill(ntracks);
  fHistNTracksVsTPCclusters->Fill(totTPCclusters,ntracks);
  fHistNTracksVsITSclusters->Fill(totITSclusters,ntracks);
  fHistITSclustersVsTPCclusters->Fill(totTPCclusters,totITSclusters);
  fHistNTracksFB4VsTPCclusters->Fill(totTPCclusters,ntracksFB4);
  fHistNTracksFB4VsITSclusters->Fill(totITSclusters,ntracksFB4);

  for (Int_t iTrack=0; iTrack < ntracks; iTrack++) {
    AliAODTrack * track = (AliAODTrack*)aod->GetTrack(iTrack);
    if (!track) continue;
    for(Int_t jvar=0; jvar<kNumOfFloatVar; jvar++) fTreeVarFloat[jvar]=-999.;
    for(Int_t jvar=0; jvar<kNumOfIntVar; jvar++) fTreeVarInt[jvar]=-999;

    fTreeVarFloat[0]=xvert;
    fTreeVarFloat[1]=yvert;
    fTreeVarFloat[2]=zvert;
    fTreeVarInt[0]=ntracks;
    fTreeVarInt[1]=ntracklets;

    Int_t trid=track->GetID();
    Double_t ydum=0.5;
    if(trid<0) ydum=-0.5;
    for(Int_t jb=0; jb<kNumOfFilterBits; jb++){
      if(track->TestFilterBit(1<<jb)) fHistFilterBits->Fill(jb,ydum);
    }

    Int_t chtrack=track->Charge();
    Double_t pttrack=track->Pt();
    Double_t etatrack=track->Eta();
    Double_t phitrack=track->Phi();
    fTreeVarFloat[3]=track->Px();
    fTreeVarFloat[4]=track->Py();
    fTreeVarFloat[5]=track->Pz();
    fTreeVarFloat[6]=track->Pt();
    fTreeVarFloat[7]=track->P();
    fTreeVarFloat[8]=track->Eta();
    fTreeVarFloat[9]=track->Phi();
    
    Double_t d0z0[2],covd0z0[3];
    Bool_t isOK=track->PropagateToDCA(vtTrc,magField,99999.,d0z0,covd0z0);
    Float_t impactXY=-999, impactZ=-999;
    if(isOK){
      impactXY=d0z0[0];
      impactZ=d0z0[1];
    }
    Bool_t itsRefit=kFALSE;
    Int_t statusTrack=track->GetStatus();
    if(statusTrack&AliESDtrack::kITSrefit) itsRefit=kTRUE; 
    Int_t nITSclus=track->GetNcls(0);
    UChar_t clumap=track->GetITSClusterMap();
    if(itsRefit && !(statusTrack&AliESDtrack::kTPCin)){
      if(statusTrack&AliESDtrack::kITSpureSA){
	fHistITSnClusITSPureSA->Fill(nITSclus);
	for(Int_t layer=0; layer<6; layer++) {
	  if(TESTBIT(clumap,layer)) fHistITSCluInLayITSPureSA->Fill(layer);
	}
      }else{
	fHistITSnClusITSsa->Fill(nITSclus);
	for(Int_t layer=0; layer<6; layer++) {
	  if(TESTBIT(clumap,layer)) fHistITSCluInLayITSsa->Fill(layer);
	}
      }
    }
    Int_t nSPDclus=0;
    if(track->HasPointOnITSLayer(0)) nSPDclus+=1;
    if(track->HasPointOnITSLayer(1)) nSPDclus+=2;
    Bool_t spdAny=kFALSE;
    if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) spdAny=kTRUE;
    Int_t nTPCclus=track->GetNcls(1);
    Double_t chi2ndf=track->Chi2perNDF();
    Double_t chi2tpc=999.;
    if(chi2ndf>0. && nTPCclus > 5){
      chi2tpc=Float_t(nTPCclus-5)*chi2ndf;
    }
    Double_t chi2clus=-1.;
    if(nTPCclus>0) chi2clus=chi2tpc/(Float_t)nTPCclus;
    Double_t goldenChi2=track->GetChi2TPCConstrainedVsGlobal();
    Float_t nCrossedRowsTPC = track->GetTPCCrossedRows();
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (track->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
    }
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(track);
    Double_t curvrelerr = TMath::Sqrt(etp.GetSigma1Pt2())/etp.OneOverPt();

    Int_t tofBC=track->GetTOFBunchCrossing(magField);
    fTreeVarFloat[10]=impactXY;
    fTreeVarFloat[11]=impactZ;
    fTreeVarFloat[12]=chi2clus;
    fTreeVarFloat[13]=goldenChi2;
    fTreeVarFloat[14]=ratioCrossedRowsOverFindableClustersTPC;

    fTreeVarInt[2]=chtrack;
    fTreeVarInt[3]=itsRefit;
    fTreeVarInt[4]=clumap;
    fTreeVarInt[5]=nTPCclus;
    fTreeVarInt[6]=tofBC;

    Int_t trlabel=track->GetLabel();
    Float_t dedx=track->GetTPCsignal();
    Int_t  filtmap=track->GetFilterMap();
    
    Double_t nSigmaTPC[AliPID::kSPECIESC];
    for(Int_t jsp=0; jsp<AliPID::kSPECIESC; jsp++) nSigmaTPC[jsp]=-999.;
    if(pidResp){
      AliPIDResponse::EDetPidStatus status = pidResp->CheckPIDStatus(AliPIDResponse::kTPC,track);
      if (status == AliPIDResponse::kDetPidOk){
	for(Int_t jsp=0; jsp<AliPID::kSPECIESC; jsp++){
	  nSigmaTPC[jsp]=pidResp->NumberOfSigmasTPC(track,(AliPID::EParticleType)jsp);
	}
      }
    }
    fTreeVarFloat[15]=dedx;
    fTreeVarFloat[16]=nSigmaTPC[0];
    fTreeVarFloat[17]=nSigmaTPC[2];
    fTreeVarFloat[18]=nSigmaTPC[3];
    fTreeVarFloat[19]=nSigmaTPC[4];
    fTreeVarInt[7]=filtmap;
    
    Float_t ptgen=-999.;
    Float_t pgen=-999.;
    Float_t pxgen=-999.;
    Float_t pygen=-999.;
    Float_t pzgen=-999.;
    Float_t etagen=-999.;
    Float_t phigen=-999.;
    Int_t hadronSpecies=-1;
    Float_t invptgen=-999.;
    Int_t isPhysPrim=-999;
    if(fReadMC){
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(TMath::Abs(trlabel)));
      ptgen=part->Pt();
      pgen=part->P();
      pxgen=part->Px();
      pygen=part->Py();
      pzgen=part->Pz();
      if(ptgen>0.) invptgen=1./ptgen;
      etagen=part->Eta();
      phigen=part->Phi();
      if(part->IsPhysicalPrimary()) isPhysPrim=1;
      else if(part->IsSecondaryFromWeakDecay()) isPhysPrim=0;
      else if(part->IsSecondaryFromMaterial()) isPhysPrim=-1;
      fTreeVarFloat[20]=pxgen;
      fTreeVarFloat[21]=pygen;
      fTreeVarFloat[22]=pzgen;
      fTreeVarFloat[23]=ptgen;
      fTreeVarFloat[24]=pgen;
      fTreeVarFloat[25]=etagen;
      fTreeVarFloat[26]=phigen;
      if (fUseMCId) {
        int pdg = TMath::Abs(part->GetPdgCode());
        for (int iS = 0; iS < AliPID::kSPECIESC; ++iS) {
          if (pdg == AliPID::ParticleCode(iS)) hadronSpecies=iS;
        }
      }
      fTreeVarInt[8]=part->Charge();
      fTreeVarInt[9]=trlabel;
      fTreeVarInt[10]=part->GetPdgCode();
      fTreeVarInt[11]=isPhysPrim;
    }

    if (fFillTree) fTrackTree->Fill();    

    for(Int_t jb=0; jb<kNumOfFilterBits; jb++){
      if(track->TestFilterBit(1<<jb)){
	fHistImpParXYPtMulFiltBit[jb]->Fill(pttrack,impactXY*10000.,ncl1);
	fHistEtaPhiPtFiltBit[jb]->Fill(etatrack,phitrack,pttrack);
	fHistITScluPtFiltBit[jb]->Fill(pttrack,nITSclus);
 	fHistSPDcluPtFiltBit[jb]->Fill(pttrack,nSPDclus);
 	fHistTPCcluPtFiltBit[jb]->Fill(pttrack,nTPCclus);
	fHistTPCcrrowsPtFiltBit[jb]->Fill(pttrack,nCrossedRowsTPC);
	fHistTPCCrowOverFindPtFiltBit[jb]->Fill(pttrack,ratioCrossedRowsOverFindableClustersTPC);
	fHistTPCChi2clusPtFiltBit[jb]->Fill(pttrack,chi2clus);
	fHistChi2TPCConstrVsGlobPtFiltBit[jb]->Fill(pttrack,goldenChi2);
	fHistSig1ptCovMatPtFiltBit[jb]->Fill(pttrack,curvrelerr);
      }
    }

    if(track->GetID()<0) continue;
    // convert to ESD track here
    if(ConvertAndSelectAODTrack(track,vESD,magField)==kFALSE) continue;

    fHistITSnClusTPCsel->Fill(nITSclus);
    for(Int_t layer=0; layer<6; layer++) {
      if(TESTBIT(clumap,layer)) fHistITSCluInLayTPCsel->Fill(layer);
    }

    fHistEtaPhiPtTPCsel->Fill(etatrack,phitrack,pttrack);
    if(chtrack>0) fHistEtaPhiPtPosChargeTPCsel->Fill(etatrack,phitrack,pttrack);
    else if(chtrack<0) fHistEtaPhiPtNegChargeTPCsel->Fill(etatrack,phitrack,pttrack);
    if(tofBC==0) fHistEtaPhiPtTPCselTOFbc->Fill(etatrack,phitrack,pttrack);
    if(itsRefit){
      fHistEtaPhiPtTPCselITSref->Fill(etatrack,phitrack,pttrack);
      if(chtrack>0) fHistEtaPhiPtPosChargeTPCselITSref->Fill(etatrack,phitrack,pttrack);
      else if(chtrack<0) fHistEtaPhiPtNegChargeTPCselITSref->Fill(etatrack,phitrack,pttrack);
      if(tofBC==0) fHistEtaPhiPtTPCselITSrefTOFbc->Fill(etatrack,phitrack,pttrack);
      if(spdAny){ 
        fHistEtaPhiPtTPCselSPDany->Fill(etatrack,phitrack,pttrack);
        if(chtrack>0) fHistEtaPhiPtPosChargeTPCselSPDany->Fill(etatrack,phitrack,pttrack);
        else if(chtrack<0) fHistEtaPhiPtNegChargeTPCselSPDany->Fill(etatrack,phitrack,pttrack);
        if(tofBC==0) fHistEtaPhiPtTPCselSPDanyTOFbc->Fill(etatrack,phitrack,pttrack);
        if ((track->GetStatus() & AliVTrack::kTOFout) &&
            (track->GetStatus() & AliVTrack::kTIME) &&
            (track->GetIntegratedLength() > 350.)) {
          fHistEtaPhiPtTPCselSPDanyTOFpid->Fill(etatrack,phitrack,pttrack);
        }
      }
    }

    fHistTPCchi2PerClusPhiPtTPCsel->Fill(chi2clus,pttrack,phitrack);
    fHistSig1ptCovMatPhiPtTPCsel->Fill(curvrelerr,pttrack,phitrack);
    if(itsRefit){
      fHistTPCchi2PerClusPhiPtTPCselITSref->Fill(chi2clus,pttrack,phitrack);
      fHistSig1ptCovMatPhiPtTPCselITSref->Fill(curvrelerr,pttrack,phitrack);
      if(spdAny){
	fHistTPCchi2PerClusPhiPtTPCselSPDany->Fill(chi2clus,pttrack,phitrack);
	fHistSig1ptCovMatPhiPtTPCselSPDany->Fill(curvrelerr,pttrack,phitrack);
      }
    }

    bool pid[AliPID::kSPECIESC] = {false};
    if (fReadMC && fUseMCId) {
      if (hadronSpecies > -1) pid[hadronSpecies] = true;
    } else {
      for (int iS = 0; iS < AliPID::kSPECIESC; ++iS)
        pid[iS] = TMath::Abs(nSigmaTPC[iS])<3;
    }
    bool isProton = pid[AliPID::kProton];
    bool isKaon = pid[AliPID::kKaon];
    bool isPion = pid[AliPID::kPion];

    if(itsRefit && spdAny){
      if(isPion) fHistImpParXYPtMulPionTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1); 
      if(isKaon) fHistImpParXYPtMulKaonTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1); 
      if(isProton) fHistImpParXYPtMulProtonTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1);
    }

    if(fReadMC && pttrack>0.){
      Double_t ptOnX=pttrack;
      if(fUseGenPt) ptOnX=ptgen;
      fHistPtResidVsPtTPCselAll->Fill(ptOnX,(pttrack-ptgen));
      fHistOneOverPtResidVsPtTPCselAll->Fill(ptOnX,pttrack*(1./pttrack-invptgen));
      fHistPzResidVsPtTPCselAll->Fill(ptOnX,(track->Pz()-pzgen));
      fHistPzResidVsEtaTPCselAll->Fill(etatrack,(track->Pz()-pzgen));
      if (itsRefit){
	fHistPtResidVsPtTPCselITSrefAll->Fill(ptOnX,(pttrack-ptgen));
	fHistOneOverPtResidVsPtTPCselITSrefAll->Fill(ptOnX,pttrack*(1./pttrack-invptgen));
	fHistPzResidVsPtTPCselITSrefAll->Fill(ptOnX,(track->Pz()-pzgen));
	fHistPzResidVsEtaTPCselITSrefAll->Fill(etatrack,(track->Pz()-pzgen));
      }
      if(isPhysPrim==1){
	fHistPtResidVsPtTPCselPrim->Fill(ptOnX,(pttrack-ptgen));
	fHistOneOverPtResidVsPtTPCselPrim->Fill(ptOnX,pttrack*(1./pttrack-invptgen));
	if (itsRefit){
	  fHistPtResidVsPtTPCselITSrefPrim->Fill(ptOnX,(pttrack-ptgen));
	  fHistOneOverPtResidVsPtTPCselITSrefPrim->Fill(ptOnX,pttrack*(1./pttrack-invptgen));
	}
      }else if(isPhysPrim==0){
	fHistPtResidVsPtTPCselSecDec->Fill(ptOnX,(pttrack-ptgen));
	fHistOneOverPtResidVsPtTPCselSecDec->Fill(ptOnX,pttrack*(1./pttrack-invptgen));
	if (itsRefit){
	  fHistPtResidVsPtTPCselITSrefSecDec->Fill(ptOnX,(pttrack-ptgen));
	  fHistOneOverPtResidVsPtTPCselITSrefSecDec->Fill(ptOnX,pttrack*(1./pttrack-invptgen));
	}
      }else if(isPhysPrim==-1){
	fHistPtResidVsPtTPCselSecMat->Fill(ptOnX,(pttrack-ptgen));
	fHistOneOverPtResidVsPtTPCselSecMat->Fill(ptOnX,pttrack*(1./pttrack-invptgen));
	if (itsRefit){
	  fHistPtResidVsPtTPCselITSrefSecMat->Fill(ptOnX,(pttrack-ptgen));
	  fHistOneOverPtResidVsPtTPCselITSrefSecMat->Fill(ptOnX,pttrack*(1./pttrack-invptgen));
	}
      }
      if(isPhysPrim==1){
	for (int iS = 0; iS < AliPID::kSPECIESC; ++iS) {
	  if (pid[iS]) {
	    Double_t ptDiff=pttrack*AliPID::ParticleCharge(iS)-ptgen;
	    Double_t oneOverPtDiff=0;
	    if(AliPID::ParticleCharge(iS)>0) oneOverPtDiff=pttrack*AliPID::ParticleCharge(iS)*(1./(pttrack*AliPID::ParticleCharge(iS))-invptgen);
	    fHistPtResidVsPtTPCsel[iS]->Fill(ptOnX*AliPID::ParticleCharge(iS),ptDiff);
	    fHistOneOverPtResidVsPtTPCsel[iS]->Fill(ptOnX*AliPID::ParticleCharge(iS),oneOverPtDiff);
	    if (itsRefit){ 
	      fHistPtResidVsPtTPCselITSref[iS]->Fill(ptOnX*AliPID::ParticleCharge(iS),ptDiff);
	      fHistOneOverPtResidVsPtTPCselITSref[iS]->Fill(ptOnX*AliPID::ParticleCharge(iS),oneOverPtDiff);
	    }
	  }
	}
      }
    
      if(trlabel>=0){
	fHistEtaPhiPtTPCselITSrefGood->Fill(etatrack,phitrack,pttrack);
	if(itsRefit && spdAny) fHistImpParXYPtMulTPCselSPDanyGood->Fill(ptOnX,impactXY*10000.,ncl1);
      }else{
	fHistEtaPhiPtTPCselITSrefFake->Fill(etatrack,phitrack,pttrack);
	if(itsRefit && spdAny) fHistImpParXYPtMulTPCselSPDanyFake->Fill(ptOnX,impactXY*10000.,ncl1);
      }
      
      if(itsRefit && spdAny){
	if(isPhysPrim==1) fHistImpParXYPtMulTPCselSPDanyPrim->Fill(ptOnX,impactXY*10000.,ncl1);
	else if(isPhysPrim==0) fHistImpParXYPtMulTPCselSPDanySecDec->Fill(ptOnX,impactXY*10000.,ncl1);
	else if(isPhysPrim==-1) fHistImpParXYPtMulTPCselSPDanySecMat->Fill(ptOnX,impactXY*10000.,ncl1);
      }
    }
  }

  Int_t nv0s = aod->GetNumberOfV0s();
  for (Int_t iV0 = 0; iV0 < nv0s; iV0++){
    AliAODv0 *v0 = aod->GetV0(iV0);
    if (!v0) continue;
    Bool_t onFlyStatus=v0->GetOnFlyStatus();
    if(onFlyStatus==kTRUE) continue;
    
    AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0); //0->Positive Daughter
    AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1); //1->Negative Daughter
    if (!pTrack || !nTrack) {
      Printf("ERROR: Could not retreive one of the daughter track");
      continue;
    }
    if(pTrack->GetID()<0 || nTrack->GetID()<0) continue;
    if (pTrack->Charge() == nTrack->Charge()) continue;

    Double_t invMassK0s = v0->MassK0Short();
    Double_t invMassLambda = v0->MassLambda();
    Double_t invMassAntiLambda = v0->MassAntiLambda();
    Double_t ptv0=v0->Pt();
    Double_t xv0=v0->Xv();
    Double_t yv0=v0->Yv();
    Double_t rv0=TMath::Sqrt(xv0*xv0+yv0*yv0);

    if(ConvertAndSelectAODTrack(pTrack,vESD,magField)==kFALSE) continue;
    if(ConvertAndSelectAODTrack(nTrack,vESD,magField)==kFALSE) continue;
    if(fRequireITSforV0dau & (1<<kBitRequireITSrefit)){
      if(!(pTrack->GetStatus() & AliESDtrack::kITSrefit)) continue;
      if(!(nTrack->GetStatus() & AliESDtrack::kITSrefit)) continue;
    }
    if(fRequireITSforV0dau & (1<<kBitRequireSPDany)){
      if(!pTrack->HasPointOnITSLayer(0) && !pTrack->HasPointOnITSLayer(1)) continue;
      if(!nTrack->HasPointOnITSLayer(0) && !nTrack->HasPointOnITSLayer(1)) continue;
    }

    Bool_t keepK0s=kTRUE;
    Bool_t keepLambda=kTRUE;
    Bool_t keepAntiLambda=kTRUE;
    if(!fReadMC){
      if(pidResp){
	Double_t nsigmap=-999.;
	if (pidResp->CheckPIDStatus(AliPIDResponse::kTPC,pTrack) == AliPIDResponse::kDetPidOk){
	  nsigmap=pidResp->NumberOfSigmasTPC(pTrack,AliPID::kProton);
	}
	Double_t nsigman=-999.;
	if (pidResp->CheckPIDStatus(AliPIDResponse::kTPC,nTrack) == AliPIDResponse::kDetPidOk){
	  nsigman=pidResp->NumberOfSigmasTPC(nTrack,AliPID::kProton);
	}
	if(TMath::Abs(nsigmap)>3) keepLambda=kFALSE;
	if(TMath::Abs(nsigman)>3) keepAntiLambda=kFALSE;
      }
    }else{
      keepK0s=kFALSE;
      keepLambda=kFALSE;
      keepAntiLambda=kFALSE;
      Int_t labelPos=TMath::Abs(pTrack->GetLabel());
      Int_t labelNeg =TMath::Abs(nTrack->GetLabel());
      AliAODMCParticle* partPos = dynamic_cast<AliAODMCParticle*>(arrayMC->At(TMath::Abs(labelPos)));
      AliAODMCParticle* partNeg = dynamic_cast<AliAODMCParticle*>(arrayMC->At(TMath::Abs(labelNeg)));
      if(partPos && partNeg){
        Int_t labelMotherPos=partPos->GetMother() ;
        Int_t labelMotherNeg=partNeg->GetMother();
	if(labelMotherPos==labelMotherNeg && labelMotherPos>-1){
	  AliAODMCParticle* partV0 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(TMath::Abs(labelMotherPos)));
 	  Int_t pdgV0=partV0->GetPdgCode();
	  if(TMath::Abs(pdgV0)==310) keepK0s=kTRUE;
	  if(pdgV0==3122) keepLambda=kTRUE;
	  if(pdgV0==-3122) keepAntiLambda=kTRUE;
	}
      }
    }
    if(TMath::Abs(invMassK0s-0.4976)<0.03){
      keepLambda=kFALSE;
      keepAntiLambda=kFALSE;
    }

    if(keepK0s) fHistInvMassK0s->Fill(invMassK0s,ptv0,rv0);
    if(keepLambda){
      fHistInvMassLambda->Fill(invMassLambda,ptv0,rv0);
    }
    if(keepAntiLambda){
      fHistInvMassAntiLambda->Fill(invMassAntiLambda,ptv0,rv0);
    }
  }
  PostData(1,fOutput);
  PostData(2,fTrackTree);
  
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskCheckAODTracks::ConvertAndSelectAODTrack(AliAODTrack* aTrack, const AliESDVertex vESD, Double_t magField)
{
  AliESDtrack esdTrack(aTrack);
  esdTrack.SetTPCClusterMap(aTrack->GetTPCClusterMap());
  esdTrack.SetTPCSharedMap(aTrack->GetTPCSharedMap());
  esdTrack.SetTPCPointsF(aTrack->GetTPCNclsF());
  esdTrack.SetTPCNcls(aTrack->GetTPCNcls());
  esdTrack.SetITSchi2(aTrack->GetITSchi2());
  Int_t nTPCclus=aTrack->GetNcls(1);
  Double_t chi2ndf=aTrack->Chi2perNDF();
  Double_t chi2tpc=999.;
  if(chi2ndf>0. && nTPCclus > 5){
    chi2tpc=Float_t(nTPCclus-5)*chi2ndf;
  }
  esdTrack.SetTPCchi2(chi2tpc);
  // needed to calculate the impact parameters
  Bool_t okDCA=esdTrack.RelateToVertex(&vESD,magField,99999.);
  if(!okDCA) return kFALSE;
  AliAODVertex* av=aTrack->GetProdVertex();
  if(av->GetType()==AliAODVertex::kKink) return kFALSE;
  if(aTrack->GetTPCsignalN()<fMinNumOfTPCPIDclu) return kFALSE;
  return fTrCutsTPC->AcceptTrack(&esdTrack);
}
//______________________________________________________________________________
void AliAnalysisTaskCheckAODTracks::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  printf("AliAnalysisTaskCheckAODTracks::Terminate --- Number of events: read = %.0f  analysed = %.0f\n",fHistNEvents->GetBinContent(1),fHistNEvents->GetBinContent(5));
  return;
}





