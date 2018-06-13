#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliPID.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include "AliAnalysisUtils.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TChain.h>
#include "AliPIDResponse.h"
#include "AliAnalysisTaskCheckESDTracks.h"


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
// Implementation of class AliAnalysisTaskCheckESDTracks
// AliAnalysisTaskSE to extract QA and performance histos for tracks
// 
//
// Authors: 
//          F. Prino, C. Zampolli
//          
//*************************************************************************

ClassImp(AliAnalysisTaskCheckESDTracks)
//______________________________________________________________________________
AliAnalysisTaskCheckESDTracks::AliAnalysisTaskCheckESDTracks() : 
  AliAnalysisTaskSE("QAofESDTracks"), 
  fOutput{nullptr},
  fHistNEvents{nullptr},
  fHistNTracks{nullptr},
  fHistNITSClu{nullptr},
  fHistCluInITSLay{nullptr},
  fHistNtracksTPCselVsV0befEvSel{nullptr},
  fHistNtracksSPDanyVsV0befEvSel{nullptr},
  fHistNtracksTPCselVsV0aftEvSel{nullptr},
  fHistNtracksSPDanyVsV0aftEvSel{nullptr},
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
  fHistEtaPhiPtInnerTPCsel{nullptr},
  fHistEtaPhiPtInnerTPCselITSref{nullptr},
  fHistEtaPhiPtInnerTPCselSPDany{nullptr},
  fHistEtaPhiPtInnerTPCselTOFbc{nullptr},
  fHistEtaPhiPtInnerTPCselITSrefTOFbc{nullptr},
  fHistEtaPhiPtInnerTPCselSPDanyTOFbc{nullptr},
  fHistNtrackeltsPtTPCsel{nullptr},
  fHistNtrackeltsPtTPCselITSref{nullptr},
  fHistNtrackeltsPtTPCselSPDany{nullptr},
  fHistNtrackeltsPtTPCselTOFbc{nullptr},
  fHistNtrackeltsPtTPCselITSrefTOFbc{nullptr},
  fHistNtrackeltsPtTPCselSPDanyTOFbc{nullptr},
  fHistTPCchi2PerClusPhiPtTPCsel{nullptr},
  fHistTPCchi2PerClusPhiPtTPCselITSref{nullptr},
  fHistTPCchi2PerClusPhiPtTPCselSPDany{nullptr},
  fHistSig1ptCovMatPhiPtTPCsel{nullptr},
  fHistSig1ptCovMatPhiPtTPCselITSref{nullptr},
  fHistSig1ptCovMatPhiPtTPCselSPDany{nullptr},
  fHistEtaPhiPtGoodHypTPCsel{nullptr},
  fHistEtaPhiPtGoodHypTPCselITSref{nullptr},
  fHistEtaPhiPtGoodHypTPCselSPDany{nullptr},
  fHistEtaPhiPtInnerGoodHypTPCsel{nullptr},
  fHistEtaPhiPtInnerGoodHypTPCselITSref{nullptr},
  fHistEtaPhiPtInnerGoodHypTPCselSPDany{nullptr},
  fHistTPCchi2PerClusPhiPtGoodHypTPCsel{nullptr},
  fHistTPCchi2PerClusPhiPtGoodHypTPCselITSref{nullptr},
  fHistTPCchi2PerClusPhiPtGoodHypTPCselSPDany{nullptr},
  fHistdEdxVsPGoodHyp{nullptr},
  fHistImpParXYPtMulGoodHypTPCselSPDany{nullptr},
  fHistEtaPhiPtBadHypTPCsel{nullptr},
  fHistEtaPhiPtBadHypTPCselITSref{nullptr},
  fHistEtaPhiPtBadHypTPCselSPDany{nullptr},
  fHistEtaPhiPtInnerBadHypTPCsel{nullptr},
  fHistEtaPhiPtInnerBadHypTPCselITSref{nullptr},
  fHistEtaPhiPtInnerBadHypTPCselSPDany{nullptr},
  fHistTPCchi2PerClusPhiPtBadHypTPCsel{nullptr},
  fHistTPCchi2PerClusPhiPtBadHypTPCselITSref{nullptr},
  fHistTPCchi2PerClusPhiPtBadHypTPCselSPDany{nullptr},
  fHistdEdxVsPBadHyp{nullptr},
  fHistImpParXYPtMulBadHypTPCselSPDany{nullptr},
  fHistPtResidVsPtTPCselAll{nullptr},
  fHistPtResidVsPtTPCselITSrefAll{nullptr},
  fHistPtResidVsPtTPCselGoodHyp{nullptr},
  fHistPtResidVsPtTPCselBadHyp{nullptr},
  fHistPtResidVsPtTPCselITSrefGoodHyp{nullptr},
  fHistPtResidVsPtTPCselITSrefBadHyp {nullptr},
  fHistOneOverPtResidVsPtTPCselAll{nullptr},
  fHistOneOverPtResidVsPtTPCselITSrefAll{nullptr},
  fHistOneOverPtResidVsPtTPCselGoodHyp{nullptr},
  fHistOneOverPtResidVsPtTPCselBadHyp{nullptr},
  fHistOneOverPtResidVsPtTPCselITSrefGoodHyp{nullptr},
  fHistOneOverPtResidVsPtTPCselITSrefBadHyp {nullptr},
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
  fHistInvMassLambdaGoodHyp{nullptr},
  fHistInvMassAntiLambdaGoodHyp{nullptr},
  fHistInvMassLambdaBadHyp{nullptr},
  fHistInvMassAntiLambdaBadHyp{nullptr},
  fFillTree(kFALSE),
  fTrackTree{nullptr},
  fTreeVarFloat{nullptr},
  fTreeVarInt{nullptr},
  fTrCutsTPC{nullptr},
  fMinNumOfTPCPIDclu(0),
  fUseTOFbcSelection(kTRUE),
  fUsePhysSel(kTRUE),
  fUsePileupCut(kTRUE),
  fTriggerMask(AliVEvent::kAnyINT),
  fNEtaBins(10),
  fNPhiBins(144),
  fNPtBins(100),
  fMinPt(0.),
  fMaxPt(25.),
  fReadMC(kFALSE),
  fUseMCId(kFALSE)
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
  for(Int_t jsp=0; jsp<9; jsp++){
    fHistdEdxVsP[jsp]=0x0;
    fHistdEdxVsPTPCsel[jsp]=0x0;
    fHistdEdxVsPTPCselITSref[jsp]=0x0;
    fHistnSigmaVsPdEdxTPCsel[jsp]=0x0;
  }
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}


//___________________________________________________________________________
AliAnalysisTaskCheckESDTracks::~AliAnalysisTaskCheckESDTracks(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;

  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistNTracks;
    delete fHistNITSClu;
    delete fHistCluInITSLay;
    delete fHistNtracksTPCselVsV0befEvSel;
    delete fHistNtracksSPDanyVsV0befEvSel;
    delete fHistNtracksTPCselVsV0aftEvSel;
    delete fHistNtracksSPDanyVsV0aftEvSel;
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
    delete fHistEtaPhiPtInnerTPCsel;
    delete fHistEtaPhiPtInnerTPCselITSref;
    delete fHistEtaPhiPtInnerTPCselSPDany;
    delete fHistEtaPhiPtInnerTPCselTOFbc;
    delete fHistEtaPhiPtInnerTPCselITSrefTOFbc;
    delete fHistEtaPhiPtInnerTPCselSPDanyTOFbc;
    delete fHistNtrackeltsPtTPCsel;
    delete fHistNtrackeltsPtTPCselITSref;
    delete fHistNtrackeltsPtTPCselSPDany;
    delete fHistNtrackeltsPtTPCselTOFbc;
    delete fHistNtrackeltsPtTPCselITSrefTOFbc;
    delete fHistNtrackeltsPtTPCselSPDanyTOFbc;
    delete fHistTPCchi2PerClusPhiPtTPCsel;
    delete fHistTPCchi2PerClusPhiPtTPCselITSref;
    delete fHistTPCchi2PerClusPhiPtTPCselSPDany;
    delete fHistSig1ptCovMatPhiPtTPCsel;
    delete fHistSig1ptCovMatPhiPtTPCselITSref;
    delete fHistSig1ptCovMatPhiPtTPCselSPDany;
    for(Int_t j=0; j<3; j++){
      delete fHistEtaPhiPtGoodHypTPCsel[j];
      delete fHistEtaPhiPtGoodHypTPCselITSref[j];
      delete fHistEtaPhiPtGoodHypTPCselSPDany[j];
      delete fHistEtaPhiPtInnerGoodHypTPCsel[j];
      delete fHistEtaPhiPtInnerGoodHypTPCselITSref[j];
      delete fHistEtaPhiPtInnerGoodHypTPCselSPDany[j];
      delete fHistTPCchi2PerClusPhiPtGoodHypTPCsel[j];
      delete fHistTPCchi2PerClusPhiPtGoodHypTPCselITSref[j];
      delete fHistTPCchi2PerClusPhiPtGoodHypTPCselSPDany[j];
      delete fHistdEdxVsPGoodHyp[j];
      delete fHistImpParXYPtMulGoodHypTPCselSPDany[j];
      delete fHistEtaPhiPtBadHypTPCsel[j];
      delete fHistEtaPhiPtBadHypTPCselITSref[j];
      delete fHistEtaPhiPtBadHypTPCselSPDany[j];
      delete fHistEtaPhiPtInnerBadHypTPCsel[j];
      delete fHistEtaPhiPtInnerBadHypTPCselITSref[j];
      delete fHistEtaPhiPtInnerBadHypTPCselSPDany[j];
      delete fHistTPCchi2PerClusPhiPtBadHypTPCsel[j];
      delete fHistTPCchi2PerClusPhiPtBadHypTPCselITSref[j];
      delete fHistTPCchi2PerClusPhiPtBadHypTPCselSPDany[j];
      delete fHistdEdxVsPBadHyp[j];
      delete fHistImpParXYPtMulBadHypTPCselSPDany[j];
    }
    delete fHistPtResidVsPtTPCselAll;
    delete fHistPtResidVsPtTPCselITSrefAll;
    delete fHistOneOverPtResidVsPtTPCselAll;
    delete fHistOneOverPtResidVsPtTPCselITSrefAll;
    for (int iS = 0; iS < AliPID::kSPECIESC;++iS) {
      delete fHistPtResidVsPtTPCselGoodHyp[iS];
      delete fHistPtResidVsPtTPCselBadHyp[iS];
      delete fHistPtResidVsPtTPCselITSrefGoodHyp[iS];
      delete fHistPtResidVsPtTPCselITSrefBadHyp[iS];
      delete fHistOneOverPtResidVsPtTPCselGoodHyp[iS];
      delete fHistOneOverPtResidVsPtTPCselBadHyp[iS];
      delete fHistOneOverPtResidVsPtTPCselITSrefGoodHyp[iS];
      delete fHistOneOverPtResidVsPtTPCselITSrefBadHyp[iS];
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
    delete fHistInvMassLambdaGoodHyp;
    delete fHistInvMassAntiLambdaGoodHyp;
    delete fHistInvMassLambdaBadHyp;
    delete fHistInvMassAntiLambdaBadHyp;

    for(Int_t jsp=0; jsp<9; jsp++){
      delete fHistdEdxVsP[jsp];
      delete fHistdEdxVsPTPCsel[jsp];
      delete fHistdEdxVsPTPCselITSref[jsp];
      delete fHistnSigmaVsPdEdxTPCsel[jsp];
    }
    delete fTrackTree;
  }
  delete fOutput;
  delete fTrCutsTPC;
  delete [] fTreeVarFloat;
  delete [] fTreeVarInt;
}
 
//___________________________________________________________________________
void AliAnalysisTaskCheckESDTracks::UserCreateOutputObjects() {
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
  floatVarName[10]="pxtpc";
  floatVarName[11]="pytpc";
  floatVarName[12]="pztpc";
  floatVarName[13]="ptttpc";
  floatVarName[14]="ptpc";
  floatVarName[15]="etatpc";
  floatVarName[16]="phitpc";
  floatVarName[17]="phipostpc";
  floatVarName[18]="d0xy";
  floatVarName[19]="d0z";
  floatVarName[20]="chi2clustpc";
  floatVarName[21]="croverfind";
  floatVarName[22]="dedxTPC";
  floatVarName[23]="nsigel";
  floatVarName[24]="nsigpi";
  floatVarName[25]="nsigk";
  floatVarName[26]="nsigp";
  floatVarName[27]="pxgen";
  floatVarName[28]="pygen";
  floatVarName[29]="pzgen";
  floatVarName[30]="ptgen";
  floatVarName[31]="pgen";
  floatVarName[32]="etagen";
  floatVarName[33]="phigen";
  Int_t usedVar=kNumOfFloatVar-7;
  if(fReadMC) usedVar=kNumOfFloatVar;
  for(Int_t ivar=0; ivar<usedVar; ivar++){
    fTrackTree->Branch(floatVarName[ivar].Data(),&fTreeVarFloat[ivar],Form("%s/F",floatVarName[ivar].Data()));
  }
  

  TString intVarName[kNumOfIntVar];
  fTreeVarInt = new Int_t[kNumOfIntVar];
  intVarName[0]="ntrack";
  intVarName[1]="ntracklets";
  intVarName[2]="charge";
  intVarName[3]="ITSrefit";
  intVarName[4]="ITSclumap";
  intVarName[5]="nTPCclu";
  intVarName[6]="TOFbc";
  intVarName[7]="trackPIDhip";
  intVarName[8]="label";
  intVarName[9]="truePID";
  usedVar=kNumOfIntVar-2;
  if(fReadMC) usedVar=kNumOfIntVar;
  for(Int_t ivar=0; ivar<usedVar; ivar++){
    fTrackTree->Branch(intVarName[ivar].Data(),&fTreeVarInt[ivar],Form("%s/I",intVarName[ivar].Data()));
  }

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",15,-0.5,14.5);
  //fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"PhysSel"); 
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Good vertex"); 
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Pass zSPD-zTrk vert sel"); 
  fHistNEvents->GetXaxis()->SetBinLabel(5,"|zvert|<10"); 
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Pileup cut"); 
  fOutput->Add(fHistNEvents);

  fHistNTracks = new TH1F("hNTracks", "Number of tracks in ESD events ; N_{tracks}",5001,-0.5,5000.5);
  fOutput->Add(fHistNTracks);
  fHistNITSClu = new TH1F("hNITSClu"," ; N_{ITS clusters}",7,-0.5,6.5);
  fOutput->Add(fHistNITSClu);
  fHistCluInITSLay = new TH1F("hCluInITSLay"," ; Layer",7,-1.5,5.5);
  fOutput->Add(fHistCluInITSLay);

  fHistNtracksTPCselVsV0befEvSel=new TH2F("hNtracksTPCselVsV0befEvSel"," ; V0 signal ; N_{tracks}",250,0.,15000.,250,0.,5000.);
  fHistNtracksSPDanyVsV0befEvSel=new TH2F("hNtracksSPDanyVsV0befEvSel"," ; V0 signal ; N_{tracks}",250,0.,15000.,250,0.,5000.);
  fHistNtracksTPCselVsV0aftEvSel=new TH2F("hNtracksTPCselVsV0aftEvSel"," ; V0 signal ; N_{tracks}",250,0.,15000.,250,0.,5000.);
  fHistNtracksSPDanyVsV0aftEvSel=new TH2F("hNtracksSPDanyVsV0aftEvSel"," ; V0 signal ; N_{tracks}",250,0.,15000.,250,0.,5000.);
  fOutput->Add(fHistNtracksTPCselVsV0befEvSel);
  fOutput->Add(fHistNtracksSPDanyVsV0befEvSel);
  fOutput->Add(fHistNtracksTPCselVsV0aftEvSel);
  fOutput->Add(fHistNtracksSPDanyVsV0aftEvSel);

  TString pNames[9]={"Elec","Muon","Pion","Kaon","Proton","Deuteron","Triton","He3","Alpha"};
  for(Int_t jsp=0; jsp<9; jsp++){ 
    fHistdEdxVsP[jsp] = new TH2F(Form("hdEdxVsP%s",pNames[jsp].Data()),"  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
    fHistdEdxVsPTPCsel[jsp] = new TH2F(Form("hdEdxVsPTPCsel%s",pNames[jsp].Data()),"  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
    fHistdEdxVsPTPCselITSref[jsp] = new TH2F(Form("hdEdxVsPTPCselITSref%s",pNames[jsp].Data()),"  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
    fHistnSigmaVsPdEdxTPCsel[jsp] = new TH2F(Form("hnSigmaVsPdEdxTPCsel%s",pNames[jsp].Data()),Form("  ; p_{TPC} (GeV/c) ; n#sigma(%s)",pNames[jsp].Data()),100,0.,5.,200,-10.,10.);
    fOutput->Add(fHistdEdxVsP[jsp]);
    fOutput->Add(fHistdEdxVsPTPCsel[jsp]);
    fOutput->Add(fHistdEdxVsPTPCselITSref[jsp]);
    fOutput->Add(fHistnSigmaVsPdEdxTPCsel[jsp]);
  }

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
  fOutput->Add(fHistEtaPhiPtTPCselTOFbc);
  fOutput->Add(fHistEtaPhiPtTPCselITSrefTOFbc);
  fOutput->Add(fHistEtaPhiPtTPCselSPDanyTOFbc);

  fHistEtaPhiPtInnerTPCsel = new TH3F("hEtaPhiPtInnerTPCsel"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtInnerTPCselITSref = new TH3F("hEtaPhiPtInnerTPCselITSref"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtInnerTPCselSPDany = new TH3F("hEtaPhiPtInnerTPCselSPDany"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fOutput->Add(fHistEtaPhiPtInnerTPCsel);
  fOutput->Add(fHistEtaPhiPtInnerTPCselITSref);
  fOutput->Add(fHistEtaPhiPtInnerTPCselSPDany);

  fHistEtaPhiPtInnerTPCselTOFbc = new TH3F("hEtaPhiPtInnerTPCselTOFbc"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtInnerTPCselITSrefTOFbc = new TH3F("hEtaPhiPtInnerTPCselITSrefTOFbc"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fHistEtaPhiPtInnerTPCselSPDanyTOFbc = new TH3F("hEtaPhiPtInnerTPCselSPDanyTOFbc"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
  fOutput->Add(fHistEtaPhiPtInnerTPCselTOFbc);
  fOutput->Add(fHistEtaPhiPtInnerTPCselITSrefTOFbc);
  fOutput->Add(fHistEtaPhiPtInnerTPCselSPDanyTOFbc);

  fHistNtrackeltsPtTPCsel = new TH2F("hNtrackeltsPtTPCsel"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,fNPtBins,fMinPt,fMaxPt);
  fHistNtrackeltsPtTPCselITSref = new TH2F("hNtrackeltsPtTPCselITSref"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,fNPtBins,fMinPt,fMaxPt);
  fHistNtrackeltsPtTPCselSPDany = new TH2F("hNtrackeltsPtTPCselSPDany"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,fNPtBins,fMinPt,fMaxPt);
  fHistNtrackeltsPtTPCselTOFbc = new TH2F("hNtrackeltsPtTPCselTOFbc"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,fNPtBins,fMinPt,fMaxPt);
  fHistNtrackeltsPtTPCselITSrefTOFbc = new TH2F("hNtrackeltsPtTPCselITSrefTOFbc"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,fNPtBins,fMinPt,fMaxPt);
  fHistNtrackeltsPtTPCselSPDanyTOFbc = new TH2F("hNtrackeltsPtTPCselSPDanyTOFbc"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,fNPtBins,fMinPt,fMaxPt);
  fOutput->Add(fHistNtrackeltsPtTPCsel);
  fOutput->Add(fHistNtrackeltsPtTPCselITSref);
  fOutput->Add(fHistNtrackeltsPtTPCselSPDany);
  fOutput->Add(fHistNtrackeltsPtTPCselTOFbc);
  fOutput->Add(fHistNtrackeltsPtTPCselITSrefTOFbc);
  fOutput->Add(fHistNtrackeltsPtTPCselSPDanyTOFbc);
 

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
 
  // Impact parameter binning
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

  for(Int_t j=0; j<3; j++){
    
    fHistEtaPhiPtGoodHypTPCsel[j] = new TH3F(Form("hEtaPhiPtGoodHyp%sTPCsel",pNames[j+2].Data())," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtGoodHypTPCselITSref[j] = new TH3F(Form("hEtaPhiPtGoodHyp%sTPCselITSref",pNames[j+2].Data())," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtGoodHypTPCselSPDany[j] = new TH3F(Form("hEtaPhiPtGoodHyp%sTPCselSPDany",pNames[j+2].Data())," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtInnerGoodHypTPCsel[j] = new TH3F(Form("hEtaPhiPtInnerGoodHyp%sTPCsel",pNames[j+2].Data())," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtInnerGoodHypTPCselITSref[j] = new TH3F(Form("hEtaPhiPtInnerGoodHyp%sTPCselITSref",pNames[j+2].Data())," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtInnerGoodHypTPCselSPDany[j] = new TH3F(Form("hEtaPhiPtInnerGoodHyp%sTPCselSPDany",pNames[j+2].Data())," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistTPCchi2PerClusPhiPtGoodHypTPCsel[j] = new TH3F(Form("hTPCchi2PerClusPhiPtGoodHyp%sTPCsel",pNames[j+2].Data())," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
    fHistTPCchi2PerClusPhiPtGoodHypTPCselITSref[j] = new TH3F(Form("hTPCchi2PerClusPhiPtGoodHyp%sTPCselITSref",pNames[j+2].Data())," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
    fHistTPCchi2PerClusPhiPtGoodHypTPCselSPDany[j] = new TH3F(Form("hTPCchi2PerClusPhiPtGoodHyp%sTPCselSPDany",pNames[j+2].Data())," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
    fHistdEdxVsPGoodHyp[j] = new TH2F(Form("hdEdxVsPGoodHyp%s",pNames[j+2].Data()),"  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
    fHistImpParXYPtMulGoodHypTPCselSPDany[j] = new TH3F(Form("hImpParXYPtMulGoodHyp%sTPCselSPDany",pNames[j+2].Data())," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip); 
    fOutput->Add(fHistEtaPhiPtGoodHypTPCsel[j]);
    fOutput->Add(fHistEtaPhiPtGoodHypTPCselITSref[j]);
    fOutput->Add(fHistEtaPhiPtGoodHypTPCselSPDany[j]);
    fOutput->Add(fHistEtaPhiPtInnerGoodHypTPCsel[j]);
    fOutput->Add(fHistEtaPhiPtInnerGoodHypTPCselITSref[j]);
    fOutput->Add(fHistEtaPhiPtInnerGoodHypTPCselSPDany[j]);
    fOutput->Add(fHistTPCchi2PerClusPhiPtGoodHypTPCsel[j]);
    fOutput->Add(fHistTPCchi2PerClusPhiPtGoodHypTPCselITSref[j]);
    fOutput->Add(fHistTPCchi2PerClusPhiPtGoodHypTPCselSPDany[j]);
    fOutput->Add(fHistdEdxVsPGoodHyp[j]);
    fOutput->Add(fHistImpParXYPtMulGoodHypTPCselSPDany[j]);

    fHistEtaPhiPtBadHypTPCsel[j] = new TH3F(Form("hEtaPhiPtBadHyp%sTPCsel",pNames[j+2].Data())," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtBadHypTPCselITSref[j] = new TH3F(Form("hEtaPhiPtBadHyp%sTPCselITSref",pNames[j+2].Data())," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtBadHypTPCselSPDany[j] = new TH3F(Form("hEtaPhiPtBadHyp%sTPCselSPDany",pNames[j+2].Data())," ; #eta ; #varphi ; p_{T} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtInnerBadHypTPCsel[j] = new TH3F(Form("hEtaPhiPtInnerBadHyp%sTPCsel",pNames[j+2].Data())," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtInnerBadHypTPCselITSref[j] = new TH3F(Form("hEtaPhiPtInnerBadHyp%sTPCselITSref",pNames[j+2].Data())," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistEtaPhiPtInnerBadHypTPCselSPDany[j] = new TH3F(Form("hEtaPhiPtInnerBadHyp%sTPCselSPDany",pNames[j+2].Data())," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",fNEtaBins,-1.,1.,fNPhiBins,0.,2*TMath::Pi(),fNPtBins,fMinPt,fMaxPt);
    fHistTPCchi2PerClusPhiPtBadHypTPCsel[j] = new TH3F(Form("hTPCchi2PerClusPhiPtBadHyp%sTPCsel",pNames[j+2].Data())," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
    fHistTPCchi2PerClusPhiPtBadHypTPCselITSref[j] = new TH3F(Form("hTPCchi2PerClusPhiPtBadHyp%sTPCselITSref",pNames[j+2].Data())," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
    fHistTPCchi2PerClusPhiPtBadHypTPCselSPDany[j] = new TH3F(Form("hTPCchi2PerClusPhiPtBadHyp%sTPCselSPDany",pNames[j+2].Data())," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
    fHistdEdxVsPBadHyp[j] = new TH2F(Form("hdEdxVsPBadHyp%s",pNames[j+2].Data()),"  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
    fHistImpParXYPtMulBadHypTPCselSPDany[j] = new TH3F(Form("hImpParXYPtMulBadHyp%sTPCselSPDany",pNames[j+2].Data())," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
    fOutput->Add(fHistEtaPhiPtBadHypTPCsel[j]);
    fOutput->Add(fHistEtaPhiPtBadHypTPCselITSref[j]);
    fOutput->Add(fHistEtaPhiPtBadHypTPCselSPDany[j]);
    fOutput->Add(fHistEtaPhiPtInnerBadHypTPCsel[j]);
    fOutput->Add(fHistEtaPhiPtInnerBadHypTPCselITSref[j]);
    fOutput->Add(fHistEtaPhiPtInnerBadHypTPCselSPDany[j]);
    fOutput->Add(fHistTPCchi2PerClusPhiPtBadHypTPCsel[j]);
    fOutput->Add(fHistTPCchi2PerClusPhiPtBadHypTPCselITSref[j]);
    fOutput->Add(fHistTPCchi2PerClusPhiPtBadHypTPCselSPDany[j]);
    fOutput->Add(fHistdEdxVsPBadHyp[j]);
    fOutput->Add(fHistImpParXYPtMulBadHypTPCselSPDany[j]);
  }

  fHistPtResidVsPtTPCselAll = new TH2F("hPtResidVsPtTPCselAll","; p_{T,rec} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefAll = new TH2F("hPtResidVsPtTPCselITSrefAll","; p_{T,rec} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselAll = new TH2F("hOneOverPtResidVsPtTPCselAll","; p_{T,rec} (GeV/c) ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
  fHistOneOverPtResidVsPtTPCselITSrefAll = new TH2F("hOneOverPtResidVsPtTPCselITSrefAll","; p_{T,rec} (GeV/c) ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
  fOutput->Add(fHistPtResidVsPtTPCselAll);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefAll);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselAll);
  fOutput->Add(fHistOneOverPtResidVsPtTPCselITSrefAll);
  for (int iS = 0; iS < AliPID::kSPECIESC; ++iS) {
    fHistPtResidVsPtTPCselGoodHyp[iS] = new TH2F(Form("hPtResidVsPtTPCselGoodHyp%s",AliPID::ParticleShortName(iS)), Form("%s ; p_{T,rec} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",AliPID::ParticleLatexName(iS)),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistPtResidVsPtTPCselBadHyp[iS] = new TH2F(Form("hPtResidVsPtTPCselBadHyp%s",AliPID::ParticleShortName(iS)),Form("%s ; p_{T,rec} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",AliPID::ParticleLatexName(iS)),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistPtResidVsPtTPCselITSrefGoodHyp[iS] = new TH2F(Form("hPtResidVsPtTPCselITSrefGoodHyp%s",AliPID::ParticleShortName(iS)),Form("%s ; p_{T,rec} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",AliPID::ParticleLatexName(iS)),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistPtResidVsPtTPCselITSrefBadHyp[iS] = new TH2F(Form("hPtResidVsPtTPCselITSrefBadHyp%s",AliPID::ParticleShortName(iS)),Form("%s ; p_{T,rec} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",AliPID::ParticleLatexName(iS)),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistOneOverPtResidVsPtTPCselGoodHyp[iS] = new TH2F(Form("hOneOverPtResidVsPtTPCselGoodHyp%s",AliPID::ParticleShortName(iS)), Form("%s ; p_{T,rec} (GeV/c) ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",AliPID::ParticleLatexName(iS)),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistOneOverPtResidVsPtTPCselBadHyp[iS] = new TH2F(Form("hOneOverPtResidVsPtTPCselBadHyp%s",AliPID::ParticleShortName(iS)),Form("%s ; p_{T,rec} (GeV/c) ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",AliPID::ParticleLatexName(iS)),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistOneOverPtResidVsPtTPCselITSrefGoodHyp[iS] = new TH2F(Form("hOneOverPtResidVsPtTPCselITSrefGoodHyp%s",AliPID::ParticleShortName(iS)),Form("%s ; p_{T,rec} (GeV/c) ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",AliPID::ParticleLatexName(iS)),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fHistOneOverPtResidVsPtTPCselITSrefBadHyp[iS] = new TH2F(Form("hOneOverPtResidVsPtTPCselITSrefBadHyp%s",AliPID::ParticleShortName(iS)),Form("%s ; p_{T,rec} (GeV/c) ; p_{T,rec}*(1/p_{T,rec}-1/p_{T,gen})",AliPID::ParticleLatexName(iS)),fNPtBins,fMinPt,fMaxPt,100,-0.5,0.5);
    fOutput->Add(fHistPtResidVsPtTPCselGoodHyp[iS]);
    fOutput->Add(fHistPtResidVsPtTPCselBadHyp[iS]);
    fOutput->Add(fHistPtResidVsPtTPCselITSrefGoodHyp[iS]);
    fOutput->Add(fHistPtResidVsPtTPCselITSrefBadHyp[iS]);
    fOutput->Add(fHistOneOverPtResidVsPtTPCselGoodHyp[iS]);
    fOutput->Add(fHistOneOverPtResidVsPtTPCselBadHyp[iS]);
    fOutput->Add(fHistOneOverPtResidVsPtTPCselITSrefGoodHyp[iS]);
    fOutput->Add(fHistOneOverPtResidVsPtTPCselITSrefBadHyp[iS]);
  }
  fHistPzResidVsPtTPCselAll = new TH2F("hPzResidVsPtTPCselAll","; p_{T,rec} (GeV/c) ; p_{z,rec}-p_{z,gen} (GeV/c)",fNPtBins,fMinPt,fMaxPt,150,-0.3,0.3);
  fHistPzResidVsPtTPCselITSrefAll = new TH2F("hPzResidVsPtTPCselITSrefAll","; p_{T,rec} (GeV/c) ; p_{z,rec}-p_{z,gen} (GeV/c)",fNPtBins,fMinPt,fMaxPt,150,-0.3,0.3);
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


  fHistImpParXYPtMulTPCselSPDanyGood = new TH3F("hImpParXYPtMulTPCselSPDanyGood"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fHistImpParXYPtMulTPCselSPDanyFake = new TH3F("hImpParXYPtMulTPCselSPDanyFake"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fHistImpParXYPtMulTPCselSPDanyPrim = new TH3F("hImpParXYPtMulTPCselSPDanyPrim"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fHistImpParXYPtMulTPCselSPDanySecDec = new TH3F("hImpParXYPtMulTPCselSPDanySecDec"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fHistImpParXYPtMulTPCselSPDanySecMat = new TH3F("hImpParXYPtMulTPCselSPDanySecMat"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",nPtBins4ip,ptBins4ip,nIPBins4ip,ipBins4ip,nMultBins4ip,multBins4ip);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyGood);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyFake);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyPrim);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanySecDec);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanySecMat);

  fHistInvMassK0s = new TH3F("hInvMassK0s"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(K0s) ; R (cm)",200,0.4,0.6,25,0.,5.,50,0.,50.);
  fHistInvMassLambda = new TH3F("hInvMassLambda"," ;Inv.Mass (GeV/c^{2}) ; p_{T}(#Lambda) ; R (cm)",200,1.0,1.2,25,0.,5.,50,0.,50.);
  fHistInvMassAntiLambda = new TH3F("hInvMassAntiLambda"," ;Inv.Mass (GeV/c^{2}) ; p_{T}(#bar{#Lambda}) ; R (cm)",200,1.0,1.2,25,0.,5.,50,0.,50.);
  fHistInvMassLambdaGoodHyp = new TH3F("hInvMassLambdaGoodHyp"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#Lambda) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fHistInvMassAntiLambdaGoodHyp = new TH3F("hInvMassAntiLambdaGoodHyp"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#bar{#Lambda}) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fHistInvMassLambdaBadHyp = new TH3F("hInvMassLambdaBadHyp"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#Lambda) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fHistInvMassAntiLambdaBadHyp = new TH3F("hInvMassAntiLambdaBadHyp"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#bar{#Lambda}) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fOutput->Add(fHistInvMassK0s);
  fOutput->Add(fHistInvMassLambda);
  fOutput->Add(fHistInvMassAntiLambda);
  fOutput->Add(fHistInvMassLambdaGoodHyp);
  fOutput->Add(fHistInvMassAntiLambdaGoodHyp);
  fOutput->Add(fHistInvMassLambdaBadHyp);
  fOutput->Add(fHistInvMassAntiLambdaBadHyp);

  PostData(1,fOutput);
  PostData(2,fTrackTree);

}
//______________________________________________________________________________
void AliAnalysisTaskCheckESDTracks::UserExec(Option_t *)
{
  //

  AliESDEvent *esd = (AliESDEvent*) (InputEvent());
  if(!esd) {
    printf("AliAnalysisTaskCheckESDTracks::UserExec(): bad ESD\n");
    return;
  } 

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
  
  AliMCEvent* mcEvent = nullptr;

  if(fReadMC){
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
  }


  fHistNEvents->Fill(0);
  if(fUsePhysSel){
    Bool_t isPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
    if(!isPhysSel) return;
  }
  fHistNEvents->Fill(1);



  Int_t ntracks = esd->GetNumberOfTracks();
  Int_t ntracklets = 0;
  Int_t ncl1 = 0;
  const AliMultiplicity* mult=esd->GetMultiplicity();
  if(mult){
    ntracklets = mult->GetNumberOfTracklets();
    ncl1 = mult->GetNumberOfITSClusters(1);
  }
  Int_t ntracksTPCsel=0;
  Int_t ntracksSPDany=0;
  for (Int_t iTrack=0; iTrack < ntracks; iTrack++) {
    AliESDtrack * track = (AliESDtrack*)esd->GetTrack(iTrack);
    if (!track) continue;
    track->SetESDEvent(esd);
    if(fTrCutsTPC->AcceptTrack(track)){
      ntracksTPCsel++;
      Int_t statusTrack=track->GetStatus();
      if(statusTrack&AliESDtrack::kITSrefit){
	if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)){
	  ntracksSPDany++;
	}
      }
    }
  }
  Double_t vZEROampl=0;
  for(Int_t i=0;i<64;i++) vZEROampl+=esd->GetVZEROData()->GetMultiplicity(i);
  fHistNtracksTPCselVsV0befEvSel->Fill(vZEROampl,ntracksTPCsel);
  fHistNtracksSPDanyVsV0befEvSel->Fill(vZEROampl,ntracksSPDany);

  const AliVVertex* vtTrc = esd->GetPrimaryVertex();
  const AliVVertex* vtSPD = esd->GetPrimaryVertexSPD();
  TString titTrc=vtTrc->GetTitle();
  if(titTrc.IsNull() || titTrc=="vertexer: 3D" || titTrc=="vertexer: Z") return;
  if (vtSPD->GetNContributors()<1) return;
  fHistNEvents->Fill(2);

  double covTrc[6],covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = vtTrc->GetZ()-vtSPD->GetZ();
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
  if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return; // bad vertexing
  fHistNEvents->Fill(3);

  Float_t xvert=vtTrc->GetX();
  Float_t yvert=vtTrc->GetY();
  Float_t zvert=vtTrc->GetZ();

  if(TMath::Abs(zvert)>10) return;
  fHistNEvents->Fill(4);

  if(fUsePileupCut){
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(5);
    utils.SetMaxPlpChi2MV(5.);
    utils.SetMinWDistMV(15.);
    utils.SetCheckPlpFromDifferentBCMV(kTRUE);
    Bool_t isPUMV = utils.IsPileUpMV(esd);
    if(isPUMV) return;
    fHistNEvents->Fill(5);
  }

  fHistNtracksTPCselVsV0aftEvSel->Fill(vZEROampl,ntracksTPCsel);
  fHistNtracksSPDanyVsV0aftEvSel->Fill(vZEROampl,ntracksSPDany);


  fHistNTracks->Fill(ntracks);


  for (Int_t iTrack=0; iTrack < ntracks; iTrack++) {
    AliESDtrack * track = esd->GetTrack(iTrack);
    if (!track) continue;
    track->SetESDEvent(esd);
    for(Int_t jvar=0; jvar<kNumOfFloatVar; jvar++) fTreeVarFloat[jvar]=-999.;
    for(Int_t jvar=0; jvar<kNumOfIntVar; jvar++) fTreeVarInt[jvar]=-999;

    fTreeVarFloat[0]=xvert;
    fTreeVarFloat[1]=yvert;
    fTreeVarFloat[2]=zvert;
    fTreeVarInt[0]=ntracks;
    fTreeVarInt[1]=ntracklets;

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
    const AliExternalTrackParam* ippar=track->GetInnerParam();
    Double_t ptrackTPC=-999.;
    Double_t pttrackTPC=-999.;
    Double_t phitrackTPC=-999.;
    Double_t etatrackTPC=-999.;
    if(ippar){
      ptrackTPC=ippar->P();
      pttrackTPC=ippar->Pt();
      phitrackTPC=ippar->Phi();
      etatrackTPC=ippar->Eta();
      fTreeVarFloat[10]=ippar->Px();
      fTreeVarFloat[11]=ippar->Py();
      fTreeVarFloat[12]=ippar->Pz();
      fTreeVarFloat[13]=ippar->Pt();
      fTreeVarFloat[14]=ippar->P();
      fTreeVarFloat[15]=ippar->Eta();
      fTreeVarFloat[16]=ippar->Phi();
      fTreeVarFloat[17]=ippar->PhiPos();
    }
    Float_t impactXY=-999, impactZ=-999;
    track->GetImpactParameters(impactXY, impactZ);
    Bool_t itsRefit=kFALSE;
    Int_t statusTrack=track->GetStatus();
    if(statusTrack&AliESDtrack::kITSrefit) itsRefit=kTRUE; 
    UChar_t clumap=track->GetITSClusterMap();
    Double_t chi2clus = track->GetTPCNcls() ? track->GetTPCchi2()/track->GetTPCNcls() : -999.;
    Double_t curvrelerr = TMath::Sqrt(track->GetSigma1Pt2())/track->OneOverPt();
    Bool_t spdAny=kFALSE;
    if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) spdAny=kTRUE;
    Int_t nTPCclus=track->GetNcls(1);
    Float_t nCrossedRowsTPC = track->GetTPCCrossedRows();
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (track->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
    }
    Bool_t tofOK=kTRUE;
    Int_t nrun=esd->GetRunNumber();
    if(nrun>=244913 && nrun<=246994){
      // more relaxed cut in PbPb with 100ns spacing 
      if(TMath::Abs(track->GetTOFExpTDiff())>30) tofOK=kFALSE;
    }else{
      if(TMath::Abs(track->GetTOFExpTDiff())>12.5) tofOK=kFALSE;
      // NOTE: track->GetTOFBunchCrossing()==0 (used in ITSqa task)
      // is equivalent to TMath::Abs(trc->GetTOFExpTDiff())<12.5 
      // since aliased to TMath::Nint(  trc->GetTOFExpTDiff()/25 )
    }
    fTreeVarFloat[18]=impactXY;
    fTreeVarFloat[19]=impactZ;
    fTreeVarFloat[20]=chi2clus;
    fTreeVarFloat[21]=ratioCrossedRowsOverFindableClustersTPC;

    fTreeVarInt[2]=chtrack;
    fTreeVarInt[3]=itsRefit;
    fTreeVarInt[4]=clumap;
    fTreeVarInt[5]=nTPCclus;
    fTreeVarInt[6]=tofOK;

    Int_t trlabel=track->GetLabel();
    Float_t dedx=track->GetTPCsignal();
    Int_t  pidtr=track->GetPIDForTracking();
    
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
    fTreeVarFloat[22]=dedx;
    fTreeVarFloat[23]=nSigmaTPC[0];
    fTreeVarFloat[24]=nSigmaTPC[2];
    fTreeVarFloat[25]=nSigmaTPC[3];
    fTreeVarFloat[26]=nSigmaTPC[4];
    fTreeVarInt[7]=pidtr;
    
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
      TParticle* part = mcEvent->Particle(TMath::Abs(trlabel));
      if (!part) continue;
      ptgen=part->Pt();
      pgen=part->P();
      pxgen=part->Px();
      pygen=part->Py();
      pzgen=part->Pz();
      if(ptgen>0.) invptgen=1./ptgen;
      etagen=part->Eta();
      phigen=part->Phi();
      if(mcEvent->IsPhysicalPrimary(TMath::Abs(trlabel))) isPhysPrim=1;
      else if(mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(trlabel))) isPhysPrim=0;
      else if(mcEvent->IsSecondaryFromMaterial(TMath::Abs(trlabel))) isPhysPrim=-1;
      fTreeVarFloat[27]=pxgen;
      fTreeVarFloat[28]=pygen;
      fTreeVarFloat[29]=pzgen;
      fTreeVarFloat[30]=ptgen;
      fTreeVarFloat[31]=pgen;
      fTreeVarFloat[32]=etagen;
      fTreeVarFloat[33]=phigen;
      if (fUseMCId) {
        int pdg = TMath::Abs(part->GetPdgCode());
        for (int iS = 0; iS < AliPID::kSPECIESC; ++iS) {
          if (pdg == AliPID::ParticleCode(iS)) hadronSpecies=iS;
        }
      }
      fTreeVarInt[8]=trlabel;
      fTreeVarInt[9]=part->GetPdgCode();
    }

    if(pidtr>=0 && pidtr<9) fHistdEdxVsP[pidtr]->Fill(ptrackTPC,dedx);
    
    if(!fTrCutsTPC->AcceptTrack(track)) continue;
    if(track->GetTPCsignalN()<fMinNumOfTPCPIDclu) continue;
    fHistNITSClu->Fill(track->GetNcls(0));
    fHistCluInITSLay->Fill(-1);
    for(Int_t iBit=0; iBit<6; iBit++){
      if(clumap&(1<<iBit)) fHistCluInITSLay->Fill(iBit);
    }

    if (fFillTree) fTrackTree->Fill();    

    fHistEtaPhiPtTPCsel->Fill(etatrack,phitrack,pttrack);
    if(chtrack>0) fHistEtaPhiPtPosChargeTPCsel->Fill(etatrack,phitrack,pttrack);
    else if(chtrack<0) fHistEtaPhiPtNegChargeTPCsel->Fill(etatrack,phitrack,pttrack);
    fHistEtaPhiPtInnerTPCsel->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
    fHistNtrackeltsPtTPCsel->Fill(ntracklets,pttrack);
    if(tofOK){
      fHistEtaPhiPtTPCselTOFbc->Fill(etatrack,phitrack,pttrack);
      fHistEtaPhiPtInnerTPCselTOFbc->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
      fHistNtrackeltsPtTPCselTOFbc->Fill(ntracklets,pttrack);
    }
    if(itsRefit){
      fHistEtaPhiPtTPCselITSref->Fill(etatrack,phitrack,pttrack);
      if(chtrack>0) fHistEtaPhiPtPosChargeTPCselITSref->Fill(etatrack,phitrack,pttrack);
      else if(chtrack<0) fHistEtaPhiPtNegChargeTPCselITSref->Fill(etatrack,phitrack,pttrack);
      fHistEtaPhiPtInnerTPCselITSref->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
      fHistNtrackeltsPtTPCselITSref->Fill(ntracklets,pttrack);
      if(tofOK){
	fHistEtaPhiPtTPCselITSrefTOFbc->Fill(etatrack,phitrack,pttrack);
	fHistEtaPhiPtInnerTPCselITSrefTOFbc->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	fHistNtrackeltsPtTPCselITSrefTOFbc->Fill(ntracklets,pttrack);
      }
      if(spdAny){ 
	fHistEtaPhiPtTPCselSPDany->Fill(etatrack,phitrack,pttrack);
	if(chtrack>0) fHistEtaPhiPtPosChargeTPCselSPDany->Fill(etatrack,phitrack,pttrack);
	else if(chtrack<0) fHistEtaPhiPtNegChargeTPCselSPDany->Fill(etatrack,phitrack,pttrack);
	fHistEtaPhiPtInnerTPCselSPDany->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	fHistNtrackeltsPtTPCselSPDany->Fill(ntracklets,pttrack);
	if(tofOK){
	  fHistEtaPhiPtTPCselSPDanyTOFbc->Fill(etatrack,phitrack,pttrack);
	  fHistEtaPhiPtInnerTPCselSPDanyTOFbc->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	  fHistNtrackeltsPtTPCselSPDanyTOFbc->Fill(ntracklets,pttrack);
	}
      }
    }


    if(fUseTOFbcSelection && !tofOK) continue;

    if(pidtr>=0 && pidtr<9) fHistdEdxVsPTPCsel[pidtr]->Fill(ptrackTPC,dedx);
    for(Int_t jsp=0; jsp<9; jsp++){
      fHistnSigmaVsPdEdxTPCsel[jsp]->Fill(ptrackTPC,nSigmaTPC[jsp]);
    }

    fHistTPCchi2PerClusPhiPtTPCsel->Fill(chi2clus,pttrack,phitrack);
    fHistSig1ptCovMatPhiPtTPCsel->Fill(curvrelerr,pttrack,phitrack);
    if(itsRefit){
      if(pidtr>=0 && pidtr<9) fHistdEdxVsPTPCselITSref[pidtr]->Fill(ptrackTPC,dedx);
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
    bool isPion = pid[AliPID::kPion];
    bool isKaon = pid[AliPID::kKaon];

    if (fReadMC && pttrack>0.) {
      fHistPtResidVsPtTPCselAll->Fill(ptgen,(pttrack-ptgen));
      fHistOneOverPtResidVsPtTPCselAll->Fill(pttrack,pttrack*(1./pttrack-invptgen));
      fHistPzResidVsPtTPCselAll->Fill(pttrack,(track->Pz()-pzgen));
      fHistPzResidVsEtaTPCselAll->Fill(etatrack,(track->Pz()-pzgen));
      if (itsRefit){
	fHistPtResidVsPtTPCselITSrefAll->Fill(ptgen,(pttrack-ptgen));
	fHistOneOverPtResidVsPtTPCselITSrefAll->Fill(pttrack,pttrack*(1./pttrack-invptgen));
	fHistPzResidVsPtTPCselITSrefAll->Fill(pttrack,(track->Pz()-pzgen));
	fHistPzResidVsEtaTPCselITSrefAll->Fill(etatrack,(track->Pz()-pzgen));
      }
      for (int iS = 0; iS < AliPID::kSPECIESC; ++iS) {
        if (pid[iS]) {
	  Double_t ptDiff=pttrack*AliPID::ParticleCharge(iS)-ptgen;
	  Double_t oneOverPtDiff=0;
	  if(AliPID::ParticleCharge(iS)>0) oneOverPtDiff=pttrack*AliPID::ParticleCharge(iS)*(1./(pttrack*AliPID::ParticleCharge(iS))-invptgen);
          if (pidtr == iS) {
            fHistPtResidVsPtTPCselGoodHyp[iS]->Fill(pttrack*AliPID::ParticleCharge(iS),ptDiff);
            fHistOneOverPtResidVsPtTPCselGoodHyp[iS]->Fill(pttrack*AliPID::ParticleCharge(iS),oneOverPtDiff);
            if (itsRefit){
	      fHistPtResidVsPtTPCselITSrefGoodHyp[iS]->Fill(pttrack*AliPID::ParticleCharge(iS),ptDiff);
	      fHistOneOverPtResidVsPtTPCselITSrefGoodHyp[iS]->Fill(pttrack*AliPID::ParticleCharge(iS),oneOverPtDiff);
	    }
          } else {
            fHistPtResidVsPtTPCselBadHyp[iS]->Fill(pttrack*AliPID::ParticleCharge(iS),ptDiff);
            fHistOneOverPtResidVsPtTPCselBadHyp[iS]->Fill(pttrack*AliPID::ParticleCharge(iS),oneOverPtDiff);
            if (itsRefit){ 
	      fHistPtResidVsPtTPCselITSrefBadHyp[iS]->Fill(pttrack*AliPID::ParticleCharge(iS),ptDiff);
	      fHistOneOverPtResidVsPtTPCselITSrefBadHyp[iS]->Fill(pttrack*AliPID::ParticleCharge(iS),oneOverPtDiff);
	    }
	  }
	}
      }

      if (trlabel >= 0) {
        fHistEtaPhiPtTPCselITSrefGood->Fill(etatrack,phitrack,pttrack);
        if (itsRefit && spdAny) fHistImpParXYPtMulTPCselSPDanyGood->Fill(pttrack,impactXY*10000.,ncl1);
      } else {
	fHistEtaPhiPtTPCselITSrefFake->Fill(etatrack,phitrack,pttrack);
        if (itsRefit && spdAny) fHistImpParXYPtMulTPCselSPDanyFake->Fill(pttrack,impactXY*10000.,ncl1);
      }

      if(itsRefit && spdAny){
	if(isPhysPrim==1) fHistImpParXYPtMulTPCselSPDanyPrim->Fill(pttrack,impactXY*10000.,ncl1);
	else if(isPhysPrim==0) fHistImpParXYPtMulTPCselSPDanySecDec->Fill(pttrack,impactXY*10000.,ncl1);
	else if(isPhysPrim==-1) fHistImpParXYPtMulTPCselSPDanySecMat->Fill(pttrack,impactXY*10000.,ncl1);
      }

    }

    for(Int_t j=0; j<3; j++){
      if(j==0 && !isPion) continue;
      if(j==1 && !isKaon) continue;
      if(j==2 && !isProton) continue;
      Bool_t goodHyp=kFALSE;
      if(pidtr==j+2) goodHyp=kTRUE;
      if(goodHyp){
	fHistdEdxVsPGoodHyp[j]->Fill(ptrackTPC,dedx);
	fHistEtaPhiPtGoodHypTPCsel[j]->Fill(etatrack,phitrack,pttrack);
	fHistEtaPhiPtInnerGoodHypTPCsel[j]->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	fHistTPCchi2PerClusPhiPtGoodHypTPCsel[j]->Fill(chi2clus,pttrack,phitrack);
	if(itsRefit){
	  fHistEtaPhiPtGoodHypTPCselITSref[j]->Fill(etatrack,phitrack,pttrack);
	  fHistEtaPhiPtInnerGoodHypTPCselITSref[j]->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	  fHistTPCchi2PerClusPhiPtGoodHypTPCselITSref[j]->Fill(chi2clus,pttrack,phitrack);
	  if(spdAny){
	    fHistEtaPhiPtGoodHypTPCselSPDany[j]->Fill(etatrack,phitrack,pttrack);
	    fHistEtaPhiPtInnerGoodHypTPCselSPDany[j]->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	    fHistTPCchi2PerClusPhiPtGoodHypTPCselSPDany[j]->Fill(chi2clus,pttrack,phitrack);
	    fHistImpParXYPtMulGoodHypTPCselSPDany[j]->Fill(pttrack,impactXY*10000.,ncl1); 
	  }
	}
      }else{
	fHistdEdxVsPBadHyp[j]->Fill(ptrackTPC,dedx);
	fHistEtaPhiPtBadHypTPCsel[j]->Fill(etatrack,phitrack,pttrack);
	fHistEtaPhiPtInnerBadHypTPCsel[j]->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	fHistTPCchi2PerClusPhiPtBadHypTPCsel[j]->Fill(chi2clus,pttrack,phitrack);
	if(itsRefit){
	  fHistEtaPhiPtBadHypTPCselITSref[j]->Fill(etatrack,phitrack,pttrack);
	  fHistEtaPhiPtInnerBadHypTPCselITSref[j]->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	  fHistTPCchi2PerClusPhiPtBadHypTPCselITSref[j]->Fill(chi2clus,pttrack,phitrack);
	  if(spdAny){
	    fHistEtaPhiPtBadHypTPCselSPDany[j]->Fill(etatrack,phitrack,pttrack);
	    fHistEtaPhiPtInnerBadHypTPCselSPDany[j]->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	    fHistTPCchi2PerClusPhiPtBadHypTPCselSPDany[j]->Fill(chi2clus,pttrack,phitrack);
	    fHistImpParXYPtMulBadHypTPCselSPDany[j]->Fill(pttrack,impactXY*10000.,ncl1);
	  }
	}	

      }
    }
  }

  Int_t nv0s = esd->GetNumberOfV0s();
  for (Int_t iV0 = 0; iV0 < nv0s; iV0++){
    AliESDv0 *v0 = esd->GetV0(iV0);
    if (!v0) continue;
    Bool_t onFlyStatus=v0->GetOnFlyStatus();
    if(onFlyStatus==kTRUE) continue;
    
    UInt_t labPos = (UInt_t)TMath::Abs(v0->GetPindex());
    UInt_t labNeg = (UInt_t)TMath::Abs(v0->GetNindex());
    AliESDtrack *pTrack=esd->GetTrack(labPos);
    AliESDtrack *nTrack=esd->GetTrack(labNeg);
    if (!pTrack || !nTrack) {
      Printf("ERROR: Could not retreive one of the daughter track");
      continue;
    }
    v0->ChangeMassHypothesis(310);
    Double_t invMassK0s = v0->GetEffMass();
    v0->ChangeMassHypothesis(3122);
    Double_t invMassLambda = v0->GetEffMass();
    v0->ChangeMassHypothesis(-3122);
    Double_t invMassAntiLambda = v0->GetEffMass();
    Double_t ptv0=v0->Pt();
    Double_t xv0=v0->Xv();
    Double_t yv0=v0->Yv();
    Double_t rv0=TMath::Sqrt(xv0*xv0+yv0*yv0);

    if(!fTrCutsTPC->AcceptTrack(pTrack)) continue;
    if(!fTrCutsTPC->AcceptTrack(nTrack)) continue;
    
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
      Int_t labbelNeg =TMath::Abs(nTrack->GetLabel());
      TParticle* partPos=mcEvent->Particle(labelPos);
      TParticle* partNeg=mcEvent->Particle(labbelNeg);
      if(partPos && partNeg){
        Int_t labelMotherPos=partPos->GetFirstMother() ;
        Int_t labelMotherNeg=partNeg->GetFirstMother();
	if(labelMotherPos==labelMotherNeg && labelMotherPos>-1){
	  TParticle* partV0=mcEvent->Particle(labelMotherPos);
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

    if(keepK0s){
      fHistInvMassK0s->Fill(invMassK0s,ptv0,rv0);
    }

    if(keepLambda){
      fHistInvMassLambda->Fill(invMassLambda,ptv0,rv0);
      if(pTrack->GetPIDForTracking()==AliPID::kProton){
	fHistInvMassLambdaGoodHyp->Fill(invMassLambda,ptv0,pTrack->GetInnerParam()->Pt());
      }else{
	fHistInvMassLambdaBadHyp->Fill(invMassLambda,ptv0,pTrack->GetInnerParam()->Pt());
      }
    }
    if(keepAntiLambda){
      fHistInvMassAntiLambda->Fill(invMassAntiLambda,ptv0,rv0);
      if(nTrack->GetPIDForTracking()==AliPID::kProton){
	fHistInvMassAntiLambdaGoodHyp->Fill(invMassAntiLambda,ptv0,nTrack->GetInnerParam()->Pt());
      }else{
	fHistInvMassAntiLambdaBadHyp->Fill(invMassAntiLambda,ptv0,nTrack->GetInnerParam()->Pt());
      }
    }
  }
  PostData(1,fOutput);
  PostData(2,fTrackTree);
  
}
//______________________________________________________________________________
void AliAnalysisTaskCheckESDTracks::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  printf("AliAnalysisTaskCheckESDTracks::Terminate --- Number of events: read = %.0f  analysed = %.0f\n",fHistNEvents->GetBinContent(1),fHistNEvents->GetBinContent(5));
  return;
}





