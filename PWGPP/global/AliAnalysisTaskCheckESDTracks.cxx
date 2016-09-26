#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliPID.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
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
  AliAnalysisTaskSE("ITSsa resolution"), 
  fOutput(0x0),
  fHistNEvents(0x0),
  fHistNTracks(0x0),
  fHistEtaPhiPtTPCsel(0x0),
  fHistEtaPhiPtTPCselITSref(0x0),
  fHistEtaPhiPtTPCselSPDany(0x0),
  fHistEtaPhiPtTPCselTOFbc(0x0),
  fHistEtaPhiPtTPCselITSrefTOFbc(0x0),
  fHistEtaPhiPtTPCselSPDanyTOFbc(0x0),
  fHistEtaPhiPtInnerTPCsel(0x0),
  fHistEtaPhiPtInnerTPCselITSref(0x0),
  fHistEtaPhiPtInnerTPCselSPDany(0x0),
  fHistEtaPhiPtInnerTPCselTOFbc(0x0),
  fHistEtaPhiPtInnerTPCselITSrefTOFbc(0x0),
  fHistEtaPhiPtInnerTPCselSPDanyTOFbc(0x0),
  fHistNtrackeltsPtTPCsel(0x0),
  fHistNtrackeltsPtTPCselITSref(0x0),
  fHistNtrackeltsPtTPCselSPDany(0x0),
  fHistNtrackeltsPtTPCselTOFbc(0x0),
  fHistNtrackeltsPtTPCselITSrefTOFbc(0x0),
  fHistNtrackeltsPtTPCselSPDanyTOFbc(0x0),
  fHistTPCchi2PerClusPhiPtTPCsel(0x0),
  fHistTPCchi2PerClusPhiPtTPCselITSref(0x0),
  fHistTPCchi2PerClusPhiPtTPCselSPDany(0x0),
  fHistEtaPhiPtGoodHypProtTPCsel(0x0),
  fHistEtaPhiPtGoodHypProtTPCselITSref(0x0),
  fHistEtaPhiPtGoodHypProtTPCselSPDany(0x0),
  fHistEtaPhiPtInnerGoodHypProtTPCsel(0x0),
  fHistEtaPhiPtInnerGoodHypProtTPCselITSref(0x0),
  fHistEtaPhiPtInnerGoodHypProtTPCselSPDany(0x0),
  fHistTPCchi2PerClusPhiPtGoodHypProtTPCsel(0x0),
  fHistTPCchi2PerClusPhiPtGoodHypProtTPCselITSref(0x0),
  fHistTPCchi2PerClusPhiPtGoodHypProtTPCselSPDany(0x0),
  fHistdEdxVsPGoodHypProt(0x0),
  fHistEtaPhiPtBadHypProtTPCsel(0x0),
  fHistEtaPhiPtBadHypProtTPCselITSref(0x0),
  fHistEtaPhiPtBadHypProtTPCselSPDany(0x0),
  fHistEtaPhiPtInnerBadHypProtTPCsel(0x0),
  fHistEtaPhiPtInnerBadHypProtTPCselITSref(0x0),
  fHistEtaPhiPtInnerBadHypProtTPCselSPDany(0x0),
  fHistTPCchi2PerClusPhiPtBadHypProtTPCsel(0x0),
  fHistTPCchi2PerClusPhiPtBadHypProtTPCselITSref(0x0),
  fHistTPCchi2PerClusPhiPtBadHypProtTPCselSPDany(0x0),
  fHistdEdxVsPBadHypProt(0x0),
  fHistImpParXYPtMulGoodHypPionTPCselSPDany(0x0),
  fHistImpParXYPtMulBadHypPionTPCselSPDany(0x0),
  fHistImpParXYPtMulGoodHypProtTPCselSPDany(0x0),
  fHistImpParXYPtMulBadHypProtTPCselSPDany(0x0),
  fHistPtResidVsPtTPCselAll(0x0),
  fHistPtResidVsPtTPCselGoodHypPion(0x0),
  fHistPtResidVsPtTPCselGoodHypProton(0x0),
  fHistPtResidVsPtTPCselBadHypPion(0x0),
  fHistPtResidVsPtTPCselBadHypProton(0x0),
  fHistPtResidVsPtTPCselITSrefAll(0x0),
  fHistPtResidVsPtTPCselITSrefGoodHypPion(0x0),
  fHistPtResidVsPtTPCselITSrefGoodHypProton(0x0),
  fHistPtResidVsPtTPCselITSrefBadHypPion(0x0),
  fHistPtResidVsPtTPCselITSrefBadHypProton(0x0),
  fHistEtaPhiPtTPCselITSrefGood(0x0),
  fHistEtaPhiPtTPCselITSrefFake(0x0),
  fHistImpParXYPtMulTPCselSPDanyGood(0x0),
  fHistImpParXYPtMulTPCselSPDanyFake(0x0),
  fHistInvMassK0s(0x0),
  fHistInvMassLambdaGoodHyp(0x0),
  fHistInvMassAntiLambdaGoodHyp(0x0),
  fHistInvMassLambdaBadHyp(0x0),
  fHistInvMassAntiLambdaBadHyp(0x0),
  fFillTree(kFALSE),
  fTrackTree(0x0),
  fTreeVarFloat(0x0),
  fTreeVarInt(0x0),
  fTrCutsTPC(0x0),
  fMinNumOfTPCPIDclu(0),
  fUseTOFbcSelection(kTRUE),
  fUsePhysSel(kTRUE),
  fTriggerMask(AliVEvent::kAnyINT),
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
}


//___________________________________________________________________________
AliAnalysisTaskCheckESDTracks::~AliAnalysisTaskCheckESDTracks(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistNTracks;
    delete fHistEtaPhiPtTPCsel;
    delete fHistEtaPhiPtTPCselITSref;
    delete fHistEtaPhiPtTPCselSPDany;
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
    delete fHistEtaPhiPtGoodHypProtTPCsel;
    delete fHistEtaPhiPtGoodHypProtTPCselITSref;
    delete fHistEtaPhiPtGoodHypProtTPCselSPDany;
    delete fHistEtaPhiPtInnerGoodHypProtTPCsel;
    delete fHistEtaPhiPtInnerGoodHypProtTPCselITSref;
    delete fHistEtaPhiPtInnerGoodHypProtTPCselSPDany;
    delete fHistTPCchi2PerClusPhiPtGoodHypProtTPCsel;
    delete fHistTPCchi2PerClusPhiPtGoodHypProtTPCselITSref;
    delete fHistTPCchi2PerClusPhiPtGoodHypProtTPCselSPDany;
    delete fHistdEdxVsPGoodHypProt;
    delete fHistEtaPhiPtBadHypProtTPCsel;
    delete fHistEtaPhiPtBadHypProtTPCselITSref;
    delete fHistEtaPhiPtBadHypProtTPCselSPDany;
    delete fHistEtaPhiPtInnerBadHypProtTPCsel;
    delete fHistEtaPhiPtInnerBadHypProtTPCselITSref;
    delete fHistEtaPhiPtInnerBadHypProtTPCselSPDany;
    delete fHistTPCchi2PerClusPhiPtBadHypProtTPCsel;
    delete fHistTPCchi2PerClusPhiPtBadHypProtTPCselITSref;
    delete fHistTPCchi2PerClusPhiPtBadHypProtTPCselSPDany;
    delete fHistdEdxVsPBadHypProt;
    delete fHistImpParXYPtMulGoodHypPionTPCselSPDany;
    delete fHistImpParXYPtMulBadHypPionTPCselSPDany;
    delete fHistImpParXYPtMulGoodHypProtTPCselSPDany;
    delete fHistImpParXYPtMulBadHypProtTPCselSPDany;
    delete fHistPtResidVsPtTPCselAll;
    delete fHistPtResidVsPtTPCselGoodHypPion;
    delete fHistPtResidVsPtTPCselGoodHypProton;
    delete fHistPtResidVsPtTPCselBadHypPion;
    delete fHistPtResidVsPtTPCselBadHypProton;
    delete fHistPtResidVsPtTPCselITSrefAll;
    delete fHistPtResidVsPtTPCselITSrefGoodHypPion;
    delete fHistPtResidVsPtTPCselITSrefGoodHypProton;
    delete fHistPtResidVsPtTPCselITSrefBadHypPion;
    delete fHistPtResidVsPtTPCselITSrefBadHypProton;
    delete fHistEtaPhiPtTPCselITSrefGood;
    delete fHistEtaPhiPtTPCselITSrefFake;
    delete fHistImpParXYPtMulTPCselSPDanyGood;
    delete fHistImpParXYPtMulTPCselSPDanyFake;
    delete fHistInvMassK0s;
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
  intVarName[2]="ITSrefit";
  intVarName[3]="ITSclumap";
  intVarName[4]="nTPCclu";
  intVarName[5]="TOFbc";
  intVarName[6]="trackPIDhip";
  intVarName[7]="label";
  intVarName[8]="truePID";
  usedVar=kNumOfIntVar-2;
  if(fReadMC) usedVar=kNumOfIntVar;
  for(Int_t ivar=0; ivar<usedVar; ivar++){
    fTrackTree->Branch(intVarName[ivar].Data(),&fTreeVarInt[ivar],Form("%s/I",intVarName[ivar].Data()));
  }
  fOutput->Add(fTrackTree);

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",15,-0.5,14.5);
  //fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"PhysSel"); 
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Good vertex"); 
  fHistNEvents->GetXaxis()->SetBinLabel(4,"|zvert|<10"); 
  fOutput->Add(fHistNEvents);

  fHistNTracks = new TH1F("hNTracks", "Number of tracks in ESD events ; N_{tracks}",5001,-0.5,5000.5);
  fOutput->Add(fHistNTracks);
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

  fHistEtaPhiPtTPCsel = new TH3F("hEtaPhiPtTPCsel"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtTPCselITSref = new TH3F("hEtaPhiPtTPCselITSref"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtTPCselSPDany = new TH3F("hEtaPhiPtTPCselSPDany"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fOutput->Add(fHistEtaPhiPtTPCsel);
  fOutput->Add(fHistEtaPhiPtTPCselITSref);
  fOutput->Add(fHistEtaPhiPtTPCselSPDany);

  fHistEtaPhiPtTPCselTOFbc = new TH3F("hEtaPhiPtTPCselTOFbc"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtTPCselITSrefTOFbc = new TH3F("hEtaPhiPtTPCselITSrefTOFbc"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtTPCselSPDanyTOFbc = new TH3F("hEtaPhiPtTPCselSPDanyTOFbc"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fOutput->Add(fHistEtaPhiPtTPCselTOFbc);
  fOutput->Add(fHistEtaPhiPtTPCselITSrefTOFbc);
  fOutput->Add(fHistEtaPhiPtTPCselSPDanyTOFbc);

  fHistEtaPhiPtInnerTPCsel = new TH3F("hEtaPhiPtInnerTPCsel"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerTPCselITSref = new TH3F("hEtaPhiPtInnerTPCselITSref"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerTPCselSPDany = new TH3F("hEtaPhiPtInnerTPCselSPDany"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fOutput->Add(fHistEtaPhiPtInnerTPCsel);
  fOutput->Add(fHistEtaPhiPtInnerTPCselITSref);
  fOutput->Add(fHistEtaPhiPtInnerTPCselSPDany);

  fHistEtaPhiPtInnerTPCselTOFbc = new TH3F("hEtaPhiPtInnerTPCselTOFbc"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerTPCselITSrefTOFbc = new TH3F("hEtaPhiPtInnerTPCselITSrefTOFbc"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerTPCselSPDanyTOFbc = new TH3F("hEtaPhiPtInnerTPCselSPDanyTOFbc"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fOutput->Add(fHistEtaPhiPtInnerTPCselTOFbc);
  fOutput->Add(fHistEtaPhiPtInnerTPCselITSrefTOFbc);
  fOutput->Add(fHistEtaPhiPtInnerTPCselSPDanyTOFbc);

  fHistNtrackeltsPtTPCsel = new TH2F("hNtrackeltsPtTPCsel"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,40,0.,4.);
  fHistNtrackeltsPtTPCselITSref = new TH2F("hNtrackeltsPtTPCselITSref"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,40,0.,4.);
  fHistNtrackeltsPtTPCselSPDany = new TH2F("hNtrackeltsPtTPCselSPDany"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,40,0.,4.);
  fHistNtrackeltsPtTPCselTOFbc = new TH2F("hNtrackeltsPtTPCselTOFbc"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,40,0.,4.);
  fHistNtrackeltsPtTPCselITSrefTOFbc = new TH2F("hNtrackeltsPtTPCselITSrefTOFbc"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,40,0.,4.);
  fHistNtrackeltsPtTPCselSPDanyTOFbc = new TH2F("hNtrackeltsPtTPCselSPDanyTOFbc"," ; N_{tracklets} ;  p_{T} (GeV/c)",200,0.,10000.,40,0.,4.);
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
  
  fHistEtaPhiPtGoodHypProtTPCsel = new TH3F("hEtaPhiPtGoodHypProtTPCsel"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtGoodHypProtTPCselITSref = new TH3F("hEtaPhiPtGoodHypProtTPCselITSref"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtGoodHypProtTPCselSPDany = new TH3F("hEtaPhiPtGoodHypProtTPCselSPDany"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerGoodHypProtTPCsel = new TH3F("hEtaPhiPtInnerGoodHypProtTPCsel"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerGoodHypProtTPCselITSref = new TH3F("hEtaPhiPtInnerGoodHypProtTPCselITSref"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerGoodHypProtTPCselSPDany = new TH3F("hEtaPhiPtInnerGoodHypProtTPCselSPDany"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistTPCchi2PerClusPhiPtGoodHypProtTPCsel = new TH3F("hTPCchi2PerClusPhiPtGoodHypProtTPCsel"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistTPCchi2PerClusPhiPtGoodHypProtTPCselITSref = new TH3F("hTPCchi2PerClusPhiPtGoodHypProtTPCselITSref"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistTPCchi2PerClusPhiPtGoodHypProtTPCselSPDany = new TH3F("hTPCchi2PerClusPhiPtGoodHypProtTPCselSPDany"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistdEdxVsPGoodHypProt = new TH2F("hdEdxVsPGoodHypProt","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
  fOutput->Add(fHistEtaPhiPtGoodHypProtTPCsel);
  fOutput->Add(fHistEtaPhiPtGoodHypProtTPCselITSref);
  fOutput->Add(fHistEtaPhiPtGoodHypProtTPCselSPDany);
  fOutput->Add(fHistEtaPhiPtInnerGoodHypProtTPCsel);
  fOutput->Add(fHistEtaPhiPtInnerGoodHypProtTPCselITSref);
  fOutput->Add(fHistEtaPhiPtInnerGoodHypProtTPCselSPDany);
  fOutput->Add(fHistTPCchi2PerClusPhiPtGoodHypProtTPCsel);
  fOutput->Add(fHistTPCchi2PerClusPhiPtGoodHypProtTPCselITSref);
  fOutput->Add(fHistTPCchi2PerClusPhiPtGoodHypProtTPCselSPDany);
  fOutput->Add(fHistdEdxVsPGoodHypProt);

  fHistEtaPhiPtBadHypProtTPCsel = new TH3F("hEtaPhiPtBadHypProtTPCsel"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtBadHypProtTPCselITSref = new TH3F("hEtaPhiPtBadHypProtTPCselITSref"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtBadHypProtTPCselSPDany = new TH3F("hEtaPhiPtBadHypProtTPCselSPDany"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerBadHypProtTPCsel = new TH3F("hEtaPhiPtInnerBadHypProtTPCsel"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerBadHypProtTPCselITSref = new TH3F("hEtaPhiPtInnerBadHypProtTPCselITSref"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtInnerBadHypProtTPCselSPDany = new TH3F("hEtaPhiPtInnerBadHypProtTPCselSPDany"," ; #eta_{TPC} ; #varphi_{TPC} ; p_{T,TPC} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistTPCchi2PerClusPhiPtBadHypProtTPCsel = new TH3F("hTPCchi2PerClusPhiPtBadHypProtTPCsel"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistTPCchi2PerClusPhiPtBadHypProtTPCselITSref = new TH3F("hTPCchi2PerClusPhiPtBadHypProtTPCselITSref"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistTPCchi2PerClusPhiPtBadHypProtTPCselSPDany = new TH3F("hTPCchi2PerClusPhiPtBadHypProtTPCselSPDany"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistdEdxVsPBadHypProt = new TH2F("hdEdxVsPBadHypProt","  ; p_{TPC} (GeV/c) ; dE/dx",100,0.,5.,100,0.,600.);
  fOutput->Add(fHistEtaPhiPtBadHypProtTPCsel);
  fOutput->Add(fHistEtaPhiPtBadHypProtTPCselITSref);
  fOutput->Add(fHistEtaPhiPtBadHypProtTPCselSPDany);
  fOutput->Add(fHistEtaPhiPtInnerBadHypProtTPCsel);
  fOutput->Add(fHistEtaPhiPtInnerBadHypProtTPCselITSref);
  fOutput->Add(fHistEtaPhiPtInnerBadHypProtTPCselSPDany);
  fOutput->Add(fHistTPCchi2PerClusPhiPtBadHypProtTPCsel);
  fOutput->Add(fHistTPCchi2PerClusPhiPtBadHypProtTPCselITSref);
  fOutput->Add(fHistTPCchi2PerClusPhiPtBadHypProtTPCselSPDany);
  fOutput->Add(fHistdEdxVsPBadHypProt);

  fHistImpParXYPtMulGoodHypPionTPCselSPDany = new TH3F("hImpParXYPtMulGoodHypPionTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fHistImpParXYPtMulBadHypPionTPCselSPDany = new TH3F("hImpParXYPtMulBadHypPionTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fHistImpParXYPtMulGoodHypProtTPCselSPDany = new TH3F("hImpParXYPtMulGoodHypProtTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fHistImpParXYPtMulBadHypProtTPCselSPDany = new TH3F("hImpParXYPtMulBadHypProtTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fOutput->Add(fHistImpParXYPtMulGoodHypPionTPCselSPDany);
  fOutput->Add(fHistImpParXYPtMulBadHypPionTPCselSPDany);
  fOutput->Add(fHistImpParXYPtMulGoodHypProtTPCselSPDany);
  fOutput->Add(fHistImpParXYPtMulBadHypProtTPCselSPDany);

  fHistPtResidVsPtTPCselAll = new TH2F("hPtResidVsPtTPCselAll"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselGoodHypPion = new TH2F("hPtResidVsPtTPCselGoodHypPion"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselGoodHypProton = new TH2F("hPtResidVsPtTPCselGoodHypProton"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselBadHypPion = new TH2F("hPtResidVsPtTPCselBadHypPion"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselBadHypProton = new TH2F("hPtResidVsPtTPCselBadHypProton"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefAll = new TH2F("hPtResidVsPtTPCselITSrefAll"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefGoodHypPion = new TH2F("hPtResidVsPtTPCselITSrefGoodHypPion"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefGoodHypProton = new TH2F("hPtResidVsPtTPCselITSrefGoodHypProton"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefBadHypPion = new TH2F("hPtResidVsPtTPCselITSrefBadHypPion"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefBadHypProton = new TH2F("hPtResidVsPtTPCselITSrefBadHypProton"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fOutput->Add(fHistPtResidVsPtTPCselAll);
  fOutput->Add(fHistPtResidVsPtTPCselGoodHypPion);
  fOutput->Add(fHistPtResidVsPtTPCselGoodHypProton);
  fOutput->Add(fHistPtResidVsPtTPCselBadHypPion);
  fOutput->Add(fHistPtResidVsPtTPCselBadHypProton);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefAll);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefGoodHypPion);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefGoodHypProton);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefBadHypPion);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefBadHypProton);
 
  fHistEtaPhiPtTPCselITSrefGood = new TH3F("hEtaPhiPtTPCselITSrefGood"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtTPCselITSrefFake = new TH3F("hEtaPhiPtTPCselITSrefFake"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fOutput->Add(fHistEtaPhiPtTPCselITSrefGood);
  fOutput->Add(fHistEtaPhiPtTPCselITSrefFake);

  fHistImpParXYPtMulTPCselSPDanyGood = new TH3F("hImpParXYPtMulTPCselSPDanyGood"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fHistImpParXYPtMulTPCselSPDanyFake = new TH3F("hImpParXYPtMulTPCselSPDanyFake"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyGood);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyFake);

  fHistInvMassK0s = new TH2F("hInvMassK0s"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(K0s) ",200,0.4,0.6,25,0.,5.);
  fHistInvMassLambdaGoodHyp = new TH3F("hInvMassLambdaGoodHyp"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#Lambda) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fHistInvMassAntiLambdaGoodHyp = new TH3F("hInvMassAntiLambdaGoodHyp"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#bar{#Lambda}) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fHistInvMassLambdaBadHyp = new TH3F("hInvMassLambdaBadHyp"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#Lambda) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fHistInvMassAntiLambdaBadHyp = new TH3F("hInvMassAntiLambdaBadHyp"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#bar{#Lambda}) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fOutput->Add(fHistInvMassK0s);
  fOutput->Add(fHistInvMassLambdaGoodHyp);
  fOutput->Add(fHistInvMassAntiLambdaGoodHyp);
  fOutput->Add(fHistInvMassLambdaBadHyp);
  fOutput->Add(fHistInvMassAntiLambdaBadHyp);

  PostData(1,fOutput);

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
  
  AliStack* stack=0;

  if(fReadMC){
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    stack = mcEvent->Stack();
    if (!stack) {
      Printf("ERROR: stack not available");
      return;
    }
  }


  fHistNEvents->Fill(0);
  if(fUsePhysSel){
    Bool_t isPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
    if(!isPhysSel) return;
  }
  fHistNEvents->Fill(1);

  const AliESDVertex *spdv=esd->GetPrimaryVertexSPD();
  if(spdv->GetNContributors()<=0) return;
  fHistNEvents->Fill(2);

  Float_t xvert=spdv->GetX();
  Float_t yvert=spdv->GetY();
  Float_t zvert=spdv->GetZ();

  if(TMath::Abs(zvert)>10) return;
  fHistNEvents->Fill(3);

  esd->InitMagneticField();

  const Int_t ntSize=33;
  Float_t xnt[ntSize];
  
  Int_t ntracks = esd->GetNumberOfTracks();
  Int_t ntracklets = 0;
  Int_t ncl1 = 0;
  const AliMultiplicity* mult=esd->GetMultiplicity();
  if(mult){
    ntracklets = mult->GetNumberOfTracklets();
    ncl1 = mult->GetNumberOfITSClusters(1);
  }
  Int_t nPureSAtracks=0;
  Int_t nITSTPCtracks=0;
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

    Double_t pttrack=track->Pt();
    Double_t ptrack=track->P();
    Double_t pxtrack=track->Px();
    Double_t pytrack=track->Py();
    Double_t pztrack=track->Pz();
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
    Double_t chi2clus=track->GetTPCchi2()/track->GetTPCNcls();
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

    fTreeVarInt[2]=itsRefit;
    fTreeVarInt[3]=clumap;
    fTreeVarInt[4]=nTPCclus;
    fTreeVarInt[5]=tofOK;

    Int_t trlabel=track->GetLabel();
    Float_t invpttrack=track->OneOverPt();
    Float_t dedx=track->GetTPCsignal();
    Int_t  pidtr=track->GetPIDForTracking();
    
    Double_t nSigmaTPC[9]={-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.};
    AliPIDResponse::EDetPidStatus status = pidResp->CheckPIDStatus(AliPIDResponse::kTPC,track);
    if (status == AliPIDResponse::kDetPidOk){
      for(Int_t jsp=0; jsp<9; jsp++){
	nSigmaTPC[jsp]=pidResp->NumberOfSigmasTPC(track,(AliPID::EParticleType)jsp);
      }
    }
    fTreeVarFloat[22]=dedx;
    fTreeVarFloat[23]=nSigmaTPC[0];
    fTreeVarFloat[24]=nSigmaTPC[2];
    fTreeVarFloat[25]=nSigmaTPC[3];
    fTreeVarFloat[26]=nSigmaTPC[4];
    fTreeVarInt[6]=pidtr;
    
    Float_t ptgen=-999.;
    Float_t pgen=-999.;
    Float_t pxgen=-999.;
    Float_t pygen=-999.;
    Float_t pzgen=-999.;
    Float_t etagen=-999.;
    Float_t phigen=-999.;
    Int_t hadronSpecie=-1;
    Float_t invptgen=-999.;
    if(fReadMC){
      TParticle* part = stack->Particle(TMath::Abs(trlabel));
      ptgen=part->Pt();
      pgen=part->P();
      pxgen=part->Px();
      pygen=part->Py();
      pzgen=part->Pz();
      if(ptgen>0.) invptgen=1./ptgen;
      etagen=part->Eta();
      phigen=part->Phi();
      fTreeVarFloat[27]=pxgen;
      fTreeVarFloat[28]=pygen;
      fTreeVarFloat[29]=pzgen;
      fTreeVarFloat[30]=ptgen;
      fTreeVarFloat[31]=pgen;
      fTreeVarFloat[32]=etagen;
      fTreeVarFloat[33]=phigen;
      if(fUseMCId){
	Int_t pdg=TMath::Abs(part->GetPdgCode());
	if(pdg==11) hadronSpecie=0;
	else if(pdg==13) hadronSpecie=1;
	else if(pdg==211) hadronSpecie=2;
	else if(pdg==321) hadronSpecie=3;
	else if(pdg==2212) hadronSpecie=4;
      }
      fTreeVarInt[7]=trlabel;
      fTreeVarInt[8]=part->GetPdgCode();
    }

    if(pidtr>=0 && pidtr<9) fHistdEdxVsP[pidtr]->Fill(ptrackTPC,dedx);
    
    if(!fTrCutsTPC->AcceptTrack(track)) continue;
    if(track->GetTPCsignalN()<fMinNumOfTPCPIDclu) continue;

    if (fFillTree) fTrackTree->Fill();    

    fHistEtaPhiPtTPCsel->Fill(etatrack,phitrack,pttrack);
    fHistEtaPhiPtInnerTPCsel->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
    fHistNtrackeltsPtTPCsel->Fill(ntracklets,pttrack);
    if(tofOK){
      fHistEtaPhiPtTPCselTOFbc->Fill(etatrack,phitrack,pttrack);
      fHistEtaPhiPtInnerTPCselTOFbc->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
      fHistNtrackeltsPtTPCselTOFbc->Fill(ntracklets,pttrack);
    }
    if(itsRefit){
      fHistEtaPhiPtTPCselITSref->Fill(etatrack,phitrack,pttrack);
      fHistEtaPhiPtInnerTPCselITSref->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
      fHistNtrackeltsPtTPCselITSref->Fill(ntracklets,pttrack);
      if(tofOK){
	fHistEtaPhiPtTPCselITSrefTOFbc->Fill(etatrack,phitrack,pttrack);
	fHistEtaPhiPtInnerTPCselITSrefTOFbc->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	fHistNtrackeltsPtTPCselITSrefTOFbc->Fill(ntracklets,pttrack);
      }
      if(spdAny){ 
	fHistEtaPhiPtTPCselSPDany->Fill(etatrack,phitrack,pttrack);
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
    if(itsRefit){
      if(pidtr>=0 && pidtr<9) fHistdEdxVsPTPCselITSref[pidtr]->Fill(ptrackTPC,dedx);
      fHistTPCchi2PerClusPhiPtTPCselITSref->Fill(chi2clus,pttrack,phitrack);
      if(spdAny){ 
	fHistTPCchi2PerClusPhiPtTPCselSPDany->Fill(chi2clus,pttrack,phitrack);
      }
    }

    Bool_t isProton=kFALSE;
    Bool_t isPion=kFALSE;
    if(fReadMC && fUseMCId){
      if(hadronSpecie==4) isProton=kTRUE;
      if(hadronSpecie==2) isPion=kTRUE;
    }else{
      if(TMath::Abs(nSigmaTPC[AliPID::kProton])<3) isProton=kTRUE;
      if(TMath::Abs(nSigmaTPC[AliPID::kPion])<3) isPion=kTRUE;
    }

    if(fReadMC){
      fHistPtResidVsPtTPCselAll->Fill(ptgen,(pttrack-ptgen));
      if(itsRefit) fHistPtResidVsPtTPCselITSrefAll->Fill(ptgen,(pttrack-ptgen));
      if(trlabel>=0){
	fHistEtaPhiPtTPCselITSrefGood->Fill(etatrack,phitrack,pttrack);
	if(itsRefit && spdAny) fHistImpParXYPtMulTPCselSPDanyGood->Fill(pttrack,impactXY*10000.,ncl1);
      }else{
	fHistEtaPhiPtTPCselITSrefFake->Fill(etatrack,phitrack,pttrack);
	if(itsRefit && spdAny) fHistImpParXYPtMulTPCselSPDanyFake->Fill(pttrack,impactXY*10000.,ncl1);
      }
    }

    if(isProton){
      if(pidtr==4){
	fHistdEdxVsPGoodHypProt->Fill(ptrackTPC,dedx);
	fHistEtaPhiPtGoodHypProtTPCsel->Fill(etatrack,phitrack,pttrack);
	fHistEtaPhiPtInnerGoodHypProtTPCsel->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	fHistTPCchi2PerClusPhiPtGoodHypProtTPCsel->Fill(chi2clus,pttrack,phitrack);
	if(itsRefit){
	  fHistEtaPhiPtGoodHypProtTPCselITSref->Fill(etatrack,phitrack,pttrack);
	  fHistEtaPhiPtInnerGoodHypProtTPCselITSref->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	  fHistTPCchi2PerClusPhiPtGoodHypProtTPCselITSref->Fill(chi2clus,pttrack,phitrack);
	  if(spdAny){
	    fHistEtaPhiPtGoodHypProtTPCselSPDany->Fill(etatrack,phitrack,pttrack);
	    fHistEtaPhiPtInnerGoodHypProtTPCselSPDany->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	    fHistTPCchi2PerClusPhiPtGoodHypProtTPCselSPDany->Fill(chi2clus,pttrack,phitrack);
	    fHistImpParXYPtMulGoodHypProtTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1); 
	  }
	}
	if(fReadMC){
	  fHistPtResidVsPtTPCselGoodHypProton->Fill(ptgen,(pttrack-ptgen));
	  if(itsRefit) fHistPtResidVsPtTPCselITSrefGoodHypProton->Fill(ptgen,(pttrack-ptgen));
	}
      }else{
	fHistdEdxVsPBadHypProt->Fill(ptrackTPC,dedx);
	fHistEtaPhiPtBadHypProtTPCsel->Fill(etatrack,phitrack,pttrack);
	fHistEtaPhiPtInnerBadHypProtTPCsel->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	fHistTPCchi2PerClusPhiPtBadHypProtTPCsel->Fill(chi2clus,pttrack,phitrack);
	if(itsRefit){
	  fHistEtaPhiPtBadHypProtTPCselITSref->Fill(etatrack,phitrack,pttrack);
	  fHistEtaPhiPtInnerBadHypProtTPCselITSref->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	  fHistTPCchi2PerClusPhiPtBadHypProtTPCselITSref->Fill(chi2clus,pttrack,phitrack);
	  if(spdAny){
	    fHistEtaPhiPtBadHypProtTPCselSPDany->Fill(etatrack,phitrack,pttrack);
	    fHistEtaPhiPtInnerBadHypProtTPCselSPDany->Fill(etatrackTPC,phitrackTPC,pttrackTPC);
	    fHistTPCchi2PerClusPhiPtBadHypProtTPCselSPDany->Fill(chi2clus,pttrack,phitrack);
	    fHistImpParXYPtMulBadHypProtTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1);
	  }
	}	
	if(fReadMC){
	  fHistPtResidVsPtTPCselBadHypProton->Fill(ptgen,(pttrack-ptgen));
	  if(itsRefit) fHistPtResidVsPtTPCselITSrefBadHypProton->Fill(ptgen,(pttrack-ptgen));
	}
      }
    }

    if(isPion){
      if(pidtr==2){
	if(itsRefit){
	  if(spdAny){
	    fHistImpParXYPtMulGoodHypPionTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1);
	  }
	}
	if(fReadMC){
	  fHistPtResidVsPtTPCselGoodHypPion->Fill(ptgen,(pttrack-ptgen));
	  if(itsRefit) fHistPtResidVsPtTPCselITSrefGoodHypPion->Fill(ptgen,(pttrack-ptgen));
	}
      }else{
	if(itsRefit){
	  if(spdAny){ 
	    fHistImpParXYPtMulBadHypPionTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1);
	  }
	}
	if(fReadMC){
	  fHistPtResidVsPtTPCselBadHypPion->Fill(ptgen,(pttrack-ptgen));
	  if(itsRefit) fHistPtResidVsPtTPCselITSrefBadHypPion->Fill(ptgen,(pttrack-ptgen));
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

    if(!fTrCutsTPC->AcceptTrack(pTrack)) continue;
    if(!fTrCutsTPC->AcceptTrack(nTrack)) continue;
    
    Bool_t keepK0s=kTRUE;
    Bool_t keepLambda=kTRUE;
    Bool_t keepAntiLambda=kTRUE;
    if(!fReadMC){
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
    }else{
      keepK0s=kFALSE;
      keepLambda=kFALSE;
      keepAntiLambda=kFALSE;
      Int_t labelPos=TMath::Abs(pTrack->GetLabel());
      Int_t labbelNeg =TMath::Abs(nTrack->GetLabel());
      TParticle* partPos=stack->Particle(labelPos);
      TParticle* partNeg=stack->Particle(labbelNeg);
      if(partPos && partNeg){
        Int_t labelMotherPos=partPos->GetFirstMother() ;
        Int_t labelMotherNeg=partNeg->GetFirstMother();
	if(labelMotherPos==labelMotherNeg && labelMotherPos>-1){
	  TParticle* partV0=stack->Particle(labelMotherPos);
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

    if(keepLambda){
      if(pTrack->GetPIDForTracking()==AliPID::kProton){
	fHistInvMassLambdaGoodHyp->Fill(invMassLambda,ptv0,pTrack->GetInnerParam()->Pt());
      }else{
	fHistInvMassLambdaBadHyp->Fill(invMassLambda,ptv0,pTrack->GetInnerParam()->Pt());
      }
    }
    if(keepAntiLambda){
      if(nTrack->GetPIDForTracking()==AliPID::kProton){
	fHistInvMassAntiLambdaGoodHyp->Fill(invMassAntiLambda,ptv0,nTrack->GetInnerParam()->Pt());
      }else{
	fHistInvMassAntiLambdaBadHyp->Fill(invMassAntiLambda,ptv0,nTrack->GetInnerParam()->Pt());
      }
    }
  }
  PostData(1,fOutput);
  
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
  printf("AliAnalysisTaskCheckESDTracks::Terminate --- Number of events: read = %.0f  analysed = %.0f\n",fHistNEvents->GetBinContent(1),fHistNEvents->GetBinContent(4));
  return;
}





