#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODEvent.h"
#include "AliStack.h"
#include "AliPID.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAODTracklets.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include <AliAODMCParticle.h>
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
  AliAnalysisTaskSE("ITSsa resolution"), 
  fOutput(0x0),
  fHistNEvents(0x0),
  fHistNTracks(0x0),
  fHistFilterBits(0x0),
  fHistEtaPhiPtTPCsel(0x0),
  fHistEtaPhiPtTPCselITSref(0x0),
  fHistEtaPhiPtTPCselSPDany(0x0),
  fHistEtaPhiPtTPCselTOFbc(0x0),
  fHistEtaPhiPtTPCselITSrefTOFbc(0x0),
  fHistEtaPhiPtTPCselSPDanyTOFbc(0x0),
  fHistTPCchi2PerClusPhiPtTPCsel(0x0),
  fHistTPCchi2PerClusPhiPtTPCselITSref(0x0),
  fHistTPCchi2PerClusPhiPtTPCselSPDany(0x0),
  fHistImpParXYPtMulPionTPCselSPDany(0x0),
  fHistImpParXYPtMulKaonTPCselSPDany(0x0),
  fHistImpParXYPtMulProtonTPCselSPDany(0x0),
  fHistPtResidVsPtTPCselAll(0x0),
  fHistPtResidVsPtTPCselPion(0x0),
  fHistPtResidVsPtTPCselKaon(0x0),
  fHistPtResidVsPtTPCselProton(0x0),
  fHistPtResidVsPtTPCselITSrefAll(0x0),
  fHistPtResidVsPtTPCselITSrefPion(0x0),
  fHistPtResidVsPtTPCselITSrefKaon(0x0),
  fHistPtResidVsPtTPCselITSrefProton(0x0),
  fHistEtaPhiPtTPCselITSrefGood(0x0),
  fHistEtaPhiPtTPCselITSrefFake(0x0), 
  fHistImpParXYPtMulTPCselSPDanyGood(0x0),
  fHistImpParXYPtMulTPCselSPDanyFake(0x0),
  fHistImpParXYPtMulTPCselSPDanyPrim(0x0),
  fHistImpParXYPtMulTPCselSPDanySec(0x0),
  fHistInvMassK0s(0x0),
  fHistInvMassLambda(0x0),
  fHistInvMassAntiLambda(0x0),
  fFillTree(kFALSE),
  fTrackTree(0x0),
  fTreeVarFloat(0x0),
  fTreeVarInt(0x0),
  fTrCutsTPC(0x0),
  fMinNumOfTPCPIDclu(0),
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
  for(Int_t jb=0; jb<kNumOfFilterBits; jb++){
    fHistImpParXYPtMulFiltBit[jb]=0x0;
    fHistEtaPhiPtFiltBit[jb]=0x0;
    fHistITScluPtFiltBit[jb]=0x0;
    fHistSPDcluPtFiltBit[jb]=0x0;
    fHistTPCcluPtFiltBit[jb]=0x0;
    fHistTPCcrrowsPtFiltBit[jb]=0x0;
    fHistTPCCrowOverFindPtFiltBit[jb]=0x0;
    fHistTPCChi2ndfPtFiltBit[jb]=0x0;
    fHistChi2TPCConstrVsGlobPtFiltBit[jb]=0x0;
  }
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskCheckAODTracks::~AliAnalysisTaskCheckAODTracks(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistNTracks;
    delete fHistFilterBits;
    delete fHistEtaPhiPtTPCsel;
    delete fHistEtaPhiPtTPCselITSref;
    delete fHistEtaPhiPtTPCselSPDany;
    delete fHistEtaPhiPtTPCselTOFbc;
    delete fHistEtaPhiPtTPCselITSrefTOFbc;
    delete fHistEtaPhiPtTPCselSPDanyTOFbc;
    delete fHistTPCchi2PerClusPhiPtTPCsel;
    delete fHistTPCchi2PerClusPhiPtTPCselITSref;
    delete fHistTPCchi2PerClusPhiPtTPCselSPDany;
    delete fHistImpParXYPtMulPionTPCselSPDany;
    delete fHistImpParXYPtMulKaonTPCselSPDany;
    delete fHistImpParXYPtMulProtonTPCselSPDany;
    delete fHistPtResidVsPtTPCselAll;
    delete fHistPtResidVsPtTPCselPion;
    delete fHistPtResidVsPtTPCselKaon;
    delete fHistPtResidVsPtTPCselProton;
    delete fHistPtResidVsPtTPCselITSrefAll;
    delete fHistPtResidVsPtTPCselITSrefPion;
    delete fHistPtResidVsPtTPCselITSrefKaon;
    delete fHistPtResidVsPtTPCselITSrefProton;
    delete fHistEtaPhiPtTPCselITSrefGood;
    delete fHistEtaPhiPtTPCselITSrefFake;
    delete fHistImpParXYPtMulTPCselSPDanyGood;
    delete fHistImpParXYPtMulTPCselSPDanyFake;
    delete fHistImpParXYPtMulTPCselSPDanyPrim;
    delete fHistImpParXYPtMulTPCselSPDanySec;
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
      delete fHistTPCChi2ndfPtFiltBit[jb];
      delete fHistChi2TPCConstrVsGlobPtFiltBit[jb];
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
  intVarName[2]="ITSrefit";
  intVarName[3]="ITSclumap";
  intVarName[4]="nTPCclu";
  intVarName[5]="TOFbc";
  intVarName[6]="filbits";
  intVarName[7]="charge";
  intVarName[8]="label";
  intVarName[9]="truePID";
  intVarName[10]="isPhysPrim";
  usedVar=kNumOfIntVar-4;
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

  fHistNTracks = new TH1F("hNTracks", "Number of tracks in AOD events ; N_{tracks}",5001,-0.5,5000.5);
  fOutput->Add(fHistNTracks);

  fHistFilterBits = new TH2F("hFilterBits", " ; Filter Bit ; Id ; N_{tracks}",10,-0.5,9.5,2,-1,1);
  fHistFilterBits->GetYaxis()->SetBinLabel(1,"Neg. ID");
  fHistFilterBits->GetYaxis()->SetBinLabel(2,"Pos. ID");
  fOutput->Add(fHistFilterBits);
  

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

  fHistTPCchi2PerClusPhiPtTPCsel = new TH3F("hTPCchi2PerClusPhiPtTPCsel"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistTPCchi2PerClusPhiPtTPCselITSref = new TH3F("hTPCchi2PerClusPhiPtTPCselITSref"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fHistTPCchi2PerClusPhiPtTPCselSPDany = new TH3F("hTPCchi2PerClusPhiPtTPCselSPDany"," ; TPC #chi^{2}/nClusters; p_{T} (GeV/c) ; #varphi",100, 0, 10, 100, 0, 10, 72, 0, 2*TMath::Pi());
  fOutput->Add(fHistTPCchi2PerClusPhiPtTPCsel);
  fOutput->Add(fHistTPCchi2PerClusPhiPtTPCselITSref);
  fOutput->Add(fHistTPCchi2PerClusPhiPtTPCselSPDany);
  
  fHistImpParXYPtMulPionTPCselSPDany = new TH3F("hImpParXYPtMulPionTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fHistImpParXYPtMulKaonTPCselSPDany = new TH3F("hImpParXYPtMulKaonTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fHistImpParXYPtMulProtonTPCselSPDany = new TH3F("hImpParXYPtMulProtonTPCselSPDany"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fOutput->Add(fHistImpParXYPtMulPionTPCselSPDany);
  fOutput->Add(fHistImpParXYPtMulKaonTPCselSPDany);
  fOutput->Add(fHistImpParXYPtMulProtonTPCselSPDany);
  for(Int_t jb=0; jb<kNumOfFilterBits; jb++){
    fHistEtaPhiPtFiltBit[jb] =  new TH3F(Form("hEtaPhiPtFiltBit%d",jb)," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
    fHistImpParXYPtMulFiltBit[jb] = new TH3F(Form("hImpParXYPtMulPionFiltBit%d",jb)," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
    fHistITScluPtFiltBit[jb] = new TH2F(Form("hITScluPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; n ITS clusters",50,0.,10.,7,-0.5,6.5);
    fHistSPDcluPtFiltBit[jb] = new TH2F(Form("hSPDcluPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; SPD clusters",50,0.,10.,4,-0.5,3.5);
    fHistSPDcluPtFiltBit[jb]->GetYaxis()->SetBinLabel(1,"kNone");
    fHistSPDcluPtFiltBit[jb]->GetYaxis()->SetBinLabel(2,"kOnlyFirst");
    fHistSPDcluPtFiltBit[jb]->GetYaxis()->SetBinLabel(3,"kOnlySecond");
    fHistSPDcluPtFiltBit[jb]->GetYaxis()->SetBinLabel(4,"kBoth");
    fHistTPCcluPtFiltBit[jb] = new TH2F(Form("hTPCcluPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; n TPC clusters",50,0.,10.,161,-0.5,160.5);  
    fHistTPCcrrowsPtFiltBit[jb] = new TH2F(Form("hTPCcrrowsPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; n TPC Crossed Rows",50,0.,10.,161,-0.5,160.5);  
    fHistTPCCrowOverFindPtFiltBit[jb] = new TH2F(Form("hTPCCrowOverFindPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; #chi^{2}/ndf",50,0.,10.,100,0.,2.);    
    fHistTPCChi2ndfPtFiltBit[jb] = new TH2F(Form("hTPCChi2ndfPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; #chi^{2}/ndf",50,0.,10.,160,0.,8.);
    fHistChi2TPCConstrVsGlobPtFiltBit[jb] = new TH2F(Form("hChi2TPCConstrVsGlobPtFiltBit%d",jb)," ; p_{T} (GeV/c) ; golden #chi^{2}",50,0.,10.,160,0.,8.);
    fOutput->Add(fHistEtaPhiPtFiltBit[jb]);
    fOutput->Add(fHistImpParXYPtMulFiltBit[jb]);
    fOutput->Add(fHistITScluPtFiltBit[jb]);
    fOutput->Add(fHistSPDcluPtFiltBit[jb]);
    fOutput->Add(fHistTPCcluPtFiltBit[jb]);
    fOutput->Add(fHistTPCcrrowsPtFiltBit[jb]);
    fOutput->Add(fHistTPCChi2ndfPtFiltBit[jb]);
    fOutput->Add(fHistChi2TPCConstrVsGlobPtFiltBit[jb]);
    fOutput->Add(fHistTPCCrowOverFindPtFiltBit[jb]);
  }

  fHistPtResidVsPtTPCselAll = new TH2F("hPtResidVsPtTPCselAll"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselPion = new TH2F("hPtResidVsPtTPCselPion"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselKaon = new TH2F("hPtResidVsPtTPCselKaon"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselProton = new TH2F("hPtResidVsPtTPCselProton"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefAll = new TH2F("hPtResidVsPtTPCselITSrefAll"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefPion = new TH2F("hPtResidVsPtTPCselITSrefPion"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefKaon = new TH2F("hPtResidVsPtTPCselITSrefKaon"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fHistPtResidVsPtTPCselITSrefProton = new TH2F("hPtResidVsPtTPCselITSrefProton"," ; p_{T,gen} (GeV/c) ; p_{T,reco}-p_{T,gen} (GeV/c)",40,0.,4.,100,-0.5,0.5);
  fOutput->Add(fHistPtResidVsPtTPCselAll);
  fOutput->Add(fHistPtResidVsPtTPCselPion);
  fOutput->Add(fHistPtResidVsPtTPCselKaon);
  fOutput->Add(fHistPtResidVsPtTPCselProton);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefAll);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefPion);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefKaon);
  fOutput->Add(fHistPtResidVsPtTPCselITSrefProton);
 
  fHistEtaPhiPtTPCselITSrefGood = new TH3F("hEtaPhiPtTPCselITSrefGood"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fHistEtaPhiPtTPCselITSrefFake = new TH3F("hEtaPhiPtTPCselITSrefFake"," ; #eta ; #varphi ; p_{T} (GeV/c)",20,-1.,1.,72,0.,2*TMath::Pi(),40,0.,4.);
  fOutput->Add(fHistEtaPhiPtTPCselITSrefGood);
  fOutput->Add(fHistEtaPhiPtTPCselITSrefFake);

  fHistImpParXYPtMulTPCselSPDanyGood = new TH3F("hImpParXYPtMulTPCselSPDanyGood"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fHistImpParXYPtMulTPCselSPDanyFake = new TH3F("hImpParXYPtMulTPCselSPDanyFake"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyGood);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyFake);

  fHistImpParXYPtMulTPCselSPDanyPrim = new TH3F("hImpParXYPtMulTPCselSPDanyPrim"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fHistImpParXYPtMulTPCselSPDanySec = new TH3F("hImpParXYPtMulTPCselSPDanySec"," ; p_{T} (GeV/c) ; d_{0}^{xy} (#mum) ; N_{CL1}",200,0.,4.,400,-1500.,1500,50,0.,10000.);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanyPrim);
  fOutput->Add(fHistImpParXYPtMulTPCselSPDanySec);


  fHistInvMassK0s = new TH2F("hInvMassK0s"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(K0s) ",200,0.4,0.6,25,0.,5.);
  fHistInvMassLambda = new TH3F("hInvMassLambda"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#Lambda) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fHistInvMassAntiLambda = new TH3F("hInvMassAntiLambda"," ; Inv.Mass (GeV/c^{2}) ; p_{T}(#bar{#Lambda}) ; p_{T,TPC}(p)",200,1.0,1.2,25,0.,5.,50,0.,5.);
  fOutput->Add(fHistInvMassK0s);
  fOutput->Add(fHistInvMassLambda);
  fOutput->Add(fHistInvMassAntiLambda);

  PostData(1,fOutput);

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

  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  if(vtx1->GetNContributors()<=0) return;
  fHistNEvents->Fill(2);

  Float_t xvert=vtx1->GetX();
  Float_t yvert=vtx1->GetY();
  Float_t zvert=vtx1->GetZ();
  if(TMath::Abs(zvert)>10) return;
  fHistNEvents->Fill(3);

  Double_t pos[3],cov[6];
  vtx1->GetXYZ(pos);
  vtx1->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);

  Int_t ntracks = aod->GetNumberOfTracks();
  Int_t ntracklets = 0;
  AliAODTracklets *mult=aod->GetTracklets();
  if(mult) ntracklets=mult->GetNumberOfTracklets();
  Int_t ncl1 = aod->GetNumberOfITSClusters(1);
  fHistNTracks->Fill(ntracks);

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
    
    Double_t d0z0[2],covd0z0[3];
    Bool_t isOK=track->PropagateToDCA(vtx1,magField,99999.,d0z0,covd0z0);
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
    Int_t nSPDclus=0;
    if(track->HasPointOnITSLayer(0)) nSPDclus+=1;
    if(track->HasPointOnITSLayer(1)) nSPDclus+=2;
    Bool_t spdAny=kFALSE;
    if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) spdAny=kTRUE;
    Int_t nTPCclus=track->GetNcls(1);
    Double_t chi2clus=track->Chi2perNDF();
    Double_t goldenChi2=track->GetChi2TPCConstrainedVsGlobal();
    Float_t nCrossedRowsTPC = track->GetTPCCrossedRows();
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (track->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
    }
    Int_t tofBC=track->GetTOFBunchCrossing(magField);
    fTreeVarFloat[10]=impactXY;
    fTreeVarFloat[11]=impactZ;
    fTreeVarFloat[12]=chi2clus;
    fTreeVarFloat[13]=goldenChi2;
    fTreeVarFloat[14]=ratioCrossedRowsOverFindableClustersTPC;

    fTreeVarInt[2]=itsRefit;
    fTreeVarInt[3]=clumap;
    fTreeVarInt[4]=nTPCclus;
    fTreeVarInt[5]=tofBC;

    Int_t trlabel=track->GetLabel();
    Float_t dedx=track->GetTPCsignal();
    Int_t  filtmap=track->GetFilterMap();
    
    Double_t nSigmaTPC[9]={-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.};
    AliPIDResponse::EDetPidStatus status = pidResp->CheckPIDStatus(AliPIDResponse::kTPC,track);
    if (status == AliPIDResponse::kDetPidOk){
      for(Int_t jsp=0; jsp<9; jsp++){
	nSigmaTPC[jsp]=pidResp->NumberOfSigmasTPC(track,(AliPID::EParticleType)jsp);
      }
    }
    fTreeVarFloat[15]=dedx;
    fTreeVarFloat[16]=nSigmaTPC[0];
    fTreeVarFloat[17]=nSigmaTPC[2];
    fTreeVarFloat[18]=nSigmaTPC[3];
    fTreeVarFloat[19]=nSigmaTPC[4];
    fTreeVarInt[6]=filtmap;
    
    Float_t ptgen=-999.;
    Float_t pgen=-999.;
    Float_t pxgen=-999.;
    Float_t pygen=-999.;
    Float_t pzgen=-999.;
    Float_t etagen=-999.;
    Float_t phigen=-999.;
    Int_t hadronSpecie=-1;
    Float_t invptgen=-999.;
    Bool_t isPhysPrim=kFALSE;
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
      isPhysPrim=part->IsPhysicalPrimary();
      fTreeVarFloat[20]=pxgen;
      fTreeVarFloat[21]=pygen;
      fTreeVarFloat[22]=pzgen;
      fTreeVarFloat[23]=ptgen;
      fTreeVarFloat[24]=pgen;
      fTreeVarFloat[25]=etagen;
      fTreeVarFloat[26]=phigen;
      if(fUseMCId){
	Int_t pdg=TMath::Abs(part->GetPdgCode());
	if(pdg==11) hadronSpecie=0;
	else if(pdg==13) hadronSpecie=1;
	else if(pdg==211) hadronSpecie=2;
	else if(pdg==321) hadronSpecie=3;
	else if(pdg==2212) hadronSpecie=4;
      }
      fTreeVarInt[7]=part->Charge();
      fTreeVarInt[8]=trlabel;
      fTreeVarInt[9]=part->GetPdgCode();
      fTreeVarInt[10]=isPhysPrim;
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
	fHistTPCChi2ndfPtFiltBit[jb]->Fill(pttrack,chi2clus);
	fHistChi2TPCConstrVsGlobPtFiltBit[jb]->Fill(pttrack,goldenChi2);
      }
    }

    if(track->GetID()<0) continue;
    // convert to ESD track here
    AliESDtrack esdTrack(track);
    // set the TPC cluster info
    esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
    esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
    esdTrack.SetTPCPointsF(track->GetTPCNclsF());
    // needed to calculate the impact parameters
    esdTrack.RelateToVertex(&vESD,0.,3.);
    if(!fTrCutsTPC->AcceptTrack(&esdTrack)) continue;
    if(track->GetTPCsignalN()<fMinNumOfTPCPIDclu) continue;

    fHistEtaPhiPtTPCsel->Fill(etatrack,phitrack,pttrack);
    if(tofBC==0) fHistEtaPhiPtTPCselTOFbc->Fill(etatrack,phitrack,pttrack);
    if(itsRefit){
      fHistEtaPhiPtTPCselITSref->Fill(etatrack,phitrack,pttrack);
      if(tofBC==0) fHistEtaPhiPtTPCselITSrefTOFbc->Fill(etatrack,phitrack,pttrack);
      if(spdAny){ 
	fHistEtaPhiPtTPCselSPDany->Fill(etatrack,phitrack,pttrack);
	if(tofBC==0) fHistEtaPhiPtTPCselSPDanyTOFbc->Fill(etatrack,phitrack,pttrack);
      }
    }

    fHistTPCchi2PerClusPhiPtTPCsel->Fill(chi2clus,pttrack,phitrack);
    if(itsRefit){
      fHistTPCchi2PerClusPhiPtTPCselITSref->Fill(chi2clus,pttrack,phitrack);
      if(spdAny) fHistTPCchi2PerClusPhiPtTPCselSPDany->Fill(chi2clus,pttrack,phitrack);
    }

    Bool_t isProton=kFALSE;
    Bool_t isKaon=kFALSE;
    Bool_t isPion=kFALSE;
    if(fReadMC && fUseMCId){
      if(hadronSpecie==4) isProton=kTRUE;
      if(hadronSpecie==3) isKaon=kTRUE;
      if(hadronSpecie==2) isPion=kTRUE;
    }else{
      if(TMath::Abs(nSigmaTPC[AliPID::kProton])<3) isProton=kTRUE;
      if(TMath::Abs(nSigmaTPC[AliPID::kKaon])<3) isKaon=kTRUE;
      if(TMath::Abs(nSigmaTPC[AliPID::kPion])<3) isPion=kTRUE;
    }

    if(itsRefit && spdAny){
      if(isPion) fHistImpParXYPtMulPionTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1); 
      if(isKaon) fHistImpParXYPtMulKaonTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1); 
      if(isProton) fHistImpParXYPtMulProtonTPCselSPDany->Fill(pttrack,impactXY*10000.,ncl1);
    }

    if(fReadMC){
      fHistPtResidVsPtTPCselAll->Fill(ptgen,(pttrack-ptgen));
      if(isPion) fHistPtResidVsPtTPCselPion->Fill(ptgen,(pttrack-ptgen));
      if(isKaon) fHistPtResidVsPtTPCselKaon->Fill(ptgen,(pttrack-ptgen));
      if(isProton) fHistPtResidVsPtTPCselProton->Fill(ptgen,(pttrack-ptgen));
      if(itsRefit){
	fHistPtResidVsPtTPCselITSrefAll->Fill(ptgen,(pttrack-ptgen));
	if(isPion) fHistPtResidVsPtTPCselITSrefPion->Fill(ptgen,(pttrack-ptgen));
	if(isKaon) fHistPtResidVsPtTPCselITSrefKaon->Fill(ptgen,(pttrack-ptgen));
	if(isProton) fHistPtResidVsPtTPCselITSrefProton->Fill(ptgen,(pttrack-ptgen));
      }
    }
    if(trlabel>=0){
      fHistEtaPhiPtTPCselITSrefGood->Fill(etatrack,phitrack,pttrack);
      if(itsRefit && spdAny) fHistImpParXYPtMulTPCselSPDanyGood->Fill(pttrack,impactXY*10000.,ncl1);
    }else{
      fHistEtaPhiPtTPCselITSrefFake->Fill(etatrack,phitrack,pttrack);
      if(itsRefit && spdAny) fHistImpParXYPtMulTPCselSPDanyFake->Fill(pttrack,impactXY*10000.,ncl1);
    }

    
    if(itsRefit && spdAny){
      if(isPhysPrim) fHistImpParXYPtMulTPCselSPDanyPrim->Fill(pttrack,impactXY*10000.,ncl1);
      else fHistImpParXYPtMulTPCselSPDanySec->Fill(pttrack,impactXY*10000.,ncl1);
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
    Double_t invMassK0s = v0->MassK0Short();
    Double_t invMassLambda = v0->MassLambda();
    Double_t invMassAntiLambda = v0->MassAntiLambda();
    Double_t ptv0=v0->Pt();

    AliESDtrack pEsdTrack(pTrack);
    pEsdTrack.SetTPCClusterMap(pTrack->GetTPCClusterMap());
    pEsdTrack.SetTPCSharedMap(pTrack->GetTPCSharedMap());
    pEsdTrack.SetTPCPointsF(pTrack->GetTPCNclsF());
    if(!fTrCutsTPC->AcceptTrack(&pEsdTrack)) continue;
    AliESDtrack nEsdTrack(pTrack);
    nEsdTrack.SetTPCClusterMap(nTrack->GetTPCClusterMap());
    nEsdTrack.SetTPCSharedMap(nTrack->GetTPCSharedMap());
    nEsdTrack.SetTPCPointsF(nTrack->GetTPCNclsF());
    if(!fTrCutsTPC->AcceptTrack(&nEsdTrack)) continue;

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

    if(keepK0s) fHistInvMassK0s->Fill(invMassK0s,ptv0);
    if(keepLambda){
      fHistInvMassLambda->Fill(invMassLambda,ptv0,pTrack->Pt());
    }
    if(keepAntiLambda){
      fHistInvMassAntiLambda->Fill(invMassAntiLambda,ptv0,nTrack->Pt());
    }
  }
  PostData(1,fOutput);
  
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
  printf("AliAnalysisTaskCheckAODTracks::Terminate --- Number of events: read = %.0f  analysed = %.0f\n",fHistNEvents->GetBinContent(1),fHistNEvents->GetBinContent(4));
  return;
}





