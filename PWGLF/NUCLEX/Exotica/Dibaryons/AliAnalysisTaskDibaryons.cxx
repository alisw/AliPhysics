#include <THashList.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TVector3.h>
#include <TChain.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliPIDResponse.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAODHeader.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDv0.h"
#include "AliAODv0.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliAnalysisTaskDibaryons.h"

using namespace std;

ClassImp(AliAnalysisTaskDibaryons)
//_______________________________________________________________________________________________
AliAnalysisTaskDibaryons::AliAnalysisTaskDibaryons():
  AliAnalysisTaskSE(),
  fAliEventCuts(),
  fAnalysisType("ESD"),
  fCollidingSystem(0),
  fkTriggerClass(AliVEvent::kINT7),
  fPIDResponse(0),
  fPileupCut(kTRUE),
  fOutput(0)
{
  // default constructor
}
//_______________________________________________________________________________________________
AliAnalysisTaskDibaryons::AliAnalysisTaskDibaryons(const char *name):
  AliAnalysisTaskSE(name),
  fAliEventCuts(),
  fAnalysisType("ESD"),
  fCollidingSystem(0),
  fkTriggerClass(AliVEvent::kINT7),
  fPIDResponse(0),
  fPileupCut(kTRUE),
  fOutput(0)
{
  // constructor

  DefineInput(0,TChain::Class());
  DefineOutput(1,THashList::Class());
}
//_______________________________________________________________________________________________
AliAnalysisTaskDibaryons::~AliAnalysisTaskDibaryons()
{
  // destructor

  delete fOutput;
}
//_______________________________________________________________________________________________
void AliAnalysisTaskDibaryons::UserCreateOutputObjects()
{
  // create output objects
  // called once

  fOutput = new THashList();
  fOutput->SetOwner(kTRUE);

  fPIDResponse = dynamic_cast<AliPIDResponse*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetPIDResponse());

  if(fkTriggerClass != AliVEvent::kINT7) {
    fAliEventCuts.OverrideAutomaticTriggerSelection(fkTriggerClass,true);
  }
  fAliEventCuts.AddQAplotsToList(fOutput);

  TH1F *hNPartStatistics = new TH1F("hNPartStatistics","Number of candidates under certain condition",10,0.5,10.5);
  hNPartStatistics->GetXaxis()->SetBinLabel(1,"p");
  hNPartStatistics->GetXaxis()->SetBinLabel(2,"#bar{p}");
  hNPartStatistics->GetXaxis()->SetBinLabel(3,"#Lambda");
  hNPartStatistics->GetXaxis()->SetBinLabel(4,"#bar{#Lambda}");
  hNPartStatistics->GetXaxis()->SetBinLabel(5,"#Xi^{-}");
  hNPartStatistics->GetXaxis()->SetBinLabel(6,"#Xi^{+}");
  hNPartStatistics->GetXaxis()->SetBinLabel(7,"#Omega^{-}");
  hNPartStatistics->GetXaxis()->SetBinLabel(8,"#Omega^{+}");

  fOutput->Add(hNPartStatistics);

  // Define histograms related to invariant mass for V0 and Cascade
  TH1F *hInvMassLambdawoCuts = new TH1F("hInvMassLambdawoCuts","Invariant mass of Lambda(p #pi-) without topological cuts;M_{p#pi^{-}} (GeV/c^{2})",400,1.0,1.2);
  TH1F *hInvMassLambdawCuts = new TH1F("hInvMassLambdawCuts","Invariant mass of Lambda(p #pi-) with topological cuts;M_{p#pi^{-}} (GeV/c^{2})",400,1.0,1.2);
  TH1F *hInvMassAntiLambdawoCuts = new TH1F("hInvMassAntiLambdawoCuts","Invariant mass of AntiLambda(a-p #pi+) without topological cuts;M_{#bar{p}#pi^{+}} (GeV/c^{2})",400,1.0,1.2);
  TH1F *hInvMassAntiLambdawCuts = new TH1F("hInvMassAntiLambdawCuts","Invariant mass of AntiLambda(a-p #pi+) with topological cuts;M_{#bar{p}#pi^{+}} (GeV/c^{2})",400,1.0,1.2);
  TH1F *hInvMassLambdaAsCascDghter = new TH1F("hInvMassLambdaAsCascDghter","Invariant mass of Lambda(p #pi-) stemming from Xi;M_{p#pi^{-}} (GeV/c^{2})",400,1.0,1.2);
  TH1F *hInvMassXimwoCuts = new TH1F("hInvMassXimwoCuts","Invariant mass of Xi-(p #pi- #pi-) without topological cuts;M_{p#pi^{-}#pi^{-}} (GeV/c^{2})",300,1.2,1.5);
  TH1F *hInvMassXimwCuts = new TH1F("hInvMassXimwCuts","Invariant mass of Xi-(p #pi- #pi-) with topological cuts;M_{p#pi^{-}#pi^{-}} (GeV/c^{2})",300,1.2,1.5);
  TH1F *hInvMassXipwoCuts = new TH1F("hInvMassXipwoCuts","Invariant mass of Xi+(a-p #pi+ #pi+) without topological cuts;M_{#bar{p}#pi^{+}#pi^{+}} (GeV/c^{2})",300,1.2,1.5);
  TH1F *hInvMassXipwCuts = new TH1F("hInvMassXipwCuts","Invariant mass of Xi+(a-p #pi+ #pi+) with topological cuts;M_{#bar{p}#pi^{+}#pi^{+}} (GeV/c^{2})",300,1.2,1.5);
  TH1F *hInvMassOmegamwoCuts = new TH1F("hInvMassOmegamwoCuts","Invariant mass of Omega-(p K- #pi-) without topological cuts;M_{pK^{-}#pi^{-}} (GeV/c^{2})",300,1.5,1.8);
  TH1F *hInvMassOmegamwCuts = new TH1F("hInvMassOmegamwCuts","Invariant mass of Omega-(p K- #pi-) with topological cuts;M_{pK^{-}#pi^{-}} (GeV/c^{2})",300,1.5,1.8);
  TH1F *hInvMassOmegapwoCuts = new TH1F("hInvMassOmegapwoCuts","Invariant mass of Omega+(a-p K+ #pi+) without topological cuts;M_{#bar{p}K^{+}#pi^{+}} (GeV/c^{2})",300,1.5,1.8);
  TH1F *hInvMassOmegapwCuts = new TH1F("hInvMassOmegapwCuts","Invariant mass of Omega+(a-p K+ #pi+) with topological cuts;M_{#bar{p}K^{+}#pi^{+}} (GeV/c^{2})",300,1.5,1.8);

  fOutput->Add(hInvMassLambdawoCuts);
  fOutput->Add(hInvMassLambdawCuts);
  fOutput->Add(hInvMassAntiLambdawoCuts);
  fOutput->Add(hInvMassAntiLambdawCuts);
  fOutput->Add(hInvMassLambdaAsCascDghter);
  fOutput->Add(hInvMassXimwoCuts);
  fOutput->Add(hInvMassXimwCuts);
  fOutput->Add(hInvMassXipwoCuts);
  fOutput->Add(hInvMassXipwCuts);
  fOutput->Add(hInvMassOmegamwoCuts);
  fOutput->Add(hInvMassOmegamwCuts);
  fOutput->Add(hInvMassOmegapwoCuts);
  fOutput->Add(hInvMassOmegapwCuts);

  // Define QA plots for topological observables
  TH2F *hProtonDCAxyDCAz = new TH2F("hProtonDCAxyDCAz","Distribution of DCAz vs DCAxy;DCA_{xy} (cm);DCA_{z} (cm)",500,-5,5,1000,-20,20);
  TH1F *hLambdaDCADaughterTracks = new TH1F("hLambdaDCADaughterTracks","Distance of closest approach of p and #pi track;DCA (cm)",100,0,10);
  TH1F *hLambdaDCAPosDaughPrimVertex = new TH1F("hLambdaDCAPosDaughPrimVertex","Distance of closest approach of V0 positive daughter track to primary vertex;DCA (cm)",500,0,100);
  TH1F *hLambdaDCANegDaughPrimVertex = new TH1F("hLambdaDCANegDaughPrimVertex","Distance of closest approach of V0 negative daughter track to primary vertex;DCA (cm)",500,0,100);
  TH1F *hLambdaTransverseRadius = new TH1F("hLambdaTransverseRadius","Transverse distance between primary vertex and V0 decay vertex;r_{xy} (cm);",400,0,200);
  TH1F *hLambdaCosPointingAngle = new TH1F("hLambdaCosPointingAngle","Cosine of pointing angle;cos(#Theta)",200,0.8,1);
  TH1F *hXiDCADaughterTracks = new TH1F("hXiDCADaughterTracks","Distance of closest approach of V0 and #pi track;DCA (cm)",400,0,2);
  TH1F *hXiDCAV0DaughterTracks = new TH1F("hXiDCAV0DaughterTracks","Distance of closest approach of p and #pi track;DCA (cm)",400,0,2);
  TH1F *hXiDCAV0PrimVertex = new TH1F("hXiDCAV0PrimVertex","Distance of closest approach of Xi V0 track to primary vertex;DCA (cm)",500,0,100);
  TH1F *hXiDCABachPrimVertex = new TH1F("hXiDCABachPrimVertex","Distance of closest approach of Xi bachelor track to primary vertex;DCA (cm)",500,0,100);
  TH1F *hXiDCAPosDaughPrimVertex = new TH1F("hXiDCAPosDaughPrimVertex","Distance of closest approach of Xi positive daughter track to primary vertex;DCA (cm)",500,0,100);
  TH1F *hXiDCANegDaughPrimVertex = new TH1F("hXiDCANegDaughPrimVertex","Distance of closest approach of Xi negative daughter track to primary vertex;DCA (cm)",500,0,100);
  TH1F *hXiTransverseRadius = new TH1F("hXiTransverseRadius","Transverse distance between primary vertex and Xi decay vertex;r_{xy} (cm);",400,0,200);
  TH1F *hXiV0TransverseRadius = new TH1F("hXiV0TransverseRadius","Transverse distance between primary vertex and V0 decay vertex;r_{xy} (cm);",400,0,200);
  TH1F *hXiCosPointingAngle = new TH1F("hXiCosPointingAngle","Cosine of pointing angle of Xi;cos(#Theta)",200,0.8,1);
  TH1F *hXiV0CosPointingAngle = new TH1F("hXiV0CosPointingAngle","Cosine of pointing angle of V0;cos(#Theta)",200,0.8,1);
  TH1F *hOmegaDCADaughterTracks = new TH1F("hOmegaDCADaughterTracks","Distance of closest approach of V0 and K track;DCA (cm)",100,0,10);
  TH1F *hOmegaDCAV0DaughterTracks = new TH1F("hOmegaDCAV0DaughterTracks","Distance of closest approach of p and #pi track;DCA (cm)",100,0,10);
  TH1F *hOmegaDCAV0PrimVertex = new TH1F("hOmegaDCAV0PrimVertex","Distance of closest approach of Omega V0 track to primary vertex;DCA (cm)",500,0,10);
  TH1F *hOmegaDCABachPrimVertex = new TH1F("hOmegaDCABachPrimVertex","Distance of closest approach of Omega bachelor track to primary vertex;DCA (cm)",500,0,10);
  TH1F *hOmegaDCAPosDaughPrimVertex = new TH1F("hOmegaDCAPosDaughPrimVertex","Distance of closest approach of Omega positive daughter track to primary vertex;DCA (cm)",500,0,10);
  TH1F *hOmegaDCANegDaughPrimVertex = new TH1F("hOmegaDCANegDaughPrimVertex","Distance of closest approach of Omega negative daughter track to primary vertex;DCA (cm)",500,0,10);
  TH1F *hOmegaTransverseRadius = new TH1F("hOmegaTransverseRadius","Transverse distance between primary vertex and Omega decay vertex;r_{xy} (cm);",400,0,200);
  TH1F *hOmegaV0TransverseRadius = new TH1F("hOmegaV0TransverseRadius","Transverse distance between primary vertex and V0 decay vertex;r_{xy} (cm);",400,0,200);
  TH1F *hOmegaCosPointingAngle = new TH1F("hOmegaCosPointingAngle","Cosine of pointing angle;cos(#Theta)",200,0.8,1);
  TH1F *hOmegaV0CosPointingAngle = new TH1F("hOmegaV0CosPointingAngle","Cosine of pointing angle;cos(#Theta)",200,0.8,1);

  fOutput->Add(hProtonDCAxyDCAz);
  fOutput->Add(hLambdaDCADaughterTracks);
  fOutput->Add(hLambdaDCAPosDaughPrimVertex);
  fOutput->Add(hLambdaDCANegDaughPrimVertex);
  fOutput->Add(hLambdaTransverseRadius);
  fOutput->Add(hLambdaCosPointingAngle);
  fOutput->Add(hXiDCADaughterTracks);
  fOutput->Add(hXiDCAV0DaughterTracks);
  fOutput->Add(hXiDCAV0PrimVertex);
  fOutput->Add(hXiDCABachPrimVertex);
  fOutput->Add(hXiDCAPosDaughPrimVertex);
  fOutput->Add(hXiDCANegDaughPrimVertex);
  fOutput->Add(hXiTransverseRadius);
  fOutput->Add(hXiV0TransverseRadius);
  fOutput->Add(hXiCosPointingAngle);
  fOutput->Add(hXiV0CosPointingAngle);
  fOutput->Add(hOmegaDCADaughterTracks);
  fOutput->Add(hOmegaDCAV0DaughterTracks);
  fOutput->Add(hOmegaDCAV0PrimVertex);
  fOutput->Add(hOmegaDCABachPrimVertex);
  fOutput->Add(hOmegaDCAPosDaughPrimVertex);
  fOutput->Add(hOmegaDCANegDaughPrimVertex);
  fOutput->Add(hOmegaTransverseRadius);
  fOutput->Add(hOmegaV0TransverseRadius);
  fOutput->Add(hOmegaCosPointingAngle);
  fOutput->Add(hOmegaV0CosPointingAngle);

  // Define QA plots for kinematic observables
  TH1F *hProtonPt = new TH1F("hProtonPt","Transverse momentum of proton;p_{T} (GeV/c)",500,0,10);
  TH1F *hProtonPhi = new TH1F("hProtonPhi","Phi angle of proton;#varphi (rad)",100,0,TMath::TwoPi());
  TH1F *hProtonEta = new TH1F("hProtonEta","Pseudorapidity of proton;#eta",200,-1,+1);
  TH1F *hAntiProtonPt = new TH1F("hAntiProtonPt","Transverse momentum of Anti-proton;p_{T} (GeV/c)",1000,0,10);
  TH1F *hLambdaPt = new TH1F("hLambdaPt","Transverse momentum of V0;p_{T} (GeV/c)",500,0,10);
  TH1F *hLambdaPosDaughPt = new TH1F("hLambdaPosDaughPt","Transverse momentum of positive daughter of V0;p_{T} (GeV/c)",500,0,10);
  TH1F *hLambdaNegDaughPt = new TH1F("hLambdaNegDaughPt","Transverse momentum of negative daughter of V0;p_{T} (GeV/c)",500,0,10);
  TH1F *hLambdaPhi = new TH1F("hLambdaPhi","Phi angle of V0;#varphi (rad)",100,0,TMath::TwoPi());
  TH1F *hLambdaEta = new TH1F("hLambdaEta","Pseudorapidity of V0;#eta",400,-2,2);
  TH1F *hAntiLambdaPt = new TH1F("hAntiLambdaPt","Transverse momentum of Anti-V0;p_{T} (GeV/c)",500,0,10);
  TH1F *hXimPt = new TH1F("hXimPt","Transverse momentum of Xi-;p_{T} (GeV/c)",500,0,10);
  TH1F *hXimBachPt = new TH1F("hXimBachPt","Transverse momentum of bachelor of Xi-;p_{T} (GeV/c)",500,0,10);
  TH1F *hXimPosDaughPt = new TH1F("hXimPosDaughPt","Transverse momentum of positive daughter of Xi-;p_{T} (GeV/c)",500,0,10);
  TH1F *hXimNegDaughPt = new TH1F("hXimNegDaughPt","Transverse momentum of negative daughter of Xi-;p_{T} (GeV/c)",500,0,10);
  TH1F *hXimPhi = new TH1F("hXimPhi","Phi angle of Xi-;#varphi (rad)",100,0,TMath::TwoPi());
  TH1F *hXimEta = new TH1F("hXimEta","Pseudorapidity of Xi-;#eta",400,-2,2);
  TH1F *hXipPt = new TH1F("hXipPt","Transverse momentum of Xi+;p_{T} (GeV/c)",500,0,10);
  TH1F *hOmegamPt = new TH1F("hOmegamPt","Transverse momentum of Omega-;p_{T} (GeV/c)",500,0,10);
  TH1F *hOmegamBachPt = new TH1F("hOmegamBachPt","Transverse momentum of bachelor of Omega-;p_{T} (GeV/c)",500,0,10);
  TH1F *hOmegamPosDaughPt = new TH1F("hOmegamPosDaughPt","Transverse momentum of positive daughter of Omega-;p_{T} (GeV/c)",500,0,10);
  TH1F *hOmegamNegDaughPt = new TH1F("hOmegamNegDaughPt","Transverse momentum of negative daughter of Omega-;p_{T} (GeV/c)",500,0,10);
  TH1F *hOmegamPhi = new TH1F("hOmegamPhi","Phi angle of Omega-;#varphi (rad)",100,0,TMath::TwoPi());
  TH1F *hOmegamEta = new TH1F("hOmegamEta","Pseudorapidity of Omega-;#eta",400,-2,2);
  TH1F *hOmegapPt = new TH1F("hOmegapPt","Transverse momentum of Omega+;p_{T} (GeV/c)",500,0,10);

  fOutput->Add(hProtonPt);
  fOutput->Add(hProtonPhi);
  fOutput->Add(hProtonEta);
  fOutput->Add(hAntiProtonPt);
  fOutput->Add(hLambdaPt);
  fOutput->Add(hLambdaPosDaughPt);
  fOutput->Add(hLambdaNegDaughPt);
  fOutput->Add(hLambdaPhi);
  fOutput->Add(hLambdaEta);
  fOutput->Add(hAntiLambdaPt);
  fOutput->Add(hXimPt);
  fOutput->Add(hXimBachPt);
  fOutput->Add(hXimPosDaughPt);
  fOutput->Add(hXimNegDaughPt);
  fOutput->Add(hXimPhi);
  fOutput->Add(hXimEta);
  fOutput->Add(hXipPt);
  fOutput->Add(hOmegamPt);
  fOutput->Add(hOmegamBachPt);
  fOutput->Add(hOmegamPosDaughPt);
  fOutput->Add(hOmegamNegDaughPt);
  fOutput->Add(hOmegamPhi);
  fOutput->Add(hOmegamEta);
  fOutput->Add(hOmegapPt);

  // Define histogrms related to TPC track info
  TH1F *hNCrossedRows = new TH1F("hNCrossedRows","Number of Crossed Rows",160,0,160);
  TH1F *hNCluster = new TH1F("hNCluster","Number of TPC Clusters",160,0,160);
  TH1F *hNSharedCluster = new TH1F("hNSharedCluster","Number of Shared Clusters",160,0,160);
  TH1F *hNFindableCluster = new TH1F("hNFindableCluster","Number of Findable Clusters",160,0,160);
  TH1F *hRatioFindableCrossed = new TH1F("hRatioFindableCrossed","Ratio of Number of Crossed Rows / Number of Findable Clusters",150,0,1.5);

  fOutput->Add(hNCrossedRows);
  fOutput->Add(hNCluster);
  fOutput->Add(hNSharedCluster);
  fOutput->Add(hNFindableCluster);
  fOutput->Add(hRatioFindableCrossed);

  // Define PID related histograms
  TH2F *hdEdxVsP = new TH2F("hdEdxVsP","dE/dx of all particles vs momentum;p (GeV/c);#frac{dE}{dx} (a.u.)",1000,0,10,200,0,200);
  TH2F *hProtonNsigmaTPC = new TH2F("hProtonNsigmaTPC","nsigma_TPC of protons;p (GeV/c);n#sigma_{TPC}",1000,0,10,200,-10,10);
  TH2F *hProtonNsigmaTOF = new TH2F("hProtonNsigmaTOF","nsigma_TOF of protons;p (GeV/c);n#sigma_{TOF}",1000,0,10,200,-10,10);
  TH2F *hProtonNsigmaCombined = new TH2F("hProtonNsigmaCombined","nsigma_combined of protons;p (GeV/c);n#sigma_{comb}=#sqrt{n#sigma_{TPC}^{2}+n#sigma_{TOF}^{2}};",1000,0,10,100,0,10);
  TH2F *hProtonNsigmaTPCwPID = new TH2F("hProtonNsigmaTPCwPID","nsigma_TPC of protons with PID;p (GeV/c);n#sigma_{TPC}",1000,0,10,200,-10,10);
  TH2F *hProtonNsigmaCombinedwPID = new TH2F("hProtonNsigmaCombinedwPID","nsigma_combined of protons with PID;p (GeV/c);n#sigma_{comb}=#sqrt{n#sigma_{TPC}^{2}+n#sigma_{TOF}^{2}};",1000,0,10,100,0,10);

  fOutput->Add(hdEdxVsP);
  fOutput->Add(hProtonNsigmaTPC);
  fOutput->Add(hProtonNsigmaTOF);
  fOutput->Add(hProtonNsigmaCombined);
  fOutput->Add(hProtonNsigmaTPCwPID);
  fOutput->Add(hProtonNsigmaCombinedwPID);

  // Define histograms related to pair analysis
  TH1F *hNLambdaLambdaPairs = new TH1F("hNLambdaLambdaPairs","Number of V0-V0 pairs in an event",100,0,100);
  TH1F *hNProtonXiPairs = new TH1F("hNProtonXiPairs","Number of Proton-Xi pairs in an event",100,0,100);
  TH1F *hNLambdaXiPairs = new TH1F("hNLambdaXiPairs","Number of V0-Xi pairs in an event",100,0,100);
  TH1F *hNPairStatistics = new TH1F("hNPairStatistics","Number of pairs under certain condition",10,0.5,10.5);
  hNPairStatistics->GetXaxis()->SetBinLabel(1,"p-#Lambda");
  hNPairStatistics->GetXaxis()->SetBinLabel(2,"#Lambda-#Lambda");
  hNPairStatistics->GetXaxis()->SetBinLabel(3,"p-#Xi^{-}");
  hNPairStatistics->GetXaxis()->SetBinLabel(4,"#Lambda-#Xi^{-}");
  hNPairStatistics->GetXaxis()->SetBinLabel(5,"#Xi^{-}-#Omega^{-}");

  fOutput->Add(hNLambdaLambdaPairs);
  fOutput->Add(hNProtonXiPairs);
  fOutput->Add(hNLambdaXiPairs);
  fOutput->Add(hNPairStatistics);

  // Define histograms related to invariant mass for dibaryons
  TH2F *hInvMassRelPLambdaLambda = new TH2F("hInvMassRelPLambdaLambda","Invariant mass vs relative momentum of Lambda-Lambda pair;M_{#Lambda#Lambda} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPLambdaLambdaRsideband = new TH2F("hInvMassRelPLambdaLambdaRsideband","Invariant mass vs relative momentum of Lambda and Lambda in right sideband;M_{#Lambda#Lambda} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPLambdaLambdaLsideband = new TH2F("hInvMassRelPLambdaLambdaLsideband","Invariant mass vs relative momentum of Lambda and Lambda in left sideband;M_{#Lambda#Lambda} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPProtonLambda = new TH2F("hInvMassRelPProtonLambda","Invariant mass vs relative momentum of proton-Lambda pair;M_{p#Lambda} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPProtonLambdaRsideband = new TH2F("hInvMassRelPProtonLambdaRsideband","Invariant mass vs relative momentum of proton and Lambda in right side band;M_{p#Lambda} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPProtonLambdaLsideband = new TH2F("hInvMassRelPProtonLambdaLsideband","Invariant mass vs relative momentum of proton and Lambda in left side band;M_{p#Lambda} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPProtonXi = new TH2F("hInvMassRelPProtonXi","Invariant mass vs relative momentum of Proton-Xi pair;M_{p#Xi} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPProtonXiRsideband = new TH2F("hInvMassRelPProtonXiRsideband","Invariant mass vs relative momentum of proton and Xi in right sideband;M_{p#Xi} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPProtonXiLsideband = new TH2F("hInvMassRelPProtonXiLsideband","Invariant mass vs relative momentum of proton and Xi in left sideband;M_{p#Xi} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPLambdaXi = new TH2F("hInvMassRelPLambdaXi","Invariant mass vs relative momentum of Lambda-Xi pair;M_{#Lambda#Xi} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPLambdaXiRsideband = new TH2F("hInvMassRelPLambdaXiRsideband","Invariant mass vs relative momentum of Lambda and Xi in right sideband;M_{#Lambda#Xi} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPLambdaXiLsideband = new TH2F("hInvMassRelPLambdaXiLsideband","Invariant mass vs relative momentum of Lambda and Xi in left sideband;M_{#Lambda#Xi} (GeV/c^{2});p_{Rel} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassRelPXiOmega = new TH2F("hInvMassRelPXiOmega","Invariant mass vs relative momentum of Xi-Omega pair;M_{#Xi#Omega} (GeV/c^{2});p_{Rel} (GeV/c)",1000,3,4,100,0,10);
  TH2F *hInvMassRelPXiOmegaRsideband = new TH2F("hInvMassRelPXiOmegaRsideband","Invariant mass vs relative momentum of Xi and Omega in right sideband;M_{#Xi#Omega} (GeV/c^{2});p_{Rel} (GeV/c)",1000,3,4,100,0,10);
  TH2F *hInvMassRelPXiOmegaLsideband = new TH2F("hInvMassRelPXiOmegaLsideband","Invariant mass vs relative momentum of Xi and Omega in left sideband;M_{#Xi#Omega} (GeV/c^{2});p_{Rel} (GeV/c)",1000,3,4,100,0,10);

  fOutput->Add(hInvMassRelPLambdaLambda);
  fOutput->Add(hInvMassRelPLambdaLambdaRsideband);
  fOutput->Add(hInvMassRelPLambdaLambdaLsideband);
  fOutput->Add(hInvMassRelPProtonLambda);
  fOutput->Add(hInvMassRelPProtonLambdaRsideband);
  fOutput->Add(hInvMassRelPProtonLambdaLsideband);
  fOutput->Add(hInvMassRelPProtonXi);
  fOutput->Add(hInvMassRelPProtonXiRsideband);
  fOutput->Add(hInvMassRelPProtonXiLsideband);
  fOutput->Add(hInvMassRelPLambdaXi);
  fOutput->Add(hInvMassRelPLambdaXiRsideband);
  fOutput->Add(hInvMassRelPLambdaXiLsideband);
  fOutput->Add(hInvMassRelPXiOmega);
  fOutput->Add(hInvMassRelPXiOmegaRsideband);
  fOutput->Add(hInvMassRelPXiOmegaLsideband);

  PostData(1,fOutput);
}

//_______________________________________________________________________________________________
void AliAnalysisTaskDibaryons::UserExec(Option_t *option)
{
  // Main loop
  // Called for each event

  AliESDEvent *esdEvent = 0x0;
  AliAODEvent *aodEvent = 0x0;

  // Check the PID response
  if(!fPIDResponse) {
    AliError("Cannot get pid response");
    return;
  }

  // Load the InputEvent and check it
  if(fAnalysisType == "ESD") {
    esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
    if(!esdEvent) {
      AliWarning("ERROR: esdEvent not available");
      return;
    }
  } else if(fAnalysisType == "AOD") {
    aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!aodEvent) {
      AliWarning("ERROR: aodEvent not available");
      return;
    }
  } else {
    AliFatal("Analysis type (ESD or AOD) not specified");
    return;
  }

  // Event selection
  if(!fAliEventCuts.AcceptEvent(InputEvent())) {
    AliInfo("This event is rejected by AliEventCuts::AcceptEvent() ... return!");
    return;
  }

  // Primary Vertex position
  Double_t primaryVtxPos[3] = {-999., -999., -999.};
  if(fAnalysisType == "ESD") {
    const AliESDVertex *esdVtx = esdEvent->GetPrimaryVertex();
    if(!esdVtx) {
      AliWarning("No prim. vertex in ESD... return!");
      return;
    }
    esdVtx->GetXYZ(primaryVtxPos);
  } else if(fAnalysisType == "AOD") {
    const AliAODVertex *aodVtx = aodEvent->GetPrimaryVertex();
    if(!aodVtx) {
      AliWarning("No prim. vertex in AOD... return!");
      return;
    }
    aodVtx->GetXYZ(primaryVtxPos);
  }

  // Magnetic field
  Double_t bz = -10.;
  if     (fAnalysisType == "ESD") bz = esdEvent->GetMagneticField();
  else if(fAnalysisType == "AOD") bz = aodEvent->GetMagneticField();

  // PDG mass
  const Double_t massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass(); 
  const Double_t massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass(); 
  const Double_t massXi     = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  const Double_t massOmega  = TDatabasePDG::Instance()->GetParticle(3334)->Mass();


  //______________________________________________________________________________
  // Loop over the reconstructed candidates

  // - Track loop

  Int_t nTrack = 0;
  if     (fAnalysisType == "ESD") nTrack = esdEvent->GetNumberOfTracks();
  else if(fAnalysisType == "AOD") nTrack = aodEvent->GetNumberOfTracks();

  TObjArray candProton(nTrack);
  TObjArray candAntiProton(nTrack);

  for(Int_t iTrack=0; iTrack < nTrack; iTrack++) {// This is the beginning of the track loop

    AliESDtrack *esdTrack = 0x0;
    AliAODTrack *aodTrack = 0x0;

    // Initialisation of the local variables that will be needed for ESD/AOD
    Double_t charge    = 0.;
    Double_t momX      = 0.;
    Double_t momY      = 0.;
    Double_t momZ      = 0.;
    Double_t transvMom = 0.;
    Double_t totMom    = 0.;
    Double_t eta       = 0.;
    Double_t phi       = 0.;
    Double_t dEdx      = 0.;

    Float_t nTPCCrossedRows  = -1.;
    UShort_t nTPCClusters    = -1;
    UShort_t nTPCSharedCls   = -1;
    UShort_t nTPCFindableCls = -1;
    Float_t ratio = -1.;

    Float_t nSigmaTPCproton      = -10.;
    Float_t nSigmaTOFproton      = -10.;
    Float_t nSigmaTPCTOFcombined = -10.;

    Float_t DCAxy = -999., DCAz = -999.;

    if(fAnalysisType == "ESD") {

      esdTrack = dynamic_cast<AliESDtrack*>(esdEvent->GetTrack(iTrack));
      if(!esdTrack) {
        AliFatal("Not a standard ESD");
        continue;
      }

      charge    = esdTrack->Charge();
      momX      = esdTrack->Px();
      momY      = esdTrack->Py();
      momZ      = esdTrack->Pz();
      transvMom = esdTrack->Pt();
      totMom    = esdTrack->P();
      eta       = esdTrack->Eta();
      phi       = esdTrack->Phi();
      dEdx      = esdTrack->GetTPCsignal();

      nTPCCrossedRows = esdTrack->GetTPCCrossedRows();
      nTPCClusters    = esdTrack->GetTPCNcls();
      nTPCSharedCls   = esdTrack->GetTPCnclsS();
      nTPCFindableCls = esdTrack->GetTPCNclsF();
      if(nTPCFindableCls > 0) ratio = nTPCCrossedRows / nTPCFindableCls;

      nSigmaTPCproton      = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
      nSigmaTOFproton      = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton);
      nSigmaTPCTOFcombined = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));

      esdTrack->GetImpactParameters(DCAxy,DCAz);

    } // end of ESD treatment

    else if(fAnalysisType == "AOD") {

      aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
      if(!aodTrack) {
        AliFatal("Not a standard AOD");
        continue;
      }

      charge    = aodTrack->Charge();
      momX      = aodTrack->Px();
      momY      = aodTrack->Py();
      momZ      = aodTrack->Pz();
      transvMom = aodTrack->Pt();
      totMom    = aodTrack->P();
      eta       = aodTrack->Eta();
      phi       = aodTrack->Phi();
      dEdx      = aodTrack->GetTPCsignal();

      nTPCCrossedRows = aodTrack->GetTPCCrossedRows();
      nTPCClusters    = aodTrack->GetTPCNcls();
      nTPCSharedCls   = aodTrack->GetTPCnclsS();
      nTPCFindableCls = aodTrack->GetTPCNclsF();
      if(nTPCFindableCls > 0) ratio = nTPCCrossedRows / nTPCFindableCls;

      nSigmaTPCproton      = fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton);
      nSigmaTOFproton      = fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kProton);
      nSigmaTPCTOFcombined = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));

      aodTrack->GetImpactParameters(DCAxy,DCAz);

    } // end of AOD treatment

    dynamic_cast<TH1F*>(fOutput->FindObject("hNCrossedRows"))        ->Fill(nTPCCrossedRows);
    dynamic_cast<TH1F*>(fOutput->FindObject("hNCluster"))            ->Fill(nTPCClusters);
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedCluster"))      ->Fill(nTPCSharedCls);
    dynamic_cast<TH1F*>(fOutput->FindObject("hNFindableCluster"))    ->Fill(nTPCFindableCls);
    dynamic_cast<TH1F*>(fOutput->FindObject("hRatioFindableCrossed"))->Fill(ratio);

    if(nTPCClusters    < 80)   continue;
    if(nTPCCrossedRows < 70)   continue;
    if(ratio           < 0.83) continue;
    if(nTPCSharedCls   > 0)    continue;

    dynamic_cast<TH2F*>(fOutput->FindObject("hdEdxVsP"))             ->Fill(totMom, dEdx);

    if(charge > 0.) {
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaTPC"))     ->Fill(totMom, nSigmaTPCproton);
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaTOF"))     ->Fill(totMom, nSigmaTOFproton);
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaCombined"))->Fill(totMom, nSigmaTPCTOFcombined);
    }

    // Proton PID
    if(     totMom < 0.75 && TMath::Abs(nSigmaTPCproton) > 3.) continue;
    else if(totMom > 0.75 && nSigmaTPCTOFcombined        > 3.) continue;

    if(charge > 0.) {
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaTPCwPID"))     ->Fill(totMom, nSigmaTPCproton);
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaCombinedwPID"))->Fill(totMom, nSigmaTPCTOFcombined);
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonDCAxyDCAz"))         ->Fill(DCAxy,DCAz);
    }

    if(TMath::Abs(eta) > 0.8) continue;
    if(transvMom < 0.5)  continue;
    if(transvMom > 4.05) continue;
    if(TMath::Abs(DCAz)  > 0.2) continue; 
    if(TMath::Abs(DCAxy) > 0.1) continue;

    // proton
    if(charge > 0.) {

      dynamic_cast<TH1F*>(fOutput->FindObject("hProtonPt")) ->Fill(transvMom);
      dynamic_cast<TH1F*>(fOutput->FindObject("hProtonPhi"))->Fill(phi);
      dynamic_cast<TH1F*>(fOutput->FindObject("hProtonEta"))->Fill(eta);
      dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(1);

      if(fAnalysisType == "ESD")      candProton.Add(esdTrack);
      else if(fAnalysisType == "AOD") candProton.Add(aodTrack);

    }
    // Anti-proton
    else if(charge < 0.) {

      dynamic_cast<TH1F*>(fOutput->FindObject("hAntiProtonPt"))->Fill(transvMom);
      dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(2);

      if(fAnalysisType == "ESD")      candAntiProton.Add(esdTrack);
      else if(fAnalysisType == "AOD") candAntiProton.Add(aodTrack);

    }

  }// end of track loop


  // - V0 loop

  Int_t nV0 = 0;
  if     (fAnalysisType == "ESD") nV0 = esdEvent->GetNumberOfV0s();
  else if(fAnalysisType == "AOD") nV0 = aodEvent->GetNumberOfV0s();

  TObjArray candLambda(nV0);
  TObjArray candAntiLambda(nV0);
  TObjArray rsideLambda(nV0);
  TObjArray rsideAntiLambda(nV0);
  TObjArray lsideLambda(nV0);
  TObjArray lsideAntiLambda(nV0);
  
  for(Int_t iV0=0; iV0<nV0; iV0++) {// This is the beginning of the V0 loop

    AliESDv0 *esdV0 = 0x0;
    AliAODv0 *aodV0 = 0x0;

    // Initialisation of the local variables that will be needed for ESD/AOD
    Double_t invMassLambda     = 0.;
    Double_t invMassAntiLambda = 0.;

    Double_t dcaV0Dghters    = -1.; 
    Double_t dcaV0ToPrimVtx  = -1.;
    Double_t cosPointAngle   = -1.;
    Double_t vtxPosV0[3]     = {-999., -999., -999.};
    Double_t radius          = -999.;
    Double_t dcaPosToPrimVtx = -1.;
    Double_t dcaNegToPrimVtx = -1.;

    Bool_t isPosProton = kFALSE;
    Bool_t isPosPion   = kFALSE;
    Bool_t isNegProton = kFALSE;
    Bool_t isNegPion   = kFALSE;

    Double_t etaV0        = -20.;
    Double_t phiV0        = 720.;
    Double_t momV0X       = 0.;
    Double_t momV0Y       = 0.;
    Double_t momV0Z       = 0.;
    Double_t transvMomV0  = 0.;
    Double_t etaPos       = -20.;
    Double_t momPosX      = 0.;
    Double_t momPosY      = 0.;
    Double_t momPosZ      = 0.;
    Double_t transvMomPos = 0.;
    Double_t etaNeg       = -20.;
    Double_t momNegX      = 0.;
    Double_t momNegY      = 0.;
    Double_t momNegZ      = 0.;
    Double_t transvMomNeg = 0.;

    if(fAnalysisType == "ESD") {
      
      esdV0 = esdEvent->GetV0(iV0);
      if(!esdV0) continue;

      if(esdV0->GetOnFlyStatus()) continue; // select offline v0

      // Get the tracks for the daughters
      AliESDtrack *pTrack = esdEvent->GetTrack(esdV0->GetPindex());
      AliESDtrack *nTrack = esdEvent->GetTrack(esdV0->GetNindex());

      if(!pTrack || !nTrack) {
        AliWarning("ERROR: Could not retrieve one of the 2 ESD daughter tracks of the V0 ...");
        continue;
      }

      // Rejection of a poor quality tracks
      if(pTrack->GetTPCNcls() < 70) continue;
      if(nTrack->GetTPCNcls() < 70) continue;
      if(TMath::Abs(pTrack->Eta()) > 0.8) continue;
      if(TMath::Abs(nTrack->Eta()) > 0.8) continue;
      if(esdV0->Pt() < 0.3) continue;

      // Daughter track PID using TPC
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton)) < 5.0) isPosProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion  )) < 5.0) isPosPion   = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton)) < 5.0) isNegProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion  )) < 5.0) isNegPion   = kTRUE;

      // Calculate the invariant mass
      esdV0->ChangeMassHypothesis(3122);
      invMassLambda     = esdV0->GetEffMass();
      esdV0->ChangeMassHypothesis(-3122);
      invMassAntiLambda = esdV0->GetEffMass();

      // Get topological values
      dcaV0Dghters    = esdV0->GetDcaV0Daughters();
      dcaPosToPrimVtx = TMath::Abs(pTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));
      dcaNegToPrimVtx = TMath::Abs(pTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));
      cosPointAngle   = esdV0->GetV0CosineOfPointingAngle(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]);
      esdV0->GetXYZ(vtxPosV0[0], vtxPosV0[1], vtxPosV0[2]);
      radius = TMath::Sqrt(vtxPosV0[0]*vtxPosV0[0] + vtxPosV0[1]*vtxPosV0[1]);

      esdV0->GetPxPyPz(momV0X, momV0Y, momV0Z);
      esdV0->GetPPxPyPz(momPosX, momPosY, momPosZ);
      esdV0->GetNPxPyPz(momNegX, momNegY, momNegZ);
      transvMomV0  = TMath::Sqrt(momV0X*momV0X + momV0Y*momV0Y);
      transvMomPos = TMath::Sqrt(momPosX*momPosX + momPosY*momPosY);
      transvMomNeg = TMath::Sqrt(momNegX*momNegX + momNegY*momNegY);
      etaV0  = esdV0->Eta();
      etaPos = pTrack->Eta();
      etaNeg = nTrack->Eta();
      phiV0  = esdV0->Phi();

    } // end of ESD treatment

    else if(fAnalysisType == "AOD") {

      aodV0 = aodEvent->GetV0(iV0);
      if(!aodV0) continue;

      if(aodV0->GetNDaughters() != 2)                 continue;
      if(aodV0->GetNProngs() != 2)                    continue;
      if(aodV0->GetCharge() != 0)                     continue;
      if(aodV0->ChargeProng(0) == aodV0->ChargeProng(1)) continue;

      if(aodV0->GetOnFlyStatus()) continue; // select offline V0

      // Get the tracks for the daughters
      AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
      AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));

      if(!pTrack || !nTrack) {
        AliWarning("ERROR: Could not retrieve one of the 2 AOD daughter tracks of the V0 ...");
        continue;
      }

      // Rejection of a poor quality tracks
      if(pTrack->GetTPCNcls() < 70) continue;
      if(nTrack->GetTPCNcls() < 70) continue;
      if(TMath::Abs(pTrack->Eta()) > 0.8) continue;
      if(TMath::Abs(nTrack->Eta()) > 0.8) continue;
      if(aodV0->Pt() < 0.3) continue;

      // Out-of-bunch pile-up removal: require either a hit in the ITS SPD or ITS SDD, or TOF in-bunch timing
      if(fPileupCut) {
        if(!(pTrack->HasPointOnITSLayer(0)) && !(pTrack->HasPointOnITSLayer(1)) &&
           !(pTrack->HasPointOnITSLayer(4)) && !(pTrack->HasPointOnITSLayer(5)) &&
           !(pTrack->GetTOFBunchCrossing() == 0)) continue;
        if(!(nTrack->HasPointOnITSLayer(0)) && !(nTrack->HasPointOnITSLayer(1)) &&
           !(nTrack->HasPointOnITSLayer(4)) && !(nTrack->HasPointOnITSLayer(5)) &&
           !(nTrack->GetTOFBunchCrossing() == 0)) continue;
      }

      // Daughter track PID using TPC
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton)) < 5.0) isPosProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion  )) < 5.0) isPosPion   = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton)) < 5.0) isNegProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion  )) < 5.0) isNegPion   = kTRUE;

      // Calculate the invariant mass
      invMassLambda     = aodV0->MassLambda();
      invMassAntiLambda = aodV0->MassAntiLambda();

      // Get topological values
      dcaV0Dghters    = aodV0->DcaV0Daughters();
      dcaPosToPrimVtx = aodV0->DcaPosToPrimVertex();
      dcaNegToPrimVtx = aodV0->DcaNegToPrimVertex();
      cosPointAngle   = aodV0->CosPointingAngle(primaryVtxPos);
      aodV0->GetXYZ(vtxPosV0);
      radius = TMath::Sqrt(vtxPosV0[0]*vtxPosV0[0] + vtxPosV0[1]*vtxPosV0[1]);

      momV0X  = aodV0->MomV0X();
      momV0Y  = aodV0->MomV0Y();
      momV0Z  = aodV0->MomV0Z();
      momPosX = aodV0->MomPosX();
      momPosY = aodV0->MomPosY();
      momPosZ = aodV0->MomPosZ();
      momNegX = aodV0->MomNegX();
      momNegY = aodV0->MomNegY();
      momNegZ = aodV0->MomNegZ();
      transvMomV0  = TMath::Sqrt(momV0X*momV0X + momV0Y*momV0Y);
      transvMomPos = TMath::Sqrt(momPosX*momPosX + momPosY*momPosY);
      transvMomNeg = TMath::Sqrt(momNegX*momNegX + momNegY*momNegY);
      phiV0   = aodV0->Phi();
      etaV0   = aodV0->Eta();
      etaPos  = pTrack->Eta();
      etaNeg  = nTrack->Eta();

    } // end of AOD treatment

    if(isPosProton && isNegPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassLambdawoCuts"))        ->Fill(invMassLambda);
      dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaDCADaughterTracks"))    ->Fill(dcaV0Dghters);
      dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaDCAPosDaughPrimVertex"))->Fill(dcaPosToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaDCANegDaughPrimVertex"))->Fill(dcaNegToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaCosPointingAngle"))     ->Fill(cosPointAngle);
      dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaTransverseRadius"))     ->Fill(radius);
    }
    else if(isNegProton && isPosPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassAntiLambdawoCuts"))->Fill(invMassAntiLambda);
    }

    // Topological cuts
    if(TMath::Abs(vtxPosV0[0]) > 100) continue;
    if(TMath::Abs(vtxPosV0[1]) > 100) continue;
    if(TMath::Abs(vtxPosV0[2]) > 100) continue;
    if(dcaV0Dghters > 1.5) continue;
    if(cosPointAngle < 0.99) continue;
    if(dcaPosToPrimVtx < 0.05) continue;
    if(dcaNegToPrimVtx < 0.05) continue;
    if(radius < 0.2 ) continue;
    if(radius > 100 ) continue;

    // Lambda
    if(isPosProton && isNegPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassLambdawCuts"))->Fill(invMassLambda);

      if(TMath::Abs(invMassLambda - massLambda) < 0.004) {
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaPt"))        ->Fill(transvMomV0);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaPosDaughPt"))->Fill(transvMomPos);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaNegDaughPt"))->Fill(transvMomNeg);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaPhi"))       ->Fill(phiV0);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaEta"))       ->Fill(etaV0);
        dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics")) ->Fill(3);

        if(fAnalysisType == "ESD")      candLambda.Add(esdV0);
        else if(fAnalysisType == "AOD") candLambda.Add(aodV0);

      } else if(invMassLambda - massLambda > 0.004) {

        if(fAnalysisType == "ESD")      rsideLambda.Add(esdV0);
        else if(fAnalysisType == "AOD") rsideLambda.Add(aodV0);

      } else if(invMassLambda - massLambda < -0.004) {

        if(fAnalysisType == "ESD")      lsideLambda.Add(esdV0);
        else if(fAnalysisType == "AOD") lsideLambda.Add(aodV0);

      }
    }
    // AntiLambda
    if(isNegProton && isPosPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassAntiLambdawCuts"))->Fill(invMassAntiLambda);

      if(TMath::Abs(invMassAntiLambda - massLambda) < 0.004) {
        dynamic_cast<TH1F*>(fOutput->FindObject("hAntiLambdaPt"))->Fill(transvMomV0);
        dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(4);

        if(fAnalysisType == "ESD")      candAntiLambda.Add(esdV0);
        else if(fAnalysisType == "AOD") candAntiLambda.Add(aodV0);

      } else if(invMassLambda - massLambda > 0.004) {

        if(fAnalysisType == "ESD")      rsideAntiLambda.Add(esdV0);
        else if(fAnalysisType == "AOD") rsideAntiLambda.Add(aodV0);

      } else if(invMassLambda - massLambda < -0.004) {

        if(fAnalysisType == "ESD")      lsideAntiLambda.Add(esdV0);
        else if(fAnalysisType == "AOD") lsideAntiLambda.Add(aodV0);

      }
    }

  }// end of V0 loop


  // - Cascade loop

  Int_t nCascade = 0;
  if     (fAnalysisType == "ESD") nCascade = esdEvent->GetNumberOfCascades();
  else if(fAnalysisType == "AOD") nCascade = aodEvent->GetNumberOfCascades();

  TObjArray candXim(nCascade);
  TObjArray candXip(nCascade);
  TObjArray rsideXim(nCascade);
  TObjArray rsideXip(nCascade);
  TObjArray lsideXim(nCascade);
  TObjArray lsideXip(nCascade);
  TObjArray candOmegam(nCascade);
  TObjArray candOmegap(nCascade);
  TObjArray rsideOmegam(nCascade);
  TObjArray rsideOmegap(nCascade);
  TObjArray lsideOmegam(nCascade);
  TObjArray lsideOmegap(nCascade);
  
  for(Int_t iXi = 0; iXi < nCascade; iXi++) {// This is the beginning of the Cascade loop

    AliESDcascade *esdXi = 0x0;
    AliAODcascade *aodXi = 0x0;

    // Initialisation of the local variables that will be needed for ESD/AOD
    Double_t invMassXiMinus    = 0.;
    Double_t invMassXiPlus     = 0.;
    Double_t invMassOmegaMinus = 0.;
    Double_t invMassOmegaPlus  = 0.;
    Double_t invMassLambdaAsCascDghter     = 0.;

    Double_t dcaXiDghters     = -1.; 
    Double_t cosPointAngleXi  = -1.;
    Double_t vtxPosXi[3]      = {-999., -999., -999.};
    Double_t radiusXi         = -999.;
    Double_t dcaV0Dghters     = -1.; 
    Double_t dcaV0ToPrimVtx   = -1.;
    Double_t cosPointAngleV0  = -1.;
    Double_t vtxPosV0[3]      = {-999., -999., -999.};
    Double_t radiusV0         = -999.;
    Double_t v0quality        = 0.;
    Double_t dcaBachToPrimVtx = -1.;
    Double_t dcaPosToPrimVtx  = -1.;
    Double_t dcaNegToPrimVtx  = -1.;

    Bool_t isBachelorKaon   = kFALSE;
    Bool_t isBachelorPion   = kFALSE;
    Bool_t isPosProton      = kFALSE;
    Bool_t isPosPion        = kFALSE;
    Bool_t isNegProton      = kFALSE;
    Bool_t isNegPion        = kFALSE;

    Double_t etaXi         = -20.;
    Double_t phiXi         = 720.;
    Double_t momXiX        = 0.;
    Double_t momXiY        = 0.;
    Double_t momXiZ        = 0.;
    Double_t transvMomXi   = 0.;
    Double_t etaBach       = -20.;
    Double_t momBachX      = 0.;
    Double_t momBachY      = 0.;
    Double_t momBachZ      = 0.;
    Double_t transvMomBach = 0.;
    Double_t etaPos        = -20.;
    Double_t momPosX       = 0.;
    Double_t momPosY       = 0.;
    Double_t momPosZ       = 0.;
    Double_t transvMomPos  = 0.;
    Double_t etaNeg        = -20.;
    Double_t momNegX       = 0.;
    Double_t momNegY       = 0.;
    Double_t momNegZ       = 0.;
    Double_t transvMomNeg  = 0.;

    Short_t chargeXi = -2;

    if(fAnalysisType == "ESD") {
      
      esdXi = esdEvent->GetCascade(iXi);
      if(!esdXi) continue;

      chargeXi = esdXi->Charge();

      // Get the tracks for the bachelor and daughters
      AliESDtrack *bachTrack = esdEvent->GetTrack(esdXi->GetBindex());
      AliESDtrack *pTrack    = esdEvent->GetTrack(esdXi->GetPindex());
      AliESDtrack *nTrack    = esdEvent->GetTrack(esdXi->GetNindex());

      if(!bachTrack || !pTrack || !nTrack) {
        AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
        continue;
      }

      // Rejection of a poor quality tracks
      if(bachTrack->GetTPCNcls() < 70) continue;
      if(pTrack   ->GetTPCNcls() < 70) continue;
      if(nTrack   ->GetTPCNcls() < 70) continue;
      if(TMath::Abs(bachTrack->Eta()) > 0.8) continue;
      if(TMath::Abs(pTrack   ->Eta()) > 0.8) continue;
      if(TMath::Abs(nTrack   ->Eta()) > 0.8) continue;
      if(bachTrack->Pt() < 0.3) continue;
      if(pTrack   ->Pt() < 0.3) continue;
      if(nTrack   ->Pt() < 0.3) continue;

      // Bachelor and daughter track PID using TPC
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bachTrack,AliPID::kKaon)) < 4.0) isBachelorKaon = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bachTrack,AliPID::kPion)) < 4.0) isBachelorPion = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton )) < 4.0) isPosProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion   )) < 4.0) isPosPion   = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton )) < 4.0) isNegProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion   )) < 4.0) isNegPion   = kTRUE;

      // Calculate the invariant mass
      invMassLambdaAsCascDghter     = esdXi->GetEffMass();
      v0quality = 0.;
      esdXi->ChangeMassHypothesis(v0quality, 3312);
      invMassXiMinus = esdXi->GetEffMassXi();
      v0quality = 0.;
      esdXi->ChangeMassHypothesis(v0quality, -3312);
      invMassXiPlus = esdXi->GetEffMassXi();
      v0quality = 0.;
      esdXi->ChangeMassHypothesis(v0quality, 3334);
      invMassOmegaMinus = esdXi->GetEffMassXi();
      v0quality = 0.;
      esdXi->ChangeMassHypothesis(v0quality, -3334);
      invMassOmegaPlus = esdXi->GetEffMassXi();

      // Get topological values
      esdXi->GetXYZcascade(vtxPosXi[0], vtxPosXi[1], vtxPosXi[2]);
      esdXi->GetXYZ(vtxPosV0[0], vtxPosV0[1], vtxPosV0[2]);
      dcaXiDghters     = esdXi->GetDcaXiDaughters();
      cosPointAngleXi  = esdXi->GetCascadeCosineOfPointingAngle(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]);
      radiusXi         = TMath::Sqrt(vtxPosXi[0]*vtxPosXi[0] + vtxPosXi[1]*vtxPosXi[1]);
      dcaV0Dghters     = esdXi->GetDcaV0Daughters();
      dcaV0ToPrimVtx   = esdXi->GetD(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]);
      cosPointAngleV0  = esdXi->GetV0CosineOfPointingAngle(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]);
      radiusV0         = TMath::Sqrt(vtxPosV0[0]*vtxPosV0[0] + vtxPosV0[1]*vtxPosV0[1]);
      dcaBachToPrimVtx = TMath::Abs(bachTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));
      dcaPosToPrimVtx  = TMath::Abs(pTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));
      dcaNegToPrimVtx  = TMath::Abs(nTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));

      esdXi->GetPxPyPz(momXiX, momXiY, momXiZ);
      esdXi->GetBPxPyPz(momBachX, momBachY, momBachZ);
      esdXi->GetPPxPyPz(momPosX, momPosY, momPosZ);
      esdXi->GetNPxPyPz(momNegX, momNegY, momNegZ);
      transvMomXi   = TMath::Sqrt(momXiX*momXiX + momXiY*momXiY);
      transvMomBach = TMath::Sqrt(momBachX*momBachX + momBachY*momBachY);
      transvMomPos  = TMath::Sqrt(momPosX*momPosX + momPosY*momPosY);
      transvMomNeg  = TMath::Sqrt(momNegX*momNegX + momNegY*momNegY);
      phiXi   = esdXi->Phi();
      etaXi   = esdXi->Eta();
      etaBach = bachTrack->Eta();
      etaPos  = pTrack->Eta();
      etaNeg  = nTrack->Eta();

    } // end of ESD treatment
    
    else if(fAnalysisType == "AOD") {

      aodXi = aodEvent->GetCascade(iXi);
      if(!aodXi) continue;

      chargeXi = aodXi->ChargeXi();

      // Get the tracks for the bachelor and daughters
      AliAODTrack *bachTrack = dynamic_cast<AliAODTrack*>(aodXi->GetDecayVertexXi()->GetDaughter(0));
      AliAODTrack *pTrack    = dynamic_cast<AliAODTrack*>(aodXi->GetDaughter(0));
      AliAODTrack *nTrack    = dynamic_cast<AliAODTrack*>(aodXi->GetDaughter(1));

      if(!bachTrack || !pTrack || !nTrack) {
        AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
        continue;
      }

      // Rejection of a poor quality tracks
      if(bachTrack->GetTPCNcls() < 70) continue;
      if(pTrack   ->GetTPCNcls() < 70) continue;
      if(nTrack   ->GetTPCNcls() < 70) continue;
      if(TMath::Abs(bachTrack->Eta()) > 0.8) continue;
      if(TMath::Abs(pTrack   ->Eta()) > 0.8) continue;
      if(TMath::Abs(nTrack   ->Eta()) > 0.8) continue;
      if(bachTrack->Pt() < 0.3) continue;
      if(pTrack   ->Pt() < 0.3) continue;
      if(nTrack   ->Pt() < 0.3) continue;

      // Out-of-bunch pile-up removal: require either a hit in the ITS SPD or ITS SDD, or TOF in-bunch timing
      if(fPileupCut) {
        if(!(bachTrack->HasPointOnITSLayer(0)) && !(bachTrack->HasPointOnITSLayer(1)) &&
           !(bachTrack->HasPointOnITSLayer(4)) && !(bachTrack->HasPointOnITSLayer(5)) &&
           !(bachTrack->GetTOFBunchCrossing() == 0)) continue;
        if(!(pTrack->HasPointOnITSLayer(0)) && !(pTrack->HasPointOnITSLayer(1)) &&
           !(pTrack->HasPointOnITSLayer(4)) && !(pTrack->HasPointOnITSLayer(5)) &&
           !(pTrack->GetTOFBunchCrossing() == 0)) continue;
        if(!(nTrack->HasPointOnITSLayer(0)) && !(nTrack->HasPointOnITSLayer(1)) &&
           !(nTrack->HasPointOnITSLayer(4)) && !(nTrack->HasPointOnITSLayer(5)) &&
           !(nTrack->GetTOFBunchCrossing() == 0)) continue;
      }

      // Bachelor and daughter track PID using TPC
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bachTrack,AliPID::kKaon)) < 4.0) isBachelorKaon = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bachTrack,AliPID::kPion)) < 4.0) isBachelorPion = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton )) < 4.0) isPosProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion   )) < 4.0) isPosPion   = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton )) < 4.0) isNegProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion   )) < 4.0) isNegPion   = kTRUE;

      // Calculate the invariant mass
      if(chargeXi < 0.) {
        invMassLambdaAsCascDghter = aodXi->MassLambda();
        invMassXiMinus    = aodXi->MassXi();
        invMassOmegaMinus = aodXi->MassOmega();
      }
      else if(chargeXi > 0.) {
        invMassLambdaAsCascDghter = aodXi->MassAntiLambda();
        invMassXiPlus     = aodXi->MassXi();
        invMassOmegaPlus  = aodXi->MassOmega();
      }

      // Get topological values
      vtxPosXi[0] = aodXi->DecayVertexXiX();
      vtxPosXi[1] = aodXi->DecayVertexXiY();
      vtxPosXi[2] = aodXi->DecayVertexXiZ();
      vtxPosV0[0] = aodXi->DecayVertexV0X();
      vtxPosV0[1] = aodXi->DecayVertexV0Y();
      vtxPosV0[2] = aodXi->DecayVertexV0Z();
      dcaXiDghters     = aodXi->DcaXiDaughters();
      cosPointAngleXi  = aodXi->CosPointingAngleXi(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]);
      radiusXi         = TMath::Sqrt(vtxPosXi[0]*vtxPosXi[0] + vtxPosXi[1]*vtxPosXi[1]);
      dcaV0Dghters     = aodXi->DcaV0Daughters();
      dcaV0ToPrimVtx   = aodXi->DcaV0ToPrimVertex();
      cosPointAngleV0  = aodXi->CosPointingAngle(primaryVtxPos);
      radiusV0         = TMath::Sqrt(vtxPosV0[0]*vtxPosV0[0] + vtxPosV0[1]*vtxPosV0[1]);
      dcaBachToPrimVtx = aodXi->DcaBachToPrimVertex();
      dcaPosToPrimVtx  = aodXi->DcaPosToPrimVertex();
      dcaNegToPrimVtx  = aodXi->DcaNegToPrimVertex();

      momXiX = aodXi->MomXiX();
      momXiY = aodXi->MomXiY();
      momXiZ = aodXi->MomXiZ();
      momBachX = aodXi->MomBachX();
      momBachY = aodXi->MomBachY();
      momBachZ = aodXi->MomBachZ();
      momPosX = aodXi->MomPosX();
      momPosY = aodXi->MomPosY();
      momPosZ = aodXi->MomPosZ();
      momNegX = aodXi->MomNegX();
      momNegY = aodXi->MomNegY();
      momNegZ = aodXi->MomNegZ();
      transvMomXi   = TMath::Sqrt(momXiX*momXiX + momXiY*momXiY);
      transvMomBach = TMath::Sqrt(momBachX*momBachX + momBachY*momBachY);
      transvMomPos  = TMath::Sqrt(momPosX*momPosX + momPosY*momPosY);
      transvMomNeg  = TMath::Sqrt(momNegX*momNegX + momNegY*momNegY);
      phiXi   = aodXi->Phi();
      etaXi   = aodXi->Eta();
      etaBach = bachTrack->Eta();
      etaPos  = pTrack->Eta();
      etaNeg  = nTrack->Eta();

    } // end of AOD treatment

    if((chargeXi<0) && isBachelorPion && isPosProton && isNegPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassXimwoCuts"))          ->Fill(invMassXiMinus);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCADaughterTracks"))       ->Fill(dcaXiDghters);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCAV0DaughterTracks"))     ->Fill(dcaV0Dghters);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCAV0PrimVertex"))         ->Fill(dcaV0ToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCABachPrimVertex"))       ->Fill(dcaBachToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCAPosDaughPrimVertex"))   ->Fill(dcaPosToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCANegDaughPrimVertex"))   ->Fill(dcaNegToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiCosPointingAngle"))        ->Fill(cosPointAngleXi);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiV0CosPointingAngle"))      ->Fill(cosPointAngleV0);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiTransverseRadius"))        ->Fill(radiusXi);
      dynamic_cast<TH1F*>(fOutput->FindObject("hXiV0TransverseRadius"))      ->Fill(radiusV0);
    }
    if((chargeXi>0) && isBachelorPion && isNegProton && isPosPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassXipwoCuts"))          ->Fill(invMassXiPlus);
    }
    if((chargeXi<0) && isBachelorKaon && isPosProton && isNegPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassOmegamwoCuts"))       ->Fill(invMassOmegaMinus);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaDCADaughterTracks"))    ->Fill(dcaXiDghters);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaDCAV0DaughterTracks"))  ->Fill(dcaV0Dghters);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaDCAV0PrimVertex"))      ->Fill(dcaV0ToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaDCABachPrimVertex"))    ->Fill(dcaBachToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaDCAPosDaughPrimVertex"))->Fill(dcaPosToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaDCANegDaughPrimVertex"))->Fill(dcaNegToPrimVtx);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaCosPointingAngle"))     ->Fill(cosPointAngleXi);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaV0CosPointingAngle"))   ->Fill(cosPointAngleV0);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaTransverseRadius"))     ->Fill(radiusXi);
      dynamic_cast<TH1F*>(fOutput->FindObject("hOmegaV0TransverseRadius"))   ->Fill(radiusV0);
    }
    if((chargeXi>0) && isBachelorKaon && isNegProton && isPosPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassOmegapwoCuts"))       ->Fill(invMassOmegaPlus);
    }

    // Topological cuts
    if(dcaXiDghters > 1.6) continue;
    if(dcaV0Dghters > 1.6) continue;
    if(dcaV0Dghters < 0.07) continue;
    if(dcaBachToPrimVtx < 0.05) continue;
    if(dcaPosToPrimVtx < 0.04) continue;
    if(dcaNegToPrimVtx < 0.04) continue;
    if(cosPointAngleXi < 0.97) continue;
    if(cosPointAngleV0 < 0.97) continue;
    if(radiusXi < 0.8) continue;
    if(radiusXi > 200) continue;
    if(radiusV0 < 1.4) continue;
    if(radiusV0 > 200) continue;

    dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassLambdaAsCascDghter"))->Fill(invMassLambdaAsCascDghter);

    if(TMath::Abs(invMassLambdaAsCascDghter - massLambda) < 0.006) {

      // Xi Minus
      if((chargeXi<0) && isBachelorPion && isPosProton && isNegPion) {

        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassXimwCuts"))->Fill(invMassXiMinus);

        if(TMath::Abs(invMassXiMinus - massXi) < 0.008) {
          dynamic_cast<TH1F*>(fOutput->FindObject("hXimPt"))        ->Fill(transvMomXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXimBachPt"))    ->Fill(transvMomBach);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXimPosDaughPt"))->Fill(transvMomPos);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXimNegDaughPt"))->Fill(transvMomNeg);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXimPhi"))       ->Fill(phiXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXimEta"))       ->Fill(etaXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(5);

          if(fAnalysisType == "ESD")      candXim.Add(esdXi);
          else if(fAnalysisType == "AOD") candXim.Add(aodXi);

        } else if(invMassXiMinus - massXi > 0.008) {

          if(fAnalysisType == "ESD")      rsideXim.Add(esdXi);
          else if(fAnalysisType == "AOD") rsideXim.Add(aodXi);

        } else if(invMassXiMinus - massXi < -0.008) {

          if(fAnalysisType == "ESD")      lsideXim.Add(esdXi);
          else if(fAnalysisType == "AOD") lsideXim.Add(aodXi);

        }
      }
      // Xi Plus
      if((chargeXi>0) && isBachelorPion && isNegProton && isPosPion) {

        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassXipwCuts"))->Fill(invMassXiPlus);

        if(TMath::Abs(invMassXiPlus - massXi) < 0.008) {
          dynamic_cast<TH1F*>(fOutput->FindObject("hXipPt"))->Fill(transvMomXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(6);

          if(fAnalysisType == "ESD")      candXip.Add(esdXi);
          else if(fAnalysisType == "AOD") candXip.Add(aodXi);

        } else if(invMassXiPlus - massXi > 0.008) {

          if(fAnalysisType == "ESD")      rsideXip.Add(esdXi);
          else if(fAnalysisType == "AOD") rsideXip.Add(aodXi);

        } else if(invMassXiPlus - massXi < -0.008) {

          if(fAnalysisType == "ESD")      lsideXip.Add(esdXi);
          else if(fAnalysisType == "AOD") lsideXip.Add(aodXi);

        }
      }
      // Omega Minus
      if((chargeXi<0) && isBachelorKaon && isPosProton && isNegPion) {

        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassOmegamwCuts"))->Fill(invMassOmegaMinus);

        if(TMath::Abs(invMassOmegaMinus - massOmega) < 0.008) {
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamPt"))        ->Fill(transvMomXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamBachPt"))    ->Fill(transvMomBach);
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamPosDaughPt"))->Fill(transvMomPos);
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamNegDaughPt"))->Fill(transvMomNeg);
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamPhi"))       ->Fill(phiXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamEta"))       ->Fill(etaXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(7);

          if(fAnalysisType == "ESD")      candOmegam.Add(esdXi);
          else if(fAnalysisType == "AOD") candOmegam.Add(aodXi);

        } else if(invMassOmegaMinus - massOmega > 0.008) {

          if(fAnalysisType == "ESD")      rsideOmegam.Add(esdXi);
          else if(fAnalysisType == "AOD") rsideOmegam.Add(aodXi);

        } else if(invMassOmegaMinus - massOmega < -0.008) {

          if(fAnalysisType == "ESD")      lsideOmegam.Add(esdXi);
          else if(fAnalysisType == "AOD") lsideOmegam.Add(aodXi);

        }
      }
      // Omega Plus
      if((chargeXi>0) && isBachelorKaon && isNegProton && isPosPion) {

        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassOmegapwCuts"))->Fill(invMassOmegaPlus);

        if(TMath::Abs(invMassOmegaPlus - massOmega) < 0.008) {
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegapPt"))->Fill(transvMomXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(8);

          if(fAnalysisType == "ESD")      candOmegam.Add(esdXi);
          else if(fAnalysisType == "AOD") candOmegam.Add(aodXi);

        } else if(invMassOmegaPlus - massOmega > 0.008) {

          if(fAnalysisType == "ESD")      rsideOmegap.Add(esdXi);
          else if(fAnalysisType == "AOD") rsideOmegap.Add(aodXi);

        } else if(invMassOmegaPlus - massOmega < -0.008) {

          if(fAnalysisType == "ESD")      lsideOmegap.Add(esdXi);
          else if(fAnalysisType == "AOD") lsideOmegap.Add(aodXi);

        }
      }
    }

  } // end of cascade loop


  //______________________________________________________________________________
  // Calculate invariant mass for dibaryons

  if(fAnalysisType == "ESD") {

    // ppK- -> proton + Lambda
    for(Int_t i=0; i<candProton.GetEntriesFast(); i++) {
      AliESDtrack *track = (AliESDtrack*)candProton.At(i);

      for(Int_t j=0; j<candLambda.GetEntriesFast(); j++) { // signal
        AliESDv0 *v0 = (AliESDv0*)candLambda.At(j);

        if(track->GetID() == v0->GetPindex()) continue;
        if(track->GetID() == v0->GetNindex()) continue;

        TLorentzVector proton, lambda, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);

        trackSum = proton + lambda;
        trackRel = proton - lambda;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(1);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonLambda"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideLambda.GetEntriesFast(); j++) { // right sideband
        AliESDv0 *v0 = (AliESDv0*)rsideLambda.At(j);

        if(track->GetID() == v0->GetPindex()) continue;
        if(track->GetID() == v0->GetNindex()) continue;

        TLorentzVector proton, lambda, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);

        trackSum = proton + lambda;
        trackRel = proton - lambda;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonLambdaRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideLambda.GetEntriesFast(); j++) { // left sideband
        AliESDv0 *v0 = (AliESDv0*)lsideLambda.At(j);

        if(track->GetID() == v0->GetPindex()) continue;
        if(track->GetID() == v0->GetNindex()) continue;

        TLorentzVector proton, lambda, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);

        trackSum = proton + lambda;
        trackRel = proton - lambda;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonLambdaLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }

    // H-Dibaryon -> Lambda + Lambda
    for(Int_t i=0; i<candLambda.GetEntriesFast(); i++) {
      AliESDv0 *v01 = (AliESDv0*)candLambda.At(i);

      for(Int_t j=i+1; j<candLambda.GetEntriesFast(); j++) { // signal
        AliESDv0 *v02 = (AliESDv0*)candLambda.At(j);

        if(v01->GetPindex() == v02->GetPindex()) continue;
        if(v01->GetNindex() == v02->GetNindex()) continue;

        TLorentzVector lambda1, lambda2, trackSum, trackRel;

        lambda1.SetXYZM(v01->Px(), v01->Py(), v01->Pz(), massLambda);
        lambda2.SetXYZM(v02->Px(), v02->Py(), v02->Pz(), massLambda);

        trackSum = lambda1 + lambda2;
        trackRel = lambda1 - lambda2;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(2);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaLambda"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideLambda.GetEntriesFast(); j++) { // right sideband
        AliESDv0 *v02 = (AliESDv0*)rsideLambda.At(j);

        if(v01->GetPindex() == v02->GetPindex()) continue;
        if(v01->GetNindex() == v02->GetNindex()) continue;

        TLorentzVector lambda1, lambda2, trackSum, trackRel;

        lambda1.SetXYZM(v01->Px(), v01->Py(), v01->Pz(), massLambda);
        lambda2.SetXYZM(v02->Px(), v02->Py(), v02->Pz(), massLambda);

        trackSum = lambda1 + lambda2;
        trackRel = lambda1 - lambda2;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaLambdaRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideLambda.GetEntriesFast(); j++) { // left sideband
        AliESDv0 *v02 = (AliESDv0*)lsideLambda.At(j);

        if(v01->GetPindex() == v02->GetPindex()) continue;
        if(v01->GetNindex() == v02->GetNindex()) continue;

        TLorentzVector lambda1, lambda2, trackSum, trackRel;

        lambda1.SetXYZM(v01->Px(), v01->Py(), v01->Pz(), massLambda);
        lambda2.SetXYZM(v02->Px(), v02->Py(), v02->Pz(), massLambda);

        trackSum = lambda1 + lambda2;
        trackRel = lambda1 - lambda2;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaLambdaLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }

    // H-Dibaryon -> proton + Xi-
    for(Int_t i=0; i<candProton.GetEntriesFast(); i++) {
      AliESDtrack *track = (AliESDtrack*)candProton.At(i);

      for(Int_t j=0; j<candXim.GetEntriesFast(); j++) { // signal
        AliESDcascade *xi = (AliESDcascade*)candXim.At(j);

        if(track->GetID() == xi->GetBindex()) continue;
        if(track->GetID() == xi->GetPindex()) continue;
        if(track->GetID() == xi->GetNindex()) continue;

        TLorentzVector proton, xim, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = proton + xim;
        trackRel = proton - xim;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(3);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonXi"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideXim.GetEntriesFast(); j++) { // right sideband
        AliESDcascade *xi = (AliESDcascade*)rsideXim.At(j);

        if(track->GetID() == xi->GetBindex()) continue;
        if(track->GetID() == xi->GetPindex()) continue;
        if(track->GetID() == xi->GetNindex()) continue;

        TLorentzVector proton, xim, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = proton + xim;
        trackRel = proton - xim;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonXiRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideXim.GetEntriesFast(); j++) { // left sideband
        AliESDcascade *xi = (AliESDcascade*)lsideXim.At(j);

        if(track->GetID() == xi->GetBindex()) continue;
        if(track->GetID() == xi->GetPindex()) continue;
        if(track->GetID() == xi->GetNindex()) continue;

        TLorentzVector proton, xim, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = proton + xim;
        trackRel = proton - xim;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonXiLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }

    // nOmega- -> Lambda + Xi-
    for(Int_t i=0; i<candLambda.GetEntriesFast(); i++) {
      AliESDv0 *v0 = (AliESDv0*)candLambda.At(i);

      for(Int_t j=0; j<candXim.GetEntriesFast(); j++) { // signal
        AliESDcascade *xi = (AliESDcascade*)candXim.At(j);

        if(v0->GetPindex() == xi->GetPindex()) continue;
        if(v0->GetNindex() == xi->GetNindex()) continue;
        if(v0->GetNindex() == xi->GetBindex()) continue;

        TLorentzVector lambda, xim, trackSum, trackRel;

        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = lambda + xim;
        trackRel = lambda - xim;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(4);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaXi"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideXim.GetEntriesFast(); j++) { // right sideband
        AliESDcascade *xi = (AliESDcascade*)rsideXim.At(j);

        if(v0->GetPindex() == xi->GetPindex()) continue;
        if(v0->GetNindex() == xi->GetNindex()) continue;
        if(v0->GetNindex() == xi->GetBindex()) continue;

        TLorentzVector lambda, xim, trackSum, trackRel;

        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = lambda + xim;
        trackRel = lambda - xim;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaXiRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideXim.GetEntriesFast(); j++) { // left sideband
        AliESDcascade *xi = (AliESDcascade*)lsideXim.At(j);

        if(v0->GetPindex() == xi->GetPindex()) continue;
        if(v0->GetNindex() == xi->GetNindex()) continue;
        if(v0->GetNindex() == xi->GetBindex()) continue;

        TLorentzVector lambda, xim, trackSum, trackRel;

        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = lambda + xim;
        trackRel = lambda - xim;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaXiLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }

    // Di-Omega -> Xi- + Omega-
    for(Int_t i=0; i<candXim.GetEntriesFast(); i++) {
      AliESDcascade *xi1 = (AliESDcascade*)candXim.At(i);

      for(Int_t j=0; j<candOmegam.GetEntriesFast(); j++) { // signal
        AliESDcascade *xi2 = (AliESDcascade*)candOmegam.At(j);

        if(xi1->GetPindex() == xi2->GetPindex()) continue;
        if(xi1->GetNindex() == xi2->GetNindex()) continue;
        if(xi1->GetNindex() == xi2->GetBindex()) continue;
        if(xi1->GetBindex() == xi2->GetNindex()) continue;

        TLorentzVector xim, omegam, trackSum, trackRel;

        xim.SetXYZM(xi1->Px(), xi1->Py(), xi1->Pz(), massXi);
        omegam.SetXYZM(xi2->Px(), xi2->Py(), xi2->Pz(), massOmega);

        trackSum = xim + omegam;
        trackRel = xim - omegam;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(5);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPXiOmega"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideOmegam.GetEntriesFast(); j++) { // right sideband
        AliESDcascade *xi2 = (AliESDcascade*)rsideOmegam.At(j);

        if(xi1->GetPindex() == xi2->GetPindex()) continue;
        if(xi1->GetNindex() == xi2->GetNindex()) continue;
        if(xi1->GetNindex() == xi2->GetBindex()) continue;
        if(xi1->GetBindex() == xi2->GetNindex()) continue;

        TLorentzVector xim, omegam, trackSum, trackRel;

        xim.SetXYZM(xi1->Px(), xi1->Py(), xi1->Pz(), massXi);
        omegam.SetXYZM(xi2->Px(), xi2->Py(), xi2->Pz(), massOmega);

        trackSum = xim + omegam;
        trackRel = xim - omegam;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPXiOmegaRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideOmegam.GetEntriesFast(); j++) { // left sideband
        AliESDcascade *xi2 = (AliESDcascade*)lsideOmegam.At(j);

        if(xi1->GetPindex() == xi2->GetPindex()) continue;
        if(xi1->GetNindex() == xi2->GetNindex()) continue;
        if(xi1->GetNindex() == xi2->GetBindex()) continue;
        if(xi1->GetBindex() == xi2->GetNindex()) continue;

        TLorentzVector xim, omegam, trackSum, trackRel;

        xim.SetXYZM(xi1->Px(), xi1->Py(), xi1->Pz(), massXi);
        omegam.SetXYZM(xi2->Px(), xi2->Py(), xi2->Pz(), massOmega);

        trackSum = xim + omegam;
        trackRel = xim - omegam;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPXiOmegaLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }
  }

  else if(fAnalysisType == "AOD") {

    // ppK- -> proton + Lambda
    for(Int_t i=0; i<candProton.GetEntriesFast(); i++) {
      AliAODTrack *track = (AliAODTrack*)candProton.At(i);

      for(Int_t j=0; j<candLambda.GetEntriesFast(); j++) { // signal
        AliAODv0 *v0 = (AliAODv0*)candLambda.At(j);

        if(track->GetID() == v0->GetPosID()) continue;
        if(track->GetID() == v0->GetNegID()) continue;

        TLorentzVector proton, lambda, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);

        trackSum = proton + lambda;
        trackRel = proton - lambda;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(1);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonLambda"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideLambda.GetEntriesFast(); j++) { // right sideband
        AliAODv0 *v0 = (AliAODv0*)rsideLambda.At(j);

        if(track->GetID() == v0->GetPosID()) continue;
        if(track->GetID() == v0->GetNegID()) continue;

        TLorentzVector proton, lambda, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);

        trackSum = proton + lambda;
        trackRel = proton - lambda;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonLambdaRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideLambda.GetEntriesFast(); j++) { // left sideband
        AliAODv0 *v0 = (AliAODv0*)lsideLambda.At(j);

        if(track->GetID() == v0->GetPosID()) continue;
        if(track->GetID() == v0->GetNegID()) continue;

        TLorentzVector proton, lambda, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);

        trackSum = proton + lambda;
        trackRel = proton - lambda;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonLambdaLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }

    // H-Dibaryon -> Lambda + Lambda
    for(Int_t i=0; i<candLambda.GetEntriesFast(); i++) {
      AliAODv0 *v01 = (AliAODv0*)candLambda.At(i);

      for(Int_t j=i+1; j<candLambda.GetEntriesFast(); j++) { // signal
        AliAODv0 *v02 = (AliAODv0*)candLambda.At(j);

        if(v01->GetPosID() == v02->GetPosID()) continue;
        if(v01->GetNegID() == v02->GetNegID()) continue;

        TLorentzVector lambda1, lambda2, trackSum, trackRel;

        lambda1.SetXYZM(v01->Px(), v01->Py(), v01->Pz(), massLambda);
        lambda2.SetXYZM(v02->Px(), v02->Py(), v02->Pz(), massLambda);

        trackSum = lambda1 + lambda2;
        trackRel = lambda1 - lambda2;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(2);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaLambda"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideLambda.GetEntriesFast(); j++) { // right sideband
        AliAODv0 *v02 = (AliAODv0*)rsideLambda.At(j);

        if(v01->GetPosID() == v02->GetPosID()) continue;
        if(v01->GetNegID() == v02->GetNegID()) continue;

        TLorentzVector lambda1, lambda2, trackSum, trackRel;

        lambda1.SetXYZM(v01->Px(), v01->Py(), v01->Pz(), massLambda);
        lambda2.SetXYZM(v02->Px(), v02->Py(), v02->Pz(), massLambda);

        trackSum = lambda1 + lambda2;
        trackRel = lambda1 - lambda2;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaLambdaRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideLambda.GetEntriesFast(); j++) { // left sideband
        AliAODv0 *v02 = (AliAODv0*)lsideLambda.At(j);

        if(v01->GetPosID() == v02->GetPosID()) continue;
        if(v01->GetNegID() == v02->GetNegID()) continue;

        TLorentzVector lambda1, lambda2, trackSum, trackRel;

        lambda1.SetXYZM(v01->Px(), v01->Py(), v01->Pz(), massLambda);
        lambda2.SetXYZM(v02->Px(), v02->Py(), v02->Pz(), massLambda);

        trackSum = lambda1 + lambda2;
        trackRel = lambda1 - lambda2;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaLambdaLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }

    // H-Dibaryon -> proton + Xi-
    for(Int_t i=0; i<candProton.GetEntriesFast(); i++) {
      AliAODTrack *track = (AliAODTrack*)candProton.At(i);

      for(Int_t j=0; j<candXim.GetEntriesFast(); j++) { // signal
        AliAODcascade *xi = (AliAODcascade*)candXim.At(j);

        if(track->GetID() == xi->GetBachID()) continue;
        if(track->GetID() == xi->GetPosID()) continue;
        if(track->GetID() == xi->GetNegID()) continue;

        TLorentzVector proton, xim, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = proton + xim;
        trackRel = proton - xim;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(3);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonXi"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideXim.GetEntriesFast(); j++) { // right sideband
        AliAODcascade *xi = (AliAODcascade*)rsideXim.At(j);

        if(track->GetID() == xi->GetBachID()) continue;
        if(track->GetID() == xi->GetPosID()) continue;
        if(track->GetID() == xi->GetNegID()) continue;

        TLorentzVector proton, xim, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = proton + xim;
        trackRel = proton - xim;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonXiRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideXim.GetEntriesFast(); j++) { // left sideband
        AliAODcascade *xi = (AliAODcascade*)lsideXim.At(j);

        if(track->GetID() == xi->GetBachID()) continue;
        if(track->GetID() == xi->GetPosID()) continue;
        if(track->GetID() == xi->GetNegID()) continue;

        TLorentzVector proton, xim, trackSum, trackRel;

        proton.SetXYZM(track->Px(), track->Py(), track->Pz(), massProton);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = proton + xim;
        trackRel = proton - xim;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPProtonXiLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }

    // nOmega- -> Lambda + Xi-
    for(Int_t i=0; i<candLambda.GetEntriesFast(); i++) {
      AliAODv0 *v0 = (AliAODv0*)candLambda.At(i);

      for(Int_t j=0; j<candXim.GetEntriesFast(); j++) {
        AliAODcascade *xi = (AliAODcascade*)candXim.At(j);

        if(v0->GetPosID() == xi->GetPosID()) continue;
        if(v0->GetNegID() == xi->GetNegID()) continue;
        if(v0->GetNegID() == xi->GetBachID()) continue;

        TLorentzVector lambda, xim, trackSum, trackRel;

        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = lambda + xim;
        trackRel = lambda - xim;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(4);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaXi"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideXim.GetEntriesFast(); j++) { // right sideband
        AliAODcascade *xi = (AliAODcascade*)rsideXim.At(j);

        if(v0->GetPosID() == xi->GetPosID()) continue;
        if(v0->GetNegID() == xi->GetNegID()) continue;
        if(v0->GetNegID() == xi->GetBachID()) continue;

        TLorentzVector lambda, xim, trackSum, trackRel;

        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = lambda + xim;
        trackRel = lambda - xim;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaXiRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideXim.GetEntriesFast(); j++) { // left sideband
        AliAODcascade *xi = (AliAODcascade*)lsideXim.At(j);

        if(v0->GetPosID() == xi->GetPosID()) continue;
        if(v0->GetNegID() == xi->GetNegID()) continue;
        if(v0->GetNegID() == xi->GetBachID()) continue;

        TLorentzVector lambda, xim, trackSum, trackRel;

        lambda.SetXYZM(v0->Px(), v0->Py(), v0->Pz(), massLambda);
        xim.SetXYZM(xi->Px(), xi->Py(), xi->Pz(), massXi);

        trackSum = lambda + xim;
        trackRel = lambda - xim;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPLambdaXiLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }

    // Di-Omega -> Xi- + Omega-
    for(Int_t i=0; i<candXim.GetEntriesFast(); i++) {
      AliAODcascade *xi1 = (AliAODcascade*)candXim.At(i);

      for(Int_t j=0; j<candOmegam.GetEntriesFast(); j++) { // signal
        AliAODcascade *xi2 = (AliAODcascade*)candOmegam.At(j);

        if(xi1->GetPosID() == xi2->GetPosID()) continue;
        if(xi1->GetNegID() == xi2->GetNegID()) continue;
        if(xi1->GetNegID() == xi2->GetBachID()) continue;
        if(xi1->GetBachID() == xi2->GetNegID()) continue;

        TLorentzVector xim, omegam, trackSum, trackRel;

        xim.SetXYZM(xi1->Px(), xi1->Py(), xi1->Pz(), massXi);
        omegam.SetXYZM(xi2->Px(), xi2->Py(), xi2->Pz(), massOmega);

        trackSum = xim + omegam;
        trackRel = xim - omegam;

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(5);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPXiOmega"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<rsideOmegam.GetEntriesFast(); j++) { // right sideband
        AliAODcascade *xi2 = (AliAODcascade*)rsideOmegam.At(j);

        if(xi1->GetPosID() == xi2->GetPosID()) continue;
        if(xi1->GetNegID() == xi2->GetNegID()) continue;
        if(xi1->GetNegID() == xi2->GetBachID()) continue;
        if(xi1->GetBachID() == xi2->GetNegID()) continue;

        TLorentzVector xim, omegam, trackSum, trackRel;

        xim.SetXYZM(xi1->Px(), xi1->Py(), xi1->Pz(), massXi);
        omegam.SetXYZM(xi2->Px(), xi2->Py(), xi2->Pz(), massOmega);

        trackSum = xim + omegam;
        trackRel = xim - omegam;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPXiOmegaRsideband"))->Fill(trackSum.M(), trackRel.P());

      }

      for(Int_t j=0; j<lsideOmegam.GetEntriesFast(); j++) { // left sideband
        AliAODcascade *xi2 = (AliAODcascade*)lsideOmegam.At(j);

        if(xi1->GetPosID() == xi2->GetPosID()) continue;
        if(xi1->GetNegID() == xi2->GetNegID()) continue;
        if(xi1->GetNegID() == xi2->GetBachID()) continue;
        if(xi1->GetBachID() == xi2->GetNegID()) continue;

        TLorentzVector xim, omegam, trackSum, trackRel;

        xim.SetXYZM(xi1->Px(), xi1->Py(), xi1->Pz(), massXi);
        omegam.SetXYZM(xi2->Px(), xi2->Py(), xi2->Pz(), massOmega);

        trackSum = xim + omegam;
        trackRel = xim - omegam;

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelPXiOmegaLsideband"))->Fill(trackSum.M(), trackRel.P());

      }
    }
  }

}
//_______________________________________________________________________________________________
void AliAnalysisTaskDibaryons::Terminate(Option_t *option)
{
  //Called once at the end of the query
  AliInfo(Form("%s is done.",GetName()));
}
//_______________________________________________________________________________________________
