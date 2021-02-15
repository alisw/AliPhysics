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
  fPIDResponse(0x0),
  fFilterBit(0),
  fPileupCut(kTRUE),
  fEventMixing(kTRUE),
  fOutput(0x0),
  fTrackArray(0x0),
  fProtonArray(0x0),
  fLambdaArray(0x0),
  fXiArray(0x0),
  fOmegaArray(0x0),
  fTrackBuffSize(2500),
  fV0BuffSize(100),
  fCascadeBuffSize(100),
  fProtonEMpool(),
  fLambdaEMpool(),
  fXiEMpool(),
  fOmegaEMpool()
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
  fPIDResponse(0x0),
  fFilterBit(0),
  fPileupCut(kTRUE),
  fEventMixing(kTRUE),
  fOutput(0x0),
  fTrackArray(0x0),
  fProtonArray(0x0),
  fLambdaArray(0x0),
  fXiArray(0x0),
  fOmegaArray(0x0),
  fTrackBuffSize(2500),
  fV0BuffSize(100),
  fCascadeBuffSize(100),
  fProtonEMpool(),
  fLambdaEMpool(),
  fXiEMpool(),
  fOmegaEMpool()
{
  // constructor

  DefineInput(0,TChain::Class());
  DefineOutput(1,THashList::Class());
}
//_______________________________________________________________________________________________
AliAnalysisTaskDibaryons::~AliAnalysisTaskDibaryons()
{
  // destructor

  delete fPIDResponse;
  delete fOutput;
  delete fTrackArray;
  delete fProtonArray;
  delete fLambdaArray;
  delete fXiArray;
  delete fOmegaArray;
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

  if(fAnalysisType == "ESD") {
    fProtonArray = new TClonesArray("AliESDtrack",fTrackBuffSize);
    fLambdaArray = new TClonesArray("AliESDv0",fV0BuffSize);
    fXiArray     = new TClonesArray("AliESDcascade",fCascadeBuffSize);
    fOmegaArray  = new TClonesArray("AliESDcascade",fCascadeBuffSize);
  }
  else if(fAnalysisType == "AOD") {
    fTrackArray  = new AliAODTrack*[fTrackBuffSize];
    fProtonArray = new TClonesArray("AliAODTrack",fTrackBuffSize);
    fLambdaArray = new TClonesArray("AliAODv0",fV0BuffSize);
    fXiArray     = new TClonesArray("AliAODcascade",fCascadeBuffSize);
    fOmegaArray  = new TClonesArray("AliAODcascade",fCascadeBuffSize);
  }

  fProtonEMpool.resize(10,std::vector<EventPool>(26));
  fLambdaEMpool.resize(10,std::vector<EventPool>(26));
  fXiEMpool.resize(10,std::vector<EventPool>(26));
  fOmegaEMpool.resize(10,std::vector<EventPool>(26));

  TH1F *hRefMultiplicity08 = new TH1F("hRefMultiplicity08","Reference multiplicity;Multipicity;Counts",600,0,600);
  TH1F* hNBinsMultMixing = new TH1F("hNBinsMultMixing","Bins in multiplicity that are used for event mixing",26,0,26);
  TH1F* hNBinsVertexMixing = new TH1F("hNBinsVertexMixing","Bins in z-Vertex that are used for event mixing",10,0,10);
  TH1F *hNPartStatistics = new TH1F("hNPartStatistics","Number of candidates under certain condition",10,0.5,10.5);
  hNPartStatistics->GetXaxis()->SetBinLabel(1,"p");
  hNPartStatistics->GetXaxis()->SetBinLabel(2,"#bar{p}");
  hNPartStatistics->GetXaxis()->SetBinLabel(3,"#Lambda");
  hNPartStatistics->GetXaxis()->SetBinLabel(4,"#bar{#Lambda}");
  hNPartStatistics->GetXaxis()->SetBinLabel(5,"#Xi^{-}");
  hNPartStatistics->GetXaxis()->SetBinLabel(6,"#Xi^{+}");
  hNPartStatistics->GetXaxis()->SetBinLabel(7,"#Omega^{-}");
  hNPartStatistics->GetXaxis()->SetBinLabel(8,"#Omega^{+}");
  hNPartStatistics->GetYaxis()->SetTitle("Counts");

  fOutput->Add(hRefMultiplicity08);
  fOutput->Add(hNBinsMultMixing);
  fOutput->Add(hNBinsVertexMixing);
  fOutput->Add(hNPartStatistics);

  // Define histograms related to invariant mass for V0 and Cascade
  TH2F *hInvMassLambda = new TH2F("hInvMassLambda","Invariant mass of p#pi^{-};M_{p#pi} (GeV/c^{2});p_{T} (GeV/c)",400,1.0,1.2,100,0,10);
  TH2F *hInvMassAntiLambda = new TH2F("hInvMassAntiLambda","Invariant mass of #bar{p}#pi^{+};M_{p#pi} (GeV/c^{2});p_{T} (GeV/c)",400,1.0,1.2,100,0,10);
  TH2F *hInvMassXim = new TH2F("hInvMassXim","Invariant mass of #pi^{-}#Lambda;M_{#pi#Lambda} (GeV/c^{2});p_{T} (GeV/c)",300,1.2,1.5,100,0,10);
  TH2F *hInvMassXip = new TH2F("hInvMassXip","Invariant mass of #pi^{+}#bar{#Lambda};M_{#pi#Lambda} (GeV/c^{2});p_{T} (GeV/c)",300,1.2,1.5,100,0,10);
  TH2F *hInvMassOmegam = new TH2F("hInvMassOmegam","Invariant mass of K^{-}#Lambda;M_{K#Lambda} (GeV/c^{2});p_{T} (GeV/c)",300,1.5,1.8,100,0,10);
  TH2F *hInvMassOmegap = new TH2F("hInvMassOmegap","Invariant mass of K^{+}#bar{#Lambda};M_{K#Lambda} (GeV/c^{2});p_{T} (GeV/c)",300,1.5,1.8,100,0,10);
  TH1F *hInvMassLambdaAsCascDghter = new TH1F("hInvMassLambdaAsCascDghter","Invariant mass of p#pi stemming from cascade;M_{p#pi} (GeV/c^{2});Counts",400,1.0,1.2);
  TH1F *hInvMassLambdawoCuts = new TH1F("hInvMassLambdawoCuts","Invariant mass of p#pi^{-} without topological cuts;M_{p#pi} (GeV/c^{2});Counts",400,1.0,1.2);
  TH1F *hInvMassAntiLambdawoCuts = new TH1F("hInvMassAntiLambdawoCuts","Invariant mass of #bar{p}#pi^{+} without topological cuts;M_{p#pi} (GeV/c^{2});Counts",400,1.0,1.2);
  TH1F *hInvMassXimwoCuts = new TH1F("hInvMassXimwoCuts","Invariant mass of #pi^{-}#Lambda without topological cuts;M_{#pi#Lambda} (GeV/c^{2});Counts",300,1.2,1.5);
  TH1F *hInvMassXipwoCuts = new TH1F("hInvMassXipwoCuts","Invariant mass of #pi^{+}#bar{#Lambda} without topological cuts;M_{#pi#Lambda} (GeV/c^{2});Counts",300,1.2,1.5);
  TH1F *hInvMassOmegamwoCuts = new TH1F("hInvMassOmegamwoCuts","Invariant mass of K^{-}#Lambda without topological cuts;M_{K#Lambda} (GeV/c^{2});Counts",300,1.5,1.8);
  TH1F *hInvMassOmegapwoCuts = new TH1F("hInvMassOmegapwoCuts","Invariant mass of K^{+}#bar{#Lambda} without topological cuts;M_{K#Lambda} (GeV/c^{2});Counts",300,1.5,1.8);
  TH1F *hInvMassMissIDK0S = new TH1F("hInvMassMissIDK0S","Invariant mass of #pi^{+}#pi^{-} under assumption of selected #Lambda;M_{#pi#pi} (GeV/c^{2})",400,0.4,0.6);
  TH1F *hInvMassMissIDOmegam = new TH1F("hInvMassMissIDOmegam","Invariant mass of K^{-}#Lambda under assumption of selected #Xi^{-};M_{K#Lambda} (GeV/c^{2})",300,1.5,1.8);
  TH1F *hInvMassMissIDXim = new TH1F("hInvMassMissIDXim","Invariant mass of #pi^{-}#Lambda under assumption of selected #Omega^{-};M_{#pi#Lambda} (GeV/c^{2})",300,1.2,1.5);

  fOutput->Add(hInvMassLambda);
  fOutput->Add(hInvMassAntiLambda);
  fOutput->Add(hInvMassXim);
  fOutput->Add(hInvMassXip);
  fOutput->Add(hInvMassOmegam);
  fOutput->Add(hInvMassOmegap);
  fOutput->Add(hInvMassLambdaAsCascDghter);
  fOutput->Add(hInvMassLambdawoCuts);
  fOutput->Add(hInvMassAntiLambdawoCuts);
  fOutput->Add(hInvMassXimwoCuts);
  fOutput->Add(hInvMassXipwoCuts);
  fOutput->Add(hInvMassOmegamwoCuts);
  fOutput->Add(hInvMassOmegapwoCuts);
  fOutput->Add(hInvMassMissIDK0S);
  fOutput->Add(hInvMassMissIDOmegam);
  fOutput->Add(hInvMassMissIDXim);

  // Define QA plots for topological observables
  TH2F *hProtonDCAxyDCAz = new TH2F("hProtonDCAxyDCAz","Proton DCAxy vs DCAz;DCA_{xy} (cm);DCA_{z} (cm)",500,-5,5,1000,-20,20);
  TH1F *hLambdaDCADaughterTracks = new TH1F("hLambdaDCADaughterTracks","DCA between #Lambda daughters;DCA (cm);Counts",200,0,10);
  TH1F *hLambdaDCAPosDaughPrimVertex = new TH1F("hLambdaDCAPosDaughPrimVertex","DCA of proton to PV;DCA (cm);Counts",500,0,100);
  TH1F *hLambdaDCANegDaughPrimVertex = new TH1F("hLambdaDCANegDaughPrimVertex","DCA of pion to PV;DCA (cm);Counts",500,0,100);
  TH1F *hLambdaTransverseRadius = new TH1F("hLambdaTransverseRadius","Transverse radius of #Lambda decay vertex;r_{xy} (cm);Counts;",400,0,200);
  TH1F *hLambdaCosPointingAngle = new TH1F("hLambdaCosPointingAngle","Cosine of #Lambda pointing angle;cos(PA);Counts",200,0.95,1);
  TH1F *hXiDCADaughterTracks = new TH1F("hXiDCADaughterTracks","DCA between #Xi daughters;DCA (cm);Counts",100,0,10);
  TH1F *hXiDCAV0DaughterTracks = new TH1F("hXiDCAV0DaughterTracks","DCA between #Lambda daughters, stemming from #Xi;DCA (cm);Counts",100,0,10);
  TH1F *hXiDCAV0PrimVertex = new TH1F("hXiDCAV0PrimVertex","DCA of #Lambda to PV, stemming from #Xi;DCA (cm);Counts",500,0,100);
  TH1F *hXiDCABachPrimVertex = new TH1F("hXiDCABachPrimVertex","DCA of bachelor pion to PV;DCA (cm);Counts",500,0,100);
  TH1F *hXiDCAPosDaughPrimVertex = new TH1F("hXiDCAPosDaughPrimVertex","DCA of daughter proton to PV, stemming from #Xi;DCA (cm);Counts",500,0,100);
  TH1F *hXiDCANegDaughPrimVertex = new TH1F("hXiDCANegDaughPrimVertex","DCA of daughter pion to PV, stemming from #Xi;DCA (cm);Counts",500,0,100);
  TH1F *hXiTransverseRadius = new TH1F("hXiTransverseRadius","Transverse radius of #Xi decay vertex;r_{xy} (cm);Counts",400,0,200);
  TH1F *hXiV0TransverseRadius = new TH1F("hXiV0TransverseRadius","Transverse rsdius of #Lambda decay vertex, stemming form #Xi;r_{xy} (cm);Counts",400,0,200);
  TH1F *hXiCosPointingAngle = new TH1F("hXiCosPointingAngle","Cosine of #Xi pointing angle;cos(PA);Counts",200,0.95,1);
  TH1F *hXiV0CosPointingAngle = new TH1F("hXiV0CosPointingAngle","Cosine of #Lambda pointing angle, stemming from #Xi;cos(PA);Counts",200,0.95,1);
  TH1F *hOmegaDCADaughterTracks = new TH1F("hOmegaDCADaughterTracks","DCA between #Omega daughters;DCA (cm);Counts",200,0,10);
  TH1F *hOmegaDCAV0DaughterTracks = new TH1F("hOmegaDCAV0DaughterTracks","DCA between #Lambda daughters, stemming from #Omega;DCA (cm);Counts",200,0,10);
  TH1F *hOmegaDCAV0PrimVertex = new TH1F("hOmegaDCAV0PrimVertex","DCA of #Lambda to PV, stemming from #Omega;DCA (cm);Counts",500,0,100);
  TH1F *hOmegaDCABachPrimVertex = new TH1F("hOmegaDCABachPrimVertex","DCA of bachelor kaon to PV;DCA (cm);Counts",500,0,100);
  TH1F *hOmegaDCAPosDaughPrimVertex = new TH1F("hOmegaDCAPosDaughPrimVertex","DCA of daughter proton to PV, stemming from #Omega;DCA (cm);Counts",500,0,100);
  TH1F *hOmegaDCANegDaughPrimVertex = new TH1F("hOmegaDCANegDaughPrimVertex","DCA of daughter pion to PV, stemming from #Omega;DCA (cm);Counts",500,0,100);
  TH1F *hOmegaTransverseRadius = new TH1F("hOmegaTransverseRadius","Transverse radius of #Omega decay vertex;r_{xy} (cm);Counts",400,0,200);
  TH1F *hOmegaV0TransverseRadius = new TH1F("hOmegaV0TransverseRadius","Transverse rsdius of #Lambda decay vertex, stemming form #Omega;r_{xy} (cm);Counts",400,0,200);
  TH1F *hOmegaCosPointingAngle = new TH1F("hOmegaCosPointingAngle","Cosine of #Omega pointing angle;cos(PA);Counts",200,0.95,1);
  TH1F *hOmegaV0CosPointingAngle = new TH1F("hOmegaV0CosPointingAngle","Cosine of #Lambda pointing angle, stemming from #Omega;cos(PA);Counts",200,0.95,1);

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
  TH1F *hProtonPt = new TH1F("hProtonPt","Transverse momentum of proton;p_{T} (GeV/c);Counts",200,0,10);
  TH2F *hProtonEtaPhi = new TH2F("hProtonEtaPhi","Pseudorapidity vs azimuthal angle of proton;#eta;#varphi (rad)",200,-1,+1,200,0,TMath::TwoPi());
  TH1F *hAntiProtonPt = new TH1F("hAntiProtonPt","Transverse momentum of Anti-proton;p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hLambdaPt = new TH1F("hLambdaPt","Transverse momentum of #Lambda;p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hLambdaPosDaughPt = new TH1F("hLambdaPosDaughPt","Transverse momentum of #Lambda daughter proton;p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hLambdaNegDaughPt = new TH1F("hLambdaNegDaughPt","Transverse momentum of #Lambda daughter pion;p_{T} (GeV/c);Counts",200,0,10);
  TH2F *hLambdaEtaPhi = new TH2F("hLambdaEtaPhi","Pseudorapidity vs azimuthal angle of #Lambda;#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH2F *hLambdaPosDaughEtaPhi = new TH2F("hLambdaPosDaughEtaPhi","Pseudorapidity vs azimuthal angle of #Lambda daughter proton;#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH2F *hLambdaNegDaughEtaPhi = new TH2F("hLambdaNegDaughEtaPhi","Pseudorapidity vs azimuthal angle of #Lambda daughter pion;#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH1F *hAntiLambdaPt = new TH1F("hAntiLambdaPt","Transverse momentum of #bar{#Lambda};p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hXimPt = new TH1F("hXimPt","Transverse momentum of #Xi^{-};p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hXimBachPt = new TH1F("hXimBachPt","Transverse momentum of bachelor pion;p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hXimPosDaughPt = new TH1F("hXimPosDaughPt","Transverse momentum of daughter proton from #Xi;p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hXimNegDaughPt = new TH1F("hXimNegDaughPt","Transverse momentum of daughter pion from #Xi;p_{T} (GeV/c);Counts",200,0,10);
  TH2F *hXimEtaPhi = new TH2F("hXimEtaPhi","Pseudorapidity vs azimuthal angle of #Xi^{-};#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH2F *hXimBachEtaPhi = new TH2F("hXimBachEtaPhi","Pseudorapidity vs azimuthal angle of #Xi^{-} bachelor pion;#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH2F *hXimPosDaughEtaPhi = new TH2F("hXimPosDaughEtaPhi","Pseudorapidity vs azimuthal angle of #Xi^{-} daughter proton;#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH2F *hXimNegDaughEtaPhi = new TH2F("hXimNegDaughEtaPhi","Pseudorapidity vs azimuthal angle of #Xi^{-} daughter pion;#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH1F *hXipPt = new TH1F("hXipPt","Transverse momentum of #Xi^{+};p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hOmegamPt = new TH1F("hOmegamPt","Transverse momentum of #Omega^{-};p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hOmegamBachPt = new TH1F("hOmegamBachPt","Transverse momentum of bachelor kaon;p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hOmegamPosDaughPt = new TH1F("hOmegamPosDaughPt","Transverse momentum of daughter proton from #Omega;p_{T} (GeV/c);Counts",200,0,10);
  TH1F *hOmegamNegDaughPt = new TH1F("hOmegamNegDaughPt","Transverse momentum of daughter pion from #Omega;p_{T} (GeV/c);Counts",200,0,10);
  TH2F *hOmegamEtaPhi = new TH2F("hOmegamEtaPhi","Pseudorapidity vs azimuthal angle of #Omega^{-};#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH2F *hOmegamBachEtaPhi = new TH2F("hOmegamBachEtaPhi","Pseudorapidity vs azimuthal angle of #Omega^{-} bachelor kaon;#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH2F *hOmegamPosDaughEtaPhi = new TH2F("hOmegamPosDaughEtaPhi","Pseudorapidity vs azimuthal angle of #Omega^{-} daughter proton;#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH2F *hOmegamNegDaughEtaPhi = new TH2F("hOmegamNegDaughEtaPhi","Pseudorapidity vs azimuthal angle of #Omega^{-} daughter pion;#eta;#varphi (rad)",200,-1,1,200,0,TMath::TwoPi());
  TH1F *hOmegapPt = new TH1F("hOmegapPt","Transverse momentum of #Omega^{+};p_{T} (GeV/c);Counts",200,0,10);

  fOutput->Add(hProtonPt);
  fOutput->Add(hProtonEtaPhi);
  fOutput->Add(hAntiProtonPt);
  fOutput->Add(hLambdaPt);
  fOutput->Add(hLambdaPosDaughPt);
  fOutput->Add(hLambdaNegDaughPt);
  fOutput->Add(hLambdaEtaPhi);
  fOutput->Add(hLambdaPosDaughEtaPhi);
  fOutput->Add(hLambdaNegDaughEtaPhi);
  fOutput->Add(hAntiLambdaPt);
  fOutput->Add(hXimPt);
  fOutput->Add(hXimBachPt);
  fOutput->Add(hXimPosDaughPt);
  fOutput->Add(hXimNegDaughPt);
  fOutput->Add(hXimEtaPhi);
  fOutput->Add(hXimBachEtaPhi);
  fOutput->Add(hXimPosDaughEtaPhi);
  fOutput->Add(hXimNegDaughEtaPhi);
  fOutput->Add(hXipPt);
  fOutput->Add(hOmegamPt);
  fOutput->Add(hOmegamBachPt);
  fOutput->Add(hOmegamPosDaughPt);
  fOutput->Add(hOmegamNegDaughPt);
  fOutput->Add(hOmegamEtaPhi);
  fOutput->Add(hOmegamBachEtaPhi);
  fOutput->Add(hOmegamPosDaughEtaPhi);
  fOutput->Add(hOmegamNegDaughEtaPhi);
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
  TH2F *hdEdxVsP = new TH2F("hdEdxVsP","TPC dE/dx vs momentum;p (GeV/c);TPC #frac{dE}{dx} (a.u.)",1000,0,10,1000,0,1000);
  TH2F *hBetaVsP = new TH2F("hBetaVsP","TOF #beta vs momentum;p (GeV/c);TOF #beta",1000,0,5,1000,0.1,1.1);
  TH2F *hProtonNsigmaTPC = new TH2F("hProtonNsigmaTPC","Nsigma TPC vs momentum;p (GeV/c);n_{#sigma,TPC}",1000,0,10,200,-10,10);
  TH2F *hProtonNsigmaTOF = new TH2F("hProtonNsigmaTOF","Nsigma TOF vs momentum;p (GeV/c);n_{#sigma,TOF}",1000,0,10,200,-10,10);
  TH2F *hProtonNsigmaTPCvsTOF = new TH2F("hProtonNsigmaTPCvsTOF","Nsigma TPC vs TOF;n_{#sigma,TPC};n_{#sigma,TOF};",200,-10,10,200,-10,10);
  TH2F *hProtonNsigmaTPCwPID = new TH2F("hProtonNsigmaTPCwPID","Nsimga TPC vs momentum after PID selection;p (GeV/c);n_{#sigma,TPC}",500,0,5,100,-5,5);
  TH2F *hProtonNsigmaTOFwPID = new TH2F("hProtonNsigmaTOFwPID","Nsimga TOF vs momentum after PID selection;p (GeV/c);n_{#sigma,TOF}",500,0,5,100,-5,5);
  TH2F *hProtonNsigmaTPCvsTOFwCuts = new TH2F("hProtonNsigmaTPCvsTOFwCuts","Nsigma TPC vs TOF with hypothesis cuts;n_{#sigma,TPC};n_{#sigma,TOF};",200,-10,10,200,-10,10);

  fOutput->Add(hdEdxVsP);
  fOutput->Add(hBetaVsP);
  fOutput->Add(hProtonNsigmaTPC);
  fOutput->Add(hProtonNsigmaTOF);
  fOutput->Add(hProtonNsigmaTPCvsTOF);
  fOutput->Add(hProtonNsigmaTPCwPID);
  fOutput->Add(hProtonNsigmaTOFwPID);
  fOutput->Add(hProtonNsigmaTPCvsTOFwCuts);

  // Define histograms related to pair analysis
  TH1F *hNPairStatistics = new TH1F("hNPairStatistics","Number of pairs under certain condition",10,0.5,10.5);
  hNPairStatistics->GetXaxis()->SetBinLabel(1,"p-#Lambda");
  hNPairStatistics->GetXaxis()->SetBinLabel(2,"#Lambda-#Lambda");
  hNPairStatistics->GetXaxis()->SetBinLabel(3,"p-#Xi^{-}");
  hNPairStatistics->GetXaxis()->SetBinLabel(4,"p-#Omega^{-}");
  hNPairStatistics->GetXaxis()->SetBinLabel(5,"#Lambda-#Xi^{-}");
  hNPairStatistics->GetXaxis()->SetBinLabel(6,"#Xi^{-}-#Omega^{-}");
  hNPairStatistics->GetXaxis()->SetBinLabel(7,"#Omega^{-}-#Omega^{-}");
  hNPairStatistics->GetYaxis()->SetTitle("Counts");
  TH1F *hNSharedTracksProtonLambda = new TH1F("hNSharedTracksProtonLambda","Number of p-#Lambda pairs which shares tracks in an event",10,0,10);
  TH1F *hNSharedTracksProtonXi = new TH1F("hNSharedTracksProtonXi","Number of p-#Xi pairs which shares tracks in an event",10,0,10);
  TH1F *hNSharedTracksProtonOmega = new TH1F("hNSharedTracksProtonOmega","Number of p-#Omega pairs which shares tracks in an event",10,0,10);
  TH1F *hNSharedTracksLambdaLambda = new TH1F("hNSharedTracksLambdaLambda","Number of #Lambda-#Lambda pairs which shares tracks in an event",10,0,10);
  TH1F *hNSharedTracksLambdaXi = new TH1F("hNSharedTracksLambdaXi","Number of #Lambda-#Xi pairs which shares tracks in an event",10,0,10);
  TH1F *hNSharedTracksLambdaOmega = new TH1F("hNSharedTracksLambdaOmega","Number of #Lambda-#Omega pairs which shares tracks in an event",10,0,10);
  TH1F *hNSharedTracksXiXi = new TH1F("hNSharedTracksXiXi","Number of #Xi-#Xi pairs which shares tracks in an event",10,0,10);
  TH1F *hNSharedTracksXiOmega = new TH1F("hNSharedTracksXiOmega","Number of #Xi-#Omega pairs which shares tracks in an event",10,0,10);
  TH1F *hNSharedTracksOmegaOmega = new TH1F("hNSharedTracksOmegaOmega","Number of #Omega-#Omega pairs which shares tracks in an event",10,0,10);

  fOutput->Add(hNPairStatistics);
  fOutput->Add(hNSharedTracksProtonLambda);
  fOutput->Add(hNSharedTracksProtonXi);
  fOutput->Add(hNSharedTracksProtonOmega);
  fOutput->Add(hNSharedTracksLambdaLambda);
  fOutput->Add(hNSharedTracksLambdaXi);
  fOutput->Add(hNSharedTracksLambdaOmega);
  fOutput->Add(hNSharedTracksXiXi);
  fOutput->Add(hNSharedTracksXiOmega);
  fOutput->Add(hNSharedTracksOmegaOmega);

  // Define histograms related to invariant mass for dibaryons
  TH1F *hInvMassProtonLambda = new TH1F("hInvMassProtonLambda","Invariant mass of p-#Lambda pair;M_{p#Lambda} (GeV/c^{2});Counts",1000,2,2.4);
  TH1F *hInvMassLambdaLambda = new TH1F("hInvMassLambdaLambda","Invariant mass of #Lambda-#Lambda pair;M_{#Lambda#Lambda} (GeV/c^{2});Counts",1000,2.2,2.3);
  TH1F *hInvMassProtonXi = new TH1F("hInvMassProtonXi","Invariant mass of p-#Xi pair;M_{p#Xi} (GeV/c^{2});Counts",1000,2.25,2.35);
  TH1F *hInvMassProtonOmega = new TH1F("hInvMassProtonOmega","Invariant mass of p-#Omega pair;M_{p#Omega} (GeV/c^{2});Counts",1000,2.6,2.7);
  TH1F *hInvMassLambdaXi = new TH1F("hInvMassLambdaXi","Invariant mass of #Lambda-#Xi pair;M_{#Lambda#Xi} (GeV/c^{2});Counts",1000,2.55,2.65);
  TH2F *hInvMassPtProtonLambda = new TH2F("hInvMassPtProtonLambda","Invariant mass vs transverse momentum of p-#Lambda pair;M_{p#Lambda} (GeV/c^{2});p_{T} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassPtLambdaLambda = new TH2F("hInvMassPtLambdaLambda","Invariant mass vs transverse momentum of #Lambda-#Lambda pair;M_{#Lambda#Lambda} (GeV/c^{2});p_{T} (GeV/c)",1000,2.2,3.2,100,0,10);
  TH2F *hInvMassPtProtonXi = new TH2F("hInvMassPtProtonXi","Invariant mass vs transverse momentum of p-#Xi pair;M_{p#Xi} (GeV/c^{2});p_{T} (GeV/c)",1000,2.2,3.2,100,0,10);
  TH2F *hInvMassPtProtonOmega = new TH2F("hInvMassPtProtonOmega","Invariant mass vs transverse momentum of p-#Omega pair;M_{p#Omega} (GeV/c^{2});p_{T} (GeV/c)",1000,2.6,3.6,100,0,10);
  TH2F *hInvMassPtLambdaXi = new TH2F("hInvMassPtLambdaXi","Invariant mass vs transverse momentum of #Lambda-#Xi pair;M_{#Lambda#Xi} (GeV/c^{2});p_{T} (GeV/c)",1000,2.4,3.4,100,0,10);
  TH2F *hInvMassPtXiOmega = new TH2F("hInvMassPtXiOmega","Invariant mass vs transverse momentum of #Xi-#Omega pair;M_{#Xi#Omega} (GeV/c^{2});p_{T} (GeV/c)",1000,2.9,3.9,100,0,10);
  TH2F *hInvMassPtOmegaOmega = new TH2F("hInvMassPtOmegaOmega","Invariant mass vs transverse momentum of #Omega-#Omega pair;M_{#Omega#Omega} (GeV/c^{2});p_{T} (GeV/c)",1000,3.3,4.3,100,0,10);
  TH2F *hInvMassRelKProtonLambda = new TH2F("hInvMassRelKProtonLambda","Invariant mass vs relative momentum of p-#Lambda pair;M_{p#Lambda} (GeV/c^{2});k^{*} (MeV/c)",1000,2,3,100,0,1000);
  TH2F *hInvMassRelKLambdaLambda = new TH2F("hInvMassRelKLambdaLambda","Invariant mass vs relative momentum of #Lambda-#Lambda pair;M_{#Lambda#Lambda} (GeV/c^{2});k^{*} (MeV/c)",1000,2.2,3.2,100,0,1000);
  TH2F *hInvMassRelKProtonXi = new TH2F("hInvMassRelKProtonXi","Invariant mass vs relative momentum of p-#Xi pair;M_{p#Xi} (GeV/c^{2});k^{*} (MeV/c)",1000,2.2,3.2,100,0,1000);
  TH2F *hInvMassRelKProtonOmega = new TH2F("hInvMassRelKProtonOmega","Invariant mass vs relative momentum of p-#Omega pair;M_{p#Omega} (GeV/c^{2});k^{*} (MeV/c)",1000,2.6,3.6,100,0,1000);
  TH2F *hInvMassRelKLambdaXi = new TH2F("hInvMassRelKLambdaXi","Invariant mass vs relative momentum of #Lambda-#Xi pair;M_{#Lambda#Xi} (GeV/c^{2});k^{*} (MeV/c)",1000,2.4,3.4,100,0,1000);
  TH2F *hInvMassRelKXiOmega = new TH2F("hInvMassRelKXiOmega","Invariant mass vs relative momentum of #Xi-#Omega pair;M_{#Xi#Omega} (GeV/c^{2});k^{*} (MeV/c)",1000,2.9,3.9,100,0,1000);
  TH2F *hInvMassRelKOmegaOmega = new TH2F("hInvMassRelKOmegaOmega","Invariant mass vs relative momentum of #Omega-#Omega pair;M_{#Omega#Omega} (GeV/c^{2});k^{*} (MeV/c)",1000,3.3,4.3,100,0,1000);
  TH2F *hInvMassPtProtonLambdaME = new TH2F("hInvMassPtProtonLambdaME","Invariant mass vs transverse momentum of p-#Lambda pair ME;M_{p#Lambda} (GeV/c^{2});p_{T} (GeV/c)",1000,2,3,100,0,10);
  TH2F *hInvMassPtLambdaLambdaME = new TH2F("hInvMassPtLambdaLambdaME","Invariant mass vs transverse momentum of #Lambda-#Lambda pair ME;M_{#Lambda#Lambda} (GeV/c^{2});p_{T} (GeV/c)",1000,2.2,3.2,100,0,10);
  TH2F *hInvMassPtProtonXiME = new TH2F("hInvMassPtProtonXiME","Invariant mass vs transverse momentum of p-#Xi pair ME;M_{p#Xi} (GeV/c^{2});p_{T} (GeV/c)",1000,2.2,3.2,100,0,10);
  TH2F *hInvMassPtProtonOmegaME = new TH2F("hInvMassPtProtonOmegaME","Invariant mass vs transverse momentum of p-#Omega pair ME;M_{p#Omega} (GeV/c^{2});p_{T} (GeV/c)",1000,2.6,3.6,100,0,10);
  TH2F *hInvMassPtLambdaXiME = new TH2F("hInvMassPtLambdaXiME","Invariant mass vs transverse momentum of #Lambda-#Xi pair ME;M_{#Lambda#Xi} (GeV/c^{2});p_{T} (GeV/c)",1000,2.4,3.4,100,0,10);
  TH2F *hInvMassPtXiOmegaME = new TH2F("hInvMassPtXiOmegaME","Invariant mass vs transverse momentum of #Xi-#Omega pair ME;M_{#Xi#Omega} (GeV/c^{2});p_{T} (GeV/c)",1000,2.9,3.9,100,0,10);
  TH2F *hInvMassPtOmegaOmegaME = new TH2F("hInvMassPtOmegaOmegaME","Invariant mass vs transverse momentum of #Omega-#Omega pair ME;M_{#Omega#Omega} (GeV/c^{2});p_{T} (GeV/c)",1000,3.3,4.3,100,0,10);
  TH2F *hInvMassRelKProtonLambdaME = new TH2F("hInvMassRelKProtonLambdaME","Invariant mass vs relative momentum of p-#Lambda pair ME;M_{p#Lambda} (GeV/c^{2});k^{*} (MeV/c)",1000,2,3,100,0,1000);
  TH2F *hInvMassRelKLambdaLambdaME = new TH2F("hInvMassRelKLambdaLambdaME","Invariant mass vs relative momentum of #Lambda-#Lambda pair ME;M_{#Lambda#Lambda} (GeV/c^{2});k^{*} (MeV/c)",1000,2.2,3.2,100,0,1000);
  TH2F *hInvMassRelKProtonXiME = new TH2F("hInvMassRelKProtonXiME","Invariant mass vs relative momentum of p-#Xi pair ME;M_{p#Xi} (GeV/c^{2});k^{*} (MeV/c)",1000,2.2,3.2,100,0,1000);
  TH2F *hInvMassRelKProtonOmegaME = new TH2F("hInvMassRelKProtonOmegaME","Invariant mass vs relative momentum of p-#Omega pair ME;M_{p#Omega} (GeV/c^{2});k^{*} (MeV/c)",1000,2.6,3.6,100,0,1000);
  TH2F *hInvMassRelKLambdaXiME = new TH2F("hInvMassRelKLambdaXiME","Invariant mass vs relative momentum of #Lambda-#Xi pair ME;M_{#Lambda#Xi} (GeV/c^{2});k^{*} (MeV/c)",1000,2.4,3.4,100,0,1000);
  TH2F *hInvMassRelKXiOmegaME = new TH2F("hInvMassRelKXiOmegaME","Invariant mass vs relative momentum of #Xi-#Omega pair ME;M_{#Xi#Omega} (GeV/c^{2});k^{*} (MeV/c)",1000,2.9,3.9,100,0,1000);
  TH2F *hInvMassRelKOmegaOmegaME = new TH2F("hInvMassRelKOmegaOmegaME","Invariant mass vs relative momentum of #Omega-#Omega pair ME;M_{#Omega#Omega} (GeV/c^{2});k^{*} (MeV/c)",1000,3.3,4.3,100,0,1000);

  fOutput->Add(hInvMassProtonLambda);
  fOutput->Add(hInvMassLambdaLambda);
  fOutput->Add(hInvMassProtonXi);
  fOutput->Add(hInvMassProtonOmega);
  fOutput->Add(hInvMassLambdaXi);
  fOutput->Add(hInvMassPtProtonLambda);
  fOutput->Add(hInvMassPtLambdaLambda);
  fOutput->Add(hInvMassPtProtonXi);
  fOutput->Add(hInvMassPtProtonOmega);
  fOutput->Add(hInvMassPtLambdaXi);
  fOutput->Add(hInvMassPtXiOmega);
  fOutput->Add(hInvMassPtOmegaOmega);
  fOutput->Add(hInvMassRelKProtonLambda);
  fOutput->Add(hInvMassRelKLambdaLambda);
  fOutput->Add(hInvMassRelKProtonXi);
  fOutput->Add(hInvMassRelKProtonOmega);
  fOutput->Add(hInvMassRelKLambdaXi);
  fOutput->Add(hInvMassRelKXiOmega);
  fOutput->Add(hInvMassRelKOmegaOmega);
  fOutput->Add(hInvMassPtProtonLambdaME);
  fOutput->Add(hInvMassPtLambdaLambdaME);
  fOutput->Add(hInvMassPtProtonXiME);
  fOutput->Add(hInvMassPtProtonOmegaME);
  fOutput->Add(hInvMassPtLambdaXiME);
  fOutput->Add(hInvMassPtXiOmegaME);
  fOutput->Add(hInvMassPtOmegaOmegaME);
  fOutput->Add(hInvMassRelKProtonLambdaME);
  fOutput->Add(hInvMassRelKLambdaLambdaME);
  fOutput->Add(hInvMassRelKProtonXiME);
  fOutput->Add(hInvMassRelKProtonOmegaME);
  fOutput->Add(hInvMassRelKLambdaXiME);
  fOutput->Add(hInvMassRelKXiOmegaME);
  fOutput->Add(hInvMassRelKOmegaOmegaME);

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

  // Reference multiplicity
  Int_t refMult08 = 0;
  if(fAnalysisType == "ESD") {
    refMult08 = AliESDtrackCuts::GetReferenceMultiplicity(esdEvent,AliESDtrackCuts::kTrackletsITSTPC,0.8,0);
  } else if(fAnalysisType == "AOD") {
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(aodEvent->GetHeader());
    if(!header) {
      AliWarning("Could not get AOD header... return!");
      return;
    }
    refMult08 = header->GetRefMultiplicityComb08();
  }
  dynamic_cast<TH1F*>(fOutput->FindObject("hRefMultiplicity08"))->Fill(refMult08);

  // Find multiplicity bin and z vertex bin for mixed event
  Int_t multBin = -1;
  for(Int_t i=0; i<25; i++) {
    if(4*i < refMult08 && refMult08 <= 4*(i+1)) {
      multBin = i;
      break;
    }
  }
  if(refMult08 > 100) multBin = 25;
  dynamic_cast<TH1F*>(fOutput->FindObject("hNBinsMultMixing"))->Fill(multBin);

  Double_t zVertex = primaryVtxPos[2];
  Int_t zBin = -1;
  for(Int_t i=0; i<10; i++) {
    if((-10+i*2) < zVertex && zVertex < (-10+(i+1)*2)) {
      zBin = i;
      break;
    }
  }
  dynamic_cast<TH1F*>(fOutput->FindObject("hNBinsVertexMixing"))->Fill(zBin);

  // Magnetic field
  Double_t bz = -10.;
  if     (fAnalysisType == "ESD") bz = esdEvent->GetMagneticField();
  else if(fAnalysisType == "AOD") bz = aodEvent->GetMagneticField();

  // PDG mass
  const Double_t massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass(); 
  const Double_t massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass(); 
  const Double_t massXi     = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  const Double_t massOmega  = TDatabasePDG::Instance()->GetParticle(3334)->Mass();

  // Speed of light
  const Double_t c = 2.99792457999999984e-02; // cm/ps

  //______________________________________________________________________________
  // Loop over the reconstructed candidates

  // - Track loop for proton and anti-proton selecton

  Int_t nTrack = 0;
  if     (fAnalysisType == "ESD") nTrack = esdEvent->GetNumberOfTracks();
  else if(fAnalysisType == "AOD") nTrack = aodEvent->GetNumberOfTracks();

  fProtonArray->Clear("C");
  Int_t countProton = 0;

  if(fAnalysisType == "AOD") { // set global track

    // Reset global track reference
    for(Int_t i=0; i < fTrackBuffSize; i++) fTrackArray[i] = 0;

    // Store global track reference
    for(Int_t iTrack=0; iTrack < nTrack; iTrack++) {

      AliAODTrack *track = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
      if(!track) {
        AliFatal("Not a standard AOD");
        continue;
      }

      const Int_t trackID = track->GetID();
      if(trackID < 0) continue;
      if(trackID >= fTrackBuffSize) {
        AliWarning(Form("Track ID too big for buffer: ID: %d, buffer %d",trackID,fTrackBuffSize));
        continue;
      }

      if (fTrackArray[trackID]) {
        if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) continue;
      }

      fTrackArray[trackID] = track; // store the pointer to the global track
    }
  }

  for(Int_t iTrack=0; iTrack < nTrack; iTrack++) {

    AliESDtrack *esdTrack = 0x0;
    AliAODTrack *aodTrack = 0x0;

    // Initialisation of the local variables
    Double_t charge = 0.;
    Double_t pt     = 0.;
    Double_t eta    = 0.;
    Double_t phi    = 0.;
    Double_t p      = 0.;
    Double_t dEdx   = 0.;
    Double_t beta   = 0.;

    Float_t nTPCCrossedRows  = -1.;
    UShort_t nTPCClusters    = -1;
    UShort_t nTPCSharedCls   = -1;
    UShort_t nTPCFindableCls = -1;
    Float_t ratio = -1.;

    Bool_t isthereTPC = kFALSE;
    Bool_t isthereTOF = kFALSE;

    Float_t nSigmaTPCproton      = -10.;
    Float_t nSigmaTPCkaon        = -10.;
    Float_t nSigmaTPCpion        = -10.;
    Float_t nSigmaTPCelectron    = -10.;
    Float_t nSigmaTOFproton      = -10.;
    Float_t nSigmaTOFkaon        = -10.;
    Float_t nSigmaTOFpion        = -10.;
    Float_t nSigmaTOFelectron    = -10.;
    Float_t nSigmaTPCTOFproton   = -10.;
    Float_t nSigmaTPCTOFkaon     = -10.;
    Float_t nSigmaTPCTOFpion     = -10.;
    Float_t nSigmaTPCTOFelectron = -10.;

    Float_t DCAxy = -999., DCAz = -999.;

    if(fAnalysisType == "ESD") {

      esdTrack = (AliESDtrack*)esdEvent->GetTrack(iTrack);
      if(!esdTrack) {
        AliFatal("Not a standard ESD");
        continue;
      }

      charge = esdTrack->Charge();
      pt     = esdTrack->Pt();
      eta    = esdTrack->Eta();
      phi    = esdTrack->Phi();
      dEdx   = esdTrack->GetTPCsignal();
      p      = esdTrack->GetTPCmomentum();

      const Double_t len = esdTrack->GetIntegratedLength();
      const Double_t tim = esdTrack->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->GetTPCmomentum());
      if(tim != 0.) beta = len / (tim * c);

      nTPCCrossedRows = esdTrack->GetTPCCrossedRows();
      nTPCClusters    = esdTrack->GetTPCNcls();
      nTPCSharedCls   = esdTrack->GetTPCnclsS();
      nTPCFindableCls = esdTrack->GetTPCNclsF();
      if(nTPCFindableCls > 0) ratio = nTPCCrossedRows / nTPCFindableCls;

      AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,esdTrack);
      AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,esdTrack);
      if(statusTPC == AliPIDResponse::kDetPidOk) isthereTPC = kTRUE;
      if(statusTOF == AliPIDResponse::kDetPidOk) isthereTOF = kTRUE;

      nSigmaTPCproton      = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
      nSigmaTPCkaon        = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon);
      nSigmaTPCpion        = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion);
      nSigmaTPCelectron    = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kElectron);
      nSigmaTOFproton      = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton);
      nSigmaTOFkaon        = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon);
      nSigmaTOFpion        = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion);
      nSigmaTOFelectron    = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kElectron);
      nSigmaTPCTOFproton   = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));
      nSigmaTPCTOFkaon     = TMath::Sqrt(pow(nSigmaTPCkaon,2.) + pow(nSigmaTOFkaon,2.));
      nSigmaTPCTOFpion     = TMath::Sqrt(pow(nSigmaTPCpion,2.) + pow(nSigmaTOFpion,2.));
      nSigmaTPCTOFelectron = TMath::Sqrt(pow(nSigmaTPCelectron,2.) + pow(nSigmaTOFelectron,2.));

      esdTrack->GetImpactParameters(DCAxy,DCAz);

    } // end of ESD treatment

    else if(fAnalysisType == "AOD") {

      aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
      if(!aodTrack) {
        AliFatal("Not a standard AOD");
        continue;
      }

      if(!((AliAODTrack*)aodTrack)->TestFilterBit(fFilterBit)) continue;

      const Int_t trackID = aodTrack->GetID();
      AliAODTrack *globalTrack;
      if(trackID < 0) {
        if(!fTrackArray[-trackID-1]) continue; // check if a global track exists for the TPC-only track
        globalTrack = fTrackArray[-trackID-1];
      } else {
        globalTrack = aodTrack;
      }

      charge = aodTrack->Charge();
      pt     = aodTrack->Pt();
      eta    = aodTrack->Eta();
      phi    = aodTrack->Phi();
      p      = globalTrack->GetTPCmomentum();
      dEdx   = globalTrack->GetTPCsignal();

      const Double_t len = globalTrack->GetIntegratedLength();
      const Double_t tim = globalTrack->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(globalTrack->GetTPCmomentum());
      if(tim != 0.) beta = len / (tim * c);

      nTPCCrossedRows = aodTrack->GetTPCCrossedRows();
      nTPCClusters    = aodTrack->GetTPCNcls();
      nTPCSharedCls   = aodTrack->GetTPCnclsS();
      nTPCFindableCls = aodTrack->GetTPCNclsF();
      if(nTPCFindableCls > 0) ratio = nTPCCrossedRows / nTPCFindableCls;

      AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,globalTrack);
      AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,globalTrack);
      if(statusTPC == AliPIDResponse::kDetPidOk) isthereTPC = kTRUE;
      if(statusTOF == AliPIDResponse::kDetPidOk) isthereTOF = kTRUE;

      nSigmaTPCproton      = fPIDResponse->NumberOfSigmasTPC(globalTrack,AliPID::kProton);
      nSigmaTPCkaon        = fPIDResponse->NumberOfSigmasTPC(globalTrack,AliPID::kKaon);
      nSigmaTPCpion        = fPIDResponse->NumberOfSigmasTPC(globalTrack,AliPID::kPion);
      nSigmaTPCelectron    = fPIDResponse->NumberOfSigmasTPC(globalTrack,AliPID::kElectron);
      nSigmaTOFproton      = fPIDResponse->NumberOfSigmasTOF(globalTrack,AliPID::kProton);
      nSigmaTOFkaon        = fPIDResponse->NumberOfSigmasTOF(globalTrack,AliPID::kKaon);
      nSigmaTOFpion        = fPIDResponse->NumberOfSigmasTOF(globalTrack,AliPID::kPion);
      nSigmaTOFelectron    = fPIDResponse->NumberOfSigmasTOF(globalTrack,AliPID::kElectron);
      nSigmaTPCTOFproton   = TMath::Sqrt(pow(nSigmaTPCproton,2.) + pow(nSigmaTOFproton,2.));
      nSigmaTPCTOFkaon     = TMath::Sqrt(pow(nSigmaTPCkaon,2.) + pow(nSigmaTOFkaon,2.));
      nSigmaTPCTOFpion     = TMath::Sqrt(pow(nSigmaTPCpion,2.) + pow(nSigmaTOFpion,2.));
      nSigmaTPCTOFelectron = TMath::Sqrt(pow(nSigmaTPCelectron,2.) + pow(nSigmaTOFelectron,2.));

      globalTrack->GetImpactParameters(DCAxy,DCAz);

    } // end of AOD treatment

    // Check track quality
    dynamic_cast<TH1F*>(fOutput->FindObject("hNCrossedRows"))        ->Fill(nTPCCrossedRows);
    dynamic_cast<TH1F*>(fOutput->FindObject("hNCluster"))            ->Fill(nTPCClusters);
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedCluster"))      ->Fill(nTPCSharedCls);
    dynamic_cast<TH1F*>(fOutput->FindObject("hNFindableCluster"))    ->Fill(nTPCFindableCls);
    dynamic_cast<TH1F*>(fOutput->FindObject("hRatioFindableCrossed"))->Fill(ratio);

    if(nTPCClusters    < 80)   continue;
    if(nTPCCrossedRows < 70)   continue;
    if(ratio           < 0.83) continue;
    if(nTPCSharedCls   > 0)    continue;

    // Fill histos related to PID
    if(isthereTPC) {
      dynamic_cast<TH2F*>(fOutput->FindObject("hdEdxVsP"))->Fill(p,dEdx);
      if(charge > 0.) dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaTPC"))->Fill(p,nSigmaTPCproton);
    }
    if(isthereTOF) {
      dynamic_cast<TH2F*>(fOutput->FindObject("hBetaVsP"))->Fill(p,beta);
      if(charge > 0.) dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaTOF"))->Fill(p,nSigmaTOFproton);
    }

    // proton PID
    if(p < 0.75) { // for p < 0.75 use TPC only
      if(!isthereTPC) continue;
      if(TMath::Abs(nSigmaTPCproton) > 3.) continue;
    }
    else if(p > 0.75) { // for p > 0.75 use TPC & TOF
      if(!isthereTPC || !isthereTOF) continue;
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaTPCvsTOF"))->Fill(nSigmaTPCproton,nSigmaTOFproton);

      // Reject if the other hypotheis fits better
      if(nSigmaTPCTOFproton > nSigmaTPCTOFkaon) continue;
      if(nSigmaTPCTOFproton > nSigmaTPCTOFpion) continue;
      if(nSigmaTPCTOFproton > nSigmaTPCTOFelectron) continue;

      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaTPCvsTOFwCuts"))->Fill(nSigmaTPCproton,nSigmaTOFproton);
      if(nSigmaTPCTOFproton > 3.) continue;
    }

    if(charge > 0.) {
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonDCAxyDCAz"))->Fill(DCAxy,DCAz);
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaTPCwPID"))->Fill(p,nSigmaTPCproton);
      if(isthereTOF) dynamic_cast<TH2F*>(fOutput->FindObject("hProtonNsigmaTOFwPID"))->Fill(p,nSigmaTOFproton);
    }

    // proton selection
    if(TMath::Abs(eta) > 0.8) continue;
    if(pt < 0.5 || 4.05 < pt) continue;
    if(TMath::Abs(DCAz)  > 0.2) continue; 
    if(TMath::Abs(DCAxy) > 0.1) continue;

    // Store candidates after selection
    if(charge > 0.) { // for proton

      dynamic_cast<TH1F*>(fOutput->FindObject("hProtonPt")) ->Fill(pt);
      dynamic_cast<TH2F*>(fOutput->FindObject("hProtonEtaPhi"))->Fill(eta,phi);

      dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(1);

      if(fAnalysisType == "ESD") {

        AliESDtrack *track = (AliESDtrack*)fProtonArray->ConstructedAt(countProton);
        esdTrack->Copy(*track);
        countProton++;
      }
      else if(fAnalysisType == "AOD") {

        AliAODTrack *track = (AliAODTrack*)fProtonArray->ConstructedAt(countProton);
        *track = *aodTrack;
        countProton++;
      }
    }
    else if(charge < 0.) { // for anti-proton

      dynamic_cast<TH1F*>(fOutput->FindObject("hAntiProtonPt"))->Fill(pt);

      dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(2);
    }

  }// end of track loop


  // - V0 loop for Lambda selection

  Int_t nV0 = 0;
  if     (fAnalysisType == "ESD") nV0 = esdEvent->GetNumberOfV0s();
  else if(fAnalysisType == "AOD") nV0 = aodEvent->GetNumberOfV0s();

  fLambdaArray->Clear("C");
  Int_t countLambda = 0;
  
  for(Int_t iV0=0; iV0<nV0; iV0++) {

    AliESDv0 *esdV0 = 0x0;
    AliAODv0 *aodV0 = 0x0;

    // Initialisation of the local variables
    Double_t invMassLambda     = 0.;
    Double_t invMassAntiLambda = 0.;
    Double_t invMassK0S        = 0.;

    Double_t dcaV0Dghters    = -1.; 
    Double_t cosPointAngle   = -1.;
    Double_t vtxPosV0[3]     = {-999., -999., -999.};
    Double_t radius          = -999.;
    Double_t dcaPosToPrimVtx = -1.;
    Double_t dcaNegToPrimVtx = -1.;

    Bool_t isPosProton = kFALSE;
    Bool_t isPosPion   = kFALSE;
    Bool_t isNegProton = kFALSE;
    Bool_t isNegPion   = kFALSE;

    Double_t ptV0   = 0.;
    Double_t ptPos  = 0.;
    Double_t ptNeg  = 0.;
    Double_t etaV0  = -20.;
    Double_t etaPos = -20.;
    Double_t etaNeg = -20.;
    Double_t phiV0  = 720.;
    Double_t phiPos = 720.;
    Double_t phiNeg = 720.;

    if(fAnalysisType == "ESD") {
      
      esdV0 = (AliESDv0*)esdEvent->GetV0(iV0);
      if(!esdV0) continue;

      if(esdV0->GetOnFlyStatus()) continue; // select offline V0

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

      // Check TPC PID status
      AliPIDResponse::EDetPidStatus statusPosTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,pTrack);
      AliPIDResponse::EDetPidStatus statusNegTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,nTrack);
      if(statusPosTPC != AliPIDResponse::kDetPidOk) continue;
      if(statusNegTPC != AliPIDResponse::kDetPidOk) continue;

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
      esdV0->ChangeMassHypothesis(310);
      invMassK0S        = esdV0->GetEffMass();

      // Get topological values
      dcaV0Dghters    = TMath::Abs(esdV0->GetDcaV0Daughters());
      dcaPosToPrimVtx = TMath::Abs(pTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));
      dcaNegToPrimVtx = TMath::Abs(pTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));
      cosPointAngle   = esdV0->GetV0CosineOfPointingAngle(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]);
      esdV0->GetXYZ(vtxPosV0[0], vtxPosV0[1], vtxPosV0[2]);
      radius = TMath::Sqrt(vtxPosV0[0]*vtxPosV0[0] + vtxPosV0[1]*vtxPosV0[1]);

      // Get kinematic values
      ptV0   = esdV0->Pt();
      ptPos  = pTrack->Pt();
      ptNeg  = nTrack->Pt();
      etaV0  = esdV0->Eta();
      etaPos = pTrack->Eta();
      etaNeg = nTrack->Eta();
      phiV0  = esdV0->Phi();
      phiPos = pTrack->Phi();
      phiNeg = nTrack->Phi();

    } // end of ESD treatment

    else if(fAnalysisType == "AOD") {

      aodV0 = (AliAODv0*)aodEvent->GetV0(iV0);
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

      // Out-of-bunch pile-up removal: require either a hit in the ITS SPD or ITS SDD, or TOF in-bunch timing
      if(fPileupCut) {
        if(!(pTrack->HasPointOnITSLayer(0)) && !(pTrack->HasPointOnITSLayer(1)) &&
           !(pTrack->HasPointOnITSLayer(4)) && !(pTrack->HasPointOnITSLayer(5)) &&
           !(pTrack->GetTOFBunchCrossing() == 0)) continue;
        if(!(nTrack->HasPointOnITSLayer(0)) && !(nTrack->HasPointOnITSLayer(1)) &&
           !(nTrack->HasPointOnITSLayer(4)) && !(nTrack->HasPointOnITSLayer(5)) &&
           !(nTrack->GetTOFBunchCrossing() == 0)) continue;
      }

      // Check TPC PID status
      AliPIDResponse::EDetPidStatus statusPosTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,pTrack);
      AliPIDResponse::EDetPidStatus statusNegTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,nTrack);
      if(statusPosTPC != AliPIDResponse::kDetPidOk) continue;
      if(statusNegTPC != AliPIDResponse::kDetPidOk) continue;

      // Daughter track PID using TPC
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton)) < 5.0) isPosProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion  )) < 5.0) isPosPion   = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton)) < 5.0) isNegProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion  )) < 5.0) isNegPion   = kTRUE;

      // Calculate the invariant mass
      invMassLambda     = aodV0->MassLambda();
      invMassAntiLambda = aodV0->MassAntiLambda();
      invMassK0S        = aodV0->MassK0Short();

      // Get topological values
      dcaV0Dghters    = TMath::Abs(aodV0->DcaV0Daughters());
      dcaPosToPrimVtx = TMath::Abs(aodV0->DcaPosToPrimVertex());
      dcaNegToPrimVtx = TMath::Abs(aodV0->DcaNegToPrimVertex());
      cosPointAngle   = aodV0->CosPointingAngle(primaryVtxPos);
      aodV0->GetXYZ(vtxPosV0);
      radius = TMath::Sqrt(vtxPosV0[0]*vtxPosV0[0] + vtxPosV0[1]*vtxPosV0[1]);

      // Get kinematic values
      ptV0   = aodV0->Pt();
      ptPos  = pTrack->Pt();
      ptNeg  = nTrack->Pt();
      etaV0  = aodV0->Eta();
      etaPos = pTrack->Eta();
      etaNeg = nTrack->Eta();
      phiV0  = aodV0->Phi();
      phiPos = pTrack->Phi();
      phiNeg = nTrack->Phi();

    } // end of AOD treatment

    if(isPosProton && isNegPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassLambdawoCuts"))->Fill(invMassLambda);
    }
    else if(isNegProton && isPosPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassAntiLambdawoCuts"))->Fill(invMassAntiLambda);
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassMissIDK0S"))->Fill(invMassK0S);
    if(0.48 < invMassK0S && invMassK0S < 0.515) continue; // reject K0Short

    // Lambda selection
    if(TMath::Abs(etaPos) > 0.8) continue;
    if(TMath::Abs(etaNeg) > 0.8) continue;
    if(ptV0 < 0.3) continue;
    if(TMath::Abs(vtxPosV0[0]) > 100) continue;
    if(TMath::Abs(vtxPosV0[1]) > 100) continue;
    if(TMath::Abs(vtxPosV0[2]) > 100) continue;
    if(dcaV0Dghters > 1.5) continue;
    if(cosPointAngle < 0.99) continue;
    if(dcaPosToPrimVtx < 0.05) continue;
    if(dcaNegToPrimVtx < 0.05) continue;
    if(radius < 0.2 || 100 < radius) continue;

    // Store candidates after selection
    if(isPosProton && isNegPion) { // for Lambda

      dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassLambda"))->Fill(invMassLambda,ptV0);

      if(TMath::Abs(invMassLambda - massLambda) < 0.004) { // mass window selection

        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaDCADaughterTracks"))    ->Fill(dcaV0Dghters);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaDCAPosDaughPrimVertex"))->Fill(dcaPosToPrimVtx);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaDCANegDaughPrimVertex"))->Fill(dcaNegToPrimVtx);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaCosPointingAngle"))     ->Fill(cosPointAngle);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaTransverseRadius"))     ->Fill(radius);

        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaPt"))            ->Fill(ptV0);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaPosDaughPt"))    ->Fill(ptPos);
        dynamic_cast<TH1F*>(fOutput->FindObject("hLambdaNegDaughPt"))    ->Fill(ptNeg);
        dynamic_cast<TH2F*>(fOutput->FindObject("hLambdaEtaPhi"))        ->Fill(etaV0,phiV0);
        dynamic_cast<TH2F*>(fOutput->FindObject("hLambdaPosDaughEtaPhi"))->Fill(etaPos,phiPos);
        dynamic_cast<TH2F*>(fOutput->FindObject("hLambdaNegDaughEtaPhi"))->Fill(etaNeg,phiNeg);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics")) ->Fill(3);

        if(fAnalysisType == "ESD") {

          AliESDv0 *v0 = (AliESDv0*)fLambdaArray->ConstructedAt(countLambda);
          esdV0->Copy(*v0);
          countLambda++;
        }
        else if (fAnalysisType == "AOD") {

          AliAODv0 *v0 = (AliAODv0*)fLambdaArray->ConstructedAt(countLambda);
          *v0 = *aodV0;
          countLambda++;
        }
      }
    }
    if(isNegProton && isPosPion) { // for Anti-Lambda

      dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassAntiLambda"))->Fill(invMassAntiLambda,ptV0);

      if(TMath::Abs(invMassAntiLambda - massLambda) < 0.004) { // mass window selection

        dynamic_cast<TH1F*>(fOutput->FindObject("hAntiLambdaPt"))->Fill(ptV0);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(4);
      }
    }

  }// end of V0 loop


  // - Cascade loop for Xi and Omega selection

  Int_t nCascade = 0;
  if     (fAnalysisType == "ESD") nCascade = esdEvent->GetNumberOfCascades();
  else if(fAnalysisType == "AOD") nCascade = aodEvent->GetNumberOfCascades();

  fXiArray->Clear("C");
  fOmegaArray->Clear("C");
  Int_t countXi = 0;
  Int_t countOmega = 0;

  for(Int_t iXi = 0; iXi < nCascade; iXi++) {

    AliESDcascade *esdXi = 0x0;
    AliAODcascade *aodXi = 0x0;

    // Initialisation of the local variables
    Double_t invMassXiMinus    = 0.;
    Double_t invMassXiPlus     = 0.;
    Double_t invMassOmegaMinus = 0.;
    Double_t invMassOmegaPlus  = 0.;
    Double_t invMassLambdaAsCascDghter = 0.;

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

    Bool_t isBachelorKaon = kFALSE;
    Bool_t isBachelorPion = kFALSE;
    Bool_t isPosProton    = kFALSE;
    Bool_t isPosPion      = kFALSE;
    Bool_t isNegProton    = kFALSE;
    Bool_t isNegPion      = kFALSE;

    Double_t ptXi    = 0.;
    Double_t ptBach  = 0.;
    Double_t ptPos   = 0.;
    Double_t ptNeg   = 0.;
    Double_t etaXi   = -20.;
    Double_t etaBach = -20.;
    Double_t etaPos  = -20.;
    Double_t etaNeg  = -20.;
    Double_t phiXi   = 720.;
    Double_t phiBach = 720.;
    Double_t phiPos  = 720.;
    Double_t phiNeg  = 720.;

    Short_t chargeXi = -2;

    Bool_t standerdXi    = kTRUE;
    Bool_t standerdOmega = kTRUE;

    if(fAnalysisType == "ESD") {
      
      esdXi = (AliESDcascade*)esdEvent->GetCascade(iXi);
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

      // Check TPC PID status
      AliPIDResponse::EDetPidStatus statusBachTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,bachTrack);
      AliPIDResponse::EDetPidStatus statusPosTPC  = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,pTrack);
      AliPIDResponse::EDetPidStatus statusNegTPC  = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,nTrack);
      if(statusBachTPC != AliPIDResponse::kDetPidOk) continue;
      if(statusPosTPC  != AliPIDResponse::kDetPidOk) continue;
      if(statusNegTPC  != AliPIDResponse::kDetPidOk) continue;

      // Bachelor and daughter track PID using TPC
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bachTrack,AliPID::kKaon)) < 4.0) isBachelorKaon = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bachTrack,AliPID::kPion)) < 4.0) isBachelorPion = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton )) < 4.0) isPosProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion   )) < 4.0) isPosPion   = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton )) < 4.0) isNegProton = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion   )) < 4.0) isNegPion   = kTRUE;

      // Calculate the invariant mass
      invMassLambdaAsCascDghter = esdXi->GetEffMass();
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
      dcaXiDghters     = TMath::Abs(esdXi->GetDcaXiDaughters());
      cosPointAngleXi  = esdXi->GetCascadeCosineOfPointingAngle(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]);
      radiusXi         = TMath::Sqrt(vtxPosXi[0]*vtxPosXi[0] + vtxPosXi[1]*vtxPosXi[1]);
      dcaV0Dghters     = TMath::Abs(esdXi->GetDcaV0Daughters());
      dcaV0ToPrimVtx   = TMath::Abs(esdXi->GetD(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]));
      cosPointAngleV0  = esdXi->GetV0CosineOfPointingAngle(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]);
      radiusV0         = TMath::Sqrt(vtxPosV0[0]*vtxPosV0[0] + vtxPosV0[1]*vtxPosV0[1]);
      dcaBachToPrimVtx = TMath::Abs(bachTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));
      dcaPosToPrimVtx  = TMath::Abs(pTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));
      dcaNegToPrimVtx  = TMath::Abs(nTrack->GetD(primaryVtxPos[0], primaryVtxPos[1], bz));

      // Get kinematic values
      ptXi    = esdXi->Pt();
      ptBach  = bachTrack->Pt();
      ptPos   = pTrack->Pt();
      ptNeg   = nTrack->Pt();
      etaXi   = esdXi->Eta();
      etaBach = bachTrack->Eta();
      etaPos  = pTrack->Eta();
      etaNeg  = nTrack->Eta();
      phiXi   = esdXi->Phi();
      phiBach = bachTrack->Phi();
      phiPos  = pTrack->Phi();
      phiNeg  = nTrack->Phi();

    } // end of ESD treatment
    
    else if(fAnalysisType == "AOD") {

      aodXi = (AliAODcascade*)aodEvent->GetCascade(iXi);
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

      // Check TPC PID status
      AliPIDResponse::EDetPidStatus statusBachTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,bachTrack);
      AliPIDResponse::EDetPidStatus statusPosTPC  = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,pTrack);
      AliPIDResponse::EDetPidStatus statusNegTPC  = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,nTrack);
      if(statusBachTPC != AliPIDResponse::kDetPidOk) continue;
      if(statusPosTPC  != AliPIDResponse::kDetPidOk) continue;
      if(statusNegTPC  != AliPIDResponse::kDetPidOk) continue;

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
      dcaXiDghters     = TMath::Abs(aodXi->DcaXiDaughters());
      cosPointAngleXi  = aodXi->CosPointingAngleXi(primaryVtxPos[0], primaryVtxPos[1], primaryVtxPos[2]);
      radiusXi         = TMath::Sqrt(vtxPosXi[0]*vtxPosXi[0] + vtxPosXi[1]*vtxPosXi[1]);
      dcaV0Dghters     = TMath::Abs(aodXi->DcaV0Daughters());
      dcaV0ToPrimVtx   = TMath::Abs(aodXi->DcaV0ToPrimVertex());
      cosPointAngleV0  = aodXi->CosPointingAngle(primaryVtxPos);
      radiusV0         = TMath::Sqrt(vtxPosV0[0]*vtxPosV0[0] + vtxPosV0[1]*vtxPosV0[1]);
      dcaBachToPrimVtx = TMath::Abs(aodXi->DcaBachToPrimVertex());
      dcaPosToPrimVtx  = TMath::Abs(aodXi->DcaPosToPrimVertex());
      dcaNegToPrimVtx  = TMath::Abs(aodXi->DcaNegToPrimVertex());

      // Get kinematic values
      ptXi    = aodXi->Pt();
      ptBach  = bachTrack->Pt();
      ptPos   = pTrack->Pt();
      ptNeg   = nTrack->Pt();
      etaXi   = aodXi->Eta();
      etaBach = bachTrack->Eta();
      etaPos  = pTrack->Eta();
      etaNeg  = nTrack->Eta();
      phiXi   = aodXi->Phi();
      phiBach = bachTrack->Phi();
      phiPos  = pTrack->Phi();
      phiNeg  = nTrack->Phi();

    } // end of AOD treatment

    if((chargeXi<0) && isBachelorPion && isPosProton && isNegPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassXimwoCuts"))   ->Fill(invMassXiMinus);
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassMissIDOmegam"))->Fill(invMassOmegaMinus);
    }
    if((chargeXi>0) && isBachelorPion && isNegProton && isPosPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassXipwoCuts"))   ->Fill(invMassXiPlus);
    }
    if((chargeXi<0) && isBachelorKaon && isPosProton && isNegPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassOmegamwoCuts"))->Fill(invMassOmegaMinus);
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassMissIDXim"))   ->Fill(invMassXiMinus);
    }
    if((chargeXi>0) && isBachelorKaon && isNegProton && isPosPion) {
      dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassOmegapwoCuts"))->Fill(invMassOmegaPlus);
    }

    // Kinematic cuts for Xi and Omega
    if(TMath::Abs(etaBach) > 0.8) continue;
    if(TMath::Abs(etaPos)  > 0.8) continue;
    if(TMath::Abs(etaNeg) > 0.8) continue;
    if(ptBach < 0.3) continue;
    if(ptPos  < 0.3) continue;
    if(ptNeg  < 0.3) continue;

    // Topological cuts for Xi
    if(dcaXiDghters > 1.6)      standerdXi = kFALSE;
    if(dcaV0Dghters > 1.6)      standerdXi = kFALSE;
    if(dcaV0ToPrimVtx < 0.07)   standerdXi = kFALSE;
    if(dcaBachToPrimVtx < 0.05) standerdXi = kFALSE;
    if(dcaPosToPrimVtx < 0.05)  standerdXi = kFALSE;
    if(dcaNegToPrimVtx < 0.05)  standerdXi = kFALSE;
    if(cosPointAngleXi < 0.98)  standerdXi = kFALSE;
    if(cosPointAngleV0 < 0.97)  standerdXi = kFALSE;
    if(radiusXi < 0.8 || 200 < radiusXi) standerdXi = kFALSE;
    if(radiusV0 < 1.4 || 200 < radiusV0) standerdXi = kFALSE;

    // Topological cuts for Omega
    if(dcaXiDghters > 0.8)      standerdOmega = kFALSE;
    if(dcaV0Dghters > 1.2)      standerdOmega = kFALSE;
    if(dcaV0ToPrimVtx < 0.06)   standerdOmega = kFALSE;
    if(dcaBachToPrimVtx < 0.05) standerdOmega = kFALSE;
    if(dcaPosToPrimVtx < 0.05)  standerdOmega = kFALSE;
    if(dcaNegToPrimVtx < 0.05)  standerdOmega = kFALSE;
    if(cosPointAngleXi < 0.995) standerdOmega = kFALSE;
    if(cosPointAngleV0 < 0.97)  standerdOmega = kFALSE;
    if(radiusXi < 0.2 || 200 < radiusXi) standerdOmega = kFALSE;
    if(radiusV0 < 1.0 || 200 < radiusV0) standerdOmega = kFALSE;

    dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassLambdaAsCascDghter"))->Fill(invMassLambdaAsCascDghter);

    // Mass window cut for V0
    if(TMath::Abs(invMassLambdaAsCascDghter - massLambda) > 0.006) {
      standerdXi    = kFALSE;
      standerdOmega = kFALSE;
    }

    // Store candidates after selection
    if(standerdXi) {

      if((chargeXi<0) && // for Xi-
         (isBachelorPion && isPosProton && isNegPion) && // TPC PID info
         (invMassOmegaMinus < 1.667 || 1.677 < invMassOmegaMinus)) { // reject Omega-

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassXim"))->Fill(invMassXiMinus,ptXi);

        if(TMath::Abs(invMassXiMinus - massXi) < 0.005) { // mass window selection

          dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCADaughterTracks"))    ->Fill(dcaXiDghters);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCAV0DaughterTracks"))  ->Fill(dcaV0Dghters);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCAV0PrimVertex"))      ->Fill(dcaV0ToPrimVtx);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCABachPrimVertex"))    ->Fill(dcaBachToPrimVtx);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCAPosDaughPrimVertex"))->Fill(dcaPosToPrimVtx);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXiDCANegDaughPrimVertex"))->Fill(dcaNegToPrimVtx);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXiCosPointingAngle"))     ->Fill(cosPointAngleXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXiV0CosPointingAngle"))   ->Fill(cosPointAngleV0);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXiTransverseRadius"))     ->Fill(radiusXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXiV0TransverseRadius"))   ->Fill(radiusV0);

          dynamic_cast<TH1F*>(fOutput->FindObject("hXimPt"))            ->Fill(ptXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXimBachPt"))        ->Fill(ptBach);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXimPosDaughPt"))    ->Fill(ptPos);
          dynamic_cast<TH1F*>(fOutput->FindObject("hXimNegDaughPt"))    ->Fill(ptNeg);
          dynamic_cast<TH2F*>(fOutput->FindObject("hXimEtaPhi"))        ->Fill(etaXi,phiXi);
          dynamic_cast<TH2F*>(fOutput->FindObject("hXimBachEtaPhi"))    ->Fill(etaBach,phiBach);
          dynamic_cast<TH2F*>(fOutput->FindObject("hXimPosDaughEtaPhi"))->Fill(etaPos,phiPos);
          dynamic_cast<TH2F*>(fOutput->FindObject("hXimNegDaughEtaPhi"))->Fill(etaNeg,phiNeg);

          dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(5);

          if(fAnalysisType == "ESD") {

            AliESDcascade *xi = (AliESDcascade*)fXiArray->ConstructedAt(countXi);
            esdXi->Copy(*xi);
            countXi++;
          }
          else if (fAnalysisType == "AOD") {

            AliAODcascade *xi = (AliAODcascade*)fXiArray->ConstructedAt(countXi);
            *xi = *aodXi;
            countXi++;
          }
        }
      }
      if((chargeXi>0) && // for Xi+
         (isBachelorPion && isNegProton && isPosPion) && // TOC PID info
         (invMassOmegaPlus < 1.667 || 1.677 < invMassOmegaPlus)) { // reject Omega+

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassXip"))->Fill(invMassXiPlus,ptXi);

        if(TMath::Abs(invMassXiPlus - massXi) < 0.005) { // mass window selection

          dynamic_cast<TH1F*>(fOutput->FindObject("hXipPt"))->Fill(ptXi);

          dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(6);
        }
      }
    }

    if(standerdOmega) {

      if((chargeXi<0) && // for Omega-
         (isBachelorKaon && isPosProton && isNegPion) && // TPC PID info
         (invMassXiMinus < 1.317 || 1.327 < invMassXiMinus)) { // reject Xi-

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassOmegam"))->Fill(invMassOmegaMinus,ptXi);

        if(TMath::Abs(invMassOmegaMinus - massOmega) < 0.005) { // mass window selection

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

          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamPt"))            ->Fill(ptXi);
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamBachPt"))        ->Fill(ptBach);
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamPosDaughPt"))    ->Fill(ptPos);
          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegamNegDaughPt"))    ->Fill(ptNeg);
          dynamic_cast<TH2F*>(fOutput->FindObject("hOmegamEtaPhi"))        ->Fill(etaXi,phiXi);
          dynamic_cast<TH2F*>(fOutput->FindObject("hOmegamBachEtaPhi"))    ->Fill(etaBach,phiBach);
          dynamic_cast<TH2F*>(fOutput->FindObject("hOmegamPosDaughEtaPhi"))->Fill(etaPos,phiPos);
          dynamic_cast<TH2F*>(fOutput->FindObject("hOmegamNegDaughEtaPhi"))->Fill(etaNeg,phiNeg);

          dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(7);

          if(fAnalysisType == "ESD") {

            AliESDcascade *xi = (AliESDcascade*)fOmegaArray->ConstructedAt(countOmega);
            esdXi->Copy(*xi);
            countOmega++;
          }
          else if (fAnalysisType == "AOD") {

            AliAODcascade *xi = (AliAODcascade*)fOmegaArray->ConstructedAt(countOmega);
            *xi = *aodXi;
            countOmega++;
          }
        }
      }
      if((chargeXi>0) && // for Omega+
         (isBachelorKaon && isNegProton && isPosPion) && // TPC PID info
         (invMassXiPlus < 1.317 || 1.327 < invMassXiPlus)) { // reject Xi+

        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassOmegap"))->Fill(invMassOmegaPlus,ptXi);

        if(TMath::Abs(invMassOmegaPlus - massOmega) < 0.005) { // mass window selection

          dynamic_cast<TH1F*>(fOutput->FindObject("hOmegapPt"))->Fill(ptXi);

          dynamic_cast<TH1F*>(fOutput->FindObject("hNPartStatistics"))->Fill(8);
        }
      }
    }

  } // end of cascade loop

  PairCleaner(); // check candidates for shared daughter tracks


  //______________________________________________________________________________
  // Calculate invariant mass for dibaryons 

  const Int_t nProton = fProtonArray->GetEntriesFast();
  const Int_t nLambda = fLambdaArray->GetEntriesFast();
  const Int_t nXi     = fXiArray->GetEntriesFast();
  const Int_t nOmega  = fOmegaArray->GetEntriesFast();

  if(fAnalysisType == "ESD") {

    // ppK- -> proton + Lambda
    for(Int_t i=0; i<nProton; i++) {

      AliESDtrack *track = (AliESDtrack*)fProtonArray->At(i);
      if(!track) continue;

      for(Int_t j=0; j<nLambda; j++) {

        AliESDv0 *v0 = (AliESDv0*)fLambdaArray->At(j);
        if(!v0) continue;

        if(track->GetID() == v0->GetPindex()) continue;
        if(track->GetID() == v0->GetNindex()) continue;

        TLorentzVector trackProton, trackLambda, trackSum;

        trackProton.SetXYZM(track->Px(),track->Py(),track->Pz(),massProton);
        trackLambda.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),massLambda);
        trackSum = trackProton + trackLambda;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackProton,trackLambda);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(1);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassProtonLambda"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtProtonLambda"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKProtonLambda"))->Fill(mass,relK);
      }
    }

    // H-Dibaryon -> Lambda + Lambda
    for(Int_t i=0; i<nLambda; i++) {

      AliESDv0 *v01 = (AliESDv0*)fLambdaArray->At(i);
      if(!v01) continue;

      for(Int_t j=i+1; j<nLambda; j++) {

        AliESDv0 *v02 = (AliESDv0*)fLambdaArray->At(j);
        if(!v02) continue;

        if(v01->GetPindex() == v02->GetPindex()) continue;
        if(v01->GetNindex() == v02->GetNindex()) continue;

        TLorentzVector trackLambda1, trackLambda2, trackSum;

        trackLambda1.SetXYZM(v01->Px(),v01->Py(),v01->Pz(),massLambda);
        trackLambda2.SetXYZM(v02->Px(),v02->Py(),v02->Pz(),massLambda);
        trackSum = trackLambda1 + trackLambda2;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackLambda1,trackLambda2);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(2);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassLambdaLambda"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtLambdaLambda"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKLambdaLambda"))->Fill(mass,relK);
      }
    }

    // H-Dibaryon -> proton + Xi-
    for(Int_t i=0; i<nProton; i++) {

      AliESDtrack *track = (AliESDtrack*)fProtonArray->At(i);
      if(!track) continue;

      for(Int_t j=0; j<nXi; j++) {

        AliESDcascade *xi = (AliESDcascade*)fXiArray->At(j);
        if(!xi) continue;

        if(track->GetID() == xi->GetBindex()) continue;
        if(track->GetID() == xi->GetPindex()) continue;
        if(track->GetID() == xi->GetNindex()) continue;

        TLorentzVector trackProton, trackXi, trackSum;

        trackProton.SetXYZM(track->Px(),track->Py(),track->Pz(),massProton);
        trackXi.SetXYZM(xi->Px(),xi->Py(),xi->Pz(),massXi);
        trackSum = trackProton + trackXi;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackProton,trackXi);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(3);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassProtonXi"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtProtonXi"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKProtonXi"))->Fill(mass,relK);
      }
    }

    // pOmega -> proton + Omega-
    for(Int_t i=0; i<nProton; i++) {

      AliESDtrack *track = (AliESDtrack*)fProtonArray->At(i);
      if(!track) continue;

      for(Int_t j=0; j<nOmega; j++) {

        AliESDcascade *xi = (AliESDcascade*)fOmegaArray->At(j);
        if(!xi) continue;

        if(track->GetID() == xi->GetBindex()) continue;
        if(track->GetID() == xi->GetPindex()) continue;
        if(track->GetID() == xi->GetNindex()) continue;

        TLorentzVector trackProton, trackOmega, trackSum;

        trackProton.SetXYZM(track->Px(),track->Py(),track->Pz(),massProton);
        trackOmega.SetXYZM(xi->Px(),xi->Py(),xi->Pz(),massOmega);
        trackSum = trackProton + trackOmega;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackProton, trackOmega);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(4);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassProtonOmega"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtProtonOmega"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKProtonOmega"))->Fill(mass,relK);
      }
    }

    // nOmega- -> Lambda + Xi-
    for(Int_t i=0; i<nLambda; i++) {

      AliESDv0 *v0 = (AliESDv0*)fLambdaArray->At(i);
      if(!v0) continue;

      for(Int_t j=0; j<nXi; j++) {

        AliESDcascade *xi = (AliESDcascade*)fXiArray->At(j);
        if(!xi) continue;

        if(v0->GetPindex() == xi->GetPindex()) continue;
        if(v0->GetNindex() == xi->GetNindex()) continue;
        if(v0->GetNindex() == xi->GetBindex()) continue;

        TLorentzVector trackLambda, trackXi, trackSum;

        trackLambda.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),massLambda);
        trackXi.SetXYZM(xi->Px(),xi->Py(),xi->Pz(),massXi);
        trackSum = trackLambda + trackXi;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackLambda,trackXi);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(5);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassLambdaXi"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtLambdaXi"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKLambdaXi"))->Fill(mass,relK);
      }
    }

    // Di-Omega -> Xi- + Omega-
    for(Int_t i=0; i<nXi; i++) {

      AliESDcascade *xi1 = (AliESDcascade*)fXiArray->At(i);
      if(!xi1) continue;

      for(Int_t j=0; j<nOmega; j++) {

        AliESDcascade *xi2 = (AliESDcascade*)fOmegaArray->At(j);
        if(!xi2) continue;

        if(xi1->GetPindex() == xi2->GetPindex()) continue;
        if(xi1->GetNindex() == xi2->GetNindex()) continue;
        if(xi1->GetNindex() == xi2->GetBindex()) continue;
        if(xi1->GetBindex() == xi2->GetNindex()) continue;

        TLorentzVector trackXi, trackOmega, trackSum;

        trackXi.SetXYZM(xi1->Px(),xi1->Py(),xi1->Pz(),massXi);
        trackOmega.SetXYZM(xi2->Px(),xi2->Py(),xi2->Pz(),massOmega);
        trackSum = trackXi + trackOmega;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackXi,trackOmega);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(6);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtXiOmega"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKXiOmega"))->Fill(mass,relK);
      }
    }

    // Di-Omega -> Omega- + Omega-
    for(Int_t i=0; i<nOmega; i++) {

      AliESDcascade *xi1 = (AliESDcascade*)fOmegaArray->At(i);
      if(!xi1) continue;

      for(Int_t j=i+1; j<nOmega; j++) {

        AliESDcascade *xi2 = (AliESDcascade*)fOmegaArray->At(j);
        if(!xi2) continue;

        if(xi1->GetPindex() == xi2->GetPindex()) continue;
        if(xi1->GetNindex() == xi2->GetNindex()) continue;
        if(xi1->GetNindex() == xi2->GetBindex()) continue;
        if(xi1->GetBindex() == xi2->GetNindex()) continue;

        TLorentzVector trackOmega1, trackOmega2, trackSum;

        trackOmega1.SetXYZM(xi1->Px(),xi1->Py(),xi1->Pz(),massOmega);
        trackOmega2.SetXYZM(xi2->Px(),xi2->Py(),xi2->Pz(),massOmega);
        trackSum = trackOmega1 + trackOmega2;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackOmega1,trackOmega2);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(7);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtOmegaOmega"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKOmegaOmega"))->Fill(mass,relK);
      }
    }

  }

  else if(fAnalysisType == "AOD") {

    // ppK- -> proton + Lambda
    for(Int_t i=0; i<nProton; i++) {

      AliAODTrack *track = (AliAODTrack*)fProtonArray->At(i);
      if(!track) continue;

      Int_t trackID = track->GetID();
      if(trackID < 0) trackID = -trackID - 1;

      for(Int_t j=0; j<nLambda; j++) {

        AliAODv0 *v0 = (AliAODv0*)fLambdaArray->At(j);
        if(!v0) continue;

        TLorentzVector trackProton, trackLambda, trackSum;

        trackProton.SetXYZM(track->Px(),track->Py(),track->Pz(),massProton);
        trackLambda.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),massLambda);
        trackSum = trackProton + trackLambda;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackProton, trackLambda);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(1);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassProtonLambda"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtProtonLambda"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKProtonLambda"))->Fill(mass,relK);
      }
    }

    // H-Dibaryon -> Lambda + Lambda
    for(Int_t i=0; i<nLambda; i++) {

      AliAODv0 *v01 = (AliAODv0*)fLambdaArray->At(i);
      if(!v01) continue;

      for(Int_t j=i+1; j<nLambda; j++) {

        AliAODv0 *v02 = (AliAODv0*)fLambdaArray->At(j);
        if(!v02) continue;

        TLorentzVector trackLambda1, trackLambda2, trackSum;

        trackLambda1.SetXYZM(v01->Px(),v01->Py(),v01->Pz(),massLambda);
        trackLambda2.SetXYZM(v02->Px(),v02->Py(),v02->Pz(),massLambda);
        trackSum = trackLambda1 + trackLambda2;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackLambda1,trackLambda2);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(2);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassLambdaLambda"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtLambdaLambda"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKLambdaLambda"))->Fill(mass,relK);
      }
    }

    // H-Dibaryon -> proton + Xi-
    for(Int_t i=0; i<nProton; i++) {

      AliAODTrack *track = (AliAODTrack*)fProtonArray->At(i);
      if(!track) continue;

      Int_t trackID = track->GetID();
      if(trackID < 0) trackID = -trackID - 1;

      for(Int_t j=0; j<nXi; j++) {

        AliAODcascade *xi = (AliAODcascade*)fXiArray->At(j);
        if(!xi) continue;

        TLorentzVector trackProton, trackXi, trackSum;

        trackProton.SetXYZM(track->Px(),track->Py(),track->Pz(),massProton);
        trackXi.SetXYZM(xi->Px(),xi->Py(),xi->Pz(),massXi);
        trackSum = trackProton + trackXi;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackProton,trackXi);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(3);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassProtonXi"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtProtonXi"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKProtonXi"))->Fill(mass,relK);
      }
    }

    // pOmega -> proton + Omega-
    for(Int_t i=0; i<nProton; i++) {

      AliAODTrack *track = (AliAODTrack*)fProtonArray->At(i);
      if(!track) continue;

      Int_t trackID = track->GetID();
      if(trackID < 0) trackID = -trackID - 1;

      for(Int_t j=0; j<nOmega; j++) {

        AliAODcascade *xi = (AliAODcascade*)fOmegaArray->At(j);
        if(!xi) continue;

        TLorentzVector trackProton, trackOmega, trackSum;

        trackProton.SetXYZM(track->Px(),track->Py(),track->Pz(),massProton);
        trackOmega.SetXYZM(xi->Px(),xi->Py(),xi->Pz(),massOmega);
        trackSum = trackProton + trackOmega;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackProton,trackOmega);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(4);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassProtonOmega"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtProtonOmega"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKProtonOmega"))->Fill(mass,relK);
      }
    }

    // nOmega- -> Lambda + Xi-
    for(Int_t i=0; i<nLambda; i++) {

      AliAODv0 *v0 = (AliAODv0*)fLambdaArray->At(i);
      if(!v0) continue;

      for(Int_t j=0; j<nXi; j++) {

        AliAODcascade *xi = (AliAODcascade*)fXiArray->At(j);
        if(!xi) continue;

        TLorentzVector trackLambda, trackXi, trackSum;

        trackLambda.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),massLambda);
        trackXi.SetXYZM(xi->Px(),xi->Py(),xi->Pz(),massXi);
        trackSum = trackLambda + trackXi;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackLambda,trackXi);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(5);
        dynamic_cast<TH1F*>(fOutput->FindObject("hInvMassLambdaXi"))->Fill(mass);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtLambdaXi"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKLambdaXi"))->Fill(mass,relK);
      }
    }

    // Di-Omega -> Xi- + Omega-
    for(Int_t i=0; i<nXi; i++) {

      AliAODcascade *xi1 = (AliAODcascade*)fXiArray->At(i);
      if(!xi1) continue;

      for(Int_t j=0; j<nOmega; j++) {

        AliAODcascade *xi2 = (AliAODcascade*)fOmegaArray->At(j);
        if(!xi2) continue;

        TLorentzVector trackXi, trackOmega, trackSum;

        trackXi.SetXYZM(xi1->Px(),xi1->Py(),xi1->Pz(),massXi);
        trackOmega.SetXYZM(xi2->Px(),xi2->Py(),xi2->Pz(),massOmega);
        trackSum = trackXi + trackOmega;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackXi,trackOmega);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(6);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtXiOmega"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKXiOmega"))->Fill(mass,relK);
      }
    }

    // Di-Omega -> Omega- + Omega-
    for(Int_t i=0; i<nOmega; i++) {

      AliAODcascade *xi1 = (AliAODcascade*)fOmegaArray->At(i);
      if(!xi1) continue;

      for(Int_t j=i+1; j<nOmega; j++) {

        AliAODcascade *xi2 = (AliAODcascade*)fOmegaArray->At(j);
        if(!xi2) continue;

        TLorentzVector trackOmega1, trackOmega2, trackSum;

        trackOmega1.SetXYZM(xi1->Px(),xi1->Py(),xi1->Pz(),massOmega);
        trackOmega2.SetXYZM(xi2->Px(),xi2->Py(),xi2->Pz(),massOmega);
        trackSum = trackOmega1 + trackOmega2;

        Double_t mass = trackSum.M();
        Double_t pt   = trackSum.Pt();
        Double_t relK = 1000*relKcalc(trackOmega1,trackOmega2);

        dynamic_cast<TH1F*>(fOutput->FindObject("hNPairStatistics"))->Fill(7);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtOmegaOmega"))->Fill(mass,pt);
        dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKOmegaOmega"))->Fill(mass,relK);
      }
    }

  }

  //______________________________________________________________________________
  // Event Mixing

  if((zBin >= 0) && (multBin >= 0) && fEventMixing) {

    EventPool& evPoolProton = fProtonEMpool[zBin][multBin];
    EventPool& evPoolLambda = fLambdaEMpool[zBin][multBin];
    EventPool& evPoolXi     = fXiEMpool[zBin][multBin];
    EventPool& evPoolOmega  = fOmegaEMpool[zBin][multBin];

    if(nProton > 0) {
      evPoolProton.push_front((TClonesArray*)fProtonArray->Clone());
      if(evPoolProton.size() > 10) evPoolProton.pop_back();
    }
    if(nLambda > 0) {
      evPoolLambda.push_front((TClonesArray*)fLambdaArray->Clone());
      if(evPoolLambda.size() > 10) evPoolLambda.pop_back();
    }
    if(nXi > 0) {
      evPoolXi.push_front((TClonesArray*)fXiArray->Clone());
      if(evPoolXi.size() > 10) evPoolXi.pop_back();
    }
    if(nOmega > 0) {
      evPoolOmega.push_front((TClonesArray*)fOmegaArray->Clone());
      if(evPoolOmega.size() > 10) evPoolOmega.pop_back();
    }

    if(fAnalysisType == "AOD") {

      // ppK- -> proton + Lambda
      for(Int_t i=0; i<nProton; i++) {

        AliAODTrack *track = (AliAODTrack*)fProtonArray->At(i);
        if(!track) continue;

        Int_t trackID = track->GetID();
        if(trackID < 0) trackID = -trackID - 1;

        for(UInt_t evnum=0; evnum<evPoolLambda.size(); evnum++) {

          if(evnum == 0 && nLambda > 0) continue; // same event

          for(Int_t j=0; j<(evPoolLambda.at(evnum))->GetEntriesFast(); j++) {

            AliAODv0 *v0 = (AliAODv0*)(evPoolLambda.at(evnum))->At(j);
            if(!v0) continue;

            TLorentzVector trackProton, trackLambda, trackSum;

            trackProton.SetXYZM(track->Px(),track->Py(),track->Pz(),massProton);
            trackLambda.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),massLambda);
            trackSum = trackProton + trackLambda;

            Double_t mass = trackSum.M();
            Double_t pt   = trackSum.Pt();
            Double_t relK = 1000*relKcalc(trackProton, trackLambda);

            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtProtonLambdaME"))->Fill(mass,pt);
            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKProtonLambdaME"))->Fill(mass,relK);
          }
        }
      }

      // H-Dibaryon -> Lambda + Lambda
      for(Int_t i=0; i<nLambda; i++) {

        AliAODv0 *v01 = (AliAODv0*)fLambdaArray->At(i);
        if(!v01) continue;

        for(UInt_t evnum=0; evnum<evPoolLambda.size(); evnum++) {

          if(evnum == 0) continue; // same event

          for(Int_t j=0; j<(evPoolLambda.at(evnum))->GetEntriesFast(); j++) {

            AliAODv0 *v02 = (AliAODv0*)(evPoolLambda.at(evnum))->At(j);
            if(!v02) continue;

            TLorentzVector trackLambda1, trackLambda2, trackSum;

            trackLambda1.SetXYZM(v01->Px(),v01->Py(),v01->Pz(),massLambda);
            trackLambda2.SetXYZM(v02->Px(),v02->Py(),v02->Pz(),massLambda);
            trackSum = trackLambda1 + trackLambda2;

            Double_t mass = trackSum.M();
            Double_t pt   = trackSum.Pt();
            Double_t relK = 1000*relKcalc(trackLambda1,trackLambda2);

            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtLambdaLambdaME"))->Fill(mass,pt);
            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKLambdaLambdaME"))->Fill(mass,relK);
          }
        }
      }

      // H-Dibaryon -> proton + Xi-
      for(Int_t i=0; i<nProton; i++) {

        AliAODTrack *track = (AliAODTrack*)fProtonArray->At(i);
        if(!track) continue;

        Int_t trackID = track->GetID();
        if(trackID < 0) trackID = -trackID - 1;

        for(UInt_t evnum=0; evnum<evPoolXi.size(); evnum++) {

          if(evnum == 0 && nXi > 0) continue; // same event

          for(Int_t j=0; j<(evPoolXi.at(evnum))->GetEntriesFast(); j++) {

            AliAODcascade *xi = (AliAODcascade*)(evPoolXi.at(evnum))->At(j);
            if(!xi) continue;

            TLorentzVector trackProton, trackXi, trackSum;

            trackProton.SetXYZM(track->Px(),track->Py(),track->Pz(),massProton);
            trackXi.SetXYZM(xi->Px(),xi->Py(),xi->Pz(),massXi);
            trackSum = trackProton + trackXi;

            Double_t mass = trackSum.M();
            Double_t pt   = trackSum.Pt();
            Double_t relK = 1000*relKcalc(trackProton,trackXi);

            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtProtonXiME"))->Fill(mass,pt);
            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKProtonXiME"))->Fill(mass,relK);
          }
        }
      }

      // pOmega -> proton + Omega-
      for(Int_t i=0; i<nProton; i++) {

        AliAODTrack *track = (AliAODTrack*)fProtonArray->At(i);
        if(!track) continue;

        Int_t trackID = track->GetID();
        if(trackID < 0) trackID = -trackID - 1;

        for(UInt_t evnum=0; evnum<evPoolOmega.size(); evnum++) {

          if(evnum == 0 && nOmega > 0) continue; // same event

          for(Int_t j=0; j<(evPoolOmega.at(evnum))->GetEntriesFast(); j++) {

            AliAODcascade *xi = (AliAODcascade*)(evPoolOmega.at(evnum))->At(j);
            if(!xi) continue;

            TLorentzVector trackProton, trackOmega, trackSum;

            trackProton.SetXYZM(track->Px(),track->Py(),track->Pz(),massProton);
            trackOmega.SetXYZM(xi->Px(),xi->Py(),xi->Pz(),massOmega);
            trackSum = trackProton + trackOmega;

            Double_t mass = trackSum.M();
            Double_t pt   = trackSum.Pt();
            Double_t relK = 1000*relKcalc(trackProton,trackOmega);

            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtProtonOmegaME"))->Fill(mass,pt);
            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKProtonOmegaME"))->Fill(mass,relK);
          }
        }
      }

      // nOmega- -> Lambda + Xi-
      for(Int_t i=0; i<nLambda; i++) {

        AliAODv0 *v0 = (AliAODv0*)fLambdaArray->At(i);
        if(!v0) continue;

        for(UInt_t evnum=0; evnum<evPoolXi.size(); evnum++) {

          if(evnum == 0 && nXi > 0) continue; // same event

          for(Int_t j=0; j<(evPoolXi.at(evnum))->GetEntriesFast(); j++) {

            AliAODcascade *xi = (AliAODcascade*)(evPoolXi.at(evnum))->At(j);
            if(!xi) continue;

            TLorentzVector trackLambda, trackXi, trackSum;

            trackLambda.SetXYZM(v0->Px(),v0->Py(),v0->Pz(),massLambda);
            trackXi.SetXYZM(xi->Px(),xi->Py(),xi->Pz(),massXi);
            trackSum = trackLambda + trackXi;

            Double_t mass = trackSum.M();
            Double_t pt   = trackSum.Pt();
            Double_t relK = 1000*relKcalc(trackLambda,trackXi);

            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtLambdaXiME"))->Fill(mass,pt);
            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKLambdaXiME"))->Fill(mass,relK);
          }
        }
      }

      // Di-Omega -> Xi- + Omega-
      for(Int_t i=0; i<nXi; i++) {

        AliAODcascade *xi1 = (AliAODcascade*)fXiArray->At(i);
        if(!xi1) continue;

        for(UInt_t evnum=0; evnum<evPoolOmega.size(); evnum++) {

          if(evnum == 0 && nOmega > 0) continue; // same event

          for(Int_t j=0; j<(evPoolOmega.at(evnum))->GetEntriesFast(); j++) {

            AliAODcascade *xi2 = (AliAODcascade*)(evPoolOmega.at(evnum))->At(j);
            if(!xi2) continue;

            TLorentzVector trackXi, trackOmega, trackSum;

            trackXi.SetXYZM(xi1->Px(),xi1->Py(),xi1->Pz(),massXi);
            trackOmega.SetXYZM(xi2->Px(),xi2->Py(),xi2->Pz(),massOmega);
            trackSum = trackXi + trackOmega;

            Double_t mass = trackSum.M();
            Double_t pt   = trackSum.Pt();
            Double_t relK = 1000*relKcalc(trackXi,trackOmega);

            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtXiOmegaME"))->Fill(mass,pt);
            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKXiOmegaME"))->Fill(mass,relK);
          }
        }
      }

      // Di-Omega -> Omega- + Omega-
      for(Int_t i=0; i<nOmega; i++) {

        AliAODcascade *xi1 = (AliAODcascade*)fOmegaArray->At(i);
        if(!xi1) continue;

        for(UInt_t evnum=0; evnum<evPoolOmega.size(); evnum++) {

          if(evnum == 0) continue; // same event

          for(Int_t j=0; j<(evPoolOmega.at(evnum))->GetEntriesFast(); j++) {

            AliAODcascade *xi2 = (AliAODcascade*)(evPoolOmega.at(evnum))->At(j);
            if(!xi2) continue;

            TLorentzVector trackOmega1, trackOmega2, trackSum;

            trackOmega1.SetXYZM(xi1->Px(),xi1->Py(),xi1->Pz(),massOmega);
            trackOmega2.SetXYZM(xi2->Px(),xi2->Py(),xi2->Pz(),massOmega);
            trackSum = trackOmega1 + trackOmega2;

            Double_t mass = trackSum.M();
            Double_t pt   = trackSum.Pt();
            Double_t relK = 1000*relKcalc(trackOmega1,trackOmega2);

            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassPtOmegaOmegaME"))->Fill(mass,pt);
            dynamic_cast<TH2F*>(fOutput->FindObject("hInvMassRelKOmegaOmegaME"))->Fill(mass,relK);
          }
        }
      }

    }
  }

}
//_______________________________________________________________________________________________
void AliAnalysisTaskDibaryons::PairCleaner()
{
  const Int_t nProton = fProtonArray->GetEntriesFast();
  const Int_t nLambda = fLambdaArray->GetEntriesFast();
  const Int_t nXi     = fXiArray->GetEntriesFast();
  const Int_t nOmega  = fOmegaArray->GetEntriesFast();

  if(fAnalysisType == "AOD") {

    Int_t nShared = 0;

    // p-Lambda
    for(Int_t i=0; i<nProton; i++) {

      AliAODTrack *track = (AliAODTrack*)fProtonArray->At(i);
      if(!track) continue;

      Int_t trackID = track->GetID();
      if(trackID < 0) trackID = -trackID - 1;

      for(Int_t j=0; j<nLambda; j++) {

        AliAODv0 *v0 = (AliAODv0*)fLambdaArray->At(j);
        if(!v0) continue;

        if(trackID == v0->GetPosID() || trackID == v0->GetNegID()) {

          fLambdaArray->RemoveAt(j);
          nShared++;
        }
      }
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedTracksProtonLambda"))->Fill(nShared);

    nShared = 0;

    // p-Xi
    for(Int_t i=0; i<nProton; i++) {

      AliAODTrack *track = (AliAODTrack*)fProtonArray->At(i);
      if(!track) continue;

      Int_t trackID = track->GetID();
      if(trackID < 0) trackID = -trackID - 1;

      for(Int_t j=0; j<nXi; j++) {

        AliAODcascade *xi = (AliAODcascade*)fXiArray->At(j);
        if(!xi) continue;

        if(trackID == xi->GetBachID() || trackID == xi->GetPosID() || trackID == xi->GetNegID()) {

          fXiArray->RemoveAt(j);
          nShared++;
        }
      }
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedTracksProtonXi"))->Fill(nShared);

    nShared = 0;

    // p-Omega
    for(Int_t i=0; i<nProton; i++) {

      AliAODTrack *track = (AliAODTrack*)fProtonArray->At(i);
      if(!track) continue;

      Int_t trackID = track->GetID();
      if(trackID < 0) trackID = -trackID - 1;

      for(Int_t j=0; j<nOmega; j++) {

        AliAODcascade *xi = (AliAODcascade*)fOmegaArray->At(j);
        if(!xi) continue;

        if(trackID == xi->GetBachID() || trackID == xi->GetPosID() || trackID == xi->GetNegID()) {

          fOmegaArray->RemoveAt(j);
          nShared++;
        }
      }
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedTracksProtonOmega"))->Fill(nShared);

    nShared = 0;

    // Lambda-Lambda
    for(Int_t i=0; i<nLambda; i++) {

      AliAODv0 *v01 = (AliAODv0*)fLambdaArray->At(i);
      if(!v01) continue;

      for(Int_t j=i+1; j<nLambda; j++) {

        AliAODv0 *v02 = (AliAODv0*)fLambdaArray->At(j);
        if(!v02) continue;

        if(v01->GetPosID() == v02->GetPosID() || v01->GetNegID() == v02->GetNegID()) {

          if(v01->DcaV0Daughters() < v02->DcaV0Daughters()) fLambdaArray->RemoveAt(j);
          else                                              fLambdaArray->RemoveAt(i);
          nShared++;
        }
      }
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedTracksLambdaLambda"))->Fill(nShared);

    nShared = 0;

    // Lambda-Xi
    for(Int_t i=0; i<nLambda; i++) {

      AliAODv0 *v0 = (AliAODv0*)fLambdaArray->At(i);
      if(!v0) continue;

      for(Int_t j=0; j<nXi; j++) {

        AliAODcascade *xi = (AliAODcascade*)fXiArray->At(j);
        if(!xi) continue;

        if(v0->GetPosID() == xi->GetPosID() || v0->GetNegID() == xi->GetNegID() || v0->GetNegID() == xi->GetBachID()) {

          fXiArray->RemoveAt(j);
          nShared++;
        }
      }
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedTracksLambdaXi"))->Fill(nShared);

    nShared = 0;

    // Lambda-Omega
    for(Int_t i=0; i<nLambda; i++) {

      AliAODv0 *v0 = (AliAODv0*)fLambdaArray->At(i);
      if(!v0) continue;

      for(Int_t j=0; j<nOmega; j++) {

        AliAODcascade *xi = (AliAODcascade*)fOmegaArray->At(j);
        if(!xi) continue;

        if(v0->GetPosID() == xi->GetPosID() || v0->GetNegID() == xi->GetNegID() || v0->GetNegID() == xi->GetBachID()) {

          fOmegaArray->RemoveAt(j);
          nShared++;
        }
      }
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedTracksLambdaOmega"))->Fill(nShared);

    nShared = 0;

    // Xi-Xi
    for(Int_t i=0; i<nXi; i++) {

      AliAODcascade *xi1 = (AliAODcascade*)fXiArray->At(i);
      if(!xi1) continue;

      for(Int_t j=i+1; j<nXi; j++) {

        AliAODcascade *xi2 = (AliAODcascade*)fXiArray->At(j);
        if(!xi2) continue;

        if(xi1->GetPosID()  == xi2->GetPosID()  || xi1->GetNegID()  == xi2->GetNegID()
        || xi1->GetNegID()  == xi2->GetBachID() || xi1->GetBachID() == xi2->GetNegID()
        || xi1->GetBachID() == xi2->GetBachID()) {

          if(xi1->DcaXiDaughters() < xi2->DcaXiDaughters()) fXiArray->RemoveAt(j);
          else                                              fXiArray->RemoveAt(i);
          nShared++;
        }
      }
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedTracksXiXi"))->Fill(nShared);

    nShared = 0;

    // Xi-Omega
    for(Int_t i=0; i<nXi; i++) {

      AliAODcascade *xi1 = (AliAODcascade*)fXiArray->At(i);
      if(!xi1) continue;

      for(Int_t j=0; j<nOmega; j++) {

        AliAODcascade *xi2 = (AliAODcascade*)fOmegaArray->At(j);
        if(!xi2) continue;

        if(xi1->GetPosID()  == xi2->GetPosID()  || xi1->GetNegID()  == xi2->GetNegID()
        || xi1->GetNegID()  == xi2->GetBachID() || xi1->GetBachID() == xi2->GetNegID()
        || xi1->GetBachID() == xi2->GetBachID()) {

          fOmegaArray->RemoveAt(j);
          nShared++;
        }
      }
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedTracksXiOmega"))->Fill(nShared);

    nShared = 0;

    // Omega-Omega
    for(Int_t i=0; i<nOmega; i++) {

      AliAODcascade *xi1 = (AliAODcascade*)fOmegaArray->At(i);
      if(!xi1) continue;

      for(Int_t j=i+1; j<nOmega; j++) {

        AliAODcascade *xi2 = (AliAODcascade*)fOmegaArray->At(j);
        if(!xi2) continue;

        if(xi1->GetPosID()  == xi2->GetPosID()  || xi1->GetNegID()  == xi2->GetNegID()
        || xi1->GetNegID()  == xi2->GetBachID() || xi1->GetBachID() == xi2->GetNegID()
        || xi1->GetBachID() == xi2->GetBachID()) {

          if(xi1->DcaXiDaughters() < xi2->DcaXiDaughters()) fOmegaArray->RemoveAt(j);
          else                                              fOmegaArray->RemoveAt(i);
          nShared++;
        }
      }
    }
    dynamic_cast<TH1F*>(fOutput->FindObject("hNSharedTracksOmegaOmega"))->Fill(nShared);

    fProtonArray->Compress();
    fLambdaArray->Compress();
    fXiArray->Compress();
    fOmegaArray->Compress();
  }
}
//_______________________________________________________________________________________________
Double_t AliAnalysisTaskDibaryons::relKcalc(TLorentzVector track1,TLorentzVector track2)
{
  //This function calculates the relative momentum k* between any particle pair

  TLorentzVector trackSum, track1cms, track2cms;
  trackSum = track1 + track2;

  Double_t beta = trackSum.Beta();
  Double_t betax = beta*cos(trackSum.Phi())*sin(trackSum.Theta());
  Double_t betay = beta*sin(trackSum.Phi())*sin(trackSum.Theta());
  Double_t betaz = beta*cos(trackSum.Theta());

  track1cms = track1;
  track2cms = track2;

  track1cms.Boost(-betax,-betay,-betaz);
  track2cms.Boost(-betax,-betay,-betaz);

  TLorentzVector trackRelK;

  trackRelK = track1cms - track2cms;
  Double_t relK = 0.5*trackRelK.P();

  return relK;
}
//_______________________________________________________________________________________________
void AliAnalysisTaskDibaryons::Terminate(Option_t *option)
{
  //Called once at the end of the query
  AliInfo(Form("%s is done.",GetName()));
}
//_______________________________________________________________________________________________
