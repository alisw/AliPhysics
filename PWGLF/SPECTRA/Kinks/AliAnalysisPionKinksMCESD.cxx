/**************************************************************************
 * Authors: Eftychios Cheiladakis, for his Master Thesis                  *
 * at the  Physics Department of Athens University                        *
 * under the supervision of  Prof. Martha Spyropoulou-Stassinaki          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//                 AliAnalysisPionKinksESDMC class
//       Example of an analysis task for kink topology study
//      pions from kink topology are 'identified' in this code
//-----------------------------------------------------------------

#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TF1.h" 
#include "TH1.h" 
#include "TH2.h" 
#include "TH3.h"
#include "TList.h"
#include "TParticle.h"

#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDkink.h"
#include "AliESDpid.h"
#include "AliPID.h"
#include "AliStack.h"

#include "AliAnalysisTask.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliPIDResponse.h" 
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisPionKinksMCESD.h"

ClassImp(AliAnalysisPionKinksMCESD)

//________________________________________________________________________
AliAnalysisPionKinksMCESD::AliAnalysisPionKinksMCESD(const char *name)
:AliAnalysisTaskSE(name),
fMaxKinkAngKmu(0), 
fMaxKinkAngPimu(0), //functions
hMCMult(0),
hMCMultPrim(0),
hMCPtAll(0),
hMCEtaAll(0),
hMCPtPrim(0), 
hMCEtaPrim(0), 
hMCPt(0), 
hMCEta(0), 
hMCPdg(0), 
hMCMultPiPlus(0), 
hMCPtPiPlus(0), 
hMCEtaPiPlus(0), 
hMCRapidityPiPlus(0),
hMCNDaughtersPlus(0),
hMCRadPiDauPlus(0), 
hMCKinkPosZPlus(0),
hMCUIDPiDauPlus(0), 
hMCPdgPiNonDecayedPlus(0), 
hMCPdgPiDauPlus(0),
hMCQtPlus(0), 
hMCKinkAnglePlus(0),
hMCPKinkAngPlus(0),
hMCPdgCodemdPlus(0), 
hMCPtmdPlus(0), 
hMCPtPimuonPlus(0),
hMCEtaPimuonPlus(0),
hMCRapidityPimuonPlus(0),
hMCQtPimuonPlus(0), 
hMCPKinkAngPimuonPlus(0),
hMCPtPiotherPlus(0),
hMCEtaPiotherPlus(0),
hMCRapidityPiotherPlus(0),
hMCQtPiotherPlus(0),
hMCPKinkAngPiotherPlus(0),
hMCMultPiMinus(0), 
hMCPtPiMinus(0), 
hMCEtaPiMinus(0), 
hMCRapidityPiMinus(0),
hMCNDaughtersMinus(0),
hMCRadPiDauMinus(0), 
hMCKinkPosZMinus(0),
hMCUIDPiDauMinus(0), 
hMCPdgPiNonDecayedMinus(0), 
hMCPdgPiDauMinus(0),
hMCQtMinus(0), 
hMCKinkAngleMinus(0),
hMCPKinkAngMinus(0),
hMCPdgCodemdMinus(0), 
hMCPtmdMinus(0), 
hMCPtPimuonMinus(0),
hMCEtaPimuonMinus(0),
hMCRapidityPimuonMinus(0),
hMCQtPimuonMinus(0), 
hMCPKinkAngPimuonMinus(0),
hMCPtPiotherMinus(0),
hMCEtaPiotherMinus(0),
hMCRapidityPiotherMinus(0),
hMCQtPiotherMinus(0),
hMCPKinkAngPiotherMinus(0),//MC histograms
hMult(0),
hAcceptedMult(0),
hMultPS(0),
hvtx(0),
hvtxy(0),
hvtyz(0),
hvtxz(0),
hMultPSV(0),
hPtAll(0),
hEtaAll(0),
hTrackPos(0),
hTrackPosxy(0),
hTrackPosyz(0),
hTrackPosxz(0),
//hTPCchi2clusters(0),
//hdcaToVertexXY(0),
//hdcaToVertexZ(0),
hMultPrim(0),
hPtPrim(0),
hEtaPrim(0),
hPrimTrackPos(0),
hPrimTrackPosxy(0),
hPrimTrackPosyz(0),
hPrimTrackPosxz(0),
hPt(0),
hEta(0),
//hRapidity(0),
hPtKink(0),
hEtaKink(0),
hRapidityKink(0),
hPmP(0),
hKinkPosRTPCclusters1(0),
hKinkPosRTPCclusters2(0),
hQt(0),
hKinkAngle(0),
hDCAkink(0),
hPmKinkAng(0),
hKinkPosXY(0),
hKinkPosZY(0),
hKinkPosZR(0),
hKinkPosR(0),
hKinkPosZ(0),
hKinkPosZMCKinkPosZ(0),
hPdgCodemd(0),
hPmd(0),
hMinvPimu(0),
hUIDKinkDau(0),
hdEdx(0),
hPtKinkFake(0),
hEtaKinkFake(0),
hRapidityKinkFake(0),
hPmPFake(0),
hKinkPosRTPCclusters1Fake(0),
hKinkPosRTPCclusters2Fake(0),
hQtFake(0),
hKinkAngleFake(0),
hDCAkinkFake(0),
hPmKinkAngFake(0),
hKinkPosXYFake(0),
hKinkPosZYFake(0),
hKinkPosZRFake(0),
hKinkPosRFake(0),
hKinkPosZFake(0),
hKinkPosZMCKinkPosZFake(0),
hPdgCodemdFake(0),
hPmdFake(0),
hMinvPimuFake(0),
hUIDKinkDauFake(0),
hdEdxFake(0),
hPtPosRSelected(0),
hPdgCodemdZRejected(0),
hPtZSelected(0),
hPdgCodemdAngRejected(0),
hPtAngSelected(0),
hPdgCodemdPmRejected(0),
hPtPmSelected(0),
hPdgCodemdQtLowRejected(0),
hPtGoodKink(0), 
hEtaGoodKink(0), 
hRapidityGoodKink(0), 
hQtGoodKink(0), 
hPmGoodKinkAng(0),  
hPdgCodemdGoodKink(0), 
hPmdGoodKink(0),
hUIDGoodKinkDau(0),
hdEdxGoodKink(0),
hUIDPiDauPlus(0), 
hMultPiPlus(0),
hPtPiPlus(0),
hEtaPiPlus(0),
hRapidityPiPlus(0),
hQtPiPlus(0),
hKinkAnglePiPlus(0),
hPmKinkAngPiPlus(0),
hKinkPosXYPiPlus(0),
hKinkPosZRPiPlus(0),
hKinkPosRPiPlus(0),
hDCAkinkPiPlus(0),
hPdgCodemdPiPlus(0),
hPmdPiPlus(0),
hdEdxPiPlus(0),
hQtPimuPlus(0),
hPmKinkAngPimuPlus(0),
hQtPiotherPlus(0),
hPmKinkAngPiotherPlus(0), 
hPdgCodemdPiotherPlus(0),
hUIDPiDauMinus(0),
hMultPiMinus(0),
hPtPiMinus(0),
hEtaPiMinus(0),
hRapidityPiMinus(0),
hQtPiMinus(0),
hKinkAnglePiMinus(0),
hPmKinkAngPiMinus(0),
hKinkPosXYPiMinus(0),
hKinkPosZRPiMinus(0),
hKinkPosRPiMinus(0),
hDCAkinkPiMinus(0),
hPdgCodemdPiMinus(0),
hPmdPiMinus(0),
hdEdxPiMinus(0),
hQtPimuMinus(0),
hPmKinkAngPimuMinus(0),
hQtPiotherMinus(0),
hPmKinkAngPiotherMinus(0), 
hPdgCodemdPiotherMinus(0),
hPdgCodemdQtRejected(0),
hPtQtSelected(0),
hPdgCodemdMaxAngRejected(0),
hPtMaxAngSelected(0),
hPdgCodemdRTPCclustersRejected(0),
hPtRTPCclustersSelected(0),
hRTPCclustersRTPCclustersSelected(0),
hPdgCodemdMinvRejected(0),
hPtSelected(0), 
hEtaSelected(0), 
hRapiditySelected(0), 
hQtSelected(0), 
hKinkAngleSelected(0),
hDCAkinkSelected(0),
hPmKinkAngSelected(0), 
hKinkPosXYSelected(0), 
hKinkPosZRSelected(0), 
hKinkPosRSelected(0), 
hPdgCodemdSelected(0), 
hPmdSelected(0),
hMinvPimuSelected(0),  
hUIDKinkDauSelected(0), 
hdEdxSelected(0), 
hPtSelectedFake(0), 
hEtaSelectedFake(0), 
hRapiditySelectedFake(0), 
hQtSelectedFake(0), 
hKinkAngleSelectedFake(0),
hDCAkinkSelectedFake(0),
hPmKinkAngSelectedFake(0),
hKinkPosXYSelectedFake(0),
hKinkPosZRSelectedFake(0),
hKinkPosRSelectedFake(0),
hPmdSelectedFake(0),
hMinvPimuSelectedFake(0), 
hdEdxSelectedFake(0),
hPdgCodemddEdxRejected(0),
hPtPiSelected(0), 
hEtaPiSelected(0), 
hRapidityPiSelected(0), 
hQtPiSelected(0), 
hKinkAnglePiSelected(0),
hDCAkinkPiSelected(0),
hPmKinkAngPiSelected(0), 
hKinkPosRTPCclusters1PiSelected(0),
hKinkPosRTPCclusters2PiSelected(0),
hKinkPosXYPiSelected(0), 
hKinkPosZRPiSelected(0), 
hKinkPosRPiSelected(0), 
hKinkPosZPiSelected(0),
hPmPPiSelected(0),
hPdgCodemdPiSelected(0), 
hPmdPiSelected(0),
hMinvPimuPiSelected(0),  
hUIDKinkDauPiSelected(0), 
hdEdxPiSelected(0), 
hPtPiSelectedPlus(0), 
hEtaPiSelectedPlus(0), 
hRapidityPiSelectedPlus(0), 
hQtPiSelectedPlus(0), 
hKinkAnglePiSelectedPlus(0),
hDCAkinkPiSelectedPlus(0),
hPmKinkAngPiSelectedPlus(0), 
hKinkPosXYPiSelectedPlus(0), 
hKinkPosZRPiSelectedPlus(0),  
hPdgCodemdPiSelectedPlus(0), 
hPmdPiSelectedPlus(0),
hMinvPimuPiSelectedPlus(0),  
hUIDPiDaumuSelectedPlus(0), 
hdEdxPiSelectedPlus(0),
hPtPrimPiKinksPlus(0), 
hEtaPrimPiKinksPlus(0), 
hRapidityPrimPiKinksPlus(0),
hPtSecondPiKinksPlus(0), 
hEtaSecondPiKinksPlus(0), 
hRapiditySecondPiKinksPlus(0),
hPtNonPiKinksPlus(0), 
hEtaNonPiKinksPlus(0), 
hRapidityNonPiKinksPlus(0), 
hPdgCodemdNonPiKinksPlus(0), 
hPtPiSelectedMinus(0),
hEtaPiSelectedMinus(0), 
hRapidityPiSelectedMinus(0), 
hQtPiSelectedMinus(0), 
hKinkAnglePiSelectedMinus(0),
hDCAkinkPiSelectedMinus(0),
hPmKinkAngPiSelectedMinus(0), 
hKinkPosXYPiSelectedMinus(0), 
hKinkPosZRPiSelectedMinus(0),  
hPdgCodemdPiSelectedMinus(0), 
hPmdPiSelectedMinus(0),
hMinvPimuPiSelectedMinus(0), 
hUIDPiDaumuSelectedMinus(0), 
hdEdxPiSelectedMinus(0),
hPtPrimPiKinksMinus(0), 
hEtaPrimPiKinksMinus(0), 
hRapidityPrimPiKinksMinus(0),
hPtSecondPiKinksMinus(0), 
hEtaSecondPiKinksMinus(0), 
hRapiditySecondPiKinksMinus(0),
hPtNonPiKinksMinus(0), 
hEtaNonPiKinksMinus(0), 
hRapidityNonPiKinksMinus(0), 
hPdgCodemdNonPiKinksMinus(0),// reconstruction histograms
fListOfHistos(0),
fLowMulcut(-1), fUpMulcut(-1), 
cLowPt(0), cRapidityLim(0),
cLowR(0), cUpR(0),
cLowZ(0), cUpZ(0),
cLowKinkAngle(0),
cLowQt(0), cUpQt(0),
cLowInvMass(0), cUpInvMass(0),
cSigmaCut(0),
cPdgKaon(321), cPdgPion(211), cPdgMuon(13), cPdgElectron(11),
cKaonMass(0), cPionMass(0), cMuonMass(0), cElectronMass(0),
nBinsMult(0), hLowMult(0), hUpMult(0),
nBinsPt(0), hLowPt(0), hUpPt(0),
nBinsEta(0), hLowEta(0), hUpEta(0),
nBinsQt(0), hLowQt(0), hUpQt(0),
nBinsPdg(0), hLowPdg(0), hUpPdg(0),
nBinsPdg2(0), hLowPdg2(0), hUpPdg2(0),
nBinsUID(0), hLowUID(0), hUpUID(0),
nBinsR(0), hLowR(0), hUpR(0),
nBinsZ(0), hLowZ(0), hUpZ(0),
nBinsXY(0), hLowXY(0), hUpXY(0),
nBinsAngle(0), hLowAngle(0), hUpAngle(0),
nBinsZV(0), hLowZV(0), hUpZV(0),
nBinsXYV(0), hLowXYV(0), hUpXYV(0),
nBinsInvMass(0), hLowInvMass(0), hUpInvMass(0),
nBinsdEdx(0), hLowdEdx(0), hUpdEdx(0), fPIDResponse(0),
fMaxDCAtoVtxCut(0), fTrackCuts(0)

{
//Multiplicity bins 
fMaxDCAtoVtxCut=new AliESDtrackCuts("fMaxDCAtoVtxCut","fMaxDCAtoVtxCut");
fMaxDCAtoVtxCut->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
fMaxDCAtoVtxCut->SetMaxChi2TPCConstrainedGlobal(36);

fTrackCuts = new AliESDtrackCuts("Multiplicity bins","Multiplicity bins");
fTrackCuts->SetMinNClustersTPC(70);
fTrackCuts->SetMaxChi2PerClusterTPC(4);
fTrackCuts->SetAcceptKinkDaughters(kFALSE); 
fTrackCuts->SetRequireTPCRefit(kTRUE);
fTrackCuts->SetRequireITSRefit(kTRUE);
fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
fTrackCuts->SetMaxDCAToVertexZ(2);
fTrackCuts->SetDCAToVertex2D(kFALSE);
fTrackCuts->SetRequireSigmaToVertex(kFALSE);
fTrackCuts->SetEtaRange(-0.8,+0.8);
fTrackCuts->SetPtRange(0.15, 1e10);

//DefineOutput(0, TList::Class());
DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisPionKinksMCESD::UserCreateOutputObjects() {
fListOfHistos=new TList();

//maximum kink angle for kaons to muons 
fMaxKinkAngKmu=new TF1("fMaxKinkAngKmu","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",1.1,10.0);
fMaxKinkAngKmu->SetParameter(0,cKaonMass); 
fMaxKinkAngKmu->SetParameter(1,0.9127037);
fMaxKinkAngKmu->SetParameter(2,TMath::Pi());

//maximum kink angle for pions to muons 
fMaxKinkAngPimu=new TF1("fMaxKinkAngPimu","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.1,10.0);
fMaxKinkAngPimu->SetParameter(0,cPionMass);
fMaxKinkAngPimu->SetParameter(1,0.2731374);
fMaxKinkAngPimu->SetParameter(2,TMath::Pi());

//Create histograms
TH1::SetDefaultSumw2();
TH2::SetDefaultSumw2();

//MC histograms
hMCMult = new TH1F("hMCMult", "MC multiplicity; Number of tracks; Number of events", 100, 0.0, 2000);
hMCMultPrim = new TH1F("hMCMultPrim", "MC primary tracks multiplicity; Number of tracks; Number of events", nBinsMult, hLowMult, hUpMult);
hMCPtAll = new TH1F("hMCPtAll", "Transverse momentum of all MC tracks; p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hMCEtaAll = new TH1F("hMCEtaAll", "Pseudorapidity of all MC tracks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCPtPrim = new TH1F("hMCPtPrim", "Transverse momentum of primary MC tracks; p_{T} (GeV/c); dN/dp_{T}",nBinsPt, hLowPt, hUpPt);
hMCEtaPrim = new TH1F("hMCEtaPrim", "Pseudorapidity of primary MC tracks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCPt = new TH1F("hMCPt", "Transverse momentum of selected MC tracks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hMCEta = new TH1F("hMCEta", "Pseudorapidity of selected MC tracks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCPdg = new TH1F("hMCPdg", "Pdg code of selected MC tracks; Pdg code; Number of particles", nBinsPdg, hLowPdg, hUpPdg);
hMCMultPiPlus = new TH1F("hMCMultPiPlus", "MC pion multiplicity; Number of pion tracks; Number of events", nBinsMult, hLowMult, hUpMult);
hMCPtPiPlus = new TH1F("hMCPtPiPlus", "Transverse momentum of selected MC pions; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hMCEtaPiPlus = new TH1F("hMCEtaPiPlus", "Pseudorapidity of selected MC piaons; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCRapidityPiPlus = new TH1F("hMCRapidityPiPlus", "Pseudorapidity of selected MC pions; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCNDaughtersPlus = new TH1F("hMCNDaughtersPlus", "Number of daughters; number of daughers; number of pions", 10,0,10);
hMCRadPiDauPlus = new TH1F("hMCRadPiDauPlus", "Radius of MC daughter generation position; R (cm); dN/dR", nBinsR, hLowR, hUpR);
hMCKinkPosZPlus = new TH1F("hMCKinkPosZPlus", "z position of MC daughter generation vertex; z (cm); dN/dz", 100, 0.0, 500.0);
hMCUIDPiDauPlus = new TH1F("hMCUIDPiDauPlus", "UID (method of production) of MC pion daughters; UID (method of production); Number of particles", 21, 0.0, 20);
hMCPdgPiNonDecayedPlus = new TH1F("hMCPdgPiNonDecayedPlus", "Pdg code of MC pion non-decayed products; Pdg code; Number of particles", nBinsPdg, hLowPdg, hUpPdg);
hMCPdgPiDauPlus = new TH1F("hMCPdgPiDauPlus", "Pdg code of MC pion daughters; Pdg code; Number of particles", nBinsPdg, hLowPdg, hUpPdg);
hMCQtPlus = new TH1F("hMCQtPlus", "Daughter's transverse momentum in mother's frame for MC pion daughters; q_{T} (GeV/c); dN/dq_{T}", nBinsQt, hLowQt, hUpQt);
hMCKinkAnglePlus = new TH1F("hMCKinkAnglePlus", "MC angle between pion mother's and daughter's momentum; #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hMCPKinkAngPlus = new TH2F("hMCPKinkAngPlus", "MC mother's P vs kink angle; P (GeV/c); #theta (#circ); dN/d#theta",   nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hMCPdgCodemdPlus = new TH2F("hMCPdgCodemdPlus", "MC mother vs daughter pdg code; Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hMCPtmdPlus = new TH2F("hMCPtmdPlus", "MC mother vs daughter transverse momentum; Mother's p_{T} (GeV/c); Daughter's p_{T} (GeV/c)",   nBinsPt, hLowPt, hUpPt,   nBinsPt, hLowPt, hUpPt);
hMCPtPimuonPlus = new TH1F("hMCPtPimuonPlus", "Transverse momentum of selected pions MC decaying to muons; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hMCEtaPimuonPlus = new TH1F("hMCEtaPimuonPlus", "Pseudorapidity; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCRapidityPimuonPlus = new TH1F("hMCRapidityPimuonPlus", "Rapidity of selected MC tracks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCQtPimuonPlus = new TH1F("hMCQtPimuonPlus", "MC daughter's (muon) transverse momentum in mother's (pion) frame; q_{T} (GeV/c); dN/dq_{T}", nBinsQt, hLowQt, hUpQt);
hMCPKinkAngPimuonPlus = new TH1F("hMCPKinkAngPimuonPlus", "MC angle between pion mother's and muon daughter's momentum; #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hMCPtPiotherPlus = new TH1F("hMCPtPiotherPlus", "Transverse momentum of selected pions MC not decaying to muons; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hMCEtaPiotherPlus = new TH1F("hMCEtaPiotherPlus", "Pseudorapidity; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCRapidityPiotherPlus = new TH1F("hMCRapidityPiotherPlus", "Rapidity of selected MC tracks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCQtPiotherPlus = new TH1F("hMCQtPiotherPlus", "MC daughter's (muon) transverse momentum in mother's (pion) frame; q_{T} (GeV/c); dN/dq_{T}", nBinsQt, hLowQt, hUpQt);
hMCPKinkAngPiotherPlus = new TH1F("hMCPKinkAngPiotherPlus", "MC angle between pion mother's and non-muon daughter's momentum; #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hMCMultPiMinus = new TH1F("hMCMultPiMinus", "MC pion multiplicity; Number of pion tracks; Number of events", nBinsMult, hLowMult, hUpMult);
hMCPtPiMinus = new TH1F("hMCPtPiMinus", "Transverse momentum of selected MC pions; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hMCEtaPiMinus = new TH1F("hMCEtaPiMinus", "Pseudorapidity of selected MC pions; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCRapidityPiMinus = new TH1F("hMCRapidityPiMinus", "Pseudorapidity of selected MC pions; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCNDaughtersMinus = new TH1F("hMCNDaughtersMinus", "Number of daughters; number of daughers; number of pions", 10,0,10);
hMCRadPiDauMinus = new TH1F("hMCRadPiDauMinus", "Radius of MC daughter generation position; R (cm); dN/dR", nBinsR, hLowR, hUpR);
hMCKinkPosZMinus = new TH1F("hMCKinkPosZMinus", "z position of MC daughter generation vertex; z (cm); dN/dz", 100, 0.0, 500.0);
hMCUIDPiDauMinus = new TH1F("hMCUIDPiDauMinus", "UID (method of production) of MC pion daughters; UID (method of production); Number of particles", 21, 0.0, 20);
hMCPdgPiNonDecayedMinus = new TH1F("hMCPdgPiNonDecayedMinus", "Pdg code of MC pion non-decayed products; Pdg code; Number of particles", nBinsPdg, hLowPdg, hUpPdg);
hMCPdgPiDauMinus = new TH1F("hMCPdgPiDauMinus", "Pdg code of MC pion daughters; Pdg code; Number of particles", nBinsPdg, hLowPdg, hUpPdg);
hMCQtMinus = new TH1F("hMCQtMinus", "Daughter's transverse momentum in mother's frame for MC pion daughters; q_{T} (GeV/c); dN/dq_{T}", nBinsQt, hLowQt, hUpQt);
hMCKinkAngleMinus = new TH1F("hMCKinkAngleMinus", "MC angle between pion mother's and daughter's momentum; #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hMCPKinkAngMinus = new TH2F("hMCPKinkAngMinus", "MC mother's P vs kink angle; P (GeV/c); #theta (#circ); dN/d#theta",   nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hMCPdgCodemdMinus = new TH2F("hMCPdgCodemdMinus", "MC mother vs daughter pdg code; Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hMCPtmdMinus = new TH2F("hMCPtmdMinus", "MC mother vs daughter transverse momentum; Mother's p_{T} (GeV/c); Daughter's p_{T} (GeV/c)",   nBinsPt, hLowPt, hUpPt,   nBinsPt, hLowPt, hUpPt);
hMCPtPimuonMinus = new TH1F("hMCPtPimuonMinus", "Transverse momentum of selected pions MC decaying to muons; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hMCEtaPimuonMinus = new TH1F("hMCEtaPimuonMinus", "Pseudorapidity; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCRapidityPimuonMinus = new TH1F("hMCRapidityPimuonMinus", "Rapidity of selected MC tracks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCQtPimuonMinus = new TH1F("hMCQtPimuonMinus", "MC daughter's (muon) transverse momentum in mother's (pion) frame; q_{T} (GeV/c); dN/dq_{T}", nBinsQt, hLowQt, hUpQt);
hMCPKinkAngPimuonMinus = new TH1F("hMCPKinkAngPimuonMinus", "MC angle betweenpion mother's and muon daughter's momentum; #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hMCPtPiotherMinus = new TH1F("hMCPtPiotherMinus", "Transverse momentum of selected pions MC not decaying to muons; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hMCEtaPiotherMinus = new TH1F("hMCEtaPiotherMinus", "Pseudorapidity; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCRapidityPiotherMinus = new TH1F("hMCRapidityPiotherMinus", "Rapidity of selected MC tracks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hMCQtPiotherMinus = new TH1F("hMCQtPiotherMinus", "MC daughter's (muon) transverse momentum in mother's (pion) frame; q_{T} (GeV/c); dN/dq_{T}", nBinsQt, hLowQt, hUpQt);
hMCPKinkAngPiotherMinus = new TH1F("hMCPKinkAngPiotherMinus", "MC angle betweenpion mother's and non-muon daughter's momentum; #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);


//Reconstruction histograms
hMult = new TH1F("hMult", "Multiplicity (unbiased); Number of tracks; Number of events", nBinsMult, hLowMult, hUpMult);
hAcceptedMult = new TH1F("hAcceptedMult", "Multiplicity (biased); Number of tracks; Number of events", nBinsMult, hLowMult, hUpMult);
hMultPS = new TH1F("hMultPS", "Multiplicity after physics selection; Number of tracks; Number of events", nBinsMult, hLowMult, hUpMult);
hvtx = new TH3F("hvtx", "Reconstructed primary vertex position; x axis; y axis; z axis", nBinsXYV, hLowXYV, hUpXYV, nBinsXYV, hLowXYV, hUpXYV, nBinsZV, hLowZV, hUpZV);
hvtxy = new TH2F("hvtxy", "Reconstructed primary vertex position in x-y plane; x axis; y axis", nBinsXYV, hLowXYV, hUpXYV, nBinsXYV, hLowXYV, hUpXYV);
hvtyz = new TH2F("hvtyz", "Reconstructed primary vertex position in y-z plane; y axis; z axis", nBinsXYV, hLowXYV, hUpXYV, nBinsZV, hLowZV, hUpZV);
hvtxz = new TH2F("hvtxz", "Reconstructed primary vertex position in x-z plane; x axis; z axis", nBinsXYV, hLowXYV, hUpXYV, nBinsZV, hLowZV, hUpZV);
hMultPSV = new TH1F("hMultPSV", "Multiplicity after physics selection & vertex cut; Number of tracks; Number of events", nBinsMult, hLowMult, hUpMult);
hPtAll = new TH1F("hPtAll", "Transverse momentum of all tracks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaAll = new TH1F("hEtaAll", "Pseudorapidity of all tracks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hTrackPos = new TH3F("hTrackPos", "Generetion position of all reconstructed tracks", nBinsXYV, hLowXYV, hUpXYV, nBinsXYV, hLowXYV, hUpXYV, nBinsZV, hLowZV, hUpZV);
hTrackPosxy = new TH2F("hTrackPosxy", "Generetion position of all reconstructed tracks in x-y plane; x axis; y axis", nBinsXYV, hLowXYV, hUpXYV, nBinsXYV, hLowXYV, hUpXYV);
hTrackPosyz = new TH2F("hTrackPosyz", "Generetion position of all reconstructed tracks in y-z plane; y axis; z axis", nBinsXYV, hLowXYV, hUpXYV, nBinsZV, hLowZV, hUpZV);
hTrackPosxz = new TH2F("hTrackPosxz", "Generetion position of all reconstructed tracks in x-z plane; x axis; z axis", nBinsXYV, hLowXYV, hUpXYV, nBinsZV, hLowZV, hUpZV);
//hTPCchi2clusters = new TH1F("hTPCchi2clusters", "#chi^{2}/Number of TPC clusters; #chi^{2}; TPC clusters)", 100, 0.0, 2.0);//
//hdcaToVertexXY = new TH1F("hdcaToVertexXY", "Track to vertex impact parameter in x-y plane; DCA_{z} (cm); dN/d(DCA_{z})", 100, 0.0, 2.0);//
//hdcaToVertexZ = new TH1F("hdcaToVertexZ", "Track to vertex impact parameter in z axis; DCA_{z} (cm); dN/d(DCA_{z})", 100, 0.0, 2.0);//
hMultPrim = new TH1F("hMultPrim", "ESD primary - supposed tracks multiplicity; Number of tracks; Number of events", nBinsMult, hLowMult, hUpMult);
hPtPrim = new TH1F("hPtPrim", "Transverse momentum of primary - supposed ESD tracks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPrim = new TH1F("hEtaPrim", "Pseudorapidity of primary - supposed ESD tracks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hPrimTrackPos = new TH3F("hPrimTrackPos", "Generetion position of selected tracks (DCA and quality cuts)", nBinsXYV, hLowXYV, hUpXYV, nBinsXYV, hLowXYV, hUpXYV, nBinsZV, hLowZV, hUpZV);
hPrimTrackPosxy = new TH2F("hPrimTrackPosxy", "Generetion position of selected tracks in x-y plane (DCA and quality cuts); x axis; y axis", nBinsXYV, hLowXYV, hUpXYV, nBinsXYV, hLowXYV, hUpXYV);
hPrimTrackPosyz = new TH2F("hPrimTrackPosyz", "Generetion position of selected tracks in y-z plane (DCA and quality cuts); y axis; z axis", nBinsXYV, hLowXYV, hUpXYV, nBinsZV, hLowZV, hUpZV);
hPrimTrackPosxz = new TH2F("hPrimTrackPosxz", "Generetion position of selected tracks in x-z plane (DCA and quality cuts); x axis; z axis", nBinsXYV, hLowXYV, hUpXYV, nBinsZV, hLowZV, hUpZV);
hPt = new TH1F("hPt", "Transverse momentum of selected ESD tracks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEta = new TH1F("hEta", "Pseudorapidity of selected ESD tracks; n; Number of tracks dN/dn", nBinsEta, hLowEta, hUpEta);
//hRapidity = new TH1F("hRapidity", "Rapidity of selected ESD tracks; n; Number of tracks dN/dn", nBinsEta, hLowEta, hUpEta);
hPtKink = new TH1F("hPtKink", "Mother's transverse momentum for all ESD kinks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaKink = new TH1F("hEtaKink", "Mother's pseudorapidity for all ESD kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityKink = new TH1F("hRapidityKink", "Mother's rapidity for all ESD kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hPmP = new TH2F("hPmP", "Mother's momentum as calculated by kink and by track; P_{kink}; P_{track}", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hKinkPosRTPCclusters1 = new TH1F("hKinkPosRTPCclusters1", " ;kinkposR; tpc clusters",100,0,1);
hKinkPosRTPCclusters2 = new TH2F("hKinkPosRTPCclusters2", " ;kinkposR; tpc clusters",nBinsR, hLowR, hUpR,100,0,200 );
hQt = new TH1F("hQt", "Daughter's transverse momentum in mother's frame for all ESD kinks; q_{T} (GeV/c); dN/dq_{T}", nBinsQt, hLowQt, hUpQt);
hKinkAngle = new TH1F("hKinkAngle", "Kink angle for all ESD kinks; #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hDCAkink = new TH1F("hDCAkink", "DCA between the two kink tracks; DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPmKinkAng = new TH2F("hPmKinkAng", "k, p_(GeV/c);  #theta (#circ); d^{2}N/dpd#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXY = new TH2F("hKinkPosXY", "X-Y Position of all kinks; X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZY = new TH2F("hKinkPosZY", "z vs y position of kinks (DCA, quality, p_{T} & y cuts); z (cm); y (cm)", nBinsZ, hLowZ, hUpZ, nBinsXY, hLowXY, hUpXY);	
hKinkPosZR = new TH2F("hKinkPosZR", "Z vs radius of all kinks position; Z (cm); R (cm)", nBinsZ, hLowZ, hUpZ, nBinsR, hLowR, hUpR);
hKinkPosR = new TH1F("hKinkPosR", "Position radius of all ESD kinks; R (cm); dN/dR", nBinsR, hLowR, hUpR);
hKinkPosZ = new TH1F("hKinkPosZ", "z position of kinks (DCA, quality, p_{T} & y cuts); z (cm); dN/dz", nBinsZ, -1, 1);
hKinkPosZMCKinkPosZ = new TH2F("hKinkPosZMCKinkPosZ", "Reconstructed vs generated z position of kinks (DCA, quality, p_{T} & y cuts); z (cm); dN/dz", nBinsZ, hLowZ, hUpZ, nBinsZ, hLowZ, hUpZ); 
hPdgCodemd = new TH2F("hPdgCodemd", "Mother vs daughter pdg code for all ESD kinks; Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPmd = new TH2F("hPmd", "ESD mother vs daughter momentum magnitude; Mother's P (GeV/c); Daughter's P (GeV/c)", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimu = new TH1F("hMinvPimu", "Invariant mass of ESD pions decaying to muons; m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hUIDKinkDau = new TH1F("hUIDKinkDau", "UID (method of production) of all ESD pion kink daughters (perfect PID); UID (method of production); Number of particles", 21, 0.0, 20);
hdEdx = new TH2F("hdEdx", "dE/dx vs mother's momentum; p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hPtKinkFake = new TH1F("hPtKinkFake", "Mother's transverse momentum for all fake kinks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaKinkFake = new TH1F("hEtaKinkFake", "Mother's pseudorapidity for all fake kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityKinkFake = new TH1F("hRapidityKinkFake", "Mother's rapidity for all fake kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hPmPFake = new TH2F("hPmPFake", "Mother's momentum magnitude calculated from kink vs calculated from track for all fake kinks; p_{kink} (GeV/c); p_{track} (GeV/c)",   nBinsPt, hLowPt, hUpPt,   nBinsPt, hLowPt, hUpPt);
hKinkPosRTPCclusters1Fake = new TH1F("hKinkPosRTPCclusters1Fake", " ;kinkposR; tpc clusters",100,0,1);
hKinkPosRTPCclusters2Fake = new TH2F("hKinkPosRTPCclusters2Fake", " ;kinkposR; tpc clusters",nBinsR, hLowR, hUpR,100,0,200 );
hQtFake = new TH1F("hQtFake", "Daughter's transverse momentum in mother's frame for all fake kinks; q_{T} (GeV/c); dN/dq_{T}", nBinsQt, hLowQt, hUpQt);
hKinkAngleFake = new TH1F("hKinkAngleFake", "Kink angle for all fake kinks; #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hDCAkinkFake = new TH1F("hDCAkinkFake", "DCA between the two fake kink tracks; DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPmKinkAngFake = new TH2F("hPmKinkAngFake", "k, p_(GeV/c);  #theta (#circ); d^{2}N/dpd#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXYFake = new TH2F("hKinkPosXYFake", "X-Y Position of all fake kinks; X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZYFake = new TH2F("hKinkPosZYFake", "z vs y position of fake kinks (DCA, quality, p_{T} & y cuts); z (cm); y (cm)", nBinsZ, hLowZ, hUpZ, nBinsXY, hLowXY, hUpXY);
hKinkPosZRFake = new TH2F("hKinkPosZRFake", "Z vs radius of all fake kinks position; Z (cm); R (cm)", nBinsZ, hLowZ, hUpZ, nBinsR, hLowR, hUpR);
hKinkPosRFake = new TH1F("hKinkPosRFake", "Position radius of all fake kinks; R (cm); dN/dR", nBinsR, hLowR, hUpR);
hKinkPosZFake = new TH1F("hKinkPosZFake", "z position of fake kinks (DCA, quality, p_{T} & y cuts); z (cm); dN/dz", nBinsZ, -1, 1);
hKinkPosZMCKinkPosZFake = new TH2F("hKinkPosZMCKinkPosZFake", "Reconstructed vs generated z position of fake kinks (DCA, quality, p_{T} & y cuts); z (cm); dN/dz", nBinsZ, hLowZ, hUpZ, nBinsZ, hLowZ, hUpZ); 
hPdgCodemdFake = new TH2F("hPdgCodemdFake", "Mother vs daughter pdg code for fakes; Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPmdFake = new TH2F("hPmdFake", "ESD mother vs daughter momentum magnitude (fake kinks); Mother's P (GeV/c); Daughter's P (GeV/c)", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuFake = new TH1F("hMinvPimuFake", "Invariant mass of ESD pions decaying to muons (fake kinks); m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hUIDKinkDauFake = new TH1F("hUIDKinkDauFake", "UID (method of production) of all ESD pion kink daughters (perfect PID); UID (method of production); Number of particles", 21, 0.0, 20);
hdEdxFake = new TH2F("hdEdxFake", "dE/dx vs mother's momentum; p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hPtPosRSelected = new TH1F("hPtPosRSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y & R cuts); R (cm); dN/dR", nBinsPt, hLowPt, hUpPt);
hPdgCodemdZRejected = new TH2F("hPdgCodemdZRejected", "Mother's vs daughter's pdg code for rejected kinks (DCA, quality, p_{T}, y, R & z cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPtZSelected = new TH1F("hPtZSelected",  "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R & z cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPdgCodemdAngRejected = new TH2F("hPdgCodemdAngRejected", "Mother's vs daughter's pdg code for rejected kinks (DCA, quality, p_{T}, y, R, z & #theta cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPtAngSelected = new TH1F("hPtAngSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z & #theta cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPdgCodemdPmRejected = new TH2F("hPdgCodemdPmRejected", "Mother's vs daughter's pdg code for rejected kinks (DCA, quality, p_{T}, y, R, z, #theta & p_{track}/p{kink} cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPtPmSelected = new TH1F("hPtPmSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z, #theta & p_{track}/p{kink} cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPdgCodemdQtLowRejected = new TH2F("hPdgCodemdQtLowRejected", "Mother's vs daughter's pdg code for rejected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink} & low q_{T} cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPtGoodKink = new TH1F("hPtGoodKink", "Mother's transverse momentum for real ESD kinks (realistic PID); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaGoodKink = new TH1F("hEtaGoodKink", "Mother's pseudorapidity for real ESD kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityGoodKink = new TH1F("hRapidityGoodKink", "Mother's rapidity for real ESD kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtGoodKink = new TH1F("hQtGoodKink", "Daughter's transverse momentum in mother's frame for real ESD kinks (realistic PID); q_{T} (GeV/c); dN/dq_{T}", 200, 0.0, 0.3);
hPmGoodKinkAng = new TH2F("hPmGoodKinkAng", "Mother's momentum magnitude vs kink angle for real ESD kinks (realistic PID); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hPdgCodemdGoodKink = new TH2F("hPdgCodemdGoodKink", "Mother vs daughter pdg code for real ESD kinks; Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPmdGoodKink = new TH2F("hPmdGoodKink", "Mother vs daughter momentum magnitude for real ESD kinks (realistic PID); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hUIDGoodKinkDau = new TH1F("hUIDGoodKinkDau", "UID (method of production); UID (method of production); Number of particles", 21, 0.0, 20);
hdEdxGoodKink = new TH2F("hdEdxGoodKink", "dE/dx vs mother's momentum; p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hUIDPiDauPlus = new TH1F("hUIDPiDauPlus", "UID (method of production) of all ESD pion kink daughters (perfect PID); UID (method of production); Number of particles", 21, 0.0, 20);
hMultPiPlus = new TH1F("hMultPiPlus", "ESD pion multiplicity in selected TPC area; Number of pions; Number of events", 100, 0.0, 200.0);
hPtPiPlus = new TH1F("hPtPiPlus", "Mother's transverse momentum for ESD pion kinks (perfect PID); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPiPlus = new TH1F("hEtaPiPlus", "Mother's pseudorapidity for ESD pion kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityPiPlus = new TH1F("hRapidityPiPlus", "Mother's rapidity for ESD pion kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtPiPlus = new TH1F("hQtPiPlus", "Daughter's transverse momentum in mother's frame for ESD pion kinks (perfect PID); q_{T} (GeV/c); dN/dq_{T}", 100, 0.0, 0.1);
hKinkAnglePiPlus = new TH1F("hKinkAnglePiPlus", "Kink angle of ESD pion kinks (perfect PID); #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hPmKinkAngPiPlus = new TH2F("hPmKinkAngPiPlus", "Mother's momentum magnitude vs kink angle for ESD pion kinks (perfect PID); p_{T} (GeV/c); #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXYPiPlus = new TH2F("hKinkPosXYPiPlus", "X-Y Position ; X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZRPiPlus = new TH2F("hKinkPosZRPiPlus", "Z vs radius of all fake kinks position; Z (cm); R (cm)", nBinsZ, hLowZ, hUpZ, nBinsR, hLowR, hUpR);
hKinkPosRPiPlus = new TH1F("hKinkPosRPiPlus", "Position radius of all fake kinks; R (cm); dN/dR", nBinsR, hLowR, hUpR);
hDCAkinkPiPlus = new TH1F("hDCAkinkPiPlus", "DCA; DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPdgCodemdPiPlus = new TH2F("hPdgCodemdPiPlus", "Mother vs daughter pdg code for pi kinks; Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPmdPiPlus = new TH2F("hPmdPiPlus", "Mother vs daughter momentum magnitude for pi kinks (realistic PID); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hdEdxPiPlus = new TH2F("hdEdxPiPlus", "dE/dx vs mother's momentum; p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hQtPimuPlus = new TH1F("hQtPimuPlus", "Daughter's transverse momentum in mother's frame for ESD pion to muon kinks (perfect PID); q_{T} (GeV/c); dN/dq_{T}", 100, 0.0, 0.1);
hPmKinkAngPimuPlus = new TH2F("hPmKinkAngPimuPlus", "Mother's momentum magnitude vs kink angle for ESD pion to muon kinks (perfect PID); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hQtPiotherPlus = new TH1F("hQtPiotherPlus", "Daughter's transverse momentum in mother's frame for ESD pion to other (not muon) kinks (perfect PID); q_{T} (GeV/c); dN/dq_{T}", 100, 0.0, 0.1);
hPmKinkAngPiotherPlus = new TH2F("hPmKinkAngPiotherPlus", "Mother's momentum magnitude vs kink angle for ESD pion to other (not muon) kinks (perfect PID); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hPdgCodemdPiotherPlus = new TH2F("hPdgCodemdPiotherPlus", "Mother vs daughter pdg code for ESD pion to other (not muon) kinks (perfect PID); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hUIDPiDauMinus = new TH1F("hUIDPiDauMinus", "UID (method of production) of all ESD pion kink daughters (perfect PID); UID (method of production); Number of particles", 21, 0.0, 20);
hMultPiMinus = new TH1F("hMultPiMinus", "ESD pion multiplicity in selected TPC area; Number of pions; Number of events", 100, 0.0, 200.0);
hPtPiMinus = new TH1F("hPtPiMinus", "Mother's transverse momentum for ESD pion kinks (perfect PID); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPiMinus = new TH1F("hEtaPiMinus", "Mother's pseudorapidity for ESD pion kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityPiMinus = new TH1F("hRapidityPiMinus", "Mother's rapidity for ESD pion kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtPiMinus = new TH1F("hQtPiMinus", "Daughter's transverse momentum in mother's frame for ESD pion kinks (perfect PID); q_{T} (GeV/c); dN/dq_{T}", 100, 0.0, 0.1);
hKinkAnglePiMinus = new TH1F("hKinkAnglePiMinus", "Kink angle of ESD pion kinks (perfect PID); #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hPmKinkAngPiMinus = new TH2F("hPmKinkAngPiMinus", "Mother's momentum magnitude vs kink angle for ESD pion kinks (perfect PID); p_{T} (GeV/c); #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXYPiMinus = new TH2F("hKinkPosXYPiMinus", "X-Y Position ; X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZRPiMinus = new TH2F("hKinkPosZRPiMinus", "Z vs radius of all fake kinks position; Z (cm); R (cm)", nBinsZ, hLowZ, hUpZ, nBinsR, hLowR, hUpR);
hKinkPosRPiMinus = new TH1F("hKinkPosRPiMinus", "Position radius of all fake kinks; R (cm); dN/dR", nBinsR, hLowR, hUpR);
hDCAkinkPiMinus = new TH1F("hDCAkinkPiMinus", "DCA; DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPdgCodemdPiMinus = new TH2F("hPdgCodemdPiMinus", "Mother vs daughter pdg code for pi kinks; Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPmdPiMinus = new TH2F("hPmdPiMinus", "Mother vs daughter momentum magnitude for pi kinks (realistic PID); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hdEdxPiMinus = new TH2F("hdEdxPiMinus", "dE/dx vs mother's momentum; p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hQtPimuMinus = new TH1F("hQtPimuMinus", "Daughter's transverse momentum in mother's frame for ESD pion to muon kinks (perfect PID); q_{T} (GeV/c); dN/dq_{T}", 100, 0.0, 0.1);
hPmKinkAngPimuMinus = new TH2F("hPmKinkAngPimuMinus", "Mother's momentum magnitude vs kink angle for ESD pion to muon kinks (perfect PID); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hQtPiotherMinus = new TH1F("hQtPiotherMinus", "Daughter's transverse momentum in mother's frame for ESD pion to other (not muon) kinks (perfect PID); q_{T} (GeV/c); dN/dq_{T}", 100, 0.0, 0.1);
hPmKinkAngPiotherMinus = new TH2F("hPmKinkAngPiotherMinus", "Mother's momentum magnitude vs kink angle for ESD pion to other (not muon) kinks (perfect PID); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hPdgCodemdPiotherMinus = new TH2F("hPdgCodemdPiotherMinus", "Mother vs daughter pdg code for ESD pion to other (not muon) kinks (perfect PID); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPdgCodemdQtRejected = new TH2F("hPdgCodemdQtRejected", "Mother's vs daughter's pdg code for rejected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink} & q^{T} cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPtQtSelected = new TH1F("hPtQtSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink} & q^{T} cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPdgCodemdMaxAngRejected = new TH2F("hPdgCodemdMaxAngRejected", "Mother's vs daughter's pdg code for rejected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink}, q^{T} & #theta_{max} cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPtMaxAngSelected = new TH1F("hPtMaxAngSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink}, q^{T} & #theta_{max} cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPdgCodemdRTPCclustersRejected = new TH2F("hPdgCodemdRTPCclustersRejected", "Mother's vs daughter's pdg code for rejected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink}, q^{T}, #theta_{max} & R/TPC clusters cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPtRTPCclustersSelected = new TH1F("hPtRTPCclustersSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink}, q^{T}, #theta_{max} & R/TPC clusters cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hRTPCclustersRTPCclustersSelected = new TH2F("hRTPCclustersRTPCclustersSelected", "Number of TPC clusters vs radius for selected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink}, q^{T}, #theta_{max} & R/TPC clusters cuts) R (cm); number of TPC clusters", nBinsR, hLowR, hUpR, 100, 0, 200); 
hPdgCodemdMinvRejected = new TH2F("hPdgCodemdMinvRejected", "Mother's vs daughter's pdg code for rejected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink}, q^{T}, #theta_{max}, R/TPC clusters & m_{inv} cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPtSelected = new TH1F("hPtSelected", "Mother's transverse momentum for selected kinks (all cuts except dE/dx); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaSelected = new TH1F("hEtaSelected", "Mother's pseudorapidity for selected kinks (all cuts except dE/dx); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapiditySelected = new TH1F("hRapiditySelected", "Mother's rapidity for selected kinks (all cuts except dE/dx); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtSelected = new TH1F("hQtSelected", "Daughter's transverse momentum in mother's frame for selected kinks (all cuts except dE/dx); q_{T} (GeV/c); dN/dq_{T}", 200, 0.0, 0.3);
hKinkAngleSelected = new TH1F("hKinkAngleSelected", "Kink angle of selected kinks (all cuts except dE/dx); #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hDCAkinkSelected = new TH1F("hDCAkinkSelected", "DCA; DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPmKinkAngSelected = new TH2F("hPmKinkAngSelected", "Mother's momentum magnitude vs kink angle for selected kinks (all cuts except dE/dx); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXYSelected = new TH2F("hKinkPosXYSelected", "X-Y Position ; X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZRSelected = new TH2F("hKinkPosZRSelected", "Z vs radius of selected kinks (all cuts except dE/dx); Z (cm); R (cm)", nBinsZ, hLowZ, hUpZ, nBinsR, hLowR, hUpR);
hKinkPosRSelected = new TH1F("hKinkPosRSelected", "Position radius of selected kinks (all cuts except dE/dx); R (cm); dN/dR", nBinsR, hLowR, hUpR);
hPdgCodemdSelected = new TH2F("hPdgCodemdSelected", "Mother vs daughter pdg code for selected kinks (all cuts except dE/dx); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPmdSelected = new TH2F("hPmdSelected", "Mother vs daughter momentum magnitude for selected kinks (all cuts except dE/dx); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuSelected = new TH1F("hMinvPimuSelected", "Invariant mass for #pi^{#pm} decaying to #mu^{#pm} for selected kinks (all cuts except dE/dx); m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hUIDKinkDauSelected = new TH1F("hUIDKinkDauSelected", "UID (method of production); UID (method of production); Number of particles", 100, 0.0, 20);
hdEdxSelected = new TH2F("hdEdxSelected", "dE/dx vs mother's momentum for selected kinks (all cuts except dE/dx); p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hPtSelectedFake = new TH1F("hPtSelectedFake", "Mother's transverse momentum for fake kinks (all cuts except dE/dx); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaSelectedFake = new TH1F("hEtaSelectedFake", "Mother's pseudorapidity for fake kinks (all cuts except dE/dx); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapiditySelectedFake = new TH1F("hRapiditySelectedFake", "Mother's rapidity ffor fake kinks (all cuts except dE/dx); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtSelectedFake = new TH1F("hQtSelectedFake", "Daughter's transverse momentum in mother's frame ffor fake kinks (all cuts except dE/dx); q_{T} (GeV/c); dN/dq_{T}", 200, 0.0, 0.3);
hKinkAngleSelectedFake = new TH1F("hKinkAngleSelectedFake", "Kink angle of fake kinks (all cuts except dE/dx); #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hDCAkinkSelectedFake = new TH1F("hDCAkinkSelectedFake", "Kink DCA for fake kinks (all cuts except dE/dx); DCA; Number of kinks", 100, 0.0, 2.0);
hPmKinkAngSelectedFake = new TH2F("hPmKinkAngSelectedFake", "Mother's momentum magnitude vs kink angle for fake kinks (all cuts except dE/dx);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXYSelectedFake = new TH2F("hKinkPosXYSelectedFake", "X-Y Position for fake kinks (all cuts except dE/dx); X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZRSelectedFake = new TH2F("hKinkPosZRSelectedFake", "Z vs radius of fake kinks (all cuts except dE/dx); Z (cm); R (cm)", nBinsZ, hLowZ, hUpZ, nBinsR, hLowR, hUpR);
hKinkPosRSelectedFake = new TH1F("hKinkPosRSelectedFake", "Position radius of fake kinks (all cuts except dE/dx); R (cm); dN/dR", nBinsR, hLowR, hUpR);
hPmdSelectedFake = new TH2F("hPmdSelectedFake", "Mother vs daughter momentum magnitude for fake kinks (all cuts except dE/dx); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuSelectedFake = new TH1F("hMinvPimuSelectedFake", "Invariant mass for #pi^{#pm} decaying to #mu^{#pm} for selected kinks (all cuts) for fake kinks (all cuts except dE/dx); m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hdEdxSelectedFake = new TH2F("hdEdxSelectedFake", "dE/dx vs mother's momentum for fake kinks (all cuts except dE/dx); p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hPdgCodemddEdxRejected = new TH2F("hPdgCodemddEdxRejected", "Mother's vs daughter's pdg code for rejected kinks (all cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPtPiSelected = new TH1F("hPtPiSelected", "Mother's transverse momentum for selected kinks (all cuts); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPiSelected = new TH1F("hEtaPiSelected", "Mother's pseudorapidity for selected kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityPiSelected = new TH1F("hRapidityPiSelected", "Mother's rapidity for selected kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtPiSelected = new TH1F("hQtPiSelected", "Daughter's transverse momentum in mother's frame for selected kinks (all cuts); q_{T} (GeV/c); dN/dq_{T}", 200, 0.0, 0.3);
hKinkAnglePiSelected = new TH1F("hKinkAnglePiSelected", "Kink angle of selected kinks (all cuts); #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hDCAkinkPiSelected = new TH1F("hDCAkinkPiSelected", "DCA; DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPmKinkAngPiSelected = new TH2F("hPmKinkAngPiSelected", "Mother's momentum magnitude vs kink angle for selected kinks (all cuts); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosRTPCclusters1PiSelected = new TH1F("hKinkPosRTPCclusters1PiSelected", " ;kinkposR; tpc clusters",100,0,1);
hKinkPosRTPCclusters2PiSelected = new TH2F("hKinkPosRTPCclusters2PiSelected", " ;kinkposR; tpc clusters",nBinsR, hLowR, hUpR,100,0,200 );
hKinkPosXYPiSelected = new TH2F("hKinkPosXYPiSelected", "X-Y Position ; X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZRPiSelected = new TH2F("hKinkPosZRPiSelected", "Z vs radius of selected kinks (all cuts); Z (cm); R (cm)", nBinsZ, hLowZ, hUpZ, nBinsR, hLowR, hUpR);
hKinkPosRPiSelected = new TH1F("hKinkPosRPiSelected", "Position radius of selected kinks (all cuts); R (cm); dN/dR", nBinsR, hLowR, hUpR);
hKinkPosZPiSelected = new TH1F("hKinkPosZPiSelected", "z position of selected kinks (all cuts); R (cm); dN/dR", nBinsZ, hLowZ, hUpZ);
hPmPPiSelected = new TH2F("hPmPPiSelected", "Mother's momentum as calculated by kink and by track of selected kinks (all cuts); R (cm); dN/dR", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hPdgCodemdPiSelected = new TH2F("hPdgCodemdPiSelected", "Mother vs daughter pdg code for selected kinks (all cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPmdPiSelected = new TH2F("hPmdPiSelected", "Mother vs daughter momentum magnitude for selected kinks (all cuts); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuPiSelected = new TH1F("hMinvPimuPiSelected", "Invariant mass for #pi^{#pm} decaying to #mu^{#pm} for selected kinks (all cuts); m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hUIDKinkDauPiSelected = new TH1F("hUIDKinkDauPiSelected", "UID (method of production); UID (method of production); Number of particles", 100, 0.0, 20);
hdEdxPiSelected = new TH2F("hdEdxPiSelected", "dE/dx vs mother's momentum for selected kinks (all cuts); p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);

hPtPiSelectedPlus = new TH1F("hPtPiSelectedPlus", "Mother's transverse momentum for selected positive kinks (all cuts); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPiSelectedPlus = new TH1F("hEtaPiSelectedPlus", "Mother's pseudorapidity for selected positive kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityPiSelectedPlus = new TH1F("hRapidityPiSelectedPlus", "Mother's rapidity for selected positive kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtPiSelectedPlus = new TH1F("hQtPiSelectedPlus", "Daughter's transverse momentum in mother's frame for selected positive kinks (all cuts); q_{T} (GeV/c); dN/dq_{T}", 200, 0.0, 0.3);
hKinkAnglePiSelectedPlus = new TH1F("hKinkAnglePiSelectedPlus", "Kink angle for selected positive kinks (all cuts); #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hDCAkinkPiSelectedPlus = new TH1F("hDCAkinkPiSelectedPlus", "DCA for selected positive kinks (all cuts); DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPmKinkAngPiSelectedPlus = new TH2F("hPmKinkAngPiSelectedPlus", "Mother's momentum magnitude vs kink angle for selected positive kinks (all cuts); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXYPiSelectedPlus = new TH2F("hKinkPosXYPiSelectedPlus", "X-Y Position of selected positive kinks (all cuts); X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZRPiSelectedPlus = new TH2F("hKinkPosZRPiSelectedPlus", "Z vs radius of selected positive kinks (all cuts); Z (cm); R (cm)", nBinsZ, hLowZ, hUpZ,nBinsR, hLowR, hUpR);
hPdgCodemdPiSelectedPlus = new TH2F("hPdgCodemdPiSelectedPlus", "Mother vs daughter pdg code for selected positive kinks (all cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPmdPiSelectedPlus = new TH2F("hPmdPiSelectedPlus", "Mother vs daughter momentum magnitude for selected positive kinks (all cuts); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuPiSelectedPlus = new TH1F("hMinvPimuPiSelectedPlus", "Invariant mass for #pi^{+} decaying to #mu^{+} for selected kinks (all cuts); m (GeV/c^{2}); dN/dm",nBinsInvMass, hLowInvMass, hUpInvMass);
hUIDPiDaumuSelectedPlus = new TH1F("hUIDPiDaumuSelectedPlus", "UID (method of production) of selected positive kinks (all cuts); UID (method of production); Number of particles", 100, 0.0, 20);
hdEdxPiSelectedPlus = new TH2F("hdEdxPiSelectedPlus", "dE/dx vs mother's momentum for selected positive kinks (all cuts); p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hPtPrimPiKinksPlus = new TH1F("hPtPrimPiKinksPlus", "Mother's transverse momentum for selected real #pi^{+} primary kinks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPrimPiKinksPlus = new TH1F("hEtaPrimPiKinksPlus", "Mother's pseudorapidity for selected real #pi^{+} primary kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityPrimPiKinksPlus = new TH1F("hRapidityPrimPiKinksPlus", "Mother's rapidity for selected real #pi^{+} primary kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hPtSecondPiKinksPlus = new TH1F("hPtSecondPiKinksPlus", "Mother's transverse momentum for selected real #pi^{+} secondary kinks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaSecondPiKinksPlus = new TH1F("hEtaSecondPiKinksPlus", "Mother's pseudorapidity for selected real #pi^{+} secondary kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapiditySecondPiKinksPlus = new TH1F("hRapiditySecondPiKinksPlus", "Mother's rapidity for selected real #pi^{+} secondary kinks; y; dN/dy", nBinsEta, hLowEta, hUpEta);
hPtNonPiKinksPlus = new TH1F("hPtNonPiKinksPlus", "Mother's transverse momentum for selected non #pi^{+} kinks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaNonPiKinksPlus = new TH1F("hEtaNonPiKinksPlus", "Mother's pseudorapidity for selected non #pi^{+} kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityNonPiKinksPlus = new TH1F("hRapidityNonPiKinksPlus", "Mother's rapidity for selected non #pi^{+} kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hPdgCodemdNonPiKinksPlus = new TH2F("hPdgCodemdNonPiKinksPlus", "Mother vs daughter pdg code for selected non #pi^{+} kinks; Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);

hPtPiSelectedMinus = new TH1F("hPtPiSelectedMinus", "Mother's transverse momentum for selected negative kinks (all cuts); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPiSelectedMinus = new TH1F("hEtaPiSelectedMinus", "Mother's pseudorapidity for selected negative kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityPiSelectedMinus = new TH1F("hRapidityPiSelectedMinus", "Mother's rapidity for selected negative kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtPiSelectedMinus = new TH1F("hQtPiSelectedMinus", "Daughter's transverse momentum in mother's frame for selected negative kinks (all cuts); q_{T} (GeV/c); dN/dq_{T}", 200, 0.0, 0.3);
hKinkAnglePiSelectedMinus = new TH1F("hKinkAnglePiSelectedMinus", "Kink angle for selected negative kinks (all cuts); #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hDCAkinkPiSelectedMinus = new TH1F("hDCAkinkPiSelectedMinus", "DCA for selected negative kinks (all cuts); DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPmKinkAngPiSelectedMinus = new TH2F("hPmKinkAngPiSelectedMinus", "Mother's momentum magnitude vs kink angle for selected negative kinks (all cuts); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXYPiSelectedMinus = new TH2F("hKinkPosXYPiSelectedMinus", "X-Y Position for selected negative kinks (all cuts); X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZRPiSelectedMinus = new TH2F("hKinkPosZRPiSelectedMinus", "Z vs radius for selected negative kinks (all cuts); Z (cm); R (cm)",nBinsZ, hLowZ, hUpZ,nBinsR, hLowR, hUpR);
hPdgCodemdPiSelectedMinus = new TH2F("hPdgCodemdPiSelectedMinus", "Mother vs daughter pdg code for selected negative kinks (all cuts); Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);
hPmdPiSelectedMinus = new TH2F("hPmdPiSelectedMinus", "Mother vs daughter momentum magnitude for selected negative kinks (all cuts); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuPiSelectedMinus = new TH1F("hMinvPimuPiSelectedMinus", "Invariant mass for #pi^{-} decaying to #mu^{-} for selected kinks (all cuts); m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hUIDPiDaumuSelectedMinus = new TH1F("hUIDPiDaumuSelectedMinus", "UID (method of production) of selected negative kinks (all cuts); UID (method of production); Number of particles", 100, 0.0, 20);
hdEdxPiSelectedMinus = new TH2F("hdEdxPiSelectedMinus", "dE/dx vs mother's momentum for selected negative kinks (all cuts); p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hPtPrimPiKinksMinus = new TH1F("hPtPrimPiKinksMinus", "Mother's transverse momentum for selected real #pi^{-} primary kinks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPrimPiKinksMinus = new TH1F("hEtaPrimPiKinksMinus", "Mother's pseudorapidity for selected real #pi^{-} primary kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityPrimPiKinksMinus = new TH1F("hRapidityPrimPiKinksMinus", "Mother's rapidity for selected real #pi^{-} primary kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hPtSecondPiKinksMinus = new TH1F("hPtSecondPiKinksMinus", "Mother's transverse momentum for selected real #pi^{-} secondary kinks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaSecondPiKinksMinus = new TH1F("hEtaSecondPiKinksMinus", "Mother's pseudorapidity for selected real #pi^{-} secondary kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapiditySecondPiKinksMinus = new TH1F("hRapiditySecondPiKinksMinus", "Mother's rapidity for selected real #pi^{-} secondary kinks; y; dN/dy", nBinsEta, hLowEta, hUpEta);
hPtNonPiKinksMinus = new TH1F("hPtNonPiKinksMinus", "Mother's transverse momentum for selected non #pi^{-} kinks; p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaNonPiKinksMinus = new TH1F("hEtaNonPiKinksMinus", "Mother's pseudorapidity for selected non #pi^{-} kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityNonPiKinksMinus = new TH1F("hRapidityNonPiKinksMinus", "Mother's rapidity for selected non #pi^{-} kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hPdgCodemdNonPiKinksMinus = new TH2F("hPdgCodemdNonPiKinksMinus", "Mother vs daughter pdg code for selected non #pi^{-} kinks; Mother's pdg code; Daughter's pdg code", nBinsPdg, hLowPdg, hUpPdg, nBinsPdg, hLowPdg, hUpPdg);



fListOfHistos->Add(hMCMult);
fListOfHistos->Add(hMCMultPrim); 
fListOfHistos->Add(hMCPtAll);
fListOfHistos->Add(hMCEtaAll);
fListOfHistos->Add(hMCPtPrim);
fListOfHistos->Add(hMCEtaPrim);
fListOfHistos->Add(hMCPt);
fListOfHistos->Add(hMCEta);
fListOfHistos->Add(hMCPdg);
fListOfHistos->Add(hMCMultPiPlus);
fListOfHistos->Add(hMCPtPiPlus);
fListOfHistos->Add(hMCEtaPiPlus);
fListOfHistos->Add(hMCRapidityPiPlus);
fListOfHistos->Add(hMCNDaughtersPlus);
fListOfHistos->Add(hMCRadPiDauPlus);
fListOfHistos->Add(hMCKinkPosZPlus);
fListOfHistos->Add(hMCUIDPiDauPlus);
fListOfHistos->Add(hMCPdgPiNonDecayedPlus);
fListOfHistos->Add(hMCPdgPiDauPlus);
fListOfHistos->Add(hMCQtPlus);
fListOfHistos->Add(hMCKinkAnglePlus);
fListOfHistos->Add(hMCPKinkAngPlus);
fListOfHistos->Add(hMCPdgCodemdPlus);
fListOfHistos->Add(hMCPtmdPlus);
fListOfHistos->Add(hMCPtPimuonPlus);
fListOfHistos->Add(hMCEtaPimuonPlus);
fListOfHistos->Add(hMCRapidityPimuonPlus);
fListOfHistos->Add(hMCQtPimuonPlus);
fListOfHistos->Add(hMCPKinkAngPimuonPlus);
fListOfHistos->Add(hMCPtPiotherPlus);
fListOfHistos->Add(hMCEtaPiotherPlus);
fListOfHistos->Add(hMCRapidityPiotherPlus);
fListOfHistos->Add(hMCQtPiotherPlus);
fListOfHistos->Add(hMCPKinkAngPiotherPlus);
fListOfHistos->Add(hMCMultPiMinus);
fListOfHistos->Add(hMCPtPiMinus);
fListOfHistos->Add(hMCEtaPiMinus);
fListOfHistos->Add(hMCRapidityPiMinus);
fListOfHistos->Add(hMCNDaughtersMinus);
fListOfHistos->Add(hMCRadPiDauMinus);
fListOfHistos->Add(hMCKinkPosZMinus);
fListOfHistos->Add(hMCUIDPiDauMinus);
fListOfHistos->Add(hMCPdgPiNonDecayedMinus);
fListOfHistos->Add(hMCPdgPiDauMinus);
fListOfHistos->Add(hMCQtMinus);
fListOfHistos->Add(hMCKinkAngleMinus);
fListOfHistos->Add(hMCPKinkAngMinus);
fListOfHistos->Add(hMCPdgCodemdMinus);
fListOfHistos->Add(hMCPtmdMinus);
fListOfHistos->Add(hMCPtPimuonMinus);
fListOfHistos->Add(hMCEtaPimuonMinus);
fListOfHistos->Add(hMCRapidityPimuonMinus);
fListOfHistos->Add(hMCQtPimuonMinus);
fListOfHistos->Add(hMCPKinkAngPimuonMinus);
fListOfHistos->Add(hMCPtPiotherMinus);
fListOfHistos->Add(hMCEtaPiotherMinus);
fListOfHistos->Add(hMCRapidityPiotherMinus);
fListOfHistos->Add(hMCQtPiotherMinus);
fListOfHistos->Add(hMCPKinkAngPiotherMinus);

fListOfHistos->Add(hMult);
fListOfHistos->Add(hAcceptedMult);
fListOfHistos->Add(hMultPS);
fListOfHistos->Add(hvtx);
fListOfHistos->Add(hvtxy);
fListOfHistos->Add(hvtyz);
fListOfHistos->Add(hvtxz);
fListOfHistos->Add(hMultPSV);
fListOfHistos->Add(hPtAll);
fListOfHistos->Add(hEtaAll);
fListOfHistos->Add(hTrackPos);
fListOfHistos->Add(hTrackPosxy);
fListOfHistos->Add(hTrackPosyz);
fListOfHistos->Add(hTrackPosxz);
//fListOfHistos->Add(hTPCchi2clusters);
//fListOfHistos->Add(hdcaToVertexXY);
//fListOfHistos->Add(hdcaToVertexZ);
fListOfHistos->Add(hMultPrim);
fListOfHistos->Add(hPtPrim);
fListOfHistos->Add(hEtaPrim);
fListOfHistos->Add(hPrimTrackPos);
fListOfHistos->Add(hPrimTrackPosxy);
fListOfHistos->Add(hPrimTrackPosyz);
fListOfHistos->Add(hPrimTrackPosxz);
fListOfHistos->Add(hPt);
fListOfHistos->Add(hEta);
//fListOfHistos->Add(hRapidity);
fListOfHistos->Add(hPtKink);
fListOfHistos->Add(hEtaKink);
fListOfHistos->Add(hRapidityKink);
fListOfHistos->Add(hPmP);
fListOfHistos->Add(hKinkPosRTPCclusters1);
fListOfHistos->Add(hKinkPosRTPCclusters2);
fListOfHistos->Add(hQt);
fListOfHistos->Add(hKinkAngle);
fListOfHistos->Add(hDCAkink);
fListOfHistos->Add(hPmKinkAng);
fListOfHistos->Add(hKinkPosXY);
fListOfHistos->Add(hKinkPosZY);
fListOfHistos->Add(hKinkPosZR);
fListOfHistos->Add(hKinkPosR);
fListOfHistos->Add(hKinkPosZ);
fListOfHistos->Add(hKinkPosZMCKinkPosZ);
fListOfHistos->Add(hPdgCodemd);
fListOfHistos->Add(hPmd);
fListOfHistos->Add(hMinvPimu);
fListOfHistos->Add(hUIDKinkDau);
fListOfHistos->Add(hdEdx);
fListOfHistos->Add(hPtKinkFake);
fListOfHistos->Add(hEtaKinkFake);
fListOfHistos->Add(hRapidityKinkFake);
fListOfHistos->Add(hPmPFake);
fListOfHistos->Add(hKinkPosRTPCclusters1Fake);
fListOfHistos->Add(hKinkPosRTPCclusters2Fake);
fListOfHistos->Add(hQtFake);
fListOfHistos->Add(hKinkAngleFake);
fListOfHistos->Add(hDCAkinkFake);
fListOfHistos->Add(hPmKinkAngFake);
fListOfHistos->Add(hKinkPosXYFake);
fListOfHistos->Add(hKinkPosZYFake);
fListOfHistos->Add(hKinkPosZRFake);
fListOfHistos->Add(hKinkPosRFake);
fListOfHistos->Add(hKinkPosZFake);
fListOfHistos->Add(hKinkPosZMCKinkPosZFake);
fListOfHistos->Add(hPdgCodemdFake);
fListOfHistos->Add(hPmdFake);
fListOfHistos->Add(hMinvPimuFake);
fListOfHistos->Add(hUIDKinkDauFake);
fListOfHistos->Add(hdEdxFake);
fListOfHistos->Add(hPtPosRSelected);
fListOfHistos->Add(hPdgCodemdZRejected);
fListOfHistos->Add(hPtZSelected);
fListOfHistos->Add(hPdgCodemdAngRejected);
fListOfHistos->Add(hPtAngSelected);
fListOfHistos->Add(hPdgCodemdPmRejected);
fListOfHistos->Add(hPtPmSelected);
fListOfHistos->Add(hPdgCodemdQtLowRejected);
fListOfHistos->Add(hPtGoodKink);
fListOfHistos->Add(hEtaGoodKink);
fListOfHistos->Add(hRapidityGoodKink);
fListOfHistos->Add(hQtGoodKink);
fListOfHistos->Add(hPmGoodKinkAng);
fListOfHistos->Add(hPdgCodemdGoodKink);
fListOfHistos->Add(hPmdGoodKink);
fListOfHistos->Add(hUIDGoodKinkDau);
fListOfHistos->Add(hdEdxGoodKink);
fListOfHistos->Add(hUIDPiDauPlus);
fListOfHistos->Add(hMultPiPlus);
fListOfHistos->Add(hPtPiPlus);
fListOfHistos->Add(hEtaPiPlus);
fListOfHistos->Add(hRapidityPiPlus);
fListOfHistos->Add(hQtPiPlus);
fListOfHistos->Add(hKinkAnglePiPlus);
fListOfHistos->Add(hPmKinkAngPiPlus);
fListOfHistos->Add(hKinkPosXYPiPlus);
fListOfHistos->Add(hKinkPosZRPiPlus);
fListOfHistos->Add(hKinkPosRPiPlus);
fListOfHistos->Add(hDCAkinkPiPlus);
fListOfHistos->Add(hPdgCodemdPiPlus);
fListOfHistos->Add(hPmdPiPlus);
fListOfHistos->Add(hdEdxPiPlus);
fListOfHistos->Add(hQtPimuPlus);
fListOfHistos->Add(hPmKinkAngPimuPlus);
fListOfHistos->Add(hQtPiotherPlus);
fListOfHistos->Add(hPmKinkAngPiotherPlus);
fListOfHistos->Add(hPdgCodemdPiotherPlus);
fListOfHistos->Add(hUIDPiDauMinus);
fListOfHistos->Add(hMultPiMinus);
fListOfHistos->Add(hPtPiMinus);
fListOfHistos->Add(hEtaPiMinus);
fListOfHistos->Add(hRapidityPiMinus);
fListOfHistos->Add(hQtPiMinus);
fListOfHistos->Add(hKinkAnglePiMinus);
fListOfHistos->Add(hPmKinkAngPiMinus);
fListOfHistos->Add(hKinkPosXYPiMinus);
fListOfHistos->Add(hKinkPosZRPiMinus);
fListOfHistos->Add(hKinkPosRPiMinus);
fListOfHistos->Add(hDCAkinkPiMinus);
fListOfHistos->Add(hPdgCodemdPiMinus);
fListOfHistos->Add(hPmdPiMinus);
fListOfHistos->Add(hdEdxPiMinus);
fListOfHistos->Add(hQtPimuMinus);
fListOfHistos->Add(hPmKinkAngPimuMinus);
fListOfHistos->Add(hQtPiotherMinus);
fListOfHistos->Add(hPmKinkAngPiotherMinus);
fListOfHistos->Add(hPdgCodemdPiotherMinus);
fListOfHistos->Add(hPdgCodemdQtRejected);
fListOfHistos->Add(hPtQtSelected);
fListOfHistos->Add(hPdgCodemdMaxAngRejected);
fListOfHistos->Add(hPtMaxAngSelected);
fListOfHistos->Add(hPdgCodemdRTPCclustersRejected);
fListOfHistos->Add(hPtRTPCclustersSelected);
fListOfHistos->Add(hRTPCclustersRTPCclustersSelected);
fListOfHistos->Add(hPdgCodemdMinvRejected);
fListOfHistos->Add(hPtSelected);
fListOfHistos->Add(hEtaSelected);
fListOfHistos->Add(hRapiditySelected);
fListOfHistos->Add(hQtSelected);
fListOfHistos->Add(hKinkAngleSelected);
fListOfHistos->Add(hDCAkinkSelected);
fListOfHistos->Add(hPmKinkAngSelected);
fListOfHistos->Add(hKinkPosXYSelected);
fListOfHistos->Add(hKinkPosZRSelected);
fListOfHistos->Add(hKinkPosRSelected);
fListOfHistos->Add(hPdgCodemdSelected);
fListOfHistos->Add(hPmdSelected);
fListOfHistos->Add(hMinvPimuSelected);
fListOfHistos->Add(hUIDKinkDauSelected);
fListOfHistos->Add(hdEdxSelected);
fListOfHistos->Add(hPtSelectedFake);
fListOfHistos->Add(hEtaSelectedFake);
fListOfHistos->Add(hRapiditySelectedFake);
fListOfHistos->Add(hQtSelectedFake);
fListOfHistos->Add(hKinkAngleSelectedFake);
fListOfHistos->Add(hDCAkinkSelectedFake);
fListOfHistos->Add(hPmKinkAngSelectedFake);
fListOfHistos->Add(hKinkPosXYSelectedFake);
fListOfHistos->Add(hKinkPosZRSelectedFake);
fListOfHistos->Add(hKinkPosRSelectedFake);
fListOfHistos->Add(hPmdSelectedFake);
fListOfHistos->Add(hMinvPimuSelectedFake);
fListOfHistos->Add(hdEdxSelectedFake);
fListOfHistos->Add(hPdgCodemddEdxRejected);
fListOfHistos->Add(hPtPiSelected);
fListOfHistos->Add(hEtaPiSelected);
fListOfHistos->Add(hRapidityPiSelected);
fListOfHistos->Add(hQtPiSelected);
fListOfHistos->Add(hKinkAnglePiSelected);
fListOfHistos->Add(hDCAkinkPiSelected);
fListOfHistos->Add(hPmKinkAngPiSelected);
fListOfHistos->Add(hKinkPosRTPCclusters1PiSelected);
fListOfHistos->Add(hKinkPosRTPCclusters2PiSelected);
fListOfHistos->Add(hKinkPosXYPiSelected);
fListOfHistos->Add(hKinkPosZRPiSelected);
fListOfHistos->Add(hKinkPosRPiSelected);
fListOfHistos->Add(hKinkPosZPiSelected);
fListOfHistos->Add(hPmPPiSelected);
fListOfHistos->Add(hPdgCodemdPiSelected);
fListOfHistos->Add(hPmdPiSelected);
fListOfHistos->Add(hMinvPimuPiSelected);
fListOfHistos->Add(hUIDKinkDauPiSelected);
fListOfHistos->Add(hdEdxPiSelected);

fListOfHistos->Add(hPtPiSelectedPlus);
fListOfHistos->Add(hEtaPiSelectedPlus);
fListOfHistos->Add(hRapidityPiSelectedPlus);
fListOfHistos->Add(hQtPiSelectedPlus);
fListOfHistos->Add(hKinkAnglePiSelectedPlus);
fListOfHistos->Add(hDCAkinkPiSelectedPlus);
fListOfHistos->Add(hPmKinkAngPiSelectedPlus);
fListOfHistos->Add(hKinkPosXYPiSelectedPlus);
fListOfHistos->Add(hKinkPosZRPiSelectedPlus);
fListOfHistos->Add(hPdgCodemdPiSelectedPlus);
fListOfHistos->Add(hPmdPiSelectedPlus);
fListOfHistos->Add(hMinvPimuPiSelectedPlus);
fListOfHistos->Add(hUIDPiDaumuSelectedPlus);
fListOfHistos->Add(hdEdxPiSelectedPlus);
fListOfHistos->Add(hPtPrimPiKinksPlus);
fListOfHistos->Add(hEtaPrimPiKinksPlus);
fListOfHistos->Add(hRapidityPrimPiKinksPlus);
fListOfHistos->Add(hPtSecondPiKinksPlus);
fListOfHistos->Add(hEtaSecondPiKinksPlus);
fListOfHistos->Add(hRapiditySecondPiKinksPlus);
fListOfHistos->Add(hPtNonPiKinksPlus);
fListOfHistos->Add(hEtaNonPiKinksPlus);
fListOfHistos->Add(hRapidityNonPiKinksPlus);
fListOfHistos->Add(hPdgCodemdNonPiKinksPlus);

fListOfHistos->Add(hPtPiSelectedMinus);
fListOfHistos->Add(hEtaPiSelectedMinus);
fListOfHistos->Add(hRapidityPiSelectedMinus);
fListOfHistos->Add(hQtPiSelectedMinus);
fListOfHistos->Add(hKinkAnglePiSelectedMinus);
fListOfHistos->Add(hDCAkinkPiSelectedMinus);
fListOfHistos->Add(hPmKinkAngPiSelectedMinus);
fListOfHistos->Add(hKinkPosXYPiSelectedMinus);
fListOfHistos->Add(hKinkPosZRPiSelectedMinus);
fListOfHistos->Add(hPdgCodemdPiSelectedMinus);
fListOfHistos->Add(hPmdPiSelectedMinus);
fListOfHistos->Add(hMinvPimuPiSelectedMinus);
fListOfHistos->Add(hUIDPiDaumuSelectedMinus);
fListOfHistos->Add(hdEdxPiSelectedMinus);
fListOfHistos->Add(hPtPrimPiKinksMinus);
fListOfHistos->Add(hEtaPrimPiKinksMinus);
fListOfHistos->Add(hRapidityPrimPiKinksMinus);
fListOfHistos->Add(hPtSecondPiKinksMinus);
fListOfHistos->Add(hEtaSecondPiKinksMinus);
fListOfHistos->Add(hRapiditySecondPiKinksMinus);
fListOfHistos->Add(hPtNonPiKinksMinus);
fListOfHistos->Add(hEtaNonPiKinksMinus);
fListOfHistos->Add(hRapidityNonPiKinksMinus);
fListOfHistos->Add(hPdgCodemdNonPiKinksMinus);

fListOfHistos->SetOwner(kTRUE);
PostData(1, fListOfHistos);
}

//________________________________________________________________________
void AliAnalysisPionKinksMCESD::UserExec(Option_t *) {
	AliVEvent *event = InputEvent();
	if (!event) {
		Printf("ERROR: Could not retrieve event");
		return;
	}

	AliMCEvent* mcEvent = MCEvent();
	if (!mcEvent) {
		Printf("ERROR: Could not retrieve MC event");
		return;
	}
	AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(event);
	if (!esdEvent) {
		Printf("ERROR: Could not retrieve esd");
		return;
	}

 
//------------------------------ MC Multiplicity unbiased------------------------//
	AliStack* MCstack = mcEvent->Stack(); //stack of MC events
	Int_t MCNTracks = mcEvent->GetNumberOfTracks(); 
//	Int_t MCNTracks2 = MCstack->GetNtrack(); 
 	hMCMult->Fill(MCNTracks);
	Int_t MCNPrimTracks = MCstack->GetNprimary(); 
	//Printf("There are %d number of primary tracks in this Monte Carlo event", MCNPrimTracks);
	hMCMultPrim->Fill(MCNPrimTracks); 
	
//------------------------ Data Multiplicity unbiased --------------------------//
	Int_t NTracks = esdEvent->GetNumberOfTracks(); //number of ESD tracks 
	hMult->Fill(NTracks);
		
//----------------------------- Accepted Multiplicity -----------------------------//
	Float_t NAcceptedTracks = fTrackCuts->CountAcceptedTracks(esdEvent); 
	if(fLowMulcut>-1) {
		if(NAcceptedTracks<fLowMulcut) return;
	}
	if(fUpMulcut>-1) {
		if(NAcceptedTracks>fUpMulcut) return;
	}
	hAcceptedMult->Fill(NAcceptedTracks); //to check if the multiplicity limits are ok


//-------------------------------- MC data Analysis -----------------------------//
 	for (Int_t iMC = 0; iMC < MCNTracks; iMC++) { //loop on all accepted MC tracks
		TParticle* MCparticle = MCstack->Particle(iMC);
		if (!MCparticle) {
			Printf("ERROR: Could not receive MC particle %d", iMC);
			continue;
		}

		Double_t MCPt = MCparticle->Pt(); 
		Double_t MCP = MCparticle->P(); 
		Double_t MCEta = MCparticle->Eta(); 
		hMCPtAll->Fill(MCPt); 
		hMCEtaAll->Fill(MCEta); 
		
		if (!MCstack->IsPhysicalPrimary(iMC)) continue; 
		hMCPtPrim->Fill(MCPt);
		hMCEtaPrim->Fill(MCEta); 

		if (MCPt<cLowPt) continue;

		hMCPt->Fill(MCPt);
		hMCEta->Fill(MCEta); 

		Int_t MCPdg = MCparticle->GetPdgCode(); 
		hMCPdg->Fill(MCPdg); 

		if (MCPdg==cPdgPion) { //positive pion selection
		Double_t MCRapidity = GetMCRapidity(MCparticle);	
	 	if (TMath::Abs(MCRapidity)>cRapidityLim) continue;
			hMCMultPiPlus->Fill(iMC); //MC pion multiplicity	
			hMCPtPiPlus->Fill(MCPt); 
			hMCEtaPiPlus->Fill(MCEta); 
			hMCRapidityPiPlus->Fill(MCRapidity);

			Int_t MCNDaughtersPlus = MCparticle->GetNDaughters();
			hMCNDaughtersPlus->Fill(MCNDaughtersPlus);
			Int_t FirstDPlus = MCparticle->GetFirstDaughter(); //first daughter's label
			Int_t LastDPlus = MCparticle->GetLastDaughter(); //last daughter's label
			if ((FirstDPlus>MCNTracks)||(LastDPlus>MCNTracks)) continue;


			for (Int_t iMCdPlus=FirstDPlus; iMCdPlus<=LastDPlus; iMCdPlus++) { //loop on pion daughters
				if (iMCdPlus<0) continue; //debug
				TParticle* MCdaughterPlus = MCstack->Particle(iMCdPlus); 

				Double_t RdPlus = MCdaughterPlus->R(); //position radius of daughter's generation vertex 
				hMCRadPiDauPlus->Fill(RdPlus); 

				Double_t MCKinkPosZPlus = MCdaughterPlus->Vz();
				hMCKinkPosZPlus->Fill(MCKinkPosZPlus);
				
				if ((TMath::Abs(MCKinkPosZPlus)<cLowZ) || (TMath::Abs(MCKinkPosZPlus)>cUpZ)) continue;

				if ((RdPlus<cLowR) || (RdPlus>cUpR)) continue; //selection of daughters that where generated in TPC
	
				Int_t dcodePlus = MCdaughterPlus->GetPdgCode(); 
						
				UInt_t MCDauUIDPlus = MCdaughterPlus->GetUniqueID(); //daughter's unique id (= method of production) 
				hMCUIDPiDauPlus->Fill(MCDauUIDPlus); 

				if (MCDauUIDPlus!=4) { //selection of decayed daughters
					hMCPdgPiNonDecayedPlus->Fill(dcodePlus); 
					continue;
				}

				hMCPdgPiDauPlus->Fill(TMath::Abs(dcodePlus));

				Double_t MCPtdPlus = MCdaughterPlus->Pt(); //daughter's transverse momentum in lab frame
				Double_t MCQtPlus = MCPQt(mcEvent, iMC, MCdaughterPlus); //daughter's transverse momentum in mother's frame (Qt)
				Double_t MCKinkAnglePlus2 = fuMCKinkAngle(mcEvent, iMC, MCdaughterPlus, kTRUE); //kink angle in degrees
				Double_t MCMaxKinkAngPimuPlus=fMaxKinkAngPimu->Eval(MCP,0.,0.,0.); //maximum decay angle in lab for pion decaying to muon
				
				// if  (MCKinkAnglePlus2>(MCMaxKinkAngPimuPlus*1.1)) continue;
				if  (MCKinkAnglePlus2>(MCMaxKinkAngPimuPlus)) continue;
			
				hMCQtPlus->Fill(MCQtPlus); 
				hMCKinkAnglePlus->Fill(MCKinkAnglePlus2);
				hMCPKinkAngPlus->Fill(MCP,MCKinkAnglePlus2);
				hMCPdgCodemdPlus->Fill(TMath::Abs(MCPdg),TMath::Abs(dcodePlus)); 
				hMCPtmdPlus->Fill(MCPt, MCPtdPlus); 	

				if ((MCQtPlus<cLowQt) || (MCQtPlus>cUpQt) || (MCKinkAnglePlus2<cLowKinkAngle)) continue; 

				if (dcodePlus==-cPdgMuon) { //muon daughters selection
					hMCPtPimuonPlus->Fill(MCPt);
					hMCEtaPimuonPlus->Fill(MCEta);
					hMCRapidityPimuonPlus->Fill(MCRapidity);
					hMCQtPimuonPlus->Fill(MCQtPlus); 
					hMCPKinkAngPimuonPlus->Fill(MCKinkAnglePlus2); 
				} else {hMCPtPiotherPlus->Fill(MCPt);
					hMCEtaPiotherPlus->Fill(MCEta);
					hMCRapidityPiotherPlus->Fill(MCRapidity);
					hMCQtPiotherPlus->Fill(MCQtPlus); 
					hMCPKinkAngPiotherPlus->Fill(MCKinkAnglePlus2);
				}
			} //end of pion daughters' loop
		} else if (MCPdg==-cPdgPion) { //negative pion selection
		Double_t MCRapidity = GetMCRapidity(MCparticle);	
	 	if (TMath::Abs(MCRapidity)>cRapidityLim) continue;
			hMCMultPiMinus->Fill(iMC); //MC pion multiplicity	
			hMCPtPiMinus->Fill(MCPt); 
			hMCEtaPiMinus->Fill(MCEta); 
			hMCRapidityPiMinus->Fill(MCRapidity);

			Int_t MCNDaughtersMinus = MCparticle->GetNDaughters();
			hMCNDaughtersMinus->Fill(MCNDaughtersMinus);
			Int_t FirstDMinus = MCparticle->GetFirstDaughter(); //first daughter's label
			Int_t LastDMinus = MCparticle->GetLastDaughter(); //last daughter's label
			if ((FirstDMinus>MCNTracks)||(LastDMinus>MCNTracks)) continue;

			for (Int_t iMCdMinus=FirstDMinus; iMCdMinus<=LastDMinus; iMCdMinus++) { //loop on pion daughters
				if (iMCdMinus<0) continue; //debug
				TParticle* MCdaughterMinus = MCstack->Particle(iMCdMinus); 

				Double_t RdMinus = MCdaughterMinus->R(); //position radius of daughter's generation vertex 
				hMCRadPiDauMinus->Fill(RdMinus); 

				Double_t MCKinkPosZMinus = MCdaughterMinus->Vz();
				hMCKinkPosZMinus->Fill(MCKinkPosZMinus);
				
				if ((TMath::Abs(MCKinkPosZMinus)<cLowZ) || (TMath::Abs(MCKinkPosZMinus)>cUpZ)) continue;

				if ((RdMinus<cLowR) || (RdMinus>cUpR)) continue; //selection of daughters that where generated in TPC
	
				Int_t dcodeMinus = MCdaughterMinus->GetPdgCode(); 
						
				UInt_t MCDauUIDMinus = MCdaughterMinus->GetUniqueID(); //daughter's unique id (= method of production) 
				hMCUIDPiDauMinus->Fill(MCDauUIDMinus); 

				if (MCDauUIDMinus!=4) { //selection of decayed daughters
					hMCPdgPiNonDecayedMinus->Fill(dcodeMinus); 
					continue;
				}

				hMCPdgPiDauMinus->Fill(TMath::Abs(dcodeMinus));

 
				Double_t MCPtdMinus = MCdaughterMinus->Pt(); //daughter's transverse momentum in lab frame
				Double_t MCQtMinus = MCPQt(mcEvent, iMC, MCdaughterMinus); //daughter's transverse momentum in mother's frame (Qt)
				Double_t MCKinkAngleMinus2 = fuMCKinkAngle(mcEvent, iMC, MCdaughterMinus, kTRUE); //kink angle in degrees
				Double_t MCMaxKinkAngPimuMinus=fMaxKinkAngPimu->Eval(MCP,0.,0.,0.); //maximum decay angle in lab for pion decaying to muon
				
				if  (MCKinkAngleMinus2>(MCMaxKinkAngPimuMinus*1.1)) continue;
			
				hMCQtMinus->Fill(MCQtMinus); 
				hMCKinkAngleMinus->Fill(MCKinkAngleMinus2);
				hMCPKinkAngMinus->Fill(MCP,MCKinkAngleMinus2);
				hMCPdgCodemdMinus->Fill(TMath::Abs(MCPdg),TMath::Abs(dcodeMinus)); 
				hMCPtmdMinus->Fill(MCPt, MCPtdMinus); 	

				if ((MCQtMinus<cLowQt) || (MCQtMinus>cUpQt) || (MCKinkAngleMinus2<cLowKinkAngle)) continue; 

				if (dcodeMinus==cPdgMuon) { //muon daughters selection
					hMCPtPimuonMinus->Fill(MCPt);
					hMCEtaPimuonMinus->Fill(MCEta);
					hMCRapidityPimuonMinus->Fill(MCRapidity);
					hMCQtPimuonMinus->Fill(MCQtMinus); 
					hMCPKinkAngPimuonMinus->Fill(MCKinkAngleMinus2); 
				} else {hMCPtPiotherMinus->Fill(MCPt);///////////////////////////
					hMCEtaPiotherMinus->Fill(MCEta);
					hMCRapidityPiotherMinus->Fill(MCRapidity);
					hMCQtPiotherMinus->Fill(MCQtMinus); 
					hMCPKinkAngPiotherMinus->Fill(MCKinkAngleMinus2);
				}
			} //end of pion daughters' loop
		}

	} // end of MC tracks loop



//----------------------------- Physics selection ------------------------------//
	Bool_t IsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB;
	if ( IsSelected ==kFALSE) return;
	
//--------------- Data Multiplicity after Physics selection --------------------//
	hMultPS->Fill(NTracks);
	
//------------------------------------ Vertex ----------------------------------//
	const AliESDVertex* vtx = GetEventVertex(esdEvent); //ESD primary vertex
	if (!vtx) return;

	Double_t vtxpos[3]; //vertex position vector
	vtx->GetXYZ(vtxpos); 
	hvtx->Fill(vtxpos[0], vtxpos[1], vtxpos[2]); //ESD primary vertex position (x-y-z)
	hvtxy->Fill(vtxpos[0], vtxpos[1]); 
	hvtyz->Fill(vtxpos[1], vtxpos[2]); 
	hvtxz->Fill(vtxpos[0], vtxpos[2]); 

	if (TMath::Abs(vtxpos[2])>10.) return;

//-------------------- Data Multiplicity after vertex cut -----------------------//
	hMultPSV->Fill(NTracks);

//------------------------------- dE/dx parameters -----------------------------//
  if(!fPIDResponse) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
  } 
	
//------------------------- Reconstructed data Analysis -----------------------//
	for (Int_t iData = 0; iData<NTracks; iData++) { //loop on all ESD tracks
		AliESDtrack *ESDtrack = esdEvent->GetTrack(iData);
		if (!ESDtrack) {
			Printf("ERROR: Could not receive ESD track %d", iData);
			continue;
		}

		Double_t Pt = ESDtrack->Pt(); 
		Double_t Eta = ESDtrack->Eta(); 

		hPtAll->Fill(Pt); 
		hEtaAll->Fill(Eta); 
		//hRapidityAll->Fill(Rapidity(ESDtrack));

		Double_t TrackPos[3]; //starting position of ESD tracks (x-y-z)
		ESDtrack->GetXYZ(TrackPos);
		hTrackPos->Fill(TrackPos[0],TrackPos[1],TrackPos[2]); 
		hTrackPosxy->Fill(TrackPos[0], TrackPos[1]); 
		hTrackPosyz->Fill(TrackPos[1], TrackPos[2]); 
		hTrackPosxz->Fill(TrackPos[0], TrackPos[2]); 

		if ((IsPrimaryTrack(ESDtrack) == kFALSE) || (IsGoodTrack(ESDtrack) == kFALSE)) continue; //reject bad & secondary tracks
		hMultPrim->Fill(NTracks); 
 		hPtPrim->Fill(Pt);
		hEtaPrim->Fill(Eta); 
		//hRapidityPrim->Fill(Rapidity(ESDtrack));

		hPrimTrackPos->Fill(TrackPos[0],TrackPos[1],TrackPos[2]); //starting position of primary - supposed ESD tracks (x-y-z)
		hPrimTrackPosxy->Fill(TrackPos[0], TrackPos[1]); 
		hPrimTrackPosyz->Fill(TrackPos[1], TrackPos[2]); 
		hPrimTrackPosxz->Fill(TrackPos[0], TrackPos[2]); 

		if (Pt<cLowPt) continue;

		hPt->Fill(Pt); 
		hEta->Fill(Eta); 
		//hRapidity->Fill(fuRapidity(ESDtrack));

		Int_t KinkIndex = ESDtrack->GetKinkIndex(0); //kink index (1st component is negative if the track is a kink candidate)
		if (KinkIndex>=0) continue; //kink selection 

		Double_t Rapidity = fuRapidity(ESDtrack);
		if (TMath::Abs(Rapidity)>cRapidityLim) continue;
		
		Double_t dEdx = ESDtrack->GetTPCsignal();
		Double_t NTPCclusters = ESDtrack->GetTPCclusters(0);

		AliESDkink *kink = esdEvent->GetKink(TMath::Abs(KinkIndex)-1);
		
		Int_t Label1 = kink->GetLabel(0); //mother's track MC label (first component is mother)
		Int_t Label2 = kink->GetLabel(1); //daughter's track MC label (second component is daughter)
		if (Label1>MCNTracks) continue;
		if (Label2>MCNTracks) continue;
		
		const TVector3 KinkPos(kink->GetPosition());
		Double_t KinkPosR = kink->GetR(); //kink's position radius 
		Double_t KinkDistance = kink->GetDistance();
		Double_t Qt = kink->GetQt(); //daughter's transverse momentum in mother's frame
		Double_t KinkAngle =TMath::RadToDeg()*kink->GetAngle(2); //kink angle in mother frame (3rd component is angle in rads)

		const TVector3 PMother(kink->GetMotherP()); //ESD mother's momentum
		Double_t Pmx = PMother.Px();
		Double_t Pmy = PMother.Py();
		Double_t Pmz = PMother.Pz();
		Double_t Pm = PMother.Mag(); //ESD mother's momentum magnitude
		Double_t PTrack[3];
		ESDtrack->GetPxPyPz(PTrack);
		TVector3 Pvector3(PTrack[0], PTrack[1], PTrack[2]);
		Double_t P = Pvector3.Mag();
		Double_t PTPC = ESDtrack->GetInnerParam()->GetP();
		
		const TVector3 PDaughter(kink->GetDaughterP()); //ESD daughter's momentum
		Double_t Pdx = PDaughter.Px();
		Double_t Pdy = PDaughter.Py();
		Double_t Pdz = PDaughter.Pz();
		Double_t Pd = PDaughter.Mag(); //ESD daugter's momentum magnitude
		Double_t Edmu = TMath::Sqrt(TMath::Power(Pd,2)+TMath::Power(cMuonMass,2)); //ESD muon daughter's energy

		Double_t DP = TMath::Sqrt((Pmx-Pdx)*(Pmx-Pdx)+(Pmy-Pdy)*(Pmy-Pdy)+(Pmz-Pdz)*(Pmz-Pdz)); //transferred momentum magnitude
		Double_t MinvPimu = TMath::Sqrt((Edmu+DP)*(Edmu+DP)-TMath::Power(Pm,2)); //pion mother's invariant mass when decaying to muon
		Double_t MaxKinkAngPimu=fMaxKinkAngPimu->Eval(P,0.,0.,0.); //maximum decay angle in lab frame for a pion decaying to muon

		TParticle *particle1 = MCstack->Particle(TMath::Abs(kink->GetLabel(0))); //mother MC particle object
		TParticle *particle2 = MCstack->Particle(TMath::Abs(kink->GetLabel(1))); //daughter MC particle object
		
		Int_t code1 = particle1->GetPdgCode(); //mother's pdg code obtained from MC mother object
		Int_t code2 = particle2->GetPdgCode(); //daughter's pdg code obtained from MC daughter object

		UInt_t DauUID = particle2->GetUniqueID(); //ESD daughter's unique id (method of production)

		Double_t MCKinkPosZ = particle2->Vz();

		hPtKink->Fill(Pt); 
		hEtaKink->Fill(Eta); 
		hRapidityKink->Fill(Rapidity);
		hPmP->Fill(Pm,P);
		hKinkPosRTPCclusters1->Fill(NTPCclusters/KinkPosR); 
		hKinkPosRTPCclusters2->Fill(KinkPosR, NTPCclusters);
		hQt->Fill(Qt);
		hKinkAngle->Fill(KinkAngle); 
		hDCAkink->Fill(KinkDistance);
		hPmKinkAng->Fill(P, KinkAngle);
		hKinkPosXY->Fill(KinkPos[0],KinkPos[1]);
		hKinkPosZY->Fill(KinkPos[2],KinkPos[1]);
		hKinkPosZR->Fill(KinkPos[2],KinkPosR);
		hKinkPosR->Fill(KinkPosR);
		hKinkPosZ->Fill(KinkPos[2]);
		hKinkPosZMCKinkPosZ->Fill(KinkPos[2], MCKinkPosZ);
 		hPdgCodemd->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
		hPmd->Fill(Pm,Pd);
		hMinvPimu->Fill(MinvPimu); 
		hUIDKinkDau->Fill(DauUID);
		hdEdx->Fill(PTPC, dEdx);

		if ( (IsRealPionKink(kink, MCstack) == kFALSE) || (Label1>MCNPrimTracks) ) {
			hPtKinkFake->Fill(Pt);
			hEtaKinkFake->Fill(Eta);
			hRapidityKinkFake->Fill(Rapidity);
			hPmPFake->Fill(Pm,P);
			hKinkPosRTPCclusters1Fake->Fill(NTPCclusters/KinkPosR);
			hKinkPosRTPCclusters2Fake->Fill(KinkPosR,NTPCclusters);
			hQtFake->Fill(Qt);
			hKinkAngleFake->Fill(KinkAngle); 
			hDCAkinkFake->Fill(KinkDistance);
			hPmKinkAngFake->Fill(P, KinkAngle);
			hKinkPosXYFake->Fill(KinkPos[0],KinkPos[1]);
			hKinkPosZYFake->Fill(KinkPos[2],KinkPos[1]);
			hKinkPosZRFake->Fill(KinkPos[2],KinkPosR);
			hKinkPosRFake->Fill(KinkPosR);
			hKinkPosZFake->Fill(KinkPos[2]);	
			hKinkPosZMCKinkPosZFake->Fill(KinkPos[2], MCKinkPosZ);
 			hPdgCodemdFake->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			hPmdFake->Fill(Pm,Pd);
			hMinvPimuFake->Fill(MinvPimu); 
			hUIDKinkDauFake->Fill(DauUID);
			hdEdxFake->Fill(PTPC, dEdx);
		}

		if ((KinkPosR<cLowR)||(KinkPosR>cUpR)) continue; //selection of kinks that are detected in the main region of TPC
 		
		hPtPosRSelected->Fill(Pt);

		if ((TMath::Abs(KinkPos[2])<0.5) || (TMath::Abs(KinkPos[2])>225)){ 
			hPdgCodemdZRejected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			continue; 
		}
		hPtZSelected->Fill(Pt);


		if  (KinkAngle<cLowKinkAngle) { //reject fake kinks (with very small kink angle)
			hPdgCodemdAngRejected->Fill(TMath::Abs(code1), TMath::Abs(code2));
			continue;
		}
		hPtAngSelected->Fill(Pt);

		if ((Pm/P<0.7) || (1.3<Pm/P)) { 
			hPdgCodemdPmRejected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			continue; 
		}//good tracks check mother track and kink if near
		hPtPmSelected->Fill(Pt);

		if (Qt<cLowQt) { //pion kinks Qt cuts
			hPdgCodemdQtLowRejected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			continue;
		}

		//Good Kinks
		hPtGoodKink->Fill(Pt);
		hEtaGoodKink->Fill(Eta);
		hRapidityGoodKink->Fill(Rapidity);
		hQtGoodKink->Fill(Qt);  
		hPmGoodKinkAng->Fill(P, KinkAngle); 
		hPdgCodemdGoodKink->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
		hPmdGoodKink->Fill(Pm,Pd); 
		hUIDGoodKinkDau->Fill(DauUID);
		hdEdxGoodKink->Fill(PTPC, dEdx);



		//---------------------- perfect PID from pdg code -------------------------//

		if (code1==cPdgPion) { //selection of pion plus kinks
			hUIDPiDauPlus->Fill(DauUID); 
			if ((Label1<=MCNPrimTracks) && (DauUID==4)) { //select primary pions that decay
				hMultPiPlus->Fill(NTracks); 
				hPtPiPlus->Fill(Pt);
				hEtaPiPlus->Fill(Eta); 
				hRapidityPiPlus->Fill(Rapidity); 
				hQtPiPlus->Fill(Qt); 
				hKinkAnglePiPlus->Fill(KinkAngle); 
				hPmKinkAngPiPlus->Fill(P, KinkAngle); 
				hKinkPosXYPiPlus->Fill(KinkPos[0],KinkPos[1]);
				hKinkPosZRPiPlus->Fill(KinkPos[2],KinkPosR);
				hKinkPosRPiPlus->Fill(KinkPosR);
				hDCAkinkPiPlus->Fill(KinkDistance);
		 		hPdgCodemdPiPlus->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
				hPmdPiPlus->Fill(Pm,Pd);
				hdEdxPiPlus->Fill(PTPC, dEdx);

				if (code2==-cPdgMuon) { //selection of muon daughters
					hQtPimuPlus->Fill(Qt); 
					hPmKinkAngPimuPlus->Fill(P, KinkAngle); 
				} else { //other daughters?
					hQtPiotherPlus->Fill(Qt); 
					hPmKinkAngPiotherPlus->Fill(P, KinkAngle); 
					hPdgCodemdPiotherPlus->Fill(TMath::Abs(code1), TMath::Abs(code2)); 
        	    		}
			}
		} else if (code1==-cPdgPion) {
			hUIDPiDauMinus->Fill(DauUID); 
			if ((Label1<=MCNPrimTracks) && (DauUID==4)) { //select primary pions that decay
				hMultPiMinus->Fill(NTracks); 
				hPtPiMinus->Fill(Pt);
				hEtaPiMinus->Fill(Eta); 
				hRapidityPiMinus->Fill(Rapidity); 
				hQtPiMinus->Fill(Qt); 
				hKinkAnglePiMinus->Fill(KinkAngle); 
				hPmKinkAngPiMinus->Fill(P, KinkAngle); 
				hKinkPosXYPiMinus->Fill(KinkPos[0],KinkPos[1]);
				hKinkPosZRPiMinus->Fill(KinkPos[2],KinkPosR);
				hKinkPosRPiMinus->Fill(KinkPosR);
				hDCAkinkPiMinus->Fill(KinkDistance);
		 		hPdgCodemdPiMinus->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
				hPmdPiMinus->Fill(Pm,Pd);
				hdEdxPiMinus->Fill(PTPC, dEdx);

				if (code2==cPdgMuon) { //selection of muon daughters
					hQtPimuMinus->Fill(Qt); 
					hPmKinkAngPimuMinus->Fill(P, KinkAngle); 
				} else { //other daughters?
					hQtPiotherMinus->Fill(Qt); 
					hPmKinkAngPiotherMinus->Fill(P, KinkAngle); 
					hPdgCodemdPiotherMinus->Fill(TMath::Abs(code1), TMath::Abs(code2)); 
        	    		}
			}
		}// end of pion selection


		//------------------------ realistic PID from physical criteria --------------------//



		if (Qt>cUpQt){ //pion kinks Qt cuts
			hPdgCodemdQtRejected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			continue;
		}
		hPtQtSelected->Fill(Pt);

		if  (KinkAngle>(MaxKinkAngPimu*1.1)) { //pion kink angle cuts
			hPdgCodemdMaxAngRejected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			continue;
		}
		hPtMaxAngSelected->Fill(Pt);

//		if ( ((NTPCclusters/KinkPosR)>0.63) || ((NTPCclusters/KinkPosR)<0.20) ) { //good TPC tracks selection
		Double_t tpcNClHigh = -51.67+ (11./12.)*KinkPosR;
		Double_t tpcNClMin  = -85.5 + (65./95.)*KinkPosR;
		if ( (NTPCclusters>tpcNClHigh) || (NTPCclusters<tpcNClMin) ){ //good TPC tracks selection  
			hPdgCodemdRTPCclustersRejected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			continue;
		} 
		hPtRTPCclustersSelected->Fill(Pt);
		hRTPCclustersRTPCclustersSelected->Fill(KinkPosR,NTPCclusters);

		if ((MinvPimu>cUpInvMass) || (MinvPimu<cLowInvMass)) { //pion kinks invariant mass cuts
			hPdgCodemdMinvRejected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			continue;
		}


/////////////////selected without dE/dx cut
		hPtSelected->Fill(Pt); 
		hEtaSelected->Fill(Eta); 
		hRapiditySelected->Fill(Rapidity);
		hQtSelected->Fill(Qt);
		hKinkAngleSelected->Fill(KinkAngle); 
		hDCAkinkSelected->Fill(KinkDistance);
		hPmKinkAngSelected->Fill(P, KinkAngle);
		hKinkPosXYSelected->Fill(KinkPos[0],KinkPos[1]);
		hKinkPosZRSelected->Fill(KinkPos[2],KinkPosR);
		hKinkPosRSelected->Fill(KinkPosR);
 		hPdgCodemdSelected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
		hPmdSelected->Fill(Pm,Pd);
		hMinvPimuSelected->Fill(MinvPimu); 
		hUIDKinkDauSelected->Fill(DauUID);
		hdEdxSelected->Fill(PTPC, dEdx);
		if ( (IsRealPionKink(kink, MCstack) == kFALSE) || (Label1>MCNPrimTracks) ) {
			hPtSelectedFake->Fill(Pt);
			hEtaSelectedFake->Fill(Eta);
			hRapiditySelectedFake->Fill(Rapidity);
			hQtSelectedFake->Fill(Qt);
			hKinkAngleSelectedFake->Fill(KinkAngle);
			hDCAkinkSelectedFake->Fill(KinkDistance);
			hPmKinkAngSelectedFake->Fill(P, KinkAngle);
			hKinkPosXYSelectedFake->Fill(KinkPos[0],KinkPos[1]);
			hKinkPosZRSelectedFake->Fill(KinkPos[2],KinkPosR);
			hKinkPosRSelectedFake->Fill(KinkPosR);
			hPmdSelectedFake->Fill(Pm,Pd); 
			hMinvPimuSelectedFake->Fill(MinvPimu); 
			hdEdxSelectedFake->Fill(PTPC, dEdx);
		}


		Double_t NSigmaTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ESDtrack, AliPID::kPion));
		if (NSigmaTPC>cSigmaCut)  { // dE/dx cut
			hPdgCodemddEdxRejected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			continue; 
		}

		Double_t Sign = ESDtrack->GetSign();

///////////////// dEdx selected
		hPtPiSelected->Fill(Pt); 
		hEtaPiSelected->Fill(Eta); 
		hRapidityPiSelected->Fill(Rapidity);
		hQtPiSelected->Fill(Qt);
		hKinkAnglePiSelected->Fill(KinkAngle); 
		hDCAkinkPiSelected->Fill(KinkDistance);
		hPmKinkAngPiSelected->Fill(P, KinkAngle);
		hKinkPosRTPCclusters1PiSelected->Fill(NTPCclusters/KinkPosR); 
		hKinkPosRTPCclusters2PiSelected->Fill(KinkPosR, NTPCclusters);
		hKinkPosXYPiSelected->Fill(KinkPos[0],KinkPos[1]);
		hKinkPosZRPiSelected->Fill(KinkPos[2],KinkPosR);
		hKinkPosRPiSelected->Fill(KinkPosR);
		hKinkPosZPiSelected->Fill(KinkPos[2]);
		hPmPPiSelected->Fill(Pm,P);
 		hPdgCodemdPiSelected->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
		hPmdPiSelected->Fill(Pm,Pd);
		hMinvPimuPiSelected->Fill(MinvPimu); 
		hUIDKinkDauPiSelected->Fill(DauUID);
		hdEdxPiSelected->Fill(PTPC, dEdx);
		
		if (Sign>0) {
			hPtPiSelectedPlus->Fill(Pt); 
			hEtaPiSelectedPlus->Fill(Eta); 
			hRapidityPiSelectedPlus->Fill(Rapidity);
			hQtPiSelectedPlus->Fill(Qt); 
			hKinkAnglePiSelectedPlus->Fill(KinkAngle);
			hDCAkinkPiSelectedPlus->Fill(KinkDistance);
			hPmKinkAngPiSelectedPlus->Fill(P, KinkAngle); 
			hKinkPosXYPiSelectedPlus->Fill(KinkPos[0],KinkPos[1]);
			hKinkPosZRPiSelectedPlus->Fill(KinkPos[2],KinkPosR);
			hPdgCodemdPiSelectedPlus->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			hPmdPiSelectedPlus->Fill(Pm,Pd); 
			hMinvPimuPiSelectedPlus->Fill(MinvPimu); 
			hUIDPiDaumuSelectedPlus->Fill(DauUID); 
			hdEdxPiSelectedPlus->Fill(PTPC, dEdx);
			
//------------------- contamination of results of realistic PID ----------------//	
			if (IsRealPionKink(kink, MCstack) == kTRUE) { //real pion kinks
				if (Label1<=MCNPrimTracks) { //primaries
					hPtPrimPiKinksPlus->Fill(Pt);
					hEtaPrimPiKinksPlus->Fill(Eta); 
					hRapidityPrimPiKinksPlus->Fill(Rapidity); 
				} else { //secondary real kinks
					hPtSecondPiKinksPlus->Fill(Pt); 
					hEtaSecondPiKinksPlus->Fill(Eta); 
					hRapiditySecondPiKinksPlus->Fill(Rapidity); 
				}
			} else { //fakes
				hPtNonPiKinksPlus->Fill(Pt);
				hEtaNonPiKinksPlus->Fill(Eta); 
				hRapidityNonPiKinksPlus->Fill(Rapidity);
				hPdgCodemdNonPiKinksPlus->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
        	}
						
		} else if (Sign<0) {
			hPtPiSelectedMinus->Fill(Pt);
			hEtaPiSelectedMinus->Fill(Eta); 
			hRapidityPiSelectedMinus->Fill(Rapidity);
			hQtPiSelectedMinus->Fill(Qt); 
			hKinkAnglePiSelectedMinus->Fill(KinkAngle);
			hDCAkinkPiSelectedMinus->Fill(KinkDistance);
			hPmKinkAngPiSelectedMinus->Fill(P, KinkAngle); 
			hKinkPosXYPiSelectedMinus->Fill(KinkPos[0],KinkPos[1]);
			hKinkPosZRPiSelectedMinus->Fill(KinkPos[2],KinkPosR);
			hPdgCodemdPiSelectedMinus->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
			hPmdPiSelectedMinus->Fill(Pm,Pd); 
			hMinvPimuPiSelectedMinus->Fill(MinvPimu);  
			hUIDPiDaumuSelectedMinus->Fill(DauUID); 
			hdEdxPiSelectedMinus->Fill(PTPC, dEdx);

//------------------- contamination of results of realistic PID ----------------//	
			if (IsRealPionKink(kink, MCstack) == kTRUE) { //real pion kinks
				if (Label1<=MCNPrimTracks) { //primaries
					hPtPrimPiKinksMinus->Fill(Pt);
					hEtaPrimPiKinksMinus->Fill(Eta); 
					hRapidityPrimPiKinksMinus->Fill(Rapidity); 
				} else { //secondary real kinks
					hPtSecondPiKinksMinus->Fill(Pt); 
					hEtaSecondPiKinksMinus->Fill(Eta); 
					hRapiditySecondPiKinksMinus->Fill(Rapidity); 
				}
			} else { //fakes
				hPtNonPiKinksMinus->Fill(Pt);
				hEtaNonPiKinksMinus->Fill(Eta); 
				hRapidityNonPiKinksMinus->Fill(Rapidity);
				hPdgCodemdNonPiKinksMinus->Fill(TMath::Abs(code1),TMath::Abs(code2)); 
        	}
		}
    } //end of ESD track loop
	PostData(1, fListOfHistos);
}   

    
//________________________________________________________________________
void AliAnalysisPionKinksMCESD::Terminate(Option_t *) {
}

//_________________________________________________________________________
const AliESDVertex* AliAnalysisPionKinksMCESD::GetEventVertex(AliESDEvent* esd) { //Gets ESD vertex and returns it if it is valid
  	const AliESDVertex* vertex = esd->GetPrimaryVertexTracks();
	if(vertex->GetStatus()==kTRUE) return vertex;
	else { 
		vertex = esd->GetPrimaryVertexSPD();
		if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>0)) return vertex;
		else return 0;
	}
}

//________________________________________________________________________
Double_t AliAnalysisPionKinksMCESD::Energy(AliESDtrack* track) const { //calculates the energy for a pion track
	Double_t TrackMom[3];
	track->GetPxPyPz(TrackMom);
	TVector3 P(TrackMom[0], TrackMom[1], TrackMom[2]);
	//Double_t EnergyK = TMath::Sqrt(P.Mag()*P.Mag()+0.493677*0.493677); //kaon's energy
	Double_t EnergyPi = TMath::Sqrt(P.Mag()*P.Mag()+0.13957018*0.13957018); //pion's energy
	return EnergyPi;
}

//________________________________________________________________________
Double_t AliAnalysisPionKinksMCESD::GetMCRapidity(TParticle* particle) const { //calculates the rapidity of a track
	TVector3 P(particle->Px(),particle->Py(),particle->Pz());
	Double_t MCRapidity = 0.5*(TMath::Log(((particle->Energy())+P[2])/((particle->Energy())-P[2])));
	return MCRapidity;
}

//________________________________________________________________________
Double_t AliAnalysisPionKinksMCESD::fuRapidity(AliESDtrack* track) const { //calculates the rapidity for a track
	Double_t TrackMom[3];
	track->GetPxPyPz(TrackMom);
	TVector3 P(TrackMom[0], TrackMom[1], TrackMom[2]);
	Double_t RapidityPi = 0.5*(TMath::Log((Energy(track)+P[2])/(Energy(track)-P[2])));
	return RapidityPi;
}

//________________________________________________________________________
Double_t AliAnalysisPionKinksMCESD::MCPQt(AliMCEvent* mcEvent, Int_t iMC, TParticle* MCdaughter) const {
	TVector3 DecayMomentumPi(0,0,0); //mother's momentum when it decays
	TClonesArray* trArray=0;
	TParticle* tempParticle=0;
	if (mcEvent->GetParticleAndTR(iMC, tempParticle, trArray) != -1) { 
    		AliTrackReference* MCtrackReference = static_cast<AliTrackReference*>(trArray->Last());
		DecayMomentumPi.SetXYZ(MCtrackReference->Px(), MCtrackReference->Py(), MCtrackReference->Pz());
	} else return 0;
  						
	const TVector3 MCP3d(MCdaughter->Px(), MCdaughter->Py(), MCdaughter->Pz()); //daughter's momentum when producedframe
	Double_t MCQt = MCP3d.Perp(DecayMomentumPi); //daughter's transverse momentum in mother's frame (Qt)
	return MCQt;
}

//________________________________________________________________________
Double_t AliAnalysisPionKinksMCESD::fuMCKinkAngle(AliMCEvent* mcEvent, Int_t iMC, TParticle* MCdaughter, Bool_t degrees) const {
	Double_t MCQt=MCPQt(mcEvent, iMC, MCdaughter);
	Double_t MCKinkAngle = TMath::ASin(MCQt/(MCdaughter->P())); //kink angle in rads
	if (degrees==kFALSE) {	
		return MCKinkAngle;
     } else {
		MCKinkAngle = TMath::RadToDeg()*MCKinkAngle; //kink angle in degrees
		return MCKinkAngle;
	}
}

//________________________________________________________________________
Bool_t AliAnalysisPionKinksMCESD::IsGoodTrack(AliESDtrack* ESDtrack) const { //Checks if a track is acceptable as good
	UInt_t status = ESDtrack->GetStatus();
	if ((status&AliESDtrack::kITSrefit)==0) return kFALSE; //cut tracks that cannot be reconstructed by ITS only
	if ((status&AliESDtrack::kTPCrefit)==0) return kFALSE; //cut tracks that cannot be reconstructed by TPC only
	
	Double_t NTPCclusters = ESDtrack->GetTPCclusters(0);
	Double_t TPCchi2 = ESDtrack->GetTPCchi2();
	if (NTPCclusters<20) return kFALSE;
	//if (NTPCclusters<30) return kFALSE; // cut sta 30 gia syst

	//hTPCchi2clusters->Fill(TPCchi2/NTPCclusters);//histo to check???????
	if (TPCchi2/NTPCclusters>3.8) return kFALSE; //cut tracks of bad quality fit with the contributing clusters

	/*Double_t ExtCov[15]; //external covariances matrix
	ESDtrack->GetExternalCovariance(ExtCov); 
	if(ExtCov[0]>2) return kFALSE; //sigma(y^2)
	if(ExtCov[2]>2) return kFALSE; //sigma(z^2)
	if(ExtCov[5]>0.5) return kFALSE; //sigma(sinphi^2) 
	if(ExtCov[9]>0.5) return kFALSE; //sigma(tanlamda^2)
	if(ExtCov[14]>2) return kFALSE; //sigma(1/Pt^2)*/
	
	return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisPionKinksMCESD::IsPrimaryTrack(AliESDtrack* ESDtrack) const { //Checks if a track is acceptable as a primary 
	Float_t ImpParam[2]; //DCA to vertex, 0->in x-y plane, 1->in z
	Float_t ImpParamCov[3]; //covariances of DCA to vertex
	ESDtrack->GetImpactParameters(ImpParam,ImpParamCov);
	//hdcaToVertexXY->Fill(ImpParam[0]);
	//hdcaToVertexZ->Fill(ImpParam[1]);
	if (ImpParamCov[0]<=0 || ImpParamCov[2]<=0) {
		AliDebug (1, "Estimated DCA covariance lower or equal zero!");
		ImpParamCov[0]=0; ImpParamCov[2]=0;
	}

	//if((TMath::Abs(ImpParam[0])>0.3) || (TMath::Abs((ImpParam[1])>2.5))) return kFALSE; //absolute DCA cut
	if (!fMaxDCAtoVtxCut->AcceptTrack(ESDtrack))	return kFALSE;
	else return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisPionKinksMCESD::IsRealPionKink(AliESDkink* kink, AliStack* MCstack) const { //checks if a reconstructed kink is a pion kink
	Int_t Label1 = kink->GetLabel(0); //mother's track MC label (first component is mother)
	Int_t Label2 = kink->GetLabel(1); //daughter's track MC label (second component is daughter)
	//if (Label1>MCNTracks) return kFALSE;
	//if (Label2>MCNTracks) return kFALSE;

	TParticle *particle1 = MCstack->Particle(TMath::Abs(Label1)); //mother MC particle object
	TParticle *particle2 = MCstack->Particle(TMath::Abs(Label2)); //daughter MC particle object
		
	Int_t code1 = particle1->GetPdgCode(); //mother's pdg code obtained from MC mother object
	Int_t code2 = particle2->GetPdgCode(); //daughter's pdg code obtained from MC daughter object

	UInt_t UID = particle2->GetUniqueID(); //ESD daughter's unique id (method of production)
	if (UID!=4) return kFALSE;

	if( (code1==cPdgPion) && (code2==-cPdgMuon) ){ 
		return kTRUE;
	} else if ( (code1==-cPdgPion) && (code2==cPdgMuon) ){
		return kTRUE;
	} else {
		return kFALSE;
	}
}

