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

#include "AliVParticle.h"
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
#include "AliESDtrackCuts.h"

#include "AliAnalysisPionKinksESD.h"

ClassImp(AliAnalysisPionKinksESD)

//________________________________________________________________________
AliAnalysisPionKinksESD::AliAnalysisPionKinksESD(const char *name)
:AliAnalysisTaskSE(name),
fMaxKinkAngKmu(0), 
fMaxKinkAngPimu(0), //functions

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
hPmd(0),
hMinvPimu(0),
hdEdx(0),
hPtPosRSelected(0),
hPtZSelected(0),
hPtAngSelected(0),
hPtPmSelected(0),
hPtGoodKink(0), 
hEtaGoodKink(0), 
hRapidityGoodKink(0), 
hQtGoodKink(0), 
hPmGoodKinkAng(0),   
hPmdGoodKink(0),
hdEdxGoodKink(0), 
hPtQtSelected(0),
hPtMaxAngSelected(0),
hPtRTPCclustersSelected(0),
hRTPCclustersRTPCclustersSelected(0),
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
hPmdSelected(0),
hMinvPimuSelected(0),  
hdEdxSelected(0), 
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
hPmdPiSelected(0),
hMinvPimuPiSelected(0),   
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
hPmdPiSelectedPlus(0),
hMinvPimuPiSelectedPlus(0),  
hdEdxPiSelectedPlus(0),
hPtPiSelectedMinus(0),
hEtaPiSelectedMinus(0), 
hRapidityPiSelectedMinus(0), 
hQtPiSelectedMinus(0), 
hKinkAnglePiSelectedMinus(0),
hDCAkinkPiSelectedMinus(0),
hPmKinkAngPiSelectedMinus(0), 
hKinkPosXYPiSelectedMinus(0), 
hKinkPosZRPiSelectedMinus(0),  
hPmdPiSelectedMinus(0),
hMinvPimuPiSelectedMinus(0), 
hdEdxPiSelectedMinus(0), // reconstruction histograms
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
void AliAnalysisPionKinksESD::UserCreateOutputObjects() {
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
hPmd = new TH2F("hPmd", "ESD mother vs daughter momentum magnitude; Mother's P (GeV/c); Daughter's P (GeV/c)", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimu = new TH1F("hMinvPimu", "Invariant mass of ESD pions decaying to muons; m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hdEdx = new TH2F("hdEdx", "dE/dx vs mother's momentum; p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hPtPosRSelected = new TH1F("hPtPosRSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y & R cuts); R (cm); dN/dR", nBinsPt, hLowPt, hUpPt);
hPtZSelected = new TH1F("hPtZSelected",  "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R & z cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPtAngSelected = new TH1F("hPtAngSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z & #theta cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPtPmSelected = new TH1F("hPtPmSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z, #theta & p_{track}/p{kink} cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPtGoodKink = new TH1F("hPtGoodKink", "Mother's transverse momentum for real ESD kinks (realistic PID); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaGoodKink = new TH1F("hEtaGoodKink", "Mother's pseudorapidity for real ESD kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityGoodKink = new TH1F("hRapidityGoodKink", "Mother's rapidity for real ESD kinks; n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtGoodKink = new TH1F("hQtGoodKink", "Daughter's transverse momentum in mother's frame for real ESD kinks (realistic PID); q_{T} (GeV/c); dN/dq_{T}", 200, 0.0, 0.3);
hPmGoodKinkAng = new TH2F("hPmGoodKinkAng", "Mother's momentum magnitude vs kink angle for real ESD kinks (realistic PID); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hPmdGoodKink = new TH2F("hPmdGoodKink", "Mother vs daughter momentum magnitude for real ESD kinks (realistic PID); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hdEdxGoodKink = new TH2F("hdEdxGoodKink", "dE/dx vs mother's momentum; p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hPtQtSelected = new TH1F("hPtQtSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink} & q^{T} cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPtMaxAngSelected = new TH1F("hPtMaxAngSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink}, q^{T} & #theta_{max} cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hPtRTPCclustersSelected = new TH1F("hPtRTPCclustersSelected", "Mother's transverse momentum for selected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink}, q^{T}, #theta_{max} & R/TPC clusters cuts); p_{T} (GeV/c); dN/dp_{T}", nBinsPt, hLowPt, hUpPt);
hRTPCclustersRTPCclustersSelected = new TH2F("hRTPCclustersRTPCclustersSelected", "Number of TPC clusters vs radius for selected kinks (DCA, quality, p_{T}, y, R, z, #theta, p_{track}/p{kink}, q^{T}, #theta_{max} & R/TPC clusters cuts); R (cm); number of TPC clusters", nBinsR, hLowR, hUpR, 100, 0, 200); 
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
hPmdSelected = new TH2F("hPmdSelected", "Mother vs daughter momentum magnitude for selected kinks (all cuts except dE/dx); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuSelected = new TH1F("hMinvPimuSelected", "Invariant mass for #pi^{#pm} decaying to #mu^{#pm} for selected kinks (all cuts except dE/dx); m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hdEdxSelected = new TH2F("hdEdxSelected", "dE/dx vs mother's momentum for selected kinks (all cuts except dE/dx); p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
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
hPmdPiSelected = new TH2F("hPmdPiSelected", "Mother vs daughter momentum magnitude for selected kinks (all cuts); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuPiSelected = new TH1F("hMinvPimuPiSelected", "Invariant mass for #pi^{#pm} decaying to #mu^{#pm} for selected kinks (all cuts); m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hdEdxPiSelected = new TH2F("hdEdxPiSelected", "dE/dx vs mother's momentum for selected kinks (all cuts); p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);
hPtPiSelectedPlus = new TH1F("hPtPiSelectedPlus", "Mother's transverse momentum for selected positive kinks (all cuts); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPiSelectedPlus = new TH1F("hEtaPiSelectedPlus", "Mother's pseudorapidity for selected positive kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityPiSelectedPlus = new TH1F("hRapidityPiSelectedPlus", "Mother's rapidity for selected positive kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtPiSelectedPlus = new TH1F("hQtKselectedPlus", "Daughter's transverse momentum in mother's frame for selected positive kinks (all cuts); q_{T} (GeV/c); dN/dq_{T}", 200, 0.0, 0.3);
hKinkAnglePiSelectedPlus = new TH1F("hKinkAnglePiSelectedPlus", "Kink angle for selected positive kinks (all cuts); #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hDCAkinkPiSelectedPlus = new TH1F("hDCAkinkPiSelectedPlus", "DCA for selected positive kinks (all cuts); DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPmKinkAngPiSelectedPlus = new TH2F("hPmKinkAngPiSelectedPlus", "Mother's momentum magnitude vs kink angle for selected positive kinks (all cuts); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXYPiSelectedPlus = new TH2F("hKinkPosXYPiSelectedPlus", "X-Y Position of selected positive kinks (all cuts); X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZRPiSelectedPlus = new TH2F("hKinkPosZRPiSelectedPlus", "Z vs radius of selected positive kinks (all cuts); Z (cm); R (cm)", nBinsZ, hLowZ, hUpZ,nBinsR, hLowR, hUpR);
hPmdPiSelectedPlus = new TH2F("hPmdPiSelectedPlus", "Mother vs daughter momentum magnitude for selected positive kinks (all cuts); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuPiSelectedPlus = new TH1F("hMinvPimuPiSelectedPlus", "Invariant mass for #pi^{+} decaying to #mu^{+} for selected kinks (all cuts); m (GeV/c^{2}); dN/dm",nBinsInvMass, hLowInvMass, hUpInvMass);
hdEdxPiSelectedPlus = new TH2F("hdEdxPiSelectedPlus", "dE/dx vs mother's momentum for selected positive kinks (all cuts); p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);

hPtPiSelectedMinus = new TH1F("hPtPiSelectedMinus", "Mother's transverse momentum for selected negative kinks (all cuts); p_{T} (GeV/c); dN/dp_{T}",   nBinsPt, hLowPt, hUpPt);
hEtaPiSelectedMinus = new TH1F("hEtaPiSelectedMinus", "Mother's pseudorapidity for selected negative kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hRapidityPiSelectedMinus = new TH1F("hRapidityPiSelectedMinus", "Mother's rapidity for selected negative kinks (all cuts); n; dN/dn", nBinsEta, hLowEta, hUpEta);
hQtPiSelectedMinus = new TH1F("hQtPiSelectedMinus", "Daughter's transverse momentum in mother's frame for selected negative kinks (all cuts); q_{T} (GeV/c); dN/dq_{T}", 200, 0.0, 0.3);
hKinkAnglePiSelectedMinus = new TH1F("hKinkAnglePiSelectedMinus", "Kink angle for selected negative kinks (all cuts); #theta (#circ); dN/d#theta", nBinsAngle, hLowAngle, hUpAngle);
hDCAkinkPiSelectedMinus = new TH1F("hDCAkinkPiSelectedMinus", "DCA for selected negative kinks (all cuts); DCA(cm); Number of kinks", 100, 0.0, 2.0);
hPmKinkAngPiSelectedMinus = new TH2F("hPmKinkAngPiSelectedMinus", "Mother's momentum magnitude vs kink angle for selected negative kinks (all cuts); p_{T} (GeV/c);  #theta (#circ); d^{2}N/dp_{T}d#theta", nBinsPt, hLowPt, hUpPt, nBinsAngle, hLowAngle, hUpAngle);
hKinkPosXYPiSelectedMinus = new TH2F("hKinkPosXYPiSelectedMinus", "X-Y Position for selected negative kinks (all cuts); X (cm); Y (cm)", nBinsXY, hLowXY, hUpXY, nBinsXY, hLowXY, hUpXY);
hKinkPosZRPiSelectedMinus = new TH2F("hKinkPosZRPiSelectedMinus", "Z vs radius for selected negative kinks (all cuts); Z (cm); R (cm)",nBinsZ, hLowZ, hUpZ,nBinsR, hLowR, hUpR);
hPmdPiSelectedMinus = new TH2F("hPmdPiSelectedMinus", "Mother vs daughter momentum magnitude for selected negative kinks (all cuts); Mother's P; Daughter's P", nBinsPt, hLowPt, hUpPt, nBinsPt, hLowPt, hUpPt);
hMinvPimuPiSelectedMinus = new TH1F("hMinvPimuPiSelectedMinus", "Invariant mass for #pi^{-} decaying to #mu^{-} for selected kinks (all cuts); m (GeV/c^{2}); dN/dm", nBinsInvMass, hLowInvMass, hUpInvMass);
hdEdxPiSelectedMinus = new TH2F("hdEdxPiSelectedMinus", "dE/dx vs mother's momentum for selected negative kinks (all cuts); p (GeV/c); dE/dx (GeV/cm)",   nBinsPt, hLowPt, hUpPt, 100, 0.0, 300.0);


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
fListOfHistos->Add(hPmd);
fListOfHistos->Add(hMinvPimu);
fListOfHistos->Add(hdEdx);
fListOfHistos->Add(hPtPosRSelected);
fListOfHistos->Add(hPtZSelected);
fListOfHistos->Add(hPtAngSelected);
fListOfHistos->Add(hPtPmSelected);
fListOfHistos->Add(hPtGoodKink);
fListOfHistos->Add(hEtaGoodKink);
fListOfHistos->Add(hRapidityGoodKink);
fListOfHistos->Add(hQtGoodKink);
fListOfHistos->Add(hPmGoodKinkAng);
fListOfHistos->Add(hPmdGoodKink);
fListOfHistos->Add(hdEdxGoodKink);
fListOfHistos->Add(hPtQtSelected);
fListOfHistos->Add(hPtMaxAngSelected);
fListOfHistos->Add(hPtRTPCclustersSelected);
fListOfHistos->Add(hRTPCclustersRTPCclustersSelected);
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
fListOfHistos->Add(hPmdSelected);
fListOfHistos->Add(hMinvPimuSelected);
fListOfHistos->Add(hdEdxSelected);
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
fListOfHistos->Add(hPmdPiSelected);
fListOfHistos->Add(hMinvPimuPiSelected);
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
fListOfHistos->Add(hPmdPiSelectedPlus);
fListOfHistos->Add(hMinvPimuPiSelectedPlus);
fListOfHistos->Add(hdEdxPiSelectedPlus);

fListOfHistos->Add(hPtPiSelectedMinus);
fListOfHistos->Add(hEtaPiSelectedMinus);
fListOfHistos->Add(hRapidityPiSelectedMinus);
fListOfHistos->Add(hQtPiSelectedMinus);
fListOfHistos->Add(hKinkAnglePiSelectedMinus);
fListOfHistos->Add(hDCAkinkPiSelectedMinus);
fListOfHistos->Add(hPmKinkAngPiSelectedMinus);
fListOfHistos->Add(hKinkPosXYPiSelectedMinus);
fListOfHistos->Add(hKinkPosZRPiSelectedMinus);
fListOfHistos->Add(hPmdPiSelectedMinus);
fListOfHistos->Add(hMinvPimuPiSelectedMinus);
fListOfHistos->Add(hdEdxPiSelectedMinus);

fListOfHistos->SetOwner(kTRUE);
PostData(1, fListOfHistos);
}

//________________________________________________________________________
void AliAnalysisPionKinksESD::UserExec(Option_t *) {
	AliVEvent *event = InputEvent();
	if (!event) {
		Printf("ERROR: Could not retrieve event");
		return;
	}

	AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(event);
	if (!esdEvent) {
		Printf("ERROR: Could not retrieve esd");
		return;
	}

	
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
		hPmd->Fill(Pm,Pd);
		hMinvPimu->Fill(MinvPimu); 
		hdEdx->Fill(PTPC, dEdx);

		if ((KinkPosR<cLowR)||(KinkPosR>cUpR)) continue; //selection of kinks that are detected in the main region of TPC
 		
		hPtPosRSelected->Fill(Pt);

		if ((TMath::Abs(KinkPos[2])<cLowZ) || (TMath::Abs(KinkPos[2])>cUpZ)) continue; 
		hPtZSelected->Fill(Pt);


		if  (KinkAngle<cLowKinkAngle) continue;
		hPtAngSelected->Fill(Pt);

		if ((Pm/P<0.7) || (1.3<Pm/P)) continue; 
		hPtPmSelected->Fill(Pt);

		if (Qt<cLowQt)continue;

		//Good Kinks
		hPtGoodKink->Fill(Pt);
		hEtaGoodKink->Fill(Eta);
		hRapidityGoodKink->Fill(Rapidity);
		hQtGoodKink->Fill(Qt);  
		hPmGoodKinkAng->Fill(P, KinkAngle); 
		hPmdGoodKink->Fill(Pm,Pd); 
		hdEdxGoodKink->Fill(PTPC, dEdx);

		//------------------------ realistic PID from physical criteria --------------------//

		if (Qt>cUpQt) continue;
		hPtQtSelected->Fill(Pt);

		if  (KinkAngle>(MaxKinkAngPimu*1.1)) continue;
		hPtMaxAngSelected->Fill(Pt);


//		if ( ((NTPCclusters/KinkPosR)>0.63) || ((NTPCclusters/KinkPosR)<0.20) ) //good TPC tracks selection
		Double_t tpcNClHigh = -51.67+ (11./12.)*KinkPosR;
		Double_t tpcNClMin  = -85.5 + (65./95.)*KinkPosR;
		if ( (NTPCclusters>tpcNClHigh) || (NTPCclusters<tpcNClMin) ) continue;//good TPC tracks selection
//		if ( (NTPCclusters>((11/12)*KinkPosR-51.67)) || (NTPCclusters<((65/95)*KinkPosR-85.5)) )  continue;
		hPtRTPCclustersSelected->Fill(Pt);
		hRTPCclustersRTPCclustersSelected->Fill(KinkPosR,NTPCclusters);

		if ((MinvPimu>cUpInvMass) || (MinvPimu<cLowInvMass)) continue;

		//selected
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
		hPmdSelected->Fill(Pm,Pd);
		hMinvPimuSelected->Fill(MinvPimu); 
		hdEdxSelected->Fill(PTPC, dEdx);


		Double_t NSigmaTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ESDtrack, AliPID::kPion));
		if (NSigmaTPC>cSigmaCut) continue; 
		Double_t Sign = ESDtrack->GetSign();

		// dEdx selected
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
		hPmdPiSelected->Fill(Pm,Pd);
		hMinvPimuPiSelected->Fill(MinvPimu); 
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
			hPmdPiSelectedPlus->Fill(Pm,Pd); 
			hMinvPimuPiSelectedPlus->Fill(MinvPimu); 
			hdEdxPiSelectedPlus->Fill(PTPC, dEdx);
			
						
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
			hPmdPiSelectedMinus->Fill(Pm,Pd); 
			hMinvPimuPiSelectedMinus->Fill(MinvPimu);  
			hdEdxPiSelectedMinus->Fill(PTPC, dEdx);

		}
    } //end of ESD track loop
	PostData(1, fListOfHistos);
}   

    
//________________________________________________________________________
void AliAnalysisPionKinksESD::Terminate(Option_t *) {
}

//_________________________________________________________________________
const AliESDVertex* AliAnalysisPionKinksESD::GetEventVertex(AliESDEvent* esd) { //Gets ESD vertex and returns it if it is valid
  	const AliESDVertex* vertex = esd->GetPrimaryVertexTracks();
	if(vertex->GetStatus()==kTRUE) return vertex;
	else { 
		vertex = esd->GetPrimaryVertexSPD();
		if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>0)) return vertex;
		else return 0;
	}
}

//________________________________________________________________________
Double_t AliAnalysisPionKinksESD::Energy(AliESDtrack* track) const { //calculates the energy for a pion track
	Double_t TrackMom[3];
	track->GetPxPyPz(TrackMom);
	TVector3 P(TrackMom[0], TrackMom[1], TrackMom[2]);
	//Double_t EnergyK = TMath::Sqrt(P.Mag()*P.Mag()+0.493677*0.493677); //kaon's energy
	Double_t EnergyPi = TMath::Sqrt(P.Mag()*P.Mag()+0.13957018*0.13957018); //pion's energy
	return EnergyPi;
}


//________________________________________________________________________
Double_t AliAnalysisPionKinksESD::fuRapidity(AliESDtrack* track) const { //calculates the rapidity for a track
	Double_t TrackMom[3];
	track->GetPxPyPz(TrackMom);
	TVector3 P(TrackMom[0], TrackMom[1], TrackMom[2]);
	Double_t RapidityK = 0.5*(TMath::Log((Energy(track)+P[2])/(Energy(track)-P[2])));
	return RapidityK;
}


//________________________________________________________________________
Bool_t AliAnalysisPionKinksESD::IsGoodTrack(AliESDtrack* ESDtrack) const { //Checks if a track is acceptable as good
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
Bool_t AliAnalysisPionKinksESD::IsPrimaryTrack(AliESDtrack* ESDtrack) const { //Checks if a track is acceptable as a primary 
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

