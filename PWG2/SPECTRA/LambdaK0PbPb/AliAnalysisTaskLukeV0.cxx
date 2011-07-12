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

/* $Id: AliAnalysisTaskLukeV0.cxx 46301 2011-01-06 14:25:27Z agheata $ */

/* AliAnalysisTaskLukeV0.cxx
 *
 * Task analysing lambda, antilambda & K0 spectra
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 * Adapted by Luke Hanratty
 *
 */

#include "AliAnalysisTaskLukeV0.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TPDGCode.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliPID.h"
#include "AliPIDResponse.h"

ClassImp(AliAnalysisTaskLukeV0)

//________________________________________________________________________
AliAnalysisTaskLukeV0::AliAnalysisTaskLukeV0() // All data members should be initialised here
   :AliAnalysisTaskSE(),
    fOutputList(0),
    fTrackCuts(0),
	fPIDResponse(0),
    fHistPt(0), 
    fHistEta(0),
	fHistLuke(0),
	fHistBetaV0(0), 
	fHistCosPA(0), 
	fHistCosTheta(0), 
	fHistCosThetaPi2(0), 
	fHistCosThetaProton2(0), 
	fHistDCAV0Daughters(0), 
	fHistDecayL(0), 
	fHistDecayLxy(0), 
	fHistDeltaTheta(0), 
	fHistdNV0sdT2(0), 
	fHistEtaNTracks(0),
	fHistEtaPTracks(0),
	fHistEtaV0s(0),
	fHistImpactxyN(0), 
	fHistImpactzN(0), 
	fHistImpactxyP(0), 
	fHistImpactzP(0), 
	fHistKinkIndexFalse(0), 
	fHistKinkIndexTrue(0),
	fHistLambdaBgRapidity(0),
	fHistLambdaBgEta(0),
	fHistLambdaBgPt(0),	
	fHistLambdaSigRapidity(0),
	fHistLambdaSigEta(0),
	fHistLambdaSigPt(0),
	fHistMagneticField(0),
	fHistMcLog(0),
	fHistMcNLambdaPrimary(0),
	fHistMcNLambda(0),
	fHistMcNAntilambda(0),
	fHistMcNKshort(0),
	fHistMcPtV0La(0),
	fHistMcPtV0Lb(0),
	fHistMcPtV0K0(0),
	fHistMcSigmaPProton(0),
	fHistMcSigmaPNProton(0),
	fHistMcSLambdaEta(0),
	fHistMcSLambdaDaughter(0),
	fHistMcSLambdaPt(0),
	fHistMcTPCTrackLength(0),
	fHistMK0(0), 
	fHistMK00(0), 
	fHistMK01(0), 
	fHistMK02(0), 
	fHistMK03(0), 
	fHistMK04(0), 
	fHistMK05(0), 
	fHistMK06(0), 
	fHistMK07(0), 
	fHistMK08(0), 
	fHistMK09(0),
	fHistMLa(0), 
	fHistMLa0(0), 
	fHistMLa1(0), 
	fHistMLa2(0), 
	fHistMLa3(0), 
	fHistMLa4(0), 
	fHistMLa5(0), 
	fHistMLa6(0), 
	fHistMLa7(0), 
	fHistMLa8(0), 
	fHistMLa9(0), 
	fHistMLb(0), 
	fHistMLb0(0), 
	fHistMLb1(0), 
	fHistMLb2(0), 
	fHistMLb3(0), 
	fHistMLb4(0), 
	fHistMLb5(0), 
	fHistMLb6(0), 
	fHistMLb7(0), 
	fHistMLb8(0), 
	fHistMLb9(0), 
	fHistNLambda(0), 
	fHistNV0(0), 
	fHistNTracks(0), 
	fHistPtV0(0), 
	fHistPVZ(0), 
	fHistTauLa(0), 
	fHistTheta(0), 
	fHistThetaPi2(0), 
	fHistThetaProton2(0), 
	fHistV0Z(0), 
	fHistZ(0),
	fHistBetaERatio(0), 
	fHistBetaPRatio(0), 
	fHistBetaPXRatio(0), 
	fHistBetheBlochITSNeg(0), 
	fHistBetheBlochITSPos(0), 
	fHistBetheBlochTPCNeg(0), 
	fHistBetheBlochTPCPos(0), 
	fHistCosPAMLa(0), 
	fHistDCAPtPSig(0),
	fHistDCAPtPBg(0),
	fHistDCAPtPbarSig(0),
	fHistDCAPtPbarBg(0),
	fHistDecayLDCA(0), 
	fHistDecayLMLa(0), 
	fHistDeltaThetaMLa(0), 
	fHistImpactxyImpactz(0), 
	fHistImpactxyMLa(0), 
	fHistImpactzMLa(0), 
	fHistMcLamdaPProductionVertex(0),
	fHistMcLamdaSProductionVertex(0),
	fHistMcLamdaSDecayVertex(0),
	fHistMcPMK0Pt(0),
	fHistMcPMLaPt(0),
	fHistMcPMLbPt(0),
	fHistMcSLambdaDaughterPairs(0),
	fHistMcV0MK0Pt(0),
	fHistMcV0MLaPt(0),
	fHistMcV0MLbPt(0),
	fHistMcV0LamdaSProductionVertex(0),
	fHistMcV0IDLamdaSProductionVertex(0),
	fHistMK0PtArm(0), 
	fHistMK0Pt(0), 
	fHistMK0R(0), 
	fHistMLaPtArm(0), 
	fHistMLaPt(0), 
	fHistMLaR(0), 
	fHistMLbPtArm(0), 
	fHistMLbPt(0), 
	fHistMLbR(0),
	fHistNegChi2PerClusterMLa(0), 
	fHistNegTPCRefitMLa(0), 
	fHistNTPCNegClustersMLa(0), 
	fHistNTPCPosClustersMLa(0), 
	fHistNV0sNTracks(0),  
	fHistPosChi2PerClusterMLa(0), 
	fHistPosTPCRefitMLa(0), 
	fHistPtArm(0), 
	fHistPtArmR(0), 
	fHistPtV0Z(0), 
	fHistRZ(0), 
	fHistTauLaMLa(0), 
	fHistXZ(0), 
	fHistYZ(0) // The last in the above list should not have a comma after it
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskLukeV0::AliAnalysisTaskLukeV0(const char *name) // All data members should be initialised here
   :AliAnalysisTaskSE(name),
    fOutputList(0),
	fTrackCuts(0),
	fPIDResponse(0),
	fHistPt(0), 
	fHistEta(0),
	fHistLuke(0),
	fHistBetaV0(0), 
	fHistCosPA(0), 
	fHistCosTheta(0), 
	fHistCosThetaPi2(0), 
	fHistCosThetaProton2(0), 
	fHistDCAV0Daughters(0), 
	fHistDecayL(0), 
	fHistDecayLxy(0), 
	fHistDeltaTheta(0), 
	fHistdNV0sdT2(0), 
	fHistEtaNTracks(0),
	fHistEtaPTracks(0),
	fHistEtaV0s(0),
	fHistImpactxyN(0), 
	fHistImpactzN(0), 
	fHistImpactxyP(0), 
	fHistImpactzP(0), 
	fHistKinkIndexFalse(0), 
	fHistKinkIndexTrue(0),
	fHistLambdaBgRapidity(0),
	fHistLambdaBgEta(0),
	fHistLambdaBgPt(0),	
	fHistLambdaSigRapidity(0),
	fHistLambdaSigEta(0),
	fHistLambdaSigPt(0),
	fHistMagneticField(0),
	fHistMcLog(0),
	fHistMcNLambdaPrimary(0),
	fHistMcNLambda(0),
	fHistMcNAntilambda(0),
	fHistMcNKshort(0),
	fHistMcPtV0La(0),
	fHistMcPtV0Lb(0),
	fHistMcPtV0K0(0),
	fHistMcSigmaPProton(0),
	fHistMcSigmaPNProton(0),
	fHistMcSLambdaEta(0),
	fHistMcSLambdaDaughter(0),
	fHistMcSLambdaPt(0),
	fHistMcTPCTrackLength(0),
	fHistMK0(0), 
	fHistMK00(0), 
	fHistMK01(0), 
	fHistMK02(0), 
	fHistMK03(0), 
	fHistMK04(0), 
	fHistMK05(0), 
	fHistMK06(0), 
	fHistMK07(0), 
	fHistMK08(0), 
	fHistMK09(0),
	fHistMLa(0), 
	fHistMLa0(0), 
	fHistMLa1(0), 
	fHistMLa2(0), 
	fHistMLa3(0), 
	fHistMLa4(0), 
	fHistMLa5(0), 
	fHistMLa6(0), 
	fHistMLa7(0), 
	fHistMLa8(0), 
	fHistMLa9(0), 
	fHistMLb(0), 
	fHistMLb0(0), 
	fHistMLb1(0), 
	fHistMLb2(0), 
	fHistMLb3(0), 
	fHistMLb4(0), 
	fHistMLb5(0), 
	fHistMLb6(0), 
	fHistMLb7(0), 
	fHistMLb8(0), 
	fHistMLb9(0), 
	fHistNLambda(0), 
	fHistNV0(0), 
	fHistNTracks(0), 
	fHistPtV0(0), 
	fHistPVZ(0), 
	fHistTauLa(0), 
	fHistTheta(0), 
	fHistThetaPi2(0), 
	fHistThetaProton2(0), 
	fHistV0Z(0), 
	fHistZ(0),
	fHistBetaERatio(0), 
	fHistBetaPRatio(0), 
	fHistBetaPXRatio(0), 
	fHistBetheBlochITSNeg(0), 
	fHistBetheBlochITSPos(0), 
	fHistBetheBlochTPCNeg(0), 
	fHistBetheBlochTPCPos(0), 
	fHistCosPAMLa(0), 
	fHistDCAPtPSig(0),
	fHistDCAPtPBg(0),
	fHistDCAPtPbarSig(0),
	fHistDCAPtPbarBg(0),
	fHistDecayLDCA(0), 
	fHistDecayLMLa(0), 
	fHistDeltaThetaMLa(0), 
	fHistImpactxyImpactz(0), 
	fHistImpactxyMLa(0), 
	fHistImpactzMLa(0), 
	fHistMcLamdaPProductionVertex(0),
	fHistMcLamdaSProductionVertex(0),
	fHistMcLamdaSDecayVertex(0),
	fHistMcPMK0Pt(0),
	fHistMcPMLaPt(0),
	fHistMcPMLbPt(0),
	fHistMcSLambdaDaughterPairs(0),
	fHistMcV0MK0Pt(0),
	fHistMcV0MLaPt(0),
	fHistMcV0MLbPt(0),
	fHistMcV0LamdaSProductionVertex(0),
	fHistMcV0IDLamdaSProductionVertex(0),
	fHistMK0PtArm(0), 
	fHistMK0Pt(0), 
	fHistMK0R(0), 
	fHistMLaPtArm(0), 
	fHistMLaPt(0), 
	fHistMLaR(0), 
	fHistMLbPtArm(0), 
	fHistMLbPt(0), 
	fHistMLbR(0),
	fHistNegChi2PerClusterMLa(0), 
	fHistNegTPCRefitMLa(0), 
	fHistNTPCNegClustersMLa(0), 
	fHistNTPCPosClustersMLa(0), 
	fHistNV0sNTracks(0), 
	fHistPosChi2PerClusterMLa(0), 
	fHistPosTPCRefitMLa(0), 
	fHistPtArm(0), 
	fHistPtArmR(0), 
	fHistPtV0Z(0), 
	fHistRZ(0), 
	fHistTauLaMLa(0), 
	fHistXZ(0), 
	fHistYZ(0)// The last in the above list should not have a comma after it
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskLukeV0::~AliAnalysisTaskLukeV0()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fOutputList;
    }
    if (fTrackCuts) delete fTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskLukeV0::UserCreateOutputObjects()
{
    // Create histograms
    // Called once (on the worker node)
        
    fOutputList = new TList();
    fOutputList->SetOwner();  // IMPORTANT!
    
	fTrackCuts = new AliESDtrackCuts();	
	
	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
	fPIDResponse = inputHandler->GetPIDResponse();

    // Create histograms - original test histograms
    Int_t ptBins = 15;
    Float_t ptLow = 0.1, ptUp = 3.1;
    fHistPt = new TH1F("fHistPt", "P_{T} distribution for reconstructed", ptBins, ptLow, ptUp);        
    Int_t etaBins = 40;
    Float_t etaLow = -2.0, etaUp = 2.0;
    fHistEta = new TH1F("fHistEta","#eta distribution for reconstructed",etaBins, etaLow, etaUp);
	fHistLuke = new TH1F("fHistLuke","Lukes Histogram",100, 0, 10);

	// lambda plot parameters
	int div = 96;
	float max = 1.2;
	float min = 1.08;
	
	// Create remaining histograms
	// TH1F first
	fHistBetaV0 = new	TH1F("fHistBetaV0","Beta of the v0 candidate; Beta; NV0s",50,0,1);
	fHistCosPA = new	TH1F("fHistCosPA", "Cosine of Pointing Angle of V0s; Cos PA; N(v0s)",202,0.8,1.01);
	fHistCosTheta = new	TH1F("fHistCosTheta","Cos Theta of decay in Lambda rest frame ;Cos Theta ;N V0s",180,-1,1);
	fHistCosThetaPi2 = new	TH1F("fHistCosThetaPi2","Cos Theta of pion decay in lab frame ;Cos Theta ;N V0s",180,-1,1);
	fHistCosThetaProton2 = new	TH1F("fHistCosThetaProton2","Cos Theta proton in lab frame ;Cos Theta ;N V0s",180,-1,1);
	fHistDCAV0Daughters = new	TH1F("fHistDCAV0Daughters", "DCA between V0 daughters; DCA (cm); N V0s", 100, 0, 2);
	fHistDecayL = new	TH1F("fHistDecayL", "Distance between V0 and PV; Distance(cm); N(v0s)",200,-0.1,30);
	fHistDecayLxy = new	TH1F("fHistDecayLxy", "Distance between V0 and PV in xy plane; Distance(cm); N(v0s)",200,-0.1,30);
	fHistDeltaTheta = new	TH1F("fHistDeltaTheta","Difference in theta of Proton and Pi in lambda rest frame; Delta Theta; NV0s",180,-7,7);
	fHistdNV0sdT2 = new	TH1F("fHistdNV0sdT2", "No of V0s over No of Tracks squared; NV0s/T^2; N", 30, 0, 2e-3);
	fHistEtaNTracks = new TH1F("fHistEtaNTracks","Eta for negative tracks; Eta; N",100,-5,5);
	fHistEtaPTracks = new TH1F("fHistEtaPTracks","Eta for positive tracks; Eta; N",100,-5,5);
	fHistEtaV0s = new TH1F("fHistEtaV0s","Eta for v0s; Eta; N",100,-5,5);
	fHistImpactxyN = new	TH1F("fHistImpactxyN", "RSM DCA between negative particle and primary vertex in xy plane; RSM DCA (cm); N(v0s)",100,0,1);
	fHistImpactzN = new	TH1F("fHistImpactzN", "RSM DCA between negative particle and primary vertex in z direction; RSM DCA (cm); N(v0s)",100,0,1);
	fHistImpactxyP = new	TH1F("fHistImpactxyP", "RSM DCA between positive particle and primary vertex in xy plane; RSM DCA (cm); N(v0s)",100,0,1);
	fHistImpactzP = new	TH1F("fHistImpactzP", "RSM DCA between positive particle and primary vertex in z direction; RSM DCA (cm); N(v0s)",100,0,1);
	fHistKinkIndexFalse = new	TH1F("fHistKinkIndexFalse","Lambda mass of non-kink candidates; M(p#pi^{-}) (GeV/c^{2})",96,1.08,1.2);
	fHistKinkIndexTrue = new	TH1F("fHistKinkIndexTrue","Lambda mass of kink candidates; M(p#pi^{-}) (GeV/c^{2})",96,1.08,1.2);
	fHistLambdaBgRapidity = new	TH1F("fHistLambdaBgRapidity","Rapidity of background V0s near Lambda mass; rapidity ",100,-1.5,1.5);
	fHistLambdaBgEta = new	TH1F("fHistLambdaBgEta","Psuedorapidity of background V0s near Lambda mass; Eta ",100,-1.5,1.5);
	fHistLambdaBgPt = new	TH1F("fHistLambdaBgPt","Transverse momentum of background V0s near Lambda mass; Pt(GeV/c) ",150,0,15);
	fHistLambdaSigRapidity = new	TH1F("fHistLambdaSigRapidity","Rapidity of signal V0s on Lambda mass; rapidity ",100,-1.5,1.5);
	fHistLambdaSigEta = new	TH1F("fHistLambdaSigEta","Psuedorapidity of signal V0s on Lambda mass; Eta ",100,-1.5,1.5);
	fHistLambdaSigPt = new	TH1F("fHistLambdaSigPt","Transverse momentum of signal V0s on Lambda mass; Pt(GeV/c) ",150,0,15);	
	fHistMagneticField = new	TH1F("fHistMagneticField", "Magnetic Field; Magnetic Field; N",100,-100,100);
	fHistMcLog = new	TH1F("fHistMcLog", "Log: 1=P ID, 2=P nID, 3=nP ID, 4=nP nID; Code; N",21,-0.25,10.25);
	fHistMcNLambdaPrimary = new	TH1F("fHistMcNLambdaPrimary","Number of primary lambdas in MC; NLambdas; i",6,-0.25,2.25);
	fHistMcNLambda = new	TH1F("fHistMcNLambda","Number of lambdas in MC; NLambdas; i",31,-0.5,30);
	fHistMcNAntilambda = new	TH1F("fHistMcNAntilambda","Number of antilambdas in MC; NAntiLambdas; i",31,-0.5,30);
	fHistMcNKshort = new	TH1F("fHistMcNKshort","Number of K0s in MC; NKshort; i",31,-0.5,30);
	fHistMcPtV0La = new	TH1F("fHistMcPtV0La","Pt distribution of V0s confirmed as lambdas; Pt (GeV/c); dN/dydPt",200,0,10);
	fHistMcPtV0Lb = new	TH1F("fHistMcPtV0Lb","Pt distribution of V0s confirmed as antilambdas; Pt (GeV/c); dN/dydPt",200,0,10);
	fHistMcPtV0K0 = new	TH1F("fHistMcPtV0K0","Pt distribution of V0s confirmed as K0s; Pt (GeV/c); dN/dydPt",200,0,10);
	fHistMcSigmaPProton = new	TH1F("fHistMcSigmaPProton","Sigma of being a proton from TPC response (MC protons); N Sigma; N",200,0,10);
	fHistMcSigmaPNProton = new	TH1F("fHistMcSigmaPNProton","Sigma of being a proton from TPC response (MC not protons); N Sigma; N",200,0,10);	
	fHistMcSLambdaEta = new TH1F("fHistMcSLambdaEta","#eta distribution for MC secondary lambdas",60, -3, 3);
	fHistMcSLambdaDaughter = new TH1F("fHistMcSLambda_Daughter","PDG code for first daughter of MC secondary lambdas",8001, -4000.5,4000.5);
	fHistMcSLambdaPt = new TH1F("fHistMcSLambdaPt","#Pt distribution for MC secondary lambdas",500, 0, 100);
	fHistMcTPCTrackLength = new TH1F("fHistMcTPCTrackLength","TPC track length for charged Lambda daughters",500, 0, 100);
	fHistMK0 = new	TH1F("fHistMK0","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK00 = new	TH1F("fHistMK00","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK01 = new	TH1F("fHistMK01","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK02 = new	TH1F("fHistMK02","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK03 = new	TH1F("fHistMK03","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK04 = new	TH1F("fHistMK04","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK05 = new	TH1F("fHistMK05","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK06 = new	TH1F("fHistMK06","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK07 = new	TH1F("fHistMK07","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK08 = new	TH1F("fHistMK08","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMK09 = new	TH1F("fHistMK09","K0Short Mass; M(#pi^{+}#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",140,0.414,0.582);
	fHistMLa = new	TH1F("fHistMLa","Lambda Mass; M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa0 = new	TH1F("fHistMLa0", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa1 = new	TH1F("fHistMLa1", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa2 = new	TH1F("fHistMLa2", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa3 = new	TH1F("fHistMLa3", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa4 = new	TH1F("fHistMLa4", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa5 = new	TH1F("fHistMLa5", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa6 = new	TH1F("fHistMLa6", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa7 = new	TH1F("fHistMLa7", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa8 = new	TH1F("fHistMLa8", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLa9 = new	TH1F("fHistMLa9", "V0 Mass; M(M(p#pi^{-}) (GeV/c^{2}); dN/dM (0.12 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb = new	TH1F("fHistMLb","AntiLambda Mass; M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb0 = new	TH1F("fHistMLb0", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb1 = new	TH1F("fHistMLb1", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb2 = new	TH1F("fHistMLb2", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb3 = new	TH1F("fHistMLb3", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb4 = new	TH1F("fHistMLb4", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb5 = new	TH1F("fHistMLb5", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb6 = new	TH1F("fHistMLb6", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb7 = new	TH1F("fHistMLb7", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb8 = new	TH1F("fHistMLb8", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistMLb9 = new	TH1F("fHistMLb9", "V0 Mass; M(M(#bar{p}#pi^{+}) (GeV/c^{2}); dN/dM (0.125 GeV/c^{2})^{-1}",div,min,max);
	fHistNLambda = new	TH1F("fHistNLambda", "Number of lambda per event; N(lambda); N(events)",50,-0.5,49.5);
	fHistNV0 = new	TH1F("fHistNV0","V0 frequency distribution; Number of V0 Candidates",1000,0,100000);
	fHistNTracks = new	TH1F("fHistNTracks", "Track frequency distribution; Number of Tracks; N Tracks", 1000, 0, 100000);
	fHistPtV0 = new	TH1F("fHistPtV0","V0 P_{T}; P_{T} (GeV/c);dN/dP_{T} (GeV/c)^{-1}",40,0.,4.);
	fHistPVZ = new	TH1F("fHistPVZ","Z primary; Z (cm); Counts",100,-10,10);
	fHistTauLa = new	TH1F("fHistTauLa", "Lifetime under Lambda mass hypothesis; Lifetime(s); N(v0s)",200,0,100);
	fHistTheta = new	TH1F("fHistTheta","Angle of decay in Lambda rest frame ;Theta ;N V0s",180,-4,4);
	fHistThetaPi2 = new	TH1F("fHistThetaPi2","Angle of pion decay in lab frame ;Theta ;N V0s",180,-4,4);
	fHistThetaProton2 = new	TH1F("fHistThetaProton2","Angle of proton decay in lab frame;Theta ;N V0s",180,-4,4);
	fHistV0Z = new	TH1F("fHistV0Z","Z decay; Z (cm); Counts",100,-10,10);
	fHistZ = new	TH1F("fHistZ","Z decay - Z primary; Z (cm); Counts",100,-10,10);
	
	//TH2F follow
	fHistBetaERatio = new	TH2F("fHistBetaERatio","Ratio of Energies of Pion & Proton ;Beta ;Ratio",100,0,1,100,0,1);
	fHistBetaPRatio = new	TH2F("fHistBetaPRatio","Ratio of Momentum of Pion & Proton ;Beta ;Ratio",100,0,1,100,0,2);
	fHistBetaPXRatio = new	TH2F("fHistBetaPXRatio","Ratio of Momentum in boost direction of Pion & Proton ;Beta ;Ratio",100,0,1,100,-10,10);
	fHistBetheBlochITSNeg = new	TH2F("fHistBetheBlochITSNeg","-dE/dX against Momentum for negative daughter from ITS; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
	fHistBetheBlochITSPos = new	TH2F("fHistBetheBlochITSPos","-dE/dX against Momentum for positive daughter from ITS; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
	fHistBetheBlochTPCNeg = new	TH2F("fHistBetheBlochTPCNeg","-dE/dX against Momentum for negative daughter from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
	fHistBetheBlochTPCPos = new	TH2F("fHistBetheBlochTPCPos","-dE/dX against Momentum for positive daughter from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
	fHistCosPAMLa = new	TH2F("fHistCosPAMLa", "Cosine of Pointing Angle of V0s; Cos PA; M(p#pi^{-}) (GeV/c^{2})",202,0.8,1.01,96,1.08,1.2);
	fHistDCAPtPSig = new	TH2F("fHistDCAPtPSig","DCA of proton daughter from Lambda in peak region to PV versus Pt; P_{perp} (GeV/c); DCA(cm)",200,0,10,200,0,20);
	fHistDCAPtPBg = new	TH2F("fHistDCAPtPBg","DCA of proton daughter from Lambda in sideband region to PV versus Pt; P_{perp} (GeV/c); DCA(cm)",200,0,10,200,0,20);
	fHistDCAPtPbarSig = new	TH2F("fHistDCAPtPbarSig","DCA of antiproton daughter from antiLambda in peak region to PV versus Pt; P_{perp} (GeV/c); DCA(cm)",200,0,10,200,0,20);
	fHistDCAPtPbarBg = new	TH2F("fHistDCAPtPbarBg","DCA of antiproton daughter from antiLambda in sideband region to PV versus Pt; P_{perp} (GeV/c); DCA(cm)",200,0,10,200,0,20);
	fHistDecayLDCA = new	TH2F("fHistDecayLDCA", "Distance between V0 and PV against DCA between daughter tracks; Distance(cm); DCA (cm) ",301,-0.1,30,100,0,2);
	fHistDecayLMLa = new	TH2F("fHistDecayLMLa", "Distance between V0 and PV; Distance(cm); M(p#pi^{-}) (GeV/c^{2})",301,-0.1,30,96,1.08,1.2);
	fHistDeltaThetaMLa = new	TH2F("fHistDeltaThetaMLa"," Delta theta against effective lambda mass; Delta Theta; MLa", 180, -3, 3, 96, 1.08, 1.2);
	fHistImpactxyImpactz = new	TH2F("fHistImpactxyImpactz", "RSM DCA between negative particle and primary vertex in xy plane; RSM DCA xy (cm); RSM DCA z (cm)",100,0,1,100,0,10);
	fHistImpactxyMLa = new	TH2F("fHistImpactxyMLa", "RSM DCA between negative particle and primary vertex in xy plane; RSM DCA (cm); M(p#pi^{-}) (GeV/c^{2})",100,0,1,96,1.08,1.2);
	fHistImpactzMLa = new	TH2F("fHistImpactzMLa", "RSM DCA between negative particle and primary vertex in xy plane; RSM DCA (cm); M(p#pi^{-}) (GeV/c^{2})",100,0,1,96,1.08,1.2);
	fHistMcLamdaPProductionVertex = new	TH2F("fHistMcLamdaPProductionVertex", "Production vertex of primary Lambdas in RZ plane; Zv; Rv",100,-100,100,100,0,100);
	fHistMcLamdaSProductionVertex = new	TH2F("fHistMcLamdaSProductionVertex", "Production vertex of secondary Lambdas in RZ plane; Zv; Rv",100,-100,100,100,0,100);
	fHistMcLamdaSDecayVertex = new	TH2F("fHistMcLamdaSDecayVertex", "Decay vertex of secondary Lambdas in RZ plane; Zv; Rv",100,-100,100,100,0,100);
	fHistMcPMK0Pt = new	TH2F("fHistMcPMK0Pt","Monte Carlo primary K0 Mass versus Pt; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPt = new	TH2F("fHistMcPMLaPt","Monte Carlo primary (& sigma0) Lambda Mass versus Pt; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPt = new	TH2F("fHistMcPMLbPt","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcSLambdaDaughterPairs = new	TH2F("fHistMcSLambdaDaughterPairs", "PDG codes of daughter particles of secondary Lambdas; Daughter1; Daughter2",21,-10.5,10.5,21,-10.5,10.5);
	fHistMcV0MK0Pt = new	TH2F("fHistMcV0MK0Pt","Monte Carlo V0s passing cuts. K0 Mass versus Pt; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcV0MLaPt = new	TH2F("fHistMcV0MLaPt","Monte Carlo V0s passing cuts. Lambda Mass versus Pt; P_{perp} (GeV/c); Lambda Mass (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcV0MLbPt = new	TH2F("fHistMcV0MLbPt","Monte Carlo V0s passing cuts. Antilambda Mass versus Pt; P_{perp} (GeV/c); Antilambda Mass (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcV0LamdaSProductionVertex = new	TH2F("fHistMcV0LamdaSProductionVertex", "Production vertex of V0 secondary Lambdas in RZ plane; Zv; Rv",100,-100,100,100,0,100);
	fHistMcV0IDLamdaSProductionVertex = new	TH2F("fHistMcV0IDLamdaSProductionVertex", "Production vertex of identified V0 secondary Lambdas in RZ plane; Zv; Rv",100,-100,100,100,0,100);	
	fHistMK0PtArm = new	TH2F("fHistMK0PtArm","K0 Mass versus PtArm; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",40,0,0.25,140,0.414,0.582);
	fHistMK0Pt = new	TH2F("fHistMK0Pt","K0 Mass versus Pt; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMK0R = new	TH2F("fHistMK0R","K0 Mass versus R; R (cm); K0 Mass (GeV/c^2)",120,0,120,140,0.414,0.582);
	fHistMLaPtArm = new	TH2F("fHistMLaPtArm","Lambda Mass versus PtArm; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",40,0,0.25,96,1.08,1.2);
	fHistMLaPt = new	TH2F("fHistMLaPt","Lambda Mass versus Pt; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLaR = new	TH2F("fHistMLaR","Lambda Mass versus R; R (cm); M(p#pi^{-}) (GeV/c^2)",120,0,120,96,1.08,1.2);
	fHistMLbPtArm = new	TH2F("fHistMLbPtArm","AntiLambda Mass versus PtArm; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",40,0,0.25,96,1.08,1.2);
	fHistMLbPt = new	TH2F("fHistMLbPt","AntiLambda Mass versus Pt; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbR = new	TH2F("fHistMLbR","AntiLambda Mass versus R; R (cm); M(#bar{p}#pi^{+}) (GeV/c^2)",120,0,120,96,1.08,1.2);
	fHistNegChi2PerClusterMLa = new	TH2F("fHistNegChi2PerClusterMLa","Chi Squared per cluster against Lambda mass (negative Track);Chi^2 per cluster; M(p#pi^{-}) (GeV/c^{2})",200,0,7,96,1.08,1.2);; 
	fHistNegTPCRefitMLa = new	TH2F("fHistNegTPCRefitMLa","TPC refit?; TPC refit?; M(p#pi^{-}) (GeV/c^{2})",200,0,100,96,1.08,1.2);; 
	fHistNTPCNegClustersMLa = new	TH2F("fHistNTPCNegClustersMLa","Number of TPC clusters per negative track against Lambda mass;Number of TPC clusters; M(p#pi^{-}) (GeV/c^{2})",200,0,200,96,1.08,1.2);; 
	fHistNTPCPosClustersMLa = new	TH2F("fHistNTPCPosClustersMLa","Number of TPC clusters per positive track against Lambda mass;Number of TPC clusters; M(p#pi^{-}) (GeV/c^{2})",200,0,200,96,1.08,1.2);; 
	fHistNV0sNTracks = new	TH2F("fHistNV0sNTracks", "Number of tracks squared against number of v0s; NTracks^2; NV0s", 1000, 0, 1e8, 1000, 0, 1.5e5);
	fHistPosChi2PerClusterMLa = new	TH2F("fHistPosChi2PerClusterMLa","Chi Squared per cluster against Lambda mass (positive Track);Chi^2 per cluster; M(p#pi^{-}) (GeV/c^{2})",200,0,7,96,1.08,1.2);; 
	fHistPosTPCRefitMLa = new	TH2F("fHistPosTPCRefitMLa","TPC refit?; TPC refit?; M(p#pi^{-}) (GeV/c^{2})",200,0,100,96,1.08,1.2);; 
	fHistPtArm = new	TH2F("fHistPtArm","Podolanski-Armenteros Plot; #alpha; P_{perp} (GeV/c)",40,-1,1,80,0,0.5);
	fHistPtArmR = new	TH2F("fHistPtArmR","PtArm versus R; R decay (cm); P_{perp}",120,0,120,50,0,0.25);
	fHistPtV0Z = new	TH2F("fHistPtV0Z","Pt of V0 vs Z position; Pt (GeV/c); Z (cm)",200,-0.1,1.9,200,-50,50);
	fHistRZ = new	TH2F("fHistRZ","R decay versus Z decay; Z (cm); R (cm)",100,-50,50,120,0,220);
	fHistTauLaMLa = new	TH2F("fHistTauLaMLa", "Lifetime under Lambda mass hypothesis; Lifetime(s); M(p#pi^{-}) (GeV/c^{2})",200,0,100,96,1.08,1.2);
	fHistXZ = new	TH2F("fHistXZ","X decay versus Z decay; Z (cm); X (cm)",100,-50,50,200,-200,200);
	fHistYZ = new	TH2F("fHistYZ","Y decay versus Z decay; Z (cm); Y (cm)",100,-50,50,200,-200,200);  
	
	
        
	// All histograms must be added to output list
        
	fOutputList->Add(fHistPt); 
    fOutputList->Add(fHistEta);
	fOutputList->Add(fHistLuke);
	fOutputList->Add(fHistBetaV0); 
	fOutputList->Add(fHistCosPA); 
	fOutputList->Add(fHistCosTheta); 
	fOutputList->Add(fHistCosThetaPi2); 
	fOutputList->Add(fHistCosThetaProton2); 
	fOutputList->Add(fHistDCAV0Daughters); 
	fOutputList->Add(fHistDecayL); 
	fOutputList->Add(fHistDecayLxy); 
	fOutputList->Add(fHistDeltaTheta); 
	fOutputList->Add(fHistdNV0sdT2); 
	fOutputList->Add(fHistEtaNTracks);
	fOutputList->Add(fHistEtaPTracks);
	fOutputList->Add(fHistEtaV0s);
	fOutputList->Add(fHistImpactxyN); 
	fOutputList->Add(fHistImpactzN); 
	fOutputList->Add(fHistImpactxyP); 
	fOutputList->Add(fHistImpactzP); 
	fOutputList->Add(fHistKinkIndexFalse); 
	fOutputList->Add(fHistKinkIndexTrue);
	fOutputList->Add(fHistLambdaBgRapidity);
	fOutputList->Add(fHistLambdaBgEta);
	fOutputList->Add(fHistLambdaBgPt);	
	fOutputList->Add(fHistLambdaSigRapidity);
	fOutputList->Add(fHistLambdaSigEta);
	fOutputList->Add(fHistLambdaSigPt);
	fOutputList->Add(fHistMagneticField);
	fOutputList->Add(fHistMcLog);
	fOutputList->Add(fHistMcNLambdaPrimary);
	fOutputList->Add(fHistMcNLambda);
	fOutputList->Add(fHistMcNAntilambda);
	fOutputList->Add(fHistMcNKshort);
	fOutputList->Add(fHistMcPtV0La);
	fOutputList->Add(fHistMcPtV0Lb);
	fOutputList->Add(fHistMcPtV0K0);
	fOutputList->Add(fHistMcSigmaPProton);
	fOutputList->Add(fHistMcSigmaPNProton);
	fOutputList->Add(fHistMcSLambdaEta);
	fOutputList->Add(fHistMcSLambdaDaughter);
	fOutputList->Add(fHistMcSLambdaPt);
	fOutputList->Add(fHistMcTPCTrackLength);
	fOutputList->Add(fHistMK0); 
	fOutputList->Add(fHistMK00); 
	fOutputList->Add(fHistMK01); 
	fOutputList->Add(fHistMK02); 
	fOutputList->Add(fHistMK03); 
	fOutputList->Add(fHistMK04); 
	fOutputList->Add(fHistMK05); 
	fOutputList->Add(fHistMK06); 
	fOutputList->Add(fHistMK07); 
	fOutputList->Add(fHistMK08); 
	fOutputList->Add(fHistMK09);
	fOutputList->Add(fHistMLa); 
	fOutputList->Add(fHistMLa0); 
	fOutputList->Add(fHistMLa1); 
	fOutputList->Add(fHistMLa2); 
	fOutputList->Add(fHistMLa3); 
	fOutputList->Add(fHistMLa4); 
	fOutputList->Add(fHistMLa5); 
	fOutputList->Add(fHistMLa6); 
	fOutputList->Add(fHistMLa7); 
	fOutputList->Add(fHistMLa8); 
	fOutputList->Add(fHistMLa9); 
	fOutputList->Add(fHistMLb); 
	fOutputList->Add(fHistMLb0); 
	fOutputList->Add(fHistMLb1); 
	fOutputList->Add(fHistMLb2); 
	fOutputList->Add(fHistMLb3); 
	fOutputList->Add(fHistMLb4); 
	fOutputList->Add(fHistMLb5); 
	fOutputList->Add(fHistMLb6); 
	fOutputList->Add(fHistMLb7); 
	fOutputList->Add(fHistMLb8); 
	fOutputList->Add(fHistMLb9); 
	fOutputList->Add(fHistNLambda); 
	fOutputList->Add(fHistNV0); 
	fOutputList->Add(fHistNTracks); 
	fOutputList->Add(fHistPtV0); 
	fOutputList->Add(fHistPVZ); 
	fOutputList->Add(fHistTauLa); 
	fOutputList->Add(fHistTheta); 
	fOutputList->Add(fHistThetaPi2); 
	fOutputList->Add(fHistThetaProton2); 
	fOutputList->Add(fHistV0Z); 
	fOutputList->Add(fHistZ);
	fOutputList->Add(fHistBetaERatio); 
	fOutputList->Add(fHistBetaPRatio); 
	fOutputList->Add(fHistBetaPXRatio); 
	fOutputList->Add(fHistBetheBlochITSNeg); 
	fOutputList->Add(fHistBetheBlochITSPos); 
	fOutputList->Add(fHistBetheBlochTPCNeg); 
	fOutputList->Add(fHistBetheBlochTPCPos); 
	fOutputList->Add(fHistCosPAMLa); 
	fOutputList->Add(fHistDCAPtPSig);
	fOutputList->Add(fHistDCAPtPBg);
	fOutputList->Add(fHistDCAPtPbarSig);
	fOutputList->Add(fHistDCAPtPbarBg);
	fOutputList->Add(fHistDecayLDCA); 
	fOutputList->Add(fHistDecayLMLa); 
	fOutputList->Add(fHistDeltaThetaMLa); 
	fOutputList->Add(fHistImpactxyImpactz); 
	fOutputList->Add(fHistImpactxyMLa); 
	fOutputList->Add(fHistImpactzMLa); 
	fOutputList->Add(fHistMcLamdaPProductionVertex);
	fOutputList->Add(fHistMcLamdaSProductionVertex);
	fOutputList->Add(fHistMcLamdaSDecayVertex);
	fOutputList->Add(fHistMcPMK0Pt);
	fOutputList->Add(fHistMcPMLaPt);
	fOutputList->Add(fHistMcPMLbPt);
	fOutputList->Add(fHistMcSLambdaDaughterPairs);
	fOutputList->Add(fHistMcV0MK0Pt);
	fOutputList->Add(fHistMcV0MLaPt);
	fOutputList->Add(fHistMcV0MLbPt);
	fOutputList->Add(fHistMcV0LamdaSProductionVertex);
	fOutputList->Add(fHistMcV0IDLamdaSProductionVertex);
	fOutputList->Add(fHistMK0PtArm); 
	fOutputList->Add(fHistMK0Pt); 
	fOutputList->Add(fHistMK0R); 
	fOutputList->Add(fHistMLaPtArm); 
	fOutputList->Add(fHistMLaPt); 
	fOutputList->Add(fHistMLaR); 
	fOutputList->Add(fHistMLbPtArm); 
	fOutputList->Add(fHistMLbPt); 
	fOutputList->Add(fHistMLbR);
	fOutputList->Add(fHistNegChi2PerClusterMLa); 
	fOutputList->Add(fHistNegTPCRefitMLa); 
	fOutputList->Add(fHistNTPCNegClustersMLa); 
	fOutputList->Add(fHistNTPCPosClustersMLa); 
	fOutputList->Add(fHistNV0sNTracks); 
	fOutputList->Add(fHistPosChi2PerClusterMLa); 
	fOutputList->Add(fHistPosTPCRefitMLa); 
	fOutputList->Add(fHistPtArm); 
	fOutputList->Add(fHistPtArmR); 
	fOutputList->Add(fHistPtV0Z); 
	fOutputList->Add(fHistRZ); 
	fOutputList->Add(fHistTauLaMLa); 
	fOutputList->Add(fHistXZ); 
	fOutputList->Add(fHistYZ);
	
	
    PostData(1, fOutputList); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskLukeV0::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
	
	// paramaters used for most cuts, to minimise editing
	double cutCosPa(0.998), cutcTau(2);
	double cutNImpact(-999), cutDCA(0.4);
	double cutBetheBloch(3);
	double cutMinNClustersTPC(70), cutMaxChi2PerClusterTPC(-999);
	double isMonteCarlo(false);
	double cutEta(0.8);
	
	//Track Cuts set here
	if(cutMinNClustersTPC != -999)
	{(fTrackCuts->SetMinNClustersTPC(int(cutMinNClustersTPC)));}
	if(cutMaxChi2PerClusterTPC != -999)
	{fTrackCuts->SetMaxChi2PerClusterTPC(cutMaxChi2PerClusterTPC);}
	fTrackCuts->SetAcceptKinkDaughters(kFALSE); 
	fTrackCuts->SetRequireTPCRefit(kTRUE);
	
        
    // Create pointer to reconstructed event

	AliVEvent *event = InputEvent();
	if (!event) { Printf("ERROR: Could not retrieve event"); return; }
	
	// create pointer to event
	AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(event);
	if (!fESD) {
		AliError("Cannot get the ESD event");
		return;
	}

	/*********************************************************************/
	// MONTE CARLO SECTION
	// This section loops over all MC tracks

	int nLambdaMC = 0;
	int nAntilambdaMC = 0;
	int nKshortMC = 0;
	
	if(isMonteCarlo) 
	{

		// If the task accesses MC info, this can be done as in the commented block below:
		
		// Create pointer to reconstructed event
		AliMCEvent *mcEvent = MCEvent();
		if (!mcEvent) 
			{ 
				Printf("ERROR: Could not retrieve MC event"); 
				//return; 
			}
		else
			{
				Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
			}
			 
		// set up a stack for use in check for primary/stable particles
		AliStack* mcStack = mcEvent->Stack();
		if( !mcStack ) { Printf( "Stack not available"); return; }
		
		AliMCVertex *mcpv = (AliMCVertex *) mcEvent->GetPrimaryVertex();
		Double_t mcpvPos[3];
		if (mcpv != 0)
		{
			mcpv->GetXYZ(mcpvPos);
		}
		else 
		{
			Printf("ERROR: Could not resolve MC primary vertex");
			return;
		}
		
		//loop over all MC tracks
		for(Int_t iMCtrack = 0; iMCtrack < mcEvent->GetNumberOfTracks(); iMCtrack++)
		{
			
			//booleans to check if track is La, Lb, K0 and primary
			bool lambdaMC = false;
			bool antilambdaMC = false;
			bool kshortMC = false;
			bool isprimaryMC = false;
			
			AliMCParticle *mcPart = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(iMCtrack));
			
			if(mcPart->PdgCode() == kLambda0)
			{
				lambdaMC = true;
				nLambdaMC++;
			}
			else if(mcPart->PdgCode() == kK0Short)
			{
				kshortMC = true;
				nKshortMC++;
			}
			else if(mcPart->PdgCode() == kLambda0Bar)
			{
				antilambdaMC = true;
				nAntilambdaMC++;
			}
			
			//if only interested in La,Lb,K0, can use this to terminate loop for other particles
			//if(lambdaMC == false && antilambdaMC == false && kshortMC == false)
			//{continue;}
			
			Int_t firstDaughterLabel = mcPart->GetFirstDaughter();
			Int_t lastDaughterLabel = mcPart->GetLastDaughter();
			Int_t motherLabel = mcPart->GetMother();
			
			if(firstDaughterLabel < 0 || lastDaughterLabel < 0)
			{continue;}
			
			AliMCParticle *mcFirstDaughter = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(firstDaughterLabel));
			AliMCParticle *mcLastDaughter = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(lastDaughterLabel));
			AliMCParticle *mcMother = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(motherLabel));
		
			double mcRadius = TMath::Sqrt((mcPart->Xv())*(mcPart->Xv())+(mcPart->Yv())*(mcPart->Yv()));
			double mcFirstDaughterRadius = TMath::Sqrt((mcFirstDaughter->Xv())*(mcFirstDaughter->Xv())+(mcFirstDaughter->Yv())*(mcFirstDaughter->Yv()));
			double mcLastDaughterRadius = TMath::Sqrt((mcLastDaughter->Xv())*(mcLastDaughter->Xv())+(mcLastDaughter->Yv())*(mcLastDaughter->Yv()));			
			
			double motherType = -1;
			if(motherLabel >= 0)
			{motherType = mcMother->PdgCode();}
			
			// this block of code is used to include Sigma0 decays as primary lambda/antilambda
			bool sigma0MC = false;
			if(motherType == 3212 || motherType == -3212)
				{
				if(mcEvent->IsPhysicalPrimary(motherLabel))
				   {sigma0MC = true;}
				}
				   
			
			if(mcEvent->IsPhysicalPrimary(iMCtrack) || sigma0MC)
			   {
				   isprimaryMC = true;
				   if(lambdaMC)
				   {
					   fHistMcNLambdaPrimary->Fill(1);
					   fHistMcLamdaPProductionVertex->Fill(mcPart->Zv(),mcRadius);
					   
					   if(TMath::Abs(mcPart->Y())<=0.5)
					   {fHistMcPMLaPt->Fill(mcPart->Pt(),mcPart->M());}
				   }
				   if(antilambdaMC)
				   {
					   if(TMath::Abs(mcPart->Y())<=0.5)
					   {fHistMcPMLbPt->Fill(mcPart->Pt(),mcPart->M());}
				   }
				   if(kshortMC)
				   {
					   if(TMath::Abs(mcPart->Y())<=0.5)
					   {fHistMcPMK0Pt->Fill(mcPart->Pt(),mcPart->M());}
				   }
				}
			else 
				{
					isprimaryMC = false;
					if(lambdaMC)
					{
						fHistMcNLambdaPrimary->Fill(2);
						
						fHistMcSLambdaDaughter->Fill(mcFirstDaughter->PdgCode());
						
						if(motherLabel >=0)
						{
							if(mcMother->PdgCode() != kLambda0)
							{fHistMcLamdaSProductionVertex->Fill(mcPart->Zv(),mcRadius);}
						}
						else
							{fHistMcLamdaSProductionVertex->Fill(mcPart->Zv(),mcRadius);}
						
						if(mcFirstDaughter->PdgCode() == kProton || mcFirstDaughter->PdgCode() == kPiMinus)
							{fHistMcLamdaSDecayVertex->Fill(mcFirstDaughter->Zv(),mcFirstDaughterRadius);}
								
						fHistMcSLambdaEta->Fill(mcPart->Eta());
						fHistMcSLambdaPt->Fill(mcPart->Pt());
					}
				}
			
		}
		

	} 


	fHistMcNLambda->Fill(nLambdaMC);
	fHistMcNAntilambda->Fill(nAntilambdaMC);
	fHistMcNKshort->Fill(nKshortMC);
	
	//END OF MONTE CARLO SECTION
	/*********************************************************************/
	
			
	

	
    // Do some fast cuts first
    // check for good reconstructed vertex
    if(!(fESD->GetPrimaryVertex()->GetStatus())) return;
    // if vertex is from spd vertexZ, require more stringent cut
    if (fESD->GetPrimaryVertex()->IsFromVertexerZ()) {
        if (fESD->GetPrimaryVertex()->GetDispersion()>0.02 ||  fESD->GetPrimaryVertex()->GetZRes()>0.25 ) return; // bad vertex from VertexerZ
    }
     
	Double_t tV0Position[3];
	Double_t tPVPosition[3];
	Double_t radius;
	
	// physics selection
	UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	if(!maskIsSelected)
    {
		//printf("Event failed physics selection\n");
		return;
    }
	
	//if additional initial cuts wanted, can set conditions here
	bool isCut = (fESD->GetNumberOfTracks()==0);
	if (isCut)
	{return;}
	
	//gets primary vertex for the event
	const AliESDVertex *kPv = ((AliESDEvent *)fESD)->GetPrimaryVertex();
	if ( kPv != 0 ) 
    {
		tPVPosition[0] = kPv->GetXv();
		tPVPosition[1] = kPv->GetYv();
		tPVPosition[2] = kPv->GetZv();
		if( tPVPosition[2] == 0. ) 
		{
			//printf("WARNING: Primary vertex a Z = 0, aborting\n");
			return;
		}
    }
	else 
    {
		//printf("ERROR: Primary vertex not found\n");
		return;
    }
	if( !kPv->GetStatus())
    {return;}
	
	

	
	int nLambda(0);
	int nV0s(0);
	double lambdaLabel[fESD->GetNumberOfV0s()];
	
	// V0 loop - runs over every v0 in the event
	for (Int_t iV0 = 0; iV0 < fESD->GetNumberOfV0s(); iV0++) 
    {
		
		AliESDv0 *v0 = fESD->GetV0(iV0);
		if (!v0) 
		{
			printf("ERROR: Could not receive v0 %d\n", iV0);
			continue;
		}
		
		bool lambdaCandidate = true;
		bool antilambdaCandidate = true;
		bool kshortCandidate = true;
		
		// keep only events of interest for fHistMLa plots
		
        if (v0->GetEffMass(4,2) < 1.08 || v0->GetEffMass(4,2) > 1.2 || TMath::Abs(v0->Y(3122))>0.5 )
		{lambdaCandidate = false;}
        if (v0->GetEffMass(2,4) < 1.08 || v0->GetEffMass(2,4) > 1.2 || TMath::Abs(v0->Y(-3122))>0.5)
		{antilambdaCandidate = false;}
        if (v0->GetEffMass(2,2) < 0.414 || v0->GetEffMass(2,2) > 0.582 || TMath::Abs(v0->Y(310))>0.5)
		{kshortCandidate = false;}
		if (v0->GetOnFlyStatus())
		{continue;}
		
		if(!isMonteCarlo) 
		{if(lambdaCandidate == false && antilambdaCandidate == false && kshortCandidate == false)
		{continue;}}
		
		
		//gets details of the v0
		v0->GetXYZ(tV0Position[0],tV0Position[1],tV0Position[2]);
		radius = TMath::Sqrt(tV0Position[0]*tV0Position[0]+tV0Position[1]*tV0Position[1]);
		
		double decayLength = (sqrt((tV0Position[0]-tPVPosition[0])*(tV0Position[0]-tPVPosition[0])+(tV0Position[1]-tPVPosition[1])*(tV0Position[1]-tPVPosition[1])+(tV0Position[2]-tPVPosition[2])*(tV0Position[2]-tPVPosition[2])));
		double cTauLa = decayLength*(v0->GetEffMass(4,2))/(v0->P());
		double cTauLb = decayLength*(v0->GetEffMass(2,4))/(v0->P());
		double cTauK0 = decayLength*(v0->GetEffMass(2,2))/(v0->P());
		
		Int_t indexP, indexN;
		indexP = TMath::Abs(v0->GetPindex());
		AliESDtrack *posTrack = ((AliESDEvent*)fESD)->GetTrack(indexP);
		indexN = TMath::Abs(v0->GetNindex());
		AliESDtrack *negTrack = ((AliESDEvent*)fESD)->GetTrack(indexN);
		
		if(!posTrack || !negTrack)
		{continue;}
		
		double pTrackMomentum[3];
		double nTrackMomentum[3];
		double pV0x, pV0y, pV0z;
		posTrack->GetConstrainedPxPyPz(pTrackMomentum);
		negTrack->GetConstrainedPxPyPz(nTrackMomentum);
		v0->GetPxPyPz(pV0x, pV0y, pV0z);

		const double kMLambda = 1.115;
		const double kMProton = 0.938;
		const double kMPi     = 0.140;
		
		double pPos2 = sqrt(pTrackMomentum[0]*pTrackMomentum[0]+pTrackMomentum[1]*pTrackMomentum[1]+pTrackMomentum[2]*pTrackMomentum[2]);
		double pNeg2 = sqrt(nTrackMomentum[0]*nTrackMomentum[0]+nTrackMomentum[1]*nTrackMomentum[1]+nTrackMomentum[2]*nTrackMomentum[2]);
		double pV02 = sqrt(pV0x*pV0x+pV0y*pV0y+pV0z*pV0z);
		
		//to prevent segfaults when ratios etc taken
		if(pV02 < 0.01 || pPos2 <0.01 || pNeg2 <0.01)
		{continue;}
		
		Float_t pImpactxy(0), pImpactz(0);
		Float_t nImpactxy(0), nImpactz(0);
		posTrack->GetImpactParameters(pImpactxy,pImpactz);
		negTrack->GetImpactParameters(nImpactxy,nImpactz);
		nImpactxy = sqrt((nImpactxy*nImpactxy));
		nImpactz  = sqrt((nImpactz *nImpactz ));
		pImpactxy = sqrt((pImpactxy*pImpactxy));
		pImpactz  = sqrt((pImpactz *pImpactz ));
		
		/*********************************************************************/
		// Cuts are implemented here.
		
		if(!(fTrackCuts->IsSelected(posTrack)) || !(fTrackCuts->IsSelected(negTrack)))
		{
			lambdaCandidate = false;
			antilambdaCandidate = false;
			kshortCandidate = false;
		}
		
		//extra cut to account for difference between p2 & p1 data
		if(nImpactxy < 0.1 || pImpactxy < 0.1)
		{
			lambdaCandidate = false;
			antilambdaCandidate = false;
			kshortCandidate = false;
		}
		
		//psuedorapidity cut
		if(cutEta != -999)
		{
			if(TMath::Abs(posTrack->Eta()) > cutEta || TMath::Abs(negTrack->Eta())  >cutEta)
			{
				lambdaCandidate = false;
				antilambdaCandidate = false;
				kshortCandidate = false;
			}
		}
		
		//pointing angle cut
		if(cutCosPa != -999)
		{
			if (v0->GetV0CosineOfPointingAngle(tPVPosition[0],tPVPosition[1],tPVPosition[2]) < cutCosPa)
			{
				lambdaCandidate = false;
				antilambdaCandidate = false;
				kshortCandidate = false;
			}
		}
		
		//lifetime cut
		if(cutcTau != -999)
		{
			if(cTauLa < cutcTau)
			{
				lambdaCandidate = false;
			}
			if(cTauLb < cutcTau)
			{
				antilambdaCandidate = false;
			}
			if(cTauK0 < cutcTau)
			{
				kshortCandidate = false;
			}
			
		}
		
		// Impact paramater cut (on neg particle)
		if(cutNImpact != -999)
		{
			if(nImpactxy < cutNImpact || nImpactz < cutNImpact)
			{
				lambdaCandidate = false;
			}
			if(pImpactxy < cutNImpact || pImpactz < cutNImpact)
			{
				antilambdaCandidate = false;
			}
		}
		

		// DCA between daughterscut
		if(cutDCA != -999)
		{
			if(v0->GetDcaV0Daughters() > cutDCA)
			{
				lambdaCandidate = false;
				antilambdaCandidate = false;
				kshortCandidate = false;
			}
		}
		
		// Bethe Bloch cut. Made sightly complicated as options for crude cuts still included. Should probably reduce to just 'official' cuts
		if(cutBetheBloch != -999)
		{ 
			if(posTrack->GetTPCsignal() <0 || negTrack->GetTPCsignal()<0)
			{continue;}
			
			if(lambdaCandidate)
			{
				if(cutBetheBloch > 0)
				{
				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kProton)) > cutBetheBloch )
				{lambdaCandidate = false;}
				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kPion)) > cutBetheBloch )
				{lambdaCandidate = false;}
				}
				
				if(cutBetheBloch == -4)  
				{				
				if(isMonteCarlo) 
					{
					double beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+0.9*0.9))),2);
					double gamma2 = 1.0/(1.0-beta2);
					if(posTrack->GetTPCsignal() < (2.0/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2))
					{lambdaCandidate = false;}
					}
				else 
					{ 
					double beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+kMProton*kMProton))),2);
					double gamma2 = 1.0/(1.0-beta2);
					if(posTrack->GetTPCsignal() < (2.3/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2))
					{lambdaCandidate = false;}
					}
				}
				
			}
			
			if(antilambdaCandidate)
			{
				if(cutBetheBloch > 0)
				{
					if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kProton)) > cutBetheBloch )
					{antilambdaCandidate = false;}
					if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kPion)) > cutBetheBloch )
					{antilambdaCandidate = false;}
				}
				
				if(cutBetheBloch == -4)  
				{ 
					if(isMonteCarlo) 
					{
						double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+0.9*0.9))),2);
						double gamma2 = 1.0/(1.0-beta2);
						if(negTrack->GetTPCsignal() < (2.0/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2))
						{antilambdaCandidate = false;}
					}
					else 
					{ 
						double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+0.9*0.9))),2);
						double gamma2 = 1.0/(1.0-beta2);
						if(negTrack->GetTPCsignal() < (2.3/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2))
						{antilambdaCandidate = false;}
					}
				}
				
			}
			
			if(kshortCandidate)
			{
				if(cutBetheBloch > 0)
				{
					if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kPion)) > cutBetheBloch )
					{kshortCandidate = false;}
					if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kPion)) > cutBetheBloch )
					{kshortCandidate = false;}
				}
				
				
				if(cutBetheBloch == -4)  
				{ 
					double par0 = 0.20;
					double par1 = 4.2;
					double par2 = 1000000;
					
					if(isMonteCarlo) 
					{ 
						double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+par0*par0))),2);
						double gamma2 = 1.0/(1.0-beta2);
						if(negTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pNeg2) > -0.6)
						{kshortCandidate = false;}
						
						beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+par0*par0))),2);
						gamma2 = 1.0/(1.0-beta2);
						if(posTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pNeg2) > -0.6)
						{kshortCandidate = false;}
					}
					else 
					{ 
						double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+par0*par0))),2);
						double gamma2 = 1.0/(1.0-beta2);
						if(negTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pNeg2) > -0.6)
						{kshortCandidate = false;}
						
						beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+par0*par0))),2);
						gamma2 = 1.0/(1.0-beta2);
						if(posTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pPos2) > -0.6)
						{kshortCandidate = false;}
					}
				}
				
			}
		}
		
		// Selection of random cuts which I've been playing with
		/*if(nImpactxy > 3 || pImpactxy > 2)
		 {				
		 lambdaCandidate = false;
		 }*/
		
		/*if(nImpactxy < 0.4 || pImpactxy < 0.4 || nImpactxy > 2.5 || pImpactxy >2.5)
		 {				
		 antilambdaCandidate = false;
		 }	*/
		
		//if(decayLength > 25 )
		//{lambdaCandidate = false;}
		
		// Cuts to look at just lowpt region of lambdas
		if(v0->Pt() < 0.3 || v0->Pt() > 0.7 || !lambdaCandidate)
		{continue;}
		
		// cuts to just look at signal/background region of lambda mass
		//if(!((v0->GetEffMass(4,2) >= 1.096 && v0->GetEffMass(4,2) < 1.106) || (v0->GetEffMass(4,2) >= 1.126 && v0->GetEffMass(4,2) < 1.136) ))
		//if(!(v0->GetEffMass(4,2) >= 1.106 && v0->GetEffMass(4,2) < 1.126 ))
		//{continue;}
		/*********************************************************************/

		/*********************************************************************/
		// MONTE CARLO SECTION 2
		// this section looks at the individual V0s
		
		bool mcLambdaCandidate = true;
		bool mcAntilambdaCandidate = true;
		bool mcK0Candidate = true;
		bool isPrimaryMC2 = false;
		bool realParticle = true;
		
		if(isMonteCarlo) 
		{
			
			AliMCEvent *mcEvent = MCEvent();
			AliStack* mcStack = mcEvent->Stack();
			if( !mcStack ) { Printf( "Stack not available"); return; }
			
			TParticle *negParticle = mcStack->Particle( TMath::Abs(negTrack->GetLabel())); 
			TParticle *posParticle = mcStack->Particle( TMath::Abs(posTrack->GetLabel())); 
			
			Int_t negParticleMotherLabel = negParticle->GetFirstMother();
			Int_t posParticleMotherLabel = posParticle->GetFirstMother();
			
			if( negParticleMotherLabel == -1 || posParticleMotherLabel == -1)
			{
			realParticle = false;
			mcLambdaCandidate = false;
			mcAntilambdaCandidate = false;
			mcK0Candidate  =false;
			}

				
			

			
			if( negParticleMotherLabel != posParticleMotherLabel)
			{
				mcLambdaCandidate = false;
				mcAntilambdaCandidate = false;
				mcK0Candidate  =false;
			}
			
			if(realParticle == true)
			{
				AliMCParticle *mcPart2 = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(negParticleMotherLabel));
				fHistMcTPCTrackLength->Fill(mcPart2->GetNumberOfTrackReferences());
				Int_t firstDaughterLabel2 = mcPart2->GetFirstDaughter();
				Int_t lastDaughterLabel2 = mcPart2->GetLastDaughter();
				AliMCParticle *mcFirstDaughter2 = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(firstDaughterLabel2));
				AliMCParticle *mcLastDaughter2 = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(lastDaughterLabel2));
				fHistMcTPCTrackLength->Fill(mcFirstDaughter2->GetNumberOfTrackReferences());
				fHistMcTPCTrackLength->Fill(mcLastDaughter2->GetNumberOfTrackReferences());
				
				if(mcPart2->PdgCode() != kLambda0)
					{mcLambdaCandidate = false;}
				if(mcPart2->PdgCode() != kLambda0Bar)
					{mcAntilambdaCandidate = false;}
				if(mcPart2->PdgCode() != kK0Short)
					{mcK0Candidate  =false;}
				
				if( mcEvent->IsPhysicalPrimary(negParticleMotherLabel) )
					{isPrimaryMC2 = true;}
				
				Int_t motherMotherLabel2 = mcPart2->GetMother();
				AliMCParticle *mcMotherMother2 = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(motherMotherLabel2));
				if(motherMotherLabel2 >= 0)
				{
					if((mcMotherMother2->PdgCode() == 3212 && mcLambdaCandidate) || (mcMotherMother2->PdgCode() == -3212 && mcAntilambdaCandidate))
					{
						if(mcEvent->IsPhysicalPrimary(motherMotherLabel2) )
						{isPrimaryMC2 = true;}
					}
				}
				
				
				
				double mcRadius = TMath::Sqrt((mcPart2->Xv())*(mcPart2->Xv())+(mcPart2->Yv())*(mcPart2->Yv()));
				
				if(mcLambdaCandidate) 
				{
					lambdaLabel[iV0] = negParticleMotherLabel;
					bool once = false;
					for(int j = 0; j < iV0; j++)
					{
						if(lambdaLabel[iV0] == lambdaLabel[j])
						{	
							if(!once)
								{
									fHistMcLog->Fill(6);
									once = true;
								}
							fHistMcLog->Fill(7);
						}	
					}
				}
				
				if(mcLambdaCandidate && lambdaCandidate)
				{
					fHistMcPtV0La->Fill(v0->Pt());
					fHistMcV0MLaPt->Fill(v0->Pt(),v0->GetEffMass(4,2)); 
				}
				if(mcAntilambdaCandidate && antilambdaCandidate)
				{
					fHistMcPtV0Lb->Fill(v0->Pt());
					fHistMcV0MLbPt->Fill(v0->Pt(),v0->GetEffMass(2,4));
				}
				if(mcK0Candidate && kshortCandidate)
				{
					fHistMcPtV0K0->Fill(v0->Pt());
					fHistMcV0MK0Pt->Fill(v0->Pt(),v0->GetEffMass(2,2));
				}
				
				if(mcLambdaCandidate && isPrimaryMC2 && lambdaCandidate)
					{fHistMcLog->Fill(1);}
				if(mcLambdaCandidate && isPrimaryMC2 && !lambdaCandidate)
					{fHistMcLog->Fill(2);}
				if(mcLambdaCandidate && !isPrimaryMC2 && lambdaCandidate)
				{
					fHistMcLog->Fill(3);
					fHistMcV0IDLamdaSProductionVertex->Fill(mcPart2->Zv(),mcRadius);
					fHistMcV0LamdaSProductionVertex->Fill(mcPart2->Zv(),mcRadius);
				}
				if(mcLambdaCandidate && !isPrimaryMC2 && !lambdaCandidate)
				{
					fHistMcLog->Fill(4);
					fHistMcV0LamdaSProductionVertex->Fill(mcPart2->Zv(),mcRadius);
				}

					
			}
			

			
			AliMCParticle *posParticle2 = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(posTrack->GetLabel())));
			AliMCParticle *negParticle2 = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(negTrack->GetLabel())));
			
			if(posParticle2->PdgCode() == kProton)
			{
				fHistMcSigmaPProton->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kProton)));
			}
			if(posParticle2->PdgCode() != kProton)
			{
				fHistMcSigmaPNProton->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kProton)));
			}
		}
	
		// END OF MONTE CARLO SECTION 2
		/*********************************************************************/

		//remove all non-candidates
		if(lambdaCandidate == false && antilambdaCandidate == false && kshortCandidate == false)
		{continue;}
		
		
		//count number of valid v0s
		nV0s+=1;
			
		/*********************************************************************/
		//This section fills histograms
		
		if(lambdaCandidate)
		{
			if((v0->GetEffMass(4,2) >= 1.096 && v0->GetEffMass(4,2) < 1.106) || (v0->GetEffMass(4,2) >= 1.126 && v0->GetEffMass(4,2) < 1.136) )
			{fHistDCAPtPBg->Fill(v0->Pt(),nImpactxy);}
			if(v0->GetEffMass(4,2) >= 1.106 && v0->GetEffMass(4,2) < 1.126 )
			{fHistDCAPtPSig->Fill(v0->Pt(),nImpactxy);}
		}
		
		if(antilambdaCandidate)
		{
			if((v0->GetEffMass(2,4) >= 1.096 && v0->GetEffMass(2,4) < 1.106) || (v0->GetEffMass(2,4) >= 1.126 && v0->GetEffMass(2,4) < 1.136) )
			{fHistDCAPtPbarBg->Fill(v0->Pt(),pImpactxy);}
			if(v0->GetEffMass(2,4) >= 1.106 && v0->GetEffMass(2,4) < 1.126 )
			{fHistDCAPtPbarSig->Fill(v0->Pt(),pImpactxy);}
		}
			
		
		
		fHistEtaPTracks->Fill(posTrack->Eta());
		fHistEtaNTracks->Fill(negTrack->Eta());
		fHistEtaV0s->Fill(v0->Eta());
		
		
		if(lambdaCandidate)
		{
			double laMass = v0->GetEffMass(4,2);
			
			if((laMass > 1.096 && laMass < 1.106) || (laMass > 1.126 && laMass < 1.136) )
			{
			fHistLambdaBgRapidity->Fill(v0->Y(3122));
			fHistLambdaBgEta->Fill(v0->Eta());
			fHistLambdaBgPt->Fill(v0->Pt());
			}
			if(laMass > 1.106 && laMass < 1.126 )
			{
			fHistLambdaSigRapidity->Fill(v0->Y(3122));
			fHistLambdaSigEta->Fill(v0->Eta());
			fHistLambdaSigPt->Fill(v0->Pt());
			}
			
			fHistNTPCPosClustersMLa->Fill((posTrack->GetTPCNcls()),(v0->GetEffMass(4,2)));
			fHistNTPCNegClustersMLa->Fill((negTrack->GetTPCNcls()),(v0->GetEffMass(4,2)));
			fHistPosChi2PerClusterMLa->Fill((posTrack->GetTPCchi2())/(posTrack->GetTPCNcls()),(v0->GetEffMass(4,2)));
			fHistNegChi2PerClusterMLa->Fill((negTrack->GetTPCchi2())/(negTrack->GetTPCNcls()),(v0->GetEffMass(4,2)));
			if((posTrack->GetKinkIndex(0) <= 0) && (negTrack->GetKinkIndex(0) <= 0))
			{fHistKinkIndexFalse->Fill((v0->GetEffMass(4,2)));}
			else
			{fHistKinkIndexTrue->Fill((v0->GetEffMass(4,2)));}
			UInt_t status = posTrack->GetStatus();
			fHistPosTPCRefitMLa->Fill((status&AliESDtrack::kTPCrefit),(v0->GetEffMass(4,2)));
			status = negTrack->GetStatus();
			fHistNegTPCRefitMLa->Fill((status&AliESDtrack::kTPCrefit),(v0->GetEffMass(4,2)));				 
		}
		
		fHistDCAV0Daughters->Fill(v0->GetDcaV0Daughters());
		fHistDecayLDCA->Fill(decayLength,v0->GetDcaV0Daughters());
		
		fHistCosPA->Fill(v0->GetV0CosineOfPointingAngle(tPVPosition[0],tPVPosition[1],tPVPosition[2]));
		fHistDecayL->Fill(decayLength);
		fHistDecayLxy->Fill(sqrt((tV0Position[0]-tPVPosition[0])*(tV0Position[0]-tPVPosition[0])+(tV0Position[1]-tPVPosition[1])*(tV0Position[1]-tPVPosition[1])));
		fHistTauLa->Fill(cTauLa);
		if(lambdaCandidate)
		{
			fHistCosPAMLa->Fill(v0->GetV0CosineOfPointingAngle(tPVPosition[0],tPVPosition[1],tPVPosition[2]),(v0->GetEffMass(4,2)));
			fHistDecayLMLa->Fill(decayLength,(v0->GetEffMass(4,2)));
			fHistTauLaMLa->Fill(cTauLa,(v0->GetEffMass(4,2)));
		}
		

		{fHistBetheBlochTPCPos->Fill(TMath::Log10(pPos2),posTrack->GetTPCsignal());	}
		{fHistBetheBlochTPCNeg->Fill(TMath::Log10(pNeg2),negTrack->GetTPCsignal());}
		
		fHistBetheBlochITSPos->Fill(TMath::Log10(pPos2),posTrack->GetITSsignal());
		fHistBetheBlochITSNeg->Fill(TMath::Log10(pNeg2),negTrack->GetITSsignal());
		
		fHistImpactxyN->Fill(nImpactxy);
		fHistImpactzN->Fill(nImpactz);
		fHistImpactxyP->Fill(pImpactxy);
		fHistImpactzP->Fill(pImpactz);
		
		if(lambdaCandidate)
		{
			fHistImpactxyMLa->Fill(nImpactxy,v0->GetEffMass(4,2));
			fHistImpactzMLa->Fill(nImpactz,v0->GetEffMass(4,2));
		}
		fHistImpactxyImpactz->Fill(nImpactxy,nImpactz);
		
		fHistV0Z->Fill(tV0Position[2]);
		fHistZ->Fill(tV0Position[2]-tPVPosition[2]);
		fHistPtArmR->Fill(radius,v0->PtArmV0());
		
		fHistRZ->Fill(tV0Position[2],radius);
		fHistPtV0Z->Fill(v0->Pt(),tV0Position[2]);
		
		fHistPtArm->Fill(v0->AlphaV0(),v0->PtArmV0());
		fHistXZ->Fill(tV0Position[2],tV0Position[0]);
		fHistYZ->Fill(tV0Position[2],tV0Position[1]);
		fHistPtArmR->Fill(radius,v0->PtArmV0());
		fHistPtV0->Fill(v0->Pt());
		
		//effective mass histograms
		
		//sets assumed particle type of pos/neg daughters. 
		// 0 = electron, 1 = Muon, 2 = pion, 3 = kaon, 4 = proton.
		int dPos = 0;
		int dNeg = 0;
		
		//    v0->ChangeMassHypothesis(kLambda0);
		dPos = 4;
		dNeg = 2;
		if(v0->GetEffMass(dPos,dNeg) > 1.11 && v0->GetEffMass(dPos,dNeg) < 1.13)
		{
			if(!(v0->GetOnFlyStatus()))
			{
				nLambda++;
			}
		}
		if(lambdaCandidate)
		{
			fHistMLa->Fill(v0->GetEffMass(dPos,dNeg));
			fHistMLaR->Fill(radius,v0->GetEffMass(dPos,dNeg));
			fHistMLaPtArm->Fill(v0->PtArmV0(),v0->GetEffMass(dPos,dNeg));   
			fHistMLaPt->Fill(v0->Pt(),v0->GetEffMass(dPos,dNeg));   
		}
		
		//    v0->ChangeMassHypothesis(kK0Short);
		dPos = 2;
		dNeg = 2;
		if(kshortCandidate)
		{
			fHistMK0->Fill(v0->GetEffMass(dPos,dNeg));
			fHistMK0R->Fill(radius,v0->GetEffMass(dPos,dNeg));
			fHistMK0PtArm->Fill(v0->PtArmV0(),v0->GetEffMass(dPos,dNeg));
			fHistMK0Pt->Fill(v0->Pt(),v0->GetEffMass(dPos,dNeg));
		}
		//    v0->ChangeMassHypothesis(kLambda0Bar);
		dPos = 2;
		dNeg = 4;
		if(antilambdaCandidate)
		{
			fHistMLb->Fill(v0->GetEffMass(dPos,dNeg));
			fHistMLbR->Fill(radius,v0->GetEffMass(dPos,dNeg));
			fHistMLbPtArm->Fill(v0->PtArmV0(),v0->GetEffMass(dPos,dNeg));
			fHistMLbPt->Fill(v0->Pt(),v0->GetEffMass(dPos,dNeg));
		}
		
		float binMin = 0.20;
		float binWidth = 0.35;
		
		//   v0->ChangeMassHypothesis(kLambda0);
		
		if(lambdaCandidate)
		{
			dPos = 4;
			dNeg = 2;
			if((v0->Pt()) < binMin || (v0->Pt())>= (binMin + 10*binWidth))
				continue;
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 1*binWidth))    
				fHistMLa0->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin  && (v0->Pt())< (binMin + 2*binWidth))    
				fHistMLa1->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 3*binWidth))    
				fHistMLa2->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 4*binWidth))    
				fHistMLa3->Fill(v0->GetEffMass(dPos,dNeg));    
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 5*binWidth))    
				fHistMLa4->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 6*binWidth))    
				fHistMLa5->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin  && (v0->Pt())< (binMin + 7*binWidth))    
				fHistMLa6->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 8*binWidth))    
				fHistMLa7->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 9*binWidth))    
				fHistMLa8->Fill(v0->GetEffMass(dPos,dNeg));    
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 10*binWidth))    
				fHistMLa9->Fill(v0->GetEffMass(dPos,dNeg));
		}
		//    v0->ChangeMassHypothesis(kLambda0Bar);
		
		if(antilambdaCandidate)
		{
			dPos = 2;
			dNeg = 4;
			if((v0->Pt()) < binMin || (v0->Pt())>= (binMin + 10*binWidth))
				continue;
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 1*binWidth))    
				fHistMLb0->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin  && (v0->Pt())< (binMin + 2*binWidth))    
				fHistMLb1->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 3*binWidth))    
				fHistMLb2->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 4*binWidth))    
				fHistMLb3->Fill(v0->GetEffMass(dPos,dNeg));    
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 5*binWidth))    
				fHistMLb4->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 6*binWidth))    
				fHistMLb5->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin  && (v0->Pt())< (binMin + 7*binWidth))    
				fHistMLb6->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 8*binWidth))    
				fHistMLb7->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 9*binWidth))    
				fHistMLb8->Fill(v0->GetEffMass(dPos,dNeg));    
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 10*binWidth))    
				fHistMLb9->Fill(v0->GetEffMass(dPos,dNeg));
		}
		
		if(kshortCandidate)
		{
			dPos = 2;
			dNeg = 2;
			if((v0->Pt()) < binMin || (v0->Pt())>= (binMin + 10*binWidth))
				continue;
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 1*binWidth))    
				fHistMK00->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin  && (v0->Pt())< (binMin + 2*binWidth))    
				fHistMK01->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 3*binWidth))    
				fHistMK02->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 4*binWidth))    
				fHistMK03->Fill(v0->GetEffMass(dPos,dNeg));    
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 5*binWidth))    
				fHistMK04->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 6*binWidth))    
				fHistMK05->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin  && (v0->Pt())< (binMin + 7*binWidth))    
				fHistMK06->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 8*binWidth))    
				fHistMK07->Fill(v0->GetEffMass(dPos,dNeg));
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 9*binWidth))    
				fHistMK08->Fill(v0->GetEffMass(dPos,dNeg));    
			else if((v0->Pt()) >= binMin && (v0->Pt())< (binMin + 10*binWidth))    
				fHistMK09->Fill(v0->GetEffMass(dPos,dNeg));
		}
		
    }
	
	// track loop, runs over all tracks in event
	int nTracks = 0;
	for (Int_t iTrack=0;  iTrack < fESD->GetNumberOfTracks(); iTrack++)
	{
		
		AliESDtrack* track=fESD->GetTrack(iTrack);
		//Float_t xy=0;
		//Float_t z=0;
		//track->GetImpactParameters(xy,z);
		//&&(xy<1.0)&&(z<1.0)) 
		if ( !fTrackCuts->IsSelected(track)) 
		{continue;}
		if (cutEta != -999 && TMath::Abs(track->Eta()) > cutEta)
		{continue;}
		
		nTracks += 1;
		fHistEta->Fill(track->Eta());
		
	}
	
	fHistPVZ->Fill(tPVPosition[2]);
	fHistNV0->Fill(nV0s);
	fHistNTracks->Fill(nTracks);	
	fHistNV0sNTracks->Fill(nTracks*nTracks,nV0s); 
	fHistdNV0sdT2->Fill( (1.0*nV0s)/(nTracks*nTracks));
	fHistMagneticField->Fill(fESD->GetMagneticField()); 
	fHistNLambda->Fill(nLambda);	
	
	
	fHistLuke->Fill(3); // test histogram
    // NEW HISTO should be filled before this point, as PostData puts the
    // information for this iteration of the UserExec in the container
	PostData(1, fOutputList);
}


//________________________________________________________________________
void AliAnalysisTaskLukeV0::Terminate(Option_t *) 
{
    // Draw result to screen, or perform fitting, normalizations
    // Called once at the end of the query
        
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if(!fOutputList) { Printf("ERROR: could not retrieve TList fOutputList"); return; }
      
 }
