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

/* $Id: AliAnalysisTaskLukeAOD.cxx 46301 2011-01-06 14:25:27Z agheata $ */

/* AliAnalysisTaskLukeAOD.cxx
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 *
 * Instructions for adding histograms can be found below, starting with NEW HISTO
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 * Edited by Luke Hanratty for AODS
 */
#include "AliAnalysisTaskLukeAOD.h"

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
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"

#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"
#include "AliCentrality.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

#include "AliInputEventHandler.h"



ClassImp(AliAnalysisTaskLukeAOD)

//________________________________________________________________________
AliAnalysisTaskLukeAOD::AliAnalysisTaskLukeAOD() // All data members should be initialised here
:AliAnalysisTaskSE(),fIsMonteCarlo(false), fcutCosPa(0.998),fcutcTauMin(-999), fcutNcTauMax(3.0), fcutBetheBloch(3.0), fcutMinNClustersTPC(70), fcutRatio(0.8), fcutEta(0.8), fcutRapidity(0.5), fcutArmenteros(0.2),
fOutput(0),
fPIDResponse(0),

fHistPt(0), 
fHistEta(0),
fHistLog(0),
fHistNV0(0),
fHistZVertex(0),
fHistMCZVertex(0),
fHistCentrality(0),

fHistBBK0Pos(0),
fHistBBK0Neg(0),
fHistBBLaPos(0),
fHistBBLaNeg(0),
fHistBBLbPos(0),
fHistBBLbNeg(0),

fHistBB3SigProton(0),
fHistMK0Pt(0), 
fHistMLaPt(0), 
fHistMLbPt(0),
fHistMcPMK0Pt(0),
fHistMcPMLaPt(0),
fHistMcPMLbPt(0),

fHistMcFMLaPt(0),

fHistMK0PtCent0005(0), 
fHistMLaPtCent0005(0), 
fHistMLbPtCent0005(0),
fHistMcPMK0PtCent0005(0),
fHistMcPMLaPtCent0005(0),
fHistMcPMLbPtCent0005(0),
fHistMcAsMK0PtCent0005(0),
fHistMcAsMLaPtCent0005(0),
fHistMcAsMLbPtCent0005(0),
fHistZVertexCent0005(0),
fHistMCZVertexCent0005(0),

fHistMK0PtCent0510(0), 
fHistMLaPtCent0510(0), 
fHistMLbPtCent0510(0),
fHistMcPMK0PtCent0510(0),
fHistMcPMLaPtCent0510(0),
fHistMcPMLbPtCent0510(0),
fHistMcAsMK0PtCent0510(0),
fHistMcAsMLaPtCent0510(0),
fHistMcAsMLbPtCent0510(0),
fHistZVertexCent0510(0),
fHistMCZVertexCent0510(0),

fHistMK0PtCent1020(0), 
fHistMLaPtCent1020(0), 
fHistMLbPtCent1020(0),
fHistMcPMK0PtCent1020(0),
fHistMcPMLaPtCent1020(0),
fHistMcPMLbPtCent1020(0),
fHistMcAsMK0PtCent1020(0),
fHistMcAsMLaPtCent1020(0),
fHistMcAsMLbPtCent1020(0),
fHistZVertexCent1020(0),
fHistMCZVertexCent1020(0),

fHistMK0PtCent2040(0), 
fHistMLaPtCent2040(0), 
fHistMLbPtCent2040(0),
fHistMcPMK0PtCent2040(0),
fHistMcPMLaPtCent2040(0),
fHistMcPMLbPtCent2040(0),
fHistMcAsMK0PtCent2040(0),
fHistMcAsMLaPtCent2040(0),
fHistMcAsMLbPtCent2040(0),
fHistZVertexCent2040(0),
fHistMCZVertexCent2040(0),

fHistMK0PtCent4060(0), 
fHistMLaPtCent4060(0), 
fHistMLbPtCent4060(0),
fHistMcPMK0PtCent4060(0),
fHistMcPMLaPtCent4060(0),
fHistMcPMLbPtCent4060(0),
fHistMcAsMK0PtCent4060(0),
fHistMcAsMLaPtCent4060(0),
fHistMcAsMLbPtCent4060(0),
fHistZVertexCent4060(0),
fHistMCZVertexCent4060(0),

fHistMK0PtCent6090(0), 
fHistMLaPtCent6090(0), 
fHistMLbPtCent6090(0),
fHistMcPMK0PtCent6090(0),
fHistMcPMLaPtCent6090(0),
fHistMcPMLbPtCent6090(0),
fHistMcAsMK0PtCent6090(0),
fHistMcAsMLaPtCent6090(0),
fHistMcAsMLbPtCent6090(0),
fHistZVertexCent6090(0),
fHistMCZVertexCent6090(0),

fHistMK0PtCent0090(0), 
fHistMLaPtCent0090(0), 
fHistMLbPtCent0090(0),
fHistMcPMK0PtCent0090(0),
fHistMcPMLaPtCent0090(0),
fHistMcPMLbPtCent0090(0),
fHistMcAsMK0PtCent0090(0),
fHistMcAsMLaPtCent0090(0),
fHistMcAsMLbPtCent0090(0),
fHistZVertexCent0090(0),
fHistMCZVertexCent0090(0),

fHistCosPaMLa(0),
fHistCosPaMLb(0),
fHistCosPaMK0(0),	
fHistMcGenCosPaMLa(0),
fHistMcGenCosPaMLb(0),
fHistMcGenCosPaMK0(0),	
fHistMcAsReconCosPaMLa(0),
fHistMcAsReconCosPaMLb(0),
fHistMcAsReconCosPaMK0(0),
fHistMcAsTruthCosPaMLa(0),
fHistMcAsTruthCosPaMLb(0),
fHistMcAsTruthCosPaMK0(0),

fHistcTauMLa(0),	
fHistcTauMLb(0),
fHistcTauMK0(0),		
fHistMcGencTauMLa(0),
fHistMcGencTauMLb(0),
fHistMcGencTauMK0(0),	
fHistMcAsReconcTauMLa(0),
fHistMcAsReconcTauMLb(0),	
fHistMcAsReconcTauMK0(0),
fHistMcAsTruthcTauMLa(0),
fHistMcAsTruthcTauMLb(0),	
fHistMcAsTruthcTauMK0(0),

fHistDcaMLa(0),	
fHistDcaMLb(0),	
fHistDcaMK0(0),	
fHistMcGenDcaMLa(0),
fHistMcGenDcaMLb(0),
fHistMcGenDcaMK0(0),	
fHistMcAsReconDcaMLa(0),	
fHistMcAsReconDcaMLb(0),	
fHistMcAsReconDcaMK0(0),
fHistMcAsTruthDcaMLa(0),	
fHistMcAsTruthDcaMLb(0),	
fHistMcAsTruthDcaMK0(0),

fHistNSigmaMLa(0),	
fHistNSigmaMLb(0),		
fHistNSigmaMK0(0),		
fHistMcGenNSigmaMLa(0),	
fHistMcGenNSigmaMLb(0),	
fHistMcGenNSigmaMK0(0),	
fHistMcAsReconNSigmaMLa(0),	
fHistMcAsReconNSigmaMLb(0),	
fHistMcAsReconNSigmaMK0(0),
fHistMcAsTruthNSigmaMLa(0),	
fHistMcAsTruthNSigmaMLb(0),	
fHistMcAsTruthNSigmaMK0(0),

fHistEtaMLa(0),	
fHistEtaMLb(0),	
fHistEtaMK0(0),		
fHistMcGenEtaMLa(0),
fHistMcGenEtaMLb(0),	
fHistMcGenEtaMK0(0),
fHistMcAsReconEtaMLa(0),	
fHistMcAsReconEtaMLb(0),	
fHistMcAsReconEtaMK0(0),
fHistMcAsTruthEtaMLa(0),	
fHistMcAsTruthEtaMLb(0),	
fHistMcAsTruthEtaMK0(0),

fHistRapMLa(0),
fHistRapMLb(0),		
fHistRapMK0(0),	
fHistMcGenRapMLa(0),	
fHistMcGenRapMLb(0),	
fHistMcGenRapMK0(0),	
fHistMcAsReconRapMLa(0),
fHistMcAsReconRapMLb(0),
fHistMcAsReconRapMK0(0),
fHistMcAsTruthRapMLa(0),
fHistMcAsTruthRapMLb(0),
fHistMcAsTruthRapMK0(0),

fHistArmPodK0(0),
fHistArmPodLa(0),
fHistArmPodLb(0),
fHistMcGenArmPodK0(0),
fHistMcGenArmPodLa(0),
fHistMcGenArmPodLb(0),
fHistMcAsReconArmPodK0(0),
fHistMcAsReconArmPodLa(0),
fHistMcAsReconArmPodLb(0),
fHistMcAsTruthArmPodK0(0),
fHistMcAsTruthArmPodLa(0),
fHistMcAsTruthArmPodLb(0)

// The last in the above list should not have a comma after it
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskLukeAOD::AliAnalysisTaskLukeAOD(const char *name) // All data members should be initialised here
:AliAnalysisTaskSE(name), fIsMonteCarlo(false), fcutCosPa(0.998),fcutcTauMin(-999), fcutNcTauMax(3.0), fcutBetheBloch(3.0), fcutMinNClustersTPC(70), fcutRatio(0.8), fcutEta(0.8), fcutRapidity(0.5), fcutArmenteros(0.2),
fOutput(0),
fPIDResponse(0),

fHistPt(0), 
fHistEta(0),
fHistLog(0),
fHistNV0(0),
fHistZVertex(0),
fHistMCZVertex(0),
fHistCentrality(0),

fHistBBK0Pos(0),
fHistBBK0Neg(0),
fHistBBLaPos(0),
fHistBBLaNeg(0),
fHistBBLbPos(0),
fHistBBLbNeg(0),

fHistBB3SigProton(0),
fHistMK0Pt(0), 
fHistMLaPt(0), 
fHistMLbPt(0),
fHistMcPMK0Pt(0),
fHistMcPMLaPt(0),
fHistMcPMLbPt(0),

fHistMcFMLaPt(0),

fHistMK0PtCent0005(0), 
fHistMLaPtCent0005(0), 
fHistMLbPtCent0005(0),
fHistMcPMK0PtCent0005(0),
fHistMcPMLaPtCent0005(0),
fHistMcPMLbPtCent0005(0),
fHistMcAsMK0PtCent0005(0),
fHistMcAsMLaPtCent0005(0),
fHistMcAsMLbPtCent0005(0),
fHistZVertexCent0005(0),
fHistMCZVertexCent0005(0),

fHistMK0PtCent0510(0), 
fHistMLaPtCent0510(0), 
fHistMLbPtCent0510(0),
fHistMcPMK0PtCent0510(0),
fHistMcPMLaPtCent0510(0),
fHistMcPMLbPtCent0510(0),
fHistMcAsMK0PtCent0510(0),
fHistMcAsMLaPtCent0510(0),
fHistMcAsMLbPtCent0510(0),
fHistZVertexCent0510(0),
fHistMCZVertexCent0510(0),

fHistMK0PtCent1020(0), 
fHistMLaPtCent1020(0), 
fHistMLbPtCent1020(0),
fHistMcPMK0PtCent1020(0),
fHistMcPMLaPtCent1020(0),
fHistMcPMLbPtCent1020(0),
fHistMcAsMK0PtCent1020(0),
fHistMcAsMLaPtCent1020(0),
fHistMcAsMLbPtCent1020(0),
fHistZVertexCent1020(0),
fHistMCZVertexCent1020(0),

fHistMK0PtCent2040(0), 
fHistMLaPtCent2040(0), 
fHistMLbPtCent2040(0),
fHistMcPMK0PtCent2040(0),
fHistMcPMLaPtCent2040(0),
fHistMcPMLbPtCent2040(0),
fHistMcAsMK0PtCent2040(0),
fHistMcAsMLaPtCent2040(0),
fHistMcAsMLbPtCent2040(0),
fHistZVertexCent2040(0),
fHistMCZVertexCent2040(0),

fHistMK0PtCent4060(0), 
fHistMLaPtCent4060(0), 
fHistMLbPtCent4060(0),
fHistMcPMK0PtCent4060(0),
fHistMcPMLaPtCent4060(0),
fHistMcPMLbPtCent4060(0),
fHistMcAsMK0PtCent4060(0),
fHistMcAsMLaPtCent4060(0),
fHistMcAsMLbPtCent4060(0),
fHistZVertexCent4060(0),
fHistMCZVertexCent4060(0),

fHistMK0PtCent6090(0), 
fHistMLaPtCent6090(0), 
fHistMLbPtCent6090(0),
fHistMcPMK0PtCent6090(0),
fHistMcPMLaPtCent6090(0),
fHistMcPMLbPtCent6090(0),
fHistMcAsMK0PtCent6090(0),
fHistMcAsMLaPtCent6090(0),
fHistMcAsMLbPtCent6090(0),
fHistZVertexCent6090(0),
fHistMCZVertexCent6090(0),

fHistMK0PtCent0090(0), 
fHistMLaPtCent0090(0), 
fHistMLbPtCent0090(0),
fHistMcPMK0PtCent0090(0),
fHistMcPMLaPtCent0090(0),
fHistMcPMLbPtCent0090(0),
fHistMcAsMK0PtCent0090(0),
fHistMcAsMLaPtCent0090(0),
fHistMcAsMLbPtCent0090(0),
fHistZVertexCent0090(0),
fHistMCZVertexCent0090(0),

fHistCosPaMLa(0),
fHistCosPaMLb(0),
fHistCosPaMK0(0),	
fHistMcGenCosPaMLa(0),
fHistMcGenCosPaMLb(0),
fHistMcGenCosPaMK0(0),	
fHistMcAsReconCosPaMLa(0),
fHistMcAsReconCosPaMLb(0),
fHistMcAsReconCosPaMK0(0),
fHistMcAsTruthCosPaMLa(0),
fHistMcAsTruthCosPaMLb(0),
fHistMcAsTruthCosPaMK0(0),

fHistcTauMLa(0),	
fHistcTauMLb(0),
fHistcTauMK0(0),		
fHistMcGencTauMLa(0),
fHistMcGencTauMLb(0),
fHistMcGencTauMK0(0),	
fHistMcAsReconcTauMLa(0),
fHistMcAsReconcTauMLb(0),	
fHistMcAsReconcTauMK0(0),
fHistMcAsTruthcTauMLa(0),
fHistMcAsTruthcTauMLb(0),	
fHistMcAsTruthcTauMK0(0),

fHistDcaMLa(0),	
fHistDcaMLb(0),	
fHistDcaMK0(0),	
fHistMcGenDcaMLa(0),
fHistMcGenDcaMLb(0),
fHistMcGenDcaMK0(0),	
fHistMcAsReconDcaMLa(0),	
fHistMcAsReconDcaMLb(0),	
fHistMcAsReconDcaMK0(0),
fHistMcAsTruthDcaMLa(0),	
fHistMcAsTruthDcaMLb(0),	
fHistMcAsTruthDcaMK0(0),

fHistNSigmaMLa(0),	
fHistNSigmaMLb(0),		
fHistNSigmaMK0(0),		
fHistMcGenNSigmaMLa(0),	
fHistMcGenNSigmaMLb(0),	
fHistMcGenNSigmaMK0(0),	
fHistMcAsReconNSigmaMLa(0),	
fHistMcAsReconNSigmaMLb(0),	
fHistMcAsReconNSigmaMK0(0),
fHistMcAsTruthNSigmaMLa(0),	
fHistMcAsTruthNSigmaMLb(0),	
fHistMcAsTruthNSigmaMK0(0),

fHistEtaMLa(0),	
fHistEtaMLb(0),	
fHistEtaMK0(0),		
fHistMcGenEtaMLa(0),
fHistMcGenEtaMLb(0),	
fHistMcGenEtaMK0(0),
fHistMcAsReconEtaMLa(0),	
fHistMcAsReconEtaMLb(0),	
fHistMcAsReconEtaMK0(0),
fHistMcAsTruthEtaMLa(0),	
fHistMcAsTruthEtaMLb(0),	
fHistMcAsTruthEtaMK0(0),

fHistRapMLa(0),
fHistRapMLb(0),		
fHistRapMK0(0),	
fHistMcGenRapMLa(0),	
fHistMcGenRapMLb(0),	
fHistMcGenRapMK0(0),	
fHistMcAsReconRapMLa(0),
fHistMcAsReconRapMLb(0),
fHistMcAsReconRapMK0(0),
fHistMcAsTruthRapMLa(0),
fHistMcAsTruthRapMLb(0),
fHistMcAsTruthRapMK0(0),

fHistArmPodK0(0),
fHistArmPodLa(0),
fHistArmPodLb(0),
fHistMcGenArmPodK0(0),
fHistMcGenArmPodLa(0),
fHistMcGenArmPodLb(0),
fHistMcAsReconArmPodK0(0),
fHistMcAsReconArmPodLa(0),
fHistMcAsReconArmPodLb(0),
fHistMcAsTruthArmPodK0(0),
fHistMcAsTruthArmPodLa(0),
fHistMcAsTruthArmPodLb(0)

// The last in the above list should not have a comma after it
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskLukeAOD::~AliAnalysisTaskLukeAOD()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fOutput;
    }
}

//________________________________________________________________________
void AliAnalysisTaskLukeAOD::UserCreateOutputObjects()
{
    // Create histograms
    // Called once (on the worker node)
	
    fOutput = new TList();
    fOutput->SetOwner();  // IMPORTANT!
	
	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
	fPIDResponse = inputHandler->GetPIDResponse();
	//Bool_t maskIsSelected = inputHandler->IsEventSelected();
    
	// Create histograms
    Int_t ptbins = 15;
    Float_t ptlow = 0.1, ptup = 3.1;
    fHistPt = new TH1F("fHistPt", "P_{T} distribution for reconstructed", ptbins, ptlow, ptup);
    fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPt->SetMarkerStyle(kFullCircle);
	
    Int_t etabins = 40;
    Float_t etalow = -2.0, etaup = 2.0;
    fHistEta = new TH1F("fHistEta","#eta distribution for reconstructed",etabins, etalow, etaup);
    fHistEta->GetXaxis()->SetTitle("#eta");
    fHistEta->GetYaxis()->SetTitle("counts");
	
	fHistLog = new TH1F("fHistLog","Log Variables",100, -0.5, 99.5);
	fHistNV0 = new TH1F("fHistNV0","Number of V0s per event",100, 0, 5000);
	fHistZVertex = new TH1F("fHistZVertex","Z coordinate of primary vertex",60, -15, 15);
	fHistMCZVertex = new TH1F("fHistMCZVertex","Z coordinate of primary vertex in MC truth",60, -15, 15);
	fHistCentrality = new	TH1F("fHistCentrality", "Centrality of Events; Centrality Percentile (%); N",110,-5,105);
	
	fHistBBK0Pos = new	TH2F("fHistBBK0Pos","PID of the positive daughter of K0 candidates; P (GeV); -dE/dx (keV/cm ?)",1000,0,10,1000,0,200);
	fHistBBK0Neg = new	TH2F("fHistBBK0Neg","PID of the negative daughter of K0 candidates; P (GeV); -dE/dx (keV/cm ?)",1000,0,10,1000,0,200);
	fHistBBLaPos = new	TH2F("fHistBBLaPos","PID of the positive daughter of La candidates; P (GeV); -dE/dx (keV/cm ?)",1000,0,10,1000,0,200);
	fHistBBLaNeg = new	TH2F("fHistBBLaNeg","PID of the negative daughter of La candidates; P (GeV); -dE/dx (keV/cm ?)",1000,0,10,1000,0,200);
	fHistBBLbPos = new	TH2F("fHistBBLbPos","PID of the positive daughter of Lb candidates; P (GeV); -dE/dx (keV/cm ?)",1000,0,10,1000,0,200);
	fHistBBLbNeg = new	TH2F("fHistBBLbNeg","PID of the negative daughter of Lb candidates; P (GeV); -dE/dx (keV/cm ?)",1000,0,10,1000,0,200);
	
	fHistBB3SigProton = new	TH2F("fHistBB3SigProton","-dE/dX against Momentum for Protons @3sigma from TPC; P (GeV); -dE/dx (keV/cm ?)",1000,0,10,1000,0,200);
	fHistMK0Pt = new	TH2F("fHistMK0Pt","K0 Mass versus Pt; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPt = new	TH2F("fHistMLaPt","Lambda Mass versus Pt; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPt = new	TH2F("fHistMLbPt","AntiLambda Mass versus Pt; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0Pt = new	TH2F("fHistMcPMK0Pt","Monte Carlo primary K0 Mass versus Pt; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPt = new	TH2F("fHistMcPMLaPt","Monte Carlo primary (& sigma0) Lambda Mass versus Pt; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPt = new	TH2F("fHistMcPMLbPt","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	
	fHistMcFMLaPt = new	TH2F("fHistMcFMLaPt","Monte Carlo Reconstructed Lambda Mass versus Pt for feeddown lambda; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	
	fHistMK0PtCent0005 = new	TH2F("fHistMK0PtCent0005","K0 Mass versus Pt for Centrality 0-5%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent0005 = new	TH2F("fHistMLaPtCent0005","Lambda Mass versus Pt for Centrality 0-5%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent0005 = new	TH2F("fHistMLbPtCent0005","AntiLambda Mass versus Pt for Centrality 0-5%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent0005 = new	TH2F("fHistMcPMK0PtCent0005","Monte Carlo primary K0 Mass versus Pt for Centrality 0-5%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent0005 = new	TH2F("fHistMcPMLaPtCent0005","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 0-5%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent0005 = new	TH2F("fHistMcPMLbPtCent0005","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 0-5%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMK0PtCent0005 = new	TH2F("fHistMcAsMK0PtCent0005","Monte Carlo associated K0 Mass versus Pt for Centrality 0-5%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcAsMLaPtCent0005 = new	TH2F("fHistMcAsMLaPtCent0005","Monte Carlo associated (& sigma0) Lambda Mass versus Pt for Centrality 0-5%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMLbPtCent0005 = new	TH2F("fHistMcAsMLbPtCent0005","Monte Carlo associated (& sigma0) AntiLambda Mass versus Pt for Centrality 0-5%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent0005 = new TH1F("fHistZVertexCent0005","Z coordinate of primary vertex for Centrality 0-5%",60, -15, 15);
	fHistMCZVertexCent0005 = new TH1F("fHistMCZVertexCent0005","Z coordinate of primary vertex in MC truth for Centrality 0-5%",60, -15, 15);
	
	fHistMK0PtCent0510 = new	TH2F("fHistMK0PtCent0510","K0 Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent0510 = new	TH2F("fHistMLaPtCent0510","Lambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent0510 = new	TH2F("fHistMLbPtCent0510","AntiLambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent0510 = new	TH2F("fHistMcPMK0PtCent0510","Monte Carlo primary K0 Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent0510 = new	TH2F("fHistMcPMLaPtCent0510","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent0510 = new	TH2F("fHistMcPMLbPtCent0510","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMK0PtCent0510 = new	TH2F("fHistMcAsMK0PtCent0510","Monte Carlo associated K0 Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcAsMLaPtCent0510 = new	TH2F("fHistMcAsMLaPtCent0510","Monte Carlo associated (& sigma0) Lambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMLbPtCent0510 = new	TH2F("fHistMcAsMLbPtCent0510","Monte Carlo associated (& sigma0) AntiLambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent0510 = new TH1F("fHistZVertexCent0510","Z coordinate of primary vertex for Centrality 5-10%",60, -15, 15);
	fHistMCZVertexCent0510 = new TH1F("fHistMCZVertexCent0510","Z coordinate of primary vertex in MC truth for Centrality 5-10%",60, -15, 15);
	
	fHistMK0PtCent1020 = new	TH2F("fHistMK0PtCent1020","K0 Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent1020 = new	TH2F("fHistMLaPtCent1020","Lambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent1020 = new	TH2F("fHistMLbPtCent1020","AntiLambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent1020 = new	TH2F("fHistMcPMK0PtCent1020","Monte Carlo primary K0 Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent1020 = new	TH2F("fHistMcPMLaPtCent1020","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent1020 = new	TH2F("fHistMcPMLbPtCent1020","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMK0PtCent1020 = new	TH2F("fHistMcAsMK0PtCent1020","Monte Carlo associated K0 Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcAsMLaPtCent1020 = new	TH2F("fHistMcAsMLaPtCent1020","Monte Carlo associated (& sigma0) Lambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMLbPtCent1020 = new	TH2F("fHistMcAsMLbPtCent1020","Monte Carlo associated (& sigma0) AntiLambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent1020 = new TH1F("fHistZVertexCent1020","Z coordinate of primary vertex for Centrality 10-20%",60, -15, 15);
	fHistMCZVertexCent1020 = new TH1F("fHistMCZVertexCent1020","Z coordinate of primary vertex in MC truth for Centrality 10-20%",60, -15, 15);
	
	fHistMK0PtCent2040 = new	TH2F("fHistMK0PtCent2040","K0 Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent2040 = new	TH2F("fHistMLaPtCent2040","Lambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent2040 = new	TH2F("fHistMLbPtCent2040","AntiLambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent2040 = new	TH2F("fHistMcPMK0PtCent2040","Monte Carlo primary K0 Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent2040 = new	TH2F("fHistMcPMLaPtCent2040","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent2040 = new	TH2F("fHistMcPMLbPtCent2040","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMK0PtCent2040 = new	TH2F("fHistMcAsMK0PtCent2040","Monte Carlo associated K0 Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcAsMLaPtCent2040 = new	TH2F("fHistMcAsMLaPtCent2040","Monte Carlo associated (& sigma0) Lambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMLbPtCent2040 = new	TH2F("fHistMcAsMLbPtCent2040","Monte Carlo associated (& sigma0) AntiLambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent2040 = new TH1F("fHistZVertexCent2040","Z coordinate of primary vertex for Centrality 20-40%",60, -15, 15);
	fHistMCZVertexCent2040 = new TH1F("fHistMCZVertexCent2040","Z coordinate of primary vertex in MC truth for Centrality 20-40%",60, -15, 15);
	
	fHistMK0PtCent4060 = new	TH2F("fHistMK0PtCent4060","K0 Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent4060 = new	TH2F("fHistMLaPtCent4060","Lambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent4060 = new	TH2F("fHistMLbPtCent4060","AntiLambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent4060 = new	TH2F("fHistMcPMK0PtCent4060","Monte Carlo primary K0 Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent4060 = new	TH2F("fHistMcPMLaPtCent4060","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent4060 = new	TH2F("fHistMcPMLbPtCent4060","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMK0PtCent4060 = new	TH2F("fHistMcAsMK0PtCent4060","Monte Carlo associated K0 Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcAsMLaPtCent4060 = new	TH2F("fHistMcAsMLaPtCent4060","Monte Carlo associated (& sigma0) Lambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMLbPtCent4060 = new	TH2F("fHistMcAsMLbPtCent4060","Monte Carlo associated (& sigma0) AntiLambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent4060 = new TH1F("fHistZVertexCent4060","Z coordinate of primary vertex for Centrality 40-60%",60, -15, 15);
	fHistMCZVertexCent4060 = new TH1F("fHistMCZVertexCent4060","Z coordinate of primary vertex in MC truth for Centrality 40-60%",60, -15, 15);
	
	fHistMK0PtCent6090 = new	TH2F("fHistMK0PtCent6090","K0 Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent6090 = new	TH2F("fHistMLaPtCent6090","Lambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent6090 = new	TH2F("fHistMLbPtCent6090","AntiLambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent6090 = new	TH2F("fHistMcPMK0PtCent6090","Monte Carlo primary K0 Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent6090 = new	TH2F("fHistMcPMLaPtCent6090","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent6090 = new	TH2F("fHistMcPMLbPtCent6090","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMK0PtCent6090 = new	TH2F("fHistMcAsMK0PtCent6090","Monte Carlo associated K0 Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcAsMLaPtCent6090 = new	TH2F("fHistMcAsMLaPtCent6090","Monte Carlo associated (& sigma0) Lambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMLbPtCent6090 = new	TH2F("fHistMcAsMLbPtCent6090","Monte Carlo associated (& sigma0) AntiLambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent6090 = new TH1F("fHistZVertexCent6090","Z coordinate of primary vertex for Centrality 60-90%",60, -15, 15);
	fHistMCZVertexCent6090 = new TH1F("fHistMCZVertexCent6090","Z coordinate of primary vertex in MC truth for Centrality 60-90%",60, -15, 15);
	
	fHistMK0PtCent0090 = new	TH2F("fHistMK0PtCent0090","K0 Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent0090 = new	TH2F("fHistMLaPtCent0090","Lambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent0090 = new	TH2F("fHistMLbPtCent0090","AntiLambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent0090 = new	TH2F("fHistMcPMK0PtCent0090","Monte Carlo primary K0 Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent0090 = new	TH2F("fHistMcPMLaPtCent0090","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent0090 = new	TH2F("fHistMcPMLbPtCent0090","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMK0PtCent0090 = new	TH2F("fHistMcAsMK0PtCent0090","Monte Carlo associated K0 Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcAsMLaPtCent0090 = new	TH2F("fHistMcAsMLaPtCent0090","Monte Carlo associated (& sigma0) Lambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcAsMLbPtCent0090 = new	TH2F("fHistMcAsMLbPtCent0090","Monte Carlo associated (& sigma0) AntiLambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent0090 = new TH1F("fHistZVertexCent0090","Z coordinate of primary vertex for Centrality 0-90%",60, -15, 15);
	fHistMCZVertexCent0090 = new TH1F("fHistMCZVertexCent0090","Z coordinate of primary vertex in MC truth for Centrality 0-90%",60, -15, 15);
	
	
	fHistCosPaMLa		 = new	TH2F("fHistCosPaMLa","	Reconstructed Mass vs CosPa for Lambda Candidates;CosPA; M(p#pi^{-}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistCosPaMLb		 = new	TH2F("fHistCosPaMLb","	Reconstructed Mass vs CosPa for AntiLambda Candidates;CosPA; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistCosPaMK0		 = new	TH2F("fHistCosPaMK0","	Reconstructed Mass vs CosPa for K0Short Candidates;CosPA; K0 Mass (GeV/c^2)",200,0.99,1.001,140,0.414,0.582);
	fHistMcGenCosPaMLa	 = new	TH2F("fHistMcGenCosPaMLa","	Reconstructed Mass vs MC-Truth CosPa for all MC primary Lambda;CosPA; M(p#pi^{-}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcGenCosPaMLb	 = new	TH2F("fHistMcGenCosPaMLb","	Reconstructed Mass vs MC-Truth CosPa for all MC primary AntiLambda;CosPA; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcGenCosPaMK0	 = new	TH2F("fHistMcGenCosPaMK0","	Reconstructed Mass vs MC-Truth CosPa for all MC primary K0Short;CosPA; K0 Mass (GeV/c^2)",200,0.99,1.001,140,0.414,0.582);
	fHistMcAsReconCosPaMLa	 = new	TH2F("fHistMcAsReconCosPaMLa","	Reconstructed Mass vs CosPa for reconstructed MC primary Lambda;CosPA; M(p#pi^{-}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcAsReconCosPaMLb	 = new	TH2F("fHistMcAsReconCosPaMLb","	Reconstructed Mass vs CosPa for reconstructed MC primary AntiLambda;CosPA; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcAsReconCosPaMK0 = new	TH2F("fHistMcAsReconCosPaMK0","	Reconstructed Mass vs CosPa for reconstructed MC primary K0Short;CosPA; K0 Mass (GeV/c^2)",200,0.99,1.001,140,0.414,0.582);
	fHistMcAsTruthCosPaMLa	 = new	TH2F("fHistMcAsTruthCosPaMLa","	Reconstructed Mass vs MC-Truth CosPa for reconstructed MC primary Lambda;CosPA; M(p#pi^{-}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcAsTruthCosPaMLb	 = new	TH2F("fHistMcAsTruthCosPaMLb","	Reconstructed Mass vs MC-Truth CosPa for reconstructed MC primary AntiLambda;CosPA; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcAsTruthCosPaMK0 = new	TH2F("fHistMcAsTruthCosPaMK0","	Reconstructed Mass vs MC-Truth CosPa for reconstructed MC primary K0Short;CosPA; K0 Mass (GeV/c^2)",200,0.99,1.001,140,0.414,0.582);
	
	fHistcTauMLa		 = new	TH2F("fHistcTauMLa","	Reconstructed Mass vs cTau for Lambda Candidates; cTau; M(p#pi^{-}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistcTauMLb		 = new	TH2F("fHistcTauMLb","	Reconstructed Mass vs cTau for AntiLambda Candidates; cTau; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistcTauMK0		 = new	TH2F("fHistcTauMK0","	Reconstructed Mass vs cTau for K0Short Candidates; cTau; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	fHistMcGencTauMLa	 = new	TH2F("fHistMcGencTauMLa","	Reconstructed Mass vs MC-Truth cTau for all MC primary Lambda; cTau; M(p#pi^{-}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcGencTauMLb	 = new	TH2F("fHistMcGencTauMLb","	Reconstructed Mass vs MC-Truth cTau for all MC primary AntiLambda; cTau; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcGencTauMK0	 = new	TH2F("fHistMcGencTauMK0","	Reconstructed Mass vs MC-Truth cTau for all MC primary K0Short; cTau; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	fHistMcAsReconcTauMLa	 = new	TH2F("fHistMcAsReconcTauMLa","	Reconstructed Mass vs cTau for reconstructed MC primary Lambda; cTau; M(p#pi^{-}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcAsReconcTauMLb	 = new	TH2F("fHistMcAsReconcTauMLb","	Reconstructed Mass vs cTau for reconstructed MC primary AntiLambda; cTau; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcAsReconcTauMK0 = new	TH2F("fHistMcAsReconcTauMK0","	Reconstructed Mass vs cTau for reconstructed MC primary K0Short; cTau; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	fHistMcAsTruthcTauMLa	 = new	TH2F("fHistMcAsTruthcTauMLa","	Reconstructed Mass vs MC-Truth cTau for reconstructed MC primary Lambda; cTau; M(p#pi^{-}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcAsTruthcTauMLb	 = new	TH2F("fHistMcAsTruthcTauMLb","	Reconstructed Mass vs MC-Truth cTau for reconstructed MC primary AntiLambda; cTau; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcAsTruthcTauMK0 = new	TH2F("fHistMcAsTruthcTauMK0","	Reconstructed Mass vs MC-Truth cTau for reconstructed MC primary K0Short; cTau; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	
	fHistDcaMLa		 = new	TH2F("fHistDcaMLa","	Reconstructed Mass vs Dca for Lambda Candidates; DCA; M(p#pi^{-}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistDcaMLb		 = new	TH2F("fHistDcaMLb","	Reconstructed Mass vs Dca for AntiLambda Candidates; DCA; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistDcaMK0		 = new	TH2F("fHistDcaMK0","	Reconstructed Mass vs Dca for K0Short Candidates; DCA; K0 Mass (GeV/c^2)",200,0,1.5,140,0.414,0.582);
	fHistMcGenDcaMLa	 = new	TH2F("fHistMcGenDcaMLa","	Reconstructed Mass vs MC-Truth Dca for all MC primary Lambda; DCA; M(p#pi^{-}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcGenDcaMLb	 = new	TH2F("fHistMcGenDcaMLb","	Reconstructed Mass vs MC-Truth Dca for all MC primary AntiLambda; DCA; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcGenDcaMK0	 = new	TH2F("fHistMcGenDcaMK0","	Reconstructed Mass vs MC-Truth Dca for all MC primary K0Short; DCA; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	fHistMcAsReconDcaMLa	 = new	TH2F("fHistMcAsReconDcaMLa","	Reconstructed Mass vs Dca for reconstructed MC primary Lambda; DCA; M(p#pi^{-}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistMcAsReconDcaMLb	 = new	TH2F("fHistMcAsReconDcaMLb","	Reconstructed Mass vs Dca for reconstructed MC primary AntiLambda; DCA; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistMcAsReconDcaMK0 = new	TH2F("fHistMcAsReconDcaMK0","	Reconstructed Mass vs Dca for reconstructed MC primary K0Short; DCA; K0 Mass (GeV/c^2)",200,0,1.5,140,0.414,0.582);
	fHistMcAsTruthDcaMLa	 = new	TH2F("fHistMcAsTruthDcaMLa","	Reconstructed Mass vs MC-Truth Dca for reconstructed MC primary Lambda; DCA; M(p#pi^{-}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistMcAsTruthDcaMLb	 = new	TH2F("fHistMcAsTruthDcaMLb","	Reconstructed Mass vs MC-Truth Dca for reconstructed MC primary AntiLambda; DCA; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistMcAsTruthDcaMK0 = new	TH2F("fHistMcAsTruthDcaMK0","	Reconstructed Mass vs MC-Truth Dca for reconstructed MC primary K0Short; DCA; K0 Mass (GeV/c^2)",200,0,1.5,140,0.414,0.582);
	
	fHistNSigmaMLa		 = new	TH2F("fHistNSigmaMLa","	Reconstructed Mass vs NSigma for Lambda Candidates; NSigma; M(p#pi^{-}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistNSigmaMLb		 = new	TH2F("fHistNSigmaMLb","	Reconstructed Mass vs NSigma for AntiLambda Candidates; NSigma; M(#bar{p}#pi^{+}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistNSigmaMK0		 = new	TH2F("fHistNSigmaMK0","	Reconstructed Mass vs NSigma for K0Short Candidates; NSigma; K0 Mass (GeV/c^2)",50,-5,5,140,0.414,0.582);
	fHistMcGenNSigmaMLa	 = new	TH2F("fHistMcGenNSigmaMLa","	Reconstructed Mass vs MC-Truth NSigma for all MC primary Lambda; NSigma; M(p#pi^{-}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcGenNSigmaMLb	 = new	TH2F("fHistMcGenNSigmaMLb","	Reconstructed Mass vs MC-Truth NSigma for all MC primary AntiLambda; NSigma; M(#bar{p}#pi^{+}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcGenNSigmaMK0	 = new	TH2F("fHistMcGenNSigmaMK0","	Reconstructed Mass vs MC-Truth NSigma for all MC primary K0Short; NSigma; K0 Mass (GeV/c^2)",50,-5,5,140,0.414,0.582);
	fHistMcAsReconNSigmaMLa	 = new	TH2F("fHistMcAsReconNSigmaMLa","	Reconstructed Mass vs NSigma for reconstructed MC primary Lambda; NSigma; M(p#pi^{-}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcAsReconNSigmaMLb	 = new	TH2F("fHistMcAsReconNSigmaMLb","	Reconstructed Mass vs NSigma for reconstructed MC primary AntiLambda; NSigma; M(#bar{p}#pi^{+}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcAsReconNSigmaMK0 = new	TH2F("fHistMcAsReconNSigmaMK0","	Reconstructed Mass vs NSigma for reconstructed MC primary K0Short; NSigma; K0 Mass (GeV/c^2)",50,-5,5,140,0.414,0.582);
	fHistMcAsTruthNSigmaMLa	 = new	TH2F("fHistMcAsTruthNSigmaMLa","	Reconstructed Mass vs MC-Truth NSigma for reconstructed MC primary Lambda; NSigma; M(p#pi^{-}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcAsTruthNSigmaMLb	 = new	TH2F("fHistMcAsTruthNSigmaMLb","	Reconstructed Mass vs MC-Truth NSigma for reconstructed MC primary AntiLambda; NSigma; M(#bar{p}#pi^{+}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcAsTruthNSigmaMK0 = new	TH2F("fHistMcAsTruthNSigmaMK0","	Reconstructed Mass vs MC-Truth NSigma for reconstructed MC primary K0Short; NSigma; K0 Mass (GeV/c^2)",50,-5,5,140,0.414,0.582);
	
	fHistEtaMLa		 = new	TH2F("fHistEtaMLa","	Reconstructed Mass vs Eta for Lambda Candidates; Eta; M(p#pi^{-}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistEtaMLb		 = new	TH2F("fHistEtaMLb","	Reconstructed Mass vs Eta for AntiLambda Candidates; Eta; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistEtaMK0		 = new	TH2F("fHistEtaMK0","	Reconstructed Mass vs Eta for K0Short Candidates; Eta; K0 Mass (GeV/c^2)",200,-3,3,140,0.414,0.582);
	fHistMcGenEtaMLa	 = new	TH2F("fHistMcGenEtaMLa","	Reconstructed Mass vs MC-Truth Eta for all MC primary Lambda; Eta; M(p#pi^{-}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcGenEtaMLb	 = new	TH2F("fHistMcGenEtaMLb","	Reconstructed Mass vs MC-Truth Eta for all MC primary AntiLambda; Eta; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcGenEtaMK0	 = new	TH2F("fHistMcGenEtaMK0","	Reconstructed Mass vs MC-Truth Eta for all MC primary K0Short; Eta; K0 Mass (GeV/c^2)",200,-3,3,140,0.414,0.582);
	fHistMcAsReconEtaMLa	 = new	TH2F("fHistMcAsReconEtaMLa","	Reconstructed Mass vs Eta for reconstructed MC primary Lambda; Eta; M(p#pi^{-}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcAsReconEtaMLb	 = new	TH2F("fHistMcAsReconEtaMLb","	Reconstructed Mass vs Eta for reconstructed MC primary AntiLambda; Eta; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcAsReconEtaMK0 = new	TH2F("fHistMcAsReconEtaMK0","	Reconstructed Mass vs Eta for reconstructed MC primary K0Short; Eta; K0 Mass (GeV/c^2)",200,-3,3,140,0.414,0.582);
	fHistMcAsTruthEtaMLa	 = new	TH2F("fHistMcAsTruthEtaMLa","	Reconstructed Mass vs MC-Truth Eta for reconstructed MC primary Lambda; Eta; M(p#pi^{-}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcAsTruthEtaMLb	 = new	TH2F("fHistMcAsTruthEtaMLb","	Reconstructed Mass vs MC-Truth Eta for reconstructed MC primary AntiLambda; Eta; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcAsTruthEtaMK0 = new	TH2F("fHistMcAsTruthEtaMK0","	Reconstructed Mass vs MC-Truth Eta for reconstructed MC primary K0Short; Eta; K0 Mass (GeV/c^2)",200,-3,3,140,0.414,0.582);
	
	fHistRapMLa		 = new	TH2F("fHistRapMLa","	Reconstructed Mass vs Rap for Lambda Candidates; Rapidity; M(p#pi^{-}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistRapMLb		 = new	TH2F("fHistRapMLb","	Reconstructed Mass vs Rap for AntiLambda Candidates; Rapidity; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistRapMK0		 = new	TH2F("fHistRapMK0","	Reconstructed Mass vs Rap for K0Short Candidates; Rapidity; K0 Mass (GeV/c^2)",200,-1,1,140,0.414,0.582);
	fHistMcGenRapMLa	 = new	TH2F("fHistMcGenRapMLa","	Reconstructed Mass vs MC-Truth Rap for all MC primary Lambda; Rapidity; M(p#pi^{-}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcGenRapMLb	 = new	TH2F("fHistMcGenRapMLb","	Reconstructed Mass vs MC-Truth Rap for all MC primary AntiLambda; Rapidity; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcGenRapMK0	 = new	TH2F("fHistMcGenRapMK0","	Reconstructed Mass vs MC-Truth Rap for all MC primary K0Short; Rapidity; K0 Mass (GeV/c^2)",200,-1,1,140,0.414,0.582);
	fHistMcAsReconRapMLa	 = new	TH2F("fHistMcAsReconRapMLa","	Reconstructed Mass vs Rap for reconstructed MC primary Lambda; Rapidity; M(p#pi^{-}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcAsReconRapMLb	 = new	TH2F("fHistMcAsReconRapMLb","	Reconstructed Mass vs Rap for reconstructed MC primary AntiLambda; Rapidity; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcAsReconRapMK0 = new	TH2F("fHistMcAsReconRapMK0","	Reconstructed Mass vs Rap for reconstructed MC primary K0Short; Rapidity; K0 Mass (GeV/c^2)",200,-1,1,140,0.414,0.582);
	fHistMcAsTruthRapMLa	 = new	TH2F("fHistMcAsTruthRapMLa","	Reconstructed Mass vs MC-Truth Rap for reconstructed MC primary Lambda; Rapidity; M(p#pi^{-}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcAsTruthRapMLb	 = new	TH2F("fHistMcAsTruthRapMLb","	Reconstructed Mass vs MC-Truth Rap for reconstructed MC primary AntiLambda; Rapidity; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcAsTruthRapMK0 = new	TH2F("fHistMcAsTruthRapMK0","	Reconstructed Mass vs MC-Truth Rap for reconstructed MC primary K0Short; Rapidity; K0 Mass (GeV/c^2)",200,-1,1,140,0.414,0.582);
	
	fHistArmPodK0 = new	TH2F("fHistArmPodK0","Armenteros plot for K0 candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistArmPodLa = new	TH2F("fHistArmPodLa","Armenteros plot for La candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistArmPodLb = new	TH2F("fHistArmPodLb","Armenteros plot for Lb candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistMcGenArmPodK0 = new	TH2F("fHistMcGenArmPodK0","Armenteros plot for MC Generated K0 candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistMcGenArmPodLa = new	TH2F("fHistMcGenArmPodLa","Armenteros plot for MC Generated La candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistMcGenArmPodLb = new	TH2F("fHistMcGenArmPodLb","Armenteros plot for MC Generated Lb candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistMcAsReconArmPodK0 = new	TH2F("fHistMcAsReconArmPodK0","Armenteros plot for MC Associated K0 candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistMcAsReconArmPodLa = new	TH2F("fHistMcAsReconArmPodLa","Armenteros plot for MC Associated La candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistMcAsReconArmPodLb = new	TH2F("fHistMcAsReconArmPodLb","Armenteros plot for MC Associated Lb candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistMcAsTruthArmPodK0 = new	TH2F("fHistMcAsTruthArmPodK0","True Armenteros plot for MC Associated K0 candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistMcAsTruthArmPodLa = new	TH2F("fHistMcAsTruthArmPodLa","True Armenteros plot for MC Associated La candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistMcAsTruthArmPodLb = new	TH2F("fHistMcAsTruthArmPodLb","True Armenteros plot for MC Associated Lb candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	
	
	
    // NEW HISTO should be defined here, with a sensible name,
	
    fOutput->Add(fHistPt);
    fOutput->Add(fHistEta);
    fOutput->Add(fHistLog);
	fOutput->Add(fHistNV0);
	fOutput->Add(fHistZVertex);
	fOutput->Add(fHistMCZVertex);
	fOutput->Add(fHistCentrality);
	
	fOutput->Add(fHistBBK0Pos);
	fOutput->Add(fHistBBK0Neg);
	fOutput->Add(fHistBBLaPos);
	fOutput->Add(fHistBBLaNeg);
	fOutput->Add(fHistBBLbPos);
	fOutput->Add(fHistBBLbNeg);
	
	fOutput->Add(fHistBB3SigProton);
	fOutput->Add(fHistMK0Pt); 
	fOutput->Add(fHistMLaPt); 
	fOutput->Add(fHistMLbPt);
	fOutput->Add(fHistMcPMK0Pt);
	fOutput->Add(fHistMcPMLaPt);
	fOutput->Add(fHistMcPMLbPt);
	
	fOutput->Add(fHistMcFMLaPt);
	
	fOutput->Add(fHistMK0PtCent0005); 
	fOutput->Add(fHistMLaPtCent0005); 
	fOutput->Add(fHistMLbPtCent0005);
	fOutput->Add(fHistMcPMK0PtCent0005);
	fOutput->Add(fHistMcPMLaPtCent0005);
	fOutput->Add(fHistMcPMLbPtCent0005);
	fOutput->Add(fHistMcAsMK0PtCent0005);
	fOutput->Add(fHistMcAsMLaPtCent0005);
	fOutput->Add(fHistMcAsMLbPtCent0005);
	fOutput->Add(fHistZVertexCent0005);
	fOutput->Add(fHistMCZVertexCent0005);
	
	fOutput->Add(fHistMK0PtCent0510); 
	fOutput->Add(fHistMLaPtCent0510); 
	fOutput->Add(fHistMLbPtCent0510);
	fOutput->Add(fHistMcPMK0PtCent0510);
	fOutput->Add(fHistMcPMLaPtCent0510);
	fOutput->Add(fHistMcPMLbPtCent0510);
	fOutput->Add(fHistMcAsMK0PtCent0510);
	fOutput->Add(fHistMcAsMLaPtCent0510);
	fOutput->Add(fHistMcAsMLbPtCent0510);
	fOutput->Add(fHistZVertexCent0510);
	fOutput->Add(fHistMCZVertexCent0510);
	
	fOutput->Add(fHistMK0PtCent1020); 
	fOutput->Add(fHistMLaPtCent1020); 
	fOutput->Add(fHistMLbPtCent1020);
	fOutput->Add(fHistMcPMK0PtCent1020);
	fOutput->Add(fHistMcPMLaPtCent1020);
	fOutput->Add(fHistMcPMLbPtCent1020);
	fOutput->Add(fHistMcAsMK0PtCent1020);
	fOutput->Add(fHistMcAsMLaPtCent1020);
	fOutput->Add(fHistMcAsMLbPtCent1020);
	fOutput->Add(fHistZVertexCent1020);
	fOutput->Add(fHistMCZVertexCent1020);
	
	fOutput->Add(fHistMK0PtCent2040); 
	fOutput->Add(fHistMLaPtCent2040); 
	fOutput->Add(fHistMLbPtCent2040);
	fOutput->Add(fHistMcPMK0PtCent2040);
	fOutput->Add(fHistMcPMLaPtCent2040);
	fOutput->Add(fHistMcPMLbPtCent2040);
	fOutput->Add(fHistMcAsMK0PtCent2040);
	fOutput->Add(fHistMcAsMLaPtCent2040);
	fOutput->Add(fHistMcAsMLbPtCent2040);
	fOutput->Add(fHistZVertexCent2040);
	fOutput->Add(fHistMCZVertexCent2040);
	
	fOutput->Add(fHistMK0PtCent4060); 
	fOutput->Add(fHistMLaPtCent4060); 
	fOutput->Add(fHistMLbPtCent4060);
	fOutput->Add(fHistMcPMK0PtCent4060);
	fOutput->Add(fHistMcPMLaPtCent4060);
	fOutput->Add(fHistMcPMLbPtCent4060);
	fOutput->Add(fHistMcAsMK0PtCent4060);
	fOutput->Add(fHistMcAsMLaPtCent4060);
	fOutput->Add(fHistMcAsMLbPtCent4060);
	fOutput->Add(fHistZVertexCent4060);
	fOutput->Add(fHistMCZVertexCent4060);
	
	fOutput->Add(fHistMK0PtCent6090); 
	fOutput->Add(fHistMLaPtCent6090); 
	fOutput->Add(fHistMLbPtCent6090);
	fOutput->Add(fHistMcPMK0PtCent6090);
	fOutput->Add(fHistMcPMLaPtCent6090);
	fOutput->Add(fHistMcPMLbPtCent6090);
	fOutput->Add(fHistMcAsMK0PtCent6090);
	fOutput->Add(fHistMcAsMLaPtCent6090);
	fOutput->Add(fHistMcAsMLbPtCent6090);
	fOutput->Add(fHistZVertexCent6090);
	fOutput->Add(fHistMCZVertexCent6090);
	
	fOutput->Add(fHistMK0PtCent0090); 
	fOutput->Add(fHistMLaPtCent0090); 
	fOutput->Add(fHistMLbPtCent0090);
	fOutput->Add(fHistMcPMK0PtCent0090);
	fOutput->Add(fHistMcPMLaPtCent0090);
	fOutput->Add(fHistMcPMLbPtCent0090);
	fOutput->Add(fHistMcAsMK0PtCent0090);
	fOutput->Add(fHistMcAsMLaPtCent0090);
	fOutput->Add(fHistMcAsMLbPtCent0090);
	fOutput->Add(fHistZVertexCent0090);
	fOutput->Add(fHistMCZVertexCent0090);
	
	fOutput->Add(fHistCosPaMLa);
	fOutput->Add(fHistCosPaMLb);
	fOutput->Add(fHistCosPaMK0);	
	fOutput->Add(fHistMcGenCosPaMLa);
	fOutput->Add(fHistMcGenCosPaMLb);
	fOutput->Add(fHistMcGenCosPaMK0);	
	fOutput->Add(fHistMcAsReconCosPaMLa);
	fOutput->Add(fHistMcAsReconCosPaMLb);
	fOutput->Add(fHistMcAsReconCosPaMK0);
	fOutput->Add(fHistMcAsTruthCosPaMLa);
	fOutput->Add(fHistMcAsTruthCosPaMLb);
	fOutput->Add(fHistMcAsTruthCosPaMK0);
	
	
	fOutput->Add(fHistcTauMLa);	
	fOutput->Add(fHistcTauMLb);
	fOutput->Add(fHistcTauMK0);		
	fOutput->Add(fHistMcGencTauMLa);
	fOutput->Add(fHistMcGencTauMLb);
	fOutput->Add(fHistMcGencTauMK0);	
	fOutput->Add(fHistMcAsReconcTauMLa);
	fOutput->Add(fHistMcAsReconcTauMLb);	
	fOutput->Add(fHistMcAsReconcTauMK0);
	fOutput->Add(fHistMcAsTruthcTauMLa);
	fOutput->Add(fHistMcAsTruthcTauMLb);	
	fOutput->Add(fHistMcAsTruthcTauMK0);
	
	fOutput->Add(fHistDcaMLa);	
	fOutput->Add(fHistDcaMLb);	
	fOutput->Add(fHistDcaMK0);	
	fOutput->Add(fHistMcGenDcaMLa);
	fOutput->Add(fHistMcGenDcaMLb);
	fOutput->Add(fHistMcGenDcaMK0);	
	fOutput->Add(fHistMcAsReconDcaMLa);	
	fOutput->Add(fHistMcAsReconDcaMLb);	
	fOutput->Add(fHistMcAsReconDcaMK0);
	fOutput->Add(fHistMcAsTruthDcaMLa);	
	fOutput->Add(fHistMcAsTruthDcaMLb);	
	fOutput->Add(fHistMcAsTruthDcaMK0);
	
	fOutput->Add(fHistNSigmaMLa);	
	fOutput->Add(fHistNSigmaMLb);		
	fOutput->Add(fHistNSigmaMK0);		
	fOutput->Add(fHistMcGenNSigmaMLa);	
	fOutput->Add(fHistMcGenNSigmaMLb);	
	fOutput->Add(fHistMcGenNSigmaMK0);	
	fOutput->Add(fHistMcAsReconNSigmaMLa);	
	fOutput->Add(fHistMcAsReconNSigmaMLb);	
	fOutput->Add(fHistMcAsReconNSigmaMK0);
	fOutput->Add(fHistMcAsTruthNSigmaMLa);	
	fOutput->Add(fHistMcAsTruthNSigmaMLb);	
	fOutput->Add(fHistMcAsTruthNSigmaMK0);
	
	fOutput->Add(fHistEtaMLa);	
	fOutput->Add(fHistEtaMLb);	
	fOutput->Add(fHistEtaMK0);		
	fOutput->Add(fHistMcGenEtaMLa);
	fOutput->Add(fHistMcGenEtaMLb);	
	fOutput->Add(fHistMcGenEtaMK0);
	fOutput->Add(fHistMcAsReconEtaMLa);	
	fOutput->Add(fHistMcAsReconEtaMLb);	
	fOutput->Add(fHistMcAsReconEtaMK0);
	fOutput->Add(fHistMcAsTruthEtaMLa);	
	fOutput->Add(fHistMcAsTruthEtaMLb);	
	fOutput->Add(fHistMcAsTruthEtaMK0);
	
	fOutput->Add(fHistRapMLa);
	fOutput->Add(fHistRapMLb);		
	fOutput->Add(fHistRapMK0);	
	fOutput->Add(fHistMcGenRapMLa);	
	fOutput->Add(fHistMcGenRapMLb);	
	fOutput->Add(fHistMcGenRapMK0);	
	fOutput->Add(fHistMcAsReconRapMLa);
	fOutput->Add(fHistMcAsReconRapMLb);
	fOutput->Add(fHistMcAsReconRapMK0);
	fOutput->Add(fHistMcAsTruthRapMLa);
	fOutput->Add(fHistMcAsTruthRapMLb);
	fOutput->Add(fHistMcAsTruthRapMK0);
	
	fOutput->Add(fHistArmPodK0);
	fOutput->Add(fHistArmPodLa);
	fOutput->Add(fHistArmPodLb);
	fOutput->Add(fHistMcGenArmPodK0);
	fOutput->Add(fHistMcGenArmPodLa);
	fOutput->Add(fHistMcGenArmPodLb);
	fOutput->Add(fHistMcAsReconArmPodK0);
	fOutput->Add(fHistMcAsReconArmPodLa);
	fOutput->Add(fHistMcAsReconArmPodLb);
	fOutput->Add(fHistMcAsTruthArmPodK0);
	fOutput->Add(fHistMcAsTruthArmPodLa);
	fOutput->Add(fHistMcAsTruthArmPodLb);
    
	// NEW HISTO added to fOutput here
    PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
/*
static Bool_t checkPrimaryStatus(const AliAODMCParticle *mcPart, double pVx, double pVy, double pVz) 
{
	double dx = pVx - mcPart->Xv();
	double dy = pVy - mcPart->Yv();
	double dz = pVz - mcPart->Zv();
	
	double prodVertex = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
	if(prodVertex <= 0.001) {return true;}
	
	return false;
}
*/

//________________________________________________________________________

static Bool_t AcceptTrack(const AliAODTrack *t, double fcutMinNClustersTPC, double fcutRatio) 
{
	if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
	//if (t->GetKinkIndex(0)>0) return kFALSE;
	
	Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
	if (nCrossedRowsTPC < fcutMinNClustersTPC && fcutMinNClustersTPC != -999) return kFALSE;
	Int_t findable=t->GetTPCNclsF();
	if (findable <= 0) return kFALSE;
	if (nCrossedRowsTPC/findable < fcutRatio && fcutRatio != -999) return kFALSE;
	
	return kTRUE;
}

//________________________________________________________________________
static Bool_t AcceptV0_general(const AliAODv0 *v1, const AliAODEvent *aod, double fcutCosPa, double fcutNImpact, double fcutDCA, double fcutEta, double fcutMinNClustersTPC, double fcutRatio) 
{
	
	if (v1->GetOnFlyStatus()) return kFALSE;
	
	int nnum = 1, pnum = 0;	
	const AliAODTrack *ntracktest=(AliAODTrack *)v1->GetDaughterLabel(nnum);
	if(ntracktest->Charge() > 0){nnum = 0; pnum = 1;}	
	
	const AliAODTrack *ntrack1=(AliAODTrack *)v1->GetDaughterLabel(nnum);
	if (!AcceptTrack(ntrack1, fcutMinNClustersTPC, fcutRatio)) return kFALSE;
	
	const AliAODTrack *ptrack1=(AliAODTrack *)v1->GetDaughterLabel(pnum);
	if (!AcceptTrack(ptrack1, fcutMinNClustersTPC, fcutRatio)) return kFALSE;
	
	Float_t impact=v1->DcaNegToPrimVertex();
	if (TMath::Abs(impact)<0.1) return kFALSE;
	if (TMath::Abs(impact)<fcutNImpact && fcutNImpact != -999) return kFALSE;
	impact=v1->DcaPosToPrimVertex();
	if (TMath::Abs(impact)<0.1) return kFALSE;
	if (TMath::Abs(impact)<fcutNImpact && fcutNImpact != -999) return kFALSE;
	
	Double_t dca=v1->DcaV0Daughters();
	if (TMath::Abs(dca)>fcutDCA && fcutDCA != -999) return kFALSE;
	
	Double_t cpa=v1->CosPointingAngle(aod->GetPrimaryVertex());
	if (cpa<fcutCosPa && fcutCosPa != -999) return kFALSE;
	
	Double_t etaN = v1->PseudoRapNeg();
	Double_t etaP = v1->PseudoRapPos();
	if ((TMath::Abs(etaN)>fcutEta || TMath::Abs(etaP)>fcutEta) && fcutEta != -999) return kFALSE;
	
	return kTRUE;
}

//________________________________________________________________________
static Bool_t AcceptV0_particle(const AliAODv0 *v1, int type,  double fcutcTauMin, double fcutRapidity, Double_t decayL, double fcutNcTauMax) 
{
	
	Double_t cTau = 0;
	if(type == 1)
	{cTau = decayL*(v1->MassLambda())/(v1->P());}
	if(type == 2)
	{cTau = decayL*(v1->MassAntiLambda())/(v1->P());}
	if(type == 0)
	{cTau = decayL*(v1->MassK0Short())/(v1->P());}
	
	if (cTau < fcutcTauMin && fcutcTauMin != -999 ) return kFALSE;
	
	if (cTau > (fcutNcTauMax*7.89) && fcutNcTauMax != -999 && (type ==1 || type ==2)) return kFALSE;
	if (cTau > (fcutNcTauMax*2.68) && fcutNcTauMax != -999 && (type ==0)) return kFALSE;
	
	
	Double_t rap = 0;
	if(type == 1 || type == 2)
	{rap = v1->RapLambda();}
	if(type == 0)
	{rap = v1->RapK0Short();}
	if (TMath::Abs(rap)>fcutRapidity && fcutRapidity != -999) return kFALSE;
	
	return kTRUE;
}

//________________________________________________________________________
static Bool_t AcceptV0_lowpt(const AliAODv0 *v1, AliPIDResponse *PIDResponse,int type, double fcutBetheBloch, bool fIsMonteCarlo) 
{
	
	int nnum = 1, pnum = 0;
	const AliAODTrack *ntracktest=(AliAODTrack *)v1->GetDaughterLabel(nnum);
	if(ntracktest->Charge() > 0){nnum = 0; pnum = 1;}	
	
	const AliAODTrack *ntrack1=(AliAODTrack *)v1->GetDaughterLabel(nnum);
	const AliAODTrack *ptrack1=(AliAODTrack *)v1->GetDaughterLabel(pnum);
	
	Double_t nsig_p = 0;
	Double_t nsig_n = 0;
	
	const AliAODPid *pid_p=ptrack1->GetDetPid();
	const AliAODPid *pid_n=ntrack1->GetDetPid();
	
	if (pid_p) 
	{
		if(type == 1)
		{
			nsig_p=PIDResponse->NumberOfSigmasTPC(ptrack1,AliPID::kProton);
			if (TMath::Abs(nsig_p) > fcutBetheBloch && fcutBetheBloch >0 && !fIsMonteCarlo  && ptrack1->P() <= 1) return kFALSE;
		}
		
		if(type == 2)
		{
			nsig_p=PIDResponse->NumberOfSigmasTPC(ptrack1,AliPID::kProton);
			if (TMath::Abs(nsig_p) <= fcutBetheBloch && fcutBetheBloch >0 && !fIsMonteCarlo  && ptrack1->P() <= 1) return kFALSE;
		}
		
	}
	
	if (pid_n) 
	{
		if(type == 2)
		{
			nsig_n=PIDResponse->NumberOfSigmasTPC(ntrack1,AliPID::kProton);
			if (TMath::Abs(nsig_n) > fcutBetheBloch && fcutBetheBloch >0 && !fIsMonteCarlo  && ntrack1->P() <= 1) return kFALSE;
		}
		
		if(type == 1)
		{
			nsig_n=PIDResponse->NumberOfSigmasTPC(ntrack1,AliPID::kProton);
			if (TMath::Abs(nsig_n) <= fcutBetheBloch && fcutBetheBloch >0 && !fIsMonteCarlo  && ntrack1->P() <= 1) return kFALSE;
		}
		
		
	}
	
	return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskLukeAOD::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
	
	fHistLog->Fill(1);
	
	// parameters used for most fcuts, to minimise editing
	//double fcutCosPa(0.998), fcutcTauMin(-999), fcutNcTauMax(3.0);
	double fcutNImpact(-999), fcutDCA(-999);
	//double fcutBetheBloch(3.0); // NOTE - BB fcut only applies to data, must be accounted for when constructing corrected yields
	//double fcutMinNClustersTPC(70), fcutRatio(0.8);
	//bool fIsMonteCarlo(false); 
	int isMCtype(0);    //1 = Pure Hijing only, 0 = Anything, -1 = Injected only
	//double fcutEta(0.8), fcutRapidity(0.5), fcutArmenteros(0.2);
	
	Double_t mcXv=0., mcYv=0., mcZv=0.;
	
    // Create pointer to reconstructed event
	AliAODEvent *aod=(AliAODEvent *)InputEvent();	
	if (!aod) 
	{
		Printf("ERROR: aod not available");
		fHistLog->Fill(98);
		return;
	}
	
	
	//Bool_t isSelected = (maskIsSelected && AliVEvent::kMB);
	/*if (!isSelected) 
	 {
	 Printf("ERROR: failed physics selection");
	 fHistLog->Fill(88);
	 return;
	 }*/
	
	AliCentrality* centrality = aod->GetCentrality();
	Double_t centPercentile = centrality->GetCentralityPercentile("V0M");
	fHistCentrality->Fill(centPercentile);
	
	/********MC Loop************************************/
	
	int isInjected = -1;
	
	if(fIsMonteCarlo) 
	{
		
		// Create pointer to reconstructed event
		TList *list = aod->GetList();
		TClonesArray *stack = 0x0;
		stack = (TClonesArray*)list->FindObject(AliAODMCParticle::StdBranchName());
		if (!stack) {
			Printf("ERROR: stack not available");
			fHistLog->Fill(84);	
		}
		else{
			
			AliAODMCHeader *mcHdr=(AliAODMCHeader*)list->FindObject(AliAODMCHeader::StdBranchName());
			
			mcXv=mcHdr->GetVtxX(); mcYv=mcHdr->GetVtxY(); mcZv=mcHdr->GetVtxZ();
			if (TMath::Abs(mcZv) > 10.)
			{
				fHistLog->Fill(76);
			}  
			else{
				fHistMCZVertex->Fill(mcZv);
				
				if(centPercentile >= 0.0001 && centPercentile <= 5.0)
				{fHistMCZVertexCent0005->Fill(mcZv);}
				if(centPercentile > 5.0 && centPercentile <= 10.0)
				{fHistMCZVertexCent0510->Fill(mcZv);}
				if(centPercentile > 10.0 && centPercentile <= 20.0)
				{fHistMCZVertexCent1020->Fill(mcZv);}
				if(centPercentile > 20.0 && centPercentile <= 40.0)
				{fHistMCZVertexCent2040->Fill(mcZv);}
				if(centPercentile > 40.0 && centPercentile <= 60.0)
				{fHistMCZVertexCent4060->Fill(mcZv);}
				if(centPercentile > 60.0 && centPercentile <= 90.0)
				{fHistMCZVertexCent6090->Fill(mcZv);}
				if(centPercentile >= 0.0001 && centPercentile <= 90.0)
				{fHistMCZVertexCent0090->Fill(mcZv);}
				
				
				Int_t ntrk=stack->GetEntriesFast(), ntrk0=ntrk;
				
				for(int iMCtrack = 0; iMCtrack<ntrk0; iMCtrack++)
					
				{	
					
					//booleans to check if track is La, Lb, K0 and primary
					bool lambdaMC = false;
					bool antilambdaMC = false;
					bool kshortMC = false;
					bool isprimaryMC = false;
					
					AliAODMCParticle *mcPart =(AliAODMCParticle*)stack->UncheckedAt(iMCtrack);
					
					if ((mcPart->GetStatus() == 21) || (mcPart->GetPdgCode() == 443 && mcPart->GetMother() == -1)) 
					{
						isInjected = iMCtrack;
					}
					
					if(isInjected >= 0 && isMCtype == 1)
					{continue;}
					if(isInjected < 0 && isMCtype == -1)
					{continue;}
					
					Int_t code=mcPart->GetPdgCode();
					if (code != kK0Short && code != kLambda0 && code != kLambda0Bar ) 
					{continue;}
					
					if(code == kLambda0)
					{
						lambdaMC = true;
					}
					else if(code == kK0Short)
					{
						kshortMC = true;
					}
					else if(code == kLambda0Bar)
					{
						antilambdaMC = true;
					}
					
					
					Int_t motherLabel = mcPart->GetMother();
					AliAODMCParticle *mcMother = (AliAODMCParticle *)stack->UncheckedAt(motherLabel);
					Int_t motherType = -1;
					if(motherLabel >= 0)
					{motherType = mcMother->GetPdgCode();}
					
					// this block of code is used to include primary  Sigma0 decays as primary lambda/antilambda
					bool sigma0MC = false;
					if(motherType == 3212 || motherType == -3212)// || motherType == 3322 || motherType == -3322 || motherType == 3312 || motherType == -3312)
					{
						if(mcMother->IsPhysicalPrimary() && (lambdaMC || antilambdaMC))
							//if(checkPrimaryStatus(mcMother, mcXv, mcYv, mcZv))
						{sigma0MC = true;}
					}
					
					
					if(mcPart->IsPhysicalPrimary() || sigma0MC)
						//if(checkPrimaryStatus(mcPart, mcXv, mcYv, mcZv) || sigma0MC)
					{
						isprimaryMC = true;
						if(lambdaMC)
						{
							if(TMath::Abs(mcPart->Y())<=fcutRapidity)
							{
								fHistMcPMLaPt->Fill(mcPart->Pt(),mcPart->M());
								
								if(centPercentile >= 0.0001 && centPercentile <= 5.0)
								{fHistMcPMLaPtCent0005->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 5.0 && centPercentile <= 10.0)
								{fHistMcPMLaPtCent0510->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 10.0 && centPercentile <= 20.0)
								{fHistMcPMLaPtCent1020->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 20.0 && centPercentile <= 40.0)
								{fHistMcPMLaPtCent2040->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 40.0 && centPercentile <= 60.0)
								{fHistMcPMLaPtCent4060->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 60.0 && centPercentile <= 90.0)
								{fHistMcPMLaPtCent6090->Fill(mcPart->Pt(),mcPart->M());}
								
								if(centPercentile >= 0.0001 && centPercentile <= 90.0)
								{fHistMcPMLaPtCent0090->Fill(mcPart->Pt(),mcPart->M());}
							}
						}
						if(antilambdaMC)
						{
							if(TMath::Abs(mcPart->Y())<=fcutRapidity)
							{
								fHistMcPMLbPt->Fill(mcPart->Pt(),mcPart->M());
								
								if(centPercentile >= 0.0001 && centPercentile <= 5.0)
								{fHistMcPMLbPtCent0005->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 5.0 && centPercentile <= 10.0)
								{fHistMcPMLbPtCent0510->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 10.0 && centPercentile <= 20.0)
								{fHistMcPMLbPtCent1020->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 20.0 && centPercentile <= 40.0)
								{fHistMcPMLbPtCent2040->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 40.0 && centPercentile <= 60.0)
								{fHistMcPMLbPtCent4060->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 60.0 && centPercentile <= 90.0)
								{fHistMcPMLbPtCent6090->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile >= 0.0001 && centPercentile <= 90.0)
								{fHistMcPMLbPtCent0090->Fill(mcPart->Pt(),mcPart->M());}
								
							}
						}
						if(kshortMC)
						{
							if(TMath::Abs(mcPart->Y())<=fcutRapidity)
							{
								fHistMcPMK0Pt->Fill(mcPart->Pt(),mcPart->M());
								
								if(centPercentile >= 0.0001 && centPercentile <= 5.0)
								{fHistMcPMK0PtCent0005->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 5.0 && centPercentile <= 10.0)
								{fHistMcPMK0PtCent0510->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 10.0 && centPercentile <= 20.0)
								{fHistMcPMK0PtCent1020->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 20.0 && centPercentile <= 40.0)
								{fHistMcPMK0PtCent2040->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 40.0 && centPercentile <= 60.0)
								{fHistMcPMK0PtCent4060->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile > 60.0 && centPercentile <= 90.0)
								{fHistMcPMK0PtCent6090->Fill(mcPart->Pt(),mcPart->M());}
								if(centPercentile >= 0.0001 && centPercentile <= 90.0)
								{fHistMcPMK0PtCent0090->Fill(mcPart->Pt(),mcPart->M());}
								
							}
						}
					}
					
					
					Int_t daughter0Label = mcPart->GetDaughterLabel(0);
					AliAODMCParticle *mcDaughter0 = (AliAODMCParticle *)stack->UncheckedAt(daughter0Label);
					Int_t daughter0Type = -1;
					if(daughter0Label >= 0)
					{daughter0Type = mcDaughter0->GetPdgCode();}
					
					Int_t daughter1Label = mcPart->GetDaughterLabel(1);
					AliAODMCParticle *mcDaughter1 = (AliAODMCParticle *)stack->UncheckedAt(daughter1Label);
					Int_t daughter1Type = -1;
					if(daughter1Label >= 1)
					{daughter1Type = mcDaughter1->GetPdgCode();}		
					
					if( isprimaryMC && kshortMC && TMath::Abs(daughter0Type) == 211 && TMath::Abs(daughter1Type) == 211)
					{
						Double_t dz=mcDaughter0->Zv() - mcZv, dy= mcDaughter0->Yv() - mcYv, dx= mcDaughter0->Xv() - mcXv;
						Double_t mcDecayLength = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
						Double_t mccTau = mcDecayLength*(mcPart->M())/(mcPart->P());
						Double_t mcCosPA = (dx*mcPart->Px()+dy*mcPart->Py()+dz*mcPart->Pz())/(mcDecayLength*mcPart->P());
						if(centPercentile >= 0.0001 && centPercentile <= 90.0)
						{
							fHistMcGenCosPaMK0->Fill(mcCosPA,mcPart->M()); 	
							fHistMcGencTauMK0->Fill(mccTau,mcPart->M()); 	
							fHistMcGenDcaMK0->Fill(1,mcPart->M()); 
							fHistMcGenNSigmaMK0->Fill(1.0,mcPart->M()); 	
							fHistMcGenEtaMK0->Fill(mcPart->Eta(),mcPart->M()); 
							fHistMcGenRapMK0->Fill(mcPart->Y(),mcPart->M()); 
							fHistMcGenArmPodK0->Fill(1.0,1.0);
						}
						
					}	
					
					
					if( isprimaryMC && lambdaMC && ((daughter0Type == -211 && daughter1Type == 2212) || (daughter1Type == -211 && daughter0Type == 2212)) )
					{
						Double_t dz=mcDaughter0->Zv() - mcZv, dy= mcDaughter0->Yv() - mcYv, dx= mcDaughter0->Xv() - mcXv;
						Double_t mcDecayLength = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
						Double_t mccTau = mcDecayLength*(mcPart->M())/(mcPart->P());
						Double_t mcCosPA = (dx*mcPart->Px()+dy*mcPart->Py()+dz*mcPart->Pz())/(mcDecayLength*mcPart->P());
						if(centPercentile >= 0.0001 && centPercentile <= 90.0)
						{
							fHistMcGenCosPaMLa->Fill(mcCosPA,mcPart->M()); 	
							fHistMcGencTauMLa->Fill(mccTau,mcPart->M()); 	
							fHistMcGenDcaMLa->Fill(1,mcPart->M()); 
							fHistMcGenNSigmaMLa->Fill(1.0,mcPart->M()); 	
							fHistMcGenEtaMLa->Fill(mcPart->Eta(),mcPart->M()); 
							fHistMcGenRapMLa->Fill(mcPart->Y(),mcPart->M()); 
							fHistMcGenArmPodLa->Fill(1.0,1.0);
						}
						
					}
					
					
					if( isprimaryMC && antilambdaMC && ((daughter0Type == 211 && daughter1Type == -2212) || (daughter1Type == 211 && daughter0Type == -2212)) )
					{
						Double_t dz=mcDaughter0->Zv() - mcZv, dy= mcDaughter0->Yv() - mcYv, dx= mcDaughter0->Xv() - mcXv;
						Double_t mcDecayLength = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
						Double_t mccTau = mcDecayLength*(mcPart->M())/(mcPart->P());
						Double_t mcCosPA = (dx*mcPart->Px()+dy*mcPart->Py()+dz*mcPart->Pz())/(mcDecayLength*mcPart->P());
						if(centPercentile >= 0.0001 && centPercentile <= 90.0)
						{
							fHistMcGenCosPaMLb->Fill(mcCosPA,mcPart->M()); 	
							fHistMcGencTauMLb->Fill(mccTau,mcPart->M()); 	
							fHistMcGenDcaMLb->Fill(1,mcPart->M()); 
							fHistMcGenNSigmaMLb->Fill(1.0,mcPart->M()); 	
							fHistMcGenEtaMLb->Fill(mcPart->Eta(),mcPart->M()); 
							fHistMcGenRapMLb->Fill(mcPart->Y(),mcPart->M()); 
							fHistMcGenArmPodLb->Fill(1.0,1.0);
						}
						
					}
				}
				
			}
		}
	} 
	
	/********End of MC************************************/
	
	// Find vertex, check if its good
	const AliAODVertex *vtx=aod->GetPrimaryVertex();
	
	if (!vtx) 
	{
		fHistLog->Fill(97);
		return;
	}
	
	
	if (vtx->GetNContributors()<3) 
	{
		fHistLog->Fill(vtx->GetNContributors()+20);
		fHistLog->Fill(97);
		return;
	}
	
	Double_t zv=vtx->GetZ(), xv=vtx->GetX(), yv=vtx->GetY();
	
	if (TMath::Abs(zv) > 10.)
	{
		fHistLog->Fill(96);
		return;
	}  
	
	fHistZVertex->Fill(zv);
	
	if(centPercentile >= 0.0001 && centPercentile <= 5.0)
	{fHistZVertexCent0005->Fill(zv);}
	if(centPercentile > 5.0 && centPercentile <= 10.0)
	{fHistZVertexCent0510->Fill(zv);}
	if(centPercentile > 10.0 && centPercentile <= 20.0)
	{fHistZVertexCent1020->Fill(zv);}
	if(centPercentile > 20.0 && centPercentile <= 40.0)
	{fHistZVertexCent2040->Fill(zv);}
	if(centPercentile > 40.0 && centPercentile <= 60.0)
	{fHistZVertexCent4060->Fill(zv);}
	if(centPercentile > 60.0 && centPercentile <= 90.0)
	{fHistZVertexCent6090->Fill(zv);}
	if(centPercentile >= 0.0001 && centPercentile <= 90.0)
	{fHistZVertexCent0090->Fill(zv);}
	
    /********V0 loop for reconstructed event************************************/
	
    Int_t nv0s = aod->GetNumberOfV0s();
	fHistNV0->Fill(nv0s);
    
	for(Int_t i = 0; i < nv0s; i++) 
	{
		fHistLog->Fill(7);
        AliAODv0 *v0=aod->GetV0(i); // pointer to reconstructed v0          
        if(!v0) 
		{ 
			//Printf("No V0 ");
			fHistLog->Fill(94);
            continue; 
		}
		
		/* V0s not consistent with K0, Lambda (La) or Antilambda (Lb) are rejected */
		
		Bool_t lambdaCandidate = kTRUE; 
		Bool_t antilambdaCandidate = kTRUE; 
		Bool_t k0Candidate = kTRUE; 
		
		if (v0->MassLambda() < 1.08 || v0->MassLambda() > 1.2)
		{lambdaCandidate = kFALSE;}
        if (v0->MassAntiLambda() < 1.08 || v0->MassAntiLambda() > 1.2)
		{antilambdaCandidate = kFALSE;}
        if (v0->MassK0Short() < 0.414 || v0->MassK0Short() > 0.582)
		{k0Candidate = kFALSE;}
		
		if(lambdaCandidate == kFALSE && antilambdaCandidate == kFALSE && k0Candidate == kFALSE)
		{continue;}
		
		Double_t cosPA=v0->CosPointingAngle(aod->GetPrimaryVertex());
		Double_t xyz[3]; 
		v0->GetSecondaryVtx(xyz);
		Double_t decayL = TMath::Sqrt((xyz[0]-xv)*(xyz[0]-xv)+(xyz[1]-yv)*(xyz[1]-yv)+(xyz[2]-zv)*(xyz[2]-zv));
		Double_t dca=v0->DcaV0Daughters();
		Double_t eta=v0->PseudoRapV0();
		
		if(!AcceptV0_general(v0,aod,fcutCosPa,fcutNImpact,fcutDCA,fcutEta,fcutMinNClustersTPC, fcutRatio))
		{ 
			fHistLog->Fill(86);
			continue; 
		}
		
		int nnum = 1, pnum = 0;
		const AliAODTrack *ntracktest=(AliAODTrack *)v0->GetDaughterLabel(nnum);
		if(ntracktest->Charge() > 0){nnum = 0; pnum = 1;}	
		
		const AliAODTrack *ntrack1=(AliAODTrack *)v0->GetDaughterLabel(nnum);
		const AliAODTrack *ptrack1=(AliAODTrack *)v0->GetDaughterLabel(pnum);
		
		if(ntrack1->Charge()>0)
		{fHistLog->Fill(55);}
		if(ntrack1->Charge()==0)
		{ fHistLog->Fill(50);}
		if(ntrack1->Charge()<0)
		{ fHistLog->Fill(45);}
		
		const AliAODPid *pid_p1=ptrack1->GetDetPid();
		const AliAODPid *pid_n1=ntrack1->GetDetPid();
		
		if(!AcceptV0_particle(v0,1,fcutcTauMin, fcutRapidity, decayL, fcutNcTauMax))
		{ lambdaCandidate = kFALSE; }
		if(!AcceptV0_particle(v0,2,fcutcTauMin, fcutRapidity, decayL, fcutNcTauMax))
		{ antilambdaCandidate = kFALSE; }
		if(!AcceptV0_particle(v0,0,fcutcTauMin, fcutRapidity, decayL, fcutNcTauMax))
		{ k0Candidate = kFALSE; }
		
		if(TMath::Sqrt(v0->Pt2V0())<2)
		{
			if(!AcceptV0_lowpt(v0,fPIDResponse,1,fcutBetheBloch,fIsMonteCarlo))
			{ lambdaCandidate = kFALSE; }
			if(!AcceptV0_lowpt(v0,fPIDResponse,2,fcutBetheBloch,fIsMonteCarlo))
			{ antilambdaCandidate = kFALSE; }
			if(!AcceptV0_lowpt(v0,fPIDResponse,0,fcutBetheBloch,fIsMonteCarlo))
			{ k0Candidate = kFALSE; }
		}
		
		if(lambdaCandidate == kFALSE && antilambdaCandidate == kFALSE && k0Candidate == kFALSE)
		{continue;}
		
		fHistLog->Fill(7);
        fHistPt->Fill(TMath::Sqrt(v0->Pt2V0()));
        fHistEta->Fill(v0->PseudoRapV0());
		
		double ArmenterosAlpha = v0->Alpha();
		double ArmenterosPt	   = v0->QtProng();
		if( ArmenterosPt <= TMath::Abs(fcutArmenteros*ArmenterosAlpha) && fcutArmenteros !=-999 )
		{k0Candidate = false;}
		
		/* MC Associated selection is performed, but candidates failing are tagged not rejected, so MC and data can be compared */
		
		bool feeddown = false;
		bool mcPrimary2 = true;
		
		double mcAsCosPa(0), mcAsMass(0), mcAscTau(0), mcAsEta(0), mcAsRap(0);
		
		fHistLog->Fill(31);
		
		Bool_t mcAslambdaCandidate = lambdaCandidate; 
		Bool_t mcAsantilambdaCandidate = antilambdaCandidate; 
		Bool_t mcAsk0Candidate = k0Candidate; 
		if(fIsMonteCarlo)
		{
			bool passedTests = false;
			TList *list = aod->GetList();
			TClonesArray *stack = 0x0;
			stack = (TClonesArray*)list->FindObject(AliAODMCParticle::StdBranchName());
			if (!stack)
			{
				Printf("ERROR: stack not available");
				fHistLog->Fill(84);	
			}
			else
				
			{
				Int_t negAssLabel = TMath::Abs(ntrack1->GetLabel());
				Int_t posAssLabel = TMath::Abs(ptrack1->GetLabel());
				fHistLog->Fill(35);
				if(negAssLabel>=0 && negAssLabel < stack->GetEntriesFast() && posAssLabel>=0 && posAssLabel < stack->GetEntriesFast() )
				{
					fHistLog->Fill(36);
					AliAODMCParticle *mcNegPart =(AliAODMCParticle*)stack->UncheckedAt(negAssLabel);
					Int_t v0Label = mcNegPart->GetMother();
					AliAODMCParticle *mcPosPart =(AliAODMCParticle*)stack->UncheckedAt(posAssLabel);
					Int_t v0PosLabel = mcPosPart->GetMother();
					if(v0Label >= 0 && v0Label < stack->GetEntriesFast() && v0Label == v0PosLabel)
					{
						fHistLog->Fill(37);
						AliAODMCParticle *mcv0 = (AliAODMCParticle *)stack->UncheckedAt(v0Label);
						passedTests = true;
						
						if ((v0Label >= isInjected && isInjected >= 0 && isMCtype == 1) || (v0Label < isInjected && isInjected >= 0 && isMCtype == -1)) 
						{
							mcAslambdaCandidate = false;
							mcAsk0Candidate = false;
							mcAsantilambdaCandidate = false;
						}
						
						if(mcv0->GetPdgCode() != kLambda0)
						{mcAslambdaCandidate = false;}
						if(mcv0->GetPdgCode() != kK0Short)
						{mcAsk0Candidate = false;}
						if(mcv0->GetPdgCode() != kLambda0Bar)
						{mcAsantilambdaCandidate = false;}
						
						Double_t mcAsdz=mcNegPart->Zv() - mcZv, mcAsdy= mcNegPart->Yv() - mcYv, mcAsdx= mcNegPart->Xv() - mcXv;
						Double_t mcAsDecayLength = TMath::Sqrt(mcAsdx*mcAsdx + mcAsdy*mcAsdy + mcAsdz*mcAsdz);
						mcAscTau = mcAsDecayLength*(mcv0->M())/(mcv0->P());
						mcAsCosPa = (mcAsdx*mcv0->Px()+mcAsdy*mcv0->Py()+mcAsdz*mcv0->Pz())/(mcAsDecayLength*mcv0->P());
						mcAsMass = mcv0->M();
						mcAsEta = mcv0->Eta();
						mcAsRap = mcv0->Y();
						
						
						Int_t motherLabel = mcv0->GetMother();
						Int_t motherType = -1;
						bool sigma0MC2 = false;
						
						if(motherLabel >= 0 && v0Label < stack->GetEntriesFast() )
						{
							AliAODMCParticle *mcMother = (AliAODMCParticle *)stack->UncheckedAt(motherLabel);
							motherType = mcMother->GetPdgCode();
							
							// this block of code is used to include primary Sigma0 decays as primary lambda/antilambda
							
							if ((motherLabel >= isInjected && isInjected >= 0 && isMCtype == 1) || (motherLabel < isInjected && isInjected >= 0 && isMCtype == -1)) 
							{
								mcAslambdaCandidate = false;
								mcAsk0Candidate = false;
								mcAsantilambdaCandidate = false;
							}
							if(motherType == 3212 || motherType == -3212)
							{
								if(mcMother->IsPhysicalPrimary() && (lambdaCandidate || antilambdaCandidate))
									//if(checkPrimaryStatus(mcMother, xv, yv, zv))
								{sigma0MC2 = true;}
							}
							if(motherType == 3322 || motherType == -3322 || motherType == 3312 || motherType == -3312 )
							{
								if(mcMother->IsPhysicalPrimary() && (lambdaCandidate || antilambdaCandidate))
									//if(checkPrimaryStatus(mcMother, xv, yv, zv))
								{feeddown = true;}
							}
						}		
						
						
						if(!sigma0MC2 && !mcv0->IsPhysicalPrimary() && !feeddown)
							//if(!sigma0MC2 && !checkPrimaryStatus(mcv0, mcXv, mcYv, mcZv) && !feeddown)
						{
							fHistLog->Fill(38);
							mcAslambdaCandidate = false;
							mcAsk0Candidate = false;
							mcAsantilambdaCandidate = false;
						}
						
						if(!(sigma0MC2 || mcv0->IsPhysicalPrimary()))
						{mcPrimary2 = false;}
					}
				}
			}
			
			if(passedTests == false)
			{
				fHistLog->Fill(39);
				mcAslambdaCandidate = false;
				mcAsk0Candidate = false;
				mcAsantilambdaCandidate = false;
			}
			
		}
		

		fHistLog->Fill(32);
		
		/* We now fill histograms, starting with all v0 candidates, and then MC particles only */
		
		if(lambdaCandidate)
		{
			
			fHistMLaPt->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());
			if(centPercentile >= 0.0001 && centPercentile <= 5.0)
			{fHistMLaPtCent0005->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
			if(centPercentile > 5.0 && centPercentile <= 10.0)
			{fHistMLaPtCent0510->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
			if(centPercentile > 10.0 && centPercentile <= 20.0)
			{fHistMLaPtCent1020->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
			if(centPercentile > 20.0 && centPercentile <= 40.0)
			{fHistMLaPtCent2040->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
			if(centPercentile > 40.0 && centPercentile <= 60.0)
			{fHistMLaPtCent4060->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
			if(centPercentile > 60.0 && centPercentile <= 90.0)
			{fHistMLaPtCent6090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
			
			
			if(centPercentile >= 0.0001 && centPercentile <= 90.0)
			{
				
				fHistMLaPtCent0090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());
				fHistBBLaPos->Fill(pid_p1->GetTPCmomentum(),pid_p1->GetTPCsignal());
				fHistBBLaNeg->Fill(pid_n1->GetTPCmomentum(),pid_n1->GetTPCsignal());
				
				
				
				fHistCosPaMLa->Fill(cosPA,v0->MassLambda()); 
				fHistcTauMLa->Fill(decayL*(TMath::Sqrt(v0->Pt2V0()))/(v0->P()),v0->MassLambda()); 
				fHistDcaMLa->Fill(dca,v0->MassLambda()); 	
				fHistNSigmaMLa->Fill(1.0,v0->MassLambda()); 	
				fHistEtaMLa->Fill(eta,v0->MassLambda()); 
				fHistRapMLa->Fill(v0->RapLambda(),v0->MassLambda());
				fHistArmPodLa->Fill(ArmenterosAlpha,ArmenterosPt);
				
				
				
				
				
			}
			
		}
		if(antilambdaCandidate)
		{
			
			fHistMLbPt->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());  
			
			if(centPercentile >= 0.0001 && centPercentile <= 5.0)
			{fHistMLbPtCent0005->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
			if(centPercentile > 5.0 && centPercentile <= 10.0)
			{fHistMLbPtCent0510->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
			if(centPercentile > 10.0 && centPercentile <= 20.0)
			{fHistMLbPtCent1020->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
			if(centPercentile > 20.0 && centPercentile <= 40.0)
			{fHistMLbPtCent2040->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
			if(centPercentile > 40.0 && centPercentile <= 60.0)
			{fHistMLbPtCent4060->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
			if(centPercentile > 60.0 && centPercentile <= 90.0)
			{fHistMLbPtCent6090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
			
			if(centPercentile >= 0.0001 && centPercentile <= 90.0)
			{
				fHistMLbPtCent0090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());
				fHistBBLbPos->Fill(pid_p1->GetTPCmomentum(),pid_p1->GetTPCsignal());
				fHistBBLbNeg->Fill(pid_n1->GetTPCmomentum(),pid_n1->GetTPCsignal());
				
				
				
				fHistCosPaMLb->Fill(cosPA,v0->MassAntiLambda()); 
				fHistcTauMLb->Fill(decayL*(v0->MassAntiLambda())/(v0->P()),v0->MassAntiLambda()); 
				fHistDcaMLb->Fill(dca,v0->MassAntiLambda()); 	
				fHistNSigmaMLb->Fill(1.0,v0->MassAntiLambda()); 	
				fHistEtaMLb->Fill(eta,v0->MassAntiLambda()); 	
				fHistRapMLb->Fill(v0->RapLambda(),v0->MassAntiLambda()); 
				fHistArmPodLb->Fill(ArmenterosAlpha,ArmenterosPt);
				
				
				
			}
		}
		if(k0Candidate){fHistLog->Fill(33);}
		if(k0Candidate)
		{
			fHistLog->Fill(34);
			fHistMK0Pt->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());  
			
			if(centPercentile >= 0.0001 && centPercentile <= 5.0)
			{fHistMK0PtCent0005->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
			if(centPercentile > 5.0 && centPercentile <= 10.0)
			{fHistMK0PtCent0510->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
			if(centPercentile > 10.0 && centPercentile <= 20.0)
			{fHistMK0PtCent1020->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
			if(centPercentile > 20.0 && centPercentile <= 40.0)
			{fHistMK0PtCent2040->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
			if(centPercentile > 40.0 && centPercentile <= 60.0)
			{fHistMK0PtCent4060->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
			if(centPercentile > 60.0 && centPercentile <= 90.0)
			{fHistMK0PtCent6090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
			
			if(centPercentile >= 0.0001 && centPercentile <= 90.0)
			{
				fHistMK0PtCent0090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());
				fHistBBK0Pos->Fill(pid_p1->GetTPCmomentum(),pid_p1->GetTPCsignal());
				fHistBBK0Neg->Fill(pid_n1->GetTPCmomentum(),pid_n1->GetTPCsignal());
				
				
				
				fHistCosPaMK0->Fill(cosPA,v0->MassK0Short()); 
				fHistcTauMK0->Fill(decayL*(v0->MassK0Short())/(v0->P()),v0->MassK0Short()); 	
				fHistDcaMK0->Fill(dca,v0->MassK0Short()); 	
				fHistNSigmaMK0->Fill(1.0,v0->MassK0Short()); 	
				fHistEtaMK0->Fill(eta,v0->MassK0Short()); 
				fHistRapMK0->Fill(v0->RapK0Short(),v0->MassK0Short()); 
				fHistArmPodK0->Fill(ArmenterosAlpha,ArmenterosPt);
				
				
			}
			
		}
		
		/* below here, non MC-As candidates are rejected, and MC-As histograms are filled. */
		
		if(!mcAslambdaCandidate && !mcAsantilambdaCandidate && !mcAsk0Candidate)
		{continue;}
		
		if(fIsMonteCarlo)
		{
			if(mcAslambdaCandidate)
			{
				
				if(!feeddown && mcPrimary2)
				{
					if(centPercentile >= 0.0001 && centPercentile <= 5.0)
					{fHistMcAsMLaPtCent0005->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
					if(centPercentile > 5.0 && centPercentile <= 10.0)
					{fHistMcAsMLaPtCent0510->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
					if(centPercentile > 10.0 && centPercentile <= 20.0)
					{fHistMcAsMLaPtCent1020->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
					if(centPercentile > 20.0 && centPercentile <= 40.0)
					{fHistMcAsMLaPtCent2040->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
					if(centPercentile > 40.0 && centPercentile <= 60.0)
					{fHistMcAsMLaPtCent4060->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
					if(centPercentile > 60.0 && centPercentile <= 90.0)
					{fHistMcAsMLaPtCent6090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
				}
				
				if(centPercentile >= 0.0001 && centPercentile <= 90.0)
				{
					if(feeddown)
					{fHistMcFMLaPt->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
					if(!feeddown && mcPrimary2)
					{
						fHistMcAsMLaPtCent0090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());
						
						fHistMcAsReconCosPaMLa->Fill(cosPA,v0->MassLambda()); 
						fHistMcAsReconcTauMLa->Fill(decayL*(TMath::Sqrt(v0->Pt2V0()))/(v0->P()),v0->MassLambda()); 
						fHistMcAsReconDcaMLa->Fill(dca,v0->MassLambda()); 	
						fHistMcAsReconNSigmaMLa->Fill(1.0,v0->MassLambda()); 	
						fHistMcAsReconEtaMLa->Fill(eta,v0->MassLambda()); 
						fHistMcAsReconRapMLa->Fill(v0->RapLambda(),v0->MassLambda());
						fHistMcAsReconArmPodLa->Fill(ArmenterosAlpha,ArmenterosPt);
						
						
						
							fHistMcAsTruthCosPaMLa->Fill(mcAsCosPa,mcAsMass); 
							fHistMcAsTruthcTauMLa->Fill(mcAscTau,mcAsMass); 
							fHistMcAsTruthDcaMLa->Fill(1.0,mcAsMass); 	
							fHistMcAsTruthNSigmaMLa->Fill(1.0,mcAsMass); 
							fHistMcAsTruthEtaMLa->Fill(mcAsEta,mcAsMass); 	
							fHistMcAsTruthRapMLa->Fill(mcAsRap,mcAsMass); 
							fHistMcAsTruthArmPodLa->Fill(1.0,1.0);
						
					}
					
				}
				
			}
			if(mcAsantilambdaCandidate)
			{
				if(!feeddown && mcPrimary2)
				{
					if(centPercentile >= 0.0001 && centPercentile <= 5.0)
					{fHistMcAsMLbPtCent0005->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
					if(centPercentile > 5.0 && centPercentile <= 10.0)
					{fHistMcAsMLbPtCent0510->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
					if(centPercentile > 10.0 && centPercentile <= 20.0)
					{fHistMcAsMLbPtCent1020->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
					if(centPercentile > 20.0 && centPercentile <= 40.0)
					{fHistMcAsMLbPtCent2040->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
					if(centPercentile > 40.0 && centPercentile <= 60.0)
					{fHistMcAsMLbPtCent4060->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
					if(centPercentile > 60.0 && centPercentile <= 90.0)
					{fHistMcAsMLbPtCent6090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());}
					
					if(centPercentile >= 0.0001 && centPercentile <= 90.0)
					{
						fHistMcAsMLbPtCent0090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassAntiLambda());
						
						fHistMcAsReconCosPaMLb->Fill(cosPA,v0->MassAntiLambda()); 
						fHistMcAsReconcTauMLb->Fill(decayL*(v0->MassAntiLambda())/(v0->P()),v0->MassAntiLambda()); 
						fHistMcAsReconDcaMLb->Fill(dca,v0->MassAntiLambda()); 	
						fHistMcAsReconNSigmaMLb->Fill(1.0,v0->MassAntiLambda()); 	
						fHistMcAsReconEtaMLb->Fill(eta,v0->MassAntiLambda()); 	
						fHistMcAsReconRapMLb->Fill(v0->RapLambda(),v0->MassAntiLambda()); 
						fHistMcAsReconArmPodLb->Fill(ArmenterosAlpha,ArmenterosPt);
						
						
							fHistMcAsTruthCosPaMLb->Fill(mcAsCosPa,mcAsMass); 
							fHistMcAsTruthcTauMLb->Fill(mcAscTau,mcAsMass); 
							fHistMcAsTruthDcaMLb->Fill(1.0,mcAsMass); 	
							fHistMcAsTruthNSigmaMLb->Fill(1.0,mcAsMass); 	
							fHistMcAsTruthEtaMLb->Fill(mcAsEta,mcAsMass); 	
							fHistMcAsTruthRapMLb->Fill(mcAsRap,mcAsMass);
							fHistMcAsTruthArmPodLb->Fill(1.0,1.0); 
						
					}
				}
			}
			if(mcAsk0Candidate && mcPrimary2)
			{
				if(centPercentile >= 0.0001 && centPercentile <= 5.0)
				{fHistMcAsMK0PtCent0005->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
				if(centPercentile > 5.0 && centPercentile <= 10.0)
				{fHistMcAsMK0PtCent0510->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
				if(centPercentile > 10.0 && centPercentile <= 20.0)
				{fHistMcAsMK0PtCent1020->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
				if(centPercentile > 20.0 && centPercentile <= 40.0)
				{fHistMcAsMK0PtCent2040->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
				if(centPercentile > 40.0 && centPercentile <= 60.0)
				{fHistMcAsMK0PtCent4060->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
				if(centPercentile > 60.0 && centPercentile <= 90.0)
				{fHistMcAsMK0PtCent6090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());}
				if(centPercentile >= 0.0001 && centPercentile <= 90.0)
				{
					fHistMcAsMK0PtCent0090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassK0Short());
					
					fHistMcAsReconCosPaMK0->Fill(cosPA,v0->MassK0Short()); 
					fHistMcAsReconcTauMK0->Fill(decayL*(v0->MassK0Short())/(v0->P()),v0->MassK0Short()); 	
					fHistMcAsReconDcaMK0->Fill(dca,v0->MassK0Short()); 	
					fHistMcAsReconNSigmaMK0->Fill(1.0,v0->MassK0Short()); 	
					fHistMcAsReconEtaMK0->Fill(eta,v0->MassK0Short()); 
					fHistMcAsReconRapMK0->Fill(v0->RapK0Short(),v0->MassK0Short()); 
					fHistMcAsReconArmPodK0->Fill(ArmenterosAlpha,ArmenterosPt);
					
					
						fHistMcAsTruthCosPaMK0->Fill(mcAsCosPa,mcAsMass); 
						fHistMcAsTruthcTauMK0->Fill(mcAscTau,mcAsMass); 
						fHistMcAsTruthDcaMK0->Fill(1.0,mcAsMass); 
						fHistMcAsTruthNSigmaMK0->Fill(1.0,mcAsMass); 
						fHistMcAsTruthEtaMK0->Fill(mcAsEta,mcAsMass); 
						fHistMcAsTruthRapMK0->Fill(mcAsRap,mcAsMass); 
						fHistMcAsTruthArmPodK0->Fill(1.0,1.0);
					
				}
				
			}
		}
		
	}
	
	/********End of V0 loop for reconstructed event*****************************/
	
	
	
	
	fHistLog->Fill(9);
	// NEW HISTO should be filled before this point, as PostData puts the
	// information for this iteration of the UserExec in the container
	PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskLukeAOD::Terminate(Option_t *) 
{
	// Draw result to screen, or perform fitting, normalizations
	// Called once at the end of the query
	
	fOutput = dynamic_cast<TList*> (GetOutputData(1));
	if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
	
	fHistPt = dynamic_cast<TH1F*> (fOutput->FindObject("fHistPt"));
	if (!fHistPt) { Printf("ERROR: could not retrieve fHistPt"); return;}
	fHistEta = dynamic_cast<TH1F*> (fOutput->FindObject("fHistEta"));
	if (!fHistEta) { Printf("ERROR: could not retrieve fHistEta"); return;}
	
	// Get the physics selection histograms with the selection statistics
	//AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	//AliESDInputHandler *inputH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
	//TH2F *histStat = (TH2F*)inputH->GetStatistics();
	
	
	// NEW HISTO should be retrieved from the TList container in the above way,
	// so it is available to draw on a canvas such as below
	
	TCanvas *c = new TCanvas("AliAnalysisTaskLukeAOD","P_{T} & #eta",10,10,1020,510);
	c->Divide(2,1);
	c->cd(1)->SetLogy();
	fHistPt->DrawCopy("E");
	c->cd(2);
	fHistEta->DrawCopy("E");
}
