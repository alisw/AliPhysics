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
:AliAnalysisTaskSE(),
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
fHistArmPodK0(0),
fHistArmPodLa(0),
fHistArmPodLb(0),

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
fHistZVertexCent0005(0),
fHistMCZVertexCent0005(0),

fHistMK0PtCent0510(0), 
fHistMLaPtCent0510(0), 
fHistMLbPtCent0510(0),
fHistMcPMK0PtCent0510(0),
fHistMcPMLaPtCent0510(0),
fHistMcPMLbPtCent0510(0),
fHistZVertexCent0510(0),
fHistMCZVertexCent0510(0),

fHistMK0PtCent1020(0), 
fHistMLaPtCent1020(0), 
fHistMLbPtCent1020(0),
fHistMcPMK0PtCent1020(0),
fHistMcPMLaPtCent1020(0),
fHistMcPMLbPtCent1020(0),
fHistZVertexCent1020(0),
fHistMCZVertexCent1020(0),

fHistMK0PtCent2040(0), 
fHistMLaPtCent2040(0), 
fHistMLbPtCent2040(0),
fHistMcPMK0PtCent2040(0),
fHistMcPMLaPtCent2040(0),
fHistMcPMLbPtCent2040(0),
fHistZVertexCent2040(0),
fHistMCZVertexCent2040(0),

fHistMK0PtCent4060(0), 
fHistMLaPtCent4060(0), 
fHistMLbPtCent4060(0),
fHistMcPMK0PtCent4060(0),
fHistMcPMLaPtCent4060(0),
fHistMcPMLbPtCent4060(0),
fHistZVertexCent4060(0),
fHistMCZVertexCent4060(0),

fHistMK0PtCent6090(0), 
fHistMLaPtCent6090(0), 
fHistMLbPtCent6090(0),
fHistMcPMK0PtCent6090(0),
fHistMcPMLaPtCent6090(0),
fHistMcPMLbPtCent6090(0),
fHistZVertexCent6090(0),
fHistMCZVertexCent6090(0),

fHistMK0PtCent0090(0), 
fHistMLaPtCent0090(0), 
fHistMLbPtCent0090(0),
fHistMcPMK0PtCent0090(0),
fHistMcPMLaPtCent0090(0),
fHistMcPMLbPtCent0090(0),
fHistZVertexCent0090(0),
fHistMCZVertexCent0090(0),

fHistCosPaLaPt(0),
fHistCosPaLbPt(0),
fHistCosPaK0Pt(0),	
fHistMcCosPaAllLaPt(0),
fHistMcCosPaAllLbPt(0),
fHistMcCosPaAllK0Pt(0),	
fHistMcCosPaFoundLaPt(0),
fHistMcCosPaFoundLbPt(0),
fHistMcCosPaAFoundK0Pt(0),

fHistcTauLaPt(0),	
fHistcTauLbPt(0),
fHistcTauK0Pt(0),		
fHistMccTauAllLaPt(0),
fHistMccTauAllLbPt(0),
fHistMccTauAllK0Pt(0),	
fHistMccTauFoundLaPt(0),
fHistMccTauFoundLbPt(0),	
fHistMccTauAFoundK0Pt(0),

fHistDcaLaPt(0),	
fHistDcaLbPt(0),	
fHistDcaK0Pt(0),	
fHistMcDcaAllLaPt(0),
fHistMcDcaAllLbPt(0),
fHistMcDcaAllK0Pt(0),	
fHistMcDcaFoundLaPt(0),	
fHistMcDcaFoundLbPt(0),	
fHistMcDcaAFoundK0Pt(0),

fHistNSigmaLaPt(0),	
fHistNSigmaLbPt(0),		
fHistNSigmaK0Pt(0),		
fHistMcNSigmaAllLaPt(0),	
fHistMcNSigmaAllLbPt(0),	
fHistMcNSigmaAllK0Pt(0),	
fHistMcNSigmaFoundLaPt(0),	
fHistMcNSigmaFoundLbPt(0),	
fHistMcNSigmaAFoundK0Pt(0),

fHistEtaLaPt(0),	
fHistEtaLbPt(0),	
fHistEtaK0Pt(0),		
fHistMcEtaAllLaPt(0),
fHistMcEtaAllLbPt(0),	
fHistMcEtaAllK0Pt(0),
fHistMcEtaFoundLaPt(0),	
fHistMcEtaFoundLbPt(0),	
fHistMcEtaAFoundK0Pt(0),

fHistRapLaPt(0),
fHistRapLbPt(0),		
fHistRapK0Pt(0),	
fHistMcRapAllLaPt(0),	
fHistMcRapAllLbPt(0),	
fHistMcRapAllK0Pt(0),	
fHistMcRapFoundLaPt(0),
fHistMcRapFoundLbPt(0),
fHistMcRapAFoundK0Pt(0)

// The last in the above list should not have a comma after it
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskLukeAOD::AliAnalysisTaskLukeAOD(const char *name) // All data members should be initialised here
:AliAnalysisTaskSE(name),
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
fHistArmPodK0(0),
fHistArmPodLa(0),
fHistArmPodLb(0),

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
fHistZVertexCent0005(0),
fHistMCZVertexCent0005(0),

fHistMK0PtCent0510(0), 
fHistMLaPtCent0510(0), 
fHistMLbPtCent0510(0),
fHistMcPMK0PtCent0510(0),
fHistMcPMLaPtCent0510(0),
fHistMcPMLbPtCent0510(0),
fHistZVertexCent0510(0),
fHistMCZVertexCent0510(0),

fHistMK0PtCent1020(0), 
fHistMLaPtCent1020(0), 
fHistMLbPtCent1020(0),
fHistMcPMK0PtCent1020(0),
fHistMcPMLaPtCent1020(0),
fHistMcPMLbPtCent1020(0),
fHistZVertexCent1020(0),
fHistMCZVertexCent1020(0),

fHistMK0PtCent2040(0), 
fHistMLaPtCent2040(0), 
fHistMLbPtCent2040(0),
fHistMcPMK0PtCent2040(0),
fHistMcPMLaPtCent2040(0),
fHistMcPMLbPtCent2040(0),
fHistZVertexCent2040(0),
fHistMCZVertexCent2040(0),

fHistMK0PtCent4060(0), 
fHistMLaPtCent4060(0), 
fHistMLbPtCent4060(0),
fHistMcPMK0PtCent4060(0),
fHistMcPMLaPtCent4060(0),
fHistMcPMLbPtCent4060(0),
fHistZVertexCent4060(0),
fHistMCZVertexCent4060(0),

fHistMK0PtCent6090(0), 
fHistMLaPtCent6090(0), 
fHistMLbPtCent6090(0),
fHistMcPMK0PtCent6090(0),
fHistMcPMLaPtCent6090(0),
fHistMcPMLbPtCent6090(0),
fHistZVertexCent6090(0),
fHistMCZVertexCent6090(0),

fHistMK0PtCent0090(0), 
fHistMLaPtCent0090(0), 
fHistMLbPtCent0090(0),
fHistMcPMK0PtCent0090(0),
fHistMcPMLaPtCent0090(0),
fHistMcPMLbPtCent0090(0),
fHistZVertexCent0090(0),
fHistMCZVertexCent0090(0),

fHistCosPaLaPt(0),
fHistCosPaLbPt(0),
fHistCosPaK0Pt(0),	
fHistMcCosPaAllLaPt(0),
fHistMcCosPaAllLbPt(0),
fHistMcCosPaAllK0Pt(0),	
fHistMcCosPaFoundLaPt(0),
fHistMcCosPaFoundLbPt(0),
fHistMcCosPaAFoundK0Pt(0),

fHistcTauLaPt(0),	
fHistcTauLbPt(0),
fHistcTauK0Pt(0),		
fHistMccTauAllLaPt(0),
fHistMccTauAllLbPt(0),
fHistMccTauAllK0Pt(0),	
fHistMccTauFoundLaPt(0),
fHistMccTauFoundLbPt(0),	
fHistMccTauAFoundK0Pt(0),

fHistDcaLaPt(0),	
fHistDcaLbPt(0),	
fHistDcaK0Pt(0),	
fHistMcDcaAllLaPt(0),
fHistMcDcaAllLbPt(0),
fHistMcDcaAllK0Pt(0),	
fHistMcDcaFoundLaPt(0),	
fHistMcDcaFoundLbPt(0),	
fHistMcDcaAFoundK0Pt(0),

fHistNSigmaLaPt(0),	
fHistNSigmaLbPt(0),		
fHistNSigmaK0Pt(0),		
fHistMcNSigmaAllLaPt(0),	
fHistMcNSigmaAllLbPt(0),	
fHistMcNSigmaAllK0Pt(0),	
fHistMcNSigmaFoundLaPt(0),	
fHistMcNSigmaFoundLbPt(0),	
fHistMcNSigmaAFoundK0Pt(0),

fHistEtaLaPt(0),	
fHistEtaLbPt(0),	
fHistEtaK0Pt(0),		
fHistMcEtaAllLaPt(0),
fHistMcEtaAllLbPt(0),	
fHistMcEtaAllK0Pt(0),
fHistMcEtaFoundLaPt(0),	
fHistMcEtaFoundLbPt(0),	
fHistMcEtaAFoundK0Pt(0),

fHistRapLaPt(0),
fHistRapLbPt(0),		
fHistRapK0Pt(0),	
fHistMcRapAllLaPt(0),	
fHistMcRapAllLbPt(0),	
fHistMcRapAllK0Pt(0),	
fHistMcRapFoundLaPt(0),
fHistMcRapFoundLbPt(0),
fHistMcRapAFoundK0Pt(0)

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
	maskIsSelected = inputHandler->IsEventSelected();
    
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
	fHistArmPodK0 = new	TH2F("fHistArmPodK0","Armenteros plot for K0 candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistArmPodLa = new	TH2F("fHistArmPodLa","Armenteros plot for La candidates; Alpha; PtArm",100,-1,1,50,0,0.5);
	fHistArmPodLb = new	TH2F("fHistArmPodLb","Armenteros plot for Lb candidates; Alpha; PtArm",100,-1,1,50,0,0.5);

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
	fHistZVertexCent0005 = new TH1F("fHistZVertexCent0005","Z coordinate of primary vertex for Centrality 0-5%",60, -15, 15);
	fHistMCZVertexCent0005 = new TH1F("fHistMCZVertexCent0005","Z coordinate of primary vertex in MC truth for Centrality 0-5%",60, -15, 15);
	
	fHistMK0PtCent0510 = new	TH2F("fHistMK0PtCent0510","K0 Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent0510 = new	TH2F("fHistMLaPtCent0510","Lambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent0510 = new	TH2F("fHistMLbPtCent0510","AntiLambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent0510 = new	TH2F("fHistMcPMK0PtCent0510","Monte Carlo primary K0 Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent0510 = new	TH2F("fHistMcPMLaPtCent0510","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent0510 = new	TH2F("fHistMcPMLbPtCent0510","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 5-10%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent0510 = new TH1F("fHistZVertexCent0510","Z coordinate of primary vertex for Centrality 5-10%",60, -15, 15);
	fHistMCZVertexCent0510 = new TH1F("fHistMCZVertexCent0510","Z coordinate of primary vertex in MC truth for Centrality 5-10%",60, -15, 15);
	
	fHistMK0PtCent1020 = new	TH2F("fHistMK0PtCent1020","K0 Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent1020 = new	TH2F("fHistMLaPtCent1020","Lambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent1020 = new	TH2F("fHistMLbPtCent1020","AntiLambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent1020 = new	TH2F("fHistMcPMK0PtCent1020","Monte Carlo primary K0 Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent1020 = new	TH2F("fHistMcPMLaPtCent1020","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent1020 = new	TH2F("fHistMcPMLbPtCent1020","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 10-20%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent1020 = new TH1F("fHistZVertexCent1020","Z coordinate of primary vertex for Centrality 10-20%",60, -15, 15);
	fHistMCZVertexCent1020 = new TH1F("fHistMCZVertexCent1020","Z coordinate of primary vertex in MC truth for Centrality 10-20%",60, -15, 15);
	
	fHistMK0PtCent2040 = new	TH2F("fHistMK0PtCent2040","K0 Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent2040 = new	TH2F("fHistMLaPtCent2040","Lambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent2040 = new	TH2F("fHistMLbPtCent2040","AntiLambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent2040 = new	TH2F("fHistMcPMK0PtCent2040","Monte Carlo primary K0 Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent2040 = new	TH2F("fHistMcPMLaPtCent2040","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent2040 = new	TH2F("fHistMcPMLbPtCent2040","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 20-40%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent2040 = new TH1F("fHistZVertexCent2040","Z coordinate of primary vertex for Centrality 20-40%",60, -15, 15);
	fHistMCZVertexCent2040 = new TH1F("fHistMCZVertexCent2040","Z coordinate of primary vertex in MC truth for Centrality 20-40%",60, -15, 15);
	
	fHistMK0PtCent4060 = new	TH2F("fHistMK0PtCent4060","K0 Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent4060 = new	TH2F("fHistMLaPtCent4060","Lambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent4060 = new	TH2F("fHistMLbPtCent4060","AntiLambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent4060 = new	TH2F("fHistMcPMK0PtCent4060","Monte Carlo primary K0 Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent4060 = new	TH2F("fHistMcPMLaPtCent4060","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent4060 = new	TH2F("fHistMcPMLbPtCent4060","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 40-60%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent4060 = new TH1F("fHistZVertexCent4060","Z coordinate of primary vertex for Centrality 40-60%",60, -15, 15);
	fHistMCZVertexCent4060 = new TH1F("fHistMCZVertexCent4060","Z coordinate of primary vertex in MC truth for Centrality 40-60%",60, -15, 15);
	
	fHistMK0PtCent6090 = new	TH2F("fHistMK0PtCent6090","K0 Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent6090 = new	TH2F("fHistMLaPtCent6090","Lambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent6090 = new	TH2F("fHistMLbPtCent6090","AntiLambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent6090 = new	TH2F("fHistMcPMK0PtCent6090","Monte Carlo primary K0 Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent6090 = new	TH2F("fHistMcPMLaPtCent6090","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent6090 = new	TH2F("fHistMcPMLbPtCent6090","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 60-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent6090 = new TH1F("fHistZVertexCent6090","Z coordinate of primary vertex for Centrality 60-90%",60, -15, 15);
	fHistMCZVertexCent6090 = new TH1F("fHistMCZVertexCent6090","Z coordinate of primary vertex in MC truth for Centrality 60-90%",60, -15, 15);
	
	fHistMK0PtCent0090 = new	TH2F("fHistMK0PtCent0090","K0 Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMLaPtCent0090 = new	TH2F("fHistMLaPtCent0090","Lambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMLbPtCent0090 = new	TH2F("fHistMLbPtCent0090","AntiLambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMK0PtCent0090 = new	TH2F("fHistMcPMK0PtCent0090","Monte Carlo primary K0 Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); K0 Mass (GeV/c^2)",200,0,10,140,0.414,0.582);
	fHistMcPMLaPtCent0090 = new	TH2F("fHistMcPMLaPtCent0090","Monte Carlo primary (& sigma0) Lambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(p#pi^{-}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistMcPMLbPtCent0090 = new	TH2F("fHistMcPMLbPtCent0090","Monte Carlo primary (& sigma0) AntiLambda Mass versus Pt for Centrality 0-90%; P_{perp} (GeV/c); M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,10,96,1.08,1.2);
	fHistZVertexCent0090 = new TH1F("fHistZVertexCent0090","Z coordinate of primary vertex for Centrality 0-90%",60, -15, 15);
	fHistMCZVertexCent0090 = new TH1F("fHistMCZVertexCent0090","Z coordinate of primary vertex in MC truth for Centrality 0-90%",60, -15, 15);
	
	
	fHistCosPaLaPt		 = new	TH2F("fHistCosPaLaPt","	Reconstructed Mass vs CosPa for Lambda Candidates; M(p#pi^{-}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistCosPaLbPt		 = new	TH2F("fHistCosPaLbPt","	Reconstructed Mass vs CosPa for AntiLambda Candidates; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistCosPaK0Pt		 = new	TH2F("fHistCosPaK0Pt","	Reconstructed Mass vs CosPa for K0Short Candidates; K0 Mass (GeV/c^2)",200,0.99,1.001,140,0.414,0.582);
	fHistMcCosPaAllLaPt	 = new	TH2F("fHistMcCosPaAllLaPt","	Reconstructed Mass vs CosPa for all MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcCosPaAllLbPt	 = new	TH2F("fHistMcCosPaAllLbPt","	Reconstructed Mass vs CosPa for all MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcCosPaAllK0Pt	 = new	TH2F("fHistMcCosPaAllK0Pt","	Reconstructed Mass vs CosPa for all MC primary K0Short; K0 Mass (GeV/c^2)",200,0.99,1.001,140,0.414,0.582);
	fHistMcCosPaFoundLaPt	 = new	TH2F("fHistMcCosPaFoundLaPt","	Reconstructed Mass vs CosPa for reconstructed MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcCosPaFoundLbPt	 = new	TH2F("fHistMcCosPaFoundLbPt","	Reconstructed Mass vs CosPa for reconstructed MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0.99,1.001,96,1.08,1.2);
	fHistMcCosPaAFoundK0Pt = new	TH2F("fHistMcCosPaAFoundK0Pt","	Reconstructed Mass vs CosPa for reconstructed MC primary K0Short; K0 Mass (GeV/c^2)",200,0.99,1.001,140,0.414,0.582);
	
	fHistcTauLaPt		 = new	TH2F("fHistcTauLaPt","	Reconstructed Mass vs cTau for Lambda Candidates; M(p#pi^{-}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistcTauLbPt		 = new	TH2F("fHistcTauLbPt","	Reconstructed Mass vs cTau for AntiLambda Candidates; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistcTauK0Pt		 = new	TH2F("fHistcTauK0Pt","	Reconstructed Mass vs cTau for K0Short Candidates; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	fHistMccTauAllLaPt	 = new	TH2F("fHistMccTauAllLaPt","	Reconstructed Mass vs cTau for all MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMccTauAllLbPt	 = new	TH2F("fHistMccTauAllLbPt","	Reconstructed Mass vs cTau for all MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMccTauAllK0Pt	 = new	TH2F("fHistMccTauAllK0Pt","	Reconstructed Mass vs cTau for all MC primary K0Short; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	fHistMccTauFoundLaPt	 = new	TH2F("fHistMccTauFoundLaPt","	Reconstructed Mass vs cTau for reconstructed MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMccTauFoundLbPt	 = new	TH2F("fHistMccTauFoundLbPt","	Reconstructed Mass vs cTau for reconstructed MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMccTauAFoundK0Pt = new	TH2F("fHistMccTauAFoundK0Pt","	Reconstructed Mass vs cTau for reconstructed MC primary K0Short; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	
	fHistDcaLaPt		 = new	TH2F("fHistDcaLaPt","	Reconstructed Mass vs Dca for Lambda Candidates; M(p#pi^{-}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistDcaLbPt		 = new	TH2F("fHistDcaLbPt","	Reconstructed Mass vs Dca for AntiLambda Candidates; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistDcaK0Pt		 = new	TH2F("fHistDcaK0Pt","	Reconstructed Mass vs Dca for K0Short Candidates; K0 Mass (GeV/c^2)",200,0,1.5,140,0.414,0.582);
	fHistMcDcaAllLaPt	 = new	TH2F("fHistMcDcaAllLaPt","	Reconstructed Mass vs Dca for all MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcDcaAllLbPt	 = new	TH2F("fHistMcDcaAllLbPt","	Reconstructed Mass vs Dca for all MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,50,96,1.08,1.2);
	fHistMcDcaAllK0Pt	 = new	TH2F("fHistMcDcaAllK0Pt","	Reconstructed Mass vs Dca for all MC primary K0Short; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	fHistMcDcaFoundLaPt	 = new	TH2F("fHistMcDcaFoundLaPt","	Reconstructed Mass vs Dca for reconstructed MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistMcDcaFoundLbPt	 = new	TH2F("fHistMcDcaFoundLbPt","	Reconstructed Mass vs Dca for reconstructed MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,0,1.5,96,1.08,1.2);
	fHistMcDcaAFoundK0Pt = new	TH2F("fHistMcDcaAFoundK0Pt","	Reconstructed Mass vs Dca for reconstructed MC primary K0Short; K0 Mass (GeV/c^2)",200,0,1.5,140,0.414,0.582);
	
	fHistNSigmaLaPt		 = new	TH2F("fHistNSigmaLaPt","	Reconstructed Mass vs NSigma for Lambda Candidates; M(p#pi^{-}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistNSigmaLbPt		 = new	TH2F("fHistNSigmaLbPt","	Reconstructed Mass vs NSigma for AntiLambda Candidates; M(#bar{p}#pi^{+}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistNSigmaK0Pt		 = new	TH2F("fHistNSigmaK0Pt","	Reconstructed Mass vs NSigma for K0Short Candidates; K0 Mass (GeV/c^2)",50,-5,5,140,0.414,0.582);
	fHistMcNSigmaAllLaPt	 = new	TH2F("fHistMcNSigmaAllLaPt","	Reconstructed Mass vs NSigma for all MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcNSigmaAllLbPt	 = new	TH2F("fHistMcNSigmaAllLbPt","	Reconstructed Mass vs NSigma for all MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcNSigmaAllK0Pt	 = new	TH2F("fHistMcNSigmaAllK0Pt","	Reconstructed Mass vs NSigma for all MC primary K0Short; K0 Mass (GeV/c^2)",50,-5,5,140,0.414,0.582);
	fHistMcNSigmaFoundLaPt	 = new	TH2F("fHistMcNSigmaFoundLaPt","	Reconstructed Mass vs NSigma for reconstructed MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcNSigmaFoundLbPt	 = new	TH2F("fHistMcNSigmaFoundLbPt","	Reconstructed Mass vs NSigma for reconstructed MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",50,-5,5,96,1.08,1.2);
	fHistMcNSigmaAFoundK0Pt = new	TH2F("fHistMcNSigmaAFoundK0Pt","	Reconstructed Mass vs NSigma for reconstructed MC primary K0Short; K0 Mass (GeV/c^2)",50,-5,5,140,0.414,0.582);
	
	fHistEtaLaPt		 = new	TH2F("fHistEtaLaPt","	Reconstructed Mass vs Eta for Lambda Candidates; M(p#pi^{-}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistEtaLbPt		 = new	TH2F("fHistEtaLbPt","	Reconstructed Mass vs Eta for AntiLambda Candidates; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistEtaK0Pt		 = new	TH2F("fHistEtaK0Pt","	Reconstructed Mass vs Eta for K0Short Candidates; K0 Mass (GeV/c^2)",200,-3,3,140,0.414,0.582);
	fHistMcEtaAllLaPt	 = new	TH2F("fHistMcEtaAllLaPt","	Reconstructed Mass vs Eta for all MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcEtaAllLbPt	 = new	TH2F("fHistMcEtaAllLbPt","	Reconstructed Mass vs Eta for all MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcEtaAllK0Pt	 = new	TH2F("fHistMcEtaAllK0Pt","	Reconstructed Mass vs Eta for all MC primary K0Short; K0 Mass (GeV/c^2)",200,-3,3,140,0.414,0.582);
	fHistMcEtaFoundLaPt	 = new	TH2F("fHistMcEtaFoundLaPt","	Reconstructed Mass vs Eta for reconstructed MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcEtaFoundLbPt	 = new	TH2F("fHistMcEtaFoundLbPt","	Reconstructed Mass vs Eta for reconstructed MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-3,3,96,1.08,1.2);
	fHistMcEtaAFoundK0Pt = new	TH2F("fHistMcEtaAFoundK0Pt","	Reconstructed Mass vs Eta for reconstructed MC primary K0Short; K0 Mass (GeV/c^2)",200,-3,3,140,0.414,0.582);
	
	fHistRapLaPt		 = new	TH2F("fHistRapLaPt","	Reconstructed Mass vs Rap for Lambda Candidates; M(p#pi^{-}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistRapLbPt		 = new	TH2F("fHistRapLbPt","	Reconstructed Mass vs Rap for AntiLambda Candidates; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistRapK0Pt		 = new	TH2F("fHistRapK0Pt","	Reconstructed Mass vs Rap for K0Short Candidates; K0 Mass (GeV/c^2)",200,-1,1,140,0.414,0.582);
	fHistMcRapAllLaPt	 = new	TH2F("fHistMcRapAllLaPt","	Reconstructed Mass vs Rap for all MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcRapAllLbPt	 = new	TH2F("fHistMcRapAllLbPt","	Reconstructed Mass vs Rap for all MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcRapAllK0Pt	 = new	TH2F("fHistMcRapAllK0Pt","	Reconstructed Mass vs Rap for all MC primary K0Short; K0 Mass (GeV/c^2)",200,-1,1,140,0.414,0.582);
	fHistMcRapFoundLaPt	 = new	TH2F("fHistMcRapFoundLaPt","	Reconstructed Mass vs Rap for reconstructed MC primary Lambda; M(p#pi^{-}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcRapFoundLbPt	 = new	TH2F("fHistMcRapFoundLbPt","	Reconstructed Mass vs Rap for reconstructed MC primary AntiLambda; M(#bar{p}#pi^{+}) (GeV/c^2)",200,-1,1,96,1.08,1.2);
	fHistMcRapAFoundK0Pt = new	TH2F("fHistMcRapAFoundK0Pt","	Reconstructed Mass vs Rap for reconstructed MC primary K0Short; K0 Mass (GeV/c^2)",200,-1,1,140,0.414,0.582);
	
	
	
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
	fOutput->Add(fHistArmPodK0);
	fOutput->Add(fHistArmPodLa);
	fOutput->Add(fHistArmPodLb);
	
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
	fOutput->Add(fHistZVertexCent0005);
	fOutput->Add(fHistMCZVertexCent0005);
	
	fOutput->Add(fHistMK0PtCent0510); 
	fOutput->Add(fHistMLaPtCent0510); 
	fOutput->Add(fHistMLbPtCent0510);
	fOutput->Add(fHistMcPMK0PtCent0510);
	fOutput->Add(fHistMcPMLaPtCent0510);
	fOutput->Add(fHistMcPMLbPtCent0510);
	fOutput->Add(fHistZVertexCent0510);
	fOutput->Add(fHistMCZVertexCent0510);
	
	fOutput->Add(fHistMK0PtCent1020); 
	fOutput->Add(fHistMLaPtCent1020); 
	fOutput->Add(fHistMLbPtCent1020);
	fOutput->Add(fHistMcPMK0PtCent1020);
	fOutput->Add(fHistMcPMLaPtCent1020);
	fOutput->Add(fHistMcPMLbPtCent1020);
	fOutput->Add(fHistZVertexCent1020);
	fOutput->Add(fHistMCZVertexCent1020);
	
	fOutput->Add(fHistMK0PtCent2040); 
	fOutput->Add(fHistMLaPtCent2040); 
	fOutput->Add(fHistMLbPtCent2040);
	fOutput->Add(fHistMcPMK0PtCent2040);
	fOutput->Add(fHistMcPMLaPtCent2040);
	fOutput->Add(fHistMcPMLbPtCent2040);
	fOutput->Add(fHistZVertexCent2040);
	fOutput->Add(fHistMCZVertexCent2040);
	
	fOutput->Add(fHistMK0PtCent4060); 
	fOutput->Add(fHistMLaPtCent4060); 
	fOutput->Add(fHistMLbPtCent4060);
	fOutput->Add(fHistMcPMK0PtCent4060);
	fOutput->Add(fHistMcPMLaPtCent4060);
	fOutput->Add(fHistMcPMLbPtCent4060);
	fOutput->Add(fHistZVertexCent4060);
	fOutput->Add(fHistMCZVertexCent4060);
	
	fOutput->Add(fHistMK0PtCent6090); 
	fOutput->Add(fHistMLaPtCent6090); 
	fOutput->Add(fHistMLbPtCent6090);
	fOutput->Add(fHistMcPMK0PtCent6090);
	fOutput->Add(fHistMcPMLaPtCent6090);
	fOutput->Add(fHistMcPMLbPtCent6090);
	fOutput->Add(fHistZVertexCent6090);
	fOutput->Add(fHistMCZVertexCent6090);
	
	fOutput->Add(fHistMK0PtCent0090); 
	fOutput->Add(fHistMLaPtCent0090); 
	fOutput->Add(fHistMLbPtCent0090);
	fOutput->Add(fHistMcPMK0PtCent0090);
	fOutput->Add(fHistMcPMLaPtCent0090);
	fOutput->Add(fHistMcPMLbPtCent0090);
	fOutput->Add(fHistZVertexCent0090);
	fOutput->Add(fHistMCZVertexCent0090);
	
	fOutput->Add(fHistCosPaLaPt);
	fOutput->Add(fHistCosPaLbPt);
	fOutput->Add(fHistCosPaK0Pt);	
	fOutput->Add(fHistMcCosPaAllLaPt);
	fOutput->Add(fHistMcCosPaAllLbPt);
	fOutput->Add(fHistMcCosPaAllK0Pt);	
	fOutput->Add(fHistMcCosPaFoundLaPt);
	fOutput->Add(fHistMcCosPaFoundLbPt);
	fOutput->Add(fHistMcCosPaAFoundK0Pt);
	
	fOutput->Add(fHistcTauLaPt);	
	fOutput->Add(fHistcTauLbPt);
	fOutput->Add(fHistcTauK0Pt);		
	fOutput->Add(fHistMccTauAllLaPt);
	fOutput->Add(fHistMccTauAllLbPt);
	fOutput->Add(fHistMccTauAllK0Pt);	
	fOutput->Add(fHistMccTauFoundLaPt);
	fOutput->Add(fHistMccTauFoundLbPt);	
	fOutput->Add(fHistMccTauAFoundK0Pt);
	
	fOutput->Add(fHistDcaLaPt);	
	fOutput->Add(fHistDcaLbPt);	
	fOutput->Add(fHistDcaK0Pt);	
	fOutput->Add(fHistMcDcaAllLaPt);
	fOutput->Add(fHistMcDcaAllLbPt);
	fOutput->Add(fHistMcDcaAllK0Pt);	
	fOutput->Add(fHistMcDcaFoundLaPt);	
	fOutput->Add(fHistMcDcaFoundLbPt);	
	fOutput->Add(fHistMcDcaAFoundK0Pt);
	
	fOutput->Add(fHistNSigmaLaPt);	
	fOutput->Add(fHistNSigmaLbPt);		
	fOutput->Add(fHistNSigmaK0Pt);		
	fOutput->Add(fHistMcNSigmaAllLaPt);	
	fOutput->Add(fHistMcNSigmaAllLbPt);	
	fOutput->Add(fHistMcNSigmaAllK0Pt);	
	fOutput->Add(fHistMcNSigmaFoundLaPt);	
	fOutput->Add(fHistMcNSigmaFoundLbPt);	
	fOutput->Add(fHistMcNSigmaAFoundK0Pt);
	
	fOutput->Add(fHistEtaLaPt);	
	fOutput->Add(fHistEtaLbPt);	
	fOutput->Add(fHistEtaK0Pt);		
	fOutput->Add(fHistMcEtaAllLaPt);
	fOutput->Add(fHistMcEtaAllLbPt);	
	fOutput->Add(fHistMcEtaAllK0Pt);
	fOutput->Add(fHistMcEtaFoundLaPt);	
	fOutput->Add(fHistMcEtaFoundLbPt);	
	fOutput->Add(fHistMcEtaAFoundK0Pt);
	
	fOutput->Add(fHistRapLaPt);
	fOutput->Add(fHistRapLbPt);		
	fOutput->Add(fHistRapK0Pt);	
	fOutput->Add(fHistMcRapAllLaPt);	
	fOutput->Add(fHistMcRapAllLbPt);	
	fOutput->Add(fHistMcRapAllK0Pt);	
	fOutput->Add(fHistMcRapFoundLaPt);
	fOutput->Add(fHistMcRapFoundLbPt);
	fOutput->Add(fHistMcRapAFoundK0Pt);
    
	// NEW HISTO added to fOutput here
    PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________

static Bool_t AcceptTrack(const AliAODTrack *t, double cutMinNClustersTPC, double cutRatio) 
{
	if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
	//if (t->GetKinkIndex(0)>0) return kFALSE;
	
	Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1); 
	if (nCrossedRowsTPC < cutMinNClustersTPC) return kFALSE;
	Int_t findable=t->GetTPCNclsF();
	if (findable <= 0) return kFALSE;
	if (nCrossedRowsTPC/findable < cutRatio) return kFALSE;
	
	return kTRUE;
}

//________________________________________________________________________
static Bool_t AcceptV0_general(const AliAODv0 *v1, const AliAODEvent *aod, double cutCosPa, double cutNImpact, double cutDCA, double cutEta, double cutMinNClustersTPC, double cutRatio) 
{
	
	if (v1->GetOnFlyStatus()) return kFALSE;
	
	int nnum = 1, pnum = 0;	
	const AliAODTrack *ntracktest=(AliAODTrack *)v1->GetDaughter(nnum);
	if(ntracktest->Charge() > 0){nnum = 0; pnum = 1;}	
	
	const AliAODTrack *ntrack1=(AliAODTrack *)v1->GetDaughter(nnum);
	if (!AcceptTrack(ntrack1, cutMinNClustersTPC, cutRatio)) return kFALSE;
	
	const AliAODTrack *ptrack1=(AliAODTrack *)v1->GetDaughter(pnum);
	if (!AcceptTrack(ptrack1, cutMinNClustersTPC, cutRatio)) return kFALSE;
	
	Float_t impact=v1->DcaNegToPrimVertex();
	if (TMath::Abs(impact)<0.1) return kFALSE;
	if (TMath::Abs(impact)<cutNImpact && cutNImpact != -999) return kFALSE;
	impact=v1->DcaPosToPrimVertex();
	if (TMath::Abs(impact)<0.1) return kFALSE;
	if (TMath::Abs(impact)<cutNImpact && cutNImpact != -999) return kFALSE;
	
	Double_t dca=v1->DcaV0Daughters();
	if (TMath::Abs(dca)>cutDCA && cutDCA != -999) return kFALSE;
	
	Double_t cpa=v1->CosPointingAngle(aod->GetPrimaryVertex());
	if (cpa<cutCosPa && cutCosPa != -999) return kFALSE;
	
	Double_t etaN = v1->PseudoRapNeg();
	Double_t etaP = v1->PseudoRapPos();
	if ((TMath::Abs(etaN)>cutEta || TMath::Abs(etaP)>cutEta) && cutEta != -999) return kFALSE;
	
	return kTRUE;
}

//________________________________________________________________________
static Bool_t AcceptV0_particle(const AliAODv0 *v1, int type,  double cutcTau, double cutRapidity, Double_t decayL) 
{
	
	Double_t cTau = 0;
	if(type == 1)
	{cTau = decayL*(v1->MassLambda())/(v1->P());}
	if(type == 2)
	{cTau = decayL*(v1->MassAntiLambda())/(v1->P());}
	if(type == 0)
	{cTau = decayL*(v1->MassK0Short())/(v1->P());}
	
	if (cTau < cutcTau && cTau != -999 ) return kFALSE;
	
	Double_t rap = 0;
	if(type == 1 || type == 2)
	{rap = v1->RapLambda();}
	if(type == 0)
	{rap = v1->RapK0Short();}
	if (TMath::Abs(rap)>cutRapidity && cutRapidity != -999) return kFALSE;
	
	return kTRUE;
}

//________________________________________________________________________
static Bool_t AcceptV0_lowpt(const AliAODv0 *v1, AliPIDResponse *PIDResponse,int type, double cutBetheBloch, Double_t decayL, bool isMonteCarlo) 
{
	
	
	Double_t cTau = 0;
	if(type == 1)
	{cTau = decayL*(v1->MassLambda())/(v1->P());}
	if(type == 2)
	{cTau = decayL*(v1->MassAntiLambda())/(v1->P());}
	if(type == 0)
	{cTau = decayL*(v1->MassK0Short())/(v1->P());}
	
	if (cTau > (3*7.89) && (type ==1 || type ==2)) return kFALSE;
	if (cTau > (3*2.68) && (type ==0)) return kFALSE;
	
	int nnum = 1, pnum = 0;
	const AliAODTrack *ntracktest=(AliAODTrack *)v1->GetDaughter(nnum);
	if(ntracktest->Charge() > 0){nnum = 0; pnum = 1;}	
	
	const AliAODTrack *ntrack1=(AliAODTrack *)v1->GetDaughter(nnum);
	const AliAODTrack *ptrack1=(AliAODTrack *)v1->GetDaughter(pnum);
	
	Double_t nsig_p = 0;
	Double_t nsig_n = 0;
	
	const AliAODPid *pid_p=ptrack1->GetDetPid();
	const AliAODPid *pid_n=ntrack1->GetDetPid();
	
	if (pid_p) 
	{
		if(type == 1)
		{
			nsig_p=PIDResponse->NumberOfSigmasTPC(ptrack1,AliPID::kProton);
			if (TMath::Abs(nsig_p) > cutBetheBloch && cutBetheBloch >0 && !isMonteCarlo  && ptrack1->P() <= 1) return kFALSE;
		}
		
		if(type == 2)
		{
			nsig_p=PIDResponse->NumberOfSigmasTPC(ptrack1,AliPID::kProton);
			if (TMath::Abs(nsig_p) <= cutBetheBloch && cutBetheBloch >0 && !isMonteCarlo  && ptrack1->P() <= 1) return kFALSE;
		}
	
	}
	
	if (pid_n) 
	{
		if(type == 2)
		{
			nsig_n=PIDResponse->NumberOfSigmasTPC(ntrack1,AliPID::kProton);
			if (TMath::Abs(nsig_n) > cutBetheBloch && cutBetheBloch >0 && !isMonteCarlo  && ntrack1->P() <= 1) return kFALSE;
		}
		
		if(type == 1)
		{
			nsig_n=PIDResponse->NumberOfSigmasTPC(ntrack1,AliPID::kProton);
			if (TMath::Abs(nsig_n) <= cutBetheBloch && cutBetheBloch >0 && !isMonteCarlo  && ntrack1->P() <= 1) return kFALSE;
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
	
	// parameters used for most cuts, to minimise editing
	double cutCosPa(0.998), cutcTau(-999);
	double cutNImpact(-999), cutDCA(-999);
	double cutBetheBloch(3); // NOTE - BB cut only applies to data, must be accounted for when constructing corrected yields
	double cutMinNClustersTPC(70), cutRatio(0.8);
	bool isMonteCarlo(false); 
	int isMCtype(0);    //1 = Pure Hijing only, 0 = Anything, -1 = Injected only
	double cutEta(0.8), cutRapidity(0.5), cutArmenteros(0.2);
	
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
	
	if(isMonteCarlo) 
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
			
			Double_t mcXv=0., mcYv=0., mcZv=0.;
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
					
					// this block of code is used to include primary Sigma0 decays as primary lambda/antilambda
					bool sigma0MC = false;
					if(motherType == 3212 || motherType == -3212)// || motherType == 3322 || motherType == -3322 || motherType == 3312 || motherType == -3312)
					{
						if(mcMother->IsPhysicalPrimary())
						{sigma0MC = true;}
					}
					
					Double_t dz=mcZv - mcPart->Zv(), dy=mcYv - mcPart->Yv(), dx=mcXv - mcPart->Xv();
					Double_t mcDecayLength = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
					Double_t mccTau = mcDecayLength*(mcPart->M())/(mcPart->P());
					Double_t mcCosPA = (dx*mcPart->Px()+dy*mcPart->Py()+dz*mcPart->Pz())/(mcDecayLength*mcPart->P());
					
					if(mcPart->IsPhysicalPrimary() || sigma0MC)
					{
						isprimaryMC = true;
						if(lambdaMC)
						{
							if(TMath::Abs(mcPart->Y())<=cutRapidity)
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
								
								fHistMcCosPaAllLaPt->Fill(mcCosPA,mcPart->M()); 
								fHistMccTauAllLaPt->Fill(mccTau,mcPart->M()); 
								fHistMcDcaAllLaPt->Fill(mcDecayLength,mcPart->M()); 
								fHistMcNSigmaAllLaPt->Fill(1.0,mcPart->M()); 
								fHistMcEtaAllLaPt->Fill(mcPart->Eta(),mcPart->M()); 
								fHistMcRapAllLaPt->Fill(mcPart->Y(),mcPart->M()); 
							}
						}
						if(antilambdaMC)
						{
							if(TMath::Abs(mcPart->Y())<=cutRapidity)
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
								
								fHistMcCosPaAllLbPt->Fill(mcCosPA,mcPart->M()); 
								fHistMccTauAllLbPt->Fill(mccTau,mcPart->M()); 
								fHistMcDcaAllLbPt->Fill(mcDecayLength,mcPart->M()); 
								fHistMcNSigmaAllLbPt->Fill(1.0,mcPart->M()); 	
								fHistMcEtaAllLbPt->Fill(mcPart->Eta(),mcPart->M()); 
								fHistMcRapAllLbPt->Fill(mcPart->Y(),mcPart->M()); 	
								
							}
						}
						if(kshortMC)
						{
							if(TMath::Abs(mcPart->Y())<=cutRapidity)
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
								
								fHistMcCosPaAllK0Pt->Fill(mcCosPA,mcPart->M()); 	
								fHistMccTauAllK0Pt->Fill(mccTau,mcPart->M()); 	
								fHistMcDcaAllK0Pt->Fill(mcDecayLength,mcPart->M()); 
								fHistMcNSigmaAllK0Pt->Fill(1.0,mcPart->M()); 	
								fHistMcEtaAllK0Pt->Fill(mcPart->Eta(),mcPart->M()); 
								fHistMcRapAllK0Pt->Fill(mcPart->Y(),mcPart->M()); 
								
							}
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
		
		if(!AcceptV0_general(v0,aod,cutCosPa,cutNImpact,cutDCA,cutEta,cutMinNClustersTPC, cutRatio))
		{ 
			fHistLog->Fill(86);
			continue; 
		}
		
		int nnum = 1, pnum = 0;
		const AliAODTrack *ntracktest=(AliAODTrack *)v0->GetDaughter(nnum);
		if(ntracktest->Charge() > 0){nnum = 0; pnum = 1;}	
		
		const AliAODTrack *ntrack1=(AliAODTrack *)v0->GetDaughter(nnum);
		const AliAODTrack *ptrack1=(AliAODTrack *)v0->GetDaughter(pnum);
		if(ntrack1->Charge()>0)
		{ 
			fHistLog->Fill(55);
		}
		if(ntrack1->Charge()==0)
		{ 
			fHistLog->Fill(50);
		}
		if(ntrack1->Charge()<0)
		{ 
			fHistLog->Fill(45);
		}
		const AliAODPid *pid_p1=ptrack1->GetDetPid();
		const AliAODPid *pid_n1=ntrack1->GetDetPid();
		
		
		/*if(peterCuts)
		{
			const AliAODTrack *ptrack1=(AliAODTrack *)v0->GetDaughter(0);
			if(ntrack1->Pt() < 0.160 || ptrack1->Pt() < 0.160) continue;
			
			Double_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];
			Double_t radius = TMath::Sqrt(r2);
			if(radius <= 0.9 || radius >= 100) continue;
			
		}*/
		
		if(!AcceptV0_particle(v0,1,cutcTau, cutRapidity, decayL))
		{ lambdaCandidate = kFALSE; }
		if(!AcceptV0_particle(v0,2,cutcTau, cutRapidity, decayL))
		{ antilambdaCandidate = kFALSE; }
		if(!AcceptV0_particle(v0,0,cutcTau, cutRapidity, decayL))
		{ k0Candidate = kFALSE; }
		
		if(TMath::Sqrt(v0->Pt2V0())<2)
		{
			if(!AcceptV0_lowpt(v0,fPIDResponse,1,cutBetheBloch,decayL,isMonteCarlo))
			{ lambdaCandidate = kFALSE; }
			if(!AcceptV0_lowpt(v0,fPIDResponse,2,cutBetheBloch,decayL,isMonteCarlo))
			{ antilambdaCandidate = kFALSE; }
			if(!AcceptV0_lowpt(v0,fPIDResponse,0,cutBetheBloch,decayL,isMonteCarlo))
			{ k0Candidate = kFALSE; }
		}
		
		if(lambdaCandidate == kFALSE && antilambdaCandidate == kFALSE && k0Candidate == kFALSE)
		{continue;}
		
		fHistLog->Fill(7);
        fHistPt->Fill(TMath::Sqrt(v0->Pt2V0()));
        fHistEta->Fill(v0->PseudoRapV0());
		
		bool feeddown = false;
		
		if(isMonteCarlo)
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
				if(negAssLabel>=0 && negAssLabel < stack->GetEntriesFast())
				{
				AliAODMCParticle *mcNegPart =(AliAODMCParticle*)stack->UncheckedAt(negAssLabel);
				Int_t v0Label = mcNegPart->GetMother();
				if(v0Label >= 0 && v0Label < stack->GetEntriesFast())
				{
					AliAODMCParticle *mcv0 = (AliAODMCParticle *)stack->UncheckedAt(v0Label);
					passedTests = true;
					
					if ((v0Label >= isInjected && isInjected >= 0 && isMCtype == 1) || (v0Label < isInjected && isInjected >= 0 && isMCtype == -1)) 
					{
						lambdaCandidate = false;
						k0Candidate = false;
						antilambdaCandidate = false;
					}
					
						if(mcv0->GetPdgCode() != kLambda0)
						{
							lambdaCandidate = false;
						}
						if(mcv0->GetPdgCode() != kK0Short)
						{
							k0Candidate = false;
						}
						if(mcv0->GetPdgCode() != kLambda0Bar)
						{
							antilambdaCandidate = false;
						}
					
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
							lambdaCandidate = false;
							k0Candidate = false;
							antilambdaCandidate = false;
						}
						if(motherType == 3212 || motherType == -3212)
						{
							if(mcMother->IsPhysicalPrimary())
							{sigma0MC2 = true;}
						}
						if(motherType == 3322 || motherType == -3322 || motherType == 3312 || motherType == -3312 )
						{
							if(mcMother->IsPhysicalPrimary())
							{feeddown = true;}
						}
					}
					
					if(!sigma0MC2 && !mcv0->IsPhysicalPrimary() && !feeddown)
					{
						lambdaCandidate = false;
						k0Candidate = false;
						antilambdaCandidate = false;
					}
				}
				}
			}
			
			if(passedTests == false)
			{
				lambdaCandidate = false;
				k0Candidate = false;
				antilambdaCandidate = false;
			}
		
		}
		
		double ArmenterosAlpha = v0->Alpha();
		double ArmenterosPt	   = v0->QtProng();
		
		if(lambdaCandidate && centPercentile >= 0.0001 && centPercentile <= 90.0 &&!feeddown )
		{ fHistArmPodLa->Fill(ArmenterosAlpha,ArmenterosPt); }
		if( antilambdaCandidate &&  centPercentile >= 0.0001 && centPercentile <= 90.0 &&!feeddown )
		{ fHistArmPodLb->Fill(ArmenterosAlpha,ArmenterosPt); }
		if(k0Candidate &&  centPercentile >= 0.0001 && centPercentile <= 90.0 &&!feeddown )
		{ fHistArmPodK0->Fill(ArmenterosAlpha,ArmenterosPt); }
		
		if( ArmenterosPt <= TMath::Abs(cutArmenteros*ArmenterosAlpha) && cutArmenteros !=-999 )
		{k0Candidate = false;}
		
		if(lambdaCandidate)
		{
			
			if(!feeddown)
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
			}
			
			if(centPercentile >= 0.0001 && centPercentile <= 90.0)
			{
				if(!feeddown)
				{
				fHistMLaPtCent0090->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());
				fHistBBLaPos->Fill(pid_p1->GetTPCmomentum(),pid_p1->GetTPCsignal());
				fHistBBLaNeg->Fill(pid_n1->GetTPCmomentum(),pid_n1->GetTPCsignal());
				}
				if(feeddown)
				{fHistMcFMLaPt->Fill(TMath::Sqrt(v0->Pt2V0()),v0->MassLambda());}
			}
			if(!feeddown)
			{
			fHistCosPaLaPt->Fill(cosPA,v0->MassLambda()); 
			fHistcTauLaPt->Fill(decayL*(TMath::Sqrt(v0->Pt2V0()))/(v0->P()),v0->MassLambda()); 
			fHistDcaLaPt->Fill(dca,v0->MassLambda()); 	
			fHistNSigmaLaPt->Fill(1.0,v0->MassLambda()); 	
			fHistEtaLaPt->Fill(eta,v0->MassLambda()); 
			fHistRapLaPt->Fill(v0->RapLambda(),v0->MassLambda()); 
			
			
			if(isMonteCarlo) 
			{
				fHistMcCosPaFoundLaPt->Fill(1.0,v0->MassLambda()); 
				fHistMccTauFoundLaPt->Fill(1.0,v0->MassLambda()); 
				fHistMcDcaFoundLaPt->Fill(1.0,v0->MassLambda()); 	
				fHistMcNSigmaFoundLaPt->Fill(1.0,v0->MassLambda()); 
				fHistMcEtaFoundLaPt->Fill(1.0,v0->MassLambda()); 	
				fHistMcRapFoundLaPt->Fill(1.0,v0->MassLambda()); 
			}
			}
			
		}
		if(antilambdaCandidate)
		{
			if(!feeddown)
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
			}
			
			
			fHistCosPaLbPt->Fill(cosPA,v0->MassAntiLambda()); 
			fHistcTauLbPt->Fill(decayL*(v0->MassAntiLambda())/(v0->P()),v0->MassAntiLambda()); 
			fHistDcaLbPt->Fill(dca,v0->MassAntiLambda()); 	
			fHistNSigmaLbPt->Fill(1.0,v0->MassAntiLambda()); 	
			fHistEtaLbPt->Fill(eta,v0->MassAntiLambda()); 	
			fHistRapLbPt->Fill(v0->RapLambda(),v0->MassAntiLambda()); 
			
			if(isMonteCarlo) 
			{
				fHistMcCosPaFoundLbPt->Fill(1.0,v0->MassAntiLambda()); 
				fHistMccTauFoundLbPt->Fill(1.0,v0->MassAntiLambda()); 	
				fHistMcDcaFoundLbPt->Fill(1.0,v0->MassAntiLambda()); 	
				fHistMcNSigmaFoundLbPt->Fill(1.0,v0->MassAntiLambda()); 	
				fHistMcEtaFoundLbPt->Fill(1.0,v0->MassAntiLambda()); 	
				fHistMcRapFoundLbPt->Fill(1.0,v0->MassAntiLambda()); 
			}
			}
		}
		if(k0Candidate)
		{
			if(!feeddown)
			{
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
			}
			
			
			fHistCosPaK0Pt->Fill(cosPA,v0->MassK0Short()); 
			fHistcTauK0Pt->Fill(decayL*(v0->MassK0Short())/(v0->P()),v0->MassK0Short()); 	
			fHistDcaK0Pt->Fill(dca,v0->MassK0Short()); 	
			fHistNSigmaK0Pt->Fill(1.0,v0->MassK0Short()); 	
			fHistEtaK0Pt->Fill(eta,v0->MassK0Short()); 
			fHistRapK0Pt->Fill(v0->RapK0Short(),v0->MassK0Short()); 
			
			if(isMonteCarlo) 
			{
				fHistMcCosPaAFoundK0Pt->Fill(1.0,v0->MassK0Short()); 
				fHistMccTauAFoundK0Pt->Fill(1.0,v0->MassK0Short()); 
				fHistMcDcaAFoundK0Pt->Fill(1.0,v0->MassK0Short()); 
				fHistMcNSigmaAFoundK0Pt->Fill(1.0,v0->MassK0Short()); 
				fHistMcEtaAFoundK0Pt->Fill(1.0,v0->MassK0Short()); 
				fHistMcRapAFoundK0Pt->Fill(1.0,v0->MassK0Short()); 
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
