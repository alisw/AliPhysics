/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskLukeAOD.h 45956 2010-12-10 12:55:37Z agheata $ */
/* AliAnalysisTaskLukeAOD.h
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 * Edited by Luke Hanratty for AODs
 */
#ifndef ALIANALYSISTASKLUKEAOD_H
#define ALIANALYSISTASKLUKEAOD_H

class TH1F;
class TH2F;
class TList;
class AliPIDResponse;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskLukeAOD : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskLukeAOD();
    AliAnalysisTaskLukeAOD(const char *name);
    virtual ~AliAnalysisTaskLukeAOD();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    
 private:
	TList           *fOutput;        // Output list
	AliPIDResponse	*fPIDResponse;	 // PID
	UInt_t			maskIsSelected; // Physics Selection
	
    TH1F            *fHistPt;        // Pt spectrum
    TH1F            *fHistEta;       // pseudorapidity spectrum
	TH1F			*fHistLog;		 // storage of log variables
	TH1F			*fHistNV0;	 // Number of Tracks per event
	TH1F			*fHistZVertex;	 //	Z coordinate of primary vertex
	TH1F			*fHistMCZVertex;	 //	Z coordinate of MC primary vertex
	TH1F			*fHistCentrality; // Centrality of Events
	
	TH2F			*fHistBBK0Pos;		//PID of the positive daughter of K0 candidates
	TH2F			*fHistBBK0Neg;		//PID of the negative daughter of K0 candidates
	TH2F			*fHistBBLaPos;		//PID of the positive daughter of lambda candidates
	TH2F			*fHistBBLaNeg;		//PID of the negative daughter of lambda candidates
	TH2F			*fHistBBLbPos;		//PID of the positive daughter of antilambda candidates
	TH2F			*fHistBBLbNeg;		//PID of the negative daughter of antilambda candidates
	TH2F			*fHistArmPodK0;		//Armenteros plot for K0 candidates.
	TH2F			*fHistArmPodLa;		//Armenteros plot for Lambda candidates.
	TH2F			*fHistArmPodLb;		//Armenteros plot for Antilambda candidates.
	
	TH2F			*fHistBB3SigProton;	//Bethe Bloch plot of protons @3sigma
	TH2F			*fHistMK0Pt;			//	Mass of 'K0' vs transverse momentum
	TH2F			*fHistMLaPt;			//	Mass of 'Lambda' vs transverse momentum
	TH2F			*fHistMLbPt;			//	Mass of 'Antilambda' vs transverse momentum
	TH2F			*fHistMcPMK0Pt;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC
	TH2F			*fHistMcPMLaPt;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC
	TH2F			*fHistMcPMLbPt;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC
	
	TH2F			*fHistMcFMLaPt;			//	Transverse momentum distribution vs reconstructed Lambda mass of feedown from Xi Lambdas that are detected in MC
	
	TH2F			*fHistMK0PtCent0005;			//	Mass of 'K0' vs transverse momentum for centrality 0-5%
	TH2F			*fHistMLaPtCent0005;			//	Mass of 'Lambda' vs transverse momentum for centrality 0-5%
	TH2F			*fHistMLbPtCent0005;			//	Mass of 'Antilambda' vs transverse momentum for centrality 0-5%
	TH2F			*fHistMcPMK0PtCent0005;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 0-5%
	TH2F			*fHistMcPMLaPtCent0005;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 0-5%
	TH2F			*fHistMcPMLbPtCent0005;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 0-5%
	TH1F			*fHistZVertexCent0005;					//	Z coordinate of primary vertex for centrality 0-5%
	TH1F			*fHistMCZVertexCent0005;				//	Z coordinate of MC primary vertex for centrality 0-5%
	
	TH2F			*fHistMK0PtCent0510;			//	Mass of 'K0' vs transverse momentum for centrality 5-10%
	TH2F			*fHistMLaPtCent0510;			//	Mass of 'Lambda' vs transverse momentum for centrality 5-10%
	TH2F			*fHistMLbPtCent0510;			//	Mass of 'Antilambda' vs transverse momentum for centrality 5-10%
	TH2F			*fHistMcPMK0PtCent0510;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 5-10%
	TH2F			*fHistMcPMLaPtCent0510;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 5-10%
	TH2F			*fHistMcPMLbPtCent0510;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 5-10%
	TH1F			*fHistZVertexCent0510;					//	Z coordinate of primary vertex for centrality 5-10%
	TH1F			*fHistMCZVertexCent0510;				//	Z coordinate of MC primary vertex for centrality 5-10%
	
	
	TH2F			*fHistMK0PtCent1020;			//	Mass of 'K0' vs transverse momentum for centrality 10-20%
	TH2F			*fHistMLaPtCent1020;			//	Mass of 'Lambda' vs transverse momentum for centrality 10-20%
	TH2F			*fHistMLbPtCent1020;			//	Mass of 'Antilambda' vs transverse momentum for centrality 10-20%
	TH2F			*fHistMcPMK0PtCent1020;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 10-20%
	TH2F			*fHistMcPMLaPtCent1020;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 10-20%
	TH2F			*fHistMcPMLbPtCent1020;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 10-20%
	TH1F			*fHistZVertexCent1020;					//	Z coordinate of primary vertex for centrality 10-20%
	TH1F			*fHistMCZVertexCent1020;				//	Z coordinate of MC primary vertex for centrality 10-20%
	
	
	TH2F			*fHistMK0PtCent2040;			//	Mass of 'K0' vs transverse momentum for centrality 20-40%
	TH2F			*fHistMLaPtCent2040;			//	Mass of 'Lambda' vs transverse momentum for centrality 20-40%
	TH2F			*fHistMLbPtCent2040;			//	Mass of 'Antilambda' vs transverse momentum for centrality 20-40%
	TH2F			*fHistMcPMK0PtCent2040;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 20-40%
	TH2F			*fHistMcPMLaPtCent2040;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 20-40%
	TH2F			*fHistMcPMLbPtCent2040;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 20-40%
	TH1F			*fHistZVertexCent2040;					//	Z coordinate of primary vertex for centrality 20-40%
	TH1F			*fHistMCZVertexCent2040;				//	Z coordinate of MC primary vertex for centrality 20-40%
	
	
	TH2F			*fHistMK0PtCent4060;			//	Mass of 'K0' vs transverse momentum for centrality 40-60%
	TH2F			*fHistMLaPtCent4060;			//	Mass of 'Lambda' vs transverse momentum for centrality 40-60%
	TH2F			*fHistMLbPtCent4060;			//	Mass of 'Antilambda' vs transverse momentum for centrality 40-60%
	TH2F			*fHistMcPMK0PtCent4060;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 40-60%
	TH2F			*fHistMcPMLaPtCent4060;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 40-60%
	TH2F			*fHistMcPMLbPtCent4060;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 40-60%
	TH1F			*fHistZVertexCent4060;					//	Z coordinate of primary vertex for centrality 40-60%
	TH1F			*fHistMCZVertexCent4060;				//	Z coordinate of MC primary vertex for centrality 40-60%
	
	
	TH2F			*fHistMK0PtCent6090;			//	Mass of 'K0' vs transverse momentum for centrality 60-90%
	TH2F			*fHistMLaPtCent6090;			//	Mass of 'Lambda' vs transverse momentum for centrality 60-90%
	TH2F			*fHistMLbPtCent6090;			//	Mass of 'Antilambda' vs transverse momentum for centrality 60-90%
	TH2F			*fHistMcPMK0PtCent6090;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 60-90%
	TH2F			*fHistMcPMLaPtCent6090;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 60-90%
	TH2F			*fHistMcPMLbPtCent6090;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 60-90%
	TH1F			*fHistZVertexCent6090;					//	Z coordinate of primary vertex for centrality 60-90%
	TH1F			*fHistMCZVertexCent6090;				//	Z coordinate of MC primary vertex for centrality 60-90%
	
	
	TH2F			*fHistMK0PtCent0090;			//	Mass of 'K0' vs transverse momentum for centrality 0-90%
	TH2F			*fHistMLaPtCent0090;			//	Mass of 'Lambda' vs transverse momentum for centrality 0-90%
	TH2F			*fHistMLbPtCent0090;			//	Mass of 'Antilambda' vs transverse momentum for centrality 0-90%
	TH2F			*fHistMcPMK0PtCent0090;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 0-90%
	TH2F			*fHistMcPMLaPtCent0090;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 0-90%
	TH2F			*fHistMcPMLbPtCent0090;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 0-90%
	TH1F			*fHistZVertexCent0090;					//	Z coordinate of primary vertex for centrality 0-90%
	TH1F			*fHistMCZVertexCent0090;				//	Z coordinate of MC primary vertex for centrality 0-90%
		
	
	TH2F			*fHistCosPaLaPt;		//	Transverse momentum distribution vs CosPa for Lambda Candidates
	TH2F			*fHistCosPaLbPt;		//	Transverse momentum distribution vs CosPa for AntiLambda Candidates
	TH2F			*fHistCosPaK0Pt;		//	Transverse momentum distribution vs CosPa for K0Short Candidates
	TH2F			*fHistMcCosPaAllLaPt;	//	Transverse momentum distribution vs CosPa for all MC primary Lambda
	TH2F			*fHistMcCosPaAllLbPt;	//	Transverse momentum distribution vs CosPa for all MC primary AntiLambda
	TH2F			*fHistMcCosPaAllK0Pt;	//	Transverse momentum distribution vs CosPa for all MC primary K0Short
	TH2F			*fHistMcCosPaFoundLaPt;	//	Transverse momentum distribution vs CosPa for reconstructed MC primary Lambda
	TH2F			*fHistMcCosPaFoundLbPt;	//	Transverse momentum distribution vs CosPa for reconstructed MC primary AntiLambda
	TH2F			*fHistMcCosPaAFoundK0Pt;//	Transverse momentum distribution vs CosPa for reconstructed MC primary K0Short
	
	TH2F			*fHistcTauLaPt;		//	Transverse momentum distribution vs cTau for Lambda Candidates
	TH2F			*fHistcTauLbPt;		//	Transverse momentum distribution vs cTau for AntiLambda Candidates
	TH2F			*fHistcTauK0Pt;		//	Transverse momentum distribution vs cTau for K0Short Candidates
	TH2F			*fHistMccTauAllLaPt;	//	Transverse momentum distribution vs cTau for all MC primary Lambda
	TH2F			*fHistMccTauAllLbPt;	//	Transverse momentum distribution vs cTau for all MC primary AntiLambda
	TH2F			*fHistMccTauAllK0Pt;	//	Transverse momentum distribution vs cTau for all MC primary K0Short
	TH2F			*fHistMccTauFoundLaPt;	//	Transverse momentum distribution vs cTau for reconstructed MC primary Lambda
	TH2F			*fHistMccTauFoundLbPt;	//	Transverse momentum distribution vs cTau for reconstructed MC primary AntiLambda
	TH2F			*fHistMccTauAFoundK0Pt;//	Transverse momentum distribution vs cTau for reconstructed MC primary K0Short
	
	TH2F			*fHistDcaLaPt;		//	Transverse momentum distribution vs Dca for Lambda Candidates
	TH2F			*fHistDcaLbPt;		//	Transverse momentum distribution vs Dca for AntiLambda Candidates
	TH2F			*fHistDcaK0Pt;		//	Transverse momentum distribution vs Dca for K0Short Candidates
	TH2F			*fHistMcDcaAllLaPt;	//	Transverse momentum distribution vs Dca for all MC primary Lambda
	TH2F			*fHistMcDcaAllLbPt;	//	Transverse momentum distribution vs Dca for all MC primary AntiLambda
	TH2F			*fHistMcDcaAllK0Pt;	//	Transverse momentum distribution vs Dca for all MC primary K0Short
	TH2F			*fHistMcDcaFoundLaPt;	//	Transverse momentum distribution vs Dca for reconstructed MC primary Lambda
	TH2F			*fHistMcDcaFoundLbPt;	//	Transverse momentum distribution vs Dca for reconstructed MC primary AntiLambda
	TH2F			*fHistMcDcaAFoundK0Pt;//	Transverse momentum distribution vs Dca for reconstructed MC primary K0Short
	
	TH2F			*fHistNSigmaLaPt;		//	Transverse momentum distribution vs NSigma for Lambda Candidates
	TH2F			*fHistNSigmaLbPt;		//	Transverse momentum distribution vs NSigma for AntiLambda Candidates
	TH2F			*fHistNSigmaK0Pt;		//	Transverse momentum distribution vs NSigma for K0Short Candidates
	TH2F			*fHistMcNSigmaAllLaPt;	//	Transverse momentum distribution vs NSigma for all MC primary Lambda
	TH2F			*fHistMcNSigmaAllLbPt;	//	Transverse momentum distribution vs NSigma for all MC primary AntiLambda
	TH2F			*fHistMcNSigmaAllK0Pt;	//	Transverse momentum distribution vs NSigma for all MC primary K0Short
	TH2F			*fHistMcNSigmaFoundLaPt;	//	Transverse momentum distribution vs NSigma for reconstructed MC primary Lambda
	TH2F			*fHistMcNSigmaFoundLbPt;	//	Transverse momentum distribution vs NSigma for reconstructed MC primary AntiLambda
	TH2F			*fHistMcNSigmaAFoundK0Pt;//	Transverse momentum distribution vs NSigma for reconstructed MC primary K0Short
	
	TH2F			*fHistEtaLaPt;		//	Transverse momentum distribution vs Eta for Lambda Candidates
	TH2F			*fHistEtaLbPt;		//	Transverse momentum distribution vs Eta for AntiLambda Candidates
	TH2F			*fHistEtaK0Pt;		//	Transverse momentum distribution vs Eta for K0Short Candidates
	TH2F			*fHistMcEtaAllLaPt;	//	Transverse momentum distribution vs Eta for all MC primary Lambda
	TH2F			*fHistMcEtaAllLbPt;	//	Transverse momentum distribution vs Eta for all MC primary AntiLambda
	TH2F			*fHistMcEtaAllK0Pt;	//	Transverse momentum distribution vs Eta for all MC primary K0Short
	TH2F			*fHistMcEtaFoundLaPt;	//	Transverse momentum distribution vs Eta for reconstructed MC primary Lambda
	TH2F			*fHistMcEtaFoundLbPt;	//	Transverse momentum distribution vs Eta for reconstructed MC primary AntiLambda
	TH2F			*fHistMcEtaAFoundK0Pt;//	Transverse momentum distribution vs Eta for reconstructed MC primary K0Short
	
	TH2F			*fHistRapLaPt;		//	Transverse momentum distribution vs Rap for Lambda Candidates
	TH2F			*fHistRapLbPt;		//	Transverse momentum distribution vs Rap for AntiLambda Candidates
	TH2F			*fHistRapK0Pt;		//	Transverse momentum distribution vs Rap for K0Short Candidates
	TH2F			*fHistMcRapAllLaPt;	//	Transverse momentum distribution vs Rap for all MC primary Lambda
	TH2F			*fHistMcRapAllLbPt;	//	Transverse momentum distribution vs Rap for all MC primary AntiLambda
	TH2F			*fHistMcRapAllK0Pt;	//	Transverse momentum distribution vs Rap for all MC primary K0Short
	TH2F			*fHistMcRapFoundLaPt;	//	Transverse momentum distribution vs Rap for reconstructed MC primary Lambda
	TH2F			*fHistMcRapFoundLbPt;	//	Transverse momentum distribution vs Rap for reconstructed MC primary AntiLambda
	TH2F			*fHistMcRapAFoundK0Pt;//	Transverse momentum distribution vs Rap for reconstructed MC primary K0Short
	
    // NEW HISTO to be declared here
    
    AliAnalysisTaskLukeAOD(const AliAnalysisTaskLukeAOD&); // not implemented
    AliAnalysisTaskLukeAOD& operator=(const AliAnalysisTaskLukeAOD&); // not implemented
    
    ClassDef(AliAnalysisTaskLukeAOD, 1); // example of analysis
};

#endif

