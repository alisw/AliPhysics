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
		
	void	SetIsMonteCarlo			(bool isMonteCarlo = false)			{fIsMonteCarlo = isMonteCarlo;}
	void	SetCutCosPa				(double cutCosPa = 0.998)			{fcutCosPa = cutCosPa;}
	void	SetCutcTauMin			(double cutcTauMin = -999)			{fcutcTauMin = cutcTauMin;}
	void	SetCutNcTauMax			(double cutNcTauMax = 3.0)			{fcutNcTauMax = cutNcTauMax;}
	void	SetCutBetheBloch		(double cutBetheBloch = 3.0)		{fcutBetheBloch = cutBetheBloch;}
	void	SetCutMinNClustersTPC	(double cutMinNClustersTPC = 70)	{fcutMinNClustersTPC = cutMinNClustersTPC;}
	void	SetCutRatio				(double cutRatio = 0.8)				{fcutRatio = cutRatio;}
	void	SetCutEta				(double cutEta = 0.8)				{fcutEta = cutEta;}
	void	SetCutRapidity			(double cutRapidity = 0.5)			{fcutRapidity = cutRapidity;}
	void	SetCutArmenteros		(double cutArmenteros = 0.2)		{fcutArmenteros = cutArmenteros;}
	
 private:
	bool	fIsMonteCarlo;
	
	double	fcutCosPa;
	double	fcutcTauMin;
	double	fcutNcTauMax;
	double	fcutBetheBloch;
	double	fcutMinNClustersTPC;
	double	fcutRatio;
	double	fcutEta;
	double	fcutRapidity;
	double	fcutArmenteros;
	
	TList           *fOutput;        // Output list
	AliPIDResponse	*fPIDResponse;	 // PID
	//UInt_t			maskIsSelected; // Physics Selection
	
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
	TH2F			*fHistMcAsMK0PtCent0005;			//	Transverse momentum distribution vs reconstructed K0 mass of primary reconstructed K0s in MC for centrality 0-5%
	TH2F			*fHistMcAsMLaPtCent0005;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary reconstructed Lambda in MC for centrality 0-5%
	TH2F			*fHistMcAsMLbPtCent0005;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary reconstructed Antilambda in MC for centrality 0-5%
	TH1F			*fHistZVertexCent0005;					//	Z coordinate of primary vertex for centrality 0-5%
	TH1F			*fHistMCZVertexCent0005;				//	Z coordinate of MC primary vertex for centrality 0-5%
	
	TH2F			*fHistMK0PtCent0510;			//	Mass of 'K0' vs transverse momentum for centrality 5-10%
	TH2F			*fHistMLaPtCent0510;			//	Mass of 'Lambda' vs transverse momentum for centrality 5-10%
	TH2F			*fHistMLbPtCent0510;			//	Mass of 'Antilambda' vs transverse momentum for centrality 5-10%
	TH2F			*fHistMcPMK0PtCent0510;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 5-10%
	TH2F			*fHistMcPMLaPtCent0510;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 5-10%
	TH2F			*fHistMcPMLbPtCent0510;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 5-10%
	TH2F			*fHistMcAsMK0PtCent0510;			//	Transverse momentum distribution vs reconstructed K0 mass of primary reconstructed K0s in MC for centrality 5-10%
	TH2F			*fHistMcAsMLaPtCent0510;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary reconstructed Lambda in MC for centrality 5-10%
	TH2F			*fHistMcAsMLbPtCent0510;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary reconstructed Antilambda in MC for centrality 5-10%
	TH1F			*fHistZVertexCent0510;					//	Z coordinate of primary vertex for centrality 5-10%
	TH1F			*fHistMCZVertexCent0510;				//	Z coordinate of MC primary vertex for centrality 5-10%
	
	
	TH2F			*fHistMK0PtCent1020;			//	Mass of 'K0' vs transverse momentum for centrality 10-20%
	TH2F			*fHistMLaPtCent1020;			//	Mass of 'Lambda' vs transverse momentum for centrality 10-20%
	TH2F			*fHistMLbPtCent1020;			//	Mass of 'Antilambda' vs transverse momentum for centrality 10-20%
	TH2F			*fHistMcPMK0PtCent1020;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 10-20%
	TH2F			*fHistMcPMLaPtCent1020;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 10-20%
	TH2F			*fHistMcPMLbPtCent1020;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 10-20%
	TH2F			*fHistMcAsMK0PtCent1020;			//	Transverse momentum distribution vs reconstructed K0 mass of primary reconstructed K0s in MC for centrality 10-20%
	TH2F			*fHistMcAsMLaPtCent1020;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary reconstructed Lambda in MC for centrality 10-20%
	TH2F			*fHistMcAsMLbPtCent1020;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary reconstructed Antilambda in MC for centrality 10-20%
	TH1F			*fHistZVertexCent1020;					//	Z coordinate of primary vertex for centrality 10-20%
	TH1F			*fHistMCZVertexCent1020;				//	Z coordinate of MC primary vertex for centrality 10-20%
	
	
	TH2F			*fHistMK0PtCent2040;			//	Mass of 'K0' vs transverse momentum for centrality 20-40%
	TH2F			*fHistMLaPtCent2040;			//	Mass of 'Lambda' vs transverse momentum for centrality 20-40%
	TH2F			*fHistMLbPtCent2040;			//	Mass of 'Antilambda' vs transverse momentum for centrality 20-40%
	TH2F			*fHistMcPMK0PtCent2040;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 20-40%
	TH2F			*fHistMcPMLaPtCent2040;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 20-40%
	TH2F			*fHistMcPMLbPtCent2040;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 20-40%
	TH2F			*fHistMcAsMK0PtCent2040;			//	Transverse momentum distribution vs reconstructed K0 mass of primary reconstructed K0s in MC for centrality 20-40%
	TH2F			*fHistMcAsMLaPtCent2040;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary reconstructed Lambda in MC for centrality 20-40%
	TH2F			*fHistMcAsMLbPtCent2040;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary reconstructed Antilambda in MC for centrality 20-40%
	TH1F			*fHistZVertexCent2040;					//	Z coordinate of primary vertex for centrality 20-40%
	TH1F			*fHistMCZVertexCent2040;				//	Z coordinate of MC primary vertex for centrality 20-40%
	
	
	TH2F			*fHistMK0PtCent4060;			//	Mass of 'K0' vs transverse momentum for centrality 40-60%
	TH2F			*fHistMLaPtCent4060;			//	Mass of 'Lambda' vs transverse momentum for centrality 40-60%
	TH2F			*fHistMLbPtCent4060;			//	Mass of 'Antilambda' vs transverse momentum for centrality 40-60%
	TH2F			*fHistMcPMK0PtCent4060;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 40-60%
	TH2F			*fHistMcPMLaPtCent4060;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 40-60%
	TH2F			*fHistMcPMLbPtCent4060;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 40-60%
	TH2F			*fHistMcAsMK0PtCent4060;			//	Transverse momentum distribution vs reconstructed K0 mass of primary reconstructed K0s in MC for centrality 40-60%
	TH2F			*fHistMcAsMLaPtCent4060;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary reconstructed Lambda in MC for centrality 40-60%
	TH2F			*fHistMcAsMLbPtCent4060;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary reconstructed Antilambda in MC for centrality 40-60%
	TH1F			*fHistZVertexCent4060;					//	Z coordinate of primary vertex for centrality 40-60%
	TH1F			*fHistMCZVertexCent4060;				//	Z coordinate of MC primary vertex for centrality 40-60%
	
	
	TH2F			*fHistMK0PtCent6090;			//	Mass of 'K0' vs transverse momentum for centrality 60-90%
	TH2F			*fHistMLaPtCent6090;			//	Mass of 'Lambda' vs transverse momentum for centrality 60-90%
	TH2F			*fHistMLbPtCent6090;			//	Mass of 'Antilambda' vs transverse momentum for centrality 60-90%
	TH2F			*fHistMcPMK0PtCent6090;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 60-90%
	TH2F			*fHistMcPMLaPtCent6090;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 60-90%
	TH2F			*fHistMcPMLbPtCent6090;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 60-90%
	TH2F			*fHistMcAsMK0PtCent6090;			//	Transverse momentum distribution vs reconstructed K0 mass of primary reconstructed K0s in MC for centrality 60-90%
	TH2F			*fHistMcAsMLaPtCent6090;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary reconstructed Lambda in MC for centrality 60-90%
	TH2F			*fHistMcAsMLbPtCent6090;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary reconstructed Antilambda in MC for centrality 60-90%
	TH1F			*fHistZVertexCent6090;					//	Z coordinate of primary vertex for centrality 60-90%
	TH1F			*fHistMCZVertexCent6090;				//	Z coordinate of MC primary vertex for centrality 60-90%
	
	
	TH2F			*fHistMK0PtCent0090;			//	Mass of 'K0' vs transverse momentum for centrality 0-90%
	TH2F			*fHistMLaPtCent0090;			//	Mass of 'Lambda' vs transverse momentum for centrality 0-90%
	TH2F			*fHistMLbPtCent0090;			//	Mass of 'Antilambda' vs transverse momentum for centrality 0-90%
	TH2F			*fHistMcPMK0PtCent0090;			//	Transverse momentum distribution vs reconstructed K0 mass of primary K0s in MC for centrality 0-90%
	TH2F			*fHistMcPMLaPtCent0090;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary Lambda in MC for centrality 0-90%
	TH2F			*fHistMcPMLbPtCent0090;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary Antilambda in MC for centrality 0-90%
	TH2F			*fHistMcAsMK0PtCent0090;			//	Transverse momentum distribution vs reconstructed K0 mass of primary reconstructed K0s in MC for centrality 0-90%
	TH2F			*fHistMcAsMLaPtCent0090;			//	Transverse momentum distribution vs reconstructed Lambda mass of primary reconstructed Lambda in MC for centrality 0-90%
	TH2F			*fHistMcAsMLbPtCent0090;			//	Transverse momentum distribution vs reconstructed Antilambd mass of primary reconstructed Antilambda in MC for centrality 0-90%
	TH1F			*fHistZVertexCent0090;					//	Z coordinate of primary vertex for centrality 0-90%
	TH1F			*fHistMCZVertexCent0090;				//	Z coordinate of MC primary vertex for centrality 0-90%
		
	
	TH2F			*fHistCosPaMLa;		//	Transverse momentum distribution vs CosPa for Lambda Candidates
	TH2F			*fHistCosPaMLb;		//	Transverse momentum distribution vs CosPa for AntiLambda Candidates
	TH2F			*fHistCosPaMK0;		//	Transverse momentum distribution vs CosPa for K0Short Candidates
	TH2F			*fHistMcGenCosPaMLa;	//	Transverse momentum distribution vs MC-Truth CosPa for all MC primary Lambda
	TH2F			*fHistMcGenCosPaMLb;	//	Transverse momentum distribution vs MC-Truth CosPa for all MC primary AntiLambda
	TH2F			*fHistMcGenCosPaMK0;	//	Transverse momentum distribution vs MC-Truth CosPa for all MC primary K0Short
	TH2F			*fHistMcAsReconCosPaMLa;	//	Transverse momentum distribution vs CosPa for reconstructed MC primary Lambda
	TH2F			*fHistMcAsReconCosPaMLb;	//	Transverse momentum distribution vs CosPa for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsReconCosPaMK0;//	Transverse momentum distribution vs CosPa for reconstructed MC primary K0Short
	TH2F			*fHistMcAsTruthCosPaMLa;	//	Transverse momentum distribution vs MC-Truth CosPa for reconstructed MC primary Lambda
	TH2F			*fHistMcAsTruthCosPaMLb;	//	Transverse momentum distribution vs MC-Truth CosPa for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsTruthCosPaMK0;//	Transverse momentum distribution vs MC-Truth CosPa for reconstructed MC primary K0Short
	
	TH2F			*fHistcTauMLa;		//	Transverse momentum distribution vs cTau for Lambda Candidates
	TH2F			*fHistcTauMLb;		//	Transverse momentum distribution vs cTau for AntiLambda Candidates
	TH2F			*fHistcTauMK0;		//	Transverse momentum distribution vs cTau for K0Short Candidates
	TH2F			*fHistMcGencTauMLa;	//	Transverse momentum distribution vs MC-Truth cTau for all MC primary Lambda
	TH2F			*fHistMcGencTauMLb;	//	Transverse momentum distribution vs MC-Truth cTau for all MC primary AntiLambda
	TH2F			*fHistMcGencTauMK0;	//	Transverse momentum distribution vs MC-Truth cTau for all MC primary K0Short
	TH2F			*fHistMcAsReconcTauMLa;	//	Transverse momentum distribution vs cTau for reconstructed MC primary Lambda
	TH2F			*fHistMcAsReconcTauMLb;	//	Transverse momentum distribution vs cTau for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsReconcTauMK0;//	Transverse momentum distribution vs cTau for reconstructed MC primary K0Short
	TH2F			*fHistMcAsTruthcTauMLa;	//	Transverse momentum distribution vs MC-Truth cTau for reconstructed MC primary Lambda
	TH2F			*fHistMcAsTruthcTauMLb;	//	Transverse momentum distribution vs MC-Truth cTau for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsTruthcTauMK0;//	Transverse momentum distribution vs MC-Truth cTau for reconstructed MC primary K0Short
	
	TH2F			*fHistDcaMLa;		//	Transverse momentum distribution vs Dca for Lambda Candidates
	TH2F			*fHistDcaMLb;		//	Transverse momentum distribution vs Dca for AntiLambda Candidates
	TH2F			*fHistDcaMK0;		//	Transverse momentum distribution vs Dca for K0Short Candidates
	TH2F			*fHistMcGenDcaMLa;	//	Transverse momentum distribution vs MC-Truth Dca for all MC primary Lambda
	TH2F			*fHistMcGenDcaMLb;	//	Transverse momentum distribution vs MC-Truth Dca for all MC primary AntiLambda
	TH2F			*fHistMcGenDcaMK0;	//	Transverse momentum distribution vs MC-Truth Dca for all MC primary K0Short
	TH2F			*fHistMcAsReconDcaMLa;	//	Transverse momentum distribution vs Dca for reconstructed MC primary Lambda
	TH2F			*fHistMcAsReconDcaMLb;	//	Transverse momentum distribution vs Dca for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsReconDcaMK0;//	Transverse momentum distribution vs Dca for reconstructed MC primary K0Short
	TH2F			*fHistMcAsTruthDcaMLa;	//	Transverse momentum distribution vs MC-Truth Dca for reconstructed MC primary Lambda
	TH2F			*fHistMcAsTruthDcaMLb;	//	Transverse momentum distribution vs MC-Truth Dca for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsTruthDcaMK0;//	Transverse momentum distribution vs MC-Truth Dca for reconstructed MC primary K0Short
	
	TH2F			*fHistNSigmaMLa;		//	Transverse momentum distribution vs NSigma for Lambda Candidates
	TH2F			*fHistNSigmaMLb;		//	Transverse momentum distribution vs NSigma for AntiLambda Candidates
	TH2F			*fHistNSigmaMK0;		//	Transverse momentum distribution vs NSigma for K0Short Candidates
	TH2F			*fHistMcGenNSigmaMLa;	//	Transverse momentum distribution vs MC-Truth NSigma for all MC primary Lambda
	TH2F			*fHistMcGenNSigmaMLb;	//	Transverse momentum distribution vs MC-Truth NSigma for all MC primary AntiLambda
	TH2F			*fHistMcGenNSigmaMK0;	//	Transverse momentum distribution vs MC-Truth NSigma for all MC primary K0Short
	TH2F			*fHistMcAsReconNSigmaMLa;	//	Transverse momentum distribution vs NSigma for reconstructed MC primary Lambda
	TH2F			*fHistMcAsReconNSigmaMLb;	//	Transverse momentum distribution vs NSigma for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsReconNSigmaMK0;//	Transverse momentum distribution vs NSigma for reconstructed MC primary K0Short
	TH2F			*fHistMcAsTruthNSigmaMLa;	//	Transverse momentum distribution vs MC-Truth NSigma for reconstructed MC primary Lambda
	TH2F			*fHistMcAsTruthNSigmaMLb;	//	Transverse momentum distribution vs MC-Truth NSigma for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsTruthNSigmaMK0;//	Transverse momentum distribution vs MC-Truth NSigma for reconstructed MC primary K0Short
	
	TH2F			*fHistEtaMLa;		//	Transverse momentum distribution vs Eta for Lambda Candidates
	TH2F			*fHistEtaMLb;		//	Transverse momentum distribution vs Eta for AntiLambda Candidates
	TH2F			*fHistEtaMK0;		//	Transverse momentum distribution vs Eta for K0Short Candidates
	TH2F			*fHistMcGenEtaMLa;	//	Transverse momentum distribution vs MC-Truth Eta for all MC primary Lambda
	TH2F			*fHistMcGenEtaMLb;	//	Transverse momentum distribution vs MC-Truth Eta for all MC primary AntiLambda
	TH2F			*fHistMcGenEtaMK0;	//	Transverse momentum distribution vs MC-Truth Eta for all MC primary K0Short
	TH2F			*fHistMcAsReconEtaMLa;	//	Transverse momentum distribution vs Eta for reconstructed MC primary Lambda
	TH2F			*fHistMcAsReconEtaMLb;	//	Transverse momentum distribution vs Eta for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsReconEtaMK0;//	Transverse momentum distribution vs Eta for reconstructed MC primary K0Short
	TH2F			*fHistMcAsTruthEtaMLa;	//	Transverse momentum distribution vs MC-Truth Eta for reconstructed MC primary Lambda
	TH2F			*fHistMcAsTruthEtaMLb;	//	Transverse momentum distribution vs MC-Truth Eta for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsTruthEtaMK0;//	Transverse momentum distribution vs MC-Truth Eta for reconstructed MC primary K0Short
	
	TH2F			*fHistRapMLa;		//	Transverse momentum distribution vs Rap for Lambda Candidates
	TH2F			*fHistRapMLb;		//	Transverse momentum distribution vs Rap for AntiLambda Candidates
	TH2F			*fHistRapMK0;		//	Transverse momentum distribution vs Rap for K0Short Candidates
	TH2F			*fHistMcGenRapMLa;	//	Transverse momentum distribution vs MC-Truth Rap for all MC primary Lambda
	TH2F			*fHistMcGenRapMLb;	//	Transverse momentum distribution vs MC-Truth Rap for all MC primary AntiLambda
	TH2F			*fHistMcGenRapMK0;	//	Transverse momentum distribution vs MC-Truth Rap for all MC primary K0Short
	TH2F			*fHistMcAsReconRapMLa;	//	Transverse momentum distribution vs Rap for reconstructed MC primary Lambda
	TH2F			*fHistMcAsReconRapMLb;	//	Transverse momentum distribution vs Rap for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsReconRapMK0;//	Transverse momentum distribution vs Rap for reconstructed MC primary K0Short
	TH2F			*fHistMcAsTruthRapMLa;	//	Transverse momentum distribution vs MC-Truth Rap for reconstructed MC primary Lambda
	TH2F			*fHistMcAsTruthRapMLb;	//	Transverse momentum distribution vs MC-Truth Rap for reconstructed MC primary AntiLambda
	TH2F			*fHistMcAsTruthRapMK0;//	Transverse momentum distribution vs MC-Truth Rap for reconstructed MC primary K0Short
	
	
	TH2F			*fHistArmPodK0;		//Armenteros plot for K0 candidates.
	TH2F			*fHistArmPodLa;		//Armenteros plot for Lambda candidates.
	TH2F			*fHistArmPodLb;		//Armenteros plot for Antilambda candidates.
	TH2F			*fHistMcGenArmPodK0;		//Armenteros plot for K0 candidates.
	TH2F			*fHistMcGenArmPodLa;		//Armenteros plot for Lambda candidates.
	TH2F			*fHistMcGenArmPodLb;		//Armenteros plot for Antilambda candidates.
	TH2F			*fHistMcAsReconArmPodK0;		//Armenteros plot for K0 candidates.
	TH2F			*fHistMcAsReconArmPodLa;		//Armenteros plot for Lambda candidates.
	TH2F			*fHistMcAsReconArmPodLb;		//Armenteros plot for Antilambda candidates.
	TH2F			*fHistMcAsTruthArmPodK0;		//Armenteros plot for K0 candidates.
	TH2F			*fHistMcAsTruthArmPodLa;		//Armenteros plot for Lambda candidates.
	TH2F			*fHistMcAsTruthArmPodLb;		//Armenteros plot for Antilambda candidates.
	
	
    // NEW HISTO to be declared here
    
    AliAnalysisTaskLukeAOD(const AliAnalysisTaskLukeAOD&); // not implemented
    AliAnalysisTaskLukeAOD& operator=(const AliAnalysisTaskLukeAOD&); // not implemented
    
    ClassDef(AliAnalysisTaskLukeAOD, 1); // example of analysis
};

#endif

