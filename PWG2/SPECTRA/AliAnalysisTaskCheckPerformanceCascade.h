#ifndef ALIANALYSISTASKCHECKPERFORMANCECASCADE_H
#define ALIANALYSISTASKCHECKPERFORMANCECASCADE_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//		 AliAnalysisTaskCheckPerformanceCascade class
//            This task is for a performance study of cascade identification.
//            It works with MC info and ESD/AOD tree.
//            Origin : A.Maire Mar2009, antonin.maire@ires.in2p3.fr
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckPerformanceCascade : public AliAnalysisTaskSE {
 public:
 
  AliAnalysisTaskCheckPerformanceCascade();
  AliAnalysisTaskCheckPerformanceCascade(const char *name );
  virtual ~AliAnalysisTaskCheckPerformanceCascade() {}
  
  //virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetDebugLevelCascade(Int_t lDebugCascade = 0)          {fDebugCascade = lDebugCascade;}
  void SetCollidingSystems (Short_t collidingSystems = 0)     {fCollidingSystems = collidingSystems;}
  void SetAnalysisType     (const char* analysisType = "ESD") {fAnalysisType = analysisType;}
  
 private:
	Int_t	fDebugCascade;			// Denug Flag for this task devoted to cascade
	TString fAnalysisType;			// "ESD" or "AOD" analysis type	
	Short_t fCollidingSystems;		// 0 = pp collisions or 1 = AA collisions
	
 	TList	*fListHistCascade;		//! List of Cascade histograms
		// - Histos 
	TH1F	*fHistMCTrackMultiplicity;	//! MC Track multiplicity
	
	
	// proton
	TH1F	*fHistEtaGenProton;   			//! MC Pseudo-rapidity of any generated p+ (no cuts in acceptance)
	TH1F	*fHistEtaGenAntiProton;   		//! MC Pseudo-rapidity of any generated p+ (no cuts in acceptance)

// Part 1 - Generated cascades
	
	//--------------
	// Xi-
	TH1F	*fHistEtaGenCascXiMinus;   		//! MC Pseudo-rapidity of any generated Xi- (no cuts in acceptance)
	
	// - Plots for |y(MC)| < 1
	TH1F	*fHistYGenCascMidRapXiMinus;		//! MC rapidity of Xi- generated within |y(MC)| < 1
	TH1F	*fHistEtaGenCascMidRapXiMinus;	//! MC Eta of Xi-      generated within |y(MC)| < 1
	TH1F	*fHistThetaGenCascMidRapXiMinus;	//! MC theta of Xi-    generated within |y(MC)| < 1
	TH1F	*fHistPtGenCascMidRapXiMinus;		//! MC Pt of Xi-       generated within |y(MC)| < 1
	
	// - Histos planned for Xi- emitted within the acceptance (cuts in theta + pt of daughters)
	// 	= findable cascades
	TH1F	*fHistThetaGenCascXiMinus;		//! MC Theta angle of the generated Xi-
	TH1F	*fHistPtFdblGenCascXiMinus;		//! MC Pt of the cascade, Xi-
	
	TH1F	*fHistThetaLambdaXiMinus;		//! MC Theta angle of the Lambda daughter of the generated Xi-
	TH1F	*fHistThetaBachXiMinus;		//! MC Theta angle of the Bachelor (pi-)
	
	TH1F	*fHistThetaMesDghterXiMinus;		//! MC Theta angle of the mesonic  V0 daughter in the generated cascade, pi-
	TH1F	*fHistThetaBarDghterXiMinus;		//! MC Theta angle of the baryonic V0 daughter in the generated cascade, p+
	
	TH1F	*fHistPtBachXiMinus;   		//! MC Pt of the Bachelor (pi-)                   (Control Plot)
	TH1F	*fHistPtMesDghterXiMinus;		//! MC Pt of the meson daughter of the 'Lambda0', pi-   (Control Plot)
	TH1F	*fHistPtBarDghterXiMinus;		//! MC Pt of the baryon daughter of the 'Lambda0', p+   (Control Plot)
	
	
	
	//--------------
	// Xi+
	TH1F	*fHistEtaGenCascXiPlus;   		//! MC Pseudo-rapidity of any generated Xi+ (no cuts in acceptance)
	
	// - Plots for |y(MC)| < 1
	TH1F	*fHistYGenCascMidRapXiPlus;		//! MC rapidity of Xi+ generated within |y(MC)| < 1
	TH1F	*fHistEtaGenCascMidRapXiPlus;		//! MC Eta of Xi+      generated within |y(MC)| < 1
	TH1F	*fHistThetaGenCascMidRapXiPlus;	//! MC theta of Xi+    generated within |y(MC)| < 1
	TH1F	*fHistPtGenCascMidRapXiPlus;		//! MC Pt of Xi+       generated within |y(MC)| < 1
		
	// - Histos planned for Xi+ emitted within the acceptance (cuts in theta + pt of daughters)
	TH1F	*fHistThetaGenCascXiPlus;		//! MC Theta angle of the generated Xi+
	TH1F	*fHistPtFdblGenCascXiPlus;		//! MC Pt of the cascade, Xi+
	
	TH1F	*fHistThetaLambdaXiPlus;		//! MC Theta angle of the anti-Lambda daughter of the generated Xi+
	TH1F	*fHistThetaBachXiPlus;		//! MC Theta angle of the Bachelor (pi+)
	
	TH1F	*fHistThetaMesDghterXiPlus;		//! MC Theta angle of the mesonic  V0 daughter in the generated cascade, pi+
	TH1F	*fHistThetaBarDghterXiPlus;		//! MC Theta angle of the baryonic V0 daughter in the generated cascade, p-
	
	TH1F	*fHistPtBachXiPlus;   		//! MC Pt of the Bachelor (pi+)                   (Control Plot)
	TH1F	*fHistPtMesDghterXiPlus;		//! MC Pt of the meson daughter of the 'Lambda0', pi+ (Control Plot)
	TH1F	*fHistPtBarDghterXiPlus;		//! MC Pt of the baryon daughter of the 'Lambda0', p- (Control Plot)
	
	
	
	//--------------
	// Omega-
	TH1F	*fHistEtaGenCascOmegaMinus;   	//! MC Pseudo-rapidity of any generated Omega- (no cuts in acceptance)
	
	// - Plots for |y(MC)| < 1
	TH1F	*fHistYGenCascMidRapOmegaMinus;	//! MC rapidity of Omega- generated within |y(MC)| < 1
	TH1F	*fHistEtaGenCascMidRapOmegaMinus;	//! MC Eta of Omega-      generated within |y(MC)| < 1
	TH1F	*fHistThetaGenCascMidRapOmegaMinus;	//! MC theta of Omega-    generated within |y(MC)| < 1
	TH1F	*fHistPtGenCascMidRapOmegaMinus;	//! MC Pt of Omega-       generated within |y(MC)| < 1
	
	// - Histos planned for Omega- emitted within the acceptance (cuts in theta + pt of daughters)
	TH1F	*fHistThetaGenCascOmegaMinus;		//! MC Theta angle of the generated Omega-
	TH1F	*fHistPtFdblGenCascOmegaMinus;	//! MC Pt of the cascade, Omega-
	
	TH1F	*fHistThetaLambdaOmegaMinus;		//! MC Theta angle of the Lambda daughter of the generated Omega-
	TH1F	*fHistThetaBachOmegaMinus;		//! MC Theta angle of the Bachelor (K-)
	
	TH1F	*fHistThetaMesDghterOmegaMinus;	//! MC Theta angle of the mesonic  V0 daughter in the generated cascade, pi-
	TH1F	*fHistThetaBarDghterOmegaMinus;	//! MC Theta angle of the baryonic V0 daughter in the generated cascade, p+
	
	TH1F	*fHistPtBachOmegaMinus;   		//! MC Pt of the Bachelor (K-)                   (Control Plot)
	TH1F	*fHistPtMesDghterOmegaMinus;		//! MC Pt of the meson daughter of the 'Lambda0', pi- (Control Plot)
	TH1F	*fHistPtBarDghterOmegaMinus;		//! MC Pt of the baryon daughter of the 'Lambda0', p+ (Control Plot)
	
	
	
	//--------------
	// Omega+
	TH1F	*fHistEtaGenCascOmegaPlus;   		//! MC Pseudo-rapidity of any generated Omega+ (no cuts in acceptance)
		
	// - Plots for |y(MC)| < 1
	TH1F	*fHistYGenCascMidRapOmegaPlus;	//! MC rapidity of Omega+ generated within |y(MC)| < 1
	TH1F	*fHistEtaGenCascMidRapOmegaPlus;	//! MC Eta of Omega+      generated within |y(MC)| < 1
	TH1F	*fHistThetaGenCascMidRapOmegaPlus;	//! MC theta of Omega+    generated within |y(MC)| < 1
	TH1F	*fHistPtGenCascMidRapOmegaPlus;	//! MC Pt of Omega+       generated within |y(MC)| < 1
	
	// - Histos planned for Omega+ emitted within the acceptance (cuts in theta + pt of daughters)
	TH1F	*fHistThetaGenCascOmegaPlus;		//! MC Theta angle of the generated Omega+
	TH1F	*fHistPtFdblGenCascOmegaPlus;		//! MC Pt of the cascade, Omega+
	
	TH1F	*fHistThetaLambdaOmegaPlus;		//! MC Theta angle of the anti-Lambda daughter of the generated Omega+
	TH1F	*fHistThetaBachOmegaPlus;		//! MC Theta angle of the Bachelor (K+)
	
	TH1F	*fHistThetaMesDghterOmegaPlus;	//! MC Theta angle of the mesonic  V0 daughter in the generated cascade, pi+
	TH1F	*fHistThetaBarDghterOmegaPlus;	//! MC Theta angle of the baryonic V0 daughter in the generated cascade, p-
	
	TH1F	*fHistPtBachOmegaPlus;   		//! MC Pt of the Bachelor (K+) (Control Plot)
	TH1F	*fHistPtMesDghterOmegaPlus;		//! MC Pt of the meson daughter of the 'Lambda0', pi+ (Control Plot)
	TH1F	*fHistPtBarDghterOmegaPlus;		//! MC Pt of the baryon daughter of the 'Lambda0', p- (Control Plot)
	
	
	
// Part 2 - Any reconstructed cascades + reconstructed cascades associated with MC
	// - Effective mass histos for all the cascade candidates
	TH1F	*fHistMassXiMinus;			//! reconstructed cascade effective mass, under Xi- hyp.
	TH1F	*fHistMassXiPlus;			//! reconstructed cascade effective mass, under Xi+ hyp.
	TH1F	*fHistMassOmegaMinus;			//! reconstructed cascade effective mass, under Omega- hyp.
	TH1F	*fHistMassOmegaPlus;			//! reconstructed cascade effective mass, under Omega+ hyp.
	
	// - Effective mass histos for the cascade candidates associated with MC
	TH1F	*fHistAsMCMassXiMinus;			//! reconstr. cascade effective mass, under Xi- hyp. for Associated cand.
	TH1F	*fHistAsMCMassXiPlus;			//! reconstr. cascade effective mass, under Xi+ hyp. for Associated cand.
	TH1F	*fHistAsMCMassOmegaMinus;		//! reconstr. cascade effective mass, under Omega- hyp. for Associated cand.
	TH1F	*fHistAsMCMassOmegaPlus;		//! reconstr. cascade effective mass, under Omega+ hyp. for Associated cand.
	
	// - Generated Pt of the cascade candidates associated with MC
	TH1F	*fHistAsMCGenPtXiMinus;			//! generated Pt from the MC Xi- associated with a reconstr. cascade
	TH1F	*fHistAsMCGenPtXiPlus;			//! generated Pt from the MC Xi+ associated with a reconstr. cascade
	TH1F	*fHistAsMCGenPtOmegaMinus;		//! generated Pt from the MC Omega- associated with a reconstr. cascade
	TH1F	*fHistAsMCGenPtOmegaPlus;		//! generated Pt from the MC Omega+ associated with a reconstr. cascade
  
	// - Generated Y of the cascade candidates associated with MC
	TH1F	*fHistAsMCGenYXiMinus;			//! generated Rap. from the MC Xi- associated with a reconstr. cascade
	TH1F	*fHistAsMCGenYXiPlus;			//! generated Rap. from the MC Xi+ associated with a reconstr. cascade
	TH1F	*fHistAsMCGenYOmegaMinus;		//! generated Rap. from the MC Omega- associated with a reconstr. cascade
	TH1F	*fHistAsMCGenYOmegaPlus;		//! generated Rap. from the MC Omega+ associated with a reconstr. cascade
	
	// - Generated Y Vs Generated Pt, for the cascade candidates associated with MC
	TH2F	*f2dHistAsMCGenYVsGenPtXiMinus;		//! gen. Rap. Vs gen. Pt from the MC Xi- associated with a reconstr. cascade
	TH2F	*f2dHistAsMCGenYVsGenPtXiPlus;		//! gen. Rap. Vs gen. Pt from the MC Xi+ associated with a reconstr. cascade
	TH2F	*f2dHistAsMCGenYVsGenPtOmegaMinus;	//! gen. Rap. Vs gen. Pt from the MC Omega- associated with a reconstr. cascade
	TH2F	*f2dHistAsMCGenYVsGenPtOmegaPlus;	//! gen. Rap. Vs gen. Pt from the MC Omega+ associated with a reconstr. cascade
	
	// - Generated Eta of the the cascade candidates associated with MC
	TH1F	*fHistAsMCGenEtaXiMinus;		//! generated Eta from the MC Xi- associated with a reconstr. cascade
	TH1F	*fHistAsMCGenEtaXiPlus;			//! generated Eta from the MC Xi+ associated with a reconstr. cascade
	TH1F	*fHistAsMCGenEtaOmegaMinus;		//! generated Eta from the MC Omega- associated with a reconstr. cascade
	TH1F	*fHistAsMCGenEtaOmegaPlus;		//! generated Eta from the MC Omega+ associated with a reconstr. cascade
	
	// - Resolution in Pt as function of generated Pt
	TH2F	*f2dHistAsMCResPtXiMinus;		//! resolution in Pt as function of gen. Pt, for Xi-
	TH2F	*f2dHistAsMCResPtXiPlus;		//! resolution in Pt as function of gen. Pt, for Xi-
	TH2F	*f2dHistAsMCResPtOmegaMinus;		//! resolution in Pt as function of gen. Pt, for Omega-
	TH2F	*f2dHistAsMCResPtOmegaPlus;		//! resolution in Pt as function of gen. Pt, for Omega+
	
	// - Resolution in R(2D) as function of generated R
	TH2F	*f2dHistAsMCResRXiMinus;		//! resolution in transv. R as function of transv. gen. R, for Xi-
	TH2F	*f2dHistAsMCResRXiPlus;			//! resolution in transv. R as function of transv. gen. R, for Xi+
	TH2F	*f2dHistAsMCResROmegaMinus;		//! resolution in transv. R as function of transv. gen. R, for Omega-
	TH2F	*f2dHistAsMCResROmegaPlus;		//! resolution in transv. R as function of transv. gen. R, for Omega+
	
	
  AliAnalysisTaskCheckPerformanceCascade(const AliAnalysisTaskCheckPerformanceCascade&);            // not implemented
  AliAnalysisTaskCheckPerformanceCascade& operator=(const AliAnalysisTaskCheckPerformanceCascade&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckPerformanceCascade, 1);
};

#endif
