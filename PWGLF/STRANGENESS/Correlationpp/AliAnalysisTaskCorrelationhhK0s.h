/*Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorrelationhhK0s_H
#define AliAnalysisTaskCorrelationhhK0s_H
class AliPIDResponse;
class AliMultSelection;
class AliAODMCParticle;
class AliCentrality;
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisCorrelationEventCollection.h"
#include "AliEventCuts.h"

class AliAnalysisTaskCorrelationhhK0s : public AliAnalysisTaskSE  
{
 public:
  AliAnalysisTaskCorrelationhhK0s();
  AliAnalysisTaskCorrelationhhK0s(const char *name);
  virtual                 ~AliAnalysisTaskCorrelationhhK0s();

  virtual void            UserCreateOutputObjects();
  virtual void            UserExec(Option_t* option);
  virtual void            Terminate(Option_t* option);
  
  void SetMinPt(Float_t ptmin) {fminPtj = ptmin;}
  void SetMaxPt(Float_t ptmax) {fmaxPtj = ptmax;}
  void SetCorr(Bool_t ishhCorr){fIshhCorr = ishhCorr;}
  void SetMC(Bool_t isMC){fReadMCTruth = isMC;}
  void SetEff(Bool_t isEff){isEfficiency = isEff;}
  void SetHybridMCTruth(Bool_t isHybridMCTr){isHybridMCTruth = isHybridMCTr;}
  void SetEvtToMix(Int_t EvtToMix){fnEventsToMix = EvtToMix;}
  void SetEtaTrigger(Float_t EtaTrigger){fEtaTrigger = EtaTrigger;}
  void SetEtahAssoc(Float_t EtahAssoc){fEtahAssoc = EtahAssoc;}
  void SetEtaV0Assoc(Float_t EtaV0Assoc){fEtaV0Assoc = EtaV0Assoc;}
  void SetFilterBit(Int_t FilterBitValue){fFilterBitValue = FilterBitValue;}
  void SetYear (Int_t year = 2010) { fYear = year;}
  void SetAssocParticle (TString AssocParticle = "K0s") { fV0 = AssocParticle;}

  Float_t GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b );

  void ProcessMCParticles(Bool_t Generated, AliAODTrack* track, Int_t& labelPrimOrSec, Float_t lPercentiles, Bool_t isV0, Double_t ZAtDCA, Float_t PtTriggMax, Bool_t fIshhCorr, Int_t PdgCode[], Int_t ParticleLabel[]);
  Double_t CalculateDeltaTheta( Double_t theta1, Double_t theta2 ); 
  Double_t CalculateDeltaPhi( Double_t phi1, Double_t phi2 ) ; 
  Double_t CalculateDeltaEta( Double_t eta1, Double_t eta2 ) ; 
  double CalculateDPhiStar(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2,Double_t rad);

  //  Double_t ThetaS( Double_t posSftR125[3] )const; 
  //Double_t EtaS( Double_t posSftR125[3] ) const ; 

  void DoPairsh1h2 ( const Float_t lcentrality, int fieldsignf, Double_t lBestPrimaryVtxPos, Double_t ptTriggerMassimo);
  //void DoPairshh ( const Float_t lcentrality, int fieldsignf);

 private:
  TString                 fAnalysisType;                  // "ESD" or "AOD" analysis type
  TString                 fCollidingSystem;               // "pp", "pPb", "PbPb" 
  AliAODEvent*            fAOD;             //! input event
  AliPIDResponse *        fPIDResponse;     //!PID response object 
  //  AliMultSelection *      fMultSelection;   //! 
  AliEventCuts            fEventCuts; //! 
  
  TList*                  fOutputList;      //! output list
  TTree*                  fSignalTree;      //! output tree
  TTree*                  fBkgTree;         //! output tree
  TList*                  fOutputList2;     //! output list
  TList*                  fOutputList3;     //! output list
  TList*                  fOutputList4;     //! output list
  
  AliMCEvent *            fMCEvent;         //!
  Bool_t                  fReadMCTruth;
  Bool_t                  fIshhCorr;
  Bool_t                  isEfficiency;
  Bool_t                  isHybridMCTruth;
  AliAnalysisCorrelationEventCollection ***fEventColl;  //!
  AliAnalysisCorrelationEvent *    fEvt;                //!

  Int_t                    fzVertexBins; 
  Int_t                    fnMultBins;	 
  Int_t                    fMaxFirstMult;
  Int_t                    fMaxSecondMult;
  Int_t                    fnEventsToMix;
  Float_t                  fEtaTrigger;
  Float_t                  fEtahAssoc;
  Float_t                  fEtaV0Assoc;
  Int_t                    fFilterBitValue;
  Int_t                    fYear;

  TH1F*                   fHistPt;                   //! 
  TH1F*                   fHistDCAxym1;              //! 
  TH1F*                   fHistDCAxym2;              //! 
  TH1F*                   fHistDCAzm1;               //! 
  TH1F*                   fHistDCAzm2;               //! 
  TH1F*                   fHistPtV0;                 //! 
  TH1F*                   fHistPthAssoc;             //! 
  TH2F*                   fHistPtTMaxBefAllCfrDataMC; //!
  TH1F*                   fHistPtTMinBefAll;          //! 
  TH1F*                   fHistPtTMinBefAllMC;        //! 
  TH2F*                   fHistPtvsMult;              //! 
  TH2F*                   fHistPtvsMultBefAll;        //! 
  TH2F*                   fHistPtMaxvsMult;           //! 
  TH2F*                   fHistPtMaxvsMultBefAll;     //! 
  TH2F*                   fHistPtMaxvsMultBefAllReco;  //! 
  TH2F*                   fHistPtMaxvsMultBefAllGen; //! 
  TH1F*                   fHistZvertex;               //!
  TH1F*                   fHistFractionSharedTPCClusters;               //!
  TH3F*                   fHistNumberChargedAllEvents; //!
  TH3F*                   fHistNumberChargedNoTrigger; //!
  TH3F*                   fHistNumberChargedTrigger; //!
  TH2F*                   fHist_eta_phi;                //!
  TH2F*                   fHist_eta_phi_PtMax;   	//!
  TH1F*                   fHist_multiplicity; 		//!
  TH1F*                   fHist_multiplicity_EvwTrigger;//!
  TH1F*                   fHistEventMult;   		//!
  TH1F*                   fHistEventV0;   		//!
  TH1F*                   fHistTrack;       		//!
  TH1F*                   fHistTOFBunchCrossing; //!
  TH2F*                   fHistLengthvsCrossedRowsAfterSel; //!
  TH2F*                   fHistLengthvsCrossedRows;       //!
  TH1F*                   fHistIsCommonParton        ;      //!                                                             
  TH3F*                   fHistCommonPartonTrueCasc ;      //!                                                              
  TH3F*                   fHistCommonParton        ;      //!            
  TH3F*                   fHistTriggerComposition;  	  //!
  TH3F*                   fHistTriggerCompositionMCTruth; //! 
  TH3F*                   fHistAssocComposition;  	  //!
  TH3F*                   fHistAssocCompositionMCTruth;   //!
  TH1F*                   fHistTrackAssoc;       	  //!
  TH1F*                   fHistPDG;         		  //!
  TH1F*                   fHistTrackBufferOverflow;       //!  	
  TH2F*                   fHistSecondParticleAll; 	  //!
  TH2F*                   fHistSecondParticleTruthAll; 	  //!
  TH2F*                   fHistSecondParticle; 		  //!
  TH2F*                   fHistSecondParticleTruth; 	  //!
  TH1F*                   fMassV0;          		  //!
  TH2F *                  fDCAxyDaughters;                //!
  TH2F *                  fDCAzDaughters;                 //!
  TH2F *                  fDCAxyDaughtersBis;                //!
  TH2F *                  fDCAzDaughtersBis;                 //!
  TH1F*  fDCAxyPos; //!
  TH1F*  fDCAzPos; //!
  TH1F*  fDCAxyNeg; //!
  TH1F*  fDCAzNeg; //!
  TH2F *                  fHistMultvsV0All; 		  //!
  TH2F *                  fHistMultvsV0AllTruth; 	  //!
  TH2F *                  fHistMultvsV0MCAll; 		  //!
  TH3F*                   fHistTriggerNotLeading; 	  //!
  TH3F*                   fHistTriggerNotLeadingMC;       //!  
  TH2F**                  fHistMassvsPt_tagli;              //!           
  TH2F*                   fHistMultvsTriggerBefAll; 	    //!
  TH2F*                   fHistMultvsTriggerMCTruthBefAll;  //!
  TH2F*                   fHistMultvsTrigger; 		    //!
  TH2F*                   fHistMultvsTriggerMCTruth; 	    //!
  TH1F*                   fHistMassPhoton;  		    //!
  TH1F*                   fHistMass2Photon;  		    //!
  TH2F*                   fHistPtArmvsAlpha;  		    //!
  TH2F*                   fHistPtArmvsAlphaAfterSelection;  //!
  TH2F*                   fHistPtArmvsAlphaAfterPhotonSelection;          //! 
  TH2F*                   fHistPtArmvsAlphaAfterLambdaRejectionSelection; //!  
  TH1F*                   fHistTrigger;                     //!
  TH1F*                   fHistTriggerMCTruth;		    //!
  TH1F*                   fHistTriggerwV0;		    //!
  TH1F*                   fHistTriggerwV0MCTruth;	    //!
  TH2F*                   fHistMultiplicityVsVertexZ; 	    //!
  TH1F*                   fHistTriggervsMult; 		    //!
  TH1F*                   fHistTriggervsMultMC; 	    //!
  TH2F*                   fHistMultiplicityOfMixedEvent;    //!
  TH3F *                  fHistGeneratedTriggerPtPhi; 	    //!
  TH3F **                 fHistSelectedTriggerPtPhi; 	    //!
  TH3F **                 fHistSelectedGenTriggerPtPhi;     //!
  TH3F **                 fHistGeneratedV0PtTMaxPhi; 	    //!
  TH3F **                 fHistCPGeneratedV0PtTMaxPhi; 	    //!
  TH3F **                 fHistSelectedV0PtTMaxPhi; 	    //!
  TH3F *                  fHistGeneratedTriggerPtEta; 	    //!
  TH3F **                 fHistSelectedTriggerPtEta; 	    //!
  TH3F **                 fHistSelectedGenTriggerPtEta;     //!
  TH3F **                 fHistGeneratedV0PtTMaxEta; 	    //!
  TH3F **                 fHistCPGeneratedV0PtTMaxEta; 	    //!
  TH3F **                 fHistSelectedV0PtTMaxEta; 	    //!
  TH3F **                 fHistGeneratedV0PtPtTMax; 	    //!
  TH3F **                 fHistCPGeneratedV0PtPtTMax; 	    //!
  TH3F **                 fHistSelectedV0PtPtTMax; 	    //!
  TH3F **                 fHistSelectedGenV0PtPtTMax; 	    //!
  TH3F **                 fHistGeneratedV0PtEta; 	    //!
  TH3F **                 fHistSelectedV0PtEta; 	    //!
  TH3F *                  fHistReconstructedV0PtMass; 	    //!
  TH3F *                  fHistSelectedV0PtMass;            //!

  TH2F*  fHistTriggerPtRecovsPtGen;         //!
  TH2F*  fHistAssocPtRecovsPtGen; 	    //!
  TH2F*  fHistTriggerPtRecovsPtGenNotPrim;  //!
  TH2F*  fHistAssocPtRecovsPtGenNotPrim;    //!
  TH2F*  fHistTriggerPtRecovsPtGenPion;     //!
  TH2F*  fHistAssocPtRecovsPtGenPion; 	    //!
  TH2F*  fHistTriggerPtRecovsPtGenProton;   //!
  TH2F*  fHistAssocPtRecovsPtGenProton;     //!
  TH2F*  fHistTriggerPtRecovsPtGenKaon;     //!
  TH2F*  fHistAssocPtRecovsPtGenKaon; 	    //!
					   
  TH2F*  fHistResolutionTriggerPt; 	    //!
  TH2F*  fHistResolutionTriggerPhi; 	    //!
  TH2F*  fHistResolutionTriggerEta; 	    //!
  TH2F*  fHistResolutionTriggerPhiPhi; 	    //!
  TH2F*  fHistResolutionTriggerPhiEta; 	    //!
  TH2F*  fHistResolutionV0Pt; 		    //!
  TH2F*  fHistResolutionV0Phi; 		    //!
  TH2F*  fHistResolutionV0PhivsPt; 	    //!
  TH2F*  fHistResolutionV0Eta; 		    //!
  TH2F*  fHistResolutionV0PtvsPt;           //!

  TH2F *** fHistPrimaryTrigger; //!
  TH3F ***  fHistPrimaryV0; //!  
 
  Float_t                 fminPtj;
  Float_t                 fmaxPtj;
  TString                 fV0;
  Float_t                 fminPtV0;
  Float_t                 fmaxPtV0;
  Float_t                 fminPthAssoc;
  Float_t                 fmaxPthAssoc;
  Int_t                   Evcounter;
  Int_t                   Evcounterczero;
  Int_t                   fmolt;
  
  Int_t *                 farrGT;                //!
  UShort_t                fTrackBufferSize;      // Size fo the above array, ~12000 for PbPb

  //tree leaf
  Double_t fTreeVariablePtTrigger;		       
  Int_t    fTreeVariableChargeTrigger;		       
  Double_t fTreeVariableEtaTrigger; 		       
  Double_t fTreeVariablePhiTrigger;		       
  Double_t fTreeVariableDCAz;			       
  Double_t fTreeVariableDCAxy;
  Int_t    fTreeVariableChargeAssoc;			       
  Bool_t   fTreeVariableSkipAssoc;			       
  Bool_t   fTreeVariableIsCommonParton;
  Int_t    fTreeVariablePdgCommonParton;
  Double_t fTreeVariableAssocDCAz;			       
  Double_t fTreeVariableAssocDCAxy;			       
  Double_t fTreeVariableRapK0Short;		       	      
  Int_t    fTreeVariableisPrimaryTrigger;
  Int_t    fTreeVariableisPrimaryV0;
  Double_t fTreeVariableDcaV0ToPrimVertex ;	       	      
  Double_t fTreeVariableDcaPosToPrimVertex;	       	      
  Double_t fTreeVariableDcaNegToPrimVertex;	       	      
  Double_t fTreeVariableV0CosineOfPointingAngle;       	      
  Double_t fTreeVariablePtV0;			       
  Double_t fTreeVariablectau;			       
  Double_t fTreeVariableInvMassK0s;		       
  Double_t fTreeVariableInvMassLambda;		       
  Double_t fTreeVariableInvMassAntiLambda;		       
  Double_t fTreeVariableEtaV0;			       
  Double_t fTreeVariablePhiV0;			       
  Double_t fTreeVariablePtArmenteros;                   
  Double_t fTreeVariableAlpha;	   
  Double_t fTreeVariableDeltaEta;			       
  Double_t fTreeVariableDeltaPhi;			       
  Double_t fTreeVariableDeltaTheta;

  Double_t fTreeVariableMultiplicity;                   
  Double_t fTreeVariableZvertex;
  Int_t    fTreeVariablePDGCodeTrigger;
  Int_t    fTreeVariablePDGCodeAssoc;

  bool FifoShiftok;	                        		      

  AliAnalysisTaskCorrelationhhK0s(const AliAnalysisTaskCorrelationhhK0s&); 
  AliAnalysisTaskCorrelationhhK0s& operator=(const AliAnalysisTaskCorrelationhhK0s&); 

  ClassDef(AliAnalysisTaskCorrelationhhK0s, 1);
};

#endif
