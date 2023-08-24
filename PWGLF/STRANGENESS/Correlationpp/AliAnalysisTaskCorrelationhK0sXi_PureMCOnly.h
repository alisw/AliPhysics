/*Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCorrelationhK0sXi_PureMCOnly_H
#define AliAnalysisTaskCorrelationhK0sXi_PureMCOnly_H
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisCorrelationEventCollection.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"

class AliAnalysisTaskCorrelationhK0sXi_PureMCOnly : public AliAnalysisTaskSE  
{
 public:
  AliAnalysisTaskCorrelationhK0sXi_PureMCOnly();
  AliAnalysisTaskCorrelationhK0sXi_PureMCOnly(const char *name);
  virtual                 ~AliAnalysisTaskCorrelationhK0sXi_PureMCOnly();

  virtual void            UserCreateOutputObjects();
  virtual void            UserExec(Option_t* option);
  virtual void            Terminate(Option_t* option);
  
  void SetMinPt(Float_t ptmin) {fminPtj = ptmin;}
  void SetMaxPt(Float_t ptmax) {fmaxPtj = ptmax;}
  void SetCorr(Bool_t ishhCorr){fIshhCorr = ishhCorr;}
  void SetEvtToMix(Int_t EvtToMix){fnEventsToMix = EvtToMix;}
  void SetEtaTrigger(Float_t EtaTrigger){fEtaTrigger = EtaTrigger;}
  void SetEtahAssoc(Float_t EtahAssoc){fEtahAssoc = EtahAssoc;}
  void SetEtaV0Assoc(Float_t EtaV0Assoc){fEtaV0Assoc = EtaV0Assoc;}
  void SetFilterBit(Int_t FilterBitValue){fFilterBitValue = FilterBitValue;}
  void SetYear (Int_t year = 2010) { fYear = year;}
  void SetInclusive (Bool_t isInclusive = 1) {fisInclusiveINELgt0 = isInclusive;}
  void SetHM (Bool_t isHM) { fisHM = isHM;}
  void SetMinimumMultPercentile (Float_t PercentilesMin) { lPercentilesMin = PercentilesMin;}
  void SetMaximumMultPercentile (Float_t PercentilesMax) { lPercentilesMax = PercentilesMax;}
  void SetAssocParticle (TString AssocParticle = "K0s") { fV0 = AssocParticle;}

  Double_t CalculateDeltaTheta( Double_t theta1, Double_t theta2 ); 
  Double_t CalculateDeltaPhi( Double_t phi1, Double_t phi2 ) ; 
  Double_t CalculateDeltaEta( Double_t eta1, Double_t eta2 ) ; 
  double CalculateDPhiStar(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2,Double_t rad);

  //  Double_t ThetaS( Double_t posSftR125[3] )const; 
  //Double_t EtaS( Double_t posSftR125[3] ) const ; 

  void DoPairsh1h2 ( const Float_t lcentrality, const Float_t NchVZEROA, const Float_t lNchVZEROC, int fieldsignf, Double_t lBestPrimaryVtxPos, Double_t ptTriggerMassimo);
  //void DoPairshh ( const Float_t lcentrality, int fieldsignf);

 private:
  TString                 fAnalysisType;                  // "ESD" or "AOD" analysis type
  TString                 fCollidingSystem;               // "pp", "pPb", "PbPb" 
  
  TList*                  fOutputList;      //! output list
  TTree*                  fSignalTree;      //! output tree
  TTree*                  fBkgTree;         //! output tree
  TList*                  fOutputList2;     //! output list
  TList*                  fOutputList3;     //! output list
  TList*                  fOutputList4;     //! output list
  
  AliMCEvent *            fMCEvent;         //!
  Bool_t                  fIshhCorr;
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
  Bool_t                   fisInclusiveINELgt0;
  Float_t                  lPercentilesMin;
  Float_t                  lPercentilesMax;
  Bool_t                   fisHM;

  TH2F*                   fHistInelCorr;                    //!
  TH2F*                   fHistEvtNoTrigger;                //! 
  TH1F*                   fHistPt;                          //! 
  TH1F*                   fHistPtV0;                        //! 
  TH2F*                   fHistPtTMaxBefAllCfrDataMC;       //!
  TH1F*                   fHistPtTMinBefAll;                //! 
  TH1F*                   fHistPtTMinBefAllMC;              //! 
  TH2F*                   fHistPtvsMult;                    //! 
  TH2F*                   fHistPtvsMultBefAll;              //! 
  TH2F*                   fHistPtMaxvsMult;                 //! 
  TH2F*                   fHistPtMaxvsMultBefAll;           //! 
  TH1F*                   fHistZvertex;                     //!
  TH2F*                   fHistMultForwardvsMidRap;         //!
  TH3F*                   fHistMultForwardvsMidRapvsPt;     //!
  TH1F*                   fHistFractionSharedTPCClusters;   //!
  TH3F*                   fHistNumberChargedAllEvents;      //!
  TH3F*                   fHistNumberChargedNoTrigger;      //!
  TH3F*                   fHistNumberChargedTrigger;        //!
  TH2F*                   fHist_eta_phi;                    //!
  TH2F*                   fHist_eta_phi_PtMax;   	    //!
  TH1F*                   fHist_multiplicityAllSelEvents;   //!
  TH1F*                   fHist_multiplicity; 		    //!
  TH1F*                   fHist_multiplicity_EvwTrigger;    //!
  TH1F*                   fHistEventMult;   		    //!
  TH1F*                   fHistEventV0;   		    //!
  TH2F*                   fHistEventV0Pt;   		    //!
  TH1F*                   fHistTrack;       		    //!
  TH1F*                   fHistIsCommonParton        ;      //!                                                             
  TH3F*                   fHistCommonParton        ;        //!            
  TH3F*                   fHistCommonPartonTrueCasc ;       //!                                                              
  TH3F*                   fHistTriggerComposition;  	    //!
  TH3F*                   fHistTriggerCompositionMCTruth;   //! 
  TH3F*                   fHistAssocComposition;  	    //!
  TH3F*                   fHistAssocCompositionMCTruth;     //!
  TH1F*                   fHistTrackAssoc;       	    //!
  TH1F*                   fHistPDG;         		    //!
  TH1F*                   fHistTrackBufferOverflow;         //!  	
  TH2F*                   fHistSecondParticleAll; 	    //!
  TH2F*                   fHistSecondParticleTruthAll; 	    //!
  TH2F*                   fHistSecondParticle; 		    //!
  TH2F*                   fHistSecondParticleTruth; 	    //!
  TH2F*                   fHistMultvsV0All; 		    //!
  TH2F*                   fHistMultvsV0AllTruth; 	    //!
  TH2F*                   fHistMultvsV0MCAll; 		    //!
  TH3F*                   fHistTriggerNotLeading; 	    //!
  TH3F*                   fHistTriggerNotLeadingMC;         //!  
  TH2F*                   fHistMultvsTriggerBefAll; 	    //!
  TH2F*                   fHistMultvsTriggerMCTruthBefAll;  //!
  TH2F*                   fHistMultvsTrigger; 		    //!
  TH2F*                   fHistMultvsTriggerMCTruth; 	    //!
  TH1F*                   fHistTrigger;                     //!
  TH1F*                   fHistTriggerMCTruth;		    //!
  TH1F*                   fHistTriggerwV0;		    //!
  TH1F*                   fHistTriggerwV0MCTruth;	    //!
  TH2F*                   fHistMultiplicityVsVertexZ; 	    //!
  TH1F*                   fHistTriggervsMult; 		    //!
  TH1F*                   fHistTriggervsMultMC; 	    //!
  TH2F*                   fHistMultiplicityOfMixedEvent;    //!
  TH2F* 	          fHistGeneratedV0Pt;               //!
  TH3F *                  fHistGeneratedTriggerPtPhi; 	    //!
  TH3F *                  fHistGeneratedV0PtTMaxPhi; 	    //!
  TH3F *                  fHistGeneratedTriggerPtEta; 	    //!
  TH3F *                  fHistGeneratedV0PtTMaxEta; 	    //!
  TH3F *                  fHistGeneratedV0PtPtTMax; 	    //!
  TH3F *                  fHistGeneratedV0PtEta; 	    //!

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
  Int_t    fTreeVariableisPrimaryTrigger;
  Int_t    fTreeVariableisPrimaryV0;
  Double_t fTreeVariableRapK0Short;		       	      
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
  Double_t fTreeVariableMultiplicityV0A;
  Double_t fTreeVariableMultiplicityV0C;
  Double_t fTreeVariableZvertex;
  Int_t    fTreeVariablePDGCodeTrigger;
  Int_t    fTreeVariablePDGCodeAssoc;
  Int_t    fTreeVariableLabelTrigger;
  Int_t    fTreeVariableLabelPos;
  Int_t    fTreeVariableLabelNeg;

  bool FifoShiftok;	                        		      

  AliAnalysisTaskCorrelationhK0sXi_PureMCOnly(const AliAnalysisTaskCorrelationhK0sXi_PureMCOnly&); 
  AliAnalysisTaskCorrelationhK0sXi_PureMCOnly& operator=(const AliAnalysisTaskCorrelationhK0sXi_PureMCOnly&); 

  ClassDef(AliAnalysisTaskCorrelationhK0sXi_PureMCOnly, 1);
};

#endif
