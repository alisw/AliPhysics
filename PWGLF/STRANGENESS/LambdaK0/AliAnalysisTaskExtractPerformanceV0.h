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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- David Dobrigkeit Chinellato
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef ALIANALYSISTASKEXTRACTPERFORMANCEV0_H
#define ALIANALYSISTASKEXTRACTPERFORMANCEV0_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskExtractPerformanceV0 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskExtractPerformanceV0();
  AliAnalysisTaskExtractPerformanceV0(const char *name);
  virtual ~AliAnalysisTaskExtractPerformanceV0();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Double_t MyRapidity(Double_t rE, Double_t rPz) const;
  void CheckChargeV0(AliESDv0 *thisv0);

  void SetIsNuclear           (Bool_t lIsNuclear   = kTRUE ) { fkIsNuclear   = lIsNuclear;   }
  void SetIsLowEnergyPP       (Bool_t lLowEnergyPP = kTRUE ) { fkLowEnergyPP = lLowEnergyPP; }
  void SetUseOnTheFly         (Bool_t lUseOnTheFly = kTRUE ) { fkUseOnTheFly = lUseOnTheFly; }
  
 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
  TList  *fListHistV0;  //! List of Cascade histograms
  TTree  *fTree;              //! Output Tree

  AliPIDResponse *fPIDResponse;     // PID response object

  //Objects Controlling Task Behaviour 
  
  Bool_t fkIsNuclear;   //if true, replace multiplicity est. by centrality (default FALSE) 
  Bool_t fkLowEnergyPP; //if true, skip FASTOnly (default FALSE)
  Bool_t fkUseOnTheFly; //if true, will use On-the-fly V0s instead of Offline V0s (default FALSE)

  //Variables for Tree
   Int_t    fTreeVariablePrimaryStatus;      //!
   Int_t    fTreeVariablePrimaryStatusMother;      //!
   Float_t fTreeVariableChi2V0;             //!
   Float_t fTreeVariableDcaV0Daughters; //!
   Float_t fTreeVariableDcaV0ToPrimVertex; //!
   Float_t fTreeVariableDcaPosToPrimVertex; //!
   Float_t fTreeVariableDcaNegToPrimVertex; //!
   Float_t fTreeVariableV0CosineOfPointingAngle; //!
   Float_t fTreeVariableV0Radius; //!
   Float_t fTreeVariablePt; //!
   Float_t fTreeVariablePtMC; //!
   Float_t fTreeVariableRapK0Short; //!
   Float_t fTreeVariableRapLambda; //!
   Float_t fTreeVariableRapMC; //!
   Float_t fTreeVariableInvMassK0s; //!
   Float_t fTreeVariableInvMassLambda; //!
   Float_t fTreeVariableInvMassAntiLambda; //!
   Float_t fTreeVariableAlphaV0; //!
   Float_t fTreeVariablePtArmV0;//!
   Float_t fTreeVariableNegTotMomentum; //!               
   Float_t fTreeVariablePosTotMomentum; //!
   Float_t fTreeVariableNegTransvMomentum; //!   
   Float_t fTreeVariablePosTransvMomentum; //!
   Float_t fTreeVariableNegTransvMomentumMC; //!   
   Float_t fTreeVariablePosTransvMomentumMC; //!
   
   Float_t fTreeVariableNSigmasPosProton; //!
   Float_t fTreeVariableNSigmasPosPion; //! 
   Float_t fTreeVariableNSigmasNegProton; //!
   Float_t fTreeVariableNSigmasNegPion; //! 

   Float_t fTreeVariablePtMother; //!
   Float_t fTreeVariableV0CreationRadius; //!
   Int_t fTreeVariablePID; //!
   Int_t fTreeVariablePIDPositive; //!
   Int_t fTreeVariablePIDNegative; //!
   Int_t fTreeVariablePIDMother; //!
   Int_t fTreeVariableIndexStatus; //!
   Int_t fTreeVariableIndexStatusMother; //!

   //Note: TDistOverTotMom needs a mass hypothesis to be converted to proper decaylength.
   Float_t fTreeVariableDistOverTotMom;//!

   Float_t fTreeVariablePosEta; //!
   Float_t fTreeVariableNegEta; //!

   Int_t fTreeVariableLeastNbrCrossedRows;//!
   Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;//!
   Int_t fTreeVariableMultiplicity;//!

   TH1F      *fHistV0MultiplicityBeforeTrigSel;              //! V0 multiplicity distribution
   TH1F      *fHistV0MultiplicityForTrigEvt;                 //! V0 multiplicity distribution
   TH1F      *fHistV0MultiplicityForSelEvt;                  //! V0 multiplicity distribution
   TH1F      *fHistV0MultiplicityForSelEvtNoTPCOnly;         //! V0 multiplicity distribution
   TH1F      *fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup; //! V0 multiplicity distribution

   TH1F      *fHistMultiplicityBeforeTrigSel;     //! multiplicity distribution      
   TH1F      *fHistMultiplicityForTrigEvt;        //! multiplicity distribution
   TH1F      *fHistMultiplicity;                  //! multiplicity distribution
   TH1F      *fHistMultiplicityNoTPCOnly;         //! multiplicity distribution
   TH1F      *fHistMultiplicityNoTPCOnlyNoPileup; //! multiplicity distribution

//---> Filled At Analysis Scope

   TH3F      *f3dHistPrimAnalysisPtVsYVsMultLambda;     //! Lambda
   TH3F      *f3dHistPrimAnalysisPtVsYVsMultAntiLambda; //! AntiLambda
   TH3F      *f3dHistPrimAnalysisPtVsYVsMultK0Short;    //! K0Short

//---> Containers for monte carlo information for calculating efficiency! 

   TH3F      *f3dHistPrimRawPtVsYVsMultLambda;     //! Lambda
   TH3F      *f3dHistPrimRawPtVsYVsMultAntiLambda; //! AntiLambda
   TH3F      *f3dHistPrimRawPtVsYVsMultK0Short;    //! K0Short

//---> Filled vs Decay Length

   TH3F      *f3dHistPrimRawPtVsYVsDecayLengthLambda;     //! Lambda
   TH3F      *f3dHistPrimRawPtVsYVsDecayLengthAntiLambda; //! AntiLambda
   TH3F      *f3dHistPrimRawPtVsYVsDecayLengthK0Short;    //! K0Short

//---> Needed for FeedDown Corrections

   TH3F      *f3dHistGenPtVsYVsMultXiMinus;      //! Generated Xi- Distrib
   TH3F      *f3dHistGenPtVsYVsMultXiPlus;       //! Generated Xi+ Distrib

   TH1F      *fHistPVx;                      //! PVx distrib
   TH1F      *fHistPVy;                      //! PVy distrib
   TH1F      *fHistPVz;                      //! PVz distrib
   TH1F      *fHistPVxAnalysis;                      //! PVx distrib
   TH1F      *fHistPVyAnalysis;                      //! PVy distrib
   TH1F      *fHistPVzAnalysis;                      //! PVz distrib
   TH1F      *fHistPVxAnalysisHasHighPtLambda;                      //! PVx distrib
   TH1F      *fHistPVyAnalysisHasHighPtLambda;                      //! PVy distrib
   TH1F      *fHistPVzAnalysisHasHighPtLambda;                      //! PVz distrib

   TH1F      *fHistSwappedV0Counter;                      //! Swapped v0 counter

   AliAnalysisTaskExtractPerformanceV0(const AliAnalysisTaskExtractPerformanceV0&);            // not implemented
   AliAnalysisTaskExtractPerformanceV0& operator=(const AliAnalysisTaskExtractPerformanceV0&); // not implemented
   
   ClassDef(AliAnalysisTaskExtractPerformanceV0, 11);
};

#endif
