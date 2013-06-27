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

#ifndef AliAnalysisTaskExtractPerformanceV0pPb_H
#define AliAnalysisTaskExtractPerformanceV0pPb_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskExtractPerformanceV0pPb : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskExtractPerformanceV0pPb();
  AliAnalysisTaskExtractPerformanceV0pPb(const char *name);
  virtual ~AliAnalysisTaskExtractPerformanceV0pPb();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Double_t MyRapidity(Double_t rE, Double_t rPz) const;
 
  void SetDiffractiveOnly ( Bool_t lDiffractiveOnly = kTRUE ) { fDiffractiveOnly = lDiffractiveOnly; }

 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
  TList  *fListHistV0;  //! List of Cascade histograms
  TTree  *fTree;              //! Output Tree, V0
  TTree  *fTreeEvents;        //! Output Tree for events

  AliPIDResponse *fPIDResponse;     // PID response object

  //Objects Controlling Task Behaviour 
  Bool_t fDiffractiveOnly; //Only look at diffractive generated events

//===========================================================================================
//   Variables for Tree, V0s
//===========================================================================================
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

   //Note: TDistOverTotMom needs a mass hypothesis to be converted to proper decaylength.
   Float_t fTreeVariableDistOverTotMom;//!

   Float_t fTreeVariablePosEta; //!
   Float_t fTreeVariableNegEta; //!

   Float_t fTreeVariableVertexZ; //!

   Int_t fTreeVariableLeastNbrCrossedRows;//!
   Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;//!
  
  Float_t fTreeVariableCentrality; //! 

  //Variables for Event Tree
  Float_t fTreeEventsCentrality; //! 

//===========================================================================================
//   Histograms
//===========================================================================================

  //Default V0A Centrality only 
   TH1F      *fHistCentralityProcessed; //! All processed
   TH1F      *fHistCentralityTrigEvt;   //! Those selected by PS / Trigger
   TH1F      *fHistCentralityHasVtx;    //! Those that have a well-established vertex
   TH1F      *fHistCentralityVtxZ;      //! Those whose vertex falls within |z|<10cm

//---> Filled At Analysis Scope

//---> Strategy for yCMS shift: do binning in "400,-1,1" combination in rapidity. 
//     This will allow for post-processing shifts with a granularity of 0.005 to 
//     be applied as needed, and note the shift we need is 0.465 so it fits this 
//     binning. 

//---> At final Analysis Event selection
   TH3F      *f3dHist_Analysis_PtVsYVsV0A_Lambda;      //! Lambda
   TH3F      *f3dHist_Analysis_PtVsYVsV0A_AntiLambda;  //! AntiLambda
   TH3F      *f3dHist_Analysis_PtVsYVsV0A_K0Short;     //! K0Short

//---> At generator level (before any event selection except centrality binning)
//     ...though note that centrality is already a selection. 
   TH3F      *f3dHist_Generated_PtVsYVsV0A_Lambda;      //! Lambda
   TH3F      *f3dHist_Generated_PtVsYVsV0A_AntiLambda;  //! AntiLambda
   TH3F      *f3dHist_Generated_PtVsYVsV0A_K0Short;     //! K0Short

//Cross-checking histograms: Charged Kaons (to compare with neutral ones at generator level) 
  TH3F      *f3dHist_Analysis_PtVsYVsV0A_KPlus;     //! Added for cross-check of any bias
  TH3F      *f3dHist_Analysis_PtVsYVsV0A_KMinus;    //! Added for cross-check of any bias
  TH3F      *f3dHist_Generated_PtVsYVsV0A_KPlus;     //! Added for cross-check of any bias
  TH3F      *f3dHist_Generated_PtVsYVsV0A_KMinus;    //! Added for cross-check of any bias

//---> Needed for FeedDown Corrections

//V0A Containers

  TH3F      *f3dHist_Analysis_PtVsYVsV0A_XiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHist_Analysis_PtVsYVsV0A_XiPlus;       //! Generated Xi+ Distrib

  TH3F      *f3dHist_Generated_PtVsYVsV0A_XiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHist_Generated_PtVsYVsV0A_XiPlus;       //! Generated Xi+ Distrib

//With YCMS shift... 

//---> At final Analysis Event selection
   TH3F      *f3dHist_Analysis_PtVsYCMSVsV0A_Lambda;      //! Lambda
   TH3F      *f3dHist_Analysis_PtVsYCMSVsV0A_AntiLambda;  //! AntiLambda
   TH3F      *f3dHist_Analysis_PtVsYCMSVsV0A_K0Short;     //! K0Short

//---> At generator level (before any event selection except centrality binning)
//     ...though note that centrality is already a selection. 
   TH3F      *f3dHist_Generated_PtVsYCMSVsV0A_Lambda;      //! Lambda
   TH3F      *f3dHist_Generated_PtVsYCMSVsV0A_AntiLambda;  //! AntiLambda
   TH3F      *f3dHist_Generated_PtVsYCMSVsV0A_K0Short;     //! K0Short

//Cross-checking histograms: Charged Kaons (to compare with neutral ones at generator level) 
  TH3F      *f3dHist_Analysis_PtVsYCMSVsV0A_KPlus;     //! Added for cross-check of any bias
  TH3F      *f3dHist_Analysis_PtVsYCMSVsV0A_KMinus;    //! Added for cross-check of any bias
  TH3F      *f3dHist_Generated_PtVsYCMSVsV0A_KPlus;     //! Added for cross-check of any bias
  TH3F      *f3dHist_Generated_PtVsYCMSVsV0A_KMinus;    //! Added for cross-check of any bias

//---> Needed for FeedDown Corrections

//V0A Containers

  TH3F      *f3dHist_Analysis_PtVsYCMSVsV0A_XiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHist_Analysis_PtVsYCMSVsV0A_XiPlus;       //! Generated Xi+ Distrib

  TH3F      *f3dHist_Generated_PtVsYCMSVsV0A_XiMinus;      //! Generated Xi- Distrib
  TH3F      *f3dHist_Generated_PtVsYCMSVsV0A_XiPlus;       //! Generated Xi+ Distrib

  
   AliAnalysisTaskExtractPerformanceV0pPb(const AliAnalysisTaskExtractPerformanceV0pPb&);            // not implemented
   AliAnalysisTaskExtractPerformanceV0pPb& operator=(const AliAnalysisTaskExtractPerformanceV0pPb&); // not implemented
   
   ClassDef(AliAnalysisTaskExtractPerformanceV0pPb, 11);
};

#endif
