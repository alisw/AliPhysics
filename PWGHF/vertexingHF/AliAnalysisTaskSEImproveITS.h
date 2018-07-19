/* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALI_ANALYSIS_TASK_SE_IMPROVE_ITS_H
#define ALI_ANALYSIS_TASK_SE_IMPROVE_ITS_H

#include "AliAnalysisTaskSE.h"

/// \class AliAnalysisTaskSEImproveITS

class TGraph;
class TList;
class AliAODTrack;
class TClonesArray;
class TObjArray;
class AliESDVertex;
class AliVVertex;

class AliAnalysisTaskSEImproveITS:public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSEImproveITS();
  AliAnalysisTaskSEImproveITS(const char *name,
                              const char *period,
                              const char *systematic,
                              Bool_t isRunInVertexing,
                              Int_t ndebug);

  virtual ~AliAnalysisTaskSEImproveITS();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  //  virtual void Init();
  //  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  //  virtual void Terminate(Option_t *option);
  void SetImproveTracks(Bool_t flag=kTRUE) { fImproveTracks=flag; return; }
  void SetUpdateSecVertCovMat(Bool_t flag=kTRUE) { fUpdateSecVertCovMat=flag; return; }
  void SetUpdateSTCovMatrix(Bool_t opt=kTRUE){fUpdateSTCovMatrix=opt;}
  void SetUpdatePulls(Bool_t opt=kTRUE){fUpdatePulls=opt;}
  void SetMimicData(Bool_t opt=kFALSE){fMimicData=opt;}
  void SetAOD(Bool_t flag=kTRUE) { fIsAOD=flag; return; }

private:
  AliAnalysisTaskSEImproveITS(const AliAnalysisTaskSEImproveITS&);
  AliAnalysisTaskSEImproveITS& operator=(const AliAnalysisTaskSEImproveITS&); 

  /// Helper functions
  Double_t EvalGraph(Double_t x,const TGraph *graph,const TGraph *graphSA=0) const; 
  void SmearTrack(AliVTrack *track, Double_t bz);
  AliESDVertex* RecalculateVertex(const AliVVertex *old,TObjArray *tracks,Double_t bField);
  Int_t PhiBin(Double_t phi) const;

  TGraph *fD0ZResPCur  ; /// old pt dep. d0 res. in z for protons
  TGraph *fD0ZResKCur  ; /// old pt dep. d0 res. in z for kaons
  TGraph *fD0ZResPiCur ; /// old pt dep. d0 res. in z for pions
  TGraph *fD0ZResECur ; /// old pt dep. d0 res. in z for electrons
  TGraph *fD0RPResPCur; /// old pt dep. d0 res. in rphi for protons
  TGraph *fD0RPResKCur; /// old pt dep. d0 res. in rphi for kaons
  TGraph *fD0RPResPiCur; /// old pt dep. d0 res. in rphi for pions
  TGraph *fD0RPResECur; /// old pt dep. d0 res. in rphi for electrons
  TGraph *fD0RPSigmaPullRatioP; /// pt dep. d0 sigma pull MC/data in rphi for protons
  TGraph *fD0RPSigmaPullRatioK; /// pt dep. d0 sigma pull MC/data in rphi for kaons
  TGraph *fD0RPSigmaPullRatioPi; /// pt dep. d0 sigma pull MC/data in rphi for pions
  TGraph *fD0RPSigmaPullRatioE; /// pt dep. d0 sigma pull MC/data in rphi for electrons
  TGraph *fD0RPMeanPCur[2][4]; /// old pt dep. d0 mean. in rphi for protons in 4 phi regions
  TGraph *fD0RPMeanKCur[2][4]; /// old pt dep. d0 mean. in rphi for kaons in 4 phi regions
  TGraph *fD0RPMeanPiCur[2][4]; /// old pt dep. d0 mean. in rphi for pions in 4 phi regions
  TGraph *fD0RPMeanECur[2][4]; /// old pt dep. d0 mean. in rphi for electrons in 4 phi regions
  TGraph *fPt1ResPCur  ; /// old pt dep. 1/pt res. for protons
  TGraph *fPt1ResKCur  ; /// old pt dep. 1/pt res. for kaons
  TGraph *fPt1ResPiCur ; /// old pt dep. 1/pt res. for pions
  TGraph *fPt1ResECur ; /// old pt dep. 1/pt res. for electrons
  TGraph *fD0ZResPUpg  ; /// new pt dep. d0 res. in z for protons
  TGraph *fD0ZResKUpg  ; /// new pt dep. d0 res. in z for kaons
  TGraph *fD0ZResPiUpg ; /// new pt dep. d0 res. in z for pions
  TGraph *fD0ZResEUpg ; /// new pt dep. d0 res. in z for electrons
  TGraph *fD0RPResPUpg; /// new pt dep. d0 res. in rphi for protons
  TGraph *fD0RPResKUpg; /// new pt dep. d0 res. in rphi for kaons
  TGraph *fD0RPResPiUpg; /// new pt dep. d0 res. in rphi for pions
  TGraph *fD0RPResEUpg; /// new pt dep. d0 res. in rphi for electrons
  TGraph *fD0RPMeanPUpg[2][4]; /// new pt dep. d0 mean in rphi for protons in 4 phi regions
  TGraph *fD0RPMeanKUpg[2][4]; /// new pt dep. d0 mean in rphi for kaons in 4 phi regions
  TGraph *fD0RPMeanPiUpg[2][4]; /// new pt dep. d0 mean in rphi for pions in 4 phi regions
  TGraph *fD0RPMeanEUpg[2][4]; /// new pt dep. d0 mean in rphi for electrons in 4 phi regions
  TGraph *fPt1ResPUpg  ; /// new pt dep. 1/pt res. for protons
  TGraph *fPt1ResKUpg  ; /// new pt dep. 1/pt res. for kaons
  TGraph *fPt1ResPiUpg ; /// new pt dep. 1/pt res. for pions
  TGraph *fPt1ResEUpg ; /// new pt dep. 1/pt res. for electrons
  TGraph *fD0ZResPCurSA  ; /// old standalone pt dep. d0 res. in z for protons
  TGraph *fD0ZResKCurSA  ; /// old standalone pt dep. d0 res. in z for kaons
  TGraph *fD0ZResPiCurSA ; /// old standalone pt dep. d0 res. in z for pions
  TGraph *fD0ZResECurSA ; /// old standalone pt dep. d0 res. in z for electrons
  TGraph *fD0RPResPCurSA ; /// old standalone pt dep. d0 res. in rphi for protons
  TGraph *fD0RPResKCurSA ; /// old standalone pt dep. d0 res. in rphi for kaons
  TGraph *fD0RPResPiCurSA; /// old standalone pt dep. d0 res. in rphi for pions
  TGraph *fD0RPResECurSA; /// old standalone pt dep. d0 res. in rphi for electrons
  TGraph *fPt1ResPCurSA  ; /// old standalone pt dep. 1/pt res. for protons
  TGraph *fPt1ResKCurSA  ; /// old standalone pt dep. 1/pt res. for kaons
  TGraph *fPt1ResPiCurSA ; /// old standalone pt dep. 1/pt res. for pions
  TGraph *fPt1ResECurSA ; /// old standalone pt dep. 1/pt res. for electrons
  TGraph *fD0ZResPUpgSA  ; /// new standalone pt dep. d0 res. in z for protons
  TGraph *fD0ZResKUpgSA  ; /// new standalone pt dep. d0 res. in z for kaons
  TGraph *fD0ZResPiUpgSA ; /// new standalone pt dep. d0 res. in z for pions
  TGraph *fD0ZResEUpgSA ; /// new standalone pt dep. d0 res. in z for electrons
  TGraph *fD0RPResPUpgSA ; /// new standalone pt dep. d0 res. in rphi for protons
  TGraph *fD0RPResKUpgSA ; /// new standalone pt dep. d0 res. in rphi for kaons
  TGraph *fD0RPResPiUpgSA; /// new standalone pt dep. d0 res. in rphi for pions
  TGraph *fD0RPResEUpgSA; /// new standalone pt dep. d0 res. in rphi for electrons
  TGraph *fPt1ResPUpgSA  ; /// new standalone pt dep. 1/pt res. for protons
  TGraph *fPt1ResKUpgSA  ; /// new standalone pt dep. 1/pt res. for kaons
  TGraph *fPt1ResPiUpgSA ; /// new standalone pt dep. 1/pt res. for pions
  TGraph *fPt1ResEUpgSA ; /// new standalone pt dep. 1/pt res. for electrons

  Bool_t fRunInVertexing; /// flag to run hybrid task before the vertexingHF task or in standard mode
  Bool_t fImproveTracks; /// this is always kTRUE. kFALSE only if re-running on already improved AODs
  Bool_t fUpdateSecVertCovMat; /// flag to swicth on/off the modification of the sec vert cov matrix
  Bool_t fUpdateSTCovMatrix; /// flag to switch on/off the update of the single track covariance matrix
  Bool_t fUpdatePulls; /// flag to switch on/off the correction of the pulls
  Bool_t fMimicData;
  Bool_t fIsAOD;          // flag to run on AOD 
  TClonesArray *fMCs;      // pointer to AOD MC info
  TList   *fDebugOutput; //!<! collection of debug output
  TNtuple *fDebugNtuple; //!<! debug send on output slot 1
  Float_t *fDebugVars;   //!<! variables to store as degug info 
  Int_t   fNDebug;       /// Max number of debug entries into Ntuple

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEImproveITS,10);
  /// \endcond
};

#endif

