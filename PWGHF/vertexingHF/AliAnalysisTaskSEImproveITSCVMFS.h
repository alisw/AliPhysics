/* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALI_ANALYSIS_TASK_SE_IMPROVE_ITS_CVMFS_H
#define ALI_ANALYSIS_TASK_SE_IMPROVE_ITS_CVMFS_H

#include "AliAnalysisTaskSE.h"

/// \class AliAnalysisTaskSEImproveITSCVMFS

class TGraph;
class TList;
class AliAODTrack;
class TClonesArray;
class TObjArray;
class AliESDVertex;
class AliVVertex;

class TNtuple;

class AliAnalysisTaskSEImproveITSCVMFS:public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSEImproveITSCVMFS();
  AliAnalysisTaskSEImproveITSCVMFS(const char *name,
                              const char *period,
                              const char *systematic,
                              Bool_t isRunInVertexing,
                              Int_t ndebug);

  virtual ~AliAnalysisTaskSEImproveITSCVMFS();

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
  void SetSmearOnlySignal(Bool_t opt=kTRUE) {fSmearOnlySignal=opt;}
  void SetImproverSuffix(TString suffix="central") {fImproverSuffix=suffix;}

private:
  AliAnalysisTaskSEImproveITSCVMFS(const AliAnalysisTaskSEImproveITSCVMFS&);
  AliAnalysisTaskSEImproveITSCVMFS& operator=(const AliAnalysisTaskSEImproveITSCVMFS&); 
  /// Helper functions
  Double_t EvalGraph(Double_t x,const TGraph *graph,const TGraph *graphSA=0) const; 
  void SmearTrack(AliVTrack *track, Double_t bz);
  AliESDVertex* RecalculateVertex(const AliVVertex *old,TObjArray *tracks,Double_t bField);
  Int_t PhiBin(Double_t phi) const;
  void OpenImproverHistos(AliVEvent *event);

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

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //  Specific stuff for PbPb 2018 periods                                   //
  //                                                                         //
  //  The correction is supposed to be different for tracks                  //
  //  satisfying the following SPD requirements:                             //
  //    1) kFirst                                                            //
  //    2) kOnlySecond                                                       //
  //  The correction for the mean is then performed in 24 bins of φ,         //
  //  instead of quarters.                                                   //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  Bool_t fIsPbPb2018;  // flag declaring wheter the PbPb 2018 periods are analysed
  // kFirst
  TGraph *fD0ZResPCur_PbPb2018_kFirst; /// 
  TGraph *fD0ZResKCur_PbPb2018_kFirst; /// 
  TGraph *fD0ZResPiCur_PbPb2018_kFirst; /// 
  TGraph *fD0ZResECur_PbPb2018_kFirst; /// 
  TGraph *fD0RPResPCur_PbPb2018_kFirst; /// 
  TGraph *fD0RPResKCur_PbPb2018_kFirst; /// 
  TGraph *fD0RPResPiCur_PbPb2018_kFirst; /// 
  TGraph *fD0RPResECur_PbPb2018_kFirst; ///
  TGraph *fD0RPSigmaPullRatioP_PbPb2018_kFirst; /// 
  TGraph *fD0RPSigmaPullRatioK_PbPb2018_kFirst; /// 
  TGraph *fD0RPSigmaPullRatioPi_PbPb2018_kFirst; /// 
  TGraph *fD0RPSigmaPullRatioE_PbPb2018_kFirst; ///
  TGraph *fD0RPMeanPCur_PbPb2018_kFirst[2][24]; ///
  TGraph *fD0RPMeanKCur_PbPb2018_kFirst[2][24]; ///
  TGraph *fD0RPMeanPiCur_PbPb2018_kFirst[2][24]; ///
  TGraph *fD0RPMeanECur_PbPb2018_kFirst[2][24]; ///
  TGraph *fPt1ResPCur_PbPb2018_kFirst; /// 
  TGraph *fPt1ResKCur_PbPb2018_kFirst; /// 
  TGraph *fPt1ResPiCur_PbPb2018_kFirst; /// 
  TGraph *fPt1ResECur_PbPb2018_kFirst; /// 
  TGraph *fD0ZResPUpg_PbPb2018_kFirst; /// 
  TGraph *fD0ZResKUpg_PbPb2018_kFirst; /// 
  TGraph *fD0ZResPiUpg_PbPb2018_kFirst; /// 
  TGraph *fD0ZResEUpg_PbPb2018_kFirst; /// 
  TGraph *fD0RPResPUpg_PbPb2018_kFirst; /// 
  TGraph *fD0RPResKUpg_PbPb2018_kFirst; /// 
  TGraph *fD0RPResPiUpg_PbPb2018_kFirst; /// 
  TGraph *fD0RPResEUpg_PbPb2018_kFirst; /// 
  TGraph *fD0RPMeanPUpg_PbPb2018_kFirst[2][24]; ///
  TGraph *fD0RPMeanKUpg_PbPb2018_kFirst[2][24]; ///
  TGraph *fD0RPMeanPiUpg_PbPb2018_kFirst[2][24]; ///
  TGraph *fD0RPMeanEUpg_PbPb2018_kFirst[2][24]; ///
  TGraph *fPt1ResPUpg_PbPb2018_kFirst; /// 
  TGraph *fPt1ResKUpg_PbPb2018_kFirst; /// 
  TGraph *fPt1ResPiUpg_PbPb2018_kFirst; ///
  TGraph *fPt1ResEUpg_PbPb2018_kFirst; /// 
  TGraph *fD0ZResPCurSA_PbPb2018_kFirst; /// 
  TGraph *fD0ZResKCurSA_PbPb2018_kFirst; ///
  TGraph *fD0ZResPiCurSA_PbPb2018_kFirst; /// 
  TGraph *fD0ZResECurSA_PbPb2018_kFirst; /// 
  TGraph *fD0RPResPCurSA_PbPb2018_kFirst; /// 
  TGraph *fD0RPResKCurSA_PbPb2018_kFirst; /// 
  TGraph *fD0RPResPiCurSA_PbPb2018_kFirst; /// 
  TGraph *fD0RPResECurSA_PbPb2018_kFirst; /// 
  TGraph *fPt1ResPCurSA_PbPb2018_kFirst; /// 
  TGraph *fPt1ResKCurSA_PbPb2018_kFirst; /// 
  TGraph *fPt1ResPiCurSA_PbPb2018_kFirst; /// 
  TGraph *fPt1ResECurSA_PbPb2018_kFirst; /// 
  TGraph *fD0ZResPUpgSA_PbPb2018_kFirst; /// 
  TGraph *fD0ZResKUpgSA_PbPb2018_kFirst; /// 
  TGraph *fD0ZResPiUpgSA_PbPb2018_kFirst; /// 
  TGraph *fD0ZResEUpgSA_PbPb2018_kFirst; /// 
  TGraph *fD0RPResPUpgSA_PbPb2018_kFirst; /// 
  TGraph *fD0RPResKUpgSA_PbPb2018_kFirst; /// 
  TGraph *fD0RPResPiUpgSA_PbPb2018_kFirst; /// 
  TGraph *fD0RPResEUpgSA_PbPb2018_kFirst; /// 
  TGraph *fPt1ResPUpgSA_PbPb2018_kFirst; ///
  TGraph *fPt1ResKUpgSA_PbPb2018_kFirst; /// 
  TGraph *fPt1ResPiUpgSA_PbPb2018_kFirst; ///
  TGraph *fPt1ResEUpgSA_PbPb2018_kFirst; /// 
  // kOnlySecond
  TGraph *fD0ZResPCur_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResKCur_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResPiCur_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResECur_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResPCur_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResKCur_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResPiCur_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResECur_PbPb2018_kOnlySecond; ///
  TGraph *fD0RPSigmaPullRatioP_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPSigmaPullRatioK_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPSigmaPullRatioPi_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPSigmaPullRatioE_PbPb2018_kOnlySecond; ///
  TGraph *fD0RPMeanPCur_PbPb2018_kOnlySecond[2][24]; ///
  TGraph *fD0RPMeanKCur_PbPb2018_kOnlySecond[2][24]; ///
  TGraph *fD0RPMeanPiCur_PbPb2018_kOnlySecond[2][24]; ///
  TGraph *fD0RPMeanECur_PbPb2018_kOnlySecond[2][24]; ///
  TGraph *fPt1ResPCur_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResKCur_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResPiCur_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResECur_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResPUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResKUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResPiUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResEUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResPUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResKUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResPiUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResEUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPMeanPUpg_PbPb2018_kOnlySecond[2][24]; ///
  TGraph *fD0RPMeanKUpg_PbPb2018_kOnlySecond[2][24]; ///
  TGraph *fD0RPMeanPiUpg_PbPb2018_kOnlySecond[2][24]; ///
  TGraph *fD0RPMeanEUpg_PbPb2018_kOnlySecond[2][24]; ///
  TGraph *fPt1ResPUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResKUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResPiUpg_PbPb2018_kOnlySecond; ///
  TGraph *fPt1ResEUpg_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResPCurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResKCurSA_PbPb2018_kOnlySecond; ///
  TGraph *fD0ZResPiCurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResECurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResPCurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResKCurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResPiCurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResECurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResPCurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResKCurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResPiCurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResECurSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResPUpgSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResKUpgSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResPiUpgSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0ZResEUpgSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResPUpgSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResKUpgSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResPiUpgSA_PbPb2018_kOnlySecond; /// 
  TGraph *fD0RPResEUpgSA_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResPUpgSA_PbPb2018_kOnlySecond; ///
  TGraph *fPt1ResKUpgSA_PbPb2018_kOnlySecond; /// 
  TGraph *fPt1ResPiUpgSA_PbPb2018_kOnlySecond; ///
  TGraph *fPt1ResEUpgSA_PbPb2018_kOnlySecond; /// 

  Bool_t fRunInVertexing; /// flag to run hybrid task before the vertexingHF task or in standard mode
  Bool_t fImproveTracks; /// this is always kTRUE. kFALSE only if re-running on already improved AODs
  Bool_t fUpdateSecVertCovMat; /// flag to swicth on/off the modification of the sec vert cov matrix
  Bool_t fUpdateSTCovMatrix; /// flag to switch on/off the update of the single track covariance matrix
  Bool_t fUpdatePulls; /// flag to switch on/off the correction of the pulls
  Bool_t fMimicData;
  Bool_t fIsAOD;          /// flag to run on AOD
  Bool_t fSmearOnlySignal;  /// flag to control whether to smear only injected signal 
  TClonesArray *fMCs;      // pointer to AOD MC info
  TList   *fDebugOutput; //!<! collection of debug output
  TNtuple *fDebugNtuple; //!<! debug send on output slot 1
  Float_t *fDebugVars;   //!<! variables to store as degug info 
  Int_t   fNDebug;       /// Max number of debug entries into Ntuple
  TString fImproverSuffix; /// suffix for path of correction file on grid 
  TString fOverridePeriodName; /// custom name for overriding auto period name
  Bool_t fFilesOpen; 	 /// check to ensure improver files from CVMFS are accessed only once
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSEImproveITSCVMFS,13);
  /// \endcond
};

#endif

