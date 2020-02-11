/* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALI_ANALYSIS_TASK_SE_IMPROVE_ITS3_H
#define ALI_ANALYSIS_TASK_SE_IMPROVE_ITS3_H

#include "AliAnalysisTaskSE.h"

class TGraph;
class TList;
class AliAODTrack;
class TClonesArray;
class TObjArray;
class AliESDVertex;
class AliVVertex;

class AliAnalysisTaskSEImproveITS3:public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSEImproveITS3();
  AliAnalysisTaskSEImproveITS3(const char *name,
                              const char *resfileCurURI,
                              const char *resfileUpgURI,
                              Bool_t isRunInVertexing,
                              Bool_t isImproveDeuteron,
                              Int_t ndebug);

  virtual ~AliAnalysisTaskSEImproveITS3();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  //  virtual void Init();
  //  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  //  virtual void Terminate(Option_t *option);
  void SetImproveTracks(Bool_t flag=kTRUE) { fImproveTracks=flag; return; }
  void SetUpdateSTCovMatrix(Bool_t opt=kTRUE){fUpdateSTCovMatrix=opt;}
  void SetUpdateSecVertCovMat(Bool_t flag=kTRUE) { fUpdateSecVertCovMat=flag; return; }
  void SetOnlyProcessFilledCand(Bool_t flag=kTRUE){ fOnlyProcessFilledCand = flag;}

private:
  AliAnalysisTaskSEImproveITS3(const AliAnalysisTaskSEImproveITS3&);
  AliAnalysisTaskSEImproveITS3& operator=(const AliAnalysisTaskSEImproveITS3&); 

  // Helper functions
  Double_t EvalGraph(Double_t x,const TGraph *graph/*,const TGraph *graphSA=0*/) const; 
  void SmearTrack(AliAODTrack *track,const TClonesArray *mcs);
  AliESDVertex* RecalculateVertex(const AliVVertex *old,TObjArray *tracks,Double_t bField);

  TGraph *fD0ZResPCur  ; // old pt dep. d0 res. in z for protons
  TGraph *fD0ZResKCur  ; // old pt dep. d0 res. in z for kaons
  TGraph *fD0ZResPiCur ; // old pt dep. d0 res. in z for pions
  TGraph *fD0ZResDCur  ; // old pt dep. d0 res. in z for deuterons
  TGraph *fD0RPResPCur ; // old pt dep. d0 res. in rphi for protons
  TGraph *fD0RPResKCur ; // old pt dep. d0 res. in rphi for kaons
  TGraph *fD0RPResPiCur; // old pt dep. d0 res. in rphi for pions
  TGraph *fD0RPResDCur ; // old pt dep. d0 res. in rphi for deuterons
  TGraph *fPt1ResPCur  ; // old pt dep. 1/pt res. for protons
  TGraph *fPt1ResKCur  ; // old pt dep. 1/pt res. for kaons
  TGraph *fPt1ResPiCur ; // old pt dep. 1/pt res. for pions
  TGraph *fPt1ResDCur  ; // old pt dep. 1/pt res. for deuterons
  TGraph *fD0ZResPUpg  ; // new pt dep. d0 res. in z for protons
  TGraph *fD0ZResKUpg  ; // new pt dep. d0 res. in z for kaons
  TGraph *fD0ZResPiUpg ; // new pt dep. d0 res. in z for pions
  TGraph *fD0ZResDUpg  ; // new pt dep. d0 res. in z for deuterons
  TGraph *fD0RPResPUpg ; // new pt dep. d0 res. in rphi for protons
  TGraph *fD0RPResKUpg ; // new pt dep. d0 res. in rphi for kaons
  TGraph *fD0RPResPiUpg; // new pt dep. d0 res. in rphi for pions
  TGraph *fD0RPResDUpg ; // new pt dep. d0 res. in rphi for deuterons
  TGraph *fPt1ResPUpg  ; // new pt dep. 1/pt res. for protons
  TGraph *fPt1ResKUpg  ; // new pt dep. 1/pt res. for kaons
  TGraph *fPt1ResPiUpg ; // new pt dep. 1/pt res. for pions
  TGraph *fPt1ResDUpg  ; // new pt dep. 1/pt res. for deuterons
/*  TGraph *fD0ZResPCurSA  ; // old standalone pt dep. d0 res. in z for protons
  TGraph *fD0ZResKCurSA  ; // old standalone pt dep. d0 res. in z for kaons
  TGraph *fD0ZResPiCurSA ; // old standalone pt dep. d0 res. in z for pions
  TGraph *fD0RPResPCurSA ; // old standalone pt dep. d0 res. in rphi for protons
  TGraph *fD0RPResKCurSA ; // old standalone pt dep. d0 res. in rphi for kaons
  TGraph *fD0RPResPiCurSA; // old standalone pt dep. d0 res. in rphi for pions
  TGraph *fPt1ResPCurSA  ; // old standalone pt dep. 1/pt res. for protons
  TGraph *fPt1ResKCurSA  ; // old standalone pt dep. 1/pt res. for kaons
  TGraph *fPt1ResPiCurSA ; // old standalone pt dep. 1/pt res. for pions
  TGraph *fD0ZResPUpgSA  ; // new standalone pt dep. d0 res. in z for protons
  TGraph *fD0ZResKUpgSA  ; // new standalone pt dep. d0 res. in z for kaons
  TGraph *fD0ZResPiUpgSA ; // new standalone pt dep. d0 res. in z for pions
  TGraph *fD0RPResPUpgSA ; // new standalone pt dep. d0 res. in rphi for protons
  TGraph *fD0RPResKUpgSA ; // new standalone pt dep. d0 res. in rphi for kaons
  TGraph *fD0RPResPiUpgSA; // new standalone pt dep. d0 res. in rphi for pions
  TGraph *fPt1ResPUpgSA  ; // new standalone pt dep. 1/pt res. for protons
  TGraph *fPt1ResKUpgSA  ; // new standalone pt dep. 1/pt res. for kaons
  TGraph *fPt1ResPiUpgSA ; // new standalone pt dep. 1/pt res. for pions
*/
  Bool_t fRunInVertexing; // flag to run hybrid task before the vertexingHF task or in standard mode
  Bool_t fImproveDeuteron; // flag to switch deuteron smearing on/off
  Bool_t fImproveTracks; // this is always kTRUE. kFALSE only if re-running on already improved AODs
  Bool_t fUpdateSTCovMatrix; /// flag to switch on/off the update of the single track covariance matrix
  Bool_t fUpdateSecVertCovMat; /// flag to swicth on/off the modification of the sec vert cov matrix
  TList   *fDebugOutput; //! collection of debug output
  TNtuple *fDebugNtuple; //! debug send on output slot 1
  Float_t *fDebugVars;   //! variables to store as degug info 
  Int_t   fNDebug;       // Max number of debug entries into Ntuple
  Bool_t fOnlyProcessFilledCand; ///Flag to only process already filled candidates and skip others

  ClassDef(AliAnalysisTaskSEImproveITS3,7);
};

#endif

