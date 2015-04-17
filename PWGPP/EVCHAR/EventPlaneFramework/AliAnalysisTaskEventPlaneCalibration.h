/*
***********************************************************
  event plane corrections framework
  contact: jaap onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  2014/12/10
  *********************************************************
*/

//#include "AliSysInfo.h"

#ifndef ALIANALYSISTASKEVENTPLANECALIBRATION_H
#define ALIANALYSISTASKEVENTPLANECALIBRATION_H

#include "TFile.h"
#include "TTree.h"
//#include "framework/AliEventPlaneCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysis;
class AliReducedEvent;
class AliEventPlaneCuts;
class AliEventPlaneBinning;
class AliEventPlaneQvector;
class AliEventPlaneVarManager;
class AliEventPlaneConfiguration;
class AliEventPlaneDetector;
class AliEventPlaneManager;
class AliEventPlaneHelper;

//_________________________________________________________
//class AliAnalysisTaskEventPlaneCalibration : public TObject {
class AliAnalysisTaskEventPlaneCalibration : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskEventPlaneCalibration();
  AliAnalysisTaskEventPlaneCalibration(const char *name);
  virtual ~AliAnalysisTaskEventPlaneCalibration(){}


  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();


  void SetRunByRunCalibration(Bool_t cal) {fCalibrateByRun = cal;}
  //Bool_t IsEventSelected(AliReducedEvent* event);
  Bool_t IsEventSelected(Float_t* values);
  //Bool_t MultiplicitySelected(AliReducedEvent* event);
  //Bool_t TriggerSelected(AliReducedEvent* event);
  AliEventPlaneManager* EventPlaneManager() {return fEventPlaneManager;}
  Bool_t IsAOD() const {return fIsAOD;}
  Bool_t IsESD() const {return fIsESD;}
  Bool_t IsReduced() const {return fIsReduced;}

  void DefineHistograms(const Char_t* histClasses, Int_t runNumber, AliEventPlaneManager* EPmanager);
  void InitializeCalibrationHistograms(Int_t currentRunNumber, AliEventPlaneManager* epconf);
  void InitEqualizationHistograms(Int_t runNo, AliEventPlaneConfiguration* epconf);
  void InitQvecCalibrationHistograms(Int_t runNo, AliEventPlaneConfiguration* epconf);
  void SetEventPlaneManager(AliEventPlaneManager* EPmanager)  {fEventPlaneManager = EPmanager;}
  void SetEventCuts(AliEventPlaneCuts* cuts)  {fEventCuts = cuts;}
  void SetInputAOD() {fIsAOD=kTRUE;}
  void SetInputESD() {fIsESD=kTRUE;}
  void SetInputReduced() {fIsReduced=kTRUE;} // input AliReducedEvent
  void SetQvectorFile(TFile* file) {fQvectorFile=file;} 
  TFile* QvectorFile() {return fQvectorFile;} 

 private:
  Bool_t fRunLightWeight;
  Bool_t fCalibrateByRun;
  Bool_t fUseFriendEvent;
  Bool_t fFillTPC;
  Bool_t fFillVZERO;
  Bool_t fFillTZERO;
  Bool_t fFillZDC;
  Bool_t fFillFMD;
  Bool_t fChannelEqualization;
  Bool_t fRecenterQvec;
  Bool_t fRotateQvec;
  Bool_t fTwistQvec;
  Bool_t fScaleQvec;
  Bool_t fInitialized;
  TTree* fTree;
  TFile* fHistosFile;
  TFile* fQvectorFile;
  TList fListHistos;                 //! List of histogram managers in the eventplane framework classes
  TList fListHistosQA;                 //! List of histogram managers in the eventplane framework classes
  AliEventPlaneManager * fEventPlaneManager;
  AliEventPlaneCuts * fEventCuts;
  Int_t fFriendEventNo;
  Int_t fOffset;
  Int_t fNfiles;
  Bool_t fIsAOD;
  Bool_t fIsESD;
  Bool_t fIsReduced;
  
  AliAnalysisTaskEventPlaneCalibration(const AliAnalysisTaskEventPlaneCalibration &c);
  AliAnalysisTaskEventPlaneCalibration& operator= (const AliAnalysisTaskEventPlaneCalibration &c);

  ClassDef(AliAnalysisTaskEventPlaneCalibration, 1);
};

#endif


