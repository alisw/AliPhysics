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
////#include "framework/AliQnCorrectionsCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysis;
//class AliReducedEvent;
class AliQnCorrectionsCuts;
class AliQnCorrectionsFillEvent;
class AliQnCorrectionsManager;
class AliQnCorrectionsHistos;
//class QnCorrectionsReducedVarManager;

//_________________________________________________________
//class AliAnalysisTaskFlowVectorCorrections : public TObject {
class AliAnalysisTaskFlowVectorCorrections : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskFlowVectorCorrections();
  AliAnalysisTaskFlowVectorCorrections(const char *name);
  virtual ~AliAnalysisTaskFlowVectorCorrections(){}


  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();


  void SetRunByRunCalibration(Bool_t cal) {fCalibrateByRun = cal;}
  void SetCalibrationHistograms(TList* input) {fListInputHistogramsQnCorrections = input;}

  //void InitializeCalibrationHistograms(TString label);
  void SetEventPlaneManager(AliQnCorrectionsManager* EPmanager)  {fEventPlaneManager = EPmanager;}
  void SetVarManager(AliQnCorrectionsFillEvent* eventfill)  {fFillEvent = eventfill;}
  //void SetVarManager(QnCorrectionsReducedVarManager* eventfill)  {fFillEvent = eventfill;}
  void SetEventCuts(AliQnCorrectionsCuts* cuts)  {fEventCuts = cuts;}
  void SetTrigger(UInt_t triggerbit) {fTriggerMask=triggerbit;}
  void AddHistogramClass(TString hist) {fQAhistograms+=hist+";";}
  void DefineInOutput();
  void FillExchangeContainerWithQvectors(Bool_t b=kTRUE) {fProvideQnVectorsList=b;}
  void SetRunListPath(TString path) {fRunListPath=path;}
  void SetCalibrationFilePath(TString path) {fCalibrationFilePath=path;}

  AliQnCorrectionsManager* EventPlaneManager() {return fEventPlaneManager;}
  AliQnCorrectionsHistos* GetHistograms() {return fEventPlaneHistos;}
  AliQnCorrectionsFillEvent* GetFillEvent() {return fFillEvent;}
  AliQnCorrectionsCuts* EventCuts()  const {return fEventCuts;}
  Int_t OutputSlotEventQA()        const {return fOutputSlotEventQA;}
  Int_t OutputSlotHistQA()        const {return fOutputSlotHistQA;}
  Int_t OutputSlotHistQn()        const {return fOutputSlotHistQn;}
  Int_t OutputSlotGetListQnVectors() const {return fOutputSlotQnVectorsList;}
  Int_t OutputSlotTree()          const {return fOutputSlotTree;}
  Bool_t IsEventSelected(Float_t* values);
  Bool_t IsFillExchangeContainerWithQvectors() const  {return fProvideQnVectorsList;}
  Bool_t IsFillEventQA() const  {return fFillEventQA;}

 private:
  Bool_t fRunLightWeight;
  Bool_t fCalibrateByRun;
  Bool_t fUseFriendEvent;
  UInt_t fTriggerMask;
  Bool_t fInitialized;
  TList* fListInputHistogramsQnCorrections;          //! List of input histograms for corrections
  TList* fEventQAList;
  AliQnCorrectionsManager * fEventPlaneManager;
  AliQnCorrectionsCuts * fEventCuts;
  //QnCorrectionsReducedVarManager* fFillEvent;
  AliQnCorrectionsFillEvent* fFillEvent;
  AliQnCorrectionsHistos* fEventPlaneHistos;
  TString fLabel;
  TString fQAhistograms;
  Bool_t fFillEventQA;
  Bool_t fProvideQnVectorsList;
  Int_t fOutputSlotEventQA;
  Int_t fOutputSlotHistQA;
  Int_t fOutputSlotHistQn;
  Int_t fOutputSlotQnVectorsList;
  Int_t fOutputSlotTree;
  TString fRunListPath;
  TString fCalibrationFilePath;
  
  AliAnalysisTaskFlowVectorCorrections(const AliAnalysisTaskFlowVectorCorrections &c);
  AliAnalysisTaskFlowVectorCorrections& operator= (const AliAnalysisTaskFlowVectorCorrections &c);

  ClassDef(AliAnalysisTaskFlowVectorCorrections, 2);
};

#endif


