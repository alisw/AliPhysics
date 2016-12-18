#ifndef ALIANALYSISTASKFLOWVECTORCORRECTION_H
#define ALIANALYSISTASKFLOWVECTORCORRECTION_H

/***************************************************************************
 * Package:       FlowVectorCorrections ALICE glue                         *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2014-2016                                                *
 ***************************************************************************/

#include <TObject.h>
#include "Rtypes.h"

#include "TFile.h"
#include "TTree.h"

#include "AliAnalysisTaskSE.h"
#include "AliQnCorrectionsFillEventTask.h"

class AliQnCorrectionsManager;
class AliQnCorrectionsCutsSet;
class AliQnCorrectionsHistos;

class AliAnalysisTaskFlowVectorCorrections : public AliQnCorrectionsFillEventTask {

public:
  /// \enum CalibrationFileSource
  /// \brief The supported sources for the calibration file
  enum CalibrationFileSource {
    CALIBSRC_local,    ///< the calibration file will be taken locally when the task object is created
    CALIBSRC_aliensingle,    ///< the calibration file, which contains correction parameters for every run, will be taken from alien on each execution node
    CALIBSRC_alienmultiple,  ///< run calibration file, which contains correction parameters only for the intended run, will be taken from alien on each execution node
    CALIBSRC_OADBsingle,    ///< the calibration file, which contains correction parameters for every run, will be taken from OADB on each execution node
    CALIBSRC_OADBmultiple,  ///< run calibration file, which contains correction parameters only for the intended run, will be taken from OADB on each execution node
  };


  AliAnalysisTaskFlowVectorCorrections();
  AliAnalysisTaskFlowVectorCorrections(const char *name);
  virtual ~AliAnalysisTaskFlowVectorCorrections(){}


  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  virtual void NotifyRun();


  void SetRunByRunCalibration(Bool_t enable) { fCalibrateByRun = enable; }
  void SetAliQnCorrectionsManager(AliQnCorrectionsManager* QnManager)  {fAliQnCorrectionsManager = QnManager;}
  void SetEventCuts(AliQnCorrectionsCutsSet *cuts)  {fEventCuts = cuts;}
  void SetFillExchangeContainerWithQvectors(Bool_t enable = kTRUE) { fProvideQnVectorsList = enable; }
  void SetFillEventQA(Bool_t enable = kTRUE) { fFillEventQA = enable; }
  void SetTrigger(UInt_t triggerbit) {fTriggerMask=triggerbit;}
  void AddHistogramClass(TString hist) {fQAhistograms+=hist+";";}
  void SetCalibrationHistogramsFile(CalibrationFileSource source, const char *filename);
  void DefineInOutput();
  void SetRunsLabels(TObjArray *runsList) { fAliQnCorrectionsManager->SetListOfProcessesNames(runsList); }

  AliQnCorrectionsManager *GetAliQnCorrectionsManager() {return fAliQnCorrectionsManager;}
  AliQnCorrectionsHistos* GetEventHistograms() {return fEventHistos;}
  AliQnCorrectionsCutsSet* GetEventCuts()  const {return fEventCuts;}
  Int_t OutputSlotEventQA()        const {return fOutputSlotEventQA;}
  Int_t OutputSlotHistQA()        const {return fOutputSlotHistQA;}
  Int_t OutputSlotHistNveQA()        const {return fOutputSlotHistNveQA;}
  Int_t OutputSlotHistQn()        const {return fOutputSlotHistQn;}
  Int_t OutputSlotGetListQnVectors() const {return fOutputSlotQnVectorsList;}
  Int_t OutputSlotTree()          const {return fOutputSlotTree;}
  Bool_t IsEventSelected(Float_t* values);
  Bool_t GetFillExchangeContainerWithQvectors() const  {return fProvideQnVectorsList;}
  Bool_t GetFillEventQA() const  {return fFillEventQA;}

private:
  Bool_t fCalibrateByRun;
  TString fCalibrationFile;                       ///< the name of the calibration file
  CalibrationFileSource fCalibrationFileSource;   ///< the source of the calibration file
  UInt_t fTriggerMask;
  TList* fEventQAList;
  AliQnCorrectionsCutsSet *fEventCuts;
  TString fLabel;
  TString fQAhistograms;
  Bool_t fFillEventQA;
  Bool_t fProvideQnVectorsList;
  Int_t fOutputSlotEventQA;
  Int_t fOutputSlotHistQA;
  Int_t fOutputSlotHistNveQA;
  Int_t fOutputSlotHistQn;
  Int_t fOutputSlotQnVectorsList;
  Int_t fOutputSlotTree;

  AliAnalysisTaskFlowVectorCorrections(const AliAnalysisTaskFlowVectorCorrections &c);
  AliAnalysisTaskFlowVectorCorrections& operator= (const AliAnalysisTaskFlowVectorCorrections &c);

  ClassDef(AliAnalysisTaskFlowVectorCorrections, 4);
};

#endif // ALIANALYSISTASKFLOWVECTORCORRECTION_H


