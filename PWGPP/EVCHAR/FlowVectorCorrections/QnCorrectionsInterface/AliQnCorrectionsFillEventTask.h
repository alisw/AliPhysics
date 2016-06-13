#ifndef ALIALIQNCORRECTIONS_FILLEVENT_H
#define ALIALIQNCORRECTIONS_FILLEVENT_H

/***************************************************************************
 * Package:       FlowVectorCorrections ALICE glue                         *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2014-2016                                                *
 ***************************************************************************/

#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskSE.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>

#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsVarManagerTask.h"
#include "AliQnCorrectionsHistos.h"

class AliESDtrack;
class AliVParticle;

class AliQnCorrectionsFillEventTask : public AliQnCorrectionsVarManagerTask {
public:

  AliQnCorrectionsFillEventTask();
  AliQnCorrectionsFillEventTask(const char *name);
  ~AliQnCorrectionsFillEventTask();

public:

  virtual void UserExec(Option_t *) = 0;
  virtual void UserCreateOutputObjects() = 0;
  virtual void FinishTaskOutput() = 0;


  void SetUseTPCStandaloneTracks(Bool_t enable = kTRUE) { fUseTPCStandaloneTracks = enable; }

protected:
  /* Fill event data methods */
  void FillEventData();

  void FillDetectors();
  void FillTPC();
  void FillEsdTPC();
  void FillAodTPC();
  void FillVZERO();
  void FillTZERO();
  void FillZDC();
  void FillFMD();
  void FillRawFMD();
  void FillSPDTracklets();

  void FillEventInfo();
  void FillTrackInfo(AliESDtrack* p);
  void FillTrackInfo(AliVParticle* p);

  void SetDetectors();

private:

  AliQnCorrectionsFillEventTask(const AliQnCorrectionsFillEventTask &c);
  AliQnCorrectionsFillEventTask& operator= (const AliQnCorrectionsFillEventTask &c);

protected:
  AliVEvent* fEvent;
  AliQnCorrectionsManager *fAliQnCorrectionsManager;
  AliQnCorrectionsHistos* fEventHistos;
  Float_t *fDataBank;                             //!<! The event variables values data bank. Transient!
private:
  static const Float_t fVZEROSignalThreshold; ///< the VZERO channel signal threshold for building a data vector
  static const Float_t fTZEROSignalThreshold; ///< the TZERO channel signal threshold for building a data vector
  static const Float_t fZDCSignalThreshold; ///< the ZDC channel signal threshold for building a data vector
  static const Float_t fFMDSignalThreshold; ///< the FMD channel signal threshold for building a data vector

  Bool_t fUseTPCStandaloneTracks;
  Bool_t fFillVZERO;
  Bool_t fFillTPC;
  Bool_t fFillZDC;
  Bool_t fFillTZERO;
  Bool_t fFillFMD;
  Bool_t fFillRawFMD;
  Bool_t fFillSPD;
  Bool_t fIsAOD;
  Bool_t fIsESD;

  ClassDef(AliQnCorrectionsFillEventTask, 1);
};

#endif
