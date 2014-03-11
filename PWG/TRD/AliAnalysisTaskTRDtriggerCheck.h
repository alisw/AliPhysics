#ifndef ALIANALYSISTASKTRDTRIGGERCHECK_H
#define ALIANALYSISTASKTRDTRIGGERCHECK_H

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "AliLog.h"

#include "AliAnalysisTaskSE.h"

#define ID(x) x, #x
#define LAB(x) x+1, #x

class TList;

class AliAnalysisTaskTRDtriggerCheck :
  public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskTRDtriggerCheck(const char *name = "trd_trg_check");
  ~AliAnalysisTaskTRDtriggerCheck();

  // analysis operations
  virtual void   UserCreateOutputObjects();
  virtual Bool_t Notify();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t *option);

  // task configuration

  // histograms
  enum Hist_t {
      kHistStat = 0,
      kHistTrgStat,
      kHistTrgStatSec,
      kHistTrgStatStackCond,
      kHistTrgStatStackCondNotFired,
      kHistLast
  };

  // statistics
  enum Stat_t {
      kStatSeen = 1,
      kStatLast
  };

  enum Trigger_t {
    kTrgHCO = 0,
    kTrgHJT,
    kTrgHSE,
    kTrgHQU,
    kTrgHEE,
    kTrgLast
  };

  enum Comb_t {
    // ok
    kNothing = 1,
    kFiredCondition,
    kTriggeredFiredCondition,
    // not ok
    kTriggered,
    kFired,
    kCondition,
    kTriggeredFired,
    kTriggeredCondition,
    kLastComb
  };

protected:
  // output objects
  TList *fOutputList;		// list of output objects

  // histogram management
  TH1  *fHist[kHistLast];	//! pointers to histogram
  const char *fShortTaskId;	//! short identifier for the task

  TH1*&  GetHistogram(Hist_t hist, const Int_t idx = 0) { return fHist[hist + idx]; }

  TH1*   AddHistogram(Hist_t hist, const char *hid, TString title,
                      Int_t xbins, Float_t xmin, Float_t xmax, Int_t binType = 1);
  TH2*   AddHistogram(Hist_t hist, const char *hid, TString title,
                      Int_t xbins, Float_t xmin, Float_t xmax,
                      Int_t ybins, Float_t ymin, Float_t ymax, Int_t binType = 1);
  TH3*   AddHistogram(Hist_t hist, const char *hid, TString title,
                      Int_t xbins, Float_t xmin, Float_t xmax,
                      Int_t ybins, Float_t ymin, Float_t ymax,
                      Int_t zbins, Float_t zmin, Float_t zmax, Int_t binType = 1);

  void    FillH1(Hist_t hist, Float_t x, Float_t weight = 1., Int_t idx = 0)
  { GetHistogram(hist, idx)->Fill(x, weight); }
  void    FillH2(Hist_t hist, Float_t x, Float_t y, Float_t weight = 1., Int_t idx = 0)
  { ((TH2*) GetHistogram(hist, idx))->Fill(x, y, weight); }
  void    FillH3(Hist_t hist, Float_t x, Float_t y, Float_t z, Float_t weight = 1., Int_t idx = 0)
  { ((TH3*) GetHistogram(hist, idx))->Fill(x, y, z, weight); }

  // internal use

  // task configuration

private:
  // not implemented
  AliAnalysisTaskTRDtriggerCheck(const AliAnalysisTaskTRDtriggerCheck &rhs);
  AliAnalysisTaskTRDtriggerCheck& operator=(const AliAnalysisTaskTRDtriggerCheck &rhs);

  ClassDef(AliAnalysisTaskTRDtriggerCheck, 1);
};

#endif
