/*
  Example of AliEffFDContainer usage.
  The PCC framework (AliMCSpectraWeights) used in AliEffFDContainer is written by Patrick Huhn.
  For postprocessing of output to get efficiencies and feed-down corrections, use the
  postprocessing macro at: https://github.com/vvislavi/Feeddown
  If used, please acknowledge the authors of AliMCSpectraWeights as well as myself
  Author: Vytautas Vislavicius
*/
#ifndef ALIANALYSISTASKEFFANDFDTEST__H
#define ALIANALYSISTASKEFFANDFDTEST__H
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "TString.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"


class TList;
class AliVTrack;
class AliVVertex;
class AliInputEventHandler;
class AliAnalysisUtils;
class AliVParticle;
class AliEffFDContainer;

class AliAnalysisTaskEffFDExample : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEffFDExample();
  AliAnalysisTaskEffFDExample(const char *name, Bool_t IsMC=kTRUE);
  virtual ~AliAnalysisTaskEffFDExample();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  void SetTriggerType(UInt_t newval) {fTriggerType = newval; };
 protected:
  AliEventCuts fEventCuts;
 private:
  AliAnalysisTaskEffFDExample(const AliAnalysisTaskEffFDExample&);
  AliAnalysisTaskEffFDExample& operator=(const AliAnalysisTaskEffFDExample&);
  Bool_t fIsMC;
  AliMCEvent *fMCEvent; //! MC event
  UInt_t fTriggerType;
  AliEffFDContainer *fEfFd;
  ClassDef(AliAnalysisTaskEffFDExample,1);
};

#endif
