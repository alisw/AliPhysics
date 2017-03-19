#ifndef __AliAnalysisTaskEventCutsValidation__
#define __AliAnalysisTaskEventCutsValidation__

#include <TF1.h>

#include "AliAnalysisTaskSE.h"
#include "AliNuclexEventCuts.h"
#include "AliESDtrackCuts.h"

class TList;
class TH2F;

class AliAnalysisTaskEventCutsValidation : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEventCutsValidation(bool storeCuts = false, TString taskname = "EventCutsValidation");
  virtual ~AliAnalysisTaskEventCutsValidation();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  AliNuclexEventCuts  fEventCut;    ///<  Event cuts
  AliESDtrackCuts     fFB32cuts;    ///<  Track cuts corresponding to filter bit 32
  AliESDtrackCuts     fTPConlyCuts; ///< Cuts to select TPC only tracks

private:
  AliAnalysisTaskEventCutsValidation (const AliAnalysisTaskEventCutsValidation &source);
  AliAnalysisTaskEventCutsValidation &operator=(const AliAnalysisTaskEventCutsValidation &source);

  TList *fList;            //!<! Output list

  TH2F  *fTOFvsFB32[2];    //!<!
  TH2F  *fTPCvsAll[2];     //!<!
  TH2F  *fMultvsV0M[2];    //!<!

  TF1* fMultTOFLowCut;     //!<!
  TF1* fMultTOFHighCut;    //!<!
  TF1* fMultCentLowCut;    //!<!

  bool                fStoreCuts; ///<  If true the EventCuts object is put in the output
  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskEventCutsValidation, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskEventCutsValidation__) */
