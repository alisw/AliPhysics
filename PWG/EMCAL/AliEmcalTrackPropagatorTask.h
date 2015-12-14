#ifndef ALIEMCALTRACKPROPAGATORTASK_H
#define ALIEMCALTRACKPROPAGATORTASK_H

#include "AliAnalysisTaskEmcal.h"

class AliEmcalTrackPropagatorTask : public AliAnalysisTaskEmcal {
 public:
  AliEmcalTrackPropagatorTask();
  AliEmcalTrackPropagatorTask(const char *name);
  virtual ~AliEmcalTrackPropagatorTask();

  void               SetDist(Double_t d)               { fDist           = d; }
  void               SetOnlyIfNotSet(Bool_t b)         { fOnlyIfNotSet   = b; }
  void               SetOnlyIfEmcal(Bool_t b)          { fOnlyIfEmcal    = b; }

 protected:
  void               ExecOnce();
  Bool_t             Run();
   
  Double_t           fDist;              // distance to surface (440cm default)
  Bool_t             fOnlyIfNotSet;      // propagate only if needed
  Bool_t             fOnlyIfEmcal;       // propagate only if it is in the EMCal acceptance

 private:
  AliEmcalTrackPropagatorTask(const AliEmcalTrackPropagatorTask&);            // not implemented
  AliEmcalTrackPropagatorTask &operator=(const AliEmcalTrackPropagatorTask&); // not implemented

  ClassDef(AliEmcalTrackPropagatorTask, 3); // Class to propagate and store track parameters at EMCAL surface
};
#endif
