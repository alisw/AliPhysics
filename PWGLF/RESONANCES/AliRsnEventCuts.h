//
// Class AliRsnEventCuts
//
// It works with ESD and AOD events.
//
// authors: P. Ganoti, F. Bellini


#ifndef ALIRSNEVENTCUTS_H
#define ALIRSNEVENTCUTS_H

#include "AliRsnCut.h"
#include "AliVEvent.h"
#include "AliEventCuts.h"

class AliVVertex;
class AliEventCuts;

class AliRsnEventCuts : public AliRsnCut {
 public:

  AliRsnEventCuts(const char *name = "EventCuts");
  virtual ~AliRsnEventCuts() {;};
  
  Bool_t         Init(TObject *object);
  Bool_t         IsSelected(TObject *object);
  void           ForceSetupPbPb2018() {fEvCuts.SetManualMode(); fEvCuts.SetupPbPb2018();};
  void           ForceSetupRun2pp() {fEvCuts.SetManualMode(); fEvCuts.OverrideCentralityFramework(1); fEvCuts.SetupRun2pp();};
  void           SetUseMultSelectionEvtSel() {fEvCuts.SetManualMode(); fEvCuts.UseMultSelectionEventSelection();};
  AliEventCuts   fEvCuts; //pointer to the AliAnalysisUtils object

  ClassDef(AliRsnEventCuts, 3)

};

#endif
