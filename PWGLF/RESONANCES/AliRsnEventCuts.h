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
  AliRsnEventCuts(const AliRsnEventCuts &copy);
  AliRsnEventCuts &operator=(const AliRsnEventCuts &copy);
  virtual ~AliRsnEventCuts() {;};
  
  Bool_t         Init(TObject *object);
  Bool_t         IsSelected(TObject *object);
  void           SetSetupPbPb2018(Bool_t doit) {fUsePbPb2018 = doit;};
  Bool_t         IsAcceptedMultSelection();

 private:  
  AliEventCuts        *fEvCuts; //pointer to the AliAnalysisUtils object
  Bool_t              fUsePbPb2018;
  
  ClassDef(AliRsnEventCuts, 1)
    
};

#endif
