#ifndef ALIANALYSISTASKTRIGGERREJECTION_H
#define ALIANALYSISTASKTRIGGERREJECTION_H

class TH1;
class TH2;
class TH3;
class TH3F;
class TProfile;
class THnSparse;
class TClonesArray;
class TArrayI;

#include <TRef.h>
#include <TBits.h>
#include <TMath.h>

#include <AliVEvent.h>

#include "AliAnalysisTaskEmcalJet.h"
#include "AliEmcalTriggerPatchInfo.h"

class AliAnalysisTaskTriggerRejection : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskTriggerRejection();
  AliAnalysisTaskTriggerRejection(const char *name);
  virtual ~AliAnalysisTaskTriggerRejection();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  //Setters
  void SetContainerFull(Int_t c)            { fContainerFull      = c;}
  void SetContainerCharged(Int_t c)         { fContainerCharged   = c;}

 protected:
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        ExtractMainPatch();

 private:
  Int_t              fContainerFull;         // number of container with full jets DET
  Int_t              fContainerCharged;      // number of container with charged jets DET
  AliEmcalTriggerPatchInfo *fMaxPatch;       // main patch
  THnSparse          *fhnTriggerInfo;        //! correlation between jets, patch energy and event observables

  AliAnalysisTaskTriggerRejection(const AliAnalysisTaskTriggerRejection&);            // not implemented
  AliAnalysisTaskTriggerRejection &operator=(const AliAnalysisTaskTriggerRejection&); // not implemented

  ClassDef(AliAnalysisTaskTriggerRejection, 1)
};
#endif


