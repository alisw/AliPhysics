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
#include "AliEMCALTriggerPatchInfo.h"

namespace JETriggerRejectionAna {
  class AliAnalysisTaskTriggerRejection : public AliAnalysisTaskEmcalJet {
  public:
    enum MainPatchType {
      kManual = 0,    //just select highest energy patch in array
      kEmcalJet = 1   //use functionality of AliAnalysisTaskEmcal
    };
    AliAnalysisTaskTriggerRejection();
    AliAnalysisTaskTriggerRejection(const char *name);
    virtual ~AliAnalysisTaskTriggerRejection();

    void                        UserCreateOutputObjects();
    void                        Terminate(Option_t *option);
    
    //Setters
    void SetContainerFull(Int_t c)            { fContainerFull      = c;}
    void SetContainerCharged(Int_t c)         { fContainerCharged   = c;}
    void SetMainPatchType(MainPatchType t)    { fMainPatchType      = t;}
    void SetMainTriggerTypeCat(TriggerCategory cat, Bool_t b) {fMainTrigCat = cat; fMainTrigSimple = b;}

  protected:
    Bool_t                      FillHistograms()   ;
    Bool_t                      Run()              ;
    void                        ExtractMainPatch();

  private:
    Int_t              fContainerFull;         // number of container with full jets DET
    Int_t              fContainerCharged;      // number of container with charged jets DET
    AliEMCALTriggerPatchInfo *fMaxPatch;       // main patch
    THnSparse         *fhnTriggerInfo;         //! correlation between jets, patch energy and event observables
    MainPatchType      fMainPatchType;         // method to select main patch
    TriggerCategory    fMainTrigCat;           // trigger category for main trigger from AliAnalysisTaskEmcal::GetMainTriggerPatch
    Bool_t             fMainTrigSimple;        // use offline trigger instead of online from AliAnalysisTaskEmcal::GetMainTriggerPatch

    AliAnalysisTaskTriggerRejection(const AliAnalysisTaskTriggerRejection&);            // not implemented
    AliAnalysisTaskTriggerRejection &operator=(const AliAnalysisTaskTriggerRejection&); // not implemented
    
    ClassDef(AliAnalysisTaskTriggerRejection, 2)
      };
}
#endif


