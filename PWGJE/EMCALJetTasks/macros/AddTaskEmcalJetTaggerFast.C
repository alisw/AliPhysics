#if !defined(__CINT__) && !defined(__CLING__)
#include "AliEmcalJetTaggerTaskFast.h"
#endif

PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast *AddTaskEmcalJetTaggerFast(const char * njetsBase,
    const char * njetsTag,
    Double_t     R,
    const char * nrhoBase,
    const char * nrhoTag,
    const char * ntracks,
    const char * nclusters,
    const char * type,
    const char * CentEst,
    Int_t        pSel,
    const char * trigClass) {
  return PWGJE::EMCALJetTasks::AliEmcalJetTaggerTaskFast::AddTaskJetTaggerFast(njetsBase, njetsTag, R, nrhoBase, nrhoTag, ntracks, nclusters, type, CentEst, pSel, trigClass);
}
