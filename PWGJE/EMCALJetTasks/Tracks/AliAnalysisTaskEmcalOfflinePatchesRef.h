#ifndef AliAnalysisTaskEmcalOfflinePatchesRef_H
#define AliAnalysisTaskEmcalOfflinePatchesRef_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliAnalysisTaskSE.h>

class AliAnalysisUtils;
class AliEMCALGeometry;
class AliEMCALTriggerPatchInfo;
class THistManager;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalOfflinePatchesRef: public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEmcalOfflinePatchesRef();
  AliAnalysisTaskEmcalOfflinePatchesRef(const char *name);
  virtual ~AliAnalysisTaskEmcalOfflinePatchesRef();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);

protected:
  AliAnalysisUtils                        *fAnalysisUtil;
  AliEMCALGeometry                        *fGeometry;
  THistManager                            *fHistos;

private:
  AliAnalysisTaskEmcalOfflinePatchesRef(const AliAnalysisTaskEmcalOfflinePatchesRef &);
  AliAnalysisTaskEmcalOfflinePatchesRef &operator=(const AliAnalysisTaskEmcalOfflinePatchesRef &);

  void FillTriggerPatchHistos(const char *patchtype, const AliEMCALTriggerPatchInfo * const recpatch, Int_t supermodule, Int_t sector, Bool_t evsel);

  ClassDef(AliAnalysisTaskEmcalOfflinePatchesRef, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* AliAnalysisTaskEmcalOfflinePatchesRef_H */
