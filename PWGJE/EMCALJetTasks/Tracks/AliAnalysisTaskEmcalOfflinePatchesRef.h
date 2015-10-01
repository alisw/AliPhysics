#ifndef AliAnalysisTaskEmcalOfflinePatchesRef_H
#define AliAnalysisTaskEmcalOfflinePatchesRef_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliAnalysisTaskSE.h>

class AliAnalysisUtils;
class AliEMCALGeometry;
class AliEmcalTriggerPatchInfo;

namespace EMCalTriggerPtAnalysis {

class AliEMCalHistoContainer;

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
  AliEMCalHistoContainer                  *fHistos;

private:
  AliAnalysisTaskEmcalOfflinePatchesRef(const AliAnalysisTaskEmcalOfflinePatchesRef &);
  AliAnalysisTaskEmcalOfflinePatchesRef &operator=(const AliAnalysisTaskEmcalOfflinePatchesRef &);

  void FillTriggerPatchHistos(const char *patchtype, const AliEmcalTriggerPatchInfo * const recpatch, Int_t supermodule, Int_t sector);

  ClassDef(AliAnalysisTaskEmcalOfflinePatchesRef, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* AliAnalysisTaskEmcalOfflinePatchesRef_H */
