#ifndef ALIANALYSISTASKEGAMONITOR_H
#define ALIANALYSISTASKEGAMONITOR_H

class AliAnalysisUtils;
class AliEMCALGeometry;
class THistManager;

#include "AliAnalysisTaskSE.h"

namespace EMCalTriggerPtAnalysis {

/**
 * @class AliAnalysisTaskEGAMonitor
 * @brief Simplistic class, checks distribution of the online EGA trigger patches
 */
class AliAnalysisTaskEGAMonitor : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEGAMonitor();
  AliAnalysisTaskEGAMonitor(const char *name);
  virtual ~AliAnalysisTaskEGAMonitor();

protected:
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);

  virtual void ExecOnce();

private:
  THistManager                            *fHistos;           //!<!  Histogram manager
  AliAnalysisUtils                        *fAnalysisUtils;    //!<!  Analysis utils
  AliEMCALGeometry                        *fGeom;             //!<!  EMCAL geometry (for mapping)
  Bool_t                                  fLocalInitialized;  //!<!  Local initialization

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEGAMonitor, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEGAMONITOR_H */
