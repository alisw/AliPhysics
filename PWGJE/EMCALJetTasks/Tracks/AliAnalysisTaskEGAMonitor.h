#ifndef ALIANALYSISTASKEGAMONITOR_H
#define ALIANALYSISTASKEGAMONITOR_H

class THistManager;

#include "AliAnalysisTaskEmcal.h"

namespace EMCalTriggerPtAnalysis {

/**
 * @class AliAnalysisTaskEGAMonitor
 * @brief Simplistic class, checks distribution of the online EGA trigger patches
 */
class AliAnalysisTaskEGAMonitor : public AliAnalysisTaskEmcal {
public:
  AliAnalysisTaskEGAMonitor();
  AliAnalysisTaskEGAMonitor(const char *name);
  virtual ~AliAnalysisTaskEGAMonitor();

  void SetUseRecalcPatches(Double_t recalcLow, Double_t recalcHigh) {
    fUseRecalcPatches = kTRUE;
    fRecalcLow = recalcLow;
    fRecalcHigh = recalcHigh;
  }

protected:
  virtual void UserCreateOutputObjects();
  virtual bool IsEventSelected();
  virtual bool Run();

private:
  THistManager                            *fHistos;           //!<!  Histogram manager

  Bool_t                                   fUseRecalcPatches; ///< Defined whether to use recalc patches
  Double_t                                 fRecalcLow;        ///< Low threshold for recalc gamma trigger
  Double_t                                 fRecalcHigh;       ///< High threshold for recalc gamma trigger

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEGAMonitor, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEGAMONITOR_H */
