#ifndef ALIANALYSISTASKEGAMONITOR_H
#define ALIANALYSISTASKEGAMONITOR_H

class THistManager;
class AliOADBContainer;

#include "AliAnalysisTaskEmcal.h"
#include <TString.h>
#include <vector>

namespace PWGJE {

namespace EMCALJetTasks {

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

  void AddMaskedFastor(int fastorAbsID) { fMaskedFastors.push_back(fastorAbsID); }
  void SetMaskedFastorOADB(TString oadbname) { fNameMaskedFastorOADB = oadbname; }
  bool IsPatchRejected(int col, int row);

protected:
  virtual void UserCreateOutputObjects();
  virtual bool IsEventSelected();
  virtual bool Run();

  virtual void ExecOnce();
  virtual void RunChanged(Int_t newrun);

private:
  THistManager                            *fHistos;                 //!<!  Histogram manager

  Bool_t                                   fUseRecalcPatches;       ///< Defined whether to use recalc patches
  Double_t                                 fRecalcLow;              ///< Low threshold for recalc gamma trigger
  Double_t                                 fRecalcHigh;             ///< High threshold for recalc gamma trigger

  TString                                  fNameMaskedFastorOADB;   ///< Name of the OADB container with the masked fastor information
  AliOADBContainer                        *fMaskedFastorOADB;       ///< OADB container with masked fastor information
  std::vector<int>                         fMaskedFastors;          ///< List of masked fastors

  ClassDef(AliAnalysisTaskEGAMonitor, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEGAMONITOR_H */
