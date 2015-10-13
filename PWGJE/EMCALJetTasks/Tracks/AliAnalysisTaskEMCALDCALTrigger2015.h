#ifndef ALIANALYSISTASKEMCALDCALTRIGGER2015_H_
#define ALIANALYSISTASKEMCALDCALTRIGGER2015_H_

#include <TString.h>

#include "AliAnalysisTaskSE.h"

class TArrayD;
class TClonesArray;

class AliEmcalTriggerPatchInfo;
class AliVCluster;

namespace EMCalTriggerPtAnalysis {

class AliEMCalHistoContainer;

class AliAnalysisTaskEMCALDCALTrigger2015 : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEMCALDCALTrigger2015();
  AliAnalysisTaskEMCALDCALTrigger2015(const char *name);
  virtual ~AliAnalysisTaskEMCALDCALTrigger2015() {}

  void UserCreateOutputObjects();
  void UserExec(Option_t * opt);

  void SetClusterContainerName(const char *name) { fClusterContainerName = name; }

protected:
  void ProcessCluster(const TString &triggerclass, const AliVCluster * const clust, bool isCalib);
  void ProcessPatch(const TString &triggerclass, const AliEmcalTriggerPatchInfo * const patch, bool isOnline);
  void CreateEnergyBinning(TArrayD& binning) const;
  void CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const;

  static const TString fkTriggerClasses[6];

  TString                       fClusterContainerName;        //
  AliEMCalHistoContainer        *fHistos;                     //!<!
  AliEMCALGeometry              *fGeometry;                   //!<!
  TClonesArray                  *fClusterContainer;           //!<!
  TClonesArray                  *fPatchContainer;             //!<!

private:
  AliAnalysisTaskEMCALDCALTrigger2015(const AliAnalysisTaskEMCALDCALTrigger2015 &);
  AliAnalysisTaskEMCALDCALTrigger2015 &operator=(const AliAnalysisTaskEMCALDCALTrigger2015 &);

  ClassDef(AliAnalysisTaskEMCALDCALTrigger2015, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif
