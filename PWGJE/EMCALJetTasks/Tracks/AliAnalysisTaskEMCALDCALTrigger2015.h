#ifndef ALIANALYSISTASKEMCALDCALTRIGGER2015_H_
#define ALIANALYSISTASKEMCALDCALTRIGGER2015_H_

#include <TString.h>

#include "AliAnalysisTaskSE.h"

class TArrayD;
class TClonesArray;
class THistManager;

class AliEMCALTriggerPatchInfo;
class AliVCluster;

namespace EMCalTriggerPtAnalysis {

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
  void ProcessPatch(const TString &triggerclass, const AliEMCALTriggerPatchInfo * const patch, bool isOnline);
  void CreateEnergyBinning(TArrayD& binning) const;
  void CreateLinearBinning(TArrayD& binning, int nbins, double min, double max) const;

  static const TString fgkTriggerClasses[11];
  static const TString fgkBeamDirs[4];

  TString                       fClusterContainerName;        //
  THistManager                  *fHistos;                     //!<!
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
