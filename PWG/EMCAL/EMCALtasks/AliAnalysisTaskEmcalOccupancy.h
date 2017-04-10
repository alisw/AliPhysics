#ifndef ALIANALYSISTASKEMCALOCCUPANCY_H
#define ALIANALYSISTASKEMCALOCCUPANCY_H

#include "AliAnalysisTaskEmcalLight.h"

#include <TString.h>
class THistManager;

/**
 * @class AliAnalysisTaskEmcalOccupancy
 * @brief Task monitoring the occupancy on cluster and cell level
 */
class AliAnalysisTaskEmcalOccupancy : public AliAnalysisTaskEmcalLight {
public:
  AliAnalysisTaskEmcalOccupancy();
  AliAnalysisTaskEmcalOccupancy(const char *name);
  virtual ~AliAnalysisTaskEmcalOccupancy();

  void SetUseCentralityPercentile(Bool_t doUse) { fUseCentrality = doUse; }
  void SetNameClusters(const char *name) { fNameClusters = name; }

  static AliAnalysisTaskEmcalOccupancy *AddOccupancyTask(const char *name);

protected:
  virtual void UserCreateOutputObjects();
  virtual void ExecOnce();
  virtual bool Run();

private:
  AliAnalysisTaskEmcalOccupancy(const AliAnalysisTaskEmcalOccupancy &);
  AliAnalysisTaskEmcalOccupancy &operator=(const AliAnalysisTaskEmcalOccupancy &);

  TString                           fNameClusters;      ///< Name of the cluster container
  THistManager                      *fHistos;           ///< Histogram container;
  Bool_t                            fUseCentrality;     ///< Switch whether to use the centrality percentile
  UChar_t                           *fCellCounter;      //!<! Counting how often a cell is fired per event

  ClassDef(AliAnalysisTaskEmcalOccupancy, 1);
};

#endif /* ALIANALYSISTASKEMCALOCCUPANCY_H */
