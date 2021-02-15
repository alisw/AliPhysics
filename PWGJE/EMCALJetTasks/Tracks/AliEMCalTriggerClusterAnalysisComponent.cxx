/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <map>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TClonesArray.h>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>

#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVEvent.h"

#include "AliEMCalTriggerAnaClassManager.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerClusterAnalysisComponent.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerClusterAnalysisComponent)

using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy (I/O) constructor, not to be used by the user
 */
AliEMCalTriggerClusterAnalysisComponent::AliEMCalTriggerClusterAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent(),
  fEnergyRange()
{
}

/**
 * Main constructor, initializes all elements with default values
 * \param name Name of the component
 */
AliEMCalTriggerClusterAnalysisComponent::AliEMCalTriggerClusterAnalysisComponent(const char* name) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fEnergyRange()
{
  fEnergyRange.SetLimits(0., 1000.);
}

/**
 * Create histos for clusters in different event categories
 */
void AliEMCalTriggerClusterAnalysisComponent::CreateHistos() {
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();

  // Create trigger definitions
  std::map<std::string, std::string> triggerCombinations;
  GetAllTriggerNamesAndTitles(triggerCombinations);

  // Create axis definitions
  const TBinning *ptbinning = fBinning->GetBinning("pt"),
      *etabinning = fBinning->GetBinning("eta"),
      *phibinning = fBinning->GetBinning("phi"),
      *vertexbinning = fBinning->GetBinning("zvertex");
  /*
  ptbinning->Print();
  etabinning->Print();
  phibinning->Print();
  vertexbinning->Print();
  */

  const TAxis *clusteraxes[5] = {
      DefineAxis("energy", *ptbinning),
      DefineAxis("eta", *etabinning),
      DefineAxis("phi", *phibinning),
      DefineAxis("zvertex", *vertexbinning),
      DefineAxis("mbtrigger", TLinearBinning(2, -0.5, 1.5))
  };

  // Build histograms
  for(std::map<std::string,std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); ++it){
    const std::string name = it->first, &title = it->second;
    fHistos->CreateTHnSparse(Form("hClusterCalibHist%s", name.c_str()), Form("Calib. cluster-based histogram for %s events", title.c_str()), 5, clusteraxes, "s");
    fHistos->CreateTHnSparse(Form("hClusterUncalibHist%s", name.c_str()), Form("Uncalib. cluster-based histogram for %s events", title.c_str()), 5, clusteraxes, "s");
  }

  for(int iaxis = 0; iaxis < 5; iaxis++) delete clusteraxes[iaxis];
}

/**
 * Run loop over calibrated and uncalibrated clusters
 * \param data All data of the event
 */
void AliEMCalTriggerClusterAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {

  // First loop over uncalibrated clusters
  AliDebug(1, Form("Number of calibrated clusters: %d", data->GetClusterContainer()->GetEntries()));

  AliVCluster *clust(NULL);
  AliVEvent *recEv = data->GetRecEvent();
  std::vector<std::string> triggerNames;
  GetMachingTriggerNames(triggerNames);
  for(int iclust = 0; iclust < recEv->GetNumberOfCaloClusters(); iclust++){
    clust = recEv->GetCaloCluster(iclust);
    if(!clust->IsEMCAL()) continue;
    if(!fEnergyRange.IsInRange(clust->E())) continue;
    for(std::vector<std::string>::iterator name = triggerNames.begin(); name != triggerNames.end(); ++name)
      FillHistogram(Form("hClusterUncalibHist%s", name->c_str()), clust, recEv, fTriggerClassManager->HasMinBiasTrigger());
  }

  // Loop also over calibrated clusters
  if(!data->GetClusterContainer())
    printf("Cluster container not found \n");
  TIter clusterIter(data->GetClusterContainer());
  while((clust = dynamic_cast<AliVCluster *>(clusterIter()))){
    if(!clust->IsEMCAL()) continue;
    if(!fEnergyRange.IsInRange(clust->E())) continue;
    for(std::vector<std::string>::iterator name = triggerNames.begin(); name != triggerNames.end(); ++name)
      FillHistogram(Form("hClusterCalibHist%s", name->c_str()), clust, recEv, fTriggerClassManager->HasMinBiasTrigger());
  }
}

/**
 * Fill Histogram for cluster
 * \param histname the histogram to fill
 * \param clust the cluster analysed
 * \param event reconstructed event information
 * \param inMB true if event fulfills min bias condition
 */
void AliEMCalTriggerClusterAnalysisComponent::FillHistogram(const TString& histname, const AliVCluster* clust, AliVEvent *ev, Bool_t inMB) {
  TLorentzVector vec;
  double xyz[3];
  ev->GetPrimaryVertex()->GetXYZ(xyz);
  clust->GetMomentum(vec, xyz);
  double infs[5] = {clust->E(), vec.Eta(), vec.Phi(), xyz[2], inMB ? 1. : 0.};
  fHistos->FillTHnSparse(histname.Data(), infs);
}
