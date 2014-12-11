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
/*
 * Analysis component for EMCal clusters. Loops over calibrated and uncalibrated clusters
 *
 *   Author: Markus Fasel
 */
#include <map>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliVEvent.h"

#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerClusterAnalysisComponent.h"


ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerClusterAnalysisComponent)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerClusterAnalysisComponent::AliEMCalTriggerClusterAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent(),
  fEnergyRange(),
  fUsePatches(kFALSE)
{
  /*
   * Dummy (I/O) constructor
   */
}

//______________________________________________________________________________
AliEMCalTriggerClusterAnalysisComponent::AliEMCalTriggerClusterAnalysisComponent(const char* name) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fEnergyRange(),
  fUsePatches(kFALSE)
{
  /*
   * Main constructor
   */
  fEnergyRange.SetLimits(0., 1000.);
}

//______________________________________________________________________________
void AliEMCalTriggerClusterAnalysisComponent::CreateHistos() {
  /*
   * Create histos for clusters in different event categories
   */
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();

  // Create trigger definitions
  std::map<std::string, std::string> triggerCombinations;
  const char *triggernames[12] = {"MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh",
      "EMCGLow", "NoEMCal", "EMCHighBoth", "EMCHighGammaOnly", "EMCHighJetOnly",
      "EMCLowBoth", "EMCLowGammaOnly", "EMCLowJetOnly"};
  // Define names and titles for different triggers in the histogram container
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[0], "min. bias events"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[1], "jet-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[2], "jet-triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[3], "gamma-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[4], "gamma-triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[5], "non-EMCal-triggered events"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[6], "jet and gamma triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[7], "exclusively gamma-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[8], "exclusively jet-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[9], "jet and gamma triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[10], "exclusively gamma-triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[11], "exclusively-triggered events (low threshold)"));

  // Create axis definitions
  const AliEMCalTriggerBinningDimension *ptbinning = fBinning->GetBinning("pt"),
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
      DefineAxis("energy", ptbinning),
      DefineAxis("eta", etabinning),
      DefineAxis("phi", phibinning),
      DefineAxis("zvertex", vertexbinning),
      DefineAxis("mbtrigger", 2, -0.5, 1.5)
  };

  // Build histograms
  for(std::map<std::string,std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); ++it){
    const std::string name = it->first, &title = it->second;
    fHistos->CreateTHnSparse(Form("hClusterCalibHist%s", name.c_str()), Form("Calib. cluster-based histogram for %s events", title.c_str()), 5, clusteraxes, "s");
    fHistos->CreateTHnSparse(Form("hClusterUncalibHist%s", name.c_str()), Form("Uncalib. cluster-based histogram for %s events", title.c_str()), 5, clusteraxes, "s");
  }

  for(int iaxis = 0; iaxis < 5; iaxis++) delete clusteraxes[iaxis];
}

//______________________________________________________________________________
void AliEMCalTriggerClusterAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  /*
   * Run loop over calibrated and uncalibrated clusters
   */

  // First loop over uncalibrated clusters
  AliVCluster *clust(NULL);
  AliVEvent *recEv = data->GetRecEvent();
  std::vector<std::string> triggerNames;
  this->GetMachingTriggerNames(triggerNames, fUsePatches);
  for(int iclust = 0; iclust < recEv->GetNumberOfCaloClusters(); iclust++){
    clust = recEv->GetCaloCluster(iclust);
    if(!clust->IsEMCAL()) continue;
    if(!fEnergyRange.IsInRange(clust->E())) continue;
    for(std::vector<std::string>::iterator name = triggerNames.begin(); name != triggerNames.end(); ++name)
      FillHistogram(Form("hClusterUncalibHist%s", name->c_str()), clust, recEv, fTriggerDecision->IsMinBias());
  }

  // Loop also over calibrated clusters
  TIter clusterIter(data->GetClusterContainer());
  while((clust = dynamic_cast<AliVCluster *>(clusterIter()))){
    if(!clust->IsEMCAL()) continue;
    if(!fEnergyRange.IsInRange(clust->E())) continue;
    for(std::vector<std::string>::iterator name = triggerNames.begin(); name != triggerNames.end(); ++name)
      FillHistogram(Form("hClusterUncalibHist%s", name->c_str()), clust, recEv, fTriggerDecision->IsMinBias());
  }
}

//______________________________________________________________________________
void AliEMCalTriggerClusterAnalysisComponent::FillHistogram(const TString& histname, const AliVCluster* clust, AliVEvent *ev, Bool_t inMB) {
  /*
   * Fill Histogram for cluster
   *
   * @param histname: the histogram to fill
   * @param clust: the cluster analysed
   * @param event: reconstructed event information
   * @param inMB: true if event fulfills min bias condition
   */
  TLorentzVector vec;
  double xyz[3];
  ev->GetPrimaryVertex()->GetXYZ(xyz);
  clust->GetMomentum(vec, xyz);
  double infs[5] = {clust->E(), vec.Eta(), vec.Phi(), xyz[2], inMB ? 1. : 0.};
  fHistos->FillTHnSparse(histname.Data(), infs);
}


} /* namespace EMCalTriggerPtAnalysis */
