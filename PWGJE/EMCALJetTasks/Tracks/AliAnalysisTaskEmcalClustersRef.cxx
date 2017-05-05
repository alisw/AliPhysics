/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <array>
#include <bitset>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <TArrayD.h>
#include <TClonesArray.h>
#include <TGrid.h>
#include <THashList.h>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TParameter.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliClusterContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalList.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliVVertex.h"

#include "AliAnalysisTaskEmcalClustersRef.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef)
/// \endcond

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef() :
    AliAnalysisTaskEmcalTriggerBase(),
    fCentralityRange(-999., 999.),
    fRequestCentrality(false),
    fEventCentrality(-1),
    fCentralityEstimator("V0M"),
    fBunchCrossingIndex(-1),
    fEnergyDefinition(kDefaultEnergy),
    fEnableSumw2(false),
    fClusterTimeRange(-50e-6, 50e-6)
{
}

AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef(const char *name) :
    AliAnalysisTaskEmcalTriggerBase(name),
    fCentralityRange(-999., 999.),
    fRequestCentrality(false),
    fEventCentrality(-1),
    fCentralityEstimator("V0M"),
    fBunchCrossingIndex(-1),
    fEnergyDefinition(kDefaultEnergy),
    fEnableSumw2(false),
    fClusterTimeRange(-50e-6, 50e-6)
{
}

AliAnalysisTaskEmcalClustersRef::~AliAnalysisTaskEmcalClustersRef() {
}


void AliAnalysisTaskEmcalClustersRef::CreateUserHistos(){

  EnergyBinning energybinning;
  TLinearBinning smbinning(21, -0.5, 20.5), etabinning(100, -0.7, 0.7), timebinning(1000, -500e-9, 500e-9), ncellbinning(101, -0.5, 100.5);
  TString optionstring = fEnableSumw2 ? "s" : "";

  /*
   * Exclusive classes are defined as triggered events
   * in a class without lower threshold classes firing.
   * This is needed to make the triggers statistically
   * independent.
   */
  std::array<Double_t, 5> encuts = {1., 2., 5., 10., 20.};
  Int_t sectorsWithEMCAL[10] = {4, 5, 6, 7, 8, 9, 13, 14, 15, 16};

  // Binnings for Multiplicity correlation
  TLinearBinning v0abinning(1000, 0., 1000.), trackletbinning(500, 0., 500.), itsclustbinning(500, 0., 500.), emcclustbinning(100, 0., 100.), emccellbinning(3000, 0., 3000.);
  const TBinning *multbinning[6] = {&v0abinning, &trackletbinning, &trackletbinning, &itsclustbinning, &emcclustbinning, &emccellbinning};
  for(auto trg : GetSupportedTriggers()){
    fHistos->CreateTH1("hEventCount" + trg, "Event count for trigger class " + trg, 1, 0.5, 1.5, optionstring);
    fHistos->CreateTH1("hEventCentrality" + trg, "Event centrality for trigger class " + trg, 103, -2., 101., optionstring);
    fHistos->CreateTH1("hVertexZ" + trg, "z-position of the primary vertex for trigger class " + trg, 200, -40., 40., optionstring);
    fHistos->CreateTHnSparse("hMultiplicityCorrelation" + trg, "Multiplicity correlation for trigger" + trg, 6, multbinning);
    fHistos->CreateTH1("hClusterEnergy" + trg, "Cluster energy for trigger class " + trg, energybinning, optionstring);
    fHistos->CreateTH1("hClusterET" + trg, "Cluster transverse energy for trigger class " + trg, energybinning, optionstring);
    fHistos->CreateTH1("hClusterEnergyFired" + trg, "Cluster energy for trigger class " + trg + ", firing the trigger", energybinning, optionstring);
    fHistos->CreateTH1("hClusterETFired" + trg, "Cluster transverse energy for trigger class " + trg + ", firing the trigger" , energybinning, optionstring);
    fHistos->CreateTH2("hClusterEnergySM" + trg, "Cluster energy versus supermodule for trigger class " + trg, smbinning, energybinning, optionstring);
    fHistos->CreateTH2("hClusterETSM" + trg, "Cluster transverse energy versus supermodule for trigger class " + trg, smbinning, energybinning, optionstring);
    fHistos->CreateTH2("hClusterEnergyFiredSM" + trg, "Cluster energy versus supermodule for trigger class " + trg + ", firing the trigger" , smbinning, energybinning, optionstring);
    fHistos->CreateTH2("hClusterETFiredSM" + trg, "Cluster transverse energy versus supermodule for trigger class " + trg + ", firing the trigger" , smbinning, energybinning, optionstring);
    fHistos->CreateTH2("hEtaEnergy" + trg, "Cluster energy vs. eta for trigger class " + trg, etabinning, energybinning, optionstring);
    fHistos->CreateTH2("hEtaET" + trg, "Cluster transverse energy vs. eta for trigger class " + trg, etabinning, energybinning, optionstring);
    fHistos->CreateTH2("hTimeEnergy" + trg, "Cluster time vs. energy for trigger class " + trg, timebinning, energybinning, optionstring);
    fHistos->CreateTH2("hNCellEnergy" + trg, "Cluster number of cells vs energy for trigger class " + trg, ncellbinning, energybinning, optionstring);
    fHistos->CreateTH2("hNCellET" + trg, "Cluster number of cells vs transverse energy for trigger class " + trg, ncellbinning, energybinning, optionstring);
    fHistos->CreateTH2("hEtaEnergyFired" + trg, "Cluster energy vs. eta for trigger class " + trg + ", firing the trigger", etabinning, energybinning, optionstring);
    fHistos->CreateTH2("hEtaETFired" + trg, "Cluster transverse energy vs. eta for trigger class " + trg + ", firing the trigger", etabinning, energybinning, optionstring);
    for(int ism = 0; ism < 20; ism++){
      fHistos->CreateTH2(TString::Format("hEtaEnergySM%d", ism) + trg, TString::Format("Cluster energy vs. eta in Supermodule %d for trigger ", ism) + trg, etabinning, energybinning, optionstring);
      fHistos->CreateTH2(TString::Format("hEtaETSM%d", ism) + trg, TString::Format("Cluster transverse energy vs. eta in Supermodule %d for trigger ", ism) + trg, etabinning, energybinning, optionstring);
      fHistos->CreateTH2(TString::Format("hEtaEnergyFiredSM%d", ism) + trg, TString::Format("Cluster energy vs. eta in Supermodule %d for trigger ", ism) + trg + ",  firing the trigger", etabinning, energybinning, optionstring);
      fHistos->CreateTH2(TString::Format("hEtaETFiredSM%d", ism) + trg, TString::Format("Cluster transverse energy vs. eta in Supermodule %d for trigger ", ism) + trg +", firing the trigger", etabinning, energybinning, optionstring);
    }
    for(int isec = 0; isec < 10; isec++){
      fHistos->CreateTH2(TString::Format("hEtaEnergySec%d", sectorsWithEMCAL[isec]) + trg, TString::Format("Cluster energy vs.eta in tracking sector %d for trigger ", sectorsWithEMCAL[isec]) + trg, etabinning, energybinning, optionstring);
      fHistos->CreateTH2(TString::Format("hEtaETSec%d", sectorsWithEMCAL[isec]) + trg, TString::Format("Cluster transverse energy vs.eta in tracking sector %d for trigger ", sectorsWithEMCAL[isec]) + trg, etabinning, energybinning, optionstring);
      fHistos->CreateTH2(TString::Format("hEtaEnergyFiredSec%d", sectorsWithEMCAL[isec]) + trg, TString::Format("Cluster energy vs.eta in tracking sector %d for trigger ", sectorsWithEMCAL[isec]) +  trg + ", firing the trigger", etabinning, energybinning, optionstring);
      fHistos->CreateTH2(TString::Format("hEtaETFiredSec%d", sectorsWithEMCAL[isec]) + trg, TString::Format("Cluster transverse energy vs.eta in tracking sector %d for trigger ", sectorsWithEMCAL[isec]) +  trg + ", firing the trigger", etabinning, energybinning, optionstring);
    }
    for(auto ien : encuts){
      fHistos->CreateTH2(TString::Format("hEtaPhi%dG", static_cast<int>(ien)) + trg, TString::Format("cluster #eta-#phi map for clusters with energy larger than %f GeV/c for trigger class ", ien) + trg, 100, -0.7, 0.7, 200, 0, 2*TMath::Pi(), optionstring);
      fHistos->CreateTH2(TString::Format("hEtaPhiFired%dG", static_cast<int>(ien)) +  trg, TString::Format("cluster #eta-#phi map for clusters fired the trigger with energy larger than %f GeV/c for trigger class", ien) + trg + ", firing the trigger", 200, -0.7, 0.7, 200, 0, 2*TMath::Pi(), optionstring);
    }
  }
}

bool AliAnalysisTaskEmcalClustersRef::IsUserEventSelected(){
  fEventCentrality = -1;
  if(fRequestCentrality){
    AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if(!mult){
      AliErrorStream() << GetName() << ": Centrality selection enabled but no centrality estimator found" << std::endl;
      return false;
    }
    if(mult->IsEventSelected()) return false;
    fEventCentrality = mult->GetEstimator(fCentralityEstimator)->GetPercentile();
    AliDebugStream(1) << GetName() << ": Centrality " <<  fEventCentrality << std::endl;
    if(!fCentralityRange.IsInRange(fEventCentrality)){
      AliDebugStream(1) << GetName() << ": reject centrality: " << fEventCentrality << std::endl;
      return false;
    } else {
      AliDebugStream(1) << GetName() << ": select centrality " << fEventCentrality << std::endl;
    }
  } else {
    AliDebugStream(1) << GetName() << ": No centrality selection applied" << std::endl;
  }

  if(fBunchCrossingIndex > -1){
    int bcindex = fInputEvent->GetHeader()->GetBunchCrossNumber() % 4;
    if(bcindex != fBunchCrossingIndex) return false;
  }
  return true;
}

bool AliAnalysisTaskEmcalClustersRef::Run(){
  AliDebugStream(1) << GetName() << ": UserExec start" << std::endl;

  TList ej1patches, dj1patches, ej2patches, dj2patches, eg1patches, dg1patches, eg2patches, dg2patches;
  FindPatchesForTrigger("EJ2", fTriggerPatchInfo, ej2patches);
  FindPatchesForTrigger("DJ2", fTriggerPatchInfo, dj2patches);
  FindPatchesForTrigger("EJ1", fTriggerPatchInfo, ej1patches);
  FindPatchesForTrigger("DJ1", fTriggerPatchInfo, dj1patches);
  FindPatchesForTrigger("EG2", fTriggerPatchInfo, eg2patches);
  FindPatchesForTrigger("DG2", fTriggerPatchInfo, dg2patches);
  FindPatchesForTrigger("EG1", fTriggerPatchInfo, eg1patches);
  FindPatchesForTrigger("DG1", fTriggerPatchInfo, dg1patches);

  Double_t energy, et, eta, phi;
  TList *selpatches(nullptr);
  for(auto clust : GetClusterContainer(fNameClusterContainer.Data())->all()){
    //AliVCluster *clust = static_cast<AliVCluster *>(*clustIter);
    if(!clust->IsEMCAL()) continue;
    if(clust->GetIsExotic()) continue;
    if(!fClusterTimeRange.IsInRange(clust->GetTOF())) continue;

    // Distinguish energy definition
    switch(fEnergyDefinition){
    case kDefaultEnergy:
    	AliDebugStream(2) << GetName() << ": Using cluster energy definition: default" << std::endl;
    	energy = clust->E();
    	break;
    case kNonLinCorrEnergy:
    	AliDebugStream(2) << GetName() << ": Using cluster energy definition: corrected for non-linearity" << std::endl;
    	energy = clust->GetNonLinCorrEnergy();
    	break;
    case kHadCorrEnergy:
    	AliDebugStream(2) << GetName() << ": Using cluster energy definition: corrected for hadronic contribution" << std::endl;
    	energy = clust->GetHadCorrEnergy();
    	break;
    };

    AliDebugStream(2) << GetName() << ": Using energy " << energy << " (def: " << clust->E()
    		<< " | NL: " << clust->GetNonLinCorrEnergy()
			<< " | HD: " << clust->GetHadCorrEnergy()
			<< ")" << std::endl;

    TLorentzVector posvec;
    clust->GetMomentum(posvec, fVertex);
    posvec.SetE(energy);        // use energy definition as selected in the task
    et = posvec.Et();
    eta = posvec.Eta();
    phi = posvec.Phi();

    // fill histograms allEta
    for(const auto & trg : fSelectedTriggers){
      selpatches = nullptr;
      if(trg.Contains("EJ2")) selpatches = &ej2patches;
      if(trg.Contains("DJ2")) selpatches = &dj2patches;
      if(trg.Contains("EJ1")) selpatches = &ej1patches;
      if(trg.Contains("DJ1")) selpatches = &dj1patches;
      if(trg.Contains("EG2")) selpatches = &eg2patches;
      if(trg.Contains("DG2")) selpatches = &dg2patches;
      if(trg.Contains("EG1")) selpatches = &eg1patches;
      if(trg.Contains("DG1")) selpatches = &dg1patches;
      FillClusterHistograms(trg.Data(), energy, et, eta, phi, clust->GetTOF(), clust->GetNCells(), nullptr);
    }
  }
  return true;
}

void AliAnalysisTaskEmcalClustersRef::FillClusterHistograms(const TString &triggerclass, double energy, double transverseenergy, double eta, double phi, double clustertime, int ncell, TList *fTriggerPatches){
  Bool_t hasTriggerPatch = fTriggerPatches  ? CorrelateToTrigger(eta, phi, fTriggerPatches) : kFALSE;
  Int_t supermoduleID = -1, sector = -1;
  Double_t weight = GetTriggerWeight(triggerclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << triggerclass << std::endl;

  fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
  fHistos->FillTH1("hClusterEnergy" + triggerclass, energy, weight);
  fHistos->FillTH1("hClusterET" + triggerclass, transverseenergy, weight);
  fHistos->FillTH2("hEtaEnergy" + triggerclass, eta, energy, weight);
  fHistos->FillTH2("hEtaET" + triggerclass, eta, transverseenergy, weight);
  fHistos->FillTH2("hTimeEnergy" + triggerclass, clustertime, energy, weight);
  fHistos->FillTH2("hNCellEnergy" + triggerclass, ncell, energy, weight);
  fHistos->FillTH2("hNCellET" + triggerclass, ncell, transverseenergy, weight);
  if(supermoduleID >= 0){
    fHistos->FillTH2("hClusterEnergySM" + triggerclass, supermoduleID, energy, weight);
    fHistos->FillTH2("hClusterETSM" + triggerclass, supermoduleID, transverseenergy, weight);
    fHistos->FillTH2(TString::Format("hEtaEnergySM%d", supermoduleID) + triggerclass, eta, energy, weight);
    fHistos->FillTH2(TString::Format("hEtaETSM%d", supermoduleID) + triggerclass, eta, transverseenergy, weight);
    if(supermoduleID < 12)
      sector = 4 + int(supermoduleID/2); // EMCAL
    else
      sector = 13 + int((supermoduleID-12)/2);  // DCAL
    fHistos->FillTH2(TString::Format("hEtaEnergySec%d", sector) + triggerclass, eta, energy, weight);
    fHistos->FillTH2(TString::Format("hEtaETSec%d", sector) + triggerclass, eta, transverseenergy, weight);
  }
  if(hasTriggerPatch){
    fHistos->FillTH1("hClusterEnergyFired" + triggerclass, energy, weight);
    fHistos->FillTH1("hClusterETFired" + triggerclass, energy, weight);
    fHistos->FillTH2("hEtaEnergyFired" + triggerclass, eta, energy, weight);
    fHistos->FillTH2("hEtaETFired" + triggerclass, eta, energy, weight);
    if(supermoduleID >= 0){
      fHistos->FillTH2("hClusterEnergyFiredSM" + triggerclass, supermoduleID, energy, weight);
      fHistos->FillTH2("hClusterETFiredSM" + triggerclass, supermoduleID, transverseenergy, weight);
      fHistos->FillTH2(TString::Format("hEtaEnergyFiredSM%d", supermoduleID) + triggerclass, eta, energy,weight);
      fHistos->FillTH2(TString::Format("hEtaETFiredSM%d", supermoduleID) + triggerclass, eta, transverseenergy, weight);
      fHistos->FillTH2(TString::Format("hEtaEnergyFiredSec%d", sector) + triggerclass, eta, energy, weight);
      fHistos->FillTH2(TString::Format("hEtaETFiredSec%d", sector) + triggerclass, eta, transverseenergy, weight);
    }
  }
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(int ien = 0; ien < 5; ien++){
    if(energy > encuts[ien]){
      fHistos->FillTH2(TString::Format("hEtaPhi%dG", static_cast<int>(encuts[ien])) + triggerclass, eta, phi, weight);
      if(hasTriggerPatch){
        fHistos->FillTH2(TString::Format("hEtaPhiFired%dG", static_cast<int>(encuts[ien])) + triggerclass, eta, phi, weight);
      }
    }
  }
}

void AliAnalysisTaskEmcalClustersRef::UserFillHistosAfterEventSelection(){
  double v0amult = fInputEvent->GetVZEROData()->GetMTotV0A(),
         trackletmult = static_cast<double>(CountTracklets(-0.8, 0.8, 0., TMath::TwoPi())),
         emctrackletmult = static_cast<double>(CountTracklets(-0.8, 0.8, 1.4, TMath::Pi())),
         itsclustermult = fInputEvent->GetMultiplicity()->GetNumberOfSPDClusters(),
         emcclustermult = static_cast<double>(CountEmcalClusters(0.5)),
         emccellocc = static_cast<double>(this->GetEMCALCellOccupancy(0.1));
  for(const auto &t : fSelectedTriggers){
    Double_t weight = GetTriggerWeight(t);
    fHistos->FillTH1("hEventCount" + t, 1, weight);
    fHistos->FillTH1("hEventCentrality" + t, fEventCentrality, weight);
    fHistos->FillTH1("hVertexZ" + t, fVertex[2], weight);

    // Multiplicity correlation (no correction for downscaling)
    double data[6] = {v0amult, trackletmult, emctrackletmult, itsclustermult, emcclustermult, emccellocc};
    fHistos->FillTHnSparse("hMultiplicityCorrelation" + t, data);
  }
}

Bool_t AliAnalysisTaskEmcalClustersRef::CorrelateToTrigger(Double_t etaclust, Double_t phiclust, TList *fTriggerPatches) const {
  Bool_t hasfound = kFALSE;
  for(TIter patchIter = TIter(fTriggerPatches).Begin(); patchIter != TIter::End(); ++patchIter){
    Double_t boundaries[4];
    GetPatchBoundaries(*patchIter, boundaries);
    Double_t etamin = TMath::Min(boundaries[0], boundaries[1]),
        etamax = TMath::Max(boundaries[0], boundaries[1]),
        phimin = TMath::Min(boundaries[2], boundaries[3]),
        phimax = TMath::Max(boundaries[2], boundaries[3]);
    if(etaclust > etamin && etaclust < etamax && phiclust > phimin && phiclust < phimax){
      hasfound = kTRUE;
      break;
    }
  }
  return hasfound;
}

void AliAnalysisTaskEmcalClustersRef::FindPatchesForTrigger(TString triggerclass, const TClonesArray * triggerPatches, TList &foundtriggers) const {
  foundtriggers.Clear();
  if(!triggerPatches) return;
  AliEmcalTriggerOfflineSelection::EmcalTriggerClass myclass = AliEmcalTriggerOfflineSelection::kTrgEL0;
  if(triggerclass == "EG1") myclass = AliEmcalTriggerOfflineSelection::kTrgEG1;
  if(triggerclass == "EG2") myclass = AliEmcalTriggerOfflineSelection::kTrgEG2;
  if(triggerclass == "EJ1") myclass = AliEmcalTriggerOfflineSelection::kTrgEJ1;
  if(triggerclass == "EJ2") myclass = AliEmcalTriggerOfflineSelection::kTrgEJ2;
  if(triggerclass == "DMC7") myclass = AliEmcalTriggerOfflineSelection::kTrgDL0;
  if(triggerclass == "DG1") myclass = AliEmcalTriggerOfflineSelection::kTrgDG1;
  if(triggerclass == "DG2") myclass = AliEmcalTriggerOfflineSelection::kTrgDG2;
  if(triggerclass == "DJ1") myclass = AliEmcalTriggerOfflineSelection::kTrgDJ1;
  if(triggerclass == "DJ2") myclass = AliEmcalTriggerOfflineSelection::kTrgDJ2;
  for(auto patchiter : *triggerPatches){
    AliEMCALTriggerPatchInfo *mypatch = static_cast<AliEMCALTriggerPatchInfo *>(patchiter);
    if(!mypatch->IsOfflineSimple()) continue;
    if(AliEmcalTriggerOfflineSelection::IsDCAL(myclass)){
      if(!mypatch->IsDCalPHOS()) continue;
    } else {
      if(mypatch->IsDCalPHOS()) continue;
    }
    if(AliEmcalTriggerOfflineSelection::IsSingleShower(myclass)){
      if(!mypatch->IsGammaLowSimple()) continue;
    } else {
      if(!mypatch->IsJetLowSimple()) continue;
    }
    double threshold = fTriggerSelection ? fTriggerSelection->GetThresholdForTrigger(myclass) : -1;
    if(mypatch->GetPatchE() > threshold) foundtriggers.Add(patchiter);
  }
}

void AliAnalysisTaskEmcalClustersRef::GetPatchBoundaries(TObject *o, Double_t *boundaries) const {
  AliEMCALTriggerPatchInfo *patch= dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  boundaries[0] = patch->GetEtaMin();
  boundaries[1] = patch->GetEtaMax();
  boundaries[2] = patch->GetPhiMin();
  boundaries[3] = patch->GetPhiMax();
}

int AliAnalysisTaskEmcalClustersRef::CountEmcalClusters(double ecut){
	int nclusters = 0;
	for(auto clust : GetClusterContainer(fNameClusterContainer.Data())->all()){
	  if(!clust->IsEMCAL()) continue;
	  if(clust->GetIsExotic()) continue;
	  if(clust->E() > ecut) nclusters++;
	  nclusters++;
	}
	return nclusters;
}

int AliAnalysisTaskEmcalClustersRef::CountTracklets(double etamin, double etamax, double phimin, double phimax){
  int ntracklets = 0;
  AliVMultiplicity *mult = fInputEvent->GetMultiplicity();
  for(int itl = 0; itl < mult->GetNumberOfTracklets(); itl++){
    double eta = mult->GetEta(itl), phi = mult->GetPhi(itl);
    if(!(eta > etamin && eta < etamax)) continue;
    if(!(phi > phimin && phi < phimax)) continue;
    ntracklets++;
  }
  return ntracklets;
}

int AliAnalysisTaskEmcalClustersRef::GetEMCALCellOccupancy(double ecut){
  std::set<int> cellIDs;
  AliVCaloCells *emccells = fInputEvent->GetEMCALCells();
  for(short icell = 0; icell < emccells->GetNumberOfCells(); icell++){
    if(emccells->GetAmplitude(icell) > ecut){
      int cellID = emccells->GetCellNumber(icell);
      if(cellIDs.find(cellID) == cellIDs.end()) cellIDs.insert(cellID);
    }
  }
  return cellIDs.size();
}

AliAnalysisTaskEmcalClustersRef *AliAnalysisTaskEmcalClustersRef::AddTaskEmcalClustersRef(const TString &nclusters, const TString &suffix){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString clusName(nclusters == "usedefault" ? AliEmcalAnalysisFactory::ClusterContainerNameFactory(mgr->GetInputEventHandler()->InheritsFrom("AliAODInputHandler")) : nclusters);

  TString taskname = "emcalClusterQA_" + suffix;

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef(taskname.Data());
  task->AddClusterContainer(clusName);
  task->SetClusterContainer(clusName);
  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ClusterQA_" + TString(suffix);
  TString containername = "ClusterResults_" + TString(suffix);
  printf("Outfile: %s, container: %s\n", outfile.Data(), containername.Data());

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(containername.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}

AliAnalysisTaskEmcalClustersRef *AliAnalysisTaskEmcalClustersRef::AddTaskEmcalClustersRefDefault(const TString &nClusters){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef("emcalClusterQA");
  mgr->AddTask(task);

  // Adding cluster container
  TString clusName(nClusters == "usedefault" ? AliEmcalAnalysisFactory::ClusterContainerNameFactory(mgr->GetInputEventHandler()->InheritsFrom("AliAODInputHandler")) : nClusters);
  task->AddClusterContainer(clusName.Data());
  task->SetClusterContainer(clusName);

  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // These numers are determined as:
  // EMC7: 3.5 GeV
  // EG1:  14 GeV
  // EG2:  8 GeV
  // EJ1:  22 GeV
  // EJ2:  12 GeV
  task->SetOfflineTriggerSelection(
      EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TriggerSelectionFactory(5, 14, 8, 22, 12)
  );

  TString outfile(mgr->GetCommonFileName());
  outfile += ":ClusterQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("ClusterResults", AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}


/**
 * Create new energy binning
 */
AliAnalysisTaskEmcalClustersRef::EnergyBinning::EnergyBinning():
  TCustomBinning()
{
  this->SetMinimum(0.);
  this->AddStep(1, 0.05);
  this->AddStep(2, 0.1);
  this->AddStep(4, 0.2);
  this->AddStep(7, 0.5);
  this->AddStep(16, 1);
  this->AddStep(32, 2);
  this->AddStep(40, 4);
  this->AddStep(50, 5);
  this->AddStep(100, 10);
  this->AddStep(200, 20);
}


} /* namespace EMCalTriggerPtAnalysis */
