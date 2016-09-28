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
#include <sstream>
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

#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliClusterContainer.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliOADBContainer.h"
#include "AliVCluster.h"
#include "AliVVertex.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"

#include "AliAnalysisTaskEmcalClustersRef.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalClustersRef)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy (I/O) constructor
 */
AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef() :
    AliAnalysisTaskEmcal(),
    fTriggerSelection(nullptr),
    fHistos(nullptr),
    fAcceptTriggers(),
    fNameClusterContainer(""),
    fRequestAnalysisUtil(kTRUE),
    fTriggerStringFromPatches(kFALSE),
    fCentralityRange(-999., 999.),
    fVertexRange(-999., 999.),
    fRequestCentrality(false),
    fNameDownscaleOADB(""),
    fDownscaleOADB(nullptr),
    fDownscaleFactors(nullptr)
{
  SetNeedEmcalGeom(true);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

/**
 * Named constructor
 * @param name Name of the task
 */
AliAnalysisTaskEmcalClustersRef::AliAnalysisTaskEmcalClustersRef(const char *name) :
    AliAnalysisTaskEmcal(name, kTRUE),
    fTriggerSelection(nullptr),
    fHistos(nullptr),
    fAcceptTriggers(),
    fNameClusterContainer(""),
    fRequestAnalysisUtil(kTRUE),
    fTriggerStringFromPatches(kFALSE),
    fCentralityRange(-999., 999.),
    fVertexRange(-999., 999.),
    fRequestCentrality(false),
    fNameDownscaleOADB(""),
    fDownscaleOADB(nullptr),
    fDownscaleFactors(nullptr)
{
  SetNeedEmcalGeom(true);
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalClustersRef::~AliAnalysisTaskEmcalClustersRef() {
  if(fTriggerSelection) delete fTriggerSelection;
}

/**
 * Creates output histograms: distribution of cluster energy for different trigger classes and number of events
 */
void AliAnalysisTaskEmcalClustersRef::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if(fRequestAnalysisUtil && !fAliAnalysisUtils){
    AliInfoStream() << "Creating histograms for task " << GetName() << std::endl;
    fAliAnalysisUtils = new AliAnalysisUtils;
  }

  EnergyBinning energybinning;
  TLinearBinning smbinning(21, -0.5, 20.5), etabinning(100, -0.7, 0.7);

  fHistos = new THistManager("Ref");
  /*
   * Exclusive classes are defined as triggered events
   * in a class without lower threshold classes firing.
   * This is needed to make the triggers statistically
   * independent.
   */
  std::array<TString, 21> triggers = {
      "MB", "EMC7", "DMC7",
      "EJ1", "EJ2", "EG1", "EG2", "DJ1", "DJ2", "DG1", "DG2",
      "EMC7excl", "DMC7excl", "EG2excl", "EJ2excl", "DG2excl", "DJ2excl",
      "EJ1excl", "DJ1excl", "EG1excl", "DG1excl"
  };
  std::array<Double_t, 5> encuts = {1., 2., 5., 10., 20.};
  Int_t sectorsWithEMCAL[10] = {4, 5, 6, 7, 8, 9, 13, 14, 15, 16};
  for(auto trg : triggers){
    fHistos->CreateTH1(Form("hEventCount%s", trg.Data()), Form("Event count for trigger class %s", trg.Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hEventCentrality%s", trg.Data()), Form("Event centrality for trigger class %s", trg.Data()), 103, -2., 101.);
    fHistos->CreateTH1(Form("hVertexZ%s", trg.Data()), Form("z-position of the primary vertex for trigger class %s", trg.Data()), 200, -40., 40.);
    fHistos->CreateTH1(Form("hClusterEnergy%s", trg.Data()), Form("Cluster energy for trigger class %s", trg.Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterET%s", trg.Data()), Form("Cluster transverse energy for trigger class %s", trg.Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterEnergyFired%s", trg.Data()), Form("Cluster energy for trigger class %s, firing the trigger", trg.Data()), energybinning);
    fHistos->CreateTH1(Form("hClusterETFired%s", trg.Data()), Form("Cluster transverse energy for trigger class %s, firing the trigger", trg.Data()), energybinning);
    fHistos->CreateTH2(Form("hClusterEnergySM%s", trg.Data()), Form("Cluster energy versus supermodule for trigger class %s", trg.Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterETSM%s", trg.Data()), Form("Cluster transverse energy versus supermodule for trigger class %s", trg.Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterEnergyFiredSM%s", trg.Data()), Form("Cluster energy versus supermodule for trigger class %s, firing the trigger", trg.Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hClusterETFiredSM%s", trg.Data()), Form("Cluster transverse energy versus supermodule for trigger class %s, firing the trigger", trg.Data()), smbinning, energybinning);
    fHistos->CreateTH2(Form("hEtaEnergy%s", trg.Data()), Form("Cluster energy vs. eta for trigger class %s", trg.Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaET%s", trg.Data()), Form("Cluster transverse energy vs. eta for trigger class %s", trg.Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaEnergyFired%s", trg.Data()), Form("Cluster energy vs. eta for trigger class %s, firing the trigger", trg.Data()), etabinning, energybinning);
    fHistos->CreateTH2(Form("hEtaETFired%s", trg.Data()), Form("Cluster transverse energy vs. eta for trigger class %s, firing the trigger", trg.Data()), etabinning, energybinning);
    for(int ism = 0; ism < 20; ism++){
      fHistos->CreateTH2(Form("hEtaEnergySM%d%s", ism, trg.Data()), Form("Cluster energy vs. eta in Supermodule %d for trigger %s", ism, trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETSM%d%s", ism, trg.Data()), Form("Cluster transverse energy vs. eta in Supermodule %d for trigger %s", ism, trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaEnergyFiredSM%d%s", ism, trg.Data()), Form("Cluster energy vs. eta in Supermodule %d for trigger %s, firing the trigger", ism, trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETFiredSM%d%s", ism, trg.Data()), Form("Cluster transverse energy vs. eta in Supermodule %d for trigger %s, firing the trigger", ism, trg.Data()), etabinning, energybinning);
    }
    for(int isec = 0; isec < 10; isec++){
      fHistos->CreateTH2(Form("hEtaEnergySec%d%s", sectorsWithEMCAL[isec], trg.Data()), Form("Cluster energy vs.eta in tracking sector %d for trigger %s", sectorsWithEMCAL[isec], trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETSec%d%s", sectorsWithEMCAL[isec], trg.Data()), Form("Cluster transverse energy vs.eta in tracking sector %d for trigger %s", sectorsWithEMCAL[isec], trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaEnergyFiredSec%d%s", sectorsWithEMCAL[isec], trg.Data()), Form("Cluster energy vs.eta in tracking sector %d for trigger %s, firing the trigger", sectorsWithEMCAL[isec], trg.Data()), etabinning, energybinning);
      fHistos->CreateTH2(Form("hEtaETFiredSec%d%s", sectorsWithEMCAL[isec], trg.Data()), Form("Cluster transverse energy vs.eta in tracking sector %d for trigger %s, firing the trigger", sectorsWithEMCAL[isec], trg.Data()), etabinning, energybinning);
    }
    for(auto ien : encuts){
      fHistos->CreateTH2(Form("hEtaPhi%dG%s", static_cast<int>(ien), trg.Data()), Form("cluster #eta-#phi map for clusters with energy larger than %f GeV/c for trigger class %s", ien, trg.Data()), 100, -0.7, 0.7, 200, 0, 2*TMath::Pi());
      fHistos->CreateTH2(Form("hEtaPhiFired%dG%s", static_cast<int>(ien), trg.Data()), Form("cluster #eta-#phi map for clusters fired the trigger with energy larger than %f GeV/c for trigger class %s", ien, trg.Data()), 200, -0.7, 0.7, 200, 0, 2*TMath::Pi());
    }
  }
  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);
  PostData(1, fOutput);
  AliDebugStream(1) << "End creating histograms" << std::endl;
}

bool AliAnalysisTaskEmcalClustersRef::IsEventSelected(){
  fAcceptTriggers.clear();
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  std::stringstream triggerdebug;
  triggerdebug << "Offline bits: " << std::bitset<sizeof(UInt_t) * 8>(selectionstatus);
  AliDebug(2, triggerdebug.str().c_str());
  TString triggerstring = "";
  if(fTriggerStringFromPatches){
    triggerstring = GetFiredTriggerClassesFromPatches(fTriggerPatchInfo);
  } else {
    triggerstring = fInputEvent->GetFiredTriggerClasses();
  }
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1"),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2"),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1"),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2"),
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("EMC7"),
      isDJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ1"),
      isDJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("DJ2"),
      isDG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG1"),
      isDG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("DG2"),
      isDMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("DMC7");
  if(fTriggerPatchInfo && fTriggerSelection){
      isEJ1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ1, fTriggerPatchInfo);
      isEJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ2, fTriggerPatchInfo);
      isEG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG1, fTriggerPatchInfo);
      isEG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG2, fTriggerPatchInfo);
      isEMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEL0, fTriggerPatchInfo);
      isDJ1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ1, fTriggerPatchInfo);
      isDJ2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDJ2, fTriggerPatchInfo);
      isDG1 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG1, fTriggerPatchInfo);
      isDG2 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDG2, fTriggerPatchInfo);
      isDMC7 &= fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgDL0, fTriggerPatchInfo);
  }
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2 || isDMC7 || isDG1 || isDG2 || isDJ1 || isDJ2)){
    AliDebugStream(1) << GetName() << ": Reject trigger" << std::endl;
    return false;
  }

  double centrality = -1;
  AliDebugStream(1) << "Event selected" << std::endl;
  if(fRequestCentrality){
    AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if(!mult){
      AliErrorStream() << GetName() << ": Centrality selection enabled but no centrality estimator found" << std::endl;
      return false;
    }
    if(mult->IsEventSelected()) return false;
    centrality = mult->GetEstimator("V0M")->GetPercentile();
    AliDebugStream(1) << GetName() << ": Centrality " <<  centrality << std::endl;
    if(!fCentralityRange.IsInRange(centrality)){
      AliDebugStream(1) << GetName() << ": reject centrality: " << centrality << std::endl;
      return false;
    } else {
      AliDebugStream(1) << GetName() << ": select centrality " << centrality << std::endl;
    }
  } else {
    AliDebugStream(1) << GetName() << ": No centrality selection applied" << std::endl;
  }
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  if(!vtx) vtx = fInputEvent->GetPrimaryVertexSPD();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1){
    AliDebug(1, Form("%s: Reject contributors\n", GetName()));
    return false;
  }
  // Fill reference distribution for the primary vertex before any z-cut
  if(fRequestAnalysisUtil){
    AliDebugStream(1) << GetName() << " : Reject analysis util" << std::endl;
    if(fInputEvent->IsA() == AliESDEvent::Class() && fAliAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return false;
    if(!fAliAnalysisUtils->IsVertexSelected2013pA(fInputEvent)) return false;       // Apply new vertex cut
    if(fAliAnalysisUtils->IsPileUpEvent(fInputEvent)) return false;       // Apply new vertex cut
  }
  // Apply vertex z cut
  if(!fVertexRange.IsInRange(vtx->GetZ())){
    AliDebugStream(1) << GetName() << ": Reject z[" << vtx->GetZ() << "]" << std::endl;
    return false;
  }

  // Do MC outlier cut
  if(fIsPythia){
    if(!CheckMCOutliers()){
      AliDebugStream(1) << GetName() << ": Reject MC outliers" << std::endl;
      return false;
    }
  }

  AliDebugStream(1) << GetName() << ": Event Selected" << std::endl;

  // Store trigger decision
  if(isMinBias) fAcceptTriggers.push_back("MB");
  if(isEMC7){
    fAcceptTriggers.push_back("EMC7");
    if(!isMinBias) fAcceptTriggers.push_back("EMC7excl");
  }
  if(isDMC7){
    fAcceptTriggers.push_back("DMC7");
    if(!isMinBias) fAcceptTriggers.push_back("DMC7excl");
  }
  if(isEG2){
    fAcceptTriggers.push_back("EG2");
    if(!(isMinBias || isEMC7)) fAcceptTriggers.push_back("EG2excl");
  }
  if(isEG1){
    fAcceptTriggers.push_back("EG1");
    if(!(isMinBias || isEMC7 || isEG2)) fAcceptTriggers.push_back("EG1excl");
  }
  if(isDG2){
    fAcceptTriggers.push_back("DG2");
    if(!(isMinBias || isDMC7)) fAcceptTriggers.push_back("DG2excl");
  }
  if(isDG1){
    fAcceptTriggers.push_back("DG1");
    if(!(isMinBias || isDMC7 || isDG2)) fAcceptTriggers.push_back("DG1excl");
  }
  if(isEJ2){
    fAcceptTriggers.push_back("EJ2");
    if(!(isMinBias || isEMC7)) fAcceptTriggers.push_back("EJ2excl");
  }
  if(isEJ1){
    fAcceptTriggers.push_back("EJ1");
    if(!(isMinBias || isEMC7 || isEJ2)) fAcceptTriggers.push_back("EJ1excl");
  }
  if(isDJ2){
    fAcceptTriggers.push_back("DJ2");
    if(!(isMinBias || isDMC7)) fAcceptTriggers.push_back("DJ2excl");
  }
  if(isDJ1){
    fAcceptTriggers.push_back("DJ1");
    if(!(isMinBias || isDMC7 || isDJ2)) fAcceptTriggers.push_back("DJ1excl");
  }

  return true;
}


/**
 *
 * @param
 */
bool AliAnalysisTaskEmcalClustersRef::Run(){
  AliDebugStream(1) << GetName() << ": UserExec start" << std::endl;

  double centrality = -1;
  AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
  if(mult && mult->IsEventSelected()) centrality = mult->GetEstimator("V0M")->GetPercentile();


  // Fill Event counter and reference vertex distributions for the different trigger classes
  for(const auto &trg : fAcceptTriggers) FillEventHistograms(trg.c_str(), centrality, fVertex[2]);


  /*
  TList *objects = fInputEvent->GetList();
  for(TIter objiter = TIter(objects).Begin(); objiter != TIter::End(); ++ objiter){
    printf("Object %s\n", (*objiter)->GetName());
  }
  */

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

    TLorentzVector posvec;
    energy = clust->GetNonLinCorrEnergy();
    et = posvec.Et();
    clust->GetMomentum(posvec, fVertex);
    eta = posvec.Eta();
    phi = posvec.Phi();

    // fill histograms allEta
    for(const auto & trg : fAcceptTriggers){
      selpatches = nullptr;
      if(trg.find("EJ2") != std::string::npos) selpatches = &ej2patches;
      if(trg.find("DJ2") != std::string::npos) selpatches = &dj2patches;
      if(trg.find("EJ1") != std::string::npos) selpatches = &ej1patches;
      if(trg.find("DJ1") != std::string::npos) selpatches = &dj1patches;
      if(trg.find("EG2") != std::string::npos) selpatches = &eg2patches;
      if(trg.find("DG2") != std::string::npos) selpatches = &dg2patches;
      if(trg.find("EG1") != std::string::npos) selpatches = &eg1patches;
      if(trg.find("DG1") != std::string::npos) selpatches = &dg1patches;
      FillClusterHistograms(trg.c_str(), energy, et, eta, phi, nullptr);
    }
  }
  return true;
}

void AliAnalysisTaskEmcalClustersRef::ExecOnce(){
  AliAnalysisTaskEmcal::ExecOnce();
  if(!fLocalInitialized) return;

  // Handle OADB container with downscaling factors
  if(fNameDownscaleOADB.Length()){
    if(fNameDownscaleOADB.Contains("alien://") && ! gGrid) TGrid::Connect("alien://");
    fDownscaleOADB = new AliOADBContainer("AliEmcalDownscaleFactors");
    fDownscaleOADB->InitFromFile(fNameDownscaleOADB.Data(), "AliEmcalDownscaleFactors");
  }
}

void AliAnalysisTaskEmcalClustersRef::RunChanged(Int_t runnumber){
 if(fDownscaleOADB){
    fDownscaleFactors = static_cast<TObjArray *>(fDownscaleOADB->GetObject(runnumber));
  }
}

/**
 * Get a trigger class dependent event weight. The weight
 * is defined as 1/downscalefactor. The downscale factor
 * is taken from the OADB. For triggers which are not downscaled
 * the weight is always 1.
 * @param[in] triggerclass Class for which to obtain the trigger.
 * @return Downscale facror for the trigger class (1 if trigger is not downscaled or no OADB container is available)
 */
Double_t AliAnalysisTaskEmcalClustersRef::GetTriggerWeight(const TString &triggerclass) const {
  if(fDownscaleFactors){
    TParameter<double> *result(nullptr);
    // Downscaling only done on MB, L0 and the low threshold triggers
    if(triggerclass.Contains("MB")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("INT7"));
    else if(triggerclass.Contains("EMC7")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EMC7"));
    else if(triggerclass.Contains("EJ2")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EJ2"));
    else if(triggerclass.Contains("EG2")) result = static_cast<TParameter<double> *>(fDownscaleFactors->FindObject("EG2"));
    if(result) return 1./result->GetVal();
  }
  return 1.;
}

void AliAnalysisTaskEmcalClustersRef::FillClusterHistograms(const TString &triggerclass, double energy, double transverseenergy, double eta, double phi, TList *fTriggerPatches){
  Bool_t hasTriggerPatch = fTriggerPatches  ? CorrelateToTrigger(eta, phi, fTriggerPatches) : kFALSE;
  Int_t supermoduleID = -1, sector = -1;
  Double_t weight = GetTriggerWeight(triggerclass);
  AliDebugStream(1) << GetName() << ": Using weight " << weight << " for trigger " << triggerclass << std::endl;

  fGeom->SuperModuleNumberFromEtaPhi(eta, phi, supermoduleID);
  fHistos->FillTH1(Form("hClusterEnergy%s", triggerclass.Data()), energy, weight);
  fHistos->FillTH1(Form("hClusterET%s", triggerclass.Data()), transverseenergy, weight);
  fHistos->FillTH2(Form("hEtaEnergy%s", triggerclass.Data()), eta, energy, weight);
  fHistos->FillTH2(Form("hEtaET%s", triggerclass.Data()), eta, transverseenergy, weight);
  if(supermoduleID >= 0){
    fHistos->FillTH2(Form("hClusterEnergySM%s", triggerclass.Data()), supermoduleID, energy, weight);
    fHistos->FillTH2(Form("hClusterETSM%s", triggerclass.Data()), supermoduleID, transverseenergy, weight);
    fHistos->FillTH2(Form("hEtaEnergySM%d%s", supermoduleID, triggerclass.Data()), eta, energy, weight);
    fHistos->FillTH2(Form("hEtaETSM%d%s", supermoduleID, triggerclass.Data()), eta, transverseenergy, weight);
    if(supermoduleID < 12)
      sector = 4 + int(supermoduleID/2); // EMCAL
    else
      sector = 13 + int((supermoduleID-12)/2);  // DCAL
    fHistos->FillTH2(Form("hEtaEnergySec%d%s", sector, triggerclass.Data()), eta, energy, weight);
    fHistos->FillTH2(Form("hEtaETSec%d%s", sector, triggerclass.Data()), eta, transverseenergy, weight);
  }
  if(hasTriggerPatch){
    fHistos->FillTH1(Form("hClusterEnergyFired%s", triggerclass.Data()), energy, weight);
    fHistos->FillTH1(Form("hClusterETFired%s", triggerclass.Data()), energy, weight);
    fHistos->FillTH2(Form("hEtaEnergyFired%s", triggerclass.Data()), eta, energy, weight);
    fHistos->FillTH2(Form("hEtaETFired%s", triggerclass.Data()), eta, energy, weight);
    if(supermoduleID >= 0){
      fHistos->FillTH2(Form("hClusterEnergyFiredSM%s", triggerclass.Data()), supermoduleID, energy, weight);
      fHistos->FillTH2(Form("hClusterETFiredSM%s", triggerclass.Data()), supermoduleID, transverseenergy, weight);
      fHistos->FillTH2(Form("hEtaEnergyFiredSM%d%s", supermoduleID, triggerclass.Data()), eta, energy,weight);
      fHistos->FillTH2(Form("hEtaETFiredSM%d%s", supermoduleID, triggerclass.Data()), eta, transverseenergy, weight);
      fHistos->FillTH2(Form("hEtaEnergyFiredSec%d%s", sector, triggerclass.Data()), eta, energy, weight);
      fHistos->FillTH2(Form("hEtaETFiredSec%d%s", sector, triggerclass.Data()), eta, transverseenergy, weight);
    }
  }
  Double_t encuts[5] = {1., 2., 5., 10., 20.};
  for(int ien = 0; ien < 5; ien++){
    if(energy > encuts[ien]){
      fHistos->FillTH2(Form("hEtaPhi%dG%s", static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi, weight);
      if(hasTriggerPatch){
        fHistos->FillTH2(Form("hEtaPhiFired%dG%s", static_cast<int>(encuts[ien]), triggerclass.Data()), eta, phi, weight);
      }
    }
  }
}

/**
 * Fill event-based histograms. Monitored are
 * - Number of events
 * - Centrality percentile (if available)
 * - z-position of the primary vertex
 * In case a downscaling correction is avaiable it is applied to all
 * histograms as a weight.
 * @param[in] triggerclass Name of the trigger class
 * @param[in] centrality Centrality percentile of the event
 * @param[in] vertexz z-position of the
 */
void AliAnalysisTaskEmcalClustersRef::FillEventHistograms(const TString &triggerclass, double centrality, double vertexz){
  Double_t weight = GetTriggerWeight(triggerclass);
  fHistos->FillTH1(Form("hEventCount%s", triggerclass.Data()), 1, weight);
  fHistos->FillTH1(Form("hEventCentrality%s", triggerclass.Data()), centrality, weight);
  fHistos->FillTH1(Form("hVertexZ%s", triggerclass.Data()), vertexz, weight);
}

/**
 * Check whether cluster is inside a trigger patch which has fired the trigger
 * @param[in] etaclust \f$ \eta \f$ of the cluster at center
 * @param[in] phiclust \f$ \phi \f$ of the cluster at center
 * @param[in] fTriggerPatches List of trigger patches which have fired the trigger
 * @return[in] True if the cluster can be correlated to a triggerpatch fired the trigger, false otherwise
 */
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

/**
 * Find all patches in an event which could have fired the trigger
 * Attention: This task groups into single shower triggers (L0, EG1, EG2) and jet triggers (EJ1 and EJ2).
 * Per convention the low threshold patch is selected. No energy cut should be applied in the trigger maker
 * @param triggerclass EMCAL trigger class firing
 * @param fTriggerPatches Trigger patches found in the event
 * @return List of patches which could have fired the trigger
 */
void AliAnalysisTaskEmcalClustersRef::FindPatchesForTrigger(TString triggerclass, const TClonesArray * fTriggerPatches, TList &foundtriggers) const {
  foundtriggers.Clear();
  if(!fTriggerPatches) return;
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
  for(TIter patchiter = TIter(fTriggerPatches).Begin(); patchiter != TIter::End(); ++patchiter){
    if(!IsOfflineSimplePatch(*patchiter)) continue;
    if(AliEmcalTriggerOfflineSelection::IsDCAL(myclass)){
      if(!SelectDCALPatch(*patchiter)) continue;
    } else {
      if(SelectDCALPatch(*patchiter)) continue;
    }
    if(AliEmcalTriggerOfflineSelection::IsSingleShower(myclass)){
      if(!SelectSingleShowerPatch(*patchiter)) continue;
    } else {
      if(!SelectJetPatch(*patchiter)) continue;
    }
    double threshold = fTriggerSelection ? fTriggerSelection->GetThresholdForTrigger(myclass) : -1;
    if(GetPatchEnergy(*patchiter) > threshold) foundtriggers.Add(*patchiter);
  }
}

/**
 * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
 * @param fTriggerPatches Trigger patches found by the trigger maker
 * @return String with EMCAL trigger decision
 */
TString AliAnalysisTaskEmcalClustersRef::GetFiredTriggerClassesFromPatches(const TClonesArray* fTriggerPatches) const {
  TString triggerstring = "";
  if(!fTriggerPatches) return triggerstring;
  Int_t nEJ1 = 0, nEJ2 = 0, nEG1 = 0, nEG2 = 0, nDJ1 = 0, nDJ2 = 0, nDG1 = 0, nDG2 = 0;
  double  minADC_J1 = 260.,
          minADC_J2 = 127.,
          minADC_G1 = 140.,
          minADC_G2 = 89.;
  for(TIter patchIter = TIter(fTriggerPatches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(patch->IsJetHighSimple() && patch->GetADCOfflineAmp() > minADC_J1){
      if(patch->IsDCalPHOS()) nDJ1++;
      else nEJ1++;
    }
    if(patch->IsJetLowSimple() && patch->GetADCOfflineAmp() > minADC_J2){
      if(patch->IsDCalPHOS()) nDJ2++;
      else nEJ2++;
    }
    if(patch->IsGammaHighSimple() && patch->GetADCOfflineAmp() > minADC_G1){
      if(patch->IsDCalPHOS()) nDG1++;
      else nEG1++;
    }
    if(patch->IsGammaLowSimple() && patch->GetADCOfflineAmp() > minADC_G2){
      if(patch->IsDCalPHOS()) nDG2++;
      else nEG2++;
    }
  }
  if(nEJ1) triggerstring += "EJ1";
  if(nEJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EJ2";
  }
  if(nEG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG1";
  }
  if(nEG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "EG2";
  }
  if(nDJ1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DJ1";
  }
  if(nDJ2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DJ2";
  }
  if(nDG1){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DG1";
  }
  if(nDG2){
    if(triggerstring.Length()) triggerstring += ",";
    triggerstring += "DG2";
  }
  return triggerstring;
}

void AliAnalysisTaskEmcalClustersRef::GetPatchBoundaries(TObject *o, Double_t *boundaries) const {
  AliEMCALTriggerPatchInfo *patch= dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  boundaries[0] = patch->GetEtaMin();
  boundaries[1] = patch->GetEtaMax();
  boundaries[2] = patch->GetPhiMin();
  boundaries[3] = patch->GetPhiMax();
}

bool AliAnalysisTaskEmcalClustersRef::IsOfflineSimplePatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->IsOfflineSimple();
}

bool AliAnalysisTaskEmcalClustersRef::SelectDCALPatch(TObject *o) const {
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->GetRowStart() >= 64;
}

bool AliAnalysisTaskEmcalClustersRef::SelectSingleShowerPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  return patch->IsGammaLowSimple();
}

bool AliAnalysisTaskEmcalClustersRef::SelectJetPatch(TObject *o) const{
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  if(!patch->IsOfflineSimple()) return false;
  return patch->IsJetLowSimple();
}

double AliAnalysisTaskEmcalClustersRef::GetPatchEnergy(TObject *o) const {
  double energy = 0.;
  AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo *>(o);
  energy = patch->GetPatchE();
  return energy;
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
