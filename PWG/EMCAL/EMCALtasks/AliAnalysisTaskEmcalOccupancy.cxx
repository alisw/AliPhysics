#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <TLinearBinning.h>
#include <THistManager.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalOccupancy.h"
#include "AliClusterContainer.h"
#include "AliEMCALGeometry.h"
#include "AliMultSelection.h"
#include "AliVCaloCells.h"

ClassImp(AliAnalysisTaskEmcalOccupancy)

AliAnalysisTaskEmcalOccupancy::AliAnalysisTaskEmcalOccupancy() :
   AliAnalysisTaskEmcalLight(),
   fNameClusters(),
   fHistos(nullptr),
   fUseCentrality(false),
   fCellCounter(nullptr)
{

}

AliAnalysisTaskEmcalOccupancy::AliAnalysisTaskEmcalOccupancy(const char *name) :
  AliAnalysisTaskEmcalLight(name, true),
  fNameClusters(),
  fHistos(nullptr),
  fUseCentrality(false),
  fCellCounter(nullptr)
{
  SetNeedEmcalGeom(true);
  this->SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalOccupancy::~AliAnalysisTaskEmcalOccupancy() {
  if(fCellCounter) delete fCellCounter;
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalOccupancy::UserCreateOutputObjects() {
  AliAnalysisTaskEmcalLight::UserCreateOutputObjects();
  fHistos = new THistManager("histos");

  TLinearBinning multbinning(101, -0.5, 100.5),
                 cellbinning(20001, -0.5, 20000.5),
                 clusterbinning(2001, -0.5, 2000.5);
  const TBinning *histbinning[5] = {&multbinning, &cellbinning, &cellbinning, &clusterbinning, &clusterbinning};
  fHistos->CreateTHnSparse("EMCALOccupancy", "EMCAL occupancy", 5, histbinning);

  for(auto h : *fHistos) fOutput->Add(h);
  PostData(1, fOutput);
}

void AliAnalysisTaskEmcalOccupancy::ExecOnce() {
  AliAnalysisTaskEmcalLight::ExecOnce();
  if(!fLocalInitialized) return;

  // Create buffer for cells (would normally be done with a std::unique_ptr
  fCellCounter = new UChar_t[fGeom->GetNCells()];
}

Bool_t AliAnalysisTaskEmcalOccupancy::Run() {
  // Reset cell counter
  memset(fCellCounter, 0, sizeof(UChar_t) * fGeom->GetNCells());

  // Get the centrality percentile - steered by AliAnalysisTaskEmcalLight
  double centralityPercentile = (this->fForceBeamType == AliAnalysisTaskEmcalOccupancy::kAA) ? fCent : 0.;

  AliVCaloCells *emccells = InputEvent()->GetEMCALCells();
  int cellcounterRaw(0), cellcounterCorr(0);
  for(int icell = 0; icell < emccells->GetNumberOfCells(); icell++) {
    Short_t absID = emccells->GetCellNumber(icell);
    Double_t amplitude = emccells->GetAmplitude(icell);
    if(amplitude > DBL_EPSILON) {
      cellcounterRaw++;
      fCellCounter[absID]++;
    }
  }

  // Filter out all non-0 cells
  cellcounterCorr = std::count_if(fCellCounter, fCellCounter + fGeom->GetNCells(), [] (UChar_t count) -> bool { return count > 0; } );

  // check on cluster level - with and without exotics
  int clustercounterRaw(0), clustercounterCorr(0);
  AliClusterContainer *clusters = GetClusterContainer(fNameClusters.Data());
  clustercounterRaw = clusters->GetNClusters();
  const AliClusterIterableContainer &clusteriter = clusters->all();
  clustercounterCorr = std::count_if(clusteriter.begin(), clusteriter.end(), [](const AliVCluster *c) -> bool { return !c->GetIsExotic(); } );

  Double_t datapoint[5] = {centralityPercentile, static_cast<double>(cellcounterRaw), static_cast<double>(cellcounterCorr), static_cast<double>(clustercounterRaw), static_cast<double>(clustercounterCorr)};

  fHistos->FillTHnSparse("EMCALOccupancy", datapoint);

  return true;
}

AliAnalysisTaskEmcalOccupancy *AliAnalysisTaskEmcalOccupancy::AddOccupancyTask(const char *name) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    std::cout << "AddOccupancyTask: Analysis Manager not available - exiting ..." << std::endl;
    return nullptr;
  }

  AliAnalysisTaskEmcalOccupancy *occupancyTask = new AliAnalysisTaskEmcalOccupancy(name);
  mgr->AddTask(occupancyTask);

  TString outname = mgr->GetCommonFileName();
  outname += TString::Format(":EMCALOccupancy");

  mgr->ConnectInput(occupancyTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(occupancyTask, 1, mgr->CreateContainer("EMCALOccupancyHists", TList::Class(), AliAnalysisManager::kOutputContainer, outname));

  return occupancyTask;
}
