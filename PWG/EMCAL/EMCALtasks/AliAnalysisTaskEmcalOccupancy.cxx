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
}

AliAnalysisTaskEmcalOccupancy::~AliAnalysisTaskEmcalOccupancy() {
  if(fCellCounter) delete fCellCounter;
}

void AliAnalysisTaskEmcalOccupancy::UserCreateOutputObjects() {
  fHistos = new THistManager("histos");

  for(auto h : *fHistos) fOutput->Add(h);
  PostData(1, fOutput);
}

void AliAnalysisTaskEmcalOccupancy::ExecOnce() {
  AliAnalysisTaskEmcalLight::ExecOnce();
  if(!fLocalInitialized) return;

  TLinearBinning multbinning(1000, 0., 100.),
                 cellbinning(25001, -0.5, 25000.5),
                 clusterbinning(1001, -0.5, 101.5);

  const TBinning *histbinning[5] = {&multbinning, &cellbinning, &cellbinning, &clusterbinning, &clusterbinning};
  fHistos->CreateTHnSparse("EMCALOccupancy", "EMCAL occupancy", 5, histbinning);

  // Create buffer for cells (would normally be done with a std::unique_ptr
  fCellCounter = new UChar_t[fGeom->GetNCells()];
}

Bool_t AliAnalysisTaskEmcalOccupancy::Run() {
  // Reset cell counter
  memset(fCellCounter, 0, sizeof(fGeom->GetNCells()));

  // Get the centrality percentile
  double centralityPercentile = 0;
  if(fUseCentrality) {
    AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if(mult) centralityPercentile = mult->GetMultiplicityPercentile("V0M");
  }

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

AliAnalysisTaskEmcalOccupancy *AddOccupancyTask(const char *name) {
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
