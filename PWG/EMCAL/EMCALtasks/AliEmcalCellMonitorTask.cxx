#include <algorithm>
#include <map>
#include <vector>
#include <TArrayD.h>
#include <THashList.h>
#include <TLinearBinning.h>
#include <THistManager.h>

#include "AliEMCALGeometry.h"
#include "AliEmcalCellMonitorTask.h"
#include "AliInputEventHandler.h"
#include "AliVCaloCells.h"
#include "AliVEvent.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCellMonitorTask)
/// \endcond

AliEmcalCellMonitorTask::AliEmcalCellMonitorTask() :
   AliAnalysisTaskSE(),
   fHistManager(nullptr),
   fGeometry(nullptr),
   fMinCellAmplitude(0),
   fRequestTrigger(AliVEvent::kAnyINT),
   fTriggerString(""),
   fNumberOfCells(12288),
   fMaskedCells()
{

}

AliEmcalCellMonitorTask::AliEmcalCellMonitorTask(const char *name) :
   AliAnalysisTaskSE(name),
   fHistManager(nullptr),
   fGeometry(nullptr),
   fMinCellAmplitude(0),
   fRequestTrigger(AliVEvent::kAnyINT),
   fTriggerString(""),
   fNumberOfCells(12288),
   fMaskedCells()
{
  DefineOutput(1, TList::Class());
}

AliEmcalCellMonitorTask::~AliEmcalCellMonitorTask() {
  if(fGeometry) delete fGeometry;
}

void AliEmcalCellMonitorTask::UserCreateOutputObjects(){
  fHistManager = new THistManager("EMCALCellMonitor");

  fHistManager->CreateTH1("cellFrequency", "Frequency of cell firing", TLinearBinning(fNumberOfCells, -0.5, fNumberOfCells - 0.5));
  fHistManager->CreateTH2("cellAmplitude", "Energy distribution per cell", TLinearBinning(fNumberOfCells, -0.5, fNumberOfCells - 0.5), AliEmcalCellMonitorAmplitudeBinning());
  fHistManager->CreateTH2("cellAmplitudeCut", "Energy distribution per cell (after energy cut)", TLinearBinning(fNumberOfCells, -0.5, fNumberOfCells - 0.5), AliEmcalCellMonitorAmplitudeBinning());
  fHistManager->CreateTH2("cellTime", "Time distribution per cell", fNumberOfCells, -0.5, fNumberOfCells - 0.5, 200, -1e-6, 1e-6);
  fHistManager->CreateTH2("cellTimeOutlier", "Outlier time distribution per cell", fNumberOfCells, -0.5, fNumberOfCells - 0.5, 100, 1e-6, 5e-5);
  fHistManager->CreateTH2("cellTimeMain", "Time distribution per cell for the main bunch", fNumberOfCells, -0.5, fNumberOfCells - 0.5, 150, -50e-9, 100e-9);
  for(int ism = 0; ism < 20; ++ism){
    fHistManager->CreateTH2(Form("cellAmpSM%d", ism), Form("Integrated cell amplitudes for SM %d; col; row", ism), 48, -0.5, 47.5, 24, -0.5, 23.5);
    fHistManager->CreateTH2(Form("cellCountSM%d", ism), Form("Count rate per cell for SM %d; col; row", ism), 48, -0.5, 47.5, 24, -0.5, 23.5);
  }

  for(int ism = 0; ism < 20; ++ism){
    fHistManager->CreateTH2(Form("cellAmpTimeCorrSM%d", ism), Form("Correlation between cell amplitude and time in Supermodule %d", ism), 1000, -5e-7, 5e-7, 1000, 0., 100.);
  }

  PostData(1, fHistManager->GetListOfHistograms());
}

void AliEmcalCellMonitorTask::UserExec(Option_t *){
  if(!fGeometry) fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(fInputEvent->GetRunNumber());

  // Check trigger
  if(!(fInputHandler->IsEventSelected() & fRequestTrigger)) return;
  if(fTriggerString.Length()){
    if(!TString(InputEvent()->GetFiredTriggerClasses()).Contains(fTriggerString)) return;
  }

  AliVCaloCells *emcalcells = fInputEvent->GetEMCALCells();

  // input data
  Short_t cellNumber;
  Double_t amplitude, celltime, efrac;
  Int_t mclabel;

  Int_t sm, mod, meta, mphi, ieta, iphi;
  for(int icell = 0; icell < emcalcells->GetNumberOfCells(); icell++){
    emcalcells->GetCell(icell, cellNumber, amplitude, celltime, mclabel, efrac);
    if(IsCellMasked(cellNumber)) continue;
    fHistManager->FillTH2("cellAmplitude", cellNumber, amplitude);
    if(amplitude < fMinCellAmplitude) continue;
    fHistManager->FillTH1("cellAmplitudeCut", cellNumber, amplitude);
    fHistManager->FillTH1("cellFrequency", cellNumber);
    fHistManager->FillTH2("cellTime", cellNumber, celltime);
    if(celltime >= 1e-6) fHistManager->FillTH2("cellTimeOutlier", cellNumber, celltime);
    if(celltime > -5e-8 && celltime < 1e-7) fHistManager->FillTH2("cellTimeMain", cellNumber, celltime);

    // Get Cell index in eta-phi of sm
    fGeometry->GetCellIndex(cellNumber, sm, mod, mphi, meta);
    fGeometry->GetCellPhiEtaIndexInSModule(sm, mod, mphi, meta, iphi, ieta);

    fHistManager->FillTH2(Form("cellCountSM%d", sm), ieta, iphi);
    fHistManager->FillTH2(Form("cellAmpSM%d", sm), ieta, iphi, amplitude);
    fHistManager->FillTH2(Form("cellAmpTimeCorrSM%d", sm), celltime, amplitude);
  }
  PostData(1, fHistManager->GetListOfHistograms());
}

void AliEmcalCellMonitorTask::SetBadCell(Int_t cellId){
  if(std::find(fMaskedCells.begin(), fMaskedCells.end(), cellId) != fMaskedCells.end()) return;
  fMaskedCells.push_back(cellId);
}

bool AliEmcalCellMonitorTask::IsCellMasked(Int_t cellId) const {
  return (std::find(fMaskedCells.begin(), fMaskedCells.end(), cellId) != fMaskedCells.end());
}

AliEmcalCellMonitorTask::AliEmcalCellMonitorAmplitudeBinning::AliEmcalCellMonitorAmplitudeBinning():
    TCustomBinning()
{
  SetMinimum(0);
  AddStep(1., 0.1);
  AddStep(10., 0.5);
  AddStep(20., 1.);
  AddStep(50., 2.);
  AddStep(100., 5.);
  AddStep(200., 10.);
}
