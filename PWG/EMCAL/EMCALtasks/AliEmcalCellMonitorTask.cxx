#include <THashList.h>
#include <THistManager.h>

#include "AliEMCALGeometry.h"
#include "AliEmcalCellMonitorTask.h"
#include "AliVCaloCells.h"
#include "AliVEvent.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCellMonitorTask)
/// \endcond

AliEmcalCellMonitorTask::AliEmcalCellMonitorTask() :
   AliAnalysisTaskSE(),
   fHistManager(nullptr),
   fGeometry(nullptr)
{

}

AliEmcalCellMonitorTask::AliEmcalCellMonitorTask(const char *name) :
   AliAnalysisTaskSE(name),
   fHistManager(nullptr),
   fGeometry(nullptr)
{
  DefineOutput(1, TList::Class());
}

AliEmcalCellMonitorTask::~AliEmcalCellMonitorTask() {
  if(fGeometry) delete fGeometry;
}

void AliEmcalCellMonitorTask::UserCreateOutputObjects(){
  fHistManager = new THistManager("EMCALCellMonitor");

  fHistManager->CreateTH2("cellAmplitude", "Energy distribution per cell", 20001, -0.5, 20000.5, 1000, 0., 100);
  fHistManager->CreateTH2("cellTime", "Time distribution per cell", 20001, -0.5, 20000.5, 1000, -1e-6, 1e-6);
  for(int ism = 0; ism < 20; ++ism){
    fHistManager->CreateTH2(Form("cellAmpSM%d", ism), Form("Integrated cell amplitudes for SM %d; col; row", ism), 48, -0.5, 47.5, 24, -0.5, 23.5);
    fHistManager->CreateTH2(Form("cellCountSM%d", ism), Form("Count rate per cell for SM %d; col; row", ism), 48, -0.5, 47.5, 24, -0.5, 23.5);
  }

  PostData(1, fHistManager->GetListOfHistograms());
}

void AliEmcalCellMonitorTask::UserExec(Option_t *){
  if(!fGeometry) fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(fInputEvent->GetRunNumber());

  AliVCaloCells *emcalcells = fInputEvent->GetEMCALCells();

  // input data
  Short_t cellNumber;
  Double_t amplitude, celltime, efrac;
  Int_t mclabel;

  Int_t sm, mod, meta, mphi, ieta, iphi;
  for(int icell = 0; icell < emcalcells->GetNumberOfCells(); icell++){
    emcalcells->GetCell(icell, cellNumber, amplitude, celltime, mclabel, efrac);
    if(amplitude <= 0) continue; // handle only cells with pos. amp
    fHistManager->FillTH2("cellAmplitude", cellNumber, amplitude);
    fHistManager->FillTH2("cellTime", cellNumber, celltime);

    // Get Cell index in eta-phi of sm
    fGeometry->GetCellIndex(cellNumber, sm, mod, mphi, meta);
    fGeometry->GetCellPhiEtaIndexInSModule(sm, mod, mphi, meta, iphi, ieta);

    fHistManager->FillTH2(Form("cellCountSM%d", sm), ieta, iphi);
    fHistManager->FillTH2(Form("cellAmpSM%d", sm), ieta, iphi, amplitude);
  }
  PostData(1, fHistManager->GetListOfHistograms());
}
