#include <THashList.h>
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
   fTriggerString("")
{

}

AliEmcalCellMonitorTask::AliEmcalCellMonitorTask(const char *name) :
   AliAnalysisTaskSE(name),
   fHistManager(nullptr),
   fGeometry(nullptr),
   fMinCellAmplitude(0),
   fRequestTrigger(AliVEvent::kAnyINT),
   fTriggerString("")
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
  fHistManager->CreateTH2("cellTimeOutlier", "Outlier time distribution per cell", 20001, -0.5, 20000.5, 100, 1e-6, 5e-5);
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
    if(amplitude <= 0) continue; // handle only cells with pos. amp
    fHistManager->FillTH2("cellAmplitude", cellNumber, amplitude);
    if(amplitude > fMinCellAmplitude){
      fHistManager->FillTH2("cellTime", cellNumber, celltime);
      if(celltime >= 1e-6) fHistManager->FillTH2("cellTimeOutlier", cellNumber, celltime);
    }

    // Get Cell index in eta-phi of sm
    fGeometry->GetCellIndex(cellNumber, sm, mod, mphi, meta);
    fGeometry->GetCellPhiEtaIndexInSModule(sm, mod, mphi, meta, iphi, ieta);

    fHistManager->FillTH2(Form("cellCountSM%d", sm), ieta, iphi);
    fHistManager->FillTH2(Form("cellAmpSM%d", sm), ieta, iphi, amplitude);
    fHistManager->FillTH2(Form("cellAmpTimeCorrSM%d", sm), celltime, amplitude);
  }
  PostData(1, fHistManager->GetListOfHistograms());
}
