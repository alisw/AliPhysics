#include <array>
#include <string>
#include <vector>

#include <THistManager.h>
#include <THashList.h>

#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskEGAMonitor.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerTypes.h"
#include "AliInputEventHandler.h"
#include "AliVVertex.h"
#include "AliVCaloTrigger.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEGAMonitor)
/// \endcond

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskEGAMonitor::AliAnalysisTaskEGAMonitor() :
    AliAnalysisTaskSE(),
    fHistos(nullptr),
    fAnalysisUtils(nullptr),
    fGeom(nullptr),
    fLocalInitialized(false)
{

}

AliAnalysisTaskEGAMonitor::AliAnalysisTaskEGAMonitor(const char *name) :
    AliAnalysisTaskSE(name),
    fHistos(nullptr),
    fAnalysisUtils(nullptr),
    fGeom(nullptr),
    fLocalInitialized(false)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEGAMonitor::~AliAnalysisTaskEGAMonitor() {
  if(fHistos) delete fHistos;
  if(fGeom) delete fGeom;
}

void AliAnalysisTaskEGAMonitor::UserCreateOutputObjects(){
  fAnalysisUtils = new AliAnalysisUtils;

  fHistos = new THistManager("EGAhistos");
  std::array<std::string, 4> triggers = {"EG1", "EG2", "DG1", "DG2"};
  for(const auto &t : triggers){
    fHistos->CreateTH2(Form("hColRowG1%s", t.c_str()), Form("Col-Row distribution of online G1 patches for trigger %s", t.c_str()), 48, -0.5, 47.5, 104, -0.5, 103.5);
    fHistos->CreateTH2(Form("hColRowG2%s", t.c_str()), Form("Col-Row distribution of online G2 patches for trigger %s", t.c_str()), 48, -0.5, 47.5, 104, -0.5, 103.5);
  }

  PostData(1, fHistos);
}

void AliAnalysisTaskEGAMonitor::UserExec(Option_t *){
  if(!fLocalInitialized){
    ExecOnce();
    fLocalInitialized = true;
  }

  if(!fAnalysisUtils->IsVertexSelected2013pA(InputEvent())) return;
  if(fAnalysisUtils->IsPileUpEvent(InputEvent())) return;
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA)) return;

  std::vector<std::string> triggers;
  if(InputEvent()->GetFiredTriggerClasses().Contains("EG1")) triggers.push_back("EG1");
  if(InputEvent()->GetFiredTriggerClasses().Contains("EG2")) triggers.push_back("EG2");
  if(InputEvent()->GetFiredTriggerClasses().Contains("DG1")) triggers.push_back("DG1");
  if(InputEvent()->GetFiredTriggerClasses().Contains("DG2")) triggers.push_back("DG2");

  AliVCaloTrigger *emctrigraw = InputEvent()->GetCaloTrigger("EMCAL");

  emctrigraw->Reset();
  Int_t col(-1), row(-1), triggerbits(0);
  while(emctrigraw->Next()){
    emctrigraw->GetTriggerBits(triggerbits);
    if(!(triggerbits & (BIT(kL1GammaHigh) | BIT(kL1GammaLow)))) continue;

    emctrigraw->GetPosition(col, row);
    if(triggerbits & BIT(kL1GammaHigh)){
      for(const auto &t : triggers) fHistos->FillTH2(Form("hColRowG1%s", t.c_str()), col, row);
    }
    if(triggerbits & BIT(kL1GammaLow)){
      for(const auto &t : triggers) fHistos->FillTH2(Form("hColRowG2%s", t.c_str()), col, row);
    }
  }

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEGAMonitor::ExecOnce(){
  if(!fGeom) fGeom = AliEMCALGeometry::GetInstanceFromRunNumber(InputEvent()->GetRunNumber());
}


} /* namespace EMCalTriggerPtAnalysis */
