///
/// \file AliAnalysisTaskFemtoNu.cxx
///

#include "TROOT.h"
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TInterpreter.h"

#include "AliESDEvent.h"

#include "AliFemtoAnalysis.h"
#include "AliAnalysisTaskFemtoNu.h"
#include "AliVHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliFemtoEventReaderAODMultSelection.h"


const unsigned AliAnalysisTaskFemtoNu::RESULT_STORAGE_OUTPUT_SLOT = 1;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskFemtoNu);
/// \endcond


AliAnalysisTaskFemtoNu::AliAnalysisTaskFemtoNu()
  : AliAnalysisTaskFemto()
  , fStorage(nullptr)
  , fUseCustomExec(false)
{
}

// Default name for the setup macro of femto analysis
AliAnalysisTaskFemtoNu::AliAnalysisTaskFemtoNu(const TString &name,
                                               const TString &mac,
                                               const Bool_t verbose)
  : AliAnalysisTaskFemtoNu(name, mac, "", verbose)
{
}

// This function MUST be defined in the separate file !!!
// extern AliFemtoManager *ConfigFemtoAnalysis();

//________________________________________________________________________
AliAnalysisTaskFemtoNu::AliAnalysisTaskFemtoNu(const TString &name,
                                               const TString &mac,
                                               const TString &par,
                                               const Bool_t verbose)
  : AliAnalysisTaskFemto(name, mac, par, verbose)
  , fStorage(nullptr)
  , fUseCustomExec(false)
{
  DefineOutput(RESULT_STORAGE_OUTPUT_SLOT, AliFemtoResultStorage::Class());
}

AliAnalysisTaskFemtoNu::~AliAnalysisTaskFemtoNu()
{
}

void
AliAnalysisTaskFemtoNu::ConnectInputData(Option_t *opt)
{
  AliAnalysisTaskFemto::ConnectInputData(opt);

  if (auto *rdr = dynamic_cast<AliFemtoEventReaderAODChain *>(fReader)) {
    rdr->SetAODSource(fAOD);
  }
}


void
AliAnalysisTaskFemtoNu::CreateOutputObjects()
{
  if (fVerbose) {
    AliInfo("Creating Nouveau output objects");
  }

  gROOT->LoadMacro(fConfigMacro);

  const TString cmd = Form("ConfigFemtoAnalysis(%s)", fConfigParams ? fConfigParams.Data() : "");
  if (fVerbose) {
    AliInfo(Form("Calling `%s`", cmd.Data()));
  }

  auto femto_mgr = reinterpret_cast<AliFemtoManager*>(gInterpreter->ProcessLine(cmd));
  if (femto_mgr == nullptr) {
    AliError(Form("Could not create AliFemtoManager from %s::%s",
                  fConfigMacro.Data(), cmd.Data()));
    return;
  }

  if (femto_mgr->AnalysisCollection()->size() == 0) {
    AliError("FemtoAnalysisManager has no FemtoAnalyses");
    delete femto_mgr;
    return;
  }

  SetFemtoManager(femto_mgr);

  fStorage = new AliFemtoResultStorage(fName, *femto_mgr);

  PostData(RESULT_STORAGE_OUTPUT_SLOT, fStorage);
}

void
AliAnalysisTaskFemtoNu::Exec(Option_t *ex)
{
  // just use the standard Exec
  if (!fUseCustomExec) {
    return AliAnalysisTaskFemto::Exec(ex);
  }

  auto *mgr = AliAnalysisManager::GetAnalysisManager();
  auto *input_handler = mgr->GetInputEventHandler();
  // auto *mc_input_handler = mgr->GetMCtruthEventHandler();

  // Task making a femtoscopic analysis.
  if (fOfflineTriggerMask) {
    Bool_t is_selected = input_handler->IsEventSelected() & fOfflineTriggerMask;
    if (!is_selected) {
      if (fVerbose) {
        cout << "AliAnalysisTaskFemto: is not selected "
             << input_handler->IsEventSelected() << " != " << fOfflineTriggerMask << "\n";
      }
      return;
    }
  }

  fManager->ProcessEvent();
}
