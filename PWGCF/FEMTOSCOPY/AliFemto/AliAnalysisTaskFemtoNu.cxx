///
/// \file AliAnalysisTaskFemtoNu.cxx
///

#include "AliAnalysisTaskFemtoNu.h"

#include "AliFemtoAnalysis.h"
#include "AliFemtoEventReaderAODMultSelection.h"

#include <TROOT.h>
#include <TChain.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TFile.h>
#include <TInterpreter.h>

#include <iostream>


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
AliAnalysisTaskFemtoNu::SetupContainers(const TString &outputfile)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  const TString dest = outputfile.IsWhitespace() ? mgr->GetCommonFileName() : outputfile.Data();

  auto *out_container = mgr->CreateContainer(fName,
                                             AliFemtoResultStorage::Class(),
                                             AliAnalysisManager::kOutputContainer,
                                             dest);
  auto *in_container = mgr->GetCommonInputContainer();

  mgr->ConnectInput(this, 0, in_container);
  mgr->ConnectOutput(this, RESULT_STORAGE_OUTPUT_SLOT, out_container);

  // this is required to silence a runtime warning about no container
  // for output slot 0 (defined by AliAnalysiTaskFemto)
  //
  // when this class gets re-written without TList output remove this
  // line
  //
  auto *list_container = mgr->CreateContainer(fName + "_list",
                                              TList::Class(),
                                              AliAnalysisManager::kOutputContainer,
                                              dest);
  mgr->ConnectOutput(this, 0, list_container);
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

  auto tlist = new TList();
  tlist->SetOwner();
  PostData(0, tlist);
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
  if (fOfflineTriggerMask &&
      (input_handler->IsEventSelected() & fOfflineTriggerMask) == 0) {
    if (fVerbose) {
      std::cout << "AliAnalysisTaskFemto: is not selected "
                << input_handler->IsEventSelected() << " != "
                << fOfflineTriggerMask << "\n";
    }
    return;
  }

  static_cast<AliFemtoEventReaderAODChain*>(fReader)->SetAODSource(fAOD);

  fManager->ProcessEvent();
}


AliAnalysisTaskFemtoNu::MacroCfg::MacroCfg()
  : TNamed("cfg", "MacroCfg")
  , macro("%%/ConfigNuFemtoAnalysisR6.C")
  , auto_directory("$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/Train/PionPionFemto")
  , output_filename(AliAnalysisManager::GetAnalysisManager()->GetCommonFileName())
  , output_container("PWG2FEMTO")
  , subwagon_array("")
  , subwagon_type("centrality")
{}


#include <regex>


AliAnalysisTask*
AliAnalysisTaskFemtoNu::BuildForMacro(TString container,
                                      TString configuration,
                                      TString params,
                                      TString subwagon_suffix)
{
  std::cout << "\n\n============== AliAnalysisTaskFemtoNu::BuildForMacro (ROOT6 Only) ===============\n";

  // Get the global analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    std::cerr << "E-AliAnalysisTaskFemtoNu::BuildForMacro:  "
              << "Could not get the global AliAnalysisManager.\n";
    return nullptr;
  }

  MacroCfg cfg;

  bool verbose = kFALSE;

  TObjArray* lines = configuration.Tokenize("\n;");

  TIter next_line(lines);
  TObject *line_obj = nullptr;

  gDirectory->Add(&cfg);

  while ((line_obj = next_line())) {
    TObjString *s = static_cast<TObjString*>(line_obj);
    TString cmd = s->String().Strip(TString::kBoth, ' ');
    if (cmd.IsWhitespace()) {
      continue;
    }

    cmd.ReplaceAll("'", '"');

    std::cmatch cm;
    {
      std::regex re(" *subwagon_type *= *\"([^\"]+)\"");
      if (std::regex_match(cmd.Data(), cm, re)) {
        cfg.subwagon_type = cm[1].str().c_str();
        continue;
      }
    }

    {
      std::regex re(" *subwagon_array *= *\"([^\"]+)\"");
      if (std::regex_match(cmd.Data(), cm, re)) {
        cfg.subwagon_array = cm[1].str().c_str();
        std::cout << "ARRAY: " << cfg.subwagon_array << "\n";
        continue;
      }
    }

    // cmd = "static_cast<AliAnalysisTaskFemtoNu::MacroCfg*>(cfg)->" + cmd + ";";
    // std::cout << "running `" << cmd  << "`\n";
    // gROOT->ProcessLine(cmd);
  }

  gDirectory->Remove(&cfg);

  // Replace %% with this directory for convenience
  cfg.macro.ReplaceAll("%%", cfg.auto_directory);

  // Dealing with subwagons
  if (!subwagon_suffix.IsWhitespace()) {
    Int_t index = subwagon_suffix.Atoi();
    TObjArray *values = cfg.subwagon_array.Tokenize(",");
    TIter next_value(values);
    if (values->GetEntries() < index) {
      std::cerr << "Could not use subwagon-index " << index << " in subwagon_array of only " << values->GetEntries() << " entries\n";
      return nullptr;
    }
    for (int i=0; i<index; ++i) {
      next_value();
    }
    TString ss = ((TObjString*)next_value())->String();
    params += ";" + cfg.subwagon_type + " = " + ss;

    container += "_" + subwagon_suffix;
  }

  std::cout << "[AddTaskPionPion]\n"
               "   container: " << container << "\n"
               "   output: '" << cfg.output_filename << "'\n"
               "   macro: '" << cfg.macro << "'\n"
               "   params: '" << params << "'\n";

  // The analysis config macro for PionPionFemto accepts a single string
  // argument, which it interprets.
  // This line escapes some escapable characters (backslash, newline, tab)
  // and wraps that string in double quotes, ensuring that the interpreter
  // reads a string when passing to the macro.
  const TString analysis_params = '"' + params.ReplaceAll("\\", "\\\\")
                                              .ReplaceAll("\n", "\\n")
                                              .ReplaceAll("\"", "\\\"")
                                              .ReplaceAll("\t", "\\t") + '"';

  AliAnalysisTaskFemto *femtotask = new AliAnalysisTaskFemtoNu(
    container,
    cfg.macro,
    analysis_params,
    verbose
  );

  mgr->AddTask(femtotask);

  const TString outputfile = cfg.GetFilename();
  auto *out_container = mgr->CreateContainer(container,
                                             AliFemtoResultStorage::Class(),
                                             AliAnalysisManager::kOutputContainer,
                                             outputfile);

  mgr->ConnectInput(femtotask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(femtotask, AliAnalysisTaskFemtoNu::RESULT_STORAGE_OUTPUT_SLOT, out_container);

  // quiet a warning
  out_container = mgr->CreateContainer(container + "_list",
                                       TList::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputfile);
  mgr->ConnectOutput(femtotask, 0, out_container);

  std::cout << "============== AddNuTaskPionPion : Done ===============\n\n";
  return femtotask;
}


AliAnalysisTaskFemtoNu::MacroCfg::~MacroCfg()
{}
