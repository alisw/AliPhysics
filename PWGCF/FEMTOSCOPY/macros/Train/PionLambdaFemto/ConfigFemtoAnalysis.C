///
/// \file PionLambdaFemto/ConfigFemtoAnalysis.C
///
/// \brief The configuration macro which sets up the pion-lambda analysis
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

#if !defined(__CINT__) || defined(__MAKECINT_)

#include "AliFemtoAnalysisPionLambda.h"

#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"

#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorV0.h"

#endif

typedef AliFemtoAnalysisPionLambda AFAPL;

void BuildConfiguration(const TString&,
                        AliFemtoAnalysisPionLambda::AnalysisParams&,
                        AliFemtoAnalysisPionLambda::CutParams&);

AliFemtoManager*
ConfigFemtoAnalysis(const TString& param_str = "")
{
  const double PionMass = 0.13956995,
             LambdaMass = 1.115683;

  AliFemtoManager *manager = new AliFemtoManager();


  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAODChain();
  rdr->SetFilterBit(7);
  rdr->SetEPVZERO(kTRUE);
  rdr->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  rdr->SetCentralityFlattening(kFALSE);
  rdr->SetReadV0(1);
  rdr->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
  manager->SetEventReader(rdr);

  // Get the default configurations
  AFAPL::AnalysisParams analysis_config = AFAPL::DefaultConfig();
  AFAPL::CutParams cut_config = AFAPL::DefaultCutConfig();

  // Read string and update configurations
  BuildConfiguration(param_str, analysis_config, cut_config);

  AliFemtoAnalysisPionLambda *analysis = new AliFemtoAnalysisPionLambda("PionLambda",
                                                                        analysis_config,
                                                                        cut_config);
    analysis->SetNumEventsToMix(10);
    analysis->SetMinSizePartCollection(1);
    analysis->SetVerboseMode(kFALSE);

  analysis->EventCut()->AddCutMonitor(new AliFemtoCutMonitorEventMult("EVPass"),
                                      new AliFemtoCutMonitorEventMult("EVFail"));

  analysis->GetPionCut()->AddCutMonitor(new AliFemtoCutMonitorParticleYPt("PionPass", PionMass),
                                        new AliFemtoCutMonitorParticleYPt("PionFail", PionMass));

  analysis->GetLambdaCut()->AddCutMonitor(new AliFemtoCutMonitorV0("LambdaPass"),
                                          new AliFemtoCutMonitorV0("LambdaFail"));

  manager->AddAnalysis(analysis);

  return manager;
}

void BuildConfiguration(const TString &text,
                        AliFemtoAnalysisPionLambda::AnalysisParams &a,
                        AliFemtoAnalysisPionLambda::CutParams &cut)
{
  std::cout << "[BuildAnalysisFromString]\n";
  std::cout << "   '" << text << "'\n";

  const TString analysis_varname = "a",
                     cut_varname = "cut";


  TObjArray* lines = text.Tokenize("\n;");

  TIter next_line(lines);
  TObject *line_obj = NULL;

  while (line_obj = next_line()) {

    const TString line = ((TObjString*)line_obj)->String().ReplaceAll(".", "_")
                                                          .Strip(TString::kBoth, ' ');
                                                          

    TString cmd("");

    switch (line[0]) {
    case '$':
      cmd = cut_varname + "." + line(1, line.Length() - 1);
      break;

    case '@':
      cmd = analysis_varname + "." + line(1, line.Length() - 1);
      break;

    default:
      continue;
    }

    cmd += ";";

    cout << "CMD: `" << cmd << "`\n";

    gROOT->ProcessLineFast(cmd);
  }
}
