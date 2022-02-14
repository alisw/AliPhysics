///
/// \file PionLambdaFemto/ConfigFemtoAnalysis.C
///
/// \brief The configuration macro which sets up the pion-lambda analysis
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

// #include <set>

#if !defined(__CINT__) || defined(__MAKECINT_)

#include "AliFemtoAnalysisPionLambda.h"

#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"

#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorParticleYPt.h"
#include "AliFemtoCutMonitorV0.h"

#include "AliFemtoAvgSepCorrFctn.h"
#include "AliFemtoCorrFctnKStar.h"
#include "AliFemtoModelCorrFctnKStar.h"

#include <TROOT.h>

#endif

typedef AliFemtoAnalysisPionLambda AFAPL;

/**
 * Parameter container for this ConfigFemtoAnalysis macro
 */
struct MacroParams {
  MacroParams(): analysis_configurations(){};
//   std::set<short> analysis_configurations;

  /// A bitmask of different pion/lambda configurations
  /// bit 0 : Pi+Lambda
  ///     1 : Pi+AntiLambda
  ///     2 : Pi-Lambda
  ///     3 : Pi+AntiLambda
  unsigned short analysis_configurations;
};

void BuildConfiguration(
  const TString&,
  AliFemtoAnalysisPionLambda::AnalysisParams&,
  AliFemtoAnalysisPionLambda::CutParams&,
  MacroParams&
);

AliFemtoManager*
ConfigFemtoAnalysis(const TString& param_str = "")
{
  std::cout << "[ConfigFemtoAnalysis]\n";

  const double PionMass = 0.13956995,
             LambdaMass = 1.115683;

  // Get the default configurations
  AFAPL::AnalysisParams analysis_config = AFAPL::DefaultConfig();
  AFAPL::CutParams cut_config = AFAPL::DefaultCutConfig();
  MacroParams macro_config;

  // Read parameter string and update configurations
  BuildConfiguration(param_str, analysis_config, cut_config, macro_config);

  // Begin to build the manager and analyses
  AliFemtoManager *manager = new AliFemtoManager();

  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAODChain();
  rdr->SetFilterBit(7);
  rdr->SetEPVZERO(kTRUE);
  rdr->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
  rdr->SetCentralityFlattening(kFALSE);
  rdr->SetReadV0(1);
  rdr->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
  rdr->SetReadMC(analysis_config.is_mc_analysis);
  manager->SetEventReader(rdr);

//   if (macro_config.analysis_configurations.empty()) {
//     // insert pi+ lambda code
//     macro_config.analysis_configurations.insert(3);
//   }

  // if no configuration was set, set the first bit
  if (macro_config.analysis_configurations == 0) {
    macro_config.analysis_configurations = 0x01;
  }

//   std::set<short>::iterator config;
//   for (config = macro_config.analysis_configurations.begin();
//        config != macro_config.analysis_configurations.end();
//        ++config) {

  for (int bit = 0; bit < 4; ++bit) {

    const int config_code = (1 << bit);

    // if bit is not set - skip it
    if (!(macro_config.analysis_configurations & config_code)) {
      continue;
    }

    // Update the particle types
    // This currently uses the fact that there are only two states : lam/antilam plus/minus
    // that happen to correlate with both 0/1 integers and true/false booleans. We take
    // advantage of this because ROOT doesn't like enums in macros. If this becomes compiled
    // it should be re-written to set to real enum values.
    // If first bit is set, use antilambda (1) else lambda (0)
    analysis_config.lambda_type = (bit & 1) == 1 ? AFAPL::kAntiLambda : AFAPL::kLambda;

    // If second bit is set, use pi- (1) else pi+ (0)
    analysis_config.pion_type = (bit & 2) == 2 ? AFAPL::kPiMinus : AFAPL::kPiPlus;

//     analysis_config.pion_type = (*config & 2) == 2 ? AFAPL::kPiPlus : AFAPL::kPiMinus;
//     analysis_config.lambda_type = (*config & 1) == 1 ? AFAPL::kLambda : AliFemtoAnalysisPionLambda::kAntiLambda;

    TString analysis_name;
    analysis_name += (analysis_config.pion_type == AFAPL::kPiPlus) ? "PiPlus" : "PiMinus";
    analysis_name += (analysis_config.lambda_type == AFAPL::kLambda) ? "Lam" : "ALam";

    AliFemtoAnalysisPionLambda *analysis = new AliFemtoAnalysisPionLambda(analysis_name,
                                                                          analysis_config,
                                                                          cut_config);

    analysis->AddStanardCutMonitors();

/*
    analysis->EventCut()->AddCutMonitor(new AliFemtoCutMonitorEventMult("EVPass"),
                                        new AliFemtoCutMonitorEventMult("EVFail"));

    analysis->GetPionCut()->AddCutMonitor(new AliFemtoCutMonitorParticleYPt("PionPass", PionMass),
                                          new AliFemtoCutMonitorParticleYPt("PionFail", PionMass));

    analysis->GetLambdaCut()->AddCutMonitor(new AliFemtoCutMonitorV0("LambdaPass"),
                                            new AliFemtoCutMonitorV0("LambdaFail"));
*/

    analysis->AddCorrFctn(new AliFemtoCorrFctnKStar("CF", 660, 0, 2.0));

    AliFemtoAvgSepCorrFctn *avgsep_cf = new AliFemtoAvgSepCorrFctn("avgsep_cf", 240, 0.0, 40.0);
    avgsep_cf->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
    analysis->AddCorrFctn(avgsep_cf);

    if (analysis_config.is_mc_analysis) {

      AliFemtoModelCorrFctnKStar *cf = new AliFemtoModelCorrFctnKStar("mc_cf", 660, 0, 2.0);
      cf->SetPairType(AliFemtoAvgSepCorrFctn::kTrackV0);
      const Int_t lam_pdg = analysis_config.lambda_type == 0 ? 3122 : -3122,
                   pi_pdg = analysis_config.pion_type == 0 ? 211 : -211;
      cf->SetExpectedPDGCodes(lam_pdg, pi_pdg);

      analysis->AddCorrFctn(cf);

    }

    manager->AddAnalysis(analysis);
  }

  return manager;
}

void BuildConfiguration(const TString &text,
                        AliFemtoAnalysisPionLambda::AnalysisParams &a,
                        AliFemtoAnalysisPionLambda::CutParams &cut,
                        MacroParams &mac
                        )
{
  std::cout << "I-BuildConfiguration:" << TBase64::Encode(text) << " \n";
  //   std::cout << "   '" << text << "'\n";

  const TString analysis_varname = "a",
                     cut_varname = "cut",
                   macro_varname = "mac";


  TObjArray* lines = text.Tokenize("\n;");

  TIter next_line(lines);
  TObject *line_obj = NULL;

  while ((line_obj = next_line())) {

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

    case '+':
    {
      const bool plus = line.Contains('p'),
                minus = line.Contains('m'),
               lambda = line.Contains('l'),
           antilambda = line.Contains('a');
      // one or the other of each
      if ((plus ^ minus) && (lambda ^ antilambda)) {
        // Code : 0b00 : minus,antilambda
        //          01 : minus,lambda
        //          10 : plus,antilambda
        //          11 : plus,lambda
        const short code = (short(plus) << 1) + short(lambda);
//         mac.analysis_configurations.insert(code);
        mac.analysis_configurations |= (1 << code);
      } else {
        std::cout << "W-BuildConfiguration: WARNING: Bad analysis_configuration '" << line << "'\n";
      }
    }
      continue;

    default:
      continue;
    }

    cmd += ";";

    cout << "I-BuildConfiguration: `" << cmd << "`\n";

    gROOT->ProcessLineFast(cmd);
  }
}
