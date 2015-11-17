///
/// \file PionPionFemto/ConfigFemtoAnalysis.C
///
/// \brief The configuration macro which sets up identical pion-pion analyses
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

#include "AliFemtoCorrFctnKStar.h"

#include <TROOT.h>

#endif

#include <vector>

// typedef std::pair<float, float> pair_of_floats;


struct AnalysisParams;
struct CutParams;
struct MacroParams;


void BuildConfiguration(const TString&, AnalysisParams&, CutParams&, MacroParams&);


AliFemtoManager*
ConfigFemtoAnalysis(const TString& param_str = "")
{
  const double PionMass = 0.13956995;

  // Get the default configurations
  AnalysisParams analysis_config;
  CutParams cut_config;
  MacroParams macro_config;

  // Read parameter string and update configurations
  BuildConfiguration(param_str, analysis_config, cut_config, macro_config);

  if (analysis_config.mult_ranges.empty()) {
    // analysis_config.mult_ranges.push_back(std::pair<float, float> (0.0, 300));
    analysis_config.mult_ranges.push_back(0.0);
    analysis_config.mult_ranges.push_back(300.0);
  }

  // Begin to build the manager and analyses
  AliFemtoManager *manager = new AliFemtoManager();

  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAODChain();
    rdr->SetFilterBit(7);
    rdr->SetEPVZERO(kTRUE);
    rdr->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
    rdr->SetCentralityFlattening(kFALSE);
    rdr->SetReadV0(0);
    rdr->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
  // rdr->SetReadMC(analysis_config.is_mc_analysis);
  manager->SetEventReader(rdr);

  int runmults[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  int multbins[11] = {0.001, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900};

  const char* chrgs[] = { "_pip_", "_pim_", "_pimpip_" };
  int runtype = 1;
  int runshlcms = 1;
  int runqinv = 1;
  int runktdep = 0;
  int aniter = -1;

  double shqmax = 1.2;
  int nbinssh = 100;
  int imult = 0;

  double ktrng[3] = {0.01, 0.7, 100.0};
  int numOfkTbins = 2;

 //  for (std::vector<std::pair<float, float>>::iterator mult_range = analysis_config.mult_ranges.begin();
//        mult_range != analysis_config.mult_ranges.end();
//        ++mult_range) {

  for (int mult_it = 0; mult_it + 1 < analysis_config.mult_ranges.size(); mult_it += 2) {

    int ichg = 0;
    aniter++;
//     const double mult_low = mult_range->first,
//                 mult_high = mult_range->second;

    const double mult_low = analysis_config.mult_ranges[mult_it],
                mult_high = analysis_config.mult_ranges[mult_it + 1];

    AliFemtoSimpleAnalysis *analysis = new AliFemtoVertexMultAnalysis(
      analysis_config.vertex_bins,
      analysis_config.vertex_min,
      analysis_config.vertex_max,

      analysis_config.mult_bins,
      mult_low,
      mult_high
    );
    analysis->SetVerboseMode(kFALSE);

    analysis->SetNumEventsToMix(10 /*analysis_config.num_events_to_mix*/);
    analysis->SetMinSizePartCollection(1 /*analysis_config.min_coll_size*/);

    AliFemtoBasicEventCut *ev_cut = new AliFemtoBasicEventCut();
    ev_cut->SetEventMult(0.001, 100000);
    ev_cut->SetVertZPos(-8, 8);

    // if (!analysis_config.is_mc_analysis) {
    //   ev_cut->SetAcceptOnlyPhysics(kTRUE);
    // }

    ev_cut->AddCutMonitor(
      new AliFemtoCutMonitorEventMult(TString::Format("cutPass%stpcM%i", chrgs[ichg], imult)),
      new AliFemtoCutMonitorEventMult(TString::Format("cutFail%stpcM%i", chrgs[ichg], imult))
    );
    ev_cut->AddCutMonitor(
      new AliFemtoCutMonitorEventVertex(TString::Format("cutPass%stpcM%i", chrgs[ichg], imult)),
      new AliFemtoCutMonitorEventVertex(TString::Format("cutFail%stpcM%i", chrgs[ichg], imult))
    );
    analysis->SetEventCut(ev_cut);

    //
    // Setup Track Cuts
    //

    AliFemtoBasicTrackCut *track_cut1 = new AliFemtoBasicTrackCut(),
                          *track_cut2 = (ichg == 2) ? new AliFemtoBasicTrackCut() : NULL;
    if (ichg == 2) {
      analysis->SetFirstParticleCut(track_cut1);
      analysis->SetSecondParticleCut(track_cut2);
    } else {
      analysis->SetFirstParticleCut(track_cut1);
      analysis->SetSecondParticleCut(track_cut1);
    }

    switch (ichg) {
    case 0:
      track_cut1->SetCharge(1.0);
      break;
    case 1:
      track_cut1->SetCharge(-1.0);
      break;
    case 2:
      track_cut1->SetCharge(-1.0);
      track_cut2->SetCharge(1.0);
    }

//     track_cut1->SetEta(-0.8, 0.8);
    track_cut1->SetMass(PionMass);
    track_cut1->SetNSigmaPion(-0.02, 0.02);
    track_cut1->SetPt(0.3, 1.0);
    track_cut1->SetRapidity(-0.8, 0.8);
    track_cut1->SetDCA(0.5, 4.0);
//     track_cut1->SetMostProbablePion();

  track_cut1->AddCutMonitor(
    new AliFemtoCutMonitorParticleYPt("PionPass", PionMass),
    new AliFemtoCutMonitorParticleYPt("PionFail", PionMass)
  );
    

    // Track quality cuts
/*
    switch (runtype) {
    case 0:
      track_cut1->SetStatus(AliESDtrack::kTPCrefit | AliESDtrack::kITSrefit);
      track_cut1->SetminTPCncls(80);
      track_cut1->SetRemoveKinks(kTRUE);
      track_cut1->SetLabel(kFALSE);
      track_cut1->SetMaxTPCChiNdof(4.0);
      track_cut1->SetMaxImpactXY(0.2);
      track_cut1->SetMaxImpactZ(0.15);
      break;
    case 1:
      track_cut1->SetStatus(AliESDtrack::kITSrefit);
      track_cut1->SetRemoveKinks(kTRUE);
      track_cut1->SetLabel(kFALSE);
      track_cut1->SetMaxImpactXY(0.2);
      track_cut1->SetMaxImpactZ(0.25);
      break;
    case 2:
      track_cut1->SetStatus(AliESDtrack::kTPCin);
      track_cut1->SetminTPCncls(80);
      track_cut1->SetRemoveKinks(kTRUE);
      track_cut1->SetLabel(kFALSE);
      track_cut1->SetMaxTPCChiNdof(4.0);
      track_cut1->SetMaxImpactXY(2.4);
      track_cut1->SetMaxImpactZ(3.0);

      if (ichg == 2) {
        track_cut2->SetStatus(AliESDtrack::kTPCin);
        track_cut2->SetminTPCncls(80);
        track_cut2->SetRemoveKinks(kTRUE);
        track_cut2->SetLabel(kFALSE);
        track_cut2->SetMaxTPCChiNdof(4.0);
        track_cut2->SetMaxImpactXY(0.25);
        track_cut2->SetMaxImpactZ(0.15);
      }
      break;
    }
*/

    AliFemtoPairCutAntiGamma *pair_cut = new AliFemtoPairCutAntiGamma();

    switch (runtype) {
    case 0:
      pair_cut->SetShareQualityMax(1.0);
      pair_cut->SetShareFractionMax(0.05);
      pair_cut->SetRemoveSameLabel(kFALSE);
      break;
    case 1:
      pair_cut->SetShareQualityMax(1.0);
      pair_cut->SetShareFractionMax(1.05);
      pair_cut->SetRemoveSameLabel(kFALSE);
      break;
    case 2:
      // pair_cut->SetUseAOD(kTRUE);
      pair_cut->SetShareQualityMax(1.0);
      pair_cut->SetShareFractionMax(0.05);
      pair_cut->SetRemoveSameLabel(kFALSE);

      if (gammacut == 0) {
        pair_cut->SetMaxEEMinv(0.0);
        pair_cut->SetMaxThetaDiff(0.0);
      } else if (gammacut == 1) {
        pair_cut->SetMaxEEMinv(0.002);
        pair_cut->SetMaxThetaDiff(0.008);
      }
      break;
    }

    analysis->SetPairCut(pair_cut);
/*
    AliFemtoCorrFctn *cf = (ichg == 2)
                         ? new AliFemtoCorrFctnNonIdDR(
                             Form("ckstar%stpcM%i", chrgs[ichg], imult),
                             nbinssh, 0.0, shqmax);
                         : new AliFemtoQinvCorrFctn(
                             Form("cqinv%stpcM%i", chrgs[ichg], imult),
                             nbinssh, 0.0, shqmax);
*/
    AliFemtoCorrFctn *cf = new AliFemtoQinvCorrFctn(
                             Form("cqinv%stpcM%i", chrgs[ichg], imult),
                             nbinssh, 0.0, shqmax);
    analysis->AddCorrFctn(cf);

    analysis->AddCorrFctn(
      new AliFemtoCorrFctnDirectYlm(
        Form("cylm%stpcM%i", chrgs[ichg], imult),
        3,
        nbinssh,
        0.0,
        shqmax,
        runshlcms
      )
    );

    if (runktdep) {
      for (int ikt = 0; ikt < numOfkTbins; ikt++) {
        int ktm = aniter * numOfkTbins + ikt;

        AliFemtoKTPairCut *kt_cut = new AliFemtoKTPairCut(ktrng[ikt], ktrng[ikt+1]);
        AliFemtoCorrFctnDirectYlm *ylm_cf = new AliFemtoCorrFctnDirectYlm(
          Form("cylm%stpcM%ikT%i", chrgs[ichg], imult, ikt),
          3,
          nbinssh,
          0.0,
          shqmax,
          0
        );

        ylm_cf->SetPairSelectionCut(kt_cut);
        analysis->AddCorrFctn(ylm_cf);


        AliFemtoCorrFctn *qcf = new AliFemtoQinvCorrFctn(
          Form("cqinv%stpcM%ikT%i", chrgs[ichg], imult, ikt),
          nbinssh,
          0.0,
          (imult > 6) ? shqmax * 2.5 : shqmax
        );

        qcf->SetPairSelectionCut(kt_cut);
        analysis->AddCorrFctn(qcf);
      }
    }

    manager->AddAnalysis(analysis);

  }

  return manager;
}

struct AnalysisParams {

  Bool_t is_mc_analysis;

  UInt_t vertex_bins;
  Float_t vertex_min,
          vertex_max;

  UInt_t mult_bins;
//   std::vector< std::pair<float, float> > mult_ranges;
  std::vector<float> mult_ranges;

  UInt_t num_events_to_mix;
  UInt_t min_coll_size;

  AnalysisParams():
    is_mc_analysis(kFALSE)
  , vertex_bins(8)
  , vertex_min(-8)
  , vertex_max(8)
  , mult_bins(4)
  , mult_ranges()
  {}

};

struct CutParams {


};

struct MacroParams {

};


void
BuildConfiguration(const TString &text,
                   AnalysisParams &a,
                   CutParams &cut,
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
//         const short code = (short(plus) << 1) + short(lambda);
//         mac.analysis_configurations.insert(code);
//         mac.analysis_configurations |= (1 << code);
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
