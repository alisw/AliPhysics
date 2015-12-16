///
/// \file PionPionFemto/ConfigFemtoAnalysis.C
///
/// \brief The configuration macro which sets up identical pion-pion analyses
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

#if !defined(__CINT__) || defined(__MAKECINT_)

#include "AliFemtoAnalysisPionPion.h"

#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"

#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorParticleYPt.h"

#include "AliFemtoCorrFctnKStar.h"

#include <TROOT.h>

#endif

typedef AliFemtoAnalysisPionPion AFAPP;

bool DEFAULT_DO_KT = kFALSE;

struct MacroParams {
  std::vector<int> centrality_ranges;
  std::vector<unsigned char> pair_codes;
  bool do_kt_q3d;
  bool do_kt_qinv;
};

void
BuildConfiguration(
  const TString&,
  AliFemtoAnalysisPionPion::AnalysisParams&,
  AliFemtoAnalysisPionPion::CutParams&,
  MacroParams&
);


AliFemtoManager*
ConfigFemtoAnalysis(const TString& param_str = "")
{
  std::cout << "[ConfigFemtoAnalysis (PionPion)]\n";

  const double PionMass = 0.13956995;

  // Get the default configurations
  AFAPP::AnalysisParams analysis_config = AFAPP::DefaultConfig();
  AFAPP::CutParams cut_config = AFAPP::DefaultCutConfig();
  MacroParams macro_config;

  macro_config.do_kt_q3d = macro_config.do_kt_qinv = DEFAULT_DO_KT;

  // Read parameter string and update configurations
  BuildConfiguration(param_str, analysis_config, cut_config, macro_config);

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

  if (macro_config.centrality_ranges.empty()) {
    // macro_config.centrality_ranges.push_back(std::pair<int, int>(0.0, 300));
    macro_config.centrality_ranges.push_back(0);
    macro_config.centrality_ranges.push_back(10.0);
  }

  // default to identical pi+ and pi- analyses
  if (macro_config.pair_codes.empty()) {
    macro_config.pair_codes.push_back(1);
    macro_config.pair_codes.push_back(0);
  }


  AFAPP::PionType PI_PLUS = AliFemtoAnalysisPionPion::kPiPlus,
                 PI_MINUS = AliFemtoAnalysisPionPion::kPiMinus,
                     NOPI = AliFemtoAnalysisPionPion::kNone;

  // loop over centrality ranges
  for (int cent_it = 0; cent_it + 1 < macro_config.centrality_ranges.size(); cent_it += 2) {

    const int mult_low = macro_config.centrality_ranges[cent_it],
             mult_high = macro_config.centrality_ranges[cent_it + 1];

    const TString mult_low_str = TString::Format("%0.2i", mult_low),
                 mult_high_str = TString::Format("%0.2i", mult_high);

    // loop over pair types
    for (int pair_it = 0; pair_it < macro_config.pair_codes.size(); ++pair_it) {

      const int pair_code = macro_config.pair_codes[pair_it];

      AFAPP::PionType type_1 = (pair_code & 0x1) ? PI_PLUS : PI_MINUS,
                      type_2 = NOPI;

      switch (pair_code & 0x6) {
      case 0:
        break;
      case 2:
        type_2 = PI_PLUS;
        break;
      case 4:
        type_2 = PI_MINUS;
        break;
      default:
        cout << "W-ConfigFemtoAnalysis: Invalid pair code " << pair_code << ". Skipping.\n";
        continue;
      }

      const TString pair_type_str =
        TString(type_1 == AliFemtoAnalysisPionPion::kPiPlus ? "pip" : "pim")
        + (type_2 == AliFemtoAnalysisPionPion::kNone ? "" : type_2 == AliFemtoAnalysisPionPion::kPiPlus ? "_pip" : "_pim");

      // build unique analysis name from centrality and pair types
      const TString analysis_name = TString::Format(
        "PiPiAnalysis_%0.2i_%0.2i_%s", mult_low, mult_high, pair_type_str.Data()
      );

      analysis_config.pion_type_1 = type_1;
      analysis_config.pion_type_2 = type_2;

      cut_config.event_CentralityMin = mult_low;
      cut_config.event_CentralityMax = mult_high;

      AliFemtoAnalysisPionPion *analysis = new AliFemtoAnalysisPionPion(analysis_name, analysis_config, cut_config);
      analysis->SetVerboseMode(kTRUE);

      analysis->AddStanardCutMonitors();

      TString cf_title = TString("_qinv_") + pair_type_str;
      AliFemtoCorrFctn *cf = new AliFemtoQinvCorrFctn(cf_title.Data(), 100, 0.0, 1.2);
      analysis->AddCorrFctn(cf);

      cf_title = TString("_q3D_") + pair_type_str;
      analysis->AddCorrFctn(new AliFemtoCorrFctn3DLCMSSym(cf_title, 72, 0.25));

      if (macro_config.do_kt_qinv) {
        AliFemtoKtBinnedCorrFunc *binned = new AliFemtoKtBinnedCorrFunc("KT_Qinv", new AliFemtoQinvCorrFctn(*(AliFemtoQinvCorrFctn*)cf));
        binned->AddKtRange(0.2, 0.4);
        binned->AddKtRange(0.4, 0.6);
        binned->AddKtRange(0.6, 1.0);
        analysis->AddCorrFctn(binned);
      }

      if (macro_config.do_kt_q3d) {
        AliFemtoKtBinnedCorrFunc *kt_q3d = new AliFemtoKtBinnedCorrFunc("KT_Q3D", new AliFemtoCorrFctn3DLCMSSym(TString("q3D_") + pair_type_str, 72, 1.1));
        kt_q3d->AddKtRange(0.2, 0.4);
        kt_q3d->AddKtRange(0.4, 0.6);
        kt_q3d->AddKtRange(0.6, 1.0);
        analysis->AddCorrFctn(kt_q3d);
      }

      // analysis->AddCorrFctn(
      //   new AliFemtoCorrFctnDirectYlm(
      //     Form("cylm%stpcM%i", chrgs[ichg], imult),
      //     3,
      //     nbinssh,
      //     0.0,
      //     shqmax,
      //     runshlcms
      //   )
      // );
      //

      manager->AddAnalysis(analysis);


      continue;

      // analysis->SetNumEventsToMix(10 /*analysis_config.num_events_to_mix*/);
      // analysis->SetMinSizePartCollection(1 /*analysis_config.min_coll_size*/);

      // AliFemtoBasicEventCut *ev_cut = new AliFemtoBasicEventCut();
      // ev_cut->SetEventMult(0.001, 100000);
      // ev_cut->SetVertZPos(-8, 8);

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

    // KtRange()
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

    }
  }

  return manager;
}


void
BuildConfiguration(const TString &text,
                   AliFemtoAnalysisPionPion::AnalysisParams &a,
                   AliFemtoAnalysisPionPion::CutParams &cut,
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

    const TString line = ((TObjString*)line_obj)->String()
                            .ReplaceAll(".", "_")
                            .Strip(TString::kBoth, ' ');

    TString cmd("");

    switch (line[0]) {
    case '$':
      cmd = cut_varname + "." + line(1, line.Length() - 1);
      break;

    case '@':
      cmd = analysis_varname + "." + line(1, line.Length() - 1);
      break;

    case '~':
      cmd = macro_varname + "." + line(1, line.Length() - 1);
      break;

    case '{':
    {
      UInt_t rangeend = text.Index("}");
      if (rangeend == -1) {
        rangeend = line.Length();
      }
      TString centrality_ranges = line(1, rangeend - 1);
      TObjArray *range_groups = centrality_ranges.Tokenize(",");

      TIter next_range_group(range_groups);
      TObjString *range_group = NULL;

      while (range_group = (TObjString*)next_range_group()) {
        TObjArray *subrange = range_group->String().Tokenize(":");
        TIter next_subrange(subrange);
        TObjString *subrange_it = (TObjString *)next_subrange();
        TString prev = TString::Format("%0.2d", subrange_it->String().Atoi());
        while (subrange_it = (TObjString *)next_subrange()) {
          TString next = TString::Format("%0.2d", subrange_it->String().Atoi());

          cmd = macro_varname + ".centrality_ranges.push_back(" + prev + ");";
          gROOT->ProcessLineFast(cmd);

          cmd = macro_varname + ".centrality_ranges.push_back(" + next + ");";
          gROOT->ProcessLineFast(cmd);
          prev = next;
        }
      }
    }
    continue;

    case '+':
      if (line == "+p") {
        cmd = macro_varname + ".pair_codes.push_back(1)";
      }
      else if (line == "+m") {
        cmd = macro_varname + ".pair_codes.push_back(0)";
      }
      else if (line == "+pp") {
        cmd = macro_varname + ".pair_codes.push_back(3)";
      }
      else if (line == "+mm") {
        cmd = macro_varname + ".pair_codes.push_back(4)";
      }
      else if (line == "+pm") {
        cmd = macro_varname + ".pair_codes.push_back(2)";
      }
      else if (line == "+mp") {
        cmd = macro_varname + ".pair_codes.push_back(5)";
      }
      else {
        continue;
      }
      break;

    default:
      continue;
    }

    cmd += ";";

    cout << "I-BuildConfiguration: `" << cmd << "`\n";
    gROOT->ProcessLineFast(cmd);
  }
}
