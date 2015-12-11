///
/// \file AliFemtoAnalysisPionPion.cxx
///

#include "AliFemtoAnalysisPionPion.h"

#include "AliFemtoPionLambdaCutMonitor.h"

#include "AliESDtrack.h"

#include "AliFemtoV0TrackCut.h"
#include "AliFemtoBasicTrackCut.h"
#include "AliFemtoESDTrackCut.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoEventCutCentrality.h"
#include "AliFemtoPairCutAntiGamma.h"


#include <TROOT.h>
#include <TInterpreter.h>

#include <cassert>

static const double PionMass = 0.13956995;

const struct { unsigned int bin_count; float min; float max; }

VertexBinning = {16, -10.0, 10.0},
MultBinning = {30, 0, 10000};

const bool default_enable_pair_monitors = kTRUE,
           default_group_output_objects = kTRUE,
                 default_is_mc_analysis = kFALSE;

const UInt_t default_num_events_to_mix = 6,
                 default_min_coll_size = 100;


const float default_pion_PtMin = 0.2,
            default_pion_PtMax = 2.0,

            default_pion_EtaMin = -0.8,
            default_pion_EtaMax = 0.8,

            default_pion_DCAMin = 0.5,
            default_pion_DCAMax = 4.0,

            default_pion_NSigmaMin = -3.0,
            default_pion_NSigmaMax = 3.0;


const float default_event_EventMultMin = 0,
            default_event_EventMultMax = 100000,

            default_event_EventCentralityMin = 0,
            default_event_EventCentralityMax = 90,

            default_event_VertZPosMin = -10.0,
            default_event_VertZPosMax = 10.0,

            default_event_EPVZEROMin = -1000.0,
            default_event_EPVZEROMax = 1000.0;

const  int  default_event_TriggerSelection = 0;
const  bool default_event_AcceptBadVertex = kFALSE;


const Bool_t  default_pair_TPCOnly = kTRUE;
const Float_t default_pair_TPCExitSepMin = -1.0,
              default_pair_MinAvgSeparationPos = 0.0,
              default_pair_MinAvgSeparationNeg = 0.0;


const AliFemtoAnalysisPionPion::PionType
  default_PionType = AliFemtoAnalysisPionPion::kNone;



void
AliFemtoAnalysisPionPion::_Init(const AliFemtoAnalysisPionPion::CutParams& p = AliFemtoAnalysisPionPion::DefaultCutConfig())
{ // This function is to be called at the end of all constructors. It will
  // create default objects for any cuts that haven't been initialized

  if (fFirstParticleCut == NULL) {
    SetFirstParticleCut(BuildPionCut1(p));
  }

  if (fPionType_2 != kNone && fSecondParticleCut == NULL) {
    SetSecondParticleCut(BuildPionCut2(p));
  } else if (fPionType_2 == kNone) {
    fSecondParticleCut = fFirstParticleCut;
  }

  if (fEventCut == NULL) {
    SetEventCut(BuildEventCut(p));
  }

  if (fPairCut == NULL) {
    SetPairCut(BuildPairCut(p));
  }
}

AliFemtoAnalysisPionPion::AliFemtoAnalysisPionPion():
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count,
                             VertexBinning.min,
                             VertexBinning.max,
                             MultBinning.bin_count,
                             MultBinning.min,
                             MultBinning.max)
  , fAnalysisName("AliFemtoAnalysisPionPion")
  , fPionType_1(default_PionType)
  , fPionType_2(default_PionType)
  , fGroupOutputObjects(default_group_output_objects)
  , fMCAnalysis(default_is_mc_analysis)
{
  SetEnablePairMonitors(default_enable_pair_monitors);
  SetNumEventsToMix(default_num_events_to_mix);
  _Init();
}

AliFemtoAnalysisPionPion::AliFemtoAnalysisPionPion(const char *name):
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count,
                             VertexBinning.min,
                             VertexBinning.max,
                             MultBinning.bin_count,
                             MultBinning.min,
                             MultBinning.max)
  , fAnalysisName(name)
  , fPionType_1(default_PionType)
  , fPionType_2(default_PionType)
  , fGroupOutputObjects(default_group_output_objects)
  , fMCAnalysis(default_is_mc_analysis)
{
  SetEnablePairMonitors(default_enable_pair_monitors);
  _Init();
}


AliFemtoAnalysisPionPion
  ::AliFemtoAnalysisPionPion(const char *name,
                             const PionType pion_1,
                             const PionType pion_2):
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count,
                             VertexBinning.min,
                             VertexBinning.max,
                             MultBinning.bin_count,
                             MultBinning.min,
                             MultBinning.max)
  , fAnalysisName(name)
  , fPionType_1(pion_1)
  , fPionType_2(pion_2)
  , fGroupOutputObjects(default_group_output_objects)
  , fMCAnalysis(default_is_mc_analysis)
{
  SetEnablePairMonitors(default_enable_pair_monitors);
  SetNumEventsToMix(default_num_events_to_mix);
  _Init();
}

AliFemtoAnalysisPionPion
  ::AliFemtoAnalysisPionPion(const char *name,
                             const PionType pion_1,
                             const PionType pion_2,
                             const CutParams& params):
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count,
                             VertexBinning.min,
                             VertexBinning.max,
                             MultBinning.bin_count,
                             MultBinning.min,
                             MultBinning.max)
  , fAnalysisName(name)
  , fPionType_1(pion_1)
  , fPionType_2(pion_2)
  , fGroupOutputObjects(default_group_output_objects)
  , fMCAnalysis(default_is_mc_analysis)
{
  SetEnablePairMonitors(default_enable_pair_monitors);
  SetNumEventsToMix(default_num_events_to_mix);
  _Init(params);
}

AliFemtoAnalysisPionPion::AliFemtoAnalysisPionPion(const char *name,
                                                   const AnalysisParams &params,
                                                   const CutParams &cut_params):
  AliFemtoVertexMultAnalysis(params.vertex_bins,
                             params.vertex_min,
                             params.vertex_max,
                             params.mult_bins,
                             params.mult_min,
                             params.mult_max)
  , fAnalysisName(name)
  , fPionType_1(params.pion_type_1)
  , fPionType_2(params.pion_type_2)
  , fGroupOutputObjects(params.group_output_objects)
  , fMCAnalysis(params.is_mc_analysis)
{
  SetEnablePairMonitors(params.enable_pair_monitors);
  SetNumEventsToMix(params.num_events_to_mix);
//   cout << "[AliFemtoAnalysisPionPion] Constructed with : " <<
//     params.vertex_bins <<
//     " " << params.vertex_min <<
//     " " << params.vertex_max <<
//     " " << params.mult_bins <<
//     " " << params.mult_min <<
//     " " << params.mult_max <<'\n';
  _Init(cut_params);
}


AliFemtoAnalysisPionPion::AnalysisParams
AliFemtoAnalysisPionPion::DefaultConfig()
{
  AliFemtoAnalysisPionPion::AnalysisParams params = {
    VertexBinning.bin_count
  , VertexBinning.min
  , VertexBinning.max

  , MultBinning.bin_count
  , MultBinning.min
  , MultBinning.max

  , default_PionType
  , default_PionType

  , default_num_events_to_mix
  , default_min_coll_size

  , default_enable_pair_monitors
  , default_group_output_objects
  , default_is_mc_analysis
  };

  assert(params.min_coll_size == default_min_coll_size);
  assert(params.is_mc_analysis == default_is_mc_analysis);

  return params;
}


AliFemtoAnalysisPionPion::CutParams
AliFemtoAnalysisPionPion::DefaultCutConfig()
{
  AliFemtoAnalysisPionPion::CutParams params = {
    // Event
    default_event_EventMultMin
  , default_event_EventMultMax
  , default_event_EventCentralityMin
  , default_event_EventCentralityMax
  , default_event_VertZPosMin
  , default_event_VertZPosMax
  , default_event_EPVZEROMin
  , default_event_EPVZEROMax
  , default_event_TriggerSelection
  , default_event_AcceptBadVertex

    // Pion 1
  , default_pion_PtMin
  , default_pion_PtMax
  , default_pion_EtaMin
  , default_pion_EtaMax
  , default_pion_DCAMin
  , default_pion_DCAMax

  , default_pion_NSigmaMin
  , default_pion_NSigmaMax

    // Pion 2
  , default_pion_PtMin
  , default_pion_PtMax
  , default_pion_EtaMin
  , default_pion_EtaMax
  , default_pion_DCAMin
  , default_pion_DCAMax

  , default_pion_NSigmaMin
  , default_pion_NSigmaMax

    // pair
  , default_pair_TPCOnly
  , default_pair_TPCExitSepMin
  , default_pair_MinAvgSeparationPos
  , default_pair_MinAvgSeparationNeg
  };

  // sanity checks
  assert(params.event_MultMin == default_event_EventMultMin);
  assert(params.pion_1_PtMin == default_pion_PtMin);
  assert(params.pion_2_PtMin == default_pion_PtMin);
  assert(params.pair_TPCOnly == default_pair_TPCOnly);
  assert(params.pair_TPCExitSepMin == default_pair_TPCExitSepMin);

  return params;
}

AliFemtoTrackCut*
AliFemtoAnalysisPionPion::BuildPionCut1(const CutParams &p) const
{
  const int charge = (fPionType_1 == kPiMinus) ? -1 :
                      (fPionType_1 == kPiPlus) ? +1 :
                                                -999;
  if (charge == -999) {
    std::cerr << "E-AliFemtoAnalysisPionPion::BuildPionCut: Invalid pion type: '" << fPionType_1 << "'\n";
  }

//   AliFemtoBasicTrackCut *cut = new AliFemtoBasicTrackCut();
//   cut->SetCharge(charge);
//   cut->SetMass(PionMass);
//   cut->SetNSigmaPion(p.pion_1_NSigmaMin, p.pion_1_NSigmaMax);
//   cut->SetPt(p.pion_1_PtMin, p.pion_1_PtMax);
//   cut->SetRapidity(p.pion_1_EtaMin, p.pion_1_EtaMax);
//   cut->SetDCA(p.pion_1_DCAMin, p.pion_1_DCAMax);
  
  AliFemtoESDTrackCut *cut = new AliFemtoESDTrackCut();
  cut->SetCharge(charge);
  cut->SetMass(PionMass);
  cut->SetEta(p.pion_1_EtaMin, p.pion_1_EtaMax);
  cut->SetMostProbablePion();
//   cut->SetStatus(AliESDtrack::kTPCrefit | AliESDtrack::kITSrefit);


  /// Settings for TPC-Inner Runmode
  cut->SetStatus(AliESDtrack::kTPCin);
  cut->SetminTPCncls(80);
  cut->SetRemoveKinks(kTRUE);
  cut->SetLabel(kFALSE);
  cut->SetMaxTPCChiNdof(0.032);
  cut->SetMaxImpactXY(2.4);
  cut->SetMaxImpactZ(3.0);
  return cut;
}


AliFemtoTrackCut*
AliFemtoAnalysisPionPion::BuildPionCut2(const CutParams &p) const
{
  const int charge = (fPionType_2 == kPiMinus) ? -1 :
                      (fPionType_2 == kPiPlus) ? +1 :
                                                -999;
  if (charge == -999) {
    std::cerr << "E-AliFemtoAnalysisPionPion::BuildPionCut: Invalid pion type: '" << fPionType_2 << "'\n";
  }

//   AliFemtoBasicTrackCut *cut = new AliFemtoBasicTrackCut();
//   cut->SetCharge(charge);
//   cut->SetMass(PionMass);
//   cut->SetNSigmaPion(p.pion_2_NSigmaMin, p.pion_2_NSigmaMax);
//   cut->SetPt(p.pion_2_PtMin, p.pion_2_PtMax);
//   cut->SetRapidity(p.pion_2_EtaMin, p.pion_2_EtaMax);
//   cut->SetDCA(p.pion_2_DCAMin, p.pion_2_DCAMax);

  AliFemtoESDTrackCut *cut = new AliFemtoESDTrackCut();
  cut->SetCharge(charge);
  cut->SetMass(PionMass);
  cut->SetEta(p.pion_2_EtaMin, p.pion_2_EtaMax);
  cut->SetMostProbablePion();

  /// Settings for TPC-Inner Runmode
  cut->SetStatus(AliESDtrack::kTPCin);
  cut->SetminTPCncls(80);
  cut->SetRemoveKinks(kTRUE);
  cut->SetLabel(kFALSE);
  cut->SetMaxTPCChiNdof(2.0);
  cut->SetMaxImpactXY(2.4);
  cut->SetMaxImpactZ(3.0);
  return cut;
}


AliFemtoEventCut*
AliFemtoAnalysisPionPion::BuildEventCut(const AliFemtoAnalysisPionPion::CutParams& p) const
{
/*
  AliFemtoBasicEventCut *cut = new AliFemtoBasicEventCut();

  cut->SetEventMult(p.event_MultMin,
                    p.event_MultMax);
  cut->SetVertZPos(p.event_VertexZMin,
                   p.event_VertexZMax);
  cut->SetEPVZERO(p.event_EP_VZeroMin,
                  p.event_EP_VZeroMax);
  cut->SetTriggerSelection(p.event_TriggerSelection);
  cut->SetAcceptBadVertex(p.event_AcceptBadVertex);
*/
  AliFemtoEventCutCentrality *cut = new AliFemtoEventCutCentrality();
  cut->SetCentralityRange(p.event_CentralityMin, p.event_CentralityMax);
  cut->SetZPosRange(p.event_VertexZMin, p.event_VertexZMax);
  cut->SetEPVZERO(p.event_EP_VZeroMin, p.event_EP_VZeroMax);
  cut->SetTriggerSelection(p.event_TriggerSelection);

//   if (!fMCAnalysis) {
//     cut->SetAcceptOnlyPhysics(kTRUE);
//   }

  return cut;
}


AliFemtoPairCutAntiGamma*
AliFemtoAnalysisPionPion::BuildPairCut(const CutParams &p) const
{
  AliFemtoPairCutAntiGamma *cut = new AliFemtoPairCutAntiGamma();

  cut->SetShareQualityMax(1.0);
  cut->SetShareFractionMax(0.05);
  cut->SetRemoveSameLabel(kFALSE);

  // cut->SetTPCOnly(p.pair_TPCOnly);
  // cut->SetTPCExitSepMinimum(p.pair_TPCExitSepMin);
  // cut->SetMinAvgSeparation(0, p.pair_MinAvgSeparationPos);
  // cut->SetMinAvgSeparation(1, p.pair_MinAvgSeparationNeg);

  return cut;
}

void AliFemtoAnalysisPionPion::AddStanardCutMonitors()
{
  if (fEventCut) {
    fEventCut->AddCutMonitor(new AliFemtoPionLambdaCutMonitor::Event(true, fMCAnalysis, false),
                             new AliFemtoPionLambdaCutMonitor::Event(false, fMCAnalysis, false));
  } else {
    std::cout << " NO fEventCut!\n";
    exit(1);
  }

  const TString p1_type_str = (fPionType_1 == kPiPlus) ? "Pi+" : (fPionType_1 == kPiMinus) ? "Pi-" : "ERROR";

  if (fFirstParticleCut) {
    fFirstParticleCut->AddCutMonitor(new AliFemtoPionLambdaCutMonitor::Pion(true, p1_type_str, fMCAnalysis),
                                     new AliFemtoPionLambdaCutMonitor::Pion(true, p1_type_str, fMCAnalysis));
  }

  if (!AnalyzeIdenticalParticles()) {
    const TString p2_type_str = (fPionType_2 == kPiPlus) ? "Pi+" : (fPionType_2 == kPiMinus) ? "Pi-" : "ERROR";
    fSecondParticleCut->AddCutMonitor(new AliFemtoPionLambdaCutMonitor::Pion(true, p2_type_str, fMCAnalysis),
                                      new AliFemtoPionLambdaCutMonitor::Pion(false, p2_type_str, fMCAnalysis));

    // Non-Identical Pion Pairs
    if (fPairCut) {
      const TString pair_type_str = p1_type_str + "/" + p2_type_str;
      fPairCut->AddCutMonitor(new AliFemtoPionLambdaCutMonitor::Pair(true, pair_type_str, fMCAnalysis),
                              new AliFemtoPionLambdaCutMonitor::Pair(false, pair_type_str, fMCAnalysis));
    }
  }
  // Identical Pion Pairs
  else if (fPairCut) {
    const TString pair_type_str = "(Identical) " + p1_type_str + "/" + p1_type_str;
    fPairCut->AddCutMonitor(new AliFemtoPionLambdaCutMonitor::Pair(true, pair_type_str, fMCAnalysis),
                            new AliFemtoPionLambdaCutMonitor::Pair(false, pair_type_str, fMCAnalysis));
  }

}

static TObjArray* GetPassFailOutputList(const TString &name,
                                        AliFemtoCutMonitorHandler *handler)
{
  AliFemtoCutMonitorCollection *p_coll = handler->PassMonitorColl(),
                               *f_coll = handler->FailMonitorColl();

  TObjArray *res = new TObjArray(),
            *passout = new TObjArray(),
            *failout = new TObjArray();

  res->SetName(name);
  passout->SetName("pass");
  failout->SetName("fail");

  res->Add(passout);
  res->Add(failout);

  for (AliFemtoCutMonitorIterator it = p_coll->begin();
                                  it != p_coll->end();
                                ++it) {
    TList *o_list = (*it)->GetOutputList();
    passout->AddAll(o_list);
    delete o_list;
  }

  for (AliFemtoCutMonitorIterator it = f_coll->begin();
                                  it != f_coll->end();
                                ++it) {
    TList *o_list = (*it)->GetOutputList();
    failout->AddAll(o_list);
    delete o_list;
  }

  return res;
}

TList* AliFemtoAnalysisPionPion::GetOutputList()
{
  TList *outputlist = NULL;
  TSeqCollection *output = NULL;

  // get "standard" outputs - that's what we output
  if (!fGroupOutputObjects) {
    outputlist = AliFemtoVertexMultAnalysis::GetOutputList();
    output = outputlist;
  } else {

    output = new TObjArray();
    output->SetName(fAnalysisName);

    output->Add(GetPassFailOutputList("Event", EventCut()));

    if (AnalyzeIdenticalParticles()) {
      output->Add(GetPassFailOutputList("Tracks", fFirstParticleCut));
    } else {
      output->Add(GetPassFailOutputList("Track1", fFirstParticleCut));
      output->Add(GetPassFailOutputList("Track2", fSecondParticleCut));
    }

    output->Add(GetPassFailOutputList("Pair", PairCut()));

    for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
                                  iter != fCorrFctnCollection->end();
                                  ++iter) {
      TList *cf_output = (*iter)->GetOutputList();
      output->AddAll(cf_output);
      delete cf_output;
    }

    outputlist = new TList();
    outputlist->Add(output);
  }

  // get the list of settings from each cut and correlation function
  {
    TString settings("== settings ==\n");

    TList *setting_list = ListSettings();
    TIter next_setting(setting_list);
    for (TObject *setting = next_setting();
                  setting != NULL;
                  setting = next_setting()) {
      TObjString *setting_str = (TObjString*)setting;
      settings += setting_str->String() + "\n";
    }
    delete setting_list;
    output->Add(new TObjString(settings));
  }
  return outputlist;
}

TList* AliFemtoAnalysisPionPion::ListSettings()
{
  TList *setting_list = new TList();


  setting_list->AddVector(

    new TObjString(
      TString::Format("AliFemtoAnalysisPionPion.piontype=%d", fPionType_1)
    ),

    new TObjString(
      TString::Format("AliFemtoAnalysisPionPion.lambdatype=%d", fPionType_2)
    ),

    new TObjString(
      TString::Format("AliFemtoAnalysisPionPion.mc_analysis=%d", fMCAnalysis)
    ),

  NULL);

  TList *parent_list = AliFemtoVertexMultAnalysis::ListSettings();

  setting_list->AddAll(parent_list);
  delete parent_list;

  setting_list->SetOwner(kTRUE);
  return setting_list;
}

void AliFemtoAnalysisPionPion::EventBegin(const AliFemtoEvent* ev)
{
  fEventCut->EventBegin(ev);
  fFirstParticleCut->EventBegin(ev);

  if (!AnalyzeIdenticalParticles()) {
    fSecondParticleCut->EventBegin(ev);
  }

  fPairCut->EventBegin(ev);
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
                                iter != fCorrFctnCollection->end();
                                ++iter) {
    (*iter)->EventBegin(ev);
  }
}

void AliFemtoAnalysisPionPion::EventEnd(const AliFemtoEvent* ev)
{
  fEventCut->EventEnd(ev);
  fFirstParticleCut->EventEnd(ev);
  if (!AnalyzeIdenticalParticles()) {
    fSecondParticleCut->EventEnd(ev);
  }
  fPairCut->EventEnd(ev);
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
                                iter != fCorrFctnCollection->end();
                                ++iter) {
    (*iter)->EventEnd(ev);
  }
}
