///
/// \file AliFemtoAnalysisPionLambda.cxx
///

#include "AliFemtoAnalysisPionLambda.h"


#include "AliESDtrack.h"

#include "AliFemtoV0TrackCut.h"
#include "AliFemtoBasicTrackCut.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoV0TrackPairCut.h"

#include "AliFemtoPionLambdaCutMonitor.h"

#include <TROOT.h>
#include <TInterpreter.h>

#include <cassert>

static const double LambdaMass = 1.115683,
                      PionMass = 0.13956995;

const struct { unsigned int bin_count; float min; float max; }

VertexBinning = {16, -10.0, 10.0},
MultBinning = {30, 0, 10000};

const bool default_group_output_objects = kTRUE;

const float default_lambda_PtMin = 0.2,
            default_lambda_PtMax = 5.0,
            default_lambda_Eta = 0.8,

            default_lambda_MaxDCAV0 = 5.0,
            default_lambda_MaxV0DecayLength = 100.0,
            default_lambda_MinCosPointingAngle = 0.990,

            default_lambda_InvariantMassMin = LambdaMass - 0.005, // 1.110683
            default_lambda_InvariantMassMax = LambdaMass + 0.005; // 1.120683

const float default_lambda_eta = 0.8,
            default_lambda_EtaDaughters = 0.8,
            default_lambda_TPCnclsDaughters = 80.0,
            default_lambda_MaxDcaV0Daughters = 1.5,

            default_lambda_PtProtonMin = 0.5,
            default_lambda_PtProtonMax = 4.0,

            default_lambda_PtPionMin = 0.5,
            default_lambda_PtPionMax = 4.0,

            default_lambda_MinToPrimVertexProton = 0.1,
            default_lambda_MinToPrimVertexPion = 0.3;


const int   default_lambda_StatusDaughters = AliESDtrack::kTPCrefit | AliESDtrack::kITSrefit,
            default_lambda_NdofDaughters = 4;

const bool  default_lambda_OnFlyStatus = kFALSE;

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

            default_event_VertZPosMin = -10.0,
            default_event_VertZPosMax = 10.0,

            default_event_EPVZEROMin = -1000.0,
            default_event_EPVZEROMax = 1000.0;

const  int  default_event_TriggerSelection = 0;
const  bool default_event_AcceptBadVertex = kFALSE;

const Bool_t  default_pair_TPCOnly = kTRUE;
const Float_t default_pair_TPCExitSepMin = -1.0;


const AliFemtoAnalysisPionLambda::PionType
  default_PionType = AliFemtoAnalysisPionLambda::kPiPlus;

const AliFemtoAnalysisPionLambda::LambdaType
  default_LambdaType = AliFemtoAnalysisPionLambda::kLambda;



void
AliFemtoAnalysisPionLambda::_Init(const AliFemtoAnalysisPionLambda::CutParams& p = AliFemtoAnalysisPionLambda::DefaultCutConfig())
{ // This function is to be called at the end of all constructors. It will
  // create default objects for any cuts that haven't been initialized

  if (fFirstParticleCut == NULL) {
    SetFirstParticleCut(BuildLambdaCut(p));
  }

  if (fSecondParticleCut == NULL) {
    SetSecondParticleCut(BuildPionCut(p));
  }

  if (fEventCut == NULL) {
    SetEventCut(BuildEventCut(p));
  }

  if (fPairCut == NULL) {
    SetPairCut(BuildPairCut(p));
  }
}

AliFemtoAnalysisPionLambda::AliFemtoAnalysisPionLambda():
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count, VertexBinning.min, VertexBinning.max,
                               MultBinning.bin_count,   MultBinning.min,   MultBinning.max),
  fAnalysisName("AliFemtoAnalysisPionLambda"),
  fPionType(default_PionType),
  fLambdaType(default_LambdaType),
  fGroupOutputObjects(default_group_output_objects)
{
  _Init();
}

AliFemtoAnalysisPionLambda::AliFemtoAnalysisPionLambda(const char *name):
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count, VertexBinning.min, VertexBinning.max,
                               MultBinning.bin_count,   MultBinning.min,   MultBinning.max),
  fAnalysisName(name),
  fPionType(default_PionType),
  fLambdaType(default_LambdaType),
  fGroupOutputObjects(default_group_output_objects)
{
  _Init();
}


AliFemtoAnalysisPionLambda::AliFemtoAnalysisPionLambda(const char *name,
                                                       PionType pion,
                                                       LambdaType lambda
):
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count, VertexBinning.min, VertexBinning.max,
                               MultBinning.bin_count,   MultBinning.min,   MultBinning.max),
  fAnalysisName(name),
  fPionType(pion),
  fLambdaType(lambda),
  fGroupOutputObjects(default_group_output_objects)
{
  _Init();
}

AliFemtoAnalysisPionLambda
  ::AliFemtoAnalysisPionLambda(const char *name,
                               const PionType pion,
                               const LambdaType lambda,
                               const CutParams& params):
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count, VertexBinning.min, VertexBinning.max,
                             MultBinning.bin_count,   MultBinning.min,   MultBinning.max),
  fAnalysisName(name),
  fPionType(pion),
  fLambdaType(lambda),
  fGroupOutputObjects(default_group_output_objects)
{
  _Init(params);
}

AliFemtoAnalysisPionLambda
  ::AliFemtoAnalysisPionLambda(const char *name,
                               const AnalysisParams &params,
                               const CutParams &cut_params):
  AliFemtoVertexMultAnalysis(params.vertex_bins, params.vertex_min, params.vertex_max,
                             params.mult_bins, params.mult_min, params.mult_max),
  fAnalysisName(name),
  fPionType(params.pion_type),
  fLambdaType(params.lambda_type),
  fGroupOutputObjects(params.group_output_objects)
{
//   cout << "[AliFemtoAnalysisPionLambda] Constructed with : " <<
//     params.vertex_bins <<
//     " " << params.vertex_min <<
//     " " << params.vertex_max <<
//     " " << params.mult_bins <<
//     " " << params.mult_min <<
//     " " << params.mult_max <<'\n';
  _Init(cut_params);
}


AliFemtoAnalysisPionLambda::AnalysisParams
AliFemtoAnalysisPionLambda::DefaultConfig()
{
  AliFemtoAnalysisPionLambda::AnalysisParams params = {
    VertexBinning.bin_count,
    VertexBinning.min,
    VertexBinning.max,

    MultBinning.bin_count,
    MultBinning.min,
    MultBinning.max,

    default_PionType,
    default_LambdaType,

    default_group_output_objects
  };

  return params;
}


AliFemtoAnalysisPionLambda::CutParams
AliFemtoAnalysisPionLambda::DefaultCutConfig()
{
  AliFemtoAnalysisPionLambda::CutParams params = {
    // Event
    default_event_EventMultMin,
    default_event_EventMultMax,
    default_event_VertZPosMin,
    default_event_VertZPosMax,
    default_event_EPVZEROMin,
    default_event_EPVZEROMax,
    default_event_TriggerSelection,
    default_event_AcceptBadVertex,

    // Pion
    default_pion_PtMin,
    default_pion_PtMax,
    default_pion_EtaMin,
    default_pion_EtaMax,
    default_pion_DCAMin,
    default_pion_DCAMax,

    default_pion_NSigmaMin,
    default_pion_NSigmaMax,

    // Lambda
    default_lambda_PtMin,
    default_lambda_PtMax,
    default_lambda_Eta,

    default_lambda_MaxDCAV0,
    default_lambda_MaxV0DecayLength,
    default_lambda_MinCosPointingAngle,

    default_lambda_InvariantMassMin,
    default_lambda_InvariantMassMax,

    default_lambda_OnFlyStatus,

    // daughter params
    default_lambda_EtaDaughters,
    default_lambda_TPCnclsDaughters,
    default_lambda_MaxDcaV0Daughters,

    default_lambda_StatusDaughters,
    default_lambda_NdofDaughters,

    // pion daughter
    default_lambda_PtPionMin,
    default_lambda_PtPionMax,
    default_lambda_MinToPrimVertexPion,

    // proton daughter
    default_lambda_PtProtonMin,
    default_lambda_PtProtonMax,
    default_lambda_MinToPrimVertexProton,

    // pair
    default_pair_TPCOnly,
    default_pair_TPCExitSepMin
  };

    
  assert(params.event_MultMin == default_event_EventMultMin);
  assert(params.pion_PtMin == default_pion_PtMin);
  assert(params.lambda_PtMin == default_lambda_PtMin);
  assert(params.lambda_daughter_StatusDaughters == default_lambda_StatusDaughters);
  assert(params.pair_TPCOnly == default_pair_TPCOnly);
  assert(params.pair_TPCExitSepMin == default_pair_TPCExitSepMin);
  return params;
}


AliFemtoV0TrackCut*
AliFemtoAnalysisPionLambda::BuildLambdaCut(const CutParams &p) const
{
  AliFemtoV0TrackCut* cut = new AliFemtoV0TrackCut();

  cut->SetMass(LambdaMass);
  cut->SetInvariantMassLambda(p.lambda_MassMin, p.lambda_MassMax);
  cut->SetPt(p.lambda_PtMin, p.lambda_PtMax);
  cut->SetEta(p.lambda_Eta);

  cut->SetOnFlyStatus(p.lambda_OnFlyStatus);
  cut->SetMaxDcaV0(p.lambda_MaxDCA);
  cut->SetMaxV0DecayLength(p.lambda_MaxDecayLength);
  cut->SetMinCosPointingAngle(p.lambda_MinCosPointingAngle);

  cut->SetEtaDaughters(p.lambda_daughter_Eta);
  cut->SetTPCnclsDaughters(p.lambda_daughter_TPCncls);
  cut->SetNdofDaughters(p.lambda_daughter_Ndof);
  cut->SetStatusDaughters(p.lambda_daughter_StatusDaughters);
  cut->SetMaxDcaV0Daughters(p.lambda_daughter_MaxDCA);

  switch (fLambdaType) {
  default:
    std::cout << "E-AliFemtoAnalysisPionLambda::BuildLambdaCut: Invalid lambda type: '" << fLambdaType << "'. Using kLambda.\n";
    const_cast<AliFemtoAnalysisPionLambda*>(this)->fLambdaType = kLambda;

  case kLambda:
    cut->SetParticleType(AliFemtoV0TrackCut::kLambda);
    cut->SetPtPosDaughter(p.lambda_daughter_proton_PtMin,
                          p.lambda_daughter_proton_PtMax);
    cut->SetPtNegDaughter(p.lambda_daughter_pion_PtMin,
                          p.lambda_daughter_pion_PtMax);
    cut->SetMinDaughtersToPrimVertex(p.lambda_daughter_pion_MinToPrimVertex,
                                     p.lambda_daughter_proton_MinToPrimVertex);
    break;

  case kAntiLambda:
    cut->SetParticleType(AliFemtoV0TrackCut::kAntiLambda);
    cut->SetPtPosDaughter(p.lambda_daughter_pion_PtMin,
                          p.lambda_daughter_pion_PtMax);
    cut->SetPtNegDaughter(p.lambda_daughter_proton_PtMin,
                          p.lambda_daughter_proton_PtMax);
    cut->SetMinDaughtersToPrimVertex(p.lambda_daughter_proton_MinToPrimVertex,
                                     p.lambda_daughter_pion_MinToPrimVertex);
    break;
  }

  return cut;
}

AliFemtoBasicTrackCut*
AliFemtoAnalysisPionLambda::BuildPionCut(const CutParams &p) const
{
  AliFemtoBasicTrackCut *cut = new AliFemtoBasicTrackCut();

  const int charge = (fPionType == kPiMinus) ? -1 :
                      (fPionType == kPiPlus) ? +1 :
                                              -999;
  if (charge == -999) {
    std::cerr << "E-AliFemtoAnalysisPionLambda::BuildPionCut: Invalid pion type: '" << fPionType << "'\n";
  }

  cut->SetCharge(charge);
  cut->SetMass(PionMass);

  cut->SetNSigmaPion(p.pion_NSigmaMin, p.pion_NSigmaMax);
  cut->SetPt(p.pion_PtMin, p.pion_PtMax);
  cut->SetRapidity(p.pion_EtaMin, p.pion_EtaMax);
  cut->SetDCA(p.pion_DCAMin, p.pion_DCAMax);

  return cut;
}

AliFemtoBasicEventCut*
AliFemtoAnalysisPionLambda::BuildEventCut(
  const AliFemtoAnalysisPionLambda::CutParams& p
) const
{
  AliFemtoBasicEventCut* cut = new AliFemtoBasicEventCut();

  cut->SetEventMult(p.event_MultMin,
                    p.event_MultMax);
  cut->SetVertZPos(p.event_VertexZMin,
                   p.event_VertexZMax);
  cut->SetEPVZERO(p.event_EP_VZeroMin,
                  p.event_EP_VZeroMax);
  cut->SetTriggerSelection(p.event_TriggerSelection);
  cut->SetAcceptBadVertex(p.event_AcceptBadVertex);

  return cut;
}


AliFemtoV0TrackPairCut*
AliFemtoAnalysisPionLambda::BuildPairCut(const CutParams &p) const
{
  AliFemtoV0TrackPairCut *cut = new AliFemtoV0TrackPairCut();

  cut->SetTPCOnly(p.pair_TPCOnly);
  cut->SetTPCExitSepMinimum(p.pair_TPCExitSepMin);

  return cut;
}

void AliFemtoAnalysisPionLambda::AddPionLambdaCutMonitor(AliFemtoCutMonitorPionLambda *monitor)
{

}


AliFemtoV0TrackCut* AliFemtoAnalysisPionLambda::GetLambdaCut()
{
  return static_cast<AliFemtoV0TrackCut*>(fFirstParticleCut);
}

AliFemtoTrackCut* AliFemtoAnalysisPionLambda::GetPionCut()
{
  return static_cast<AliFemtoTrackCut*>(fSecondParticleCut);
}

const AliFemtoV0TrackCut* AliFemtoAnalysisPionLambda::GetLambdaCut() const
{
  return static_cast<const AliFemtoV0TrackCut*>(fFirstParticleCut);
}

const AliFemtoTrackCut* AliFemtoAnalysisPionLambda::GetPionCut() const
{
  return static_cast<const AliFemtoTrackCut*>(fSecondParticleCut);
}


void AliFemtoAnalysisPionLambda::AddStanardCutMonitors()
{
  if (fEventCut) {
    fEventCut->AddCutMonitor(new AliFemtoPionLambdaCutMonitor::Event(true),
                             new AliFemtoPionLambdaCutMonitor::Event(false));
  } else {
    std::cout << " NO fEventCut!\n";
    exit(1);
  }

  if (fFirstParticleCut) {
    const TString ltype = (fLambdaType == kLambda) ? "Lambda" : "AntiLambda";
    fFirstParticleCut->AddCutMonitor(new AliFemtoPionLambdaCutMonitor::Lambda(true, ltype, fLambdaType),
                                     new AliFemtoPionLambdaCutMonitor::Lambda(false, ltype, fLambdaType));
  }

  if (fSecondParticleCut) {
    const TString ptype = (fPionType == kPiPlus) ? "Pi+" : "Pi-";
    fSecondParticleCut->AddCutMonitor(new AliFemtoPionLambdaCutMonitor::Pion(true, ptype),
                                      new AliFemtoPionLambdaCutMonitor::Pion(false, ptype));
  }

  if (fPairCut) {
    const TString ltype = (fLambdaType == kLambda) ? "Lam" : "ALam",
                  ptype = (fPionType == kPiPlus) ? "Pi+" : "Pi-",
                 pair_t = ptype + "/" + ltype;
    fPairCut->AddCutMonitor(new AliFemtoPionLambdaCutMonitor::Pair(true, pair_t),
                            new AliFemtoPionLambdaCutMonitor::Pair(false, pair_t));
  }

}



TList* AliFemtoAnalysisPionLambda::GetOutputList()
{
  TList *outputlist = NULL;
  TSeqCollection *output = NULL;
  // get "standard" outputs - that's what we output
  if (!fGroupOutputObjects) {
    outputlist = AliFemtoVertexMultAnalysis::GetOutputList();
    output = outputlist;
  } else {
    outputlist = new TList();
    output = new TObjArray();
    outputlist->Add(output);

    output->SetName(fAnalysisName);

    TList *super_objects = AliFemtoVertexMultAnalysis::GetOutputList();

    TIter next(super_objects);
    TObject *obj = NULL;
    while (obj = next()) {
      output->Add(obj);
    }

    delete super_objects;
  }

  // get the analysis report
  output->Add(new TObjString(TString("  REPORT  \n") + AliFemtoVertexMultAnalysis::Report()));

  // get the list of settings from each cut and correlation function
  {
    TObjArray *settings = new TObjArray();
    settings->SetName("settings");
    TList *setting_list = ListSettings();
    TIter next_setting(setting_list);
    TObject *setting = NULL;
    while (setting = next_setting()) {
      settings->Add(setting);
    }
    delete setting_list;
    settings->SetOwner(kTRUE);
    output->Add(settings);
 }

  return outputlist;
}

TList* AliFemtoAnalysisPionLambda::ListSettings()
{
  TList *setting_list = new TList();
  
  setting_list->Add(new TObjString(
    TString::Format("AliFemtoAnalysisPionLambda.piontype=%d", fPionType)
  ));
  
  setting_list->Add(new TObjString(
    TString::Format("AliFemtoAnalysisPionLambda.lambdatype=%d", fLambdaType)
  ));
  
  TList *parent_list = AliFemtoVertexMultAnalysis::ListSettings();
  
  setting_list->AddAll(parent_list);
  delete parent_list;
 
  
  return setting_list;
}

