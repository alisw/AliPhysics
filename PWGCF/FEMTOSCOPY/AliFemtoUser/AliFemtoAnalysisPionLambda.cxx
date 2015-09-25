///
/// \file AliFemtoAnalysisPionLambda.cxx
///

#include "AliFemtoAnalysisPionLambda.h"

#include "AliESDtrack.h"

#include "AliFemtoV0TrackCut.h"
#include "AliFemtoBasicTrackCut.h"
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoV0TrackPairCut.h"

static const double LambdaMass = 1.115683,
                    PionMass = 0.139;

const struct { unsigned int bin_count; float min; float max; }

VertexBinning = {16, -10.0, 10.0},
MultBinning = {30, 0, 10000};


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
            default_lambda_PtProtonMax = 0.4,

            default_lambda_PtPionMin = 0.5,
            default_lambda_PtPionMax = 4.0,

            default_lambda_MinToPrimVertexProton = 0.1,
            default_lambda_MinToPrimVertexPion = 0.3;


const int   default_lambda_StatusDaughters = AliESDtrack::kTPCrefit,
            default_lambda_NdofDaughters = 4;

const bool  default_lambda_OnFlyStatus = kFALSE;

const float default_pion_PtMin = 0.2,
            default_pion_PtMax = 2.0,
            default_pion_EtaMin = 0.8,
            default_pion_EtaMax = -0.8,
            default_pion_DCAMin = 0.5,
            default_pion_DCAMax = 4.0;

void
AliFemtoAnalysisPionLambda::_Init()
{ // This function is to be called at the end of all constructors. It will
  // create default objects for any cuts that haven't been initialized

  if (fFirstParticleCut == NULL) {
    SetFirstParticleCut((AliFemtoParticleCut*)BuildDefaultLambdaCut());
  }

  if (fSecondParticleCut == NULL) {
    SetSecondParticleCut((AliFemtoParticleCut*)BuildDefaultPionCut());
  }

  if (fEventCut == NULL) {
    SetEventCut((AliFemtoEventCut*)BuildDefaultEventCut());
  }

  if (fPairCut == NULL) {
    SetPairCut(BuildDefaultPairCut());
  }
}

AliFemtoAnalysisPionLambda::AliFemtoAnalysisPionLambda():
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count, VertexBinning.min, VertexBinning.max,
                               MultBinning.bin_count,   MultBinning.min,   MultBinning.max),
  fPionType(kPiPlus),
  fLambdaType(kLambda)
{
  _Init();
}

AliFemtoAnalysisPionLambda::AliFemtoAnalysisPionLambda(const char *name):
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count, VertexBinning.min, VertexBinning.max,
                               MultBinning.bin_count,   MultBinning.min,   MultBinning.max),
  fPionType(kPiPlus),
  fLambdaType(kLambda)
{
  _Init();
}


AliFemtoAnalysisPionLambda::AliFemtoAnalysisPionLambda(const char *name,
                                                       PionType pion,
                                                       LambdaType lambda
):
  AliFemtoVertexMultAnalysis(VertexBinning.bin_count, VertexBinning.min, VertexBinning.max,
                               MultBinning.bin_count,   MultBinning.min,   MultBinning.max),
  fPionType(pion),
  fLambdaType(lambda)
{
  _Init();
}

AliFemtoV0TrackCut*
AliFemtoAnalysisPionLambda::BuildDefaultLambdaCut()
{
  AliFemtoV0TrackCut* cut = new AliFemtoV0TrackCut();
  cut->SetMass(LambdaMass);
  cut->SetPt(default_lambda_PtMin, default_lambda_PtMax);
  cut->SetEta(default_lambda_Eta);

  cut->SetEtaDaughters(default_lambda_EtaDaughters);
  cut->SetTPCnclsDaughters(default_lambda_EtaDaughters);
  cut->SetNdofDaughters(default_lambda_NdofDaughters);
  cut->SetStatusDaughters(default_lambda_StatusDaughters);
  cut->SetMaxDcaV0Daughters(default_lambda_MaxDcaV0Daughters);

  cut->SetOnFlyStatus(default_lambda_OnFlyStatus);
  cut->SetMaxDcaV0(default_lambda_MaxDCAV0);
  cut->SetMaxV0DecayLength(default_lambda_MaxV0DecayLength);
  cut->SetMinCosPointingAngle(default_lambda_MinCosPointingAngle);
  cut->SetInvariantMassLambda(default_lambda_InvariantMassMin,
                              default_lambda_InvariantMassMax);

  switch (fLambdaType) {
  case kLambda:
    cut->SetParticleType(AliFemtoV0TrackCut::kLambda);
    cut->SetPtPosDaughter(default_lambda_PtProtonMin,
                          default_lambda_PtProtonMax);
    cut->SetPtNegDaughter(default_lambda_PtPionMin,
                          default_lambda_PtPionMax);
    cut->SetMinDaughtersToPrimVertex(default_lambda_MinToPrimVertexPion,
                                     default_lambda_MinToPrimVertexProton);

  case kAntiLambda:
    cut->SetParticleType(AliFemtoV0TrackCut::kAntiLambda);
    cut->SetPtPosDaughter(default_lambda_PtPionMin,
                          default_lambda_PtPionMax);
    cut->SetPtNegDaughter(default_lambda_PtProtonMin,
                          default_lambda_PtProtonMax);
    cut->SetMinDaughtersToPrimVertex(default_lambda_MinToPrimVertexPion,
                                     default_lambda_MinToPrimVertexProton);
  }

  return cut;
}

const float default_pion_NSigmaPionMin = -3.0,
            default_pion_NSigmaPionMax = 3.0;

AliFemtoBasicTrackCut*
AliFemtoAnalysisPionLambda::BuildDefaultPionCut()
{
  AliFemtoBasicTrackCut *cut = new AliFemtoBasicTrackCut();
  const int charge = (fPionType == kPiMinus) ? -1 :
                      (fPionType == kPiPlus) ? +1 :
                                              -999;
  if (charge == -999) {
    std::cerr << "E-AliFemtoAnalysisPionLambda::BuildDefaultPionCut: Invalid pion type: '" << fPionType << "'\n";
  }

  cut->SetCharge(charge);
  cut->SetMass(PionMass);
  cut->SetNSigmaPion(default_pion_NSigmaPionMin, default_pion_NSigmaPionMax);
  cut->SetPt(default_pion_PtMin, default_pion_PtMax);
  cut->SetRapidity(default_pion_EtaMin, default_pion_EtaMax);
  cut->SetDCA(default_pion_DCAMin, default_pion_DCAMax);

  return cut;
}


const float default_event_EventMultMin = 0,
            default_event_EventMultMax = 100000,

            default_event_VertZPosMin = -10.0,
            default_event_VertZPosMax = 10.0,

            default_event_EPVZEROMin = -1000.0,
            default_event_EPVZEROMax = 1000.0;

const  int  default_event_TriggerSelection = 0;
const  bool default_event_AcceptBadVertex = kFALSE;

AliFemtoBasicEventCut*
AliFemtoAnalysisPionLambda::BuildDefaultEventCut()
{
  AliFemtoBasicEventCut* cut = new AliFemtoBasicEventCut();
  cut->SetEventMult(default_event_EventMultMin,
                    default_event_EventMultMax);
  cut->SetVertZPos(default_event_VertZPosMin,
                   default_event_VertZPosMax);
  cut->SetEPVZERO(default_event_EPVZEROMin,
                  default_event_EPVZEROMax);
  cut->SetTriggerSelection(default_event_TriggerSelection);
  cut->SetAcceptBadVertex(default_event_AcceptBadVertex);

  return cut;
}


const Bool_t default_pair_TPCOnly = kTRUE;
const Float_t default_pair_TPCExitSepMin = -1.0;

AliFemtoV0TrackPairCut*
AliFemtoAnalysisPionLambda::BuildDefaultPairCut()
{
  AliFemtoV0TrackPairCut* cut = new AliFemtoV0TrackPairCut();
  cut->SetTPCOnly(default_pair_TPCOnly);
  cut->SetTPCExitSepMinimum(default_pair_TPCExitSepMin);
  return cut;
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
