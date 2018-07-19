///
/// \file AliFemtoAnalysisPionLambda.h
///

#pragma once

#ifndef ALIFEMTOANALYSIS_PION_LAMBDA_H_
#define ALIFEMTOANALYSIS_PION_LAMBDA_H_

class AliFemtoV0TrackCut;
class AliFemtoBasicTrackCut;
class AliFemtoV0TrackPairCut;
class AliFemtoBasicEventCut;

class AliFemtoCutMonitorPionLambda;

class TList;

#include "AliFemtoVertexMultAnalysis.h"
#include "AliFemtoV0TrackPairCut.h"


/// \class AliFemtoAnalysisPionLambda
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///
/// \brief A simple analysis for studying femtoscopy between charged
///        pions and lambda baryons.
///
/// This class reduces boilerplate code required to cut on specific
/// variations of pions and lambdas.
///
/// The 'first' and 'second' particle cuts are aliased to 'lambda'
/// and 'pion' cuts with specific types :ref:AliFemtoV0TrackCut and
/// :ref:AliFemtoTrackCut, instead of the general AliFemtoParticleCut.
/// These are accessable by the GetLambdaCut and GetPionCut.
///
/// To remove the hassle of managing which daughter particle should
/// be the proton and which the pion when looking at lambdas vs
/// antilambdas, aliases for daughter cut parameters have been added.
/// These rely on the "type" of the lambda and pion being specified
/// in the analysis' constructor via the enumerated particle types.
///
class AliFemtoAnalysisPionLambda : public AliFemtoVertexMultAnalysis {
public:
  enum PionType {kPiPlus, kPiMinus};
  enum LambdaType {kLambda, kAntiLambda};

  struct CutParams;
  struct AnalysisParams;

public:
  /// default constructor
  AliFemtoAnalysisPionLambda();

  /// Construct the analysis with only a name.
  /// The particle types default to PiPlus and Lambda
  ///
  AliFemtoAnalysisPionLambda(const char *name);

  /// Construct the analysis with a name and specify the pion and
  /// lambda types to analyze.
  ///
  AliFemtoAnalysisPionLambda(const char *name,
                             const PionType pion,
                             const LambdaType lambda);

  /// Construct the analysis with a name, the pion and lambda types,
  /// and cut-configuration object.
  ///
  AliFemtoAnalysisPionLambda(const char *name,
                             const PionType pion,
                             const LambdaType lambda,
                             const CutParams& params);


  /// Construct the analysis with a name, the pion and lambda types,
  /// and cut-configuration object.
  ///
  AliFemtoAnalysisPionLambda(const char *name,
                             const AnalysisParams&,
                             const CutParams&);

  virtual void EventBegin(const AliFemtoEvent*);
  virtual void EventEnd(const AliFemtoEvent*);

  /// Create a Cut-Parameter object with default options.
  static CutParams DefaultCutConfig();

  /// Create an analysis-parameter object with default options.
  /// This is used internally when using the default constructor.
  /// It is recommended to call this method to construct the
  /// AnalysisParams object, and modify that to create a custom
  /// AliFemtoAnalysisPionLambda.
  ///
  static AnalysisParams DefaultConfig();
  static AnalysisParams DefaultConfigWithTypes(PionType, LambdaType);

  /// Return the particle cut which selects the lambda-like particle
  /// (implemented as ParticleCut1 - but could change)
  AliFemtoV0TrackCut* GetLambdaCut();
  const AliFemtoV0TrackCut* GetLambdaCut() const;

  /// Return the particle cut which selects the charged pion-like
  /// track (implemented as ParticleCut2 - but could change)
  AliFemtoTrackCut* GetPionCut();
  const AliFemtoTrackCut* GetPionCut() const;

  AliFemtoV0TrackCut* BuildLambdaCut(const CutParams&) const;
  AliFemtoBasicTrackCut* BuildPionCut(const CutParams&) const;
  AliFemtoBasicEventCut* BuildEventCut(const CutParams&) const;
  AliFemtoV0TrackPairCut* BuildPairCut(const CutParams&) const;

  /// For each cut already set (i.e. pointer is not NULL) this
  /// function creates two "standard" cut monitor (for pass/fail)
  /// and adds them to the cut.
  /// These cut monitors are defined with the
  /// AliFemtoPionLambdaCutMonitor class.
  ///
  void AddStanardCutMonitors();

  /// Returns a TList of all objects
  virtual TList *GetOutputList();

  /// Returns a TList of TObjStrings containing the settings of the
  /// analysis.
  virtual TList *ListSettings();

protected:

  /// The name of the analysis to identify in the output list
  TString fAnalysisName;

  /// The type of Pion particles this analysis will search for
  PionType fPionType;

  /// The type of Lambda particles this analysis will search for.
  /// This determines which V0 daughter (pos or neg) will pass proton
  /// vs pion daughter cuts
  LambdaType fLambdaType;

  /// Output objects in multiple TObjArrays instead of one TList
  Bool_t fGroupOutputObjects;

  /// This is a Monte Carlo analysis
  Bool_t fMCAnalysis;

};

/// \class AliFemtoAnalysisPionLambda::AnalysisParams
struct AliFemtoAnalysisPionLambda::AnalysisParams {
  UInt_t vertex_bins;
  Float_t vertex_min;
  Float_t vertex_max;

  UInt_t mult_bins;
  Float_t mult_min;
  Float_t mult_max;

  PionType pion_type;
  LambdaType lambda_type;

  UInt_t num_events_to_mix;
  UInt_t min_coll_size;

  Bool_t verbose;
  Bool_t enable_pair_monitors;
  Bool_t group_output_objects;
  Bool_t is_mc_analysis;
};

/// \class AliFemtoAnalysisPionLambda::CutParams
/// \brief *not yet implemented* Structure containing all fit
///        parameters for 'easy' setting of custom fit parameters
///        via one command.
///
struct AliFemtoAnalysisPionLambda::CutParams {
  // EVENT
  Float_t event_MultMin,
          event_MultMax;

  Float_t event_VertexZMin,
          event_VertexZMax;

  Float_t event_EP_VZeroMin,
          event_EP_VZeroMax;

  Int_t   event_TriggerSelection;
  Bool_t  event_AcceptBadVertex;

  // PION
  Float_t pion_PtMin,
          pion_PtMax;

  Float_t pion_EtaMin,
          pion_EtaMax;

  Float_t pion_DCAMin,
          pion_DCAMax;

  Float_t pion_NSigmaMin,
          pion_NSigmaMax;

  //  LAMBDA
  Float_t lambda_PtMin,
          lambda_PtMax,
          lambda_Eta,
          lambda_MaxDCA,
          lambda_MaxDecayLength,
          lambda_MinCosPointingAngle,

          lambda_MassMin,
          lambda_MassMax;

  Bool_t lambda_OnFlyStatus;

    // Lambda daughter parameters
    Float_t lambda_daughter_Eta,
            lambda_daughter_TPCncls,
            lambda_daughter_MaxDCA;

    Int_t lambda_daughter_StatusDaughters,
          lambda_daughter_Ndof;


      Float_t lambda_daughter_pion_PtMin,
              lambda_daughter_pion_PtMax,
              lambda_daughter_pion_MinToPrimVertex;

      Float_t lambda_daughter_proton_PtMin,
              lambda_daughter_proton_PtMax,
              lambda_daughter_proton_MinToPrimVertex;

  // PAIR
  Bool_t pair_TPCOnly;
  Float_t pair_TPCExitSepMin;
  Float_t pair_MinAvgSeparationPos;
  Float_t pair_MinAvgSeparationNeg;
};

#endif
