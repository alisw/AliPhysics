///
/// \file AliFemtoAnalysisPionLambda.h
///


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
/// \brief A simple analysis for studying femtoscopy between charged pions and
///        lambda baryons.
///
/// This class reduces boilerplate code required to cut on specific variations
/// of pions and lambdas.
///
/// The 'first' and 'second' particle cuts are aliased to 'lambda' and 'pion'
/// cuts with specific types :ref:AliFemtoV0TrackCut and :ref:AliFemtoTrackCut,
/// instead of the general AliFemtoParticleCut. These are accessable by the
/// GetLambdaCut and GetPionCut.
///
/// To remove the hassle of managing which daughter particle should be the
/// proton and which the pion when looking at lambdas vs antilambdas, aliases
/// for daughter cut parameters have been added. These rely on the "type" of
/// the lambda and pion being specified in the analysis' constructor via the
/// enumerated particle types.
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

  /**
   * Construct the analysis with only a name.
   * The particle types default to PiPlus and Lambda
   */
  AliFemtoAnalysisPionLambda(const char *name);

  /**
   * Construct the analysis with a name and specify the pion and lambda types
   * to analyze.
   */
  AliFemtoAnalysisPionLambda(const char *name,
                             const PionType pion,
                             const LambdaType lambda);

  /**
   * Construct the analysis with a name, the pion and lambda types, and
   * cut-configuration object.
   */
  AliFemtoAnalysisPionLambda(const char *name,
                             const PionType pion,
                             const LambdaType lambda,
                             const CutParams& params);


  /**
   * Construct the analysis with a name, the pion and lambda types, and
   * cut-configuration object.
   */
  AliFemtoAnalysisPionLambda(const char *name,
                             const AnalysisParams&,
                             const CutParams&);

  /**
   * Create a Cut-Parameter object with default options.
   */
  static CutParams DefaultCutConfig();
  static AnalysisParams DefaultConfig();

  AliFemtoV0TrackCut* GetLambdaCut();
  AliFemtoTrackCut* GetPionCut();

  const AliFemtoV0TrackCut* GetLambdaCut() const;
  const AliFemtoTrackCut* GetPionCut() const;

  AliFemtoV0TrackCut* BuildLambdaCut(const CutParams&) const;
  AliFemtoBasicTrackCut* BuildPionCut(const CutParams&) const ;
  AliFemtoBasicEventCut* BuildEventCut(const CutParams&) const;
  AliFemtoV0TrackPairCut* BuildPairCut(const CutParams&) const;

  // Adds a AliFemtoCutMonitorPionLambda to the cuts
  void AddPionLambdaCutMonitor(AliFemtoCutMonitorPionLambda*);

  void AddStanardCutMonitors();

  virtual TList *GetOutputList();

  virtual TList *ListSettings();

protected:

  /// The name of the analysis to identify in the output list
  TString fAnalysisName;
  
  /// The type of Pion particles this analysis will search for
  PionType fPionType;
  
  /// The type of Lambda particles this analysis will search for.
  /// This determines which V0 daughter (pos or neg) will pass proton vs pion daughter cuts
  LambdaType fLambdaType;

  Bool_t fGroupOutputObjects;

private:

  void _Init(const CutParams& params);

#ifndef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoAnalysisPionLambda, 0);
  /// \endcond
#endif

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

  Bool_t group_output_objects;
};

/// \class AliFemtoAnalysisPionLambda::CutParams
/// \brief *not yet implemented* Structure containing all fit parameters for
///        'easy' setting of custom fit parameters via one command.
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

    // Lambda daugther parameters
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
};





#endif
