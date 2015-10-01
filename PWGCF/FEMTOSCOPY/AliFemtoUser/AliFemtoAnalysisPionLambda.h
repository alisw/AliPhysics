///
/// \file AliFemtoAnalysisPionLambda.h
///


#ifndef ALIFEMTOANALYSIS_PION_LAMBDA_H_
#define ALIFEMTOANALYSIS_PION_LAMBDA_H_

#include "AliFemtoVertexMultAnalysis.h"

class AliFemtoV0TrackCut;
class AliFemtoBasicTrackCut;
class AliFemtoV0TrackPairCut;
class AliFemtoBasicEventCut;



struct AnalysisParams_Type {
  struct {
    UInt_t bins;
    Float_t min;
    Float_t max;
  } vertex;

  float pion_type;
  float lambda_type;
};



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
/// cuts with specific types :ref:AliFemtoV0TrackCut and :ref:AliFemtoTrackCut, instead
/// of the general AliFemtoParticleCut. These are accessable by the
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


protected:

  PionType fPionType;
  LambdaType fLambdaType;

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
};

/// \class AliFemtoAnalysisPionLambda::CutParams
/// \brief *not yet implemented* Structure containing all fit parameters for
///        'easy' setting of custom fit parameters via one command.
///
struct AliFemtoAnalysisPionLambda::CutParams {
  // EVENT
  struct {
    Float_t MultMin,
            MultMax;
    Float_t VertexZMin,
            VertexZMax;
    Float_t EP_VZeroMin,
            EP_VZeroMax;
    Int_t   TriggerSelection;
    Bool_t  AcceptBadVertex;
  } event;

  // PION
  struct {
    Float_t PtMin,
            PtMax;

    Float_t EtaMin,
            EtaMax;

    Float_t DCAMin,
            DCAMax;

    Float_t NSigmaPionMin,
            NSigmaPionMax;
  } pion;

  //  LAMBDA
  struct {
    Float_t PtMin,
            PtMax,
            Eta,
            MaxDCA,
            MaxDecayLength,
            MinCosPointingAngle;

    Float_t MassMin,
            MassMax;

    Bool_t OnFlyStatus;

    /// Lambda daugther parameters
    struct {

      Float_t Eta,
              TPCncls,
              MaxDCA;

      Int_t StatusDaughters,
            Ndof;

      struct {
        Float_t PtMin,
                PtMax,
                MinToPrimVertex;
      } pion,
        proton;

    } daughter;

  } lambda;

  struct {
    Bool_t TPCOnly;
    Float_t TPCExitSepMin;
  } pair;

};





#endif
