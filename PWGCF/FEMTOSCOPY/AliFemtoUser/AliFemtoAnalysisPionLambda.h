///
/// \file AliFemtoAnalysisPionLambda.h
///

#pragma once

#ifndef ALIFEMTOANALYSIS_PION_LAMBDA_H_
#define ALIFEMTOANALYSIS_PION_LAMBDA_H_

#include "AliFemtoVertexMultAnalysis.h"

class AliFemtoV0TrackCut;
class AliFemtoBasicTrackCut;
class AliFemtoV0TrackPairCut;
class AliFemtoBasicEventCut;

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
                             PionType pion,
                             LambdaType lambda);


  AliFemtoV0TrackCut* GetLambdaCut();
  AliFemtoTrackCut* GetPionCut();

  const AliFemtoV0TrackCut* GetLambdaCut() const;
  const AliFemtoTrackCut* GetPionCut() const;

  AliFemtoV0TrackCut* BuildDefaultLambdaCut();
  AliFemtoBasicTrackCut* BuildDefaultPionCut();
  AliFemtoBasicEventCut* BuildDefaultEventCut();
  AliFemtoV0TrackPairCut* BuildDefaultPairCut();



protected:

  PionType fPionType;
  LambdaType fLambdaType;

private:

  void _Init();

#ifndef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoAnalysisPionLambda, 0);
  /// \endcond
#endif

};

/// \class AliFemtoAnalysisPionLambda::CutParams
/// \brief *not yet implemented* Structure containing all fit parameters for
///        'easy' setting of custom fit parameters via one command.
///
struct AliFemtoAnalysisPionLambda::CutParams {
  struct {
    float PtMin;
    float PtMax;
  } lambda;

  struct {
    float PtMin;
    float PtMax;
  } pion;
};





#endif
