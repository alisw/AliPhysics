///
/// \file AliFemtoAnalysisPionPion.h
///


#ifndef ALIFEMTOANALYSIS_PIONPION_H_
#define ALIFEMTOANALYSIS_PIONPION_H_

class AliFemtoBasicTrackCut;
class AliFemtoPairCutAntiGamma;
class AliFemtoBasicEventCut;
class AliFemtoPairCutDetaDphi;

class TList;


#include "AliFemtoVertexMultAnalysis.h"


/// \class AliFemtoAnalysisPionPion
/// \brief A simple analysis for studying femtoscopic systems of charged pions
///
/// This class reduces boilerplate code required to cut on specific variations
/// of pion-pion femtoscopy.
///
/// \author Andrew Kubera, Ohio State University <andrew.kubera@cern.ch>
///
class AliFemtoAnalysisPionPion : public AliFemtoVertexMultAnalysis {
public:
  /// Used to name pion types
  enum PionType {
    /// Special value to allow identical particle femtoscopy (the second type is none)
    kNone = -1,

    /// identifies pi+ tracks
    kPiPlus,

    /// identifies pi- tracks
    kPiMinus
  };

  struct CutParams;
  struct AnalysisParams;

public:
  /// default constructor
  AliFemtoAnalysisPionPion();

  /// Construct the analysis with only a name.
  /// The particle types default to PiPlus
  ///
  AliFemtoAnalysisPionPion(const char *name);

  /// Construct the analysis with a name and specify the pion types
  /// to analyze.
  ///
  AliFemtoAnalysisPionPion(const char *name,
                           const PionType pion_1,
                           const PionType pion_2);

  /// Construct the analysis with a name, the pion types, and specialized
  /// cut-configuration object.
  ///
  AliFemtoAnalysisPionPion(const char *name,
                           const PionType pion_1,
                           const PionType pion_2,
                           const CutParams& params);


  /// Construct the analysis with a name, analysis parameters, and
  /// cut-configuration object.
  ///
  AliFemtoAnalysisPionPion(const char *name,
                           const AnalysisParams&,
                           const CutParams&);

  virtual void EventBegin(const AliFemtoEvent*);
  virtual void EventEnd(const AliFemtoEvent*);

  /// Create a Cut-Parameter object with default options.
  ///
  static CutParams DefaultCutConfig();

  /**
   * Create an analysis-parameter object with default options. This is used
   * internally when using the default constructor. It is recommended to call
   * this method to constrcut the AnalysisParams object, and modify that to
   * create a custom AliFemtoAnalysisPionPion.
   */
  static AnalysisParams DefaultConfig();

  AliFemtoTrackCut* BuildPionCut1(const CutParams&) const;
  AliFemtoTrackCut* BuildPionCut2(const CutParams&) const;
  AliFemtoEventCut* BuildEventCut(const CutParams&) const;
  AliFemtoPairCut* BuildPairCut(const CutParams&) const;

  /**
   * For each cut already set (i.e. pointer is not NULL) this function
   * creates two "standard" cut monitor (for pass/fail) and adds them to the
   * cut. These cut monitors are found within the AliFemtoCutMonitorPionPion
   * nam.
   */
  void AddStanardCutMonitors();

  /**
   * Returns a TList of all objects
   */
  virtual TList* GetOutputList();

  /**
   * Returns a TList of TObjStrings containing the settings of the analysis.
   */
  virtual TList* ListSettings();

protected:

  /// The name of the analysis to identify in the output list
  TString fAnalysisName;

  /// The type of Pion particles this analysis will search for
  PionType fPionType_1;

  /// The type of Lambda particles this analysis will search for.
  /// This determines which V0 daughter (pos or neg) will pass proton vs pion daughter cuts
  PionType fPionType_2;

  /// Output objects in multiple TObjArrays instead of one TList
  Bool_t fGroupOutputObjects;

  /// This is a Monte Carlo analysis
  Bool_t fMCAnalysis;

};

/// \class AliFemtoAnalysisPionPion::AnalysisParams
struct AliFemtoAnalysisPionPion::AnalysisParams {

  UInt_t vertex_bins;
  Float_t vertex_min,
          vertex_max;

  UInt_t mult_bins;
  Float_t mult_min,
          mult_max;

  PionType pion_type_1;
  PionType pion_type_2;

  UInt_t num_events_to_mix;
  UInt_t min_coll_size;

  Bool_t verbose;
  Bool_t enable_pair_monitors;
  Bool_t group_output_objects;
  Bool_t is_mc_analysis;
/*
  AnalysisParams():
    vertex_bins(8)
  , vertex_min(-8)
  , vertex_max(8)

  , mult_bins(4)
  , mult_min(0)
  , mult_max(10000)

  , pion_type_1(kNone)
  , pion_type_2(kNone)

  , num_events_to_mix(10)
  , min_coll_size(100)
  , enable_pair_monitors(kTRUE)
  , group_output_objects(kTRUE)
  , is_mc_analysis(kFALSE)
  {}
*/
};

/// \class AliFemtoAnalysisPionPion::CutParams
/// \brief Structure containing all cut parameters for use with
///        constructing an AliFemtoAnalysisPionPion object.
///
/// The expected way to use this class
///
struct AliFemtoAnalysisPionPion::CutParams {
  // EVENT
  Float_t event_MultMin,
          event_MultMax;

  Float_t event_CentralityMin,
          event_CentralityMax;

  Float_t event_VertexZMin,
          event_VertexZMax;

  Float_t event_EP_VZeroMin,
          event_EP_VZeroMax;

  Int_t   event_TriggerSelection;
  Bool_t  event_AcceptBadVertex;


  // PION - 1
  Float_t pion_1_PtMin,
          pion_1_PtMax;

  Float_t pion_1_EtaMin,
          pion_1_EtaMax;

  Float_t pion_1_DCAMin,
          pion_1_DCAMax;

  Float_t pion_1_NSigmaMin
        , pion_1_NSigmaMax
        ;

  Float_t pion_1_max_impact_xy
        , pion_1_max_impact_z
        , pion_1_max_tpc_chi_ndof
        , pion_1_max_its_chi_ndof
        ;

  UInt_t pion_1_min_tpc_ncls;
  Bool_t pion_1_remove_kinks,
         pion_1_set_label;


  // PION - 2
  Float_t pion_2_PtMin,
          pion_2_PtMax;

  Float_t pion_2_EtaMin,
          pion_2_EtaMax;

  Float_t pion_2_DCAMin,
          pion_2_DCAMax;

  Float_t pion_2_NSigmaMin,
          pion_2_NSigmaMax;

  Float_t pion_2_max_impact_xy
        , pion_2_max_impact_z
        , pion_2_max_tpc_chi_ndof
        , pion_2_max_its_chi_ndof
        ;

  UInt_t pion_2_min_tpc_ncls;
  Bool_t pion_2_remove_kinks,
         pion_2_set_label;


  // PAIR
  Bool_t pair_TPCOnly;
  // Float_t pair_TPCExitSepMin;
  // Float_t pair_MinAvgSeparationPos;
  // Float_t pair_MinAvgSeparationNeg;

  Float_t pair_delta_eta_min
        , pair_delta_phi_min
        , pair_phi_star_radius
        ;

  Float_t pair_max_share_quality,
          pair_max_share_fraction;
  Bool_t pair_remove_same_label;

};

#endif
