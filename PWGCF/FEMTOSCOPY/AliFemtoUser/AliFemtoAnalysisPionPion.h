///
/// \file AliFemtoAnalysisPionPion.h
///


#ifndef ALIFEMTOANALYSIS_PIONPION_H_
#define ALIFEMTOANALYSIS_PIONPION_H_

class AliFemtoEventReader;
class TList;



#include "AliFemtoConfigObject.h"

#include "AliFemtoVertexMultAnalysis.h"


#include <TClass.h>
#include <TNamed.h>


/// \class AliFemtoAnalysisPionPion
/// \brief A simple analysis for studying femtoscopic systems of
///        charged pions
///
/// This class reduces boilerplate code required to cut on specific
/// variations of pion-pion femtoscopy.
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

  /// Construct from configuration object
  ///
  AliFemtoAnalysisPionPion(const AliFemtoConfigObject &);

  virtual void EventBegin(const AliFemtoEvent*);
  virtual void EventEnd(const AliFemtoEvent*);

  /// Create a Cut-Parameter object with default options.
  static CutParams DefaultCutConfig();

  /// Create an analysis-parameter object with default options.
  /// This is used internally when using the default constructor.
  /// It is recommended to call this method to constrcut the
  /// AnalysisParams object, and modify that to create a custom
  /// AliFemtoAnalysisPionPion.
  static AnalysisParams DefaultConfig();

  AliFemtoTrackCut* BuildPionCut1(const CutParams&) const;
  AliFemtoTrackCut* BuildPionCut2(const CutParams&) const;
  AliFemtoEventCut* BuildEventCut(const CutParams&) const;
  AliFemtoPairCut* BuildPairCut(const CutParams&) const;

  /// For each cut already set (i.e. pointer is not NULL) this
  /// function creates two "standard" cut monitor (for pass/fail) and
  /// adds them to the cut.
  /// These cut monitors are found within the AliFemtoCutMonitorPionPion
  /// namespace.
  void AddStanardCutMonitors();

  /// Store configuration of eventreader
  void StoreEventReaderConfiguration(const AliFemtoEventReader &);

  /// Return configuration object describing analysis
  virtual AliFemtoConfigObject GetConfiguration() const;

  /// Returns a TList of all objects
  virtual TList* GetOutputList();

  /// Returns a TList of TObjStrings containing the settings of the
  /// analysis.
  virtual TList* ListSettings();

  /// Construct analysis based on the input from config object
  ///
  /// This static method is only using the class as a namespace (this
  /// should be changed eventually)
  ///
  static AliFemtoAnalysis* BuildAnalysisFromConfiguration(AliFemtoConfigObject);

  /// Construct an Event Reader from config object
  ///
  /// This static method is only using the class as a namespace (this
  /// should be changed eventually)
  ///
  static AliFemtoEventReader* ConstructEventReader(AliFemtoConfigObject cfg);

  /// Construct event cut via configuration
  static AliFemtoEventCut* ConstructEventCut(AliFemtoConfigObject);

  /// Construct pion cut via configuration
  static AliFemtoParticleCut* ConstructParticleCut(AliFemtoConfigObject);

  /// Construct pair cut via configuration
  static AliFemtoPairCut* ConstructPairCut(AliFemtoConfigObject);

  /// Construct a CorrelationFunction from config object
  static AliFemtoCorrFctn* ConstructCorrelationFunction(AliFemtoConfigObject);

  template <typename T>
  static AliFemtoConfigObject GetConfigurationOf(const T&);

  static TString make_random_string(const TString &prefix="");

protected:

  /// The name of this analysis used for identification in the output list
  TString fAnalysisName;

  /// The type of Pion particles this analysis will search for
  PionType fPionType_1;

  /// The 2nd type of pion particles this analysis will search for
  /// If kNone, this is an identical-pion analysis.
  PionType fPionType_2;

  /// Output objects in multiple TObjArrays instead of one TList
  Bool_t fGroupOutputObjects;

  /// If false, will not including settings string in the output list
  Bool_t fOutputSettings;

  /// This is a Monte Carlo analysis
  Bool_t fMCAnalysis;

  /// Configuration of event reader
  AliFemtoConfigObject fEventReaderCfg;
};

/// \class AliFemtoAnalysisPionPion::AnalysisParams
struct AliFemtoAnalysisPionPion::AnalysisParams : public TNamed {

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
  Bool_t output_settings;
  Bool_t is_mc_analysis;

  // Bool_t auto_mult_bin;

  /// get multiplicty from cut params
  void calc_automult(const AliFemtoAnalysisPionPion::CutParams &);

  /// Default Values
  AnalysisParams();

  ClassDef(AnalysisParams, 1);
};

/// \class AliFemtoAnalysisPionPion::CutParams
/// \brief Structure containing all cut parameters for use with
///        constructing an AliFemtoAnalysisPionPion object.
///
/// The expected way to use this class
///
struct AliFemtoAnalysisPionPion::CutParams : public TNamed {

  Bool_t cuts_use_attrs;

  Bool_t mc_pion_only;
  Bool_t mc_nonpion_only;

  Bool_t event_use_basic;

  // EVENT
  std::pair<int, int> event_mult;

  std::pair<double, double> event_centrality,
                            event_vertex_z,
                            event_ep_vzero;

  Int_t event_TriggerSelection;
  Bool_t event_AcceptBadVertex;
  Bool_t event_AcceptOnlyPhysics;
  Int_t event_zdc_part;

  // PION - 1
  std::pair<double, double> pion_1_pt,
                            pion_1_eta,
                            pion_1_DCA;

  ULong_t pion_1_status;

  Float_t pion_1_tpc_sigma;
  Float_t pion_1_tof_sigma;
  Bool_t pion_1_tof_required;

  Float_t pion_1_tof_kreject_sigma;
  Float_t pion_1_tof_kreject_pmin;

  Float_t pion_1_tof_preject_sigma;
  Float_t pion_1_tof_preject_pmin;

  Float_t pion_1_tpc_ereject_sigma;

  Float_t pion_1_max_impact_xy,
          pion_1_max_impact_z,
          pion_1_min_tpc_chi_ndof,
          pion_1_max_tpc_chi_ndof,
          pion_1_max_its_chi_ndof;
  UInt_t pion_1_min_tpc_ncls;
  UInt_t pion_1_min_its_ncls;

  Bool_t pion_1_remove_kinks,
         pion_1_rm_neg_lbl,
         pion_1_use_tpctof;

  ULong_t pion_1_ideal_pid;

  // PION - 2
  Float_t pion_2_PtMin,
          pion_2_PtMax;

  Float_t pion_2_EtaMin,
          pion_2_EtaMax;

  Float_t pion_2_DCAMin,
          pion_2_DCAMax;

  Float_t pion_2_NSigmaMin,
          pion_2_NSigmaMax;

  Float_t pion_2_max_impact_xy,
          pion_2_max_impact_z,
          pion_2_max_tpc_chi_ndof,
          pion_2_max_its_chi_ndof;

  UInt_t pion_2_min_tpc_ncls;
  Bool_t pion_2_remove_kinks,
         pion_2_set_label;


  // PAIR
  Bool_t pair_use_avgsep;
  Bool_t pair_TPCOnly;
  // Float_t pair_TPCExitSepMin;
  Float_t pair_min_avgsep;
  // Float_t pair_MinAvgSeparationNeg;

  Float_t pair_delta_eta_min,
          pair_delta_phi_min,
          pair_phi_star_radius;

  Float_t pair_ee_min;

  Float_t pair_max_share_quality,
          pair_max_share_fraction;
  Bool_t pair_remove_same_label;
  Int_t pair_algorithm;

  /// Default Values
  CutParams();


  ClassDef(CutParams, 0);
};

template <typename T>
AliFemtoConfigObject
AliFemtoAnalysisPionPion::GetConfigurationOf(const T &cut)
{
  auto *cls = TClass::GetClass(typeid(cut));
  if (cls) {
    return AliFemtoConfigObject("");
  }
  AliFemtoConfigObject::MapValue_t result;
  result["_class"] = "AliFemtoSomething";

  return AliFemtoConfigObject(result);
}

#endif
