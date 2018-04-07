///
/// \file AliFemtoUser/AliFemtoAnalysisPionPionObjectConstructor.h
///
/// This header file depends on C++11, and thus can be depended on at compile-time,
/// but not runtime (i.e. .cxx files may #include it, but not .h files)
///
/// This is a strange template header file that checks for include-guards already defined by other files
/// to avoid importing EVERY CLASS DEFINITION for every cxx file that includes it. This means that this
/// must be the last header included.
///


#ifndef ALIFEMTOANALYSISPIONPIONOBJECTCONSTRUCTORS_H
#define ALIFEMTOANALYSISPIONPIONOBJECTCONSTRUCTORS_H

#if __cplusplus < 201103L
#error This files requires C++11
#endif

#include <TInterpreter.h>

#include "AliFemtoConfigObject.h"

/// Generic Configuration struct used to convert AliFemtoConfigObjects
/// into actual AliFemto objects
///
/// The template is the "target" type, and SHOULD implement the method
/// `void Configure(T&) const` used by other classes to actually setup
/// an object.
///
/// These methods should have a constructor taking an AliFemtoConfigObject
/// by non-const reference
///
template <typename T>
struct Configuration;


/// Specific case of a Configuration class; one that wraps an abstract-
/// base-class used by the others for specifying possible pointer
/// transformations.
///
/// Use their virtual functions to create objects they
///
/// For example:
///
/// If `AliFemtoBasicEventCut` and `AliFemtoEventCutCentrality` are
/// subclasses of the same `AliFemtoEventCut` base class, then their
/// associated Configuration classes (`Configuration<AliFemtoBasicEventCut>`
/// & `Configuration<AliFemtoEventCutCentrality>`) should both inherit
/// from `AbstractConfiguration<AliFemtoEventCut>` and implement the
/// `operator AliFemtoEventCut*() const` method.
///
/// TODO: Rethink using the name 'Abstract'
///
template <typename T>
struct AbstractConfiguration {
  virtual operator T*() const = 0;
};

#endif



//---------------------------------------------------------
//
//   Event Readers
//
//------------------

#if defined(ALIFEMTOEVENTREADERAOD_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOEVENTREADERAOD_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOEVENTREADERAOD_H
///
template <>
struct Configuration<AliFemtoEventReaderAOD> : AbstractConfiguration<AliFemtoEventReader> {
  int filter_bit { 7 };
  AliFemtoEventReaderAOD::EstEventMult multiplicity { AliFemtoEventReaderAOD::kCentrality };

  bool read_v0 { false };
  bool read_mc { false };
  bool evp_zero { false };
  bool global_tack_DCA { true };
  bool centrality_flattening { false };

  Configuration(AliFemtoConfigObject &cfg)
  {
    std::string mult_str;
    int m_int;
    // try loading multiplicity as a string
    if (cfg.pop_and_load("multiplicity", mult_str)) {
      std::cout << "Setting multiplicity " << mult_str << "\n";
      gInterpreter->ProcessLine(Form("std::cout << '~' << %s << '\\n';", mult_str.c_str()));
    }

    else if (cfg.pop_and_load("multiplicity", m_int)) {
      std::cout << "Setting multiplicity " << m_int << "\n";
      multiplicity = static_cast<AliFemtoEventReaderAOD::EstEventMult>(m_int);

      // NOTE : These should be the bounds of the enum
      if (m_int < AliFemtoEventReaderAOD::kCentrality || AliFemtoEventReaderAOD::kCentralityNPA < m_int) {
        std::cerr << "Invalid multiplicity " << m_int << "\n";
      }
    }

    cfg.pop_all()
      ("filterbit", filter_bit)
      ("evp_zero", evp_zero)
      ("global_tack_DCA", global_tack_DCA)
      ("centrality_flattening", centrality_flattening)
      ("readV0", read_v0)
      ("read_mc", read_mc)
      .WarnOfRemainingItems();
  }

  void Configure(AliFemtoEventReaderAOD &obj) const
  {
    obj.SetFilterBit(filter_bit);
    obj.SetEPVZERO(evp_zero);
    obj.SetUseMultiplicity(multiplicity);
    obj.SetCentralityFlattening(centrality_flattening);
    obj.SetReadV0(read_v0);
    // obj.SetPrimaryVertexCorrectionTPCPoints(kTRUE);
    obj.SetDCAglobalTrack(global_tack_DCA);
    obj.SetReadMC(read_mc);
  }

  operator AliFemtoEventReader*() const
  {
    AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAOD();
    Configure(*rdr);
    return rdr;
  }

};
#endif


#if defined(ALIFEMTOEVENTREADERAODCHAIN_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOEVENTREADERAODCHAIN_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOEVENTREADERAODCHAIN_H
template <>
struct Configuration<AliFemtoEventReaderAODChain> : Configuration<AliFemtoEventReaderAOD> {

  Configuration(AliFemtoConfigObject &cfg):
    Configuration<AliFemtoEventReaderAOD>(cfg)
  {}

  void Configure(AliFemtoEventReaderAODChain &obj) const
  {
    Configuration<AliFemtoEventReaderAOD>::Configure(obj);
  }

  virtual operator AliFemtoEventReader*() const
  {
    auto *rdr = new AliFemtoEventReaderAODChain();
    Configure(*rdr);
    return rdr;
  }
};
#endif



#if defined(ALIFEMTOEVENTREADER_MULTSELECTION_H_) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOEVENTREADER_MULTSELECTION_H_)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOEVENTREADER_MULTSELECTION_H_
template <>
struct Configuration<AliFemtoEventReaderAODMultSelection> : Configuration<AliFemtoEventReaderAODChain> {
  AliFemtoEventReaderAOD::EstEventMult multiplicity { AliFemtoEventReaderAOD::kCentrality };

  bool read_v0 { false };
  bool read_mc { false };
  bool evp_zero { false };
  bool global_tack_DCA { true };
  bool centrality_flattening { false };

  void Configure(AliFemtoEventReaderAODMultSelection &obj) const
  {
    Configuration<AliFemtoEventReaderAODChain>::Configure(obj);
    obj.SetFilterBit(filter_bit);
    obj.SetEPVZERO(evp_zero);
    obj.SetUseMultiplicity(multiplicity);
    obj.SetCentralityFlattening(centrality_flattening);
    obj.SetReadV0(read_v0);
    // obj.SetPrimaryVertexCorrectionTPCPoints(kTRUE);
    obj.SetDCAglobalTrack(global_tack_DCA);
    obj.SetReadMC(read_mc);

  }

  Configuration(AliFemtoConfigObject &cfg):
    Configuration<AliFemtoEventReaderAODChain>(cfg)
  {
    std::string mult_str;
    int m_int;
    // try loading multiplicity as a string
    if (cfg.pop_and_load("multiplicity", mult_str)) {
      std::cout << "Setting multiplicity " << mult_str << "\n";
      gInterpreter->ProcessLine(Form("std::cout << '~' << %s << '\\n';", mult_str.c_str()));
    }

    else if (cfg.pop_and_load("multiplicity", m_int)) {
      std::cout << "Setting multiplicity " << m_int << "\n";
      multiplicity = static_cast<AliFemtoEventReaderAOD::EstEventMult>(m_int);

      // NOTE : These should be the bounds of the enum
      if (m_int < AliFemtoEventReaderAOD::kCentrality || AliFemtoEventReaderAOD::kCentralityNPA < m_int) {
        std::cerr << "Invalid multiplicity " << m_int << "\n";
      }
    }

    cfg.pop_all()
      ("filterbit", filter_bit)
      ("evp_zero", evp_zero)
      ("global_tack_DCA", global_tack_DCA)
      ("centrality_flattening", centrality_flattening)
      ("readV0", read_v0)
      ("read_mc", read_mc)
      .WarnOfRemainingItems();
  }

  virtual operator AliFemtoEventReader*() const
  {
    AliFemtoEventReaderAODMultSelection *rdr = new AliFemtoEventReaderAODMultSelection();
    Configure(*rdr);
    return rdr;
  }
};
#endif




//------------------------
//
//   Event Cuts
//
//------------------------


#if defined(ALIFEMTOBASICEVENTCUT_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOBASICEVENTCUT_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOBASICEVENTCUT_H

template <>
struct Configuration<AliFemtoBasicEventCut> : AbstractConfiguration<AliFemtoEventCut> {
  using RangeF_t = AliFemtoConfigObject::RangeValue_t;

  std::pair<int, int> multiplicity = {0, 1000000};
  RangeF_t centrality = {0, 90},
           vertex_z = {-10.0f, 10.0f},
           EP_VZero = {-1000.0, 1000.0};

  Int_t trigger_selection = 0;
    Bool_t accept_bad_vertex = false;
    Bool_t accept_only_physics = true;

    UInt_t min_coll_size = 10;

    /// default constructor required to use default initialized members
    Configuration(){};

    /// Templated member for constructing AliFemtoEventCut objects from
    /// these parameters.
    operator AliFemtoEventCut*() const
    {
      auto ptr = new AliFemtoBasicEventCut();
      Configure(*ptr);
      return ptr;
    }

    void Configure(AliFemtoBasicEventCut &cut) const {
      cut.SetTriggerSelection(trigger_selection);
    }

    /// Build from config object
    Configuration(AliFemtoConfigObject obj) {

      obj.pop_all()
        ("multiplicity", multiplicity)
        ("centrality", centrality)
        ("vertex_z", vertex_z)
        ("ep_v0", EP_VZero)
        ("trigger", trigger_selection)
        ("accept_bad_vertex", accept_bad_vertex)
        ("accept_only_physics", accept_only_physics)
        ("min_collection_size", min_coll_size)
        .WarnOfRemainingItems();
    }

    /// Construct a config object with this object's properties
    operator AliFemtoConfigObject() const {
        return AliFemtoConfigObject::BuildMap()
          ("multiplicity", multiplicity)
          ("centrality", centrality)
          ("vertex_z", vertex_z)
          ("ep_v0", EP_VZero)
          ("trigger", trigger_selection)
          ("accept_bad_vertex", accept_bad_vertex)
          ("accept_only_physics", accept_only_physics)
          ("min_collection_size", min_coll_size);
    }
};
#endif

/*
#if defined(AliFemtoEventCut_hh) && !defined(ALIFEMTOCONFSTRUCTOR_AliFemtoEventCut_hh)
#define ALIFEMTOCONFSTRUCTOR_AliFemtoEventCut_hh
template <>
struct AbstractConfiguration<AliFemtoEventCut> {
  operator AliFemtoe
};
#endif
*/

#if defined(ALIFEMTOEVENTCUTCENTRALITY_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOEVENTCUTCENTRALITY_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOEVENTCUTCENTRALITY_H

template<>
struct Configuration<AliFemtoEventCutCentrality> : AbstractConfiguration<AliFemtoEventCut> {
  using Super = AbstractConfiguration<AliFemtoEventCut>;

  using RangeF_t = AliFemtoConfigObject::RangeValue_t;

  RangeF_t centrality, z, epvzero;
  int trigger_selection;

  Configuration(AliFemtoConfigObject &cfg)
  {
    cfg.pop_all()
      ("centrality", centrality)
      ("ep_v0", epvzero)
      ;
  }

  void Configure(AliFemtoEventCutCentrality &cut) const
  {
    //Super::Configure(cut);
    cut.SetCentralityRange(centrality.first, centrality.second);
    cut.SetZPosRange(z.first, z.second);
    cut.SetEPVZERO(epvzero.first, epvzero.second);
    cut.SetTriggerSelection(trigger_selection);
  }

  virtual operator AliFemtoEventCut*() const {
    AliFemtoEventCutCentrality *cut = new AliFemtoEventCutCentrality();
    Configure(*cut);
    return cut;
  }
};
#endif



//------------------------
//
//   Particle Cuts
//
//------------------------


#if defined(ALIFEMTOPARTICLECUT_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOPARTICLECUT_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOPARTICLECUT_H
template<>
struct AbstractConfiguration<AliFemtoParticleCut>
{
  Double_t mass;
  void Configure(AliFemtoParticleCut &cut) const {
    cut.SetMass(mass);
  }

  virtual operator AliFemtoParticleCut*() const = 0;
};
#endif


#if defined(ALIFEMTOESDTRACKCUT_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOESDTRACKCUT_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOESDTRACKCUT_H
template<>
struct Configuration<AliFemtoESDTrackCut> : AbstractConfiguration<AliFemtoParticleCut>
{
  using Super = AbstractConfiguration<AliFemtoParticleCut>;
  using RangeF_t = std::pair<float, float>;
  RangeF_t pt = {0.2, 2.0},
           rapidity = {-2.0, 2.0},
            eta = {-0.8, 0.8},
            DCA = {0.5, 4.0},
            nSigma = {-3.0, 3.0};

  Int_t charge = 0,
        label = 0,
        min_tpc_ncls = 80;

  Float_t max_impact_xy = 2.4,
          min_impact_xy = 0,
          max_impact_z = 3.0,
          max_tpc_chi_ndof = 0.032,
          max_its_chi_ndof = 0.032;
  Float_t sigma = 3.;
  Bool_t electron_rejection = true,
         remove_kinks = true;

  // Pion
  Int_t most_probable { 2 };

  void Configure(AliFemtoESDTrackCut &cut) const {
    Super::Configure(cut);
    cut.SetCharge(charge);
    cut.SetPt(pt.first, pt.second);
    cut.SetEta(eta.first, eta.second);
    cut.SetRapidity(rapidity.first, rapidity.second);
    cut.SetMostProbablePion();

    /// Settings for TPC-Inner Runmode
    cut.SetStatus(AliESDtrack::kTPCin);
    cut.SetminTPCncls(min_tpc_ncls);
    cut.SetLabel(label);
    cut.SetMaxTPCChiNdof(max_tpc_chi_ndof);
    cut.SetMaxITSChiNdof(max_its_chi_ndof);

    cut.SetMaxImpactZ(max_impact_z);
    cut.SetMaxImpactXY(max_impact_xy);
    cut.SetMinImpactXY(min_impact_xy);
    cut.SetElectronRejection(electron_rejection);
    cut.SetRemoveKinks(remove_kinks);
  }

  virtual operator AliFemtoParticleCut*() const {
    AliFemtoESDTrackCut *cut = new AliFemtoESDTrackCut();
    Configure(*cut);
    return cut;
  }

};
#endif


//------------------------
//
//   Pair Cuts
//
//------------------------


#if defined(ALIFEMTOSHAREQUALITYPAIRCUT_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOSHAREQUALITYPAIRCUT_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOSHAREQUALITYPAIRCUT_H
template<>
struct Configuration<AliFemtoShareQualityPairCut> : AbstractConfiguration<AliFemtoPairCut> {
  Bool_t remove_same_label { false };
  Double_t max_share_quality { 1.0 },
           max_share_fraction { 0.05 };


  Configuration(AliFemtoConfigObject &cfg)
  {
    cfg.pop_all()
      ("max_share_quality", max_share_quality)
      ("max_share_fraction", max_share_fraction)
      ("remove_same_label", remove_same_label)
      ;
  }

  virtual operator AliFemtoPairCut*() const {
    AliFemtoShareQualityPairCut *result = new AliFemtoShareQualityPairCut();
    Configure(*result);
    return result;
  }

  void Configure(AliFemtoShareQualityPairCut &cut) const {
    cut.SetShareQualityMax(max_share_quality);
    cut.SetShareFractionMax(max_share_fraction);
    cut.SetRemoveSameLabel(remove_same_label);
  }

  static void ReadConfigurationInto(AliFemtoConfigObject &dest, const AliFemtoShareQualityPairCut &cut)
  {
    dest.Update(AliFemtoConfigObject::BuildMap()
                ("max_share_quality", cut.GetShareQualityMax())
                ("max_share_fraction", cut.GetShareFractionMax())
                ("remove_same_label", cut.GetRemoveSameLabel()));
  }

  static AliFemtoConfigObject GetConfigurationOf(const AliFemtoShareQualityPairCut &cut)
  {
    auto result = AliFemtoConfigObject::BuildMap();
    ReadConfigurationInto(result, cut);
    return result;
  }
};
#endif


#if defined(ALIFEMTOPAIRCUTDETADPHI_H_) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOPAIRCUTDETADPHI_H_)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOPAIRCUTDETADPHI_H_
template<>
struct Configuration<AliFemtoPairCutDetaDphi> : Configuration<AliFemtoShareQualityPairCut> {
  using Super = Configuration<AliFemtoShareQualityPairCut>;

  Double_t delta_eta_min { 0.0 },
           delta_phi_min { 0.0 },
           radius { 1.2 };
  Bool_t use_quad { true };

  Configuration(AliFemtoConfigObject &cfg):
    Super(cfg)
  {
    cfg.pop_all()
      ("use_quad", use_quad)
      ("radius", radius)
      ("delta_phi_min", delta_phi_min)
      ("delta_eta_min", delta_eta_min);

  }

  virtual operator AliFemtoPairCut*() const {
    AliFemtoPairCutDetaDphi *result = new AliFemtoPairCutDetaDphi();
    Configure(*result);
    return result;
  }

  void Configure(AliFemtoPairCutDetaDphi &cut) const
  {
    Super::Configure(cut);
    cut.SetCutTechnique(use_quad ? AliFemtoPairCutDetaDphi::Quad : AliFemtoPairCutDetaDphi::Simple);
    cut.SetR(radius);
    cut.SetMinEta(delta_eta_min);
    cut.SetMinPhi(delta_phi_min);
  }

  static AliFemtoConfigObject GetConfigOf(const AliFemtoPairCutDetaDphi &cut)
  {
    return AliFemtoConfigObject::BuildMap()
      ("class", "AliFemtoPairCutDetaDphi")
      ("share_fraction_max", cut.GetShareFractionMax())
      ("share_quality_max", cut.GetShareQualityMax())
      ("remove_same_label", cut.GetRemoveSameLabel())
      ("use_quad", cut.GetUsesQuadratureTechnique())
      ("radius", cut.GetRadius())
      ("delta_phi_min", cut.GetMinDeltaPhi())
      ("delta_eta_min", cut.GetMinDeltaEta());
  }
};
#endif

#if defined(ALIFEMTOPAIRCUTANTIGAMMA_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOPAIRCUTANTIGAMMA_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOPAIRCUTANTIGAMMA_H
template<>
struct Configuration<AliFemtoPairCutAntiGamma> : Configuration<AliFemtoShareQualityPairCut> {
  using Super = Configuration<AliFemtoShareQualityPairCut>;

  Double_t minv_ee_max { 0.0 },
           dtheta_max { 0.0 },
           dtpc_min { 0.0 },
           min_avg_sep { 0.0 };

  AliFemtoPairCut::AliFemtoDataType data_type { AliFemtoPairCut::kAOD };

  Configuration(AliFemtoConfigObject &cfg):
    Super(cfg)
  {
    std::cout << "Configuration<AliFemtoPairCutAntiGamma>!!\n";
    cfg.pop_and_load("min_avg_sep", min_avg_sep);
  }

  virtual operator AliFemtoPairCut*() const
  {
    AliFemtoPairCutAntiGamma *result = new AliFemtoPairCutAntiGamma();
    Configure(*result);
    return result;
  }

  void Configure(AliFemtoShareQualityPairCut &cut) const
  {
    Super::Configure(cut);
  }
};

#endif


//---------------------------------------------------------
//
//   Correlation Functions
//
//------------------------


#if defined(ALIFEMTOCORRFCTN_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOCORRFCTN_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOCORRFCTN_H
template<>
struct AbstractConfiguration<AliFemtoCorrFctn> {

  AliFemtoConfigObject pair_cut_cfg;

  AbstractConfiguration(AliFemtoConfigObject &obj)
  {
    if (obj.pop_and_load("pair_cut", pair_cut_cfg)) {
      if (!pair_cut_cfg.is_map()) {
        std::cerr << "'pair_cut' object should be a map describing a pair cut;"
                     " instead got " << pair_cut_cfg.name_of_type() << "\n";
      }
    }
  }

  void Configure(AliFemtoCorrFctn &cf) const
  {
    if (!pair_cut_cfg.is_empty()) {
      cf.SetPairSelectionCut(AliFemtoAnalysisPionPion::ConstructPairCut(pair_cut_cfg));
    }
  }
  virtual operator AliFemtoCorrFctn*() const = 0;
};
#endif



#if defined(ALIFEMTOAVGSEPCORRFCTN_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOAVGSEPCORRFCTN_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOAVGSEPCORRFCTN_H
template<>
struct Configuration<AliFemtoAvgSepCorrFctn> : AbstractConfiguration<AliFemtoCorrFctn> {
  using Super = AbstractConfiguration<AliFemtoCorrFctn>;

  std::string name;
  Int_t nbins { 100 };
  Double_t low { 0.0 },
           high { 40.0 };

  Configuration(AliFemtoConfigObject &obj)
  : Super(obj)
  {
    obj.pop_all()
      ("name", name)
      ("bin_count", nbins)
      ("qmin", low)
      ("qmax", high)
      .WarnOfRemainingItems();
  }

  virtual operator AliFemtoCorrFctn*() const {
    AliFemtoAvgSepCorrFctn *ptr = new AliFemtoAvgSepCorrFctn(name.c_str(), nbins, low, high);
    Super::Configure(*ptr);
    return ptr;
  }
};
#endif


#if defined(ALIFEMTOMODELCORRFCTN_TRUEQ3D_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOMODELCORRFCTN_TRUEQ3D_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOMODELCORRFCTN_TRUEQ3D_H
template<>
struct Configuration<AliFemtoModelCorrFctnTrueQ3D> : AbstractConfiguration<AliFemtoCorrFctn> {
  using Super = AbstractConfiguration<AliFemtoCorrFctn>;

  AliFemtoModelCorrFctnTrueQ3D::Parameters params;

  Configuration(AliFemtoConfigObject &obj)
  : Super(obj)
  {
    obj.pop_all()
      ("title", params.title)
      ("bin_count", params.bin_count)
      ("qmin", params.qmin)
      ("qmax", params.qmax)
      .WarnOfRemainingItems();
  }

  virtual operator AliFemtoCorrFctn*() const {
    AliFemtoModelCorrFctnTrueQ3D *ptr = new AliFemtoModelCorrFctnTrueQ3D(params);
    Super::Configure(*ptr);
    return ptr;
  }
};
#endif


#if defined(ALIFEMTOMODELCORRFCTN_TRUEQ3D_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOMODELCORRFCTN_TRUEQ3D_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOMODELCORRFCTN_TRUEQ3D_H
template<>
struct Configuration<AliFemtoModelCorrFctnDEtaDPhiStar> : AbstractConfiguration<AliFemtoCorrFctn> {
  using Super = AbstractConfiguration<AliFemtoCorrFctn>;

  AliFemtoModelCorrFctnTrueQ3D::Parameters params;

  Configuration(AliFemtoConfigObject &obj)
  : Super(obj)
  {
    obj.pop_all()
      ("title", params.title)
      ("bin_count", params.bin_count)
      ("qmin", params.qmin)
      ("qmax", params.qmax)
      .WarnOfRemainingItems();
  }

  virtual operator AliFemtoCorrFctn*() const {
    AliFemtoModelCorrFctnDEtaDPhiStar *ptr = new AliFemtoModelCorrFctnDEtaDPhiStar(params);
    Super::Configure(*ptr);
    return ptr;
  }
};
#endif


#if defined(ALIFEMTOCORRFCTN3DLCMS_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOCORRFCTN3DLCMS_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOCORRFCTN3DLCMS_H
template<>
struct Configuration<AliFemtoCorrFctn3DLCMSSym> : AbstractConfiguration<AliFemtoCorrFctn> {
  using Super = AbstractConfiguration<AliFemtoCorrFctn>;

  Bool_t use_LCMS { false };

  TString title;
  Int_t nbins { 100 };
  Float_t qmax { 1.0 };

  Configuration(AliFemtoConfigObject &obj)
  : Super(obj)
  {
    obj.pop_all()
      ("title", title)
      ("bin_count", nbins)
      ("qmax", qmax)
      .WarnOfRemainingItems();
  }

  virtual operator AliFemtoCorrFctn*() const {
    AliFemtoCorrFctn3DLCMSSym *ptr = new AliFemtoCorrFctn3DLCMSSym(title, nbins, qmax);
    Super::Configure(*ptr);
    ptr->SetUseLCMS(use_LCMS);
    return ptr;
  }
};
#endif


#if defined(ALIFEMTOCORRFCTNDIRECTYLM_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOCORRFCTNDIRECTYLM_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOCORRFCTNDIRECTYLM_H
template<>
struct Configuration<AliFemtoCorrFctnDirectYlm> : AbstractConfiguration<AliFemtoCorrFctn> {
  using Super = AbstractConfiguration<AliFemtoCorrFctn>;

  TString name { "cf" };

  int maxL { 4 },
      ibin { 30 };
  Double_t vmin { 0.0 },
           vmax { 0.3 };

  Bool_t use_LCMS { false };

  Configuration(AliFemtoConfigObject &obj)
  : Super(obj)
  {
    obj.pop_all()
      ("title", name)
      ("maxL", maxL)
      ("bin_count", ibin)
      ("vmin", vmin)
      ("vmax", vmax)
      ("use_LCMS", use_LCMS)
      .WarnOfRemainingItems();
  }

  virtual operator AliFemtoCorrFctn*() const {
    AliFemtoCorrFctnDirectYlm *ptr = new AliFemtoCorrFctnDirectYlm(name, maxL, ibin, vmin, vmax, use_LCMS);
    Super::Configure(*ptr);
    return ptr;
  }
};
#endif


#if defined(ALIFEMTOCORRFCTNDPHISTARDETA_H) && !defined(ALIFEMTOCONSTRUCTOR_ALIFEMTOCORRFCTNDPHISTARDETA_H)
#define ALIFEMTOCONSTRUCTOR_ALIFEMTOCORRFCTNDPHISTARDETA_H
template<>
struct Configuration<AliFemtoCorrFctnDPhiStarDEta> : AbstractConfiguration<AliFemtoCorrFctn> {
  using Super = AbstractConfiguration<AliFemtoCorrFctn>;

  TString name { "cf" };

  AliFemtoConfigObject::RangeValue_t eta_range {-0.1, 0.1},
                                     phi_range {-0.1, 0.1};

  Double_t radius { 0.8 };
  Int_t eta_bins { 50 },
        phi_bins { 50 };


  Configuration(AliFemtoConfigObject &obj)
  : Super(obj)
  {
    obj.pop_all()
      ("title", name)
      ("radius", radius)
      ("eta_bins", eta_bins)
      ("phi_bins", phi_bins)
      ("eta_range", eta_range)
      ("phi_range", phi_range)
      .WarnOfRemainingItems();
  }

  virtual operator AliFemtoCorrFctn*() const {
    AliFemtoCorrFctnDPhiStarDEta *ptr = new AliFemtoCorrFctnDPhiStarDEta(name,
                                                                         radius,
                                                                         eta_bins,
                                                                         eta_range.first,
                                                                         eta_range.second,
                                                                         phi_bins,
                                                                         phi_range.first,
                                                                         phi_range.second);
    Super::Configure(*ptr);
    return ptr;
  }

  static AliFemtoConfigObject GetConfigOf(const AliFemtoCorrFctnDPhiStarDEta &cut)
  {
    return AliFemtoConfigObject::BuildMap()
      ("class", "AliFemtoCorrFctnDPhiStarDEta")
      ("share_fraction_max", cut.GetShareFractionMax())
      ("share_quality_max", cut.GetShareQualityMax())
      ("remove_same_label", cut.GetRemoveSameLabel());
  }
};
#endif
