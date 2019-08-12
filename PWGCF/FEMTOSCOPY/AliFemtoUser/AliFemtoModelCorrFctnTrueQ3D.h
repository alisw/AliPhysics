///
/// \file AliFemtoModelCorrFctnTrueQ3D.h
///

#pragma once

#ifndef ALIFEMTOMODELCORRFCTN_TRUEQ3D_H
#define ALIFEMTOMODELCORRFCTN_TRUEQ3D_H


#include "AliFemtoCorrFctn.h"

class AliFemtoPair;
class AliFemtoModelManager;
class TH3F;
class TRandom;

/// \class AliFemtoModelCorrFctnTrueQ3D
/// \brief Correlation function storing true momentum in LCMS frame.
///
///
/// Build with the parameters class
///
/// ```cpp
/// AliFemtoModelCorrFctnTrueQ3D *cf = AliFemtoModelCorrFctnTrueQ3D::Build()
///                                        .NamePrefix("MRC")
///                                        .QRange(-.2, .2)
///                                        .BinCount(87)
///                                        .Manager(mc_manager);
///
/// ```
///
/// \author Andrew Kubera, The Ohio State University <andrew.michael.kubera@cern.ch>
///
class AliFemtoModelCorrFctnTrueQ3D : public AliFemtoCorrFctn {
public:

  struct Parameters {
    UInt_t bin_count;
    Double_t qmin;
    Double_t qmax;
    TString prefix;
    bool enable_extra_hists;
    bool enable_extra_denoms;
    AliFemtoModelManager *mc_manager;

    /// Build Parameters object with default values
    Parameters()
      : bin_count(56)
      , qmin(-0.14)
      , qmax(0.14)
      , prefix("")
      , enable_extra_hists(false)
      , enable_extra_denoms(false)
      , mc_manager(NULL)
      {}

    Parameters(const Parameters &orig)
      : bin_count(orig.bin_count)
      , qmin(orig.qmin)
      , qmax(orig.qmax)
      , prefix(orig.prefix)
      , enable_extra_hists(orig.enable_extra_hists)
      , enable_extra_denoms(orig.enable_extra_denoms)
      , mc_manager(orig.mc_manager)
      { }

    Parameters& operator=(const Parameters &rhs)
      {
        bin_count = rhs.bin_count;
        qmin = rhs.qmin;
        qmax = rhs.qmax;
        prefix = rhs.prefix;
        enable_extra_hists = rhs.enable_extra_hists;
        enable_extra_denoms = rhs.enable_extra_denoms;
        mc_manager = rhs.mc_manager;
        return *this;
      }

    static Parameters Default()
    {
      return Parameters().NamePrefix("CF_TrueQ3D");
    }

    Parameters NamePrefix(const TString& prefix) const {
      Parameters p(*this);
      p.prefix = prefix;
      return p;
    }
    Parameters Title(const TString& title) const {
      Parameters p(*this);
      p.prefix = title;
      return p;
    }
    Parameters EnableExtraHists(bool t) const {
      Parameters p(*this);
      p.enable_extra_hists = t;
      return p;
    }
    Parameters EnableWeightedDenominators(bool t) const {
      Parameters p(*this);
      p.enable_extra_denoms = t;
      return p;
    }
    Parameters Manager(AliFemtoModelManager *manager) const {
      Parameters p(*this);
      p.mc_manager = manager;
      return p;
    }
    Parameters QRange(Double_t max) const {
      Parameters p(*this);
      p.qmax = std::abs(max);
      p.qmin = -p.qmax;
      return p;
    }
    Parameters QRange(Double_t min, Double_t max) const {
      Parameters p(*this);
      p.qmax = max;
      p.qmin = min;
      return p;
    }
    Parameters BinCount(UInt_t nbins) const {
      Parameters p(*this);
      p.bin_count = nbins;
      return p;
    }
    Parameters NBin(UInt_t nbins) const {
      Parameters p(*this);
      p.bin_count = nbins;
      return p;
    }
    Parameters AxisInfo(UInt_t nbins, Double_t q_max) const {
      Parameters p(*this);
      p.bin_count = nbins;
      p.qmax = std::abs(q_max);
      p.qmin = -p.qmax;
      return p;
    }
    Parameters AxisInfo(UInt_t nbins, Double_t q_min, Double_t q_max) const {
      Parameters p(*this);
      p.bin_count = nbins;
      p.qmin = q_min;
      p.qmax = q_max;
      return p;
    }

    operator AliFemtoModelCorrFctnTrueQ3D*() const {
      return into_ptr();
    }

    AliFemtoModelCorrFctnTrueQ3D* into_ptr() const {
      return new AliFemtoModelCorrFctnTrueQ3D(*this);
    }
  };

  /// Deafult parameters
  ///
  /// - Prefix: "CF_TrueQ3D"
  /// - Binning Paramters: (56, -0.14, 0.14)
  ///
  AliFemtoModelCorrFctnTrueQ3D();

  /// Construct with name-prefix
  ///
  /// Use default binning parameters
  AliFemtoModelCorrFctnTrueQ3D(const char *prefix);

  /// Symmetric constructor
  ///
  /// Construct with nbins from -qmax to qmax in both directions
  AliFemtoModelCorrFctnTrueQ3D(const char *prefix, UInt_t nbins, Double_t qmax);

  /// Custom constructor
  ///
  AliFemtoModelCorrFctnTrueQ3D(const TString &prefix, UInt_t nbins, Double_t aQinvLo, Double_t aQinvHi);

  /// Big Constructor
  ///
  /// enable_extra_hists -- adds four complimentary weighted/unweighted histograms
  ///
  AliFemtoModelCorrFctnTrueQ3D(const TString &prefix,
                               UInt_t nbins,
                               Double_t aQinvLo,
                               Double_t aQinvHi,
                               Bool_t enable_extra_hists,
                               Bool_t enable_extra_denoms);

  /// Build with variable bins
  ///
  AliFemtoModelCorrFctnTrueQ3D(const TString &prefix,
                               const std::vector<double> &obins,
                               const std::vector<double> &sbins,
                               const std::vector<double> &lbins,
                               AliFemtoModelManager *mc_manager=nullptr);

  /// Unnamed, q-symmetric constructor
  ///
  AliFemtoModelCorrFctnTrueQ3D(UInt_t nbins, Double_t qmax);

  /// Unnamed constructor
  ///
  AliFemtoModelCorrFctnTrueQ3D(UInt_t nbins, Double_t qmin, Double_t qmax);

  /// Construct from parameter object
  ///
  AliFemtoModelCorrFctnTrueQ3D(const Parameters &params);

  /// copy constructor - Structure preserved - no histogram contents copied
  AliFemtoModelCorrFctnTrueQ3D(const AliFemtoModelCorrFctnTrueQ3D& aCorrFctn);

  /// Assignment Operator - unused
  //   AliFemtoModelCorrFctnTrueQ3D& operator=(const AliFemtoModelCorrFctnTrueQ3D& aCorrFctn);

  /// Destructor - histograms destroyed, ModelManager is NOT
  virtual ~AliFemtoModelCorrFctnTrueQ3D();

  /// assignement operator - copy any histograms
  AliFemtoModelCorrFctnTrueQ3D& operator=(const AliFemtoModelCorrFctnTrueQ3D&);

  /// Set the MC model manager
  virtual void SetManager(AliFemtoModelManager *manager)
    { fManager = manager; }

  virtual AliFemtoString Report();

  /// Add pair of particles from same event
  virtual void AddRealPair(AliFemtoPair* aPair);

  /// Add pair of particles from differnt events
  virtual void AddMixedPair(AliFemtoPair* aPair);

  /// no-op
  virtual void Finish();

  virtual TList* GetOutputList();
  virtual void AddOutputObjectsTo(TCollection &);
  virtual TList* AppendOutputList(TList &);

  /// Return copied corr fctn
  virtual AliFemtoCorrFctn* Clone() const
    { return new AliFemtoModelCorrFctnTrueQ3D(*this); }

  /// Get a builder-pattern constructor object
  static Parameters Build()
    { return Parameters(); }

protected:
  AliFemtoModelManager *fManager; //!<! Link back to the manager to retrieve weights

  /// Numerator made with pairs using monte-carlo generated momentum
  TH3F *fNumeratorGenerated;

  /// Numerator made with pairs from the same event
  TH3F *fNumeratorReconstructed;

  /// Numerator made with pairs from the same event - no femto weighting
  TH3F *fNumeratorGenUnweighted;
  TH3F *fNumeratorRecUnweighted;

  /// Denominator with the monte-carlo generated momentum
  TH3F *fDenominatorGenerated;

  /// Denominator with reconstructed data
  TH3F *fDenominatorReconstructed;

  TH3F *fDenominatorGenWeighted;
  TH3F *fDenominatorRecWeighted;

};

#endif
