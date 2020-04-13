///
/// \file AliFemtoModelCorrFctnTrueQ3DByParent.h
///

#pragma once

#ifndef ALIFEMTOMODELCORRFCTNQ3DBYPARENT_H
#define ALIFEMTOMODELCORRFCTNQ3DBYPARENT_H


#include "AliFemtoCorrFctn.h"

#include <THnSparse.h>

class AliFemtoPair;
class AliFemtoModelManager;
class TRandom;


/// \class AliFemtoModelCorrFctnTrueQ3DByParent
/// \brief Correlation function storing true momentum in LCMS frame.
///
/// \author Andrew Kubera, The Ohio State University <andrew.michael.kubera@cern.ch>
///
class AliFemtoModelCorrFctnTrueQ3DByParent : public AliFemtoCorrFctn {
public:

  struct Parameters {
    Int_t bin_count;
    Double_t qmin;
    Double_t qmax;
    TString prefix;
    AliFemtoModelManager *mc_manager;

    Parameters()
      : bin_count(57)
      , qmin(-0.1425)
      , qmax(0.1425)
      , prefix("")
      , mc_manager(NULL)
      { }

    Parameters(const Parameters &orig)
      : bin_count(orig.bin_count)
      , qmin(orig.qmin)
      , qmax(orig.qmax)
      , prefix(orig.prefix)
      , mc_manager(orig.mc_manager)
      { }

    Parameters& operator=(const Parameters &rhs)
      {
        bin_count = rhs.bin_count;
        qmin = rhs.qmin;
        qmax = rhs.qmax;
        prefix = rhs.prefix;
        mc_manager = rhs.mc_manager;
        return *this;
      }

    /// Build Parameters object with default values
    static Parameters Default()
      {
        return Parameters().NamePrefix("CF_Q3DByParent");
      }

    #define ImplSetter(__name, __type, __target) \
      Parameters __name(__type var) const  \
        { Parameters p; p.__target = var; return p; }

    ImplSetter(NamePrefix, const TString&, prefix);
    ImplSetter(Manager, AliFemtoModelManager *, mc_manager);
    ImplSetter(NBin, UInt_t, bin_count);

    #undef ImplSetter

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

    AliFemtoModelCorrFctnTrueQ3DByParent into()
      { return AliFemtoModelCorrFctnTrueQ3DByParent(*this); }

    AliFemtoModelCorrFctnTrueQ3DByParent* into_ptr()
      { return new AliFemtoModelCorrFctnTrueQ3DByParent(*this); }

    operator AliFemtoModelCorrFctnTrueQ3DByParent*()
      { return into_ptr(); }
  };

  /// Deafult parameters
  ///
  /// - Name: "CF_Q3DByParent"
  /// - Binning Paramters: (56, -0.14, 0.14)
  ///
  AliFemtoModelCorrFctnTrueQ3DByParent();

  /// Custom title
  ///
  /// Use default binning parameters
  AliFemtoModelCorrFctnTrueQ3DByParent(const char *title);

  /// Symmetric constructor
  ///
  /// Construct with nbins from -qmax to qmax in both directions
  AliFemtoModelCorrFctnTrueQ3DByParent(const char *prefix, UInt_t nbins, Double_t qmax);

  /// Custom constructor
  ///
  AliFemtoModelCorrFctnTrueQ3DByParent(const char *prefix, UInt_t nbins, Double_t qlo, Double_t qhigh);

  /// Construct with standard name
  AliFemtoModelCorrFctnTrueQ3DByParent(Int_t nbins, Double_t qmax);

  /// Construct with explicit q-range
  AliFemtoModelCorrFctnTrueQ3DByParent(Int_t nbins, Double_t qlow, Double_t qhigh);

  /// Construct from parameter object
  ///
  AliFemtoModelCorrFctnTrueQ3DByParent(const Parameters &params);

  /// copy constructor - Structure preserved - no histogram contents copied
  AliFemtoModelCorrFctnTrueQ3DByParent(const AliFemtoModelCorrFctnTrueQ3DByParent& aCorrFctn);

  /// Assignment Operator - unused
  //   AliFemtoModelCorrFctnTrueQ3DByParent& operator=(const AliFemtoModelCorrFctnTrueQ3DByParent& aCorrFctn);

  /// Destructor - histograms destroyed, ModelManager is NOT
  virtual ~AliFemtoModelCorrFctnTrueQ3DByParent();


  AliFemtoModelCorrFctnTrueQ3DByParent& operator=(const AliFemtoModelCorrFctnTrueQ3DByParent&);

  /// Set the MC model manager
  virtual void SetManager(AliFemtoModelManager *);

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  /// Finish Data
  virtual void Finish();

  virtual TList* GetOutputList();
  virtual TList* AppendOutputList(TList &);
  virtual void AddOutputObjectsTo(TCollection &);

  virtual AliFemtoCorrFctn* Clone() const;

  Double_t GetQinvTrue(AliFemtoPair*);

  //Special MC analysis for K selected by PDG code -->
  void SetKaonPDG(Bool_t aSetKaonAna);


  static Parameters Build()
    { return Parameters::Default(); }

protected:
  AliFemtoModelManager *fManager; //!<! Link back to the manager to retrieve weights

  /// Numerator made with pairs using monte-carlo generated momentum
  THnSparseF *fNumeratorGenerated;

  /// Numerator made with pairs from the same event
  THnSparseF *fNumeratorReconstructed;

  /// Denominator with the monte-carlo generated momentum
  THnSparseF *fDenominatorGenerated;

  /// Denominator with reconstructed data
  THnSparseF *fDenominatorReconstructed;
};

inline
void
AliFemtoModelCorrFctnTrueQ3DByParent::SetManager(AliFemtoModelManager *manager)
{
  fManager = manager;
}

inline
void
AliFemtoModelCorrFctnTrueQ3DByParent::Finish()
{
  // no-op
}

inline
AliFemtoCorrFctn*
AliFemtoModelCorrFctnTrueQ3DByParent::Clone() const
{
  return new AliFemtoModelCorrFctnTrueQ3DByParent(*this);
}

#endif
