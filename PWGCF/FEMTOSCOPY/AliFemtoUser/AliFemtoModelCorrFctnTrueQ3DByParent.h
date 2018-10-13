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
    TString title;
    AliFemtoModelManager *mc_manager;

    /// Build Parameters object with default values
    static Parameters Default()
    {
      return {
        56, -0.14, 0.14,  // histogram bin-count & range
        "CF_Q3DByParent", // title
        NULL              // pointer to MC manager
      };
    }

    #define ImplSetter(__name, __type, __target) \
      Parameters __name(__type var) const  \
        { Parameters p; p.__target = var; return p; }

    ImplSetter(Title, const TString&, title);
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

    /*
    operator AliFemtoModelCorrFctnTrueQ3DByParent()
      { return into(); }
    */

    operator AliFemtoModelCorrFctnTrueQ3DByParent*()
      { return new AliFemtoModelCorrFctnTrueQ3DByParent(*this); }
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
  AliFemtoModelCorrFctnTrueQ3DByParent(const char *title, UInt_t nbins, Double_t qmax);

  /// Custom constructor
  ///
  AliFemtoModelCorrFctnTrueQ3DByParent(const char *title, UInt_t nbins, Double_t aQinvLo, Double_t aQinvHi);

  /// Construct
  AliFemtoModelCorrFctnTrueQ3DByParent(Int_t nbins, Double_t aQinvLo, Double_t aQinvHi);

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

  virtual AliFemtoCorrFctn* Clone() const;

  Double_t GetQinvTrue(AliFemtoPair*);

  //Special MC analysis for K selected by PDG code -->
  void SetKaonPDG(Bool_t aSetKaonAna);

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

  /// Random number generator used for randomizing order of pair momentums
  TRandom *fRng;
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
