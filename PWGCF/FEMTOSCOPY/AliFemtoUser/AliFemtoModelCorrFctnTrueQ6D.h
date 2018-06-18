///
/// \file AliFemtoModelCorrFctnTrueQ6D.h
///

#pragma once

#ifndef ALIFEMTOMODELCORRFCTN_TRUEQ6D_H
#define ALIFEMTOMODELCORRFCTN_TRUEQ6D_H


#include "AliFemtoCorrFctn.h"
#include <THnSparse.h>
#include "AliFemtoConfigObject.h"

// #include "AliFemtoCorr.h"

class AliFemtoPair;
class AliFemtoModelManager;
class TRandom;


/// \class AliFemtoModelCorrFctnTrueQ6D
/// \brief Correlation function storing map from true MonteCarlo generated
///        momentum to "measured" momentum in LCMS frame.
///
/// \author Andrew Kubera, The Ohio State University <andrew.michael.kubera@cern.ch>
///
class AliFemtoModelCorrFctnTrueQ6D : public AliFemtoCorrFctn {
public:

  struct Builder {
    UInt_t bin_count;
    Double_t qmin;
    Double_t qmax;
    TString title;
    AliFemtoModelManager *mc_manager;

    Builder()
      : bin_count(56)
      , qmin(-0.24)
      , qmax(0.24)
      , title("CF_TrueQ6D")
      , mc_manager(NULL)
    {
    }

    Builder(AliFemtoConfigObject);

    Builder Title(const TString& title) const
    {
      Builder b(*this);
      b.title = title;
      return b;
    }

    Builder Manager(AliFemtoModelManager *manager) const
    {
      Builder b(*this);
      b.mc_manager = manager;
      return b;
    }

    Builder Q(double q) const
    {
      Builder b(*this);
      b.qmax = std::abs(q);
      b.qmin = -b.qmax;
      return b;
    }

    Builder Qmin(double q) const
    {
      Builder b(*this);
      b.qmin = q;
      return b;
    }

    Builder Qmax(double q) const
    {
      Builder b(*this);
      b.qmax = q;
      return b;
    }

    operator AliFemtoModelCorrFctnTrueQ6D*()
    {
      return new AliFemtoModelCorrFctnTrueQ6D(*this);
    }
  };

  static Builder Build()
  {
    return Builder();
  }

  /// Deafult parameters
  ///
  /// - Name: "CF_TrueQ6D"
  /// - Binning Paramters: (56, -0.14, 0.14)
  ///
  AliFemtoModelCorrFctnTrueQ6D();

  /// Custom title
  ///
  /// Use default binning parameters
  AliFemtoModelCorrFctnTrueQ6D(const TString &title);

  /// Symmetric constructor
  ///
  /// Construct with nbins from -qmax to qmax in both directions
  AliFemtoModelCorrFctnTrueQ6D(const TString &title,
                               UInt_t nbins,
                               Double_t qmax);

  /// Custom constructor
  ///
  AliFemtoModelCorrFctnTrueQ6D(const TString &title,
                               UInt_t nbins,
                               Double_t aQinvLo,
                               Double_t aQinvHi);

  /// Construct from parameter object
  ///
  AliFemtoModelCorrFctnTrueQ6D(const Builder &);

  /// copy constructor - Structure preserved - no histogram contents copied
  AliFemtoModelCorrFctnTrueQ6D(const AliFemtoModelCorrFctnTrueQ6D&);

  /// Assignment Operator - unused
  //   AliFemtoModelCorrFctnTrueQ6D& operator=(const AliFemtoModelCorrFctnTrueQ6D& aCorrFctn);

  /// Destructor - histograms destroyed, ModelManager is NOT
  virtual ~AliFemtoModelCorrFctnTrueQ6D();

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

  /// Histogram of data
  THnSparseS *fHistogram;

  /// Random number generator used for randomizing order of pair momentums
  TRandom *fRng;
};

inline
void
AliFemtoModelCorrFctnTrueQ6D::SetManager(AliFemtoModelManager *manager)
{
  fManager = manager;
  std::cout << "fManager set to " << fManager << "\n";
}

inline
void
AliFemtoModelCorrFctnTrueQ6D::Finish()
{
  // no-op
}

inline
AliFemtoCorrFctn*
AliFemtoModelCorrFctnTrueQ6D::Clone() const
{
  return new AliFemtoModelCorrFctnTrueQ6D(*this);
}

#endif
