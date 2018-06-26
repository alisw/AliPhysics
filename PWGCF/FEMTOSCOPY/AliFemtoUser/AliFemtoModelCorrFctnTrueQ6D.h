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

  typedef THnSparseI HistType;

  enum BinMethod {
    kGenRecLSO,
    kGenRecOSL,
    kGenLSORecOSL,
    kRecGenLSO,
    kRecGenOSL,
    kRecLSOGenOSL,
    kGroupedAxisOSL,
    kGroupedAxisLSO
  };

  struct Builder {
    UInt_t bin_count;
    Double_t qmin;
    Double_t qmax;
    BinMethod bin_method;
    TString title;
    AliFemtoModelManager *mc_manager;

    Double_t qout_range_min;
    Double_t qout_range_max;
    Double_t qside_range_min;
    Double_t qside_range_max;
    Double_t qlong_range_min;
    Double_t qlong_range_max;

    bool ignore_zeromass;

    Builder()
      : bin_count(120)
      , qmin(-0.3)
      , qmax(0.3)
      , bin_method(kGenRecOSL)
      , title("CF_TrueQ6D")
      , mc_manager(NULL)
      , qout_range_min(0.0)
      , qout_range_max(0.0)
      , qside_range_min(0.0)
      , qside_range_max(0.0)
      , qlong_range_min(0.0)
      , qlong_range_max(0.0)
      , ignore_zeromass(true)
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

    #define CREATE_SETTER_METHOD(__name, __target)    \
      Builder __name(double low, double high) const { \
        Builder b(*this);                             \
        b. __target ## _range_min = low;              \
        b. __target ## _range_max = high;             \
        return b; }

    CREATE_SETTER_METHOD(QoutRange, qout);
    CREATE_SETTER_METHOD(QsideRange, qside);
    CREATE_SETTER_METHOD(QlongRange, qlong);

    #undef CREATE_SETTER_METHOD


    #define CREATE_SETTER_METHOD(__name, __type, __target) \
      Builder __name(__type x) const {                     \
        Builder b(*this); b. __target = x; return b; }

    CREATE_SETTER_METHOD(NBins, Int_t, bin_count);
    CREATE_SETTER_METHOD(IgnoreZeroMass, Bool_t, ignore_zeromass);
    CREATE_SETTER_METHOD(Binning, BinMethod, bin_method);

    #undef CREATE_SETTER_METHOD


    AliFemtoModelCorrFctnTrueQ6D* IntoCF()
    {
      return new AliFemtoModelCorrFctnTrueQ6D(*this);
    }

    operator AliFemtoModelCorrFctnTrueQ6D*()
    {
      return IntoCF();

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

  /// Construct from data
  AliFemtoModelCorrFctnTrueQ6D(const HistType &, AliFemtoModelManager *m=nullptr);

  /// Building using pointer
  AliFemtoModelCorrFctnTrueQ6D(HistType *&, AliFemtoModelManager *m=nullptr);

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

  /// set q range
  void SetQrange(const double out[2], const double s[2], const double l[2]);

protected:
  AliFemtoModelManager *fManager; //!<! Link back to the manager to retrieve weights

  /// Histogram of data
  HistType *fHistogram;

  /// Random number generator used for randomizing order of pair momentums
  TRandom *fRng;

  BinMethod fBinMethod;
  Bool_t fIgnoreZeroMassParticles;

  std::pair<double, double> fQlimits[3];

  void UpdateQlimits();

private:
  void AddPair(const AliFemtoParticle &, const AliFemtoParticle &);
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
