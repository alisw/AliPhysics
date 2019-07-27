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

    Double_t qout_range_min;
    Double_t qout_range_max;
    Double_t qside_range_min;
    Double_t qside_range_max;
    Double_t qlong_range_min;
    Double_t qlong_range_max;

    bool ignore_zeromass;

    Builder()
      : bin_count(59)
      , qmin(-0.295)
      , qmax(0.295)
      , bin_method(kRecLSOGenOSL)
      , title("CF_TrueQ6D")
      , qout_range_min(0.0)
      , qout_range_max(0.0)
      , qside_range_min(0.0)
      , qside_range_max(0.0)
      , qlong_range_min(0.0)
      , qlong_range_max(0.0)
      , ignore_zeromass(true)
    {
    }

    Builder(const Builder&);
    Builder& operator=(const Builder&);
    Builder(AliFemtoConfigObject);

    Builder Title(const TString& title) const
    {

      Builder b(*this);
      b.title = title;
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
        return b; }                                   \
      Builder __name(double r[2]) const {             \
        Builder b(*this);                             \
        b. __target ## _range_min = r[0];             \
        b. __target ## _range_max = r[1];             \
        return b; }

    CREATE_SETTER_METHOD(QoutRange, qout);
    CREATE_SETTER_METHOD(QsideRange, qside);
    CREATE_SETTER_METHOD(QlongRange, qlong);

    #undef CREATE_SETTER_METHOD


    #define CREATE_SETTER_METHOD(__name, __type, __target) \
      Builder __name(__type x) const                       \
        { Builder b(*this); b. __target = x; return b; }

    CREATE_SETTER_METHOD(NBins, Int_t, bin_count);
    CREATE_SETTER_METHOD(IgnoreZeroMass, Bool_t, ignore_zeromass);
    CREATE_SETTER_METHOD(Binning, BinMethod, bin_method);

    #undef CREATE_SETTER_METHOD

    AliFemtoModelCorrFctnTrueQ6D* IntoCF()
      { return new AliFemtoModelCorrFctnTrueQ6D(*this); }

    operator AliFemtoModelCorrFctnTrueQ6D*()
      { return IntoCF(); }
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
  AliFemtoModelCorrFctnTrueQ6D(const HistType &);

  /// Building using pointer
  AliFemtoModelCorrFctnTrueQ6D(HistType *&);

  /// Custom title
  ///
  /// Use default binning parameters
  ///
  AliFemtoModelCorrFctnTrueQ6D(const TString &title);

  /// Symmetric constructor
  ///
  /// Construct with nbins from -qmax to qmax in both directions
  ///
  AliFemtoModelCorrFctnTrueQ6D(const TString &title,
                               UInt_t nbins,
                               Double_t qmax);

  /// Custom constructor
  ///
  AliFemtoModelCorrFctnTrueQ6D(const TString &prefix,
                               UInt_t nbins,
                               Double_t aQinvLo,
                               Double_t aQinvHi,
                               BinMethod binning=kRecGenOSL);

  /// Asymmetric constructor
  ///
  AliFemtoModelCorrFctnTrueQ6D(const TString &prefix,
                               Int_t nbins_out,
                               Double_t qout_lo,
                               Double_t qout_hi,
                               Int_t nbins_side,
                               Double_t qside_lo,
                               Double_t qside_hi,
                               Int_t nbins_long,
                               Double_t qlong_lo,
                               Double_t qlong_hi,
                               BinMethod binning=kRecGenOSL);

  /// Big Asymmetric constructor
  ///
  AliFemtoModelCorrFctnTrueQ6D(const TString &prefix,
                               Int_t nbins_out,
                               Double_t qout_lo,
                               Double_t qout_hi,
                               Int_t nbins_side,
                               Double_t qside_lo,
                               Double_t qside_hi,
                               Int_t nbins_long,
                               Double_t qlong_lo,
                               Double_t qlong_hi,
                               Int_t nbins_out_true,
                               Double_t qout_lo_true,
                               Double_t qout_hi_true,
                               Int_t nbins_side_true,
                               Double_t qside_lo_true,
                               Double_t qside_hi_true,
                               Int_t nbins_long_true,
                               Double_t qlong_lo_true,
                               Double_t qlong_hi_true,
                               BinMethod binning=kRecGenOSL);

  AliFemtoModelCorrFctnTrueQ6D(const TString &prefix,
                               const std::vector<double> &obins,
                               const std::vector<double> &sbins,
                               const std::vector<double> &lbins,
                               BinMethod binning=kRecGenOSL);

  /// Construct from parameter object
  ///
  AliFemtoModelCorrFctnTrueQ6D(const Builder &);

  /// copy constructor - Structure preserved - no histogram contents copied
  AliFemtoModelCorrFctnTrueQ6D(const AliFemtoModelCorrFctnTrueQ6D&);

  /// Assignment Operator - unused
  AliFemtoModelCorrFctnTrueQ6D& operator=(const AliFemtoModelCorrFctnTrueQ6D&);

  /// Destructor - histograms destroyed, ModelManager is NOT
  virtual ~AliFemtoModelCorrFctnTrueQ6D();

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

  /// try to guess how bins are layed out by axis name
  static BinMethod GuessBinMethod(const HistType &hist);

protected:

  /// Histogram of data
  HistType *fHistogram;

  BinMethod fBinMethod;
  Bool_t fIgnoreZeroMassParticles;

  std::pair<double, double> fQlimits[3];

  void UpdateQlimits();

private:
  void AddPair(const AliFemtoParticle &, const AliFemtoParticle &);
};

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
