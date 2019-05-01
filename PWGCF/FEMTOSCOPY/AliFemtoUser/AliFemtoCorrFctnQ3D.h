///
/// \file AliFemtoCorrFctnQ3D.h
///

#pragma once


#ifndef ALIFEMTOCORRFCTNQ3D_H
#define ALIFEMTOCORRFCTNQ3D_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"

#include <TH3F.h>
#include <TProfile3D.h>

#if __cplusplus >= 201103L
#include <tuple>
#define USE_TUPLE
#endif


/// \class AliFemtoCorrFctnQ3D
/// \brief A class to calculate 3D correlation functions for pairs
///        of identical particles vs. Bertsh-Pratt coordinates
///
/// It is recommended that you use one of the concrete classes which
/// bind this templated class to a frame-of-reference.
///
///   * AliFemtoCorrFctnQ3DLCMS
///   * AliFemtoCorrFctnQ3DBF
///   * AliFemtoCorrFctnQ3DPF
///
/// This is templated to enable custom calculation of the out, side,
/// long coordinates by writing a class (struct) with the static
/// methods `GetQ` and `FrameName` (look to the femtoQ3D namespace for
/// example implementation).
/// I would rather use std::tuple to return the three components,
/// (expanded with std::tie) but at this time we cannot support
/// c++11/ROOT6, so I use mutable references.
///
/// \author Andrew Kubera, Ohio State University <andrew.kubera@cern.ch>
///
template <typename Frame_t>
class AliFemtoCorrFctnQ3D : public AliFemtoCorrFctn {
public:

  /// Build parameter object for easy & explicit construction
  struct Build {
    TString name;
    UInt_t nbins;
    Float_t qmax;
    Bool_t use_simple_name;
    Bool_t use_tprofile;
    Bool_t single_qinv_hist;

    #define IMPL(__name, __type, __param) \
      Build __name(const __type p) const  \
        { Build b(*this); b.__param = p; return b; }

    IMPL(Name, TString, name);
    IMPL(Title, TString, name);
    IMPL(Bins, UInt_t, nbins);
    IMPL(QMax, float, qmax);
    IMPL(UseSimpleNaming, bool, use_simple_name);
    IMPL(UseTProfile, bool, use_tprofile);
    IMPL(SingleQinvHist, bool, single_qinv_hist);

    #undef IMPL

    Build()
      : name("Q3D")
      , nbins(59)
      , qmax(.295)
      , use_simple_name(true)
      , use_tprofile(false)
      , single_qinv_hist(true)
      { }

    AliFemtoCorrFctnQ3D<Frame_t> into() const
      { return AliFemtoCorrFctnQ3D<Frame_t>(name, nbins, qmax, use_simple_name, single_qinv_hist, use_tprofile); }

    AliFemtoCorrFctnQ3D<Frame_t>* into_ptr() const
      { return new AliFemtoCorrFctnQ3D<Frame_t>(name, nbins, qmax, use_simple_name, single_qinv_hist, use_tprofile); }

    operator AliFemtoCorrFctnQ3D<Frame_t>() const
      { return into(); }

    operator AliFemtoCorrFctnQ3D<Frame_t>*() const
      { return into_ptr(); }
  };

  /// Build the correlation function with parameters.
  ///
  /// \param title The title with which to give the output
  /// \param nbins The number of bins in each direction of q{out,side,long}
  /// \param QHi The limit of each axis
  /// \param simple_name - do not append title to the beginning of the name
  ///
  AliFemtoCorrFctnQ3D(const char* title,
                      const int nbins,
                      const float QHi,
                      const bool simple_name=true,
                      const bool single_qinv_hist=true,
                      const bool use_tprofile=false);

  /// Copy Constructor
  AliFemtoCorrFctnQ3D(const AliFemtoCorrFctnQ3D<Frame_t>&);

  /// Assignment - clone histograms
  AliFemtoCorrFctnQ3D& operator=(const AliFemtoCorrFctnQ3D<Frame_t>&);

  /// Deletes histograms
  virtual ~AliFemtoCorrFctnQ3D();

  virtual AliFemtoCorrFctn* Clone() const
    { return new AliFemtoCorrFctnQ3D<Frame_t>(*this); }

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* pair)
    { AddRealPair(const_cast<const AliFemtoPair&>(*pair)); }

  virtual void AddMixedPair(AliFemtoPair *pair)
    { AddMixedPair(const_cast<const AliFemtoPair&>(*pair)); }

  void AddRealPair(const AliFemtoPair &pair)
    { AddPair(pair, *fNumerator, *fNumeratorW); }

  void AddMixedPair(const AliFemtoPair &pair)
    { AddPair(pair, *fDenominator, fDenominatorW ? *fDenominatorW : *fNumeratorW); }

  /// Remove underflow-overflow contents to improve compressed file size
  virtual void Finish()
    { // no-op
    }

  /// Return denominator
  TH3& Numerator()
    { return *fNumerator; }

  /// Return denominator
  TH3& Denominator()
    { return *fDenominator; }

  /// Return bins weighed by qinv (Numerator + Denominator)
  /// -- NULL if using two histograms
  TH3* QinvW()
    { return fDenominatorW == nullptr ? fNumeratorW : nullptr; }

  /// Return numerator weighed by qinv -- NULL if using one histogram
  TH3* NumeratorW()
    { return fDenominatorW == nullptr ? nullptr : fNumeratorW; }

  /// Return denominator weighed by qinv -- NULL if using one histogram
  TH3* DenominatorW()
    { return fDenominatorW; }

  virtual TList* GetOutputList()
  {
    TList *list = new TList();
    AddOutputObjectsTo(*list);
    return list;
  }

  virtual void AddOutputObjectsTo(TCollection &dest)
  {
    dest.Add(fNumerator);
    dest.Add(fDenominator);

    dest.Add(fNumeratorW);

    if (fDenominatorW) {
      dest.Add(fDenominatorW);
    }
  }

  /// Load 3D q-vector components into variables
  static void GetQ(const AliFemtoPair &pair, double &out, double &side, double &lon);

#ifdef USE_TUPLE
  /// Return 3D q-vector components
  std::tuple<double, double, double> GetQ(const AliFemtoPair &pair) const
  {
    double x, y, z;
    GetQ(pair, x, y, z);
    return std::make_tuple(x, y, z);
  }
#endif

protected:

  void AddPair(const AliFemtoPair &pair, TH3& dest, TH3& qinv)
  {
    if (fPairCut && !fPairCut->Pass(&pair)) {
      return;
    }

    // auto [qout, qside, qlong] = Frame_t::GetQ(pair); // maybe someday...
    double qout, qside, qlong;
    Frame_t::GetQ(pair, qout, qside, qlong);

    Int_t bin = dest.FindBin(qout, qlong, qside);
    if (!(dest.IsBinOverflow(bin) or dest.IsBinUnderflow(bin))) {
      dest.Fill(qout, qside, qlong);
      qinv.Fill(qout, qside, qlong, pair.QInv());
    }
  }

  TH3F* fNumerator;     ///<!< Numerator
  TH3F* fDenominator;   ///<!< Denominator
  TH3* fNumeratorW;     ///<!< Qinv-Weighted numerator
  TH3* fDenominatorW;   ///<!< Qinv-Weighted denominator
};


template <typename T>
AliFemtoCorrFctnQ3D<T>::AliFemtoCorrFctnQ3D(const char* title,
                                            const int nbins,
                                            const float QHi,
                                            const bool simple_name,
                                            const bool single_qinv_hist,
                                            const bool use_tprofile)
  : AliFemtoCorrFctn()
  , fNumerator(nullptr)
  , fDenominator(nullptr)
  , fNumeratorW(nullptr)
  , fDenominatorW(nullptr)
{
  TString hist_title = TString::Format("%s (Frame=%s); q_{out} (GeV); q_{side} (GeV); q_{long} (GeV)", title, T::FrameName());
  TString hist_name = simple_name ? "" : title;

  auto new_3d_hist = [&] (const TString name, const TString htitle)
    {
      return new TH3F(simple_name ? name : name + " " + title,
                      htitle + hist_title,
                      nbins, -QHi, QHi,
                      nbins, -QHi, QHi,
                      nbins, -QHi, QHi);
    };

  fNumerator = new_3d_hist("Num", "Numerator");
  fDenominator = new_3d_hist("Den", "Denominator");

  auto qinv_member_builder = [&] (const TString name, const TString htitle) -> TH3*
    {
      if (!use_tprofile) {
        return new_3d_hist(name, htitle);
      }

      return new TProfile3D(simple_name ? name : name + " " + title,
                            htitle + hist_title,
                            nbins, -QHi, QHi,
                            nbins, -QHi, QHi,
                            nbins, -QHi, QHi);
    };

  // note: non-weighted histograms do not need Sumw2 - save space and time by not enabling
  if (single_qinv_hist) {
    fNumeratorW = qinv_member_builder("QinvW", "Q_{inv} Weights (divide by Num + Den)");
    fNumeratorW->Sumw2();
  } else {
    fNumeratorW = qinv_member_builder("NumQinvW", "Q_{inv} Weighted Numerator");
    fDenominatorW = qinv_member_builder("DenQinvW", "Q_{inv} Weighted Denominator");
    fNumeratorW->Sumw2();
    fDenominatorW->Sumw2();
  }
}

template <typename T>
AliFemtoCorrFctnQ3D<T>::AliFemtoCorrFctnQ3D(const AliFemtoCorrFctnQ3D<T>& orig)
  : AliFemtoCorrFctn(orig)
  , fNumerator(new TH3F(*orig.fNumerator))
  , fDenominator(new TH3F(*orig.fDenominator))
  , fNumeratorW(static_cast<TH3*>(orig.fNumeratorW->Clone()))
  , fDenominatorW(static_cast<TH3*>(orig.fDenominatorW ? orig.fDenominatorW->Clone() : nullptr))
{
}

template <typename T>
AliFemtoCorrFctnQ3D<T>&
AliFemtoCorrFctnQ3D<T>::operator=(const AliFemtoCorrFctnQ3D<T> &rhs)
{
  if (this != &rhs) {
    AliFemtoCorrFctn::operator=(rhs);
    *fNumerator = *rhs.fNumerator;
    *fDenominator = *rhs.fDenominator;
    delete fNumeratorW;
    fNumeratorW = static_cast<TH3*>(rhs.fNumeratorW->Clone());

    delete fDenominatorW;
    fDenominatorW = rhs.fDenominatorW ? static_cast<TH3*>(rhs.fDenominatorW->Clone()) : nullptr;
  }

  return *this;
}

template <typename T>
AliFemtoCorrFctnQ3D<T>::~AliFemtoCorrFctnQ3D()
{
  delete fNumerator;
  delete fDenominator;
  delete fNumeratorW;
  delete fDenominatorW;
}

template <typename T>
AliFemtoString
AliFemtoCorrFctnQ3D<T>::Report()
{
  // Construct the report
  AliFemtoString report
    = AliFemtoString("Bertsch-Pratt 3D Correlation Function")
    + Form(" (Frame = %s) Report:\n", T::FrameName())
    + Form("Number of entries in numerator:\t%E\n", fNumerator->GetEntries())
    + Form("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());

  if (fPairCut) {
    report += "Here is the PairCut specific to this CorrFctn\n";
    report += fPairCut->Report();
  } else {
    report += "No PairCut specific to this CorrFctn\n";
  }

  return report;
}

#undef SINGLE_WQINV

/*
#ifdef USE_TUPLE
template <typename T>
void
AliFemtoCorrFctnQ3D<T>::GetQ(const AliFemtoPair &pair, double &x, double &y, double &z)
{
  std::tie(x, y, z) = static_cast<T*>(this)->GetQ(pair);
}
#endif
*/

struct AliFemtoCorrFctnQ3DLCMS : public AliFemtoCorrFctnQ3D<AliFemtoCorrFctnQ3DLCMS> {

  typedef AliFemtoCorrFctnQ3D<AliFemtoCorrFctnQ3DLCMS> Super;
  typedef Super::Build Build;

  AliFemtoCorrFctnQ3DLCMS(const Super::Build &b)
    : Super(b.into())
    { }

  AliFemtoCorrFctnQ3DLCMS(const char *name, int x, float y, bool b=true)
    : Super(name, x, y, b)
    { }

  static void GetQ(const AliFemtoPair &pair, double &x, double &y, double &z)
  {
    x = pair.QOutCMS();
    y = pair.QSideCMS();
    z = pair.QLongCMS();
  }

  static const char* FrameName()
    { return "LCMS"; }
};

struct AliFemtoCorrFctnQ3DPF : public AliFemtoCorrFctnQ3D<AliFemtoCorrFctnQ3DPF> {

  typedef AliFemtoCorrFctnQ3D<AliFemtoCorrFctnQ3DPF> Super;
  typedef Super::Build Build;

  AliFemtoCorrFctnQ3DPF(const Super::Build &b)
    : Super(b.into())
    { }

  AliFemtoCorrFctnQ3DPF(const char *name, int x, float y, bool b=true)
    : Super(name, x, y, b)
    { }

  static void GetQ(const AliFemtoPair &pair, double &x, double &y, double &z)
  {
    x = pair.QOutPf();
    y = pair.QSidePf();
    z = pair.QLongPf();
  }

  static const char* FrameName()
    { return "PF"; }
};

struct AliFemtoCorrFctnQ3DBF : public AliFemtoCorrFctnQ3D<AliFemtoCorrFctnQ3DBF> {

  typedef AliFemtoCorrFctnQ3D<AliFemtoCorrFctnQ3DBF> Super;
  typedef Super::Build Build;

  AliFemtoCorrFctnQ3DBF(const Super::Build &b)
    : Super(b.into())
    { }

  AliFemtoCorrFctnQ3DBF(const char *name, int x, float y, bool b=true)
    : Super(name, x, y, b)
    { }

  static void GetQ(const AliFemtoPair &pair, double &x, double &y, double &z)
  {
    x = pair.QOutBf();
    y = pair.QSideBf();
    z = pair.QLongBf();
  }

  static const char* FrameName()
    { return "BF"; }
};

#undef USE_TUPLE

#endif
