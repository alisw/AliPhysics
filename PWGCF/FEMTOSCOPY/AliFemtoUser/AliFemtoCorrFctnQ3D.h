///
/// \file AliFemtoCorrFctnQ3D.h
///

#pragma once


#ifndef ALIFEMTOCORRFCTNQ3D_H
#define ALIFEMTOCORRFCTNQ3D_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"

#include <TH3I.h>
#include <TH3F.h>


// preprocessor flag to enable using ONE histogram to store
#define SINGLE_WQINV


/// namespace for storing BP (Q{out,side,long}) calculator classes
namespace femtoQ3D {

  /// \class FrameLCMS
  /// \breif Calculate
  struct FrameLCMS {
    static void GetQ(const AliFemtoPair &pair,
                     double &x, double &y, double &z)
    {
      x = pair.QOutCMS();
      y = pair.QSideCMS();
      z = pair.QLongCMS();
    }

    static const char* FrameName()
      { return "LCMS"; }
  };


  /// \class FramePF
  struct FramePF {
    static void GetQ(const AliFemtoPair &pair,
                     double &x, double &y, double &z)
    {
      x = pair.QOutPf();
      y = pair.QSidePf();
      z = pair.QLongPf();
    }

    static const char* FrameName()
      { return "PF"; }
  };

  /// \class FrameBF
  struct FrameBF {
    static void GetQ(const AliFemtoPair &pair,
                     double &x, double &y, double &z)
    {
      x = pair.QOutBf();
      y = pair.QSideBf();
      z = pair.QLongBf();
    }

    static const char* FrameName()
      { return "BF"; }
  };
}


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

    #define IMPL(__name, __type, __param) \
      Build __name(const __type p) const  \
        { Build b(*this); b.__param = p; return b; }

    IMPL(Name, TString, name);
    IMPL(Title, TString, name);
    IMPL(Bins, UInt_t, nbins);
    IMPL(QMax, float, qmax);
    IMPL(UseSimpleNaming, bool, use_simple_name);

    #undef IMPL

    Build()
      : name("Q3D")
      , nbins(59)
      , qmax(.295)
      , use_simple_name(true)
      { }

    AliFemtoCorrFctnQ3D<Frame_t> into() const
      { return AliFemtoCorrFctnQ3D<Frame_t>(name, nbins, qmax, use_simple_name); }

    AliFemtoCorrFctnQ3D<Frame_t>* into_ptr() const
      { return new AliFemtoCorrFctnQ3D<Frame_t>(name, nbins, qmax, use_simple_name); }

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
                      const bool simple_name=true);

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

#ifdef SINGLE_WQINV
  void AddRealPair(const AliFemtoPair &pair)
    { AddPair(pair, *fNumerator, *fQinvW); }

  void AddMixedPair(const AliFemtoPair &pair)
    { AddPair(pair, *fDenominator, *fQinvW); }
#else
  void AddRealPair(const AliFemtoPair &pair)
    { AddPair(pair, *fNumerator, *fNumeratorW); }

  void AddMixedPair(const AliFemtoPair &pair)
    { AddPair(pair, *fDenominator, *fDenominatorW); }
#endif

  /// No-op
  virtual void Finish()
    { }

  /// Return denominator
  TH3& Numerator()
    { return *fNumerator; }

  /// Return denominator
  TH3& Denominator()
    { return *fDenominator; }

#ifdef SINGLE_WQINV

  /// Return weighed by qinv
  TH3& QinvW()
    { return *fQinvW; }

#else

  /// Return numerator weighed by qinv
  TH3& NumeratorW()
    { return *fNumeratorW; }

  /// Return denominator weighed by qinv
  TH3& DenominatorW()
    { return *fDenominatorW; }
#endif

  virtual TList* GetOutputList()
  {
    TList *list = new TList();
    list->Add(fNumerator);
    list->Add(fDenominator);

#ifdef SINGLE_WQINV
    list->Add(fQinvW);
#else
    list->Add(fNumeratorW);
    list->Add(fDenominatorW);
#endif

    return list;
  }

protected:

  void AddPair(const AliFemtoPair &pair, TH3& dest, TH3& qinv)
  {
    if (fPairCut && !fPairCut->Pass(&pair)) {
      return;
    }

    // auto [qout, qside, qlong] = Frame_t::GetQ(pair); // maybe someday...
    double qout, qside, qlong;
    Frame_t::GetQ(pair, qout, qside, qlong);
    dest.Fill(qout, qside, qlong);
    qinv.Fill(qout, qside, qlong, pair.QInv());
  }

  TH3I* fNumerator;     ///<!< Numerator
  TH3I* fDenominator;   ///<!< Denominator
#ifdef SINGLE_WQINV
  TH3F* fQinvW;         ///<!< Qinv-Weighted histogram
#else
  TH3F* fNumeratorW;    ///<!< Qinv-Weighted numerator
  TH3F* fDenominatorW;  ///<!< Qinv-Weighted denominator
#endif
};

template <typename T>
AliFemtoCorrFctnQ3D<T>::AliFemtoCorrFctnQ3D(const char* title,
                                            const int nbins,
                                            const float QHi,
                                            const bool simple_name)
  : AliFemtoCorrFctn()
  , fNumerator(nullptr)
  , fDenominator(nullptr)
#ifdef SINGLE_WQINV
  , fQinvW(nullptr)
#else
  , fNumeratorW(nullptr)
  , fDenominatorW(nullptr)
#endif
{
  TString hist_title = TString::Format("%s (Frame=%s); q_{out} (GeV); q_{side} (GeV); q_{long} (GeV)", title, T::FrameName());
  TString hist_name = simple_name ? "" : title;

  fNumerator = new TH3I(simple_name ? "Num" : (TString("Num") + title).Data(),
                        "Numerator " + hist_title,
                        nbins, -QHi, QHi,
                        nbins, -QHi, QHi,
                        nbins, -QHi, QHi);

  fDenominator = new TH3I(simple_name ? "Den" : (TString("Den") + title).Data(),
                          "Denominator " + hist_title,
                          nbins, -QHi, QHi,
                          nbins, -QHi, QHi,
                          nbins, -QHi, QHi);

#ifdef SINGLE_WQINV
  fQinvW = new TH3F(simple_name ? "QinvW" : (TString("QinvW") + title).Data(),
                         "Q_{inv} Weights (divide by Num + Den)" + hist_title,
                         nbins, -QHi, QHi,
                         nbins, -QHi, QHi,
                         nbins, -QHi, QHi);

  fQinvW->Sumw2();
#else
  fNumeratorW = new TH3F(simple_name ? "NumWqinv" : (TString("NumWqinv") + title).Data(),
                         "Q_{inv} Weighted Numerator " + hist_title,
                         nbins, -QHi, QHi,
                         nbins, -QHi, QHi,
                         nbins, -QHi, QHi);

  fDenominatorW = new TH3F(simple_name ? "DenWqinv" : (TString("DenWqinv") + title).Data(),
                           "Q_{inv} Weighted Denominator " + hist_title,
                           nbins, -QHi, QHi,
                           nbins, -QHi, QHi,
                           nbins, -QHi, QHi);

  // note: non-weighted histograms do not need Sumw2 - save space and time by not enabling
  fNumeratorW->Sumw2();
  fDenominatorW->Sumw2();
#endif
}

template <typename T>
AliFemtoCorrFctnQ3D<T>::AliFemtoCorrFctnQ3D(const AliFemtoCorrFctnQ3D<T>& orig)
  : AliFemtoCorrFctn(orig)
  , fNumerator(new TH3I(*orig.fNumerator))
  , fDenominator(new TH3I(*orig.fDenominator))
#ifdef SINGLE_WQINV
  , fQinvW(new TH3F(*orig.fQinvW))
#else
  , fNumeratorW(new TH3F(*orig.fNumeratorW))
  , fDenominatorW(new TH3F(*orig.fDenominatorW))
#endif
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
#ifdef SINGLE_WQINV
    *fQinvW = *rhs.fQinvW;
#else
    *fNumeratorW = *rhs.fNumeratorW;
    *fDenominatorW = *rhs.fDenominatorW;
#endif
  }

  return *this;
}

template <typename T>
AliFemtoCorrFctnQ3D<T>::~AliFemtoCorrFctnQ3D()
{
  delete fNumerator;
  delete fDenominator;
#ifdef SINGLE_WQINV
  delete fQinvW;
#else
  delete fNumeratorW;
  delete fDenominatorW;
#endif
}

template <typename T>
AliFemtoString
AliFemtoCorrFctnQ3D<T>::Report()
{
  // Construct the report
  TString report = TString::Format("Bertsch-Pratt 3D Correlation Function (Frame = %s) Report:\n", T::FrameName())
                 + Form("Number of entries in numerator:\t%E\n", fNumerator->GetEntries())
                 + Form("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());

  if (fPairCut) {
    report += "Here is the PairCut specific to this CorrFctn\n";
    report += fPairCut->Report();
  } else {
    report += "No PairCut specific to this CorrFctn\n";
  }

  return AliFemtoString(report.Data());
}

#undef SINGLE_WQINV


typedef AliFemtoCorrFctnQ3D<femtoQ3D::FrameLCMS> AliFemtoCorrFctnQ3DLCMS;
typedef AliFemtoCorrFctnQ3D<femtoQ3D::FrameBF> AliFemtoCorrFctnQ3DBF;
typedef AliFemtoCorrFctnQ3D<femtoQ3D::FramePF> AliFemtoCorrFctnQ3DPF;

#endif
