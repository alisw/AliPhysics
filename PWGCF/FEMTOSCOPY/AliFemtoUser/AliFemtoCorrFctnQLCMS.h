///
/// \file AliFemtoCorrFctnQLCMS.h
///

#pragma once


#ifndef ALIFEMTOCORRFCTNQLCMS_H
#define ALIFEMTOCORRFCTNQLCMS_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"

#include <TH1F.h>
#include <TProfile.h>


/// \class AliFemtoCorrFctnQLCMS
/// \brief A class to calculate correlation function of LCMS
///        three-momentum difference
///
/// $ |q_{lcms}| = \sqrt{(p_{1x} - p_{2x})^2 + (p_{1y} - p_{2y})^2 + q_{long,LCMS}^2} $
///
/// \author Andrew Kubera, Ohio State University <andrew.kubera@cern.ch>
///
class AliFemtoCorrFctnQLCMS : public AliFemtoCorrFctn {
public:

  /// Build parameter object for easy & explicit construction
  struct Parameters {
    TString prefix;
    TString suffix;
    UInt_t nbins;
    Float_t qmin;
    Float_t qmax;

    #define IMPL(__name, __type, __param) \
      Parameters __name(const __type p) const  \
        { Parameters b(*this); b.__param = p; return b; }

    IMPL(Prefix, TString, prefix);
    IMPL(Suffix, TString, suffix);
    IMPL(Bins, UInt_t, nbins);
    IMPL(QMin, float, qmin);
    IMPL(QMax, float, qmax);

    #undef IMPL

    Parameters()
      : prefix("")
      , suffix("QLCMS")
      , nbins(72)
      , qmin(0.0)
      , qmax(1.0)
      { }

    AliFemtoCorrFctnQLCMS into() const
      { return AliFemtoCorrFctnQLCMS(*this); }

    AliFemtoCorrFctnQLCMS* into_ptr() const
      { return new AliFemtoCorrFctnQLCMS(*this); }

  };

  /// Build the correlation function with parameters.
  ///
  /// \param title The title with which to give the output
  /// \param nbins The number of bins in each direction of q{out,side,long}
  /// \param QHi The limit of each axis
  /// \param simple_name - do not append title to the beginning of the name
  ///
  AliFemtoCorrFctnQLCMS(const char* prefix,
                        const char* suffix,
                        const int nbins,
                        const float QHi);

  /// Build the correlation function with parameters.
  ///
  /// \param title The title with which to give the output
  /// \param nbins The number of bins in each direction of q{out,side,long}
  /// \param QHi The limit of each axis
  /// \param simple_name - do not append title to the beginning of the name
  ///
  AliFemtoCorrFctnQLCMS(const char* prefix,
                        const char* suffix,
                        const int nbins,
                        float QLo,
                        float QHi);

  /// Build with parameter struct
  AliFemtoCorrFctnQLCMS(const Parameters &);

  /// Copy Constructor
  AliFemtoCorrFctnQLCMS(const AliFemtoCorrFctnQLCMS&);

  /// Assignment - clone histograms
  AliFemtoCorrFctnQLCMS& operator=(const AliFemtoCorrFctnQLCMS&);

  /// Deletes histograms
  virtual ~AliFemtoCorrFctnQLCMS();

  virtual AliFemtoCorrFctnQLCMS* Clone() const
    { return new AliFemtoCorrFctnQLCMS(*this); }

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* pair)
    { AddRealPair(const_cast<const AliFemtoPair&>(*pair)); }

  virtual void AddMixedPair(AliFemtoPair *pair)
    { AddMixedPair(const_cast<const AliFemtoPair&>(*pair)); }

  void AddRealPair(const AliFemtoPair &pair)
    { AddPair(pair, *fNumerator, *fQinv); }

  void AddMixedPair(const AliFemtoPair &pair)
    { AddPair(pair, *fDenominator, *fQinv); }

  /// Remove underflow-overflow contents to improve compressed file size
  virtual void Finish()
    { // no-op
    }

  /// Return numerator histogram
  TH1& Numerator()
    { return *fNumerator; }

  /// Return denominator histogram
  TH1& Denominator()
    { return *fDenominator; }

  /// Return qinv profile
  TH1& Qinv()
    { return *fQinv; }

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
      dest.Add(fQinv);
    }

  /// Return default parameter object
  static Parameters Build()
    { return Parameters(); }

protected:

  void AddPair(const AliFemtoPair &pair, TH1& dest, TProfile& qinv)
    {
      if (fPairCut && !fPairCut->Pass(&pair)) {
        return;
      }

      const double
        qlong = pair.QLongCMS(),
        qt_sqrd = pair.FourMomentumDiff().Perp2(),
        qlcms = std::sqrt(qt_sqrd + qlong * qlong);


      dest.Fill(qlcms);
      qinv.Fill(qlcms, pair.QInv());
    }

  TH1F* fNumerator;    ///<!< Numerator
  TH1F* fDenominator;  ///<!< Denominator
  TProfile* fQinv;     ///<!< Qinv vs QLCMS profile
};


#endif
