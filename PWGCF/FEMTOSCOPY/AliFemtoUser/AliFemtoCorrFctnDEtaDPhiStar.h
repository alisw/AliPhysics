///
/// \file AliFemtoCorrFctnDEtaDPhiStar.h
///

#pragma once

#ifndef ALIFEMTOCORRFCTNDETADPHISTAR_H
#define ALIFEMTOCORRFCTNDETADPHISTAR_H


#include "AliFemtoCorrFctn.h"

#include <TH2F.h>

/// \class AliFemtoCorrFctnDEtaDPhiStar
/// \breif Correlation function for dphistar
class AliFemtoCorrFctnDEtaDPhiStar : public AliFemtoCorrFctn {
public:

  struct Params {
    TString suffix;
    double radius,
           phi_max,
           eta_max;
    UInt_t nbins_phi,
           nbins_eta;

    Params()
      : suffix("DEtaDPhistar")
      , radius(1.2)
      , phi_max(0.2)
      , eta_max(0.2)
      , nbins_phi(144)
      , nbins_eta(144)
      { }

    #define IMPL_SETTER(__name, __type, __dest) \
      Params __name(const __type &p) \
        { Params res(*this); res. __dest = p; return res; }

    IMPL_SETTER(Suffix, TString, suffix);
    IMPL_SETTER(Radius, double, radius);
    IMPL_SETTER(Phi, double, phi_max);
    IMPL_SETTER(Eta, double, eta_max);
    IMPL_SETTER(NbinsPhi, UInt_t, nbins_phi);
    IMPL_SETTER(NbinsEta, UInt_t, nbins_eta);

    #undef IMPL_SETTER

    operator AliFemtoCorrFctnDEtaDPhiStar() const
      { return AliFemtoCorrFctnDEtaDPhiStar(*this); }

    AliFemtoCorrFctnDEtaDPhiStar into() const
      { return AliFemtoCorrFctnDEtaDPhiStar(*this); }

    AliFemtoCorrFctnDEtaDPhiStar* into_ptr() const
      { return new AliFemtoCorrFctnDEtaDPhiStar(*this); }
  };

  AliFemtoCorrFctnDEtaDPhiStar();
  AliFemtoCorrFctnDEtaDPhiStar(const char *suffix, double radius=1.2);
  AliFemtoCorrFctnDEtaDPhiStar(const Params &);
  AliFemtoCorrFctnDEtaDPhiStar(const AliFemtoCorrFctnDEtaDPhiStar &);

  AliFemtoCorrFctnDEtaDPhiStar& operator=(const AliFemtoCorrFctnDEtaDPhiStar &);

  virtual ~AliFemtoCorrFctnDEtaDPhiStar();

  AliFemtoString Report()
    { return ""; }

  virtual void Finish()
    { }

  virtual void EventBegin(const AliFemtoEvent *ev)
    {
      fMagneticField = ev->MagneticField();
      if (fabs(fMagneticField) < 1e-10) {
        fMagneticField *= 1e13;
      }

      AliFemtoCorrFctn::EventBegin(ev);
    }

  virtual TList* GetOutputList();

  virtual void AddRealPair(AliFemtoPair *);
  virtual void AddMixedPair(AliFemtoPair *);

  virtual AliFemtoCorrFctn* Clone() const
    { return new AliFemtoCorrFctnDEtaDPhiStar(*this); }

protected:

  void AddPair(TH2&, const AliFemtoPair&);

  /// Numerator of dEta dPhi* function
  TH2F *fDPhiDEtaNumerator;

  /// Denominator of dEta dPhi* function
  TH2F *fDPhiDEtaDenominator;

  double fRadius;

  double fMagneticField;
};


#endif
