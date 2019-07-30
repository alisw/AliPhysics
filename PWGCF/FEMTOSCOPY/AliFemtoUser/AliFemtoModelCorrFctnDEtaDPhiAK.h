///
/// \file AliFemtoModelCorrFctnDEtaDPhiAK.h
///

#ifndef ALIFEMTOMODELCORRFCTNDETADPHIAK_H
#define ALIFEMTOMODELCORRFCTNDETADPHIAK_H

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelCorrFctn.h"

/// \class AliFemtoModelCorrFctnDEtaDPhiAK
/// \brief Correlation Function solving problems
///
class AliFemtoModelCorrFctnDEtaDPhiAK : public AliFemtoModelCorrFctn {
public:

  /// Builder-pattern struct
  struct Params {
    TString suffix;
    UInt_t nbins_phi,
           nbins_eta;

    Params()
      : suffix("")
      , nbins_phi(50)
      , nbins_eta(50)
      { }

    AliFemtoModelCorrFctnDEtaDPhiAK*
    into_ptr() const
      { return new AliFemtoModelCorrFctnDEtaDPhiAK(*this); }

    operator AliFemtoModelCorrFctnDEtaDPhiAK() const
      { return AliFemtoModelCorrFctnDEtaDPhiAK(*this); }

    AliFemtoModelCorrFctnDEtaDPhiAK
    into() const
      { return AliFemtoModelCorrFctnDEtaDPhiAK(*this); }

    #define IMPL_SETTER(__name, __type, __dest) \
      Params __name(const __type &p) \
        { Params res(*this); res. __dest = p; return res; }

    IMPL_SETTER(Suffix, char*, suffix);
    IMPL_SETTER(NBinsPhi, UInt_t, nbins_phi);
    IMPL_SETTER(NBinsEta, UInt_t, nbins_eta);

    #undef IMPL_SETTER
  };

  /// Construct with name suffix and bin-count
  AliFemtoModelCorrFctnDEtaDPhiAK(const char* suffix, const UInt_t phibins, const UInt_t etabins);
  AliFemtoModelCorrFctnDEtaDPhiAK(const Params &);
  AliFemtoModelCorrFctnDEtaDPhiAK(const AliFemtoModelCorrFctnDEtaDPhiAK&);
  virtual ~AliFemtoModelCorrFctnDEtaDPhiAK();

  AliFemtoModelCorrFctnDEtaDPhiAK& operator=(const AliFemtoModelCorrFctnDEtaDPhiAK&);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair *pair)
    { AddRealPair(*pair); }
  virtual void AddMixedPair(AliFemtoPair *pair)
    { AddMixedPair(*pair); }

  virtual void AddRealPair(const AliFemtoPair &);
  virtual void AddMixedPair(const AliFemtoPair &);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoModelCorrFctn* Clone() const
    { return new AliFemtoModelCorrFctnDEtaDPhiAK(*this); }

protected:

  TH2D *fDPhiDEtaNumeratorTrue;  //< Numerator of dEta dPhi true function
  TH2D *fDPhiDEtaNumeratorFake;  //< Numerator of dEta dPhi fake function
  TH2D *fDPhiDEtaDenominator;    //< Denominator of dEta dPhi function

  TH2D *fDPhiDEtaColNumerator;   //< Numerator of colinear dEta dPhi function
  TH2D *fDPhiDEtaColDenominator; //< Denominator of colinear dEta dPhi function

  TH1D *fDPhiNumeratorTrue;      //< Numerator of dPhi true correlation
  TH1D *fDPhiNumeratorFake;      //< Numerator of dPhi fake correlation
  TH1D *fDPhiDenominator;        //< Denominator of dPhi correlation

  TH1D *fDCosNumeratorTrue;      //< Numerator of colinearity true correlation
  TH1D *fDCosNumeratorFake;      //< Numerator of colinearity fake correlation
  TH1D *fDCosDenominator;        //< Denominator of colinearity correlation

  TH2D *fDPhiPtNumerator;        //< Numerator of dPhi correlation vs. Pt min
  TH2D *fDPhiPtDenominator;      //< Denominator of dPhi correlation vs. Pt min

  TH2D *fDCosPtNumerator;        //< Numerator of colinearity correlation vs. Pt min
  TH2D *fDCosPtDenominator;      //< Denominator of colinearity correlation vs. Pt min
};


#endif
