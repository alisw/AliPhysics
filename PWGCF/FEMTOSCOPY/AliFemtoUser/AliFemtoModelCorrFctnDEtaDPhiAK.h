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

class AliFemtoModelCorrFctnDEtaDPhiAK : public AliFemtoModelCorrFctn {
public:

  /// Construct with name suffix and bin-count
  AliFemtoModelCorrFctnDEtaDPhiAK(const char* suffix, const UInt_t phibins, const UInt_t etabins);
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

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnDEtaDPhiAK, 1)
#endif
};


#endif
