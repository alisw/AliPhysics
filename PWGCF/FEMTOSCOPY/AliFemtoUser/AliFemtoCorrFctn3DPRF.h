///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoCorrFctn3DPRF: a class to calculate 3D correlation            //
// for pairs of particles in PRF                                         //
//  (by Ashutosh Kumar Pandey)                                           //
//                                                                       //
//////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTN3DPRF_H
#define ALIFEMTOCORRFCTN3DPRF_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3F.h"

class AliFemtoCorrFctn3DPRF : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctn3DPRF(const char* title, const int& nbins, const float& QHi);
  AliFemtoCorrFctn3DPRF(const AliFemtoCorrFctn3DPRF& aCorrFctn);
  virtual ~AliFemtoCorrFctn3DPRF();

  AliFemtoCorrFctn3DPRF& operator=(const AliFemtoCorrFctn3DPRF& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair( AliFemtoPair* aPair);
  virtual void AddMixedPair( AliFemtoPair* aPair);

  virtual void Finish();

  TH3F* Numerator();
  TH3F* Denominator();
  TH3F* NumeratorW();//Weighted by |k*|
  TH3F* DenominatorW();


  void WriteOutHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctn3DPRF(*this); }

private:

  TH3F* fNumerator;         // numerator
  TH3F* fDenominator;       // denominator
  TH3F* fNumeratorW;         // numerator
  TH3F* fDenominatorW;       // denominator
#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctn3DPRF, 1)
#endif
};

inline  TH3F* AliFemtoCorrFctn3DPRF::Numerator(){return fNumerator;}
inline  TH3F* AliFemtoCorrFctn3DPRF::Denominator(){return fDenominator;}
inline  TH3F* AliFemtoCorrFctn3DPRF::NumeratorW(){return fNumeratorW;}
inline  TH3F* AliFemtoCorrFctn3DPRF::DenominatorW(){return fDenominatorW;}
#endif

