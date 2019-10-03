///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoCorrFctn3DPRF_qosl_q: a class to calculate 3D correlation        //
// for pairs of particles in PRF
//Written by Ashutosh Kumar Pandey (ashutosh.kumar.pandey@cern.ch) for Source Imaging//
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoCorrFctn3DPRF_qosl_q_H
#define AliFemtoCorrFctn3DPRF_qosl_q_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"
#include <THnSparse.h>
#include "TH3F.h"

class AliFemtoCorrFctn3DPRF_qosl_q : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctn3DPRF_qosl_q(const char* title, const int& nbins, const float& QHi);
  AliFemtoCorrFctn3DPRF_qosl_q(const AliFemtoCorrFctn3DPRF_qosl_q& aCorrFctn);
  virtual ~AliFemtoCorrFctn3DPRF_qosl_q();

  AliFemtoCorrFctn3DPRF_qosl_q& operator=(const AliFemtoCorrFctn3DPRF_qosl_q& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair( AliFemtoPair* aPair);
  virtual void AddMixedPair( AliFemtoPair* aPair);

  virtual void Finish();

  TH3F* Numerator();
  TH3F* Denominator();
  TH3F* NumeratorW();//Weighed by qinv
  TH3F* DenominatorW();
  THnSparse* NumeratorS(); //numerator
  THnSparse* DenominatorS(); //denominator

  void WriteOutHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctn3DPRF_qosl_q(*this); }

private:

  TH3F* fNumerator;         // numerator
  TH3F* fDenominator;       // denominator
  TH3F* fNumeratorW;         // numerator
  TH3F* fDenominatorW;       // denominator
  THnSparse* fNumS; //numerator
  THnSparse* fDenS; //denominator
#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctn3DPRF_qosl_q, 1)
#endif
};

inline  TH3F* AliFemtoCorrFctn3DPRF_qosl_q::Numerator(){return fNumerator;}
inline  TH3F* AliFemtoCorrFctn3DPRF_qosl_q::Denominator(){return fDenominator;}
inline  TH3F* AliFemtoCorrFctn3DPRF_qosl_q::NumeratorW(){return fNumeratorW;}
inline  TH3F* AliFemtoCorrFctn3DPRF_qosl_q::DenominatorW(){return fDenominatorW;}
inline  THnSparse* AliFemtoCorrFctn3DPRF_qosl_q::NumeratorS(){return fNumS;}
inline  THnSparse* AliFemtoCorrFctn3DPRF_qosl_q::DenominatorS(){return fDenS;}
#endif

