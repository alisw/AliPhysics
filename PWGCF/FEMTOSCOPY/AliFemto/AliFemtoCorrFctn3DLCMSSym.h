///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoCorrFctn3DLCMSSym: a class to calculate 3D correlation        //
// for pairs of identical particles vs. Bertsh-Pratt coordinates.        //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTN3DLCMS_H
#define ALIFEMTOCORRFCTN3DLCMS_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3F.h"

class AliFemtoCorrFctn3DLCMSSym : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctn3DLCMSSym(char* title, const int& nbins, const float& QHi);
  AliFemtoCorrFctn3DLCMSSym(const AliFemtoCorrFctn3DLCMSSym& aCorrFctn);
  virtual ~AliFemtoCorrFctn3DLCMSSym();

  AliFemtoCorrFctn3DLCMSSym& operator=(const AliFemtoCorrFctn3DLCMSSym& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair( AliFemtoPair* aPair);
  virtual void AddMixedPair( AliFemtoPair* aPair);

  virtual void Finish();

  TH3F* Numerator();
  TH3F* Denominator();

  void WriteOutHistos();
  virtual TList* GetOutputList();

private:

  TH3F* fNumerator;         // numerator
  TH3F* fDenominator;       // denominator

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctn3DLCMSSym, 1)
#endif
};

inline  TH3F* AliFemtoCorrFctn3DLCMSSym::Numerator(){return fNumerator;}
inline  TH3F* AliFemtoCorrFctn3DLCMSSym::Denominator(){return fDenominator;}

#endif

