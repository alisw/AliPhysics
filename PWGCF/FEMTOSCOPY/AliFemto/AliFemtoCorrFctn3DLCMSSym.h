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
//#include "TArrayD.h"

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
  TH3F* NumeratorW();//Weighed by qinv
  TH3F* DenominatorW();


  void WriteOutHistos();
  virtual TList* GetOutputList();

  void SetUseLCMS(int);
  int  GetUseLCMS();

private:

  TH3F* fNumerator;         // numerator
  TH3F* fDenominator;       // denominator
  TH3F* fNumeratorW;         // numerator
  TH3F* fDenominatorW;       // denominator

  int    fUseLCMS;             // 0 - Use PRF, 1 - Use LCMS

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctn3DLCMSSym, 1);
  /// \endcond
#endif
};

inline  TH3F* AliFemtoCorrFctn3DLCMSSym::Numerator(){return fNumerator;}
inline  TH3F* AliFemtoCorrFctn3DLCMSSym::Denominator(){return fDenominator;}
inline  TH3F* AliFemtoCorrFctn3DLCMSSym::NumeratorW(){return fNumeratorW;}
inline  TH3F* AliFemtoCorrFctn3DLCMSSym::DenominatorW(){return fDenominatorW;}
#endif
