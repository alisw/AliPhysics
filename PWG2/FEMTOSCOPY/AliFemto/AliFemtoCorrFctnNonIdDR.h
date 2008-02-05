////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnNonIdDR - correlation function for non-identical particles //
// uses k* as a function variable. Stores the correlation function separately //
// for positive and negative signs of k* projections into out, side and long  //
// directions, enabling the calculations of double ratios                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOCORRFCTNNONIDDR_H
#define ALIFEMTOCORRFCTNNONIDDR_H

#include "TH1D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnNonIdDR : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnNonIdDR(char* title, const int& nbins, const float& QinvLo, const float& QinvHi);
  AliFemtoCorrFctnNonIdDR(const AliFemtoCorrFctnNonIdDR& aCorrFctn);
  virtual ~AliFemtoCorrFctnNonIdDR();

  AliFemtoCorrFctnNonIdDR& operator=(const AliFemtoCorrFctnNonIdDR& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  virtual TList* GetOutputList();
  void Write();

private:
  TH1D *fNumOutP;     // Numerator for pair with positive k*out
  TH1D *fNumOutN;     // Numerator for pair with negative k*out
  TH1D *fNumSideP;    // Numerator for pair with positive k*side
  TH1D *fNumSideN;    // Numerator for pair with negative k*side
  TH1D *fNumLongP;    // Numerator for pair with positive k*long
  TH1D *fNumLongN;    // Numerator for pair with negative k*long

  TH1D *fDenOutP;     // Denominator for pair with positive k*out
  TH1D *fDenOutN;     // Denominator for pair with negative k*out
  TH1D *fDenSideP;    // Denominator for pair with positive k*side
  TH1D *fDenSideN;    // Denominator for pair with negative k*side
  TH1D *fDenLongP;    // Denominator for pair with positive k*long
  TH1D *fDenLongN;    // Denominator for pair with negative k*long

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnNonIdDR, 1)
#endif
};


#endif

