////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoModelCorrFctnNonIdDR - correlation function for non-identical      //
// particles which uses k* as a function variable. Stores the correlation     //
// function separately for positive and negative signs of k* projections into //
// out, side and long directions, enabling the calculations of double ratios  //
// Uses pair weight to simulate the model correlation function.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELCORRFCTNNONIDDR_H
#define ALIFEMTOMODELCORRFCTNNONIDDR_H

#include "TH1D.h"
#include "AliFemtoModelCorrFctn.h"

class AliFemtoModelCorrFctnNonIdDR : public AliFemtoModelCorrFctn {
public:
  AliFemtoModelCorrFctnNonIdDR(char* title, const int& nbins, const float& QinvLo, const float& QinvHi);
  AliFemtoModelCorrFctnNonIdDR(const AliFemtoModelCorrFctnNonIdDR& aCorrFctn);
  virtual ~AliFemtoModelCorrFctnNonIdDR();

  AliFemtoModelCorrFctnNonIdDR& operator=(const AliFemtoModelCorrFctnNonIdDR& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  virtual AliFemtoModelCorrFctn* Clone();

  virtual TList* GetOutputList();
  void Write();

private:
  TH1D *fNumTOutP;     // True numerator for pair with positive k*out
  TH1D *fNumTOutN;     // True numerator for pair with negative k*out
  TH1D *fNumTSideP;    // True numerator for pair with positive k*side
  TH1D *fNumTSideN;    // True numerator for pair with negative k*side
  TH1D *fNumTLongP;    // True numerator for pair with positive k*long
  TH1D *fNumTLongN;    // True numerator for pair with negative k*long

  TH1D *fNumFOutP;     // Fake numerator for pair with positive k*out
  TH1D *fNumFOutN;     // Fake numerator for pair with negative k*out
  TH1D *fNumFSideP;    // Fake numerator for pair with positive k*side
  TH1D *fNumFSideN;    // Fake numerator for pair with negative k*side
  TH1D *fNumFLongP;    // Fake numerator for pair with positive k*long
  TH1D *fNumFLongN;    // Fake numerator for pair with negative k*long

  TH1D *fDenOutP;     // Denominator for pair with positive k*out
  TH1D *fDenOutN;     // Denominator for pair with negative k*out
  TH1D *fDenSideP;    // Denominator for pair with positive k*side
  TH1D *fDenSideN;    // Denominator for pair with negative k*side
  TH1D *fDenLongP;    // Denominator for pair with positive k*long
  TH1D *fDenLongN;    // Denominator for pair with negative k*long

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnNonIdDR, 1)
#endif
};


#endif

