////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctnTrueQ - the class for correlation function which    ///
/// uses the model framework and weight generation and saves the correlation ///
/// function in true qinv                                                    ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELCORRFCTNTRUEQ_H
#define ALIFEMTOMODELCORRFCTNTRUEQ_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelCorrFctn.h"

class AliFemtoModelCorrFctnTrueQ: public AliFemtoModelCorrFctn {

public:
  AliFemtoModelCorrFctnTrueQ();
  AliFemtoModelCorrFctnTrueQ(const char *title, Int_t aNbins, Double_t aQinvHi);
  AliFemtoModelCorrFctnTrueQ(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi);
  AliFemtoModelCorrFctnTrueQ(const AliFemtoModelCorrFctnTrueQ& aCorrFctn);
  virtual ~AliFemtoModelCorrFctnTrueQ();

  AliFemtoModelCorrFctnTrueQ& operator=(const AliFemtoModelCorrFctnTrueQ& aCorrFctn);

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPir);

  virtual void Write();
  virtual TList* GetOutputList();
  virtual AliFemtoModelCorrFctn* Clone() const;

protected:

  TH1D *fTrueNum;           // Numerator in true q
  TH1D *fTrueDen;           // Denominator in true q

private:

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnTrueQ, 1)
#endif
};

#endif
