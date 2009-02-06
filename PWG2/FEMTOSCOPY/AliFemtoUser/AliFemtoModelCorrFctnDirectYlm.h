////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoModelCorrFctnDirectYlm - the class for correlation function which  //
// uses the model framework and weight generation and saves the correlation   //
// function directly in spherical harmonics                                   //
// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELCORRFCTNDIRECTYLM_H
#define ALIFEMTOMODELCORRFCTNDIRECTYLM_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoCorrFctnDirectYlm.h"

class AliFemtoModelCorrFctnDirectYlm: public AliFemtoModelCorrFctn {

public:
  AliFemtoModelCorrFctnDirectYlm();
  AliFemtoModelCorrFctnDirectYlm(const char *title, Int_t aMaxL, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi);
  AliFemtoModelCorrFctnDirectYlm(const AliFemtoModelCorrFctnDirectYlm& aCorrFctn);
  virtual ~AliFemtoModelCorrFctnDirectYlm();
  
  AliFemtoModelCorrFctnDirectYlm& operator=(const AliFemtoModelCorrFctnDirectYlm& aCorrFctn);

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPir);

  virtual void Finish();
  virtual void Write();
  virtual TList* GetOutputList();

  virtual AliFemtoModelCorrFctn* Clone();

protected:

  AliFemtoCorrFctnDirectYlm* fCYlmTrue;     // True Correlation function in spherical harmonics
  AliFemtoCorrFctnDirectYlm* fCYlmFake;     // Fake Correlation function in spherical harmonics

private:

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnDirectYlm, 1)
#endif
};

#endif
