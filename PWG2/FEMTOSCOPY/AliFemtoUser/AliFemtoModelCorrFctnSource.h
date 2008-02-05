////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctnSource - the class for correlation function which   ///
/// uses the model framework and weight generation and saves the generated   ///
/// emission source                                                          ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELCORRFCTNSOURCE_H
#define ALIFEMTOMODELCORRFCTNSOURCE_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelCorrFctn.h"

class AliFemtoModelCorrFctnSource: public AliFemtoModelCorrFctn {

public:
  AliFemtoModelCorrFctnSource();
  AliFemtoModelCorrFctnSource(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi);
  AliFemtoModelCorrFctnSource(const AliFemtoModelCorrFctnSource& aCorrFctn);
  virtual ~AliFemtoModelCorrFctnSource();
  
  AliFemtoModelCorrFctnSource& operator=(const AliFemtoModelCorrFctnSource& aCorrFctn);

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPir);

  virtual void Write();
  virtual TList* GetOutputList();

  virtual AliFemtoModelCorrFctn* Clone();

protected:

  TH1D *fHistROut;     // Distribution of Rout
  TH1D *fHistRSide;    // Distribution of Rside
  TH1D *fHistRLong;    // Distribution of Rlong
  TH1D *fHistRStar;    // Distribution of RStar
  TH1D *fHistdNdR;     // Distribution of RStar weighted by Jacobian 

private:

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnSource, 1)
#endif
};

#endif
