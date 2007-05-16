////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctnSource - the class for correlation function which   ///
/// uses the model framework and weight generation and saves the generated   ///
/// emission source                                                          ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoModelCorrFctnSource_hh
#define AliFemtoModelCorrFctnSource_hh

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

  virtual AliFemtoModelCorrFctn* Clone();

protected:

  TH1D *fHistROut;
  TH1D *fHistRSide;
  TH1D *fHistRLong;
  TH1D *fHistRStar;
  TH1D *fHistdNdR;

private:

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnSource, 1)
#endif
};

#endif
