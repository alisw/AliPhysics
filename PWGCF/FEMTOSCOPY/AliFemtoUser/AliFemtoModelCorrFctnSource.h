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

#include "TH2D.h"
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
  virtual AliFemtoModelCorrFctn* Clone() const { return new AliFemtoModelCorrFctnSource(*this); }

  void SetUseRPSelection(unsigned short aRPSel);
protected:

  TH1D *fHistROut;     // Distribution of Rout
  TH1D *fHistRSide;    // Distribution of Rside
  TH1D *fHistRLong;    // Distribution of Rlong
  TH1D *fHistRStar;    // Distribution of RStar
  TH1D *fHistdNdR;     // Distribution of RStar weighted by Jacobian 
  TH2D *fHistNumWS;    // Weight spread for numerator
  TH2D *fHistDenWS;    // Weight spread for denominator

private:

  unsigned short fUseRPSelection;  // The pair cut uses RP selection

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnSource, 1)
#endif
};

#endif
