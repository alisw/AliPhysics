////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctn - the base class for correlation function which    ///
/// uses the model framework and weight generation                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELCORRFCTN_H
#define ALIFEMTOMODELCORRFCTN_H

#include "AliFemtoCorrFctn.h"
class AliFemtoPair;
class AliFemtoModelManager;
class TH1D;

class AliFemtoModelCorrFctn: public AliFemtoCorrFctn {

public:
  AliFemtoModelCorrFctn();
  AliFemtoModelCorrFctn(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi);
  AliFemtoModelCorrFctn(const AliFemtoModelCorrFctn& aCorrFctn);
  virtual ~AliFemtoModelCorrFctn();

  AliFemtoModelCorrFctn& operator=(const AliFemtoModelCorrFctn& aCorrFctn);

  virtual void ConnectToManager(AliFemtoModelManager *aManager);

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPir);

  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);
  virtual void Finish();

  virtual void Write();

  virtual AliFemtoModelCorrFctn* Clone();

  AliFemtoAnalysis* HbtAnalysis(){return fyAnalysis;};
  void SetAnalysis(AliFemtoAnalysis* aAnalysis);

protected:
  AliFemtoModelManager *fManager; // Link back to the managet to get the weights
  
  TH1D *fNumeratorTrue; // Numerator made with pairs from the same event
  TH1D *fNumeratorFake; // Numerator made with pairs from different events (mixed pairs)
  TH1D *fDenominator;   // Denominator made with mixed pairs

private:

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctn, 1)
#endif
};

#endif
