////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctnWithWeights - the base class for correlation function which    ///
/// uses the model framework and weight generation                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELCORRFCTNWITHWEIGHTS_H
#define ALIFEMTOMODELCORRFCTNWITHWEIGHTS_H

#include "AliFemtoCorrFctn.h"
class AliFemtoPair;
class AliFemtoModelManager;
class TH1D;
class TH2D;

class AliFemtoModelCorrFctnWithWeights: public AliFemtoCorrFctn {

public:
  AliFemtoModelCorrFctnWithWeights(TH2D *filter1, TH2D *filter2);
  AliFemtoModelCorrFctnWithWeights(const char *title, TH2D *filter1, TH2D *filter2, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi);
  AliFemtoModelCorrFctnWithWeights(const AliFemtoModelCorrFctnWithWeights& aCorrFctn);
  virtual ~AliFemtoModelCorrFctnWithWeights();

  AliFemtoModelCorrFctnWithWeights& operator=(const AliFemtoModelCorrFctnWithWeights& aCorrFctn);

  virtual void ConnectToManager(AliFemtoModelManager *aManager);

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);
  virtual void Finish();

  virtual TList* GetOutputList();
  virtual void Write();

  virtual AliFemtoModelCorrFctnWithWeights* Clone() const;

  Double_t GetQinvTrue(AliFemtoPair*);

protected:
  TH2D *filterHist1;
  TH2D *filterHist2;

  AliFemtoModelManager *fManager; // Link back to the manager to get the weights

  TH1D *fNumeratorTrue; // Numerator made with pairs from the same event
  TH1D *fNumeratorFake; // Numerator made with pairs from different events (mixed pairs)
  TH1D *fDenominator;   // Denominator made with mixed pairs

  TH1D *fNumeratorTrueIdeal; // Numerator made with pairs (true qinv) from the same event
  TH1D *fNumeratorFakeIdeal; // Numerator made with pairs (true qinv) from different events (mixed pairs)
  TH1D *fDenominatorIdeal;   // Denominator made with mixed pairs (true qinv)

  TH2D *fQgenQrec; // Qinv true (generated) vs. Qinv reconstructed


private:

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelCorrFctnWithWeights, 1);
  /// \endcond
#endif
};

#endif
