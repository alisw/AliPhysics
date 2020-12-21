/******************************************************************************/
/*                                                                            */
/*  AliFemtoModelCorrFctnKK -   numerator and denominator                     */
/*      histograms from the true Monte Carlo KchKch tracks                    */
/*         Konstantin.Mikhaylov@cern.ch                                       */
/*                                                                            */
/******************************************************************************/

#ifndef ALIFEMTOMODELCORRFCTNKK_H
#define ALIFEMTOMODELCORRFCTNKK_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelCorrFctn.h"
#include "TH2D.h"
#include "TH1D.h"

class AliFemtoModelCorrFctnKK: public AliFemtoModelCorrFctn {

public:
  // default constructor:
  AliFemtoModelCorrFctnKK();
  // constructor with some parameters(histogram bins,low,high and so on):
  AliFemtoModelCorrFctnKK(const char *title,
			  Int_t aNbins,
			  Double_t aQinvLo,
			  Double_t aQinvHi);
  // copy constructor:
  AliFemtoModelCorrFctnKK(const AliFemtoModelCorrFctnKK& aCorrFctn);
  // destructor:
  virtual ~AliFemtoModelCorrFctnKK();

  AliFemtoModelCorrFctnKK& operator=(const AliFemtoModelCorrFctnKK& aCorrFctn);

  virtual void ConnectToManager(AliFemtoModelManager *aManager);

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

 
  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);
  virtual void Finish();

  virtual TList* GetOutputList();
  virtual void Write();

  virtual AliFemtoModelCorrFctn* Clone() const;
  
  //void SetFillkT(bool fillkT){fFillkT = fillkT;} //i do not need this
   
  Double_t GetQinvTrue(AliFemtoPair*);
  
  //Special MC analysis for K selected by PDG code -->
  void SetKaonPDG(Bool_t aSetKaonAna);

protected:
  
  AliFemtoModelManager *fManager; // Link back to the manager to get the weights

  TH1D *fNumeratorTrue; // Numerator made with pairs from the same event
  TH1D *fNumeratorFake; // Numerator made with pairs from different events (mixed pairs)
  TH1D *fDenominator;   // Denominator made with mixed pairs

  TH1D *fNumeratorTrueIdeal; // Numerator made with pairs (true qinv) from the same event
  TH1D *fNumeratorFakeIdeal; // Numerator made with pairs (true qinv) from different events (mixed pairs)
  TH1D *fDenominatorIdeal;   // Denominator made with mixed pairs (true qinv)

  TH2D *fQgenQrec; // Qinv true (generated) vs. Qinv reconstructed

  TH1D *fdP;//delta_p/p
  TH1D *fdPt;//delta_p_T/p_T
  TH2D *fdPtvsPt; //delta_p_T/p_T versus p_T

private:
  
  //Special MC analysis for K selected by PDG code -->
  Bool_t fKaonPDG;

  Double_t fP1x, fP2x;//do not take the same particles
  /* 
    bool fFillkT;
    int fNbbPairs = 21;
    TH1D *fkTdists[21]; // histograms with kT distributions for different BB pairs
    double GetParentsKt(AliFemtoPair *pair);
  */
    int GetPairNumber(AliFemtoPair *pair); // returns pair code
  
#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelCorrFctnKK, 1);
  /// \endcond
#endif
};

#endif
