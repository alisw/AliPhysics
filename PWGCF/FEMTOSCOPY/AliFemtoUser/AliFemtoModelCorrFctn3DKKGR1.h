////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelCorrFctn3DKKGR1 - class for 3D KK correlation function       ///
/// gleb.romanenko@cern.ch                                                   ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELCORRFCTN3DKKGR1_H
#define ALIFEMTOMODELCORRFCTN3DKKGR1_H

//#include "AliFemtoCorrFctn.h"
#include "AliFemtoModelCorrFctn.h"
class AliFemtoPair;
class AliFemtoModelManager;
class TH1D;
class TH2D;
class TH3D;

class AliFemtoModelCorrFctn3DKKGR1: public AliFemtoModelCorrFctn {

public:
  AliFemtoModelCorrFctn3DKKGR1();
  AliFemtoModelCorrFctn3DKKGR1(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi);
  AliFemtoModelCorrFctn3DKKGR1(const AliFemtoModelCorrFctn3DKKGR1& aCorrFctn);
  virtual ~AliFemtoModelCorrFctn3DKKGR1();

  AliFemtoModelCorrFctn3DKKGR1& operator=(const AliFemtoModelCorrFctn3DKKGR1& aCorrFctn);

  virtual void ConnectToManager(AliFemtoModelManager *aManager);

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void EventBegin(const AliFemtoEvent* aEvent);
  virtual void EventEnd(const AliFemtoEvent* aEvent);
  virtual void Finish();

  virtual void Write();
  virtual TList* GetOutputList();
  virtual AliFemtoModelCorrFctn3DKKGR1* Clone() const { return new AliFemtoModelCorrFctn3DKKGR1(*this); }

  void SetFillkT(bool fillkT){fFillkT = fillkT;}

  void SetFillQgenQrecInv(bool FillQgenQrecInv){fFillQgenQrecInv = FillQgenQrecInv;}
  void SetFillQgenQrecOutSideLong(bool FillQgenQrecOutSideLong){fFillQgenQrecOutSideLong = FillQgenQrecOutSideLong;}
  void SetFillDpTpTvspT(bool FillDpTpTvspT){fFillDpTpTvspT = FillDpTpTvspT;}
  
  void SetSpecificPairCut(AliFemtoPairCut* aCut);
  
  Double_t GetQinvTrue(AliFemtoPair*);
  Double_t GetQoutTrue(AliFemtoPair*);
  Double_t GetQsideTrue(AliFemtoPair*);
  Double_t GetQlongTrue(AliFemtoPair*);

  //Special MC analysis for K selected by PDG code -->
  void SetKaonPDG(Bool_t aSetKaonAna);

protected:
  AliFemtoModelManager *fManager; // Link back to the manager to get the weights

  TH3D *fNumeratorTrue; // Numerator made with pairs from the same event
  TH3D *fNumeratorFake; // Numerator made with pairs from different events (mixed pairs)
  TH3D *fDenominator;   // Denominator made with mixed pairs

  TH3D *fNumeratorTrueIdeal; // Numerator made with pairs (true qinv) from the same event
  TH3D *fNumeratorFakeIdeal; // Numerator made with pairs (true qinv) from different events (mixed pairs)
  TH3D *fDenominatorIdeal;   // Denominator made with mixed pairs (true qinv)

  TH2D *fQgenQrec; // Qinv true (generated) vs. Qinv reconstructed
  TH2D *fQgenQrecOut; // Qout true (generated) vs. Qout reconstructed
  TH2D *fQgenQrecSide; // Qside true (generated) vs. Qside reconstructed
  TH2D *fQgenQrecLong; // Qlong true (generated) vs. Qlong reconstructed

  TH2D *fDpTpTvspT; // (pT_sim - pT_gen)/pT_sim vs. pT_sim


private:

  //Special MC analysis for K selected by PDG code -->
  Bool_t fKaonPDG;
  
  bool fFillkT;
  bool fFillQgenQrecInv;
  bool fFillQgenQrecOutSideLong;
  bool fFillDpTpTvspT;
  const int fNbbPairs = 21;
  TH1D *fkTdists[21]; // histograms with kT distributions for different BB pairs
  double GetParentsKt(AliFemtoPair *pair);
  int GetPairNumber(AliFemtoPair *pair); // returns pair code

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelCorrFctn3DKKGR1, 1);
  /// \endcond
#endif
};

#endif
