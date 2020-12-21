////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelBPLCMSCorrFctn - the class for correlation function which   ///
/// uses the model framework and weight generation and calculated the 3D     ///
/// correlation function in the Bertsh-Pratt LCMS system                     ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELBPLCMSCORRFCTNKK_H
#define ALIFEMTOMODELBPLCMSCORRFCTNKK_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3D.h"

class AliFemtoModelBPLCMSCorrFctnKK : public AliFemtoModelCorrFctn {
 public:
  AliFemtoModelBPLCMSCorrFctnKK()  :
    AliFemtoModelCorrFctn(),
    fNumerator3DTrue(0),
    fNumerator3DFake(0),
    fDenominator3D(0),
    fQinvHisto(0),
    fPairCut(0),
    fUseRPSelection(0),
    fNumerator3DTrueIdeal(0),
    fNumerator3DFakeIdeal(0),
    fDenominator3DIdeal(0){}
  AliFemtoModelBPLCMSCorrFctnKK(const char* title, const int& nbins, const float& QLo, const float& QHi);
  AliFemtoModelBPLCMSCorrFctnKK(const AliFemtoModelBPLCMSCorrFctnKK& aCorrFctn);
  virtual ~AliFemtoModelBPLCMSCorrFctnKK();

  AliFemtoModelBPLCMSCorrFctnKK& operator=(const AliFemtoModelBPLCMSCorrFctnKK& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* pair);
  virtual void AddMixedPair(AliFemtoPair* pair);

  virtual void Finish();

  virtual void Write();
  virtual TList* GetOutputList();
  virtual AliFemtoModelCorrFctn* Clone() const { return new AliFemtoModelBPLCMSCorrFctnKK(*this); }

  void SetSpecificPairCut(AliFemtoPairCut* aCut);
  void SetUseRPSelection(unsigned short aRPSel);

  //ml
  Double_t GetQinvTrue(AliFemtoPair*);
  Double_t GetQoutTrue(AliFemtoPair*);
  Double_t GetQsideTrue(AliFemtoPair*);
  Double_t GetQlongTrue(AliFemtoPair*);


  virtual AliFemtoModelCorrFctn* Clone();

protected:
  TH3D* fNumerator3DTrue;            // 3D Numerator with pairs from same event only
  TH3D* fNumerator3DFake;            // 3D Numerator with pairs from mixed events
  TH3D* fDenominator3D;              // 3D Denominator with the weight of 1.0

  TH3D* fQinvHisto;                  // Averag qinv histogram


  //ml

  AliFemtoPairCut* fPairCut;    //! this is a PairCut specific to THIS CorrFctn, not the Analysis
  unsigned short fUseRPSelection;  // The pair cut uses RP selection

  TH3D *fNumerator3DTrueIdeal; // Numerator made with pairs (true qosl) from the same event
  TH3D *fNumerator3DFakeIdeal; // Numerator made with pairs (true qosl) from different events (mixed pairs)
  TH3D *fDenominator3DIdeal;   // Denominator made with mixed pairs (true qosl)




#ifdef __ROOT__
  ClassDef(AliFemtoModelBPLCMSCorrFctnKK, 1)
#endif
};

#endif
