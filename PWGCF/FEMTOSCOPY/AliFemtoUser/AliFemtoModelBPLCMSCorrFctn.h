////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelBPLCMSCorrFctn - the class for correlation function which   ///
/// uses the model framework and weight generation and calculated the 3D     ///
/// correlation function in the Bertsh-Pratt LCMS system                     ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELBPLCMSCORRFCTN_H
#define ALIFEMTOMODELBPLCMSCORRFCTN_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3D.h"

class AliFemtoModelBPLCMSCorrFctn : public AliFemtoModelCorrFctn {
 public:
  AliFemtoModelBPLCMSCorrFctn()  :
    AliFemtoModelCorrFctn(),
    fNumerator3DTrue(0),
    fNumerator3DFake(0),
    fDenominator3D(0),
    fQinvHisto(0),
    fPairCut(0),
    fUseRPSelection(0){}
  AliFemtoModelBPLCMSCorrFctn(const char* title, const int& nbins, const float& QLo, const float& QHi);
  AliFemtoModelBPLCMSCorrFctn(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn);
  virtual ~AliFemtoModelBPLCMSCorrFctn();

  AliFemtoModelBPLCMSCorrFctn& operator=(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* pair);
  virtual void AddMixedPair(AliFemtoPair* pair);

  virtual void Finish();

  virtual void Write();
  virtual TList* GetOutputList();

  void SetSpecificPairCut(AliFemtoPairCut* aCut);
  void SetUseRPSelection(unsigned short aRPSel);

  virtual AliFemtoModelCorrFctn* Clone();

protected:
  TH3D* fNumerator3DTrue;            // 3D Numerator with pairs from same event only
  TH3D* fNumerator3DFake;            // 3D Numerator with pairs from mixed events
  TH3D* fDenominator3D;              // 3D Denominator with the weight of 1.0

  TH3D* fQinvHisto;                  // Averag qinv histogram

  AliFemtoPairCut* fPairCut;    //! this is a PairCut specific to THIS CorrFctn, not the Analysis

  unsigned short fUseRPSelection;  // The pair cut uses RP selection
#ifdef __ROOT__
  ClassDef(AliFemtoModelBPLCMSCorrFctn, 1)
#endif
};

#endif

