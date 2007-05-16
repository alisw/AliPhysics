////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelBPLCMSCorrFctn - the class for correlation function which   ///
/// uses the model framework and weight generation and calculated the 3D     ///
/// correlation function in the Bertsh-Pratt LCMS system                     ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef AliFemtoModelBPLCMS3DCorrFctn_hh
#define AliFemtoModelBPLCMS3DCorrFctn_hh

#include "AliFemtoCorrFctn.h"
#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3D.h"

class AliFemtoModelBPLCMS3DCorrFctn : public AliFemtoModelCorrFctn {
 public:
  AliFemtoModelBPLCMSCorrFctn();
  AliFemtoModelBPLCMSCorrFctn(char* title, const int& nbins, const float& QLo, const float& QHi);
  AliFemtoModelBPLCMSCorrFctn(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn);
  virtual ~AliFemtoModelBPLCMSCorrFctn();

  AliFemtoModelBPLCMSCorrFctn& operator=(const AliFemtoModelBPLCMSCorrFctn& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair*);
  virtual void AddMixedPair(AliFemtoPair*);

  virtual void Finish();

  virtual void Write();

  virtual AliFemtoModelCorrFctnSource* Clone();

  TH3D* Numerator();
  TH3D* Denominator();
  TH3D* QinvHisto();

private:
  TH3D* fNumerator;
  TH3D* fDenominator;

  TH3D* fQinvHisto;

#ifdef __ROOT__
  ClassDef(AliFemtoModelBPLCMSCorrFctn, 1)
#endif
};

inline  TH3D* AliFemtoModelBPLCMSCorrFctn::Numerator(){return fNumerator;}
inline  TH3D* AliFemtoModelBPLCMSCorrFctn::Denominator(){return fDenominator;}
inline  TH3D* AliFemtoModelBPLCMSCorrFctn::QinvHisto(){return fQinvHisto;}

#endif

