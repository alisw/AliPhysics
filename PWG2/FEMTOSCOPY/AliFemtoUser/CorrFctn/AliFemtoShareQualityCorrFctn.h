////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoShareQualityCorrFctn - A correlation function that saves the     ///
/// amount of sharing and splitting hits per pair as a function of qinv      ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoShareQualityCorrFctn_hh
#define AliFemtoShareQualityCorrFctn_hh

#include "TH1D.h"
#include "TH2D.h"
#include "Base/AliFemtoCorrFctn.h"

class AliFemtoShareQualityCorrFctn : public AliFemtoCorrFctn {
public:
  AliFemtoShareQualityCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi);
  AliFemtoShareQualityCorrFctn(const AliFemtoShareQualityCorrFctn& aCorrFctn);
  virtual ~AliFemtoShareQualityCorrFctn();

  AliFemtoShareQualityCorrFctn& operator=(const AliFemtoShareQualityCorrFctn& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(const AliFemtoPair*);
  virtual void AddMixedPair(const AliFemtoPair*);

  virtual void Finish();

  void WriteHistos();
private:
  
  TH2D *fShareNumerator;
  TH2D *fShareDenominator;

  TH2D *fQualityNumerator;
  TH2D *fQualityDenominator;

#ifdef __ROOT__
  ClassDef(AliFemtoShareQualityCorrFctn, 1)
#endif
};


#endif

