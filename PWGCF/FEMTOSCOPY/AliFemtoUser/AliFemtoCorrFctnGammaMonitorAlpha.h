////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnGammaMonitor - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNGAMMAMONITORALPHA_H
#define ALIFEMTOCORRFCTNGAMMAMONITORALPHA_H

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnGammaMonitorAlpha : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnGammaMonitorAlpha(const char* title, const int& aMinvBins, const int& aDAlphaBins);
  AliFemtoCorrFctnGammaMonitorAlpha(const AliFemtoCorrFctnGammaMonitorAlpha& aCorrFctn);
  virtual ~AliFemtoCorrFctnGammaMonitorAlpha();

  AliFemtoCorrFctnGammaMonitorAlpha& operator=(const AliFemtoCorrFctnGammaMonitorAlpha& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnGammaMonitorAlpha(*this); }

private:

  TH2D *fNumPMinvDAlpha;        // Numerator Minv vs. DAlpha Positive kSide
  TH2D *fDenPMinvDAlpha;        // Denominator Minv vs. DAlpha Positive kSide

  TH2D *fNumNMinvDAlpha;        // Numerator Minv vs. DAlpha Negative kSide
  TH2D *fDenNMinvDAlpha;        // Denominator Minv vs. DAlpha Negative kSide

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnGammaMonitorAlpha, 1)
#endif
};


#endif

