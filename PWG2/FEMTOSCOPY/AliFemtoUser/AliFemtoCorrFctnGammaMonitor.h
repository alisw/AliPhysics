////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnGammaMonitor - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNGAMMAMONITOR_H
#define ALIFEMTOCORRFCTNGAMMAMONITOR_H

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnGammaMonitor : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnGammaMonitor(char* title, const int& aMinvBins, const int& aDThetaBins);
  AliFemtoCorrFctnGammaMonitor(const AliFemtoCorrFctnGammaMonitor& aCorrFctn);
  virtual ~AliFemtoCorrFctnGammaMonitor();

  AliFemtoCorrFctnGammaMonitor& operator=(const AliFemtoCorrFctnGammaMonitor& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
private:
  
  TH2D *fNumPMinvDTheta;        // Numerator Minv vs. DTheta Positive kSide
  TH2D *fDenPMinvDTheta;        // Denominator Minv vs. DTheta Positive kSide

  TH2D *fNumNMinvDTheta;        // Numerator Minv vs. DTheta Negative kSide
  TH2D *fDenNMinvDTheta;        // Denominator Minv vs. DTheta Negative kSide

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnGammaMonitor, 1)
#endif
};


#endif

