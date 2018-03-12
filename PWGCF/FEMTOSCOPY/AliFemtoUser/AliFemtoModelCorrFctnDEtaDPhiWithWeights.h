////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhi - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
// it uses histograms to provide basic acceptance filter                      //
// Authors: Rafal Maselek rmaselek@cern.ch                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMODELCORRFCTNDETADPHIWITHWEIGHTS_H
#define ALIFEMTOMODELCORRFCTNDETADPHIWITHWEIGHTS_H

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelCorrFctn.h"

class AliFemtoModelCorrFctnDEtaDPhiWithWeights : public AliFemtoModelCorrFctn {
public:
  AliFemtoModelCorrFctnDEtaDPhiWithWeights(const char* title,  TH2D* filter1,  TH2D* filter2, const int& aPhiBins, const int& aEtaBins);
  AliFemtoModelCorrFctnDEtaDPhiWithWeights(const AliFemtoModelCorrFctnDEtaDPhiWithWeights& aCorrFctn);
  virtual ~AliFemtoModelCorrFctnDEtaDPhiWithWeights();

  AliFemtoModelCorrFctnDEtaDPhiWithWeights& operator=(const AliFemtoModelCorrFctnDEtaDPhiWithWeights& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
private:
  TH2D *filterHist1; //filter hisotgram used for the first particle
  TH2D *filterHist2; //filter hisotgram used for the second particle

  TH2D *fDPhiDEtaNumeratorTrue;      // Numerator of dEta dPhi true function
  TH2D *fDPhiDEtaNumeratorFake;      // Numerator of dEta dPhi fake function
  TH2D *fDPhiDEtaDenominator;        // Denominator of dEta dPhi function

  TH2D *fDPhiDEtaColNumerator;       // Numerator of colinear dEta dPhi function
  TH2D *fDPhiDEtaColDenominator;     // Denominator of colinear dEta dPhi function

  TH1D *fDPhiNumeratorTrue;          // Numerator of dPhi true correlation
  TH1D *fDPhiNumeratorFake;          // Numerator of dPhi fake correlation
  TH1D *fDPhiDenominator;            // Denominator of dPhi correlation

  TH1D *fDCosNumeratorTrue;           // Numerator of colinearity true correlation
  TH1D *fDCosNumeratorFake;           // Numerator of colinearity fake correlation
  TH1D *fDCosDenominator;            // Denominator of colinearity correlation

  TH2D *fDPhiPtNumerator;            // Numerator of dPhi correlation vs. Pt min
  TH2D *fDPhiPtDenominator;          // Denominator of dPhi correlation vs. Pt min

  TH2D *fDCosPtNumerator;            // Numerator of colinearity correlation vs. Pt min
  TH2D *fDCosPtDenominator;          // Denominator of colinearity correlation vs. Pt min

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnDEtaDPhiWithWeights, 1)
#endif
};


#endif

