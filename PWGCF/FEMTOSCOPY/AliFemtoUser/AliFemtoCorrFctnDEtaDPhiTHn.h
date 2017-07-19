////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhiTHn - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Malgorzata Janik majanik@cern.ch                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDETADPHITHN_H
#define ALIFEMTOCORRFCTNDETADPHITHN_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnDEtaDPhiTHn : public AliFemtoCorrFctn {
 public:
  AliFemtoCorrFctnDEtaDPhiTHn(char* title, const int& aPhiBins, const int& aEtaBins, const int &pT1Bins, const double &pT1min, const double &pT1max,  const int &pT2Bins, const double &pT2min, const double &pT2max,  const int &zvtxBins, const double &zvtxmin,  const double &zvtxmax, const int &multBins, const int &multmin, const int &multmax );
  AliFemtoCorrFctnDEtaDPhiTHn(const AliFemtoCorrFctnDEtaDPhiTHn& aCorrFctn);
  virtual ~AliFemtoCorrFctnDEtaDPhiTHn();

  
  AliFemtoCorrFctnDEtaDPhiTHn& operator=(const AliFemtoCorrFctnDEtaDPhiTHn& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);
  virtual void Finish();
  void WriteHistos();
  virtual TList* GetOutputList();

 private:
  THnSparseF *fDPhiDEtaNum;  // Numerator of dEta dPhi function
  THnSparseF *fDPhiDEtaDen;  // Denominator of dEta dPhi function
  double fphiL;
  double fphiT;

  double fPt1Min;
  double fPt1Max;
  double fPt2Min;
  double fPt2Max;
  double fZvtxMin;
  double fZvtxMax;
  double fMultMin;
  double fMultMax;

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDEtaDPhiTHn, 1)
#endif
    };


#endif

