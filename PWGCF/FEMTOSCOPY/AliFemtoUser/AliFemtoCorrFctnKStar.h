// file AliFemtoCorrFctnKStar.h

/*
 *  class AliFemtoCorrFctnKStar.h
 *  Based off AliFemtoQinvCorrFctn and Kubera's AliFemtoCorrFctnKStar
 *  brief A simple KStar correlation function
 *
 *  authors:  Andrew Kubera, Ohio State University, <andrew.kubera@cern.ch>
 *            Jesse Buxton, Ohio State University, <jesse.thomas.buxton@cern.ch>
 *
 */

#ifndef ALIFEMTOCORRFCTNKSTAR_H
#define ALIFEMTOCORRFCTNKSTAR_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"


#include "AliFemtoCorrFctn.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoCorrFctnKStar : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnKStar();
  AliFemtoCorrFctnKStar(const char* title, const int& nbins, const float& KStarLo, const float& KStarHi); //If non-default values are desired for other binning,
                                                                                                          //use methods SetKStarVskTBins, SetKStarVsmTBins, etc.
  AliFemtoCorrFctnKStar(const AliFemtoCorrFctnKStar& aCorrFctn);
  virtual ~AliFemtoCorrFctnKStar();

  AliFemtoCorrFctnKStar& operator=(const AliFemtoCorrFctnKStar& aCorrFctn);

  virtual AliFemtoString Report();
  virtual TList* GetOutputList();
  virtual void Finish();
  void Write();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  void FillDEtaDPhiS(TH2D* aHist, AliFemtoPair* aPair);

  //TODO check these
  void SetkTMonitorBins(int aNbinskT, double akTMin, double akTMax);
  void SetKStarVskTBins(int aNbinsKStar, double aKStarMin, double aKStarMax,
                        int aNbinskT,    double akTMin,    double akTMax);
  void SetKStarVsmTBins(int aNbinsKStar, double aKStarMin, double aKStarMax,
                        int aNbinsmT,    double amTMin,    double amTMax);

  void SetDEtaDPhiSBins(int aNbinsDEta,  double aDEtaMin,  double aDEtaMax,
                        int aNbinsDPhis, double aDPhiSMin, double aDPhiSMax);

  void Set3dBins(int aNbinsKStarOut,  double aKStarOutMin,  double aKStarOutMax,
                 int aNbinsKStarSide, double aKStarSideMin, double aKStarSideMax,
                 int aNbinsKStarLong, double aKStarLongMin, double aKStarLongMax);

  float CalcMt(const AliFemtoPair* aPair);

  //inline functions
  void SetCalculateDetaDphis(Bool_t, Double_t);
  void SetCalculatePairKinematics(Bool_t);

  void SetBuildkTBinned(Bool_t aBuild);
  void SetBuildmTBinned(Bool_t aBuild);
  void SetBuild3d(Bool_t aBuild);

  TH1D* Numerator();
  TH1D* Denominator();
  TH1D* Ratio();

protected:
  TH1D* fNumerator;          // numerator - real pairs
  TH1D* fDenominator;        // denominator - mixed pairs
  TH1D* fRatio;              // unnormalized ratio num/den
  TH1D* fkTMonitor;          // Monitor the kT of pairs in the function

  Bool_t fDetaDphiscal;
  Bool_t fPairKinematics;

  Double_t fRaddedps;
  TH2D* fNumDEtaDPhiS;
  TH2D* fDenDEtaDPhiS;

  TNtuple* fPairKStar; //PairReader for CorrFit

  // kT binned k* Cf
  Bool_t fBuildkTBinned;
  TH2F *fNumerator_kT;
  TH2F *fDenominator_kT;

  // mT binned k* Cf
  Bool_t fBuildmTBinned;
  TH2F *fNumerator_mT;
  TH2F *fDenominator_mT;

  // 3D k*_out, _side, _long Cf
  Bool_t fBuild3d;
  TH3F *fNumerator3d;
  TH3F *fDenominator3d;



#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnKStar, 1)
#endif
};

inline void AliFemtoCorrFctnKStar::SetCalculateDetaDphis(Bool_t dedpsc, Double_t rad) {fDetaDphiscal = dedpsc; fRaddedps = rad;}
inline void AliFemtoCorrFctnKStar::SetCalculatePairKinematics(Bool_t pk) {fPairKinematics = pk;}

inline void AliFemtoCorrFctnKStar::SetBuildkTBinned(Bool_t aBuild) {fBuildkTBinned = aBuild;}
inline void AliFemtoCorrFctnKStar::SetBuildmTBinned(Bool_t aBuild) {fBuildmTBinned = aBuild;}
inline void AliFemtoCorrFctnKStar::SetBuild3d(Bool_t aBuild) {fBuild3d = aBuild;}

inline TH1D* AliFemtoCorrFctnKStar::Numerator(){return fNumerator;}
inline TH1D* AliFemtoCorrFctnKStar::Denominator(){return fDenominator;}
inline TH1D* AliFemtoCorrFctnKStar::Ratio(){return fRatio;}


#endif
