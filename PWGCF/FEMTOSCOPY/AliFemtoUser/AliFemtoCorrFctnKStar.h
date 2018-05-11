///
/// \file AliFemtoCorrFctnKStar.h
///

#pragma once

#ifndef ALIFEMTOCORRFCTNKSTAR_H
#define ALIFEMTOCORRFCTNKSTAR_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TString.h"

#include "AliFemtoCorrFctn.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"


/// \class AliFemtoCorrFctnKStar
/// \brief A simple KStar correlation function
///
/// Based off AliFemtoQinvCorrFctn and Kubera's AliFemtoCorrFctnKStar
///
///  \authors  Andrew Kubera, Ohio State University, <andrew.kubera@cern.ch>
///            Jesse Buxton, Ohio State University, <jesse.thomas.buxton@cern.ch>
///
class AliFemtoCorrFctnKStar : public AliFemtoCorrFctn {
public:

  /// Default Constructor
  ///
  /// Title is set to "CorrFctnKStar", consult source file for other default values.
  ///
  AliFemtoCorrFctnKStar();

  /// Constructor
  ///
  /// If non-default values are desired for other binning,
  /// use methods SetKStarVskTBins, SetKStarVsmTBins, etc.
  AliFemtoCorrFctnKStar(const char* title,
                        const int nbins,
                        const float KStarLo,
                        const float KStarHi);

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

  void RotateThreeVecBy180InTransversePlane(AliFemtoThreeVector &a3Vec);
  bool PassPairCut_RotatePar2(const AliFemtoPair* aPair);
  double CalcKStar_RotatePar2(const AliFemtoPair* aPair);    //Rotate the second particle in the pair by 180 degrees
                                                             // about the z-axis

  float CalcMt(const AliFemtoPair* aPair);
  float CalcMtv2(const AliFemtoPair* aPair);  //TODO testing effect of m_reduced vs 0.5(m1+m2) in calculation


  void SetCalculateDetaDphis(Bool_t, Double_t);
  void SetCalculatePairKinematics(Bool_t);

  void SetBuildkTBinned(Bool_t aBuild);
  void SetBuildmTBinned(Bool_t aBuild);
  void SetBuild3d(Bool_t aBuild);

  //inline functions
  TH1D* Numerator();
  TH1D* Denominator();
  TH1D* Ratio();

  virtual AliFemtoCorrFctn* Clone() const  { return new AliFemtoCorrFctnKStar(*this); }

protected:
  TString fTitle;
  int fNbinsKStar;
  double fKStarLow, fKStarHigh;
  TH1D* fNumerator;          // numerator - real pairs
  TH1D* fDenominator;        // denominator - mixed pairs
  TH1D* fRatio;              // unnormalized ratio num/den
  TH1D* fkTMonitor;          // Monitor the kT of pairs in the function

  TH1D* fNumerator_RotatePar2;

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
  //TODO testing effect of m_reduced vs 0.5(m1+m2) in calculation
  TH2F *fNumeratorv2_mT;
  TH2F *fDenominatorv2_mT;

  // 3D k*_out, _side, _long Cf
  Bool_t fBuild3d;
  TH3F *fNumerator3d;
  TH3F *fDenominator3d;
};

inline TH1D* AliFemtoCorrFctnKStar::Numerator(){return fNumerator;}
inline TH1D* AliFemtoCorrFctnKStar::Denominator(){return fDenominator;}
inline TH1D* AliFemtoCorrFctnKStar::Ratio(){return fRatio;}


#endif
