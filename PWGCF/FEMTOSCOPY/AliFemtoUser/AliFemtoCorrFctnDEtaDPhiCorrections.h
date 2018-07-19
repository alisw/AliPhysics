////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhiCorrections - A correlation function that analyzes //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDETADPHICORRECTIONS_H
#define ALIFEMTOCORRFCTNDETADPHICORRECTIONS_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "THnSparse.h"
#include "TFile.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnDEtaDPhiCorrections : public AliFemtoCorrFctn {
public:
  enum CorrectionType {kNone=0, kPt=1, kEta=2};
  enum ParticleType {kNoCorrection=0, kPion=1, kKaon=2, kProton=3, kAll=4, kPionMinus=5, kKaonMinus=6, kProtonMinus=7};
  typedef enum CorrectionType ReadCorrectionType;

  AliFemtoCorrFctnDEtaDPhiCorrections(const char* title, const int& aPhiBins, const int& aEtaBins);
  AliFemtoCorrFctnDEtaDPhiCorrections(const AliFemtoCorrFctnDEtaDPhiCorrections& aCorrFctn);
  virtual ~AliFemtoCorrFctnDEtaDPhiCorrections();

  AliFemtoCorrFctnDEtaDPhiCorrections& operator=(const AliFemtoCorrFctnDEtaDPhiCorrections& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();
  void SetDoFullAnalysis(Bool_t do2d);
  void SetCalculatePairPurity(Bool_t dopp);
  double CalculateCorrectionWeight(double pT1, double pT2);
  double CalculateCorrectionWeight(double pT1);
  double CalculateCorrectionWeight(double pT1, double pT2, double eta1, double eta2, double phi1, double phi2, double zvert1, double zvert2);
  void LoadCorrectionTabFromROOTFile1D(const char *file, ParticleType partType1, ParticleType partType2);
  void LoadCorrectionTabFromROOTFile(const char *file, ParticleType partType1, ParticleType partType2, bool doPtCorr, bool doEtaCorr, bool doPhiCorr, bool doZVertCorr);
  //  void LoadCorrectionTabFromFile(const char *pTtab, const char *corrTab); // Not implemented
  void SetCorrectionTab(ParticleType partType);
  double GetPurity(double pT1, int n=1); //n == 1 - first particle, n == 2 - second particle

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnDEtaDPhiCorrections(*this); }

private:

  TH2D *fDPhiDEtaNumerator;          // Numerator of dEta dPhi function
  TH2D *fDPhiDEtaDenominator;        // Denominator of dEta dPhi function

  TH1D *fDPhiNumerator;              // Numerator of dPhi correlation
  TH1D *fDPhiDenominator;            // Denominator of dPhi correlation

  TH1D *fDCosNumerator;              // Numerator of colinearity correlation
  TH1D *fDCosDenominator;            // Denominator of colinearity correlation

  Bool_t   fDoFullAnalysis;               // set to 1 to do 2D Pt analysis
  Bool_t   fCalculatePairPurity;          // set to 1 to calculate pair purity

  TH1D *fPhi;
  TH1D *fEta;
  TH1D *fPtSumDist;

  TH2D *fYtYtNumerator;
  TH2D *fYtYtDenominator;

  TH2F *fPairPurity;
  TH2F *fDPhiDEtaNumeratorNoCorr;

  bool fIfCorrection;
  THnSparseF *fPtCorrectionsNum;
  THnSparseF *fPtCorrectionsDen;

  THnSparseF *fEtaCorrectionsNum;
  THnSparseF *fEtaCorrectionsDen;

  double* fCorrFactorTab;
  double* fpTab;
  ParticleType fPartType; // particle type for calculations of correction factor

  double fphiL;
  double fphiT;

  TFile *ifileCorrTab;
  bool fdoPtCorr;
  bool fdoEtaCorr;
  bool fdoPhiCorr;
  bool fdoZVertCorr;
  int fpartType1;
  int fpartType2;

  THnT<float>* fhntReco1;
  THnT<float>* fhntReco2;
  TH1F *fh1Reco1;
  TH1F *fh1Reco2;
  TH2F *fh2Reco1;
  TH2F *fh2Reco2;
  TH3F *fh3Reco1;
  TH3F *fh3Reco2;
  TH1D *fhCont1;
  TH1D *fhCont2;

  TH1F *fSinglePurity1;
  TH1F *fSinglePurity2;

  bool fCorr1D;

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDEtaDPhiCorrections, 1)
#endif
};


#endif

