/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Class for TPC PID
// Does electron selection based on dE/dx
// For more information please check the implementation file
//
#ifndef ALIHFEPIDTPC_H
#define ALIHFEPIDTPC_H

#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

#ifndef ALIPID_H
#include "AliPID.h"
#endif

class TList;
class TF1;
class TH2D;
class AliAODTrack;
class AliAODMCParticle;
class AliESDtrack;
class AliMCParticle;
class AliVParticle;
class AliHFEcollection;
class AliHFEpidQAmanager;

class AliHFEpidTPC : public AliHFEpidBase{
  public:
    AliHFEpidTPC();
    AliHFEpidTPC(const Char_t *name);
    AliHFEpidTPC(const AliHFEpidTPC &ref);
    AliHFEpidTPC &operator=(const AliHFEpidTPC &ref);
    virtual ~AliHFEpidTPC();
    
    virtual Bool_t InitializePID(Int_t /*run*/);
    virtual Int_t IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const;

    void AddTPCdEdxLineCrossing(Int_t species, Double_t sigma);
    Bool_t HasAsymmetricSigmaCut() const { return TestBit(kAsymmetricSigmaCut);}
    Bool_t HasParticleRejection() const { return TestBit(kRejection); }
    void SetTPCnSigma(Short_t nSigma) { fNsigmaTPC = nSigma; };
    void SetUseOnlyOROC(Bool_t useOnlyOROC) { fUseOnlyOROC = useOnlyOROC; };
    inline void SetAsymmetricTPCsigmaCut(Float_t pmin, Float_t pmax, Float_t sigmaMin, Float_t sigmaMax);
    inline void SetRejectParticle(Int_t species, Float_t pmin, Float_t sigmaMin, Float_t pmax, Float_t sigmaMax);

    void SetUpperSigmaCutDefault(const TF1 * const model) { fkUpperSigmaCut[0] = model; fHasCutModel = kTRUE; }
    void SetUpperSigmaCutCentrality(const TF1 * const model, Int_t centralityBin) { if(centralityBin < 11) fkUpperSigmaCut[centralityBin+1] = model; fHasCutModel = kTRUE; }
    void SetLowerSigmaCutDefault(const TF1 * const model) { fkLowerSigmaCut[0] = model; fHasCutModel = kTRUE; }
    void SetLowerSigmaCutCentrality(const TF1 * const model, Int_t centralityBin) { if(centralityBin < 11) fkLowerSigmaCut[centralityBin+1] = model; fHasCutModel = kTRUE; }
    void SetEtaCorrection(const TF1 *const param) { fkEtaCorrection = param; }
    void SetCentralityCorrection(const TF1 *const param){ fkCentralityCorrection = param; }
    void SetEtaCorrections(const TF1 *const mean, const TF1 *const wdth) { fkEtaMeanCorrection = mean; fkEtaWidthCorrection = wdth; }
   void SetMomentumCorrections(const TF1 *const mean, const TF1 *const wdth) { fkPMeanCorrection = mean; fkPWidthCorrection = wdth; }
    void SetCentralityCorrections(const TF1 *const mean, const TF1 *const wdth) { fkCentralityMeanCorrection = mean; fkCentralityWidthCorrection = wdth; }
    void SetJpsiCorrections(const TH2D *const mean, const TH2D *const wdth) { fkCentralityEtaCorrectionMeanJpsi = mean; fkCentralityEtaCorrectionWidthJpsi = wdth; }
    void UsedEdx() { fUsedEdx = kTRUE; }
    void UseNSigma() { fUsedEdx = kFALSE; }
    Bool_t HasEtaCorrection() const { return fkEtaCorrection != NULL; }
    Bool_t HasCentralityCorrection() const { return fkCentralityCorrection != NULL; } 
    Bool_t IsUsingdEdx() const { return fUsedEdx; }

    Double_t GetP(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anaType) const;
    void ApplyEtaCorrection(AliVTrack *track, AliHFEpidObject::AnalysisType_t anatype) const;
    void ApplyCentralityCorrection(AliVTrack *track, Double_t centralityEstimator, AliHFEpidObject::AnalysisType_t anatype) const;
    Double_t GetCorrectedTPCnSigma(Double_t eta, Double_t centralityEstimator, Double_t tpcNsigma, Double_t mom = 0) const;
    Double_t GetCorrectedTPCnSigmaJpsi(Double_t eta, Double_t centralityEstimator, Double_t tpcNsigma) const;
    void UseOROC(AliVTrack *track, AliHFEpidObject::AnalysisType_t anatype) const;

  protected:
    void Copy(TObject &o) const;
    Int_t Reject(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anaType) const;

    Bool_t CutSigmaModel(const AliHFEpidObject *anaType) const;

  private:
    enum{
      kAsymmetricSigmaCut = BIT(20),
      kRejection = BIT(21)
    };
    Double_t fLineCrossingSigma[AliPID::kSPECIES];          // with of the exclusion point
    UChar_t fLineCrossingsEnabled;                          // Bitmap showing which line crossing is set
    const TF1 *fkUpperSigmaCut[12];                         // Upper Sigma Cut
    const TF1 *fkLowerSigmaCut[12];                         // Lower Sigma Cut
    const TF1 *fkEtaCorrection;                             // Correction for the eta dependence
    const TF1 *fkCentralityCorrection;                      // Correction for the centrality dependence
    const TF1 *fkEtaMeanCorrection;                         // Correct eta dependence of the mean of the TPC n sigma
    const TF1 *fkEtaWidthCorrection;                        // Correct eta dependence of the width of the TPC n sigma
   const TF1 *fkPMeanCorrection;                         // Correct momentum dependence of the mean of the TPC n sigma
   const TF1 *fkPWidthCorrection;                        // Correct momentum dependence of the width of the TPC n sigma
    const TF1 *fkCentralityMeanCorrection;                  // Correct centrality dependence of the mean of the TPC n sigma
    const TF1 *fkCentralityWidthCorrection;                 // Correct centrality dependence of the width of the TPC n sigma
    const TH2D *fkCentralityEtaCorrectionMeanJpsi;          // Correction from J/psi group for the mean
    const TH2D *fkCentralityEtaCorrectionWidthJpsi;         // Correction from J/psi group for the width
    Bool_t fHasCutModel;                                    // Has cut model functions
    Bool_t fUseOnlyOROC;                                    // Use only OROC
    Float_t fPAsigCut[2];                                   // Momentum region where to perform asymmetric sigma cut
    Float_t fNAsigmaTPC[2];                                 // Asymmetric TPC Sigma band        
    Short_t fNsigmaTPC;                                     // TPC sigma band
    Float_t fRejection[4*AliPID::kSPECIES];                 // All informations for Particle Rejection, order pmin, sigmin, pmax, sigmax
    UChar_t fRejectionEnabled;                              // Bitmap for enabled particle rejection
    Bool_t  fUsedEdx;                                       // Apply cut on dE/dx instead of number of sigmas

  ClassDef(AliHFEpidTPC, 3)   // TPC Electron ID class
};

inline void AliHFEpidTPC::SetAsymmetricTPCsigmaCut(Float_t pmin, Float_t pmax, Float_t sigmaMin, Float_t sigmaMax) { 
  fPAsigCut[0] = pmin; 
  fPAsigCut[1] = pmax; 
  fNAsigmaTPC[0] = sigmaMin; 
  fNAsigmaTPC[1] = sigmaMax; 
  SetBit(kAsymmetricSigmaCut, kTRUE);
}

inline void AliHFEpidTPC::SetRejectParticle(Int_t species, Float_t pmin, Float_t sigmaMin, Float_t pmax, Float_t sigmaMax){
  if(species < 0 || species >= AliPID::kSPECIES) return;
  fRejection[4*species]   = pmin;
  fRejection[4*species+1] = sigmaMin;
  fRejection[4*species+2] = pmax;
  fRejection[4*species+3] = sigmaMax;
  SETBIT(fRejectionEnabled, species);
  SetBit(kRejection, kTRUE);
}
 
#endif
