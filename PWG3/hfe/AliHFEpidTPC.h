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
    inline void SetAsymmetricTPCsigmaCut(Float_t pmin, Float_t pmax, Float_t sigmaMin, Float_t sigmaMax);
    inline void SetRejectParticle(Int_t species, Float_t pmin, Float_t sigmaMin, Float_t pmax, Float_t sigmaMax);

    void SetUpperSigmaCutDefault(const TF1 * const model) { fkUpperSigmaCut[0] = model; fHasCutModel = kTRUE; }
    void SetUpperSigmaCutCentrality(const TF1 * const model, Int_t centralityBin) { if(centralityBin < 11) fkUpperSigmaCut[centralityBin+1] = model; fHasCutModel = kTRUE; }
    void SetLowerSigmaCutDefault(const TF1 * const model) { fkLowerSigmaCut[0] = model; fHasCutModel = kTRUE; }
    void SetLowerSigmaCutCentrality(const TF1 * const model, Int_t centralityBin) { if(centralityBin < 11) fkLowerSigmaCut[centralityBin+1] = model; fHasCutModel = kTRUE; }

    Double_t GetP(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anaType) const;

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
    Bool_t fHasCutModel;                                    // Has cut model functions
    Float_t fPAsigCut[2];                                   // Momentum region where to perform asymmetric sigma cut
    Float_t fNAsigmaTPC[2];                                 // Asymmetric TPC Sigma band        
    Short_t fNsigmaTPC;                                     // TPC sigma band
    Float_t fRejection[4*AliPID::kSPECIES];                 // All informations for Particle Rejection, order pmin, sigmin, pmax, sigmax
    UChar_t fRejectionEnabled;                              // Bitmap for enabled particle rejection

  ClassDef(AliHFEpidTPC, 1)   // TPC Electron ID class
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
