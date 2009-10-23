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
class AliTPCpidESD;
class AliVParticle;

class AliHFEpidTPC : public AliHFEpidBase{
  public:
    AliHFEpidTPC(const Char_t *name);
    AliHFEpidTPC(const AliHFEpidTPC &ref);
    AliHFEpidTPC &operator=(const AliHFEpidTPC &ref);
    virtual ~AliHFEpidTPC();
    
    virtual Bool_t InitializePID();
    virtual Int_t IsSelected(AliHFEpidObject *track);
    virtual Bool_t HasQAhistos() const { return kTRUE; };

    void AddTPCdEdxLineCrossing(Int_t species, Double_t sigma);
    Bool_t HasAsymmetricSigmaCut() const { return TestBit(kAsymmetricSigmaCut);}
    Bool_t HasParticleRejection() const { return TestBit(kRejection); }
    void SetTPCnSigma(Short_t nSigma) { fNsigmaTPC = nSigma; };
    inline void SetAsymmetricTPCsigmaCut(Float_t pmin, Float_t pmax, Float_t sigmaMin, Float_t sigmaMax);
    inline void SetRejectParticle(Int_t species, Float_t pmin, Float_t sigmaMin, Float_t pmax, Float_t sigmaMax);

  protected:
    void Copy(TObject &o) const;
    void AddQAhistograms(TList *qaList);
    void FillTPChistograms(const AliESDtrack *track, const AliMCParticle *mctrack);
    Int_t MakePIDaod(AliAODTrack *aodTrack, AliAODMCParticle *mcTrack);
    Int_t MakePIDesd(AliESDtrack *esdTrack, AliMCParticle *mcTrack);
    Int_t Reject(AliESDtrack *track);
    Double_t Likelihood(const AliESDtrack *track, Int_t species, Float_t rsig = 2.);
    Double_t Suppression(const AliESDtrack *track, Int_t species);

  private:
    typedef enum{
      kHistTPCelectron = 0,
      kHistTPCpion = 1,
      kHistTPCmuon = 2,
      kHistTPCkaon = 3,
      kHistTPCproton = 4,
      kHistTPCothers = 5,
      kHistTPCall = 6,
      kHistTPCselected = 7,
      kHistTPCprobEl = 8,
      kHistTPCprobPi = 9,
      kHistTPCprobMu = 10,
      kHistTPCprobKa = 11,
      kHistTPCprobPro = 12,
      kHistTPCprobOth = 13,
      kHistTPCprobAll = 14,
      kHistTPCsuppressPi = 15,
      kHistTPCsuppressMu = 16,
      kHistTPCsuppressKa = 17,
      kHistTPCsuppressPro = 18,
      kHistTPCenhanceElPi = 19,
      kHistTPCenhanceElMu = 20,
      kHistTPCenhanceElKa = 21,
      kHistTPCenhanceElPro = 22,
      kHistTPCElprobPi = 23,
      kHistTPCElprobMu = 24,
      kHistTPCElprobKa = 25,
      kHistTPCElprobPro = 26
    } QAHist_t;
    enum{
      kAsymmetricSigmaCut = BIT(20),
      kRejection = BIT(21)
    };
    Double_t fLineCrossingSigma[AliPID::kSPECIES];          // with of the exclusion point
    UChar_t fLineCrossingsEnabled;                          // Bitmap showing which line crossing is set
    Float_t fPAsigCut[2];                                   // Momentum region where to perform asymmetric sigma cut
    Float_t fNAsigmaTPC[2];                                 // Asymmetric TPC Sigma band        
    Short_t fNsigmaTPC;                                     // TPC sigma band
    Float_t fRejection[4*AliPID::kSPECIES];                 // All informations for Particle Rejection, order pmin, sigmin, pmax, sigmax
    UChar_t fRejectionEnabled;                              // Bitmap for enabled particle rejection
    AliPID *fPID;                                           //! PID Object
    AliTPCpidESD *fPIDtpcESD;                               //! TPC PID object
    TList *fQAList;                                         //! QA histograms

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
