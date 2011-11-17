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
// TRD PID Class
// Does PID either on a x% electron efficiency basis or on dE/dx
// For more information please check the implementation file
//
#ifndef ALIHFEPIDTRD_H
#define ALIHFEPIDTRD_H

 #ifndef ALIHFEPIDBASE_H
 #include "AliHFEpidBase.h"
 #endif

class AliAODTrack;
class AliAODMCParticle;
class AliESDtrack;
class AliHFEcollection;
class AliMCParticle;
class AliOADBContainer;
class AliVParticle;
class AliVTrack;
class TList;
class TH2F;

class AliHFEpidTRD : public AliHFEpidBase{
  public:
    typedef enum{
      kLQ = 0,
      kNN = 1
    } PIDMethodTRD_t;
    enum{
      kThreshParams = 4
    };
    enum{
      kHistTRDlikeBefore = 0,
      kHistTRDlikeAfter = 1,
      kHistTRDthresholds = 2,
      kHistTRDSigV1 = 3,
      kHistTRDSigV2 = 4,
      kHistOverallSpecies = 5
    };
    AliHFEpidTRD();
    AliHFEpidTRD(const Char_t *name);
    AliHFEpidTRD(const AliHFEpidTRD &ref);
    AliHFEpidTRD& operator=(const AliHFEpidTRD &ref);
    virtual ~AliHFEpidTRD();

    virtual Bool_t InitializePID(Int_t run);
    virtual Int_t IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const;

    Double_t GetTRDSignalV1(const AliESDtrack *track, Float_t truncation = 0.7) const;
    Double_t GetTRDSignalV2(const AliESDtrack *track, Float_t trucation = 0.7) const;

    Bool_t IsCalculateTRDSignals() const { return TestBit(kTRDsignals); }
    Bool_t IsRenormalizeElPi() const { return TestBit(kTRDrenormalize); }
    void SelectCutOnTheFly(Bool_t onFly = kTRUE) { if(onFly) SetBit(kSelectCutOnTheFly, kTRUE); else SetBit(kSelectCutOnTheFly, kFALSE);}
    void SetOADBThresholds(AliOADBContainer *cont) { fOADBThresholds = cont; }
    void SetTotalChargeInSlice0() { fTotalChargeInSlice0 = kTRUE; }
    void SetRenormalizeElPi(Bool_t doRenorm = kTRUE) { if(doRenorm) SetBit(kTRDrenormalize, kTRUE); else SetBit(kTRDrenormalize, kFALSE);}
    void SetElectronEfficiency(Double_t electronEfficiency) { fElectronEfficiency = electronEfficiency; }
    void SetNTracklets(Int_t nTracklets) { fNTracklets = nTracklets; }
    void SetMinP(Double_t p) { fMinP = p; }
    void CalculateTRDSignals(Bool_t docalc) { SetBit(kTRDsignals, docalc); } 

    Double_t GetElectronLikelihood(const AliVTrack *track, AliHFEpidObject::AnalysisType_t anaType) const;
    Int_t    GetNTracklets() const { return fNTracklets; }
    void     GetTRDmomenta(const AliVTrack *track, Double_t *mom) const;
    Double_t GetP(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anaType) const;
    Double_t GetTRDthresholds(Double_t p) const;
    Double_t GetTRDthresholds(Double_t p, UInt_t nTracklets) const;
    Double_t GetChargeLayer(const AliVParticle *track, UInt_t layer, AliHFEpidObject::AnalysisType_t anatype) const;
  protected:
    enum{
      kTRDsignals = BIT(16),
      kThresholdsInitialized = BIT(17),
      kTRDrenormalize = BIT(18),
      kSelectCutOnTheFly = BIT(19)
    };
    void Copy(TObject &ref) const;
    Bool_t InitParamsFromOADB(Int_t run);
    void RenormalizeElPi(const Double_t * const likein, Double_t * const likeout) const;

  private:
    AliOADBContainer *fOADBThresholds;                      // OADBContainer with thresholds
    Double_t fMinP;                                         // Minimum momentum above which TRD PID is applied
    Int_t    fNTracklets;                                   // Select cut for the number of tracklets
    Int_t    fRunNumber;                                    // Run number
    Double_t fElectronEfficiency;                           // Cut on electron efficiency
    Double_t fThreshParams[kThreshParams];                  // Threshold parametrisation
    Bool_t fTotalChargeInSlice0;                            // Flag for foreward/backward compatibility for the TRD total charge
  ClassDef(AliHFEpidTRD, 1)     // TRD electron ID class
};
#endif

