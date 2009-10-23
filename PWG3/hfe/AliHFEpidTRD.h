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
#ifndef ALIHFEPIDTRD_H
#define ALIHFEPIDTRD_H

 #ifndef ALIHFEPIDBASE_H
 #include "AliHFEpidBase.h"
 #endif

class AliAODTrack;
class AliAODMCParticle;
class AliESDtrack;
class AliMCParticle;
class AliVParticle;
class TList;
class TH2F;

class AliHFEpidTRD : public AliHFEpidBase{
  public:
    typedef enum{
      kLQ = 0,
      kNN = 1
    } PIDMethodTRD_t;
    enum{
      kThreshParams = 24
    };
    enum{
      kHistTRDlikeBefore = 0,
      kHistTRDlikeAfter = 1,
      kHistTRDthresholds = 2,
      kHistTRDSigV1 = 3,
      kHistTRDSigV2 = 4,
      kHistOverallSpecies = 5
    };
    AliHFEpidTRD(const Char_t *name);
    AliHFEpidTRD(const AliHFEpidTRD &ref);
    AliHFEpidTRD& operator=(const AliHFEpidTRD &ref);
    virtual ~AliHFEpidTRD();
    
    virtual Bool_t InitializePID();
    virtual Int_t IsSelected(AliHFEpidObject *track);
    virtual Bool_t HasQAhistos() const { return kTRUE; };

    Double_t GetTRDSignalV1(AliESDtrack *track, Int_t mcPID);
    Double_t GetTRDSignalV2(AliESDtrack *track, Int_t mcPID);

    void SetPIDMethod(PIDMethodTRD_t method) { fPIDMethod = method; };
  protected:
    void Copy(TObject &ref) const;
    Int_t MakePIDesd(AliESDtrack *esdTrack, AliMCParticle *mcTrack);
    Int_t MakePIDaod(AliAODTrack *aofTrack, AliAODMCParticle *aodMC);
    Double_t GetTRDthresholds(Double_t electronEff, Double_t p);
    Int_t GetMCpid(AliESDtrack *track);
    void InitParameters();
    virtual void AddQAhistograms(TList *l);
    void GetParameters(Double_t electronEff, Double_t *parameters);

    void FillHistogramsLikelihood(Int_t whenFilled, Float_t p, Float_t elProb);
    void FillHistogramsTRDSignalV1(Double_t signal, Double_t p, Int_t species);
    void FillHistogramsTRDSignalV2(Double_t signal, Double_t p, Int_t species);
  private:
    static const Double_t fgkVerySmall;                       // Check for 0
    PIDMethodTRD_t fPIDMethod;                              // PID Method: 2D Likelihood or Neural Network
    Double_t fThreshParams[kThreshParams];                  // Threshold parametrisation
    TList *fContainer;                                      // QA  Histogram Container
  ClassDef(AliHFEpidTRD, 1)     // TRD electron ID class
};

#endif
