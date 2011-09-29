#ifndef AliAnalysisTaskPi0Reconstruction_cxx
#define AliAnalysisTaskPi0Reconstruction_cxx

#include "TH1.h"
#include "AliESDEvent.h"

#include "AliAODConversionMother.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionPhoton.h"
#include "AliConversionAODBGHandlerRP.h"

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliKFParticle.h"
#include "AliAnalysisTaskSE.h"
#include "TRandom3.h"
#include "TH1.h"
#include "AliLog.h"
#include "AliV0ReaderV1.h"

class AliAnalysisTaskPi0Reconstruction : public AliV0ReaderV1{
public:
    AliAnalysisTaskPi0Reconstruction(const char *name);
    virtual ~AliAnalysisTaskPi0Reconstruction();
  
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    // public getter functions

    TClonesArray *GetPi0Candidates(){return fPi0Candidates;}
    TClonesArray *GetBGPi0s(){return fBGPi0s;}
    TClonesArray *GetTrueMCPi0s(){return fMCTruePi0s;}


    // public Set Functions
    void SetBGHandler(AliConversionAODBGHandlerRP *bgHandler);
    void SetDefaultBGHandler();
    void SetDeltaAODBranchName(TString string) { fDeltaAODBranchName = string;AliInfo(Form("Set DeltaAOD BranchName to: %s",fDeltaAODBranchName.Data()));}
    void SetPi0MassRange(Double_t low,Double_t up){Pi0MassRange[0]=low;Pi0MassRange[1]=up;};
    void SetEventMixing(Bool_t k){kEventMixing=k;};
    void SetIsHeavyIon(Bool_t k){fIsHeavyIon=k;};
    void SetUseSatelliteAODs(Bool_t k){kUseSatelliteAODs=k;};
    void SetInvMassRange(Double_t range[2]){Pi0MassRange[0]=range[0];Pi0MassRange[1]=range[1];AliInfo(Form("Set Invariant Mass Range: %f - %f",Pi0MassRange[0],Pi0MassRange[1]));};
    void SetRapidityCut(Double_t maxrapidity){fRapidityCut=maxrapidity;AliInfo(Form(" Rapidity Range for Pi0: %f - %f",-fRapidityCut,fRapidityCut));};
    void SetAlphaCut(Double_t maxalpha){fAlphaCut=maxalpha;AliInfo(Form("Alpha Cut: %f",fAlphaCut));};
    void SetUseOnlyTaggedPhotons(Bool_t ktag){kUseOnlyTaggedPhotons=ktag;}
    void SetNRotations(Int_t nrot){fNRandomEventsForBGRotation=nrot;}

protected:

    Bool_t IsMCMesonInReconstructionAcceptance(TParticle *fMCMother);

    void FindDeltaAODBranchName();
    void InitializeBGHandler();

    void CalculatePi0Candidates();
    void CalculateBackground();

    void RotateParticle(AliAODConversionPhoton *gamma=0x0,Double_t angle=0);

    Bool_t GetConversionGammas();
    Bool_t GetPHOSGammas();
    Bool_t GetEMCALGammas();

    AliVTrack *GetTrack(Int_t label=-1);
    Bool_t IsTruePi0(AliAODConversionMother *);
    Bool_t CheckGamma(AliAODConversionPhoton *gamma=0x0);
    Bool_t CheckPi0(AliAODConversionMother *pi0=0x0,Bool_t IsSignal=kTRUE);

    void ProcessMCMesons();

    Bool_t IsMCPhotonReconstructed(TParticle *MDaughter);

    TString     fDeltaAODBranchName;//! File where Gamma Conv AOD is located, if not in default AOD
    AliConversionAODBGHandlerRP *fBGHandler;

    Double_t *Pi0MassRange;
    Bool_t kEventMixing;
    Bool_t fIsHeavyIon;
    TRandom3 *fRandomizer;
    Bool_t kUseSatelliteAODs;
    TClonesArray *fPHOSGammas;
    TClonesArray *fEMCALGammas;
    TClonesArray *fPi0Candidates;
    TClonesArray *fBGPi0s;
    TClonesArray *fMCTruePi0s;
    Double_t fRapidityCut;
    Double_t fAlphaCut;
    Bool_t kUseOnlyTaggedPhotons;
    Int_t fNRandomEventsForBGRotation;

    // Histograms

    TH1F *hPi0Cuts;
    TH1F *hPi0BGCuts;
    TH1F *hPi0Alpha;
    TH1F *hPi0AlphaRejected;
    TH1F *hPi0OpeningAngle;
    TH1F *hPi0OpeningAngleRejected;
    TH1F *hPi0Rapidity;
    TH1F *hPi0RapidityRejected;
    TH1F **hPi0TRUE;
    TH2F **hPi0RECOTRUE;

    TH3F *hPool;

    AliAnalysisTaskPi0Reconstruction(const AliAnalysisTaskPi0Reconstruction&); // not implemented
    AliAnalysisTaskPi0Reconstruction& operator=(const AliAnalysisTaskPi0Reconstruction&); // not implemented
  
    ClassDef(AliAnalysisTaskPi0Reconstruction, 1); // example of analysis
};

#endif

