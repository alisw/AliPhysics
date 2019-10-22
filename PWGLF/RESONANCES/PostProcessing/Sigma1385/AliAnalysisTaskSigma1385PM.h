#ifndef AliAnalysisTaskSigma1385PM_H
#define AliAnalysisTaskSigma1385PM_H

#include <THnSparse.h>
#include <TNtupleD.h>
#include <deque>
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
class THistManager;
class AliPIDResponse;
class AliESDtrackCuts;

class AliAnalysisTaskSigma1385PM : public AliAnalysisTaskSE {
   public:
    AliAnalysisTaskSigma1385PM();
    AliAnalysisTaskSigma1385PM(const char* name, Bool_t MCcase);
    virtual ~AliAnalysisTaskSigma1385PM();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* option);
    virtual void Terminate(Option_t* option);
    void SetMixing(Bool_t setmixing) { fsetmixing = setmixing; }
    void SetnMix(Int_t nMix) { fnMix = nMix; }
    void SetIsPrimaryMC(Bool_t isprimarymc) { IsPrimaryMC = isprimarymc; }
    void SetFillnTuple(Bool_t fillntuple) { fFillnTuple = fillntuple; }

    // Setter for cut variables
    void SetFilterbitSigmaStarPion(Double_t lParameter) {
        fFilterBit = lParameter;
    }
    void SetMaxNsigSigmaStarPion(Double_t lParameter) {
        fTPCNsigSigmaStarPionCut = lParameter;
    }
    void SetMaxEtaSigmaStarPion(Double_t lParameter) {
        fSigmaStarPionEtaCut = lParameter;
    }
    void SetMaxVertexZSigmaStarPion(Double_t lParameter) {
        fSigmaStarPionZVertexCut = lParameter;
    }
    void SetMaxVertexXYsigSigmaStarPion(Double_t lParameter) {
        fSigmaStarPionXYVertexSigmaCut = lParameter;
    }

    void SetMaxNsigV0Proton(Double_t lParameter) {
        fTPCNsigLambdaProtonCut = lParameter;
    }
    void SetMaxNsigV0Pion(Double_t lParameter) {
        fTPCNsigLambdaPionCut = lParameter;
    }
    void SetMaxDCAV0daughters(Double_t lParameter) {
        fDCADistLambdaDaughtersCut = lParameter;
    }
    void SetMaxDCAPVV0(Double_t lParameter) {
        fDCArDistLambdaPVCut = lParameter;
    }
    void SetMinCPAV0(Double_t lParameter) {
        fV0CosineOfPointingAngleCut = lParameter;
    }
    void SetMaxRapidity0(Double_t lParameter) {
        fMaxLambdaRapidity = lParameter;
    }
    void SetLowRadiusV0(Double_t lParameter) {
        fLambdaLowRadius = lParameter;
    }
    void SetHighRadiusV0(Double_t lParameter) {
        fLambdaHighRadius = lParameter;
    }
    void SetLifetimeV0(Double_t lParameter) {
        fLambdaLifetime = lParameter;
    }
    void SetMaxMassWindowV0(Double_t lParameter) {
        fV0MassWindowCut = lParameter;
    }

    void SetSigmaStarRapidityCutHigh(Double_t lParameter) {
        fSigmaStarYCutHigh = lParameter;
    }
    void SetSigmaStarRapidityCutLow(Double_t lParameter) {
        fSigmaStarYCutLow = lParameter;
    }

    Bool_t GoodTracksSelection();
    Bool_t GoodV0Selection();
    void FillTracks();
    void FillNtuples();
    void FillMCinput(AliMCEvent* fMCEvent);
    void FillTrackToEventPool();
    Bool_t IsTrueSigmaStar(UInt_t v0, UInt_t pion);
    double GetTPCnSigma(AliVTrack* track, AliPID::EParticleType type);
    void GetImpactParam(AliVTrack* track, Float_t p[2], Float_t cov[3]);
    void SetCutOpen();

    // helper
    THnSparse* CreateTHnSparse(TString name,
                               TString title,
                               Int_t ndim,
                               std::vector<TAxis> bins,
                               Option_t* opt = "");
    THnSparse* CreateTHnSparse(TString name,
                               TString title,
                               TString templ,
                               Option_t* opt = "");
    Long64_t FillTHnSparse(TString name,
                           std::vector<Double_t> x,
                           Double_t w = 1.);
    Long64_t FillTHnSparse(THnSparse* h,
                           std::vector<Double_t> x,
                           Double_t w = 1.);
    TAxis AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax);
    TAxis AxisVar(TString name, std::vector<Double_t> bin);
    TAxis AxisStr(TString name, std::vector<TString> bin);

    AliEventCuts fEventCuts;  // Event cuts

   private:
    typedef std::vector<AliVTrack*> tracklist;
    typedef std::deque<tracklist> eventpool;
    typedef std::vector<std::vector<eventpool>> mixingpool;

    AliESDtrackCuts* fTrackCuts = nullptr;   //!
    AliPIDResponse* fPIDResponse = nullptr;  //!

    AliVEvent* fEvt = nullptr;        //!
    AliMCEvent* fMCEvent = nullptr;   //!
    THistManager* fHistos = nullptr;  //!
    AliAODVertex* vertex = nullptr;   //!
    Bool_t fsetmixing = kFALSE;
    Bool_t IsMC = kFALSE;
    Bool_t IsPrimaryMC = kFALSE;
    Bool_t fFillnTuple = kFALSE;
    Bool_t IsNano = kFALSE;
    TNtupleD* fNtupleSigma1385;        //! Ntuple for the analysis
    TClonesArray* fMCArray = nullptr;  //!
    mixingpool fEMpool;                //!
    TAxis binCent;                     //!
    TAxis binZ;                        //!
    Double_t lPosPV[3];

    Double_t fCent = -1;
    Int_t fnMix = 10;
    Int_t centbin = -1;
    Int_t zbin = -1;

    // Pion cuts
    UInt_t fFilterBit = 32.0;
    Double_t fTPCNsigSigmaStarPionCut = 3.0;
    Double_t fSigmaStarPionEtaCut = 0.8;
    Double_t fSigmaStarPionZVertexCut = 2.0;  // 2.0
    Double_t fSigmaStarPionXYVertexSigmaCut = 7.0;

    // Lambda cuts
    Double_t fTPCNsigLambdaProtonCut = 3.0;
    Double_t fTPCNsigLambdaPionCut = 3.0;
    Double_t fDCADistLambdaDaughtersCut = 0.5;    // 0.5
    Double_t fDCArDistLambdaPVCut = 0.3;          // 0.3
    Double_t fV0CosineOfPointingAngleCut = 0.99;  // 0.99
    Double_t fMaxLambdaRapidity = 0.5;
    Double_t fLambdaLowRadius = 0.5;
    Double_t fLambdaHighRadius = 200.0;
    Double_t fLambdaLifetime = 30.0;
    Double_t fV0MassWindowCut = 0.01;

    // Sigma Star cut
    Double_t fSigmaStarYCutHigh = 0.5;
    Double_t fSigmaStarYCutLow = -0.5;

    std::vector<UInt_t> goodtrackindices;  //!
    std::vector<std::vector<UInt_t>> goodv0indices;  //!

    ClassDef(AliAnalysisTaskSigma1385PM, 4);
    // Add rapidity/radius/Lifetime/Y cut of lambda
    // Add NanoOption
    // 4: Add GetImpactParm function for nano
};

#endif