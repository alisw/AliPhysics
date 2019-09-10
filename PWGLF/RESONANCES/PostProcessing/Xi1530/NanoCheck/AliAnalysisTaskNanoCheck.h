#ifndef AliAnalysisTaskNanoCheck_H
#define AliAnalysisTaskNanoCheck_H

#include <THnSparse.h>
#include <TNtupleD.h>
#include <deque>
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
class THistManager;
class AliPIDResponse;
class AliESDtrackCuts;

class AliAnalysisTaskNanoCheck : public AliAnalysisTaskSE {
   public:
    AliAnalysisTaskNanoCheck();
    AliAnalysisTaskNanoCheck(const char* name, Bool_t MCcase);
    virtual ~AliAnalysisTaskNanoCheck();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* option);
    virtual void Terminate(Option_t* option);

    // Setter for cut variables

    Bool_t GoodTracksSelection();
    Bool_t GoodV0Selection();
    void FillTracks();
    double GetTPCnSigma(AliVTrack* track, AliPID::EParticleType type);

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
    AliESDtrackCuts* fTrackCuts = nullptr;   //!
    AliPIDResponse* fPIDResponse = nullptr;  //!

    AliVEvent* fEvt = nullptr;        //!
    AliMCEvent* fMCEvent = nullptr;   //!
    THistManager* fHistos = nullptr;  //!
    AliAODVertex* vertex = nullptr;   //!
    Bool_t IsMC = kFALSE;
    TClonesArray* fMCArray = nullptr;  //!
    TAxis binCent;                     //!
    TAxis binZ;                        //!
    Double_t lPosPV[3];
    Double_t fZ = 999;

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

    ClassDef(AliAnalysisTaskNanoCheck, 2);
    // Add rapidity/radius/Lifetime/Y cut of lambda
};

#endif