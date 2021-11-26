#ifndef AliAnalysisTaskXi1530_H
#define AliAnalysisTaskXi1530_H
//
// Class AliAnalysisTaskXi1530
//
// AliAnalysisTaskXi1530
//  author: Bong-Hwi Lim (bong-hwi.lim@cern.ch)
//        , Beomkyu  KIM (kimb@cern.ch)
//

#include <AliAnalysisTaskSE.h>
#include <THnSparse.h>
#include <deque>
#include "AliEventCuts.h"

class AliAnalysisTask;
class AliESDtrackCuts;
class AliESDEvent;
class AliMCEvent;
class AliAODEvent;
class AliPIDResponse;
class AliPIDCombined;
class THistManager;

class AliAnalysisTaskXi1530 : public AliAnalysisTaskSE {
   public:
    enum {
        kSD = 0,
        kDD,
        kND,
        kCD,
        kAllProc,
        kXiStarCode = 3324,  // Xi(1530)^0 MC code
        kXiCode = 3312,      // Xi- MC code
        kLambdaCode = 3122,  // Lambda MC code
        kProtonCode = 2212,  // Proton+ MC code
        kPionCode = 211
    };  // Pion+ MC code
    // PN = unlike sign, PP and NN are like signs

    AliAnalysisTaskXi1530();
    AliAnalysisTaskXi1530(const char* name, const char* option);
    AliAnalysisTaskXi1530(const AliAnalysisTaskXi1530& ap);
    AliAnalysisTaskXi1530& operator=(const AliAnalysisTaskXi1530& ap);
    ~AliAnalysisTaskXi1530();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t*);
    virtual void Terminate(Option_t*);

    void SetOption(char* option) { fOption = option; }
    void SetFilterBit(UInt_t filterbit) { fFilterBit = filterbit; }
    void SetMixing(Bool_t setmixing) { fsetmixing = setmixing; }
    void SetIsAA(Bool_t isaa) { IsAA = isaa; }
    void SetIsMC(Bool_t ismc) { IsMC = ismc; }
    void SetnMix(Int_t nMix) { fnMix = nMix; }
    void SetHighMult(Bool_t highmult) { IsHighMult = highmult; }
    void SetIsPrimaryMC(Bool_t isprimarymc) { IsPrimaryMC = isprimarymc; }
    void SetQA(Bool_t QAinput) { fQA = QAinput; }
    void SetQAEvent(Bool_t QAinput) { IsQAEvent = QAinput; }
    void SetQAPID(Bool_t QAinput) { IsQAPID = QAinput; }
    void SetQAInvMass(Bool_t QAinput) { IsQAInvMass = QAinput; }

    void SetXiCPAstudy(Bool_t input) { fCPAstudy = input; }
    void SetINEL(Bool_t input) { isINEL = input; }
    // Set Functions for the cut study & Systematic study
    void SetTPCNsigXi1530PionCut(Int_t fSysOption, Double_t nXi1530PionCut) {
        if (fSysOption == 0)
            fTPCNsigXi1530PionCut = nXi1530PionCut;
        if (fSysOption == 1)
            fTPCNsigXi1530PionCut_loose = nXi1530PionCut;
        if (fSysOption == 2)
            fTPCNsigXi1530PionCut_tight = nXi1530PionCut;
    }
    void SetTPCNsigLambdaProtonCut(Int_t fSysOption,
                                   Double_t nLambdaProtonCut) {
        if (fSysOption == 0)
            fTPCNsigLambdaProtonCut = nLambdaProtonCut;
        if (fSysOption == 1)
            fTPCNsigLambdaProtonCut_loose = nLambdaProtonCut;
        if (fSysOption == 2)
            fTPCNsigLambdaProtonCut_tight = nLambdaProtonCut;
    }
    void SetTPCNsigLambdaPionCut(Int_t fSysOption, Double_t nLambdaPionCut) {
        if (fSysOption == 0)
            fTPCNsigLambdaPionCut = nLambdaPionCut;
        if (fSysOption == 1)
            fTPCNsigLambdaPionCut_loose = nLambdaPionCut;
        if (fSysOption == 2)
            fTPCNsigLambdaPionCut_tight = nLambdaPionCut;
    }
    void SetTPCNsigBachelorPionCut(Int_t fSysOption,
                                   Double_t nBachelorPionCut) {
        if (fSysOption == 0)
            fTPCNsigBachelorPionCut = nBachelorPionCut;
        if (fSysOption == 1)
            fTPCNsigBachelorPionCut_loose = nBachelorPionCut;
        if (fSysOption == 2)
            fTPCNsigBachelorPionCut_tight = nBachelorPionCut;
    }
    void SetTPCNsigXiCut(Int_t fSysOption, Double_t nXiCut) {
        if (fSysOption == 0) {
            fTPCNsigLambdaProtonCut = nXiCut;
            fTPCNsigLambdaPionCut = nXiCut;
            fTPCNsigBachelorPionCut = nXiCut;
        }
        if (fSysOption == 1) {
            fTPCNsigLambdaProtonCut_loose = nXiCut;
            fTPCNsigLambdaPionCut_loose = nXiCut;
            fTPCNsigBachelorPionCut_loose = nXiCut;
        }
        if (fSysOption == 2) {
            fTPCNsigLambdaProtonCut_tight = nXiCut;
            fTPCNsigLambdaPionCut_tight = nXiCut;
            fTPCNsigBachelorPionCut_tight = nXiCut;
        }
    }
    void SetXi1530PionEtaCut(Double_t nXi1530PionEtaCut) {
        fXi1530PionEtaCut = nXi1530PionEtaCut;
    }
    void SetXiEtaCut(Double_t nXiEtaCut) { fXiEtaCut = nXiEtaCut; }
    void SetXi1530PionZVertexCut(Int_t fSysOption,
                                 Double_t nXi1530PionZVertexCut) {
        if (fSysOption == 0)
            fXi1530PionZVertexCut = nXi1530PionZVertexCut;
        if (fSysOption == 1)
            fXi1530PionZVertexCut_loose = nXi1530PionZVertexCut;
        if (fSysOption == 2)
            fXi1530PionZVertexCut_tight = nXi1530PionZVertexCut;
    }
    void SetDCADist_LambdaDaughtersCut(Int_t fSysOption,
                                       Double_t nDCADist_LambdaDaughtersCut) {
        if (fSysOption == 0)
            fDCADist_LambdaDaughtersCut = nDCADist_LambdaDaughtersCut;
        if (fSysOption == 1)
            fDCADist_LambdaDaughtersCut_loose = nDCADist_LambdaDaughtersCut;
        if (fSysOption == 2)
            fDCADist_LambdaDaughtersCut_tight = nDCADist_LambdaDaughtersCut;
    }
    void SetDCADist_XiDaughtersCut(Int_t fSysOption,
                                   Double_t nDCADist_XiDaughtersCut) {
        if (fSysOption == 0)
            fDCADist_XiDaughtersCut = nDCADist_XiDaughtersCut;
        if (fSysOption == 1)
            fDCADist_XiDaughtersCut_loose = nDCADist_XiDaughtersCut;
        if (fSysOption == 2)
            fDCADist_XiDaughtersCut_tight = nDCADist_XiDaughtersCut;
    }

    void SetDCADist_Lambda_PVCut(Int_t fSysOption,
                                 Double_t nDCADist_Lambda_PVCut) {
        if (fSysOption == 0)
            fDCADist_Lambda_PVCut = nDCADist_Lambda_PVCut;
        if (fSysOption == 1)
            fDCADist_Lambda_PVCut_loose = nDCADist_Lambda_PVCut;
        if (fSysOption == 2)
            fDCADist_Lambda_PVCut_tight = nDCADist_Lambda_PVCut;
    }

    void SetV0CosineOfPointingAngleCut(Int_t fSysOption,
                                       Double_t nV0CosineOfPointingAngleCut) {
        if (fSysOption == 0)
            fV0CosineOfPointingAngleCut = nV0CosineOfPointingAngleCut;
        if (fSysOption == 1)
            fV0CosineOfPointingAngleCut_loose = nV0CosineOfPointingAngleCut;
        if (fSysOption == 2)
            fV0CosineOfPointingAngleCut_tight = nV0CosineOfPointingAngleCut;
    }
    void SetCascadeCosineOfPointingAngleCut(
        Int_t fSysOption,
        Double_t nCascadeCosineOfPointingAngleCut) {
        if (fSysOption == 0)
            fCascadeCosineOfPointingAngleCut = nCascadeCosineOfPointingAngleCut;
        if (fSysOption == 1)
            fCascadeCosineOfPointingAngleCut_loose =
                nCascadeCosineOfPointingAngleCut;
        if (fSysOption == 2)
            fCascadeCosineOfPointingAngleCut_tight =
                nCascadeCosineOfPointingAngleCut;
    }
    void SetXiMassWindowCut(Int_t fSysOption, Double_t nXiMassWindowCut) {
        if (fSysOption == 0)
            fXiMassWindowCut = nXiMassWindowCut;
        if (fSysOption == 1)
            fXiMassWindowCut_loose = nXiMassWindowCut;
        if (fSysOption == 2)
            fXiMassWindowCut_tight = nXiMassWindowCut;
    }
    void SetXi1530RapidityCut_high(Double_t nXi1530RapidityCut) {
        fXi1530RapidityCut_high = nXi1530RapidityCut;
    }
    void SetXi1530RapidityCut_low(Double_t nXi1530RapidityCut) {
        fXi1530RapidityCut_low = nXi1530RapidityCut;
    }
    void SetXi1530MultiplicityEstimator(TString fEstimator) {
        MultiplicityEstimator = fEstimator;
    }
    void SetXiSysTrackCut(Bool_t cutoption) { fsetXiSysTrackCut = cutoption; }
    void SetSystematics(Bool_t fSystematics) { fsetsystematics = fSystematics; }

    void SetExoticFinder(Bool_t exotic) {
        fExoticFinder = exotic;
        fXiMassWindowCut = 0.015;
        fTPCNsigXi1530PionCut = 1.5;
    }
    void SetExoticFinder2(Bool_t exotic) {
        fExoticFinder2 = exotic;
    }
    void SetUseSimpleEventCut(Bool_t fInput) { fSimplieEventCut = fInput; };

    Bool_t GoodTracksSelection();
    Bool_t GoodCascadeSelection();
    void FillTracks();
    void FillTracksAOD();

    Bool_t IsTrueXi1530(AliESDcascade* Xi, AliVTrack* pion);
    Bool_t IsTrueXi1530AOD(AliAODcascade* Xi, AliVTrack* pion);
    Bool_t IsTrueXi(UInt_t xiIndex);
    void FillMCinput(AliMCEvent* fMCEvent, Int_t check);
    void FillMCinputAOD(AliMCEvent* fMCEvent, Int_t check);
    void FillMCinputdXi(AliMCEvent* fMCEvent, Int_t check);
    void FillMCinputdXiAOD(AliMCEvent* fMCEvent, Int_t check);
    double GetTPCnSigma(AliVTrack* track, AliPID::EParticleType type);
    void GetImpactParam(AliVTrack* track, Float_t p[2], Float_t cov[3]);

    TAxis AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax);
    TAxis AxisVar(TString name, std::vector<Double_t> bin);
    TAxis AxisLog(TString name,
                  int nbin,
                  Double_t xmin,
                  Double_t xmax,
                  Double_t xmin0);
    TAxis AxisStr(TString name, std::vector<TString> bin);
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

    AliEventCuts fEventCuts;  // Event cuts

   private:
    typedef std::vector<AliVTrack*> tracklist;
    typedef std::deque<tracklist> eventpool;
    typedef std::vector<std::vector<eventpool>> mixingpool;

    TString fOption;

    AliESDtrackCuts* fTrackCuts = nullptr;   //!
    AliVEvent* fEvt = nullptr;               //!
    UInt_t fFilterBit = 32;

    Double_t fCent = -1;
    Double_t ftrackmult = -1;
    std::vector<UInt_t> goodtrackindices;    //!
    std::vector<UInt_t> goodtrackfullindices;    //!
    std::vector<UInt_t> goodcascadeindices;  //!

    AliPIDResponse* fPIDResponse = nullptr;  //!
    AliPIDCombined* fPIDCombined = nullptr;  //!
    // Histograms below are main

    mixingpool fEMpool;  //!
    TAxis binCent;       //!
    TAxis binZ;          //!
    Int_t fnMix = 10;
    Int_t centbin = -1;
    Int_t zbin = -1;
    std::vector<TString> SysCheck;
    TAxis binSystematics;

    Double_t fTPCNsigXi1530PionCut_loose = 3.5;
    Double_t fTPCNsigXi1530PionCut = 3.0;
    Double_t fTPCNsigXi1530PionCut_tight = 2.5;

    Double_t fTPCNsigLambdaProtonCut_loose = 3.5;
    Double_t fTPCNsigLambdaProtonCut = 3.0;
    Double_t fTPCNsigLambdaProtonCut_tight = 2.5;

    Double_t fTPCNsigLambdaPionCut_loose = 3.5;
    Double_t fTPCNsigLambdaPionCut = 3.0;
    Double_t fTPCNsigLambdaPionCut_tight = 2.5;

    Double_t fTPCNsigBachelorPionCut_loose = 3.5;
    Double_t fTPCNsigBachelorPionCut = 3.0;
    Double_t fTPCNsigBachelorPionCut_tight = 2.5;

    Double_t fXi1530PionEtaCut = 0.8;
    Double_t fXiEtaCut = 0.8;

    Double_t fXi1530PionZVertexCut_loose = 2.5;
    Double_t fXi1530PionZVertexCut = 2.0;
    Double_t fXi1530PionZVertexCut_tight = 1.5;

    Double_t fDCADist_LambdaDaughtersCut_loose = 1.5;
    Double_t fDCADist_LambdaDaughtersCut = 1.4;
    Double_t fDCADist_LambdaDaughtersCut_tight = 1.3;

    Double_t fDCADist_XiDaughtersCut_loose = 2.0;
    Double_t fDCADist_XiDaughtersCut = 1.6;
    Double_t fDCADist_XiDaughtersCut_tight = 1.2;

    Double_t fDCADist_Lambda_PVCut_loose = 0.05;
    Double_t fDCADist_Lambda_PVCut = 0.07;
    Double_t fDCADist_Lambda_PVCut_tight = 0.09;

    Double_t fV0CosineOfPointingAngleCut_loose = 0.96;
    Double_t fV0CosineOfPointingAngleCut = 0.97;
    Double_t fV0CosineOfPointingAngleCut_tight = 0.98;

    Double_t fCascadeCosineOfPointingAngleCut_loose = 0.96;
    Double_t fCascadeCosineOfPointingAngleCut = 0.97;
    Double_t fCascadeCosineOfPointingAngleCut_tight = 0.98;

    Double_t fXiMassWindowCut_loose = 0.008;
    Double_t fXiMassWindowCut = 0.007;
    Double_t fXiMassWindowCut_tight = 0.006;

    Double_t fXi1530RapidityCut_low = -0.5;
    Double_t fXi1530RapidityCut_high = 0.5;

    Bool_t fsetXiSysTrackCut = kFALSE;
    Bool_t fsetsystematics = kFALSE;
    Bool_t fsetmixing = kTRUE;
    Bool_t IsAA = kFALSE;
    Bool_t IsMC = kFALSE;
    Bool_t IsPS = kFALSE;
    Bool_t IsINEL0Rec = kFALSE;
    Bool_t IsINEL0True = kFALSE;
    Bool_t IsHighMult = kFALSE;
    Bool_t IsPrimaryMC = kTRUE;
    Bool_t IsQAEvent = kTRUE;
    Bool_t IsQAPID = kTRUE;
    Bool_t IsQAInvMass = kTRUE;
    Bool_t fQA = kTRUE;
    Bool_t fExoticFinder = kFALSE;
    Bool_t fExoticFinder2 = kFALSE;
    Bool_t fSimplieEventCut = kFALSE;
    Bool_t IsNano = kFALSE;
    Bool_t IsAOD = kFALSE;
    Bool_t fCPAstudy = kFALSE;
    Bool_t isINEL = kFALSE;
    TString MultiplicityEstimator = "V0M";
    THistManager* fHistos = nullptr;   //!
    TClonesArray* fMCArray = nullptr;  //!
    AliMCEvent* fMCEvent = nullptr;    //!

    Double_t lPosPV[3];
    Double_t bField = 999;
    ClassDef(AliAnalysisTaskXi1530, 33);
    // 1: Frist version
    // 2: Add Track cut2 for the Xi daughter particles
    // 3: Add FillMixingPool function
    // 4: Add Cut parameters to header and add "Set" fuction for cut
    // study&Systematic study
    // 5: include AliAnalysisTaskSE.h to avoid compile problem.
    // 6: Add IsPrimaryMC option for MC study
    // 7: Add Systematics option for systematics study
    // 8: Change default systmatics on TPC PID for Xi daughter,
    // Add Setter function for that.
    // 9: Hot fix for the AliPhysics building
    // 10: Add NoQA option to reduce output file size
    // 11: Not using AliStack informed by DPG and BTG coordination
    // 12: Add Exotic option and change default systematic vaule
    // 13: Change Lambda PV Cut default values(reversed)
    // 14: Modify Fill MC input function for the Signal Loss
    // 15: Separate Rapidity cut option for pA analysis.
    // 16: Add option for using different Multiplicity Estimators(Default: V0M).
    // 17: Add Simple event cut option to use AliMultSelection cut only
    // 18: Enable AOD functionality
    // 19: Update default filterbit
    // 20: Remove RunTable Class.
    // 21: Use AliEventCuts
    // 22: Remove some function thanks to AliEventCuts
    // 23: Integrate true Xi checker(ESD/AOD)
    // 24: Add IsNano boolean variable
    // 25: Remove Fill event mix function
    // 26: try to update the version only.
    // 27: add Xi CPA study option // Change PV variable type
    // 28: Add INEL variable
    // 29: Add Special finder
    // 30: Remove unused variables and reformatting
    // 31: Add option to remove PIDqa plots
    // 32: Add goodtrackfullindices for further trackloop
    // 33: Modify for the nanoCheck
};

#endif
