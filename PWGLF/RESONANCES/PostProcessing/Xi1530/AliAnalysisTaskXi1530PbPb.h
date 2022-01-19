#ifndef AliAnalysisTaskXi1530PBPB_H
#define AliAnalysisTaskXi1530PBPB_H

#include <THnSparse.h>
#include <TNtupleD.h>

#include <deque>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskTrackMixer.h"
#include "AliEventCuts.h"
class THistManager;
class AliPIDResponse;
class AliESDtrackCuts;

class AliAnalysisTaskXi1530PbPb : public AliAnalysisTaskSE {
   public:
    AliAnalysisTaskXi1530PbPb();
    AliAnalysisTaskXi1530PbPb(const char *name, Bool_t MCcase);
    virtual ~AliAnalysisTaskXi1530PbPb();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    void SetFillQAPlot(Bool_t input = kTRUE) { fFillQAPlot = input; }
    void SetMixerTask(AliAnalysisTaskTrackMixer *mixerTask) {
        fMixingPool = mixerTask;
    }
    void SetUseBuiltinMixer(Bool_t builtinmixer) {
        fUseBuiltinMixer = builtinmixer;
    }
    void SetMixing(Bool_t setmixing) { fSetMixing = setmixing; }
    void SetnMix(Int_t nMix) { fnMix = nMix; }
    void SetIsPrimaryMC(Bool_t isprimarymc) { fIsPrimaryMC = isprimarymc; }
    void SetINEL(Bool_t input) { fIsINEL = input; }
    void SetHighMult(Bool_t input) { fIsHM = input; }
    void SetFillnTuple(Bool_t fillntuple = kTRUE) { fFillnTuple = fillntuple; }
    void SetTurnOnExoFinder(Bool_t input = kTRUE) { fExoticFinder = input; }
    void SetLambdaCPAtoXi(Bool_t input = kTRUE) { fLambdaCPAtoXi = input; }

    // Setter for cut variables
    void SetFilterbitXi1530Pion(Double_t lParameter) { fFilterBit = lParameter; }
    void SetMaxNsigXi1530Pion(Double_t lParameter) { fTPCNsigXi1530PionCut = lParameter; }
    void SetMaxEtaXi1530Pion(Double_t lParameter) { fXi1530PionEtaCut = lParameter; }
    void SetMaxVertexZXi1530Pion(Double_t lParameter) { fXi1530PionPVzCut = lParameter; }
    void SetMaxPVrSigmaXi1530Pion(Double_t lParameter) { fXi1530PionPVrSigmaCut = lParameter; }

    void SetMaxNsigV0Proton(Double_t lParameter) { fTPCNsigLambdaProtonCut = lParameter; }
    void SetMaxNsigV0Pion(Double_t lParameter) { fTPCNsigLambdaPionCut = lParameter; }
    void SetMaxNsigBachelorPion(Double_t lParameter) { fTPCNsigBachelorPionCut = lParameter; }
    void SetMaxDCAV0daughters(Double_t lParameter) { fDCALambdaDaughtersCut = lParameter; }
    void SetMaxDCAXidaughters(Double_t lParameter) { fDCAXiDaughtersCut = lParameter; }
    void SetMinDCAPVV0(Double_t lParameter) { fDCALambdaPVCut = lParameter; }
    void SetMinDCAPVV0Proton(Double_t lParameter) { fDCALambdaProtonPVCut = lParameter; }
    void SetMinDCAPVV0Pion(Double_t lParameter) { fDCALambdaPionPVCut = lParameter; }
    void SetMinDCAPVBachelorPion(Double_t lParameter) { fDCABachelorPionPVCut = lParameter; }
    void SetMinDCAPVXi(Double_t lParameter) { fDCAXiPVCut = lParameter; }
    void SetMinCPAV0(Double_t lParameter) { fV0CPACut = lParameter; }
    void SetMinCPAXi(Double_t lParameter) { fXiCPACut = lParameter; }
    void SetMaxEtaXi(Double_t lParameter) { fXiEtaCut = lParameter; }
    void SetLowRadiusV0(Double_t lParameter) { fLambdaLowRadius = lParameter; }
    void SetHighRadiusV0(Double_t lParameter) { fLambdaHighRadius = lParameter; }
    void SetLowRadiusXi(Double_t lParameter) { fXiLowRadius = lParameter; }
    void SetHighRadiusXi(Double_t lParameter) { fXiHighRadius = lParameter; }
    void SetMaxMassWindowV0(Double_t lParameter) { fV0MassWindowCut = lParameter; }
    void SetMaxMassWindowXi(Double_t lParameter) { fXiMassWindowCut = lParameter; }
    void SetExoticMaxOpenAngle(Double_t lParameter) { fExoticMaxOpenAngle = lParameter; }

    void SetXi1530RapidityCutHigh(Double_t lParameter) {
        fXi1530YCutHigh = lParameter;
    }
    void SetXi1530RapidityCutLow(Double_t lParameter) {
        fXi1530YCutLow = lParameter;
    }

    Bool_t GoodTracksSelection();
    Bool_t GoodCascadeSelection();
    void FillTracks();
    void FillNtuples();
    void FillMCinput(AliMCEvent *fMCEvent, int Fillbin = 0);
    void FillTrackToEventPool();
    Bool_t IsTrueXi1530(UInt_t v0, UInt_t pion);
    double GetTPCnSigma(AliVTrack *track, AliPID::EParticleType type);
    void GetImpactParam(AliVTrack *track, Float_t p[2], Float_t cov[3]);

    Bool_t IsSelectedTPCGeoCut(AliAODTrack *track);
    Bool_t IsSelectedTPCGeoCut(AliESDtrack *track);
    void SetCutOpen();

    // helper
    THnSparse *CreateTHnSparse(TString name,
                               TString title,
                               Int_t ndim,
                               std::vector<TAxis> bins,
                               Option_t *opt = "");
    THnSparse *CreateTHnSparse(TString name,
                               TString title,
                               TString templ,
                               Option_t *opt = "");
    Long64_t FillTHnSparse(TString name,
                           std::vector<Double_t> x,
                           Double_t w = 1.);
    Long64_t FillTHnSparse(THnSparse *h,
                           std::vector<Double_t> x,
                           Double_t w = 1.);
    TAxis AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax);
    TAxis AxisVar(TString name, std::vector<Double_t> bin);
    TAxis AxisStr(TString name, std::vector<TString> bin);

    AliEventCuts fEventCuts;  // Event cuts
    // TPC GeoCut
    Bool_t fCheckTPCGeo;                   //
    Double_t fTPCActiveLengthCutDeltaY;    //
    Double_t fTPCActiveLengthCutDeltaZ;    //
    Double_t fRequireCutGeoNcrNclLength;   //
    Double_t fRequireCutGeoNcrNclGeom1Pt;  //
    Double_t fCutGeoNcrNclFractionNcr;     //
    Double_t fCutGeoNcrNclFractionNcl;     //

   private:
    typedef std::vector<AliVTrack *> tracklist;
    typedef std::deque<tracklist> eventpool;
    typedef std::vector<std::vector<eventpool> > mixingpool;

    AliESDtrackCuts *fTrackCuts;             //!
    AliPIDResponse *fPIDResponse;            //!
    AliAnalysisTaskTrackMixer *fMixingPool;  //!

    AliVEvent *fEvt;          //!
    AliMCEvent *fMCEvent;     //!
    THistManager *fHistos;    //!
    AliAODVertex *fVertex;    //!
    TNtupleD *fNtupleXi1530;  //!
    TClonesArray *fMCArray;   //!

    Bool_t fIsAOD;            //!
    Bool_t fIsNano;           //!
    Bool_t fSetMixing;        //
    Bool_t fUseBuiltinMixer;  //
    Bool_t fFillQAPlot;       //
    Bool_t fIsMC;             //
    Bool_t fIsPrimaryMC;      //
    Bool_t fFillnTuple;       //
    Bool_t fIsINEL;           //
    Bool_t fIsHM;             //
    Bool_t fLambdaCPAtoXi;    //
    Bool_t fExoticFinder;     //

    mixingpool fEMpool;  //!
    TAxis fBinCent;      //!
    TAxis fBinZ;         //!
    Double_t fPosPV[3];  //!
    Double_t fMagField;  //!

    Double_t fCent;  //!
    Int_t fnMix;     //!
    Int_t fCentBin;  //!
    Int_t fZbin;     //!

    // Pion cuts
    UInt_t fFilterBit;                //
    Double_t fTPCNsigXi1530PionCut;   //
    Double_t fXi1530PionEtaCut;       //
    Double_t fXi1530PionPVzCut;       //
    Double_t fXi1530PionPVrSigmaCut;  //

    // Xi cuts
    Double_t fTPCNsigLambdaProtonCut;  //
    Double_t fTPCNsigLambdaPionCut;    //
    Double_t fTPCNsigBachelorPionCut;  //
    Double_t fDCALambdaDaughtersCut;   //
    Double_t fDCAXiDaughtersCut;       //
    Double_t fDCALambdaPVCut;          //
    Double_t fDCAXiPVCut;              //
    Double_t fDCALambdaProtonPVCut;    //
    Double_t fDCALambdaPionPVCut;      //
    Double_t fDCABachelorPionPVCut;    //
    Double_t fV0CPACut;                //
    Double_t fXiCPACut;                //
    Double_t fXiEtaCut;                //
    Double_t fLambdaLowRadius;         //
    Double_t fLambdaHighRadius;        //
    Double_t fXiLowRadius;             //
    Double_t fXiHighRadius;            //
    Double_t fV0MassWindowCut;         //
    Double_t fXiMassWindowCut;         //

    // Xi1530 cut
    Double_t fExoticMaxOpenAngle;  //
    Double_t fXi1530YCutHigh;      //
    Double_t fXi1530YCutLow;       //

    // Good track/cascade vector array
    std::vector<UInt_t> fGoodTrackArray;  //
    std::vector<UInt_t> fGoodXiArray;     //

    ClassDef(AliAnalysisTaskXi1530PbPb, 2);
    // Code beautify, fix minor bugs
};

#endif