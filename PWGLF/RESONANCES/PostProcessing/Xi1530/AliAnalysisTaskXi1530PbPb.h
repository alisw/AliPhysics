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

struct StructXi1530PbPb
{
    Double32_t TPCNsigXi1530Pion;
    Double32_t TPCNsigLambdaProton;
    Double32_t TPCNsigLambdaPion;
    Double32_t TPCNsigBachelorPion;
    Double32_t Xi1530PionEta;
    Double32_t Xi1530PionPVz;
    Double32_t Xi1530PionXYVertexSigma;
    Double32_t DCALambdaDaughters;
    Double32_t DCAXiDaughters;
    Double32_t DCALambdaPV;
    Double32_t DCALambdaProtonPV;
    Double32_t DCALambdaPionPV;
    Double32_t DCABachelorPionPV;
    Double32_t DCAXiPV;
    Double32_t V0CPA;
    Double32_t XiCPA;
    Double32_t XiEta;
    Double32_t LambdaRadius;
    Double32_t XiRadius;
    Double32_t V0Mass;
    Double32_t XiMass;
    Double32_t Xi1530Pt;
    Double32_t Xi1530Mass;
    Double32_t Xi1530Radius;
    Double32_t Xi1530DCA;
    Double32_t Xi1530CPA;
    Double32_t Centrality;
    Double32_t zVertex;
    Int_t Antiflag;
};

struct StructXi1530PbPbMC : public StructXi1530PbPb
{
    float MCflag;
    float ptMC;
    Int_t Antiflag;
    bool isReconstructed;
};

class AliAnalysisTaskXi1530PbPb : public AliAnalysisTaskSE
{
public:
    AliAnalysisTaskXi1530PbPb();
    AliAnalysisTaskXi1530PbPb(const char *name, Bool_t MCcase);
    virtual ~AliAnalysisTaskXi1530PbPb();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    void SetFillQAPlot(Bool_t input = kTRUE) { fFillQAPlot = input; }
    void SetMixerTask(AliAnalysisTaskTrackMixer *mixerTask)
    {
        fMixingPool = mixerTask;
    }
    void SetUseBuiltinMixer(Bool_t builtinmixer)
    {
        fUseBuiltinMixer = builtinmixer;
    }
    void SetMixing(Bool_t setmixing) { fSetMixing = setmixing; }
    void SetnMix(Int_t nMix) { fnMix = nMix; }
    void SetIsPrimaryMC(Bool_t isprimarymc) { fIsPrimaryMC = isprimarymc; }
    void SetINEL(Bool_t input) { fIsINEL = input; }
    void SetHighMult(Bool_t input) { fIsHM = input; }
    void SetFillnTuple(Bool_t fillntuple = kTRUE) { fFillTree = fillntuple; }
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
    void SetMinDCAPVV0(Double_t lParameter) { fMinDCALambdaPVCut = lParameter; }
    void SetMaxDCAPVV0(Double_t lParameter) { fMaxDCALambdaPVCut = lParameter; }
    void SetMinDCAPVV0Proton(Double_t lParameter) { fMinDCALambdaProtonPVCut = lParameter; }
    void SetMinDCAPVV0Pion(Double_t lParameter) { fMinDCALambdaPionPVCut = lParameter; }
    void SetMinDCAPVBachelorPion(Double_t lParameter) { fMinDCABachelorPionPVCut = lParameter; }
    void SetMaxDCAPVV0Proton(Double_t lParameter) { fMaxDCALambdaProtonPVCut = lParameter; }
    void SetMaxDCAPVV0Pion(Double_t lParameter) { fMaxDCALambdaPionPVCut = lParameter; }
    void SetMaxDCAPVBachelorPion(Double_t lParameter) { fMaxDCABachelorPionPVCut = lParameter; }
    void SetMinDCAPVXi(Double_t lParameter) { fMinDCAXiPVCut = lParameter; }
    void SetMaxDCAPVXi(Double_t lParameter) { fMaxDCAXiPVCut = lParameter; }
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
    void SetSkipFillingHistogram(Double_t lParameter) { fSkipFillingHisto = lParameter; }

    void SetXi1530RapidityCutHigh(Double_t lParameter)
    {
        fXi1530YCutHigh = lParameter;
    }
    void SetXi1530RapidityCutLow(Double_t lParameter)
    {
        fXi1530YCutLow = lParameter;
    }
    void SetXi1530MassCutHigh(Double_t lParameter)
    {
        fXi1530MassHigh = lParameter;
    }
    void SetXi1530MassCutLow(Double_t lParameter)
    {
        fXi1530MassLow = lParameter;
    }

    Bool_t GoodTracksSelection();
    Bool_t GoodCascadeSelection();
    void FillTracks();
    void FillTree();
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

    AliEventCuts fEventCuts; // Event cuts
    // TPC GeoCut
    Bool_t fCheckTPCGeo;                  //
    Double_t fTPCActiveLengthCutDeltaY;   //
    Double_t fTPCActiveLengthCutDeltaZ;   //
    Double_t fRequireCutGeoNcrNclLength;  //
    Double_t fRequireCutGeoNcrNclGeom1Pt; //
    Double_t fCutGeoNcrNclFractionNcr;    //
    Double_t fCutGeoNcrNclFractionNcl;    //

private:
    typedef std::vector<AliVTrack *> tracklist;
    typedef std::deque<tracklist> eventpool;
    typedef std::vector<std::vector<eventpool>> mixingpool;

    AliESDtrackCuts *fTrackCuts;            //!
    AliPIDResponse *fPIDResponse;           //!
    AliAnalysisTaskTrackMixer *fMixingPool; //!

    AliVEvent *fEvt;                            //!
    AliMCEvent *fMCEvent;                       //!
    THistManager *fHistos;                      //!
    AliAODVertex *fVertex;                      //!
    TTree *fTree = nullptr;                     //!<! Tree for Xi1530
    StructXi1530PbPb *fTreeXi1530Rec = nullptr; //!<! Transient fRecXi1530
    StructXi1530PbPbMC fTreeXi1530Gen;
    TClonesArray *fMCArray; //!

    Bool_t fIsAOD;            //!
    Bool_t fIsNano;           //!
    Bool_t fSetMixing;        //
    Bool_t fUseBuiltinMixer;  //
    Bool_t fFillQAPlot;       //
    Bool_t fIsMC;             //
    Bool_t fIsPrimaryMC;      //
    Bool_t fFillTree;         //
    Bool_t fIsINEL;           //
    Bool_t fIsHM;             //
    Bool_t fLambdaCPAtoXi;    //
    Bool_t fExoticFinder;     //
    Bool_t fSkipFillingHisto; //

    mixingpool fEMpool; //!
    TAxis fBinCent;     //!
    TAxis fBinZ;        //!
    Double_t fPosPV[3]; //!
    Double_t fMagField; //!

    Double_t fCent; //!
    Int_t fnMix;    //!
    Int_t fCentBin; //!
    Int_t fZbin;    //!

    // Pion cuts
    UInt_t fFilterBit;               //
    Double_t fTPCNsigXi1530PionCut;  //
    Double_t fXi1530PionEtaCut;      //
    Double_t fXi1530PionPVzCut;      //
    Double_t fXi1530PionPVrSigmaCut; //

    // Xi cuts
    Double_t fTPCNsigLambdaProtonCut;  //
    Double_t fTPCNsigLambdaPionCut;    //
    Double_t fTPCNsigBachelorPionCut;  //
    Double_t fDCALambdaDaughtersCut;   //
    Double_t fDCAXiDaughtersCut;       //
    Double_t fMinDCALambdaPVCut;       //
    Double_t fMaxDCALambdaPVCut;       //
    Double_t fMinDCAXiPVCut;           //
    Double_t fMaxDCAXiPVCut;           //
    Double_t fMinDCALambdaProtonPVCut; //
    Double_t fMinDCALambdaPionPVCut;   //
    Double_t fMinDCABachelorPionPVCut; //
    Double_t fMaxDCALambdaProtonPVCut; //
    Double_t fMaxDCALambdaPionPVCut;   //
    Double_t fMaxDCABachelorPionPVCut; //
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
    Double_t fExoticMaxOpenAngle; //
    Double_t fXi1530YCutHigh;     //
    Double_t fXi1530YCutLow;      //
    Double_t fXi1530MassHigh;     //
    Double_t fXi1530MassLow;      //

    // Good track/cascade vector array
    std::vector<UInt_t> fGoodTrackArray; //
    std::vector<UInt_t> fGoodXiArray;    //

    ClassDef(AliAnalysisTaskXi1530PbPb, 4);
    // 2. Code beautify, fix minor bugs
    // 3. Add secondary vertex option
    // 4. Update nTuples to Tree
};

#endif