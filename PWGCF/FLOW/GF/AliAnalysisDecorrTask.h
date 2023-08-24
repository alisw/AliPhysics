#ifndef AliAnalysisDecorrTask_H
#define AliAnalysisDecorrTask_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliGFWWeights.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "TComplex.h"
#include "AliVParticle.h"
#include "TAxis.h"
#include "AliDecorrFlowCorrTask.h"
#include <vector>

class Taxis;
class AliVEvent;
class TProfile;
class TObjArray;
class AliDecorrFlowCorrTask;

class AliAnalysisDecorrTask : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisDecorrTask();
                                AliAnalysisDecorrTask(const char *name, bool IsMC, bool isOnTheFly);
        virtual                 ~AliAnalysisDecorrTask();

        virtual void            UserCreateOutputObjects();
        virtual void            NotifyRun();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        //Analysis setters
        void                    SetSampling(Bool_t sample, Int_t iNum, Bool_t tracks = kFALSE) { fSampling = sample; fNumSamples = iNum; fRedTracks = tracks; }      //Use jack-knife resampling
        void                    SetFillQA(Bool_t fill = kTRUE) { fFillQA = fill; }
        //event selection
        void                    SetTriggerType(UInt_t newval) { fTrigger = newval; }
        void                    SetRejectAddPileUp(Bool_t use = kTRUE) { fEventRejectAddPileUp = use; }
        void                    SetPileupCut(Int_t cut) { fCentralPileupCut = cut; }
        void                    SetCentralityEst(TString est){ fCentEstimator = est; }
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        void                    SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
        void                    SetCentBin(Int_t nbins, Double_t *bins) { fCentAxis->Set(nbins,bins); }
        void                    SetCentLim(Double_t min, Double_t max) { fCentMin = min; fCentMax = max; } //Not used yet
        void                    SetPtBins(Int_t nbins, Double_t *bins) { fPtAxis->Set(nbins, bins); }
        void                    SetRequireHighPtTracks(Int_t Ntracks, Double_t ptcut) { fRequireHighPtTracks = true; fNHighPtTracks = Ntracks; fHighPtCut = ptcut; }
        void                    SetRequirePrimariesAndCh(Bool_t newval) { bOnlyPrimariesAndCh = newval; }
        AliEventCuts            fEventCuts;
        //track selection
        void                    SetDCAzMax(Double_t dcaz) {  fCutDCAzMax = dcaz; }
        void                    SetDCAxyMax(Double_t dcaxy) {  fCutDCAxyMax = dcaxy; }
        void                    SetChi2(Double_t chi2) { fChi2Cut = chi2; }
        void                    SetNumTPCclsMin(UShort_t tpcCls) { fCutNumTPCclsMin = tpcCls; }
        void                    SetUseLikeSign(Bool_t use, Int_t sign) { bUseLikeSign = use; iSign = sign; }
        void                    SetChargedTrackFilterBit(UInt_t filter) { fCutChargedTrackFilterBit = filter; } //Not implemented
        void                    SetTrackPercentage(Double_t percent) { fTrackprevent = percent; }
        //Flow selection
        void                    AddCorr(std::vector<Int_t> harms, std::vector<Double_t> gaps = std::vector<Double_t>(), Bool_t doRef = kTRUE, Bool_t doDiff = kTRUE, Bool_t doPtA = kFALSE, Bool_t doPtRef = kFALSE, Bool_t doPtB = kFALSE) { fVecCorrTask.push_back(new AliDecorrFlowCorrTask(doRef, doDiff, doPtA, doPtRef, doPtB, harms, gaps)); }
        void                    SetPOIsPt(Double_t min, Double_t max) { fPOIsPtmin = min; fPOIsPtmax = max; }
        void                    SetRFPsPt(Double_t min, Double_t max) { fRFPsPtMin = min; fRFPsPtMax = max; }
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }
        void                    SetEtaBins(Int_t bins) { fEtaBinNum = bins; }
        void                    SetPhiBins(Int_t bins) { fPhiBinNum = bins; }
        void                    SetFillWeights(Bool_t fill) { fFillWeights = fill; }    //Fill histograms for weights calculations
        void                    SetUseWeightsOne(Bool_t useOne) { fUseWeightsOne = useOne; }
        void                    SetCurrSystFlag(int sys) { fCurrSystFlag = sys; }
        void                    SetRequireTwoPart(Bool_t req) { fRequireTwoPart = req; }
        void                    OverrideMCFlag(Bool_t newval) { fIsMC = newval; };
        void                    OverrideCentrality(double val) { fOverrideCentrality = val; }
        void                    DebugTrain(bool dbg) { fDebugTrain = dbg; }
        struct StoreRef
        {
            double val;
        };
        struct StorePtA
        {
            double val[3][10];
        };
        struct StorePtB
        {
            double val[10][10];
        };

    private:
        static const Int_t      fNumHarms = 8;             // Used to be 13 for 4-pc without gap. Reduced for cpu. Maximum harmonics length of flow vector array
        static const Int_t      fNumPowers = 2;             // maximum weight power length of flow vector array
        static const Int_t      NcentBinMax = 11;           //
        static const Int_t      NPtBinMax = 30;             //

        TList*                  fFlowList;                //! output list
        //TList*                fFlowWeights;             //! 
        TList*                  fQA;                      //!
        TList*                  lSampleList;                //!


        AliMCEvent*           getMCEvent();  
        bool                  InitTask();
        bool                  LoadWeights();
        double                GetWeights(double dPhi, double dEta, double dVz);
        void                  CreateProfiles(AliDecorrFlowCorrTask* task, double* lCentEdges, double* lPtEdges);
        bool                  IsEventSelected(Double_t lCent);
        bool                  CheckTrigger(double lCent);
        bool                  IsEventRejectedAddPileUp(const int fPileupCut) const;
        bool                  IsTrackSelected(const AliAODTrack* track) const;
        int                   GetSamplingIndex() const;
        double                GetCentralityFromImpactParameter();
        double*               GetBinsFromAxis(TAxis *inax);
        void                  InitTypes(const AliDecorrFlowCorrTask* const task);

        std::vector<StoreRef> refData;
        std::vector<StorePtA> ptaData;
        std::vector<StorePtB> ptbData;

        //Weights and QA hists
        AliGFWWeights*          fWeights;                   //!
        TList*                  fWeightList;                //!
        TH1D*                   fV0MMulti;                  //!
        TH1D*                   fMulti;                     //!
        TH1I*                   hITSclsB;                   //!
        TH1I*                   hTPCclsB;                   //!
        TH1D*                   hTPCchi2B;                  //!    
        TH2D*                   hDCAxyB;                      //! 
        TH1I*                   hITSclsA;                   //!
        TH1I*                   hTPCclsA;                   //!
        TH1D*                   hTPCchi2A;                  //!    
        TH2D*                   hDCAxyA;                      //!    
        TH3D*                   hPtPhiEtaB;                 //!
        TH3D*                   hPtPhiEtaA;                 //!
        TH1D*                   hNumTracksB;                //!
        TH1D*                   hNumTracksA;                //!
        TH1D*                   hNumHighPtTracksA;          //!
        TH1I*                   hCharge;                    //!
        TH1D*                   fhEventSel;                 //!

        void                    FillWeights();   
        double                  getCentrality();     
        
        //Flow methods
        bool                    IsWithinRP(const AliVParticle* track) const;
        bool                    IsWithinPOI(const AliVParticle* track) const;
        void                    FillRPvectors(const AliDecorrFlowCorrTask* const task);
        int                     FillPOIvectors(const AliDecorrFlowCorrTask* const task, const double dPtLow, const double dPtHigh); 
        void                    FillPtBvectors(const AliDecorrFlowCorrTask* const task, const double dPtLow, const double dPtHigh); 
        void                    CalculateCorrelations(const AliDecorrFlowCorrTask* const task, const double &lCent, const double &lpta, const double &lptb, Bool_t bRef, Bool_t bDiff, Bool_t bPtA, Bool_t bPtRef, Bool_t bPtB);
        void                    CalculateGFW(TComplex c[2], const AliDecorrFlowCorrTask* const task, int type, bool gap, bool switchreg = false);
        //void                    RunningCovariance(double x, double y, const double &lCent, const double &lpta);
        void                    FillProfiles(TComplex c[2], const double &lCent, TString suffix);
        void                    FillProfiles(TComplex c[2], const double &lCent, const double &lpta, TString suffix);
        void                    FillProfiles(TComplex c[2], const double &lCent, const double &lpta, const double lptb, TString suffix);

        //Flow vectors
        TComplex pvector[fNumHarms][fNumPowers];
        TComplex pvector10M[fNumHarms][fNumPowers];
        TComplex pvector10P[fNumHarms][fNumPowers];
        TComplex qvector[fNumHarms][fNumPowers];
        TComplex pvectorPtB[fNumHarms][fNumPowers];
        TComplex qvectorPtB[fNumHarms][fNumPowers];
        TComplex pvectorPtB10M[fNumHarms][fNumPowers];
        TComplex pvectorPtB10P[fNumHarms][fNumPowers];
        TComplex Qvector[fNumHarms][fNumPowers];
        TComplex Qvector10M[fNumHarms][fNumPowers];
        TComplex Qvector10P[fNumHarms][fNumPowers];

        TComplex Q(int n, int p);
        TComplex QGap10M(int n, int p);
        TComplex QGap10P(int n, int p);
        TComplex p(int n, int p);
        TComplex pGap10M(int n, int p);
        TComplex pGap10P(int n, int p);
        TComplex pPtBGap10M(int n, int p);
        TComplex pPtBGap10P(int n, int p);
        TComplex q(int n, int p);
        TComplex qGap10M(int n, int p);
        TComplex qGap10P(int n, int p);
        TComplex pPtA(int n, int p);
        TComplex pPtB(int n, int p);
        TComplex qPtA(int n, int p);
        TComplex qPtB(int n, int p);
        void ResetFlowVector(TComplex (&array)[fNumHarms][fNumPowers]);

        TComplex Two(int n1, int n2);
        TComplex TwoGap10(int n1, int n2);
        TComplex TwoDiff(int n1, int n2);
        TComplex TwoDiffGap10M(int n1, int n2);
        TComplex TwoDiffGap10P(int n1, int n2);
        TComplex TwoDiff_Pt(int n1, int n2);
        TComplex TwoDiffGap10_Pt(int n1, int n2);
        TComplex TwoDiffGap10_PtB(int n1, int n2);
        TComplex TwoDiff_PtA(int n1, int n2);
        TComplex TwoDiff_PtB(int n1, int n2);
        //Subevent methods
        TComplex Two_SubP(int n1, int n2);
        TComplex Two_SubM(int n1, int n2);
        TComplex TwoDiff_SubP(int n1, int n2);
        TComplex TwoDiff_SubM(int n1, int n2);
        TComplex TwoDiff_SubM_PtA(int n1, int n2);
        TComplex TwoDiff_SubP_PtA(int n1, int n2);
        TComplex TwoDiff_SubM_PtB(int n1, int n2);
        TComplex TwoDiff_SubP_PtB(int n1, int n2);
        TComplex TwoDiff_SubP_PtA_PtB(int n1, int n2);
        TComplex TwoDiff_SubM_PtA_PtB(int n1, int n2);

        TComplex TwoDiffGap10M_PtA(int n1, int n2);
        TComplex TwoDiffGap10P_PtB(int n1, int n2);
        TComplex TwoDiff_PtA_PtB(int n1, int n2);
        TComplex TwoDiffGap10_PtA_PtB(int n1, int n2);
        TComplex Three(int n1, int n2, int n3);
        TComplex ThreeGapP(int n1, int n2, int n3);
        TComplex ThreeGapM(int n1, int n2, int n3);
        TComplex ThreeDiff(int n1, int n2, int n3);
        TComplex ThreeDiffGapP(int n1, int n2, int n3);
        TComplex ThreeDiffGapM(int n1, int n2, int n3);
        TComplex Four(int n1, int n2, int n3, int n4);
        TComplex FourGap10(int n1, int n2, int n3, int n4);
        TComplex FourDiff(int n1, int n2, int n3, int n4);
        TComplex Four_2Diff_2Ref(int n1, int n2, int n3, int n4);
        TComplex FourGapM_2Diff_2Ref(int n1, int n2, int n3, int n4);
        TComplex FourGapP_2Diff_2Ref(int n1, int n2, int n3, int n4);
        TComplex FourGap_2Diff_2Ref_OS(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10P(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10M(int n1, int n2, int n3, int n4);
        TComplex FourDiff_PtA_PtA(int n1, int n2, int n3, int n4);
        TComplex FourDiff_PtA_PtB(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10_PtA_PtA(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10M_PtA_PtB(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10P_PtA_PtB(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10_OS_PtA_PtB(int n1, int n2, int n3, int n4);
        TComplex Five(int n1, int n2, int n3, int n4, int n5);
        TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6);
        TComplex SixDiff(int n1, int n2, int n3, int n4, int n5, int n6);
        TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
        TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);

    

        AliAnalysisDecorrTask(const AliAnalysisDecorrTask&); // not implemented
        AliAnalysisDecorrTask& operator=(const AliAnalysisDecorrTask&); // not implemented

        //Array lengths and constants
        Int_t                   fIndexSampling;
        Int_t                   fTaskCounter;
        AliAODEvent*            fAOD;                       //! input event
        Bool_t                  fIsMC;
        AliMCEvent*             fMCEvent;                   //! MC event
        AliVEvent*              fEvent;                     //!
        Bool_t                  fInitTask;                  //Initialization
        Bool_t                  fOnTheFly;                  
        Bool_t                  fDebugTrain;
        Double_t                fImpactParameterMC;
        std::vector<AliDecorrFlowCorrTask*>    fVecCorrTask;   //
        TString                 fCorrName;                     //
        TString                 fCorrLabel;                    //
        int                     fCorrOrder;                    //
        std::map<double,double>      centralitymap;  

        
        //cuts & selection: Analysis
        Bool_t                  fSampling;      //Bootstrapping sampling
        Bool_t                  fRedTracks;
        Double_t                fTrackprevent;
        Bool_t                  fFillQA;        //Fill QA histograms
        //cuts & selection: events
        UInt_t                  fTrigger;
        Bool_t                  fEventRejectAddPileUp;
        Int_t                   fCentralPileupCut;
        Int_t                   fDefaultPileupCut;
        TString                 fCentEstimator;
        UInt_t                  fFilterBit;
        TAxis*                  fPtAxis;                    //
        TAxis*                  fCentAxis;                  //
        Int_t                   NcentBin;
        Int_t                   NPtBin;
        Double_t                fCentMin;
        Double_t                fCentMax;
        double                  fOverrideCentrality;
        Double_t                fPVtxCutX;
        Double_t                fPVtxCutY;
        Double_t                fPVtxCutZ;
        Bool_t                  fRequireHighPtTracks;
        Int_t                   fNHighPtTracks;
        Double_t                fHighPtCut;
        //cuts & selection: tracks
        UInt_t                  fCutChargedTrackFilterBit; // (-) tracks filter bit
        UShort_t                fCutNumTPCclsMin;  // (-) Minimal number of TPC clusters used for track reconstruction
        Double_t                fCutDCAzMax; // (cm) Maximal DCA-z cuts for tracks (pile-up rejection suggested for LHC16)
        Double_t                fCutDCAxyMax; // (cm) Maximal DCA-xy cuts for tracks (pile-up rejection suggested for LHC16)
        Double_t                fChi2Cut;
        Bool_t                  bUseLikeSign;  //Select same charge particle tracks
        Int_t                   iSign;         //+1 or -1
        //cuts & selection: flow
        Double_t                fAbsEtaMax;
        //Double_t                dEtaGap; //Now gotten from CorrTask
        Int_t                   fEtaBinNum;
        Int_t                   fPhiBinNum;
        Bool_t                  fUseWeightsOne;
        int                     fCurrSystFlag;
        Bool_t                  fFillWeights;
        Int_t                   fNumSamples;        //Number of samples for bootstrapping
        //Bool_t                  bHasGap; //Also gotten from CorrTask
        Bool_t                  bDiff;
        Bool_t                  bRef;
        Bool_t                  bPtA;
        Bool_t                  bPtRef;
        Bool_t                  bPtB;
        Double_t                fPOIsPtmax;
        Double_t                fPOIsPtmin;
        Double_t                fRFPsPtMax;
        Double_t                fRFPsPtMin;
        Bool_t                  fRequireTwoPart;
        Bool_t                  bEqualPt;
        Bool_t                  bOnlyPrimariesAndCh;
        //QA
        TH2D*                   fhQAEventsfMult32vsCentr;   //!
        TH2D*                   fhQAEventsMult128vsCentr;   //!
        TH2D*                   fhQAEventsfMultTPCvsTOF;    //!
        TH2D*                   fhQAEventsfMultTPCvsESD;    //!

        ClassDef(AliAnalysisDecorrTask, 1);
};

#endif
