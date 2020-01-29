#ifndef AliAnalysisDecorrTask_H
#define AliAnalysisDecorrTask_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliGFWWeights.h"
#include "AliAODTrack.h"
#include "TComplex.h"
#include "TAxis.h"
#include "AliUniFlowCorrTask.h"

class Taxis;
class AliVEvent;
class TProfile;
class AliUniFlowCorrTask;

class AliAnalysisDecorrTask : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisDecorrTask();
                                AliAnalysisDecorrTask(const char *name);
        virtual                 ~AliAnalysisDecorrTask();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        //Analysis setters
        void                    SetSampling(Bool_t sample, Int_t iNum) { fSampling = sample; fNumSamples = iNum; }      //Use jack-knife resampling
        void                    SetFillQA(Bool_t fill = kTRUE) { fFillQA = fill; }
        void                    SetSmallSystem(Bool_t small = kTRUE) { fSmallSystem = kTRUE; }
        //event selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetRejectAddPileUp(Bool_t use = kTRUE) { fEventRejectAddPileUp = use; }
        void                    SetCentralityEst(TString est){ fCentEstimator = est; }
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        void                    SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
        void                    SetCentBin(Int_t nbins, Double_t *bins) { fCentAxis->Set(nbins,bins); }
        void                    SetCentLim(Double_t min, Double_t max) { fCentMin = min; fCentMax = max; } //Not used yet
        void                    SetPtBins(Int_t nbins, Double_t *bins) { fPtAxis->Set(nbins, bins); }
        AliEventCuts            fEventCuts;
        //track selection
        void                    SetDCAzMax(Double_t dcaz) {  fCutDCAzMax = dcaz; }
        void                    SetDCAxyMax(Double_t dcaxy) {  fCutDCAxyMax = dcaxy; }
        void                    SetNumTPCclsMin(UShort_t tpcCls) { fCutNumTPCclsMin = tpcCls; }
        void                    SetUseLikeSign(Bool_t use, Int_t sign) { bUseLikeSign = use; iSign = sign; }
        void                    SetChargedTrackFilterBit(UInt_t filter) { fCutChargedTrackFilterBit = filter; } //Not implemented
        //Flow selection
        void                    AddCorr(std::vector<Int_t> harms, std::vector<Double_t> gaps = std::vector<Double_t>(), Bool_t doRFPs = kTRUE, Bool_t doPOIs = kTRUE) { fVecCorrTask.push_back(new AliUniFlowCorrTask(doRFPs, doPOIs, harms, gaps)); }
        void                    SetPOIsPt(Double_t min, Double_t max) { fPOIsPtmin = min; fPOIsPtmax = max; }
        void                    SetRFPsPt(Double_t min, Double_t max) { fRFPsPtMin = min; fRFPsPtMax = max; }
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }
        void                    SetEtaBins(Int_t bins) { fEtaBinNum = bins; }
        void                    SetPhiBins(Int_t bins) { fPhiBinNum = bins; }
        void                    SetEtaGap(double etaGap) { dEtaGap = etaGap; }
        void                    SetUseWeights3D(Bool_t use) { fUseWeights3D = use; }    //Use 3D weights (phi, eta Vz)
        void                    SetFillWeights(Bool_t fill) { fFillWeights = fill; }    //Only fill histograms for weights calculations
        Bool_t                  GetUseWeights3D() { return fUseWeights3D; }             //Check if 3D weights are used for macro path to weights
        void                    HasGap(Bool_t hasGap) { bHasGap = hasGap; } 
        void                    CalculateHigherOrderVn(Bool_t calc) { bHigherOrder = calc; }   //Calculate higher order particle correlation differential vn with jack-knife resampling

        //Observable selection
        void                    DoRFPs(Bool_t ref) { bRef = ref; }              //Calculate integrad flow
        void                    DoDiff(Bool_t diff) { bDiff = diff; }           //Calculate pt differential flow
        void                    DoPtB(Bool_t ptb) { bPtB = ptb; }               //Calculate flow with particles from different pt bins
        void                    DoSC(Bool_t integrated, Bool_t singlediff) { fInt = integrated, fSingle = singlediff; }

    
    private:
        static const Int_t      fNumHarms = 13;             // maximum harmonics length of flow vector array
        static const Int_t      fNumPowers = 9;             // maximum weight power length of flow vector array
        static const Int_t      NcentBinMax = 11;           //
        static const Int_t      NPtBinMax = 30;             //

        TList*                  fFlowList;                //! output list
        TList*                  fFlowWeights;             //! 
        TList*                  fQA;                      //!


        Bool_t                  InitTask();
        Bool_t                  LoadWeights();
        double                  GetWeights(double dPhi, double dEta, double dVz);
        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;
        Int_t                   GetSamplingIndex() const;
        //Weights
        AliGFWWeights*          fWeights;                   //!
        TList*                  fWeightList;                //!
        TH2D*                   fh2Weights;                 //!
        TH3D*                   fh3Weights;                 //!
        void                    FillWeights();
        
        //Flow methods
        bool                    IsWithinRP(const AliAODTrack* track) const;
        bool                    IsWithinPOI(const AliAODTrack* track) const;
        void                    FillRPvectors(double dEtaLimit);
        void                    FillPOIvectors(const double dEtaLimit, const double dPtLow, const double dPtHigh); 
        void                    FillPtBvectors(const double dEtaLimit, const double dPtLow, const double dPtHigh); 
        void                    CalculateCorrelations(double centrality, double dPtA, double dPtB, Bool_t doRef, Bool_t doDiff, Bool_t doPtB);

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
        TComplex FourDiffGap10P(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10M(int n1, int n2, int n3, int n4);
        TComplex FourDiff_PtA_PtA(int n1, int n2, int n3, int n4);
        TComplex FourDiff_PtA_PtB(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10_PtA_PtB(int n1, int n2, int n3, int n4);
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
        AliAODEvent*            fAOD;                       //! input event
        Bool_t                  fInitTask;                  //Initialization
        
        std::vector<AliUniFlowCorrTask*>    fVecCorrTask;   //
        
        //cuts & selection: Analysis
        Bool_t                  fSampling;      //Bootstrapping sampling
        Bool_t                  fFillQA;        //Fill QA histograms
        Bool_t                  fSmallSystem;   //Analyse small system
        //cuts & selection: events
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fEventRejectAddPileUp;
        TString                 fCentEstimator;
        UInt_t                  fFilterBit;
        TAxis*                  fPtAxis;                    //
        TAxis*                  fCentAxis;                  //
        Int_t                   NcentBin;
        Int_t                   NPtBin;
        Double_t                centEdges[NcentBinMax+1];
        Double_t                PtEdges[NPtBinMax+1];
        Double_t                fCentMin;
        Double_t                fCentMax;
        Double_t                fPVtxCutZ;
        //cuts & selection: tracks
        UInt_t                  fCutChargedTrackFilterBit; // (-) tracks filter bit
        UShort_t                fCutNumTPCclsMin;  // (-) Minimal number of TPC clusters used for track reconstruction
        Double_t                fCutDCAzMax; // (cm) Maximal DCA-z cuts for tracks (pile-up rejection suggested for LHC16)
        Double_t                fCutDCAxyMax; // (cm) Maximal DCA-xy cuts for tracks (pile-up rejection suggested for LHC16)
        Bool_t                  bUseLikeSign;  //Select same charge particle tracks
        Int_t                   iSign;         //+1 or -1
        //cuts & selection: flow
        Double_t                fAbsEtaMax;
        Double_t                dEtaGap;
        Int_t                   fEtaBinNum;
        Int_t                   fPhiBinNum;
        Bool_t                  fUseWeights3D;
        Bool_t                  fFillWeights;
        Int_t                   fNumSamples;        //Number of samples for bootstrapping
        Bool_t                  bHasGap;
        Bool_t                  bDiff;
        Bool_t                  bRef;
        Bool_t                  bPtB;
        Bool_t                  fInt;
        Bool_t                  fSingle;
        Bool_t                  bHigherOrder;
        Double_t                fPOIsPtmax;
        Double_t                fPOIsPtmin;
        Double_t                fRFPsPtMax;
        Double_t                fRFPsPtMin;

        ClassDef(AliAnalysisDecorrTask, 1);
};

#endif
