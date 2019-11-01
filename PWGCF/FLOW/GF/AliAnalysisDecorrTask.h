#ifndef AliAnalysisDecorrTask_H
#define AliAnalysisDecorrTask_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliGFWWeights.h"
#include "AliAODTrack.h"
#include "TComplex.h"

class AliAnalysisDecorrTask : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisDecorrTask();
                                AliAnalysisDecorrTask(const char *name);
        virtual                 ~AliAnalysisDecorrTask();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        //event and track selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetRejectAddPileUp(Bool_t use = kTRUE) { fEventRejectAddPileUp = use; }
        void                    SetCentralityEst(TString est){ fCentEstimator = est; }
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        AliEventCuts            fEventCuts;
        void                    SetPtRange(Double_t min, Double_t max) {fPtMin = min; fPtMax = max; }
        void                    SetPOIsPtRange(Double_t min, Double_t max) { fPOIsPtmin = min; fPOIsPtmax = max; }
        void                    SetRPsPtRange(Double_t min, Double_t max) { fRPsPtmin = min; fRPsPtmax = max; }
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }
        void                    SetEtaGap(double etaGap) { dEtaGap = etaGap; }
        void                    SetUseWeights3D(Bool_t useWeights3D) { fUseWeights3D = useWeights3D; }
        void                    HasGap(Bool_t hasGap) { bHasGap = hasGap; }
        void                    SetCentMax(Int_t centmax);
        //Observable selection
        void                    SetDiff(Bool_t diff) { bDiff = diff; }
        void                    SetPtB(Bool_t ptb) { bPtB = ptb; }
        TH2F*                   nuacentral;


    private:
        static const Int_t      fNumHarms = 13;             // maximum harmonics length of flow vector array
        static const Int_t      fNumPowers = 9;             // maximum weight power length of flow vector array
        static const Int_t      fHarmPlots = 3;             // Number of harmonics for plotting (v2, v3, v4)
        static const Int_t      NcentBinMax = 11;           //
        AliAODEvent*            fAOD;                       //! input event
        TList*                  fOutputList;                //! output list
        TList*                  fReferenceFlowList;         //! Output list for reference flow observables
        TList*                  fDifferentialFlowList;      //! Output list for differential flow observables
        TList*                  fDifferentialPtAList;       //! Output list for differential same bin
        TList*                  fPtA_PtB_List;              //! Output list for different pt bins
        TList*                  fObservablesList;           //! List of common observables
        //Weights
        AliGFWWeights*          fWeights;                   //!
        TList*                  fWeightList;                //!
        TH2D*                   fh2Weights;                 //!

        
        //output histograms                
        TH2F*                   fHistPhiEta;                //!
        TH2F*                   fHistPhiPt;                 //!
        TH2F*                   fHistEtaPt;                 //!
        TH1F*                   fHistEta;                   //!
        TH1F*                   fHistPhi;                   //!
        TH1F*                   fHistPt;                    //!
        TProfile*               fHistCent;                  //!
        TH1F*                   fHistVz;                    //!
        
        //Reference flow histograms
        TProfile*               fHistcn[fHarmPlots][5];     //!
        TProfile*               fHistGap[fHarmPlots];       //!
        //Differential flow histograms
        TProfile*               fHistdnNeg[fHarmPlots][NcentBinMax];  //!
        TProfile*               fHistcnPtA[fHarmPlots][NcentBinMax];  //!
        TProfile*               fHistdnPos[fHarmPlots][NcentBinMax];  //!

        //PtA PtB histograms
        TProfile*               fHistPtA_PtB[fHarmPlots][NcentBinMax][22];      //!
        TProfile*               fHistPtA_PtB_LS[fHarmPlots][NcentBinMax][22];   //!
        TProfile*               fHistPtA_PtB_OS[fHarmPlots][NcentBinMax][22];   //!

        //Flow method
        bool                    IsWithinRP(const AliAODTrack* track) const;
        bool                    IsWithinPOI(const AliAODTrack* track) const;
        void                    FillRPvectors(AliAODEvent *fAOD, double dEtaLimit);
        void                    FillPOIvectors(AliAODEvent* fAOD, const double dEtaLimit, const double dPtLow, const double dPtHigh); 
        void                    FillPtBvectors(AliAODEvent* fAOD, const double dEtaLimit, const double dPtLow, const double dPtHigh); 
        void                    CalculateCorrelations(double centrality, double dPt);
        void                    CalculatePtBCorrelations(double centrality, double dPtB, int iPtA);

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
        TComplex FourDiffGap10P(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10M(int n1, int n2, int n3, int n4);
        TComplex FourDiff_PtA_PtA(int n1, int n2, int n3, int n4);
        TComplex FourDiff_PtA_PtB(int n1, int n2, int n3, int n4);
        TComplex FourDiffGap10_PtA_PtB(int n1, int n2, int n3, int n4);
        TComplex Five(int n1, int n2, int n3, int n4, int n5);
        TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6);
        TComplex SixDiff(int n1, int n2, int n3, int n4, int n5, int n6);
        TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
        TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);

    

        AliAnalysisDecorrTask(const AliAnalysisDecorrTask&); // not implemented
        AliAnalysisDecorrTask& operator=(const AliAnalysisDecorrTask&); // not implemented

        //event and track selection
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fEventRejectAddPileUp;
        UInt_t                  fFilterBit;
        Double_t                fPtMin;
        Double_t                fPtMax;
        Double_t                fAbsEtaMax;
        Double_t                dEtaGap;
        Bool_t                  fUseWeights3D;
        TString                 fCentEstimator;
        Bool_t                  InitTask();
        Bool_t                  LoadWeights();
        double                  GetWeights(double dPhi, double dEta, double dVz);
        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;
        Bool_t                  bHasGap;
        Bool_t                  bDiff;
        Bool_t                  bRef;
        Bool_t                  bPtB;
        Double_t                fPOIsPtmax;
        Double_t                fPOIsPtmin;
        Double_t                fRPsPtmax;
        Double_t                fRPsPtmin;
        Int_t                   NcentBin;
        Double_t                centEdges[NcentBinMax+1]; 

        ClassDef(AliAnalysisDecorrTask, 1);
};

#endif
