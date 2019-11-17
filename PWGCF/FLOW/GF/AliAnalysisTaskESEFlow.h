#ifndef AliAnalysisTaskESEFlow_H
#define AliAnalysisTaskESEFlow_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliGFWWeights.h"

class AliAnalysisTaskESEFlow : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskESEFlow();
                                AliAnalysisTaskESEFlow(const char *name);
        virtual                 ~AliAnalysisTaskESEFlow();

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
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }
        void                    SetUseWeightsRunByRun(Bool_t bRunByRun) { fFlowRunByRunWeights = bRunByRun; }

        void                    SetEtaGap(Double_t val) { dGap = val; }

        void                    SetqSelectionRun(Bool_t actqRun = kFALSE) { fqRun = actqRun; }

        void                    SetReadMC(Bool_t activate = kFALSE) { fReadMC = activate;}
        void                    SetFlowRFPsPt(Double_t min, Double_t max) { fFlowRFPsPtMin = min; fFlowRFPsPtMax = max; }
        void                    SetFlowPOIsPt(Double_t min, Double_t max) { fFlowPOIsPtMin = min; fFlowPOIsPtMax = max; } 

        void                    SetWeights(Bool_t kOwn) { bUseOwnWeights = kOwn; }
        

        //runAnalysis inputs
        Bool_t                  fFlowRunByRunWeights;
        Bool_t                  bUseOwnWeights;

        Double_t                dGap;

    private:

    static const Int_t      fNumHarms = 10; // maximum harmonics length of flow vector array
    static const Int_t      fNumPowers = 10; // maximum weight power length of flow vector array
    static const Int_t      fNumHarmHists = 3; // how many harmonics hists
    static const Int_t      fNumCentHists = 7; // how many cent hists should there be
    static const Int_t      fnqCuts = 11;        // number of q selection cuts from 0-10 in 10% intervals --- probably not used anywhere?

    Bool_t                  fInit; // ini check
    Bool_t                  fqRun;

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TList*                  fObservables;   //!
        TList*                  fCorrDist;      //!
        TList*                  fpTDiff;        //!
        TList*                  fqnDist;        //!
        TList*                  fpTDiffqselec;  //!
        TList*                  fcnqselec;      //!

        //list in lists WIP

        TList*                  fFlowWeightsList; //! list of weights
        TList*                  fqCutsList;     //! list of q selection cuts
        TH1F*                   fHistqCuts[3][fnqCuts];     //!
        AliGFWWeights*          fWeights;           //!
        //output histograms
        TH2F*                   fHistPhiEta;    //!
        TH1F*                   fHistPhi;       //!
        TH1F*                   fHistEta;       //!
        TH1F*                   fHistPt;        //!
        TH1F*                   fHistZVertex;   //!

        TProfile*               fProfcn_2gap[fNumHarmHists]; //!

        TProfile*               fProfPTdn_2gap[fNumHarmHists][fNumCentHists];    //!
        TH1F*                   fHistqn_reduced[fNumHarmHists][fNumCentHists];     //!
        TProfile*               fProfPTdn_2gap_large[fNumHarmHists][fNumCentHists];   //!
        TProfile*               fProfPTdn_2gap_small[fNumHarmHists][fNumCentHists];   //!
        TProfile*               fProfcn_2gap_largeqn[fNumHarmHists][fNumCentHists]; //! used for normalizing
        TProfile*               fProfcn_2gap_smallqn[fNumHarmHists][fNumCentHists]; //! used for normalizing

        TProfile*               fProfPTdn_2gap_B[fNumHarmHists][fNumCentHists];    //!

        TH2D*                   fh2Weights; //! contains phi eta weights

        TH1F*                   fHistPDG; //! histogram of pdg codes for MC

        /////////////////// Work in progress //////////////////////////
        TProfile*               fProfcn_2gap_qn[3][fNumHarmHists][fNumCentHists][10]; //!        
        ///////////////////////////////////////////////////////////////

        TH1F*                   fHistq2_red_cent_30_31;                     //!
        TH1F*                   fHistqn_red_cent_0_1[2];                     //!

        // fill dphi/deta/dpt histograms for weights
        void FillObsDistributions(const Int_t iTracks, const AliAODEvent* fAOD);
        // Calculate flow vectors for reference and POIs
        void RFPVectors(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD);
        void POIVectors(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD);
        void ReducedRFPVectors(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD);
        void FillRFP(const Float_t centrality, const int nHarm);
        void Filldn_2(const Float_t centrality, const double dPt, const int nHarm);
        void Fillqnreduced(const Float_t centrality);
        void FillPOI(const Double_t dPtL, const Double_t dPtLow, const Double_t dPtHigh, const float dVz, const Int_t iTracks);
        void FillIntegratecqnCut(const Float_t centrality, const int nHarm, const double c, const double c_weight, const int q_i);
        Bool_t WithinRFP(const AliVParticle* track) const;
        Bool_t WithinPOI(const AliVParticle* track) const;
        Bool_t ProcessMCParticles();

        Bool_t InitializeTask();
        Bool_t LoadWeights(); // load weights histograms
        Double_t GetFlowWeight(const AliAODTrack* track, const float dVz) const;
        Bool_t LoadqCuts();
        Double_t GetqSelectionCut(Int_t nHarm, Int_t CentRange, Int_t Entry);
        //############ GENERIC FRAMEWORK ############# MODIFIED WITH ESE //

        double GetWeight(double phi, double eta, double vz,  double runNumber);
        double GetPtWeight(double pt, double eta, float vz,  double runNumber);
        void ResetFlowVector(TComplex (&array)[fNumHarms][fNumPowers]); // set values to TComplex(0,0,0) for given array
        void ResetReducedFlowVector(Double_t (&array)[fNumHarms][fNumPowers]);

        Bool_t sortPt(const AliAODTrack* t1, const AliAODTrack* t2) { return (t1->Pt() < t2->Pt()); } // function for std::sort
        
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
        Double_t qvector_red[fNumHarms][fNumPowers];
        Double_t sumCos[fNumHarms][fNumPowers];
        Double_t sumSin[fNumHarms][fNumPowers];

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

        //############# END #############//

        AliAnalysisTaskESEFlow(const AliAnalysisTaskESEFlow&); // not implemented
        AliAnalysisTaskESEFlow& operator=(const AliAnalysisTaskESEFlow&); // not implemented

        //event and track selection
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fEventRejectAddPileUp;
        UInt_t                  fFilterBit;
        Double_t                fPtMin;
        Double_t                fPtMax;
        Double_t                fAbsEtaMax;
        TString                 fCentEstimator;
        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;
        Bool_t                  fReadMC;
        AliMCEvent*             fMCEvent;       //! corresponding MC event
        int                     N_counter;
        Float_t                 fCentInterval[12];
        Double_t                fFlowRFPsPtMin; // [0.2] (GeV/c) min pT treshold for RFPs particle for reference flow
        Double_t                fFlowRFPsPtMax; // [5.0] (GeV/c) max pT treshold for RFPs particle for reference flow
        Double_t                fFlowPOIsPtMin; // [0] (GeV/c) min pT treshold for POIs for differential flow
        Double_t                fFlowPOIsPtMax; // [10] (GeV/c) max pT treshold for POIs for differential flow


        ClassDef(AliAnalysisTaskESEFlow, 1);
};

#endif
