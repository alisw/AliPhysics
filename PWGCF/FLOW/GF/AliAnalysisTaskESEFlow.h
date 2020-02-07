#ifndef AliAnalysisTaskESEFlow_H
#define AliAnalysisTaskESEFlow_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliGFWWeights.h"
#include "TComplex.h"
#include "TF1.h"
#include "TSpline.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"

#include "AliUniFlowCorrTask.h"

class AliAnalysisTaskESEFlow : public AliAnalysisTaskSE
{
    public:
        enum    ColSystem {kPP = 0, kPPb, kPbPb}; // tag for collisional system
                                AliAnalysisTaskESEFlow();
                                AliAnalysisTaskESEFlow(const char *name, ColSystem colSys, Bool_t bUseV0Calibration = kFALSE);
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
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }
        void                    SetUseWeightsRunByRun(Bool_t bRunByRun) { fFlowRunByRunWeights = bRunByRun; }

        void                    SetVtxZCut(Double_t zCut) { fVtxZCuts = zCut; }
        void                    SetEtaGap(Double_t val) { dEtaGap = val; }
        void                    SetHasEtaGap( Bool_t fEtaGap) { bHasGap = fEtaGap; }
        void                    SetChargedNumTPCclsMin(UShort_t tpcCls) { fCutChargedNumTPCclsMin = tpcCls; }

        void                    SetMakeqSelectionRun(Bool_t actqRun) { fMakeqSelectionRun = actqRun; }
        void                    SetMakeRBRweights(Bool_t actRBRrun) { fMakeRBRweightsRun = actRBRrun; }
        void                    SetUseV0Calibration(Bool_t actV0run) { fV0RunByRunCalibration = actV0run; }

        void                    SetReadMC(Bool_t activate) { fReadMC = activate;}
        void                    SetFlowRFPsPt(Double_t min, Double_t max) { fFlowRFPsPtMin = min; fFlowRFPsPtMax = max; }
        void                    SetFlowPOIsPt(Double_t min, Double_t max) { fFlowPOIsPtMin = min; fFlowPOIsPtMax = max; } 

        void                    SetTPCEse(Bool_t actTPCEse) { fTPCEse = actTPCEse; }
        void                    SetV0CEse(Bool_t actV0CEse) { fV0CEse = actV0CEse; }
        void                    SetV0AEse(Bool_t actV0AEse) { fV0AEse = actV0AEse; }

        void                    AddCorr(std::vector<Int_t> harms, std::vector<Double_t> gaps = std::vector<Double_t>(), Bool_t doRFPs = kTRUE, Bool_t doPOIs = kTRUE) { fVecCorrTask.push_back(new AliUniFlowCorrTask(doRFPs, doPOIs, harms, gaps)); }

        void                    SetWeights(Bool_t kOwn) { bUseOwnWeights = kOwn; }

        void                    SetSampling(Bool_t sample, Int_t iNum) { fSampling = sample; fNumSamples = iNum; }

        void                    SetCentBin(Int_t nbins, Double_t *bins) { fCentAxis->Set(nbins,bins); }
        void                    SetPtBins(Int_t nbins, Double_t *bins) { fPtAxis->Set(nbins, bins); }

        void                    SetSPAnalyzer(Bool_t ActSPAna) { fSPAnalysis = ActSPAna; }
        

    private:
        Bool_t                  fFlowRunByRunWeights;
        Bool_t                  fV0RunByRunCalibration;
        Bool_t                  bUseOwnWeights;
        Double_t                dEtaGap;
        Bool_t                  bHasGap;
        Bool_t                  fSampling;      //Bootstrapping sampling

        static const Int_t      fNumHarms = 13; // maximum harmonics length of flow vector array
        static const Int_t      fNumPowers = 9; // maximum weight power length of flow vector array

        static const Int_t      nCentBinMax = 11;           // maximum number of centrality bins
        static const Int_t      nPtBinMax = 30;             // maximum number of pt bins

        Bool_t                  fInit; // initilization check
        Bool_t                  fMakeqSelectionRun; // make q-selections also used for V0 Calibration runs
        Bool_t                  fMakeRBRweightsRun;

        AliAODEvent*            fAOD;           //!
        TList*                  fOutputList;    //!
        TList*                  fObservables;   //!
        TList*                  fCorrDist;      //!
        TList*                  fpTDiff;        //!
        TList*                  fqnDist;        //!
        TList*                  fpTDiffESETPC;  //!
        TList*                  fcnESETPC;      //!
        TList*                  fpTDiffESEV0C;  //!
        TList*                  fcnESEV0C;      //!
        TList*                  fpTDiffESEV0A;  //!
        TList*                  fcnESEV0A;      //!
        TList*                  SPFlowList;     //!
        TList*                  SPFlowEseList;     //!
        TList*                  fQAEvents;      //!

        TList*                  fFlowWeightsList; //! 
        AliGFWWeights*          fWeights;           //!
        TList*                  fV0CalibList;   //!
        TList*                  fqSelList;   //!
        //output histograms
        TH3F*                   fHistPhiEtaVz;    //!
        TH1F*                   fHistPhi;       //!
        TH1F*                   fHistEta;       //!
        TH1F*                   fHistPt;        //!
        TH1F*                   fHistZVertex;   //!

        TSpline3*               fSplq2TPC[90];  // q2 TPC splines
        TSpline3*               fSplq3TPC[90];  // q3 TPC splines
        TSpline3*               fSplq2V0C[90];  // q2 V0C splines
        TSpline3*               fSplq3V0C[90];  // q3 V0C splines
        TSpline3*               fSplq2V0A[90];  // q2 V0A splines
        TSpline3*               fSplq3V0A[90];  // q3 V0A splines

        TH3F*                   fh3Weights; //!
        TH1F*                   fhV0Calib;  //!
        TH1F*                   fHistPDG; //!
        

        /////////////////////////// CALIBRATION QA HISTOGRAMS ////////////////////////////////////
        TH2D*                   fq2TPC;    //!
        TH2D*                   fq3TPC;    //!
        TH2D*                   fq2V0C;    //!
        TH2D*                   fq3V0C;    //!
        TH2D*                   fq2V0A;    //!
        TH2D*                   fq3V0A;    //!

        TH2F*                   fQnxV0C[2];    //!
        TH2F*                   fQnyV0C[2];    //!
        TH2F*                   fQnxV0A[2];    //!
        TH2F*                   fQnyV0A[2];    //!
        TH2F*                   fQnxTPC[2];    //!
        TH2F*                   fQnyTPC[2];    //!

        
        TH1F*                   fQnxV0Cm[2];    //!
        TH1F*                   fQnyV0Cm[2];    //!
        TH1F*                   fQnxV0Am[2];    //!
        TH1F*                   fQnyV0Am[2];    //!
        TH1F*                   fQnxTPCm[2];    //!
        TH1F*                   fQnyTPCm[2];    //!

        TH1F*                   fQnxV0Cs[2];    //!
        TH1F*                   fQnyV0Cs[2];    //!
        TH1F*                   fQnxV0As[2];    //!
        TH1F*                   fQnyV0As[2];    //!

        TH2F*                   fQnxV0CEse[2];    //!
        TH2F*                   fQnyV0CEse[2];    //!
        TH2F*                   fQnxV0AEse[2];    //!
        TH2F*                   fQnyV0AEse[2];    //!
        TH2F*                   fQnxTPCEse[2];    //!
        TH2F*                   fQnyTPCEse[2];    //!

        TH2F*                   fQnxV0CCor[2];    //!
        TH2F*                   fQnyV0CCor[2];    //!
        TH2F*                   fQnxV0ACor[2];    //!
        TH2F*                   fQnyV0ACor[2];    //! 
        ////////////////////////// end /////////////////////////////////////////

        //// SCALAR-PRODUCT UNIT VECTOR FLOW /////

        //Event-plane nHarm=2
        TH1F*                   fhEvPlPsi_2V0C;  //!
        TH1F*                   fhEvPlPsi_2V0A;  //!


        TProfile*               fProfNPar; //!
        TH2F*                   fhV0Multiplicity;    //!
        TH2F*                   fhV0CorrMult;       //!
        TH2F*                   fhqnTPCvqnV0C[2];  //!
        TH2F*                   fhqnV0CvqnV0A[2];  //!
        TH2F*                   fhqnTPCvqnV0A[2];  //!

        void CorrelationTask(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz, Int_t fSpCent);
        void SPVienna(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz, Int_t q2ESECodeV0C);
        void SPFillPOI(const Double_t dPtL, const Double_t dPtLow, const Double_t dPtHigh, const float dVz, Int_t iTracks,const Int_t CenterCode);
        void FillCorrelation(const AliUniFlowCorrTask* const task, const Float_t centrality, const Double_t dPt, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C, Int_t q2ESECodeV0A, Int_t q3ESECodeV0A, Bool_t doRef, Bool_t doDiff);
        void FillObsDistributions(const Int_t iTracks, const AliAODEvent* fAOD, const float dVz, const Float_t centrality);
        // Calculate flow vectors for reference and POIs
        void RFPVectors(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz);
        void POIVectors(const Int_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz, Int_t q2ESECodeTPC, Int_t q3ESECodeTPC, Int_t q2ESECodeV0C, Int_t q3ESECodeV0C,Int_t q2ESECodeV0A, Int_t q3ESECodeV0A);
        void ReducedqVectorsTPC(const Float_t centrality, const Int_t iTracks, const AliAODEvent* fAOD, const float dVz, const Int_t SPCode);
        void ReducedqVectorsV0(const Float_t centrality, const AliAODEvent* fAOD, const Int_t SPCode);
        
        void FillqnRedTPC(const Float_t centrality);
        void FillqnRedV0(const Float_t centrality, TString V0type);
        void FillPOI(const Double_t dPtL, const Double_t dPtLow, const Double_t dPtHigh, const float dVz, Int_t iTracks);

        Int_t GetSamplingIndex() const;
        
        Int_t GetCentralityCode(const Float_t centrality);
        Int_t GetEsePercentileCode(Double_t qPerc) const;
        Bool_t WithinRFP(const AliVParticle* track) const;
        Bool_t WithinPOI(const AliVParticle* track) const;
        Bool_t LoadqSelection();
        Bool_t ProcessMCParticles();

        Bool_t InitializeTask();
        Bool_t LoadWeights(); // load weights histograms
        Bool_t LoadV0Calibration();
        Double_t GetFlowWeight(const AliAODTrack* track, const float dVz) const;
        //############ GENERIC FRAMEWORK ############# MODIFIED WITH ESE //

        double GetWeight(double phi, double eta, double vz,  double runNumber);
        double GetPtWeight(double pt, double eta, float vz,  double runNumber);
        void ResetFlowVector(TComplex (&array)[fNumHarms][fNumPowers]); // set values to TComplex(0,0,0) for given array
        void ResetReducedqVector(Double_t (&array)[2]);

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

        Double_t qnTPC[2];
        Double_t QxnTPC[2];
        Double_t QynTPC[2];

        Double_t qnV0C[2];
        Double_t QxnV0C[2];
        Double_t QynV0C[2];

        Double_t qnV0A[2];
        Double_t QxnV0A[2];
        Double_t QynV0A[2];

        Double_t QxnTPCEse[2]; // after ese correction
        Double_t QynTPCEse[2];

        Double_t QxnV0CEse[2];
        Double_t QynV0CEse[2];

        Double_t QxnV0AEse[2];
        Double_t QynV0AEse[2];

        Double_t QxnV0CCorr[2];
        Double_t QynV0CCorr[2];

        Double_t QxnV0ACorr[2];
        Double_t QynV0ACorr[2];

        Double_t vnSPV0A;
        Double_t vnSPV0C;
        

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
        TComplex SixGap10(int n1, int n2, int n3, int n4, int n5, int n6);
        TComplex ThreePos(int n1, int n2, int n3);
        TComplex ThreeNeg(int n1, int n2, int n3);
        TComplex SixDiffGap10M(int n1, int n2, int n3, int n4, int n5, int n6);
        TComplex ThreeDiffNeg(int n1, int n2, int n3);
        TComplex ThreeDiffPos(int n1, int n2, int n3);

        //############# END #############//

        AliAnalysisTaskESEFlow(const AliAnalysisTaskESEFlow&); // not implemented
        AliAnalysisTaskESEFlow& operator=(const AliAnalysisTaskESEFlow&); // not implemented

        //event and track selection
        ColSystem               fColSystem; // collisional system
        AliVEvent::EOfflineTriggerTypes    fTrigger;
        Bool_t                  fEventRejectAddPileUp;
        UInt_t                  fFilterBit;
        Double_t                fAbsEtaMax;
        Double_t                fVtxZCuts;
        TString                 fCentEstimator;
        UShort_t                fCutChargedNumTPCclsMin;
        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;

        Bool_t                  fReadMC;
        AliMCEvent*             fMCEvent;       //! corresponding MC event
        Double_t                fFlowRFPsPtMin; // [0.2] (GeV/c) min pT treshold for RFPs particle for reference flow
        Double_t                fFlowRFPsPtMax; // [5.0] (GeV/c) max pT treshold for RFPs particle for reference flow
        Double_t                fFlowPOIsPtMin; // [0] (GeV/c) min pT treshold for POIs for differential flow
        Double_t                fFlowPOIsPtMax; // [10] (GeV/c) max pT treshold for POIs for differential flow

        TAxis*                  fPtAxis;
        TAxis*                  fCentAxis;
        Int_t                   nCentBin;
        Int_t                   nPtBin;
        Double_t                CentEdges[nCentBinMax+1];
        Double_t                PtEdges[nPtBinMax+1];

    
        Bool_t                  fTPCEse;
        Bool_t                  fV0CEse;
        Bool_t                  fV0AEse;

        Int_t                   fIndexSampling;
        Int_t                   fNumSamples;

        Bool_t                  fSPAnalysis;

        std::vector<AliUniFlowCorrTask*>    fVecCorrTask;


        ClassDef(AliAnalysisTaskESEFlow, 1);
};

#endif
