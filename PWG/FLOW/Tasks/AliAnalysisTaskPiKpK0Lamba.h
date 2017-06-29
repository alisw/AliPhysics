#ifndef ALIANALYSISTASKDATA_H
#define ALIANALYSISTASKDATA_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TString.h>
#include <THnSparse.h>

// AliRoot includes
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliPIDResponse.h"


class AliAnalysisTaskPiKpK0Lamba : public AliAnalysisTaskSE {
    public:
        AliAnalysisTaskPiKpK0Lamba();
        AliAnalysisTaskPiKpK0Lamba(const char *name);

        virtual ~AliAnalysisTaskPiKpK0Lamba();

        virtual void   UserCreateOutputObjects();
        virtual void   UserExec(Option_t *option);
        virtual void   Terminate(Option_t *);


        virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
        virtual void  SetFilterbit(UInt_t filterbit){fFilterbit = filterbit;}
        virtual void  SetNoClus(Int_t noclus){fNoClus = noclus;}
        virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
        virtual void  SetMinPt(Double_t minPt){fMinPt = minPt;}
        virtual void  SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
        virtual void  SetFlagPileUp(Bool_t IsPileUP){fPileUp = IsPileUP;}
        virtual void  SetFlagPileUpTOF(Bool_t IsPileUPTOF){fPileUpTOF = IsPileUPTOF;}
        virtual void  SetCentFlag(Short_t nCent){fCent = nCent;}
        virtual void  SetNHarmonic(Double_t nHarm){fNHarm = nHarm;}
        virtual void  SetNsigCut(Double_t nSigCut){fNsigCut = nSigCut;}
        virtual void  SetFlagQA(Bool_t isQA){fQA = isQA;}
        virtual void  SetFlagPhiCut(Bool_t isPhiCut){fPhiCut = isPhiCut;}
        virtual void  SetMinPiPtHCut(Double_t minPi){fMinPiCut = minPi;}
        virtual void  SetMaxPiPtHCut(Double_t maxPi){fMaxPiCut = maxPi;}
        virtual void  SetMinPPtHCut(Double_t minP){fMinPCut = minP;}
        virtual void  SetMaxPPtHCut(Double_t maxP){fMaxPCut = maxP;}
        virtual void  SetFlagQPos(Bool_t isQpos){fQPos = isQpos;}
        virtual void  SetFlagQNeg(Bool_t isQneg){fQNeg = isQneg;}
        virtual void  SetFlagExclPID(Bool_t isExcl){fExclPID = isExcl;}
        virtual void  SetNoClusPID(Int_t nocluspid){fNoClusPid = nocluspid;}
        virtual void  SetFlagQAV0(Bool_t isQAV0){fQAV0 = isQAV0;}
        virtual void  SetFlagQAOutl(Bool_t isQAOutl){fQAOutl = isQAOutl;}
        virtual void  SetFlagBin1Cent(Bool_t isBin1Cent){fBin1Cent = isBin1Cent;}
        virtual void  SetFlagRemoveCh46V0A(Bool_t isRemChV0A){fRemChV0A = isRemChV0A;}
        virtual void  SetFlagEtaRange(Short_t etaRange){fEtaRange = etaRange;}
        virtual void  SetRemovePhiReg(Bool_t remPhiReg){fRemPhiReg = remPhiReg;}
        virtual void  SetCutMultESDdif(Float_t cutMultESDdif){fCutMultESDdif = cutMultESDdif;}
        virtual void  SetCrsRowsFrcShClsCuts(Bool_t crsRowsFrcShCls){fCrsRowsFrcShCls = crsRowsFrcShCls;}
        virtual void  SetFlagPsi42A(Bool_t flagPsi42A){flagPsi42A = flagPsi42A;}
        virtual void  SetNPtBins(Int_t nPtB){fNPtBins = nPtB;}
        virtual void  SetNcrFind(Float_t ncrFind){fNcrFind = ncrFind;}
        virtual void  SetDCADghtPV(Float_t dcaDghtPV){fDCADghtPV = dcaDghtPV;}
        virtual void  SetMaxDCADght(Float_t maxDCADght){fMaxDCADght = maxDCADght;}
        virtual void  SetCosPA(Float_t cosPA){fCosPA = cosPA;}
        virtual void  SetMinRad(Float_t minRad){fMinRad = minRad;}
        virtual void  SetMaxRad(Float_t maxRad){fMaxRad = maxRad;}
        virtual void  SetFlagArmPodCut(Bool_t isArmPod){fArmPodCut = isArmPod;}
        virtual void  SetFlagMinPtCutDght(Bool_t isMinPtDght){fMinPtDght = isMinPtDght;}
        virtual void  SetPtBins(Double_t ptBins[18]) {for (Int_t i = 0; i < 18; i++) fPtBins[i] = ptBins[i];}


    private:
        virtual Float_t GetVertex(AliAODEvent* aod) const;
        virtual void Analyze(AliAODEvent* aod, Float_t vtxZ);
        void OpenInfoCalbration(Int_t run);
        Double_t GetWDist(const AliVVertex* v0, const AliVVertex* v1);
        Bool_t plpMV(const AliVEvent *event);
        Short_t FindMinNSigma(Double_t nSpi, Double_t nSk, Double_t nSp);
        Bool_t GetDoubleCountingPi(Double_t nSpi, Short_t minNSigma);
        Bool_t GetDoubleCountingK(Double_t nSk, Short_t minNSigma);
        Bool_t GetDoubleCountingP(Double_t nSp, Short_t minNSigma);
        Double_t GetRapidity(Double_t mass, Double_t Pt, Double_t Eta);
        Int_t GetPtBin(Double_t valPt) const;


        AliAODEvent* fAOD;                //! AOD object
        AliPIDResponse*  fPIDResponse;    //! pid response object

        Int_t        fRun;                // run number - for calibration
        TH1D*        fMultV0;             // profile from V0 multiplicity

        TH1D*        fQxnmV0A;            // <Qx2> V0A
        TH1D*        fQynmV0A;            // <Qy2> V0A
        TH1D*        fQxnsV0A;            // sigma Qx2 V0A
        TH1D*        fQynsV0A;            // sigma Qy2 V0A
        TH1D*        fQxnmV0C;            // <Qx2> V0C
        TH1D*        fQynmV0C;            // <Qy2> V0C
        TH1D*        fQxnsV0C;            // sigma Qx2 V0C
        TH1D*        fQynsV0C;            // sigma Qy2 V0C


        //
        // Cuts and options
        //
        TF1*         fLowCut;             // cut low for centrality outliers
        TF1*         fHighCut;            // cut high for centrality outliers
        Double_t     fVtxCut;             // Vtx cut on z position in cm
        UInt_t       fFilterbit;          // filter bit
        Double_t     fEtaCut;             // Eta cut used to select particles
        Int_t        fNoClus;	          // No of TPC clusters
        Double_t     fMinPt;              // Min pt - for histogram limits
        Double_t     fMaxPt;              // Max pt - for histogram limits
        Bool_t       fPileUp;             // flag for pileup
        Bool_t       fPileUpTOF;          // flag for pileup
        Short_t      fCent;               // centrality flag
        Double_t     fNHarm;              // harmonic number
        Double_t     fNsigCut;            // combined sigma cut value
        Bool_t       fQA;                 // qa flag
        Bool_t       fPhiCut;             // flag to remove tracks close to the TPC sectors
        TF1*         fPhiCutLow;          // phi cut low to remove tracks close to the TPC sectors
        TF1*         fPhiCutHigh;         // phi cut high to remove tracks close to the TPC sectors
        TF1*         fMultTOFLowCut;      // cut low for TOF multiplicity outliers
        TF1*         fMultTOFHighCut;     // cut high for TOF multiplicity outliers
        Double_t     fMinPiCut;           // min cut for Pi TPC
        Double_t     fMaxPiCut;           // max cut for Pi TPC
        Double_t     fMinPCut;            // min cut for P TPC
        Double_t     fMaxPCut;            // max cut for P TPC
        Bool_t       fQPos;               // flag for positive particles
        Bool_t       fQNeg;               // flag for negative particles
        Bool_t       fExclPID;            // flag for exclusive PID
        Int_t        fNoClusPid;          // no of TPC clusters used for PID
        Bool_t       fQAV0;               // qa V0 flag
        TF1*         fMultCentLowCut;     // cut low for multiplicity centrality outliers
        Bool_t       fQAOutl;             // flag for QA outliers
        Bool_t       fBin1Cent;           // flag to run for bin 1 centrality
        Bool_t       fRemChV0A;           // flag to remove channel 46 from V0A
        Short_t      fEtaRange;           // switch to use only pos or neg eta tracks
        Bool_t       fRemPhiReg;          // remove 1.7<phi<2.7 region
        Float_t      fCutMultESDdif;      // cut for multESDdif pileup
        Bool_t       fCrsRowsFrcShCls;    // cuts for new cross rows and fraction of shared clusters
        Bool_t       flagPsi42A;          // flag for Psi42 analysis
        Int_t        fNPtBins;            // number of pt bins
        Float_t      fNcrFind;            // number of cross rows over findable clusters
        Float_t      fDCADghtPV;          // DCA daughters to primary vertex
        Float_t      fMaxDCADght;         // DCA daughters
        Float_t      fCosPA;              // cos pointing angle
        Float_t      fMinRad;             // V0 radius cut low
        Float_t      fMaxRad;             // V0 radius cut low
        Bool_t       fArmPodCut;          // flag for Arm-Pod cut
        Bool_t       fMinPtDght;          // flag for min pT cut for daughters
        Double_t     fPtBins[18];         // pt bins edges




        //
        // Output objects
        //
        TList*        fListOfObjects;     //! Output list of objects

        TH1I*         fVtx;               //! Event vertex info
        TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
        TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts

        TH2F*         fMultvsCentBef;        //! Multiplicity vs centrality
        TH2F*         fMultvsCentAft;        //! Multiplicity vs centrality

        TH1F*         fCentrBef;             //! centrality distribution
        TH1F*         fCentrAft;             //! centrality distribution

        TH2F*         fCenCL0vsV0MBef;       //! CL0-V0M correlation before
        TH2F*         fCenCL1vsV0MBef;       //! CL0-V0M correlation before
        TH2F*         fCenCL0vsCL1Bef;       //! CL0-CL1 correlation before
        TH2F*         fCenCL0vsV0MAft;       //! CL0-V0M correlation before
        TH2F*         fCenCL1vsV0MAft;       //! CL0-V0M correlation before
        TH2F*         fCenCL0vsCL1Aft;       //! CL0-CL1 correlation before


        TH2I* fMultFBvsMultFBTOFBef;         //! multiplicity FB vs FB+TOF before
        TH2I* fMultFBvsMultFBTOFAft;         //! multiplicity FB vs FB+TOF after

        TH2D* fMultESDDifvsMultTPCBef;       //! multiplicity dif ESD-3.38*TPC vs TPC before
        TH2D* fMultESDDifvsMultTPCAft;       //! multiplicity dif ESD-3.38*TPC vs TPC after

        TH2I* fSPclsvsSPDtrksBef;            //! SPD clusters vs tracklets before
        TH2I* fSPclsvsSPDtrksAft;            //! SPD clusters vs tracklets after

        TH2F* fMultV0vsMultTPCoutBef;        //! multV0-multTPCout before
        TH2F* fMultV0vsMultTPCoutAft;        //! multV0-multTPCout before

        THnSparseF*    fPidQA;             //! qa histo for tracks
        THnSparseF*    fAllQA;             //! qa histo for tracks
        THnSparseF*    fV0QA;              //! qa histo for V0

        THnSparseF*    fPidQAB1C;             //! qa histo for tracks
        THnSparseF*    fAllQAB1C;             //! qa histo for tracks

        TH2D*         fQxavsV0Bef;          //! corrected Qx V0A cent
        TH2D*         fQyavsV0Bef;          //! corrected Qy V0A cent
        TH2D*         fQxcvsV0Bef;          //! corrected Qx V0C cent
        TH2D*         fQycvsV0Bef;          //! corrected Qx V0C cent

        TH2D*         fQxavsVtxZBef;          //! corrected Qx V0A vtx
        TH2D*         fQyavsVtxZBef;          //! corrected Qy V0A vtx
        TH2D*         fQxcvsVtxZBef;          //! corrected Qx V0C vtx
        TH2D*         fQycvsVtxZBef;          //! corrected Qx V0C vtx


        TH2D*         fQxavsV0Aft;          //! corrected Qx V0A cent
        TH2D*         fQyavsV0Aft;          //! corrected Qy V0A cent
        TH2D*         fQxcvsV0Aft;          //! corrected Qx V0C cent
        TH2D*         fQycvsV0Aft;          //! corrected Qx V0C cent

        TH2D*         fQxavsVtxZAft;          //! corrected Qx V0A vtx
        TH2D*         fQyavsVtxZAft;          //! corrected Qy V0A vtx
        TH2D*         fQxcvsVtxZAft;          //! corrected Qx V0C vtx
        TH2D*         fQycvsVtxZAft;          //! corrected Qx V0C vtx


        TProfile*     fV0AV0Cvn;         //! correlation V0A-V0C for resolution
        TProfile*     fV0ATPCvn;         //! correlation V0A-TPC for resolution
        TProfile*     fV0CTPCvn;         //! correlation V0C-TPC for resolution

        TProfile*     fV0AV0CvnB1C;         //! correlation V0A-V0C for resolution
        TProfile*     fV0ATPCvnB1C;         //! correlation V0A-TPC for resolution
        TProfile*     fV0CTPCvnB1C;         //! correlation V0C-TPC for resolution

        TProfile*     fV0AV0Cvnsq;         //! correlation V0A-V0C for resolution
        TProfile*     fV0ATPCvnsq;         //! correlation V0A-TPC for resolution
        TProfile*     fV0CTPCvnsq;         //! correlation V0C-TPC for resolution

        TProfile*     fV0AV0CvnB1Csq;         //! correlation V0A-V0C for resolution
        TProfile*     fV0ATPCvnB1Csq;         //! correlation V0A-TPC for resolution
        TProfile*     fV0CTPCvnB1Csq;         //! correlation V0C-TPC for resolution


        TProfile*     fVnAllA[10];             //! vn V0A all
        TProfile*     fVnPiA[10];              //! vn V0A pi+ + pi-
        TProfile*     fVnKA[10];               //! vn V0A K+ + K-
        TProfile*     fVnAntiPA[10];           //! vn V0A anti-p
        TProfile*     fVnPihighPtA[10];        //! vn V0A pi+ + pi-
        TProfile*     fVnPhighPtA[10];         //! vn V0A anti-p + p

        TProfile*     fVnAllC[10];             //! vn V0C all
        TProfile*     fVnPiC[10];              //! vn V0C pi+ + pi-
        TProfile*     fVnKC[10];               //! vn V0C K+ + K-
        TProfile*     fVnAntiPC[10];           //! vn V0C anti-p
        TProfile*     fVnPihighPtC[10];        //! vn V0C pi+ + pi-
        TProfile*     fVnPhighPtC[10];         //! vn V0C anti-p + p


        TH1D*         fPsiA[10];               //! Psi V0A
        TH1D*         fPsiC[10];               //! Psi V0C
        TH2D*         fPsiAvsPsiC[10];         //! correlation PsiA-PsiC

        TProfile*     fSinTrkCosV0A[10];       //! sin(Trk)cos(V0A) for debug
        TProfile*     fCosTrkSinV0A[10];       //! cos(Trk)sin(V0A) for debug

        TProfile*     fSinTrkCosV0C[10];       //! sin(Trk)cos(V0C) for debug
        TProfile*     fCosTrkSinV0C[10];       //! cos(Trk)sin(V0C) for debug

        TProfile*     fSinTrkSinV0A[10];       //! sin(Trk)sin(V0A) for debug
        TProfile*     fCosTrkCosV0A[10];       //! cos(Trk)cos(V0A) for debug

        TProfile*     fSinTrkSinV0C[10];       //! sin(Trk)sin(V0C) for debug
        TProfile*     fCosTrkCosV0C[10];       //! cos(Trk)cos(V0C) for debug



        TH1D*         fInvMassK0[17][10];     //! inv mass K0
        TProfile*     fVnK0A[17][10];          //! vn V0A K0
        TProfile*     fVnK0C[17][10];          //! vn V0C K0

        TH1D*         fInvMassL[17][10];     //! inv mass Lambda
        TProfile*     fVnLA[17][10];          //! vn V0A Lambda
        TProfile*     fVnLC[17][10];          //! vn V0A Lambda

        TH1D*         fInvMassK0B1C[17][1];     //! inv mass K0
        TProfile*     fVnK0AB1C[17][1];          //! vn V0A K0
        TProfile*     fVnK0CB1C[17][1];          //! vn V0C K0

        TH1D*         fInvMassLB1C[17][1];     //! inv mass Lambda
        TProfile*     fVnLAB1C[17][1];          //! vn V0A Lambda
        TProfile*     fVnLCB1C[17][1];          //! vn V0A Lambda



        TProfile*     fVnAllAB1C[90];             //! vn V0A all
        TProfile*     fVnPiAB1C[90];              //! vn V0A pi+ + pi-
        TProfile*     fVnKAB1C[90];               //! vn V0A K+ + K-
        TProfile*     fVnAntiPAB1C[90];           //! vn V0A anti-p
        TProfile*     fVnPihighPtAB1C[90];        //! vn V0A pi+ + pi-
        TProfile*     fVnPhighPtAB1C[90];         //! vn V0A anti-p + p

        TProfile*     fVnAllCB1C[90];             //! vn V0C all
        TProfile*     fVnPiCB1C[90];              //! vn V0C pi+ + pi-
        TProfile*     fVnKCB1C[90];               //! vn V0C K+ + K-
        TProfile*     fVnAntiPCB1C[90];           //! vn V0C anti-p
        TProfile*     fVnPihighPtCB1C[90];        //! vn V0C pi+ + pi-
        TProfile*     fVnPhighPtCB1C[90];         //! vn V0C anti-p + p


        TH1D*         fPsiAB1C[90];               //! Psi V0A
        TH1D*         fPsiCB1C[90];               //! Psi V0C
        TH2D*         fPsiAvsPsiCB1C[90];         //! correlation PsiA-PsiC

        TProfile*     fSinTrkCosV0AB1C[90];       //! sin(Trk)cos(V0A) for debug
        TProfile*     fCosTrkSinV0AB1C[90];       //! cos(Trk)sin(V0A) for debug

        TProfile*     fSinTrkCosV0CB1C[90];       //! sin(Trk)cos(V0C) for debug
        TProfile*     fCosTrkSinV0CB1C[90];       //! cos(Trk)sin(V0C) for debug

        TProfile*     fSinTrkSinV0AB1C[90];       //! sin(Trk)sin(V0A) for debug
        TProfile*     fCosTrkCosV0AB1C[90];       //! cos(Trk)cos(V0A) for debug

        TProfile*     fSinTrkSinV0CB1C[90];       //! sin(Trk)sin(V0C) for debug
        TProfile*     fCosTrkCosV0CB1C[90];       //! cos(Trk)cos(V0C) for debug


        AliAnalysisTaskPiKpK0Lamba(const AliAnalysisTaskPiKpK0Lamba&);      // not implemented
        AliAnalysisTaskPiKpK0Lamba& operator=(const AliAnalysisTaskPiKpK0Lamba&);   // not implemented

        ClassDef(AliAnalysisTaskPiKpK0Lamba, 1);    //Analysis task for high pt analysis 
};

#endif
