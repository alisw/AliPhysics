#ifndef ALIANALYSISTASKCMEESE_H
#define ALIANALYSISTASKCMEESE_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TSpline.h>
#include <THnSparse.h>

// AliRoot includes
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"


class AliAnalysisTaskCmeEse : public AliAnalysisTaskSE {
    public:
        AliAnalysisTaskCmeEse();
        AliAnalysisTaskCmeEse(const char *name);

        virtual ~AliAnalysisTaskCmeEse();

        virtual void   UserCreateOutputObjects();
        virtual void   UserExec(Option_t *option);
        virtual void   Terminate(Option_t *); 

        virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
        virtual void  SetFilterbit(UInt_t filterbit){fFilterbit = filterbit;}	
        virtual void  SetNoClus(Int_t noclus){fNoClus = noclus;}	
        virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
        virtual void  SetMinPt(Double_t minPt){fMinPt = minPt;}
        virtual void  SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
        virtual void  SetFlagLHC10h(Bool_t IsLHC10h){fLHC10h = IsLHC10h;}
        virtual void  SetFlagPileUp(Bool_t IsPileUP){fPileUp = IsPileUP;}
        virtual void  SetCentFlag(Short_t nCent){fCent = nCent;}
        virtual void  SetFlagRecEff(Bool_t IsRec){fRecEff = IsRec;}
        virtual void  SetFlagQAV0(Bool_t isQAV0){fQAV0 = isQAV0;}
        virtual void  SetFlagQATrk(Bool_t isTrkQA){fTrkQA = isTrkQA;}


    private:
        virtual Float_t GetVertex(AliAODEvent* aod) const;
        virtual void Analyze(AliAODEvent* aod, Float_t vtxZ);
        void OpenInfoCalbration(Int_t run);
        Short_t GetPercCode(Double_t perc) const;
        virtual Int_t GetTPCMult(AliVEvent* ev) const;
        virtual Int_t GetGlobalMult(AliVEvent* ev) const;
        Double_t GetWDist(const AliVVertex* v0, const AliVVertex* v1);
        Bool_t plpMV(const AliVEvent *event);


        AliAODEvent* fAOD;                //! AOD object

        Int_t        fRun;                // run number - for calibration
        TH1D*        fMultV0;             // profile from V0 multiplicity

        TH1D*        fQx2mV0A;            // <Qx2> V0A
        TH1D*        fQy2mV0A;            // <Qy2> V0A
        TH1D*        fQx2sV0A;            // sigma Qx2 V0A
        TH1D*        fQy2sV0A;            // sigma Qy2 V0A
        TH1D*        fQx2mV0C;            // <Qx2> V0C
        TH1D*        fQy2mV0C;            // <Qy2> V0C
        TH1D*        fQx2sV0C;            // sigma Qx2 V0C
        TH1D*        fQy2sV0C;            // sigma Qy2 V0C

        // Cuts and options
        TF1*         fFitRecLow;          // rec eff fit 1
        TF1*         fFitRecHigh;         // rec eff fit 2
        Bool_t       fRecEff;             // rec eff flag
        Double_t     fVtxCut;             // Vtx cut on z position in cm
        UInt_t       fFilterbit;          // filter bit
        Double_t     fEtaCut;             // Eta cut used to select particles
        Int_t        fNoClus;	          // No of TPC clusters
        Double_t     fMinPt;              // Min pt - for histogram limits
        Double_t     fMaxPt;              // Max pt - for histogram limits
        Bool_t       fLHC10h;             // flag to LHC10h data
        Bool_t       fPileUp;             // flag for pileup
        Short_t      fCent;               // centrality flag
        Bool_t       fQAV0;               // qa V0 flag
        Bool_t       fTrkQA;              // qa for tracks
        TSpline3*    fSplQ2c[90];         // splines for q2 V0C

        // Output objects
        TList*        fListOfObjects;     //! Output list of objects
        TH1I*         fVtx;               //! Event vertex info
        TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
        TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts  
        TH2F*         fMultCorBeforeCuts; //! correlation between TPC and global track multiplicity before cuts
        TH2F*         fMultCorAfterCuts;  //! correlation between TPC and global track multiplicity after cuts

        TH2D*         fPercqc2;          //! percentile distribution

        THnSparseF*    fAllQA;               //! qa histo for tracks

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


        TProfile*     fV2A;              //! v2 V0A  
        TProfile*     fV0AV0Cv2;         //! correlation V0A-V0C for resolution 
        TProfile*     fV0ATPCv2;         //! correlation V0A-TPC for resolution 
        TProfile*     fV0CTPCv2;         //! correlation V0C-TPC for resolution 

        TProfile*     fCMESameQPos;  //! cme same charge pos
        TProfile*     fCMESameQNeg;  //! cme same charge neg
        TProfile*     fCMEOppQ;      //! cme opp charge

        TProfile*     fCMESameQPosCos;  //! cme same charge pos <cos*cos>
        TProfile*     fCMESameQPosSin;  //! cme same charge pos <sin*sin>
        TProfile*     fCMESameQPosv1;   //! cme same charge pos v1

        TProfile*     fCMESameQNegCos;  //! cme same charge neg <cos*cos>
        TProfile*     fCMESameQNegSin;  //! cme same charge neg <sin*sin>
        TProfile*     fCMESameQNegv1;   //! cme same charge neg v1

        TProfile*     fCMEOppQCos;  //! cme opp charge <cos*cos>
        TProfile*     fCMEOppQSin;  //! cme opp charge <sin*sin>
        TProfile*     fCMEOppQv1;   //! cme opp charge v1



        TProfile*     fV2Aqc2[10];        //! v2 V0A ESE
        TProfile*     fV0AV0Cv2qc2[10];   //! correlation V0A-V0C for resolution ESE
        TProfile*     fV0ATPCv2qc2[10];   //! correlation V0A-TPC for resolution ESE
        TProfile*     fV0CTPCv2qc2[10];   //! correlation V0C-TPC for resolution ESE

        TProfile*     fCMESameQPosqc2[10];  //! cme same charge pos ESE
        TProfile*     fCMESameQNegqc2[10];  //! cme same charge neg ESE
        TProfile*     fCMEOppQqc2[10];      //! cme opp charge ESE

        TProfile*     fCMESameQPosCosqc2[10];  //! cme same charge pos <cos*cos>
        TProfile*     fCMESameQPosSinqc2[10];  //! cme same charge pos <sin*sin>
        TProfile*     fCMESameQPosv1qc2[10];   //! cme same charge pos v1

        TProfile*     fCMESameQNegCosqc2[10];  //! cme same charge neg <cos*cos>
        TProfile*     fCMESameQNegSinqc2[10];  //! cme same charge neg <sin*sin>
        TProfile*     fCMESameQNegv1qc2[10];   //! cme same charge neg v1

        TProfile*     fCMEOppQCosqc2[10];  //! cme opp charge <cos*cos>
        TProfile*     fCMEOppQSinqc2[10];  //! cme opp charge <sin*sin>
        TProfile*     fCMEOppQv1qc2[10];   //! cme opp charge v1


        TProfile*     fV2Pt[9];              //! v2 pt unbiased

        TProfile*     fCMESameQPosEta[9];    //! cme same charge pos eta
        TProfile*     fCMESameQPosPtDif[9];  //! cme same charge pos pt dif
        TProfile*     fCMESameQPosPtSum[9];  //! cme same charge pos pt sum

        TProfile*     fCMESameQNegEta[9];    //! cme same charge neg eta
        TProfile*     fCMESameQNegPtDif[9];  //! cme same charge neg pt dif
        TProfile*     fCMESameQNegPtSum[9];  //! cme same charge neg pt sum

        TProfile*     fCMEOppQEta[9];        //! cme opp charge eta
        TProfile*     fCMEOppQPtDif[9];      //! cme opp charge pt dif
        TProfile*     fCMEOppQPtSum[9];      //! cme opp charge pt sum


        TH1D*         fPsiA[10];               //! Psi V0A
        TH1D*         fPsiC[10];               //! Psi V0C
        TH2D*         fPsiAvsPsiC[10];         //! correlation PsiA-PsiC


        TProfile*     fV2Ptqc2[9][10];            //! v2 pt ESE

        TProfile*     fCMESameQPosEtaqc2[9][10];  //! cme same charge pos ESE
        TProfile*     fCMESameQPosPtDifqc2[9][10];  //! cme same charge pos ESE
        TProfile*     fCMESameQPosPtSumqc2[9][10];  //! cme same charge pos ESE

        TProfile*     fCMESameQNegEtaqc2[9][10];  //! cme same charge neg ESE
        TProfile*     fCMESameQNegPtDifqc2[9][10];  //! cme same charge neg ESE
        TProfile*     fCMESameQNegPtSumqc2[9][10];  //! cme same charge neg ESE

        TProfile*     fCMEOppQEtaqc2[9][10];      //! cme opp charge ESE
        TProfile*     fCMEOppQPtDifqc2[9][10];      //! cme opp charge ESE
        TProfile*     fCMEOppQPtSumqc2[9][10];      //! cme opp charge ESE


        AliAnalysisTaskCmeEse(const AliAnalysisTaskCmeEse&);
        AliAnalysisTaskCmeEse& operator=(const AliAnalysisTaskCmeEse&);

        ClassDef(AliAnalysisTaskCmeEse, 1);
};

#endif
