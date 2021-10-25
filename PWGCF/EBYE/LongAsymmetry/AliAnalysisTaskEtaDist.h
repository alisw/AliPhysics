////////////////////////////////////////////////////////
// AliAnalysisTaskEtaDist:
// Author: Raquel Quishpe (raquel.quishpe@cern.ch)
////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskEtaDist_H
#define AliAnalysisTaskEtaDist_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

#define PI 3.1415927

class AliAnalysisTaskEtaDist : public AliAnalysisTaskSE
{
    public:
                        AliAnalysisTaskEtaDist();
                        AliAnalysisTaskEtaDist(const char *name);
        virtual         ~AliAnalysisTaskEtaDist();

        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t* option);
        virtual void    Terminate(Option_t* option);

        void            SetMCRead(Bool_t flag) { fIsMC = flag; }
        void            SetChi2DoF(Double_t Chi2DoF) { fChi2DoF = Chi2DoF; }
        void            SetNclTPC(Int_t ncl) { fTPCNcls = ncl; }
        void            SetPtLimits(Double_t ptmin, Double_t ptmax) { fPtmin = ptmin; fPtmax=ptmax; }
        void            SetEtaLimit(Double_t etalimit) { fEta = etalimit; }
        void            SetFilterBit(Int_t filterbit) { fBit = filterbit; }
        void            SetPileUpRead(Bool_t flag) {fIsPileUpCuts = flag;}

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TTree*                fNt;
        
        Bool_t fIsMC;
        Double_t fChi2DoF;
        Int_t fTPCNcls;
        Double_t fPtmin;
        Double_t fPtmax;
        Double_t fEta;
        Int_t fBit; //filter bit - 96 default
        Bool_t fIsPileUpCuts;
        AliEventCuts fEventCuts;

        Float_t mCentV0M;
        Float_t mCentCL0;
        Float_t mCentCL1;
        Float_t mPVz;
        Float_t mZDCN1;
        Float_t mZDCN2;
        Float_t mNTLs;
        Float_t mNGlob;
        Float_t mNHyb;
        Float_t mSpTGlob;
        Float_t mSpTHyb;
        Float_t mNMC;
        Float_t mSpTMC;
    
        Float_t aEtaPos[16];
        Float_t aEtaNeg[16];
        Float_t aEtaPosMC[16];
        Float_t aEtaNegMC[16];

    
        AliAnalysisTaskEtaDist(const AliAnalysisTaskEtaDist&); // not implemented
        AliAnalysisTaskEtaDist& operator=(const AliAnalysisTaskEtaDist&); // not implemented
    
        ClassDef(AliAnalysisTaskEtaDist, 2);
};

#endif
