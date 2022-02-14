////////////////////////////////////////////////////////
// AliAnalysisTaskLongFluctuations:
// Description: Analysis task for 2PC for eta1, eta2
// Author: Raquel Quishpe (raquel.quishpe@cern.ch)
////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskLongFluctuations_H
#define AliAnalysisTaskLongFluctuations_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

#define PI 3.1415927

class AliAnalysisTaskLongFluctuations : public AliAnalysisTaskSE
{
    public:
                        AliAnalysisTaskLongFluctuations();
                        AliAnalysisTaskLongFluctuations(const char *name);
        virtual         ~AliAnalysisTaskLongFluctuations();

        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t* option);
        virtual void    Terminate(Option_t* option);
        void            SetMCRead(Bool_t bs);
        void            SetChi2DoF(Double_t Chi2DoF);
        void            SetNclTPC(Int_t ncl);
        void            SetPtLimits(Double_t ptmin, Double_t ptmax);
        void            SetEtaLimit(Double_t etalimit);
        void            SetPileUpRead(Bool_t ps);

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

    
        AliAnalysisTaskLongFluctuations(const AliAnalysisTaskLongFluctuations&); // not implemented
        AliAnalysisTaskLongFluctuations& operator=(const AliAnalysisTaskLongFluctuations&); // not implemented
    
        ClassDef(AliAnalysisTaskLongFluctuations, 2);
};

#endif
