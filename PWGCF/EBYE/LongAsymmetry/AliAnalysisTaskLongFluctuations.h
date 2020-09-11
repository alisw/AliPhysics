////////////////////////////////////////////////////////
// AliAnalysisTaskLongFluctuations:
// Description: Analysis task for 2PC for eta1, eta2
// Author: Raquel Quishpe (raquel.quishpe@cern.ch)
////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskLongFluctuations_H
#define AliAnalysisTaskLongFluctuations_H

#include "AliAnalysisTaskSE.h"

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


    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TNtuple*                fNt;
        
        TH3D*                   fEtaGlobCentPVz;
        TH3D*                   fEtaTPCCentPVz;
        TH3D*                   fEtaHybCentPVz;
        TH2D*                   fEtaMCCent;
    
        Bool_t fIsMC;
        Double_t fChi2DoF;
        Int_t fTPCNcls;
        Double_t fPtmin;
        Double_t fPtmax;
        Double_t fEta;

        AliAnalysisTaskLongFluctuations(const AliAnalysisTaskLongFluctuations&); // not implemented
        AliAnalysisTaskLongFluctuations& operator=(const AliAnalysisTaskLongFluctuations&); // not implemented

        Float_t GetAnCoeff(Int_t order, TH1D *hist);
        Float_t LegPol(Int_t order, Double_t x);
    
        ClassDef(AliAnalysisTaskLongFluctuations, 1);
};

#endif
