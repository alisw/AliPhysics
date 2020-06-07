////////////////////////////////////////////////////////
// AliAnalysisTaskLongFluctuations2PC:
// Description: Analysis task for 2PC for eta1, eta2
// Author: Raquel Quishpe (raquel.quishpe@cern.ch)
////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskLongFluctuations2PC_H
#define AliAnalysisTaskLongFluctuations2PC_H

#include "AliAnalysisTaskSE.h"

#define PI 3.1415927

class AliAnalysisTaskLongFluctuations2PC : public AliAnalysisTaskSE
{
    public:
                        AliAnalysisTaskLongFluctuations2PC();
                        AliAnalysisTaskLongFluctuations2PC(const char *name);
        virtual         ~AliAnalysisTaskLongFluctuations2PC();

        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t* option);
        virtual void    Terminate(Option_t* option);
        void            SetMCRead(Bool_t bs);
        void            SetCentrality(TString Cent);
        void            SetChi2DoF(Double_t Chi2DoF);
        void            SetNclTPC(Int_t ncl);
        void            SetPtLimits(Double_t ptmin, Double_t ptmax);
        void            SetEtaLimit(Double_t etalimit);

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TH2D*                   fHistBgTPC[9][17];        //! 2PC background histogram (eta,phi) TPC
        TH3D*                   fHistSigTPC[9][17];        //! 2PC signal histogram (eta1, eta2,d,phi) TPC
        TH2D*                   fHistBgGlob[9][17];        //! 2PC background histogram with gap (eta,phi) Globals
        TH3D*                   fHistSigGlob[9][17];       //! 2PC signal histogram with gap (eta1, eta2,dphi) Globals

        TH2D*                   fHistBgMC[9];        //! 2PC MC background histogram (eta)
        TH3D*                   fHistSigMC[9];        //! 2PC MC signal histogram (eta1, eta2)
        TH1D*                   fMultD;
        TProfile*               fMultDpT;
        TH1D*                   fMultDMC;
        TProfile*               fMultDpTMC;
        TH2D*                   fMultDReMC;
        TH2D*                   fPVzCentNevents;
        TH2D*                   fCentMultQA;
        TH3D*                   fAllCentQA;
    
        Bool_t fIsMC;
        TString fCentrality;
        Double_t fChi2DoF;
        Int_t fTPCNcls;
        Double_t fPtmin;
        Double_t fPtmax;
        Double_t fEta;

        AliAnalysisTaskLongFluctuations2PC(const AliAnalysisTaskLongFluctuations2PC&); // not implemented
        AliAnalysisTaskLongFluctuations2PC& operator=(const AliAnalysisTaskLongFluctuations2PC&); // not implemented

        ClassDef(AliAnalysisTaskLongFluctuations2PC, 3);
};

#endif
