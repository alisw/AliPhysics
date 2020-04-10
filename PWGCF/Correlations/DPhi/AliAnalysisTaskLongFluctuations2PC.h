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
        void            SetFilterBit(Int_t fb);
        void            SetCentrality(TString Cent);
        void            SetChi2DoF(Double_t Chi2DoF);
        void            SetNclTPC(Int_t ncl);
        void            SetPtLimits(Double_t ptmin, Double_t ptmax);
        void            SetEtaLimit(Double_t etalimit);

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TH1D*                   fHistBg[9][17];        //! 2PC background histogram (eta)
        TH2D*                   fHistSig[9][17];        //! 2PC signal histogram (eta1, eta2)
        TH2D*                   fHistBgGap[9][17];        //! 2PC background histogram with gap (eta,phi)
        TH2D*                   fHistSigGap[9][17];       //! 2PC signal histogram with gap (eta1, eta2)

        TH1D*                   fHistBgMC[9];        //! 2PC MC background histogram (eta)
        TH2D*                   fHistSigMC[9];        //! 2PC MC signal histogram (eta1, eta2)
        TH2D*                   fHistBgGapMC[9];        //! 2PC MC background histogram with gap (eta,phi)
        TH2D*                   fHistSigGapMC[9];       //! 2PC MC signal histogram with gap (eta1, eta2)
        TH1D*                   fMultD;
        TProfile*               fMultDpT;
        TH1D*                   fMultDMC;
        TProfile*               fMultDpTMC;
        TH2D*                   fMultDReMC;
        TH2D*                   fPVzCentNevents;
        TH2D*                   fCentMultQA;
        TH3D*                   fAllCentQA;
    
        Bool_t fIsMC;
        Int_t fFB;
        TString fCentrality;
        Double_t fChi2DoF;
        Int_t fTPCNcls;
        Double_t fPtmin;
        Double_t fPtmax;
        Double_t fEta;

        AliAnalysisTaskLongFluctuations2PC(const AliAnalysisTaskLongFluctuations2PC&); // not implemented
        AliAnalysisTaskLongFluctuations2PC& operator=(const AliAnalysisTaskLongFluctuations2PC&); // not implemented

        ClassDef(AliAnalysisTaskLongFluctuations2PC, 2);
};

#endif
