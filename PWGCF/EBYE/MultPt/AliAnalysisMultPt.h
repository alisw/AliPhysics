#ifndef AliAnalysisMultPt_H
#define AliAnalysisMultPt_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#define PI 3.1415927


class AliAnalysisMultPt : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisMultPt();
                                AliAnalysisMultPt(const char *name);
        virtual                 ~AliAnalysisMultPt();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        
        void            SetMCRead(Bool_t flag);
        void            SetPtLimits(Double_t ptmin, Double_t ptmax) { fPtmin = ptmin; fPtmax=ptmax; }
        void            SetEtaLimit(Double_t etalimit) { fEta = etalimit; }
        void            SetFilterBit(Int_t filterbit) { fBit = filterbit; }
        
 
    private:
        AliAODEvent             *fAOD;             //! input event
        TList                   *fOutputList;      //! output list
        TH2D                    *fHistMultPt;      //! pT-Mult histogram
        TH1D                    *fHistMult;        //! Mult histogram
        TH2D                    *fHistMultPtMC;    //! pT-Mult MC histogram
        TH2D                    *fHistMultRatio;   //! Ratio of Gen mult vs rec mult histogram
        Bool_t   fIsMC;                            //! MC flag
        Double_t fPtmin;                           //! min PT
        Double_t fPtmax;                           //! max PT
        Double_t fEta;                             //! max eta
        Int_t    fBit;                             //! filter bit 

 

    AliAnalysisMultPt(const AliAnalysisMultPt&); // not implemented
    AliAnalysisMultPt& operator=(const AliAnalysisMultPt&); // not implemented

    ClassDef(AliAnalysisMultPt, 2);
};

#endif
