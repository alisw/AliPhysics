/* AliAnaysisTaskPtN
 * Maintainer: Yifan Zhang
 * calculating the mean Pt and fluctuations with respect to Nch
 */

#ifndef AliAnalysisPtN_cxx
#define AliAnalysisPtN_cxx

#include "AliAnalysisTaskSE.h"
#include <TList.h>
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "AliMultSelectionTask.h"
// #include "CorrelationCalculator.h"
#include "AliAODEvent.h"
#include "AliEventCuts.h"
#include <TComplex.h>
#include <TObject.h>
# include <TRandom3.h>

// class TList;
// class TH1F;
// class TH3F;
// class TProfile;
// class TComplex;
// class AliAODEvent;

//class AliPIDResponse;

class AliAnalysisPtN : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisPtN();
                                AliAnalysisPtN(const char *name);
        virtual                 ~AliAnalysisPtN();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        Double_t GetWeightNUE(double pt);

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TList*                  fWeightsListNUE;
        TList*                  fBstList;
        TH1D*                   fWeightNUE;
        TH1F*                   fTestNonWeight;
        TH2F*                   fPtNch;
        TProfile*               dPtNch;
        TProfile*               dPt2Nch;
        TProfile*               dPt3Nch;

        TProfile*               dPtNch_0;
        TProfile*               dPtNch_1;
        TProfile*               dPtNch_2;
        TProfile*               dPtNch_3;
        TProfile*               dPtNch_4;
        TProfile*               dPtNch_5;
        TProfile*               dPtNch_6;
        TProfile*               dPtNch_7;
        TProfile*               dPtNch_8;
        TProfile*               dPtNch_9;
        TProfile*               dPt2Nch_0;
        TProfile*               dPt2Nch_1;
        TProfile*               dPt2Nch_2;
        TProfile*               dPt2Nch_3;
        TProfile*               dPt2Nch_4;
        TProfile*               dPt2Nch_5;
        TProfile*               dPt2Nch_6;
        TProfile*               dPt2Nch_7;
        TProfile*               dPt2Nch_8;
        TProfile*               dPt2Nch_9;
        TProfile*               dPt3Nch_0;
        TProfile*               dPt3Nch_1;
        TProfile*               dPt3Nch_2;
        TProfile*               dPt3Nch_3;
        TProfile*               dPt3Nch_4;
        TProfile*               dPt3Nch_5;
        TProfile*               dPt3Nch_6;
        TProfile*               dPt3Nch_7;
        TProfile*               dPt3Nch_8;
        TProfile*               dPt3Nch_9;
        
        TH2F*                   fPtNchUCC;
        TProfile*               TestPtCtr;
        TProfile*               TestPt2Ctr;
        TProfile*               TestPt3Ctr;
        TProfile*               TestNchCtr;
        TProfile*               TestNchSelectedCtr;
        
        AliMultSelection*       multSelection;  //centrality class

        AliEventCuts	  fEventCuts;					// Event cuts
        TRandom3          radm;
        //CorrelationCalculator  correlator;
       
   //     AliPIDResponse*         fPIDResponse; //! pid response objec

        AliAnalysisPtN(const AliAnalysisPtN&); // not implemented
        AliAnalysisPtN& operator=(const AliAnalysisPtN&); // not implemented
        ClassDef(AliAnalysisPtN, 1);
};

#endif
