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

        virtual void   SetMinPt(Double_t minPt){ fMinPt = minPt; }
		virtual void   SetMaxPt(Double_t maxPt){ fMaxPt = maxPt; }
        virtual void   SetVz(Double_t Vz){ fVz = Vz; }
        virtual void   SetDCAz(Double_t dcaz){ fDCAz = dcaz; }
        virtual void   SetDCAxy(Double_t dcaxy){ fDCAxy = dcaxy; }
        virtual void   SetTPCcls(Int_t tpc){ fTPCcls = tpc;}
        virtual void   SetFilterBit(Int_t filterbit){ filterBit = filterbit; }
        virtual void   SetTrigger(Int_t trig){ fTrigger = trig; }
        virtual void   SetCtrType(TString ctrtype){ ctrType = ctrtype; }
        virtual void   SetPeriod(TString period) { fPeriod = period; }
        virtual void   SetNUE(TString NUE) { fNUE = NUE; }

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        Double_t GetWeightNUE(double pt);

       

        AliAODEvent*            fAOD;           //! input event

        TList*                  fOutputList;    //! output list
        TList*                  fBstList;
        TList*                  fWeightsListNUE;
        TH1D*                   fWeightNUE;
                
        Double_t		fMinPt;					// Min pt - for histogram limits
		Double_t		fMaxPt;					// Max pt - for histogram limits
        Double_t        fVz;
        Double_t        fDCAz;
        Double_t        fDCAxy;
        Int_t           fTPCcls;
        Int_t           filterBit;

		Int_t			fTrigger;				// flag for trigger
        TString         ctrType;
        TString			fPeriod;				// period
        TString         fNUE;
        
        TH1F*                   fTestNonWeight;
        TH1F*                   fNchDistri;
        TH2F*                   fPtNch;
        TProfile*               dPtNch;
        TProfile*               dPt2Nch;
        TProfile*               dPt3Nch;
        TProfile*               dPt4Nch;
        TProfile*               dPtNchNB;

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

        TProfile*               dPt4Nch_0;
        TProfile*               dPt4Nch_1;
        TProfile*               dPt4Nch_2;
        TProfile*               dPt4Nch_3;
        TProfile*               dPt4Nch_4;
        TProfile*               dPt4Nch_5;
        TProfile*               dPt4Nch_6;
        TProfile*               dPt4Nch_7;
        TProfile*               dPt4Nch_8;
        TProfile*               dPt4Nch_9;
        
        TH2F*                   fPtNchUCC;
        TProfile*               TestPtCtr;
        TProfile*               TestPt2Ctr;
        TProfile*               TestPt3Ctr;
        TProfile*               TestPt4Ctr;
        TProfile*               TestNchCtr;
        TProfile*               TestNchSelectedCtr;

        TProfile*               TestPtCtr_0;
        TProfile*               TestPtCtr_1;
        TProfile*               TestPtCtr_2;
        TProfile*               TestPtCtr_3;
        TProfile*               TestPtCtr_4;
        TProfile*               TestPtCtr_5;
        TProfile*               TestPtCtr_6;
        TProfile*               TestPtCtr_7;
        TProfile*               TestPtCtr_8;
        TProfile*               TestPtCtr_9;

        TProfile*               TestPt2Ctr_0;
        TProfile*               TestPt2Ctr_1;
        TProfile*               TestPt2Ctr_2;
        TProfile*               TestPt2Ctr_3;
        TProfile*               TestPt2Ctr_4;
        TProfile*               TestPt2Ctr_5;
        TProfile*               TestPt2Ctr_6;
        TProfile*               TestPt2Ctr_7;
        TProfile*               TestPt2Ctr_8;
        TProfile*               TestPt2Ctr_9;

        TProfile*               TestPt3Ctr_0;
        TProfile*               TestPt3Ctr_1;
        TProfile*               TestPt3Ctr_2;
        TProfile*               TestPt3Ctr_3;
        TProfile*               TestPt3Ctr_4;
        TProfile*               TestPt3Ctr_5;
        TProfile*               TestPt3Ctr_6;
        TProfile*               TestPt3Ctr_7;
        TProfile*               TestPt3Ctr_8;
        TProfile*               TestPt3Ctr_9;

        TProfile*               TestPt4Ctr_0;
        TProfile*               TestPt4Ctr_1;
        TProfile*               TestPt4Ctr_2;
        TProfile*               TestPt4Ctr_3;
        TProfile*               TestPt4Ctr_4;
        TProfile*               TestPt4Ctr_5;
        TProfile*               TestPt4Ctr_6;
        TProfile*               TestPt4Ctr_7;
        TProfile*               TestPt4Ctr_8;
        TProfile*               TestPt4Ctr_9;
        
        TProfile*               TestNchCtr_0;
        TProfile*               TestNchCtr_1;
        TProfile*               TestNchCtr_2;
        TProfile*               TestNchCtr_3;
        TProfile*               TestNchCtr_4;
        TProfile*               TestNchCtr_5;
        TProfile*               TestNchCtr_6;
        TProfile*               TestNchCtr_7;
        TProfile*               TestNchCtr_8;
        TProfile*               TestNchCtr_9;

        TProfile*               TestNchSelectedCtr_0;
        TProfile*               TestNchSelectedCtr_1;
        TProfile*               TestNchSelectedCtr_2;
        TProfile*               TestNchSelectedCtr_3;
        TProfile*               TestNchSelectedCtr_4;
        TProfile*               TestNchSelectedCtr_5;
        TProfile*               TestNchSelectedCtr_6;
        TProfile*               TestNchSelectedCtr_7;
        TProfile*               TestNchSelectedCtr_8;
        TProfile*               TestNchSelectedCtr_9;

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
