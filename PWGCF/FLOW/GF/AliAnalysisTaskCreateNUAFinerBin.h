/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCreateNUAFinerBin_H
#define AliAnalysisTaskCreateNUAFinerBin_H

#include "AliAnalysisTaskSE.h"
#include "AliGFWWeightsFinerBin.h"
#include "AliGFWXXCuts.h"
#include "AliGFWPbpass23Cuts.h"
#include "AliEventCuts.h"
#include <vector>

class AliAnalysisTaskCreateNUAFinerBin : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskCreateNUAFinerBin();
                                AliAnalysisTaskCreateNUAFinerBin(const char *name);
        virtual                 ~AliAnalysisTaskCreateNUAFinerBin();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void   SetPeriod(TString period) { fPeriod = period; }
        virtual void   SetSystFlag(int flag) { fCurrSystFlag = flag; }
        virtual int    GetSystFlag() { return fCurrSystFlag; }
        virtual void   SetMinPt(Double_t minPt){fMinPt = minPt;}
		virtual void   SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
		virtual void   SetUseHM(Bool_t useHM){fUseHM = useHM;}
		virtual void   SetUseCuts(Bool_t useCuts){fUseCuts = useCuts;}

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        AliGFWWeightsFinerBin*          weights;        //!
        std::vector<int>        RunNumber_LHC16;//!
        std::vector<int>        RunNumber_LHC17;//!
        std::vector<int>        RunNumber_LHC18;//!
        std::vector<TString>    Period_LHC16;//!
        std::vector<TString>    Period_LHC17;//!
        std::vector<TString>    Period_LHC18;//!
        std::vector<int>        RunNumber_LHC15o;//!
        std::vector<int>        RunNumber_LHC15opass2;//!
        std::vector<int>        RunNumber_LHC18qrpass3;//!
        std::vector<int>        RunNumber_LHC16qt;//!
        std::vector<int>        RunNumber_LHC17n;//!
        Int_t                   Last_RunNumer;  //!
        Int_t                   Last_Position;  //!

        TString			        fPeriod;		// period
        Int_t                   fCurrSystFlag;  // Systematics flag
        Bool_t                  fIsMC;          //
        Bool_t                  fUseHM;         // if use High Multiplicity Trigger
        Bool_t                  fUseCuts;       // if use Cuts in code
        Bool_t                  fPbPb;          // if use Pb-Pb dataset

        //Standard Cuts and Systematics
        AliEventCuts	fEventCuts;					// Event cuts
		AliGFWXXCuts*     fGFWSelection;                                  //!
    AliGFWPbpass23Cuts*     fGFWSelection15o;                                  //!


        //cuts
		Double_t		fMinPt;					// Min pt - for histogram limits
		Double_t		fMaxPt;					// Max pt - for histogram limits

        //Events histogram
        TH1D*			hEventCount;			//! counting events passing given event cuts

        AliAnalysisTaskCreateNUAFinerBin(const AliAnalysisTaskCreateNUAFinerBin&); // not implemented
        AliAnalysisTaskCreateNUAFinerBin& operator=(const AliAnalysisTaskCreateNUAFinerBin&); // not implemented
        const char* ReturnPPperiod(const Int_t runNumber) const; //!
        const Int_t ReturnPosi_WeightList(const Int_t RunNumber, const Int_t sysflag); //!
        Bool_t AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp);//!

        ClassDef(AliAnalysisTaskCreateNUAFinerBin, 1);
};

#endif