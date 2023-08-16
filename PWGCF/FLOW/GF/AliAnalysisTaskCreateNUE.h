/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCreateNUE_H
#define AliAnalysisTaskCreateNUE_H

#include "AliAnalysisTaskSE.h"
#include "AliGFWWeights.h"
#include "AliGFWMCuts.h"
#include "AliGFWNFCuts.h"
#include "AliEventCuts.h"
#include "AliAODMCParticle.h"
#include <vector>

class AliAnalysisTaskCreateNUE : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskCreateNUE();
                                AliAnalysisTaskCreateNUE(const char *name);
        virtual                 ~AliAnalysisTaskCreateNUE();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void            SetPeriod(TString period) { fPeriod = period; }
        virtual void            SetSystFlag(int flag) { fCurrSystFlag = flag; }
        virtual int             GetSystFlag() { return fCurrSystFlag; }
        virtual void            SetMinPt(Double_t minPt){fMinPt = minPt;}
		    virtual void            SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
		    virtual void            SetUseHM(Bool_t useHM){fUseHM = useHM;}
		    virtual void            SetUseCuts(Bool_t useCuts){fUseCuts = useCuts;}

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list

        TString			            fPeriod;		    // period
        Int_t                   fCurrSystFlag;  // Systematics flag
        Bool_t                  fIsMC;          //
        Bool_t                  fUseHM;         // if use High Multiplicity Trigger
        Bool_t                  fUseCuts;       // if use Cuts in code

        //Standard Cuts and Systematics
        AliEventCuts	    fEventCuts;					// Event cuts
		    AliGFWMCuts*      fGFWSelection;                                  //!
        AliGFWNFCuts*     fGFWSelection15o;                                  //!


        //cuts
		    Double_t		fMinPt;					// Min pt - for histogram limits
		    Double_t		fMaxPt;					// Max pt - for histogram limits
        Double_t    fEtaCut;        // Eta Range (0.8)

        //Events histogram
        TH1D*			hPtMCGen;			    //! MC Gen Pt distribution
        TH1D*			hPtMCRec;			    //! MC Rec Pt distribution
        TH1D*			hEventCount;			//! counting events passing given event cuts

        AliAnalysisTaskCreateNUE(const AliAnalysisTaskCreateNUE&); // not implemented
        AliAnalysisTaskCreateNUE& operator=(const AliAnalysisTaskCreateNUE&); // not implemented
        Bool_t AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp);//!
  Bool_t AcceptMCTruthTrack(AliAODMCParticle *mtr); //!

        ClassDef(AliAnalysisTaskCreateNUE, 1);
};

#endif
