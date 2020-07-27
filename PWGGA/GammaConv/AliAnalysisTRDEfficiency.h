/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTRDEfficiency_H
#define AliAnalysisTRDEfficiency_H

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "AliV0ReaderV1.h"
#include "AliTRDonlineTrackMatching.h"
//#include "AliAnalysisTaskESDfilter.h"
#include "AliConversionSelection.h"


class AliPIDResponse;

class AliAnalysisTRDEfficiency : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTRDEfficiency();
                                AliAnalysisTRDEfficiency(const char *name);
        virtual                 ~AliAnalysisTRDEfficiency();
        
        //virtual Bool_t           UserNotify();
        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual Bool_t          checkPi0(TClonesArray* lst, AliAODConversionMother* pi0, Double_t tmp[16], Int_t lbl);
        virtual Int_t           GetEventCuts(AliESDTrdTrack* trdtrack, TString clss);
        //virtual Bool_t          GetTrackCuts(AliESDtrack* track);
        virtual Double_t        GetSagitta(AliESDTrdTrack* trdtrack);
        virtual Double_t        GetRating(AliESDv0 *v0, AliESDtrack *track, AliESDTrdTrack *trdtrack);
        virtual void            Terminate(Option_t* option);

        void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
        
    private:
        AliESDEvent*            fESD;           //! input event
        TList*                  fOutputList;    //! output list
        AliV0ReaderV1*          fV0Reader;      //
        TString                 fV0ReaderName;
        
        THnSparse*              fHistGamma;        //! tracking photons histogram
        THnSparse*              fHistPi0;          //! tracking pi0 histogram
        THnSparse*              fHistPi0bkg;       //! tracking background pi0 histogram
        THnSparse*              fHistEvent;        //! tracking events histogram
        TH1F*                   fHistEventTrigger; //! tracks all event triggers

        AliPIDResponse*         fPIDResponse;    //!
        AliTRDonlineTrackMatching* fOnline;      //! For rating trdtrack to esdtrack matching
        //AliAnalysisTaskESDfilter*  esdfilter;
        AliConversionSelection* fConvsel;        //! to get 'make' pi0

        AliAnalysisTRDEfficiency(const AliAnalysisTRDEfficiency&); // not implemented
        AliAnalysisTRDEfficiency& operator=(const AliAnalysisTRDEfficiency&); // not implemented

        ClassDef(AliAnalysisTRDEfficiency, 2);
};

#endif
