/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisESDsTask_H
#define AliAnalysisESDsTask_H

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "AliV0ReaderV1.h"
#include "AliTRDonlineTrackMatching.h"
#include "AliAnalysisTaskESDfilter.h"


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
        virtual Bool_t          GetTrackCuts(AliESDtrack* track);
        virtual Double_t        GetSagitta(AliESDTrdTrack* trdtrack);
        virtual Double_t        GetRating(AliESDv0 *v0, AliESDtrack *track, AliESDTrdTrack *trdtrack);
        virtual void            Terminate(Option_t* option);

    private:
        AliESDEvent*            fESD;           //! input event
        TList*                  fOutputList;    //! output list
        TH1F*                   fHistPt;        //! dummy histogram

        AliV0ReaderV1*          v0reader;
        TH1F*                   fhBhqu;
        TH1F*                   fhpt;
        TH1F*                   fhsag1;
        TH1F*                   fhsag2;
        TH1F*                   fhlabel;
        
        THnSparse*              fhgevent1;
        THnSparse*              fhgevent2;
        THnSparse*              fhgevent3;
        THnSparse*              fhgevent4;
        THnSparse*              fhgevent5;
        THnSparse*              fhgevent6;
        THnSparse*              fhgevent7;
        THnSparse*              fhgevent8;
        THnSparse*              fhgevent9;
        
        THnSparse*              fhg;
        
        THnSparse*              fhgdghtr;
        THnSparse*              fhgtest;

        Int_t                   eventNumber = 0;
        TObjArray               *flst;    
        AliPIDResponse*         fPIDResponse;
        AliTRDonlineTrackMatching* online;
        AliAnalysisTaskESDfilter*  esdfilter;

        AliAnalysisTRDEfficiency(const AliAnalysisTRDEfficiency&); // not implemented
        AliAnalysisTRDEfficiency& operator=(const AliAnalysisTRDEfficiency&); // not implemented

        ClassDef(AliAnalysisTRDEfficiency, 1);
};

#endif
