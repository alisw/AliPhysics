/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTRDEfficiency_H
#define AliAnalysisTRDEfficiency_H

#include "AliAnalysisTaskSE.h"
//#include "AliKFConversionPhoton.h"
#include "THnSparse.h"

class AliAnalysisTRDEfficiency : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTRDEfficiency();
                                AliAnalysisTRDEfficiency(const char *name);
        virtual                 ~AliAnalysisTRDEfficiency();

        virtual void            UserCreateOutputObjects();
        virtual Bool_t          GetAODConversionGammas(AliAODEvent* fAOD);
        virtual void            UserExec(Option_t* option);
        virtual void            Photons(AliAODEvent* fAOD);
        //virtual TFile*          OpenDigitsFile(TString* inputfile, String* digfile, TString* opt);
        //virtual void          Tracks(AliAODEvent* fAOD);
        virtual void            Terminate(Option_t* option);
        
        
    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TH1F*                   fHistPt;        //! dummy histogram
        
        TFile*                  file;
        TList*                  fConversionGammas;
        
        TH1F*                   fhm1pt2;
        TH1F*                   fhNpttr;
        TH1F*                   fhNptun;
        TH1F*                   fhfA;
        //include filter bits
        
        // v0
        // all tracks
        TH1F*                   fhR;
        TH2F*                   fhRpt;
        TH1F*                   fhMv;
        TH1F*                   fhvpost;
        // hqu tracks
        TH1F*                   fhRhqu;
        TH2F*                   fhRpthqu;
        TH1F*                   fhMhqu;
        TH1F*                   fhvposthqu;

        //tracks
        TH1F*                   fhtxv;           // track->Xv()
        TH1F*                   fhtyv;
        TH1F*                   fhtzv;
        
        // gamma
        TH1F*                   fhgpt;
        TH1F*                   fhgpttrd;
        TH2F*                   fhgRpt;
        TH2F*                   fhgRpttrd;      // actually hqu but too late to change
        
        TH2F*                   fhgMinvM;       // comparing the GetPhotonMass and the GetInvMass functions
        TH1F*                   fhgR;
        TH2F*                   fhgptM;
        TH2F*                   fhgptMhqu;
        
        TH2F*                   fhgptQ;     // pt vs photon Quality
        TH2F*                   fhgptQhqu;
        TH2F*                   fhgetaphi;
        TH2F*                   fhgetaphihqu;
        
        TH2F*                   fhgxy;
        TH2F*                   fhgxyhqu;
        
        // v0 duagther particles
        TH1F*                   fhdn;
        TH1F*                   fhdpt;
        
        // n dimensional
        THnSparse*              fhna;
        THnSparse*              fhnp;
        THnSparse*              fhnhqu;
        
        // all of these events have a reconstructed photon
        THnSparse*              fhgevent2;      // events with 
        THnSparse*              fhgevent3;
        THnSparse*              fhgevent4;
        THnSparse*              fhgevent5;
        THnSparse*              fhgevent6;
        THnSparse*              fhgevent7;
        THnSparse*              fhgevent8;
        THnSparse*              fhgevent9;
        
        //event counter
        TH1F*                   fhgevent;    // event counter
        TH1F*                   fhevent;
        
        //track the events
        THnSparse*              fhtrckvnt;   // track event
        THnSparse*              fhtrckvnthqu;// track hqu event
        TH1*                    fhtrvnt;
        TH1*                    fhtrvnthqu;
        TList*                  lsttrckvnt;
        TList*                  lsttrckvnthqu;
        //AliConversionPhotonCuts fConversionCuts;
        
        TH2F*                   fhgetaphi1;
        TH1F*                   fhgR1;
        TH1F*                   fhgpt1;
        TH2F*                   fhgetaphi5;
        TH2F*                   fhgetaphi8;
        TH2F*                   fhgetaphi9;
        
        AliAnalysisTRDEfficiency(const AliAnalysisTRDEfficiency&); // not implemented
        AliAnalysisTRDEfficiency& operator=(const AliAnalysisTRDEfficiency&); // not implemented

        ClassDef(AliAnalysisTRDEfficiency, 1);
};

#endif
