/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Kaon2PC:
// Description: Analysis task to calculate Two-Particle Angular Correlation Functions of Neutral and Charged Kaons
// Author: Anjaly Sasikumar Menon
// (anjaly.sasikumar.menon@cern.ch)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef Kaon2PC_H
#define Kaon2PC_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "TMacro.h"

class AliAODEvent;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliPIDResponse;
class Kaon2PC : public AliAnalysisTaskSE
{
    public:
                                Kaon2PC();
                                Kaon2PC(const char *name);
        virtual                 ~Kaon2PC();

        virtual void            UserCreateOutputObjects();
        virtual void     SetTrackCuts(Double_t c1, Double_t c2, Double_t c3, Double_t c4);
        virtual void     SetV0TrackCuts(Double_t c5, Double_t c6, Double_t c7, Double_t c8, Double_t c9, Double_t c10, Double_t c11, Double_t c12, Double_t c13, Double_t c14, Double_t c15);
        Bool_t AcceptTrack(const AliAODTrack* Trk);
        Bool_t AcceptPosTrack(const AliAODTrack* Trk);
        Bool_t AcceptNegTrack(const AliAODTrack* Trk);
        Bool_t AcceptV0(const AliAODv0 *v0, Double_t *vertex);
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        AliEventCuts            fEventCuts; // event cuts

    private:
       AliAODEvent*            fAOD;           //! input event
       TList*                  fOutputList;    //! output list
       AliPIDResponse*         fPIDResponse;  //! pid response object’
       TH1F*                   fHistMK0;        //! dummy histogram
       TH1F*                   fHistMK0Cuts;        //! dummy histogram
       TH1F*                   fHistNV0;        //! dummy histogram
       TH1F*                   fHistKChPt;     //! dummy histogram
       TH1F*                   fHistK0Pt;     //! dummy histogram
       TH1F*                   fKChPhi;     //! dummy histogram
       TH1F*                   fK0Phi;     //! dummy histogram 
       TH1F*                   fHistKpPhi;     //! dummy histogram 
       TH1F*                   fHistKnPhi;     //! dummy histogram 
       TH1F*                   fVtx;        //! dummy histogram
       TH1F*                   fClusters;    //! dummy histogram
       TH2F*                   fPID;        //! dummy histogram
       TH2F*                   fPIDKa;      //! dummy histogram
       TH2F*                   fPIDKaon;        //! dummy histogram
       TH2F*                   fPIDK;       //! dummy histogram
       TH2F*                   fPIDKeCut;       //! dummy histogram
       TH2F*                   fPIDKpiCut;      //! dummy histogram
       TH3F*                   fHistK0PhiEta;       //! dummy histogram
       TH2F*                   fHistK0Phi;          //! dummy histogram
       TH2F*                   fHistK0Eta;      //! dummy histogram
       TH2F*                   fHistChPhi;      //! dummy histogram
       TH2F*                   fHistChEta;      //! dummy histogram
       TH2F*                   fHistChRap;      //! dummy histogram
       TH2F*                   fHistPosPhi;     //! dummy histogram
       TH2F*                   fHistPosEta;     //! dummy histogram
       TH3F*                   fHistPosPhiEta;  //! dummy histogram
       TH2F*                   fHistPosRap;     //! dummy histogram
       TH2F*                   fHistNegPhi;     //! dummy histogram
       TH2F*                   fHistNegEta;     //! dummy histogram
       TH3F*                   fHistNegPhiEta;  //! dummy histogram
       TH2F*                   fHistNegRap;     //! dummy histogram
       TH2F*                   fnsigmakaon;     //! dummy histogram
       TH2F*                   fNsigmaKaon;     //! dummy histogram
       TH2F*                   fNsigmaTOFK;     //! dummy histogram
       TH2F*                   fNsigmaTOFKaon;      //! dummy histogram
       TH2F*                   fNsigmaTPCTOFK;      //! dummy histogram
       TH1F*                   fHistNEvents;        //! dummy histogram
       TH1F*                   fHistEta;        //! dummy histogram
       TH1F*                   fHistDEta;       //! dummy histogram
       TH1F*                   fHistPhi;        //! dummy histogram
       TH1F*                   fHistDPhi;       //! dummy histogram
       TH1F*                   fHistMult;        //! dummy histogram
       TH1F*                   fHistCent;       //! dummy histogram
       TH1F*                   fHistSigCent;    //! dummy histogram
       TH2F*                   fHistInvCent;    //! dummy histogram
       TH2F*                   fHistCFPhi;      //! dummy histogram
       TH2F*                   fHistCFEta;      //! dummy histogram
       TH3F*                   fHistCF;     //! dummy histogram
       TH3F*                   fHistKChCh;      //! dummy histogram
       TH3F*                   fHistKPosKPos;       //! dummy histogram
       TH3F*                   fHistKPosKNeg;       //! dummy histogram
       TH3F*                   fHistKNegKNeg;       //! dummy histogram
       TH3F*                   fHistK0K0;       //! dummy histogram
       TH1F*                   fHistPPionPhi;       //! dummy histogram
       TH1F*                   fHistNPionPhi;       //! dummy histogram
       Double_t        PVz;
       Double_t        fLpTCut;        //not a pointer???
       Double_t        fUpTCut;
       Double_t        fEtaCut;
       Double_t        fSigCut;
       Double_t        fDecayLv0Cut;
       Double_t        fLpTv0Cut;
       Double_t        fUpTv0Cut;
       Double_t        fEtav0Cut;
       Double_t        fDcaPosToPrimVtxv0Cut;
       Double_t        fDcaNegToPrimVtxv0Cut;
       Double_t        fEtaPosv0Cut;
       Double_t        fEtaNegv0Cut;
       Double_t        fCosPACut;
       Double_t        fSigPosv0Cut;
       Double_t        fSigNegv0Cut;


        Kaon2PC(const Kaon2PC&); // not implemented
        Kaon2PC& operator=(const Kaon2PC&); // not implemented

        ClassDef(Kaon2PC, 1);
};

#endif
