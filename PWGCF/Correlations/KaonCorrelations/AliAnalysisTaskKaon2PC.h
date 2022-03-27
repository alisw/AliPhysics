 #ifndef AliAnalysisTaskKaon2PC_H
#define AliAnalysisTaskKaon2PC_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAODEvent;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliPIDResponse;
class AliAnalysisTaskKaon2PC : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskKaon2PC();
                                AliAnalysisTaskKaon2PC(const char *name);
        virtual                 ~AliAnalysisTaskKaon2PC();

        virtual void            UserCreateOutputObjects();
        virtual void     SetTrackCuts(Double_t c1, Double_t c2, Double_t c3, Double_t c4);
        virtual void     SetV0TrackCuts(Double_t c5, Double_t c6, Double_t c7, Double_t c8, Double_t c9, Double_t c10, Double_t c11, Double_t c12, Double_t c13, Double_t c14, Double_t c15);
        Bool_t AcceptTrack(const AliAODTrack* Trk);
        Bool_t AcceptV0(const AliAODv0 *v0, Double_t *vertex);
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        AliEventCuts            fEventCuts; // event cuts

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        AliPIDResponse*         fPIDResponse;  //! pid response object’
        TH1F*                   fHistMK0;
        TH1F*                   fHistPt;
        TH2F*                   fPID;
        TH2F*                   fPIDKa;
        TH2F*                   fPIDKaon;
        TH2F*                   fPIDK;
        TH2F*                   fHistK0Phi;
        TH2F*                   fHistK0Eta;
        TH2F*                   fHistChPhi;
        TH2F*                   fHistChEta;
        TH2F*                   fHistChRap;
        TH2F*                   fHistPosPhi;
        TH2F*                   fHistPosEta;
        TH2F*                   fHistPosRap;
        TH2F*                   fHistNegPhi;
        TH2F*                   fHistNegEta;
        TH2F*                   fHistNegRap;
        TH2F*                   fnsigmakaon;
        TH2F*                   fNsigmaKaon;
        TH1F*                   fHistNEvents;        //! dummy histogram
        TH1F*                   fHistEta;        //! dummy histogram
        TH1F*                   fHistDEta;
        TH1F*                   fHistPhi;        //! dummy histogram
        TH1F*                   fHistDPhi;
        TH1F*                   fHistMult;        //! dummy histogram
        TH1F*                   fHistCent;
        TH2F*                   fHistCFPhi;
        TH2F*                   fHistCFPhiCuts;
        TH2F*                   fHistCFPhiLCuts;
        TH2F*                   fHistCFEta;
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

        AliAnalysisTaskKaon2PC(const AliAnalysisTaskKaon2PC&); // not implemented
        AliAnalysisTaskKaon2PC& operator=(const AliAnalysisTaskKaon2PC&); // not implemented

        ClassDef(AliAnalysisTaskKaon2PC, 1);
};

#endif
