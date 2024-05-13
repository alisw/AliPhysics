#ifndef AliAnalysisTaskJetAngCorrelations_cxx
#define AliAnalysisTaskJetAngCorrelations_cxx

#include "AliAODTrackSelection.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"

//___________________________________________________________________________________________________________________________________________
class AliAnalysisTaskJetAngCorrelations : public AliAnalysisTaskSE {

public:
    AliAnalysisTaskJetAngCorrelations();
    AliAnalysisTaskJetAngCorrelations(const char *name);
    virtual ~AliAnalysisTaskJetAngCorrelations();

    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    void RunData();
    void SetPreliminaryTrackCuts(Bool_t kinkDaugh, Bool_t TPCRefit, Bool_t ITSRefit, Float_t TPCnCr, Float_t TPCnCrnCl, Int_t ITSnCl, Float_t ITSchi2, Float_t TPCchi2, Double_t etaCut, Float_t DCAxy, Float_t DCAz) {
        fKinkDaugh  = kinkDaugh;
        fTPCRefit   = TPCRefit;
        fITSRefit   = ITSRefit;
        fTPCnCr     = TPCnCr;
        fTPCnCrnCl  = TPCnCrnCl;
        fITSnCl     = ITSnCl;
        fITSchi2    = ITSchi2;
        fTPCchi2    = TPCchi2;
        fEtaCut     = etaCut;
        fDCAxy      = DCAxy;
        fDCAz       = DCAz;
    };
    void SetJetCuts(Float_t jetRadius, Float_t minJetPt, Float_t minJetParticlePt, Float_t minLeadingPt) {
        fJetRadius          = jetRadius;
        fMinJetPt           = minJetPt;
        fMinJetParticlePt   = minJetParticlePt;
        fMinLeadingPt       = minLeadingPt;
    }
    void SetProtonCuts(Float_t protonDCAxy, Float_t protonDCAz, Float_t protonITSTOFpT, Float_t protonITSnsig, Float_t protonTPCnsigITS, Float_t protonTPCnsigTOF, Float_t protonTOFnsig) {
        fProtonDCAxy        = protonDCAxy;
        fProtonDCAz         = protonDCAz;
        fProtonITSTOFpT     = protonITSTOFpT;
        fProtonITSnsig      = protonITSnsig;
        fProtonTPCnsigITS   = protonTPCnsigTOF;
        fProtonTOFnsig      = protonTOFnsig;
    }
    void SetAntiprotonCuts(Float_t antiprotonDCAxy, Float_t antiprotonDCAz, Float_t antiprotonITSTOFpT, Float_t antiprotonITSnsig, Float_t antiprotonTPCnsigITS, Float_t antiprotonTPCnsigTOF, Float_t antiprotonTOFnsig) {
        fAntiprotonDCAxy        = antiprotonDCAxy;
        fAntiprotonDCAz         = antiprotonDCAz;
        fAntiprotonITSTOFpT     = antiprotonITSTOFpT;
        fAntiprotonITSnsig      = antiprotonITSnsig;
        fAntiprotonTPCnsigITS   = antiprotonTPCnsigTOF;
        fAntiprotonTOFnsig      = antiprotonTOFnsig;
    }
    void SetDeuteronCuts(Float_t deuteronDCAxy, Float_t deuteronDCAz, Float_t deuteronITSTOFpT, Float_t deuteronITSnsig, Float_t deuteronTPCnsigITS, Float_t deuteronTPCnsigTOF, Float_t deuteronTOFnsig) {
        fDeuteronDCAxy        = deuteronDCAxy;
        fDeuteronDCAz         = deuteronDCAz;
        fDeuteronITSTOFpT     = deuteronITSTOFpT;
        fDeuteronITSnsig      = deuteronITSnsig;
        fDeuteronTPCnsigITS   = deuteronTPCnsigTOF;
        fDeuteronTOFnsig      = deuteronTOFnsig;
    }
    void SetAntideuteronCuts(Float_t antideuteronDCAxy, Float_t antideuteronDCAz, Float_t antideuteronITSTOFpT, Float_t antideuteronITSnsig, Float_t antideuteronTPCnsigITS, Float_t antideuteronTPCnsigTOF, Float_t antideuteronTOFnsig) {
        fAntideuteronDCAxy        = antideuteronDCAxy;
        fAntideuteronDCAz         = antideuteronDCAz;
        fAntideuteronITSTOFpT     = antideuteronITSTOFpT;
        fAntideuteronITSnsig      = antideuteronITSnsig;
        fAntideuteronTPCnsigITS   = antideuteronTPCnsigTOF;
        fAntideuteronTOFnsig      = antideuteronTOFnsig;
    }

    //User Functions
    Bool_t   GetEvent                   ();
    Bool_t   IsProton                   (AliAODTrack  *track);
    Bool_t   IsAntiproton               (AliAODTrack  *track);
    Bool_t   IsDeuteron                 (AliAODTrack  *track);
    Bool_t   IsAntideuteron             (AliAODTrack  *track);
    Bool_t   PassedTrackSelection       (AliAODTrack  *track);
    Double_t GetDCAtoPrimaryVertex      (AliAODTrack  *track, Int_t index);

private:
    AliAODEvent     *fAODEvent;     //!
    AliPIDResponse  *fPIDResponse;  //!
    AliAODTrackSelection *fAODTrackCuts; //!
    AliESDtrackCuts *fESDTrackCuts; //!
    TList           *fOutputData;   //!
    TList           *fQAList;       //!

    //### Configurable Variables ###
    //Preliminary
    Bool_t          fKinkDaugh;     //
    Bool_t          fTPCRefit;      //
    Bool_t          fITSRefit;      //
    Float_t         fTPCnCr;        //
    Float_t         fTPCnCrnCl;     //
    Int_t           fITSnCl;        //
    Float_t         fITSchi2;       //
    Float_t         fTPCchi2;       //
    Double_t        fEtaCut;        //
    Float_t         fDCAxy;         //
    Float_t         fDCAz;          //
    //Jet
    Float_t         fJetRadius;        //
    Float_t         fMinJetPt;         //
    Float_t         fMinJetParticlePt; //
    Float_t         fMinLeadingPt;     //
    //Proton
    Float_t         fProtonDCAxy;      //
    Float_t         fProtonDCAz;       //
    Float_t         fProtonITSTOFpT;   //
    Float_t         fProtonITSnsig;    //
    Float_t         fProtonTPCnsigITS; //
    Float_t         fProtonTPCnsigTOF; //
    Float_t         fProtonTOFnsig;    //
    //Antiproton
    Float_t         fAntiprotonDCAxy;      //
    Float_t         fAntiprotonDCAz;       //
    Float_t         fAntiprotonITSTOFpT;   //
    Float_t         fAntiprotonITSnsig;    //
    Float_t         fAntiprotonTPCnsigITS; //
    Float_t         fAntiprotonTPCnsigTOF; //
    Float_t         fAntiprotonTOFnsig;    //
    //Deuteron
    Float_t         fDeuteronDCAxy;      //
    Float_t         fDeuteronDCAz;       //
    Float_t         fDeuteronITSTOFpT;   //
    Float_t         fDeuteronITSnsig;    //
    Float_t         fDeuteronTPCnsigITS; //
    Float_t         fDeuteronTPCnsigTOF; //
    Float_t         fDeuteronTOFnsig;    //
    //Antideuteron
    Float_t         fAntideuteronDCAxy;      //
    Float_t         fAntideuteronDCAz;       //
    Float_t         fAntideuteronITSTOFpT;   //
    Float_t         fAntideuteronITSnsig;    //
    Float_t         fAntideuteronTPCnsigITS; //
    Float_t         fAntideuteronTPCnsigTOF; //
    Float_t         fAntideuteronTOFnsig;    //

    //### Histograms ###
    TH1F *hNumberOfEvents;          //!
    TH1F *hNumberOfJets;            //!
    TH1F *hEventProtocol;           //!
    TH1F *hTrackProtocol;           //!
    TH1I *hNumberOfParticlesInJet;  //!
    TH1I *hNumberOfJetsInEvent;     //!

    TH1F *hEtaFullEvent;            //!
    TH1D *hJetRapidity;             //!

    TH1D *hPtFullEvent;             //!
    TH1D *hFullConePt;              //!
    TH1D *hPtJetParticle;           //!
    TH1D *hPtTotalJet;              //!
    TH1D *hPtSubtractedJet;         //!
    TH2D *hPtJetProtonDeuteron;     //!
    TH1D *hPtJetDeuteron;           //!
    TH1D *hPtDiff;                  //!

    TH1D *hJetConeRadius;           //!

    TH2F *hITSSignal;               //!
    TH2F *hTPCSignal;               //!
    TH2F *hTOFSignal;               //!
    TH2F *hTPCnsigma;               //!
    TH2F *hTPCnsigmaProton;         //!
    TH2F *hTPCnsigmaDeuteron;       //!
    TH2F *hTOFnsigma;               //!
    TH2F *hTOFnsigmaProton;         //!
    TH2F *hTOFnsigmaDeuteron;       //!
    TH2F *hITSnsigma;               //!
    TH2F *hITSnsigmaProton;         //!
    TH2F *hITSnsigmaDeuteron;       //!

    TH1F *hDCAxyFullEvent;          //!
    TH1F *hDCAzFullEvent;           //!
    TH1F *hDCAzJetParticle;         //!

    AliAnalysisTaskJetAngCorrelations(const AliAnalysisTaskJetAngCorrelations&);
    AliAnalysisTaskJetAngCorrelations& operator=(const AliAnalysisTaskJetAngCorrelations&);

    ClassDef(AliAnalysisTaskJetAngCorrelations, 1);

};
//___________________________________________________________________________________________________________________________________________
#endif
