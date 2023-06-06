/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskKaon2PC:
// Description: Analysis task to calculate Two-Particle Angular Correlation Functions of Neutral and Charged Kaons
// Author: Anjaly Sasikumar Menon
// (anjaly.sasikumar.menon@cern.ch)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskKaon2PC_H
#define AliAnalysisTaskKaon2PC_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "TMacro.h"
#include "AliEventPoolManager.h" //required in the header

class AliAODEvent;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliPIDResponse;
class AliMCParticle;
class THnSparse;
class AliAODv0;
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
        Bool_t AcceptPosTrack(const AliAODTrack* Trk);
        Bool_t AcceptNegTrack(const AliAODTrack* Trk);
        Double_t Beta(const AliAODTrack *track);
        Bool_t AcceptV0(const AliAODv0 *v0, Double_t *vertex);
        Bool_t SelectK0TracksMC(AliMCParticle *mcTrack);
        Bool_t SelectKPosTracksMC(AliMCParticle *mcTrack);
        Bool_t SelectKNegTracksMC(AliMCParticle *mcTrack);
        Bool_t SelectKchTracksMC(AliMCParticle *mcTrack);
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void            Fill2DHist(Double_t DPhi, Double_t DEta, TH2F* hist);
        virtual void            FillDPhiHist(Double_t DPhi, TH2F* hist, Double_t fWeight);
        virtual void            Fill2DHistMCTruth(Double_t DPhi, Double_t DEta, TH2F* hist);
        virtual void            Fill3DHist(Double_t DPhi, Double_t DEta, Double_t PVz, TH3F* hist);
        AliEventCuts            fEventCuts; // event cuts

        void                    SetMCRead(Bool_t flag) {fAnalysisMC = flag;}
        void                    SetFilterBit(Int_t filterbit) {fBit = filterbit;}
        void                    SetCentLimit(Double_t CentMin, Double_t CentMax) {fCentMin = CentMin; fCentMax = CentMax; }
        void                    SetPtLimits(Double_t ptmin, Double_t ptmax) { fLpTCut = ptmin; fUpTCut=ptmax; }
        //void                  SetEtaLimit(Double_t etalimit) { fEta = etalimit; }
        //mixing
        //void                  SetNofSamples(Int_t n) { fNOfSamples = n; } //sampling setter
        void                    SetCentBinsForMixing(Int_t nofBins, std::vector<Double_t> bins) { fNCentBins = nofBins; fCentBins = bins; }
        
    private:
       void RunData();
       void RunMC();

       AliAODEvent*            fAOD;           //! input event
       TList*                  fOutputList;    //! output list
       AliPIDResponse*         fPIDResponse;   //! pid response object’

       TH1F*                   fEnergy;        //! dummy histogram
       TH1F*                   fEnergyCuts;    //! dummy histogram
       TH2F*                   fPID;            //! dummy histogram
       TH2F*                   fPIDKaon;        //! dummy histogram
       TH2F*                   fPIDK;           //! dummy histogram
       TH2F*                   fPIDKpiCut;      //! dummy histogram
       TH2F*                   fnsigmakaon;     //! dummy histogram
       TH2F*                   fNsigmaKaon;     //! dummy histogram
       TH2F*                   fNsigmaTOFK;     //! dummy histogram
       TH2F*                   fNsigmaTOFKaon;  //! dummy histogram
       TH2F*                   fNsigmaTPCTOFK;  //! dummy histogram
       TH2F*                   fHistTOFKch;     //! dummy histogram

       TH1F*                   fVtx;            //! dummy histogram
       TH1F*                   fClusters;       //! dummy histogram
       TH1F*                   fHistPVz;        //! dummy histogram
       TH1F*                   fHistNEvents;    //! dummy histogram
       TH1F*                   fHistNV0;        //! dummy histogram
       TH1F*                   fHistEta;        //! dummy histogram
       TH1F*                   fHistDEta;       //! dummy histogram
       TH1F*                   fHistPhi;        //! dummy histogram
       TH1F*                   fHistDPhi;       //! dummy histogram
       TH1F*                   fHistMult;       //! dummy histogram
       TH1F*                   fHistCent;       //! dummy histogram
       TH2F*                   fHistTPCTracksVsClusters; //! dummy histogram

       TH1F*                   fHistMK0;       //! dummy histogram
       TH1F*                   fHistMK0Cuts;   //! dummy histogram
       TH1F*                   fHistKChPt;     //! dummy histogram
       TH1F*                   fHistK0Pt;      //! dummy histogram
       TH1F*                   fHistKChPhi;    //! dummy histogram
       TH1F*                   fHistK0Phi;     //! dummy histogram 
       TH1F*                   fHistKpPhi;     //! dummy histogram 
       TH1F*                   fHistKnPhi;     //! dummy histogram 
       TH1F*                   fHistPPionPhi;  //! dummy histogram
       TH1F*                   fHistNPionPhi;  //! dummy histogram
       
       TH2F*                   f2DHistK0Phi;      //! dummy histogram
       TH2F*                   f2DHistK0Eta;      //! dummy histogram
       TH2F*                   f2DHistChPhi;      //! dummy histogram
       TH2F*                   f2DHistChEta;      //! dummy histogram
       TH2F*                   f2DHistChRap;      //! dummy histogram
       TH2F*                   f2DHistPosPhi;     //! dummy histogram
       TH2F*                   f2DHistPosEta;     //! dummy histogram
       TH2F*                   f2DHistPosRap;     //! dummy histogram
       TH2F*                   f2DHistNegPhi;     //! dummy histogram
       TH2F*                   f2DHistNegEta;     //! dummy histogram
       TH2F*                   f2DHistNegRap;     //! dummy histogram

       TH2F*                   fHistK0PhiEta;    //! dummy histogram
       TH2F*                   fHistPosPhiEta;   //! dummy histogram
       TH2F*                   fHistNegPhiEta;   //! dummy histogram

       TH2F*                   fHistCFPhi;       //! dummy histogram
       TH2F*                   fHistCFEta;       //! dummy histogram
       TH2F*                   fHistKChKChPhi;   //! dummy histogram
       TH2F*                   fHistKPosKNegPhi; //! dummy histogram
       TH2F*                   fHistCF;          //! dummy histogram
       TH3F*                   fHistCFz;         //! dummy histogram
       TH2F*                   fHistKChKCh;      //! dummy histogram
       TH2F*                   fHistKPosKNeg;    //! dummy histogram
       TH3F*                   fHistKPosKNegz;   //! dummy histogram

       TH1F*                   hPt;
       TH1F*                   hPt_kPos;
       TH2F*                   fHistCF_Bg;
       TH2F*                   fHistCF_KpKn_Bg;
       TH3F*                   fHistCF_Bgz;
       TH3F*                   fHistCF_KpKn_Bgz;

       AliMCEvent*             fmcEvent;
       THnSparse*              fMCK0;
       THnSparse*              fMCKpos;
       THnSparse*              fMCKneg;
       TH2F*                   fHistKpKnMC;
       TH2F*                   fHistK0KchMC;
       TH1D*                   fHistGenMultiplicity;

       Double_t        fPV[3];

       Bool_t                  fAnalysisMC; // enable MC study
       Bool_t                  fRejectEventPileUp; // enable to use Pile-up cuts
       Bool_t                  fKpKnCorr;
       Bool_t                  fK0KchCorr;

       Double_t        PVx;
       Double_t        PVy;
       Double_t        PVz;
       Double_t        fBit;
       Double_t        fCentMin;
       Double_t        fCentMax;
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


       // mixing

        TObjArray*              fSelectedKCh; //!
        TObjArray*              fSelectedK0s; //!
        TObjArray*              fSelectedKpos; //!
        TObjArray*              fSelectedKneg; //!
        AliEventPoolManager*    fPoolMgr;  //!  event pool manager for Event Mixing
        // Int_t                   fNOfSamples;
        // std::vector<Double_t>   fsampleBins; //sampling

        Int_t                   fPoolMaxNEvents;   // maximum number of events in the pool
        Int_t                   fPoolMinNTracks;   // minimum number of tracks in the pool
        Int_t                   fMinEventsToMix;   // minimum number of events for mixing
        Int_t                   fNzVtxBins; // number of PV z bins
        Int_t                   fNCentBins; // number of centrality bins
        std::vector<Double_t>   fzVtxBins;
        std::vector<Double_t>   fCentBins;
        //Double_t                fMergingCut; // [0.02] cut for track spliting/merging
        //Int_t                   fSampleIndex;
        
        //TString                 fSystematicsFlag; // ""
        //Int_t                   fEtaPolarity; // trig particles in one half of detector, -1,0,1

        // Double_t                fPVzCut;
        // Double_t                fCutDCAz;
        // Double_t                fCutDCAxySigma;
        // Double_t                fCutTPCchi2pCl;
        // Double_t                fTPCclMin;


       AliAnalysisTaskKaon2PC(const AliAnalysisTaskKaon2PC&); // not implemented
       AliAnalysisTaskKaon2PC& operator=(const AliAnalysisTaskKaon2PC&); // not implemented

       ClassDef(AliAnalysisTaskKaon2PC, 1);
};

#endif