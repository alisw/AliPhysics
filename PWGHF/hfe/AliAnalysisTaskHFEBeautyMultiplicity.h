/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice                               */
/* $Id$                                                                   */

#ifndef AliAnalysisTaskHFEBeautyMultiplicity_H
#define AliAnalysisTaskHFEBeautyMultiplicity_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "TProfile.h"

class TH1F;
class AliAODEvent;
class AliMultSelection;
class AliAODMCParticle;
class AliAODMCHeader;


class AliAnalysisTaskHFEBeautyMultiplicity : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskHFEBeautyMultiplicity();
                                AliAnalysisTaskHFEBeautyMultiplicity(const char *name);
        virtual                 ~AliAnalysisTaskHFEBeautyMultiplicity();
        //---- called once at beginning of runtime ----//
        virtual void            UserCreateOutputObjects();
    
        //---- called for each event ----//
        virtual void            UserExec(Option_t* option);
        virtual void            GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
        //virtual void            FindPatches(Bool_t &hasfiredEG1, Bool_t &hasfiredEG2, Double_t emceta, Double_t emcphi);
        virtual void            SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC, Double_t TrkPt, Double_t DCAxy, Int_t Bsign);
    
        //---- called for end of analysis ----//
        virtual void            Terminate(Option_t* option);



        void SetAODAnalysis()   { SetBit(kAODanalysis, kTRUE); };
        void SetESDAnalysis()   { SetBit(kAODanalysis, kFALSE);};
        Bool_t IsAODanalysis()  const { return TestBit(kAODanalysis); };

        void SetCentralityEstimator(const char *estimator) {fCentralityEstimator = estimator; };

	void SetMultiProfileLHC16qt(TProfile *hprof) { fMultiEstimatorAvg = new TProfile(*hprof); }

	TProfile* GetEstimatorHistogram(const AliAODEvent *fAOD);
	Double_t Nref;
	void SetNref(Double_t nref) { Nref = nref; };

	Double_t GetCorrectedNtrackletsD(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult);
    
    
        //---- MC analysis ----//
        virtual void            FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom);
        Bool_t                  IsPdecay(int mpid);
        Bool_t                  IsDdecay(int mpid);
        Bool_t                  IsBdecay(int mpid);
        virtual void            CheckMCgen(AliAODMCHeader* fMCheader,Double_t CutEta);

	TProfile* GetEstimatorHistogramMC(const AliAODEvent *fAOD);
    
    private:
        enum{
            kAODanalysis = BIT(20),
            };

        AliAODEvent*        fAOD;                   // input event
        TList*              fOutputList;            // output list
        AliPIDResponse*     fpidResponse;           // pid response(+)
        AliVEvent*          fVevent;                // event object
        AliMultSelection*   fMultSelection;
        TString             fCentralityEstimator;
    
        //---- Tender ----//
        Bool_t              fUseTender;             // switch to add tender
        TClonesArray*       fTracks_tender;         // Tender tracks
        TClonesArray*       fCaloClusters_tender;   // Tender cluster
        
    
        TH1F* fNevents;
        TH1F* fCent;
        TH2F* fMult;
        TH1F* fVtxZ;                  // Zvertex position
        TH1F* fVtxZ_2;                // Zvertex position after event cut
        TH1F* fVtxX;                  // Xvertex position
        TH1F* fVtxY;                  // Yvertex position

        TH2F* fVtxCorrelation;        // Primary Zvertex vs. SPD Zvertex
	TH2F* fNcont;		      // Number of contribution
	
    	TH2F* fZvtx_Ntrklet;	      // Z vertex vs No. of tracklets
    	TH2F* fZvtx_Ntrklet_Corr;     // Z vertex vs No. of tracklets (corrected)


        TH2F* fEMCClsEtaPhi;          // EMCal Clusters (Eta and Phi)
        TH2F* fHistNCells;            // No of EMCal Cells in a cluster
        TH2F* fHistCalCells;          // EMCal cells in a cluster
        TH1F* fHistClustE;            // EMCal cluster energy distributioon
        TH1F* fHistNClsE;             // No of EMCal cluster in the event (All)
        TH1F* fHistNClsE1;            // No of EMCal cluster in the event (E > 0.1 GeV)
        TH1F* fHistNClsE2;            // No of EMCal cluster in the event (E > 0.2 GeV)
        TH1F* fHistNClsE3;            // No of EMCal cluster in the event (E > 0.5 GeV)
        TH1F* fAllTrkPt;              // pT distribution (All)
        TH1F* fEMCTrkPt;              // pT distribution (EMCal matched)
        
        TH1F* fAllTrkEta;             // eta distribution (All)
        TH1F* fEMCTrkEta;             // eta distribution (EMCal matched)
        TH1F* fAllTrkPhi;             // phi distribution (All)
        TH1F* fEMCTrkPhi;             // phi distribution (EMCal matched)
    	TH2F* fPhiEta;		      // Trk Phi vs. Eta
    	TH1F* fTPCCrossedRow;         // TPC CrossedRows

        TH2F* fdEdx;                  // dE/dx distribution (electron)
        TH2F* fTPCnsig;               // TPC Nsigma distribution (electron)
        TH2F* fTPCnsig_Pi;            // TPC Nsigma distribution (pion)
        TH2F* fTOFnsig;               // TOF Nsigma distribution (electron)
        TH2F* fITSnsig;               // ITS Nsigma distribution (electron)
        TH2F* fTPCnsigEta0;           // TPC Nsigma vs. eta (pT>2 geV/c)
        TH2F* fTPCnsigEta1;           // TPC Nsigma vs. eta (pT>3 geV/c)
        TH2F* fTPCnsigEta2;           // TPC Nsigma vs. eta (pT>5 geV/c)
        TH2F* fClsEtaPhiAftMatch;     // cluster eta&phi (after track matching)
        TH2F* fClsEtaPhiAftMatchEMCin;// cluster eta&phi (after track matching inside EMCal)
        TH2F* fClsEtaPhiAftMatchEMCout;// cluster eta&phi (after track matching outside EMCal)
        TH1F* fHistEMCTrkMatch_Eta;   // distance of EMCal cluster to its closest track (Eta)
        TH1F* fHistEMCTrkMatch_Phi;   // distance of EMCal cluster to its closest track (Phi)
        TH2F* fEMCTrkMatch_EtaPhi;
        TH2F* fEMCTrkMatch_EtaPhi_AfterCut;

	
    	TH1F* fTrkPt_2;		      // track pT (after track cut)
    	TH1F* fTrkEta_2;	      // track Eta (after track cut)
    	TH1F* fTrkPhi_2;	      // track Phi (after track cut)
    	TH2F* fdEdx_2;		      // dE/dx (after track cut)
    	TH2F* fTPCnsig_2;	      // TPC Nsigma (after track cut)
    	TH2F* fTOFnsig_2;	      // TOF Nsigma (after track cut)
    	TH2F* fITSnsig_2;	      // ITS Nsigma (after track cut)
    	TH1F* fTPCCrossedRow_2;       // TPC CrossedRows (after track cut)
    
        TH2F* fM02_1;
        TH2F* fM20_1;
        TH2F* fM02_2;
        TH2F* fM20_2;
        TH1F* fNtracks;
    	TH2F* fTrkEtaPhi_AfterCut;  // Track Eta vs. Phi (after cut)_
    
        TH1F* fHistEopAll;          // E/p (all)
        TH2F* fEopElectron1;        // pT vs E/p (electron)
        TH2F* fEopHadron1;          // pT vs E/p (electron)
    
        TH2F* fInvmassLS;           // Invariant mass vs. pT (Like-sign)
        TH2F* fInvmassULS;          // Invariant mass vs. pT (Unlike-sign)
    
        TH2F* fDCAxy_Ele_1;         // DCA Electron (DCA*charge*Bsign)
        TH2F* fDCAxy_Had_1;         // DCA Hadron (DCA*charge*Bsign)
        TH2F* fDCAxy_LS_1;          // DCA Like-sign (DCA*charge*Bsign)
        TH2F* fDCAxy_ULS_1;         // DCA Unlike-sign (DCA*charge*Bsign)

        TH2F* fDCAxy_Ele_2;         // DCA Electron (DCA*charge)
        TH2F* fDCAxy_Had_2;         // DCA Hadron (DCA*charge)
        TH2F* fDCAxy_LS_2;          // DCA Like-sign (DCA*charge)
        TH2F* fDCAxy_ULS_2;         // DCA Unlike-sign (DCA*charge)

        TH2F* fDCAxy_Ele_3;         // DCA Electron (DCA)
        TH2F* fDCAxy_Had_3;         // DCA Hadron (DCA)
        TH2F* fDCAxy_LS_3;          // DCA Like-sign (DCA)
        TH2F* fDCAxy_ULS_3;         // DCA Unlike-sign (DCA)

        TH2F* fDCAxy_Ele_4;

        TH1F* fEopElectron2;        // electron pT (tight)
        TH1F* fEopElectron3;        // Except photonic (invariant mass)
        TH1F* fEopHadron2;          // hadron pT (tight)

	TH2F* fHistConv_R;
    	TH2F* fElectronEtaPhi;      // eta vs. phi (electron)
    	TH2F* fHadronEtaPhi;	    // eta vs. phi (hadron)
        
    
    
    
        //---- MC output ----//
        TH1F*               fMCcheckMother;
        AliAODMCParticle*   fMCparticle;            // MC track particle
        AliAODMCParticle*   fMCTrackpart;           // MC track particle
        TClonesArray*       fMCarray;               // MC array
        AliAODMCHeader*     fMCheader;
    
        TH1F*               fHistPho_Reco0;
        TH1F*               fHistPho_Reco1;
        TH1F*               fHistPho_Reco2;
    
        Int_t               NembMCpi0;
        Int_t               NembMCeta;
        Int_t               NpureMCproc;
        Int_t               NpureMC;
        Int_t               Nch;
	Bool_t		    iBevt;
	TH1F*		    fNoB;
	TH1F*		    fNoD;
    
        TH1F*               fCheckEtaMC;
        TH2F*               fHistMCorg_Pi0;
        TH2F*               fHistMCorg_Eta;
        TH1F*               fHistMCorg_D;
        TH1F*               fHistMCorg_BD;
        TH1F*               fHistMCorg_B;
        TH1F*               fHistMCorg_Lc;
        TH2F*               fPt_Btoe;
        TH1F*               fHistPt_HFE_MC_B;
        TH1F*               fHistPt_HFE_MC_D;
        TH1F*               fHistPt_HFE_MC_Lc;

	TH2F*		    fDCAxy_MC_B;
	TH2F*		    fDCAxy_MC_D;
	TH2F*		    fDCAxy_MC_Dpm;
	TH2F*		    fDCAxy_MC_D0;
	TH2F*		    fDCAxy_MC_Ds;
	TH2F*		    fDCAxy_MC_Lc;

	TH2F*		    fDCAxy_MC_ele;
	TH2F*		    fDCAxy_MC_Phot;
    	
	TH1F*		    fHistPt_B_TrkCut0;
	TH1F*		    fHistPt_B_TrkCut1;
	TH1F*		    fHistPt_B_TrkCut2;
	TH1F*		    fHistPt_B_TrkCut3;
	TH1F*		    fHistPt_B_TrkCut4;
	TH1F*		    fHistPt_B_TrkCut5;
	TH1F*		    fHistPt_B_TrkCut6;
	TH1F*		    fHistPt_B_TrkCut7;
	TH1F*		    fHistPt_B_TrkCut8;
	TH1F*		    fHistPt_B_TrkCut9;

	TH1F*		    fHistPt_D_TrkCut0;
	TH1F*		    fHistPt_D_TrkCut1;
	TH1F*		    fHistPt_D_TrkCut2;
	TH1F*		    fHistPt_D_TrkCut3;
	TH1F*		    fHistPt_D_TrkCut4;
	TH1F*		    fHistPt_D_TrkCut5;
	TH1F*		    fHistPt_D_TrkCut6;
	TH1F*		    fHistPt_D_TrkCut7;
	TH1F*		    fHistPt_D_TrkCut8;
	TH1F*		    fHistPt_D_TrkCut9;
    
    
    

        AliAnalysisTaskHFEBeautyMultiplicity(const AliAnalysisTaskHFEBeautyMultiplicity&);                // not implemented
        AliAnalysisTaskHFEBeautyMultiplicity& operator = (const AliAnalysisTaskHFEBeautyMultiplicity&);   // not implemented

	TProfile*	fMultiEstimatorAvg;

        ClassDef(AliAnalysisTaskHFEBeautyMultiplicity, 1);
};

#endif
