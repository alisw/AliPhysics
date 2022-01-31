/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice                               */
/* $Id$                                                                   */
#ifndef AliAnalysisTaskHFEBeautyMultiplicity_H
#define AliAnalysisTaskHFEBeautyMultiplicity_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "TProfile.h"

class TH1F;
class TGraphErrors;
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
	virtual void		SelectTrack_DCA(Int_t itrack, AliAODTrack *track, Double_t TPC_CrossedRows, Double_t ITS_chi2, Double_t TPC_chi2NDF, Double_t dca_xy, Double_t dca_z, Double_t track_Eta, Bool_t &fFlagTrkSelect_forDCA);
        //---- called for end of analysis ----//
        virtual void            Terminate(Option_t* option);



        void SetAODAnalysis()   { SetBit(kAODanalysis, kTRUE); };
        void SetESDAnalysis()   { SetBit(kAODanalysis, kFALSE);};
        Bool_t IsAODanalysis()  const { return TestBit(kAODanalysis); };

        void SetCentralityEstimator(const char *estimator) {fCentralityEstimator = estimator; };


	void SetMultiProfileLHC16qt(TProfile *hprof)    { fMultiEstimatorAvg[0] = new TProfile(*hprof); }
	void SetMultiProfileLHC16qt_MC(TProfile *hprof) { fMultiEstimatorAvg[1] = new TProfile(*hprof); }

	TProfile* GetEstimatorHistogram(const AliAODEvent *fAOD, Bool_t iData);
	Double_t Nref;
	Double_t MinNtrklet;
	Double_t MaxNtrklet;

	void SetNref(Double_t nref) { Nref = nref; };
	void SetNtrkletMin(Double_t minNtrklet) { MinNtrklet = minNtrklet; };
	void SetNtrkletMax(Double_t maxNtrklet) { MaxNtrklet = maxNtrklet; };
	void SetWeightNtrklet(TH1D* hWeight){ fweightNtrkl = new TH1D(*hWeight); };

	void SetMCtype(Bool_t iMCtype){ iGPMC = iMCtype; };

	void SetWeightDmeson(TGraphErrors* WeightPt_Dmeson){ pTWeight_D = WeightPt_Dmeson; };
	void SetWeightLc(TGraphErrors* WeightPt_Lc){ pTWeight_Lc = WeightPt_Lc; };
	void SetWeightBmeson(TGraphErrors* WeightPt_Bmeson){ pTWeight_B = WeightPt_Bmeson; };
	//void SetWeightPi0(TF1* WeightPt_Pi0){ pTWeight_Pi0 = WeightPt_Pi0; };
	//void SetWeightEta(TF1* WeightPt_Eta){ pTWeight_Eta = WeightPt_Eta; };
	

	void SetTrackEta(Double_t etaMin, Double_t etaMax) {TrackEtaMin = etaMin, TrackEtaMax = etaMax;};
	void SetNsigma(Double_t NsigMin, Double_t NsigMax, Double_t HadNsig) {NsigmaMin = NsigMin, NsigmaMax = NsigMax, HadNsigma = HadNsig;};
	void SetM20(Double_t m20Min, Double_t m20Max) {M20Min = m20Min, M20Max = m20Max;};
	void SetEop(Double_t eopMin, Double_t eopMax) {EopMin = eopMin, EopMax = eopMax;};
	void SetDCA(Double_t xy, Double_t z) {DCAxy = xy, DCAz = z;};
	void SetTrackClust(Int_t TPC, Int_t ITS, Int_t Crossed, Double_t dEdx) {NTPCClust = TPC, NITSClust = ITS, NCrossedRow = Crossed, TPCdEdx = dEdx;};
	void SetDiff(Double_t diff) {EtaPhiDiff = diff;};
	void SetMass(Double_t invmass, Double_t photPt) {PhotInvMass = invmass, PhotMinPt = photPt;};
	



        //---- MC analysis ----//
        virtual void            FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom);
        Bool_t                  IsPdecay(int mpid);
        Bool_t                  IsDdecay(int mpid);
        Bool_t                  IsBdecay(int mpid);
        virtual void            CheckMCgen(AliAODMCHeader* fMCheader,Double_t CutEta);

    
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

	//---- Cut Parameter ----//
	
	Double_t TrackEtaMin;
	Double_t TrackEtaMax;
	Double_t NsigmaMin;
	Double_t NsigmaMax;
	Double_t HadNsigma;
	Double_t M20Min;
	Double_t M20Max;
	Double_t EopMin;
	Double_t EopMax;
	Int_t NTPCClust;
	Int_t NITSClust;
	Double_t TPCdEdx;
	Double_t DCAxy;
	Double_t DCAz;
	Int_t NCrossedRow;
	Double_t EtaPhiDiff;
	Double_t PhotInvMass;
	Double_t PhotMinPt;


    
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
        TH2F* fEMCTrkMatch_EtaPhi;    // EMCAL track matching Eta&Phi difference
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
        
    	TH1F* fHist_Tracklet;

	TH2F *fNsigma_Electron;
	TH2F *fNsigma_Hadron;

	TH1F *fHistPt_BeforePID;
    	TH2F *fdEdx_BeforePID;
    	TH2F *fTPCnsig_BeforePID;
    
    
        //---- MC output ----//
        TH1F*               fMCcheckMother;
        AliAODMCParticle*   fMCparticle;            // MC track particle
        AliAODMCParticle*   fMCTrackpart;           // MC track particle
        TClonesArray*       fMCarray;               // MC array
        AliAODMCHeader*     fMCheader;
    
        TH1F*	fHistPho_Reco0;
        TH1F*	fHistPho_Reco1;
        TH1F*	fHistPho_Reco2;
        TH1F*	fHistPho_Reco0_Pi0;
        TH1F*	fHistPho_Reco1_Pi0;
        TH1F*	fHistPho_Reco2_Pi0;
        TH1F*	fHistPho_Reco0_Eta;
        TH1F*	fHistPho_Reco1_Eta;
        TH1F*	fHistPho_Reco2_Eta;
    
        Int_t	NembMCpi0;
        Int_t	NembMCeta;
        Int_t	NpureMCproc;
        Int_t	NpureMC;
        Int_t	Nch;
	Int_t	Nmc;
	Bool_t	iGPMC;
	Bool_t	iBevt;
	TH1F*	fNoB;
	TH1F*	fNoD;
    
        TH1F*	fCheckEtaMC;
        TH2F*	fHistMCorg_Pi0;
        TH2F*	fHistMCorg_Eta;
        TH1F*	fHistMCorg_D;
        TH1F*	fHistMCorg_BD;
        TH1F*	fHistMCorg_B;
        TH1F*	fHistMCorg_Lc;
        TH2F*	fPt_Btoe;
        TH1F*	fHistPt_HFE_MC_B;
        TH1F*	fHistPt_HFE_MC_D;
        TH1F*	fHistPt_HFE_MC_Lc;

	TH2F*	fDCAxy_MC_B;
	TH2F*	fDCAxy_MC_B_weight;
	TH2F*	fDCAxy_MC_D;
	TH2F*	fDCAxy_MC_Dpm;
	TH2F*	fDCAxy_MC_Dpm_weight;
	TH2F*	fDCAxy_MC_D0;
	TH2F*	fDCAxy_MC_D0_weight;
	TH2F*	fDCAxy_MC_Ds;
	TH2F*	fDCAxy_MC_Ds_weight;
	TH2F*	fDCAxy_MC_Lc;
	TH2F*	fDCAxy_MC_Lc_weight;

	TH2F*	fDCA_B_total;
	TH2F*	fDCA_B_total_weight;
	TH2F*	fDCA_D_total;
	TH2F*	fDCA_Dpm_total;
	TH2F*	fDCA_Dpm_total_weight;
	TH2F*	fDCA_D0_total;
	TH2F*	fDCA_D0_total_weight;
	TH2F*	fDCA_Ds_total;
	TH2F*	fDCA_Ds_total_weight;
	TH2F*	fDCA_Lc_total;
	TH2F*	fDCA_Lc_total_weight;

	TH2F*	fDCAxy_MC_ele;
	TH2F*	fDCAxy_MC_Phot;
    	
	TH1F*	fHistPt_B_TrkCut;
	TH1F*	fHistPt_B_TrkCut0;
	TH1F*	fHistPt_B_TrkCut1;
	TH1F*	fHistPt_B_TrkCut2;
	TH1F*	fHistPt_B_TrkCut3;
	TH1F*	fHistPt_B_TrkCut4;
	TH1F*	fHistPt_B_TrkCut5;
	TH1F*	fHistPt_B_TrkCut6;
	TH1F*	fHistPt_B_TrkCut7;
	TH1F*	fHistPt_B_TrkCut8;
	TH1F*	fHistPt_B_TrkCut9;
	TH1F*	fHistPt_B_TrkCut10;
	TH1F*	fHistPt_B_TrkCut11;
	TH1F*	fHistPt_B_TrkCut12;
	
	TH1F*	fHistPt_Tracking_B;
	TH1F*	fHistPt_Tracking_B_2;
	TH1F*	fHistPt_EMCmatch_B;

	TH1F*	fHistPt_D_TrkCut;
	TH1F*	fHistPt_D_TrkCut0;
	TH1F*	fHistPt_D_TrkCut1;
	TH1F*	fHistPt_D_TrkCut2;
	TH1F*	fHistPt_D_TrkCut3;
	TH1F*	fHistPt_D_TrkCut4;
	TH1F*	fHistPt_D_TrkCut5;
	TH1F*	fHistPt_D_TrkCut6;
	TH1F*	fHistPt_D_TrkCut7;
	TH1F*	fHistPt_D_TrkCut8;
	TH1F*	fHistPt_D_TrkCut9;
	TH1F*	fHistPt_D_TrkCut10;
	TH1F*	fHistPt_D_TrkCut11;
	TH1F*	fHistPt_D_TrkCut12;

	TH1F*	fHistPt_Tracking_D;
	TH1F*	fHistPt_Tracking_D_2;
	TH1F*	fHistPt_EMCmatch_D;

	TH1F*	fHistEta_forTrackingEff;

	TH2F*	fNtrkletNch;
	TH2F*	fNtrkletNch_Corr;
	TH1F*	fNtrklet_Corr;

	TH2F*	fPhot_InvMass_vs_DCA;
	TH2F*	fPhot_InvMass_vs_DCA2;
	TH2F*	fPhot_InvMass_vs_DCA3;
	TH2F*	fPhot_InvMass_vs_DCA_data;
	TH2F*	fPhot_InvMass_vs_DCA_data2;
	TH2F*	fPhot_InvMass_vs_DCA_data3;

	TH1F*	fHistOrg_B;
	TH1F*	fHistOrg_D;
	TH1F*	fHistOrg_Dpm;
	TH1F*	fHistOrg_D0;
	TH1F*	fHistOrg_Ds;
	TH1F*	fHistOrg_Lc;

	TH1F*	fHistMCorg_Pi0_Enhance;
	TH1F*	fHistMCorg_Pi0_True;
	TH1F*	fHistMCorg_Eta_Enhance;
	TH1F*	fHistMCorg_Eta_True;

	TH2F* 	fHistPt_ele_vs_D;
	TH2F* 	fHistPt_ele_vs_BtoD;
	TH2F* 	fHistPt_ele_vs_B;
	TH2F* 	fHistPt_ele_vs_Lc;
    

        AliAnalysisTaskHFEBeautyMultiplicity(const AliAnalysisTaskHFEBeautyMultiplicity&);                // not implemented
        AliAnalysisTaskHFEBeautyMultiplicity& operator = (const AliAnalysisTaskHFEBeautyMultiplicity&);   // not implemented

	TProfile*	fMultiEstimatorAvg[2];
	TH1D*		fweightNtrkl;
	TGraphErrors* pTWeight_D;
	TGraphErrors* pTWeight_Lc;
	TGraphErrors* pTWeight_B;
	//TF1* pTWeight_B;
	TF1* pTWeight_Pi0;
	TF1* pTWeight_Eta;

        ClassDef(AliAnalysisTaskHFEBeautyMultiplicity, 1);
};

#endif
