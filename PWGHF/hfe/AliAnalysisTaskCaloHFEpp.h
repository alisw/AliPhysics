/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCaloHFEpp_H
#define AliAnalysisTaskCaloHFEpp_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
class AliAODMCParticle;
class AliAODMCHeader;
class AliMultSelection;

class AliAnalysisTaskCaloHFEpp : public AliAnalysisTaskSE  
{
	public:
		AliAnalysisTaskCaloHFEpp();
		AliAnalysisTaskCaloHFEpp(const char *name);
		//AliAnalysisTaskCaloHFEpp(const char *name,
		//								      Bool_t flagEG1,
		//						          Bool_t flagEG2,
		//						          Bool_t SetFlagClsTypeEMC,
		//						          Bool_t SetFlagClsTypeDCAL);
		virtual                 ~AliAnalysisTaskCaloHFEpp();

		virtual void            UserCreateOutputObjects();
		virtual void            UserExec(Option_t* option);
		virtual void            Terminate(Option_t* option);
		virtual void            SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC, Double_t TrkPt);
		virtual void            CheckMCgen(AliAODMCHeader* fMCheader);
		virtual void            FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom);
    virtual void            SetEtaRange(Int_t etarange){fetarange = etarange;};
    Bool_t GetEMCalTriggerEG1() { return fEMCEG1; };
    Bool_t GetEMCalTriggerEG2() { return fEMCEG2; };
    Bool_t GetDCalTriggerDG1()  { return fDCDG1; };
    Bool_t GetDCalTriggerDG2()  { return fDCDG2; };
    //void SetEMCalTriggerEG1(Bool_t flagTr1) { fEMCEG1=flagTr1; fEMCEG2=kFALSE;};
    //void SetEMCalTriggerEG2(Bool_t flagTr2) { fEMCEG2=flagTr2; fEMCEG1=kFALSE;};
		Bool_t                  IsPdecay(int mpid);
		Bool_t                  IsDdecay(int mpid);
		Bool_t                  IsBdecay(int mpid);
    Bool_t                  fEMCEG1;//EMcal Threshold EG1
    Bool_t                  fEMCEG2;//EMcal Threshold EG2
    Bool_t                  fDCDG1;//EMcal Threshold EG2
    Bool_t                  fDCDG2;//EMcal Threshold EG2
    void SetEG1(Bool_t flagEG1) { fEMCEG1= flagEG1;};
    void SetEG2(Bool_t flagEG2) { fEMCEG2= flagEG2;};
    void SetDG1(Bool_t flagDG1) { fDCDG1= flagDG1;};
    void SetDG2(Bool_t flagDG2) { fDCDG2= flagDG2;};
    Bool_t fFlagClsTypeEMC;//switch to select EMC clusters
    Bool_t fFlagClsTypeDCAL;//switch to select DCAL clusters
		void SetfFlagClsTypeEMC(Bool_t fEMC){fFlagClsTypeDCAL = fEMC;};
		void SetfFlagClsTypeDCAL(Bool_t fDCAL){fFlagClsTypeDCAL = fDCAL;};

	private:
		AliAODEvent*            fAOD;           //! input event
		TList*                  fOutputList;    //! output list
		AliVEvent   *fVevent;  //!event object
    TF1*                    fPi010;
    TF1*                    fEta010;
		TH1F*                   fHist_trackPt;        //! dummy histogram
		TH1F*                   fHistMatchPt;        
		TH1F*                   fHistSelectPt;        
		TH1F*                   fHist_ClustE;        //! dummy histogram
		TH1F*                   fHist_SelectClustE;
		TH1F*                   fHistMatchE;
		TH1F*                   fHistEta_track;        //! dummy histogram
		TH1F*                   fHistPhi_track;        //! dummy histogram
		TH1F*                   fHistEta_EMcal;        //! dummy histogram
		TH1F*                   fHistPhi_EMcal;        //! dummy histogram
		TH1F*                   fHist_VertexZ;        //! dummy histogram
		TH1F*                   fHist_Centrality;        //! dummy histogram
		TH1F*                   fNevents;
		TH1F*                   fTPCNcls;
		TH1F*                   fITSNcls;
		TH1F*                   fInvmassLS;
		TH1F*                   fInvmassULS;
		TH1F*                   fHistPt_Inc;
		TH1F*                   fEtadiff;
		TH1F*                   fPhidiff;
		TH1F*                   fEop_electron;
		TH1F*                   fEop_hadron;
		TH1F*                   fMCcheckMother;
		TH1F*                   fCheckEtaMC;
		TH1F*                   fHistMCorgD;
		TH1F*                   fHistMCorgB;
		TH1D*               		fHistPhoPi0;
		TH1D*               		fHistPhoPi1;
		TH1D*               		fHistPhoEta0;
		TH1D*               		fHistPhoEta1;
		TH1F*               		fHistPt_HFE_MC_D;
		TH1F*               		fHistPt_HFE_MC_B;
		//TH1F*               		fHist_eff_pretrack;
		//TH1F*               		fHist_eff_posttrack;
		TH1F*               		fHist_eff_HFE;
		TH1F*               		fHist_eff_match;
		TH1F*               		fHist_eff_TPC;

		TH2F*                   fHistScatter_EMcal;        //! dummy histogram
		TH2F*                   fHistScatter_EMcal_aftMatch;        //! dummy histogram
		TH2F*                   fHist_Mult;        //! dummy histogram
		TH2F*                   fdEdx;
		TH2F*                   fTPCnsig;
		TH2F*                   fHistNsigEop;
		TH2F*                   fM02;
		TH2F*                   fM20;
		TH2F*                   fM02_2;
		TH2F*                   fM20_2;
		TH2F*                   fEopPt_ele_loose;
		TH2F*                   fEopPt_ele_tight;
		TH2F*                   fEopPt_had;
		TH2F*                   fInv_pT_ULS;
		TH2F*                   fInv_pT_LS;
		TH2F*                   fHistMCorgPi0;
		TH2F*                   fHistMCorgEta;
		TH2F*                   fTrigMulti;
		TH1D*                   fHistPhoReco0;
		TH1D*                   fHistPhoReco1;
		TH1D*                   fHistPhoReco2;
                TH2F*                   fHistoNCells;

		AliAnalysisTaskCaloHFEpp(const AliAnalysisTaskCaloHFEpp&); // not implemented
		AliAnalysisTaskCaloHFEpp& operator=(const AliAnalysisTaskCaloHFEpp&); // not implemented
		AliMultSelection *fMultSelection;
		AliPIDResponse *fpidResponse; //!pid response  
		Bool_t       fUseTender;// switch to add tender
		TClonesArray  *fMCarray;//! MC array
		TClonesArray  *fTracks_tender;
		TClonesArray  *fCaloClusters_tender;
		AliAODMCParticle  *fMCparticle;
		AliAODMCHeader *fMCheader;
    Int_t NembMCpi0; // # of process in MC (no GEANT process)
    Int_t NembMCeta; // # of process in MC (no GEANT process)
    Int_t NpureMCproc; // # of process in MC (no GEANT process)
    Int_t fetarange;  


		ClassDef(AliAnalysisTaskCaloHFEpp, 1);
};

#endif
