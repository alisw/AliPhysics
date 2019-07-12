/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskFlowTPCEMCalRun2_H
#define AliAnalysisTaskFlowTPCEMCalRun2_H

//#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliQnCorrectionsManager.h"
#include "AliRDHFCuts.h"
#include "AliAODEvent.h"
#include "AliHFQnVectorHandler.h"
#include "AliAnalysisVertexingHF.h"


//class AliAODMCParticle;
//class AliAODMCHeader;

class AliAnalysisTaskFlowTPCEMCalRun2 : public AliAnalysisTaskSE  
{
	public:

		enum EventPlaneMeth{kTPC,kTPCVZERO,kVZERO,kVZEROA,kVZEROC,kPosTPCVZERO,kNegTPCVZERO};//Event Plane to be calculated in the task
		enum FlowMethod{kEP,kSP,kEvShapeEP,kEvShapeSP,kEPVsMass,kEvShapeEPVsMass};//Event Plane, Scalar Product or Event Shape Engeneering methods
		enum EventPlaneDet{kNone=-1,kFullTPC,kPosTPC,kNegTPC,kFullV0,kV0A,kV0C};
		enum q2Method{kq2TPC,kq2PosTPC,kq2NegTPC,kq2VZERO,kq2VZEROA,kq2VZEROC};

		AliAnalysisTaskFlowTPCEMCalRun2();
		//AliAnalysisTaskFlowTPCEMCalRun2(const char *name, AliRDHFCuts *rdCuts);
		AliAnalysisTaskFlowTPCEMCalRun2(const char *name);

		virtual                 ~AliAnalysisTaskFlowTPCEMCalRun2();

		virtual void            UserCreateOutputObjects();
		virtual void            UserExec(Option_t* option);
		virtual void            Terminate(Option_t* option);
		void            FindMother(AliAODMCParticle* part,int &label, int &pid, double &ptmom);

		//		void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
		//      	void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
		//　　　　　　　Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };
		void SetTenderTaskName(TString name) {fTenderTaskName = name;}

		void SetHarmonic(int harmonic)       {fHarmonic = harmonic;}
		void SetCalibraionType(int calibtype) {fCalibType = calibtype;}
		void SetNormMethod(int norm) {fNormMethod = norm;}
		void SetOADBFileName(TString filename) {fOADBFileName = filename;}

		void SetFlowMethod(int meth) {fFlowMethod = meth;}
		void SetEventPlaneDetector(int det) {fEvPlaneDet = det;}
		void SetSubEventDetectors(int detsubA, int detsubB) {fSubEvDetA = detsubA; fSubEvDetB = detsubB;}
		void SetqnMethod(int qnmethod) {fqnMeth = qnmethod;}
		void SetqnPercentileSelection(TString splinesfilepath) {fPercentileqn = true; fqnSplineFileName = splinesfilepath;}
		void SetTPCHalvesEtaGap(double etagap = 0.2) {fEtaGapInTPCHalves=etagap;}

		void GetTrkClsEtaPhiDiff(AliVTrack *t,AliVCluster *v,Double_t &phidiff, Double_t &etadiff);
		void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Double_t TrkPt, Double_t DCAxy, Int_t Bsign);
		void CheckMCgen(AliAODMCHeader* fMCheader,Double_t CutEta);
		void SetDCA(Double_t xy, Double_t z){DCAxy = xy, DCAz = z;};

		void SetMinCentrality(float mincentr=0.) {fMinCentr = mincentr;}
		void SetMaxCentrality(float maxcentr=100.) {fMaxCentr = maxcentr;}

		//virtual void LocalInit();

		Bool_t IsPdecay(int mpid);
		Bool_t IsDdecay(int mpid);
		Bool_t IsBdecay(int mpid);


	private:
		AliAODEvent* fAOD;           //! input event
		TList* fOutputList;    //! output list
		AliVEvent *fVevent;  //!event object
		AliPIDResponse *fpidResponse;

		//Bool_t fUseTender;

		double GetDeltaPsiSubInRange(double psi1, double psi2);

		void GetMainQnVectorInfo(double &mainPsin, double &mainMultQn, double mainQn[2], double &SubAPsin, double &SubAMultQn, double SubAQn[2], double &SunBPsin, double &SubBMultQn,  double SubBQn[2],AliHFQnVectorHandler* HFQnVectorHandler);

		bool LoadSplinesForqnPercentile();

		//cut parameter
		Double_t DCAxy, DCAz;
		Double_t massMin;
		Int_t Nch;

		//		TClonesArray *fTracks_tender;//Tender tracks

		TH1F* fHistPt; //! dummy histogram
		TH1F* fNevents; //no of events
		TH1F* fCent; //centrality
		TH1F* fVtxZ; //!Vertex z
		TH1F* fVtxX; //!Vertex x
		TH1F* fVtxY; //!Vertex y
		TH1F* fTrkPt;
		TH1F* fTrketa;    
		TH1F* fTrkphi;
		TH1F* fTrketa2;
		TH2F* fdEdx;
		TH1F* fTrkP;
		TH1F* fHistClustE;
		TH2F* fEMCClsEtaPhi;
		TH2F* fHistNoCells;
		TH2F* fHistCalCell;
		TH1F* fHistNCls;
		TH1F* fHistNClsE1;
		TH1F* fHistNClsE2;
		TH1F* fHistNClsE3;
		TH1F* fEMCTrkPt;
		TH1F* fEMCTrketa;
		TH1F* fEMCTrkphi;
		TH2F* fClsEtaPhiAftMatch;
		TH2F* fTPCnsig;
		TH1F* fHisteop;
		TH2F* fM20;
		TH2F* fM02;
		TH2F* fHistNsigEop;
		TH1F* fEMCTrkMatchPhi;
		TH1F* fEMCTrkMatchEta;
		TH2F* fHistelectron;
		TH2F* fHisthadron;
		TH1F* fInvmassLS;
		TH1F* fInvmassULS;
		TH2F* fInvmassLS_2D;
		TH2F* fInvmassULS_2D;

		//=====MC output=======
		AliAODMCParticle* fMCparticle;
		TH1F* fMCcheckMother;
		TClonesArray* fMCarray;
		AliAODMCParticle* fMCTrackpart;
		AliAODMCHeader* fMCheader;
		TH1F* fCheckEtaMC;
		TH2F* fHistMCorgPi0;
		TH2F* fHistMCorgEta;
		TH1F* fHistMCorgD;
		TH1F* fHistMCorgB;

		Int_t NpureMC;
		Int_t NpureMCproc;
		Int_t NembMCpi0;
		Int_t NembMCeta;

		TH2F* fPt_Btoe;
		TH1F* fNDB;
		TH1F* fHist_eff_HFE;
		TH1F* fHist_eff_TPC;

		TF1* fPi000;
		TF1* fEta000;

		TH1F* fHistPhoReco0;
		TH1F* fHistPhoReco1;
		TH1F* fHistPhoReco2;

		TH1F* fHistPhoPi0;
		TH1F* fHistPhoPi01;
		TH1F* fHistPhoEta;
		TH1F* fHistPhoEta1;

		TH1F* flPercentile;

		//AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask;
		//AliQnCorrectionsManager* fFlowQnVectorMgr;

		//TH1F* fmyEventPlane;
		//TH1F* fmyEventPlane2;
		//TH1F* fmyEventPlane3;
		TH2F* fEPcorV0AC;
		TH2F* fEPcorV0ATPC;
		TH2F* fEPcorV0CTPC;
		//TH1F* fTrkPhiEPFullTPC;
		//TH1F* fTrkPhiEPFullV0;
		//TH2F* fTrkPhiEPFullTPC_Pt;
		TH2F* fTrkPhiEPV0A_Pt;
		TH2F* fTrkPhiEPV0A_Pt_ele;
		//TH1F* fTrkPhiEP2;
		//TH2F* fTrkPhiEP_Pt;
		//TH2F* fTrkPhiEP2_Pt;
		//TH2F* fTrkPhicos2;
		//TH1F* fInplane;
		//TH1F* fOutplane;

		TH2F* fDCAxy_Pt_ele;
		TH2F* fDCAxy_Pt_had;
		TH2F* fDCAxy_Pt_LS;
		TH2F* fDCAxy_Pt_ULS;

		TH2F* fsubV0ACcos2;
		TH2F* fsubV0ATPCcos2;
		TH2F* fsubV0CTPCcos2;

		TH2F* fcorTrkPhicent_charge;
		TH2F* fcorTrkPhicent_ele;
		TH2F* fcorcentcos2_charge;
		TH2F* fcorcentInplane;
		TH2F* fcorcentOutplane;

		//TH2F* fQx1;
		//TH2F* fQy1;
		//TH2F* fQx2;
		//TH2F* fQy2;
		//TH2F* fQx3;
		//TH2F* fQy3;

		TH1F* fHistEvPlane[6]; 
		TH2F* fHistEvPlaneQncorr[6]; //histograms for EP angle vs. centrality 
		//TH3F* fHistEvPlaneQncorr[6]; //histograms for EP angle vs. centrality vs.qn
		TH2F* fHistqnVsCentrPercCalib[6]; //histograms for qn percentile calibration with qn fine binning
		TH2F* fHistEPResolVsCentrVsqn[3]; //histos needed to compute EP resolution vs. centrality vs. qn
		//TH3F* fHistEPResolVsCentrVsqn[3]; //histos needed to compute EP resolution vs. centrality vs. qn
		//AliRDHFCuts* fRDCuts;

		TString fTenderTaskName; //name of tender task needed to get the calibrated Qn vectors

		Int_t fHarmonic;
		Int_t fCalibType;
		Int_t fNormMethod;
		TString fOADBFileName;

		Int_t fFlowMethod;
		Int_t fEvPlaneDet;
		Int_t fSubEvDetA;
		Int_t fSubEvDetB;
		Int_t fqnMeth;

		bool fPercentileqn;
		TString fqnSplineFileName;		
		bool fLoadedSplines;
		double fEtaGapInTPCHalves;
		double fScalProdLimit;

		float fMinCentr;
		float fMaxCentr;

		TList* fqnSplinesList[6];

		AliAnalysisTaskFlowTPCEMCalRun2(const AliAnalysisTaskFlowTPCEMCalRun2&); // not implemented
		AliAnalysisTaskFlowTPCEMCalRun2& operator=(const AliAnalysisTaskFlowTPCEMCalRun2&); // not implemented

		//ClassDef(AliAnalysisTaskFlowTPCEMCalRun2, 1);
		ClassDef(AliAnalysisTaskFlowTPCEMCalRun2, 5);
};

#endif
