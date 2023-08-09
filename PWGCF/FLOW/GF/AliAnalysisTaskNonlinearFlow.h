#ifndef ALIANALYSISTASKNONLINEARFLOW_H
#define ALIANALYSISTASKNONLINEARFLOW_H
#include "AliAnalysisTaskSE.h"
#include "AliGFWMCuts.h"
#include "AliGFWNFCuts.h"
#include "AliGFWWeights.h"
#include "CorrelationCalculator.h"
#include "AliEventCuts.h"
#include <TComplex.h>

#include <TObject.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TComplex.h>
#include <TBits.h>
#include <TRandom3.h>
// AliRoot includes
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"
#include "AliGenHijingEventHeader.h"

#ifndef __CINT__
// ROOT includes
// # include <TList.h>
// # include <TH1.h>
// # include <TH2.h>
// # include <TH3.h>
// # include <TProfile.h>
// # include <TComplex.h>
// # include <TBits.h>
# include <TRandom3.h>
// AliRoot includes
// # include "AliESDEvent.h"
// # include "AliAODEvent.h"
// # include "AliVEvent.h"
// # include "AliVTrack.h"
// # include "AliVVertex.h"
// # include "AliAnalysisFilter.h"
// # include "AliESDtrackCuts.h"
#else
#endif
class TList;
class TF1;
class TH1;
class TH2;
class TH3F;
class TH1F;
class TH2F;
class TH3D;
class TProfile;
class TProfile2D;
class TComplex;
class AliESDEvent;
class AliAODEvent;
class AliVEvent;
class AliVTrack;
class AliVVertex;
class AliESDtrackCuts;
class AliAODITSsaTrackCuts;
class AliInputEventHandler;

#include <THnSparse.h>

class PhysicsProfile : public TObject {
	public:
		PhysicsProfile();
		PhysicsProfile(const PhysicsProfile&);
		// Physics profiles
		TProfile*	 fChsc4242;		             	//! SC(4,2)
		TProfile*	 fChsc4242_Gap0;			//! SC(4,2) |#Delta#eta| > 0.0
		TProfile*	 fChsc4242_Gap2;			//! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap4;                        //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap6;                        //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap8;                        //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap10;                       //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*	 fChsc4242_3sub;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subMMLRA;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subMMLRB;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subLLMRA;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subLLMRB;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subRRMLA;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4242_3subRRMLB;			//! SC(4,2)_A 3subevent method
		TProfile*	 fChsc4224_3sub;			//! SC(4,2)_B 3subevent method
		TProfile*	 fChsc4242_3subGap2;			//! SC(4,2)_A 3subevent method |#Delta#eta| > 0.2
		TProfile*	 fChsc4224_3subGap2;			//! SC(4,2)_B 3subevent method |#Delta#eta| > 0.2
		TProfile*	 fChsc3232;				//! SC(3,2)
		TProfile*	 fChsc3232_Gap0;			//! SC(3,2) |#Delta#eta| > 0.0
		TProfile*	 fChsc3232_Gap2;			//! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap4;                        //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap6;                        //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap8;                        //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap10;                       //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*     	 fChsc3232_3sub;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subMMLRA;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subMMLRB;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subLLMRA;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subLLMRB;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subRRMLA;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3232_3subRRMLB;			//! SC(3,2)_A 3subevent method
		TProfile*     	 fChsc3223_3sub;			//! SC(3,2)_B 3subevent method
		TProfile*     	 fChsc3232_3subGap2;			//! SC(3,2)_A 3subevent method |#Delta#eta| > 0.2
    TProfile*     	 fChsc3223_3subGap2;			//! SC(3,2)_B 3subevent method |#Delta#eta| > 0.2

		// Standard correlation profiles for different harmonics
		TProfile*	 fChc422;          //!
		TProfile*        fChc532;          //!
		TProfile*        fChc422_Gap0A;    //!
		TProfile*        fChc422_Gap0B;    //!
		TProfile*        fChc532_Gap0A;    //!
		TProfile*        fChc532_Gap0B;    //!
		TProfile*        fChc422_Gap2A;    //!
		TProfile*        fChc422_Gap2B;    //!
		TProfile*        fChc532_Gap2A;    //!
		TProfile*        fChc532_Gap2B;    //!
		TProfile*        fChc422_Gap4A;    //!
		TProfile*        fChc422_Gap4B;    //!
		TProfile*        fChc532_Gap4A;    //!
		TProfile*        fChc532_Gap4B;    //!
		TProfile*        fChc422_Gap6A;    //!
		TProfile*        fChc422_Gap6B;    //!
		TProfile*        fChc532_Gap6A;    //!
		TProfile*        fChc532_Gap6B;    //!
		TProfile*        fChc422_Gap8A;    //!
		TProfile*        fChc422_Gap8B;    //!
		TProfile*        fChc532_Gap8A;    //!
		TProfile*        fChc532_Gap8B;    //!
		TProfile*        fChc422_Gap10A;    //!
		TProfile*        fChc422_Gap10B;    //!
		TProfile*        fChc532_Gap10A;    //!
		TProfile*        fChc532_Gap10B;    //!
		TProfile*        fChc422_3subL;    //!
		TProfile*        fChc422_3subM;    //!
		TProfile*        fChc422_3subR;    //!
		TProfile*        fChc532_3subLA;    //!
		TProfile*        fChc532_3subLB;    //!
		TProfile*        fChc532_3subMA;    //!
		TProfile*        fChc532_3subMB;    //!
		TProfile*        fChc532_3subRA;    //!
		TProfile*        fChc532_3subRB;    //!

    TProfile*  fMeanPt;         //! Average of Pt
    TProfile*  fc22w;           //! vn^2 with event weight
    TProfile*  fPcc;            //! v2^2-pt
    TProfile*  fc22nw;           //! vn^2 without event weight
    TProfile*  fc24nw;           //! vn^4 without event weight
    TProfile*  fPtVariancea;                  //! for variance of pt
    TProfile*  fPtVarianceb;                  //! for variance of pt
    TProfile*  fPtVariancec;                  //! for variance of pt

		TProfile*	 fChcn2[6]; 			//! <<2>> in unit bins of Ntrks
		TProfile*    	 fChcn2_Gap0[6];  		//! <<2>> |#Delta#eta| > 0.0
		TProfile*	 fChcn2_Gap2[6];  		//! <<2>> |#Delta#eta| > 0.2
		TProfile*	 fChcn2_Gap4[6];  		//! <<2>> |#Delta#eta| > 0.4
		TProfile*        fChcn2_Gap6[6];                //! <<2>> |#Delta#eta| > 0.6
		TProfile*	 fChcn2_Gap8[6];  		//! <<2>> |#Delta#eta| > 0.8
		TProfile*	 fChcn2_Gap10[6];  		//! <<2>> |#Delta#eta| > 1.0
		TProfile*	 fChcn2_Gap14[6];  	        //! <<2>> |#Delta#eta| > 1.4
		TProfile*	 fChcn2_Gap16[6];  	        //! <<2>> |#Delta#eta| > 1.6
		TProfile*	 fChcn2_Gap18[6];  	        //! <<2>> |#Delta#eta| > 1.8

		TProfile*	 fChcn2_3subLM[6];  		//! <<2>> left vs. middle subevent
		TProfile*	 fChcn2_3subRM[6];  		//! <<2>> middle vs. right subevent
		TProfile*	 fChcn2_3subLR[6];  		//! <<2>> middle vs. right subevent
		TProfile*	 fChcn2_3subGap2LM[6];          //! <<2>> left vs. middle subevent |#Delta#eta| > 0.2
		TProfile*	 fChcn2_3subGap2RM[6];          //! <<2>> middle vs. right subevent |#Delta#eta| > 0.2

		TProfile*	 fChcn4[6];  			//! <<4>> in unit bins of Ntrks
		TProfile*	 fChcn4_Gap0[6];  		//! <<4>> |#Delta#eta| > 0.0
		TProfile*	 fChcn4_Gap2[6];  		//! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap4[6];                //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap6[6];                //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap8[6];                //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap10[6];               //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_3subMMLR[6];            //! <<4>> 3subevent method
		TProfile*        fChcn4_3subLLMR[6];            //! <<4>> 3subevent method
		TProfile*        fChcn4_3subRRML[6];            //! <<4>> 3subevent method
		TProfile*        fChcn4_3subGap2[6];            //! <<4>> 3subevent method |#Delta#eta| > 0.2

                TProfile*	 fChcn6[6];  			//! <<6>> in unit bins of Ntrks
		TProfile*	 fChcn6_Gap0[6];  		//! <<6>> |#Delta#eta| > 0.0
                TProfile*	 fChcn8[6];  			//! <<8>> in unit bins of Ntrks
		TProfile*	 fChcn8_Gap0[6];  		//! <<8>> |#Delta#eta| > 0.0

	private:
		ClassDef(PhysicsProfile, 8);    //Analysis task
};

class AliAnalysisTaskNonlinearFlow : public AliAnalysisTaskSE {
	public:

		enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kCharUnidentified, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest; NB: kUknown last as counter

		enum  NonflowSupress {knStandard = 1 << 0, kn0Gap = 1 << 1, knLargeGap = 1 << 2, knThreeSub = 1 << 3, knGapScan = 1 << 4}; // list of nonflow supression method



		// const unsigned int usev2345flag = 1 << 0;
		// const unsigned int usev678flag = 1 << 1;

		AliAnalysisTaskNonlinearFlow();
		AliAnalysisTaskNonlinearFlow(const char *name);
		AliAnalysisTaskNonlinearFlow(const char *name, int NUA, int NUE, TString fPeriod);

		virtual ~AliAnalysisTaskNonlinearFlow();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t* option);
                virtual void   NotifyRun();
		virtual void   Terminate(Option_t* );

		virtual void   SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void   SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void   SetMinPt(Double_t minPt){fMinPt = minPt;}
		virtual void   SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
		virtual void   SetTrigger(Int_t trig){fTrigger = trig;}
		virtual void   SetNUEFlag(Bool_t NUE){fNUE = NUE;}
		virtual void   SetNUA(Bool_t NUA){fNUA = NUA;}
		virtual void   SetIsMC(Bool_t isMC){fIsMC = isMC;}
		virtual void   SetNtrksName(TString ntrksname){fNtrksName = ntrksname;}
		virtual void   SetPeriod(TString period) { fPeriod = period; }
		virtual void   SetSystFlag(int flag) { fCurrSystFlag = flag; }
		virtual int    GetSystFlag() { return fCurrSystFlag; }
		virtual void   SetSpringMode(bool flag = true) { fSpringMode = flag; }
		virtual void   SetLowMultiplicityMode(bool flag = true) {fLowMultiplicityMode = flag;}
		virtual void   SetAdditionalTPCPileupCuts(bool flag = true) {fAddTPCPileupCuts = flag;}
		virtual void   SetESDvsTPConlyLinearCut(double cut = 15000) {fESDvsTPConlyLinearCut = cut;}
    virtual void   SetUseOutOfBunchPileupCut(double flag = true) {fUseOutOfBunchPileupCut = flag;}
		virtual void   SetUseCorrectedNTracks(bool flag = true) {fUseCorrectedNTracks = flag;}
    virtual void   SetBinningFactor(double factor = 1.0) {binning_factor = factor;}
		virtual void   SetUseNarrowBin(bool flag = true) {fUseNarrowBin = flag;}
		virtual void   SetExtremeEfficiency(int flag = 0) {fExtremeEfficiency = flag;}
		virtual void   SetTPCchi2perCluster(double fchi2 = 4) {fTPCchi2perCluster = fchi2;}
		virtual void   SetUseAdditionalDCACut(double flag = true) {fUseAdditionalDCACut = flag;}
		virtual void   SetUseDefaultWeight(double flag = true) {fUseDefaultWeight = flag;}
    virtual void   SetUseLikeSign(int sign = 0) {bUseLikeSign = sign; iSign = sign;}
    virtual void   SetExtendV0MAcceptance(double flag = true) {fExtendV0MAcceptance = flag;}
    virtual void   SetV0MRatioCut(double ratio=5) {fV0MRatioCut = ratio;}
    virtual void   SetEtaGap3Sub(Double_t feta1 = 0.4, Double_t feta2 = 0.4) {fEtaGap3Sub1 = feta1; fEtaGap3Sub2 = feta2;}
		virtual void   SetCentralityCut(Double_t cent = 100) {fCentralityCut = cent;}
    virtual void   SetOnTheFly(Bool_t flag=false) {fOnTheFly = flag;} 
		// unsigned fgFlowHarmonics = 0;        calculate v2, v3, v4, v5
		// unsigned fgFlowHarmonicsHigher = 0;  calculate v6, v7, v8 ..
		// unsigned fgFlowHarmonicsMult = 0;    calculate v2{4} // yet v2{6}, v2{8}
		// unsigned fgNonlinearFlow = 0;        calculate v_4,22, v_5,32
		// unsigned fgSymmetricCumulants = 0;   calculate SC(3,2), SC(4,2)
		virtual void SetCalculateFlowHarmonics(unsigned flag)       { fgFlowHarmonics = flag; }
		virtual void SetCalculateFlowHarmonicsHigher(unsigned flag) { fgFlowHarmonicsHigher = flag; }
		virtual void SetCalculateFlowHarmonicsMult(unsigned flag)   { fgFlowHarmonicsMult = flag; }
		virtual void SetCalculateNonlinearFlow(unsigned flag)       { fgNonlinearFlow = flag; }
		virtual void SetCalculateSymmetricCumulants(unsigned flag)  { fgSymmetricCumulants = flag; }
    virtual void SetCalculateVnPtCorr(unsigned flag)  { fgVnPtCorr = flag; }


	private:
		AliAnalysisTaskNonlinearFlow(const AliAnalysisTaskNonlinearFlow&);
		AliAnalysisTaskNonlinearFlow& operator=(const AliAnalysisTaskNonlinearFlow&);

    AliMCEvent*     getMCEvent();
		virtual void		AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus);
		virtual void		AnalyzeMCTruth(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus);
    virtual void		AnalyzeMCOnTheFly(AliMCEvent* event);
		virtual void            NTracksCalculation(AliVEvent* aod);
		Bool_t                  AcceptAOD(AliAODEvent *inEv);
		Bool_t                  AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp);
		Bool_t                  AcceptMCTruthTrack(AliAODMCParticle *mtr);

    Bool_t                  LoadWeightsSystematics();
		Bool_t                  LoadWeightsKatarina();
		Bool_t                  LoadPtWeights();
		Bool_t                  LoadPtWeightsKatarina();

		Double_t GetWeightKatarina(double phi, double eta, double vz);
		Double_t GetPtWeightKatarina(double pt, double eta, double vz);
		Double_t GetFlowWeightSystematics(const AliVParticle* track, double fVtxZ, const PartSpecies species);
    double 			GetPtWeight(double pt, double eta, float vz, double runNumber);
    int GetEtaPtFlag(double eta);

		const char* ReturnPPperiod(const Int_t runNumber) const;
		const char* ReturnPPperiodMC(const Int_t runNumber) const;

		AliEventCuts	  fEventCuts;					// Event cuts
		AliGFWMCuts*     fGFWSelection;                                  //!
		AliGFWNFCuts*   fGFWSelection15o;                               //!
		AliAODEvent*    fAOD;                                           //! AOD object

		// Cuts and options
		Double_t		fEtaCut;				// Eta cut used to select particles
		Double_t		fVtxCut;				// Vtx cut on z position in cm
		Double_t		fMinPt;					// Min pt - for histogram limits
		Double_t		fMaxPt;					// Max pt - for histogram limits
		Int_t			  fTrigger;				// flag for trigger
		Int_t			  fAliTrigger;		// name for trigger
		Bool_t			fNUE;					  // flag for NUE correction
		Bool_t			fNUA;					  // 0: no NUA correction, 1: NUA correction
		bool        fIsMC;          // The observable for MonteCarlo truth
		TString     fNtrksName;     // Cent or Mult
		TString			fPeriod;				// period
		Int_t                   fCurrSystFlag;                          // Systematics flag
		Bool_t                  fSpringMode;                            // The mode with spring cuts.
		Bool_t                  fLowMultiplicityMode;                   // The mode to consider low-multiplicity region 
		Bool_t                  fAddTPCPileupCuts;                      // Additional TPC pileup cuts
    Double_t                fESDvsTPConlyLinearCut;                 // ESDvsTPConlyLinearCut : default = 15000
    Bool_t                  fUseOutOfBunchPileupCut;                // Out of bunch pileup cut
		Bool_t                  fUseCorrectedNTracks;                   // Use corrected Ntracks in the filling of xbins;
    Double_t                binning_factor;                         // A factor to account efficiency corrected binning
		Double_t                fCentralityCut;                         // Apply an extra centrality cut.
    Bool_t                  fUseNarrowBin;                          // Use Narrow bin
		Int_t                   fExtremeEfficiency;                     // The flag to set extreme efficiency
		Double_t                fTPCchi2perCluster;                     // Additional cuts for TPC chi2 / cluster
		Bool_t                  fUseAdditionalDCACut;                   // Additianal cuts for dca: < 1 cm
		Bool_t                  fUseDefaultWeight;                      // Force to use the default weight 
    Bool_t                  bUseLikeSign;                           // Flag to use like sign tracks
    Int_t                   iSign;                                  // Sign of selected tracks
    Bool_t                  fExtendV0MAcceptance;                   // Use V0M centrality cut 0-100%
    Double_t                fV0MRatioCut;                           // Cut on V0M / <V0M>
		Double_t                fEtaGap3Sub1;                           // The Eta Gap for 3 sub sample (Left most gap), the default is 0.4
    Double_t                fEtaGap3Sub2;                           // The Eta Gap for 3 sub sample (Middle gap), the default is 0.4
    Bool_t                  fOnTheFly;                              // flag to tune on on-the-fly

		// Output objects
		TList*			fListOfObjects;			//! Output list of objects
		TList*			fListOfProfile;			//! Output list of objects
		TList*			fListOfProfiles[30];		//! Output list of objects

		// NUE
		TH3F*			hTrackEfficiencyRun;            //! histogram with tracking efficiency (Katarina's format)

		// NUA
		TList*                  fFlowWeightsList;               //! flowWightsList
		TList*                  fFlowPtWeightsList;             //! PtflowWightsList
		TList*                  fFlowFeeddownList;              //! FeeddownList
		TList*			            fPhiWeight;	                    //! file with phi weights
		TFile*			            fPhiWeightFile;	                //! file with phi weights
		TH2D*                   fh2Weights[kUnknown];           //! container for GF weights (phi,eta,pt) (2D)
		TH3D*                   fh3Weights[kUnknown];           //! container for GF weights (phi,eta,pt)
		AliGFWWeights*          fWeightsSystematics;            //! Weights for systematics
		TH1D*                   fPtWeightsSystematics;          //! PtWeights for systematics
		TH1D*                   fPtWeightsFeeddown;             //! Feeddown for systematics
    TH1D*                   fEtaPtWeightsSystematics[8];          //! eta dependent PtWeights for systematics for pPb
    TH1D*                   fEtaPtWeightsFeeddown[8];             //! eta dependent Feeddown for systematics for pPb


		TH3F*			hPhiWeightRun;			//! 3D weight run-by-run for pPb 5TeV LHC16q
		TH1F*			hPhiWeight1D;			//! 1D weight in one MC case (maybe need to redo to 3D weight)

		// Event histograms
		TH1D*			hEventCount;			//! counting events passing given event cuts
		TH1F*			hMult;				//! multiplicity distribution
		TH1F*			fVtxAfterCuts;			//! Vertex z dist after cuts
		TH1F*			fCentralityDis;			//! distribution of centrality percentile using V0M estimator
    TH1F*			fV0CentralityDis;		//! distribution of centrality percentile using V0M estimator
    TH1F*			fV0CentralityDisNarrow;	//! distribution centrality percentile using V0M estimator
    TH1F*     fV0MMultiplicity;       //! V0M multiplicity
    TH1F*     fV0MRatio;              //! V0M multiplicity ratio: V0M/<V0M>

		// Track histograms
		TH1D*				fPhiDis1D;		//! phi dis 1D
		TH1D*				fPhiDis1DBefore;		//! phi dis 1D before track cuts
		TH3D*				fPhiDis;		//! phi dist
		TH1D*				fEtaDis;		//! eta dist
		TH1D*				fEtaBefore;		//! eta dist before track cuts
		TH1D*				fPtDis;			//! pt dist
		TH1D*				fPtBefore;		//! pt dist before track cuts
		TH2D*				hDCAxyBefore; 		//!
		TH1D*				hDCAzBefore; 		//!
		TH1F*				hITSclustersBefore; 	//!
		TH1D*				hChi2Before; 		//!
    TH1D*				hnTPCClu;  		//!
		TH2D*				hDCAxy; 		//!
		TH1D*				hDCAz; 			//!
		TH1F*				hITSclusters; 		//!
		TH1D*				hChi2; 			//!

		TH2D*                           hTracksCorrection2d;    //! Corrected Tracks - v.s. uncorrected tracks
		TProfile*                       hnCorrectedTracks;      //! Averaged number of corrected tracks in a specific bin;

    // TH2D* QDis[10];        // QDistribution for No gap
    // TH2D* QDisGap0P[10];        // QDistribution for gap 0
    // TH2D* QDisGap0M[10];        // QDistribution for gap 0
    // TH2D* QDisGap10P[10];        // QDistribution for gap 10
    // TH2D* QDisGap10M[10];        // QDistribution for gap 10
    // TH2D* QDisGap14P[10];        // QDistribution for gap 14
    // TH2D* QDisGap14M[10];        // QDistribution for gap 14
    // TH2D* QDis3subL[10];        // QDistribution for 3sub
    // TH2D* QDis3subM[10];        // QDistribution for 3sub
    // TH2D* QDis3subR[10];        // QDistribution for 3sub

    AliMCEvent *fMCEvent;           //! MC event

		// Global variables
		double NtrksCounter = 0;       //!
		double NTracksCorrected = 0;   //!
		double NTracksUncorrected = 0; //!
		int NtrksAfter = 0;            //!
		int NtrksAfterGap0M = 0;       //!
		int NtrksAfterGap0P = 0;       //!
		int NtrksAfterGap2M = 0;       //!
		int NtrksAfterGap2P = 0;       //!
		int NtrksAfterGap4M = 0;       //!
		int NtrksAfterGap4P = 0;       //!
		int NtrksAfterGap6M = 0;       //!
		int NtrksAfterGap6P = 0;       //!
		int NtrksAfterGap8M = 0;       //!
		int NtrksAfterGap8P = 0;       //!
		int NtrksAfterGap10M = 0;      //!
		int NtrksAfterGap10P = 0;      //!
		int NtrksAfterGap14M = 0;      //!
		int NtrksAfterGap14P = 0;      //!
		int NtrksAfter3subL = 0;       //!
		int NtrksAfter3subM = 0;       //!
		int NtrksAfter3subR = 0;       //!

		int lastRunNumber = 0;         //!

		PhysicsProfile multProfile;    //!
		PhysicsProfile multProfile_bin[30]; //!

    CorrelationCalculator correlator; //!
    Double_t sumPtw; //!
    Double_t sumPtw2; //!
    Double_t sumPt2w2; //!
    Double_t sumWeight; //!
    Double_t sumWeight2; //!
    Double_t eventWeight; //!
    Double_t eventWeight2; //!
    TRandom3 rand;         //!
    Int_t bootstrap_value = -1; //!


		unsigned fgFlowHarmonics = 0;        // calculate v2, v3, v4, v5
		unsigned fgFlowHarmonicsHigher = 0;  // calculate v6, v7, v8 ..
		unsigned fgFlowHarmonicsMult = 0;    // calculate v2{4} // yet v2{6}, v2{8}
    unsigned fgNonlinearFlow = 0;        // calculate v_4,22, v_5,32
		unsigned fgSymmetricCumulants = 0;   // calculate SC(3,2), SC(4,2)
    unsigned fgVnPtCorr = 0;             // calculate <v2^2-[pt]>

		unsigned fgTwoParticleCorrelation = 0;       //!
		unsigned fgTwoParticleCorrelationHigher = 0; //!
		unsigned fgThreeParticleCorrelation = 0;     //!
		unsigned fgFourParticleCorrelation = 0;      //!
		unsigned fgSixParticleCorrelation = 0;       //!
		unsigned fgEightParticleCorrelation = 0;     //!

		bool fuTwoParticleCorrelationStandard = 0;        //!
		bool fuTwoParticleCorrelation0Gap = 0;            //!
		bool fuTwoParticleCorrelationLargeGap = 0;        //!
		bool fuTwoParticleCorrelationThreeSub = 0;        //!
		bool fuTwoParticleCorrelationGapScan = 0;         //!
		bool fuTwoParticleCorrelationHigherStandard = 0;  //!
		bool fuTwoParticleCorrelationHigher0Gap = 0;      //!
		bool fuTwoParticleCorrelationHigherLargeGap = 0;  //!
		bool fuTwoParticleCorrelationHigherThreeSub = 0;  //!
		bool fuTwoParticleCorrelationHigherGapScan = 0;   //!
		bool fuThreeParticleCorrelationStandard = 0;      //!
		bool fuThreeParticleCorrelation0Gap = 0;          //!
		bool fuThreeParticleCorrelationLargeGap = 0;      //!
		bool fuThreeParticleCorrelationThreeSub = 0;      //!
		bool fuThreeParticleCorrelationGapScan = 0;       //!
		bool fuFourParticleCorrelationStandard = 0;       //!
		bool fuFourParticleCorrelation0Gap = 0;           //!
		bool fuFourParticleCorrelationLargeGap = 0;       //!
		bool fuFourParticleCorrelationThreeSub = 0;       //!
		bool fuFourParticleCorrelationGapScan = 0;        //!
		bool fuSixParticleCorrelationStandard = 0;        //!
		bool fuSixParticleCorrelation0Gap = 0;            //!
		bool fuEightParticleCorrelationStandard = 0;      //!
		bool fuEightParticleCorrelation0Gap = 0;          //!

		bool fuQStandard = 0; //!
		bool fuQ0Gap     = 0; //!
		bool fuQLargeGap = 0; //!
		bool fuQThreeSub = 0; //!
		bool fuQGapScan  = 0; //!

    Double_t fImpactParameterMC;                       //! Impact parameter

		double xbins[3000+10] = {}; //!
		int nn = 0; //!
		void CalculateProfile(PhysicsProfile& profile, double Ntrks);
		void InitProfile(PhysicsProfile& profile, TString name, TList* listOfProfile);

		ClassDef(AliAnalysisTaskNonlinearFlow, 27);    //Analysis task
};

#endif
// Local Variables:
//  mode: C++
// End:
