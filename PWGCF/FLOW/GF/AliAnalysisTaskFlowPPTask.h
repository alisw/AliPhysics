/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskFlowPPTask_H
#define AliAnalysisTaskFlowPPTask_H

#include "AliAnalysisTaskSE.h"
#include "AliGFWXXCuts.h"
#include "AliGFWPbpass23Cuts.h"
#include "AliGFWWeights.h"
#include "CorrelationCalculator.h"
#include "AliEventCuts.h"
#include <TComplex.h>
#include "AliMultSelection.h"

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
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"
#include <THnSparse.h>
// Preload to boost program
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


class PhysicsProfilePPTask : public TObject {
	public:
        PhysicsProfilePPTask();
		// Physics profiles
		TProfile*	 fChsc4242;		             	//! SC(4,2)
		TProfile*	 fChsc4242_Gap0;			//! SC(4,2) |#Delta#eta| > 0.0
		TProfile*	 fChsc4242_Gap2;			//! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap4;                        //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap6;                        //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap8;                        //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap10;                       //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap12;                       //! SC(4,2) |#Delta#eta| > 0.2
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
		TProfile*        fChsc3232_Gap12;                       //! SC(3,2) |#Delta#eta| > 0.2
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
		TProfile*        fChc422_Gap12A;    //!
		TProfile*        fChc422_Gap12B;    //!
		TProfile*        fChc532_Gap12A;    //!
		TProfile*        fChc532_Gap12B;    //!

		TProfile*	 fChcn2[5]; 			//! <<2>> in unit bins of Ntrks
		TProfile*    fChcn2_Gap0[5];  		//! <<2>> |#Delta#eta| > 0.0
		TProfile*	 fChcn2_Gap2[5];  		//! <<2>> |#Delta#eta| > 0.2
		TProfile*	 fChcn2_Gap4[5];  		//! <<2>> |#Delta#eta| > 0.4
		TProfile*    fChcn2_Gap6[5];        //! <<2>> |#Delta#eta| > 0.6
		TProfile*	 fChcn2_Gap8[5];  		//! <<2>> |#Delta#eta| > 0.8
		TProfile*	 fChcn2_Gap10[5];  		//! <<2>> |#Delta#eta| > 1.0
		TProfile*	 fChcn2_Gap12[5];  		//! <<2>> |#Delta#eta| > 1.2
		TProfile*	 fChcn2_Gap14[5];  	        //! <<2>> |#Delta#eta| > 1.4
		TProfile*	 fChcn2_Gap16[5];  	        //! <<2>> |#Delta#eta| > 1.6
		TProfile*	 fChcn2_Gap18[5];  	        //! <<2>> |#Delta#eta| > 1.8

		TProfile*	 fChcn2_3subLM[5];  		//! <<2>> left vs. middle subevent
		TProfile*	 fChcn2_3subRM[5];  		//! <<2>> middle vs. right subevent
		TProfile*	 fChcn2_3subLR[5];  		//! <<2>> middle vs. right subevent
		TProfile*	 fChcn2_3subGap2LM[5];          //! <<2>> left vs. middle subevent |#Delta#eta| > 0.2
		TProfile*	 fChcn2_3subGap2RM[5];          //! <<2>> middle vs. right subevent |#Delta#eta| > 0.2

		TProfile*	 fChcn4[5];  			//! <<4>> in unit bins of Ntrks
		TProfile*	 fChcn4_Gap0[5];  		//! <<4>> |#Delta#eta| > 0.0
		TProfile*	 fChcn4_Gap2[5];  		//! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap4[5];                //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap6[5];                //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap8[5];                //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap10[5];               //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_Gap12[5];               //! <<4>> |#Delta#eta| > 0.2
		TProfile*        fChcn4_3subMMLR[5];            //! <<4>> 3subevent method
		TProfile*        fChcn4_3subLLMR[5];            //! <<4>> 3subevent method
		TProfile*        fChcn4_3subRRML[5];            //! <<4>> 3subevent method
		TProfile*        fChcn4_3subGap2[5];            //! <<4>> 3subevent method |#Delta#eta| > 0.2

	
		//5,6 particle correlation
		//TProfile* fChc5_A42222;  		//! <<5>> 
		//TProfile* fChc5_A52322;  		//! <<5>> 
		//TProfile* fChc6_222222;  		//! <<6>> 
		//TProfile* fChc6_322322;  		//! <<6>> 

		//Addtional 6 particle correlation
		//TProfile* fChsc6222_Gap0;  		//! <<4>> |#Delta#eta| > 1.0
		//TProfile* fChsc6222_Gap10;  		//! <<4>> |#Delta#eta| > 1.0
		//TProfile* fChsc633_Gap0A;  		//! <<3>> |#Delta#eta| > 1.0
		//TProfile* fChsc633_Gap10A;  		//! <<3>> |#Delta#eta| > 1.0
		TProfile*	fChcn6[6];  			//! <<6>> in unit bins of Ntrks
		TProfile*   fChcn6_Gap10[6];               //! <<6>> |#Delta#eta| > 1.0
		TProfile*   fChcn6_Gap0[6];               //! <<6>> |#Delta#eta| > 0.

		// 8 particle correlation
		TProfile*	fChcn8[6];  			//! <<8>> in unit bins of Ntrks
		TProfile*   fChcn8_Gap0[6];               //! <<8>> |#Delta#eta| > 1.0

		private:
		ClassDef(PhysicsProfilePPTask, 1);    //Analysis task
};


class AliAnalysisTaskFlowPPTask : public AliAnalysisTaskSE  
{
    public:
        enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kCharUnidentified, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest; NB: kUknown last as counter
                                AliAnalysisTaskFlowPPTask();
                                AliAnalysisTaskFlowPPTask(const char *name);
        virtual                 ~AliAnalysisTaskFlowPPTask();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
		virtual void   			NotifyRun();
        virtual void            Terminate(Option_t* option);

        
        virtual void   SetFilterbit(Int_t bit){fFilterbit = bit;}
		virtual void   SetFilterbitDefault(Int_t bit){fFilterbitDefault = bit;}
		virtual void   SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void   SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void   SetVtxCutDefault(Double_t vtxCut){fVtxCutDefault = vtxCut;} // The vtxCut for NtrksCounter
		virtual void   SetMinPt(Double_t minPt){fMinPt = minPt;}
		virtual void   SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
		virtual void   SetTPCclusters(Int_t tpcClus){fTPCclusters = tpcClus;}
		virtual void   SetTPCclustersDefault(Int_t tpcClus){fTPCclustersDefault = tpcClus;} // The tpcCluster for NtrksCounter 
        virtual void   SetChi2PerTPCcluster(Double_t chi2){fChi2PerTPCcluster = chi2;}	    // max. chi2 per TPC cluster
		virtual void	SetEventWeightSetOne(Bool_t useOne){
			fEventWeightSetToOne=useOne;	//Set Event Weight to 1
		}
		virtual void	SetAddTPCPileupCuts(Bool_t usePileupCuts){fAddTPCPileupCuts=usePileupCuts;} // use TPC pile up Cuts
		virtual void   SetESDvsTPConlyLinearCut(double cut = 15000) {fESDvsTPConlyLinearCut = cut;}
		virtual void   SetMinITSClusters(Int_t minClus){fMinITSClus = minClus;}
		virtual void   SetMaxChi2(Double_t maxChi){fMaxChi2 = maxChi;}
		virtual void   SetUseDCAzCut(Bool_t usedcaz){fUseDCAzCut = usedcaz;}
		virtual void   SetDCAzCut(Double_t dcaz){fDCAz = dcaz;}
		virtual void   SetDCAzCutDefault(Double_t dcaz){fDCAzDefault = dcaz;} // The dcaz cut for NtrksCounter
		virtual void   SetUseDCAxyCut(Bool_t usedcaxy){fUseDCAxyCut = usedcaxy;}
		virtual void   SetDCAxyCut(Double_t dcaxy){fDCAxy = dcaxy;}
		virtual void   SetDCAxyCutDefault(Double_t dcaxy){fDCAxyDefault = dcaxy;} // The dcaxy cut for NtrksCounter
		virtual void   SetUseCL1Cent(Bool_t CL1cent){fUseCL1Centrality = CL1cent;}
		virtual void   SetIsSample(Int_t IsSample){fSample = IsSample;}
		virtual void   SetCentFlag(Short_t nCent){fCentFlag = nCent;}
		virtual void   SetTrigger(Int_t trig){fTrigger = trig;}
		virtual void   SetLSFlag(Bool_t LS){fLS = LS;}
		virtual void   SetNUEFlag(Bool_t NUE){fNUE = NUE;}
		virtual void   SetNUA(Bool_t NUA){fNUA = NUA;}
		virtual void   SetNtrksName(TString ntrksname){fNtrksName = ntrksname;}
		virtual void   SetUseWeigthsRunByRun(Bool_t bRunByRun = kTRUE) { fFlowRunByRunWeights = bRunByRun; }
		virtual void   SetUsePeriodWeigths(Bool_t weight = kTRUE) { fFlowPeriodWeights = weight; }
		virtual void   SetUseWeights3D(Bool_t use = kTRUE) { fFlowUse3Dweights = use; }
		virtual void   SetPeriod(TString period) { fPeriod = period; }
        virtual void   SetSystFlag(int flag) { fCurrSystFlag = flag; }
        virtual int    GetSystFlag() { return fCurrSystFlag; }
		virtual void	SetOnlineTrackCorrection(Bool_t UseCorrectedNTracks){fUseCorrectedNTracks = UseCorrectedNTracks;}
		virtual void	SetAdditionalDCACut(Bool_t UseAdditionalDCACut){fUseAdditionalDCACut = UseAdditionalDCACut;}
        
        //===================================================================================

    private:

        AliAnalysisTaskFlowPPTask(const AliAnalysisTaskFlowPPTask&); // not implemented
        AliAnalysisTaskFlowPPTask& operator=(const AliAnalysisTaskFlowPPTask&); // not implemented

        
        virtual void		AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus);
		virtual void            NTracksCalculation(AliVEvent* aod);
		Bool_t			CheckTrigger();
		double 			GetPtWeight(double pt, double eta, float vz, double runNumber);
		Bool_t                  LoadWeightsSystematics();

		Double_t GetFlowWeightSystematics(const AliVParticle* track, double fVtxZ, const PartSpecies species);
        

		AliEventCuts	fEventCuts;					// Event cuts
		AliGFWXXCuts*     fGFWSelection;                                  //!
		AliGFWPbpass23Cuts*     fGFWSelectionPbPb;                                  //!
		AliAODEvent* fAOD;                                              //! AOD object

		// Cuts and options
		Int_t			fFilterbit;				// flag for filter bit
		Int_t			fFilterbitDefault;			// flag for filter bit (for NtrksCounter)
		Double_t		fEtaCut;				// Eta cut used to select particles
		Double_t		fVtxCut;				// Vtx cut on z position in cm
		Double_t		fVtxCutDefault;				// Vtx cut on z position in cm (for NtrksCounter)
		Double_t		fMinPt;					// Min pt - for histogram limits
		Double_t		fMaxPt;					// Max pt - for histogram limits
		Int_t			fTPCclusters;				// min. TPC clusters
		Int_t			fTPCclustersDefault;			// min. TPC clusters (for NtrksCounter)
		Int_t			fChi2PerTPCcluster;			// max. chi2 per TPC cluster
		Bool_t			fEventWeightSetToOne;		// Set Event weight to 1
		Bool_t			fAddTPCPileupCuts;			// Pile up Cuts in TPC
		Double_t		fESDvsTPConlyLinearCut;     // ESDvsTPConlyLinearCut : default = 15000
		Int_t			fMinITSClus;				// min ITS clusters, LHC15ijl
		Double_t		fMaxChi2;				// max chi2 per ITS cluster, LHC15ijl
		Bool_t			fUseDCAzCut;				// switch to choose whether I want to use DCAz cut or not (for systematic studies, otherwise it is in FB selection by default)
		Bool_t			fUseCL1Centrality;				// switch to choose whether I want to use DCAz cut or not (for systematic studies, otherwise it is in FB selection by default)
		Double_t		fDCAz;					// max DCAz, for systematics
		Double_t		fDCAzDefault;				// max DCAz, (for NtrksCounter)
		Bool_t			fUseDCAxyCut;				// the same switch as for DCAxy
		Double_t		fDCAxy;					// max DCAxy, for systematics
		Double_t		fDCAxyDefault;				// max DCAxy, (for NtrksCounter)
		Int_t			fSample;				// number of sample
		Short_t			fCentFlag;				// centrality flag
		Int_t			fTrigger;				// flag for trigger
		Int_t			fAliTrigger;				// name for trigger
		Bool_t			fLS;					// charge, 1:all, 2:pp,  3: mm
		Bool_t			fNUE;					// flag for NUE correction
		Bool_t			fNUA;					// 0: no NUA correction, 1: NUA correction
		TString                 fNtrksName;                             // Cent or Mult
		TString			fPeriod;				// period
                Int_t                   fCurrSystFlag;                          // Systematics flag
		Bool_t			fUseCorrectedNTracks;	// flag for online track correction
		Bool_t			fUseAdditionalDCACut;	// flag for online track correction


		Double_t			fCurrCentrality;	//! Centrality of Ongoing Event

		// Output objects
		TList*			fListOfObjects;			//! Output list of objects

		// NUE
		TList*			fTrackEfficiency;		//! file with tracking efficiency
		//TH1D*			hTrackEfficiency;		//! histogram with tracking efficiency
		TH1D*			hTrackEfficiencyRun;            //! histogram with tracking efficiency

		// NUA
		bool fFlowRunByRunWeights;                              // flag of whether get the Run by run weight 
		bool fFlowPeriodWeights;                                // flag of whether to use period weight
		bool fFlowUse3Dweights;                                 // flag of whether to use 3d weight

		//
		TList*                  fFlowWeightsList;               //! flowWightsList
		AliGFWWeights*          fWeightsSystematics;            //! Weights for systematics


		TH3F*			hPhiWeight;			//! 3D weight for all periods except LHC15ijl

		// Event histograms
		TH1D*			hEventCount;			//! counting events passing given event cuts
		TH2D*			hTracksCorrection2d;	//! Correlation table for number of tracks table
		TProfile*		hnCorrectedTracks;		//! Number of corrected tracks in a ntracks bin
		TH1F*			hMult;				//! multiplicity distribution
		TH2F*			hMultCent;				//! multiplicity-centrality distribution
		TH1F*			fCentralityDis;			//! distribution of centrality percentile using V0M estimator
		TH1F*			hVtz;			//! distribution of Vertex


		// // Track histograms
		TH2F*				hDCAxyBefore; 		//!
		TH2F*				hDCAzBefore; 		//!
		TH2F*				hDCAxy; 		//!
		TH2F*				hDCAz; 			//!
		TH1D*				hChi2Before; 			//!
		TH1D*				hChi2; 			//!
		TH1D*			hPhi;			//! distribution of Phi
		TH1D*			hPhiBefore;			//! distribution of Phi
		TH1D*			fEtaDis;			//! distribution of Phi
		TH1D*			fPtDis;			//! distribution of Phi
		TH1D*			hnTPCClu;			//! distribution of Phi

		// Global variables
		double NtrksCounter = 0;               //!
		double NTracksCorrected = 0;			//!
		double NTracksUncorrected = 0;			//!
		int NtrksAfter = 0;                 //!
		int NtrksAfterGap0M = 0;           //!
		int NtrksAfterGap0P = 0;           //!
		int NtrksAfterGap2M = 0;           //!
		int NtrksAfterGap2P = 0;           //!
		int NtrksAfterGap4M = 0;           //!
		int NtrksAfterGap4P = 0;           //!
		int NtrksAfterGap6M = 0;           //!
		int NtrksAfterGap6P = 0;           //!
		int NtrksAfterGap8M = 0;           //!
		int NtrksAfterGap8P = 0;           //!
		int NtrksAfterGap10M = 0;           //!
		int NtrksAfterGap10P = 0;           //!
		int NtrksAfterGap12M = 0;           //!
		int NtrksAfterGap12P = 0;           //!
		int NtrksAfterGap14M = 0;           //!
		int NtrksAfterGap14P = 0;           //!
		int NtrksAfter3subL = 0;            //!
		int NtrksAfter3subM = 0;            //!
		int NtrksAfter3subR = 0;            //!
		
        
		PhysicsProfilePPTask multProfile;          //!
		PhysicsProfilePPTask multProfile_bin[30];  //!
        
		TRandom3 rand;                      //!
		Int_t bootstrap_value;              //!

		CorrelationCalculator correlator;   //!

		double xbins[3000+10];              //!
		int nn;                             //!
		//TH1F* MyEventNumber;					//!
        
		Bool_t AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp);
		void CalculateProfile(PhysicsProfilePPTask& profile, double Ntrks);
		void InitProfile(PhysicsProfilePPTask& profile, TString);
		Bool_t AcceptAODEvent(AliAODEvent* aliev);
        
        //===================================================================================

        ClassDef(AliAnalysisTaskFlowPPTask, 1);
};

#endif
