#ifndef ALIANALYSISTASKNONLINEARFLOW_H
#define ALIANALYSISTASKNONLINEARFLOW_H
#include "AliAnalysisTaskSE.h"
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
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

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
		// Physics profiles
		TProfile*		fChsc4242;									//! SC(4,2)
		TProfile*		fChsc4242_Gap0;							//! SC(4,2) |#Delta#eta| > 0.0
		TProfile*		fChsc4242_Gap2;							//! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap4;                            //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap6;                            //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap8;                            //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*        fChsc4242_Gap10;                            //! SC(4,2) |#Delta#eta| > 0.2
		TProfile*		fChsc4242_3sub;							//! SC(4,2)_A 3subevent method
		TProfile*		fChsc4242_3subMMLRA;							//! SC(4,2)_A 3subevent method
		TProfile*		fChsc4242_3subMMLRB;							//! SC(4,2)_A 3subevent method
		TProfile*		fChsc4242_3subLLMRA;							//! SC(4,2)_A 3subevent method
		TProfile*		fChsc4242_3subLLMRB;							//! SC(4,2)_A 3subevent method
		TProfile*		fChsc4242_3subRRMLA;							//! SC(4,2)_A 3subevent method
		TProfile*		fChsc4242_3subRRMLB;							//! SC(4,2)_A 3subevent method
		TProfile*		fChsc4224_3sub;							//! SC(4,2)_B 3subevent method
		TProfile*		fChsc4242_3subGap2;					//! SC(4,2)_A 3subevent method |#Delta#eta| > 0.2
		TProfile*		fChsc4224_3subGap2;					//! SC(4,2)_B 3subevent method |#Delta#eta| > 0.2
		TProfile*		fChsc3232;									//! SC(3,2)
		TProfile*		fChsc3232_Gap0;							//! SC(3,2) |#Delta#eta| > 0.0
		TProfile*		fChsc3232_Gap2;							//! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap4;                            //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap6;                            //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap8;                            //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*        fChsc3232_Gap10;                            //! SC(3,2) |#Delta#eta| > 0.2
		TProfile*		fChsc3232_3sub;							//! SC(3,2)_A 3subevent method
		TProfile*		fChsc3232_3subMMLRA;							//! SC(3,2)_A 3subevent method
		TProfile*		fChsc3232_3subMMLRB;							//! SC(3,2)_A 3subevent method
		TProfile*		fChsc3232_3subLLMRA;							//! SC(3,2)_A 3subevent method
		TProfile*		fChsc3232_3subLLMRB;							//! SC(3,2)_A 3subevent method
		TProfile*		fChsc3232_3subRRMLA;							//! SC(3,2)_A 3subevent method
		TProfile*		fChsc3232_3subRRMLB;							//! SC(3,2)_A 3subevent method
		TProfile*		fChsc3223_3sub;							//! SC(3,2)_B 3subevent method
		TProfile*		fChsc3232_3subGap2;					//! SC(3,2)_A 3subevent method |#Delta#eta| > 0.2
		TProfile*		fChsc3223_3subGap2;					//! SC(3,2)_B 3subevent method |#Delta#eta| > 0.2

		TProfile*		fsc4242[12];								//! the same but for different fBin (sampling)
		TProfile*		fsc4242Gap0[12];						//!
		TProfile*		fsc4242Gap2[12];						//!
		TProfile*        fsc4242Gap4[12];                        //!
		TProfile*        fsc4242Gap6[12];                        //!
		TProfile*        fsc4242Gap8[12];                        //!
		TProfile*        fsc4242Gap10[12];                        //!
		TProfile*		fsc4242_3sub[12];						//!
		TProfile*		fsc4224_3sub[12];						//!
		TProfile*		fsc4242_3subGap2[12];				//!
		TProfile*		fsc4224_3subGap2[12];				//!
		TProfile*		fsc3232[12];								//!
		TProfile*		fsc3232Gap0[12];						//!
		TProfile*		fsc3232Gap2[12];						//!
		TProfile*        fsc3232Gap4[12];                        //!
		TProfile*        fsc3232Gap6[12];                        //!
		TProfile*        fsc3232Gap8[12];                        //!
		TProfile*        fsc3232Gap10[12];                        //!
		TProfile*		fsc3232_3sub[12];						//!
		TProfile*		fsc3223_3sub[12];						//!
		TProfile*		fsc3232_3subGap2[12];				//!
		TProfile*		fsc3223_3subGap2[12];				//!

		// Standard correlation profiles for different harmonics
		TProfile*	     fChc422; //!
		TProfile*        fChc532;//!
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


		TProfile*		fChcn2[6];  				//! <<2>> in unit bins of Ntrks
		TProfile*		fChcn2_Gap0[6];  		//! <<2>> |#Delta#eta| > 0.0
		TProfile*		fChcn2_Gap2[6];  		//! <<2>> |#Delta#eta| > 0.2
		TProfile*		fChcn2_Gap4[6];  		//! <<2>> |#Delta#eta| > 0.4
		TProfile*           fChcn2_Gap6[6];          //! <<2>> |#Delta#eta| > 0.6
		TProfile*		fChcn2_Gap8[6];  		//! <<2>> |#Delta#eta| > 0.8
		TProfile*		fChcn2_Gap10[6];  		//! <<2>> |#Delta#eta| > 1.0
		TProfile*		fChcn2_Gap14[6];  	//! <<2>> |#Delta#eta| > 1.4
		TProfile*		fChcn2_Gap16[6];  	//! <<2>> |#Delta#eta| > 1.6
		TProfile*		fChcn2_Gap18[6];  	//! <<2>> |#Delta#eta| > 1.8

		TProfile*		fChcn2_3subLM[6];  		//! <<2>> left vs. middle subevent
		TProfile*		fChcn2_3subRM[6];  		//! <<2>> middle vs. right subevent
		TProfile*		fChcn2_3subLR[6];  		//! <<2>> middle vs. right subevent
		TProfile*		fChcn2_3subGap2LM[6];  //! <<2>> left vs. middle subevent |#Delta#eta| > 0.2
		TProfile*		fChcn2_3subGap2RM[6];  //! <<2>> middle vs. right subevent |#Delta#eta| > 0.2

		TProfile*		fChcn4[6];  				//! <<4>> in unit bins of Ntrks
		TProfile*		fChcn4_Gap0[6];  		//! <<4>> |#Delta#eta| > 0.0
		TProfile*		fChcn4_Gap2[6];  		//! <<4>> |#Delta#eta| > 0.2
		TProfile*      fChcn4_Gap4[6];          //! <<4>> |#Delta#eta| > 0.2
		TProfile*      fChcn4_Gap6[6];          //! <<4>> |#Delta#eta| > 0.2
		TProfile*      fChcn4_Gap8[6];          //! <<4>> |#Delta#eta| > 0.2
		TProfile*      fChcn4_Gap10[6];          //! <<4>> |#Delta#eta| > 0.2
		TProfile*     fChcn4_3subMMLR[6];         //! <<4>> 3subevent method
		TProfile*     fChcn4_3subLLMR[6];         //! <<4>> 3subevent method
		TProfile*     fChcn4_3subRRML[6];         //! <<4>> 3subevent method
		TProfile*     fChcn4_3subGap2[6];//! <<4>> 3subevent method |#Delta#eta| > 0.2
		private:
		ClassDef(PhysicsProfile, 1);    //Analysis task
};

class AliAnalysisTaskNonlinearFlow : public AliAnalysisTaskSE {
	public:

      enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kCharUnidentified, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest; NB: kUknown last as counter

		AliAnalysisTaskNonlinearFlow();
		AliAnalysisTaskNonlinearFlow(const char *name);

		virtual ~AliAnalysisTaskNonlinearFlow();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t* option);
		virtual void   Terminate(Option_t* );

		virtual void   SetFilterbit(Int_t bit){fFilterbit = bit;}
		virtual void   SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void   SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
		virtual void   SetMinPt(Double_t minPt){fMinPt = minPt;}
		virtual void   SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
		virtual void   SetTPCclusters(Int_t tpcClus){fTPCclusters = tpcClus;}
		virtual void   SetMinITSClusters(Int_t minClus){fMinITSClus = minClus;}
		virtual void   SetMaxChi2(Double_t maxChi){fMaxChi2 = maxChi;}
		virtual void   SetUseDCAzCut(Bool_t usedcaz){fUseDCAzCut = usedcaz;}
		virtual void   SetDCAzCut(Double_t dcaz){fDCAz = dcaz;}
		virtual void   SetUseDCAxyCut(Bool_t usedcaxy){fUseDCAxyCut = usedcaxy;}
		virtual void   SetDCAxyCut(Double_t dcaxy){fDCAxy = dcaxy;}
		virtual void   SetIsSample(Int_t IsSample){fSample = IsSample;}
		virtual void   SetCentFlag(Short_t nCent){fCentFlag = nCent;}
		virtual void   SetTrigger(Int_t trig){fTrigger = trig;}
		virtual void   SetLSFlag(Bool_t LS){fLS = LS;}
		virtual void   SetNUEFlag(Bool_t NUE){fNUE = NUE;}
		virtual void   SetNUA(Bool_t NUA){fNUA = NUA;}
		virtual void   SetNtrksName(TString ntrksname){fNtrksName = ntrksname;}
		void           SetUseWeigthsRunByRun(Bool_t bRunByRun = kTRUE) { fFlowRunByRunWeights = bRunByRun; }
		void           SetUsePeriodWeigths(Bool_t weight = kTRUE) { fFlowPeriodWeights = weight; }
		void           SetUseWeights3D(Bool_t use = kTRUE) { fFlowUse3Dweights = use; }
		virtual void   SetPeriod(TString period){fPeriod = period;}

	private:
		AliAnalysisTaskNonlinearFlow(const AliAnalysisTaskNonlinearFlow&);
		AliAnalysisTaskNonlinearFlow& operator=(const AliAnalysisTaskNonlinearFlow&);

		virtual void		AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus);
		Short_t			GetCentrCode(AliVEvent* ev);
		bool 			CheckPrimary(AliVEvent *aod, double label);
		bool			IsGoodPSEvent(AliVEvent *aod);
		bool			IsSPDClusterVsTrackletBG(const AliVEvent* event, bool fillHist);
		bool			IsV0C012vsTklBG         (const AliVEvent* event, bool fillHist);
		bool			IsV0Casym               (const AliVEvent* event, bool fillHist);
		bool			IsV0MOnVsOfPileup       (const AliVEvent* event, bool fillHist);
		bool			IsSPDOnVsOfPileup       (const AliVEvent* event, bool fillHist);
		bool			IsV0PFPileup            (const AliVEvent* event);
		int 			GetRunPart(int run);
		double 			GetWeight(double phi, double eta, double pt, int run, bool fPlus, double vz, double runNumber);
		double 			GetPtWeight(double pt, double eta, float vz, double runNumber);
		Bool_t                  LoadWeights();

		Double_t GetFlowWeight(const AliVParticle* track, double fVtxZ, const PartSpecies species);
                const char* ReturnPPperiod(const Int_t runNumber) const;
                const char* GetSpeciesName(const PartSpecies species) const;

		AliEventCuts	fEventCuts;					// Event cuts
		AliAODEvent* fAOD;                //! AOD object
		AliAODITSsaTrackCuts* fitssatrackcuts; //itssatrackcuts object

		// Cuts and options
		Int_t				fFilterbit;					// flag for filter bit
		Double_t		fEtaCut;						// Eta cut used to select particles
		Double_t		fVtxCut;						// Vtx cut on z position in cm
		Double_t		fMinPt;							// Min pt - for histogram limits
		Double_t		fMaxPt;							// Max pt - for histogram limits
		Int_t				fTPCclusters;				// min. TPC clusters
		Int_t				fMinITSClus;				// min ITS clusters, LHC15ijl
		Double_t		fMaxChi2;						// max chi2 per ITS cluster, LHC15ijl
		Bool_t			fUseDCAzCut;				// switch to choose whether I want to use DCAz cut or not (for systematic studies, otherwise it is in FB selection by default)
		Double_t		fDCAz;							// max DCAz, for systematics
		Bool_t			fUseDCAxyCut;				// the same switch as for DCAxy
		Double_t		fDCAxy;							// max DCAxy, for systematics
		Int_t				fSample;						// number of sample
		Short_t			fCentFlag;					// centrality flag
		Int_t				fTrigger;						// flag for trigger
		Bool_t			fLS;								// charge, 1:all, 2:pp,  3: mm
		Bool_t			fNUE;								// flag for NUE correction
		Bool_t			fNUA;								// 0: no NUA correction, 1: NUA correction
		TString                     fNtrksName;
		//....
		TString			fPeriod;						// period

		// Output objects
		TList*			fListOfObjects;			//! Output list of objects

		// Cut functions for LHC15o
		TF1*				fMultTOFLowCut;			// cut low for TOF multiplicity outliers
		TF1*				fMultTOFHighCut;		// cut high for TOF multiplicity outliers
		TF1*				fMultCentLowCut;		// cut low for multiplicity centrality outliers

		// NUE
		TFile*			fTrackEfficiency;		//! file with tracking efficiency
		TH3F*				hTrackEfficiency;		//! histogram with tracking efficiency
		TH3F*				hTrackEfficiencyRun;//! histogram with tracking efficiency

		// NUA
		bool fFlowRunByRunWeights;            // flag of whether get the Run by run weight 
		bool fFlowPeriodWeights;              // flag of whether to use period weight
		bool fFlowUse3Dweights;               // flag of whether to use 3d weight
		//
		TList*                  fFlowWeightsList;     //! flowWightsList
		TList*			fPhiWeight;	      //! file with phi weights
		TList*			fPhiWeightPlus;	      //! file with phi weights
		TList*			fPhiWeightMinus;      //! file with phi weights
                TH2D*                   fh2Weights[kUnknown]; //! container for GF weights (phi,eta,pt) (2D)
                TH3D*                   fh3Weights[kUnknown]; //! container for GF weights (phi,eta,pt)
                TH2D*                   fh2AfterWeights[kUnknown]; //! distribution after applying GF weights - lightweight QA (phi)
                TH3D*                   fh3AfterWeights[kUnknown]; //! distribution after applying GF weights - full QA (phi,eta,pt)


		TH3F*				hPhiWeight;					//! 3D weight for all periods except LHC15ijl
		TH3F*				hPhiWeightRun;			//! 3D weight run-by-run for pPb 5TeV LHC16q
		TH1F*				hPhiWeight1D;				//! 1D weight in one MC case (maybe need to redo to 3D weight)
		TH3D*				hPhiWeight_LHC15i_part1;	//! LHC15i, part1 runs
		TH3D*				hPhiWeight_LHC15i_part2;	//! LHC15i, part2 runs
		TH3D*				hPhiWeight_LHC15j_part1;	//! LHC15j, part1 runs
		TH3D*				hPhiWeight_LHC15j_part2;	//! LHC15j, part2 runs
		TH3D*				hPhiWeight_LHC15l_part1;	//! LHC15l, part1 runs
		TH3D*				hPhiWeight_LHC15l_part2;	//! LHC15l, part2 runs
		TH3D*				hPhiWeight_LHC15l_part3;	//! LHC15l, part3 runs

		TH3D*				hPhiWeightPlus_LHC15i_part1;	//! LHC15i, part1 runs
		TH3D*				hPhiWeightPlus_LHC15i_part2;	//! LHC15i, part2 runs
		TH3D*				hPhiWeightPlus_LHC15j_part1;	//! LHC15j, part1 runs
		TH3D*				hPhiWeightPlus_LHC15j_part2;	//! LHC15j, part2 runs
		TH3D*				hPhiWeightPlus_LHC15l_part1;	//! LHC15l, part1 runs
		TH3D*				hPhiWeightPlus_LHC15l_part2;	//! LHC15l, part2 runs
		TH3D*				hPhiWeightPlus_LHC15l_part3;	//! LHC15l, part3 runs

		TH3D*				hPhiWeightMinus_LHC15i_part1;	//! LHC15i, part1 runs
		TH3D*				hPhiWeightMinus_LHC15i_part2;	//! LHC15i, part2 runs
		TH3D*				hPhiWeightMinus_LHC15j_part1;	//! LHC15j, part1 runs
		TH3D*				hPhiWeightMinus_LHC15j_part2;	//! LHC15j, part2 runs
		TH3D*				hPhiWeightMinus_LHC15l_part1;	//! LHC15l, part1 runs
		TH3D*				hPhiWeightMinus_LHC15l_part2;	//! LHC15l, part2 runs
		TH3D*				hPhiWeightMinus_LHC15l_part3;	//! LHC15l, part3 runs

		// Event histograms
		TH1D*				hEventCount;								//! counting events passing given event cuts
		TH1F*				hMult;											//! multiplicity distribution
		TH1F*				hMultfBin[12]; 							//! multiplicity distribution in fBin
		TH1F*				fVtxAfterCuts;							//! Vertex z dist after cuts
		TH1F*				fCentralityDis;							//! distribution of centrality percentile using V0M estimator
		TH1F*				fV0CentralityDis;						//! distribution of V0M/<V0M>
		TH2F*				hMultV0vsNtrksAfterCuts;		//! Number of tracks vs. V0M/<V0M>
		TH2F*				hMultSPDvsNtrksAfterCuts;		//! Number of tracks vs. SPD/<SPD>
		TH2F*				hNtrksVSmultPercentile; 		//! Number of tracks vs. percentile using V0M estimator
		TH2F*				fCentralityV0MCL1;					//! LHC15o: V0M vs. CL1 percentile
		TH2F*				fCentralityV0MCL0;					//! LHC15o: V0M vs. CL0 percentile
		TH2F*				fCentralityCL0CL1;					//! LHC15o: CL0 vs. CL1 percentile
		TH2F*				fMultvsCentr;								//! LHC15o: Number of tracks vs. percentile
		TH2F*				fMult128vsCentr;						//! LHC15o: Number of FB128 tracks vs. percentile
		TH2F*				fMultTPCvsTOF;							//! LHC15o: Number of TPC tracks vs. ToF tracks
		TH2F*				fMultTPCvsESD;							//! LHC15o: Number of TPC tracks vs. ESD tracks

		TH2D*				hSPDClsVsTrk;								//! SPD clusters vs. tracklets without any cuts
		TH2D*				hV0C012vsTkl;								//! V0C mult. in 0,1,2nd ring vs. SPD tracklets without any cuts
		TH2D*				hV0C012vsV0C3;							//! V0C mult. in 0,1,2nd ring vs. V0C mult. in 3rd ring without any cuts
		TH2D*				hV0MOnVsOf;									//! V0M amplitude online vs. offline without any cuts
		TH2D*				hSPDOnVsOf;									//! SPD amplitude online vs. offline without anycuts


		// Track histograms
		TH1F*				fPhiDis1D;									//! phi dis 1D
		TH3F*				fPhiDis;										//! phi dist
		TH1F*				fEtaDis;										//! eta dist
		TH1F*				fEtaBefore;									//! eta dist before track cuts
		TH1F*				fPtDis;											//! pt dist
		TH1F*				fPtBefore;									//! pt dist before track cuts
		TH1F*				hDCAxyBefore; 							//!
		TH1F*				hDCAzBefore; 								//!
		TH1F*				hITSclustersBefore; 				//!
		TH1F*				hChi2Before; 								//!
		TH1F*				hDCAxy; 										//!
		TH1F*				hDCAz; 											//!
		TH1F*				hITSclusters; 							//!
		TH1F*				hChi2; 											//!

		// Correlation histograms
		TH2F*				hNtrksPt0530Pt0230; 				//! for correction of Xaxis of results with 0.5<pt<3.0 to 0.2<pt<3.0
		TH2F*				hNtrksPt0730Pt0230; 				//! for correction of Xaxis of results with 0.7<pt<3.0 to 0.2<pt<3.0
		TH2F*				hNtrksEta09Eta10;						//! for correction of Xaxis of results with |eta|<0.9 to |eta|<1.0
		TH2F*				hNtrksEta08Eta10;						//! for correction of Xaxis of results with |eta|<0.8 to |eta|<1.0
		TH2F*				hNtrksAllNtrksLS; 					//! for correction of Xaxis of LS results to all charged ptcls results
		TH2F*				hNtrksNoGapGap0;						//! for correction of Xaxis of Gap0 to NoGap results
		TH2F*				hNtrksNoGapGap2;						//! for correction of Xaxis of Gap2 to NoGap results
		TH2F*				hNtrksNoGapGap4;						//! for correction of Xaxis of Gap4 to NoGap results
		TH2F*               hNtrksNoGapGap6;                        //! for correction of Xaxis of Gap6 to NoGap results
		TH2F*				hNtrksNoGapGap8;						//! for correction of Xaxis of Gap8 to NoGap results
		TH2F*				hNtrksNoGapGap;							//! for correction of Xaxis of Gap10 to NoGap results
		TH2F*				hNtrksNoGapGap14;						//! for correction of Xaxis of Gap14 to NoGap results
		TH2F*				hNtrksNoGapGap16;						//! for correction of Xaxis of Gap16 to NoGap results
		TH2F*				hNtrksNoGapGap18;						//! for correction of Xaxis of Gap18 to NoGap results
		TH2F*				hNtrksNoGap3sub;						//! for correction of Xaxis of 3sub to NoGap results
		TH2F*				hNtrksNoGap3subGap;					//! for correction of Xaxis of 3sub Gap2 to NoGap results

		// Global variables
		int NtrksAfter = 0;
		int NtrksAfterGap10M = 0;
		int NtrksAfterGap10P = 0;
		int NtrksAfterGap14M = 0;
		int NtrksAfterGap14P = 0;
		int NtrksAfter3subL = 0;
		int NtrksAfter3subM = 0;
		int NtrksAfter3subR = 0;

		PhysicsProfile multProfile;
		PhysicsProfile multProfile_bin[10];

		//......
		// MC
		// Event histograms
		TH1F*				hMultMC;										//! Multiplicity distribution
		TH1F*				fPhiDisTruth;								//! Phi distribution
		TH1F*				fEtaDisTruth;								//! Eta distribution
		TH1F*				fPtDisTruth;								//! Pt distribution

		TH3F*				hReco;											//! pt, eta, Vz of reconstructed tracks
		TH3F*				hRecoPion;									//! pt, eta, Vz of reconstructed pions
		TH3F*				hRecoKaon;									//! pt, eta, Vz of reconstructed kaons
		TH3F*				hRecoProton;								//! pt, eta, Vz of reconstructed protons
		TH3F*				hRecoElectron;							//! pt, eta, Vz of reconstructed electrons
		TH3F*				hRecoMuon;									//! pt, eta, Vz of reconstructed muons
		TH3F*				hRecoLSplus;								//! pt, eta, Vz of positive reconstructed tracks
		TH3F*				hRecoLSminus;								//! pt, eta, Vz of negative reconstructed tracks

		TH2F*				hPtRecoNtrks;								//! reconstructed Pt vs. Number of generated particles
		TH2F*				hEtaRecoNtrks;							//! reconstructed Eta vs. Number of generated particles
		TH2F*				hVzRecoNtrks;								//! reconstructed Vz vs. Number of generated particles
		TH2F*				hPtRecoNtrksReco;						//! reconstructed Pt vs. Number of reconstructed tracks
		TH2F*				hEtaRecoNtrksReco;					//! reconstructed Eta vs. Number of reconstructed tracks
		TH2F*				hVzRecoNtrksReco;						//! reconstructed Vz vs. Number of reconstructed tracks

		TH3F*				hTruth;											//! pt, eta, Vz of generated particles
		TH3F*				hTruthPion;									//! pt, eta, Vz of generated pions
		TH3F*				hTruthKaon;									//! pt, eta, Vz of generated kaons
		TH3F*				hTruthProton;								//! pt, eta, Vz of generated protons
		TH3F*				hTruthElectron;							//! pt, eta, Vz of generated electrons
		TH3F*				hTruthMuon;									//! pt, eta, Vz of generated muons
		TH3F*				hTruthLSplus;								//! pt, eta, Vz of positive generated particles
		TH3F*				hTruthLSminus;							//! pt, eta, Vz of negative generated particles

		TH2F*				hPtTruthNtrks;							//! generated Pt vs. Number of generated particles
		TH2F*				hEtaTruthNtrks;							//! generated Eta vs. Number of generated particles
		TH2F*				hVzTruthNtrks;							//! generated Vs vs. Number of generated particles
		TH2F*				hPtTruthNtrksReco;					//! generated Pt vs. Number of reconstructed tracks
		TH2F*				hEtaTruthNtrksReco;					//! generated Eta vs. Number of reconstructed tracks
		TH2F*				hVzTruthNtrksReco;					//! generated Vz vs. Number of reconstructed tracks

		TH1D*				hPrimary;										//! number of primary particles
		TH1D*				hPions;											//! number of pions

		TH2F*				hDCAptMC;										//! DCA vs. Pt
		TH2F*				hDCAptMC_material;					//! DCA vs. Pt of decays from material interactions
		TH2F*				hDCAptMC_weak;							//! DCA vs. Pt of weak decays

		TH2F*				hNtrksRecoNtrksTruth; 			//! Number of generated particles vs. number of reconstructed tracks
		TH2F*				hNtrksRecoCorrNtrksTruth; 	//! Number of generated particles vs. number of reconstructed tracks weighted to match the difference between data and MC seen in LHC15ijl

		TH3F*				hRecoMult[8]; 							//!
		TH3F*				hTruthMult[8]; 							//!

		TH1F*				hTruthMultfBin[12]; 				//!

		TRandom3 rand;
		Int_t bootstrap_value;

		CorrelationCalculator correlator;

		double xbins[3000+10];
		int nn;
		void CalculateProfile(PhysicsProfile& profile, double Ntrks);
		void InitProfile(PhysicsProfile& profile, TString);

		ClassDef(AliAnalysisTaskNonlinearFlow, 1);    //Analysis task
};

#endif
// Local Variables:
//  mode: C++
// End:
