#ifndef ALIANALYSISTASKNONLINEARFLOW_H
#define ALIANALYSISTASKNONLINEARFLOW_H
#include "AliAnalysisTaskSE.h"
#include "AliGFWCuts.h"
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
		private:
		ClassDef(PhysicsProfile, 4);    //Analysis task
};

class AliAnalysisTaskNonlinearFlow : public AliAnalysisTaskSE {
	public:

        enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kCharUnidentified, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest; NB: kUknown last as counter
        
        enum  NonflowSupress {knStandard = 1 << 0, kn0Gap = 1 << 1, knLargeGap = 1 << 2, knThreeSub = 1 << 3, knGapScan = 1 << 4}; // list of nonflow supression method

      

                // const unsigned int usev2345flag = 1 << 0;
	        // const unsigned int usev678flag = 1 << 1;

		AliAnalysisTaskNonlinearFlow();
		AliAnalysisTaskNonlinearFlow(const char *name);
		AliAnalysisTaskNonlinearFlow(const char *name, int NUA, int NUE);

		virtual ~AliAnalysisTaskNonlinearFlow();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t* option);
		virtual void   Terminate(Option_t* );

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
		virtual void   SetMinITSClusters(Int_t minClus){fMinITSClus = minClus;}
		virtual void   SetMaxChi2(Double_t maxChi){fMaxChi2 = maxChi;}
		virtual void   SetUseDCAzCut(Bool_t usedcaz){fUseDCAzCut = usedcaz;}
		virtual void   SetDCAzCut(Double_t dcaz){fDCAz = dcaz;}
		virtual void   SetDCAzCutDefault(Double_t dcaz){fDCAzDefault = dcaz;} // The dcaz cut for NtrksCounter
		virtual void   SetUseDCAxyCut(Bool_t usedcaxy){fUseDCAxyCut = usedcaxy;}
		virtual void   SetDCAxyCut(Double_t dcaxy){fDCAxy = dcaxy;}
		virtual void   SetDCAxyCutDefault(Double_t dcaxy){fDCAxyDefault = dcaxy;} // The dcaxy cut for NtrksCounter
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



	private:
		AliAnalysisTaskNonlinearFlow(const AliAnalysisTaskNonlinearFlow&);
		AliAnalysisTaskNonlinearFlow& operator=(const AliAnalysisTaskNonlinearFlow&);

		virtual void		AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus);
		virtual void            NTracksCalculation(AliVEvent* aod);
                Bool_t                  AcceptAOD(AliAODEvent *inEv);
                Bool_t                  AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp);
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
		Bool_t                  LoadWeightsKatarina();
		Bool_t                  LoadPtWeights();
		Bool_t                  LoadPtWeightsKatarina();
		Bool_t                  LoadWeightsSystematics();

                Double_t GetWeightKatarina(double phi, double eta, double vz);
                Double_t GetPtWeightKatarina(double pt, double eta, double vz);
		Double_t GetFlowWeight(const AliVParticle* track, double fVtxZ, const PartSpecies species);
		Double_t GetFlowWeightSystematics(const AliVParticle* track, double fVtxZ, const PartSpecies species);
                const char* ReturnPPperiod(const Int_t runNumber) const;
                const char* GetSpeciesName(const PartSpecies species) const;

		AliEventCuts	fEventCuts;					// Event cuts
                AliGFWCuts*     fGFWSelection;                                  //!
		AliAODEvent*    fAOD;                                           //! AOD object
		AliAODITSsaTrackCuts* fitssatrackcuts;                          //! itssatrackcuts object

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
		Int_t			fMinITSClus;				// min ITS clusters, LHC15ijl
		Double_t		fMaxChi2;				// max chi2 per ITS cluster, LHC15ijl
		Bool_t			fUseDCAzCut;				// switch to choose whether I want to use DCAz cut or not (for systematic studies, otherwise it is in FB selection by default)
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

		// Output objects
		TList*			fListOfObjects;			//! Output list of objects
		TList*			fListOfProfile;			//! Output list of objects
		TList*			fListOfProfiles[30];		//! Output list of objects

		// Cut functions for LHC15o
		TF1*			fMultTOFLowCut;			// cut low for TOF multiplicity outliers
		TF1*			fMultTOFHighCut;		// cut high for TOF multiplicity outliers
		TF1*			fMultCentLowCut;		// cut low for multiplicity centrality outliers

		// NUE
		TFile*			fTrackEfficiency;		//! file with tracking efficiency
		TH3F*			hTrackEfficiency;		//! histogram with tracking efficiency
		TH3F*			hTrackEfficiencyRun;            //! histogram with tracking efficiency

		// NUA
		bool fFlowRunByRunWeights;                              // flag of whether get the Run by run weight 
		bool fFlowPeriodWeights;                                // flag of whether to use period weight
		bool fFlowUse3Dweights;                                 // flag of whether to use 3d weight

		//
		TList*                  fFlowWeightsList;               //! flowWightsList
		TList*                  fFlowPtWeightsList;             //! PtflowWightsList
		TFile*                  fFlowPtWeightsFile;             //! PtflowWightsList
		TList*			fPhiWeight;	                //! file with phi weights
		TFile*			fPhiWeightFile;	                //! file with phi weights
		TList*			fPhiWeightPlus;	                //! file with phi weights
		TList*			fPhiWeightMinus;                //! file with phi weights
		TH2D*                   fh2Weights[kUnknown];           //! container for GF weights (phi,eta,pt) (2D)
		TH3D*                   fh3Weights[kUnknown];           //! container for GF weights (phi,eta,pt)
		TH2D*                   fh2AfterWeights[kUnknown];      //! distribution after applying GF weights - lightweight QA (phi)
                TH3D*                   fh3AfterWeights[kUnknown];      //! distribution after applying GF weights - full QA (phi,eta,pt)
		AliGFWWeights*          fWeightsSystematics;            //! Weights for systematics
		TH1D*                   fPtWeightsSystematics;          //! PtWeights for systematics


		TH3F*			hPhiWeight;			//! 3D weight for all periods except LHC15ijl
		TH3F*			hPhiWeightRun;			//! 3D weight run-by-run for pPb 5TeV LHC16q
		TH1F*			hPhiWeight1D;			//! 1D weight in one MC case (maybe need to redo to 3D weight)

		// Event histograms
		TH1D*			hEventCount;			//! counting events passing given event cuts
		TH1F*			hMult;				//! multiplicity distribution
		TH1F*			hMultfBin[12]; 			//! multiplicity distribution in fBin
		TH1F*			fVtxAfterCuts;			//! Vertex z dist after cuts
		TH1F*			fCentralityDis;			//! distribution of centrality percentile using V0M estimator
		TH1F*			fV0CentralityDis;		//! distribution of V0M/<V0M>
		TH2F*			hMultV0vsNtrksAfterCuts;	//! Number of tracks vs. V0M/<V0M>
		TH2F*			hMultSPDvsNtrksAfterCuts;	//! Number of tracks vs. SPD/<SPD>
		TH2F*			hNtrksVSmultPercentile; 	//! Number of tracks vs. percentile using V0M estimator
		TH2F*			fCentralityV0MCL1;		//! LHC15o: V0M vs. CL1 percentile
		TH2F*			fCentralityV0MCL0;		//! LHC15o: V0M vs. CL0 percentile
		TH2F*			fCentralityCL0CL1;		//! LHC15o: CL0 vs. CL1 percentile
		TH2F*			fMultvsCentr;	  		//! LHC15o: Number of tracks vs. percentile
		TH2F*			fMult128vsCentr;  		//! LHC15o: Number of FB128 tracks vs. percentile
		TH2F*			fMultTPCvsTOF;	  		//! LHC15o: Number of TPC tracks vs. ToF tracks
		TH2F*			fMultTPCvsESD;	  		//! LHC15o: Number of TPC tracks vs. ESD tracks

		TH2D*			hSPDClsVsTrk;			//! SPD clusters vs. tracklets without any cuts
		TH2D*			hV0C012vsTkl;			//! V0C mult. in 0,1,2nd ring vs. SPD tracklets without any cuts
		TH2D*			hV0C012vsV0C3;			//! V0C mult. in 0,1,2nd ring vs. V0C mult. in 3rd ring without any cuts
		TH2D*			hV0MOnVsOf;			//! V0M amplitude online vs. offline without any cuts
		TH2D*			hSPDOnVsOf;			//! SPD amplitude online vs. offline without anycuts


		// Track histograms
		TH1F*				fPhiDis1D;		//! phi dis 1D
		TH3F*				fPhiDis;		//! phi dist
		TH1F*				fEtaDis;		//! eta dist
		TH1F*				fEtaBefore;		//! eta dist before track cuts
		TH1F*				fPtDis;			//! pt dist
		TH1F*				fPtBefore;		//! pt dist before track cuts
		TH1F*				hDCAxyBefore; 		//!
		TH1F*				hDCAzBefore; 		//!
		TH1F*				hITSclustersBefore; 	//!
		TH1F*				hChi2Before; 		//!
		TH1F*				hDCAxy; 		//!
		TH1F*				hDCAz; 			//!
		TH1F*				hITSclusters; 		//!
		TH1F*				hChi2; 			//!

		// Global variables
		int NtrksCounter = 0;        //!
		int NtrksAfter = 0;          //!
		int NtrksAfterGap0M = 0;     //!
		int NtrksAfterGap0P = 0;     //!
		int NtrksAfterGap2M = 0;     //!
		int NtrksAfterGap2P = 0;     //!
		int NtrksAfterGap4M = 0;     //!
		int NtrksAfterGap4P = 0;     //!
		int NtrksAfterGap6M = 0;     //!
		int NtrksAfterGap6P = 0;     //!
		int NtrksAfterGap8M = 0;     //!
		int NtrksAfterGap8P = 0;     //!
		int NtrksAfterGap10M = 0;    //!
		int NtrksAfterGap10P = 0;    //!
		int NtrksAfterGap14M = 0;    //!
		int NtrksAfterGap14P = 0;    //!
		int NtrksAfter3subL = 0;     //!
		int NtrksAfter3subM = 0;     //!
		int NtrksAfter3subR = 0;     //!

                int lastRunNumber = 0;       //!

		PhysicsProfile multProfile; //!
		PhysicsProfile multProfile_bin[30]; //!

		CorrelationCalculator correlator; //!
		TRandom3 rand;         //!
		Int_t bootstrap_value = -1; //!


                unsigned fgFlowHarmonics = 0;        //! calculate v2, v3, v4, v5
                unsigned fgFlowHarmonicsHigher = 0;  //! calculate v6, v7, v8 ..
                unsigned fgFlowHarmonicsMult = 0;    //! calculate v2{4} // yet v2{6}, v2{8}
                unsigned fgNonlinearFlow = 0;        //! calculate v_4,22, v_5,32
                unsigned fgSymmetricCumulants = 0;   //! calculate SC(3,2), SC(4,2)

                unsigned fgTwoParticleCorrelation = 0;       //!
                unsigned fgTwoParticleCorrelationHigher = 0; //!
                unsigned fgThreeParticleCorrelation = 0;     //!
                unsigned fgFourParticleCorrelation = 0;      //! 

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

                bool fuQStandard = 0; //!
                bool fuQ0Gap     = 0; //!
                bool fuQLargeGap = 0; //!
                bool fuQThreeSub = 0; //!
                bool fuQGapScan  = 0; //!

		double xbins[300] = {}; //!
		int nn = 0; //!
		void CalculateProfile(PhysicsProfile& profile, double Ntrks);
		void InitProfile(PhysicsProfile& profile, TString name, TList* listOfProfile);

		ClassDef(AliAnalysisTaskNonlinearFlow, 4);    //Analysis task
};

#endif
// Local Variables:
//  mode: C++
// End:
