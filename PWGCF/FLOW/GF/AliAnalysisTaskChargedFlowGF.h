#ifndef ALIANALYSISTASKCHARGEDFLOWGF_H
#define ALIANALYSISTASKCHARGEDFLOWGF_H
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include <TComplex.h>
#ifndef __CINT__
// ROOT includes
// # include <TList.h>
// # include <TH1.h>
// # include <TH2.h>
// # include <TH3.h>
// # include <TProfile.h>
// # include <TComplex.h>
// # include <TBits.h> 
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


class AliAnalysisTaskChargedFlowGF : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskChargedFlowGF();
  AliAnalysisTaskChargedFlowGF(const char *name);

  virtual ~AliAnalysisTaskChargedFlowGF();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t* option);
  virtual void   Terminate(Option_t* ); 

	virtual void	SetFilterbit(Int_t bit){fFilterbit = bit;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetMinPt(Double_t minPt){fMinPt = minPt;}
  virtual void  SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
	virtual void	SetTPCclusters(Int_t tpcClus){fTPCclusters = tpcClus;}
	virtual void	SetUseDCAzCut(Bool_t usedcaz){fUseDCAzCut = usedcaz;}
	virtual void	SetDCAzCut(Double_t dcaz){fDCAz = dcaz;}
	virtual void	SetUseDCAxyCut(Bool_t usedcaxy){fUseDCAxyCut = usedcaxy;}
	virtual void	SetDCAxyCut(Double_t dcaxy){fDCAxy = dcaxy;}
  virtual void  SetIsSample(Int_t IsSample){fSample = IsSample;}
	virtual void	SetTrigger(Int_t trig){fTrigger = trig;}
	virtual void	SetNUEFlag(Bool_t NUE){fNUE = NUE;}
	virtual void	SetNUA(Bool_t NUA){fNUA = NUA;}
		//..MC
	virtual void	SetMCFlag(Bool_t mc){fIsMC = mc;}
		//....
	virtual void	SetPeriod(TString period){fPeriod = period;}
		
 private:
  AliAnalysisTaskChargedFlowGF(const AliAnalysisTaskChargedFlowGF&);
  AliAnalysisTaskChargedFlowGF& operator=(const AliAnalysisTaskChargedFlowGF&);

  virtual void		AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float fVtxZ);  
	virtual double	ProcessMCTruth(AliVEvent* aod, float VtxZ);

  double 					GetWeight(double phi, double eta, double vz, double runNumber);
	double 					GetPtWeight(double pt, double eta, float vz, double runNumber);

  TComplex Qvector[20][20];
  TComplex QvectorM[20][20];
  TComplex QvectorP[20][20];
  TComplex Qvector0M[20][20];
  TComplex Qvector0P[20][20];
  TComplex Qvector2M[20][20];
  TComplex Qvector2P[20][20];
  TComplex Qvector4M[20][20];
  TComplex Qvector4P[20][20];
  TComplex Qvector8M[20][20];
  TComplex Qvector8P[20][20];
  TComplex Qvector14M[20][20];
  TComplex Qvector14P[20][20];
	TComplex QvectorSubLeft[20][20];
	TComplex QvectorSubMiddle[20][20];
	TComplex QvectorSubRight[20][20];
	
  TComplex Q(int n, int p);
  TComplex QGapM(int n, int p);
  TComplex QGapP(int n, int p);
  TComplex QGap0M(int n, int p);
  TComplex QGap0P(int n, int p);
  TComplex QGap2M(int n, int p);
  TComplex QGap2P(int n, int p);
  TComplex QGap4M(int n, int p);
  TComplex QGap4P(int n, int p);
  TComplex QGap8M(int n, int p);
  TComplex QGap8P(int n, int p);
  TComplex QGap14M(int n, int p);
  TComplex QGap14P(int n, int p);
  TComplex QGap16M(int n, int p);
  TComplex QGap16P(int n, int p);
	TComplex QsubLeft(int n, int p);
	TComplex QsubMiddle(int n, int p);
	TComplex QsubRight(int n, int p);

  TComplex Two(int n1, int n2);
  TComplex TwoGap(int n1, int n2);
  TComplex TwoGap0(int n1, int n2);
  TComplex TwoGap2(int n1, int n2);
  TComplex TwoGap4(int n1, int n2);
  TComplex TwoGap8(int n1, int n2);
  TComplex TwoGap14(int n1, int n2);
  TComplex Two_3SubLM(int n1, int n2);
  TComplex Two_3SubRM(int n1, int n2);
  TComplex Two_3SubLM_perm2(int n1, int n2);
  TComplex Two_3SubLR_perm2(int n1, int n2);
  TComplex Two_3SubRM_perm3(int n1, int n2);
  TComplex Two_3SubRL_perm3(int n1, int n2);

  TComplex Three(int n1, int n2, int n3);
  TComplex ThreeGap0M(int n1, int n2, int n3);
  TComplex ThreeGap0P(int n1, int n2, int n3);
  TComplex Four(int n1, int n2, int n3, int n4);
	TComplex FourGap0(int n1, int n2, int n3, int n4);
	TComplex FourGap0M(int n1, int n2, int n3, int n4);
	TComplex FourGap0P(int n1, int n2, int n3, int n4);
	TComplex Four_3SubEvts(int n1, int n2, int n3, int n4);
	TComplex Four_3SubEvts_perm2(int n1, int n2, int n3, int n4);
	TComplex Four_3SubEvts_perm3(int n1, int n2, int n3, int n4);
	TComplex Five(int n1, int n2, int n3, int n4, int n5);
	TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6);
	TComplex SixGap0(int n1, int n2, int n3, int n4, int n5, int n6);
	TComplex Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7);
	TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
	TComplex EightGap0(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
  
	AliEventCuts	fEventCuts;					// Event cuts
  AliAODEvent* fAOD;                //! AOD object

  // Cuts and options
	Int_t				fFilterbit;					// flag for filter bit
	Double_t		fEtaCut;						// Eta cut used to select particles
	Double_t		fVtxCut;						// Vtx cut on z position in cm
	Double_t		fMinPt;							// Min pt - for histogram limits
	Double_t		fMaxPt;							// Max pt - for histogram limits
	Int_t				fTPCclusters;				// min. TPC clusters
	Bool_t			fUseDCAzCut;				// switch to choose whether I want to use DCAz cut or not (for systematic studies, otherwise it is in FB selection by default)
	Double_t		fDCAz;							// max DCAz, for systematics
	Bool_t			fUseDCAxyCut;				// the same switch as for DCAxy
	Double_t		fDCAxy;							// max DCAxy, for systematics
	Int_t				fSample;						// number of sample
	Int_t				fTrigger;						// flag for trigger
	Bool_t			fNUE;								// flag for NUE correction
	Bool_t			fNUA;								// 0: no NUA correction, 1: NUA correction
		//..MC
	Bool_t			fIsMC;							// flag to do MonteCarlo part
		//....
	TString			fPeriod;						// period
 
  // Output objects
	TList*			fListOfObjects;			//! Output list of objects
	TList*			fListOfObjectsMC;		//! Output list of objects for MonteCarlo

	// Cut functions for LHC15o
	TF1*				fMultTOFLowCut;			// cut low for TOF multiplicity outliers
	TF1*				fMultTOFHighCut;		// cut high for TOF multiplicity outliers
	TF1*				fMultCentLowCut;		// cut low for multiplicity centrality outliers

	// NUE
	TFile*			fTrackEfficiency;		//! file with tracking efficiency
	TH3F*				hTrackEfficiency;		//! histogram with tracking efficiency
	TH3F*				hTrackEfficiencyRun;//! histogram with tracking efficiency

	// NUA
	TFile*			fPhiWeight;					//! file with phi weights
	TH3F*				hPhiWeight;					//! 3D weight for all periods except LHC15ijl
	TH3F*				hPhiWeightRun;			//! 3D weight run-by-run for pPb 5TeV LHC16q

	// Event histograms
  TH1D*				hEventCount;								//! counting events passing given event cuts
  TH1F*				hMult;											//! multiplicity distribution
	TH1F*				fVtxAfterCuts;							//! Vertex z dist after cuts
	TH1F*				fCentralityDis;							//! distribution of centrality percentile using V0M estimator
	TH1F*				fV0CentralityDis;						//! distribution of V0M/<V0M>
  TH2F*				hMultV0vsNtrksAfterCuts;		//! Number of tracks vs. V0M/<V0M>
	TH2F*				hNtrksVSmultPercentile; 		//! Number of tracks vs. percentile using V0M estimator
	TH2F*				fCentralityV0MCL1;					//! LHC15o: V0M vs. CL1 percentile
	TH2F*				fCentralityV0MCL0;					//! LHC15o: V0M vs. CL0 percentile
	TH2F*				fCentralityCL0CL1;					//! LHC15o: CL0 vs. CL1 percentile
  TH2F*				fMultvsCentr;								//! LHC15o: Number of tracks vs. percentile
  TH2F*				fMult128vsCentr;						//! LHC15o: Number of FB128 tracks vs. percentile
	TH2F*				fMultTPCvsTOF;							//! LHC15o: Number of TPC tracks vs. ToF tracks
	TH2F*				fMultTPCvsESD;							//! LHC15o: Number of TPC tracks vs. ESD tracks

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

	// Physics profiles
	TProfile*		fChsc4242;									//! SC(4,2)
	TProfile*		fChsc4242Gap0;							//! SC(4,2) |#Delta#eta| > 0.0
	TProfile*		fChsc4242_3sub;							//! SC(4,2)_A 3subevent method
	TProfile*		fChsc4224_3sub;							//! SC(4,2)_B 3subevent method
	TProfile*		fChsc4242_3sub_perm2;				//! SC(4,2)_A 3subevent method
	TProfile*		fChsc4224_3sub_perm2;				//! SC(4,2)_B 3subevent method
	TProfile*		fChsc4242_3sub_perm3;				//! SC(4,2)_A 3subevent method
	TProfile*		fChsc4224_3sub_perm3;				//! SC(4,2)_B 3subevent method

	TProfile*		fChsc3232;									//! SC(3,2)
	TProfile*		fChsc3232Gap0;							//! SC(3,2) |#Delta#eta| > 0.0
	TProfile*		fChsc3232_3sub;							//! SC(3,2)_A 3subevent method
	TProfile*		fChsc3223_3sub;							//! SC(3,2)_B 3subevent method
	TProfile*		fChsc3232_3sub_perm2;				//! SC(3,2)_A 3subevent method
	TProfile*		fChsc3223_3sub_perm2;				//! SC(3,2)_B 3subevent method
	TProfile*		fChsc3232_3sub_perm3;				//! SC(3,2)_A 3subevent method
	TProfile*		fChsc3223_3sub_perm3;				//! SC(3,2)_B 3subevent method

	TProfile*		fsc4242[12];								//! the same but for different fBin (sampling)
	TProfile*		fsc4242Gap0[12];						//!
	TProfile*		fsc4242_3sub[12];						//!	
	TProfile*		fsc4224_3sub[12];						//!
	TProfile*		fsc4242_3sub_perm2[12];			//!	
	TProfile*		fsc4224_3sub_perm2[12];			//!
	TProfile*		fsc4242_3sub_perm3[12];			//!	
	TProfile*		fsc4224_3sub_perm3[12];			//!

	TProfile*		fsc3232[12];								//!
	TProfile*		fsc3232Gap0[12];						//!
	TProfile*		fsc3232_3sub[12];						//!
	TProfile*		fsc3223_3sub[12];						//!
	TProfile*		fsc3232_3sub_perm2[12];			//!
	TProfile*		fsc3223_3sub_perm2[12];			//!
	TProfile*		fsc3232_3sub_perm3[12];			//!
	TProfile*		fsc3223_3sub_perm3[12];			//!

	// Standard correlation profiles for different harmonics
  TProfile*		fChcn2Ntrks1bin[6];  				//! <<2>> in unit bins of Ntrks
  TProfile*		fChcn2Gap0Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.0
  TProfile*		fChcn2Gap2Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.2
  TProfile*		fChcn2Gap4Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.4
  TProfile*		fChcn2Gap8Ntrks1bin[6];  		//! <<2>> |#Delta#eta| > 0.8
  TProfile*		fChcn2GapNtrks1bin[6];  		//! <<2>> |#Delta#eta| > 1.0
  TProfile*		fChcn2Gap14Ntrks1bin[6];  	//! <<2>> |#Delta#eta| > 1.4

  TProfile*		fChcn2_3subLMNtrks1bin[6];  		//! <<2>> left vs. middle subevent
  TProfile*		fChcn2_3subRMNtrks1bin[6];  		//! <<2>> middle vs. right subevent
  TProfile*		fChcn2_3subLM_perm2Ntrks1bin[6];//! <<2>> left vs. middle subevent
  TProfile*		fChcn2_3subLR_perm2Ntrks1bin[6];//! <<2>> left vs. right subevent
  TProfile*		fChcn2_3subRM_perm3Ntrks1bin[6];//! <<2>> right vs. middle subevent
  TProfile*		fChcn2_3subRL_perm3Ntrks1bin[6];//! <<2>> right vs. left subevent

  TProfile*		fChcn4Ntrks1bin[6];  				//! <<4>> in unit bins of Ntrks
  TProfile*		fChcn4Gap0Ntrks1bin[6];  		//! <<4>> |#Delta#eta| > 0.0
	TProfile*		fChcn4_3subNtrks1bin[6]; 		//! <<4>> 3subevent method
	TProfile*		fChcn4_3sub_perm2Ntrks1bin[6];		//!
	TProfile*		fChcn4_3sub_perm3Ntrks1bin[6];		//!

  TProfile*		fChcn6Ntrks1bin[6];  				//! <<6>> in unit bins of Ntrks
  TProfile*		fChcn6Gap0Ntrks1bin[6];  		//! <<6>> |#Delta#eta| > 0.0

  TProfile*		fChcn8Ntrks1bin[6];  				//! <<8>> in unit bins of Ntrks
  TProfile*		fChcn8Gap0Ntrks1bin[6];  		//! <<8>> |#Delta#eta| > 0.0

	// the same profiles but for different fBin (sampling)
  TProfile*		fcn2Ntrks1bin[6][12];  			//!
  TProfile*		fcn2Gap0Ntrks1bin[6][12];		//!
  TProfile*		fcn2Gap2Ntrks1bin[6][12];		//!
  TProfile*		fcn2Gap4Ntrks1bin[6][12];		//!
  TProfile*		fcn2Gap8Ntrks1bin[6][12];		//!
  TProfile*		fcn2GapNtrks1bin[6][12];		//!
  TProfile*		fcn2Gap14Ntrks1bin[6][12];	//!

  TProfile*		fcn2_3subLMNtrks1bin[6][12];  		//!
  TProfile*		fcn2_3subRMNtrks1bin[6][12];  		//!
  TProfile*		fcn2_3subLM_perm2Ntrks1bin[6][12];//!
  TProfile*		fcn2_3subLR_perm2Ntrks1bin[6][12];//!
  TProfile*		fcn2_3subRM_perm3Ntrks1bin[6][12];//!
  TProfile*		fcn2_3subRL_perm3Ntrks1bin[6][12];//!

  TProfile*		fcn4Ntrks1bin[6][12];				//!
  TProfile*		fcn4Gap0Ntrks1bin[6][12];		//!
	TProfile*		fcn4_3subNtrks1bin[6][12];	//!
	TProfile*		fcn4_3sub_perm2Ntrks1bin[6][12];	//!
	TProfile*		fcn4_3sub_perm3Ntrks1bin[6][12];	//!

  TProfile*		fcn6Ntrks1bin[6][12];				//!
  TProfile*		fcn6Gap0Ntrks1bin[6][12];		//!

  TProfile*		fcn8Ntrks1bin[6][12];				//!
  TProfile*		fcn8Gap0Ntrks1bin[6][12];		//!

	//......
	// MC
	// Event histograms
	TH1F*				hMultMC;										//! Multiplicity distribution
	TH1F*				fPhiDisTruth;								//! Phi distribution
	TH1F*				fEtaDisTruth;								//! Eta distribution
	TH1F*				fPtDisTruth;								//! Pt distribution
	
	TH3F*				hReco;											//! pt, eta, Vz of reconstructed tracks

	TH3F*				hTruth;											//! pt, eta, Vz of generated particles

	TH1D*				hPrimary;										//! number of primary particles
	TH1D*				hPions;											//! number of pions

	TH2F*				hDCAptMC;										//! DCA vs. Pt
	TH2F*				hDCAptMC_material;					//! DCA vs. Pt of decays from material interactions
	TH2F*				hDCAptMC_weak;							//! DCA vs. Pt of weak decays

	TH2F*				hNtrksRecoNtrksTruth; 			//! Number of generated particles vs. number of reconstructed tracks

	// Physics profiles	
	TProfile*		fChMCsc4242;								//! same profiles as for real data (or reconstructed part of MC)
	TProfile*		fChMCsc4242Gap0;						//!
	TProfile*		fChMCsc4242_3sub;						//!
	TProfile*		fChMCsc4224_3sub;						//!
	TProfile*		fChMCsc3232;								//!
	TProfile*		fChMCsc3232Gap0;						//!
	TProfile*		fChMCsc3232_3sub;						//!
	TProfile*		fChMCsc3223_3sub;						//!

	TProfile*		fMCsc4242[12];							//!
	TProfile*		fMCsc4242Gap0[12];					//!
	TProfile*		fMCsc4242_3sub[12];					//!
	TProfile*		fMCsc4224_3sub[12];					//!
	TProfile*		fMCsc3232[12];							//!
	TProfile*		fMCsc3232Gap0[12];					//!
	TProfile*		fMCsc3232_3sub[12];					//!
	TProfile*		fMCsc3223_3sub[12];					//!

  TProfile*		fChMCcn2Ntrks1bin[6];  			//!
  TProfile*		fChMCcn2Gap0Ntrks1bin[6];		//!
  TProfile*		fChMCcn2Gap2Ntrks1bin[6];		//!
  TProfile*		fChMCcn2Gap4Ntrks1bin[6];		//!
  TProfile*		fChMCcn2Gap8Ntrks1bin[6];		//!
  TProfile*		fChMCcn2GapNtrks1bin[6];  	//!
  TProfile*		fChMCcn2Gap14Ntrks1bin[6];	//!

  TProfile*		fChMCcn2_3subLMNtrks1bin[6];			//!
  TProfile*		fChMCcn2_3subRMNtrks1bin[6];			//!

  TProfile*		fChMCcn4Ntrks1bin[6];				//!
  TProfile*		fChMCcn4Gap0Ntrks1bin[6];		//!
	TProfile*		fChMCcn4_3subNtrks1bin[6];	//!

  TProfile*		fChMCcn6Ntrks1bin[6];				//!
  TProfile*		fChMCcn6Gap0Ntrks1bin[6];		//!

  TProfile*		fChMCcn8Ntrks1bin[6];				//!
  TProfile*		fChMCcn8Gap0Ntrks1bin[6];		//!


  TProfile*		fMCcn2Ntrks1bin[6][12];			//!
  TProfile*		fMCcn2Gap0Ntrks1bin[6][12];	//!
  TProfile*		fMCcn2Gap2Ntrks1bin[6][12];	//!
  TProfile*		fMCcn2Gap4Ntrks1bin[6][12];	//!
  TProfile*		fMCcn2Gap8Ntrks1bin[6][12];	//!
  TProfile*		fMCcn2GapNtrks1bin[6][12];	//!
  TProfile*		fMCcn2Gap14Ntrks1bin[6][12];//!

  TProfile*		fMCcn2_3subLMNtrks1bin[6][12];			//!
  TProfile*		fMCcn2_3subRMNtrks1bin[6][12];			//!

  TProfile*		fMCcn4Ntrks1bin[6][12];			//!
  TProfile*		fMCcn4Gap0Ntrks1bin[6][12];	//!
	TProfile*		fMCcn4_3subNtrks1bin[6][12];//!

  TProfile*		fMCcn6Ntrks1bin[6][12];			//!
  TProfile*		fMCcn6Gap0Ntrks1bin[6][12];	//!

  TProfile*		fMCcn8Ntrks1bin[6][12];			//!
  TProfile*		fMCcn8Gap0Ntrks1bin[6][12];	//!

  ClassDef(AliAnalysisTaskChargedFlowGF, 1);    //Analysis task
};

#endif
// Local Variables: 
//  mode: C++
// End:
