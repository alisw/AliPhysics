#ifndef ALIANALYSISTASKCORRFORNONLINEARFLOW_H
#define ALIANALYSISTASKCORRFORNONLINEARFLOW_H

#include "AliAnalysisTaskSE.h"
#include "AliGFWCuts.h"
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
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

#include "AliEventPoolManager.h"
#include "AliTHn.h"

#include "AliAODForwardMult.h"
#include "AliPartSimpleForCorr.h"

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

class AliAnalysisTaskCorrForNonlinearFlow : public AliAnalysisTaskSE {
public:

  enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kCharUnidentified, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest; NB: kUknown last as counter

  AliAnalysisTaskCorrForNonlinearFlow();
  AliAnalysisTaskCorrForNonlinearFlow(const char *name);
  AliAnalysisTaskCorrForNonlinearFlow(const char *name, int NUA, int NUE, TString fPeriod);

  virtual ~AliAnalysisTaskCorrForNonlinearFlow();

  virtual void  UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
  virtual void Terminate(Option_t *);

  virtual void   SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void   SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void   SetVtxCutDefault(Double_t vtxCut){fVtxCutDefault = vtxCut;} // The vtxCut for NtrksCounter
  virtual void   SetMinPt(Double_t minPt){fMinPt = minPt;}
  virtual void   SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
  virtual void   SetMinPtAss(Double_t minPt){fPtMinAss = minPt;}
  virtual void   SetMaxPtAss(Double_t maxPt){fPtMaxAss = maxPt;}
  virtual void   SetMinPtTrig(Double_t minPt){fPtMinTrig = minPt;}
  virtual void   SetMaxPtTrig(Double_t maxPt){fPtMaxTrig = maxPt;}
  virtual void   SetMaxCent(Double_t maxCent){fCentMax = maxCent;}
  virtual void   SetMinCent(Double_t maxCent){fCentMin = maxCent;}
  
  virtual void   SetIsSample(Int_t IsSample){fSample = IsSample;}
  virtual void   SetIsMC(Bool_t isMC, Bool_t useTPCTruth, Bool_t useFMDTruth){fIsMC = isMC; fUseTPCTruth = useTPCTruth; fUseFMDTruth = useFMDTruth;}
  virtual void   SetNtrksName(TString ntrksname){fNtrksName = ntrksname;}
  virtual void   SetUseWeigthsRunByRun(Bool_t bRunByRun = kTRUE) { fFlowRunByRunWeights = bRunByRun; }
  virtual void   SetUsePeriodWeigths(Bool_t weight = kTRUE) { fFlowPeriodWeights = weight; }
  virtual void   SetUseWeights3D(Bool_t use = kTRUE) { fFlowUse3Dweights = use; }
  virtual void   SetSpringMode(bool flag = true) { fSpringMode = flag; }
  virtual void   SetLowMultiplicityMode(bool flag = true) {fLowMultiplicityMode = flag;}
  virtual void   SetAdditionalTPCPileupCuts(bool flag = true) {fAddTPCPileupCuts = flag;}
  virtual void   SetESDvsTPConlyLinearCut(double cut = 15000) {fESDvsTPConlyLinearCut = cut;}
  virtual void   SetUseCorrectedNTracks(bool flag = true) {fUseCorrectedNTracks = flag;}
  virtual void   SetUseFlippedEta(bool flag = true) {fUseFlippedEta = flag;}
  virtual void   SetUseNarrowBin(bool flag = true) {fUseNarrowBin = flag;}
  virtual void   SetExtremeEfficiency(int flag = 0) {fExtremeEfficiency = flag;}
  virtual void   SetTPCchi2perCluster(double fchi2 = 4) {fTPCchi2perCluster = fchi2;}
  // virtual void   SetUseAdditionalDCACut(double flag = true) {fUseAdditionalDCACut = flag;}
  // virtual void   SetUseDefaultWeight(double flag = true) {fUseDefaultWeight = flag;}
  virtual void   SetAnaType(TString type) {anaType = type;}
  virtual void   SetTrigger(Int_t trig){fTrigger = trig;}
  virtual void   SetNUEFlag(Bool_t NUE){fNUE = NUE;}
  virtual void   SetNUA(Bool_t NUA){fNUA = NUA;}

  virtual void   SetPeriod(TString period){fPeriod = period;}
  virtual void   SetSystFlag(int syst){fCurrSystFlag = syst;} 
  virtual int    GetSystFlag(){return fCurrSystFlag;}
  virtual void   UseBootstrap(bool ftest = true){fBootstrapStat = ftest;}
  virtual void   SetBootstrapRange(int low=-100, int high = 100){sampleLow = low; sampleHigh = high;}
  virtual void   SetFastMode(bool mode = true){fFastMode = mode;}
  virtual void   SetAlternativePhiBinning(bool mode = true){fAlternativePhiBinning = mode;}

  void SetUsePhiStarCut(Bool_t cut = kTRUE) { fUsePhiStarCut = cut; }
  void SetUseFMDcut(Bool_t cut = kTRUE) { fUseFMDcut = cut; }
  void SetFMDcutParameters(Double_t par0a, Double_t par1a, Double_t par0c, Double_t par1c) { fFMDcutapar0 = par0a; fFMDcutapar1 = par1a; fFMDcutcpar0 = par0c; fFMDcutcpar1 = par1c; }
  virtual void   SetFMDacceptanceCuts(double AL, double AH, double CL, double CH) { fFMDAacceptanceCutLower = AL; fFMDAacceptanceCutUpper = AH; fFMDCacceptanceCutLower = CL; fFMDCacceptanceCutUpper = CH; }
  virtual void   SetnSamples(int nsamples) {nSamples = nsamples;}
  virtual void   SetBinningMethod(unsigned method) {fBinMethod = method;}
  virtual void   SetMixParameter(int n1, int n2, int n3) {fPoolMaxNEvents = n1; fPoolMinNTracks = n2; fMinEventsToMix = n3;}
  Double_t RangePhi(Double_t DPhi);
  Double_t  RangePhiFMD(Double_t DPhi);
  Double_t GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius);

private:
  AliAnalysisTaskCorrForNonlinearFlow(const AliAnalysisTaskCorrForNonlinearFlow&);
  AliAnalysisTaskCorrForNonlinearFlow& operator=(const AliAnalysisTaskCorrForNonlinearFlow&);

  virtual void NTracksCalculation(AliVEvent* aod);
  Bool_t PrepareTPCFMDTracks();
  Bool_t PrepareTPCFMDTracksMCTruth();
  virtual void FillCorrelations();
  virtual void FillCorrelationsMixed();
  Bool_t AcceptAOD(AliAODEvent *inEv);
  Bool_t AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp);
  Bool_t AcceptMCTruthTrack(AliAODMCParticle *mtrk);

  double 			GetWeight(double phi, double eta, double pt, int run, bool fPlus, double vz, double runNumber);
  double GetPtWeight(double pt, double eta, float vz, double runNumber);

  Bool_t                  LoadWeightsKatarina();
  Bool_t                  LoadPtWeights();
  Bool_t                  LoadPtWeightsKatarina();
  Bool_t                  LoadWeightsSystematics();

  Double_t GetWeightKatarina(double phi, double eta, double vz);
  Double_t GetPtWeightKatarina(double pt, double eta, double vz);
  Double_t GetFlowWeight(const AliVParticle* track, double fVtxZ, const PartSpecies species);
  Double_t GetFlowWeightSystematics(const AliVParticle* track, double fVtxZ, const PartSpecies species);
  const char* ReturnPPperiod(const Int_t runNumber) const;
  const char* ReturnPPperiodMC(const Int_t runNumber) const;
  const char* GetSpeciesName(const PartSpecies species) const;

  AliEventCuts	fEventCuts;					// Event cuts
  AliGFWCuts*     fGFWSelection;                                  //!
  AliGFWNFCuts*   fGFWSelection15o;                               //!
  AliAODEvent*    fAOD;                                           //! AOD object
  AliAODITSsaTrackCuts* fitssatrackcuts;                          //! itssatrackcuts object

  // Cuts and options
  Double_t		fEtaCut;				// Eta cut used to select particles
  Double_t		fVtxCut;				// Vtx cut on z position in cm
  Double_t		fVtxCutDefault;				// Vtx cut on z position in cm (for NtrksCounter)
  Double_t		fMinPt;					// Min pt - for histogram limits
  Double_t		fMaxPt;					// Max pt - for histogram limits
  Double_t                fPtMinAss;                              // Min pt - for Associate particle 
  Double_t                fPtMaxAss;                              // Min pt - for Associate particle 
  Double_t                fPtMinTrig;                             // Min pt - for trigger particle 
  Double_t                fPtMaxTrig;                             // Min pt - for trigger particle 
  Double_t                fCentMin;                               // Min Centrality
  Double_t                fCentMax;                               // Max Centrality
  Int_t			fSample;				// number of sample
  Int_t			fTrigger;				// flag for trigger
  Int_t			fAliTrigger;				// name for trigger
  // Bool_t		fLS;					// charge, 1:all, 2:pp,  3: mm
  Bool_t			fNUE;					// flag for NUE correction
  Bool_t			fNUA;					// 0: no NUA correction, 1: NUA correction
  bool                    fIsMC;                                  // The observable for MonteCarlo truth
  bool                    fUseTPCTruth;                           // Use the MC truth particles of TPC
  bool                    fUseFMDTruth;                           // Use the MC truth particles of TPC
  TString                 fNtrksName;                             // Cent or Mult
  TString                 anaType;                                // TPC-TPC or TPC-FMD
  TString					fPeriod;								// period
  Int_t                   fCurrSystFlag;                          // Systematics flag
  Bool_t                  fSpringMode;                            // The mode with spring cuts.
  Bool_t                  fLowMultiplicityMode;                   // The mode to consider low-multiplicity region 
  Bool_t                  fAddTPCPileupCuts;                      // Additional TPC pileup cuts
  Double_t                fESDvsTPConlyLinearCut;                 // ESDvsTPConlyLinearCut : default = 15000
  Bool_t                  fUseCorrectedNTracks;                   // Use corrected Ntracks in the filling of xbins;
  Bool_t                  fUseFlippedEta;                         // Flip the eta region to merge the pPb and Pbp sample;
  Bool_t                  fUseNarrowBin;                          // Use Narrow bin
  Int_t                   fExtremeEfficiency;                     // The flag to set extreme efficiency
  Double_t                fTPCchi2perCluster;                     // Additional cuts for TPC chi2 / cluster
  Bool_t                  fUseAdditionalDCACut;                   // Additianal cuts for dca: < 1 cm
  Bool_t                  fUseDefaultWeight;                      // Force to use the default weight 
  Double_t                fEtaGap3Sub;                            // The Eta Gap for 3 sub sample, the default is 0.4
  Int_t                   fPoolMaxNEvents;                        // Maximum number of events in a pool
  Int_t                   fPoolMinNTracks;                        // Minimum number of tracks to mix
  Int_t                   fMinEventsToMix;                        // Minimum numver of events to mix
  Bool_t                  fBootstrapStat;                         // Flag to calculate statistical uncertainty with bootstrap
  Bool_t                  fUsePhiStarCut;                         // [kTRUE]
  Bool_t                  fUseFMDcut;                             // [kTRUE]
  Double_t                fFMDcutapar0;                           // [1.64755]
  Double_t                fFMDcutapar1;                           // [119.602]
  Double_t                fFMDcutcpar0;                           // [2.73426]
  Double_t                fFMDcutcpar1;                           // [150.31]
  Double_t                fFMDAacceptanceCutLower;                // FMDCut
  Double_t                fFMDAacceptanceCutUpper;                // FMDCut
  Double_t                fFMDCacceptanceCutLower;                // FMDCut
  Double_t                fFMDCacceptanceCutUpper;                // FMDCut
  unsigned                fBinMethod;                             // fBinMethod
  int                     nSamples;                               // Number of bootstrap samples
  int                     sampleLow;                              // lower bound of bootstrap sample
  int                     sampleHigh;                             // higher bound of bootstrap sample
  Bool_t                  fFastMode;                                // Clear the deque after correlation;
  Bool_t                  fAlternativePhiBinning;                   // Alternative Phi Binning method

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
  TList*                  fFlowFeeddownList;              //! FeeddownList
  TFile*                  fFlowPtWeightsFile;             //! PtflowWightsList
  TList*			fPhiWeight;	                //! file with phi weights
  TFile*			fPhiWeightFile;	                //! file with phi weights
  TList*			fPhiWeightPlus;	                //! file with phi weights
  TList*			fPhiWeightMinus;                //! file with phi weights
  TH2D*                   hWeight2D;                      //! phi weights 2d 80 bins
  TH2D*                   fh2Weights[kUnknown];           //! container for GF weights (phi,eta,pt) (2D)
  TH3D*                   fh3Weights[kUnknown];           //! container for GF weights (phi,eta,pt)
  TH2D*                   fh2AfterWeights[kUnknown];      //! distribution after applying GF weights - lightweight QA (phi)
  TH3D*                   fh3AfterWeights[kUnknown];      //! distribution after applying GF weights - full QA (phi,eta,pt)
  AliGFWWeights*          fWeightsSystematics;            //! Weights for systematics
  TH1D*                   fPtWeightsSystematics;          //! PtWeights for systematics
  TH1D*                   fPtWeightsFeeddown;             //! Feeddown for systematics


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
  TH1D*				fPhiDis1D;		    //! phi dis 1D
  TH1D*				fPhiDis1DBefore;    //! phi dis 1D before track cuts
  TH3D*				fPhiDis;		    //! phi dist
  TH1D*				fEtaTriDis;		    //! eta dist
  TH1D*				fEtaTriDisBefore;	//! eta dist before track cuts
  TH1D*				fEtaAssDis;		    //! eta dist
  TH1D*				fEtaAssDisBefore;	//! eta dist before track cuts
  TH1D*				fPhiTriDis;		    //! phi dist
  TH1D*				fPhiTriDisBefore;	//! phi dist before track cuts
  TH1D*				fPhiAssDis;		    //! phi dist
  TH1D*				fPhiAssDisBefore;	//! phi dist before track cuts
  TH1D*				fPtTriDis;			//! pt dist
  TH1D*				fPtTriDisBefore;	//! pt dist before track cuts
  TH1D*				fPtAssDis;			//! pt dist
  TH1D*				fPtAssDisBefore;	//! pt dist before track cuts
  TH1F*				hDCAxyBefore; 		//!
  TH1F*				hDCAzBefore; 		//!
  TH1F*				hITSclustersBefore; //!
  TH1F*				hChi2Before; 		//!
  TH1F*				hDCAxy; 	       	//!
  TH1F*				hDCAz; 			    //!
  TH1F*				hITSclusters; 		//!
  TH1F*				hChi2; 			    //!
  TH2D*               hFMDAvsV0;          //!
  TH2D*               hFMDCvsV0;          //!


    
  TObjArray*                   fTracksTrigCharged;     //! List of charged tracks
  TObjArray*                   fTracksAss;             //! List of associate tracks
  AliTHn*                      fhChargedSE;            //!
  AliTHn*                      fhChargedME;            //!
  AliEventPoolManager*         fPoolMgr;               //!  event pool manager for Event Mixing
  TH2D*                        fhTracksTrigCent;       //! Trigger particle histogram

  // Global variables
  TRandom3 rand;                 //!
  double NtrksCounter = 0;       //!
  double NTracksCorrected = 0;   //!
  double NTracksUncorrected = 0; //!
  int NtrksAfter = 0;            //!

  int bootstrap_value;           //!
  int lastRunNumber   = 0;       //!
  double fPVz;                   //!
  double fCentrality;            //!
  Double_t fbSign;               //!

  ClassDef(AliAnalysisTaskCorrForNonlinearFlow, 11); // Analysis task
};

#endif
