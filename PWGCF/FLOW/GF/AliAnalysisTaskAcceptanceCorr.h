#ifndef ALIANALYSISTASKACCEPTANCECORR_H
#define ALIANALYSISTASKACCEPTANCECORR_H

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

class AliAnalysisTaskAcceptanceCorr : public AliAnalysisTaskSE {
	public:

		enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kCharUnidentified, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest; NB: kUknown last as counter

		AliAnalysisTaskAcceptanceCorr();
		AliAnalysisTaskAcceptanceCorr(const char *name);
		AliAnalysisTaskAcceptanceCorr(const char *name, int NUA, int NUE, TString fPeriod);

		virtual ~AliAnalysisTaskAcceptanceCorr();

		virtual void  UserCreateOutputObjects();
		virtual void UserExec(Option_t *option);
		virtual void NotifyRun();
		virtual void Terminate(Option_t *);

		virtual void   SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
		virtual void   SetMinPt(Double_t minPt){fMinPt = minPt;}
		virtual void   SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
	
		virtual void   SetNtrksName(TString ntrksname){fNtrksName = ntrksname;}
		virtual void   SetAdditionalTPCPileupCuts(bool flag = true) {fAddTPCPileupCuts = flag;}
		virtual void   SetESDvsTPConlyLinearCut(double cut = 15000) {fESDvsTPConlyLinearCut = cut;}
		virtual void   SetTPCchi2perCluster(double fchi2 = 4) {fTPCchi2perCluster = fchi2;}
		virtual void   SetTrigger(Int_t trig){fTrigger = trig;}

		virtual void   SetPeriod(TString period){fPeriod = period;}
		virtual void   SetSystFlag(int syst){fCurrSystFlag = syst;}
		virtual int    GetSystFlag(){return fCurrSystFlag;}

        void SetUseFMDcut(Bool_t cut = kTRUE) { fUseFMDcut = cut; }
        void SetFMDcutParameters(Double_t par0a, Double_t par1a, Double_t par0c, Double_t par1c) { fFMDcutapar0 = par0a; fFMDcutapar1 = par1a; fFMDcutcpar0 = par0c; fFMDcutcpar1 = par1c; }
		virtual void   SetFMDacceptanceCuts(double AL, double AH, double CL, double CH) { fFMDAacceptanceCutLower = AL; fFMDAacceptanceCutUpper = AH; fFMDCacceptanceCutLower = CL; fFMDCacceptanceCutUpper = CH; }

	private:
		AliAnalysisTaskAcceptanceCorr(const AliAnalysisTaskAcceptanceCorr&);
		AliAnalysisTaskAcceptanceCorr& operator=(const AliAnalysisTaskAcceptanceCorr&);

                Bool_t PrepareTPCFMDTracks();
		Bool_t AcceptAOD(AliAODEvent *inEv);
		Bool_t AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp);

		AliEventCuts	fEventCuts;					// Event cuts
		AliGFWCuts*     fGFWSelection;                                  //!
		AliGFWNFCuts*   fGFWSelection15o;                               //!
		AliAODEvent*    fAOD;                                           //! AOD object

		// Cuts and options
		Double_t		fEtaCut;				// Eta cut used to select particles
		Double_t		fMinPt;					// Min pt - for histogram limits
		Double_t		fMaxPt;					// Max pt - for histogram limits
		Int_t			fTrigger;				// flag for trigger
		Int_t			fAliTrigger;				// name for trigger
		// Bool_t		fLS;					// charge, 1:all, 2:pp,  3: mm
		Bool_t			fNUE;					// flag for NUE correction
		Bool_t			fNUA;					// 0: no NUA correction, 1: NUA correction
		bool                    fIsMC;                                  // The observable for MonteCarlo truth
		TString                 fNtrksName;                             // Cent or Mult
		TString                 anaType;                                // TPC-TPC or TPC-FMD
		TString					fPeriod;								// period
		Int_t                   fCurrSystFlag;                          // Systematics flag
		Bool_t                  fAddTPCPileupCuts;                      // Additional TPC pileup cuts

  Int_t                   fESDvsTPConlyLinearCut;                   // ESD-TPC linear cut
		Double_t                fTPCchi2perCluster;                     // Additional cuts for TPC chi2 / cluster
		Bool_t                  fUseAdditionalDCACut;                   // Additianal cuts for dca: < 1 cm

        Bool_t                  fUseFMDcut;                             // [kTRUE]
        Double_t                fFMDcutapar0;                           // [1.64755]
        Double_t                fFMDcutapar1;                           // [119.602]
        Double_t                fFMDcutcpar0;                           // [2.73426]
        Double_t                fFMDcutcpar1;                           // [150.31]
        Double_t                fFMDAacceptanceCutLower;                // FMDCut
        Double_t                fFMDAacceptanceCutUpper;                // FMDCut
        Double_t                fFMDCacceptanceCutLower;                // FMDCut
        Double_t                fFMDCacceptanceCutUpper;                // FMDCut

		// Output objects
		TList*			fListOfObjects;			//! Output list of objects
  TH2D* hEtaPhiDist;               // eta-phi distribution;

		// Cut functions for LHC15o
		TF1*			fMultTOFLowCut;			// cut low for TOF multiplicity outliers
		TF1*			fMultTOFHighCut;		// cut high for TOF multiplicity outliers
		TF1*			fMultCentLowCut;		// cut low for multiplicity centrality outliers


		// Event histograms
		TH1D*			hEventCount;			//! counting events passing given event cuts
		TH1F*			hMult;				//! multiplicity distribution
		TH1F*			fVtxAfterCuts;			//! Vertex z dist after cuts
		TH1F*			fCentralityDis;			//! distribution of centrality percentile using V0M estimator

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


		double fPVz;                   //!
		double fCentrality;            //!

		ClassDef(AliAnalysisTaskAcceptanceCorr, 0); // Analysis task
};

#endif
