#ifndef AliAnaTaskSEXic0SL_H
#define AliAnaTaskSEXic0SL_H

#include "AliAnalysisTaskSE.h"

class AliAODcascade;
class AliAODMCParticle;
class AliAODTrack;
class AliESDtrackCuts;
class AliInputEventHandler;
class AliMCEvent;
class AliMultSelection;
class AliNormalizationCounter;
class AliPIDResponse;
class AliRDHFCutsXictoeleXifromAODtracks;
class AliVEvent;
class AliVVertex;
class THistManager;

class TString;
class TTree;

class AliAnaTaskSEXic0SL : public AliAnalysisTaskSE
{
	public:

		AliAnaTaskSEXic0SL();
		AliAnaTaskSEXic0SL(const char* name, const char* option);
		virtual ~AliAnaTaskSEXic0SL();

		virtual void Terminate(Option_t* option);
		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t* option);

		//Options invoked via AddTask
		void SetMC(bool set) { IsMC = set; };
		void SetPA(bool set) { IsPA = set; };
		void SetTrigMB(bool set) { TrigOnMB = set; }
		void SetTrigHMV0(bool set) { TrigOnHMV0 = set; }

		//User defined functions
		int CheckOrigin(AliMCEvent* MCEvt, AliAODMCParticle *MCPart, bool SearchUpToQuark); //<=0:err, 4:c, 5:b
		bool FilterTrack(AliAODTrack* Trk, const AliVVertex* Vtx);
		bool FilterTrackElectron(AliAODTrack* Trk, AliPIDResponse* PID);
		bool FilterCascade(AliAODcascade* Casc, const AliVVertex* Vtx, AliPIDResponse* PID);

		void ControlAnaObjects(int option); //0: set, 1: delete
		void ControlOutputContainers(int option); //0: define, 1: post
		void ControlOutputTree(TTree* T, bool isMC, bool readOnly = false);
		void ResetTreeVariables(void);

		void ExecuteOffSel(TTree* T, TString Cond); //For offline selection of Grid/LT output

	private:

		//Analysis objects
		AliAODcascade*                      fCasc = nullptr;
		AliAODMCParticle*                   fMCPart = nullptr;
		AliAODTrack*                        fTrk = nullptr;
		AliESDtrackCuts*                    fTrkCuts = nullptr;
		AliInputEventHandler*               fInputHandler = nullptr;
		AliMCEvent*                         fMCEvt  = nullptr;
		AliMultSelection*                   fMultSel = nullptr;
		AliNormalizationCounter*            fANC_MB_0to100 = nullptr;
		AliNormalizationCounter*            fANC_MB_30to100 = nullptr;
		AliNormalizationCounter*            fANC_MB_0p1to30 = nullptr;
		AliNormalizationCounter*            fANC_HMV0_0to0p1 = nullptr;
		AliNormalizationCounter*            fANC_INEL0_MB_0to100 = nullptr;
		AliNormalizationCounter*            fANC_INEL0_MB_30to100 = nullptr;
		AliNormalizationCounter*            fANC_INEL0_MB_0p1to30 = nullptr;
		AliNormalizationCounter*            fANC_INEL0_HMV0_0to0p1 = nullptr;
		AliPIDResponse*                     fPID = nullptr;
		AliRDHFCutsXictoeleXifromAODtracks* fEvtCutsMB = nullptr;
		AliRDHFCutsXictoeleXifromAODtracks* fEvtCutsHMV0 = nullptr;
		AliVEvent*                          fEvt = nullptr;
		const AliVVertex*                   fEvtVtx = nullptr;
		THistManager*                       fHisto = nullptr;
		TTree*                              fTree = nullptr;

		//Options
		TString fTaskOpt;
		bool IsMC = false;
		bool IsPA = false;
		bool TrigOnMB = false;
		bool TrigOnHMV0 = false;
		std::vector<UInt_t> TrigStore; //Container for triggers of interest

		//-----------------------------------------------------------

		//Constants
		static const int MaxNXic0       = 20; //max. # of Xic0 per event (MC truth)
		static const int MaxNEle        = 20; //After all cut: typically 1-2 per event
		static const int MaxNCasc       = 20;
		static const int PDGCode_e      = 11;
		static const int PDGCode_Lambda = 3122; //Lambda0, 1115.683 +- 0.006 (MeV)
		static const int PDGCode_Omega  = 3334; //Omega-, 1672.45 +- 0.29 (MeV)
		static const int PDGCode_Xi     = 3312; //Xi- (strange baryon), 1321.71 +- 0.07 (MeV)
		static const int PDGCode_Xic0   = 4132; //2470.87 +0.28 -0.31 (MeV)
		const double MassEle = 0.51099895 * 1.E-3; //Electron mass in GeV
		const double MassLmb = 1115.683 * 1.E-3; //Lambda0 mass in GeV
		const double MassXi  = 1321.71 * 1.E-3; //Xi mass in GeV

		//Eventwise cut
		const Int_t   cut_runNoLo = 252000;
		const Int_t   cut_runNoUp = 295000;
		const Int_t   cut_vtxNContributors = 1;
		const Float_t cut_bfield = 0.001;
		const Float_t cut_eta    = 1.0;
		const Float_t cut_vtxZ   = 10.0;

		//Trackwise cut, common
		const Int_t   cut_minNClustersITS = 2;
		const Int_t   cut_TPCsignalN      = 50; //fSetProdTrackTPCNclsPID in old code
		const Float_t cut_maxChi2PerClusterITS = 36.;
		const Float_t cut_maxDCAToVertexXY     = 1.0;
		const Float_t cut_maxDCAToVertexZ      = 2.0;
		const Float_t cut_trkEta               = 0.8; //For daughter particles
		const Float_t cut_trkPt                = 0.5; //Lower limit of electron

		//Trackwise cut, electron
		const Float_t cutEle_massConv         = 0.05; //GeV, max. conversion mass
		const Float_t cutEle_nSigmaTOFAbs     = 3.0;
		const Float_t cutEle_nSigmaTPCAbsConv = 5.0; //Used only to check conversion mass (FilterTrackElectron)
		const Float_t cutEle_nSigmaTPCMax     = 3.0; //Cf. Min: pT dependent (FilterTrackElectron)

		//Trackwise cut, cascade
		const Float_t cutCasc_massTolLambda = 0.008;
		const Float_t cutCasc_massTolOmega  = 0.3 * 1.E-3; //ckim, NOT used for now
		const Float_t cutCasc_massTolXi     = 0.03; //!! Old code: online 0.01, offline selection 0.008

		const Float_t cutCasc_nSigmaTPCAbs    = 4.0;
		const Float_t cutCasc_minDecayLenXi   = 0.2; //Decay length btw PV to cascade
		const Float_t cutCasc_minDecayLenV0   = 0.2; //Decay length btw cascade to V0
		const Float_t cutCasc_minDcaBachToPV  = 0.01; //DCA of bachelor track to primary vertex
		const Float_t cutCasc_minDcaV0ToPV    = 0.01; //DCA of V0 to PV
		const Float_t cutCasc_maxDcaXiDau     = 1.68; //DCA of Cascade (Xi) to its daughters
		const Float_t cutCasc_maxDcaV0Dau     = 1.68; //DCA of V0 to its daughters
		const Float_t cutCasc_minDcaV0PosToPV = 0.05; //DCA of V0 daughter (positive) to PV
		const Float_t cutCasc_minDcaV0NegToPV = 0.05;
		const Float_t cutCasc_minCosPAngleXi  = 0.98;
		const Float_t cutCasc_minCosPAngleV0  = 0.98;

		//-----------------------------------------------------------

		//Tree variables for events
		UInt_t   fEvtID = 0;
		UInt_t   fEvtTrig;
		Int_t    fEvtRunNo;
		Float_t  fEvtMult;
		Double_t fEvtVtxZ;
		Bool_t   fEvtINELLgt0;
		Bool_t   fEvtPileupMB;
		Bool_t   fEvtPileupHMV0;

		Int_t    fMCNum;
		Int_t    fMCOrig  [MaxNXic0];
		Double_t fMCXic0Pt[MaxNXic0];
		Double_t fMCXic0Y [MaxNXic0];
		Double_t fMCCascPt[MaxNXic0];
		Double_t fMCCascY [MaxNXic0];
		Double_t fMCElePt [MaxNXic0];
		Double_t fMCEleY  [MaxNXic0];

		//Tree variables for tracks
		Int_t    fEleNum;
		Int_t    fEleChg      [MaxNEle];
		Int_t    fEleITSNcls  [MaxNEle]; //Previous notation: ITS
		Float_t  fEleMinMassLS[MaxNEle]; //Minimum mass of e+e- suspect from photon conversion, likesign
		Float_t  fEleMinMassUS[MaxNEle]; //Minimum mass of e+e- suspect from photon conversion, unlikesign
		Float_t  fEleNSigmaTOF[MaxNEle];
		Float_t  fEleNSigmaTPC[MaxNEle];
		Double_t fEleEta      [MaxNEle];
		Double_t fElePhi      [MaxNEle];
		Double_t fElePt       [MaxNEle];
		Double_t fElePx       [MaxNEle];
		Double_t fElePy       [MaxNEle];
		Double_t fElePz       [MaxNEle];
		Double_t fEleY        [MaxNEle];
		UShort_t fEleTPCNsig  [MaxNEle]; //Previous notation: TPCPID
		UShort_t fEleTPCNxedR [MaxNEle]; //Previous notation: e_crossedrows
		UShort_t fEleTPCNclsF [MaxNEle]; //Previous notation: e_findable

		Int_t    fCascNum;
		Int_t    fCascChgXi      [MaxNCasc];
		Double_t fCascCosPAXi    [MaxNCasc]; //Cosine of pointing angle
		Double_t fCascCosPAV0    [MaxNCasc];
		Double_t fCascDcaBachToPV[MaxNCasc]; //DCA of Bachelor track to Primary Vertex
		Double_t fCascDcaV0ToPV  [MaxNCasc];
		Double_t fCascDcaXiDau   [MaxNCasc]; //DCA of Xi daughters
		Double_t fCascDcaV0Dau   [MaxNCasc];
		Double_t fCascDcaPosToPV [MaxNCasc]; //DCA of Positive V0 daughter to PV
		Double_t fCascDcaNegToPV [MaxNCasc];
		//
		Double_t fCascDecayLenXi   [MaxNCasc]; //Decay length, Xi to PV
		Double_t fCascDecayLenXiOld[MaxNCasc]; //Decay length (in truth, radial length) at the old code
		Double_t fCascDecayLenV0   [MaxNCasc]; //Decay length, V0 ti Xi
		Double_t fCascDecayLenV0Old[MaxNCasc]; //Decay length (in truth, radial length) at the old code
		Double_t fCascMassLmb      [MaxNCasc]; //Lambda0
		Double_t fCascMassLmbAnti  [MaxNCasc]; //Lambda0_bar
		Double_t fCascMassOmega    [MaxNCasc];
		Double_t fCascMassXi       [MaxNCasc];
		Double_t fCascPtXi         [MaxNCasc]; //pT of Xi
		Double_t fCascPxXi         [MaxNCasc];
		Double_t fCascPyXi         [MaxNCasc];
		Double_t fCascPzXi         [MaxNCasc];
		//
		Double_t fCascPt_BachPi      [MaxNCasc];
		Double_t fCascPt_V0dPos      [MaxNCasc];
		Double_t fCascPt_V0dNeg      [MaxNCasc];
		UShort_t fCascTPCNxedR_BachPi[MaxNCasc]; //TPCNcrossedRows, previously bpion_crossedrows or crossedratio
		UShort_t fCascTPCNxedR_V0dPos[MaxNCasc]; //V0 daughter, positive
		UShort_t fCascTPCNxedR_V0dNeg[MaxNCasc]; //V0 daughter, negative
		UShort_t fCascTPCNclsF_BachPi[MaxNCasc]; //Previously bpion_findable
		UShort_t fCascTPCNclsF_V0dPos[MaxNCasc];
		UShort_t fCascTPCNclsF_V0dNeg[MaxNCasc];

		ClassDef(AliAnaTaskSEXic0SL, 1);
};

#endif //AliAnaTaskSEXic0SL_H
