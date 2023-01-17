#ifndef AliAnalysisTaskSEXic0SL_H
#define AliAnalysisTaskSEXic0SL_H

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

class AliAnalysisTaskSEXic0SL : public AliAnalysisTaskSE
{
	public:

		AliAnalysisTaskSEXic0SL():AliAnalysisTaskSE("AliAnalysisTaskSEXic0SL") {}
		AliAnalysisTaskSEXic0SL(const char* name, const char* option);
		virtual ~AliAnalysisTaskSEXic0SL();

		virtual void Terminate(Option_t* option);
		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t* option);

		//Options invoked via AddTask
		void SetMC(bool set) { IsMC = set; };
		void SetPA(bool set) { IsPA = set; };
		void SetTrigMB(bool set) { TrigOnMB = set; }
		void SetTrigHMV0(bool set) { TrigOnHMV0 = set; }
		void SetCutsLegacy(bool set) { IsLegacy = set; }
		void SetCutsByFile(bool set) { IsCutsByFile = set; }
		void SetValidEvtOnly(bool set) { ValidEvtOnly = set; }

		//User defined functions
		int CheckOrigin(AliMCEvent* MCEvt, AliAODMCParticle *MCPart); //<0:no_quark, 4:c, 5:b
		int GetCascLabel(AliMCEvent* MCEvt, AliAODcascade* Casc, bool getLabelXic0); //true for Xic0, false for Xi

		bool FilterTrack(AliAODTrack* Trk, const AliVVertex* Vtx);
		bool FilterTrackElectron(AliAODTrack* Trk, AliPIDResponse* PID);
		bool FilterCascade(AliAODcascade* Casc, const AliVVertex* Vtx, AliPIDResponse* PID);

		void ControlAnaObjects(int option); //0: set, 1: delete
		void ControlOutputContainers(int option); //0: define, 1: post
		void ControlOutputTree(TTree* T, bool isMC, bool readOnly = false);
		void ResetTreeVariables(void);
		void SetConstants(void);

		void ExecuteOffSel(TTree* T, TString Cond); //For offline selection of Grid/LT output

	private:

		TString fTaskOpt;

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
		bool IsMC;         //= false;
		bool IsPA;         //= false;
		bool TrigOnMB;     //= false;
		bool TrigOnHMV0;   //= false;
		bool IsLegacy;     //= false;
		bool IsCutsByFile; //= false;
		bool ValidEvtOnly; //= false;
		std::vector<UInt_t> TrigStore; //Container for triggers of interest

		//-----------------------------------------------------------

		//Constants
		int MaxNTruth;      //= 20; //max. # of generated particle w/ eXi pair per event (MC truth)
		int MaxNEle;        //= 20; //max. # of electrons per event, for both MC truth and reco
		int MaxNCasc;       //= 20; //max. # of Xi per event, for both MC truth and reco
		int PDGCode_e;      //= 11;
		int PDGCode_Lambda; //= 3122; //Lambda0, 1115.683 +- 0.006 (MeV)
		int PDGCode_Omega;  //= 3334; //Omega-, 1672.45 +- 0.29 (MeV)
		int PDGCode_Xi;     //= 3312; //Xi- (strange baryon), 1321.71 +- 0.07 (MeV)
		int PDGCode_Xistm;  //= 3314; //Xi*- (Xi 1530), 1535.0 +- 0.6 (MeV)
		int PDGCode_Xist0;  //= 3324; //Xi*0 (Xi 1530), 1531.80 +- 0.32 (MeV)
		int PDGCode_Xic0;   //= 4132; //2470.87 +0.28 -0.31 (MeV)
		int PDGCode_Xicp;   //= 4232; //Xic+, 2467.93 +- 0.18 (MeV)
		double MassEle;     //= 0.51099895 * 1.E-3; //Electron mass in GeV
		double MassLmb;     //= 1115.683 * 1.E-3; //Lambda0 mass in GeV
		double MassXi;      //= 1321.71 * 1.E-3; //Xi mass in GeV

		//Eventwise cut
		Int_t   cut_runNoLo;          //= 252000; 
		Int_t   cut_runNoUp;          //= 295000;
		Int_t   cut_vtxNContributors; //= 1;
		Float_t cut_bfield;           //= 0.001;
		Float_t cut_eta;              //= 1.0;
		Float_t cut_vtxZ;             //= 10.0;

		//Trackwise cut, common
		Int_t   cut_minNClustersITS;      //= 2;
		Int_t   cut_TPCsignalN;           //= 50; //fSetProdTrackTPCNclsPID in old code
		Float_t cut_maxChi2PerClusterITS; //= 36.;
		Float_t cut_maxDCAToVertexXY;     //= 1.0;
		Float_t cut_maxDCAToVertexZ;      //= 2.0;
		Float_t cut_trkEta;               //= 0.8; //For daughter particles
		Float_t cut_trkPt;                //= 0.5; //Lower limit of electron

		//Trackwise cut, electron
		Float_t cutEle_massConv;         //= 0.05; //GeV, max. conversion mass
		Float_t cutEle_nSigmaTOFAbs;     //= 3.0;
		Float_t cutEle_nSigmaTPCAbsConv; //= 5.0; //Used only to check conversion mass (FilterTrackElectron)
		Float_t cutEle_nSigmaTPCMax;     //= 3.0; //Cf. Min: pT dependent (FilterTrackElectron)

		//Trackwise cut, cascade
		Float_t cutCasc_massTolLambda;   //= 0.008;
		Float_t cutCasc_massTolOmega;    //= 0.3 * 1.E-3; //ckim, NOT used for now
		Float_t cutCasc_massTolXi;       //= 0.03; //!! Old code: online 0.01, offline selection 0.008
		Float_t cutCasc_nSigmaTPCAbs;    //= 4.0;
		Float_t cutCasc_minDecayLenXi;   //= 0.2; //Decay length btw PV to cascade
		Float_t cutCasc_minDecayLenV0;   //= 0.2; //Decay length btw cascade to V0
		Float_t cutCasc_minDcaBachToPV;  //= 0.01; //DCA of bachelor track to primary vertex
		Float_t cutCasc_minDcaV0ToPV;    //= 0.01; //DCA of V0 to PV
		Float_t cutCasc_maxDcaXiDau;     //= 1.68; //DCA of Cascade (Xi) to its daughters
		Float_t cutCasc_maxDcaV0Dau;     //= 1.68; //DCA of V0 to its daughters
		Float_t cutCasc_minDcaV0PosToPV; //= 0.05; //DCA of V0 daughter (positive) to PV
		Float_t cutCasc_minDcaV0NegToPV; //= 0.05;
		Float_t cutCasc_minCosPAngleXi;  //= 0.98;
		Float_t cutCasc_minCosPAngleV0;  //= 0.98;

		//-----------------------------------------------------------

		//Tree variables for events
		UInt_t   fEvtID;
		UInt_t   fEvtTrig;
		Int_t    fEvtRunNo;
		Float_t  fEvtMult;
		Double_t fEvtVtxZ;
		Bool_t   fEvtGoodMB;
		Bool_t   fEvtGoodHMV0;
		Bool_t   fEvtINELLgt0;

		Int_t    fMCNum;
		Int_t    fMCLabel[20]; //20: MaxNTruth
		Int_t    fMCNDau [20];
		Int_t    fMCOrig [20];
		Int_t    fMCPDG  [20];
		Double_t fMCPt   [20];
		Double_t fMCY    [20];
		Double_t fMCElePt[20];
		Double_t fMCEleY [20];
		Double_t fMCXiPt [20];
		Double_t fMCXiY  [20];

		//Tree variables for tracks
		Int_t    fEleNum;
		Int_t    fEleChg      [20]; //20: MaxNEle
		Int_t    fEleITSNcls  [20]; //Previous notation: ITS
		Float_t  fEleMinMassLS[20]; //Minimum mass of e+e- suspect from photon conversion, likesign
		Float_t  fEleMinMassUS[20]; //Minimum mass of e+e- suspect from photon conversion, unlikesign
		Float_t  fEleNSigmaTOF[20];
		Float_t  fEleNSigmaTPC[20];
		Double_t fEleEta      [20];
		Double_t fElePhi      [20];
		Double_t fElePt       [20];
		Double_t fElePx       [20];
		Double_t fElePy       [20];
		Double_t fElePz       [20];
		Double_t fEleY        [20];
		UShort_t fEleTPCNsig  [20]; //Previous notation: TPCPID
		UShort_t fEleTPCNxedR [20]; //Previous notation: e_crossedrows
		UShort_t fEleTPCNclsF [20]; //Previous notation: e_findable
		//
		Int_t fEleLabel   [20]; //MC only, electron candidate's label, to check if it's negative or not
		Int_t fElePDG     [20]; //MC only, electron candidate's PDG code
		Int_t fEleMomLabel[20]; //MC only, mother particle's (Xic0) label
		Int_t fEleMomPDG  [20]; //MC only, mother particle's (Xic0) PDG code

		//Tree variables for cascades
		Int_t    fCascNum;
		Int_t    fCascChgXi      [20]; //20: MaxNCasc
		Double_t fCascCosPAXi    [20]; //Cosine of pointing angle
		Double_t fCascCosPAV0    [20];
		Double_t fCascDcaBachToPV[20]; //DCA of Bachelor track to Primary Vertex
		Double_t fCascDcaV0ToPV  [20];
		Double_t fCascDcaXiDau   [20]; //DCA of Xi daughters
		Double_t fCascDcaV0Dau   [20];
		Double_t fCascDcaPosToPV [20]; //DCA of Positive V0 daughter to PV
		Double_t fCascDcaNegToPV [20];
		//
		Double_t fCascDecayLenXi   [20]; //Decay length, Xi to PV
		Double_t fCascDecayLenXiOld[20]; //Decay length (in truth, radial length) at the old code
		Double_t fCascDecayLenV0   [20]; //Decay length, V0 ti Xi
		Double_t fCascDecayLenV0Old[20]; //Decay length (in truth, radial length) at the old code
		Double_t fCascMassLmb      [20]; //Lambda0
		Double_t fCascMassLmbAnti  [20]; //Lambda0_bar
		Double_t fCascMassOmega    [20];
		Double_t fCascMassXi       [20];
		Double_t fCascPtXi         [20]; //pT of Xi
		Double_t fCascPxXi         [20];
		Double_t fCascPyXi         [20];
		Double_t fCascPzXi         [20];
		//
		Double_t fCascPt_BachPi      [20];
		Double_t fCascPt_V0dPos      [20];
		Double_t fCascPt_V0dNeg      [20];
		UShort_t fCascTPCNxedR_BachPi[20]; //TPCNcrossedRows, previously bpion_crossedrows or crossedratio
		UShort_t fCascTPCNxedR_V0dPos[20]; //V0 daughter, positive
		UShort_t fCascTPCNxedR_V0dNeg[20]; //V0 daughter, negative
		UShort_t fCascTPCNclsF_BachPi[20]; //Previously bpion_findable
		UShort_t fCascTPCNclsF_V0dPos[20];
		UShort_t fCascTPCNclsF_V0dNeg[20];
		//
		Int_t fCascPDG     [20]; //MC only, Xi candidate's PDG code
		Int_t fCascMomLabel[20]; //MC only, mother particle's (Xic0) label
		Int_t fCascMomPDG  [20]; //MC only, mother particle's (Xic0) PDG code

		ClassDef(AliAnalysisTaskSEXic0SL, 1);
};

#endif //AliAnalysisTaskSEXic0SL_H
