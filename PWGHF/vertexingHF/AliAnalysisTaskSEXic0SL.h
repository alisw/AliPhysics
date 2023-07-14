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

		AliAnalysisTaskSEXic0SL();
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

		bool FilterTrack(AliAODTrack* Trk, const AliVVertex* Vtx, float cutPt);
		bool FilterTrackElectron(AliAODTrack* Trk, AliPIDResponse* PID);
		bool FilterCascade(AliAODcascade* Casc, const AliVVertex* Vtx, AliPIDResponse* PID);

		void ControlAnaObjects(int option); //0: set, 1: delete
		void ControlOutputContainers(int option); //0: define, 1: post
		void ControlOutputTree(TTree* T, bool isMC, bool readOnly = false);
		void DeleteTreeVariables(void);
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
		Float_t cut_maxChi2PerClusterTPC; //= 4; //Updated July 6, 2023
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
		Float_t  fEvtNSPDtl; //NSPDTracklets, July 14 (2023)
		Double_t fEvtVtxZ;
		Bool_t   fEvtGoodMB;
		Bool_t   fEvtGoodHMV0;
		Bool_t   fEvtINELLgt0;

		//Tree variables for MC truth
		Int_t     fMCNum;
		Int_t*    fMCLabel; //[fMCNum]
		Int_t*    fMCOrig;  //[fMCNum]
		Int_t*    fMCPDG;   //[fMCNum]
		Double_t* fMCPt;    //[fMCNum]
		Double_t* fMCY;     //[fMCNum]
		Double_t* fMCElePt; //[fMCNum]
		Double_t* fMCEleY;  //[fMCNum]
		Double_t* fMCXiPt;  //[fMCNum]
		Double_t* fMCXiY;   //[fMCNum]
		Int_t*    fMCXiMomLabel; //[fMCNum]
		Int_t*    fMCXiMomPDG;   //[fMCNum]

		//Tree variables for tracks
		Int_t     fEleNum;
		Int_t*    fEleChg;       //[fEleNum]
		Int_t*    fEleITSNcls;   //[fEleNum]
		Float_t*  fEleMinMassLS; //[fEleNum]
		Float_t*  fEleMinMassUS; //[fEleNum]
		Float_t*  fEleNSigmaTOF; //[fEleNum]
		Float_t*  fEleNSigmaTPC; //[fEleNum]
		Double_t* fEleDCAd;      //[fEleNum]
		Double_t* fEleDCAz;      //[fEleNum]
		Double_t* fEleEta;       //[fEleNum]
		Double_t* fElePhi;       //[fEleNum]
		Double_t* fElePt;        //[fEleNum]
		Double_t* fElePx;        //[fEleNum]
		Double_t* fElePy;        //[fEleNum]
		Double_t* fElePz;        //[fEleNum]
		Double_t* fEleY;         //[fEleNum]
		UShort_t* fEleTPCNsig;   //[fEleNum]
		UShort_t* fEleTPCNxedR;  //[fEleNum]
		UShort_t* fEleTPCNclsF;  //[fEleNum]
		///Tree variables for tracks, MC only heareafter
		Int_t* fEleLabel;    //[fEleNum]
		Int_t* fElePDG;      //[fEleNum]
		Int_t* fEleMomLabel; //[fEleNum]
		Int_t* fEleMomPDG;   //[fEleNum]

		//Tree variables for cascades
		Int_t     fCascNum;
		Int_t*    fCascChgXi;       //[fCascNum]
		Double_t* fCascCosPAXi;     //[fCascNum]
		Double_t* fCascCosPAV0;     //[fCascNum]
		Double_t* fCascDcaBachToPV; //[fCascNum]
		Double_t* fCascDcaV0ToPV;   //[fCascNum]
		Double_t* fCascDcaXiDau;    //[fCascNum]
		Double_t* fCascDcaV0Dau;    //[fCascNum]
		Double_t* fCascDcaPosToPV;  //[fCascNum]
		Double_t* fCascDcaNegToPV;  //[fCascNum]
		//
		Double_t* fCascDecayLenXi;    //[fCascNum]
		Double_t* fCascDecayLenXiOld; //[fCascNum]
		Double_t* fCascDecayLenV0;    //[fCascNum]
		Double_t* fCascDecayLenV0Old; //[fCascNum]
		Double_t* fCascMassLmb;       //[fCascNum]
		Double_t* fCascMassLmbAnti;   //[fCascNum]
		Double_t* fCascMassOmega;     //[fCascNum]
		Double_t* fCascMassXi;        //[fCascNum]
		Double_t* fCascMassXi1530;    //[fCascNum]
		Double_t* fCascPtXi;          //[fCascNum]
		Double_t* fCascPxXi;          //[fCascNum]
		Double_t* fCascPyXi;          //[fCascNum]
		Double_t* fCascPzXi;          //[fCascNum]
		//
		Double_t* fCascPt_BachPi;       //[fCascNum]
		Double_t* fCascPt_V0dPos;       //[fCascNum]
		Double_t* fCascPt_V0dNeg;       //[fCascNum]
		UShort_t* fCascTPCNxedR_BachPi; //[fCascNum]
		UShort_t* fCascTPCNxedR_V0dPos; //[fCascNum]
		UShort_t* fCascTPCNxedR_V0dNeg; //[fCascNum]
		UShort_t* fCascTPCNclsF_BachPi; //[fCascNum]
		UShort_t* fCascTPCNclsF_V0dPos; //[fCascNum]
		UShort_t* fCascTPCNclsF_V0dNeg; //[fCascNum]
		//Tree variables for cascades, MC only hereafter
		Int_t* fCascPDG;      //[fCascNum]
		Int_t* fCascMomLabel; //[fCascNum]
		Int_t* fCascMomPDG;   //[fCascNum]

		ClassDef(AliAnalysisTaskSEXic0SL, 4);
};

#endif //AliAnalysisTaskSEXic0SL_H
