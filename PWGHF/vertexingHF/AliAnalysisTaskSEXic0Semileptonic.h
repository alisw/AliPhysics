/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#ifndef AliAnalysisTaskSEXic0Semileptonic_H
#define AliAnalysisTaskSEXic0Semileptonic_H

#include "TROOT.h"
#include "TVector.h"
#include "TVector2.h"
#include "TSystem.h"
#include "TProfile.h"
#include "THistManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliRDHFCuts.h"
#include "AliAODRecoDecayHF.h"
#include "AliNormalizationCounter.h"
#include "AliRDHFCutsXictoeleXifromAODtracks.h"
#include <vector>

//kimc
#include <iostream>
using namespace std;

class AliNormalizationCounter;

class AliAnalysisTaskSEXic0RunTable
{
	public:

		enum {kPP, kPA, kAA, kUnknownCollType};

		AliAnalysisTaskSEXic0RunTable();
		AliAnalysisTaskSEXic0RunTable(Int_t runnumber);
		~AliAnalysisTaskSEXic0RunTable();

		//Bool_t IsPP() { return fCollisionType==kPP; };
		//Bool_t IsPA() { return fCollisionType==kPA; };
		//Bool_t IsAA() { return fCollisionType==kAA; };

	private:

		Int_t fCollisionType; //! Is proton-proton collisions?
};

class AliAnalysisTaskSEXic0Semileptonic : public AliAnalysisTaskSE
{
	public:

		enum {kSD=0, kDD, kND, kCD, kAllProc}; //PN = unlike sign, PP and NN are like signs

		AliAnalysisTaskSEXic0Semileptonic();
		AliAnalysisTaskSEXic0Semileptonic( const char *name, const char *option );
		AliAnalysisTaskSEXic0Semileptonic( const AliAnalysisTaskSEXic0Semileptonic& ap );
		AliAnalysisTaskSEXic0Semileptonic& operator = ( const AliAnalysisTaskSEXic0Semileptonic& ap );
		~AliAnalysisTaskSEXic0Semileptonic();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *);
		virtual void FinishTaskOutput();
		virtual void Terminate(Option_t *);

		void SetMC(Bool_t ismc) { IsMC = ismc; };
		void SetPA(Bool_t ispa) { IsPA = ispa; };
		void SetRunOffset(Int_t RunOffset) { fRunOffset = RunOffset; };
		void SetHighMultiplicity(Bool_t IsHM) { IsHighMul = IsHM; };
		Int_t GetRunOffset() { return fRunOffset; };

		void DefineMCXicTree();
		void DefinePaireXiTree();
		void DefineMCPaireXiTree();
		void DefineEventTree();

		Bool_t FilterTrack(AliAODTrack *trk, Int_t NoFillHistos);
		Bool_t FilterElectron(AliAODTrack *trk, Double_t &mass, Double_t &samesign_mass,
				              Bool_t IsSameSign, Bool_t IsLoose, Int_t NoFillHistos);
		Bool_t FilterCascade(AliAODcascade *casc);
		Bool_t FillMCXib(AliAODMCParticle *mcpart);

		void FillBottomContribution(AliAODMCParticle *mcpart, AliAODcascade *casc, AliAODTrack *trk);
		void FillPairEleXi(AliAODcascade *casc, AliAODTrack *trk);
		void FillMCXic0(AliAODMCParticle *mcpart);
		void FillXiHistFromPromptNonPrompt(Bool_t ismc, AliAODcascade *casc);

		Int_t MatchToMCXic0(AliAODcascade *casc, AliAODTrack *trk);
		Int_t MatchToMCXib(AliAODcascade *casc, AliAODTrack *trk);
		Int_t MatchToMCXic0(AliAODcascade *casc);
		Int_t MatchToMCXi(AliAODcascade *casc);
		Int_t MatchToMCele(AliAODTrack *trk);

		//Type: 1=loose, 2=standard, 3=tight
		Bool_t StandardCutFlag(AliAODTrack *track, AliAODcascade *casc,
				               Bool_t e_reco, Bool_t e_pid, Bool_t Xi_reco, Bool_t Xi_pid);

		unsigned int GetEvID();

		//kimc
		//*************************************************

		//Add high multiplicity triggers in event selection
		void UseTrig_kINT7(void) { fUsekINT7 = kTRUE; }
		void UseTrig_kHMV0(void) { fUsekHMV0 = kTRUE; }
		void UseTrig_kHMSPD(void) { fUsekHMSPD = kTRUE; }

	private:

		AliAnalysisTaskSEXic0RunTable      *fRunTable = nullptr; //!
		AliESDtrackCuts                    *fTrackCuts = nullptr; //!
		AliMCEvent                         *fMC = nullptr; //!
		AliNormalizationCounter            *fCounter = nullptr; //!<! Counter for normalization
		AliPIDResponse                     *fPIDResponse = nullptr; //!
		AliRDHFCutsXictoeleXifromAODtracks *fEvtCuts = nullptr;
		AliVEvent                          *fEvt = nullptr; //!

		TString fOption; //!
		Bool_t  IsMC = kTRUE;
		Bool_t  IsHighMul = kFALSE; //high multiplicity condition
		Bool_t  IsPA = kFALSE; // to use V0A for pPb

		TList        *fOutput = nullptr; //!
		THistManager *fHistos = nullptr; //!
		TTree        *fMCXicTree = nullptr; //!
		TTree        *fPaireXiTree = nullptr; //!
		TTree        *fMCTree = nullptr; //!
		TTree        *fEventTree = nullptr; //!

		Float_t* fMCXicTreeVariable = nullptr; //!
		Float_t* fPaireXiTreeVariable = nullptr; //!
		Float_t* fMCTreeVariable = nullptr; //!
		Float_t* fEventTreeVariable = nullptr; //!

		//Cut Values
		//---------------------------------------------

		Int_t   fRunNumber = 0; //!
		Int_t   fNSPDTracklets = 0; //kimc
		Int_t   fNeXiPair      = 0; //kimc
		Float_t fBzkG = 0; //! Magnetic filed for of event
		Float_t fCentrality = 9999; //!
		Float_t fEtaCut = 0.8; //for daugther particle
		Float_t fPtCut = 0.5; //lower limit of electron
		Float_t fRunOffset = 0; //!
		Float_t fCentralSPD = 9999; //kimc
		Float_t fVtxZ       = 9999; //kimc

		Double_t MassTolLambda = 0.008;
		Double_t MassTolXi = 0.01;
		Double_t DCAV0PrToPrimVertexMin = 0.05;
		Double_t DCAV0PiToPrimVertexMin = 0.05;
		Double_t DCABachToPrimVertexMin = 0.01;
		Double_t DCAV0ToPrimVertexMin = 0.01;
		Double_t V0CosineOfPoiningAngleXiMin = 0.98;
		Double_t CascDecayLengthMin = 0.2;
		Double_t DecayLengthV0 = 0.2;

		Double_t fSetProdTrackTPCNclsPID = 50;
		Double_t fNClustersITSMin = 2;
		Double_t fNXiCrossedRowsMin = 70;
		Double_t fNXiCrossedRowsOverFindalbeRatioMin = 0.77;
		//Double_t fProdTrackTPCNclsRatioMin = 0.6;
		//Double_t fNClustersTPCMin = 70;

		//kimc
		//*************************************************

		UInt_t fEventTreeVarTrig = 0; //Dummy index for saving triggerbit to tree

        vector<UInt_t> fTargetTriggers; //Container for triggerbits
        Bool_t fUsekINT7 = kFALSE;
        Bool_t fUsekHMV0 = kFALSE;
        Bool_t fUsekHMSPD = kFALSE;

		ClassDef(AliAnalysisTaskSEXic0Semileptonic, 1);
};

#endif //AliAnalysisTaskSEXic0Semileptonic
