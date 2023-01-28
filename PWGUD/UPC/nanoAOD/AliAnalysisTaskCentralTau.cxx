/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

// c++ headers
#include <iostream>

// root headers
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1I.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDZDC.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"
#include "AliTimeRangeCut.h"
#include "AliTOFTriggerMask.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"

// my headers
#include "AliAnalysisTaskCentralTau.h"

ClassImp(AliAnalysisTaskCentralTau);

using std::cout;
using std::endl;

//analysis skeleton of UPC nano AODs,

//_____________________________________________________________________________
AliAnalysisTaskCentralTau::AliAnalysisTaskCentralTau()
		: AliAnalysisTaskSE(),
		  fPIDResponse(nullptr), fTrackCutsBit0(nullptr), fTrackCutsBit1(nullptr), fTrackCutsBit4(nullptr), fOutputList(nullptr),
		  fOutputPID(nullptr), tTwoTracks(nullptr), tPID(nullptr), hTriggerCounter(nullptr), hParticleTypeCounter(nullptr), fESDtracks(nullptr),
			fSign(0), fZNAenergy(0), fZNCenergy(0), fRunNumber(0), fADAdecision(0), fADCdecision(0), fV0Adecision(0), fV0Cdecision(0)
{
//Dummy constructor

}//AliAnalysisTaskCentralTau


//_____________________________________________________________________________
AliAnalysisTaskCentralTau::AliAnalysisTaskCentralTau(const char *name)
		: AliAnalysisTaskSE(name),
		  fPIDResponse(nullptr), fTrackCutsBit0(nullptr), fTrackCutsBit1(nullptr), fTrackCutsBit4(nullptr), fOutputList(nullptr),
		  fOutputPID(nullptr), tTwoTracks(nullptr),tPID(nullptr), hTriggerCounter(nullptr), hParticleTypeCounter(nullptr), fESDtracks(nullptr),
			fSign(0), fZNAenergy(0), fZNCenergy(0), fRunNumber(0), fADAdecision(0), fADCdecision(0), fV0Adecision(0), fV0Cdecision(0)
{
	for(bool & fTrigger : fTriggers)          fTrigger = false;
	for(bool & fTriggerClas : fTriggerClass) fTriggerClas = false;
	for (int it(0);it<2;it++) {
		fTrackPIDid[it] = -1;
		fPIDpt[it] = -999;
		fPIDmomentum[it] = -999;
		fTPCsignal[it] = -999;
		fTOFsignal[it] = -999;
		fTPCmostProbableTrackType[it] = -1;
		fTOFmostProbableTrackType[it] = -1;
		for (int typ(0);typ<5;typ++){
			fPIDTPC[typ][it] = -999.;
			fPIDTOF[typ][it] = -999.;
		}
	}
	for (int it(0);it<4;it++){
		fZNAtime[it] = 0.;
		fZNCtime[it] = 0.;
	}
	DefineOutput(1, TList::Class());
	DefineOutput(2, TList::Class());

}//AliAnalysisTaskCentralTau

//_____________________________________________________________________________
AliAnalysisTaskCentralTau::~AliAnalysisTaskCentralTau()
{
	// Destructor

	// Destructor
	if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
		delete fOutputList;
		fOutputList = nullptr;
		delete fOutputPID;
		fOutputPID = nullptr;
		delete fESDtracks;
		fESDtracks = nullptr;
	}

}//~AliAnalysisTaskCentralTau


//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::UserCreateOutputObjects()
{

	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	auto *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
	fPIDResponse = inputHandler->GetPIDResponse();

	fTrackCutsBit0 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	fTrackCutsBit0->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
	fTrackCutsBit1 = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kFALSE,kTRUE);
	fTrackCutsBit4 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);

	fESDtracks = new TClonesArray("AliESDtrack", 2);

	fOutputList = new TList();
	fOutputList ->SetOwner();

	tTwoTracks = new TTree("tTwoTracks", "tTwoTracks");
	tTwoTracks ->Branch("fESDtracks", &fESDtracks);
	tTwoTracks ->Branch("fPIDTPC", &fPIDTPC[0][0],"fPIDTPC[5][2]/F");
	tTwoTracks ->Branch("fPIDTOF", &fPIDTOF[0][0],"fPIDTOF[5][2]/F");
	tTwoTracks ->Branch("fZNAenergy", &fZNAenergy,"fZNAenergy/D");
	tTwoTracks ->Branch("fZNCenergy", &fZNCenergy,"fZNCenergy/D");
	tTwoTracks ->Branch("fZNAtime", &fZNAtime[0],"fZNAtime[4]/D");
	tTwoTracks ->Branch("fZNCtime", &fZNCtime[0],"fZNCtime[4]/D");
	tTwoTracks ->Branch("fSign", &fSign, "fSign/I");
	tTwoTracks ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
	tTwoTracks ->Branch("fTriggers", &fTriggers, Form("fTriggers[%i]/O",(Int_t)(sizeof(fTriggers)/sizeof(fTriggers[0]))));
	tTwoTracks ->Branch("fADAdecision", &fADAdecision, "fADAdecision/I");
	tTwoTracks ->Branch("fADCdecision", &fADCdecision, "fADCdecision/I");
	tTwoTracks ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
	tTwoTracks ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");
	tTwoTracks ->Branch("fPIDpt", &fPIDpt[0], "fPIDpt[2]/D");
	tTwoTracks ->Branch("fPIDmomentum", &fPIDmomentum[0], "fPIDmomentum[2]/D");
	tTwoTracks ->Branch("fTPCsignal", &fTPCsignal[0], "fTPCsignal[2]/D");
	tTwoTracks ->Branch("fTOFsignal", &fTOFsignal[0], "fTOFsignal[2]/D");
	tTwoTracks ->Branch("fTPCmostProbableTrackType", &fTPCmostProbableTrackType[0], "fTPCmostProbableTrackType[2]/I");
	tTwoTracks ->Branch("fTOFmostProbableTrackType", &fTOFmostProbableTrackType[0], "fTOFmostProbableTrackType[2]/I");
	fOutputList->Add(tTwoTracks);

	hTriggerCounter = new TH2I("hTriggerCounter","Number of analyzed UPC triggers per run",3,1,4,3000,295000,298000);
	fOutputList->Add(hTriggerCounter);
	hParticleTypeCounter = new TH1I("hParticleTypeCounter","Electron, Muon, Pion, Kaon, Proton",6,-0.5,5.5);
	fOutputList->Add(hParticleTypeCounter);

	fOutputPID = new TList();
	fOutputPID ->SetOwner(); // @suppress("Ambiguous problem")

	tPID = new TTree("tPID", "tPID");
	tPID ->Branch("fPIDpt", &fPIDpt[0], "fPIDpt[2]/D");
	tPID ->Branch("fPIDmomentum", &fPIDmomentum[0], "fPIDmomentum[2]/D");
	tPID ->Branch("fTPCsignal", &fTPCsignal[0], "fTPCsignal[2]/D");
	tPID ->Branch("fTOFsignal", &fTOFsignal[0], "fTOFsignal[2]/D");
	tPID ->Branch("fTPCmostProbableTrackType", &fTPCmostProbableTrackType[0], "fTPCmostProbableTrackType[2]/I");
	tPID ->Branch("fTOFmostProbableTrackType", &fTOFmostProbableTrackType[0], "fTOFmostProbableTrackType[2]/I");
	fOutputPID->Add(tPID);


	PostData(1, fOutputList);
	PostData(2, fOutputPID);

}//UserCreateOutputObjects


//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::UserExec(Option_t *)
{
	enum MyParticle { P_ELECTRON = 0, P_MUON = 1, P_PION = 2, P_KAON = 3, P_PROTON = 4};
	//
	// META
	//
	AliVEvent *fEvent = InputEvent();
	if(!fEvent) return;

	fRunNumber = fEvent->GetRunNumber();

	AliTimeRangeCut fTimeRangeCut;
	fTimeRangeCut.InitFromEvent(InputEvent());
	if(fTimeRangeCut.CutEvent(InputEvent()))return;

	//
	// TRIGGER
	//
	TString trigger = fEvent->GetFiredTriggerClasses();

	if(fRunNumber>=295881 && fRunNumber<296594)	if(!trigger.Contains("CCUP29-B-SPD2-CENTNOTRD") && !trigger.Contains("CCUP30-B-SPD2-CENTNOTRD") && !trigger.Contains("CCUP31-B-SPD2-CENTNOTRD"))return;
	if(fRunNumber>=296594) 					          	if(!trigger.Contains("CCUP29-U-SPD2-CENTNOTRD") && !trigger.Contains("CCUP30-B-SPD2-CENTNOTRD") && !trigger.Contains("CCUP31-B-SPD2-CENTNOTRD"))return;
	if(fRunNumber< 295881)  						        if(!trigger.Contains("CCUP29-B-NOPF-CENTNOTRD") && !trigger.Contains("CCUP30-B-NOPF-CENTNOTRD") && !trigger.Contains("CCUP31-B-NOPF-CENTNOTRD"))return;


	UInt_t fL0inputs = fEvent->GetHeader()->GetL0TriggerInputs();
	fTriggers[0] = trigger.Contains("CCUP29-B-SPD2-CENTNOTRD");
	fTriggers[1] = trigger.Contains("CCUP29-B-NOPF-CENTNOTRD");
	fTriggers[2] = trigger.Contains("CCUP29-U-SPD2-CENTNOTRD");
	fTriggers[3] = trigger.Contains("CCUP30-B-NOPF-CENTNOTRD");
	fTriggers[4] = trigger.Contains("CCUP30-B-SPD2-CENTNOTRD");
	fTriggers[5] = trigger.Contains("CCUP31-B-NOPF-CENTNOTRD");
	fTriggers[6] = trigger.Contains("CCUP31-B-SPD2-CENTNOTRD");
	fTriggers[7] =  fL0inputs & (1 << 11);//OM2
	fTriggers[8] =  fL0inputs & (1 << 12);//OMU

	if(trigger.Contains("CCUP29-B") || trigger.Contains("CCUP29-U")) hTriggerCounter->Fill(1,fRunNumber);
	if(trigger.Contains("CCUP30-B"))                                    hTriggerCounter->Fill(2,fRunNumber);
	if(trigger.Contains("CCUP31-B"))                                    hTriggerCounter->Fill(3,fRunNumber);

	//
	// VERTEX CUT
	//
	const AliVVertex *fVertex = fEvent->GetPrimaryVertex();
	if(fVertex->GetNContributors()<2) return;
	if(TMath::Abs(fVertex->GetZ())>15)return;

	//
	// OFFLINE VETOS INFORMATION
	//
	AliVVZERO *fV0data = fEvent->GetVZEROData();
	AliVAD *fADdata = fEvent->GetADData();

	auto *fZDCdata = (AliESDZDC*)fEvent->GetZDCData();
	fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
	fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
	Int_t detChZNA  = fZDCdata->GetZNATDCChannel();
	Int_t detChZNC  = fZDCdata->GetZNCTDCChannel();
	if (fEvent->GetRunNumber()>=245726 && fEvent->GetRunNumber()<=245793) detChZNA = 10;
	for (Int_t i=0;i<4;i++){
		fZNAtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNA,i);
		fZNCtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNC,i);
	}

	fV0Adecision = fV0data->GetV0ADecision();
	fV0Cdecision = fV0data->GetV0CDecision();

	fADAdecision = fADdata->GetADADecision();
	fADCdecision = fADdata->GetADCDecision();

	//
	// TRACKS INFO
	//
	UInt_t nGoodTracksTPC=0;
	UInt_t nGoodTracksSPD=0;
	Int_t TrackIndexTPC[5] = {-1,-1,-1,-1,-1};

	//
	// LOOP OVER TRACKS
	//
	if (fEvent->GetNumberOfTracks() < 1) return;
	for(Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++){
		Bool_t goodTPCTrack = kTRUE;
		auto *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(iTrack));
		if( !trk ) continue;
		if(!fTrackCutsBit4->AcceptTrack(trk)) goodTPCTrack = kFALSE;
		else{ if(trk->HasPointOnITSLayer(0) && trk->HasPointOnITSLayer(1))  nGoodTracksSPD++;}
		if(goodTPCTrack){
			TrackIndexTPC[nGoodTracksTPC] = iTrack;
			nGoodTracksTPC++;
		}
		if (nGoodTracksTPC == 5) {
			Printf("***** WARNING: Event has 5 TPC good tracks. Loop over tracks terminated.");
			break;
		}
	}//Track loop

	//
	// TRIGGER INFO
	//
	Bool_t isCUP29, isCUP30, isCUP31;
	isCUP29 = fTriggers[0] || fTriggers[1] || fTriggers[2];
	isCUP30 = fTriggers[3] || fTriggers[4];
	isCUP31 = fTriggers[5] || fTriggers[6];
	fTriggerClass[0] = isCUP29; fTriggerClass[1] = isCUP30; fTriggerClass[2] = isCUP31;

	//
	// SETUP STG
	//
	Int_t crossedFO[4];
	TBits fFOCrossedChips(1200);
	const AliVMultiplicity *mult = fEvent->GetMultiplicity();
	TBits fFOFiredChips = mult->GetFastOrFiredChips();

	//
	// TRACKS SETTING
	//
	Short_t qTrack[2];

	//
	// SELECT TWO TRACKS
	//
	fESDtracks->Clear("C");
	if(nGoodTracksTPC == 2 && nGoodTracksSPD == 2){
		fFOCrossedChips.ResetAllBits(kFALSE);

		//
		// LOOP OVER TWO TRACKS
		//
		for(Int_t iTrack=0; iTrack<2; iTrack++) {

			//
			// Fast-OR chips crossed STG
			//
			auto *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(TrackIndexTPC[iTrack]));
			if(trk->Pt()>10.) return; // just skip this event completelly
			crossedFO[0] = trk->GetITSModuleIndex(0);
			crossedFO[1] = trk->GetITSModuleIndex(1);
			crossedFO[2] = trk->GetITSModuleIndex(6);
			crossedFO[3] = trk->GetITSModuleIndex(7);
			SetCrossed(crossedFO, fFOCrossedChips);

			// PID standalone analysis
			TPCandTOFsignalInfo(trk, iTrack);

			//
			// PID AND KINEMATIC INFO
			//

			qTrack[iTrack] = trk->Charge();

			new((*fESDtracks)[iTrack]) AliESDtrack(*trk);

			fPIDTPC[P_ELECTRON][iTrack] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron));
			fPIDTPC[P_MUON][iTrack] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon));
			fPIDTPC[P_PION][iTrack] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion));
			fPIDTPC[P_KAON][iTrack] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon));
			fPIDTPC[P_PROTON][iTrack] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton));

			fPIDTOF[P_ELECTRON][iTrack] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron));
			fPIDTOF[P_MUON][iTrack] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon));
			fPIDTOF[P_PION][iTrack] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion));
			fPIDTOF[P_KAON][iTrack] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon));
			fPIDTOF[P_PROTON] [iTrack]= TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton));

			fTPCsignal[iTrack] = trk->GetTPCsignal();
			fTOFsignal[iTrack] = trk->GetTOFsignal();
			fTrackPIDid[iTrack] = TestPIDhypothesis(trk);
			fPIDpt[iTrack] = trk->Pt();
			fPIDmomentum[iTrack] = trk->P();
			hParticleTypeCounter->Fill(fTrackPIDid[iTrack]);

		}//Two tracks loop
		
		//
		// SAVE STG DECISION
		//
		TBits fFOCrossFiredChips = fFOCrossedChips & fFOFiredChips;
		fTriggers[9] = IsSTGFired(fFOCrossFiredChips,fRunNumber >= 295753 ? 9 : 3);

		//
		// CHARGE? (OPPOSITE, LIKE)
		//
		if(qTrack[0]*qTrack[1]<0)fSign = -1;
		if(qTrack[0]*qTrack[1]>0)fSign = 1;

		tPID->Fill();
		tTwoTracks->Fill();

	}//Two good tracks if

	PostData(1, fOutputList);
	PostData(2, fOutputPID);

}//UserExec

//_____________________________________________________________________________
Int_t AliAnalysisTaskCentralTau::TestPIDhypothesis(AliESDtrack *trk)
// Choose, which particle it is accroding to PID
{

	Float_t PIDTPC[5];
	PIDTPC[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron));
	PIDTPC[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon));
	PIDTPC[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion));
	PIDTPC[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon));
	PIDTPC[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton));
	Int_t idTPC = std::distance(std::begin(PIDTPC),std::min_element(std::begin(PIDTPC),std::end(PIDTPC)));
	Float_t PIDTPCpick[3] = {PIDTPC[0],PIDTPC[1],PIDTPC[2]};
	Int_t idTPCpick = std::distance(std::begin(PIDTPCpick),std::min_element(std::begin(PIDTPCpick),std::end(PIDTPCpick)));

	Float_t PIDTOF[5];
	PIDTOF[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron));
	PIDTOF[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon));
	PIDTOF[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion));
	PIDTOF[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon));
	PIDTOF[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton));
	Int_t idTOF = std::distance(std::begin(PIDTOF),std::min_element(std::begin(PIDTOF),std::end(PIDTOF)));

	if      (idTPC == 0 || idTPC == 1 || idTPC == 2){
		if      (idTOF == 3) return 3; // probably kaon
		else if (idTOF == 4) return 4; // probably proton
		else {
			if      (idTPC == 0) return 0; // probably electron
			else if (idTPC == 1) return 1; // probably muon
			else                 return 2; // probably pion
		}
	}
	else if (idTPC == 3){
		if      (idTOF == 3) return 3; // probably kaon
		else if (idTOF == 4) return 4; // probably proton
		else {
			if      (idTPCpick == 0) return 0; // probably misidentified electron
			else if (idTPCpick == 1) return 1; // probably misidentified muon
			else if (idTPCpick == 2) return 2; // probably misidentified pion
		}
	}
	else {
		if      (idTOF == 3) return 3; // probably kaon
		else if (idTOF == 4) return 4; // probably proton
		else {
			if      (idTPCpick == 0) return 0; // probably misidentified electron
			else if (idTPCpick == 1) return 1; // probably misidentified muon
			else if (idTPCpick == 2) return 2; // probably misidentified pion
		}
	}

	Printf("Something is wrong with your brilliant logic, genius!");
	return -1;
}

//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::TPCandTOFsignalInfo(AliESDtrack *trk, Int_t trkID){

	fPIDpt[trkID] = trk->Pt();
	fTPCsignal[trkID] = trk->GetTPCsignal();
	fTOFsignal[trkID] = trk->GetTOFsignal();

	Float_t PIDTPC[5];
	PIDTPC[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron));
	PIDTPC[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon));
	PIDTPC[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion));
	PIDTPC[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon));
	PIDTPC[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton));

	fTPCmostProbableTrackType[trkID] = std::distance(std::begin(PIDTPC),std::min_element(std::begin(PIDTPC),std::end(PIDTPC)));

	Float_t PIDTOF[5];
	PIDTOF[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron));
	PIDTOF[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon));
	PIDTOF[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion));
	PIDTOF[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon));
	PIDTOF[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton));

	fTOFmostProbableTrackType[trkID] = std::distance(std::begin(PIDTOF),std::min_element(std::begin(PIDTOF),std::end(PIDTOF)));
}

//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::SetCrossed(Int_t spd[4], TBits &crossed){

	Int_t chipId2;
	for(Int_t iLayer = 0; iLayer<4 ;iLayer++)
		if(spd[iLayer]>0) { crossed.SetBitNumber(GetChipId(spd[iLayer],chipId2)); crossed.SetBitNumber(chipId2); }
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskCentralTau::GetChipId(Int_t index, Int_t &chipId2, Bool_t debug){
//  Int_t status   = (index%1000000)/100000;
	Int_t iModule  = index/1000000;           // 0 - 239
	Int_t iPhi     = iModule/4;               // 0-19 - inner, 20-59 outer
//  Int_t iModuleZ = iModule%4;               // 0-3
	Int_t iSign    = (index%100000)/10000;    // 1-4
	Int_t signZ    = iPhi<20 ? (iSign%2==1 ? 1 : -1) : (iSign%2==0 ? 1 : -1); // 1 or -1
	Int_t iX       = (index%10000)/100;       // ??
	Int_t iZ       = index%100;               // 0-36 [mm]
	Int_t signZiZ  = (36-signZ*iZ);
	Int_t chipId   = iModule*5+signZiZ*5/72;
	if (chipId<0) return 1200;
	if (chipId>=1200) return 1201;
	if (signZiZ<0) return 1202;
	if (signZiZ>72) return 1203;
	if (signZiZ==72 && chipId%20==0 && chipId>=400) return 1204;
	chipId2=chipId;

	if (signZiZ==0  && chipId%20!=0)  chipId2=chipId-1;
	if (signZiZ==72 && chipId%20!=19) chipId2=chipId+1;
	if (signZiZ==13)  chipId2=chipId+1;
	if (signZiZ==14)  chipId2=chipId+1;
	if (signZiZ==15)  chipId2=chipId-1;
	if (signZiZ==16)  chipId2=chipId-1;
	if (signZiZ==27)  chipId2=chipId+1;
	if (signZiZ==28)  chipId2=chipId+1;
	if (signZiZ==29)  chipId2=chipId-1;
	if (signZiZ==30)  chipId2=chipId-1;
	if (signZiZ==42)  chipId2=chipId+1;
	if (signZiZ==43)  chipId2=chipId+1;
	if (signZiZ==44)  chipId2=chipId-1;
	if (signZiZ==45)  chipId2=chipId-1;
	if (signZiZ==56)  chipId2=chipId+1;
	if (signZiZ==57)  chipId2=chipId+1;
	if (signZiZ==58)  chipId2=chipId-1;
	if (signZiZ==59)  chipId2=chipId-1;
	if (debug) printf("%4i %4i %3i %3i %3i\n",chipId,chipId2,iX,signZiZ,iSign);
	return chipId;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCentralTau::IsSTGFired(const TBits& bits, Int_t dphiMin, Int_t dphiMax, Bool_t tolerance){
	UInt_t n1 = bits.CountBits(400);
	UInt_t n0 = bits.CountBits()-n1;
	//cout<<n0<<" "<<n1<<endl;
	if (n0<1 || n1<1) return false;
	Bool_t stg = false;
	Bool_t l0[20]={false};
	Bool_t l1[40]={false};
	Bool_t phi[20]={false};
	for (Int_t i=0;   i< 400; ++i) if (bits.TestBitNumber(i)) l0[      i/20] = true;
	for (Int_t i=400; i<1200; ++i) if (bits.TestBitNumber(i)) l1[(i-400)/20] = true;
	for (Int_t i=0; i<20; ++i) {
		if (tolerance) phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40] | l1[(2*i+2)%40] | l1[(2*i+39)%40]);
		else           phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40]);
	}
	for (Int_t dphi=dphiMin;dphi<=dphiMax;dphi++)
		for (Int_t i=0; i<20; ++i) stg |= phi[i] & phi[(i+dphi)%20];
	return stg;
}




//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::Terminate(Option_t *)
{
	cout<<"Analysis complete."<<endl;
}//Terminate

