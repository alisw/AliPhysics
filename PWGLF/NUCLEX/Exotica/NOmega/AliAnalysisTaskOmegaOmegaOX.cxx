#include "TSystem.h"
#include "TROOT.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliStack.h"
#include "AliAnalysisTask.h"

#include "AliCentrality.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliNeutralTrackParam.h"
#include "AliNormalizationCounter.h"
#include "AliPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliInputEventHandler.h"
#include "AliVertexerTracks.h"
#include "AliEventPoolManager.h"

#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliAODPidHF.h"

#include "AliESD.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliCascadeVertexer.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskOmegaOmegaOX.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskOmegaOmegaOX)

//______________________________________________________________________________________________________
AliAnalysisTaskOmegaOmegaOX::AliAnalysisTaskOmegaOmegaOX() : 
	AliAnalysisTaskSE(),
	fESDEvent(0x0),
	fESDtrackCutsV0(0x0),
	fESDtrackCuts(0x0),
	fUseMCInfo(kFALSE),
	fPIDResponse(0x0),
	fIsEventSelected(kFALSE),
	fMixedEvent(kFALSE),
	fParametersTree(0x0),
	fVariablesTree(0x0),
	fParameters(),
	fCandidateVariables(),
	fVtx1(0x0),
	fBzkG(0),
	fRecoTypeDB(1),// Reconstruction type of DiBaryon (0:All, 1:OmOm, 2:OmXi, 3:XiOm, 4:XiXi)
	fLikeSignDB(0),// Like-sign of DB (0:ALL, 1:OO, 2:(Obar)(Obar))
	fRecoSelfCasc(1),// Cascade reconstruction is made by (0:ESD class, 1:by myself(with ESD), 2:by myself(with string))
	fReqSigmaTPC(3.0),
	fReqClustersTPC(80),
	fReqSigmaTOF(3.0),
	fReqPseudoRap(0.9),
	fProtonPMax(999.),// Max momentum of proton
	fPionPMax(1.5),//!!CHECK ALSO ESDtrackCuts!! Max momentum of pion
	fKaonPMax(1.5),//0.6!!CHECK ALSO ESDtrackCuts!! Max momentum of kaon
	fCPADibaryon(0.9),//0.99875
	fDCADibaryon(1.0),
	fMassWinCascade(999.),
	fCsReqChi2max(33.),
	fCsReqDV0min(0.2),//0.01
	fCsReqMassWinLambda(0.008),
	fCsReqDBachMin(0.1),//0.01
	fCsReqDCAmax(1.0),//2.0
	fCsReqCPAmin(0.5),//0.98
	fCsReqRmin(0.5),//0.2
	fCsReqRmax(100.),
	fCentrality(0),
	fCountEvent(0),
	fIsDCA(0),
	fNAll(0)
{
}
//______________________________________________________________________________________________________
AliAnalysisTaskOmegaOmegaOX::AliAnalysisTaskOmegaOmegaOX(const Char_t* name) :
	AliAnalysisTaskSE(name),
	fESDEvent(0x0),
	fESDtrackCutsV0(0x0),
//	fESDCutsV0(0),
	fESDtrackCuts(0x0),
	fUseMCInfo(kFALSE),
	fPIDResponse(0x0),
	fIsEventSelected(kFALSE),
	fMixedEvent(kFALSE),
	fParametersTree(0),
	fVariablesTree(0),
	fParameters(),
	fCandidateVariables(),
	fVtx1(0),
	fBzkG(0),
	fRecoTypeDB(1),// Reconstruction type of DiBaryon (0:All, 1:OmOm, 2:OmXi, 3:XiOm, 4:XiXi)
	fLikeSignDB(0),// Like-sign of DB (0:ALL, 1:OO, 2:(Obar)(Obar))
	fRecoSelfCasc(1),// Cascade reconstruction is made by (0:ESD class, 1:by myself(with ESD), 2:by myself(with string))
	fReqSigmaTPC(3.0),
	fReqClustersTPC(80),
	fReqSigmaTOF(3.0),
	fReqPseudoRap(0.9),
	fProtonPMax(999.),// Max momentum of proton
	fPionPMax(1.5),//!!CHECK ALSO ESDtrackCuts!! Max momentum of pion
	fKaonPMax(1.5),//0.6!!CHECK ALSO ESDtrackCuts!! Max momentum of kaon
	fCPADibaryon(0.9),//0.99875
	fDCADibaryon(1.0),
	fMassWinCascade(999.),
	fCsReqChi2max(33.),
	fCsReqDV0min(0.2),//0.01
	fCsReqMassWinLambda(0.008),
	fCsReqDBachMin(0.1),//0.01
	fCsReqDCAmax(1.0),//2.0
	fCsReqCPAmin(0.5),//0.98
	fCsReqRmin(0.5),//0.2
	fCsReqRmax(100.),
	fCentrality(0),
	fCountEvent(0),
	fIsDCA(0),
	fNAll(0)
{
	//
	// Constructor. Initialization of Inputs and Outputs
	//
	Info("AliAnalysisTaskOmegaOmegaOX","Calling Constructor");

	DefineOutput(1,TTree::Class());  //My private output
	DefineOutput(2,TTree::Class());  //My private output


	//V0 cuts
	fESDtrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
	fESDtrackCutsV0->SetMinNClustersTPC(80);
	fESDtrackCutsV0->SetAcceptKinkDaughters(kFALSE);
	fESDtrackCutsV0->SetMaxChi2PerClusterTPC(5);
	fESDtrackCutsV0->SetRequireTPCRefit(kTRUE);
	fESDtrackCutsV0->SetEtaRange(-0.9,0.9);
	fESDtrackCutsV0->SetPtRange(0.2,1.5);//pi
//	fESDtrackCutsV0->SetPtRange(0.2,0.6);//K
	fESDtrackCutsV0->SetMinDCAToVertexXY(2); //war inzwischen 1 & 3
	fESDtrackCutsV0->SetMinDCAToVertexZ(2); //war inzwischen 1 & 3

//	fESDCutsV0 = new AliESDv0Cuts("AliESDCutsV0","AliESDCutsV0");
//	fESDCutsV0->SetMaxDcaV0Daughters(1.0);
//	fESDCutsV0->SetMinDcaNegToVertex(2); //1.5
//	fESDCutsV0->SetMinDcaPosToVertex(2); //1.5

	//ESD Track cuts
	fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
	fESDtrackCuts->SetMinNClustersTPC(80);
	fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
	fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
	fESDtrackCuts->SetRequireTPCRefit(kTRUE);
	fESDtrackCuts->SetEtaRange(-0.9,0.9);

}
//______________________________________________________________________________________________________
AliAnalysisTaskOmegaOmegaOX::~AliAnalysisTaskOmegaOmegaOX() {
	//
	// destructor
	//
	Info("~AliAnalysisTaskOmegaOmegaOX","Calling Destructor");
  
	if (fPIDResponse) {
		delete  fPIDResponse;
	}

	if (fParametersTree) {
		delete fParametersTree;
		fParametersTree = 0;
	}

	if (fVariablesTree) {
		delete fVariablesTree;
		fVariablesTree = 0;
	}

	if(fESDtrackCutsV0) delete fESDtrackCutsV0;
//	if(fESDCutsV0) delete fESDCutsV0;
	if(fESDtrackCuts) delete fESDtrackCuts;

}
//______________________________________________________________________________________________________
void AliAnalysisTaskOmegaOmegaOX::Init() {
	//
	// Initialization
	//
	//

	//Copied from $ALICE_ROOT/PWGHF/vertexingHF/ConfigVertexingHF.C

	fIsEventSelected=kFALSE;

	if (fDebug > 1) AliInfo("Init");

	return;
}
//______________________________________________________________________________________________________
void AliAnalysisTaskOmegaOmegaOX::UserExec(Option_t *)
{
// user exec
	//----------------------
	// Check the PIDresponse
	//----------------------
	if(!fPIDResponse) {
		AliError("No pid response");
		return;
	}

	//---------------------
	// Check the InputEvent
	//---------------------
	fESDEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
	if (!fESDEvent) {
		AliError("NO EVENT FOUND!");
		return;
	}

	//------------------------------------------------
	// First check if the event has proper vertex and B
	//------------------------------------------------
	AliESDVertex *esdTrkPV = (AliESDVertex*)fESDEvent->GetPrimaryVertexTracks();
	Double_t PosTrkPV[3];
	esdTrkPV->GetXYZ(PosTrkPV);

	AliESDVertex *esdPV = (AliESDVertex*)fESDEvent->GetPrimaryVertex();
	Double_t posPV[3];
	esdPV->GetXYZ(posPV);

	if (TMath::Abs(PosTrkPV[2])>10.) {
		AliError("Z of primary vertex is not within +-10cm !");
		return;
	}

	Int_t runNumber = 0;
	runNumber = (InputEvent())->GetRunNumber();

	fBzkG = (Double_t)fESDEvent->GetMagneticField(); 
	AliKFParticle::SetField(fBzkG);
	if (TMath::Abs(fBzkG)<0.001) {
		AliError("No magnet field !");
		return;
	}

	//------------------------------------------------
  // MC analysis setting
  //------------------------------------------------
  TClonesArray *mcArray = 0;
  Int_t nSelectedAnal = 0;

	//------------------------------------------------
  // Main analysis done in this function
  //------------------------------------------------
	MakeAnalysis(mcArray, fESDEvent);

	PostData(1,fParametersTree);
	PostData(2,fVariablesTree);

  fIsEventSelected=kFALSE;

  return;
}

//______________________________________________________________________________________________________
void AliAnalysisTaskOmegaOmegaOX::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //AliInfo("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  return;
}

//______________________________________________________________________________________________________
void AliAnalysisTaskOmegaOmegaOX::UserCreateOutputObjects() 
{ 
  // output
  AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));

	DefineTreeVariables();
	PostData(1,fParametersTree);
	PostData(2,fVariablesTree);

	// I don't want to use the PID through the cut object, 
	// but I will use the PID response directly!!!
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  return;
}
//______________________________________________________________________________________________________
void AliAnalysisTaskOmegaOmegaOX::MakeAnalysis(TClonesArray *mcArray,AliESDEvent *fESDEvent)
{


  //------------------------------------------------------------------------------------------
  // version OO2-0-10 (2016/04/13)
  // Reconstruct 2 cascade by myself
  // and calculate the invariant mass of Omega-Omega
  //------------------------------------------------------------------------------------------

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Set cut parameter
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// Others ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Int_t setStartNumber = 0;
//	fMixedEvent = kTRUE;
	fMixedEvent = kFALSE;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Initialization
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// Get or calculate constant ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Double_t mElectronPDG = TDatabasePDG::Instance()->GetParticle(11)->Mass();//0.000511
	Double_t mPionPDG     = TDatabasePDG::Instance()->GetParticle(211)->Mass();//0.139570
	Double_t mKaonPDG     = TDatabasePDG::Instance()->GetParticle(321)->Mass();//0.493677
	Double_t mProtonPDG   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();//0.938272
	Double_t mLambdaPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();//1.115680
  Double_t mXiPDG       = TDatabasePDG::Instance()->GetParticle(3312)->Mass();//1.321710
//  Double_t mXi1530PDG   = TDatabasePDG::Instance()->GetParticle(3314)->Mass();//1.535000
  Double_t mOmegaPDG    = TDatabasePDG::Instance()->GetParticle(3334)->Mass();//1.672450

  const Int_t nTracks = fESDEvent->GetNumberOfTracks();
  if (nTracks==0) {
    return;
  }

  const Int_t nV0s = fESDEvent->GetNumberOfV0s();
  if (nV0s==0) {
    return;
  }

  const Int_t nCascades = fESDEvent->GetNumberOfCascades();

//	printf("### nTracks:%6d, nV0s:%6d, nCascades:%6d\n",nTracks,nV0s,nCascades);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Output cut parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fParameters[ 0] = fRecoTypeDB;
	fParameters[ 1] = fLikeSignDB;
	fParameters[ 2] = fReqSigmaTPC;
	fParameters[ 3] = fReqClustersTPC;
	fParameters[ 4] = fReqSigmaTOF;
	fParameters[ 5] = fReqPseudoRap;
	fParameters[ 6] = fCPADibaryon;
	fParameters[ 7] = fMassWinCascade;
  fParameters[ 8] = fCsReqChi2max;
  fParameters[ 9] = fCsReqDV0min;
  fParameters[10] = fCsReqMassWinLambda;
  fParameters[11] = fCsReqDBachMin;
  fParameters[12] = fCsReqDCAmax; 
  fParameters[13] = fCsReqCPAmin;
  fParameters[14] = fCsReqRmin;
  fParameters[15] = fCsReqRmax;
	fParameters[16] = fRecoSelfCasc;

  fParametersTree->Fill();
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Get or calculate constant for this event ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const AliESDVertex *esdTrkPV = (AliESDVertex*)fESDEvent->GetPrimaryVertexTracks();
	Double_t PosTrkPV[3];
	esdTrkPV->GetXYZ(PosTrkPV);

	const AliESDVertex *esdPV = (AliESDVertex*)fESDEvent->GetPrimaryVertex();
	Double_t PosPV[3];
	esdPV->GetXYZ(PosPV);

	Double_t vPVX = PosPV[0];
	Double_t vPVY = PosPV[1];
	Double_t vPVZ = PosPV[2];

	Int_t runNumber;
	AliCentrality *cent = fESDEvent->GetCentrality();
	fCentrality = cent->GetCentralityClass10("V0M");
	runNumber = fESDEvent->GetRunNumber();

	Int_t eventNumber = fESDEvent->GetEventNumberInFile();
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // START ANALYSIS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------
	// Cut for event
	//------------------------------------------------------------------------------------------

	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
	Bool_t isSelectedCn = (inputHandler->IsEventSelected()&AliVEvent::kCentral);
	Bool_t isSelectedSC = (inputHandler->IsEventSelected()&AliVEvent::kSemiCentral);
	Bool_t isSelectedMB = (inputHandler->IsEventSelected()&AliVEvent::kMB);
	Int_t triggerType = 0;
	if ( isSelectedCn ) triggerType = 1;
	if ( isSelectedSC ) triggerType = 2;
	if ( isSelectedMB ) triggerType = 3;

//	if (fCountEvent<setStartNumber) return;
//	if (fCountEvent>10) return;


	//------------------------------------------------------------------------------------------
	// Get the number of cascades (Reconstruct cascade by myself)
	//------------------------------------------------------------------------------------------

	Int_t nCascadeOri = fESDEvent->GetNumberOfCascades();
	Int_t nCasc;

	const Int_t nCascComb = 1;
//	printf("### number of str:%d(V0:%6d,nTrk:%6d)\n",nCascComb,nV0s,nTracks);
	Int_t v0IDCascCount[nCascComb];
	Int_t trkIDCascCount[nCascComb];

	if ( fRecoSelfCasc==1 ) {//Reconstruction of cascade is done by myself(with ESD)
		if (!ReconstructCascade()) return;
		nCasc = 1;
	}
	if ( fRecoSelfCasc==2 ) {//Reconstruction of cascade is done by myself(with string)
		// Count the number of reconstructed cascades
		return;
		//FIXME
		if (!ReconstructCascadeString(kFALSE,nCasc,v0IDCascCount,trkIDCascCount)) return;
	}

	const Int_t nStrCasc = nCasc;
	Int_t v0IDCasc[nStrCasc];
	Int_t trkIDCasc[nStrCasc];

	if ( fRecoSelfCasc==2 ) {//Reconstruction of cascade is done by myself(with string)
		if (!ReconstructCascadeString(kTRUE,nCasc,v0IDCasc,trkIDCasc)) return;
	}

	if ( fRecoSelfCasc==0 || fRecoSelfCasc==1 ) {//with ESD
		nCasc = fESDEvent->GetNumberOfCascades();
	}

  if (nCasc==0) {
    return;
  }
//	printf("\n### reconstructed nCascade:%d, ESDCascade:%d\n\n",nCasc,nCascadeOri);

	//------------------------------------------------------------------------------------------
	// Cascade loop 1 (To find Omega1) (START)
	//------------------------------------------------------------------------------------------

	for (Int_t iCasc1=0; iCasc1<nCasc; iCasc1++) {
		AliESDcascade *casc1 = fESDEvent->GetCascade(iCasc1);

		// Get track information
		AliESDtrack *trkP1 = fESDEvent->GetTrack(casc1->GetPindex());
		AliESDtrack *trkN1 = fESDEvent->GetTrack(casc1->GetNindex());
		AliESDtrack *trkB1 = fESDEvent->GetTrack(casc1->GetBindex());

		Double_t trkP1Charge = trkP1->GetSign();
		Double_t trkN1Charge = trkN1->GetSign();
		Double_t trkB1Charge = trkB1->GetSign();
		if ( trkP1Charge>0 && trkN1Charge<0 ) {
			// Nothing to do
		} else if ( trkP1Charge<0 && trkN1Charge>0 ) {
			// P <-> N
			trkP1 = fESDEvent->GetTrack(casc1->GetNindex());
			trkN1 = fESDEvent->GetTrack(casc1->GetPindex());
		} else continue;

		Int_t trkP1ID = trkP1->GetID();
		Int_t trkN1ID = trkN1->GetID();
		Int_t trkB1ID = trkB1->GetID();

		if( !fESDtrackCuts->AcceptTrack(trkP1) ) continue;
		if( !fESDtrackCuts->AcceptTrack(trkN1) ) continue;
		if( !fESDtrackCuts->AcceptTrack(trkB1) ) continue;

		// Pre-track cut (TPCsigma) (p+,pi-,K- or p-,pi+,K+)
		Double_t trkP1TPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkP1, AliPID::kProton));
		Double_t trkP1TPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkP1, AliPID::kPion  ));
		Double_t trkN1TPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkN1, AliPID::kProton));
		Double_t trkN1TPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkN1, AliPID::kPion  ));
		Double_t trkB1TPCKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkB1, AliPID::kKaon  ));
		Double_t trkB1TPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkB1, AliPID::kPion  ));

//		if ( trkP1TPCProton>fReqSigmaTPC && trkP1TPCPion>fReqSigmaTPC ) continue;
//		if ( trkN1TPCProton>fReqSigmaTPC && trkN1TPCPion>fReqSigmaTPC ) continue;
		if ( trkB1TPCKaon>fReqSigmaTPC   && trkB1TPCPion>fReqSigmaTPC ) continue;
		if ( fRecoTypeDB==1||fRecoTypeDB==2 ) {//Omega-???
			if ( trkB1TPCKaon>fReqSigmaTPC ) continue;
		} else if ( fRecoTypeDB==3||fRecoTypeDB==4 ) {//Xi-???
			if ( trkB1TPCPion>fReqSigmaTPC ) continue;
		}
		if ( !(   (trkP1TPCProton<=fReqSigmaTPC && trkN1TPCPion<=fReqSigmaTPC  )||
		          (trkP1TPCPion<=fReqSigmaTPC   && trkN1TPCProton<=fReqSigmaTPC)   ) ) {
			continue;
		}


	//------------------------------------------------------------------------------------------
	// Cascade loop 2 (To find Omega2) (START)
	//------------------------------------------------------------------------------------------

		for (Int_t iCasc2=iCasc1+1; iCasc2<nCasc; iCasc2++) {
			AliESDcascade *casc2 = fESDEvent->GetCascade(iCasc2);

			// Get track information
			AliESDtrack *trkP2 = fESDEvent->GetTrack(casc2->GetPindex());
			AliESDtrack *trkN2 = fESDEvent->GetTrack(casc2->GetNindex());
			AliESDtrack *trkB2 = fESDEvent->GetTrack(casc2->GetBindex());

			Double_t trkP2Charge = trkP2->GetSign();
			Double_t trkN2Charge = trkN2->GetSign();
			Double_t trkB2Charge = trkB2->GetSign();
			if ( trkP2Charge>0 && trkN2Charge<0 ) {
				// Nothing to do
			} else if ( trkP2Charge<0 && trkN2Charge>0 ) {
				// P <-> N
				trkP2 = fESDEvent->GetTrack(casc2->GetNindex());
				trkN2 = fESDEvent->GetTrack(casc2->GetPindex());
			} else continue;

			Int_t trkP2ID = trkP2->GetID();
			Int_t trkN2ID = trkN2->GetID();
			Int_t trkB2ID = trkB2->GetID();

			if( !fESDtrackCuts->AcceptTrack(trkP2) ) continue;
			if( !fESDtrackCuts->AcceptTrack(trkN2) ) continue;
			if( !fESDtrackCuts->AcceptTrack(trkB2) ) continue;

			// Pre-track cut (TPCsigma) (p+,pi-,K- or p-,pi+,K+)
			Double_t trkP2TPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkP2, AliPID::kProton));
			Double_t trkP2TPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkP2, AliPID::kPion  ));
			Double_t trkN2TPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkN2, AliPID::kProton));
			Double_t trkN2TPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkN2, AliPID::kPion  ));
			Double_t trkB2TPCKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkB2, AliPID::kKaon  ));
			Double_t trkB2TPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkB2, AliPID::kPion  ));
//			if ( trkP2TPCProton>fReqSigmaTPC && trkP2TPCPion>fReqSigmaTPC ) continue;
//			if ( trkN2TPCProton>fReqSigmaTPC && trkN2TPCPion>fReqSigmaTPC ) continue;
			if ( trkB2TPCKaon>fReqSigmaTPC   && trkB2TPCPion>fReqSigmaTPC ) continue;
			if ( fRecoTypeDB==1||fRecoTypeDB==3 ) {//???-Omega
				if ( trkB2TPCKaon>fReqSigmaTPC ) continue;
			} else if ( fRecoTypeDB==2||fRecoTypeDB==4 ) {//???-Xi
				if ( trkB2TPCPion>fReqSigmaTPC ) continue;
			}
			if ( !(   (trkP2TPCProton<=fReqSigmaTPC && trkN2TPCPion<=fReqSigmaTPC  )||
			          (trkP2TPCPion<=fReqSigmaTPC   && trkN2TPCProton<=fReqSigmaTPC)   ) ) continue;

	//------------------------------------------------------------------------------------------
	// Preparation to reconstruct Omega-Omega
	//------------------------------------------------------------------------------------------
//			printf("### loop casc1:%6d, casc2:%6d\n",iCasc1,iCasc2);

			if ( trkP1ID==trkN1ID || trkN1ID==trkB1ID || trkB1ID==trkP1ID ) continue;
			if ( trkP2ID==trkN2ID || trkN2ID==trkB2ID || trkB2ID==trkP2ID ) continue;

			if ( trkP1ID==trkP2ID || trkP1ID==trkN2ID || trkP1ID==trkB2ID ) continue;
			if ( trkN1ID==trkP2ID || trkN1ID==trkN2ID || trkN1ID==trkB2ID ) continue;
			if ( trkB1ID==trkP2ID || trkB1ID==trkN2ID || trkB1ID==trkB2ID ) continue;

//			if ( trkP1TPCProton>fReqSigmaTPC && trkN1TPCProton>fReqSigmaTPC ) continue;
//			if ( trkP1TPCPion  >fReqSigmaTPC && trkN1TPCPion  >fReqSigmaTPC ) continue;
//			if ( trkP2TPCProton>fReqSigmaTPC && trkN2TPCProton>fReqSigmaTPC ) continue;
//			if ( trkP2TPCPion  >fReqSigmaTPC && trkN2TPCPion  >fReqSigmaTPC ) continue;

////			if( !fESDtrackCuts->AcceptTrack(trkP1) ) continue;
////			if( !fESDtrackCuts->AcceptTrack(trkP2) ) continue;
////			if( !fESDtrackCuts->AcceptTrack(trkN1) ) continue;
////			if( !fESDtrackCuts->AcceptTrack(trkN2) ) continue;
////			if( !fESDtrackCuts->AcceptTrack(trkB1) ) continue;
////			if( !fESDtrackCuts->AcceptTrack(trkB2) ) continue;

	//------------------------------------------------------------------------------------------
	// Get Cascade momenta and positions
	//------------------------------------------------------------------------------------------
//casc1->Print();
//casc2->Print();

			Double_t PosCasc1[3], PosCasc2[3];
			Double_t MomCasc1[3], MomCasc2[3];
			Double_t PosV01[3], PosV02[3];
			Double_t MomP1[3], MomP2[3];
			Double_t MomN1[3], MomN2[3];
			Double_t MomB1[3], MomB2[3];
			Double_t MomL1[3], MomL2[3];
			Double_t MomO1[3], MomO2[3];

			casc1->XvYvZv(PosCasc1);
			casc2->XvYvZv(PosCasc2);
			casc1->PxPyPz(MomCasc1);
			casc2->PxPyPz(MomCasc2);
			casc1->GetXYZ(PosV01[0],PosV01[1],PosV01[2]);
			casc2->GetXYZ(PosV02[0],PosV02[1],PosV02[2]);
//			printf("### cascpos: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",PosCasc1[0],PosCasc1[1],PosCasc1[2],PosCasc2[0],PosCasc2[1],PosCasc2[2]);
//			printf("### cascmom: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",MomCasc1[0],MomCasc1[1],MomCasc1[2],MomCasc2[0],MomCasc2[1],MomCasc2[2]);
//			printf("### v0  pos: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",PosV01[0],PosV01[1],PosV01[2],PosV02[0],PosV02[1],PosV02[2]);

			casc1->GetPPxPyPz(MomP1[0],MomP1[1],MomP1[2]);
			casc2->GetPPxPyPz(MomP2[0],MomP2[1],MomP2[2]);
			casc1->GetNPxPyPz(MomN1[0],MomN1[1],MomN1[2]);
			casc2->GetNPxPyPz(MomN2[0],MomN2[1],MomN2[2]);
			casc1->GetBPxPyPz(MomB1[0],MomB1[1],MomB1[2]);
			casc2->GetBPxPyPz(MomB2[0],MomB2[1],MomB2[2]);

			for (Int_t i=0; i<3; i++) {
				MomL1[i] = MomP1[i] + MomN1[i];
				MomL2[i] = MomP2[i] + MomN2[i];
				MomO1[i] = MomL1[i] + MomB1[i];
				MomO2[i] = MomL2[i] + MomB2[i];
			}

			Double_t MomP1P = GetPaFromPxPyPz(MomP1);
			Double_t MomP2P = GetPaFromPxPyPz(MomP2);
			Double_t MomN1P = GetPaFromPxPyPz(MomN1);
			Double_t MomN2P = GetPaFromPxPyPz(MomN2);
			Double_t MomB1P = GetPaFromPxPyPz(MomB1);
			Double_t MomB2P = GetPaFromPxPyPz(MomB2);

			Double_t Casc1toV01[3], Casc2toV02[3];
			for (Int_t i=0; i<3; i++) {
				Casc1toV01[i] = PosV01[i] - PosCasc1[i];
				Casc2toV02[i] = PosV02[i] - PosCasc2[i];
			}
			Double_t rV01    = GetPaFromPxPyPz(Casc1toV01); 
			Double_t rV02    = GetPaFromPxPyPz(Casc2toV02); 
			Double_t ctauV01 = rV01*mLambdaPDG/GetPaFromPxPyPz(MomL1);
			Double_t ctauV02 = rV02*mLambdaPDG/GetPaFromPxPyPz(MomL2);

//			printf("### P   mom: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",MomP1[0],MomP1[1],MomP1[2],MomP2[0],MomP2[1],MomP2[2]);
//			printf("### N   mom: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",MomN1[0],MomN1[1],MomN1[2],MomN2[0],MomN2[1],MomN2[2]);
//			printf("### B   mom: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",MomB1[0],MomB1[1],MomB1[2],MomB2[0],MomB2[1],MomB2[2]);
//			printf("### L   mom: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",MomL1[0],MomL1[1],MomL1[2],MomL2[0],MomL2[1],MomL2[2]);
//			printf("\n#############################################\n");
//			printf("### O   mom: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",MomCasc1[0],MomCasc1[1],MomCasc1[2],MomCasc2[0],MomCasc2[1],MomCasc2[2]);
//			printf("### Ori mom: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",MomO1[0],MomO1[1],MomO1[2],MomO2[0],MomO2[1],MomO2[2]);
//			printf("### O   pos: 1(%10f,%10f,%10f), 2(%10f,%10f,%10f)\n",PosCasc1[0],PosCasc1[1],PosCasc1[2],PosCasc2[0],PosCasc2[1],PosCasc2[2]);

			// calculate DCA of two daughters of dibaryon
//			AliExternalTrackParam *trkLikeCasc1 = new AliExternalTrackParam();
//			AliExternalTrackParam *trkLikeCasc2 = new AliExternalTrackParam();

			Double_t cvDef[21];
			for (Int_t i=0; i<21; i++) cvDef[i] = 0.;
			cvDef[ 0] = 1.;
			cvDef[ 2] = 1.;
			cvDef[ 5] = 1.;
			cvDef[ 9] = 1.;
			cvDef[14] = 1.;
			cvDef[20] = 1.;

			AliExternalTrackParam trkLikeCasc1(PosCasc1,MomCasc1,cvDef,static_cast<Short_t>(trkB1Charge));
			AliExternalTrackParam trkLikeCasc2(PosCasc2,MomCasc2,cvDef,static_cast<Short_t>(trkB2Charge));

			Double_t xDcaCasc1, xDcaCasc2;
			Double_t DCADB; 
			DCADB = trkLikeCasc1.GetDCA(&trkLikeCasc2,fBzkG,xDcaCasc1,xDcaCasc2);
			if ( DCADB > fDCADibaryon ) continue;
//			printf("##### DCADB:%10f(<%10f)\n",DCADB,fDCADibaryon);

			trkLikeCasc1.PropagateTo(xDcaCasc1,fBzkG);
			trkLikeCasc2.PropagateTo(xDcaCasc2,fBzkG);

			Double_t MomCasc1DB[3], MomCasc2DB[3];
			trkLikeCasc1.GetPxPyPz(MomCasc1DB);
			trkLikeCasc2.GetPxPyPz(MomCasc2DB);

			Double_t PosCasc1DB[3], PosCasc2DB[3];
			trkLikeCasc1.GetXYZ(PosCasc1DB);
			trkLikeCasc2.GetXYZ(PosCasc2DB);

			// get dibaryon momentum and position
			Double_t MomDB[3], PosDB[3], PVtoDB[3];
			for (Int_t i=0; i<3; i++) {
				MomDB[i]  = MomCasc1DB[i] + MomCasc2DB[i];
				PosDB[i]  = (PosCasc1DB[i] + PosCasc2DB[i]) / 2.;
				PVtoDB[i] = PosDB[i] - PosPV[i];
			}

			Double_t rDB    = GetPaFromPxPyPz(PVtoDB);
			Double_t ctauDB = rDB/(2*mOmegaPDG)/GetPaFromPxPyPz(MomDB);
			Double_t cpaDB  = (PVtoDB[0]*MomDB[0]+PVtoDB[1]*MomDB[1]+PVtoDB[2]*MomDB[2])/rDB/GetPaFromPxPyPz(MomDB);
			if ( cpaDB < fCPADibaryon ) continue;
//			printf("##### CPADB:%10f\n",cpaDB);

			Double_t DBtoCasc1[3], DBtoCasc2[3];
			for (Int_t i=0; i<3; i++) {
				DBtoCasc1[i] = PosCasc1[i] - PosDB[i];
				DBtoCasc2[i] = PosCasc2[i] - PosDB[i];
			}
			Double_t rCasc1 = GetPaFromPxPyPz(DBtoCasc1); 
			Double_t rCasc2 = GetPaFromPxPyPz(DBtoCasc2); 
			Double_t ctauO1 = rCasc1*mOmegaPDG/GetPaFromPxPyPz(MomCasc1);
			Double_t ctauO2 = rCasc2*mOmegaPDG/GetPaFromPxPyPz(MomCasc2);
//			printf("##### CTAUO1:%10f, CTAUO2:%10f\n",ctauO1,ctauO2);

//			Double_t CPACasc1Test = casc1->GetCascadeCosineOfPointingAngle(PosDB[0],PosDB[1],PosDB[2]);
//			Double_t CPACasc2Test = casc2->GetCascadeCosineOfPointingAngle(PosDB[0],PosDB[1],PosDB[2]);
			Double_t CPACasc1Test = 0.;
			Double_t CPACasc2Test = 0.;
			Double_t CPACasc1 = (DBtoCasc1[0]*MomCasc1[0]+DBtoCasc1[1]*MomCasc1[1]+DBtoCasc1[2]*MomCasc1[2])/rCasc1/GetPaFromPxPyPz(MomCasc1); 
			Double_t CPACasc2 = (DBtoCasc2[0]*MomCasc2[0]+DBtoCasc2[1]*MomCasc2[1]+DBtoCasc2[2]*MomCasc2[2])/rCasc2/GetPaFromPxPyPz(MomCasc2);
//			printf("### CPAcasc1:%10f, ESD:%10f, CPAcasc2:%10f, ESD:%10f\n",CPACasc1,CPACasc1Test,CPACasc2,CPACasc2Test);


	//------------------------------------------------------------------------------------------
	// PID information
	//------------------------------------------------------------------------------------------

			// proLambda or antiLambda
			Int_t proLambda1  = 0;
			Int_t proLambda2  = 0;

			if ( trkP1TPCProton<fReqSigmaTPC && trkN1TPCPion  <fReqSigmaTPC ) proLambda1 =  1;
			else if ( trkP1TPCPion  <fReqSigmaTPC && trkN1TPCProton<fReqSigmaTPC ) proLambda1 = -1;
			if ((trkP1TPCProton<fReqSigmaTPC && trkN1TPCPion  <fReqSigmaTPC) &&
					(trkP1TPCPion  <fReqSigmaTPC && trkN1TPCProton<fReqSigmaTPC)) {
				Double_t checkProLambda1 = TMath::Sqrt(trkP1TPCProton*trkP1TPCProton +
																							 trkN1TPCPion  *trkN1TPCPion   );
				Double_t checkAntLambda1 = TMath::Sqrt(trkP1TPCPion  *trkP1TPCPion   +
																							 trkN1TPCProton*trkN1TPCProton );
				if ( checkProLambda1 < checkAntLambda1 ) proLambda1 = 1;
				else if ( checkProLambda1 >= checkAntLambda1 ) proLambda1 = -1;
			}

			if ( trkP2TPCProton<fReqSigmaTPC && trkN2TPCPion  <fReqSigmaTPC ) proLambda2 =  1;
			else if ( trkP2TPCPion  <fReqSigmaTPC && trkN2TPCProton<fReqSigmaTPC ) proLambda2 = -1;
			if ((trkP2TPCProton<fReqSigmaTPC && trkN2TPCPion  <fReqSigmaTPC) &&
					(trkP2TPCPion  <fReqSigmaTPC && trkN2TPCProton<fReqSigmaTPC)) {
				Double_t checkProLambda2 = TMath::Sqrt(trkP2TPCProton*trkP2TPCProton +
																							 trkN2TPCPion  *trkN2TPCPion   );
				Double_t checkAntLambda2 = TMath::Sqrt(trkP2TPCPion  *trkP2TPCPion   +
																							 trkN2TPCProton*trkN2TPCProton );
				if ( checkProLambda2 < checkAntLambda2 ) proLambda2 = 1;
				else if ( checkProLambda2 >= checkAntLambda2 ) proLambda2 = -1;
			}

			// Kaon or pion for Bachelor
			Int_t kaonBach1  = 0;//1: kaon, -1:pion
			Int_t kaonBach2  = 0;//1: kaon, -1:pion

			if ( trkB1TPCKaon<fReqSigmaTPC ) kaonBach1 =  1;
			else if ( trkB1TPCPion<fReqSigmaTPC ) kaonBach1 = -1;
			if ( trkB1TPCKaon<fReqSigmaTPC && trkB1TPCPion<fReqSigmaTPC ) {
				if ( trkB1TPCKaon<trkB1TPCPion ) kaonBach1 = 1;
				else if ( trkB1TPCKaon>=trkB1TPCPion ) kaonBach1 = -1;
			}

			if ( trkB2TPCKaon<fReqSigmaTPC ) kaonBach2 =  1;
			else if ( trkB2TPCPion<fReqSigmaTPC ) kaonBach2 = -1;
			if ( trkB2TPCKaon<fReqSigmaTPC && trkB2TPCPion<fReqSigmaTPC ) {
				if ( trkB2TPCKaon<trkB2TPCPion ) kaonBach2 = 1;
				else if ( trkB2TPCKaon>=trkB2TPCPion ) kaonBach2 = -1;
			}

			// Xi or Omega
			Int_t typeOfCasc1 = 0;//1: Omega-(LK-), 2:Xi-(Lpi-), 3:LK+, 4:Lpi+, negative:Anti-
			Int_t typeOfCasc2 = 0;//1: Omega-(LK-), 2:Xi-(Lpi-), 3:LK+, 4:Lpi+, negative:Anti-

			if ( proLambda1 == 1 ) {//pro-Lambda
				if ( trkB1Charge == -1 ) {//minus bachelor
					if ( kaonBach1 ==  1 ) typeOfCasc1 = 1;//Omega-
					else if ( kaonBach1 == -1 ) typeOfCasc1 = 2;//Xi-
					else typeOfCasc1 = -99;
				} else if ( trkB1Charge ==  1 ) {//plus bachelor
					if ( kaonBach1 ==  1 ) typeOfCasc1 = 3;//LK+
					else if ( kaonBach1 == -1 ) typeOfCasc1 = 4;//Lpi+
					else typeOfCasc1 = -99;
				} else typeOfCasc1 = -99;
			} else if ( proLambda1 == -1 ) {//anti-Lambda
				if ( trkB1Charge == -1 ) {//minus bachelor
					if ( kaonBach1 ==  1 ) typeOfCasc1 = -3;//LbarK-
					else if ( kaonBach1 == -1 ) typeOfCasc1 = -4;//Lbarpi-
					else typeOfCasc1 = -99;
				} else if ( trkB1Charge ==  1 ) {//plus bachelor
					if ( kaonBach1 ==  1 ) typeOfCasc1 = -1;//Omega+
					else if ( kaonBach1 == -1 ) typeOfCasc1 = -2;//Xi+
					else typeOfCasc1 = -99;
				} else typeOfCasc1 = -99;
			} else typeOfCasc1 = -99;

			if ( proLambda2 == 1 ) {//pro-Lambda
				if ( trkB2Charge == -1 ) {//minus bachelor
					if ( kaonBach2 ==  1 ) typeOfCasc2 = 1;//Omega-
					else if ( kaonBach2 == -1 ) typeOfCasc2 = 2;//Xi-
					else typeOfCasc2 = -99;
				} else if ( trkB2Charge ==  1 ) {//plus bachelor
					if ( kaonBach2 ==  1 ) typeOfCasc2 = 3;//LK+
					else if ( kaonBach2 == -1 ) typeOfCasc2 = 4;//Lpi+
					else typeOfCasc2 = -99;
				} else typeOfCasc2 = -99;
			} else if ( proLambda2 == -1 ) {//anti-Lambda
				if ( trkB2Charge == -1 ) {//minus bachelor
					if ( kaonBach2 ==  1 ) typeOfCasc2 = -3;//LbarK-
					else if ( kaonBach2 == -1 ) typeOfCasc2 = -4;//Lbarpi-
					else typeOfCasc2 = -99;
				} else if ( trkB2Charge ==  1 ) {//plus bachelor
					if ( kaonBach2 ==  1 ) typeOfCasc2 = -1;//Omega+
					else if ( kaonBach2 == -1 ) typeOfCasc2 = -2;//Xi+
					else typeOfCasc2 = -99;
				} else typeOfCasc2 = -99;
			} else typeOfCasc2 = -99;

			// select dibaryon
			if ( fRecoTypeDB==0 ) {//Default (Nothing to do)
			} else if ( fRecoTypeDB==1 ) {//Omega-Omega
				if ( !(  (typeOfCasc1==1||typeOfCasc1==-1||typeOfCasc1==3||typeOfCasc1==-3)&&
				         (typeOfCasc2==1||typeOfCasc2==-1||typeOfCasc2==3||typeOfCasc2==-3)  ) ) continue;
			} else if ( fRecoTypeDB==2 ) {//Omega-Xi
				if ( !(  (typeOfCasc1==1||typeOfCasc1==-1||typeOfCasc1==3||typeOfCasc1==-3)&&
				         (typeOfCasc2==2||typeOfCasc2==-2||typeOfCasc2==4||typeOfCasc2==-4)  ) ) continue;
			} else if ( fRecoTypeDB==3 ) {//Xi-Omega
				if ( !(  (typeOfCasc1==2||typeOfCasc1==-2||typeOfCasc1==4||typeOfCasc1==-4)&&
				         (typeOfCasc2==1||typeOfCasc2==-1||typeOfCasc2==3||typeOfCasc2==-3)  ) ) continue;
			} else if ( fRecoTypeDB==4 ) {//Omega-Omega
				if ( !(  (typeOfCasc1==2||typeOfCasc1==-2||typeOfCasc1==4||typeOfCasc1==-4)&&
				         (typeOfCasc2==2||typeOfCasc2==-2||typeOfCasc2==4||typeOfCasc2==-4)  ) ) continue;
			} else {
				printf("\n\n################################################################################\n");
				printf("### ERROR!! what is the reconstruct type???                                  ###\n");
				printf("################################################################################\n");
			}
/*
			printf("\n### 1:P(%2d):P=%8f, pi=%8f, N(%2d):P=%8f, pi=%8f, B(%2d):K=%8f, pi=%8f\n",static_cast<Int_t>(trkP1Charge),trkP1TPCProton,trkP1TPCPion,static_cast<Int_t>(trkN1Charge),trkN1TPCProton,trkN1TPCPion,static_cast<Int_t>(trkB1Charge),trkB1TPCKaon,trkB1TPCPion);
			printf("### -> This cascade is ...");
			if (typeOfCasc1== 0) printf(" None\n");
			if (typeOfCasc1== 1) printf(" Omega-\n");
			if (typeOfCasc1== 2) printf(" Xi-\n");
			if (typeOfCasc1== 3) printf(" Lambda K+\n");
			if (typeOfCasc1== 4) printf(" Lambda pi+\n");
			if (typeOfCasc1==-1) printf(" Omega+\n");
			if (typeOfCasc1==-2) printf(" Xi+\n");
			if (typeOfCasc1==-3) printf(" anti-Lambda K-\n");
			if (typeOfCasc1==-4) printf(" anti-Lambda pi-\n");
			if (typeOfCasc1==-99) printf(" error\n");
			printf("### 2:P(%2d):P=%8f, pi=%8f, N(%2d):P=%8f, pi=%8f, B(%2d):K=%8f, pi=%8f\n",static_cast<Int_t>(trkP2Charge),trkP2TPCProton,trkP2TPCPion,static_cast<Int_t>(trkN2Charge),trkN2TPCProton,trkN2TPCPion,static_cast<Int_t>(trkB2Charge),trkB2TPCKaon,trkB2TPCPion);
			printf("### -> This cascade is ...");
			if (typeOfCasc2== 0) printf(" None\n");
			if (typeOfCasc2== 1) printf(" Omega-\n");
			if (typeOfCasc2== 2) printf(" Xi-\n");
			if (typeOfCasc2== 3) printf(" Lambda K+\n");
			if (typeOfCasc2== 4) printf(" Lambda pi+\n");
			if (typeOfCasc2==-1) printf(" Omega+\n");
			if (typeOfCasc2==-2) printf(" Xi+\n");
			if (typeOfCasc2==-3) printf(" anti-Lambda K-\n");
			if (typeOfCasc2==-4) printf(" anti-Lambda pi-\n");
			if (typeOfCasc2==-99) printf(" error\n");
*/

	//------------------------------------------------------------------------------------------
	// Invariant mass
	//------------------------------------------------------------------------------------------

			TLorentzVector vP1,vN1,vB1,vL1,vC1;
			TLorentzVector vP2,vN2,vB2,vL2,vC2;
			TLorentzVector vLm1,vLm2,vCs1,vCs1L,vCs2,vCs2L;
			TLorentzVector vCsCs,vCsCsLL,vCsCsCC;
			TLorentzVector vC1DB,vC1DBL,vC1DBC,vC2DB,vC2DBL,vC2DBC;
			TLorentzVector vDB,vDBLL,vDBCC;
			if      (proLambda1== 1) vP1.SetXYZM(MomP1[0],MomP1[1],MomP1[2],mProtonPDG);
			else if (proLambda1==-1) vP1.SetXYZM(MomP1[0],MomP1[1],MomP1[2],mPionPDG  );
			else                     vP1.SetXYZM(0.,0.,0.,0.);
			if      (proLambda1== 1) vN1.SetXYZM(MomN1[0],MomN1[1],MomN1[2],mPionPDG  );
			else if (proLambda1==-1) vN1.SetXYZM(MomN1[0],MomN1[1],MomN1[2],mProtonPDG);
			else                     vN1.SetXYZM(0.,0.,0.,0.);
			if      (kaonBach1== 1)  vB1.SetXYZM(MomB1[0],MomB1[1],MomB1[2],mKaonPDG);
			else if (kaonBach1==-1)  vB1.SetXYZM(MomB1[0],MomB1[1],MomB1[2],mPionPDG);
			else                     vB1.SetXYZM(0.,0.,0.,0.);
			                         vL1.SetXYZM(MomL1[0],MomL1[1],MomL1[2],mLambdaPDG);
			if      (kaonBach1== 1)  vC1.SetXYZM(MomCasc1[0],MomCasc1[1],MomCasc1[2],mOmegaPDG);
			else if (kaonBach1==-1)  vC1.SetXYZM(MomCasc1[0],MomCasc1[1],MomCasc1[2],mXiPDG   );
			else                     vC1.SetXYZM(0.,0.,0.,0.);

			if      (proLambda2== 1) vP2.SetXYZM(MomP2[0],MomP2[1],MomP2[2],mProtonPDG);
			else if (proLambda2==-1) vP2.SetXYZM(MomP2[0],MomP2[1],MomP2[2],mPionPDG  );
			else                     vP2.SetXYZM(0.,0.,0.,0.);
			if      (proLambda2== 1) vN2.SetXYZM(MomN2[0],MomN2[1],MomN2[2],mPionPDG  );
			else if (proLambda2==-1) vN2.SetXYZM(MomN2[0],MomN2[1],MomN2[2],mProtonPDG);
			else                     vN2.SetXYZM(0.,0.,0.,0.);
			if      (kaonBach2== 1)  vB2.SetXYZM(MomB2[0],MomB2[1],MomB2[2],mKaonPDG);
			else if (kaonBach2==-1)  vB2.SetXYZM(MomB2[0],MomB2[1],MomB2[2],mPionPDG);
			else                     vB2.SetXYZM(0.,0.,0.,0.);
			                         vL2.SetXYZM(MomL2[0],MomL2[1],MomL2[2],mLambdaPDG);
			if      (kaonBach2== 1)  vC2.SetXYZM(MomCasc2[0],MomCasc2[1],MomCasc2[2],mOmegaPDG);
			else if (kaonBach2==-1)  vC2.SetXYZM(MomCasc2[0],MomCasc2[1],MomCasc2[2],mXiPDG   );
			else                     vC2.SetXYZM(0.,0.,0.,0.);

			if      (kaonBach1== 1)  vC1DBC.SetXYZM(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2],mOmegaPDG);
			else if (kaonBach1==-1)  vC1DBC.SetXYZM(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2],mXiPDG   );
			else                     vC1DBC.SetXYZM(0.,0.,0.,0.);
			if      (kaonBach2== 1)  vC2DBC.SetXYZM(MomCasc2DB[0],MomCasc2DB[1],MomCasc2DB[2],mOmegaPDG);
			else if (kaonBach2==-1)  vC2DBC.SetXYZM(MomCasc2DB[0],MomCasc2DB[1],MomCasc2DB[2],mXiPDG   );
			else                     vC2DBC.SetXYZM(0.,0.,0.,0.);

			vLm1    = vP1    + vN1;
			vLm2    = vP2    + vN2;
			vCs1    = vLm1   + vB1;
			vCs1L   = vL1    + vB1;
			vCs2    = vLm2   + vB2;
			vCs2L   = vL2    + vB2;
			vCsCs   = vCs1   + vCs2;
			vCsCsLL = vCs1L  + vCs2L;
			vCsCsCC = vC1    + vC2;
			vDBCC   = vC1DBC + vC2DBC;

			vC1DBL.SetXYZM(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2],vCs1L.M());
			vC2DBL.SetXYZM(MomCasc2DB[0],MomCasc2DB[1],MomCasc2DB[2],vCs2L.M());
			vC1DB.SetXYZM(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2],vCs1.M());
			vC2DB.SetXYZM(MomCasc2DB[0],MomCasc2DB[1],MomCasc2DB[2],vCs2.M());

			vDBLL   = vC1DBL + vC2DBL;
			vDB     = vC1DB  + vC2DB;


	//------------------------------------------------------------------------------------------
	// Other information
	//------------------------------------------------------------------------------------------

			Double_t CPACasc1PV, CPACasc2PV;
			Double_t IPcasc1, IPcasc2;
			Double_t DCAOmDa1, DCAOmDa2;
			Double_t CPAV01, CPAV02;
			Double_t DCAV0Da1, DCAV0Da2;
			CPACasc1PV = casc1->GetCascadeCosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2]);
			CPACasc2PV = casc2->GetCascadeCosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2]);
			IPcasc1  = casc1->GetDcascade(PosPV[0],PosPV[1],PosPV[2]);
			IPcasc2  = casc2->GetDcascade(PosPV[0],PosPV[1],PosPV[2]);
			DCAOmDa1 = casc1->GetDcaXiDaughters();
			DCAOmDa2 = casc2->GetDcaXiDaughters();
			CPAV01   = casc1->GetV0CosineOfPointingAngle(PosCasc1[0],PosCasc1[1],PosCasc1[2]);
			CPAV02   = casc2->GetV0CosineOfPointingAngle(PosCasc2[0],PosCasc2[1],PosCasc2[2]);
			DCAV0Da1 = casc1->GetDcaV0Daughters();
			DCAV0Da2 = casc2->GetDcaV0Daughters();

			Double_t AlphaV01 = 0.;
			Double_t PtArmV01 = 0.;
			TVector3 momentumPosV01(MomP1[0],MomP1[1],MomP1[2]);
			TVector3 momentumNegV01(MomN1[0],MomN1[1],MomN1[2]);
			TVector3 momentumTotV01(MomP1[0]+MomN1[0],MomP1[1]+MomN1[1],MomP1[2]+MomN1[2]);
			Double_t qPosV01 = momentumPosV01.Dot(momentumTotV01)/momentumTotV01.Mag();
			Double_t qNegV01 = momentumNegV01.Dot(momentumTotV01)/momentumTotV01.Mag();
			AlphaV01 = (qPosV01-qNegV01)/(qPosV01+qNegV01);
			PtArmV01 = momentumNegV01.Perp(momentumTotV01);

			Double_t AlphaV02 = 0.;
			Double_t PtArmV02 = 0.;
			TVector3 momentumPosV02(MomP2[0],MomP2[1],MomP2[2]);
			TVector3 momentumNegV02(MomN2[0],MomN2[1],MomN2[2]);
			TVector3 momentumTotV02(MomP2[0]+MomN2[0],MomP2[1]+MomN2[1],MomP2[2]+MomN2[2]);
			Double_t qPosV02 = momentumPosV02.Dot(momentumTotV02)/momentumTotV02.Mag();
			Double_t qNegV02 = momentumNegV02.Dot(momentumTotV02)/momentumTotV02.Mag();
			AlphaV02 = (qPosV02-qNegV02)/(qPosV02+qNegV02);
			PtArmV02 = momentumNegV02.Perp(momentumTotV02);

			Double_t AlphaDB = 0.;
			Double_t PtArmDB = 0.;
			TVector3 momentumDB1(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2]);
			TVector3 momentumDB2(MomCasc2DB[0],MomCasc2DB[1],MomCasc2DB[2]);
			TVector3 momentumTotDB(MomCasc1DB[0]+MomCasc2DB[0],MomCasc1DB[1]+MomCasc2DB[1],MomCasc1DB[2]+MomCasc2DB[2]);
			Double_t qDB1 = momentumDB1.Dot(momentumTotDB)/momentumTotDB.Mag();
			Double_t qDB2 = momentumDB2.Dot(momentumTotDB)/momentumTotDB.Mag();
			AlphaDB = (qDB1-qDB2)/(qDB1+qDB2);
			PtArmDB = momentumDB2.Perp(momentumTotDB);


	//------------------------------------------------------------------------------------------
	// Output
	//------------------------------------------------------------------------------------------

			fCandidateVariables[ 0] = vLm1.M();
			fCandidateVariables[ 1] = vLm2.M();
			fCandidateVariables[ 2] = vCs1.M();
			fCandidateVariables[ 3] = vCs1L.M();
			fCandidateVariables[ 4] = vCs2.M();
			fCandidateVariables[ 5] = vCs2L.M();
			fCandidateVariables[ 6] = vCsCs.M();
			fCandidateVariables[ 7] = vCsCsLL.M();
			fCandidateVariables[ 8] = vCsCsCC.M();
			fCandidateVariables[ 9] = vDB.M();
			fCandidateVariables[10] = vDBLL.M();
			fCandidateVariables[11] = vDBCC.M();
//			fCandidateVariables[12] = rCasc1;
//			fCandidateVariables[13] = rCasc2;
			fCandidateVariables[12] = ctauO1;
			fCandidateVariables[13] = ctauO2;
//			fCandidateVariables[16] = rV01;
//			fCandidateVariables[17] = rV02;
			fCandidateVariables[14] = ctauV01;
			fCandidateVariables[15] = ctauV02;
			fCandidateVariables[16] = CPACasc1;
			fCandidateVariables[17] = CPACasc2;
			fCandidateVariables[18] = IPcasc1;
			fCandidateVariables[19] = IPcasc2;
			fCandidateVariables[20] = DCAOmDa1;
			fCandidateVariables[21] = DCAOmDa2;
//			fCandidateVariables[22] = CPAV01;
//			fCandidateVariables[23] = CPAV02;
			fCandidateVariables[22] = DCAV0Da1;
			fCandidateVariables[23] = DCAV0Da2;
			fCandidateVariables[24] = DCADB;
			fCandidateVariables[25] = ctauDB;
			fCandidateVariables[26] = cpaDB;
			fCandidateVariables[27] = typeOfCasc1;
			fCandidateVariables[28] = typeOfCasc2;
//			fCandidateVariables[35] = MomP1P;
//			fCandidateVariables[36] = MomP2P;
//			fCandidateVariables[37] = MomN1P;
//			fCandidateVariables[38] = MomN2P;
//			fCandidateVariables[39] = MomB1P;
//			fCandidateVariables[40] = MomB2P;
			fCandidateVariables[29] = static_cast<Float_t>(triggerType);
//			fCandidateVariables[42] = AlphaV01;
//			fCandidateVariables[43] = PtArmV01;
//			fCandidateVariables[44] = AlphaV02;
//			fCandidateVariables[45] = PtArmV02;
//			fCandidateVariables[46] = AlphaDB;
//			fCandidateVariables[47] = PtArmDB;
			fCandidateVariables[30] = CPACasc1PV;
			fCandidateVariables[31] = CPACasc2PV;

			fVariablesTree->Fill();

		}//cascade loop 2

	}//cascade loop 1

//	printf("### Number of V0s (END)  : %d\n",nV01Test);

/*
	Double_t ratioLambda = static_cast<Double_t>(nAntiLambda)/static_cast<Double_t>(nLambda);
	Double_t ratioLambdaAll = static_cast<Double_t>(fCountAntiLambda)/static_cast<Double_t>(fCountLambda);
	printf("#######################################################################################################################\n");
	printf("#####               Lambda: Pro:%7d->%7d,  Anti:%7d->%7d  (ratio:%6.4f->%6.4f) ######################\n",nLambda,fCountLambda,nAntiLambda,fCountAntiLambda,ratioLambda,ratioLambdaAll);
	printf("#######################################################################################################################\n");
	printf("##### Analyzed events: %7d, Run Number: %10d, Centrality: %5d, nV0s: %5d ################################\n",fCountEvent,runNumber,fCentrality,nV0s);
	printf("#######################################################################################################################\n\n\n");
*/
	fCountEvent++;

}
//______________________________________________________________________________________________________
void AliAnalysisTaskOmegaOmegaOX::DefineTreeVariables() {
  //
  // This is to define tree variables
  //

  const char* nameoutput1 = GetOutputSlot(1)->GetContainer()->GetName();
  fParametersTree = new TTree(nameoutput1,"Parameters tree");
  Int_t nVar1 = 17;
  fParameters = new Float_t [nVar1];
  TString * fParameterNames = new TString[nVar1];

	fParameterNames[ 0]="RecoTypeDB";// = fRecoTypeDB;
	fParameterNames[ 1]="LikeSignDB";// = fLikeSignDB;
	fParameterNames[ 2]="TPCSigma";// = fReqSigmaTPC;
	fParameterNames[ 3]="TPCClusters";// = fReqClustersTPC;
	fParameterNames[ 4]="TOFSigma";// = fReqSigmaTOF;
	fParameterNames[ 5]="PseudoRap";// = fReqPseudoRap;
	fParameterNames[ 6]="CPADB";// = fCPADibaryon;
	fParameterNames[ 7]="WinCascade";// = fMassWinCascade;
  fParameterNames[ 8]="Chi2max";// = fCsReqChi2max;
  fParameterNames[ 9]="DV0min";// = fCsReqDV0min;
  fParameterNames[10]="WinLambda";// = fCsReqMassWinLambda;
  fParameterNames[11]="DBachMin";// = fCsReqDBachMin;
  fParameterNames[12]="DCAmax";// = fCsReqDCAmax; 
  fParameterNames[13]="CPAmin";// = fCsReqCPAmin;
  fParameterNames[14]="Rmin";// = fCsReqRmin;
  fParameterNames[15]="Rmax";// = fCsReqRmax;
	fParameterNames[16]="RecoSelfCasc";// = fRecoSelfCasc;

  for (Int_t ivar=0; ivar<nVar1; ivar++) {
    fParametersTree->Branch(fParameterNames[ivar].Data(),&fParameters[ivar],Form("%s/f",fParameterNames[ivar].Data()));
  }


  const char* nameoutput2 = GetOutputSlot(2)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput2,"Candidates variables tree");
  Int_t nVar2 = 32;
  fCandidateVariables = new Float_t [nVar2];
  TString * fCandidateVariableNames = new TString[nVar2];

	fCandidateVariableNames[ 0]="IML1";// = vLm1.M();
	fCandidateVariableNames[ 1]="IML2";// = vLm2.M();
	fCandidateVariableNames[ 2]="IMC1";// = vCs1.M();
	fCandidateVariableNames[ 3]="IMC1L";// = vCs1L.M();
	fCandidateVariableNames[ 4]="IMC2";// = vCs2.M();
	fCandidateVariableNames[ 5]="IMC2L";// = vCs2L.M();
	fCandidateVariableNames[ 6]="IMCC";// = vCsCs.M();
	fCandidateVariableNames[ 7]="IMCCL";// = vCsCsLL.M();
	fCandidateVariableNames[ 8]="IMCCC";// = vCsCsCC.M();
	fCandidateVariableNames[ 9]="IMDB";// = vDB.M();
	fCandidateVariableNames[10]="IMDBL";// = vDBLL.M();
	fCandidateVariableNames[11]="IMDBC";// = vDBCC.M();
//	fCandidateVariableNames[12]="rCs1";// = rCasc1;
//	fCandidateVariableNames[13]="rCs2";// = rCasc2;
	fCandidateVariableNames[12]="ctauO1";// = ctauO1;
	fCandidateVariableNames[13]="ctauO2";// = ctauO2;
//	fCandidateVariableNames[16]="rV01";// = rV01;
//	fCandidateVariableNames[17]="rV02";// = rV02;
	fCandidateVariableNames[14]="ctauL1";// = ctauV01;
	fCandidateVariableNames[15]="ctauL2";// = ctauV02;
	fCandidateVariableNames[16]="CPACs1";// = CPACasc1;
	fCandidateVariableNames[17]="CPACs2";// = CPACasc2;
	fCandidateVariableNames[18]="IPCs1";// = IPcasc1;
	fCandidateVariableNames[19]="IPCs2";// = IPcasc2;
	fCandidateVariableNames[20]="DCAOmDa1";// = DCAOmDa1;
	fCandidateVariableNames[21]="DCAOmDa2";// = DCAOmDa2;
//	fCandidateVariableNames[22]="CPAV01";// = CPAV01;
//	fCandidateVariableNames[23]="CPAV02";// = CPAV02;
	fCandidateVariableNames[22]="DCAV0Da1";// = DCAV0Da1;
	fCandidateVariableNames[23]="DCAV0Da2";// = DCAV0Da2;
	fCandidateVariableNames[24]="DCADB";// = DCADB;
	fCandidateVariableNames[25]="ctauDB";// = ctauDB;
	fCandidateVariableNames[26]="CPADB";// = cpaDB;
	fCandidateVariableNames[27]="t1";// = typeOfCasc1;
	fCandidateVariableNames[28]="t2";// = typeOfCasc2;
//	fCandidateVariableNames[35]="P1P";// = MomP1P;
//	fCandidateVariableNames[36]="P2P";// = MomP2P;
//	fCandidateVariableNames[37]="P1N";// = MomN1P;
//	fCandidateVariableNames[38]="P2N";// = MomN2P;
//	fCandidateVariableNames[39]="P1B";// = MomB1P;
//	fCandidateVariableNames[40]="P2B";// = MomB2P;
	fCandidateVariableNames[29]="TT";// = static_cast<Float_t>(triggerType);
//	fCandidateVariableNames[42]="AlphaV01";// = AlphaV01;
//	fCandidateVariableNames[43]="PtArmV01";// = PtArmV01;
//	fCandidateVariableNames[44]="AlphaV02";// = AlphaV02;
//	fCandidateVariableNames[45]="PtArmV02";// = PtArmV02;
//	fCandidateVariableNames[46]="AlphaDB";// = AlphaDB;
//	fCandidateVariableNames[47]="PtArmDB";// = PtArmDB;
	fCandidateVariableNames[30]="CPACs1PV";// = CPACasc1PV;
	fCandidateVariableNames[31]="CPACs2PV";// = CPACasc2PV;

  for (Int_t ivar=0; ivar<nVar2; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;

}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskOmegaOmegaOX::InvMassLambda(Double_t MomPos[3], Double_t MomNeg[3], AliESDtrack *pos, AliESDtrack *neg, Double_t v0Return[1]) {

	//------------------------------------------------------------------------------------------
	// version 1.20 (2016/03/02)
	// Input:  MomPos & MomNeg: Array of momentum (positive and negative tracks) (0:Px, 1:Py, 2:Pz)
	//         pos & neg: AliESDtrack information for each track that has MomPos/MomNeg momentum
	// Return: Invariant mass of Lambda (Proton(+/-) + Pion(-/+))
	//         v0Return[0]: 1:Lambda, -1:LambdaBar, 0:Not Lambda(Bar)
	//------------------------------------------------------------------------------------------

	// Set cut parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Double_t cutProtonMomTPC = fProtonPMax;
	Double_t cutPionMomTPC   = fPionPMax;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Prepare constant and parameter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Double_t mProtonPDG = TDatabasePDG::Instance()->GetParticle(2212)->Mass();//0.938272
	Double_t mPionPDG     = TDatabasePDG::Instance()->GetParticle(211)->Mass();//0.139570

	AliExternalTrackParam posEtp(*pos);
	AliExternalTrackParam negEtp(*neg);

	v0Return[0] = 0;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// TPC sigma cut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// For positive track
	Double_t posTPCProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos, AliPID::kProton  ));
	Double_t posTPCPion     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos, AliPID::kPion    ));
	if (posTPCProton>fReqSigmaTPC && posTPCPion>fReqSigmaTPC) return 0.;

	// For negative track
	Double_t negTPCProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kProton  ));
	Double_t negTPCPion     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kPion    ));
	if (negTPCProton>fReqSigmaTPC && negTPCPion>fReqSigmaTPC) return 0.;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// TOF sigma cut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// For positive track
	Double_t posTOFPion   = -999.;
	Double_t posTOFProton = -999.;
	Bool_t   posTOFOn     = kFALSE;
	ULong_t posStatus = pos->GetStatus();
	if(posStatus&AliESDtrack::kTOFpid) {
		posTOFProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pos, AliPID::kProton  ));
		posTOFPion     = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pos, AliPID::kPion    ));
//		if (posTOFProton>fReqSigmaTOF && posTOFPion>fReqSigmaTOF) return 0.;
		if (posTOFProton>fReqSigmaTOF) return 0.;
		posTOFOn       = kTRUE;
	}

	// For negative track
	Double_t negTOFPion   = -999.;
	Double_t negTOFProton = -999.;
	Bool_t   negTOFOn     = kFALSE;
	ULong_t negStatus = neg->GetStatus();
	if(negStatus&AliESDtrack::kTOFpid) {
		negTOFProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(neg, AliPID::kProton  ));
		negTOFPion     = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(neg, AliPID::kPion    ));
		if (negTOFProton>fReqSigmaTOF && negTOFPion>fReqSigmaTOF) return 0.;
		negTOFOn       = kTRUE;
	}
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Calculate invariant mass ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// PID for each track
	Bool_t posProton = kFALSE;
	Bool_t posPion   = kFALSE;
	Bool_t negProton = kFALSE;
	Bool_t negPion   = kFALSE;
	if (posTOFOn) {// PID with TOF information
		if (posTOFProton<fReqSigmaTOF) posProton = kTRUE;
		if (posTOFPion<fReqSigmaTOF)   posPion   = kTRUE;
	} else {// PID with TPC information
		if (posTPCProton<fReqSigmaTPC && posEtp.GetP()<cutProtonMomTPC) posProton = kTRUE;
		if (posTPCPion<fReqSigmaTPC   && posEtp.GetP()<cutPionMomTPC)   posPion   = kTRUE;
	}
	if (negTOFOn) {// PID with TOF information
		if (negTOFProton<fReqSigmaTOF) negProton = kTRUE;
		if (negTOFPion<fReqSigmaTOF)   negPion   = kTRUE;
	} else {// PID with TPC information
		if (negTPCProton<fReqSigmaTPC && negEtp.GetP()<cutProtonMomTPC) negProton = kTRUE;
		if (negTPCPion<fReqSigmaTPC   && negEtp.GetP()<cutPionMomTPC)   negPion   = kTRUE;
	}

	// PID combination
	Bool_t Lambda    = kFALSE;
	Bool_t LambdaBar = kFALSE;
	if (posProton==kTRUE && negPion==kTRUE) Lambda    = kTRUE;
	if (negProton==kTRUE && posPion==kTRUE) LambdaBar = kTRUE;

	// Avoid double count
	if (Lambda==kTRUE && LambdaBar==kTRUE) {
		Double_t sigmaLambda    = -999.;
		Double_t sigmaLambdaBar = -999.;
		if (posTOFOn) {
			sigmaLambda    = posTOFProton * posTOFProton;
			sigmaLambdaBar = posTOFPion   * posTOFPion;
		} else {
			sigmaLambda    = posTPCProton * posTPCProton;
			sigmaLambdaBar = posTPCPion   * posTPCPion;
		}
		if (negTOFOn) {
			sigmaLambda    = sigmaLambda    + negTOFPion   * negTOFPion;
			sigmaLambdaBar = sigmaLambdaBar + negTOFProton * negTOFProton;
		} else {
			sigmaLambda    = sigmaLambda    + negTPCPion   * negTPCPion;
			sigmaLambdaBar = sigmaLambdaBar + negTPCProton * negTPCProton;
		}
		sigmaLambda    = TMath::Sqrt(sigmaLambda);
		sigmaLambdaBar = TMath::Sqrt(sigmaLambdaBar);

		if (sigmaLambda<=sigmaLambdaBar) LambdaBar = kFALSE;
		else Lambda = kFALSE;
	}

	if (Lambda) v0Return[0] = 1.;
	if (LambdaBar) v0Return[0] = -1.;

	TLorentzVector vLambda,vPos,vNeg;
	if (Lambda) {
		vPos.SetXYZM(MomPos[0],MomPos[1],MomPos[2],mProtonPDG);
		vNeg.SetXYZM(MomNeg[0],MomNeg[1],MomNeg[2],mPionPDG);
	} else if (LambdaBar) {
		vPos.SetXYZM(MomPos[0],MomPos[1],MomPos[2],mPionPDG);
		vNeg.SetXYZM(MomNeg[0],MomNeg[1],MomNeg[2],mProtonPDG);
	} else {
		vPos.SetXYZM(0.,0.,0.,0.);
		vNeg.SetXYZM(0.,0.,0.,0.);
	}
	vLambda = vPos + vNeg;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//	printf("\n### TPC:  POS Pr:%10f, POS Pi:%10f,\tNEG Pr:%10f, NEG Pi:%10f\n",posTPCProton,posTPCPion,negTPCProton,negTPCPion);
//	printf("!! Bool:  POS Pr:%10d, POS Pi:%10d,\tNEG Pr:%10d, NEG Pi:%10d,\tL:%d, A:%d\n",posProton,posPion,negProton,negPion,Lambda,LambdaBar);
//	if (posTOFOn||negTOFOn) { 
//	  printf("### TOF:  POS Pr:%10f, POS Pi:%10f,\tNEG Pr:%10f, NEG Pi:%10f\n",posTOFProton,posTOFPion,negTOFProton,negTOFPion);
//	} 
//	printf("!!! THEREFORE: %f (1:Lambda, -1:Anti-Lambda, 0:Other)\n",v0Return[0]);


  return vLambda.M();

}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskOmegaOmegaOX::InvMassLambdaStar(Double_t MomPos[3], Double_t MomNeg[3], AliESDtrack *pos, AliESDtrack *neg, Double_t v0Return[1]) {

	//------------------------------------------------------------------------------------------
	// version 1.20 (2016/03/02)
	// Input:  MomPos & MomNeg: Array of momentum (positive and negative tracks) (0:Px, 1:Py, 2:Pz)
	//         pos & neg: AliESDtrack information for each track that has MomPos/MomNeg momentum
	// Return: Invariant mass of Lambda(1520) (Proton(+/-) + Kaon(-/+))
	//         v0Return[0]: 1:LambdaStar, -1:LambdaStarBar, 0:Not LambdaStar(Bar)
	//------------------------------------------------------------------------------------------

	// Set cut parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Double_t cutProtonMomTPC = fProtonPMax;
	Double_t cutKaonMomTPC   = fKaonPMax;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Prepare constant and parameter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Double_t mProtonPDG = TDatabasePDG::Instance()->GetParticle(2212)->Mass();//0.938272
	Double_t mKaonPDG   = TDatabasePDG::Instance()->GetParticle(321)->Mass();//0.493677

	AliExternalTrackParam posEtp(*pos);
	AliExternalTrackParam negEtp(*neg);

	v0Return[0] = 0;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// TPC sigma cut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// For positive track
	Double_t posTPCProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos, AliPID::kProton  ));
	Double_t posTPCKaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos, AliPID::kKaon    ));
	if (posTPCProton>fReqSigmaTPC && posTPCKaon>fReqSigmaTPC) return 0.;

	// For negative track
	Double_t negTPCProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kProton  ));
	Double_t negTPCKaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kKaon    ));
	if (negTPCProton>fReqSigmaTPC && negTPCKaon>fReqSigmaTPC) return 0.;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// TOF sigma cut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// For positive track
	Double_t posTOFKaon   = -999.;
	Double_t posTOFProton = -999.;
	Bool_t   posTOFOn     = kFALSE;
	ULong_t posStatus = pos->GetStatus();
	if(posStatus&AliESDtrack::kTOFpid) {
		posTOFProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pos, AliPID::kProton  ));
		posTOFKaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pos, AliPID::kKaon    ));
		if (posTOFProton>fReqSigmaTOF && posTOFKaon>fReqSigmaTOF) return 0.;
		posTOFOn       = kTRUE;
	}

	// For negative track
	Double_t negTOFKaon   = -999.;
	Double_t negTOFProton = -999.;
	Bool_t   negTOFOn     = kFALSE;
	ULong_t negStatus = neg->GetStatus();
	if(negStatus&AliESDtrack::kTOFpid) {
		negTOFProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(neg, AliPID::kProton  ));
		negTOFKaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(neg, AliPID::kKaon    ));
		if (negTOFProton>fReqSigmaTOF && negTOFKaon>fReqSigmaTOF) return 0.;
		negTOFOn       = kTRUE;
	}
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Calculate invariant mass ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// PID for each track
	Bool_t posProton = kFALSE;
	Bool_t posKaon   = kFALSE;
	Bool_t negProton = kFALSE;
	Bool_t negKaon   = kFALSE;
	if (posTOFOn) {// PID with TOF information
		if (posTOFProton<fReqSigmaTOF) posProton = kTRUE;
		if (posTOFKaon<fReqSigmaTOF)   posKaon   = kTRUE;
	} else {// PID with TPC information
		if (posTPCProton<fReqSigmaTPC && posEtp.GetP()<cutProtonMomTPC) posProton = kTRUE;
		if (posTPCKaon<fReqSigmaTPC   && posEtp.GetP()<cutKaonMomTPC)   posKaon   = kTRUE;
	}
	if (negTOFOn) {// PID with TOF information
		if (negTOFProton<fReqSigmaTOF) negProton = kTRUE;
		if (negTOFKaon<fReqSigmaTOF)   negKaon   = kTRUE;
	} else {// PID with TPC information
		if (negTPCProton<fReqSigmaTPC && negEtp.GetP()<cutProtonMomTPC) negProton = kTRUE;
		if (negTPCKaon<fReqSigmaTPC   && negEtp.GetP()<cutKaonMomTPC)   negKaon   = kTRUE;
	}

	// PID combination
	Bool_t LambdaStar    = kFALSE;
	Bool_t LambdaStarBar = kFALSE;
	if (posProton==kTRUE && negKaon==kTRUE) LambdaStar    = kTRUE;
	if (negProton==kTRUE && posKaon==kTRUE) LambdaStarBar = kTRUE;

	// Avoid double count
	if (LambdaStar==kTRUE && LambdaStarBar==kTRUE) {
		Double_t sigmaLambdaStar    = -999.;
		Double_t sigmaLambdaStarBar = -999.;
		if (posTOFOn) {
			sigmaLambdaStar    = posTOFProton * posTOFProton;
			sigmaLambdaStarBar = posTOFKaon   * posTOFKaon;
		} else {
			sigmaLambdaStar    = posTPCProton * posTPCProton;
			sigmaLambdaStarBar = posTPCKaon   * posTPCKaon;
		}
		if (negTOFOn) {
			sigmaLambdaStar    = sigmaLambdaStar    + negTOFKaon   * negTOFKaon;
			sigmaLambdaStarBar = sigmaLambdaStarBar + negTOFProton * negTOFProton;
		} else {
			sigmaLambdaStar    = sigmaLambdaStar    + negTPCKaon   * negTPCKaon;
			sigmaLambdaStarBar = sigmaLambdaStarBar + negTPCProton * negTPCProton;
		}
		sigmaLambdaStar    = TMath::Sqrt(sigmaLambdaStar);
		sigmaLambdaStarBar = TMath::Sqrt(sigmaLambdaStarBar);

		if (sigmaLambdaStar<=sigmaLambdaStarBar) LambdaStarBar = kFALSE;
		else LambdaStar = kFALSE;
	}

	if (LambdaStar) v0Return[0] = 1.;
	if (LambdaStarBar) v0Return[0] = -1.;

	TLorentzVector vLambdaStar,vPos,vNeg;
	if (LambdaStar) {
		vPos.SetXYZM(MomPos[0],MomPos[1],MomPos[2],mProtonPDG);
		vNeg.SetXYZM(MomNeg[0],MomNeg[1],MomNeg[2],mKaonPDG);
	} else if (LambdaStarBar) {
		vPos.SetXYZM(MomPos[0],MomPos[1],MomPos[2],mKaonPDG);
		vNeg.SetXYZM(MomNeg[0],MomNeg[1],MomNeg[2],mProtonPDG);
	} else {
		vPos.SetXYZM(0.,0.,0.,0.);
		vNeg.SetXYZM(0.,0.,0.,0.);
	}
	vLambdaStar = vPos + vNeg;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  return vLambdaStar.M();

}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskOmegaOmegaOX::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 2x2 determinant
  //--------------------------------------------------------------------
  return a00*a11 - a01*a10;
}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskOmegaOmegaOX::Det(Double_t a00,Double_t a01,Double_t a02,
         Double_t a10,Double_t a11,Double_t a12,
         Double_t a20,Double_t a21,Double_t a22) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 3x3 determinant
  //--------------------------------------------------------------------
  return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskOmegaOmegaOX::PropagateToDCA(AliESDv0 *v, AliExternalTrackParam *t, Double_t b) {
  //--------------------------------------------------------------------
  // This function returns the DCA between the V0 and the track
  //--------------------------------------------------------------------
  Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
  Double_t r[3]; t->GetXYZ(r);
  Double_t x1=r[0], y1=r[1], z1=r[2];
  Double_t p[3]; t->GetPxPyPz(p);
  Double_t px1=p[0], py1=p[1], pz1=p[2];

  Double_t x2,y2,z2;     // position and momentum of V0
  Double_t px2,py2,pz2;

  v->GetXYZ(x2,y2,z2);
  v->GetPxPyPz(px2,py2,pz2);

// calculation dca

  Double_t dd= Det(x2-x1,y2-y1,z2-z1,px1,py1,pz1,px2,py2,pz2);
  Double_t ax= Det(py1,pz1,py2,pz2);
  Double_t ay=-Det(px1,pz1,px2,pz2);
  Double_t az= Det(px1,py1,px2,py2);

  Double_t dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);
  if (dca > fCsReqDCAmax) return 1.e+33;

//points of the DCA
  Double_t t1 = Det(x2-x1,y2-y1,z2-z1,px2,py2,pz2,ax,ay,az)/
                Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);

  x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;

  if (x1*x1+y1*y1 > (fCsReqRmax+5.)*(fCsReqRmax+5.)) return 1.e+33;

  //propagate track to the points of DCA

  x1=x1*cs1 + y1*sn1;

  if (!t->PropagateTo(x1,b)) {
    AliError("Propagation failed");
    //    AliErrorF("Propagation failed for X=%f | V0: %f %f %f",x1,x2,y2,z2);
    //    t->Print();
    //
    return 1.e+33;
  }

  return dca;
}
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskOmegaOmegaOX::ReconstructCascadeString(Bool_t allocateStr, Int_t &nCascade, Int_t v0IDCasc[], Int_t trkIDCasc[]) {

	//------------------------------------------------------------------------------------------
	// version 1.30 (2016/04/13)
	// Reconstruction of Cascade by myself without cascade class (from AliCascadeVertexer.cxx)
	// Return: Number of Reconstructed cascades
	//------------------------------------------------------------------------------------------

  const Int_t nCascades = fESDEvent->GetNumberOfCascades();

//	Int_t nKill = 100;
	Int_t nKill = nCascades;

	// Get or calculate constant ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Double_t mElectronPDG = TDatabasePDG::Instance()->GetParticle(11)->Mass();//0.000511
	Double_t mPionPDG     = TDatabasePDG::Instance()->GetParticle(211)->Mass();//0.139570
	Double_t mKaonPDG     = TDatabasePDG::Instance()->GetParticle(321)->Mass();//0.493677
	Double_t mProtonPDG   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();//0.938272
	Double_t mLambdaPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();//1.115680
  Double_t mXiPDG       = TDatabasePDG::Instance()->GetParticle(3312)->Mass();//1.321710
  Double_t mOmegaPDG    = TDatabasePDG::Instance()->GetParticle(3334)->Mass();//1.672450

  const Int_t nTracks = fESDEvent->GetNumberOfTracks();
  if (nTracks==0) {
    return kFALSE;
  }

  const Int_t nV0s = fESDEvent->GetNumberOfV0s();
  if (nV0s==0) {
    return kFALSE;
  }

  Double_t fCsReqRmin2 = fCsReqRmin*fCsReqRmin;
  Double_t fCsReqRmax2 = fCsReqRmax*fCsReqRmax;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Get or calculate constant for this event ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const AliESDVertex *esdTrkPV = (AliESDVertex*)fESDEvent->GetPrimaryVertexTracks();
	Double_t PosTrkPV[3];
	esdTrkPV->GetXYZ(PosTrkPV);

	const AliESDVertex *esdPV = (AliESDVertex*)fESDEvent->GetPrimaryVertex();
	Double_t PosPV[3];
	esdPV->GetXYZ(PosPV);

	Double_t vPVX = PosPV[0];
	Double_t vPVY = PosPV[1];
	Double_t vPVZ = PosPV[2];
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//------------------------------------------------------------------------------------------
	// Stores relevant V0s and tracks in array
	//------------------------------------------------------------------------------------------

	Int_t nCasc = fESDEvent->GetNumberOfCascades();
//	printf("\n\n### Stored nCascade:%d, nV0s:%d, nTracks:%d\n",nCasc,nV0s,nTracks);
	if ( fRecoSelfCasc == 1 ) {
		fESDEvent->ResetCascades();
		nCasc = fESDEvent->GetNumberOfCascades();
//		printf("### Clear nCascade:%d\n",nCasc);
	}

	// V0s
	Int_t v0Array[nV0s];
	Int_t nV0Survived=0;
	for (Int_t iV0=0; iV0<nV0s; iV0++) {
		AliESDv0 *v0rel = fESDEvent->GetV0(iV0);
		if (v0rel->GetOnFlyStatus()) continue;
		if (v0rel->GetD(PosPV[0],PosPV[1],PosPV[2])<fCsReqDV0min) continue;
		AliESDtrack *trkPrel = fESDEvent->GetTrack(v0rel->GetPindex());
		AliESDtrack *trkNrel = fESDEvent->GetTrack(v0rel->GetNindex());
		if( !fESDtrackCuts->AcceptTrack(trkPrel) ) continue;
		if( !fESDtrackCuts->AcceptTrack(trkNrel) ) continue;
		Double_t trkPrelTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkPrel, AliPID::kProton));
		Double_t trkPrelTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkPrel, AliPID::kPion  ));
		if ( trkPrelTPCProton>fReqSigmaTPC && trkPrelTPCPion>fReqSigmaTPC ) continue;
		Double_t trkNrelTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkNrel, AliPID::kProton));
		Double_t trkNrelTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkNrel, AliPID::kPion  ));
		if ( trkNrelTPCProton>fReqSigmaTPC && trkNrelTPCPion>fReqSigmaTPC ) continue;
		v0Array[nV0Survived++] = iV0;
	}

	// tracks
	Int_t trkArray[nTracks];
	Int_t nTrkSurvived=0;
	for (Int_t iTrk=0; iTrk<nTracks; iTrk++) {
		AliESDtrack *trkrel = fESDEvent->GetTrack(iTrk);
		ULong_t status = trkrel->GetStatus();
		if (status&AliESDtrack::kITSpureSA) continue;
		if ((status&AliESDtrack::kITSrefit)==0)
			if ((status&AliESDtrack::kTPCrefit)==0) continue;
		if (TMath::Abs(trkrel->GetD(PosPV[0],PosPV[1],fBzkG))<fCsReqDBachMin) continue;
		if( !fESDtrackCuts->AcceptTrack(trkrel) ) continue;
		Double_t trkrelTPCKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkrel, AliPID::kKaon));
		Double_t trkrelTPCPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkrel, AliPID::kPion));
		if ( trkrelTPCKaon>fReqSigmaTPC && trkrelTPCPion>fReqSigmaTPC ) continue;
		trkArray[nTrkSurvived++] = iTrk;
	}


	//------------------------------------------------------------------------------------------
	// Reconstruction of cascades
	//------------------------------------------------------------------------------------------

	Int_t nCascSurvived = 0;
	Int_t nCascExist = 0;

	// Candidate of cascade
	if (fLikeSignDB==0||fLikeSignDB==1) {//including Omega or Xi
		for (Int_t iV0=0; iV0<nV0Survived; iV0++) {//V0 loop
			Int_t indexV0 = v0Array[iV0];
			AliESDv0 *v0point = fESDEvent->GetV0(indexV0);
			AliESDv0 v0Casc(*v0point);
			Int_t trkPind = v0Casc.GetPindex();
			Int_t trkNind = v0Casc.GetNindex();
			v0Casc.ChangeMassHypothesis(kLambda0);
			if (TMath::Abs(v0Casc.GetEffMass()-mLambdaPDG)>fCsReqMassWinLambda) continue;

			for (Int_t iTrk=0; iTrk<nTrkSurvived; iTrk++) {//bachelor loop
				Int_t indexTrk = trkArray[iTrk];
				if (indexTrk==trkPind) continue;
				if (indexTrk==trkNind) continue;

				AliESDtrack *trkCasc = fESDEvent->GetTrack(indexTrk);
				if (trkCasc->GetSign()>0) continue;

				AliESDv0 *v0Propagated = &v0Casc;
				AliExternalTrackParam etpTrkCasc(*trkCasc);
				AliExternalTrackParam *trkPropagated = &etpTrkCasc;

				Double_t dcaCasc = PropagateToDCA(v0Propagated,trkPropagated,fBzkG);
				if (dcaCasc > fCsReqDCAmax) continue;

				AliESDcascade cascade(*v0Propagated,*trkPropagated,indexTrk);

				Double_t x,y,z;
				cascade.GetXYZcascade(x,y,z);
				Double_t r2 = x*x + y*y;
				if (r2 > fCsReqRmax2) continue;
				if (r2 < fCsReqRmin2) continue;

				Double_t pxV0,pyV0,pzV0;
				v0Propagated->GetPxPyPz(pxV0,pyV0,pzV0);
				if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

				Double_t x1,y1,z1;
				v0Propagated->GetXYZ(x1,y1,z1);
				if (r2 > (x1*x1+y1*y1)) continue;

				Double_t CPACascPV = cascade.GetCascadeCosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2]);
				if (CPACascPV<fCsReqCPAmin) continue; //condition on the cascade pointing angle 

				Double_t IPcasc  = cascade.GetDcascade(PosPV[0],PosPV[1],PosPV[2]);
				Double_t DCAOmDa = cascade.GetDcaXiDaughters();
				Double_t CPAV0   = cascade.GetV0CosineOfPointingAngle(x,y,z);
				Double_t DCAV0Da = cascade.GetDcaV0Daughters();

				Double_t pxB,pyB,pzB;
				cascade.GetBPxPyPz(pxB,pyB,pzB);
				TLorentzVector vCascade,vLambda,vBachelor;
				vLambda.SetXYZM(pxV0,pyV0,pzV0,mLambdaPDG);
				vBachelor.SetXYZM(pxB,pyB,pzB,mKaonPDG);
				vCascade = vLambda + vBachelor;
				Double_t massOmega = vCascade.M();
				Double_t massDiffOmega = TMath::Abs(massOmega-mOmegaPDG);
				vBachelor.SetXYZM(pxB,pyB,pzB,mPionPDG);
				vCascade = vLambda + vBachelor;
				Double_t massXi = vCascade.M();
				Double_t massDiffXi = TMath::Abs(massXi-mXiPDG);
				if (massDiffOmega>fMassWinCascade&&massDiffXi>fMassWinCascade) continue;
				if (fRecoTypeDB==1&&massDiffOmega>fMassWinCascade) continue;
				if (fRecoTypeDB==4&&massDiffXi>fMassWinCascade) continue;

//				if (nCascSurvived>=nKill) continue;
				if (allocateStr) {
					v0IDCasc[nCascSurvived]  = indexV0;
					trkIDCasc[nCascSurvived] = indexTrk;
				}

				nCascSurvived++;
//				printf("%d ",nCascSurvived);

			}//bachelor loop

		}//V0 loop

	}//including Omega or Xi

	// Candidate of anti-cascade
	if (fLikeSignDB==0||fLikeSignDB==2) {//including Omegabar or Xibar
		for (Int_t iV0=0; iV0<nV0Survived; iV0++) {//V0 loop
			Int_t indexV0 = v0Array[iV0];
			AliESDv0 *v0point = fESDEvent->GetV0(indexV0);
			AliESDv0 v0Casc(*v0point); 
			Int_t trkPind = v0Casc.GetPindex();
			Int_t trkNind = v0Casc.GetNindex();
			v0Casc.ChangeMassHypothesis(kLambda0Bar); 
			if (TMath::Abs(v0Casc.GetEffMass()-mLambdaPDG)>fCsReqMassWinLambda) continue;

			for (Int_t iTrk=0; iTrk<nTrkSurvived; iTrk++) {//bachelor loop
				Int_t indexTrk = trkArray[iTrk];
				if (indexTrk==trkPind) continue;
				if (indexTrk==trkNind) continue;

				AliESDtrack *trkCasc = fESDEvent->GetTrack(indexTrk);
				if (trkCasc->GetSign()<0) continue;

				AliESDv0 *v0Propagated = &v0Casc;
				AliExternalTrackParam etpTrkCasc(*trkCasc);
				AliExternalTrackParam *trkPropagated = &etpTrkCasc;

				Double_t dcaCasc = PropagateToDCA(v0Propagated,trkPropagated,fBzkG);
				if (dcaCasc > fCsReqDCAmax) continue;

				AliESDcascade cascade(*v0Propagated,*trkPropagated,indexTrk);

				Double_t x,y,z;
				cascade.GetXYZcascade(x,y,z);
				Double_t r2 = x*x + y*y;
				if (r2 > fCsReqRmax2) continue;
				if (r2 < fCsReqRmin2) continue;

				Double_t pxV0,pyV0,pzV0;
				v0Propagated->GetPxPyPz(pxV0,pyV0,pzV0);
				if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

				Double_t x1,y1,z1;
				v0Propagated->GetXYZ(x1,y1,z1);
				if (r2 > (x1*x1+y1*y1)) continue;

				Double_t CPACascPV = cascade.GetCascadeCosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2]);
				if (CPACascPV<fCsReqCPAmin) continue; //condition on the cascade pointing angle 

				Double_t IPcasc  = cascade.GetDcascade(PosPV[0],PosPV[1],PosPV[2]);
				Double_t DCAOmDa = cascade.GetDcaXiDaughters();
				Double_t CPAV0   = cascade.GetV0CosineOfPointingAngle(x,y,z);
				Double_t DCAV0Da = cascade.GetDcaV0Daughters();

				Double_t pxB,pyB,pzB;
				cascade.GetBPxPyPz(pxB,pyB,pzB);
				TLorentzVector vCascade,vLambda,vBachelor;
				vLambda.SetXYZM(pxV0,pyV0,pzV0,mLambdaPDG);
				vBachelor.SetXYZM(pxB,pyB,pzB,mKaonPDG);
				vCascade = vLambda + vBachelor;
				Double_t massOmega = vCascade.M();
				Double_t massDiffOmega = TMath::Abs(massOmega-mOmegaPDG);
				vBachelor.SetXYZM(pxB,pyB,pzB,mPionPDG);
				vCascade = vLambda + vBachelor;
				Double_t massXi = vCascade.M();
				Double_t massDiffXi = TMath::Abs(massXi-mXiPDG);
				if (massDiffOmega>fMassWinCascade&&massDiffXi>fMassWinCascade) continue;
				if (fRecoTypeDB==1&&massDiffOmega>fMassWinCascade) continue;
				if (fRecoTypeDB==4&&massDiffXi>fMassWinCascade) continue;

//				if (nCascSurvived>=nKill) continue;
				if (allocateStr) {
					v0IDCasc[nCascSurvived]  = indexV0;
					trkIDCasc[nCascSurvived] = indexTrk;
				}

				nCascSurvived++;
//				printf("%d ",nCascSurvived);

			}//bachelor loop

		}//V0 loop

	}//including Omegabar or Xibar

//	printf("### looped nCascade:%d\n",nCascSurvived);

	nCascade = nCascSurvived;

	return kTRUE;

}
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskOmegaOmegaOX::ReconstructCascade() {

	//------------------------------------------------------------------------------------------
	// version 1.00 (2016/04/01)
	// Reconstruction of Cascade by myself (from AliCascadeVertexer.cxx)
	// Return: Number of Reconstructed cascades
	//------------------------------------------------------------------------------------------

	// Get or calculate constant ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Double_t mElectronPDG = TDatabasePDG::Instance()->GetParticle(11)->Mass();//0.000511
	Double_t mPionPDG     = TDatabasePDG::Instance()->GetParticle(211)->Mass();//0.139570
	Double_t mKaonPDG     = TDatabasePDG::Instance()->GetParticle(321)->Mass();//0.493677
	Double_t mProtonPDG   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();//0.938272
	Double_t mLambdaPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();//1.115680
	Double_t mXiPDG       = TDatabasePDG::Instance()->GetParticle(3312)->Mass();//1.321710
	Double_t mOmegaPDG    = TDatabasePDG::Instance()->GetParticle(3334)->Mass();//1.672450

	const Int_t nTracks = fESDEvent->GetNumberOfTracks();
	if (nTracks==0) {
		return kFALSE;
	}

	const Int_t nV0s = fESDEvent->GetNumberOfV0s();
	if (nV0s==0) {
		return kFALSE;
	}

	Double_t fCsReqRmin2 = fCsReqRmin*fCsReqRmin;
	Double_t fCsReqRmax2 = fCsReqRmax*fCsReqRmax;

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Get or calculate constant for this event ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const AliESDVertex *esdTrkPV = (AliESDVertex*)fESDEvent->GetPrimaryVertexTracks();
	Double_t PosTrkPV[3];
	esdTrkPV->GetXYZ(PosTrkPV);

	const AliESDVertex *esdPV = (AliESDVertex*)fESDEvent->GetPrimaryVertex();
	Double_t PosPV[3];
	esdPV->GetXYZ(PosPV);

	Double_t vPVX = PosPV[0];
	Double_t vPVY = PosPV[1];
	Double_t vPVZ = PosPV[2];
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//------------------------------------------------------------------------------------------
	// Stores relevant V0s and tracks in array
	//------------------------------------------------------------------------------------------

	Int_t nCasc = fESDEvent->GetNumberOfCascades();
//	printf("\n\n### Stored nCascade:%d, nV0s:%d, nTracks:%d\n",nCasc,nV0s,nTracks);
	fESDEvent->ResetCascades();
	nCasc = fESDEvent->GetNumberOfCascades();
//	printf("### Clear nCascade:%d\n",nCasc);

	// V0s
	Int_t v0Array[nV0s];
	Int_t nV0Survived=0;
	for (Int_t iV0=0; iV0<nV0s; iV0++) {
		AliESDv0 *v0rel = fESDEvent->GetV0(iV0);
		if (v0rel->GetOnFlyStatus()) continue;
		if (v0rel->GetD(PosPV[0],PosPV[1],PosPV[2])<fCsReqDV0min) continue;
		AliESDtrack *trkPrel = fESDEvent->GetTrack(v0rel->GetPindex());
		AliESDtrack *trkNrel = fESDEvent->GetTrack(v0rel->GetNindex());
		if( !fESDtrackCuts->AcceptTrack(trkPrel) ) continue;
		if( !fESDtrackCuts->AcceptTrack(trkNrel) ) continue;
		Double_t trkPrelTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkPrel, AliPID::kProton));
		Double_t trkPrelTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkPrel, AliPID::kPion  ));
		if ( trkPrelTPCProton>fReqSigmaTPC && trkPrelTPCPion>fReqSigmaTPC ) continue;
		Double_t trkNrelTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkNrel, AliPID::kProton));
		Double_t trkNrelTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkNrel, AliPID::kPion  ));
		if ( trkNrelTPCProton>fReqSigmaTPC && trkNrelTPCPion>fReqSigmaTPC ) continue;
		v0Array[nV0Survived++] = iV0;
	}

	// tracks
	Int_t trkArray[nTracks];
	Int_t nTrkSurvived=0;
	for (Int_t iTrk=0; iTrk<nTracks; iTrk++) {
		AliESDtrack *trkrel = fESDEvent->GetTrack(iTrk);
		ULong_t status = trkrel->GetStatus();
		if (status&AliESDtrack::kITSpureSA) continue;
		if ((status&AliESDtrack::kITSrefit)==0)
			if ((status&AliESDtrack::kTPCrefit)==0) continue;
		if (TMath::Abs(trkrel->GetD(PosPV[0],PosPV[1],fBzkG))<fCsReqDBachMin) continue;
		if( !fESDtrackCuts->AcceptTrack(trkrel) ) continue;
		Double_t trkrelTPCKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkrel, AliPID::kKaon));
		Double_t trkrelTPCPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkrel, AliPID::kPion));
		if ( trkrelTPCKaon>fReqSigmaTPC && trkrelTPCPion>fReqSigmaTPC ) continue;
		trkArray[nTrkSurvived++] = iTrk;
	}



	//------------------------------------------------------------------------------------------
	// Reconstruction of cascades
	//------------------------------------------------------------------------------------------

	Int_t nCascSurvived = 0;
	Int_t nCascExist = 0;

  // Candidate of cascade
	if (fLikeSignDB==0||fLikeSignDB==1) {//including Omega or Xi
		for (Int_t iV0=0; iV0<nV0Survived; iV0++) {//V0 loop
			Int_t indexV0 = v0Array[iV0];
			AliESDv0 *v0point = fESDEvent->GetV0(indexV0);
			AliESDv0 v0Casc(*v0point);
			v0Casc.ChangeMassHypothesis(kLambda0);
			if (TMath::Abs(v0Casc.GetEffMass()-mLambdaPDG)>fCsReqMassWinLambda) continue;

			for (Int_t iTrk=0; iTrk<nTrkSurvived; iTrk++) {//bachelor loop
				Int_t indexTrk = trkArray[iTrk];
				if (indexTrk==v0Casc.GetIndex(0)) continue;
				if (indexTrk==v0Casc.GetIndex(1)) continue;

				AliESDtrack *trkCasc = fESDEvent->GetTrack(indexTrk);
				if (trkCasc->GetSign()>0) continue;

				AliESDv0 *v0Propagated = &v0Casc;
				AliExternalTrackParam etpTrkCasc(*trkCasc);
				AliExternalTrackParam *trkPropagated = &etpTrkCasc;

				Double_t dcaCasc = PropagateToDCA(v0Propagated,trkPropagated,fBzkG);
				if (dcaCasc > fCsReqDCAmax) continue;

				AliESDcascade cascade(*v0Propagated,*trkPropagated,indexTrk);

				Double_t x,y,z;
				cascade.GetXYZcascade(x,y,z);
				Double_t r2 = x*x + y*y;
				if (r2 > fCsReqRmax2) continue;
				if (r2 < fCsReqRmin2) continue;

				Double_t pxV0,pyV0,pzV0;
				v0Propagated->GetPxPyPz(pxV0,pyV0,pzV0);
				if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

				Double_t x1,y1,z1;
				v0Propagated->GetXYZ(x1,y1,z1);
				if (r2 > (x1*x1+y1*y1)) continue;

				if (cascade.GetCascadeCosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2])<fCsReqCPAmin) continue; //condition on the cascade pointing angle 

				Double_t pxB,pyB,pzB;
				cascade.GetBPxPyPz(pxB,pyB,pzB);
				TLorentzVector vCascade,vLambda,vBachelor;
				vLambda.SetXYZM(pxV0,pyV0,pzV0,mLambdaPDG);
				vBachelor.SetXYZM(pxB,pyB,pzB,mKaonPDG);
				vCascade = vLambda + vBachelor;
				Double_t massOmega = vCascade.M();
				vBachelor.SetXYZM(pxB,pyB,pzB,mPionPDG);
				vCascade = vLambda + vBachelor;
				Double_t massXi = vCascade.M();
				if (TMath::Abs(massOmega-mOmegaPDG)>fMassWinCascade&&
				    TMath::Abs(massXi-mXiPDG)>fMassWinCascade        ) continue;
				if (fRecoTypeDB==1&&TMath::Abs(massOmega-mOmegaPDG)>fMassWinCascade) continue;
				if (fRecoTypeDB==4&&TMath::Abs(massXi-mXiPDG)>fMassWinCascade) continue;

				cascade.SetDcaXiDaughters(dcaCasc);
				fESDEvent->AddCascade(&cascade);
				nCascSurvived++;

			}//bachelor loop

		}//V0 loop

	}//including Omega or Xi

	// Candidate of anti-cascade
	if (fLikeSignDB==0||fLikeSignDB==2) {//including Omegabar or Xibar
		for (Int_t iV0=0; iV0<nV0Survived; iV0++) {//V0 loop
			Int_t indexV0 = v0Array[iV0];
			AliESDv0 *v0point = fESDEvent->GetV0(indexV0);
			AliESDv0 v0Casc(*v0point);
			v0Casc.ChangeMassHypothesis(kLambda0Bar);
			if (TMath::Abs(v0Casc.GetEffMass()-mLambdaPDG)>fCsReqMassWinLambda) continue;

			for (Int_t iTrk=0; iTrk<nTrkSurvived; iTrk++) {//bachelor loop
				Int_t indexTrk = trkArray[iTrk];
				if (indexTrk==v0Casc.GetIndex(0)) continue;
				if (indexTrk==v0Casc.GetIndex(1)) continue;

				AliESDtrack *trkCasc = fESDEvent->GetTrack(indexTrk);
				if (trkCasc->GetSign()<0) continue;

				AliESDv0 *v0Propagated = &v0Casc;
				AliExternalTrackParam etpTrkCasc(*trkCasc);
				AliExternalTrackParam *trkPropagated = &etpTrkCasc;

				Double_t dcaCasc = PropagateToDCA(v0Propagated,trkPropagated,fBzkG);
				if (dcaCasc > fCsReqDCAmax) continue;

				AliESDcascade cascade(*v0Propagated,*trkPropagated,indexTrk);

				Double_t x,y,z;
				cascade.GetXYZcascade(x,y,z);
				Double_t r2 = x*x + y*y;
				if (r2 > fCsReqRmax2) continue;
				if (r2 < fCsReqRmin2) continue;

				Double_t pxV0,pyV0,pzV0;
				v0Propagated->GetPxPyPz(pxV0,pyV0,pzV0);
				if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

				Double_t x1,y1,z1;
				v0Propagated->GetXYZ(x1,y1,z1);
				if (r2 > (x1*x1+y1*y1)) continue;

				if (cascade.GetCascadeCosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2])<fCsReqCPAmin) continue; //condition on the cascade pointing angle 

				Double_t pxB,pyB,pzB;
				cascade.GetBPxPyPz(pxB,pyB,pzB);
				TLorentzVector vCascade,vLambda,vBachelor;
				vLambda.SetXYZM(pxV0,pyV0,pzV0,mLambdaPDG);
				vBachelor.SetXYZM(pxB,pyB,pzB,mKaonPDG);
				vCascade = vLambda + vBachelor;
				Double_t massOmega = vCascade.M();
				vBachelor.SetXYZM(pxB,pyB,pzB,mPionPDG);
				vCascade = vLambda + vBachelor;
				Double_t massXi = vCascade.M();
				if (TMath::Abs(massOmega-mOmegaPDG)>fMassWinCascade&&
				    TMath::Abs(massXi-mXiPDG)>fMassWinCascade        ) continue;
				if (fRecoTypeDB==1&&TMath::Abs(massOmega-mOmegaPDG)>fMassWinCascade) continue;
				if (fRecoTypeDB==4&&TMath::Abs(massXi-mXiPDG)>fMassWinCascade) continue;

				cascade.SetDcaXiDaughters(dcaCasc);
				fESDEvent->AddCascade(&cascade);
				nCascSurvived++;

			}//bachelor loop

		}//V0 loop

	}//including Omegabar or Xibar

//  printf("### looped nCascade:%d\n",nCascSurvived);

	return kTRUE;

}
//______________________________________________________________________________________________________
//__________________________________________________________________________
