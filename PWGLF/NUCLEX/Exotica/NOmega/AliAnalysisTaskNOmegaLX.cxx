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
#include "AliAnalysisTaskNOmegaLX.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNOmegaLX)

//______________________________________________________________________________________________________
AliAnalysisTaskNOmegaLX::AliAnalysisTaskNOmegaLX() : 
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
	fRecoTypeDB(0),// Reconstruction type of DiBaryon (0:All, 1:OmL)
	fLikeSignDB(1),// Like-sign of DB (0:ALL, 1:XL, 2:(Xbar)(Lbar), 3:X(Lbar), 4:(Xbar)L)
	fRecoSelfCasc(0),// Cascade reconstruction is made by (0:ESD class, 1:by myself)
	fReqSigmaTPC(3.0),
	fReqClustersTPC(80),
	fReqSigmaTOF(3.0),
	fReqPseudoRap(0.9),
	fProtonPMax(999.),// Max momentum of proton
	fPionPMax(1.5),//!!CHECK ALSO ESDtrackCuts!! Max momentum of pion
	fKaonPMax(0.6),//!!CHECK ALSO ESDtrackCuts!! Max momentum of kaon
	fTrackPMin(0.0),//Min momentum of track
	fCPADibaryon(0.9),//0.99875
	fDCADibaryon(2.0),
	fMassWinCascade(999.),
	fCsChi2max(33.),
	fCsDV0min(0.01),
	fCsMassWinLambda(0.0045),
	fCsDBachMin(0.01),
	fCsDCAmax(2.0), 
	fCsCPAmin(0.0),//0.98
	fCsRmin(0.2),
	fCsRmax(100.),
	fFiducialVolMin2(2.0),// Min radius of the fiducial volume
	fFiducialVolMax2(200.),// Max radius of the fiducial volume
	fPosDCAToPVMin2(0.0),//2.0 Min length of impact parameter for the positive track
	fNegDCAToPVMin2(0.0),//2.0 Min length of impact parameter for the negative track
	fDCADaughterMax2(1.0),// Max DCA between the daughter tracks
	fCPAMin2(-1.),// Min cosine of V0's pointing angle to PV
	fDCAToPVMin2(0.0),// Min DCA V0 to PV
	fCOADaughterMin2(0.0),// Min cosine between the daughter tracks
	fDCAZDaughterMax2(2.0),// Max DCAZ V0 to PV
	fWindowV02(0.015),// Mass window cut for Lambda
	fCentrality(0),
	fCountEvent(0),
	fIsDCA(0),
	fNAll(0)
{
}
//______________________________________________________________________________________________________
AliAnalysisTaskNOmegaLX::AliAnalysisTaskNOmegaLX(const Char_t* name) :
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
	fRecoTypeDB(0),// Reconstruction type of DiBaryon (0:All, 1:OmL)
	fLikeSignDB(1),// Like-sign of DB (0:ALL, 1:XL, 2:(Xbar)(Lbar), 3:X(Lbar), 4:(Xbar)L)
	fRecoSelfCasc(0),// Cascade reconstruction is made by (0:ESD class, 1:by myself)
	fReqSigmaTPC(3.0),
	fReqClustersTPC(80),
	fReqSigmaTOF(3.0),
	fReqPseudoRap(0.9),
	fProtonPMax(999.),// Max momentum of proton
	fPionPMax(1.5),//!!CHECK ALSO ESDtrackCuts!! Max momentum of pion
	fKaonPMax(0.6),//!!CHECK ALSO ESDtrackCuts!! Max momentum of kaon
	fTrackPMin(0.0),//Min momentum of track
	fCPADibaryon(0.9),//0.99875
	fDCADibaryon(2.0),
	fMassWinCascade(999.),
	fCsChi2max(33.),
	fCsDV0min(0.01),
	fCsMassWinLambda(0.0045),
	fCsDBachMin(0.01),
	fCsDCAmax(2.0), 
	fCsCPAmin(0.0),//0.98
	fCsRmin(0.2),
	fCsRmax(100.),
	fFiducialVolMin2(2.0),// Min radius of the fiducial volume
	fFiducialVolMax2(200.),// Max radius of the fiducial volume
	fPosDCAToPVMin2(0.0),//2.0 Min length of impact parameter for the positive track
	fNegDCAToPVMin2(0.0),//2.0 Min length of impact parameter for the negative track
	fDCADaughterMax2(1.0),// Max DCA between the daughter tracks
	fCPAMin2(-1.),// Min cosine of V0's pointing angle to PV
	fDCAToPVMin2(0.0),// Min DCA V0 to PV
	fCOADaughterMin2(0.0),// Min cosine between the daughter tracks
	fDCAZDaughterMax2(2.0),// Max DCAZ V0 to PV
	fWindowV02(0.015),// Mass window cut for Lambda
	fCentrality(0),
	fCountEvent(0),
	fIsDCA(0),
	fNAll(0)
{
	//
	// Constructor. Initialization of Inputs and Outputs
	//
	Info("AliAnalysisTaskNOmegaLX","Calling Constructor");

	DefineOutput(1,TTree::Class());  //My private output
	DefineOutput(2,TTree::Class());  //My private output


	//V0 cuts
	fESDtrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
	fESDtrackCutsV0->SetMinNClustersTPC(80);
	fESDtrackCutsV0->SetAcceptKinkDaughters(kFALSE);
	fESDtrackCutsV0->SetMaxChi2PerClusterTPC(5);
	fESDtrackCutsV0->SetRequireTPCRefit(kTRUE);
	fESDtrackCutsV0->SetEtaRange(-0.9,0.9);
	fESDtrackCutsV0->SetPtRange(fTrackPMin,1.5);//pi
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
AliAnalysisTaskNOmegaLX::~AliAnalysisTaskNOmegaLX() {
	//
	// destructor
	//
	Info("~AliAnalysisTaskNOmegaLX","Calling Destructor");
  
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
void AliAnalysisTaskNOmegaLX::Init() {
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
void AliAnalysisTaskNOmegaLX::UserExec(Option_t *)
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
void AliAnalysisTaskNOmegaLX::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //AliInfo("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  return;
}

//______________________________________________________________________________________________________
void AliAnalysisTaskNOmegaLX::UserCreateOutputObjects() 
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
void AliAnalysisTaskNOmegaLX::MakeAnalysis(TClonesArray *mcArray,AliESDEvent *fESDEvent)
{

  //------------------------------------------------------------------------------------------
  // version NO2-0-2 (2016/03/03)
  // Reconstruct 1 cascade and 1 V0 (by myself)
  // and calculate the invariant mass of N-Omega
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
  fParameters[ 8] = fCsChi2max;
  fParameters[ 9] = fCsDV0min;
  fParameters[10] = fCsMassWinLambda;
  fParameters[11] = fCsDBachMin;
  fParameters[12] = fCsDCAmax; 
  fParameters[13] = fCsCPAmin;
  fParameters[14] = fCsRmin;
  fParameters[15] = fCsRmax;
	fParameters[16] = fRecoSelfCasc;
	fParameters[17] = fProtonPMax;
	fParameters[18] = fPionPMax;
	fParameters[19] = fKaonPMax;
	fParameters[20] = fDCADibaryon;
	fParameters[21] = fFiducialVolMin2;
	fParameters[22] = fFiducialVolMax2;
	fParameters[23] = fPosDCAToPVMin2;
	fParameters[24] = fNegDCAToPVMin2;
	fParameters[25] = fDCADaughterMax2;
	fParameters[26] = fCPAMin2;
	fParameters[27] = fDCAToPVMin2;
	fParameters[28] = fCOADaughterMin2;
	fParameters[29] = fDCAZDaughterMax2;
	fParameters[30] = fWindowV02;
	fParameters[31] = fTrackPMin;

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
	// Stores relevant V0s and tracks in array
	//------------------------------------------------------------------------------------------

  // V0s
  Int_t v0Array[nV0s];
  Int_t nV0Survived=0;
  for (Int_t iV0=0; iV0<nV0s; iV0++) {
    AliESDv0 *v0rel = fESDEvent->GetV0(iV0);
    AliESDtrack *trkPrel = fESDEvent->GetTrack(v0rel->GetPindex());
    AliESDtrack *trkNrel = fESDEvent->GetTrack(v0rel->GetNindex());
    if ( trkPrel->GetSign()>0 && trkNrel->GetSign()<0 ) {
    } else if ( trkPrel->GetSign()<0 && trkNrel->GetSign()>0 ) {
      trkPrel = fESDEvent->GetTrack(v0rel->GetNindex());
      trkNrel = fESDEvent->GetTrack(v0rel->GetPindex());
    } else continue;
    if( !fESDtrackCuts->AcceptTrack(trkPrel) ) continue;
    if( !fESDtrackCuts->AcceptTrack(trkNrel) ) continue;

		Double_t v0pos[3], v0neg[3];
		if ( trkPrel->GetSign()>0 && trkNrel->GetSign()<0 ) {
			v0rel->GetPPxPyPz(v0pos[0],v0pos[1],v0pos[2]);
			v0rel->GetNPxPyPz(v0neg[0],v0neg[1],v0neg[2]);
		} else if ( trkPrel->GetSign()<0 && trkNrel->GetSign()>0 ) {
			v0rel->GetPPxPyPz(v0neg[0],v0neg[1],v0neg[2]);
			v0rel->GetNPxPyPz(v0pos[0],v0pos[1],v0pos[2]);
		}

		Double_t invmassL = -9.;
		Double_t returnv0[1];
		returnv0[0] = 0.;
		invmassL = InvMassLambda(v0pos,v0neg,trkPrel,trkNrel,returnv0);
		if ( TMath::Abs(invmassL-mLambdaPDG)>fWindowV02 ) continue;


		if ( fLikeSignDB == 0 ) {//All
			if ( returnv0[0]==0 ) continue;
		} else if ( fLikeSignDB==1 || fLikeSignDB==4 ) {//???-Lambda
			if ( returnv0[0]!=1 ) continue;
		} else if ( fLikeSignDB==2 || fLikeSignDB==3 ) {//???-Lambdabar
			if ( returnv0[0]!=-1 ) continue;
		}

		v0Array[nV0Survived++] = iV0;
  }


	//------------------------------------------------------------------------------------------
	// Cascade loop (To find Xi) (START)
	//------------------------------------------------------------------------------------------

	if ( fRecoSelfCasc==1 ) {//Reconstruction of cascade is done by myself
		if (!ReconstructCascade()) return;
	}
	Int_t nCasc = fESDEvent->GetNumberOfCascades();
  if (nCasc==0) {
    return;
  }
//	printf("\n### reconstructed nCascade:%d\n\n",nCasc);

	for (Int_t iCasc1=0; iCasc1<nCasc; iCasc1++) {
		AliESDcascade *casc1 = fESDEvent->GetCascade(iCasc1);

		// Get track information
		AliESDtrack *trkP1 = fESDEvent->GetTrack(casc1->GetPindex());
		Double_t trkP1Charge = trkP1->GetSign();
		AliESDtrack *trkN1 = fESDEvent->GetTrack(casc1->GetNindex());
		Double_t trkN1Charge = trkN1->GetSign();
		AliESDtrack *trkB1 = fESDEvent->GetTrack(casc1->GetBindex());
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
		Double_t trkB1TPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkB1, AliPID::kPion  ));

		if ( trkB1TPCPion>fReqSigmaTPC ) continue;
		if ( !(   (trkP1TPCProton<=fReqSigmaTPC && trkN1TPCPion<=fReqSigmaTPC  )||
		          (trkP1TPCPion<=fReqSigmaTPC   && trkN1TPCProton<=fReqSigmaTPC)   ) ) continue;


	//------------------------------------------------------------------------------------------
	// V0 loop (To find Lambda) (START)
	//------------------------------------------------------------------------------------------

	for (Int_t iV02S=0; iV02S<nV0Survived; iV02S++) {
		Int_t iV02 = v0Array[iV02S];
		AliESDv0 *esdV02 = fESDEvent->GetV0(iV02);
		Bool_t IsOnFly = esdV02->GetOnFlyStatus();
      
		AliESDtrack *ptrk = fESDEvent->GetTrack(esdV02->GetPindex());
		AliESDtrack *ntrk = fESDEvent->GetTrack(esdV02->GetNindex());
		if( !fESDtrackCuts->AcceptTrack(ptrk) ) continue;
		if( !fESDtrackCuts->AcceptTrack(ntrk) ) continue;
      
		Double_t ptrkCharge = ptrk->GetSign();
		Double_t ntrkCharge = ntrk->GetSign();
		if        ( ptrkCharge>0 && ntrkCharge>0 ) {//(+,+)->continue
			continue;
		} else if ( ptrkCharge<0 && ntrkCharge<0 ) {//(-,-)->continue
			continue;
		} else if ( ptrkCharge>0 && ntrkCharge<0 ) {//(+,-)->ok
		} else if ( ptrkCharge<0 && ntrkCharge>0 ) {//(-,+)->(+,-)
			ptrk = fESDEvent->GetTrack(esdV02->GetNindex());
			ntrk = fESDEvent->GetTrack(esdV02->GetPindex());
		} 
      
		Double_t ptrkTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrk, AliPID::kProton));
		Double_t ptrkTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrk, AliPID::kPion  ));
		Double_t ntrkTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrk, AliPID::kProton));
		Double_t ntrkTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrk, AliPID::kPion  ));

		Int_t pos2ID = ptrk->GetID();
		Int_t neg2ID = ntrk->GetID();

		if (pos2ID==trkP1ID || pos2ID==trkN1ID || pos2ID==trkB1ID) continue;
		if (neg2ID==trkP1ID || neg2ID==trkN1ID || neg2ID==trkB1ID) continue;

		Double_t MomV02[3], PosV02[3];
		Double_t MomV02Pos[3], MomV02Neg[3];
		esdV02->PxPyPz(MomV02);
		esdV02->XvYvZv(PosV02);
		if ( ptrkCharge>0 && ntrkCharge<0 ) {//(+,-)->ok
			esdV02->GetPPxPyPz(MomV02Pos[0],MomV02Pos[1],MomV02Pos[2]);
			esdV02->GetNPxPyPz(MomV02Neg[0],MomV02Neg[1],MomV02Neg[2]);
		} else if ( ptrkCharge<0 && ntrkCharge>0 ) {//(-,+)->(+,-)
			esdV02->GetPPxPyPz(MomV02Neg[0],MomV02Neg[1],MomV02Neg[2]);
			esdV02->GetNPxPyPz(MomV02Pos[0],MomV02Pos[1],MomV02Pos[2]);
		}
		Double_t dca2 = esdV02->GetDcaV0Daughters();
		if ( dca2>fDCADaughterMax2 ) continue;

		Double_t cpa2ESD  = esdV02->GetV0CosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2]);
		Double_t ip2ESD   = esdV02->GetD(PosPV[0],PosPV[1],PosPV[2]);
		if ( cpa2ESD < fCPAMin2 ) continue;
		if ( ip2ESD < fDCAToPVMin2 ) continue;

    Double_t Mom2PosPt  = GetPtFromPxPyPz(MomV02Pos);
    Double_t Mom2NegPt  = GetPtFromPxPyPz(MomV02Neg);
    Double_t Mom2PosP   = GetPaFromPxPyPz(MomV02Pos);
    Double_t Mom2NegP   = GetPaFromPxPyPz(MomV02Neg);
    Double_t coa2PosNeg = (MomV02Pos[0]*MomV02Neg[0]+MomV02Pos[1]*MomV02Neg[1]+MomV02Pos[2]*MomV02Neg[2])/Mom2PosP/Mom2NegP;
		if ( coa2PosNeg < fCOADaughterMin2 ) continue;

    Double_t pDca = ptrk->GetD(PosPV[0],PosPV[1],fBzkG);
    Double_t nDca = ntrk->GetD(PosPV[0],PosPV[1],fBzkG);
    if ( TMath::Abs(pDca)<fPosDCAToPVMin2 ) continue;
    if ( TMath::Abs(nDca)<fNegDCAToPVMin2 ) continue;

    Double_t InvMassLambdaV02 = -9.;
    Double_t L2ReturnV0[1];
    L2ReturnV0[0] = 0.;
    InvMassLambdaV02     = InvMassLambda(MomV02Pos,MomV02Neg,ptrk,ntrk,L2ReturnV0);
		Int_t proLambda2;
		proLambda2  = static_cast<Int_t>(L2ReturnV0[0]);

    if ( TMath::Abs(InvMassLambdaV02-mLambdaPDG)>fWindowV02 ) continue;


	//------------------------------------------------------------------------------------------
	// Get Cascade/V0 momenta and positions
	//------------------------------------------------------------------------------------------

			Double_t PosCasc1[3];
			casc1->XvYvZv(PosCasc1);

			Double_t PVtoCasc1[3], PVtoV02[3] ;
			for (Int_t i=0; i<3; i++) {
				PVtoCasc1[i] = PosCasc1[i] - PosPV[i];
				PVtoV02[i]   = PosV02[i]   - PosPV[i];
			}
			Double_t rCasc1 = GetPaFromPxPyPz(PVtoCasc1); 
			Double_t rV02   = GetPaFromPxPyPz(PVtoV02); 
			Double_t ctauX1 = rCasc1*mXiPDG/casc1->P();
			Double_t ctauL2 = rV02*mLambdaPDG/GetPaFromPxPyPz(MomV02);

			Double_t MomCasc1[3];
			casc1->PxPyPz(MomCasc1);

			Double_t PosV01[3];
			casc1->GetXYZ(PosV01[0],PosV01[1],PosV01[2]);

			Double_t MomP1[3];
			Double_t MomN1[3];
			Double_t MomB1[3];
			Double_t MomL1[3];
			Double_t MomX1[3];
			casc1->GetPPxPyPz(MomP1[0],MomP1[1],MomP1[2]);
			casc1->GetNPxPyPz(MomN1[0],MomN1[1],MomN1[2]);
			casc1->GetBPxPyPz(MomB1[0],MomB1[1],MomB1[2]);
			for (Int_t i=0; i<3; i++) {
				MomL1[i] = MomP1[i] + MomN1[i];
				MomX1[i] = MomL1[i] + MomB1[i];
			}

			Double_t MomP1P = GetPaFromPxPyPz(MomP1);
			Double_t MomN1P = GetPaFromPxPyPz(MomN1);
			Double_t MomB1P = GetPaFromPxPyPz(MomB1);

			Double_t Casc1toV01[3];
			for (Int_t i=0; i<3; i++) {
				Casc1toV01[i] = PosV01[i] - PosCasc1[i];
			}
			Double_t rV01    = GetPaFromPxPyPz(Casc1toV01); 
			Double_t ctauV01 = rV01*mLambdaPDG/GetPaFromPxPyPz(MomL1);

	//------------------------------------------------------------------------------------------
	// Get Dibaryon momenta and positions
	//------------------------------------------------------------------------------------------

			AliESDv0 v0Casc(*esdV02);
			AliESDv0 *v0Propagated = &v0Casc;
			AliExternalTrackParam *trkLikeCasc1 = new AliExternalTrackParam();

			Double_t cvDef[21];
			for (Int_t i=0; i<21; i++) cvDef[i] = 0.;
			cvDef[ 0] = 1.;
			cvDef[ 2] = 1.;
			cvDef[ 5] = 1.;
			cvDef[ 9] = 1.;
			cvDef[14] = 1.;
			cvDef[20] = 1.;

			trkLikeCasc1->Set(PosCasc1,MomCasc1,cvDef,static_cast<Short_t>(trkB1Charge));
			//trkLikeCasc1->Print();

			Double_t DCADB = 999.;
			DCADB = PropagateToDCA(v0Propagated,trkLikeCasc1,fBzkG);
			if ( DCADB>fDCADibaryon ) continue;

			Int_t fakeBachId = -1;
			AliESDcascade *cascade = new AliESDcascade(*v0Propagated,*trkLikeCasc1,fakeBachId);

			Double_t PosDB[3];
			cascade->GetXYZcascade(PosDB[0],PosDB[1],PosDB[2]);

			Double_t MomDB[3];
			cascade->GetPxPyPz(MomDB[0],MomDB[1],MomDB[2]);

			Double_t MomV02DB[3];
			v0Propagated->GetPxPyPz(MomV02DB[0],MomV02DB[1],MomV02DB[2]);

			Double_t PosV02DB[3];
			v0Propagated->GetXYZ(PosV02DB[0],PosV02DB[1],PosV02DB[2]);

			Double_t MomCasc1DB[3];
			cascade->GetBPxPyPz(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2]);

			Double_t PosCasc1DB[3];
			trkLikeCasc1->GetXYZ(PosCasc1DB);

//			printf("\n### Pos DB(%10f,%10f,%10f) SUM(%10f,%10f,%10f)\n",PosDB[0],PosDB[1],PosDB[2],(PosCasc1DB[0]+PosV02DB[0])/2.,(PosCasc1DB[1]+PosV02DB[1])/2.,(PosCasc1DB[2]+PosV02DB[2])/2.);
//			printf("### Pos V0(%10f,%10f,%10f) Cas(%10f,%10f,%10f)\n",PosV02DB[0],PosV02DB[1],PosV02DB[2],PosCasc1DB[0],PosCasc1DB[1],PosCasc1DB[2]);
//			printf("### Mom DB(%10f,%10f,%10f) SUM(%10f,%10f,%10f)\n",MomDB[0],MomDB[1],MomDB[2],MomCasc1DB[0]+MomV02DB[0],MomCasc1DB[1]+MomV02DB[1],MomCasc1DB[2]+MomV02DB[2]);

			// get dibaryon momentum and position
			Double_t DBtoCasc1[3], DBtoV02[3];
			Double_t PVtoDB[3];
			for (Int_t i=0; i<3; i++) {
				DBtoCasc1[i] = PosCasc1[i] - PosDB[i];
				DBtoV02[i]   = PosV02[i] - PosDB[i];
				PVtoDB[i]    = PosDB[i] - PosPV[i];
			}

			Double_t rCasc1toDB = GetPaFromPxPyPz(DBtoCasc1);
			Double_t rV02toDB   = GetPaFromPxPyPz(DBtoV02);
			Double_t rDB        = GetPaFromPxPyPz(PVtoDB);
			Double_t ctauXitoDB = rCasc1toDB/mXiPDG/GetPaFromPxPyPz(MomCasc1);
			Double_t ctauLtoDB  = rV02toDB/mLambdaPDG/GetPaFromPxPyPz(MomV02);
			Double_t ctauDB     = rDB/(mOmegaPDG+mProtonPDG)/GetPaFromPxPyPz(MomDB);
			Double_t cpaXitoDB  = (DBtoCasc1[0]*MomCasc1DB[0]+DBtoCasc1[1]*MomCasc1DB[1]+DBtoCasc1[2]*MomCasc1DB[2])/rCasc1toDB/GetPaFromPxPyPz(MomCasc1DB);
			Double_t cpaLtoDB   = (DBtoV02[0]*MomV02DB[0]+DBtoV02[1]*MomV02DB[1]+DBtoV02[2]*MomV02DB[2])/rV02toDB/GetPaFromPxPyPz(MomV02DB);
			Double_t cpaDB      = (PVtoDB[0]*MomDB[0]+PVtoDB[1]*MomDB[1]+PVtoDB[2]*MomDB[2])/rDB/GetPaFromPxPyPz(MomDB);
			if ( cpaXitoDB < 0. ) continue;
			if ( cpaLtoDB < 0. ) continue;
			if ( cpaDB < fCPADibaryon ) continue;

			if(trkLikeCasc1) delete trkLikeCasc1;
			if(cascade) delete cascade;

			Double_t coaCasc1V02 = (MomCasc1DB[0]*MomV02DB[0]+MomCasc1DB[1]*MomV02DB[1]+MomCasc1DB[2]*MomV02DB[2])/GetPaFromPxPyPz(MomCasc1DB)/GetPaFromPxPyPz(MomV02DB);
			Double_t caPVCasc1   = (PosCasc1[0]*PosDB[0]+PosCasc1[1]*PosDB[1]+PosCasc1[2]*PosDB[2])/GetPaFromPxPyPz(PosCasc1)/GetPaFromPxPyPz(PosDB);
			Double_t caPVV02     = (PosV02[0]*PosDB[0]+PosV02[1]*PosDB[1]+PosV02[2]*PosDB[2])/GetPaFromPxPyPz(PosV02)/GetPaFromPxPyPz(PosDB);
	//------------------------------------------------------------------------------------------
	// PID information
	//------------------------------------------------------------------------------------------

			// proLambda or antiLambda
			Int_t proLambda1  = 0;

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

			// Xi or Anti-Xi
			Int_t typeOfCasc1 = 0;//2:Xi-(Lpi-), -2:Xi+(Lbarpi+)

			if ( proLambda1 == 1 ) {//pro-Lambda
				if ( trkB1Charge == -1 ) {//minus bachelor
					typeOfCasc1 = 2;//Xi-
				} else if ( trkB1Charge ==  1 ) {//plus bachelor
					typeOfCasc1 = 4;//Lpi+
				} else typeOfCasc1 = -99;
			} else if ( proLambda1 == -1 ) {//anti-Lambda
				if ( trkB1Charge == -1 ) {//minus bachelor
					typeOfCasc1 = -4;//Lbarpi-
				} else if ( trkB1Charge ==  1 ) {//plus bachelor
					typeOfCasc1 = -2;//Xi+
				} else typeOfCasc1 = -99;
			} else typeOfCasc1 = -99;

			if ( fLikeSignDB==1 || fLikeSignDB==3 ) {//(Xi-)-???
				if (typeOfCasc1!=2) continue;
			} else if ( fLikeSignDB==2 || fLikeSignDB==4 ) {//(Xi+)-???
				if (typeOfCasc1!=-2) continue;
			}

/*			printf("\n### 1:P(%2d):P=%8f, pi=%8f, N(%2d):P=%8f, pi=%8f, B(%2d):K=%8f, pi=%8f\n",static_cast<Int_t>(trkP1Charge),trkP1TPCProton,trkP1TPCPion,static_cast<Int_t>(trkN1Charge),trkN1TPCProton,trkN1TPCPion,static_cast<Int_t>(trkB1Charge),trkB1TPCKaon,trkB1TPCPion);
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
			TLorentzVector vP2,vN2,vL2;
			TLorentzVector vLm1,vLm2,vCs1,vCs1L;
			TLorentzVector vCsV0,vCsV0LL,vCsV0CL;
			TLorentzVector vC1DB,vC1DBL,vC1DBC,vL2DB,vL2DBL;
			TLorentzVector vDB,vDBLL,vDBCL;
			if      (proLambda1== 1) vP1.SetXYZM(MomP1[0],MomP1[1],MomP1[2],mProtonPDG);
			else if (proLambda1==-1) vP1.SetXYZM(MomP1[0],MomP1[1],MomP1[2],mPionPDG  );
			else                     vP1.SetXYZM(0.,0.,0.,0.);
			if      (proLambda1== 1) vN1.SetXYZM(MomN1[0],MomN1[1],MomN1[2],mPionPDG  );
			else if (proLambda1==-1) vN1.SetXYZM(MomN1[0],MomN1[1],MomN1[2],mProtonPDG);
			else                     vN1.SetXYZM(0.,0.,0.,0.);
			                         vB1.SetXYZM(MomB1[0],MomB1[1],MomB1[2],mPionPDG);
			                         vL1.SetXYZM(MomL1[0],MomL1[1],MomL1[2],mLambdaPDG);
			                         vC1.SetXYZM(MomCasc1[0],MomCasc1[1],MomCasc1[2],mXiPDG   );

			if      (proLambda2== 1) vP2.SetXYZM(MomV02Pos[0],MomV02Pos[1],MomV02Pos[2],mProtonPDG);
			else if (proLambda2==-1) vP2.SetXYZM(MomV02Pos[0],MomV02Pos[1],MomV02Pos[2],mPionPDG  );
			else                     vP2.SetXYZM(0.,0.,0.,0.);
			if      (proLambda2== 1) vN2.SetXYZM(MomV02Neg[0],MomV02Neg[1],MomV02Neg[2],mPionPDG  );
			else if (proLambda2==-1) vN2.SetXYZM(MomV02Neg[0],MomV02Neg[1],MomV02Neg[2],mProtonPDG);
			else                     vN2.SetXYZM(0.,0.,0.,0.);
			                         vL2.SetXYZM(MomV02[0],MomV02[1],MomV02[2],mLambdaPDG);

			                         vC1DBC.SetXYZM(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2],mXiPDG   );
			                         vL2DBL.SetXYZM(MomV02DB[0],MomV02DB[1],MomV02DB[2],mLambdaPDG   );

			vLm1    = vP1    + vN1;
			vLm2    = vP2    + vN2;
			vCs1    = vLm1   + vB1;
			vCs1L   = vL1    + vB1;
			vCsV0   = vCs1   + vLm2;
			vCsV0LL = vCs1L  + vL2;
			vCsV0CL = vC1    + vL2;
			vDBCL   = vC1DBC + vL2DBL;

			vC1DBL.SetXYZM(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2],vCs1L.M());
			vC1DB.SetXYZM(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2],vCs1.M());
			vL2DB.SetXYZM(MomV02DB[0],MomV02DB[1],MomV02DB[2],vLm2.M());

			vDBLL   = vC1DBL + vL2DBL;
			vDB     = vC1DB  + vL2DB;


	//------------------------------------------------------------------------------------------
	// Other information
	//------------------------------------------------------------------------------------------

			Double_t CPACasc1 = casc1->GetCascadeCosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2]);
			Double_t IPcasc1  = casc1->GetDcascade(PosPV[0],PosPV[1],PosPV[2]);
			Double_t DCAXiDa1 = casc1->GetDcaXiDaughters();

			Double_t CPAV01   = casc1->GetV0CosineOfPointingAngle(PosCasc1[0],PosCasc1[1],PosCasc1[2]);
			Double_t DCAV0Da1 = casc1->GetDcaV0Daughters();

			// Armenteros-Podolanski (self)
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
			TVector3 momentumPosV02(MomV02Pos[0],MomV02Pos[1],MomV02Pos[2]);
			TVector3 momentumNegV02(MomV02Neg[0],MomV02Neg[1],MomV02Neg[2]);
			TVector3 momentumTotV02(MomV02Pos[0]+MomV02Neg[0],MomV02Pos[1]+MomV02Neg[1],MomV02Pos[2]+MomV02Neg[2]);
			Double_t qPosV02 = momentumPosV02.Dot(momentumTotV02)/momentumTotV02.Mag();
			Double_t qNegV02 = momentumNegV02.Dot(momentumTotV02)/momentumTotV02.Mag();
			AlphaV02 = (qPosV02-qNegV02)/(qPosV02+qNegV02);
			PtArmV02 = momentumNegV02.Perp(momentumTotV02);

			Double_t AlphaDB = 0.;
			Double_t PtArmDB = 0.;
			TVector3 momentumDB1(MomCasc1DB[0],MomCasc1DB[1],MomCasc1DB[2]);
			TVector3 momentumDB2(MomV02DB[0],MomV02DB[1],MomV02DB[2]);
			TVector3 momentumTotDB(MomCasc1DB[0]+MomV02DB[0],MomCasc1DB[1]+MomV02DB[1],MomCasc1DB[2]+MomV02DB[2]);
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
			fCandidateVariables[ 4] = vCsV0.M();
			fCandidateVariables[ 5] = vCsV0LL.M();
			fCandidateVariables[ 6] = vCsV0CL.M();
			fCandidateVariables[ 7] = vDB.M();
			fCandidateVariables[ 8] = vDBLL.M();
			fCandidateVariables[ 9] = vDBCL.M();
			fCandidateVariables[10] = rCasc1;
			fCandidateVariables[11] = ctauX1;
			fCandidateVariables[12] = ctauL2;
			fCandidateVariables[13] = rV01;
			fCandidateVariables[14] = rV02;
			fCandidateVariables[15] = ctauV01;
			fCandidateVariables[16] = CPACasc1;
			fCandidateVariables[17] = cpa2ESD;
			fCandidateVariables[18] = IPcasc1;
			fCandidateVariables[19] = ip2ESD;
			fCandidateVariables[20] = DCAXiDa1;
			fCandidateVariables[21] = dca2;
			fCandidateVariables[22] = CPAV01;
			fCandidateVariables[23] = DCAV0Da1;
			fCandidateVariables[24] = DCADB;
			fCandidateVariables[25] = ctauDB;
			fCandidateVariables[26] = cpaDB;
			fCandidateVariables[27] = typeOfCasc1;
			fCandidateVariables[28] = proLambda2;
			fCandidateVariables[29] = static_cast<Float_t>(triggerType);
			fCandidateVariables[30] = AlphaV01;
			fCandidateVariables[31] = PtArmV01;
			fCandidateVariables[32] = AlphaV02;
			fCandidateVariables[33] = PtArmV02;
			fCandidateVariables[34] = AlphaDB;
			fCandidateVariables[35] = PtArmDB;
			fCandidateVariables[36] = ctauXitoDB;
			fCandidateVariables[37] = ctauLtoDB;
			fCandidateVariables[38] = cpaXitoDB;
			fCandidateVariables[39] = cpaLtoDB;
			fCandidateVariables[40] = rCasc1toDB;
			fCandidateVariables[41] = rV02;
			fCandidateVariables[42] = coaCasc1V02;
			fCandidateVariables[43] = caPVCasc1;
			fCandidateVariables[44] = caPVV02;
			fCandidateVariables[45] = pDca;
			fCandidateVariables[46] = nDca;

			fVariablesTree->Fill();

		}//v0 loop

	}//cascade loop

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
void AliAnalysisTaskNOmegaLX::DefineTreeVariables() {
  //
  // This is to define tree variables
  //

  const char* nameoutput1 = GetOutputSlot(1)->GetContainer()->GetName();
  fParametersTree = new TTree(nameoutput1,"Parameters tree");
  Int_t nVar1 = 32;
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
  fParameterNames[ 8]="Chi2max";// = fCsChi2max;
  fParameterNames[ 9]="DV0min";// = fCsDV0min;
  fParameterNames[10]="WinLambda1";// = fCsMassWinLambda;
  fParameterNames[11]="DBachMin";// = fCsDBachMin;
  fParameterNames[12]="DCAmax";// = fCsDCAmax; 
  fParameterNames[13]="CPAmin";// = fCsCPAmin;
  fParameterNames[14]="Rmin";// = fCsRmin;
  fParameterNames[15]="Rmax";// = fCsRmax;
	fParameterNames[16]="RecoSelfCasc";// = fRecoSelfCasc;
	fParameterNames[17]="PProtonMax";// = fProtonPMax;
	fParameterNames[18]="PPionMax";// = fPionPMax;
	fParameterNames[19]="PKaonMax";// = fKaonPMax;
	fParameterNames[20]="DCADB";// = fDCADibaryon;
	fParameterNames[21]="FVMin2";// = fFiducialVolMin2;
	fParameterNames[22]="FVMax2";// = fFiducialVolMax2;
	fParameterNames[23]="IPPosMin2";// = fPosDCAToPVMin2;
	fParameterNames[24]="IPNegMin2";// = fNegDCAToPVMin2;
	fParameterNames[25]="DCADauMax2";// = fDCADaughterMax2;
	fParameterNames[26]="CPAMin2";// = fCPAMin2;
	fParameterNames[27]="DCAPVMin2";// = fDCAToPVMin2;
	fParameterNames[28]="COADauMin2";// = fCOADaughterMin2;
	fParameterNames[29]="DCAZDauMax2";// = fDCAZDaughterMax2;
	fParameterNames[30]="WinLambda2";// = fWindowV02;
	fParameterNames[31]="PTrackMin";// = fTrackPMin;

  for (Int_t ivar=0; ivar<nVar1; ivar++) {
    fParametersTree->Branch(fParameterNames[ivar].Data(),&fParameters[ivar],Form("%s/f",fParameterNames[ivar].Data()));
  }


  const char* nameoutput2 = GetOutputSlot(2)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput2,"Candidates variables tree");
  Int_t nVar2 = 47;
  fCandidateVariables = new Float_t [nVar2];
  TString * fCandidateVariableNames = new TString[nVar2];

	fCandidateVariableNames[ 0]="IML1";// = vLm1.M();
	fCandidateVariableNames[ 1]="IML2";// = vLm2.M();
	fCandidateVariableNames[ 2]="IMC1";// = vCs1.M();
	fCandidateVariableNames[ 3]="IMC1L";// = vCs1L.M();
	fCandidateVariableNames[ 4]="IMCV";// = vCsV0.M();
	fCandidateVariableNames[ 5]="IMCVL";// = vCsV0LL.M();
	fCandidateVariableNames[ 6]="IMCVC";// = vCsV0CL.M();
	fCandidateVariableNames[ 7]="IMDB";// = vDB.M();
	fCandidateVariableNames[ 8]="IMDBL";// = vDBLL.M();
	fCandidateVariableNames[ 9]="IMDBC";// = vDBCL.M();
	fCandidateVariableNames[10]="rCs1PV";// = rCasc1;
	fCandidateVariableNames[11]="ctauX1PV";// = ctauX1;
	fCandidateVariableNames[12]="ctauL2PV";// = ctauL2;
	fCandidateVariableNames[13]="rV01";// = rV01;
	fCandidateVariableNames[14]="rV02PV";// = rV02;
	fCandidateVariableNames[15]="ctauL1";// = ctauV01;
	fCandidateVariableNames[16]="CPACs1";// = CPACasc1;
	fCandidateVariableNames[17]="CPAV02";// = cpa2ESD;
	fCandidateVariableNames[18]="IPCs1";// = IPcasc1;
	fCandidateVariableNames[19]="IPV02";// = ip2ESD;
	fCandidateVariableNames[20]="DCAXiDa1";// = DCAXiDa1;
	fCandidateVariableNames[21]="DCALDa2";// = dca2;
	fCandidateVariableNames[22]="CPAV01";// = CPAV01;
	fCandidateVariableNames[23]="DCAV0Da1";// = DCAV0Da1;
	fCandidateVariableNames[24]="DCADB";// = DCADB;
	fCandidateVariableNames[25]="ctauDB";// = ctauDB;
	fCandidateVariableNames[26]="CPADB";// = cpaDB;
	fCandidateVariableNames[27]="t1";// = typeOfCasc1;
	fCandidateVariableNames[28]="proL2";// = proLambda2;
	fCandidateVariableNames[29]="TT";// = static_cast<Float_t>(triggerType);
	fCandidateVariableNames[30]="AlphaV01";// = AlphaV01;
	fCandidateVariableNames[31]="PtArmV01";// = PtArmV01;
	fCandidateVariableNames[32]="AlphaV02";// = AlphaV02;
	fCandidateVariableNames[33]="PtArmV02";// = PtArmV02;
	fCandidateVariableNames[34]="AlphaDB";// = AlphaDB;
	fCandidateVariableNames[35]="PtArmDB";// = PtArmDB;
	fCandidateVariableNames[36]="ctauXi";// = ctauXitoDB;
	fCandidateVariableNames[37]="ctauL2";// = ctauLtoDB;
	fCandidateVariableNames[38]="CPAXi";// = cpaXitoDB;
	fCandidateVariableNames[39]="CPAL2";// = cpaLtoDB;
	fCandidateVariableNames[40]="rXi";// = rCasc1toDB;
	fCandidateVariableNames[41]="rL2";// = rV02;
	fCandidateVariableNames[42]="COACsV0";// = coaCasc1V02;
	fCandidateVariableNames[43]="CACs1";// = caPVCasc1;
	fCandidateVariableNames[44]="CAV02";// = caPVV02;
	fCandidateVariableNames[45]="DcaPos2";// = pDca;
	fCandidateVariableNames[46]="DcaNeg2";// = nDca;

  for (Int_t ivar=0; ivar<nVar2; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;

}
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskNOmegaLX::PreTrackCut(AliESDtrack *track) {

	//------------------------------------------------------------------------------------------
	// version 0.11 (2015/10/23)
	// First track cut for speed up (reduce the number of roops)
	// Input:  AliESDtrack: track
	// Return: kTRUE: Survived
	//        kFALSE: Excluded
	//------------------------------------------------------------------------------------------

	Double_t trackCharge = track->GetSign();

	if ( track->GetTPCNcls()<fReqClustersTPC ) return kFALSE;

	if ( !track->IsOn(AliESDtrack::kTPCrefit) ) return kFALSE;

	if ( track->GetKinkIndex(0)>0 ) return kFALSE;

	Double_t trackTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
	Double_t trackTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion  ));
	Double_t trackTPCKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon  ));
//	if ( trackTPCProton>fReqSigmaTPC && trackTPCPion>fReqSigmaTPC && trackTPCKaon>fReqSigmaTPC ) return kFALSE;
//	if ( trackTPCProton>fReqSigmaTPC && trackTPCPion>fReqSigmaTPC ) return kFALSE;
	if ( trackTPCProton>fReqSigmaTPC && trackTPCKaon>fReqSigmaTPC ) return kFALSE;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//	if (TMath::Abs(pos1DCAPV)>0.0182+0.035/MomPos01Pt) continue;
//	if (TMath::Abs(neg1DCAPV)>0.0182+0.035/MomNeg01Pt) continue;

	Double_t Mom[3];
	track->GetPxPyPz(Mom);
	TVector3 vMom;
	vMom.SetXYZ(Mom[0],Mom[1],Mom[2]);
	if (TMath::Abs(vMom.PseudoRapidity())>fReqPseudoRap) return kFALSE;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	return kTRUE;

}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskNOmegaLX::InvMassLambda(Double_t MomPos[3], Double_t MomNeg[3], AliESDtrack *pos, AliESDtrack *neg, Double_t v0Return[1]) {

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
		if (posTOFProton>fReqSigmaTOF && posTOFPion>fReqSigmaTOF) return 0.;
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
Double_t AliAnalysisTaskNOmegaLX::InvMassLambdaStar(Double_t MomPos[3], Double_t MomNeg[3], AliESDtrack *pos, AliESDtrack *neg, Double_t v0Return[1]) {

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
Double_t AliAnalysisTaskNOmegaLX::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 2x2 determinant
  //--------------------------------------------------------------------
  return a00*a11 - a01*a10;
}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskNOmegaLX::Det(Double_t a00,Double_t a01,Double_t a02,
         Double_t a10,Double_t a11,Double_t a12,
         Double_t a20,Double_t a21,Double_t a22) const {
  //--------------------------------------------------------------------
  // This function calculates locally a 3x3 determinant
  //--------------------------------------------------------------------
  return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskNOmegaLX::PropagateToDCA(AliESDv0 *v, AliExternalTrackParam *t, Double_t b) {
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
  if (dca > fCsDCAmax) return 1.e+33;

//points of the DCA
  Double_t t1 = Det(x2-x1,y2-y1,z2-z1,px2,py2,pz2,ax,ay,az)/
                Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);

  x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;

  if (x1*x1+y1*y1 > (fCsRmax+5.)*(fCsRmax+5.)) return 1.e+33;

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
Bool_t AliAnalysisTaskNOmegaLX::ReconstructCascade() {

	//------------------------------------------------------------------------------------------
	// version 0.00 (2016/02/24)
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
//  Double_t mXi1530PDG   = TDatabasePDG::Instance()->GetParticle(3314)->Mass();//1.535000
  Double_t mOmegaPDG    = TDatabasePDG::Instance()->GetParticle(3334)->Mass();//1.672450

  const Int_t nTracks = fESDEvent->GetNumberOfTracks();
  if (nTracks==0) {
    return kFALSE;
  }

  const Int_t nV0s = fESDEvent->GetNumberOfV0s();
  if (nV0s==0) {
    return kFALSE;
  }

  Double_t fCsRmin2 = fCsRmin*fCsRmin;
  Double_t fCsRmax2 = fCsRmax*fCsRmax;

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
		if (v0rel->GetD(PosPV[0],PosPV[1],PosPV[2])<fCsDV0min) continue;
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
		if (TMath::Abs(trkrel->GetD(PosPV[0],PosPV[1],fBzkG))<fCsDBachMin) continue;
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
			if (TMath::Abs(v0Casc.GetEffMass()-mLambdaPDG)>fCsMassWinLambda) continue;

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
				if (dcaCasc > fCsDCAmax) continue;

				AliESDcascade cascade(*v0Propagated,*trkPropagated,indexTrk);

				Double_t x,y,z;
				cascade.GetXYZcascade(x,y,z);
				Double_t r2 = x*x + y*y;
				if (r2 > fCsRmax2) continue;
				if (r2 < fCsRmin2) continue;

				Double_t pxV0,pyV0,pzV0;
				v0Propagated->GetPxPyPz(pxV0,pyV0,pzV0);
				if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

				Double_t x1,y1,z1;
				v0Propagated->GetXYZ(x1,y1,z1);
				if (r2 > (x1*x1+y1*y1)) continue;

				if (cascade.GetCascadeCosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2])<fCsCPAmin) continue; //condition on the cascade pointing angle 

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
			if (TMath::Abs(v0Casc.GetEffMass()-mLambdaPDG)>fCsMassWinLambda) continue;

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
				if (dcaCasc > fCsDCAmax) continue;

				AliESDcascade cascade(*v0Propagated,*trkPropagated,indexTrk);

				Double_t x,y,z;
				cascade.GetXYZcascade(x,y,z);
				Double_t r2 = x*x + y*y;
				if (r2 > fCsRmax2) continue;
				if (r2 < fCsRmin2) continue;

				Double_t pxV0,pyV0,pzV0;
				v0Propagated->GetPxPyPz(pxV0,pyV0,pzV0);
				if (x*pxV0+y*pyV0+z*pzV0 < 0) continue; //causality

				Double_t x1,y1,z1;
				v0Propagated->GetXYZ(x1,y1,z1);
				if (r2 > (x1*x1+y1*y1)) continue;

				if (cascade.GetCascadeCosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2])<fCsCPAmin) continue; //condition on the cascade pointing angle 

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

//	printf("### looped nCascade:%d\n",nCascSurvived);

	return kTRUE;

}
//______________________________________________________________________________________________________
//__________________________________________________________________________
