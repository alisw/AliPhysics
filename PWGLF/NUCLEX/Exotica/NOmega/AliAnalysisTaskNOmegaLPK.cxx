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
#include "AliPhysicsSelection.h"
#include "AliVEvent.h"
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

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskNOmegaLPK.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNOmegaLPK)

//______________________________________________________________________________________________________
AliAnalysisTaskNOmegaLPK::AliAnalysisTaskNOmegaLPK() :
	AliAnalysisTaskSE(),
	fESDEvent(0x0),
	fVEvent(0x0),
	fESDtrackCutsNeg(0x0),
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
	fCentrality(0),
	fCountEvent(0),
	fIsDCA(0),
	fNAll(0),
	fRecoTypeV0(0),// Reconstruction type of V0 (0:N-Omega, 1:H-dibaryon)
	fLikeSignDB(1),// Like-sign of DB (0:ALL, 1:LL, 2:(Lbar)(Lbar), 3:L(Lbar), 4:(Lbar)L)
	fReqSigmaTPC(3.0),// TPC PIDcut sigma
	fReqClustersTPC(80),// TPC number of clusters
	fReqSigmaTOF(3.0),// TOF PIDcut sigma
	fReqPseudoRap(0.9),// PseudoRapidity
	fProtonPMax(999.),// Max momentum of proton
	fPionPMax(1.5),//!!CHECK ALSO ESDtrackCuts!! Max momentum of pion
	fKaonPMax(1.5),//0.6!!CHECK ALSO ESDtrackCuts!! Max momentum of kaon
	fTrackPMin(0.0),//Min momentum of track
	fFiducialVolMin1(2.0),// Min radius of the fiducial volume
	fFiducialVolMax1(200.),// Max radius of the fiducial volume
	fPosDCAToPVMin1(0.1),//2.0 Min length of impact parameter for the positive track
	fNegDCAToPVMin1(0.1),//2.0 Min length of impact parameter for the negative track
	fDCADaughterMax1(0.5),//1.0 Max DCA between the daughter tracks
	fCPAMin1(-1.),// Min cosine of V0's pointing angle to PV
	fDCAToPVMin1(0.0),// Min DCA V0 to PV
	fCOADaughterMin1(0.0),// Min cosine between the daughter tracks
	fDCAZDaughterMax1(2.0),// Max DCAZ V0 to PV
	fWindowV01(0.0045),// Mass window cut for Lambda
	fMassGammaMin1(0.05),// Min mass of gamma conversion
	fFiducialVolMin2(2.0),// Min radius of the fiducial volume
	fFiducialVolMax2(200.),// Max radius of the fiducial volume
	fPosDCAToPVMin2(0.5),//2.0 Min length of impact parameter for the positive track
	fNegDCAToPVMin2(0.5),//2.0 Min length of impact parameter for the negative track
	fDCADaughterMax2(1.0),// Max DCA between the daughter tracks
	fCPAMin2(-1.),// Min cosine of V0's pointing angle to PV
	fDCAToPVMin2(0.0),// Min DCA V0 to PV
	fCOADaughterMin2(0.0),// Min cosine between the daughter tracks
	fDCAZDaughterMax2(2.0),// Max DCAZ V0 to PV
	fWindowV02(0.0045),// Mass window cut for Lambda
	fCPAV01toV02(0.99),// Min cosine of V02's pointing angle to V01
	fCPADibaryonNO(0.99)//0.99875 Min cosine of dibaryon's pointing angle
{
}
//______________________________________________________________________________________________________
AliAnalysisTaskNOmegaLPK::AliAnalysisTaskNOmegaLPK(const Char_t* name) :
	AliAnalysisTaskSE(name),
	fESDEvent(0x0),
	fVEvent(0x0),
	fESDtrackCutsNeg(0x0),
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
	fCentrality(0),
	fCountEvent(0),
	fIsDCA(0),
	fNAll(0),
	fRecoTypeV0(0),// Reconstruction type of V0 (0:N-Omega, 1:H-dibaryon)
	fLikeSignDB(1),// Like-sign of DB (0:ALL, 1:LL, 2:(Lbar)(Lbar), 3:L(Lbar), 4:(Lbar)L)
	fReqSigmaTPC(3.0),// TPC PIDcut sigma
	fReqClustersTPC(80),// TPC number of clusters
	fReqSigmaTOF(3.0),// TOF PIDcut sigma
	fReqPseudoRap(0.9),// PseudoRapidity
	fProtonPMax(999.),// Max momentum of proton
	fPionPMax(1.5),//!!CHECK ALSO ESDtrackCuts!! Max momentum of pion
	fKaonPMax(1.5),//!!CHECK ALSO ESDtrackCuts!! Max momentum of kaon
	fTrackPMin(0.0),//Min momentum of track
	fFiducialVolMin1(2.0),// Min radius of the fiducial volume
	fFiducialVolMax1(200.),// Max radius of the fiducial volume
	fPosDCAToPVMin1(0.1),//2.0 Min length of impact parameter for the positive track
	fNegDCAToPVMin1(0.1),//2.0 Min length of impact parameter for the negative track
	fDCADaughterMax1(0.5),//1.0 Max DCA between the daughter tracks
	fCPAMin1(-1.),// Min cosine of V0's pointing angle to PV
	fDCAToPVMin1(0.0),// Min DCA V0 to PV
	fCOADaughterMin1(0.0),// Min cosine between the daughter tracks
	fDCAZDaughterMax1(2.0),// Max DCAZ V0 to PV
	fWindowV01(0.0045),// Mass window cut for Lambda
	fMassGammaMin1(0.05),// Min mass of gamma conversion
	fFiducialVolMin2(2.0),// Min radius of the fiducial volume
	fFiducialVolMax2(200.),// Max radius of the fiducial volume
	fPosDCAToPVMin2(0.5),//2.0 Min length of impact parameter for the positive track
	fNegDCAToPVMin2(0.5),//2.0 Min length of impact parameter for the negative track
	fDCADaughterMax2(1.0),// Max DCA between the daughter tracks
	fCPAMin2(-1.),// Min cosine of V0's pointing angle to PV
	fDCAToPVMin2(0.0),// Min DCA V0 to PV
	fCOADaughterMin2(0.0),// Min cosine between the daughter tracks
	fDCAZDaughterMax2(2.0),// Max DCAZ V0 to PV
	fWindowV02(0.0045),// Mass window cut for Lambda
	fCPAV01toV02(0.99),// Min cosine of V02's pointing angle to V01
	fCPADibaryonNO(0.99)//0.99875 Min cosine of dibaryon's pointing angle
{
	//
	// Constructor. Initialization of Inputs and Outputs
	//
	Info("AliAnalysisTaskNOmegaLPK","Calling Constructor");

	DefineOutput(1,TTree::Class());  //My private output
	DefineOutput(2,TTree::Class());  //My private output

	//ESD Track cuts for Negative track
	fESDtrackCutsNeg = new AliESDtrackCuts("AliESDtrackCutsNeg","AliESDtrackCutsNeg");
	fESDtrackCutsNeg->SetMinNClustersTPC(80);
	fESDtrackCutsNeg->SetAcceptKinkDaughters(kFALSE);
	fESDtrackCutsNeg->SetMaxChi2PerClusterTPC(5);
	fESDtrackCutsNeg->SetRequireTPCRefit(kTRUE);
	fESDtrackCutsNeg->SetEtaRange(-0.9,0.9);
	fESDtrackCutsNeg->SetPtRange(fTrackPMin,1.5);

	//ESD Track cuts
	fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
	fESDtrackCuts->SetMinNClustersTPC(80);
	fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
	fESDtrackCuts->SetMaxChi2PerClusterTPC(5);
	fESDtrackCuts->SetRequireTPCRefit(kTRUE);
	fESDtrackCuts->SetEtaRange(-0.9,0.9);

}
//______________________________________________________________________________________________________
AliAnalysisTaskNOmegaLPK::~AliAnalysisTaskNOmegaLPK() {
	//
	// destructor
	//
	Info("~AliAnalysisTaskNOmegaLPK","Calling Destructor");
  
	if(fESDtrackCutsNeg) delete fESDtrackCutsNeg;
//	if(fESDCutsV0) delete fESDCutsV0;
	if (fESDtrackCuts) delete fESDtrackCuts;
  if (fPIDResponse) delete  fPIDResponse;

	if (fParametersTree) {
		delete fParameters;
		fParametersTree = 0;
	}

	if (fVariablesTree) {
		delete fVariablesTree;
		fVariablesTree = 0;
	}

}
//______________________________________________________________________________________________________
void AliAnalysisTaskNOmegaLPK::Init() {
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
void AliAnalysisTaskNOmegaLPK::UserExec(Option_t *)
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
//	AliAODMCHeader *mcHeader=0;
/*
	if (fUseMCInfo) {
		// MC array need for maching
		mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
		if (!mcArray) {
			AliError("Could not find Monte-Carlo in AOD");
			return;
		}

		// load MC header
		mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
		if (!mcHeader) {
			AliError("AliAnalysisTaskNOmegaLPK::UserExec: MC header branch not found!\n");
			return;
		}

		Double_t zMCVertex = mcHeader->GetVtxZ();
	}
*/

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
void AliAnalysisTaskNOmegaLPK::Terminate(Option_t*)
{    
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
  
	//AliInfo("Terminate","");
	AliAnalysisTaskSE::Terminate();
  
	return;
}

//______________________________________________________________________________________________________
void AliAnalysisTaskNOmegaLPK::UserCreateOutputObjects() 
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
void AliAnalysisTaskNOmegaLPK::MakeAnalysis(TClonesArray *mcArray,AliESDEvent *fESDEvent)
{

  //------------------------------------------------------------------------------------------
  // version PO1-0-10 (2016/03/03)
  // Reconstruct dibaryon from two V0s by ESD class and myself
  // and calculate the invariant mass of dibaryon
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
//  Double_t mXiPDG       = TDatabasePDG::Instance()->GetParticle(3312)->Mass();//1.321710
  Double_t mOmegaPDG    = TDatabasePDG::Instance()->GetParticle(3334)->Mass();//1.672450

  const Int_t nTracks = fESDEvent->GetNumberOfTracks();
  if (nTracks==0) {
    return;
  }

	const Int_t nV0s = fESDEvent->GetNumberOfV0s();
  if (nV0s==0) {
    return;
  }

  // Initialization of parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Int_t posID[nV0s];
	Int_t negID[nV0s];
	Double_t infoV0Tracks[nV0s][15];

	Int_t nLambda     = 0;
	Int_t nAntiLambda = 0;

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Output cut parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	fParameters[ 0] = fRecoTypeV0;
	fParameters[ 1] = fReqSigmaTPC;
	fParameters[ 2] = fReqClustersTPC;
	fParameters[ 3] = fReqSigmaTOF;
	fParameters[ 4] = fReqPseudoRap;
	fParameters[ 5] = fFiducialVolMin1;
	fParameters[ 6] = fFiducialVolMax1;
	fParameters[ 7] = fPosDCAToPVMin1;
	fParameters[ 8] = fNegDCAToPVMin1;
	fParameters[ 9] = fDCADaughterMax1;
	fParameters[10] = fCPAMin1;
	fParameters[11] = fDCAToPVMin1;
	fParameters[12] = fCOADaughterMin1;
	fParameters[13] = fDCAZDaughterMax1;
	fParameters[14] = fWindowV01;
	fParameters[15] = fFiducialVolMin2;
	fParameters[16] = fFiducialVolMax2;
	fParameters[17] = fPosDCAToPVMin2;
	fParameters[18] = fNegDCAToPVMin2;
	fParameters[19] = fDCADaughterMax2;
	fParameters[20] = fCPAMin2;
	fParameters[21] = fDCAToPVMin2;
	fParameters[22] = fCOADaughterMin2;
	fParameters[23] = fDCAZDaughterMax2;
	fParameters[24] = fWindowV02;
	fParameters[25] = fCPAV01toV02;
	fParameters[26] = fCPADibaryonNO;
	fParameters[27] = fProtonPMax;
	fParameters[28] = fPionPMax;
	fParameters[29] = fKaonPMax;
	fParameters[30] = fLikeSignDB;
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
	AliCentrality *cent = (AliCentrality*)fESDEvent->GetCentrality();
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
//	printf("triggerType:%d (Cn:%d, SC:%d, MB:%d)\n",triggerType,static_cast<Int_t>(isSelectedCn),static_cast<Int_t>(isSelectedSC),static_cast<Int_t>(isSelectedMB));

//	if (fCountEvent<setStartNumber) return;
//	if (fCountEvent>2) return;


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

  // tracks
  Int_t trkPArray[nTracks];
  Int_t trkNArray[nTracks];
  Int_t nTrkPSurvived=0; 
  Int_t nTrkNSurvived=0; 
  for (Int_t iTrk=0; iTrk<nTracks; iTrk++) {
    AliESDtrack *trkrel = fESDEvent->GetTrack(iTrk);
		Double_t trkrelCharge = trkrel->GetSign();
    if( !fESDtrackCuts->AcceptTrack(trkrel) ) continue;
    if( trkrelCharge<0 ) {
			if (!fESDtrackCutsNeg->AcceptTrack(trkrel))  continue;
		}

		Double_t trkrelDCAPV = trkrel->GetD(PosPV[0],PosPV[1],fBzkG);
		if ( TMath::Abs(trkrelDCAPV)<fPosDCAToPVMin1 && TMath::Abs(trkrelDCAPV)<fNegDCAToPVMin1 ) continue;

    Double_t trkrelTPCKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkrel, AliPID::kKaon));
    Double_t trkrelTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkrel, AliPID::kPion));
    Double_t trkrelTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trkrel, AliPID::kProton));

		if ( fRecoTypeV0 == 0 ) {//N-Omega
			if ( fLikeSignDB == 0 ) {//All
				if ( trkrelTPCProton>fReqSigmaTPC && trkrelTPCKaon>fReqSigmaTPC ) continue;
				if ( trkrelTPCProton>fReqSigmaTPC && trkrel->P()>fKaonPMax      ) continue;
			} else if ( fLikeSignDB==1 || fLikeSignDB==3 ) {//LambdaStar-???
				if ( trkrelCharge>0 && (trkrelTPCProton>fReqSigmaTPC||trkrel->P()>fProtonPMax) ) continue;
				if ( trkrelCharge<0 && (trkrelTPCKaon>fReqSigmaTPC  ||trkrel->P()>fKaonPMax  ) ) continue;
			} else if ( fLikeSignDB==2 || fLikeSignDB==4 ) {//LambdaStarbar-???
				if ( trkrelCharge>0 && (trkrelTPCKaon>fReqSigmaTPC  ||trkrel->P()>fKaonPMax  ) ) continue;
				if ( trkrelCharge<0 && (trkrelTPCProton>fReqSigmaTPC||trkrel->P()>fProtonPMax) ) continue;
			} 
		} else if ( fRecoTypeV0 == 1) {//H-dibaryon
			if ( fLikeSignDB == 0 ) {//All
				if ( trkrelTPCProton>fReqSigmaTPC && trkrelTPCPion>fReqSigmaTPC ) continue;
				if ( trkrelTPCProton>fReqSigmaTPC && trkrel->P()>fPionPMax      ) continue;
			} else if ( fLikeSignDB==1 || fLikeSignDB==3 ) {//Lambda-???
				if ( trkrelCharge>0 && (trkrelTPCProton>fReqSigmaTPC||trkrel->P()>fProtonPMax) ) continue;
				if ( trkrelCharge<0 && (trkrelTPCPion>fReqSigmaTPC  ||trkrel->P()>fPionPMax  ) ) continue;
			} else if ( fLikeSignDB==2 || fLikeSignDB==4 ) {//Lambdabar-???
				if ( trkrelCharge>0 && (trkrelTPCPion>fReqSigmaTPC  ||trkrel->P()>fPionPMax  ) ) continue;
				if ( trkrelCharge<0 && (trkrelTPCProton>fReqSigmaTPC||trkrel->P()>fProtonPMax) ) continue;
			} 
		}
    if ( trkrelCharge>0 ) trkPArray[nTrkPSurvived++] = iTrk;
    if ( trkrelCharge<0 ) trkNArray[nTrkNSurvived++] = iTrk;
  }

//  printf("\n\n### nTracks:%d, pSurvive:%d, nSurvive:%d\n",nTracks,nTrkPSurvived,nTrkNSurvived);
//  printf("### nV0s:%d, nV0Survived:%d\n",nV0s,nV0Survived);


	//------------------------------------------------------------------------------------------
	// Track loop 1 (To find V01) (START)
	//------------------------------------------------------------------------------------------

	for (Int_t iTrk1=0; iTrk1<nTrkPSurvived; iTrk1++) {
		Int_t iPos = trkPArray[iTrk1];
		AliESDtrack *pos1Trk=fESDEvent->GetTrack(iPos);
		Double_t trk1Charge = pos1Trk->GetSign();
		if (trk1Charge<0) printf("### BUG!!! at P\n");

		Double_t pos1TPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos1Trk, AliPID::kProton));
		Double_t pos1TPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos1Trk, AliPID::kPion  ));
		Double_t pos1TPCKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos1Trk, AliPID::kKaon  ));
	//------------------------------------------------------------------------------------------
	// Track loop 2 (To find V01) (START)
	//------------------------------------------------------------------------------------------

		for (Int_t iTrk2=0; iTrk2<nTrkNSurvived; iTrk2++) {
			Int_t iNeg = trkNArray[iTrk2];
			AliESDtrack *neg1Trk=fESDEvent->GetTrack(iNeg);
			Double_t trk2Charge = neg1Trk->GetSign();
			if (trk2Charge>0) printf("### BUG!!! at N\n");

			if ( trk1Charge*trk2Charge > 0 ) printf("### BUG!!! at PN\n");

			Double_t neg1TPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg1Trk, AliPID::kProton));
			Double_t neg1TPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg1Trk, AliPID::kPion  ));
			Double_t neg1TPCKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg1Trk, AliPID::kKaon  ));

	//------------------------------------------------------------------------------------------
	// Preparation to reconstruct V01
	//------------------------------------------------------------------------------------------

			Double_t pos1DCAPV = pos1Trk->GetD(PosPV[0],PosPV[1],fBzkG);
			Double_t neg1DCAPV = neg1Trk->GetD(PosPV[0],PosPV[1],fBzkG);
			if ( TMath::Abs(pos1DCAPV)<fPosDCAToPVMin1 ) continue;
			if ( TMath::Abs(neg1DCAPV)<fNegDCAToPVMin1 ) continue;

			if ( fRecoTypeV0 == 0 ) {//N-Omega
				if ( !(   (pos1TPCProton<=fReqSigmaTPC && neg1TPCKaon<=fReqSigmaTPC  )||
				          (pos1TPCKaon<=fReqSigmaTPC   && neg1TPCProton<=fReqSigmaTPC)   ) ) continue;
			} else if ( fRecoTypeV0 == 1) {//H-dibaryon
				if ( !(   (pos1TPCProton<=fReqSigmaTPC && neg1TPCPion<=fReqSigmaTPC  )||
				          (pos1TPCPion<=fReqSigmaTPC   && neg1TPCProton<=fReqSigmaTPC)   ) ) continue;
			}

			Int_t pos1ID = pos1Trk->GetID();
			Int_t neg1ID = neg1Trk->GetID();


	//------------------------------------------------------------------------------------------
	// V0 reconstruction at Primary vertex (V01)
	//------------------------------------------------------------------------------------------

			Double_t MomPosOri1[3], MomNegOri1[3];
			pos1Trk->GetPxPyPz(MomPosOri1);
			neg1Trk->GetPxPyPz(MomNegOri1);

			Double_t MomPosOri1P  = GetPaFromPxPyPz(MomPosOri1);
			Double_t MomPosOri1Pt = GetPtFromPxPyPz(MomPosOri1);
			Double_t MomNegOri1P  = GetPaFromPxPyPz(MomNegOri1);
			Double_t MomNegOri1Pt = GetPtFromPxPyPz(MomNegOri1);

			Double_t InvMassLambdaOri1 = -9.;
			Double_t L1ReturnOri[1];
			L1ReturnOri[0] = 0.;
			if ( fRecoTypeV0 == 0 ) {//N-Omega
				InvMassLambdaOri1 = InvMassLambdaStar(MomPosOri1,MomNegOri1,pos1Trk,neg1Trk,L1ReturnOri);
			} else if ( fRecoTypeV0 == 1) {//H-dibaryon
				InvMassLambdaOri1 = InvMassLambda(    MomPosOri1,MomNegOri1,pos1Trk,neg1Trk,L1ReturnOri);
			}


	//------------------------------------------------------------------------------------------
	// Get V01 Momentumas
	//------------------------------------------------------------------------------------------

			Double_t ReturnV01[6];//0:DCAxy, 1:CPA, 2:COA, 3:DCAz, 4:V0xy, 5:V0r
			Double_t ReturnV01ESD[6];//0:DCAxy, 1:CPA, 2:COA, 3:DCAz, 4:V0xy, 5:V0r
			Double_t MomV01[3], PosV01[3], MomV01Pos[3], MomV01Neg[3];
			Double_t MomV01ESD[3], PosV01ESD[3], MomV01PosESD[3], MomV01NegESD[3];
			Bool_t IsV01    = kFALSE;
			Bool_t IsESDv01 = kFALSE;
			IsV01    = GetSelfV0Momentum(1,pos1Trk,neg1Trk,ReturnV01,MomV01,PosV01,MomV01Pos,MomV01Neg);
//			IsESDv01 = GetEsdV0Momentum(1,pos1ID,neg1ID,ReturnV01ESD,MomV01ESD,PosV01ESD,MomV01PosESD,MomV01NegESD);

			Double_t InvMassLambdaV01 = -9.;
			Double_t L1ReturnV0[1];
			L1ReturnV0[0] = 0.;
			Double_t dcaV01       = -999.;
			Double_t cpaV01       = -999.;
			Double_t coaV01       = -999.;
			Double_t dcaZV01      = -999.;
			Double_t xyV01        = -999.;
			Double_t rV01         = -999.;
			Double_t ctauV01      = -999.;

			if ( fRecoTypeV0 == 0 ) {//N-Omega
				InvMassLambdaV01 = InvMassLambdaStar(MomV01Pos,MomV01Neg,pos1Trk,neg1Trk,L1ReturnV0);
			} else if ( fRecoTypeV0 == 1) {//H-dibaryon
				InvMassLambdaV01 = InvMassLambda(    MomV01Pos,MomV01Neg,pos1Trk,neg1Trk,L1ReturnV0);
			}
			dcaV01  = ReturnV01[0];
			cpaV01  = ReturnV01[1];
			coaV01  = ReturnV01[2];
			dcaZV01 = ReturnV01[3];
			xyV01   = ReturnV01[4];
			rV01    = ReturnV01[5];
			ctauV01 = ReturnV01[5]*(mOmegaPDG+mProtonPDG)/GetPaFromPxPyPz(MomV01);

			Double_t antiLambda1 = 0.;
			if ( L1ReturnV0[0]==1. ) {
				antiLambda1 = -1.;//1: Anti, -1: Not Anti
			} else if ( L1ReturnV0[0]==-1. ) {
				antiLambda1 = 1.;//1: Anti, -1: Not Anti
			}

			// Invariant mass of ee
			TLorentzVector vE1,vE2,vCon;
			vE1.SetXYZM(MomV01Pos[0],MomV01Pos[1],MomV01Pos[2],mElectronPDG);
			vE2.SetXYZM(MomV01Neg[0],MomV01Neg[1],MomV01Neg[2],mElectronPDG);
			vCon = vE1 + vE2;
			if (vCon.M()<fMassGammaMin1) continue;

			// Armenteros-Podolanski (self)
			Double_t AlphaV01 = 0.;
			Double_t PtArmV01 = 0.;
			TVector3 momentumPosV01(MomV01Pos[0],MomV01Pos[1],MomV01Pos[2]);
			TVector3 momentumNegV01(MomV01Neg[0],MomV01Neg[1],MomV01Neg[2]);
			TVector3 momentumTotV01(MomV01Pos[0]+MomV01Neg[0],MomV01Pos[1]+MomV01Neg[1],MomV01Pos[2]+MomV01Neg[2]);
			Double_t qPosV01 = momentumPosV01.Dot(momentumTotV01)/momentumTotV01.Mag();
			Double_t qNegV01 = momentumNegV01.Dot(momentumTotV01)/momentumTotV01.Mag();
			AlphaV01 = (qPosV01-qNegV01)/(qPosV01+qNegV01);
			PtArmV01 = momentumNegV01.Perp(momentumTotV01);

	//------------------------------------------------------------------------------------------
	// V0 Cut (V01)
	//------------------------------------------------------------------------------------------

			if (InvMassLambdaV01<1.) continue;
			if ( !IsV01 ) continue;


	//------------------------------------------------------------------------------------------
	// V0 loop (To find V02) (START)
	//------------------------------------------------------------------------------------------

	for (Int_t iV02S=0; iV02S<nV0Survived; iV02S++) {
		Int_t iV02 = v0Array[iV02S];
		AliESDv0 *esdV02 = fESDEvent->GetV0(iV02);
    Bool_t IsOnFly = esdV02->GetOnFlyStatus();
//		printf("### loop Pos:%5d/%5d, Neg:%5d/%5d, V0:%5d/%5d\n",iTrk1,nTrkPSurvived-1,iTrk2,nTrkNSurvived-1,iV02S,nV0Survived);

		AliESDtrack *ptrk = fESDEvent->GetTrack(esdV02->GetPindex());
		AliESDtrack *ntrk = fESDEvent->GetTrack(esdV02->GetNindex());
		if( !fESDtrackCuts->AcceptTrack(ptrk) ) continue;
		if( !fESDtrackCuts->AcceptTrack(ntrk) ) continue;

		Double_t ptrkCharge = ptrk->GetSign();
		Double_t ntrkCharge = ntrk->GetSign();
		if        ( ptrkCharge>0 && ntrkCharge>0 ) {//(3,4)=(+,+)->continue
			continue;
		} else if ( ptrkCharge<0 && ntrkCharge<0 ) {//(3,4)=(-,-)->continue
			continue;
		} else if ( ptrkCharge>0 && ntrkCharge<0 ) {//(3,4)=(+,-)->ok
		} else if ( ptrkCharge<0 && ntrkCharge>0 ) {//(3,4)=(-,+)->(+,-)
			ptrk = fESDEvent->GetTrack(esdV02->GetNindex());
			ntrk = fESDEvent->GetTrack(esdV02->GetPindex());
		}

		Double_t ptrkTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrk, AliPID::kProton));
		Double_t ptrkTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrk, AliPID::kPion  ));
		Double_t ntrkTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrk, AliPID::kProton));
		Double_t ntrkTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrk, AliPID::kPion  ));

		Int_t pos2ID = ptrk->GetID();
		Int_t neg2ID = ntrk->GetID();

		if (pos1ID==pos2ID) continue;
		if (neg1ID==neg2ID) continue;

		Double_t MomV02[3], PosV02[3];
		Double_t MomV02Pos[3], MomV02Neg[3];
		esdV02->PxPyPz(MomV02);
		esdV02->XvYvZv(PosV02);
		if ( ptrkCharge>0 && ntrkCharge<0 ) {//(3,4)=(+,-)->ok
			esdV02->GetPPxPyPz(MomV02Pos[0],MomV02Pos[1],MomV02Pos[2]);
			esdV02->GetNPxPyPz(MomV02Neg[0],MomV02Neg[1],MomV02Neg[2]);
		} else if ( ptrkCharge<0 && ntrkCharge>0 ) {//(3,4)=(-,+)->(+,-)
			esdV02->GetPPxPyPz(MomV02Neg[0],MomV02Neg[1],MomV02Neg[2]);
			esdV02->GetNPxPyPz(MomV02Pos[0],MomV02Pos[1],MomV02Pos[2]);
		}
		Double_t dca2 = esdV02->GetDcaV0Daughters();
		if ( dca2>fDCADaughterMax2 ) continue;

		Double_t xV01toV02  = PosV02[0] - PosV01[0];
		Double_t yV01toV02  = PosV02[1] - PosV01[1];
		Double_t zV01toV02  = PosV02[2] - PosV01[2];
		Double_t xyV01toV02 = TMath::Sqrt(xV01toV02*xV01toV02+yV01toV02*yV01toV02);
		Double_t rV01toV02  = TMath::Sqrt(xV01toV02*xV01toV02+yV01toV02*yV01toV02+zV01toV02*zV01toV02);
		Double_t MomV02Pt  = GetPtFromPxPyPz(MomV02);
		Double_t MomV02P   = GetPaFromPxPyPz(MomV02);
		Double_t ctau1to2L  = rV01toV02*mLambdaPDG/MomV02P;
		Double_t cpa2 = (xV01toV02*MomV02[0]+yV01toV02*MomV02[1]+zV01toV02*MomV02[2])/rV01toV02/MomV02P;
		if ( cpa2<fCPAV01toV02 ) continue;
		Double_t cpa2ESD  = esdV02->GetV0CosineOfPointingAngle(PosPV[0],PosPV[1],PosPV[2]);

		Double_t Mom2PosPt  = GetPtFromPxPyPz(MomV02Pos);
		Double_t Mom2NegPt  = GetPtFromPxPyPz(MomV02Neg);
		Double_t Mom2PosP   = GetPaFromPxPyPz(MomV02Pos);
		Double_t Mom2NegP   = GetPaFromPxPyPz(MomV02Neg);
		Double_t coa2PosNeg = (MomV02Pos[0]*MomV02Neg[0]+MomV02Pos[1]*MomV02Neg[1]+MomV02Pos[2]*MomV02Neg[2])/Mom2PosP/Mom2NegP;
		if ( coa2PosNeg<fCOADaughterMin2 ) continue;

		Double_t pDca = ptrk->GetD(PosV01[0],PosV01[1],fBzkG);
		Double_t nDca = ntrk->GetD(PosV01[0],PosV01[1],fBzkG);
		if ( TMath::Abs(pDca)<fPosDCAToPVMin2 ) continue;
		if ( TMath::Abs(nDca)<fNegDCAToPVMin2 ) continue;

		Double_t InvMassLambdaV02 = -9.;
		Double_t L2ReturnV0[1];
		L2ReturnV0[0] = 0.;
		InvMassLambdaV02     = InvMassLambda(MomV02Pos,MomV02Neg,ptrk,ntrk,L2ReturnV0);
		Double_t antiLambda2 = 0.;
		if ( L2ReturnV0[0]==1. ) {
			antiLambda2 = -1.;//1: Anti, -1: Not Anti
		} else if ( L2ReturnV0[0]==-1. ) {
			antiLambda2 = 1.;//1: Anti, -1: Not Anti
		}

		if ( TMath::Abs(InvMassLambdaV02-mLambdaPDG)>fWindowV02 ) continue;

		TLorentzVector vV01,vV01Ori,vV02,vV02L,vDB,vDBL,vDBOri,vDBOriL;
		vV01.SetXYZM(MomV01[0],MomV01[1],MomV01[2],InvMassLambdaV01);
		vV02.SetXYZM(MomV02[0],MomV02[1],MomV02[2],InvMassLambdaV02);
		vV02L.SetXYZM(MomV02[0],MomV02[1],MomV02[2],mLambdaPDG);
		vV01Ori.SetXYZM(MomPosOri1[0]+MomNegOri1[0],MomPosOri1[1]+MomNegOri1[1],MomPosOri1[2]+MomNegOri1[2],InvMassLambdaOri1);
		vDB     = vV01     + vV02;
		vDBL    = vV01     + vV02L;
		vDBOri  = vV01Ori  + vV02;
		vDBOriL = vV01Ori  + vV02L;
		Double_t cpaDB = (vDB.Px()*(PosV01[0]-PosPV[0])+vDB.Py()*(PosV01[1]-PosPV[1])+vDB.Pz()*(PosV01[2]-PosPV[2]))/vDB.P()/rV01;
		if ( cpaDB<fCPADibaryonNO ) continue;


	//------------------------------------------------------------------------------------------
	// Output
	//------------------------------------------------------------------------------------------

				fCandidateVariables[ 0] = antiLambda1;
				fCandidateVariables[ 1] = MomPosOri1Pt;
				fCandidateVariables[ 2] = MomNegOri1Pt;
				fCandidateVariables[ 3] = pos1DCAPV;
				fCandidateVariables[ 4] = neg1DCAPV;
				fCandidateVariables[ 5] = InvMassLambdaV01;
				fCandidateVariables[ 6] = dcaV01;
				fCandidateVariables[ 7] = cpaV01;
				fCandidateVariables[ 8] = coaV01;
				fCandidateVariables[ 9] = dcaZV01;
				fCandidateVariables[10] = xyV01;
				fCandidateVariables[11] = rV01;
				fCandidateVariables[12] = ctauV01;
				fCandidateVariables[13] = fCentrality;
				fCandidateVariables[14] = cpa2;
				fCandidateVariables[15] = InvMassLambdaV02;
				fCandidateVariables[16] = vDB.M();
				fCandidateVariables[17] = cpaDB;
				fCandidateVariables[18] = static_cast<Float_t>(IsOnFly);
				fCandidateVariables[19] = pDca;
				fCandidateVariables[20] = nDca;
				fCandidateVariables[21] = coa2PosNeg;
				fCandidateVariables[22] = dca2;
				fCandidateVariables[23] = xyV01toV02;
				fCandidateVariables[24] = rV01toV02;
				fCandidateVariables[25] = antiLambda2;
				fCandidateVariables[26] = InvMassLambdaOri1;
				fCandidateVariables[27] = vDBOri.M();
				fCandidateVariables[28] = static_cast<Float_t>(IsESDv01);
				fCandidateVariables[29] = cpa2ESD;
				fCandidateVariables[30] = AlphaV01;
				fCandidateVariables[31] = PtArmV01;
				fCandidateVariables[32] = vCon.M();
				fCandidateVariables[33] = vDBL.M();
				fCandidateVariables[34] = vDBOriL.M();
				fCandidateVariables[35] = ctau1to2L;
				fCandidateVariables[36] = static_cast<Float_t>(triggerType);

				fVariablesTree->Fill();

			}//v0 loop

		}//track loop 2

	}//track loop 1


/*
	Double_t ratioLambda = static_cast<Double_t>(nAntiLambda)/static_cast<Double_t>(nLambda);
	Double_t ratioLambdaAll = static_cast<Double_t>(fCountAntiLambda)/static_cast<Double_t>(fCountLambda);
	printf("#######################################################################################################################\n");
	printf("#####               Lambda: Pro:%7d->%7d,  Anti:%7d->%7d  (ratio:%6.4f->%6.4f) ######################\n",nLambda,fCountLambda,nAntiLambda,fCountAntiLambda,ratioLambda,ratioLambdaAll);
	printf("#######################################################################################################################\n");
	printf("##### Analyzed events: %7d, Run Number: %10d, Centrality: %d(%d), nV0s: %5d ################################\n",fCountEvent,runNumber,fCentrality,fESDEvent->GetEventSpecie(),nV0s);
	printf("#######################################################################################################################\n\n\n");

	fCountEvent++;
*/
}
//______________________________________________________________________________________________________
void AliAnalysisTaskNOmegaLPK::DefineTreeVariables() {
  //
  // This is to define tree variables
  //

  const char* nameoutput1 = GetOutputSlot(1)->GetContainer()->GetName();
  fParametersTree = new TTree(nameoutput1,"Parameters tree");
  Int_t nVar1 = 32;
  fParameters = new Float_t [nVar1];
  TString * fParameterNames = new TString[nVar1];

	fParameterNames[ 0]="RecoTypeV0";// = fRecoTypeV0;
	fParameterNames[ 1]="TPCSigma";// = fReqSigmaTPC;
	fParameterNames[ 2]="TPCClusters";// = fReqClustersTPC;
	fParameterNames[ 3]="TOFSigma";// = fReqSigmaTOF;
	fParameterNames[ 4]="PseudoRap";// = fReqPseudoRap;
	fParameterNames[ 5]="FVMin1";// = fFiducialVolMin1;
	fParameterNames[ 6]="FVMax1";// = fFiducialVolMax1;
	fParameterNames[ 7]="IPPOsMin1";// = fPosDCAToPVMin1;
	fParameterNames[ 8]="IPNegMin1";// = fNegDCAToPVMin1;
	fParameterNames[ 9]="DCADauMax1";// = fDCADaughterMax1;
	fParameterNames[10]="CPAMin1";// = fCPAMin1;
	fParameterNames[11]="DCAPVMin1";// = fDCAToPVMin1;
	fParameterNames[12]="COADauMin1";// = fCOADaughterMin1;
	fParameterNames[13]="DCAZDauMax1";// = fDCAZDaughterMax1;
	fParameterNames[14]="V0Win1";// = fWindowV01;
	fParameterNames[15]="FVMin2";// = fFiducialVolMin2;
	fParameterNames[16]="FVMax2";// = fFiducialVolMax2;
	fParameterNames[17]="IPPosMin2";// = fPosDCAToPVMin2;
	fParameterNames[18]="IPNegMin2";// = fNegDCAToPVMin2;
	fParameterNames[19]="DCADauMax2";// = fDCADaughterMax2;
	fParameterNames[20]="CPAMin2";// = fCPAMin2;
	fParameterNames[21]="DCAPVMin2";// = fDCAToPVMin2;
	fParameterNames[22]="COADauMin2";// = fCOADaughterMin2;
	fParameterNames[23]="DCAZDauMax2";// = fDCAZDaughterMax2;
	fParameterNames[24]="V0Win2";// = fWindowV02;
	fParameterNames[25]="CPA2to1";// = fCPAV01toV02;
	fParameterNames[26]="CPADB";// = fCPADibaryon;
	fParameterNames[27]="PProtonMax";// = fProtonPMax;
	fParameterNames[28]="PPionMax";// = fPionPMax;
	fParameterNames[29]="PKaonMax";// = fKaonPMax;
	fParameterNames[30]="LikeSsingDB";// = fLikeSignDB;
	fParameterNames[31]="PTrackMin";// = fTrackPMin;

  for (Int_t ivar=0; ivar<nVar1; ivar++) {
    fParametersTree->Branch(fParameterNames[ivar].Data(),&fParameters[ivar],Form("%s/f",fParameterNames[ivar].Data()));
  }


  const char* nameoutput2 = GetOutputSlot(2)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput2,"Candidates variables tree");
  Int_t nVar2 = 37;
  fCandidateVariables = new Float_t [nVar2];
  TString * fCandidateVariableNames = new TString[nVar2];

	fCandidateVariableNames[ 0]="antiL1";// = antiLambda1;
	fCandidateVariableNames[ 1]="PtPos1";// = MomPosOri1Pt;
	fCandidateVariableNames[ 2]="PtNeg1";// = MomNegOri1Pt;
	fCandidateVariableNames[ 3]="DcaPos1";// = posDCAPV;
	fCandidateVariableNames[ 4]="DcaNeg1";// = negDCAPV;
	fCandidateVariableNames[ 5]="IMLV01";// = InvMassLambdaV01;
	fCandidateVariableNames[ 6]="Dca1";// = dcaV01;
	fCandidateVariableNames[ 7]="Cpa1";// = cpaV01;
	fCandidateVariableNames[ 8]="Coa1";// = coaV01;
	fCandidateVariableNames[ 9]="DcaZ1";// = dcaZV01;
	fCandidateVariableNames[10]="xyDB";// = xyV01;
	fCandidateVariableNames[11]="rDB";// = rV01;
	fCandidateVariableNames[12]="ctauDB";// = ctauV01;
	fCandidateVariableNames[13]="Centrality";// = fCentrality;
	fCandidateVariableNames[14]="Cpa2";// = cpa2;
	fCandidateVariableNames[15]="IMLV02";// = InvMassLambdaV02;
	fCandidateVariableNames[16]="IMDB";// = vDB.M();
	fCandidateVariableNames[17]="CpaDB";// = cpaDB;
	fCandidateVariableNames[18]="OnFly2";// = static_cast<Float_t>(IsOnFly);
	fCandidateVariableNames[19]="DcaPos2";// = pDca;
	fCandidateVariableNames[20]="DcaNeg2";// = nDca;
	fCandidateVariableNames[21]="Coa2";// = coa2PosNeg;
	fCandidateVariableNames[22]="Dca2";// = dca2;
	fCandidateVariableNames[23]="xy2";// = xyV01toV02;
	fCandidateVariableNames[24]="r2";// = rV01toV02;
	fCandidateVariableNames[25]="antiL2";// = antiLambda2;
	fCandidateVariableNames[26]="IMLOri1";// = InvMassLambdaOri1;
	fCandidateVariableNames[27]="IMDBOri";// = vDBOri.M();
	fCandidateVariableNames[28]="ESD";// = static_cast<Float_t>(IsESDv01);
	fCandidateVariableNames[29]="Cpa2PV";// = cpa2ESD;
	fCandidateVariableNames[30]="Alpha";// = AlphaV01;
	fCandidateVariableNames[31]="PtArm";// = PtArmV01;
	fCandidateVariableNames[32]="IMG1";// = vCon.M;
	fCandidateVariableNames[33]="IMDBL";// = vDBL.M();
	fCandidateVariableNames[34]="IMDBOriL";// = vDBOriL.M();
	fCandidateVariableNames[35]="ctau2L";// = ctau1to2L;
	fCandidateVariableNames[36]="TT";// = static_cast<Float_t>(triggeType);

  for (Int_t ivar=0; ivar<nVar2; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;

}
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskNOmegaLPK::GetEsdV0Momentum(Int_t typeV0, Int_t id1, Int_t id2, Double_t ReturnV0[6], Double_t MomV0[3], Double_t PosV0[3], Double_t MomV0Pos[3], Double_t MomV0Neg[3]) {

	//------------------------------------------------------------------------------------------
	// version 1.10 (2016/02/20)
	// Check if the set of 2 tracks are assigned as ESD V0 and return ESD V0 momentum and position
	// Input:  V0 type (1:V01, 2:V02), Track IDs of 2 tracks
	// Return: kTRUE: Found
	//        kFALSE: Not found
	//      ReturnV0: 0:DCAxy, 1:CPA, 2:COA, 3:DCAz, 4:V0xy, 5:V0r
	//         V0Mom: Momentum of V0
	//         V0Pos: Position of V0
	//------------------------------------------------------------------------------------------

//	Bool_t recut = kTRUE;
	Bool_t recut = kFALSE;

	for ( Int_t i=0; i<6; i++ ) {
		ReturnV0[i] = -999.;
	}
 
	const AliESDVertex *esdPV = (AliESDVertex*)fESDEvent->GetPrimaryVertex();
  Double_t PosPV[3];
  esdPV->GetXYZ(PosPV);

	for ( Int_t i=0; i<3; i++ ) {
		MomV0[i] = -999.;
		PosV0[i] = -999.;
	}

	if ( id1==id2 ) return kFALSE;

	const Int_t nV0s = fESDEvent->GetNumberOfV0s();
	if (nV0s==0) {
		return kFALSE;
	}

	for ( Int_t iV0 = 0; iV0<nV0s; iV0++ ) {
		AliESDv0 *esdV0 = fESDEvent->GetV0(iV0);
		Bool_t IsOnFly = esdV0->GetOnFlyStatus();
		if ( IsOnFly ) continue;// select only offline reconstruction

		AliESDtrack *ptrk = fESDEvent->GetTrack(esdV0->GetPindex());
		AliESDtrack *ntrk = fESDEvent->GetTrack(esdV0->GetNindex());

		Double_t ptrkCharge = ptrk->GetSign();
		Double_t ntrkCharge = ntrk->GetSign();
		if        ( ptrkCharge>0 && ntrkCharge>0 ) {//(3,4)=(+,+)->continue
			return kFALSE;
		} else if ( ptrkCharge<0 && ntrkCharge<0 ) {//(3,4)=(-,-)->continue
			return kFALSE;
		} else if ( ptrkCharge>0 && ntrkCharge<0 ) {//(3,4)=(+,-)->ok
		} else if ( ptrkCharge<0 && ntrkCharge>0 ) {//(3,4)=(-,+)->(+,-)
			ptrk = fESDEvent->GetTrack(esdV0->GetNindex());
			ntrk = fESDEvent->GetTrack(esdV0->GetPindex());
		}

		Int_t pid = ptrk->GetID();
		Int_t nid = ntrk->GetID();

		if ( (id1==pid && id2==nid) || (id2==pid && id1==nid) ) {
			esdV0->PxPyPz(MomV0);
			esdV0->XvYvZv(PosV0);
			esdV0->GetPPxPyPz(MomV0Pos[0],MomV0Pos[1],MomV0Pos[2]);
			esdV0->GetNPxPyPz(MomV0Neg[0],MomV0Neg[1],MomV0Neg[2]);
			Double_t dca = esdV0->GetDcaV0Daughters();

			Double_t xPVtoV0  = PosV0[0] - PosPV[0];
			Double_t yPVtoV0  = PosV0[1] - PosPV[1];
			Double_t zPVtoV0  = PosV0[2] - PosPV[2];
			Double_t xyPVtoV0 = TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0);
			Double_t rPVtoV0  = TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0+zPVtoV0*zPVtoV0);
			Double_t MomV0Pt  = GetPtFromPxPyPz(MomV0);
			Double_t MomV0P   = GetPaFromPxPyPz(MomV0);
			Double_t cpa = (xPVtoV0*MomV0[0]+yPVtoV0*MomV0[1]+zPVtoV0*MomV0[2])/rPVtoV0/MomV0P;

			Double_t MomPosPt  = GetPtFromPxPyPz(MomV0Pos);
			Double_t MomNegPt  = GetPtFromPxPyPz(MomV0Neg);
			Double_t MomPosP   = GetPaFromPxPyPz(MomV0Pos);
			Double_t MomNegP   = GetPaFromPxPyPz(MomV0Neg);
			Double_t coaPosNeg = (MomV0Pos[0]*MomV0Neg[0]+MomV0Pos[1]*MomV0Neg[1]+MomV0Pos[2]*MomV0Neg[2])/MomPosP/MomNegP;

			Double_t pDca = ptrk->GetD(PosPV[0],PosPV[1],fBzkG);
			Double_t nDca = ntrk->GetD(PosPV[0],PosPV[1],fBzkG);

			Double_t xPos;
			Double_t xNeg;
			Double_t dcaCut;
			dcaCut = ptrk->GetDCA(ntrk,fBzkG,xPos,xNeg);

			AliExternalTrackParam posEtp(*ptrk);
			AliExternalTrackParam negEtp(*ntrk);

			Bool_t corrected = kFALSE;
			if ((posEtp.GetX() > 3.) && (xPos < 3.)) {
				corrected = kTRUE;//correct for the beam pipe material
			}
			if ((negEtp.GetX() > 3.) && (xNeg < 3.)) {
				corrected = kTRUE;//correct for the beam pipe material
			}
			if (corrected) {
				dcaCut = posEtp.GetDCA(&negEtp,fBzkG,xPos,xNeg);
			}

			Double_t r2V0 = PosV0[0]*PosV0[0] + PosV0[1]*PosV0[1];
			Double_t rV0  = TMath::Sqrt(r2V0);

			Double_t dcaV0ToPV = TMath::Abs(MomV0[0]/MomV0Pt*yPVtoV0-MomV0[1]/MomV0Pt*xPVtoV0);

			posEtp.PropagateTo(xPos,fBzkG);
			negEtp.PropagateTo(xNeg,fBzkG);
			Double_t PosPos[3], PosNeg[3];
			posEtp.GetXYZ(PosPos);
			negEtp.GetXYZ(PosNeg);
			Double_t dcaZ = PosPos[2] - PosNeg[2];

			ReturnV0[0] = dca;
			ReturnV0[1] = cpa;
			ReturnV0[2] = coaPosNeg;
			ReturnV0[3] = dcaZ;
			ReturnV0[4] = xyPVtoV0;
			ReturnV0[5] = rPVtoV0;

			if ( recut ) {
//				printf("##### ESD pDca: %10f(%10f-%10f), nDca: %10f(%10f-%10f)\n",pDca,fPosDCAToPVMin1,fFiducialVolMax1,nDca,fNegDCAToPVMin1,fFiducialVolMax1);

				if ( typeV0==1 ) {
					if ( TMath::Abs(pDca)<fPosDCAToPVMin1 ) return kFALSE;
					if ( TMath::Abs(pDca)>fFiducialVolMax1 ) return kFALSE;
					if ( TMath::Abs(nDca)<fNegDCAToPVMin1 ) return kFALSE;
					if ( TMath::Abs(nDca)>fFiducialVolMax1 ) return kFALSE;
					if ( dcaCut>fDCADaughterMax1 ) return kFALSE;
					if ( (xPos+xNeg)>2*fFiducialVolMax1 ) return kFALSE;
					if ( (xPos+xNeg)<2*fFiducialVolMin1 ) return kFALSE;
					if ( coaPosNeg<fCOADaughterMin1 ) return kFALSE;
					if ( rV0<fFiducialVolMin1) return kFALSE;
					if ( rV0>fFiducialVolMax1) return kFALSE;
					if ( cpa<fCPAMin1 ) return kFALSE;
					if ( dcaV0ToPV<fDCAToPVMin1 ) return kFALSE;
					if ( TMath::Abs(dcaZ)>fDCAZDaughterMax1 ) return kFALSE;
				} else if ( typeV0==2 ) {
					if ( TMath::Abs(pDca)<fPosDCAToPVMin2 ) return kFALSE;
					if ( TMath::Abs(pDca)>fFiducialVolMax2 ) return kFALSE;
					if ( TMath::Abs(nDca)<fNegDCAToPVMin2 ) return kFALSE;
					if ( TMath::Abs(nDca)>fFiducialVolMax2 ) return kFALSE;
					if ( dcaCut>fDCADaughterMax2 ) return kFALSE;
					if ( (xPos+xNeg)>2*fFiducialVolMax2 ) return kFALSE;
					if ( (xPos+xNeg)<2*fFiducialVolMin2 ) return kFALSE;
					if ( coaPosNeg<fCOADaughterMin2 ) return kFALSE;
					if ( rV0<fFiducialVolMin2) return kFALSE;
					if ( rV0>fFiducialVolMax2) return kFALSE;
					if ( cpa<fCPAMin2 ) return kFALSE;
					if ( dcaV0ToPV<fDCAToPVMin2 ) return kFALSE;
					if ( TMath::Abs(dcaZ)>fDCAZDaughterMax2 ) return kFALSE;
				} else return kFALSE;

			}

			return kTRUE;

		}

	}

	return kFALSE;

}
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskNOmegaLPK::GetSelfV0Momentum(Int_t typeV0, AliESDtrack *posTrk, AliESDtrack *negTrk,Double_t ReturnV0[6], Double_t MomV0[3], Double_t PosV0[3], Double_t MomPos[3], Double_t MomNeg[3]) {

	//------------------------------------------------------------------------------------------
	// version 1.10 (2016/02/20)
	// Copy from AliV0vertexer
	// reconstruct V0 from 2 tracks
	// Input:  V0 type (1:V01, 2:V02), 2 tracks of V0 positive and negative tracks
	// Return: kTRUE: Success
	//        kFALSE: Failed
	//      ReturnV0: 0:DCAxy, 1:CPA, 2:COA, 3:DCAz, 4:V0xy, 5:V0r
	//         V0Mom: Momentum of V0
	//         V0Pos: Position of V0
	//------------------------------------------------------------------------------------------

	for ( Int_t i=0; i<6; i++ ) {
		ReturnV0[i] = -999.;
	}
 
	AliESDEvent *fESDEvent = dynamic_cast<AliESDEvent*>(fInputEvent);

  const AliESDVertex *esdPV = (AliESDVertex*)fESDEvent->GetPrimaryVertex();
  Double_t PosPV[3];
  esdPV->GetXYZ(PosPV);

	Int_t eventNumber = fESDEvent->GetEventNumberInFile();

	for (Int_t i=0; i<3; i++) {
	  MomV0[i] = -999.;
	  PosV0[i] = -999.;
	  MomPos[i] = -999.;
	  MomNeg[i] = -999.;
	}

	//------------------------------------------------------------------------------------------
	// Track cuts
	//------------------------------------------------------------------------------------------

	// positive track cut
	Double_t posDCAPV = posTrk->GetD(PosPV[0],PosPV[1],fBzkG);
	if (typeV0==1) {
		if (TMath::Abs(posDCAPV)<fPosDCAToPVMin1) return kFALSE;
		if (TMath::Abs(posDCAPV)>fFiducialVolMax1) return kFALSE;
	} else if (typeV0==2) {
		if (TMath::Abs(posDCAPV)<fPosDCAToPVMin2) return kFALSE;
		if (TMath::Abs(posDCAPV)>fFiducialVolMax2) return kFALSE;
	} else return kFALSE;

	// negative track cut
	Double_t negDCAPV = negTrk->GetD(PosPV[0],PosPV[1],fBzkG);
	if (typeV0==1) {
		if (TMath::Abs(negDCAPV)<fNegDCAToPVMin1) return kFALSE;
		if (TMath::Abs(negDCAPV)>fFiducialVolMax1) return kFALSE;
	} else if (typeV0==2) {
		if (TMath::Abs(negDCAPV)<fNegDCAToPVMin2) return kFALSE;
		if (TMath::Abs(negDCAPV)>fFiducialVolMax2) return kFALSE;
	} else return kFALSE;
//	printf("#### Self pDca: %10f(%10f-%10f), nDca: %10f(%10f-%10f)\n",posTrk->GetD(PosPV[0],PosPV[1],fBzkG),fPosDCAToPVMin1,fFiducialVolMax1,negTrk->GetD(PosPV[0],PosPV[1],fBzkG),fNegDCAToPVMin1,fFiducialVolMax1);


	//------------------------------------------------------------------------------------------
	// Reconstruct V0 and cut V0
	//------------------------------------------------------------------------------------------

	// get DCA, propagate tracks, set V0
	Double_t xPos;
	Double_t xNeg;
	Double_t dcaV0Daughters;
	dcaV0Daughters = posTrk->GetDCA(negTrk,fBzkG,xPos,xNeg);
	ReturnV0[0] = dcaV0Daughters;
	if (typeV0==1) {
		if (dcaV0Daughters > fDCADaughterMax1) return kFALSE;
		if ((xPos+xNeg) > 2*fFiducialVolMax1) return kFALSE;
		if ((xPos+xNeg) < 2*fFiducialVolMin1) return kFALSE;
	} else if (typeV0==2) {
		if (dcaV0Daughters > fDCADaughterMax2) return kFALSE;
		if ((xPos+xNeg) > 2*fFiducialVolMax2) return kFALSE;
		if ((xPos+xNeg) < 2*fFiducialVolMin2) return kFALSE;
	} else return kFALSE;

	AliExternalTrackParam posEtp(*posTrk);
	AliExternalTrackParam negEtp(*negTrk);

	Bool_t corrected = kFALSE;
	if ((posEtp.GetX() > 3.) && (xPos < 3.)) {
		corrected = kTRUE;//correct for the beam pipe material
	}
	if ((negEtp.GetX() > 3.) && (xNeg < 3.)) {
		corrected = kTRUE;//correct for the beam pipe material
	}
	if (corrected) {
		dcaV0Daughters = posEtp.GetDCA(&negEtp,fBzkG,xPos,xNeg);
		ReturnV0[0] = dcaV0Daughters;
		if (typeV0==1) {
			if (dcaV0Daughters > fDCADaughterMax1) return kFALSE;
			if ((xPos+xNeg) > 2*fFiducialVolMax1) return kFALSE;
			if ((xPos+xNeg) < 2*fFiducialVolMin1) return kFALSE;
		} else if (typeV0==2) {
			if (dcaV0Daughters > fDCADaughterMax2) return kFALSE;
			if ((xPos+xNeg) > 2*fFiducialVolMax2) return kFALSE;
			if ((xPos+xNeg) < 2*fFiducialVolMin2) return kFALSE;
		} else return kFALSE;
	}

	// get the momentum and position of tracks propagated to V0 (DCA of 2 tracks)
	posEtp.PropagateTo(xPos,fBzkG);
	negEtp.PropagateTo(xNeg,fBzkG);

	Double_t chargeV0 = posEtp.GetSign() + negEtp.GetSign();
	posEtp.GetPxPyPz(MomPos);
	negEtp.GetPxPyPz(MomNeg);
	Double_t MomPosP  = GetPaFromPxPyPz(MomPos);
	Double_t MomPosPt = GetPtFromPxPyPz(MomPos);
	Double_t MomNegP  = GetPaFromPxPyPz(MomNeg);
	Double_t MomNegPt = GetPtFromPxPyPz(MomNeg);
	Double_t PosPos[3], PosNeg[3];
	posEtp.GetXYZ(PosPos);
	negEtp.GetXYZ(PosNeg);
	Double_t dcaZ = PosPos[2] - PosNeg[2];
	ReturnV0[3] = dcaZ;
	if (typeV0==1) {
		if ( TMath::Abs(dcaZ)>fDCAZDaughterMax1 ) return kFALSE;
	} else if (typeV0==2) {
		if ( TMath::Abs(dcaZ)>fDCAZDaughterMax2 ) return kFALSE;
	} else return kFALSE;

	// get V0 momentum and position
	for (Int_t i=0; i<3; i++) {
	  MomV0[i] = MomPos[i] + MomNeg[i];
	  PosV0[i] = (PosPos[i] + PosNeg[i]) / 2.;
	}
	Double_t MomV0P  = GetPaFromPxPyPz(MomV0);
	Double_t MomV0Pt = GetPtFromPxPyPz(MomV0);

	Double_t coaV0Daughters   = (MomPos[0]*MomNeg[0]+MomPos[1]*MomNeg[1]+MomPos[2]*MomNeg[2])
	                            /MomPosP/MomNegP;
	Double_t coaV0Daughters2D = (MomPos[0]*MomNeg[0]+MomPos[1]*MomNeg[1])
	                            /MomPosPt/MomNegPt;
	ReturnV0[2] = coaV0Daughters;
	if (typeV0==1) {
		if (coaV0Daughters2D < fCOADaughterMin1) return kFALSE;
	} else if (typeV0==2) {
		if (coaV0Daughters2D < fCOADaughterMin2) return kFALSE;
	} else return kFALSE;

	Double_t r2V0 = PosV0[0]*PosV0[0] + PosV0[1]*PosV0[1];
	Double_t rV0  = TMath::Sqrt(r2V0);
	if (typeV0==1) {
		if (rV0 < fFiducialVolMin1) return kFALSE;
		if (rV0 > fFiducialVolMax1) return kFALSE;
	} else if (typeV0==2) {
		if (rV0 < fFiducialVolMin2) return kFALSE;
		if (rV0 > fFiducialVolMax2) return kFALSE;
	} else return kFALSE;

	Double_t xPVtoV0  = PosV0[0] - PosPV[0];
	Double_t yPVtoV0  = PosV0[1] - PosPV[1];
	Double_t zPVtoV0  = PosV0[2] - PosPV[2];
	Double_t cpaPVtoV0 = (xPVtoV0*MomV0[0]+yPVtoV0*MomV0[1]+zPVtoV0*MomV0[2])/MomV0P
	                      /TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0+zPVtoV0*zPVtoV0);
	ReturnV0[1] = cpaPVtoV0;
	if (typeV0==1) {
		if (cpaPVtoV0 < fCPAMin1) return kFALSE;
	} else if (typeV0==2) {
		if (cpaPVtoV0 < fCPAMin2) return kFALSE;
	} else return kFALSE;

	Double_t xyPVtoV0 = TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0);
	ReturnV0[4] = xyPVtoV0;
	Double_t rPVtoV0  = TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0+zPVtoV0*zPVtoV0);
	ReturnV0[5] = rPVtoV0;

	Double_t dcaV0ToPV = TMath::Abs(MomV0[0]/MomV0Pt*yPVtoV0-MomV0[1]/MomV0Pt*xPVtoV0);
	if (typeV0==1) {
  	if (dcaV0ToPV < fDCAToPVMin1) return kFALSE;
	} else if (typeV0==2) {
  	if (dcaV0ToPV < fDCAToPVMin2) return kFALSE;
	} else return kFALSE;
//	printf("1,2:posDCA    , 3,4:posDCA    , 5:dcaV0Daugh, 6,7:xPos+xNeg , 8:coaV0Daugh, 9,10:rV0       , 11:cpaPVtoV0 , 12:dcaV0ToPV , 13:dcaZ     \n",TMath::Abs(posDCAPV),TMath::Abs(negDCAPV),dcaV0Daughters,(xPos+xNeg),coaV0Daughters2D,rV0,cpaPVtoV0,dcaV0ToPV,TMath::Abs(dcaZ));
//	printf("1,2:%10f, 3,4:%10f, 5:%10f, 6,7:%10f, 8:%10f, 9,10:%10f, 11:%10f, 12:%10f, 13:%10f\n",TMath::Abs(posDCAPV),TMath::Abs(negDCAPV),dcaV0Daughters,(xPos+xNeg),coaV0Daughters2D,rV0,cpaPVtoV0,dcaV0ToPV,TMath::Abs(dcaZ));

	return kTRUE; 


}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskNOmegaLPK::InvMassLambda(Double_t MomPos[3], Double_t MomNeg[3], AliESDtrack *pos, AliESDtrack *neg, Double_t v0Return[1]) {

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
//		printf("### TOF:  POS Pr:%10f, POS Pi:%10f,\tNEG Pr:%10f, NEG Pi:%10f\n",posTOFProton,posTOFPion,negTOFProton,negTOFPion);
//	}
//	printf("!!! THEREFORE: %f (1:Lambda, -1:Anti-Lambda, 0:Other)\n",v0Return[0]);
	

  return vLambda.M();

}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskNOmegaLPK::InvMassLambdaStar(Double_t MomPos[3], Double_t MomNeg[3], AliESDtrack *pos, AliESDtrack *neg, Double_t v0Return[1]) {

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
//__________________________________________________________________________
