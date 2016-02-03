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
#include "AliPIDCombined.h"
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

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskNOmegaLPK.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNOmegaLPK)

//______________________________________________________________________________________________________
AliAnalysisTaskNOmegaLPK::AliAnalysisTaskNOmegaLPK() : 
	AliAnalysisTaskSE(),
	fESDtrackCutsV0(0),
//	fESDCutsV0(0),
	fESDtrackCuts(0),
  fUseMCInfo(kFALSE),
  fPIDResponse(0),
  fPIDCombined(0),
  fIsEventSelected(kFALSE),
  fMixedEvent(kFALSE),
  fVariablesTree1(0),
  fVariablesTree2(0),
  fCandidateVariables1(),
  fCandidateVariables2(),
	fVtx1(0),
  fBzkG(0),
	fSigmaTPC(0),
	fClustersTPC(0),
	fSigmaTOF(0),
	fPseudoRap(0),
  fCentrality(0),
  fCountMatch(0),
  fCountOnlySelf(0),
  fCountOnlyV0(0),
  fCountLambda(0),
  fCountAntiLambda(0),
	fCountLambdaSelf(0),
	fCountAntiLambdaSelf(0),
	fCountMatchLambda(0),
	fCountMatchAntiLambda(0),
	fCountEvent(0),
	fIsDCA(0),
	fNAll(0),
	fFiducialVolMin1(2.0),
	fFiducialVolMax1(200.),
	fPosDCAToPVMin1(2.0),
	fNegDCAToPVMin1(2.0),
	fDCADaughterMax1(1.0),
	fCPAMin1(-1.),
	fDCAToPVMin1(0.0),
	fCOADaughterMin1(-1.0),
	fDCAZDaughterMax1(2.0),
	fWindowLambda(0.0045),
	fCPADibaryon(0.99875)
{
}
//______________________________________________________________________________________________________
AliAnalysisTaskNOmegaLPK::AliAnalysisTaskNOmegaLPK(const Char_t* name) :
  AliAnalysisTaskSE(name),
	fESDtrackCutsV0(0),
//	fESDCutsV0(0),
	fESDtrackCuts(0),
  fUseMCInfo(kFALSE),
  fPIDResponse(0),
  fPIDCombined(0),
  fIsEventSelected(kFALSE),
  fMixedEvent(kFALSE),
  fVariablesTree1(0),
  fVariablesTree2(0),
  fCandidateVariables1(),
  fCandidateVariables2(),
  fVtx1(0),
  fBzkG(0),
	fSigmaTPC(0),
	fClustersTPC(0),
	fSigmaTOF(0),
	fPseudoRap(0),
  fCentrality(0),
  fCountMatch(0),
  fCountOnlySelf(0),
  fCountOnlyV0(0),
  fCountLambda(0),
  fCountAntiLambda(0),
	fCountLambdaSelf(0),
	fCountAntiLambdaSelf(0),
	fCountMatchLambda(0),
	fCountMatchAntiLambda(0),
	fCountEvent(0),
	fIsDCA(0),
	fNAll(0),
	fFiducialVolMin1(2.0),
	fFiducialVolMax1(200.),
	fPosDCAToPVMin1(2.0),
	fNegDCAToPVMin1(2.0),
	fDCADaughterMax1(1.0),
	fCPAMin1(-1.),
	fDCAToPVMin1(0.0),
	fCOADaughterMin1(-1.0),
	fDCAZDaughterMax1(2.0),
	fWindowLambda(0.0045),
	fCPADibaryon(0.99875)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskNOmegaLPK","Calling Constructor");

	DefineOutput(1,TTree::Class());  //My private output
	DefineOutput(2,TTree::Class());  //My private output


	//V0 cuts
	fESDtrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0","AliESDtrackCutsV0");
	fESDtrackCutsV0->SetMinNClustersTPC(80);
	fESDtrackCutsV0->SetAcceptKinkDaughters(kFALSE);
	fESDtrackCutsV0->SetMaxChi2PerClusterTPC(5);
	fESDtrackCutsV0->SetRequireTPCRefit(kTRUE);
	fESDtrackCutsV0->SetEtaRange(-0.9,0.9);
	fESDtrackCutsV0->SetPtRange(0.2,1.5);
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
AliAnalysisTaskNOmegaLPK::~AliAnalysisTaskNOmegaLPK() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskNOmegaLPK","Calling Destructor");
  
  if (fPIDResponse) {
    delete  fPIDResponse;
  }

  if (fPIDCombined) {
    delete  fPIDCombined;
  }

  if (fVariablesTree1) {
    delete fVariablesTree1;
    fVariablesTree1 = 0;
  }

  if (fVariablesTree2) {
    delete fVariablesTree2;
    fVariablesTree2 = 0;
  }

	if(fESDtrackCutsV0) delete fESDtrackCutsV0;
//	if(fESDCutsV0) delete fESDCutsV0;
	if(fESDtrackCuts) delete fESDtrackCuts;

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
	AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
	if (!esdEvent) {
 	  AliError("NO EVENT FOUND!");
 		return;
  }

  //------------------------------------------------
  // First check if the event has proper vertex and B
  //------------------------------------------------
	AliESDVertex *esdTrkPV = (AliESDVertex*)esdEvent->GetPrimaryVertexTracks();
	Double_t PosTrkPV[3];
	esdTrkPV->GetXYZ(PosTrkPV);

	AliESDVertex *esdPV = (AliESDVertex*)esdEvent->GetPrimaryVertex();
	Double_t posPV[3];
	esdPV->GetXYZ(posPV);

	if (TMath::Abs(PosTrkPV[2])>10.) {
 	  AliError("Z of primary vertex is not within +-10cm !");
		return;
	}

	Int_t runNumber = 0;
  runNumber = (InputEvent())->GetRunNumber();

	fBzkG = (Double_t)esdEvent->GetMagneticField(); 
	AliKFParticle::SetField(fBzkG);
  if (TMath::Abs(fBzkG)<0.001) {
 	  AliError("No magnet field !");
    return;
  }

	//------------------------------------------------
  // MC analysis setting
  //------------------------------------------------
  TClonesArray *mcArray = 0;
//  AliAODMCHeader *mcHeader=0;
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
  Int_t nSelectedAnal = 0;

	//------------------------------------------------
  // Main analysis done in this function
  //------------------------------------------------
	MakeAnalysis(mcArray, nSelectedAnal, esdEvent);

	PostData(1,fVariablesTree1);
	PostData(2,fVariablesTree2);

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
	PostData(1,fVariablesTree1);
	PostData(2,fVariablesTree2);

	// I don't want to use the PID through the cut object, 
	// but I will use the PID response directly!!!
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  // Setting properties of PID
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);

  return;
}
//______________________________________________________________________________________________________
void AliAnalysisTaskNOmegaLPK::MakeAnalysis(TClonesArray *mcArray,Int_t &nSelectedAnal,AliESDEvent *esdEvent)
{

  //------------------------------------------------------------------------------------------
  // version 4-0-0 (2016/02/01) from 3-3-1
	// revised 2016/02/03 for train
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

	// For track cut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	fSigmaTPC          = 3.;//TPC PIDcut sigma
	fClustersTPC       = 80;//TPC number of clusters
	fSigmaTOF          = 3.;//TOF PIDcut sigma
	fPseudoRap         = 0.9;//PseudoRapidity
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// For V01 cut (Inner 2 tracks)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
///	Double_t fFiducialVolMin1  = 2.;//0.2;// Min radius of the fiducial volume
///	Double_t fFiducialVolMax1  = 200.;// Max radius of the fiducial volume
///	Double_t fPosDCAToPVMin1   = 2.;//0.05;// Min length of impact parameter for the positive track
///	Double_t fNegDCAToPVMin1   = 2.;//0.05;// Min length of impact parameter for the negative track
///	Double_t fDCADaughterMax1  = 1.0;//2.0;// Max DCA between the daughter tracks
///	Double_t fCPAMin1          = -1.0//0.0;//0.5;// Min cosine of V0's pointing angle
///	Double_t fDCAToPVMin1      = 0.0;// Min DCA V0 to PV
///	Double_t fCOADaughterMin1  = 0.0;// Min cosine between the daughter tracks
///	Double_t fDCAZDaughterMax1 = 2.0;// Max DCAZ V0 to PV
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// For dibaryon cut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
///	Double_t fWindowLambda = 0.0045;//Mass window cut for Lambda
///	Double_t fCPADibaryon  = 0.99875;//0.9;//Mass window cut for dibaryon
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
//  Double_t mXi1530PDG   = TDatabasePDG::Instance()->GetParticle(3314)->Mass();//1.535000
//  Double_t mOmegaPDG    = TDatabasePDG::Instance()->GetParticle(3334)->Mass();//1.672450

  const Int_t nTracks = esdEvent->GetNumberOfTracks();
  if (nTracks==0) {
    return;
  }

	const Int_t nV0s = esdEvent->GetNumberOfV0s();

  // Initialization of parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Int_t posID[nV0s];
	Int_t negID[nV0s];
	Double_t infoV0Tracks[nV0s][15];

	Int_t nLambda     = 0;
	Int_t nAntiLambda = 0;

	Double_t cutV01Parameters[9];
	cutV01Parameters[0] = fFiducialVolMin1;
	cutV01Parameters[1] = fFiducialVolMax1;
	cutV01Parameters[2] = fPosDCAToPVMin1;
	cutV01Parameters[3] = fNegDCAToPVMin1;
	cutV01Parameters[4] = fDCADaughterMax1;
	cutV01Parameters[5] = fCPAMin1;
	cutV01Parameters[6] = fDCAToPVMin1;
	cutV01Parameters[7] = fCOADaughterMin1;
	cutV01Parameters[8] = fDCAZDaughterMax1;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Get or calculate constant for this event ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const AliESDVertex *esdTrkPV = (AliESDVertex*)esdEvent->GetPrimaryVertexTracks();
	Double_t PosTrkPV[3];
	esdTrkPV->GetXYZ(PosTrkPV);

	const AliESDVertex *esdPV = (AliESDVertex*)esdEvent->GetPrimaryVertex();
	Double_t PosPV[3];
	esdPV->GetXYZ(PosPV);

	Double_t vPVX = PosPV[0];
	Double_t vPVY = PosPV[1];
	Double_t vPVZ = PosPV[2];

	Int_t runNumber;
	AliCentrality *cent = esdEvent->GetCentrality();
	fCentrality = cent->GetCentralityClass10("V0M");//GetCentralityPercentile????
//	fCentrality = cent->GetCentralityPercentile("V0M");//GetCentralityPercentile????
	runNumber = esdEvent->GetRunNumber();

	Int_t eventNumber = esdEvent->GetEventNumberInFile();
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // START ANALYSIS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//------------------------------------------------------------------------------------------
	// Cut for event
	//------------------------------------------------------------------------------------------

//	if (fCountEvent<setStartNumber) return;
//	if (fCountEvent>5) return;

	//------------------------------------------------------------------------------------------
	// Track loop 1 (To find V01) (START)
	//------------------------------------------------------------------------------------------

	Int_t nV01Test = 0;

	for (Int_t iTrk1=0; iTrk1<nTracks; iTrk1++) {
		AliESDtrack *trk1=esdEvent->GetTrack(iTrk1);
		Double_t trk1Charge = trk1->GetSign();

		if( !fESDtrackCuts->AcceptTrack(trk1) ) continue;

		// Pre-track cut (TPCrefit, TPCclusters, TPCsigma, RejectKink, PseudoRapidity)
		if( !PreTrackCut(trk1) ) continue;

		// Impact parameter cut
		Double_t trk1DCAPV = trk1->GetD(PosPV[0],PosPV[1],fBzkG);
		if ( TMath::Abs(trk1DCAPV)<fPosDCAToPVMin1 && TMath::Abs(trk1DCAPV)<fNegDCAToPVMin1 ) continue;

	//------------------------------------------------------------------------------------------
	// Track loop 2 (To find V01) (START)
	//------------------------------------------------------------------------------------------

		for (Int_t iTrk2=iTrk1+1; iTrk2<nTracks; iTrk2++) {
			AliESDtrack *trk2=esdEvent->GetTrack(iTrk2);
			Double_t trk2Charge = trk2->GetSign();

			if( !fESDtrackCuts->AcceptTrack(trk2) ) continue;

			// Pre-track cut (TPCrefit, TPCclusters, TPCsigma, RejectKink, PseudoRapidity)
			if( !PreTrackCut(trk2) ) continue;

			// Impact parameter cut
			Double_t trk2DCAPV = trk2->GetD(PosPV[0],PosPV[1],fBzkG);
			if ( TMath::Abs(trk2DCAPV)<fPosDCAToPVMin1 && TMath::Abs(trk2DCAPV)<fNegDCAToPVMin1 ) continue;

	//------------------------------------------------------------------------------------------
	// Preparation to reconstruct V01
	//------------------------------------------------------------------------------------------

			// Ordering positive -> negative track 
			AliESDtrack *pos1Trk, *neg1Trk;
			if        ( trk1Charge>0 && trk2Charge>0 ) {//(1,2)=(+,+)->continue
				continue;
			} else if ( trk1Charge<0 && trk2Charge<0 ) {//(1,2)=(-,-)->continue
				continue;
			} else if ( trk1Charge>0 && trk2Charge<0 ) {//(1,2)=(+,-)->ok
				pos1Trk=esdEvent->GetTrack(iTrk1);
				neg1Trk=esdEvent->GetTrack(iTrk2);
			} else if ( trk1Charge<0 && trk2Charge>0 ) {//(1,2)=(-,+)->(+,-)
				pos1Trk=esdEvent->GetTrack(iTrk2);
				neg1Trk=esdEvent->GetTrack(iTrk1);
			}

			if( !fESDtrackCutsV0->AcceptTrack(neg1Trk) ) continue;

			Double_t pos1DCAPV = pos1Trk->GetD(PosPV[0],PosPV[1],fBzkG);
			Double_t neg1DCAPV = neg1Trk->GetD(PosPV[0],PosPV[1],fBzkG);
			if ( TMath::Abs(pos1DCAPV)<fPosDCAToPVMin1 ) continue;
			if ( TMath::Abs(neg1DCAPV)<fNegDCAToPVMin1 ) continue;

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
//			Double_t InvMassLambdaStarOri1 = -9.;
			Double_t L1ReturnOri[1];
			L1ReturnOri[0] = 0.;
//			Double_t LS1ReturnOri[1];
//			LS1ReturnOri[0] = 0.;
			InvMassLambdaOri1     = InvMassLambda(    MomPosOri1,MomNegOri1,pos1Trk,neg1Trk,L1ReturnOri);
//			InvMassLambdaStarOri1 = InvMassLambdaStar(MomPosOri1,MomNegOri1,pos1Trk,neg1Trk,LS1ReturnOri);

/*
			Double_t antiLambdaOri1 = 0.;
			if ( L1ReturnOri[0]==1. ) {
				antiLambdaOri1 = -1.;//1: Anti, -1: Not Anti
				nAntiLambda++;
				fCountAntiLambda++;
			} else if ( L1ReturnOri[0]==-1. ) {
				antiLambdaOri1 = 1.;//1: Anti, -1: Not Anti
				nLambda++;
				fCountLambda++;
			}
			Double_t antiLambdaStarOri1 = 0.;
			if ( LS1ReturnOri[0]==1. ) {
				antiLambdaStarOri1 = -1.;//1: Anti, -1: Not Anti
			} else if ( LS1ReturnOri[0]==-1. ) {
				antiLambdaStarOri1 = 1.;//1: Anti, -1: Not Anti
			}
*/

	//------------------------------------------------------------------------------------------
	// Get V01 Momentumas
	//------------------------------------------------------------------------------------------

//			Double_t ReturnV01[6],
			Double_t ReturnV01Self[6];//, ReturnV01Cross[6];//0:DCAxy, 1:CPA, 2:COA, 3:DCAz, 4:V0xy, 5:V0r
//			Double_t MomV01[3], PosV01[3], MomV01Pos[3], MomV01Neg[3];
			Double_t MomV01Self[3], PosV01Self[3], MomV01PosSelf[3], MomV01NegSelf[3];
//			Double_t MomV01Cross[3], PosV01Cross[3], MomV01PosCross[3], MomV01NegCross[3];
//			Bool_t IsV01,
			Bool_t IsV01Self;//, IsV01Cross;
//			IsV01 = GetEsdV0Momentum(esdEvent,pos1ID,neg1ID,cutV01Parameters,ReturnV01,MomV01,PosV01,MomV01Pos,MomV01Neg);
			IsV01Self = GetSelfV0Momentum(pos1Trk,neg1Trk,cutV01Parameters,ReturnV01Self,MomV01Self,PosV01Self,MomV01PosSelf,MomV01NegSelf);
//			IsV01Cross = GetSelfV0CrossMomentum(pos1Trk,neg1Trk,cutV01Parameters,ReturnV01Cross,MomV01Cross,PosV01Cross,MomV01PosCross,MomV01NegCross);

//			Double_t InvMassLambdaV01 = -9.;
			Double_t InvMassLambdaV01Self = -9.;
//			Double_t InvMassLambdaV01Cross = -9.;
//			Double_t InvMassLambdaStarV01 = -9.;
//			Double_t InvMassLambdaStarV01Self = -9.;
//			Double_t InvMassLambdaStarV01Cross = -9.;
//			Double_t L1ReturnV0[1],
				Double_t L1ReturnV0Self[1];//, L1ReturnV0Cross[1];
//			L1ReturnV0[0] = 0.;
			L1ReturnV0Self[0] = 0.;
//			L1ReturnV0Cross[0] = 0.;
//			Double_t LS1ReturnV0[1],
//			Double_t LS1ReturnV0Self[1];//, LS1ReturnV0Cross[1];
//			LS1ReturnV0[0] = 0.;
//			LS1ReturnV0Self[0] = 0.;
//			LS1ReturnV0Cross[0] = 0.;
//			Double_t dcaV01       = -999.;
			Double_t dcaV01Self   = -999.;
//			Double_t dcaV01Cross  = -999.;
//			Double_t cpaV01       = -999.;
			Double_t cpaV01Self   = -999.;
//			Double_t cpaV01Cross  = -999.;
//			Double_t coaV01       = -999.;
			Double_t coaV01Self   = -999.;
//			Double_t coaV01Cross  = -999.;
//			Double_t dcaZV01      = -999.;
			Double_t dcaZV01Self  = -999.;
//			Double_t dcaZV01Cross = -999.;
//			Double_t xyV01        = -999.;
			Double_t xyV01Self    = -999.;
//			Double_t xyV01Cross   = -999.;
//			Double_t rV01         = -999.;
			Double_t rV01Self     = -999.;
//			Double_t rV01Cross    = -999.;
//			Double_t ctauV01      = -999.;
			Double_t ctauV01Self  = -999.;
//			Double_t ctauV01Cross = -999.;

/*
			if ( IsV01 ) {
				nV01Test++;
				InvMassLambdaV01     = InvMassLambda(    MomV01Pos,MomV01Neg,pos1Trk,neg1Trk,L1ReturnV0);
				InvMassLambdaStarV01 = InvMassLambdaStar(MomV01Pos,MomV01Neg,pos1Trk,neg1Trk,LS1ReturnV0);
				dcaV01  = ReturnV01[0];
				cpaV01  = ReturnV01[1];
				coaV01  = ReturnV01[2];
				dcaZV01 = ReturnV01[3];
				xyV01   = ReturnV01[4];
				rV01    = ReturnV01[5];
				ctauV01 = ReturnV01[5]*mLambdaPDG/GetPaFromPxPyPz(MomV01);
			}
*/
			if ( IsV01Self ) {
				InvMassLambdaV01Self     = InvMassLambda(    MomV01PosSelf,MomV01NegSelf,pos1Trk,neg1Trk,L1ReturnV0Self);
//				InvMassLambdaStarV01Self = InvMassLambdaStar(MomV01PosSelf,MomV01NegSelf,pos1Trk,neg1Trk,LS1ReturnV0Self);
				dcaV01Self  = ReturnV01Self[0];
				cpaV01Self  = ReturnV01Self[1];
				coaV01Self  = ReturnV01Self[2];
				dcaZV01Self = ReturnV01Self[3];
				xyV01Self   = ReturnV01Self[4];
				rV01Self    = ReturnV01Self[5];
				ctauV01Self = ReturnV01Self[5]*mLambdaPDG/GetPaFromPxPyPz(MomV01Self);
			}
/*
			if ( IsV01Cross ) {
				InvMassLambdaV01Cross     = InvMassLambda(    MomV01PosCross,MomV01NegCross,pos1Trk,neg1Trk,L1ReturnV0Cross);
				InvMassLambdaStarV01Cross = InvMassLambdaStar(MomV01PosCross,MomV01NegCross,pos1Trk,neg1Trk,LS1ReturnV0Cross);
				dcaV01Cross  = ReturnV01Cross[0];
				cpaV01Cross  = ReturnV01Cross[1];
				coaV01Cross  = ReturnV01Cross[2];
				dcaZV01Cross = ReturnV01Cross[3];
				xyV01Cross   = ReturnV01Cross[4];
				rV01Cross    = ReturnV01Cross[5];
				ctauV01Cross = ReturnV01Cross[5]*mLambdaPDG/GetPaFromPxPyPz(MomV01Cross);
			}
*/

	//------------------------------------------------------------------------------------------
	// V0 Cut (V01)
	//------------------------------------------------------------------------------------------

//			if (InvMassLambdaV01<1. && InvMassLambdaStarV01<1. ) {
//			if (InvMassLambdaV01Self<1. && InvMassLambdaStarV01Self<1. ) {
			if (InvMassLambdaV01Self<1.) {
//			if (InvMassLambdaV01Cross<1. && InvMassLambdaStarV01Cross<1. ) {
						 continue;
//					}
				}
//			}
//			if ( !IsV01 && !IsV01Self && !IsV01Cross ) continue;
			if ( !IsV01Self ) continue;
//			if ( IsV01Cross ) {
//			printf("### (ESD,Self,Cross) = (%d,%d,%d)\n",static_cast<Int_t>(IsV01),static_cast<Int_t>(IsV01Self),static_cast<Int_t>(IsV01Cross));
//			}

	//------------------------------------------------------------------------------------------
	// Count Lambda (V01)
	//------------------------------------------------------------------------------------------

			Double_t antiLambdaOri1 = 0.;
			if ( L1ReturnOri[0]==1. ) {
				antiLambdaOri1 = -1.;//1: Anti, -1: Not Anti
				nAntiLambda++;
				fCountAntiLambda++;
			} else if ( L1ReturnOri[0]==-1. ) {
				antiLambdaOri1 = 1.;//1: Anti, -1: Not Anti
				nLambda++;
				fCountLambda++;
			}
//			Double_t antiLambdaStarOri1 = 0.;
//			if ( LS1ReturnOri[0]==1. ) {
//				antiLambdaStarOri1 = -1.;//1: Anti, -1: Not Anti
//			} else if ( LS1ReturnOri[0]==-1. ) {
//				antiLambdaStarOri1 = 1.;//1: Anti, -1: Not Anti
//			}

	//------------------------------------------------------------------------------------------
	// Mixed Event (Test)
	//------------------------------------------------------------------------------------------

/*
	Int_t maxEvts       = 40;
	Int_t minNTracks    = 10;
	Int_t nMultBins     = 9;
	Double_t multbins[] = {0,6,12,19,28,39,49,59,71,82};
	Int_t nZvtxBins     = 7;
	Double_t zvtxbins[] = {-10,-5,-2.5,0,2.5,5,10};
  if(fMixedEvent){
		AliEventPoolManager *fPoolMgr = new AliEventPoolManager(maxEvts,minNTracks,nMultBins,multbins,nZvtxBins,zvtxbins);
// 		AliEventPool *fPool= fPoolMgr->GetEventPool(fCentrality, PosPV[2]);
 		AliEventPool *fPool= fPoolMgr->GetEventPool(nTracks, PosPV[2]);
    if(!fPool){
//      AliInfo(Form("No pool found for Event: multiplicity = %f, zVtx = %f cm", fCentrality, PosPV[2]));
      AliInfo(Form("No pool found for Event: nTracks = %f, zVtx = %f cm", nTracks, PosPV[2]));
      return;
    }
  }
*/

	//------------------------------------------------------------------------------------------
	// V0 loop (To find V02) (START)
	//------------------------------------------------------------------------------------------

	for (Int_t iV02=0; iV02<nV0s; iV02++) {
		AliESDv0 *esdV02 = esdEvent->GetV0(iV02);
    Bool_t IsOnFly = esdV02->GetOnFlyStatus();

		AliESDtrack *ptrk = esdEvent->GetTrack(esdV02->GetPindex());
		AliESDtrack *ntrk = esdEvent->GetTrack(esdV02->GetNindex());
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
			ptrk = esdEvent->GetTrack(esdV02->GetNindex());
			ntrk = esdEvent->GetTrack(esdV02->GetPindex());
		}

		Int_t pos2ID = ptrk->GetID();
		Int_t neg2ID = ntrk->GetID();

		if (pos1ID==pos2ID) continue;
		if (neg1ID==neg2ID) continue;

		Double_t MomV02[3], PosV02[3];
		Double_t MomV02Pos[3], MomV02Neg[3];
		esdV02->PxPyPz(MomV02);
		esdV02->XvYvZv(PosV02);
		esdV02->GetPPxPyPz(MomV02Pos[0],MomV02Pos[1],MomV02Pos[2]);
		esdV02->GetNPxPyPz(MomV02Neg[0],MomV02Neg[1],MomV02Neg[2]);
		Double_t dca2 = esdV02->GetDcaV0Daughters();
		if ( dca2>1. ) continue;

//////////			Double_t MomV01Self[3], PosV01Self[3], MomV01PosSelf[3], MomV01NegSelf[3];
		Double_t xV02toV01  = PosV02[0] - PosV01Self[0];
		Double_t yV02toV01  = PosV02[1] - PosV01Self[1];
		Double_t zV02toV01  = PosV02[2] - PosV01Self[2];
		Double_t xyV02toV01 = TMath::Sqrt(xV02toV01*xV02toV01+yV02toV01*yV02toV01);
		Double_t rPVtoV0  = TMath::Sqrt(xV02toV01*xV02toV01+yV02toV01*yV02toV01+zV02toV01*zV02toV01);
		Double_t MomV02Pt  = GetPtFromPxPyPz(MomV02);
		Double_t cpa2 = (xV02toV01*MomV02[0]+yV02toV01*MomV02[1])/xyV02toV01/MomV02Pt;

		Double_t Mom2PosPt  = GetPtFromPxPyPz(MomV02Pos);
		Double_t Mom2NegPt  = GetPtFromPxPyPz(MomV02Neg);
		Double_t coa2PosNeg = (MomV02Pos[0]*MomV02Neg[0]+MomV02Pos[1]*MomV02Neg[1])/Mom2PosPt/Mom2NegPt;

		Double_t pDca = ptrk->GetD(PosV01Self[0],PosV01Self[1],fBzkG);
		Double_t nDca = ntrk->GetD(PosV01Self[0],PosV01Self[1],fBzkG);
		if ( pDca<2. ) continue;
		if ( nDca<2. ) continue;

		Double_t InvMassLambdaV02 = -9.;
		Double_t L2ReturnV0[1];
		L2ReturnV0[0] = 0.;
		InvMassLambdaV02     = InvMassLambda(    MomV02Pos,MomV02Neg,ptrk,ntrk,L2ReturnV0);

		if ( TMath::Abs(InvMassLambdaV02-mLambdaPDG)>fWindowLambda ) continue;

		TLorentzVector vV01,vV02,vH;
		vV01.SetXYZM( MomV01Self[0],MomV01Self[1],MomV01Self[2],InvMassLambdaV01Self);
		vV02.SetXYZM( MomV02[0],MomV02[1],MomV02[2],InvMassLambdaV02);
		vH = vV01  + vV02;
		Double_t cpaH = (vH.Px()*(PosV01Self[0]-PosPV[0])+vH.Py()*(PosV01Self[1]-PosPV[1]))/vH.Pt()/xyV01Self;
//		printf("!!!!! cpaH=%10f\n",cpaH);
		if ( cpaH<fCPADibaryon ) continue;

	//------------------------------------------------------------------------------------------
	// Output
	//------------------------------------------------------------------------------------------

				fCandidateVariables2[ 0] = antiLambdaOri1;
//				fCandidateVariables2[ 1] = antiLambdaStarOri1;
				fCandidateVariables2[ 1] = MomPosOri1Pt;
				fCandidateVariables2[ 2] = MomNegOri1Pt;
				fCandidateVariables2[ 3] = pos1DCAPV;
				fCandidateVariables2[ 4] = neg1DCAPV;
//				fCandidateVariables2[ 5] = InvMassLambdaV01; 
//				fCandidateVariables2[ 7] = InvMassLambdaStarV01;
//				fCandidateVariables2[ 8] = dcaV01;
//				fCandidateVariables2[ 9] = cpaV01;
//				fCandidateVariables2[10] = coaV01;
//				fCandidateVariables2[11] = dcaZV01;
//				fCandidateVariables2[12] = xyV01;
//				fCandidateVariables2[13] = rV01;
//				fCandidateVariables2[14] = ctauV01;
				fCandidateVariables2[ 5] = InvMassLambdaV01Self;
//				fCandidateVariables2[16] = InvMassLambdaStarV01Self;
				fCandidateVariables2[ 6] = dcaV01Self;
				fCandidateVariables2[ 7] = cpaV01Self;
				fCandidateVariables2[ 8] = coaV01Self;
				fCandidateVariables2[ 9] = dcaZV01Self;
				fCandidateVariables2[10] = xyV01Self;
				fCandidateVariables2[11] = rV01Self;
				fCandidateVariables2[12] = ctauV01Self;
//				fCandidateVariables2[24] = InvMassLambdaV01Cross;
//				fCandidateVariables2[25] = InvMassLambdaStarV01Cross;
//				fCandidateVariables2[26] = dcaV01Cross;
//				fCandidateVariables2[27] = cpaV01Cross;
//				fCandidateVariables2[28] = coaV01Cross;
//				fCandidateVariables2[29] = dcaZV01Cross;
//				fCandidateVariables2[30] = xyV01Cross;
//				fCandidateVariables2[31] = rV01Cross;
//				fCandidateVariables2[32] = ctauV01Cross;
				fCandidateVariables2[13] = fCentrality;
				fCandidateVariables2[14] = cpa2;
				fCandidateVariables2[15] = InvMassLambdaV02;
				fCandidateVariables2[16] = vH.M();
				fCandidateVariables2[17] = cpaH;
				fCandidateVariables2[18] = static_cast<Float_t>(IsOnFly);
				fCandidateVariables2[19] = pDca;
				fCandidateVariables2[20] = nDca;
				fCandidateVariables2[21] = coa2PosNeg;
				fCandidateVariables2[22] = dca2;
				fCandidateVariables2[23] = xyV02toV01;
				fCandidateVariables2[24] = rPVtoV0;

				fVariablesTree2->Fill();

			}//v0 loop

		}//track loop 2

	}//track loop 1

//	printf("### Number of V0s (END)  : %d\n",nV01Test);


	Double_t ratioLambda = static_cast<Double_t>(nAntiLambda)/static_cast<Double_t>(nLambda);
	Double_t ratioLambdaAll = static_cast<Double_t>(fCountAntiLambda)/static_cast<Double_t>(fCountLambda);
	printf("#######################################################################################################################\n");
	printf("#####               Lambda: Pro:%7d->%7d,  Anti:%7d->%7d  (ratio:%6.4f->%6.4f) ######################\n",nLambda,fCountLambda,nAntiLambda,fCountAntiLambda,ratioLambda,ratioLambdaAll);
	printf("#######################################################################################################################\n");
	printf("##### Analyzed events: %7d, Run Number: %10d, Centrality: %5d, nV0s: %5d ################################\n",fCountEvent,runNumber,fCentrality,nV0s);
	printf("#######################################################################################################################\n\n\n");

	fCountEvent++;

}
//______________________________________________________________________________________________________
void AliAnalysisTaskNOmegaLPK::DefineTreeVariables() {
  //
  // This is to define tree variables
  //

  const char* nameoutput1 = GetOutputSlot(1)->GetContainer()->GetName();
  fVariablesTree1 = new TTree(nameoutput1,"Candidates variables tree");
  Int_t nVar1 = 11;
  fCandidateVariables1 = new Float_t [nVar1];
  TString * fCandidateVariableNames1 = new TString[nVar1];

	fCandidateVariableNames1[ 0]="antiLambda";//= antiLambda;
	fCandidateVariableNames1[ 1]="antiLambdaStar";//= antiLambdaStar;
	fCandidateVariableNames1[ 2]="IMLV0";//= InvMassLambdaV0;
	fCandidateVariableNames1[ 3]="IMLSV0";//= InvMassLambdaStarV0;
	fCandidateVariableNames1[ 4]="DcaPos";//= testDcaPos;
	fCandidateVariableNames1[ 5]="DcaNeg";//= testDcaNeg;
	fCandidateVariableNames1[ 6]="DcaV0";//= testDcaV0;
	fCandidateVariableNames1[ 7]="FV";//= testFV;
	fCandidateVariableNames1[ 8]="R";//= testR;
	fCandidateVariableNames1[ 9]="CpaV0";//= testCpa;
	fCandidateVariableNames1[10]="CoaV0";//= CoaV0Daughters2D;

  for (Int_t ivar=0; ivar<nVar1; ivar++) {
    fVariablesTree1->Branch(fCandidateVariableNames1[ivar].Data(),&fCandidateVariables1[ivar],Form("%s/f",fCandidateVariableNames1[ivar].Data()));
  }


  const char* nameoutput2 = GetOutputSlot(2)->GetContainer()->GetName();
  fVariablesTree2 = new TTree(nameoutput2,"Candidates variables tree");
  Int_t nVar2 = 25;
  fCandidateVariables2 = new Float_t [nVar2];
  TString * fCandidateVariableNames2 = new TString[nVar2];

	fCandidateVariableNames2[ 0]="antiL1";//= antiLambdaOri1;
//	fCandidateVariableNames2[ 1]="antiLS1";//= antiLambdaStarOri1;
	fCandidateVariableNames2[ 1]="PtPos1";//= MomPosOri1Pt;
	fCandidateVariableNames2[ 2]="PtNeg1";//= MomNegOri1Pt;
	fCandidateVariableNames2[ 3]="DcaPos1";//= posDCAPV;
	fCandidateVariableNames2[ 4]="DcaNeg1";//= negDCAPV;
//	fCandidateVariableNames2[ 5]="IMLV01";//= InvMassLambdaV01; 
//	fCandidateVariableNames2[ 7]="IMLSV01";//= InvMassLambdaStarV01;
//	fCandidateVariableNames2[ 8]="Dca1";//= dcaV01;
//	fCandidateVariableNames2[ 9]="Cpa1";//= cpaV01;
//	fCandidateVariableNames2[10]="Coa1";//= coaV01;
//	fCandidateVariableNames2[11]="DcaZ1";//= dcaZV01;
//	fCandidateVariableNames2[12]="xy1";//= xyV01;
//	fCandidateVariableNames2[13]="r1";//= rV01;
//	fCandidateVariableNames2[14]="ctau1";//= ctauV01;
	fCandidateVariableNames2[ 5]="IMLV01Self";//= InvMassLambdaV01Self;
//	fCandidateVariableNames2[16]="IMLSV01Self";//= InvMassLambdaStarV01Self;
	fCandidateVariableNames2[ 6]="Dca1Self";//= dcaV01Self;
	fCandidateVariableNames2[ 7]="Cpa1Self";//= cpaV01Self;
	fCandidateVariableNames2[ 8]="Coa1Self";//= coaV01Self;
	fCandidateVariableNames2[ 9]="DcaZ1Self";//= dcaZV01Self;
	fCandidateVariableNames2[10]="xy1Self";//= xyV01Self;
	fCandidateVariableNames2[11]="r1Self";//= rV01Self;
	fCandidateVariableNames2[12]="ctau1Self";//= ctauV01Self;
//	fCandidateVariableNames2[24]="IMLV01Cross";//= InvMassLambdaV01Cross;
//	fCandidateVariableNames2[25]="IMLSV01Cross";//= InvMassLambdaStarV01Cross;
//	fCandidateVariableNames2[26]="Dca1Cross";//= dcaV01Cross;
//	fCandidateVariableNames2[27]="Cpa1Cross";//= cpaV01Cross;
//	fCandidateVariableNames2[28]="Coa1Cross";//= coaV01Cross;
//	fCandidateVariableNames2[29]="Dca1Cross";//= dcaZV01Cross;
//	fCandidateVariableNames2[30]="xy1Cross";//= xyV01Cross;
//	fCandidateVariableNames2[31]="r1Cross";//= rV01Cross;
//	fCandidateVariableNames2[32]="ctau1Cross";//= ctauV01Cross;
	fCandidateVariableNames2[13]="Centrality";//= fCentrality;
	fCandidateVariableNames2[14]="Cpa2";// = cpa2;
	fCandidateVariableNames2[15]="IMLV02";// = InvMassLambdaV02;
	fCandidateVariableNames2[16]="IMH";// = vH.M();
	fCandidateVariableNames2[17]="CpaH";// = cpaH;
	fCandidateVariableNames2[18]="OnFly2";// = static_cast<Float_t>(IsOnFly);
	fCandidateVariableNames2[19]="DcaPos2";// = pDca;
	fCandidateVariableNames2[20]="DcaNeg2";// = nDca;
	fCandidateVariableNames2[21]="Coa2";// = coa2PosNeg;
	fCandidateVariableNames2[22]="Dca2";// = dca2;
	fCandidateVariableNames2[23]="xy2";// = xyV02toV01;
	fCandidateVariableNames2[24]="r2";// = rPVtoV0;

  for (Int_t ivar=0; ivar<nVar2; ivar++) {
    fVariablesTree2->Branch(fCandidateVariableNames2[ivar].Data(),&fCandidateVariables2[ivar],Form("%s/f",fCandidateVariableNames2[ivar].Data()));
  }

  return;

}
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskNOmegaLPK::PreTrackCut(AliESDtrack *track) {

	//------------------------------------------------------------------------------------------
	// version 0.11 (2015/10/23)
	// First track cut for speed up (reduce the number of roops)
	// Input:  AliESDtrack: track
	// Return: kTRUE: Survived
	//        kFALSE: Excluded
	//------------------------------------------------------------------------------------------

	Double_t trackCharge = track->GetSign();

	if ( track->GetTPCNcls()<fClustersTPC ) return kFALSE;

	if ( !track->IsOn(AliESDtrack::kTPCrefit) ) return kFALSE;

	if ( track->GetKinkIndex(0)>0 ) return kFALSE;

	Double_t trackTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
	Double_t trackTPCPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion  ));
	Double_t trackTPCKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon  ));
//	if ( trackTPCProton>fSigmaTPC && trackTPCPion>fSigmaTPC && trackTPCKaon>fSigmaTPC ) return kFALSE;
	if ( trackTPCProton>fSigmaTPC && trackTPCPion>fSigmaTPC ) return kFALSE;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//	if (TMath::Abs(pos1DCAPV)>0.0182+0.035/MomPos01Pt) continue;
//	if (TMath::Abs(neg1DCAPV)>0.0182+0.035/MomNeg01Pt) continue;

	Double_t Mom[3];
	track->GetPxPyPz(Mom);
	TVector3 vMom;
	vMom.SetXYZ(Mom[0],Mom[1],Mom[2]);
	if (TMath::Abs(vMom.PseudoRapidity())>fPseudoRap) return kFALSE;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	return kTRUE;

}
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskNOmegaLPK::GetEsdV0Momentum(AliESDEvent *esdEvent, Int_t id1, Int_t id2, Double_t cutV0Parameters[9], Double_t ReturnV0[6], Double_t MomV0[3], Double_t PosV0[3], Double_t MomV0Pos[3], Double_t MomV0Neg[3]) {

	//------------------------------------------------------------------------------------------
	// version 0.30 (2016/02/01)
	// Check if the set of 2 tracks are assigned as ESD V0 and return ESD V0 momentum and position
	// Input:  Track IDs of 2 tracks, cut parameters
	// Return: kTRUE: Found
	//        kFALSE: Not found
	//      ReturnV0: 0:DCAxy, 1:CPA, 2:COA, 3:DCAz, 4:V0xy, 5:V0r
	//         V0Mom: Momentum of V0
	//         V0Pos: Position of V0
	//------------------------------------------------------------------------------------------

	Bool_t recut = kTRUE;

	for ( Int_t i=0; i<6; i++ ) {
		ReturnV0[i] = -999.;
	}
 
	const AliESDVertex *esdPV = (AliESDVertex*)esdEvent->GetPrimaryVertex();
  Double_t PosPV[3];
  esdPV->GetXYZ(PosPV);

	for ( Int_t i=0; i<3; i++ ) {
		MomV0[i] = -999.;
		PosV0[i] = -999.;
	}

	if ( id1==id2 ) return kFALSE;

	const Int_t nV0s = esdEvent->GetNumberOfV0s();
	if (nV0s==0) {
		return kFALSE;
	}

	for ( Int_t iV0 = 0; iV0<nV0s; iV0++ ) {
		AliESDv0 *esdV0 = esdEvent->GetV0(iV0);

		AliESDtrack *ptrk = esdEvent->GetTrack(esdV0->GetPindex());
		AliESDtrack *ntrk = esdEvent->GetTrack(esdV0->GetNindex());

		Double_t ptrkCharge = ptrk->GetSign();
		Double_t ntrkCharge = ntrk->GetSign();
		if        ( ptrkCharge>0 && ntrkCharge>0 ) {//(3,4)=(+,+)->continue
			return kFALSE;
		} else if ( ptrkCharge<0 && ntrkCharge<0 ) {//(3,4)=(-,-)->continue
			return kFALSE;
		} else if ( ptrkCharge>0 && ntrkCharge<0 ) {//(3,4)=(+,-)->ok
		} else if ( ptrkCharge<0 && ntrkCharge>0 ) {//(3,4)=(-,+)->(+,-)
			ptrk = esdEvent->GetTrack(esdV0->GetNindex());
			ntrk = esdEvent->GetTrack(esdV0->GetPindex());
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
			Double_t cpa = (xPVtoV0*MomV0[0]+yPVtoV0*MomV0[1])/xyPVtoV0/MomV0Pt;

			Double_t MomPosPt  = GetPtFromPxPyPz(MomV0Pos);
			Double_t MomNegPt  = GetPtFromPxPyPz(MomV0Neg);
			Double_t coaPosNeg = (MomV0Pos[0]*MomV0Neg[0]+MomV0Pos[1]*MomV0Neg[1])/MomPosPt/MomNegPt;

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
/*
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
*/
				Int_t debugESD = 0;
				if ( TMath::Abs(pDca)<fPosDCAToPVMin1 ) debugESD = 1;
				if ( TMath::Abs(pDca)>fFiducialVolMax1 ) debugESD = 2;
				if ( TMath::Abs(nDca)<fNegDCAToPVMin1 ) debugESD = 3;
				if ( TMath::Abs(nDca)>fFiducialVolMax1 ) debugESD = 4;

				if ( dcaCut>fDCADaughterMax1 ) debugESD = 5;
				if ( (xPos+xNeg)>2*fFiducialVolMax1 ) debugESD = 6;
				if ( (xPos+xNeg)<2*fFiducialVolMin1 ) debugESD = 7;

				if ( coaPosNeg<fCOADaughterMin1 ) debugESD = 8;

				if ( rV0<fFiducialVolMin1) debugESD = 9;
				if ( rV0>fFiducialVolMax1) debugESD = 10;

				if ( cpa<fCPAMin1 ) debugESD = 11;

				if ( dcaV0ToPV<fDCAToPVMin1 ) debugESD = 12;

				if ( TMath::Abs(dcaZ)>fDCAZDaughterMax1 ) debugESD = 13;
//				printf("####### Cut at #%d\n",debugESD);
				if ( debugESD > 0 ) return kFALSE;

			}

			return kTRUE;

		}

	}

	return kFALSE;

}
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskNOmegaLPK::GetSelfV0Momentum(AliESDtrack *posTrk, AliESDtrack *negTrk, Double_t cutV0Parameters[9],Double_t ReturnV0[6], Double_t MomV0[3], Double_t PosV0[3], Double_t MomPos[3], Double_t MomNeg[3]) {

	//------------------------------------------------------------------------------------------
	// version 0.30 (2016/02/01)
	// Copy from AliV0vertexer
	// reconstruct V0 from 2 tracks
	// Input:  2 tracks of V0 positive and negative tracks, cut parameters
	// Return: kTRUE: Success
	//        kFALSE: Failed
	//      ReturnV0: 0:DCAxy, 1:CPA, 2:COA, 3:DCAz, 4:V0xy, 5:V0r
	//         V0Mom: Momentum of V0
	//         V0Pos: Position of V0
	//------------------------------------------------------------------------------------------

	for ( Int_t i=0; i<6; i++ ) {
		ReturnV0[i] = -999.;
	}
 
	AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);

  const AliESDVertex *esdPV = (AliESDVertex*)esdEvent->GetPrimaryVertex();
  Double_t PosPV[3];
  esdPV->GetXYZ(PosPV);

	Int_t eventNumber = esdEvent->GetEventNumberInFile();

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
	if (TMath::Abs(posDCAPV)<fPosDCAToPVMin1) return kFALSE;
	if (TMath::Abs(posDCAPV)>fFiducialVolMax1) return kFALSE;

	// negative track cut
	Double_t negDCAPV = negTrk->GetD(PosPV[0],PosPV[1],fBzkG);
	if (TMath::Abs(negDCAPV)<fNegDCAToPVMin1) return kFALSE;
	if (TMath::Abs(negDCAPV)>fFiducialVolMax1) return kFALSE;
//	printf("#### Self pDca: %10f(%10f-%10f), nDca: %10f(%10f-%10f)\n",posTrk->GetD(PosPV[0],PosPV[1],fBzkG),fPosDCAToPVMin1,fFiducialVolMax1,negTrk->GetD(PosPV[0],PosPV[1],fBzkG),fNegDCAToPVMin1,fFiducialVolMax1);


	//------------------------------------------------------------------------------------------
	// Reconstruct V01 and cut V01
	//------------------------------------------------------------------------------------------

	// get DCA, propagate tracks, set V0
	Double_t xPos;
	Double_t xNeg;
	Double_t dcaV0Daughters;
	dcaV0Daughters = posTrk->GetDCA(negTrk,fBzkG,xPos,xNeg);
	ReturnV0[0] = dcaV0Daughters;
	if (dcaV0Daughters > fDCADaughterMax1) return kFALSE;
	if ((xPos+xNeg) > 2*fFiducialVolMax1) return kFALSE;
	if ((xPos+xNeg) < 2*fFiducialVolMin1) return kFALSE;

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
		if (dcaV0Daughters > fDCADaughterMax1) return kFALSE;
		if ((xPos+xNeg) > 2*fFiducialVolMax1) return kFALSE;
		if ((xPos+xNeg) < 2*fFiducialVolMin1) return kFALSE;
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
	if ( TMath::Abs(dcaZ)>fDCAZDaughterMax1 ) return kFALSE;

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
	ReturnV0[2] = coaV0Daughters2D;
	if (coaV0Daughters2D < fCOADaughterMin1) return kFALSE;

	Double_t r2V0 = PosV0[0]*PosV0[0] + PosV0[1]*PosV0[1];
	Double_t rV0  = TMath::Sqrt(r2V0);
	if (rV0 < fFiducialVolMin1) return kFALSE;
	if (rV0 > fFiducialVolMax1) return kFALSE;

	Double_t xPVtoV0  = PosV0[0] - PosPV[0];
	Double_t yPVtoV0  = PosV0[1] - PosPV[1];
	Double_t zPVtoV0  = PosV0[2] - PosPV[2];
	Double_t cpaPVtoV0 = (xPVtoV0*MomV0[0]+yPVtoV0*MomV0[1]+zPVtoV0*MomV0[2])/MomV0P
	                      /TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0+zPVtoV0*zPVtoV0);
	ReturnV0[1] = cpaPVtoV0;
	if (cpaPVtoV0 < fCPAMin1) return kFALSE;

	Double_t xyPVtoV0 = TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0);
	ReturnV0[4] = xyPVtoV0;
	Double_t rPVtoV0  = TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0+zPVtoV0*zPVtoV0);
	ReturnV0[5] = rPVtoV0;

	Double_t dcaV0ToPV = TMath::Abs(MomV0[0]/MomV0Pt*yPVtoV0-MomV0[1]/MomV0Pt*xPVtoV0);
  if (dcaV0ToPV < fDCAToPVMin1) return kFALSE;
//	printf("1,2:posDCA    , 3,4:posDCA    , 5:dcaV0Daugh, 6,7:xPos+xNeg , 8:coaV0Daugh, 9,10:rV0       , 11:cpaPVtoV0 , 12:dcaV0ToPV , 13:dcaZ     \n",TMath::Abs(posDCAPV),TMath::Abs(negDCAPV),dcaV0Daughters,(xPos+xNeg),coaV0Daughters2D,rV0,cpaPVtoV0,dcaV0ToPV,TMath::Abs(dcaZ));
//	printf("1,2:%10f, 3,4:%10f, 5:%10f, 6,7:%10f, 8:%10f, 9,10:%10f, 11:%10f, 12:%10f, 13:%10f\n",TMath::Abs(posDCAPV),TMath::Abs(negDCAPV),dcaV0Daughters,(xPos+xNeg),coaV0Daughters2D,rV0,cpaPVtoV0,dcaV0ToPV,TMath::Abs(dcaZ));

	return kTRUE; 


}
//______________________________________________________________________________________________________
Bool_t AliAnalysisTaskNOmegaLPK::GetSelfV0CrossMomentum(AliESDtrack *posTrk, AliESDtrack *negTrk,Double_t cutV0Parameters[9], Double_t ReturnV0[6], Double_t MomV0[3], Double_t PosV0[3], Double_t MomPos[3], Double_t MomNeg[3]) {

	//------------------------------------------------------------------------------------------
	// version 1.30 (2016/02/01) !!! There is a BUG in getting zVtx for ver0.xx !!!
	// reconstruct V0 without ExternalTrackParam class (by hand)
	// Input:  2 tracks of V0 positive and negative tracks, cut parameters
	// Return: kTRUE: Success
	//        kFALSE: Failed
	//      ReturnV0: 0:DCAxy, 1:CPA, 2:COA, 3:DCAz, 4:V0xy, 5:V0r
	//         V0Mom: Momentum of V0
	//         V0Pos: Position of V0
	//------------------------------------------------------------------------------------------

//	printf("\n#########################################################\n");
	// Select cross mode ON/OFF (kTRUE:all cross V0 accepted, kFALSE:cross V0 with DCA cut accepted)
	Bool_t modeCross = kTRUE;

	Int_t debugCross = 0;

	for ( Int_t i=0; i<6; i++ ) {
		ReturnV0[i] = -999.;
	}
 
	AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);

  const AliESDVertex *esdPV = (AliESDVertex*)esdEvent->GetPrimaryVertex();
  Double_t PosPV[3];
  esdPV->GetXYZ(PosPV);

	// Charge of track (-1 or +1)
	Double_t posCharge = posTrk->Charge();
	Double_t negCharge = negTrk->Charge();
	if ( posCharge==-1. || negCharge==1. ) {
		printf("#####\n##### ERROR Get V0 because posCharge<0 or negCharge>0 !!! \n#####\n");
//		return kFALSE;
		debugCross = 1;
	}

	// Position(x,y,z) and Momentum(Px,Py,Pz) in global coordinates
	Double_t XYZTrk1[3];
	Double_t MomTrk1[3];
	Double_t XYZTrk2[3];
	Double_t MomTrk2[3];
	Double_t XYZVtx[3];
	Double_t MomVtx[3];

	posTrk->GetXYZ(XYZTrk1);
	posTrk->GetPxPyPz(MomTrk1);// This function returns the global track momentum components
	negTrk->GetXYZ(XYZTrk2);
	negTrk->GetPxPyPz(MomTrk2);// This function returns the global track momentum components

	const Double_t xTrk1Ori = XYZTrk1[0];
	const Double_t yTrk1Ori = XYZTrk1[1];
	const Double_t zTrk1Ori = XYZTrk1[2];
	TVector2 vMomTrk1;
	vMomTrk1.Set(MomTrk1[0],MomTrk1[1]);
	const Double_t modMomTrk1 = vMomTrk1.Mod();
	vMomTrk1 = vMomTrk1.Unit();

	const Double_t xTrk2Ori = XYZTrk2[0];
	const Double_t yTrk2Ori = XYZTrk2[1];
	const Double_t zTrk2Ori = XYZTrk2[2];
	TVector2 vMomTrk2;
	vMomTrk2.Set(MomTrk2[0],MomTrk2[1]);
	const Double_t modMomTrk2 = vMomTrk2.Mod();
	vMomTrk2 = vMomTrk2.Unit();

	// calculate tangent of dip angle of track
	Double_t tanDipAngleMom1 = MomTrk1[2]/modMomTrk1;
	Double_t tanDipAngleMom2 = MomTrk2[2]/modMomTrk2;

	// calculate curvature of track
	const Double_t crv1  = fBzkG*kB2C/modMomTrk1;//kB2C=0.299792458e-3
	const Double_t crv2  = fBzkG*kB2C/modMomTrk2;//kB2C=0.299792458e-3

	// calculate radius of track
	const Double_t radius1 = 1./crv1;
	const Double_t radius2 = 1./crv2;

	// calculate center of helix
	TVector2 vCircleCenter1, vCircleCenter2;
	vCircleCenter1.Set(xTrk1Ori-vMomTrk1.Y()*radius1* 1.,yTrk1Ori+vMomTrk1.X()*radius1* 1.);// Pos track
	vCircleCenter2.Set(xTrk2Ori-vMomTrk2.Y()*radius2*-1.,yTrk2Ori+vMomTrk2.X()*radius2*-1.);// Neg track
	Double_t distanceTest = TMath::Sqrt(TMath::Power(vCircleCenter1.X()-vCircleCenter2.X(),2)+TMath::Power(vCircleCenter1.Y()-vCircleCenter2.Y(),2));//

	// calculate the angle from X-axis to center1->center2 in global coordinates
	TVector2 vCent1toCent2;
	vCent1toCent2.Set(vCircleCenter2.X()-vCircleCenter1.X(),vCircleCenter2.Y()-vCircleCenter1.Y());

//	if (vCent1toCent2.Mod()==0.) return kFALSE;
	if (vCent1toCent2.Mod()==0.) debugCross = 2;
	Double_t cosCent1toCent2 = vCent1toCent2.X()/vCent1toCent2.Mod();
	Double_t sinCent1toCent2 = vCent1toCent2.Y()/vCent1toCent2.Mod();

	const Double_t phiCent1toCent2 = TVector2::Phi_mpi_pi(vCent1toCent2.Phi());
	vCent1toCent2 = vCent1toCent2.Rotate(-phiCent1toCent2);

	// rotate center in original global coordinates (direction of center1->center2 is rotated to parallel to X-axis)
	vCircleCenter2.Set(vCircleCenter2.X()-vCircleCenter1.X(),vCircleCenter2.Y()-vCircleCenter1.Y());
	vCircleCenter2 = vCircleCenter2.Rotate(-phiCent1toCent2);
	Double_t dCenters = vCircleCenter2.X();
//	printf("! rotated (xCen2,yCen2) = (%11f,%11f)) (=(x,0), =%11f(distance))\n",vCircleCenter2.X(),vCircleCenter2.Y(),distanceTest);

	//////////////////////////
	// select vertex candidate
	Double_t SignedDca = 0.;
	Double_t dcaXY     = 0.;
	Double_t xVtx      = 0.;
	Double_t yVtx      = 0.;
	Double_t xVtxDCA   = 0.;
	Double_t yVtxDCA   = 0.;
	Bool_t crossTrack  = kFALSE;

	Double_t yVtxP = 0.;
	Double_t yVtxN = 0.;
	SignedDca =  dCenters - (radius1+radius2);
	//NOTE!! SignedDca<0 means crossing (or one including another)! abs(signedDCA) is not real DCA !!
	dcaXY = TMath::Abs(SignedDca);
	xVtx  = radius1 + SignedDca/2.;
	xVtxDCA = xVtx;
	if ( modeCross && SignedDca<0 ) {//crossing tracks with crossMode=ON
		crossTrack = kTRUE;
		dcaXY =  0.;
		xVtx = (radius1*radius1 - radius2*radius2 + dCenters*dCenters)/(2.*dCenters);// x of crossing point
		yVtxP =  TMath::Sqrt(radius1*radius1-xVtx*xVtx);
		yVtxN = -TMath::Sqrt(radius1*radius1-xVtx*xVtx);
	}
	ReturnV0[0] = dcaXY;
	if ( dcaXY>fDCADaughterMax1 ) return kFALSE;
	printf("#### SignedDca = %10f\n",SignedDca);
//	if ( dcaXY>fDCADaughterMax1 ) debugCross = 3;
//	printf("!!! xVtx: %11f=%11f\n",radius1+SignedDca/2.,(radius1*radius1-radius2*radius2+dCenters*dCenters)/(2.*dCenters));
//	printf("!!! xVtx via Centers:%11f=%11f, Dca/2. via Cos:%11f=%11f\n",radius1+SignedDca/2.,dCenters-radius2-SignedDca/2.,radius1-(radius1*radius1-radius2*radius2+dCenters*dCenters)/(2.*dCenters),dCenters-radius2-(radius1*radius1-radius2*radius2+dCenters*dCenters)/(2.*dCenters));
//	printf("!!! Sin via Centers:%11f=%11f, via Cos:%11f=%11f\n",TMath::Sqrt(1.-TMath::Power((radius1+SignedDca/2.)/radius1,2))*radius1,TMath::Sqrt(1.-TMath::Power(((radius1+SignedDca/2.-dCenters)/radius2),2))*radius2,TMath::Sqrt(1.-TMath::Power((radius1*radius1-radius2*radius2+dCenters*dCenters)/(2.*radius1*dCenters),2))*radius1,TMath::Sqrt(1.-TMath::Power(((radius1*radius1-radius2*radius2+dCenters*dCenters)/(2.*dCenters)-dCenters)/radius2,2))*radius2);
//	printf("!!! radius1=%11f, radius2=%11f, dCenters=%11f\n",radius1,radius2,dCenters);

	// rotate others in original global coordinates (direction of center1->center2 is rotated to parallel to X-axis)
	TVector2 vTrk1, vTrk2;
	vTrk1.Set(xTrk1Ori-vCircleCenter1.X(),yTrk1Ori-vCircleCenter1.Y());
	vTrk2.Set(xTrk2Ori-vCircleCenter1.X(),yTrk2Ori-vCircleCenter1.Y());
	vTrk1 = vTrk1.Rotate(-phiCent1toCent2);
	vTrk2 = vTrk2.Rotate(-phiCent1toCent2);
///	printf("# rotated (xTrk1,yTrk1) = (%11f,%11f))\n",vTrk1.X(),vTrk1.Y());
///	printf("# rotated (xTrk2,yTrk2) = (%11f,%11f))\n",vTrk2.X(),vTrk2.Y());

	Double_t cosTrk1Cross, cosTrk2Cross;
	TVector2 vTrk1CrossP, vTrk2CrossP, vTrk1CrossN, vTrk2CrossN;
	if ( modeCross && crossTrack ) {
		cosTrk1Cross = xVtx/radius1;
		cosTrk2Cross = (xVtx-dCenters)/radius2;

//		if ( TMath::Abs(cosTrk1Cross)>1. || TMath::Abs(cosTrk2Cross)>1. ) return kFALSE;//one circle including another
		if ( TMath::Abs(cosTrk1Cross)>1. || TMath::Abs(cosTrk2Cross)>1. ) debugCross = 4;//one circle including another

		vTrk1CrossP.Set(cosTrk1Cross, TMath::Sqrt(1.-cosTrk1Cross*cosTrk1Cross));
		vTrk1CrossN.Set(cosTrk1Cross,-TMath::Sqrt(1.-cosTrk1Cross*cosTrk1Cross));
		vTrk2CrossP.Set(cosTrk2Cross, TMath::Sqrt(1.-cosTrk2Cross*cosTrk2Cross));
		vTrk2CrossN.Set(cosTrk2Cross,-TMath::Sqrt(1.-cosTrk2Cross*cosTrk2Cross));
	}

	const Double_t phiTrk1CrossP = TVector2::Phi_mpi_pi(vTrk1CrossP.Phi());
	const Double_t phiTrk1CrossN = TVector2::Phi_mpi_pi(vTrk1CrossN.Phi());
	const Double_t phiTrk2CrossP = TVector2::Phi_mpi_pi(vTrk2CrossP.Phi());
	const Double_t phiTrk2CrossN = TVector2::Phi_mpi_pi(vTrk2CrossN.Phi());

	vMomTrk1 = vMomTrk1.Rotate(-phiCent1toCent2);
	TVector2 vTrk1ToRotate, vTrk1ToRotateP, vTrk1ToRotateN;
	vTrk1ToRotate = vMomTrk1;
	vTrk1ToRotate = vTrk1ToRotate.Rotate(-TMath::Pi()/2.);
	if ( modeCross && crossTrack ) {
		vTrk1ToRotateP = vTrk1ToRotate.Rotate(-phiTrk1CrossP);
		vTrk1ToRotateN = vTrk1ToRotate.Rotate(-phiTrk1CrossN);
	}
	const Double_t phiTrk1ToRotate0 = TVector2::Phi_mpi_pi(-vTrk1ToRotate.Phi());
	const Double_t phiTrk1ToRotateP = TVector2::Phi_mpi_pi(-vTrk1ToRotateP.Phi());
	const Double_t phiTrk1ToRotateN = TVector2::Phi_mpi_pi(-vTrk1ToRotateN.Phi());

	vMomTrk2 = vMomTrk2.Rotate(-phiCent1toCent2);
	TVector2 vTrk2ToRotate, vTrk2ToRotateP, vTrk2ToRotateN;
	vTrk2ToRotate = vMomTrk2;
	vTrk2ToRotate = vTrk2ToRotate.Rotate(TMath::Pi()/2.);
	if ( modeCross && crossTrack ) {
		vTrk2ToRotateP = vTrk2ToRotate.Rotate(-phiTrk2CrossP);
		vTrk2ToRotateN = vTrk2ToRotate.Rotate(-phiTrk2CrossN);
	}
	vTrk2ToRotate = vTrk2ToRotate.Rotate(-TMath::Pi());
	const Double_t phiTrk2ToRotate0 = TVector2::Phi_mpi_pi(-vTrk2ToRotate.Phi());
	const Double_t phiTrk2ToRotateP = TVector2::Phi_mpi_pi(-vTrk2ToRotateP.Phi());
	const Double_t phiTrk2ToRotateN = TVector2::Phi_mpi_pi(-vTrk2ToRotateN.Phi());

	if ( modeCross && crossTrack ) {
//		printf("### Angle: Trk1: P:%11f, 0:%11f, N:%11f, Trk2: P:%11f, 0:%11f, N:%11f\n",TVector2::Phi_mpi_pi(phiTrk1CrossP),TVector2::Phi_mpi_pi((vMomTrk1.Rotate(-TMath::Pi()/2.).Phi())),TVector2::Phi_mpi_pi(phiTrk1CrossN),TVector2::Phi_mpi_pi(phiTrk2CrossP),TVector2::Phi_mpi_pi((vMomTrk2.Rotate(TMath::Pi()/2.).Phi())),TVector2::Phi_mpi_pi(phiTrk2CrossN));
//		printf("### Angle: Trk1: P:%11f, 0:%11f, N:%11f, Trk2: P:%11f, 0:%11f, N:%11f\n",TVector2::Phi_mpi_pi(phiTrk1ToRotateP-phiTrk1ToRotate0),TVector2::Phi_mpi_pi(-phiTrk1ToRotate0),TVector2::Phi_mpi_pi(phiTrk1ToRotateN-phiTrk1ToRotate0),TVector2::Phi_mpi_pi(TMath::Pi()+phiTrk2ToRotateP-phiTrk2ToRotate0),TVector2::Phi_mpi_pi(TMath::Pi()-phiTrk2ToRotate0),TVector2::Phi_mpi_pi(TMath::Pi()+phiTrk2ToRotateN-phiTrk2ToRotate0));
	}

	//-------------------------------------------------------------------------------------------------------------------
	// calculate zVtx
	//   for DCA
	Double_t dRtTrk1, dZTrk1, dRtTrk2, dZTrk2;
	Double_t zTrk1, zTrk2, dcaZ, zVtx, dca3D;
	Double_t zTrk1P, zTrk2P, dcaZP, zVtxP;
	Double_t zTrk1N, zTrk2N, dcaZN, zVtxN;
	dRtTrk1 = phiTrk1ToRotate0*radius1;
	dZTrk1  = dRtTrk1*tanDipAngleMom1*1.;// Pos track
	zTrk1   = zTrk1Ori + dZTrk1;
	dRtTrk2 = phiTrk2ToRotate0*radius2;
	dZTrk2  = dRtTrk2*tanDipAngleMom2*-1.;// Neg track
	zTrk2   = zTrk2Ori + dZTrk2;
	dcaZ    =   zTrk1 - zTrk2;
	zVtx    = ( zTrk1 + zTrk2 ) / 2.;
	dca3D   = TMath::Sqrt(SignedDca*SignedDca+dcaZ*dcaZ);

	//   for yVtxP
	if ( modeCross && crossTrack ) {
		dRtTrk1 = phiTrk1ToRotateP*radius1;
		dZTrk1  = dRtTrk1*tanDipAngleMom1*1.;// Pos track
		zTrk1P  = zTrk1Ori + dZTrk1;
		dRtTrk2 = phiTrk2ToRotateP*radius2;
		dZTrk2  = dRtTrk2*tanDipAngleMom2*-1.;// Neg track
		zTrk2P  = zTrk2Ori + dZTrk2;
		dcaZP   =   zTrk1P - zTrk2P;
		zVtxP   = ( zTrk1P + zTrk2P ) / 2.;
	}

	//   for yVtxN
		if ( modeCross && crossTrack ) {
		dRtTrk1 = phiTrk1ToRotateN*radius1;
		dZTrk1  = dRtTrk1*tanDipAngleMom1*1.;// Pos track
		zTrk1N  = zTrk1Ori + dZTrk1;
		dRtTrk2 = phiTrk2ToRotateN*radius2;
		dZTrk2  = dRtTrk2*tanDipAngleMom2*-1.;// Neg track
		zTrk2N  = zTrk2Ori + dZTrk2;
		dcaZN   =   zTrk1N - zTrk2N;
		zVtxN   = ( zTrk1N + zTrk2N ) / 2.;
	}

	if ( modeCross && crossTrack ) {
//		if ( (zTrk1P < zTrk1 && zTrk1 < zTrk1N) || (zTrk1P > zTrk1 && zTrk1 > zTrk1N) ) return kFALSE;//for debug
//		if ( (zTrk2P < zTrk2 && zTrk2 < zTrk1N) || (zTrk2P > zTrk2 && zTrk2 > zTrk2N) ) return kFALSE;//for debug
//		printf("\n### cos sin: Trk1: 1=%11f (%11f,%11f), Trk2: 1=%11f (%11f,%11f)\n",vTrk1CrossP.X()*vTrk1CrossP.X()+vTrk1CrossP.Y()*vTrk1CrossP.Y(),vTrk2CrossP.X(),vTrk2CrossP.Y(),vTrk2CrossP.X()*vTrk2CrossP.X()+vTrk2CrossP.Y()*vTrk2CrossP.Y(),vTrk2CrossP.X(),vTrk2CrossP.Y());
		printf("### Z: Trk1: P:%11f 0:%11f, N:%11f, Trk2: P:%11f 0:%11f, N:%11f\n",zTrk1P,zTrk1,zTrk1N,zTrk2P,zTrk2,zTrk2N);
		printf("### DCA  VtxP:%11f, Vtx0:%11f, VtxN:%11f  (DCAxy:%11f, sDCAxy:%11f, DCAz:%11f)\n",dcaZP,dca3D,dcaZN,dcaXY,SignedDca,dcaZ);
	}

	Int_t v0Position = 0;
	if ( modeCross && crossTrack ) {
		if ( TMath::Abs(dcaZP) < TMath::Abs(dcaZN) ) {//P or 0
			if ( TMath::Abs(dcaZP) < dca3D ) {//P
				v0Position = 1;
			} else {//0
				v0Position = 0;
			}
		} else {//N or 0
			if ( TMath::Abs(dcaZN) < dca3D ) {//N
				v0Position = -1;
			} else {//0
				v0Position = 0;
			}
		}
	} else {//no cross mode
		v0Position = 0;
	}

	Double_t phiTrk1ToRotate;
	Double_t phiTrk2ToRotate;
	if ( v0Position == 0 ) {
		yVtx = 0;
		dcaZ = dca3D;
		phiTrk1ToRotate = phiTrk1ToRotate0;
		phiTrk2ToRotate = phiTrk2ToRotate0;
///		printf("###### V0 position is determined from DCA\n");
	} else if ( v0Position == 1 ) {
		yVtx = yVtxP;
		zVtx = zVtxP;
		dcaZ = dcaZP;
		phiTrk1ToRotate = phiTrk1ToRotateP;
		phiTrk2ToRotate = phiTrk2ToRotateP;
//		printf("###### V0 position is determined from CrossP\n");
	} else if ( v0Position == -1 ) {
		yVtx = yVtxN;
		zVtx = zVtxN;
		dcaZ = dcaZN;
		phiTrk1ToRotate = phiTrk1ToRotateN;
		phiTrk2ToRotate = phiTrk2ToRotateN;
//		printf("###### V0 position is determined from CrossN\n");
//	} else return kFALSE;
	} else debugCross = 5;



//	if ( TMath::Abs(dcaZ)>fDCAZDaughterMax1 ) return kFALSE;
	if ( TMath::Abs(dcaZ)>fDCAZDaughterMax1 ) debugCross = 6;
//		printf("\n###### V0 position is determined from %d (1:P, 0:DCA, -1:N)\n",v0Position);
//		printf("### Z: Trk1: P:%11f 0:%11f, N:%11f, Trk2: P:%11f 0:%11f, N:%11f\n",zTrk1P,zTrk1,zTrk1N,zTrk2P,zTrk2,zTrk2N);
//		printf("### DCA  VtxP:%11f, Vtx0:%11f, VtxN:%11f  (DCAxy:%11f, sDCAxy:%11f, DCAz:%11f)\n",dcaZP,dca3D,dcaZN,dcaXY,SignedDca,dcaZ);

	//-------------------------------------------------------------------------------------------------------------------

	// propagate to DCA
	vTrk1 = vTrk1.Rotate(phiTrk1ToRotate);
	vTrk2.Set(vTrk2.X()-dCenters,vTrk2.Y());
	vTrk2 = vTrk2.Rotate(phiTrk2ToRotate);
	vTrk2.Set(vTrk2.X()+dCenters,vTrk2.Y());
//	printf("! Propagated pos position: (%11f,%11f) (=(x,0))\n",vTrk1.X(),vTrk1.Y());
//	printf("! Propagated neg position: (%11f,%11f) (=(x,0))\n",vTrk2.X(),vTrk2.Y());
//	printf("! (pos+neg)/2 position: (%11f,%11f) (=(x,0))\n",(vTrk1.X()+vTrk2.X())/2., (vTrk1.Y()+vTrk2.Y())/2.);

	vMomTrk1 = vMomTrk1.Rotate(phiTrk1ToRotate);
//	printf("! Propagated pos momentum: (%11f,%11f) (=(0,+1))\n",vMomTrk1.X(),vMomTrk1.Y());

	vMomTrk2 = vMomTrk2.Rotate(phiTrk2ToRotate);
//	printf("! Propagated neg momentum: (%11f,%11f) (=(0,+1))\n",vMomTrk2.X(),vMomTrk2.Y());

	// re-rotate in original global coordinates
	TVector2 vVtx;
	vVtx.Set(xVtx,yVtx);
	vTrk1 = vTrk1.Rotate(phiCent1toCent2);
	vTrk2 = vTrk2.Rotate(phiCent1toCent2);
	vVtx  = vVtx.Rotate(phiCent1toCent2);
	vTrk1.Set(vTrk1.X()+vCircleCenter1.X(),vTrk1.Y()+vCircleCenter1.Y());
	vTrk2.Set(vTrk2.X()+vCircleCenter1.X(),vTrk2.Y()+vCircleCenter1.Y());
	vVtx.Set(vVtx.X()+vCircleCenter1.X(),vVtx.Y()+vCircleCenter1.Y());
//	printf("# Propagated (xTrk1,yTrk1) = (%11f,%11f))\n",vTrk1.X(),vTrk1.Y());
//	printf("# Propagated (xTrk2,yTrk2) = (%11f,%11f))\n",vTrk2.X(),vTrk2.Y());
//	printf("# Propagated (xVtx, yVtx ) = (%11f,%11f))\n",vVtx.X(), vVtx.Y());

	vMomTrk1 = vMomTrk1.Rotate(phiCent1toCent2);
	vMomTrk2 = vMomTrk2.Rotate(phiCent1toCent2);
	vMomTrk1 = vMomTrk1*modMomTrk1;
	vMomTrk2 = vMomTrk2*modMomTrk2;

	Double_t r2Cross  = vVtx.X()*vVtx.X()+vVtx.Y()*vVtx.Y();
//	if ( TMath::Sqrt(r2Cross)<fFiducialVolMin1 || TMath::Sqrt(r2Cross)>fFiducialVolMax1 ) return kFALSE;
	if ( TMath::Sqrt(r2Cross)<fFiducialVolMin1 || TMath::Sqrt(r2Cross)>fFiducialVolMax1 ) debugCross = 7;

	// output
	PosV0[0]  = vVtx.X();
	PosV0[1]  = vVtx.Y();
	PosV0[2]  = zVtx;
	MomV0[0]  = vMomTrk1.X() + vMomTrk2.X();
	MomV0[1]  = vMomTrk1.Y() + vMomTrk2.Y();
	MomV0[2]  = MomTrk1[2] + MomTrk2[2];
	MomPos[0] = vMomTrk1.X();
	MomPos[1] = vMomTrk1.Y();
	MomPos[2] = MomTrk1[2];//No change
	MomNeg[0] = vMomTrk2.X();
	MomNeg[1] = vMomTrk2.Y();
	MomNeg[2] = MomTrk2[2];//No change

	Double_t xPVtoV0  = PosV0[0] - PosPV[0];
	Double_t yPVtoV0  = PosV0[1] - PosPV[1];
	Double_t zPVtoV0  = PosV0[2] - PosPV[2];
	Double_t xyPVtoV0 = TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0);
	ReturnV0[4] = xyPVtoV0;
	Double_t rPVtoV0  = TMath::Sqrt(xPVtoV0*xPVtoV0+yPVtoV0*yPVtoV0+zPVtoV0*zPVtoV0);
	ReturnV0[5] = rPVtoV0;
	Double_t MomV0Pt  = GetPtFromPxPyPz(MomV0);
	Double_t cpa = (xPVtoV0*MomV0[0]+yPVtoV0*MomV0[1])/xyPVtoV0/MomV0Pt;
	ReturnV0[1] = cpa;
//	if ( cpa<fCPAMin1 ) return kFALSE;
	if ( cpa<fCPAMin1 ) debugCross = 8;

	Double_t coaPos    = (vMomTrk1.X()*MomV0[0]+vMomTrk1.Y()*MomV0[1])/vMomTrk1.Mod()/MomV0Pt;
	Double_t coaNeg    = (vMomTrk2.X()*MomV0[0]+vMomTrk2.Y()*MomV0[1])/vMomTrk2.Mod()/MomV0Pt;
	Double_t coaPosNeg = (vMomTrk1.X()*vMomTrk2.X()+vMomTrk1.Y()*vMomTrk2.Y())/vMomTrk1.Mod()/vMomTrk2.Mod();
	ReturnV0[2] = coaPosNeg;
//	if ( coaPosNeg<fCOADaughterMin1 ) return kFALSE;
	if ( coaPosNeg<fCOADaughterMin1 ) debugCross = 9;
//	printf("########## Cut at #%d (Cross)\n",debugCross);



//	printf("######## dcaZ:%10f(%10f), r2:%10f(%10f-%10f), cpa:%10f(%10f)\n",TMath::Abs(dcaZ),fDCAZDaughterMax1,TMath::Sqrt(r2Cross),fFiducialVolMin1,fFiducialVolMax1,cpa,fCPAMin1);



	if ( debugCross > 0 ) return kFALSE;

	// Summary
/*
	printf("\n################# V0 found!! DCAxy=%9f, DCAz=%9f, charge: %2.0f(trk1),%2.0f(trk2) ##################\n",dcaXY,dcaZ,posCharge,negCharge);
	printf("### Trk1Ori XYZ:(%11f,%11f,%11f), P:%11f (%11f,%11f,%11f)\n",xTrk1Ori,yTrk1Ori,zTrk1Ori,GetPaFromPxPyPz(MomTrk1),MomTrk1[0],MomTrk1[1],MomTrk1[2]);
	printf("### Trk2Ori XYZ:(%11f,%11f,%11f), P:%11f (%11f,%11f,%11f)\n",xTrk2Ori,yTrk2Ori,zTrk2Ori,GetPaFromPxPyPz(MomTrk2),MomTrk2[0],MomTrk2[1],MomTrk2[2]);
	printf("### Trk1Prp XYZ:(%11f,%11f,%11f), P:%11f (%11f,%11f,%11f)\n",vTrk1.X(),vTrk1.Y(),zTrk1,GetPaFromPxPyPz(MomPos),MomPos[0],MomPos[1],MomPos[2]);
	printf("### Trk2Prp XYZ:(%11f,%11f,%11f), P:%11f (%11f,%11f,%11f)\n",vTrk2.X(),vTrk2.Y(),zTrk2,GetPaFromPxPyPz(MomNeg),MomNeg[0],MomNeg[1],MomNeg[2]);
	printf("### Vtx     XYZ:(%11f,%11f,%11f), P:%11f (%11f,%11f,%11f)\n",PosV0[0],PosV0[1],PosV0[2],GetPaFromPxPyPz(MomV0),MomV0[0],MomV0[1],MomV0[2]);
	printf("##########################################################################################################\n");
*/
	return kTRUE;

}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskNOmegaLPK::InvMassLambda(Double_t MomPos[3], Double_t MomNeg[3], AliESDtrack *pos, AliESDtrack *neg, Double_t v0Return[1]) {

	//------------------------------------------------------------------------------------------
	// version 1.11 (2015/10/23)
	// Input:  MomPos & MomNeg: Array of momentum (positive and negative tracks) (0:Px, 1:Py, 2:Pz)
	//         pos & neg: AliESDtrack information for each track that has MomPos/MomNeg momentum
	// Return: Invariant mass of Lambda (Proton(+/-) + Pion(-/+))
	//         v0Return[0]: 1:Lambda, -1:LambdaBar, 0:Not Lambda(Bar)
	//------------------------------------------------------------------------------------------

	// Set cut parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Double_t cutProtonMomTPC = 999.;
	Double_t cutPionMomTPC   = 1.5;
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
//	if (posTPCProton>fSigmaTPC && posTPCPion>fSigmaTPC) return 0.;
	if (posTPCProton>fSigmaTPC) return 0.;

	// For negative track
	Double_t negTPCProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kProton  ));
	Double_t negTPCPion     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kPion    ));
//	if (negTPCProton>fSigmaTPC && negTPCPion>fSigmaTPC) return 0.;
	if (negTPCPion>fSigmaTPC) return 0.;
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
//		if (posTOFProton>fSigmaTOF && posTOFPion>fSigmaTOF) return 0.;
		if (posTOFProton>fSigmaTOF) return 0.;
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
//		if (negTOFProton>fSigmaTOF && negTOFPion>fSigmaTOF) return 0.;
		if (negTOFPion>fSigmaTOF) return 0.;
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
		if (posTOFProton<fSigmaTOF) posProton = kTRUE;
		if (posTOFPion<fSigmaTOF)   posPion   = kTRUE;
	} else {// PID with TPC information
		if (posTPCProton<fSigmaTPC && posEtp.GetP()<cutProtonMomTPC) posProton = kTRUE;
		if (posTPCPion<fSigmaTPC   && posEtp.GetP()<cutPionMomTPC)   posPion   = kTRUE;
	}
	if (negTOFOn) {// PID with TOF information
		if (negTOFProton<fSigmaTOF) negProton = kTRUE;
		if (negTOFPion<fSigmaTOF)   negPion   = kTRUE;
	} else {// PID with TPC information
		if (negTPCProton<fSigmaTPC && negEtp.GetP()<cutProtonMomTPC) negProton = kTRUE;
		if (negTPCPion<fSigmaTPC   && negEtp.GetP()<cutPionMomTPC)   negPion   = kTRUE;
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
//	if (Lambda) {
		vPos.SetXYZM(MomPos[0],MomPos[1],MomPos[2],mProtonPDG);
		vNeg.SetXYZM(MomNeg[0],MomNeg[1],MomNeg[2],mPionPDG);
//	} else if (LambdaBar) {
//		vPos.SetXYZM(MomPos[0],MomPos[1],MomPos[2],mPionPDG);
//		vNeg.SetXYZM(MomNeg[0],MomNeg[1],MomNeg[2],mProtonPDG);
//	} else {
//		vPos.SetXYZM(0.,0.,0.,0.);
//		vNeg.SetXYZM(0.,0.,0.,0.);
//	}
	vLambda = vPos + vNeg;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  return vLambda.M();

}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskNOmegaLPK::InvMassLambdaStar(Double_t MomPos[3], Double_t MomNeg[3], AliESDtrack *pos, AliESDtrack *neg, Double_t v0Return[1]) {

	//------------------------------------------------------------------------------------------
	// version 1.11 (2015/10/23)
	// Input:  MomPos & MomNeg: Array of momentum (positive and negative tracks) (0:Px, 1:Py, 2:Pz)
	//         pos & neg: AliESDtrack information for each track that has MomPos/MomNeg momentum
	// Return: Invariant mass of Lambda(1520) (Proton(+/-) + Kaon(-/+))
	//         v0Return[0]: 1:LambdaStar, -1:LambdaStarBar, 0:Not LambdaStar(Bar)
	//------------------------------------------------------------------------------------------

	// Set cut parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	Double_t cutProtonMomTPC = 1.1;
	Double_t cutProtonMomTPC = 999.;
	Double_t cutKaonMomTPC   = 0.6;
//	Double_t cutKaonMomTPC   = 999.;;
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
//	if (posTPCProton>fSigmaTPC && posTPCKaon>fSigmaTPC) return 0.;
	if (posTPCProton>fSigmaTPC) return 0.;

	// For negative track
	Double_t negTPCProton   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kProton  ));
	Double_t negTPCKaon     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kKaon    ));
//	if (negTPCProton>fSigmaTPC && negTPCKaon>fSigmaTPC) return 0.;
	if (negTPCKaon>fSigmaTPC) return 0.;
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
//		if (posTOFProton>fSigmaTOF && posTOFKaon>fSigmaTOF) return 0.;
		if (posTOFProton>fSigmaTOF) return 0.;
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
//		if (negTOFProton>fSigmaTOF && negTOFKaon>fSigmaTOF) return 0.;
		if (negTOFKaon>fSigmaTOF) return 0.;
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
		if (posTOFProton<fSigmaTOF) posProton = kTRUE;
		if (posTOFKaon<fSigmaTOF)   posKaon   = kTRUE;
	} else {// PID with TPC information
		if (posTPCProton<fSigmaTPC && posEtp.GetP()<cutProtonMomTPC) posProton = kTRUE;
		if (posTPCKaon<fSigmaTPC   && posEtp.GetP()<cutKaonMomTPC)   posKaon   = kTRUE;
	}
	if (negTOFOn) {// PID with TOF information
		if (negTOFProton<fSigmaTOF) negProton = kTRUE;
		if (negTOFKaon<fSigmaTOF)   negKaon   = kTRUE;
	} else {// PID with TPC information
		if (negTPCProton<fSigmaTPC && negEtp.GetP()<cutProtonMomTPC) negProton = kTRUE;
		if (negTPCKaon<fSigmaTPC   && negEtp.GetP()<cutKaonMomTPC)   negKaon   = kTRUE;
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
//	if (LambdaStar) {
		vPos.SetXYZM(MomPos[0],MomPos[1],MomPos[2],mProtonPDG);
		vNeg.SetXYZM(MomNeg[0],MomNeg[1],MomNeg[2],mKaonPDG);
//	} else if (LambdaStarBar) {
//		vPos.SetXYZM(MomPos[0],MomPos[1],MomPos[2],mKaonPDG);
//		vNeg.SetXYZM(MomNeg[0],MomNeg[1],MomNeg[2],mProtonPDG);
//	} else {
//		vPos.SetXYZM(0.,0.,0.,0.);
//		vNeg.SetXYZM(0.,0.,0.,0.);
//	}
	vLambdaStar = vPos + vNeg;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  return vLambdaStar.M();

}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskNOmegaLPK::GetPaFromPxPyPz(Double_t Momentum[3]) { 

	//------------------------------------------------------------------------------------------
	// version 0.00 (2015/07/31)
	// Input:  Momentum: Array of momentum (0:Px, 1:Py, 2:Pz)
	// Return: Norm of momentum (Sqrt(Px*Px+Py*Py+Pz*Pz))
	//------------------------------------------------------------------------------------------

	return TMath::Sqrt( Momentum[0]*Momentum[0] +	Momentum[1]*Momentum[1] + Momentum[2]*Momentum[2] );

}
//______________________________________________________________________________________________________
Double_t AliAnalysisTaskNOmegaLPK::GetPtFromPxPyPz(Double_t Momentum[3]) { 

	//------------------------------------------------------------------------------------------
	// version 0.00 (2015/07/31)
	// Input:  Momentum: Array of momentum (0:Px, 1:Py, 2:Pz)
	// Return: Pt (Sqrt(Px*Px+Py*Py))
	//------------------------------------------------------------------------------------------

	return TMath::Sqrt( Momentum[0]*Momentum[0] +	Momentum[1]*Momentum[1] );

}
//______________________________________________________________________________________________________
//__________________________________________________________________________
