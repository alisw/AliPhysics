#include "TChain.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskLambdaNNRun2.h"
#include "AliPIDResponse.h"
#include "AliInputEventHandler.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliExternalTrackParam.h"
#include "AliESDv0.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"
#include <iostream>

class AliAnalysisTaskLambdaNNRun2;
using namespace std;

ClassImp(AliAnalysisTaskLambdaNNRun2)

const double triton_mass = 2.808921;
const double pion_mass = 0.13957;

AliAnalysisTaskLambdaNNRun2::AliAnalysisTaskLambdaNNRun2() : AliAnalysisTaskSE(),
fOutputTree(0), fMC(kTRUE), fAOD(0), fMCEvent(0), fEventCut(0), fPID(0), fOutputList(0)
{}
//_____________________________________________________________________________
AliAnalysisTaskLambdaNNRun2::AliAnalysisTaskLambdaNNRun2(const char* name, Bool_t isMC) : AliAnalysisTaskSE(name),
fOutputTree(0), fMC(isMC), fAOD(0), fMCEvent(0), fEventCut(0), fPID(0), fOutputList(0)
{
	// constructor
	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
	DefineOutput(2, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskLambdaNNRun2::~AliAnalysisTaskLambdaNNRun2()
{
	// destructor
	if (fOutputList) {
		delete fOutputList;
	}
}
//_____________________________________________________________________________
void AliAnalysisTaskLambdaNNRun2::UserCreateOutputObjects()
{
	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);

	fOutputEvent = new AnalysisEvent();
	fOutputTree = new TTree("OutputTree", "Output Tree");
	fOutputTree->Branch("event", &fOutputEvent);

	fEventCut.AddQAplotsToList(fOutputList);

	PostData(1, fOutputList);
	PostData(2, fOutputTree);
}
//_____________________________________________________________________________
void AliAnalysisTaskLambdaNNRun2::UserExec(Option_t *)
{
	// For the output tree
	fOutputEvent->Init();

	// Event selection via AliEventCuts
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!fAOD) return;

	TClonesArray* AODMCParticlesArray = NULL;
	if(fMC){
		// MCEvent
		fMCEvent = MCEvent();
		if (!fMCEvent) return;

		// process MC particles
		AODMCParticlesArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
		if (AODMCParticlesArray == NULL) return;
	}

	bool isEventGood = fEventCut.AcceptEvent(fAOD);
	if (!isEventGood) {
		PostData(1, fOutputList);
		return;
	}
	fOutputEvent->RunNumber = fAOD->GetRunNumber();

	// Centrality selection in Pb-Pb uses the percentile from V0
	fOutputEvent->Centrality = fEventCut.GetCentrality();

	// Get PID
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
	fPID = handl->GetPIDResponse();
	if (!fPID) return;

	// V0 loop
	Double_t vertex[3] = { -100.0, -100.0, -100.0 };
	const AliAODVertex *vertexAOD = fAOD->GetPrimaryVertex();
	vertexAOD->GetXYZ(vertex);
	for (Int_t ivertex = 0; ivertex < fAOD->GetNumberOfV0s(); ivertex++) {
		AliAODv0 * v0 = fAOD->GetV0(ivertex);

		AliAODTrack * track0 = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0)); //pos
		AliAODTrack * track1 = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1)); //neg

		if ( fAOD->GetMagneticField() * ( track0->Px()*track1->Py() - track0->Py()*track1->Px() ) < 0 )
		  fAnalysis_V0.isCowboy = kTRUE;
		else
		  fAnalysis_V0.isCowboy = kFALSE;
		Int_t label0 = track0->GetLabel();
		Int_t label1 = track1->GetLabel();

		if (!track0 || !track1) continue;

		// Cut on filter bit
		if (track0->TestFilterBit(1) == false || track1->TestFilterBit(1) == false) continue;

		// Loose track cuts
		if (LooseTrackCuts(track0) == false || LooseTrackCuts(track1) == false) continue;

		// Loose cut on CosPointingAngle
		if (v0->CosPointingAngle(vertex) < 0.95) continue;

		// Loose cut on dca
		if (v0->DcaV0Daughters() > 1.0) continue;

		// Online Vertex finder
		//if (!v0->GetOnFlyStatus()) continue;

		////////////////////////////
		fAnalysis_V0.Reset();
		fAnalysis_V0.onlineV0 = v0->GetOnFlyStatus();
		fAnalysis_V0.dca = v0->DcaV0Daughters();
		fAnalysis_V0.cosPointing = v0->CosPointingAngle(vertex);
		fAnalysis_V0.decayRadius = v0->DecayLengthV0(vertex);

		Int_t pdgcode0 = 0;
		Int_t pdgcode1 = 0;
		Int_t pdgcodeMother0 = 0;
		Int_t pdgcodeMother1 = 0;
		Int_t labelMother0 = 0;
		Int_t labelMother1 = 0;

		if(fMC){
			AliAODMCParticle* mcTrack0 = static_cast<AliAODMCParticle*>(AODMCParticlesArray->At(TMath::Abs(label0)));
			AliAODMCParticle* mcTrack1 = static_cast<AliAODMCParticle*>(AODMCParticlesArray->At(TMath::Abs(label1)));
			labelMother0 = mcTrack0->GetMother();
			labelMother1 = mcTrack1->GetMother();
			AliAODMCParticle* mcMother0 = static_cast<AliAODMCParticle*>(AODMCParticlesArray->At(TMath::Abs(labelMother0)));
			AliAODMCParticle* mcMother1 = static_cast<AliAODMCParticle*>(AODMCParticlesArray->At(TMath::Abs(labelMother1)));
			pdgcode0 = mcTrack0->GetPdgCode();
			pdgcode1 = mcTrack1->GetPdgCode();
			pdgcodeMother0 = mcMother0->GetPdgCode();
			pdgcodeMother1 = mcMother1->GetPdgCode();
		}

		Bool_t isTriton[2] = {kFALSE, kFALSE};
		isTriton[0] = fabs(fPID->NumberOfSigmasTPC(track0, AliPID::kTriton)) < 4.;
		isTriton[1] = fabs(fPID->NumberOfSigmasTPC(track1, AliPID::kTriton)) < 4.;

		// // for p (measured by TPC) > 1.5 I also look at TOF
		// if (isTriton[0] && track0->GetTPCmomentum() > 1.5) {
		// 	if ( fabs(fPID->NumberOfSigmasTOF(track0, AliPID::kTriton)) >= 4.) continue;
		// }
		// else if (isTriton[1] && track1->GetTPCmomentum() > 1.5) {
		// 	if ( fabs(fPID->NumberOfSigmasTOF(track1, AliPID::kTriton)) >= 4.) continue;
		// }

		Bool_t isPion[2] = {kFALSE, kFALSE};
		isPion[0] = fabs(fPID->NumberOfSigmasTPC(track0, AliPID::kPion)) < 4.;
		isPion[1] = fabs(fPID->NumberOfSigmasTPC(track1, AliPID::kPion)) < 4.;

		// AntiLambdaNN t- p+ ////////
		if (isTriton[0] == kTRUE &&
				track0->Charge() < 0 &&
		                track0->GetTPCNcls()>=90 &&
				isPion[1] == kTRUE &&
				track1->Charge() > 0) {
					if(!fMC)
						FillEvent(AnalysisV0::Type::AntiLambdaNN, track0, track1);
					else if (	pdgcode0 == -1000010030 			&&
										pdgcode1 ==  211 							&&
										labelMother0 == labelMother1 	&&
										pdgcodeMother0 == -1010000020	&&
										pdgcodeMother1 == -1010000020	)
									FillEvent(AnalysisV0::Type::AntiLambdaNN, track0, track1);
				}

		// AntiLambdaNN p+ t- ////////
		if (isTriton[1] == kTRUE &&
				track1->Charge() < 0 &&
   		                track1->GetTPCNcls()>=90 &&
				isPion[0] == kTRUE &&
				track0->Charge() > 0) {
				if(!fMC)
					FillEvent(AnalysisV0::Type::AntiLambdaNN, track1, track0);
				else if (	pdgcode0 ==  211 							&&
									pdgcode1 == -1000010030 			&&
									labelMother0 == labelMother1 	&&
									pdgcodeMother0 == -1010000020	&&
									pdgcodeMother1 == -1010000020	)
								FillEvent(AnalysisV0::Type::AntiLambdaNN, track1, track0);
			}

		// LambdaNN d+ p- ////////
		if (isTriton[0] == kTRUE &&
				track0->Charge() > 0 &&
		                track0->GetTPCNcls()>=90 &&
				isPion[1] == kTRUE &&
				track1->Charge() < 0) {
				if(!fMC)
					FillEvent(AnalysisV0::Type::LambdaNN, track0, track1);
				else if (	pdgcode0 ==  1000010030 			&&
									pdgcode1 == -211 							&&
									labelMother0 == labelMother1 	&&
									pdgcodeMother0 ==  1010000020	&&
									pdgcodeMother1 ==  1010000020	)
								FillEvent(AnalysisV0::Type::LambdaNN, track0, track1);
			}

		// LambdaNN p- d+ ////////
		if (isTriton[1] == kTRUE &&
				track1->Charge() > 0 &&
		                track1->GetTPCNcls()>=90 &&
				isPion[0] == kTRUE &&
				track0->Charge() < 0) {
				if(!fMC)
					FillEvent(AnalysisV0::Type::LambdaNN, track1, track0);
				else if (	pdgcode0 == -211 							&&
									pdgcode1 ==  1000010030 			&&
									labelMother0 == labelMother1 	&&
									pdgcodeMother0 ==  1010000020	&&
									pdgcodeMother1 ==  1010000020	)
								FillEvent(AnalysisV0::Type::LambdaNN, track1, track0);
			}

	}

	if (fOutputEvent->Size() == 0) return;
	fOutputEvent->NV0s = fOutputEvent->Size();
	fOutputTree->Fill();
	PostData(2, fOutputTree);

}
//_____________________________________________________________________________
bool AliAnalysisTaskLambdaNNRun2::LooseTrackCuts(AliAODTrack *track) {

	// |eta| < 0.8
	if (fabs(track->Eta()) >= 0.8) return false;

	// TPC clusters > 70 (80 analysis cut)
	if (track->GetTPCNcls() <= 70) return false;

	// chi2 per TPC clusters < 6 (5 analysis cut)                                                                                                                                                
        if (track->Chi2perNDF() >= 6) return false;

	// Max signal in TPC
	if (track->GetTPCNcls() >= 1000) return false;

	// Kink daughter reject
	AliAODVertex *vtx = (AliAODVertex*)track->GetProdVertex();
	if (!vtx) return false;
	if (vtx->GetType() == AliAODVertex::kKink) return false;

	// TPC refit true
	ULong64_t status = track->GetStatus();
	if (!(status & AliAODTrack::kTPCrefit)) return false;

	return true;
}
//_____________________________________________________________________________
void AliAnalysisTaskLambdaNNRun2::FillEvent(AnalysisV0::Type etype, AliAODTrack* triton, AliAODTrack* pion) {
	TLorentzVector triton_v4, pion_v4, LambdaNN_v4;
	triton_v4.SetPxPyPzE(triton->Px(), triton->Py(), triton->Pz(), TMath::Sqrt(triton_mass * triton_mass + triton->P() * triton->P()));
	pion_v4.SetPxPyPzE(pion->Px(), pion->Py(), pion->Pz(), TMath::Sqrt(pion_mass * pion_mass + pion->P() * pion->P()));
	LambdaNN_v4 = triton_v4 + pion_v4;

	if (LambdaNN_v4.M() > 2.8 && LambdaNN_v4.M() < 4.2 && triton->P()<3.) {
		fAnalysis_V0.mass = LambdaNN_v4.M();
		
		fAnalysis_V0.topology = etype;
		fAnalysis_V0.triton_Charge = triton->Charge();
		fAnalysis_V0.triton_P = triton->P();
		fAnalysis_V0.triton_PTPC = triton->GetTPCmomentum();
		fAnalysis_V0.triton_Px = triton->Px();
		fAnalysis_V0.triton_Py = triton->Py();
		fAnalysis_V0.triton_Pz = triton->Pz();
		fAnalysis_V0.triton_Eta = triton->Eta();
		fAnalysis_V0.triton_NSigmaTPC = fPID->NumberOfSigmasTPC(triton, AliPID::kTriton);
		fAnalysis_V0.triton_NSigmaTOF = fPID->NumberOfSigmasTOF(triton, AliPID::kTriton);
		fAnalysis_V0.triton_TPCsignal = triton->GetTPCsignal();
		fAnalysis_V0.triton_TPCNcls = triton->GetTPCNcls();
		fAnalysis_V0.triton_Chi2perNDF = triton->Chi2perNDF();
		fAnalysis_V0.triton_TPCExpSignal = fPID->GetTPCResponse().GetExpectedSignal(triton->GetTPCmomentum(), AliPID::kTriton);
		fAnalysis_V0.triton_TOFExpSignal = fPID->GetTOFResponse().GetExpectedSignal(triton, AliPID::kTriton);

		fAnalysis_V0.pion_Charge = pion->Charge();
		fAnalysis_V0.pion_P = pion->P();
		fAnalysis_V0.pion_PTPC = pion->GetTPCmomentum();
		fAnalysis_V0.pion_Px = pion->Px();
		fAnalysis_V0.pion_Py = pion->Py();
		fAnalysis_V0.pion_Pz = pion->Pz();
		fAnalysis_V0.pion_Eta = pion->Eta();
		fAnalysis_V0.pion_NSigmaTPC = fPID->NumberOfSigmasTPC(pion, AliPID::kPion);
		fAnalysis_V0.pion_NSigmaTOF = fPID->NumberOfSigmasTOF(pion, AliPID::kPion);
		fAnalysis_V0.pion_TPCsignal = pion->GetTPCsignal();
		fAnalysis_V0.pion_TPCNcls = pion->GetTPCNcls();
		fAnalysis_V0.pion_Chi2perNDF = pion->Chi2perNDF();
		fAnalysis_V0.pion_TPCExpSignal = fPID->GetTPCResponse().GetExpectedSignal(pion->GetTPCmomentum(), AliPID::kPion);
		fAnalysis_V0.triton_TOFExpSignal = fPID->GetTOFResponse().GetExpectedSignal(pion, AliPID::kPion);

		fOutputEvent->AddV0(fAnalysis_V0);
	}
}
//_____________________________________________________________________________
void AliAnalysisTaskLambdaNNRun2::Terminate(Option_t *)
{}
//_____________________________________________________________________________
