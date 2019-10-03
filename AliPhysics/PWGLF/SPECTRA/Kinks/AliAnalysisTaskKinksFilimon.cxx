#include "AliAnalysisTaskKinksFilimon.h"
#include <limits>
#include <TCanvas.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TMCProcess.h>
#include <TPDGCode.h>
#include <TTree.h>
#include <TDirectory.h>
#include <AliAnalysisManager.h>
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliESDtrackCuts.h>
#include <AliCFParticleGenCuts.h>
#include <AliPIDResponse.h>
#include <AliCFParticleGenCuts.h>
#include <AliESDpidCuts.h>
#include <AliESDkink.h>
#include <AliStack.h>
#include <AliMCParticle.h>
#include <AliMCVertex.h>
#include <AliCFTrackKineCuts.h>
#include <AliKineTrackCuts.h>
#include <AliAODEvent.h>
#include <AliAODInputHandler.h>
#include <AliAODMCHeader.h>
#include <AliAODTrack.h>
#include <AliAODMCParticle.h>
#include <AliMultiplicity.h>
#include <AliMixInputEventHandler.h>
#include <AliLog.h>

// Kink and resonance analysis task
// Author: Filimon Roukoutakis, University of Athens

ClassImp(AliAnalysisTaskKinksFilimon)

//________________________________________________________________________
AliAnalysisTaskKinksFilimon::AliAnalysisTaskKinksFilimon(const char *name) 
  : AliAnalysisTaskSE(name), fOutputContClass(0), fOutputCont(0), fMainInputHandler(0), fIsAOD(kFALSE), fUseMC(kFALSE), fFillCutHist(kFALSE), fEventSelectionCutsMC(kFALSE), fCollisionType(kPP), fOfflineTriggerType(AliVEvent::kMB), fPartIdDedx(AliPID::kElectron), fPIDResponse(0), fRecEventCuts(0), fRecCentCuts(0), fCommonKineTrackCuts(0), fEsdTrackCutsKinkMother(0), fEsdTrackCutsKinkDaughter(0), fEsdTrackCutsPartner(0), fEsdPIDCutsArray(), fRecKinkCutsKaon(0), fAODFilterMap(0), /*fRecKinkCutsPion(0),*/ /*fMcKinkCutsPion(0),*/ fHistPtYTemplate(0), fHistRPhiZTemplate(0), fKaonMass(0), fMuonMass(0), fPionMass(0), fProtonMass(0), fElectronMass(0), fPhiMassRangeMin(0), fPhiMassRangeMax(0), fKstarMassRangeMin(0), fKstarMassRangeMax(0), fLambdaMassRangeMin(0), fLambdaMassRangeMax(0),
  fhRecMultUnbiased(0), fhRecMult(0), fhRecMultCentralityUnbiased(0), fhRecMultCentrality(0), fhRecAllPtY(0), fhESDnKinks(0), fhESDKinkQt(0), fhESDKinkAngle(0), fhESDKinkPAngle(0), fhESDKinkDCA(0), fhESDAllKinkPtY(0), fhESDUnknownKinkPtY(0), fhESDKaonKinkPtY(0), fhESDPionKinkPtY(0), fhESDKaonKinkPtEta(0), fhESDKPlusKinkPtY(0), fhESDKMinusKinkPtY(0), fhESDKPlusKinkPt(0), fhESDKMinusKinkPt(0), fhESDPiPlusKinkPtY(0), fhESDPiMinusKinkPtY(0), fhESDKaonKinkPtYCent(0), fhESDPionKinkPtYCent(0), fhESDKPlusKinkPtYCent(0), fhESDPiPlusKinkPtYCent(0), fhESDKMinusKinkPtYCent(0), fhESDPiMinusKinkPtYCent(0), fhRecPrimaryVertexRPhiZ(0), fhESDAllKinkDecaysRPhiZ(0), fhESDKaonKinkDecaysRPhiZ(0), fhESDPionKinkDecaysRPhiZ(0), fhESDKaonKinkQt(0), fhESDKaonKinkAngle(0), fhESDPionKinkQt(0), fhESDPionKinkAngle(0), fhESDMomentumTPCSignalKaonKinks(0), fhESDKaonKinksMotherAndDaughterTPCncls(0), fhESDKaonKinksMotherVSDaughterTPCncls(0), fhESDMomentumTPCSignalPionKinks(0), fhESDPionKinksMotherAndDaughterTPCncls(0), fhESDPionKinksMotherVSDaughterTPCncls(0), fhESDKaonOverPionEbE(0), fhESDInvMassKinkDecayMuNu(0),
  fhMCmainVertexRPhiZ(0), fhMCKaonKinkDecaysRPhiZ(0), fhMCmult(0), fhMCmultPrim(0), fhMCPtYall(0), fhMCPtYprim(0), fhMCpdg(0), fhMCPtYprimKaon(0), fhMCPtYprimKPlus(0), fhMCPtYprimKMinus(0), fhMCPtprimKPlus(0), fhMCPtprimKMinus(0), fhMCPtYCentprimKaon(0), fhMCPtYCentprimKPlus(0), fhMCPtYCentprimKMinus(0), fhMCPtYkinkKaonFiducial(0), fhMCPtYkinkKPlusFiducial(0), fhMCPtYkinkKMinusFiducial(0), fhMCKinkKaonQt(0), fhMCKinkKaonAngle(0), fhMCKaonKinkDecaysTrackLength(0), fhMCKaonKinkDecaysTrackTime(0), fhMCkaonMotherPdg(0), fhMCpionMotherPdg(0), fhMCPhiPtY(0), fhMCPtYprimPhi(0), fhMCPhiInvMassPtYtest(0), fhMCInvMassPtprimPhi(0), fhMCInvMassPtCentprimPhi(0), fhMCKstarPtY(0), fhMCPtYprimKstar(0), fhMCInvMassPtprimKstar(0), fhMCInvMassPtCentprimKstar(0), fhMCLambdaPtY(0), fhMCPtYprimLambda(0), fhMCInvMassPtprimLambda(0), fhMCInvMassPtCentprimLambda(0), fhMCPhi2KaonPtY(0), fhMCKstar2KaonPtY(0), fhMCKaon2MuonPtY(0), fhMCKaon2PionPtY(0), fhMCPtYprimPion(0), fhMCPtYprimPiPlus(0), fhMCPtYprimPiMinus(0), fhMCPtYCentprimPion(0), fhMCPtYCentprimPiPlus(0), fhMCPtYCentprimPiMinus(0), fhMCPionKinkDecaysRPhiZ(0), fhMCPion2MuonPtY(0), fhMCPtYkinkPionFiducial(0), fhMCPtYkinkPiPlusFiducial(0), fhMCPtYkinkPiMinusFiducial(0), fhMCKinkPionQt(0), fhMCKinkPionAngle(0), fhMCMotherDaughterPdg(0), fhMCKaonOverPionEbE(0), fhMCPhiDecayOpeningAngle(0), fhMCKstarDecayOpeningAngle(0), fhMCLambdaDecayOpeningAngle(0), 
  fhESDKinkQtTrueMC(0), fhESDKinkAngleTrueMC(0), fhESDKaonKinkPtYTrueMC(0), fhESDKaonKinkPtYFakeMC(0), fhESDKaonKinkPtYTrueMCsecondary(0), fhESDKPlusKinkPtYTrueMC(0), fhESDKPlusKinkPtYFakeMC(0), fhESDKPlusKinkPtYTrueMCsecondary(0), fhESDKMinusKinkPtYTrueMC(0), fhESDKMinusKinkPtYFakeMC(0), fhESDKMinusKinkPtYTrueMCsecondary(0), fhESDPionKinkPtYTrueMC(0), fhESDPionKinkPtYFakeMC(0), fhESDPionKinkPtYTrueMCsecondary(0), fhESDPiPlusKinkPtYTrueMC(0), fhESDPiPlusKinkPtYFakeMC(0), fhESDPiPlusKinkPtYTrueMCsecondary(0), fhESDPiMinusKinkPtYTrueMC(0), fhESDPiMinusKinkPtYFakeMC(0), fhESDPiMinusKinkPtYTrueMCsecondary(0), fhESDPhi2KaonPtYTrueMC(0), fhESDPhiPtYtestTrueMC(0), fhESDPhiPtYtestFakeMC(0), fhESDPhiPtYtestLikeSignTrueMC(0), fhESDPhiPtYtestLikeSignFakeMC(0), fhESDKstar2KaonPtYTrueMC(0), fhESDKstarPtYtestTrueMC(0), fhESDKstarPtYtestFakeMC(0), fhESDKstarPtYtestLikeSignTrueMC(0), fhESDKstarPtYtestLikeSignFakeMC(0), fhESDLambda2KaonPtYTrueMC(0),
  fhAODnVertices(0), fhAODnKinkVertices(0), fhAODnKinkDaughters(0), fhAODAllKinkDecaysRPhiZ(0), fhAODKinkKaonQt(0), fhAODKinkKaonAngle(0), fhAODKinkMothervsDaughterPIDcheck(0), fhAODKaonKinkPID(0), fhAODKaonKinkPIDFake(0), fhAODmotherVertexRPhiZ(0), fhAODmotherVertexDCARPhiZ(0),
  fhRecPhiInvMassPtYtest(0), fhRecPhiInvMassPtYtestLikeSign(0), fhRecPhiPtYtest(0), fhRecPhiPtYtestLikeSign(0), fhRecPhiInvMassPt(0), fhRecPhiInvMassPtLikeSign(0), fhRecPhiInvMassBothKinks(0), fhRecPhiInvMassBothKinksLikeSign(0), fhRecKstarPtYtest(0), fhRecKstarPtYtestLikeSign(0), fhRecKstarInvMassPt(0), fhRecKstarInvMassPtLikeSign(0), fhRecKstarInvMassBothKinks(0), fhRecKstarInvMassBothKinksLikeSign(0), fhRecLambdaPtYtest(0), fhRecLambdaPtYtestLikeSign(0), fhRecLambdaInvMassPt(0), fhRecLambdaInvMassPtLikeSign(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskKinksFilimon::AliAnalysisTaskKinksFilimon(const char *name, const AliVEvent::EOfflineTriggerTypes offlineTriggerType, AliESDkinkCuts* esdKinkCuts, AliESDtrackCuts* esdTrackCutsKinkMother, AliESDtrackCuts* esdTrackCutsKinkDaughter, AliESDtrackCuts* esdTrackCutsPartner, TArrayF esdPIDCutsResonances, /*const AliCFParticleGenCuts* genTrackCuts, const TH1* histPtTemplate=0x0,*/ const Bool_t useMC, const Bool_t fillCutHist, const TH2* histPtYTemplate, const THnSparseF* histRPhiZTemplate, const ECollisionType collisionType, TClass* outputContClass, const UInt_t aodFilterMap)
  : AliAnalysisTaskSE(name), fOutputContClass(outputContClass), fOutputCont(0), fMainInputHandler(0), fIsAOD(kFALSE), fUseMC(useMC), fFillCutHist(fillCutHist), fEventSelectionCutsMC(kFALSE), fCollisionType(collisionType), fOfflineTriggerType(offlineTriggerType), fPartIdDedx(AliPID::kElectron), fPIDResponse(0), fRecEventCuts(0), fRecCentCuts(0), fCommonKineTrackCuts(0), fEsdTrackCutsKinkMother(esdTrackCutsKinkMother), fEsdTrackCutsKinkDaughter(esdTrackCutsKinkDaughter), fEsdTrackCutsPartner(esdTrackCutsPartner), fEsdPIDCutsArray(esdPIDCutsResonances), fRecKinkCutsKaon(esdKinkCuts), fAODFilterMap(aodFilterMap), /*fRecKinkCutsPion(0),*/ /*fMcKinkCutsPion(0),*/ fHistPtYTemplate(histPtYTemplate), fHistRPhiZTemplate(histRPhiZTemplate), fKaonMass(0), fMuonMass(0), fPionMass(0), fProtonMass(0), fElectronMass(0), fPhiMassRangeMin(0), fPhiMassRangeMax(0), fKstarMassRangeMin(0), fKstarMassRangeMax(0), fLambdaMassRangeMin(0), fLambdaMassRangeMax(0),
  fhRecMultUnbiased(0), fhRecMult(0), fhRecMultCentralityUnbiased(0), fhRecMultCentrality(0), fhRecAllPtY(0), fhESDnKinks(0), fhESDKinkQt(0), fhESDKinkAngle(0), fhESDKinkPAngle(0), fhESDKinkDCA(0), fhESDAllKinkPtY(0), fhESDUnknownKinkPtY(0), fhESDKaonKinkPtY(0), fhESDPionKinkPtY(0), fhESDKaonKinkPtEta(0), fhESDKPlusKinkPtY(0), fhESDKMinusKinkPtY(0), fhESDKPlusKinkPt(0), fhESDKMinusKinkPt(0), fhESDPiPlusKinkPtY(0), fhESDPiMinusKinkPtY(0), fhESDKaonKinkPtYCent(0), fhESDPionKinkPtYCent(0), fhESDKPlusKinkPtYCent(0), fhESDPiPlusKinkPtYCent(0), fhESDKMinusKinkPtYCent(0), fhESDPiMinusKinkPtYCent(0), fhRecPrimaryVertexRPhiZ(0), fhESDAllKinkDecaysRPhiZ(0), fhESDKaonKinkDecaysRPhiZ(0), fhESDPionKinkDecaysRPhiZ(0), fhESDKaonKinkQt(0), fhESDKaonKinkAngle(0), fhESDPionKinkQt(0), fhESDPionKinkAngle(0), fhESDMomentumTPCSignalKaonKinks(0), fhESDKaonKinksMotherAndDaughterTPCncls(0), fhESDKaonKinksMotherVSDaughterTPCncls(0), fhESDMomentumTPCSignalPionKinks(0), fhESDPionKinksMotherAndDaughterTPCncls(0), fhESDPionKinksMotherVSDaughterTPCncls(0), fhESDKaonOverPionEbE(0), fhESDInvMassKinkDecayMuNu(0),
  fhMCmainVertexRPhiZ(0), fhMCKaonKinkDecaysRPhiZ(0), fhMCmult(0), fhMCmultPrim(0), fhMCPtYall(0), fhMCPtYprim(0), fhMCpdg(0), fhMCPtYprimKaon(0), fhMCPtYprimKPlus(0), fhMCPtYprimKMinus(0), fhMCPtprimKPlus(0), fhMCPtprimKMinus(0), fhMCPtYCentprimKaon(0), fhMCPtYCentprimKPlus(0), fhMCPtYCentprimKMinus(0), fhMCPtYkinkKaonFiducial(0), fhMCPtYkinkKPlusFiducial(0), fhMCPtYkinkKMinusFiducial(0), fhMCKinkKaonQt(0), fhMCKinkKaonAngle(0), fhMCKaonKinkDecaysTrackLength(0), fhMCKaonKinkDecaysTrackTime(0), fhMCkaonMotherPdg(0), fhMCpionMotherPdg(0), fhMCPhiPtY(0), fhMCPtYprimPhi(0), fhMCPhiInvMassPtYtest(0), fhMCInvMassPtprimPhi(0), fhMCInvMassPtCentprimPhi(0), fhMCKstarPtY(0), fhMCPtYprimKstar(0), fhMCInvMassPtprimKstar(0), fhMCInvMassPtCentprimKstar(0), fhMCLambdaPtY(0), fhMCPtYprimLambda(0), fhMCInvMassPtprimLambda(0), fhMCInvMassPtCentprimLambda(0), fhMCPhi2KaonPtY(0), fhMCKstar2KaonPtY(0), fhMCKaon2MuonPtY(0), fhMCKaon2PionPtY(0), fhMCPtYprimPion(0), fhMCPtYprimPiPlus(0), fhMCPtYprimPiMinus(0), fhMCPtYCentprimPion(0), fhMCPtYCentprimPiPlus(0), fhMCPtYCentprimPiMinus(0), fhMCPionKinkDecaysRPhiZ(0), fhMCPion2MuonPtY(0), fhMCPtYkinkPionFiducial(0), fhMCPtYkinkPiPlusFiducial(0), fhMCPtYkinkPiMinusFiducial(0), fhMCKinkPionQt(0), fhMCKinkPionAngle(0), fhMCMotherDaughterPdg(0), fhMCKaonOverPionEbE(0), fhMCPhiDecayOpeningAngle(0), fhMCKstarDecayOpeningAngle(0), fhMCLambdaDecayOpeningAngle(0), 
  fhESDKinkQtTrueMC(0), fhESDKinkAngleTrueMC(0), fhESDKaonKinkPtYTrueMC(0), fhESDKaonKinkPtYFakeMC(0), fhESDKaonKinkPtYTrueMCsecondary(0), fhESDKPlusKinkPtYTrueMC(0), fhESDKPlusKinkPtYFakeMC(0), fhESDKPlusKinkPtYTrueMCsecondary(0), fhESDKMinusKinkPtYTrueMC(0), fhESDKMinusKinkPtYFakeMC(0), fhESDKMinusKinkPtYTrueMCsecondary(0), fhESDPionKinkPtYTrueMC(0), fhESDPionKinkPtYFakeMC(0), fhESDPionKinkPtYTrueMCsecondary(0), fhESDPiPlusKinkPtYTrueMC(0), fhESDPiPlusKinkPtYFakeMC(0), fhESDPiPlusKinkPtYTrueMCsecondary(0), fhESDPiMinusKinkPtYTrueMC(0), fhESDPiMinusKinkPtYFakeMC(0), fhESDPiMinusKinkPtYTrueMCsecondary(0), fhESDPhi2KaonPtYTrueMC(0), fhESDPhiPtYtestTrueMC(0), fhESDPhiPtYtestFakeMC(0), fhESDPhiPtYtestLikeSignTrueMC(0), fhESDPhiPtYtestLikeSignFakeMC(0), fhESDKstar2KaonPtYTrueMC(0), fhESDKstarPtYtestTrueMC(0), fhESDKstarPtYtestFakeMC(0), fhESDKstarPtYtestLikeSignTrueMC(0), fhESDKstarPtYtestLikeSignFakeMC(0), fhESDLambda2KaonPtYTrueMC(0),
  fhAODnVertices(0), fhAODnKinkVertices(0), fhAODnKinkDaughters(0), fhAODAllKinkDecaysRPhiZ(0), fhAODKinkKaonQt(0), fhAODKinkKaonAngle(0), fhAODKinkMothervsDaughterPIDcheck(0), fhAODKaonKinkPID(0), fhAODKaonKinkPIDFake(0), fhAODmotherVertexRPhiZ(0), fhAODmotherVertexDCARPhiZ(0),
  fhRecPhiInvMassPtYtest(0), fhRecPhiInvMassPtYtestLikeSign(0), fhRecPhiPtYtest(0), fhRecPhiPtYtestLikeSign(0), fhRecPhiInvMassPt(0), fhRecPhiInvMassPtLikeSign(0), fhRecPhiInvMassBothKinks(0), fhRecPhiInvMassBothKinksLikeSign(0), fhRecKstarPtYtest(0), fhRecKstarPtYtestLikeSign(0), fhRecKstarInvMassPt(0), fhRecKstarInvMassPtLikeSign(0), fhRecKstarInvMassBothKinks(0), fhRecKstarInvMassBothKinksLikeSign(0), fhRecLambdaPtYtest(0), fhRecLambdaPtYtestLikeSign(0), fhRecLambdaInvMassPt(0), fhRecLambdaInvMassPtLikeSign(0)
{
  // Constructor

	if (!fEsdTrackCutsKinkMother) AliWarning("Kink mother ESD track cuts  = 0x0, skipping ESD track cuts in analysis");
	if (!fEsdTrackCutsKinkDaughter) AliInfo("Kink daughter ESD track cuts = 0x0"); //AliFatal("Kink daughter ESD track cuts required");
	if (!fEsdTrackCutsPartner) AliInfo("Resonance partner ESD track cuts = 0x0, skipping ESD resonance analysis"); //AliFatal("Resonance partner ESD track cuts required");
  if(!fHistPtYTemplate) AliFatal("histPtYTemplate required");
  if(!fOutputContClass) {
		AliInfo("outputContClass missing, using default container TList");
		fOutputContClass = TList::Class();
	}

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, fOutputContClass);
}

//________________________________________________________________________
AliAnalysisTaskKinksFilimon::~AliAnalysisTaskKinksFilimon() {

   // Clean-up the output container, but not the histograms that are put inside
   // (the container is owner and will clean-up these histograms). Protect in PROOF case.
   if (fOutputCont && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
     delete fOutputCont;
   }
}

//________________________________________________________________________
void AliAnalysisTaskKinksFilimon::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
	//TDatabasePDG::Instance()->GetParticle("K+")->PdgCode(), TDatabasePDG::Instance()->GetParticle("pi+")->PdgCode(), TDatabasePDG::Instance()->GetParticle("mu-")->PdgCode(), TDatabasePDG::Instance()->GetParticle("e-")->PdgCode();
	 fKaonMass = TDatabasePDG::Instance()->GetParticle("K+")->Mass();
	 fMuonMass = TDatabasePDG::Instance()->GetParticle("mu+")->Mass();
	 fPionMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
	 fProtonMass = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
	 fElectronMass = TDatabasePDG::Instance()->GetParticle("e-")->Mass();
	 Double_t fPhiMass = TDatabasePDG::Instance()->GetParticle("phi")->Mass();
	 Double_t fKstarMass = TDatabasePDG::Instance()->GetParticle("K*0")->Mass();
	 Double_t fLambdaMass = 1.5195;//TDatabasePDG::Instance()->GetParticle("Lambda1520")->Mass();
	 Float_t nFWHM = 1.3; // WARNING!!! Very important parameter for the resonances
	 fPhiMassRangeMin = fPhiMass-nFWHM*TDatabasePDG::Instance()->GetParticle("phi")->Width();
	 fPhiMassRangeMax = fPhiMass+nFWHM*TDatabasePDG::Instance()->GetParticle("phi")->Width();
	 fKstarMassRangeMin = fKstarMass-nFWHM*TDatabasePDG::Instance()->GetParticle("K*0")->Width();
	 fKstarMassRangeMax = fKstarMass+nFWHM*TDatabasePDG::Instance()->GetParticle("K*0")->Width();
	 fLambdaMassRangeMin = fLambdaMass-nFWHM*0.0156;//TDatabasePDG::Instance()->GetParticle("Lambda1520")->Width();
	 fLambdaMassRangeMax = fLambdaMass+nFWHM*0.0156;//TDatabasePDG::Instance()->GetParticle("Lambda1520")->Width();
	
	fMainInputHandler = fMultiInputHandler ? fMultiInputHandler->GetFirstInputEventHandler() : fInputHandler;
	fIsAOD = dynamic_cast<AliAODInputHandler*>(fMainInputHandler);
	/*!fMultiInputHandler ? dynamic_cast<AliAODInputHandler*>(fInputHandler) : dynamic_cast<AliAODInputHandler*>(fMultiInputHandler->GetFirstInputEventHandler());*/
	fUseMC = fUseMC && (fIsAOD || fMCEventHandler);
  fPIDResponse = fMainInputHandler->GetPIDResponse();
	if(!fPIDResponse) AliWarning("No PID response used");
	/*if (genTrackCuts) genTrackCuts->GetPtRange(fRecKinkCutsKaon->fMinPt, tempMax);
	else AliFatal("Primary track cuts required for MC");*/
	fCommonKineTrackCuts = new AliKineTrackCuts("fCommonKineTrackCuts", "fCommonKineTrackCuts"); // Temporarily here to understand how to use in both MC and ESD
	Float_t tempMin = 0, tempMax = 0;
	fRecEventCuts = new AliRecEventCuts("fRecEventCuts", "Reconstructed event cuts", fOfflineTriggerType, kTRUE, -1, 1e10, -10, 10, -1, 101); // Need to have this init here for the moment as AliRecEventCuts is a "local" class
	if (fEsdTrackCutsKinkMother) {
	  fEsdTrackCutsKinkMother->GetPtRange(tempMin, tempMax);
	  fCommonKineTrackCuts->SetPtRange(tempMin, tempMax);
	  fEsdTrackCutsKinkMother->GetRapRange(tempMin, tempMax);
	  fCommonKineTrackCuts->SetRapRange(tempMin, tempMax);
	  fEsdTrackCutsKinkMother->GetEtaRange(tempMin, tempMax);
	  fCommonKineTrackCuts->SetEtaRange(tempMin, tempMax);
	  fEsdTrackCutsKinkMother->GetPtRange(tempMin, tempMax);
	}
	else {
	  fCommonKineTrackCuts->SetPtRange();
	  fCommonKineTrackCuts->SetRapRange();
	  fCommonKineTrackCuts->SetEtaRange();
	}
	
	if (!fRecKinkCutsKaon) AliFatal("Rec Kaon Kink Cuts required");
	fRecKinkCutsKaon->SetPIDResponse(fPIDResponse);
	//TH1::SetDefaultSumw2();
	//TH2::SetDefaultSumw2();
	//TDirectory* histDirectory = gDirectory->mkdir("KinksFilimon");
	//histDirectory->cd();
	//fOutputCont = histDirectory->GetList();
  fOutputCont = dynamic_cast<TCollection*>(static_cast<TObject*>(fOutputContClass->New()));
	//OpenFile(1); // May save memory WARNING! not fully checked functionality
	fOutputCont->SetOwner(kTRUE);
	//TH1::AddDirectory(kTRUE);
	//TH2::AddDirectory(kTRUE);
	//THnSparseF::AddDirectory(kTRUE);
	//TDirectory::AddDirectory(kTRUE);
	
	Int_t bins[3] = {250, 40, 300};
  Double_t xmin[3] = {0., -4.0, -300.};
  Double_t xmax[3] = {500., 4.0, 300.};
  fHistRPhiZTemplate = new THnSparseF("histRPhiZTemplate", "histRPhiZTemplate", 3, bins, xmin, xmax); // Temporary internal definition
  if(!fHistRPhiZTemplate) AliFatal("fHistRPhiZTemplate required");
	
        Float_t hLoPt = fHistPtYTemplate->GetXaxis()->GetXmin(), hHiPt = fHistPtYTemplate->GetXaxis()->GetXmax()/*, hLoY = fHistPtYTemplate->GetYaxis()->GetXmin(), hHiY = fHistPtYTemplate->GetYaxis()->GetXmax()*//*, hLoEta = -1.5, hHiEta = 1.5*/;//, fLowX=0.6, fHighX=1.8;
	Int_t nBinsPt = fHistPtYTemplate->GetNbinsX(), nBinsPtInvMassSlice = nBinsPt, nBinsY = fHistPtYTemplate->GetNbinsY();/*, nBinsEta = TMath::Nint(10*(hHiEta-hLoEta))*/
	const Double_t* binsPtCArray = 0x0;
	const Double_t* binsYCArray = 0x0;
	if (fHistPtYTemplate->GetXaxis()->IsVariableBinSize()) binsPtCArray = fHistPtYTemplate->GetXaxis()->GetXbins()->GetArray();
	else {
		binsPtCArray = new Double_t[nBinsPt+1];
		GetCustomBins(nBinsPt, hLoPt, hHiPt, const_cast<Double_t* const>(binsPtCArray));
	}
	if (!binsPtCArray) AliFatal("Cannot determine Pt binning");
	if (fHistPtYTemplate->GetYaxis()->IsVariableBinSize()) binsYCArray = fHistPtYTemplate->GetYaxis()->GetXbins()->GetArray();
	else {
		/*binsYCArray = new Double_t[nBinsY+1];
		GetCustomBins(nBinsY, hLoY, hHiY, const_cast<Double_t* const>(binsYCArray));*/
		Double_t stdYbins[] = { -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6 }; // Temporary until GetCustomBins is debugged
		binsYCArray = stdYbins;
	}
	if (!binsYCArray) AliFatal("Cannot determine Y binning");
  Int_t nBinsMult = (fCollisionType == kPP)	? 50 : 1000;
	Int_t hMaxMult =  (fCollisionType == kPP)	? 500 : 10000;
	const Double_t binsCentCArray[] = { 0, 5, 10, 20, 30, 40, 50, 60, 80, 100 };
	Int_t nBinsCent = sizeof(binsCentCArray)/sizeof(Double_t)-1;
	//AliInfo(Form("%f %f %f %f %f %f", nBinsPt, hLoPt, hHiPt, nBinsY, hLoY, hHiY));
	
	Bool_t isESD = dynamic_cast<AliESDInputHandler*>(fMainInputHandler);
	
	// Rec (Common ESD/AOD)
	if (fIsAOD || isESD ) {
  fhRecMultUnbiased = new TH1F("fhRecMultUnbiased", "Reconstructed unbiased multiplicity (before PS); Number of tracks; Number of events", nBinsMult, 0, hMaxMult);
  fhRecMult = new TH1F("fhRecMult", "Reconstructed multiplicity (after PS); Number of tracks; Number of events", nBinsMult, 0, hMaxMult);
	if ( fCollisionType != kPP ) {
		fhRecMultCentralityUnbiased = new TH2F("fhRecMultCentralityUnbiased", "Reconstructed unbiased multiplicity (before PS) vs centrality; Number of tracks; Centrality", nBinsMult, 0, hMaxMult, nBinsCent, binsCentCArray);
		fhRecMultCentralityUnbiased->SetOption("colz");
	  fhRecMultCentrality= new TH2F("fhRecMultCentrality", "Reconstructed multiplicity (after PS) vs centrality; Number of tracks; Centrality", nBinsMult, 0, hMaxMult, nBinsCent, binsCentCArray);
		fhRecMultCentrality->SetOption("colz");
	}
	fhRecPrimaryVertexRPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhRecPrimaryVertexRPhiZ"));
	fhRecPrimaryVertexRPhiZ->GetAxis(0)->Set(200, 0, 20);
	fhRecPrimaryVertexRPhiZ->GetAxis(2)->Set(200, -50, 50);
	fhRecPrimaryVertexRPhiZ->SetTitle("Reconstructed primary vertex position;R (cm); phi (rad); Z (cm)");
  fhRecAllPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhRecAllPtY"));
  fhRecAllPtY->SetTitle("Transverse momentum vs Rapidity of all reconstructed tracks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");

	fOutputCont->Add(fhRecMultUnbiased);
	fOutputCont->Add(fhRecMult);
	fOutputCont->Add(fhRecMultCentralityUnbiased);
	fOutputCont->Add(fhRecMultCentrality);
	fOutputCont->Add(fhRecPrimaryVertexRPhiZ);
	fOutputCont->Add(fhRecAllPtY);
	
	/*}
	
	// ESD
	if (!fIsAOD && isESD) {*/
	fhESDAllKinkPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDAllKinkPtY"));
	fhESDAllKinkPtY->SetTitle("Transverse momentum vs Rapidity of all ESD Kinks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDUnknownKinkPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDUnknownKinkPtY"));
	fhESDUnknownKinkPtY->SetTitle("Transverse momentum vs Rapidity of all ESD Unknown Kinks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKaonKinkPtEta = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKaonKinkPtEta"));
	fhESDKaonKinkPtEta->SetTitle("Transverse momentum vs Pseudorapidity of all ESD K^{#pm} Kinks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKaonKinkPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKaonKinkPtY"));
	fhESDKaonKinkPtY->SetTitle("Transverse momentum vs Rapidity of all ESD K^{#pm} Kinks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPionKinkPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPionKinkPtY"));
	fhESDPionKinkPtY->SetTitle("Transverse momentum vs Rapidity of all ESD #pi^{#pm} Kinks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKPlusKinkPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKPlusKinkPtY"));
	fhESDKPlusKinkPtY->SetTitle("Transverse momentum vs Rapidity of all ESD K^{+} Kinks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKMinusKinkPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKMinusKinkPtY"));
	fhESDKMinusKinkPtY->SetTitle("Transverse momentum vs Rapidity of all ESD K^{-} Kinks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPiPlusKinkPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPiPlusKinkPtY"));
	fhESDPiPlusKinkPtY->SetTitle("Transverse momentum vs Rapidity of all ESD #pi^{+} Kinks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPiMinusKinkPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPiMinusKinkPtY"));
	fhESDPiMinusKinkPtY->SetTitle("Transverse momentum vs Rapidity of all ESD #pi^{-} Kinks; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKPlusKinkPt = fhESDKPlusKinkPtY->ProjectionX("fhESDKPlusKinkPt");
	fhESDKPlusKinkPt->SetTitle("Transverse momentum of of all ESD K^{+} Kinks; #it{p}_{T} (GeV/#it{c});d#it{N}/(d#it{p}_{T}) [(GeV/#it{c})]^{-2}");
	fhESDKMinusKinkPt = fhESDKMinusKinkPtY->ProjectionX("fhESDKMinusKinkPt");
	fhESDKMinusKinkPt->SetTitle("Transverse momentum of of all ESD K^{-} Kinks; #it{p}_{T} (GeV/#it{c});d#it{N}/(d#it{p}_{T}) [(GeV/#it{c})]^{-2}");
	fhESDAllKinkDecaysRPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhESDAllKinkDecaysRPhiZ"));
	fhESDAllKinkDecaysRPhiZ->SetTitle("fhESDAllKinkDecaysRPhiZ");
	fhESDKaonKinkDecaysRPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhESDKaonKinkDecaysRPhiZ"));
	fhESDKaonKinkDecaysRPhiZ->SetTitle("fhESDKaonKinkDecaysRPhiZ");
	fhESDPionKinkDecaysRPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhESDPionKinkDecaysRPhiZ"));
	fhESDPionKinkDecaysRPhiZ->SetTitle("fhESDPionKinkDecaysRPhiZ");
	if ( fCollisionType != kPP ) { // d^{3}N/(dp_{T}dydc) (c/GeV)
		fhESDKaonKinkPtYCent = new TH3F("fhESDKaonKinkPtYCent", "Transverse momentum vs Rapidity vs Centrality of all ESD K^{#pm} Kinks; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhESDPionKinkPtYCent = new TH3F("fhESDPionKinkPtYCent", "Transverse momentum vs Rapidity vs Centrality of all ESD #pi^{#pm} Kinks; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhESDKPlusKinkPtYCent = new TH3F("fhESDKPlusKinkPtYCent", "Transverse momentum vs Rapidity vs Centrality of all ESD K^{+} Kinks; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhESDPiPlusKinkPtYCent = new TH3F("fhESDPiPlusKinkPtYCent", "Transverse momentum vs Rapidity vs Centrality of all ESD #pi^{+} Kinks; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhESDKMinusKinkPtYCent = new TH3F("fhESDKMinusKinkPtYCent", "Transverse momentum vs Rapidity vs Centrality of all ESD K^{-} Kinks; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhESDPiMinusKinkPtYCent = new TH3F("fhESDPiMinusKinkPtYCent", "Transverse momentum vs Rapidity vs Centrality of all ESD #pi^{-} Kinks; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
	}
  fhESDnKinks = new TH1F("fhESDnKinks", "ESD kinks multiplicity (after PS); Number of kinks; Number of events", nBinsMult, 0, nBinsMult);
	fhESDKinkQt = new TH1F("fhESDKinkQt", "ESD kink q_{T} distribution; q_{T} (GeV/c); Number of events", nBinsPt, 0, 0.5);
	fhESDKinkAngle = new TH1F("fhESDKinkAngle", "ESD kink angle distribution; angle (deg); Number of events", nBinsPt, 0, 30);
	fhESDKaonKinkQt = new TH1F("fhESDKaonKinkQt", "ESD K^{#pm} kink q_{T} distribution; q_{T} (GeV/c); Number of events", nBinsPt, 0, 0.5);
	fhESDKaonKinkAngle = new TH1F("fhESDKaonKinkAngle", "ESD K^{#pm} kink angle distribution; angle (deg); Number of events", nBinsPt, 0, 30);
	fhESDPionKinkQt = new TH1F("fhESDPionKinkQt", "ESD #pi^{#pm} kink q_{T} distribution; q_{T} (GeV/c); Number of events", nBinsPt, 0, 0.5);
	fhESDPionKinkAngle = new TH1F("fhESDPionKinkAngle", "ESD #pi^{#pm} kink angle distribution; angle (deg); Number of events", nBinsPt, 0, 30);
	fhESDKinkPAngle = new TH2F("fhESDKinkPAngle", "ESD kink angle vs momentum distribution; p (GeV/c); angle (deg)", nBinsPt, hLoPt, hHiPt, nBinsPt, 0, 30);
	fhESDKinkPAngle->SetOption("colz");
	fhESDKinkDCA = new TH1F("fhESDKinkDCA", "ESD kink DCA distribution; DCA (cm); Number of events", nBinsPt, 0, 10);
	fhESDKaonOverPionEbE = new TH1F("fhESDKaonOverPionEbE", "ESD kink K/#pi distribution; q_{T} (GeV/c); Number of events", nBinsPt, 0, 0.05);
	fhESDMomentumTPCSignalKaonKinks = new TH2F("fhESDMomentumTPCSignalKaonKinks", "K^{#pm} Kinks momentum vs TPC de/dx; Mother momentum (GeV/c); TPC dE/dx central value", 10*nBinsPt, hLoPt, hHiPt+5, 150, 0, 300); // MSS compatible binning
	fhESDMomentumTPCSignalKaonKinks->SetOption("colz");
	fhESDKaonKinksMotherVSDaughterTPCncls = new TH2F("fhESDKaonKinksMotherVSDaughterTPCncls", "K^{#pm} Kinks mother vs daughter TPC nlcs; Mother ncls; Daughter ncls", 70, 0, 140, 70, 0, 140);
	fhESDKaonKinksMotherVSDaughterTPCncls->SetOption("colz");
	fhESDKaonKinksMotherAndDaughterTPCncls = new TH1F("fhESDKaonKinksMotherAndDaughterTPCncls", "K^{#pm} Kinks mother plus daughter TPC nlcs; Mother+Daughter ncls; Number of events", 100, 0, 200);
	fhESDMomentumTPCSignalPionKinks = new TH2F("fhESDMomentumTPCSignalPionKinks", "#pi^{#pm} Kinks momentum vs TPC de/dx; Mother momentum (GeV/c); TPC dE/dx central value", 10*nBinsPt, hLoPt, hHiPt+5, 150, 0, 300); // MSS compatible binning
	fhESDMomentumTPCSignalPionKinks->SetOption("colz");
	fhESDPionKinksMotherVSDaughterTPCncls = new TH2F("fhESDPionKinksMotherVSDaughterTPCncls", "#pi^{#pm} Kinks mother vs daughter TPC nlcs; Mother ncls; Daughter ncls", 70, 0, 140, 70, 0, 140);
	fhESDPionKinksMotherVSDaughterTPCncls->SetOption("colz");
	fhESDPionKinksMotherAndDaughterTPCncls = new TH1F("fhESDPionKinksMotherAndDaughterTPCncls", "#pi^{#pm} Kinks mother plus daughter TPC nlcs; Mother+Daughter ncls; Number of events", 100, 0, 200);
	fhESDInvMassKinkDecayMuNu = new TH1F("fhESDInvMassKinkDecayMuNu","Invariant mass of kink->mu+nu decays", 1000, 0.0, 1.0);
	
	fOutputCont->Add(fRecEventCuts->GetHist());
  fOutputCont->Add(fhESDnKinks);
	fOutputCont->Add(fhESDAllKinkDecaysRPhiZ);
	fOutputCont->Add(fhESDKinkQt);
	fOutputCont->Add(fhESDKinkAngle);
	fOutputCont->Add(fhESDKinkPAngle);
	fOutputCont->Add(fhESDKinkDCA);
	fOutputCont->Add(fhESDAllKinkPtY);
	fOutputCont->Add(fhESDUnknownKinkPtY);
	fOutputCont->Add(fhESDKaonKinkDecaysRPhiZ);
  fOutputCont->Add(fhESDKaonKinkPtEta);
  fOutputCont->Add(fhESDKaonKinkPtY);
  fOutputCont->Add(fhESDKPlusKinkPtY);
  fOutputCont->Add(fhESDKMinusKinkPtY);
  fOutputCont->Add(fhESDKPlusKinkPt);
  fOutputCont->Add(fhESDKMinusKinkPt);
	fOutputCont->Add(fhESDKaonKinkQt);
	fOutputCont->Add(fhESDKaonKinkAngle);
	fOutputCont->Add(fhESDPionKinkDecaysRPhiZ);
  fOutputCont->Add(fhESDPionKinkPtY);
  fOutputCont->Add(fhESDPiPlusKinkPtY);
  fOutputCont->Add(fhESDPiMinusKinkPtY);
	fOutputCont->Add(fhESDPionKinkQt);
	fOutputCont->Add(fhESDPionKinkAngle);
	fOutputCont->Add(fhESDMomentumTPCSignalKaonKinks);
	fOutputCont->Add(fhESDKaonKinksMotherVSDaughterTPCncls);
	fOutputCont->Add(fhESDKaonKinksMotherAndDaughterTPCncls);
	fOutputCont->Add(fhESDMomentumTPCSignalPionKinks);
	fOutputCont->Add(fhESDPionKinksMotherVSDaughterTPCncls);
	fOutputCont->Add(fhESDPionKinksMotherAndDaughterTPCncls);
	if ( fCollisionType != kPP ) {
		fOutputCont->Add(fhESDKaonKinkPtYCent);
		fOutputCont->Add(fhESDPionKinkPtYCent);
		fOutputCont->Add(fhESDKPlusKinkPtYCent);
		fOutputCont->Add(fhESDPiPlusKinkPtYCent);
		fOutputCont->Add(fhESDKMinusKinkPtYCent);
		fOutputCont->Add(fhESDPiMinusKinkPtYCent); 
	}
	fOutputCont->Add(fhESDKaonOverPionEbE);
	fOutputCont->Add(fRecKinkCutsKaon->fhCutFiducial);
	fOutputCont->Add(fRecKinkCutsKaon->fhCutFiducialKine);
	//fOutputCont->Add(fRecKinkCutsPion->fhCutFiducial);
	//fOutputCont->Add(fRecKinkCutsPion->fhCutFiducialKine);
	fOutputCont->Add(fhESDInvMassKinkDecayMuNu);
	}

	// MC
	if (fUseMC) {
 	fhMCmult = new TH1F("fhMCmult", "MC multiplicity; Number of tracks; Number of events", nBinsMult, 0, hMaxMult);  	
	fhMCmultPrim = new TH1F("fhMCmultPrim", "MC primary tracks multiplicity; Number of tracks; Number of events", nBinsMult, 0, hMaxMult);
  fhMCpdg = new TH1F("fhMCpdg", "fhMCpdg", 7000, -3500, 3500); 
  fhMCMotherDaughterPdg = new TH2F("fhMCMotherDaughterPdg", "fhMCMotherDaughterPdg", 7000, -3500, 3500, 7000, -3500, 3500);
	fhMCMotherDaughterPdg->SetDrawOption("text");
  fhMCkaonMotherPdg = new TH1F("fhMCkaonMotherPdg", "fhMCkaonMotherPdg", 7000, -3500, 3500); 
  fhMCpionMotherPdg = new TH1F("fhMCpionMotherPdg", "fhMCpionMotherPdg", 7000, -3500, 3500); 
	fhMCKaonKinkDecaysRPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhMCKaonKinkDecaysRPhiZ"));
	fhMCKaonKinkDecaysRPhiZ->SetTitle("fhMCKaonKinkDecaysRPhiZ");
	fhMCmainVertexRPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhMCmainVertexRPhiZ"));
	fhMCmainVertexRPhiZ->GetAxis(0)->Set(200, 0, 20);
	fhMCmainVertexRPhiZ->GetAxis(2)->Set(200, -50, 50);
	fhMCmainVertexRPhiZ->SetTitle("MC main vertex position;R (cm); phi (rad); Z (cm)");
	fhMCPtYall = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYall"));
	fhMCPtYall->SetTitle("Transverse momentum vs Rapidity of all generated particles; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYprim = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprim"));
	fhMCPtYprim->SetTitle("Transverse momentum vs Rapidity of primary generated particles; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYprimKaon = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprimKaon"));
	fhMCPtYprimKaon->SetTitle("Transverse momentum vs Rapidity of primary generated K^{#pm}; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYkinkKaonFiducial = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYkinkKaonFiducial"));
	fhMCPtYkinkKaonFiducial->SetTitle("Transverse momentum vs Rapidity of primary generated K^{#pm} Kink decays in fiducial volume; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYprimKPlus = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprimKPlus"));
	fhMCPtYprimKPlus->SetTitle("Transverse momentum vs Rapidity of primary generated K^{+}; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYkinkKPlusFiducial = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYkinkKPlusFiducial"));
	fhMCPtYkinkKPlusFiducial->SetTitle("Transverse momentum vs Rapidity of primary generated K^{+} Kink decays in fiducial volume; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYprimKMinus = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprimKMinus"));
	fhMCPtYprimKMinus->SetTitle("Transverse momentum vs Rapidity of primary generated K^{-}; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYkinkKMinusFiducial = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYkinkKMinusFiducial"));
	fhMCPtYkinkKMinusFiducial->SetTitle("Transverse momentum vs Rapidity of primary generated K^{-} Kink decays in fiducial volume; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtprimKPlus = fhMCPtYprimKPlus->ProjectionX("fhMCPtprimKPlus");
	fhMCPtprimKPlus->SetTitle("Transverse momentum of primary generated K^{+}; #it{p}_{T} (GeV/#it{c});d#it{N}/(d#it{p}_{T}) [(GeV/#it{c})]^{-2}");
	fhMCPtprimKMinus = fhMCPtYprimKMinus->ProjectionX("fhMCPtprimKMinus");
	fhMCPtprimKMinus->SetTitle("Transverse momentum of primary generated K^{-}; #it{p}_{T} (GeV/#it{c});d#it{N}/(d#it{p}_{T}) [(GeV/#it{c})]^{-2}");
	fhMCPhiPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPhiPtY"));
	fhMCPhiPtY->SetTitle("fhMCPhiPtY; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYprimPhi = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprimPhi"));
	fhMCPtYprimPhi->SetTitle("fhMCPtYprimPhi; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCInvMassPtprimPhi = new TH2F("fhMCInvMassPtprimPhi", "fhMCInvMassPtprimPhi; #it{M}_{KK} (GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c})", 200, 0.95, 1.15, nBinsPtInvMassSlice, binsPtCArray);
	fhMCInvMassPtprimPhi->SetOption("colz");
	fhMCPhiInvMassPtYtest = new TH3F("fhMCPhiInvMassPtYtest", "fhMCPhiInvMassPtYtest; #it{M}_{KK} (GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}", fhMCInvMassPtprimPhi->GetNbinsX(), fhMCInvMassPtprimPhi->GetXaxis()->GetXbins()->GetArray(),  
	nBinsPtInvMassSlice, binsPtCArray, nBinsY, binsYCArray);
	fhMCKstarPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCKstarPtY"));
	fhMCKstarPtY->SetTitle("fhMCKstarPtY; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYprimKstar = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprimKstar"));
	fhMCPtYprimKstar->SetTitle("fhMCPtYprimKstar; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCInvMassPtprimKstar = new TH2F("fhMCInvMassPtprimKstar", "fhMCInvMassPtprimKstar; #it{M}_{K#pi} (GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c})", 50, 0.6, 1.1, nBinsPtInvMassSlice, binsPtCArray);
	fhMCInvMassPtprimKstar->SetOption("colz");
	fhMCLambdaPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCLambdaPtY"));
	fhMCLambdaPtY->SetTitle("fhMCLambdaPtY; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYprimLambda = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprimLambda"));
	fhMCPtYprimLambda->SetTitle("fhMCPtYprimLambda; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCInvMassPtprimLambda = new TH2F("fhMCInvMassPtprimLambda", "fhMCInvMassPtprimLambda; #it{M}_{K#pi} (GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c})", 100, 1.4, 1.8, nBinsPtInvMassSlice, binsPtCArray);
	fhMCInvMassPtprimKstar->SetOption("colz");
	fhMCPhi2KaonPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPhi2KaonPtY"));
	fhMCPhi2KaonPtY->SetTitle("Transverse momentum vs Rapidity of primary generated Phi decaying to Kaons; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCKstar2KaonPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCKstar2KaonPtY"));
	fhMCKstar2KaonPtY->SetTitle("Transverse momentum vs Rapidity of primary generated Kstar decaying to Kaons; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCKaon2MuonPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCKaon2MuonPtY"));
	fhMCKaon2MuonPtY->SetTitle("Transverse momentum vs Rapidity of primary generated Kaons decaying to Muons; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCKaon2PionPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCKaon2PionPtY"));
	fhMCKaon2PionPtY->SetTitle("Transverse momentum vs Rapidity of primary generated Kaons decaying to Pions; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYprimPion = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprimPion"));
	fhMCPtYprimPion->SetTitle("Transverse momentum vs Rapidity of primary generated #pi^{#pm}; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPionKinkDecaysRPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhMCPionKinkDecaysRPhiZ"));
	fhMCPionKinkDecaysRPhiZ->SetTitle("fhMCPionKinkDecaysRPhiZ");
	fhMCPion2MuonPtY = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPion2MuonPtY"));
	fhMCPion2MuonPtY->SetTitle("fhMCPion2MuonPtY");
	fhMCPtYkinkPionFiducial = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYkinkPionFiducial"));
	fhMCPtYkinkPionFiducial->SetTitle("Transverse momentum vs Rapidity of primary generated #pi^{#pm} Kink decays in fiducial volume; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCKinkKaonQt = new TH1F("fhMCKinkKaonQt", "MC kink q_{T} distribution; q_{T} (GeV/c); Number of events", nBinsPt, 0, 0.5);
	fhMCKinkKaonAngle = new TH1F("fhMCKinkKaonAngle", "MC kink angle distribution; angle (deg); Number of events", nBinsPt, 0, 30);
	fhMCKaonKinkDecaysTrackLength = new TH1F("fhMCKaonKinkDecaysTrackLength", "fhMCKaonKinkDecaysTrackLength; Track length (cm); Number of events", 250, 0, 500);
	fhMCKaonKinkDecaysTrackTime = new TH1F("fhMCKaonKinkDecaysTrackTime", "fhMCKaonKinkDecaysTrackTime; Track time (ns); Number of events", 1000, 0, 1);
	
	fhMCPtYprimPiPlus = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprimPiPlus"));
	fhMCPtYprimPiPlus->SetTitle("Transverse momentum vs Rapidity of primary generated #pi^{+}; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYkinkPiPlusFiducial = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYkinkPiPlusFiducial"));
	fhMCPtYkinkPiPlusFiducial->SetTitle("Transverse momentum vs Rapidity of primary generated #pi^{+} Kink decays in fiducial volume; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYprimPiMinus = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYprimPiMinus"));
	fhMCPtYprimPiMinus->SetTitle("Transverse momentum vs Rapidity of primary generated #pi^{-}; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCPtYkinkPiMinusFiducial = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhMCPtYkinkPiMinusFiducial"));
	fhMCPtYkinkPiMinusFiducial->SetTitle("Transverse momentum vs Rapidity of primary generated #pi^{-} Kink decays in fiducial volume; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhMCKinkPionQt = new TH1F("fhMCKinkPionQt", "MC Pion kink q_{T} distribution; q_{T} (GeV/c); Number of events", nBinsPt, 0, 0.5);
	fhMCKinkPionAngle = new TH1F("fhMCKinkPionAngle", "MC Pion kink angle distribution; angle (deg); Number of events", nBinsPt, 0, 30);
	if ( fCollisionType != kPP ) { // d^{3}N/(dp_{T}dydc) (c/GeV)
		fhMCPtYCentprimKaon = new TH3F("fhMCPtYCentprimKaon", "Transverse momentum vs Rapidity vs Centrality of primary generated K^{#pm}; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhMCPtYCentprimKPlus = new TH3F("fhMCPtYCentprimKPlus", "Transverse momentum vs Rapidity vs Centrality of primary generated K^{+}; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhMCPtYCentprimKMinus = new TH3F("fhMCPtYCentprimKMinus", "Transverse momentum vs Rapidity vs Centrality of primary generated K^{-}; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhMCPtYCentprimPion = new TH3F("fhMCPtYCentprimPion", "Transverse momentum vs Rapidity vs Centrality of primary generated #pi^{#pm}; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhMCPtYCentprimPiPlus = new TH3F("fhMCPtYCentprimPiPlus", "Transverse momentum vs Rapidity vs Centrality of primary generated #pi^{+}; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		fhMCPtYCentprimPiMinus = new TH3F("fhMCPtYCentprimPiMinus", "Transverse momentum vs Rapidity vs Centrality of primary generated #pi^{-}; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		//TODO://fhMCInvMassPtCentprimPhi = new TH3F("fhMCInvMassPtCentprimPhi", "Invariant Mass vs Transverse momentum vs Centrality of primary generated #Phi; #it{p}_{T} (GeV/#it{c}); y; centrality", nBinsPt, binsPtCArray, nBinsY, binsYCArray, nBinsCent, binsCentCArray);
		//fhMCInvMassPtCentprimKstar
		//fhMCInvMassPtCentprimLambda
	}
	fhMCPhiDecayOpeningAngle = new TH2F("fhMCPhiDecayOpeningAngle", "Transverse momentum vs cos of decay opening angle of primary generated Phi decaying to Kaons; #it{p}_{T} (GeV/#it{c}); cos(#theta); d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}", nBinsPt, binsPtCArray, 20, -1, 1);
	fhMCPhiDecayOpeningAngle->SetOption("colz");
	fhMCKstarDecayOpeningAngle = new TH2F("fhMCKstarDecayOpeningAngle", "Transverse momentum vs cos of decay opening angle of primary generated Kstar decaying to Kaons; #it{p}_{T} (GeV/#it{c}); cos(#theta); d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}", nBinsPt, binsPtCArray, 20, -1, 1); 
	fhMCKstarDecayOpeningAngle->SetOption("colz");
	fhMCLambdaDecayOpeningAngle = new TH2F("fhMCLambdaDecayOpeningAngle", "Transverse momentum vs cos of decay opening angle of primary generated Lambda decaying to Kaons; #it{p}_{T} (GeV/#it{c}); cos(#theta); d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}", nBinsPt, binsPtCArray, 20, -1, 1); 
	fhMCLambdaDecayOpeningAngle->SetOption("colz");
	
	fhMCKaonOverPionEbE = new TH1F("fhMCKaonOverPionEbE", "MC kink K/#pi distribution; q_{T} (GeV/c); Number of events", nBinsPt, 0, 0.05);

	fOutputCont->Add(fhMCmainVertexRPhiZ);
  fOutputCont->Add(fhMCmult);
	fOutputCont->Add(fhMCmultPrim);
	fOutputCont->Add(fhMCPtYall);
	fOutputCont->Add(fhMCPtYprim);
  fOutputCont->Add(fhMCpdg);
	fOutputCont->Add(fhMCPhiPtY);
	fOutputCont->Add(fhMCPtYprimPhi);
	fOutputCont->Add(fhMCPhiInvMassPtYtest);
	fOutputCont->Add(fhMCInvMassPtprimPhi);
	fOutputCont->Add(fhMCKstarPtY);
	fOutputCont->Add(fhMCPtYprimKstar);
	fOutputCont->Add(fhMCInvMassPtprimKstar);
	fOutputCont->Add(fhMCLambdaPtY);
	fOutputCont->Add(fhMCPtYprimLambda);
	fOutputCont->Add(fhMCMotherDaughterPdg);
	fOutputCont->Add(fhMCkaonMotherPdg);
	fOutputCont->Add(fhMCpionMotherPdg);
  fOutputCont->Add(fhMCPtYprimKaon);
  fOutputCont->Add(fhMCPtYprimKPlus);
  fOutputCont->Add(fhMCPtYprimKMinus);
  fOutputCont->Add(fhMCPtprimKPlus);
  fOutputCont->Add(fhMCPtprimKMinus);
	fOutputCont->Add(fhMCPtYkinkKaonFiducial);
	fOutputCont->Add(fhMCPtYkinkKPlusFiducial);
	fOutputCont->Add(fhMCPtYkinkKMinusFiducial);
	fOutputCont->Add(fhMCKaonKinkDecaysTrackLength);
	fOutputCont->Add(fhMCKaonKinkDecaysTrackTime);
	fOutputCont->Add(fhMCKinkKaonQt);
	fOutputCont->Add(fhMCKinkKaonAngle);
	fOutputCont->Add(fhMCPhi2KaonPtY);
	fOutputCont->Add(fhMCKstar2KaonPtY);
	fOutputCont->Add(fhMCKaon2MuonPtY);
	fOutputCont->Add(fhMCKaon2PionPtY);
	fOutputCont->Add(fhMCKaonKinkDecaysRPhiZ);
	fOutputCont->Add(fhMCPtYprimPion);
  fOutputCont->Add(fhMCPtYprimPiPlus);
  fOutputCont->Add(fhMCPtYprimPiMinus);
	fOutputCont->Add(fhMCPtYkinkPionFiducial);
	fOutputCont->Add(fhMCPtYkinkPiPlusFiducial);
	fOutputCont->Add(fhMCPtYkinkPiMinusFiducial);
	//fOutputCont->Add(fhMCPionKinkDecaysTrackLength);
	//fOutputCont->Add(fhMCPionKinkDecaysTrackTime);
	fOutputCont->Add(fhMCKinkPionQt);
	fOutputCont->Add(fhMCKinkPionAngle);
	fOutputCont->Add(fhMCPionKinkDecaysRPhiZ);
	fOutputCont->Add(fhMCPion2MuonPtY);
	fOutputCont->Add(fhMCKaonOverPionEbE);
	fOutputCont->Add(fhMCPhiDecayOpeningAngle);
	fOutputCont->Add(fhMCKstarDecayOpeningAngle);
	fOutputCont->Add(fhMCLambdaDecayOpeningAngle);
	if ( fCollisionType != kPP ) {
		fOutputCont->Add(fhMCPtYCentprimKaon);
		fOutputCont->Add(fhMCPtYCentprimKPlus);
		fOutputCont->Add(fhMCPtYCentprimKMinus);
		fOutputCont->Add(fhMCPtYCentprimPion);
		fOutputCont->Add(fhMCPtYCentprimPiPlus);
		fOutputCont->Add(fhMCPtYCentprimPiMinus);
		//fOutputCont->Add(fhMCInvMassPtCentprimPhi);
		//fOutputCont->Add(fhMCInvMassPtCentprimKstar);
		//fOutputCont->Add(fhMCInvMassPtCentprimLambda);
	}
	//fOutputCont->Add(fRecKinkCutsKaon->fhCutFiducial);
	//fOutputCont->Add(fRecKinkCutsKaon->fhCutFiducialKine);
	
	// MC truth check
	if (fIsAOD || isESD) {
	fhESDKinkQtTrueMC = new TH1F("fhESDKinkQtTrueMC", "MC truth kink q_{T} distribution; q_{T} (GeV/c); Number of events", nBinsPt, 0, 0.5);
	fhESDKinkAngleTrueMC = new TH1F("fhESDKinkAngleTrueMC", "MC truth kink angle distribution; angle (deg); Number of events", nBinsPt, 0, 30);
	fhESDKaonKinkPtYTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKaonKinkPtYTrueMC"));
	fhESDKaonKinkPtYTrueMC->SetTitle("fhESDKaonKinkPtYTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKaonKinkPtYFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKaonKinkPtYFakeMC"));
	fhESDKaonKinkPtYFakeMC->SetTitle("fhESDKaonKinkPtYFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKaonKinkPtYTrueMCsecondary = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKaonKinkPtYTrueMCsecondary"));
	fhESDKaonKinkPtYTrueMCsecondary->SetTitle("fhESDKaonKinkPtYTrueMCsecondary; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKPlusKinkPtYTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKPlusKinkPtYTrueMC"));
	fhESDKPlusKinkPtYTrueMC->SetTitle("fhESDKPlusKinkPtYTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKPlusKinkPtYFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKPlusKinkPtYFakeMC"));
	fhESDKPlusKinkPtYFakeMC->SetTitle("fhESDKPlusKinkPtYFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKPlusKinkPtYTrueMCsecondary = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKPlusKinkPtYTrueMCsecondary"));
	fhESDKPlusKinkPtYTrueMCsecondary->SetTitle("fhESDKPlusKinkPtYTrueMCsecondary; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKMinusKinkPtYTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKMinusKinkPtYTrueMC"));
	fhESDKMinusKinkPtYTrueMC->SetTitle("fhESDKMinusKinkPtYTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKMinusKinkPtYFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKMinusKinkPtYFakeMC"));
	fhESDKMinusKinkPtYFakeMC->SetTitle("fhESDKMinusKinkPtYFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKMinusKinkPtYTrueMCsecondary = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKMinusKinkPtYTrueMCsecondary"));
	fhESDKMinusKinkPtYTrueMCsecondary->SetTitle("fhESDKMinusKinkPtYTrueMCsecondary; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPionKinkPtYTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPionKinkPtYTrueMC"));
	fhESDPionKinkPtYTrueMC->SetTitle("fhESDPionKinkPtYTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPionKinkPtYFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPionKinkPtYFakeMC"));
	fhESDPionKinkPtYFakeMC->SetTitle("fhESDPionKinkPtYFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPionKinkPtYTrueMCsecondary = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPionKinkPtYTrueMCsecondary"));
	fhESDPionKinkPtYTrueMCsecondary->SetTitle("fhESDPionKinkPtYTrueMCsecondary; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPiPlusKinkPtYTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPiPlusKinkPtYTrueMC"));
	fhESDPiPlusKinkPtYTrueMC->SetTitle("fhESDPiPlusKinkPtYTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPiPlusKinkPtYFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPiPlusKinkPtYFakeMC"));
	fhESDPiPlusKinkPtYFakeMC->SetTitle("fhESDPiPlusKinkPtYFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPiPlusKinkPtYTrueMCsecondary = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPiPlusKinkPtYTrueMCsecondary"));
	fhESDPiPlusKinkPtYTrueMCsecondary->SetTitle("fhESDPiPlusKinkPtYTrueMCsecondary; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPiMinusKinkPtYTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPiMinusKinkPtYTrueMC"));
	fhESDPiMinusKinkPtYTrueMC->SetTitle("fhESDPiMinusKinkPtYTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPiMinusKinkPtYFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPiMinusKinkPtYFakeMC"));
	fhESDPiMinusKinkPtYFakeMC->SetTitle("fhESDPiMinusKinkPtYFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPiMinusKinkPtYTrueMCsecondary = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPiMinusKinkPtYTrueMCsecondary"));
	fhESDPiMinusKinkPtYTrueMCsecondary->SetTitle("fhESDPiMinusKinkPtYTrueMCsecondary; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPhi2KaonPtYTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPhi2KaonPtYTrueMC"));
	fhESDPhi2KaonPtYTrueMC->SetTitle("fhESDPhi2KaonPtYTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKstar2KaonPtYTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKstar2KaonPtYTrueMC"));
	fhESDKstar2KaonPtYTrueMC->SetTitle("fhESDKstar2KaonPtYTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDLambda2KaonPtYTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDLambda2KaonPtYTrueMC"));
	fhESDLambda2KaonPtYTrueMC->SetTitle("fhESDLambda2KaonPtYTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");

	fhESDPhiPtYtestTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPhiPtYtestTrueMC"));
	fhESDPhiPtYtestTrueMC->SetTitle("fhESDPhiPtYtestTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPhiPtYtestFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPhiPtYtestFakeMC"));
	fhESDPhiPtYtestFakeMC->SetTitle("fhESDPhiPtYtestFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPhiPtYtestLikeSignTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPhiPtYtestLikeSignTrueMC"));
	fhESDPhiPtYtestLikeSignTrueMC->SetTitle("fhESDPhiPtYtestLikeSignTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDPhiPtYtestLikeSignFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDPhiPtYtestLikeSignFakeMC"));
	fhESDPhiPtYtestLikeSignFakeMC->SetTitle("fhESDPhiPtYtestLikeSignFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKstarPtYtestTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKstarPtYtestTrueMC"));
	fhESDKstarPtYtestTrueMC->SetTitle("fhESDKstarPtYtestTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKstarPtYtestFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKstarPtYtestFakeMC"));
	fhESDKstarPtYtestFakeMC->SetTitle("fhESDKstarPtYtestFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKstarPtYtestLikeSignTrueMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKstarPtYtestLikeSignTrueMC"));
	fhESDKstarPtYtestLikeSignTrueMC->SetTitle("fhESDKstarPtYtestLikeSignTrueMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhESDKstarPtYtestLikeSignFakeMC = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhESDKstarPtYtestLikeSignFakeMC"));
	fhESDKstarPtYtestLikeSignFakeMC->SetTitle("fhESDKstarPtYtestLikeSignFakeMC; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");

	fOutputCont->Add(fhESDKinkQtTrueMC);
	fOutputCont->Add(fhESDKinkAngleTrueMC);
	fOutputCont->Add(fhESDKaonKinkPtYTrueMC);
	fOutputCont->Add(fhESDKaonKinkPtYFakeMC);
	fOutputCont->Add(fhESDKaonKinkPtYTrueMCsecondary);
	fOutputCont->Add(fhESDKPlusKinkPtYTrueMC);
	fOutputCont->Add(fhESDKPlusKinkPtYFakeMC);
	fOutputCont->Add(fhESDKPlusKinkPtYTrueMCsecondary);
	fOutputCont->Add(fhESDKMinusKinkPtYTrueMC);
	fOutputCont->Add(fhESDKMinusKinkPtYFakeMC);
	fOutputCont->Add(fhESDKMinusKinkPtYTrueMCsecondary);
	fOutputCont->Add(fhESDPionKinkPtYTrueMC);
	fOutputCont->Add(fhESDPionKinkPtYFakeMC);
	fOutputCont->Add(fhESDPionKinkPtYTrueMCsecondary);
	fOutputCont->Add(fhESDPiPlusKinkPtYTrueMC);
	fOutputCont->Add(fhESDPiPlusKinkPtYFakeMC);
	fOutputCont->Add(fhESDPiPlusKinkPtYTrueMCsecondary);
	fOutputCont->Add(fhESDPiMinusKinkPtYTrueMC);
	fOutputCont->Add(fhESDPiMinusKinkPtYFakeMC);
	fOutputCont->Add(fhESDPiMinusKinkPtYTrueMCsecondary);
	fOutputCont->Add(fhESDPhi2KaonPtYTrueMC);
	fOutputCont->Add(fhESDKstar2KaonPtYTrueMC);
	fOutputCont->Add(fhESDLambda2KaonPtYTrueMC);
	fOutputCont->Add(fhESDPhiPtYtestTrueMC);
	fOutputCont->Add(fhESDPhiPtYtestFakeMC);
	fOutputCont->Add(fhESDPhiPtYtestLikeSignTrueMC);
	fOutputCont->Add(fhESDPhiPtYtestLikeSignFakeMC);
	fOutputCont->Add(fhESDKstarPtYtestTrueMC);
	fOutputCont->Add(fhESDKstarPtYtestFakeMC);
	fOutputCont->Add(fhESDKstarPtYtestLikeSignTrueMC);
	fOutputCont->Add(fhESDKstarPtYtestLikeSignFakeMC);
	}
	}
	
	// AOD test
	if (fIsAOD) {
	fhAODnVertices = new TH1F("fhAODnVertices", "AOD vertex multiplicity; Number of vertices; Number of events", nBinsMult, 0, hMaxMult);
	//fhAODnKinkVertices = new TH1F("fhAODnKinkVertices", "AOD kink vertex multiplicity; Number of kink vertices; Number of events", nBinsMult, 0, hMaxMult/10);
	fhAODnKinkDaughters = new TH1F("fhAODnKinkDaughters", "AOD kink daughters multiplicity; Number of kink daughters; Number of events", nBinsMult, 0, hMaxMult);
	fhAODAllKinkDecaysRPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhAODAllKinkDecaysRPhiZ"));
	fhAODAllKinkDecaysRPhiZ->SetTitle("fhAODAllKinkDecaysRPhiZ");
	fhAODKinkMothervsDaughterPIDcheck = new TH2F("fhAODKinkMothervsDaughterPIDcheck", "fhAODKinkMothervsDaughterPIDcheck", 14, -2, 12, 14, -2, 12);
	TAxis* axis = fhAODKinkMothervsDaughterPIDcheck->GetXaxis();
	axis->SetBinLabel(AliAODTrack::kElectron, "Electron");
	axis->SetBinLabel(AliAODTrack::kMuon, "Muon");
	axis->SetBinLabel(AliAODTrack::kPion, "Pion");
	axis->SetBinLabel(AliAODTrack::kKaon, "Kaon");
	axis->SetBinLabel(AliAODTrack::kProton, "Proton");
	TH2* hAODPtYTemplate = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("hAODPtYTemplate"));
	fhAODKaonKinkPID = hAODPtYTemplate->ProjectionX("fhAODKaonKinkPID");
	fhAODKaonKinkPID->SetTitle("fhAODKaonKinkPID; #it{p}_{T} (GeV/#it{c}); dN/dp_{T} (c/GeV)");
	fhAODKaonKinkPIDFake = hAODPtYTemplate->ProjectionX("fhAODKaonKinkPIDFake");
	fhAODKaonKinkPIDFake->SetTitle("fhAODKaonKinkPIDFake; #it{p}_{T} (GeV/#it{c}); dN/dp_{T} (c/GeV)");
	fhAODmotherVertexRPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhAODmotherVertexRPhiZ"));
	fhAODmotherVertexRPhiZ->SetTitle("fhAODmotherVertexRPhiZ");
	fhAODmotherVertexDCARPhiZ = dynamic_cast<THnSparseF*>(fHistRPhiZTemplate->Clone("fhAODmotherVertexDCARPhiZ"));
	fhAODmotherVertexDCARPhiZ->SetTitle("fhAODmotherVertexDCARPhiZ");
	//fhAODKinkKaonQt = new TH1F("fhAODKinkKaonQt", "AOD Kaon kink q_{T} distribution; q_{T} (GeV/c); Number of events", nBinsPt, 0, 0.5);
	//fhAODKinkKaonAngle = new TH1F("fhAODKinkKaonAngle", "AOD Kaon kink angle distribution; angle (deg); Number of events", nBinsPt, 0, 30);
	
	fOutputCont->Add(fhAODnVertices);
	//fOutputCont->Add(fhAODnKinkVertices); //fhESDnKinks
	fOutputCont->Add(fhAODnKinkDaughters);
	//fOutputCont->Add(fhAODAllKinkDecaysRPhiZ); //fhESDAllKinkDecaysRPhiZ
	fOutputCont->Add(fhAODKinkMothervsDaughterPIDcheck);
	fOutputCont->Add(fhAODKaonKinkPID);
	fOutputCont->Add(fhAODKaonKinkPIDFake);
	fOutputCont->Add(fhAODmotherVertexRPhiZ);	
	fOutputCont->Add(fhAODmotherVertexDCARPhiZ);
	//fOutputCont->Add(fhAODKinkKaonQt); //fhESDKinkQt
	//fOutputCont->Add(fhAODKinkKaonAngle); //fhESDKinkAngle
	}

	// Rec Resonances
	if (fEsdTrackCutsPartner) {
	fhRecPhiInvMassPt = new TH2F("fhRecPhiInvMassPt", "fhRecPhiInvMassPt; #it{M}_{KK} (GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c})", 200, 0.95, 1.15, nBinsPtInvMassSlice, binsPtCArray);
	fhRecPhiInvMassPt->SetOption("colz");
	fhRecPhiInvMassPtLikeSign = dynamic_cast<TH2*>(fhRecPhiInvMassPt->Clone("fhRecPhiInvMassPtLikeSign"));
	fhRecPhiInvMassPtLikeSign->SetTitle("fhRecPhiInvMassPtLikeSign");
	fhRecPhiInvMassPtYtest = new TH3F("fhRecPhiInvMassPtYtest", "fhRecPhiInvMassPtYtest; #it{M}_{KK} (GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}", fhRecPhiInvMassPt->GetNbinsX(), fhRecPhiInvMassPt->GetXaxis()->GetXbins()->GetArray(),  
	nBinsPtInvMassSlice, binsPtCArray, nBinsY, binsYCArray);
	fhRecPhiInvMassPtYtestLikeSign = dynamic_cast<TH3*>(fhRecPhiInvMassPtYtest->Clone("fhRecPhiInvMassPtYtestLikeSign"));
	fhRecPhiInvMassPtYtestLikeSign->SetTitle("fhRecPhiInvMassPtYtestLikeSign");
	fhRecPhiInvMassBothKinks = fhRecPhiInvMassPt->ProjectionX("fhRecPhiInvMassBothKinks");//new TH1F("fhRecPhiInvMassBothKinks", "fhRecPhiInvMassBothKinks; #it{M}_{KK} (GeV/#it{c^2})", fhRecPhiInvMassPt->GetNbinsX(), fhRecPhiInvMassPt->GetXaxis()->GetXbins()->GetArray());
	fhRecPhiInvMassBothKinks->SetTitle("fhRecPhiInvMassBothKinks; #it{M}_{KK} (GeV/#it{c^2})");
	//fhRecPhiInvMassBothKinks->Reset();
	fhRecPhiInvMassBothKinksLikeSign = dynamic_cast<TH1*>(fhRecPhiInvMassBothKinks->Clone("fhRecPhiInvMassBothKinksLikeSign"));//new TH1F("fhRecPhiInvMassBothKinksLikeSign", "fhRecPhiInvMassBothKinksLikeSign; #it{M}_{KK} (GeV/#it{c^2})", fhRecPhiInvMassPt->GetNbinsX(), fhRecPhiInvMassPt->GetXaxis()->GetXbins()->GetArray());
	fhRecPhiInvMassBothKinksLikeSign->SetTitle("fhRecPhiInvMassBothKinksLikeSign; #it{M}_{KK} (GeV/#it{c^2})");
	fhRecKstarInvMassPt = new TH2F("fhRecKstarInvMassPt", "fhRecKstarInvMassPt; #it{M}_{K#pi} (GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c})", 50, 0.6, 1.1, nBinsPtInvMassSlice, binsPtCArray);
	fhRecKstarInvMassPt->SetOption("colz");
	fhRecKstarInvMassPtLikeSign = dynamic_cast<TH2*>(fhRecKstarInvMassPt->Clone("fhRecKstarInvMassPtLikeSign"));
	fhRecKstarInvMassPtLikeSign->SetTitle("fhRecKstarInvMassPtLikeSign");
	fhRecKstarInvMassBothKinks = fhRecKstarInvMassPt->ProjectionX("fhRecKstarInvMassBothKinks");//new TH1F("fhRecKstarInvMassBothKinks", "fhRecKstarInvMassBothKinks; #it{M}_{K#pi} (GeV/#it{c^2})", fhRecKstarInvMassPt->GetNbinsX(), fhRecKstarInvMassPt->GetXaxis()->GetXbins()->GetArray());
	fhRecKstarInvMassBothKinks->SetTitle("fhRecKstarInvMassBothKinks; #it{M}_{K#pi} (GeV/#it{c^2})");
	//fhRecKstarInvMassBothKinks->Reset();
	fhRecKstarInvMassBothKinksLikeSign = dynamic_cast<TH1*>(fhRecKstarInvMassBothKinks->Clone("fhRecKstarInvMassBothKinksLikeSign"));//new TH1F("fhRecKstarInvMassBothKinksLikeSign", "fhRecKstarInvMassBothKinksLikeSign; #it{M}_{K#pi} (GeV/#it{c^2})", fhRecKstarInvMassPt->GetNbinsX(), fhRecKstarInvMassPt->GetXaxis()->GetXbins()->GetArray());
	fhRecKstarInvMassBothKinksLikeSign->SetTitle("fhRecKstarInvMassBothKinksLikeSign; #it{M}_{K#pi} (GeV/#it{c^2})");
	fhRecLambdaPtYtest = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhRecLambdaPtYtest"));
	fhRecLambdaPtYtest->SetTitle("fhRecLambdaPtYtest; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhRecLambdaPtYtestLikeSign = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhRecLambdaPtYtestLikeSign"));
	fhRecLambdaPtYtestLikeSign->SetTitle("fhRecLambdaPtYtestLikeSign; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");	
	fhRecLambdaInvMassPt = new TH2F("fhRecLambdaInvMassPt", "fhESDInvariantMassLambda; #it{M}_{Kp} (GeV/#it{c^2})", 100, 1.4, 1.8, nBinsPtInvMassSlice, binsPtCArray);
	fhRecLambdaInvMassPt->SetOption("colz");
	fhRecLambdaInvMassPtLikeSign = new TH2F("fhRecLambdaInvMassPtLikeSign", "fhESDInvariantMassLambdaLikeSign; #it{M}_{Kp} (GeV/#it{c^2})", 100, 1.4, 1.8, nBinsPtInvMassSlice, binsPtCArray);
	fhRecLambdaInvMassPtLikeSign->SetOption("colz");
	fhRecPhiPtYtest = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhRecPhiPtYtest"));
	fhRecPhiPtYtest->SetTitle("fhRecPhiPtYtest; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhRecPhiPtYtestLikeSign = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhRecPhiPtYtestLikeSign"));
	fhRecPhiPtYtestLikeSign->SetTitle("fhRecPhiPtYtestLikeSign; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhRecKstarPtYtest = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhRecKstarPtYtest"));
	fhRecKstarPtYtest->SetTitle("fhRecKstarPtYtest; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");
	fhRecKstarPtYtestLikeSign = dynamic_cast<TH2*>(fHistPtYTemplate->Clone("fhRecKstarPtYtestLikeSign"));
	fhRecKstarPtYtestLikeSign->SetTitle("fhRecKstarPtYtestLikeSign; #it{p}_{T} (GeV/#it{c}); y; d^2#it{N}/(d#it{p}_{T}dy) [(GeV/#it{c})]^{-2}");

  fOutputCont->Add(fhRecPhiInvMassPt);
  fOutputCont->Add(fhRecPhiInvMassPtLikeSign);
  fOutputCont->Add(fhRecPhiInvMassBothKinks);
	fOutputCont->Add(fhRecPhiInvMassBothKinksLikeSign);
	fOutputCont->Add(fhRecPhiPtYtest);
	fOutputCont->Add(fhRecPhiPtYtestLikeSign);
	fOutputCont->Add(fhRecPhiInvMassPtYtest);
	fOutputCont->Add(fhRecPhiInvMassPtYtestLikeSign);
  fOutputCont->Add(fhRecKstarInvMassPt);
  fOutputCont->Add(fhRecKstarInvMassPtLikeSign);
  fOutputCont->Add(fhRecKstarInvMassBothKinks);
	fOutputCont->Add(fhRecKstarInvMassBothKinksLikeSign);
	fOutputCont->Add(fhRecKstarPtYtest);
	fOutputCont->Add(fhRecKstarPtYtestLikeSign);
	fOutputCont->Add(fhRecLambdaInvMassPt);
	fOutputCont->Add(fhRecLambdaInvMassPtLikeSign);
	fOutputCont->Add(fhRecLambdaPtYtest);
	fOutputCont->Add(fhRecLambdaPtYtestLikeSign);
  }
  // Post output data.
  PostData(1, fOutputCont);
}

//________________________________________________________________________
void AliAnalysisTaskKinksFilimon::UserExec(Option_t *) 
{
	if (!fMultiInputHandler) KinkAnalysis();
  // Post output data.
  PostData(1, fOutputCont);
}

//________________________________________________________________________
void AliAnalysisTaskKinksFilimon::UserExecMix(Option_t *) 
{
	if (fMultiInputHandler) KinkAnalysis();
  // Post output data.
  PostData(1, fOutputCont);
}

//________________________________________________________________________
Int_t AliAnalysisTaskKinksFilimon::KinkAnalysis(/*Option_t **/)
{
  // Main loop
  // Called for each event

	//AliDebug(AliLog::kDebug, Form("Multi=%d MultiMix=%d", esdEvent->GetNumberOfTracks(), esdEventMix->GetNumberOfTracks()));

  AliVEvent* const inputEvent = GetMainEvent();
	if ( !inputEvent ) return -1; // Should return in order to correctly calculate PE
	AliESDEvent* esdEvent = 0x0;
	const AliAODEvent* aodEvent = 0x0;
  if ( !fIsAOD ) esdEvent = dynamic_cast<AliESDEvent*>(inputEvent);
	else aodEvent = dynamic_cast<const AliAODEvent* const>(inputEvent);
	if ( !(esdEvent || aodEvent) ) AliWarning("No rec event found");
	Float_t centralityF = GetRecCentrality(inputEvent, "V0M");
  //if ( (centralityF < 40) || (centralityF > 50) ) return -1;
	
  AliMCEvent* const mcEvent = ( fUseMC && !fIsAOD ) ? MCEvent() : 0;
	const TClonesArray* const arrayAODMC = ( fUseMC && fIsAOD && aodEvent && aodEvent->GetList() ) ? dynamic_cast<const TClonesArray* const>(aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName())) : 0;
	if ( fUseMC && !(mcEvent || arrayAODMC) ) AliWarning("No MC event found");
	//fdotheMCLoopAfterEventCuts TODO: handle this to have MC loop only on ESD/AOD selected events
	
  // Event selection before MC
  if ( fEventSelectionCutsMC && (EventSelection(inputEvent, esdEvent, aodEvent, centralityF) == -1) ) return -1;
	
	// Process MC truth
	const AliVVertex* mcVertex = 0x0;
	Double_t mcVertexXYZArray[3] = { 0, 0, 0 };
	Int_t mcNtracks = -1, mcNprimaries = -1;
	if (mcEvent) {
	  const AliStack* const mcStack = mcEvent->Stack(); //stack of MC events
	  if ( !mcStack ) return -1;
    //AliCFEventGenCuts
	  mcVertex = dynamic_cast<const AliVVertex* const>(mcEvent->GetPrimaryVertex());
	  if (mcVertex) mcVertex->GetXYZ(mcVertexXYZArray);
	  mcNtracks = mcEvent->GetNumberOfTracks(); 
	  mcNprimaries = mcEvent->GetNumberOfPrimaries(); 
	}
	else if (arrayAODMC) {
		mcNtracks = arrayAODMC->GetEntries();
	  //mcNprimaries = mcEvent->GetNumberOfPrimaries();
		const AliAODMCHeader* const mcAODHeader = ( fUseMC && fIsAOD && aodEvent && aodEvent->GetList() ) ? dynamic_cast<const AliAODMCHeader* const >(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName())) : 0;
	  if (mcAODHeader) mcAODHeader->GetVertex(mcVertexXYZArray);
	}
	if (fUseMC) {
	TVector3 mcVertexXYZ(mcVertexXYZArray);
	Double_t mcVertexRPhiZarray[3] = {mcVertexXYZ.Perp(), mcVertexXYZ.Phi(), mcVertexXYZ.Z()};
	fhMCmainVertexRPhiZ->Fill(mcVertexRPhiZarray);
 	fhMCmult->Fill(mcNtracks); 
	fhMCmultPrim->Fill(mcNprimaries); 
	Int_t nMCKaonKinks = 0, nMCPionKinks = 0, nDaughters = -1, MCPdg = -1;
	TLorentzVector mcTrackMomentum;
	const Double_t resMCabsYmax = 1.5;
	Double_t MCPt = -1, MCY = -1, MCEta = -1;
	Bool_t /*isPrimary = kFALSE,*/ isPhysicalPrimary = kFALSE;
	AliVParticle* vMCParticle = 0x0;
 	for (Int_t iMCtrack = 0; iMCtrack < mcNtracks; ++iMCtrack) { //loop on all accepted MC tracks
		if (mcEvent) {
		  AliMCParticle* mcTrack = static_cast<AliMCParticle*>(mcEvent->GetTrack(iMCtrack));
		  if ( !mcTrack ) continue;
		  TParticle* mcParticle = mcTrack->Particle();//mcStack->Particle(iMCtrack);
		  if ( !mcParticle ) continue;
		  //if ( mcParticle != mcStack->Particle(iMCtrack) ) AliFatal("AliMCEvent-AliStack error");
			isPhysicalPrimary = mcEvent->IsPhysicalPrimary(iMCtrack);
		  nDaughters = mcParticle->GetNDaughters();
			vMCParticle = mcTrack;
		}
		else if (arrayAODMC) {
		  /*AliVParticle*/AliAODMCParticle* mcTrack = static_cast<AliAODMCParticle*>(arrayAODMC->At(iMCtrack));
			isPhysicalPrimary = mcTrack->IsPhysicalPrimary();
		  nDaughters = mcTrack->GetNDaughters();
			vMCParticle = mcTrack;
		}
		if ( mcEvent || arrayAODMC ) {
			if ( vMCParticle->M() < fElectronMass ) continue; // Try to avoid vMCParticle->Y() FPE
		  MCY = vMCParticle->Y();//( vMCParticle->E() > TMath::Abs(vMCParticle->Pz()) ) ? mcTrackMomentum.Rapidity() : -999;//Rapidity(mcTrackMomentum);//vMCParticle->Y();//mcTrackMomentum.Rapidity(); //( TMath::Abs(vMCParticle->E()-vMCParticle->Pz()) > std::numeric_limits<Float_t>::epsilon() ) ? vMCParticle->Y() : -999;
		  if (MCY == -999) continue;
		  MCPdg = vMCParticle->PdgCode(); 
		  fhMCpdg->Fill(MCPdg);
		  //if ( !vMCParticle->Charge() ) continue; // Skip neutrals WARNING! also skips neutral resonances!
			mcTrackMomentum.SetXYZT(vMCParticle->Px(), vMCParticle->Py(), vMCParticle->Pz(), vMCParticle->E());
		  MCPt = vMCParticle->Pt();//vMCParticle->Pt();//mcTrackMomentum.Pt();
		  //if (MCPt < fRecKinkCutsKaon->fMinPt) continue; // Try to avoid vMCParticle->Y() FPE
		  MCEta = vMCParticle->Eta(); 
		  //if (MCEta > 20) continue; // Try to avoid vMCParticle->Y() FPE
		  fhMCPtYall->Fill(MCPt, MCY);
		  if ( isPhysicalPrimary ) {
			  fhMCPtYprim->Fill(MCPt, MCY);
		  }
		//else { // All secondaries
		//}
		}
		//if ( MCPt < 0 ) continue;
		switch ( MCPdg ) {
			case kPhi: {
				fhMCPhiPtY->Fill(MCPt, MCY);
		    if ( (MCPt < fRecKinkCutsKaon->fMinPt) || (TMath::Abs(MCY) > resMCabsYmax) ) continue; // Temporary cut 20130420
			  //if ( /*!mcParticle->IsPrimary()*/ /*!mcEvent->IsPhysicalPrimary(iMCtrack) ||*/ !fRecKinkCutsKaon->IsInFiducialKine(mcTrackMomentum) ) continue;
			  fhMCPtYprimPhi->Fill(MCPt, MCY);
				fhMCInvMassPtprimPhi->Fill(mcTrackMomentum.M(), MCPt);
				fhMCPhiInvMassPtYtest->Fill(mcTrackMomentum.M(), MCPt, MCY);
				//if ( fCollisionType != kPP ) fhMCInvMassPtCentprimPhi->Fill(mcTrackMomentum.M(), MCPt, centralityF);
			  //MCPdg == kKPlus ? fhMCPtYprimKPlus->Fill(MCPt, MCY) : fhMCPtYprimKMinus->Fill(MCPt, MCY);
				//Int_t nPhiDaughters = mcParticle->GetNDaughters();
				if (nDaughters != 2) continue;
				fhMCPhi2KaonPtY->Fill(MCPt, MCY);
				AliVParticle* daughter1 = static_cast<AliVParticle*>(mcEvent ? mcEvent->GetTrack(static_cast<AliMCParticle*>(vMCParticle)->GetDaughterFirst()) : arrayAODMC->At(static_cast<AliAODMCParticle*>(vMCParticle)->GetDaughterLabel(0)));
				AliVParticle* daughter2 = static_cast<AliVParticle*>(mcEvent ? mcEvent->GetTrack(static_cast<AliMCParticle*>(vMCParticle)->GetDaughterLast()) : arrayAODMC->At(static_cast<AliAODMCParticle*>(vMCParticle)->GetDaughterLabel(1)));
				Double_t openingAngle = TVector3(daughter1->Px(), daughter1->Py(), daughter1->Pz()).Angle(TVector3(daughter2->Px(), daughter2->Py(), daughter2->Pz()));
				fhMCPhiDecayOpeningAngle->Fill(MCPt, TMath::Cos(openingAngle));
			}
				break;
			case kKstar0: {
				fhMCKstarPtY->Fill(MCPt, MCY);
		    if ( (MCPt < fRecKinkCutsKaon->fMinPt) || (TMath::Abs(MCY) > resMCabsYmax) ) continue; // Temporary cut 20130420
			  //if ( !(  /*mcParticle->IsPrimary()*/ /*mcEvent->IsPhysicalPrimary(iMCtrack) &&*/ fRecKinkCutsKaon->IsInFiducialKine(mcTrackMomentum) )  ) continue;
			  fhMCPtYprimKstar->Fill(MCPt, MCY);
				fhMCInvMassPtprimKstar->Fill(mcTrackMomentum.M(), MCPt);
				//if ( fCollisionType != kPP ) fhMCInvMassPtCentprimKstar->Fill(mcTrackMomentum.M(), MCPt, centralityF);
			  //MCPdg == kKPlus ? fhMCPtYprimKPlus->Fill(MCPt, MCY) : fhMCPtYprimKMinus->Fill(MCPt, MCY);
				//Int_t nKstarDaughters = mcParticle->GetNDaughters();
				if (nDaughters != 2) continue;
				fhMCKstar2KaonPtY->Fill(MCPt, MCY);
				AliVParticle* daughter1 = static_cast<AliVParticle*>(mcEvent ? mcEvent->GetTrack(static_cast<AliMCParticle*>(vMCParticle)->GetDaughterFirst()) : arrayAODMC->At(static_cast<AliAODMCParticle*>(vMCParticle)->GetDaughterLabel(0)));
				AliVParticle* daughter2 = static_cast<AliVParticle*>(mcEvent ? mcEvent->GetTrack(static_cast<AliMCParticle*>(vMCParticle)->GetDaughterLast()) : arrayAODMC->At(static_cast<AliAODMCParticle*>(vMCParticle)->GetDaughterLabel(1)));
				Double_t openingAngle = TVector3(daughter1->Px(), daughter1->Py(), daughter1->Pz()).Angle(TVector3(daughter2->Px(), daughter2->Py(), daughter2->Pz()));
				fhMCKstarDecayOpeningAngle->Fill(MCPt, TMath::Cos(openingAngle));
			}
				break;
			case kLambda1520: {
				fhMCLambdaPtY->Fill(MCPt, MCY);
		    if ( (MCPt < fRecKinkCutsKaon->fMinPt) || (TMath::Abs(MCY) > resMCabsYmax) ) continue; // Temporary cut 20130420
			  //if ( !(  /*mcParticle->IsPrimary()*/ /*mcEvent->IsPhysicalPrimary(iMCtrack) &&*/ fRecKinkCutsKaon->IsInFiducialKine(mcTrackMomentum) )  ) continue;
			  fhMCPtYprimLambda->Fill(MCPt, MCY);
				fhMCInvMassPtprimLambda->Fill(mcTrackMomentum.M(), MCPt);
				//if ( fCollisionType != kPP ) fhMCInvMassPtCentprimLambda->Fill(mcTrackMomentum.M(), MCPt, centralityF);
			  //MCPdg == kKPlus ? fhMCPtYprimKPlus->Fill(MCPt, MCY) : fhMCPtYprimKMinus->Fill(MCPt, MCY);
				//Int_t nLambdaDaughters = mcParticle->GetNDaughters();
				if (nDaughters != 2) continue;
				//fhMCLambda2KaonPtY->Fill(MCPt, MCY);
				AliVParticle* daughter1 = static_cast<AliVParticle*>(mcEvent ? mcEvent->GetTrack(static_cast<AliMCParticle*>(vMCParticle)->GetDaughterFirst()) : arrayAODMC->At(static_cast<AliAODMCParticle*>(vMCParticle)->GetDaughterLabel(0)));
				AliVParticle* daughter2 = static_cast<AliVParticle*>(mcEvent ? mcEvent->GetTrack(static_cast<AliMCParticle*>(vMCParticle)->GetDaughterLast()) : arrayAODMC->At(static_cast<AliAODMCParticle*>(vMCParticle)->GetDaughterLabel(1)));
				Double_t openingAngle = TVector3(daughter1->Px(), daughter1->Py(), daughter1->Pz()).Angle(TVector3(daughter2->Px(), daughter2->Py(), daughter2->Pz()));
				fhMCLambdaDecayOpeningAngle->Fill(MCPt, TMath::Cos(openingAngle));
			}
				break;
			case kKPlus: case kKMinus:
	  { // Study Kaons
		  //if ( !isPhysicalPrimary ||/*!mcEvent->IsPhysicalPrimary(iMCtrack) ||*/ (MCPt < fRecKinkCutsKaon->fMinPt) || (TMath::Abs(MCY) > fRecKinkCutsKaon->fMaxAbsY) ) continue; // Temporary cut 20130420
			//if ( /*!mcParticle->IsPrimary()*/ !mcEvent->IsPhysicalPrimary(iMCtrack) || !fRecKinkCutsKaon->IsInFiducialKine(mcTrackMomentum) ) continue;
			if ( !isPhysicalPrimary || !fRecKinkCutsKaon->IsInFiducialKine(mcTrackMomentum) ) continue;
			//if ( mcEvent->IsPhysicalPrimary(iMCtrack) /*mcParticle->IsPrimary()*/ /*&& fRecKinkCutsKaon->IsInFiducialKine(mcTrackMomentum)*/ ) {
			fhMCPtYprimKaon->Fill(MCPt, MCY);
			MCPdg == kKPlus ? fhMCPtYprimKPlus->Fill(MCPt, MCY) : fhMCPtYprimKMinus->Fill(MCPt, MCY);
			MCPdg == kKPlus ? fhMCPtprimKPlus->Fill(MCPt) : fhMCPtprimKMinus->Fill(MCPt);
			if ( fCollisionType != kPP ) {
				fhMCPtYCentprimKaon->Fill(MCPt, MCY, centralityF);
			  MCPdg == kKPlus ? fhMCPtYCentprimKPlus->Fill(MCPt, MCY, centralityF) : fhMCPtYCentprimKMinus->Fill(MCPt, MCY, centralityF);
			}
			//}
			if (arrayAODMC) continue; // No trackrefs available in AOD
			AliMCParticle* mcTrack = static_cast<AliMCParticle*>(mcEvent->GetTrack(iMCtrack));
			// Study Kaon mother (eventually move this analysis as resonance daughter)
			Int_t kaonMotherIdx = mcTrack->GetMother();
			if ( (kaonMotherIdx != -1 ) && (kaonMotherIdx < mcNprimaries) ) {
				AliVParticle* kaonMother = mcEvent->GetTrack(kaonMotherIdx);
				Int_t kaonMotherPdg = kaonMother->PdgCode();
				fhMCkaonMotherPdg->Fill(kaonMotherPdg);
				/*if ( TMath::Abs(kaonMotherPdg) == kPhi ) fhMCPhi2KaonPtY->Fill(kaonMother->Pt(), kaonMother->Y());
				else if ( TMath::Abs(kaonMotherPdg) == kKstar0 ) fhMCKstar2KaonPtY->Fill(kaonMother->Pt(), kaonMother->Y());*/
			}
			// Study Kaon daughters and decays // This does not work properly yet
			//if (mcParticle->GetNDaughters != 1) continue;
			AliTrackReference* lastTrackReference = 0x0;
			TVector3 decay3Momentum;
			TVector3 decay3Position;
			Int_t nTrackRefs = GetLastTrackRef(mcTrack, lastTrackReference, decay3Momentum, decay3Position);
			if ( (nTrackRefs < 1) || !lastTrackReference ) continue; // At least 1 trackref for decay
			Double_t rPhiZ[3] = {decay3Position.Perp(), decay3Position.Phi(), decay3Position.Z()};
			fhMCKaonKinkDecaysRPhiZ->Fill(rPhiZ);
		  //if ( decay3Position.Y() > 0 ) continue; // WARNING!!! MSS request to check LHC12f1a symmetry
			fhMCKaonKinkDecaysTrackLength->Fill(lastTrackReference->GetLength());
			fhMCKaonKinkDecaysTrackTime->Fill(1000000 * lastTrackReference->GetTime());
			Int_t firstDaughterIdx=0, lastDaughterIdx=0, nDaughtersCalc=0, nChargedDaughters=0;
			if ( ( (nDaughtersCalc=GetDaughterIdx(mcTrack, mcNtracks, firstDaughterIdx, lastDaughterIdx)) < 1 ) /*|| (nDaughtersCalc > 3 )*/ ) continue;
			for (Int_t iDaughterIdx = firstDaughterIdx; iDaughterIdx <= lastDaughterIdx; ++iDaughterIdx) {
				AliMCParticle* daughterMC = static_cast<AliMCParticle*>(mcEvent->GetTrack(iDaughterIdx));
				if ( !daughterMC ) continue;
				if ( daughterMC->Charge() ) ++nChargedDaughters;
				//if ( nChargedDaughters > 1 ) break;
			  TParticle* daughterPart = daughterMC->Particle();
        if ( !daughterPart || (daughterPart->GetUniqueID() != kPDecay) ) break;
				Int_t daughterPdg = daughterMC->PdgCode();
				if ( !fRecKinkCutsKaon->IsInFiducial(decay3Position) ) continue;
				if ( TMath::Abs(daughterPdg) == kMuonMinus ) fhMCKaon2MuonPtY->Fill(MCPt, MCY);
				else if ( TMath::Abs(daughterPdg) == kPiPlus ) fhMCKaon2PionPtY->Fill(MCPt, MCY);
				TVector3 daughterMC3Momentum(daughterMC->Px(), daughterMC->Py(), daughterMC->Pz());
				fhMCKinkKaonQt->Fill(daughterMC3Momentum.Perp(decay3Momentum));
				fhMCKinkKaonAngle->Fill(daughterMC3Momentum.Angle(decay3Momentum));
			//if ( !firstDaughter /*|| (firstDaughter->GetUniqueID() != kPDecay)*/ || ( TMath::Abs(firstDaughterPdg) != kMuonMinus ) || ( TMath::Abs(firstDaughterPdg) != kPiPlus )) continue; // Kink decays
		  }
			if ( /*( nChargedDaughters == 1 ) &&*/ mcEvent->IsPhysicalPrimary(iMCtrack) && fRecKinkCutsKaon->IsInFiducial(decay3Position) ) {
				++nMCKaonKinks;
				fhMCPtYkinkKaonFiducial->Fill(MCPt, MCY); // GA, RE, PE can eventually be computed with this
			  MCPdg == kKPlus ? fhMCPtYkinkKPlusFiducial->Fill(MCPt, MCY) : fhMCPtYkinkKMinusFiducial->Fill(MCPt, MCY);
			}
		} // Study Kaons end
		break;
			case kPiPlus: case kPiMinus:
		{ // Study Pions
		  //if ( !isPhysicalPrimary /*!mcEvent->IsPhysicalPrimary(iMCtrack)*/ || (MCPt < fRecKinkCutsKaon->fMinPt) || (TMath::Abs(MCY) > fRecKinkCutsKaon->fMaxAbsY) ) continue; // Temporary cut 20130420
			//if ( /*!mcParticle->IsPrimary()*/ !mcEvent->IsPhysicalPrimary(iMCtrack) || !fRecKinkCutsKaon->IsInFiducialKine(mcTrackMomentum) ) continue;
			if ( !isPhysicalPrimary || !fRecKinkCutsKaon->IsInFiducialKine(mcTrackMomentum) ) continue;
			//if ( mcEvent->IsPhysicalPrimary(iMCtrack) /*mcParticle->IsPrimary()*/ /*&& fRecKinkCutsKaon->IsInFiducialKine(mcTrackMomentum)*/ ) {
			fhMCPtYprimPion->Fill(MCPt, MCY);
			MCPdg == kPiPlus ? fhMCPtYprimPiPlus->Fill(MCPt, MCY) : fhMCPtYprimPiMinus->Fill(MCPt, MCY);
			if ( fCollisionType != kPP ) {
				fhMCPtYCentprimPion->Fill(MCPt, MCY, centralityF);
			  MCPdg == kPiPlus ? fhMCPtYCentprimPiPlus->Fill(MCPt, MCY, centralityF) : fhMCPtYCentprimPiMinus->Fill(MCPt, MCY, centralityF);
			}
			//}
			if (arrayAODMC) continue; // No trackrefs available in AOD
			AliMCParticle* mcTrack = static_cast<AliMCParticle*>(mcEvent->GetTrack(iMCtrack));
			// Study Pion mother (eventually move this analysis as resonance daughter)
			Int_t pionMotherIdx = mcTrack->GetMother();
			if ( (pionMotherIdx != -1 ) && (pionMotherIdx < mcNprimaries) ) {
				AliVParticle* pionMother = mcEvent->GetTrack(pionMotherIdx);
				Int_t pionMotherPdg = pionMother->PdgCode();
				fhMCpionMotherPdg->Fill(pionMotherPdg);
				/*if ( TMath::Abs(pionMotherPdg) == kPhi ) fhMCPhi2PionPtY->Fill(pionMother->Pt(), pionMother->Y());
				else if ( TMath::Abs(pionMotherPdg) == kKstar0 ) fhMCKstar2PionPtY->Fill(pionMother->Pt(), pionMother->Y());*/
			}
			// Study Pion daughters and decays // This does not work properly yet
			//if (mcParticle->GetNDaughters != 1) continue;
			AliTrackReference* lastTrackReference = 0x0;
			TVector3 decay3Momentum;
			TVector3 decay3Position;
			Int_t nTrackRefs = GetLastTrackRef(mcTrack, lastTrackReference, decay3Momentum, decay3Position);
			if ( (nTrackRefs < 1) || !lastTrackReference ) continue; // At least 1 trackref for decay
			Double_t rPhiZ[3] = {decay3Position.Perp(), decay3Position.Phi(), decay3Position.Z()};
			fhMCPionKinkDecaysRPhiZ->Fill(rPhiZ);
		  //if ( decay3Position.Y() > 0 ) continue; // WARNING!!! MSS request to check LHC12f1a symmetry
			//fhMCPionKinkDecaysTrackLength->Fill(lastTrackReference->GetLength());
			//fhMCPionKinkDecaysTrackTime->Fill(1000000 * lastTrackReference->GetTime());
			Int_t firstDaughterIdx=0, lastDaughterIdx=0, nDaughtersCalc=0, nChargedDaughters=0;
			if ( (nDaughtersCalc=GetDaughterIdx(mcTrack, mcNtracks, firstDaughterIdx, lastDaughterIdx)) < 1 ) continue;
			for (Int_t iDaughterIdx = firstDaughterIdx; iDaughterIdx <= lastDaughterIdx; ++iDaughterIdx) {
				AliMCParticle* daughterMC = static_cast<AliMCParticle*>(mcEvent->GetTrack(iDaughterIdx));
				if ( !daughterMC ) continue;
				if ( daughterMC->Charge() ) ++nChargedDaughters;
				//if ( nChargedDaughters > 1 ) break;
			  TParticle* daughterPart = daughterMC->Particle();
        if ( !daughterPart || (daughterPart->GetUniqueID() != kPDecay) ) break;
				Int_t daughterPdg = daughterMC->PdgCode();
				if ( !fRecKinkCutsKaon->IsInFiducial(decay3Position) ) continue;
				if ( TMath::Abs(daughterPdg) == kMuonMinus ) fhMCPion2MuonPtY->Fill(MCPt, MCY);
				//else if ( TMath::Abs(daughterPdg) == kPiPlus ) fhMCPion2PionPtY->Fill(MCPt, MCY);
				TVector3 daughterMC3Momentum(daughterMC->Px(), daughterMC->Py(), daughterMC->Pz());
				fhMCKinkPionQt->Fill(daughterMC3Momentum.Perp(decay3Momentum));
				fhMCKinkPionAngle->Fill(daughterMC3Momentum.Angle(decay3Momentum));
		  }
			if ( /*( nChargedDaughters == 1 ) &&*/ mcEvent->IsPhysicalPrimary(iMCtrack) && fRecKinkCutsKaon->IsInFiducial(decay3Position) ) {
				++nMCPionKinks;
				fhMCPtYkinkPionFiducial->Fill(MCPt, MCY); // GA, RE, PE can eventually be computed with this
			  MCPdg == kPiPlus ? fhMCPtYkinkPiPlusFiducial->Fill(MCPt, MCY) : fhMCPtYkinkPiMinusFiducial->Fill(MCPt, MCY);
			}
		} // Study Pions end
		break;
		} // PDG switch
	} // MC particle loop
	if (nMCPionKinks) fhMCKaonOverPionEbE->Fill(static_cast<Float_t>(nMCKaonKinks)/nMCPionKinks);
	} // MC analysis
  
  // Event selection after MC
  if ( !fEventSelectionCutsMC && (EventSelection(inputEvent, esdEvent, aodEvent, centralityF) == -1) ) return -1;
  
	// Rec
	TLorentzVector resonance4Momentum, mother4Momentum, mother4MomentumPartnerKink;
	Int_t nESDKaonKinks = 0, nESDPionKinks = 0, nESDUnknownKinks = 0;
	
	if (esdEvent) { // ESD analysis

  /*
  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < esdEvent->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = esdEvent->GetTrack(iTracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }

    fhRecAllPtY->Fill(track->Pt());
  } //track loop */
	
	//const AliMultiplicity* const esdMultiplicity = esdEvent->GetMultiplicity();
	Int_t nKinks = esdEvent->GetNumberOfKinks();
	fhESDnKinks->Fill(nKinks);
	for (Int_t iKink = 0; iKink < nKinks; ++iKink) { // ESD kinks loop
	  AliESDkink* kink = esdEvent->GetKink(iKink);
		TVector3 kinkVertex(kink->GetPosition());
		Double_t rPhiZ[3] = {kinkVertex.Perp(), kinkVertex.Phi(), kinkVertex.Z()};
		fhESDAllKinkDecaysRPhiZ->Fill(rPhiZ);
		//if ( kinkVertex.Y() > 0 ) continue; // WARNING!!! MSS request to check LHC12f1a symmetry
		fhESDKinkQt->Fill(kink->GetQt());
		fhESDKinkAngle->Fill(kink->GetAngle(2));
		fhESDKinkDCA->Fill(kink->GetDistance());
		if (!fRecKinkCutsKaon->IsInFiducial(kinkVertex)) continue;
		AliVTrack* kinkMother = 0x0;
		AliVTrack* kinkDaughter = 0x0;
		Int_t kinkMotherPDG = 0;
		AliPID::EParticleType kinkMotherPID = AliPID::kUnknown;
		Double_t kinkMotherSign = 0;
		Double_t kinkMotherPt = 0;
		Double_t kinkMotherY = 0;
		//Double_t kinkMotherEta = 0;
		/*Int_t FoundGoodKaonKink(AliVTrack* const kinkMother, TLorentzVector& mother4Momentum, const Float_t centralityF) {
		}*/
		if ( fRecKinkCutsKaon->IsGoodKaonKink(kink, esdEvent, kinkMother, kinkDaughter) /*fRecKinkCutsKaon->IsInQtLimits(kink->GetQt())*/) {
			/*fRecKinkCutsKaon->IsGoodKaonKink(kink, esdEvent, kinkMother, kinkDaughter);
			if ( !kinkMother || !kinkDaughter) continue;*/
		  ++nESDKaonKinks;
			kinkMotherPID = AliPID::kKaon;
		  kinkMotherSign = kinkMother->Charge();//GetSign();
	    kinkMotherPDG = (kinkMotherSign > 0 ) ? kKPlus : kKMinus;
		  mother4Momentum.SetVectM(TVector3(kinkMother/*TPCParam*/->Px(), kinkMother/*TPCParam*/->Py(), kinkMother/*TPCParam*/->Pz()), fKaonMass); // Kaon kink hypothesis
		  kinkMotherPt = mother4Momentum.Pt();//kinkMother->Pt();//mother4Momentum.Pt();
		  kinkMotherY = mother4Momentum.Rapidity();//kinkMother->Y();//mother4Momentum.Rapidity();
		  if ( !fRecKinkCutsKaon->IsInFiducialKine(mother4Momentum)) continue; // WARNING! Enabling this cut deviates from MSS
		  //if ( (kinkMotherPt < fRecKinkCutsKaon->fMinPt) || ( TMath::Abs(kinkMotherY) > fRecKinkCutsKaon->fMaxAbsY) ) continue; // Temporary cut 20130420, deviates from MSS, checked OK 20130620 after fixing kinkMotherPt/Y bug 
		  fhESDKaonKinkDecaysRPhiZ->Fill(rPhiZ);
	    fhESDKaonKinkPtY->Fill(kinkMotherPt, kinkMotherY);
			fhESDKaonKinkQt->Fill(kink->GetQt());
			fhESDKaonKinkAngle->Fill(kink->GetAngle(2));
		  (kinkMotherSign > 0 ) ? fhESDKPlusKinkPtY->Fill(kinkMotherPt, kinkMotherY) : fhESDKMinusKinkPtY->Fill(kinkMotherPt, kinkMotherY);
			(kinkMotherSign > 0 ) ? fhESDKPlusKinkPt->Fill(kinkMotherPt) : fhESDKMinusKinkPt->Fill(kinkMotherPt);
			if ( fCollisionType != kPP ) {
				fhESDKaonKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF);
				(kinkMotherSign > 0 ) ? fhESDKPlusKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF) : fhESDKMinusKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF);
			}
	    fhESDKaonKinkPtEta->Fill(mother4Momentum.Pt(), mother4Momentum.Eta());
		  fhESDMomentumTPCSignalKaonKinks->Fill(kinkMother->P()/*mother4Momentum.P()*/, kinkMother->GetTPCsignal());
		  fhESDKaonKinksMotherVSDaughterTPCncls->Fill(kinkMother->GetTPCNcls(), kinkDaughter->GetTPCNcls());
		  fhESDKaonKinksMotherAndDaughterTPCncls->Fill(kinkMother->GetTPCNcls() + kinkDaughter->GetTPCNcls());
		}
		else if ( fRecKinkCutsKaon->IsGoodPionKink(kink, esdEvent, kinkMother, kinkDaughter) ) { // First look for pion kink that has narrow cut and also contaminates kaons
			++nESDPionKinks;
			kinkMotherPID = AliPID::kPion;
		  kinkMotherSign = kinkMother->Charge();//GetSign();
			kinkMotherPDG = (kinkMotherSign > 0 ) ? kPiPlus : kPiMinus;
		  mother4Momentum.SetVectM(TVector3(kinkMother/*TPCParam*/->Px(), kinkMother/*TPCParam*/->Py(), kinkMother/*TPCParam*/->Pz()), fPionMass); // Pion kink hypothesis
		  kinkMotherPt = mother4Momentum.Pt();//kinkMother->Pt();//mother4Momentum.Pt();
		  kinkMotherY = mother4Momentum.Rapidity();//kinkMother->Y();//mother4Momentum.Rapidity();
		  if ( !fRecKinkCutsKaon->IsInFiducialKine(mother4Momentum)) continue; // WARNING! Enabling this cut deviates from MSS
		  //if ( (kinkMotherPt < fRecKinkCutsKaon->fMinPt) || ( TMath::Abs(kinkMotherY) > fRecKinkCutsKaon->fMaxAbsY) ) continue; // Temporary cut 20130420, deviates from MSS, checked OK 20130620 after fixing kinkMotherPt/Y bug 
		  fhESDPionKinkDecaysRPhiZ->Fill(rPhiZ);
			fhESDPionKinkPtY->Fill(kinkMotherPt, kinkMotherY);
			fhESDPionKinkQt->Fill(kink->GetQt());
			fhESDPionKinkAngle->Fill(kink->GetAngle(2));
		  (kinkMotherSign > 0 ) ? fhESDPiPlusKinkPtY->Fill(kinkMotherPt, kinkMotherY) : fhESDPiMinusKinkPtY->Fill(kinkMotherPt, kinkMotherY);
			if ( fCollisionType != kPP ) {
				fhESDPionKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF);
				(kinkMotherSign > 0 ) ? fhESDPiPlusKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF) : fhESDPiMinusKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF);
			}
	    //fhESDPionKinkPtEta->Fill(mother4Momentum.Pt(), mother4Momentum.Eta());
		  fhESDMomentumTPCSignalPionKinks->Fill(kinkMother->P()/*mother4Momentum.P()*/, kinkMother->GetTPCsignal());
		  fhESDPionKinksMotherVSDaughterTPCncls->Fill(kinkMother->GetTPCNcls(), kinkDaughter->GetTPCNcls());
		  fhESDPionKinksMotherAndDaughterTPCncls->Fill(kinkMother->GetTPCNcls() + kinkDaughter->GetTPCNcls());
		}
		else {
			++nESDUnknownKinks;
			kinkMotherPID = AliPID::kUnknown;
			if ( !kinkMother ) continue;
		  kinkMotherSign = kinkMother->Charge();//GetSign();
			//kinkMotherPDG = 0;
			fhESDUnknownKinkPtY->Fill(kinkMother->Pt(), kinkMother->Y());
			/*if (mcEvent) {
			 * Used to define proper PIDeff
			}*/
			continue; // Neither Pion or Kaon "good" kink found
		}
		fhESDAllKinkPtY->Fill(kinkMother->Pt(), kinkMother->Y());
		for (Int_t jKink = iKink + 1; jKink < nKinks; ++jKink) { // Look for phi/kstar->2kinks
	    AliESDkink* partnerKink = esdEvent->GetKink(jKink);
			AliVTrack* partnerKinkMother = 0x0;
			AliVTrack* partnerKinkDaughter = 0x0;
		  if ( !fRecKinkCutsKaon->IsInFiducial(TVector3(partnerKink->GetPosition())) ) continue;
			if ( kinkMotherPID == AliPID::kPion ) { // Current kink is Pion, look for Kaon kink to form Kstar 
		    if ( !fRecKinkCutsKaon->IsGoodKaonKink(partnerKink, esdEvent, partnerKinkMother, partnerKinkDaughter) ) continue;
			  if ( (kinkMother == partnerKinkMother) || (kinkDaughter == partnerKinkDaughter) ) continue;
			  mother4MomentumPartnerKink.SetVectM(TVector3(partnerKinkMother/*TPCParam*/->Px(), partnerKinkMother/*TPCParam*/->Py(), partnerKinkMother/*TPCParam*/->Pz()), fKaonMass); // Kaon kink hypothesis
			  //if (mother4Momentum.Angle(mother4MomentumPartnerKink.Vect()) < 0.6) continue; // Apply phi->K+K- "magic" angle cut (P-dependent?)
			  resonance4Momentum = mother4Momentum + mother4MomentumPartnerKink;
		    //if ( !fRecKinkCutsKaon->IsInFiducialKine(mother4MomentumPartnerKink)) continue; // WARNING! Enabling this cut deviates from MSS
			  //if ( resonance4Momentum.Pt() < fRecKinkCutsKaon->fMinPt ) continue; // Remove
			  if ( kinkMotherSign * partnerKinkMother->Charge()/*GetSign()*/ < 0) fhRecKstarInvMassBothKinks->Fill(resonance4Momentum.M()); // kstar->2kinks
			  else fhRecKstarInvMassBothKinksLikeSign->Fill(resonance4Momentum.M()); // kstar->2kinks likesign bgr
			}
			else if ( kinkMotherPID == AliPID::kKaon ) {  // Current kink is Kaon, look for Kaon or Pion kink to form Phi or Kstar 
		    if ( fRecKinkCutsKaon->IsGoodPionKink(partnerKink, esdEvent, partnerKinkMother, partnerKinkDaughter) ) {
			    if ( (kinkMother == partnerKinkMother) || (kinkDaughter == partnerKinkDaughter) ) continue;
			    mother4MomentumPartnerKink.SetVectM(TVector3(partnerKinkMother/*TPCParam*/->Px(), partnerKinkMother/*TPCParam*/->Py(), partnerKinkMother/*TPCParam*/->Pz()), fPionMass); // Pion kink hypothesis
			    //if (mother4Momentum.Angle(mother4MomentumPartnerKink.Vect()) < 0.6) continue; // Apply phi->K+K- "magic" angle cut (P-dependent?)
			    resonance4Momentum = mother4Momentum + mother4MomentumPartnerKink;
		      //if ( !fRecKinkCutsKaon->IsInFiducialKine(mother4MomentumPartnerKink)) continue; // WARNING! Enabling this cut deviates from MSS
			    //if ( resonance4Momentum.Pt() < fRecKinkCutsKaon->fMinPt ) continue; // Remove
			    if ( kinkMotherSign * partnerKinkMother->Charge()/*GetSign()*/ < 0) fhRecKstarInvMassBothKinks->Fill(resonance4Momentum.M()); // kstar->2kinks
			    else fhRecKstarInvMassBothKinksLikeSign->Fill(resonance4Momentum.M()); // kstar->2kinks likesign bgr
				}
				else if ( fRecKinkCutsKaon->IsGoodKaonKink(partnerKink, esdEvent, partnerKinkMother, partnerKinkDaughter) ) {
			    if ( (kinkMother == partnerKinkMother) || (kinkDaughter == partnerKinkDaughter) ) continue;
			    mother4MomentumPartnerKink.SetVectM(TVector3(partnerKinkMother/*TPCParam*/->Px(), partnerKinkMother/*TPCParam*/->Py(), partnerKinkMother/*TPCParam*/->Pz()), fKaonMass); // Kaon kink hypothesis
			    //if (mother4Momentum.Angle(mother4MomentumPartnerKink.Vect()) < 0.6) continue; // Apply phi->K+K- "magic" angle cut (P-dependent?)
			    resonance4Momentum = mother4Momentum + mother4MomentumPartnerKink;
		      //if ( !fRecKinkCutsKaon->IsInFiducialKine(mother4MomentumPartnerKink)) continue; // WARNING! Enabling this cut deviates from MSS
			    //if ( resonance4Momentum.Pt() < fRecKinkCutsKaon->fMinPt ) continue; // Remove
			    if ( kinkMotherSign * partnerKinkMother->Charge()/*GetSign()*/ < 0) fhRecPhiInvMassBothKinks->Fill(resonance4Momentum.M()); // phi->2kinks
			    else fhRecPhiInvMassBothKinksLikeSign->Fill(resonance4Momentum.M()); // phi->2kinks likesign bgr
				}
			}
		}
		if (mcEvent) { // Check MC truth
			Int_t motherMClabel = TMath::Abs(kinkMother->GetLabel());
			AliMCParticle* kinkMotherMC = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(motherMClabel));
			if ( kinkMotherMC ) {
				CheckMCtruth(kinkMother, mother4Momentum, kinkMotherPID, kinkMotherMC, mcEvent->IsPhysicalPrimary(motherMClabel));
#if 1 // TODO: Move into CheckMCtruth or make new function
			  // Study Kaon mother
			  Int_t kaonMotherIdx = kinkMotherMC->GetMother();
			  if ( (kaonMotherIdx != -1 ) && (kaonMotherIdx < mcEvent->GetNumberOfPrimaries()) ) {
				  AliVParticle* kaonMother = mcEvent->GetTrack(kaonMotherIdx);
				  Int_t kaonMotherPdg = kaonMother->PdgCode();
				  //fhMCkaonMotherPdg->Fill(kaonMotherPdg);
				  switch ( TMath::Abs(kaonMotherPdg) ) {
						case kPhi:
							fhESDPhi2KaonPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
							break;
						case kKstar0:
							if ( kinkMotherPID == AliPID::kKaon ) fhESDKstar2KaonPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
							break;
						case kLambda1520:
							fhESDLambda2KaonPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
							break;
						//default:
							//fhMCFakeResonance2KaonPtYFakeMC->Fill(kinkMotherPt, kinkMotherY);
				  }
			  }
#endif	
			} // kinkMotherMC
		} // MC truth check
		if (!fEsdTrackCutsPartner) continue; // Skip resonance analysis if no ESD track cuts provided for partner
		KinkResonanceAnalysis(esdEvent, kinkMother, mother4Momentum, kinkMotherPID, 0x0, mcEvent);
	} // ESD kink loop
	if (nESDPionKinks) fhESDKaonOverPionEbE->Fill(static_cast<Float_t>(nESDKaonKinks)/nESDPionKinks);
	} // ESD analysis
	else if (aodEvent) { // AOD analysis
	// Do "physics selection" here?
	//Float_t aodMultiplicity =	aodEvent->GetVZEROEqMultiplicity(0);
	Int_t aodNvertices = aodEvent->GetNumberOfVertices();
	fhAODnVertices->Fill(aodNvertices);
	Int_t nKinkVertices = 0;
	TLorentzVector /*ESD resonance4Momentum, mother4Momentum, mother4MomentumPartnerKink,*/ partnerTrack4MomentumKaonHypothesis, partnerTrack4MomentumPionHypothesis, partnerTrack4MomentumProtonHypothesis;
	for (Int_t iVertex = 0; iVertex < aodNvertices; ++iVertex) {
		AliAODVertex* aodVertex = aodEvent->GetVertex(iVertex);
		if ( !aodVertex || (aodVertex->GetType() != AliAODVertex::kKink) ) continue;
		++nKinkVertices;
		Int_t nKinkDaughters = aodVertex->GetNDaughters();
		fhAODnKinkDaughters->Fill(nKinkDaughters);
		//if (nKinkDaughters != 1) continue; // Remove, not needed for kKink
		TVector3 kinkVertex(aodVertex->GetX(), aodVertex->GetY(), aodVertex->GetZ());
		Double_t rPhiZ[3] = {kinkVertex.Perp(), kinkVertex.Phi(), kinkVertex.Z()};
		fhESDAllKinkDecaysRPhiZ->Fill(rPhiZ);
		//fhAODKinkQt->Fill(kink->GetQt());
		//fhAODKinkAngle->Fill(kink->GetAngle(2));
		//fhAODKinkDCA->Fill(kink->GetDistance());
		if (!fRecKinkCutsKaon->IsInFiducial(kinkVertex)) continue;
if (!fMultiInputHandler) { // This part crashes with AliMultiInputEventHandler
		AliAODTrack* kinkMother = dynamic_cast<AliAODTrack*>(aodVertex->GetParent());
		if (!kinkMother) continue;
		AliAODTrack* kinkDaughter = dynamic_cast<AliAODTrack*>(aodVertex->GetDaughter(0));
		if (!kinkDaughter) continue;
		if ( !kinkMother->IsPrimaryCandidate() || kinkDaughter->IsPrimaryCandidate() ) continue;
		// TODO: Put standard track BIT cuts here!
		//if ( !kinkMother->TestFilterBit(AliAODTrack::kTrkTPCOnly) ) continue; // AliAODTrack::kTrkGlobal not good for kinks, AliAODTrack::kTrkTPCOnly no effect with respect to PassesESDTrackCuts
		if ( !kinkMother->TestFilterMask(fAODFilterMap) ) continue;	
		//if (!PassesESDTrackCuts(kinkMother, fEsdTrackCutsKinkMother)) continue;
		Double_t b = aodEvent->GetMagneticField();
		if ( b == -999.0) continue;
		Double_t kinkMotherPt = kinkMother->Pt(); //mother4Momentum.Pt()
		//Double_t kinkMotherY = kinkMother->Y(); //mother4Momentum.Rapidity()
		// Kaon kink analysis
		TVector3 motherVertex(kinkMother->Xv(), kinkMother->Yv(), kinkMother->Zv());
		Double_t rPhiZmotherVertex[3] = {motherVertex.Perp(), motherVertex.Phi(), motherVertex.Z()};
		fhAODmotherVertexRPhiZ->Fill(rPhiZmotherVertex);
		TVector3 mother3Momentum(kinkMother->Px(), kinkMother->Py(), kinkMother->Pz());
		/*TVector3 motherVertexDCA(kinkMother->XAtDCA(), kinkMother->YAtDCA(), kinkMother->ZAtDCA());
		Double_t rPhiZmotherVertexDCA[3] = {motherVertexDCA.Perp(), motherVertexDCA.Phi(), motherVertexDCA.Z()};
		fhAODmotherVertexDCARPhiZ->Fill(rPhiZmotherVertexDCA); // To be removed!!! */
		// Calculate DCA to the kink vertex x-plane and get mother and daughter momentum there
		//TVector3 mother3MomentumDCA(kinkMother->PxAtDCA(), kinkMother->PyAtDCA(), kinkMother->PzAtDCA());
    AliESDtrack kinkMotherESD(kinkMother);
    AliESDtrack kinkDaughterESD(kinkDaughter);
    Double_t xKinkMother, xKinkDaughter, minDist = kinkMotherESD.GetDCA(&kinkDaughterESD, b, xKinkMother, xKinkDaughter);
		if ( minDist > 2.0 ) continue;
		Double_t mother3MomentumKinkVtxArray[3] = {0, 0, 0};
		Double_t daughter3MomentumKinkVtxArray[3] = {0, 0, 0};
		if ( !kinkMotherESD.GetPxPyPzAt(xKinkMother, b, mother3MomentumKinkVtxArray) ) continue;
		if ( !kinkDaughterESD.GetPxPyPzAt(xKinkDaughter, b, daughter3MomentumKinkVtxArray) ) continue;
		TVector3 mother3MomentumDCA(mother3MomentumKinkVtxArray);
		TVector3 daughter3MomentumDCA(daughter3MomentumKinkVtxArray);
		//if ( !IsGoodKaonKink(kink, esdEvent, kinkMother, kinkDaughter) ) continue;
		//mother4Momentum.SetVectM(TVector3(kinkMother/*TPCParam*/->Px(), kinkMother/*TPCParam*/->Py(), kinkMother/*TPCParam*/->Pz()), fKaonMass); // Kaon kink hypothesis
		//TVector3 daughter3Momentum(daughter3MomentumKinkVtxArray/*kinkDaughter->Px(), kinkDaughter->Py(), kinkDaughter->Pz()*/);
		Double_t qt = daughter3MomentumDCA.Perp(mother3MomentumDCA), kinkAngleRad = daughter3MomentumDCA.Angle(mother3MomentumDCA), kinkAngleDeg = TMath::RadToDeg() * kinkAngleRad;
		fhESDKinkQt->Fill(qt);
		fhESDKinkAngle->Fill(kinkAngleDeg);
		fhESDKinkDCA->Fill(minDist);
		AliPID::EParticleType kinkMotherPID = AliPID::kUnknown;
		// Temp kink PID
		/*if (  ) kinkMotherPID = AliPID::kKaon;
		else kinkMotherPID = AliPID::kPion;*/
		// Copied from ESD
	  //ESD AliESDkink* kink = esdEvent->GetKink(iKink);
		//ESD TVector3 kinkVertex(kink->GetPosition());
		//ESD Double_t rPhiZ[3] = {kinkVertex.Perp(), kinkVertex.Phi(), kinkVertex.Z()};
		//ESD fhESDAllKinkDecaysRPhiZ->Fill(rPhiZ);
		//if ( kinkVertex.Y() > 0 ) continue; // WARNING!!! MSS request to check LHC12f1a symmetry
		//ESD fhESDKinkQt->Fill(kink->GetQt());
		//ESD fhESDKinkAngle->Fill(kink->GetAngle(2));
		//ESD fhESDKinkDCA->Fill(kink->GetDistance());
		if (!fRecKinkCutsKaon->IsInFiducial(kinkVertex)) continue;
		//ESD AliESDtrack* kinkMother = 0x0;
		//ESD AliESDtrack* kinkDaughter = 0x0;
		Int_t kinkMotherPDG = 0;
		//ESD AliPID::EParticleType kinkMotherPID = AliPID::kUnknown;
		Double_t kinkMotherSign = kinkMother->Charge();
		//ESD Double_t kinkMotherPt = 0;
		Double_t kinkMotherY = 0;
		//Double_t kinkMotherEta = 0;
		Double_t maxDecAngKmu = fRecKinkCutsKaon->fMaxDecayAngleCurveKmu->Eval(mother3MomentumDCA.Mag(), 0., 0., 0.);
		Double_t maxDecAngpimu = fRecKinkCutsKaon->fMaxDecayAngleCurvePimu->Eval(mother3MomentumDCA.Mag(), 0., 0., 0.);
		//if ( kinkAngle < maxDecAngpimu ) return kFALSE; // Possibly pion, reject above theoretical pion kink curve
		//TLorentzVector mother4MomentumDCA, daughter4MomentumDCA;
		//mother4MomentumDCA.SetVectM(mother3MomentumDCA, fKaonMass);
		//daughter4MomentumDCA.SetVectM(daughter3MomentumDCA, fMuonMass);
		
		TVector3 transferedMom = mother3MomentumDCA-daughter3MomentumDCA;
		Float_t energyDaughterMu = TMath::Sqrt(daughter3MomentumDCA.Mag()*daughter3MomentumDCA.Mag()+fMuonMass*fMuonMass);
    Float_t invMassKmu = (energyDaughterMu+transferedMom.Mag())*(energyDaughterMu+transferedMom.Mag())-mother3MomentumDCA.Mag()*mother3MomentumDCA.Mag();
		if (invMassKmu > 0) invMassKmu = TMath::Sqrt(invMassKmu);
		else continue;
     /*Double_t tpcNClMax = -51.67 + (11./12.) * kinkVertex.Perp();
     Double_t tpcNClMin  = -85.5 + (65./95.) * kinkVertex.Perp();
		if ( (kinkMother->GetTPCNcls() < tpcNClMin ) || (kinkMother->GetTPCNcls() > tpcNClMax ) ) continue;*/ // MSS compatibility, no real effect observed

		//Double_t invMassKmu = (mother4MomentumDCA+daughter4MomentumDCA).M();
		fhESDInvMassKinkDecayMuNu->Fill(invMassKmu);
		if ( fRecKinkCutsKaon->IsInQtLimits(qt) && ( kinkAngleDeg > maxDecAngpimu * 1.2) && (kinkAngleDeg > fRecKinkCutsKaon->fMinAngleKaon) && ( (mother3MomentumDCA.Mag() > fRecKinkCutsKaon->fMinKaonMomentum) ? (kinkAngleDeg < maxDecAngKmu * 0.98) : kTRUE ) && (invMassKmu > 0.3) && (invMassKmu < fRecKinkCutsKaon->fInvMassKmuMaxAngle) /*&& (kinkMother->GetMostProbablePID() == AliAODTrack::kKaon)*/ && fRecKinkCutsKaon->IsSelected(kinkMother, AliPID::kKaon)/*fRecKinkCutsKaon->IsGoodKaonKink(kink, esdEvent, kinkMother, kinkDaughter)*/ ) {
		  ++nESDKaonKinks;
			kinkMotherPID = AliPID::kKaon;
		  //ESD kinkMotherSign = kinkMother->GetSign();
	    kinkMotherPDG = (kinkMotherSign > 0 ) ? kKPlus : kKMinus;
		  mother4Momentum.SetVectM(TVector3(kinkMother/*TPCParam*/->Px(), kinkMother/*TPCParam*/->Py(), kinkMother/*TPCParam*/->Pz()), fKaonMass); // Kaon kink hypothesis
		  kinkMotherPt = kinkMother->Pt();//kinkMother->Pt();//mother4Momentum.Pt();
		  kinkMotherY = mother4Momentum.Rapidity();//kinkMother->Y();//mother4Momentum.Rapidity();
		  //if ( !fRecKinkCutsKaon->IsInFiducialKine(mother4Momentum)) continue; // WARNING! Enabling this cut deviates from MSS
		  if ( (kinkMotherPt < fRecKinkCutsKaon->fMinPt) || ( TMath::Abs(kinkMotherY) > fRecKinkCutsKaon->fMaxAbsY) ) continue; // Temporary cut 20130420, deviates from MSS, checked OK 20130620 after fixing kinkMotherPt/Y bug 
		  fhESDKaonKinkDecaysRPhiZ->Fill(rPhiZ);
	    fhESDKaonKinkPtY->Fill(kinkMotherPt, kinkMotherY);
			fhESDKaonKinkQt->Fill(qt);
			fhESDKaonKinkAngle->Fill(kinkAngleDeg);
		  (kinkMotherSign > 0 ) ? fhESDKPlusKinkPtY->Fill(kinkMotherPt, kinkMotherY) : fhESDKMinusKinkPtY->Fill(kinkMotherPt, kinkMotherY);
			(kinkMotherSign > 0 ) ? fhESDKPlusKinkPt->Fill(kinkMotherPt) : fhESDKMinusKinkPt->Fill(kinkMotherPt);
			if ( fCollisionType != kPP ) {
				fhESDKaonKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF);
				(kinkMotherSign > 0 ) ? fhESDKPlusKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF) : fhESDKMinusKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF);
			}
	    fhESDKaonKinkPtEta->Fill(mother4Momentum.Pt(), mother4Momentum.Eta());
		  fhESDMomentumTPCSignalKaonKinks->Fill(kinkMother->P()/*mother4Momentum.P()*/, kinkMother->GetTPCsignal());
		  fhESDKaonKinksMotherVSDaughterTPCncls->Fill(kinkMother->GetTPCNcls(), kinkDaughter->GetTPCNcls());
		  fhESDKaonKinksMotherAndDaughterTPCncls->Fill(kinkMother->GetTPCNcls() + kinkDaughter->GetTPCNcls());
		}
#if 1
		else if ( fRecKinkCutsKaon->IsInQtLimitsPion(qt) && (kinkAngleDeg < maxDecAngpimu) && (kinkAngleDeg > fRecKinkCutsKaon->fMinAnglePion) && (invMassKmu < 0.2) && fRecKinkCutsKaon->IsSelected(kinkMother, AliPID::kPion)/*fRecKinkCutsKaon->IsGoodPionKink(kink, esdEvent, kinkMother, kinkDaughter)*/ ) { // First look for pion kink that has narrow cut and also contaminates kaons
			++nESDPionKinks;
			kinkMotherPID = AliPID::kPion;
		  //ESD kinkMotherSign = kinkMother->GetSign();
			kinkMotherPDG = (kinkMotherSign > 0 ) ? kPiPlus : kPiMinus;
		  mother4Momentum.SetVectM(TVector3(kinkMother/*TPCParam*/->Px(), kinkMother/*TPCParam*/->Py(), kinkMother/*TPCParam*/->Pz()), fPionMass); // Pion kink hypothesis
		  kinkMotherPt = mother4Momentum.Pt();//kinkMother->Pt();//mother4Momentum.Pt();
		  kinkMotherY = mother4Momentum.Rapidity();//kinkMother->Y();//mother4Momentum.Rapidity();
		  if ( !fRecKinkCutsKaon->IsInFiducialKine(mother4Momentum)) continue; // WARNING! Enabling this cut deviates from MSS
		  //if ( (kinkMotherPt < fRecKinkCutsKaon->fMinPt) || ( TMath::Abs(kinkMotherY) > fRecKinkCutsKaon->fMaxAbsY) ) continue; // Temporary cut 20130420, deviates from MSS, checked OK 20130620 after fixing kinkMotherPt/Y bug 
		  fhESDPionKinkDecaysRPhiZ->Fill(rPhiZ);
			fhESDPionKinkPtY->Fill(kinkMotherPt, kinkMotherY);
			fhESDPionKinkQt->Fill(qt);
			fhESDPionKinkAngle->Fill(kinkAngleDeg);
		  (kinkMotherSign > 0 ) ? fhESDPiPlusKinkPtY->Fill(kinkMotherPt, kinkMotherY) : fhESDPiMinusKinkPtY->Fill(kinkMotherPt, kinkMotherY);
			if ( fCollisionType != kPP ) {
				fhESDPionKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF);
				(kinkMotherSign > 0 ) ? fhESDPiPlusKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF) : fhESDPiMinusKinkPtYCent->Fill(kinkMotherPt, kinkMotherY, centralityF);
			}
	    //fhESDPionKinkPtEta->Fill(mother4Momentum.Pt(), mother4Momentum.Eta());
		  fhESDMomentumTPCSignalPionKinks->Fill(kinkMother->P()/*mother4Momentum.P()*/, kinkMother->GetTPCsignal());
		  fhESDPionKinksMotherVSDaughterTPCncls->Fill(kinkMother->GetTPCNcls(), kinkDaughter->GetTPCNcls());
		  fhESDPionKinksMotherAndDaughterTPCncls->Fill(kinkMother->GetTPCNcls() + kinkDaughter->GetTPCNcls());
		}
		else {
			++nESDUnknownKinks;
			kinkMotherPID = AliPID::kUnknown;
			if ( !kinkMother ) continue;
		  //ESD kinkMotherSign = kinkMother->GetSign();
			//kinkMotherPDG = 0;
			fhESDUnknownKinkPtY->Fill(kinkMother->Pt(), kinkMother->Y());
			/*if (mcEvent) {
			 * Used to define proper PIDeff
			}*/
			continue; // Neither Pion or Kaon "good" kink found
		}
#endif
		fhESDAllKinkPtY->Fill(kinkMother->Pt(), kinkMother->Y());
		// End copied from ESD
		fhAODKinkMothervsDaughterPIDcheck->Fill(kinkMother->GetMostProbablePID(), kinkDaughter->GetMostProbablePID());
		if (kinkMother->GetMostProbablePID() == AliAODTrack::kKaon) {
			kinkMotherPID = AliPID::kKaon;
			fhAODKaonKinkPID->Fill(kinkMotherPt/*, kinkMotherY*/);
		}
		else {
			fhAODKaonKinkPIDFake->Fill(kinkMotherPt/*, kinkMotherY*/);
			//continue;
		}
	  if (arrayAODMC) { // Check MC truth
			Int_t motherMClabel = TMath::Abs(kinkMother->GetLabel());
			AliAODMCParticle* kinkMotherMC = static_cast<AliAODMCParticle*>(arrayAODMC->At(motherMClabel));
			if ( kinkMotherMC ) {
				CheckMCtruth(kinkMother, mother4Momentum, kinkMotherPID, kinkMotherMC, kinkMotherMC->IsPhysicalPrimary());
			} // kinkMotherMC
	  }
		if (!fEsdTrackCutsPartner) continue; // Skip resonance analysis if no ESD track cuts provided for partner
		KinkResonanceAnalysis(aodEvent, kinkMother, mother4Momentum, kinkMotherPID, 0x0, 0x0);
} // This part crashes with AliMultiInputEventHandler
    //delete kinkMotherAtKinkVtx;
	} // AOD vertex-kink loop
	fhESDnKinks->Fill(nKinkVertices);
	} // AOD analysis
	return 0;
}

//________________________________________________________________________
Float_t AliAnalysisTaskKinksFilimon::GetRecCentrality(AliVEvent* const inputEvent, const char* centEstimator) const {

	Float_t centralityF = -1;
	if (inputEvent) {
    AliCentrality* recCentrality = ( fCollisionType != kPP ) ? inputEvent->GetCentrality() : 0;
    recCentrality ? centralityF = recCentrality->GetCentralityPercentile(centEstimator) : -1;
	}
	return centralityF;
}

//________________________________________________________________________
Int_t AliAnalysisTaskKinksFilimon::EventSelection(AliVEvent* const inputEvent, AliESDEvent* const esdEvent, const AliAODEvent* const aodEvent, const Float_t centralityF/*Option_t **/) {

	if (!inputEvent) return(-1);
  // Centrality?
  //AliEventplane* recEventPlane = inputEvent->GetEventplane();
	// Multiplicity
	const Int_t recNtracksUnbiased = inputEvent->GetNumberOfTracks();
	fhRecMultUnbiased->Fill(recNtracksUnbiased);
	if ( fCollisionType != kPP ) fhRecMultCentralityUnbiased->Fill(recNtracksUnbiased, centralityF);
	fRecEventCuts->GetHist()->Fill(AliRecEventCuts::kProcessedEvents);
	// Physics selection
	//AliInputEventHandler* fInputHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if ( esdEvent && !( fMainInputHandler && (fMainInputHandler->IsEventSelected()&fOfflineTriggerType)) ) return -1; // Physics selection offline (WARNING ESD only!
	fRecEventCuts->GetHist()->Fill(AliRecEventCuts::kPhysSelEvents);
	//hMultPS->Fill(NTracks); // Data Multiplicity after Physics selection
	const Int_t recNtracks = inputEvent->GetNumberOfTracks();
	fhRecMult->Fill(recNtracks);
	if ( fCollisionType != kPP ) fhRecMultCentrality->Fill(recNtracks, centralityF);
  // Primary vertex
	const AliVVertex* recPrimaryVertex = fIsAOD ? GetPrimaryVertex(aodEvent) : GetPrimaryVertex(esdEvent); // Rec primary vertex
	if ( !(recPrimaryVertex && GetPrimaryVertexQuality(recPrimaryVertex, 1, 1e10) ) ) return(-1);
	Double_t recPrimaryVertexXYZArray[3] = { 0, 0, 0 };
	recPrimaryVertex->GetXYZ(recPrimaryVertexXYZArray);
	TVector3 recPrimaryVertexXYZ(recPrimaryVertexXYZArray);
	Double_t recPrimaryVertexRPhiZarray[3] = {recPrimaryVertexXYZ.Perp(), recPrimaryVertexXYZ.Phi(), recPrimaryVertexXYZ.Z()};
	fhRecPrimaryVertexRPhiZ->Fill(recPrimaryVertexRPhiZarray);
	if (TMath::Abs(recPrimaryVertexXYZ.Z()) > 10.) return -1; // Temporarily cut here all events even with |z-vtx| > 10cm
	fRecEventCuts->GetHist()->Fill(AliRecEventCuts::kVtxRange);
	/*hvtx->Fill(recPrimaryVertexXYZArray[0], recPrimaryVertexXYZArray[1], recPrimaryVertexXYZArray[2]); //ESD primary vertex position (x-y-z)
	hvtxy->Fill(recPrimaryVertexXYZArray[0], recPrimaryVertexXYZArray[1]); 
	hvtyz->Fill(recPrimaryVertexXYZArray[1], recPrimaryVertexXYZArray[2]); 
	hvtxz->Fill(recPrimaryVertexXYZArray[0], recPrimaryVertexXYZArray[2]); 
	if (TMath::Abs(recPrimaryVertexXYZArray[2])>10.) return;*/
	//hMultPSV->Fill(NTracks); // Data Multiplicity after vertex cut
	
	// AliAnalysisTaskKinksFilimon
	/*
	 * Internally use an AliKineTrackCuts (concretely AliCFParticleGenCuts, AliCFTrackKineCuts ) for MC!!!
	 * and one for main vertex cuts
	Bool_t IsPrimary(const AliESDtrack* esdTrack) { return fPrimaryTrackCuts->IsSelected(esdTack); };
	SetAbsKinkEta, SetAbsKinkY, SetMinPt, SetNSigmaXXX, SetMinQt, SetMinkKinkR, SetMaxKinkR;
	*/
	/*AliCFParticleGenCuts* genPrimaryChargedTrackCut=new AliESDtrackCuts("genPrimaryChargedTrackCut", "genPrimaryChargedTrackCut");
	AliCFParticleGenCuts* genFiducialVolumeKinkDecayCut=new AliESDtrackCuts("genFiducialVolumeKinkDecayCut", "genFiducialVolumeKinkDecayCut");

	delete genPrimaryChargedTrackCut;
	delete genFiducialVolumeKinkDecayCut;*/

	return(0);
}

Int_t AliAnalysisTaskKinksFilimon::KinkResonanceAnalysis(const AliVEvent* const recEvent, const AliVParticle* const kinkMother, const TLorentzVector& mother4Momentum, const AliPID::EParticleType kinkMotherPID, const AliVEvent* const mixEvent, AliMCEvent* const mcEvent) const { // TODO: Unify ESD/AOD 
		//const AliESDEvent* esdEvent = 0x0;
		//Double_t kinkMotherSign = 0;
		//const AliAODEvent* aodEvent = 0x0;
		Short_t kinkMotherCharge = kinkMother->Charge();
#if 0
	  if (!fIsAOD) { // ESD
      //esdEvent = dynamic_cast<const AliESDEvent*>(mixEvent ? mixEvent : recEvent); // If mixEvent exists, use it to generate background only
	    //if (!esdEvent) return(-1);
      const AliESDtrack* esdKinkMother = dynamic_cast<const AliESDtrack*>(kinkMother);
	    if (!esdKinkMother) return(-1);
		  kinkMotherSign = esdKinkMother->GetSign();
		}
		else { // AOD
      //aodEvent = dynamic_cast<const AliAODEvent*>(mixEvent ? mixEvent : recEvent); // If mixEvent exists, use it to generate background only
	    //if (!aodEvent) return(-1);
      const AliAODTrack* aodKinkMother = dynamic_cast<const AliAODTrack*>(kinkMother);
	    if (!aodKinkMother) return(-1);
		  kinkMotherCharge = aodKinkMother->Charge(); 
		}
#endif
		/*AliExternalTrackParam* partnerTrackTPCParam = (AliExternalTrackParam *)partnerTrack->GetTPCInnerParam();
    if (!partnerTrackTPCParam) continue;*/ // Remove
		TVector3 partnerTrack3Momentum;
	  TLorentzVector resonance4Momentum, partnerTrack4MomentumKaonHypothesis, partnerTrack4MomentumPionHypothesis, partnerTrack4MomentumProtonHypothesis;
	  Int_t recNtracks = recEvent->GetNumberOfTracks();
		//if(!IsKink(esd, partnerTrack->GetKinkIndex(0), partnerTrack3Momentum)) continue; // Remove
		Bool_t isLikeSign = kFALSE;
	  const Double_t resRecabsYmax = 1.5;
		for (Int_t iPartnerTrack = 0; iPartnerTrack < recNtracks; ++iPartnerTrack) { // Look for phi,K*,Lambda->Kink+partner
			//if ( (iPartnerTrack == TMath::Abs(motherIndex)) || (iPartnerTrack == TMath::Abs(daughterIndex)) ) continue; // Remove
			AliVParticle* partnerTrack = recEvent->GetTrack(iPartnerTrack);
			if ( kinkMother == partnerTrack) continue;
#if 1			
			if (!fIsAOD) { // ESD
			  AliESDtrack* partnerTrackESD = static_cast<AliESDtrack*>(partnerTrack);//esdEvent->GetTrack(iPartnerTrack);
		    if ( !fEsdTrackCutsPartner->IsSelected(partnerTrackESD) ) continue;
				//isLikeSign = (kinkMotherSign * partnerTrackESD->GetSign() > 0) ? kTRUE : kFALSE;
			}
			else { // AOD
			  AliAODTrack* partnerTrackAOD = static_cast<AliAODTrack*>(partnerTrack);//aodEvent->GetTrack(iPartnerTrack);
		    if (!PassesESDTrackCuts(partnerTrackAOD, fEsdTrackCutsPartner)) continue;
		    //if ( !partnerTrackAOD->IsPrimaryCandidate() ) continue;
				//isLikeSign = (kinkMotherCharge * partnerTrackAOD->Charge() > 0) ? kTRUE : kFALSE;
			}
#endif
		  isLikeSign = (kinkMotherCharge * partnerTrack->Charge() > 0) ? kTRUE : kFALSE;
			partnerTrack3Momentum.SetXYZ(partnerTrack/*TPCParam*/->Px(), partnerTrack/*TPCParam*/->Py(), partnerTrack/*TPCParam*/->Pz());
			partnerTrack4MomentumKaonHypothesis.SetVectM(partnerTrack3Momentum, fKaonMass);
			partnerTrack4MomentumPionHypothesis.SetVectM(partnerTrack3Momentum, fPionMass);
			partnerTrack4MomentumProtonHypothesis.SetVectM(partnerTrack3Momentum, fProtonMass);
			// Note 20130328 Need to replace the nSigma tests with max_probability(kaon, pion, proton);
			if ( !isLikeSign && !mixEvent) { // Opposite sign - good resonance candidate
			  resonance4Momentum = mother4Momentum + partnerTrack4MomentumKaonHypothesis;
			  /*if ( (kinkMotherPID == AliPID::kPion) &&
					   (mother4Momentum.Pt() < 1.5) &&
					   (fRecKinkCutsKaon->IsInFiducialKine(resonance4Momentum)) && 
					   //(resonance4Momentum.Pt() > fRecKinkCutsKaon->fMinPt) &&
					   //(TMath::Abs(resonance4Momentum.Y()) < fRecKinkCutsKaon->fMaxAbsY) && 
					   ( IsSelected(partnerTrack, AliPID::kKaon) ) ) { // Kstar
			    fhRecKstarInvMassPt->Fill(resonance4Momentum.M(), resonance4Momentum.Pt());
				  if ( (resonance4Momentum.M() > fKstarMassRangeMin) && (resonance4Momentum.M() < fKstarMassRangeMax) ) fhRecKstarPtYtest->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
					continue;
				}*/
			  if ( (kinkMotherPID == AliPID::kKaon) &&
					   IsInFiducialKine(resonance4Momentum, fRecKinkCutsKaon->fMinPt, resRecabsYmax) &&
					   //(fRecKinkCutsKaon->IsInFiducialKine(resonance4Momentum)) && 
					   //(resonance4Momentum.Pt() > fRecKinkCutsKaon->fMinPt) && 
					   //(TMath::Abs(resonance4Momentum.Y()) < fRecKinkCutsKaon->fMaxAbsY) && 
					   ( /*partnerTrack->GetTPCsignalSigma()*/ IsSelected(partnerTrack, AliPID::kKaon) ) &&
					   ( TMath::Cos(mother4Momentum.Vect().Angle(partnerTrack4MomentumKaonHypothesis.Vect())) > 0.4 * mother4Momentum.Pt() - 1.4 ) ) { // Phi
			    fhRecPhiInvMassPt->Fill(resonance4Momentum.M(), resonance4Momentum.Pt());
				  if ( (resonance4Momentum.M() > fPhiMassRangeMin) && (resonance4Momentum.M() < fPhiMassRangeMax) ) {
						fhRecPhiPtYtest->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
						fhRecPhiInvMassPtYtest->Fill(resonance4Momentum.M(), resonance4Momentum.Pt(), resonance4Momentum.Y());
					}
				}
				resonance4Momentum = mother4Momentum + partnerTrack4MomentumPionHypothesis;
			  if ( (kinkMotherPID == AliPID::kKaon) && 
					   IsInFiducialKine(resonance4Momentum, fRecKinkCutsKaon->fMinPt, resRecabsYmax) &&
					   //(fRecKinkCutsKaon->IsInFiducialKine(resonance4Momentum)) && 
					   //(resonance4Momentum.Pt() > fRecKinkCutsKaon->fMinPt) && 
					   //(TMath::Abs(resonance4Momentum.Y()) < fRecKinkCutsKaon->fMaxAbsY) && 
					   ( IsSelected(partnerTrack, AliPID::kPion) ) &&
					   ( TMath::Cos(mother4Momentum.Vect().Angle(partnerTrack4MomentumPionHypothesis.Vect())) > 0.4 * mother4Momentum.Pt() - 1.4 ) ) { // Kstar
			    fhRecKstarInvMassPt->Fill(resonance4Momentum.M(), resonance4Momentum.Pt());
				  if ( (resonance4Momentum.M() > fKstarMassRangeMin) && (resonance4Momentum.M() < fKstarMassRangeMax) ) fhRecKstarPtYtest->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
				}
				resonance4Momentum = mother4Momentum + partnerTrack4MomentumProtonHypothesis;
			  if ( (kinkMotherPID == AliPID::kKaon) && 
					   IsInFiducialKine(resonance4Momentum, fRecKinkCutsKaon->fMinPt, resRecabsYmax) &&
					   //(fRecKinkCutsKaon->IsInFiducialKine(resonance4Momentum)) && 
					   //(resonance4Momentum.Pt() > fRecKinkCutsKaon->fMinPt) && 
					   //(TMath::Abs(resonance4Momentum.Y()) < fRecKinkCutsKaon->fMaxAbsY) && 
					   ( IsSelected(partnerTrack, AliPID::kProton) ) &&
					   ( TMath::Cos(mother4Momentum.Vect().Angle(partnerTrack4MomentumProtonHypothesis.Vect())) > 0.333 * mother4Momentum.Pt() - 1.0 ) ) { // Lambda
			    fhRecLambdaInvMassPt->Fill(resonance4Momentum.M(), resonance4Momentum.Pt());
				  if ( (resonance4Momentum.M() > fLambdaMassRangeMin) && (resonance4Momentum.M() < fLambdaMassRangeMax) ) fhRecLambdaPtYtest->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
				}
				if (mcEvent) { // Check MC truth
			    Int_t motherMClabel = TMath::Abs(kinkMother->GetLabel());
			    AliMCParticle* kinkMotherMC = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(motherMClabel));
			    if ( kinkMotherMC ) {
			      Int_t kaonMotherIdx = kinkMotherMC->GetMother();
			      if ( (kaonMotherIdx != -1 ) && (kaonMotherIdx < mcEvent->GetNumberOfPrimaries()) ) {
				    AliVParticle* kaonMother = mcEvent->GetTrack(kaonMotherIdx);
				    Int_t kaonMotherPdg = kaonMother->PdgCode();
				    //fhMCkaonMotherPdg->Fill(kaonMotherPdg);
				    switch ( TMath::Abs(kaonMotherPdg) ) {
						  case kPhi:
							  fhESDPhiPtYtestTrueMC->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y()); // FIXME: Put same conditions as above!!!
							  break;
						  case kKstar0:
							  fhESDKstarPtYtestTrueMC->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
							  break;
#if 0
						  case kLambda1520:
							  fhESDPhiPtYtestTrueMC->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
							  break;
#endif
						  //default:
							  //fhMCFakeResonance2KaonPtYFakeMC->Fill(kinkMotherPt, kinkMotherY);
				    }
			     } // kaonMotherIdx good 
			    } // kinkMotherMC
				}
			} // Opposite sign analysis end
			else { // Likesign - background
			  resonance4Momentum = mother4Momentum + partnerTrack4MomentumKaonHypothesis;
			  /*if ( (kinkMotherPID == AliPID::kPion) &&
					   (mother4Momentum.Pt() < 1.5) &&
					   (fRecKinkCutsKaon->IsInFiducialKine(resonance4Momentum)) && 
					   //(resonance4Momentum.Pt() > fRecKinkCutsKaon->fMinPt) && 
					   //(TMath::Abs(resonance4Momentum.Y()) < fRecKinkCutsKaon->fMaxAbsY) && 
					   ( IsSelected(partnerTrack, AliPID::kKaon) ) ) { // Kstar
			    fhRecKstarInvMassPtLikeSign->Fill(resonance4Momentum.M(), resonance4Momentum.Pt());
				  if ( (resonance4Momentum.M() > fKstarMassRangeMin) && (resonance4Momentum.M() < fKstarMassRangeMax) ) fhRecKstarPtYtestLikeSign->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
					continue;
				}*/
			  if ( (kinkMotherPID == AliPID::kKaon) && 
					   IsInFiducialKine(resonance4Momentum, fRecKinkCutsKaon->fMinPt, resRecabsYmax) &&
					   //(fRecKinkCutsKaon->IsInFiducialKine(resonance4Momentum)) && 
					   //(resonance4Momentum.Pt() > fRecKinkCutsKaon->fMinPt) && 
					   //(TMath::Abs(resonance4Momentum.Y()) < fRecKinkCutsKaon->fMaxAbsY) && 
					   ( IsSelected(partnerTrack, AliPID::kKaon) ) &&
					   ( TMath::Cos(mother4Momentum.Vect().Angle(partnerTrack4MomentumKaonHypothesis.Vect())) > 0.4 * mother4Momentum.Pt() - 1.4 ) ) { // Phi
			    fhRecPhiInvMassPtLikeSign->Fill(resonance4Momentum.M(), resonance4Momentum.Pt());
				  if ( (resonance4Momentum.M() > fPhiMassRangeMin) && (resonance4Momentum.M() < fPhiMassRangeMax) ) {
						fhRecPhiPtYtestLikeSign->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
						fhRecPhiInvMassPtYtestLikeSign->Fill(resonance4Momentum.M(), resonance4Momentum.Pt(), resonance4Momentum.Y());
					}
				}
				resonance4Momentum = mother4Momentum + partnerTrack4MomentumPionHypothesis;
			  if ( (kinkMotherPID == AliPID::kKaon) && 
					   IsInFiducialKine(resonance4Momentum, fRecKinkCutsKaon->fMinPt, resRecabsYmax) &&
					   //(fRecKinkCutsKaon->IsInFiducialKine(resonance4Momentum)) && 
					   //(resonance4Momentum.Pt() > fRecKinkCutsKaon->fMinPt) && 
					   //(TMath::Abs(resonance4Momentum.Y()) < fRecKinkCutsKaon->fMaxAbsY) && 
					   ( IsSelected(partnerTrack, AliPID::kPion) ) &&
					   ( TMath::Cos(mother4Momentum.Vect().Angle(partnerTrack4MomentumPionHypothesis.Vect())) > 0.4 * mother4Momentum.Pt() - 1.4 ) ) { // Kstar
			    fhRecKstarInvMassPtLikeSign->Fill(resonance4Momentum.M(), resonance4Momentum.Pt());
				  if ( (resonance4Momentum.M() > fKstarMassRangeMin) && (resonance4Momentum.M() < fKstarMassRangeMax) ) fhRecKstarPtYtestLikeSign->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
				}
				resonance4Momentum = mother4Momentum + partnerTrack4MomentumProtonHypothesis;
			  if ( (kinkMotherPID == AliPID::kKaon) && 
					   IsInFiducialKine(resonance4Momentum, fRecKinkCutsKaon->fMinPt, resRecabsYmax) &&
					   //(fRecKinkCutsKaon->IsInFiducialKine(resonance4Momentum)) && 
					   //(resonance4Momentum.Pt() > fRecKinkCutsKaon->fMinPt) && 
					   //(TMath::Abs(resonance4Momentum.Y()) < fRecKinkCutsKaon->fMaxAbsY) && 
					   ( IsSelected(partnerTrack, AliPID::kProton) ) &&
					   ( TMath::Cos(mother4Momentum.Vect().Angle(partnerTrack4MomentumProtonHypothesis.Vect())) > 0.333 * mother4Momentum.Pt() - 1.0 ) ) { // Lambda
			    fhRecLambdaInvMassPtLikeSign->Fill(resonance4Momentum.M(), resonance4Momentum.Pt());
				  if ( (resonance4Momentum.M() > fLambdaMassRangeMin) && (resonance4Momentum.M() < fLambdaMassRangeMax) ) fhRecLambdaPtYtestLikeSign->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
				}
				if (mcEvent) { // Check MC truth
			    Int_t motherMClabel = TMath::Abs(kinkMother->GetLabel());
			    AliMCParticle* kinkMotherMC = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(motherMClabel));
			    if ( kinkMotherMC ) {
			      Int_t kaonMotherIdx = kinkMotherMC->GetMother();
			      if ( (kaonMotherIdx != -1 ) && (kaonMotherIdx < mcEvent->GetNumberOfPrimaries()) ) {
				    AliVParticle* kaonMother = mcEvent->GetTrack(kaonMotherIdx);
				    Int_t kaonMotherPdg = kaonMother->PdgCode();
				    //fhMCkaonMotherPdg->Fill(kaonMotherPdg);
				    switch ( TMath::Abs(kaonMotherPdg) ) {
						  case kPhi:
							  fhESDPhiPtYtestLikeSignTrueMC->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
							  break;
						  case kKstar0:
							  fhESDKstarPtYtestLikeSignTrueMC->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
							  break;
#if 0
						  case kLambda1520:
							  fhESDPhiPtYtestTrueMC->Fill(resonance4Momentum.Pt(), resonance4Momentum.Y());
							  break;
#endif
						  //default:
							  //fhMCFakeResonance2KaonPtYFakeMC->Fill(kinkMotherPt, kinkMotherY);
				    }
			     } // kaonMotherIdx good 
			    } // kinkMotherMC
				}
			} // Likesign - background end
		} // Rec resonance partner track loop
		return(0);
}

//________________________________________________________________________
AliVEvent* AliAnalysisTaskKinksFilimon::GetMixedEvent() const {
	
  AliMixInputEventHandler* mixInputHandler = dynamic_cast<AliMixInputEventHandler*>(fMultiInputHandler->GetFirstMultiInputHandler());
  if (!mixInputHandler) return 0x0;
  if (mixInputHandler->CurrentBinIndex() < 0) {
    AliDebug(AliLog::kDebug + 1, "Current event mixInputHandler->CurrentEntry() == -1");
    return 0x0;
  }
  AliDebug(AliLog::kDebug, Form("Mixing %lld %d [%lld,%lld] %d", mixInputHandler->CurrentEntry(), mixInputHandler->CurrentBinIndex(), mixInputHandler->CurrentEntryMain(), mixInputHandler->CurrentEntryMix(), mixInputHandler->NumberMixed()));
  AliMultiInputEventHandler* inEvHMixedCurrent = mixInputHandler->GetFirstMultiInputHandler(); // for buffer = 1
  AliInputEventHandler* ihMixedCurrent = inEvHMixedCurrent ? inEvHMixedCurrent->GetFirstInputEventHandler() : 0x0;
  return ihMixedCurrent ? ihMixedCurrent->GetEvent() : 0x0;
} 

//________________________________________________________________________
Int_t AliAnalysisTaskKinksFilimon::BookCentralityHist(const char* name, const char* title, const TH2* const histPtYtemplate, const Int_t nCentBins, const Float_t* centBins, TObjArray* const histArray, TCollection* const outputCont) {
       TString titleBase(title);
       titleBase = titleBase(0, titleBase.First(';'));
       TString titleAxes(title);
       titleAxes = titleAxes(titleAxes.First(';'), titleAxes.Length());
               for ( Int_t iCent = 0; iCent < nCentBins; ++iCent ) {
                       Float_t centMin = centBins ? centBins[iCent] : TMath::Nint(static_cast<Float_t>(iCent)/nCentBins*100);
                       Float_t centMax = centBins ? centBins[iCent+1] : TMath::Nint(static_cast<Float_t>(iCent+1)/nCentBins*100);
                       TString histName(Form("%s_cent%2.1f-%2.1f", name, centMin, centMax));
                       TString histTitle(Form("%s, Centrality %2.1f-%2.1f%%;%s", titleBase.Data(), centMin, centMax, titleAxes.Data()));
           TH2* tmp = dynamic_cast<TH2*>(histPtYtemplate->Clone(histName));
           tmp->SetTitle(histTitle);
                       histArray->Add(tmp);
                       outputCont->Add(tmp);
               }
               return(0);
 }
 
//________________________________________________________________________
Int_t AliAnalysisTaskKinksFilimon::CheckMCtruth(AliVTrack* const /*kinkMotherRec*/, const TLorentzVector& mother4Momentum, const AliPID::EParticleType kinkMotherPID, AliVParticle* const kinkMotherMC, const Bool_t isPhysicalPrimary) const {

		  Double_t kinkMotherPt = mother4Momentum.Pt();//kinkMother->Pt();//mother4Momentum.Pt();
		  Double_t kinkMotherY = mother4Momentum.Rapidity();//kinkMother->Y();//mother4Momentum.Rapidity();
				Int_t kinkMotherMCPdg = kinkMotherMC->PdgCode();
				//TParticle* kinkMotherMCparticle = kinkMotherMC->Particle();
				switch (kinkMotherPID) {
					case AliPID::kKaon: {
				  if ( TMath::Abs(kinkMotherMCPdg) == kKPlus ) {
						if ( isPhysicalPrimary ) {
							fhESDKaonKinkPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
							kinkMotherMCPdg == kKPlus ? fhESDKPlusKinkPtYTrueMC->Fill(kinkMotherPt, kinkMotherY) : fhESDKMinusKinkPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
						}
						else {
							fhESDKaonKinkPtYTrueMCsecondary->Fill(kinkMotherPt, kinkMotherY);
							kinkMotherMCPdg == kKPlus ? fhESDKPlusKinkPtYTrueMCsecondary->Fill(kinkMotherPt, kinkMotherY) : fhESDKMinusKinkPtYTrueMCsecondary->Fill(kinkMotherPt, kinkMotherY);
						}
					}
			    else {
						fhESDKaonKinkPtYFakeMC->Fill(kinkMotherPt, kinkMotherY);
						kinkMotherMCPdg == kKPlus ? fhESDKPlusKinkPtYFakeMC->Fill(kinkMotherPt, kinkMotherY) : fhESDKMinusKinkPtYFakeMC->Fill(kinkMotherPt, kinkMotherY);
					}
				break;
				}
					case AliPID::kPion: {
				  if ( TMath::Abs(kinkMotherMCPdg) == kPiPlus ) {
						if ( isPhysicalPrimary ) {
							fhESDPionKinkPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
							kinkMotherMCPdg == kPiPlus ? fhESDPiPlusKinkPtYTrueMC->Fill(kinkMotherPt, kinkMotherY) : fhESDPiMinusKinkPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
						}
						else {
							fhESDPionKinkPtYTrueMCsecondary->Fill(kinkMotherPt, kinkMotherY);
							kinkMotherMCPdg == kPiPlus ? fhESDPiPlusKinkPtYTrueMCsecondary->Fill(kinkMotherPt, kinkMotherY) : fhESDPiMinusKinkPtYTrueMCsecondary->Fill(kinkMotherPt, kinkMotherY);
						}
					}
			    else {
						fhESDPionKinkPtYFakeMC->Fill(kinkMotherPt, kinkMotherY);
						kinkMotherMCPdg == kPiPlus ? fhESDPiPlusKinkPtYFakeMC->Fill(kinkMotherPt, kinkMotherY) : fhESDPiMinusKinkPtYFakeMC->Fill(kinkMotherPt, kinkMotherY);
					}
				break;
				}
					default: break;
				};
#if 0
			  // Study Kaon mother
			  Int_t kaonMotherIdx = kinkMotherMC->GetMother();
			  if ( (kaonMotherIdx != -1 ) && (kaonMotherIdx < mcEvent->GetNumberOfPrimaries()) ) {
				  AliVParticle* kaonMother = mcEvent->GetTrack(kaonMotherIdx);
				  Int_t kaonMotherPdg = kaonMother->PdgCode();
				  //fhMCkaonMotherPdg->Fill(kaonMotherPdg);
				  switch ( TMath::Abs(kaonMotherPdg) ) {
						case kPhi:
							fhESDPhi2KaonPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
							break;
						case kKstar0:
							if ( kinkMotherPID == AliPID::kKaon ) fhESDKstar2KaonPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
							break;
						case kLambda1520:
							fhESDLambda2KaonPtYTrueMC->Fill(kinkMotherPt, kinkMotherY);
							break;
						//default:
							//fhMCFakeResonance2KaonPtYFakeMC->Fill(kinkMotherPt, kinkMotherY);
				  }
			  }
#endif	
  return(0);
}

//________________________________________________________________________
void AliAnalysisTaskKinksFilimon::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputCont = dynamic_cast<TCollection*> (GetOutputData(1));
  if (!fOutputCont) {
    Printf("ERROR: Output list not available");
    return;
  }

}
