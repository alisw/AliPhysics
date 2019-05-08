#ifndef AliAnalysisTaskKinksFilimon_cxx
#define AliAnalysisTaskKinksFilimon_cxx

// Kink and resonance analysis task
// Author: Filimon Roukoutakis, University of Athens
// Uses AliESDkinkCuts helper class to identify K and pi from their kink decays to mu. These are fed to a function that searches for K0*(892) and phi(1020) resonances by combining with "partner" tracks (unlike and like sign). Works on MC, ESD and AOD (AOD currently limited).

#include <AliAnalysisTaskSE.h>
#include <AliVEvent.h>
#include <AliInputEventHandler.h>
#include <AliMultiInputEventHandler.h>
#include <AliPID.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TVector3.h>
#include <AliAnalysisCuts.h>
#include <AliTrackReference.h>
#include <AliMCParticle.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliESDpidCuts.h>
#include <AliESDkink.h>
#include <AliPIDResponse.h>
#include "AliESDkinkCuts.h"
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliESDtrackCuts.h>

//class TH1F;
class TH3;
class AliESDtrackCuts;
class AliCFParticleGenCuts;
class AliKineTrackCuts;
class AliVVertex;
class AliESDtrack;
//class TF1;
//class AliESDkink;
//class THnSparseF;
//class AliESDEvent;
//class AliPIDResponse;
//class AliESDpidCuts;

class AliAnalysisTaskKinksFilimon : public AliAnalysisTaskSE {
 public:
  enum EResonanceType { kPhi=333, kKstar0=313, kLambda1520=3124 };
  enum ECollisionType { kPP, kPbPb, kPPb, kAA };
  AliAnalysisTaskKinksFilimon() : AliAnalysisTaskSE(), fOutputContClass(0), fOutputCont(0), fMainInputHandler(0), fIsAOD(kFALSE), fUseMC(kFALSE), fFillCutHist(kFALSE), fEventSelectionCutsMC(kFALSE), fCollisionType(kPP), fOfflineTriggerType(AliVEvent::kMB), fPartIdDedx(AliPID::kElectron), fPIDResponse(0), fRecEventCuts(0), fRecCentCuts(0), fCommonKineTrackCuts(0), fEsdTrackCutsKinkMother(0), fEsdTrackCutsKinkDaughter(0), fEsdTrackCutsPartner(0), fEsdPIDCutsArray(), fRecKinkCutsKaon(0), fAODFilterMap(0), /*fRecKinkCutsPion(0),*/ /*fMcKinkCutsPion(0),*/ fHistPtYTemplate(0), fHistRPhiZTemplate(0), fKaonMass(0), fMuonMass(0), fPionMass(0), fProtonMass(0), fElectronMass(0), fPhiMassRangeMin(0), fPhiMassRangeMax(0), fKstarMassRangeMin(0), fKstarMassRangeMax(0), fLambdaMassRangeMin(0), fLambdaMassRangeMax(0),
  fhRecMultUnbiased(0), fhRecMult(0), fhRecMultCentralityUnbiased(0), fhRecMultCentrality(0), fhRecAllPtY(0), fhESDnKinks(0), fhESDKinkQt(0), fhESDKinkAngle(0), fhESDKinkPAngle(0), fhESDKinkDCA(0), fhESDAllKinkPtY(0), fhESDUnknownKinkPtY(0), fhESDKaonKinkPtY(0), fhESDPionKinkPtY(0), fhESDKaonKinkPtEta(0), fhESDKPlusKinkPtY(0), fhESDKMinusKinkPtY(0), fhESDKPlusKinkPt(0), fhESDKMinusKinkPt(0), fhESDPiPlusKinkPtY(0), fhESDPiMinusKinkPtY(0), fhESDKaonKinkPtYCent(0), fhESDPionKinkPtYCent(0), fhESDKPlusKinkPtYCent(0), fhESDPiPlusKinkPtYCent(0), fhESDKMinusKinkPtYCent(0), fhESDPiMinusKinkPtYCent(0), fhRecPrimaryVertexRPhiZ(0), fhESDAllKinkDecaysRPhiZ(0), fhESDKaonKinkDecaysRPhiZ(0), fhESDPionKinkDecaysRPhiZ(0), fhESDKaonKinkQt(0), fhESDKaonKinkAngle(0), fhESDPionKinkQt(0), fhESDPionKinkAngle(0), fhESDMomentumTPCSignalKaonKinks(0), fhESDKaonKinksMotherAndDaughterTPCncls(0), fhESDKaonKinksMotherVSDaughterTPCncls(0), fhESDMomentumTPCSignalPionKinks(0), fhESDPionKinksMotherAndDaughterTPCncls(0), fhESDPionKinksMotherVSDaughterTPCncls(0), fhESDKaonOverPionEbE(0), fhESDInvMassKinkDecayMuNu(0),
  fhMCmainVertexRPhiZ(0), fhMCKaonKinkDecaysRPhiZ(0), fhMCmult(0), fhMCmultPrim(0), fhMCPtYall(0), fhMCPtYprim(0), fhMCpdg(0), fhMCPtYprimKaon(0), fhMCPtYprimKPlus(0), fhMCPtYprimKMinus(0), fhMCPtprimKPlus(0), fhMCPtprimKMinus(0), fhMCPtYCentprimKaon(0), fhMCPtYCentprimKPlus(0), fhMCPtYCentprimKMinus(0), fhMCPtYkinkKaonFiducial(0), fhMCPtYkinkKPlusFiducial(0), fhMCPtYkinkKMinusFiducial(0), fhMCKinkKaonQt(0), fhMCKinkKaonAngle(0), fhMCKaonKinkDecaysTrackLength(0), fhMCKaonKinkDecaysTrackTime(0), fhMCkaonMotherPdg(0), fhMCpionMotherPdg(0), fhMCPhiPtY(0), fhMCPtYprimPhi(0), fhMCPhiInvMassPtYtest(0), fhMCInvMassPtprimPhi(0), fhMCInvMassPtCentprimPhi(0), fhMCKstarPtY(0), fhMCPtYprimKstar(0), fhMCInvMassPtprimKstar(0), fhMCInvMassPtCentprimKstar(0), fhMCLambdaPtY(0), fhMCPtYprimLambda(0), fhMCInvMassPtprimLambda(0), fhMCInvMassPtCentprimLambda(0), fhMCPhi2KaonPtY(0), fhMCKstar2KaonPtY(0), fhMCKaon2MuonPtY(0), fhMCKaon2PionPtY(0), fhMCPtYprimPion(0), fhMCPtYprimPiPlus(0), fhMCPtYprimPiMinus(0), fhMCPtYCentprimPion(0), fhMCPtYCentprimPiPlus(0), fhMCPtYCentprimPiMinus(0), fhMCPionKinkDecaysRPhiZ(0), fhMCPion2MuonPtY(0), fhMCPtYkinkPionFiducial(0), fhMCPtYkinkPiPlusFiducial(0), fhMCPtYkinkPiMinusFiducial(0), fhMCKinkPionQt(0), fhMCKinkPionAngle(0), fhMCMotherDaughterPdg(0), fhMCKaonOverPionEbE(0), fhMCPhiDecayOpeningAngle(0), fhMCKstarDecayOpeningAngle(0), fhMCLambdaDecayOpeningAngle(0), 
  fhESDKinkQtTrueMC(0), fhESDKinkAngleTrueMC(0), fhESDKaonKinkPtYTrueMC(0), fhESDKaonKinkPtYFakeMC(0), fhESDKaonKinkPtYTrueMCsecondary(0), fhESDKPlusKinkPtYTrueMC(0), fhESDKPlusKinkPtYFakeMC(0), fhESDKPlusKinkPtYTrueMCsecondary(0), fhESDKMinusKinkPtYTrueMC(0), fhESDKMinusKinkPtYFakeMC(0), fhESDKMinusKinkPtYTrueMCsecondary(0), fhESDPionKinkPtYTrueMC(0), fhESDPionKinkPtYFakeMC(0), fhESDPionKinkPtYTrueMCsecondary(0), fhESDPiPlusKinkPtYTrueMC(0), fhESDPiPlusKinkPtYFakeMC(0), fhESDPiPlusKinkPtYTrueMCsecondary(0), fhESDPiMinusKinkPtYTrueMC(0), fhESDPiMinusKinkPtYFakeMC(0), fhESDPiMinusKinkPtYTrueMCsecondary(0), fhESDPhi2KaonPtYTrueMC(0), fhESDPhiPtYtestTrueMC(0), fhESDPhiPtYtestFakeMC(0), fhESDPhiPtYtestLikeSignTrueMC(0), fhESDPhiPtYtestLikeSignFakeMC(0), fhESDKstar2KaonPtYTrueMC(0), fhESDKstarPtYtestTrueMC(0), fhESDKstarPtYtestFakeMC(0), fhESDKstarPtYtestLikeSignTrueMC(0), fhESDKstarPtYtestLikeSignFakeMC(0), fhESDLambda2KaonPtYTrueMC(0),
  fhAODnVertices(0), fhAODnKinkVertices(0), fhAODnKinkDaughters(0), fhAODAllKinkDecaysRPhiZ(0), fhAODKinkKaonQt(0), fhAODKinkKaonAngle(0), fhAODKinkMothervsDaughterPIDcheck(0), fhAODKaonKinkPID(0), fhAODKaonKinkPIDFake(0), fhAODmotherVertexRPhiZ(0), fhAODmotherVertexDCARPhiZ(0),
  fhRecPhiInvMassPtYtest(0), fhRecPhiInvMassPtYtestLikeSign(0), fhRecPhiPtYtest(0), fhRecPhiPtYtestLikeSign(0), fhRecPhiInvMassPt(0), fhRecPhiInvMassPtLikeSign(0), fhRecPhiInvMassBothKinks(0), fhRecPhiInvMassBothKinksLikeSign(0), fhRecKstarPtYtest(0), fhRecKstarPtYtestLikeSign(0), fhRecKstarInvMassPt(0), fhRecKstarInvMassPtLikeSign(0), fhRecKstarInvMassBothKinks(0), fhRecKstarInvMassBothKinksLikeSign(0), fhRecLambdaPtYtest(0), fhRecLambdaPtYtestLikeSign(0), fhRecLambdaInvMassPt(0), fhRecLambdaInvMassPtLikeSign(0)
	{}
  AliAnalysisTaskKinksFilimon(const char *name);
  AliAnalysisTaskKinksFilimon(const char *name, const AliVEvent::EOfflineTriggerTypes offlineTriggerType, AliESDkinkCuts* esdKinkCuts, AliESDtrackCuts* esdTrackCutsKinkMother, AliESDtrackCuts* esdTrackCutsKinkDaughter, AliESDtrackCuts* esdTrackCutsPartner, TArrayF esdPIDcutsResonances, /*const AliCFParticleGenCuts* genTrackCuts, const TH1* histPtTemplate=0,*/ const Bool_t useMC=kTRUE, const Bool_t fillCutHist=kFALSE, const TH2* histPtYTemplate=0, const THnSparseF* histRPhiZTemplate=0, const ECollisionType collisionType=kPP, TClass* outputContClass=TList::Class(), const UInt_t aodFilterMap=AliAODTrack::kTrkTPCOnly);
  virtual ~AliAnalysisTaskKinksFilimon();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   UserExecMix(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
	
	class AliRecEventCuts : public AliAnalysisCuts {
	public:
    enum { kProcessedEvents = 1, kPhysSelEvents/*, kAcceptedEvents*/ /* Centrality, multiplicity, vtx combined cut */, kVtxRange, kVtxCentral, kVtxNoEvent, kQVector, kTPCasPV, kZeroCont, kNVtxCuts=4 }; // Consult AliSpectraBothEventCuts
		AliRecEventCuts() : AliAnalysisCuts(), fFillCutHist(kFALSE), fOfflineTriggerType(AliVEvent::kMB), fMinMult(0), fMaxMult(0), fMinVtxZ(0), fMaxVtxZ(0), fMinCent(0), fMaxCent(0), fhEventCuts(0x0) {};
		AliRecEventCuts(const char* name, const char* title) : AliAnalysisCuts(name, title), fFillCutHist(kFALSE), fOfflineTriggerType(AliVEvent::kMB), fMinMult(0), fMaxMult(0), fMinVtxZ(0), fMaxVtxZ(0), fMinCent(0), fMaxCent(0), fhEventCuts(0x0) {};
    AliRecEventCuts(const char* name, const char* title, const AliVEvent::EOfflineTriggerTypes offlineTriggerTypes, const Bool_t fillCutHist = kTRUE, const Double_t minMult = -1, const Double_t maxMult = 1e10, const Double_t minVtxZ = -10, const Double_t maxVtxZ = 10, const Double_t minCent = -1, const Double_t maxCent = 101) : AliAnalysisCuts(name, title), fFillCutHist(fillCutHist), fOfflineTriggerType(offlineTriggerTypes), fMinMult(minMult), fMaxMult(maxMult), fMinVtxZ(minVtxZ), fMaxVtxZ(maxVtxZ), fMinCent(minCent), fMaxCent(maxCent), fhEventCuts(0x0) {
	    fhEventCuts = new TH1F(Form("%s::fhEventCuts", GetName()), Form("%s; Cut type; Number of events", GetTitle()), kNVtxCuts, kProcessedEvents, kNVtxCuts);
			fhEventCuts->SetOption("text0");
	    TAxis* axis = fhEventCuts->GetXaxis();
	    axis->SetBinLabel(kProcessedEvents, "processed");
	    axis->SetBinLabel(kPhysSelEvents, "physics sel");
	    //axis->SetBinLabel(kAcceptedEvents, "accepted");
	    axis->SetBinLabel(kVtxRange, "vtx range");
		};
		virtual Bool_t IsSelected(TObject*) { return kFALSE; };		
		virtual Bool_t IsSelected(TList*) { return kFALSE; };		
		inline Bool_t IsAccepted(const AliVEvent* /*recEvent*/) const {
		if (fFillCutHist) fhEventCuts->Fill(kProcessedEvents);
		//if (TMath::Abs(position.Z()) > fMaxAbsKinkZ) return kFALSE;
		if (fFillCutHist) fhEventCuts->Fill(kPhysSelEvents);
		//if (position.Perp() < fMinKinkR) return kFALSE;
		if (fFillCutHist) fhEventCuts->Fill(kVtxRange);
		//if (position.Perp() > fMaxKinkR) return kFALSE; 
		return kTRUE;
		};
		inline TH1* GetHist() const { return(fhEventCuts); };
	private:
		Bool_t fFillCutHist;
    AliVEvent::EOfflineTriggerTypes fOfflineTriggerType;
		Double_t fMinMult;
		Double_t fMaxMult;
		Double_t fMinVtxZ;
		Double_t fMaxVtxZ;
		Double_t fMinCent;
		Double_t fMaxCent;
		TH1* fhEventCuts;
    AliRecEventCuts(const AliRecEventCuts& rhs);
    AliRecEventCuts& operator=(const AliRecEventCuts& rhs);
  //ClassDef(AliRecEventCuts, 1); // example of analysis
	};

	class AliRecCentCuts : public AliAnalysisCuts {
	public:
    //enum { kProcessedEvents = 1, kPhysSelEvents/*, kAcceptedEvents*/ /* Centrality, multiplicity, vtx combined cut */, kVtxRange, kVtxCentral, kVtxNoEvent, kQVector, kTPCasPV, kZeroCont, kNVtxCuts=4 }; // Consult AliSpectraBothEventCuts
		AliRecCentCuts() : AliAnalysisCuts(), fFillCutHist(kFALSE), fCentEstimator(0x0), fhCentTemplate(0x0), /*fOfflineTriggerType(AliVEvent::kMB), fMinMult(0), fMaxMult(0), fMinVtxZ(0), fMaxVtxZ(0), fMinCent(0), fMaxCent(0),*/ fhEventCuts(0x0) {};
		AliRecCentCuts(const char* name, const char* title) : AliAnalysisCuts(name, title), fFillCutHist(kFALSE), fCentEstimator(0x0), fhCentTemplate(0x0), /*fOfflineTriggerType(AliVEvent::kMB), fMinMult(0), fMaxMult(0), fMinVtxZ(0), fMaxVtxZ(0), fMinCent(0), fMaxCent(0),*/ fhEventCuts(0x0) {};
    AliRecCentCuts(const char* name, const char* title, const char* centEstimator, TH1* const histCentTemplate, const Bool_t fillCutHist = kTRUE/*const AliVEvent::EOfflineTriggerTypes offlineTriggerTypes, const Double_t minMult = -1, const Double_t maxMult = 1e10, const Double_t minVtxZ = -10, const Double_t maxVtxZ = 10, const Double_t minCent = -1, const Double_t maxCent = 101*/) : AliAnalysisCuts(name, title), fFillCutHist(fillCutHist), fCentEstimator(centEstimator), fhCentTemplate(histCentTemplate), /*fOfflineTriggerType(offlineTriggerTypes), fMinMult(minMult), fMaxMult(maxMult), fMinVtxZ(minVtxZ), fMaxVtxZ(maxVtxZ), fMinCent(minCent), fMaxCent(maxCent),*/ fhEventCuts(0x0) {
#if 0
	    fhEventCuts = new TH1F(Form("%s::fhEventCuts", GetName()), Form("%s; Cut type; Number of events", GetTitle()), kNVtxCuts, kProcessedEvents, kNVtxCuts);
			fhEventCuts->SetOption("text0");
	    TAxis* axis = fhEventCuts->GetXaxis();
	    axis->SetBinLabel(kProcessedEvents, "processed");
	    axis->SetBinLabel(kPhysSelEvents, "physics sel");
	    //axis->SetBinLabel(kAcceptedEvents, "accepted");
	    axis->SetBinLabel(kVtxRange, "vtx range");
#endif
		};
		virtual Bool_t IsSelected(TObject*) { return kFALSE; };		
		virtual Bool_t IsSelected(TList*) { return kFALSE; };
#if 0		
		inline Bool_t IsAccepted(const AliVEvent* /*recEvent*/) const {
		if (fFillCutHist) fhEventCuts->Fill(kProcessedEvents);
		//if (TMath::Abs(position.Z()) > fMaxAbsKinkZ) return kFALSE;
		if (fFillCutHist) fhEventCuts->Fill(kPhysSelEvents);
		//if (position.Perp() < fMinKinkR) return kFALSE;
		if (fFillCutHist) fhEventCuts->Fill(kVtxRange);
		//if (position.Perp() > fMaxKinkR) return kFALSE; 
		return kTRUE;
		};
#endif
		inline TH1* GetHist() const { return(fhEventCuts); };
	private:
		Bool_t fFillCutHist;
		const char* fCentEstimator;
		TH1* /*const*/ fhCentTemplate;
    /*AliVEvent::EOfflineTriggerTypes fOfflineTriggerType;
		Double_t fMinMult;
		Double_t fMaxMult;
		Double_t fMinVtxZ;
		Double_t fMaxVtxZ;
		Double_t fMinCent;
		Double_t fMaxCent;*/
		TH1* fhEventCuts;
    AliRecCentCuts(const AliRecCentCuts& rhs);
    AliRecCentCuts& operator=(const AliRecCentCuts& rhs);
  //ClassDef(AliRecCentCuts, 1); // example of analysis
	};

	inline const AliVVertex* GetPrimaryVertex(const AliESDEvent* const esdEvent) const;
	inline const AliVVertex* GetPrimaryVertex(const AliAODEvent* const aodEvent) const;
	inline Bool_t GetPrimaryVertexQuality(const AliVVertex* const vertex, const Int_t nContrib, const Double_t chi2PerNDF) const;
	inline Bool_t GetKinkPID(AliESDkink* /*kink*/, AliESDEvent* /*esdEvent*/, /* Future general function AliESDkinkCuts* esdKinkCuts,*/ AliESDtrack*& /*kinkMother*/, AliESDtrack*& /*kinkDaughter*/, Int_t& /*kinkPID*/) const { return kFALSE; };
	inline Bool_t IsGoodKink(const AliMCParticle* const mcParticle) const;
	Int_t EventSelection(AliVEvent* const inputEvent, AliESDEvent* const esdEvent, const AliAODEvent* const aodEvent, const Float_t centralityF/*Option_t **/);
	Float_t GetRecCentrality(AliVEvent* const inputEvent, const char* centEstimator="V0M"/*AliAnalysisEvent::kCentEst_V0M*/) const;
	//inline AliVVertex* GetRecPrimaryVertex(const AliESDEvent* const esdEvent, const AliAODEvent* const aodEvent, const Int_t nContrib) const {   return(esdEvent ? GetESDVertex(esdEvent, nContrib) : (aodEvent ? aodEvent->GetPrimaryVertex() : 0 )); };
	Int_t KinkAnalysis(/*Option_t **/);
	Int_t KinkResonanceAnalysis(const AliVEvent* const recEvent, const AliVParticle* const kinkMother, const TLorentzVector& mother4Momentum, const AliPID::EParticleType kinkMotherPID, const AliVEvent* const mixEvent /*= 0x0*/, /*or Bool_t eventMixing = kFALSE if more than 1 event needed?*/AliMCEvent* const mcEvent /*= 0x0*/) const;
  inline Int_t GetDaughterIdx(const AliMCParticle* const mcTrack, const Int_t mcNtracks, Int_t& firstDaughterIdx, Int_t& lastDaughterIdx) const;
	inline Int_t GetLastTrackRef(AliMCParticle* const mcTrack, AliTrackReference*& lastTrackReference, TVector3& decay3Momentum, TVector3& decay3Position);
	void PostProcessHist() {};
	Bool_t IsPrimary(const AliESDtrack* /*esdTrack*/) const { return kFALSE; };
	Bool_t IsPrimary(const AliVParticle* /*particle*/) const { return kFALSE; };
  virtual Bool_t IsSelected(const AliVParticle *track, const AliPID::EParticleType type) const { return( (fPIDResponse && (fEsdPIDCutsArray.GetSize() == AliPID::kSPECIES) ) ? (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, type)) < fEsdPIDCutsArray[type]) : kTRUE); };
	Int_t BookCentralityHist(const char* name, const char* title, const TH2* const histPtYtemplate, const Int_t nCentBins, const Float_t* centBins, TObjArray* const histArray, TCollection* const outputCont);
	inline AliVEvent* GetMainEvent() const;
	AliVEvent* GetMixedEvent() const;
	Int_t CheckMCtruth(AliVTrack* const kinkMotherRec, const TLorentzVector& mother4Momentum, const AliPID::EParticleType kinkMotherPID, AliVParticle* const kinkMotherMC, const Bool_t isPhysicalPrimary) const;
	
  TClass* fOutputContClass;
	TCollection* fOutputCont; //! Output list
	AliVEventHandler* fMainInputHandler;
  Bool_t fIsAOD;
  Bool_t fUseMC;
	const Bool_t fFillCutHist;
	const Bool_t fEventSelectionCutsMC;
	const ECollisionType fCollisionType;
  const AliVEvent::EOfflineTriggerTypes fOfflineTriggerType;
	const AliPID::EParticleType fPartIdDedx;
	AliPIDResponse* fPIDResponse;
	AliRecEventCuts* fRecEventCuts;
	AliRecCentCuts* fRecCentCuts;
	AliKineTrackCuts* fCommonKineTrackCuts; //AliKineTrackCuts*/
	AliESDtrackCuts* fEsdTrackCutsKinkMother;
	AliESDtrackCuts* fEsdTrackCutsKinkDaughter;
	AliESDtrackCuts* fEsdTrackCutsPartner;
	const TArrayF fEsdPIDCutsArray;
	AliESDkinkCuts* fRecKinkCutsKaon;
	UInt_t fAODFilterMap;
	//const AliESDkinkCuts* fRecKinkCutsPion;
	//AliESDkinkCuts* fRecKinkCutsKaon;
	//const AliESDkinkCuts* fMcKinkCutsPion;
	const TH2* fHistPtYTemplate;
	const THnSparseF* fHistRPhiZTemplate;
	Float_t fKaonMass;
	Float_t fMuonMass;
	Float_t fPionMass;
	Float_t fProtonMass;
	Float_t fElectronMass;
	Float_t fPhiMassRangeMin;
	Float_t fPhiMassRangeMax;
	Float_t fKstarMassRangeMin;
	Float_t fKstarMassRangeMax;
	Float_t fLambdaMassRangeMin;
	Float_t fLambdaMassRangeMax;
	
	// Rec hist
	TH1* fhRecMultUnbiased;
	TH1* fhRecMult;
	TH2* fhRecMultCentralityUnbiased;
	TH2* fhRecMultCentrality;
	TH2* fhRecAllPtY;
	// ESD kinks
	TH1* fhESDnKinks;
	TH1* fhESDKinkQt;
	TH1* fhESDKinkAngle;
	TH2* fhESDKinkPAngle;
	TH1* fhESDKinkDCA;
	TH2* fhESDAllKinkPtY;
	TH2* fhESDUnknownKinkPtY;
	TH2* fhESDKaonKinkPtY;
	TH2* fhESDPionKinkPtY;
	TH2* fhESDKaonKinkPtEta;
	//TH2* fhESDPionKinkPtEta;
	TH2* fhESDKPlusKinkPtY;
	TH2* fhESDKMinusKinkPtY;
	TH1* fhESDKPlusKinkPt;
	TH1* fhESDKMinusKinkPt;
	TH2* fhESDPiPlusKinkPtY;
	TH2* fhESDPiMinusKinkPtY;
	TH3* fhESDKaonKinkPtYCent;
	TH3* fhESDPionKinkPtYCent;
	TH3* fhESDKPlusKinkPtYCent;
	TH3* fhESDPiPlusKinkPtYCent;
	TH3* fhESDKMinusKinkPtYCent;
	TH3* fhESDPiMinusKinkPtYCent;
	THnSparseF* fhRecPrimaryVertexRPhiZ;
	THnSparseF* fhESDAllKinkDecaysRPhiZ;
	THnSparseF* fhESDKaonKinkDecaysRPhiZ;
	THnSparseF* fhESDPionKinkDecaysRPhiZ;
	TH1* fhESDKaonKinkQt;
	TH1* fhESDKaonKinkAngle;
	TH1* fhESDPionKinkQt;
	TH1* fhESDPionKinkAngle;
	TH2* fhESDMomentumTPCSignalKaonKinks;
	TH1* fhESDKaonKinksMotherAndDaughterTPCncls;
	TH2* fhESDKaonKinksMotherVSDaughterTPCncls;
	TH2* fhESDMomentumTPCSignalPionKinks;
	TH1* fhESDPionKinksMotherAndDaughterTPCncls;
	TH2* fhESDPionKinksMotherVSDaughterTPCncls;
	TH1* fhESDKaonOverPionEbE;
	TH1* fhESDInvMassKinkDecayMuNu;
	// MC hist
	//THnSparseF* fhMCKaonKinkDecaysXYZ;
	THnSparseF* fhMCmainVertexRPhiZ;
	THnSparseF* fhMCKaonKinkDecaysRPhiZ;
  TH1* fhMCmult;
	TH1* fhMCmultPrim;
	TH2* fhMCPtYall;
	TH2* fhMCPtYprim;
	TH1* fhMCpdg;
	TH2* fhMCPtYprimKaon;
	TH2* fhMCPtYprimKPlus;
	TH2* fhMCPtYprimKMinus;
	TH1* fhMCPtprimKPlus;
	TH1* fhMCPtprimKMinus;
	TH3* fhMCPtYCentprimKaon;
	TH3* fhMCPtYCentprimKPlus;
	TH3* fhMCPtYCentprimKMinus;
	TH2* fhMCPtYkinkKaonFiducial;
	TH2* fhMCPtYkinkKPlusFiducial;
	TH2* fhMCPtYkinkKMinusFiducial;
	TH1* fhMCKinkKaonQt;
	TH1* fhMCKinkKaonAngle;
	TH1* fhMCKaonKinkDecaysTrackLength;
	TH1* fhMCKaonKinkDecaysTrackTime;
	TH1* fhMCkaonMotherPdg;
	TH1* fhMCpionMotherPdg;
	TH2* fhMCPhiPtY;
	TH2* fhMCPtYprimPhi;
	TH3* fhMCPhiInvMassPtYtest;
	TH2* fhMCInvMassPtprimPhi;
	TH3* fhMCInvMassPtCentprimPhi;
	TH2* fhMCKstarPtY;
	TH2* fhMCPtYprimKstar;
	TH2* fhMCInvMassPtprimKstar;
	TH3* fhMCInvMassPtCentprimKstar;
	TH2* fhMCLambdaPtY;
	TH2* fhMCPtYprimLambda;
	TH2* fhMCInvMassPtprimLambda;
	TH3* fhMCInvMassPtCentprimLambda;
	TH2* fhMCPhi2KaonPtY;
	TH2* fhMCKstar2KaonPtY;
	TH2* fhMCKaon2MuonPtY;
	TH2* fhMCKaon2PionPtY;
	TH2* fhMCPtYprimPion;
	TH2* fhMCPtYprimPiPlus;
	TH2* fhMCPtYprimPiMinus;
	TH3* fhMCPtYCentprimPion;
	TH3* fhMCPtYCentprimPiPlus;
	TH3* fhMCPtYCentprimPiMinus;
	THnSparseF* fhMCPionKinkDecaysRPhiZ;
	TH2* fhMCPion2MuonPtY;
	TH2* fhMCPtYkinkPionFiducial;
	TH2* fhMCPtYkinkPiPlusFiducial;
	TH2* fhMCPtYkinkPiMinusFiducial;
	TH1* fhMCKinkPionQt;
	TH1* fhMCKinkPionAngle;
	TH2* fhMCMotherDaughterPdg;
	TH1* fhMCKaonOverPionEbE;
	TH2* fhMCPhiDecayOpeningAngle;
	TH2* fhMCKstarDecayOpeningAngle;
	TH2* fhMCLambdaDecayOpeningAngle;
	// MC truth
	TH1* fhESDKinkQtTrueMC;
	TH1* fhESDKinkAngleTrueMC;
	TH2* fhESDKaonKinkPtYTrueMC;
	TH2* fhESDKaonKinkPtYFakeMC;
	TH2* fhESDKaonKinkPtYTrueMCsecondary;
	TH2* fhESDKPlusKinkPtYTrueMC;
	TH2* fhESDKPlusKinkPtYFakeMC;
	TH2* fhESDKPlusKinkPtYTrueMCsecondary;
	TH2* fhESDKMinusKinkPtYTrueMC;
	TH2* fhESDKMinusKinkPtYFakeMC;
	TH2* fhESDKMinusKinkPtYTrueMCsecondary;
	TH2* fhESDPionKinkPtYTrueMC;
	TH2* fhESDPionKinkPtYFakeMC;
	TH2* fhESDPionKinkPtYTrueMCsecondary;
	TH2* fhESDPiPlusKinkPtYTrueMC;
	TH2* fhESDPiPlusKinkPtYFakeMC;
	TH2* fhESDPiPlusKinkPtYTrueMCsecondary;
	TH2* fhESDPiMinusKinkPtYTrueMC;
	TH2* fhESDPiMinusKinkPtYFakeMC;
	TH2* fhESDPiMinusKinkPtYTrueMCsecondary;
	TH2* fhESDPhi2KaonPtYTrueMC;
	TH2* fhESDPhiPtYtestTrueMC;
	TH2* fhESDPhiPtYtestFakeMC;
	TH2* fhESDPhiPtYtestLikeSignTrueMC;
	TH2* fhESDPhiPtYtestLikeSignFakeMC;
	TH2* fhESDKstar2KaonPtYTrueMC;
	TH2* fhESDKstarPtYtestTrueMC;
	TH2* fhESDKstarPtYtestFakeMC;
	TH2* fhESDKstarPtYtestLikeSignTrueMC;
	TH2* fhESDKstarPtYtestLikeSignFakeMC;
	TH2* fhESDLambda2KaonPtYTrueMC;
	// AOD test
	TH1* fhAODnVertices;
	TH1* fhAODnKinkVertices;
	TH1* fhAODnKinkDaughters;
	THnSparseF* fhAODAllKinkDecaysRPhiZ;
	TH1* fhAODKinkKaonQt;
	TH1* fhAODKinkKaonAngle;
	TH2* fhAODKinkMothervsDaughterPIDcheck;
	TH1* fhAODKaonKinkPID;
	TH1* fhAODKaonKinkPIDFake;
	THnSparseF* fhAODmotherVertexRPhiZ;
	THnSparseF* fhAODmotherVertexDCARPhiZ;
	// Rec resonances
	TH3* fhRecPhiInvMassPtYtest;
	TH3* fhRecPhiInvMassPtYtestLikeSign;
	TH2* fhRecPhiPtYtest;
	TH2* fhRecPhiPtYtestLikeSign;
	TH2* fhRecPhiInvMassPt;
	TH2* fhRecPhiInvMassPtLikeSign;
	TH1* fhRecPhiInvMassBothKinks;
	TH1* fhRecPhiInvMassBothKinksLikeSign;
	TH2* fhRecKstarPtYtest;
	TH2* fhRecKstarPtYtestLikeSign;
	TH2* fhRecKstarInvMassPt;
	TH2* fhRecKstarInvMassPtLikeSign;
	TH1* fhRecKstarInvMassBothKinks;
	TH1* fhRecKstarInvMassBothKinksLikeSign;
	TH2* fhRecLambdaPtYtest;
	TH2* fhRecLambdaPtYtestLikeSign;
	TH2* fhRecLambdaInvMassPt;
	TH2* fhRecLambdaInvMassPtLikeSign;
  
  AliAnalysisTaskKinksFilimon(const AliAnalysisTaskKinksFilimon&); // not implemented
  AliAnalysisTaskKinksFilimon& operator=(const AliAnalysisTaskKinksFilimon&); // not implemented

  ClassDef(AliAnalysisTaskKinksFilimon, 1);
};

//________________________________________________________________________
AliVEvent* AliAnalysisTaskKinksFilimon::GetMainEvent() const {

	return fMainInputHandler->GetEvent();//fMultiInputHandler ? fMultiInputHandler->GetFirstInputEventHandler()->GetEvent() : InputEvent();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskKinksFilimon::IsGoodKink(const AliMCParticle* const mcParticle) const {
	
	if ( !mcParticle ) return(kFALSE);
	TParticle* particle = mcParticle->Particle();
	if ( !particle ) return(kFALSE);
	return kTRUE;
}

//________________________________________________________________________
const AliVVertex* AliAnalysisTaskKinksFilimon::GetPrimaryVertex(const AliESDEvent* const esdEvent) const
{
	
	if ( !esdEvent ) return(0x0);
  const AliESDVertex* esdVertex = esdEvent->GetPrimaryVertex();
  if ( esdVertex && esdVertex->GetStatus() ) return esdVertex;
  esdVertex = esdEvent->GetPrimaryVertexTracks();
  if ( esdVertex && esdVertex->GetStatus() ) return esdVertex;
  esdVertex = esdEvent->GetPrimaryVertexSPD();
  if ( esdVertex && esdVertex->GetStatus() ) return esdVertex;
	/*esdVertex = esdEvent->GetPrimaryVertexTPC();
  if ( esdVertex && esdVertex->GetStatus() ) return esdVertex;*/
  return(0x0);
}

//________________________________________________________________________
const AliVVertex* AliAnalysisTaskKinksFilimon::GetPrimaryVertex(const AliAODEvent* const aodEvent) const
{
	
	if ( !aodEvent ) return(0x0);
  const AliAODVertex* aodVertex = aodEvent->GetPrimaryVertex();
  if ( aodVertex /*&& aodVertex->GetStatus()*/ ) return aodVertex;
  /*aodVertex = aodEvent->GetPrimaryVertexTracks();
  if ( aodVertex && aodVertex->GetStatus() ) return aodVertex;*/
  aodVertex = aodEvent->GetPrimaryVertexSPD();
  if ( aodVertex /*&& aodVertex->GetStatus()*/ ) return aodVertex;
	/*aodVertex = aodEvent->GetPrimaryVertexTPC();
  if ( aodVertex && aodVertex->GetStatus() ) return aodVertex;*/
  return(0x0);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskKinksFilimon::GetPrimaryVertexQuality(const AliVVertex* const vertex, const Int_t nContrib, const Double_t chi2PerNDF) const
{
	if (!vertex) return kFALSE;
	if (vertex->GetNContributors() < nContrib) return kFALSE;
	if (vertex->GetChi2perNDF() > chi2PerNDF) return kFALSE;
	return kTRUE;
}

//________________________________________________________________________
Int_t AliAnalysisTaskKinksFilimon::GetDaughterIdx(const AliMCParticle* const mcTrack, const Int_t mcNtracks, Int_t& firstDaughterIdx, Int_t& lastDaughterIdx) const {
	
  if ( !(mcNtracks && mcTrack) ) return(-1);
	firstDaughterIdx = mcTrack->GetDaughterFirst();
	lastDaughterIdx = mcTrack->GetDaughterLast();
	if ( (firstDaughterIdx < 0) || (lastDaughterIdx < 0) || (firstDaughterIdx > mcNtracks) || (lastDaughterIdx > mcNtracks) ) return(-1); // Invalid decay
  return(lastDaughterIdx-firstDaughterIdx+1);
	
}

//________________________________________________________________________
Int_t AliAnalysisTaskKinksFilimon::GetLastTrackRef(AliMCParticle* const mcTrack, AliTrackReference*& lastTrackReference, TVector3& decay3Momentum, TVector3& decay3Position) {
	
	if ( !mcTrack ) return(0);
	Int_t nTrackRefs = mcTrack->GetNumberOfTrackReferences();
	if ( nTrackRefs < 1 ) return(nTrackRefs); // At least 1 trackref for decay
	lastTrackReference = mcTrack->GetTrackReference(nTrackRefs-1);
	if ( !lastTrackReference ) return(0);
	decay3Momentum.SetXYZ(lastTrackReference->Px(), lastTrackReference->Py(), lastTrackReference->Pz());
	decay3Position.SetXYZ(lastTrackReference->X(), lastTrackReference->Y(), lastTrackReference->Z());
	return(nTrackRefs);
}

//________________________________________________________________________
inline Double_t Rapidity(const AliVParticle* const part) {
  Double_t e  = part->E();
  Double_t pz = part->Pz(); 
  if (e > TMath::Abs(pz)) return 0.5*TMath::Log((e+pz)/(e-pz));
  else return -999.;
}

//________________________________________________________________________
inline Double_t Rapidity(const TLorentzVector& part) {
  Double_t e  = part.E();
  Double_t pz = part.Pz(); 
  if (e > TMath::Abs(pz)) return 0.5*TMath::Log((e+pz)/(e-pz));
  else return -999.;
}

//________________________________________________________________________
inline Int_t GetCustomBins(const Int_t nBins, const Double_t lowBin, const Double_t upBin, Double_t* const binArray) {
	/*Double_t binWidth = (upBin-lowBin)/nBins, currentBin = lowBin;
  for ( Int_t iBin = 0; iBin < nBins; ++iBin) {
		binArray[iBin] = currentBin;
		currentBin+=binWidth;
	}*/
	TH1* templateHist = new TH1F("templateHist", "templateHist", nBins, lowBin, upBin);
	if (!templateHist) return(-1);
  for ( Int_t iBin = 0; iBin < nBins; ++iBin) {
		binArray[iBin] = templateHist->GetXaxis()->GetBinLowEdge(iBin+1);
	}
	delete templateHist;
	return(0);
}

//________________________________________________________________________
inline Bool_t IsInFiducialKine(const TLorentzVector& momentum, const Double_t minPt, const Double_t maxAbsY) {
	return( (momentum.Pt() > minPt) && (TMath::Abs(momentum.Rapidity()) < maxAbsY) );
}

//________________________________________________________________________
inline Bool_t PassesESDTrackCuts(AliAODTrack* const aodTrack, AliESDtrackCuts* const esdTrackCuts) {
	Float_t tempMin = 0, tempMax = 0;
	esdTrackCuts->GetPtRange(tempMin, tempMax);
	if ( (aodTrack->Pt() < tempMin) || (aodTrack->Pt() > tempMax) ) return kFALSE;
	esdTrackCuts->GetRapRange(tempMin, tempMax);
	if ( (aodTrack->Y() < tempMin) || (aodTrack->Y() > tempMax) ) return kFALSE;
	if ( aodTrack->GetTPCNcls() < esdTrackCuts->GetMinNClusterTPC() ) return kFALSE;
	if ( aodTrack->Chi2perNDF() > esdTrackCuts->GetMaxChi2PerClusterTPC() ) return kFALSE;
	//if (esdTrackCuts->GetMaxChi2TPCConstrainedGlobal() < aodTrack->Chi2perNDF() ) return kFALSE;
	if ( esdTrackCuts->GetRequireTPCRefit() && !aodTrack->IsOn(AliVTrack::kTPCrefit) ) return kFALSE;
	if ( esdTrackCuts->GetRequireITSRefit() && !aodTrack->IsOn(AliVTrack::kITSrefit) ) return kFALSE;
	//if ( !esdTrackCuts->GetAcceptKinkDaughters() && aodTrack->GetProdVertex() && (aodTrack->GetProdVertex()->GetType() == AliAODVertex::kKink) ) return kFALSE;
	return kTRUE;
}

#endif
