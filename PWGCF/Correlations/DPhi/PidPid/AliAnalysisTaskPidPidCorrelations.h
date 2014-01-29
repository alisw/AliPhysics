#ifndef ALIANALYSISTASKPIDPIDCORRELATION_H
#define ALIANALYSISTASKPIDPIDCORRELATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TAxis.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TObjArray.h"

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliAODHeader.h"
#include "AliVParticle.h"
#include "AliVVertex.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliInputEventHandler.h" // event mixing
#include "AliEventPoolManager.h"  // event mixing
#include "AliLog.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"

#include "AliCFContainer.h"
// #include "AliCFGridSparse.h"
// #include "AliCFEffGrid.h"
// #include "AliCFManager.h"
// #include "AliCFCutBase.h"

#include "TH1F.h"
#include "TH2F.h"

#include "THnSparse.h"
// #include "AliTHn.h"

class TList;
class TH1F;
class TH2F;
// class AliTHn;
class AliCFContainer;
// class AliCFGridSparse;
// class AliCFEffGrid;
// class AliCFCutBase;
// class AliCFManager;
class AliPIDResponse;

#include <vector>
using std::vector;

// const Int_t kTrackVariablesSingle = 4+1;       // track variables in histogram (particle type, pTtrig, vertexZ, centrality)
const Int_t kTrackVariablesPair   = 5;       // track variables in histogram (particle type trig, particle type assoc, pTtrig, ptAssociated, dPhi, dEta, vertexZ, centrality)

namespace myAliPID {
  
  //enum PIDType { kMyNSigmaITS=0, kMyNSigmaTPC, kMyNSigmaTOF, kMyNSigmaHMP, kMyNSigmaTPCTOF, kMyNSigmaPIDType = kMyNSigmaTPCTOF };
  enum PIDType { kMyNSigmaITS=0, kMyNSigmaTPC, kMyNSigmaTOF, kMyNSigmaHMP, kMyNSigmaPIDType };
  enum AliParticleSpeciesQA { fPartElectronQA=0, fPartMuonQA, fPartPionQA, fPartKaonQA, fPartProtonQA, fPartNSpeciesQA, fPartUndefinedQA = 999 };
  enum AliParticleSpecies { fPartPionMinus=1, fPartPionPlus, fPartKaonMinus, fPartKaonPlus, fPartProtonMinus, fPartProtonPlus, fPartHadronMinus, fPartHadronPlus, fPartNSpecies, fPartUndefined = 999 };
}

using namespace myAliPID;

class AliAnalysisTaskPidPidCorrelations : public AliAnalysisTaskSE { //

  public:
  
  AliAnalysisTaskPidPidCorrelations();
  AliAnalysisTaskPidPidCorrelations(const Char_t* name/* = "AliAnalysisTaskPidPidCorrelations"*/);
  virtual ~AliAnalysisTaskPidPidCorrelations();
  
  //___ Implementation of interface methods
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t* option);
  virtual void  Terminate(Option_t* option);

  // Setters/Getters
  void 	SetMaxVertexZ(Double_t maxVertexZ) { fVzMax = maxVertexZ; }
  void	SetTriggerMask(Int_t triggerType)   { fTriggerType = triggerType; }  
//   void	SetFilterBit(Int_t fbit)		{ fFilterBit = fbit; }
  void	SetCentralityEstimator(const Char_t* centralityestimator) { fCentralityEstimator = centralityestimator; }
  void 	SetCentralityRange(Float_t min, Float_t max)    { fCentralityPercentileMin = min; fCentralityPercentileMax = max; }
//   void	SetEtaRange(Double_t etaMin, Double_t etaMax) { fTrackEtaMin = etaMin; fTrackEtaMax = etaMax; }
  void	SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) { fVxMax = vx; fVyMax = vy; fVzMax = vz;}  
  void	SetKinematicsCutsAOD(Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax) { fTrackPtMin = ptmin;  fTrackPtMax = ptmax; fTrackEtaMin = etamin; fTrackEtaMax = etamax; }
  //void SetFilterMask(UInt_t i,Int_t iType = 0){fFilterMask = i;fFilterType = iType;}

  void	SetMixingTracks(Int_t tracks, Int_t poolsize) { fMinNumTrack = tracks; fPoolSize = poolsize;}
  void	SetFillpT(Bool_t flag) { fFillpT = flag; }
  void	SetTwoTrackEfficiencyCut(Float_t value = 0.02, Float_t min = 0.8) { twoTrackEfficiencyCutValue = value; fTwoTrackCutMinRadius = min; }
  void	SetPairCuts(Bool_t conversions, Bool_t resonances) { fCutConversions = conversions; fCutResonances = resonances; }
  void	SetWeightPerEvent(Bool_t flag = kTRUE)   { fWeightPerEvent = flag; }
  void	SetPtOrder(Bool_t flag) { fPtOrder = flag; }
  void	SetEtaOrdering(Bool_t flag) { fEtaOrdering = flag; }
  void	SetEventMixing(Bool_t flag) { fFillMixed = flag; }
  void	SetRejectResonanceDaughters(Int_t value) { fRejectResonanceDaughters = value; }
  void	SetSelectCharge(Int_t selectCharge) { fSelectCharge = selectCharge; }
  void	SetSelectTriggerCharge(Int_t selectCharge) { fTriggerSelectCharge = selectCharge; }
  void	SetSelectAssociatedCharge(Int_t selectCharge) { fAssociatedSelectCharge = selectCharge; }
  void	SetTriggerRestrictEta(Float_t eta) { fTriggerRestrictEta = eta; }
  void	SetOnlyOneEtaSide(Int_t flag)     { fOnlyOneEtaSide = flag; }

  void SetPIDsToCorrelate(Int_t triggerPID, Int_t assocPID) { fPIDtrigType = triggerPID; fPIDassocType = assocPID; }
  
  void 		SetMC(Bool_t isMC) {fUseMC = isMC;}
  Bool_t		GetMC() const {return fUseMC;}

  void 	SetCentBinning(Int_t nBins, Double_t* limits);
  Int_t	GetCentBin(Double_t centrality);
  void 	SetZvtxBinning(Int_t nBins, Double_t* limits);
  Int_t	GetZvtxBin(Double_t zvtx);
  void 	SetPtBinning(Int_t nBins, Double_t* limits);
  Int_t	GetPtBin(Double_t valPt);
  void	SetEtaBinning(Int_t nBins, Double_t *limits);
  Int_t	GetEtaBin(Double_t valEta);  
    
  TString	GetCentralityEstimator() const { return fCentralityEstimator; }
  Double_t	GetMaxVertexZ() const		{ return fVzMax; }
//   Int_t		GetMyFilterBit()		{ return fFilterBit; }
//  Double_t* GetNSigmas(AliParticleSpecies species) { return fnsigmas[species]; }	// get nsigma[ipart][idet], calculated in CalculateNSigmas(trk)
/*
  //_______ AliTHn
  AliTHn*	GetHistCorrSingle() {return fHistCorrSingle;}
  void 		SetHistCorrSingle(AliTHn *gHist) { fHistCorrSingle = gHist; }//fHistP->FillParent(); fHistP->DeleteContainers();}
  AliTHn*	GetHistCorrPairSame() {return fHistCorrPairSame;}
  void 		SetHistCorrPairSame(AliTHn *gHist) { fHistCorrPairSame = gHist; }//fHistN->FillParent(); fHistN->DeleteContainers();}
  AliTHn*	GetHistCorrPairMixed() {return fHistCorrPairMixed;}
  void 		SetHistCorrPairMixed(AliTHn *gHist) { fHistCorrPairMixed = gHist; }//fHistN->FillParent(); fHistN->DeleteContainers();}
*/
  //AliTHn* MakeHistCorr(const Char_t* name) const;

  void		UseMomentumDifferenceCut(Bool_t fqcut = kFALSE,Double_t gDeltaPtCutMin = 0.01) { fQCut = fqcut; fDeltaPtMin = fqcut ? gDeltaPtCutMin : 0.; }
  void 		SetupForMixing();
  void  		AddSettingsTree();
  void  		Analyse();
	
  void		FillCorrelations(TObjArray* particles, TObjArray* particlesMixed, Double_t centrality, Double_t zVtx, Double_t bSign, Bool_t twoTrackEfficiencyCut, Double_t twoTrackEfficiencyCutValue, /*Int_t step,*/ Double_t weight);
  Bool_t 	CheckMcDistributions(TClonesArray* arrayMC, AliAODMCHeader* mcHeader);
  TString 	GetGenerator(Int_t label, AliAODMCHeader* MCheader);
  Bool_t 	IsFromHijingBg(Int_t mcTrackLabel);
  void 		FillMcGeneratorHistos(TString genLabel);
  Bool_t  	VertexSelection(TObject* obj, Int_t ntracks, Int_t centBin, Double_t gVxMax, Double_t gVyMax, Double_t gVzMax);
  void 		CleanUp(TObjArray* tracks, TObject* mcObj);
  TObjArray* 	AcceptTracks(Int_t centBin, TObject* arrayMC, /*Bool_t onlyprimaries,*/ Bool_t useCuts);
  TObjArray*    	AcceptMcTracks(Int_t centBin, Bool_t onlyprimaries, Bool_t useCuts);
  TObjArray* 	AcceptMcRecoMachedTracks(Int_t centBin, Bool_t onlyprimaries, Bool_t useCuts);
  Double_t* 	GetBinning(const Char_t* configuration, const Char_t* tag, Int_t& nBins);
  Bool_t 	CheckTrack(AliAODTrack* track);
  void 		CalculateNSigmas(AliAODTrack* track, Int_t centBin, Bool_t* pidFlag, Bool_t fillQA);
  Int_t 		FindNSigma(AliAODTrack* track);
  Int_t 		GetParticleID(AliVTrack* trk, Int_t centbin, Bool_t fillQA); // DATA and MC-reco
  Int_t 		GetParticleIDMC(AliVTrack* trk, Int_t centbin, Bool_t fillQA); // MC-truth
  Double_t 	MakeTPCPID(AliAODTrack* track, Double_t* nSigma);
  Double_t 	MakeTOFPID(AliAODTrack* track, Double_t* nSigma);
  Bool_t 	HasTPCPID(AliAODTrack* track) const;
  Bool_t 	HasTOFPID(AliAODTrack* track) const;
  Double_t 	GetBeta(AliAODTrack* track);
  void 		RemoveDuplicates(TObjArray* tracks);
  void 		RemoveWeakDecays(TObjArray* tracks, TObject* mcObj);
  Double_t	DeltaPhi(Double_t Dphi) const;
  TH2F* 		GetHisto2D(const Char_t* name);
  Double_t 	GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius, Double_t bSign);
  Float_t 	GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
  Float_t	GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
  void 		PrintPoolManagerContents();
  TObjArray* 	CloneAndReduceTrackList(TObjArray* tracks);
  TString 	GetOutputListName() const;

  //___________ Correction
  enum ECorrStep { kStepGen = 0, kStepRec, kStepRecMatch, kNsteps };
//   enum ECorrVars { kVarEta, kVarPhi, kVarPt, kVarPID, kVarCent, kVarZvtx, kNvars };
//   enum ECorrVars { kVarEta, kVarPhi, kVarPt, kVarZvtx, kNvars };
  enum ECorrVars { kVarEta, kVarPt, kVarPID, kVarCent, kVarZvtx, kNvars };
  void FillCFcontainers(TObjArray* mca, TObjArray* rec, TObjArray* recmatch, Double_t cent/*, Double_t zvtx*/);

  private:

  AliAnalysisTaskPidPidCorrelations(const AliAnalysisTaskPidPidCorrelations&);            // not implemented
  AliAnalysisTaskPidPidCorrelations& operator=(const AliAnalysisTaskPidPidCorrelations&); // not implemented

  Bool_t fUseMC;

  AliAODEvent* 		fMyAODEvent; //!
  AliAODHeader*		fMyAODHeader; //!
  AliAODTrack*		fMyAODTrack; //!
  AliPIDResponse*	fPIDResponse; //! PID response object
  AliAODVertex*		fMyPrimVtx; //!
  TClonesArray*		fMyMcArray; //!
  AliAODMCHeader*   	fMyMcHeader; //!
  AliMCEventHandler*	fMyMcHandler; //!
  AliEventPoolManager*	fPoolMgr; //!
  AliCFContainer* 	fMyCFCont; //!
  
  TObjArray* 		fMyprimRecoTracksPID; //!
  TObjArray* 		fMyprimMCParticlesPID; //!
  TObjArray* 		fMyprimRecoTracksMatchedPID;  //!
  
  
//   Bool_t		fWriteCorrTree;
//   TTree* 	fVariablesTreeCorr;
//   Float_t*	fCorrVariables;
  
  Int_t 		fTriggerType; 	//  sets trigger -> AliVEvent::kMB, AliVEvent::kHighMult
  Int_t		fMyMcType ;//
//   Int_t 		fFilterBit; 		// track selection cuts
  UInt_t  	fFilterType;    // filter type 0 = all, 1 = ITSTPC, 2 = TPC
  Double_t 	fCentrality;//
  Double_t 	fCentralityPercentileMin;
  Double_t 	fCentralityPercentileMax;
  Int_t 		fNbinsCent;
  TAxis*		fCentAxis; //
  Int_t		fNbinsZvtx;
  TAxis*		fZvtxAxis; //
  Int_t		fNbinsPt;
  TAxis*		fPtAxis; //
  Int_t 		fNbinsEta;
  TAxis*		fEtaAxis; //
  TString 	fCentralityEstimator; // Method to determine the centrality, e.g. "V0M","TRK","TKL","ZDC","FMD", "CL0", "CL1"
  Double_t 	fTrackEtaMin;
  Double_t 	fTrackEtaMax;
  Double_t 	fTrackPtMin;
  Double_t	fTrackPtMax;
  UInt_t		fTrackStatus;          // if non-0, the bits set in this variable are required for each track
  Int_t		fnTracksVertex;        // QA tracks pointing to principal vertex
  Bool_t		fRejectZeroTrackEvents;  // reject events which have no tracks (using the eta, pT cuts defined)
  Double_t 	fEtaCut;
  Double_t 	fVxMax;
  Double_t 	fVyMax;
  Double_t 	fVzMax;
  Bool_t 	fRemoveWeakDecays;   // remove secondaries from weak decays from tracks and particles
  Bool_t 	fRemoveDuplicates;   // remove particles with the same label (double reconstruction)

  Double_t 	fDeltaEtaMax;				// maximum delta eta for output THnSparse

  Int_t		fSelectCharge;           // (un)like sign selection when building correlations: 0: no selection; 1: unlike sign; 2: like sign
  Int_t 		fTriggerSelectCharge;    // select charge of trigger particle: 1: positive; -1 negative
  Int_t 		fAssociatedSelectCharge; // select charge of associated particle: 1: positive; -1 negative
  Float_t 	fTriggerRestrictEta;   // restrict eta range for trigger particle (default: -1 [off])
  Bool_t 	fCutConversions;        // cut on conversions (inv mass)
  Bool_t 	fCutResonances;         // cut on resonances (inv mass)
  Int_t 		fRejectResonanceDaughters; // reject all daughters of all resonance candidates (1: test method (cut at m_inv=0.9); 2: k0; 3: lambda)
  Int_t 		fOnlyOneEtaSide;       // decides that only trigger particle from one eta side are considered (0 = all; -1 = negative, 1 = positive)
  Bool_t 	fWeightPerEvent;	   // weight with the number of trigger particles per event
  Bool_t 	fPtOrder;		   			// apply pT,a < pt,t condition; default: kTRUE
  Bool_t 	fQCut;						//cut on momentum difference to suppress femtoscopic effect correlations
  Double_t	fDeltaPtMin;				//delta pt cut: minimum value
  
  Int_t		fPIDtrigType;	// PID trigger type to correlate
  Int_t		fPIDassocType;	// PID assoc type to correlate

  TString 	fCustomBinning;//for setting customized binning
  TString 	fBinningString;//final binning string	

  Int_t 		fMinNumTrack; // AliEventPoolManager(), Size of track buffer for event mixing (number of tracks to fill the pool)
  Int_t 		fPoolSize; // AliEventPoolManager(), max number of event to mix
  Int_t 		fMinNEventsToMix; //

  Bool_t 	fFillpT;
  Float_t	fTwoTrackEfficiencyCut;   // enable two-track efficiency cut
  Float_t	twoTrackEfficiencyCutValue; // set it 0.02 as default in dphicorrelations
  Float_t	fTwoTrackCutMinRadius;    // minimum radius for two-track efficiency cut
  Bool_t		fEtaOrdering;           // eta ordering, see AliUEHistograms.h for documentation
  Bool_t		fFillMixed;		// enable event mixing (default: ON)

  static const Int_t fNMaxBinsCentrality = 1;
  static const Int_t fNMaxBinsZvtx = 15;
  static const Int_t fNMaxBinsPt = 9;
  static const Int_t fNMaxBinsEta = 25;
  
  TList*		fList; //
  
  TH1F*		fHistQA[fNMaxBinsCentrality][14]; //!
  TH1I*		fHistNev;//!
  TH1F*		fHistTriggerStats; //!
  TH1F*		fHistTriggerRun;  //!
  TH1F*		fHistEventStat; //!
  TH2F*		fHistRefTracks; //!
  TH2F*		fHistRefTracksCent[fNMaxBinsCentrality][6];//!
  TH2F*		fHistCentStats; //!
  TH1F*		fHistCentralityPercentile;        //! centrality class
  TH1F*		fHistCentralityClass10;           //! centrality class by 10
  TH1F*		fHistCentralityClass5;            //! centrality class by 5
  TH2F*		fHistV0M; //!
  
//   TH1F*		fHistDeltaPhi[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt]; //!
//   TH1F*		fHistDeltaPhiMix[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt]; //!
//   TH2F*		fHistDphiDeta[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt]; //!
//   TH2F*		fHistDphiDetaMix[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt]; //!
  TH2F*		fHistTracksEtaTrigVsEtaAssoc[fNMaxBinsCentrality]; //!
  TH2F*		fHistTracksEtaTrigVsEtaAssocMixed[fNMaxBinsCentrality]; //!
  TH1F*		fHistSingleHadronsPt[fNMaxBinsCentrality]; //!
//   TH1F*		fHistSingleHadronsPt_Mixed[fNMaxBinsCentrality]; //!
  TH2F*		fHistSingleHadronsEtaPt[fNMaxBinsCentrality]; //!
//   TH2F*		fHistSingleHadronsEtaPt_Mixed[fNMaxBinsCentrality]; //!

  //____ MC
  TH1F*		fHistMcGenerator; //!
  TH1F*		fHist_HijingBg; //!
  TH1F*		fHistNumOfPartPerEvt; //!
  TH1F*		fHistMcStats; //!
  TH1F*		fHistMcAllPt; //!
  TH1F*		fHistMcAllPt_Hijing; //!
  TH1F*		fHistMcAllPt_Dec; //!
  TH1F*		fHistMcAllPt_Inj; //!
  TH1F*		fHistMcAllEta_NotHijing; //!
  TH1F*		fHistMcAllEta_Hijing; //!
  TH1F*		fHistMcAllEta; //!

  //____ Correlation
//   AliTHn* fHistCorrSingle;
//   AliTHn* fHistCorrPairSame;
//   AliTHn* fHistCorrPairMixed;
//   THnSparseD* fHistCorrSingle;
  THnSparseD* fHistCorrPair[2]; //!
  
//   THnSparseD* fHistPoolMgrQA;


  TH2F*		fHistControlConvResoncances; //!
  
  TH1D*		fHistTriggerWeighting; //!
  TAxis*		fTriggerWeighting;

  TH2F*		fHistTwoTrackDistancePt[fNMaxBinsPt][2];
  
  TH2F*		fHistHBTbefore; //! Delta Eta vs. Delta Phi before HBT inspired cuts
  TH2F*		fHistHBTafter; //! Delta Eta vs. Delta Phi after HBT inspired cuts
//   TH2F*		fHistConversionbefore; //! Delta Eta vs. Delta Phi before Conversion cuts
//   TH2F*		fHistConversionafter; //! Delta Eta vs. Delta Phi before Conversion cuts
//   TH2F*		fHistPsiMinusPhi;//! psi - phi QA histogram
//   TH3D*		fHistResonancesBefore; //! 3D histogram (Deta,Dphi,Invmass) before resonance cuts
//   TH3D*		fHistResonancesRho;    //! 3D histogram (Deta,Dphi,Invmass) after removing rho 
//   TH3D*		fHistResonancesK0;     //! 3D histogram (Deta,Dphi,Invmass) after removing rho, K0 
//   TH3D*		fHistResonancesLambda; //! 3D histogram (Deta,Dphi,Invmass) after removing rho, K0, and Lambda
//   TH3D*		fHistQbefore; //! Delta Eta vs. Delta Phi before cut on momentum difference
//   TH3D*		fHistQafter; //! Delta Eta vs. Delta Phi after cut on momentum difference


  TH2F* 		fHistoNSigma[fNMaxBinsCentrality]; //!
  Double_t	nsigmaITS[fPartNSpeciesQA];
  Double_t	nsigmaTPC[fPartNSpeciesQA];
  Double_t	nsigmaTOF[fPartNSpeciesQA];
  Double_t	nsigmaHMPID[fPartNSpeciesQA];

  Double_t 	fnsigmas[fPartNSpeciesQA][kMyNSigmaPIDType]; // nsigma values

  TH2F*		fHistTPCdEdx[fNMaxBinsCentrality];		//! TPC dEdx
  TH2F*		fHistTOFbeta[fNMaxBinsCentrality];		//! TOF beta
  TH2F*		fHistTPCdEdx_selected[fNMaxBinsCentrality];	//! TPC dEdx
  TH2F*		fHistTOFbeta_selected[fNMaxBinsCentrality];	//! TOF beta

  TH2F*		fHistNSigmaTPC[fNMaxBinsCentrality][AliPID::kSPECIES];	//! nsigma TPC
  TH2F*		fHistNSigmaTOF[fNMaxBinsCentrality][AliPID::kSPECIES];	//! nsigma TOF
 
  ClassDef(AliAnalysisTaskPidPidCorrelations, 1);
};

//_____ Reduced Tracks -- contains only quantities requires for this analysis to reduce memory consumption for event mixing
class AliPidPidCorrelationReducedTrack : public AliVParticle // TObject
{
  public:
    AliPidPidCorrelationReducedTrack(Int_t partID, Double_t eta, Double_t phi, Double_t pt, Short_t charge)
    : fParticleIDReduced(partID), fEtaReduced(eta), fPhiReduced(phi), fPtReduced(pt), fChargeReduced(charge)
    {
    }
    ~AliPidPidCorrelationReducedTrack() {}
   
  // AliVParticle functions
  virtual Double_t	Px() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t	Py() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t	Pz() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t	Pt() const { return fPtReduced; }
  virtual Double_t	P() const	{ AliFatal("Not implemented"); return 0; }
  virtual Bool_t		PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
  virtual Double_t 	Xv() const	{ AliFatal("Not implemented"); return 0; }
  virtual Double_t 	Yv() const	{ AliFatal("Not implemented"); return 0; }
  virtual Double_t 	Zv() const	{ AliFatal("Not implemented"); return 0; }
  virtual Bool_t   	XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
  virtual Double_t 	OneOverPt()	const	{ AliFatal("Not implemented"); return 0; }
  virtual Double_t	Phi() const	{ return fPhiReduced; }
  virtual Double_t 	Theta() const	{ AliFatal("Not implemented"); return 0; }
  virtual Double_t 	E() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t 	M() const	{ AliFatal("Not implemented"); return 0; }
  virtual Double_t 	Eta() const	{ return fEtaReduced; }
  virtual Double_t 	Y() const	{ AliFatal("Not implemented"); return 0; }
  virtual Short_t 	Charge() const	{ return fChargeReduced; }
    
  // void Print() { Printf(Form("Reduced track, eta: %lf phi: %lf pt: %lf p: %lf", fEtaReduced, fPhiReduced, fPtReduced, fPReduced)); }

  //________ PID
  Int_t				GetMyPartID() const { return fParticleIDReduced; }
  virtual Bool_t			IsEqual(const TObject* obj) const { return (obj->GetUniqueID() == GetUniqueID()); }
  virtual Int_t			GetLabel() const { AliFatal("Not implemented"); return 0; }
  virtual Int_t			PdgCode() const	{ AliFatal("Not implemented"); return 0; }
  virtual const Double_t*	PID() const { AliFatal("Not implemented"); return 0; }
    
  private:
    
  Int_t		fParticleIDReduced;	// particle ID 
  Double_t	fEtaReduced;    		// eta
  Double_t	fPhiReduced;     	// phi
  Double_t	fPtReduced;      	// pT
  Short_t	fChargeReduced;  	// charge

  ClassDef(AliPidPidCorrelationReducedTrack, 1); // reduced track class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};

#endif
