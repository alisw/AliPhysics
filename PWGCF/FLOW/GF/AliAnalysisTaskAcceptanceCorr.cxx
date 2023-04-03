/* -------------------------------------------
 * Maintainer: Mingrui Zhao
 */
#include "AliAnalysisTaskAcceptanceCorr.h"
#include "AliGFWCuts.h"
#include "AliGFWNFCuts.h"
#include "AliGFWWeights.h"
#include "CorrelationCalculator.h"

#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TComplex.h>
#include <TBits.h>
// AliRoot includes
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TComplex.h>

// AliRoot includes
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODITSsaTrackCuts.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"

// STL includes
#include <iostream>
using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskAcceptanceCorr)
	// ---------------------------------------------------------------------------------
	AliAnalysisTaskAcceptanceCorr::AliAnalysisTaskAcceptanceCorr():
		AliAnalysisTaskSE(),
		fEventCuts(),
		fGFWSelection(NULL),
		fGFWSelection15o(NULL),
		fAOD(0),
		fEtaCut(0.8),
		fMinPt(0.2),
		fMaxPt(3.0),
		fTrigger(0),
		fAliTrigger(0),
		fNtrksName("Mult"),
		//....
		fPeriod("LHC15o"),
		fCurrSystFlag(0),
		fAddTPCPileupCuts(false),
		fESDvsTPConlyLinearCut(15000.),
		fTPCchi2perCluster(4.0),
		fUseAdditionalDCACut(false),

		fListOfObjects(0),

		hEventCount(0),
		hMult(0),
		fVtxAfterCuts(0),
		fCentralityDis(0),
		fUseFMDcut(kTRUE),

		fFMDcutapar0(1.64755),
		fFMDcutapar1(119.602),
		fFMDcutcpar0(2.73426),
		fFMDcutcpar1(150.31),
		fFMDAacceptanceCutLower(1.8),
		fFMDAacceptanceCutUpper(4.8),
		fFMDCacceptanceCutLower(-3.2),
		fFMDCacceptanceCutUpper(-1.8)
{
}
//______________________________________________________________________________
AliAnalysisTaskAcceptanceCorr::AliAnalysisTaskAcceptanceCorr(const char *name):
  AliAnalysisTaskSE(name),
  fEventCuts(),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fAOD(0),
  fEtaCut(0.8),
  fMinPt(0.2),
  fMaxPt(3.0),
  fTrigger(0),
  fAliTrigger(0),
  fNtrksName("Mult"),
  //....
  fPeriod("LHC15o"),
  fCurrSystFlag(0),
  fAddTPCPileupCuts(false),
  fESDvsTPConlyLinearCut(15000.),
  fTPCchi2perCluster(4.0),
  fUseAdditionalDCACut(false),

	fListOfObjects(0),
	
	hEventCount(0),
	hMult(0),
	fVtxAfterCuts(0),
	fCentralityDis(0),
	fUseFMDcut(kTRUE),
	fFMDcutapar0(1.64755),
	fFMDcutapar1(119.602),
	fFMDcutcpar0(2.73426),
	fFMDcutcpar1(150.31),
  fFMDAacceptanceCutLower(1.8),
  fFMDAacceptanceCutUpper(4.8),
  fFMDCacceptanceCutLower(-3.2),
  fFMDCacceptanceCutUpper(-1.8)
{

	// Output slot #1 writes into a TList
	DefineOutput(1, TList::Class());
}


// ---------------------------------------------------------------------------------
AliAnalysisTaskAcceptanceCorr::~AliAnalysisTaskAcceptanceCorr() {
	// Destructor
}

// ---------------------------------------------------------------------------------
void AliAnalysisTaskAcceptanceCorr::UserCreateOutputObjects() {
	cout << "UserCreateOutputObjects" << endl;
	// Create output objects
	fListOfObjects = new TList();
	fListOfObjects->SetOwner();

	// Setting for AliEventCuts:
	fEventCuts.AddQAplotsToList(fListOfObjects);

	if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n") ) {
		// Only for LHC15o pass1
		fGFWSelection15o = new AliGFWNFCuts();
		fGFWSelection15o->PrintSetup();
	} else {
		fGFWSelection = new AliGFWCuts();
		fGFWSelection->PrintSetup();
	}

	hEventCount = new TH1D("hEventCount", "; centrality;;", 1, 0, 1);
	fListOfObjects->Add(hEventCount);

  Double_t binning_eta[] = {-4, -3.95, -3.9, -3.85, -3.8, -3.75, -3.7, -3.65, -3.6, -3.55, -3.5, -3.45, -3.4, -3.35, -3.3, -3.25, -3.2, -3.15, -3.1, -3.05,
    -3, -2.95, -2.9, -2.85, -2.8, -2.75, -2.7, -2.65, -2.6, -2.55, -2.5, -2.45, -2.4, -2.35, -2.3, -2.25, -2.2, -2.15, -2.1, -2.05,
    -2, -1.95, -1.9, -1.85, -1.8, -1.75, -1.7, -1.65, -1.6, -1.55, -1.5, -1.45, -1.4, -1.35, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05,
    -1, -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05,
    0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
    1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95,
    2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95,
    3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95,
    4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95,
    5, 5.05, 5.1, 5.15, 5.2, 5.25, 5.3, 5.35, 5.4, 5.45, 5.5, 5.55, 5.6, 5.65, 5.7, 5.75, 5.8, 5.85, 5.9, 5.95, 6};

  hEtaPhiDist = new TH2D("hEtaPhi", "#eta-#phi distribution", 200, binning_eta, 80, 0, 2 * TMath::Pi());
  fListOfObjects->Add(hEtaPhiDist);
	PostData(1, fListOfObjects);
}

// ---------------------------------------------------------------------------------
void AliAnalysisTaskAcceptanceCorr::NotifyRun() {
	if (fAddTPCPileupCuts) {
		Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
		fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
		fEventCuts.fESDvsTPConlyLinearCut[0] = fESDvsTPConlyLinearCut;
	}
}

// ---------------------------------------------------------------------------------
void AliAnalysisTaskAcceptanceCorr::UserExec(Option_t *) {
	// Check if it can pass the trigger
	//..apply physics selection
	UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	Bool_t isTrigselected = false;
	if (fTrigger == 0) {
		isTrigselected = fSelectMask&AliVEvent::kINT7;
		fAliTrigger = AliVEvent::kINT7;
	} else if (fTrigger == 1) {
		isTrigselected = fSelectMask&AliVEvent::kHighMultV0;
		fAliTrigger = AliVEvent::kHighMultV0;
	}
	if(isTrigselected == false) return;

	//..check if I have AOD
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
	if(!fAOD) {
		Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
		return;
	}

	// Check if it passed the standard AOD selection
	if (!AcceptAOD(fAOD) ) {
		PostData(1, fListOfObjects);
		return;
	}
	hEventCount->Fill("after fEventCuts", 1.);
  // cout << "Fill after event selection" << endl;

	if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n") ) { // Only for LHC15o pass1
		fGFWSelection15o->ResetCuts();
	} else {
		fGFWSelection->ResetCuts();
	}
	//..filling Vz distribution
	AliVVertex *vtx = fAOD->GetPrimaryVertex();
	float fVtxZ = vtx->GetZ();
	fPVz = fVtxZ;

	if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n") ) { // Only for LHC15o pass1
		if (!fGFWSelection15o->AcceptVertex(fAOD)) {
			PostData(1, fListOfObjects);
			return;
		}
	} else {
		if (!fGFWSelection->AcceptVertex(fAOD)) {
			PostData(1, fListOfObjects);
			return;
		}
	}

	AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
	if(!multSelection) { return; }
	fCentrality = multSelection->GetMultiplicityPercentile("V0M");

	// Here we calcuate the multiplicity distribution

	PrepareTPCFMDTracks();

	PostData(1, fListOfObjects);
	return;
}

//________________________________________________________________________
void AliAnalysisTaskAcceptanceCorr::Terminate(Option_t *)
{
}

Bool_t AliAnalysisTaskAcceptanceCorr::PrepareTPCFMDTracks() {


  // cout << "Preparing tracks" << endl;
  Int_t nTracks = fInputEvent->GetNumberOfTracks();
  //..for DCA
  double pos[3], vz, vx, vy;
  vz = fInputEvent->GetPrimaryVertex()->GetZ();
  vx = fInputEvent->GetPrimaryVertex()->GetX();
  vy = fInputEvent->GetPrimaryVertex()->GetY();
  double vtxp[3] = {vx, vy, vz};

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(iTrack));

    // Require track to be existing and pass the track selection
    if (!track) continue;

    track->GetXYZ(pos);
    if (!AcceptAODTrack(track, pos,vtxp)) continue;

    hEtaPhiDist->Fill(track->Eta(), track->Phi());
  }

  AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
  const TH2D& d2Ndetadphi = aodForward->GetHistogram();
  int nEta = d2Ndetadphi.GetXaxis()->GetNbins();
  int nPhi = d2Ndetadphi.GetYaxis()->GetNbins();


  for (int iEta = 1; iEta <= nEta; iEta++) {
    double eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta); 
    // cout << iEta << " " << eta - 0.025 << endl;
    for (int iPhi = 1; iPhi <= nPhi; iPhi++) {
      double phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi); 
      // cout << iPhi << " " << phi/TMath::Pi()*80 << endl;


      double etaAccepted = false;
      if ( fFMDAacceptanceCutLower < eta && eta < fFMDAacceptanceCutUpper) etaAccepted = true;
      if ( fFMDCacceptanceCutLower < eta && eta < fFMDCacceptanceCutUpper) etaAccepted = true;
      if (!etaAccepted) continue;

      double mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
      if (mostProbableN > 0) {
        hEtaPhiDist->Fill(eta, phi, mostProbableN);
      }
    }
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskAcceptanceCorr::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp) {
	// Pt cut
	if(mtr->Pt() < fMinPt) return kFALSE;
	if(mtr->Pt() > fMaxPt) return kFALSE;

	// DCA cut
	if(ltrackXYZ && vtxp) {
		mtr->GetXYZ(ltrackXYZ);
		ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
		ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
		ltrackXYZ[2] = abs(ltrackXYZ[2]-vtxp[2]);
	} else return kFALSE; //DCA cut is a must for now

	// Additional cut for TPCchi2perCluster
	if (mtr->GetTPCchi2perCluster()>fTPCchi2perCluster) return kFALSE;

	if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n") ) { // Only for LHC15o pass1 and LHC17n
		return fGFWSelection15o->AcceptTrack(mtr,ltrackXYZ,0,kFALSE);
	} else {
		return fGFWSelection->AcceptTrack(mtr,ltrackXYZ,0,kFALSE);
	}
}

Bool_t AliAnalysisTaskAcceptanceCorr::AcceptAOD(AliAODEvent *inEv) {
	// LHC15i, LHC15l, LHC16, LHC17, LHC18: means: pp sample
	if (fPeriod.EqualTo("LHC15i") ||
			fPeriod.EqualTo("LHC15l") ||
			fPeriod.EqualTo("LHC16Preview") ||
			fPeriod.EqualTo("LHC17Preview") ||
			fPeriod.EqualTo("LHC18Preview") || 
			fPeriod.EqualTo("LHC16") ||
			fPeriod.EqualTo("LHC17") ||
			fPeriod.EqualTo("LHC18") || 
			fPeriod.EqualTo("LHC16ZM") ||
			fPeriod.EqualTo("LHC17ZM") ||
			fPeriod.EqualTo("LHC18ZM") ) {
		fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
	}


	if(!fEventCuts.AcceptEvent(inEv)) return false;

	// Primary vertex
	const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertex());
	if(!vtx || vtx->GetNContributors() < 1)
		return kFALSE;

	// SPD Vertex
	const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertexSPD());
	Double_t dMaxResol = 0.25; // suggested from DPG
	Double_t cov[6] = {0};
	vtxSPD->GetCovarianceMatrix(cov);
	Double_t zRes = TMath::Sqrt(cov[5]);
	if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;

	if (fPeriod.EqualTo("LHC15o") ||
			fPeriod.EqualTo("LHC15o_pass2") ||
			fPeriod.EqualTo("LHC18qr_pass3") ||
			fPeriod.EqualTo("LHC16qt") ||
			fPeriod.EqualTo("LHC17n") ||
			fPeriod.EqualTo("LHC15oKatarina")) {
		// return false;
	} else {
		// if(fAOD->IsPileupFromSPDInMultBins() ) { return false; }

		// AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
		// if (!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return false; }

		// if(!multSelection->GetThisEventIsNotPileup() || !multSelection->GetThisEventIsNotPileupInMultBins() || !multSelection->GetThisEventHasNoInconsistentVertices() || !multSelection->GetThisEventPassesTrackletVsCluster()) { return false; }

		// Int_t nTracksPrim = fAOD->GetPrimaryVertex()->GetNContributors();
		// if(nTracksPrim < 0.5) { return false; }
	}

	// Vertex Z
	const Double_t aodVtxZ = vtx->GetZ();
	if(TMath::Abs(aodVtxZ) > 10) return kFALSE;


	// FMD cut should be checked **before** the Preparation of FMD tracks
	if(fUseFMDcut){

		double nFMD_fwd_hits = 0;
		double nFMD_bwd_hits = 0;
		AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
		const TH2D& d2Ndetadphi = aodForward->GetHistogram();
		int nEta = d2Ndetadphi.GetXaxis()->GetNbins();
		int nPhi = d2Ndetadphi.GetYaxis()->GetNbins();

		for (int iEta = 1; iEta <= nEta; iEta++) {
			for (int iPhi = 1; iPhi <= nPhi; iPhi++) {
				double eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta); 
				double phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi); 

				double etaAccepted = false;
				if ( fFMDAacceptanceCutLower < eta && eta < fFMDAacceptanceCutUpper) etaAccepted = true;
				if ( fFMDCacceptanceCutLower < eta && eta < fFMDCacceptanceCutUpper) etaAccepted = true;
				if (!etaAccepted) continue;

				double mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
				if (mostProbableN > 0) {
					if (eta > 0) nFMD_fwd_hits += mostProbableN;
					else nFMD_bwd_hits += mostProbableN;
				}
			}
		}
		if(nFMD_fwd_hits==0 || nFMD_bwd_hits==0) {
			return kFALSE;
		}

		AliAODVZERO *fvzero = fAOD->GetVZEROData();
		if(!fvzero) { AliError("Problem with VZEROData, terminating!"); return kFALSE; }
		Float_t nV0A_hits = fvzero->GetMTotV0A();
		Float_t nV0C_hits = fvzero->GetMTotV0C();

    // cout << nV0A_hits << endl;
		// cout << nV0C_hits << endl;
		// cout << nFMD_fwd_hits << endl;
		// cout << nFMD_bwd_hits << endl;

		if((nV0A_hits<(fFMDcutapar0*nFMD_fwd_hits-fFMDcutapar1)) || (nV0C_hits<(fFMDcutcpar0*nFMD_bwd_hits-fFMDcutcpar1))){
			return kFALSE;
		}
		// hFMDAvsV0->Fill(nFMD_fwd_hits, nV0A_hits);
		// hFMDCvsV0->Fill(nFMD_bwd_hits, nV0C_hits);
	}

	return kTRUE;
}
