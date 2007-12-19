/**************************************************************************
 * Author: Boris Hippolyte.                                               *
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

//-----------------------------------------------------------------
//                 AliAnalysisTaskESDCheckV0 class
//            This task is for QAing the dE/dx from the ESD
//              Origin: B.H. Nov2007, hippolyt@in2p3.fr
//-----------------------------------------------------------------

//-----------------------------------------------------------------
//    ---> Next upgrades
//      1) extent the use of fMidPseudoRapidityFlag to primaries
//      2) discuss the GetSigmaToVertex(AliESDtrack*)
//      3) decide the kind of refit to be used
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVector3.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliESDVertex.h"
#include "AliESDv0.h"

#include "AliAnalysisTaskESDDedx.h"

ClassImp(AliAnalysisTaskESDDedx)

//________________________________________________________________________
  AliAnalysisTaskESDDedx::AliAnalysisTaskESDDedx(const char    *rName,
						 const Bool_t   rAllConstrainedFlag,
						 const Bool_t   rMidPseudoRapidityFlag,
						 const Bool_t   rSelTrackRemoveKink,
						 const Bool_t   rSelTrackWithOnTheFlyV0,
						 const Int_t    rSelTrackMinClustersTPC,
						 const Int_t    rSelTrackMinClustersITS,
						 const Float_t  rSelTrackMaxChi2PerClusterTPC,
						 const Float_t  rSelTrackMaxChi2PerClusterITS,
						 const Double_t rSelTrackMaxCov11,
						 const Double_t rSelTrackMaxCov22,
						 const Double_t rSelTrackMaxCov33,
						 const Double_t rSelTrackMaxCov44,
						 const Double_t rSelTrackMaxCov55,
						 const Double_t rSelV0MaxDcaDaughters,
						 const Double_t rSelV0MinDecayLength)

    : AliAnalysisTask(rName, ""), fESD(0), fListHist(), fHistPtot(0),

    fHistMultiplicity(0), fHistTPCDedxVsMomentum(0), fHistITSDedxVsMomentum(0),
    fHistMassK0(0), fHistMassLambda(0), fHistMassAntiLambda(0),
    fHistTPCDedxVsMomPosK0(0), fHistTPCDedxVsMomNegK0(0),
    fHistTPCDedxVsMomPosLambda(0), fHistTPCDedxVsMomNegLambda(0),
    fHistTPCDedxVsMomPosAntiLambda(0), fHistTPCDedxVsMomNegAntiLambda(0),
    fHistDiffInOutMomentum(0), fHistDiffPrimOutMomentum(0),
    fHistDiffPrimMeanMomentum(0), fHistPercPrimMeanMomentum(0), fHistPrimEta(0),
    fHistPercPrimMeanMomentumVsEta(0), fHistPercPrimMeanMomentumVsPrim(0),

    fHistMultiplicityCuts(0), fHistTPCDedxVsMomentumCuts(0), fHistITSDedxVsMomentumCuts(0),
    fHistMassK0Cuts(0), fHistMassLambdaCuts(0), fHistMassAntiLambdaCuts(0), 
    fHistTPCDedxVsMomPosK0Cuts(0), fHistTPCDedxVsMomNegK0Cuts(0), 
    fHistTPCDedxVsMomPosLambdaCuts(0), fHistTPCDedxVsMomNegLambdaCuts(0),
    fHistTPCDedxVsMomPosAntiLambdaCuts(0), fHistTPCDedxVsMomNegAntiLambdaCuts(0),
    fHistDiffInOutMomentumCuts(0), fHistDiffPrimOutMomentumCuts(0),
    fHistDiffPrimMeanMomentumCuts(0), fHistPercPrimMeanMomentumCuts(0), fHistPrimEtaCuts(0),
    fHistPercPrimMeanMomentumVsEtaCuts(0), fHistPercPrimMeanMomentumVsPrimCuts(0),

    fAllConstrainedFlag(rAllConstrainedFlag), fMidPseudoRapidityFlag(rMidPseudoRapidityFlag),

    fSelTrackRemoveKink(rSelTrackRemoveKink), fSelTrackWithOnTheFlyV0(rSelTrackWithOnTheFlyV0),
    fSelTrackMinClustersTPC(rSelTrackMinClustersTPC), fSelTrackMinClustersITS(rSelTrackMinClustersITS),
    fSelTrackMaxChi2PerClusterTPC(rSelTrackMaxChi2PerClusterTPC), fSelTrackMaxChi2PerClusterITS(rSelTrackMaxChi2PerClusterITS),
    fSelTrackMaxCov11(rSelTrackMaxCov11), fSelTrackMaxCov22(rSelTrackMaxCov22), fSelTrackMaxCov33(rSelTrackMaxCov33),
    fSelTrackMaxCov44(rSelTrackMaxCov44), fSelTrackMaxCov55(rSelTrackMaxCov55),
    fSelV0MaxDcaDaughters(rSelV0MaxDcaDaughters), fSelV0MinDecayLength(rSelV0MinDecayLength)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TH1F::Class());
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskESDDedx::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("fTracks.*", kTRUE);
    tree->SetBranchStatus("fV0s.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskESDDedx::CreateOutputObjects()
{
  // Create histograms
  // Called once
  fListHist = new TList();
  if (!fHistPtot) {
    fHistPtot = new TH1F("fHistPtot", "P_{tot} distribution;P_{tot} (GeV/c);dN/dP_{tot} (c/GeV)", 15, 0.1, 3.1);
    fHistPtot->SetMarkerStyle(kFullCircle);
    fListHist->Add(fHistPtot);
  }
  if (!fHistMultiplicity) {
    fHistMultiplicity = new TH1F("fHistMultiplicity", "Multiplicity distribution;Number of tracks;Events", 250, 0, 250);
    fHistMultiplicity->SetMarkerStyle(kFullCircle);
    fListHist->Add(fHistMultiplicity);
  }
  if (!fHistTPCDedxVsMomentum) {
    fHistTPCDedxVsMomentum = new TH2F("h2TPCDedxVsMomentum","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomentum);
  }
  if (!fHistITSDedxVsMomentum) {
    fHistITSDedxVsMomentum = new TH2F("h2ITSDedxVsMomentum","Bethe-Bloch Distribution for ITS;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistITSDedxVsMomentum);
  }
  if (!fHistMassK0) {
    fHistMassK0 = new TH1F("h1MassK0","K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts",100,0.4,0.6);
    fListHist->Add(fHistMassK0);
  }
  if (!fHistMassLambda) {
    fHistMassLambda = new TH1F("h1MassLambda","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts",75,1.05,1.2);
    fListHist->Add(fHistMassLambda);
  }
  if (!fHistMassAntiLambda) {
    fHistMassAntiLambda = new TH1F("h1MassAntiLambda","#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts",75,1.05,1.2);
    fListHist->Add(fHistMassAntiLambda);
  }
  if (!fHistTPCDedxVsMomPosK0) {
    fHistTPCDedxVsMomPosK0 = new TH2F("h2TPCDedxVsMomPosK0","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomPosK0);
  }
  if (!fHistTPCDedxVsMomNegK0) {
    fHistTPCDedxVsMomNegK0 = new TH2F("h2TPCDedxVsMomNegK0","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomNegK0);
  }
  if (!fHistTPCDedxVsMomPosLambda) {
    fHistTPCDedxVsMomPosLambda = new TH2F("h2TPCDedxVsMomPosLambda","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomPosLambda);
  }
  if (!fHistTPCDedxVsMomNegLambda) {
    fHistTPCDedxVsMomNegLambda = new TH2F("h2TPCDedxVsMomNegLambda","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomNegLambda);
  }
  if (!fHistTPCDedxVsMomPosAntiLambda) {
    fHistTPCDedxVsMomPosAntiLambda = new TH2F("h2TPCDedxVsMomPosAntiLambda","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomPosAntiLambda);
  }
  if (!fHistTPCDedxVsMomNegAntiLambda) {
    fHistTPCDedxVsMomNegAntiLambda = new TH2F("h2TPCDedxVsMomNegAntiLambda","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomNegAntiLambda);
  }
  if (!fHistDiffInOutMomentum) {
    fHistDiffInOutMomentum = new TH1F("h1DiffInOutMomentum","Momentum Difference In-Out;Momentum (GeV/c);Counts",200,-0.1,0.1); 
    fListHist->Add(fHistDiffInOutMomentum);
  }
  if (!fHistDiffPrimOutMomentum) {
    fHistDiffPrimOutMomentum = new TH1F("h1DiffPrimOutMomentum","Momentum Difference Prim-Out;Momentum (GeV/c);Counts",200,-0.1,0.1); 
    fListHist->Add(fHistDiffPrimOutMomentum);
  }
  if (!fHistDiffPrimMeanMomentum) {
    fHistDiffPrimMeanMomentum = new TH1F("h1DiffPrimMeanMomentum","Momentum Difference Prim-Mean;Momentum (GeV/c);Counts",200,-0.1,0.1); 
    fListHist->Add(fHistDiffPrimMeanMomentum);
  }
  if (!fHistPercPrimMeanMomentum) {
    fHistPercPrimMeanMomentum = new TH1F("h1PercPrimMeanMomentum","Momentum Percentage (Prim-Mean)/Prim;(Primary-Mean)/Primary (%);Counts",200,-10,10); 
    fListHist->Add(fHistPercPrimMeanMomentum);
  }
  if (!fHistPrimEta) {
    fHistPrimEta = new TH1F("h1PrimEta","Pseudorapidity Distribution (Primaries);#eta;Counts",300,-1.5,1.5); 
    fListHist->Add(fHistPrimEta);
  }
  if (!fHistPercPrimMeanMomentumVsEta) {
    fHistPercPrimMeanMomentumVsEta = new TH2F("h2PercPrimMeanMomentumVsEta","Momentum Percentage (Prim-Mean)/Prim vs #eta;#eta;(Primary-Mean)/Primary (%)",100,-1.5,1.5,100,-10,10); 
    fListHist->Add(fHistPercPrimMeanMomentumVsEta);
  }
  if (!fHistPercPrimMeanMomentumVsPrim) {
    fHistPercPrimMeanMomentumVsPrim = new TH2F("h2PercPrimMeanMomentumVsPrim","Momentum Percentage (Prim-Mean)/Prim vs Prim;Momentum (GeV/c);(Primary-Mean)/Primary (%)",1500,0,15,100,-10,10); 
    fListHist->Add(fHistPercPrimMeanMomentumVsPrim);
  }
  // Histograms after selections
  if (!fHistMultiplicityCuts) {
    fHistMultiplicityCuts = new TH1F("h1HistMultiplicityCuts","Number of Good tracks;Number of tracks;Event",201,-0.5,200.5);
    fListHist->Add(fHistMultiplicityCuts);
  }
  if (!fHistTPCDedxVsMomentumCuts) {
    fHistTPCDedxVsMomentumCuts = new TH2F("h2TPCDedxVsMomentumCuts","Bethe-Bloch Distribution for TPC w/cuts;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomentumCuts);
  }
  if (!fHistITSDedxVsMomentumCuts) {
    fHistITSDedxVsMomentumCuts = new TH2F("h2ITSDedxVsMomentumCuts","Bethe-Bloch Distribution for ITS w/cuts;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistITSDedxVsMomentumCuts);
  }
  if (!fHistMassK0Cuts) {
    fHistMassK0Cuts = new TH1F("h1MassK0Cuts","K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts",100,0.4,0.6);
    fListHist->Add(fHistMassK0Cuts);
  }
  if (!fHistMassLambdaCuts) {
    fHistMassLambdaCuts = new TH1F("h1MassLambdaCuts","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts",75,1.05,1.2);
    fListHist->Add(fHistMassLambdaCuts);
  }
  if (!fHistMassAntiLambdaCuts) {
    fHistMassAntiLambdaCuts = new TH1F("h1MassAntiLambdaCuts","#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts",75,1.05,1.2);
    fListHist->Add(fHistMassAntiLambdaCuts);
  }
  if (!fHistTPCDedxVsMomPosK0Cuts) {
    fHistTPCDedxVsMomPosK0Cuts = new TH2F("h2TPCDedxVsMomPosK0Cuts","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomPosK0Cuts);
  }
  if (!fHistTPCDedxVsMomNegK0Cuts) {
    fHistTPCDedxVsMomNegK0Cuts = new TH2F("h2TPCDedxVsMomNegK0Cuts","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomNegK0Cuts);
  }
  if (!fHistTPCDedxVsMomPosLambdaCuts) {
    fHistTPCDedxVsMomPosLambdaCuts = new TH2F("h2TPCDedxVsMomPosLambdaCuts","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomPosLambdaCuts);
  }
  if (!fHistTPCDedxVsMomNegLambdaCuts) {
    fHistTPCDedxVsMomNegLambdaCuts = new TH2F("h2TPCDedxVsMomNegLambdaCuts","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomNegLambdaCuts);
  }
  if (!fHistTPCDedxVsMomPosAntiLambdaCuts) {
    fHistTPCDedxVsMomPosAntiLambdaCuts = new TH2F("h2TPCDedxVsMomPosAntiLambdaCuts","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomPosAntiLambdaCuts);
  }
  if (!fHistTPCDedxVsMomNegAntiLambdaCuts) {
    fHistTPCDedxVsMomNegAntiLambdaCuts = new TH2F("h2TPCDedxVsMomNegAntiLambdaCuts","Bethe-Bloch Distribution for TPC;Momentum (GeV/c);dE/dx (A.U.)",1500,0,15,100,0,100);
    fListHist->Add(fHistTPCDedxVsMomNegAntiLambdaCuts);
  }
  if (!fHistDiffInOutMomentumCuts) {
    fHistDiffInOutMomentumCuts = new TH1F("h1DiffInOutMomentumCuts","Momentum Difference In-Out;Momentum (GeV/c);Counts",200,-0.1,0.1); 
    fListHist->Add(fHistDiffInOutMomentumCuts);
  }
  if (!fHistDiffPrimOutMomentumCuts) {
    fHistDiffPrimOutMomentumCuts = new TH1F("h1DiffPrimOutMomentumCuts","Momentum Difference Prim-Out;Momentum (GeV/c);Counts",200,-0.1,0.1); 
    fListHist->Add(fHistDiffPrimOutMomentumCuts);
  }
  if (!fHistDiffPrimMeanMomentumCuts) {
    fHistDiffPrimMeanMomentumCuts = new TH1F("h1DiffPrimMeanMomentumCuts","Momentum Difference Prim-Mean;Momentum (GeV/c);Counts",200,-0.1,0.1); 
    fListHist->Add(fHistDiffPrimMeanMomentumCuts);
  }
  if (!fHistPercPrimMeanMomentumCuts) {
    fHistPercPrimMeanMomentumCuts = new TH1F("h1PercPrimMeanMomentumCuts","Momentum Percentage (Prim-Mean)/Prim;(Primary-Mean)/Primary (%);Counts",200,-10,10); 
    fListHist->Add(fHistPercPrimMeanMomentumCuts);
  }
  if (!fHistPrimEtaCuts) {
    fHistPrimEtaCuts = new TH1F("h1PrimEtaCuts","Pseudorapidity Distribution (Primaries);#eta;Counts",300,-1.5,1.5); 
    fListHist->Add(fHistPrimEtaCuts);
  }
  if (!fHistPercPrimMeanMomentumVsEtaCuts) {
    fHistPercPrimMeanMomentumVsEtaCuts = new TH2F("h2PercPrimMeanMomentumVsEtaCuts","Momentum Percentage (Prim-Mean)/Prim vs #eta;#eta;(Primary-Mean)/Primary (%)",100,-1.5,1.5,100,-10,10); 
    fListHist->Add(fHistPercPrimMeanMomentumVsEtaCuts);
  }
  if (!fHistPercPrimMeanMomentumVsPrimCuts) {
    fHistPercPrimMeanMomentumVsPrimCuts = new TH2F("h2PercPrimMeanMomentumVsPrimCuts","Momentum Percentage (Prim-Mean)/Prim vs Prim;Momentum (GeV/c);(Primary-Mean)/Primary (%)",1500,0,15,100,-10,10); 
    fListHist->Add(fHistPercPrimMeanMomentumVsPrimCuts);
  }
}

//________________________________________________________________________
void AliAnalysisTaskESDDedx::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  Int_t    nGetTracks        = fESD->GetNumberOfTracks();
  //  Printf("There are %d tracks in this event", nGetTracks);

  Int_t    nAllTracks        = 0, nGoodTracks       = 0;
  Double_t lTPCDedx          = 0, lITSDedx          = 0;

  // Track loop
  for (Int_t iTracks = 0; iTracks < nGetTracks; iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not retrieve track %d", iTracks);
      continue;
    }
    fHistPtot->Fill(track->P()); // Fill a ptot spectrum before any selection

    Double_t lPrimvtxMomemtum[3];
    Bool_t lIsTrackConstrained = track->GetConstrainedPxPyPz(lPrimvtxMomemtum);
    if (fAllConstrainedFlag && (!lIsTrackConstrained)) continue; // Constrained or not to Prim Vertex.

    const AliExternalTrackParam *lInnerETP = track->GetInnerParam();
    const AliExternalTrackParam *lOuterETP = track->GetOuterParam();
    if ((!lInnerETP)||(!lOuterETP)) continue;  // No inner nor outer params for this track.

    Double_t lInnerMomentum = 0, lOuterMomentum = 0, lMeanMomentum = 0;
    if (lInnerETP) lInnerMomentum = lInnerETP->GetP(); // The "inner" momentum.
    if (lOuterETP) lOuterMomentum = lOuterETP->GetP(); // The "outer" momentum.
    if (lInnerETP&&lOuterETP) lMeanMomentum = (lInnerMomentum+lOuterMomentum)/2.; // The "mean" momentum 

    lTPCDedx = track->GetTPCsignal();
    lITSDedx = track->GetITSsignal();
    // This needs to be discussed and implemented
    // nSigmaToVertex = GetSigmaToVertex(track);

    Float_t lPrimPseudoRap = 0, lPrimMomentum = 0;
    if (lIsTrackConstrained){
      TVector3 pPrimvtx(lPrimvtxMomemtum);
      lPrimPseudoRap = pPrimvtx.Eta();
      lPrimMomentum = pPrimvtx.Mag();
    }

    if (
	((!fAllConstrainedFlag)||(TMath::Abs(lPrimPseudoRap) < 0.5)) &&
	(lTPCDedx)
	)
      // fAllConstrained = 0 -> !fAllConstrained = 1 so ok anyway
      // fAllConstrained = 1 -> !fAllConstrained = 0 so it depends on eta condition 
      {
	nAllTracks++;
	fHistTPCDedxVsMomentum->Fill(lInnerMomentum,lTPCDedx);
	fHistDiffInOutMomentum->Fill(lInnerMomentum-lOuterMomentum);
	if (lIsTrackConstrained){
	  fHistITSDedxVsMomentum->Fill(lPrimMomentum,lITSDedx);
 	  fHistDiffPrimOutMomentum->Fill(lPrimMomentum-lOuterMomentum);
 	  fHistDiffPrimMeanMomentum->Fill(lPrimMomentum-lMeanMomentum);
	  if (lPrimMomentum){
 	    fHistPercPrimMeanMomentum->Fill(100.*(lPrimMomentum-lMeanMomentum)/lPrimMomentum);
 	    fHistPercPrimMeanMomentumVsEta->Fill(lPrimPseudoRap,100.*(lPrimMomentum-lMeanMomentum)/lPrimMomentum);
 	    fHistPercPrimMeanMomentumVsPrim->Fill(lPrimMomentum,100.*(lPrimMomentum-lMeanMomentum)/lPrimMomentum);
	  }
 	  fHistPrimEta->Fill(lPrimPseudoRap);
	}
	if(
	   // (!fAllConstrainedFlag || (nSigmaToVertex < 3)) &&
	   (IsAccepted(track))
	   ){
 	  fHistTPCDedxVsMomentumCuts->Fill(lInnerMomentum,lTPCDedx);
 	  fHistDiffInOutMomentumCuts->Fill(lInnerMomentum-lOuterMomentum);
	  if (lIsTrackConstrained){
 	    fHistITSDedxVsMomentumCuts->Fill(lPrimMomentum,lITSDedx);
 	    fHistDiffPrimOutMomentumCuts->Fill(lPrimMomentum-lOuterMomentum);
 	    fHistDiffPrimMeanMomentumCuts->Fill(lPrimMomentum-lMeanMomentum);
	    if (lPrimMomentum){
 	      fHistPercPrimMeanMomentumCuts->Fill(100.*(lPrimMomentum-lMeanMomentum)/lPrimMomentum);
 	      fHistPercPrimMeanMomentumVsEtaCuts->Fill(lPrimPseudoRap,100.*(lPrimMomentum-lMeanMomentum)/lPrimMomentum);
 	      fHistPercPrimMeanMomentumVsPrimCuts->Fill(lPrimMomentum,100.*(lPrimMomentum-lMeanMomentum)/lPrimMomentum);
	    }
 	    fHistPrimEtaCuts->Fill(lPrimPseudoRap);
	  }
	  nGoodTracks++;
	}
      }
  } //track loop 
  fHistMultiplicity->Fill(nGetTracks);
  fHistMultiplicityCuts->Fill(nGoodTracks);

  Int_t nv0s = 0;
  nv0s = fESD->GetNumberOfV0s();

  Int_t    lIndexTrackPos       = 0,  lIndexTrackNeg       = 0;
  Double_t lTPCDedxPos          = 0,  lTPCDedxNeg          = 0;
  Float_t lInvMassK0 = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) 
    {// This is the V0 loop
      AliESDv0 *v0 = fESD->GetV0(iV0);
      if (!v0) continue;
      if (v0->GetOnFlyStatus() != fSelTrackWithOnTheFlyV0) continue;

      const AliESDVertex *lPrimaryVertex = fESD->GetVertex();
      Double_t tPositionPrimaryVertex[3]; lPrimaryVertex->GetXYZ(tPositionPrimaryVertex);
      Double_t v0PositionX = 0, v0PositionY = 0, v0PositionZ = 0; v0->GetXYZ(v0PositionX,v0PositionY,v0PositionZ);
      Double_t v0DecayLength = TMath::Sqrt(
					   (v0PositionX-tPositionPrimaryVertex[0])*(v0PositionX-tPositionPrimaryVertex[0])+
					   (v0PositionY-tPositionPrimaryVertex[1])*(v0PositionY-tPositionPrimaryVertex[1])+
					   (v0PositionZ-tPositionPrimaryVertex[2])*(v0PositionZ-tPositionPrimaryVertex[2])
					   );
      if (v0->GetDcaV0Daughters() > fSelV0MaxDcaDaughters) continue;
      if (v0DecayLength < fSelV0MinDecayLength) continue;

      // Getting invariant mass infos directly from ESD
      v0->ChangeMassHypothesis(310);
      lInvMassK0 = v0->GetEffMass();
      v0->ChangeMassHypothesis(3122);
      lInvMassLambda = v0->GetEffMass();
      v0->ChangeMassHypothesis(-3122);
      lInvMassAntiLambda = v0->GetEffMass();

      // Filling invariant mass histos for all candidates
      fHistMassK0->Fill(lInvMassK0);
      fHistMassLambda->Fill(lInvMassLambda);
      fHistMassAntiLambda->Fill(lInvMassAntiLambda);

      // Accessing the daughter track infos
      lIndexTrackPos = TMath::Abs(v0->GetPindex());
      lIndexTrackNeg = TMath::Abs(v0->GetNindex());
      AliESDtrack* trackPos = fESD->GetTrack(lIndexTrackPos);
      AliESDtrack* trackNeg = fESD->GetTrack(lIndexTrackNeg);
      
      Double_t lPrimvtxMomPos[3];
      Bool_t lIsPosConstrained = trackPos->GetConstrainedPxPyPz(lPrimvtxMomPos);
      if (fAllConstrainedFlag && (!lIsPosConstrained)) continue; // Constrained or not to Prim Vertex.

      Double_t lPrimvtxMomNeg[3];
      Bool_t lIsNegConstrained = trackNeg->GetConstrainedPxPyPz(lPrimvtxMomNeg);
      if (fAllConstrainedFlag && (!lIsNegConstrained)) continue; // Constrained or not to Prim Vertex.

      const AliExternalTrackParam *lInnerPosETP = trackPos->GetInnerParam();
      const AliExternalTrackParam *lInnerNegETP = trackNeg->GetInnerParam();
      if ((!lInnerPosETP)||(!lInnerNegETP)) continue; // No inner params for at least one track.

      Double_t lInnerMomPos = 0, lInnerMomNeg = 0;
      if (lInnerPosETP) lInnerMomPos = lInnerPosETP->GetP();
      if (lInnerNegETP) lInnerMomNeg = lInnerNegETP->GetP();

      lTPCDedxPos = trackPos->GetTPCsignal();
      lTPCDedxNeg = trackNeg->GetTPCsignal();

      Float_t lPrimPseudoRapPos = 0, lPrimPseudoRapNeg = 0;
      if (lIsPosConstrained && lIsNegConstrained){
	TVector3 pPrimPos(lPrimvtxMomPos);
	TVector3 pPrimNeg(lPrimvtxMomNeg);
	lPrimPseudoRapPos = pPrimPos.Eta();
	lPrimPseudoRapNeg = pPrimNeg.Eta();
      }

      if (TMath::Abs(lInvMassK0-0.497)<0.01) {
	fHistMassK0Cuts->Fill(lInvMassK0);
	fHistTPCDedxVsMomPosK0->Fill(lInnerMomPos,lTPCDedxPos);
	fHistTPCDedxVsMomNegK0->Fill(lInnerMomNeg,lTPCDedxNeg);
	if (
	    (!fMidPseudoRapidityFlag) || ( TMath::Abs(lPrimPseudoRapPos) < 0.5) &&
	    (IsAccepted(trackPos))
	    )
	  fHistTPCDedxVsMomPosK0Cuts->Fill(lInnerMomPos,lTPCDedxPos);
	if (
	    (!fMidPseudoRapidityFlag) || ( TMath::Abs(lPrimPseudoRapNeg) < 0.5) &&
	    (IsAccepted(trackNeg))
	    )
	  fHistTPCDedxVsMomNegK0Cuts->Fill(lInnerMomNeg,lTPCDedxNeg);
      }
      if (TMath::Abs(lInvMassLambda-1.115)<0.01) {
	fHistMassLambdaCuts->Fill(lInvMassLambda);
	fHistTPCDedxVsMomPosLambda->Fill(lInnerMomPos,lTPCDedxPos);
	fHistTPCDedxVsMomNegLambda->Fill(lInnerMomNeg,lTPCDedxNeg);
	if (
	    (!fMidPseudoRapidityFlag) || ( TMath::Abs(lPrimPseudoRapPos) < 0.5) &&
	    (IsAccepted(trackPos))
	    )
	  fHistTPCDedxVsMomPosLambdaCuts->Fill(lInnerMomPos,lTPCDedxPos);
	if (
	    (!fMidPseudoRapidityFlag) || ( TMath::Abs(lPrimPseudoRapNeg) < 0.5) &&
	    (IsAccepted(trackNeg))
	    )
	  fHistTPCDedxVsMomNegLambdaCuts->Fill(lInnerMomNeg,lTPCDedxNeg);
      }
      if (TMath::Abs(lInvMassAntiLambda-1.115)<0.01) {
	fHistMassAntiLambdaCuts->Fill(lInvMassAntiLambda);
	fHistTPCDedxVsMomPosAntiLambda->Fill(lInnerMomPos,lTPCDedxPos);
	fHistTPCDedxVsMomNegAntiLambda->Fill(lInnerMomNeg,lTPCDedxNeg);
	if (
	    (!fMidPseudoRapidityFlag) || ( TMath::Abs(lPrimPseudoRapPos) < 0.5) &&
	    (IsAccepted(trackPos))
	    )
	  fHistTPCDedxVsMomPosAntiLambdaCuts->Fill(lInnerMomPos,lTPCDedxPos);
	if (
	    (!fMidPseudoRapidityFlag) || ( TMath::Abs(lPrimPseudoRapNeg) < 0.5) &&
	    (IsAccepted(trackNeg))
	    )
	  fHistTPCDedxVsMomNegAntiLambdaCuts->Fill(lInnerMomNeg,lTPCDedxNeg);
      }
    }//V0 loop 
  
  // Post output data.
  PostData(0, fHistPtot);
  PostData(1, fListHist);
}      

//________________________________________________________________________
void AliAnalysisTaskESDDedx::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  Printf("Reminder for options and selections: fAllConstrainedFlag = %d fMidPseudoRapidityFlag = %d fSelTrackRemoveKink = %d fSelTrackWithOnTheFlyV0 = %d fSelTrackMinClustersTPC = %d fSelTrackMinClustersITS = %d fSelTrackMaxChi2PerClusterTPC = %.2f fSelTrackMaxChi2PerClusterITS = %.2f fSelTrackMaxCov11 = %.2f fSelTrackMaxCov22 = %.2f fSelTrackMaxCov33 = %.2f fSelTrackMaxCov44 = %.2f fSelTrackMaxCov55 = %.2f fSelV0MaxDcaDaughters = %.2f fSelV0MinDecayLength= %.2f",
	 fAllConstrainedFlag,
	 fMidPseudoRapidityFlag,fSelTrackRemoveKink,fSelTrackWithOnTheFlyV0,
	 fSelTrackMinClustersTPC,fSelTrackMinClustersITS,
	 fSelTrackMaxChi2PerClusterTPC,fSelTrackMaxChi2PerClusterITS,
	 fSelTrackMaxCov11,fSelTrackMaxCov22,fSelTrackMaxCov33,fSelTrackMaxCov44,fSelTrackMaxCov55,
	 fSelV0MaxDcaDaughters,fSelV0MinDecayLength);

  fHistPtot = dynamic_cast<TH1F*> (GetOutputData(0));
  if (!fHistPtot) {
    Printf("ERROR: fHistPtot not available");
    return;
  }
   
  TCanvas *cTaskDedx = new TCanvas("AliAnalysisTaskESDDedx","Ptot",10,10,510,510);
  cTaskDedx->cd(1)->SetLogy();
  fHistPtot->DrawCopy("E");
}
//____________________________________________________________________//
Bool_t AliAnalysisTaskESDDedx::IsAccepted(AliESDtrack* track) {
  // Checks if the track is excluded from the cuts

  Int_t   fIdxInt[200];
  Int_t   nClustersITS = track->GetITSclusters(fIdxInt);
  Int_t   nClustersTPC = track->GetTPCclusters(fIdxInt);

  Float_t chi2PerClusterITS = -1;
  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

  Double_t extCov[15];
  track->GetExternalCovariance(extCov);

  if (fSelTrackRemoveKink)
    if (track->GetKinkIndex(0)>0) return kFALSE;
  if (fSelTrackMinClustersTPC)
    if (nClustersTPC < fSelTrackMinClustersTPC) return kFALSE;
  if (fSelTrackMinClustersITS)
    if (nClustersITS < fSelTrackMinClustersITS) return kFALSE;
  if (fSelTrackMaxChi2PerClusterTPC)
    if (chi2PerClusterTPC > fSelTrackMaxChi2PerClusterTPC) return kFALSE; 
  if (fSelTrackMaxChi2PerClusterITS)
    if (chi2PerClusterITS > fSelTrackMaxChi2PerClusterITS) return kFALSE; 

  if (fSelTrackMaxCov11)
    if (extCov[0]  > fSelTrackMaxCov11) return kFALSE;
  if (fSelTrackMaxCov22)
    if (extCov[2]  > fSelTrackMaxCov22) return kFALSE;
  if (fSelTrackMaxCov33)
    if (extCov[5]  > fSelTrackMaxCov33) return kFALSE;
  if (fSelTrackMaxCov44)
    if (extCov[9]  > fSelTrackMaxCov44) return kFALSE;
  if(fSelTrackMaxCov55)
    if (extCov[14] > fSelTrackMaxCov55) return kFALSE;
  /*
  if(fMaxSigmaToVertexFlag)
    if(GetSigmaToVertex(track) > fMaxSigmaToVertex) return kFALSE;
  if(fITSRefitFlag)
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) return kFALSE;
  if(fTPCRefitFlag)
    if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) return kFALSE;
  if((Pt < fMinPt) || (Pt > fMaxPt)) return kFALSE;
  */
  return kTRUE;
}
