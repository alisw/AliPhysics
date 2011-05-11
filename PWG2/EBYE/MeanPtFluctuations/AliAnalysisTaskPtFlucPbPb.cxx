#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"
#include "TList.h"
#include "iostream"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"

#include "AliLog.h"

#include "AliAnalysisTaskPtFlucPbPb.h"


using namespace std;

// Analysis of Pt FLuctuations (PbPb)
// Author: Stefan Heckel
// Version of PbPb task: 5.0, 18.04.2011


ClassImp(AliAnalysisTaskPtFlucPbPb)

//________________________________________________________________________
AliAnalysisTaskPtFlucPbPb::AliAnalysisTaskPtFlucPbPb(const char *name)
  :AliAnalysisTaskSE(name),
  fESD(0),
  fOutputList(0),
  fPtSpec(0),
  fMult(0),
  fMultSum(0),
  fMultSumPt(0),
  fMultNrPairs(0),
  fCent(0),
  fCentSum(0),
  fCentSumPt(0),
  fCentNrPairs(0),
  fEta(0),
  fEtaPhiPlus(0),
  fEtaPhiMinus(0),
  fVtxZ(0),
  fVtxZCut(0),
  fVtxZCont(0),
  fVtxZTrackCuts(0),
  fEventMeanPt(0),
  fEventMeanPtSq(0),
  fMultEventMeanPt(0),
  fMultEventMeanPtSq(0),
  fCentEventMeanPt(0),
  fCentEventMeanPtSq(0),
  fTwoPartCorrEv(0),
  fTwoPartCorrEvSq(0),
  fTwoPartCorrEvCent(0),
  fTwoPartCorrEvCentSq(0),
  fESDTrackCuts(0),
  fMaxVertexZ(0),
  fNContributors(0),
  fUseCentrality(0),
  fMC(0)
{
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskPtFlucPbPb::~AliAnalysisTaskPtFlucPbPb()
{
  if(fOutputList) delete fOutputList;  fOutputList =0;
}

//________________________________________________________________________
void AliAnalysisTaskPtFlucPbPb::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(1, "RECREATE");
  fOutputList = new TList();
  fOutputList->SetOwner();


  fPtSpec = new TH1F("fPtSpec","Pt spectrum",50,0,2.5);

  fMult = new TH1F("fMult","Multiplicity distribution",30,0,3000);

  fMultSum = new TH1F("fMultSum","Sum of nrtracks of all events in multiplicity bins",30,0,3000);

  fMultSumPt = new TH1F("fMultSumPt","Sum of pTs of all events in multiplicity bins",30,0,3000);

  fMultNrPairs = new TH1F("fMultNrPairs","Sum of number of pairs in multiplicity bins",30,0,3000);

  fCent = new TH1F("fCent","Centrality distribution",21,0,105);

  fCentSum = new TH1F("fCentSum","Sum of nrtracks of all events in centrality bins",21,0,105);

  fCentSumPt = new TH1F("fCentSumPt","Sum of pTs of all events in centrality bins",21,0,105);

  fCentNrPairs = new TH1F("fCentNrPairs","Sum of number of pairs in centrality bins",21,0,105);


  fEta = new TH1F("fEta","Eta distribution",80,-2,2);

  fEtaPhiPlus = new TH1F("fEtaPhiPlus","Phi distribution for positive eta",62,0,6.2);

  fEtaPhiMinus = new TH1F("fEtaPhiMinus","Phi distribution for negative eta",62,0,6.2);

  fVtxZ = new TH1F("fVtxZ","Vertex Z distribution before cuts",100,-20,20);

  fVtxZCut = new TH1F("fVtxZCut","Vertex Z distribution after vtxZ cut",110,-11,11);

  fVtxZCont = new TH1F("fVtxZCont","Vertex Z distribution after nCont cut",110,-11,11);

  fVtxZTrackCuts = new TH1F("fVtxZTrackCuts","Vertex Z distribution after track cuts",110,-11,11);


  fEventMeanPt = new TH1F("fEventMeanPt","Mean-Pt event by event",50,0,2.5);

  fEventMeanPtSq = new TH1F("fEventMeanPtSq","Mean-Pt event by event squared",100,0,5);

  fMultEventMeanPt = new TH1F("fMultEventMeanPt","Mean-Pt event by event vs. multiplicity",30,0,3000);

  fMultEventMeanPtSq = new TH1F("fMultEventMeanPtSq","Mean-Pt event by event squared vs. multiplicity",30,0,3000);

  fCentEventMeanPt = new TH1F("fCentEventMeanPt","Mean-Pt event by event vs. centrality",21,0,105);

  fCentEventMeanPtSq = new TH1F("fCentEventMeanPtSq","Mean-Pt event by event squared vs. centrality",21,0,105);


  fTwoPartCorrEv = new TH1F("fTwoPartCorrEv","Two-particle correlator vs. multiplicity",30,0,3000);

  fTwoPartCorrEvSq = new TH1F("fTwoPartCorrEvSq","Two-particle correlator squared vs. multiplicity",30,0,3000);

  fTwoPartCorrEvCent = new TH1F("fTwoPartCorrEvCent","Two-particle correlator vs. centrality",21,0,105);

  fTwoPartCorrEvCentSq = new TH1F("fTwoPartCorrEvCentSq","Two-particle correlator squared vs. centrality",21,0,105);


  // Add histograms to the output list
  fOutputList->Add(fPtSpec);
  fOutputList->Add(fMult);
  fOutputList->Add(fMultSum);
  fOutputList->Add(fMultSumPt);
  fOutputList->Add(fMultNrPairs);
  fOutputList->Add(fCent);
  fOutputList->Add(fCentSum);
  fOutputList->Add(fCentSumPt);
  fOutputList->Add(fCentNrPairs);
  fOutputList->Add(fEta);
  fOutputList->Add(fEtaPhiPlus);
  fOutputList->Add(fEtaPhiMinus);
  fOutputList->Add(fVtxZ);
  fOutputList->Add(fVtxZCut);
  fOutputList->Add(fVtxZCont);
  fOutputList->Add(fVtxZTrackCuts);
  fOutputList->Add(fEventMeanPt);
  fOutputList->Add(fEventMeanPtSq);
  fOutputList->Add(fMultEventMeanPt);
  fOutputList->Add(fMultEventMeanPtSq);
  fOutputList->Add(fCentEventMeanPt);
  fOutputList->Add(fCentEventMeanPtSq);
  fOutputList->Add(fTwoPartCorrEv);
  fOutputList->Add(fTwoPartCorrEvSq);
  fOutputList->Add(fTwoPartCorrEvCent);
  fOutputList->Add(fTwoPartCorrEvCentSq);
}

//________________________________________________________________________
void AliAnalysisTaskPtFlucPbPb::UserExec(Option_t *)
{
  // Main loop
  // Called for each event


  // --- ESD event handler ---

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
  }
  else {
    fESD = esdH->GetEvent();
  }

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  if(!fESDTrackCuts) Printf("ERROR: No esd track cut");

  // --- End event handler ---


  // --- Vertex cuts ---

//   // TPC+ITS vertex
//   const AliESDVertex* vtxESD = fESD->GetPrimaryVertexTracks();
//   Double_t vtxZ = vtxESD->GetZv();
//   Double_t vtxNCont = vtxESD->GetNContributors();

  // TPConly vertex
  const AliESDVertex* vtxESDTPC = fESD->GetPrimaryVertexTPC();
  Double_t vtxZ = vtxESDTPC->GetZv();
  Double_t vtxNCont = vtxESDTPC->GetNContributors();

  fVtxZ->Fill(vtxZ); // VtxZ before cuts

  // Event cut on the z-Position of the vertex
  if(vtxZ > fMaxVertexZ || vtxZ < (-1.*fMaxVertexZ)) {
//     Printf("VertexZ out of range, Zv = %f",vtxZ);
    return;
  }

  fVtxZCut->Fill(vtxZ); // VtxZ after cut on vtxZ

  // Event cut on the number of contributors
  if(vtxNCont < fNContributors) {
//     Printf("Vertex has no contributors");
    return;
  }
//   else {
//     Printf("Vertex nContributors = %i",vtxNCont);
//   }

  fVtxZCont->Fill(vtxZ); // VtxZ after cut on nContributors

  // --- End vertex cuts ---


  Int_t nrTracks = 0;		// count for number of tracks which pass the track cuts
  Int_t nrESDTracks = 0;	// number of tracks in the ESD file
  Double_t sumPt = 0.;		// sum of track Pts
  Double_t twoPartCorrEv = 0.;	// Two-particle correlator of one event for multiplicity bin analysis
  Double_t twoPartCorrEvCent = 0.; // Two-particle correlator of one event for centrality bin analysis

  Double_t tracks[3000] = {0.};	// array of track Pts, needed for the two-particle correlator

  Double_t trackPt = 0., trackEta = 0., trackPhi = 0.;
  Double_t eventMeanPt = 0., eventMeanPtSq = 0., evMptMult = 0.;
  Double_t twoPartCorrPair = 0., twoPartCorrEvSq = 0.;
  Double_t twoPartCorrPairCent = 0., twoPartCorrEvCentSq = 0.;
  Double_t nrPairs = 0.;

  Double_t *nbins = 0x0; // Mean pT values for multiplicity bin analysis
  Double_t *centbins = 0x0; // Mean pT values for centrality bin analysis

  Double_t centralityVZERO = 0.;
  Int_t centralityVZEROBin = 0.;

  Int_t centBin = 0.;
  Double_t evMptCent = 0.;


//   Printf("MC: %d",fMC);


// --- Mean pT values ---

// !!!!! Have to be calculated in a first run - for each sample, which should be analysed, separately! !!!!!


// -- Mean pT values for multiplicity bins (first bin is for nrTracks = 0 and has to be 0) --

  // MC PbPb 2.76 ATeV 11a10a (Hijing) (Pt-Range: 0.15 - 2)
  Double_t nbinsMC276[32] = {0.000, 0.525, 0.535, 0.539, 0.541, 0.543, 0.544, 0.545, 0.545, 0.546, 0.547, 0.548, 0.548, 0.549, 0.549, 0.550, 0.551, 0.551, 0.552, 0.552, 0.553, 0.554, 0.554, 0.555, 0.556, 0.556, 0.557, 0.558, 0.549, 0.000, 0.000, 0.000};

  // Data PbPb 2.76 ATeV 10h.pass1 (Pt-Range: 0.15 - 2)
  Double_t nbinsData276[32] = {0.000, 0.580, 0.605, 0.618, 0.626, 0.632, 0.636, 0.640, 0.642, 0.645, 0.647, 0.648, 0.649, 0.650, 0.651, 0.652, 0.653, 0.653, 0.654, 0.654, 0.654, 0.654, 0.654, 0.654, 0.655, 0.656, 0.658, 0.660, 0.661, 0.659, 0.660, 0.656};

// -- End mean pT values for multiplicity bins --


// -- Mean pT values for centrality bins --

  // MC PbPb 2.76 ATeV 11a10a (Hijing) (Pt-Range: 0.15 - 2)
  Double_t centbinsMC276[21] = {0.555, 0.553, 0.551, 0.549, 0.548, 0.547, 0.545, 0.544, 0.543, 0.541, 0.539, 0.537, 0.536, 0.533, 0.530, 0.528, 0.524, 0.518, 0.505, 0.000, 0.534};

  // Data PbPb 2.76 ATeV 10h.pass1 (Pt-Range: 0.15 - 2)
  Double_t centbinsData276[21] = {0.654, 0.654, 0.652, 0.650, 0.647, 0.644, 0.640, 0.635, 0.629, 0.623, 0.616, 0.609, 0.601, 0.594, 0.586, 0.579, 0.568, 0.551, 0.535, 0.000, 0.621};

// -- End mean pT values for centrality bins --


// -- Selection of MC/Data; whole sample values --

if (fMC) { // - MC -

//   Printf(" -- MC, 2.76 ATeV -- ");

  nbins = nbinsMC276;
  centbins = centbinsMC276;

} // - End MC -
else { // - Data -

//   Printf(" -- Data, 2.76 ATeV -- ");

  nbins = nbinsData276;
  centbins = centbinsData276;

} // - End data -

// -- End selection of MC/Data; whole sample values --

// --- End mean pT values ---


  nrESDTracks = fESD->GetNumberOfTracks();
//   if( nrESDTracks ) { printf("Found event with %i tracks.\n",nrESDTracks); }


  if(fUseCentrality != 0) {

    AliCentrality *esdCentrality = fESD->GetCentrality();

    //V0
    if (fUseCentrality == 1){

      centralityVZERO = esdCentrality->GetCentralityPercentile("V0M");

      if      ( centralityVZERO >   0. && centralityVZERO <   5.) centralityVZEROBin =  0;
      else if ( centralityVZERO >=  5. && centralityVZERO <  10.) centralityVZEROBin =  5;
      else if ( centralityVZERO >= 10. && centralityVZERO <  15.) centralityVZEROBin = 10;
      else if ( centralityVZERO >= 15. && centralityVZERO <  20.) centralityVZEROBin = 15;
      else if ( centralityVZERO >= 20. && centralityVZERO <  25.) centralityVZEROBin = 20;
      else if ( centralityVZERO >= 25. && centralityVZERO <  30.) centralityVZEROBin = 25;
      else if ( centralityVZERO >= 30. && centralityVZERO <  35.) centralityVZEROBin = 30;
      else if ( centralityVZERO >= 35. && centralityVZERO <  40.) centralityVZEROBin = 35;
      else if ( centralityVZERO >= 40. && centralityVZERO <  45.) centralityVZEROBin = 40;
      else if ( centralityVZERO >= 45. && centralityVZERO <  50.) centralityVZEROBin = 45;
      else if ( centralityVZERO >= 50. && centralityVZERO <  55.) centralityVZEROBin = 50;
      else if ( centralityVZERO >= 55. && centralityVZERO <  60.) centralityVZEROBin = 55;
      else if ( centralityVZERO >= 60. && centralityVZERO <  65.) centralityVZEROBin = 60;
      else if ( centralityVZERO >= 65. && centralityVZERO <  70.) centralityVZEROBin = 65;
      else if ( centralityVZERO >= 70. && centralityVZERO <  75.) centralityVZEROBin = 70;
      else if ( centralityVZERO >= 75. && centralityVZERO <  80.) centralityVZEROBin = 75;
      else if ( centralityVZERO >= 80. && centralityVZERO <  85.) centralityVZEROBin = 80;
      else if ( centralityVZERO >= 85. && centralityVZERO <  90.) centralityVZEROBin = 85;
      else if ( centralityVZERO >= 90. && centralityVZERO <  95.) centralityVZEROBin = 90;
      else if ( centralityVZERO >= 95. && centralityVZERO <  99.) centralityVZEROBin = 95;
      else if ( centralityVZERO >= 99. ) centralityVZEROBin = 100;
      else if ( centralityVZERO <= 0.  ) centralityVZEROBin = 100;

      centBin = centralityVZEROBin / 5;
      evMptCent = centbins[centBin];
    }

//     Printf("\nRaw mult: %i   CentVZERO: %f   CentVZERObin: %i   centBin: %i   evMptCent: %.3f ",nrESDTracks,centralityVZERO,centralityVZEROBin,centBin,evMptCent);

  }


  // --- Loop over all tracks of one event ---

  for (Int_t iTracks = 0; iTracks < nrESDTracks; iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }

      if(!fESDTrackCuts->AcceptTrack(track))continue;

	trackPt = track->Pt();
	fPtSpec->Fill(trackPt);
	tracks[nrTracks] = trackPt;

	trackEta = track->Eta();
	fEta->Fill(trackEta);

	trackPhi = track->Phi();

	if (trackEta > 0.) {
	  fEtaPhiPlus->Fill(trackPhi);
	}
	else if (trackEta < 0.) {
	  fEtaPhiMinus->Fill(trackPhi);
	}

	sumPt = sumPt + trackPt;
	nrTracks++;

// 	printf("Track Pt = %.3f   Track Eta = %.3f   Track Phi = %3f\n",trackPt,trackEta,trackPhi);

  } // --- End track loop ---


  // --- Calculation of various values and filling of histograms (for remaining events with N_acc > 0) ---

  if(nrTracks != 0) {

    // VtxZ after all track cuts
    fVtxZTrackCuts->Fill(vtxZ);

    // Multiplicity distributions
    fMult->Fill(nrTracks);
    fMultSum->Fill(nrTracks,nrTracks);

    // Number of pairs in event
    nrPairs = 0.5 * nrTracks * (nrTracks-1);
    fMultNrPairs->Fill(nrTracks,nrPairs);

    // Calculation of mean Pt and mean Pt Squared
    eventMeanPt = sumPt / nrTracks;
    eventMeanPtSq = eventMeanPt * eventMeanPt;

    // Mean-Pt and Mean-Pt squared
    fEventMeanPt->Fill(eventMeanPt);
    fEventMeanPtSq->Fill(eventMeanPtSq);

    // Mean-Pt and Mean-Pt squared depending on multiplicity
    fMultEventMeanPt->Fill(nrTracks,eventMeanPt);
    fMultEventMeanPtSq->Fill(nrTracks,eventMeanPtSq);
    fMultSumPt->Fill(nrTracks,sumPt);

//     printf("nrTracks: %i   sumPt: %.8f   meanPt: %.8f   meanPtSq: %.8f\n",nrTracks,sumPt,eventMeanPt,eventMeanPtSq);


    // Centrality V0
    if (fUseCentrality == 1) {

      // Centrality distributions
      fCent->Fill(centralityVZEROBin);
      fCentSum->Fill(centralityVZEROBin,nrTracks);

      // Number of pairs in event
      fCentNrPairs->Fill(centralityVZEROBin,nrPairs);

      // Mean-Pt and Mean-Pt squared depending on centrality
      fCentEventMeanPt->Fill(centralityVZEROBin,eventMeanPt);
      fCentEventMeanPtSq->Fill(centralityVZEROBin,eventMeanPtSq);
      fCentSumPt->Fill(centralityVZEROBin,sumPt);

//       printf("CentV0bin: %i   NrTracks: %i   NrPairs: %.0f   SumPt: %.1f\n",centralityVZEROBin,nrTracks,nrPairs,sumPt);
    }


    // --- Two-particle correlator ---

    // -- Analysis in multiplicity bins --

    for (int k=1; k<32; k++) {
      if (nrTracks > 100 * (k-1) && nrTracks <= 100 * k) {
	evMptMult = nbins[k];
	break;
      }
    }
//     printf("nrTracks = %3i   evMptMult = %.3f \n",nrTracks,evMptMult);

    for (int i=0; i<nrTracks; i++) {
      for (int j=i+1; j<nrTracks; j++) {
	twoPartCorrPair = (tracks[i] - evMptMult) * (tracks[j] - evMptMult);
	twoPartCorrEv = twoPartCorrEv + twoPartCorrPair;
      }
    }

    twoPartCorrEvSq = twoPartCorrEv * twoPartCorrEv;
    fTwoPartCorrEv->Fill(nrTracks,twoPartCorrEv);
    fTwoPartCorrEvSq->Fill(nrTracks,twoPartCorrEvSq);

//     printf("twoPartCorrEv = %.3f   twoPartCorrEvSq = %.3f \n\n",twoPartCorrEv,twoPartCorrEvSq);

    // -- End analysis in multiplicity bins --


    // -- Analysis in centrality bins --

    // Centrality V0
    if (fUseCentrality == 1) {

//       printf("CentV0bin: %i  ",centralityVZEROBin);

      if (centralityVZEROBin < 91) {

// 	printf("CentV0bin: %i   ",centralityVZEROBin);

	for (int i=0; i<nrTracks; i++) {
	  for (int j=i+1; j<nrTracks; j++) {
	    twoPartCorrPairCent = (tracks[i] - evMptCent) * (tracks[j] - evMptCent);
	    twoPartCorrEvCent = twoPartCorrEvCent + twoPartCorrPairCent;
	  }
	}

	twoPartCorrEvCentSq = twoPartCorrEvCent * twoPartCorrEvCent;
	fTwoPartCorrEvCent->Fill(centralityVZEROBin,twoPartCorrEvCent);
	fTwoPartCorrEvCentSq->Fill(centralityVZEROBin,twoPartCorrEvCentSq);

// 	printf("CentV0bin: %i   evMptCent: %.3f   CorrEvCent: %.2f   CorrEvCentSq: %.2f\n",centralityVZEROBin,evMptCent,twoPartCorrEvCent,twoPartCorrEvCentSq);
      }
    }

    // -- End analysis in centrality bins --

    // --- End two-particle correlator ---


  } // --- End calculation of various values and filling of histograms ---


  // Post output data
  PostData(1, fOutputList);
}

void AliAnalysisTaskPtFlucPbPb::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    Printf("ERROR: fOutputList not available\n");
    return;
  }

}
