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

#include "AliLog.h"

#include "AliAnalysisTaskPtFluc.h"


using namespace std;

// Analysis of Pt FLuctuations (pp)
// Author: Stefan Heckel
// Version of pp task:   8.0, 19.04.2011


ClassImp(AliAnalysisTaskPtFluc)

//________________________________________________________________________
AliAnalysisTaskPtFluc::AliAnalysisTaskPtFluc(const char *name)
  :AliAnalysisTaskSE(name),
  fESD(0),
  fOutputList(0),
  fPtSpec(0),
  fMult(0),
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
  fTwoPartCorrEv(0),
  fTwoPartCorrEvSq(0),
  fTwoPartCorrEvSample(0),
  fTwoPartCorrEvSampleSq(0),
  fESDTrackCuts(0),
  fMaxVertexZ(0),
  fNContributors(0),
  fMC(0)
{
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskPtFluc::~AliAnalysisTaskPtFluc()
{
  if(fOutputList) delete fOutputList;  fOutputList =0;
}

//________________________________________________________________________
void AliAnalysisTaskPtFluc::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(1, "RECREATE");
  fOutputList = new TList();
  fOutputList->SetOwner();


  fPtSpec = new TH1F("fPtSpec","Pt spectrum",50,0,2.5);

  fMult = new TH1F("fMult","Multiplicity distribution",120,-0.5,119.5);

  fEta = new TH1F("fEta","Eta distribution",80,-2,2);

  fEtaPhiPlus = new TH1F("fEtaPhiPlus","Phi distribution for positive eta",62,0,6.2);

  fEtaPhiMinus = new TH1F("fEtaPhiMinus","Phi distribution for negative eta",62,0,6.2);

  fVtxZ = new TH1F("fVtxZ","Vertex Z distribution before cuts",100,-20,20);

  fVtxZCut = new TH1F("fVtxZCut","Vertex Z distribution after vtxZ cut",110,-11,11);

  fVtxZCont = new TH1F("fVtxZCont","Vertex Z distribution after nCont cut",110,-11,11);

  fVtxZTrackCuts = new TH1F("fVtxZTrackCuts","Vertex Z distribution after track cuts",110,-11,11);


  fEventMeanPt = new TH1F("fEventMeanPt","Mean-Pt event by event",50,0,2.5);

  fEventMeanPtSq = new TH1F("fEventMeanPtSq","Mean-Pt event by event squared",100,0,5);

  fMultEventMeanPt = new TH1F("fMultEventMeanPt","Mean-Pt event by event vs. Multiplicity",100,-0.5,99.5);

  fMultEventMeanPtSq = new TH1F("fMultEventMeanPtSq","Mean-Pt event by event squared vs. Multiplicity",100,-0.5,99.5);


  fTwoPartCorrEv = new TH1F("fTwoPartCorrEv","Two-particle correlator vs. Multiplicity",100,-0.5,99.5);

  fTwoPartCorrEvSq = new TH1F("fTwoPartCorrEvSq","Two-particle correlator squared vs. Multiplicity",100,-0.5,99.5);

  fTwoPartCorrEvSample = new TH1F("fTwoPartCorrEvSample","Two-particle correlator vs. Multiplicity",100,-0.5,99.5);

  fTwoPartCorrEvSampleSq = new TH1F("fTwoPartCorrEvSampleSq","Two-particle correlator squared vs. Multiplicity",100,-0.5,99.5);


  // Add histograms to the output list
  fOutputList->Add(fPtSpec);
  fOutputList->Add(fMult);
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
  fOutputList->Add(fTwoPartCorrEv);
  fOutputList->Add(fTwoPartCorrEvSq);
  fOutputList->Add(fTwoPartCorrEvSample);
  fOutputList->Add(fTwoPartCorrEvSampleSq);
}

//________________________________________________________________________
void AliAnalysisTaskPtFluc::UserExec(Option_t *)
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

  // TPC+ITS vertex
  const AliESDVertex* vtxESD = fESD->GetPrimaryVertexTracks();
  Double_t vtxZ = vtxESD->GetZv();
  Double_t vtxNCont = vtxESD->GetNContributors();

//   // TPConly vertex
//   const AliESDVertex* vtxESDTPC = fESD->GetPrimaryVertexTPC();
//   Double_t vtxZ = vtxESDTPC->GetZv();
//   Double_t vtxNCont = vtxESDTPC->GetNContributors();

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
  Double_t twoPartCorrEvSample = 0.; // Two-part. corr. of one event for whole sample analysis

  Double_t tracks[120] = {0.};	// array of track Pts, needed for the two-particle correlator

  Double_t trackPt, trackEta, trackPhi;
  Double_t eventMeanPt, eventMeanPtSq, evMptMult;
  Double_t twoPartCorrPair, twoPartCorrEvSq;
  Double_t twoPartCorrPairSample, twoPartCorrEvSampleSq;

  Double_t *nbins = 0x0; // Mean pT values for multiplicity bin analysis
  Double_t evMptSample;  // Mean pT value for whole sample analysis

  Double_t energy; // beam energy

  energy = fESD->GetBeamEnergy();


//   Printf("MC: %d,  %f",fMC,energy);


// --- Mean pT values ---

// !!!!! Have to be calculated in a first run - for each sample, which should be analysed, separately! !!!!!


// -- Mean pT values for the multiplicity bins (starting with N = 0, therefore the first value has to be 0.000) --

// MC 900 GeV 10e13 (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2)
Double_t nbinsMC900[100] = {0.000, 0.456, 0.475, 0.490, 0.507, 0.523, 0.535, 0.546, 0.554, 0.561, 0.568, 0.572, 0.578, 0.582, 0.586, 0.590, 0.593, 0.596, 0.602, 0.601, 0.606, 0.605, 0.610, 0.618, 0.610, 0.620, 0.633, 0.615, 0.609, 0.620, 0.632, 0.634, 0.661, 0.647, 0.560, 0.671, 0.660, 0.678, 0.000, 0.612, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

// MC 7 TeV 10f6a (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2)
Double_t nbinsMC7[100] = {0.000, 0.453, 0.475, 0.498, 0.523, 0.546, 0.563, 0.577, 0.587, 0.595, 0.603, 0.610, 0.616, 0.622, 0.627, 0.632, 0.636, 0.640, 0.644, 0.648, 0.651, 0.655, 0.658, 0.662, 0.664, 0.667, 0.669, 0.671, 0.674, 0.677, 0.679, 0.682, 0.683, 0.686, 0.688, 0.690, 0.692, 0.695, 0.696, 0.696, 0.699, 0.700, 0.701, 0.705, 0.709, 0.705, 0.701, 0.705, 0.714, 0.718, 0.727, 0.702, 0.725, 0.738, 0.724, 0.728, 0.710, 0.729, 0.711, 0.748, 0.702, 0.708, 0.732, 0.000, 0.683, 0.000, 0.000, 0.675, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};


// // MC 900 GeV 10e12 (Phojet) (Pt-Range: 0.15 - 2)
// Double_t nbinsMC900[100] = {0.000, 0.478, 0.483, 0.488, 0.494, 0.498, 0.502, 0.505, 0.507, 0.509, 0.513, 0.516, 0.520, 0.522, 0.526, 0.528, 0.531, 0.535, 0.539, 0.540, 0.544, 0.540, 0.548, 0.548, 0.549, 0.556, 0.562, 0.549, 0.557, 0.550, 0.561, 0.559, 0.554, 0.547, 0.539, 0.564, 0.546, 0.603, 0.532, 0.583, 0.656, 0.531, 0.458, 0.547, 0.000, 0.633, 0.000, 0.603, 0.588, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

// // MC 7 TeV 10f6 (Phojet) (Pt-Range: 0.15 - 2)
// Double_t nbinsMC7[100] = {0.000, 0.488, 0.499, 0.512, 0.526, 0.536, 0.544, 0.549, 0.553, 0.557, 0.560, 0.563, 0.566, 0.569, 0.572, 0.574, 0.577, 0.579, 0.580, 0.582, 0.584, 0.585, 0.586, 0.587, 0.588, 0.589, 0.589, 0.590, 0.591, 0.591, 0.590, 0.591, 0.592, 0.591, 0.589, 0.591, 0.590, 0.588, 0.590, 0.589, 0.590, 0.588, 0.589, 0.589, 0.591, 0.588, 0.582, 0.580, 0.582, 0.586, 0.575, 0.575, 0.560, 0.570, 0.569, 0.604, 0.616, 0.549, 0.614, 0.611, 0.480, 0.628, 0.568, 0.000, 0.677, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};


// Data 900 GeV 10c.pass3 (Pt-Range: 0.15 - 2)
Double_t nbinsData900[100] = {0.000, 0.490, 0.491, 0.495, 0.503, 0.512, 0.521, 0.529, 0.536, 0.542, 0.548, 0.551, 0.556, 0.561, 0.564, 0.567, 0.570, 0.573, 0.575, 0.578, 0.580, 0.582, 0.583, 0.590, 0.591, 0.591, 0.590, 0.598, 0.593, 0.608, 0.609, 0.609, 0.587, 0.615, 0.603, 0.613, 0.560, 0.583, 0.614, 0.649, 0.584, 0.584, 0.594, 0.000, 0.622, 0.000, 0.000, 0.703, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

// Data 2.76 TeV 11a.pass1 (Pt-Range: 0.15 - 2)
Double_t nbinsData276[100] = {0.000, 0.491, 0.500, 0.511, 0.524, 0.536, 0.548, 0.558, 0.567, 0.574, 0.581, 0.586, 0.592, 0.596, 0.601, 0.605, 0.609, 0.613, 0.616, 0.619, 0.622, 0.625, 0.628, 0.629, 0.632, 0.634, 0.637, 0.639, 0.641, 0.645, 0.645, 0.646, 0.649, 0.652, 0.655, 0.658, 0.656, 0.655, 0.659, 0.668, 0.660, 0.672, 0.672, 0.665, 0.669, 0.662, 0.666, 0.669, 0.677, 0.673, 0.595, 0.642, 0.700, 0.623, 0.670, 0.704, 0.694, 0.000, 0.000, 0.667, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

// Data 7 TeV 10d.pass2 + 10e.pass2 (Pt-Range: 0.15 - 2)
Double_t nbinsData7[100] = {0.000, 0.495, 0.498, 0.505, 0.516, 0.528, 0.539, 0.550, 0.559, 0.566, 0.573, 0.579, 0.584, 0.589, 0.593, 0.597, 0.601, 0.605, 0.608, 0.611, 0.614, 0.617, 0.619, 0.622, 0.624, 0.626, 0.629, 0.631, 0.633, 0.634, 0.636, 0.638, 0.640, 0.642, 0.643, 0.645, 0.646, 0.648, 0.649, 0.651, 0.652, 0.654, 0.654, 0.656, 0.657, 0.658, 0.659, 0.660, 0.662, 0.662, 0.663, 0.665, 0.664, 0.665, 0.668, 0.667, 0.670, 0.666, 0.670, 0.670, 0.677, 0.673, 0.672, 0.679, 0.680, 0.667, 0.682, 0.687, 0.675, 0.658, 0.677, 0.663, 0.686, 0.671, 0.682, 0.683, 0.667, 0.000, 0.000, 0.662, 0.720, 0.683, 0.720, 0.622, 0.000, 0.729, 0.000, 0.603, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

// -- End mean pT values for the multiplicity bins --


// -- Selection of MC/Data and energy; whole sample values --

if (fMC) { // - MC -

  if (energy < 1000.) {

//     Printf(" -- MC, 900 GeV -- ");

    evMptSample = 0.5307; // MC 900 GeV 10e13 (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2)
//     evMptSample = 0.5044; // MC 900 GeV 10e12 (Phojet) (Pt-Range: 0.15 - 2)

    nbins = nbinsMC900;
  }
  else {

//     Printf(" -- MC, 7 TeV -- ");

    evMptSample = 0.5847; // MC 7 TeV 10f6a (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2)
//     evMptSample = 0.5517; // MC 7 TeV 10f6 (Phojet) (Pt-Range: 0.15 - 2)

    nbins = nbinsMC7;
  }

} // - End MC -
else { // - Data -

  if (energy < 1000.) {

//     Printf(" -- Data, 900 GeV -- ");

    evMptSample = 0.5267; // Data 900 GeV 10c.pass3 (Pt-Range: 0.15 - 2)
    nbins = nbinsData900;
  }
  else if (energy > 1000. && energy < 3000.) {

//     Printf(" -- Data, 2.76 TeV -- ");

    evMptSample = 0.5622; // Data 2.76 TeV 11a.pass1 (Pt-Range: 0.15 - 2)
    nbins = nbinsData276;
  }
  else {

//     Printf(" -- Data, 7 TeV -- ");

    evMptSample = 0.5743; // Data 7 TeV 10d.pass2 + 10e.pass2 (Pt-Range: 0.15 - 2)
    nbins = nbinsData7;
  }

} // - End data -

// -- End selection of MC/Data and energy; whole sample values --

// --- End mean pT values ---


  nrESDTracks = fESD->GetNumberOfTracks();
//   if( nrESDTracks ) { printf("Found event with %i tracks.\n",nrESDTracks); }


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

    // Multiplicity distribution
    fMult->Fill(nrTracks);

    // Calculation of mean Pt and mean Pt Squared
    eventMeanPt = sumPt / nrTracks;
    eventMeanPtSq = eventMeanPt * eventMeanPt;

    // Mean-Pt and Mean-Pt squared
    fEventMeanPt->Fill(eventMeanPt);
    fEventMeanPtSq->Fill(eventMeanPtSq);

    // Mean-Pt and Mean-Pt squared depending on multiplicity
    fMultEventMeanPt->Fill(nrTracks,eventMeanPt);
    fMultEventMeanPtSq->Fill(nrTracks,eventMeanPtSq);

//     printf("nrTracks: %i   sumPt: %.8f   meanPt: %.8f   meanPtSq: %.8f\n",nrTracks,sumPt,eventMeanPt,eventMeanPtSq);


    // --- Two-particle correlator ---

    evMptMult = nbins[nrTracks];
//     printf("nrTracks = %3i   evMptMult = %.3f   evMptSample = %.3f \n\n",nrTracks,evMptMult,evMptSample);

    for (int i=0; i<nrTracks; i++) {

//       printf("--- TrackPt = %.3f ",tracks[i]);

      for (int j=i+1; j<nrTracks; j++) {

	twoPartCorrPair = (tracks[i] - evMptMult) * (tracks[j] - evMptMult); // Multiplicity bins
	twoPartCorrEv = twoPartCorrEv + twoPartCorrPair;

	twoPartCorrPairSample = (tracks[i] - evMptSample) * (tracks[j] - evMptSample); // Whole sample
	twoPartCorrEvSample = twoPartCorrEvSample + twoPartCorrPairSample;

// 	printf("twoPartCorrPair = %.3f   twoPartCorrPairSample = %.3f \n",twoPartCorrPair,twoPartCorrPairSample);
      }
//       printf("twoPartCorrEv = %.3f   twoPartCorrEvSample = %.3f \n",twoPartCorrEv,twoPartCorrEvSamples);
    }
//     printf("Total Event: twoPartCorrEv = %.3f   twoPartCorrEvSample = %.3f \n\n",twoPartCorrEv,twoPartCorrEvSample);

    // Multiplicity bins
    twoPartCorrEvSq = twoPartCorrEv * twoPartCorrEv;
    fTwoPartCorrEv->Fill(nrTracks,twoPartCorrEv);
    fTwoPartCorrEvSq->Fill(nrTracks,twoPartCorrEvSq);

    // Whole sample
    twoPartCorrEvSampleSq = twoPartCorrEvSample * twoPartCorrEvSample;
    fTwoPartCorrEvSample->Fill(nrTracks,twoPartCorrEvSample);
    fTwoPartCorrEvSampleSq->Fill(nrTracks,twoPartCorrEvSampleSq);

    // --- End two-particle correlator ---


  } // --- End calculation of various values and filling of histograms ---


  // Post output data
  PostData(1, fOutputList);
}

void AliAnalysisTaskPtFluc::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    Printf("ERROR: fOutputList not available\n");
    return;
  }

}
