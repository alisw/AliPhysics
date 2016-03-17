#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TRandom3.h"
#include "iostream"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskPtFluc.h"


using namespace std;

// Analysis of Pt Fluctuations (pp)
// Author: Stefan Heckel
// Version of pp task:   11.4, 11.06.2012


ClassImp(AliAnalysisTaskPtFluc)

//________________________________________________________________________
AliAnalysisTaskPtFluc::AliAnalysisTaskPtFluc(const char *name)
  :AliAnalysisTaskSE(name),
  fESD(0),
  fMCev(0),
  fRandom3(0),
  fOutputList(0),
  fPtSpec(0),
  fPtSpec2(0),
  fMult(0),
  fEta(0),
  fEtaPhiPlus(0),
  fEtaPhiMinus(0),
  fVtxZ(0),
  fVtxZCut(0),
  fVtxZCont(0),
  fVtxZCutDiff(0),
  fVtxZPileup(0),
  fVtxZTrackCuts(0),
  fVtxZDiff1(0),
  fVtxZDiff2(0),
  fVtxZDiff3(0),
  fVtxZDiff1b(0),
  fVtxZDiff2b(0),
  fVtxZDiff3b(0),
  fEventMeanPt(0),
  fEventMeanPtSq(0),
  fEventMeanPtMult(0),
  fMultEventMeanPt(0),
  fMultEventMeanPtSq(0),
  fTwoPartCorrEv(0),
  fTwoPartCorrEvSq(0),
  fTwoPartCorrEvSample(0),
  fTwoPartCorrEvSampleSq(0),
  fESDTrackCuts(0),
  fMaxVertexZ(0),
  fMaxVertexZDiff1(0),
  fNContributors(0),
  fMC(0),
  fMCType(0)
{
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskPtFluc::~AliAnalysisTaskPtFluc()
{
  if (fRandom3) delete fRandom3; fRandom3 = 0;

  if (fOutputList) {
    fOutputList->Clear();
    delete fOutputList;
  }
  fOutputList = 0;

  if (fESDTrackCuts) delete fESDTrackCuts; fESDTrackCuts = 0;

}

//________________________________________________________________________
void AliAnalysisTaskPtFluc::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(1, "RECREATE");
  fOutputList = new TList();
  fOutputList->SetOwner();


  fPtSpec = new TH1F("fPtSpec","Pt spectrum",100,0,2.5);
  fPtSpec2 = new TH1F("fPtSpec2","Pt spectrum 2 - MC ESD",100,0,2.5);
  fMult = new TH1F("fMult","Multiplicity distribution",120,-0.5,119.5);

  fEta = new TH1F("fEta","Eta distribution",80,-2,2);
  fEtaPhiPlus = new TH1F("fEtaPhiPlus","Phi distribution for positive eta",62,0,6.2);
  fEtaPhiMinus = new TH1F("fEtaPhiMinus","Phi distribution for negative eta",62,0,6.2);

  fVtxZ = new TH1F("fVtxZ","Vertex Z distribution before cuts",100,-20,20);
  fVtxZCut = new TH1F("fVtxZCut","Vertex Z distribution after vtxZ cut",110,-11,11);
  fVtxZCont = new TH1F("fVtxZCont","Vertex Z distribution after nCont cut",110,-11,11);
  fVtxZCutDiff = new TH1F("fVtxZCutDiff","Vertex Z distribution after cut on vtx Z Diff",110,-11,11);
  fVtxZPileup = new TH1F("fVtxZPileup","Vertex Z distribution after pileup rejection",110,-11,11);
  fVtxZTrackCuts = new TH1F("fVtxZTrackCuts","Vertex Z distribution after track cuts",110,-11,11);

  fVtxZDiff1 = new TH1F("fVtxZDiff1","Difference 1 between vertex Z distributions",100,-5,5);
  fVtxZDiff2 = new TH1F("fVtxZDiff2","Difference 2 between vertex Z distributions",100,-5,5);
  fVtxZDiff3 = new TH1F("fVtxZDiff3","Difference 3 between vertex Z distributions",100,-5,5);
  fVtxZDiff1b = new TH1F("fVtxZDiff1b","Difference 1 between vertex Z distributions after all cuts",100,-5,5);
  fVtxZDiff2b = new TH1F("fVtxZDiff2b","Difference 2 between vertex Z distributions after all cuts",100,-5,5);
  fVtxZDiff3b = new TH1F("fVtxZDiff3b","Difference 3 between vertex Z distributions after all cuts",100,-5,5);


  fEventMeanPt = new TH1F("fEventMeanPt","Mean-Pt distribution",50,0,2.5);
  fEventMeanPtSq = new TH1F("fEventMeanPtSq","Mean-Pt squared distribution",100,0,5);
  fEventMeanPtMult = new TH2F("fEventMeanPtMult","Mean-Pt for single events vs. multiplicity",20,0,100,200,0.,2.);

  fMultEventMeanPt = new TH1F("fMultEventMeanPt","Mean-Pt event by event vs. Multiplicity",100,-0.5,99.5);
  fMultEventMeanPtSq = new TH1F("fMultEventMeanPtSq","Mean-Pt event by event squared vs. Multiplicity",100,-0.5,99.5);


  fTwoPartCorrEv = new TH1F("fTwoPartCorrEv","Two-particle correlator vs. Multiplicity",100,-0.5,99.5);
  fTwoPartCorrEvSq = new TH1F("fTwoPartCorrEvSq","Two-particle correlator squared vs. Multiplicity",100,-0.5,99.5);
  fTwoPartCorrEvSample = new TH1F("fTwoPartCorrEvSample","Two-particle correlator vs. Multiplicity",100,-0.5,99.5);
  fTwoPartCorrEvSampleSq = new TH1F("fTwoPartCorrEvSampleSq","Two-particle correlator squared vs. Multiplicity",100,-0.5,99.5);


  // Add histograms to the output list
  fOutputList->Add(fPtSpec);
  fOutputList->Add(fPtSpec2);
  fOutputList->Add(fMult);
  fOutputList->Add(fEta);
  fOutputList->Add(fEtaPhiPlus);
  fOutputList->Add(fEtaPhiMinus);
  fOutputList->Add(fVtxZ);
  fOutputList->Add(fVtxZCut);
  fOutputList->Add(fVtxZCont);
  fOutputList->Add(fVtxZCutDiff);
  fOutputList->Add(fVtxZPileup);
  fOutputList->Add(fVtxZTrackCuts);
  fOutputList->Add(fVtxZDiff1);
  fOutputList->Add(fVtxZDiff2);
  fOutputList->Add(fVtxZDiff3);
  fOutputList->Add(fVtxZDiff1b);
  fOutputList->Add(fVtxZDiff2b);
  fOutputList->Add(fVtxZDiff3b);
  fOutputList->Add(fEventMeanPt);
  fOutputList->Add(fEventMeanPtSq);
  fOutputList->Add(fEventMeanPtMult);
  fOutputList->Add(fMultEventMeanPt);
  fOutputList->Add(fMultEventMeanPtSq);
  fOutputList->Add(fTwoPartCorrEv);
  fOutputList->Add(fTwoPartCorrEvSq);
  fOutputList->Add(fTwoPartCorrEvSample);
  fOutputList->Add(fTwoPartCorrEvSampleSq);

  // Post output data (if histograms are not used later, PostData is at least called here)
  PostData(1, fOutputList);

}

//________________________________________________________________________
void AliAnalysisTaskPtFluc::UserExec(Option_t *)
{
  // Main loop
  // Called for each event


  // --- ESD event handler ---
  if (!fMC || fMCType < 3) {

//     Printf("ESD handler");

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    }
    else {
      fESD = (AliESDEvent*)esdH->GetEvent();
    }

    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }

    if(!fESDTrackCuts) Printf("ERROR: No esd track cut");

  } // --- End ESD event handler ---


  // --- MC event handler ---

  AliStack *stack = NULL;

  if (fMC) {

//     Printf("MC handler");

    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

    if (!mcH) {
      Printf("ERROR: Could not get MCInputHandler");
      return;
    }
    else {
     fMCev = mcH->MCEvent();
    }

    if (!fMCev) {
      Printf("ERROR: fMCev not available");
      return;
    }

   stack = fMCev->Stack();

  }

  // --- End MC event handler ---


  // --- Vertex cuts ---

  Double_t vtxZGlobal=0., vtxZSPD=0., vtxZTPC=0.;
  Double_t vtxNContGlobal=0.;
//   Double_t  vtxNContSPD=0., vtxNContTPC=0.;
  Double_t vtxZ=0., vtxNCont=0.;
  Double_t vtxZdiff1=0., vtxZdiff2=0., vtxZdiff3=0.;

 if (!fMC || fMCType < 3) {

  // Global vertex
  const AliESDVertex* vtxESD = fESD->GetPrimaryVertexTracks();
  vtxZGlobal = vtxESD->GetZ();
  vtxNContGlobal = vtxESD->GetNContributors();

  // SPD vertex
  const AliESDVertex* vtxESDSPD = fESD->GetPrimaryVertexSPD();
  vtxZSPD = vtxESDSPD->GetZ();
//   vtxNContSPD = vtxESDSPD->GetNContributors();

  // TPC vertex
  const AliESDVertex* vtxESDTPC = fESD->GetPrimaryVertexTPC();
  vtxZTPC = vtxESDTPC->GetZ();
//   vtxNContTPC = vtxESDTPC->GetNContributors();

  vtxZ = vtxZGlobal;
  vtxNCont = vtxNContGlobal;

  fVtxZ->Fill(vtxZ); // VtxZ before cuts

  // Differences between vertices
  vtxZdiff1 = vtxZTPC - vtxZGlobal;
  vtxZdiff2 = vtxZTPC - vtxZSPD;
  vtxZdiff3 = vtxZGlobal - vtxZSPD;
  fVtxZDiff1->Fill(vtxZdiff1);
  fVtxZDiff2->Fill(vtxZdiff2);
  fVtxZDiff3->Fill(vtxZdiff3);

  // Event cut on the z-position of the vertex
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

  fVtxZCont->Fill(vtxZ); // VtxZ after cut on nContributors

  // Event cut on the difference of the z-positions of the vertices (TPC - global)
  if (fMaxVertexZDiff1 > 0.) {
    if(vtxZdiff1 > fMaxVertexZDiff1 || vtxZdiff1 < (-1.*fMaxVertexZDiff1)) {
//       Printf("VertexZ Diff (TPC - global) out of range, vtxZdiff1 = %f",vtxZ);
      return;
    }
  }

  fVtxZCutDiff->Fill(vtxZ); // VtxZ after cut on vtxZDiff

  } // --- End vertex cuts ---


  // --- Pileup rejection ---
  if (!fMC || fMCType < 3) {

    if (fESD->IsPileupFromSPD(3,0.8,3.,2.,5.)) {
//       Printf("IsPileupFromSPD");
      return;
    }

    fVtxZPileup->Fill(vtxZ); // VtxZ after pileup rejection

  } // --- End pileup rejection ---


  Int_t nrTracks = 0;		// count for number of tracks which pass the track cuts
//   Int_t nrTracks2 = 0;		// count for number of tracks which pass the track cuts - for MC ESD vs. MC truth
  Int_t nrESDTracks = 0;	// number of tracks in the ESD file
  Int_t nrMCTracks = 0;		// number of tracks in the MC file
  Double_t sumPt = 0.;		// sum of track Pts
  Double_t twoPartCorrEv = 0.;	// Two-particle correlator of one event for multiplicity bin analysis
  Double_t twoPartCorrEvSample = 0.; // Two-part. corr. of one event for whole sample analysis

  Double_t tracks[120] = {0.};	// array of track Pts, needed for the two-particle correlator

  Double_t trackPt=0., trackEta=0., trackPhi=0.;
  Double_t trackPt2=0.; // for MC ESD
  Double_t eventMeanPt=0., eventMeanPtSq=0., evMptMult=0.;
  Double_t twoPartCorrPair=0., twoPartCorrEvSq=0.;
  Double_t twoPartCorrPairSample=0., twoPartCorrEvSampleSq=0.;
//   Double_t nrPairs=0.;

  Double_t *nbins = 0x0; // Mean pT values for multiplicity bin analysis
  Double_t evMptSample=0.;  // Mean pT value for whole sample analysis

  Double_t random=0., trackPt40=0., ptRatio=0.;
  Int_t trackPtBin=0;

  fRandom3 = new TRandom3();

  if (fMCType == 2 || fMCType == 4) { // MC, Type = mod. MC truth
   fRandom3->SetSeed(0);
  }


  Double_t energy = 0.; // beam energy ---> has to be set by the user (in GeV) for MCType >= 3 (w/o ESD) <---

  if (!fMC || fMCType < 3) {
    energy = fESD->GetBeamEnergy();
  }


// --- Mean pT values ---

// !!!!! Have to be calculated in a first run - for each sample, which should be analysed, separately! !!!!!


// -- Mean pT values for the multiplicity bins (starting with N = 0, therefore the first value has to be 0.000) --

// MC 900 GeV 10e13 (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2) ---> TPC standalone <------> values for MC ESD <---
Double_t nbinsMC900[100] = {0.000, 0.452, 0.468, 0.482, 0.497, 0.511, 0.524, 0.535, 0.544, 0.550, 0.557, 0.563, 0.567, 0.571, 0.576, 0.580, 0.583, 0.586, 0.590, 0.592, 0.595, 0.599, 0.600, 0.602, 0.604, 0.611, 0.604, 0.618, 0.614, 0.612, 0.632, 0.623, 0.617, 0.612, 0.640, 0.632, 0.649, 0.609, 0.658, 0.615, 0.677, 0.614, 0.646, 0.000, 0.000, 0.000, 0.000, 0.000, 0.672, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

// MC 2.76 TeV 11b10a (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2) ---> TPC standalone <------> values for MC ESD <---
Double_t nbinsMC276[100] = {0.000, 0.451, 0.467, 0.485, 0.505, 0.523, 0.540, 0.551, 0.563, 0.569, 0.577, 0.583, 0.588, 0.594, 0.598, 0.603, 0.608, 0.611, 0.615, 0.617, 0.622, 0.622, 0.626, 0.630, 0.633, 0.632, 0.637, 0.633, 0.641, 0.643, 0.650, 0.649, 0.655, 0.650, 0.669, 0.664, 0.653, 0.671, 0.668, 0.683, 0.693, 0.662, 0.765, 0.643, 0.677, 0.660, 0.000, 0.689, 0.708, 0.653, 0.635, 0.000, 0.719, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

// MC 7 TeV 10f6a (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2) ---> TPC standalone <------> values for MC ESD <---
Double_t nbinsMC7[100] = {0.000, 0.447, 0.465, 0.485, 0.508, 0.530, 0.548, 0.562, 0.573, 0.581, 0.589, 0.595, 0.601, 0.607, 0.612, 0.617, 0.621, 0.625, 0.629, 0.633, 0.637, 0.640, 0.643, 0.646, 0.649, 0.652, 0.654, 0.657, 0.659, 0.662, 0.664, 0.666, 0.668, 0.670, 0.672, 0.674, 0.677, 0.677, 0.680, 0.682, 0.684, 0.685, 0.687, 0.688, 0.693, 0.692, 0.692, 0.693, 0.696, 0.698, 0.696, 0.700, 0.702, 0.697, 0.701, 0.708, 0.712, 0.714, 0.705, 0.717, 0.702, 0.701, 0.716, 0.716, 0.761, 0.706, 0.750, 0.715, 0.731, 0.728, 0.711, 0.733, 0.000, 0.714, 0.704, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.652, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};
// MC 7 TeV 10f6a (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2) ---> values for MC truth <---
// Double_t nbinsMC7[100] = {0.000, 0.441, 0.455, 0.473, 0.493, 0.515, 0.534, 0.550, 0.561, 0.570, 0.578, 0.585, 0.591, 0.597, 0.602, 0.607, 0.612, 0.616, 0.621, 0.624, 0.628, 0.631, 0.635, 0.638, 0.641, 0.644, 0.646, 0.649, 0.652, 0.654, 0.656, 0.659, 0.661, 0.663, 0.665, 0.666, 0.669, 0.670, 0.672, 0.674, 0.676, 0.678, 0.679, 0.683, 0.683, 0.685, 0.686, 0.686, 0.692, 0.693, 0.690, 0.691, 0.690, 0.692, 0.701, 0.701, 0.700, 0.704, 0.701, 0.704, 0.709, 0.698, 0.707, 0.711, 0.692, 0.739, 0.689, 0.715, 0.691, 0.720, 0.790, 0.670, 0.000, 0.717, 0.833, 0.716, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.619, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};


// MC 7 TeV 10f6 (Phojet) (Pt-Range: 0.15 - 2) ---> TPC standalone <------> values for MC ESD <---
// Double_t nbinsMC7[100] = {0.000, 0.483, 0.490, 0.502, 0.515, 0.527, 0.535, 0.541, 0.545, 0.548, 0.552, 0.555, 0.558, 0.561, 0.563, 0.566, 0.568, 0.571, 0.573, 0.574, 0.576, 0.578, 0.579, 0.581, 0.582, 0.583, 0.584, 0.585, 0.585, 0.586, 0.586, 0.587, 0.588, 0.589, 0.588, 0.588, 0.588, 0.588, 0.587, 0.589, 0.588, 0.588, 0.587, 0.587, 0.586, 0.583, 0.588, 0.587, 0.589, 0.583, 0.586, 0.582, 0.584, 0.583, 0.585, 0.582, 0.587, 0.576, 0.580, 0.588, 0.566, 0.581, 0.594, 0.585, 0.531, 0.000, 0.580, 0.499, 0.534, 0.534, 0.580, 0.432, 0.674, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.655, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};
// MC 7 TeV 10f6 (Phojet) (Pt-Range: 0.15 - 2) ---> values for MC truth <---
// Double_t nbinsMC7[100] = {0.000, 0.479, 0.483, 0.493, 0.506, 0.518, 0.527, 0.533, 0.537, 0.540, 0.543, 0.546, 0.549, 0.551, 0.554, 0.556, 0.558, 0.560, 0.562, 0.564, 0.566, 0.568, 0.569, 0.570, 0.571, 0.572, 0.573, 0.574, 0.575, 0.576, 0.576, 0.576, 0.577, 0.577, 0.577, 0.579, 0.577, 0.577, 0.577, 0.576, 0.577, 0.577, 0.577, 0.575, 0.575, 0.575, 0.573, 0.575, 0.571, 0.575, 0.568, 0.569, 0.570, 0.576, 0.566, 0.570, 0.568, 0.575, 0.565, 0.562, 0.564, 0.545, 0.563, 0.527, 0.551, 0.540, 0.571, 0.561, 0.562, 0.553, 0.000, 0.586, 0.000, 0.603, 0.480, 0.489, 0.566, 0.513, 0.715, 0.000, 0.000, 0.000, 0.441, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};


// Data 900 GeV 10c.pass3 (Pt-Range: 0.15 - 2) ---> TPC standalone <---
Double_t nbinsData900[100] = {0.000, 0.488, 0.486, 0.489, 0.495, 0.504, 0.513, 0.521, 0.528, 0.535, 0.541, 0.545, 0.549, 0.554, 0.558, 0.561, 0.565, 0.568, 0.570, 0.573, 0.576, 0.577, 0.580, 0.583, 0.582, 0.587, 0.587, 0.589, 0.591, 0.594, 0.601, 0.599, 0.598, 0.601, 0.608, 0.609, 0.600, 0.606, 0.611, 0.605, 0.612, 0.603, 0.576, 0.606, 0.617, 0.609, 0.572, 0.607, 0.546, 0.000, 0.604, 0.000, 0.689, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

// Data 2.76 TeV 11a.pass2 w/oSDD (Pt-Range: 0.15 - 2) ---> TPC standalone <---
Double_t nbinsData276[100] = {0.000, 0.492, 0.491, 0.496, 0.503, 0.513, 0.523, 0.533, 0.541, 0.549, 0.555, 0.560, 0.565, 0.570, 0.574, 0.578, 0.582, 0.585, 0.589, 0.591, 0.594, 0.597, 0.599, 0.601, 0.604, 0.606, 0.608, 0.610, 0.612, 0.613, 0.615, 0.618, 0.619, 0.621, 0.622, 0.624, 0.625, 0.626, 0.626, 0.628, 0.632, 0.633, 0.633, 0.635, 0.635, 0.638, 0.634, 0.641, 0.642, 0.644, 0.638, 0.646, 0.644, 0.647, 0.644, 0.646, 0.640, 0.650, 0.650, 0.655, 0.666, 0.650, 0.659, 0.644, 0.646, 0.657, 0.629, 0.626, 0.602, 0.652, 0.629, 0.000, 0.659, 0.580, 0.000, 0.000, 0.000, 0.000, 0.000, 0.706, 0.000, 0.596, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

// Data 7 TeV 10d.pass2 + 10e.pass2 (Pt-Range: 0.15 - 2) ---> TPC standalone <---
Double_t nbinsData7[100] = {0.000, 0.492, 0.492, 0.497, 0.505, 0.516, 0.527, 0.538, 0.547, 0.555, 0.562, 0.568, 0.573, 0.578, 0.583, 0.587, 0.591, 0.595, 0.598, 0.601, 0.604, 0.607, 0.610, 0.612, 0.615, 0.617, 0.619, 0.621, 0.623, 0.625, 0.627, 0.629, 0.631, 0.632, 0.634, 0.636, 0.637, 0.639, 0.640, 0.642, 0.643, 0.645, 0.646, 0.647, 0.648, 0.650, 0.650, 0.652, 0.653, 0.654, 0.655, 0.657, 0.657, 0.659, 0.659, 0.659, 0.660, 0.663, 0.663, 0.662, 0.663, 0.665, 0.668, 0.665, 0.667, 0.668, 0.672, 0.671, 0.669, 0.672, 0.672, 0.671, 0.675, 0.666, 0.683, 0.670, 0.673, 0.676, 0.684, 0.672, 0.673, 0.675, 0.689, 0.678, 0.678, 0.720, 0.675, 0.690, 0.664, 0.622, 0.000, 0.698, 0.678, 0.686, 0.751, 0.657, 0.000, 0.663, 0.000, 0.677};

// -- End mean pT values for the multiplicity bins --


// -- Selection of MC/Data and energy; whole sample values --

if (fMC) { // - MC -

  if (energy < 1000.) {

//     Printf(" -- MC, 900 GeV -- ");

    evMptSample = 0.5281; // MC 900 GeV 10e13 (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2) ---> TPC standalone <---

    nbins = nbinsMC900;
  }
  else if (energy > 1000. && energy < 3000.) {

//     Printf(" -- MC, 2.76 TeV -- ");

    evMptSample = 0.5574; // MC 2.76 TeV 11b10a (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2) ---> TPC standalone <---

    nbins = nbinsMC276;
  }
  else {

//     Printf(" -- MC, 7 TeV -- ");

    evMptSample = 0.5809; // MC 7 TeV 10f6a (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2) ---> TPC standalone <---
//     evMptSample = 0.5732; // MC 7 TeV 10f6a (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2) ---> value for MC truth <---

//     evMptSample = 0.5481; // MC 7 TeV 10f6 (Phojet) (Pt-Range: 0.15 - 2) ---> TPC standalone <---
//     evMptSample = 0.5409; // MC 7 TeV 10f6 (Phojet) (Pt-Range: 0.15 - 2) ---> value for MC truth <---

    nbins = nbinsMC7;
  }

} // - End MC -
else { // - Data -

  if (energy < 1000.) {

//     Printf(" -- Data, 900 GeV -- ");

    evMptSample = 0.5263; // Data 900 GeV 10c.pass3 (Pt-Range: 0.15 - 2) ---> TPC standalone <---

    nbins = nbinsData900;
  }
  else if (energy > 1000. && energy < 3000.) {

//     Printf(" -- Data, 2.76 TeV -- ");

    evMptSample = 0.5548; // Data 2.76 TeV 11a.pass2 w/oSDD (Pt-Range: 0.15 - 2) ---> TPC standalone <---

    nbins = nbinsData276;
  }
  else {

//     Printf(" -- Data, 7 TeV -- ");

    evMptSample = 0.5738; // Data 7 TeV 10d.pass2 + 10e.pass2 (Pt-Range: 0.15 - 2) ---> TPC standalone <---

    nbins = nbinsData7;
  }

} // - End data -

// -- End selection of MC/Data and energy; whole sample values --

// --- End mean pT values ---


// --- Ratio of pT distributions MC ESD / MC truth ---

// MC 7 TeV 10f6a (Pythia6 - Perugia0) (Pt-Range: 0.15 - 2, array-range: 0 - 2.1 GeV/c)
Double_t nbinsPtRatio[84] = {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.823, 0.908, 0.908, 0.901, 0.903, 0.910, 0.918, 0.931, 0.941, 0.947, 0.951, 0.953, 0.956, 0.959, 0.959, 0.960, 0.961, 0.962, 0.963, 0.963, 0.963, 0.964, 0.964, 0.964, 0.963, 0.963, 0.963, 0.964, 0.964, 0.962, 0.963, 0.964, 0.963, 0.963, 0.962, 0.962, 0.961, 0.963, 0.964, 0.964, 0.965, 0.965, 0.965, 0.964, 0.964, 0.964, 0.964, 0.964, 0.963, 0.964, 0.961, 0.957, 0.958, 0.959, 0.959, 0.958, 0.958, 0.957, 0.959, 0.956, 0.956, 0.954, 0.954, 0.956, 0.955, 0.954, 0.954, 0.954, 0.955, 0.954, 0.951, 0.952, 0.951, 0.949, 0.000, 0.000, 0.000, 0.000};

// --- End ratio of pT distributions MC ESD / MC truth ---


  if (!fMC || fMCType < 3) {
    nrESDTracks = fESD->GetNumberOfTracks();
  }

  if (fMC) {
    nrMCTracks = stack->GetNtrack();
  }

//   Printf("\n\n nrESDTracks: %i   nrMCTracks : %i \n",nrESDTracks,nrMCTracks);


if (fMC) { // - MC -

 if (fMCType == 0) { // MC, Type = ESD


  // --- Loop over all tracks of one ESD event ---
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

// 	Printf("Track Pt = %.3f   Track Eta = %.3f   Track Phi = %3f\n",trackPt,trackEta,trackPhi);

  } // --- End track loop ---


 } // End MC, Type = ESD
 else if (fMCType > 0) { // MC, Type = MC truth


  // --- Loop over all tracks of one MC event - MC truth / stack ---
  for (Int_t iTracks = 0; iTracks < nrMCTracks; iTracks++) {

    TParticle *track = stack->Particle(iTracks);

    if (!track) {
      Printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }

	// only charged particles
	TParticlePDG* pdg = track->GetPDG();
	if(!pdg) continue;

	Double_t charge = pdg->Charge()/3.;
	if ( TMath::Abs(charge) < 0.001 ) continue;

	// only physical primary
	Bool_t prim = stack->IsPhysicalPrimary(iTracks);
	if(!prim) continue;

	// get pt, eta and phi and cut in eta and pt
	trackPt = track->Pt();
	trackEta = track->Eta();
	trackPhi = track->Phi();

// 	Printf("Track Pt = %.3f   Track Eta = %.3f   Track Phi = %3f  ",trackPt,trackEta,trackPhi);

	if (trackEta > 0.8) continue;
	if (trackEta < -0.8) continue;

	if (trackPt < 0.15) continue;
	if (trackPt > 2.) continue;

// 	Printf("Track Pt = %.3f   Track Eta = %.3f   Track Phi = %3f  ",trackPt,trackEta,trackPhi);
// 	Printf("-> Track accepted!");


	// -- Reject tracks according to ratio of pT distribution MC ESD / MC truth --
	if (fMCType == 2 || fMCType == 4) { // MC, Type = mod. MC truth

	  trackPt40 = trackPt*40.;
	  trackPtBin = floor(trackPt40);

	  ptRatio = nbinsPtRatio[trackPtBin];

// 	  Printf("trackPt: %f   trackPt40: %f   trackPtBin: %i   ptRatio: %.3f   ",trackPt,trackPt40,trackPtBin,ptRatio);

	  random = fRandom3->Rndm();
// 	  Printf("Random: %f  \n",random);

	  if (random > ptRatio) continue;
// 	  Printf("accepted random: %f \n",random);

	} // -- End reject tracks according to ratio of pT distribution MC ESD / MC truth --


	fPtSpec->Fill(trackPt);
	tracks[nrTracks] = trackPt;

	fEta->Fill(trackEta);

	if (trackEta > 0.) {
	  fEtaPhiPlus->Fill(trackPhi);
	}
	else if (trackEta < 0.) {
	  fEtaPhiMinus->Fill(trackPhi);
	}

	sumPt = sumPt + trackPt;
	nrTracks++;

// 	Printf("Track Pt = %.3f   Track Eta = %.3f   Track Phi = %3f\n",trackPt,trackEta,trackPhi);

  } // --- End stack track loop ---

//   Printf("nrTracks MC truth: %i \n",nrTracks);


  // MC, Type = MC truth with corresponding ESD
  if (fMCType < 3) {

    // --- Loop over all tracks of one MC event - MC ESD ---
    for (Int_t jTracks = 0; jTracks < nrESDTracks; jTracks++) {
      AliESDtrack* track = fESD->GetTrack(jTracks);
      if (!track) {
        Printf("ERROR: Could not receive track %d\n", jTracks);
        continue;
      }

        if(!fESDTrackCuts->AcceptTrack(track))continue;

	  trackPt2 = track->Pt();
	  fPtSpec2->Fill(trackPt2);

// 	  nrTracks2++;

// 	  Printf("Track Pt ESD = %.3f\n\n",trackPt2);

    } // --- End ESD track loop ---

//     Printf("nrTracks ESD: %i \n\n",nrTracks2);

  } // End MC, Tpye = MC truth with corresponding ESD

 } // End MC, Type = MC truth


} // - End MC -
else { // - Data -


  // --- Loop over all tracks of one ESD event ---
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

// 	Printf("Track Pt = %.3f   Track Eta = %.3f   Track Phi = %3f\n",trackPt,trackEta,trackPhi);

  } // --- End track loop ---

} // - End data -


  // --- Calculation of various values and filling of histograms (for remaining events with N_acc > 0) ---

  if(nrTracks != 0) {

    if (!fMC || fMCType < 3) {

      // VtxZ after all track cuts
      fVtxZTrackCuts->Fill(vtxZ);

      // Differences between vertices after all cuts
      fVtxZDiff1b->Fill(vtxZdiff1);
      fVtxZDiff2b->Fill(vtxZdiff2);
      fVtxZDiff3b->Fill(vtxZdiff3);

    }

    // Multiplicity distribution
    fMult->Fill(nrTracks);

//     // Number of pairs in event
//     nrPairs = 0.5 * nrTracks * (nrTracks-1);

    // Calculation of mean Pt and mean Pt Squared
    eventMeanPt = sumPt / nrTracks;
    eventMeanPtSq = eventMeanPt * eventMeanPt;

    // Mean-Pt and Mean-Pt squared
    fEventMeanPt->Fill(eventMeanPt);
    fEventMeanPtSq->Fill(eventMeanPtSq);
    fEventMeanPtMult->Fill(nrTracks,eventMeanPt);

    // Mean-Pt and Mean-Pt squared depending on multiplicity
    fMultEventMeanPt->Fill(nrTracks,eventMeanPt);
    fMultEventMeanPtSq->Fill(nrTracks,eventMeanPtSq);

//     Printf("nrTracks: %i   sumPt: %.8f   meanPt: %.8f   meanPtSq: %.8f\n",nrTracks,sumPt,eventMeanPt,eventMeanPtSq);


    // --- Two-particle correlator ---

    evMptMult = nbins[nrTracks];
//     Printf("nrTracks = %3i   evMptMult = %.3f   evMptSample = %.3f \n\n",nrTracks,evMptMult,evMptSample);

    for (int i=0; i<nrTracks; i++) {

      for (int j=i+1; j<nrTracks; j++) {

	twoPartCorrPair = (tracks[i] - evMptMult) * (tracks[j] - evMptMult); // Multiplicity bins
	twoPartCorrEv = twoPartCorrEv + twoPartCorrPair;

	twoPartCorrPairSample = (tracks[i] - evMptSample) * (tracks[j] - evMptSample); // Whole sample
	twoPartCorrEvSample = twoPartCorrEvSample + twoPartCorrPairSample;

      }
    }

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

  // Delet pointer
  if (fRandom3) delete fRandom3; fRandom3 = 0;

}

void AliAnalysisTaskPtFluc::Terminate(Option_t *)
{
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    Printf("ERROR: fOutputList not available\n");
    return;
  }

}
