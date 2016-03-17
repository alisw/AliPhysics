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
#include "AliCentrality.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskPtFlucPbPb.h"


using namespace std;

// Analysis of Pt Fluctuations (PbPb)
// Author: Stefan Heckel
// Version of PbPb task:  9.2, 03.07.2012


ClassImp(AliAnalysisTaskPtFlucPbPb)

//________________________________________________________________________
AliAnalysisTaskPtFlucPbPb::AliAnalysisTaskPtFlucPbPb(const char *name)
  :AliAnalysisTaskSE(name),
  fESD(0),
  fMCev(0),
  fRandom3(0),
  fOutputList(0),
  fPtSpec(0),
  fPtSpec2(0),
  fMult(0),
  fMultNbins(0),
  fMultSum(0),
  fMultSumPt(0),
  fMultNrPairs(0),
  fMult1(0),
  fMultSum1(0),
  fMultSumPt1(0),
  fMultNrPairs1(0),
  fMult10(0),
  fMultSum10(0),
  fMultSumPt10(0),
  fMultNrPairs10(0),
  fMult80(0),
  fMultSum80(0),
  fMultSumPt80(0),
  fMultNrPairs80(0),
  fMult801(0),
  fMultSum801(0),
  fMultSumPt801(0),
  fMultNrPairs801(0),
  fMult810(0),
  fMultSum810(0),
  fMultSumPt810(0),
  fMultNrPairs810(0),
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
  fVtxZCutDiff(0),
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
  fMultEventMeanPtNbins(0),
  fMultEventMeanPtSqNbins(0),
  fCentEventMeanPt(0),
  fCentEventMeanPtSq(0),
  fEventMeanPtCent05(0),
  fEventMeanPtCent2030(0),
  fEventMeanPtCent7080(0),
  fTwoPartCorrEv(0),
  fTwoPartCorrEvSq(0),
  fTwoPartCorrEv1(0),
  fTwoPartCorrEvSq1(0),
  fTwoPartCorrEv10(0),
  fTwoPartCorrEvSq10(0),
  fTwoPartCorrEv80(0),
  fTwoPartCorrEvSq80(0),
  fTwoPartCorrEv801(0),
  fTwoPartCorrEvSq801(0),
  fTwoPartCorrEv810(0),
  fTwoPartCorrEvSq810(0),
  fTwoPartCorrEvCent(0),
  fTwoPartCorrEvCentSq(0),
  fESDTrackCuts(0),
  fMaxVertexZ(0),
  fMaxVertexZDiff1(0),
  fNContributors(0),
  fUseCentrality(0),
  fMC(0),
  fMCType(0),
  fMCAMPT(0)
{
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskPtFlucPbPb::~AliAnalysisTaskPtFlucPbPb()
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
void AliAnalysisTaskPtFlucPbPb::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(1, "RECREATE");
  fOutputList = new TList();
  fOutputList->SetOwner();


  fPtSpec = new TH1F("fPtSpec","Pt spectrum",100,0,2.5);
  fPtSpec2 = new TH1F("fPtSpec2","Pt spectrum 2 - MC ESD",100,0,2.5);

  fMult = new TH1F("fMult","Multiplicity distribution",30,0,3000);
  fMultNbins = new TH1F("fMultNbins","Multiplicity distribution",3000,0,3000);
  fMultSum = new TH1F("fMultSum","Sum of nrtracks of all events in multiplicity bins",30,0,3000);
  fMultSumPt = new TH1F("fMultSumPt","Sum of pTs of all events in multiplicity bins",30,0,3000);
  fMultNrPairs = new TH1F("fMultNrPairs","Sum of number of pairs in multiplicity bins",30,0,3000);

  fMult1 = new TH1F("fMult1","Multiplicity distribution",4,0,100);
  fMultSum1 = new TH1F("fMultSum1","Sum of nrtracks of all events in multiplicity bins",4,0,100);
  fMultSumPt1 = new TH1F("fMultSumPt1","Sum of pTs of all events in multiplicity bins",4,0,100);
  fMultNrPairs1 = new TH1F("fMultNrPairs1","Sum of number of pairs in multiplicity bins",4,0,100);

  fMult10 = new TH1F("fMult10","Multiplicity distribution",5,0,50);
  fMultSum10 = new TH1F("fMultSum10","Sum of nrtracks of all events in multiplicity bins",5,0,50);
  fMultSumPt10 = new TH1F("fMultSumPt10","Sum of pTs of all events in multiplicity bins",5,0,50);
  fMultNrPairs10 = new TH1F("fMultNrPairs10","Sum of number of pairs in multiplicity bins",5,0,50);

  fMult80 = new TH1F("fMult80","Multiplicity distribution",30,0,3000);
  fMultSum80 = new TH1F("fMultSum80","Sum of nrtracks of all events in multiplicity bins",30,0,3000);
  fMultSumPt80 = new TH1F("fMultSumPt80","Sum of pTs of all events in multiplicity bins",30,0,3000);
  fMultNrPairs80 = new TH1F("fMultNrPairs80","Sum of number of pairs in multiplicity bins",30,0,3000);

  fMult801 = new TH1F("fMult801","Multiplicity distribution",4,0,100);
  fMultSum801 = new TH1F("fMultSum801","Sum of nrtracks of all events in multiplicity bins",4,0,100);
  fMultSumPt801 = new TH1F("fMultSumPt801","Sum of pTs of all events in multiplicity bins",4,0,100);
  fMultNrPairs801 = new TH1F("fMultNrPairs801","Sum of number of pairs in multiplicity bins",4,0,100);

  fMult810 = new TH1F("fMult810","Multiplicity distribution",5,0,50);
  fMultSum810 = new TH1F("fMultSum810","Sum of nrtracks of all events in multiplicity bins",5,0,50);
  fMultSumPt810 = new TH1F("fMultSumPt810","Sum of pTs of all events in multiplicity bins",5,0,50);
  fMultNrPairs810 = new TH1F("fMultNrPairs810","Sum of number of pairs in multiplicity bins",5,0,50);


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
  fVtxZCutDiff = new TH1F("fVtxZCutDiff","Vertex Z distribution after cut on vtx Z Diff",110,-11,11);
  fVtxZTrackCuts = new TH1F("fVtxZTrackCuts","Vertex Z distribution after track cuts",110,-11,11);

  fVtxZDiff1 = new TH1F("fVtxZDiff1","Difference 1 between vertex Z distributions",100,-5,5);
  fVtxZDiff2 = new TH1F("fVtxZDiff2","Difference 2 between vertex Z distributions",100,-5,5);
  fVtxZDiff3 = new TH1F("fVtxZDiff3","Difference 3 between vertex Z distributions",100,-5,5);
  fVtxZDiff1b = new TH1F("fVtxZDiff1b","Difference 1 between vertex Z distributions after all cuts",100,-5,5);
  fVtxZDiff2b = new TH1F("fVtxZDiff2b","Difference 2 between vertex Z distributions after all cuts",100,-5,5);
  fVtxZDiff3b = new TH1F("fVtxZDiff3b","Difference 3 between vertex Z distributions after all cuts",100,-5,5);


  fEventMeanPt = new TH1F("fEventMeanPt","Mean-Pt distribution",250,0,2.5);
  fEventMeanPtSq = new TH1F("fEventMeanPtSq","Mean-Pt squared distribution",500,0,5);
  fEventMeanPtMult = new TH2F("fEventMeanPtMult","Mean-Pt for single events vs. multiplicity",30,0,3000,200,0.,2.);

  fMultEventMeanPt = new TH1F("fMultEventMeanPt","Mean-Pt event by event vs. multiplicity",30,0,3000);
  fMultEventMeanPtSq = new TH1F("fMultEventMeanPtSq","Mean-Pt event by event squared vs. multiplicity",30,0,3000);

  fMultEventMeanPtNbins = new TH1F("fMultEventMeanPtNbins","Mean-Pt event by event vs. multiplicity",3000,0,3000);
  fMultEventMeanPtSqNbins = new TH1F("fMultEventMeanPtSqNbins","Mean-Pt event by event squared vs. multiplicity",3000,0,3000);

  fCentEventMeanPt = new TH1F("fCentEventMeanPt","Mean-Pt event by event vs. centrality",21,0,105);
  fCentEventMeanPtSq = new TH1F("fCentEventMeanPtSq","Mean-Pt event by event squared vs. centrality",21,0,105);

  fEventMeanPtCent05 = new TH1F("fEventMeanPtCent05","Mean-Pt distribution for centrality 0-5%",500,0.4,0.9);
  fEventMeanPtCent2030 = new TH1F("fEventMeanPtCent2030","Mean-Pt distribution for centrality 20-30%",500,0.4,0.9);
  fEventMeanPtCent7080 = new TH1F("fEventMeanPtCent7080","Mean-Pt distribution for centrality 70-80%",1500,0,1.5);


  fTwoPartCorrEv = new TH1F("fTwoPartCorrEv","Two-particle correlator vs. multiplicity",30,0,3000);
  fTwoPartCorrEvSq = new TH1F("fTwoPartCorrEvSq","Two-particle correlator squared vs. multiplicity",30,0,3000);

  fTwoPartCorrEv1 = new TH1F("fTwoPartCorrEv1","Two-particle correlator vs. multiplicity",4,0,100);
  fTwoPartCorrEvSq1 = new TH1F("fTwoPartCorrEvSq1","Two-particle correlator squared vs. multiplicity",4,0,100);

  fTwoPartCorrEv10 = new TH1F("fTwoPartCorrEv10","Two-particle correlator vs. multiplicity",5,0,50);
  fTwoPartCorrEvSq10 = new TH1F("fTwoPartCorrEvSq10","Two-particle correlator squared vs. multiplicity",5,0,50);

  fTwoPartCorrEv80 = new TH1F("fTwoPartCorrEv80","Two-particle correlator vs. multiplicity",30,0,3000);
  fTwoPartCorrEvSq80 = new TH1F("fTwoPartCorrEvSq80","Two-particle correlator squared vs. multiplicity",30,0,3000);

  fTwoPartCorrEv801 = new TH1F("fTwoPartCorrEv801","Two-particle correlator vs. multiplicity",4,0,100);
  fTwoPartCorrEvSq801 = new TH1F("fTwoPartCorrEvSq801","Two-particle correlator squared vs. multiplicity",4,0,100);

  fTwoPartCorrEv810 = new TH1F("fTwoPartCorrEv810","Two-particle correlator vs. multiplicity",5,0,50);
  fTwoPartCorrEvSq810 = new TH1F("fTwoPartCorrEvSq810","Two-particle correlator squared vs. multiplicity",5,0,50);

  fTwoPartCorrEvCent = new TH1F("fTwoPartCorrEvCent","Two-particle correlator vs. centrality",21,0,105);
  fTwoPartCorrEvCentSq = new TH1F("fTwoPartCorrEvCentSq","Two-particle correlator squared vs. centrality",21,0,105);


  // Add histograms to the output list
  fOutputList->Add(fPtSpec);
  fOutputList->Add(fPtSpec2);
  fOutputList->Add(fMult);
  fOutputList->Add(fMultNbins);
  fOutputList->Add(fMultSum);
  fOutputList->Add(fMultSumPt);
  fOutputList->Add(fMultNrPairs);
  fOutputList->Add(fMult1);
  fOutputList->Add(fMultSum1);
  fOutputList->Add(fMultSumPt1);
  fOutputList->Add(fMultNrPairs1);
  fOutputList->Add(fMult10);
  fOutputList->Add(fMultSum10);
  fOutputList->Add(fMultSumPt10);
  fOutputList->Add(fMultNrPairs10);
  fOutputList->Add(fMult80);
  fOutputList->Add(fMultSum80);
  fOutputList->Add(fMultSumPt80);
  fOutputList->Add(fMultNrPairs80);
  fOutputList->Add(fMult801);
  fOutputList->Add(fMultSum801);
  fOutputList->Add(fMultSumPt801);
  fOutputList->Add(fMultNrPairs801);
  fOutputList->Add(fMult810);
  fOutputList->Add(fMultSum810);
  fOutputList->Add(fMultSumPt810);
  fOutputList->Add(fMultNrPairs810);
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
  fOutputList->Add(fVtxZCutDiff);
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
  fOutputList->Add(fMultEventMeanPtNbins);
  fOutputList->Add(fMultEventMeanPtSqNbins);
  fOutputList->Add(fCentEventMeanPt);
  fOutputList->Add(fCentEventMeanPtSq);
  fOutputList->Add(fEventMeanPtCent05);
  fOutputList->Add(fEventMeanPtCent2030);
  fOutputList->Add(fEventMeanPtCent7080);
  fOutputList->Add(fTwoPartCorrEv);
  fOutputList->Add(fTwoPartCorrEvSq);
  fOutputList->Add(fTwoPartCorrEv1);
  fOutputList->Add(fTwoPartCorrEvSq1);
  fOutputList->Add(fTwoPartCorrEv10);
  fOutputList->Add(fTwoPartCorrEvSq10);
  fOutputList->Add(fTwoPartCorrEv80);
  fOutputList->Add(fTwoPartCorrEvSq80);
  fOutputList->Add(fTwoPartCorrEv801);
  fOutputList->Add(fTwoPartCorrEvSq801);
  fOutputList->Add(fTwoPartCorrEv810);
  fOutputList->Add(fTwoPartCorrEvSq810);
  fOutputList->Add(fTwoPartCorrEvCent);
  fOutputList->Add(fTwoPartCorrEvCentSq);

  // Post output data (if histograms are not used later, PostData is at least called here)
  PostData(1, fOutputList);

}

//________________________________________________________________________
void AliAnalysisTaskPtFlucPbPb::UserExec(Option_t *)
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


  Int_t nrTracks = 0;		// count for number of tracks which pass the track cuts
//   Int_t nrTracks2 = 0;		// count for number of tracks which pass the track cuts - for MC ESD vs. MC truth
  Int_t nrESDTracks = 0;	// number of tracks in the ESD file
  Int_t nrMCTracks = 0;		// number of tracks in the MC file
  Double_t sumPt = 0.;		// sum of track Pts
  Double_t twoPartCorrEv = 0.;	// Two-particle correlator of one event for multiplicity bin analysis
  Double_t twoPartCorrEvCent = 0.; // Two-particle correlator of one event for centrality bin analysis

  Double_t tracks[3000] = {0.};	// array of track Pts, needed for the two-particle correlator

  Double_t trackPt=0., trackEta=0., trackPhi=0.;
  Double_t trackPt2=0.; // for MC ESD
  Double_t eventMeanPt=0., eventMeanPtSq=0., evMptMult=0.;
  Double_t twoPartCorrPair=0., twoPartCorrEvSq=0.;
  Double_t twoPartCorrPairCent=0., twoPartCorrEvCentSq=0.;
  Double_t nrPairs=0.;

  Double_t *nbins = 0x0; // Mean pT values for multiplicity bin analysis
  Double_t *centbins = 0x0; // Mean pT values for centrality bin analysis
  Double_t *sNbinsMC276 = 0x0; // Mean pT values for mult. with Delta N_acc = 1, MC

  Double_t centralityVZERO=0.;
  Int_t centralityVZEROBin=0;

  Int_t centBin=0;
  Double_t evMptCent=0.;

  Int_t maxNrTracksMC=0, minNrBinMC=0; // special values, which have to be set differently for AMPT;

  Double_t random=0., trackPt40=0., ptRatio=0.;
  Int_t trackPtBin=0;

  fRandom3 = new TRandom3();

  if (fMCType == 2 || fMCType == 4) { // MC, Type = mod. MC truth
   fRandom3->SetSeed(0);
  }


// --- Mean pT values ---

// !!!!! Have to be calculated in a first run - for each sample, which should be analysed, separately! !!!!!


// -- Mean pT values for multiplicity bins (first bin is for nrTracks = 0 and has to be 0) --

  // MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2) ---> values for MC ESD <---
  Double_t nbinsMC276H[32] = {0.0000, 0.5251, 0.5350, 0.5385, 0.5405, 0.5420, 0.5435, 0.5444, 0.5452, 0.5459, 0.5467, 0.5473, 0.5478, 0.5485, 0.5490, 0.5496, 0.5500, 0.5505, 0.5513, 0.5519, 0.5525, 0.5531, 0.5540, 0.5547, 0.5554, 0.5560, 0.5570, 0.5581, 0.5578, 0.5693, 0.0000, 0.0000};
  // MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2) ---> values for MC truth <---
//   Double_t nbinsMC276H[32] = {0.000, 0.515, 0.524, 0.528, 0.530, 0.531, 0.532, 0.533, 0.534, 0.534, 0.535, 0.535, 0.536, 0.536, 0.536, 0.536, 0.537, 0.537, 0.537, 0.537, 0.538, 0.538, 0.538, 0.539, 0.539, 0.539, 0.540, 0.541, 0.541, 0.541, 0.543, 0.543};
  // MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2) ---> values for mod. MC truth <---
//   Double_t nbinsMC276H[32] = {0.000, 0.528, 0.537, 0.540, 0.542, 0.544, 0.545, 0.546, 0.546, 0.547, 0.547, 0.548, 0.548, 0.548, 0.549, 0.549, 0.549, 0.550, 0.550, 0.550, 0.551, 0.551, 0.552, 0.552, 0.552, 0.553, 0.554, 0.554, 0.554, 0.558, 0.565, 0.000};


  // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) ---> values for MC ESD <---
//   Double_t nbinsMC276A[32] = {0.0000, 0.5274, 0.5365, 0.5398, 0.5433, 0.5449, 0.5440, 0.5416, 0.5398, 0.5383, 0.5371, 0.5353, 0.5340, 0.5324, 0.5304, 0.5293, 0.5284, 0.5276, 0.5267, 0.5258, 0.5249, 0.5242, 0.5238, 0.5239, 0.5241, 0.5247, 0.5258, 0.5269, 0.0000, 0.0000, 0.0000, 0.0000};
  // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) ---> values for MC ESD <------> merged NEW <---
  Double_t nbinsMC276A[32] = {0.0000, 0.5275, 0.5380, 0.5418, 0.5432, 0.5429, 0.5421, 0.5411, 0.5398, 0.5383, 0.5371, 0.5354, 0.5341, 0.5327, 0.5310, 0.5299, 0.5290, 0.5279, 0.5269, 0.5260, 0.5253, 0.5245, 0.5242, 0.5241, 0.5246, 0.5252, 0.5266, 0.5266, 0.0000, 0.0000, 0.0000, 0.0000};

  // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) ---> values for MC truth <---
//   Double_t nbinsMC276A[32] = {0.0000, 0.5072, 0.5179, 0.5219, 0.5234, 0.5237, 0.5227, 0.5221, 0.5205, 0.5189, 0.5177, 0.5158, 0.5146, 0.5127, 0.5114, 0.5099, 0.5090, 0.5075, 0.5061, 0.5055, 0.5045, 0.5038, 0.5038, 0.5045, 0.5037, 0.5059, 0.5032, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
  // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) ---> values for mod. MC truth <---
//   Double_t nbinsMC276A[32] = {0.0000, 0.5205, 0.5312, 0.5350, 0.5363, 0.5359, 0.5353, 0.5337, 0.5318, 0.5306, 0.5287, 0.5272, 0.5254, 0.5237, 0.5222, 0.5211, 0.5196, 0.5185, 0.5174, 0.5166, 0.5170, 0.5174, 0.5162, 0.5191, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};


  // Data PbPb 2.76 ATeV 10h.pass2 (Pt-Range: 0.15 - 2)
  Double_t nbinsData276[32] = {0.0000, 0.5785, 0.6037, 0.6159, 0.6240, 0.6300, 0.6345, 0.6381, 0.6409, 0.6432, 0.6452, 0.6467, 0.6480, 0.6492, 0.6501, 0.6508, 0.6515, 0.6521, 0.6525, 0.6528, 0.6530, 0.6532, 0.6532, 0.6533, 0.6538, 0.6550, 0.6563, 0.6579, 0.6595, 0.6599, 0.6607, 0.6549};

// -- End mean pT values for multiplicity bins --


// -- Mean pT values for centrality bins --

// - 5% centrality bins -

  // MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2) ---> values for MC ESD <---
  Double_t centbinsMC276H[21] = {0.5541, 0.5517, 0.5499, 0.5484, 0.5473, 0.5460, 0.5448, 0.5436, 0.5419, 0.5404, 0.5389, 0.5369, 0.5350, 0.5320, 0.5299, 0.5268, 0.5232, 0.5179, 0.5150, 0.0000, 0.5072};
  // MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2) ---> values for MC truth <---
//   Double_t centbinsMC276H[21] = {0.538, 0.537, 0.536, 0.535, 0.534, 0.534, 0.533, 0.531, 0.530, 0.529, 0.527, 0.526, 0.523, 0.521, 0.519, 0.516, 0.512, 0.507, 0.504, 0.000, 0.538};
  // MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2) ---> values for mod. MC truth <---
//   Double_t centbinsMC276H[21] = {0.551, 0.549, 0.548, 0.547, 0.547, 0.546, 0.545, 0.543, 0.542, 0.541, 0.539, 0.537, 0.535, 0.533, 0.530, 0.528, 0.524, 0.518, 0.516, 0.000, 0.550};


  // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) ---> values for MC ESD <---
//   Double_t centbinsMC276A[21] = {0.0000, 0.5230, 0.5239, 0.5262, 0.5287, 0.5326, 0.5361, 0.5398, 0.5431, 0.5450, 0.5439, 0.5417, 0.5387, 0.5369, 0.5344, 0.5304, 0.5264, 0.5216, 0.5169, 0.0000, 0.5309};
  // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) ---> values for MC ESD <------> merged NEW <---
  Double_t centbinsMC276A[21] = {0.0000, 0.5241, 0.5242, 0.5263, 0.5294, 0.5329, 0.5362, 0.5397, 0.5422, 0.5432, 0.5434, 0.5428, 0.5407, 0.5378, 0.5346, 0.5303, 0.5258, 0.5205, 0.5143, 0.0000, 0.5316};


  // Data PbPb 2.76 ATeV 10h.pass2 (Pt-Range: 0.15 - 2)
  Double_t centbinsData276[21] = {0.6534, 0.6526, 0.6510, 0.6488, 0.6460, 0.6425, 0.6383, 0.6334, 0.6277, 0.6215, 0.6147, 0.6075, 0.6001, 0.5928, 0.5855, 0.5778, 0.5677, 0.5510, 0.5382, 0.0000, 0.5652};

// - End 5% centrality bins -


// - 10% centrality bins -

  // MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2) ---> values for MC ESD <---
//   Double_t centbinsMC276H[11] = {0.553, 0.549, 0.547, 0.544, 0.541, 0.538, 0.534, 0.529, 0.521, 0.517, 0.509};


  // Data PbPb 2.76 ATeV 10h.pass2 (Pt-Range: 0.15 - 2)
//   Double_t centbinsData276H[11] = {0.653, 0.650, 0.644, 0.636, 0.625, 0.612, 0.597, 0.582, 0.562, 0.539, 0.548};

// - End 10% centrality bins -

// -- End mean pT values for centrality bins --


// -- Mean pT values for single multiplicity bins ---> the first value has to be 0.000 <---

  // MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2) (up to N_acc = 500) ---> values for MC ESD <---
  Double_t sNbinsMC276H[501] = {0.000, 0.480, 0.481, 0.483, 0.484, 0.489, 0.493, 0.498, 0.501, 0.504, 0.506, 0.508, 0.510, 0.511, 0.512, 0.513, 0.514, 0.515, 0.516, 0.517, 0.517, 0.518, 0.519, 0.519, 0.520, 0.520, 0.521, 0.521, 0.521, 0.522, 0.522, 0.523, 0.523, 0.523, 0.523, 0.524, 0.524, 0.524, 0.525, 0.525, 0.525, 0.525, 0.525, 0.526, 0.526, 0.526, 0.526, 0.526, 0.527, 0.527, 0.527, 0.527, 0.527, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543};

  // MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2) (up to N_acc = 500) ---> values for MC truth <---
//   Double_t sNbinsMC276H[501] = {0.000, 0.477, 0.472, 0.475, 0.474, 0.481, 0.483, 0.486, 0.489, 0.495, 0.495, 0.495, 0.503, 0.504, 0.504, 0.505, 0.506, 0.506, 0.507, 0.507, 0.508, 0.508, 0.509, 0.509, 0.510, 0.510, 0.510, 0.511, 0.511, 0.512, 0.512, 0.512, 0.512, 0.513, 0.513, 0.513, 0.513, 0.514, 0.514, 0.514, 0.514, 0.515, 0.515, 0.515, 0.515, 0.515, 0.516, 0.516, 0.516, 0.516, 0.516, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.521, 0.521, 0.521, 0.521, 0.521, 0.521, 0.521, 0.521, 0.521, 0.521, 0.521, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.523, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532};


  // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) (up to N_acc = 200) ---> values for MC ESD <---
//   Double_t sNbinsMC276A[201] = {0.000, 0.474, 0.499, 0.492, 0.489, 0.489, 0.495, 0.489, 0.492, 0.495, 0.498, 0.500, 0.502, 0.503, 0.505, 0.506, 0.508, 0.509, 0.510, 0.511, 0.512, 0.513, 0.513, 0.514, 0.515, 0.516, 0.516, 0.517, 0.517, 0.518, 0.519, 0.519, 0.520, 0.520, 0.520, 0.521, 0.521, 0.522, 0.522, 0.522, 0.523, 0.523, 0.523, 0.524, 0.524, 0.524, 0.525, 0.525, 0.525, 0.525, 0.526, 0.526, 0.526, 0.526, 0.527, 0.527, 0.527, 0.527, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538};
  // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) (up to N_acc = 300) ---> values for MC ESD <------> merged NEW <---
  Double_t sNbinsMC276A[301] = {0.000, 0.476, 0.497, 0.489, 0.491, 0.490, 0.498, 0.491, 0.494, 0.496, 0.499, 0.501, 0.502, 0.504, 0.505, 0.507, 0.508, 0.509, 0.510, 0.511, 0.512, 0.513, 0.514, 0.514, 0.515, 0.516, 0.516, 0.517, 0.518, 0.518, 0.519, 0.519, 0.520, 0.520, 0.521, 0.521, 0.521, 0.522, 0.522, 0.523, 0.523, 0.523, 0.524, 0.524, 0.524, 0.525, 0.525, 0.525, 0.525, 0.526, 0.526, 0.526, 0.527, 0.527, 0.527, 0.527, 0.528, 0.528, 0.528, 0.528, 0.528, 0.529, 0.529, 0.529, 0.529, 0.529, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.535, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.537, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.539, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.540, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.541, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.542, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543, 0.543};

//   // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) (up to N_acc = 200) ---> values for MC truth <---
//   Double_t sNbinsMC276A[201] = {0.000, 0.404, 0.404, 0.426, 0.439, 0.448, 0.455, 0.461, 0.465, 0.469, 0.472, 0.475, 0.477, 0.479, 0.481, 0.483, 0.484, 0.485, 0.487, 0.488, 0.489, 0.490, 0.491, 0.492, 0.493, 0.494, 0.494, 0.495, 0.496, 0.496, 0.497, 0.498, 0.498, 0.499, 0.499, 0.500, 0.500, 0.501, 0.501, 0.501, 0.502, 0.502, 0.503, 0.503, 0.503, 0.504, 0.504, 0.504, 0.505, 0.505, 0.505, 0.506, 0.506, 0.506, 0.506, 0.507, 0.507, 0.507, 0.507, 0.508, 0.508, 0.508, 0.508, 0.508, 0.509, 0.509, 0.509, 0.509, 0.509, 0.510, 0.510, 0.510, 0.510, 0.510, 0.511, 0.511, 0.511, 0.511, 0.511, 0.511, 0.511, 0.512, 0.512, 0.512, 0.512, 0.512, 0.512, 0.512, 0.513, 0.513, 0.513, 0.513, 0.513, 0.513, 0.513, 0.513, 0.514, 0.514, 0.514, 0.514, 0.514, 0.514, 0.514, 0.514, 0.514, 0.515, 0.515, 0.515, 0.515, 0.515, 0.515, 0.515, 0.515, 0.515, 0.515, 0.516, 0.516, 0.516, 0.516, 0.516, 0.516, 0.516, 0.516, 0.516, 0.516, 0.516, 0.516, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.517, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.518, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.519, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.520, 0.521, 0.521, 0.521, 0.521, 0.521};

  // MC PbPb 2.76 ATeV AMPT (Pt-Range: 0.15 - 2) (up to N_acc = 200) ---> values for mod. MC truth <---
//   Double_t sNbinsMC276A[201] = {0.000, 0.479, 0.480, 0.481, 0.482, 0.482, 0.483, 0.485, 0.486, 0.488, 0.490, 0.491, 0.493, 0.494, 0.496, 0.497, 0.498, 0.499, 0.500, 0.502, 0.502, 0.503, 0.504, 0.505, 0.506, 0.507, 0.507, 0.508, 0.509, 0.509, 0.510, 0.511, 0.511, 0.512, 0.512, 0.513, 0.513, 0.514, 0.514, 0.515, 0.515, 0.515, 0.516, 0.516, 0.517, 0.517, 0.517, 0.518, 0.518, 0.518, 0.519, 0.519, 0.519, 0.519, 0.520, 0.520, 0.520, 0.520, 0.521, 0.521, 0.521, 0.521, 0.522, 0.522, 0.522, 0.522, 0.523, 0.523, 0.523, 0.523, 0.523, 0.524, 0.524, 0.524, 0.524, 0.524, 0.524, 0.525, 0.525, 0.525, 0.525, 0.525, 0.525, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.527, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.528, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.529, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.530, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.531, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.532, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.534, 0.534, 0.534, 0.534, 0.534, 0.534};


  // Data PbPb 2.76 ATeV 10h.pass2 (Pt-Range: 0.15 - 2) (up to N_acc = 1000)
  Double_t sNbinsData276[1001] = {0.000, 0.487, 0.488, 0.491, 0.497, 0.504, 0.510, 0.515, 0.522, 0.527, 0.531, 0.536, 0.539, 0.543, 0.546, 0.550, 0.551, 0.554, 0.556, 0.558, 0.560, 0.561, 0.563, 0.564, 0.566, 0.567, 0.568, 0.568, 0.570, 0.571, 0.571, 0.572, 0.573, 0.574, 0.574, 0.575, 0.575, 0.576, 0.577, 0.577, 0.578, 0.578, 0.578, 0.580, 0.580, 0.580, 0.580, 0.581, 0.581, 0.582, 0.582, 0.583, 0.583, 0.583, 0.583, 0.584, 0.584, 0.585, 0.585, 0.586, 0.586, 0.586, 0.586, 0.586, 0.586, 0.587, 0.587, 0.587, 0.588, 0.588, 0.588, 0.589, 0.589, 0.589, 0.589, 0.589, 0.590, 0.590, 0.590, 0.591, 0.591, 0.590, 0.591, 0.592, 0.592, 0.592, 0.592, 0.592, 0.593, 0.593, 0.593, 0.593, 0.593, 0.594, 0.594, 0.594, 0.594, 0.595, 0.594, 0.594, 0.595, 0.596, 0.596, 0.596, 0.596, 0.596, 0.597, 0.597, 0.596, 0.597, 0.597, 0.598, 0.597, 0.598, 0.598, 0.598, 0.598, 0.599, 0.598, 0.599, 0.598, 0.599, 0.599, 0.600, 0.600, 0.600, 0.600, 0.600, 0.600, 0.601, 0.601, 0.601, 0.601, 0.601, 0.602, 0.602, 0.602, 0.602, 0.602, 0.602, 0.602, 0.603, 0.603, 0.602, 0.603, 0.603, 0.603, 0.604, 0.604, 0.604, 0.604, 0.605, 0.604, 0.604, 0.605, 0.605, 0.604, 0.605, 0.605, 0.605, 0.605, 0.605, 0.606, 0.606, 0.606, 0.606, 0.606, 0.606, 0.607, 0.607, 0.606, 0.607, 0.607, 0.608, 0.608, 0.607, 0.607, 0.608, 0.608, 0.608, 0.608, 0.608, 0.609, 0.608, 0.609, 0.608, 0.610, 0.609, 0.609, 0.609, 0.609, 0.610, 0.609, 0.610, 0.610, 0.610, 0.610, 0.610, 0.610, 0.610, 0.611, 0.611, 0.611, 0.611, 0.610, 0.611, 0.611, 0.611, 0.612, 0.612, 0.612, 0.612, 0.612, 0.612, 0.613, 0.613, 0.613, 0.613, 0.613, 0.612, 0.613, 0.612, 0.613, 0.613, 0.613, 0.613, 0.614, 0.613, 0.614, 0.613, 0.613, 0.614, 0.614, 0.614, 0.614, 0.614, 0.615, 0.615, 0.615, 0.615, 0.615, 0.615, 0.615, 0.615, 0.615, 0.616, 0.615, 0.616, 0.616, 0.616, 0.616, 0.615, 0.616, 0.616, 0.616, 0.617, 0.617, 0.617, 0.617, 0.616, 0.617, 0.617, 0.617, 0.618, 0.618, 0.618, 0.617, 0.618, 0.617, 0.617, 0.618, 0.618, 0.618, 0.618, 0.619, 0.618, 0.618, 0.619, 0.619, 0.618, 0.619, 0.619, 0.619, 0.619, 0.619, 0.619, 0.619, 0.619, 0.620, 0.620, 0.619, 0.620, 0.620, 0.620, 0.620, 0.621, 0.621, 0.620, 0.621, 0.621, 0.620, 0.621, 0.620, 0.620, 0.621, 0.621, 0.621, 0.621, 0.621, 0.621, 0.621, 0.621, 0.622, 0.621, 0.622, 0.621, 0.622, 0.622, 0.621, 0.621, 0.622, 0.622, 0.622, 0.622, 0.622, 0.622, 0.622, 0.623, 0.623, 0.622, 0.623, 0.622, 0.623, 0.623, 0.623, 0.623, 0.624, 0.623, 0.623, 0.623, 0.624, 0.623, 0.624, 0.624, 0.624, 0.623, 0.624, 0.624, 0.623, 0.624, 0.624, 0.624, 0.624, 0.624, 0.625, 0.625, 0.624, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.626, 0.625, 0.625, 0.625, 0.626, 0.626, 0.626, 0.625, 0.626, 0.626, 0.626, 0.626, 0.626, 0.626, 0.627, 0.626, 0.626, 0.627, 0.627, 0.627, 0.627, 0.627, 0.626, 0.627, 0.627, 0.627, 0.627, 0.627, 0.627, 0.627, 0.627, 0.627, 0.627, 0.628, 0.627, 0.628, 0.628, 0.628, 0.628, 0.627, 0.628, 0.628, 0.628, 0.628, 0.628, 0.628, 0.628, 0.628, 0.628, 0.628, 0.628, 0.628, 0.629, 0.628, 0.628, 0.629, 0.629, 0.630, 0.629, 0.629, 0.629, 0.629, 0.629, 0.630, 0.629, 0.630, 0.629, 0.630, 0.630, 0.629, 0.630, 0.630, 0.629, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.631, 0.630, 0.631, 0.630, 0.630, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.631, 0.632, 0.631, 0.632, 0.631, 0.632, 0.632, 0.631, 0.631, 0.632, 0.631, 0.632, 0.632, 0.632, 0.633, 0.632, 0.632, 0.632, 0.632, 0.633, 0.632, 0.632, 0.632, 0.632, 0.632, 0.633, 0.633, 0.633, 0.633, 0.633, 0.633, 0.633, 0.633, 0.633, 0.632, 0.633, 0.633, 0.633, 0.633, 0.633, 0.633, 0.633, 0.633, 0.634, 0.633, 0.634, 0.633, 0.634, 0.634, 0.634, 0.634, 0.634, 0.633, 0.633, 0.634, 0.634, 0.634, 0.633, 0.634, 0.633, 0.634, 0.634, 0.634, 0.634, 0.634, 0.634, 0.634, 0.634, 0.634, 0.634, 0.635, 0.635, 0.634, 0.634, 0.635, 0.635, 0.634, 0.634, 0.635, 0.634, 0.635, 0.635, 0.634, 0.635, 0.635, 0.635, 0.635, 0.635, 0.635, 0.635, 0.635, 0.635, 0.636, 0.636, 0.635, 0.635, 0.635, 0.636, 0.636, 0.635, 0.636, 0.636, 0.635, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.637, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.637, 0.637, 0.637, 0.636, 0.637, 0.636, 0.636, 0.636, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.637, 0.638, 0.638, 0.637, 0.638, 0.637, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.638, 0.639, 0.638, 0.638, 0.638, 0.639, 0.638, 0.639, 0.638, 0.639, 0.638, 0.638, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.638, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.639, 0.640, 0.639, 0.640, 0.640, 0.640, 0.639, 0.639, 0.639, 0.639, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.640, 0.641, 0.640, 0.640, 0.641, 0.641, 0.641, 0.641, 0.641, 0.641, 0.641, 0.641, 0.641, 0.640, 0.640, 0.640, 0.641, 0.641, 0.641, 0.641, 0.641, 0.641, 0.642, 0.641, 0.642, 0.642, 0.641, 0.641, 0.641, 0.642, 0.641, 0.642, 0.641, 0.641, 0.641, 0.642, 0.641, 0.642, 0.642, 0.641, 0.642, 0.642, 0.641, 0.642, 0.641, 0.641, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.642, 0.643, 0.642, 0.642, 0.642, 0.642, 0.642, 0.643, 0.642, 0.642, 0.643, 0.643, 0.643, 0.642, 0.642, 0.643, 0.642, 0.643, 0.643, 0.643, 0.642, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.642, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.644, 0.643, 0.643, 0.643, 0.643, 0.643, 0.643, 0.644, 0.643, 0.643, 0.643, 0.643, 0.644, 0.643, 0.644, 0.644, 0.643, 0.643, 0.643, 0.644, 0.643, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.643, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.645, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.644, 0.645, 0.644, 0.645, 0.644, 0.644, 0.645, 0.644, 0.644, 0.645, 0.645, 0.644, 0.644, 0.644, 0.644, 0.645, 0.645, 0.645, 0.645, 0.644, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.646, 0.646, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.645, 0.646, 0.645, 0.646, 0.646, 0.645, 0.646, 0.646, 0.646, 0.645, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646, 0.645, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646, 0.645, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646};

// -- End mean pT values for single multiplicity bins --


// -- Selection of MC/Data --

if (fMC) { // - MC -

  if (fMCAMPT) {

//     Printf(" -- MC, 2.76 ATeV - AMPT -- ");

    nbins = nbinsMC276A;
    centbins = centbinsMC276A;
    sNbinsMC276 = sNbinsMC276A;

    maxNrTracksMC = 300;
    minNrBinMC = 4;

  }
  else {

//     Printf(" -- MC, 2.76 ATeV - HIIJNG -- ");

    nbins = nbinsMC276H;
    centbins = centbinsMC276H;
    sNbinsMC276 = sNbinsMC276H;

    maxNrTracksMC = 500;
    minNrBinMC = 6;

  }

} // - End MC -
else { // - Data -

//   Printf(" -- Data, 2.76 ATeV -- ");

  nbins = nbinsData276;
  centbins = centbinsData276;

} // - End data -

// -- End selection of MC/Data --

// --- End mean pT values ---


// --- Ratio of pT distributions MC ESD / MC truth ---

// MC PbPb 2.76 ATeV 11a10a_bis (Hijing) (Pt-Range: 0.15 - 2, array-range: 0 - 2.1 GeV/c)
Double_t nbinsPtRatio[84] = {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.667, 0.850, 0.871, 0.880, 0.884, 0.890, 0.899, 0.913, 0.923, 0.931, 0.935, 0.937, 0.941, 0.944, 0.944, 0.945, 0.946, 0.947, 0.947, 0.948, 0.948, 0.949, 0.949, 0.950, 0.949, 0.949, 0.948, 0.950, 0.949, 0.949, 0.951, 0.953, 0.951, 0.951, 0.950, 0.950, 0.950, 0.952, 0.954, 0.953, 0.955, 0.954, 0.955, 0.955, 0.955, 0.957, 0.955, 0.955, 0.956, 0.956, 0.952, 0.951, 0.950, 0.948, 0.950, 0.951, 0.949, 0.948, 0.945, 0.946, 0.943, 0.943, 0.942, 0.940, 0.941, 0.939, 0.936, 0.936, 0.936, 0.932, 0.932, 0.927, 0.928, 0.927, 0.000, 0.000, 0.000, 0.000};

// --- End ratio of pT distributions MC ESD / MC truth ---


  if (!fMC || fMCType < 3) {
    nrESDTracks = fESD->GetNumberOfTracks();
  }

  if (fMC) {
    nrMCTracks = stack->GetNtrack();
  }

//   Printf("\n\n nrESDTracks: %i   nrMCTracks : %i \n",nrESDTracks,nrMCTracks);


  // Get event centrality
  if (!fMC || fMCType < 3) {

   if(fUseCentrality != 0) {

    AliCentrality *esdCentrality = fESD->GetCentrality();

    //V0
    if (fUseCentrality == 1) {

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

//       if      ( centralityVZERO >   0. && centralityVZERO <  10.) centralityVZEROBin =  0;
//       else if ( centralityVZERO >= 10. && centralityVZERO <  20.) centralityVZEROBin = 10;
//       else if ( centralityVZERO >= 20. && centralityVZERO <  30.) centralityVZEROBin = 20;
//       else if ( centralityVZERO >= 30. && centralityVZERO <  40.) centralityVZEROBin = 30;
//       else if ( centralityVZERO >= 40. && centralityVZERO <  50.) centralityVZEROBin = 40;
//       else if ( centralityVZERO >= 50. && centralityVZERO <  60.) centralityVZEROBin = 50;
//       else if ( centralityVZERO >= 60. && centralityVZERO <  70.) centralityVZEROBin = 60;
//       else if ( centralityVZERO >= 70. && centralityVZERO <  80.) centralityVZEROBin = 70;
//       else if ( centralityVZERO >= 80. && centralityVZERO <  90.) centralityVZEROBin = 80;
//       else if ( centralityVZERO >= 90. && centralityVZERO <  99.) centralityVZEROBin = 90;
//       else if ( centralityVZERO >= 99. ) centralityVZEROBin = 100;
//       else if ( centralityVZERO <= 0.  ) centralityVZEROBin = 100;
// 
//       centBin = centralityVZEROBin / 10;
//       evMptCent = centbins[centBin];
    }

//     Printf("\nRaw mult: %i   CentVZERO: %f   CentVZERObin: %i   centBin: %i   evMptCent: %.3f ",nrESDTracks,centralityVZERO,centralityVZEROBin,centBin,evMptCent);

   }
  } // End get event centrality


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

    // Multiplicity distributions
    fMult->Fill(nrTracks);
    fMultNbins->Fill(nrTracks);
    fMultSum->Fill(nrTracks,nrTracks);
    fMultSumPt->Fill(nrTracks,sumPt);
    // Same for first bin divided in 4
    if (nrTracks < 101) {
      fMult1->Fill(nrTracks);
      fMultSum1->Fill(nrTracks,nrTracks);
      fMultSumPt1->Fill(nrTracks,sumPt);
    }
    // Same for five bins in 0 < Nacc < 50
    if (nrTracks < 51) {
      fMult10->Fill(nrTracks);
      fMultSum10->Fill(nrTracks,nrTracks);
      fMultSumPt10->Fill(nrTracks,sumPt);
    }
    // Same for events with centrality < 80%
    if (centralityVZEROBin < 80) {
	fMult80->Fill(nrTracks);
	fMultSum80->Fill(nrTracks,nrTracks);
	fMultSumPt80->Fill(nrTracks,sumPt);
      // Same for first bin divided in 4 only for events with centrality < 80%
      if (nrTracks < 101) {
	fMult801->Fill(nrTracks);
	fMultSum801->Fill(nrTracks,nrTracks);
	fMultSumPt801->Fill(nrTracks,sumPt);
      }
      // Same for five bins in 0 < Nacc < 50 only for events with centrality < 80%
      if (nrTracks < 51) {
	fMult810->Fill(nrTracks);
	fMultSum810->Fill(nrTracks,nrTracks);
	fMultSumPt810->Fill(nrTracks,sumPt);
      }
    }

    // Number of pairs in event
    nrPairs = 0.5 * nrTracks * (nrTracks-1);
    fMultNrPairs->Fill(nrTracks,nrPairs);
    // Same for first bin divided in 4
    if (nrTracks < 101) {
      fMultNrPairs1->Fill(nrTracks,nrPairs);
    }
    // Same for five bins in 0 < Nacc < 50
    if (nrTracks < 51) {
      fMultNrPairs10->Fill(nrTracks,nrPairs);
    }
    // Same for events with centrality < 80%
    if (centralityVZEROBin < 80) {
	fMultNrPairs80->Fill(nrTracks,nrPairs);
      // Same for first bin divided in 4 only for events with centrality < 80%
      if (nrTracks < 101) {
	fMultNrPairs801->Fill(nrTracks,nrPairs);
      }
      // Same for five bins in 0 < Nacc < 50 only for events with centrality < 80%
      if (nrTracks < 51) {
	fMultNrPairs810->Fill(nrTracks,nrPairs);
      }
    }

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
    fMultEventMeanPtNbins->Fill(nrTracks,eventMeanPt);
    fMultEventMeanPtSqNbins->Fill(nrTracks,eventMeanPtSq);

//     Printf("nrTracks: %i   sumPt: %.8f   meanPt: %.8f   meanPtSq: %.8f\n",nrTracks,sumPt,eventMeanPt,eventMeanPtSq);


    // Centrality V0
    if (!fMC || fMCType < 3) {

     if (fUseCentrality == 1) {

      // Centrality distributions
      fCent->Fill(centralityVZEROBin);
      fCentSum->Fill(centralityVZEROBin,nrTracks);
      fCentSumPt->Fill(centralityVZEROBin,sumPt);

      // Number of pairs in event
      fCentNrPairs->Fill(centralityVZEROBin,nrPairs);

      // Mean-Pt and Mean-Pt squared depending on centrality
      fCentEventMeanPt->Fill(centralityVZEROBin,eventMeanPt);
      fCentEventMeanPtSq->Fill(centralityVZEROBin,eventMeanPtSq);

      if (centBin == 0) {
	fEventMeanPtCent05->Fill(eventMeanPt);
// 	Printf("centralityV0 = %.3f   --   0-5% -- \n",centralityVZERO);
      }
      else if (centBin == 4 || centBin == 5) {
	fEventMeanPtCent2030->Fill(eventMeanPt);
// 	Printf("centralityV0 = %.3f   -- 20-30% -- \n",centralityVZERO);
      }
      else if (centBin == 14 || centBin == 15) {
	fEventMeanPtCent7080->Fill(eventMeanPt);
// 	Printf("centralityV0 = %.3f   -- 70-80% -- \n",centralityVZERO);
      }

//       Printf("CentV0bin: %i   NrTracks: %i   NrPairs: %.0f   SumPt: %.1f\n",centralityVZEROBin,nrTracks,nrPairs,sumPt);
     }
    }


    // --- Setting of mean pT values

    if (fMC) { // - MC -

//       if (nrTracks <= 500) {
//       if (nrTracks <= 200) {
      if (nrTracks <= maxNrTracksMC) {
	evMptMult = sNbinsMC276[nrTracks];
	evMptCent = evMptMult;
      }
      else {
// 	for (int k=6; k<32; k++) {
// 	for (int k=3; k<32; k++) {
	for (int k=minNrBinMC; k<32; k++) {
	  if (nrTracks > 100 * (k-1) && nrTracks <= 100 * k) {
	    evMptMult = nbins[k];
	    break;
	  }
	}
      }

    } // - End MC -
    else { // - Data -

      if (nrTracks <= 1000) {
	evMptMult = sNbinsData276[nrTracks];
	evMptCent = evMptMult;
      }
      else {
	for (int k=11; k<32; k++) {
	  if (nrTracks > 100 * (k-1) && nrTracks <= 100 * k) {
	    evMptMult = nbins[k];
	    break;
	  }
	}
      }

    } // - End data -


//     Printf("maxNrTracksMC = %3i   minNrBinMC = %3i   ",maxNrTracksMC,minNrBinMC);
//     Printf("nrTracks = %3i   evMptMult = %.3f   evMptCent = %.3f \n",nrTracks,evMptMult,evMptCent);

    // --- End setting of mean pT values


    // --- Two-particle correlator ---

    // -- Analysis in multiplicity bins --

    for (int i=0; i<nrTracks; i++) {
      for (int j=i+1; j<nrTracks; j++) {
	twoPartCorrPair = (tracks[i] - evMptMult) * (tracks[j] - evMptMult);
	twoPartCorrEv = twoPartCorrEv + twoPartCorrPair;
      }
    }

    twoPartCorrEvSq = twoPartCorrEv * twoPartCorrEv;
    fTwoPartCorrEv->Fill(nrTracks,twoPartCorrEv);
    fTwoPartCorrEvSq->Fill(nrTracks,twoPartCorrEvSq);
    // Same for first bin divided in 4
    if (nrTracks < 101) {
      fTwoPartCorrEv1->Fill(nrTracks,twoPartCorrEv);
      fTwoPartCorrEvSq1->Fill(nrTracks,twoPartCorrEvSq);
    }
    // Same for five bins in 0 < Nacc < 50
    if (nrTracks < 51) {
      fTwoPartCorrEv10->Fill(nrTracks,twoPartCorrEv);
      fTwoPartCorrEvSq10->Fill(nrTracks,twoPartCorrEvSq);
    }
    // Same for events with centrality < 80%
    if (centralityVZEROBin < 80) {
	fTwoPartCorrEv80->Fill(nrTracks,twoPartCorrEv);
	fTwoPartCorrEvSq80->Fill(nrTracks,twoPartCorrEvSq);
      // Same for first bin divided in 4 only for events with centrality < 80%
      if (nrTracks < 101) {
	fTwoPartCorrEv801->Fill(nrTracks,twoPartCorrEv);
	fTwoPartCorrEvSq801->Fill(nrTracks,twoPartCorrEvSq);
      }
      // Same for five bins in 0 < Nacc < 50 only for events with centrality < 80%
      if (nrTracks < 51) {
	fTwoPartCorrEv810->Fill(nrTracks,twoPartCorrEv);
	fTwoPartCorrEvSq810->Fill(nrTracks,twoPartCorrEvSq);
      }
    }

    // -- End analysis in multiplicity bins --


    // -- Analysis in centrality bins --
    if (!fMC || fMCType < 3) {

     // Centrality V0
     if (fUseCentrality == 1) {

      if (centralityVZEROBin < 91) {

	for (int i=0; i<nrTracks; i++) {
	  for (int j=i+1; j<nrTracks; j++) {
	    twoPartCorrPairCent = (tracks[i] - evMptCent) * (tracks[j] - evMptCent);
	    twoPartCorrEvCent = twoPartCorrEvCent + twoPartCorrPairCent;
	  }
	}

	twoPartCorrEvCentSq = twoPartCorrEvCent * twoPartCorrEvCent;
	fTwoPartCorrEvCent->Fill(centralityVZEROBin,twoPartCorrEvCent);
	fTwoPartCorrEvCentSq->Fill(centralityVZEROBin,twoPartCorrEvCentSq);

      }
     }
    } // -- End analysis in centrality bins --

    // --- End two-particle correlator ---


  } // --- End calculation of various values and filling of histograms ---


  // Post output data
  PostData(1, fOutputList);

  // Delet pointer
  if (fRandom3) delete fRandom3; fRandom3 = 0;

}

void AliAnalysisTaskPtFlucPbPb::Terminate(Option_t *)
{
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    Printf("ERROR: fOutputList not available\n");
    return;
  }

}
