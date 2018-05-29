#include "AliAnalysisTaskJetLikeCorrelation.h"

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "TList.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "TChain.h"
#include "THnSparse.h"
#include "TGrid.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

#include "AliAnalysisManager.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliCollisionGeometry.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliEventplane.h"
#include "AliMultSelection.h"
#include "AliEventPoolManager.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "TString.h"



using namespace std;

ClassImp(AliExtractedTrack)
ClassImp(AliExtractedEvent)
ClassImp(MixedParticle)
ClassImp(AliAnalysisTaskJetLikeCorrelation)

AliExtractedTrack::AliExtractedTrack()
  :
    fPt(0),
    fEta(0),
    fPhi(0),
    fCharge(0), 
    fFilterBit(0),
    kIsPrimary(0) {

}

AliExtractedTrack::AliExtractedTrack(const AliExtractedTrack &track) 
 : fPt(track.fPt),
   fEta(track.fEta),
   fPhi(track.fPhi),
   fCharge(track.fCharge),
   fFilterBit(track.fFilterBit),
   kIsPrimary(track.kIsPrimary) {
}

AliExtractedTrack &AliExtractedTrack::operator=(const AliExtractedTrack &track) {
  if (this != &track) {
    this->fPt = track.fPt;
    this->fEta = track.fEta;
    this->fPhi = track.fPhi;
    this->fCharge = track.fCharge;
    this->fFilterBit = track.fFilterBit;
    this->kIsPrimary = track.kIsPrimary;
  }
  return *this;
}

AliExtractedTrack::AliExtractedTrack(AliExtractedTrack &&track)
  : fPt(std::move(track.fPt)),
  fEta(std::move(track.fEta)), 
  fPhi(std::move(track.fPhi)), 
  fCharge(std::move(track.fCharge)),
  fFilterBit(std::move(track.fFilterBit)),
  kIsPrimary(std::move(track.kIsPrimary))
{
  track.fPt = 0;
  track.fEta = 0;
  track.fPhi = 0;
  track.fCharge = 0;
  track.fFilterBit = 0;
  track.kIsPrimary = 0;
}

AliExtractedTrack &AliExtractedTrack::operator=(AliExtractedTrack &&track)
{
  if (this != &track) {
    fPt = std::move(track.fPt);
    fEta = std::move(track.fEta); 
    fPhi = std::move(track.fPhi); 
    fCharge = std::move(track.fCharge); 
    fFilterBit = std::move(track.fFilterBit);
    kIsPrimary = std::move(track.kIsPrimary);
  }
  return *this;

}

AliExtractedTrack::~AliExtractedTrack() {
}


AliExtractedEvent::AliExtractedEvent()
  :
    fRunNumber(0), fEventID(0), fBSign(0), fCentrality(0), fZVertex(0), 
    fEventPlane(0), fEventPlaneV0A(0), fEventPlaneV0C(0), fEventPlaneTPC(0)
{

}

AliExtractedEvent::AliExtractedEvent(const AliExtractedEvent &event)  :
  fRunNumber(event.fRunNumber), 
  fEventID(event.fEventID), 
  fBSign(event.fBSign),
  fCentrality(event.fCentrality), 
  fZVertex(event.fZVertex), 
  fTracks(event.fTracks),
  fEventPlane(event.fEventPlane),
  fEventPlaneV0A(event.fEventPlaneV0A),
  fEventPlaneV0C(event.fEventPlaneV0C),
  fEventPlaneTPC(event.fEventPlaneTPC)
{
}

AliExtractedEvent &AliExtractedEvent::operator=(const AliExtractedEvent &event) 
{
  if (this != &event) {
    this->fRunNumber = event.fRunNumber;
    this->fEventID = event.fEventID; 
    this->fCentrality = event.fCentrality; 
    this->fZVertex = event.fZVertex; 
    this->fTracks = event.fTracks;
    this->fEventPlane = event.fEventPlane;
    this->fEventPlaneV0A = event.fEventPlaneV0A;
    this->fEventPlaneV0C = event.fEventPlaneV0C;
    this->fEventPlaneTPC = event.fEventPlaneTPC;
  }
  return *this;
}

AliExtractedEvent::AliExtractedEvent(AliExtractedEvent &&event)
  : fRunNumber(std::move(event.fRunNumber)),
  fEventID(std::move(event.fEventID)), 
  fBSign(std::move(event.fBSign)),
  fCentrality(std::move(event.fCentrality)), 
  fZVertex(std::move(event.fZVertex)), 
  fTracks(std::move(event.fTracks)),
  fEventPlane(std::move(event.fEventPlane)),
  fEventPlaneV0A(std::move(event.fEventPlaneV0A)),
  fEventPlaneV0C(std::move(event.fEventPlaneV0C)),
  fEventPlaneTPC(std::move(event.fEventPlaneTPC))
{
}

AliExtractedEvent &AliExtractedEvent::operator=(AliExtractedEvent &&event)
{
  if (this != &event) {
   this->fRunNumber = std::move(event.fRunNumber);
   this->fEventID = std::move(event.fEventID); 
   this->fBSign = std::move(event.fBSign); 
   this->fCentrality = std::move(event.fCentrality); 
   this->fZVertex = std::move(event.fZVertex); 
   this->fTracks = std::move(event.fTracks);
   this->fEventPlane = std::move(event.fEventPlane);
   this->fEventPlaneV0A = std::move(event.fEventPlaneV0A);
   this->fEventPlaneV0C = std::move(event.fEventPlaneV0C);
   this->fEventPlaneTPC = std::move(event.fEventPlaneTPC);
  }
  return *this;
}

AliExtractedEvent::~AliExtractedEvent() {
}

AliAnalysisTaskJetLikeCorrelation::AliAnalysisTaskJetLikeCorrelation() :
  AliAnalysisTaskSE(),
  fInputHandler(0), fMCHandler(0), fMCEvent(0), fMCHeader(0), 
  fHistZVertex(0), fHistPt(0), fHistPhi(0),fHistCent(0),
  fHistEta(0), fHistPos(0), fHistNeg(0), fHistEventPlaneV0A(0), fHistEventPlaneV0C(0),
  fHistEventPlaneTPC(0), fPoolMgr(0), fInputRoot(0), fInputTree(0),
  fMCCorrection(0), fMCTruth(0), bUseMixingPool(0), fMixingPoolSize(100), fDebug(0), fMinNumTrack(500),
  fEtaCut(0.9), fPhiCut(TMath::TwoPi()), fMinPtTrigCut(3.0), fTwoTrackEffCut(0.02), 
  fFilterBit(128), fTrackDepth(5000), fEventID(0), fIsFirstEvent(1),
  fMinimumPtABinMerging(5), fResonancesVCut(0.02), fConversionsVCut(0.04),
  fHistNevtSame(0),
  fHistEtaSparse(0), fHistV2(0),
  fHistContamination(0), fHighPtMixing(0), fTreeStart(0), fTreeEnd(0),
  fCentPercentile(0), fZVertex(0), fEventPlane(0), fEventPlaneV0A(0),
  fEventPlaneV0C(0), fEventPlaneTPC(0), pi(TMath::Pi()), fCollision(kPbPb), fPeriod(k10h),
  flowQnVectorTask(0), fFlowQnVectorMgr(0)
{

}

AliAnalysisTaskJetLikeCorrelation::AliAnalysisTaskJetLikeCorrelation(const char *name) :
  AliAnalysisTaskSE(name),
  fInputHandler(0), fMCHandler(0), fMCEvent(0), fMCHeader(0), 
  fHistZVertex(0), fHistPt(0), fHistPhi(0),fHistCent(0),
  fHistEta(0), fHistPos(0), fHistNeg(0), fHistEventPlaneV0A(0), fHistEventPlaneV0C(0),
  fHistEventPlaneTPC(0), fPoolMgr(0), fInputRoot(0), fInputTree(0),
  fMCCorrection(0), fMCTruth(0), bUseMixingPool(0), fMixingPoolSize(100), fDebug(0), fMinNumTrack(500),
  fEtaCut(0.9), fPhiCut(TMath::TwoPi()), fMinPtTrigCut(3.0), fTwoTrackEffCut(0.02), 
  fFilterBit(128), fTrackDepth(5000), fEventID(0), fIsFirstEvent(1),
  fMinimumPtABinMerging(5), fResonancesVCut(0.02), fConversionsVCut(0.04),
  fHistNevtSame(0), 
  fHistEtaSparse(0), fHistV2(0),
  fHistContamination(0), fHighPtMixing(0), fTreeStart(0), fTreeEnd(0),
  fCentPercentile(0), fZVertex(0), fEventPlane(0), fEventPlaneV0A(0),
  fEventPlaneV0C(0), fEventPlaneTPC(0), pi(TMath::Pi()), fCollision(kPbPb), fPeriod(k10h),
  flowQnVectorTask(0), fFlowQnVectorMgr(0)
{
  DefineInput(0, TChain::Class());
//  DefineInput(1, TList::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
  DefineOutput(7, TList::Class());
//  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
//    DefineOutput(iplane+1, TList::Class());    //Incl, In, Out, M1, M2, M3, M4 : total 7
//    cout << "DefineOutput " << iplane << endl;
//  }
}

AliAnalysisTaskJetLikeCorrelation::~AliAnalysisTaskJetLikeCorrelation() {

  for (int iin = 0; iin < fNumberOfPlanes; iin++) {
    if (fOutput[iin] && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput[iin];
    }
  }
}

void AliAnalysisTaskJetLikeCorrelation::Terminate(Option_t *option) {

  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    fOutput[iplane] = dynamic_cast<TList*>(GetOutputData(iplane+1));
    if (NULL == fOutput[iplane]) {
      AliFatal(Form("fOutput%d == NULL", iplane));
      return;
    }
  }
  
}


void AliAnalysisTaskJetLikeCorrelation::UserCreateOutputObjects() {
//  fOutput = new TList();
//  fOutput->SetOwner(kTRUE);
  fIsFirstEvent = 1;

  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    fOutput[iplane] = new TList();
    fOutput[iplane]->SetOwner(kTRUE);
  }

  fHistZVertex = new TH1D("fHistZVertex", "fHistZVertex", 400, -20, 20);
  fHistPt = new TH1D("fHistPt", "fHistPt", 200, 0, 20);
  fHistPhi = new TH1D("fHistPhi", "fHistPhi", 315*2 , 0, 3.15*2);
  fHistCent = new TH1D("fHistCent", "fHistCent", 100, 0, 100);
  fHistEta = new TH1D("fHistEta", "fHistEta", 200, -1, 1);
  fHistPos = new TH1D("fHistPos", "fHistPos", 200, 0, 20);
  fHistNeg = new TH1D("fHistNeg", "fHistNeg", 200, 0, 20);
  fHistEventPlaneV0A = new TH1D("fHistEventPlaneV0A", "fHistEventPlaneV0A", 315 , 0, 3.15);
  fHistEventPlaneV0C = new TH1D("fHistEventPlaneV0C", "fHistEventPlaneV0C", 315 , 0, 3.15);
  fHistEventPlaneTPC = new TH1D("fHistEventPlaneTPC", "fHistEventPlaneTPC", 315 , 0, 3.15);

  fHistV2 = new TProfile("fHistV2", "fHistV2", 10, 0, 100);
  fHistResolutionV2[0] = new TProfile("fHistResolutionV2_01", "fHistResolutionV2_01", 10, 0, 100);
  fHistResolutionV2[1] = new TProfile("fHistResolutionV2_02", "fHistResolutionV2_02", 10, 0, 100);
  fHistResolutionV2[2] = new TProfile("fHistResolutionV2_12", "fHistResolutionV2_12", 10, 0, 100);

  int nbins[6] = {2, 5, 20, 9, 20, 250};
  // Contamination/Primary, In/out(5), Cent, ZVtx, Eta, Pt // total 6
  double maxbin[6] = {1.5, 4.5, 100, 9,1,25 };
  double minbin[6] = {-0.5, -0.5, 0, -9, -1, 0};
  if (fMCCorrection) {
    fHistContamination = new THnSparseD("fHistContamination", "fHistContamination", 6, nbins, minbin, maxbin);
    fHistContamination->Sumw2();
  }
  InitHistograms();

  
  fOutput[0]->Add(fHistZVertex);
  fOutput[0]->Add(fHistPt);
  fOutput[0]->Add(fHistPhi);
  fOutput[0]->Add(fHistCent);
  fOutput[0]->Add(fHistEta);
  fOutput[0]->Add(fHistPos);
  fOutput[0]->Add(fHistNeg);
  fOutput[0]->Add(fHistEtaSparse);
  fOutput[0]->Add(fHistEventPlaneV0A);
  fOutput[0]->Add(fHistEventPlaneV0C);
  fOutput[0]->Add(fHistEventPlaneTPC);
  fOutput[0]->Add(fHistV2);
  fOutput[0]->Add(fHistResolutionV2[0]);
  fOutput[0]->Add(fHistResolutionV2[1]);
  fOutput[0]->Add(fHistResolutionV2[2]);
  fOutput[0]->Add(fHistNevtSame);
  if (fMCCorrection)
    fOutput[0]->Add(fHistContamination);

  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    fOutput[iplane]->Add(fHistPtSame[iplane]);
  }

//  int lMinPtBin = GetPtBin(fMinPtTrigCut);
  

  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    for (int icent = 0; icent < fCentArray.GetSize() - 1; icent++) {
      for (int izvertex = 0; izvertex < fZVertexArray.GetSize() - 1; izvertex++) {
        for (int iptt = 0; iptt < fPttArray.GetSize() - 1; iptt++) {
          for (int ipta = 0; ipta < fPtaArray.GetSize() - 1; ipta++) {
            if (fPttArray[iptt] + 0.1 < fPtaArray[ipta] ) continue;

            fOutput[iplane]->Add(fHistdEtadPhiSame[iplane][icent][izvertex][iptt][ipta]);
            fOutput[iplane]->Add(fHistdEtadPhiMixed[iplane][icent][izvertex][iptt][ipta]);
            if (fMCCorrection) {
              fOutput[iplane]->Add(fHistdEtadPhiSameMCCorrPrim[iplane][icent][izvertex][iptt][ipta]);
              fOutput[iplane]->Add(fHistdEtadPhiSameMCCorrCont[iplane][icent][izvertex][iptt][ipta]);
            }
          }
        }
      }
    }
  }
    
  if (fCollision == kPbPb) {
    flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections*>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
    //    cout << flowQnVectorTask << endl;
    //    AliInfo("FlowQnVectorTask");
    if (flowQnVectorTask != NULL) {
      /* AliQnCorrectionsManager *fFlowQnVectorMgr; shall be a member of the user's analysis task */
      /* And store the framework manager */
      fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
    } else {
      AliFatal("Flow Qn vector corrections framework needed but it's not present. ABORTING!");
    }
  }

  SetupMixing();
 
//  ConnectInputHandler();

  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    PostData(iplane+1, fOutput[iplane]);
  }


}

void AliAnalysisTaskJetLikeCorrelation::ConnectInputHandler() {
  fInputHandler = (AliInputEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  fMCHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}

void AliAnalysisTaskJetLikeCorrelation::UserExec(Option_t *option) {

  // Main loop executed for each event
  // Create Pointer to reconstructed event
  AliVEvent *event = InputEvent();
  if (!event) { 
    AliError("ERROR: No Event! End Event");
    return;
  }

  AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*> (event);
  if (!aodEvent) {
    AliError("ERROR: No AODEvent Available");
    return;
  }

  AliAODHeader *aodHeader = static_cast<AliAODHeader*> (aodEvent->GetHeader());
  if (!aodHeader) {
    AliError("ERROR: No Header Available");
    return;
  }

  fEventID = GetEventIdAsLong(aodHeader);
  int runnumber = aodHeader->GetRunNumber();

  if (fIsFirstEvent && bUseMixingPool) {
    TGrid::Connect("alien:");
    string  filename;
    if (fCollision == kPbPb)
      filename = Form("alien:///alice/cern.ch/user/h/hyeonjoo/MixingPool_10h/output/000%d/AnalysisResults.root", runnumber); // pbpb
    else 
      filename = Form("alien:///alice/cern.ch/user/h/hyeonjoo/MixingPool_11a/output/000%d/AnalysisResults.root", runnumber); // pbpb
    //    string filename = Form("alien:///alice/cern.ch/user/h/hyeonjoo/Hybrid_11aMixing/output/000%d/AnalysisResults.root", runnumber); // pp
    //    string filename = Form("./146858.root");
    fInputRoot = TFile::Open(filename.c_str());
    if (!fInputRoot) {
      cout << "File : " << filename << " doesn't exist!" << endl;
      return;
    }
    fInputTree = (TTree*)fInputRoot->Get("cTree");
    fInputTree->SetMaxVirtualSize(1024*1024*1024);
    fInputTree->SetBranchAddress("fEvent", &fHighPtMixing);
    if (fInputTree->GetEntries() > fMixingPoolSize) {       // for pp, 1M, PbPb : 0.1M
      TRandom r;
      int temp_value = int(r.Uniform(0,fInputTree->GetEntries() - fMixingPoolSize));
      fTreeStart = temp_value;
      fTreeEnd = temp_value + fMixingPoolSize;
    } else {
      fTreeStart = 0;
      fTreeEnd = fInputTree->GetEntries();
    }
    fIsFirstEvent = 0;
    cout << "fIsFirstEvent = 0" << endl;
    for (int icent = 0; icent < fCentArray.GetSize() -1; icent++) {
      for (int izvertex = 0; izvertex < fZVertexArray.GetSize() -1; izvertex++) {
        fPoolSize[icent][izvertex] = 0;
      }
    }
    //    TStopwatch stwa;
    //    stwa.Start();
    //    cout << "Start time : " << endl;
    for (int ientry = fTreeStart; ientry < fTreeEnd; ientry++) {
      fInputTree->GetEntry(ientry);
      if (!fHighPtMixing) {
        cout << "No High Pt Mixing! " << endl;
        return;
      }
      fPoolSize[GetCentBin(fHighPtMixing->GetCentrality())][GetZVertexBin(fHighPtMixing->GetZVertex())]++;
    }
  }

  if (fMCCorrection || fMCTruth) 
  {
    fMCEvent = MCEvent();
    fMCHeader = (AliAODMCHeader*) aodEvent->FindListObject(AliAODMCHeader::StdBranchName());
  }

  UInt_t mask;
  UInt_t bTargetEvents;

  // Event Selection Cut

  if (137161 <= runnumber && runnumber <= 139507) {// PbPb, 10h
    mask = static_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsEventSelected();
    bTargetEvents = mask & (AliVEvent::kMB);
    if (!bTargetEvents) return;
    fPeriod = k10h;
    fCollision = kPbPb;
  } else if (139510 <= runnumber && 137366 >= runnumber) { // PbPb, 11h
    mask = static_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsEventSelected();
    bTargetEvents = mask & (AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB);
    if (!bTargetEvents) return;
    fPeriod = k11h;
    fCollision = kPbPb;
  }
  else if ( runnumber <= 146860 && runnumber >= 146746) { // pp, 11a, 2.76TeV
    mask = static_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsEventSelected();
    bTargetEvents = mask & (AliVEvent::kMB);
    if (!bTargetEvents) return;
    fPeriod = k11a;
    fCollision = kpp;
  } else if (runnumber >= 244917 && runnumber <= 246994) {  // for Run2 dataset
    mask = static_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsEventSelected();
    bTargetEvents = mask & (AliVEvent::kINT7);
    if (!bTargetEvents) return;
    fPeriod = k15o;
    fCollision = kPbPb;
  } else {
    mask = static_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsEventSelected();
    bTargetEvents = mask & (AliVEvent::kINT7);
    if (!bTargetEvents) return;
    fPeriod = k17p;
    fCollision = kpp;

  }

  fEventPlane = -999;
  fEventPlaneV0A = -999;
  fEventPlaneV0C = -999;
  fEventPlaneTPC = -999;

  if (GetEventInformation(aodEvent)) return;

  AliCollisionGeometry *mcCollGeo;
  int nmctracks, label = -999;
  if (fMCCorrection || fMCTruth) {
    nmctracks = fMCEvent->GetNumberOfTracks();
  }

  TObjArray *lTrackArray = new TObjArray();
  lTrackArray->SetOwner(kTRUE);

  float ldeltaEventPlane = -999;
  float lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso = -999;
  float lChargeTrig, lChargeAsso = -999;
  int lpTTrigBin, lpTAssoBin = -999;
  int lCentBin, lZVertexBin = -1;
  //  int lTreeEntry = 0;
  float deltaEta, deltaPhi = -999;
  float ldPhiStar = -999;
  int lBSign = 0;

  lBSign = (event->GetMagneticField() > 0) ? 1 : -1;

  int ntracks = aodEvent->GetNumberOfTracks();
  short kPrimaryTrig, kPrimaryAsso = 0;

  int cont_inout = 0;
  double contvars[6] = {0,};
  double etavars[5] = {0,};

  etavars[1] = fCentPercentile;
  etavars[2] = fZVertex;

  fHistCent->Fill(fCentPercentile);
  fHistZVertex->Fill(fZVertex);
  fHistEventPlaneV0A->Fill(fEventPlaneV0A);
  fHistEventPlaneV0C->Fill(fEventPlaneV0C);
  fHistEventPlaneTPC->Fill(fEventPlaneTPC);
  lCentBin = GetCentBin(fCentPercentile);
  lZVertexBin = GetZVertexBin(fZVertex);

    if (fEntry%1000 == 0 ) cout << "Entry : " << fEntry << " Centrality : " << fCentPercentile << " ZVertex : " << fZVertex << " Evp : " << fEventPlaneV0A << endl;
  if (lCentBin == -1 || lZVertexBin == -1) return;   //To reduce time for mixing

  if (!fMCTruth) {
    if (fCollision == kPbPb && bUseMixingPool) {
      for (int itrack = 0; itrack < ntracks; itrack++) {
        AliAODTrack *aodTrack = static_cast<AliAODTrack*>(aodEvent->GetTrack(itrack));
        AliAODMCParticle *lmcTrack, *lmcTrackAsso;

        if (!aodTrack) {
          AliError(Form("ERROR: No track pointer for track %d", itrack));
          continue;
        }

        if (fFilterBit == 768){
          if (!(aodTrack->IsHybridGlobalConstrainedGlobal())) continue;  // new version of hybrid selection
          //      if (!aodTrack->TestFilterBit(fFilterBit)) continue;
        } else {
          //          if (aodTrack->GetType() != AliAODTrack::kPrimary) continue; // Check whether it's primary -> Leads to non flat phi distribution in hybrid tracks!
          if (!aodTrack->TestFilterBit(fFilterBit)) {
            //     std::cout << "filterbit no" << std::endl;
            continue; //TPC Only Track cut : 128 
          }
        }

        lChargeTrig = aodTrack->Charge();

        lpTTrig = aodTrack->Pt();

        ldeltaEventPlane = lPhiTrig - fEventPlaneV0C;

        if (lpTTrig > 0.5 && lpTTrig < 5) {
          fHistV2->Fill(fCentPercentile, TMath::Cos(2*ldeltaEventPlane));
        }

        if (lpTTrig < 0.8) continue;
        lEtaTrig = aodTrack->Eta();
        lPhiTrig = aodTrack->Phi();


        if (TMath::Abs(lEtaTrig) >= 0.8) continue;   

        fHistPhi->Fill(lPhiTrig);
        fHistPt->Fill(lpTTrig);
        fHistEta->Fill(lEtaTrig);
        if (lChargeTrig < 0) fHistPos->Fill(lpTTrig);
        if (lChargeTrig > 0) fHistNeg->Fill(lpTTrig);

        //      ldeltaEventPlane = lPhiTrig - fEventPlaneV0A;
        if (ldeltaEventPlane < 0) ldeltaEventPlane += TMath::TwoPi();
        if (ldeltaEventPlane > TMath::TwoPi()) ldeltaEventPlane -= TMath::TwoPi();

        if (fMCCorrection) {
          label = TMath::Abs(aodTrack->GetLabel());
          if (label > nmctracks) continue;
          lmcTrack = (AliAODMCParticle*) fMCEvent->GetTrack(label);
          kPrimaryTrig = lmcTrack->IsPhysicalPrimary();
        }

        etavars[0] = 0;
        etavars[3] = lEtaTrig;
        etavars[4] = lpTTrig;

        fHistPtSame[0]->Fill(lCentBin, lZVertexBin, lpTTrig);

        lTrackArray->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));

        lpTTrigBin = GetPttBin(lpTTrig);
        if (lpTTrigBin == -1) continue;

        fHistEtaSparse->Fill(etavars);

        cont_inout = CheckInOut(ldeltaEventPlane);

        fHistEtaSparse->Fill(etavars);

        if (fMCCorrection) {
          // Inclusive & All Reco particles (contaminated + primary)
          fHistContamination->Fill(contvars);
          if (kPrimaryTrig) {
            contvars[0] = 1;
            fHistContamination->Fill(contvars);
          } // In, out, m1, m2

          contvars[1] = 0; // Inclusive
          contvars[0] = 0; // contaminated + primary
          fHistContamination->Fill(contvars);
          if (kPrimaryTrig) {
            contvars[0] = 1;
            fHistContamination->Fill(contvars);
          }
        }

        for (int jtrack = 0; jtrack < ntracks; jtrack++) {
          if (itrack == jtrack) continue;
          AliAODTrack *aodTrackAsso = static_cast<AliAODTrack*>(aodEvent->GetTrack(jtrack));
          if (!aodTrackAsso) {
            continue;
          }

          if (fFilterBit == 768){
            if (!(aodTrackAsso->IsHybridGlobalConstrainedGlobal())) continue;  // new version of hybrid selection
            //        if (!(aodTrackAsso->TestFilterBit(fFilterBit))) continue;
          } else {
            //            if (aodTrackAsso->GetType() != AliAODTrack::kPrimary) continue;
            if (!aodTrackAsso->TestFilterBit(fFilterBit)) {
              //     std::cout << "filterbit no" << std::endl;
              continue; //TPC Only Track cut : 128 
            }
          }

          lpTAsso = aodTrackAsso->Pt();
          if (lpTAsso < 0.8 || lpTTrig < lpTAsso) continue;
          lpTAssoBin = GetPtaBin(lpTAsso);
          if (lpTAssoBin == -1) continue;

          lEtaAsso = aodTrackAsso->Eta();
          if (TMath::Abs(lEtaAsso) >= 0.8) continue;
          lPhiAsso = aodTrackAsso->Phi();
          lChargeAsso = aodTrackAsso->Charge();


          if (fConversionsVCut > 0 && lChargeAsso * lChargeTrig < 0)
          {
            Float_t mass = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.510e-3, 0.510e-3);

            if (mass < fConversionsVCut * 5)
            {
              mass = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.510e-3, 0.510e-3);

              if (mass < fConversionsVCut*fConversionsVCut) 
                continue;
            }
          }

          // K0s
          if (fResonancesVCut > 0 && lChargeAsso * lChargeTrig < 0)
          {
            Float_t mass = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.1396);

            const Float_t kK0smass = 0.4976;

            if (TMath::Abs(mass - kK0smass*kK0smass) < fResonancesVCut * 5)
            {
              mass = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.1396);

              if (mass > (kK0smass-fResonancesVCut)*(kK0smass-fResonancesVCut) && mass < (kK0smass+fResonancesVCut)*(kK0smass+fResonancesVCut))
                continue;
            }
          }


          // Lambda
          if (fResonancesVCut > 0 && lChargeAsso * lChargeTrig < 0)
          {
            Float_t mass1 = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.9383);
            Float_t mass2 = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.9383, 0.1396);

            const Float_t kLambdaMass = 1.115;

            if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
            {
              mass1 = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.9383);

              if (mass1 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass1 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                continue;
            }
            if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
            {
              mass2 = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.9383, 0.1396);

              if (mass2 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass2 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                continue;
            }
          }

          deltaEta = lEtaTrig - lEtaAsso;
          deltaPhi = lPhiTrig - lPhiAsso;
          if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
          if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

          ldPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, lChargeTrig, lChargeAsso, lpTTrig, lpTAsso, lBSign);
          if (TMath::Abs(ldPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
            continue;

          if (fMCCorrection) {
            label = TMath::Abs(aodTrackAsso->GetLabel());
            if (label > nmctracks) continue;
            lmcTrackAsso = (AliAODMCParticle*) fMCEvent->GetTrack(label);
            kPrimaryAsso = lmcTrackAsso->IsPhysicalPrimary();


            fHistdEtadPhiSameMCCorrCont[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            if (cont_inout > 0)
              fHistdEtadPhiSameMCCorrPrim[cont_inout][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            if (kPrimaryAsso && kPrimaryTrig) {
              fHistdEtadPhiSameMCCorrPrim[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
              if (cont_inout > 0)
                fHistdEtadPhiSameMCCorrPrim[cont_inout][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            }
          }

          fHistdEtadPhiSame[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);

          if (cont_inout > 0) fHistdEtadPhiSame[cont_inout][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);

        }

      }
    } else if (fCollision == kPbPb && !bUseMixingPool) {
      for (int itrack = 0; itrack < ntracks; itrack++) {
        AliAODTrack *aodTrack = static_cast<AliAODTrack*>(aodEvent->GetTrack(itrack));
        AliAODMCParticle *lmcTrack, *lmcTrackAsso;

        if (!aodTrack) {
          AliError(Form("ERROR: No track pointer for track %d", itrack));
          continue;
        }

        if (fFilterBit == 768){
          if (!(aodTrack->IsHybridGlobalConstrainedGlobal())) continue;  // new version of hybrid selection
          //      if (!aodTrack->TestFilterBit(fFilterBit)) continue;
        } else {
          //          if (aodTrack->GetType() != AliAODTrack::kPrimary) continue; // Check whether it's primary -> Leads to non flat phi distribution in hybrid tracks!
          if (!aodTrack->TestFilterBit(fFilterBit)) {
            //     std::cout << "filterbit no" << std::endl;
            continue; //TPC Only Track cut : 128 
          }
        }

        lChargeTrig = aodTrack->Charge();

        lpTTrig = aodTrack->Pt();

        ldeltaEventPlane = lPhiTrig - fEventPlaneV0C;

        if (lpTTrig > 0.5 && lpTTrig < 5) {
          fHistV2->Fill(fCentPercentile, TMath::Cos(2*ldeltaEventPlane));
        }
        if (lpTTrig < 0.8) continue;
        lpTTrigBin = GetPttBin(lpTTrig);
        lTrackArray->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));
        lEtaTrig = aodTrack->Eta();
        lPhiTrig = aodTrack->Phi();


        if (TMath::Abs(lEtaTrig) >= 0.8) continue;   

        fHistPhi->Fill(lPhiTrig);
        fHistPt->Fill(lpTTrig);
        fHistEta->Fill(lEtaTrig);
        if (lChargeTrig < 0) fHistPos->Fill(lpTTrig);
        if (lChargeTrig > 0) fHistNeg->Fill(lpTTrig);

        //      ldeltaEventPlane = lPhiTrig - fEventPlaneV0A;
        if (ldeltaEventPlane < 0) ldeltaEventPlane += TMath::TwoPi();
        if (ldeltaEventPlane > TMath::TwoPi()) ldeltaEventPlane -= TMath::TwoPi();

        if (fMCCorrection) {
          label = TMath::Abs(aodTrack->GetLabel());
          if (label > nmctracks) continue;
          lmcTrack = (AliAODMCParticle*) fMCEvent->GetTrack(label);
          kPrimaryTrig = lmcTrack->IsPhysicalPrimary();
        }
        etavars[0] = 0;

        etavars[3] = lEtaTrig;
        etavars[4] = lpTTrig;

        fHistPtSame[0]->Fill(lCentBin, lZVertexBin, lpTTrig);
        fHistEtaSparse->Fill(etavars);

        cont_inout = CheckInOut(ldeltaEventPlane);
        if (cont_inout > 0) {
          contvars[1] = cont_inout;
          etavars[0] = cont_inout;
          fHistPtSame[cont_inout]->Fill(lCentBin, lZVertexBin, lpTTrig);
        }

        fHistEtaSparse->Fill(etavars);

        if (fMCCorrection) {
          // Inclusive & All Reco particles (contaminated + primary)
          contvars[0] = 0;
          contvars[2] = fCentPercentile;
          contvars[3] = fZVertex;
          contvars[4] = lEtaTrig;
          contvars[5] = lpTTrig;

          fHistContamination->Fill(contvars);
          if (kPrimaryTrig) {
            contvars[0] = 1;
            fHistContamination->Fill(contvars);
          } // In, out, m1, m2

          contvars[1] = 0; // Inclusive
          contvars[0] = 0; // contaminated + primary
          fHistContamination->Fill(contvars);
          if (kPrimaryTrig) {
            contvars[0] = 1;
            fHistContamination->Fill(contvars);
          }
        }

        if (lpTTrigBin == -1) continue;

        for (int jtrack = 0; jtrack < ntracks; jtrack++) {
          if (itrack == jtrack) continue;
          AliAODTrack *aodTrackAsso = static_cast<AliAODTrack*>(aodEvent->GetTrack(jtrack));
          if (!aodTrackAsso) {
            continue;
          }

          if (fFilterBit == 768){
            if (!(aodTrackAsso->IsHybridGlobalConstrainedGlobal())) continue;  // new version of hybrid selection
            //        if (!(aodTrackAsso->TestFilterBit(fFilterBit))) continue;
          } else {
            //            if (aodTrackAsso->GetType() != AliAODTrack::kPrimary) continue;
            if (!aodTrackAsso->TestFilterBit(fFilterBit)) {
              //     std::cout << "filterbit no" << std::endl;
              continue; //TPC Only Track cut : 128 
            }
          }

          lpTAsso = aodTrackAsso->Pt();
          if (lpTAsso < 0.8 || lpTTrig < lpTAsso) continue;
          lpTAssoBin = GetPtaBin(lpTAsso);
          if (lpTAssoBin == -1) continue;

          lEtaAsso = aodTrackAsso->Eta();
          if (TMath::Abs(lEtaAsso) >= 0.8) continue;
          lPhiAsso = aodTrackAsso->Phi();
          lChargeAsso = aodTrackAsso->Charge();


          if (fConversionsVCut > 0 && lChargeAsso * lChargeTrig < 0)
          {
            Float_t mass = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.510e-3, 0.510e-3);

            if (mass < fConversionsVCut * 5)
            {
              mass = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.510e-3, 0.510e-3);

              if (mass < fConversionsVCut*fConversionsVCut) 
                continue;
            }
          }

          // K0s
          if (fResonancesVCut > 0 && lChargeAsso * lChargeTrig < 0)
          {
            Float_t mass = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.1396);

            const Float_t kK0smass = 0.4976;

            if (TMath::Abs(mass - kK0smass*kK0smass) < fResonancesVCut * 5)
            {
              mass = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.1396);

              if (mass > (kK0smass-fResonancesVCut)*(kK0smass-fResonancesVCut) && mass < (kK0smass+fResonancesVCut)*(kK0smass+fResonancesVCut))
                continue;
            }
          }


          // Lambda
          if (fResonancesVCut > 0 && lChargeAsso * lChargeTrig < 0)
          {
            Float_t mass1 = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.9383);
            Float_t mass2 = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.9383, 0.1396);

            const Float_t kLambdaMass = 1.115;

            if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
            {
              mass1 = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.9383);

              if (mass1 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass1 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                continue;
            }
            if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
            {
              mass2 = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.9383, 0.1396);

              if (mass2 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass2 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                continue;
            }
          }

          deltaEta = lEtaTrig - lEtaAsso;
          deltaPhi = lPhiTrig - lPhiAsso;
          if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
          if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

          ldPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, lChargeTrig, lChargeAsso, lpTTrig, lpTAsso, lBSign);
          if (TMath::Abs(ldPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
            continue;

//          cout << fMCCorrection << endl;
          if (fMCCorrection) {
            label = TMath::Abs(aodTrackAsso->GetLabel());
            nmctracks = fMCEvent->GetNumberOfTracks();
            if (label > nmctracks) continue;
            lmcTrackAsso = (AliAODMCParticle*) fMCEvent->GetTrack(label);
            kPrimaryAsso = lmcTrackAsso->IsPhysicalPrimary();
            

            fHistdEtadPhiSameMCCorrCont[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            if (cont_inout > 0) fHistdEtadPhiSameMCCorrCont[cont_inout][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);

            if (kPrimaryAsso && kPrimaryTrig) {
              fHistdEtadPhiSameMCCorrPrim[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
              if (cont_inout > 0) fHistdEtadPhiSameMCCorrPrim[cont_inout][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            }
          }

          fHistdEtadPhiSame[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          if (cont_inout > 0) 
            fHistdEtadPhiSame[cont_inout][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);

        }

      }

    } else {
      for (int itrack = 0; itrack < ntracks; itrack++) {
        AliAODTrack *aodTrack = static_cast<AliAODTrack*>(aodEvent->GetTrack(itrack));
        AliAODMCParticle *lmcTrack, *lmcTrackAsso;

        if (!aodTrack) {
          AliError(Form("ERROR: No track pointer for track %d", itrack));
          continue;
        }

        if (fFilterBit == 768){
          if (!(aodTrack->IsHybridGlobalConstrainedGlobal())) {
            continue; 
          }// new version of hybrid selection
          //      if (!aodTrack->TestFilterBit(fFilterBit)) continue;
        } else {
          //          if (aodTrack->GetType() != AliAODTrack::kPrimary) continue; // Check whether it's primary -> Leads to non flat phi distribution in hybrid tracks!
          if (!aodTrack->TestFilterBit(fFilterBit)) {
            //     std::cout << "filterbit no" << std::endl;
            continue; //TPC Only Track cut : 128 
          }
        }

        lChargeTrig = aodTrack->Charge();

        lpTTrig = aodTrack->Pt();
        if (lpTTrig < 0.8) continue;
        lEtaTrig = aodTrack->Eta();
        lPhiTrig = aodTrack->Phi();
        lpTTrigBin = GetPttBin(lpTTrig);
        lTrackArray->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));
        if (lpTTrigBin == -1) continue;


        if (TMath::Abs(lEtaTrig) >= 0.8) continue;   

        fHistPhi->Fill(lPhiTrig);
        fHistPt->Fill(lpTTrig);
        fHistEta->Fill(lEtaTrig);
        if (lChargeTrig < 0) fHistPos->Fill(lpTTrig);
        if (lChargeTrig > 0) fHistNeg->Fill(lpTTrig);

        if (fMCCorrection) {
          label = TMath::Abs(aodTrack->GetLabel());
          nmctracks = fMCEvent->GetNumberOfTracks();
          if (label > nmctracks) continue;
          lmcTrack = (AliAODMCParticle*) fMCEvent->GetTrack(label);
          kPrimaryTrig = lmcTrack->IsPhysicalPrimary();
        }

        etavars[0] = 0;
        etavars[3] = lEtaTrig;
        etavars[4] = lpTTrig;

        fHistPtSame[0]->Fill(lCentBin, lZVertexBin, lpTTrig);
        fHistEtaSparse->Fill(etavars);

        if (fMCCorrection) {
          // Inclusive & All Reco particles (contaminated + primary)
          contvars[0] = 0;
          contvars[2] = fCentPercentile;
          contvars[3] = fZVertex;
          contvars[4] = lEtaTrig;
          contvars[5] = lpTTrig;

          fHistContamination->Fill(contvars);
          if (kPrimaryTrig) {
            contvars[0] = 1;
            fHistContamination->Fill(contvars);
          } // In, out, m1, m2

          contvars[1] = 0; // Inclusive
          contvars[0] = 0; // contaminated + primary
          fHistContamination->Fill(contvars);
          if (kPrimaryTrig) {
            contvars[0] = 1;
            fHistContamination->Fill(contvars);
          }
        }


        for (int jtrack = 0; jtrack < ntracks; jtrack++) {
          if (itrack == jtrack) continue;
          AliAODTrack *aodTrackAsso = static_cast<AliAODTrack*>(aodEvent->GetTrack(jtrack));
          if (!aodTrackAsso) {
            continue;
          }

          if (fFilterBit == 768){
            if (!(aodTrackAsso->IsHybridGlobalConstrainedGlobal())) continue;  // new version of hybrid selection
            //        if (!(aodTrackAsso->TestFilterBit(fFilterBit))) continue;
          } else {
            //            if (aodTrackAsso->GetType() != AliAODTrack::kPrimary) continue;
            if (!aodTrackAsso->TestFilterBit(fFilterBit)) {
              //     std::cout << "filterbit no" << std::endl;
              continue; //TPC Only Track cut : 128 
            }
          }

          lpTAsso = aodTrackAsso->Pt();
          if (lpTAsso < 0.8 || lpTTrig < lpTAsso) continue;
          lpTAssoBin = GetPtaBin(lpTAsso);
          if (lpTAssoBin == -1) continue;

          lEtaAsso = aodTrackAsso->Eta();
          if (TMath::Abs(lEtaAsso) >= 0.8) continue;
          lPhiAsso = aodTrackAsso->Phi();
          lChargeAsso = aodTrackAsso->Charge();


          if (fConversionsVCut > 0 && lChargeAsso * lChargeTrig < 0)
          {
            Float_t mass = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.510e-3, 0.510e-3);

            if (mass < fConversionsVCut * 5)
            {
              mass = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.510e-3, 0.510e-3);

              if (mass < fConversionsVCut*fConversionsVCut) 
                continue;
            }
          }

          // K0s
          if (fResonancesVCut > 0 && lChargeAsso * lChargeTrig < 0)
          {
            Float_t mass = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.1396);

            const Float_t kK0smass = 0.4976;

            if (TMath::Abs(mass - kK0smass*kK0smass) < fResonancesVCut * 5)
            {
              mass = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.1396);

              if (mass > (kK0smass-fResonancesVCut)*(kK0smass-fResonancesVCut) && mass < (kK0smass+fResonancesVCut)*(kK0smass+fResonancesVCut))
                continue;
            }
          }


          // Lambda
          if (fResonancesVCut > 0 && lChargeAsso * lChargeTrig < 0)
          {
            Float_t mass1 = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.9383);
            Float_t mass2 = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.9383, 0.1396);

            const Float_t kLambdaMass = 1.115;

            if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
            {
              mass1 = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.9383);

              if (mass1 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass1 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                continue;
            }
            if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
            {
              mass2 = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.9383, 0.1396);

              if (mass2 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass2 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                continue;
            }
          }



          deltaEta = lEtaTrig - lEtaAsso;
          deltaPhi = lPhiTrig - lPhiAsso;
          if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
          if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

          ldPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, lChargeTrig, lChargeAsso, lpTTrig, lpTAsso, lBSign);
          if (TMath::Abs(ldPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
            continue;

          if (fMCCorrection) {
            label = TMath::Abs(aodTrackAsso->GetLabel());
            nmctracks = fMCEvent->GetNumberOfTracks();
            if (label > nmctracks) continue;
            lmcTrackAsso = (AliAODMCParticle*) fMCEvent->GetTrack(label);
            kPrimaryAsso = lmcTrackAsso->IsPhysicalPrimary();


            fHistdEtadPhiSameMCCorrCont[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            if (kPrimaryAsso && kPrimaryTrig) {
              fHistdEtadPhiSameMCCorrPrim[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            }
          }

          fHistdEtadPhiSame[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        }

      }
    }  // fCollision = pp
  } else {

    for (int itrack = 0; itrack < nmctracks; itrack++) {
      AliAODMCParticle *lmcTrack = (AliAODMCParticle*) fMCEvent->GetTrack(itrack);
      if (!lmcTrack) continue;
      lChargeTrig = lmcTrack->Charge();
      if (!lChargeTrig) continue;
      int pdg = TMath::Abs(lmcTrack->GetPdgCode());
      if ( !((pdg >= 211 && pdg <= 533) || (pdg > 1000 && pdg < 6000))) continue; 
      if (! (lmcTrack->IsPhysicalPrimary()) ) continue;

      lEtaTrig = lmcTrack->Eta();
      lpTTrig = lmcTrack->Pt();
      lPhiTrig = lmcTrack->Phi();
      ldeltaEventPlane = lPhiTrig - fEventPlaneV0C;

      if (TMath::Abs(lEtaTrig) >= 0.8) continue;
      lpTTrigBin = GetPttBin(lpTTrig);
      lTrackArray->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));
      if (lpTTrigBin == -1) continue;

      fHistPhi->Fill(lPhiTrig);
      fHistPt->Fill(lpTTrig);
      fHistEta->Fill(lEtaTrig);
      if (ldeltaEventPlane < 0) ldeltaEventPlane += TMath::TwoPi();
      if (ldeltaEventPlane > TMath::TwoPi()) ldeltaEventPlane -= TMath::TwoPi();

      fHistPtSame[0]->Fill(lCentBin, lZVertexBin, lpTTrig);
      cont_inout = CheckInOut(ldeltaEventPlane);
      if (cont_inout > 0) {
        contvars[1] = cont_inout;
        etavars[0] = cont_inout;
        fHistPtSame[cont_inout]->Fill(lCentBin, lZVertexBin, lpTTrig);
      }


      for (int jtrack = 0; jtrack < nmctracks; jtrack++) {
        if (itrack == jtrack) continue;
        AliAODMCParticle *lmcTrackAsso = dynamic_cast<AliAODMCParticle*>(fMCEvent->GetTrack(jtrack));
        if (!lmcTrackAsso) continue;

        if (!(lmcTrackAsso->IsPhysicalPrimary())) continue;
        int pdgAsso = TMath::Abs(lmcTrackAsso->GetPdgCode());
        if ( !((pdgAsso >= 211 && pdgAsso <= 533) || (pdgAsso > 1000 && pdgAsso < 6000))) continue; 
        lChargeAsso = lmcTrackAsso->Charge();
        if (!lChargeAsso) continue;
        lEtaAsso = lmcTrackAsso->Eta();
        lpTAsso = lmcTrackAsso->Pt();
        lPhiAsso = lmcTrackAsso->Phi();
        if (lpTTrig < lpTAsso) continue;

        if (TMath::Abs(lEtaAsso) >= 0.8) continue;
        lpTAssoBin = GetPtaBin(lpTAsso);
        if (lpTAssoBin == -1) continue;

        if (fConversionsVCut > 0 && lChargeAsso * lChargeTrig < 0)
        {
          Float_t mass = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.510e-3, 0.510e-3);

          if (mass < fConversionsVCut * 5)
          {
            mass = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.510e-3, 0.510e-3);

            if (mass < fConversionsVCut*fConversionsVCut) 
              continue;
          }
        }

        // K0s
        if (fResonancesVCut > 0 && lChargeAsso * lChargeTrig < 0)
        {
          Float_t mass = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.1396);

          const Float_t kK0smass = 0.4976;

          if (TMath::Abs(mass - kK0smass*kK0smass) < fResonancesVCut * 5)
          {
            mass = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.1396);

            if (mass > (kK0smass-fResonancesVCut)*(kK0smass-fResonancesVCut) && mass < (kK0smass+fResonancesVCut)*(kK0smass+fResonancesVCut))
              continue;
          }
        }


        // Lambda
        if (fResonancesVCut > 0 && lChargeAsso * lChargeTrig < 0)
        {
          Float_t mass1 = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.9383);
          Float_t mass2 = GetInvMassSquaredCheap(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.9383, 0.1396);

          const Float_t kLambdaMass = 1.115;

          if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
          {
            mass1 = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.1396, 0.9383);

            if (mass1 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass1 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
              continue;
          }
          if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
          {
            mass2 = GetInvMassSquared(lpTTrig, lEtaTrig, lPhiTrig, lpTAsso, lEtaAsso, lPhiAsso, 0.9383, 0.1396);

            if (mass2 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass2 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
              continue;
          }
        }

        deltaEta = lEtaTrig - lEtaAsso;
        deltaPhi = lPhiTrig - lPhiAsso;
        if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
        if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

        ldPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, lChargeTrig, lChargeAsso, lpTTrig, lpTAsso, lBSign);
        if (TMath::Abs(ldPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
          continue;

        //        cout << lCentBin << " " << lZVertexBin << " " << lpTTrigBin << " " << lpTAssoBin << endl;
        fHistdEtadPhiSame[0][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        if (cont_inout > 0)
          fHistdEtadPhiSame[cont_inout][lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
      } // jtrack mc


    }  //itrack mc

  } // fMCTruth

  //  cout << "End of Event" << endl;
  DoMixing(fCentPercentile, fZVertex, lTrackArray, fEventPlaneV0C, lBSign);
  fHistNevtSame->Fill(lCentBin, lZVertexBin);


  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    PostData(iplane+1, fOutput[iplane]);
  }
}

int AliAnalysisTaskJetLikeCorrelation::SetupMixing() {
  double *centralityarr = fCentArray.GetArray();
  double *zVertexarr = fZVertexArray.GetArray();
  int nCentBin = fCentArray.GetSize();
  int nZVertBin = fZVertexArray.GetSize();

  if (fCollision == kPbPb) {
    int poolsize = 100;
    fPoolMgr = new AliEventPoolManager(poolsize, fTrackDepth, nCentBin, centralityarr, nZVertBin, zVertexarr);
  } else if (fCollision == kpPb) {
    int poolsize = 150;
    fPoolMgr = new AliEventPoolManager(poolsize, fTrackDepth, nCentBin, centralityarr, nZVertBin, zVertexarr);
  } else if ( fCollision == kpp) {
    int poolsize = 300;
    cout << "Collision pp, making pool of size : " << poolsize << endl;
    fPoolMgr = new AliEventPoolManager(poolsize, fTrackDepth, nCentBin, centralityarr, nZVertBin, zVertexarr);
  } else {
    AliError("[SetupMixing] Incorrect Collision Type! Aborting..");
    return 1;
  }
  return 0;

}

int AliAnalysisTaskJetLikeCorrelation::DoMixing(float lCentrality, float fZVertex, TObjArray *lTracksTrig, float eventPlane, int lbSign) {

  //  cout << "Starting Do Mixing.. " << endl;
  int lMinimumPtTBinMerging = GetPttBin(6.00);

  if (!fMCTruth) {
    if (bUseMixingPool) {
      //    lTracksMixing->SetOwner();
      //    int lnEvent = fTreeEnd - fTreeStart + 1; // = pool->GetCurrentNEvents();

      int lCentBin = GetCentBin(lCentrality);
      int lZVertexBin = GetZVertexBin(fZVertex);
      int lnEvent = fPoolSize[lCentBin][lZVertexBin];
      long lnTracksMixing = 0;
      long lnTracksTrig = 0;

      float trigpT, trigEta, trigPhi, trigCharge = -999;
      float mixpT, mixEta, mixPhi, mixCharge = -999;
      unsigned long mixFilterBit = 0;
      float deltaPhi, deltaEta = -999;
      float dPhiStar = -999;
      float ldeltaEventPlane = -999;
      int trigpTBin, mixpTBin;
      int evp_inout = 0;


      double tempcent, tempzvertex;
      unsigned long long tempeventid;

      vector<AliExtractedTrack> ltracksmixing;
      if (!fInputTree){
        cout << "No Input Tree! returning! " << endl;
        return 999;
      }

      for (int ientry = fTreeStart; ientry < fTreeEnd; ientry++) {
        fInputTree->GetEntry(ientry);
        if (!fHighPtMixing) {
          cout << "No High Pt Mixing! " << endl;
          return 998;
        }
        tempeventid = fHighPtMixing->GetEventID();
        if (tempeventid == fEventID) lnEvent -= 1;
      }


      if (fCollision == kPbPb) {      
        for (int ientry = fTreeStart; ientry < fTreeEnd; ientry++) {
          fInputTree->GetEntry(ientry);

          tempcent = fHighPtMixing->GetCentrality();
          tempzvertex = fHighPtMixing->GetZVertex();
          tempeventid = fHighPtMixing->GetEventID();


          if (GetCentBin(tempcent) != lCentBin || GetZVertexBin(tempzvertex) != lZVertexBin || tempeventid == fEventID) continue;
          ltracksmixing = fHighPtMixing->GetTracks();
          lnTracksMixing = ltracksmixing.size();
          //    if (lnTracksMixing < 1) continue;
          lnTracksTrig = lTracksTrig->GetEntriesFast();
          //    cout << lCentBin << " " << lZVertexBin << " " << lnTracksMixing << endl;

          for (int itrack = 0; itrack < lnTracksTrig; itrack++) {
            if (!(lTracksTrig->At(itrack))) return 20;
            const MixedParticle *lParticleTrig = dynamic_cast<MixedParticle*>(lTracksTrig->At(itrack));
            trigpT = lParticleTrig->Pt();
            trigEta = lParticleTrig->Eta();
            if (TMath::Abs(trigEta) >= 0.8) continue;
            trigPhi = lParticleTrig->Phi();
            trigCharge = lParticleTrig->Charge();
            //        cout << ldeltaEventPlane << "\n";
            trigpTBin = GetPttBin(trigpT);
            if (trigpTBin == -1) continue;

            ldeltaEventPlane = trigPhi - eventPlane;
            if (ldeltaEventPlane < 0) ldeltaEventPlane += TMath::TwoPi();
            if (ldeltaEventPlane > TMath::TwoPi()) ldeltaEventPlane -= TMath::TwoPi();

            evp_inout = CheckInOut(ldeltaEventPlane);

            for (int jtrack = 0; jtrack < lnTracksMixing; jtrack++) {
              mixFilterBit = ltracksmixing[jtrack].GetFilterBit();
              if (!(mixFilterBit & fFilterBit)) continue;       // Test filter bit
//              if (mixFilterBit == 128 && !(ltracksmixing[jtrack].IsPrimary())) continue;     // Only for TPC only track primary check
              mixpT = ltracksmixing[jtrack].GetPt();
              if (mixpT < 3.0 ||  trigpT < mixpT) continue;  // for pbpb
              //          if (trigpT < mixpT) continue;  // for pp
              mixEta = ltracksmixing[jtrack].GetEta();
              if (TMath::Abs(mixEta) >= 0.8) continue;
              mixPhi = ltracksmixing[jtrack].GetPhi();
              mixCharge = ltracksmixing[jtrack].GetCharge();
              mixpTBin = GetPtaBin(mixpT);

              if (mixpTBin == -1) continue;

              if (fConversionsVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.510e-3, 0.510e-3);
                if (mass < fConversionsVCut * 5)
                {
                  mass = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.510e-3, 0.510e-3);
                  if (mass < fConversionsVCut*fConversionsVCut) 
                    continue;
                }
              }

              // K0s
              if (fResonancesVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.1396);

                const Float_t kK0smass = 0.4976;

                if (TMath::Abs(mass - kK0smass*kK0smass) < fResonancesVCut * 5)
                {
                  mass = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.1396);

                  if (mass > (kK0smass-fResonancesVCut)*(kK0smass-fResonancesVCut) && mass < (kK0smass+fResonancesVCut)*(kK0smass+fResonancesVCut))
                    continue;
                }
              }


              // Lambda
              if (fResonancesVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass1 = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.9383);
                Float_t mass2 = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.9383, 0.1396);

                const Float_t kLambdaMass = 1.115;

                if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
                {
                  mass1 = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.9383);

                  if (mass1 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass1 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                    continue;
                }
                if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
                {
                  mass2 = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.9383, 0.1396);

                  if (mass2 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass2 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                    continue;
                }
              }

              deltaPhi = trigPhi - mixPhi;
              deltaEta = trigEta - mixEta;

              if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
              if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

              dPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, trigCharge, mixCharge, trigpT, mixpT, lbSign);
              if (TMath::Abs(dPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
                continue;

              if (trigpTBin >= lMinimumPtTBinMerging && mixpTBin >= fMinimumPtABinMerging) {
                fHistdEtadPhiMixed[0][lCentBin][lZVertexBin][trigpTBin][fMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
                if (evp_inout > 0) fHistdEtadPhiMixed[evp_inout][lCentBin][lZVertexBin][trigpTBin][fMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

              } else {
                fHistdEtadPhiMixed[0][lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
                if (evp_inout > 0) fHistdEtadPhiMixed[evp_inout][lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

              }
            }
          }
        }
      }// fCollision == kPbPb 
      else // fCollision == pp
      {
        for (int ientry = fTreeStart; ientry < fTreeEnd; ientry++) {
          fInputTree->GetEntry(ientry);

          tempcent = fHighPtMixing->GetCentrality();
          tempzvertex = fHighPtMixing->GetZVertex();
          tempeventid = fHighPtMixing->GetEventID();

          if (GetCentBin(tempcent) != lCentBin || GetZVertexBin(tempzvertex) != lZVertexBin || tempeventid == fEventID) continue;
          ltracksmixing = fHighPtMixing->GetTracks();
          lnTracksMixing = ltracksmixing.size();
          //    if (lnTracksMixing < 1) continue;
          lnTracksTrig = lTracksTrig->GetEntriesFast();
          //    cout << lCentBin << " " << lZVertexBin << " " << lnTracksMixing << endl;
          for (int itrack = 0; itrack < lnTracksTrig; itrack++) {
            if (!(lTracksTrig->At(itrack))) return 20;
            const MixedParticle *lParticleTrig = dynamic_cast<MixedParticle*>(lTracksTrig->At(itrack));
            trigpT = lParticleTrig->Pt();
            trigEta = lParticleTrig->Eta();
            if (TMath::Abs(trigEta) >= 0.8) continue;
            trigPhi = lParticleTrig->Phi();
            trigCharge = lParticleTrig->Charge();

            //        cout << ldeltaEventPlane << "\n";
            trigpTBin = GetPttBin(trigpT);
            if (trigpTBin == -1) continue;

            for (int jtrack = 0; jtrack < lnTracksMixing; jtrack++) {
              mixFilterBit = ltracksmixing[jtrack].GetFilterBit();
              if (!(mixFilterBit & fFilterBit)) continue;       // Test filter bit
              if (mixFilterBit == 128 && !(ltracksmixing[jtrack].IsPrimary())) continue;     // Only for TPC only track primary check
              mixpT = ltracksmixing[jtrack].GetPt();
              //                    if (mixpT < 4.0 ||  trigpT < mixpT) continue;  // for pbpb
              if (trigpT < mixpT || mixpT < 0.8) continue;  // for pp
              mixEta = ltracksmixing[jtrack].GetEta();
              if (TMath::Abs(mixEta) >= 0.8) continue;
              mixPhi = ltracksmixing[jtrack].GetPhi();
              mixCharge = ltracksmixing[jtrack].GetCharge();
              mixpTBin = GetPtaBin(mixpT);

              if (mixpTBin == -1) continue;

              if (fConversionsVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.510e-3, 0.510e-3);
                if (mass < fConversionsVCut * 5)
                {
                  mass = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.510e-3, 0.510e-3);
                  if (mass < fConversionsVCut*fConversionsVCut) 
                    continue;
                }
              }

              // K0s
              if (fResonancesVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.1396);

                const Float_t kK0smass = 0.4976;

                if (TMath::Abs(mass - kK0smass*kK0smass) < fResonancesVCut * 5)
                {
                  mass = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.1396);

                  if (mass > (kK0smass-fResonancesVCut)*(kK0smass-fResonancesVCut) && mass < (kK0smass+fResonancesVCut)*(kK0smass+fResonancesVCut))
                    continue;
                }
              }


              // Lambda
              if (fResonancesVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass1 = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.9383);
                Float_t mass2 = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.9383, 0.1396);

                const Float_t kLambdaMass = 1.115;

                if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
                {
                  mass1 = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.9383);

                  if (mass1 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass1 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                    continue;
                }
                if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
                {
                  mass2 = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.9383, 0.1396);

                  if (mass2 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass2 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                    continue;
                }
              }

              deltaPhi = trigPhi - mixPhi;
              deltaEta = trigEta - mixEta;

              if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
              if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

              dPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, trigCharge, mixCharge, trigpT, mixpT, lbSign);
              if (TMath::Abs(dPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
                continue;

              if (trigpTBin >= lMinimumPtTBinMerging && mixpTBin >= fMinimumPtABinMerging) {
                fHistdEtadPhiMixed[0][lCentBin][lZVertexBin][trigpTBin][fMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
              } else { 
                fHistdEtadPhiMixed[0][lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
              }

            }
          }
        }
      }
    } else {
      AliEventPool *pool = fPoolMgr->GetEventPool(lCentrality, fZVertex);
      if (!pool) return 1;

      TObjArray *lTracksMixing;
      int lnEvent = pool->GetCurrentNEvents();

      int lCentBin = GetCentBin(lCentrality);
      int lZVertexBin = GetZVertexBin(fZVertex);
      long lnTracksMixing = 0;
      long lnTracksTrig = 0;

      float trigpT, trigEta, trigPhi, trigCharge = -999;
      float mixpT, mixEta, mixPhi, mixCharge = -999;
      float deltaPhi, deltaEta = -999;
      float dPhiStar = -999;
      float ldeltaEventPlane = -999;
      int trigpTBin, mixpTBin;
      int evp_inout = 0;

      //  cout << "pool : " << lnEvent << " " << pool->IsReady() << endl;

      if (pool->IsReady() || pool->GetCurrentNEvents() >= 5) {
        if (fDebug > 4) pool->PrintInfo();

        for (int ievent = 0; ievent < lnEvent; ievent++) {
          lTracksMixing = pool->GetEvent(ievent);
          //      cout << "lNEvt : " << lnEvent << endl;
          if (!lTracksMixing) continue;
          lnTracksMixing = lTracksMixing->GetEntriesFast();
          lnTracksTrig = lTracksTrig->GetEntriesFast();
          //      cout << "lNEvt : " << lnEvent << " MixingT " << lnTracksMixing << " Trig " << lnTracksTrig << endl;

          for (int itrack = 0; itrack < lnTracksTrig; itrack++) {
            if (!(lTracksTrig->At(itrack))) return 1;
            const MixedParticle *lParticleTrig = dynamic_cast<MixedParticle*>(lTracksTrig->At(itrack));
            trigpT = lParticleTrig->Pt();
            if (trigpT < fMinPtTrigCut) continue;
            trigEta = lParticleTrig->Eta();
            if (TMath::Abs(trigEta) >= 0.8) continue;
            trigPhi = lParticleTrig->Phi();
            trigCharge = lParticleTrig->Charge();
            ldeltaEventPlane = trigPhi - eventPlane;
            if (ldeltaEventPlane < 0) ldeltaEventPlane += TMath::TwoPi();
            if (ldeltaEventPlane > TMath::TwoPi()) ldeltaEventPlane -= TMath::TwoPi();
            //        cout << ldeltaEventPlane << "\n";
            trigpTBin = GetPttBin(trigpT);

            if (trigpTBin == -1) continue;
            evp_inout = CheckInOut(ldeltaEventPlane);

            for (int jtrack = 0; jtrack < lnTracksMixing; jtrack++) {
              if (!(lTracksMixing->At(jtrack))) return 1;
              const MixedParticle *lParticleMixing = dynamic_cast<MixedParticle*>(lTracksMixing->At(jtrack));
              //          cout << "PROBLEM1 ?" << endl;
              mixpT = lParticleMixing->Pt();
              if (trigpT < mixpT) continue;
              mixEta = lParticleMixing->Eta();
              if (TMath::Abs(mixEta) >= 0.8) continue;
              mixPhi = lParticleMixing->Phi();
              mixCharge = lParticleMixing->Charge();
              mixpTBin = GetPtaBin(mixpT);

              if (mixpTBin == -1) continue;

              if (fConversionsVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.510e-3, 0.510e-3);

                if (mass < fConversionsVCut * 5)
                {
                  mass = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.510e-3, 0.510e-3);

                  if (mass < fConversionsVCut*fConversionsVCut) 
                    continue;
                }
              }

              // K0s
              if (fResonancesVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.1396);

                const Float_t kK0smass = 0.4976;

                if (TMath::Abs(mass - kK0smass*kK0smass) < fResonancesVCut * 5)
                {
                  mass = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.1396);

                  if (mass > (kK0smass-fResonancesVCut)*(kK0smass-fResonancesVCut) && mass < (kK0smass+fResonancesVCut)*(kK0smass+fResonancesVCut))
                    continue;
                }
              }


              // Lambda
              if (fResonancesVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass1 = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.9383);
                Float_t mass2 = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.9383, 0.1396);

                const Float_t kLambdaMass = 1.115;

                if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
                {
                  mass1 = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.9383);

                  if (mass1 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass1 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                    continue;
                }
                if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
                {
                  mass2 = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.9383, 0.1396);

                  if (mass2 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass2 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                    continue;
                }
              }


              deltaPhi = trigPhi - mixPhi;
              deltaEta = trigEta - mixEta;

//                        cout << evp_inout << " " << mixpTBin << " and " << trigpTBin << "  " << jtrack << "/" << lnTracksMixing << " and " << itrack << "/" << lnTracksTrig  <<  endl;

              if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
              if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

              dPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, trigCharge, mixCharge, trigpT, mixpT, lbSign);
              if (TMath::Abs(dPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
                continue;

              if (trigpTBin >= lMinimumPtTBinMerging && mixpTBin >= fMinimumPtABinMerging) {
                fHistdEtadPhiMixed[0][lCentBin][lZVertexBin][trigpTBin][fMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
                if (evp_inout > 0) fHistdEtadPhiMixed[evp_inout][lCentBin][lZVertexBin][trigpTBin][fMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

              } else {
                fHistdEtadPhiMixed[0][lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
                if (evp_inout > 0)fHistdEtadPhiMixed[evp_inout][lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

              }
              //          cout << "Ho ! " << mixPhi <<"  " <<  mixCharge << "  " << mixEta << "  "<< endl;

            } // jtrack

          } // itrack
        }  // ievent
        //    cout << "problem!" << endl;
        pool->UpdatePool(lTracksTrig);
      } else {
        pool->UpdatePool(lTracksTrig);
        //    cout << lTracksTrig->GetEntriesFast() << endl;
      } // if pool is ready

    } 
  } else {

      AliEventPool *pool = fPoolMgr->GetEventPool(lCentrality, fZVertex);
      if (!pool) return 1;

      TObjArray *lTracksMixing;
      int lnEvent = pool->GetCurrentNEvents();

      int lCentBin = GetCentBin(lCentrality);
      int lZVertexBin = GetZVertexBin(fZVertex);
      long lnTracksMixing = 0;
      long lnTracksTrig = 0;

      float trigpT, trigEta, trigPhi, trigCharge = -999;
      float mixpT, mixEta, mixPhi, mixCharge = -999;
      float deltaPhi, deltaEta = -999;
      float dPhiStar = -999;
      float ldeltaEventPlane = -999;
      int trigpTBin, mixpTBin;
      int evp_inout = 0;

//        cout << "pool : " << lnEvent << " " << pool->IsReady() << endl;

      if (pool->IsReady() || pool->GetCurrentNEvents() >= 5) {
        if (fDebug > 4) pool->PrintInfo();

        for (int ievent = 0; ievent < lnEvent; ievent++) {
          lTracksMixing = pool->GetEvent(ievent);
          if (!lTracksMixing) continue;
          lnTracksMixing = lTracksMixing->GetEntriesFast();
          lnTracksTrig = lTracksTrig->GetEntriesFast();
          if (!lnTracksTrig) return 1;
//          cout << "lNEvt : " << lnEvent << " MixingT " << lnTracksMixing << " Trig " << lnTracksTrig << " " << lCentBin << " " << lZVertexBin << endl;

          for (int itrack = 0; itrack < lnTracksTrig; itrack++) {
            if (!(lTracksTrig->At(itrack))) return 1;
            const MixedParticle *lParticleTrig = dynamic_cast<MixedParticle*>(lTracksTrig->At(itrack));
            trigpT = lParticleTrig->Pt();
            if (trigpT < fMinPtTrigCut) continue;
            trigEta = lParticleTrig->Eta();
            if (TMath::Abs(trigEta) >= 0.8) continue;
            trigPhi = lParticleTrig->Phi();
            trigCharge = lParticleTrig->Charge();
            ldeltaEventPlane = trigPhi - eventPlane;
            if (ldeltaEventPlane < 0) ldeltaEventPlane += TMath::TwoPi();
            if (ldeltaEventPlane > TMath::TwoPi()) ldeltaEventPlane -= TMath::TwoPi();
            //        cout << ldeltaEventPlane << "\n";
            trigpTBin = GetPttBin(trigpT);

            if (trigpTBin == -1) continue;
            evp_inout = CheckInOut(ldeltaEventPlane);

            for (int jtrack = 0; jtrack < lnTracksMixing; jtrack++) {
              if (!(lTracksMixing->At(jtrack))) return 1;
              const MixedParticle *lParticleMixing = dynamic_cast<MixedParticle*>(lTracksMixing->At(jtrack));
              //          cout << "PROBLEM1 ?" << endl;
              mixpT = lParticleMixing->Pt();
              if (trigpT < mixpT) continue;
              mixEta = lParticleMixing->Eta();
              if (TMath::Abs(mixEta) >= 0.8) continue;
              mixPhi = lParticleMixing->Phi();
              mixCharge = lParticleMixing->Charge();
              mixpTBin = GetPtaBin(mixpT);

              if (mixpTBin == -1) continue;

              if (fConversionsVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.510e-3, 0.510e-3);

                if (mass < fConversionsVCut * 5)
                {
                  mass = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.510e-3, 0.510e-3);

                  if (mass < fConversionsVCut*fConversionsVCut) 
                    continue;
                }
              }

              // K0s
              if (fResonancesVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.1396);

                const Float_t kK0smass = 0.4976;

                if (TMath::Abs(mass - kK0smass*kK0smass) < fResonancesVCut * 5)
                {
                  mass = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.1396);

                  if (mass > (kK0smass-fResonancesVCut)*(kK0smass-fResonancesVCut) && mass < (kK0smass+fResonancesVCut)*(kK0smass+fResonancesVCut))
                    continue;
                }
              }


              // Lambda
              if (fResonancesVCut > 0 && mixCharge * trigCharge < 0)
              {
                Float_t mass1 = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.9383);
                Float_t mass2 = GetInvMassSquaredCheap(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.9383, 0.1396);

                const Float_t kLambdaMass = 1.115;

                if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
                {
                  mass1 = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.1396, 0.9383);

                  if (mass1 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass1 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                    continue;
                }
                if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < fResonancesVCut * 5)
                {
                  mass2 = GetInvMassSquared(trigpT, trigEta, trigPhi, mixpT, mixEta, mixPhi, 0.9383, 0.1396);

                  if (mass2 > (kLambdaMass-fResonancesVCut)*(kLambdaMass-fResonancesVCut) && mass2 < (kLambdaMass+fResonancesVCut)*(kLambdaMass+fResonancesVCut))
                    continue;
                }
              }


              deltaPhi = trigPhi - mixPhi;
              deltaEta = trigEta - mixEta;

//              cout << mixpTBin << " and " << trigpTBin << "  " << jtrack << "/" << lnTracksMixing << " and " << itrack << "/" << lnTracksTrig  << "cent " << lCentBin " " << lZVertexBin << endl;

              if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
              if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

              dPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, trigCharge, mixCharge, trigpT, mixpT, lbSign);
              if (TMath::Abs(dPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
                continue;

              if (trigpTBin >= lMinimumPtTBinMerging && mixpTBin >= fMinimumPtABinMerging) {
                fHistdEtadPhiMixed[0][lCentBin][lZVertexBin][trigpTBin][fMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
                if (evp_inout > 0) fHistdEtadPhiMixed[evp_inout][lCentBin][lZVertexBin][trigpTBin][fMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

              } else {
                fHistdEtadPhiMixed[0][lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
                if (evp_inout > 0)
                fHistdEtadPhiMixed[evp_inout][lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

              }
              //          cout << "Ho ! " << mixPhi <<"  " <<  mixCharge << "  " << mixEta << "  "<< endl;

            } // jtrack

          } // itrack
        }  // ievent
        //    cout << "problem!" << endl;
        pool->UpdatePool(lTracksTrig);
      } else {
        pool->UpdatePool(lTracksTrig);
        //    cout << lTracksTrig->GetEntriesFast() << endl;
      } // if pool is ready
  }

//  cout << "Domixing return" << endl;


  return 0;
}

void AliAnalysisTaskJetLikeCorrelation::InitHistogramVectors() {

  fHistdEtadPhiSame = vector <vector < vector <vector< vector<TH2D*> > > > > (fNumberOfPlanes, vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPttArray.GetSize() - 1, vector<TH2D*> (fPtaArray.GetSize() - 1) ) ) ) );
  fHistdEtadPhiMixed = vector <vector < vector <vector< vector<TH2D*> > > > > (fNumberOfPlanes,  vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPttArray.GetSize() - 1, vector<TH2D*> (fPtaArray.GetSize() - 1) ) ) ) );

  fHistdEtadPhiSameMCCorrCont = vector <vector < vector <vector< vector<TH2D*> > > > > (fNumberOfPlanes, vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPttArray.GetSize() - 1, vector<TH2D*> (fPtaArray.GetSize() - 1) ) ) ) );
  
  fHistdEtadPhiSameMCCorrPrim = vector <vector < vector <vector< vector<TH2D*> > > > > (fNumberOfPlanes, vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPttArray.GetSize() - 1, vector<TH2D*> (fPtaArray.GetSize() - 1) ) ) ) );

  fPoolSize = vector< vector <int> > (fCentArray.GetSize() - 1, vector< int > (fZVertexArray.GetSize() - 1 ) );
}


void AliAnalysisTaskJetLikeCorrelation::InitHistograms() {

  InitHistogramVectors();
  int lPtBin = 250;
  int lPtMin = 0;
  int lPtMax = 25;
//  int lMinPtBin = GetPtBin(fMinPtTrigCut);

  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    for (int icent = 0; icent < fCentArray.GetSize() - 1; icent++) {
      for (int izvertex = 0; izvertex < fZVertexArray.GetSize() - 1; izvertex++) {
        for (int iptt = 0; iptt < fPttArray.GetSize() - 1; iptt++) {
//          if (iptt < lMinPtBin) continue;
          for (int ipta = 0; ipta < fPtaArray.GetSize() - 1; ipta++) {
            if (fPttArray[iptt] + 0.1 < fPtaArray[ipta] ) continue;
            fHistdEtadPhiSame[iplane][icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSame%02d%02d%02d%02d%02d", iplane, icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSame%02d%02d%02d%02d%02d", iplane, icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*8/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
            fHistdEtadPhiMixed[iplane][icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiMixed%02d%02d%02d%02d%02d",iplane, icent,izvertex, iptt, ipta), Form("fHistdEtadPhiMixed%02d%02d%02d%02d%02d", iplane, icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*8/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);

            fHistdEtadPhiSame[iplane][icent][izvertex][iptt][ipta]->Sumw2();
            fHistdEtadPhiMixed[iplane][icent][izvertex][iptt][ipta]->Sumw2();

            if (fMCCorrection) {
              fHistdEtadPhiSameMCCorrCont[iplane][icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrCont%02d%02d%02d%02d%02d", iplane, icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrCont%02d%02d%02d%02d%02d", iplane, icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*8/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);

              fHistdEtadPhiSameMCCorrCont[iplane][icent][izvertex][iptt][ipta]->Sumw2();

              fHistdEtadPhiSameMCCorrPrim[iplane][icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrPrim%02d%02d%02d%02d%02d", iplane, icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrPrim%02d%02d%02d%02d%02d", iplane, icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*8/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);

              fHistdEtadPhiSameMCCorrPrim[iplane][icent][izvertex][iptt][ipta]->Sumw2();

            }
          } // ipta
        } // iptt
      } //izvertex 
    }
  } // iplane
  
  fHistNevtSame = new TH2D(Form("fHistNevtSame"), Form("fHistNevtSame"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
  fHistNevtSame->GetXaxis()->SetTitle("CentBin");
  fHistNevtSame->GetYaxis()->SetTitle("ZVertBin");

  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    fHistPtSame[iplane] = new TH3D(Form("fHistPtSame%02d", iplane), Form("fHistPtSame%02d", iplane), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5, lPtBin, lPtMin, lPtMax);
    fHistPtSame[iplane]->GetXaxis()->SetTitle("CentBin");
    fHistPtSame[iplane]->GetYaxis()->SetTitle("ZVertBin");
    fHistPtSame[iplane]->GetZaxis()->SetTitle("Pt(GeV/c)");

  }
//  float lEtaMin = -0.9;
//  float lEtaMax = 0.9;

  int nbins[5] = {8, 20, 9, 20, 50};
  // In/out(5), Cent, ZVtx, Eta, Pt // total 5
  double maxbin[5] = {7.5, 100, 9,1,25 };
  double minbin[5] = { -0.5, 0, -9, -1, 0};
  fHistEtaSparse = new THnSparseD("fHistEtaSparse", "fHistEtaSparse", 5, nbins, minbin, maxbin);
  fHistEtaSparse->Sumw2();

}

float AliAnalysisTaskJetLikeCorrelation::CalculatedPhiStar (float dPhi, float dEta, float charge1, float charge2, float pT1, float pT2, float bSign) {
  float dPhiStartmp = 0;
  float dPhiStartmpabs = 0;
  float dPhiStar = 999;
  float dPhiStarabs = 999;
  float Bze = 0.075;
  float radius = 0;

  const float kLimit = fTwoTrackEffCut*3;
  static const double kPi = TMath::Pi();

  if (TMath::Abs(dEta) < fTwoTrackEffCut * 2.5 * 3) {
    float dPhiStarInner = dPhi + charge1 * bSign * TMath::ASin( (Bze * 0.8)/(2*pT1) ) - charge2 * bSign * TMath::ASin( (Bze * 0.8)/(2*pT2) );
    float dPhiStarOuter = dPhi + charge1 * bSign * TMath::ASin( (Bze * 2.5)/(2*pT1) ) - charge2 * bSign * TMath::ASin( (Bze * 2.5)/(2*pT2) );

    if (dPhiStarInner > kPi) dPhiStarInner = kPi * 2 - dPhiStarInner;
    if (dPhiStarInner < -kPi) dPhiStarInner = -kPi * 2 - dPhiStarInner;
    if (dPhiStarInner > kPi) dPhiStarInner = kPi * 2 - dPhiStarInner; // might look funny but is needed

    if (dPhiStarOuter > kPi) dPhiStarOuter = kPi * 2 - dPhiStarOuter;
    if (dPhiStarOuter < -kPi) dPhiStarOuter = -kPi * 2 - dPhiStarOuter;
    if (dPhiStarOuter > kPi) dPhiStarOuter = kPi * 2 - dPhiStarOuter; // might look funny but is needed

    if (TMath::Abs(dPhiStarInner) < kLimit || TMath::Abs(dPhiStarOuter) < kLimit || dPhiStarInner * dPhiStarOuter < 0) {

      //Find the smallest dPhiStar

      for (int ir = 80; ir < 251; ir++) {  // Radius of TPC : 0.8 - 2.5 (m)
        radius = ir*0.01;
        dPhiStartmp = dPhi + charge1 * bSign * TMath::ASin( (Bze * radius)/(2*pT1) ) - charge2 * bSign * TMath::ASin( (Bze * radius) / (2*pT2));
        dPhiStartmpabs = TMath::Abs(dPhiStartmp);

        if (ir == 80) { dPhiStar = dPhiStartmp; 
          dPhiStarabs = dPhiStartmpabs; }
        if (dPhiStartmp < dPhiStar) { dPhiStar = dPhiStartmp;
          dPhiStarabs = dPhiStartmpabs; }
      }
    }

    if (dPhiStar > kPi) dPhiStar = kPi * 2 - dPhiStar;
    if (dPhiStar < -1*kPi) dPhiStar = -2*kPi - dPhiStar;
    if (dPhiStar > kPi) dPhiStar = kPi * 2 - dPhiStar;

  }
  return dPhiStar;

}

int AliAnalysisTaskJetLikeCorrelation::GetCentBin(double cent) {
    for (int icentbin = 0; icentbin < fCentArray.GetSize() - 1 ; icentbin++) {
          if (cent >= fCentArray[icentbin] && cent < fCentArray[icentbin+1])
                  return icentbin;
            }
      return -1;
}

int AliAnalysisTaskJetLikeCorrelation::GetZVertexBin(double zvertex) {
    for (int izvertexbin = 0; izvertexbin < fZVertexArray.GetSize() - 1; izvertexbin++) {
          if (zvertex >= fZVertexArray[izvertexbin] && zvertex  < fZVertexArray[izvertexbin+1])
                  return izvertexbin;
            }
      return -1;
}

int AliAnalysisTaskJetLikeCorrelation::GetPttBin(double pt) {
  for (int iptbin = 0; iptbin < fPttArray.GetSize() - 1 ; iptbin++) {
    if (pt >= fPttArray[iptbin] && pt < fPttArray[iptbin+1])
      return iptbin;
  }
  return -1;
}

int AliAnalysisTaskJetLikeCorrelation::GetPtaBin(double pt) {
  for (int iptbin = 0; iptbin < fPtaArray.GetSize() - 1 ; iptbin++) {
    if (pt >= fPtaArray[iptbin] && pt < fPtaArray[iptbin+1])
      return iptbin;
  }
  return -1;
}

int AliAnalysisTaskJetLikeCorrelation::GetEventInformation(AliAODEvent *aodevent) {

  if (fMCTruth) {
    fZVertex = fMCHeader->GetVtxZ();
    AliGenHijingEventHeader *lMCHeaderHijing = dynamic_cast<AliGenHijingEventHeader*>(fMCEvent->GenEventHeader());
    if (lMCHeaderHijing) {
      double lIP = lMCHeaderHijing->ImpactParameter();
      fCentPercentile = GetCentralityFromIP(lIP);
      if (fCentPercentile - 0.0001 < 0.0001) return 20;

      float rpAngle = -999;
      if (lMCHeaderHijing) {
        rpAngle = lMCHeaderHijing->ReactionPlaneAngle();
        if (rpAngle < 0) rpAngle += TMath::Pi();   // V0 EP angle varies from -pi/2 to +pi/2
        fEventPlane = rpAngle;
        fEventPlaneV0A = rpAngle;
        fEventPlaneV0C = rpAngle;
        fEventPlaneTPC = rpAngle;
      }
    }
    if (fCollision == kpp) fCentPercentile = 1;

    return 0;

  }

  
  const AliVVertex *vtxTrk = aodevent->GetPrimaryVertex();
  const AliVVertex *vtxSPD = aodevent->GetPrimaryVertexSPD();

  if (fPeriod > k15a) { // Additional vertex cuts, for run2
    if (!(vtxTrk) || vtxTrk->GetNContributors() < 1 || !(vtxSPD)) return 10;
    double covTrk[6],covSPD[6];
    vtxTrk->GetCovarianceMatrix(covTrk);
    vtxSPD->GetCovarianceMatrix(covSPD);
    double dz = vtxTrk->GetZ()-vtxSPD->GetZ();
    double errTot = TMath::Sqrt(covTrk[5]+covSPD[5]);
    double errTrk = TMath::Sqrt(covTrk[5]);
    double nsigTot = TMath::Abs(dz)/errTot, nsigTrk = TMath::Abs(dz)/errTrk;
    if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrk>20){
      // reject, bad reconstructed track vertex
      return 10;
    }
  } else {
    if (!vtxTrk || vtxTrk->GetNContributors() < 1) return 10;
  } // Bad ZVertex : error code 10's

  fZVertex = vtxTrk->GetZ();
  
  fCentPercentile = aodevent->GetCentrality()->GetCentralityPercentile("V0M");
  if (fCollision == kPbPb) {
    if (fCentPercentile > 110 || fCentPercentile < 0) return 20;  // Centrality error : 20's


    const AliQnCorrectionsQnVector *myV0QnVector;
    const AliQnCorrectionsQnVector *myV0AQnVector;
    const AliQnCorrectionsQnVector *myV0CQnVector;
    const AliQnCorrectionsQnVector *myTPCQnVector;
    myV0QnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZERO");
    if (myV0QnVector != NULL) {
      fEventPlane = myV0QnVector->EventPlane(2);
    }
    myV0AQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA");
    if (myV0AQnVector != NULL){
      fEventPlaneV0A = myV0AQnVector->EventPlane(2);
    }
    myV0CQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC");
    if (myV0CQnVector != NULL){
      fEventPlaneV0C = myV0CQnVector->EventPlane(2);
    }
    myTPCQnVector = fFlowQnVectorMgr->GetDetectorQnVector("TPC");
    if (myTPCQnVector != NULL) {
      fEventPlaneTPC = myTPCQnVector->EventPlane(2);
    }

    if (fEventPlaneTPC >=0 && fEventPlaneV0A >= 0 && fEventPlaneV0C >= 0) {
      fHistResolutionV2[0]->Fill(fCentPercentile, TMath::Cos(2*(fEventPlaneV0C - fEventPlaneTPC)));
      fHistResolutionV2[1]->Fill(fCentPercentile, TMath::Cos(2*(fEventPlaneV0C - fEventPlaneV0A)));
      fHistResolutionV2[2]->Fill(fCentPercentile, TMath::Cos(2*(fEventPlaneTPC - fEventPlaneV0A)));
    }
  
  }
  else if (fCollision == kpp) {
    fCentPercentile = 1;
    fEventPlane = 0;
    fEventPlaneV0A =0;
    fEventPlaneV0C =0;
    fEventPlaneTPC =0;
  } else {

  } // for pPb's case

  return 0;

}

double AliAnalysisTaskJetLikeCorrelation::GetCentralityFromIP(double ip) {

	//https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies
	static double bmin[12] = {0.0,1.60,2.27,3.72,5.23,7.31,8.88,10.20,11.38,12.47,14.51,100};
	static double centmean[12] = {0.5,1.5,3.5,7.5,15,25,35,45,55,65,75,90};
	for(UInt_t i = 0; i < 11; i++){
		if(bmin[i+1] > ip)
			return centmean[i];
	}
	return 0.0;

}

int AliAnalysisTaskJetLikeCorrelation::FillPt(double deltaevp, int centbin, int zvertbin, double trigpt) {

  if (deltaevp < 1./12*pi && deltaevp >= 0) {
    fHistPtSame[1]->Fill(centbin, zvertbin, trigpt); 
    return 1;
  } 
  else if (deltaevp < 2./12*pi && deltaevp >= 1./12*pi) {
    fHistPtSame[3]->Fill(centbin, zvertbin, trigpt); 
    return 3;
  } 
  else if (deltaevp < 3./12*pi && deltaevp >= 2./12*pi) { 
    fHistPtSame[4]->Fill(centbin, zvertbin, trigpt); 
    return 4;
  } 
  else if (deltaevp < 4./12*pi && deltaevp >= 3./12*pi) { 
    fHistPtSame[5]->Fill(centbin, zvertbin, trigpt); 
    return 5;
  } 
  else if (deltaevp < 5./12*pi && deltaevp >= 4./12*pi) { 
    fHistPtSame[6]->Fill(centbin, zvertbin, trigpt); 
    return 6;
  } 
  else if (deltaevp < 7./12*pi && deltaevp >= 5./12*pi) { 
    fHistPtSame[2]->Fill(centbin, zvertbin, trigpt); 
    return 2;
  } 
  else if (deltaevp < 8./12*pi && deltaevp >= 7./12*pi) { 
    fHistPtSame[6]->Fill(centbin, zvertbin, trigpt); 
    return 6;
  } 
  else if (deltaevp < 9./12*pi && deltaevp >= 8./12*pi) { 
    fHistPtSame[5]->Fill(centbin, zvertbin, trigpt); 
    return 5;
  } 
  else if (deltaevp < 10./12*pi && deltaevp >= 9./12*pi) {
    fHistPtSame[4]->Fill(centbin, zvertbin, trigpt); 
    return 4;
  } 
  else if (deltaevp < 11./12*pi && deltaevp >= 10./12*pi) { 
    fHistPtSame[3]->Fill(centbin, zvertbin, trigpt); 
    return 3;
  } 
  else if (deltaevp < 13./12*pi && deltaevp >= 11./12*pi) { 
    fHistPtSame[1]->Fill(centbin, zvertbin, trigpt); 
    return 1;
  } 
  else if (deltaevp < 14./12*pi && deltaevp >= 13./12*pi) { 
    fHistPtSame[3]->Fill(centbin, zvertbin, trigpt); 
    return 3;
  } 
  else if (deltaevp < 15./12*pi && deltaevp >= 14./12*pi) { 
    fHistPtSame[4]->Fill(centbin, zvertbin, trigpt); 
    return 4;
  } 
  else if (deltaevp < 16./12*pi && deltaevp >= 15./12*pi) { 
    fHistPtSame[5]->Fill(centbin, zvertbin, trigpt); 
    return 5;
  } 
  else if (deltaevp < 17./12*pi && deltaevp >= 16./12*pi) { 
    fHistPtSame[6]->Fill(centbin, zvertbin, trigpt); 
    return 6;
  } 
  else if (deltaevp < 19./12*pi && deltaevp >= 17./12*pi) { 
    fHistPtSame[2]->Fill(centbin, zvertbin, trigpt); 
    return 2;
  } 
  else if (deltaevp < 20./12*pi && deltaevp >= 19./12*pi) { 
    fHistPtSame[6]->Fill(centbin, zvertbin, trigpt); 
    return 6;
  } 
  else if (deltaevp < 21./12*pi && deltaevp >= 20./12*pi) { 
    fHistPtSame[5]->Fill(centbin, zvertbin, trigpt); 
    return 5;
  } 
  else if (deltaevp < 22./12*pi && deltaevp >= 21./12*pi) { 
    fHistPtSame[4]->Fill(centbin, zvertbin, trigpt); 
    return 4;
  } 
  else if (deltaevp < 23./12*pi && deltaevp >= 22./12*pi) { 
    fHistPtSame[3]->Fill(centbin, zvertbin, trigpt); 
    return 3;
  } 
  else if (deltaevp < 24./12*pi && deltaevp >= 23./12*pi) { 
    fHistPtSame[1]->Fill(centbin, zvertbin, trigpt); 
    return 1;
  } 

  return -9;
}

int AliAnalysisTaskJetLikeCorrelation::CheckInOut(double deltaevp) {

  if (fNumberOfPlanes > 2) {
    if (deltaevp < 1./12*pi && deltaevp >= 0) {
      return 1;
    } 
    else if (deltaevp < 2./12*pi && deltaevp >= 1./12*pi) {
      return 3;
    } 
    else if (deltaevp < 3./12*pi && deltaevp >= 2./12*pi) { 
      return 4;
    } 
    else if (deltaevp < 4./12*pi && deltaevp >= 3./12*pi) { 
      return 5;
    } 
    else if (deltaevp < 5./12*pi && deltaevp >= 4./12*pi) { 
      return 6;
    } 
    else if (deltaevp < 7./12*pi && deltaevp >= 5./12*pi) { 
      return 2;
    } 
    else if (deltaevp < 8./12*pi && deltaevp >= 7./12*pi) { 
      return 6;
    } 
    else if (deltaevp < 9./12*pi && deltaevp >= 8./12*pi) { 
      return 5;
    } 
    else if (deltaevp < 10./12*pi && deltaevp >= 9./12*pi) {
      return 4;
    } 
    else if (deltaevp < 11./12*pi && deltaevp >= 10./12*pi) { 
      return 3;
    } 
    else if (deltaevp < 13./12*pi && deltaevp >= 11./12*pi) { 
      return 1;
    } 
    else if (deltaevp < 14./12*pi && deltaevp >= 13./12*pi) { 
      return 3;
    } 
    else if (deltaevp < 15./12*pi && deltaevp >= 14./12*pi) { 
      return 4;
    } 
    else if (deltaevp < 16./12*pi && deltaevp >= 15./12*pi) { 
      return 5;
    } 
    else if (deltaevp < 17./12*pi && deltaevp >= 16./12*pi) { 
      return 6;
    } 
    else if (deltaevp < 19./12*pi && deltaevp >= 17./12*pi) { 
      return 2;
    } 
    else if (deltaevp < 20./12*pi && deltaevp >= 19./12*pi) { 
      return 6;
    } 
    else if (deltaevp < 21./12*pi && deltaevp >= 20./12*pi) { 
      return 5;
    } 
    else if (deltaevp < 22./12*pi && deltaevp >= 21./12*pi) { 
      return 4;
    } 
    else if (deltaevp < 23./12*pi && deltaevp >= 22./12*pi) { 
      return 3;
    } 
    else if (deltaevp < 24./12*pi && deltaevp >= 23./12*pi) { 
      return 1;
    } 
    return -9;
  } 
  else if (fNumberOfPlanes == 2) {
    if (deltaevp < 1./12*pi && deltaevp >= 0) {
      return 1;
    } 
    else if (deltaevp < 7./12*pi && deltaevp >= 5./12*pi) { 
      return 2;
    } 
    else if (deltaevp < 13./12*pi && deltaevp >= 11./12*pi) { 
      return 1;
    } 
    else if (deltaevp < 19./12*pi && deltaevp >= 17./12*pi) { 
      return 2;
    } 
    else if (deltaevp < 24./12*pi && deltaevp >= 23./12*pi) { 
      return 1;
    } 

    return -9;

  }
  else {
    return -9;
  }
}
