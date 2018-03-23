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
  fOutputInc(0), fOutputIn(0), fOutputOut(0), fOutputM1(0), fOutputM2(0),
  fInputHandler(0), fMCHandler(0), fHistZVertex(0), fHistPt(0), fHistPhi(0),fHistCent(0),
  fHistEta(0), fHistPos(0), fHistNeg(0), fHistEventPlaneV0A(0), fHistEventPlaneV0C(0),
  fHistEventPlaneTPC(0), fPoolMgr(0), fPoolMgr_Highpt(0), fInputRoot(0), fInputTree(0),
  fMCCorrection(0), bUseMixingPool(0), fMixingPoolSize(100), fDebug(0), fMinNumTrack(500),
  fEtaCut(0.9), fPhiCut(TMath::TwoPi()), fMinPtTrigCut(3.0), fTwoTrackEffCut(0.02), 
  fFilterBit(128), fTrackDepth(5000), fEventID(0), fIsFirstEvent(1),
  fResonancesVCut(0.02), fConversionsVCut(0.04),
  fHistNevtSame(0), fHistPtSame(0),  fHistPtSameIn(0), fHistPtSameOut(0),
  fHistPtSameM1(0), fHistPtSameM2(0), fHistEtaSparse(0), fHistV2(0),
  fHistContamination(0), fHighPtMixing(0), fTreeStart(0), fTreeEnd(0),
  fCentPercentile(0), fZVertex(0), fEventPlane(0), fEventPlaneV0A(0),
  fEventPlaneV0C(0), fEventPlaneTPC(0), fCollision(kPbPb), fPeriod(k10h),
  flowQnVectorTask(0), fFlowQnVectorMgr(0)
{

}

AliAnalysisTaskJetLikeCorrelation::AliAnalysisTaskJetLikeCorrelation(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputInc(0), fOutputIn(0), fOutputOut(0), fOutputM1(0), fOutputM2(0),
  fInputHandler(0), fMCHandler(0), fHistZVertex(0), fHistPt(0), fHistPhi(0),fHistCent(0),
  fHistEta(0), fHistPos(0), fHistNeg(0), fHistEventPlaneV0A(0), fHistEventPlaneV0C(0),
  fHistEventPlaneTPC(0), fPoolMgr(0), fPoolMgr_Highpt(0), fInputRoot(0), fInputTree(0),
  fMCCorrection(0), bUseMixingPool(0), fMixingPoolSize(100), fDebug(0), fMinNumTrack(500),
  fEtaCut(0.9), fPhiCut(TMath::TwoPi()), fMinPtTrigCut(3.0), fTwoTrackEffCut(0.02), 
  fFilterBit(128), fTrackDepth(5000), fEventID(0), fIsFirstEvent(1),
  fResonancesVCut(0.02), fConversionsVCut(0.04),
  fHistNevtSame(0), fHistPtSame(0),  fHistPtSameIn(0), fHistPtSameOut(0),
  fHistPtSameM1(0), fHistPtSameM2(0), fHistEtaSparse(0), fHistV2(0),
  fHistContamination(0), fHighPtMixing(0), fTreeStart(0), fTreeEnd(0),
  fCentPercentile(0), fZVertex(0), fEventPlane(0), fEventPlaneV0A(0),
  fEventPlaneV0C(0), fEventPlaneTPC(0), fCollision(kPbPb), fPeriod(k10h),
  flowQnVectorTask(0), fFlowQnVectorMgr(0)
{
  DefineInput(0, TChain::Class());
//  DefineInput(1, TList::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

AliAnalysisTaskJetLikeCorrelation::~AliAnalysisTaskJetLikeCorrelation() {

  if (fOutputInc && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputInc;
  }
  if (fOutputIn && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputIn;
  }
  if (fOutputM1 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputM1;
  }
  if (fOutputM2 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputM2;
  }
  if (fOutputOut && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputOut;
  }
}

void AliAnalysisTaskJetLikeCorrelation::Terminate(Option_t *option) {
  fOutputInc = dynamic_cast<TList*>(GetOutputData(1));
  if (NULL == fOutputInc) {
    AliFatal("fOutputIncList == NULL");
    return;
  }
  
  fOutputIn = dynamic_cast<TList*>(GetOutputData(2));
  if (NULL == fOutputIn) {
    AliFatal("fOutputInList == NULL");
    return;
  }
  fOutputOut = dynamic_cast<TList*>(GetOutputData(3));
  if (NULL == fOutputOut) {
    AliFatal("fOutputOutList == NULL");
    return;
  }
  fOutputM1 = dynamic_cast<TList*>(GetOutputData(4));
  if (NULL == fOutputM1) {
    AliFatal("fOutputM1List == NULL");
    return;
  }
  fOutputM2 = dynamic_cast<TList*>(GetOutputData(5));
  if (NULL == fOutputM2) {
    AliFatal("fOutputM2List == NULL");
    return;
  }

}


void AliAnalysisTaskJetLikeCorrelation::UserCreateOutputObjects() {
//  fOutput = new TList();
//  fOutput->SetOwner(kTRUE);
  fIsFirstEvent = 1;
  fOutputInc = new TList();
  fOutputInc->SetOwner(kTRUE);
  fOutputIn = new TList();
  fOutputIn->SetOwner(kTRUE);
  fOutputOut = new TList();
  fOutputOut->SetOwner(kTRUE);
  fOutputM1 = new TList();
  fOutputM1->SetOwner(kTRUE);
  fOutputM2 = new TList();
  fOutputM2->SetOwner(kTRUE);
  
//  fOutput->Add(fOutputInc);
//  fOutput->Add(fOutputIn);
//  fOutput->Add(fOutputOut);
//  fOutput->Add(fOutputM1);
//  fOutput->Add(fOutputM2);

  fHistZVertex = new TH1D("fHistZVertex", "fHistZVertex", 400, -20, 20);
  fHistPt = new TH1D("fHistPt", "fHistPt", 20000, 0, 20);
  fHistPhi = new TH1D("fHistPhi", "fHistPhi", 315*2 , 0, 3.15*2);
  fHistCent = new TH1D("fHistCent", "fHistCent", 100, 0, 100);
  fHistEta = new TH1D("fHistEta", "fHistEta", 200, -1, 1);
  fHistPos = new TH1D("fHistPos", "fHistPos", 20000, 0, 20);
  fHistNeg = new TH1D("fHistNeg", "fHistNeg", 20000, 0, 20);
  fHistEventPlaneV0A = new TH1D("fHistEventPlaneV0A", "fHistEventPlaneV0A", 315 , 0, 3.15);
  fHistEventPlaneV0C = new TH1D("fHistEventPlaneV0C", "fHistEventPlaneV0C", 315 , 0, 3.15);
  fHistEventPlaneTPC = new TH1D("fHistEventPlaneTPC", "fHistEventPlaneTPC", 315 , 0, 3.15);

  fHistV2 = new TProfile("fHistV2", "fHistV2", 10, 0, 100);
  fHistResolutionV2[0] = new TProfile("fHistResolutionV2_01", "fHistResolutionV2_01", 10, 0, 100);
  fHistResolutionV2[1] = new TProfile("fHistResolutionV2_02", "fHistResolutionV2_02", 10, 0, 100);
  fHistResolutionV2[2] = new TProfile("fHistResolutionV2_12", "fHistResolutionV2_12", 10, 0, 100);

  int nbins[6] = {2, 5, 20, 9, 200, 250};
  // Contamination/Primary, In/out(5), Cent, ZVtx, Eta, Pt // total 6
  double maxbin[6] = {1.5, 4.5, 100, 9,1,25 };
  double minbin[6] = {-0.5, -0.5, 0, -9, -1, 0};
  if (fMCCorrection) {
    fHistContamination = new THnSparseD("fHistContamination", "fHistContamination", 6, nbins, minbin, maxbin);
    fHistContamination->Sumw2();
  }
  InitHistograms();

  
  fOutputInc->Add(fHistZVertex);
  fOutputInc->Add(fHistPt);
  fOutputInc->Add(fHistPhi);
  fOutputInc->Add(fHistCent);
  fOutputInc->Add(fHistEta);
  fOutputInc->Add(fHistPos);
  fOutputInc->Add(fHistNeg);
  fOutputInc->Add(fHistEventPlaneV0A);
  fOutputInc->Add(fHistEventPlaneV0C);
  fOutputInc->Add(fHistEventPlaneTPC);
  fOutputInc->Add(fHistNevtSame);
  if (fMCCorrection)
    fOutputInc->Add(fHistContamination);

  fOutputInc->Add(fHistPtSame);
  fOutputIn->Add(fHistPtSameIn);
  fOutputOut->Add(fHistPtSameOut);
  fOutputM1->Add(fHistPtSameM1);
  fOutputM2->Add(fHistPtSameM2);

  int lMinPtBin = GetPtBin(fMinPtTrigCut);
  

  for (int icent = 0; icent < fCentArray.GetSize() - 1; icent++) {
    for (int izvertex = 0; izvertex < fZVertexArray.GetSize() - 1; izvertex++) {
      for (int iptt = 0; iptt < fPtArray.GetSize() - 1; iptt++) {
        if (iptt < lMinPtBin) continue;
        for (int ipta = 0; ipta < fPtArray.GetSize() - 1; ipta++) {
          if (iptt < ipta) continue;

          fOutputInc->Add(fHistdEtadPhiSame[icent][izvertex][iptt][ipta]);
          fOutputInc->Add(fHistdEtadPhiMixed[icent][izvertex][iptt][ipta]);
          if (fMCCorrection) {
            fOutputInc->Add(fHistdEtadPhiSameMCCorrPrim[icent][izvertex][iptt][ipta]);
            fOutputInc->Add(fHistdEtadPhiSameMCCorrCont[icent][izvertex][iptt][ipta]);
          }
        }
      }
    }
  }
  for (int icent = 0; icent < fCentArray.GetSize() - 1; icent++) {
    for (int izvertex = 0; izvertex < fZVertexArray.GetSize() - 1; izvertex++) {
      for (int iptt = 0; iptt < fPtArray.GetSize() - 1; iptt++) {
        if (iptt < lMinPtBin) continue;
        for (int ipta = 0; ipta < fPtArray.GetSize() - 1; ipta++) {
          if (iptt < ipta) continue;

          fOutputIn->Add(fHistdEtadPhiSameIn[icent][izvertex][iptt][ipta]);
          fOutputIn->Add(fHistdEtadPhiMixedIn[icent][izvertex][iptt][ipta]);
          if (fMCCorrection) {
            fOutputIn->Add(fHistdEtadPhiSameMCCorrPrimIn[icent][izvertex][iptt][ipta]);
            fOutputIn->Add(fHistdEtadPhiSameMCCorrContIn[icent][izvertex][iptt][ipta]);
          }

        }
      }
    }
  }
  for (int icent = 0; icent < fCentArray.GetSize() - 1; icent++) {
    for (int izvertex = 0; izvertex < fZVertexArray.GetSize() - 1; izvertex++) {
      for (int iptt = 0; iptt < fPtArray.GetSize() - 1; iptt++) {
        if (iptt < lMinPtBin) continue;
        for (int ipta = 0; ipta < fPtArray.GetSize() - 1; ipta++) {
          if (iptt < ipta) continue;

          fOutputM1->Add(fHistdEtadPhiSameM1[icent][izvertex][iptt][ipta]);
          fOutputM1->Add(fHistdEtadPhiMixedM1[icent][izvertex][iptt][ipta]);

          if (fMCCorrection) {
            fOutputM1->Add(fHistdEtadPhiSameMCCorrPrimM1[icent][izvertex][iptt][ipta]);
            fOutputM1->Add(fHistdEtadPhiSameMCCorrContM1[icent][izvertex][iptt][ipta]);
          }
        }
      }
    }
  }
  for (int icent = 0; icent < fCentArray.GetSize() - 1; icent++) {
    for (int izvertex = 0; izvertex < fZVertexArray.GetSize() - 1; izvertex++) {
      for (int iptt = 0; iptt < fPtArray.GetSize() - 1; iptt++) {
        if (iptt < lMinPtBin) continue;
        for (int ipta = 0; ipta < fPtArray.GetSize() - 1; ipta++) {
          if (iptt < ipta) continue;

          fOutputM2->Add(fHistdEtadPhiSameM2[icent][izvertex][iptt][ipta]);
          fOutputM2->Add(fHistdEtadPhiMixedM2[icent][izvertex][iptt][ipta]);

          if (fMCCorrection) {
            fOutputM2->Add(fHistdEtadPhiSameMCCorrPrimM2[icent][izvertex][iptt][ipta]);
            fOutputM2->Add(fHistdEtadPhiSameMCCorrContM2[icent][izvertex][iptt][ipta]);
          }
        }
      }
    }
  }
  for (int icent = 0; icent < fCentArray.GetSize() - 1; icent++) {
    for (int izvertex = 0; izvertex < fZVertexArray.GetSize() - 1; izvertex++) {
      for (int iptt = 0; iptt < fPtArray.GetSize() - 1; iptt++) {
        if (iptt < lMinPtBin) continue;
        for (int ipta = 0; ipta < fPtArray.GetSize() - 1; ipta++) {
          if (iptt < ipta) continue;

          fOutputOut->Add(fHistdEtadPhiSameOut[icent][izvertex][iptt][ipta]);
          fOutputOut->Add(fHistdEtadPhiMixedOut[icent][izvertex][iptt][ipta]);

          if (fMCCorrection) {
            fOutputOut->Add(fHistdEtadPhiSameMCCorrPrimOut[icent][izvertex][iptt][ipta]);
            fOutputOut->Add(fHistdEtadPhiSameMCCorrContOut[icent][izvertex][iptt][ipta]);
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

  PostData(1, fOutputInc);
  PostData(2, fOutputIn);
  PostData(3, fOutputOut);
  PostData(4, fOutputM1);
  PostData(5, fOutputM2);

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

  AliAODMCHeader *aodmcHeader;

  if (fMCCorrection) 
  {
    fMCEvent = MCEvent();
    aodmcHeader = (AliAODMCHeader*)aodEvent->FindListObject(AliAODMCHeader::StdBranchName());
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
  }

  fEventPlane = -999;
  fEventPlaneV0A = -999;
  fEventPlaneV0C = -999;
  fEventPlaneTPC = -999;

  if (GetEventInformation(aodEvent)) return;

  AliCollisionGeometry *mcCollGeo;
  int nmctracks, label = -999;
  if (fMCCorrection) {
    mcCollGeo= dynamic_cast<AliCollisionGeometry*> (fMCEvent->GenEventHeader());
    float rpAngle = -999;
    if (mcCollGeo) {
      rpAngle = mcCollGeo->ReactionPlaneAngle();
      if (rpAngle < 0) rpAngle += TMath::Pi();   // V0 EP angle varies from -pi/2 to +pi/2
    }
    nmctracks = fMCEvent->GetNumberOfTracks();
  }

  TObjArray *lTrackArray = new TObjArray();
  lTrackArray->SetOwner(kTRUE);
  TObjArray *lTrackArray_HighPt = new TObjArray();
  lTrackArray_HighPt->SetOwner(kTRUE);

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
  
  if (fEntry%1000 == 0) cout << "Entry : " << fEntry << " Centrality : " << fCentPercentile << " ZVertex : " << fZVertex << " Evp : " << fEventPlaneV0A << endl;
  if (lCentBin == -1 || lZVertexBin == -1) return;   //To reduce time for mixing
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
        if (aodTrack->GetType() != AliAODTrack::kPrimary) continue; // Check whether it's primary -> Leads to non flat phi distribution in hybrid tracks!
        if (!aodTrack->TestFilterBit(fFilterBit)) {
          //     std::cout << "filterbit no" << std::endl;
          continue; //TPC Only Track cut : 128 
        }
      }

      lChargeTrig = aodTrack->Charge();

      lpTTrig = aodTrack->Pt();

      ldeltaEventPlane = lPhiTrig - fEventPlaneV0A;

      if (lpTTrig > 0.5 && lpTTrig < 5) {
        fHistV2->Fill(fCentPercentile, TMath::Cos(2*ldeltaEventPlane));
      }

      if (lpTTrig < 0.8) continue;
      lEtaTrig = aodTrack->Eta();
      lPhiTrig = aodTrack->Phi();
      lpTTrigBin = GetPtBin(lpTTrig);
      if (lpTTrigBin == -1) continue;


      if (TMath::Abs(lEtaTrig) >= 0.8) continue;   

      fHistPhi->Fill(lPhiTrig);
      fHistPt->Fill(lpTTrig);
      fHistEta->Fill(lEtaTrig);
      if (lChargeTrig < 0) fHistPos->Fill(lpTTrig);
      if (lChargeTrig > 0) fHistNeg->Fill(lpTTrig);

//      ldeltaEventPlane = lPhiTrig - fEventPlaneV0A;
      if (ldeltaEventPlane < 0) ldeltaEventPlane += TMath::TwoPi();

      if (fMCCorrection) {
        label = TMath::Abs(aodTrack->GetLabel());
        if (label > nmctracks) continue;
        lmcTrack = (AliAODMCParticle*) fMCEvent->GetTrack(label);
        kPrimaryTrig = lmcTrack->IsPhysicalPrimary();
      }

      etavars[3] = lEtaTrig;
      etavars[4] = lpTTrig;

      fHistPtSame->Fill(lCentBin, lZVertexBin, lpTTrig);
      fHistEtaSparse->Fill(etavars);

      if ( (0 <= ldeltaEventPlane) && (ldeltaEventPlane < 1./6*TMath::Pi()) ) {
        fHistPtSameIn->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 1;
        etavars[0] = 1;
        cont_inout = 1;
      } else if ( (5./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 7./6*TMath::Pi()) ) {
        fHistPtSameIn->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 1;
        etavars[0] = 1;
        cont_inout = 1;
      } else if ( (11./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 12./6*TMath::Pi()) ) {
        fHistPtSameIn->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 1;
        etavars[0] = 1;
        cont_inout = 1;
      } else if ( (1./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 2./6*TMath::Pi()) ) {
        fHistPtSameM1->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 3;
        etavars[0] = 3;
        cont_inout = 3;
      } else if ( (7./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 8./6*TMath::Pi()) ) {
        fHistPtSameM1->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 3;
        etavars[0] = 3;
        cont_inout = 3;
      } else if ( (4./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 5./6*TMath::Pi()) ) {
        fHistPtSameM2->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 4;
        etavars[0] = 4;
        cont_inout = 4;
      } else if ( (10./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 11./6*TMath::Pi()) ) {
        fHistPtSameM2->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 4;
        etavars[0] = 4;
        cont_inout = 4;
      } else if ( (2./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 4./6*TMath::Pi()) ) {
        fHistPtSameOut->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 2;
        etavars[0] = 2;
        cont_inout = 2;
      } else if ( (8./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 10./6*TMath::Pi()) ) {
        fHistPtSameOut->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 2;
        etavars[0] = 2;
        cont_inout = 2;
      }

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
      if (lpTTrig < fMinPtTrigCut) { 
        //      lTrackArray->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));
        continue;
      }
      else 
        lTrackArray_HighPt->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));

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
          if (aodTrackAsso->GetType() != AliAODTrack::kPrimary) continue;
          if (!aodTrackAsso->TestFilterBit(fFilterBit)) {
            //     std::cout << "filterbit no" << std::endl;
            continue; //TPC Only Track cut : 128 
          }
        }

        lpTAsso = aodTrackAsso->Pt();
        if (lpTAsso < 0.8 || lpTTrig < lpTAsso) continue;
        lpTAssoBin = GetPtBin(lpTAsso);
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


          fHistdEtadPhiSameMCCorrCont[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          if ( cont_inout == 1 ) fHistdEtadPhiSameMCCorrContIn[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          else if ( cont_inout == 2 ) fHistdEtadPhiSameMCCorrContOut[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          else if ( cont_inout == 3 ) fHistdEtadPhiSameMCCorrContM1[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          else if (cont_inout == 4 ) fHistdEtadPhiSameMCCorrContM2[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          if (kPrimaryAsso && kPrimaryTrig) {
            fHistdEtadPhiSameMCCorrPrim[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            if ( cont_inout == 1 ) fHistdEtadPhiSameMCCorrPrimIn[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            else if ( cont_inout == 2 ) fHistdEtadPhiSameMCCorrPrimOut[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            else if ( cont_inout == 3 ) fHistdEtadPhiSameMCCorrPrimM1[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            else if (cont_inout == 4 ) fHistdEtadPhiSameMCCorrPrimM2[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          }
        }

        fHistdEtadPhiSame[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        if ( cont_inout == 1 ) fHistdEtadPhiSameIn[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        else if ( cont_inout == 2 ) fHistdEtadPhiSameOut[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        else if ( cont_inout == 3 ) fHistdEtadPhiSameM1[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        else if (cont_inout == 4 ) fHistdEtadPhiSameM2[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);

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
        if (aodTrack->GetType() != AliAODTrack::kPrimary) continue; // Check whether it's primary -> Leads to non flat phi distribution in hybrid tracks!
        if (!aodTrack->TestFilterBit(fFilterBit)) {
          //     std::cout << "filterbit no" << std::endl;
          continue; //TPC Only Track cut : 128 
        }
      }

      lChargeTrig = aodTrack->Charge();

      lpTTrig = aodTrack->Pt();
      
      ldeltaEventPlane = lPhiTrig - fEventPlaneV0A;

      if (lpTTrig > 0.5 && lpTTrig < 5) {
        fHistV2->Fill(fCentPercentile, TMath::Cos(2*ldeltaEventPlane));
      }
      if (lpTTrig < 0.8) continue;
      lpTTrigBin = GetPtBin(lpTTrig);
      if (lpTTrigBin == -1) continue;
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

      if (fMCCorrection) {
        label = TMath::Abs(aodTrack->GetLabel());
        if (label > nmctracks) continue;
        lmcTrack = (AliAODMCParticle*) fMCEvent->GetTrack(label);
        kPrimaryTrig = lmcTrack->IsPhysicalPrimary();
      }

      etavars[3] = lEtaTrig;
      etavars[4] = lpTTrig;

      fHistPtSame->Fill(lCentBin, lZVertexBin, lpTTrig);
      fHistEtaSparse->Fill(etavars);

      if ( (0 <= ldeltaEventPlane) && (ldeltaEventPlane < 1./6*TMath::Pi()) ) {
        fHistPtSameIn->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 1;
        etavars[0] = 1;
        cont_inout = 1;
      } else if ( (5./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 7./6*TMath::Pi()) ) {
        fHistPtSameIn->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 1;
        etavars[0] = 1;
        cont_inout = 1;
      } else if ( (11./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 12./6*TMath::Pi()) ) {
        fHistPtSameIn->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 1;
        etavars[0] = 1;
        cont_inout = 1;
      } else if ( (1./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 2./6*TMath::Pi()) ) {
        fHistPtSameM1->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 3;
        etavars[0] = 3;
        cont_inout = 3;
      } else if ( (7./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 8./6*TMath::Pi()) ) {
        fHistPtSameM1->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 3;
        etavars[0] = 3;
        cont_inout = 3;
      } else if ( (4./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 5./6*TMath::Pi()) ) {
        fHistPtSameM2->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 4;
        etavars[0] = 4;
        cont_inout = 4;
      } else if ( (10./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 11./6*TMath::Pi()) ) {
        fHistPtSameM2->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 4;
        etavars[0] = 4;
        cont_inout = 4;
      } else if ( (2./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 4./6*TMath::Pi()) ) {
        fHistPtSameOut->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 2;
        etavars[0] = 2;
        cont_inout = 2;
      } else if ( (8./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 10./6*TMath::Pi()) ) {
        fHistPtSameOut->Fill(lCentBin, lZVertexBin, lpTTrig);
        contvars[1] = 2;
        etavars[0] = 2;
        cont_inout = 2;
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
//      if (lpTTrig < fMinPtTrigCut) { 
////        lTrackArray->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));
//        continue;
//      }
//      else 
        lTrackArray_HighPt->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));
        if (lpTTrig < fMinPtTrigCut) continue;

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
          if (aodTrackAsso->GetType() != AliAODTrack::kPrimary) continue;
          if (!aodTrackAsso->TestFilterBit(fFilterBit)) {
            //     std::cout << "filterbit no" << std::endl;
            continue; //TPC Only Track cut : 128 
          }
        }

        lpTAsso = aodTrackAsso->Pt();
        if (lpTAsso < 0.8 || lpTTrig < lpTAsso) continue;
        lpTAssoBin = GetPtBin(lpTAsso);
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


          fHistdEtadPhiSameMCCorrCont[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          if ( cont_inout == 1 ) fHistdEtadPhiSameMCCorrContIn[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          else if ( cont_inout == 2 ) fHistdEtadPhiSameMCCorrContOut[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          else if ( cont_inout == 3 ) fHistdEtadPhiSameMCCorrContM1[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          else if (cont_inout == 4 ) fHistdEtadPhiSameMCCorrContM2[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          if (kPrimaryAsso && kPrimaryTrig) {
            fHistdEtadPhiSameMCCorrPrim[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            if ( cont_inout == 1 ) fHistdEtadPhiSameMCCorrPrimIn[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            else if ( cont_inout == 2 ) fHistdEtadPhiSameMCCorrPrimOut[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            else if ( cont_inout == 3 ) fHistdEtadPhiSameMCCorrPrimM1[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
            else if (cont_inout == 4 ) fHistdEtadPhiSameMCCorrPrimM2[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          }
        }

        fHistdEtadPhiSame[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        if ( cont_inout == 1 ) fHistdEtadPhiSameIn[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        else if ( cont_inout == 2 ) fHistdEtadPhiSameOut[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        else if ( cont_inout == 3 ) fHistdEtadPhiSameM1[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
        else if (cont_inout == 4 ) fHistdEtadPhiSameM2[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);

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
        if (aodTrack->GetType() != AliAODTrack::kPrimary) continue; // Check whether it's primary -> Leads to non flat phi distribution in hybrid tracks!
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
      lpTTrigBin = GetPtBin(lpTTrig);
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
      
      etavars[3] = lEtaTrig;
      etavars[4] = lpTTrig;

      fHistPtSame->Fill(lCentBin, lZVertexBin, lpTTrig);
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

//      if (lpTTrig < fMinPtTrigCut) { 
//        //      lTrackArray->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));
//        continue;
//      }
//      else 
        lTrackArray_HighPt->Add(new MixedParticle(lpTTrig, lEtaTrig, lPhiTrig, lChargeTrig));

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
          if (aodTrackAsso->GetType() != AliAODTrack::kPrimary) continue;
          if (!aodTrackAsso->TestFilterBit(fFilterBit)) {
            //     std::cout << "filterbit no" << std::endl;
            continue; //TPC Only Track cut : 128 
          }
        }

        lpTAsso = aodTrackAsso->Pt();
        if (lpTAsso < 0.8 || lpTTrig < lpTAsso) continue;
        lpTAssoBin = GetPtBin(lpTAsso);
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


          fHistdEtadPhiSameMCCorrCont[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          if (kPrimaryAsso && kPrimaryTrig) {
            fHistdEtadPhiSameMCCorrPrim[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
          }
        }

        fHistdEtadPhiSame[lCentBin][lZVertexBin][lpTTrigBin][lpTAssoBin]->Fill(deltaEta, deltaPhi);
      }

    }
  }  // fCollision = pp


  DoMixing(fCentPercentile, fZVertex, lTrackArray_HighPt, lTrackArray, fEventPlaneV0A, lBSign);
  fHistNevtSame->Fill(lCentBin, lZVertexBin);

  PostData(1, fOutputInc);
  PostData(2, fOutputIn);
  PostData(3, fOutputOut);
  PostData(4, fOutputM1);
  PostData(5, fOutputM2);
}

int AliAnalysisTaskJetLikeCorrelation::SetupMixing() {
  double *centralityarr = fCentArray.GetArray();
  double *zVertexarr = fZVertexArray.GetArray();
  int nCentBin = fCentArray.GetSize();
  int nZVertBin = fZVertexArray.GetSize();

  if (fCollision == kPbPb) {
    int poolsize = 100;
    fPoolMgr = new AliEventPoolManager(poolsize, fTrackDepth, nCentBin, centralityarr, nZVertBin, zVertexarr);
    //    fPoolMgr_Highpt = new AliEventPoolManager(1000, fTrackDepth, nCentBin, centralityarr, nZVertBin, zVertexarr);
  } else if (fCollision == kpPb) {
    int poolsize = 150;
    fPoolMgr = new AliEventPoolManager(poolsize, fTrackDepth, nCentBin, centralityarr, nZVertBin, zVertexarr);
  } else if ( fCollision == kpp) {
    int poolsize = 300;
    fPoolMgr = new AliEventPoolManager(poolsize, fTrackDepth, nCentBin, centralityarr, nZVertBin, zVertexarr);
  } else {
    AliError("[SetupMixing] Incorrect Collision Type! Aborting..");
    return 1;
  }
  return 0;

}

int AliAnalysisTaskJetLikeCorrelation::DoMixing(float lCentrality, float fZVertex, TObjArray *lTracksTrig, TObjArray *lTracks_LowPt, float eventPlane, int lbSign) {

  //  cout << "Starting Do Mixing.. " << endl;
  int lMinimumPtABinMerging = GetPtBin(3.00);
  int lMinimumPtTBinMerging = GetPtBin(6.00);

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
          trigpTBin = GetPtBin(trigpT);
          if (trigpTBin == -1) continue;

          ldeltaEventPlane = trigPhi - eventPlane;
          if (ldeltaEventPlane < 0) ldeltaEventPlane += TMath::TwoPi();

          if ( (0 <= ldeltaEventPlane) && (ldeltaEventPlane < 1./6*TMath::Pi()) ) {
            evp_inout = 1;
          } else if ( (5./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 7./6*TMath::Pi()) ) {
            evp_inout = 1;
          } else if ( (11./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 12./6*TMath::Pi()) ) {
            evp_inout = 1;
          } else if ( (1./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 2./6*TMath::Pi()) ) {
            evp_inout = 3;
          } else if ( (7./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 8./6*TMath::Pi()) ) {
            evp_inout = 3;
          } else if ( (4./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 5./6*TMath::Pi()) ) {
            evp_inout = 4;
          } else if ( (10./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 11./6*TMath::Pi()) ) {
            evp_inout = 4;
          } else if ( (2./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 4./6*TMath::Pi()) ) {
            evp_inout = 2;
          } else if ( (8./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 10./6*TMath::Pi()) ) {
            evp_inout = 2;
          }

          for (int jtrack = 0; jtrack < lnTracksMixing; jtrack++) {
            mixFilterBit = ltracksmixing[jtrack].GetFilterBit();
            if (!(mixFilterBit & fFilterBit)) continue;       // Test filter bit
            if (mixFilterBit == 128 && !(ltracksmixing[jtrack].IsPrimary())) continue;     // Only for TPC only track primary check
            mixpT = ltracksmixing[jtrack].GetPt();
            if (mixpT < 3.0 ||  trigpT < mixpT) continue;  // for pbpb
            //          if (trigpT < mixpT) continue;  // for pp
            mixEta = ltracksmixing[jtrack].GetEta();
            if (TMath::Abs(mixEta) >= 0.8) continue;
            mixPhi = ltracksmixing[jtrack].GetPhi();
            mixCharge = ltracksmixing[jtrack].GetCharge();
            mixpTBin = GetPtBin(mixpT);

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

            if (trigpTBin >= lMinimumPtTBinMerging && mixpTBin >= lMinimumPtABinMerging) {
              fHistdEtadPhiMixed[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
              if ( evp_inout == 1 ) fHistdEtadPhiMixedIn[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
              else if ( evp_inout == 2 ) fHistdEtadPhiMixedOut[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
              else if ( evp_inout == 3 ) fHistdEtadPhiMixedM1[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
              else if ( evp_inout == 4) fHistdEtadPhiMixedM2[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

            } else {
              fHistdEtadPhiMixed[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

              if ( evp_inout == 1 ) fHistdEtadPhiMixedIn[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
              else if ( evp_inout == 2 ) fHistdEtadPhiMixedOut[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
              else if ( evp_inout == 3 ) fHistdEtadPhiMixedM1[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
              else if ( evp_inout == 4) fHistdEtadPhiMixedM2[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
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
          trigpTBin = GetPtBin(trigpT);
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
            mixpTBin = GetPtBin(mixpT);

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

            if (trigpTBin >= lMinimumPtTBinMerging && mixpTBin >= lMinimumPtABinMerging) {
                fHistdEtadPhiMixed[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
            } else { 
                fHistdEtadPhiMixed[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
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
        //        cout << ldeltaEventPlane << "\n";
        trigpTBin = GetPtBin(trigpT);
        
        if (trigpTBin == -1) continue;
        if ( (0 <= ldeltaEventPlane) && (ldeltaEventPlane < 1./6*TMath::Pi()) ) {
          evp_inout = 1;
        } else if ( (5./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 7./6*TMath::Pi()) ) {
          evp_inout = 1;
        } else if ( (11./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 12./6*TMath::Pi()) ) {
          evp_inout = 1;
        } else if ( (1./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 2./6*TMath::Pi()) ) {
          evp_inout = 3;
        } else if ( (7./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 8./6*TMath::Pi()) ) {
          evp_inout = 3;
        } else if ( (4./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 5./6*TMath::Pi()) ) {
          evp_inout = 4;
        } else if ( (10./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 11./6*TMath::Pi()) ) {
          evp_inout = 4;
        } else if ( (2./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 4./6*TMath::Pi()) ) {
          evp_inout = 2;
        } else if ( (8./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 10./6*TMath::Pi()) ) {
          evp_inout = 2;
        }

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
          mixpTBin = GetPtBin(mixpT);

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

//          cout << mixpTBin << " and " << trigpTBin << "  " << jtrack << "/" << lnTracksMixing << " and " << itrack << "/" << lnTracksTrig  <<  endl;

          if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
          if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

          dPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, trigCharge, mixCharge, trigpT, mixpT, lbSign);
          if (TMath::Abs(dPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
            continue;

          if (trigpTBin >= lMinimumPtTBinMerging && mixpTBin >= lMinimumPtABinMerging) {
            fHistdEtadPhiMixed[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
            if ( evp_inout == 1 ) fHistdEtadPhiMixedIn[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
            else if ( evp_inout == 2 ) fHistdEtadPhiMixedOut[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
            else if ( evp_inout == 3 ) fHistdEtadPhiMixedM1[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
            else if ( evp_inout == 4) fHistdEtadPhiMixedM2[lCentBin][lZVertexBin][trigpTBin][lMinimumPtABinMerging]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

          } else {
            fHistdEtadPhiMixed[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

            if ( evp_inout == 1 ) fHistdEtadPhiMixedIn[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
            else if ( evp_inout == 2 ) fHistdEtadPhiMixedOut[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
            else if ( evp_inout == 3 ) fHistdEtadPhiMixedM1[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
            else if ( evp_inout == 4) fHistdEtadPhiMixedM2[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
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

  /*
  AliEventPool *pool_lowpt = fPoolMgr->GetEventPool(lCentrality, fZVertex);
  if (!pool_lowpt) return 1;

  lnEvent = pool_lowpt->GetCurrentNEvents();
//  cout << "pool_lowpt : " << lnEvent << endl;

  if (pool_lowpt->IsReady() || pool_lowpt->NTracksInPool() > fMinNumTrack || pool_lowpt->GetCurrentNEvents() >= 5) {
    if (fDebug > 4) pool_lowpt->PrintInfo();

    lCentBin = GetCentBin(lCentrality);
    lZVertexBin = GetZVertexBin(fZVertex);
    lnTracksMixing = 0;
    lnTracksTrig = 0;

     evp_inout = 0;

    for (int ievent = 0; ievent < lnEvent; ievent++) {
      lTracksMixing = pool_lowpt->GetEvent(ievent);
//      cout << "lNEvt : " << lnEvent << endl;
      if (!lTracksMixing) continue;
      lnTracksMixing = lTracksMixing->GetEntriesFast();
      lnTracksTrig = lTracksTrig->GetEntriesFast();

//      cout << "lNEvt : " << lnEvent << " MixingT " << lnTracksMixing << " Trig " << lnTracksTrig << endl;

      for (int itrack = 0; itrack < lnTracksTrig; itrack++) {
        if (!(lTracksTrig->At(itrack))) return 1;
        const MixedParticle *lParticleTrig = dynamic_cast<MixedParticle*>(lTracksTrig->At(itrack));
        trigpT = lParticleTrig->Pt();
        trigEta = lParticleTrig->Eta();
        if (TMath::Abs(trigEta) >= 0.8) continue;
        trigPhi = lParticleTrig->Phi();
        trigCharge = lParticleTrig->Charge();
        ldeltaEventPlane = trigPhi - fEventPlane;
        if (ldeltaEventPlane < 0) ldeltaEventPlane += TMath::TwoPi();
        //        cout << ldeltaEventPlane << "\n";
        trigpTBin = GetPtBin(trigpT);
        
        if (trigpTBin == -1) continue;
        if ( (0 <= ldeltaEventPlane) && (ldeltaEventPlane < 1./6*TMath::Pi()) ) {
          evp_inout = 1;
        } else if ( (5./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 7./6*TMath::Pi()) ) {
          evp_inout = 1;
        } else if ( (11./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 12./6*TMath::Pi()) ) {
          evp_inout = 1;
        } else if ( (1./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 2./6*TMath::Pi()) ) {
          evp_inout = 3;
        } else if ( (7./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 8./6*TMath::Pi()) ) {
          evp_inout = 3;
        } else if ( (4./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 5./6*TMath::Pi()) ) {
          evp_inout = 4;
        } else if ( (10./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 11./6*TMath::Pi()) ) {
          evp_inout = 4;
        } else if ( (2./6*TMath::Pi()<= ldeltaEventPlane) && (ldeltaEventPlane < 4./6*TMath::Pi()) ) {
          evp_inout = 2;
        } else if ( (8./6*TMath::Pi() <= ldeltaEventPlane) && (ldeltaEventPlane < 10./6*TMath::Pi()) ) {
          evp_inout = 2;
        }

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
          mixpTBin = GetPtBin(mixpT);

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

//          cout << mixpTBin << " and " << trigpTBin << "  " << jtrack << "/" << lnTracksMixing << " and " << itrack << "/" << lnTracksTrig  <<  endl;

          if (deltaPhi < -0.5*TMath::Pi()) deltaPhi += TMath::TwoPi();
          if (deltaPhi > 1.5*TMath::Pi()) deltaPhi -= TMath::TwoPi();

          dPhiStar = CalculatedPhiStar (deltaPhi, deltaEta, trigCharge, mixCharge, trigpT, mixpT, lbSign);
          if (TMath::Abs(dPhiStar) < 0.02 && TMath::Abs(deltaEta) < 0.02) 
            continue;

          fHistdEtadPhiMixed[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
//          cout << fHistdEtadPhiMixed[lCentBin][lZVertexBin][trigpTBin][mixpTBin] << " & " << fHistdEtadPhiMixed[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->GetEntries() << "\n";

          if ( evp_inout == 1 ) fHistdEtadPhiMixedIn[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
          else if ( evp_inout == 2 ) fHistdEtadPhiMixedOut[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
          else if ( evp_inout == 3 ) fHistdEtadPhiMixedM1[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);
          else if ( evp_inout == 4) fHistdEtadPhiMixedM2[lCentBin][lZVertexBin][trigpTBin][mixpTBin]->Fill(deltaEta, deltaPhi, 1.0/lnEvent);

//          cout << "Ho ! " << mixPhi <<"  " <<  mixCharge << "  " << mixEta << "  "<< endl;

        }

      }
    }
//    cout << "problem!" << endl;
    pool_lowpt->UpdatePool(lTracks_LowPt);

  } else 
    pool_lowpt->UpdatePool(lTracks_LowPt);
*/       // Low PT part!

  }
  return 0;
}

void AliAnalysisTaskJetLikeCorrelation::InitHistogramVectors() {

  fHistdEtadPhiSame = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiMixed = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameIn = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiMixedIn = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameM1 = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiMixedM1 = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameM2 = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiMixedM2 = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameOut = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiMixedOut = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );

  fHistdEtadPhiSameMCCorrCont = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameMCCorrContIn = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameMCCorrContM1 = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameMCCorrContM2 = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameMCCorrContOut = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  
  fHistdEtadPhiSameMCCorrPrim = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameMCCorrPrimIn = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameMCCorrPrimM1 = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameMCCorrPrimM2 = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fHistdEtadPhiSameMCCorrPrimOut = vector< vector< vector< vector<TH2D*> > > > (fCentArray.GetSize() - 1, vector< vector< vector<TH2D*> > > (fZVertexArray.GetSize() - 1, vector< vector<TH2D*> > (fPtArray.GetSize() - 1, vector<TH2D*> (fPtArray.GetSize() - 1) ) ) );
  fPoolSize = vector< vector <int> > (fCentArray.GetSize() - 1, vector< int > (fZVertexArray.GetSize() - 1 ) );
}


void AliAnalysisTaskJetLikeCorrelation::InitHistograms() {

  InitHistogramVectors();
  int lPtBin = 25000;
  int lPtMin = 0;
  int lPtMax = 25;
  int lMinPtBin = GetPtBin(fMinPtTrigCut);

  for (int icent = 0; icent < fCentArray.GetSize() - 1; icent++) {
    for (int izvertex = 0; izvertex < fZVertexArray.GetSize() - 1; izvertex++) {
      for (int iptt = 0; iptt < fPtArray.GetSize() - 1; iptt++) {
        if (iptt < lMinPtBin) continue;
        for (int ipta = 0; ipta < fPtArray.GetSize() - 1; ipta++) {
          if (iptt < ipta) continue;
          fHistdEtadPhiSame[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSame%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSame%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiMixed[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiMixed%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiMixed%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameIn[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameIn%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameIn%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiMixedIn[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiMixedIn%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiMixedIn%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameM1[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameM1%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameM1%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiMixedM1[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiMixedM1%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiMixedM1%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameM2[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameM2%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameM2%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiMixedM2[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiMixedM2%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiMixedM2%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameOut[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameOut%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameOut%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiMixedOut[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiMixedOut%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiMixedOut%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);

          fHistdEtadPhiSame[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameIn[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameOut[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameM1[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameM2[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiMixed[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiMixedIn[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiMixedOut[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiMixedM1[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiMixedM2[icent][izvertex][iptt][ipta]->Sumw2();
          
          if (fMCCorrection) {
          fHistdEtadPhiSameMCCorrCont[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrCont%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrCont%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameMCCorrContIn[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrContIn%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrContIn%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameMCCorrContM1[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrContM1%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrContM1%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameMCCorrContM2[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrContM2%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrContM2%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameMCCorrContOut[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrContOut%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrContOut%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);

          fHistdEtadPhiSameMCCorrCont[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameMCCorrContIn[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameMCCorrContOut[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameMCCorrContM1[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameMCCorrContM2[icent][izvertex][iptt][ipta]->Sumw2();
          
          fHistdEtadPhiSameMCCorrPrim[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrPrim%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrPrim%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameMCCorrPrimIn[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrPrimIn%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrPrimIn%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameMCCorrPrimM1[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrPrimM1%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrPrimM1%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameMCCorrPrimM2[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrPrimM2%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrPrimM2%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);
          fHistdEtadPhiSameMCCorrPrimOut[icent][izvertex][iptt][ipta] = new TH2D(Form("fHistdEtadPhiSameMCCorrPrimOut%02d%02d%02d%02d", icent,izvertex, iptt, ipta), Form("fHistdEtadPhiSameMCCorrPrimOut%02d%02d%02d%02d", icent,izvertex, iptt, ipta), int((fEtaCut+0.001)*4/0.1), -2*fEtaCut, 2*fEtaCut, 64, -0.25*fPhiCut, 0.75*fPhiCut);

          fHistdEtadPhiSameMCCorrPrim[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameMCCorrPrimIn[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameMCCorrPrimOut[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameMCCorrPrimM1[icent][izvertex][iptt][ipta]->Sumw2();
          fHistdEtadPhiSameMCCorrPrimM2[icent][izvertex][iptt][ipta]->Sumw2();

          }
        } // ipta
      } // iptt
    } //izvertex 
  }
  
//  fHistNevtMixed = new TH2D(Form("fHistNevtMixed"), Form("fHistNevtMixed"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
//  fHistNevtMixed->GetXaxis()->SetTitle("CentBin");
//  fHistNevtMixed->GetYaxis()->SetTitle("ZVertBin");
  fHistNevtSame = new TH2D(Form("fHistNevtSame"), Form("fHistNevtSame"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
  fHistNevtSame->GetXaxis()->SetTitle("CentBin");
  fHistNevtSame->GetYaxis()->SetTitle("ZVertBin");
  
//  fHistNevtMixedIn = new TH2D(Form("fHistNevtMixedIn"), Form("fHistNevtMixedIn"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
//  fHistNevtMixedIn->GetXaxis()->SetTitle("CentBin");
//  fHistNevtMixedIn->GetYaxis()->SetTitle("ZVertBin");
//  fHistNevtSameIn = new TH2D(Form("fHistNevtSameIn"), Form("fHistNevtSameIn"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
//  fHistNevtSameIn->GetXaxis()->SetTitle("CentBin");
//  fHistNevtSameIn->GetYaxis()->SetTitle("ZVertBin");
  //fHistNevtMixedM1 = new TH2D(Form("fHistNevtMixedM1"), Form("fHistNevtMixedM1"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
  //fHistNevtMixedM1->GetXaxis()->SetTitle("CentBin");
  //fHistNevtMixedM1->GetYaxis()->SetTitle("ZVertBin");
//fHistNevtSameM1 = new TH2D(Form("fHistNevtSameM1"), Form("fHistNevtSameM1"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
//fHistNevtSameM1->GetXaxis()->SetTitle("CentBin");
//fHistNevtSameM1->GetYaxis()->SetTitle("ZVertBin");
  //fHistNevtMixedM2 = new TH2D(Form("fHistNevtMixedM2"), Form("fHistNevtMixedM2"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
  //fHistNevtMixedM2->GetXaxis()->SetTitle("CentBin");
  //fHistNevtMixedM2->GetYaxis()->SetTitle("ZVertBin");
//fHistNevtSameM2 = new TH2D(Form("fHistNevtSameM2"), Form("fHistNevtSameM2"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
//fHistNevtSameM2->GetXaxis()->SetTitle("CentBin");
//fHistNevtSameM2->GetYaxis()->SetTitle("ZVertBin");
  //fHistNevtMixedOut = new TH2D(Form("fHistNevtMixedOut"), Form("fHistNevtMixedOut"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
  //fHistNevtMixedOut->GetXaxis()->SetTitle("CentBin");
  //fHistNevtMixedOut->GetYaxis()->SetTitle("ZVertBin");
//fHistNevtSameOut = new TH2D(Form("fHistNevtSameOut"), Form("fHistNevtSameOut"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
//fHistNevtSameOut->GetXaxis()->SetTitle("CentBin");
//fHistNevtSameOut->GetYaxis()->SetTitle("ZVertBin");

  /*
  fHistPtMixed = new TH2I(Form("fHistPtMixed"), Form("fHistPtMixed"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5);
  fHistPtMixed->GetXaxis()->SetTitle("CentBin");
  fHistPtMixed->GetYaxis()->SetTitle("ZVertBin"); */
  fHistPtSame = new TH3D(Form("fHistPtSame"), Form("fHistPtSame"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5, lPtBin, lPtMin, lPtMax);
  fHistPtSame->GetXaxis()->SetTitle("CentBin");
  fHistPtSame->GetYaxis()->SetTitle("ZVertBin");
  fHistPtSame->GetZaxis()->SetTitle("Pt(GeV/c)");
  
  
  fHistPtSameIn = new TH3D(Form("fHistPtSameIn"), Form("fHistPtSameIn"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5, lPtBin, lPtMin, lPtMax);
  fHistPtSameIn->GetXaxis()->SetTitle("CentBin");
  fHistPtSameIn->GetYaxis()->SetTitle("ZVertBin");
  fHistPtSameIn->GetZaxis()->SetTitle("Pt(GeV/c)");
  fHistPtSameM1 = new TH3D(Form("fHistPtSameM1"), Form("fHistPtSameM1"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5, lPtBin, lPtMin, lPtMax);
  fHistPtSameM1->GetXaxis()->SetTitle("CentBin");
  fHistPtSameM1->GetYaxis()->SetTitle("ZVertBin");
  fHistPtSameM1->GetZaxis()->SetTitle("Pt(GeV/c)");
  fHistPtSameM2 = new TH3D(Form("fHistPtSameM2"), Form("fHistPtSameM2"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5, lPtBin, lPtMin, lPtMax);
  fHistPtSameM2->GetXaxis()->SetTitle("CentBin");
  fHistPtSameM2->GetYaxis()->SetTitle("ZVertBin");
  fHistPtSameM2->GetZaxis()->SetTitle("Pt(GeV/c)");
  fHistPtSameOut = new TH3D(Form("fHistPtSameOut"), Form("fHistPtSameOut"), fCentArray.GetSize() - 1, -0.5, fCentArray.GetSize() - 1.5, fZVertexArray.GetSize() - 1, -0.5, fZVertexArray.GetSize() - 1.5, lPtBin, lPtMin, lPtMax);
  fHistPtSameOut->GetXaxis()->SetTitle("CentBin");
  fHistPtSameOut->GetYaxis()->SetTitle("ZVertBin");
  fHistPtSameOut->GetZaxis()->SetTitle("Pt(GeV/c)");

//  float lEtaMin = -0.9;
//  float lEtaMax = 0.9;

  int nbins[5] = {5, 20, 9, 200, 250};
  // In/out(5), Cent, ZVtx, Eta, Pt // total 5
  double maxbin[5] = {4.5, 100, 9,1,25 };
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

int AliAnalysisTaskJetLikeCorrelation::GetPtBin(double pt) {
    for (int iptbin = 0; iptbin < fPtArray.GetSize() - 1 ; iptbin++) {
          if (pt >= fPtArray[iptbin] && pt < fPtArray[iptbin+1])
                  return iptbin;
            }
      return -1;
}

int AliAnalysisTaskJetLikeCorrelation::GetEventInformation(AliAODEvent *aodevent) {

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
 //   cout << fFlowQnVectorMgr << endl;
    myV0QnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZERO");
    if (myV0QnVector != NULL) {
      fEventPlane = myV0QnVector->EventPlane(2);
      if (fEventPlane < 0) fEventPlane += TMath::Pi();
    }
    myV0AQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA");
    if (myV0AQnVector != NULL){
      fEventPlaneV0A = myV0AQnVector->EventPlane(2);
      if (fEventPlaneV0A < 0) fEventPlaneV0A += TMath::Pi();
    }
    myV0CQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC");
    if (myV0CQnVector != NULL){
      fEventPlaneV0C = myV0CQnVector->EventPlane(2);
      if (fEventPlaneV0C < 0) fEventPlaneV0C += TMath::Pi();
    }
    myTPCQnVector = fFlowQnVectorMgr->GetDetectorQnVector("TPC");
    if (myTPCQnVector != NULL) {
      fEventPlaneTPC = myTPCQnVector->EventPlane(2);
      if (fEventPlaneTPC < 0) fEventPlaneTPC += TMath::Pi();
    }

    if (fEventPlaneTPC >=0 && fEventPlaneV0A >= 0 && fEventPlaneV0C >= 0) {
      fHistResolutionV2[0]->Fill(fCentPercentile, TMath::Cos(2*(fEventPlaneV0A - fEventPlaneTPC)));
      fHistResolutionV2[1]->Fill(fCentPercentile, TMath::Cos(2*(fEventPlaneV0A - fEventPlaneV0C)));
      fHistResolutionV2[2]->Fill(fCentPercentile, TMath::Cos(2*(fEventPlaneTPC - fEventPlaneV0C)));
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
