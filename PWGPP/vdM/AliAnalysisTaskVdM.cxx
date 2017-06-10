#include <memory>

#include <RConfig.h> // SafeDelete
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TString.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliAnalysisManager.h"
#include "AliTriggerIR.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliTriggerIR.h"

#include "AliAnalysisTaskVdM.h"
#include "AliRawEventHeaderBase.h"
#include "AliVVZERO.h"
#include "AliVAD.h"

ClassImp(AliAnalysisTaskVdM);
ClassImp(AliAnalysisTaskVdM::TreeData);

// code from RS for obtaining an unconstrained vertex
Bool_t revertex(AliESDEvent* esdEv, Bool_t kVtxConstr=kTRUE, Int_t algo=6, const Double_t* cuts=0, Bool_t useSA=kTRUE);
Bool_t revertex(AliESDEvent* esdEv, Bool_t kVtxConstr, Int_t algo, const Double_t *cuts, Bool_t useSA)
{
  // Refit ESD VertexTracks and redo tracks->RelateToVertex
  // Default vertexin algorithm is 6 (multivertexer). To use old vertexed, use algo=1
  //
  static AliVertexerTracks* vtFinder = 0;
  static int currRun = 0;
  static int defAlgo = -1;
  static double bkgauss = 0;
  //
  if (!vtFinder) { // create vertexer
    vtFinder = new AliVertexerTracks(esdEv->GetMagneticField());
    printf("Initialized vertexer for VertexTracks refit with field %f kG\n",esdEv->GetMagneticField());
    //
    vtFinder->SetITSMode();
    vtFinder->SetConstraintOff();
  }
  //
  vtFinder->SetITSpureSA(useSA);
  //
  if ( (cuts && algo>11) || algo!=defAlgo) {
    // if cuts array is provided, then interpret algo as the number of parameters in the cuts.
    // otherwise, interpret it as an algorithm ID for hardwired cuts below
    if (cuts) {
      vtFinder->SetCuts((double*)cuts,algo);
      defAlgo = (Int_t)(cuts[10]);
    }
    else {
      const int kNCuts=21;
      double vtCuts[kNCuts] =
        {1.00e-01,1.00e-01,5.00e-01,3.00e+00,1.00e+00,3.00e+00,1.00e+02,
         1.00e+03,3.00e+00,3.00e+01,6.00e+00,4.00e+00,7.00e+00,1.00e+03,
         5.00e+00,5.00e-02,1.00e-03,2.00e+00,1.00e+01,1.00e+00,5.00e+01};
      //
      vtCuts[10] = algo;
      defAlgo = algo;
      vtFinder->SetCuts(vtCuts,kNCuts);
      printf("Setting vertexing algorithm to %d\n",defAlgo);
    }
  }
  if (defAlgo<0 || defAlgo>AliVertexerTracks::kMultiVertexer) {
    printf("Vertexer algorithms 0:%d are supported... \n",defAlgo);
    return kFALSE;
  }
  //
  if (currRun!=esdEv->GetRunNumber() && kVtxConstr) { // update diamond for this run
    double pos[3]={esdEv->GetDiamondX(),esdEv->GetDiamondY(),0};
    Float_t diamondcovxy[3]={0};
    esdEv->GetDiamondCovXY(diamondcovxy);
    Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,7.*7.};
    AliESDVertex initVertex(pos,cov,1.,1);
    vtFinder->SetVtxStart(&initVertex);
    bkgauss = esdEv->GetMagneticField();
    vtFinder->SetFieldkG(bkgauss);
    currRun = esdEv->GetRunNumber();
    printf("Imposed Vtx constraint for run %d\n",currRun);
    initVertex.Print();
  }
  //
  // reset old vertex info
  if (esdEv->GetPileupVerticesTracks()) esdEv->GetPileupVerticesTracks()->Clear();
  ((AliESDVertex*)esdEv->GetPrimaryVertexTracks())->SetNContributors(-1);
  //
  AliESDVertex *pvtx=vtFinder->FindPrimaryVertex(esdEv);
  if (pvtx) {
    if (pvtx->GetStatus()) {
      esdEv->SetPrimaryVertexTracks(pvtx);
      for (Int_t i=esdEv->GetNumberOfTracks(); i--;) {
        AliESDtrack *t = esdEv->GetTrack(i);
        Double_t x[3]; t->GetXYZ(x);
        t->RelateToVertex(pvtx, bkgauss, kVeryBig);
      }
    }
    delete pvtx;
  }
  else return kFALSE;
  //
  return kTRUE;
}

void AliAnalysisTaskVdM::EventInfo::Fill(const AliVEvent* vEvent) {
  const AliVHeader *vHeader = vEvent->GetHeader();
  if (!vHeader) // this is already dealt with in UserExec
    return;

  fClassMask       = vHeader->GetTriggerMask();
  fClassMaskNext50 = vHeader->GetTriggerMaskNext50();
  fRunNumber       = vEvent->GetRunNumber();

  fL0Inputs        = vHeader->GetL0TriggerInputs();
  fL1Inputs        = vHeader->GetL1TriggerInputs();
  fL2Inputs        = vHeader->GetL2TriggerInputs();

  fBCID            = vHeader->GetBunchCrossNumber();
  fOrbitID         = vHeader->GetOrbitNumber();
  fPeriod          = vHeader->GetPeriodNumber();
  fTimeStamp       = vHeader->GetTimeStamp();
}

void AliAnalysisTaskVdM::ADV0::FillInvalid() {
  fTime[0] = fTime[1] = -10240.0f;
  fBB[0] = fBG[0] = fBB[1] = fBG[1] = -1;
  for (Int_t bc=0; bc<21; ++bc)
    fPFBBA[bc] = fPFBBC[bc] = fPFBGA[bc] = fPFBGC[bc] = 0;
}

void AliAnalysisTaskVdM::ADV0::FillAD(const AliVEvent *vEvent, AliTriggerAnalysis &trigAna) {
  fDecisionOnline[0]  = trigAna.ADTrigger(vEvent, AliTriggerAnalysis::kCSide, kFALSE);
  fDecisionOnline[1]  = trigAna.ADTrigger(vEvent, AliTriggerAnalysis::kASide, kFALSE);

  fDecisionOffline[0] = trigAna.ADTrigger(vEvent, AliTriggerAnalysis::kCSide, kTRUE);
  fDecisionOffline[1] = trigAna.ADTrigger(vEvent, AliTriggerAnalysis::kASide, kTRUE);

  const AliVAD *esdAD = vEvent->GetADData();
  if (!esdAD) {
    FillInvalid();
    return;
  }
  fTime[0] = esdAD->GetADCTime();
  fTime[1] = esdAD->GetADATime();

  fBB[0] = fBB[1] = fBG[0] = fBG[1] = 0;
  for (Int_t ch=0; ch<4; ++ch) {
    fBB[0] += (esdAD->GetBBFlag(ch  ) && esdAD->GetBBFlag(ch+ 4));
    fBB[1] += (esdAD->GetBBFlag(ch+8) && esdAD->GetBBFlag(ch+12));
    fBG[0] += (esdAD->GetBGFlag(ch  ) && esdAD->GetBGFlag(ch+ 4));
    fBG[1] += (esdAD->GetBGFlag(ch+8) && esdAD->GetBGFlag(ch+12));
  }

  for (Int_t bc=0; bc<21; ++bc) {
    fPFBBA[bc] = fPFBBC[bc] = fPFBGA[bc] = fPFBGC[bc] = 0;
    for (Int_t ch=0; ch<4; ++ch) {
      fPFBBC[bc] += (esdAD->GetPFBBFlag(ch, bc) && esdAD->GetPFBBFlag(ch+4, bc));
      fPFBGC[bc] += (esdAD->GetPFBGFlag(ch, bc) && esdAD->GetPFBGFlag(ch+4, bc));

      fPFBBA[bc] += (esdAD->GetPFBBFlag(ch+8, bc) && esdAD->GetPFBBFlag(ch+12, bc));
      fPFBGA[bc] += (esdAD->GetPFBGFlag(ch+8, bc) && esdAD->GetPFBGFlag(ch+12, bc));
    }
  }
}
void AliAnalysisTaskVdM::ADV0::FillV0(const AliVEvent *vEvent, AliTriggerAnalysis &trigAna) {
  fDecisionOnline[0]  = trigAna.V0Trigger(vEvent, AliTriggerAnalysis::kCSide, kFALSE);
  fDecisionOnline[1]  = trigAna.V0Trigger(vEvent, AliTriggerAnalysis::kASide, kFALSE);

  fDecisionOffline[0] = trigAna.V0Trigger(vEvent, AliTriggerAnalysis::kCSide, kTRUE);
  fDecisionOffline[1] = trigAna.V0Trigger(vEvent, AliTriggerAnalysis::kASide, kTRUE);

  const AliVVZERO *esdV0 = vEvent->GetVZEROData();
  if (!esdV0) {
    FillInvalid();
    return;
  }

  fTime[0] = esdV0->GetV0CTime();
  fTime[1] = esdV0->GetV0ATime();

  fBB[0] = fBB[1] = fBG[0] = fBG[1] = 0;
  for (Int_t ch=0; ch<64; ++ch) {
    fBB[ch/32] += esdV0->GetBBFlag(ch);
    fBG[ch/32] += esdV0->GetBGFlag(ch);
  }

  for (Int_t bc=0; bc<21; ++bc) {
    fPFBBA[bc] = fPFBBC[bc] = fPFBGA[bc] = fPFBGC[bc] = 0;
    for (Int_t ch=0; ch<32; ++ch) {
      fPFBBC[bc] += esdV0->GetPFBBFlag(ch,    bc);
      fPFBGC[bc] += esdV0->GetPFBGFlag(ch,    bc);
      fPFBBA[bc] += esdV0->GetPFBBFlag(ch+32, bc);
      fPFBGA[bc] += esdV0->GetPFBGFlag(ch+32, bc);
    }
  }
}

AliAnalysisTaskVdM::AliAnalysisTaskVdM(const char *name)
  : AliAnalysisTaskSE(name)
  , fTreeBranchNames("")
  , fTriggerAnalysis()
  , fList(nullptr)
  , fTE(nullptr)
  , fIR1InteractionMap()
  , fIR2InteractionMap()
  , fVertexSPD()
  , fVertexTPC()
  , fVertexTracks()
  , fVertexTracksUnconstrained()
  , fTriggerIRs("AliTriggerIR", 3)
  , fFiredTriggerClasses()
  , fTreeData()
{
  for (Int_t i=0; i<kNHist;++i) {
    fHist[i] = nullptr;
  }
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskVdM::~AliAnalysisTaskVdM()
{
  const AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man && man->GetAnalysisType() == AliAnalysisManager::kProofAnalysis)
    return;

  fTriggerIRs.Delete();

  SafeDelete(fList);
  SafeDelete(fTE);
}

void AliAnalysisTaskVdM::SetBranches(TTree* t) {
  t->Branch("AliAnalysisTaskVdM::TreeData", &fTreeData);
  if (fTreeBranchNames.Contains("IRMap")) {
    t->Branch("IR1InteractionMap", &fIR1InteractionMap, 32000, 0);
    t->Branch("IR2InteractionMap", &fIR2InteractionMap, 32000, 0);
  }
  if (fTreeBranchNames.Contains("VertexSPD"))
    t->Branch("VertexSPD", &fVertexSPD, 32000, 0);
  if (fTreeBranchNames.Contains("VertexTPC"))
    t->Branch("VertexTPC", &fVertexTPC, 32000, 0);
  if (fTreeBranchNames.Contains("VertexTracks"))
    t->Branch("VertexTracks", &fVertexTracks, 32000, 0);
  if (fTreeBranchNames.Contains("VertexTracksUnconstrained"))
      t->Branch("VertexTracksUnconstrained", &fVertexTracksUnconstrained, 32000, 0);
  if (fTreeBranchNames.Contains("TriggerIR"))
    t->Branch("TriggerIRs", &fTriggerIRs, 32000, 0);
}

void AliAnalysisTaskVdM::UserCreateOutputObjects()
{
  fList = new TList;
  fList->SetOwner(kTRUE);
  fList->SetName(GetListName());
  fHist[kHistTrig] = new TH1D("HTrig", ";trigger class index", 102, -1.5, 100.5);
  fHist[kHistTrig]->SetStats(0);
  fList->Add(fHist[kHistTrig]);
  PostData(1, fList);

  TDirectory *owd = gDirectory;
  TFile *fSave = OpenFile(1);
  fTE = new TTree(GetTreeName(), "");
  SetBranches(fTE);
  PostData(2, fTE);
  owd->cd();
}

void AliAnalysisTaskVdM::NotifyRun()
{
  AliDebugF(5, "run %d", fCurrentRunNumber);
}

class TClonesArrayGuard {
public:
  TClonesArrayGuard(TClonesArray &a)
    : fA(a) {}
  ~TClonesArrayGuard() {
    fA.Delete();
  }
private:
  TClonesArrayGuard(const TClonesArrayGuard&);
  TClonesArrayGuard& operator=(const TClonesArrayGuard&);
  TClonesArray& fA;
} ;

void AliAnalysisTaskVdM::FillTriggerIR(const AliESDHeader* esdHeader) {
  if (!esdHeader)
    return;
  // store trigger IR for up to +-1 orbits around the event
  for (Int_t i=0,j=0,n=esdHeader->GetTriggerIREntries(); i<n; ++i) {
    const AliTriggerIR *ir = esdHeader->GetTriggerIR(i);
    if (!ir || TMath::Abs(Int_t(ir->GetOrbit()&0xFFFF) - Int_t(fTreeData.fEventInfo.fOrbitID)) > 1)
      continue;
    new(fTriggerIRs[j++]) AliTriggerIR(*ir);
  }
}
void AliAnalysisTaskVdM::UserExec(Option_t *)
{
  AliVEvent *event = InputEvent();
  if (!event) {
    AliFatal("nullptr == event");
    return;
  }

  AliVEvent* vEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!vEvent) {
    AliFatal("nullptr == vEvent");
    return;
  }
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(vEvent);
  const Bool_t isESD(esdEvent);

  if (!isESD)
    AliFatal("this task can run on ESD data only");

  // input handler
  const AliAnalysisManager* man(AliAnalysisManager::GetAnalysisManager());
  if (!man) {
    AliFatal("nullptr == man");
    return;
  }

  const AliVHeader *vHeader = vEvent->GetHeader();
  if (!vHeader) {
    AliFatal("nullptr == vHeader");
    return;
  }

  if (vHeader->GetEventType() != AliRawEventHeaderBase::kPhysicsEvent)
    return;

  fHist[kHistTrig]->Fill(-1); // # all analyzed events in underflow bin
  for (Int_t i=0; i<50; ++i) {
    const ULong64_t mask(1ULL<<i);
    if ((vHeader->GetTriggerMask() & mask) == mask)
      fHist[kHistTrig]->Fill(i);
    if ((vHeader->GetTriggerMaskNext50() & mask) == mask)
      fHist[kHistTrig]->Fill(50+i);
  }
  PostData(1, fList);

  fFiredTriggerClasses = vEvent->GetFiredTriggerClasses();

  fTreeData.fEventInfo.Fill(vEvent);
  fTreeData.fIsIncompleteDAQ = vEvent->IsIncompleteDAQ();
  fTreeData.fV0Info.FillV0(vEvent, fTriggerAnalysis);
  fTreeData.fADInfo.FillAD(vEvent, fTriggerAnalysis);

  fVertexSPD    = *esdEvent->GetPrimaryVertexSPD();
  fVertexTPC    = *esdEvent->GetPrimaryVertexTPC();
  fVertexTracks = *esdEvent->GetPrimaryVertexTracks();
  revertex(esdEvent, kFALSE);
  fVertexTracksUnconstrained = *esdEvent->GetPrimaryVertexTracks();

  TClonesArrayGuard guardTriggerIR(fTriggerIRs);
  FillTriggerIR(dynamic_cast<const AliESDHeader*>(vHeader));

  fIR1InteractionMap = vHeader->GetIRInt1InteractionMap();
  fIR2InteractionMap = vHeader->GetIRInt2InteractionMap();

  fTE->Fill();
  PostData(2, fTE);
}

void AliAnalysisTaskVdM::Terminate(Option_t*)
{
  fList  = dynamic_cast<TList*>(GetOutputData(1));
  if (!fList)
    Error("Terminate","fList is not available");

  fTE  = dynamic_cast<TTree*>(GetOutputData(2));
  if (!fTE)
    Error("Terminate","fTE is not available");
}
