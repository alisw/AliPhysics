#include <memory>
#include <array>
#include <functional>
#include <sstream>

#include <TRandom.h>
#include <TTree.h>
#include <TList.h>
#include <TH1.h>
#include <TH3.h>
#include <TFile.h>
#include <TString.h>
#include <TLorentzVector.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliTriggerIR.h"
#include "AliESDVertex.h"
#include "AliTriggerIR.h"

#include "AliAnalysisTaskDG.h"
#include "AliRawEventHeaderBase.h"
#include "AliVVZERO.h"
#include "AliVAD.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliITSsegmentationSPD.h"

ClassImp(AliAnalysisTaskDG);
ClassImp(AliAnalysisTaskDG::TreeData);
ClassImp(AliAnalysisTaskDG::TrackData);
ClassImp(AliAnalysisTaskDG::SPD_0STG);

void AliAnalysisTaskDG::EventInfo::Fill(const AliVEvent* vEvent) {
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

void AliAnalysisTaskDG::ADV0::FillInvalid() {
  fTime[0] = fTime[1] = -10240.0f;
  fBB[0] = fBG[0] = fBB[1] = fBG[1] = -1;
  for (Int_t bc=0; bc<21; ++bc)
    fPFBBA[bc] = fPFBBC[bc] = fPFBGA[bc] = fPFBGC[bc] = 0;
}

void AliAnalysisTaskDG::ADV0::FillAD(const AliVEvent *vEvent, AliTriggerAnalysis &trigAna) {
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
void AliAnalysisTaskDG::ADV0::FillV0(const AliVEvent *vEvent, AliTriggerAnalysis &trigAna) {
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

void AliAnalysisTaskDG::FMD::Fill(const AliVEvent *vEvent, AliTriggerAnalysis &trigAna) {
  fA = trigAna.FMDTrigger(vEvent, AliTriggerAnalysis::kASide);
  fC = trigAna.FMDTrigger(vEvent, AliTriggerAnalysis::kCSide);
}

template<typename A>
std::string array2string(const A& a, std::string fill="") {
  std::ostringstream oss;
  std::for_each(a.begin(), a.end(), [&](typename A::value_type e){ oss << e << fill; });
  return oss.str();
}

const TBits& AliAnalysisTaskDG::SPD_0STG::Fill(const TBits &bits) {
  std::array<Int_t, 20> l0;
  std::array<Int_t, 40> l1;

  l0.fill(0);
  l1.fill(0);

  for (Int_t i=0; i<400; ++i)
    l0[i/20] += bits.TestBitNumber(i);
  for (Int_t i=400; i<1200; ++i)
    l1[(i-400)/20] += bits.TestBitNumber(i);

  AliDebugF(5, "l0: %s", array2string(l0, " ").c_str());
  AliDebugF(5, "l1: %s", array2string(l1).c_str());

  fNPseudoTracklets = 0;
  std::array<Bool_t, 20> tr;
  for (Int_t i=0; i<20; ++i) {
    tr[i] = l0[i] & (l1[2*i] | l1[(2*i+1)%40] | l1[(2*i+2)%40] | l1[(2*i+39)%40]);
    fNPseudoTracklets += tr[i];
  }
  AliDebugF(5, "fNPseudoTracklets = %d (%s)", fNPseudoTracklets, array2string(tr).c_str());

  std::array<Bool_t, 10> match;
  for (Int_t j=0; j<10; ++j) {
    match[j] = kFALSE;
    for (Int_t i=0; i<20 && !match[j]; ++i) {
      match[j] = tr[i] & tr[(i+1+j)%20];
    }
  }
  fMaxDeltaPhi = -1;
  fMinDeltaPhi = -1;
  for (Int_t i=0; i<10; ++i) {
    if (match[i]) {
      fMinDeltaPhi = (fMinDeltaPhi <= 0 ? i+1 : fMinDeltaPhi);
      fMaxDeltaPhi = i+1;
    }
  }
  AliDebugF(5, "fNPseudoTracklets = %d (%s)", fNPseudoTracklets, array2string(tr).c_str());

  return bits;
}

void AliAnalysisTaskDG::FindChipKeys(AliESDtrack *tr, Short_t chipKeys[2], Int_t status[2]) {
  chipKeys[0] = chipKeys[1] = -1;
  status[0]   = status[1]   = -1;
  if (!tr)
    return;

  Int_t   idet=0;
  Float_t xloc=0, zloc=0;
  const AliITSsegmentationSPD seg;
  for (Int_t layer=0; layer<2; ++layer) {
    chipKeys[layer] = -1;
    status[layer]   = -1;
    const Int_t module = tr->GetITSModuleIndex(layer);
    if (module < 0)
      continue;

    tr->GetITSModuleIndexInfo(layer, idet, status[layer], xloc, zloc);
    if (status[layer] == 0)
      continue;

    Int_t off = seg.GetChipFromLocal(xloc, zloc);
    if (off < 0)
      continue;

    off = (layer==0 ? 4-off : off);
    AliDebugClassF(5, "layer=%d module=%10d idet=%3d xloc=%5.1f zloc=%5.1f off=%d status=%d chipKey=%4d",
		   layer, module, idet, xloc, zloc, off, status[layer], 5*idet + off);
    chipKeys[layer] = 5*idet + off;
  }
}

Int_t GetTrackSign(const AliVTrack *tr) {
  const AliAODTrack *aodTrack = dynamic_cast<const AliAODTrack*>(tr);
  return (aodTrack
	  ? aodTrack->Charge()
	  : tr->GetSign());
}
void AliAnalysisTaskDG::TrackData::Fill(AliVTrack *tr, AliPIDResponse *pidResponse=nullptr) {
  if (!tr || !pidResponse) {
    AliErrorF("tr=%p pidResponse=%p", tr, pidResponse);
    return;
  }
  fSign = GetTrackSign(tr);
  fPx   = tr->Px();
  fPy   = tr->Py();
  fPz   = tr->Pz();
  fITSsignal = tr->GetITSsignal();
  fTPCsignal = tr->GetTPCsignal();
  fTOFsignal = tr->GetTOFsignal();

  fPIDStatus[0] = pidResponse->CheckPIDStatus(AliPIDResponse::kITS, tr);
  fPIDStatus[1] = pidResponse->CheckPIDStatus(AliPIDResponse::kTPC, tr);
  fPIDStatus[2] = pidResponse->CheckPIDStatus(AliPIDResponse::kTOF, tr);

  for (Int_t i=0; i<AliPID::kSPECIES; ++i) {
    const AliPID::EParticleType particleType(static_cast<const AliPID::EParticleType>(i));
    fNumSigmaITS[i] = pidResponse->NumberOfSigmasITS(tr, particleType);
    fNumSigmaTPC[i] = pidResponse->NumberOfSigmasTPC(tr, particleType);
    fNumSigmaTOF[i] = pidResponse->NumberOfSigmasTOF(tr, particleType);
    AliDebugF(5, "%d %f %f %f", i, fNumSigmaITS[i], fNumSigmaTPC[i], fNumSigmaTOF[i]);
  }

  AliAnalysisTaskDG::FindChipKeys(dynamic_cast<AliESDtrack*>(tr), fChipKey, fStatus);

  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(tr);
  if (aodTrack)
    fFilterMap = aodTrack->GetFilterMap();
}

void AliAnalysisTaskDG::SetMaxTracksSave(Int_t m) {
  fTrackData.Expand(m);
  fMaxTracksSave = m;
}

AliAnalysisTaskDG::AliAnalysisTaskDG(const char *name)
  : AliAnalysisTaskSE(name)
  , fIsMC(kFALSE)
  , fTreeBranchNames("")
  , fTrackCutType("TPCOnly")
  , fTrackFilterMask(0)
  , fTriggerSelection("")
  , fTriggerSelectionSPD("")
  , fMaxTracksSave(4)
  , fTriggerAnalysis()
  , fAnalysisUtils()
  , fList(nullptr)
  , fTE(nullptr)
  , fIR1InteractionMap()
  , fIR2InteractionMap()
  , fFastOrMap()
  , fFiredChipMap()
  , fFiredTriggerClasses()
  , fVertexSPD()
  , fVertexTPC()
  , fVertexTracks()
  , fTOFHeader()
  , fTriggerIRs("AliTriggerIR", 3)
  , fSPD_0STG_Online()
  , fSPD_0STG_Offline()
  , fTrackData("AliAnalysisTaskDG::TrackData", fMaxTracksSave)
  , fMCTracks("TLorentzVector", 2)
  , fTrackCuts(nullptr)
{
  for (Int_t i=0; i<kNHist;++i) {
    fHist[i] = nullptr;
  }
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskDG::~AliAnalysisTaskDG()
{
  const AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man && man->GetAnalysisType() == AliAnalysisManager::kProofAnalysis)
    return;

  fTriggerIRs.Delete();
  fTrackData.Delete();
  fMCTracks.Delete();

  SafeDelete(fList);
  SafeDelete(fTE);
  SafeDelete(fTrackCuts);
}

void AliAnalysisTaskDG::SetBranches(TTree* t, Bool_t isAOD) {
  t->Branch("AliAnalysisTaskDG::TreeData",            &fTreeData);

  if (fTreeBranchNames.Contains("IRMap")) {
    t->Branch("IR1InteractionMap", &fIR1InteractionMap, 32000, 0);
    t->Branch("IR2InteractionMap", &fIR2InteractionMap, 32000, 0);
  }
  if (fTreeBranchNames.Contains("SPDMaps")) {
    t->Branch("FastOrMap",    &fFastOrMap,    32000, 0);
    t->Branch("FiredChipMap", &fFiredChipMap, 32000, 0);
  }
  if (fTreeBranchNames.Contains("0STG")) {
    t->Branch("0STG_Online",  &fSPD_0STG_Online,  32000, 0);
    t->Branch("0STG_Offline", &fSPD_0STG_Offline, 32000, 0);
  }
  if (fTreeBranchNames.Contains("VertexSPD")) {
    if (isAOD)
      t->Branch("VertexSPD", &fVertexSPD.second, 32000, 0);
    else
      t->Branch("VertexSPD", &fVertexSPD.first,  32000, 0);
  }
  if (fTreeBranchNames.Contains("VertexTPC")) {
    if (isAOD)
      t->Branch("VertexTPC", &fVertexTPC.second, 32000, 0);
    else
      t->Branch("VertexTPC", &fVertexTPC.first,  32000, 0);
  }
  if (fTreeBranchNames.Contains("VertexTracks")) {
    if (isAOD)
      t->Branch("VertexTracks", &fVertexTracks.second, 32000, 0);
    else
      t->Branch("VertexTracks", &fVertexTracks.first,  32000, 0);
  }
  if (fTreeBranchNames.Contains("TOFHeader")) {
    t->Branch("TOFHeader", &fTOFHeader, 32000, 0);
  }
  if (fTreeBranchNames.Contains("TriggerIR")) {
    t->Branch("TriggerIRs", &fTriggerIRs, 32000, 0);
  }

  if (fTreeBranchNames.Contains("Tracks")) {
    t->Branch("AliAnalysisTaskDG::TrackData", &fTrackData);
  }

  if (fIsMC) {
    t->Branch("TLorentzVector", &fMCTracks, 32000, 0);
  }

}

void AliAnalysisTaskDG::UserCreateOutputObjects()
{
  if (fTrackCutType == "ITSPureSA") {
    fTrackCuts  = AliESDtrackCuts::GetStandardITSPureSATrackCuts2010(kTRUE, kFALSE);
    fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kBoth);
  } else if (fTrackCutType == "TPCOnly") {
    fTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  } else if (fTrackCutType == "ITSTPC2011") {
    fTrackCuts  = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1);
  } else if (fTrackCutType == "ITSTPC2011_SPDboth") {
    fTrackCuts  = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1);
    fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kBoth);
  }

  if (fTrackCutType.BeginsWith("AODFilterMask")) {
    std::unique_ptr<const TObjArray> split(fTrackCutType.Tokenize("_"));
    if (split->GetEntries() != 2)
      AliFatalF("malformed track cut specification: '%s'", fTrackCutType.Data());
    fTrackFilterMask = atoi(split->At(1)->GetName());
  }
  if (!fTrackCuts && !fTrackFilterMask)
    AliFatal("no track cuts specified");


  fList = new TList;
  fList->SetOwner(kTRUE);
  fList->SetName(GetListName());
  fHist[kHistTrig] = new TH1D("HTrig", ";trigger class index", 102, -1.5, 100.5);
  fHist[kHistTrig]->SetStats(0);
  fList->Add(fHist[kHistTrig]);

  fHist[kHistSPDFiredTrk] = new TH3D("HSPDFiredTrk", fTriggerSelectionSPD+";chip key;BCmod4;mult",
				     1200, -0.5, 1199.5, 4, -0.5, 3.5, 10, -0.5, 9.5);
  fHist[kHistSPDFiredTrk]->SetStats(0);
  fList->Add(fHist[kHistSPDFiredTrk]);

  fHist[kHistSPDFOTrk] = new TH3D("HSPDFOTrk", fTriggerSelectionSPD+";chip key;BCmod4;mult",
				  1200, -0.5, 1199.5, 4, -0.5, 3.5, 10, -0.5, 9.5);
  fHist[kHistSPDFOTrk]->SetStats(0);
  fList->Add(fHist[kHistSPDFOTrk]);

  fHist[kHistSPDFiredTrkVsMult] = new TH3D("HSPDFiredTrkVsMult", fTriggerSelectionSPD+";chip key;BCmod4;log_{10}(number of tracklets)",
					   1200, -0.5, 1199.5, 4, -0.5, 3.5, 25, 0.0, 5.0);
  fHist[kHistSPDFiredTrkVsMult]->SetStats(0);
  fList->Add(fHist[kHistSPDFiredTrkVsMult]);

  fHist[kHistSPDFOTrkVsMult] = new TH3D("HSPDFOTrkVsMult", fTriggerSelectionSPD+";chip key;BCmod4;log_{10}(number of tracklets)",
					1200, -0.5, 1199.5, 4, -0.5, 3.5, 25, 0.0, 5.0);
  fHist[kHistSPDFOTrkVsMult]->SetStats(0);
  fList->Add(fHist[kHistSPDFOTrkVsMult]);

  fHist[kHistSPDFiredVsMult] = new TH3D("HSPDFiredVsMult", fTriggerSelectionSPD+";chip key;BCmod4;log_{10}(number of tracklets)",
					1200, -0.5, 1199.5, 4, -0.5, 3.5, 25, 0.0, 5.0);
  fHist[kHistSPDFiredVsMult]->SetStats(0);
  fList->Add(fHist[kHistSPDFiredVsMult]);

  fHist[kHistSPDFOVsMult] = new TH3D("HSPDFOVsMult", fTriggerSelectionSPD+";chip key;BCmod4;log_{10}(number of tracklets)",
				     1200, -0.5, 1199.5, 4, -0.5, 3.5, 25, 0.0, 5.0);
  fHist[kHistSPDFOVsMult]->SetStats(0);
  fList->Add(fHist[kHistSPDFOVsMult]);

  PostData(1, fList);

  TDirectory *owd = gDirectory;
  TFile *fSave = OpenFile(1);
  fTE = new TTree(GetTreeName(), "");
  SetBranches(fTE, fTrackFilterMask != 0);
  PostData(2, fTE);
  owd->cd();
}

void AliAnalysisTaskDG::NotifyRun()
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

void AliAnalysisTaskDG::FillTH3(Int_t idx, Double_t x, Double_t y, Double_t z, Double_t w) {
  if (idx < 0 || idx >= kNHist)
    AliFatalF("idx=%d", idx);
  TH3 *h = dynamic_cast<TH3*>(fHist[idx]);
  if (!h)
    AliFatal("h==nullptr");
  h->Fill(x, y, z, w);
}
void AliAnalysisTaskDG::FillSPDFOEffiencyHistograms(const AliESDEvent *esdEvent)
{
  if (!esdEvent)
    return;

  const AliESDHeader *esdHeader = esdEvent->GetHeader();
  const AliMultiplicity *mult = esdEvent->GetMultiplicity();

  Bool_t selectedForSPD = (fTriggerSelectionSPD == "");
  if (!selectedForSPD) { // trigger selection for SPD efficiency studies
    TString tcName = "";
    Ssiz_t  from   = 0;
    while (fTriggerSelectionSPD.Tokenize(tcName, from, "|") && !selectedForSPD)
      selectedForSPD = esdEvent->GetFiredTriggerClasses().Contains(tcName);
  }
  AliInfoF("selectedForSPD = %d", selectedForSPD);
  if (selectedForSPD) { // PF protection
    const AliVAD    *esdAD = esdEvent->GetADData();
    const AliVVZERO *esdV0 = esdEvent->GetVZEROData();
    Int_t nBB=0;
    for (Int_t bc=3; bc<=17 && !nBB; ++bc) {
      if (bc == 10)
	continue;
      for (Int_t ch=0; ch<4; ++ch) {
	nBB += (esdAD->GetPFBBFlag(ch,   bc) && esdAD->GetPFBBFlag(ch+ 4, bc));
	nBB += (esdAD->GetPFBBFlag(ch+8, bc) && esdAD->GetPFBBFlag(ch+12, bc));
      }
      for (Int_t ch=0; ch<64; ++ch)
	nBB += esdV0->GetPFBBFlag(ch, bc);
    }
    if (!nBB) {
      Int_t matched[1200] = { 0 };
      for (Int_t l=0; l<1200; ++l)
	matched[l] = 0;
      std::unique_ptr<AliESDtrackCuts> tc(AliESDtrackCuts::GetStandardITSPureSATrackCuts2010(kTRUE, kFALSE));
      std::unique_ptr<const TObjArray> oa(tc->GetAcceptedTracks(esdEvent));
      for (Int_t i=0, n=oa->GetEntries(); i<n; ++i) {
	AliESDtrack *tr = dynamic_cast<AliESDtrack*>(oa->At(i));
	Short_t chipKeys[2] = { -1, -1};
	Int_t   status[2]   = { -1, -1};
	AliAnalysisTaskDG::FindChipKeys(tr, chipKeys, status);
	for (Int_t layer=0; layer<2; ++layer) {
	  if (chipKeys[layer] >= 0 && chipKeys[layer]<1200 && status[layer] == 1)
	    matched[chipKeys[layer]] += 1;
	}
      }
      const Int_t    bcMod4         = (esdHeader->GetBunchCrossNumber() % 4);
      const Double_t log10Tracklets = (mult->GetNumberOfTracklets() > 0
				       ? TMath::Log10(mult->GetNumberOfTracklets())
				       : -1.0);
      for (Int_t chipKey=0; chipKey<1200; ++chipKey) {
	if (mult->TestFiredChipMap(chipKey)) {
	  FillTH3(kHistSPDFiredTrk,       chipKey, bcMod4, matched[chipKey]);
	  FillTH3(kHistSPDFiredTrkVsMult, chipKey, bcMod4, log10Tracklets, (matched[chipKey]>0));
	  FillTH3(kHistSPDFiredVsMult,    chipKey, bcMod4, log10Tracklets);
	}
	if (mult->TestFastOrFiredChips(chipKey)) {
	  FillTH3(kHistSPDFOTrk,       chipKey, bcMod4, matched[chipKey]);
	  FillTH3(kHistSPDFOTrkVsMult, chipKey, bcMod4, log10Tracklets, (matched[chipKey]>0));
	  FillTH3(kHistSPDFOVsMult,    chipKey, bcMod4, log10Tracklets);
	}
      }
    }
  }
}

void AliAnalysisTaskDG::FillTriggerIR(const AliESDHeader* esdHeader) {
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
void AliAnalysisTaskDG::UserExec(Option_t *)
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
  AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*>(vEvent);
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(vEvent);
  const Bool_t isESD(esdEvent);

  // input handler
  const AliAnalysisManager* man(AliAnalysisManager::GetAnalysisManager());
  if (!man) {
    AliFatal("nullptr == man");
    return;
  }

  AliESDInputHandler* inputHandlerESD(dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler()));
  AliAODInputHandler* inputHandlerAOD(dynamic_cast<AliAODInputHandler*>(man->GetInputEventHandler()));
  if (!inputHandlerESD && !inputHandlerAOD) {
    AliFatal("nullptr == iH");
    return;
  }

  AliPIDResponse* pidResponse = (inputHandlerESD
				 ? inputHandlerESD->GetPIDResponse()
				 : inputHandlerAOD->GetPIDResponse());
  if (!pidResponse) {
    AliFatal("nullptr == pidResponse");
    return;
  }

  const AliVHeader *vHeader = vEvent->GetHeader();
  if (!vHeader) {
    AliFatal("nullptr == vHeader");
    return;
  }

  if (kFALSE == fIsMC && vHeader->GetEventType() != AliRawEventHeaderBase::kPhysicsEvent)
    return;

  const AliVMultiplicity *mult = vEvent->GetMultiplicity();
  if (!mult) {
    AliFatal("nullptr == mult");
    return;
  }

  fHist[kHistTrig]->Fill(-1); // # analyzed events in underflow bin
  for (Int_t i=0; i<50; ++i) {
    const ULong64_t mask(1ULL<<i);
    if ((vHeader->GetTriggerMask() & mask))
      fHist[kHistTrig]->Fill(i);
    if ((vHeader->GetTriggerMaskNext50() & mask))
      fHist[kHistTrig]->Fill(50+i);
  }

  FillSPDFOEffiencyHistograms(dynamic_cast<AliESDEvent*>(vEvent));

  PostData(1, fList);

  Bool_t cutNotV0           = kFALSE;
  Bool_t useOnly2Trk        = kFALSE;
  Bool_t requireAtLeast2Trk = kFALSE;

  Bool_t selected = (fTriggerSelection == "");

  if (!selected) {
    // fTriggerSelection can be "CLASS1|CLASS2&NotV0|CLASS3&Only2Trk|CLASS4&AtLeast2Trk"
    Int_t sumCutNotV0(0);
    Int_t sumUseOnly2Trk(0);
    Int_t sumRequireAtLeast2Trk(0);

    Int_t   counter     = 0;
    TString tcName      = "";
    Ssiz_t  from_tcName = 0;
    while (fTriggerSelection.Tokenize(tcName, from_tcName, "|")) {
      TString tok      = "";
      Ssiz_t  from_tok = 0;
      if (!tcName.Tokenize(tok, from_tok, "&")) continue;
      if (!vEvent->GetFiredTriggerClasses().Contains(tok)) continue;
      if ( tcName.Tokenize(tok, from_tok, "&")) {
	sumCutNotV0           += (tok == "NotV0");
	sumUseOnly2Trk        += (tok == "Only2Trk");
	sumRequireAtLeast2Trk += (tok == "AtLeast2Trk");
	++counter;
      }
    }

    selected           = (counter != 0);
    cutNotV0           = (counter == sumCutNotV0);
    useOnly2Trk        = (counter == sumUseOnly2Trk);
    requireAtLeast2Trk = (counter == sumRequireAtLeast2Trk);
  }

  AliDebugF(5, "selected: %d (%d,%d,%d) %s ", selected,
	    cutNotV0, useOnly2Trk, requireAtLeast2Trk,
	    vEvent->GetFiredTriggerClasses().Data());
  if (!selected)
    return;

  fTreeData.fEventInfo.Fill(vEvent);

  fTreeData.fIsIncompleteDAQ          = vEvent->IsIncompleteDAQ();
  fTreeData.fIsSPDClusterVsTrackletBG = fAnalysisUtils.IsSPDClusterVsTrackletBG(vEvent);
  fTreeData.fIskMB                    = (inputHandlerESD
					 ? (inputHandlerESD->IsEventSelected() & AliVEvent::kMB)
					 : kFALSE);
  fTreeData.fV0Info.FillV0(vEvent, fTriggerAnalysis);
  fTreeData.fADInfo.FillAD(vEvent, fTriggerAnalysis);

  if (cutNotV0 &&
      (fTreeData.fV0Info.fDecisionOnline[0] != 0 ||
       fTreeData.fV0Info.fDecisionOnline[1] != 0))
    return;

  if (isESD) {
    fVertexSPD.first    = *esdEvent->GetPrimaryVertexSPD();
    fVertexTPC.first    = *esdEvent->GetPrimaryVertexTPC();
    fVertexTracks.first = *esdEvent->GetPrimaryVertexTracks();
  } else {
    fVertexSPD.second    = *aodEvent->GetPrimaryVertexSPD();
    fVertexTPC.second    = *aodEvent->GetPrimaryVertexTPC();
    fVertexTracks.second = *aodEvent->GetPrimaryVertex();
  }
  fTOFHeader    = *(vEvent->GetTOFHeader());

  TClonesArrayGuard guardTriggerIR(fTriggerIRs);
  FillTriggerIR(dynamic_cast<const AliESDHeader*>(vHeader));

  fFastOrMap    = fSPD_0STG_Online.Fill(mult->GetFastOrFiredChips());
  fFiredChipMap = fSPD_0STG_Offline.Fill(mult->GetFiredChipMap());

  fIR1InteractionMap = vHeader->GetIRInt1InteractionMap();
  fIR2InteractionMap = vHeader->GetIRInt2InteractionMap();

  for (Int_t i=0; i<4; ++i)
    fTreeData.fEventInfo.fnTrklet[i] = 0;
  for (Int_t i=0, n=mult->GetNumberOfTracklets(); i<n; ++i) {
    const Double_t eta = -TMath::Log(TMath::Tan(0.5*mult->GetTheta(i)));
    fTreeData.fEventInfo.fnTrklet[0] += 1;           // all tracklets
    fTreeData.fEventInfo.fnTrklet[1] += (eta <  -0.9);
    fTreeData.fEventInfo.fnTrklet[2] += (eta >= -0.9 && eta <= +0.9);
    fTreeData.fEventInfo.fnTrklet[3] += (eta >  +0.9);
  }
  fTreeData.fEventInfo.fnFO[1] = fFastOrMap.CountBits(400);
  fTreeData.fEventInfo.fnFO[0] = fFastOrMap.CountBits() - fTreeData.fEventInfo.fnFO[1];
  std::unique_ptr<TObjArray> oa;
  if (esdEvent) { // ESD
    oa = std::unique_ptr<TObjArray>(fTrackCuts->GetAcceptedTracks(esdEvent));
  } else { // AOD
    oa = std::unique_ptr<TObjArray>(new TObjArray);
    const AliAODEvent *aodEvent = dynamic_cast<const AliAODEvent*>(vEvent);
    for (Int_t i=0, n=aodEvent->GetNumberOfTracks(); i<n; ++i) {
      AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
      if (aodTrack &&
	  (aodTrack->GetFilterMap() & fTrackFilterMask) != 0                     &&
	  aodTrack->HasPointOnITSLayer(0) && aodTrack->HasPointOnITSLayer(1)     &&
	  (aodTrack->GetStatus() & AliVTrack::kITSrefit) == AliVTrack::kITSrefit &&
	  (aodTrack->GetStatus() & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit)
	oa->Add(aodTrack);
    }
  }
  fTreeData.fEventInfo.fnTrk   = oa->GetEntries();
  fTreeData.fEventInfo.fCharge = 0;

  if (useOnly2Trk && fTreeData.fEventInfo.fnTrk != 2)
    return;

  if (requireAtLeast2Trk && fTreeData.fEventInfo.fnTrk < 2)
    return;

  for (Int_t i=0, n=oa->GetEntries(); i<n; ++i)
    fTreeData.fEventInfo.fCharge += (GetTrackSign(dynamic_cast<AliVTrack*>(oa->At(i))) > 0
				     ? +1
				     : -1);

  TClonesArrayGuard guardTrackData(fTrackData);
  if (oa->GetEntries() <= fMaxTracksSave)  {
    for (Int_t i=0, n=TMath::Min(oa->GetEntries(), fMaxTracksSave); i<n; ++i)
      new(fTrackData[i]) TrackData(dynamic_cast<AliVTrack*>(oa->At(i)), pidResponse);
  }

  TClonesArrayGuard guardMCTracks(fMCTracks);
  if (fIsMC) {
    AliMCEvent *mcEvent = MCEvent();
    if (!mcEvent)
      AliFatal("nullptr ==mcEvent");

    Int_t counter = 0;
    for(Int_t i=0, n=mcEvent->GetNumberOfTracks(); i<n && counter<2; ++i) {
      AliMCParticle *p = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(i));
      if (!p) continue;
      new(fMCTracks[counter]) TLorentzVector;
      TLorentzVector *v = dynamic_cast<TLorentzVector*>(fMCTracks.At(counter));
      p->Particle()->Momentum(*v);
      ++counter;
    }
  } // fIsMC

  fTE->Fill();
  PostData(2, fTE);
}

void AliAnalysisTaskDG::Terminate(Option_t*)
{
  fList  = dynamic_cast<TList*>(GetOutputData(1));
  if (!fList)
    Error("Terminate","fList is not available");

  fTE  = dynamic_cast<TTree*>(GetOutputData(2));
  if (!fTE)
    Error("Terminate","fTE is not available");

}
