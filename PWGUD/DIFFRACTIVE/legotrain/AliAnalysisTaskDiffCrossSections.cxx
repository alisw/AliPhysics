#include <memory>
#include <cmath>
#include <vector>

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TDecompChol.h>
#include <TRandom.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliRawEventHeaderBase.h"
#include "AliESDVZERO.h"
#include "AliESDFMD.h"
#include "AliESDAD.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include "AliAnalysisTaskDiffCrossSections.h"

ClassImp(AliAnalysisTaskDiffCrossSections);
ClassImp(AliAnalysisTaskDiffCrossSections::TreeData);
ClassImp(AliAnalysisTaskDiffCrossSections::MCInfo);

void AliAnalysisTaskDiffCrossSections::MCInfo::Fill(const AliMCEvent* mcEvent, TString mcType) {
  fEventType = kInvalid;
  if (!mcEvent)
    AliFatal("NULL == mcEvent");

  const AliGenEventHeader* h(mcEvent->GenEventHeader());
  if (!h)
    AliFatal("NULL == h");

  Int_t pdgDiff_p = -1;
  if (mcType.Contains("Pythia")) {
    pdgDiff_p = 9902210;
    const AliGenPythiaEventHeader* ph = dynamic_cast<const AliGenPythiaEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case 92:
      fEventType = kSDR;
      break;
    case 93:
      fEventType = kSDL;
      break;
    case 94:
      fEventType = kDD;
      break;
    case 91:
      fEventType = kElastic;
      break;
    default:
      fEventType = kND;
   }
  }
  if (mcType.Contains("Pythia8")) {
    pdgDiff_p = 9902210;
    const AliGenPythiaEventHeader* ph = dynamic_cast<const AliGenPythiaEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case 103:
      fEventType = kSDR;
      break;
    case 104:
      fEventType = kSDL;
      break;
    case 105:
      fEventType = kDD;
      break;
    case 106:
      fEventType = kCD;
      break;
    case 102:
      fEventType = kElastic;
      break;
    default:
      fEventType = kND;
    }
  }
  if (mcType.Contains("PHOJET")) {
    const AliGenDPMjetEventHeader* ph = dynamic_cast<const AliGenDPMjetEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case  5:
      fEventType = kSDR;
      break;
    case  6:
      fEventType = kSDL;
      break;
    case  7:
      fEventType = kDD;
      break;
    case  4:
      fEventType = kCD;
      break;
    case  2:
      fEventType = kElastic;
      break;
    default:
      fEventType = kND;
    }
  }

  const AliStack* stack  = dynamic_cast<const AliStack *>(const_cast<AliMCEvent*>(mcEvent)->Stack());
  if (!stack)
    AliFatal("NULL == stack");
  const Int_t npart = stack->GetNprimary();

  TLorentzVector v;
  fDiffMass[0] = fDiffMass[1] = -1.0f;
  for (Int_t i=0, counter=0, n=stack->GetNprimary(); i<n; ++i) {
    const TParticle *p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p)
      continue;
    if (pdgDiff_p > 0 && p->GetPdgCode() == pdgDiff_p) {
      p->Momentum(v);
      AliInfoF("found diff(%d): m=%f eta=%f", counter++, v.M(), v.Eta());
      fDiffMass[v.Eta()>0] = v.M(); // eta>0,<0 -> SDL,SDR to be verified
    }
  }
}

void AliAnalysisTaskDiffCrossSections::EventInfo::Fill(const AliESDEvent* esdEvent) {
  const AliESDHeader *esdHeader = esdEvent->GetHeader();
  if (NULL == esdHeader) // this is already dealt with in UserExec
    return;

  fClassMask       = esdHeader->GetTriggerMask();
  fClassMaskNext50 = esdHeader->GetTriggerMaskNext50();

  fRunNumber       = esdEvent->GetRunNumber();

  fL0Inputs        = esdHeader->GetL0TriggerInputs();
  fL1Inputs        = esdHeader->GetL1TriggerInputs();
  fL2Inputs        = esdHeader->GetL2TriggerInputs();

  fBCID            = esdHeader->GetBunchCrossNumber();
  fOrbitID         = esdHeader->GetOrbitNumber();
  fPeriod          = esdHeader->GetPeriodNumber();
  fTimeStamp       = esdHeader->GetTimeStamp();
}

void AliAnalysisTaskDiffCrossSections::VtxInfo::Fill(const AliESDVertex *vtx) {
  if (vtx) {
    fZ      = vtx->GetZ();
    fNcontr = vtx->GetNContributors();
    AliInfoClass("NULL != vtx");
  } else {
    fZ      =  0.0f;
    fNcontr = -4;
    AliErrorClass("NULL == vtx");
  }
}

void AliAnalysisTaskDiffCrossSections::ADV0::FillInvalid() {
  fTime[0] = fTime[1] = -10240.0f;
  fCharge[0] = fCharge[1] = fBB[0] = fBG[0] = fBB[1] = fBG[1] = -1.0f;
}

void AliAnalysisTaskDiffCrossSections::ADV0::FillAD(const AliESDEvent *esdEvent, AliTriggerAnalysis &trigAna) {
  fDecisionOnline[0]  = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kCSide, kFALSE);
  fDecisionOnline[1]  = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kASide, kFALSE);

  fDecisionOffline[0] = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kCSide, kTRUE);
  fDecisionOffline[1] = trigAna.ADTrigger(esdEvent, AliTriggerAnalysis::kASide, kTRUE);

  const AliESDAD *esdAD = esdEvent->GetADData();
  if (NULL == esdAD) {
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
  fCharge[0] = fCharge[1] = 0.0;
  for (Int_t ch=0; ch<16; ++ch)
    fCharge[ch/8] += esdAD->GetMultiplicity(ch);
}
void AliAnalysisTaskDiffCrossSections::ADV0::FillV0(const AliESDEvent *esdEvent, AliTriggerAnalysis &trigAna) {
  fDecisionOnline[0]  = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kCSide, kFALSE);
  fDecisionOnline[1]  = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kASide, kFALSE);

  fDecisionOffline[0] = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kCSide, kTRUE);
  fDecisionOffline[1] = trigAna.V0Trigger(esdEvent, AliTriggerAnalysis::kASide, kTRUE);

  const AliESDVZERO *esdV0 = esdEvent->GetVZEROData();
  if (NULL == esdV0) {
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
  fCharge[0] = fCharge[1] = 0.0;
  for (Int_t ch=0; ch<64; ++ch)
    fCharge[ch/32] += esdV0->GetMultiplicity(ch);
}
void AliAnalysisTaskDiffCrossSections::FMDInfo::Fill(const AliESDEvent *esdEvent) {
  for (Int_t i=0; i<5; ++i)
    fMult[i] = 0.0f;

  const AliESDFMD* esdFMD = esdEvent->GetFMDData();
  const Float_t fFMDLowCut = 0.2;
  const Int_t idx[4][2] = {
    {0, 0},
    {0, 0},
    {1, 2},
    {3, 4}
  };
  for (UShort_t d=1; d<4; ++d) {
    const UShort_t nRng(d == 1 ? 1 : 2);
    for (UShort_t i=0; i<nRng; ++i) {
      const Char_t   r    = (i == 0 ? 'I' : 'O'); // ring
      const UShort_t nSec = (i == 0 ?  20 :  40); // sector
      const UShort_t nStr = (i == 0 ? 512 : 256); // strip
      for (UShort_t s=0; s<nSec; ++s) {
	for (UShort_t t=0; t<nStr; ++t) {
	  const Float_t fmdMult = esdFMD->Multiplicity(d, r, s, t);
	  if (fmdMult < fFMDLowCut || fmdMult == AliESDFMD::kInvalidMult)
	    continue;
	  fMult[idx[d][i]] += 1;
	}
      }
    }
  }
}


AliAnalysisTaskDiffCrossSections::AliAnalysisTaskDiffCrossSections(const char *name)
  : AliAnalysisTaskSE(name)
  , fIsMC(kFALSE)
  , fMCType("")
  , fTriggerSelection("")
  , fDetectorsUsed("ADA ADC V0 FMD SPD")
  , fUseBranch("")
  , fTriggerAnalysis()
  , fAnalysisUtils()
  , fTE(NULL)
  , fFastOrMap()
  , fFiredChipMap()
  , fTreeData()
  , fMCInfo()
  , fMeanVtxPos(3)
  , fMeanVtxCov(3,3)
  , fMeanVtxU(3,3)
{
  fTriggerAnalysis.SetAnalyzeMC(fIsMC);

  DefineOutput(1, TTree::Class());
}

AliAnalysisTaskDiffCrossSections::~AliAnalysisTaskDiffCrossSections()
{
  const AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (NULL != man && man->GetAnalysisType() == AliAnalysisManager::kProofAnalysis)
    return;

  if (NULL != fTE)
    delete fTE;
  fTE = NULL;
}

void AliAnalysisTaskDiffCrossSections::SetBranches(TTree* t) {
  t->Branch("AliAnalysisTaskDG::TreeData", &fTreeData);
  if (fUseBranch.Contains("FastOrMap"))
    t->Branch("FastOrMap",    &fFastOrMap,    32000, 0);
  if (fUseBranch.Contains("FiredChipMap"))
    t->Branch("FiredChipMap", &fFiredChipMap, 32000, 0);
  if (fIsMC)
    t->Branch("AliAnalysisTaskDG::MCInfo", &fMCInfo);

  t->Branch("eventType",    &fEventType);
  t->Branch("etaL",         &fEtaL);
  t->Branch("etaR",         &fEtaR);
  t->Branch("etaGap",       &fEtaGap);
  t->Branch("etaGapCenter", &fEtaGapCenter);
}

void AliAnalysisTaskDiffCrossSections::UserCreateOutputObjects()
{
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTE = new TTree;
  fTE->SetName(GetTreeName());
  SetBranches(fTE);
  PostData(1, fTE);
  owd->cd();
}

void AliAnalysisTaskDiffCrossSections::NotifyRun() {
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    const Int_t runs[] = { // taken from PWGPP/AliTaskCDBconnect.cxx
      139674, // < 2010
      170718, // < 2011
      194479, // < 2012
      199999, // < 2013
      208504, // < 2014
      247170, // < 2015
      999999  // < 2016
    };
    Int_t year=0;
    for (Int_t i=0, n=sizeof(runs)/sizeof(Int_t); i<n; ++i) {
      if (fCurrentRunNumber < runs[i]) {
	year = 2010+i;
	break;
      }
    }
    man->SetDefaultStorage(Form("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/%d/OCDB/", year));
  }
  man->SetRun(fCurrentRunNumber);

  AliCDBEntry *entry = man->Get("GRP/Calib/MeanVertex");
  AliESDVertex *vtx = (AliESDVertex*)entry->GetObject();

  Double_t pos[3] = { 0,0,0 };
  vtx->GetXYZ(pos);
  fMeanVtxPos.SetElements(pos);

  Double_t cov[6] = { 0,0,0,0,0,0 }; // fCovXX,fCovXY,fCovYY,fCovXZ,fCovYZ,fCovZZ;
  vtx->GetCovMatrix(cov);
  fMeanVtxCov(0,0) = cov[0];
  fMeanVtxCov(0,1) = fMeanVtxCov(1,0) = cov[1];
  fMeanVtxCov(1,1) = cov[2];
  fMeanVtxCov(0,2) = fMeanVtxCov(2,0) = cov[3];
  fMeanVtxCov(1,2) = fMeanVtxCov(2,1) = cov[4];
  fMeanVtxCov(2,2) = cov[5];

  TDecompChol l(fMeanVtxCov);
  const Bool_t success = l.Decompose();
  if (!success)
    AliFatal("TDecompChol::Decompose() failed");
  fMeanVtxU = l.GetU();
  fMeanVtxU.Transpose(fMeanVtxU);
}

TVector3 AliAnalysisTaskDiffCrossSections::GetRandomVtxPosition() const {
  // see http://pdg.lbl.gov/2013/reviews/rpp2013-rev-monte-carlo-techniques.pdf
  const Double_t rnd[3] = {
    gRandom->Gaus(),
    gRandom->Gaus(),
    gRandom->Gaus()
  };
  return TVector3((fMeanVtxPos + fMeanVtxU*TVectorD(3, rnd)).GetMatrixArray());
}
TVector3 GetV0PseudoTrack(Int_t ch) {
  const Double_t r1[2][4] = {
    { 4.6, 7.1, 11.7, 19.5 }, // C-side ring 0-4
    { 4.3, 7.7, 13.9, 22.8 }  // A-side ring 0-4
  };
  const Double_t r2[2][4] = {
    { 7.1, 11.6, 19.1, 32.0 }, // C-side ring 0-4
    { 7.5, 13.7, 22.6, 41.2 }, // A-side ring 0-4
  };
  const Double_t z[2] = {
    -87.15, // C_side
    329.00  // A-side
  };

  const Bool_t isForAcceptance = ch<0;
  ch = (ch<0 ? -ch : ch);
  const Int_t side   = ch/32;
  const Int_t ring   = (ch%32)/8;
  const Int_t sector = (ch%8);

  if (isForAcceptance)
    return TVector3(0, r1[1][ring], z[side]);

  const Double_t r   = gRandom->Uniform(r1[side][ring], r2[side][ring]);
  const Double_t phi = gRandom->Uniform(sector*TMath::Pi()/4, (1+sector)*TMath::Pi()/4);
  return TVector3(r*TMath::Cos(phi),
		  r*TMath::Sin(phi),
		  z[side]);
}

TVector3 GetADPseudoTrack(Int_t ch) {
  const Double_t z[2] = {
    -1954.4, // C-side
    +1696.7  // A-side
  };
  const Double_t innerRadius[2] = {
    3.6, // C-side
    6.2  // A-side
  };
  const Int_t signX[4] = { +1,-1,-1,+1 };
  const Int_t signY[4] = { +1,+1,-1,-1 };

  const Bool_t isForAcceptance = ch<0;
  ch = (ch<0 ? -ch : ch);
  const Int_t side   = ch/8;
  const Int_t module = (ch%4);

  if (isForAcceptance)
    return TVector3(0, innerRadius[side], z[side]);

  TVector2 v(0,0);
  const Int_t maxIter=100*1000;
  Int_t counter=0;
  for (; v.Mod() < innerRadius[side] && counter<maxIter; ++counter) {
    v.Set(gRandom->Uniform(0.54, 18.1+0.54),
	  gRandom->Uniform(0.0,  21.4));
  }
  // if (counter == maxIter)
  //   AliFatalGeneral("AD", "counter == maxIter");
  return TVector3(signX[module]*v.X(),
		  signY[module]*v.Y(),
		  z[side]);
}

void AliAnalysisTaskDiffCrossSections::UserExec(Option_t *)
{
  AliVEvent *event = InputEvent();
  if (NULL == event) {
    AliFatal("NULL == event");
    return;
  }

  // ESD Event
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (NULL == esdEvent) {
    AliFatal("NULL == esdEvent");
    return;
  }

  const AliAnalysisManager* man(AliAnalysisManager::GetAnalysisManager());
  if (NULL == man) {
    AliFatal("NULL == man");
    return;
  }

  AliESDInputHandler* inputHandler(dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler()));
  if (NULL == inputHandler) {
    AliFatal("NULL == inputHandler");
    return;
  }

  const AliESDHeader *esdHeader = esdEvent->GetHeader();
  if (NULL == esdHeader) {
    AliFatal("NULL == esdHeader");
    return;
  }

  if (kFALSE == fIsMC && esdHeader->GetEventType() != AliRawEventHeaderBase::kPhysicsEvent)
    return;

  const AliMultiplicity *mult = esdEvent->GetMultiplicity();
  if (NULL == mult) {
    AliFatal("NULL == mult");
    return;
  }

  if (!fIsMC) {
    std::unique_ptr<const TObjArray> split(fTriggerSelection.Tokenize("|"));
    Bool_t selected = kFALSE;
    for (Int_t i=0, n=split->GetEntries(); i<n && !selected; ++i)
      selected = esdEvent->GetFiredTriggerClasses().Contains(split->At(i)->GetName());

    AliInfo(Form("selected: %d %s", selected, esdEvent->GetFiredTriggerClasses().Data()));
    if (!selected)
      return;
  }

  fTreeData.fEventInfo.Fill(esdEvent);

  fTreeData.fPhysSelBits              = inputHandler->IsEventSelected();
  fTreeData.fIsIncompleteDAQ          = esdEvent->IsIncompleteDAQ();
  fTreeData.fIsSPDClusterVsTrackletBG = fAnalysisUtils.IsSPDClusterVsTrackletBG(esdEvent);
  fTreeData.fV0Info.FillV0(esdEvent, fTriggerAnalysis);
  fTreeData.fFMDInfo.Fill(esdEvent);
  fTreeData.fADInfo.FillAD(esdEvent, fTriggerAnalysis);

  fTreeData.fVtxInfo.Fill(esdEvent->GetPrimaryVertexSPD());

  fFastOrMap    = mult->GetFastOrFiredChips();
  fFiredChipMap = mult->GetFiredChipMap();

  fTreeData.fEventInfo.fnTrklet         = mult->GetNumberOfTracklets();
  fTreeData.fEventInfo.fnSPDClusters[0] = mult->GetNumberOfITSClusters(0);
  fTreeData.fEventInfo.fnSPDClusters[1] = mult->GetNumberOfITSClusters(1);

  if (fIsMC)
    fMCInfo.Fill(MCEvent(), fMCType);

  const AliESDVertex *esdVertex = esdEvent->GetPrimaryVertex();
  TVector3 vertexPosition(esdVertex->GetX(),
			  esdVertex->GetY(),
			  esdVertex->GetZ());
  if (esdVertex->GetNContributors() < 1)
    vertexPosition = GetRandomVtxPosition();

  std::vector<Double_t> eta;
  std::vector<Double_t> phi;

  // SPD
  if (fDetectorsUsed.Contains("SPD") && mult) {
    for (Int_t i=0, n=mult->GetNumberOfTracklets(); i<n; ++i) {
      phi.push_back(mult->GetPhi(i));
      eta.push_back(mult->GetEta(i));
    }
  }
  // VZERO
  const AliESDVZERO *esdVZERO = esdEvent->GetVZEROData();
  if (fDetectorsUsed.Contains("V0") && esdVZERO) {
    for (Int_t ch=0; ch<64; ++ch) {
      if (esdVZERO->GetBBFlag(ch)) {
	const TVector3 v = GetV0PseudoTrack(ch)-vertexPosition;
	phi.push_back(v.Phi());
	eta.push_back(v.Eta());
      }
    }
  }
  // FMD
  const AliESDFMD* esdFMD = esdEvent->GetFMDData();
  if (fDetectorsUsed.Contains("FMD") && esdFMD) {
    const Float_t fFMDLowCut = 0.2;
    for (UShort_t d=1; d<4; ++d) {
      const UShort_t nRng(d == 1 ? 1 : 2);
      for (UShort_t i=0; i<nRng; ++i) {
	const Char_t   r    = (i == 0 ? 'I' : 'O'); // ring
	const UShort_t nSec = (i == 0 ?  20 :  40); // sector
	const UShort_t nStr = (i == 0 ? 512 : 256); // strip
	for (UShort_t s=0; s<nSec; ++s) {
	  for (UShort_t t=0; t<nStr; ++t) {
	    const Float_t fmdMult = esdFMD->Multiplicity(d, r, s, t);
	    if (fmdMult < fFMDLowCut || fmdMult == AliESDFMD::kInvalidMult)
	      continue;
	    phi.push_back(esdFMD->Phi(d, r, s, t)*TMath::Pi()/180.0);
	    eta.push_back(esdFMD->Eta(d, r, s, t));
	  }
	}
      }
    }
  }
  // AD
  const AliESDAD *esdAD = esdEvent->GetADData();
  if (fDetectorsUsed.Contains("AD") && esdAD) {
    for (Int_t side=0; side<2; ++side) {
      if (side==0 && !fDetectorsUsed.Contains("ADC"))
	continue;
      if (side==1 && !fDetectorsUsed.Contains("ADA"))
	continue;
      for (Int_t module=0; module<4; ++module) {
	const Int_t ch = 8*side + module;
	if (esdAD->GetBBFlag(ch) &&
	    esdAD->GetBBFlag(ch+4)) {
	  const TVector3 v = GetADPseudoTrack(ch)-vertexPosition;
	  phi.push_back(v.Phi());
	  eta.push_back(v.Eta());
	}
      }
    }
  }

  // default is for V0:
  Double_t etaAccL = -3.7;
  Double_t etaAccR = +5.1;
  if (fDetectorsUsed.Contains("AD")) {
    const Double_t phiVertex = vertexPosition.Phi();
    if (fDetectorsUsed.Contains("ADC")) {
      TVector3 vL = GetADPseudoTrack(-1);
      vL.RotateZ(-TMath::Pi()/2 + phiVertex);
      etaAccL = (vL - vertexPosition).Eta();
    }
    if (fDetectorsUsed.Contains("ADA")) {
      TVector3 vR = GetADPseudoTrack(-9);
      vR.RotateZ(-TMath::Pi()/2 + phiVertex);
      etaAccR = (vR - vertexPosition).Eta();
    }
  }

  fEtaL = fEtaR = fEtaGap = fEtaGapCenter = 0.0;
  fEventType = 0;
  if (!eta.empty()) {
    std::sort(eta.begin(), eta.end());
    std::sort(phi.begin(), phi.end());
    fEtaL = eta.front();
    fEtaR = eta.back();

    const Double_t phiMin = phi.front();
    const Double_t phiMax = phi.back();
    // the following is OK for SPD,FMD,V0 but presumably needs to be generalized for AD
    const Bool_t isOneTrack = (fEtaR - fEtaL < 0.5 &&
			       TMath::Abs(std::fmod(phiMax-phiMin, TMath::Pi())) < TMath::Pi()/2);
    if (isOneTrack) {
      const Double_t etaCenter = 0.5*(fEtaL + fEtaR);
      fEventType = 2*(etaCenter > 0.0) - 1;
    } else {
      const Double_t dL = fEtaL - etaAccL;
      const Double_t dR = etaAccR - fEtaR;
      fEtaGap = 0.0;
      for (std::vector<double>::const_iterator i=eta.begin(), j=i+1, jend=eta.end(); j!=jend; ++i, ++j) {
	if (*j - *i > fEtaGap) {
	  fEtaGap       =      *j - *i;
	  fEtaGapCenter = 0.5*(*j + *i);
	}
      }
      if (fEtaGap > dL && fEtaGap > dR) {
	fEventType = 2;
      } else if ((TMath::Abs(fEtaL) < 1 &&
		  TMath::Abs(fEtaR) < 1)) {
	fEventType = 2;
      } else if (fEtaR < +1) {
	fEventType = -1;
      } else if (fEtaL > -1) {
	fEventType = +1;
      } else {
	fEventType = 2;
      }
    }
  }

  fTE->Fill();
  PostData(1, fTE);
}

void AliAnalysisTaskDiffCrossSections::Terminate(Option_t*)
{
  fTE  = dynamic_cast<TTree*>(GetOutputData(1));
  if (NULL == fTE)
    Error("Terminate","fTE is not available");
}
