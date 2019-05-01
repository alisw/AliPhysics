// -*- C++ indent-tabs-mode:nil; -*-
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
#include "AliStack.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliRawEventHeaderBase.h"
#include "AliESDFMD.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"

#include "AliAnalysisTaskDiffCrossSections.h"

ClassImp(AliAnalysisTaskDiffCrossSections);
ClassImp(AliAnalysisTaskDiffCrossSections::TreeData);
ClassImp(AliAnalysisTaskDiffCrossSections::MCInfo);
ClassImp(AliAnalysisTaskDiffCrossSections::PseudoTrack);
ClassImp(AliAnalysisTaskDiffCrossSections::PseudoTracks);

Int_t AliAnalysisTaskDiffCrossSections::PseudoTrack::Compare(const TObject *obj) const {
  if (this == obj)
    return 0;
  const PseudoTrack *p = (const PseudoTrack*)obj;
  if (fEta == p->Eta())
    return 0;
  return (*this < *p ? -1 : +1);
}
void AliAnalysisTaskDiffCrossSections::PseudoTrack::Print(Option_t* ) const {
  Printf("PT (eta=%5.2f phi=%5.2f charge=%8.2f charge2=%8.2f detFlags=0x%08X)", fEta, fPhi, fCharge, fCharge2, fDetFlags);
}
AliAnalysisTaskDiffCrossSections::PseudoTracks::PseudoTracks()
  : TObject()
  , fTracks("AliAnalysisTaskDiffCrossSections::PseudoTrack", 10000) {
  //
}
AliAnalysisTaskDiffCrossSections::PseudoTracks::~PseudoTracks() {
  //
}

void AliAnalysisTaskDiffCrossSections::PseudoTracks::AddTrack(const PseudoTrack& t) {
  new (fTracks[fTracks.GetEntriesFast()]) PseudoTrack(t);
}
void  AliAnalysisTaskDiffCrossSections::PseudoTracks::Clear(Option_t *opt) {
  fTracks.Clear(opt);
}
void  AliAnalysisTaskDiffCrossSections::PseudoTracks::Delete(Option_t *opt) {
  fTracks.Delete(opt);
}
Int_t AliAnalysisTaskDiffCrossSections::PseudoTracks::RemoveTracks(UInt_t mask) {
  Int_t nRemoved = 0;
  if (mask) {
    for (Int_t i=0, nt=fTracks.GetEntriesFast(); i<nt; ++i) {
      if (GetTrackAt(i).DetFlags() & mask) {
        fTracks.RemoveAt(i);
        ++nRemoved;
      }
    }
  }
  if (nRemoved) {
    // avoid calling TClonesArray::Compress() which is expensive
    fTracks.Sort();
  }
  return nRemoved;
}
const AliAnalysisTaskDiffCrossSections::PseudoTrack& AliAnalysisTaskDiffCrossSections::PseudoTracks::operator[](Int_t idx) const {
  const PseudoTrack* p = reinterpret_cast<const PseudoTrack* >(fTracks.UncheckedAt(idx));
  // if (!p)
  //   AliFatalF("idx=%d is outside of the allowed range [0,%d)", idx, fTracks.GetSize());
  return *p;
}
void AliAnalysisTaskDiffCrossSections::PseudoTracks::SortIfNeeded() {
  Float_t etaLast=0;
  Bool_t isSorted = kTRUE;
  for (Int_t i=0, n=fTracks.GetEntriesFast(); i<n && isSorted; ++i) {
    const PseudoTrack &t = GetTrackAt(i);
    isSorted = (i ? t.Eta() >= etaLast : kTRUE);
    etaLast = t.Eta();
  }
  if (!isSorted)
    fTracks.Sort();
}

TVector3 AliAnalysisTaskDiffCrossSections::GetADPseudoTrack(Int_t ch) {
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

void AliAnalysisTaskDiffCrossSections::PseudoTracks::FindAcceptance(UInt_t mask, Float_t &etaAccL, Float_t &etaAccR, const TVector3 &vertexPosition) const {
  // default is for V0:
  etaAccL = -3.7; // ADC: -7.15
  etaAccR = +5.1; // ADA: +6.40
  if ((mask & kAD) != 0) {
    const Double_t phiVertex = vertexPosition.Phi();
    if ((mask & kADC) != 0) {
      TVector3 vL = GetADPseudoTrack(-1);
      vL.RotateZ(-TMath::Pi()/2 + phiVertex);
      etaAccL = (vL - vertexPosition).Eta();
    }
    if ((mask & kADA) != 0) {
      TVector3 vR = GetADPseudoTrack(-9);
      vR.RotateZ(-TMath::Pi()/2 + phiVertex);
      etaAccR = (vR - vertexPosition).Eta();
    }
  }
}
void AliAnalysisTaskDiffCrossSections::PseudoTracks::FindAcceptance(UInt_t mask, Float_t &etaAccL, Float_t &etaAccR) const {
  etaAccL = ((mask & kADC) ? -7.15 : -3.7);
  etaAccR = ((mask & kADA) ? +6.40 : +5.1);
}

Bool_t SelectAll(const AliAnalysisTaskDiffCrossSections::PseudoTrack &t) {
  return kTRUE;
}

Int_t AliAnalysisTaskDiffCrossSections::PseudoTracks::ClassifyEventBits(Int_t &iEtaL, Int_t &iEtaR, Float_t &etaGap, Float_t &etaGapCenter,
                                                                        UInt_t mask, const TBits& bits) const {
  iEtaL = iEtaR = -1;
  etaGap = etaGapCenter = 0.0f;

  Int_t eventType = 0; // -1 -> 1L, +1 -> 1R, +2 -> 2A
  const Int_t nt = fTracks.GetEntriesFast();
  AliDebugF(5, "ClassifyEvent: nT=%d mask=0x%08X", nt, mask);

  if (nt == 0) // no pseudo-tracks
    return eventType;

  Float_t etaAccL, etaAccR;
  FindAcceptance(mask, etaAccL, etaAccR);
  // boundaries of acceptance:
  AliDebugF(5, "etaAccLR= %.3f %.3f", etaAccL, etaAccR);

  // fTracks.Sort(); already sorted

  // find etaL,etaR indices
  for (Int_t i=0; i<nt; ++i) {
    const PseudoTrack &t = GetTrackAt(i);
    if ((t.DetFlags() & mask) == 0 || !bits.TestBitNumber(i))
      continue;

    if (AliDebugLevelClass() >= 5)
      t.Print();

    iEtaL = (iEtaL == -1 ? i : iEtaL);
    iEtaR = i;
  }

  AliDebugF(5, "ClassifyEvent: iL,iR=%d,%d", iEtaL, iEtaR);
  if (iEtaL == -1) // no pseudo-tracks passing cuts
    return eventType;

  const PseudoTrack &tL = GetTrackAt(iEtaL);
  const PseudoTrack &tR = GetTrackAt(iEtaR);
  const Float_t etaR=tR.Eta(), etaL=tL.Eta();
  const Bool_t isOneTrack = (etaR - etaL < 0.5 &&
                             TMath::Abs(std::fmod(tR.Phi()-tL.Phi(), TMath::Pi())) < TMath::Pi()/2);
  if (isOneTrack) {
    const Float_t etaCenter = 0.5*(etaL + etaR);
    eventType = 2*(etaCenter > 0.0f) - 1;
    AliDebugF(5, "eventType=%2d center=%6.3f etaLR=%6.3f,%6.3f",
           eventType, etaCenter, etaL, etaR);
  } else {
    const Float_t dL = etaL - etaAccL;
    const Float_t dR = etaAccR - etaR;
    etaGap = etaGapCenter = 0.0;
    for (Int_t i=0; i<nt; ++i) {
      const PseudoTrack &t1 = GetTrackAt(i);
      if ((t1.DetFlags() & mask) == 0 || !bits.TestBitNumber(i))
        continue;
      const Float_t eta1 = t1.Eta();
      for (Int_t j=i+1; j<nt; ++j) {
        const PseudoTrack &t2 = GetTrackAt(j);
        if ((t2.DetFlags() & mask) != 0 && bits.TestBitNumber(j)) {
          const Float_t eta2 = t2.Eta();
          if (eta2 - eta1 > etaGap) {
            etaGap       =      eta2 - eta1;
            etaGapCenter = 0.5*(eta2 + eta1);
          }
          break;
        }
      }
    }
    if (etaGap > dL && etaGap > dR) {
      eventType =  2;
    } else if ((TMath::Abs(etaL) < 1.0f &&
                TMath::Abs(etaR) < 1.0f)) {
      eventType =  2;
    } else if (etaR < +1.0f) {
      eventType = -1;
    } else if (etaL > -1.0f) {
      eventType = +1;
    } else {
      eventType =  2;
    }
    AliDebugF(5, "eventType=%2d gap=%6.3f etaLR=%6.3f,%6.3f dLR=%6.3f,%6.3f",
           eventType, etaGap, etaL, etaR, dL, dR);
  }
  return eventType;
}

void AliAnalysisTaskDiffCrossSections::MCInfo::Fill(const AliMCEvent* mcEvent, TString mcType) {
  fEventType = fEventTypeGen = kInvalid;
  if (!mcEvent)
    AliFatal("NULL == mcEvent");

  const AliGenEventHeader* h(mcEvent->GenEventHeader());
  if (!h)
    AliFatal("NULL == h");

  h->Print();
  Int_t pdgDiff_p = -1;
  if (mcType.Contains("Pythia6")) {
    pdgDiff_p = 9902210;
    const AliGenPythiaEventHeader* ph = dynamic_cast<const AliGenPythiaEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case 92:
      fEventTypeGen = kSDR;
      break;
    case 93:
      fEventTypeGen = kSDL;
      break;
    case 94:
      fEventTypeGen = kDD;
      break;
    case 91:
      fEventTypeGen = kElastic;
      break;
    default:
      fEventTypeGen = kND;
   }
  }
  if (mcType.Contains("Pythia8")) {
    pdgDiff_p = 9902210;
    const AliGenPythiaEventHeader* ph = dynamic_cast<const AliGenPythiaEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case 103:
      fEventTypeGen = kSDR;
      break;
    case 104:
      fEventTypeGen = kSDL;
      break;
    case 105:
      fEventTypeGen = kDD;
      break;
    case 106:
      fEventTypeGen = kCD;
      break;
    case 102:
      fEventTypeGen = kElastic;
      break;
    default:
      fEventTypeGen = kND;
    }
  }
  if (mcType.Contains("PHOJET")) {
    const AliGenDPMjetEventHeader* ph = dynamic_cast<const AliGenDPMjetEventHeader* >(h);
    if (!ph)
      AliFatal("NULL == ph");
    switch (ph->ProcessType()) {
    case  5:
      fEventTypeGen = kSDR;
      break;
    case  6:
      fEventTypeGen = kSDL;
      break;
    case  7:
      fEventTypeGen = kDD;
      break;
    case  4:
      fEventTypeGen = kCD;
      break;
    case  2:
      fEventTypeGen = kElastic;
      break;
    default:
      fEventTypeGen = kND;
    }
  }

  fEventType = fEventTypeGen;

  AliStack* stack  = dynamic_cast<AliStack *>(const_cast<AliMCEvent*>(mcEvent)->Stack());
  if (!stack)
    AliFatal("NULL == stack");

  for (Int_t i=0; i<2; ++i)
    fDiffMass[i] = fDiffMassGen[i] = -1.0;

  // determine SD diffractive mass from MC generator
  TLorentzVector v;
  for (Int_t i=0, counter=0, n=stack->GetNprimary(); i<n; ++i) {
    const TParticle *p = const_cast<AliStack*>(stack)->Particle(i);
    if (!p)
      continue;
    if (pdgDiff_p > 0 && p->GetPdgCode() == pdgDiff_p) {
      p->Momentum(v);
      AliInfoF("found diff(%d): m=%f eta=%f", counter++, v.M(), v.Eta());
      fDiffMassGen[v.Eta()>0] = v.M(); // eta>0,<0 <-> SDL,SDR
    }
  }

  // determine SD mass without coorperation from MC generator
  for (Int_t i=0; i<2; ++i)
    fDiffMass[i] = fDiffMassGen[i];

  Int_t    side =  0;
  Double_t mass = -1;
  const Bool_t isSD = FindSingleDiffraction(stack, mcType, side, mass);
  AliInfoF("    isSD=%d side=%d mass=%6.1f type=%.0f mL=%6.1f mR=%6.1f",
           isSD, side, mass,
           fEventType, fDiffMass[0], fDiffMass[1]);
  if (isSD) {
    fEventType         = (side==1 ? kSDR : kSDL);
    fDiffMass[side==1] = mass;
    fDiffMass[side!=1] = -1.0;
    AliInfoF("*** isSD=%d side=%d mass=%6.1f type=%.0f mL=%6.1f mR=%6.1f",
             isSD, side, mass,
             fEventType, fDiffMass[0], fDiffMass[1]);
  }
}

TLorentzVector GetV(const TParticle* part) { // recomputes the energy of the 4-momentum
  TLorentzVector v;
  part->Momentum(v);
  v.SetXYZM(v.Px(), v.Py(), v.Pz(), part->GetMass());
  return v;
}

// based on code by M. Poghosyan
Bool_t AliAnalysisTaskDiffCrossSections::MCInfo::FindSingleDiffraction(AliStack *stack, TString mcType,
                                                                       Int_t &side, Double_t &mass) const
{
  side =  0;
  mass = -1;
  const Int_t pdgStable[20] = {
    22,   // Photon
    11,   // Electron
    12,   // Electron Neutrino
    13,   // Muon
    14,   // Muon Neutrino
    15,   // Tau
    16,   // Tau Neutrino
    211,  // Pion
    321,  // Kaon
    311,  // K0
    130,  // K0s
    310,  // K0l
    2212, // Proton
    2112, // Neutron
    3122, // Lambda_0
    3112, // Sigma Minus
    3222, // Sigma Plus
    3312, // Xsi Minus
    3322, // Xsi0
    3334  // Omega
  };
  const Int_t kNstable=sizeof(pdgStable)/sizeof(Int_t);

  Int_t  idx[2]  = {   -1,   -1 }; // index of particles with yMin, yMax
  Double_t y[2]  = {  0.0,  0.0 }; // yMin, yMax
  Double_t E_cms = 0.0;
  for (Int_t i0=(mcType.Contains("PHOJET") ? 2 : 0), i=i0, n=stack->GetNtrack(); i<n; i++) {
    const TParticle *part = stack->Particle(i);
    if (!part)
      continue;

    const Bool_t isCollProton = (i-i0 < 2 && part->GetPdgCode() == 2212);
    if (isCollProton)
      E_cms += GetV(part).E(); //part->Energy();

    if (part->GetStatusCode() != 1)
      continue;

    const Int_t pdg = TMath::Abs(part->GetPdgCode());
    Bool_t isStable = kFALSE;
    for (Int_t i1=0; i1<kNstable && !isStable; ++i1)
      isStable = (pdg == pdgStable[i1]);

    if (!isStable)
      continue;

    const Double_t yp = GetV(part).Rapidity(); //part->Y();
    if (idx[0] == -1) {
      y[0]   = yp;
      idx[0] = i;
    } else if (yp < y[0]) {
      y[0]   = yp;
      idx[0] = i;
    }
    if (idx[1] == -1) {
      y[1]   = yp;
      idx[1] = i;
    } else if (yp > y[1]) {
      y[1]   = yp;
      idx[1] = i;
    }
  } // next particle

  if (idx[0] < 0 && idx[1] < 0)
    return kFALSE;

  y[0] = TMath::Abs(y[0]);
  y[1] = TMath::Abs(y[1]);

  const TParticle *p[2] = {
    stack->Particle(idx[0]),
    stack->Particle(idx[1])
  };
  const Int_t pdg[2] = {
    p[0]->GetPdgCode(),
    p[1]->GetPdgCode()
  };

  Int_t i=-1;
  if (pdg[0]==2212 && pdg[1]==2212) {
    i = idx[y[1]>y[0]];
    if (y[1] == y[0])
      i = idx[gRandom->Uniform(0,1) > 0.5];
  } else if (pdg[0]==2212) {
    i = idx[0];
  } else if (pdg[1]==2212) {
    i = idx[1];
  }
  if (i < 0)
    return kFALSE;

  AliDebugF(5, "i=%4d idx= %4d,%4d y=%6.3f,%6.3f pdg=%4i,%4i", i, idx[0], idx[1], y[0], y[1], pdg[0], pdg[1]);

  const TParticle  *psel = stack->Particle(i);
  const TLorentzVector v = GetV(psel);
  const Double_t E = v.E(); //psel->Energy();
  const Double_t P = v.P(); //psel->P();

  const Double_t M2 = (E_cms-E-P)*(E_cms-E+P);
  mass = (M2 > 0.0 ? TMath::Sqrt(M2) : -2.0);
  if (i == idx[0])
    side = +1;

  if (i == idx[1])
    side = -1;

  return kTRUE;
}

void AliAnalysisTaskDiffCrossSections::EventInfo::Fill(const AliESDEvent* esdEvent, const char* inputFileName) {
  const AliESDHeader *esdHeader = esdEvent->GetHeader();
  if (NULL == esdHeader) // this case is already handled in UserExec
    return;

  fClassMask         = esdHeader->GetTriggerMask();
  fClassMaskNext50   = esdHeader->GetTriggerMaskNext50();

  fRunNumber         = esdEvent->GetRunNumber();
  fEventNumberInFile = esdEvent->GetEventNumberInFile();

  fL0Inputs          = esdHeader->GetL0TriggerInputs();
  fL1Inputs          = esdHeader->GetL1TriggerInputs();
  fL2Inputs          = esdHeader->GetL2TriggerInputs();

  fBCID              = esdHeader->GetBunchCrossNumber();
  fOrbitID           = esdHeader->GetOrbitNumber();
  fPeriod            = esdHeader->GetPeriodNumber();
  fTimeStamp         = esdHeader->GetTimeStamp();

  fInputFileName    = inputFileName;
}

void AliAnalysisTaskDiffCrossSections::VtxInfo::Fill(const AliESDVertex *vtx) {
  if (vtx) {
    fX      = vtx->GetX();
    fY      = vtx->GetY();
    fZ      = vtx->GetZ();
    fNcontr = vtx->GetNContributors();
  } else {
    fX      =  0.0f;
    fY      =  0.0f;
    fZ      =  0.0f;
    fNcontr = -4;
    AliErrorClass("NULL == vtx");
  }
}

void AliAnalysisTaskDiffCrossSections::ADV0::FillInvalid() {
  for (Int_t i=0; i<2; ++i) {
    fTime[i] = -10240.0f;
    fCharge[i] = fBBOnline[i] = fBBOffline[i] = fBGOnline[i] = fBGOffline[i] = -1;
  }
  for (Int_t ch=0; ch<64; ++ch) {
    fBBFlagsOnline[ch] = fBBFlagsOffline[ch] = fCharges[ch] = -1;
    fTimes[ch] = -10240.0f;
  }
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

  fBBOnline[0]  = fBBOnline[1]  = fBGOnline[0]  = fBGOnline[1]  = 0;
  fBBOffline[0] = fBBOffline[1] = fBGOffline[0] = fBGOffline[1] = 0;
  for (Int_t ch=0; ch<16; ++ch)
    fBBFlagsOnline[ch] = fBBFlagsOffline[ch] = fCharges[ch] = 0;

  for (Int_t ch=0; ch<4; ++ch) {
    fBBOnline[0]  += (esdAD->GetBBFlag(ch  )  && esdAD->GetBBFlag(ch+ 4));
    fBBOnline[1]  += (esdAD->GetBBFlag(ch+8)  && esdAD->GetBBFlag(ch+12));
    fBBOffline[0] += (esdAD->BBTriggerADC(ch) && esdAD->BBTriggerADC(ch+4));
    fBBOffline[1] += (esdAD->BBTriggerADA(ch) && esdAD->BBTriggerADA(ch+4));
    fBGOnline[0]  += (esdAD->GetBGFlag(ch  )  && esdAD->GetBGFlag(ch+ 4));
    fBGOnline[1]  += (esdAD->GetBGFlag(ch+8)  && esdAD->GetBGFlag(ch+12));
    fBGOffline[0] += (esdAD->BGTriggerADC(ch) && esdAD->BGTriggerADC(ch+4));
    fBGOffline[1] += (esdAD->BGTriggerADA(ch) && esdAD->BGTriggerADA(ch+4));
  }
  for (Int_t ch=0; ch<16; ++ch) {
      fBBFlagsOnline [ch] = esdAD->GetBBFlag(ch);
      fBBFlagsOffline[ch] = (ch < 8
                             ? esdAD->BBTriggerADC(ch)
                             : esdAD->BBTriggerADA(ch-8));
      fCharges[ch] = TMath::Nint(esdAD->GetAdc(ch));
      fTimes[ch]   = esdAD->GetTime(ch);
  }
  fCharge[0] = fCharge[1] = 0.0;
  for (Int_t ch=0; ch<16; ++ch)
    fCharge[ch/8] += esdAD->GetAdc(ch);
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

  fBBOnline[0]  = fBBOnline[1]  = fBGOnline[0]  = fBGOnline[1] = 0;
  fBBOffline[0] = fBBOffline[1] = fBGOffline[0] = fBGOffline[1] = 0;
  fCharge[0] = fCharge[1] = 0.0;
  for (Int_t ch=0; ch<16; ++ch)
    fBBFlagsOnline[ch] = fBBFlagsOffline[ch] = fCharges[ch] = 0;

  for (Int_t ch=0; ch<64; ++ch) {
    fBBOnline[ch/32]   += esdV0->GetBBFlag(ch);
    fBGOnline[ch/32]   += esdV0->GetBGFlag(ch);
    fBBOffline[ch/32]  += (ch < 32 ? esdV0->BBTriggerV0C(ch) : esdV0->BBTriggerV0A(ch-32));
    fBGOffline[ch/32]  += (ch < 32 ? esdV0->BGTriggerV0C(ch) : esdV0->BGTriggerV0A(ch-32));
    fCharge[ch/32]     += esdV0->GetAdc(ch);
    fBBFlagsOnline [ch] = esdV0->GetBBFlag(ch);
    fBBFlagsOffline[ch] = (ch < 32 ? esdV0->BBTriggerV0C(ch) : esdV0->BBTriggerV0A(ch-32));
    fCharges[ch]        = TMath::Nint(esdV0->GetAdc(ch));
    fTimes[ch]          = esdV0->GetTime(ch);
  }
}
void AliAnalysisTaskDiffCrossSections::FMDInfo::Fill(const AliESDEvent *esdEvent, Float_t fmdMultLowCut) {
  for (Int_t i=0; i<5; ++i)
    fMult[i] = 0.0f;

  const AliESDFMD* esdFMD = esdEvent->GetFMDData();
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
          if (fmdMult < fmdMultLowCut || fmdMult == AliESDFMD::kInvalidMult)
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
  , fFMDMultLowCut(0.3)
  , fTimeChargeCuts(16)
  , fTriggerAnalysis()
  , fAnalysisUtils()
  , fTE(NULL)
  , fFastOrMap()
  , fFiredChipMap()
  , fESDVZERO()
  , fESDAD()
  , fTreeData()
  , fMCInfo()
  , fIR1InteractionMap()
  , fIR2InteractionMap()
  , fMeanVtxPos(3)
  , fMeanVtxCov(3,3)
  , fMeanVtxU(3,3)
  , fEventType(0)
  , fIdxEtaL(-1)
  , fIdxEtaR(-1)
  , fEtaL(0)
  , fEtaR(0)
  , fEtaGap(0)
  , fEtaGapCenter(0)
{
  fTimeChargeCuts.SetOwner(kTRUE);
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

Bool_t AliAnalysisTaskDiffCrossSections::DoTimeChargeCut(Int_t ch, Float_t time, Float_t charge) const
{
  const TCutG *cut = dynamic_cast<const TCutG*>(fTimeChargeCuts.At(ch));
  // it no cut is set, require a time measurement (AliADReconstructor::kInvalidTime=-1024)
  return (cut ? cut->IsInside(time, charge) : (time > -1023.0f));
}

void AliAnalysisTaskDiffCrossSections::SetBranches(TTree* t) {
  t->Branch("TreeData", &fTreeData);
  if (fUseBranch.Contains("FastOrMap"))
    t->Branch("FastOrMap",    &fFastOrMap,    32000, 0);
  if (fUseBranch.Contains("FiredChipMap"))
    t->Branch("FiredChipMap", &fFiredChipMap, 32000, 0);
  if (fUseBranch.Contains("ESDVZERO"))
    t->Branch("ESDVZERO", &fESDVZERO, 32000, 0);
  if (fUseBranch.Contains("ESDAD"))
    t->Branch("ESDAD", &fESDAD, 32000, 0);
  if (fUseBranch.Contains("IRMap")) {
    t->Branch("IR1InteractionMap", &fIR1InteractionMap, 32000, 0);
    t->Branch("IR2InteractionMap", &fIR2InteractionMap, 32000, 0);
  }
  if (fIsMC)
    t->Branch("MCInfo", &fMCInfo);

  t->Branch("eventType",    &fEventType);
  t->Branch("idxEtaL",      &fIdxEtaL);
  t->Branch("idxEtaR",      &fIdxEtaR);
  t->Branch("etaL",         &fEtaL);
  t->Branch("etaR",         &fEtaR);
  t->Branch("etaGap",       &fEtaGap);
  t->Branch("etaGapCenter", &fEtaGapCenter);
}

void AliAnalysisTaskDiffCrossSections::UserCreateOutputObjects()
{
  TDirectory *owd = gDirectory;
  OpenFile(1);
  fTE = new TTree(GetTreeName(), GetTreeName());
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
TVector3 AliAnalysisTaskDiffCrossSections::GetV0PseudoTrack(Int_t ch) {
  // https://arxiv.org/pdf/1306.3130v2.pdf
  const Double_t r1[2][4] = {
    { 4.5, 7.3, 11.9, 19.5 }, // C-side ring 0-4
    { 4.3, 7.7, 13.9, 22.8 }  // A-side ring 0-4
  };
  const Double_t r2[2][4] = {
    { 7.1, 11.7, 19.3, 32.0 }, // C-side ring 0-4
    { 7.5, 13.7, 22.6, 41.2 }, // A-side ring 0-4
  };
  const Double_t z[2][4] = {
    { -86, -87, -88, -88 }, // C-side ring 0-4
    { 329, 329, 329, 329 }  // A-side ring 0-4
  };

  const Bool_t isForAcceptance = ch<0;
  ch = (ch<0 ? -ch : ch);
  const Int_t side   = ch/32;
  const Int_t ring   = (ch%32)/8;
  const Int_t sector = (ch%8);

  if (isForAcceptance)
    return TVector3(0, r1[1][ring], z[side][ring]);

  const Double_t r   = gRandom->Uniform(r1[side][ring], r2[side][ring]);
  const Double_t phi = gRandom->Uniform(sector*TMath::Pi()/4, (1+sector)*TMath::Pi()/4);
  return TVector3(r*TMath::Cos(phi),
                  r*TMath::Sin(phi),
                  z[side][ring]);
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
    TString tok = "";;
    Ssiz_t from = 0;
    Bool_t selected = kFALSE;
    while (fTriggerSelection.Tokenize(tok, from, "|") && !selected) {
      selected |= esdEvent->GetFiredTriggerClasses().Contains(tok);
    }

    AliInfo(Form("selected: %d %s", selected, esdEvent->GetFiredTriggerClasses().Data()));
    if (!selected)
      return;
  }

  fTreeData.fEventInfo.Fill(esdEvent, inputHandler->GetInputFileName());

  fTreeData.fPhysSelBits              = inputHandler->IsEventSelected();
  fTreeData.fIsIncompleteDAQ          = esdEvent->IsIncompleteDAQ();
  fTreeData.fIsSPDClusterVsTrackletBG = fAnalysisUtils.IsSPDClusterVsTrackletBG(esdEvent);
  fTreeData.fV0Info.FillV0(esdEvent, fTriggerAnalysis);
  fTreeData.fFMDInfo.Fill(esdEvent, fFMDMultLowCut);
  fTreeData.fADInfo.FillAD(esdEvent, fTriggerAnalysis);

  fFastOrMap    = mult->GetFastOrFiredChips();
  fFiredChipMap = mult->GetFiredChipMap();

  fIR1InteractionMap = esdHeader->GetIRInt1InteractionMap();
  fIR2InteractionMap = esdHeader->GetIRInt2InteractionMap();

  fTreeData.fEventInfo.fnTrklet         = mult->GetNumberOfTracklets();
  fTreeData.fEventInfo.fnSPDClusters[0] = mult->GetNumberOfITSClusters(0);
  fTreeData.fEventInfo.fnSPDClusters[1] = mult->GetNumberOfITSClusters(1);

  if (fIsMC)
    fMCInfo.Fill(MCEvent(), fMCType);

  const AliESDVertex *esdVertex = esdEvent->GetPrimaryVertexSPD();
  fTreeData.fVtxInfo.Fill(esdVertex);

  TVector3 vertexPosition(esdVertex->GetX(),
                          esdVertex->GetY(),
                          esdVertex->GetZ());
  if (esdVertex->GetNContributors() < 1) {
    vertexPosition = GetRandomVtxPosition();
    fTreeData.fVtxInfo.fX = vertexPosition.X();
    fTreeData.fVtxInfo.fY = vertexPosition.Y();
    fTreeData.fVtxInfo.fZ = vertexPosition.Z();
  }

  UInt_t mask = 0; // for selecting pseudo-tracks
  fTreeData.fPseudoTracks.Clear();

  std::vector<Double_t> eta;
  std::vector<Double_t> phi;

  // SPD
  if (fDetectorsUsed.Contains("SPD") && mult) {
    mask |= PseudoTracks::kSPD;
    for (Int_t i=0, n=mult->GetNumberOfTracklets(); i<n; ++i) {
      phi.push_back(mult->GetPhi(i));
      eta.push_back(mult->GetEta(i));
      fTreeData.fPseudoTracks.AddTrack(PseudoTrack(mult->GetEta(i),
						   mult->GetPhi(i),
						   1.0f, // charge=1 for SPD
						   0.0f,
						   PseudoTracks::kSPD |
						   PseudoTracks::kOffline));
    }
  }
  // VZERO
  const AliESDVZERO *esdVZERO = esdEvent->GetVZEROData();
  if (fDetectorsUsed.Contains("V0") && esdVZERO) {
    fESDVZERO = *esdVZERO;
    mask |= fDetectorsUsed.Contains("V0C")*PseudoTracks::kV0C;
    mask |= fDetectorsUsed.Contains("V0A")*PseudoTracks::kV0A;
    for (Int_t ch=0; ch<64; ++ch) {
      if (esdVZERO->GetBBFlag(ch)) {
	const TVector3 v = GetV0PseudoTrack(ch)-vertexPosition;
	phi.push_back(v.Phi());
	eta.push_back(v.Eta());
	const Bool_t isOffline = (ch < 32
				  ? esdVZERO->BBTriggerV0C(ch%32)
				  : esdVZERO->BBTriggerV0A(ch%32));
	fTreeData.fPseudoTracks.AddTrack(PseudoTrack(v.Eta(),
						     v.Phi(),
						     esdVZERO->GetAdc(ch),
						     0.0f,
						     (ch< 32)*PseudoTracks::kV0C |
						     (ch>=32)*PseudoTracks::kV0A |
						     PseudoTracks::kOnline | (isOffline*PseudoTracks::kOffline)));
      }
    }
  }
  // FMD
  const AliESDFMD* esdFMD = esdEvent->GetFMDData();
  if (fDetectorsUsed.Contains("FMD") && esdFMD) {
    mask |= PseudoTracks::kFMD1;
    mask |= PseudoTracks::kFMD2i;
    mask |= PseudoTracks::kFMD2o;
    mask |= PseudoTracks::kFMD3i;
    mask |= PseudoTracks::kFMD3o;
    for (UShort_t d=1; d<4; ++d) {
      const UShort_t nRng(d == 1 ? 1 : 2);
      for (UShort_t i=0; i<nRng; ++i) {
	const Char_t   r    = (i == 0 ? 'I' : 'O'); // ring
	const UShort_t nSec = (i == 0 ?  20 :  40); // sector
	const UShort_t nStr = (i == 0 ? 512 : 256); // strip
	for (UShort_t s=0; s<nSec; ++s) {
	  for (UShort_t t=0; t<nStr; ++t) {
	    const Float_t fmdMult = esdFMD->Multiplicity(d, r, s, t);
	    if (fmdMult < fFMDMultLowCut || fmdMult == AliESDFMD::kInvalidMult)
	      continue;
	    phi.push_back(esdFMD->Phi(d, r, s, t)*TMath::Pi()/180.0);
	    eta.push_back(esdFMD->Eta(d, r, s, t));
	    fTreeData.fPseudoTracks.AddTrack(PseudoTrack(esdFMD->Eta(d, r, s, t),
							 esdFMD->Phi(d, r, s, t)*TMath::Pi()/180.0,
							 fmdMult,
							 Float_t( UInt_t(t) + (UInt_t(s) << 9)), // encode s,t
							 (d==1)       *PseudoTracks::kFMD1  |
							 (d==2)*(i==0)*PseudoTracks::kFMD2i |
							 (d==2)*(i==1)*PseudoTracks::kFMD2o |
							 (d==3)*(i==0)*PseudoTracks::kFMD3i |
							 (d==3)*(i==1)*PseudoTracks::kFMD3o |
							 PseudoTracks::kOffline));
	  }
	}
      }
    }
  }
  // AD
  const AliESDAD *esdAD = esdEvent->GetADData();
  if (fDetectorsUsed.Contains("AD") && esdAD) {
    fESDAD = *esdAD;
    mask |= fDetectorsUsed.Contains("ADC")*PseudoTracks::kADC;
    mask |= fDetectorsUsed.Contains("ADA")*PseudoTracks::kADA;
    for (Int_t side=0; side<2; ++side) {
      if (side==0 && !fDetectorsUsed.Contains("ADC"))
	continue;
      if (side==1 && !fDetectorsUsed.Contains("ADA"))
	continue;
      for (Int_t module=0; module<4; ++module) {
	const Int_t ch = 8*side + module;
        const Bool_t flag1 = esdAD->GetBBFlag(ch)   || DoTimeChargeCut(ch,   esdAD->GetTime(ch),   esdAD->GetAdc(ch));
        const Bool_t flag2 = esdAD->GetBBFlag(ch+4) || DoTimeChargeCut(ch+4, esdAD->GetTime(ch+4), esdAD->GetAdc(ch+4));
	if (flag1 && flag2) {
	  const TVector3 v = GetADPseudoTrack(ch)-vertexPosition;
	  phi.push_back(v.Phi());
	  eta.push_back(v.Eta());
	  const Bool_t isOffline = (side == 0
				    ? (esdAD->BBTriggerADC(module) && esdAD->BBTriggerADC(module+4))
				    : (esdAD->BBTriggerADA(module) && esdAD->BBTriggerADA(module+4)));
          const Bool_t isFlagNotBB = (!esdAD->GetBBFlag(ch) || !esdAD->GetBBFlag(ch+4));
	  fTreeData.fPseudoTracks.AddTrack(PseudoTrack(v.Eta(),
						       v.Phi(),
						       esdAD->GetAdc(ch),
						       esdAD->GetAdc(ch+4),
						       (ch< 8)*PseudoTracks::kADC |
						       (ch>=8)*PseudoTracks::kADA |
						       (PseudoTracks::kOnline)                |
                                                       (PseudoTracks::kOffline   * isOffline) |
                                                       (PseudoTracks::kFlagNotBB * isFlagNotBB)));
	}
      }
    }
  }

  fEventType = fTreeData.fPseudoTracks.ClassifyEvent(fIdxEtaL, fIdxEtaR, fEtaGap, fEtaGapCenter, mask, SelectAll);
  fEtaL = (fIdxEtaL >= 0 ? fTreeData.fPseudoTracks[fIdxEtaL].Eta() : -999.0f);
  fEtaR = (fIdxEtaR >= 0 ? fTreeData.fPseudoTracks[fIdxEtaR].Eta() : -999.0f);

  fTE->Fill();
  PostData(1, fTE);
}

void AliAnalysisTaskDiffCrossSections::Terminate(Option_t*)
{
  fTE  = dynamic_cast<TTree*>(GetOutputData(1));
  if (NULL == fTE)
    Error("Terminate","fTE is not available");
}
