//
// Implementation file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#include "Riostream.h"

#include "TH1.h"
#include "TTree.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliMCEvent.h"

#include "AliRsnAnalysisPhi900GeV.h"


AliRsnAnalysisPhi900GeV::AliRsnAnalysisPhi900GeV(const char *name) :
  AliAnalysisTaskSE(name),
  fUseMC(kFALSE),
  fPDG(0),
  fIM(0.0),
  fPt(0.0),
  fY(0.0),
  fEta(0.0),
  fDCAr(1E6),
  fDCAz(1E6),
  fChi2(1E6),
  fNTPC(0),
  fMinTPC(-1E6),
  fMaxTPC( 1E6),
  fOutList(0x0),
  fHEvents(0x0),
  fTOFESD(kFALSE),
  fTOFSigma(210.0),
  fTOFmaker(0x0),
  fTOFSettings(AliRsnTOFT0maker::kNone)
{
//
// Constructor
//

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TList::Class());
}


AliRsnAnalysisPhi900GeV::AliRsnAnalysisPhi900GeV(const AliRsnAnalysisPhi900GeV& copy) :
  AliAnalysisTaskSE(copy),
  fUseMC(copy.fUseMC),
  fPDG(0),
  fIM(0.0),
  fPt(0.0),
  fY(0.0),
  fEta(0.0),
  fDCAr(copy.fDCAr),
  fDCAz(copy.fDCAz),
  fChi2(copy.fChi2),
  fNTPC(copy.fNTPC),
  fMinTPC(copy.fMinTPC),
  fMaxTPC(copy.fMaxTPC),
  fOutList(0x0),
  fHEvents(0x0),
  fTOFESD(copy.fTOFESD),
  fTOFSigma(copy.fTOFSigma),
  fTOFmaker(0x0),
  fTOFSettings(copy.fTOFSettings)
{
//
// Copy constructor
//
}


AliRsnAnalysisPhi900GeV& AliRsnAnalysisPhi900GeV::operator=(const AliRsnAnalysisPhi900GeV& copy)
{
//
// Assignment operator
//

  fUseMC = copy.fUseMC;

  fDCAr = copy.fDCAr;
  fDCAz = copy.fDCAz;
  fChi2 = copy.fChi2;
  fNTPC = copy.fNTPC;

  fMinTPC = copy.fMinTPC;
  fMaxTPC = copy.fMaxTPC;

  fTOFESD      = copy.fTOFESD;
  fTOFSigma    = copy.fTOFSigma;
  fTOFSettings = copy.fTOFSettings;

  return (*this);
}


AliRsnAnalysisPhi900GeV::~AliRsnAnalysisPhi900GeV()
{
//
// Destructor
//

  if (fOutTree[0]) delete fOutTree[0];
  if (fOutTree[1]) delete fOutTree[1];
}


void AliRsnAnalysisPhi900GeV::UserCreateOutputObjects()
{
//
// Create the output data container
//

  // setup TOF maker
  fTOFmaker = new AliRsnTOFT0maker;
  fTOFmaker->SetTimeResolution(fTOFSigma * 1E-12);
  fTOFmaker->SetESDdata(fTOFESD);
  fTOFmaker->fSettings = fTOFSettings;
  AliInfo(Form("TOF sigma    = %f", fTOFSigma));
  AliInfo(Form("TOF ESD      = %s", (fTOFESD ? "YES" : "NO")));
  AliInfo(Form("TOF settings = %s", fTOFmaker->Settings().Data()));
  
  // load dead channel map
  fTOFmaker->LoadChannelMap("tofmap.root");
  fTOFmaker->SetMaskOffChannel();
  
  // initialize random
  gRandom->SetSeed(0);

  // create output trees
  OpenFile(1);
  fOutTree[0] = new TTree("rsnTree", "Pairs");

  fOutTree[0]->Branch("pdg", &fPDG, "pdg/S");
  fOutTree[0]->Branch("im" , &fIM , "im/F" );
  fOutTree[0]->Branch("y"  , &fY  , "y/F"  );
  fOutTree[0]->Branch("pt" , &fPt , "pt/F" );
  fOutTree[0]->Branch("eta", &fEta, "eta/F");

  OpenFile(2);
  fOutTree[1] = new TTree("rsnTrue", "True pairs");

  fOutTree[1]->Branch("im" , &fIM , "im/F" );
  fOutTree[1]->Branch("y"  , &fY  , "y/F"  );
  fOutTree[1]->Branch("pt" , &fPt , "pt/F" );
  fOutTree[1]->Branch("eta", &fEta, "eta/F");

  OpenFile(3);
  fOutList = new TList;
  fHEvents = new TH1I("hEvents", "Event details", 4, 0, 4);
  fOutList->Add(fHEvents);
}


void AliRsnAnalysisPhi900GeV::UserExec(Option_t *)
{
//
// Main execution function.
// Fills the fHEvents data member with the following legenda:
// 0 -- event OK, prim vertex with tracks
// 1 -- event OK, prim vertex with SPD
// 2 -- event OK but vz large
// 3 -- event bad
//

  // retrieve ESD event and related stack (if available)
  AliESDEvent *esd   = dynamic_cast<AliESDEvent*>(fInputEvent);
  AliStack    *stack = (fMCEvent ? fMCEvent->Stack() : 0x0);
  
  // get the best primary vertex:
  // first try the one with tracks
  Int_t type = 0;
  const AliESDVertex *v = esd->GetPrimaryVertexTracks();
  if(v->GetNContributors() < 1)
  {
    // if not good, try SPD vertex
    type = 1;
    v = esd->GetPrimaryVertexSPD();
    
    // if this is not good skip this event
    if (v->GetNContributors() < 1)
    {
      fHEvents->Fill(3);
      PostData(3, fOutList);
      return;
    }
  }

  // if the Z position is larger than 10, skip this event
  if (TMath::Abs(v->GetZv()) > 10.0)
  {
    fHEvents->Fill(2);
    PostData(3, fOutList);
    return;
  }

  // use the type to fill the histogram
  fHEvents->Fill(type);

  // smear TOF times in case of MC
  if (stack) RemakeTOFtimeMC(esd);

  // get time zero for TOF
  Double_t *tof = fTOFmaker->RemakePID(esd);

  ProcessESD(esd, v, tof[0], stack);
  ProcessMC(stack);

  PostData(3, fOutList);
}


void AliRsnAnalysisPhi900GeV::Terminate(Option_t *)
{
//
// Terminate
//
}


void AliRsnAnalysisPhi900GeV::ProcessESD
(AliESDEvent *esd, const AliESDVertex *v, Double_t time0, AliStack *stack)
{
//
// This function works with the ESD object
//

  // prepare to look on all tracks to select the ones
  // which pass all the cuts
  Int_t   ntracks = esd->GetNumberOfTracks();
  TArrayI pos(ntracks);
  TArrayI neg(ntracks);

  // define fixed functions for TOF compatibility range
  Double_t a1 = 0.01, a2 = -0.03;
  Double_t b1 = 0.25, b2 =  0.25;
  Double_t c1 = 0.05, c2 = -0.03;
  Double_t ymax, ymin;

  // loop on all tracks
  Int_t    i, charge, nSPD, npos = 0, nneg = 0;
  Float_t  chi2, b[2], bCov[3];
  Double_t tpc, bb, mom, tofTime, tofRef, tofRel, times[10];
  Bool_t   okTOF;
  for (i = 0; i < ntracks; i++)
  {
    AliESDtrack *track = esd->GetTrack(i);
    if (!track) continue;

    // skip if it has not the required flags
    if (!track->IsOn(AliESDtrack::kTPCin)) continue;
    if (!track->IsOn(AliESDtrack::kTPCrefit)) continue;
    if (!track->IsOn(AliESDtrack::kITSrefit)) continue;

    // skip if it has not the TPC inner wall projection
    if (!track->GetInnerParam()) continue;
    
    // skip kink daughters
    if ((Int_t)track->GetKinkIndex(0) > 0) continue;

    // check clusters in TPC
    if (track->GetTPCclusters(0) < fNTPC) continue;

    // check chi2
    chi2  = (Float_t)track->GetTPCchi2();
    chi2 /= (Float_t)track->GetTPCclusters(0);
    if (chi2 > fChi2) continue;

    // check that has at least 1 SPD cluster
    nSPD = 0;
    if (track->HasPointOnITSLayer(0)) nSPD++;
    if (track->HasPointOnITSLayer(1)) nSPD++;
    if (nSPD < 1) continue;

    // check primary by reverting to vertex
    // and checking DCA
    if (!track->RelateToVertex(v, esd->GetMagneticField(), kVeryBig)) continue;
    track->GetImpactParameters(b, bCov);
    if (b[0] > fDCAr) continue;
    if (b[1] > fDCAz) continue;

    // check TPC dE/dx
    AliExternalTrackParam trackIn(*track->GetInnerParam());
    mom = trackIn.P();
    tpc = (Double_t)track->GetTPCsignal();
    bb  = AlephBB(mom);
    tpc = (tpc - bb) / bb;
    if (tpc < fMinTPC || tpc > fMaxTPC) continue;

    // if possible, check TOF
    okTOF = kTRUE;
    if (track->IsOn(AliESDtrack::kTOFpid))
    {
      mom = track->P();
      if (mom <= 0.26)
        okTOF = kTRUE;
      else
      {
        track->GetIntegratedTimes(times);
        tofTime = (Double_t)track->GetTOFsignal() - time0;
        tofRef  = times[AliPID::kKaon];
        tofRel  = (tofTime - tofRef) / tofRef;
        ymax    = a1 / (mom - b1) + c1;
        ymin    = a2 / (mom - b2) + c2;
        okTOF   = (tofRel >= ymin && tofRel <= ymax);
      }
    }
    if (!okTOF) continue;

    // if we arrive here, all cuts were passed
    // and we add the track to one array depending on charge
    charge = (Int_t)track->Charge();
    if (charge > 0)
      pos[npos++] = i;
    else if (charge < 0)
      neg[nneg++] = i;
  }
  
  // resize arrays accordingly
  pos.Set(npos);
  neg.Set(nneg);

  // loop to compute invariant mass
  Int_t           ip, in, lp, ln;
  AliPID          pid;
  Double_t        kmass = pid.ParticleMass(AliPID::kKaon);
  Double_t        phimass = 1.019455;
  TParticle      *partp = 0x0, *partn = 0x0;
  AliESDtrack    *tp = 0x0, *tn = 0x0;
  TLorentzVector  vp, vn, vsum, vref;
  for (ip = 0; ip < npos; ip++)
  {
    tp = esd->GetTrack(pos[ip]);
    lp = TMath::Abs(tp->GetLabel());
    if (stack) partp = stack->Particle(lp);

    for (in = 0; in < nneg; in++)
    {
      if (pos[ip] == neg[in]) continue;
      tn = esd->GetTrack(neg[in]);
      ln = TMath::Abs(tn->GetLabel());
      if (stack) partn = stack->Particle(ln);

      fPDG = 0;
      if (partp && partn)
      {
        if (partp->GetFirstMother() == partn->GetFirstMother())
        {
          if (partp->GetFirstMother() > 0)
          {
            TParticle *mum = stack->Particle(partp->GetFirstMother());
            fPDG = mum->GetPdgCode();
          }
        }
      }
      fPDG = TMath::Abs(fPDG);

      vp.SetXYZM(tp->Px(), tp->Py(), tp->Pz(), kmass);
      vn.SetXYZM(tn->Px(), tn->Py(), tn->Pz(), kmass);
      vsum = vp + vn;
      vref.SetXYZM(vsum.X(), vsum.Y(), vsum.Z(), phimass);

      fIM  = (Float_t)vsum.M();
      fPt  = (Float_t)vsum.Perp();
      fEta = (Float_t)vsum.Eta();
      fY   = (Float_t)vref.Rapidity();

      fOutTree[0]->Fill();
    }
  }

  PostData(1, fOutTree[0]);
}


void AliRsnAnalysisPhi900GeV::ProcessMC(AliStack *stack)
{
//
// Function to process stack only
//

  if (!stack) return;
  Int_t nPart = stack->GetNtrack();

  // loop to compute invariant mass
  Int_t           ip, in;
  AliPID          pid;
  Double_t        kmass = pid.ParticleMass(AliPID::kKaon);
  Double_t        phimass = 1.019455;
  TParticle      *partp = 0x0, *partn = 0x0;
  TLorentzVector  vp, vn, vsum, vref;

  for (ip = 0; ip < nPart; ip++)
  {
    partp = stack->Particle(ip);
    if (partp->GetPdgCode() != 321) continue;

    for (in = 0; in < nPart; in++)
    {
      partn = stack->Particle(in);
      if (partn->GetPdgCode() != -321) continue;

      fPDG = 0;
      if (partp->GetFirstMother() == partn->GetFirstMother())
      {
        if (partp->GetFirstMother() > 0)
        {
          TParticle *mum = stack->Particle(partp->GetFirstMother());
          fPDG = mum->GetPdgCode();
        }
      }
      fPDG = TMath::Abs(fPDG);
      if (fPDG != 333) continue;

      vp.SetXYZM(partp->Px(), partp->Py(), partp->Pz(), kmass);
      vn.SetXYZM(partn->Px(), partn->Py(), partn->Pz(), kmass);
      vsum = vp + vn;
      vref.SetXYZM(vsum.X(), vsum.Y(), vsum.Z(), phimass);

      fIM  = (Float_t)vsum.M();
      fPt  = (Float_t)vsum.Perp();
      fEta = (Float_t)vsum.Eta();
      fY   = (Float_t)vref.Rapidity();

      fOutTree[1]->Fill();
    }
  }

  PostData(2, fOutTree[1]);
}


void AliRsnAnalysisPhi900GeV::SetTPCparams(Bool_t isMC)
{
//
// Set TPC bethe-bloch parameters
//

  if (!isMC)
  {
    fTPCpar[0] = 1.41543;
    fTPCpar[1] = 2.63394E1;
    fTPCpar[2] = 5.0411E-11;
    fTPCpar[3] = 2.12543;
    fTPCpar[4] = 4.88663;
  }
  else
  {
    fTPCpar[0] = 2.15898;
    fTPCpar[1] = 1.75295E1;
    fTPCpar[2] = 3.40030E-9;
    fTPCpar[3] = 1.96178;
    fTPCpar[4] = 3.91720;
  }
}


Double_t AliRsnAnalysisPhi900GeV::AlephBB(Double_t p, Double_t mass) const
{
//
// Compute expected Bethe-Bloch for that momentum and mass
//

  if (mass < 1E-6) return 0.0;

  Double_t aa, bb, out, beta, bg;
  
  bg   = p / mass;
  beta = bg / TMath::Sqrt(1.0 + bg * bg);
  aa   = TMath::Power(beta, fTPCpar[3]);
  bb   = TMath::Power(1./bg, fTPCpar[4]);
  bb   = TMath::Log(fTPCpar[2] + bb);
  out  = (fTPCpar[1] - aa - bb) * fTPCpar[0]/aa;

  return out;
}


Double_t AliRsnAnalysisPhi900GeV::RemakeTOFtimeMC(AliESDEvent *& esd)
{
//
// Smear initial time for TOF in order to reproduce data resolution
//

  Double_t t0 = gRandom->Gaus(0,135.); //spread in ps
  Int_t ntracks = esd->GetNumberOfTracks();
  Double_t sigmaStandard = 80.; // 80 ps from TOF
  
  while (ntracks--) 
  {
    AliESDtrack *t = esd->GetTrack(ntracks);
    if ((t->GetStatus()&AliESDtrack::kTOFout) == 0 || (t->GetStatus()&AliESDtrack::kTIME)==0) continue;
    Double_t time = t->GetTOFsignal();
    if(fTOFSigma > sigmaStandard)
    {
      Double_t sigmaAdded = TMath::Sqrt(fTOFSigma*fTOFSigma - sigmaStandard*sigmaStandard);
      Double_t timerandomtrack = gRandom->Gaus(0, sigmaAdded); //spread in ps
      time += timerandomtrack;
   }
   time += t0;
   t->SetTOFsignal(time);
 }
 return t0;
} 
