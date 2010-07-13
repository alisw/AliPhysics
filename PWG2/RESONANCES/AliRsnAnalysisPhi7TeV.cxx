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
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"

#include "AliRsnAnalysisPhi7TeV.h"

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::AliRsnAnalysisPhi7TeV(const char *name) :
  AliAnalysisTaskSE(name),
  fUseMC(kFALSE),
  fPDG(0),
  fIM(0.0),
  fPt(0.0),
  fY(0.0),
  fEta(0.0),
  fMaxDCAr(1E6),
  fMaxDCAz(1E6),
  fMaxChi2(1E6),
  fMinNTPC(0),
  fMinTPCband(-1E6),
  fMaxTPCband( 1E6),
  fRsnTreeComp(0x0),
  fRsnTreeTrue(0x0),
  fOutList(0x0),
  fHEvents(0x0),
  fHCuts(0x0),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(kFALSE),
  fTOFcorrectTExp(kFALSE),
  fTOFuseT0(kFALSE),
  fTOFtuneMC(kFALSE),
  fTOFresolution(0.0)
  
{
//
// Constructor
//

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TList::Class());
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::AliRsnAnalysisPhi7TeV(const AliRsnAnalysisPhi7TeV& copy) :
  AliAnalysisTaskSE(copy),
  fUseMC(copy.fUseMC),
  fPDG(0),
  fIM(0.0),
  fPt(0.0),
  fY(0.0),
  fEta(0.0),
  fMaxDCAr(copy.fMaxDCAr),
  fMaxDCAz(copy.fMaxDCAz),
  fMaxChi2(copy.fMaxChi2),
  fMinNTPC(copy.fMinNTPC),
  fMinTPCband(copy.fMinTPCband),
  fMaxTPCband(copy.fMaxTPCband),
  fRsnTreeComp(0x0),
  fRsnTreeTrue(0x0),
  fOutList(0x0),
  fHEvents(0x0),
  fHCuts(0x0),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(kFALSE),
  fTOFcorrectTExp(kFALSE),
  fTOFuseT0(kFALSE),
  fTOFtuneMC(kFALSE),
  fTOFresolution(0.0)
{
//
// Copy constructor
//
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV& AliRsnAnalysisPhi7TeV::operator=(const AliRsnAnalysisPhi7TeV& copy)
{
//
// Assignment operator
//

  fUseMC = copy.fUseMC;

  fMaxDCAr = copy.fMaxDCAr;
  fMaxDCAz = copy.fMaxDCAz;
  fMaxChi2 = copy.fMaxChi2;
  fMinNTPC = copy.fMinNTPC;

  fMinTPCband = copy.fMinTPCband;
  fMaxTPCband = copy.fMaxTPCband;
  
  fTOFcalibrateESD = copy.fTOFcalibrateESD;
  fTOFcorrectTExp = copy.fTOFcorrectTExp;
  fTOFuseT0 = copy.fTOFuseT0;
  fTOFtuneMC = copy.fTOFtuneMC;
  fTOFresolution = copy.fTOFresolution;

  return (*this);
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::~AliRsnAnalysisPhi7TeV()
{
//
// Destructor
//

  if (fRsnTreeComp) delete fRsnTreeComp;
  if (fRsnTreeTrue) delete fRsnTreeTrue;
  if (fHEvents)     delete fHEvents;
  if (fHCuts)       delete fHCuts;
  if (fESDpid)      delete fESDpid;
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::UserCreateOutputObjects()
{
//
// Create the output data container
//

  // setup TPC response
  fESDpid = new AliESDpid;
  fESDpid->GetTPCResponse().SetBetheBlochParameters(fTPCpar[0],fTPCpar[1],fTPCpar[2],fTPCpar[3],fTPCpar[4]);

  // setup TOF maker & calibration
  fTOFcalib = new AliTOFcalib;
  fTOFmaker = new AliTOFT0maker(fESDpid, fTOFcalib);
  fTOFmaker->SetTimeResolution(fTOFresolution);
  
  // initialize random
  gRandom->SetSeed(0);

  // create output trees
  OpenFile(1);
  fRsnTreeComp = new TTree("rsnTree", "Pairs");

  fRsnTreeComp->Branch("pdg", &fPDG, "pdg/S");
  fRsnTreeComp->Branch("im" , &fIM , "im/F" );
  fRsnTreeComp->Branch("y"  , &fY  , "y/F"  );
  fRsnTreeComp->Branch("pt" , &fPt , "pt/F" );
  fRsnTreeComp->Branch("eta", &fEta, "eta/F");

  OpenFile(2);
  fRsnTreeTrue = new TTree("rsnTrue", "True pairs");

  fRsnTreeTrue->Branch("im" , &fIM , "im/F" );
  fRsnTreeTrue->Branch("y"  , &fY  , "y/F"  );
  fRsnTreeTrue->Branch("pt" , &fPt , "pt/F" );
  fRsnTreeTrue->Branch("eta", &fEta, "eta/F");

  OpenFile(3);
  fOutList = new TList;
  fHEvents = new TH1I("hEvents", "Event details", 4, 0, 4);
  fHCuts   = new TH1I("hCuts", "Cuts not passed", 11, 0, 11);
  fHCuts->GetXaxis()->SetBinLabel( 1, "Flags");
  fHCuts->GetXaxis()->SetBinLabel( 2, "No TPC inner");
  fHCuts->GetXaxis()->SetBinLabel( 3, "Kink daughter");
  fHCuts->GetXaxis()->SetBinLabel( 4, "Few clusters in TPC");
  fHCuts->GetXaxis()->SetBinLabel( 5, "Large #chi^{2}");
  fHCuts->GetXaxis()->SetBinLabel( 6, "No SPD clusters");
  fHCuts->GetXaxis()->SetBinLabel( 7, "Failed revert to vertex");
  fHCuts->GetXaxis()->SetBinLabel( 8, "Too large DCA");
  fHCuts->GetXaxis()->SetBinLabel( 9, "Failed TPC PID");
  fHCuts->GetXaxis()->SetBinLabel(10, "Failed TOF PID");
  fHCuts->GetXaxis()->SetBinLabel(11, "Track OK");
  fOutList->Add(fHEvents);
  fOutList->Add(fHCuts);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::UserExec(Option_t *)
{
//
// Main execution function.
// Fills the fHEvents data member with the following legenda:
// 0 -- event OK, prim vertex with tracks
// 1 -- event OK, prim vertex with SPD
// 2 -- event OK but vz large
// 3 -- event bad
//

  static Int_t evNum = 0;
  evNum++;

  // retrieve ESD event and related stack (if available)
  AliESDEvent *esd   = dynamic_cast<AliESDEvent*>(fInputEvent);
  AliStack    *stack = (fMCEvent ? fMCEvent->Stack() : 0x0);
  //cout << "NTRACKS: " << esd->GetNumberOfTracks() << endl;
  
  // get the best primary vertex:
  // first try the one with tracks
  Int_t type = 0;
  const AliESDVertex *v = esd->GetPrimaryVertexTracks();
  //cout << "[ev " << evNum << "] vertex tracks: " << v << " contrib = " << v->GetNContributors() << " -- status = " << v->GetStatus() << endl;
  if(v->GetNContributors() < 1)
  {
    // if not good, try SPD vertex
    type = 1;
    v = esd->GetPrimaryVertexSPD();
    //cout << "[ev " << evNum << "] vertex SPD: " << v << " contrib = " << v->GetNContributors() << " -- status = " << v->GetStatus() << endl;
    
    // if this is not good skip this event
    if (v->GetNContributors() < 1)
    {
      fHEvents->Fill(3);
      PostData(3, fOutList);
      return;
    }
  }

  // if the Z position is larger than 10, skip this event
  //cout << "[ev " << evNum << "] vertex Z = " << v->GetZv() << endl;
  if (TMath::Abs(v->GetZv()) > 10.0)
  {
    fHEvents->Fill(2);
    PostData(3, fOutList);
    return;
  }
  
  //cout << "EVENT " << evNum << " is OK" << endl;

  // use the type to fill the histogram
  fHEvents->Fill(type);

  ProcessESD(esd, v, stack);
  ProcessMC(stack);

  PostData(3, fOutList);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::Terminate(Option_t *)
{
//
// Terminate
//
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::ProcessESD
(AliESDEvent *esd, const AliESDVertex *v, AliStack *stack)
{
//
// This function works with the ESD object
//

  // TOF stuff #1: init OCDB
  Int_t run = esd->GetRunNumber();
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(run);
  // TOF stuff #2: init calibration
  fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
  fTOFcalib->Init();
  // TOF stuff #3: calibrate
  if (fTOFcalibrateESD) fTOFcalib->CalibrateESD(esd);
  if (fTOFtuneMC) fTOFmaker->TuneForMC(esd);
  if (fTOFuseT0) 
  {
    fTOFmaker->ComputeT0TOF(esd);
    fTOFmaker->ApplyT0TOF(esd);
    fESDpid->MakePID(esd, kFALSE, 0.);
  }

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
  Int_t    i, icut, charge, nSPD, npos = 0, nneg = 0;
  Float_t  chi2, b[2], bCov[3];
  Double_t mom, tpcNSigma, tpcMaxNSigma, tofTime, tofSigma, tofRef, tofRel, times[10];
  Bool_t   okTOF;
  for (i = 0; i < ntracks; i++)
  {
    AliESDtrack *track = esd->GetTrack(i);
    if (!track) {fHCuts->Fill(icut); continue;}

    // skip if it has not the required flags
    icut = 0;
    if (!track->IsOn(AliESDtrack::kTPCin)) {fHCuts->Fill(icut); continue;}
    if (!track->IsOn(AliESDtrack::kTPCrefit)) {fHCuts->Fill(icut); continue;}
    if (!track->IsOn(AliESDtrack::kITSrefit)) {fHCuts->Fill(icut); continue;}

    // skip if it has not the TPC inner wall projection
    icut = 1;
    if (!track->GetInnerParam()) {fHCuts->Fill(icut); continue;}
    
    // skip kink daughters
    icut = 2;
    if ((Int_t)track->GetKinkIndex(0) > 0) {fHCuts->Fill(icut); continue;}

    // check clusters in TPC
    icut = 3;
    if (track->GetTPCclusters(0) < fMinNTPC) {fHCuts->Fill(icut); continue;}

    // check chi2
    icut = 4;
    chi2  = (Float_t)track->GetTPCchi2();
    chi2 /= (Float_t)track->GetTPCclusters(0);
    if (chi2 > fMaxChi2) {fHCuts->Fill(icut); continue;}

    // check that has at least 1 SPD cluster
    icut = 5;
    nSPD = 0;
    if (track->HasPointOnITSLayer(0)) nSPD++;
    if (track->HasPointOnITSLayer(1)) nSPD++;
    if (nSPD < 1) {fHCuts->Fill(icut); continue;}

    // check primary by reverting to vertex
    // and checking DCA
    icut = 6;
    if (!track->RelateToVertex(v, esd->GetMagneticField(), kVeryBig)) {fHCuts->Fill(icut); continue;}
    icut = 7;
    track->GetImpactParameters(b, bCov);
    if (b[0] > fMaxDCAr) {fHCuts->Fill(icut); continue;}
    if (b[1] > fMaxDCAz) {fHCuts->Fill(icut); continue;}

    // check TPC dE/dx
    icut = 8;
    //cout << "TPC SIGMA: " << fESDpid->NumberOfSigmasTPC(track, AliPID::kKaon) << endl;
    tpcNSigma = TMath::Abs(fESDpid->NumberOfSigmasTPC(track, AliPID::kKaon));
    if (track->GetInnerParam()->P() > fTPCpLimit) tpcMaxNSigma = fMinTPCband; else tpcMaxNSigma = fMaxTPCband;
    if (tpcNSigma > tpcMaxNSigma) {fHCuts->Fill(icut); continue;}

    // if possible, check TOF
    icut = 9;
    okTOF = kTRUE;
    if (track->IsOn(AliESDtrack::kTOFout | AliESDtrack::kTIME))
    {
      mom = track->P();
      if (mom <= 0.26)
        okTOF = kTRUE;
      else
      {
        track->GetIntegratedTimes(times);
        tofTime  = (Double_t)track->GetTOFsignal();
        tofSigma = fTOFmaker->GetExpectedSigma(mom, times[AliPID::kKaon], AliPID::ParticleMass(AliPID::kKaon));
        tofRef   = times[AliPID::kKaon];
        tofRel   = (tofTime - tofRef) / tofRef;
        ymax     = a1 / (mom - b1) + c1;
        ymin     = a2 / (mom - b2) + c2;
        okTOF    = (tofRel >= ymin && tofRel <= ymax);
      }
    }
    if (!okTOF) {fHCuts->Fill(icut); continue;}

    // if we arrive here, all cuts were passed
    // and we add the track to one array depending on charge
    icut = 10;
    charge = (Int_t)track->Charge();
    if (charge > 0)
      pos[npos++] = i;
    else if (charge < 0)
      neg[nneg++] = i;
    
    // fill the bin corresponding to passed cuts
    fHCuts->Fill(icut);
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

      fRsnTreeComp->Fill();
    }
  }

  PostData(1, fRsnTreeComp);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::ProcessMC(AliStack *stack)
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

      fRsnTreeTrue->Fill();
    }
  }

  PostData(2, fRsnTreeTrue);
}
