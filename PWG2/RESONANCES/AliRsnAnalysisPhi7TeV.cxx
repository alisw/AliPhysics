//
// Implementation file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#include "Riostream.h"
#include <iomanip>

#include "TH1.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "AliLog.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include "AliITSPIDResponse.h"

#include "AliRsnEvent.h"
#include "AliRsnTarget.h"
#include "AliRsnDaughter.h"
#include "AliRsnCutESD2010.h"

#include "AliRsnAnalysisPhi7TeV.h"

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::AliRsnAnalysisPhi7TeV(const char *name, Bool_t isMC) :
  AliAnalysisTaskSE(name),
  fUseMC(kFALSE),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fAddITSSA(kFALSE),
  fMaxVz(1E6),
  fMaxITSband(1E6),
  fMaxITSmom(0.0),
  fTPCpLimit(0.35),
  fMinTPCband(-1E6),
  fMaxTPCband( 1E6),
  fMinTOF(-3.0),
  fMaxTOF( 3.0),
  fOutList(0x0),
  fUnlike(0x0),
  fLikePP(0x0),
  fLikeMM(0x0),
  fTrues(0x0),
  fHEvents(0x0),
  fESDtrackCutsTPC(),
  fESDtrackCutsITS(),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(!isMC),
  fTOFcorrectTExp(kTRUE),
  fTOFuseT0(kTRUE),
  fTOFtuneMC(isMC),
  fTOFresolution(100.0),
  fDaughter(),
  fRsnCuts()
{
//
// Constructor
//

  SetUseMC(isMC);
  DefineOutput(1, TList::Class());
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::AliRsnAnalysisPhi7TeV(const AliRsnAnalysisPhi7TeV& copy) :
  AliAnalysisTaskSE(copy),
  fUseMC(copy.fUseMC),
  fCheckITS(copy.fCheckITS),
  fCheckTPC(copy.fCheckTPC),
  fCheckTOF(copy.fCheckTOF),
  fAddITSSA(copy.fAddITSSA),
  fMaxVz(copy.fMaxVz),
  fMaxITSband(copy.fMaxITSband),
  fMaxITSmom(copy.fMaxITSmom),
  fTPCpLimit(copy.fTPCpLimit),
  fMinTPCband(copy.fMinTPCband),
  fMaxTPCband(copy.fMaxTPCband),
  fMinTOF(copy.fMinTOF),
  fMaxTOF(copy.fMaxTOF),
  fOutList(0x0),
  fUnlike(0x0),
  fLikePP(0x0),
  fLikeMM(0x0),
  fTrues(0x0),
  fHEvents(0x0),
  fESDtrackCutsTPC(copy.fESDtrackCutsTPC),
  fESDtrackCutsITS(copy.fESDtrackCutsITS),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(kFALSE),
  fTOFcorrectTExp(kFALSE),
  fTOFuseT0(kFALSE),
  fTOFtuneMC(kFALSE),
  fTOFresolution(0.0),
  fDaughter(),
  fRsnCuts()
{
//
// Copy constructor
//

  SetUseMC(copy.fUseMC);
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV& AliRsnAnalysisPhi7TeV::operator=(const AliRsnAnalysisPhi7TeV& copy)
{
//
// Assignment operator
//

  fUseMC = copy.fUseMC;
  fCheckITS = copy.fCheckITS;
  fCheckTPC = copy.fCheckTPC;
  fCheckTOF = copy.fCheckTOF;
  fAddITSSA = copy.fAddITSSA;

  fMaxVz   = copy.fMaxVz;
  fMaxITSband = copy.fMaxITSband;
  fMaxITSmom  = copy.fMaxITSmom;
  
  fTPCpLimit  = copy.fTPCpLimit;
  fMinTPCband = copy.fMinTPCband;
  fMaxTPCband = copy.fMaxTPCband;
  
  fESDtrackCutsTPC = copy.fESDtrackCutsTPC;
  fESDtrackCutsITS = copy.fESDtrackCutsITS;
  
  fTOFcalibrateESD = copy.fTOFcalibrateESD;
  fTOFcorrectTExp = copy.fTOFcorrectTExp;
  fTOFuseT0 = copy.fTOFuseT0;
  fTOFtuneMC = copy.fTOFtuneMC;
  fTOFresolution = copy.fTOFresolution;
  
  SetUseMC(copy.fUseMC);

  return (*this);
}

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::~AliRsnAnalysisPhi7TeV()
{
//
// Destructor
//
}

//_________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::SetUseMC(Bool_t isMC)
{
//
// Sets some aspects of cuts depending on the fact that runs on MC or not
//

  fUseMC = isMC;
  
  if (isMC)
  {
    SetTPCpar(2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720);
    SetTOFcalibrateESD(kFALSE);
    SetTOFtuneMC(kTRUE);
  }
  else
  {
    SetTPCpar(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);
    SetTOFcalibrateESD(kTRUE);
    SetTOFtuneMC(kFALSE);
  }
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::UserCreateOutputObjects()
{
//
// Create the output data container
//
  
  // initialize random
  //gRandom->SetSeed(0);

  // create output trees
  OpenFile(1);
  fOutList    = new TList;
  fHEvents    = new TH1I("hEvents", "Event details", 6, 0, 6);
  fVertexX[0] = new TH1F("hVertexTracksX", "X position of primary vertex (tracks)", 200,  -2,  2);
  fVertexY[0] = new TH1F("hVertexTracksY", "Y position of primary vertex (tracks)", 200,  -2,  2);
  fVertexZ[0] = new TH1F("hVertexTracksZ", "Z position of primary vertex (tracks)", 400, -40, 40);
  fVertexX[1] = new TH1F("hVertexSPDX", "X position of primary vertex (SPD)", 200,  -2,  2);
  fVertexY[1] = new TH1F("hVertexSPDY", "Y position of primary vertex (SPD)", 200,  -2,  2);
  fVertexZ[1] = new TH1F("hVertexSPDZ", "Z position of primary vertex (SPD)", 400, -40, 40);
  //fUnlike     = new TH3F("hPM", "+-   pairs", 500, 0.9, 1.4, 100, 0.0, 10.0, 24, -1.2, 1.2);
  //fLikePP     = new TH3F("hPP", "++   pairs", 500, 0.9, 1.4, 100, 0.0, 10.0, 24, -1.2, 1.2);
  //fLikeMM     = new TH3F("hMM", "--   pairs", 500, 0.9, 1.4, 100, 0.0, 10.0, 24, -1.2, 1.2);
  //fTrues      = new TH3F("hTR", "True pairs", 500, 0.9, 1.4, 100, 0.0, 10.0, 24, -1.2, 1.2);
  TFile *ffile = TFile::Open("template.root");
  if (ffile)
  {
    TH3F *tmp = (TH3F*)ffile->Get("template");
    fUnlike = (TH3F*)tmp->Clone("hPM");
    fLikePP = (TH3F*)tmp->Clone("hPP");
    fLikeMM = (TH3F*)tmp->Clone("hMM");
    fTrues  = (TH3F*)tmp->Clone("hTR");
  }
  
  fHEvents->GetXaxis()->SetBinLabel(1, "Good vertex with tracks");
  fHEvents->GetXaxis()->SetBinLabel(2, "Good vertex with SPD");
  fHEvents->GetXaxis()->SetBinLabel(3, "Far vertex with tracks");
  fHEvents->GetXaxis()->SetBinLabel(4, "Far vertex with SPD");
  fHEvents->GetXaxis()->SetBinLabel(5, "No good vertex");
  fHEvents->GetXaxis()->SetBinLabel(6, "Empty event");

  fOutList->Add(fHEvents);
  fOutList->Add(fVertexX[0]);
  fOutList->Add(fVertexY[0]);
  fOutList->Add(fVertexZ[0]);
  fOutList->Add(fVertexX[1]);
  fOutList->Add(fVertexY[1]);
  fOutList->Add(fVertexZ[1]);
  fOutList->Add(fUnlike);
  fOutList->Add(fLikePP);
  fOutList->Add(fLikeMM);
  fOutList->Add(fTrues);
  
  // setup RSN-related objects
  
  fRsnCuts.SetMC       (fUseMC);
  fRsnCuts.SetCheckITS (fCheckITS);
  fRsnCuts.SetCheckTPC (fCheckTPC);
  fRsnCuts.SetCheckTOF (fCheckTOF);
  fRsnCuts.SetUseITSTPC(kTRUE);
  fRsnCuts.SetUseITSSA (fAddITSSA);
  fRsnCuts.SetPID      (AliPID::kKaon);
  fRsnCuts.SetMaxITSPIDmom(fMaxITSmom);
  fRsnCuts.SetITSband(fMaxITSband);
  fRsnCuts.SetTPCpLimit(fTPCpLimit);
  fRsnCuts.SetTPCrange(fMinTPCband, fMaxTPCband);
  fRsnCuts.SetTPCpar(fTPCpar[0],fTPCpar[1],fTPCpar[2],fTPCpar[3],fTPCpar[4]);
  //fRsnCuts.SetTOFcalibrateESD(fTOFcalibrateESD);
  fRsnCuts.SetTOFcorrectTExp (fTOFcorrectTExp);
  fRsnCuts.SetTOFuseT0       (fTOFuseT0);
  fRsnCuts.SetTOFtuneMC      (fTOFtuneMC);
  fRsnCuts.SetTOFresolution  (fTOFresolution);
  fRsnCuts.SetTOFrange       (fMinTOF, fMaxTOF);
  fRsnCuts.GetCutsTPC()->Copy(fESDtrackCutsTPC);
  fRsnCuts.GetCutsITS()->Copy(fESDtrackCutsITS);
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

  // retrieve ESD event and related stack (if available)
  AliESDEvent *esd   = dynamic_cast<AliESDEvent*>(fInputEvent);
  AliStack    *stack = (fMCEvent ? fMCEvent->Stack() : 0x0);
  
  // check the event
  Int_t eval = EventEval(esd);
  fHEvents->Fill(eval);
  
  // if the event is good for analysis, process it
  if (eval == kGoodTracksPrimaryVertex || eval == kGoodSPDPrimaryVertex)
  {
    ProcessESD(esd, stack);
  }
  
  // update histogram container
  PostData(1, fOutList);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::Terminate(Option_t *)
{
//
// Terminate
//
}

//__________________________________________________________________________________________________
Int_t AliRsnAnalysisPhi7TeV::EventEval(AliESDEvent *esd)
{
//
// Checks if the event is good for analysis.
// Returns one of the flag values defined in the header
//

  static Int_t evNum = 0;
  evNum++;
  
  // reject empty events
  Int_t ntracks = esd->GetNumberOfTracks();
  if (!ntracks) return kEmptyEvent;
  
  // get the best primary vertex:
  // first try the one with tracks
  const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
  const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
  Double_t            vzTrk = 1000000.0;
  Double_t            vzSPD = 1000000.0;
  Int_t               ncTrk = -1;
  Int_t               ncSPD = -1;
  if (vTrk) ncTrk = (Int_t)vTrk->GetNContributors();
  if (vSPD) ncSPD = (Int_t)vSPD->GetNContributors();
  if (vTrk) vzTrk = TMath::Abs(vTrk->GetZv());
  if (vSPD) vzSPD = TMath::Abs(vSPD->GetZv());
  if(vTrk && ncTrk > 0)
  {
    // fill the histograms
    fVertexX[0]->Fill(vTrk->GetXv());
    fVertexY[0]->Fill(vTrk->GetYv());
    fVertexZ[0]->Fill(vTrk->GetZv());
    
    // check VZ position
    if (vzTrk <= fMaxVz)
      return kGoodTracksPrimaryVertex;
    else
      return kFarTracksPrimaryVertex;
  }
  else if (vSPD && ncSPD > 0)
  {
    // fill the histograms
    fVertexX[1]->Fill(vSPD->GetXv());
    fVertexY[1]->Fill(vSPD->GetYv());
    fVertexZ[1]->Fill(vSPD->GetZv());
    
    // check VZ position
    if (vzSPD <= fMaxVz)
      return kGoodSPDPrimaryVertex;
    else
      return kFarSPDPrimaryVertex;
  }
  else
    return kNoGoodPrimaryVertex;
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::ProcessESD
(AliESDEvent *esd, AliStack *stack)
{
//
// This function works with the ESD object
//

  static Int_t lastRun = -1;
  static Int_t evnum = 0;
  evnum++;
  
  // get current run
  Int_t run = esd->GetRunNumber();
  
  // if absent, initialize ESD pid response
  if (!fESDpid)
  {
    AliITSPIDResponse itsresponse(fUseMC);
    
    fESDpid = new AliESDpid;
    fESDpid->GetTPCResponse().SetBetheBlochParameters(fTPCpar[0],fTPCpar[1],fTPCpar[2],fTPCpar[3],fTPCpar[4]);
    fESDpid->GetITSResponse() = itsresponse;
  }

  // if the run number has changed, update it now and give a message
  if (run != lastRun)
  {
    AliInfo("============================================================================================");
    AliInfo(Form("*** CHANGING RUN NUMBER: PREVIOUS = %d --> CURRENT = %d ***", lastRun, run));
    AliInfo("============================================================================================");
    lastRun = run;
  
    AliCDBManager::Instance()->SetDefaultStorage("raw://");
    AliCDBManager::Instance()->SetRun(lastRun);
    
    if (fTOFmaker) delete fTOFmaker;
    if (fTOFcalib) delete fTOFcalib;
    
    fTOFcalib = new AliTOFcalib();
    if (fTOFcorrectTExp) fTOFcalib->SetCorrectTExp(kTRUE);
    fTOFcalib->Init();
    
    fTOFmaker = new AliTOFT0maker(fESDpid, fTOFcalib);
    fTOFmaker->SetTimeResolution(fTOFresolution);
  }

  if (fTOFcalibrateESD) fTOFcalib->CalibrateESD(esd);
  if (fTOFtuneMC) fTOFmaker->TuneForMC(esd);
  if (fTOFuseT0)
  {
    fTOFmaker->ComputeT0TOF(esd);
    fTOFmaker->ApplyT0TOF(esd);
    fESDpid->MakePID(esd, kFALSE, 0.);
  }
  
  // RSN event
  AliRsnEvent event;
  event.SetRef(esd);
  event.SetRefMC(fMCEvent);
  //fRsnCuts.ProcessEvent(esd);
  
  // loop on all tracks
  Int_t           i1, i2, label, pdg, ntracks = esd->GetNumberOfTracks();
  Bool_t          okTrack;
  AliPID          pid;
  Double_t        kmass = pid.ParticleMass(AliPID::kKaon);
  Double_t        phimass = 1.019455;
  Double_t        angle, cosangle;
  AliMCParticle  *p1 = 0x0, *p2 = 0x0;
  AliESDtrack    *t1 = 0x0, *t2 = 0x0;
  TLorentzVector  v1, v2, vsum, vref;
  
  // external loop (T1)
  for (i1 = 0; i1 < ntracks; i1++)
  {
    t1 = esd->GetTrack(i1);
    if (!t1) continue;
    
    // setup RSN
    fDaughter.SetRef(t1);
    if (stack)
    {
      p1 = (AliMCParticle*)fMCEvent->GetTrack(TMath::Abs(t1->GetLabel()));
      if (p1) fDaughter.SetRefMC(p1);
    }
    fDaughter.SetMass(kmass);
    v1.SetXYZM(t1->Px(), t1->Py(), t1->Pz(), kmass);
    
    // check track
    okTrack = OkTrack(t1, AliPID::kKaon);
    
    // internal loop (T2)
    for (i2 = i1+1; i2 < ntracks; i2++)
    {
      t2 = esd->GetTrack(i2);
      if (!t2) continue;
      
      // check track
      if (!OkTrack(t2, AliPID::kKaon)) continue;
        
      // if unlike pair, check if it is true
      pdg = 0;
      if ((t1->Charge() == t2->Charge()) && p1 && p2)
      {
        TParticle *part1 = p1->Particle();
        TParticle *part2 = p2->Particle();
        if (TMath::Abs(part1->GetPdgCode()) == 321 && TMath::Abs(part2->GetPdgCode()) == 321 && part1->GetFirstMother() == part2->GetFirstMother())
        {
          label = part1->GetFirstMother();
          if (label >= 0 && label < stack->GetNtrack())
          {
            TParticle *mum = stack->Particle(label);
            pdg = mum->GetPdgCode();
          }
        }
      }
      pdg = TMath::Abs(pdg);

      // combine momenta
      v1.SetXYZM(t1->Px(), t1->Py(), t1->Pz(), kmass);
      v2.SetXYZM(t2->Px(), t2->Py(), t2->Pz(), kmass);
      angle = v1.Angle(v2.Vect());
      cosangle = TMath::Abs(TMath::Cos(angle));
      if (cosangle <= 0.02) continue;
      vsum = v1 + v2;
      vref.SetXYZM(vsum.X(), vsum.Y(), vsum.Z(), phimass);

      // fill appropriate histogram
      if (t1->Charge() != t2->Charge())
      {
        fUnlike->Fill(vsum.M(), vsum.Perp(), vref.Rapidity());
        if (pdg == 333) fTrues->Fill(vsum.M(), vsum.Perp(), vref.Rapidity());
      }
      else if (t1->Charge() > 0)
      {
        fLikePP->Fill(vsum.M(), vsum.Perp(), vref.Rapidity());
      }
      else
      {
        fLikeMM->Fill(vsum.M(), vsum.Perp(), vref.Rapidity());
      }
    }
  }

  PostData(1, fOutList);
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhi7TeV::OkQuality(AliESDtrack *track)
{
//
// Check track quality parameters.
// Rejects all tracks which are not either TPC+ITS nor ITS standalone.
// If tracks of any type are not flagged to be used, they are rejected anyway.
//

  if (IsITSTPC(track)) 
    return fESDtrackCutsTPC.IsSelected(track);
  else if (IsITSSA (track))
  {
    if (fAddITSSA) 
      return fESDtrackCutsITS.IsSelected(track);
    else
      return kFALSE;
  }
  else
    return kFALSE;
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhi7TeV::OkITSPID (AliESDtrack *track, AliPID::EParticleType pid)
{
//
// Check ITS particle identification with 3sigma cut
//

  // reject not ITS standalone tracks
  if (!IsITSSA(track)) return kFALSE;
  
  // count PID layers and reject if they are too few
  Int_t   k, nITSpidLayers = 0;
  UChar_t itsCluMap = track->GetITSClusterMap();
  for(k = 2; k < 6; k++) if(itsCluMap & (1 << k)) ++nITSpidLayers;
  if (nITSpidLayers < 3)
  {
    AliDebug(AliLog::kDebug+2, "Rejecting track with too few ITS pid layers");
    return kFALSE;
  }
  
  // check the track type (ITS+TPC or ITS standalone)
  // and reject it if it is of none of the allowed types
  Bool_t isSA = kFALSE;
  if (IsITSTPC(track)) isSA = kFALSE;
  else if (IsITSSA(track)) isSA = kTRUE;
  else
  {
    AliWarning("Track is neither ITS+TPC nor ITS standalone");
    return kFALSE;
  }
  
  // create the PID response object and compute nsigma
  AliITSPIDResponse &itsrsp = fESDpid->GetITSResponse();
  Double_t mom    = track->P();
  Double_t nSigma = itsrsp.GetNumberOfSigmas(mom, track->GetITSsignal(), pid, nITSpidLayers, isSA);
  
  // evaluate the cut
  Bool_t ok = (TMath::Abs(nSigma) <= fMaxITSband);
  
  // debug message
  AliDebug(AliLog::kDebug + 2, Form("ITS nsigma = %f -- max = %f -- cut %s", nSigma, fMaxITSband, (ok ? "passed" : "failed")));
  
  // outcome
  return ok;
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhi7TeV::OkTPCPID (AliESDtrack *track, AliPID::EParticleType pid)
{
//
// Check TPC particle identification with {3|5}sigmacut,
// depending on the track total momentum.
//

  // reject not TPC tracks
  if (!IsITSTPC(track)) return kFALSE;

  // setup TPC PID response
  AliTPCPIDResponse &tpcrsp = fESDpid->GetTPCResponse();
  tpcrsp.SetBetheBlochParameters(fTPCpar[0],fTPCpar[1],fTPCpar[2],fTPCpar[3],fTPCpar[4]);
  
  // get momentum and number of sigmas and choose the reference band
  Double_t mom       = track->GetInnerParam()->P();
  Double_t nSigma    = tpcrsp.GetNumberOfSigmas(mom, track->GetTPCsignal(), track->GetTPCsignalN(), pid);
  Double_t maxNSigma = fMaxTPCband;
  if (mom < fTPCpLimit) maxNSigma = fMinTPCband;
  
  // evaluate the cut
  Bool_t ok = (TMath::Abs(nSigma) <= maxNSigma);
  
  // debug message
  AliDebug(AliLog::kDebug + 2, Form("TPC nsigma = %f -- max = %f -- cut %s", nSigma, maxNSigma, (ok ? "passed" : "failed")));
  
  // outcome
  return ok;
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhi7TeV::OkTOFPID (AliESDtrack *track, AliPID::EParticleType pid)
{
//
// Check TOF particle identification if matched there.
//

  // reject not TOF-matched tracks
  if (!MatchTOF(track)) return kFALSE;
  
  // setup TOF PID response
  AliTOFPIDResponse &tofrsp = fESDpid->GetTOFResponse();
  
  // get info for computation
  Double_t momentum = track->P();
  Double_t time     = track->GetTOFsignal();
  Double_t timeint[AliPID::kSPECIES];
  tofrsp.GetStartTime(momentum);
  track->GetIntegratedTimes(timeint);

  // check the cut
  Double_t timeDiff = time - timeint[(Int_t)pid];
  Double_t sigmaRef = tofrsp.GetExpectedSigma(momentum, timeint[(Int_t)pid], AliPID::ParticleMass(pid));
  Double_t nSigma   = timeDiff / sigmaRef;
  
  // evaluate the cut
  Bool_t ok = (nSigma >= fMinTOF && nSigma <= fMaxTOF);
  
  // debug message
  AliDebug(AliLog::kDebug + 2, Form("TOF nsigma = %f -- range = %f - %f -- cut %s", nSigma, fMinTOF, fMaxTOF, (ok ? "passed" : "failed")));
  //if (print) cout << Form("**PHI** TOF nsigma = %f -- range = %f - %f -- cut %s", nSigma, fMinTOF, fMaxTOF, (ok ? "passed" : "failed")) << endl;
  
  // outcome
  return ok;
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisPhi7TeV::OkTrack(AliESDtrack *track, AliPID::EParticleType pid)
{
//
// Global track cut check
//

  Bool_t okITS, okTPC, okTOF;

  // check quality and track type and reject tracks not passing this step
  if (!OkQuality(track))
  {
    AliDebug(AliLog::kDebug+2, "Failed quality cut");
    //printf("[PHI] Track %4d: REJECTED\n", i);
    return kFALSE;
  }
  
  // ITS PID can be checked always
  // if PID is not required, the flag is sed as
  // if the cut was alsways passed 
  okITS = OkITSPID(track, AliPID::kKaon);
  if (!fCheckITS) okITS = kTRUE;
  
  // TPC PID can be checked only for TPC+ITS tracks
  // if PID is not required, the flag is sed as
  // if the cut was alsways passed
  okTPC = kFALSE;
  if (IsITSTPC(track)) okTPC = OkTPCPID(track, AliPID::kKaon);
  if (!fCheckTPC) okTPC = kTRUE;
  
  // TOF PID can be checked only if TOF is matched
  // if PID is not required, the flag is sed as
  // if the cut was alsways passed
  okTOF = kFALSE;
  if (IsITSTPC(track) && MatchTOF(track)) okTOF = OkTOFPID(track, AliPID::kKaon);
  if (!fCheckTOF) okTOF = kTRUE;
  
  // now combine all outcomes according to the different possibilities:
  // -- ITS standalone:
  //    --> only ITS PID, always
  // -- ITS + TPC:
  //    --> ITS PID, only for momenta lower than 'fMaxITSPIDmom' and when the ITSpid flag is active
  //    --> TPC PID, always --> MASTER (first to be checked, if fails, track is rejected)
  //    --> TOF PID, only if matched
  if (IsITSSA(track))
  {
    if (!okITS)
    {
      AliDebug(AliLog::kDebug+2, "ITS standalone track --> ITS PID failed");
      //printf("[PHI] Track %4d: REJECTED\n", i);
      return kFALSE;
    }
  }
  else // checking IsITSTPC() is redundant due to OkQuality() cut check
  {
    if (!okTPC)
    {
      AliDebug(AliLog::kDebug+2, "ITS+TPC track --> TPC PID failed");
      //printf("[PHI] Track %4d: REJECTED\n", i);
      return kFALSE;
    }
    else if (MatchTOF(track) && !okTOF)
    {
      AliDebug(AliLog::kDebug+2, "ITS+TPC track --> TOF matched but TOF PID failed");
      //printf("[PHI] Track %4d: REJECTED\n", i);
      return kFALSE;
    }
    else if (track->IsOn(AliESDtrack::kITSpid) && track->P() <= fMaxITSmom && !okITS)
    {
      AliDebug(AliLog::kDebug+2, Form("ITS+TPC track --> Momentum lower than limit (%.2f) and ITS PID failed", fMaxITSmom));
      //printf("[PHI] Track %4d: REJECTED\n", i);
      return kFALSE;
    }
  }
  
  // this point is reached only if function didn't exit before due to somecheck not passed
  return kTRUE;
}
