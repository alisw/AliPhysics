//
// Implementation file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#include "Riostream.h"
#include <iomanip>

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
#include "AliITSPIDResponse.h"

#include "AliRsnAnalysisPhi7TeV.h"

//__________________________________________________________________________________________________
AliRsnAnalysisPhi7TeV::AliRsnAnalysisPhi7TeV(const char *name) :
  AliAnalysisTaskSE(name),
  fUseMC(kFALSE),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fPDG(0),
  fCh(0),
  fIM(0.0),
  fPt(0.0),
  fY(0.0),
  fEta(0.0),
  fMaxVz(1E6),
  fMaxITSband(1E6),
  fTPCpLimit(0.35),
  fMinTPCband(-1E6),
  fMaxTPCband( 1E6),
  fRsnTreeComp(0x0),
  fRsnTreeTrue(0x0),
  fOutList(0x0),
  fHEvents(0x0),
  fESDtrackCutsTPC(),
  fESDtrackCutsITS(),
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
  fCheckITS(copy.fCheckITS),
  fCheckTPC(copy.fCheckTPC),
  fCheckTOF(copy.fCheckTOF),
  fPDG(0),
  fCh(0),
  fIM(0.0),
  fPt(0.0),
  fY(0.0),
  fEta(0.0),
  fMaxVz(copy.fMaxVz),
  fMaxITSband(copy.fMaxITSband),
  fTPCpLimit(copy.fTPCpLimit),
  fMinTPCband(copy.fMinTPCband),
  fMaxTPCband(copy.fMaxTPCband),
  fRsnTreeComp(0x0),
  fRsnTreeTrue(0x0),
  fOutList(0x0),
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
  fCheckITS = copy.fCheckITS;
  fCheckTPC = copy.fCheckTPC;
  fCheckTOF = copy.fCheckTOF;

  fMaxVz   = copy.fMaxVz;
  fMaxITSband = copy.fMaxITSband;
  
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

  fRsnTreeComp->Branch("pdg", &fPDG, "pdg/S"   );
  fRsnTreeComp->Branch("ch" , &fCh , "ch/S"    );
  fRsnTreeComp->Branch("im" , &fIM , "im/F"    );
  fRsnTreeComp->Branch("y"  , &fY  , "y/F"     );
  fRsnTreeComp->Branch("pt" , &fPt , "pt/F"    );
  fRsnTreeComp->Branch("eta", &fEta, "eta/F"   );
  fRsnTreeComp->Branch("its", &fITS, "its[2]/S");

  OpenFile(2);
  fRsnTreeTrue = new TTree("rsnTrue", "True pairs");

  fRsnTreeTrue->Branch("im" , &fIM , "im/F" );
  fRsnTreeTrue->Branch("y"  , &fY  , "y/F"  );
  fRsnTreeTrue->Branch("pt" , &fPt , "pt/F" );
  fRsnTreeTrue->Branch("eta", &fEta, "eta/F");

  OpenFile(3);
  fOutList    = new TList;
  fHEvents    = new TH1I("hEvents", "Event details", 5, 0, 5);
  fVertexX[0] = new TH1F("hVertexTracksX", "X position of primary vertex (tracks)", 200,  -2,  2);
  fVertexY[0] = new TH1F("hVertexTracksY", "Y position of primary vertex (tracks)", 200,  -2,  2);
  fVertexZ[0] = new TH1F("hVertexTracksZ", "Z position of primary vertex (tracks)", 400, -40, 40);
  fVertexX[1] = new TH1F("hVertexSPDX", "X position of primary vertex (SPD)", 200,  -2,  2);
  fVertexY[1] = new TH1F("hVertexSPDY", "Y position of primary vertex (SPD)", 200,  -2,  2);
  fVertexZ[1] = new TH1F("hVertexSPDZ", "Z position of primary vertex (SPD)", 400, -40, 40);
  
  fHEvents->GetXaxis()->SetBinLabel(1, "Good vertex with tracks");
  fHEvents->GetXaxis()->SetBinLabel(2, "Good vertex with SPD");
  fHEvents->GetXaxis()->SetBinLabel(3, "Far vertex with tracks");
  fHEvents->GetXaxis()->SetBinLabel(4, "Far vertex with SPD");
  fHEvents->GetXaxis()->SetBinLabel(5, "No good vertex");

  fOutList->Add(fHEvents);
  fOutList->Add(fVertexX[0]);
  fOutList->Add(fVertexY[0]);
  fOutList->Add(fVertexZ[0]);
  fOutList->Add(fVertexX[1]);
  fOutList->Add(fVertexY[1]);
  fOutList->Add(fVertexZ[1]);
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
  
  // check the event
  Int_t eval = EventEval(esd);
  fHEvents->Fill(eval);
  
  // if the event is good for analysis, process it
  if (eval == kGoodTracksPrimaryVertex || eval == kGoodSPDPrimaryVertex)
  {
    ProcessESD(esd, stack);
    ProcessMC(stack);
  }
  
  // update histogram container
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
Int_t AliRsnAnalysisPhi7TeV::EventEval(AliESDEvent *esd)
{
//
// Checks if the event is good for analysis.
// Returns one of the flag values defined in the header
//

  static Int_t evNum = 0;
  evNum++;
  
  // get the best primary vertex:
  // first try the one with tracks
  const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
  const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
  Double_t            vzTrk = 1000.0;
  Double_t            vzSPD = 1000.0;
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

  // ITS stuff #1 create the response function
  Bool_t isMC = (stack != 0x0);
  AliITSPIDResponse itsrsp(isMC);

  // TOF stuff #1: init OCDB
  Int_t run = esd->GetRunNumber();
  if (run != lastRun)
  {
    cout << "Run = " << run << " -- LAST = " << lastRun << endl;
    lastRun = run;
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("raw://");
    cdb->SetRun(run);
    fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
    fTOFcalib->Init();
  }
  //AliCDBManager *cdb = AliCDBManager::Instance();
  //cdb->SetDefaultStorage("raw://");
  //cdb->SetRun(run);
  // TOF stuff #2: init calibration
  //fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
  //fTOFcalib->Init();
  // TOF stuff #3: calibrate
  if (fTOFcalibrateESD) fTOFcalib->CalibrateESD(esd);
  if (fTOFtuneMC) fTOFmaker->TuneForMC(esd);
  if (fTOFuseT0) 
  {
    fTOFmaker->ComputeT0TOF(esd);
    fTOFmaker->ApplyT0TOF(esd);
    fESDpid->MakePID(esd, kFALSE, 0.);
  }
  // TOF stuff #4: define fixed function for compatibility range
  Double_t a1 = 0.01, a2 = -0.03;
  Double_t b1 = 0.25, b2 =  0.25;
  Double_t c1 = 0.05, c2 = -0.03;
  Double_t ymax, ymin;

  // prepare to look on all tracks to select the ones
  // which pass all the cuts
  Int_t   ntracks = esd->GetNumberOfTracks();
  TArrayI pos(ntracks);
  TArrayI neg(ntracks);
  TArrayI itspos(ntracks);
  TArrayI itsneg(ntracks);
  
  // loop on all tracks
  ULong_t  status;
  Int_t    i, k, charge, npos = 0, nneg = 0, nITS;
  Double_t times[10], tpcNSigma, tpcMaxNSigma, itsSignal, itsNSigma, mom, tofTime, tofSigma, tofRef, tofDiff, tofRel;
  Bool_t   okQuality, okTOF, okTPC, okITS, okTrack, isTPC, isITSSA;
  UChar_t  itsCluMap;
  for (i = 0; i < ntracks; i++)
  {
    AliESDtrack *track = esd->GetTrack(i);
    if (!track) continue;
    
    // get commonly used variables
    status  = (ULong_t)track->GetStatus();
    mom     = track->P();
    isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
    isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
    //cout << isTPC << ' ' << isITSSA << endl;
    
    // accept only tracks which are TPC+ITS or ITS standalone
    if (!isTPC && !isITSSA) continue;
    
    // define selection properties depending on track type
    // it track is standard TPC+ITS+TOF track, check standard cuts and TOF
    // if the track is an ITS standalone, check its specific cuts only
    okTrack = kTRUE;
    
    if (isTPC)
    {
      // check standard ESD cuts
      okQuality = fESDtrackCutsTPC.IsSelected(track);
      //cout << "GLOBAL -- quality = " << (okQuality ? "GOOD" : "BAD") << endl;
      if (!okQuality) continue;
      
      // check TPC dE/dx
      if (fCheckTPC)
      {
        tpcNSigma = TMath::Abs(fESDpid->NumberOfSigmasTPC(track, AliPID::kKaon));
        if (track->GetInnerParam()->P() > fTPCpLimit) tpcMaxNSigma = fMinTPCband; else tpcMaxNSigma = fMaxTPCband;
        okTPC = (tpcNSigma <= tpcMaxNSigma);
        //cout << "ALT -- TPC    -- nsigma = " << tpcNSigma << ", max = " << tpcMaxNSigma << " --> " << (okTPC ? "OK" : "FAILED") << endl;
        //cout << "ALTTPC -- " << fTPCpar[0] << ' ' << fTPCpar[1] << ' ' << fTPCpar[2] << ' ' << fTPCpar[3] << ' ' << fTPCpar[4] << endl;
      }
      else
      {
        okTPC = kTRUE;
      }
      
      // check TOF (only if momentum is large than function asymptote and flags are OK)
      if (fCheckTOF)
      {
        if (((status & AliESDtrack::kTOFout) != 0) && ((status & AliESDtrack::kTIME) != 0) && mom > TMath::Max(b1, b2))
        {
          track->GetIntegratedTimes(times);
          tofTime  = (Double_t)track->GetTOFsignal();
          tofSigma = fTOFmaker->GetExpectedSigma(mom, times[AliPID::kKaon], AliPID::ParticleMass(AliPID::kKaon));
          tofRef   = times[AliPID::kKaon];
          if (tofRef > 0.0)
          {
            tofDiff  = (tofTime - tofRef);
            tofRel   = (tofTime - tofRef) / tofRef;
            ymax     = a1 / (mom - b1) + c1;
            ymin     = a2 / (mom - b2) + c2;
            okTOF    = (tofRel >= ymin && tofRel <= ymax);
            //cout << "TOF    -- diff = " << tofDiff << ", rel diff = " << tofRel << ", range = " << ymin << " to " << ymax << ", sigma = " << tofSigma << " --> " << (okTOF ? "OK" : "FAILED") << endl;
          }
          else
          {
            okTOF = kTRUE;
            //cout << "TOF    -- not checked due to ZERO reference time" << endl;
          }
        }
        else
        {
          okTOF = kTRUE;
          //cout << "TOF    -- not checked because TOF pid absent" << endl;
        }
      }
      else
      {
        okTOF = kTRUE;
      }
      
      okTrack = okQuality && okTPC && okTOF;
      //cout << "GLOBAL -- overall = " << (okTrack ? "ACCEPTED" : "REJECTED") << endl;
    }
    else
    {
      // check standard ESD cuts
      okQuality = fESDtrackCutsITS.IsSelected(track);
      //cout << "ITSSA  -- quality = " << (okQuality ? "GOOD" : "BAD") << endl;
      if (!okQuality) continue;
      
      // check dE/dx
      if (fCheckITS)
      {
        itsSignal = track->GetITSsignal();
        itsCluMap = track->GetITSClusterMap();
        nITS      = 0;
        for(k = 2; k < 6; k++) if(itsCluMap & (1 << k)) ++nITS;
        if (nITS < 3) 
        {
          okITS = kFALSE;
          //cout << "ITS    -- not checked due to too few PID clusters" << endl;
        }
        else
        {
          itsNSigma = itsrsp.GetNumberOfSigmas(mom, itsSignal, AliPID::kKaon, nITS, kTRUE);
          okITS = (TMath::Abs(itsNSigma) <= fMaxITSband);
          //cout << "ITS    -- nsigma = " << itsNSigma << ", max = " << fMaxITSband << " --> " << (okITS ? "OK" : "FAILED") << endl;
        }
      }
      else
      {
        okITS = kTRUE;
      }
      
      okTrack = okQuality && okITS;
      //cout << "ITSSA  -- overall = " << (okTrack ? "ACCEPTED" : "REJECTED") << endl;
    }
    
    // skip tracks not passing cuts
    if (!okTrack) continue;
    
    // if all checks are passed, add the track index in one of the
    // charged tracks arrays
    charge = (Int_t)track->Charge();
    if (charge > 0)
    {
      pos[npos] = i;
      if (isITSSA) itspos[npos] = 1; else itspos[npos] = 0;
      npos++;
    }
    else if (charge < 0)
    {
      neg[nneg] = i;
      if (isITSSA) itsneg[nneg] = 1; else itsneg[nneg] = 0;
      nneg++;
    }
  }
  
  // resize arrays accordingly
  pos.Set(npos);
  neg.Set(nneg);
  itspos.Set(npos);
  itsneg.Set(nneg);
  
  /*
  // fill unlike-sign
  Int_t i1, i2;
  for (i1 = 0; i1 < npos; i1++)
    for (i2 = 0; i2 < nneg; i2++)
      AddEntryFromESD(esd, i1, i2, itspos[i1], itsneg[i2], 0, stack);
      
  // fill ++ like-sign
  for (i1 = 0; i1 < npos; i1++)
    for (i2 = i1 + 1; i2 < npos; i2++)
      AddEntryFromESD(esd, i1, i2, itspos[i1], itspos[i2], 1, 0x0);
      
  // fill -- like-sign
  for (i1 = 0; i1 < nneg; i1++)
    for (i2 = i1 + 1; i2 < nneg; i2++)
      AddEntryFromESD(esd, i1, i2, itsneg[i1], itsneg[i2], -1, 0x0);
      */
      
  // loop on unlike-sign pairs to compute invariant mass signal
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
      if (pos[ip] == neg[in]) 
      {
        AliError("POS = NEG");
        continue;
      }
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

      fCh     = 0;
      fIM     = (Float_t)vsum.M();
      fPt     = (Float_t)vsum.Perp();
      fEta    = (Float_t)vsum.Eta();
      fY      = (Float_t)vref.Rapidity();
      fITS[0] = itspos[ip];
      fITS[1] = itsneg[in];

      //if (fIM < 0.9 || fIM >  1.4) continue;
      //if (fPt < 0.0 || fPt > 20.0) continue;
      
      fRsnTreeComp->Fill();
    }
  }
  
  // loop on like-sign pairs to compute invariant mass background
  Int_t           i1, i2;
  AliESDtrack    *t1 = 0x0, *t2 = 0x0;
  TLorentzVector  v1, v2;
  
  // pos-pos
  for (i1 = 0; i1 < npos; i1++)
  {
    t1 = esd->GetTrack(pos[i1]);

    for (i2 = i1+1; i2 < npos; i2++)
    {
      t2 = esd->GetTrack(pos[i2]);

      v1.SetXYZM(t1->Px(), t1->Py(), t1->Pz(), kmass);
      v2.SetXYZM(t2->Px(), t2->Py(), t2->Pz(), kmass);
      vsum = v1 + v2;
      vref.SetXYZM(vsum.X(), vsum.Y(), vsum.Z(), phimass);

      fPDG    = 0;
      fCh     = 1;
      fIM     = (Float_t)vsum.M();
      fPt     = (Float_t)vsum.Perp();
      fEta    = (Float_t)vsum.Eta();
      fY      = (Float_t)vref.Rapidity();
      fITS[0] = itspos[i1];
      fITS[1] = itspos[i2];

      //if (fIM < 0.9 || fIM >  1.4) continue;
      //if (fPt < 0.0 || fPt > 20.0) continue;
      
      fRsnTreeComp->Fill();
    }
  }
  // neg-neg
  for (i1 = 0; i1 < nneg; i1++)
  {
    t1 = esd->GetTrack(neg[i1]);

    for (i2 = i1+1; i2 < nneg; i2++)
    {
      t2 = esd->GetTrack(neg[i2]);

      v1.SetXYZM(t1->Px(), t1->Py(), t1->Pz(), kmass);
      v2.SetXYZM(t2->Px(), t2->Py(), t2->Pz(), kmass);
      vsum = v1 + v2;
      vref.SetXYZM(vsum.X(), vsum.Y(), vsum.Z(), phimass);

      fPDG    = 0;
      fCh     = -1;
      fIM     = (Float_t)vsum.M();
      fPt     = (Float_t)vsum.Perp();
      fEta    = (Float_t)vsum.Eta();
      fY      = (Float_t)vref.Rapidity();
      fITS[0] = itsneg[i1];
      fITS[1] = itsneg[i2];

      //if (fIM < 0.9 || fIM >  1.4) continue;
      //if (fPt < 0.0 || fPt > 20.0) continue;
      
      fRsnTreeComp->Fill();
    }
  }

  PostData(1, fRsnTreeComp);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisPhi7TeV::AddEntryFromESD(AliESDEvent *esd, Int_t i1, Int_t i2, Int_t its1, Int_t its2, Short_t charge, AliStack *stack)
{
//
// Add an entry to the output TTree computed from two tracks in the ESD
//

  // loop on unlike-sign pairs to compute invariant mass signal
  AliPID          pid;
  Double_t        kmass = pid.ParticleMass(AliPID::kKaon);
  Double_t        phimass = 1.019455;
  TLorentzVector  v1, v2, vsum, vref;
  
  AliESDtrack *t1 = esd->GetTrack(i1);
  AliESDtrack *t2 = esd->GetTrack(i2);

  v1.SetXYZM(t1->Px(), t1->Py(), t1->Pz(), kmass);
  v2.SetXYZM(t2->Px(), t2->Py(), t2->Pz(), kmass);
  vsum = v1 + v2;
  vref.SetXYZM(vsum.X(), vsum.Y(), vsum.Z(), phimass);
  
  // if stack is present, search for true pair
  fPDG = 0;
  Int_t l1 = TMath::Abs(t1->GetLabel());
  Int_t l2 = TMath::Abs(t2->GetLabel());
  if (stack) 
  {
    TParticle *part1 = stack->Particle(l1);
    TParticle *part2 = stack->Particle(l2);
    if (part1 && part2)
    {
      if (part1->GetFirstMother() == part2->GetFirstMother())
      {
        if (part1->GetFirstMother() > 0)
        {
          TParticle *mum = stack->Particle(part1->GetFirstMother());
          fPDG = mum->GetPdgCode();
        }
      }
    }
  }
  
  fPDG = TMath::Abs(fPDG);
  fCh     = charge;
  fIM     = (Float_t)vsum.M();
  fPt     = (Float_t)vsum.Perp();
  fEta    = (Float_t)vsum.Eta();
  fY      = (Float_t)vref.Rapidity();
  fITS[0] = its1;
  fITS[1] = its2;

  //if (fIM < 0.9 || fIM >  1.4) return;
  //if (fPt < 0.0 || fPt > 20.0) return;
      
  fRsnTreeComp->Fill();
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
