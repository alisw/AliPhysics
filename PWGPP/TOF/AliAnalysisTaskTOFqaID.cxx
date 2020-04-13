/*  created by fbellini@cern.ch on 29/04/2013 */
/*  last modified by fbellini   on 19/08/2013 */

#ifndef ALIANALYSISTASKTOFQAID_CXX
#define ALIANALYSISTASKTOFQAID_CXX

#include "AliAnalysisFilter.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliTOFChannelOnlineStatusArray.h"
#include "AliTOFGeometry.h"
#include "AliTOFPIDParams.h"
#include "AliTOFRawStream.h"
#include "AliTOFT0maker.h"
#include "AliTOFT0v1.h"
#include "AliTOFcalib.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THashList.h"
#include "TProfile.h"
#include "TTree.h"
//
#include "AliAnalysisTaskTOFqaID.h"

ClassImp(AliAnalysisTaskTOFqaID);

//________________________________________________________________________
AliAnalysisTaskTOFqaID::AliAnalysisTaskTOFqaID()
    : AliAnalysisTaskTOFqaID("AliAnalysisTaskTOFqaID")
{
}

//________________________________________________________________________
AliAnalysisTaskTOFqaID::AliAnalysisTaskTOFqaID(const char* name)
    : AliAnalysisTaskSE(name)
    , fVariableBinsPt(300 + 1)
    , fVariableBinsMult(100 + 1)
    , fRunNumber(0)
    , fESD(0x0)
    , fMCevent(0x0)
    , fTrackFilter(0x0)
    , fVertex(0x0)
    , fESDpid(new AliESDpid())
    , fTOFHeader(0x0)
    , fEnableAdvancedCheck(kFALSE)
    , fEnableChargeSplit(kFALSE)
    , fExpTimeBinWidth(24.4)
    , fExpTimeRangeMin(-25010.)
    , fExpTimeRangeMax(25010.)
    , fExpTimeSmallRangeMin(-5002.)
    , fExpTimeSmallRangeMax(5002.)
    , fnExpTimeBins(1)
    , fnExpTimeSmallBins(1)
    , fMyTimeZeroTOF(3)
    , fMyTimeZeroTOFsigma(3)
    , fMyTimeZeroTOFtracks(3)
    , fMyTimeZeroTOFstatus(kFALSE)
    , fIsMC(kFALSE)
    , fVerbose(kFALSE)
    , fUseTOFT0CalibMode(kFALSE)
    , fSelectedPdg(0)
    , fP(1e10)
    , fPt(1e10)
    , fEta(1e10)
    , fPhi(1e10)
    , fTPCOuterPhi(1e10)
    , fL(1e10)
    , fMatchingMomCut(1.0)
    , fMatchingEtaCut(0.8)
    , fTof(1e10)
    , fMCTOFMatch(-1.)
    , fOCDBLocation("local://$ALICE_PHYSICS/OCDB")
    , fChannelArray(0x0)
    , fCalib(0x0)
    , fResolutionMinP(.95)
    , fResolutionMaxP(1.05)
    , fHlist(0x0)
    , fHlistTimeZero(0x0)
    , fHlistPID(0x0)
    , fHlistTRD(0x0)
    , fHlistTrigger(0x0)
{
  // Constructor
  // Define input and output slots here
  Info("AliAnalysisTaskTOFqaID", "Calling Constructor");

  for (Int_t j = 0; j < 5; j++) {
    if (j < 3) {
      fT0[j] = 0.0;
      fNTOFtracks[j] = 0;
    }
    fSigmaSpecie[j] = 0.0;
    fTrkExpTimes[j] = 0.0;
    fThExpTimes[j] = 0.0;
  }
  //
  fMyTimeZeroTOF.Reset(1e20);
  fMyTimeZeroTOFsigma.Reset(1e20);
  fMyTimeZeroTOFtracks.Reset(-1);
  //
  fVariableBinsPt.Reset(0.);
  fVariableBinsMult.Reset(0.);
  //
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());

  // Output slot #0 writes into a TH1 container
  // Output slot #1 writes into a user defined  container
  DefineOutput(1, THashList::Class());
  DefineOutput(2, THashList::Class());
  DefineOutput(3, THashList::Class());
  DefineOutput(4, THashList::Class());
  DefineOutput(5, THashList::Class());
}

//________________________________________________________________________
#define cpVar(var) var(copy.var)
AliAnalysisTaskTOFqaID::AliAnalysisTaskTOFqaID(const AliAnalysisTaskTOFqaID& copy)
    : AliAnalysisTaskSE()
    , cpVar(fVariableBinsPt)
    , cpVar(fVariableBinsMult)
    , cpVar(fRunNumber)
    , cpVar(fESD)
    , cpVar(fMCevent)
    , cpVar(fTrackFilter)
    , cpVar(fVertex)
    , cpVar(fESDpid)
    , cpVar(fTOFHeader)
    , cpVar(fEnableAdvancedCheck)
    , cpVar(fEnableChargeSplit)
    , cpVar(fExpTimeBinWidth)
    , cpVar(fExpTimeRangeMin)
    , cpVar(fExpTimeRangeMax)
    , cpVar(fExpTimeSmallRangeMin)
    , cpVar(fExpTimeSmallRangeMax)
    , cpVar(fnExpTimeBins)
    , cpVar(fnExpTimeSmallBins)
    , cpVar(fMyTimeZeroTOF)
    , cpVar(fMyTimeZeroTOFsigma)
    , cpVar(fMyTimeZeroTOFtracks)
    , cpVar(fMyTimeZeroTOFstatus)
    , cpVar(fIsMC)
    , cpVar(fVerbose)
    , cpVar(fUseTOFT0CalibMode)
    , cpVar(fSelectedPdg)
    , cpVar(fP)
    , cpVar(fPt)
    , cpVar(fEta)
    , cpVar(fPhi)
    , cpVar(fTPCOuterPhi)
    , cpVar(fL)
    , cpVar(fMatchingMomCut)
    , cpVar(fMatchingEtaCut)
    , cpVar(fTof)
    , cpVar(fMCTOFMatch)
    , cpVar(fOCDBLocation)
    , cpVar(fChannelArray)
    , cpVar(fCalib)
    , cpVar(fResolutionMinP)
    , cpVar(fResolutionMaxP)
    , cpVar(fHlist)
    , cpVar(fHlistTimeZero)
    , cpVar(fHlistPID)
    , cpVar(fHlistTRD)
    , cpVar(fHlistTrigger)
{
  // Copy constructor
  for (Int_t j = 0; j < 5; j++) {
    if (j < 3) {
      fT0[j] = copy.fT0[j];
      fNTOFtracks[j] = copy.fNTOFtracks[j];
    }
    fSigmaSpecie[j] = copy.fSigmaSpecie[j];
    fTrkExpTimes[j] = copy.fTrkExpTimes[j];
    fThExpTimes[j] = copy.fThExpTimes[j];
  }
  //
}
#undef cpVar

//___________________________________________________________________________
#define cpVar(var) var = copy.var;
AliAnalysisTaskTOFqaID& AliAnalysisTaskTOFqaID::operator=(const AliAnalysisTaskTOFqaID& copy)
{
  //
  // Assignment operator
  //
  if (this != &copy) {
    AliAnalysisTaskSE::operator=(copy);
    cpVar(fRunNumber);
    cpVar(fVariableBinsPt);
    cpVar(fVariableBinsMult);
    cpVar(fESD);
    cpVar(fMCevent);
    cpVar(fTrackFilter);
    cpVar(fVertex);
    cpVar(fESDpid);
    cpVar(fTOFHeader);
    cpVar(fEnableAdvancedCheck);
    cpVar(fEnableChargeSplit);
    cpVar(fExpTimeBinWidth);
    cpVar(fExpTimeRangeMin);
    cpVar(fExpTimeRangeMax);
    cpVar(fExpTimeSmallRangeMin);
    cpVar(fExpTimeSmallRangeMax);
    cpVar(fnExpTimeBins);
    cpVar(fnExpTimeSmallBins);
    cpVar(fMyTimeZeroTOF);
    cpVar(fMyTimeZeroTOFsigma);
    cpVar(fMyTimeZeroTOFtracks);
    cpVar(fMyTimeZeroTOFstatus);
    cpVar(fIsMC);
    cpVar(fVerbose);
    cpVar(fUseTOFT0CalibMode);
    cpVar(fSelectedPdg);
    cpVar(fP);
    cpVar(fPt);
    cpVar(fEta);
    cpVar(fPhi);
    cpVar(fTPCOuterPhi);
    cpVar(fL);
    cpVar(fMatchingMomCut);
    cpVar(fMatchingEtaCut);
    cpVar(fTof);
    cpVar(fMCTOFMatch);
    cpVar(fOCDBLocation);
    cpVar(fChannelArray);
    cpVar(fCalib);
    cpVar(fResolutionMinP);
    cpVar(fResolutionMaxP);
    cpVar(fHlist);
    cpVar(fHlistTimeZero);
    cpVar(fHlistPID);
    cpVar(fHlistTRD);
    cpVar(fHlistTrigger);
    for (Int_t j = 0; j < 5; j++) {
      if (j < 3) {
        fT0[j] = copy.fT0[j];
        fNTOFtracks[j] = copy.fNTOFtracks[j];
      }
      fSigmaSpecie[j] = copy.fSigmaSpecie[j];
      fTrkExpTimes[j] = copy.fTrkExpTimes[j];
      fThExpTimes[j] = copy.fThExpTimes[j];
    }
    //
  }
  return *this;
}
#undef cpVar
//___________________________________________________________________________
AliAnalysisTaskTOFqaID::~AliAnalysisTaskTOFqaID()
{
  //
  //destructor
  //

  Info("~AliAnalysisTaskTOFqaID", "Calling Destructor");
  if (fESDpid)
    delete fESDpid;
  if (fTOFHeader)
    delete fTOFHeader;
  if (fVertex)
    delete fVertex;
  if (fTrackFilter)
    delete fTrackFilter;
  if (fChannelArray)
    delete fChannelArray;
  if (fCalib)
    delete fCalib;
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    return;

  if (fHlist) {
    delete fHlist;
    fHlist = 0;
  }
  if (fHlistTimeZero) {
    delete fHlistTimeZero;
    fHlistTimeZero = 0;
  }
  if (fHlistPID) {
    delete fHlistPID;
    fHlistPID = 0;
  }
  if (fHlistTRD) {
    delete fHlistTRD;
    fHlistTRD = 0;
  }
  if (fHlistTrigger) {
    delete fHlistTrigger;
    fHlistTrigger = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskTOFqaID::UserCreateOutputObjects()
{
  //
  //Define output objects and histograms
  //
  //retrieve PID response object
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  if (!man)
    AliFatal("Analysis manager needed");
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler)
    AliFatal("Input handler needed");
  //pid response object
  fESDpid = (AliESDpid*)inputHandler->GetPIDResponse();
  if (!fESDpid)
    AliError("PIDResponse object was not created");
  //fESDpid->SetOADBPath("$ALICE_PHYSICS/OADB");

  Info("CreateOutputObjects", "CreateOutputObjects (THashList) of task %s", GetName());
  OpenFile(1);

  //define variable binning for pt and p before creating histos
  SetVariableBinning();

  if (!fHlist)
    fHlist = new THashList();
  fHlist->SetOwner(kTRUE);
  fHlist->SetName("base");

  if (!fHlistTimeZero)
    fHlistTimeZero = new THashList();
  fHlistTimeZero->SetOwner(kTRUE);
  fHlistTimeZero->SetName("startTime");

  if (!fHlistPID)
    fHlistPID = new THashList();
  fHlistPID->SetOwner(kTRUE);
  fHlistPID->SetName("pid");

  if (!fHlistTRD)
    fHlistTRD = new THashList();
  fHlistTRD->SetOwner(kTRUE);
  fHlistTRD->SetName("TRD");

  if (!fHlistTrigger)
    fHlistTrigger = new THashList();
  fHlistTrigger->SetOwner(kTRUE);
  fHlistTrigger->SetName("trigger");

  if (fExpTimeRangeMax < fExpTimeRangeMin) {
    SetExpTimeHistoRange(-25010., 25010.);
  }
  fnExpTimeBins = TMath::Nint((fExpTimeRangeMax - fExpTimeRangeMin) / fExpTimeBinWidth); //ps
  fExpTimeRangeMax = fExpTimeRangeMin + fnExpTimeBins * fExpTimeBinWidth;                //ps

  if (fExpTimeSmallRangeMax < fExpTimeSmallRangeMin) {
    SetExpTimeHistoSmallRange(-5002., 5002.);
  }
  fnExpTimeSmallBins = TMath::Nint((fExpTimeSmallRangeMax - fExpTimeSmallRangeMin) / fExpTimeBinWidth); //ps
  fExpTimeSmallRangeMax = fExpTimeSmallRangeMin + fnExpTimeSmallBins * fExpTimeBinWidth;                //ps

  //add plots for start time QA
  AddStartTimeHisto(fHlistTimeZero, "");

  //add plots for base TOF quantities
  if (fEnableChargeSplit) {
    AddTofBaseHisto(fHlist, 1, "");
    AddTofBaseHisto(fHlist, -1, "");
  } else {
    AddTofBaseHisto(fHlist, 0, "");
  }
  //add plots for matching efficiency
  if (fEnableChargeSplit) {
    AddMatchingEffHisto(fHlist, 1, "");
    AddMatchingEffHisto(fHlist, -1, "");
    if (fIsMC) {
      AddMatchingEffHisto(fHlist, 1, "_TrueMatch");
      AddMatchingEffHisto(fHlist, -1, "_TrueMatch");
    }
  } else {
    AddMatchingEffHisto(fHlist, 0, "");
    if (fIsMC)
      AddMatchingEffHisto(fHlist, 0, "_TrueMatch");
  }
  //add plots for PID checks
  if (fEnableChargeSplit) {
    AddPidHisto(fHlistPID, 1, "");
    AddPidHisto(fHlistPID, -1, "");
  } else {
    AddPidHisto(fHlistPID, 0, "");
  }
  //add trd plots
  if (fEnableAdvancedCheck) {
    AddTrdHisto();
  }
  //Add trigger plots
  AddTofTrgHisto("");

  PostData(1, fHlist);
  PostData(2, fHlistTimeZero);
  PostData(3, fHlistPID);
  PostData(4, fHlistTRD);
  PostData(5, fHlistTrigger);
}
//________________________________________________________________________
void AliAnalysisTaskTOFqaID::UserExec(Option_t*)
{
  /* Main - executed for each event.
    It extracts event information and track information after selecting
    primary tracks via standard cuts. */
  /*
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    AliError("ERROR: Could not get ESDInputHandler");
    return;
  } else {
    fESD = (AliESDEvent*) esdH->GetEvent();
  }

  */
  fESD = (AliESDEvent*)InputEvent();
  if (!fESD) {
    AliError("fESD event not available");
    return;
  }

  if (!fESDpid) {
    AliError("PID object fESDpid not available");
    return;
  }

  //retrieve default start time type from PIDresponse
  AliPIDResponse::EStartTimeType_t startTimeMethodDefault = AliPIDResponse::kBest_T0;
  if (fESDpid->GetTOFPIDParams()) { // during reconstruction OADB not yet available
    startTimeMethodDefault = ((AliTOFPIDParams*)fESDpid->GetTOFPIDParams())->GetStartTimeMethod();
  }

  //access MC event handler for MC truth information
  if (fIsMC) {
    AliMCEventHandler* mcH = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcH) {
      AliError("Cannot get MCeventHandler");
      return;
    } else {
      fMCevent = (AliMCEvent*)mcH->MCEvent();
      if (!fMCevent) {
        AliError("Trying to retrieve an invalid MC event.");
        return;
      }
      fESDpid->SetCurrentMCEvent(fMCevent);
    }
  }

  // get run number
  Int_t runNb = fESD->GetRunNumber();
  if (runNb != fRunNumber) {
    //retrieve maps of channels from OCDB
    fRunNumber = runNb;
    LoadChannelMapsFromOCDB();
  }

  //reset matched track counters
  for (Int_t j = 0; j < 3; j++) {
    fNTOFtracks[j] = 0;
  }

  //Get vertex info and apply vertex cut
  if (!IsEventSelected(fESD))
    return;

  //set response tof_t0 for all other checks
  fESDpid->SetTOFResponse(fESD, AliESDpid::kTOF_T0); //(fill_t0, tof_t0, t0_t0, best_t0)

  AliDebug(3, Form("Momentum cut for eta and phi distributions set: Pt>%3.2f", fMatchingMomCut));

  //check existence of track filter
  if (!fTrackFilter) {
    AliInfo("No track filter found, skipping the track loop");
    return;
  }

  //Computing our own TOF_T0 once per event
  ComputeTimeZeroByTOF1GeV();

  // loop over ESD tracks
  AliESDtrack* track = 0x0;
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    track = fESD->GetTrack(iTracks);
    if (!track) {
      AliInfo(Form("Cannot receive track %d", iTracks));
      continue;
    }

    //primary tracks selection: kTPCrefit and std cuts
    if (!fTrackFilter->IsSelected(track))
      continue;

    //select specie if MC
    if (fIsMC && (!SelectMCspecies(fMCevent, track))) {
      AliDebug(4, Form("MC tracks selection: Track=%i  label=%i  Not Accepted", iTracks, track->GetLabel()));
      continue;
    }

    //apply cut for eta acceptance
    fEta = track->Eta();
    if (TMath::Abs(fEta) > fMatchingEtaCut)
      continue;

    //get other track variables
    fP = track->P();
    fPt = track->Pt();
    fPhi = track->Phi() * TMath::RadToDeg();
    fTPCOuterPhi = GetPhiAtTPCouterRadius(track);
    fL = track->GetIntegratedLength();
    track->GetIntegratedTimes(fTrkExpTimes);

    Int_t charge = 0;
    if (fEnableChargeSplit)
      charge = track->Charge();

    //Fill histograms for primary particles
    FillPrimaryTrkHisto(charge, "");

    //Fill histrograms for tracks matched at TOF (no check on true matching)
    if (IsTPCTOFMatched(track, 0)) {
      fTof = track->GetTOFsignal() * 1E-3; //in ps
      //increment track counters
      fNTOFtracks[0]++;
      if (charge > 0)
        fNTOFtracks[1]++;
      if (charge < 0)
        fNTOFtracks[2]++;
      //fill histos
      FillTofBaseHisto(track, charge, "");
      FillMatchedTrkHisto(charge, "");
      FillPidHisto(track, charge, "");
    }

    //Fill histograms for tracks truly matched at TOF (based on the matching label)
    if (fIsMC) {
      FillPrimaryTrkHisto(charge, "_TrueMatch");
      if (IsTPCTOFMatched(track, 1))
        FillMatchedTrkHisto(charge, "_TrueMatch");
    }

    if (fEnableAdvancedCheck)
      FillTrdHisto(track, charge);
  } //end loop on tracks

  //fill time zero histos
  FillStartTimeHisto("");
  if (fEnableChargeSplit) {
    ((TH1F*)fHlist->FindObject("hTOFmulti_pos"))->Fill(fNTOFtracks[1]);
    ((TH1F*)fHlist->FindObject("hTOFmulti_neg"))->Fill(fNTOFtracks[2]);
  } else {
    ((TH1F*)fHlist->FindObject("hTOFmulti_all"))->Fill(fNTOFtracks[0]);
  }
  //fill TOF trg histos from infos in TOF header
  fTOFHeader = (AliTOFHeader*)fESD->GetTOFHeader();
  if (!fTOFHeader) {
    AliWarning("Cannot get TOF header: no TOF trigger info available");
  } else {
    FillTofTrgHisto("");
  }

  //restore value set by AliPIDResponseTask for subsequent wagons
  fESDpid->SetTOFResponse(fESD, startTimeMethodDefault);

  PostData(1, fHlist);
  PostData(2, fHlistTimeZero);
  PostData(3, fHlistPID);
  PostData(4, fHlistTRD);
  PostData(5, fHlistTrigger);
}

//________________________________________________________________________
void AliAnalysisTaskTOFqaID::Terminate(Option_t*)
{
  //check on output validity
  fHlist = dynamic_cast<THashList*>(GetOutputData(1));
  if (!fHlist) {
    AliError("Base histograms list not available");
    return;
  }

  // TH1F*hDummy = ((TH1F*)fHlist->FindObject("hTOFmatchedESDPt"));
  // TH1F*hMatchingEff = (TH1F*) hDummy->Clone("hMatchingEff");
  // hMatchingEff->SetTitle("Matching efficiency");
  // hMatchingEff->Divide((TH1F*) fHlist->FindObject("hESDprimaryTrackPt"));
  // TCanvas *c1 = new TCanvas("AliAnalysisTaskTOFqaID","Matching vs Pt",10,10,510,510);
  // c1->cd(1)->SetLogy();
  // hMatchingEff->DrawCopy("E");
  // fHlist->AddLast(hMatchingEff);
  ComputeMatchingEfficiency(fHlist, "pt");
  ComputeMatchingEfficiency(fHlist, "eta");
  ComputeMatchingEfficiency(fHlist, "phi");

  PostData(1, fHlist);
}

//---------------------------------------------------------------
Int_t AliAnalysisTaskTOFqaID::GetStripIndex(const Int_t* in)
{
  /* return tof strip index between 0 and 91 */

  Int_t nStripA = AliTOFGeometry::NStripA();
  Int_t nStripB = AliTOFGeometry::NStripB();
  Int_t nStripC = AliTOFGeometry::NStripC();

  Int_t iplate = in[1];
  Int_t istrip = in[2];

  Int_t stripOffset = 0;
  switch (iplate) {
  case 0:
    stripOffset = 0;
    break;
  case 1:
    stripOffset = nStripC;
    break;
  case 2:
    stripOffset = nStripC + nStripB;
    break;
  case 3:
    stripOffset = nStripC + nStripB + nStripA;
    break;
  case 4:
    stripOffset = nStripC + nStripB + nStripA + nStripB;
    break;
  default:
    stripOffset = -1;
    break;
  };

  if (stripOffset < 0 || stripOffset > 92)
    return -1;
  else
    return (stripOffset + istrip);
}

//-----------------------------------------------------------------
Double_t AliAnalysisTaskTOFqaID::GetPhiAtTPCouterRadius(AliESDtrack* track)
{
  //get track phi at TPC outer radius
  if (!track)
    return 1e10;
  Double_t tpcoutcoord[3] = { 0., 0., 0. };
  track->GetOuterXYZ(tpcoutcoord);
  Double_t phiOuterTPC = TMath::ATan2(tpcoutcoord[1], tpcoutcoord[0]) * TMath::RadToDeg();
  if (phiOuterTPC < 0)
    phiOuterTPC += (2 * TMath::Pi() * TMath::RadToDeg());
  return phiOuterTPC;
}
//-----------------------------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::IsEventSelected(AliESDEvent* event)
{
  //select event based on primary vertex
  if (!event) {
    AliError("Invalid ESD event");
    return kFALSE;
  }
  fVertex = (AliESDVertex*)event->GetPrimaryVertexTracks();
  if (!fVertex || fVertex->GetNContributors() < 1) {
    // SPD vertex
    fVertex = (AliESDVertex*)event->GetPrimaryVertexSPD();
    if (!fVertex || fVertex->GetNContributors() < 1)
      fVertex = 0x0;
  }
  if (!fVertex)
    return kFALSE;
  if (TMath::Abs(fVertex->GetZ()) < 10.0)
    return kTRUE;
  else
    return kFALSE;
}

//-----------------------------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::IsTPCTOFMatched(AliESDtrack* track, Bool_t checkMatchLabel)
{
  //defines TOF matching
  if (!track) {
    AliWarning("Invalid track object");
    return kFALSE;
  }

  Bool_t isMatched = kFALSE;
  if ((track->IsOn(AliESDtrack::kTOFout)) && (track->IsOn(AliESDtrack::kTIME)) && (track->IsOn(AliESDtrack::kTPCout)))
    isMatched = kTRUE;

  if (checkMatchLabel && fIsMC) {
    const Int_t TrkLabel = track->GetLabel(); /*The Get*Label() getters return the label of the associated MC particle. The absolute value of this label is the index of the particle within the MC fMCStack. If the label is negative, this track was assigned a certain number of clusters that did not in fact belong to this track. */
    const Int_t AbsTrkLabel = TMath::Abs(TrkLabel);
    Int_t TOFTrkLabel[3] = { -1 }; //This can contain three particles wich occupy the same cluster
    // Int_t mfl, uniqueID;
    //Gets the labels of the tracks matched to the TOF,
    //this can be used to remove the mismatch and to compute the efficiency!
    //The label to check is the first one, the others can come from different tracks
    track->GetTOFLabel(TOFTrkLabel);
    if (TOFTrkLabel[0] == -1) {
      fMCTOFMatch = -1;
      isMatched &= kFALSE;
    } // Track was not matched to any TOF hit.
    else if (AbsTrkLabel == TOFTrkLabel[0]) {
      fMCTOFMatch = 0;
      isMatched &= kTRUE;
    } // Track was correctly matched to a TOF hit.
    else {
      fMCTOFMatch = 1;
      isMatched &= kFALSE;
    } // Track was matched to a TOF hit but comes from mismatch!
  }
  return isMatched;
}

//-----------------------------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::IsInTRD(AliESDtrack* track)
{
  //defines cut to select particles in/out TRD
  if (!track) {
    AliWarning("Invalid track object");
    return kFALSE;
  }

  if (track->IsOn(AliESDtrack::kTPCout)
      && track->IsOn(AliESDtrack::kTRDout))
    return kTRUE;
  else
    return kFALSE;
}
//-----------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillStartTimeMaskHisto(TString suffix)
{
  /* set pid response to use best_T0 and for each
     accepted track fills the histogram with the
     used start time
  */

  //set response best_t0
  //fESDpid->SetTOFResponse(fESD,AliESDpid::kBest_T0);

  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      AliInfo(Form("Cannot receive track %d", iTracks));
      continue;
    }
    //primary tracks selection: kTPCrefit and std cuts
    if (fTrackFilter) {
      if (!fTrackFilter->IsSelected(track))
        continue;
    } else {
      AliInfo("No track filter found, skipping the track loop");
      break;
    }
    if (TMath::Abs(track->Eta()) > fMatchingEtaCut)
      continue; //cut for acceptance

    Int_t StartTimeBit = fESDpid->GetTOFResponse().GetStartTimeMask(track->P());
    ((TH2F*)fHlistTimeZero->FindObject(Form("hStartTimeMask%s", suffix.Data())))->Fill(track->P(), StartTimeBit);

    //matched tracks selection: kTOFout and kTIME
    if ((track->IsOn(AliESDtrack::kTOFout)) && (track->IsOn(AliESDtrack::kTIME)) && (track->IsOn(AliESDtrack::kTPCout))) {
      ((TH2F*)fHlistTimeZero->FindObject(Form("hStartTimeMaskMatched%s", suffix.Data())))->Fill(track->P(), StartTimeBit);
    }
  }
  ((TH2F*)fHlistTimeZero->FindObject(Form("hStartTimeMaskvsTOFmulti%s", suffix.Data())))->Fill(fNTOFtracks[0], fESDpid->GetTOFResponse().GetStartTimeMask(10.));
  return;
}

//----------------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::ComputeTimeZeroByTOF1GeV()
{
  /* compute T0-TOF for tracks within momentum range [fResolutionMinP, fResolutionMaxP] */
  /* init T0-TOF */
  AliTOFT0v1 TOFT0v1(fESDpid); // TOF-T0 v1
  TOFT0v1.Init(fESD);
  TOFT0v1.DefineT0("all", fResolutionMinP, fResolutionMaxP);
  fMyTimeZeroTOF.SetAt(-1000. * TOFT0v1.GetResult(0), 0);
  fMyTimeZeroTOFsigma.SetAt(1000. * TOFT0v1.GetResult(1), 0);
  fMyTimeZeroTOFtracks.SetAt(TOFT0v1.GetResult(3), 0);
  if (fUseTOFT0CalibMode) {
    //Get TOFT0 with first half of the tracks
    Int_t mode = 1;
    TOFT0v1.DefineT0("all", fResolutionMinP, fResolutionMaxP, mode);
    fMyTimeZeroTOF.SetAt(-1000. * TOFT0v1.GetResult(0), mode);
    fMyTimeZeroTOFsigma.SetAt(1000. * TOFT0v1.GetResult(1), mode);
    fMyTimeZeroTOFtracks.SetAt(TOFT0v1.GetResult(3), mode);
    //Get TOFT0 with second half of the tracks
    mode = 2;
    TOFT0v1.DefineT0("all", fResolutionMinP, fResolutionMaxP, mode);
    fMyTimeZeroTOF.SetAt(-1000. * TOFT0v1.GetResult(0), mode);
    fMyTimeZeroTOFsigma.SetAt(1000. * TOFT0v1.GetResult(1), mode);
    fMyTimeZeroTOFtracks.SetAt(TOFT0v1.GetResult(3), mode);
  }
  fMyTimeZeroTOFstatus = kFALSE;
  /* check T0-TOF sigma */
  if (fMyTimeZeroTOFsigma.At(0) < 250.)
    fMyTimeZeroTOFstatus = kTRUE;
  return fMyTimeZeroTOFstatus;
}

//------------------------------------------------------
TString AliAnalysisTaskTOFqaID::GetSpeciesName(Int_t absPdgCode)
{
  //returns name of selected specie
  TString name;
  switch (absPdgCode) {
  case 11:
    name = "electron";
    break;
  case 13:
    name = "muon";
    break;
  case 211:
    name = "pion";
    break;
  case 321:
    name = "kaon";
    break;
  case 2212:
    name = "proton";
    break;
  default:
    name = "noPID";
    break;
  }
  return name.Data();
}

//-----------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::SelectMCspecies(AliMCEvent* ev, AliESDtrack* track)
{
  //
  //Retrieves particle true ID from MC and selects the desired species
  //
  if ((!ev) || (!track)) {
    AliError("SelectMCspecies - Invalid object set as argument");
    return kFALSE;
  }

  if (fSelectedPdg == 0)
    return kTRUE; //if fSelectedPdg==0, no species selection is applied

  Long_t label = track->GetLabel();
  if (label < 0)
    return kFALSE;

  // get number of particles
  Long_t nMC = ev->GetNumberOfTracks();
  // if label too large --> failed
  if (label >= nMC) {
    AliWarning(Form("Stack overflow: track label = %li -- stack maximum = %li", label, nMC));
    return kFALSE;
  }
  // retrieve particle
  AliMCParticle* mcPart = (AliMCParticle*)ev->GetTrack(label);
  if (!mcPart) { // if particle = NULL --> failed
    AliWarning(Form("Stack discontinuity: label %li refers to a NULL object", label));
    return kFALSE;
  }

  Int_t pdgCode = mcPart->PdgCode();
  if (!(TMath::Abs(pdgCode) == fSelectedPdg))
    return kFALSE;
  else
    return kTRUE;
}

//----------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::ComputeMatchingEfficiency(THashList* list, TString variable)
{
  //computes matching efficiency from previously filled histos
  // to be called in terminate function
  if (!list)
    return kFALSE;

  TString matchedName, primaryName, xAxisTitle;
  if (variable.Contains("pt")) {
    matchedName = "hTOFmatchedESDPt";
    primaryName = "hESDprimaryTrackPt";
    xAxisTitle = "#it{p}_{T} (GeV/#it{c})";
  }
  if (variable.Contains("eta")) {
    matchedName = "hTOFmatchedESDeta";
    primaryName = "hTOFprimaryESDeta";
    xAxisTitle = "#eta";
  }
  if (variable.Contains("phi")) {
    matchedName = "hTOFmatchedESDphi";
    primaryName = "hTOFprimaryESDphi";
    xAxisTitle = "#phi_vtx (deg)";
  }

  TH1F* hDummy = ((TH1F*)list->FindObject(matchedName.Data()));
  if (!hDummy)
    return 0;

  TH1F* hMatchingEff = (TH1F*)hDummy->Clone("hMatchingEff");
  hMatchingEff->SetNameTitle(Form("hMatchingEff_%s", variable.Data()), Form("Matching efficiency vs %s", variable.Data()));
  hMatchingEff->Divide((TH1F*)list->FindObject(primaryName.Data()));
  hMatchingEff->GetXaxis()->SetTitle(xAxisTitle.Data());
  hMatchingEff->GetYaxis()->SetRangeUser(0.0, 1.0);
  hMatchingEff->GetYaxis()->SetTitle("#epsilon_{match}");
  list->AddLast(hMatchingEff);
  return 1;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::HistogramMakeUp(TH1* hist, Color_t color, Int_t markerStyle)
{
  //set histogram style and axes style at once
  if (!hist)
    return;
  if (color >= 0) {
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
  }
  if (markerStyle >= 0)
    hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(0.7);
  //hist->Sumw2();
  return;
}

#define GetArrayBinning(arr) arr.GetSize() - 1, arr.GetArray()
#define CreateH(Hname, Type, title, ...)                                                                                           \
  Type* Hname = new Type(Form("%s%s_%s", #Hname, suffix.Data(), cLabel.Data()), Form("%s %s", cLabel.Data(), title), __VA_ARGS__); \
  HistogramMakeUp(Hname, ((charge > 0) ? kRed + 2 : kBlue + 2), 1);                                                                \
  list->AddLast(Hname);

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::AddTofBaseHisto(THashList* list, Int_t charge, TString suffix)
{
  //Creates histograms for monitoring TOF signal, time alignement and matching-related quantities
  if (!list) {
    AliError("Invalid list passed as argument.");
    return;
  }

  TString cLabel = "all";
  if (charge < 0)
    cLabel.Form("neg");
  else if (charge > 0)
    cLabel.Form("pos");

  CreateH(hTOFmulti, TH1I, Form("matched trk per event (|#eta|#leq%3.2f, #it{p}_{T}#geq0.3 GeV/#it{c});N;events", fMatchingEtaCut), GetArrayBinning(fVariableBinsMult));

  CreateH(hTime, TH1F, "matched trk TOF signal;t (ns);tracks", 250, 0., 610.);

  CreateH(hRawTime, TH1F, "matched trk TOF raw signal;t_{raw} (ns);tracks", 250, 0., 610.);

  CreateH(hTot, TH1F, "matched trk ToT;ToT (ns);tracks", 50, 0., 50.);

  CreateH(hMatchedL, TH1F, "matched trk lenght;L (cm);tracks", 900, -100., 800);

  CreateH(hMatchedDxVsPt, TH2F, "matched trk dx vs. #it{p}_{T};#it{p}_{T} (GeV/#it{c});dx (cm)", GetArrayBinning(fVariableBinsPt), 200, -10., 10.);

  CreateH(hMatchedDzVsStrip, TH2F, "matched trk dz vs. strip (#eta);strip index;dz (cm)", 92, 0., 92., 200, -10., 10.);

  CreateH(hMatchedDxVsCh, TProfile, "matched trk dx vs. channel;channel index;dx (cm)", 157248., 0., 157248.);

  CreateH(hMatchedDzVsCh, TProfile, "matched trk dz vs. channel;channel index;dz (cm)", 157248., 0., 157248.);

  return;
}

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::AddMatchingEffHisto(THashList* list, Int_t charge, TString suffix)
{
  if (!list) {
    AliError("Invalid list passed as argument.");
    return;
  }

  TString cLabel = "all";
  if (charge < 0)
    cLabel = "neg";
  else if (charge > 0)
    cLabel = "pos";

  CreateH(hMatchedP, TH1F, "matched trk p;p (GeV/#it{c});tracks", GetArrayBinning(fVariableBinsPt));

  CreateH(hMatchedPt, TH1F, "matched trk #it{p}_{T};#it{p}_{T} (GeV/#it{c});tracks", GetArrayBinning(fVariableBinsPt));

  CreateH(hMatchedPhi, TH1F, "matched trk #phi_{vtx};#phi_{vtx} (deg);tracks", fnBinsPhi, fBinsPhi[0], fBinsPhi[1]);

  CreateH(hMatchedPtVsOutPhi, TH2F, "matched trk #it{p}_{T} vs. #phi_{TPC out};#phi_{TPC out} (deg);#it{p}_{T} (GeV/#it{c})", fnBinsPhi, fBinsPhi[0], fBinsPhi[1], GetArrayBinning(fVariableBinsPt));

  CreateH(hMatchedEtaVsOutPhi, TH2F, "matched trk #eta vs. #phi_{TPC out};#phi_{TPC out} (deg);#eta", fnBinsPhi, fBinsPhi[0], fBinsPhi[1], fnBinsEta, fBinsEta[0], fBinsEta[1]);

  CreateH(hPrimaryP, TH1F, "primary trk p;p (GeV/#it{c});tracks", GetArrayBinning(fVariableBinsPt));

  CreateH(hPrimaryPt, TH1F, "primary trk #it{p}_{T};#it{p}_{T} (GeV/#it{c});tracks", GetArrayBinning(fVariableBinsPt));

  CreateH(hPrimaryPhi, TH1F, "primary trk #phi_{vtx};#phi_{vtx} (deg);tracks", fnBinsPhi, fBinsPhi[0], fBinsPhi[1]);

  CreateH(hPrimaryPtVsOutPhi, TH2F, "primary trk #it{p}_{T} vs. #phi_{TPC out};#phi_{TPC out} (deg);#it{p}_{T} (GeV/#it{c})", fnBinsPhi, fBinsPhi[0], fBinsPhi[1], GetArrayBinning(fVariableBinsPt));

  CreateH(hPrimaryEtaVsOutPhi, TH2F, "primary trk #eta vs. #phi_{TPC out}; #eta; #phi_{TPC out};#phi_{TPC out} (deg);#eta", fnBinsPhi, fBinsPhi[0], fBinsPhi[1], fnBinsEta, fBinsEta[0], fBinsEta[1]);
}

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::AddPidHisto(THashList* list, Int_t charge, TString suffix)
{
  //Creates histograms for monitoring TOF PID
  if (!list) {
    AliError("Invalid list passed as argument.");
    return;
  }

  TString cLabel = "all";
  if (charge < 0)
    cLabel.Form("neg");
  else if (charge > 0)
    cLabel.Form("pos");

  CreateH(hMatchedBetaVsP, TH2F, "matched trk #beta vs. p;#it{p} (GeV/#it{c});#beta", GetArrayBinning(fVariableBinsPt), 150, 0., 1.5);

  CreateH(hMatchedMass, TH1F, "matched p.le M;M (GeV/#it{c}^{2});entries", 500, 0., 5.);

  CreateH(hMatchedMass2, TH1F, "matched p.le M^{2};M^{2} (GeV^{2}/c^{4});entries", 500, 0., 10.);

  CreateH(hExpTimePiVsStrip, TH2F, "matched trk t_{TOF}-t_{#pi,exp} vs strip;strip (#eta);t_{TOF}-t_{#pi,exp} (ps)", 92, 0, 92, fnExpTimeSmallBins, fExpTimeSmallRangeMin, fExpTimeSmallRangeMax);

  CreateH(hExpTimePiT0Sub1GeV, TH2F, Form("trk (%.2f#leq #it{p}#leq %.2f GeV/#it{c}) t_{TOF}-t_{#pi,exp}-t_{0}^{TOF};n. tracks used for t_{0}^{TOF};t_{TOF}-t_{#pi,exp}-t_{0}^{TOF}", fResolutionMinP, fResolutionMaxP), GetArrayBinning(fVariableBinsMult), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);

  CreateH(hExpTimePiFillSub, TH1F, "trk t_{TOF}-t_{#pi,exp}-t_{0,fill};t_{TOF}-t_{#pi,exp} -t_{0,fill} (ps);entries", 6150, -75030., 75030.);

  CreateH(hExpTimePi, TH1F, "matched trk t_{TOF}-t_{#pi,exp};t_{TOF}-t_{#pi,exp} (ps);tracks", 6150, -75030., 75030.);

  CreateH(hExpTimePiVsP, TH2F, "matched trk t_{TOF}-t_{#pi,exp};#it{p} (GeV/#it{c});t_{TOF}-t_{#pi,exp} (ps)", GetArrayBinning(fVariableBinsPt), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimePiVsP, kRed + 2);

  CreateH(hExpTimeKaVsP, TH2F, "matched trk t_{TOF}-t_{K,exp};#it{p} (GeV/#it{c});t_{TOF}-t_{K,exp} (ps)", GetArrayBinning(fVariableBinsPt), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimeKaVsP, kBlue + 2);

  CreateH(hExpTimeProVsP, TH2F, "matched trk t_{TOF}-t_{p,exp};#it{p} (GeV/#it{c});t_{TOF}-t_{p,exp} (ps)", GetArrayBinning(fVariableBinsPt), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimeProVsP, kGreen + 2);

  CreateH(hTOFpidSigmaPi, TH2F, "trk n#sigma^{TOF}_{#pi} vs #it{p}_{T};#it{p} (GeV/#it{c});n#sigma_{#pi,exp} (ps)", 500, 0., 5., 200, -10., 10.);
  HistogramMakeUp(hTOFpidSigmaPi, kRed + 2);

  CreateH(hTOFpidSigmaKa, TH2F, "trk n#sigma^{TOF}_{K} vs #it{p}_{T};#it{p} (GeV/#it{c});n#sigma_{K,exp} (ps)", 500, 0., 5., 200, -10., 10.);
  HistogramMakeUp(hTOFpidSigmaKa, kBlue + 2);

  CreateH(hTOFpidSigmaPro, TH2F, "trk TOF n#sigma^{TOF}_{p} vs #it{p}_{T};#it{p} (GeV/#it{c});n#sigma_{p,exp} (ps)", 500, 0., 5., 200, -10., 10.);
  HistogramMakeUp(hTOFpidSigmaPro, kGreen + 2);

  CreateH(hExpTimePiT0SubVsP, TH2F, "trk t_{TOF}-t_{#pi,exp}-t_{0}^{TOF};#it{p} (GeV/#it{c});t_{TOF}-t_{#pi,exp}-t_{0}^{TOF}", GetArrayBinning(fVariableBinsPt), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimePiT0SubVsP, kRed + 2);

  CreateH(hExpTimeKaT0SubVsP, TH2F, "trk t_{TOF}-t_{K,exp}-t_{0}^{TOF};#it{p} (GeV/#it{c});t_{TOF}-t_{K,exp}-t_{0}^{TOF}", GetArrayBinning(fVariableBinsPt), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimeKaT0SubVsP, kBlue + 2);

  CreateH(hExpTimeProT0SubVsP, TH2F, "trk t_{TOF}-t_{p,exp}-t_{0}^{TOF};#it{p} (GeV/#it{c});t_{TOF}-t_{p,exp}-t_{0}^{TOF}", GetArrayBinning(fVariableBinsPt), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimeProT0SubVsP, kGreen + 2);

  CreateH(hExpTimePiVsOutPhi, TH2F, "matched trk t_{TOF}-t_{#pi,exp} vs #phi_{TPC out};#phi_{TPC out} (deg);t_{TOF}-t_{#pi,exp} (ps)", fnBinsPhi, fBinsPhi[0], fBinsPhi[1], fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimePiVsOutPhi, kRed + 2);

  CreateH(hExpTimeKaVsOutPhi, TH2F, "matched trk t_{TOF}-t_{K,exp} vs #phi_{TPC out};#phi_{TPC out} (deg);t_{TOF}-t_{K,exp} (ps)", fnBinsPhi, fBinsPhi[0], fBinsPhi[1], fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimeKaVsOutPhi, kBlue + 2);

  CreateH(hExpTimeProVsOutPhi, TH2F, "matched trk t_{TOF}-t_{p,exp} vs #phi_{TPC out};#phi_{TPC out} (deg);t_{TOF}-t_{p,exp} (ps)", fnBinsPhi, fBinsPhi[0], fBinsPhi[1], fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimeProVsOutPhi, kGreen + 2);

  CreateH(hExpTimePiVsPgoodCh, TH2F, "matched trk t_{TOF}-t_{#pi,exp} - good channels;#it{p} (GeV/#it{c});t_{TOF}-t_{#pi,exp} (ps)", GetArrayBinning(fVariableBinsPt), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimePiVsPgoodCh, kRed + 2);

  CreateH(hExpTimeKaVsPgoodCh, TH2F, "matched trk t_{TOF}-t_{K,exp} - good channels;#it{p} (GeV/#it{c});t_{TOF}-t_{K,exp} (ps)", GetArrayBinning(fVariableBinsPt), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimeKaVsPgoodCh, kBlue + 2);

  CreateH(hExpTimeProVsPgoodCh, TH2F, "matched trk t_{TOF}-t_{p,exp} - good channels;#it{p} (GeV/#it{c});t_{TOF}-t_{p,exp} (ps)", GetArrayBinning(fVariableBinsPt), fnExpTimeBins, fExpTimeRangeMin, fExpTimeRangeMax);
  HistogramMakeUp(hExpTimeProVsPgoodCh, kGreen + 2);

  return;
}
#undef CreateH
#define CreateH(Hname, Type, title, ...)                                           \
  Type* Hname = new Type(Form("%s%s", #Hname, suffix.Data()), title, __VA_ARGS__); \
  list->AddLast(Hname);

#define SetBinLabels(H)                                \
  if (labels_arr->GetEntries() != H->GetNbinsY())      \
    AliFatal("Not the right number of labels");        \
  for (Int_t i = 0; i < labels_arr->GetEntries(); i++) \
    H->GetYaxis()->SetBinLabel(i + 1, labels_arr->At(i)->GetName());
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::AddStartTimeHisto(THashList* list, TString suffix)
{
  //Creates histograms for monitoring T0 signal and start-time related quantities
  if (!list) {
    AliError("Invalid list passed as argument.");
    return;
  }
  CreateH(hT0AC, TH1F, "Event timeZero from T0A&C; t_{0,AC} (ps); events", 1000, -12500., 12500.);
  HistogramMakeUp(hT0AC, kRed + 2, 20);

  CreateH(hT0A, TH1F, "Event timeZero from T0A; t_{0,A} (ps); events", 1000, -12500., 12500.);
  HistogramMakeUp(hT0A, kBlue + 2, 25);

  CreateH(hT0C, TH1F, "Event timeZero from T0C; t_{0,C} (ps); events", 1000, -12500., 12500.);
  HistogramMakeUp(hT0C, kGreen + 2, 28);

  CreateH(hT0DetRes, TH1F, "T0 detector (T0A-T0C)/2; (T0A-T0C)/2 (ps); events", 200, -500., 500.);
  HistogramMakeUp(hT0DetRes, kMagenta + 1, 1);

  CreateH(hEventT0MeanVsVtx, TH2F, "T0 detector: mean vs vertex ; (t0_{A}-t0_{C})/2 [ns]; (t0_{A}+t0_{C})/2 [ns]; events", 50, -25., 25., 50, -25., 25.);
  HistogramMakeUp(hEventT0MeanVsVtx, kBlue + 2, 1);

  CreateH(hEventV0MeanVsVtx, TH2F, "V0 detector: mean vs vertex ; (V0_{A}-V0_{C})/2 [ns]; (V0_{A}+V0_{C})/2 [ns]; events", 50, -25., 25., 50, -25., 25.);
  HistogramMakeUp(hEventV0MeanVsVtx, kBlack, 1);

  TString labels = "best_t0,fill_t0,tof_t0,T0AC,T0A,T0C";
  TObjArray* labels_arr = labels.Tokenize(",");

  CreateH(hStartTime, TH2F, "Start time for each method (mask); start time (ps); method;", fnBinsT0, fBinsT0[0], fBinsT0[1], 6, -1., 5.);
  SetBinLabels(hStartTime);

  CreateH(hStartTimeRes, TH2F, "Start time resolution for each method (mask); resolution (ps); method;", 300, 0., 300., 6, -1., 5.);
  SetBinLabels(hStartTimeRes);

  CreateH(hT0TOFvsNtrk, TH2F, "Event timeZero estimated by TOF vs. TOF-matching tracks; N_{TOF}; t0 (ps)", GetArrayBinning(fVariableBinsMult), fnBinsT0, fBinsT0[0], fBinsT0[1]);
  HistogramMakeUp(hT0TOFvsNtrk, kTeal - 5, 1);

  CreateH(hT0ACvsNtrk, TH2F, "Event timeZero estimated by T0AC vs. TOF-matching tracks; N_{TOF}; t0 (ps)", GetArrayBinning(fVariableBinsMult), fnBinsT0, fBinsT0[0], fBinsT0[1]);
  HistogramMakeUp(hT0ACvsNtrk, kRed + 2, 1);

  CreateH(hT0AvsNtrk, TH2F, "Event timeZero estimated by T0A vs. TOF-matching tracks; N_{TOF}; t0 (ps)", GetArrayBinning(fVariableBinsMult), fnBinsT0, fBinsT0[0], fBinsT0[1]);
  HistogramMakeUp(hT0AvsNtrk, kBlue + 2, 1);

  CreateH(hT0CvsNtrk, TH2F, "Event timeZero estimated by T0C vs. TOF-matching tracks; N_{TOF}; t0 (ps)", GetArrayBinning(fVariableBinsMult), fnBinsT0, fBinsT0[0], fBinsT0[1]);
  HistogramMakeUp(hT0CvsNtrk, kGreen + 2, 1);

  labels = "fill_t0,tof_t0,T0A,T0A & tof_t0,T0C,T0C & tof_t0,T0AC,T0AC & tof_t0";
  labels_arr = labels.Tokenize(",");
  const Double_t startTimeMomBins[13] = { 0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.5, 2., 3., 10. };
  CreateH(hStartTimeMaskMatched, TH2F, "Start Time Mask vs p bin for matched tracks; p(GeV/#it{c});", 12, startTimeMomBins, 8, 0., 8.);
  SetBinLabels(hStartTimeMaskMatched);
  HistogramMakeUp(hStartTimeMaskMatched, kRed + 2, 1);

  CreateH(hStartTimeMask, TH2F, "Start Time Mask vs p bin for primary tracks; p(GeV/#it{c});", 12, startTimeMomBins, 8, 0., 8.);
  SetBinLabels(hStartTimeMask);
  HistogramMakeUp(hStartTimeMask, kRed + 2, 1);

  CreateH(hStartTimeMaskvsTOFmulti, TH2F, "Start Time Mask vs TOF hit multiplicity; TOF hit multiplicity;", GetArrayBinning(fVariableBinsMult), 8, 0., 8.);
  SetBinLabels(hStartTimeMaskvsTOFmulti);
  HistogramMakeUp(hStartTimeMaskvsTOFmulti, -1, -1);

  if (fUseTOFT0CalibMode) {
    CreateH(hT0TOFdiffvsNtrk, TH2F, "Event timeZero estimated by TOF (first half - second half) vs. TOF-matching tracks; n. tracks used for t_{0}^{TOF} (average f. h. and s. h.); TOFt0_{f. h.} - TOFt0_{s. h.} (ps)", GetArrayBinning(fVariableBinsMult), fnBinsT0, fBinsT0[0], fBinsT0[1]);
    HistogramMakeUp(hT0TOFdiffvsNtrk, -1, -1);

    CreateH(hT0TOFdiffNormvsNtrk, TH2F, "Event timeZero estimated by TOF (first half - second half) normalized to reso. vs. TOF-matching tracks; n. tracks used for t_{0}^{TOF} (average f. h. and s. h.); (TOFt0_{f. h.} - TOFt0_{s. h.})/#sqrt{#sigma_{TOFt0_{f. h.}}^{2} -#sigma_{TOFt0_{s. h.}}^{2}}", GetArrayBinning(fVariableBinsMult), fnBinsT0, -10, 10);
    HistogramMakeUp(hT0TOFdiffNormvsNtrk, -1, -1);
  }
}
#undef SetBinLabels
#undef CreateH
#undef GetArrayBinning

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::AddTrdHisto()
{
  //Creates histograms for monitoring TOF base quantities wrt TRD/no TRD selection
  if (!fHlistTRD) {
    AliError("Invalid TRD list");
    return;
  }

  if (fEnableChargeSplit) {
    AddMatchingEffHisto(fHlistTRD, 1, "_noTrd");
    AddMatchingEffHisto(fHlistTRD, -1, "_noTrd");
    AddMatchingEffHisto(fHlistTRD, 1, "_Trd");
    AddMatchingEffHisto(fHlistTRD, -1, "_Trd");

    AddPidHisto(fHlistTRD, 1, "_noTrd");
    AddPidHisto(fHlistTRD, -1, "_noTrd");
    AddPidHisto(fHlistTRD, 1, "_Trd");
    AddPidHisto(fHlistTRD, -1, "_Trd");
  } else {
    AddMatchingEffHisto(fHlistTRD, 0, "_noTrd");
    AddMatchingEffHisto(fHlistTRD, 0, "_Trd");
    AddPidHisto(fHlistTRD, 0, "_noTrd");
    AddPidHisto(fHlistTRD, 0, "_Trd");
  }

  return;
}

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::AddTofTrgHisto(TString suffix)
{
  //defines histo with trigger info
  if (!fHlistTrigger) {
    AliError("Invalid TOF trigger list");
    return;
  }

  TH1I* hFiredMaxipad = new TH1I("hFiredMaxipad" + suffix, "Fired maxipad per event;N_{maxipad};events", 1584, 0, 1584);
  HistogramMakeUp(hFiredMaxipad, kBlue + 2, 1);
  fHlistTrigger->AddLast(hFiredMaxipad);

  TH1I* hFiredReadoutPad = new TH1I("hFiredReadoutPad" + suffix, "Fired readout pad per event;N_{pad};events", 153000, 0, 153000);
  HistogramMakeUp(hFiredReadoutPad, kRed + 2, 1);
  fHlistTrigger->AddLast(hFiredReadoutPad);

  TH1I* hFiredReadoutTrgPad = new TH1I("hFiredReadoutTrgPad" + suffix, "Fired readout pad in trg window;N_{pad} in trg window;events", 153000, 0, 153000);
  HistogramMakeUp(hFiredReadoutTrgPad, kBlack, 1);
  fHlistTrigger->AddLast(hFiredReadoutTrgPad);

  TH2I* hFiredMaxipadVsTrgPad = new TH2I("hFiredMaxipadVsTrgPad" + suffix, "Fired maxipad vs pads in trg window per event;N_{pad} in trg window;N_{maxipad}", 100, 0, 100, 100, 0, 100);
  HistogramMakeUp(hFiredMaxipadVsTrgPad, kBlue + 2, 1);
  fHlistTrigger->AddLast(hFiredMaxipadVsTrgPad);

  TH2I* hTrgMap = new TH2I("hTrgMap" + suffix, "Map of fired maxipads;LTM;maxipad", 72, 0, 72, 23, 0, 23);
  HistogramMakeUp(hTrgMap, kBlue + 2, 1);
  fHlistTrigger->AddLast(hTrgMap);

  return;
}

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillTofBaseHisto(AliESDtrack* track, Int_t charge, TString suffix)
{
  //fill histo with TOF base quantities
  if (!track)
    return;

  //  Double_t tofTime=track->GetTOFsignal();//in ps
  Double_t tofTimeRaw = track->GetTOFsignalRaw(); //in ps
  Double_t tofToT = track->GetTOFsignalToT();     //in ps
  Int_t channel = track->GetTOFCalChannel();
  Int_t volId[5]; //(sector, plate,strip,padZ,padX)
  AliTOFGeometry::GetVolumeIndices(channel, volId);

  TString cLabel = "all";
  if (charge < 0)
    cLabel.Form("neg");
  else if (charge > 0)
    cLabel.Form("pos");

  ((TH1F*)fHlist->FindObject(Form("hTime%s_%s", suffix.Data(), cLabel.Data())))->Fill(fTof);                 //ns
  ((TH1F*)fHlist->FindObject(Form("hRawTime%s_%s", suffix.Data(), cLabel.Data())))->Fill(tofTimeRaw * 1E-3); //ns
  ((TH1F*)fHlist->FindObject(Form("hTot%s_%s", suffix.Data(), cLabel.Data())))->Fill(tofToT);
  ((TH1F*)fHlist->FindObject(Form("hMatchedL%s_%s", suffix.Data(), cLabel.Data())))->Fill(fL);
  ((TH2F*)fHlist->FindObject(Form("hMatchedDxVsPt%s_%s", suffix.Data(), cLabel.Data())))->Fill(fPt, track->GetTOFsignalDx());
  ((TH2F*)fHlist->FindObject(Form("hMatchedDzVsStrip%s_%s", suffix.Data(), cLabel.Data())))->Fill((Int_t)GetStripIndex(volId), track->GetTOFsignalDz());
  ((TProfile*)fHlist->FindObject(Form("hMatchedDxVsCh%s_%s", suffix.Data(), cLabel.Data())))->Fill(channel, track->GetTOFsignalDx());
  ((TProfile*)fHlist->FindObject(Form("hMatchedDzVsCh%s_%s", suffix.Data(), cLabel.Data())))->Fill(channel, track->GetTOFsignalDz());

  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillPrimaryTrkHisto(Int_t charge, TString suffix)
{
  // fill histos with primary tracks info
  // => denominator for matching efficiency

  TString cLabel = "all";
  if (charge < 0)
    cLabel.Form("neg");
  else if (charge > 0)
    cLabel.Form("pos");

  THashList* theL = (suffix.Contains("Trd") ? fHlistTRD : fHlist);
  ((TH1F*)theL->FindObject(Form("hPrimaryP%s_%s", suffix.Data(), cLabel.Data())))->Fill(fP);
  ((TH1F*)theL->FindObject(Form("hPrimaryPt%s_%s", suffix.Data(), cLabel.Data())))->Fill(fPt);
  ((TH2F*)theL->FindObject(Form("hPrimaryPtVsOutPhi%s_%s", suffix.Data(), cLabel.Data())))->Fill(fTPCOuterPhi, fPt);
  if (fPt >= fMatchingMomCut) {
    ((TH1F*)theL->FindObject(Form("hPrimaryPhi%s_%s", suffix.Data(), cLabel.Data())))->Fill(fPhi);
    ((TH2F*)theL->FindObject(Form("hPrimaryEtaVsOutPhi%s_%s", suffix.Data(), cLabel.Data())))->Fill(fTPCOuterPhi, fEta);
  }

  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillMatchedTrkHisto(Int_t charge, TString suffix)
{
  //get matched tracks variables (matching cut to be applied externally)
  //=> numerator for matching efficiency

  TString cLabel = "all";
  if (charge < 0)
    cLabel.Form("neg");
  else if (charge > 0)
    cLabel.Form("pos");

  THashList* theL = (suffix.Contains("Trd") ? fHlistTRD : fHlist);
  ((TH1F*)theL->FindObject(Form("hMatchedP%s_%s", suffix.Data(), cLabel.Data())))->Fill(fP);
  ((TH1F*)theL->FindObject(Form("hMatchedPt%s_%s", suffix.Data(), cLabel.Data())))->Fill(fPt);
  ((TH2F*)theL->FindObject(Form("hMatchedPtVsOutPhi%s_%s", suffix.Data(), cLabel.Data())))->Fill(fTPCOuterPhi, fPt);
  if (fPt >= fMatchingMomCut) {
    ((TH1F*)theL->FindObject(Form("hMatchedPhi%s_%s", suffix.Data(), cLabel.Data())))->Fill(fPhi);
    ((TH2F*)theL->FindObject(Form("hMatchedEtaVsOutPhi%s_%s", suffix.Data(), cLabel.Data())))->Fill(fTPCOuterPhi, fEta);
  }
  return;
}

//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillPidHisto(AliESDtrack* track, Int_t charge, TString suffix)
{
  //basic PID performance check
  if (fTof <= 0) {
    if (fVerbose)
      printf("WARNING: track with negative TOF time found! Skipping this track for PID checks\n");
    return;
  }
  if (fL <= 0) {
    if (fVerbose)
      printf("WARNING: track with negative length found!Skipping this track for PID checks\n");
    return;
  }
  if (!track)
    return;

  TString cLabel = "all";
  if (charge < 0)
    cLabel.Form("neg");
  else if (charge > 0)
    cLabel.Form("pos");

  //calculate beta
  Double_t c = TMath::C() * 1.E-9; // m/ns
  Double_t mass = 0.;              //GeV
  Double_t length = fL * 0.01;     // in meters
  Double_t tof = fTof * c;
  Double_t beta = length / tof;
  Double_t fact = (tof / length) * (tof / length) - 1.;
  Double_t fP2 = fP * fP;

  if (fact <= 0) {
    mass = -fP * TMath::Sqrt(-fact);
  } else {
    mass = fP * TMath::Sqrt(fact);
  }

  THashList* theL = (suffix.Contains("Trd") ? fHlistTRD : fHlistPID);
  ((TH2F*)theL->FindObject(Form("hMatchedBetaVsP%s_%s", suffix.Data(), cLabel.Data())))->Fill(fP, beta);
  ((TH1F*)theL->FindObject(Form("hMatchedMass%s_%s", suffix.Data(), cLabel.Data())))->Fill(mass);
  ((TH1F*)theL->FindObject(Form("hMatchedMass2%s_%s", suffix.Data(), cLabel.Data())))->Fill(mass * mass);

  //PID sigmas
  Bool_t isValidBeta[AliPID::kSPECIES] = { 0, 0, 0, 0, 0 };
  for (Int_t specie = 0; specie < AliPID::kSPECIES; specie++) {
    fSigmaSpecie[specie] = fESDpid->GetTOFResponse().GetExpectedSigma(fP, fTrkExpTimes[specie], AliPID::ParticleMass(specie));
    beta = 1 / TMath::Sqrt(1 + AliPID::ParticleMass(specie) * AliPID::ParticleMass(specie) / (fP2));
    if (beta > 0) {
      fThExpTimes[specie] = length * 1.E3 / (beta * c); //ps
      isValidBeta[specie] = kTRUE;
    } else {
      fThExpTimes[specie] = 1E-10;
      isValidBeta[specie] = kFALSE;
    }
  }
  Float_t timeZeroTOF = (Float_t)fESDpid->GetTOFResponse().GetStartTime(fPt);
  Double_t tofps = fTof * 1E3; //ps for t-texp
  Int_t channel = track->GetTOFCalChannel();
  Int_t volId[5]; //(sector, plate,strip,padZ,padX)
  AliTOFGeometry::GetVolumeIndices(channel, volId);
  Char_t partName[3][4] = { "Pi", "Ka", "Pro" };

  //fill histos for pion only
  ((TH2F*)theL->FindObject(Form("hExpTimePiVsStrip%s_%s", suffix.Data(), cLabel.Data())))->Fill((Int_t)GetStripIndex(volId), tofps - fTrkExpTimes[AliPID::kPion]); //ps
  ((TH1F*)theL->FindObject(Form("hExpTimePi%s_%s", suffix.Data(), cLabel.Data())))->Fill(tofps - fTrkExpTimes[AliPID::kPion]);                                     //ps
  if (fMyTimeZeroTOFstatus && (fP > fResolutionMinP) && (fP < fResolutionMaxP)) {
    ((TH2F*)theL->FindObject(Form("hExpTimePiT0Sub1GeV%s_%s", suffix.Data(), cLabel.Data())))->Fill(fMyTimeZeroTOFtracks.At(0), tofps - fMyTimeZeroTOF.At(0) - fTrkExpTimes[AliPID::kPion]);
  }
  //fill sigmas and deltas for each specie
  for (Int_t specie = AliPID::kPion; specie <= AliPID::kProton; specie++) {
    if (isValidBeta[specie]) {
      ((TH2F*)theL->FindObject(Form("hExpTime%sVsP%s_%s", partName[specie - 2], suffix.Data(), cLabel.Data())))->Fill(fP, tofps - fTrkExpTimes[specie]);
      ((TH2F*)theL->FindObject(Form("hTOFpidSigma%s%s_%s", partName[specie - 2], suffix.Data(), cLabel.Data())))->Fill(fPt, (tofps - fTrkExpTimes[specie]) / fSigmaSpecie[specie]);
      ((TH2F*)theL->FindObject(Form("hExpTime%sT0SubVsP%s_%s", partName[specie - 2], suffix.Data(), cLabel.Data())))->Fill(fP, tofps - fTrkExpTimes[specie] - timeZeroTOF);
      ((TH2F*)theL->FindObject(Form("hExpTime%sVsOutPhi%s_%s", partName[specie - 2], suffix.Data(), cLabel.Data())))->Fill(fTPCOuterPhi, tofps - fTrkExpTimes[specie] - timeZeroTOF);

    } // end check on beta
  }

  //re-set response kFILL_T0 to check post-alignment wih OADB
  fESDpid->SetTOFResponse(fESD, AliESDpid::kFILL_T0);                                                                                                 //(fill_t0, tof_t0, t0_t0, best_t0)
  Float_t startTimeFill = fESDpid->GetTOFResponse().GetStartTime(fP);                                                                                 //timeZero for bin pT>10GeV/#it{c}
  ((TH1F*)theL->FindObject(Form("hExpTimePiFillSub%s_%s", suffix.Data(), cLabel.Data())))->Fill(tofps - fTrkExpTimes[AliPID::kPion] - startTimeFill); //ps

  // if (fEnableAdvancedCheck && (fPt<1.)) {
  //   Double_t pos[3]={0.,0.,0.};
  //   track->GetXYZAt(378.,5.,pos);
  //   if ((pos[0]==0.)&&(pos[1]==0.)&&(pos[2]==0.))continue;

  //   Double_t phiTOF=TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
  //   if (phiTOF<0) phiTOF+= (2*TMath::Pi()*TMath::RadToDeg());

  //   //fill t-texp vs phi@TOF
  //    if ((phiOuterTPC<=75) || ((phiOuterTPC>=125)&&(phiOuterTPC<=235)) || (phiOuterTPC>=305) ) { //TRD sectors 2012
  //    if ( ((phiOuterTPC>=85)&&(phiOuterTPC<=115)) || ((phiOuterTPC>=245)&&(phiOuterTPC<=295)) ) {//no TRD sectors 2012
  // }
  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillStartTimeHisto(TString suffix)
{
#define FillH1(Hname, X1) \
  static_cast<TH1F*>(fHlistTimeZero->FindObject(Form("%s%s", Hname, suffix.Data())))->Fill(X1);
#define FillH2(Hname, ...) \
  static_cast<TH2F*>(fHlistTimeZero->FindObject(Form("%s%s", Hname, suffix.Data())))->Fill(__VA_ARGS__);

  //fill start time histo
  if (!fESD) {
    AliError("Invalid event object");
    return;
  }
  // info from V0 detector QA
  AliESDVZERO* vzero = fESD->GetVZEROData();
  Float_t V0Atime = vzero->GetV0ATime();
  Float_t V0Ctime = vzero->GetV0CTime();
  FillH2("hEventV0MeanVsVtx", (V0Atime - V0Ctime) * 0.5, (V0Atime + V0Ctime) * 0.5);

  // info from T0 detector QA
  for (Int_t j = 0; j < 3; j++) {
    fT0[j] = (Float_t)fESD->GetT0TOF(j); //ps
    if (fT0[j] > 90000.)
      fT0[j] = 99999.; //fix old default values to the new one
  }

  Float_t t0cut = 90000.;
  //Float_t t0cut =3 * t0spread; //use this cut to check t0 used in tof response
  // if(t0cut < 500) t0cut = 500;

  if (TMath::Abs(fT0[1]) < t0cut && TMath::Abs(fT0[2]) < t0cut) {
    //&& TMath::Abs(fT0[2]-fT0[1]) < 500)  //add this condition to check t0 used in tof response
    FillH1("hT0DetRes", (fT0[2] - fT0[1]) * 0.5);
    FillH1("hT0AC", fT0[0]);
    FillH2("hEventT0MeanVsVtx", (fT0[2] - fT0[1]) * 0.5e-3, (fT0[2] + fT0[1]) * 0.5e-3);
  }
  if (TMath::Abs(fT0[1]) < t0cut) {
    FillH1("hT0A", fT0[1]);
  }
  if (TMath::Abs(fT0[2]) < t0cut) {
    FillH1("hT0C", fT0[2]);
  }

  //Fill T0fill plots
  (fESDpid->GetTOFResponse()).ResetT0info();
  fESDpid->SetTOFResponse(fESD, AliESDpid::kFILL_T0);
  Double_t timeZero = fESDpid->GetTOFResponse().GetStartTime(10.);
  Double_t timeZeroRes = fESDpid->GetTOFResponse().GetStartTimeRes(10.);
  FillH2("hStartTime", timeZero, "fill_t0", 1);       //will be 0 by definition
  FillH2("hStartTimeRes", timeZeroRes, "fill_t0", 1); //taken from the header

  //fill TOF_T0 plots
  (fESDpid->GetTOFResponse()).ResetT0info();
  fESDpid->SetTOFResponse(fESD, AliESDpid::kTOF_T0);
  Int_t startTimeMaskBin = fESDpid->GetTOFResponse().GetT0binMask(9);
  if (startTimeMaskBin > 0) { //timeZeroMask for 10th bin, corresponding to p>10GeV/#it{c}
    //fill plot only when there is a start time other than fill_t0
    timeZero = fESDpid->GetTOFResponse().GetStartTime(10.);       //timeZero for bin p>10GeV/#it{c}
    timeZeroRes = fESDpid->GetTOFResponse().GetStartTimeRes(10.); //timeZero for bin p>10GeV/#it{c}
    FillH2("hStartTime", timeZero, "tof_t0", 1);
    FillH2("hStartTimeRes", timeZeroRes, "tof_t0", 1);
    FillH2("hT0TOFvsNtrk", fNTOFtracks[0], timeZero);
  }
  if (fUseTOFT0CalibMode
      && (fMyTimeZeroTOFtracks.At(1) > 0 || fMyTimeZeroTOFtracks.At(2) > 0)
      && TMath::Abs((fMyTimeZeroTOFtracks.At(1) - fMyTimeZeroTOFtracks.At(2)) / (fMyTimeZeroTOFtracks.At(1) + fMyTimeZeroTOFtracks.At(2))) < 0.2) {
    FillH2("hT0TOFdiffvsNtrk", (fMyTimeZeroTOFtracks.At(1) + fMyTimeZeroTOFtracks.At(2)) / 2., (fMyTimeZeroTOF.At(1) - fMyTimeZeroTOF.At(2)));
    if (fMyTimeZeroTOFsigma.At(1) > 0 || fMyTimeZeroTOFsigma.At(2) > 0)
      FillH2("hT0TOFdiffNormvsNtrk", (fMyTimeZeroTOFtracks.At(1) + fMyTimeZeroTOFtracks.At(2)) / 2., (fMyTimeZeroTOF.At(1) - fMyTimeZeroTOF.At(2)) / TMath::Sqrt(fMyTimeZeroTOFsigma.At(1) * fMyTimeZeroTOFsigma.At(1) + fMyTimeZeroTOFsigma.At(2) * fMyTimeZeroTOFsigma.At(2)));
  }

  //fill T0_T0 plots
  (fESDpid->GetTOFResponse()).ResetT0info();
  fESDpid->SetTOFResponse(fESD, AliESDpid::kT0_T0);
  startTimeMaskBin = fESDpid->GetTOFResponse().GetT0binMask(9);
  if (startTimeMaskBin > 0) { //timeZeroMask for 10th bin, corresponding to p>10GeV/#it{c}
    //fill plot only when there is a start time other than fill_t0
    timeZero = fESDpid->GetTOFResponse().GetStartTime(10.);       //timeZero for bin p>10GeV/#it{c}
    timeZeroRes = fESDpid->GetTOFResponse().GetStartTimeRes(10.); //timeZero for bin p>10GeV/#it{c}
    if (startTimeMaskBin == 2) {
      FillH2("hStartTime", timeZero, "T0A", 1);
      FillH2("hStartTimeRes", timeZeroRes, "T0A", 1);
      FillH2("hT0AvsNtrk", fNTOFtracks[0], timeZero);
    } else if (startTimeMaskBin == 4) {
      FillH2("hStartTime", timeZero, "T0C", 1);
      FillH2("hStartTimeRes", timeZeroRes, "T0C", 1);
      FillH2("hT0CvsNtrk", fNTOFtracks[0], timeZero);
    } else if (startTimeMaskBin == 6) {
      FillH2("hStartTime", timeZero, "T0AC", 1);
      FillH2("hStartTimeRes", timeZeroRes, "T0AC", 1);
      FillH2("hT0ACvsNtrk", fNTOFtracks[0], timeZero);
    }
  }

  //response set to best_t0
  (fESDpid->GetTOFResponse()).ResetT0info();
  fESDpid->SetTOFResponse(fESD, AliESDpid::kBest_T0);
  //fill plot only when there is a start time other than fill_t0
  timeZero = fESDpid->GetTOFResponse().GetStartTime(10.);       //timeZero for bin p>10GeV/#it{c}
  timeZeroRes = fESDpid->GetTOFResponse().GetStartTimeRes(10.); //timeZero for bin p>10GeV/#it{c}
  FillH2("hStartTime", timeZero, "best_t0", 1);
  FillH2("hStartTimeRes", timeZeroRes, "best_t0", 1);

  //Fill mask when best_t0 is set
  FillStartTimeMaskHisto(suffix.Data());

#undef FillH1
#undef FillH2
  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillTrdHisto(AliESDtrack* track, Int_t charge)
{
  //fill histograms for TRD/noTRD
  if (!track) {
    AliError("Invalid track object");
    return;
  }
  const TString suffix = IsInTRD(track) ? "_Trd" : "_noTrd";
  FillPrimaryTrkHisto(charge, suffix);
  if (IsTPCTOFMatched(track, 0)) {
    FillMatchedTrkHisto(charge, suffix);
    FillPidHisto(track, charge, suffix);
  }
  return;
}
//----------------------------------------------------------------------------------
void AliAnalysisTaskTOFqaID::FillTofTrgHisto(TString suffix)
{
  //fills histo with trigger info
  if (!fHlistTrigger) {
    AliError("Invalid TOF trigger list");
    return;
  }
  if (!fTOFHeader) {
    AliWarning("Invalid AliTOFHeader object - cannot fill trg histo");
    return;
  }

  Int_t nPad = fTOFHeader->GetNumberOfTOFclusters();    //all fired readout pads
  Int_t nTrgPad = fTOFHeader->GetNumberOfTOFtrgPads();  // fired readout pads in the trigger window
  Int_t nMaxiPad = fTOFHeader->GetNumberOfTOFmaxipad(); //fired maxipads

  // update histo with fired macropads
  AliTOFTriggerMask* trgMask = fTOFHeader->GetTriggerMask();
  for (Int_t j = 0; j < 72; j++) {
    for (Int_t i = 22; i >= 0; i--) {
      if (trgMask->IsON(j, i))
        ((TH1I*)fHlistTrigger->FindObject(Form("hTrgMap%s", suffix.Data())))->Fill(j + 1, i + 1);
    }
  }
  ((TH1I*)fHlistTrigger->FindObject(Form("hFiredMaxipad%s", suffix.Data())))->Fill(nMaxiPad);
  ((TH1I*)fHlistTrigger->FindObject(Form("hFiredReadoutPad%s", suffix.Data())))->Fill(nPad);
  ((TH1I*)fHlistTrigger->FindObject(Form("hFiredReadoutTrgPad%s", suffix.Data())))->Fill(nTrgPad);
  ((TH2I*)fHlistTrigger->FindObject(Form("hFiredMaxipadVsTrgPad%s", suffix.Data())))->Fill(nTrgPad, nMaxiPad);

  return;
}

//-----------------------------------------------------------
void AliAnalysisTaskTOFqaID::LoadChannelMapsFromOCDB()
{
  //method to get the channel maps from the OCDB
  // it is called at the CreatedOutputObject stage
  // to comply with the CAF environment

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(fOCDBLocation.Data());
  cdb->SetRun(fRunNumber);

  if (!cdb) {
    AliWarning("No CDB MANAGER, maps can not be loaded");
    return;
  }
  AliCDBPath cdbPath("TOF/Calib/Status");
  AliCDBEntry* cdbe = cdb->Get(cdbPath, fRunNumber);
  if (!cdbe) {
    printf("invalid CDB entry\n");
    fChannelArray = 0x0;
    return;
  }

  fChannelArray = (AliTOFChannelOnlineStatusArray*)cdbe->GetObject();
  fCalib = new AliTOFcalib();
  fCalib->Init();
  return;
}

//-----------------------------------------------------------
Bool_t AliAnalysisTaskTOFqaID::IsChannelGood(Int_t channel)
{
  // methos to check the given channel
  // against the status maps
  // taken from the OCDB

  if (channel < 0)
    return kFALSE;
  if (!fChannelArray || !fCalib) {
    AliError("Array with TOF channel status from OCDB not available.");
    return kFALSE;
  }
  if ((fChannelArray->GetNoiseStatus(channel) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad) || fCalib->IsChannelProblematic(channel))
    return kFALSE;

  return kTRUE;
}

//-----------------------------------------------------------
void AliAnalysisTaskTOFqaID::SetVariableBinning()
{
  //
  // defines an array with variable binning in p and pT
  //
  Int_t size = fVariableBinsPt.GetSize();
  Double_t delta = size > 1 ? fVariableBinsPt.At(0) - fVariableBinsPt.At(1) : 0.;
  if (size > 1 && TMath::Abs(delta) < 0.01) { //Check to see if it was set externally
    for (Int_t j = 0; j < size; j++) {
      if (j < 200)
        fVariableBinsPt.SetAt(0.025 * j, j);
      else
        fVariableBinsPt.SetAt(5.0 + (j - 200) * 0.050, j);
    }
  }
  //Checking that the binning is in order
  for (Int_t j = 0; j < size - 1; j++)
    if (fVariableBinsPt.GetAt(j) > fVariableBinsPt.GetAt(j + 1))
      AliFatal(Form("Binning fVariableBinsPt is well ordered at position %i: %f > %f", j, fVariableBinsPt.GetAt(j), fVariableBinsPt.GetAt(j + 1)));

  //
  // defines an array with variable binning in multiplicity (e.g. TOF hits)
  //
  size = fVariableBinsMult.GetSize();
  delta = size > 1 ? fVariableBinsMult.At(0) - fVariableBinsMult.At(1) : 0.;
  if (size > 1 && TMath::Abs(delta) < 0.01) { //Check to see if it was set externally
    fVariableBinsMult.SetAt(1., 0);
    for (Int_t j = 1; j < fVariableBinsMult.GetSize(); j++)
      fVariableBinsMult.SetAt(fVariableBinsMult.At(j - 1) + TMath::Ceil(TMath::Sqrt(fVariableBinsMult.At(j - 1))), j);
  }
  //Checking that the binning is in order
  for (Int_t j = 0; j < size - 1; j++)
    if (fVariableBinsMult.GetAt(j) > fVariableBinsMult.GetAt(j + 1))
      AliFatal(Form("Binning fVariableBinsMult is well ordered at position %i: %f > %f", j, fVariableBinsMult.GetAt(j), fVariableBinsMult.GetAt(j + 1)));

  return;
}

//-----------------------------------------------------------
const Double_t AliAnalysisTaskTOFqaID::fBinsEta[2] = { -1.0, 1.0 };

//-----------------------------------------------------------
const Double_t AliAnalysisTaskTOFqaID::fBinsPhi[2] = { 0.0, 360.0 };

//-----------------------------------------------------------
const Double_t AliAnalysisTaskTOFqaID::fBinsT0[2] = { -700., 700. };

#endif
