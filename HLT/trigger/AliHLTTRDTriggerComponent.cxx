// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Jochen Thaeder <jochen@thaeder.de>                    *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTRDTriggerComponent.cxx
/// @author Felix Rettig, Stefan Kirsch
/// @date   2012-08-16
/// @brief

#include <cstdlib>
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "TH1I.h"
#include "TH2I.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTLogging.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliESDtrackCuts.h"
//#include "AliHLTESDTrackCuts.h"
#include "AliGeomManager.h"
#include "AliHLTComponentBenchmark.h"
#include "AliHLTTRDTriggerComponent.h"
#include "AliTRDonlineTrackingDataContainer.h"
#include "AliTRDpadPlane.h"

ClassImp(AliHLTTRDTriggerComponent)

#define LogError( ... ) { HLTError(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("ERROR", __VA_ARGS__); } }
#define LogInfo( ... ) { HLTInfo(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("INFO", __VA_ARGS__); } }
#define LogInspect( ... ) { HLTDebug(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("INSPECT", __VA_ARGS__); } }
#define LogDebug( ... ) { if (fDebugLevel >= 1) { HLTInfo(__VA_ARGS__); DbgLog("DEBUG", __VA_ARGS__); } }

AliHLTTRDTriggerComponent::AliHLTTRDTriggerComponent() :
  AliHLTTrigger(),
  fName(),
  fRefTrackSelectionEtaLimit(0.9),
  fRefTrackSelectionVertexXYLimit(20.),
  fRefTrackSelectionVertexZLimit(30.),
  fRefTrackSelectionPtThreshold(0.7),
  fMatchRatingThreshold(0.25),
  fElectronTriggerPtThresholdHSE(3.),
  fElectronTriggerPIDThresholdHSE(144),
  fElectronTriggerPtThresholdHQU(2.),
  fElectronTriggerPIDThresholdHQU(164),
  fApplyRefTrackCuts(kFALSE),
  fElectronTriggerOnL1TrgOnly(kFALSE),
  fHistoMode(1),
  fDebugLevel(0),
  fExtendedHistos(kFALSE),
  fEventRendering(kFALSE),
  fPushHistos(kFALSE),
  fWriteHistos(kFALSE),
  fEventId(fgkInvalidEventId),
  fRunNumber(-1),
  fChunkId(NULL),
  fSectorsWithData(0),
  fIsMinBiasEvent(kFALSE),
  fIsTRDElectronEvent(kFALSE),
  fESDtracksPresent(kFALSE),
  fHLTtracksPresent(kFALSE),
  fTRDGeometry(NULL),
  fEsdEvent(NULL),
  fTrackingData(NULL),
  fHLTTracks(NULL),
  fRefTrackCuts(NULL),
#ifdef __TRDHLTDEBUG
  fEventDisplay(NULL),
  fBenchmark(NULL),
#endif
  fHistArray(NULL),
  fHistMatchRating(NULL),
  fHistMatchRatingByPt(NULL),
  fHistMatchRatingByPid(NULL),
  fHistTrackPt(NULL),
  fHistTrackPtMatched(NULL),
  fHistTrackPtCorr(NULL),
  fHistTrackPid(NULL),
  fHistTrackPidMatched(NULL),
  fHistElectronCandidatePt(NULL),
  fHistElectronCandidateMatchedPt(NULL),
  fHistElectronCandidatePid(NULL),
  fHistElectronCandidateMatchedPid(NULL),
  fHistRefTrackPid(NULL),
  fHistMatchedRefTrackPid(NULL),
  fHistPIDvsTruncPID(NULL),
  fHistElectronFalsePIDvsTruncPID(NULL),
  fHistElectronConfirmedPIDvsTruncPID(NULL),
  fHistTrackMatchingCombinations(NULL),
  fHistElectronTriggerBaseMinBias(NULL),
  fHistElectronTriggerBaseTrdL1(NULL)
{
}

const char* AliHLTTRDTriggerComponent::fgkDefaultOCDBEntry = "HLT/ConfigHLT/TRDTrigger";
const char* AliHLTTRDTriggerComponent::fgkTriggerDecisionElectronHSE = "TRD-ELECTRON-HSE";
const char* AliHLTTRDTriggerComponent::fgkTriggerDecisionElectronHQU = "TRD-ELECTRON-HQU";

AliHLTTRDTriggerComponent::~AliHLTTRDTriggerComponent()
{
}

const char* AliHLTTRDTriggerComponent::GetTriggerName() const
{
  if (!fName.IsNull())
    return fName.Data();
  else
    return "TRDTriggerComponent";
}

AliHLTComponent* AliHLTTRDTriggerComponent::Spawn()
{
  return new AliHLTTRDTriggerComponent;
}

int AliHLTTRDTriggerComponent::DoInit(int argc, const char** argv)
{
  int iResult = 0;

  do {

    fChunkId = new TString("XXXXX");
#ifdef __TRDHLTDEBUG
    if (fChunkId){
      // chunk identification for debug output
      *fChunkId = gSystem->WorkingDirectory();
      fChunkId->Remove(0, fChunkId->Last('/') + 1);
      if (fChunkId->Contains("hlt_trd_comp"))
	*fChunkId = "L";
    } else {
      iResult = -ENOMEM;
      break;
    }
#endif

    fTRDGeometry = new AliTRDgeometry();
    if (!fTRDGeometry) {
      iResult = -ENOMEM;
      break;
    }

    fTrackingData = new AliTRDonlineTrackingDataContainer();
    if (!fTrackingData) {
      iResult = -ENOMEM;
      break;
    }

    fHLTTracks = new vector<AliHLTGlobalBarrelTrack>;
    if (!fHLTTracks) {
      iResult = -ENOMEM;
      break;
    }

    fRefTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "No track cuts");
    if (fRefTrackCuts){
      // fRefTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 0);
      fRefTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      fRefTrackCuts->SetEtaRange(-0.8, 0.8);
    } else {
      iResult = -ENOMEM;
      break;
    }

    fHistArray = new TObjArray(25);
    if(!fHistArray)
      return -ENOMEM;
    fHistArray->SetOwner(kTRUE);

    fHistMatchRating = new TH1I("trdtrg_match_rating", "Match rating distribution HLT tracks vs. GTU track;match rating;abundance;",
				1./0.02, 0., 1.);
    fHistArray->AddLast(fHistMatchRating);

    fHistMatchRatingByPt = new TH2I("trdtrg_match_rating_pt", "Match rating distribution HLT tracks vs. GTU track by p_{T};match rating m;p_{T}^{TRD,online} (GeV/c)",
				    1./0.02, 0., 1., 20./0.02, 0., 20.);
    fHistArray->AddLast(fHistMatchRatingByPt);

    fHistMatchRatingByPid = new TH2I("trdtrg_match_rating_pid", "Match rating distribution HLT tracks vs. GTU track by PID;match rating m;PID (a.u.)",
				     1./0.02, 0., 1., 256, 0, 256);
    fHistArray->AddLast(fHistMatchRatingByPid);

    fHistTrackPt = new TH1I("trdtrg_track_pt", "GTU track p_{T};GTU track p_{T} (GeV/c);abundance", 200, -20, 20.);
    fHistArray->AddLast(fHistTrackPt);

    fHistTrackPtMatched = new TH1I("trdtrg_matched_track_pt", "matched GTU track p_{T};GTU track p_{T} (GeV/c);abundance", 200, -20, 20.);
    fHistArray->AddLast(fHistTrackPtMatched);

    fHistTrackPtCorr = new TH2I("trdtrg_track_pt_corr", "HLT vs. GTU track p_{T};p_{T}^{HLT} (GeV/c);p_{T}^{GTU} (GeV/c)",
	       400, -40., 40., 400, -40., 40.);
    fHistArray->AddLast(fHistTrackPtCorr);

    fHistTrackPid = new TH1I("trdtrg_track_pid", "GTU track PID;GTU track PID (a.u.);abundance", 256, 0, 256);
    fHistArray->AddLast(fHistTrackPid);

    fHistTrackPidMatched = new TH1I("trdtrg_matched_track_pid", "matched GTU track PID;GTU track PID (a.u.);abundance", 256, 0, 256);
    fHistArray->AddLast(fHistTrackPidMatched);

    fHistTrackMatchingCombinations = new TH2I("trdtrg_matching_combinations", "HLT-GTU track matching combinations;set;combinations to check;",
					      10, 0, 10, 2000, 0, 10000);
    fHistArray->AddLast(fHistTrackMatchingCombinations);

    fHistElectronTriggerBaseMinBias = new TH1I("trdtrg_electron_trigger_base_mb", "min. bias base numbers for electron trigger analysis;set;abundance;",
					      10, 0, 10);
    fHistElectronTriggerBaseMinBias->GetXaxis()->SetBinLabel(1, "min. bias events");
    fHistElectronTriggerBaseMinBias->GetXaxis()->SetBinLabel(2, "TRD L1 electron triggered");
    fHistElectronTriggerBaseMinBias->GetXaxis()->SetBinLabel(3, "TRD HLT electron triggered");
    fHistArray->AddLast(fHistElectronTriggerBaseMinBias);

    fHistElectronTriggerBaseTrdL1 = new TH1I("trdtrg_electron_trigger_base_l1", "TRD L1 base numbers for electron trigger analysis;set;abundance;",
					      10, 0, 10);
    fHistElectronTriggerBaseTrdL1->GetXaxis()->SetBinLabel(1, "TRD L1 electron triggered");
    fHistElectronTriggerBaseTrdL1->GetXaxis()->SetBinLabel(2, "TRD HLT electron triggered");
    fHistArray->AddLast(fHistElectronTriggerBaseTrdL1);

    fHistElectronCandidatePt = new TH1I("trdtrg_electron_candidate_pt", "GTU electron candidate p_{T};GTU track p_{T} (GeV/c);abundance", 200, -20, 20.);
    fHistArray->AddLast(fHistElectronCandidatePt);

    fHistElectronCandidateMatchedPt = new TH1I("trdtrg_electron_candidate_matched_pt", "matching GTU electron candidate p_{T};GTU track p_{T} (GeV/c);abundance", 200, -20, 20.);
    fHistArray->AddLast(fHistElectronCandidateMatchedPt);

    fHistElectronCandidatePid = new TH1I("trdtrg_electron_candidate_pid", "GTU electron candidate PID;GTU track PID (a.u.);abundance", 256, 0, 256);
    fHistArray->AddLast(fHistElectronCandidatePid);

    fHistElectronCandidateMatchedPid = new TH1I("trdtrg_electron_candidate_matched_pid", "matching GTU electron candidate PID;GTU track PID (a.u.);abundance", 256, 0, 256);
    fHistArray->AddLast(fHistElectronCandidateMatchedPid);

  } while (0);

  if (iResult < 0){

    if (fHistArray) delete fHistArray;
    fHistArray = NULL;

    if (fTRDGeometry) delete fTRDGeometry;
    fTRDGeometry = NULL;

    if (fRefTrackCuts) delete fRefTrackCuts;
    fRefTrackCuts = NULL;

    if (fHLTTracks) delete fHLTTracks;
    fHLTTracks = NULL;

    if (fTrackingData) delete fTrackingData;
    fTrackingData = NULL;

    if (fChunkId) delete fChunkId;
    fChunkId = NULL;

  }

  // check if the -triggername argument is used, the trigger name determines the following initialization
  vector<const char*> remainingArgs;
  for (int i = 0; i < argc; ++i) {
    if (strcmp(argv[i], "-triggername") == 0) {
      if (++i < argc){
	fName = argv[i];
      } else {
	LogError("invalid parameter for argument '-triggername', string expected");
	return -EINVAL;
      }
      continue;
    }
    remainingArgs.push_back(argv[i]);
  }

  TString cdbPath;
  if (!fName.IsNull()) {
    cdbPath = "HLT/ConfigHLT/";
    cdbPath += fName;
  } else {
    cdbPath = fgkDefaultOCDBEntry;
  }

  LogInfo("cdbPath = <%s>", cdbPath.Data());
  iResult = ConfigureFromCDBObject(cdbPath);

  if (iResult >= 0 && argc > 0)
    iResult = ConfigureFromArgumentString(remainingArgs.size(), &(remainingArgs[0]));

  return iResult;
}

int AliHLTTRDTriggerComponent::DoDeinit()
{

#ifdef __TRDHLTDEBUG
  if (fEventDisplay) delete fEventDisplay;
  fEventDisplay = NULL;
#endif

  if ((fHistoMode == 1) && (fWriteHistos)){
    TFile out("trg_out/trg_hists.root", "RECREATE");
    if (!out.IsZombie()) {
      out.cd();
      UInt_t numHists = fHistArray->GetEntries();
      for (UInt_t iHist = 0; iHist < numHists; ++iHist)
	if (fHistArray->At(iHist))
	  fHistArray->At(iHist)->Write();
      out.Close();
    }
  }

  if (fHistArray) delete fHistArray;
  fHistArray = NULL;

  if (fRefTrackCuts) delete fRefTrackCuts;
  fRefTrackCuts = NULL;

  if (!fHLTTracks) delete fHLTTracks;
  fHLTTracks = NULL;

  if (fTrackingData) delete fTrackingData;
  fTrackingData = NULL;

  if (fTRDGeometry) delete fTRDGeometry;
  fTRDGeometry = NULL;

  if (fChunkId) delete fChunkId;
  fChunkId = NULL;

  return 0;
}

int AliHLTTRDTriggerComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  TString cdbPath;
  if (!cdbEntry || cdbEntry[0] == 0) {
    if (!fName.IsNull()) {
      cdbPath = "HLT/ConfigHLT/";
      cdbPath += fName;
    } else {
      cdbPath = fgkDefaultOCDBEntry;
    }
  } else
    cdbPath = cdbEntry;

  return ConfigureFromCDBObject(cdbPath);
}

int AliHLTTRDTriggerComponent::ReadPreprocessorValues(const char* /*modules*/)
{
  return 0;
}

Int_t AliHLTTRDTriggerComponent::ConfigureFromCDBObject(TString /*cdbPath*/)
{
  return 0;
}

int AliHLTTRDTriggerComponent::ScanConfigurationArgument(int argc, const char** argv)
{

  if (argc <= 0)
    return 0;

  UShort_t iArg = 0;
  TString argument(argv[iArg]);

  if (!argument.CompareTo("-match-sel-eta")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fRefTrackSelectionEtaLimit = argument.Atof();
    LogInfo("ref track selection eta limit set to %.1f.", fRefTrackSelectionEtaLimit);
    return 2;
  }

  if (!argument.CompareTo("-match-sel-vxy")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fRefTrackSelectionVertexXYLimit = argument.Atof();
    LogInfo("ref track selection vertex xy limit set to %.1f.", fRefTrackSelectionVertexXYLimit);
    return 2;
  }

  if (!argument.CompareTo("-match-sel-vz")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fRefTrackSelectionVertexZLimit = argument.Atof();
    LogInfo("ref track selection vertex z limit set to %.1f.", fRefTrackSelectionVertexZLimit);
    return 2;
  }

  if (!argument.CompareTo("-match-sel-pt")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fRefTrackSelectionPtThreshold = argument.Atof();
    LogInfo("ref track selection pt threshold set to %.1f GeV/c.", fRefTrackSelectionPtThreshold);
    return 2;
  }

  if (!argument.CompareTo("-match-rating")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fMatchRatingThreshold = argument.Atof();
    LogInfo("match rating threshold set to %.2f GeV/c.", fMatchRatingThreshold);
    return 2;
  }

  if (!argument.CompareTo("-histo-mode")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fHistoMode = argument.Atoi();
    LogInfo("histo mode set to %d", fHistoMode);
    return 2;
  }

  if (!argument.CompareTo("-trghse")){
    LogInfo("TRD HLT electron trigger HSE enabled.");
    return 1;
  }

  if (!argument.CompareTo("-hse-pt")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fElectronTriggerPtThresholdHSE = argument.Atof();
    LogInfo("pt threshold for HSE trigger set to %.1f GeV/c.", fElectronTriggerPtThresholdHSE);
    return 2;
  }

  if (!argument.CompareTo("-hse-pid")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fElectronTriggerPIDThresholdHSE = argument.Atof();
    LogInfo("PID threshold for HSE trigger set to %d.", fElectronTriggerPIDThresholdHSE);
    return 2;
  }

  if (!argument.CompareTo("-trghqu")){
    LogInfo("TRD HLT electron trigger HQU enabled.");
    return 1;
  }

  if (!argument.CompareTo("-hqu-pt")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fElectronTriggerPtThresholdHQU = argument.Atof();
    LogInfo("pt threshold for HQU trigger set to %.1f GeV/c.", fElectronTriggerPtThresholdHQU);
    return 2;
  }

  if (!argument.CompareTo("-hqu-pid")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fElectronTriggerPIDThresholdHQU = argument.Atof();
    LogInfo("PID threshold for HQU trigger set to %.1f GeV/c.", fElectronTriggerPIDThresholdHQU);
    return 2;
  }

  if (!argument.CompareTo("-l1-only")){
    LogInfo("evaluation of electron trigger only for events with TRD L1 electron trigger enabled.");
    fElectronTriggerOnL1TrgOnly = kTRUE;
    return 1;
  }

  if (!argument.CompareTo("-histo-ext")){
    LogInfo("extended histogramming enabled.");
    fExtendedHistos = kTRUE; // enable extended histogramming, for debugging/tuning only!

    if (!fHistRefTrackPid){
      fHistRefTrackPid = new TH2I("trdtrg_reftrack_pid", "TPC dE/dx by p_{T} for matching GTU tracks;p_{T} (GeV/c);dE/dx (a. u.)",
				  500, 0., 10., 500, 0., 500);
      fHistArray->AddLast(fHistRefTrackPid);
    }

    if (!fHistMatchedRefTrackPid){
      fHistMatchedRefTrackPid = new TH2I("trdtrg_matched_reftrack_pid", "TPC dE/dx by p_{T} for matching GTU tracks;p_{T} (GeV/c);dE/dx (a. u.)",
					 500, 0., 10., 500, 0., 500);
      fHistArray->AddLast(fHistMatchedRefTrackPid);
    }

    if (!fHistPIDvsTruncPID){
      fHistPIDvsTruncPID = new TH2I("trdtrg_pid_trunc_pid", "GTU track PID vs. truncated PID;PID (a.u.);truncated PID (a. u.)",
				    256, 0., 256., 256, 0., 256.);
      fHistArray->AddLast(fHistPIDvsTruncPID);
    }

    if (!fHistElectronFalsePIDvsTruncPID){
      fHistElectronFalsePIDvsTruncPID = new TH2I("trdtrg_electron_false_pid_trunc_pid", "false electron PID vs. truncated PID;PID (a.u.);truncated PID (a. u.)",
						 256, 0., 256., 256, 0., 256.);
      fHistArray->AddLast(fHistElectronFalsePIDvsTruncPID);
    }

    if (!fHistElectronConfirmedPIDvsTruncPID){
      fHistElectronConfirmedPIDvsTruncPID = new TH2I("trdtrg_electron_confirmed_pid_trunc_pid", "confirmed electron PID vs. truncated PID;PID (a.u.);truncated PID (a. u.)",
						     256, 0., 256., 256, 0., 256.);
      fHistArray->AddLast(fHistElectronConfirmedPIDvsTruncPID);
    }

    return 1;
  }

  if (!argument.CompareTo("-ref-cuts")){
    LogInfo("ref track cuts for matching enabled.");
    fApplyRefTrackCuts = kTRUE;
    return 1;
  }

  if (!argument.CompareTo("-push-histograms")){
    LogInfo("pushing of histograms enabled.");
    fPushHistos = kTRUE; // enable histogram pushing, for debugging/tuning only!
    return 1;
  }

  if (!argument.CompareTo("-write-histograms")){
    LogInfo("writing of histograms enabled.");
    fWriteHistos = kTRUE; // enable histogram writing, for debugging/tuning only!
    return 1;
  }

  if (!argument.CompareTo("-render")){
#ifdef __TRDHLTDEBUG
    LogInfo("rendering of interesting events enabled.");
    if (!fEventDisplay)
      fEventDisplay = new AliTRDtrackingEventDisplay();
      LogInfo("event rendering activated. this is for debugging only!");
    fEventRendering = kTRUE; // enable event rendering, for debugging/tuning only!
#endif
    return 1;
  }

  if (!argument.CompareTo("-debug")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fDebugLevel = argument.Atoi();
    LogInfo("debug level set to %d.", fDebugLevel);
    return 2;
  }

  if (!argument.CompareTo("-ref-pt")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    Float_t minPt, maxPt;
    fRefTrackCuts->GetPtRange(minPt, maxPt);
    maxPt = argument.Atof();
    fRefTrackCuts->SetPtRange(minPt, maxPt);
    LogInfo("ref track pt range set to %.1f .. %.1f GeV/c.", minPt, maxPt);
    return 2;
  }

  if (!argument.CompareTo("-ref-tdca")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fRefTrackCuts->SetMaxDCAToVertexXY(argument.Atof());
    LogInfo("ref track DCA transverse threshold set to %.3f", argument.Atof());
    return 2;
  }

  if (!argument.CompareTo("-ref-ldca")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fRefTrackCuts->SetMaxDCAToVertexZ(argument.Atof());
    LogInfo("ref track longitudinal DCA threshold set to %.3f", argument.Atof());
    return 2;
  }

  return 1;
}

void AliHLTTRDTriggerComponent::ScanTriggerClasses(const char* firedTriggerClasses) {

  fIsMinBiasEvent = kFALSE;

  TString trg(firedTriggerClasses);
  if (trg.Index("CINT7WU-S-NOPF-ALL") >= 0)
    fIsMinBiasEvent = kTRUE;

  if (trg.Index("CINT8WU-S-NOPF-ALL") >= 0)
    fIsMinBiasEvent = kTRUE;

}

Bool_t AliHLTTRDTriggerComponent::CheckRefTrackCuts(AliESDtrack* track){

  // cuts by cut class
  if (fApplyRefTrackCuts)
    return fRefTrackCuts->AcceptTrack(track) ? kTRUE : kFALSE;

  // simple custom cuts
  Float_t dcaToVertexXY, dcaToVertexZ;
  track->GetImpactParametersTPC(dcaToVertexXY, dcaToVertexZ);
  // LogDebug("IMPACT TPC %.4f  %.4f", dcaToVertexXY, dcaToVertexZ);

  if ((dcaToVertexXY > fRefTrackSelectionVertexXYLimit) || (dcaToVertexZ >= fRefTrackSelectionVertexZLimit))
    return kFALSE;

  if (TMath::Abs(track->Eta()) > fRefTrackSelectionEtaLimit)
    return kFALSE;

  return kTRUE;
}


void AliHLTTRDTriggerComponent::DbgLog(const char* prefix, ...){
#ifdef __TRDHLTDEBUG
  AliHLTEventID_t eventNumber = fEventId;
  printf("TRDHLTGM %s-%s-%s: [TRG] %s",
 	 (fRunNumber >= 0) ? Form("%06d", fRunNumber) : "XXXXXX",
 	 fChunkId->Data(),
 	 (eventNumber != fgkInvalidEventId) ? Form("%05llu", eventNumber) : "XXXXX",
 	 (strlen(prefix) > 0) ? Form("<%s> ", prefix) : "");
#endif
  va_list args;
  va_start(args, prefix);
  char* fmt = va_arg(args, char*);
  vprintf(fmt, args);
  printf("\n");
  va_end(args);
}

int AliHLTTRDTriggerComponent::PrepareESDData(){

  int result = 0;
  fESDtracksPresent = kFALSE;

  // check access to ESD event data
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  fEsdEvent = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));

  if (fEsdEvent) {
    fEsdEvent->GetStdContent();
    fRunNumber = fEsdEvent->GetRunNumber();
    fESDtracksPresent = kTRUE;
    unsigned int numTracks = fEsdEvent->GetNumberOfTracks();

    // process trigger classes
    ScanTriggerClasses(fEsdEvent->GetFiredTriggerClasses().Data());

    LogDebug("ESD event data: %d ESD tracks [event num: %d] [minbias: %d (fired:%s)]",
		    numTracks, fEsdEvent->GetEventNumberInFile(), fIsMinBiasEvent, fEsdEvent->GetFiredTriggerClasses().Data());

    // process ESD tracks
    if ((fDebugLevel >= 2) || fExtendedHistos){
      AliESDtrack* esdTrack;
      TString paramStr("");
      for (unsigned int iTrack = 0; iTrack < numTracks; ++iTrack) {
	esdTrack = fEsdEvent->GetTrack(iTrack);

	if (fExtendedHistos){
	  fHistRefTrackPid->Fill(esdTrack->Pt(), esdTrack->GetTPCsignal());
	}

	if (fDebugLevel >= 2){
	  paramStr = "";
	  if (esdTrack){

	    if (esdTrack->GetInnerParam())
	      paramStr += "I";

	    if (esdTrack->GetOuterParam())
	      paramStr += "O";

	    LogDebug("ESD track %4d - pt: %+.2fGeV/c  [params: S%s]", iTrack, esdTrack->GetSignedPt(), paramStr.Data());
	  } else
	    LogError("ESD data for track %d invalid", iTrack);
	}

      } // loop over ESD tracks
    }

    fTrackingData->SetGtuPtMultiplierFromMagField(fEsdEvent->GetMagneticField()); // used for sign correction

    result = 1;

  }

  return result;

}

int AliHLTTRDTriggerComponent::PrepareHLTData(){

  int iResult = 0;
  UInt_t numHLTTracks = 0;
  fHLTtracksPresent = kFALSE;

  if (!fHLTTracks){
    LogError("HLT track vector instance not available.");
    return 0;
  }

  fHLTTracks->clear();

  for (const AliHLTComponentBlockData *pBlock = GetFirstInputBlock(kAliHLTDataTypeTrack | kAliHLTDataOriginTPC);
       pBlock != NULL; pBlock = GetNextInputBlock()) {
    LogDebug("#hlttrk - data block with HLT raw tracks received");
    fHLTtracksPresent = kTRUE;

    // vector<AliHLTGlobalBarrelTrack> hltTracks;
    if ((iResult = AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, *fHLTTracks)) >= 0) {
      for (vector<AliHLTGlobalBarrelTrack>::iterator element = fHLTTracks->begin();
	   element != fHLTTracks->end(); element++) {

	numHLTTracks++;
	//## AliHLTGlobalBarrelTrack -> AliKalmanTrack -> AliExternalTrackParam

      } // loop over HLT tracks
    }
  } // loop over data blocks

  LogDebug("#hlttrk - %d HLT raw tracks found", numHLTTracks);

  return iResult;

}

int AliHLTTRDTriggerComponent::PrepareTRDData() {

  int result = 1;

  fSectorsWithData = 0;
  fTrackingData->Clear();
  for (const AliHLTComponentBlockData* datablock = GetFirstInputBlock(AliHLTTRDDefinitions::fgkOnlineDataType);
       datablock != NULL;
       datablock = GetNextInputBlock())
    {
      fSectorsWithData |= datablock->fSpecification;
      fTrackingData->Decompress(datablock->fPtr, datablock->fSize, kTRUE);
    }

  fTrackingData->SetGtuPtMultiplierFromMagField(GetBz()); // used for sign correction

  fTrackingData->PrintSummary("trigger component");

  return result;

}

int AliHLTTRDTriggerComponent::MatchTRDTracksESD(){

  if (!fEsdEvent) {
    LogError("ESD event data not available in MatchTRDTracks().");
    return 0;
  }

  if (fHistoMode == 0){
    UInt_t numHists = fHistArray->GetEntries();
    for (UInt_t iHist = 0; iHist < numHists; ++iHist)
      if (fHistArray->At(iHist))
	((TH1*) (fHistArray->At(iHist)))->Reset();
  }

  int result = 1;
  unsigned int numRefTracks = fEsdEvent->GetNumberOfTracks();
  unsigned int numGtuTracks;
  Int_t refTrackIndices[fkMaxRefTracksPerStack];
  UInt_t refTrackCount = 0;
  Bool_t isRelevant;
  Double_t distY, distZ;
  Double_t matchRating;
  Double_t bestMatchRating;
  Int_t bestMatchRefIndex;
  AliESDtrack* refTrack = NULL;
  UInt_t numComparisonsDone = 0;
  UInt_t numUnmatchedTracks = 0;
  UInt_t numMatchedTracks = 0;
  Double_t magField = fEsdEvent->GetMagneticField();

  for (UShort_t iStack = 0; iStack < fkTRDStacks; ++iStack) {
    numGtuTracks = fTrackingData->GetNumTracks(iStack);
    refTrackCount = 0;

    // preselect ESD relevant ESD tracks
    for (UInt_t iRefTrack = 0; iRefTrack < numRefTracks; ++iRefTrack) {
      refTrack = fEsdEvent->GetTrack(iRefTrack);
      isRelevant = (refTrack->Pt() >= fRefTrackSelectionPtThreshold); //## implement properly
      if (isRelevant){
	if (refTrackCount < fkMaxRefTracksPerStack)
	  refTrackIndices[refTrackCount++] = iRefTrack;
	else {
	  LogError("number of HLT tracks exceeding limit of %d. Skipping some tracks.", fkMaxRefTracksPerStack);
	  break;
	}
      }
    }

    // try to match GTU track with ref tracks
    for (UInt_t iGtuTrack = 0; iGtuTrack < numGtuTracks; ++iGtuTrack) {
      bestMatchRating = 0.;
      bestMatchRefIndex = -1;
      Double_t gpt = fTrackingData->GetTrackPt(iStack, iGtuTrack);
      UShort_t gpid = fTrackingData->GetTrackPID(iStack, iGtuTrack);
      UShort_t layerMask = fTrackingData->GetTrackLayerMask(iStack, iGtuTrack);

      Float_t trklLocalY[fkTRDLayers];
      Int_t trklBinZ[fkTRDLayers];
      for (UShort_t iLayer = 0; iLayer < fkTRDLayers; ++iLayer){
	if ((layerMask >> iLayer) & 1){
	  trklLocalY[iLayer] = fTrackingData->GetTrackTrackletLocalY(iStack, iGtuTrack, iLayer);
	  trklBinZ[iLayer] = fTrackingData->GetTrackTrackletBinZ(iStack, iGtuTrack, iLayer);
	}
      }

      for (UInt_t iRefTrack = 0; iRefTrack < refTrackCount; ++iRefTrack) {
	refTrack = fEsdEvent->GetTrack(refTrackIndices[iRefTrack]);

// 	if (!CheckRefTrackCuts(refTrack))
// 	  continue;

	numComparisonsDone++;

// 	if ((!refTrack->GetOuterParam()) && (!refTrack->GetInnerParam()))  // use tracks with TPC outer param only (all HLT tracks have it anyways)
// 	  continue;

	Int_t distRes = EstimateTrackDistance(refTrack, iStack, layerMask, trklLocalY, trklBinZ,
					      magField, &distY, &distZ);

	if (fDebugLevel >= 3){
	  printf("CHECKMATCH = %i   distY %.2f   distZ %.2f   pt: %.2f   %.2f\n",
		distRes, distY, distZ, gpt, refTrack->GetSignedPt());
	}

	if (distRes == 0){
	  matchRating = RateTrackMatch(distY, distZ, gpt, refTrack->GetSignedPt());
	} else {
	  matchRating = 0.;
	}

	if ((matchRating >= fMatchRatingThreshold) && (matchRating > bestMatchRating)) {
	  bestMatchRefIndex = refTrackIndices[iRefTrack];
	  bestMatchRating = matchRating;
	} else if (matchRating > bestMatchRating)
	  bestMatchRating = matchRating;

// 	DbgLog("", Form("#match - comparing GTU track %d in S%02d-%d with ref track %d: [gpt: %+5.2f rpt: %+5.2f] dy: %.1f dz: %.1f rating: %.1f",
// 			iGtuTrack, iStack/5, iStack%5, refTrackIndices[iRefTrack],
// 			fTrackingData->GetTrackPt(iStack, iGtuTrack), refTrack->GetSignedPt(),
// 			distY, distZ, matchRating));

      } // loop over ref tracks in stack

      fHistMatchRating->Fill(bestMatchRating);
      fHistMatchRatingByPt->Fill(bestMatchRating, TMath::Abs(gpt));
      fHistMatchRatingByPid->Fill(bestMatchRating, gpid);
      fHistTrackPt->Fill(gpt);
      fHistTrackPid->Fill(gpid);

      if (fExtendedHistos){
	fHistMatchedRefTrackPid->Fill(refTrack->Pt(), refTrack->GetTPCsignal());
      }

      if (bestMatchRefIndex >= 0){
	// GTU track has matching reference track
	refTrack = fEsdEvent->GetTrack(bestMatchRefIndex);
	Double_t rpt = refTrack->GetSignedPt();
	LogDebug("#match - match found: rating %.2f, gpt: %+5.2f, rpt: %+5.2f (%i) -> diff: %.2f%%",
			    bestMatchRating,gpt, rpt, bestMatchRefIndex, (gpt - rpt)/rpt*100.);
	fTrackingData->SetTrackAddInfo(iStack, iGtuTrack, bestMatchRefIndex);
	if (fDebugLevel >= 3){
	  LogDebug("#match-info rating: %.2f gpt: %+.2f  matchref: %d  rpt: %+.2f	gpid: %d  S%02d-%d %d	rpid: %.3f",
			  bestMatchRating, gpt, bestMatchRefIndex, rpt, gpid,
			  iStack/5, iStack%5, iGtuTrack, refTrack->GetTPCsignal());
	}
	numMatchedTracks++;

	fHistTrackPtMatched->Fill(gpt);
	fHistTrackPtCorr->Fill(rpt, gpt);
	fHistTrackPidMatched->Fill(gpid);

      } else {
	if (fDebugLevel >= 3){
	  LogDebug("#match-info rating: %.2f gpt: %+.2f  no ref matching  gpid: %d  S%02d-%d %d",
			  bestMatchRating, gpt, gpid,
			  iStack/5, iStack%5, iGtuTrack);
	}
	numUnmatchedTracks++;
      }

    } // loop over gtu tracks in stack

  } // loop over stacks

  LogInfo("#match - %d matched GTU tracks, %d unmatched",
		      numMatchedTracks, numUnmatchedTracks);

  fHistTrackMatchingCombinations->Fill(0., (Double_t)(numRefTracks * fTrackingData->GetNumTracks()));  // index 0: full combinatorics
  fHistTrackMatchingCombinations->Fill(1., numComparisonsDone);  // index 1: track matching comparisons actually done

  return result;

}

int AliHLTTRDTriggerComponent::MatchTRDTracksHLT(){

  if (!fHLTTracks){
    LogError("HLT track data available.");
    return 0;
  }

  Double_t magField = GetBz();
  unsigned int numGtuTracks;
  unsigned int numHLTTracks = 0;
  Bool_t isRelevant;
  Double_t distY, distZ;
  Double_t matchRating;
  Double_t bestMatchRating;
  Int_t bestMatchRefIndex;
  UInt_t numComparisonsDone = 0;
  UInt_t numUnmatchedTracks = 0;
  UInt_t numMatchedTracks = 0;

  Bool_t hltTrackPreSel[fkMaxRefTracks];
  UInt_t iHltTrack;

  for (UShort_t iStack = 0; iStack < fkTRDStacks; ++iStack) {
    numGtuTracks = fTrackingData->GetNumTracks(iStack);

    // preselect relevant HLT tracks
    iHltTrack = 0;
    for (vector<AliHLTGlobalBarrelTrack>::iterator element = fHLTTracks->begin();
	 element != fHLTTracks->end(); element++) {
      numHLTTracks++;
      isRelevant = (element->Pt() >= fRefTrackSelectionPtThreshold); //## implement properly

      //## use cuts here
      hltTrackPreSel[iHltTrack++] = isRelevant;
      if (iHltTrack >= fkMaxRefTracks) {
	LogError("maximum number of HLT tracks exceeded.");
	break;
      }
    } // loop over HLT tracks;

    // search for matching HLT track for each GTU track
    for (UInt_t iGtuTrack = 0; iGtuTrack < numGtuTracks; ++iGtuTrack) {
      bestMatchRating = 0.;
      bestMatchRefIndex = -1;
      Double_t gpt = fTrackingData->GetTrackPt(iStack, iGtuTrack);
      UShort_t gpid = fTrackingData->GetTrackPID(iStack, iGtuTrack);
      UShort_t layerMask = fTrackingData->GetTrackLayerMask(iStack, iGtuTrack);

      Float_t trklLocalY[fkTRDLayers];
      Int_t trklBinZ[fkTRDLayers];
      for (UShort_t iLayer = 0; iLayer < fkTRDLayers; ++iLayer){
	if ((layerMask >> iLayer) & 1){
	  trklLocalY[iLayer] = fTrackingData->GetTrackTrackletLocalY(iStack, iGtuTrack, iLayer);
	  trklBinZ[iLayer] = fTrackingData->GetTrackTrackletBinZ(iStack, iGtuTrack, iLayer);
	}
      }

      iHltTrack = 0;
      for (vector<AliHLTGlobalBarrelTrack>::iterator element = fHLTTracks->begin();
	   element != fHLTTracks->end(); element++) {
	if (!hltTrackPreSel[iHltTrack]){
	  iHltTrack++;
	  continue;
	}

	// compare GTU track and relevant HLT track
	numComparisonsDone++;

	AliHLTGlobalBarrelTrack hltPar(*element);
	Int_t distRes = EstimateTrackDistance(&hltPar, iStack, layerMask, trklLocalY, trklBinZ,
					       magField, &distY, &distZ);

	if (fDebugLevel >= 3){
	  printf("CHECKMATCH = %i   distY %.2f   distZ %.2f   pt: %.2f   %.2f\n",
		 distRes, distY, distZ, gpt, element->GetSignedPt());
	}

	if (distRes == 0){
	  matchRating = RateTrackMatch(distY, distZ, gpt, element->GetSignedPt());
	} else {
	  matchRating = 0.;
	}

	if ((matchRating >= fMatchRatingThreshold) && (matchRating > bestMatchRating)) {
	  bestMatchRefIndex = iHltTrack;
	  bestMatchRating = matchRating;
	} else if (matchRating > bestMatchRating)
	  bestMatchRating = matchRating;


	iHltTrack++;
      } // loop over HLT tracks;

      fHistMatchRating->Fill(bestMatchRating);
      fHistMatchRatingByPt->Fill(bestMatchRating, TMath::Abs(gpt));
      fHistMatchRatingByPid->Fill(bestMatchRating, gpid);
      fHistTrackPt->Fill(gpt);
      fHistTrackPid->Fill(gpid);

      if (bestMatchRefIndex >= 0){
	// GTU track has matching reference track
 	Double_t rpt = fHLTTracks->at(bestMatchRefIndex).GetSignedPt();
 	LogDebug("#match - match found: rating %.2f, gpt: %+5.2f, rpt: %+5.2f (%i) -> diff: %.2f%%",
 			    bestMatchRating,gpt, rpt, bestMatchRefIndex, (gpt - rpt)/rpt*100.);
	fTrackingData->SetTrackAddInfo(iStack, iGtuTrack, bestMatchRefIndex);

// 	if (fExtendedHistos){
// 	  fHistMatchedRefTrackPid->Fill(element->Pt(), element->GetTPCsignal());
// 	}
//
// 	if (fDebugLevel >= 3){
// 	  LogDebug("#match-info rating: %.2f gpt: %+.2f  matchref: %d  rpt: %+.2f	gpid: %d  S%02d-%d %d	rpid: %.3f",
// 			  bestMatchRating, gpt, bestMatchRefIndex, rpt, gpid,
// 			  iStack/5, iStack%5, iGtuTrack, refTrack->GetTPCsignal());
// 	}
	numMatchedTracks++;

	fHistTrackPtMatched->Fill(gpt);
	fHistTrackPtCorr->Fill(rpt, gpt);
	fHistTrackPidMatched->Fill(gpid);

      } else {
	if (fDebugLevel >= 3){
	  LogDebug("#match-info rating: %.2f gpt: %+.2f  no ref matching  gpid: %d  S%02d-%d %d",
			  bestMatchRating, gpt, gpid,
			  iStack/5, iStack%5, iGtuTrack);
	}
	numUnmatchedTracks++;
      }

    } // loop over gtu tracks
  } // loop over stacks

  LogInfo("#match - %d matched GTU tracks, %d unmatched (%d comparisons)",
	  numMatchedTracks, numUnmatchedTracks, numComparisonsDone);

  fHistTrackMatchingCombinations->Fill(0., (Double_t)(numHLTTracks * fTrackingData->GetNumTracks()));  // index 0: full combinatorics
  fHistTrackMatchingCombinations->Fill(1., numComparisonsDone);  // index 1: track matching comparisons actually done

  return kTRUE;
}

int AliHLTTRDTriggerComponent::MatchTRDTracks(){

  if (!fEsdEvent) {
    LogError("ESD event data not available in MatchTRDTracks().");
    return 0;
  }

  if (fHistoMode == 0){
    UInt_t numHists = fHistArray->GetEntries();
    for (UInt_t iHist = 0; iHist < numHists; ++iHist)
      if (fHistArray->At(iHist))
	((TH1*) (fHistArray->At(iHist)))->Reset();
  }

  int result = 1;
  unsigned int numRefTracks = fEsdEvent->GetNumberOfTracks();
  unsigned int numGtuTracks;
  Int_t refTrackIndices[fkMaxRefTracksPerStack];
  UInt_t refTrackCount = 0;
  Bool_t isRelevant;
  Double_t distY, distZ;
  Double_t matchRating;
  Double_t bestMatchRating;
  Int_t bestMatchRefIndex;
  AliESDtrack* refTrack;
  UInt_t numComparisonsDone = 0;
  UInt_t numUnmatchedTracks = 0;
  UInt_t numMatchedTracks = 0;

  for (UShort_t iStack = 0; iStack < fkTRDStacks; ++iStack) {
    numGtuTracks = fTrackingData->GetNumTracks(iStack);
    refTrackCount = 0;

    // preselect ESD relevant ESD tracks
    for (UInt_t iRefTrack = 0; iRefTrack < numRefTracks; ++iRefTrack) {
      refTrack = fEsdEvent->GetTrack(iRefTrack);
      isRelevant = (refTrack->Pt() >= fRefTrackSelectionPtThreshold);
      if (isRelevant){
	if (refTrackCount < fkMaxRefTracksPerStack)
	  refTrackIndices[refTrackCount++] = iRefTrack;
	else {
	  LogError("number of HLT tracks exceeding limit of %d. Skipping some tracks.", fkMaxRefTracksPerStack);
	  break;
	}
      }
    }

    // try to match GTU track with ref tracks
    for (UInt_t iGtuTrack = 0; iGtuTrack < numGtuTracks; ++iGtuTrack) {
      bestMatchRating = 0.;
      bestMatchRefIndex = -1;
      Double_t gpt = fTrackingData->GetTrackPt(iStack, iGtuTrack);
      UShort_t gpid = fTrackingData->GetTrackPID(iStack, iGtuTrack);
      UShort_t layerMask = fTrackingData->GetTrackLayerMask(iStack, iGtuTrack);

      for (UInt_t iRefTrack = 0; iRefTrack < refTrackCount; ++iRefTrack) {
	refTrack = fEsdEvent->GetTrack(refTrackIndices[iRefTrack]);

	if (!CheckRefTrackCuts(refTrack))
	  continue;

	numComparisonsDone++;

// 	if ((!refTrack->GetOuterParam()) && (!refTrack->GetInnerParam()))  // use tracks with TPC outer param only (all HLT tracks have it anyways)
// 	  continue;

	Float_t trklLocalY[fkTRDLayers];
	Int_t trklBinZ[fkTRDLayers];
	for (UShort_t iLayer = 0; iLayer < fkTRDLayers; ++iLayer){
	  if ((layerMask >> iLayer) & 1){
	    trklLocalY[iLayer] = fTrackingData->GetTrackTrackletLocalY(iStack, iGtuTrack, iLayer);
	    trklBinZ[iLayer] = fTrackingData->GetTrackTrackletBinZ(iStack, iGtuTrack, iLayer);
	  }
	}

	Int_t distRes = EstimateTrackDistance(refTrack, iStack, layerMask, trklLocalY, trklBinZ,
					      fEsdEvent->GetMagneticField(), &distY, &distZ);

	if (distRes == 0){
	  matchRating = RateTrackMatch(distY, distZ, gpt, refTrack->GetSignedPt());
	} else {
	  matchRating = 0.;
	}

	if ((matchRating >= fMatchRatingThreshold) && (matchRating > bestMatchRating)) {
	  bestMatchRefIndex = refTrackIndices[iRefTrack];
	  bestMatchRating = matchRating;
	} else if (matchRating > bestMatchRating)
	  bestMatchRating = matchRating;

// 	DbgLog("", Form("#match - comparing GTU track %d in S%02d-%d with ref track %d: [gpt: %+5.2f rpt: %+5.2f] dy: %.1f dz: %.1f rating: %.1f",
// 			iGtuTrack, iStack/5, iStack%5, refTrackIndices[iRefTrack],
// 			fTrackingData->GetTrackPt(iStack, iGtuTrack), refTrack->GetSignedPt(),
// 			distY, distZ, matchRating));

      } // loop over ref tracks in stack

      fHistMatchRating->Fill(bestMatchRating);
      fHistMatchRatingByPt->Fill(bestMatchRating, TMath::Abs(gpt));
      fHistMatchRatingByPid->Fill(bestMatchRating, gpid);
      fHistTrackPt->Fill(gpt);
      fHistTrackPid->Fill(gpid);

      if (bestMatchRefIndex >= 0){
	// GTU track has matching reference track
	refTrack = fEsdEvent->GetTrack(bestMatchRefIndex);
	Double_t rpt = refTrack->GetSignedPt();
	LogDebug("#match - match found: rating %.2f, gpt: %+5.2f, rpt: %+5.2f (%i) -> diff: %.2f%%",
			    bestMatchRating,gpt, rpt, bestMatchRefIndex, (gpt - rpt)/rpt*100.);
	fTrackingData->SetTrackAddInfo(iStack, iGtuTrack, bestMatchRefIndex);

	if (fExtendedHistos){
	  fHistMatchedRefTrackPid->Fill(refTrack->Pt(), refTrack->GetTPCsignal());
	}

	if (fDebugLevel >= 3){
	  LogDebug("#match-info rating: %.2f gpt: %+.2f  matchref: %d  rpt: %+.2f	gpid: %d  S%02d-%d %d	rpid: %.3f",
			  bestMatchRating, gpt, bestMatchRefIndex, rpt, gpid,
			  iStack/5, iStack%5, iGtuTrack, refTrack->GetTPCsignal());
	}
	numMatchedTracks++;

	fHistTrackPtMatched->Fill(gpt);
	fHistTrackPtCorr->Fill(rpt, gpt);
	fHistTrackPidMatched->Fill(gpid);

      } else {
	if (fDebugLevel >= 3){
	  LogDebug("#match-info rating: %.2f gpt: %+.2f  no ref matching  gpid: %d  S%02d-%d %d",
			  bestMatchRating, gpt, gpid,
			  iStack/5, iStack%5, iGtuTrack);
	}
	numUnmatchedTracks++;
      }

    } // loop over gtu tracks in stack

  } // loop over stacks

  LogInfo("#match - %d matched GTU tracks, %d unmatched",
		      numMatchedTracks, numUnmatchedTracks);

  fHistTrackMatchingCombinations->Fill(0., (Double_t)(numRefTracks * fTrackingData->GetNumTracks()));  // index 0: full combinatorics
  fHistTrackMatchingCombinations->Fill(1., numComparisonsDone);  // index 1: track matching comparisons actually done

  return result;
}

void AliHLTTRDTriggerComponent::DumpTrackingData(){

  if (fTrackingData->GetNumTracklets() + fTrackingData->GetNumTracks() == 0)
    return;

  TString trklStr("");
  TString matchStr("");
  UShort_t layerMask;

//  for (UShort_t iSector = 0; iSector < 18; ++iSector){
//    if (fTrackingData->GetSectorTrgWord(iSector) != ((UInt_t)0x2345352 ^ ((UInt_t)iSector + 34)))
//      LogError("invalid sector trigger word in sector %02d: trg word is 0x%08x, should be 0x%08x",
// 	       iSector, fTrackingData->GetSectorTrgWord(iSector), ((UInt_t)0x2345352 ^ ((UInt_t)iSector + 34)));
//    for (UShort_t iStack = 0; iStack < 5; ++iStack){
//      ULong64_t dwl = ((ULong64_t)0xaffe << 16) | (ULong64_t)iSector << 8 | iStack;
//      ULong64_t dw = (((ULong64_t)0xbead << 16) | (ULong64_t)iSector << 8 | iStack) | ((ULong64_t)dwl << 32);
//      if (fTrackingData->GetStackTrgWord(iSector, iStack) != dw)
// 	LogError("stack %02d-%d trg word is 0x%016llx, should be 0x%016llx",
// 		 iSector, iStack, fTrackingData->GetStackTrgWord(iSector, iStack), dw);
//    }
//  }

  for (UShort_t iStack = 0; iStack < fkTRDStacks; ++iStack){
    for (Int_t iTrk = 0; iTrk < fTrackingData->GetNumTracks(iStack); ++iTrk){

      layerMask = fTrackingData->GetTrackLayerMask(iStack, iTrk);

      trklStr = Form("trkl: ");
      for (Short_t iLayer = 5; iLayer >= 0; --iLayer){
	if ((layerMask >> iLayer) & 1)
	  trklStr += Form("0x%08x (%+8.3f)  ",
			  fTrackingData->GetTrackTrackletWord(iStack, iTrk, iLayer),
			  fTrackingData->GetTrackTrackletLocalY(iStack, iTrk, iLayer));
	else
	  trklStr += "---------------------  ";
      } // loop over layers
      trklStr.Remove(trklStr.Length() - 2, 2);

      if (fTrackingData->GetTrackAddInfo(iStack, iTrk) >= 0){
	Double_t rpt = (fESDtracksPresent) ?  fEsdEvent->GetTrack(fTrackingData->GetTrackAddInfo(iStack, iTrk))->GetSignedPt() :
	  fHLTTracks->at(fTrackingData->GetTrackAddInfo(iStack, iTrk)).GetSignedPt();
	matchStr = Form("mpt: %+7.2f", rpt);
      } else
	matchStr = "unmatched";

      if (fDebugLevel >= 3){

	printf("###DOTDA EV%04llu  GTU TRACK - S%02d-%d	 pt: %+7.2f  pid: %3d  lm: 0x%02x %s  %s\n",
	       fEventId,
	       iStack/5, iStack%5, fTrackingData->GetTrackPt(iStack, iTrk), fTrackingData->GetTrackPID(iStack, iTrk),
	       layerMask, trklStr.Data(),
	       matchStr.Data());

	printf("###DOTDB EV%04llu  GTU TRACK - S%02d-%d	 pt: %+7.2f  pid: %3d  lm: 0x%02x %s\n",
	       fEventId,
	       iStack/5, iStack%5, fTrackingData->GetTrackPt(iStack, iTrk), fTrackingData->GetTrackPID(iStack, iTrk),
	       layerMask, trklStr.Data());
      }


      // paranoia checks
      for (Short_t iLayer = 5; iLayer >= 0; --iLayer){
	if (((layerMask >> iLayer) & 1) && (fTrackingData->GetTrackTrackletWord(iStack, iTrk, iLayer) == 0x10001000))
	  LogError("invalid layer mask / tracklet value combination A");

	if ((((layerMask >> iLayer) & 1) == 0) && (fTrackingData->GetTrackTrackletWord(iStack, iTrk, iLayer) != 0x10001000))
	  LogError("invalid layer mask / tracklet value combination B");
      }

    } // loop over tracks in stack
  } // loop over stacks

}

void AliHLTTRDTriggerComponent::AssignTrackInfo(TString* infoStr, const UInt_t stack, const UInt_t trackIndex, const char* flagStr) {

  TString flags(flagStr);
  Bool_t appendRefData = kFALSE;

  if (fTrackingData->GetTrackAddInfo(stack, trackIndex) >= 0){
    appendRefData = kTRUE;
    if (flags.First('M') < 0)
      flags += "M";
  }

  *infoStr = Form("EXCH-TRK-INFO DS<%s> CH<%s> EV%05llu SEC%02d-%d-%d  gpt: %+6.1f  gpid: %3d  flags: %s",
		  "TRDHLTTRG",
		  fChunkId->Data(), fEventId,
		  stack/5, stack%5, trackIndex,
		  fTrackingData->GetTrackPt(stack, trackIndex),
		  fTrackingData->GetTrackPID(stack, trackIndex), flags.Data());

  if (appendRefData){
    Double_t rpt = (fESDtracksPresent) ? fEsdEvent->GetTrack(fTrackingData->GetTrackAddInfo(stack, trackIndex))->GetSignedPt() :
      fHLTTracks->at(fTrackingData->GetTrackAddInfo(stack, trackIndex)).GetSignedPt();
    *infoStr += Form("  rpt: %+6.1f", rpt);
  }

}

#ifdef __TRDHLTDEBUG
void AliHLTTRDTriggerComponent::RenderEvent(const Bool_t showGtuTracks, const Bool_t showTracklets, const Bool_t showRefTracks) {

  if ((!fEventDisplay) || (!fEsdEvent))
    return;

  LogDebug("rendering event");

  const Float_t refTrackPtDisplayThreshold = 0.7;
  const Float_t trackPtEmphasizeThreshold = 1.8;

  fEventDisplay->Reset();
  fEventDisplay->SetMagField(fEsdEvent->GetMagneticField());
  for (UInt_t iDet = 0; iDet < fkTRDStacks*fkTRDLayers; ++iDet)
    fEventDisplay->SetChamberState(iDet/fkTRDLayers, iDet%6, ((fSectorsWithData >> (iDet/30)) & 1) );

  if (showTracklets){
    for (UShort_t iDet = 0; iDet < fkTRDStacks*fkTRDLayers; ++iDet){
      UInt_t numTrkl = fTrackingData->GetNumTracklets(iDet);
      for (UInt_t iTrkl = 0; iTrkl < numTrkl; ++iTrkl){
	fEventDisplay->AddTracklet(fTrackingData->GetTracklet(iDet, iTrkl));
      }
    }
  }

  if (showGtuTracks){
    for (UShort_t iStack = 0; iStack < fkTRDStacks; ++iStack){
      for (Int_t iTrk = 0; iTrk < fTrackingData->GetNumTracks(iStack); ++iTrk){
	AliESDTrdTrack* track = fTrackingData->GetTrack(iStack, iTrk);
	Int_t trkIndex = fEventDisplay->AddTrack(track, (TMath::Abs(track->Pt()) >= trackPtEmphasizeThreshold) ? (kMagenta + 2) : kMagenta, 1, 1);
	if (fTrackingData->GetTrackAddInfo(iStack, iTrk) >= 0)
	  fEventDisplay->SetTrackTrackletStyle(trkIndex, kRed, 12);
	else
	  fEventDisplay->SetTrackTrackletStyle(trkIndex, kViolet - 5, 12);
      }
    }
  }

  unsigned int numRefTracks = fEsdEvent->GetNumberOfTracks();
  unsigned int numRefTracksRendered = 0;
  if (showRefTracks){
    // check for some marks for rendering
    UShort_t marks[40000];
    memset(marks, 0, sizeof(UShort_t)*40000);
    for (UShort_t iStack = 0; iStack < fkTRDStacks; ++iStack){
      for (Int_t iTrk = 0; iTrk < fTrackingData->GetNumTracks(iStack); ++iTrk){
	if (fTrackingData->GetTrackAddInfo(iStack, iTrk) >= 0){
	  marks[fTrackingData->GetTrackAddInfo(iStack, iTrk)] |= 1;
	  if (
	      (fTrackingData->GetTrackPID(iStack, iTrk) >= fElectronTriggerPIDThresholdHQU) &&
	      (TMath::Abs(fTrackingData->GetTrackPt(iStack, iTrk)) >= fElectronTriggerPtThresholdHQU))
	    marks[fTrackingData->GetTrackAddInfo(iStack, iTrk)] |= 2;
	}
      }
    }

    // add to rendering
    for (unsigned int iTrack = 0; iTrack < numRefTracks; ++iTrack){
      AliESDtrack* track = fEsdEvent->GetTrack(iTrack);
      if ((track->Pt() >= refTrackPtDisplayThreshold) || (marks[iTrack] != 0)){
	Color_t color = (track->Pt() >= trackPtEmphasizeThreshold) ? (kGray + 2) : kGray;
	UShort_t width = 1;
	UShort_t style = 1;
	if (marks[iTrack] & 0x2){
	  color = kRed + 1;
	  width = 6;
	}
	if (!CheckRefTrackCuts(track))
	  style = 2;

	fEventDisplay->AddTrack(track, color, width, style);
	numRefTracksRendered++;
      }
    }
  }

  if (!fEventDisplay->IsEmpty()){
    fEventDisplay->SetTitle("");
    fEventDisplay->SetSetupText("HLT", Form("%.1f#scale[0.5]{ }/#scale[0.5]{ }%.1f",
					  refTrackPtDisplayThreshold, trackPtEmphasizeThreshold));
    fEventDisplay->SetBottomText(Form("%05i-%s-%05llu",
				      fRunNumber, fChunkId->Data(), fEventId));
    fEventDisplay->SetBottomTextRight(Form("%d (%d) HLT tracks, %d tracklets, %d GTU tracks",
					   numRefTracks,
					   numRefTracksRendered,
					   fTrackingData->GetNumTracklets(),
					   fTrackingData->GetNumTracks()));
    fEventDisplay->SetLook(AliTRDtrackingEventDisplay::dmMediumLight);
    fEventDisplay->SetDisplayMode(AliTRDtrackingEventDisplay::dmFullXY);
    fEventDisplay->SaveToFile(Form("display/event-%s-%05llu.eps",
				   fChunkId->Data(), fEventId));
  }

}
#endif

Bool_t AliHLTTRDTriggerComponent::TRDElectronTrigger(const char *ident, const Double_t minPt, const UShort_t minPID){

  LogDebug("Electron trigger processing (%s: pt>=%.1f, pid>=%d)...", ident, minPt, minPID);

  UInt_t numTracks;
  Bool_t highPtElectronSeenGTU = kFALSE;
  Bool_t highPtElectronSeen = kFALSE;
  UInt_t truncMeanPID = 0;
  TString trackExchangeInfo("");
  TString flags("");

  UInt_t trdSectorTrgContribs = 0;
  for (UShort_t iSector = 0; iSector < fkTRDSectors; ++iSector)
    trdSectorTrgContribs |= fTrackingData->GetSectorTrgContribs(iSector);
  if ((trdSectorTrgContribs >> 5) & 0x3)
    fIsTRDElectronEvent = kTRUE;
  else
    fIsTRDElectronEvent = kFALSE;

  if (fIsMinBiasEvent)
    fHistElectronTriggerBaseMinBias->Fill(0.);

  if (fIsTRDElectronEvent)
    fHistElectronTriggerBaseTrdL1->Fill(0.);

  if ((fElectronTriggerOnL1TrgOnly) && (!fIsTRDElectronEvent)){
    // evaluate trigger for events with TRD L1 electron trigger fired
    DbgLog("skipping %s electron trigger evaluation for event, because no TRD trigger flag set (ctbs: 0x%02x)",
	   ident, trdSectorTrgContribs);
    return 0;
  }

  for (UShort_t iStack = 0; iStack < fkTRDStacks; ++iStack){
    numTracks = fTrackingData->GetNumTracks(iStack);
    for (UInt_t iTrk = 0; iTrk < numTracks; ++iTrk){
      Double_t gpt = fTrackingData->GetTrackPt(iStack, iTrk);
      Bool_t trdElectronCandidate = kFALSE;
      Bool_t hltElectronCandidate = kFALSE;

      // re-evaluate GTU only decision for comparison
      if (
	  (TMath::Abs(gpt) >= minPt) &&
	  (fTrackingData->GetTrackPID(iStack, iTrk) >= minPID)
	  ) {
	trdElectronCandidate = kTRUE;
	highPtElectronSeenGTU = kTRUE;
	fHistElectronCandidatePt->Fill(gpt);
	fHistElectronCandidatePid->Fill(fTrackingData->GetTrackPID(iStack, iTrk));

	if (fExtendedHistos){

	  // idea to be checked
	  truncMeanPID = 0;
	  Short_t minLyr = -1;
	  UShort_t minValue = 255;
	  Short_t maxLyr = -1;
	  UShort_t maxValue = 0;
	  UShort_t lyrCount = 0;
	  UInt_t layerMask = fTrackingData->GetTrackLayerMask(iStack, iTrk);

	  // scan for min & max values
	  for (UShort_t iLayer = 0; iLayer < fkTRDLayers; ++iLayer){
	    if ((layerMask >> iLayer) & 1){
	      if (fTrackingData->GetTrackTrackletPID(iStack, iTrk, iLayer) < minValue){
		minValue = fTrackingData->GetTrackTrackletPID(iStack, iTrk, iLayer);
		minLyr = iLayer;
	      }
	      if (fTrackingData->GetTrackTrackletPID(iStack, iTrk, iLayer) > maxValue){
		maxValue = fTrackingData->GetTrackTrackletPID(iStack, iTrk, iLayer);
		maxLyr = iLayer;
	      }
	    }
	  } // loop over layers

	  // calculate trunc mean
	  for (UShort_t iLayer = 0; iLayer < fkTRDLayers; ++iLayer){
	    if (((layerMask >> iLayer) & 1) && (iLayer != minLyr) && (iLayer != maxLyr)){
	      truncMeanPID += fTrackingData->GetTrackTrackletPID(iStack, iTrk, iLayer);
	      lyrCount++;
	    }
	  } // loop over layers
	  truncMeanPID = TMath::Nint((Double_t)(truncMeanPID)/(Double_t)(lyrCount));

	  fHistPIDvsTruncPID->Fill(fTrackingData->GetTrackPID(iStack, iTrk), truncMeanPID);
	  if (fTrackingData->GetTrackAddInfo(iStack, iTrk) < 0)
	    fHistElectronFalsePIDvsTruncPID->Fill(fTrackingData->GetTrackPID(iStack, iTrk), truncMeanPID);

	}

	LogInspect("#hlt-trd-trg - GTU flagged %s %s high-pt electron seen in S%02d-%d: pt: %+6.1f  pid: %d  [id: S%02d-%d-%d] [trunc-pid: %d]",
		   (fTrackingData->GetTrackAddInfo(iStack, iTrk) < 0) ? "unmatched" : "matched", ident,
		   iStack/5, iStack%5,
		   gpt, fTrackingData->GetTrackPID(iStack, iTrk),
		   iStack/5, iStack%5, iTrk,
		   truncMeanPID);
      }

      // evaluate HLT-level trigger decision using match information
      if (
	  (fTrackingData->GetTrackAddInfo(iStack, iTrk) >= 0) &&
	  (TMath::Abs(gpt) >= minPt) &&
	  (fTrackingData->GetTrackPID(iStack, iTrk) >= minPID)
	  ) {
	hltElectronCandidate = kTRUE;
	highPtElectronSeen = kTRUE;

	fHistElectronCandidateMatchedPt->Fill(gpt);
	fHistElectronCandidateMatchedPid->Fill(fTrackingData->GetTrackPID(iStack, iTrk));

	if (fExtendedHistos){
	  fHistElectronConfirmedPIDvsTruncPID->Fill(fTrackingData->GetTrackPID(iStack, iTrk), truncMeanPID);
	}

	Double_t rpt = (fESDtracksPresent) ? fEsdEvent->GetTrack(fTrackingData->GetTrackAddInfo(iStack, iTrk))->GetSignedPt() :
	  fHLTTracks->at(fTrackingData->GetTrackAddInfo(iStack, iTrk)).GetSignedPt();
	LogInspect("#hlt-trd-trg - HLT matched %s high-pt electron seen in S%02d-%d: gpt: %+6.1f  gpid: %d   rpt: %+6.1f  [id: S%02d-%d-%d]",
		   ident, iStack/5, iStack%5, fTrackingData->GetTrackPt(iStack, iTrk), fTrackingData->GetTrackPID(iStack, iTrk),
		   rpt, iStack/5, iStack%5, iTrk);
      }

      // log output for subsequent offline analysis
      if ((fDebugLevel >= 3) && (trdElectronCandidate || hltElectronCandidate)){
	trackExchangeInfo = "";
	flags = "";
	if (trdElectronCandidate)
	  flags += "G";
	if (hltElectronCandidate)
	  flags += "H";
	AssignTrackInfo(&trackExchangeInfo, iStack, iTrk, flags.Data());
	LogDebug("%s\n", trackExchangeInfo.Data());
      }

    } // loop over tracks in stack
  } // loop over stacks

  if (highPtElectronSeenGTU || highPtElectronSeen){
    LogInspect("#hlt-trd-trg - event triggered by %s electron trigger (TRD L1: %d, TRD HLT: %d)",
	       ident, highPtElectronSeenGTU, highPtElectronSeen);
  }

  if (highPtElectronSeenGTU){
    if (fIsMinBiasEvent)
      fHistElectronTriggerBaseMinBias->Fill(1.);
  }

  if (highPtElectronSeen){
    if (fIsMinBiasEvent)
      fHistElectronTriggerBaseMinBias->Fill(2.);
    if (fIsTRDElectronEvent)
      fHistElectronTriggerBaseTrdL1->Fill(1.);
  }

  return (highPtElectronSeen) ? kTRUE : kFALSE;
}


int AliHLTTRDTriggerComponent::DoTrigger()
{

  fEsdEvent = NULL;
  fRunNumber = -1;
  int iResult = 0;
  UShort_t firedTriggers = 0;

  const AliHLTComponentEventData* hltEventData = GetEventData();
  fEventId = hltEventData->fEventID;
  LogDebug("### START DoTrigger [event id: %llu, %d blocks, size: %d]",
		      fEventId, hltEventData->fBlockCnt, hltEventData->fStructSize);

//  LogDebug("### START DoTrigger [event id: %llu, %d blocks, size: %d]",
// 		      fEventId, hltEventData->fBlockCnt, hltEventData->fStructSize);

  if (!IsDataEvent()) {  // process data events only
    IgnoreEvent();
    LogDebug("### END   DoTrigger [event id: %llu, %d blocks, size: %d] (skipped: no data event)",
	     fEventId, hltEventData->fBlockCnt, hltEventData->fStructSize);
    return iResult;
  }

  fTrackingData->SetLogPrefix(Form("TRDHLTGM XXXXXX-%05llu: [TRG] {TrkDat} ", fEventId));

  do {

    // access to TRD specific data from AliHLTTRDPreprocessorComponent
    if (!PrepareTRDData()){
      LogError("access to TRD data failed. Skipping event...");
      break;
    }

    if (fTrackingData->GetNumTracks() + fTrackingData->GetNumTracklets() == 0) {
      LogDebug("no trigger-relevant TRD information, skipping further event processing");
      break;
    }

    // access to ESD data
    if (!PrepareESDData()){
      LogInfo("access to ESD event data failed.");
    }

    // access to alternative HLT data
    if (!fESDtracksPresent){
      if (!PrepareHLTData()){
	LogError("access to HLT event data failed.");
      }
    }

    // match TRD and HLT tracks
    if (fESDtracksPresent){
      if (!MatchTRDTracksESD()){
	LogError("matching TRD tracks to ESD tracks failed. Skipping event...");
	break;
      }
    } else if (fHLTtracksPresent){
      if (!MatchTRDTracksHLT()){
	LogError("matching TRD tracks to HLT tracks failed. Skipping event...");
	break;
      }
    } else {
      LogError("No HLT track information available. Skipping event...");
      break;
    }

//    if (!MatchTRDTracks()){
//      LogError("matching TRD tracks to TPC tracks failed. Skipping event...");
//      break;
//    }

    if (fDebugLevel >= 1)
      DumpTrackingData();

    // evaluate electron trigger conditions
    if (TRDElectronTrigger("HSE", fElectronTriggerPtThresholdHSE, fElectronTriggerPIDThresholdHSE))
      firedTriggers |= fkElectronTriggerHSE;

    if (TRDElectronTrigger("HQU", fElectronTriggerPtThresholdHQU, fElectronTriggerPIDThresholdHQU))
      firedTriggers |= fkElectronTriggerHQU;

    break;

  } while (1);


  // trigger decision
  TString description("");
  if (firedTriggers & fkElectronTriggerHSE){
    if (description.Length() > 0)
      description += " ";
    description += fgkTriggerDecisionElectronHSE;
  }

  if (firedTriggers & fkElectronTriggerHQU){
    if (description.Length() > 0)
      description += " ";
    description += fgkTriggerDecisionElectronHQU;
  }

  SetDescription(description.Data());
  AliHLTTriggerDecision decision((firedTriggers) ? kTRUE : kFALSE,
				 GetTriggerName(),
				 GetReadoutList(),
				 GetDescription()
				 );
  TriggerEvent(&decision, kAliHLTDataTypeTObject | kAliHLTDataOriginOut);

  if (firedTriggers){
    LogInspect("TRD HLT trigger fired for event: description: >%s<, flags: 0x%04x",
	       description.Data(), firedTriggers);
#ifdef __TRDHLTDEBUG
    if (fEventRendering)
      RenderEvent();
#endif
  } else {
    LogDebug("TRD HLT trigger did not fire for event");
  }

  if (fPushHistos){
    PushBack(fHistArray, (kAliHLTDataTypeTObjArray | kAliHLTDataOriginTRD), 0x3fffff);
  }

  LogDebug("### END   DoTrigger [event id: %llu, %d blocks, size: %d]",
		      fEventId, hltEventData->fBlockCnt, hltEventData->fStructSize);

  return iResult;
}

Bool_t AliHLTTRDTriggerComponent::TrackPlaneIntersect(AliExternalTrackParam *trk, Double_t pnt[3], Double_t norm[3], Double_t mag){

  UInt_t its = 0;
  Double_t r = 290.;
  Double_t dist = 99999, dist_prev = 99999;
  Double_t x[3] = {0., 0., 0.};
  Bool_t ret = kTRUE;

  dist = (x[0] - pnt[0]) * norm[0] + (x[1] - pnt[1]) *norm[1] + (x[2] - pnt[2]) * norm[2];

  while(TMath::Abs(dist) > 0.1) {

    trk->GetXYZAt(r, mag, x);

    if ((x[0] * x[0] + x[1] * x[1]) < 100.) {  // extrapolation to radius failed
      ret = kFALSE;
      break;
    }

    //distance between current track position and plane
    dist_prev = TMath::Abs(dist);
    dist = (x[0] - pnt[0]) * norm[0] + (x[1] - pnt[1]) * norm[1];
    r -= dist;
    its++;
    if(TMath::Abs(dist) >= dist_prev ||
       (r > 380.) || (r < 100.)){
      ret = kFALSE;
      break;
    }
  }

  for (Int_t i=0; i<3; i++){
    if(ret)
      pnt[i] = x[i];
    else
      pnt[i] = 0.;
  }

  return kTRUE;
}


Double_t AliHLTTRDTriggerComponent::RateTrackMatch(Double_t distY, Double_t distZ, Double_t rpt, Double_t gpt){

  // maximum limits for spatial distance
  if ((distY > 5.) || (distZ > 20.))
    return 0.;

  // same pt sign required
  if ((rpt * gpt) < 0.)
    return 0.;

  Double_t rating_distY = -0.1 * distY + 1.;
  Double_t rating_distZ = -0.025 * distZ + 1.;
  Double_t rating_ptDiff = 1. - TMath::Abs((TMath::Abs(rpt) > 0.000001) ? ((gpt-rpt)/rpt) : 0.);

  if (rating_ptDiff <  0.)
    rating_ptDiff = 0.2;

  Double_t total = rating_distY * rating_distZ * rating_ptDiff;

//  DbgLog("", Form("#matching: rating:   dy: %.3f   dz: %.3f   dpt: %.3f     -> total: %.3f",
// 		    rating_distY, rating_distZ, rating_ptDiff, total));

  if (total > 1.)
    LogError("track match rating exceeds limit of 1.0: %.3f", total);

  return total;
}

Int_t AliHLTTRDTriggerComponent::EstimateTrackDistance(AliESDtrack *esd_track,
							const UShort_t stack,
							const UShort_t layerMask,
							const Float_t trklLocalY[6], const Int_t trklBinZ[6],
							Double_t mag, Double_t *ydist, Double_t *zdist){
  if (!esd_track)
    return -3;

  AliExternalTrackParam* refParam = NULL;
  if (esd_track->GetOuterParam())
    refParam = new AliExternalTrackParam(*(esd_track->GetOuterParam()));
  if (!refParam)
    refParam = new AliExternalTrackParam(*(esd_track));

  Int_t res = EstimateTrackDistance(refParam, stack, layerMask, trklLocalY, trklBinZ, mag, ydist, zdist);

  if (refParam)
    delete refParam;

  return res;

}

Int_t AliHLTTRDTriggerComponent::EstimateTrackDistance(AliExternalTrackParam *refParam,
							const UShort_t stack,
							const UShort_t layerMask,
							const Float_t trklLocalY[6], const Int_t trklBinZ[6],
							Double_t mag, Double_t *ydist, Double_t *zdist){

  Float_t diff_y = 0;
  Float_t diff_z = 0;
  Int_t nLayers = 0;
  Double_t xtrkl[3];
  Double_t ptrkl[3];
  Double_t ptrkl2[3];
  UInt_t trklDet;
  UShort_t trklLayer;
  // UInt_t stack_gtu;
  UShort_t stackInSector;

  AliTRDpadPlane* padPlane;

  for (UShort_t iLayer = 0; iLayer < 6; iLayer++){
    if ((layerMask >> iLayer) & 1){
      trklDet = stack*6 + iLayer;
      trklLayer = iLayer;
      // stack_gtu = stack;
      stackInSector = stack % 5;

      // local coordinates of the outer end point of the tracklet
      xtrkl[0] = AliTRDgeometry::AnodePos();
      xtrkl[1] = trklLocalY[iLayer];

      padPlane = fTRDGeometry->GetPadPlane(trklLayer, stackInSector);
      if(stackInSector == 2){ // corrected version by Felix Muecke
	xtrkl[2] = padPlane->GetRowPos(trklBinZ[iLayer]) - (padPlane->GetRowSize(trklBinZ[iLayer]))/2. - padPlane->GetRowPos(6);
      } else {
	xtrkl[2] = padPlane->GetRowPos(trklBinZ[iLayer]) - (padPlane->GetRowSize(trklBinZ[iLayer]))/2. - padPlane->GetRowPos(8);
      }

      // transform to global coordinates
      TGeoHMatrix *matrix = fTRDGeometry->GetClusterMatrix(trklDet);
      if (!matrix){
	LogError("invalid TRD cluster matrix in EstimateTrackDistance for detector %i", trklDet);
	return -5;
      }
      matrix->LocalToMaster(xtrkl, ptrkl);
      fTRDGeometry->RotateBack((stack/5) * 30, ptrkl, ptrkl2);  // ptrkl2 now contains the global position of the outer end point of the tracklet

      // calculate parameterization of plane representing the tracklets layer
      Double_t layer_zero_local[3] = {0., 0.,	0.};
      Double_t layer_zero_global[3], layer_zero_global2[3];

      matrix->LocalToMaster(layer_zero_local, layer_zero_global);
      fTRDGeometry->RotateBack(trklDet, layer_zero_global, layer_zero_global2); // layer_zero_global2 points to chamber origin in global coords

      Double_t layer_ref_local[3] = {AliTRDgeometry::AnodePos(), 0.,  0.};
      Double_t layer_ref_global[3], layer_ref_global2[3];

      matrix->LocalToMaster(layer_ref_local, layer_ref_global);
      fTRDGeometry->RotateBack(trklDet, layer_ref_global, layer_ref_global2); // layer_ref_global2 points to center anode pos within plane in global coords

      Double_t n0[3] = {layer_ref_global2[0]-layer_zero_global2[0],
			layer_ref_global2[1]-layer_zero_global2[1],
			layer_ref_global2[2]-layer_zero_global2[2]};

      Double_t n_len = TMath::Sqrt(n0[0]*n0[0] + n0[1]*n0[1] + n0[2]*n0[2]);
      if (n_len == 0.){ // This should never happen
	//printf("<ERROR> divison by zero in estimate_track_distance!");
	n_len = 1.;
      }
      Double_t n[3] = {n0[0]/n_len, n0[1]/n_len, n0[2]/n_len}; // normal vector of plane

      Bool_t isects = TrackPlaneIntersect(refParam, layer_ref_global2, n, mag); // find intersection point between track and TRD layer

      if (isects == kFALSE){ // extrapolation fails, because track never reaches the TRD radius
	return -1;
      }

      Double_t m[2] = {ptrkl2[0] - layer_ref_global2[0], ptrkl2[1] - layer_ref_global2[1]};
      Double_t len_m = TMath::Sqrt(m[0]*m[0] + m[1]*m[1]);
      diff_y += len_m;
      diff_z += TMath::Abs(ptrkl2[2] - layer_ref_global2[2]);
      nLayers++;
    }
  }

  if (nLayers > 0){
    *ydist = diff_y / nLayers;
    *zdist = diff_z / nLayers;
    return 0;
  } else {
    LogError("invalid number of contributing layers (%d) in EstimateTrackDistance()", nLayers);
    return -4;
  }
}
