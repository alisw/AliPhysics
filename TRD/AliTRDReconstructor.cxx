/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class for TRD reconstruction                                              //
//                                                                           //   
// For the special options which can be used during reconstruction and their //
//  default values pls. see function SetOption().                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObjString.h>
#include <TObjArray.h>
#include <TTreeStream.h>
#include <TDirectory.h>
#include <TRef.h>

#include "AliRawReader.h"
#include "AliLog.h"

#include "AliTRDReconstructor.h"
#include "AliTRDclusterizer.h"
#include "AliTRDrawData.h"
#include "AliTRDrawStream.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDtrackerV1.h"
#include "AliESDEvent.h"
#include "AliESDTrdTrack.h"
#include "AliESDTrdTracklet.h"
#include "AliESDTrdTrigger.h"
#include "AliTRDtrackletWord.h"
#include "AliTRDtrackletMCM.h"
#include "AliTRDonlineTrackMatching.h"

#define SETFLG(n,f) ((n) |= f)
#define CLRFLG(n,f) ((n) &= ~f)

ClassImp(AliTRDReconstructor)

AliESDTrdTrigger AliTRDReconstructor::fgTriggerFlags;
AliTRDonlineTrackMatching AliTRDReconstructor::fgOnlineTrackMatcher;
TClonesArray *AliTRDReconstructor::fgClusters = NULL;
TClonesArray *AliTRDReconstructor::fgTracklets = NULL;
TClonesArray *AliTRDReconstructor::fgTracks = NULL;
Char_t const * AliTRDReconstructor::fgSteerNames[kNsteer] = {
  "DigitsConversion       "
 ,"Write Clusters         "
 ,"Write Online Tracklets "
 ,"Stand Alone Tracking   "
 ,"HLT Mode              "
 ,"Process Online Trklts  "
 ,"Debug Streaming       "
 ,"Cl. Radial Correction  "
};
Char_t const * AliTRDReconstructor::fgSteerFlags[kNsteer] = {
  "dc"// digits conversion [false]
 ,"cw"// write clusters [true]
 ,"tw"// write online tracklets [false]
 ,"sa"// track seeding (stand alone tracking) [true]
 ,"hlt"// HLT reconstruction [false]
 ,"tp"// also use online tracklets for reconstruction [false]
 ,"deb"// Write debug stream [false]
 ,"cc" // Cluster radial correction during reconstruction [false]
};
Char_t const * AliTRDReconstructor::fgTaskNames[AliTRDrecoParam::kTRDreconstructionTasks] = {
  "Clusterizer"
 ,"Tracker"
 ,"PID"
};
Char_t const * AliTRDReconstructor::fgTaskFlags[AliTRDrecoParam::kTRDreconstructionTasks] = {
  "cl"
 ,"tr"
 ,"pd"
};
Int_t AliTRDReconstructor::fgNTimeBins = -1;
const  Float_t  AliTRDReconstructor::fgkMinClustersInTrack =  0.5;  //
const  Float_t  AliTRDReconstructor::fgkLabelFraction      =  0.8;  //
const  Double_t AliTRDReconstructor::fgkMaxChi2            = 12.0;  //
const  Double_t AliTRDReconstructor::fgkMaxSnp             =  0.95; // Maximum local sine of the azimuthal angle
const  Double_t AliTRDReconstructor::fgkMaxStep            =  2.0;  // Maximal step size in propagation
const Double_t  AliTRDReconstructor::fgkEpsilon            = 1.e-5;                  // Precision of radial coordinate

//_____________________________________________________________________________
AliTRDReconstructor::AliTRDReconstructor()
  :AliReconstructor()
  ,fSteerParam(0)
  ,fClusterizer(NULL)
{
  // setting default "ON" steering parameters
  // owner of debug streamers 
  SETFLG(fSteerParam, kOwner);
  // write clusters [cw]
  SETFLG(fSteerParam, kWriteClusters);
  // track seeding (stand alone tracking) [sa]
  //SETFLG(fSteerParam, kSeeding);
  // Cluster radial correction during reconstruction [cc]
  //SETFLG(fSteerParam, kClRadialCorr);
  memset(fDebugStream, 0, sizeof(TTreeSRedirector *) * AliTRDrecoParam::kTRDreconstructionTasks);
}

//_____________________________________________________________________________
AliTRDReconstructor::~AliTRDReconstructor()
{
  //
  // Destructor
  //

  if(fClusterizer){
    delete fClusterizer;
    fClusterizer = NULL;
  }
  if(fgClusters) {
    fgClusters->Delete();
    delete fgClusters;
    fgClusters = NULL;
  }
  if(fgTracklets) {
    fgTracklets->Delete();
    delete fgTracklets;
    fgTracklets = NULL;
  }
  if(fgTracks) {
    fgTracks->Delete();
    delete fgTracks;
    fgTracks = NULL;
  }
  if(fSteerParam&kOwner){
    for(Int_t itask = 0; itask < AliTRDrecoParam::kTRDreconstructionTasks; itask++)
      if(fDebugStream[itask]) delete fDebugStream[itask];
  }
}


//_____________________________________________________________________________
void AliTRDReconstructor::Init(){
  //
  // Init Options
  //
  SetOption(GetOption());
  Options(fSteerParam);

  if(!fClusterizer){
    fClusterizer = new AliTRDclusterizer(fgTaskNames[AliTRDrecoParam::kClusterizer], fgTaskNames[AliTRDrecoParam::kClusterizer]);
    fClusterizer->SetReconstructor(this);
    AliInfo(Form("Build TRD clusterizer[%p]", (void*)fClusterizer));
  }
  
  // Make Debug Streams when Debug Streaming
  if(IsDebugStreaming()){
    for(Int_t task = 0; task < AliTRDrecoParam::kTRDreconstructionTasks; task++){
      TDirectory *savedir = gDirectory;
      fDebugStream[task] = new TTreeSRedirector(Form("TRD.Debug%s.root", fgTaskNames[task]));
      savedir->cd();
      SETFLG(fSteerParam, kOwner);
    }
  }
}

//_____________________________________________________________________________
void AliTRDReconstructor::ConvertDigits(AliRawReader *rawReader
              , TTree *digitsTree) const
{
  //
  // Convert raw data digits into digit objects in a root tree
  //

  AliDebug(1, "Convert raw data digits into digit objects [RawReader -> Digit TTree]");
  AliDebug(2, Form("clusters[%p] tracklets[%p] tracks[%p]", (void*)fgClusters, (void*)fgTracklets, (void*)fgTracks));

  rawReader->Reset();
  rawReader->Select("TRD");
  ResetContainers();

  AliTRDrawData rawData;

  AliTRDdigitsManager *manager = rawData.Raw2Digits(rawReader);
  manager->MakeBranch(digitsTree);
  manager->WriteDigits();
  delete manager;

  for (Int_t iSector = 0; iSector < 18; iSector++) fgTriggerFlags.SetFlags(iSector, rawData.GetTriggerFlags(iSector));
}

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(AliRawReader *rawReader
                                    , TTree *clusterTree) const
{
  //
  // Reconstruct clusters
  //

  AliDebug(1, "Reconstruct TRD clusters from RAW data [RawReader -> Cluster TTree]");
  AliDebug(2, Form("clusters[%p] tracklets[%p] tracks[%p]", (void*)fgClusters, (void*)fgTracklets, (void*)fgTracks));
  if(!fClusterizer){
    AliFatal("Clusterizer not available!");
    return;
  }
  rawReader->Reset();
  rawReader->Select("TRD");
  ResetContainers();
  fClusterizer->OpenOutput(clusterTree);
  fClusterizer->SetUseLabels(kFALSE);
  fClusterizer->SetStoreRawSignals(kTRUE);
  fClusterizer->ResetRecPoints();
  fClusterizer->Raw2ClustersChamber(rawReader);
  fgNTimeBins = fClusterizer->GetNTimeBins();
  for (Int_t iSector = 0; iSector < 18; iSector++) fgTriggerFlags.SetFlags(iSector, fClusterizer->GetTriggerFlags(iSector));
}

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(TTree *digitsTree
                                    , TTree *clusterTree) const
{
  //
  // Reconstruct clusters
  //

  AliDebug(1, "Reconstruct TRD clusters from Digits [Digit TTree -> Cluster TTree]");
  AliDebug(2, Form("Start :: clusters[%p] tracklets[%p] tracks[%p]", (void*)fgClusters, (void*)fgTracklets, (void*)fgTracks));
  if(!fClusterizer){
    AliFatal("Clusterizer not available!");
    return;
  }

  ResetContainers();
  AliTRDclusterizer clusterizer(fgTaskNames[AliTRDrecoParam::kClusterizer], fgTaskNames[AliTRDrecoParam::kClusterizer]);
  clusterizer.SetReconstructor(this);
  clusterizer.SetUseLabels(kTRUE);
  clusterizer.SetStoreRawSignals(kTRUE);
  clusterizer.OpenOutput(clusterTree);
  clusterizer.ReadDigits(digitsTree);
  clusterizer.ReadTracklets();
  clusterizer.ReadTracks();
  clusterizer.MakeClusters();
  fgNTimeBins = clusterizer.GetNTimeBins();
}

//_____________________________________________________________________________
AliTracker *AliTRDReconstructor::CreateTracker() const
{
  //
  // Create a TRD tracker
  //

  //return new AliTRDtracker(NULL);
  AliTRDtrackerV1 *tracker = new AliTRDtrackerV1();
  tracker->SetReconstructor(this);
  AliInfo(Form("Build TRD tracker[%p]", (void*)tracker));
  return tracker;

}

//_____________________________________________________________________________
void AliTRDReconstructor::FillESD(TTree* /*digitsTree*/
        , TTree* /*clusterTree*/
        , AliESDEvent* esd) const
{
  //
  // Fill ESD
  //

  // ----- filling tracklets -----
  AliDebug(1, Form("Loading onl.tracklets(%i) to ESD", fgTracklets ? fgTracklets->GetEntriesFast() : 0));
  Int_t trackletIndex[1080] = { 0 };
  TList trklList;
  AliTRDrawStream::SortTracklets(fgTracklets, trklList, trackletIndex);
  TIter trackletIter(&trklList);
  while (AliTRDtrackletBase* tracklet = (AliTRDtrackletBase*) trackletIter()) {
    Int_t label = -2; // mark raw tracklets with label -2
    if (AliTRDtrackletMCM *trklMCM = dynamic_cast<AliTRDtrackletMCM*> (tracklet)) label = trklMCM->GetLabel();
    esd->AddTrdTracklet(new AliESDTrdTracklet(tracklet->GetTrackletWord(), tracklet->GetHCId(), label));
  }

  // ----- filling GTU tracks -----
  AliDebug(1, Form("Loading gtu.tracks(%i) to ESD", fgTracks ? fgTracks->GetEntriesFast() : 0));
  if (fgTracks) {
    for (Int_t iTrack = 0; iTrack < fgTracks->GetEntriesFast(); iTrack++) {
      AliESDTrdTrack *trdTrack = (AliESDTrdTrack*) ((*fgTracks)[iTrack]);

      UInt_t stack = trdTrack->GetStack();

      Int_t refIndex[6];
      AliTRDrawStream::AssignTracklets(trdTrack, trackletIndex, refIndex);

      for (Int_t iLayer = 0; iLayer < 6; ++iLayer) {
        Int_t det = trdTrack->GetSector()*30 + stack*6 + iLayer;
        AliESDTrdTracklet *trkl = refIndex[iLayer] > -1 ? esd->GetTrdTracklet(refIndex[iLayer]) : 0x0;
        if (trkl) {
          AliDebug(5, Form("adding tracklet with index %i: 0x%08x",
              refIndex[iLayer], trkl->GetTrackletWord()));
          if (trkl->GetDetector() != det)
            AliError(Form("inconsistent assignment of tracklet 0x%08x in det %i to track in %i",
              trkl->GetTrackletWord(), trkl->GetDetector(), det));
          trdTrack->AddTrackletReference(trkl, iLayer);
        }
      }
      // only add the track when it's complete (including tracklet references)
      esd->AddTrdTrack(trdTrack);
    }
  }
  esd->SetTrdTrigger(&fgTriggerFlags);

  // ----- matching GTU tracks to global tracks -----
  AliDebug(1, Form("TRD track matching with %i ESD, %i TRD tracks",
		   esd->GetNumberOfTracks(), esd->GetNumberOfTrdTracks()));
  fgOnlineTrackMatcher.ProcessEvent(esd);
}

//_____________________________________________________________________________
void AliTRDReconstructor::SetOption(Option_t *opt)
{
  //
  // Read option string into the steer param.
  //
  // The following string options are available during reconstruction.
  // In square brackets the default values are given.
  //   "dc"  : digits conversion [false]
  //   "cw"  : write clusters [true]
  //   "tw"  : write online tracklets [false]
  //   "sa"  : track seeding (stand alone tracking) [true]
  //   "hlt" : HLT reconstruction [false]
  //   "tp"  : also use online tracklets for reconstruction [false]
  //   "deb" : Write debug stream [false]
  //   "cc"  : Cluster radial correction during reconstruction [false]
  //
  // To check the actual options used during reconstruction include the following line in your rec.C script
  // AliLog::SetClassDebugLevel("AliTRDReconstructor", 1);

  AliReconstructor::SetOption(opt);

  TString s(opt);
  TObjArray *opar = s.Tokenize(",");
  for(Int_t ipar=0; ipar<opar->GetEntriesFast(); ipar++){
    Bool_t processed = kFALSE;
    TString sopt(((TObjString*)(*opar)[ipar])->String());
    for(Int_t iopt=0; iopt<kNsteer; iopt++){
      if(!sopt.Contains(fgSteerFlags[iopt])) continue;
      SETFLG(fSteerParam, BIT(iopt));
      if(sopt.Contains("!")) CLRFLG(fSteerParam, BIT(iopt));
      processed = kTRUE;
      break;	
    }
    if(processed) continue;

    AliWarning(Form("Unknown option flag %s.", sopt.Data()));
  }
  opar->Delete();
  delete opar;
}

//_____________________________________________________________________________
void AliTRDReconstructor::Options(UInt_t steer)
{
  //
  // Print the options
  //

  for(Int_t iopt=0; iopt<kNsteer; iopt++){
    AliDebugGeneral("AliTRDReconstructor", 1, Form(" %s[%s]%s", fgSteerNames[iopt], fgSteerFlags[iopt], steer ?(((steer>>iopt)&1)?" : ON":" : OFF"):""));
  }
}


//_____________________________________________________________________________
TClonesArray* AliTRDReconstructor::GetClusters()
{
// Build/ Retrieve cluster array
  if(!fgClusters){
    fgClusters = new TClonesArray("AliTRDcluster", Int_t(GetRecoParam()->GetNClusters()));
    fgClusters->SetOwner();
    AliInfoGeneral("AliTRDReconstructor", Form("Allocate cluster array @ %p", (void*)fgClusters));
  }
  return fgClusters;
}

//_____________________________________________________________________________
TClonesArray* AliTRDReconstructor::GetTracklets(const char *trkltype)
{
// Build/ Retrieve online tracklets array

  if (trkltype && strlen(trkltype) > 0) {
    if(fgTracklets && (TClass::GetClass(trkltype) != fgTracklets->GetClass())){
      fgTracklets->Delete();
      delete fgTracklets;
    }
    if (!fgTracklets) {
      fgTracklets = new TClonesArray(trkltype, 200);
      fgTracklets->SetOwner(kTRUE);
      AliInfoGeneral("AliTRDReconstructor", Form("Allocate online tracklets[%s] array @ %p", trkltype, (void*)fgTracklets));
    }
  }
  return fgTracklets;
}

//_____________________________________________________________________________
TClonesArray* AliTRDReconstructor::GetTracks()
{
// Build/ Retrieve cluster array
  if(!fgTracks){
    fgTracks = new TClonesArray("AliESDTrdTrack", 100);
    fgTracks->SetOwner();
    AliInfoGeneral("AliTRDReconstructor", Form("Allocate online tracks array @ %p", (void*)fgTracks));
  }
  return fgTracks;
}

//_____________________________________________________________________________
void AliTRDReconstructor::ResetContainers() const
{
// prepare data containers for a new event

  if(fgClusters){
    AliDebug(1, Form("Removing %5d clusters @ %p", fgClusters->GetEntriesFast(), (void*)fgClusters));
    fgClusters->Clear();
  }
  if(fgTracklets){
    AliDebug(1, Form("Removing %3d online tracklets @ %p", fgTracklets->GetEntriesFast(), (void*)fgTracklets));
    fgTracklets->Clear();
  }
  if(fgTracks){
    AliDebug(1, Form("Removing %3d online tracks @ %p", fgTracks->GetEntriesFast(), (void*)fgTracks));
    fgTracks->Clear();
  }
}
