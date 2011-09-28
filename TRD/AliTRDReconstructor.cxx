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

#define SETFLG(n,f) ((n) |= f)
#define CLRFLG(n,f) ((n) &= ~f)

ClassImp(AliTRDReconstructor)

AliESDTrdTrigger AliTRDReconstructor::fgTriggerFlags;
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
  if(fClusterizer){
    delete fClusterizer;
    fClusterizer = NULL;
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

  //AliInfo("Convert raw data digits into digit objects [RawReader -> Digit TTree]");

  AliTRDrawData rawData;
  rawReader->Reset();
  rawReader->Select("TRD");
  AliTRDdigitsManager *manager = rawData.Raw2Digits(rawReader);
  manager->MakeBranch(digitsTree);
  manager->WriteDigits();
  delete manager;

  // take over ownership of online tracklets
  fgTracklets = rawData.TrackletsArray();
  rawData.SetTrackletsOwner(0x0);

  // take over GTU tracks
  fgTracks = rawData.TracksArray();
  rawData.SetTracksOwner(0x0);

  for (Int_t iSector = 0; iSector < 18; iSector++) {
    fgTriggerFlags.SetFlags(iSector, rawData.GetTriggerFlags(iSector));
  }
}

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(AliRawReader *rawReader
                                    , TTree *clusterTree) const
{
  //
  // Reconstruct clusters
  //

  //AliInfo("Reconstruct TRD clusters from RAW data [RawReader -> Cluster TTree]");


  rawReader->Reset();
  rawReader->Select("TRD");

  if(!fClusterizer){
    AliFatal("Clusterizer not available!");
    return;
  }

  fClusterizer->ResetRecPoints();

  fClusterizer->OpenOutput(clusterTree);
  fClusterizer->SetUseLabels(kFALSE);
  fClusterizer->Raw2ClustersChamber(rawReader);
  
  fgNTimeBins = fClusterizer->GetNTimeBins();
  
  // take over ownership of online tracklets
  fgTracklets = fClusterizer->TrackletsArray();
  fClusterizer->SetTrackletsOwner(kFALSE);

  // take over GTU tracks
  fgTracks = fClusterizer->TracksArray();
  fClusterizer->SetTracksOwner(kFALSE);

  for (Int_t iSector = 0; iSector < 18; iSector++) {
    fgTriggerFlags.SetFlags(iSector, fClusterizer->GetTriggerFlags(iSector));
  }

  if(IsWritingClusters()) return;

  // take over ownership of clusters
  fgClusters = fClusterizer->RecPoints();
  fClusterizer->SetClustersOwner(kFALSE);
}

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(TTree *digitsTree
                                    , TTree *clusterTree) const
{
  //
  // Reconstruct clusters
  //

  //AliInfo("Reconstruct TRD clusters from Digits [Digit TTree -> Cluster TTree]");
  
  AliTRDclusterizer clusterer(fgTaskNames[AliTRDrecoParam::kClusterizer], fgTaskNames[AliTRDrecoParam::kClusterizer]);
  clusterer.SetReconstructor(this);
  clusterer.OpenOutput(clusterTree);
  clusterer.ReadDigits(digitsTree);
  clusterer.MakeClusters();

  // read tracklets and tracks if not done during reading of raw data
  if (!fgTracklets) {
    clusterer.ReadTracklets();
    fgTracklets = clusterer.TrackletsArray();
    clusterer.SetTrackletsOwner(kFALSE);
  }
  if (!fgTracks) {
    clusterer.ReadTracks();
    fgTracks = clusterer.TracksArray();
    clusterer.SetTracksOwner(kFALSE);
  }

  fgNTimeBins = clusterer.GetNTimeBins();

  if(IsWritingClusters()) return;

  // take over ownership of clusters
  fgClusters = clusterer.RecPoints();
  clusterer.SetClustersOwner(kFALSE);
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
  Int_t trackletIndex[1080] = { 0 };
  AliDebug(1, Form("Filling tracklets from %p (%i)",
		   fgTracklets, fgTracklets ? fgTracklets->GetEntriesFast() : 0));
  if (fgTracklets) {
    Int_t nTracklets = fgTracklets->GetEntriesFast();

    Int_t lastHC = -1;
    for (Int_t iTracklet = 0; iTracklet < nTracklets; iTracklet++) {
      AliTRDtrackletBase *trkl = (AliTRDtrackletBase*) ((*fgTracklets)[iTracklet]);
      Int_t hc = trkl->GetHCId();
      //      hc += trkl->GetY() > 0 ? 1 : 0;
      if ((hc < 0) || (hc >= 1080)) {
	AliError(Form("HC for tracklet: 0x%08x out of range: %i", trkl->GetTrackletWord(), trkl->GetHCId()));
	continue;
      }
      AliDebug(5, Form("hc: %4i : 0x%08x z: %2i", hc, trkl->GetTrackletWord(), trkl->GetZbin()));
      if (hc != lastHC) {
	AliDebug(2, Form("set tracklet index for HC %i to %i", hc, iTracklet));
	trackletIndex[hc] = iTracklet + 1;
	lastHC = hc;
      }
    }

    for (Int_t iDet = 0; iDet < 540; iDet++) {
      Int_t trklIndexA = trackletIndex[2*iDet + 0] - 1;
      Int_t trklIndexB = trackletIndex[2*iDet + 1] - 1;
      Int_t trklIndex  = esd->GetNumberOfTrdTracklets();
      AliTRDtrackletBase *trklA = trklIndexA > -1 ? (AliTRDtrackletBase*) ((*fgTracklets)[trklIndexA]) : 0x0;
      AliTRDtrackletBase *trklB = trklIndexB > -1 ? (AliTRDtrackletBase*) ((*fgTracklets)[trklIndexB]) : 0x0;
      while (trklA != 0x0 || trklB != 0x0) {
	AliDebug(5, Form("det %i - A: %i/%i -> %p, B: %i/%i -> %p",
			 iDet, trklIndexA, nTracklets, trklA, trklIndexB, nTracklets, trklB));
	if (trklA == 0x0) {
	  AliESDTrdTracklet esdTracklet(trklB->GetTrackletWord(), trklB->GetHCId());
	  esd->AddTrdTracklet(&esdTracklet);
	  trklIndexB++;
	  trklB = trklIndexB < nTracklets ? (AliTRDtrackletBase*) ((*fgTracklets)[trklIndexB]) : 0x0;
	  if (trklB && trklB->GetHCId() != 2*iDet + 1)
	    trklB = 0x0;
	}
	else if (trklB == 0x0) {
	  AliESDTrdTracklet esdTracklet(trklA->GetTrackletWord(), trklA->GetHCId());
	  esd->AddTrdTracklet(&esdTracklet);
	  trklIndexA++;
	  trklA = trklIndexA < nTracklets ? (AliTRDtrackletBase*) ((*fgTracklets)[trklIndexA]) : 0x0;
	  if (trklA && trklA->GetHCId() != 2*iDet)
	    trklA = 0x0;
	}
	else {
	  if ((trklA->GetZbin() < trklB->GetZbin()) ||
	      ((trklA->GetZbin() == trklB->GetZbin()) && (trklA->GetYbin() < trklB->GetYbin()))) {
	    AliESDTrdTracklet esdTracklet(trklA->GetTrackletWord(), trklA->GetHCId());
	    esd->AddTrdTracklet(&esdTracklet);
	    trklIndexA++;
	    trklA = trklIndexA < nTracklets ? (AliTRDtrackletBase*) ((*fgTracklets)[trklIndexA]) : 0x0;
	    if (trklA && trklA->GetHCId() != 2*iDet)
	      trklA = 0x0;
	  }
	  else {
	    AliESDTrdTracklet esdTracklet(trklB->GetTrackletWord(), trklB->GetHCId());
	    esd->AddTrdTracklet(&esdTracklet);
	    trklIndexB++;
	    trklB = trklIndexB < nTracklets ? (AliTRDtrackletBase*) ((*fgTracklets)[trklIndexB]) : 0x0;
	    if (trklB && trklB->GetHCId() != 2*iDet + 1)
	      trklB = 0x0;
	  }
	}
      }
      // updating tracklet indices as in ESD
      if (esd->GetNumberOfTrdTracklets() != trklIndex) {
	trackletIndex[2*iDet + 0] = trackletIndex[2*iDet + 1] = trklIndex;
      }
      else
	trackletIndex[2*iDet + 0] = trackletIndex[2*iDet + 1] = -1;
    }
  }

  // ----- filling GTU tracks -----
  AliDebug(1, Form("Now filling ESD with GTU tracks from %p (%i)",
		   fgTracks, fgTracks ? fgTracks->GetEntriesFast() : 0));
  if (fgTracks) {
    for (Int_t iTrack = 0; iTrack < fgTracks->GetEntriesFast(); iTrack++) {
      AliESDTrdTrack *trdTrack = (AliESDTrdTrack*) ((*fgTracks)[iTrack]);

      UInt_t mask  = trdTrack->GetLayerMask();
      UInt_t stack = trdTrack->GetStack();

      for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
	if (mask & (1 << iLayer)) {

	  Int_t det = trdTrack->GetSector()*30 + stack*6 + iLayer;
	  Int_t idx = trdTrack->GetTrackletIndex(iLayer);

	  if ((det < 0) || (det > 539)) {
	    AliError(Form("Invalid detector no. from track: %i", 2*det));
	    continue;
	  }
	  if (trackletIndex[2*det] >= 0) {
	    if ((trackletIndex[2*det] + idx > -1) &&
		(trackletIndex[2*det] + idx < esd->GetNumberOfTrdTracklets())) {
	      AliESDTrdTracklet *trkl = esd->GetTrdTracklet(trackletIndex[2*det] + idx);
	      if (trkl) {
		AliDebug(5, Form("adding tracklet: 0x%08x", trkl->GetTrackletWord()));
		trdTrack->AddTrackletReference(trkl, iLayer);
	      }
	    }
	  }
	}
      }
      // only add the track when it's complete (including tracklet references)
      esd->AddTrdTrack(trdTrack);
    }
  }

  esd->SetTrdTrigger(&fgTriggerFlags);

  // clearing variables for next event
  fgTracklets = 0x0;
  fgTracks = 0x0;
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

