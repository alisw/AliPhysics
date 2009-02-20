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
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TClonesArray.h>

#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliESDTrdTrack.h"
#include "AliESDEvent.h"

#include "AliTRDReconstructor.h"
#include "AliTRDclusterizer.h"
#include "AliTRDtracker.h"
#include "AliTRDpidESD.h"
#include "AliTRDgtuTrack.h"
#include "AliTRDrawData.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDrecoParam.h"

#include "TTreeStream.h"

#define SETFLG(n,f) ((n) |= f)
#define CLRFLG(n,f) ((n) &= ~f)

ClassImp(AliTRDReconstructor)

TClonesArray *AliTRDReconstructor::fgClusters = 0x0;
//_____________________________________________________________________________
AliTRDReconstructor::AliTRDReconstructor()
  :AliReconstructor()
  ,fSteerParam(0)
{
  // setting default "ON" steering parameters
  // owner of debug streamers 
  SETFLG(fSteerParam, kOwner);
  // write clusters [cw]
  SETFLG(fSteerParam, kWriteClusters);
  // track seeding (stand alone tracking) [sa]
  SETFLG(fSteerParam, kSeeding);
  // PID method in reconstruction (NN) [nn]
  SETFLG(fSteerParam, kSteerPID);
  // number of dEdx slices in the ESD track [8s]
  //SETFLG(fSteerParam, kEightSlices);

  memset(fStreamLevel, 0, kNtasks*sizeof(UChar_t));
  memset(fDebugStream, 0, sizeof(TTreeSRedirector *) * kNtasks);
  // Xe tail cancellation parameters
  fTCParams[0] = 1.156; // r1
  fTCParams[1] = 0.130; // r2
  fTCParams[2] = 0.114; // c1
  fTCParams[3] = 0.624; // c2
  // Ar tail cancellation parameters
  fTCParams[4] = 6.;    // r1
  fTCParams[5] = 0.62;  // r2
  fTCParams[6] = 0.0087;// c1
  fTCParams[7] = 0.07;  // c2
}

//_____________________________________________________________________________
AliTRDReconstructor::AliTRDReconstructor(const AliTRDReconstructor &r)
  :AliReconstructor(r)
  ,fSteerParam(r.fSteerParam)
{
  memcpy(fStreamLevel, r.fStreamLevel, kNtasks*sizeof(UChar_t));
  memcpy(fTCParams, r.fTCParams, 8*sizeof(Double_t));
  memcpy(fDebugStream, r.fDebugStream, sizeof(TTreeSRedirector *) *kNtasks);
  // ownership of debug streamers is not taken
  CLRFLG(fSteerParam, kOwner);
}

//_____________________________________________________________________________
AliTRDReconstructor::~AliTRDReconstructor()
{
  if(fgClusters) {
    fgClusters->Delete(); delete fgClusters;
  }
  if(TESTBIT(fSteerParam, kOwner)){
    for(Int_t itask = 0; itask < kNtasks; itask++)
      if(fDebugStream[itask]) delete fDebugStream[itask];
  }
}


//_____________________________________________________________________________
void AliTRDReconstructor::Init(){
  //
  // Init Options
  //
  SetOption(GetOption());

  AliInfo(Form("\tDigitsConversion       [dc] : %s", fSteerParam&kDigitsConversion?"yes":"no"));
  AliInfo(Form("\tWrite Clusters         [cw] : %s", fSteerParam&kWriteClusters?"yes":"no"));
  AliInfo(Form("\tWrite Online Tracklets [tw] : %s", fSteerParam&kWriteTracklets?"yes":"no"));
  AliInfo(Form("\tDrift Gas Argon        [ar] : %s", fSteerParam&kDriftGas?"yes":"no"));
  AliInfo(Form("\tStand Alone Tracking   [sa] : %s", fSteerParam&kSeeding?"yes":"no"));
  AliInfo(Form("\tHLT         Tracking  [hlt] : %s", fSteerParam&kHLT?"yes":"no"));
  AliInfo(Form("\tCosmic Reconstruction [cos] : %s", fSteerParam&kCosmic?"yes":"no"));
  AliInfo(Form("\tNN PID                 [nn] : %s", fSteerParam&kSteerPID?"yes":"no"));
  AliInfo(Form("\t8 dEdx slices in ESD   [8s] : %s", fSteerParam&kEightSlices?"yes":"no"));
  AliInfo(Form("\tStreaming Levels            : Clusterizer[%d] Tracker[%d] PID[%d]", fStreamLevel[kClusterizer], fStreamLevel[kTracker], fStreamLevel[kPID]));
}

//_____________________________________________________________________________
void AliTRDReconstructor::ConvertDigits(AliRawReader *rawReader
              , TTree *digitsTree) const
{
  //
  // Convert raw data digits into digit objects in a root tree
  //
  AliInfo("Feature not available for the moment."); return;

  AliInfo("Convert raw data digits into digit objects [RawReader -> Digit TTree]");

  AliTRDrawData rawData;
  rawReader->Reset();
  rawReader->Select("TRD");
  AliTRDdigitsManager *manager = rawData.Raw2Digits(rawReader);
  manager->MakeBranch(digitsTree);
  manager->WriteDigits();
  delete manager;

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

  // New (fast) cluster finder
  AliTRDclusterizer clusterer("clusterer","TRD clusterizer",this);
  //clusterer.SetReconstructor(this);                 //     ^| "this" tells the digitsmanager that we are reading raw files
  clusterer.OpenOutput(clusterTree);                  //        it is not strictly necessaray but will give a speed up
  clusterer.SetAddLabels(kFALSE);
  clusterer.Raw2ClustersChamber(rawReader);
  
  if(IsWritingClusters()) return;

  // take over ownership of clusters
  fgClusters = clusterer.RecPoints();
  clusterer.SetClustersOwner(kFALSE);
}

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(TTree *digitsTree
                                    , TTree *clusterTree) const
{
  //
  // Reconstruct clusters
  //

  //AliInfo("Reconstruct TRD clusters from Digits [Digit TTree -> Cluster TTree]");

  AliTRDclusterizer clusterer("clusterer","TRD clusterizer");
  clusterer.SetReconstructor(this);                  //    ^| no this, because we are reading from digitsTree
  clusterer.OpenOutput(clusterTree);                 //       it is necessary to NOT have the "this" here!
  clusterer.ReadDigits(digitsTree);
  clusterer.MakeClusters();

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
        , AliESDEvent* /*esd*/) const
{
  //
  // Fill ESD
  //

}


//_____________________________________________________________________________
void AliTRDReconstructor::SetOption(Option_t *opt)
{
// Read option string into the steer param.
//
// Default steer param values
//
// digits conversion [dc] = false
// drift gas [ar] = false - do not update the number of exponentials in the TC !
// write clusters [cw] = true
// write online tracklets [tw] = false
// track seeding (stand alone tracking) [sa] = true
// PID method in reconstruction (NN) [nn] = true
// 8 dEdx slices in ESD [8s] = false 
// HLT tracking [hlt] = false
// Cosmic Reconstruction [cos] = false
//

  AliReconstructor::SetOption(opt);

  TString s(opt);
  TObjArray *opar = s.Tokenize(",");
  for(Int_t ipar=0; ipar<opar->GetEntriesFast(); ipar++){
    TString sopt(((TObjString*)(*opar)[ipar])->String());
    if(sopt.Contains("dc")){
      SETFLG(fSteerParam, kDigitsConversion);
      if(sopt.Contains("!")) CLRFLG(fSteerParam, kDigitsConversion);
      continue;	
    } else if(sopt.Contains("cw")){ 
      SETFLG(fSteerParam, kWriteClusters);
      if(sopt.Contains("!")) CLRFLG(fSteerParam, kWriteClusters);
      continue;
    } else if(sopt.Contains("sa")){
      SETFLG(fSteerParam, kSeeding);
      if(sopt.Contains("!")) CLRFLG(fSteerParam, kSeeding);
      continue;
    } else if(sopt.Contains("nn")){
      SETFLG(fSteerParam, kSteerPID);
      if(sopt.Contains("!")) CLRFLG(fSteerParam, kSteerPID);
      continue;
    } else if(sopt.Contains("8s")){
      SETFLG(fSteerParam, kEightSlices);
      if(sopt.Contains("!")) CLRFLG(fSteerParam, kEightSlices);
      continue;
    } else if(sopt.Contains("tw")){
      SETFLG(fSteerParam, kWriteTracklets);
      if(sopt.Contains("!")) CLRFLG(fSteerParam, kWriteTracklets);
      continue;	
    } else if(sopt.Contains("ar")){
      SETFLG(fSteerParam, kDriftGas);
      if(sopt.Contains("!")) CLRFLG(fSteerParam, kDriftGas);
      continue;	
    } else if(sopt.Contains("hlt")){
      SETFLG(fSteerParam, kHLT);
      if(sopt.Contains("!")) CLRFLG(fSteerParam, kHLT);
      continue;	
    } else if(sopt.Contains("cos")){
      SETFLG(fSteerParam, kCosmic);
      if(sopt.Contains("!")) CLRFLG(fSteerParam, kCosmic);
    } else if(sopt.Contains("sl")){
      TObjArray *stl = sopt.Tokenize("_");
      if(stl->GetEntriesFast() < 3) continue;
      TString taskstr(((TObjString*)(*stl)[1])->String());
      TString levelstring(((TObjString*)(*stl)[2])->String());
      // Set the stream Level
      Int_t level = levelstring.Atoi();
      ETRDReconstructorTask task = kTracker;
      if(taskstr.CompareTo("raw") == 0) task = kRawReader;	
      else if(taskstr.CompareTo("cl") == 0) task = kClusterizer;	
      else if(taskstr.CompareTo("tr") == 0) task = kTracker;
      else if(taskstr.CompareTo("pi") == 0) task = kPID;
      SetStreamLevel(level, task);
      continue;
    }
  }
}

//_____________________________________________________________________________
void AliTRDReconstructor::SetStreamLevel(Int_t level, ETRDReconstructorTask task){
  //
  // Set the Stream Level for one of the tasks Clusterizer, Tracker or PID
  //
  TString taskname[4] = {"RawReader", "Clusterizer", "Tracker", "PID"};
  const Int_t minLevel[4] = {1, 1, 2, 1}; // the minimum debug level upon which a debug stream is created for different tasks
  //AliInfo(Form("Setting Stream Level for Task %s to %d", taskname.Data(),level));
  fStreamLevel[(Int_t)task] = level;
  // Initialize DebugStreamer if not yet done
  if(level >= minLevel[task] && !fDebugStream[task]){
    TDirectory *savedir = gDirectory;
    fDebugStream[task] = new TTreeSRedirector(Form("TRD.%sDebug.root", taskname[task].Data()));
    savedir->cd();
    SETFLG(fSteerParam, kOwner);
  }
}
