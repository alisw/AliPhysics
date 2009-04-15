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
Char_t* AliTRDReconstructor::fgSteerNames[kNsteer] = {
  "DigitsConversion       "
 ,"Tail Cancellation      "
 ,"Clusters LUT           "
 ,"Clusters GAUSS         "
 ,"Clusters Sharing       "
 ,"NN PID                 "
 ,"8 dEdx slices in ESD   "
 ,"Write Clusters         "
 ,"Write Online Tracklets "
 ,"Drift Gas Argon        "
 ,"Stand Alone Tracking   "
 ,"Vertex Constrain       "
 ,"Tracklet Improve       "
 ,"HLT Mode              "
 ,"Cosmic Reconstruction "
};
Char_t* AliTRDReconstructor::fgSteerFlags[kNsteer] = {
  "dc"// digits conversion [false]
 ,"tc"// apply tail cancellation [true]
 ,"lut"// look-up-table for cluster shape in the r-phi direction
 ,"gs"// gauss cluster shape in the r-phi direction
 ,"sh"// cluster sharing between tracks
 ,"nn"// PID method in reconstruction (NN) [true]
 ,"8s"// 8 dEdx slices in ESD [true] 
 ,"cw"// write clusters [true]
 ,"tw"// write online tracklets [false]
 ,"ar"// drift gas [false] - do not update the number of exponentials in the TC !
 ,"sa"// track seeding (stand alone tracking) [true]
 ,"vc"// vertex constrain on stand alone track finder [false]
 ,"ti"// improve tracklets in stand alone track finder [true]
 ,"hlt"// HLT reconstruction [false]
 ,"cos"// Cosmic Reconstruction [false]
};
Char_t* AliTRDReconstructor::fgTaskNames[kNtasks] = {
  "RawReader"
 ,"Clusterizer"
 ,"Tracker"
 ,"PID"
};
Char_t* AliTRDReconstructor::fgTaskFlags[kNtasks] = {
  "rr"
 ,"cl"
 ,"tr"
 ,"pd"
};

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
  SETFLG(fSteerParam, kEightSlices);
  // vertex constrain for stand alone track finder
  SETFLG(fSteerParam, kVertexConstrained);
  // improve tracklets for stand alone track finder
  SETFLG(fSteerParam, kImproveTracklet);
  // use look up table for cluster r-phi position
  SETFLG(fSteerParam, kLUT);
  // use tail cancellation
  SETFLG(fSteerParam, kTC);

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
  if(fSteerParam&kOwner){
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
  Options(fSteerParam, fStreamLevel);
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
  AliTRDclusterizer clusterer(fgTaskNames[kClusterizer], fgTaskNames[kClusterizer]);
  clusterer.SetReconstructor(this);
  clusterer.OpenOutput(clusterTree);
  clusterer.SetUseLabels(kFALSE);
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

  AliTRDclusterizer clusterer(fgTaskNames[kClusterizer], fgTaskNames[kClusterizer]);
  clusterer.SetReconstructor(this);
  clusterer.OpenOutput(clusterTree);
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

  AliReconstructor::SetOption(opt);

  TString s(opt);
  TObjArray *opar = s.Tokenize(",");
  for(Int_t ipar=0; ipar<opar->GetEntriesFast(); ipar++){
    Bool_t PROCESSED = kFALSE;
    TString sopt(((TObjString*)(*opar)[ipar])->String());
    for(Int_t iopt=0; iopt<kNsteer; iopt++){
      if(!sopt.Contains(fgSteerFlags[iopt])) continue;
      SETFLG(fSteerParam, BIT(iopt));
      if(sopt.Contains("!")) CLRFLG(fSteerParam, BIT(iopt));
      PROCESSED = kTRUE;
      break;	
    }
    // extra rules
    if(sopt.Contains("gs") && !sopt.Contains("!")){
      CLRFLG(fSteerParam, kLUT); PROCESSED = kTRUE;
    }

    if(PROCESSED) continue;

    if(sopt.Contains("sl")){
      TObjArray *stl = sopt.Tokenize("_");
      if(stl->GetEntriesFast() < 3) continue;
      TString taskstr(((TObjString*)(*stl)[1])->String());
      TString levelstring(((TObjString*)(*stl)[2])->String());
      Int_t level = levelstring.Atoi();

      // Set the stream Level
      PROCESSED = kFALSE;
      for(Int_t it=0; it<kNtasks; it++){
        if(taskstr.CompareTo(fgTaskFlags[it]) != 0) continue;
        SetStreamLevel(level, ETRDReconstructorTask(it));
        PROCESSED = kTRUE;
      }
    } 
    if(PROCESSED) continue;

    AliWarning(Form("Unknown option flag %s.", sopt.Data()));
  }
}

//_____________________________________________________________________________
void AliTRDReconstructor::SetStreamLevel(Int_t level, ETRDReconstructorTask task){
  //
  // Set the Stream Level for one of the tasks Clusterizer, Tracker or PID
  //
  const Int_t minLevel[4] = {1, 1, 2, 1}; // the minimum debug level upon which a debug stream is created for different tasks
  //AliInfo(Form("Setting Stream Level for Task %s to %d", taskname.Data(),level));
  fStreamLevel[(Int_t)task] = level;
  // Initialize DebugStreamer if not yet done
  if(level >= minLevel[task] && !fDebugStream[task]){
    TDirectory *savedir = gDirectory;
    fDebugStream[task] = new TTreeSRedirector(Form("TRD.Debug%s.root", fgTaskNames[task]));
    savedir->cd();
    SETFLG(fSteerParam, kOwner);
  }
}

//_____________________________________________________________________________
void AliTRDReconstructor::Options(UInt_t steer, UChar_t *stream)
{
  for(Int_t iopt=0; iopt<kNsteer; iopt++){
    AliInfoGeneral("AliTRDReconstructor", Form(" %s[%s]%s", fgSteerNames[iopt], fgSteerFlags[iopt], steer ?(((steer>>iopt)&1)?" : ON":" : OFF"):""));
  }
  AliInfoGeneral("AliTRDReconstructor", " Debug Streaming"); 
  for(Int_t it=0; it<kNtasks; it++) 
    AliInfoGeneral("AliTRDReconstructor", Form(" %s [sl_%s] %d", fgTaskNames[it], fgTaskFlags[it], stream ? stream[it] : 0));
}

