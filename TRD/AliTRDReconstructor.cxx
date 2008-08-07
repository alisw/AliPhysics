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

ClassImp(AliTRDReconstructor)


//_____________________________________________________________________________
AliTRDReconstructor::AliTRDReconstructor()
  :AliReconstructor()
  ,fSteerParam(0x00000007)
{
  memset(fStreamLevel, 0, 5*sizeof(UChar_t));
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
  memcpy(fStreamLevel, r.fStreamLevel, 5*sizeof(UChar_t));
  memcpy(fTCParams, r.fTCParams, 8*sizeof(Double_t));
}

//_____________________________________________________________________________
void AliTRDReconstructor::Init(){
	//
	// Init Options
	//
	SetOption(GetOption());
}

//_____________________________________________________________________________
void AliTRDReconstructor::ConvertDigits(AliRawReader *rawReader
				      , TTree *digitsTree) const
{
  //
  // Convert raw data digits into digit objects in a root tree
  //

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

  AliInfo("Reconstruct TRD clusters from RAW data [RawReader -> Cluster TTree]");


  rawReader->Reset();
  rawReader->Select("TRD");

  // New (fast) cluster finder
  AliTRDclusterizer clusterer("clusterer","TRD clusterizer");
  clusterer.SetReconstructor(this);
  clusterer.OpenOutput(clusterTree);
  clusterer.SetAddLabels(kFALSE);
  clusterer.Raw2ClustersChamber(rawReader);

}

//_____________________________________________________________________________
void AliTRDReconstructor::Reconstruct(TTree *digitsTree
                                    , TTree *clusterTree) const
{
  //
  // Reconstruct clusters
  //

  AliInfo("Reconstruct TRD clusters from Digits [Digit TTree -> Cluster TTree]");

  AliTRDclusterizer clusterer("clusterer","TRD clusterizer");
  clusterer.SetReconstructor(this);
  clusterer.OpenOutput(clusterTree);
  clusterer.ReadDigits(digitsTree);
  clusterer.MakeClusters();

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
// write clusters [cw] = true
// track seeding (stand alone tracking) [sa] = true
// PID method in reconstruction (NN) [nn] = true
// write online tracklets [tw] = false
// drift gas [ar] = false
//
  fSteerParam = 0x00000007;

  TString s(opt);
  TObjArray *opar = s.Tokenize(",");
  for(Int_t ipar=0; ipar<opar->GetEntriesFast(); ipar++){
    TString sopt(((TObjString*)(*opar)[ipar])->String());
    if(sopt.Contains("!cw")){ 
      fSteerParam &= ~kWriteClusters;
      continue;
    } else if(sopt.Contains("!sa")){
      fSteerParam &= ~kSeeding;
      continue;
    } else if(sopt.Contains("!nn")){
      fSteerParam &= ~kSteerPID;
      continue;
    } else if(sopt.Contains("tw")){
      fSteerParam |= kWriteTracklets;
      continue;	
    } else if(sopt.Contains("ar")){
      fSteerParam |= kDriftGas;
      continue;	
    } else if(sopt.Contains("sl")){
    	TObjArray *stl = sopt.Tokenize("_");
    	if(stl->GetEntriesFast() < 3) continue;
			TString taskstr(((TObjString*)(*stl)[1])->String());
			TString levelstring(((TObjString*)(*stl)[2])->String());
    	// Set the stream Level
    	Int_t level = levelstring.Atoi();
    	AliTRDReconstructorTask task = kTracker;
    	if(taskstr.CompareTo("cl") == 0) task = kClusterizer;	
    	else if(taskstr.CompareTo("tr") == 0) task = kTracker;
    	else if(taskstr.CompareTo("pi") == 0) task = kPID;
    	SetStreamLevel(level, task);
    	continue;
		}
  }
}

//_____________________________________________________________________________
void AliTRDReconstructor::SetStreamLevel(Int_t level, AliTRDReconstructorTask task){
	//
	// Set the Stream Level for one of the tasks Clusterizer, Tracker or PID
	//
	TString taskname;
	switch(task){
		case kClusterizer: taskname = "Clusterizer";
											 break;
		case kTracker: taskname = "Tracker";
									 break;
		case kPID: taskname = "PID";
							 break;
	}
 	AliInfo(Form("Setting Stream Level for Task %s to %d", taskname.Data(),level));
	fStreamLevel[(Int_t)task] = level;
}
