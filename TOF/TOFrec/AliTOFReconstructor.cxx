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

/* $Id: AliTOFReconstructor.cxx 59948 2012-12-12 11:05:59Z fnoferin $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include "TObjArray.h"
#include "TString.h"

#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliRawReader.h"
#include "AliTOFHeader.h"

#include "AliTOFClusterFinder.h"
#include "AliTOFClusterFinderV1.h"
#include "AliTOFcalib.h"
#include "AliTOFtrackerMI.h"
#include "AliTOFtracker.h"
#include "AliTOFtrackerV1.h"
#include "AliTOFtrackerV2.h"
#include "AliTOFT0maker.h"
#include "AliTOFReconstructor.h"
#include "AliTOFTriggerMask.h"
#include "AliTOFTrigger.h"
#include "AliCTPTimeParams.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

class TTree;

ClassImp(AliTOFReconstructor)

Double_t AliTOFReconstructor::fgExtraTolerance = 0;
Int_t AliTOFReconstructor::fgCTPtriggerLatency = -1;

 //____________________________________________________________________
AliTOFReconstructor::AliTOFReconstructor() :
  AliReconstructor(),
  fTOFcalib(0),
  /*fTOFT0maker(0),*/
  fNumberOfTofClusters(0),
  fNumberOfTofTrgPads(0),
  fClusterFinder(0),
  fClusterFinderV1(0)
{
//
// ctor
//
  
  //Retrieving the TOF calibration info  
  fTOFcalib = new AliTOFcalib();
  fTOFcalib->Init();
  fClusterFinder = new AliTOFClusterFinder(fTOFcalib);
  fClusterFinderV1 = new AliTOFClusterFinderV1(fTOFcalib);

  TString optionString = GetOption();
  if (optionString.Contains("DecoderV0")) {
    fClusterFinder->SetDecoderVersion(0);
    fClusterFinderV1->SetDecoderVersion(0);
  }
  else if (optionString.Contains("DecoderV1")) {
    fClusterFinder->SetDecoderVersion(1);
    fClusterFinderV1->SetDecoderVersion(1);
  }
  else {
    fClusterFinder->SetDecoderVersion(2);
    fClusterFinderV1->SetDecoderVersion(2);
  }



#if 0
  fTOFcalib->CreateCalObjects();

  if(!fTOFcalib->ReadParOnlineDelayFromCDB("TOF/Calib",-1)) {AliFatal("Exiting, no CDB object found!!!");exit(0);}  
  if(!fTOFcalib->ReadParOnlineStatusFromCDB("TOF/Calib",-1)) {AliFatal("Exiting, no CDB object found!!!");exit(0);}  

  if(!fTOFcalib->ReadParOfflineFromCDB("TOF/Calib",-1)) {AliFatal("Exiting, no CDB object found!!!");exit(0);}  


  if(!fTOFcalib->ReadDeltaBCOffsetFromCDB("TOF/Calib",-1)) {AliFatal("Exiting, no CDB object found!!!");exit(0);}  
  if(!fTOFcalib->ReadCTPLatencyFromCDB("TOF/Calib",-1)) {AliFatal("Exiting, no CDB object found!!!");exit(0);}  
  if(!fTOFcalib->ReadT0FillFromCDB("TOF/Calib",-1)) {AliFatal("Exiting, no CDB object found!!!");exit(0);}  
  if(!fTOFcalib->ReadRunParamsFromCDB("TOF/Calib",-1)) {AliFatal("Exiting, no CDB object found!!!");exit(0);}  
#endif

}

//_____________________________________________________________________________
AliTOFReconstructor::~AliTOFReconstructor() 
{
//
// dtor
//

  delete fTOFcalib;

  //delete fTOFT0maker;
  fNumberOfTofClusters = 0;
  fNumberOfTofTrgPads = 0;

  delete fClusterFinder;
  delete fClusterFinderV1;
}
void AliTOFReconstructor::GetPidSettings(AliESDpid *esdPID){
  Float_t tofResolution = GetRecoParam()->GetTimeResolution();// TOF time resolution in ps
  AliInfo(Form(" TOF resolution set in PIDResponse to %f ps",tofResolution));
  
  esdPID->GetTOFResponse().SetTimeResolution(tofResolution);

  return;
}

//_____________________________________________________________________________
void AliTOFReconstructor::Reconstruct(AliRawReader *rawReader,
                                      TTree *clustersTree) const
{
  AliInfo(Form("TOF Reconstruct"));
 
  //
  // reconstruct clusters from Raw Data
  //

  TString optionString = GetOption();

  // use V1 cluster finder if selected
  if (optionString.Contains("ClusterizerV1")) {
    /*
    AliTOFClusterFinderV1 tofClus(fTOFcalib);

    // decoder version option
    if (optionString.Contains("DecoderV0")) {
      tofClus.SetDecoderVersion(0);
    }
    else if (optionString.Contains("DecoderV1")) {
      tofClus.SetDecoderVersion(1);
    }
    else {
      tofClus.SetDecoderVersion(2);
    }

    tofClus.Digits2RecPoints(rawReader, clustersTree);
    */

    fClusterFinderV1->Digits2RecPoints(rawReader, clustersTree);    
  }
  else {
    /*
    AliTOFClusterFinder tofClus(fTOFcalib);
      
    // decoder version option
    if (optionString.Contains("DecoderV0")) {
      tofClus.SetDecoderVersion(0);
    }
    else if (optionString.Contains("DecoderV1")) {
      tofClus.SetDecoderVersion(1);
    }
    else {
      tofClus.SetDecoderVersion(2);
    }

    tofClus.Digits2RecPoints(rawReader, clustersTree);

    */

    fClusterFinder->Digits2RecPoints(rawReader, clustersTree);
  }


  if(fgCTPtriggerLatency < 0){ // read from OCDB
    AliCDBManager *man = AliCDBManager::Instance();
    Int_t run = man->GetRun();

    if(run > 256144){
      fgCTPtriggerLatency = 11600; // run-2 value 2016
    }
     else if(run > 244335){
      fgCTPtriggerLatency = 12800; // run-2 value
    }
    else
      fgCTPtriggerLatency = 13600; // run-1 value

    AliInfo(Form("CTP latency used for run %i to select bunch ID = %i",run,fgCTPtriggerLatency));
  }


  AliTOFTrigger::PrepareTOFMapFromRaw(rawReader,fgCTPtriggerLatency); // 13600 +/- 400 is the value to select the richt bunch crossing (in future from OCDB)
}

//_____________________________________________________________________________
void AliTOFReconstructor::Reconstruct(TTree *digitsTree,
                                      TTree *clustersTree) const
{
  //
  // reconstruct clusters from digits
  //

  AliDebug(2,Form("Global Event loop mode: Creating Recpoints from Digits Tree")); 

  TString optionString = GetOption();
  // use V1 cluster finder if selected
  if (optionString.Contains("ClusterizerV1")) {
    /*
    AliTOFClusterFinderV1 tofClus(fTOFcalib);

    // decoder version option
    if (optionString.Contains("DecoderV0")) {
      tofClus.SetDecoderVersion(0);
    }
    else if (optionString.Contains("DecoderV1")) {
      tofClus.SetDecoderVersion(1);
    }
    else {
      tofClus.SetDecoderVersion(2);
    }
    
    tofClus.Digits2RecPoints(digitsTree, clustersTree);
    */
    fClusterFinderV1->Digits2RecPoints(digitsTree, clustersTree);
  }
  else {
    /*
    AliTOFClusterFinder tofClus(fTOFcalib);

    // decoder version option
    if (optionString.Contains("DecoderV0")) {
      tofClus.SetDecoderVersion(0);
    }
    else if (optionString.Contains("DecoderV1")) {
      tofClus.SetDecoderVersion(1);
    }
    else {
      tofClus.SetDecoderVersion(2);
    }
    
    tofClus.Digits2RecPoints(digitsTree, clustersTree);
    */

    fClusterFinder->Digits2RecPoints(digitsTree, clustersTree);

  }
  AliTOFTrigger::PrepareTOFMapFromDigit(digitsTree);

}
//_____________________________________________________________________________
  void AliTOFReconstructor::ConvertDigits(AliRawReader* reader, TTree* digitsTree) const
{
// reconstruct clusters from digits

  AliDebug(2,Form("Global Event loop mode: Converting Raw Data to a Digits Tree")); 

  TString optionString = GetOption();
  // use V1 cluster finder if selected
  if (optionString.Contains("ClusterizerV1")) {
    /*
    AliTOFClusterFinderV1 tofClus(fTOFcalib);

    // decoder version option
    if (optionString.Contains("DecoderV0")) {
      tofClus.SetDecoderVersion(0);
    }
    else if (optionString.Contains("DecoderV1")) {
      tofClus.SetDecoderVersion(1);
    }
    else {
      tofClus.SetDecoderVersion(2);
    }
    
    tofClus.Raw2Digits(reader, digitsTree);
    */

    fClusterFinderV1->Digits2RecPoints(reader, digitsTree);
  }
  else {
    /*
    AliTOFClusterFinder tofClus(fTOFcalib);

    // decoder version option
    if (optionString.Contains("DecoderV0")) {
      tofClus.SetDecoderVersion(0);
    }
    else if (optionString.Contains("DecoderV1")) {
      tofClus.SetDecoderVersion(1);
    }
    else {
      tofClus.SetDecoderVersion(2);
    }
    
    tofClus.Raw2Digits(reader, digitsTree);
    */

    fClusterFinder->Digits2RecPoints(reader, digitsTree);

  }

}

//_____________________________________________________________________________
AliTracker* AliTOFReconstructor::CreateTracker() const
{

  // 
  // create a TOF tracker using 
  // TOF Reco Param collected by STEER
  //

  TString selectedTracker = GetOption();
 
  AliTracker *tracker;
  // use MI tracker if selected
  if (selectedTracker.Contains("TrackerMI")) {
    tracker = new AliTOFtrackerMI();
  }
  // use V1 tracker if selected
  else if (selectedTracker.Contains("TrackerV1")) {
    tracker =  new AliTOFtrackerV1();
  }
  else if (selectedTracker.Contains("TrackerV2")) {
    tracker =  new AliTOFtrackerV2();
  }
  else {
    tracker = new AliTOFtracker();
  }
  return tracker;

}

//_____________________________________________________________________________
void AliTOFReconstructor::FillEventTimeWithTOF(AliESDEvent *event, AliESDpid *esdPID)
{
  //
  // Fill AliESDEvent::fTOFHeader variable
  // It contains the event_time estiamted by the TOF combinatorial algorithm
  //


  // Set here F. Noferini
  AliTOFTriggerMask *mapTrigger = AliTOFTrigger::GetTOFTriggerMap();


  TString optionString = GetOption();
  if (optionString.Contains("ClusterizerV1")) {
    fNumberOfTofClusters=fClusterFinderV1->GetNumberOfTOFclusters();
    fNumberOfTofTrgPads=fClusterFinderV1->GetNumberOfTOFtrgPads();
    AliInfo(Form(" Number of TOF cluster readout = %d ",fNumberOfTofClusters));
    AliInfo(Form(" Number of TOF cluster readout in trigger window = %d ",fNumberOfTofTrgPads));
  } else {
    fNumberOfTofClusters=fClusterFinder->GetNumberOfTOFclusters();
    fNumberOfTofTrgPads=fClusterFinder->GetNumberOfTOFtrgPads();
    AliInfo(Form(" Number of TOF cluster readout = %d ",fNumberOfTofClusters));
    AliInfo(Form(" Number of TOF cluster readout in trigger window = %d ",fNumberOfTofTrgPads));
  }

  if (!GetRecoParam()) AliFatal("cannot get TOF RECO params");

  AliTOFT0maker *tofT0maker = new AliTOFT0maker(esdPID);
  //Float_t tofResolution = GetRecoParam()->GetTimeResolution();// TOF time resolution in ps
  //tofT0maker->SetTimeResolution(tofResolution); // setting performed now in GetPidSetting in AliReconstructor
  tofT0maker->ComputeT0TOF(event);
  tofT0maker->WriteInESD(event);
  tofT0maker->~AliTOFT0maker();
  delete tofT0maker;

  esdPID->SetTOFResponse(event,(AliESDpid::EStartTimeType_t)GetRecoParam()->GetStartTimeType());


  event->GetTOFHeader()->SetNumberOfTOFclusters(fNumberOfTofClusters);
  event->GetTOFHeader()->SetNumberOfTOFtrgPads(fNumberOfTofTrgPads);
  if(mapTrigger)  event->GetTOFHeader()->SetTriggerMask(mapTrigger);
  AliInfo(Form(" Number of readout cluster in trigger window = %d ; number of trgPads from Trigger map = %d",
	       event->GetTOFHeader()->GetNumberOfTOFtrgPads(),
	       event->GetTOFHeader()->GetNumberOfTOFmaxipad()));

  fClusterFinderV1->ResetDigits();
  fClusterFinderV1->ResetRecpoint();
  fClusterFinder->ResetRecpoint();
  fClusterFinderV1->Clear();
  fClusterFinder->Clear();

}

//_____________________________________________________________________________
void 
AliTOFReconstructor::FillESD(TTree *, TTree *, AliESDEvent * /*esdEvent*/) const
{
  //
  // correct Texp 
  // 
  //

  //  fTOFcalib->CalibrateTExp(esdEvent);
}
