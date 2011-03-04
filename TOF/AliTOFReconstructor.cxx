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
#include "AliTOFT0maker.h"
#include "AliTOFReconstructor.h"

class TTree;

ClassImp(AliTOFReconstructor)

 //____________________________________________________________________
AliTOFReconstructor::AliTOFReconstructor() 
  : AliReconstructor(),
    fTOFcalib(0)/*,
		  fTOFT0maker(0)*/
{
//
// ctor
//
  
  //Retrieving the TOF calibration info  
  fTOFcalib = new AliTOFcalib();
  fTOFcalib->Init();

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

}

//_____________________________________________________________________________
void AliTOFReconstructor::Reconstruct(AliRawReader *rawReader,
                                      TTree *clustersTree) const
{
  //
  // reconstruct clusters from Raw Data
  //

  TString optionString = GetOption();

  // use V1 cluster finder if selected
  if (optionString.Contains("ClusterizerV1")) {
    static AliTOFClusterFinderV1 tofClus(fTOFcalib);

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
  }
  else {
    static AliTOFClusterFinder tofClus(fTOFcalib);
    
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
  }

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
    static AliTOFClusterFinderV1 tofClus(fTOFcalib);

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
  }
  else {
    static AliTOFClusterFinder tofClus(fTOFcalib);

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
  }

}
//_____________________________________________________________________________
  void AliTOFReconstructor::ConvertDigits(AliRawReader* reader, TTree* digitsTree) const
{
// reconstruct clusters from digits

  AliDebug(2,Form("Global Event loop mode: Converting Raw Data to a Digits Tree")); 

  TString optionString = GetOption();
  // use V1 cluster finder if selected
  if (optionString.Contains("ClusterizerV1")) {
    static AliTOFClusterFinderV1 tofClus(fTOFcalib);

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
  }
  else {
    static AliTOFClusterFinder tofClus(fTOFcalib);

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

  if (!GetRecoParam()) AliFatal("cannot get TOF RECO params");

  Float_t tofResolution = GetRecoParam()->GetTimeResolution();// TOF time resolution in ps
  AliTOFT0maker *tofT0maker = new AliTOFT0maker(esdPID,fTOFcalib);
  //AliTOFT0maker tofT0maker = AliTOFT0maker(esdPID,fTOFcalib);
  tofT0maker->SetTimeResolution(tofResolution);
  tofT0maker->ComputeT0TOF(event);
  tofT0maker->WriteInESD(event);
  tofT0maker->~AliTOFT0maker();
  delete tofT0maker;

}

//_____________________________________________________________________________
void 
AliTOFReconstructor::FillESD(TTree *, TTree *, AliESDEvent *esdEvent) const
{
  //
  // correct Texp 
  // 
  //

  //  fTOFcalib->CalibrateTExp(esdEvent);
}
