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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Base class for the calibration components using 
//  as input TPCseeds and ESDs
//  Event loop outside of the component
//
//
// Base functionality to be implemeneted by component 
/* 
   //In some cases only one of this function to be implemented
   virtual void     Process(AliESDEvent *event)
   virtual void     Process(AliTPCseed *track)
   //
   virtual Long64_t Merge(TCollection *li);
   virtual void     Analyze()
   void             Terminate();
*/
// Functionality provided by base class for Algorith debuging:
//  TTreeSRedirector * cstream =  GetDebugStreamer() - get debug streamer which can be use for numerical debugging
//                      



//  marian.ivanov@cern.ch
// 
#include "AliTPCcalibBase.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include "TTimeStamp.h"
#include "AliESDEvent.h"


ClassImp(AliTPCcalibBase)

AliTPCcalibBase::AliTPCcalibBase():
    TNamed(),
    fDebugStreamer(0),
    fStreamLevel(0),   
    fRun(0),                  //!  current Run number
    fEvent(0),                //!  current Event number
    fTime(0),                 //!  current Time
    fTrigger(0),              //! current trigger type
    fMagF(0),                 //! current magnetic field
    fTriggerMaskReject(-1),   //trigger mask - reject trigger
    fTriggerMaskAccept(-1),   //trigger mask - accept trigger
    fHasLaser(kFALSE),                    //flag the laser is overlayed with given event 
    fRejectLaser(kTRUE),                 //flag- reject laser
    fDebugLevel(0)
{
  //
  // Constructor
  //
}

AliTPCcalibBase::AliTPCcalibBase(const AliTPCcalibBase&calib):
  TNamed(calib),
  fDebugStreamer(0),
  fStreamLevel(calib.fStreamLevel),
  fRun(0),                  //!  current Run number
  fEvent(0),                //!  current Event number
  fTime(0),                 //!  current Time
  fTrigger(0),              //! current trigger type
  fMagF(0),                 //! current magnetic field
  fTriggerMaskReject(calib.fTriggerMaskReject),   //trigger mask - reject trigger
  fTriggerMaskAccept(calib.fTriggerMaskAccept),   //trigger mask - accept trigger
  fHasLaser(calib.fHasLaser),                    //flag the laser is overlayed with given event
  fRejectLaser(calib.fRejectLaser),                 //flag- reject laser
  fDebugLevel(calib.fDebugLevel)
{
  //
  // copy constructor
  //
}

AliTPCcalibBase &AliTPCcalibBase::operator=(const AliTPCcalibBase&calib){
  //
  //
  //
  ((TNamed *)this)->operator=(calib);
  fDebugStreamer=0;
  fStreamLevel=calib.fStreamLevel;
  fDebugLevel=calib.fDebugLevel;
  return *this;
}


AliTPCcalibBase::~AliTPCcalibBase() {
  //
  // destructor
  //
  if (fDebugLevel>0) printf("AliTPCcalibBase::~AliTPCcalibBase\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer=0;
}

void  AliTPCcalibBase::Terminate(){
  //
  //
  //
  if (fDebugLevel>0) printf("AliTPCcalibBase::Terminate\n");
  if (fDebugStreamer) delete fDebugStreamer;
  fDebugStreamer = 0;
  return;
}

TTreeSRedirector *AliTPCcalibBase::GetDebugStreamer(){
  //
  // Get Debug streamer
  // In case debug streamer not yet initialized and StreamLevel>0 create new one
  //
  if (fStreamLevel==0) return 0;
  if (fDebugStreamer) return fDebugStreamer;
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ",""); 
  fDebugStreamer = new TTreeSRedirector(dsName.Data());
  return fDebugStreamer;
}


void    AliTPCcalibBase::UpdateEventInfo(AliESDEvent * event){
  //
  //
  //
  fRun     = event->GetRunNumber();
  fEvent   = event->GetEventNumberInFile();
  fTime    = event->GetTimeStamp();
  fTrigger = event->GetTriggerMask();
  fMagF    = event->GetMagneticField();
  fHasLaser = HasLaser(event); 
}


Bool_t AliTPCcalibBase::HasLaser(AliESDEvent *event){
  //
  //
  //
  // Thresholds more than 8 tracks with small dip angle
  
  const Int_t kMinLaserTracks = 8;
  const Float_t kThrLaser       = 0.3;
  const Float_t kLaserTgl       = 0.01;

  Int_t ntracks = event->GetNumberOfTracks();
  if (ntracks<kMinLaserTracks) return kFALSE;
  Float_t nlaser=0;
  Float_t nall=0;
  for (Int_t i=0;i<ntracks;++i) {
    AliESDtrack *track=event->GetTrack(i);
    if (!track) continue;
    if (track->GetTPCNcls()<=0) continue; 
    nall++;
    if (TMath::Abs(track->GetTgl())<kLaserTgl) nlaser++;
  }
  if (nlaser>kMinLaserTracks) return kTRUE;
  if (nall>0 && nlaser/nall>kThrLaser) return kTRUE;
  return kFALSE;
}



Bool_t AliTPCcalibBase::AcceptTrigger(){
  //
  // Apply trigger mask - Don't do calibration for non proper triggers
  // 
  if (fTriggerMaskReject==(Int_t)fTrigger) return kFALSE;
  if (fTriggerMaskAccept>0 && fTriggerMaskAccept!=(Int_t)fTrigger) return kFALSE;
  if (fHasLaser && fRejectLaser) return kFALSE;
  return kTRUE;
}


void AliTPCcalibBase::RegisterDebugOutput(const char *path){
  //
  // store  - copy debug output to the destination position
  // currently ONLY for local copy
  if (fDebugLevel>0) printf("AliTPCcalibBase::RegisterDebugOutput(%s)\n",path);
  if (fStreamLevel==0) return;
  TString dsName;
  dsName=GetName();
  dsName+="Debug.root";
  dsName.ReplaceAll(" ",""); 
  TString dsName2=path;
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+=gSystem->HostName();
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+="/";
  TTimeStamp s;
  dsName2+=Int_t(s.GetNanoSec());
  dsName2+="/";
  gSystem->MakeDirectory(dsName2.Data());
  dsName2+=dsName;
  AliInfo(Form("copy %s\t%s\n",dsName.Data(),dsName2.Data()));
  printf("copy %s\t%s\n",dsName.Data(),dsName2.Data());
  TFile::Cp(dsName.Data(),dsName2.Data());
}
