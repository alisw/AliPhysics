// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Felix Rettig, Stefan Kirsch                           *
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

/// @file   AliHLTTRDGlobalMonitorComponent.cxx
/// @author Felix Rettig, Stefan Kirsch
/// @date   2011-08-02
/// @brief  A processing component for TRD tracking/trigger data on CN-level
///

#include "AliHLTTRDGlobalMonitorComponent.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTRDDefinitions.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TObjString.h"
#include "TObjArray.h"

#define TRDMODULES 18
#define TRDMAXDDLSIZE 10000

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDGlobalMonitorComponent)

AliHLTTRDGlobalMonitorComponent::AliHLTTRDGlobalMonitorComponent()
  :
    fHistArray(NULL)
  , fHistTrackletY(NULL)
  , fHistTrackletDy(NULL)
  , fHistTrackletZ(NULL)
  , fHistTrackletPID(NULL)
  , fHistTrackletYDy(NULL)
  , fHistTrackletHC(NULL)
  , fHistTrackletBadY(NULL)
  , fHistTrackletBadPID(NULL)
  , fHistFirstTrackletTime(NULL)
  , fHistLastTrackletTime(NULL)
  , fHistTmuTime(NULL)
  , fHistSmuTime(NULL)
  , fHistTrackPt(NULL)
  , fHistTrackPID(NULL)
  , fHistTrackStack(NULL)
  , fHistTrackletsTrack(NULL)
  , fHistTrackletsTrackHpt(NULL)
  , fHistTriggerContribs(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTRDGlobalMonitorComponent::~AliHLTTRDGlobalMonitorComponent()
{
  // see header file for class documentation
}

const char* AliHLTTRDGlobalMonitorComponent::GetComponentID()
{
  // see header file for class documentation
  return "TRDGlobalMonitorComponent";
}

void AliHLTTRDGlobalMonitorComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  // list.push_back(AliHLTTRDDefinitions::fgkSimpleIntegerDataType);
  list.push_back(kAliHLTDataTypeTObject | kAliHLTDataOriginTRD);
}

AliHLTComponentDataType AliHLTTRDGlobalMonitorComponent::GetOutputDataType()
{
  // see header file for class documentation
  return (kAliHLTDataTypeTObjArray | kAliHLTDataOriginTRD);
}

void AliHLTTRDGlobalMonitorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase = 10000000;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTRDGlobalMonitorComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTRDGlobalMonitorComponent;
}

int AliHLTTRDGlobalMonitorComponent::DoInit( int argc, const char** argv )
{
  fHistArray = new TObjArray(25);
  if(!fHistArray)
    return -ENOMEM;
  fHistArray->SetOwner(kTRUE);

  fHistTrackletY = new TH1I("hist_tracklets_y", "Y-Position of online tracklets", 256, -4096, 4096);
  fHistArray->AddLast(fHistTrackletY);

  fHistTrackletDy = new TH1I("hist_tracklets_dy", "Deflections of online tracklets", 128, -64, 64);
  fHistArray->AddLast(fHistTrackletDy);

  fHistTrackletZ = new TH1I("hist_tracklets_z", "Z-Position of online tracklets", 12, 0, 12);
  fHistArray->AddLast(fHistTrackletZ);

  fHistTrackletPID = new TH1I("hist_tracklets_pid", "PID of online tracklets", 256, 0, 256);
  fHistArray->AddLast(fHistTrackletPID);

  fHistTrackletYDy = new TH2I("hist_tracklets_y_dy", "Tracklet deflection vs. tracklet position", 256, -4096, 4096, 64, -64, 64);
  fHistArray->AddLast(fHistTrackletYDy);

  fHistTrackletHC = new TH2I("hist_tracklets_hc", "Number of online tracklets by HC", 18, 0, 18, 60, 0, 60);
  fHistArray->AddLast(fHistTrackletHC);

  fHistTrackletBadY = new TH2I("hist_tracklets_bad_y", "Number of online tracklets with bad y-position by stack", 18, 0, 18, 5, 0, 5);
  fHistArray->AddLast(fHistTrackletBadY);

  fHistTrackletBadPID = new TH2I("hist_tracklets_bad_pid", "Number of online tracklets with bad PID value by stack", 18, 0, 18, 5, 0, 5);
  fHistArray->AddLast(fHistTrackletBadPID);

  fHistFirstTrackletTime = new TH1F("hist_first_tracklet_time", "Arrival time of first tracklet", 160, 0., 8.);
  fHistArray->AddLast(fHistFirstTrackletTime);

  fHistLastTrackletTime = new TH1F("hist_last_tracklet_time", "Arrival time of last tracklet", 160, 0., 8.);
  fHistArray->AddLast(fHistLastTrackletTime);

  fHistTmuTime = new TH1F("hist_tmu_time", "Tracking done time TMU-level", 160, 0., 8.);
  fHistArray->AddLast(fHistTmuTime);

  fHistSmuTime = new TH1F("hist_smu_time", "Tracking done time SMU-level", 160, 0., 8.);
  fHistArray->AddLast(fHistSmuTime);

  fHistTrackPt = new TH1F("hist_tracks_pt", "Transverse momentum of GTU  tracks", 100, 0., 20.);
  fHistArray->AddLast(fHistTrackPt);

  fHistTrackPID = new TH1I("hist_tracks_pid", "PID value of GTU  tracks", 256, 0, 256);
  fHistArray->AddLast(fHistTrackPID);

  fHistTrackStack = new TH2I("hist_tracks_stack", "Number of GTU tracks by stack", 18, 0, 18, 5, 0, 5);
  fHistArray->AddLast(fHistTrackStack);

  fHistTrackletsTrack = new TH1I("hist_tracklets_track", "Tracklets per GTU track", 7, 0, 7);
  fHistArray->AddLast(fHistTrackletsTrack);

  fHistTrackletsTrackHpt = new TH1I("hist_tracklets_track_hpt", "Tracklets per high-pt GTU track", 7, 0, 7);
  fHistArray->AddLast(fHistTrackletsTrackHpt);

  fHistTriggerContribs = new TH2I("hist_trigger_contribs", "Trigger contributions by segment", 18, 0, 18, 12, 0, 12);
  fHistArray->AddLast(fHistTriggerContribs);

  return 0;
}

int AliHLTTRDGlobalMonitorComponent::DoDeinit()
{
  if(fHistArray){
    fHistArray->Clear();
    delete fHistArray;
    fHistArray=NULL;
  }
  return 0;
}

int AliHLTTRDGlobalMonitorComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
 
  int iResult=0;
  //  Int_t totaleventsize=0;

  /*
  //Loop over integers
  for (const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(AliHLTTRDDefinitions::fgkSimpleIntegerDataType);
       pBlock!=NULL && iResult>=0;
       pBlock=GetNextInputBlock()) {
    // extract DDL id from specification
    int ddlnum=-1;
    for (unsigned pos=0; pos<8*sizeof(AliHLTUInt32_t); pos++) {
      if (pBlock->fSpecification & (0x1<<pos)) {
	if (ddlnum>=0) {
	  HLTWarning("Can not uniquely identify DDL number from specification, skipping data block %s 0x%08x",
		     DataType2Text(pBlock->fDataType).c_str(),
		     pBlock->fSpecification);
	  ddlnum=-1;
	  break;
	}
	ddlnum=pos;
      }
    }
    if (ddlnum<0 || ddlnum >= TRDMODULES) continue;

    if ( pBlock->fSize < sizeof(Int_t) ) {
      HLTWarning("block size (%d) smaller than expected size (%d)!", pBlock->fSize, sizeof(Int_t) );
      continue;
    } 

    Int_t* eventsize = (Int_t *) pBlock->fPtr;
    if(*eventsize > TRDMAXDDLSIZE)
      HLTWarning("eventsize (%d) larger than histogram range (%d)", *eventsize, TRDMAXDDLSIZE);
    // Fill histo of ddl
    TH1F* h = dynamic_cast<TH1F*>(fHistoArray->At(ddlnum));
    h->Fill((Double_t)*eventsize);

    totaleventsize += *eventsize;
  }
  */

  // reset the high-level histograms, the accumulation over events is done in the lower level
  fHistTrackletY->Reset();
  fHistTrackletDy->Reset();
  fHistTrackletZ->Reset();
  fHistTrackletPID->Reset();
  fHistTrackletYDy->Reset();
  fHistTrackletHC->Reset();
  fHistTrackletBadY->Reset();
  fHistTrackletBadPID->Reset();
  fHistFirstTrackletTime->Reset();
  fHistLastTrackletTime->Reset();
  fHistTmuTime->Reset();
  fHistSmuTime->Reset();
  fHistTrackPt->Reset();
  fHistTrackPID->Reset();
  fHistTrackStack->Reset();
  fHistTrackletsTrack->Reset();
  fHistTrackletsTrackHpt->Reset();
  fHistTriggerContribs->Reset();

  // loop over TObject-based input data
  for (const TObject* obj = 
	 GetFirstInputObject(kAliHLTDataTypeTObject | kAliHLTDataOriginTRD, "TObjArray");
       obj!=NULL && iResult>=0;
       obj=GetNextInputObject()) {
    // extract DDL id from specification
    int ddlnum=-1;
    AliHLTUInt32_t specification = GetSpecification(obj);
    for (unsigned pos=0; pos<8*sizeof(AliHLTUInt32_t); pos++) {
      if (specification & (0x1<<pos)) {
        if (ddlnum>=0) {
          HLTWarning("Can not uniquely identify DDL number from specification, skipping data block %s 0x%08x",
        	     DataType2Text(GetDataType(obj)).c_str(),
        	     specification);
          ddlnum=-1;
          break;
        }
	      ddlnum=pos;
      }
    }
    if (ddlnum<0 || ddlnum >= TRDMODULES) continue;
    
    // input object is a TObjArray containing the actual data containers
    const TObjArray* histArray = dynamic_cast<const TObjArray*>(obj);
    if (!histArray){
      HLTWarning("Received object was not a TObjAarray");
      continue;
    }

    // extract data containers and process the data

    TH1I *hist1 = (TH1I*)histArray->FindObject("hist_tracklets_y");
    if (hist1){
      fHistTrackletY->Add(hist1);
    } else
      HLTWarning("Tracklet y-position histogram not found!");

    hist1 = (TH1I*)histArray->FindObject("hist_tracklets_dy");
    if (hist1){
      fHistTrackletDy->Add(hist1);
    } else
      HLTWarning("Tracklet deflection histogram not found!");

    hist1 = (TH1I*)histArray->FindObject("hist_tracklets_z");
    if (hist1){
      fHistTrackletZ->Add(hist1);
    } else
      HLTWarning("Tracklet z-position histogram not found!");

    hist1 = (TH1I*)histArray->FindObject("hist_tracklets_pid");
    if (hist1){
      fHistTrackletPID->Add(hist1);
    } else
      HLTWarning("Tracklet PID histogram not found!");

    TH2I *hist2 = (TH2I*)histArray->FindObject("hist_tracklets_y_dy");
    if (hist2){
      fHistTrackletYDy->Add(hist2);
    } else
      HLTWarning("Tracklet y vs. dy histogram not found!");

    hist2 = (TH2I*)histArray->FindObject("hist_tracklets_hc");
    if (hist2){
      fHistTrackletHC->Add(hist2);
    } else
      HLTWarning("Tracklet number histogram not found!");

    hist2 = (TH2I*)histArray->FindObject("hist_tracklets_bad_y");
    if (hist2){
      fHistTrackletBadY->Add(hist2);
    } else
      HLTWarning("Tracklet number with bad y-position histogram not found!");

    hist2 = (TH2I*)histArray->FindObject("hist_tracklets_bad_pid");
    if (hist2){
      fHistTrackletBadPID->Add(hist2);
    } else
      HLTWarning("Tracklet number with bad y-position histogram not found!");

    TH1F *hist3 = (TH1F*)histArray->FindObject("hist_first_tracklet_time");
    if (hist3){
      fHistFirstTrackletTime->Add(hist3);
    } else
      HLTWarning("First tracklet arrival  histogram not found!");

    hist3 = (TH1F*)histArray->FindObject("hist_last_tracklet_time");
    if (hist3){
      fHistLastTrackletTime->Add(hist3);
    } else
      HLTWarning("Last tracklet arrival time histogram not found!");

    hist3 = (TH1F*)histArray->FindObject("hist_tmu_time");
    if (hist3){
      fHistTmuTime->Add(hist3);
    } else
      HLTWarning("TMU-level tracking done time histogram not found!");

    hist3 = (TH1F*)histArray->FindObject("hist_smu_time");
    if (hist3){
      fHistSmuTime->Add(hist3);
    } else
      HLTWarning("SMU-level tracking done time histogram not found!");

    hist3 = (TH1F*)histArray->FindObject("hist_tracks_pt");
    if (hist3){
      fHistTrackPt->Add(hist3);
    } else
      HLTWarning("Track pt histogram not found!");

    hist1 = (TH1I*)histArray->FindObject("hist_tracks_pid");
    if (hist1){
      fHistTrackPID->Add(hist1);
    } else
      HLTWarning("Tracklet PID histogram not found!");

    hist2 = (TH2I*)histArray->FindObject("hist_tracks_stack");
    if (hist2){
      fHistTrackStack->Add(hist2);
    } else
      HLTWarning("Track number histogram not found!");

    hist1 = (TH1I*)histArray->FindObject("hist_tracklets_track");
    if (hist1){
      fHistTrackletsTrack->Add(hist1);
    } else
      HLTWarning("Tracklets per GTU track histogram not found!");

    hist1 = (TH1I*)histArray->FindObject("hist_tracklets_track_hpt");
    if (hist1){
      fHistTrackletsTrackHpt->Add(hist1);
    } else
      HLTWarning("Tracklets per high-pt GTU track histogram not found!");

    hist2 = (TH2I*)histArray->FindObject("hist_trigger_contribs");
    if (hist2){
      fHistTriggerContribs->Add(hist2);
    } else
      HLTWarning("Trigger contribution histogram not found!");

  }

  iResult = PushBack(fHistArray, 
		     (kAliHLTDataTypeTObjArray | kAliHLTDataOriginTRD), 
		     //Specification: all 18 bits (links) set
		     (0xFFFFFFFF >> (32 - TRDMODULES) )
		     );
		       
  return iResult;
}

int AliHLTTRDGlobalMonitorComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;

  return iResult;
}

int AliHLTTRDGlobalMonitorComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // see header file for class documentation
  int iResult=0;

  return iResult;
}
