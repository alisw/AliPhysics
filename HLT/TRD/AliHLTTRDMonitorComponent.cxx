// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTRDMonitorComponent.cxx
/// @author Felix Rettig, Stefan Kirsch
/// @date   2011-08-02
/// @brief  A processing component for TRD tracking/trigger data on FEP-level
/// @ingroup alihlt_trd_components

#include "AliHLTTRDMonitorComponent.h"
#include "AliLog.h"
#include "AliHLTDataTypes.h"
#include "AliRawReaderMemory.h"
#include "AliHLTTRDDefinitions.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDrawStream.h"
#include "AliTRDtrackletWord.h"
#include "AliESDTrdTracklet.h"
#include "AliESDTrdTrack.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1F.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDMonitorComponent)

AliHLTTRDMonitorComponent::AliHLTTRDMonitorComponent()
  : AliHLTProcessor()
  , fTrackletArray(NULL)
  , fGtuTrackArray(NULL)
  , fRawReaderMem(NULL)
  , fDigitsManagerTrd(NULL)
  , fRawReaderTrd(NULL)
  , fHistArray(NULL)
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
  , fHistTrackStack(NULL)
  , fHistTrackletsTrack(NULL)
  , fHistTrackletsTrackHpt(NULL)
  , fHistTrackPID(NULL)
  , fHistTriggerContribs(NULL)
{
  // constructor
}

AliHLTTRDMonitorComponent::~AliHLTTRDMonitorComponent()
{
  // destructor
}

const char* AliHLTTRDMonitorComponent::GetComponentID()
{ 
  // component property: id
  return "TRDMonitorComponent";
}

void AliHLTTRDMonitorComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // component property: list of input data types
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTRD);
}

AliHLTComponentDataType AliHLTTRDMonitorComponent::GetOutputDataType()
{
  // component property: output data type
  return kAliHLTMultipleDataType;
}

int AliHLTTRDMonitorComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)

{
  // see header file for class documentation
  tgtList.clear();
  //  tgtList.push_back(AliHLTTRDDefinitions::fgkSimpleIntegerDataType);
  tgtList.push_back(kAliHLTDataTypeTObject | kAliHLTDataOriginTRD);
  return tgtList.size();
}

void AliHLTTRDMonitorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // component property: output size estimator
  constBase = 5000000;
  inputMultiplier = 0;
}

void AliHLTTRDMonitorComponent::GetOCDBObjectDescription( TMap* const /*targetMap*/)
{
  // Get a list of OCDB object description.
  // The list of objects is provided in a TMap
  // - key: complete OCDB path, e.g. GRP/GRP/Data
  // - value: short description why the object is needed
  // Key and value objects created inside this class go into ownership of
  // target TMap.
  //if (!targetMap) return;
  //targetMap->Add(new TObjString("HLT/ConfigSample/SampleRawAnalysis"),
  //               new TObjString("configuration object"));
}

AliHLTComponent* AliHLTTRDMonitorComponent::Spawn()
{
  // Spawn function, return new class instance
  return new AliHLTTRDMonitorComponent;
}

int AliHLTTRDMonitorComponent::DoInit( int /*argc*/, const char** /*argv*/ )
{

  int iResult=0;

  // init stage 1: default values for all data members

  // init stage 2: read configuration object
  // ScanConfigurationArgument() needs to be implemented
  // HLT component configuration objects are located in the HLT/ConfigDET
  // of the OCDB. TObjString configuration objects can be generated with
  // the macro HLT/exa/makeComponentConfigurationObject.C, e.g.
  // aliroot -b -q -l $ALICE_ROOT/HLT/exa/makeComponentConfigurationObject.C'("HLT/ConfigSample/SampleRawAnalysis", "")'
  //TString cdbPath="HLT/ConfigSample/";
  //cdbPath+=GetComponentID();
  //iResult=ConfigureFromCDBTObjString(cdbPath);

  // init stage 3: read the component arguments
  //if (iResult>=0) {
  //  iResult=ConfigureFromArgumentString(argc, argv);
  //}

  // implement the component initialization
  // Matthias 2011-08-24: this has to go into ScanConfigurationArgument
  do {

    fRawReaderMem = new AliRawReaderMemory; 
    if (!fRawReaderMem) {
      iResult=-ENOMEM;
      break;
    }

    fTrackletArray = new TClonesArray("AliTRDtrackletWord", 500);
    if (!fTrackletArray) {
      iResult=-ENOMEM;
      break;
    }

    fGtuTrackArray = new TClonesArray("AliESDTrdTrack", 50);
    if (!fGtuTrackArray){
      iResult=-ENOMEM;
      break;
    }

    fDigitsManagerTrd = new AliTRDdigitsManager();
    if (!fDigitsManagerTrd) {
      iResult=-ENOMEM;
      break;
    }
    fDigitsManagerTrd->CreateArrays();

    fRawReaderTrd = new AliTRDrawStream(fRawReaderMem);
    if (!fRawReaderTrd) {
      iResult=-ENOMEM;
      break;
    }

    fRawReaderTrd->SetDigitsManager(fDigitsManagerTrd);
    // fRawReaderTrd->SetDigitsManager(NULL);  // FIXME may be used to improve performance, needs a fix to be committed by TRD
    fRawReaderTrd->SetTrackletArray(fTrackletArray);
    fRawReaderTrd->SetTrackArray(fGtuTrackArray);
  
    // Disable raw reader error messages that could flood HLT logbook
    fRawReaderTrd->SetErrorDebugLevel( AliTRDrawStream::kLinkMonitor, 1 );

    fHistArray = new TObjArray(25);
    if(!fHistArray){
      return -ENOMEM;
      break;
    }
    fHistArray->SetOwner(kTRUE);

    fHistTrackletY = new TH1I("hist_tracklets_y", "Y-Position of online tracklets", 256, -4096, 4096);
    if (!fHistTrackletY){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletY);

    fHistTrackletDy = new TH1I("hist_tracklets_dy", "Deflection of online tracklets", 128, -64, 64);
    if (!fHistTrackletDy){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletDy);

    fHistTrackletZ = new TH1I("hist_tracklets_z", "Z-Position of online tracklets", 12, 0, 12);
    if (!fHistTrackletZ){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletZ);

    fHistTrackletPID = new TH1I("hist_tracklets_pid", "PID of online tracklets", 256, 0, 256);
    if (!fHistTrackletPID){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletPID);

    fHistTrackletYDy = new TH2I("hist_tracklets_y_dy", "Tracklet deflection vs. tracklet position", 256, -4096, 4096, 64, -64, 64);
    if (!fHistTrackletYDy){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletYDy);

    fHistTrackletHC = new TH2I("hist_tracklets_hc", "Number of online tracklets by HC", 18, 0, 18, 60, 0, 60);
    if (!fHistTrackletHC){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletHC);

    fHistTrackletBadY = new TH2I("hist_tracklets_bad_y", "Number of online tracklets with bad y-position by HC", 18, 0, 18, 5, 0, 5);
    if (!fHistTrackletBadY){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletBadY);

    fHistTrackletBadPID = new TH2I("hist_tracklets_bad_pid", "Number of online tracklets with bad y-position by HC", 18, 0, 18, 5, 0, 5);
    if (!fHistTrackletBadPID){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletBadPID);

    fHistFirstTrackletTime = new TH1F("hist_first_tracklet_time", "Arrival time of first tracklet", 160, 0., 8.);
    if (!fHistFirstTrackletTime){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistFirstTrackletTime);

    fHistLastTrackletTime = new TH1F("hist_last_tracklet_time", "Arrival time of last tracklet", 160, 0., 8.);
    if (!fHistLastTrackletTime){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistLastTrackletTime);

    fHistTmuTime = new TH1F("hist_tmu_time", "Tracking done time at TMU-level", 160, 0., 8.);
    if (!fHistTmuTime){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTmuTime);

    fHistSmuTime = new TH1F("hist_smu_time", "Tracking done time at SMU-level", 160, 0., 8.);
    if (!fHistSmuTime){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistSmuTime);

    fHistTrackPt = new TH1F("hist_tracks_pt", "Transverse momentum of GTU tracks", 100, 0., 20.);
    if (!fHistTrackPt){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackPt);

    fHistTrackPID = new TH1I("hist_tracks_pid", "PID of online tracks", 256, 0, 256);
    if (!fHistTrackPID){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackPID);

    fHistTrackStack = new TH2I("hist_tracks_stack", "Number of GTU tracks by stack", 18, 0, 18, 5, 0, 5);
    if (!fHistTrackStack){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackStack);

    fHistTrackletsTrack = new TH1I("hist_tracklets_track", "Tracklets per GTU track", 7, 0, 7);
    if (!fHistTrackletsTrack){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletsTrack);

    fHistTrackletsTrackHpt = new TH1I("hist_tracklets_track_hpt", "Tracklets per high-pt GTU track", 7, 0, 7);
    if (!fHistTrackletsTrackHpt){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTrackletsTrackHpt);

    fHistTriggerContribs = new TH2I("hist_trigger_contribs", "Trigger contributions by segment", 18, 0, 18, 12, 0, 12);
    if (!fHistTriggerContribs){
      return -ENOMEM;
      break;
    }
    fHistArray->AddLast(fHistTriggerContribs);

  } while (0);

  if (iResult<0) {
    // implement cleanup
    if (fRawReaderMem) delete fRawReaderMem;
    fRawReaderMem=NULL;

    if (fRawReaderTrd) delete fRawReaderTrd;
    fRawReaderTrd=NULL;

    if (fTrackletArray) delete fTrackletArray;
    fTrackletArray = NULL;

    if (fGtuTrackArray) delete fGtuTrackArray;
    fGtuTrackArray = NULL;

    if (fHistArray) delete fHistArray;
    fHistArray = NULL;

    fHistTrackletY = NULL;
    fHistTrackletDy = NULL;
    fHistTrackletZ = NULL;
    fHistTrackletPID = NULL;
    fHistTrackletYDy = NULL;
    fHistTrackletHC = NULL;
    fHistTrackletBadY = NULL;
    fHistTrackletBadPID = NULL;
    fHistTrackPt = NULL;
    fHistTrackStack = NULL;
    fHistFirstTrackletTime = NULL;
    fHistLastTrackletTime = NULL;
    fHistTmuTime = NULL;
    fHistSmuTime = NULL;
    fHistTrackPt = NULL;
    fHistTrackStack = NULL;
    fHistTrackletsTrack = NULL;
    fHistTrackletsTrackHpt = NULL;
    fHistTrackPID = NULL;
    fHistTriggerContribs = NULL;
  }

  return iResult;
}

int AliHLTTRDMonitorComponent::ScanConfigurationArgument(int /*argc*/, const char** /*argv*/)
{
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
  // The AliHLTComponent base class implements a parsing loop for argument strings and
  // arrays of strings which is invoked by ConfigureFromArgumentString/ConfigureFromCDBTObjString
  // The component needs to implement ScanConfigurationArgument in order to decode the arguments.

  /*
  int i=0;
  TString argument=argv[i];

  if (argument.IsNull()) return 0;

  // -mandatory1 arg
  if (argument.CompareTo("-mandatory1")==0) {
    if (++i>=argc) return -EINVAL;
    HLTInfo("got \'-mandatory1\' argument: %s", argv[i]);
    return 2; // keyword + 1 argument
  }

  // -optional1 arg
  if (argument.CompareTo("-optional1")==0) {
    if (++i>=argc) return -EINVAL;
    HLTInfo("got \'-optional1\' argument: %s", argv[i]);
    return 2; // keyword + 1 argument
  }

  // -verbose
  if (argument.CompareTo("-verbose")==0) {
    fVerbosity=1;
    return 1; // only keyword
  }
  */
  return 0;
}

int AliHLTTRDMonitorComponent::DoDeinit()
{
  // component cleanup, delete all instances of helper classes here
  if (fRawReaderMem) delete fRawReaderMem;
  fRawReaderMem=NULL;

  if (fRawReaderTrd) delete fRawReaderTrd;
  fRawReaderTrd=NULL;

  if (fTrackletArray) delete fTrackletArray;
  fTrackletArray = NULL;

  if (fGtuTrackArray) delete fGtuTrackArray;
  fGtuTrackArray = NULL;
  
  if (fHistArray) delete fHistArray;
  fHistArray = NULL;

  fHistTrackletY = NULL;
  fHistTrackletDy = NULL;
  fHistTrackletY = NULL;
  fHistTrackletDy = NULL;
  fHistTrackletZ = NULL;
  fHistTrackletPID = NULL;
  fHistTrackletHC = NULL;
  fHistTrackletBadY = NULL;
  fHistTrackletBadPID = NULL;
  fHistFirstTrackletTime = NULL;
  fHistLastTrackletTime = NULL;
  fHistTmuTime = NULL;
  fHistSmuTime = NULL;
  fHistTrackPt = NULL;
  fHistTrackStack = NULL;
  fHistTrackletsTrack = NULL;
  fHistTrackletsTrackHpt = NULL;
  fHistTrackPID = NULL;
  fHistTriggerContribs = NULL;

  return 0;
}

int AliHLTTRDMonitorComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					      AliHLTComponentTriggerData& /*trigData*/)
{
  // event processing function
  int iResult=0;

  // check if this is a data event, there are a couple of special events
  // which should be ignored for normal processing
  if (!IsDataEvent()) return 0;

  // #FIXME: Also take care of SOR, EOR, etc...

  // loop over the raw input data blocks and set up the rawreader
  for (const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTRD);
       pBlock!=NULL && iResult>=0;
       pBlock=GetNextInputBlock()) {
    // extract DDL id from specification
    int ddlnum=-1;
    for (unsigned pos=0; pos<8*sizeof(AliHLTUInt32_t); pos++) {
      if (pBlock->fSpecification & (0x1<<pos)) {
	      if (ddlnum>=0) {
	        // this is just an example, please avoid warnings in every event since those will
	        // saturate the logging system. Consider AliHLTErrorGuard for such cases, e.g.
	        // ALIHLTERRORGUARD(5, "nasty error, first occurence in event %d", event);
	        // this will show the error 5 times and then a summary at the end
	        HLTWarning("Can not uniquely identify DDL number from specification, skipping data block %s 0x%08x",
	      	     DataType2Text(pBlock->fDataType).c_str(),
	      	     pBlock->fSpecification);
	        ddlnum=-1;
	        break;
	      }
	      ddlnum=pos;
      }
    }
    if (ddlnum<0) continue;
    ddlnum += 1024;

    // add data block to rawreader
    if(!fRawReaderMem->AddBuffer((UChar_t*) pBlock->fPtr, pBlock->fSize, ddlnum)){
      HLTError("Could not add buffer of data block  %s, 0x%08x to rawreader",
	       DataType2Text(pBlock->fDataType).c_str(),
	       pBlock->fSpecification);
      continue;
    }

    // read and process TRD tracklet and GTU tracks from event
    fTrackletArray->Clear();
    fGtuTrackArray->Clear();
    fRawReaderTrd->ReadEvent();

    // read and process tracking/trigger flags
    UInt_t iSector = ddlnum-1024;
    UInt_t trgflags = fRawReaderTrd->GetTriggerFlags(iSector);
    //UInt_t done_tmu = (trgflags >> 27) & 0x1f;
    //UInt_t done_smu = (trgflags >> 22) & 0x1f;
    Float_t smu_timing = ((trgflags >> 12) & 0x3ff) * 0.0083333; // tracking finished after L0, in us
    UInt_t ctb_fired =  trgflags & 0xfff;

    /*
      printf("Trigger flags sector %02d: 0x%08x - ctb=0x%03x  tmu-done=0x%x smu-done=0x%x  smu-time: %.2fus\n", 
	   iSector, trgflags,
	   ctb_fired,
	   done_tmu,
	   done_smu,
	   smu_timing
	   );
    */
    fHistSmuTime->Fill(smu_timing);

    for (int iCtb=0; iCtb<12; iCtb++){
      if ((ctb_fired >> iCtb) & 1)
        fHistTriggerContribs->Fill(iSector, iCtb, 1);
    }

    Float_t first_tracklet_timing[5]={0};
    Float_t last_tracklet_timing[5]={0};
    Float_t tmu_timing[5]={0};
    ULong64_t trkflags;

    for (int iStack=0; iStack<5; iStack++){
      trkflags = fRawReaderTrd->GetTrkFlags(iSector, iStack);
      tmu_timing[iStack]=(trkflags & 0x3ff) * 0.02;
      first_tracklet_timing[iStack]=((trkflags >> 20) & 0x3ff) * 0.008;
      last_tracklet_timing[iStack]=((trkflags >> 10) & 0x3ff) * 0.008;
      /*
      printf("  Stack %02d_%d: 0x%016llx - first tracklet: %.2fus  last tracklet: %.2fus  tmu timing: %.2fus\n",
	     iSector, iStack, trkflags, 
	     first_tracklet_timing[iStack], last_tracklet_timing[iStack], tmu_timing[iStack]
	     );
      */
      fHistFirstTrackletTime->Fill(first_tracklet_timing[iStack]);
      fHistLastTrackletTime->Fill(last_tracklet_timing[iStack]);
      fHistTmuTime->Fill(tmu_timing[iStack]);
    }

    for (int iTracklet = 0; iTracklet < fTrackletArray->GetEntriesFast(); iTracklet++) {
      AliTRDtrackletWord *trackletWord = (AliTRDtrackletWord*) ((*fTrackletArray)[iTracklet]);
      AliESDTrdTracklet tracklet(trackletWord->GetTrackletWord(), trackletWord->GetHCId());
      /*
      printf("TRDMSG: TRD tracklet found: 0x%08x  - y=%+5d  dy=%+3d  pid=%3d\n",
	     tracklet.GetTrackletWord(),
	     tracklet.GetBinY(),
	     tracklet.GetBinDy(),
	     tracklet.GetPID());
      */
      // fill some basic histograms right here
      fHistTrackletY->Fill(tracklet.GetBinY());
      fHistTrackletDy->Fill(tracklet.GetBinDy());
      fHistTrackletZ->Fill(tracklet.GetBinZ());
      fHistTrackletPID->Fill(tracklet.GetPID());
      fHistTrackletYDy->Fill(tracklet.GetBinY(), tracklet.GetBinDy());
      fHistTrackletHC->Fill(tracklet.GetHCId()/60, tracklet.GetHCId()%60);

      if (TMath::Abs(tracklet.GetBinY()) >= 3682) 
        fHistTrackletBadY->Fill(tracklet.GetHCId()/60, (tracklet.GetHCId()%60)/12);

      if (tracklet.GetPID() < 40)
        fHistTrackletBadPID->Fill(tracklet.GetHCId()/60, (tracklet.GetHCId()%60)/12);

    }

    for (int iTrack = 0; iTrack < fGtuTrackArray->GetEntriesFast(); iTrack++) {
      AliESDTrdTrack *trdTrack = (AliESDTrdTrack*) ((*fGtuTrackArray)[iTrack]);
      /*
      printf("TRDMSG: GTU track found: 0x%016llx - Stack %02d_%d  pt=%.3fGeV/c\n",
	     trdTrack->GetTrackWord(0),
	     trdTrack->GetSector(), trdTrack->GetStack(),
	     trdTrack->Pt());
      */
      fHistTrackPt->Fill(trdTrack->Pt());
      fHistTrackStack->Fill(trdTrack->GetSector(), trdTrack->GetStack());
      fHistTrackPID->Fill(trdTrack->GetPID());

      Int_t layers=0;
      Int_t layer_mask = trdTrack->GetLayerMask();
      for (Int_t iLayer=0; iLayer<6; iLayer++)
       if ((layer_mask >> iLayer) & 1)
         layers++;
      fHistTrackletsTrack->Fill(layers);
      if (TMath::Abs(trdTrack->Pt()) >= 3.)
        fHistTrackletsTrackHpt->Fill(layers);
    }

    // do more complex processing here, preferably in a dedicated class

    // push  TObject-based data to output
    iResult = PushBack(fHistArray, 
		       (kAliHLTDataTypeTObject | kAliHLTDataOriginTRD), 
    		       pBlock->fSpecification);

    if (iResult < 0)                                  
      break;
		       
    // clear the rawreader
    fRawReaderMem->ClearBuffers();    
  }

  return iResult;
}

int AliHLTTRDMonitorComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/)
{
  // reconfigure the component from the specified CDB entry, or default CDB entry
  // function is invoked by the framework if a reconfigure command was received.
  // 
  int iResult=0;
  /*
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigSample/";
    cdbPath+=GetComponentID();
  }
  AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);
  */
  return iResult;
}

int AliHLTTRDMonitorComponent::ReadPreprocessorValues(const char* /*modules*/)
{
  // read the preprocessor values for the detectors in the modules list
  // function is invoked by the framework if the pendolino indivates an update
  // of online calibration objects, e.g temperature and pressure measurements.
  int iResult=0;
  /*
  TString detectors(modules!=NULL?modules:"");
  AliInfoClass(Form("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data()));
  */
  return iResult;
}
