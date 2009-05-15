// Author: Cvetan Cheshkov 29/01/2008

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// alirawdump                                                           //
//                                                                      //
// Program can used to dump the raw-data files in ROOT format.          //
// It dumps the event,sub-event,equipment and common-data header.       //
// Additional application of the program is to check if the CDHs        //
// of different raw-data payloads are compatible. In this sense         //
// it replaces the DAQ online checks in case the DAQ is running         //
// UNCHECKED partition.                                                 //
//                                                                      //
// Written by: Cvetan Cheshkov, 29/01/2008.                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>

#include "AliRawVEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawVEquipment.h"
#include "AliRawEquipmentHeader.h"
#include "AliRawDataHeader.h"
#include "AliRawData.h"
#include "AliDAQ.h"

#include <Riostream.h>

static Int_t miniEventIDOffset[AliDAQ::kNDetectors] = {3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565,3565};
static Bool_t detTriggerClasses[AliDAQ::kNDetectors] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

//______________________________________________________________________________
static void Usage(const char *prognam)
{
  // Prints the usage
  // of the alirawdump program
      fprintf(stderr, "Usage: %s <raw_data_root_file>\n",
              prognam);
      fprintf(stderr, " <raw_data_root_file> = file with ROOT formatted raw data\n");
}

//______________________________________________________________________________
static bool DumpCDH(AliRawDataHeader *cdh)
{
  // Dumps the CDH
  // ...
  cout << "        Size: " << cdh->fSize << endl;
  cout << "        Version: " << (Int_t)cdh->GetVersion() << endl;
  cout << "        Orbit: " << cdh->GetEventID2() << " Bunch-crossing: " << cdh->GetEventID1() << endl;
  cout << "        L1 trigger message: " << (UInt_t)cdh->GetL1TriggerMessage() << endl;
  cout << "        Participating sub-detectors: " << cdh->GetSubDetectors() << endl;
  cout << "        Block attributes: " << (Int_t)cdh->GetAttributes() << endl;
  cout << "        Status: " << cdh->GetStatus() << endl;
  cout << "        Mini event ID: " << cdh->GetMiniEventID() << endl;
  cout << "        Trigger classes: " << cdh->GetTriggerClasses() << endl;
  cout << "        ROI: " << cdh->GetROI() << endl;

  return true;
}

//______________________________________________________________________________
static bool CheckCDH(AliRawDataHeader *cdhRef,AliRawDataHeader *cdh)
{
  // Check the consistency of the CDHs
  // ...
  bool iserror = false;
  if ((cdhRef->GetEventID1() != cdh->GetEventID1())) {
    cout << "ERROR: CDH mismatch detected in EventID1: " <<  cdhRef->GetEventID1() << " != " << cdh->GetEventID1() << endl;
    iserror = true;
  }
//   if ((cdhRef->GetVersion() != cdh->GetVersion())) {
//     cout << "ERROR: CDH mismatch detected in Version: " <<  (Int_t)cdhRef->GetVersion() << " != " << (Int_t)cdh->GetVersion() << endl;
//     iserror = true;
//   }
  if ((cdhRef->GetEventID2() != cdh->GetEventID2())) {
    cout << "ERROR: CDH mismatch detected in EventID2: " <<  cdhRef->GetEventID2() << " != " << cdh->GetEventID2() << endl;
    iserror = true;
  }
//   if ((cdhRef->GetMiniEventID() != cdh->GetMiniEventID())) {
//     cout << "ERROR: CDH mismatch detected in MiniEventID: " <<  cdhRef->GetMiniEventID() << " != " << cdh->GetMiniEventID() << endl;
//     iserror = true;
//   }
//   if ((cdhRef->GetTriggerClasses() != cdh->GetTriggerClasses())) {
//     cout << "ERROR: CDH mismatch detected in TriggerClasses: " <<  cdhRef->GetTriggerClasses() << " != " << cdh->GetTriggerClasses() << endl;
//     iserror = true;
//   }

//   if ((cdhRef->GetL1TriggerMessage() != cdh->GetL1TriggerMessage())) {
//     cout << "ERROR: CDH mismatch detected in L1TriggerMessage: " <<  (Int_t)cdhRef->GetL1TriggerMessage() << " != " << (Int_t)cdh->GetL1TriggerMessage() << endl;
//     iserror = true;
//   }
  if ((cdhRef->GetSubDetectors() != cdh->GetSubDetectors())) {
    cout << "ERROR: CDH mismatch detected in ParticipatingSubDetectors: " <<  cdhRef->GetSubDetectors() << " != " << cdh->GetSubDetectors() << endl;
    iserror = true;
  }

  if (iserror) return false;
  else return true;
}

//______________________________________________________________________________
static bool DumpEvent(const char *progname, AliRawVEvent *rawEvent)
{
  // Dumps and checks one
  // raw-data event
  AliRawEventHeaderBase *rawEventHeader = rawEvent->GetHeader();

  if (rawEventHeader->GetMagic() != 0xDA1E5AFE) {
    Error(progname,"Wrong magic number ( 0x%x != 0xDA1E5AFE )",rawEventHeader->GetMagic());
    return false;
  }

  cout << "  *********** Event header ***********" << endl;
  rawEventHeader->Print();

  AliRawDataHeader *cdhRef = NULL;

  for(Int_t iSubEvent=0; iSubEvent < rawEvent->GetNSubEvents(); iSubEvent++) {
    AliRawVEvent *rawSubEvent = rawEvent->GetSubEvent(iSubEvent);
    AliRawEventHeaderBase *rawSubEventHeader = rawSubEvent->GetHeader();
    cout << "    *********** Sub-event header ***********" << endl;
    rawSubEventHeader->Print("  ");

    for(Int_t iEquipment=0; iEquipment < rawSubEvent->GetNEquipments(); iEquipment++) {
      AliRawVEquipment *rawEquip = rawSubEvent->GetEquipment(iEquipment);
      AliRawEquipmentHeader *rawEquipHeader = rawEquip->GetEquipmentHeader();
      cout << "      *********** Equipment event header ***********" << endl;
      rawEquipHeader->Print("    ");
      cout << "        *********** Common Data Header ***********" << endl;
      AliRawData *rawData = rawEquip->GetRawData();
      AliRawDataHeader *cdh = (AliRawDataHeader*)rawData->GetBuffer();

      Int_t ddlID;
      Int_t detID = AliDAQ::DetectorIDFromDdlID(rawEquipHeader->GetId(),ddlID);
      Int_t idOffset = cdh->GetMiniEventID() - cdh->GetEventID1();
      if (idOffset < 0) idOffset += 3564;
      if (miniEventIDOffset[detID] == 3565) {
	miniEventIDOffset[detID] = idOffset;
	cout << "MiniEvenID offset for detector " << AliDAQ::DetectorName(detID) << " is set to " << idOffset << endl;
      }
      else {
	if (miniEventIDOffset[detID] != idOffset) {
	  cout << "ERROR: MiniEventID offset for detector " << AliDAQ::DetectorName(detID) << " has changed ( " << idOffset << " != " << miniEventIDOffset[detID] << " )" << endl;
	}
      }

      // TPC is using version 1
      if ((cdh->GetVersion() != 2) && (detID != 3))
	cout << "ERROR: Bad CDH version: " << (Int_t)cdh->GetVersion() << endl;

      if (cdh->GetTriggerClasses() == 0) {
	if (detTriggerClasses[detID])
	  cout << "Empty trigger class mask for detector " << AliDAQ::DetectorName(detID) << endl;
	detTriggerClasses[detID] = false;
      }

      if (!DumpCDH(cdh)) return false;
      // check the CDH consistency
      if (cdhRef == NULL) {
	cdhRef = cdh;
      }
      else {
	// TPC L1 trigger message is shifted by 2 bits??
	UShort_t l1Message = cdh->GetL1TriggerMessage();
	UShort_t l1MessageRef = cdhRef->GetL1TriggerMessage();

	if (l1Message != l1MessageRef)
	  cout << "ERROR: CDH mismatch detected in L1TriggerMessage for detector " << AliDAQ::DetectorName(detID) << ": " << (Int_t)l1MessageRef << " ( " << (Int_t)cdhRef->GetL1TriggerMessage() << " ) " << " != " << (Int_t)l1Message << " ( " << (Int_t)cdh->GetL1TriggerMessage() << " )" << endl;

	if ((cdhRef->GetTriggerClasses() == 0) && (cdh->GetTriggerClasses() != 0)) {
	  // update the reference trigger class mask
	  cdhRef->fTriggerClassLow = cdh->fTriggerClassLow;     
	  cdhRef->fROILowTriggerClassHigh = (((cdhRef->fROILowTriggerClassHigh >> 28) & 0xF) << 28) | (cdh->fROILowTriggerClassHigh & 0x1FFFF);
	}
	if (cdh->GetTriggerClasses() != 0) {
	  if (cdhRef->GetTriggerClasses() != cdh->GetTriggerClasses()) {
	    cout << "ERROR: CDH mismatch detected in TriggerClasses: " <<  cdhRef->GetTriggerClasses() << " != " << cdh->GetTriggerClasses() << endl;
	  }
	}

	CheckCDH(cdhRef,cdh);
	//	if (!CheckCDH(cdhRef,cdh)) return false;
      }
    }
  }

  return true;
}

//______________________________________________________________________________
int main(int argc, char **argv)
{
  // Dumps a ROOT formatted
  // raw-data file

  gROOT->SetBatch();
  
  if ((argc == 2 && (!strcmp(argv[1], "-?") || !strcmp(argv[1], "-help"))) || argc != 2) {
      Usage(argv[0]);
      return 1;
  }

  TString str = argv[1];
  if (str.BeginsWith("alien://"))
    TGrid::Connect("alien://");

  TFile *rawFile = TFile::Open(argv[1],"READ");
  if (!rawFile) {
    Error(argv[0],"Raw data file %s can not be opened!",argv[1]);
    return 1;
  }

  TTree *rawTree=(TTree *)rawFile->Get("RAW");
  if(!rawTree) {
    Error(argv[0],"Error getting RAW tree from file %s",argv[1]);
    return 1;
  }

  AliRawVEvent *rawEvent=NULL;
 
  rawTree->SetBranchAddress("rawevent", &rawEvent);

  Int_t nEvents = rawTree->GetEntries();

  cout << "*******************************************" << endl;
  cout << "File: " << argv[1] << endl;
  cout << "GUID: " << rawFile->GetUUID().AsString() << endl;
  cout << "Total number of events: " << nEvents << endl;
  cout << "*******************************************" << endl;

  for(Int_t iEvent=0; iEvent < nEvents; iEvent++) {
    rawEvent=NULL;
    rawTree->GetEntry(iEvent);
    cout << "  *********** Event " << iEvent << " *******" << endl;
    DumpEvent(argv[0],rawEvent);
    delete rawEvent;
  }

  cout << "*******************************************" << endl;
  cout << "EOF" << endl;
  cout << "*******************************************" << endl;
  delete rawTree;
  rawFile->Close();
}
