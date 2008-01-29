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

#include "AliRawEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawEquipment.h"
#include "AliRawEquipmentHeader.h"
#include "AliRawDataHeader.h"
#include "AliRawData.h"

#include <Riostream.h>

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
static bool DumpEvent(const char *progname, AliRawEvent *rawEvent)
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

  for(Int_t iSubEvent=0; iSubEvent < rawEvent->GetNSubEvents(); iSubEvent++) {
    AliRawEvent *rawSubEvent = rawEvent->GetSubEvent(iSubEvent);
    AliRawEventHeaderBase *rawSubEventHeader = rawSubEvent->GetHeader();
    cout << "    *********** Sub-event header ***********" << endl;
    rawSubEventHeader->Print("  ");

    for(Int_t iEquipment=0; iEquipment < rawSubEvent->GetNEquipments(); iEquipment++) {
      AliRawEquipment *rawEquip = rawSubEvent->GetEquipment(iEquipment);
      AliRawEquipmentHeader *rawEquipHeader = rawEquip->GetEquipmentHeader();
      cout << "      *********** Equipment event header ***********" << endl;
      rawEquipHeader->Print("    ");
      cout << "        *********** Common Data Header ***********" << endl;
      AliRawData *rawData = rawEquip->GetRawData();
      AliRawDataHeader *cdh = (AliRawDataHeader*)rawData->GetBuffer();
      if (!DumpCDH(cdh)) return false;
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

  AliRawEvent *rawEvent=NULL;
 
  rawTree->SetBranchAddress("rawevent", &rawEvent);

  Int_t nEvents = rawTree->GetEntries();

  cout << "*******************************************" << endl;
  cout << "File: " << argv[1] << endl;
  cout << "GUID: " << rawFile->GetUUID().AsString() << endl;
  cout << "Total number of events: " << nEvents << endl;
  cout << "*******************************************" << endl;

  for(Int_t iEvent=0; iEvent < nEvents; iEvent++) {
    rawEvent=new AliRawEvent;
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
