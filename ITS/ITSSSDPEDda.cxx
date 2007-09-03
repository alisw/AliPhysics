
#include <iostream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "daqDA.h"
//#include "AliITSChannelDaSSD.h" 
#include "AliITSHandleDaSSD.h" 

using namespace std;

Bool_t GetRunSettings (const char *datafilename, Long_t &eventsnumber, Long_t &stripsnumber);

int main( int argc, char** argv )
{
  AliITSHandleDaSSD  *ssddaldc;
  ostringstream       feefname;
  Int_t               status;
  string              dafname;
  Char_t             *dadaqdir = NULL;

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  char *datafilename = argv[1];

  Long_t strn, evn;
  if (!GetRunSettings (datafilename, evn, strn))
  {
    cout << "Error GetRunSettings (datafilename, evn, strn)!" << endl;
    return -1;
  }
  cout << "Pysics events : " << evn << ";   Total strip number : " << strn 
       << ";   Modules number: " << strn / AliITSModuleDaSSD::GetStripsPerModuleConst() << endl;
  
  ssddaldc = new AliITSHandleDaSSD();
//  if (!ssddaldc->SetNumberOfModules((Int_t)(strn / AliITSModuleDaSSD::GetStripsPerModuleConst())))
  if (!ssddaldc->SetNumberOfModules(AliITSHandleDaSSD::GetNumberOfSSDModulesConst()))
  {
     cout << "Error ssddaldc->SetNumberOfModules" << endl;
     delete ssddaldc;
     return -1;
  }  
  if (!ssddaldc->ReadCalibrationDataFile(datafilename, evn))
  {
     cout << "Error !ssddaldc->ReadCalibrationDataFile" << endl;
     delete ssddaldc;
     return -1;
  }  
  daqDA_progressReport(30);
//  if (daqDA_checkShutdown() == 1) {
//    cout << "Shutdown has been requested!" << endl;
//    delete ssddaldc;
//    return -1;
//  }

  if (!ssddaldc->CalculatePedestal()) {
    cout << "Error, ssddaldc->CalculatePedestal()";
    return 1; 
  } 
  daqDA_progressReport(50);
  if (!ssddaldc->CalculateNoiseCM()) {
    cout << "Error, ssddaldc->CalculateNoiseCM()";
    return 2; 
  } 
  ssddaldc->DeleteSignal();
  daqDA_progressReport(90);
  dadaqdir = getenv ("DAQDA_TEST_DIR");
  if (dadaqdir) {
    dafname = dadaqdir;
    if (!(ssddaldc->SaveCalibrationSSDLDC(dafname))) 
      cout << "Error saving DA data to the file! Probably $DAQDA_TEST_DIR defined incorrectly!" << endl;
    else cout << "SSDDA data are saved in " << dafname << endl;
      feefname << dadaqdir << "/ssddaldc.root";
    cout << "Saving feessdda data in " << feefname.str() << endl;
    TFile *fileRun = new TFile (feefname.str().data(),"RECREATE");
    ssddaldc->Write();
    fileRun->Close();
    delete fileRun;
    status = daqDA_FES_storeFile(dafname.data(),"DASSD_DB_results");
    if (status) printf("Failed to export file : %d\n",status);
  }
  else cout << "Error: DAQDA_TEST_DIR is not defined, DA data are not saved!" << endl;
  delete ssddaldc;
  daqDA_progressReport(100);
  return 0;
}



Bool_t GetRunSettings (const char *datafilename, Long_t &eventsnumber, Long_t &stripsnumber)
{
  Long_t physeventind = 0, strn = 0, strneq = 0;
  AliRawReaderDate  *rawreaderdate = NULL;
  Int_t *data = NULL;
  Long_t datasize = 0, eqdatasize = 0, eqbelsize = 1;
  rawreaderdate = new AliRawReaderDate(datafilename, 0);
  if (!rawreaderdate) {
    cout << "GetRunSettings : Error  new DARawReader(datafilename, 0);" << endl;
    return kFALSE;
  }  
  rawreaderdate->SelectEvents(PHYSICS_EVENT);
  while (rawreaderdate->NextEvent())
  { 
    physeventind += 1;
    datasize = 0;
    strn = 0;
    while (rawreaderdate->ReadNextData((UChar_t*&)data)) {
      eqdatasize = rawreaderdate->GetDataSize();
      eqbelsize = rawreaderdate->GetEquipmentElementSize();
      if ( (eqdatasize % eqbelsize) || (eqbelsize != sizeof(long32)) ) {
        cout << "Error ReadCalibrationDataFile: equipment event data size " << eqdatasize  
	     << " is not an integer of equipment data size "      << eqbelsize << endl;
	rawreaderdate->DumpData();
        return kFALSE;
      }
      strneq = eqdatasize / eqbelsize;
      datasize += eqdatasize;
      strn += strneq;
    }  
    if ((strn * eqbelsize != datasize)) {
      if (physeventind != 1) {
        cout << "Something is wrong with data file, strip number changes from event to event! Ev = " << physeventind << endl;
        rawreaderdate->DumpData();
        return kFALSE;
      }
      if ((datasize % eqbelsize)) {
        cout << "Wrong number :  (datasize % eqbelsize) != 0" << endl;
        rawreaderdate->DumpData();
        return kFALSE;
      }
      strneq = datasize / eqbelsize;
    }
  }
  delete rawreaderdate;
  if ((physeventind > 0) && (strn > 0))
  {
    eventsnumber = physeventind;
    stripsnumber = strn;
    return kTRUE;
  }
  return kFALSE;
}
