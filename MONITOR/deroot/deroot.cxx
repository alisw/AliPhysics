// Author: Filimon Roukoutakis 02/08/2006

/******************************************************************************
  MOOD - Monitor Of On-line Data and Detector Debugger for ALICE Experiment
******************************************************************************/

#include <RConfig.h>
#include <TError.h>
#include <TSystem.h>
#include <TSysEvtHandler.h>
#include <TGrid.h>
#include "deroot.h"
#include "AliRawReader.h"

using std::cerr;
using std::endl;

int deroot(const char *rootFileName, const char *dateFileName, const char *ddlFilesFolder);

int deroot(const char *rootFileName, const char *dateFileName, const char *ddlFilesFolder) {

  AliRawReader* rawReader = AliRawReader::Create(rootFileName);
  
  if (!rawReader)
  {
    cerr << "Error creating raw reader" << endl;
    return(1);
  }
  
  std::cout << "Getting number of input events..." << std::flush;
  
  Int_t nEntries = rawReader->GetNumberOfEvents();
  
  std::cout << nEntries << std::endl;
  
  Bool_t removeEOREvents = kFALSE;

  if ( TString(rawReader->ClassName()).Contains("Chain") )
  {
	  // reading from a (local) chain.
	  // that's a special mode where EOR does not make a lot of sense, so remove those
	  removeEOREvents = kTRUE;
  }

 FILE *dateFile;
#if defined(R__SEEK64)
 if(!(dateFile=fopen64(dateFileName, "wb"))) {
#else
 if(!(dateFile=fopen(dateFileName, "wb"))) {
#endif
  cerr << "Error opening DATE file" << endl;
  return(1);
 }
 
 UInt_t eventSize = 10000000; // 10MB by default
 unsigned char *dateEvent = new unsigned char[eventSize];
 Long_t gdcCounter(0);
 ULong_t totalSize(0);
   
 while ( rawReader->NextEvent() )
 {
   
   rawReader->Reset();
   
   AliRawVEvent* rootEvent = const_cast<AliRawVEvent*>(rawReader->GetEvent());
   
    if ( removeEOREvents && TString(rootEvent->GetHeader()->GetTypeName()).Contains("END_OF_RUN") )
   {
	   continue;
   }

  if (rootEvent->GetHeader()->GetEventSize() > eventSize) {
    delete [] dateEvent;
    eventSize = (UInt_t)(1.05*rootEvent->GetHeader()->GetEventSize());
    dateEvent = new unsigned char[eventSize];
  }

  size_t gdcSize;
  if (ddlFilesFolder) {
    TString eventDir;
    
    eventDir.Form("%s/raw%ld",ddlFilesFolder,gdcCounter);
    
    TString cmd;
    
    cmd.Form("rm -rf %s",eventDir.Data());
    
    gSystem->Exec(cmd.Data());

    if (gSystem->MakeDirectory(eventDir.Data()) < 0) {
      cerr << "Can not create directory " << eventDir.Data() << endl;
      fclose(dateFile);
      delete [] dateEvent;
      return(1);
    }
    gdcSize=Root2Date(rootEvent, dateEvent, eventDir.Data());
  }
  else {
	  gdcSize=Root2Date(rootEvent, dateEvent, NULL);
  }

   ++gdcCounter;
   
   totalSize += gdcSize;
   cout << Form("\r    \r %10ld events written (%9.2f MB)",gdcCounter,totalSize/1024.0/1024.0);
   
   if (nEntries>0)
   {
     cout << Form(" [ %4.1f %% ]",100.0*gdcCounter/nEntries);
   }
   
    size_t written = fwrite(dateEvent, gdcSize, 1, dateFile);

    if (!written) break;
 }

 // Cleanup resources
 
 cout << "\r     \r";
 cout.flush();
 fclose(dateFile);
 delete [] dateEvent;
 delete rawReader;
 
 return(0);

}

int main(int argc, char **argv) {

  if (argc != 3 && argc != 4) {
    cerr << "Usage: deroot <input_root_file> <output_date_file> [<optional_folder_for_ddl_files>]" << endl;
    return 1;
  }

  if (argc ==3)
    deroot(argv[1], argv[2], NULL);
  else
    deroot(argv[1], argv[2], argv[3]);

 return(0);

}

