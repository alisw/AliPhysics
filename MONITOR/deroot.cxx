// Author: Filimon Roukoutakis 02/08/2006

/******************************************************************************
  MOOD - Monitor Of On-line Data and Detector Debugger for ALICE Experiment
******************************************************************************/

#include <TError.h>
#include <TSystem.h>
#include <TSysEvtHandler.h>
#include <TGrid.h>
#include "deroot.h"

int deroot(const char *rootFileName, const char *dateFileName, const char *ddlFilesFolder);

int deroot(const char *rootFileName, const char *dateFileName, const char *ddlFilesFolder) {

 TString str = rootFileName;
 if (str.BeginsWith("alien://"))
   TGrid::Connect("alien://");

 TFile *rootFile = TFile::Open(rootFileName,"READ");
 if (!rootFile) {
   cerr << "Raw data file can not be opened" << endl;
   return(1);
 }

 TTree *t=(TTree *)rootFile->Get("RAW");
 if(!t) {
  cerr << "Error getting RAW tree" << endl;
  return(1);
 }
 AliRawVEvent *rootEvent=NULL;
 
 t->SetBranchAddress("rawevent", &rootEvent);

 FILE *dateFile;
#ifdef __APPLE__
 if(!(dateFile=fopen(dateFileName, "wb"))) {
#else
 if(!(dateFile=fopen64(dateFileName, "wb"))) {
#endif
  cerr << "Error opening DATE file" << endl;
  return(1);
 }
 
 UInt_t eventSize = 10000000; // 10MB by default
 unsigned char *dateEvent = new unsigned char[eventSize];
 for(Long_t gdcCounter=0; gdcCounter<t->GetEntries(); gdcCounter++) {
  rootEvent=NULL;
  t->GetEntry(gdcCounter);
  if (rootEvent->GetHeader()->GetEventSize() > eventSize) {
    delete [] dateEvent;
    eventSize = (UInt_t)(1.05*rootEvent->GetHeader()->GetEventSize());
    dateEvent = new unsigned char[eventSize];
  }

  size_t gdcSize;
  if (ddlFilesFolder) {
    char command[256];
    sprintf(command, "rm -rf %s/raw%ld", ddlFilesFolder, gdcCounter);
    gSystem->Exec(command);
    sprintf(command, "%s/raw%ld", ddlFilesFolder, gdcCounter);
    if (gSystem->MakeDirectory(command) < 0) {
      cerr << "Can not create directory " << command << endl;
      return(1);
    }
    gdcSize=Root2Date(rootEvent, dateEvent, command);
  }
  else
    gdcSize=Root2Date(rootEvent, dateEvent, NULL);

  delete rootEvent;
  cerr << "\r     \r" << setprecision(3) << 100*(float)(gdcCounter+1)/t->GetEntries() << "% ";
  fwrite(dateEvent, gdcSize, 1, dateFile);
 }

 // Cleanup resources
 
 cerr << "\r     \r";
 cerr.flush();
 delete t;
 rootFile->Close();
 fclose(dateFile);
 delete [] dateEvent;
 
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

