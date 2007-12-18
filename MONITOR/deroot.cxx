// Author: Filimon Roukoutakis 02/08/2006

/******************************************************************************
  MOOD - Monitor Of On-line Data and Detector Debugger for ALICE Experiment
******************************************************************************/

#include <TError.h>
#include <TSysEvtHandler.h>
#ifdef ALI_DATE
#include "deroot.h"

int deroot(const char *rootFileName, const char *dateFileName);

int deroot(const char *rootFileName, const char *dateFileName) {

 unsigned char *dateEvent=new unsigned char [1000000000];

 FILE *dateFile;
 size_t gdcCounter, gdcSize;
 TFile rootFile(rootFileName);
 TTree *t=(TTree *)rootFile.Get("RAW");
 if(!t) {
  cerr << "Error getting RAW tree" << endl;
  return(1);
 }
 AliRawEvent *rootEvent=NULL;
 
 t->SetBranchAddress("rawevent", &rootEvent);

 if(!(dateFile=fopen(dateFileName, "wb"))) {
  cerr << "Error opening DATE file" << endl;
  return(1);
 }
 
 for(gdcCounter=0; gdcCounter<t->GetEntries(); gdcCounter++) {
  rootEvent=new AliRawEvent;
  t->GetEntry(gdcCounter);
  gdcSize=Root2Date(rootEvent, dateEvent);
  delete rootEvent;
  cerr << "\r     \r" << setprecision(3) << 100*(float)(gdcCounter+1)/t->GetEntries() << "% ";
  fwrite(dateEvent, gdcSize, 1, dateFile);
 }

 // Cleanup resources
 
 cerr << "\r     \r";
 cerr.flush();
 delete t;
 rootFile.Close();
 fclose(dateFile);
 delete [] dateEvent;
 
 return(0);

}

int main(int argc, char **argv) {

  if (argc != 3) {
    cerr << "Usage: deroot <input_root_file> <output_date_file>" << endl;
    return 1;
  }

 deroot(argv[1], argv[2]);

 return(0);

}

#else
int main(int /*argc*/, char** /*argv*/)
{
  ::Error("main", "this program was compiled without DATE");

  return 1;
}
#endif
