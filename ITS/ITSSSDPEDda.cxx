/**************************************************************************
- "Contact": - Oleksandr Borysov, aborysov@ts.infnf.it
- "Link:" - link to test files: /afs/infn.it/ts/user/aborysov/public/C23_run387.000.raw
- "Run Type:" - run type (exactly as defined in the ECS)
- "DA Type:" -  LDC
- "Number of events needed:"   at least 500

- "Input Files:" - config file:    ssdpeddaconfig 
                   previous result files:  
		   data source:   raw data file on LDC
		    
- "Output Files:" -  local names $DA_TEST_DIR/ssddaldc_<LDCID>_<RunID>.root 
                     FXS name: ITSSSDda_<LDCID>_<RunID>.root, 
                     local files are persistent over runs: data source
- "Trigger types used:"
 **************************************************************************/


#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "daqDA.h"
#include "AliITSHandleDaSSD.h" 

using namespace std;


int main( int argc, char** argv )
{
  AliITSHandleDaSSD  *ssddaldc;
  TString             feefname, cmddbsave;
  Int_t               status;
  Char_t             *dafname = NULL, *dadaqdir = NULL;

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  char *datafilename = argv[1];

  ssddaldc = new AliITSHandleDaSSD(datafilename);
  if (ssddaldc->IsZombie()) return -1;
  if (!ssddaldc->ProcessRawData())
  {
     cout << "Error !ssddaldc->ProcessRawData()" << endl;
     delete ssddaldc;
     return -1;
  }  
  daqDA_progressReport(90);

  if (!system(NULL)) {
    cout << "Error: the call system(NULL) in main() returned NULL!" << endl;
    return -1;
  }
  dadaqdir = getenv ("DA_TEST_DIR");
  if (dadaqdir) {
    dafname = dadaqdir;
    if (!(ssddaldc->SaveCalibrationSSDLDC(dafname))) 
      cout << "Error saving DA data to the file! Probably $DA_TEST_DIR defined incorrectly!" << endl;
    else cout << "SSDDA data are saved in " << dafname << endl;
    feefname = Form("%s/ssddaldc_%i_%i.root", dadaqdir, ssddaldc->GetLdcId(), ssddaldc->GetRunId());
    cout << "Saving feessdda data in " << feefname << endl;
    TFile *fileRun = new TFile (feefname.Data(),"RECREATE");
    ssddaldc->Write();
    fileRun->Close();
    delete fileRun;
    status = daqDA_FES_storeFile(dafname, "DASSD_DB_results");
    if (status) printf("Failed to export file : %d\n",status);

    if (getenv("DATE_DB_DIR")) {
      cmddbsave = Form("$DATE_DB_DIR/daqDetDB_store ssddaldc.root %s", feefname.Data());
      status = system(cmddbsave.Data());
      if (status) printf("Failed to export file to the detector db: %d\n",status);
    } else cout << "Error main(): $DATE_DB_DIR is not defined!" << endl;
  }
  else cout << "Error: DA_TEST_DIR is not defined, DA data are not saved!" << endl;
  delete ssddaldc;
  daqDA_progressReport(100);
  return 0;
}