/*

  FMD DA for online calibration of conditions

  Contact:                 canute@nbi.dk
  Link:                    fmd.nbi.dk/fmd/offline
  Run Type:                GAIN
  DA Type:                 LDC
  Number of events needed: usually 102400
  Input Files:             raw data 
  Output Files:            gains.csv
  Trigger types used:      GAIN
*/
#include <AliFMDGainDA.h>
#include <AliFMDParameters.h>
#include "FMDUtilda.h"
#include <TROOT.h>
#include <iostream>
#include <unistd.h>

int main(int argc, char **argv) 
{
  AliFMDBaseDA::Runner r;
  
  Int_t ret = r.Init(argc, argv);
  if (ret < 0) return -ret;
  if (ret > 0) return 0;

  AliFMDGainDA gainDA;
#ifdef ALI_AMORE
  gainDA.SetMakeSummaries(kTRUE);
#endif
  std::cout << "Running Gain DA" << std::endl;
  Bool_t success = r.Exec(gainDA);
  if (!success) { 
    std::cout << "Failed to execute Gain DA" 
	      << std::endl;
    ret = 1;
  }
  
  if (success && r.fUpload) {
    const char* files[] = { "conditions.csv", 
			    "gains.csv", 
			    0 }; 
    const char* ids[] = { AliFMDParameters::Instance()->GetConditionsShuttleID(), 
			  AliFMDParameters::Instance()->GetGainShuttleID(),
			  0 };
    ret = UploadFiles(files, ids, r.fOwnUpload);
    if (ret > 0) std::cerr << "Upload of gains failed" << std::endl;
  }

  if (success) 
    PostSummaries(gainDA, "gain", r.RunNumber());


  gROOT->SetMustClean(false);
  if (r.fFast) {
    std::cout << "Sleep 1 second before getting out" << std::endl;
    gSystem->Sleep(1000);
    std::cout << "Now calling _exit(" << ret << ") to finish NOW!" << std::endl;
    _exit(ret);
  }

  std::cout << "End of FMD-Gain, return " << ret << std::endl;
  return ret;
}
//
// EOF
//

