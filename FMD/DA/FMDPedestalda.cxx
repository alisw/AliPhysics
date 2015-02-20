/*
  FMD DA for online calibration of conditions

  Contact:                 canute@nbi.dk
  Link:                    fmd.nbi.dk/fmd/offline
  Run Type:                PEDESTAL
  DA Type:                 LDC
  Number of events needed: 1000
  Input Files:             raw data 
  Output Files:            peds.csv
  Trigger types used:      PEDESTAL
*/
#include <AliFMDPedestalDA.h>
#include <AliFMDParameters.h>
#include <TROOT.h>
#include <TApplication.h>
#include "FMDUtilda.h"
#include <iostream>
#include <unistd.h>

int main(int argc, char **argv) 
{
  AliFMDBaseDA::Runner r;

  Int_t ret = r.Init(argc, argv);
  if (ret < 0) return -ret;
  if (ret > 0) return 0;

  AliFMDPedestalDA pedDA;
#ifdef ALI_AMORE
  pedDA.SetMakeSummaries(kTRUE);
#endif
  std::cout << "Executing pedestal DA" << std::endl;
  Bool_t success = r.Exec(pedDA);
  if (!success) {
    std::cout << "Failed to execute Pedestal DA" 
	      << std::endl;
    ret = 1;
  }
  
  if (success && r.fUpload) {
    const char* files[] = { "conditions.csv", 
			    "peds.csv", 
			    0 };
    const char* ids[] = { AliFMDParameters::Instance()->GetConditionsShuttleID(),
			  AliFMDParameters::Instance()->GetPedestalShuttleID(), 
			  0 };
    ret = UploadFiles(files, ids, r.fOwnUpload);
    if (ret > 0) std::cerr << "Upload of pedestals failed" << std::endl;
  }

  if (success) 
    PostSummaries(pedDA, "ped", r.RunNumber());

  gROOT->SetMustClean(false);

  if (r.fFast) {
    std::cout << "Sleep 1 second before getting out" << std::endl;
    gSystem->Sleep(1000);
    std::cout << "Now calling _exit(" << ret << ") to finish NOW!" << std::endl;
    _exit(ret);
  }

  std::cout << "End of FMD-Pedestal, return " << ret << std::endl;
  return ret;
}
//
// EOF
//

