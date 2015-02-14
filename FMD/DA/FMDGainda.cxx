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
  r.Exec(gainDA);

  const char* files[] = { "conditions.csv", 
			  "gains.csv", 
			  0 }; 
  const char* ids[] = { AliFMDParameters::Instance()->GetConditionsShuttleID(), 
			AliFMDParameters::Instance()->GetGainShuttleID(),
			0 };
  ret = UploadFiles(files, ids);

  if(ret > 0) std::cerr << "Gain DA failed" << std::endl;

  PostSummaries(gainDA, "gain", r.RunNumber());

  std::cout << "End of FMD-Gain, return " << ret << std::endl;

  gROOT->SetMustClean(false);

  std::cout << "Now calling _Exit(" << ret << ") to finish NOW!" << std::endl;
  _exit(ret);

  return ret;
}
//
// EOF
//

