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
  std::cout << "Eecuting pedestal DA" << std::endl;
  r.Exec(pedDA);

  const char* files[] = { "conditions.csv", 
			  "peds.csv", 
			  0 };
  const char* ids[] = { AliFMDParameters::Instance()->GetConditionsShuttleID(),
			AliFMDParameters::Instance()->GetPedestalShuttleID(), 
			0 };
  ret = UploadFiles(files, ids);

  if(ret > 0) std::cerr << "Pedestal DA failed" << std::endl;

  PostSummaries(pedDA, "ped", r.RunNumber());

  std::cout << "End of FMD-Gain, return " << ret << std::endl;
  gROOT->SetMustClean(false);

  std::cout << "Now calling _Exit(" << ret << ") to finish NOW!" << std::endl;
  _exit(ret);

  return ret;
}
//
// EOF
//

