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
#include "FMDUtilda.h"
#include <iostream>

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
  
  return ret;
}
//
// EOF
//

