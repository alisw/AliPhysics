/*
  FMD DA for online calibration of conditions

  Contact:                 christian.holm.christensen@cern.ch
  Link:                    fmd.nbi.dk/fmd/offline
  Run Type:                PHYSICS
  DA Type:                 MON
  Number of events needed: depending on the run, being run-level
  Input Files:             raw data 
  Output Files:            conditions.csv
  Trigger types used:      PHYSICS_EVENT
*/
#include <cstdlib>
#include <Riostream.h>
#include "monitor.h"
#include "event.h"
#include <AliLog.h>
#include <TSystem.h>
#include <TString.h>
#include <AliFMDParameters.h>
#include <AliRawReader.h>
#include <TStopwatch.h>
#include <AliFMDBaseDA.h>
#include <AliRawReaderDate.h>
#include <AliRawReaderRoot.h>
#include "daqDA.h"
#include "TROOT.h"
#include "TPluginManager.h"

void
usage(std::ostream& o, const char* progname)
{
  o << "Usage: " << progname << " FILE [OPTIONS]\n\n"
    << "Options:\n"
    << "\t-h,--help         Show this help\n"
    << "\t-d,--diagnostics  Create diagnostics\n"
    << "\t-D,--debug LEVEL  Set the debug level\n"
    << std::endl;
}

int main(int argc, char **argv) 
{

#if 0
  /* magic line from Rene - for future reference! */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");
#endif
  AliFMDBaseDA::Runner r;

  Int_t ret = r.Init(argc, argv, false);
  if (ret < 0) return -ret;
  if (ret > 0) return 0;

  
  
  const Char_t* tableSOD[]  = {"ALL", "no", "SOD", "all", NULL, NULL};

  monitorDeclareTable(const_cast<char**>(tableSOD));

  int status = monitorSetDataSource(r.fSource.Data());
  if (status!=0) {
    printf("monitorSetDataSource() failed for %s: %s\n",
	   r.fSource.Data(), monitorDecodeError(status));
    return -1;
  }
  
  /* declare monitoring program */
  status = monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }
  
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);

  AliFMDBaseDA  baseDA;
  Int_t         iev     = 0;
  Bool_t        sodSeen = false;
  Bool_t        success = false;
  while(!sodSeen && iev<1000) {
    
    /* check shutdown condition */
    if (daqDA_checkShutdown()) break;
    
    /* get next event (blocking call until timeout) */
    struct eventHeaderStruct *event = 0;    
    status = monitorGetEventDynamic((void **)&event);
    if (status == MON_ERR_EOF) {
      printf ("End of File detected\n");
      break; /* end of monitoring file has been reached */
    }
    
    if (status != 0) {
      printf("monitorGetEventDynamic() failed : %s\n",
	     monitorDecodeError(status));
      break;
    }
    
    /* retry if got no event */
    if (!event) continue;
    
    iev++; 
    
    switch (event->eventType) {
    case START_OF_DATA: {
      std::cout << "Got START OF DATA event" << std::endl;
      AliRawReader* reader = new AliRawReaderDate((void*)event);
      if (!(success = baseDA.Run(reader, r.fAppendRun, true)))     
	std::cout << "Base DA failed" << std::endl;
      sodSeen = kTRUE;
    }
      break;
    default:
      break;
    
    }
  }
  int retval = success ? 0 : 1;
  if (r.fUpload && success) {
    std::cout << "Pushing to FXS" << std::endl;
    retval = 
      daqDA_FES_storeFile("conditions.csv", 
			  AliFMDParameters::Instance()
			  ->GetConditionsShuttleID());
    if (retval != 0) std::cerr << "Base DA failed" << std::endl;
      
  }
  std::cout << "End of FMD-Base - return " << retval << std::endl;
   
  return retval;
}
//
// EOF
// 
