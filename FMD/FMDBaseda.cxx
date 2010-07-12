/*

  FMD DA for online calibration of conditions

  Contact:                 canute@nbi.dk
  Link:                    fmd.nbi.dk/fmd/offline
  Run Type:                PHYSICS
  DA Type:                 MON
  Number of events needed: depending on the run, being run-level
  Input Files:             raw data 
  Output Files:            conditions.csv
  Trigger types used:      PHYSICS_EVENT
*/
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
  
  
  const Char_t* tableSOD[]  = {"ALL", "no", "SOD", "all", NULL, NULL};

  Bool_t old = kTRUE;
  monitorDeclareTable(const_cast<char**>(tableSOD));

  AliFMDParameters::Instance()->Init(kFALSE,0);
  AliFMDParameters::Instance()->UseCompleteHeader(old);
  struct eventHeaderStruct *event;
  int status;
  /* define data source : this is argument 1 */  
  if (argc < 2) { 
    std::cout << "No monitor source set" << std::endl;
    return -1;
  }
  status=monitorSetDataSource( argv[1] );
  if (status!=0) {
    printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
    return -1;
  }
  
  if (argc > 2) AliLog::SetModuleDebugLevel("FMD", atoi(argv[2]));
  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }
  
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);

  AliRawReader *reader = 0;
  AliFMDBaseDA baseDA;
  Int_t  retval = 0;
  Int_t iev = 0;
  Bool_t SODread = kFALSE;
  while(!SODread && iev<1000) {
    
    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}
    
    /* get next event (blocking call until timeout) */
    status=monitorGetEventDynamic((void **)&event);
    if (status==MON_ERR_EOF) {
      printf ("End of File detected\n");
      break; /* end of monitoring file has been reached */
    }
    
    if (status!=0) {
      printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
      break;
    }
    
    /* retry if got no event */
    if (event==NULL) {
      continue;
    }
    
    iev++; 
    
    switch (event->eventType){
      
    case START_OF_DATA:
      reader = new AliRawReaderDate((void*)event);
      baseDA.Run(reader);
      SODread = kTRUE;
      retval = 
	daqDA_FES_storeFile("conditions.csv", 
			    AliFMDParameters::Instance()->GetConditionsShuttleID());
      if (retval != 0) std::cerr << "Base DA failed" << std::endl;
      
      break;
    
    default:
      break;
    
    }
  }
 
   
  return retval;
}
