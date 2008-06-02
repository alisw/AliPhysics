/*
  EMCAL DA for online calibration
  
  Contact: silvermy@ornl.gov
  Run Type: PEDESTAL or PHYSICS or STANDALONE
  DA Type: LDC
  Number of events needed: ~100
  Input Files: argument list
  Output Files: RESULT_FILE=EMCALda1.root, to be exported to the DAQ FXS
  fileId:  FILE_ID=EMCALda1    
  Trigger types used: CALIBRATION_EVENT (temporarily also PHYSICS_EVENT to start with)

*/
/*
  This process reads RAW data from the files provided as command line arguments
  and save results (class itself) in a file (named from RESULT_FILE define - see below).
*/

#define RESULT_FILE  "EMCALda.root"
#define FILE_ID "EMCALda"
#define AliDebugLevel() -1
#define ClassName "emcCalibPedestal"

extern "C" {
#include <daqDA.h>
}
#include "event.h"
#include "monitor.h"

#include "stdio.h"
#include "stdlib.h"

//
//AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliRawEventHeaderBase.h"
#include "AliCaloRawStream.h"
#include "AliLog.h"

//
// EMC calibration-helper algorithm includes
//
#include "AliCaloCalibPedestal.h"

/*
  Main routine, EMC pedestal detector algorithm to be run on EMC LDC
  Arguments: list of DATE raw data files
*/

int main(int argc, char **argv) {

  AliLog::SetClassDebugLevel("AliCaloRawStream",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);

  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  int i, status;

  /* log start of process */
  printf("EMCAL DA started - %s\n",__FILE__);

  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }

  AliCaloCalibPedestal * calibPedestal = new 
    AliCaloCalibPedestal(AliCaloCalibPedestal::kEmCal); // pedestal and noise calibration

  /* loop over RAW data files */
  int nevents=0;
  for ( i=1; i<argc; i++ ) {

    /* define data source : this is argument i */
    printf("Processing file %s\n", argv[i]);

    AliRawReader *rawReader = new AliRawReaderDate(argv[i]);
    AliCaloRawStream *in = new AliCaloRawStream(rawReader,"EMCAL");
    AliRawEventHeaderBase *aliHeader=NULL;

    /* read until EOF */
    while ( rawReader->NextEvent() ) {

      /* check shutdown condition */
      if (daqDA_checkShutdown()) {break;}

      aliHeader = (AliRawEventHeaderBase*) rawReader->GetEventHeader();
      calibPedestal->SetRunNumber( aliHeader->Get("RunNb") ); // just for fun; keep info on last run looked at

      // select physics and calibration events now (only calibration in future)
      if ( aliHeader->Get("Type") == AliRawEventHeaderBase::kPhysicsEvent || 
	   aliHeader->Get("Type") == AliRawEventHeaderBase::kCalibrationEvent  ) {

	nevents++;

	//  Pedestal calibration
	calibPedestal->ProcessEvent(in);
      }
    } // loop over all events in file
    /* cleanup the reading handles */
    delete in;
    delete rawReader;    
  } // loop over files

  //
  // write class to rootfile
  //

  printf ("%d physics/calibration events processed.\n",nevents);

  TFile f(RESULT_FILE, "recreate");
  if (!f.IsZombie()) { 
    f.cd();
    calibPedestal->Write(ClassName);
    f.Close();
    printf("Object saved to file \"%s\" as \"%s\".\n", RESULT_FILE, ClassName); 
  } 
  else {
    printf("Could not save the object to file \"%s\".\n", RESULT_FILE);
  }

  /* store the result file on FES */
  status = daqDA_FES_storeFile(RESULT_FILE, FILE_ID);

  // see if we can delete our analysis helper also
  delete calibPedestal;

  return status;
}
