/*
TPC DA for online calibration

Contact: Haavard.Helstrup@cern.ch
Link:
Run Type: CALIBRATION_PULSER
DA Type: LDC
Number of events needed: 100
Input Files: 
Output Files: tpcPulser.root, to be exported to the DAQ FXS
fileId:   pulser
Trigger types used: CALIBRATION_EVENT

*/

/*

TPCda_pulser.cxx - calibration algorithm for TPC pulser events

10/06/2007  sylvain.chapeland@cern.ch  :  first version - clean skeleton based on DAQ DA case1
30/09/2007  haavard.helstrup@cern.ch   :  created pulser DA based on pedestal code
19/09/2008  J.Wiechula@gsi.de:            Added export of the calibration data to the AMORE data base.
                                          Added support for configuration files.
23/04/2011  Christian.Lippmann@cern.ch :  Added output of acsii files for online

 contact: marian.ivanov@cern.ch


This process reads RAW data from the files provided as command line arguments
and save results in a file (named from RESULT_FILE define - see below).

*/

#define RESULT_FILE "tpcPulser.root"
#define FILE_ID "pulser"
#define MAPPING_FILE "tpcMapping.root"
#define CONFIG_FILE "TPCPULSERda.conf"
#define Q_FILE "tpcPulserQ.data"
#define DEAD_FILE "tpcDeadChannelsPulser.data"
#define AliDebugLevel() -1


#include <daqDA.h>
#include "event.h"
#include "monitor.h"
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

//
//Root includes
//
#include <TFile.h>
#include "TROOT.h"
#include "TPluginManager.h"
#include "TString.h"
#include "TObjString.h"
#include "TDatime.h"
//
//AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliTPCmapper.h"
#include "AliTPCRawStream.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliMathBase.h"
#include "TTreeStream.h"
#include "AliLog.h"
#include "TSystem.h"
#include "AliTPCConfigDA.h"
//
//AMORE
//
#include <AmoreDA.h>
//
// TPC calibration algorithm includes
//
#include "AliTPCCalibPulser.h"

/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  /* log start of process */
  printf("TPC Pulser DA started - %s\n",__FILE__);

  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  AliLog::SetClassDebugLevel("AliTPCRawStream",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");


  TString daterolename(gSystem->Getenv("DATE_ROLE_NAME"));
  if ( daterolename == "" ) {
    printf("Error: Variable DATE_ROLE_NAME not defined! Exiting ...\n");
    return -1;
  }
  bool inner;
  if      ( daterolename.EndsWith("-0") ) inner = true;
  else if ( daterolename.EndsWith("-1") ) inner = false;
  else {
    printf("Error: Variable DATE_ROLE_NAME neither ends with -0 nor -1 (E.g. ldc-TPC-C12-1)! Exiting ...\n");
    return -1;
  }

  /* declare monitoring program */
  int i,status;
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }

  // variables 
  AliTPCmapper *mapping = 0;   // The TPC mapping
  char localfile[255];
  unsigned long32 runNb=0;      //run number
  // configuration options 
  // if  test setup get parameters from $DAQDA_TEST_DIR 
   
  if (!mapping){
    /* copy locally the mapping file from daq detector config db */
    sprintf(localfile,"./%s",MAPPING_FILE);
    status = daqDA_DB_getFile(MAPPING_FILE,localfile);
    if (status) {
      printf("Failed to get mapping file (%s) from DAQdetDB, status=%d\n", MAPPING_FILE, status);
      return -1;
    }

    /* open the mapping file and retrieve mapping object */
    TFile *fileMapping = new TFile(MAPPING_FILE, "read");
    mapping = (AliTPCmapper*) fileMapping->Get("tpcMapping");
    delete fileMapping;
  }

  if (mapping == 0) {
    printf("Failed to get mapping object from %s.  ...\n", MAPPING_FILE);
    return -1;
  } else {
    printf("Got mapping object from %s\n", MAPPING_FILE);
  }

  //
  // DA configuration from configuration file
  //
 //retrieve configuration file
  sprintf(localfile,"./%s",CONFIG_FILE);
  status = daqDA_DB_getFile(CONFIG_FILE,localfile);
  if (status) {
    printf("Failed to get configuration file (%s) from DAQdetDB, status=%d\n", CONFIG_FILE, status);
    return -1;
  }
  AliTPCConfigDA config(CONFIG_FILE);
  // check configuration options
  Bool_t  skipAmore=kFALSE;

  //skip the amore part
  if ( config.GetConfigurationMap()->GetValue("SkipAmore") ) {
    skipAmore=((TObjString*)config.GetConfigurationMap()->GetValue("SkipAmore"))->GetString().Atoi();
    printf("TPCPULSERda: Skip Amore set in config\n");
  }

  // create calibration object
  AliTPCCalibPulser calibPulser(config.GetConfigurationMap());   // pulser calibration algorithm
  calibPulser.SetAltroMapping(mapping->GetAltroMapping()); // Use altro mapping we got from daqDetDb

  //===========================//
  // loop over RAW data files //
  //==========================//
  int nevents=0;
  for(i=1;i<argc;i++) {

    /* define data source : this is argument i */
    printf("Processing file %s\n", argv[i]);
    status=monitorSetDataSource( argv[i] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* read until EOF */
    while (true) {
      struct eventHeaderStruct *event;

      /* check shutdown condition */
      if (daqDA_checkShutdown()) {break;}

      /* get next event (blocking call until timeout) */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) {
        printf ("End of File %d detected\n",i);
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
      nevents++;
      // get the run number
      runNb = event->eventRunNb;
      //  Pulser calibration
      calibPulser.ProcessEvent(event);

      /* free resources */
      free(event);
    }
  }

  //
  // Analyse pulser data and write them to rootfile
  //
  calibPulser.Analyse();
  printf ("%d events processed\n",nevents);

  TFile * fileTPC = new TFile (RESULT_FILE,"recreate");
  calibPulser.Write("tpcCalibPulser");
  delete fileTPC;
  printf("Wrote %s\n",RESULT_FILE);

  /* store the result file on FES */

  status=daqDA_FES_storeFile(RESULT_FILE,FILE_ID);
  if (status) {
    status = -2;
  }
  //
  //Send objects to the AMORE DB
  //
  if (!skipAmore){
    printf ("AMORE part\n");
    const char *amoreDANameorig=gSystem->Getenv("AMORE_DA_NAME");
    //cheet a little -- temporary solution (hopefully)
    //
    //currently amoreDA uses the environment variable AMORE_DA_NAME to create the mysql
    //table in which the calib objects are stored. This table is dropped each time AmoreDA
    //is initialised. This of course makes a problem if we would like to store different
    //calibration entries in the AMORE DB. Therefore in each DA which writes to the AMORE DB
    //the AMORE_DA_NAME env variable is overwritten.

    //find processed sector
    Char_t sideName='A';
    Int_t sector = -1;
    for ( Int_t roc = 0; roc < 72; roc++ ) {
      if ( !calibPulser.GetCalRocT0(roc) ) continue;
      if (mapping->GetSideFromRoc(roc)==1) sideName='C';
      sector = mapping->GetSectorFromRoc(roc);
    }
    //gSystem->Setenv("AMORE_DA_NAME",Form("TPC-%c%02d-%s",sideName,sector,FILE_ID));
    gSystem->Setenv("AMORE_DA_NAME",Form("%s-%s", gSystem->Getenv("DATE_ROLE_NAME"), FILE_ID));
    //
    // end cheet
    if (sector>-1){
      TDatime time;
      TObjString info(Form("Run: %u; Date: %s",runNb,time.AsSQLString()));

      amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
      Int_t statusDA=0;
      statusDA+=amoreDA.Send("PulserT0",calibPulser.GetCalPadT0());
      statusDA+=amoreDA.Send("PulserQ",calibPulser.GetCalPadQ());
      statusDA+=amoreDA.Send("PulserRMS",calibPulser.GetCalPadRMS());
      statusDA+=amoreDA.Send("arrayTmean",calibPulser.GetMeanTimeSectorArray());
      statusDA+=amoreDA.Send("Info",&info);
      if ( statusDA!=0 )
        printf("Waring: Failed to write one of the calib objects to the AMORE database\n");
    } else {
      printf("Waring: No data found!\n");
    }
    // reset env var
    if (amoreDANameorig) gSystem->Setenv("AMORE_DA_NAME",amoreDANameorig);
  }
  
  //
  // Now prepare ASCII files for local ALTRO configuration through DDL. 
  //
  ofstream deadchannelfile;
  ofstream qfile;

  qfile.open(Q_FILE);
  deadchannelfile.open(DEAD_FILE);

  qfile           << 19 << std::endl; // Pulser Q
  deadchannelfile << 14 << std::endl; // Mark file to contain NOISY or DEAD CHANNELS 

  Int_t ctr_channel = 0;
  Int_t ctr_dead= 0;

  // inner==true : calROC from ldc-0 contains: rcus 0,1 for IROC and rcu 2 for OROC 
  // inner==false: calROC from ldc-1 contains: nothingfor IROC and rcus 3,4,5 for OROC
  for ( Int_t roc = 0; roc < 72; roc++ ) {
    if ( !calibPulser.GetCalRocQ(roc) ) continue;
    bool isIROC= mapping->IsIROC(roc);
    Int_t side = mapping->GetSideFromRoc(roc);
    Int_t sec= mapping->GetSectorFromRoc(roc);
    Int_t minrcu, maxrcu;
    if( isIROC ) { minrcu=0; maxrcu=1; }
    else if ( inner) { minrcu=2; maxrcu=2; }
    else { minrcu=3; maxrcu=5; }
    for ( int rcu = minrcu; rcu <= maxrcu; rcu++ ) {
      //Int_t patch = mapping->IsIROC(roc) ? rcu : rcu+2; 
      for ( int branch = 0; branch < 2; branch++ ) {
	for ( int fec = 0; fec < mapping->GetNfec(rcu, branch); fec++ ) {
	  for ( int altro = 0; altro < 8; altro++ ) {
	    for ( int channel = 0; channel < 16; channel++ ) {
	      Int_t hwadd = mapping->CodeHWAddress(branch, fec, altro, channel);
	      Int_t row = mapping->GetPadRow(rcu, hwadd);              // row in a ROC
	      Int_t globalrow = mapping->GetGlobalPadRow(rcu, hwadd);  // row in full sector
	      Int_t pad = mapping->GetPad(rcu, hwadd);
	      // skip edge pads
	      if ( (pad<1) || (pad>mapping->GetNpads(globalrow)-2) ) continue;
	      Float_t Q = calibPulser.GetCalRocQ(roc)->GetValue(row,pad);
	      qfile << ctr_channel++ << "\t" << side << "\t" << sec << "\t" << rcu << "\t" << hwadd << "\t" << Q << std::endl;
	      if ( fabs(Q) < 2) { // dead channel
		deadchannelfile << ctr_dead++ << "\t" << side << "\t" << sec << "\t"
				<< rcu << "\t" << hwadd << "\t" << Q << std::endl;
	      }
	    } // end channel for loop 
	  } // end altro for loop 
	} // end fec for loop 
      } // end branch for loop
    } // end rcu for loop 
  } // end roc loop 
  
  qfile.close();
  deadchannelfile.close();

  printf("Wrote ASCII file. Found %d dead channels.\n", ctr_dead);

return status;
}
