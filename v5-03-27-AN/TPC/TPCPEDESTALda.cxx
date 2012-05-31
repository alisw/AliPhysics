/*
TPC DA for online calibration

Contact: Haavard.Helstrup@cern.ch
Link:
Run Type: PEDESTAL
DA Type: LDC
Number of events needed: 100
Input Files: 
Output Files: tpcPedestal.root, to be exported to the DAQ FXS
fileId:   pedestals    
Trigger types used: CALIBRATION_EVENT

*/

/*

TPCda_pedestal.cxx - calibration algorithm for TPC pedestal runs

10/06/2007  sylvain.chapeland@cern.ch :  first version - clean skeleton based on DAQ DA case1
19/10/2007  christian.lippmann@cern.ch :  Possibility to write output to ASCII file
24/10/2007  christian.lippmann@cern.ch :  Including pedestal calibration for time bins
23/11/2007  christian.lippmann@cern.ch :  Fix in order to avoid streamer problems in case of
                                          invalid ROOTSTYS. The famous magic line provided by Rene.
28/11/2007  christian.lippmann@cern.ch :  TPC mapping file is read from DaqDetDB
18/09/2008  christian.lippmann@cern.ch :  Noisy channels are output to ASCII file. Use max noise in ALTRO.
19/09/2008  J.Wiechula@gsi.de:            Added export of the calibration data to the AMORE data base.
                                          Added support for configuration files.
31/01/2011  Christian.Lippmann@cern.ch :  Updates for changed setup at P2 with 2 LDCs per sector

contact: marian.ivanov@cern.ch

This process reads RAW data from the files provided as command line arguments
and save results in a file (named from RESULT_FILE define - see below).

*/

#define RESULT_FILE  "tpcPedestal.root"
#define FILE_ID "pedestals"
#define MAPPING_FILE "tpcMapping.root"
#define CONFIG_FILE "TPCPEDESTALda.conf"
#define PED_FILE "tpcPedestals.data"
#define NOISE_FILE "tpcNoise.data"
#define PEDMEM_FILE "tpcPedestalMem.data"
#define NOISY_FILE "tpcNoisyChannels.data"
#define VERYNOISY_FILE "tpcVeryNoisyChannels.data"
#define DEAD_FILE "tpcDeadChannels.data"
#define AliDebugLevel() -1

extern "C" {
#include <daqDA.h>
}
#include "event.h"
#include "monitor.h"

#include "stdio.h"
#include "stdlib.h"
#include <fstream>

//
//Root includes
//
#include "TFile.h"
#include "TArrayF.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TSystem.h"
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
#include "AliTPCConfigDA.h"

//
//AMORE
//
#include <AmoreDA.h>

//
// TPC calibration algorithm includes
//
#include "AliTPCCalibPedestal.h"

/*
  Main routine, TPC pedestal detector algorithm to be run on TPC LDC
  Arguments: list of DATE raw data files
*/

int main(int argc, char **argv) {
  //
  // Main for TPC pedestal detector algorithm
  //
  /* log start of process */
  printf("TPC DA started - %s\n",__FILE__);
  
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

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

  AliLog::SetClassDebugLevel("AliTPCRawStream",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);

  /* magic line */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");

  /* declare monitoring program */
  int i, status;
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }
  
  // variables
  AliTPCmapper *mapping = 0;   // The TPC mapping
  char localfile[255];
  unsigned long32 runNb=0;     // run number
  // configuration options 
  Bool_t timeAnalysis = kTRUE;
  Bool_t fastDecoding = kFALSE;

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
  // retrieve configuration file
  sprintf(localfile,"./%s",CONFIG_FILE);
  status = daqDA_DB_getFile(CONFIG_FILE,localfile);
  if (status) {
    printf("Failed to get configuration file (%s) from DAQdetDB, status=%d\n", CONFIG_FILE, status);
    return -1;
  }
  AliTPCConfigDA config(CONFIG_FILE);
  // check configuration
  Bool_t  skipAmore=kFALSE;

  if ( (Int_t)config.GetValue("NoTimeAnalysis") == 1 ) {
    printf("WARNING: Time analysis was switched off in the configuration file!\n");
    timeAnalysis=kFALSE;
  }

  if ( (Int_t)config.GetValue("UseFastDecoder") == 1 ){
    printf("Info: The fast decoder will be used for the processing.\n");
    fastDecoding=kTRUE;
  }

  if ( config.GetConfigurationMap()->GetValue("SkipAmore") ) {
    skipAmore=((TObjString*)config.GetConfigurationMap()->GetValue("SkipAmore"))->GetString().Atoi();
    printf("TPCPEDESTALda: Skip Amore set in config\n");
  }


  // create calibration object
  AliTPCCalibPedestal calibPedestal(config.GetConfigurationMap()); // pedestal and noise calibration
  calibPedestal.SetAltroMapping(mapping->GetAltroMapping()); // Use altro mapping we got from daqDetDb
  calibPedestal.SetTimeAnalysis(timeAnalysis);               // pedestal(t) calibration 
  
  //===========================//
  // loop over RAW data files //
  //==========================//
  int nevents=0;
  for ( i=1; i<argc; i++ ) {
    
    /* define data source : this is argument i */
    printf("Processing file %s\n", argv[i]);
    status=monitorSetDataSource( argv[i] );
    if (status!=0) {
      printf("monitorSetDataSource() failed. Error=%s. Exiting ...\n", monitorDecodeError(status));
      return -1;
    }
    
    /* read until EOF */
    for ( ; ; ) {
      struct eventHeaderStruct *event;
      
      /* check shutdown condition */
      if (daqDA_checkShutdown()) {break;}
      
      /* get next event (blocking call until timeout) */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) {
        printf ("End of File %d (%s) detected\n", i, argv[i]);
        break; /* end of monitoring file has been reached */
      }
      
      if (status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        break;
      }

      /* retry if got no event */
      if (event==NULL)
        continue;
      
      /* skip start/end of run events */
      if ( (event->eventType != physicsEvent) && (event->eventType != calibrationEvent) )
        continue;

      nevents++;
      // get the run number
      runNb = event->eventRunNb;
      //  Pedestal calibration
      calibPedestal.ProcessEvent(event);

      /* free resources */
      free(event);
    }
  }

  //
  // Analyse pedestals and write them to rootfile
  //
  calibPedestal.Analyse();
  calibPedestal.AnalyseTime(nevents);
  printf ("%d physics/calibration events processed.\n",nevents);

  TFile *fileTPC = new TFile(RESULT_FILE, "recreate");
  calibPedestal.Write("tpcCalibPedestal");
  delete fileTPC;
  printf("Wrote %s.\n",RESULT_FILE);

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
      if ( !calibPedestal.GetCalRocPedestal(roc) ) continue;
      if (mapping->GetSideFromRoc(roc)==1) sideName='C';
      sector = mapping->GetSectorFromRoc(roc);
    }
  //   gSystem->Setenv("AMORE_DA_NAME",Form("TPC-%c%02d-%s",sideName,sector,FILE_ID));
    gSystem->Setenv("AMORE_DA_NAME",Form("%s-%s",gSystem->Getenv("DATE_ROLE_NAME"),FILE_ID));

    //
    // end cheet
    if (sector>-1){
      TDatime time;
      TObjString info(Form("Run: %u; Date: %s",runNb,time.AsSQLString()));

      amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
      Int_t statusDA=0;
      statusDA+=amoreDA.Send("Pedestals",calibPedestal.GetCalPadPedestal());
      statusDA+=amoreDA.Send("Noise",calibPedestal.GetCalPadRMS());
      statusDA+=amoreDA.Send("Info",&info);
      if ( statusDA )
        printf("Warning: Failed to write one of the calib objects to the AMORE database\n");
    }  else {
      printf("Warning: No data found!\n");
    }
    // reset env var
    if (amoreDANameorig) gSystem->Setenv("AMORE_DA_NAME",amoreDANameorig);
  }
  
  //
  // Now prepare ASCII files for local ALTRO configuration through DDL.
  //
  ofstream pedfile;
  ofstream noisefile;
  ofstream pedmemfile;
  ofstream noisychannelfile;
  ofstream verynoisychannelfile;
  ofstream deadchannelfile;

  pedfile.open(PED_FILE);
  noisefile.open(NOISE_FILE);
  pedmemfile.open(PEDMEM_FILE);
  noisychannelfile.open(NOISY_FILE);
  verynoisychannelfile.open(VERYNOISY_FILE);
  deadchannelfile.open(DEAD_FILE);

  TArrayF **timePed = calibPedestal.GetTimePedestals();  // pedestal values for each time bin

  Int_t ctr_channel = 0;
  Int_t ctr_altro   = 0;
  Int_t ctr_pattern = 0;
  Int_t ctr_noisy   = 0;
  Int_t ctr_vnoisy  = 0;
  Int_t ctr_dead    = 0;

  pedfile              << 10 << std::endl; // Mark file to contain PEDESTALS per channel
  noisefile            << 11 << std::endl; // Mark file to contain NOISE per altro
  pedmemfile           << 12 << std::endl; // Mark file to contain PEDESTALs per time bin
  noisychannelfile     << 14 << std::endl; // Mark file to contain NOISY or DEAD CHANNELS
  verynoisychannelfile << 14 << std::endl; // Mark file to contain NOISY or DEAD CHANNELS
  deadchannelfile      << 14 << std::endl; // Mark file to contain NOISY or DEAD CHANNELS

  // inner==true : calROC from ldc-0 contains: rcus 0,1 for IROC and rcu 2 for OROC 
  // inner==false: calROC from ldc-1 contains: nothing  for IROC and rcus 3,4,5 for OROC 
  for ( Int_t roc = 0; roc < 72; roc++ ) {
    if ( !calibPedestal.GetCalRocPedestal(roc) ) continue;
    bool isIROC  = mapping->IsIROC(roc);
    Int_t side   = mapping->GetSideFromRoc(roc);
    Int_t sector = mapping->GetSectorFromRoc(roc);
    Int_t minrcu, maxrcu;
    if      ( isIROC ) { minrcu=0; maxrcu=1; }
    else if ( inner  ) { minrcu=2; maxrcu=2; }
    else               { minrcu=3; maxrcu=5; }
    for ( int rcu = minrcu; rcu <= maxrcu; rcu++ ) {
      //Int_t patch = mapping->IsIROC(roc) ? rcu : rcu+2;
      for ( int branch = 0; branch < 2; branch++ ) {
        for ( int fec = 0; fec < mapping->GetNfec(rcu, branch); fec++ ) {
          for ( int altro = 0; altro < 8; altro++ ) {
            Float_t rms = 0.;
            Float_t maxrms = 0.;
            Float_t ctr_altrochannel = 0.;
            for ( int channel = 0; channel < 16; channel++ ) {
              Int_t hwadd     = mapping->CodeHWAddress(branch, fec, altro, channel);
              Int_t row       = mapping->GetPadRow(rcu, hwadd);        // row in a ROC
              Int_t globalrow = mapping->GetGlobalPadRow(rcu, hwadd);  // row in full sector
              Int_t pad       = mapping->GetPad(rcu, hwadd);
              Float_t ped     = calibPedestal.GetCalRocPedestal(roc)->GetValue(row,pad);
              // fixed pedestal
              pedfile << ctr_channel++ << "\t" << side << "\t" << sector << "\t" << rcu << "\t"
		      << hwadd << "\t" << ped << std::endl;
              // pedestal(t)=pedestal memories
              if ( timePed && fabs(timePed[globalrow][pad].GetSum()) > 1e-10 ) {
                pedmemfile << ctr_pattern++ << "\t" << side << "\t" << sector << "\t" << rcu
                    << "\t" << hwadd;
                for ( Int_t timebin = 0; timebin < 1024; timebin++ )
                  pedmemfile << "\t" << timePed[globalrow][pad].At(timebin);
                pedmemfile << std::endl;
              }
              // rms=noise
              Float_t rms2 = calibPedestal.GetCalRocRMS(roc)->GetValue(row,pad);
              if ( fabs(ped) < 1.e-10 ) {                        // dead channel
                deadchannelfile << ctr_dead++ << "\t" << side << "\t" << sector << "\t"
				<< rcu << "\t" << hwadd << "\t" << rms2 << std::endl;
              } else if ( (ped > 1.e-10) && (rms2 > 1.e-10) ) {  // not dead
		// Find noisy channels
                if ( rms2 > 6.0 ) { // VERY noisy
                  verynoisychannelfile << ctr_vnoisy++ << "\t" << side << "\t" << sector << "\t"
				       << rcu << "\t" << hwadd << "\t" << rms2 << std::endl;
		  
		} else if ( ((roc<36)             && (rms2 > 2.0))  ||  // IROC
			    ((roc>35) && (row<65) && (rms2 > 2.0))  ||  // OROC, small pads
			    ((roc>35) && (row>64) && (rms2 > 3.0)) ) {  // OROC, large pads (50% more signal)
                  noisychannelfile << ctr_noisy++ << "\t" << side << "\t" << sector << "\t"
				   << rcu << "\t" << hwadd << "\t" << rms2 << std::endl;
		} else {
		  // Not noisy. Get average and maximum noise in this ALTRO
		  rms += rms2;
		  ctr_altrochannel += 1.;
		  if (rms2 > maxrms) maxrms = rms2;
		} // end if noisy
              } // end if some signal
            } // end channel for loop
            Int_t hwadd = mapping->CodeHWAddress(branch, fec, altro, 0);   // ALTRO address
	    // Noise data (rms) averaged over all channels in this ALTRO.
            if ( ctr_altrochannel > 1.e-10 ) {
	      /*
	      // average noise of this ALTRO (excluding high-noise channels)
              noisefile << ctr_altro << "\t" << side << "\t" << sector << "\t" << rcu << "\t"
              << hwadd << "\t" << rms/ctr_altrochannel << std::endl;
	      */
	      // maximum noise of this ALTRO (excluding high-noise channels)
              noisefile << ctr_altro << "\t" << side << "\t" << sector << "\t" << rcu << "\t"
			<< hwadd << "\t" << maxrms << std::endl;
              ctr_altro++;
            }
          } // end altro for loop
        } // end fec for loop
      } // end branch for loop
    } // end rcu for loop
  } // end roc loop
  
  pedfile.close();
  noisefile.close();
  pedmemfile.close();
  noisychannelfile.close();
  verynoisychannelfile.close();
  deadchannelfile.close();
  printf("Wrote ASCII files. Found %d noisy, %d very noisy and %d dead channels.\n", ctr_noisy, ctr_vnoisy, ctr_dead);

  return status;

}
