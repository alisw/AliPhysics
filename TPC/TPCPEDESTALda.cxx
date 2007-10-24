/*

TPCda_pedestal.cxx - calibration algorithm for TPC pedestal runs

10/06/2007  sylvain.chapeland@cern.ch :  first version - clean skeleton based on DAQ DA case1

contact: marian.ivanov@cern.ch


This process reads RAW data from the files provided as command line arguments
and save results in a file (named from RESULT_FILE define - see below).

*/

#define RESULT_FILE "tpcPedestal.root"


extern "C" {
#include <daqDA.h>
}
#include "event.h"
#include "monitor.h"
#include <stdio.h>
#include <stdlib.h>
#include <fstream.h>

//
//Root includes
//
#include <TFile.h>
#include <TArrayF.h>

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
  Bool_t timeAnalysis = kTRUE;

  int i,status;
  AliTPCCalibPedestal calibPedestal;           // pedestal and noise calibration
  calibPedestal.SetTimeAnalysis(timeAnalysis); // pedestal(t) calibration

  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  /* log start of process */
  printf("TPC DA started - %s\n",__FILE__);

  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }

  /* loop over RAW data files */
  int nevents=0;
  for ( i=1; i<argc; i++ ) {

    /* define data source : this is argument i */
    printf("Processing file %s\n", argv[i]);
    status=monitorSetDataSource( argv[i] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
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

      /* skip start/end of run events */
      if ( (event->eventType != physicsEvent) && (event->eventType != calibrationEvent) )
	continue;

      /* retry if got no event */
      if (event==NULL) {
	continue;
      }

      nevents++;

      //  Pedestal calibration
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      calibPedestal.ProcessEvent(rawReader);
      delete rawReader;

      /* free resources */
      free(event);
    }
  }

  calibPedestal.Analyse(); 
  calibPedestal.AnalyseTime(nevents); 
  printf ("%d physics/calibration events processed\n",nevents);

  TFile *fileTPC = new TFile(RESULT_FILE, "recreate");
  calibPedestal.Write("calibPedestal");
  delete fileTPC;
  printf("Wrote %s\n",RESULT_FILE);

  //
  // Now prepare ASCII files for local ALTRO configuration through DDL
  //

  ofstream pedfile;
  ofstream noisefile;
  ofstream pedmemfile;
  char filename[255];
  sprintf(filename,"tpcPedestals.data");
  pedfile.open(filename);
  sprintf(filename,"tpcNoise.data");
  noisefile.open(filename);
  sprintf(filename,"tpcPedestalMem.data");
  pedmemfile.open(filename);

  TArrayF **timePed = calibPedestal.GetTimePedestals();
  AliTPCmapper mapping;

  Int_t ctr_channel = 0;
  Int_t ctr_altro = 0;
  Int_t ctr_pattern = 0;

  pedfile    << 10 << std::endl; // PEDESTALS per channel
  noisefile  << 11 << std::endl; // NOISE per altro
  pedmemfile << 12 << std::endl; // PEDESTALs per time bin

  for ( Int_t roc = 0; roc < 72; roc++ ) {
    if ( !calibPedestal.GetCalRocPedestal(roc) ) continue;
    Int_t side   = mapping.GetSideFromRoc(roc);
    Int_t sector = mapping.GetSectorFromRoc(roc);
    //printf("Analysing ROC %d (side %d, sector %d) ...\n", roc, side, sector);
    Int_t nru = mapping.IsIROC(roc) ? 2 : 4;
    for ( int rcu = 0; rcu < nru; rcu++ ) {
      Int_t patch = mapping.IsIROC(roc) ? rcu : rcu+2;
      for ( int branch = 0; branch < 2; branch++ ) {
	for ( int fec = 0; fec < mapping.GetNfec(patch, branch); fec++ ) {
	  for ( int altro = 0; altro < 8; altro++ ) {
	    Float_t rms = 0.;
	    Float_t ctr = 0.;
	    for ( int channel = 0; channel < 16; channel++ ) {
	      Int_t hwadd = mapping.CodeHWAddress(branch, fec, altro, channel);
	      Int_t row   = mapping.GetPadRow(patch, hwadd);
	      Int_t pad   = mapping.GetPad(patch, hwadd);
	      Float_t ped = calibPedestal.GetCalRocPedestal(roc)->GetValue(row,pad);
	      // fixed pedestal
	      if ( ped > 1.e-10 ) {
		pedfile << ctr_channel << "\t" << side << "\t" << sector << "\t" << patch << "\t" << hwadd << "\t" << ped << std::endl;
		ctr_channel++;
	      }
	      // pedestal(t)
	      if ( timePed && fabs(timePed[row][pad].GetSum()) > 1e-10 ) {
		pedmemfile << ctr_pattern << "\t" << side << "\t" << sector << "\t" << patch << "\t" << hwadd;
		for ( Int_t timebin = 0; timebin < 1024; timebin++ )
		  pedmemfile << "\t" << timePed[row][pad].At(timebin);
		pedmemfile << std::endl;
		ctr_pattern++;
	      }
	      // rms=noise
	      Float_t rms2 = calibPedestal.GetCalRocRMS(roc)->GetValue(row,pad);
	      if ( rms2 > 1.e-10 ) { rms += rms2; ctr += 1.; }
	    } // end channel for loop
	    // noise data (rms) averaged over all channels in this ALTRO.
	    Int_t hwadd = mapping.CodeHWAddress(branch, fec, altro, 0);
	    if ( ctr > 1.e-10 ) {
	      noisefile << ctr_altro << "\t" << side << "\t" << sector << "\t" << patch << "\t" << hwadd << "\t" << rms/ctr << std::endl;
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

  printf("Wrote ASCII files\n");

  return status;
}
