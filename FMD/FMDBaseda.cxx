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


  monitorDeclareTable(const_cast<char**>(tableSOD));

  
  Char_t* fileName = argv[1];
  
  Bool_t old = kTRUE;
  
  AliFMDParameters::Instance()->Init(kFALSE,0);
  AliFMDParameters::Instance()->UseCompleteHeader(!old);
  
  AliRawReader *reader = 0;
  TString fileNam(fileName);
  if (fileNam.EndsWith(".root")) reader = new AliRawReaderRoot(fileName);
  else reader = new AliRawReaderDate(fileName);
  if (!reader) { 
    std::cerr << "Don't know how to make reader for " << fileNam 
	      << std::endl;
    return -2;
  }
  
  TStopwatch timer;
  timer.Start();
  AliFMDBaseDA baseDA;
  
  baseDA.Run(reader);
  
  timer.Stop();
  timer.Print();

  Int_t  retval = daqDA_FES_storeFile("conditions.csv", AliFMDParameters::Instance()->GetConditionsShuttleID());
  if (retval != 0) std::cerr << "Base DA failed" << std::endl;
  
  return retval;
}
