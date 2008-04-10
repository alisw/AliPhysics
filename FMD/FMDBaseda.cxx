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
#include <TSystem.h>
#include <AliCDBManager.h>
#include <AliFMDParameters.h>
#include <AliRawReader.h>
#include <TStopwatch.h>
#include <AliFMDBaseDA.h>
#include <AliRawReaderDate.h>
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
  
  
  Char_t* fileName = argv[1];
  
  Bool_t old = kTRUE;
  
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(0);
  cdb->SetDefaultStorage("");
  
  AliFMDParameters::Instance()->Init(kFALSE,0);
  AliFMDParameters::Instance()->UseRcuTrailer(!old);
  AliFMDParameters::Instance()->UseCompleteHeader(!old);
  
  AliRawReader *reader = new AliRawReaderDate(fileName);
  TStopwatch timer;
  timer.Start();
  AliFMDBaseDA baseDA;
  
  baseDA.Run(reader);
  
  timer.Stop();
  timer.Print();

  
  
  
}
