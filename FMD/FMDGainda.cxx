/*

  FMD DA for online calibration of conditions

  Contact:                 canute@nbi.dk
  Link:                    fmd.nbi.dk/fmd/offline
  Run Type:                GAIN
  DA Type:                 GAIN
  Number of events needed: usually 102400
  Input Files:             raw data 
  Output Files:            gains.csv
  Trigger types used:      GAIN
*/
#include <TSystem.h>
#include <AliFMDParameters.h>
#include <AliRawReader.h>
#include <TStopwatch.h>
#include <AliFMDGainDA.h>
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
    
  AliFMDParameters::Instance()->Init(kFALSE,0);
  AliFMDParameters::Instance()->SetSampleRate(4);
  AliFMDParameters::Instance()->UseRcuTrailer(!old);

  //This will only work for FDR 1 data. When newer data becomes available the ! must be removed!
  AliFMDParameters::Instance()->UseCompleteHeader(!old);
  
  
  
  
  AliRawReader *reader = new AliRawReaderDate(fileName);
  TStopwatch timer;
  timer.Start();
  AliFMDGainDA gainDA;
  
  gainDA.Run(reader);
  
  timer.Stop();
  timer.Print();

  
  
  
}
