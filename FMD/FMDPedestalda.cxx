/*

  FMD DA for online calibration of conditions

  Contact:                 canute@nbi.dk
  Link:                    fmd.nbi.dk/fmd/offline
  Run Type:                PEDESTAL
  DA Type:                 Pedestal
  Number of events needed: 1000
  Input Files:             raw data 
  Output Files:            peds.csv
  Trigger types used:      PEDESTAL
*/
#include <TSystem.h>
#include <AliFMDParameters.h>
#include <AliRawReader.h>
#include <TStopwatch.h>
#include <AliFMDPedestalDA.h>
#include <AliRawReaderDate.h>
#include "TROOT.h"
#include "TPluginManager.h"



int main(int argc, char **argv) 
{

  //#if 0
  /* magic line from Rene - for future reference! */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");
  //#endif
  
  
  Char_t* fileName = argv[1];
  
  Bool_t old = kTRUE;


  AliFMDParameters::Instance()->Init(kFALSE,0);
  AliFMDParameters::Instance()->SetSampleRate(4);
  AliFMDParameters::Instance()->UseRcuTrailer(!old);
  AliFMDParameters::Instance()->UseCompleteHeader(old);

  AliRawReader *reader = new AliRawReaderDate(fileName);
  TStopwatch timer;
  timer.Start();
  AliFMDPedestalDA pedDA;
  //pedDA.SetSaveDiagnostics(kTRUE);
  pedDA.Run(reader);
  
  timer.Stop();
  timer.Print();

  
  
  
}
