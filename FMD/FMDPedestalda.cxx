/*

  FMD DA for online calibration of conditions

  Contact:                 canute@nbi.dk
  Link:                    fmd.nbi.dk/fmd/offline
  Run Type:                PEDESTAL
  DA Type:                 LDC
  Number of events needed: 1000
  Input Files:             raw data 
  Output Files:            peds.csv
  Trigger types used:      PEDESTAL
*/
#include <iostream>
#include <TSystem.h>
#include <TString.h>
#include <AliFMDParameters.h>
#include <AliRawReader.h>
#include <TStopwatch.h>
#include <AliFMDPedestalDA.h>
#include <AliRawReaderDate.h>
#include <AliRawReaderRoot.h>
#include "daqDA.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include <AliLog.h>


int main(int argc, char **argv) 
{
  /* magic line from Rene - for future reference! */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
  					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");
  
  
  Bool_t diagnostics = kFALSE;
  Char_t* fileName = argv[1];
  TString secondArgument(argv[2]);
  
  if(secondArgument.Contains("--diagnostics=true"))
    diagnostics = kTRUE;
  if(secondArgument.Contains("--help")) {
    std::cout<<"Usage: filename --diagnostics=true/false . --help this help"<<std::endl;
    return 0;
  }
  
  if(!secondArgument.IsWhitespace() && !secondArgument.Contains("--help") 
     && !secondArgument.Contains("--diagnostics=true")) {
    std::cout<<"Second argument wrong. Use --help"<<std::endl;
    return -1;
  }
  Bool_t old = kTRUE;

  AliLog::EnableDebug(kFALSE);
  AliFMDParameters::Instance()->Init(kFALSE,0);
  AliFMDParameters::Instance()->UseCompleteHeader(old);

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
  AliFMDPedestalDA pedDA;
  pedDA.SetSaveDiagnostics(diagnostics);
  pedDA.Run(reader);
  
  timer.Stop();
  timer.Print();
 
  Int_t  retvalConditions = daqDA_FES_storeFile("conditions.csv", AliFMDParameters::Instance()->GetConditionsShuttleID());
  Int_t  retvalPeds = daqDA_FES_storeFile("peds.csv", AliFMDParameters::Instance()->GetPedestalShuttleID());

  if(retvalConditions!=0 || retvalPeds!=0)
    std::cerr << "Pedestal DA failed" << std::endl;
  
  if(retvalPeds != 0) return retvalPeds;
  else return retvalConditions;
  
}
