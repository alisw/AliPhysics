/*

  FMD DA for online calibration of conditions

  Contact:                 canute@nbi.dk
  Link:                    fmd.nbi.dk/fmd/offline
  Run Type:                GAIN
  DA Type:                 LDC
  Number of events needed: usually 102400
  Input Files:             raw data 
  Output Files:            gains.csv
  Trigger types used:      GAIN
*/
#include <TSystem.h>
#include <TString.h>
#include <AliFMDParameters.h>
#include <AliRawReader.h>
#include <TStopwatch.h>
#include <AliFMDGainDA.h>
#include <AliRawReaderDate.h>
#include <AliRawReaderRoot.h>
#include <AliLog.h>
#include "daqDA.h"
#include "TROOT.h"
#include "TPluginManager.h"



int main(int argc, char **argv) 
{

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
  if(!secondArgument.IsWhitespace()&& !secondArgument.Contains("--help") 
     && !secondArgument.Contains("--diagnostics=true")) {
    std::cout<<"Second argument wrong. Use --help"<<std::endl;
    return -1;
  }
    
    Bool_t old = kTRUE;
    
  AliFMDParameters::Instance()->Init(kFALSE,0);

  //This will only work for FDR 1 data. When newer data becomes available the ! must be removed!
  AliFMDParameters::Instance()->UseCompleteHeader(old);
  
  AliLog::SetModuleDebugLevel("FMD", 1);
  
  AliRawReader *reader = 0;
  TString fileNam(fileName);
  if (fileNam.EndsWith(".root")) 
    reader = new AliRawReaderRoot(fileName);
  else reader = new AliRawReaderDate(fileName);
  if (!reader) { 
    std::cerr << "Don't know how to make reader for " << fileNam 
	      << std::endl;
    return -2;
  }

  
  TStopwatch timer;
  timer.Start();
  AliFMDGainDA gainDA;
  gainDA.SetSaveDiagnostics(diagnostics);
  gainDA.Run(reader);
  
  timer.Stop();
  timer.Print();
  
  Int_t  retvalConditions = daqDA_FES_storeFile("conditions.csv", AliFMDParameters::Instance()->GetConditionsShuttleID());
  Int_t  retvalGain = daqDA_FES_storeFile("gains.csv", AliFMDParameters::Instance()->GetGainShuttleID());

  if(retvalConditions!=0 || retvalGain!=0)
    std::cerr << "Pedestal DA failed" << std::endl;
  
  if(retvalGain != 0) return retvalGain;
  else return retvalConditions;
  
  
  
  
}
