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
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", "Minuit", 
					"TMinuitMinimizer",
					"Minuit", 
					"TMinuitMinimizer(const char *)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"GSLMultiMin", 
					"ROOT::Math::GSLMinimizer",
					"MathMore", 
					"GSLMinimizer(const char *)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"GSLMultiFit", 
					"ROOT::Math::GSLNLSMinimizer",
					"MathMore", "GSLNLSMinimizer(int)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"GSLSimAn", 
					"ROOT::Math::GSLSimAnMinimizer",
					"MathMore", 
					"GSLSimAnMinimizer(int)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"Linear", 
					"TLinearMinimizer",
					"Minuit", 
					"TLinearMinimizer(const char *)");
  gROOT->GetPluginManager()->AddHandler("ROOT::Math::Minimizer", 
					"Fumili", 
					"TFumiliMinimizer",
					"Fumili", 
					"TFumiliMinimizer(int)");


  
  
  Bool_t diagnostics = kFALSE;
  if (argc < 2) { 
    std::cerr << "No input file given" << std::endl;
    return 1;
  }
  TString fileName(argv[1]);
  if (fileName.Contains("--help")) { 
    std::cout << "Usage: " << argv[0] << " FILENAME [OPTIONS]\n\n"
	      << "Options:\n" 
	      << "    --diagnostics=BOOL Make diagnostics ROOT file\n"
	      << std::endl;
    return 0;
  }
      
  for (int i = 2; i < argc; i++) { 
    TString arg(argv[i]);
    if      (arg.Contains("--diagnostics=true")) diagnostics = kTRUE;
    else if (arg.Contains("--help")) { 
      std::cout << "Usage: " << argv[0] << " FILENAME [OPTIONS]\n\n"
		<< "Options:\n" 
		<< "    --diagnostics=BOOL Make diagnostics ROOT file\n"
		<< std::endl;
      return 0;
    }
    else { 
      std::cerr << "Unknown option: " << arg << "\n"
		<< "Try '" << argv[0] << " --help" << std::endl;
      return 1;
    }
  }
  Bool_t old = kTRUE;

  AliLog::EnableDebug(kFALSE);
  AliFMDParameters::Instance()->Init(kFALSE,0);
  AliFMDParameters::Instance()->UseCompleteHeader(old);

  AliRawReader *reader = 0;
  if (fileName.EndsWith(".root")) 
    reader = new AliRawReaderRoot(fileName.Data());
  else 
    reader = new AliRawReaderDate(fileName.Data());
  if (!reader) { 
    std::cerr << "Don't know how to make reader for " << fileName
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
 
  Int_t  retvalConditions = 
    daqDA_FES_storeFile("conditions.csv", 
			AliFMDParameters::Instance()->GetConditionsShuttleID());
  Int_t  retvalPeds = 
    daqDA_FES_storeFile("peds.csv", 
			AliFMDParameters::Instance()->GetPedestalShuttleID());

  if(retvalConditions!=0 || retvalPeds!=0)
    std::cerr << "Pedestal DA failed" << std::endl;
  
  if(retvalPeds != 0) return retvalPeds;
  else return retvalConditions;
  
}
