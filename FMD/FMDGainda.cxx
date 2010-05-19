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
#ifdef ALI_AMORE
# include <AmoreDA.h>
# include <TH2.h>
#endif



int main(int argc, char **argv) 
{

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
  Char_t* fileName = argv[1];
  TString secondArgument(argv[2]);
  
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
    
  AliFMDParameters::Instance()->Init(kFALSE,0);

  //This will only work for FDR 1 data. When newer data becomes available the ! must be removed!
  AliFMDParameters::Instance()->UseCompleteHeader(old);
  
  AliLog::EnableDebug(kFALSE);
  
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
#ifdef ALI_AMORE
  gainDA.SetMakeSummaries(kTRUE);
#endif
  gainDA.Run(reader);
  
  timer.Stop();
  timer.Print();
  
  Int_t  retvalConditions = 
    daqDA_FES_storeFile("conditions.csv", 
			AliFMDParameters::Instance()->GetConditionsShuttleID());
  Int_t  retvalGain = 
    daqDA_FES_storeFile("gains.csv", 
			AliFMDParameters::Instance()->GetGainShuttleID());

  if(retvalConditions!=0 || retvalGain!=0)
    std::cerr << "Pedestal DA failed" << std::endl;
  
#ifdef ALI_AMORE
  try { 
    amore::da::AmoreDA myAmore(amore::da::AmoreDA::kSender);

    UShort_t det = 0;
    for (det = 1; det <= 3; det++) 
      if (gainDA.HasSeenDetector(det)) break;
    if (det >= 1 && det <= 3) { 
      TObject* runNo = new TObject;
      runNo->SetUniqueID(reader->GetRunNumber());
      myAmore.Send(Form("gainRunNoFMD%d", det), runNo);
    }
		   
    TIter     next(&gainDA.GetSummaries());
    TObject*  obj = 0;
    while ((obj = next())) 
      myAmore.Send(obj->GetName(), obj);
    
  }
  catch (std::exception& e) {
    std::cerr << "Failed to make AMORE instance: " << e.what() << std::endl;
  }
			       
#endif

  if(retvalGain != 0) return retvalGain;
  return retvalConditions;

}
//
// EOF
//

