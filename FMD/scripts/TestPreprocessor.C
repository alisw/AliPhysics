/* $Id$ */

//====================================================================
//
// Helper classes
//
#include <iostream>
#include <fstream>
#ifndef __CINT__
# include <TString.h>
# include <AliPreprocessor.h>
# include <FMD/AliFMDPreprocessor.h>
# include <AliShuttleInterface.h>
# include <AliCDBStorage.h>
# include <AliCDBEntry.h>
# include <AliCDBManager.h>
# include <AliCDBId.h>
# include <AliCDBMetaData.h>
# include <AliDCSValue.h>
# include <AliLog.h>
// # include <../FMD/AliFMDCalibZeroSuppression.h>
// # include <../FMD/AliFMDCalibDeadMap.h>
# include <FMD/AliFMDParameters.h>
# include <SHUTTLE/TestShuttle/AliTestShuttle.h>
# include <TRandom.h>
# include <TSystem.h>
# include <TObjString.h>
# include <TMap.h>
# include <TError.h>
# include <iostream>
#endif

namespace { 
  //====================================================================
  //
  // Helper functions
  //
  Float_t Hardware2Ped(int ddl, int board, int chip, int, int) 
  {
    return ((chip & 0xf) | ((board & 0x1F) << 4) | (ddl << 9) & 0x3);
  }
  Float_t Hardware2Noise(int ddl, int board, int chip, int channel, int strip)
  {
    return ((strip & 0x7f) | ((channel & 0xf) << 7));
  }
  Float_t Hardware2Gain(int ddl, int board, int chip, int channel, int strip) 
  {
    return (((strip   & 0x7f) <<  0) | 
	    ((channel & 0x0f) <<  7) | 
	    ((chip    & 0x07) << 11) |
	    ((board   & 0x1f) << 14) | 
	    ((ddl     & 0x03) << 19));
  }

  //====================================================================
  //
  // Helper classes
  //
  //__________________________________________________________________
  class CreateDummyDaData
  {
  public:
    CreateDummyDaData(const char* output,
		      int firstDDL,   int lastDDL, 
		      int firstStrip, int lastStrip)
      : fOutput(output),
	fFirstDDL(firstDDL),     fLastDDL(lastDDL), 
	fFirstStrip(firstStrip), fLastStrip(lastStrip)
    {}
    void Exec()
    {
      std::cout << "Will write on " << fOutput << std::endl;
      std::ofstream file(fOutput.Data());
      if (file.bad()) { 
	std::cerr << "Failed to open output file " << fOutput << std::endl;
	return;
      }
      Header(file);
      for (int ddl = fFirstDDL; ddl <= fLastDDL; ddl++) { 
	int  boards[] = { 0, 16, (ddl==1 ? -1 : 1), (ddl==1 ? -1 : 17), -1};
	int* bptr     = boards;
	int  board    = -1;
	while ((board = (*bptr++)) >= 0) { 
	  for (int chip = 0; chip < 3; chip++) { 
	    for (int channel = 0; channel < (chip == 1 ? 8 : 16); channel++) { 
	      for (int strip = fFirstStrip; strip <= fLastStrip; strip++) {
		Output(file, ddl, board, chip, channel, strip);
	      }
	    } // for channel 
	  } // for chip
	} // while board 
      } // for ddl
      file.close();
    }
    virtual void Header(std::ostream& file) = 0;
    virtual void Output(std::ostream& file, int ddl, int board, int chip, 
			int channel, int strip) = 0;
  protected:
    TString fOutput;
    int fFirstDDL;
    int fLastDDL;
    int fFirstStrip;
    int fLastStrip;
  };

  //__________________________________________________________________
  class CreateDummyPeds : public  CreateDummyDaData
  {
  public:
    CreateDummyPeds(const char* out="peds.csv", 
		    int overSampling=4,
		    int firstDDL=0,   int lastDDL=2, 
		    int firstStrip=0, int lastStrip=127)
      : CreateDummyDaData(out, firstDDL, lastDDL, firstStrip, lastStrip), 
	fOverSampling(overSampling)
    {}
    void Output(std::ostream& file, int ddl, int board, int chip, int channel, 
		int strip)
    {
      // Format is
      //   ddl,board,chip,channel,strip,sample,ped,noise,mu,sigma,chi2
      for (int sample = 0; sample < fOverSampling; sample++) {
	Float_t ped   = Hardware2Ped(ddl, board, chip, channel, strip);
	Float_t noise = Hardware2Noise(ddl, board, chip, channel, strip);
	file // << ddl     << ","
	     << board   << "," 
	     << chip    << ","
	     << channel << ","
	     << strip   << ","
	     << sample  << ","
	     << ped     << ","  
	     << noise   << ","  
	     << chip    << ","  // Predictable mu
	     << channel << ","  // Predictable sigma 
	     << strip           // Predictable chi2/ndf
	     << endl;
      }
    }
    void Header(std::ostream& file) 
    { 
      file << "# Pedestals\n" 
	   << "# ddl,board,chip,channel,strip,sample,mean,noise,mu,sigma,chi"
	   << std::endl;
    }
  protected:
    int fOverSampling;
  };

  //__________________________________________________________________
  class CreateDummyGains : public  CreateDummyDaData
  {
  public:
    CreateDummyGains(const char* out="gains.csv", 
		     int firstDDL=0,   int lastDDL=2, 
		     int firstStrip=0, int lastStrip=127)
      : CreateDummyDaData(out, firstDDL, lastDDL, firstStrip, lastStrip) 
    {}
    void Output(std::ostream& file, int ddl, int board, int chip, int channel, 
		int strip)
    {
      // Format is
      //   ddl,board,chip,channel,strip,gain,error,chi2
      Float_t gain = Hardware2Gain(ddl, board, chip, channel, strip);
      file // << ddl     << ","
	   << board   << "," 
	   << chip    << ","
	   << channel << ","
	   << strip   << ","
	   << gain    << ","  // Predictable gain
	   << board   << ","  // Predictable error
	   << strip           // Predictable chi2/ndf
	   << endl;
    }
    void Header(std::ostream& file) 
    { 
      file << "# Gains\n" 
	   << "# ddl,board,chip,channel,strip,gain,errorchi"
	   << std::endl;
    }
  };
}

//====================================================================
//
// Read back the calibrations written, and check the values 
//
void ReadBack(const char* dbBase="local://$ALICE_ROOT/OCDB/FMD/")
{
  // AliLog::SetModuleDebugLevel("FMD", 1);
  // Set specific storage of FMD ALTRO map 
  AliCDBManager::Instance()->SetDefaultStorage(Form("%s/TestCDB", dbBase));
  AliCDBManager::Instance()->SetSpecificStorage("FMD/Calib/AltroMap",
						"local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  
  AliFMDParameters* param = AliFMDParameters::Instance();
  std::cout << "Getting the stuff via AliFMDParameters ... " << std::flush;
  param->Init(kTRUE,AliFMDParameters::kPulseGain|
	      AliFMDParameters::kPedestal|
	      AliFMDParameters::kAltroMap);
  std::cout << "done" << std::endl;
  // param->Print("FMD1I[*,*]");
  
  for (UShort_t det = 3; det <= 3; det++) { 
    Char_t  rs[] = { 'I', (det==1 ? '\0' : 'O'), '\0' };
    Char_t* pr   = rs;
    Char_t  rng  = '\0';
    while ((rng = *(pr++)) != '\0') { 
      UShort_t nsec = (rng == 'I' ?  20 :  40);
      UShort_t nstr = (rng == 'I' ? 512 : 256);
      for (UShort_t sec = 0; sec < nsec; sec++) { 
	for (UShort_t str = 0; str < nstr; str++) { 
	  Float_t  gain    = param->GetPulseGain(det,rng,sec,str);
	  Float_t  ped     = param->GetPedestal(det,rng,sec,str);
	  Float_t  noise   = param->GetPedestalWidth(det,rng,sec,str);
	  UInt_t   ddl, board, chip, channel;
	  param->Detector2Hardware(det,rng,sec,str,ddl,board,chip,channel);
	  UShort_t strip   = str % 128;
	  Float_t  eped    = Hardware2Ped(ddl, board, chip, channel, strip);
	  Float_t  enoise  = Hardware2Noise(ddl, board, chip, channel, strip);
	  Float_t  egain   = Hardware2Gain(ddl, board, chip, channel, strip);
	  if (ped  != eped) 
	    Error(Form("FMD%d%c[%2d,%3d] (%d,%2d,%1d,%2d)",
		       det,rng,sec,str,ddl,board,chip,channel),
		  "pedestal=%14.7f != %14.7f", ped, eped);
	  if (noise  != enoise) 
	    Error(Form("FMD%d%c[%2d,%3d] (%d,%2d,%1d,%2d)",
		       det,rng,sec,str,ddl,board,chip,channel),
		  "noise=%14.7f != %14.7f", noise, enoise);
#if 0
	  // Will fail due to rounding errors.
	  if (gain  != egain) 
	    Error(Form("FMD%d%c[%2d,%3d] (%d,%2d,%1d,%2d)",
		       det,rng,sec,str,ddl,board,chip,channel),
		  "gain=%14.7f != %14.7f", gain, egain);
#endif
	}
      }
    }
  }
}

 
//====================================================================
//
// This script runs the test preprocessor. It uses AliTestShuttle to
// simulate a full Shuttle process
//
// The input data is created in the functions
//
//   CreateDCSAliasMap()   creates input that would in the same way come
//                         from DCS 
//   ReadDCSAliasMap()     reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be
//                         accessed by the shuttle 
//
void TestPreprocessor(const char* runType="PEDESTAL", 
		      Bool_t createDummies=kTRUE,
		      const char* dbBase="local://$ALICE_ROOT/OCDB/FMD/")
{
  // Dummy data
  if (createDummies) { 
    CreateDummyPeds pedMaker;
    pedMaker.Exec();
    CreateDummyGains gainMaker;
    gainMaker.Exec();
  }
  
  // load library - needs to be built using make
  gSystem->Load("libTestSHUTTLE.so"); 

   // create AliTestShuttle instance
  // The parameters are run, startTime, endTime
  AliTestShuttle* shuttle = new AliTestShuttle(0, 0, 1);
  AliTestShuttle::SetMainCDB("local://$ALICE_ROOT/OCDB/FMD/TestCDB");
  AliTestShuttle::SetLocalCDB(Form("%s/TestCDB", dbBase));
  AliTestShuttle::SetMainRefStorage(Form("%s/TestReference", dbBase));

  std::cout << "Test OCDB storage URI: " << AliShuttleInterface::GetMainCDB()
	    << "\n"
	    << "Test Reference storage Uri: "
	    << AliShuttleInterface::GetMainRefStorage().Data()
	    << std::endl;
  
  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "FMD", "pedestal", 
			"source1", "peds.csv");
  shuttle->AddInputFile(AliShuttleInterface::kDAQ, "FMD", "gain", 
			"source2", "gains.csv");
  shuttle->SetInputRunType(runType);
  // shuttle->SetInputRunType("PEDESTAL PULSER PHYSICS");

  new AliFMDPreprocessor(shuttle);
  // Test the preprocessor
  shuttle->Process();


  // Read back 
  ReadBack(dbBase);
}
      
//____________________________________________________________________
// We do not use this functions .....yet
TMap* CreateDCSAliasMap()
{
  // Creates a DCS structure
  // The structure is the following:
  // 
  //   TMap (key --> value)
  //     <DCSAlias> --> <valueList>
  //     <DCSAlias> is a string
  //     <valueList> is a TObjArray of AliDCSValue
  //     An AliDCSValue consists of timestamp and a value in form of a
  //     AliSimpleValue 
  // 
  // In this example 6 aliases exists: DCSAlias1 ... DCSAlias6
  // Each contains 1000 values randomly generated by TRandom::Gaus +
  // 5*nAlias 
  TRandom random;
  TMap*   aliasMap = new TMap;
  aliasMap->SetOwner(1);

  for(int nAlias=0;nAlias<6;nAlias++) {
    TObjArray* valueSet = new TObjArray;
    valueSet->SetOwner(1);

    TString aliasName="DCSAlias";
    aliasName += nAlias;
    //printf("\n\n alias: %s\n\n",aliasName.Data());
    for (int timeStamp = 0; timeStamp < 1000; timeStamp += 10) {
      Float_t      x      =  Float_t(random.Gaus()+5*nAlias);
      AliDCSValue* dcsVal = new AliDCSValue(x, timeStamp);
      valueSet->Add(dcsVal);
    }
    aliasMap->Add(new TObjString(aliasName), valueSet);
  }

  return aliasMap;
}

//____________________________________________________________________
TMap* ReadDCSAliasMap()
{
  // Open a file that contains DCS input data
  // 
  // The CDB framework is used to open the file, this means the file
  // is located in 
  //  
  //   $ALICE_ROOT/FMD/TestCDB/<detector>/DCS/Data
  // 
  // The file contains an AliCDBEntry that contains a TMap with the
  // DCS structure.  An explanation of the structure can be found in
  // CreateDCSAliasMap() 
  AliCDBEntry *entry = 
    AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
    ->Get("DET/DCS/Data", 0);
  return dynamic_cast<TMap*> (entry->GetObject());
}

//____________________________________________________________________
void WriteDCSAliasMap()
{
  // This writes the output from CreateDCSAliasMap to a CDB file

  TMap*          dcsAliasMap = CreateDCSAliasMap();
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("Responsible person");
  metaData.SetComment("Test object for TestPreprocessor.C");

  AliCDBId id("DET/DCS/Data", 0, 0);
  
  // look into AliTestShuttle's CDB main folder
  AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
    ->Put(dcsAliasMap, id, &metaData);
}

#ifndef __CINT__
int
main(int argc, char** argv)
{
  Bool_t createDummies = kTRUE;
  TString dbBase   = "local://$ALICE_ROOT/OCDB/FMD/";
  for (int i = 1; i < argc; i++) { 
    if (argv[i][0] == '-') { 
      switch (argv[i][1])  {
      case 'h': 
	std::cout << "Usage: " << argv[0] << " [OPTIONS]\n\n"
		  << "Options:\n"
		  << "\t-h\tThis help\n"
		  << "\t-d\tToggle dummies\n" 
		  << "\t-b DIR\tSet database dir\n" 
		  << std::endl;
	return 0;
      case 'd':  createDummies = !createDummies; break;
      case 'b':  dbBase = argv[++i]; break;
      }
    }
  }
  
  TestPreprocessor(createDummies, dbBase);
  return 0;
}

#endif
  
//____________________________________________________________________
//
// EOF
//
