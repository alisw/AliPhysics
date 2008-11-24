/*
MTR DA for online

Contact: Franck Manso <manso@clermont.in2p3.fr>
Link: http://aliceinfo.cern.ch/static/Offline/dimuon/muon_html/README_mtrda.html
Reference run: 61898, 63698 (dead channels), 63701 (noisy channels)
Run Type:  PHYSICS (noisy channels), STANDALONE (dead channels)
DA Type: MON
Number of events needed: 1000 events for noisy and dead channels
Input Files: Rawdata file (DATE format)
Input Files from DB:
MtgGlobalCrate-<version>.dat
MtgRegionalCrate-<version>.dat
MtgLocalMask-<version>.dat
MtgLocalLut-<version>.dat
MtgCurrent.dat

Output Files: local dir (not persistent) 
ExportedFiles.dat
*/

//////////////////////////////////////////////////////////////////////////////
// Detector Algorithm for the MUON trigger configuration.                   //
//                                                                          //
// Calculates masks for the global trigger input, by looking at dead        //
// channels in calibration runs and at noisy channels in physics runs.      //
// Transfers trigger configuration files to the File Exchange Server and    //
// writes them (if modified) into the trigger configuration data base.      //
//                                                                          //
// Authors:                                                                 //
// Christian Fink (formerly at Subatech, Nantes)                            //
// Franck Manso (LPC Clermont-Ferrand, manso@clermont.in2p3.fr)             //
// Bogdan Vulpescu (LPC Clermont-Ferrand, vulpescu@clermont.in2p3.fr)       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include "AliRawReaderDate.h"

#include "AliMpConstants.h"
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONDDLTrigger.h"
#include "AliMUONVStore.h"
#include "AliMUON1DArray.h"
#include "AliMUONTriggerIO.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONTriggerCrateConfig.h"

//ROOT
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TPluginManager.h"

// structure for DA run parameters and DA working space
struct DAConfig {

  const TString fDAConfigFileName;
  const TString fCurrentFileName;
  const TString fLastCurrentFileName;

  TString fSodName;
  Int_t   fSodFlag;

  TString fDAName;
  Int_t   fDAFlag;

  TString fGlobalFileName;
  TString fRegionalFileName;
  TString fLocalMaskFileName;
  TString fLocalLutFileName;
  TString fSignatureFileName;

  Int_t   fGlobalFileVersion;
  Int_t   fRegionalFileVersion;
  Int_t   fLocalMaskFileVersion;
  Int_t   fLocalLutFileVersion;
  Int_t   fSignatureFileVersion;
  
  Int_t   fGlobalFileLastVersion;
  Int_t   fRegionalFileLastVersion;
  Int_t   fLocalMaskFileLastVersion;
  Int_t   fLocalLutFileLastVersion;

  Int_t   fNEventsN;
  Int_t   fNEventsD;

  Int_t   fPrintLevel;

  AliMUONVStore*                fLocalMasks;
  AliMUONRegionalTriggerConfig* fRegionalMasks;
  AliMUONGlobalCrateConfig*     fGlobalMasks;    

  AliMUONTriggerIO *fTriggerIO;

  Bool_t fAlgoNoisyInput;
  Bool_t fAlgoDeadInput;

  Float_t fThrN;
  Float_t fThrD;
  Int_t   fMinEvents;
  Int_t   fSkipEvents;
  Int_t   fMaxEvents;
  Bool_t  fWithWarnings;

  const Int_t fNLocalBoard;

  enum { kGlobalInputs = 4, kGlobalInputLength = 32 };

  Int_t fAccGlobalInputN[kGlobalInputs][kGlobalInputLength];
  Int_t fAccGlobalInputD[kGlobalInputs][kGlobalInputLength];

  DAConfig() : 
    fDAConfigFileName("DAConfig.txt"),
    fCurrentFileName("MtgCurrent.dat"), 
    fLastCurrentFileName("MtgLastCurrent.dat"), 
    fSodName(""),
    fSodFlag(0),
    fDAName(""),
    fDAFlag(0),
    fGlobalFileName(""),
    fRegionalFileName(""),
    fLocalMaskFileName(""),
    fLocalLutFileName(""),
    fSignatureFileName(""),
    fGlobalFileVersion(0),
    fRegionalFileVersion(0),
    fLocalMaskFileVersion(0),
    fLocalLutFileVersion(0),
    fSignatureFileVersion(0),
    fGlobalFileLastVersion(0),
    fRegionalFileLastVersion(0),
    fLocalMaskFileLastVersion(0),
    fLocalLutFileLastVersion(0),
    fNEventsN(0),
    fNEventsD(0),
    fPrintLevel(0),
    fLocalMasks(0x0),
    fRegionalMasks(0x0),
    fGlobalMasks(0x0),
    fTriggerIO(new AliMUONTriggerIO),
    fAlgoNoisyInput(false),
    fAlgoDeadInput(false),
    fThrN(0.1),
    fThrD(0.9),
    fMinEvents(10),
    fSkipEvents(0),
    fMaxEvents(65535),
    fWithWarnings(false),
    fNLocalBoard(AliMpConstants::TotalNofLocalBoards()+1)
  {
    for (Int_t ii = 0; ii < kGlobalInputs; ii++) {
      for (Int_t il = 0; il < kGlobalInputLength; il++) {
	fAccGlobalInputN[ii][il] = 0;
	fAccGlobalInputD[ii][il] = 0;
      }
    }
  }

  DAConfig (const DAConfig& cfg);
  DAConfig& operator=(const DAConfig& cfg);

};

//__________________________________________________________________
Bool_t ReadDAConfig(DAConfig& cfg)
{
  // read run parameters for the DA

    char line[80];

    TString file;
    file = cfg.fDAConfigFileName;
    std::ifstream in(gSystem->ExpandPathName(file.Data()));
    if (!in.good()) {
      printf("Cannot open DA configuration file %s ; use default values.\n",file.Data());
      return true;
    }

    TString tmp;
    Int_t pos;
    
    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.fPrintLevel = tmp.Atoi();
    
    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.fThrN = tmp.Atof();

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.fThrD = tmp.Atof();

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.fMinEvents = tmp.Atoi();

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.fSkipEvents = tmp.Atoi();

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.fMaxEvents = tmp.Atoi();

    in.getline(line,80);  
    tmp = line;
    pos = tmp.First(" ");
    tmp = tmp(0,pos);
    cfg.fWithWarnings = tmp.Atoi();
    
    return true;

}

//__________________________________________________________________
void WriteLastCurrentFile(DAConfig& cfg, TString currentFile)
{

    // write last current file
    ofstream out;
    TString file;
    file = currentFile;
    out.open(file.Data());
    out << cfg.fSodName << " " << cfg.fSodFlag << endl;
    out << cfg.fDAName  << " " << cfg.fDAFlag  << endl;

    out << cfg.fGlobalFileName    << " " << cfg.fGlobalFileVersion    << endl;
    out << cfg.fRegionalFileName  << " " << cfg.fRegionalFileVersion  << endl;
    out << cfg.fLocalMaskFileName << " " << cfg.fLocalMaskFileVersion << endl;
    out << cfg.fLocalLutFileName  << " " << cfg.fLocalLutFileVersion  << endl;
    out << cfg.fSignatureFileName << " " << cfg.fSignatureFileVersion << endl;

    out.close();
}

//___________________________________________________________________________________________
Bool_t ReadCurrentFile(DAConfig& cfg, TString currentFile, Bool_t lastCurrentFlag = false)
{

    // read last current file name and version
    char line[80];
    char name[80];

    TString file;
    file = currentFile;
    std::ifstream in(gSystem->ExpandPathName(file.Data()));
    if (!in.good()) {
      printf("Cannot open last current file %s\n",currentFile.Data());
      return false;
    }
    
    // read SOD 
    in.getline(line,80);  
    sscanf(line, "%s %d", name, &cfg.fSodFlag);
    cfg.fSodName = name;
    if (cfg.fPrintLevel) printf("Sod Flag %d\n", cfg.fSodFlag);

    //read DA
    in.getline(line,80);    
    sscanf(line, "%s %d", name, &cfg.fDAFlag);
    cfg.fDAName = name;
    if (cfg.fPrintLevel) printf("DA Flag: %d\n", cfg.fDAFlag);

    // read global
    in.getline(line,80);    
    TString tmp(line);
    Int_t pos =  tmp.First(" ");
    cfg.fGlobalFileName = tmp(0, pos);
    
    if (!lastCurrentFlag) {
	cfg.fGlobalFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (cfg.fPrintLevel) printf("Global File Name: %s version: %d\n", 
			    cfg.fGlobalFileName.Data(), cfg.fGlobalFileVersion);
    } else {
	cfg.fGlobalFileLastVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (cfg.fPrintLevel) printf("Global File Name: %s last version: %d\n", 
				cfg.fGlobalFileName.Data(), cfg.fGlobalFileLastVersion);
    }

    // read regional
    in.getline(line,80);
    tmp = line;
    pos = tmp.First(" ");
    cfg.fRegionalFileName = tmp(0, pos);

    if (!lastCurrentFlag) {
	cfg.fRegionalFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (cfg.fPrintLevel) printf("Regional File Name: %s version: %d\n", 
				cfg.fRegionalFileName.Data(), cfg.fRegionalFileVersion);

    } else {
	cfg.fRegionalFileLastVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (cfg.fPrintLevel) printf("Regional File Name: %s last version: %d\n", 
				cfg.fRegionalFileName.Data(), cfg.fRegionalFileLastVersion);
    }

    // read mask
    in.getline(line,80);    
    tmp = line;
    pos = tmp.First(" ");
    cfg.fLocalMaskFileName = tmp(0, pos);

    if (!lastCurrentFlag) {
      cfg.fLocalMaskFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
      if (cfg.fPrintLevel) printf("Mask File Name: %s version: %d\n", 
			    cfg.fLocalMaskFileName.Data(), cfg.fLocalMaskFileVersion);
    } else {
      cfg.fLocalMaskFileLastVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
      if (cfg.fPrintLevel) printf("Mask File Name: %s last version: %d\n", 
			    cfg.fLocalMaskFileName.Data(), cfg.fLocalMaskFileLastVersion);
    }
    // read Lut
    in.getline(line,80);    
    tmp = line;
    pos = tmp.First(" ");
    cfg.fLocalLutFileName = tmp(0, pos);

    if (!lastCurrentFlag) {
	cfg.fLocalLutFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (cfg.fPrintLevel) printf("Lut File Name: %s version: %d\n", 
				cfg.fLocalLutFileName.Data(), cfg.fLocalLutFileVersion);
    } else {
	cfg.fLocalLutFileLastVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (cfg.fPrintLevel) printf("Lut File Name: %s last version: %d\n", 
				cfg.fLocalLutFileName.Data(), cfg.fLocalLutFileLastVersion);
    }

    in.getline(line,80);    
    tmp = line;
    pos = tmp.First(" ");
    cfg.fSignatureFileName = tmp(0, pos);
    cfg.fSignatureFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
    if (cfg.fPrintLevel) printf("Lut File Name: %s version: %d\n", 
			    cfg.fSignatureFileName.Data(), cfg.fSignatureFileVersion);

    return true;
}

//_____________
void ReadFileNames(DAConfig& cfg)
{

    // if last current file does not exist than read current file
    if (!ReadCurrentFile(cfg,cfg.fLastCurrentFileName, true)) 
    {
      ReadCurrentFile(cfg,cfg.fCurrentFileName, true);
    } 

    // any case read current file
    ReadCurrentFile(cfg,cfg.fCurrentFileName);

}

//__________________
Bool_t ExportFiles(DAConfig& cfg)
{

    // Export files to FES
    // env variables have to be set (suppose by ECS ?)
    // setenv DATE_FES_PATH
    // setenv DATE_RUN_NUMBER
    // setenv DATE_ROLE_NAME
    // setenv DATE_DETECTOR_CODE

    // offline:
    //gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/infoLogger");

    // offline: use a dummy FES (local directory)
    //gSystem->Setenv("DAQDA_TEST_DIR", "/alisoft/FES");

    // update files
    Int_t status = 0;

    Bool_t modified = false;

    ofstream out;
    TString fileExp("ExportedFiles.dat");
    TString file;

    out.open(fileExp.Data());
    if (!out.good()) {
	printf("Failed to create file: %s\n",file.Data());
	return false;
    }      

    // check if MtgLastCurrent.dat exists
    // if not, do initial export of all files
    Bool_t initFES = false;
    if (gSystem->AccessPathName("MtgLastCurrent.dat"))
      initFES = true;
    if (initFES) printf("Copy all configuration files to the FES.\n");

    file = cfg.fGlobalFileName;
    if ((cfg.fGlobalFileLastVersion != cfg.fGlobalFileVersion) || initFES) {
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",cfg.fGlobalFileName.Data());
	return false;
      }
      if(cfg.fPrintLevel) printf("Export file: %s\n",cfg.fGlobalFileName.Data());
      out << cfg.fGlobalFileName.Data() << endl;
    }

    file = cfg.fLocalMaskFileName;  
    if ((cfg.fLocalMaskFileLastVersion != cfg.fLocalMaskFileVersion) || initFES) {
      modified = true;
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",cfg.fLocalMaskFileName.Data());
	return false;
      }
      if(cfg.fPrintLevel) printf("Export file: %s\n",cfg.fLocalMaskFileName.Data());
      out << cfg.fLocalMaskFileName.Data() << endl;
    }

    file = cfg.fLocalLutFileName;
    if ((cfg.fLocalLutFileLastVersion != cfg.fLocalLutFileVersion) || initFES) {
      modified = true;
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",cfg.fLocalLutFileName.Data());
	return false;
      }
      if(cfg.fPrintLevel) printf("Export file: %s\n",cfg.fLocalLutFileName.Data());
      out << cfg.fLocalLutFileName.Data() << endl;

    }

    // exported regional file whenever mask or/and Lut are modified
    file = cfg.fRegionalFileName;
    if ( (cfg.fRegionalFileLastVersion != cfg.fRegionalFileVersion) || modified || initFES) {
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",cfg.fRegionalFileName.Data());
	return false;
      }
      if(cfg.fPrintLevel) printf("Export file: %s\n",cfg.fRegionalFileName.Data());
      out << cfg.fRegionalFileName.Data() << endl;
    }

    out.close();

    // export Exported file to FES anyway
    status = daqDA_FES_storeFile(fileExp.Data(), fileExp.Data());
    if (status) {
      printf("Failed to export file: %s\n", fileExp.Data());
      return false;
    }
    if(cfg.fPrintLevel) printf("Export file: %s\n",fileExp.Data());

    // write last current file
    WriteLastCurrentFile(cfg,cfg.fLastCurrentFileName);

    return true;
}

//__________________
Bool_t ImportFiles(DAConfig& cfg)
{
    // copy locally a file from daq detector config db 
    // The current detector is identified by detector code in variable
    // DATE_DETECTOR_CODE. It must be defined.
    // If environment variable DAQDA_TEST_DIR is defined, files are copied from DAQDA_TEST_DIR
    // instead of the database. The usual environment variables are not needed.

    Int_t status = 0;

    // offline:
    //gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/db");

    status = daqDA_DB_getFile(cfg.fDAConfigFileName.Data(), cfg.fDAConfigFileName.Data());
    if (status) {
      printf("Failed to get DA config file from DB: %s\n",cfg.fDAConfigFileName.Data());
      return false;
    }
 
    ReadDAConfig(cfg);

    status = daqDA_DB_getFile(cfg.fCurrentFileName.Data(), cfg.fCurrentFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.fCurrentFileName.Data());
      return false;
    }
 
    ReadFileNames(cfg);

    status = daqDA_DB_getFile(cfg.fGlobalFileName.Data(), cfg.fGlobalFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n", cfg.fGlobalFileName.Data());
      return false;
    }

    status = daqDA_DB_getFile(cfg.fRegionalFileName.Data(), cfg.fRegionalFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.fRegionalFileName.Data());
      return false;
    }

    status = daqDA_DB_getFile(cfg.fLocalMaskFileName.Data(), cfg.fLocalMaskFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.fLocalMaskFileName.Data());
      return false;
    }

    status = daqDA_DB_getFile(cfg.fLocalLutFileName.Data(), cfg.fLocalLutFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",cfg.fLocalLutFileName.Data());
      return false;
    }
 
    return true;
}

//_____________
void ReadMaskFiles(DAConfig& cfg)
{

    // read mask files
    cfg.fLocalMasks    = new AliMUON1DArray(cfg.fNLocalBoard);
    cfg.fRegionalMasks = new AliMUONRegionalTriggerConfig();
    cfg.fGlobalMasks   = new AliMUONGlobalCrateConfig();

    TString localFile    = cfg.fLocalMaskFileName;
    TString regionalFile = cfg.fRegionalFileName;
    TString globalFile   = cfg.fGlobalFileName;

    cfg.fTriggerIO->ReadConfig(localFile.Data(), regionalFile.Data(), globalFile.Data(),
			 cfg.fLocalMasks, cfg.fRegionalMasks, cfg.fGlobalMasks);			
}

//______________________________________________________________
UInt_t GetFetMode(DAConfig& cfg)
{
  // FET mode = 3 to run algorithm for dead global inputs
  // 0x3 prepulse
  // 0x0 internal

  return cfg.fGlobalMasks->GetFetRegister(3);

}

//______________________________________________________________
void StoreGlobalInput(DAConfig& cfg, const UInt_t * const globalInput) 
{
  // accumulate and build statistics of global input values
  
  for (Int_t ii = 0; ii < cfg.kGlobalInputs; ii++) {
    for (Int_t ib = 0; ib < cfg.kGlobalInputLength; ib++) {
      // lsb -> msb
      if (cfg.fAlgoNoisyInput)
	cfg.fAccGlobalInputN[ii][ib] += (globalInput[ii] >> ib) & 0x1;
      if (cfg.fAlgoDeadInput)
	cfg.fAccGlobalInputD[ii][ib] += (globalInput[ii] >> ib) & 0x1;
    }
  }

}

//______________________________________________________________
void UpdateGlobalMasks(DAConfig& cfg) 
{
  // update the global masks
  
  // offline:
  //gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/db");

  Float_t rateN = 0.0, rateD = 0.0;
  UInt_t gmask[4], omask;
  Bool_t noise, deadc, withEvN, withEvD, updated = false;

  for (Int_t ii = 0; ii < cfg.kGlobalInputs; ii++) {
    gmask[ii] = 0;

    for (Int_t ib = 0; ib < cfg.kGlobalInputLength; ib++) {
      // lsb -> msb
      noise = false;
      deadc = false;
      withEvN = false;
      withEvD = false;
      if (cfg.fNEventsN > cfg.fMinEvents) {
	rateN = (Float_t)cfg.fAccGlobalInputN[ii][ib]/(Float_t)cfg.fNEventsN;
	noise = (rateN > cfg.fThrN);	
	withEvN = true;
      }
      if (cfg.fNEventsD > cfg.fMinEvents) {
	rateD = (Float_t)cfg.fAccGlobalInputD[ii][ib]/(Float_t)cfg.fNEventsD;
	deadc = (rateD < cfg.fThrD);
	withEvD = true;
      }
      if (!withEvN && !withEvD) {
	// - copy the bit from the old mask
	gmask[ii] |= ((cfg.fGlobalMasks->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	printf("Mask not changed (just copy the old values)\n");
      }
      if (!withEvN && withEvD) {
	if (!deadc) {
	  // - create a new mask, set the bit to 1
	  //   not allowed!
	  //gmask[ii] |= 0x1 << ib;
	  // - copy the bit from the old mask
	  gmask[ii] |= ((cfg.fGlobalMasks->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	} else {
	  // - create a new mask, set the bit to 0
	  gmask[ii] |= 0x0 << ib;
	  printf("Found dead  channel %1d:%02d (%4.2f) \n",ii,ib,rateD);
	}
      }
      if (withEvN && !withEvD) {
	if (!noise) {
	  // - create a new mask, set the bit to 1
	  //   not allowed!
	  //gmask[ii] |= 0x1 << ib;
	  // - copy the bit from the old mask
	  gmask[ii] |= ((cfg.fGlobalMasks->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	} else {
	  // - create a new mask, set the bit to 0
	  gmask[ii] |= 0x0 << ib;
	  printf("Found noisy channel %1d:%02d (%4.2f) \n",ii,ib,rateN);
	}
      }
      if (withEvN && withEvD) {
	if (!noise && !deadc) {
	  // - create a new mask, set the bit to 1
	  //   not allowed!
	  //gmask[ii] |= 0x1 << ib;
	  // - copy the bit from the old mask
	  gmask[ii] |= ((cfg.fGlobalMasks->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	} else {
	  // - create a new mask, set the bit to 0
	  gmask[ii] |= 0x0 << ib;
	  if (noise)
	    printf("Found noisy channel %1d:%02d (%4.2f) \n",ii,ib,rateN);
	  if (deadc)
	    printf("Found dead  channel %1d:%02d (%4.2f) \n",ii,ib,rateD);
	}
      }
    }
  }

  // check if at least one mask value has been changed from previous version
  for (Int_t ii = 0; ii < cfg.kGlobalInputs; ii++) {
    printf("Global mask [%1d] %08x \n",ii,gmask[ii]);
    omask = cfg.fGlobalMasks->GetGlobalMask(ii);
    if (gmask[ii] != omask) {
      updated = true;
      cfg.fGlobalMasks->SetGlobalMask(ii,gmask[ii]);
    }
  }

  Int_t status = 0;
  if (updated) {
    
    // update version
    cfg.fGlobalFileVersion++;
    
    // don't change the file version ("-x.dat")
    
    cfg.fTriggerIO->WriteGlobalConfig(cfg.fGlobalFileName,cfg.fGlobalMasks);
    
    // write last current file
    WriteLastCurrentFile(cfg,cfg.fCurrentFileName);

    status = daqDA_DB_storeFile(cfg.fGlobalFileName.Data(), cfg.fGlobalFileName.Data());
    if (status) {
      printf("Failed to export file to DB: %s\n",cfg.fGlobalFileName.Data());
      return;
    }
    
    status = daqDA_DB_storeFile(cfg.fCurrentFileName.Data(), cfg.fCurrentFileName.Data());
    if (status) {
      printf("Failed to export file to DB: %s\n",cfg.fCurrentFileName.Data());
      return;
    }

  }
  
}

//*************************************************************//

// main routine
int main(Int_t argc, Char_t **argv) 
{
  
    // needed for streamer application
    gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo", "*", "TStreamerInfo", "RIO", "TStreamerInfo()"); 

    /* check that we got some arguments = list of files */
    if (argc<2) {
      printf("Wrong number of arguments\n");
      return -1;
    }

    DAConfig cfg;

    Char_t inputFile[256] = "";
    inputFile[0] = 0;
    if (argc > 1)
      if (argv[1] != NULL)
        strncpy(inputFile, argv[1], 256);
      else {
        printf("MUONTRGda : No input File !\n");
        return -1;
      }

    // decoding the events
  
    Int_t status;
    Int_t nDateEvents = 0;

    void* event;

    // containers
    AliMUONDDLTrigger*       ddlTrigger  = 0x0;
    AliMUONDarcHeader*       darcHeader  = 0x0;

    TStopwatch timers;

    timers.Start(kTRUE); 

    // comment out, since we do not retrieve files from database
    if (!ImportFiles(cfg)) {
      printf("Import from DB failed\n");
      printf("For local test set DAQDA_TEST_DIR to the local directory where the Mtg files are located \n");
      return -1;
    }
    
    ReadMaskFiles(cfg);

    // offline: the run number extracted from the file name
    //TString tmp(inputFile);
    //Int_t pos = tmp.First("daq");
    //tmp = tmp(pos+3,5);
    //gSystem->Setenv("DATE_RUN_NUMBER",tmp.Data());
    //gSystem->Exec("echo \"DATE_RUN_NUMBER = \" $DATE_RUN_NUMBER");
    
    if(!ExportFiles(cfg)) {
      printf("ExportFiles failed\n");
      return -1;
    }

    // FET is triggered by CTP
    if (GetFetMode(cfg) != 3) {
      printf("FET is not in mode 3\n");
      return -1;
    }

    // All 5 global cards are controlled by the Mts proxy
    if (cfg.fGlobalMasks->GetGlobalCrateEnable() != 0x1F) {
      printf("The MTS proxy does not control all global cards\n");
      return -1;
    }

    // The global cards are ON (active on the global inputs)
    if (!cfg.fGlobalMasks->GetMasksOn()) {
      printf("Global masks are not ON\n");
      return -1;
    }
  
    // make sure to catch the "rare" calib events (1 every 50s in physics)
    const Char_t* tableSOD[]  = {"ALL", "yes", "CAL", "all", NULL, NULL};
    monitorDeclareTable(const_cast<char**>(tableSOD));

    status = monitorSetDataSource(inputFile);
    if (status) {
      cerr << "ERROR : monitorSetDataSource status (hex) = " << hex << status
	   << " " << monitorDecodeError(status) << endl;
      return -1;
    }
    status = monitorDeclareMp("MUON Trigger monitoring");
    if (status) {
      cerr << "ERROR : monitorDeclareMp status (hex) = " << hex << status
	   << " " << monitorDecodeError(status) << endl;
      return -1;
    }

    /* define wait event timeout - 1s max */
    monitorSetNowait();
    monitorSetNoWaitNetworkTimeout(1000);

    cout << "MUONTRGda : Reading data from file " << inputFile <<endl;

    UInt_t *globalInput;
    Bool_t doUpdate = false;
    Int_t runNumber = 0;
    Int_t nEvents = 0;

    while(1) 
    {
      if (nEvents >= cfg.fMaxEvents) break;
      if (nEvents && nEvents % 100 == 0) 	
	  cout<<"Cumulated events " << nEvents << endl;

      // check shutdown condition 
      if (daqDA_checkShutdown()) 
	  break;

      // Skip Events if needed
      while (cfg.fSkipEvents) {
	status = monitorGetEventDynamic(&event);
	cfg.fSkipEvents--;
      }

      // starts reading
      status = monitorGetEventDynamic(&event);
      if (status < 0)  {
	cout << "MUONTRGda : EOF found" << endl;
	break;
      }

      nDateEvents++;

      // decoding rawdata headers
      AliRawReader *rawReader = new AliRawReaderDate(event);
 
      Int_t eventType = rawReader->GetType();
      runNumber = rawReader->GetRunNumber();
    
      // L1Swc1
      // CALIBRATION_EVENT 
      // SYSTEM_SOFTWARE_TRIGGER_EVENT
      // DETECTOR_SOFTWARE_TRIGGER_EVENT
      cfg.fAlgoNoisyInput = false;
      cfg.fAlgoDeadInput  = false;
      if (eventType == PHYSICS_EVENT) {
	cfg.fAlgoNoisyInput = true;
	doUpdate = true;
	cfg.fNEventsN++;
      } else if (eventType == CALIBRATION_EVENT) {
	cfg.fAlgoDeadInput  = true;
	doUpdate = true;
	cfg.fNEventsD++;
      } else {
	continue;
      }
      
      nEvents++;
      if (cfg.fPrintLevel) printf("\nEvent # %d\n",nEvents);

      // decoding MUON payload
      AliMUONRawStreamTrigger* rawStream  = new AliMUONRawStreamTrigger(rawReader);
      // ... without warnings from the decoder !!!
      if (!cfg.fWithWarnings)
	rawStream->DisableWarnings();

      // loops over DDL 
      while((status = rawStream->NextDDL())) {

	if (cfg.fPrintLevel) printf("iDDL %d\n", rawStream->GetDDL());

	ddlTrigger = rawStream->GetDDLTrigger();
	darcHeader = ddlTrigger->GetDarcHeader();

	if (rawStream->GetDDL() == 0) {
	  if (cfg.fPrintLevel) printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());
	  globalInput = darcHeader->GetGlobalInput();
	  StoreGlobalInput(cfg,globalInput);
	}

      } // NextDDL

      delete rawReader;
      delete rawStream;

    } // while (1)

    // update configuration files ifrequested event types were found
    if (doUpdate && cfg.fDAFlag) 
      UpdateGlobalMasks(cfg);

    timers.Stop();

    cout << "MUONTRGda: DA enable: \t" << cfg.fDAFlag << endl;
    cout << "MUONTRGda: Run number: \t" << runNumber << endl;
    cout << "MUONTRGda: Nb of DATE events: \t" << nDateEvents << endl;
    cout << "MUONTRGda: Nb of events used: \t" << nEvents << endl;
    cout << "MUONTRGda: Nb of events used (noise): \t" << cfg.fNEventsN << endl;
    cout << "MUONTRGda: Nb of events used (deadc): \t" << cfg.fNEventsD << endl;
    cout << "MUONTRGda: Minumum nr of events for rate calculation: \t" << cfg.fMinEvents << endl;
    cout << "MUONTRGda: Maximum nr of analyzed events: \t" << cfg.fMaxEvents << endl;
    cout << "MUONTRGda: Skip events from start: \t" << cfg.fSkipEvents << endl;
    cout << "MUONTRGda: Threshold for noisy inputs: \t" << 100*cfg.fThrN << "%" << endl;
    cout << "MUONTRGda: Threshold for dead inputs: \t" << 100*cfg.fThrD << "%" << endl;
    cout << "MUONTRGda: Print level: \t" << cfg.fPrintLevel << endl;
    cout << "MUONTRGda: Show decoder warnings: \t" << cfg.fWithWarnings << endl;

    printf("MUONTRGda: Execution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

    delete cfg.fLocalMasks;
    delete cfg.fRegionalMasks;
    delete cfg.fGlobalMasks; 
    delete cfg.fTriggerIO;

    return status;

}

