/*
MTR DA for online

Contact: Franck Manso <manso@clermont.in2p3.fr>
Link: http://aliceinfo.cern.ch/static/Offline/dimuon/muon_html/README_mtrda.html
Reference run: 61898, 61945 (dead channels), 61963 (noisy channels)
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


/*

Version 2 for MUONTRGda MUON trigger
Working version for: 

 reading back raw data

 Versionning of the Mtg file

 DA for ELECTRONICS_CALIBRATION_RUN (calib)
   checking dead channels

 DA for DETECTOR_CALIBRATION_RUN (ped)
   checking the noisy channels

 Interfaced with online database and file exchange server

November 2007
Ch. Finck

To be checked:
 Writing into the online database (need update of daqDAlib)
To be done:
 Looking at scalers outputs
*/

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
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUON1DArray.h"
#include "AliMUONTriggerIO.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONTriggerCrateConfig.h"

//ROOT
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TFile.h"
#include "TH1F.h"
#include "TArrayI.h"
#include "TArrayS.h"

// global variables
const Int_t gkNLocalBoard = AliMpConstants::TotalNofLocalBoards()+1;

TString gCurrentFileName("MtgCurrent.dat");
TString gLastCurrentFileName("MtgLastCurrent.dat");

TString gSodName;
Int_t   gSodFlag = 0;

TString gDAName;
Int_t   gDAFlag = 0;

TString gGlobalFileName;
TString gRegionalFileName;
TString gLocalMaskFileName;
TString gLocalLutFileName;
TString gSignatureFileName;

Int_t gGlobalFileVersion;
Int_t gRegionalFileVersion;
Int_t gLocalMaskFileVersion;
Int_t gLocalLutFileVersion;
Int_t gSignatureFileVersion;

Int_t gGlobalFileLastVersion;
Int_t gRegionalFileLastVersion;
Int_t gLocalMaskFileLastVersion;
Int_t gLocalLutFileLastVersion;

UInt_t gRunNumber = 0;
Int_t  gNEvents = 0;
Int_t  gNEventsN = 0;
Int_t  gNEventsD = 0;

Int_t gPrintLevel = 0;

AliMUONVStore* gLocalMasks    = 0x0;
AliMUONRegionalTriggerConfig* gRegionalMasks = 0x0;
AliMUONGlobalCrateConfig* gGlobalMasks = 0x0;    

AliMUONTriggerIO gTriggerIO;

Bool_t gAlgoNoisyInput = false;
Bool_t gAlgoDeadInput  = false;

Int_t   gkGlobalInputs = 4;
Int_t   gkGlobalInputLength = 32;
Float_t gkThrN = 0.1;
Float_t gkThrD = 0.1;
Int_t   gkMinEvents = 10;

Int_t gAccGlobalInputN[4][32] = {0};
Int_t gAccGlobalInputD[4][32] = {0};

//__________________________________________________________________
void WriteLastCurrentFile(TString currentFile = gLastCurrentFileName)
{

    // write last current file
    ofstream out;
    TString file;
    file = currentFile;
    out.open(file.Data());
    out << gSodName << " " << gSodFlag << endl;
    out << gDAName  << " " << gDAFlag  << endl;

    out << gGlobalFileName    << " " << gGlobalFileVersion    << endl;
    out << gRegionalFileName  << " " << gRegionalFileVersion  << endl;
    out << gLocalMaskFileName << " " << gLocalMaskFileVersion << endl;
    out << gLocalLutFileName  << " " << gLocalLutFileVersion  << endl;
    out << gSignatureFileName << " " << gSignatureFileVersion << endl;

    out.close();
}

//___________________________________________________________________________________________
Bool_t ReadCurrentFile(TString currentFile = gCurrentFileName, Bool_t lastCurrentFlag = false)
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
    sscanf(line, "%s %d", name, &gSodFlag);
    gSodName = name;
    if (gPrintLevel) printf("Sod Flag %d\n", gSodFlag);

    //read DA
    in.getline(line,80);    
    sscanf(line, "%s %d", name, &gDAFlag);
    gDAName = name;
    if (gPrintLevel) printf("DA Flag: %d\n", gDAFlag);

    // read global
    in.getline(line,80);    
    TString tmp(line);
    Int_t pos =  tmp.First(" ");
    gGlobalFileName = tmp(0, pos);
    
    if (!lastCurrentFlag) {
	gGlobalFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (gPrintLevel) printf("Global File Name: %s version: %d\n", 
			    gGlobalFileName.Data(), gGlobalFileVersion);
    } else {
	gGlobalFileLastVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (gPrintLevel) printf("Global File Name: %s last version: %d\n", 
				gGlobalFileName.Data(), gGlobalFileLastVersion);
    }

    // read regional
    in.getline(line,80);
    tmp = line;
    pos = tmp.First(" ");
    gRegionalFileName = tmp(0, pos);

    if (!lastCurrentFlag) {
	gRegionalFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (gPrintLevel) printf("Regional File Name: %s version: %d\n", 
				gRegionalFileName.Data(), gRegionalFileVersion);

    } else {
	gRegionalFileLastVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (gPrintLevel) printf("Regional File Name: %s last version: %d\n", 
				gRegionalFileName.Data(), gRegionalFileLastVersion);
    }

    // read mask
    in.getline(line,80);    
    tmp = line;
    pos = tmp.First(" ");
    gLocalMaskFileName = tmp(0, pos);

    if (!lastCurrentFlag) {
      gLocalMaskFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
      if (gPrintLevel) printf("Mask File Name: %s version: %d\n", 
			    gLocalMaskFileName.Data(), gLocalMaskFileVersion);
    } else {
      gLocalMaskFileLastVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
      if (gPrintLevel) printf("Mask File Name: %s last version: %d\n", 
			    gLocalMaskFileName.Data(), gLocalMaskFileLastVersion);
    }
    // read Lut
    in.getline(line,80);    
    tmp = line;
    pos = tmp.First(" ");
    gLocalLutFileName = tmp(0, pos);

    if (!lastCurrentFlag) {
	gLocalLutFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (gPrintLevel) printf("Lut File Name: %s version: %d\n", 
				gLocalLutFileName.Data(), gLocalLutFileVersion);
    } else {
	gLocalLutFileLastVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
	if (gPrintLevel) printf("Lut File Name: %s last version: %d\n", 
				gLocalLutFileName.Data(), gLocalLutFileLastVersion);
    }

    in.getline(line,80);    
    tmp = line;
    pos = tmp.First(" ");
    gSignatureFileName = tmp(0, pos);
    gSignatureFileVersion = atoi(tmp(pos+1, tmp.Length()-pos).Data());
    if (gPrintLevel) printf("Lut File Name: %s version: %d\n", 
			    gSignatureFileName.Data(), gSignatureFileVersion);

    return true;
}

//_____________
void ReadFileNames()
{
    // if last current file does not exist than read current file
    if (!ReadCurrentFile(gLastCurrentFileName, true)) 
    {
      ReadCurrentFile(gCurrentFileName, true);
    } 

    // any case read current file
    ReadCurrentFile();

}

//__________________
Bool_t ExportFiles()
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

    file = gGlobalFileName;
    if ((gGlobalFileLastVersion != gGlobalFileVersion) || initFES) {
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",gGlobalFileName.Data());
	return false;
      }
      if(gPrintLevel) printf("Export file: %s\n",gGlobalFileName.Data());
      out << gGlobalFileName.Data() << endl;
    }

    file = gLocalMaskFileName;  
    if ((gLocalMaskFileLastVersion != gLocalMaskFileVersion) || initFES) {
      modified = true;
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",gLocalMaskFileName.Data());
	return false;
      }
      if(gPrintLevel) printf("Export file: %s\n",gLocalMaskFileName.Data());
      out << gLocalMaskFileName.Data() << endl;
    }

    file = gLocalLutFileName;
    if ((gLocalLutFileLastVersion != gLocalLutFileVersion) || initFES) {
      modified = true;
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",gLocalLutFileName.Data());
	return false;
      }
      if(gPrintLevel) printf("Export file: %s\n",gLocalLutFileName.Data());
      out << gLocalLutFileName.Data() << endl;

    }

    // exported regional file whenever mask or/and Lut are modified
    file = gRegionalFileName;
    if ( (gRegionalFileLastVersion != gRegionalFileVersion) || modified || initFES) {
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",gRegionalFileName.Data());
	return false;
      }
      if(gPrintLevel) printf("Export file: %s\n",gRegionalFileName.Data());
      out << gRegionalFileName.Data() << endl;
    }

    out.close();

    // export Exported file to FES anyway
    status = daqDA_FES_storeFile(fileExp.Data(), fileExp.Data());
    if (status) {
      printf("Failed to export file: %s\n", fileExp.Data());
      return false;
    }
    if(gPrintLevel) printf("Export file: %s\n",fileExp.Data());

    // write last current file
    WriteLastCurrentFile();

    return true;
}

//__________________
Bool_t ImportFiles()
{
    // copy locally a file from daq detector config db 
    // The current detector is identified by detector code in variable
    // DATE_DETECTOR_CODE. It must be defined.
    // If environment variable DAQDA_TEST_DIR is defined, files are copied from DAQDA_TEST_DIR
    // instead of the database. The usual environment variables are not needed.

    Int_t status = 0;

    // offline:
    //gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/db");

    // offline: take the first two files from the local db
    //gSystem->Unsetenv("DAQDA_TEST_DIR");

    status = daqDA_DB_getFile(gCurrentFileName.Data(), gCurrentFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",gCurrentFileName.Data());
      return false;
    }
 
    ReadFileNames();

    status = daqDA_DB_getFile(gGlobalFileName.Data(), gGlobalFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n", gGlobalFileName.Data());
      return false;
    }

    // offline: take the other files from the special local storage (not db)
    //gSystem->Setenv("DAQDA_TEST_DIR", "/alisoft/Mts-files");

    status = daqDA_DB_getFile(gRegionalFileName.Data(), gRegionalFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",gRegionalFileName.Data());
      return false;
    }

    status = daqDA_DB_getFile(gLocalMaskFileName.Data(), gLocalMaskFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",gLocalMaskFileName.Data());
      return false;
    }

    status = daqDA_DB_getFile(gLocalLutFileName.Data(), gLocalLutFileName.Data());
    if (status) {
      printf("Failed to get current config file from DB: %s\n",gLocalLutFileName.Data());
      return false;
    }
 
    return true;
}

//_____________
void ReadMaskFiles()
{
    // read mask files
    gLocalMasks    = new AliMUON1DArray(gkNLocalBoard);
    gRegionalMasks = new AliMUONRegionalTriggerConfig();
    gGlobalMasks   = new AliMUONGlobalCrateConfig();

    TString localFile    = gLocalMaskFileName;
    TString regionalFile = gRegionalFileName;
    TString globalFile   = gGlobalFileName;

    gTriggerIO.ReadConfig(localFile.Data(), regionalFile.Data(), globalFile.Data(),
			 gLocalMasks, gRegionalMasks, gGlobalMasks);			
}

//______________________________________________________________
UInt_t GetFetMode()
{
  // FET mode = 3 to run algorithm for dead global inputs
  // 0x3 prepulse
  // 0x0 internal

  return gGlobalMasks->GetFetRegister(3);

}

//______________________________________________________________
void StoreGlobalInput(UInt_t *globalInput) 
{
  // accumulate and build statistics of global input values
  
  for (Int_t ii = 0; ii < gkGlobalInputs; ii++) {
    for (Int_t ib = 0; ib < gkGlobalInputLength; ib++) {
      // lsb -> msb
      if (gAlgoNoisyInput)
	gAccGlobalInputN[ii][ib] += (globalInput[ii] >> ib) & 0x1;
      if (gAlgoDeadInput)
	gAccGlobalInputD[ii][ib] += (globalInput[ii] >> ib) & 0x1;
    }
  }

}

//______________________________________________________________
void UpdateGlobalMasks() 
{
  // update the global masks
  
  // offline:
  //gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/db");
  
  Float_t rateN = 0.0, rateD = 0.0;
  UInt_t gmask[4], omask;
  Bool_t noise, deadc, withEvN, withEvD, updated = false;

  for (Int_t ii = 0; ii < gkGlobalInputs; ii++) {
    gmask[ii] = 0;

    for (Int_t ib = 0; ib < gkGlobalInputLength; ib++) {
      // lsb -> msb
      noise = false;
      deadc = false;
      withEvN = false;
      withEvD = false;
      if (gNEventsN > gkMinEvents) {
	rateN = (Float_t)gAccGlobalInputN[ii][ib]/(Float_t)gNEventsN;
	noise = (rateN > gkThrN);	
	withEvN = true;
      }
      if (gNEventsD > gkMinEvents) {
	rateD = (Float_t)gAccGlobalInputD[ii][ib]/(Float_t)gNEventsD;
	deadc = (rateD < (1.0-gkThrD));
	withEvD = true;
      }
      if (!withEvN && !withEvD) {
	// - copy the bit from the old mask
	gmask[ii] |= ((gGlobalMasks->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	printf("Mask not changed (just copy the old values)\n");
      }
      if (!withEvN && withEvD) {
	if (!deadc) {
	  // - create a new mask, set the bit to 1
	  //   not allowed!
	  //gmask[ii] |= 0x1 << ib;
	  // - copy the bit from the old mask
	  gmask[ii] |= ((gGlobalMasks->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	} else {
	  // - create a new mask, set the bit to 0
	  gmask[ii] |= 0x0 << ib;
	  printf("Found dead channel %1d:%02d \n",ii,ib);
	}
      }
      if (withEvN && !withEvD) {
	if (!noise) {
	  // - create a new mask, set the bit to 1
	  //   not allowed!
	  //gmask[ii] |= 0x1 << ib;
	  // - copy the bit from the old mask
	  gmask[ii] |= ((gGlobalMasks->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	} else {
	  // - create a new mask, set the bit to 0
	  gmask[ii] |= 0x0 << ib;
	  printf("Found noisy channel %1d:%02d \n",ii,ib);
	}
      }
      if (withEvN && withEvD) {
	if (!noise && !deadc) {
	  // - create a new mask, set the bit to 1
	  //   not allowed!
	  //gmask[ii] |= 0x1 << ib;
	  // - copy the bit from the old mask
	  gmask[ii] |= ((gGlobalMasks->GetGlobalMask(ii) >> ib) & 0x1) << ib;
	} else {
	  // - create a new mask, set the bit to 0
	  gmask[ii] |= 0x0 << ib;
	  if (noise)
	    printf("Found noisy channel %1d:%02d \n",ii,ib);
	  if (deadc)
	    printf("Found dead channel %1d:%02d \n",ii,ib);
	}
      }
    }
    printf("Global mask [%1d] %08x \n",ii,gmask[ii]);
  }

  // check if at least one mask value has been changed from previous version
  for (Int_t ii = 0; ii < gkGlobalInputs; ii++) {
    omask = gGlobalMasks->GetGlobalMask(ii);
    if (gmask[ii] != omask) {
      updated = true;
      gGlobalMasks->SetGlobalMask(ii,gmask[ii]);
    }
  }

  Int_t status = 0;
  if (updated) {
    
    // update version
    gGlobalFileVersion++;
    
    // don't change the file version ("-x.dat")
    
    gTriggerIO.WriteGlobalConfig(gGlobalFileName,gGlobalMasks);
    
    // write last current file
    WriteLastCurrentFile(gCurrentFileName);

    status = daqDA_DB_storeFile(gGlobalFileName.Data(), gGlobalFileName.Data());
    if (status) {
      printf("Failed to export file to DB: %s\n",gGlobalFileName.Data());
      return;
    }
    
    status = daqDA_DB_storeFile(gCurrentFileName.Data(), gCurrentFileName.Data());
    if (status) {
      printf("Failed to export file to DB: %s\n",gCurrentFileName.Data());
      return;
    }

  }
  
}

//*************************************************************//

// main routine
int main(Int_t argc, Char_t **argv) 
{
  
    // needed for streamer application
    gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo", "*", "TStreamerInfo",
					  "RIO", "TStreamerInfo()"); 

    Int_t skipEvents = 0;
    Int_t maxEvents  = 1000000;
    Int_t withWarnings = 0;
    Char_t inputFile[256] = "";
    inputFile[0] = 0;
    if (argc > 1)
      if (argv[1] != NULL)
        strncpy(inputFile, argv[1], 256);
      else {
        printf("MUONTRGda : No input File !\n");
        return -1;
      }
    TString flatOutputFile;

// option handler

    // decode the input line
    for (Int_t i = 1; i < argc; i++) // argument 0 is the executable name
    {
      Char_t* arg;
      
      arg = argv[i];
      if (arg[0] != '-') continue;
      switch (arg[1])
      {
      case 'f' : 
	  i++;
	  sprintf(inputFile,argv[i]);
	  break;
      case 'p' : 
	  i++;
          gkThrD = atof(argv[i]);
	  break;
      case 'c' : 
	  i++;
          gkThrN = atof(argv[i]);
	  break;
      case 'l' :
	  i++; 
	  gPrintLevel=atoi(argv[i]);
	  break;
      case 's' :
	  i++; 
	  skipEvents=atoi(argv[i]);
	  break;
      case 'n' :
	  i++; 
	  sscanf(argv[i],"%d",&maxEvents);
	  break;
      case 'w' :
	  i++; 
	  sscanf(argv[i],"%d",&withWarnings);
	  break;
      case 'h' :
	  i++;
	  printf("\n******************* %s usage **********************",argv[0]);
	  printf("\n%s -options, the available options are :",argv[0]);
	  printf("\n-h help                   (this screen)");
	  printf("\n");
	  printf("\n Input");
	  printf("\n-f <raw data file>"); 
	  printf("\n");
	  printf("\n Output");
	  printf("\n");
	  printf("\n Options");
          printf("\n-p <thr value ped (deadc)> (default = %3.1f)",gkThrD);
          printf("\n-c <thr value cal (noise)> (default = %3.1f)",gkThrN);
	  printf("\n-l <print level>           (default = %d)",gPrintLevel);
	  printf("\n-s <skip events>           (default = %d)",skipEvents);
	  printf("\n-n <max events>            (default = %d)",maxEvents);
	  printf("\n-w <decoder warnings>      (default = %d)",withWarnings);

	  printf("\n\n");
	  exit(-1);
      default :
	  printf("%s : bad argument %s (please check %s -h)\n",argv[0],argv[i],argv[0]);
	  argc = 2; exit(-1); // exit if error
      } // end of switch  
    } // end of for i  

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
    if (!ImportFiles()) {
      printf("Import from DB failed\n");
      printf("For local test set DAQDA_TEST_DIR to the local directory where the Mtg files are located \n");
      return -1;
    }

    ReadMaskFiles();

    // offline: the run number extracted from the file name
    //TString tmp(inputFile);
    //Int_t pos = tmp.First("daq");
    //tmp = tmp(pos+3,5);
    //gSystem->Setenv("DATE_RUN_NUMBER",tmp.Data());
    //gSystem->Exec("echo \"DATE_RUN_NUMBER = \" $DATE_RUN_NUMBER");

    if(!ExportFiles()) {
      printf("ExportFiles failed\n");
      return -1;
    }

    // FET is triggered by CTP
    if (GetFetMode() != 3) {
      printf("FET is not in mode 3\n");
      return -1;
    }

    // All 5 global cards are controlled by the Mts proxy
    if (gGlobalMasks->GetGlobalCrateEnable() != 0x1F) {
      printf("The MTS proxy does not control all global cards\n");
      return -1;
    }

    // The global cards are ON (active on the global inputs)
    if (!gGlobalMasks->GetMasksOn()) {
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
    while(1) 
    {
      if (gNEvents >= maxEvents) break;
      if (gNEvents && gNEvents % 100 == 0) 	
	  cout<<"Cumulated events " << gNEvents << endl;

      // check shutdown condition 
      if (daqDA_checkShutdown()) 
	  break;

      // Skip Events if needed
      while (skipEvents) {
	status = monitorGetEventDynamic(&event);
	skipEvents--;
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
      gRunNumber = rawReader->GetRunNumber();
    
      // L1Swc1
      // CALIBRATION_EVENT 
      // SYSTEM_SOFTWARE_TRIGGER_EVENT
      // DETECTOR_SOFTWARE_TRIGGER_EVENT
      gAlgoNoisyInput = false;
      gAlgoDeadInput  = false;
      if (eventType == PHYSICS_EVENT) {
	gAlgoNoisyInput = true;
	doUpdate = true;
	gNEventsN++;
      } else if (eventType == CALIBRATION_EVENT) {
	gAlgoDeadInput  = true;
	doUpdate = true;
	gNEventsD++;
	if (gRunNumber == 61963) {   // FET overwrite, not used for dead channels
	  gAlgoNoisyInput = true;
	  gNEventsN++;
	  gAlgoDeadInput  = false;
	  gNEventsD--;
	}
      } else {
	continue;
      }
      
      gNEvents++;
      if (gPrintLevel) printf("\nEvent # %d\n",gNEvents);

      // decoding MUON payload
      AliMUONRawStreamTrigger* rawStream  = new AliMUONRawStreamTrigger(rawReader);
      // ... without warnings from the decoder !!!
      if (!withWarnings)
	rawStream->DisableWarnings();

      // loops over DDL 
      while((status = rawStream->NextDDL())) {

	if (gPrintLevel) printf("iDDL %d\n", rawStream->GetDDL());

	ddlTrigger = rawStream->GetDDLTrigger();
	darcHeader = ddlTrigger->GetDarcHeader();

	if (rawStream->GetDDL() == 0) {
	  if (gPrintLevel) printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());
	  globalInput = darcHeader->GetGlobalInput();
	  StoreGlobalInput(globalInput);
	}

      } // NextDDL

      delete rawReader;
      delete rawStream;

    } // while (1)

    // update configuration files ifrequested event types were found
    if (doUpdate && gDAFlag) 
      UpdateGlobalMasks();

    timers.Stop();

    cout << "MUONTRGda : Run number                    : " << gRunNumber << endl;
    cout << "MUONTRGda : Nb of DATE events     = "         << nDateEvents << endl;
    cout << "MUONTRGda : Nb of events used     = "         << gNEvents    << endl;
    cout << "MUONTRGda : Nb of events used (noise)    = "  << gNEventsN   << endl;
    cout << "MUONTRGda : Nb of events used (deadc)    = "  << gNEventsD   << endl;
    cout << "Threshold for noisy inputs : " << gkThrN << " of the total number of used events (noise)" << endl;
    cout << "Threshold for dead inputs : " << gkThrD << " of the total number of used events (deadc)" << endl;

    printf("Execution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

    delete gLocalMasks;
    delete gRegionalMasks;
    delete gGlobalMasks; 

    return status;

}

