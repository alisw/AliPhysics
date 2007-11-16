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

To be done:
 Writing into the online database (need update of daqDAlib)
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
const Int_t gkNLocalBoard = AliMpConstants::NofLocalBoards();

TString gCommand("ped");

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

Int_t gPrintLevel = 0;

AliMUONVStore* gLocalMasks    = 0x0;
AliMUONVStore* gRegionalMasks = 0x0;
AliMUONVCalibParam* gGlobalMasks = 0x0;

AliMUONTriggerIO gTriggerIO;

AliMUONVStore* gPatternStore =  new AliMUON1DArray(gkNLocalBoard+9);

Char_t gHistoFileName[256];

Float_t gkThreshold = 0.2;

//__________________________________________________________________
void UpdateLocalMask(Int_t localBoardId, Int_t connector, Int_t strip)
{

    // update local mask
    AliMUONVCalibParam* localMask = 
	static_cast<AliMUONVCalibParam*>(gLocalMasks->FindObject(localBoardId));

    UShort_t mask = localMask->ValueAsInt(connector,0); 

    mask ^= (0x1 << strip); // set strip mask to zero

    localMask->SetValueAsInt(connector,0, mask);  
}

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
      WriteLastCurrentFile();
    } 

    // any case read current file
    ReadCurrentFile();

}

//__________________
Bool_t ExportFiles()
{

    // Export files to FES
    // Export files to DB not yet done, waiting for a version > 1.2 of daqDAlib
    // env variables have to be set (suppose by ECS ?)
    // setenv DATE_FES_PATH
    // setenv DATE_RUN_NUMBER
    // setenv DATE_ROLE_NAME
    // setenv DATE_DETECTOR_CODE

    // to be sure that env variable is set
    gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/infoLogger");

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

    if (gGlobalFileLastVersion != gGlobalFileVersion) {
      file = gGlobalFileName.Data();
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",gGlobalFileName.Data());
	return false;
      }
      out << gGlobalFileName.Data() << endl;
      if(gPrintLevel) printf("Export file: %s\n",gGlobalFileName.Data());
    }

    if (gLocalMaskFileLastVersion != gLocalMaskFileVersion) {
      modified = true;
      file = gLocalMaskFileName;
      status = daqDA_FES_storeFile(file.Data(), file.Data());
      if (status) {
	printf("Failed to export file: %s\n",gLocalMaskFileName.Data());
	return false;
      }
      if(gPrintLevel) printf("Export file: %s\n",gLocalMaskFileName.Data());
      out << gLocalMaskFileName.Data() << endl;
    }

    if (gLocalLutFileLastVersion != gLocalLutFileVersion) {
      file = gLocalLutFileName;
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
    if ( (gRegionalFileLastVersion != gRegionalFileVersion) || modified) {
      file = gRegionalFileName;
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

    // to be sure that env variable is set
    gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/db");

    Int_t status = 0;

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
    gLocalMasks    = new AliMUON1DArray(gkNLocalBoard+9);
    gRegionalMasks = new AliMUON1DArray(16);
    gGlobalMasks   = new AliMUONCalibParamNI(1,2,1,0,0);

    TString localFile    = gLocalMaskFileName;
    TString regionalFile = gRegionalFileName;
    TString globalFile   = gGlobalFileName;

    gTriggerIO.ReadMasks(localFile.Data(), regionalFile.Data(), globalFile.Data(),
			 gLocalMasks, gRegionalMasks, gGlobalMasks, false);			
}
//__________
void MakePattern(Int_t localBoardId, TArrayS& xPattern,  TArrayS& yPattern)
{

    // calculate the hit map for each strip in x and y direction
    AliMUONVCalibParam* pat = 
	static_cast<AliMUONVCalibParam*>(gPatternStore->FindObject(localBoardId));

    if (!pat) {
      pat = new AliMUONCalibParamND(2, 64, localBoardId, 0,0.); // put default wise 0.
      gPatternStore->Add(pat);	
    }

    for (Int_t i = 0; i < 4; ++i) {
      for (Int_t j = 0; j < 16; ++j) {
	
	Int_t xMask = xPattern[i];
	Int_t yMask = yPattern[i];

	Int_t index = 16*i + j;
	Double_t patOcc = 0.;

	if ( (xMask >> j ) & 0x1 ) {
	    patOcc  = pat->ValueAsDouble(index, 0) + 1.;
	    pat->SetValueAsDouble(index, 0, patOcc);
	}
	if ( (yMask >> j ) & 0x1 ) {
	    patOcc  = pat->ValueAsDouble(index, 0) + 1.;
	    pat->SetValueAsDouble(index, 0, patOcc);
	}
      }
    }

}

//__________
void MakePatternStore(Bool_t pedestal = true)
{

    // calculates the occupancy (option: store in a root file)
    // check noisy strip (pedestal true, software trigger)
    // check dead channel (pesdetal false, FET trigger)

    Int_t localBoardId = 0;
    Bool_t updated = false;

    // histo

    Char_t name[255];
    Char_t title[255];

    TH1F*  xOccHisto[243];
    TH1F*  yOccHisto[243];
    TH1F*  xPatOccHisto = 0x0;
    TH1F*  yPatOccHisto = 0x0;

    TFile*  histoFile = 0x0;

    if (gHistoFileName[0] != 0) {
      histoFile = new TFile(gHistoFileName,"RECREATE","MUON Tracking pedestals");

      sprintf(name,"pat_x");
      sprintf(title,"Occupancy for x strip");
      Int_t nx = 200;
      Float_t xmin = -0.2;
      Float_t xmax = 1.2; 
      xPatOccHisto = new TH1F(name,title,nx,xmin,xmax);
      xPatOccHisto ->SetDirectory(histoFile);

      sprintf(name,"pat_y");
      sprintf(title,"Occupancy for y strip");
      yPatOccHisto = new TH1F(name,title,nx,xmin,xmax);
      yPatOccHisto->SetDirectory(histoFile);
    
    }

    // iterator over pedestal
    TIter next(gPatternStore->CreateIterator());
    AliMUONVCalibParam* pat;
  
    while ( ( pat = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
    {
      localBoardId  = pat->ID0();

      if (gHistoFileName[0] != 0) {

	Int_t nx = 64;
	Float_t xmin = 0;
	Float_t xmax = 64;

	sprintf(name,"pat_x_%d",localBoardId);
	sprintf(title,"Occupancy for x strip, board %d",localBoardId);
	xOccHisto[localBoardId] = new TH1F(name,title,nx,xmin,xmax);

	sprintf(name,"pat_y_%d",localBoardId);
	sprintf(title,"Occupancy for y strip, board %d",localBoardId);
	yOccHisto[localBoardId] = new TH1F(name,title,nx,xmin,xmax);

      }

      for (Int_t index = 0; index < pat->Size() ; ++index) {// 64 bits for X and 64 bits for Y strips

	Double_t patXOcc  = pat->ValueAsDouble(index, 0)/(Double_t)gNEvents;
	Double_t patYOcc  = pat->ValueAsDouble(index, 1)/(Double_t)gNEvents;

	pat->SetValueAsDouble(index, 0, patXOcc);
	pat->SetValueAsDouble(index, 1, patYOcc);


	// check for x strip
	if ( (patXOcc > gkThreshold && pedestal) || (patXOcc < 1.- gkThreshold && !pedestal) ) {
	  UShort_t strip  = index % 16;
	  Int_t connector = index/16;
	  UpdateLocalMask(localBoardId, connector, strip); 
	  updated = true;
	}

	// check for y strip
	if ( (patYOcc > gkThreshold && pedestal) || (patYOcc < 1.- gkThreshold && !pedestal) ) {
	  UShort_t strip  = index % 16;
	  Int_t connector = index/16 + 4;
	  UpdateLocalMask(localBoardId, connector, strip);
	  updated = true;

	}

	if (gHistoFileName[0] != 0)  {
	  xPatOccHisto->Fill(patXOcc);
	  yPatOccHisto->Fill(patYOcc);
	  xOccHisto[localBoardId]->Fill(index, patXOcc);
	  yOccHisto[localBoardId]->Fill(index, patYOcc);

	}	
      }
    }

    if (gHistoFileName[0] != 0) {
      histoFile->Write();
      histoFile->Close();
    }


    if (updated) {

      // update version
      gLocalMaskFileVersion++;

      TString tmp(gLocalMaskFileName);
      Int_t pos = tmp.First("-");
      gLocalMaskFileName = tmp(0,pos+1) + Form("%d",gLocalMaskFileVersion) + ".dat"; 

      // write last current file
      WriteLastCurrentFile();

      gTriggerIO.WriteMasks(gLocalMaskFileName, gRegionalFileName, gGlobalFileName, gLocalMasks, gRegionalMasks, gGlobalMasks);
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
    Char_t inputFile[256];
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
      case 't' : 
	  i++;
          gkThreshold = atof(argv[i]);
	  break;
      case 'e' : 
	  i++;
	  gCommand = argv[i];
	  break;
      case 'd' :
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
     case 'r' :
	  i++; 
	  sscanf(argv[i],"%s",gHistoFileName);
	  break;
      case 'h' :
	  i++;
	  printf("\n******************* %s usage **********************",argv[0]);
	  printf("\n%s -options, the available options are :",argv[0]);
	  printf("\n-h help                   (this screen)");
	  printf("\n");
	  printf("\n Input");
	  printf("\n-f <raw data file>        (default = %s)",inputFile); 
	  printf("\n");
	  printf("\n output");
	  printf("\n-r <root file>            (default = %s)",gHistoFileName); 
	  printf("\n");
	  printf("\n Options");
          printf("\n-t <threshold values>     (default = %3.1f)",gkThreshold);
	  printf("\n-d <print level>          (default = %d)",gPrintLevel);
	  printf("\n-s <skip events>          (default = %d)",skipEvents);
	  printf("\n-n <max events>           (default = %d)",maxEvents);
	  printf("\n-e <execute ped/calib>    (default = %s)",gCommand.Data());

	  printf("\n\n");
	  exit(-1);
      default :
	  printf("%s : bad argument %s (please check %s -h)\n",argv[0],argv[i],argv[0]);
	  argc = 2; exit(-1); // exit if error
      } // end of switch  
    } // end of for i  

    // set command to lower case
    gCommand.ToLower();

    // decoding the events
  
    Int_t status;
    Int_t nDateEvents = 0;

    void* event;

    // containers
    AliMUONDDLTrigger*       ddlTrigger  = 0x0;
    AliMUONDarcHeader*       darcHeader  = 0x0;
    AliMUONRegHeader*        regHeader   = 0x0;
    AliMUONLocalStruct*      localStruct = 0x0;

    TStopwatch timers;

    timers.Start(kTRUE); 

    // comment out, since we do not retrieve files from database
    if (!ImportFiles()) {
      printf("Import from DB failed\n");
      printf("For local test set DAQDA_TEST_DIR to the local directory where the Mtg files are located \n");
      return -1;
    }

    if (!gDAFlag) {
      if(!ExportFiles()) return -1;
      return 0;
    }

    ReadMaskFiles();

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

    cout << "MUONTRKda : Reading data from file " << inputFile <<endl;

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
	cout<<"EOF found"<<endl;
	break;
      }

      nDateEvents++;

      // decoding rawdata headers
      AliRawReader *rawReader = new AliRawReaderDate(event);
 
      Int_t eventType = rawReader->GetType();
      gRunNumber = rawReader->GetRunNumber();
    

      if (eventType != PHYSICS_EVENT) 
	  continue; // for the moment

      gNEvents++;
      if (gPrintLevel) printf("\nEvent # %d\n",gNEvents);

      // decoding MUON payload
      AliMUONRawStreamTrigger* rawStream  = new AliMUONRawStreamTrigger(rawReader);
      //rawStream->SetMaxReg(1);

      Int_t index = 0;
      // loops over DDL 
      while((status = rawStream->NextDDL())) {

	if (gPrintLevel) printf("iDDL %d\n", rawStream->GetDDL());

	ddlTrigger = rawStream->GetDDLTrigger();
	darcHeader = ddlTrigger->GetDarcHeader();

	if (gPrintLevel) printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());

	// loop over regional structures
	Int_t nReg = darcHeader->GetRegHeaderEntries();
	for(Int_t iReg = 0; iReg < nReg; ++iReg){   //REG loop

	  if (gPrintLevel) printf("RegionalId %d\n", iReg);

	  regHeader =  darcHeader->GetRegHeaderEntry(iReg);

	  // loop over local structures
	  Int_t nLocal = regHeader->GetLocalEntries();
	  for(Int_t iLocal = 0; iLocal < nLocal; ++iLocal) {  

	    localStruct = regHeader->GetLocalEntry(iLocal);

	    Int_t localBoardId = gTriggerIO.LocalBoardId(index++);
	    if (gPrintLevel) printf("local %d\n",  localBoardId );

	    TArrayS xPattern(4);
	    TArrayS yPattern(4);
	    localStruct->GetXPattern(xPattern);
	    localStruct->GetYPattern(yPattern);
	    MakePattern(localBoardId, xPattern, yPattern);

	    if (gPrintLevel) printf("X pattern %x %x %x %x, Y pattern %x %x %x %x\n", 
				   localStruct->GetX1(), localStruct->GetX2(),localStruct->GetX3(),localStruct->GetX4(),
				   localStruct->GetY1(), localStruct->GetY2(),localStruct->GetY3(),localStruct->GetY4());
	   
	  } // iLocal
	} // iReg
      } // NextDDL

      delete rawReader;
      delete rawStream;

    } // while (1)

    if (gCommand.Contains("ped")) 
	MakePatternStore();

    if (gCommand.Contains("cal"))
      printf("Options %s disabled",  gCommand.Data());
    //	MakePatternStore(false);

    if (!ExportFiles())
	return -1;

    timers.Stop();

    cout << "MUONTRKda : Run number                    : " << gRunNumber << endl;
    cout << "MUONTRKda : Histo file generated          : " << gHistoFileName  << endl;
    cout << "MUONTRKda : Nb of DATE events     = "         << nDateEvents    << endl;
    cout << "MUONTRKda : Nb of events used     = "         << gNEvents        << endl;

    printf("Execution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

    delete gLocalMasks;
    delete gRegionalMasks;
    delete gGlobalMasks; // in case
    delete gPatternStore;

    return status;
}
