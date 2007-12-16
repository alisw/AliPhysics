/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///
///  @file   AliHLTMUONTriggerCalibratorComponent.cxx
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   
///  @brief  Implementation of the AliHLTMUONTriggerCalibratorComponent class.
///

#include "AliHLTMUONTriggerCalibratorComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"

ClassImp(AliHLTMUONTriggerCalibratorComponent);

///////////////////////////////////////////////////////////////////////////////
// The code from here on was copied from MUONTRGda.cxx and addapted to the
// HLT framework.
//TODO: test that any of this actually works and clean up the code.
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
using namespace std;

#include <cstdio>
#include <cstdlib>

//AliRoot
#include "AliRawDataHeader.h"
#include "AliRawReaderMemory.h"
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

namespace
{

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

AliMUONVStore* gLocalMasks    = NULL;
AliMUONVStore* gRegionalMasks = NULL;
AliMUONVCalibParam* gGlobalMasks = NULL;

AliMUONTriggerIO gTriggerIO;

AliMUONVStore* gPatternStore = NULL;

TString gHistoFileName = "";

Float_t gkThreshold = 0.2;

}; // end of namespace

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
	if (gSodName.IsNull()) gSodName = "SOD";
	if (gDAName.IsNull()) gDAName = "MTG";
	if (gGlobalFileName.IsNull()) gGlobalFileName = "global.dat";
	if (gRegionalFileName.IsNull()) gRegionalFileName = "regional.dat";
	if (gLocalMaskFileName.IsNull()) gLocalMaskFileName = "localMask.dat";
	if (gLocalLutFileName.IsNull()) gLocalLutFileName = "localLut.dat";
	if (gSignatureFileName.IsNull()) gSignatureFileName = "signature.dat";

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
//      status = daqDA_FES_storeFile(file.Data(), file.Data());
      status = 0;
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
//      status = daqDA_FES_storeFile(file.Data(), file.Data());
      status = 0;
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
//      status = daqDA_FES_storeFile(file.Data(), file.Data());
      status = 0;
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
//      status = daqDA_FES_storeFile(file.Data(), file.Data());
      status = 0;
      if (status) {
	printf("Failed to export file: %s\n",gRegionalFileName.Data());
	return false;
      }
      if(gPrintLevel) printf("Export file: %s\n",gRegionalFileName.Data());
      out << gRegionalFileName.Data() << endl;
    }

    out.close();

    // export Exported file to FES anyway
//    status = daqDA_FES_storeFile(fileExp.Data(), fileExp.Data());
    status = 0;
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

//    status = daqDA_DB_getFile(gCurrentFileName.Data(), gCurrentFileName.Data());
    status = 0;
    if (status) {
      printf("Failed to get current config file from DB: %s\n",gCurrentFileName.Data());
      return false;
    }
 
    ReadFileNames();

//    status = daqDA_DB_getFile(gGlobalFileName.Data(), gGlobalFileName.Data());
    status = 0;
    if (status) {
      printf("Failed to get current config file from DB: %s\n", gGlobalFileName.Data());
      return false;
    }

//    status = daqDA_DB_getFile(gRegionalFileName.Data(), gRegionalFileName.Data());
    status = 0;
    if (status) {
      printf("Failed to get current config file from DB: %s\n",gRegionalFileName.Data());
      return false;
    }

//    status = daqDA_DB_getFile(gLocalMaskFileName.Data(), gLocalMaskFileName.Data());
    status = 0;
    if (status) {
      printf("Failed to get current config file from DB: %s\n",gLocalMaskFileName.Data());
      return false;
    }

//    status = daqDA_DB_getFile(gLocalLutFileName.Data(), gLocalLutFileName.Data());
    status = 0;
    if (status) {
      printf("Failed to get current config file from DB: %s\n",gLocalLutFileName.Data());
      return false;
    }
 
    return true;
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


AliHLTMUONTriggerCalibratorComponent::AliHLTMUONTriggerCalibratorComponent() :
	AliHLTCalibrationProcessor(),
	fSkipEvents(0),
	fMaxEvents(1000000)
{
	/// Default contructor.
}


AliHLTMUONTriggerCalibratorComponent::~AliHLTMUONTriggerCalibratorComponent()
{
	/// Default destructor.
}


const char* AliHLTMUONTriggerCalibratorComponent::GetComponentID()
{
	/// Inherited from AliHLTComponent.
	/// Returns the component ID string for this component type.
	
	return AliHLTMUONConstants::TriggerCalibratorId();
}


void AliHLTMUONTriggerCalibratorComponent::GetInputDataTypes(
		vector<AliHLTComponentDataType>& list
	)
{
	/// Inherited from AliHLTComponent.
	/// Returns the list of input block types expected by this component.

	list.clear();
	list.push_back( AliHLTMUONConstants::DDLRawDataType() );
}


AliHLTComponentDataType AliHLTMUONTriggerCalibratorComponent::GetOutputDataType()
{
	/// Inherited from AliHLTComponent.
	/// Returns the type of output block generated by this component.

	//TODO: fix.
	return AliHLTMUONConstants::DDLRawDataType();
}


void AliHLTMUONTriggerCalibratorComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	/// Inherited from AliHLTComponent.
	/// Returns an estimate of the expected output data size.

	constBase = 0;
	inputMultiplier = 2.;  //TODO: is this the correct estimate.
}


AliHLTComponent* AliHLTMUONTriggerCalibratorComponent::Spawn()
{
	/// Inherited from AliHLTComponent.
	/// Creates a new instance of AliHLTMUONTriggerCalibratorComponent.

	return new AliHLTMUONTriggerCalibratorComponent();
}


Int_t AliHLTMUONTriggerCalibratorComponent::ScanArgument(int argc, const char** argv)
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Parses the command line parameters.
	
	TString arg = argv[0];
	
	if (arg.CompareTo("-h") == 0)
	{
		HLTInfo("******************* usage **********************");
		HLTInfo("The available options are :");
		HLTInfo("-h help                   (this screen)");
		HLTInfo("");
		HLTInfo(" Output");
		HLTInfo("-r <root file>            (default = %s)", gHistoFileName.Data());
		HLTInfo("");
		HLTInfo(" Options");
		HLTInfo("-t <threshold values>     (default = %3.1f)", gkThreshold);
		HLTInfo("-d <print level>          (default = %d)", gPrintLevel);
		HLTInfo("-s <skip events>          (default = %d)", fSkipEvents);
		HLTInfo("-n <max events>           (default = %d)", fMaxEvents);
		HLTInfo("-e <execute ped/calib>    (default = %s)", gCommand.Data());
		return 0;  // Zero parameters parsed.
	}
	
	if (arg.CompareTo("-t") == 0)
	{
		if (argc < 2) return -EPROTO;
		gkThreshold = atof(argv[1]);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-e") == 0)
	{
		if (argc < 2) return -EPROTO;
		gCommand = argv[1];
		gCommand.ToLower();  // set command to lower case
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-d") == 0)
	{
		if (argc < 2) return -EPROTO;
		gPrintLevel = atoi(argv[1]);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-s") == 0)
	{
		if (argc < 2) return -EPROTO;
		fSkipEvents = atoi(argv[1]);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-n") == 0)
	{
		if (argc < 2) return -EPROTO;
		fMaxEvents = atoi(argv[1]);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-r") == 0)
	{
		if (argc < 2) return -EPROTO;
		gHistoFileName = argv[1];
		return 1;  // 1 parameter parsed.
	}
	
	// Do not know what this argument is so return an error code.
	HLTError("Bad argument %s (please check with -h)", arg.Data());
	return -EINVAL;
}


Int_t AliHLTMUONTriggerCalibratorComponent::InitCalibration()
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Initialise the calibration component.
	
	// comment out, since we do not retrieve files from database
	if (!ImportFiles())
	{
		HLTError("Import from DB failed");
		HLTError("For local test set DAQDA_TEST_DIR to the local directory where the Mtg files are located.");
		return -1;
	}
	/*
	if (!gDAFlag)
	{
		if(!ExportFiles()) return -1;
		return 0;
	}
	*/

	// read mask files
	gLocalMasks    = new AliMUON1DArray(gkNLocalBoard+9);
	gRegionalMasks = new AliMUON1DArray(16);
	gGlobalMasks   = new AliMUONCalibParamNI(1,2,1,0,0);
	gPatternStore  = new AliMUON1DArray(gkNLocalBoard+9);
	
	TString localFile    = gLocalMaskFileName;
	TString regionalFile = gRegionalFileName;
	TString globalFile   = gGlobalFileName;
	
	gTriggerIO.ReadMasks(localFile.Data(), regionalFile.Data(), globalFile.Data(),
				gLocalMasks, gRegionalMasks, gGlobalMasks, false);

	return 0;
}


Int_t AliHLTMUONTriggerCalibratorComponent::DeinitCalibration()
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Cleanup the calibration component releasing allocated memory.
	
	delete gLocalMasks;
	delete gRegionalMasks;
	delete gGlobalMasks; // in case
	delete gPatternStore;
	gLocalMasks = NULL;
	gRegionalMasks = NULL;
	gGlobalMasks = NULL;
	gPatternStore = NULL;
	
	return 0;
}


Int_t AliHLTMUONTriggerCalibratorComponent::ProcessCalibration(
		const AliHLTComponentEventData& /*evtData*/,
		AliHLTComponentTriggerData& /*trigData*/
	)
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Perform a calibration procedure on the new event.
	
	// Skip Events if needed
	if (fSkipEvents > 0)
	{
		fSkipEvents--;
		return 0;
	}
	
	// Do not process more than fMaxEvents.
	if (gNEvents >= fMaxEvents) return 0;
	if (gNEvents && gNEvents % 100 == 0)
		HLTInfo("Cumulated events %d", gNEvents);
	
	gNEvents++;
	
	gRunNumber = GetRunNo();
	
	/*
	Int_t eventType = GetRunType();
	if (eventType != 7)  // PHYSICS_EVENT - from event.h (DATE software)
	{
		HLTWarning("This event is not a physics event");
		return 0;
	}
	*/
	
	Int_t status;
	
	// containers
	AliMUONDDLTrigger*       ddlTrigger  = 0x0;
	AliMUONDarcHeader*       darcHeader  = 0x0;
	AliMUONRegHeader*        regHeader   = 0x0;
	AliMUONLocalStruct*      localStruct = 0x0;

	const AliHLTComponentBlockData* iter = NULL;
	
	// Loop over all DDL raw data input blocks and decode the event.
	iter = GetFirstInputBlock( AliHLTMUONConstants::DDLRawDataType() );
	while (iter != NULL)
	{
		// Make sure we have the correct muon trigger DDL type.
		if (not AliHLTMUONUtils::IsTriggerDDL(iter->fSpecification))
		{
			iter = GetNextInputBlock();
			continue;
		}
		
		// decoding rawdata headers
		AliRawReaderMemory* rawReader = new AliRawReaderMemory(
				reinterpret_cast<UChar_t*>(iter->fPtr), iter->fSize
			);
		rawReader->SetEquipmentID(AliHLTMUONUtils::SpecToEquipId(iter->fSpecification));
		
		// decoding MUON payload
		AliMUONRawStreamTrigger* rawStream  = new AliMUONRawStreamTrigger(rawReader);
		//rawStream->SetMaxReg(1);
		
		Int_t index = 0;
		// loops over DDL 
		while((status = rawStream->NextDDL()))
		{
			if (gPrintLevel) HLTInfo("iDDL %d", rawStream->GetDDL());
		
			ddlTrigger = rawStream->GetDDLTrigger();
			darcHeader = ddlTrigger->GetDarcHeader();
		
			if (gPrintLevel) HLTInfo("Global output %x", (Int_t)darcHeader->GetGlobalOutput());
		
			// loop over regional structures
			Int_t nReg = darcHeader->GetRegHeaderEntries();
			for(Int_t iReg = 0; iReg < nReg; ++iReg)
			{   //REG loop
				if (gPrintLevel) HLTInfo("RegionalId %d", iReg);
			
				regHeader =  darcHeader->GetRegHeaderEntry(iReg);
			
				// loop over local structures
				Int_t nLocal = regHeader->GetLocalEntries();
				for(Int_t iLocal = 0; iLocal < nLocal; ++iLocal)
				{
					localStruct = regHeader->GetLocalEntry(iLocal);
				
					Int_t localBoardId = gTriggerIO.LocalBoardId(index++);
					if (gPrintLevel) HLTInfo("local %d",  localBoardId );
				
					TArrayS xPattern(4);
					TArrayS yPattern(4);
					localStruct->GetXPattern(xPattern);
					localStruct->GetYPattern(yPattern);
					MakePattern(localBoardId, xPattern, yPattern);
				
					if (gPrintLevel)
					{
						HLTInfo("X pattern %x %x %x %x, Y pattern %x %x %x %x",
							localStruct->GetX1(), localStruct->GetX2(),localStruct->GetX3(),localStruct->GetX4(),
							localStruct->GetY1(), localStruct->GetY2(),localStruct->GetY3(),localStruct->GetY4()
						);
					}
				} // iLocal
			} // iReg
		} // NextDDL
		
		delete rawReader;
		delete rawStream;
	
		// Get next DDL raw data input block, with the same specification as defined in GetFirstInputBlock().
		iter = GetNextInputBlock();
	}
	
	if (gCommand.Contains("ped")) 
		MakePatternStore();
	
	if (gCommand.Contains("cal"))
		HLTWarning("Options %s disabled", gCommand.Data());
	//	MakePatternStore(false);
	
	if (!ExportFiles())
		return -1;
	
	HLTInfo("MUONTRKda : Run number                    : %d", gRunNumber);
	HLTInfo("MUONTRKda : Histo file generated          : %s", gHistoFileName.Data());
	HLTInfo("MUONTRKda : Nb of events used     = %d", gNEvents);
	
	//TODO:
	// PushBack data to shared memory ...
	//PushBack(.....);

	return 0;
}


Int_t AliHLTMUONTriggerCalibratorComponent::ShipDataToFXS(
		const AliHLTComponentEventData& /*evtData*/,
		AliHLTComponentTriggerData& /*trigData*/
	)
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Push the data to the FXS to ship off to the offline CDB.
	
	//TODO:
	// PushBack data to FXS ...
	//PushToFXS( ..... ) ;
	
	return 0;
}
