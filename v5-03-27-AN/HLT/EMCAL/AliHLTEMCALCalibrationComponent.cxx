/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTEMCALCalibrationComponent.cxx
    @author Mateusz Ploskon
    @date   
    @brief  A EMCALCalibration processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include "AliHLTEMCALCalibrationComponent.h"
#include "AliHLTEMCALDefinitions.h"

#include "AliCDBManager.h"
#include "AliRawReaderMemory.h"

#include <cstdlib>
#include <cerrno>
#include <string>

// this is a global object used for automatic component registration, do not use this
AliHLTEMCALCalibrationComponent gAliHLTEMCALCalibrationComponent;

ClassImp(AliHLTEMCALCalibrationComponent);
   
AliHLTEMCALCalibrationComponent::AliHLTEMCALCalibrationComponent()
  : AliHLTCalibrationProcessor()
  , fOutputPercentage(100) // By default we copy to the output exactly what we got as input  
  , fStrorageDBpath("local://$ALICE_ROOT/OCDB")
  , fCDB(NULL)
{
  // Default constructor
}

AliHLTEMCALCalibrationComponent::~AliHLTEMCALCalibrationComponent()
{
  // Destructor
  ;
}

const char* AliHLTEMCALCalibrationComponent::GetComponentID()
{
  // Return the component ID const char *
  return "EMCALCalibration"; // The ID of this component
}

void AliHLTEMCALCalibrationComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data
  list.clear(); // We do not have any requirements for our input data type(s).
  //list.push_back( AliHLTEMCALDefinitions::fgkDDLRawDataType );
  list.push_back( AliHLTEMCALDefinitions::fgkClusterDataType );
  //list.push_back( AliHLTEMCALDefinitions::fgkEMCALESDDataType );
}

AliHLTComponent_DataType AliHLTEMCALCalibrationComponent::GetOutputDataType()
{
  // Get the output data type
  return AliHLTEMCALDefinitions::fgkCalibrationDataType;
}

void AliHLTEMCALCalibrationComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = 0;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

AliHLTComponent* AliHLTEMCALCalibrationComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTEMCALCalibrationComponent;
};

Int_t AliHLTEMCALCalibrationComponent::ScanArgument( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  fOutputPercentage = 100;
  int i = 0;
  char* cpErr;
  while ( i < argc )
    {
      Logging( kHLTLogDebug, "HLT::EMCALCalibration::ScanArgument", "Arguments", "argv[%d] == %s", i, argv[i] );
      if ( !strcmp( argv[i], "output_percentage" ) )
	{
	  if ( i+1>=argc )
	    {
	      Logging(kHLTLogError, "HLT::EMCALCalibration::ScanArgument", "Missing Argument", "Missing output_percentage parameter");
	      return ENOTSUP;
	    }
	  Logging( kHLTLogDebug, "HLT::EMCALCalibration::ScanArgument", "Arguments", "argv[%d+1] == %s", i, argv[i+1] );
	  fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	  if ( *cpErr )
	    {
	      Logging(kHLTLogError, "HLT::EMCALCalibration::ScanArgument", "Wrong Argument", "Cannot convert output_percentage parameter '%s'", argv[i+1] );
	      return EINVAL;
	    }
	  Logging( kHLTLogInfo, "HLT::EMCALCalibration::ScanArgument", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
	  i += 2;
	  continue;
	}

      if ( strcmp( argv[i], "-cdb" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      Logging(kHLTLogError, "HLT::EMCALCalibration::ScanArgument", "Missing Argument", "Missing -cdb argument");
	      return ENOTSUP;	      
	    }
	  fStrorageDBpath = argv[i+1];
	  Logging( kHLTLogInfo, "HLT::EMCALCalibration::ScanArgument", "DB storage set", "DB storage is %s", fStrorageDBpath.c_str() );	  
	  i += 2;
	  continue;
	}      

      Logging(kHLTLogError, "HLT::EMCALCalibration::ScanArgument", "Unknown Option", "Unknown option '%s'", argv[i] );
      return EINVAL;
    }
  return 0;
}

Int_t AliHLTEMCALCalibrationComponent::InitCalibration()
{
  //init the calibration
  fCDB = AliCDBManager::Instance();
  if (!fCDB)
    {
      Logging(kHLTLogError, "HLT::EMCALCalibration::InitCalibration", "Could not get CDB instance", "fCDB 0x%x", fCDB);
    }
  else
    {
      fCDB->SetRun(0); // THIS HAS TO BE RETRIEVED !!!
      fCDB->SetDefaultStorage(fStrorageDBpath.c_str());
      Logging(kHLTLogDebug, "HLT::EMCALCalibration::InitCalibration", "CDB instance", "fCDB 0x%x", fCDB);
    }
  return 0;
}

Int_t AliHLTEMCALCalibrationComponent::DeinitCalibration()
{
  // Deinitialization of the component
  if (fCDB)
    {
      Logging( kHLTLogDebug, "HLT::EMCALCalibration::DeinitCalibration", "destroy", "fCDB");
      fCDB->Destroy();
      fCDB = 0;
    }
  return 0;
}

Int_t AliHLTEMCALCalibrationComponent::ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData)
{
  // Process an event
  // Logging( kHLTLogInfo, "HLT::EMCALCalibration::ProcessCalibration", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
  Logging( kHLTLogDebug, "HLT::EMCALCalibration::ProcessCalibration", "BLOCKS", "NofBlocks %lu", evtData.fBlockCnt );
  // Process an event
  // unsigned long totalSize = 0;

  // implement a usage of the following
  //   AliHLTUInt32_t triggerDataStructSize = trigData.fStructSize;
  //   AliHLTUInt32_t triggerDataSize = trigData.fDataSize;
  //   void *triggerData = trigData.fData;
  Logging( kHLTLogDebug, "HLT::EMCALCalibration::ProcessCalibration", "Trigger data received", 
	   "Struct size %d Data size %d Data location 0x%x", 
	   trigData.fStructSize, trigData.fDataSize, (UInt_t*)trigData.fData);

  // Loop over all input blocks in the event
  int ibForce = 0;
  TObject *tobjin = (TObject *)GetFirstInputObject(AliHLTEMCALDefinitions::fgkClusterDataType , "TTree", ibForce);
  Logging( kHLTLogInfo, "HLT::EMCALCalibration::ProcessCalibration", "1stBLOCK", "Pointer = 0x%x", tobjin);
  while (tobjin)
    {
      tobjin = (TObject *)GetNextInputObject( ibForce );
      Logging( kHLTLogInfo, "HLT::EMCALCalibration::ProcessCalibration", "nextBLOCK", "Pointer = 0x%x", tobjin);
    }

  return 0;
}

Int_t AliHLTEMCALCalibrationComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  //Int_t PushToFXS(TObject* pObject, const char* pDetector, const char* pFileID, const char* pDDLNumber = "");
  // what we will actually push depends on the calibration procedure
  //ireturn = PushToFXS(object, "EMCAL ", "EMCALCalib", "1024 ");
  Logging( kHLTLogDebug, "HLT::EMCALCalibration::ProcessCalibration", "Shipping data", 
	   "Nothing serious");
  Int_t ireturn = 0;
  return ireturn;
}
