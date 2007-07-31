// $Id$

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

/** @file   AliHLTTRDClusterizerComponentV2.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  A TRDClusterizer processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTRDClusterizerComponentV2.h"
#include "AliHLTTRDDefinitions.h"

#include "AliCDBManager.h"
#include "AliTRDclusterizerV2HLT.h"
#include "AliRawReaderMemory.h"

#include "TTree.h"
#include "TBranch.h"

#include <cstdlib>
#include <cerrno>
#include <string>

// this is a global object used for automatic component registration, do not use this
AliHLTTRDClusterizerComponentV2 gAliHLTTRDClusterizerComponentV2;

ClassImp(AliHLTTRDClusterizerComponentV2);
   
AliHLTTRDClusterizerComponentV2::AliHLTTRDClusterizerComponentV2()
  : AliHLTProcessor()
  , fOutputPercentage(100) // By default we copy to the output exactly what we got as input  
  , fStrorageDBpath("local://$ALICE_ROOT")
  , fClusterizer(NULL)
  , fCDB(NULL)
  , fMemReader(NULL)
{
  // Default constructor
}

AliHLTTRDClusterizerComponentV2::~AliHLTTRDClusterizerComponentV2()
{
  // Destructor
  ;
}

const char* AliHLTTRDClusterizerComponentV2::GetComponentID()
{
  // Return the component ID const char *
  return "TRDClusterizerV2"; // The ID of this component
}

void AliHLTTRDClusterizerComponentV2::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back( AliHLTTRDDefinitions::fgkDDLRawDataType );
}

AliHLTComponent_DataType AliHLTTRDClusterizerComponentV2::GetOutputDataType()
{
  // Get the output data type
  return AliHLTTRDDefinitions::fgkClusterDataType;
}

void AliHLTTRDClusterizerComponentV2::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = 0;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

AliHLTComponent* AliHLTTRDClusterizerComponentV2::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDClusterizerComponentV2;
};

int AliHLTTRDClusterizerComponentV2::DoInit( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  fOutputPercentage = 100;
  Int_t fRawDataVersion = 2;
  int i = 0;
  char* cpErr;
  while ( i < argc )
    {
      Logging( kHLTLogDebug, "HLT::TRDClusterizer::DoInit", "Arguments", "argv[%d] == %s", i, argv[i] );
      if ( !strcmp( argv[i], "output_percentage" ) )
	{
	  if ( i+1>=argc )
	    {
	      Logging(kHLTLogError, "HLT::TRDClusterizer::DoInit", "Missing Argument", "Missing output_percentage parameter");
	      return ENOTSUP;
	    }
	  Logging( kHLTLogDebug, "HLT::TRDClusterizer::DoInit", "Arguments", "argv[%d+1] == %s", i, argv[i+1] );
	  fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	  if ( *cpErr )
	    {
	      Logging(kHLTLogError, "HLT::TRDClusterizer::DoInit", "Wrong Argument", "Cannot convert output_percentage parameter '%s'", argv[i+1] );
	      return EINVAL;
	    }
	  Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoInit", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
	  i += 2;
	  continue;
	}

      if ( strcmp( argv[i], "-cdb" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      Logging(kHLTLogError, "HLT::TRDClusterizer::DoInit", "Missing Argument", "Missing -cdb argument");
	      return ENOTSUP;	      
	    }
	  fStrorageDBpath = argv[i+1];
	  Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoInit", "DB storage set", "DB storage is %s", fStrorageDBpath.c_str() );	  
	  i += 2;
	  continue;
	}      

      if ( strcmp( argv[i], "-rawver" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      Logging(kHLTLogError, "HLT::TRDClusterizer::DoInit", "Missing Argument", "Missing -rawver argument");
	      return ENOTSUP;	      
	    }
	  fRawDataVersion = atoi( argv[i+1] );
	  Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoInit", "Raw Data", "Version is %d", fRawDataVersion );	  
	  i += 2;
	  continue;
	}      

      Logging(kHLTLogError, "HLT::TRDClusterizer::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
      return EINVAL;
    }

  fCDB = AliCDBManager::Instance();
  if (!fCDB)
    {
      Logging(kHLTLogError, "HLT::TRDCalibration::DoInit", "Could not get CDB instance", "fCDB 0x%x", fCDB);
    }
  else
    {
      fCDB->SetRun(0); // THIS HAS TO BE RETRIEVED !!!
      fCDB->SetDefaultStorage(fStrorageDBpath.c_str());
      Logging(kHLTLogDebug, "HLT::TRDCalibration::DoInit", "CDB instance", "fCDB 0x%x", fCDB);
    }

  fMemReader = new AliRawReaderMemory;
  fClusterizer = new AliTRDclusterizerV2HLT("TRDCclusterizer", "TRDCclusterizer");
  fClusterizer->SetRawVersion(fRawDataVersion);
  fClusterizer->InitClusterTree();
  return 0;
}

int AliHLTTRDClusterizerComponentV2::DoDeinit()
{
  // Deinitialization of the component
  delete fMemReader;
  fMemReader = 0;
  delete fClusterizer;
  fClusterizer = 0;
  return 0;

  if (fCDB)
    {
      Logging( kHLTLogDebug, "HLT::TRDCalibration::DoDeinit", "destroy", "fCDB");
      fCDB->Destroy();
      fCDB = 0;
    }
}

int AliHLTTRDClusterizerComponentV2::DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
					    AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
					    AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks )
{
  // Process an event
  Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoEvent", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
  Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoEvent", "BLOCKS", "NofBlocks %lu", evtData.fBlockCnt );
  // Process an event
  unsigned long totalSize = 0;
  AliHLTUInt32_t fDblock_Specification = 0;

  // Loop over all input blocks in the event
  for ( unsigned long i = 0; i < evtData.fBlockCnt; i++ )
    {
      char tmp1[14], tmp2[14];
      DataType2Text( blocks[i].fDataType, tmp1 );
      DataType2Text( AliHLTTRDDefinitions::fgkDDLRawDataType, tmp2 );      
      Logging( kHLTLogDebug, "HLT::TRDClusterizer::DoEvent", "Event received", 
	       "Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	       evtData.fEventID, evtData.fEventID, tmp1, tmp2 );

      if ( blocks[i].fDataType != AliHLTTRDDefinitions::fgkDDLRawDataType ) 
	{
 	  Logging (kHLTLogError, "HLT::TRDClusterizer::DoEvent", "COMPARE FAILED", "type=%d is type=%d",
 		   blocks[i].fDataType, AliHLTTRDDefinitions::fgkDDLRawDataType);
	  continue;
	}
      fDblock_Specification = blocks[i].fSpecification;
      unsigned long blockSize = blocks[i].fSize;
      totalSize += blockSize;
    }

  void *memBufIn = calloc(totalSize, 1);
  AliHLTUInt8_t *pBuf = (AliHLTUInt8_t *)memBufIn;
  if (memBufIn == NULL)
    {
      Logging( kHLTLogError, "HLT::TRDClusterizer::DoEvent", "MEMORY", "Unable to allocate %lu bytes", totalSize);
      return -1;
    }

  // Make the memory continuous
  unsigned long copied = 0;
  for ( unsigned long i = 0; i < evtData.fBlockCnt; i++ )
    {
      if ( blocks[i].fDataType != AliHLTTRDDefinitions::fgkDDLRawDataType ) 
	continue;

      void *pos = (void*)(pBuf + copied);
      void *copyret = memcpy(pos, blocks[i].fPtr, blocks[i].fSize);
      if (copyret < 0)
	{
	  Logging( kHLTLogError, "HLT::TRDClusterizer::DoEvent", "MEMORY", "Unable to copy %lu bytes", blocks[i].fSize);
	  return -1;
	}
      copied += blocks[i].fSize;
    }

  Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoEvent", "COPY STATS", "total=%lu copied=%lu", totalSize, copied);

  fMemReader->Reset();
  fMemReader->SetMemory((UChar_t*)memBufIn, totalSize);
  //fMemReader->Reset();
  Bool_t ihead = fMemReader->ReadHeader();
  if (ihead == kTRUE)
    {
      Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoEvent", "HEADER", "Header read successfully");
    }
  else
    {
      Logging( kHLTLogError, "HLT::TRDClusterizer::DoEvent", "HEADER", "Header read ERROR");
      //return -1; -- not FATAL
    }

  fClusterizer->ResetTree();
//   Bool_t ireadD = fClusterizer->ReadDigits(fMemReader);
//   if (ireadD == kTRUE)
//     {
//       Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoEvent", "DIGITS", "Digits read successfully");
//     }
//   else
//     {
//       Logging( kHLTLogError, "HLT::TRDClusterizer::DoEvent", "DIGITS", "Digits read ERROR");
//       return -1;
//     }

//   Bool_t iclustered = fClusterizer->MakeClusters();
//   if (iclustered == kTRUE)
//     {
//       Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoEvent", "CLUSTERS", "Clustered successfully");
//     }
//   else
//     {
//       Logging( kHLTLogError, "HLT::TRDClusterizer::DoEvent", "CLUSTERS", "Clustering ERROR");
//       return -1;
//     }

  Bool_t iclustered = fClusterizer->Raw2ClustersChamber(fMemReader);
  if (iclustered == kTRUE)
    {
      Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoEvent", "CLUSTERS", "Clustered successfully");
    }
  else
    {
      Logging( kHLTLogError, "HLT::TRDClusterizer::DoEvent", "CLUSTERS", "Clustering ERROR");
      return -1;
    }

  free(memBufIn);

  //UInt_t memBufOutSize = 0;
  //   void *memBufOut = fClusterizer->WriteClustersToMemory(memBufOut, memBufOutSize);
  Int_t iNclusters = fClusterizer->GetNclusters();
  Logging( kHLTLogInfo, "HLT::TRDClusterizer::DoEvent", "COUNT", "N of Clusters = %d", iNclusters);

  // put the tree into output blocks of TObjArrays
  TTree *fcTree = fClusterizer->GetClusterTree();
  TList *lt = (TList*)fcTree->GetListOfBranches();
  TIter it(lt);
  it.Reset();
  TBranch *tb = 0;
  while ((tb = (TBranch*)it.Next()) != 0)
    {
      TObjArray *clusters = 0;
      tb->SetAddress(&clusters);
      for (Int_t icb = 0; icb < tb->GetEntries(); icb++)
	{
	  tb->GetEntry(icb);
	  PushBack(clusters, AliHLTTRDDefinitions::fgkClusterDataType, fDblock_Specification);
	}
    }

  return 0;
}
