//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Sergey Gorbunov                                       *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

// @file   AliHLTTPCDataCompressionIDMap.cxx
// @author Sergey Gorbunov
// @date   2013-14-05
// @brief  Map HLT cluster ID to offline index (sector,row,number) after HLT compression

#include "AliHLTTPCDataCompressionIDMap.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTComponent.h"

#include <vector>
#include <algorithm>
#include "Riostream.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCDataCompressionIDMap)

const UInt_t AliHLTTPCDataCompressionIDMap::fkNSlices = 36;
const UInt_t AliHLTTPCDataCompressionIDMap::fkNPatches = 6;
const UInt_t AliHLTTPCDataCompressionIDMap::fkNPatchesTotal = fkNSlices*fkNPatches;

AliHLTTPCDataCompressionIDMap::AliHLTTPCDataCompressionIDMap() 
:
  fPatchClusterBlocks(NULL),
  fInternalMemory(NULL),  
  fIDMap(NULL),
  fPatchFirstID(NULL),  
  fPatchAllHltIDsPresent(NULL)
{  
  //
  // Default constructor 
  //
  fPatchClusterBlocks = new const AliHLTComponentBlockData* [fkNPatchesTotal];
  fPatchFirstID = new AliHLTUInt32_t [fkNPatchesTotal+1];
  fPatchAllHltIDsPresent = new Bool_t [fkNPatchesTotal];

  if( !fPatchClusterBlocks || !fPatchFirstID || !fPatchAllHltIDsPresent ){
    HLTError("Can not allocate memory block of %d bytes!!", fkNPatchesTotal*(sizeof(AliHLTUInt32_t) + sizeof( Bool_t) + sizeof(AliHLTComponentBlockData*)));
  }

  for( unsigned int i=0; i<fkNPatchesTotal; i++ ) fPatchClusterBlocks[i] = NULL;
  for( unsigned int i=0; i<=fkNPatchesTotal; i++ ) fPatchFirstID[i] = 0;
  for( unsigned int i=0; i<fkNPatchesTotal; i++ ) fPatchAllHltIDsPresent[i] = kFALSE;
}


AliHLTTPCDataCompressionIDMap::~AliHLTTPCDataCompressionIDMap()
{
  //
  // Default destructor.
  //  
  delete fInternalMemory;
  delete[] fPatchAllHltIDsPresent;
  delete[] fPatchFirstID;
  delete[] fPatchClusterBlocks;
}

void AliHLTTPCDataCompressionIDMap::StartFilling()
{
  //
  // Start filling of IDs
  //

  delete fInternalMemory;
  fInternalMemory = new std::vector< AliIDMapEntry >;
  if( !fInternalMemory ){
    HLTError("Can not allocate ID map, no memory left!!!");
  }
  fIDMap = NULL;
  for( unsigned int i=0; i<fkNPatchesTotal; i++ ) fPatchClusterBlocks[i] = NULL;
  for( unsigned int i=0; i<=fkNPatchesTotal; i++ ) fPatchFirstID[i] = 0;
  for( unsigned int i=0; i<fkNPatchesTotal; i++ ) fPatchAllHltIDsPresent[i] = kFALSE;
}

void AliHLTTPCDataCompressionIDMap::SetClusterBlock( const AliHLTComponentBlockData* block)
{
  //
  // Set pointer to HLT clustrers 
  //
  UInt_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*block); 
  UInt_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*block);
  UInt_t maxSlice     = AliHLTTPCDefinitions::GetMaxSliceNr(*block); 
  UInt_t maxPartition = AliHLTTPCDefinitions::GetMaxPatchNr(*block);
  if( maxSlice!=minSlice || maxPartition!=minPartition ){
    HLTError("Unexpected input block (slices %d-%d, patches %d-%d): only one TPC partition per block is expected", minSlice, maxSlice, minPartition, maxPartition );
    return ;
  }

  if( minSlice>= fkNSlices || minPartition >= fkNPatches ){
    HLTError("Wrong Slice/Patch number of input block: slice %d, patch %d", minSlice, minPartition );
    return ;
  }

  UInt_t iSlicePatch = minSlice*fkNPatches + minPartition;
  fPatchClusterBlocks[iSlicePatch] = block;
}

void AliHLTTPCDataCompressionIDMap::FillHLTID( AliHLTUInt32_t hltID)
{
  //
  // fill cluster HLT ID
  //

  if( fIDMap || !fInternalMemory ){
    HLTError("Internal error: FillHLTID() called before StartFilling()");
    return;
  }

  UInt_t slice = (UInt_t) AliHLTTPCGeometry::CluID2Slice(hltID);
  UInt_t patch = (UInt_t) AliHLTTPCGeometry::CluID2Partition(hltID);
  UInt_t number = (UInt_t) AliHLTTPCGeometry::CluID2Index(hltID);
  UInt_t slicepatch = (UInt_t) ( slice*fkNPatches + patch );

  if( slice>= fkNSlices || patch >= fkNPatches ){
    HLTError("Wrong Slice/Patch number of input cluster ID: slice %d, patch %d", slice, patch);
    return ;
  }
  
  if( !fPatchClusterBlocks[slicepatch] ){
    HLTError("Cluster block is not set for slice %d, patch %d", slice, patch);
  }

  const AliHLTTPCRawClusterData *data = reinterpret_cast<const AliHLTTPCRawClusterData*>(fPatchClusterBlocks[slicepatch]->fPtr);
  if( !data ){
    HLTError("no data found in datablock slice %d, patch %d",slice,patch);
    return;
  }

  if( number>=data->fCount ){
    HLTError("cluster number %d is outside of data block length %d",number,data->fCount);
    return;
  }
  
  Int_t sliceRow = data->fClusters[number].GetPadRow();
  sliceRow+=AliHLTTPCGeometry::GetFirstRow(patch);

  Int_t sector=0, row=0;
  Bool_t ok = AliHLTTPCGeometry::Slice2Sector( slice, sliceRow, sector, row);  
  if( !ok ) {
    HLTError("Can not transform HLT slice %d, row %d to offline numeration",slice,sliceRow);
    return;
  }
 
  AliHLTUInt32_t usec = ( (AliHLTUInt32_t) sector )&0x7F;
  AliHLTUInt32_t urow = ( (AliHLTUInt32_t) row )&0x7F;
 
  AliIDMapEntry entry;
  entry.fOfflineID = (usec<<24) + (urow<<16);
  entry.fHltID = hltID;  

  fInternalMemory->push_back(entry);
}


void AliHLTTPCDataCompressionIDMap::EndFilling()
{
  //
  // End filling of HLT IDs
  //
  
  if( fIDMap || !fInternalMemory ){
    HLTError("Internal error: EndFilling() called before StartFilling()");
    return;
  }
 
  for( unsigned int i=0; i<fkNPatchesTotal; i++ ) fPatchClusterBlocks[i] = NULL;

  const int nSec = AliHLTTPCGeometry::GetNSector();
  int nRows = AliHLTTPCGeometry::GetNRowUp1() + AliHLTTPCGeometry::GetNRowUp2();
  if( nRows < AliHLTTPCGeometry::GetNRowLow() ) nRows = AliHLTTPCGeometry::GetNRowLow();
  
  int nRowsTotal = nSec*nRows;
  AliHLTUInt32_t  *nRowClusters = new AliHLTUInt32_t  [nRowsTotal];
  if( !nRowClusters ){
    HLTError("Cannot allocate memory block of %d bytes", nRowsTotal*sizeof(AliHLTUInt32_t) );
    return;
  }
  
  for( int i=0; i<nRowsTotal; i++ ) nRowClusters[i] = 0;

  for( unsigned int i=0; i<fInternalMemory->size(); i++){
    AliIDMapEntry &entry = fInternalMemory->at(i);
    int sec = entry.fOfflineID >> 24;
    int row = (entry.fOfflineID >> 16 )&0x7F;
    int secrow = sec*nRows+row;
    if( sec>=nSec || row>=nRows ){
      delete[] nRowClusters;
      HLTError("Wrong numeration: sector %d, row %d", sec, row);
      return;
    }
    entry.fOfflineID = ( entry.fOfflineID&0xFFFF0000 ) + ( nRowClusters[secrow]&0x0000FFFF );
    nRowClusters[secrow]++;
  }
  delete[] nRowClusters;

  std::sort(fInternalMemory->begin(), fInternalMemory->end(), CompareIDs );

  fIDMap = &fInternalMemory->at(0); 

  SetPatchIndices( fInternalMemory->size() );
}

void AliHLTTPCDataCompressionIDMap::SetPatchIndices(UInt_t nIDs)
{
  //
  // Set pach indices
  //
  if( !fIDMap ) return;
  for( unsigned int i=0; i<fkNPatchesTotal; i++ ) fPatchAllHltIDsPresent[i] = kTRUE;

  UInt_t currentPatch = 0;
  UInt_t nPatchClusters=0;
  fPatchFirstID[0]=0;
  int nOK=0;
  for( unsigned int i=0; i<nIDs; i++ ){
    const AliIDMapEntry &entry = fIDMap[i];
    UInt_t slice = AliHLTTPCGeometry::CluID2Slice(entry.fHltID);
    UInt_t patch = AliHLTTPCGeometry::CluID2Partition(entry.fHltID);
    UInt_t number = AliHLTTPCGeometry::CluID2Index(entry.fHltID);
    UInt_t slicepatch = slice*fkNPatches + patch;

    if( slice>= fkNSlices || patch >= fkNPatches || slicepatch < currentPatch ){
      HLTWarning("Corrupted cluster ID block: wrong numeration: slice %d, patch %d", slice, patch );
      break;
    }        
    if( slicepatch>currentPatch ){
      for(UInt_t j = currentPatch+1; j<=slicepatch; j++) fPatchFirstID[j]=i;
      currentPatch = slicepatch;
      nPatchClusters=0;
    }
    if( number!= nPatchClusters ) {
      if( fPatchAllHltIDsPresent[slicepatch] ){
	HLTDebug("Non-continious cluster numeration for slice %d, patch %d", slice, patch);
      }
      fPatchAllHltIDsPresent[slicepatch]=kFALSE;
    }
    nPatchClusters++;
    nOK++;
  }
  
  for( UInt_t j = currentPatch+1; j<=fkNPatchesTotal; j++) fPatchFirstID[j] = nOK;

}


Bool_t AliHLTTPCDataCompressionIDMap::GetOfflineID( AliHLTUInt32_t hltID, Int_t &sector, Int_t &row, UInt_t &ind )
{
  //
  // Get cluster ID
  //

  AliHLTUInt32_t offlineID = 0;
  if( !GetOfflineID(hltID, offlineID) ){
    sector = 0;
    row = 0;
    ind = 0;
    return 0;
  }
  sector = offlineID >> 24;
  row = (offlineID >>16 )&0x7F;
  ind = offlineID &0x0000FFFF;
  return 1;
}

Bool_t AliHLTTPCDataCompressionIDMap::GetOfflineID( AliHLTUInt32_t hltID, AliHLTUInt32_t &offlineID )
{
  //
  // Get cluster ID
  //

  offlineID = 0;

  UInt_t slice = AliHLTTPCGeometry::CluID2Slice(hltID);
  UInt_t patch = AliHLTTPCGeometry::CluID2Partition(hltID);
  UInt_t number = AliHLTTPCGeometry::CluID2Index(hltID);
  UInt_t slicepatch = slice*fkNPatches + patch;
  if( !fIDMap || slice>= fkNSlices || patch >= fkNPatches ){
    HLTWarning("Wrong HLT ID: slice %d, patch %d", slice, patch);
    return 0;
  }
  const AliIDMapEntry *entry = NULL;

  if( fPatchAllHltIDsPresent[slicepatch] ){
    // look at cluster position 
    if( fPatchFirstID[slicepatch]+number<fPatchFirstID[slicepatch+1] ){
      entry = fIDMap + fPatchFirstID[slicepatch]+number;
    }
  } else { 
    // do search
    const AliIDMapEntry *pEnd = fIDMap+fPatchFirstID[fkNPatchesTotal];
    AliIDMapEntry tmp;
    tmp.fHltID = hltID;
    tmp.fOfflineID = 0;
    entry = std::lower_bound(fIDMap, pEnd, tmp, CompareIDs );
    if( entry==pEnd || entry->fHltID!=hltID ){
      entry = NULL;
    }
  }
 
  if( !entry ){
    HLTWarning("Cluster ID is not present in the map: slice %d, patch %d, cluster %d", slice, patch, number);
    return 0;
  }
  
  offlineID = entry->fOfflineID;
  return 1;
}


void AliHLTTPCDataCompressionIDMap::SetIDMap( const AliHLTUInt8_t *data, AliHLTUInt32_t sizeBytes )
{
  //
  // Set cluster ID map from input bufer
  //

  delete fInternalMemory;
  fInternalMemory = 0;
  
  fIDMap = reinterpret_cast< const AliIDMapEntry * > (data);  
  if( !fIDMap ) return;
  int nIDs = sizeBytes / sizeof( AliIDMapEntry );
  SetPatchIndices( nIDs );
}

int AliHLTTPCDataCompressionIDMap::WriteIDMap( AliHLTUInt8_t* output, AliHLTUInt32_t& sizeBytes )
{
  //
  // Write cluster ID map to output buffer
  //
  AliHLTUInt32_t maxSize = sizeBytes;
  sizeBytes = 0;
  if( !fIDMap ) return sizeBytes;
  int nIDs = fPatchFirstID[fkNPatchesTotal];
  AliHLTUInt32_t size =  nIDs * sizeof(AliIDMapEntry);
  if( size > maxSize ) return sizeBytes;
  sizeBytes = size;
  memcpy( output, (const void*)fIDMap, sizeBytes);
  return sizeBytes;
}

