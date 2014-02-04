#ifndef ALIHLTTPCDATACOMPRESSIONIDMAP_H
#define ALIHLTTPCDATACOMPRESSIONIDMAP_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

// @file   AliHLTTPCDataCompressionIDMap.h
// @author Sergey Gorbunov
// @date   2013-14-05
// @brief  Map HLT cluster ID to offline index (sector,row,number) after HLT compression


#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
#include <vector>

class AliHLTComponentBlockData;

class AliHLTTPCDataCompressionIDMap: public AliHLTLogging{

 public:
  AliHLTTPCDataCompressionIDMap();
  ~AliHLTTPCDataCompressionIDMap();


  void StartFilling();
  void SetClusterBlock( const AliHLTComponentBlockData* clusters );
  void FillHLTID( AliHLTUInt32_t hltID );
  void EndFilling();

  Int_t GetNumberOfIDs(){
    return fPatchFirstID[fkNPatchesTotal];
  }

  Bool_t GetOfflineID( AliHLTUInt32_t hltID, Int_t &sector, Int_t &row, UInt_t &ind );
  Bool_t GetOfflineID( AliHLTUInt32_t hltID, AliHLTUInt32_t &offlineID);

  void SetIDMap( const AliHLTUInt8_t *data, AliHLTUInt32_t sizeBytes );
  int  WriteIDMap( AliHLTUInt8_t* output, AliHLTUInt32_t& sizeBytes );  

  // Must be public for dictionary 
  struct AliIDMapEntry{
    AliHLTUInt32_t fHltID; // HLT ID of a cluster
    AliHLTUInt32_t fOfflineID; // offline ID of a cluster
  };
 private:
  
  static Bool_t CompareIDs( const AliIDMapEntry &a, const AliIDMapEntry &b ){ return a.fHltID < b.fHltID; }
  
  AliHLTTPCDataCompressionIDMap(const AliHLTTPCDataCompressionIDMap& src);
  AliHLTTPCDataCompressionIDMap& operator=(const AliHLTTPCDataCompressionIDMap& src);

  void SetPatchIndices( UInt_t nIDs );

  static const UInt_t fkNSlices;  // !transient
  static const UInt_t fkNPatches; // !transient
  static const UInt_t fkNPatchesTotal; // NSlices*NPatches
  
  const AliHLTComponentBlockData **fPatchClusterBlocks; // pointers to data blocks (needed to read row numbers of HLT clusters)
  std::vector< AliIDMapEntry > *fInternalMemory; // memory block to fill IDs
 

  const AliIDMapEntry *fIDMap; // pointer to ID map; can point to fInternalMemory or to external memory block
  AliHLTUInt32_t *fPatchFirstID; // index of the first ID of a patch in the IDMap [size of fkNPatches+1]
  Bool_t *fPatchAllHltIDsPresent; // are all the HLT IDs of the patch present in the map [size of fkNPatches]
  
};


#endif 
