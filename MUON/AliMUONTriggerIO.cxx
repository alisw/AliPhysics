/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONTriggerIO.h"

#include "AliLog.h"
#include "AliMpDDLStore.h"
#include "AliMpTriggerCrate.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONVStore.h"
#include <Riostream.h>
#include <TSystem.h>

/// \class AliMUONTriggerIO
///
/// Handles read/write of masks and LUT to/from online files, 
/// to be used by Shuttle and Trigger DA.
/// 
/// \author Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMUONTriggerIO)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerIO::AliMUONTriggerIO() :
 TObject(), fLocalBoardIds(), fNofLocalBoards(0)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONTriggerIO::~AliMUONTriggerIO()
{
  /// dtor
}

//_____________________________________________________________________________
Int_t 
AliMUONTriggerIO::LocalBoardId(Int_t index) const
{
  /// Return the i-th localBoardId, or -1 if index is out of bounds
  if ( index >= 0 && index < fNofLocalBoards ) 
  {
    return fLocalBoardIds[index];
  }
  return -1;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerIO::ReadLocalMasks(const char* localFile, AliMUONVStore& localMasks) const
{
  /// Fills the local masks store from file
  
  if ( !NofLocalBoards() )
  {
    AliError("No local board to read");
    return 0;
  }
  
  FILE* fp = fopen(gSystem->ExpandPathName(localFile),"r");
  if (!fp)
  {
    AliError(Form("Could not read file %s",localFile));
    return 0;
  }
  
  UShort_t maskBuffer[8];
  
  Int_t nLocalBoards(0);
  
  while ( fread ( maskBuffer, 2, 8, fp ) )
  {
    Int_t localBoardId = LocalBoardId(nLocalBoards);
    AliDebug(1,Form("LB %03d X1 %4x X2 %4x X3 %4x X4 %4x "
                    "Y1 %4x Y2 %4x Y3 %4x Y4 %4x",
                    localBoardId,
                    maskBuffer[0],
                    maskBuffer[1],
                    maskBuffer[2],
                    maskBuffer[3],
                    maskBuffer[4],
                    maskBuffer[5],
                    maskBuffer[6],
                    maskBuffer[7]));
    
    if ( localBoardId ) 
    {
      AliMUONVCalibParam* localBoard = new AliMUONCalibParamNI(1,8,localBoardId,0,0);
      for ( Int_t x = 0; x < 2; ++x )
      {
        for ( Int_t y = 0; y < 4; ++y )
        {
          Int_t index = x*4+y;
          localBoard->SetValueAsInt(index,0,maskBuffer[index]);
        }
      }
      localMasks.Add(localBoard);
    }
    
    ++nLocalBoards;
  }
  
  if ( nLocalBoards != NofLocalBoards() ) 
  {
    AliError(Form("Read %d out of %d local boards",
                  nLocalBoards, NofLocalBoards()));
  }
  
  fclose(fp);
  
  return nLocalBoards;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::ReadMasks(const char* localFile,
                            const char* regionalFile,
                            const char* globalFile,
                            AliMUONVStore* localMasks,
                            AliMUONVStore* regionalMasks,
                            AliMUONVCalibParam* globalMasks)
{
  /// Fill the various masks store from files
  
  if ( !regionalFile ) 
  {
    AliError("Must have a regional file name to proceeed");
    return kFALSE;
  }
  
  Int_t nCrates = ReadRegional(regionalFile,regionalMasks);
  
  if (!nCrates) return kFALSE;
  
  if (localMasks && localFile)
  {
    Int_t nLocal = ReadLocalMasks(localFile,*localMasks);
    AliDebug(1,Form("Read masks for %d local boards",nLocal));
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerIO::ReadRegional(const char* regionalFile, AliMUONVStore* regionalMasks)
{
  /// Read regional file to fill the regional mask store *AND* 
  /// determine the order in which local boards will appear in local 
  /// and lut files.
  
  fLocalBoardIds.Reset();
  fNofLocalBoards = 0;
  
  std::ifstream in(gSystem->ExpandPathName(regionalFile));
  if (!in.good()) 
  {
    AliError(Form("Cannot read file %s",regionalFile));
    return 0;
  }

  char name[80];
  char line[80];
  
  Int_t nCrates(0);
  
  while (!in.eof())
  {
    in.getline(name,80);
    
    if (!strlen(name)) break;

    AliDebug(1,Form("Looking for crate %s",name));
    
    AliMpTriggerCrate* triggerCrate = AliMpDDLStore::Instance()->GetTriggerCrate(name);
    
    if (!triggerCrate)
    {
      AliError(Form("Mapping error : could not get crate %s",name));
      return 0;
    }
    
    ++nCrates;
    
    UShort_t id, mask;
    Int_t mode, coincidence;
    
    in.getline(line,80);    
    sscanf(line,"%hx",&id);

    in.getline(line,80);
    sscanf(line,"%d",&mode);
    
    in.getline(line,80);
    sscanf(line,"%d",&coincidence);
    
    in.getline(line,80);
    sscanf(line,"%hx",&mask);

    if ( regionalMasks ) 
    {
      AliMUONVCalibParam* regionalBoard = new AliMUONCalibParamNI(1,16,id,0,0);
      regionalBoard->SetValueAsInt(0,0,mask);
      regionalMasks->Add(regionalBoard);
      //FIXME: lines below should not be needed, as regionalMask should be only 1 16 bits word, not 16 16 bits word...
      for ( Int_t j = 1; j < 16; ++j )
      {
        regionalBoard->SetValueAsInt(j,0,0x3F);
      }      
    }
    
    AliDebug(1,Form("Name %s ID %x Mode %d Coin %d Mask %x",
                    name,id,mode,coincidence,mask));
    
    for ( Int_t i = 0; i < 16; ++i ) 
    {
      if ( (mask >> i ) & 0x1 )
      {          
        in.getline(line,80);
        char localBoardName[20];
        int j,localBoardId,switches;
        sscanf(line,"%02d %s %03d %03x",&j,localBoardName,&localBoardId,&switches);
        AliDebug(1,Form("%02d %s %03d %03x",j,localBoardName,localBoardId,switches));
        fLocalBoardIds.Set(fNofLocalBoards+1);
        fLocalBoardIds[fNofLocalBoards] = localBoardId;
        ++fNofLocalBoards;
      }
    }
  }
  
  return nCrates;  
}
