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
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpTriggerCrate.h"
#include "AliMUONTriggerLut.h"
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
/// \author Bogdan Vulpescu, LPC Clermont-Ferrand

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
AliMUONTriggerIO::AliMUONTriggerIO(const char* regionalFileToRead) :
TObject(), fLocalBoardIds(), fNofLocalBoards(0)
{
  /// ctor
  ReadRegional(regionalFileToRead,0);
}

//_____________________________________________________________________________
AliMUONTriggerIO::~AliMUONTriggerIO()
{
  /// dtor
}

//_____________________________________________________________________________
void 
AliMUONTriggerIO::DeCompAddress(UChar_t &ypos, UChar_t &ytri, UChar_t &xdev, UChar_t &xpos, 
                                UShort_t address) const
{  
  /// decompose the 15-bits address
  
  UChar_t bitsYpos = 4;
  UChar_t bitsYtri = 1;
  UChar_t bitsXdev = 5;
  //  UChar_t bitsXpos = 5;
  
  UShort_t maskYpos = 0x000F; // ...0 00001111
  UShort_t maskYtri = 0x0001; // ...0 00000001
  UShort_t maskXdev = 0x001F; // ...0 00011111
  UShort_t maskXpos = 0x001F; // ...0 00011111
  
  ypos =  address                                  & maskYpos;
  ytri = (address >>  bitsYpos)                    & maskYtri;
  xdev = (address >> (bitsYpos+bitsYtri))          & maskXdev;
  xpos = (address >> (bitsYpos+bitsYtri+bitsXdev)) & maskXpos;
}

//_____________________________________________________________________________
void 
AliMUONTriggerIO::FillLut(AliMUONTriggerLut& lut,
                          Int_t icirc, UChar_t istripX, UChar_t idev,  
                          Int_t lutLpt[16][2], Int_t lutHpt[16][2]) 
{
  /// Fill the LUT histograms
  
  if (icirc == 0 && istripX == 0 && idev == 0) 
  {
    AliDebug(1,"Copy board, not filled ...");
    return;
  }
  
  Short_t iLptPlus, iLptMinu, iLptUnde;
  Short_t iHptPlus, iHptMinu, iHptUnde;

  iLptPlus = iLptMinu = iLptUnde = 0;
  iHptPlus = iHptMinu = iHptUnde = 0;
  
  for (Int_t istripY=0; istripY<16; istripY++) 
  {
    if (lutLpt[istripY][1] == 0 && lutLpt[istripY][0] ==1)
      iLptMinu=iLptMinu+(1 << istripY);
    if (lutLpt[istripY][1] == 1 && lutLpt[istripY][0] ==0)
      iLptPlus=iLptPlus+(1 << istripY);
    if (lutLpt[istripY][1] == 1 && lutLpt[istripY][0] ==1)
      iLptUnde=iLptUnde+(1 << istripY);
    
    if (lutHpt[istripY][1] == 0 && lutHpt[istripY][0] ==1)
      iHptMinu=iHptMinu+(1 << istripY);
    if (lutHpt[istripY][1] == 1 && lutHpt[istripY][0] ==0)
      iHptPlus=iHptPlus+(1 << istripY);
    if (lutHpt[istripY][1] == 1 && lutHpt[istripY][0] ==1)
      iHptUnde=iHptUnde+(1 << istripY);
    
  } // loop on istripY
  
  lut.SetContent("LptMinu",icirc,istripX,idev,iLptMinu);
  lut.SetContent("LptUnde",icirc,istripX,idev,iLptUnde);
  lut.SetContent("LptPlus",icirc,istripX,idev,iLptPlus);

  lut.SetContent("HptMinu",icirc,istripX,idev,iLptMinu);
  lut.SetContent("HptUnde",icirc,istripX,idev,iLptUnde);
  lut.SetContent("HptPlus",icirc,istripX,idev,iLptPlus);
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
void
AliMUONTriggerIO::ReadLocalLUT(AliMUONTriggerLut& lut,
                               Int_t localBoardId,
                               FILE* flut)
{
  /// Read the LUT for one local board from an online file

  UShort_t address;
  
  UChar_t buffer;
  UChar_t mask1 = 0xF0;
  UChar_t mask2 = 0x0F;
  UChar_t maskLpt = 0x0C;
  UChar_t maskHpt = 0x03;
  UChar_t lh, lpt, hpt;
  
  UChar_t xpos, xdev, ypos, ytri;
  
  Int_t lutLpt[16][2], lutHpt[16][2];

  Int_t boardnr = localBoardId;
  
  AliDebug(1,Form("Reading LUT values for local board %d",boardnr));
  
  Int_t ny = 0;
  
  // create the 32767 addresses for the 4-bits lpt and hpt half-bytes
  for (UShort_t ilut = 0; ilut < 0x7FFF; ilut += 2) 
  {
    // read two lut addresses at once
    fread(&buffer,1,1,flut);
    
    // 1st 4-bits half-byte
    address = ilut;   
    lh = (buffer & mask1) >> 4;
    
    // Lpt and Hpt response
    lpt = (lh & maskLpt) >> 2;
    hpt =  lh & maskHpt;
    
    // decompose the 15-bits address
    DeCompAddress(ypos,ytri,xdev,xpos,address);
    
    // calculate group of y-strips
    if (ny < 16) 
    {
      lutLpt[ny][0] =  lpt & 1;
      lutLpt[ny][1] = (lpt & 2) >> 1;
      lutHpt[ny][0] =  hpt & 1;
      lutHpt[ny][1] = (hpt & 2) >> 1;
      ny++;
      if (ny == 16) 
      {
        ny = 0;
        // ytri == 1 means no trigger in y-direction
        if (ytri == 0) 
        {
          FillLut(lut,boardnr,xpos,xdev,lutLpt,lutHpt);
        }
      }
    }
    
    // 2nd 4-bits half-byte
    address = ilut+1; 
    lh = (buffer & mask2);
    
    // Lpt and Hpt response
    lpt = (lh & maskLpt) >> 2;
    hpt =  lh & maskHpt;
    
    // decompose the 15-bits address
    DeCompAddress(ypos,ytri,xdev,xpos,address);
    
    // calculate group of y-strips
    if (ny < 16) 
    {
      lutLpt[ny][0] =  lpt & 1;
      lutLpt[ny][1] = (lpt & 2) >> 1;
      lutHpt[ny][0] =  hpt & 1;
      lutHpt[ny][1] = (hpt & 2) >> 1;
      ny++;
      if (ny == 16) 
      {
        ny = 0;
        // ytri == 1 means no trigger in y-direction
        if (ytri == 0) 
        {
          FillLut(lut,boardnr,xpos,xdev,lutLpt,lutHpt);
        }
      }
    }
  }
}

//_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::ReadLUT(const char* lutFileToRead, AliMUONTriggerLut& lut)
{
  /// Fill the LUT object from online file
  
  if ( !NofLocalBoards() )
  {
    AliError("No local board id defined. Must read a regional file first");
    return kFALSE;
  }
  
  FILE* flut = fopen(gSystem->ExpandPathName(lutFileToRead),"rb");
  if (!flut) 
  {
    AliError(Form("Could not read LUT file %s",lutFileToRead));
    return kFALSE;
  }   
  
  for ( Int_t i = 0; i < NofLocalBoards(); ++i ) 
  {
    ReadLocalLUT(lut,LocalBoardId(i),flut);
  }
  
  fclose(flut);
  
  return kTRUE;
  
}

//_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::ReadMasks(const char* localFile,
                            const char* regionalFile,
                            const char* /* globalFile */,
                            AliMUONVStore* localMasks,
                            AliMUONVStore* regionalMasks,
                            AliMUONVCalibParam* /* globalMasks */)
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
  
  if (!AliMpDDLStore::Instance(kFALSE))
  {
    AliMpCDB::LoadDDLStore();
  }
  
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

//_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::WriteLUT(const AliMUONTriggerLut& lut,
                           const char* lutFileToWrite)
{
  /// Convert an offline lut into an online (binary) lut file
  
  if ( !NofLocalBoards() )
  {
    AliError("No local board id defined. Must read a regional file first");
    return kFALSE;
  }
  
  FILE* flut = fopen(gSystem->ExpandPathName(lutFileToWrite),"wb");
  if (!flut) 
  {
    AliError(Form("Could not create output LUT file %s",lutFileToWrite));
    return kFALSE;
  }   
  
  for ( Int_t i = 0; i < NofLocalBoards(); ++i ) 
  {
    WriteLocalLUT(lut,LocalBoardId(i),flut);
  }
  
  fclose(flut);
  
  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUONTriggerIO::WriteLocalLUT(const AliMUONTriggerLut& lut,
                                Int_t localBoardId,
                                FILE* flut)
{
  /// loop over the address for the 4-bits lpt and hpt decisions

  const Int_t kMaskYpos = 0x0F;
  const Int_t kMaskYtri = 0x01;
  const Int_t kMaskXdev = 0x1F;
  const Int_t kMaskXpos = 0x1F;
  
  for (Int_t i = 0; i < 32768; ++i) 
  {
    Int_t lutLpt[2] = { 0 };
    Int_t lutHpt[2] = { 0 };
    
    // decompose address
    Int_t iYpos =   i                    & kMaskYpos;	
    Int_t iYtri = ( i >>   4           ) & kMaskYtri;
    Int_t iXdev = ( i >> ( 4 + 1 )     ) & kMaskXdev;
    Int_t iXpos = ( i >> ( 4 + 1 + 5 ) ) & kMaskXpos;
    
    // iYtri == 1 means no trigger in y-direction
    if (iYtri == 0) 
    {
      lut.GetLutOutput(localBoardId,iXpos,iXdev,iYpos,lutLpt,lutHpt);
    }
    
    UChar_t buffer;
    
    // fill byte
    if (i%2 == 0) 
    {
      // upper half-byte
      buffer = 0;	    
      buffer += lutLpt[1] << 7;
      buffer += lutLpt[0] << 6;
      buffer += lutHpt[1] << 5;
      buffer += lutHpt[0] << 4;
    } else {
      // lower half-byte
      buffer += lutLpt[1] << 3;
      buffer += lutLpt[0] << 2;
      buffer += lutHpt[1] << 1;
      buffer += lutHpt[0] << 0;
      fwrite(&buffer,1,1,flut);
    }
  }
}  
  
