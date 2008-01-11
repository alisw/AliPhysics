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
#include "AliMUONTriggerLut.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONVStore.h"

#include "AliMpCDB.h"
#include "AliMpHelper.h"
#include "AliMpConstants.h"
#include "AliMpDDL.h"
#include "AliMpFiles.h"
#include "AliMpDDLStore.h"
#include "AliMpLocalBoard.h"
#include "AliMpTriggerCrate.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TSystem.h>

/// \class AliMUONTriggerIO
///
/// Handles read/write of masks and LUT to/from online files, 
/// to be used by Shuttle and Trigger DA.
/// 
/// \author Laurent Aphecetche, Christian Finck Subatech
/// \author Bogdan Vulpescu, LPC Clermont-Ferrand

/// \cond CLASSIMP
ClassImp(AliMUONTriggerIO)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerIO::AliMUONTriggerIO() 
    : TObject(), 
      fRegionalTrigger(),
      fGlobalCrate()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONTriggerIO::AliMUONTriggerIO(const char* regionalFileToRead) 
    :TObject(), 
     fRegionalTrigger(),
     fGlobalCrate()
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
Bool_t 
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

  // convert deviation format
  // online: sign 1bit , dev 4bit
  // sign    dev    trigger
  // 0       1-15   mu-
  // 1       1-15   mu+
  // 0       0      mu+, mu- infinite momentum (unde)
  // 1       0      no x-trigger
  // offline: dev 5bit
  // sign    dev    trigger
  // -        0-14  mu-
  // -       16-31  mu+
  // -       15     mu+, mu- infinite momentum (unde)

  Int_t iXdevOff, iXdevOn, iXdev, sign;
  Bool_t trigx;

  iXdev = xdev;

  iXdevOn = sign = 0;
  iXdevOn +=  iXdev & 0x0F;
  sign    += (iXdev >> 4) & 0x01;
  if (iXdevOn == 0) {
    if (sign == 0) {
      iXdevOff = 15;
      trigx = kTRUE;
    } else {
      iXdevOff = 15;
      trigx = kFALSE;
    }
  } else {
    trigx = kTRUE;
    if (sign == 0) {
      iXdevOff = - iXdevOn + 15;  // gives range  0-14
    } else {
      iXdevOff = + iXdevOn + 15;  // gives range 16-30 !
    }
  }

  xdev = iXdevOff;

  return trigx;

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

  lut.SetContent("HptMinu",icirc,istripX,idev,iHptMinu);
  lut.SetContent("HptUnde",icirc,istripX,idev,iHptUnde);
  lut.SetContent("HptPlus",icirc,istripX,idev,iHptPlus);
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
      for ( Int_t index = 0; index < 8; ++index )
      {
	localBoard->SetValueAsInt(index,0,maskBuffer[index]);
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
  UChar_t maskHpt = 0x0C;
  UChar_t maskLpt = 0x03;
  UChar_t lh, lpt, hpt;
  
  UChar_t xpos, xdev, ypos, ytri;
  
  Int_t lutLpt[16][2], lutHpt[16][2];

  Int_t boardnr = localBoardId;
  
  AliDebug(1,Form("Reading LUT values for local board %d",boardnr));
  
  Int_t ny = 0;
  Bool_t trigx = kFALSE;
  
  // create the 32767 addresses for the 4-bits lpt and hpt half-bytes
  for (UShort_t ilut = 0; ilut < 0x7FFF; ilut += 2) 
  {
    // read two lut addresses at once
    fread(&buffer,1,1,flut);
    
    // 1st 4-bits half-byte
    address = ilut;   
    lh = (buffer & mask1) >> 4;
    
    // Lpt and Hpt response
    hpt = (lh & maskHpt) >> 2;
    lpt =  lh & maskLpt;
    
    // decompose the 15-bits address
    trigx = DeCompAddress(ypos,ytri,xdev,xpos,address);
    
    // calculate group of y-strips
    if (trigx && (ny < 16)) 
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
    hpt = (lh & maskHpt) >> 2;
    lpt =  lh & maskLpt;
    
    // decompose the 15-bits address
    trigx = DeCompAddress(ypos,ytri,xdev,xpos,address);
    
    // calculate group of y-strips
    if (trigx && (ny < 16)) 
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
                            const char* globalFile,
                            AliMUONVStore* localMasks,
                            AliMUONVStore* regionalMasks,
                            AliMUONVCalibParam* globalMasks,
                            Bool_t warn)
{
  /// Fill the various masks store from files
  
  if ( !regionalFile ) 
  {
    AliError("Must have a regional file name to proceeed");
    return kFALSE;
  }
  
  Int_t nCrates = ReadRegional(regionalFile,regionalMasks, warn);
  
  if (!nCrates) return kFALSE;
  
  if (localMasks && localFile)
  {
    Int_t nLocal = ReadLocalMasks(localFile,*localMasks);
    AliDebug(1,Form("Read masks for %d local boards",nLocal));
  }
  
  Int_t nDarc = ReadGlobal(globalFile, globalMasks);
  AliDebug(1,Form("Read disable for %d DARC boards",nDarc));
  
  if (!nDarc) return kFALSE;
  
  return kTRUE;
}

//_____________________________________________________________________________
 Int_t 
 AliMUONTriggerIO::ReadGlobal(const char* globalFile, AliMUONVCalibParam* globalMasks)
{
  /// read the global crate file and file corresponding mask
  /// the masks are disable bit for each crate, 8 per darc board
  /// bit value 0 means enable, 1 means disable                                                 * 
  
  Int_t nDarc = 0;
  if ( ! fGlobalCrate.ReadData(globalFile) ) return 0;
  
  UChar_t mask    = fGlobalCrate.GetFirstDarcDisable();
  ULong_t vmeAddr = fGlobalCrate.GetFirstDarcVmeAddr();   
  if (vmeAddr) nDarc++;
  globalMasks->SetValueAsInt(0,0,mask);
  
  mask    = fGlobalCrate.GetSecondDarcDisable();
  vmeAddr = fGlobalCrate.GetSecondDarcVmeAddr();    
  if (vmeAddr) nDarc++;
  globalMasks->SetValueAsInt(1,0,mask);
  
  return nDarc;
}
 
//_____________________________________________________________________________
Int_t
AliMUONTriggerIO::ReadRegional(const char* regionalFile, AliMUONVStore* regionalMasks, Bool_t warn)
{
  /// Read regional file to fill the regional mask store *AND* 
  /// determine the order in which local boards will appear in local 
  /// and lut files.
  
  if ( ! fRegionalTrigger.ReadData(regionalFile) ) return 0;

  Int_t nCrates(0);
  
  for (Int_t iCrate = 0; iCrate < fRegionalTrigger.GetNofTriggerCrates(); ++iCrate) 
  {

    AliMpTriggerCrate* crate = fRegionalTrigger.GetTriggerCrateFast(iCrate);
    
    if (warn)
    {
      AliMpTriggerCrate* triggerCrate = AliMpDDLStore::Instance()->GetTriggerCrate(crate->GetName());
          
      if (!triggerCrate)
      {
        AliError(Form("Mapping error : could not get crate %s", crate->GetName()));
	return 0;
      }
    }
    nCrates++;
    
    UShort_t masks = 0;
    if (regionalMasks != 0x0) 
    {
      masks = crate->GetMask();
      
      AliMUONVCalibParam* regionalBoard = new AliMUONCalibParamNI(1, 1, crate->GetId(), 0, 0);
      regionalBoard->SetValueAsInt(0, 0, masks);
      regionalMasks->Add(regionalBoard);
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
Bool_t
AliMUONTriggerIO::WriteMasks(const char* localFile,
			     const char* regionalFile,
			     const char* globalFile,
			     AliMUONVStore* localMasks,
			     AliMUONVStore* regionalMasks,
			     AliMUONVCalibParam* globalMasks) const
{
    /// write mask files
    Bool_t ok;
    ok  = WriteLocalMasks(localFile, *localMasks);
    ok &= WriteRegional(regionalFile, regionalMasks);
    ok &= WriteGlobal(globalFile, globalMasks);
    
    return ok;
}

 //_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::WriteGlobal(const char* globalFile, AliMUONVCalibParam* globalMasks) const
{
    /// write global file
    /// if no global masks defined take the one of configuration

  ofstream out;
  Int_t disable = 0;
  
  out.open(globalFile);
  if (!out.good())
  {
    AliError(Form("Could not create output regional file %s", globalFile));
    return kFALSE;
  }
   
  out << fGlobalCrate.GetName() << endl;

  // Jtag
  out << fGlobalCrate.GetJtagName() << endl;
  out << Form("0x%08x", fGlobalCrate.GetJtagVmeAddr()) << endl;
  out << Form("%d %d %d", fGlobalCrate.GetJtagClockDiv(), 
              fGlobalCrate.GetJtagRxPhase(), fGlobalCrate.GetJtagRdDelay()) << endl;
 
  for (Int_t i = 0; i < fGlobalCrate.GetJtagNofLines(); ++i)
    out << Form("%d ", fGlobalCrate.GetEnableJtag(i));
  out << endl;

  
  for (Int_t i = 0; i < fGlobalCrate.GetJtagNofLines(); ++i)
  {
    out << i << endl;
    for (Int_t j = 0; j < fGlobalCrate.GetJtagNofLines(); ++j)
      out << Form(" %s", fGlobalCrate.GetJtagCrateName(i,j).Data()) << endl;
  }
  
  // first darc board
  out << fGlobalCrate.GetFirstDarcName() << endl;
  out << Form("0x%08x", fGlobalCrate.GetFirstDarcVmeAddr()) << endl;
  out << fGlobalCrate.GetFirstDarcType() << endl;
  if (globalMasks != 0x0)
    disable = globalMasks->ValueAsInt(0);
  else
    disable = fGlobalCrate.GetFirstDarcDisable();
  out << Form("0x%02x", disable) << endl;
  out << Form("0x%x", fGlobalCrate.GetFirstDarcL0Delay()) << endl;
  out << Form("0x%x", fGlobalCrate.GetFirstDarcL1TimeOut()) << endl;
  
  // second darc board
  out << fGlobalCrate.GetSecondDarcName() << endl;
  out << Form("0x%08x", fGlobalCrate.GetSecondDarcVmeAddr()) << endl;
  out << fGlobalCrate.GetSecondDarcType() << endl;
  if (globalMasks != 0x0)
    disable = globalMasks->ValueAsInt(1);
  else
    disable = fGlobalCrate.GetSecondDarcDisable();
  out << Form("0x%02x", disable) << endl;
  out << Form("0x%x", fGlobalCrate.GetSecondDarcL0Delay()) << endl;
  out << Form("0x%x", fGlobalCrate.GetSecondDarcL1TimeOut()) << endl; 
  
  // global board
  out << fGlobalCrate.GetGlobalName() << endl;
  out << Form("0x%08x", fGlobalCrate.GetGlobalVmeAddr()) << endl;
  for (Int_t i = 0; i < fGlobalCrate.GetGlobalNofRegisters(); ++i)
    out << Form("0x%x", fGlobalCrate.GetGlobalRegister(i)) << endl;
  
  // Fet board
  out << fGlobalCrate.GetFetName() << endl;
  out << Form("0x%08x", fGlobalCrate.GetFetVmeAddr()) << endl;
  for (Int_t i = 0; i < fGlobalCrate.GetFetNofRegisters(); ++i)
    out << Form("0x%x", fGlobalCrate.GetFetRegister(i)) << endl;
  
  return kTRUE;
}
 

//_____________________________________________________________________________
Bool_t
AliMUONTriggerIO::WriteRegional(const char* regionalFile, AliMUONVStore* regionalMasks) const
{

    /// write regional mask with the current configuration
   /// if regional masks not defined, take the one from current configuration

    ofstream out;
    out.open(regionalFile);
          
    if (!out.good())
    {
      AliError(Form("Could not create output regional file %s",regionalFile));
      return kFALSE;
    }

    for (Int_t iCrate = 0; iCrate < fRegionalTrigger.GetNofTriggerCrates(); ++iCrate) 
    {
      AliMpTriggerCrate* crate = fRegionalTrigger.GetTriggerCrateFast(iCrate);

      out << crate->GetName()  << endl;
      out << Form("%02x", crate->GetId())   << endl;
      out << crate->GetMode()  << endl;
      out << crate->GetCoinc() << endl;
      
      UShort_t masks = 0;
      if (regionalMasks != 0x0) 
      {
        AliMUONVCalibParam* maskParam = 
            static_cast<AliMUONVCalibParam*>(regionalMasks->FindObject(crate->GetId()));
        masks = maskParam->ValueAsInt(0,0);
      } 
      else
      {
        masks = crate->GetMask();
      } 
      
      out << Form("%04x", masks) << endl;
      
      for (Int_t iLocal = 0; iLocal < crate->GetNofLocalBoards(); ++iLocal) 
      {
	Int_t localBoardId = crate->GetLocalBoardId(iLocal);

	AliMpLocalBoard* board = fRegionalTrigger.FindLocalBoard(localBoardId);

	out << Form("%02d ", board->GetSlot())  
	    << board->GetName() 
	    << Form(" %03d ", localBoardId) 
	    << Form("%03x", board->GetSwitch()) 
	    << endl;
 
        out << " ";
        for (Int_t i = 0; i < board->GetNofDEs(); ++i)
          out << Form("%4d ", board->GetDEId(i));
        out << endl;
          
          // print copy card numbers
        out << Form(" %4d %4d", board->GetInputXfrom(), board->GetInputXto());
        out << Form(" %4d %4d", board->GetInputYfrom(), board->GetInputYto());
        out << Form(" %4d",     board->GetTC()) << endl;
      }
    }
    out.close();
    
    return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::WriteLocalMasks(const char* localFile, AliMUONVStore& localMasks) const
{
    /// write local mask
    /// removing/adding enable for a local board need a update of the configuration 
    /// before calling this method
    /// mask are written for all boards including the copy card (Ch.F.)

    FILE* fp = fopen(gSystem->ExpandPathName(localFile),"wb");
    if (!fp) 
    {
      AliError(Form("Could not create output local file %s",localFile));
      return kFALSE;
    }   

    UShort_t maskBuffer[8];

    for (Int_t iCrate = 0; iCrate < fRegionalTrigger.GetNofTriggerCrates(); ++iCrate) 
    {
      AliMpTriggerCrate* crate = fRegionalTrigger.GetTriggerCrateFast(iCrate);
      
      UShort_t mask = crate->GetMask(); // getting mask from current config

      for (Int_t iLocal = 0; iLocal < crate->GetNofLocalBoards(); ++iLocal) 
      {
	Int_t localBoardId = crate->GetLocalBoardId(iLocal);

	if ( (mask >> iLocal ) & 0x1 ) 
	{
	  AliMUONVCalibParam* localMask = 
	      static_cast<AliMUONVCalibParam*>(localMasks.FindObject(localBoardId));

	  for (Int_t index = 0; index < 8; ++index) 
	  {
	    maskBuffer[index] = localMask->ValueAsInt(index,0); 
	  }

	  fwrite ( maskBuffer, 2, 8, fp); 
	}

      }
    }

    fclose(fp);

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
    
    // convert deviation format
    // online: sign 1bit , dev 4bit
    // sign    dev    trigger
    // 0       1-15   mu-
    // 1       1-15   mu+
    // 0       0      mu+, mu- infinite momentum (unde)
    // 1       0      no x-trigger
    // offline: dev 5bit
    // sign    dev    trigger
    // -        0-14  mu-
    // -       16-31  mu+
    // -       15     mu+, mu- infinite momentum (unde)
    Int_t iXdevOn  = 0;
    Int_t iXdevOff = 0;
    Int_t sign     = 0;
    Bool_t trigx = kFALSE;
    iXdevOn +=  iXdev & 0x0F;
    sign    += (iXdev >> 4) & 0x01;
    if (iXdevOn == 0) {
      if (sign == 0) {
	iXdevOff = 15;
	trigx = kTRUE;
      } else {
	iXdevOff = 15;
	trigx = kFALSE;
      }
    } else {
      trigx = kTRUE;
      if (sign == 0) {
	iXdevOff = - iXdevOn + 15;  // gives range  0-14
      } else {
	iXdevOff = + iXdevOn + 15;  // gives range 16-30 !
      }
    }
    iXdev = iXdevOff;

    // iYtri == 1 means no trigger in y-direction
    if (iYtri == 0 && trigx) 
    {
      lut.GetLutOutput(localBoardId,iXpos,iXdev,iYpos,lutLpt,lutHpt);
    }
    
    UChar_t buffer;
    
    // fill byte
    if (i%2 == 0) 
    {
      // upper half-byte
      buffer = 0;	    
      buffer += lutHpt[1] << 7;
      buffer += lutHpt[0] << 6;
      buffer += lutLpt[1] << 5;
      buffer += lutLpt[0] << 4;
    } else {
      // lower half-byte
      buffer += lutHpt[1] << 3;
      buffer += lutHpt[0] << 2;
      buffer += lutLpt[1] << 1;
      buffer += lutLpt[0] << 0;
      fwrite(&buffer,1,1,flut);
    }
  }
}  

//_____________________________________________________________________________
Int_t 
AliMUONTriggerIO::LocalBoardId(Int_t index) const
{  
  /// Return the i-th localBoardId, or -1 if index is out of bounds

  AliMpLocalBoard* board = fRegionalTrigger.GetLocalBoard(index);
  if ( ! board ) return -1;
  
  return board->GetId(); 
}

//_____________________________________________________________________________
void 
AliMUONTriggerIO::UpdateMapping(Bool_t writeFile) const
{
/// check if mapping in database different from current Mtg configuration
/// Update mapping in databse and read regional crate file in repository (ext .out
/// to avoid overwriting). This case has a low probability to happen.

    // Assuming that crates do not change

    if (!AliMpDDLStore::Instance(kFALSE))
    {
      AliMpCDB::LoadDDLStore();
    }

    Bool_t modified = false;

    TExMapIter itr = AliMpDDLStore::Instance()->GetLocalBoardItr();

    Long_t key, value;

    while(itr.Next(key, value))
    {
      AliMpLocalBoard* boardMapping =  reinterpret_cast<AliMpLocalBoard*>(value);

      Int_t localBoardId = boardMapping->GetId();
      AliMpLocalBoard* board = fRegionalTrigger.FindLocalBoard(localBoardId);
      if ( ! board ) {
        AliFatal("Board found in mapping but not in regional trigger");
        return;
      }  

      if ( board->GetCrate().CompareTo(boardMapping->GetCrate()) != 0 ) 
      {
	AliWarning(Form("Crate Name different for board %d (%s %s)", localBoardId, boardMapping->GetCrate().Data(), 
			board->GetCrate().Data()));
	boardMapping->SetCrate( board->GetCrate() );
	modified = true;
      }

      if ((board->GetSlot()) != boardMapping->GetSlot()) 
      {
	AliWarning(Form("Slot different for board %d (%d %d)", localBoardId, boardMapping->GetSlot(), board->GetSlot()+1));
	boardMapping->SetSlot(board->GetSlot());
	modified = true;
      }
	  
      if (boardMapping->GetSwitch() != board->GetSwitch()) {
	AliWarning(Form("Switch different for board %d (0x%x 0x%x)", localBoardId, 
			boardMapping->GetSwitch(), board->GetSwitch()));
	boardMapping->SetSwitch(board->GetSwitch());
	modified = true;
      }
    }
    
    if (modified) 
    {
      AliMpDDLStore::Instance()->SetRegionalTrigger(fRegionalTrigger);
      AliMpCDB::WriteDDLStore(false);
      AliWarning("Wrote new version of mapping in databse");
      if (writeFile) 
      {
          TString file = AliMpFiles::LocalTriggerBoardMapping();
          file += ".out";
          WriteRegional(file.Data(), 0x0);
          AliWarning(Form("Wrote regional file %s", file.Data()));

      }
    }

}
