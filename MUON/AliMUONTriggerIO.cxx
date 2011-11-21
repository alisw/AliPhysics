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
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONTriggerCrateConfig.h"

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


const UInt_t AliMUONTriggerIO::fgkLocalLutSize = 1 << 14; // 16384

//_____________________________________________________________________________
AliMUONTriggerIO::AliMUONTriggerIO() 
    : TObject(), 
      fRegionalTrigger()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONTriggerIO::AliMUONTriggerIO(const char* regionalFileToRead) 
    :TObject(), 
     fRegionalTrigger()
{
  /// ctor
  ReadRegionalConfig(regionalFileToRead,0);
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
  
  Int_t localBoardIndex(0);
    
  while ( fread ( maskBuffer, 2, 8, fp ) )
  {
    Int_t localBoardId = fRegionalTrigger.LocalBoardId(localBoardIndex);
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
    
    if ( localBoardId > 0 ) 
    {
      AliMUONVCalibParam* localBoard = new AliMUONCalibParamNI(1,8,localBoardId,0,0);
      for ( Int_t index = 0; index < 8; ++index )
      {
        localBoard->SetValueAsInt(index,0,maskBuffer[index]);
      }
      localMasks.Add(localBoard);
    }
    else
    {
      AliError(Form("Oups. Got localBoardId=%d for index=%d",localBoardId,localBoardIndex));
    }
    
    ++localBoardIndex;
  }
  
  if ( localBoardIndex != NofLocalBoards() ) 
  {
    AliError(Form("Read %d out of %d local boards",
                  localBoardIndex, NofLocalBoards()));
  }
  
  fclose(fp);
  
  return localBoardIndex+1;
}

//_____________________________________________________________________________
void
AliMUONTriggerIO::ReadLocalLUT(AliMUONTriggerLut& lut,
                               Int_t localBoardId,
                               FILE* flut)
{
  /// Read the LUT for one local board from an online file

  UShort_t address;
  
  UChar_t buffer[fgkLocalLutSize];   // 32768 hpt/lpt addresses divided by two
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
  
  // read two lut addresses at once, 32768/2=16384 times
  if (fread(buffer,fgkLocalLutSize,1,flut) == 0) {
    AliWarning("Error reading the LUT file");
    return;
  }

  // create the 32767 addresses for the 4-bits lpt and hpt half-bytes
  for (UShort_t ilut = 0; ilut < fgkLocalLutSize*2; ilut += 2) 
  {
    
    // 1st 4-bits half-byte
    address = ilut;   
    lh = (buffer[ilut/2] & mask1) >> 4;
    
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
    lh = (buffer[ilut/2] & mask2);
    
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
    ReadLocalLUT(lut,fRegionalTrigger.LocalBoardId(i),flut);
  }
  
  fclose(flut);
  
  return kTRUE;
  
}

//_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::ReadConfig(const char* localFile,
                             const char* regionalFile,
                             const char* globalFile,
                             AliMUONVStore* localMasks,
                             AliMUONRegionalTriggerConfig* regionalConfig,
                             AliMUONGlobalCrateConfig* globalConfig)
{
  /// Fill the various masks store from files
  
  if ( !regionalConfig || !regionalFile || strlen(regionalFile)==0 ) 
  {
    AliError("Must have a regional file name to proceeed");
    return kFALSE;
  }
  
  AliDebug(1,Form("regionalConfig=%p",regionalConfig));
  
  Int_t nCrates = ReadRegionalConfig(regionalFile, regionalConfig);

  if (!nCrates) 
  {
    AliError("nCrates=0 !");
    return kFALSE;
  }
  
  if (localMasks && localFile && strlen(localFile) > 0 )
  {
    Int_t nLocal = ReadLocalMasks(localFile,*localMasks);
    AliDebug(1,Form("Read masks for %d local boards",nLocal));
  }
  
  Int_t nDarc = ReadGlobalConfig(globalFile, globalConfig);
  AliDebug(1,Form("Read config for %d DARC boards",nDarc));
  
  if (!nDarc) return kFALSE;
  
  return kTRUE;
}



//_____________________________________________________________________________
 Int_t 
 AliMUONTriggerIO::ReadGlobalConfig(const char* globalFile, AliMUONGlobalCrateConfig* globalConfig) const
{
  /// read the global crate file
  /// the masks are disable bit for each crate, 8 per darc board
  /// bit value 0 means enable, 1 means disable                                                 * 
  
  Int_t nDarc = 0;
  if ( !(nDarc = globalConfig->ReadData(globalFile)) ) return 0;
  
  return nDarc;
}
 
//_____________________________________________________________________________
Int_t
AliMUONTriggerIO::ReadRegionalConfig(const char* regionalFile, AliMUONRegionalTriggerConfig* regionalConfig)
{
  /// Read regional file to fill  
  
  AliDebug(1,Form("regionalConfig=%p",regionalConfig));
  
  Int_t nCrates = 0;
  if ( !(nCrates = regionalConfig->ReadData(regionalFile)) ) return 0;

  // read the mapping file also
  if ( ! fRegionalTrigger.ReadData(regionalFile) ) return 0;

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
    WriteLocalLUT(lut,fRegionalTrigger.LocalBoardId(i),flut);
  }
  
  fclose(flut);
  
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::WriteConfig(const char* localFile,
			     const char* regionalFile,
			     const char* globalFile,
			     const AliMUONVStore* localMasks,
                    AliMUONRegionalTriggerConfig* regionalConfig,
                    AliMUONGlobalCrateConfig* globalConfig) const
{
/// write config files

    Bool_t ok;
    ok  = WriteLocalMasks(localFile, *localMasks);
    ok &= WriteRegionalConfig(regionalFile, regionalConfig);
    ok &= WriteGlobalConfig(globalFile, globalConfig);
    
    return ok;


}

 
 //_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::WriteGlobalConfig(const char* globalFile, AliMUONGlobalCrateConfig* globalConfig) const
{
    /// write global config

  ofstream out;
  Int_t disable = 0;
  
  out.open(globalFile);
  if (!out.good())
  {
    AliError(Form("Could not create output global file %s", globalFile));
    return kFALSE;
  }
   
  out << globalConfig->GetName() << endl;
  out << Form("0x%x",globalConfig->GetGlobalCrateEnable()) << endl;
  
  // Jtag
  out << globalConfig->GetJtagName() << endl;
  out << Form("0x%08lx", globalConfig->GetJtagVmeAddr()) << endl;
  out << Form("%d %d %d", globalConfig->GetJtagClockDiv(), 
              globalConfig->GetJtagRxPhase(), globalConfig->GetJtagRdDelay()) << endl;
 
  for (Int_t i = 0; i < globalConfig->GetJtagNofLines(); ++i)
    out << Form("%d ", globalConfig->GetEnableJtag(i));
  out << endl;

  
  for (Int_t i = 0; i < globalConfig->GetJtagNofLines(); ++i)
  {
    out << i << endl;
    for (Int_t j = 0; j < globalConfig->GetJtagNofLines(); ++j)
      out << Form(" %s", globalConfig->GetJtagCrateName(i,j).Data()) << endl;
  }
  
  // first darc board
  out << globalConfig->GetFirstDarcName() << endl;
  out << Form("0x%08lx", globalConfig->GetFirstDarcVmeAddr()) << endl;
  out << globalConfig->GetFirstDarcType() << endl;
  disable = globalConfig->GetFirstDarcDisable();
  out << Form("0x%02x", disable) << endl;
  out << Form("0x%x", globalConfig->GetFirstDarcL0Delay()) << endl;
  out << Form("0x%x", globalConfig->GetFirstDarcL1TimeOut()) << endl;
  out << Form("0x%x", globalConfig->GetFirstDarcGlobalL0()) << endl;
  out << Form("0x%x", globalConfig->GetFirstDarcConfig()) << endl;
  
  // second darc board
  out << globalConfig->GetSecondDarcName() << endl;
  out << Form("0x%08lx", globalConfig->GetSecondDarcVmeAddr()) << endl;
  out << globalConfig->GetSecondDarcType() << endl;
  disable = globalConfig->GetSecondDarcDisable();
  out << Form("0x%02x", disable) << endl;
  out << Form("0x%x", globalConfig->GetSecondDarcL0Delay()) << endl;
  out << Form("0x%x", globalConfig->GetSecondDarcL1TimeOut()) << endl; 
  out << Form("0x%x", globalConfig->GetSecondDarcGlobalL0()) << endl; 
  out << Form("0x%x", globalConfig->GetSecondDarcConfig()) << endl; 
  
  // global board
  out << globalConfig->GetGlobalName() << endl;
  out << Form("0x%08lx", globalConfig->GetGlobalVmeAddr()) << endl;
  for (Int_t i = 0; i < globalConfig->GetGlobalNofRegisters(); ++i)
    out << Form("0x%x", globalConfig->GetGlobalRegister(i)) << endl;
  
  // Fet board
  out << globalConfig->GetFetName() << endl;
  out << Form("0x%08lx", globalConfig->GetFetVmeAddr()) << endl;
  for (Int_t i = 0; i < globalConfig->GetFetNofRegisters(); ++i)
    out << Form("0x%x", globalConfig->GetFetRegister(i)) << endl;
  
  return kTRUE;
}
 
//_____________________________________________________________________________
Bool_t
AliMUONTriggerIO::WriteRegionalConfig(const char* regionalFile, AliMUONRegionalTriggerConfig* regionalConfig) const
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

    Int_t nCrate = fRegionalTrigger.GetNofTriggerCrates();
    if (!nCrate)
    {
      AliError("Could not write regional no configuration in memory");
      return kFALSE;
    }

    Int_t nofDDLs = 0;
    TString name;
    AliMpTriggerCrate* crate;
    for (Int_t ddlId = 0; ddlId < 2; ddlId++) // right & left side            
      {
	for (Int_t crateId = 0; crateId < 8; crateId++) // 8 crates/regional boards for each side.
	  {
	    
	    name = AliMpTriggerCrate::GenerateName(crateId, ddlId, nofDDLs);
	    
	    crate = fRegionalTrigger.FindTriggerCrate(name, false);
	    
	    AliMUONTriggerCrateConfig* crateConfig = regionalConfig->FindTriggerCrate(crate->GetName());
	    if (!crateConfig) 
	      {
		AliError(Form("Cannot find crate %s in CDB", crate->GetName()));
		return kFALSE;
	      }
	    
	    out << crate->GetName()  << endl;
	    out << Form("%02x", crate->GetId())   << endl;
	    out << crateConfig->GetMode()  << endl;
	    out << crateConfig->GetCoinc() << endl;
	    out << Form("%04x", crateConfig->GetMask()) << endl;
	    out << Form("%02d",crate->GetNofLocalBoards()) << endl;
	    
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
		
		if (board->IsNotified()) {
		  for (Int_t i = 0; i < board->GetNofDEs(); ++i)
		    out << Form("%4d ", board->GetDEId(i));
		} else {
		  out << Form("%4d ", 0);
		}
		out << endl;
		
		// print copy card numbers & TC
		out << Form(" %4d %4d", board->GetInputXfrom(), board->GetInputXto());
		out << Form(" %4d %4d", board->GetInputYfrom(), board->GetInputYto());
		out << Form(" %4d",     board->GetTC()) << endl;
	      }
	  }
      }

    out.close();
    
    return kTRUE;
}


//_____________________________________________________________________________
Bool_t 
AliMUONTriggerIO::WriteLocalMasks(const char* localFile, const AliMUONVStore& localMasks) const
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
    Int_t localBoardIndex(0);
    while (localBoardIndex < NofLocalBoards()) {

      Int_t localBoardId = fRegionalTrigger.LocalBoardId(localBoardIndex);

      AliMUONVCalibParam* localMask = 
	static_cast<AliMUONVCalibParam*>(localMasks.FindObject(localBoardId));

      for (Int_t index = 0; index < 8; ++index) 
	{
	  maskBuffer[index] = localMask->ValueAsInt(index,0); 
	}
      
      fwrite ( maskBuffer, 2, 8, fp); 

      ++localBoardIndex;

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

  UChar_t buffer[fgkLocalLutSize];  // 32768 hpt/lpt addresses divided by two
  Int_t bc = 0;
  
  for (UInt_t i = 0; i < fgkLocalLutSize*2; ++i) 
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
    
    // fill byte
    if (i%2 == 0) 
    {
      // upper half-byte
      buffer[bc] = 0;	    
      buffer[bc] += lutHpt[1] << 7;
      buffer[bc] += lutHpt[0] << 6;
      buffer[bc] += lutLpt[1] << 5;
      buffer[bc] += lutLpt[0] << 4;
    } else {
      // lower half-byte
      buffer[bc] += lutHpt[1] << 3;
      buffer[bc] += lutHpt[0] << 2;
      buffer[bc] += lutLpt[1] << 1;
      buffer[bc] += lutLpt[0] << 0;
      bc++;
    }
  }
  fwrite(&buffer,bc,1,flut);
}  

//_____________________________________________________________________________
Int_t 
AliMUONTriggerIO::LocalBoardId(Int_t index) const
{  
  /// Return the i-th localBoardId, or -1 if index is out of bounds

  return fRegionalTrigger.LocalBoardId(index);
}


//______________________________________________________________________________

Int_t AliMUONTriggerIO::LocalBoardId(Int_t ddlId, Int_t crateId, Int_t localId) const
{
    /// Return local board id from crate and local indexes.
    
    Int_t nofDDLs = 0;
    TString name = AliMpTriggerCrate::GenerateName(crateId, ddlId, nofDDLs);

    AliMpTriggerCrate* crate = fRegionalTrigger.FindTriggerCrate(name, false);
    return crate->GetLocalBoardId(localId);
}
