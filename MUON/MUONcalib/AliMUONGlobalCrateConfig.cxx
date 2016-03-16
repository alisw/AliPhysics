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
// $MpId: AliMpTrigger.cxx,v 1.4 2006/05/24 13:58:52 ivana Exp $

//-----------------------------------------------------------------------------
// Class AliMUONGlobalCrateConfig
// --------------------
// The class defines the configuration of trigger crate
// Author: Ch. Finck, Subatech Nantes
//-----------------------------------------------------------------------------

#include "AliMUONGlobalCrateConfig.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpHelper.h"

#include "AliLog.h"

#include <TArrayI.h>
#include <Riostream.h>
#include <TSystem.h>

using std::ifstream;
using std::endl;
using std::ios;
/// \cond CLASSIMP
ClassImp(AliMUONGlobalCrateConfig)
/// \endcond

const Char_t*  AliMUONGlobalCrateConfig::fgkJtagName       = "JtagBoard";
const Char_t*  AliMUONGlobalCrateConfig::fgkFirstDarcName  = "LeftDarcBoard";
const Char_t*  AliMUONGlobalCrateConfig::fgkSecondDarcName = "RightDarcBoard";
const Char_t*  AliMUONGlobalCrateConfig::fgkGlobalName     = "GlobalBoard";  
const Char_t*  AliMUONGlobalCrateConfig::fgkFetName        = "FetBoard"; 

const Int_t   AliMUONGlobalCrateConfig::fgkGlobalNofRegisters =  13; 
const Int_t   AliMUONGlobalCrateConfig::fgkFetNofRegisters    =  7; 
const Int_t   AliMUONGlobalCrateConfig::fgkJtagNofLines       =  4; 
const Int_t   AliMUONGlobalCrateConfig::fgkDarcNofLines       =  8; 
//______________________________________________________________________________
AliMUONGlobalCrateConfig::AliMUONGlobalCrateConfig()
  : TNamed("GlobalCrate", "mapping trigger global crate"),
    fGlobalCrateEnable(0x0), 
    fJtagVmeAddr(0x0),
    fJtagClockDiv(0),
    fJtagRxPhase(0),
    fJtagRdDelay(0),
    fEnableJtag(0),
    fJtagCrateName(),                  
    fFirstDarcCrateName(),                  
    fSecondDarcCrateName(),                  
    fFirstDarcVmeAddr(0x0),
    fFirstDarcType(0),
    fFirstDarcDisable(0),
    fFirstDarcL0Delay(0),
    fFirstDarcL1TimeOut(0),
    fFirstDarcGlobalL0(0),
    fFirstDarcConfig(0), 
    fSecondDarcVmeAddr(0x0),
    fSecondDarcType(0),
    fSecondDarcDisable(0),
    fSecondDarcL0Delay(0),
    fSecondDarcL1TimeOut(0),
    fSecondDarcGlobalL0(0),
    fSecondDarcConfig(0),
    fGlobalVmeAddr(0x0),
    fFetVmeAddr(0x0),
    fEnableFirstDarc(0),
    fEnableSecondDarc(0)
{
/// Standard constructor

  for (Int_t i = 0; i < fgkGlobalNofRegisters; ++i)
    fGlobalRegisters[i] = 0;
  
  for (Int_t j = 0; j < fgkFetNofRegisters; ++j)
    fFetRegisters[j] = 0;
}

//______________________________________________________________________________
AliMUONGlobalCrateConfig::~AliMUONGlobalCrateConfig()
{
/// Destructor
}

//______________________________________________________________________________
Int_t AliMUONGlobalCrateConfig::ReadData(const TString& fileName)
{
    /// Fill trigger global crate object from ascii file
    /// put the method static to be used by other class w/o initializing object
  
    TString inFileName(fileName);
    if ( inFileName == "" )
      inFileName = AliMpFiles::GlobalTriggerBoardMapping();
    
    inFileName = gSystem->ExpandPathName(inFileName.Data());

    ifstream in(inFileName.Data(), ios::in);

    if (!in) {
      AliErrorStream()
         << "Global Trigger Board Mapping File " << fileName.Data() << " not found" << endl;
      return kFALSE;
    }

    Int_t nDarc = 0;

    ULong_t addr;
    UInt_t item;
    TArrayI list;
    Bool_t oldFormat = kFALSE;

    char line[255];
    in.getline(line, 255);
    TString tmp(AliMpHelper::Normalize(line));

    if (!tmp.Contains(GetName()))
        AliWarning("Wrong Global Crate File");

    // enable
    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    UInt_t en = 0;
    sscanf(tmp.Data(), "%x", &en);
    SetGlobalCrateEnable(en);

    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    
    if (tmp.Contains(GetJtagName())) {

      // old version, old regional boards

      oldFormat = kTRUE;
      AliInfo(Form("Old format of global config"));

    } else {

      if (tmp.Contains(GetFirstDarcName())) {
	
	// new version, new regional boards

	oldFormat = kFALSE;
	AliInfo(Form("New format of global config"));

      }

    }

    if (oldFormat) {

      // vme addr
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      ULong_t addr;
      sscanf(tmp.Data(), "%lx", &addr);
      SetJtagVmeAddr(addr);
      AliDebug(1, Form("Jtag Vme Address: 0x%lx", addr));
      
      // clk div, rx phase, read delay
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      TArrayI list;
      AliMpHelper::DecodeName(line, ' ', list);
      SetJtagClockDiv(list[0]);
      SetJtagRxPhase(list[1]);
      SetJtagRdDelay(list[2]);
      AliDebug(1, Form("Jtag Clock Div: %d, Rx Phase: %d, Read Delay %d", list[0], list[1], list[2]));
      
      // enable
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      AliMpHelper::DecodeName(line, ' ', list);
      UChar_t enable = 0;
      for (Int_t i = 0; i < GetJtagNofLines(); ++i)
	enable |= (list[i] << i);
      SetEnableJtag(enable);
      AliDebug(1, Form("Jtag Enable: 0x%x", enable));
      
      for (Int_t i = 0; i < GetJtagNofLines(); ++i) {
	in.getline(line, 255);
	for (Int_t j = 0; j < GetJtagNofLines(); ++j) {
	  in.getline(line, 255);
	  tmp = AliMpHelper::Normalize(line);
	  SetJtagCrateName(i*GetJtagNofLines() + j, tmp);
	  //AliDebug(1, Form("Jtag Crate Name: %s", tmp.Data()));
	}
      }
      
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);

    } // end old format
    
    // vme addr
    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    sscanf(tmp.Data(), "%lx", &addr);
    if (addr) nDarc++;
    SetFirstDarcVmeAddr(addr);
    AliDebug(1, Form("First Darc Vme Address: 0x%lx", addr));
    
    // type
    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    SetFirstDarcType(tmp.Atoi());
    AliDebug(1, Form("First Darc Type: %d", tmp.Atoi()));
    
    // disable
    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    sscanf(tmp.Data(), "%x", &item);
    SetFirstDarcDisable(item);
    AliDebug(1, Form("First Darc Disable: 0x%x", item));
    
    // L0
    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    sscanf(tmp.Data(), "%x", &item);
    SetFirstDarcL0Delay(item);
    AliDebug(1, Form("First Darc L0 Delay: 0x%x", item));
    
    // L1
    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    sscanf(tmp.Data(), "%x", &item);
    SetFirstDarcL1TimeOut(item);
    AliDebug(1, Form("First Darc L1 Time Out: 0x%x", item));
    
    // Global L0 delay
    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    sscanf(tmp.Data(), "%x", &item);
    SetFirstDarcGlobalL0(item);
    AliDebug(1, Form("First Darc Global L0 delay: 0x%x", item));
    
    // Trigger configuration
    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    sscanf(tmp.Data(), "%x", &item);
    SetFirstDarcConfig(item);
    AliDebug(1, Form("First Darc Config: 0x%x", item));    
    
    if (!oldFormat) {

      // enable
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      AliMpHelper::DecodeName(line, ' ', list);
      UChar_t enable = 0;
      for (Int_t i = 0; i < GetDarcNofLines(); ++i)
	enable |= (list[i] << i);
      SetEnableFirstDarc(enable);
      AliDebug(1, Form("First Darc Enable: 0x%x", enable));
      
      for (Int_t i = 0; i < GetDarcNofLines(); ++i) {
	in.getline(line, 255);
	tmp = AliMpHelper::Normalize(line);
	SetFirstDarcCrateName(i, tmp);
	//AliInfo(Form("First Darc Crate Name: %s", tmp.Data()));
      }

    } // end new format First Darc

    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    if (tmp.Contains(GetSecondDarcName())) {

      // vme addr
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      sscanf(tmp.Data(), "%lx", &addr);
      if (addr) nDarc++;
      SetSecondDarcVmeAddr(addr);
      AliDebug(1, Form("Second Darc Vme Address: 0x%lx", addr));
      
      // type
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      SetSecondDarcType(tmp.Atoi());
      AliDebug(1, Form("Second Darc Type: %d", tmp.Atoi()));
      
      // enable
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      sscanf(tmp.Data(), "%x", &item);
      SetSecondDarcDisable(item);
      AliDebug(1, Form("Second Darc Disable: 0x%x", item));
      
      // L0
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      sscanf(tmp.Data(), "%x", &item);
      SetSecondDarcL0Delay(item);
      AliDebug(1, Form("Second Darc L0 Delay: 0x%x", item));
      
      // L1
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      sscanf(tmp.Data(), "%x", &item);
      SetSecondDarcL1TimeOut(item);
      AliDebug(1, Form("Second Darc L1 Time Out: 0x%x", item));
      
      // Global L0 delay
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      sscanf(tmp.Data(), "%x", &item);
      SetSecondDarcGlobalL0(item);
      AliDebug(1, Form("Second Darc Global L0 delay: 0x%x", item));
      
      // Trigger configuration
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      sscanf(tmp.Data(), "%x", &item);
      SetSecondDarcConfig(item);
      AliDebug(1, Form("Second Darc Config: 0x%x", item));    
    
    } // end Second Darc old and new format

    if (!oldFormat) {

      // enable
      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      AliMpHelper::DecodeName(line, ' ', list);
      UChar_t enable = 0;
      for (Int_t i = 0; i < GetDarcNofLines(); ++i)
	enable |= (list[i] << i);
      SetEnableSecondDarc(enable);
      AliDebug(1, Form("Second Darc Enable: 0x%x", enable));
      
      for (Int_t i = 0; i < GetDarcNofLines(); ++i) {
	in.getline(line, 255);
	tmp = AliMpHelper::Normalize(line);
	SetSecondDarcCrateName(i, tmp);
	//AliInfo(Form("Second Darc Crate Name: %s", tmp.Data()));
      }

    } // end new format Second Darc

    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    if (tmp.Contains(GetGlobalName())) {

      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      sscanf(tmp.Data(), "%lx", &addr);
      SetGlobalVmeAddr(addr);
      AliDebug(1, Form("Global Vme Address: 0x%lx", addr));
      
      for (Int_t i = 0; i < GetGlobalNofRegisters(); ++i) {
	in.getline(line, 255);
	tmp = AliMpHelper::Normalize(line);
	UInt_t reg;
	sscanf(tmp.Data(), "%x", &reg);
	SetGlobalRegister(i, reg);
	AliDebug(1, Form("Global Register %d: 0x%x", i, reg));
      }
    
    } // end Global board old and new format

    in.getline(line, 255);
    tmp = AliMpHelper::Normalize(line);
    if (tmp.Contains(GetFetName())) {

      in.getline(line, 255);
      tmp = AliMpHelper::Normalize(line);
      sscanf(tmp.Data(), "%lx", &addr);
      SetFetVmeAddr(addr);
      AliDebug(1, Form("Fet Vme Address: 0x%lx", addr));
      
      for (Int_t i = 0; i < GetFetNofRegisters(); ++i) {
	in.getline(line, 255);
	tmp = AliMpHelper::Normalize(line);
	UInt_t reg;
	sscanf(tmp.Data(), "%x", &reg);
	SetFetRegister(i, reg);
	AliDebug(1, Form("Fet Register %d: 0x%x", i, reg));
      }

    } // end Fet board old and new format

    return nDarc;
}

//______________________________________________________________________________
Bool_t AliMUONGlobalCrateConfig::GetEnableJtag(Int_t index) const 
{
  /// returns enable mask for a given Jtag line

  if (index > fgkJtagNofLines) {
    AliWarning("Index size too big for Jtag line");
    return kFALSE;
  }
  return ((fEnableJtag >> index) & 0x1);

}

//______________________________________________________________________________
void AliMUONGlobalCrateConfig::SetJtagCrateName(Int_t index, TString name) 
{
  /// Set Jtag crate name for a given index 
  if (index > AliMpConstants::LocalBoardNofChannels()) {
    AliWarning("Index size too big for Jtag line");
    return;
  }                                                 
  fJtagCrateName[index] = name;
}

//______________________________________________________________________________
TString AliMUONGlobalCrateConfig::GetJtagCrateName(Int_t jtagLine, Int_t index) const 
{ 
  /// Get the crate name for a given line and a given index 
  if (jtagLine > AliMpConstants::LocalBoardNofChannels() || index > AliMpConstants::LocalBoardNofChannels())
    return "";
  else                                       
    return fJtagCrateName[jtagLine*fgkJtagNofLines + index];
}

//______________________________________________________________________________
Bool_t AliMUONGlobalCrateConfig::GetEnableFirstDarc(Int_t index) const 
{
  /// returns enable mask for a given First Darc line

  if (index >= fgkDarcNofLines) {
    AliWarning("Index size too big for First Darc line");
    return kFALSE;
  }
  return ((fEnableFirstDarc >> index) & 0x1);

}

//______________________________________________________________________________
void AliMUONGlobalCrateConfig::SetFirstDarcCrateName(Int_t index, TString name) 
{
  /// Set First Darc crate name for a given index 
  if (index >= AliMpConstants::LocalBoardNofChannels()/2) {
    AliWarning("Index size too big for First Darc line");
    return;
  }                                                 
  fFirstDarcCrateName[index] = name;
}

//______________________________________________________________________________
TString AliMUONGlobalCrateConfig::GetFirstDarcCrateName(Int_t index) const 
{ 
  /// Get the First Darc crate name for a given index
  if (index >= AliMpConstants::LocalBoardNofChannels()/2)
    return "";
  else                                       
    return fFirstDarcCrateName[index];
}

//______________________________________________________________________________
Bool_t AliMUONGlobalCrateConfig::GetEnableSecondDarc(Int_t index) const 
{
  /// returns enable mask for a given Second Darc line

  if (index >= fgkDarcNofLines) {
    AliWarning("Index size too big for Second Darc line");
    return kFALSE;
  }
  return ((fEnableSecondDarc >> index) & 0x1);

}

//______________________________________________________________________________
void AliMUONGlobalCrateConfig::SetSecondDarcCrateName(Int_t index, TString name) 
{
  /// Set Second Darc crate name for a given index 
  if (index >= AliMpConstants::LocalBoardNofChannels()/2) {
    AliWarning("Index size too big for Second Darc line");
    return;
  }                                                 
  fSecondDarcCrateName[index] = name;
}

//______________________________________________________________________________
TString AliMUONGlobalCrateConfig::GetSecondDarcCrateName(Int_t index) const 
{ 
  /// Get the Second Darc crate name for a given index
  if (index >= AliMpConstants::LocalBoardNofChannels()/2)
    return "";
  else                                       
    return fSecondDarcCrateName[index];
}

//______________________________________________________________________________
UInt_t AliMUONGlobalCrateConfig::GetGlobalRegister(Int_t index) const       
{
  /// return global register for a given index
  if (index >= fgkGlobalNofRegisters) {
    AliWarning("Index size too big for Global Register");
    return 0;
  } else 
    return fGlobalRegisters[index];
}

//______________________________________________________________________________
void AliMUONGlobalCrateConfig::SetGlobalRegister(Int_t index, UInt_t reg) 
{
  /// set Global register for a given index
  if (index >= fgkGlobalNofRegisters) {
    AliWarning("Index size too big for Global Register");
    return;
  } 
  fGlobalRegisters[index] = reg;
}
   
//______________________________________________________________________________
void AliMUONGlobalCrateConfig::SetGlobalMask(Int_t index, UInt_t mask)  
{
  /// set one word of the global mask

  if (index >= 0 && index < 4) {
    SetGlobalRegister(index,mask);
  } else {
    AliWarning(Form("Check register number of the mask (%d) \n",index));
  }

}

//______________________________________________________________________________
UInt_t AliMUONGlobalCrateConfig::GetGlobalMask(Int_t index) const       
{
  /// return one word of the global mask
  if (index >= 0 && index < 4) {
    return fGlobalRegisters[index];
  } else {
    AliWarning(Form("Check register number of the mask (%d) \n",index));
    return 0;
  }
}

//______________________________________________________________________________
Bool_t AliMUONGlobalCrateConfig::GetMasksOn() const       
{
  /// indicates if global masks are active on global inputs

  // test 7th lsb
  if (fGlobalRegisters[4] & 0x40) return kTRUE;

  return kFALSE;

}

//______________________________________________________________________________
UInt_t AliMUONGlobalCrateConfig::GetFetRegister(Int_t index) const       
{
  /// return global register for a given index
  if (index >= fgkFetNofRegisters) {
    AliWarning("Index size too big for Fet Register");
    return 0;
  } else 
    return fFetRegisters[index];
}

//______________________________________________________________________________
void AliMUONGlobalCrateConfig::SetFetRegister(Int_t index, UInt_t reg) 
{
  /// set Global register for a given index
  if (index >= fgkFetNofRegisters) {
    AliWarning("Index size too big for Global Register");
    return;
  } 
  fFetRegisters[index] = reg;
}
