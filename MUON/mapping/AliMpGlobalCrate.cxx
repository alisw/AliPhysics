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
// Class AliMpGlobalCrate
// --------------------
// The class defines the properties of trigger crate
// Author: Ch. Finck, Subatech Nantes
//-----------------------------------------------------------------------------

#include "AliMpGlobalCrate.h"
#include "AliMpDEManager.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpGlobalCrate)
/// \endcond

const Char_t*  AliMpGlobalCrate::fgkJtagName       = "JtagBoard";
const Char_t*  AliMpGlobalCrate::fgkFirstDarcName  = "LeftDarcBoard";
const Char_t*  AliMpGlobalCrate::fgkSecondDarcName = "RightDarcBoard";
const Char_t*  AliMpGlobalCrate::fgkGlobalName     = "GlobalBoard";  
const Char_t*  AliMpGlobalCrate::fgkFetName        = "FetBoard"; 

const Int_t   AliMpGlobalCrate::fgkGlobalNofRegisters =  13; 
const Int_t   AliMpGlobalCrate::fgkFetNofRegisters    =  7; 
const Int_t   AliMpGlobalCrate::fgkJtagNofLines       =  4; 
//______________________________________________________________________________
AliMpGlobalCrate::AliMpGlobalCrate()
  : TNamed("GlobalCrate", "mapping trigger global crate"),
    fJtagVmeAddr(0x0),
    fJtagClockDiv(0),
    fJtagRxPhase(0),
    fJtagRdDelay(0),
    fEnableJtag(0),
    fJtagCrateName(),                  
    fFirstDarcVmeAddr(0x0),
    fFirstDarcType(0),
    fFirstDarcDisable(0),
    fFirstDarcL0Delay(0),
    fFirstDarcL1TimeOut(0),
    fSecondDarcVmeAddr(0x0),
    fSecondDarcType(0),
    fSecondDarcDisable(0),
    fSecondDarcL0Delay(0),
    fSecondDarcL1TimeOut(0),
    fGlobalVmeAddr(0x0),
    fFetVmeAddr(0x0)
{
/// Standard constructor

  for (Int_t i = 0; i < fgkGlobalNofRegisters; ++i)
    fGlobalRegisters[i] = 0;
  
  for (Int_t j = 0; j < fgkFetNofRegisters; ++j)
    fFetRegisters[j] = 0;
}

//______________________________________________________________________________
AliMpGlobalCrate::AliMpGlobalCrate(TRootIOCtor* /*ioCtor*/)
  : TNamed("GlobalCrate", "mapping trigger global crate"),
    fJtagVmeAddr(0x0),
    fJtagClockDiv(0),
    fJtagRxPhase(0),
    fJtagRdDelay(0),
    fEnableJtag(0),
    fJtagCrateName(),                  
    fFirstDarcVmeAddr(0x0),
    fFirstDarcType(0),
    fFirstDarcDisable(0),
    fFirstDarcL0Delay(0),
    fFirstDarcL1TimeOut(0),
    fSecondDarcVmeAddr(0x0),
    fSecondDarcType(0),
    fSecondDarcDisable(0),
    fSecondDarcL0Delay(0),
    fSecondDarcL1TimeOut(0),
    fGlobalVmeAddr(0x0),
    fGlobalRegisters(),
    fFetVmeAddr(0x0),
    fFetRegisters()
{
/// Root IO constructor
}

//______________________________________________________________________________
AliMpGlobalCrate::~AliMpGlobalCrate()
{
/// Destructor
}

//______________________________________________________________________________
Bool_t AliMpGlobalCrate::GetEnableJtag(Int_t index) const 
{
  /// returns enable mask for a given Jtag line

  if (index > fgkJtagNofLines) {
    AliWarning("Index size too big for Jtag line");
    return kFALSE;
  }
  return ((fEnableJtag >> index) & 0x1);

}

//______________________________________________________________________________
void AliMpGlobalCrate::SetJtagCrateName(Int_t index, TString name) 
{
  /// Get Jtag crate name for a given index 
  if (index > AliMpConstants::LocalBoardNofChannels()) {
    AliWarning("Index size too big for Jtag line");
    return;
  }                                                 
  fJtagCrateName[index] = name;
}

//______________________________________________________________________________
TString AliMpGlobalCrate::GetJtagCrateName(Int_t jtagLine, Int_t index) const 
{ 
  /// Get the crate name for a given line and a given index 
  if (jtagLine > AliMpConstants::LocalBoardNofChannels() || index > AliMpConstants::LocalBoardNofChannels())
    return 0x0;
  else                                       
    return fJtagCrateName[jtagLine*fgkJtagNofLines + index];
}

//______________________________________________________________________________
UInt_t AliMpGlobalCrate::GetGlobalRegister(Int_t index) const       
{
  /// return global register for a given index
  if (index > fgkGlobalNofRegisters) {
    AliWarning("Index size too big for Global Register");
    return 0;
  } else 
    return fGlobalRegisters[index];
}

//______________________________________________________________________________
void AliMpGlobalCrate::SetGlobalRegister(Int_t index, UInt_t reg) 
{
  /// set Global register for a given index
  if (index > fgkGlobalNofRegisters) {
    AliWarning("Index size too big for Global Register");
    return;
  } 
  fGlobalRegisters[index] = reg;
}
   
//______________________________________________________________________________
UInt_t AliMpGlobalCrate::GetFetRegister(Int_t index) const       
{
  /// return global register for a given index
  if (index > fgkFetNofRegisters) {
    AliWarning("Index size too big for Fet Register");
    return 0;
  } else 
    return fFetRegisters[index];
}

//______________________________________________________________________________
void AliMpGlobalCrate::SetFetRegister(Int_t index, UInt_t reg) 
{
  /// set Global register for a given index
  if (index > fgkFetNofRegisters) {
    AliWarning("Index size too big for Global Register");
    return;
  } 
  fFetRegisters[index] = reg;
}
