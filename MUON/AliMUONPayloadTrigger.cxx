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

//-----------------------------------------------------------------------------
/// \class AliMUONPayloadTrigger
/// Class Payload
///
/// Decodes rawdata from buffer and stores in TClonesArray.
/// 
/// First version implement for Trigger
///
/// \author Christian Finck
//-----------------------------------------------------------------------------

#include "AliMUONPayloadTrigger.h"

#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"
#include "AliMUONLogger.h"

#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONPayloadTrigger)
/// \endcond

AliMUONPayloadTrigger::AliMUONPayloadTrigger()
  : TObject(),
    fMaxReg(8),
    fMaxLoc(16),
    fDDLTrigger(new AliMUONDDLTrigger()),
    fRegHeader(new AliMUONRegHeader()), 
    fLocalStruct(new AliMUONLocalStruct()),
    fLog(new AliMUONLogger(1000)),
    fDarcEoWErrors(0),
    fGlobalEoWErrors(0),
    fRegEoWErrors(0),
    fLocalEoWErrors(0),
    fWarnings(kTRUE),
    fNofRegSet(kFALSE)
{
  ///
  /// create an object to read MUON raw digits
  /// Default ctor for monitoring purposes
  ///

}

//___________________________________
AliMUONPayloadTrigger::~AliMUONPayloadTrigger()
{
  ///
  /// clean up
  ///
  delete fDDLTrigger;
  delete fLocalStruct;
  delete fRegHeader;
  delete fLog;
}


//______________________________________________________
Bool_t AliMUONPayloadTrigger::Decode(UInt_t *buffer, Bool_t scalerEvent)
{
  /// decode trigger DDL
  /// store only notified cards

 // reading DDL for trigger

  AliMUONDarcHeader* darcHeader = fDDLTrigger->GetDarcHeader();

  static Int_t kGlobalHeaderSize   = darcHeader->GetGlobalHeaderLength(); 
  static Int_t kDarcHeaderSize     = darcHeader->GetDarcHeaderLength(); 
  static Int_t kRegHeaderSize      = fRegHeader->GetHeaderLength();
  static Int_t kRegEmptySize       = fRegHeader->GetHeaderLength()+1 + 16*(fLocalStruct->GetLength()+1);
  static Int_t kRegEmptyScalerSize = fRegHeader->GetHeaderLength() + fRegHeader->GetScalerLength() + 1 +
                                      16*(fLocalStruct->GetLength() + fLocalStruct->GetScalerLength() + 1);

  Int_t index = 0;

  memcpy(darcHeader->GetHeader(), &buffer[index], (kDarcHeaderSize)*4); 
  index += kDarcHeaderSize;


  // darc type vadorh
  if (darcHeader->GetDarcType() == darcHeader->GetDarcVadohrType())
      fMaxReg = 1;
    
  // darc type def.
  if (darcHeader->GetDarcType() == darcHeader->GetDarcDefaultType())
      fMaxReg = 8;
      
  if(darcHeader->GetEventType() == scalerEvent) 
      if (fWarnings) AliWarning("Wrong event type obtained from the Darc header, take the one of CDH");


  if(scalerEvent) {
    // 6 DARC scaler words
    memcpy(darcHeader->GetDarcScalers(), &buffer[index], darcHeader->GetDarcScalerLength()*4);
    index += darcHeader->GetDarcScalerLength();
  }

  if (buffer[index++] != darcHeader->GetEndOfDarc()) {

      const Char_t* msg = Form("Wrong end of Darc word %x instead of %x\n",
		    buffer[index-1], darcHeader->GetEndOfDarc());
      if (fWarnings) AliWarning(msg);
      AddErrorMessage(msg);
      fDarcEoWErrors++;
  }
  // 4 words of global board input + Global board output
  memcpy(darcHeader->GetGlobalInput(), &buffer[index], (kGlobalHeaderSize)*4); 
  index += kGlobalHeaderSize; 

  if(scalerEvent) {
    // 10 Global scaler words
    memcpy(darcHeader->GetGlobalScalers(), &buffer[index], darcHeader->GetGlobalScalerLength()*4);
    index += darcHeader->GetGlobalScalerLength();
  }

  if (buffer[index++] != darcHeader->GetEndOfGlobal()) {

      const Char_t* msg = Form("Wrong end of Global word %x instead of %x\n",
		      buffer[index-1], darcHeader->GetEndOfGlobal());
      if (fWarnings) AliWarning(msg);
      AddErrorMessage(msg);
      fGlobalEoWErrors++;
  }
  // 8 regional boards
  for (Int_t iReg = 0; iReg < fMaxReg; iReg++) {           //loop over regeonal card

    // skip empty regaional board (not connected or with error reading)
    if (buffer[index] == fRegHeader->GetErrorWord()) {
      fDDLTrigger->AddRegHeader(*fRegHeader);
      if (scalerEvent)
        index += kRegEmptyScalerSize;
      else 
        index += kRegEmptySize;
      continue;
    }
    memcpy(fRegHeader->GetHeader(), &buffer[index], kRegHeaderSize*4);
    index += kRegHeaderSize;

    fDDLTrigger->AddRegHeader(*fRegHeader);
    // 11 regional scaler word
    if(scalerEvent) {
      memcpy(fRegHeader->GetScalers(), &buffer[index], fRegHeader->GetScalerLength()*4);
      index += fRegHeader->GetScalerLength();
    }

    if (buffer[index++] != fRegHeader->GetEndOfReg()) {

      const Char_t* msg = Form("Wrong end of Regional word %x instead of %x\n",
		    buffer[index-1], fRegHeader->GetEndOfReg());
      if (fWarnings) AliWarning(msg);
      AddErrorMessage(msg);
      fRegEoWErrors++;
    }
    // 16 local cards per regional board
    for (Int_t iLoc = 0; iLoc < fMaxLoc; iLoc++) {         //loop over local card
	  
      Int_t dataSize = fLocalStruct->GetLength();;

      // 5 word trigger information
      memcpy(fLocalStruct->GetData(), &buffer[index], dataSize*4); 
      index += dataSize;	 

      // 45 regional scaler word
      if(scalerEvent) {
	memcpy(fLocalStruct->GetScalers(), &buffer[index], fLocalStruct->GetScalerLength()*4);
	index += fLocalStruct->GetScalerLength();
      }

      if (buffer[index++] != fLocalStruct->GetEndOfLocal()) {

        const Char_t* msg = Form("Wrong end of Local word %x instead of %x\n",
                                 buffer[index-1], fLocalStruct->GetEndOfLocal());
        
        if (fWarnings) AliWarning(msg);
        AddErrorMessage(msg);
	fLocalEoWErrors++;
      }
      // fill only if card notified
      if (fLocalStruct->GetData(0) == fLocalStruct->GetDisableWord())
	  continue;

      fDDLTrigger->AddLocStruct(*fLocalStruct, iReg);

    } // local card loop
	
  } // regional card loop
      

  return kTRUE;
}

//______________________________________________________
void AliMUONPayloadTrigger::ResetDDL()
{
  /// reseting TClonesArray
  /// after each DDL
  ///
  AliMUONDarcHeader* darcHeader = fDDLTrigger->GetDarcHeader();
  darcHeader->GetRegHeaderArray()->Delete();
  fDarcEoWErrors   = 0;
  fGlobalEoWErrors = 0;
  fRegEoWErrors    = 0;
  fLocalEoWErrors  = 0;
}

//______________________________________________________
void AliMUONPayloadTrigger::SetMaxReg(Int_t reg) 
{
  /// set regional card number
  if (reg > 8) reg = 8;
   fMaxReg = reg;
   
  fNofRegSet = kTRUE;
}

//______________________________________________________
void AliMUONPayloadTrigger::SetMaxLoc(Int_t loc) 
{
  /// set local card number
  if (loc > 16) loc = 16;
  fMaxLoc = loc;
}

//______________________________________________________
void AliMUONPayloadTrigger::AddErrorMessage(const Char_t* msg)
{
/// adding message to logger
 
  TString tmp(msg);
  
  Int_t pos = tmp.First("\n");
  tmp[pos] = 0;
    
  fLog->Log(tmp.Data());
}

