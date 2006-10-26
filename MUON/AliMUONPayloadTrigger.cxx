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


/// \class AliMUONPayloadTrigger
/// Class Payload
///
/// Decodes rawdata from buffer and stores in TClonesArray.
/// 
/// First version implement for Trigger
///
/// \author Christian Finck

#include "AliMUONPayloadTrigger.h"

#include "AliRawReader.h"
#include "AliRawDataHeader.h"

#ifndef DATE_SYS
#include "AliLog.h"
#endif

#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"

/// \cond CLASSIMP
ClassImp(AliMUONPayloadTrigger)
/// \endcond

AliMUONPayloadTrigger::AliMUONPayloadTrigger()
  : TObject(),
    fMaxReg(8),
    fMaxLoc(16),
    fDDLTrigger(new AliMUONDDLTrigger()),
    fRegHeader(new AliMUONRegHeader()), 
    fLocalStruct(new AliMUONLocalStruct())
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
}


//______________________________________________________
Bool_t AliMUONPayloadTrigger::Decode(UInt_t *buffer)
{
  /// decode trigger DDL
  /// store only non-empty structures (TrigY ==0)

 // reading DDL for trigger

  AliMUONDarcHeader* darcHeader = fDDLTrigger->GetDarcHeader();

  static Int_t kGlobalHeaderSize = darcHeader->GetGlobalHeaderLength(); 
  static Int_t kDarcHeaderSize   = darcHeader->GetDarcHeaderLength(); 
  static Int_t kRegHeaderSize    = fRegHeader->GetHeaderLength();

  Bool_t scalerEvent = kFALSE;
  
  Int_t index = 0;

  memcpy(darcHeader->GetHeader(), &buffer[index], (kDarcHeaderSize)*4); 
  index += kDarcHeaderSize;

  if(darcHeader->GetEventType() == 0) {
    scalerEvent = kTRUE;
  } else
    scalerEvent = kFALSE;

  if(scalerEvent) {
    // 6 DARC scaler words
    memcpy(darcHeader->GetDarcScalers(), &buffer[index], darcHeader->GetDarcScalerLength()*4);
    index += darcHeader->GetDarcScalerLength();
  }

  if (buffer[index++] != darcHeader->GetEndOfDarc())
#ifndef DATE_SYS
    AliWarning(Form("Wrong end of Darc word %x instead of %x\n",buffer[index-1], darcHeader->GetEndOfDarc()));
#else 
  printf("Wrong end of Darc word %x instead of %x\n",buffer[index-1], darcHeader->GetEndOfDarc());
#endif

  // 4 words of global board input + Global board output
  memcpy(darcHeader->GetGlobalInput(), &buffer[index], (kGlobalHeaderSize)*4); 
  index += kGlobalHeaderSize; 

  if(scalerEvent) {
    // 10 Global scaler words
    memcpy(darcHeader->GetGlobalScalers(), &buffer[index], darcHeader->GetGlobalScalerLength()*4);
    index += darcHeader->GetGlobalScalerLength();
  }

  if (buffer[index++] != darcHeader->GetEndOfGlobal())
#ifndef DATE_SYS
    AliWarning(Form("Wrong end of Global word %x instead of %x\n",buffer[index-1], darcHeader->GetEndOfGlobal()));
#else 
  printf("Wrong end of Global word %x instead of %x\n",buffer[index-1], darcHeader->GetEndOfGlobal());
#endif
 
  // 8 regional boards
  for (Int_t iReg = 0; iReg < fMaxReg; iReg++) {           //loop over regeonal card

    memcpy(fRegHeader->GetHeader(), &buffer[index], kRegHeaderSize*4);
    index += kRegHeaderSize;

    fDDLTrigger->AddRegHeader(*fRegHeader);
    // 11 regional scaler word
    if(scalerEvent) {
      memcpy(fRegHeader->GetScalers(), &buffer[index], fRegHeader->GetScalerLength()*4);
      index += fRegHeader->GetScalerLength();
    }

    if (buffer[index++] != fRegHeader->GetEndOfReg())
#ifndef DATE_SYS
      AliWarning(Form("Wrong end of Reg word %x instead of %x\n",buffer[index-1], fRegHeader->GetEndOfReg()));
#else
      printf("Wrong end of Reg word %x instead of %x\n",buffer[index-1], fRegHeader->GetEndOfReg());
#endif

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

      if (buffer[index++] != fLocalStruct->GetEndOfLocal())
#ifndef DATE_SYS
	AliWarning(Form("Wrong end of local word %x instead of %x\n",buffer[index-1], fLocalStruct->GetEndOfLocal()));
#else
      printf("Wrong end of local word %x instead of %x\n",buffer[index-1], fLocalStruct->GetEndOfLocal());
#endif
	  
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
}

//______________________________________________________
void AliMUONPayloadTrigger::SetMaxReg(Int_t reg) 
{
  /// set regional card number
  if (reg > 8) reg = 8;
  fMaxReg = reg;
}

//______________________________________________________
void AliMUONPayloadTrigger::SetMaxLoc(Int_t loc) 
{
  /// set local card number
  if (loc > 16) loc = 16;
  fMaxLoc = loc;
}
