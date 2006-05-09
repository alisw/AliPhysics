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


///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to MUON digits in raw data.
///
/// It loops over all MUON digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE(under develpment).
/// It can loop also over DDL and store the decoded rawdata in TClonesArrays
/// 
/// First version implement for Trigger
///
///////////////////////////////////////////////////////////////////////////////

#include "AliMUONRawStreamTrigger.h"

#include "AliRawReader.h"
#include "AliRawDataHeader.h"
#include "AliLog.h"

#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"

ClassImp(AliMUONRawStreamTrigger)

AliMUONRawStreamTrigger::AliMUONRawStreamTrigger()
  : TObject(),
    fRawReader(0x0),
    fDDL(0),
    fSubEntries(0),
    fNextDDL(kTRUE),
    fMaxDDL(2),
    fMaxReg(8),
    fMaxLoc(16)
{
  //
  // create an object to read MUON raw digits
  // Default ctor for monitoring purposes
  //

  fDDLTrigger  = new AliMUONDDLTrigger();
  fRegHeader   = new AliMUONRegHeader();  
  fLocalStruct = new AliMUONLocalStruct();

}

//_________________________________________________________________
AliMUONRawStreamTrigger::AliMUONRawStreamTrigger(AliRawReader* rawReader)
  : TObject(),
    fDDL(0),
    fSubEntries(0),
    fNextDDL(kTRUE),
    fMaxDDL(2),
    fMaxReg(8),
    fMaxLoc(16)
{
  //
  // ctor with AliRawReader as argument
  // for reconstruction purpose
  //

  fRawReader   = rawReader;

  fDDLTrigger  = new AliMUONDDLTrigger();
  fRegHeader   = new AliMUONRegHeader();
  fLocalStruct = new AliMUONLocalStruct();
  
}

//_________________________________________________________________
AliMUONRawStreamTrigger::AliMUONRawStreamTrigger(const AliMUONRawStreamTrigger& stream) :
  TObject(stream)
{
  //
  // copy ctor
  //
  AliFatal("copy constructor not implemented");
}

//______________________________________________________________________
AliMUONRawStreamTrigger& AliMUONRawStreamTrigger::operator = (const AliMUONRawStreamTrigger& 
					      /* stream */)
{ 
  // 
  // assignment operator
  //
  AliFatal("assignment operator not implemented");
  return *this;
}

//___________________________________
AliMUONRawStreamTrigger::~AliMUONRawStreamTrigger()
{
  //
  // clean up
  //
  delete fDDLTrigger;
  delete fLocalStruct;
  delete fRegHeader;
}

//_____________________________________________________________
Bool_t AliMUONRawStreamTrigger::Next()
{
// read the next raw digit (buspatch structure)
// returns kFALSE if there is no digit left

//   if (fNextDDL){
//     if(!NextDDL()) return kFALSE;
//   }
//   Int_t nEntries = fDDLTrigger->GetBusPatchEntries();

//   if (fSubEntries < nEntries) {
//     fLocalStruct =  (AliMUONLocalStruct*)fDDLTrigger->GetBusPatchEntry(fSubEntries);
//     fSubEntries++;
//     fNextDDL = kFALSE;
//     return kTRUE;
//   } else {
//     fDDLTrigger->GetBusPatchArray()->Delete();
//     fSubEntries = 0;
//     fNextDDL = kTRUE;
//     return Next(); 
//   }

  return kFALSE;
}

//______________________________________________________
Bool_t AliMUONRawStreamTrigger::NextDDL()
{
  // reading tracker DDL
  // store buspatch info into Array
  // store only non-empty structures (buspatch info with datalength !=0)

 // reading DDL for trigger

  AliMUONDarcHeader* darcHeader = fDDLTrigger->GetDarcHeader();

  Int_t kDarcHeaderSize = darcHeader->GetHeaderLength(); 
  Int_t kRegHeaderSize  = fRegHeader->GetHeaderLength() ;
  Bool_t scalerEvent    = kFALSE;

  // reset TClones
  darcHeader->GetRegHeaderArray()->Delete();
  //darcHeader->GetRegHeaderArray()->Clear("C");

  // loop over the two ddl's
  if (fDDL >= fMaxDDL) {
    fDDL = 0;
    return kFALSE;
  }

  fRawReader->Reset();
  fRawReader->Select(0XA,fDDL,fDDL);  //Select the DDL file to be read  

  fRawReader->ReadHeader();

  Int_t totalDataWord = fRawReader->GetDataSize(); // in bytes
  UInt_t *buffer = new UInt_t[totalDataWord/4];

  fRawReader->ReadNext((UChar_t*)buffer, totalDataWord); 
  
  Int_t index = 0;
  darcHeader->SetWord(buffer[index++]);

  if(darcHeader->GetEventType() == 2) {
    scalerEvent = kTRUE;
  } else
    scalerEvent = kFALSE;

  if(scalerEvent) {
    // 6 DARC scaler words
    memcpy(darcHeader->GetDarcScalers(), &buffer[index], darcHeader->GetDarcScalerLength()*4);
    index += darcHeader->GetDarcScalerLength();
  }

  if (buffer[index++] != darcHeader->GetEndOfDarc())
    AliWarning(Form("Wrong end of Darc word %d instead of %d\n",buffer[index-1], darcHeader->GetEndOfDarc())); 

  // 4 words of global board input + Global board output
  memcpy(darcHeader->GetGlobalInput(), &buffer[index], (kDarcHeaderSize-1)*4); 
  index += kDarcHeaderSize- 1; // kind tricky cos scaler info in-between Darc header

  if(scalerEvent) {
    // 10 Global scaler words
    memcpy(darcHeader->GetGlobalScalers(), &buffer[index], darcHeader->GetGlobalScalerLength()*4);
    index += darcHeader->GetGlobalScalerLength();
  }

  if (buffer[index++] != darcHeader->GetEndOfGlobal())
    AliWarning(Form("Wrong end of Global word %d instead of %d\n",buffer[index-1], darcHeader->GetEndOfGlobal()));
    
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
      AliWarning(Form("Wrong end of Reg word %d instead of %d\n",buffer[index-1], fRegHeader->GetEndOfReg()));

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
	AliWarning(Form("Wrong end of local word %d instead of %d\n",buffer[index-1], fLocalStruct->GetEndOfLocal()));
	  
      fDDLTrigger->AddLocStruct(*fLocalStruct, iReg);

    } // local card loop
	
  } // regional card loop
      
  delete [] buffer;


  fDDL++;

  return kTRUE;
}

//______________________________________________________
void AliMUONRawStreamTrigger::ResetDDL()
{
  // reseting TClonesArray
  // after each DDL
  //
  AliMUONDarcHeader* darcHeader = fDDLTrigger->GetDarcHeader();
  darcHeader->GetRegHeaderArray()->Clear("C");
}

//______________________________________________________
void AliMUONRawStreamTrigger::SetMaxDDL(Int_t ddl) 
{
  // set DDL number
  if (ddl > 2) ddl = 2;
  fMaxDDL = ddl;
}

//______________________________________________________
void AliMUONRawStreamTrigger::SetMaxReg(Int_t reg) 
{
  // set regional card number
  if (reg > 8) reg = 8;
  fMaxReg = reg;
}

//______________________________________________________
void AliMUONRawStreamTrigger::SetMaxLoc(Int_t loc) 
{
  // set local card number
  if (loc > 16) loc = 16;
  fMaxLoc = loc;
}
