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

/* $Id $ */

//-----------------------------------------------------------------------------
/// \class AliMUONRawStreamTrigger
/// This class provides access to MUON digits in raw data.
///
/// It loops over all MUON digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE(under develpment).
/// It can loop also over DDL and store the decoded rawdata in TClonesArrays
/// in payload class.
/// 
/// First version implement for Trigger
/// \author Christian Finck
//-----------------------------------------------------------------------------

#include "AliMUONRawStreamTrigger.h"

#include "AliRawReader.h"
#include "AliRawDataHeader.h"
#include "AliDAQ.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONRawStreamTrigger)
/// \endcond

const Int_t AliMUONRawStreamTrigger::fgkMaxDDL = 2;

//___________________________________________
AliMUONRawStreamTrigger::AliMUONRawStreamTrigger()
  : TObject(),
    fRawReader(0x0),
    fPayload(new AliMUONPayloadTrigger()),
    fDDL(0),
    fSubEntries(0),
    fNextDDL(kTRUE),
    fEnableErrorLogger(kFALSE)
{
  ///
  /// create an object to read MUON raw digits
  /// Default ctor for monitoring purposes
  ///


}

//_________________________________________________________________
AliMUONRawStreamTrigger::AliMUONRawStreamTrigger(AliRawReader* rawReader)
  : TObject(),
    fRawReader(rawReader),
    fPayload(new AliMUONPayloadTrigger()),
    fDDL(0),
    fSubEntries(0),
    fNextDDL(kTRUE),
    fEnableErrorLogger(kFALSE)
{
  ///
  /// ctor with AliRawReader as argument
  /// for reconstruction purpose
  ///

}

//___________________________________
AliMUONRawStreamTrigger::~AliMUONRawStreamTrigger()
{
  ///
  /// clean up
  ///
  delete fPayload;
}

//_____________________________________________________________
Bool_t AliMUONRawStreamTrigger::Next()
{
/// read the next raw digit (buspatch structure)
/// returns kFALSE if there is no digit left

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
  /// reading tracker DDL
  /// store buspatch info into Array
  /// store only non-empty structures (buspatch info with datalength !=0)

  // reset TClones
  fPayload->ResetDDL();


  // loop over the two ddl's

  while ( fDDL < fgkMaxDDL ) {
    fRawReader->Reset();
    fRawReader->Select("MUONTRG", fDDL, fDDL);  //Select the DDL file to be read  
    if (fRawReader->ReadHeader()) break;
    AliDebug(3,Form("Skipping DDL %d which does not seem to be there",fDDL));
    ++fDDL;
  }

  if (fDDL >= fgkMaxDDL) {
    fDDL = 0;
    return kFALSE;
  }

  AliDebug(3, Form("DDL Number %d\n", fDDL ));

  Int_t totalDataWord = fRawReader->GetDataSize(); // in bytes

  UInt_t *buffer = new UInt_t[totalDataWord/4];

  // check not necessary yet, but for future developments
  if (!fRawReader->ReadNext((UChar_t*)buffer, totalDataWord)) return kFALSE; 
  
  fPayload->Decode(buffer);
  if (fEnableErrorLogger) AddErrorMessage();

  fDDL++;

  delete [] buffer;


  return kTRUE;
}

//______________________________________________________
void AliMUONRawStreamTrigger::SetMaxReg(Int_t reg) 
{
  /// set regional card number
  fPayload->SetMaxReg(reg);
}

//______________________________________________________
void AliMUONRawStreamTrigger::SetMaxLoc(Int_t loc) 
{
  /// set local card number
  fPayload->SetMaxLoc(loc);
}

//______________________________________________________
void AliMUONRawStreamTrigger::AddErrorMessage()
{
/// add message into logger of AliRawReader per event

    for (Int_t i = 0; i < fPayload->GetDarcEoWErrors(); ++i)
	fRawReader->AddMajorErrorLog(kDarcEoWErr, "Wrong end of Darc word structure");

   for (Int_t i = 0; i < fPayload->GetGlobalEoWErrors(); ++i)
	fRawReader->AddMajorErrorLog(kGlobalEoWErr, "Wrong end of Global word structure");

   for (Int_t i = 0; i < fPayload->GetRegEoWErrors(); ++i)
	fRawReader->AddMajorErrorLog(kRegEoWErr, "Wrong end of Regional word structure");

   for (Int_t i = 0; i < fPayload->GetLocalEoWErrors(); ++i)
	fRawReader->AddMajorErrorLog(kLocalEoWErr, "Wrong end of Local word structure");

}
