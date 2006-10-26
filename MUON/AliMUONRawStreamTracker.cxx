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
/// it returns kFALSE (under develpment)
/// It can loop also over DDL and store the decoded rawdata in TClonesArray
/// in Payload class.
/// 
/// First version implement for Tracker
///
///////////////////////////////////////////////////////////////////////////////

#include "AliMUONRawStreamTracker.h"

#include "AliRawReader.h"
#include "AliRawDataHeader.h"
#include "AliDAQ.h"
#include "AliLog.h"

#include "AliMpBusPatch.h"

/// \cond CLASSIMP
ClassImp(AliMUONRawStreamTracker)
/// \endcond

AliMUONRawStreamTracker::AliMUONRawStreamTracker()
  : TObject(),
    fRawReader(0x0),
    fDDL(0),
    fBusPatchId(0),
    fDspId(0),
    fBlkId(0),
    fNextDDL(kTRUE),
    fMaxDDL(20),
    fPayload(new AliMUONPayloadTracker())
{
  ///
  /// create an object to read MUON raw digits
  /// Default ctor for monitoring purposes
  ///


}

//_________________________________________________________________
AliMUONRawStreamTracker::AliMUONRawStreamTracker(AliRawReader* rawReader)
  : TObject(),
    fRawReader(rawReader),
    fDDL(0),
    fBusPatchId(0),
    fDspId(0),
    fBlkId(0),
    fNextDDL(kTRUE),
    fMaxDDL(20),
    fPayload(new AliMUONPayloadTracker())

{
  ///
  /// ctor with AliRawReader as argument
  /// for reconstruction purpose
  ///


}

//___________________________________
AliMUONRawStreamTracker::~AliMUONRawStreamTracker()
{
  ///
  /// clean up
  ///
  delete fPayload;

}

//_____________________________________________________________
Bool_t AliMUONRawStreamTracker::Next()
{
  ///
  /// read the next raw digit (buspatch structure)
  /// returns kFALSE if there is no digit left
  /// (under development)

//      AliMUONDDLTracker*       ddlTracker = 0x0;
//      AliMUONBlockHeader*      blkHeader  = 0x0;
//      AliMUONDspHeader*        dspHeader  = 0x0;
//      Int_t nBusPatch;
//      Int_t nDsp;
//      Int_t nBlock;

//  next:  
//      if (fNextDDL){
//        printf("iDDL %d\n", fDDL+1);
//        fBlkId = 0;
//        fDspId = 0;
//        fBusPatchId = 0;
//        if(!NextDDL()) 
// 	 return kFALSE;
//      }
//      fNextDDL = kFALSE;

//      ddlTracker = GetDDLTracker();

//      nBlock = ddlTracker->GetBlkHeaderEntries();
//      if (fBlkId <  nBlock) {

//        blkHeader = ddlTracker->GetBlkHeaderEntry(fBlkId);
//        nDsp      = blkHeader->GetDspHeaderEntries();

//        if( fDspId < nDsp) {
// 	 dspHeader = blkHeader->GetDspHeaderEntry(fDspId);
// 	 nBusPatch = dspHeader->GetBusPatchEntries();

// 	 if (fBusPatchId < nBusPatch) {
// 	   fBusStructPtr = dspHeader->GetBusPatchEntry(fBusPatchId++);
// 	   return kTRUE;

// 	 } else {// iBusPatch
// 	   fDspId++;
// 	   fBusPatchId = 0;
// 	   goto next;
// 	   //	Next();
// 	 }

//        } else {// iDsp
// 	 fBlkId++;
// 	 fDspId = 0;
// 	 fBusPatchId = 0;
// 	 goto next;
// 	 //      Next();
//        }

//      } else {// iBlock
//        fBlkId = 0;
//        fDspId = 0;
//        fBusPatchId = 0;
//        fNextDDL = kTRUE;
//        //return kTRUE;
//        goto next; 
//      }

     return kFALSE;
}

//______________________________________________________
Bool_t AliMUONRawStreamTracker::NextDDL()
{
  /// reading tracker DDL

  fPayload->ResetDDL();


  if (fDDL >= 20) {
    fDDL = 0;
    return kFALSE;
  }
  AliDebug(3, Form("DDL Number %d\n", fDDL ));

  fRawReader->Reset();
  fRawReader->Select("MUONTRK", fDDL, fDDL);  //Select the DDL file to be read  

  fRawReader->ReadHeader();

  Int_t totalDataWord  = fRawReader->GetDataSize(); // in bytes

  UInt_t *buffer = new UInt_t[totalDataWord/4];

  fRawReader->ReadNext((UChar_t*)buffer, totalDataWord); 

  fPayload->Decode(buffer, totalDataWord/4);

  delete[] buffer;

  fDDL++;

  return kTRUE;
}

//______________________________________________________
void AliMUONRawStreamTracker::SetMaxDDL(Int_t ddl) 
{
  /// set DDL number
  if (ddl > 20) ddl = 20;
  fMaxDDL = ddl;

}

//______________________________________________________
void AliMUONRawStreamTracker::SetMaxBlock(Int_t blk) 
{
  /// set regional card number
  fPayload->SetMaxBlock(blk);
}
