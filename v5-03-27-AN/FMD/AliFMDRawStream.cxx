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

/* $Id$ */

//____________________________________________________________________
//                                                                          
// Buffer to read RAW ALTRO FMD format from a AliRawReader 
// 
// This class derives from AliAltroBuffer, but overloads the memer
// function Next to do some extra processing.  In particular, it tries
// to autodetect the sample rate.  If zero-suppression was used when
// writing the raw data, then the automatic discovery will not work,
// and the sample rate should be set explicitly. 
//
#include "AliFMDRawStream.h"		// ALIFMDRAWSTREAM_H
// #include <AliRawReader.h>		// ALIRAWREADER_H
#include "AliFMDParameters.h"
// #include <AliLog.h>
#include "AliFMDDebug.h" // Better debug macros
// #include <iomanip>
// #include <iostream>
#include "AliRawReader.h"
#include <climits>

//____________________________________________________________________
ClassImp(AliFMDRawStream)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDRawStream::AliFMDRawStream(AliRawReader* reader) 
  : AliAltroRawStream(reader)
{
  // CTOR 
  reader->Reset();
  // Select FMD DDL's 
  SelectRawData("FMD");
}

//_____________________________________________________________________________
Bool_t 
AliFMDRawStream::ReadChannel(UInt_t& ddl, UInt_t& addr, 
			     UInt_t& len, volatile UShort_t* data)
{
  // Read one channel and return.   Returns 0 when there's no more
  // data. 
  Int_t        l         = 0;
  static Int_t last      = 0xFFFF; // 0xFFFF means signal is used
  Bool_t       next      = kTRUE;
  do {
    Int_t signal = last;
    if (last > 0x3FF) {
      AliFMDDebug(30, ("Last is 0x%x, so reading a new word", last));
      next   = Next();
      if(!next){
	AliFMDDebug(15, ("Read word # %d (!next)", l));
	addr = GetPrevHWAddress();
	ddl  = (GetPrevDDLNumber() < 0 ? UINT_MAX: UInt_t(GetPrevDDLNumber()));
	len  = l+1; // Need to add one - l points to last valid index
	last = signal;
	break;
      }
      signal = GetSignal();
      if (GetHWAddress() != GetPrevHWAddress() && GetPrevHWAddress() >= 0) {
	AliFMDDebug(15, ("New hardware address, was 0x%x, now 0x%x", 
			  GetPrevHWAddress(), GetHWAddress()));
	addr = GetPrevHWAddress();
	ddl  = (GetPrevDDLNumber() < 0 ? UINT_MAX : UInt_t(GetPrevDDLNumber()));
	len  = l+1; // Need to add one - l points to last valid index
	last = signal;
	break;
      }
    }
    // Sanity check - if the total bunch length is less than 1, then
    // read until we get the next bunch. 
    Int_t b  = GetTimeLength();
    if (b < 1) { 
      AliWarning(Form("Bunch length %0d is less than 0 for "
		      "DDL %4d address 0x%03x", 
		      b, ddl, addr));
      last = 0xFFFF;
      continue;
    }

    // Sanity check - if the current time is less than 0, then read
    // until we get a new bunch. 
    Int_t t  = GetTime();
    if (t < 0) {
      AliWarning(Form("Time %0d is less than 0 for DDL %4d address 0x%03x", 
		      t, ddl, addr));
      last = 0xFFFF;
      continue;
    }
    l        = TMath::Max(l, t);
    data[t]  = signal;
    last     = 0xFFFF;
#if 0
    AliFMDDebug(signal > 512 ? 1 : 0, ("Signal @ %d (%d) is %d", 
				       time, t, data[t]));
#endif
  } while (next);
  return next;
}


//_____________________________________________________________________________
//
// EOF
//
