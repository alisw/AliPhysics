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
/** @file    AliFMDRawReader.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:45:23 2006
    @brief   Class to read raw data 
    @ingroup FMD_rec
*/
//____________________________________________________________________
//
// Class to read ADC values from a AliRawReader object. 
//
// This class uses the AliFMDRawStreamer class to read the ALTRO
// formatted data. 
// 
//          +-------+
//          | TTask |
//          +-------+
//              ^
//              |
//      +-----------------+  <<references>>  +--------------+
//      | AliFMDRawReader |<>----------------| AliRawReader |
//      +-----------------+                  +--------------+
//              |                                  ^
//              | <<uses>>                         |
//              V                                  |
//      +-----------------+      <<uses>>          |
//      | AliFMDRawStream |------------------------+
//      +-----------------+
//              |
//              V
//      +----------------+
//      | AliAltroStream |
//      +----------------+
//
// #include <AliLog.h>		// ALILOG_H
#include "AliFMDDebug.h" // Better debug macros
#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDRawStream.h"	// ALIFMDRAWSTREAM_H 
// #include "AliRawReader.h"	// ALIRAWREADER_H 
#include "AliFMDRawReader.h"	// ALIFMDRAWREADER_H 
// #include "AliFMDAltroIO.h"	// ALIFMDALTROIO_H 
// #include <TArrayI.h>		// ROOT_TArrayI
#include <TTree.h>		// ROOT_TTree
#include <TClonesArray.h>	// ROOT_TClonesArray
// #include <iostream>
// #include <iomanip>

//____________________________________________________________________
ClassImp(AliFMDRawReader)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDRawReader::AliFMDRawReader(AliRawReader* reader, TTree* tree) 
  : TTask("FMDRawReader", "Reader of Raw ADC values from the FMD"),
    fTree(tree),
    fReader(reader), 
    fSampleRate(1)
{
  // Default CTOR
}

//____________________________________________________________________
void
AliFMDRawReader::Exec(Option_t*) 
{
  // Read the data 
  TClonesArray* array = new TClonesArray("AliFMDDigit");
  if (!fTree) {
    AliError("No tree");
    return;
  }
  fTree->Branch("FMD", &array);
  ReadAdcs(array);
  Int_t nWrite = fTree->Fill();
  AliFMDDebug(1, ("Got a grand total of %d digits, wrote %d bytes to tree", 
		   array->GetEntries(), nWrite));
}


#if 1
//____________________________________________________________________
Bool_t
AliFMDRawReader::ReadAdcs(TClonesArray* array) 
{
  // Read raw data into the digits array, using AliFMDAltroReader. 
  if (!array) {
    AliError("No TClonesArray passed");
    return kFALSE;
  }
  //  if (!fReader->ReadHeader()) {
  //    AliError("Couldn't read header");
  //    return kFALSE;
  //  }
  // Get sample rate 
  AliFMDParameters* pars = AliFMDParameters::Instance();
  AliFMDRawStream input(fReader);

  UShort_t stripMin = 0;
  UShort_t stripMax = 127;
  UShort_t preSamp  = 0;
  
  UInt_t ddl    = 0;
  UInt_t rate   = 0;
  UInt_t last   = 0;
  UInt_t hwaddr = 0;
  // Data array is approx twice the size needed. 
  UShort_t data[2048];

  Bool_t isGood = kTRUE;
  while (isGood) {
    isGood = input.ReadChannel(ddl, hwaddr, last, data);

    AliFMDDebug(5, ("Read channel 0x%x of size %d", hwaddr, last));
    UShort_t det, sec, str;
    Char_t   ring;
    if (!pars->Hardware2Detector(ddl, hwaddr, det, ring, sec, str)) {
      AliError(Form("Failed to get detector id from DDL %d "
		    "and hardware address 0x%x", ddl, hwaddr));
      continue;
    }
    rate     = pars->GetSampleRate(det, ring, sec, str);
    stripMin = pars->GetMinStrip(det, ring, sec, str);
    stripMax = pars->GetMaxStrip(det, ring, sec, str);
    AliFMDDebug(5, ("DDL 0x%04x, address 0x%03x maps to FMD%d%c[%2d,%3d]", 
		       ddl, hwaddr, det, ring, sec, str));
    
    // Loop over the `timebins', and make the digits
    for (size_t i = 0; i < last; i++) {
      if (i < preSamp) continue;
      Int_t    n      = array->GetEntriesFast();
      UShort_t curStr = str + stripMin + i / rate;
      if ((curStr-str) > stripMax) {
	AliError(Form("Current strip is %d but DB says max is %d", 
		      curStr, stripMax));
      }
      AliFMDDebug(5, ("making digit for FMD%d%c[%2d,%3d] from sample %4d", 
		       det, ring, sec, curStr, i));
      new ((*array)[n]) AliFMDDigit(det, ring, sec, curStr, data[i], 
				    (rate >= 2 ? data[i+1] : 0),
				    (rate >= 3 ? data[i+2] : 0));
      if (rate >= 2) i++;
      if (rate >= 3) i++;
    }
  }
  return kTRUE;
}
#else
//____________________________________________________________________
Bool_t
AliFMDRawReader::ReadAdcs(TClonesArray* array) 
{
  // Read raw data into the digits array, using AliFMDAltroReader. 
  if (!array) {
    AliError("No TClonesArray passed");
    return kFALSE;
  }
  //  if (!fReader->ReadHeader()) {
  //    AliError("Couldn't read header");
  //    return kFALSE;
  //  }
  // Get sample rate 
  AliFMDParameters* pars = AliFMDParameters::Instance();

  // Select FMD DDL's 
  fReader->Select("FMD");

  UShort_t stripMin = 0;
  UShort_t stripMax = 127;
  UShort_t preSamp  = 0;
  
  do {
    UChar_t* cdata;
    if (!fReader->ReadNextData(cdata)) break;
    size_t   nchar = fReader->GetDataSize();
    UShort_t ddl   = fReader->GetDDLID();
    UShort_t rate  = 0;
    AliFMDDebug(1, ("Reading %d bytes (%d 10bit words) from %d", 
		     nchar, nchar * 8 / 10, ddl));
    // Make a stream to read from 
    std::string str((char*)(cdata), nchar);
    std::istringstream s(str);
    // Prep the reader class.
    AliFMDAltroReader r(s);
    // Data array is approx twice the size needed. 
    UShort_t data[2048], hwaddr, last;
    while (r.ReadChannel(hwaddr, last, data) > 0) {
      AliFMDDebug(5, ("Read channel 0x%x of size %d", hwaddr, last));
      UShort_t det, sec, str;
      Char_t   ring;
      if (!pars->Hardware2Detector(ddl, hwaddr, det, ring, sec, str)) {
	AliError(Form("Failed to detector id from DDL %d "
		      "and hardware address 0x%x", ddl, hwaddr));
	continue;
      }
      rate     = pars->GetSampleRate(det, ring, sec, str);
      stripMin = pars->GetMinStrip(det, ring, sec, str);
      stripMax = pars->GetMaxStrip(det, ring, sec, str);
      AliFMDDebug(5, ("DDL 0x%04x, address 0x%03x maps to FMD%d%c[%2d,%3d]", 
		       ddl, hwaddr, det, ring, sec, str));

      // Loop over the `timebins', and make the digits
      for (size_t i = 0; i < last; i++) {
	if (i < preSamp) continue;
	Int_t    n      = array->GetEntries();
	UShort_t curStr = str + stripMin + i / rate;
	if ((curStr-str) > stripMax) {
	  AliError(Form("Current strip is %d but DB says max is %d", 
			curStr, stripMax));
	}
	AliFMDDebug(5, ("making digit for FMD%d%c[%2d,%3d] from sample %4d", 
			 det, ring, sec, curStr, i));
	new ((*array)[n]) AliFMDDigit(det, ring, sec, curStr, data[i], 
				      (rate >= 2 ? data[i+1] : 0),
				      (rate >= 3 ? data[i+2] : 0));
	if (rate >= 2) i++;
	if (rate >= 3) i++;
	}
	if (r.IsBof()) break;
    }
  } while (true);
  return kTRUE;
}

  

// This is the old method, for comparison.   It's really ugly, and far
// too convoluted. 
//____________________________________________________________________
void
AliFMDRawReader::Exec(Option_t*) 
{
  // Read raw data into the digits array
  //  if (!fReader->ReadHeader()) {
  //    Error("ReadAdcs", "Couldn't read header");
  //    return;
  //  }

  Int_t n = 0;
  TClonesArray* array = new TClonesArray("AliFMDDigit");
  fTree->Branch("FMD", &array);

  // Get sample rate 
  AliFMDParameters* pars = AliFMDParameters::Instance();
  fSampleRate = pars->GetSampleRate(0);

  // Use AliAltroRawStream to read the ALTRO format.  No need to
  // reinvent the wheel :-) 
  AliFMDRawStream input(fReader, fSampleRate);
  // Select FMD DDL's 
  fReader->Select("FMD");
  
  Int_t    oldDDL      = -1;
  Int_t    count       = 0;
  UShort_t detector    = 1; // Must be one here
  UShort_t oldDetector = 0;
  Bool_t   next        = kTRUE;

  // local Cache 
  TArrayI counts(10);
  counts.Reset(-1);
  
  // Loop over data in file 
  while (next) {
    next = input.Next();

    count++; 
    Int_t ddl = fReader->GetDDLID();
    AliFMDDebug(10, ("Current DDL is %d", ddl));
    if (ddl != oldDDL || input.IsNewStrip() || !next) {
      // Make a new digit, if we have some data (oldDetector == 0,
      // means that we haven't really read anything yet - that is,
      // it's the first time we get here). 
      if (oldDetector > 0) {
	// Got a new strip. 
	AliFMDDebug(10, ("Add a new strip: FMD%d%c[%2d,%3d] "
			  "(current: FMD%d%c[%2d,%3d])", 
			  oldDetector, input.PrevRing(), 
			  input.PrevSector() , input.PrevStrip(),
			  detector , input.Ring(), input.Sector(), 
			  input.Strip()));
	new ((*array)[n]) AliFMDDigit(oldDetector, 
				      input.PrevRing(), 
				      input.PrevSector(), 
				      input.PrevStrip(), 
				      counts[0], counts[1], counts[2]);
	n++;
#if 0
	AliFMDDigit* digit = 
	  static_cast<AliFMDDigit*>(fFMD->Digits()->
				    UncheckedAt(fFMD->GetNdigits()-1));
#endif 
      }
	
      if (!next) { 
	AliFMDDebug(10, ("Read %d channels for FMD%d", 
			  count + 1, detector));
	break;
      }
    
    
      // If we got a new DDL, it means we have a new detector. 
      if (ddl != oldDDL) {
	if (detector != 0) 
	  AliFMDDebug(10, ("Read %d channels for FMD%d", count + 1, detector));
	// Reset counts, and update the DDL cache 
	count       = 0;
	oldDDL      = ddl;
	// Check that we're processing a FMD detector 
	Int_t detId = fReader->GetDetectorID();
	if (detId != (AliDAQ::DetectorID("FMD"))) {
	  AliError(Form("Detector ID %d != %d",
			detId, (AliDAQ::DetectorID("FMD"))));
	  break;
	}
	// Figure out what detector we're deling with 
	oldDetector = detector;
	switch (ddl) {
	case 0: detector = 1; break;
	case 1: detector = 2; break;
	case 2: detector = 3; break;
	default:
	  AliError(Form("Unknown DDL 0x%x for FMD", ddl));
	  return;
	}
	AliFMDDebug(10, ("Reading ADCs for 0x%x  - That is FMD%d",
			  fReader->GetEquipmentId(), detector));
      }
      counts.Reset(-1);
    }
    
    counts[input.Sample()] = input.Count();
    
    AliFMDDebug(10, ("ADC of FMD%d%c[%2d,%3d] += %d",
		      detector, input.Ring(), input.Sector(), 
		      input.Strip(), input.Count()));
    oldDetector = detector;
  }
  fTree->Fill();
  return;

}
#endif

//____________________________________________________________________
// 
// EOF
//
