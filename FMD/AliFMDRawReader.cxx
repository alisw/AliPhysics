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
#include "AliFMDSDigit.h"	// ALIFMDSDIGIT_H
#include "AliFMDRawStream.h"	// ALIFMDRAWSTREAM_H 
#include "AliRawReader.h"	// ALIRAWREADER_H 
#include "AliFMDRawReader.h"	// ALIFMDRAWREADER_H 
#include "AliFMDDebug.h"
#include "AliFMDCalibSampleRate.h"
#include "AliFMDCalibStripRange.h"
#include "AliFMDAltroMapping.h"
// #include "AliFMDAltroIO.h"	// ALIFMDALTROIO_H 
#include <TArrayS.h>		// ROOT_TArrayS
#include <TTree.h>		// ROOT_TTree
#include <TClonesArray.h>	// ROOT_TClonesArray
#include <TString.h>
#include <iostream>
#include <climits>
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
    // fSampleRate(1),
    fData(0),
    fNbytes(0), 
    fSeen(0)
{
  // Default CTOR
  for (Int_t i = 0; i < 3; i++) { 
    fSampleRate[i]   = 0;
    fZeroSuppress[i] = kFALSE;
    fNoiseFactor[i]  = 1;
  }
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
		   array->GetEntriesFast(), nWrite));
}

//____________________________________________________________________
Bool_t
AliFMDRawReader::NextSample(UShort_t& det, Char_t&   rng, UShort_t& sec, 
			    UShort_t& str, UShort_t& sam, UShort_t& rat, 
			    Short_t&  adc, Bool_t&   zs,  UShort_t& fac)
{
  // Scan current event for next signal.   It returns kFALSE when
  // there's no more data in the event. 
  static AliAltroRawStream   stream(fReader); //    = 0;
  static AliFMDParameters*   pars     = 0;
  static AliFMDAltroMapping* map      = 0;
  static Int_t               ddl      = -1;
  // static UInt_t              rate     = 0;
  static UShort_t            tdet     = 0;
  static Char_t              trng     = '\0';
  static UShort_t            tsec     = 0;
  static Short_t             tstr     = 0;   
  static Short_t             bstr     = -1;
  static Short_t             tsam     = -1;   
  static UInt_t              trate    = 0;
  static Int_t               hwaddr   = -1;
  static UShort_t            stripMin = 0;
  static UShort_t            stripMax = 0; // 127;
  static UShort_t            preSamp  = 0; // 14+5;
  if (stream.GetDDLNumber() < 0) { 
    fReader->Select("FMD");

    // Reset "seen" array
    // const UShort_t kUShortMax = (1 << 16) - 1;
    // fSeen.Reset(kUShortMax);

    pars   = AliFMDParameters::Instance();
    map    = pars->GetAltroMap();
    // stream = new AliAltroRawStream(fReader);

    AliFMDDebug(5, ("Setting %d word headers", 
		    pars->HasCompleteHeader() ? 8 : 7));
    stream.SetShortDataHeader(!pars->HasCompleteHeader());

    // Reset variables
    ddl    = -1;  
    trate  = 0;   
    tdet   = 0;   
    trng   = '\0';
    tsec   = 0;   
    tstr   = 0;  
    tsam   = -1;
    hwaddr = -1;
  }
  do { 
    Bool_t next = stream.Next();
    if (!next) { 
      // if (stream) delete stream;
      // stream = 0;
      return kFALSE;
    }
    Int_t thisDDL = stream.GetDDLNumber();
    AliFMDDebug(10, ("RCU @ DDL %d", thisDDL));
    if (thisDDL != ddl) { 
      ddl   = thisDDL;
      fZeroSuppress[ddl] = zs    = stream.GetZeroSupp();
      fNoiseFactor[ddl]  = fac   = stream.GetNPostsamples();
      fSampleRate[ddl]   = trate = 0; // stream.GetNPresamples();
      tdet  = map->DDL2Detector(ddl);
      AliFMDDebug(10, ("RCU @ DDL %d zero suppression: %s",ddl, zs?"yes":"no"));
      AliFMDDebug(10, ("RCU @ DDL %d noise factor: %d", ddl,fac));
      AliFMDDebug(10, ("RCU @ DDL %d sample rate: %d", ddl, trate));
    }
    Int_t thisAddr = stream.GetHWAddress();
    AliFMDDebug(10, ("RCU @ DDL %d, Address 0x%03x", ddl, thisAddr));    
    if (thisAddr != hwaddr) { 
      hwaddr = thisAddr;
      UShort_t  board, chip, channel;
      map->ChannelAddress(hwaddr, board, chip, channel);
      map->Channel2StripBase(board, chip, channel, trng, tsec, bstr);
      AliFMDDebug(10, ("0x%04x/0x%03x maps to FMD%d%c[%2d]-%3d", 
		      ddl, hwaddr, tdet, trng, tsec, bstr));
      
      stripMin = pars->GetMinStrip(tdet, trng, tsec, bstr);
      stripMax = pars->GetMaxStrip(tdet, trng, tsec, bstr);
      preSamp  = pars->GetPreSamples(tdet, trng, tsec, bstr);
      if (trate == 0) 
	fSampleRate[ddl] = trate = pars->GetSampleRate(tdet, trng, tsec, bstr);
      AliFMDDebug(10, ("RCU @ DDL %d, Address 0x%03x sample rate: %d", 
		      ddl, hwaddr, trate));
      
      Int_t nChAddrMismatch = stream.GetNChAddrMismatch();
      Int_t nChLenMismatch  = stream.GetNChLengthMismatch();
      if (nChAddrMismatch != 0) 
	AliWarning(Form("Got %d channels with address mis-matches for 0x%03x",
			nChAddrMismatch, hwaddr));
      if (nChLenMismatch != 0) 
	AliWarning(Form("Got %d channels with length mis-matches for 0x%03x",
			nChLenMismatch, hwaddr));
    }
    // Get the signal
    adc = stream.GetSignal();

    // Sanity check - if the total bunch length is less than 1, then
    // read until we get the next bunch. 
    Int_t b  = stream.GetTimeLength();
    if (b < 1) { 
      AliWarning(Form("Bunch length %0d is less than 0 for "
		      "DDL %4d address 0x%03x", b, ddl, hwaddr));
      continue;
    }
    // Sanity check - if the current time is less than 0, then read
    // until we get a new bunch. 
    Int_t t  = stream.GetTime();
    if (t < 0) {
      AliWarning(Form("Time %0d is less than 0 for DDL %4d address 0x%03x", 
		      t, ddl, hwaddr));
      continue;
    }

    Short_t  strOff = 0;
    UShort_t samp   = 0;
    map->Timebin2Strip(tsec, t, preSamp, trate, strOff, samp);
    tstr = bstr + strOff;
    tsam = samp;
    AliFMDDebug(20, ("0x%04x/0x%03x/%04d maps to FMD%d%c[%2d,%3d]-%d", 
		     ddl, hwaddr, t, tdet, trng, tsec, tstr, samp));
    
    // Local strip number to channel
    Short_t l = (tstr > bstr ? tstr - bstr : bstr - tstr);
    AliFMDDebug(10, ("Checking if strip %d in range [%d,%d]", 
		     l, stripMin, stripMax));
    if (l < stripMin || l > stripMax) { 
      AliFMDDebug(5, ("VA channel %3d (t: %4d) of DDL %4d address 0x%03x "
		      "is out of range (%3d->%3d)", 
		      l, t, ddl, hwaddr, stripMin, stripMax));
      continue;
    }


    det = tdet;
    rng = trng;
    sec = tsec;
    str = tstr;
    sam = tsam;
    rat = trate;
    // adc = stream.GetSignal();
    
    break;
  } while (true);
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDRawReader::NextSignal(UShort_t& det, Char_t&   rng, 
			    UShort_t& sec, UShort_t& str, 
			    Short_t&  adc, Bool_t&   zs, 
			    UShort_t& fac)
{
  
  do { 
    UShort_t samp, rate;
    if (!NextSample(det, rng, sec, str, samp, rate, adc, zs, fac)) 
      return kFALSE;

    Bool_t take = kFALSE;
    switch (rate) { 
    case 1:                      take = kTRUE; break;
    case 2:  if (samp == 1)      take = kTRUE; break;
    case 3:  if (samp == 1)      take = kTRUE; break; 
    case 4:  if (samp == 2)      take = kTRUE; break;
    default: if (samp == rate-2) take = kTRUE; break;
    }
    if (!take) continue;
    break;
  } while (true);
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDRawReader::SelectSample(UShort_t samp, UShort_t rate) 
{
  Bool_t take = kFALSE;
  switch (rate) { 
  case 1:                      take = kTRUE; break;
  case 2:  if (samp == 1)      take = kTRUE; break;
  case 3:  if (samp == 1)      take = kTRUE; break; 
  case 4:  if (samp == 2)      take = kTRUE; break;
  default: if (samp == rate-2) take = kTRUE; break;
  }
  
  return take;
}
  
#if 1
//____________________________________________________________________
Bool_t
AliFMDRawReader::ReadAdcs(TClonesArray* array) 
{
  
  fReader->Select("FMD");
  // Read raw data into the digits array, using AliFMDAltroReader. 
  if (!array) {
    AliError("No TClonesArray passed");
    return kFALSE;
  }
  const UShort_t kUShortMax = (1 << 16) - 1;
  fSeen.Reset(kUShortMax);
  
  //  if (!fReader->ReadHeader()) {
  //    AliError("Couldn't read header");
  //    return kFALSE;
  //  }
  // Get sample rate 
  AliFMDParameters*     pars = AliFMDParameters::Instance();
  AliFMDAltroMapping*   map  = pars->GetAltroMap();
  AliFMDRawStream       input(fReader);
  AliFMDDebug(5, ("Setting 7 word headers"));
  input.SetShortDataHeader(!pars->HasCompleteHeader());

  UShort_t stripMin = 0;
  UShort_t stripMax = 0; // 127;
  UShort_t preSamp  = 0; // 14+5;

  Int_t  oldddl = -1;
  UInt_t ddl    = 0;
  // UInt_t rate   = 0;
  UInt_t last   = 0;
  UInt_t hwaddr = 0;
  // Data array is approx twice the size needed. 
  UShort_t data[2048];
  for (size_t i = 0; i < 2048; i++) data[i] = 0; // kUShortMax;

  Bool_t isGood = kTRUE;
  while (isGood) {
    isGood = input.ReadChannel(ddl, hwaddr, last, data);
    // if (!isGood) break;
    if (ddl >= UINT_MAX /* UInt_t(-1) */) { 
      AliFMDDebug(5, ("At end of event with %d digits", 
		      array->GetEntriesFast()));
      break;
    }
    if (UInt_t(oldddl) != ddl) { 
      fZeroSuppress[ddl] = input.GetZeroSupp();
      AliFMDDebug(20, ("RCU @ DDL %d zero suppression: %s", 
		       ddl, (fZeroSuppress[ddl] ? "yes" : "no")));

      // WARNING: We store the noise factor in the 2nd baseline
      // filters excluded post samples, since we'll never use that
      // mode. 
      fNoiseFactor[ddl]  = input.GetNPostsamples();
      AliFMDDebug(20, ("RCU @ DDL %d noise factor: %d", ddl,fNoiseFactor[ddl]));

      // WARNING: We store the noise factor in the 2nd baseline
      // filters excluded post samples, since we'll never use that
      // mode. 
      fSampleRate[ddl]     = input.GetNPretriggerSamples();
      AliFMDDebug(20, ("RCU @ DDL %d Sample rate: %d", ddl,fNoiseFactor[ddl]));

      Int_t nChAddrMismatch = input.GetNChAddrMismatch();
      Int_t nChLenMismatch  = input.GetNChLengthMismatch();
      if (nChAddrMismatch != 0) 
	AliWarning(Form("Got %d channels with address mis-matches for 0x%03x",
			nChAddrMismatch, hwaddr));
      if (nChLenMismatch != 0) 
	AliWarning(Form("Got %d channels with length mis-matches for 0x%03x",
			nChLenMismatch, hwaddr));
      oldddl = ddl;
    }
    // AliFMDDebug(5, ("Read channel 0x%x of size %d", hwaddr, last));

    UShort_t det, sec, samp, board, chip, channel;
    Short_t strbase;
    Char_t   ring;
    
    
    if (map->DDL2Detector(ddl) < 0) break;
    det = map->DDL2Detector(ddl);
    map->ChannelAddress(hwaddr, board, chip, channel);
    if (!map->Channel2StripBase(board, chip, channel, ring, sec, strbase)) {
      AliError(Form("Failed to get detector id from DDL %d, "
		    "hardware address 0x%03x", ddl, hwaddr));
      continue;
    }
    AliFMDDebug(5, ("Board: 0x%02x, Altro: 0x%x, Channel: 0x%x, Length: %4d", 
		    board, chip, channel, last));

    stripMin = pars->GetMinStrip(det, ring, sec, strbase);
    stripMax = pars->GetMaxStrip(det, ring, sec, strbase);
    preSamp  = pars->GetPreSamples(det, ring, sec, strbase);
    if (fSampleRate[ddl] == 0) 
      fSampleRate[ddl] = pars->GetSampleRate(det, ring, sec, strbase);
    
    // Loop over the `timebins', and make the digits
    for (size_t i = 0; i < last; i++) {
      // if (i < preSamp) continue;
      AliFMDDebug(15, ("0x%04x/0x%03x/%04d %4d", ddl, hwaddr, i, data[i]));

      Short_t  stroff = 0;
      map->Timebin2Strip(sec, i, preSamp, fSampleRate[ddl], stroff, samp);
      Short_t  str    = strbase + stroff;
      
      AliFMDDebug(10, ("0x%04x/0x%03x/%04d maps to FMD%d%c[%2d,%3d]-%d", 
		      ddl, hwaddr, i, det, ring, sec, str, samp));
      if (str < 0) { 
	AliFMDDebug(8, ("Got presamples at timebin %d", i));
	data[i] = 0; // Reset cache 
	continue;
      }
      
      Short_t lstrip = (i - preSamp) / fSampleRate[ddl] + stripMin;
      
      AliFMDDebug(15, ("Checking if strip %d (%d) in range [%d,%d]", 
		      lstrip, str, stripMin, stripMax));
      if (lstrip < stripMin || lstrip > stripMax) {
	AliFMDDebug(5, ("FMD%d%c[%02d,%03d]-%d out of range (%3d->%3d)", 
			det, ring, sec, samp, str, stripMin, stripMax));
	data[i] = 0; // Reset cache 
	continue;
      }
      // Possibly do pedestal subtraction of signal 
      Int_t counts = data[i];
	
      
      // Check the cache of indicies
      Int_t idx = fSeen(det, ring, sec, str);
      if (idx == kUShortMax) { 
	// We haven't seen this strip yet. 
	fSeen(det, ring, sec, str) = idx = array->GetEntriesFast();
	AliFMDDebug(7,("making digit for FMD%d%c[%2d,%3d]-%d from timebin %4d", 
		       det, ring, sec, str, samp, i));
	new ((*array)[idx]) AliFMDDigit(det, ring, sec, str);
      }
      AliFMDBaseDigit* digit = static_cast<AliFMDBaseDigit*>(array->At(idx));
      AliFMDDebug(10,
		  ("Setting FMD%d%c[%2d,%3d]-%d from timebin %4d=%4d (%4d)", 
		   det, ring, sec, str, samp, i, counts, data[i]));
      digit->SetCount(samp, counts);
      data[i] = 0; // Reset cache 
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
Bool_t AliFMDRawReader::ReadSODevent(AliFMDCalibSampleRate* sampleRate, 
				     AliFMDCalibStripRange* stripRange, 
				     TArrayS &pulseSize, 
				     TArrayS &pulseLength) 
{

  AliFMDDebug(0, ("Start of SOD/EOD"));
  
  UInt_t shift_clk[18];
  UInt_t sample_clk[18];
  UInt_t strip_low[18];
  UInt_t strip_high[18];
  UInt_t pulse_size[18];
  UInt_t pulse_length[18];  
  AliFMDParameters*   param = AliFMDParameters::Instance();
  AliFMDAltroMapping* map   = param->GetAltroMap();
  
  while(fReader->ReadNextData(fData)) {
    
    Int_t ddl   = fReader->GetDDLID();
    Int_t detID = fReader->GetDetectorID();
    AliFMDDebug(0, (" From reader: DDL number is %d , det ID is %d",ddl,detID));
    
    fNbytes = fReader->GetDataSize();
    
    ULong_t nWords      = GetNwords();
    UInt_t  trailerLast = Get32bitWord(nWords);
    AliFMDDebug(20, (" # Bytes: %d, # Words: %d, Last word: 0x%08x", 
		     fNbytes, nWords, trailerLast));
    if ((trailerLast & 0xFFFF0000) != 0xAAAA0000) {
      AliWarning(Form("Last word 0x%08x does not match RCU II trailer",
		      trailerLast));
      return kFALSE;
    }
    
    // 7 last bits is size of trailer
    ULong_t nTrailerWords = trailerLast & 0x7f; 
    ULong_t nPayloadWords = Get32bitWord(nWords - nTrailerWords+1);

    AliFMDDebug(20, (" # trailer words: %d, # payload words: %d", 
		     nTrailerWords, nPayloadWords));
    
    
    for (ULong_t i = 1; i <= nPayloadWords ; i++) {
      UInt_t payloadWord = Get32bitWord(i);
      
      // address is only 24 bit
      UInt_t address       = (0xffffff & payloadWord);
      UInt_t type          = ((address >> 21) & 0xf);
      UInt_t error         = ((address >> 20) & 0x1);
      UInt_t bcast         = ((address >> 18) & 0x1);
      UInt_t bc_not_altro  = ((address >> 17) & 0x1);
      UInt_t board         = ((address >> 12) & 0x1f);
      UInt_t instruction   = 0;
      UInt_t chip          = 0;
      UInt_t channel       = 0;
      if(bc_not_altro)
	instruction        = address & 0xfff;
      else {
	chip               = ((address >> 9) & 0x7);
	channel            = ((address >> 5) & 0x5);
	instruction        = (address & 0x1f);
      }
	
      Bool_t readDataWord = kFALSE;
      switch(type) {
      case 0x0: // Fec read
	readDataWord = kTRUE;  
      case 0x1: // Fec cmd
      case 0x2: // Fec write
	i++;  
	break;
      case 0x4: // Loop
      case 0x5: // Wait
	break;
      case 0x6: // End sequence
      case 0x7: // End Mem
	i = nPayloadWords + 1;
	break;
      default:    
	break;
      }
	
      //Don't read unless we have a FEC_RD
      if(!readDataWord)  continue;

      UInt_t dataWord      = Get32bitWord(i);
      UInt_t data          = (0xFFFFF & dataWord) ;
      //UInt_t data          = (0xFFFF & dataWord) ;
	
      if(error) {
	AliWarning(Form("error bit detected at Word 0x%06x; "
			"error % d, type %d, bc_not_altro %d, "
			"bcast %d, board 0x%02x, chip 0x%x, "
			"channel 0x%02x, instruction 0x%03x",
			address, error, type, bc_not_altro, 
			bcast,board,chip,channel,instruction));
	//process error
	continue;
      }
	
	
      switch(instruction) {
	  
      case 0x01: break;  // First ADC T           
      case 0x02: break; // I  3.3 V              
      case 0x03: break; // I  2.5 V altro digital
      case 0x04: break; // I  2.5 V altro analog 
      case 0x05: break; // I  2.5 V VA           
      case 0x06: break; // First ADC T           
      case 0x07: break; // I  3.3 V              
      case 0x08: break; // I  2.5 V altro digital
      case 0x09: break; // I  2.5 V altro analog 
      case 0x0A: break; // I  2.5 V VA           
      case 0x2D: break; // Second ADC T           
      case 0x2E: break; // I  1.5 V VA            
      case 0x2F: break; // I -2.0 V               
      case 0x30: break; // I -2.0 V VA            
      case 0x31: break; //    2.5 V Digital driver
      case 0x32: break; // Second ADC T           
      case 0x33: break; // I  1.5 V VA            
      case 0x34: break; // I -2.0 V               
      case 0x35: break; // I -2.0 V VA            
      case 0x36: break; //    2.5 V Digital driver
      case 0x37: break; // Third ADC T             
      case 0x38: break; // Temperature sens. 1     
      case 0x39: break; // Temperature sens. 2     
      case 0x3A: break; // U  2.5 altro digital (m)
      case 0x3B: break; // U  2.5 altro analog (m) 
      case 0x3C: break; // Third ADC T             
      case 0x3D: break; // Temperature sens. 1     
      case 0x3E: break; // Temperature sens. 2     
      case 0x3F: break; // U  2.5 altro digital (m)
      case 0x40: break; // U  2.5 altro analog (m) 
      case 0x41: break; // Forth ADC T  
      case 0x42: break; // U  2.5 VA (m)
      case 0x43: break; // U  1.5 VA (m)
      case 0x44: break; // U -2.0 VA (m)
      case 0x45: break; // U -2.0 (m)   
      case 0x46: break; // Forth ADC T  
      case 0x47: break; // U  2.5 VA (m)
      case 0x48: break; // U  1.5 VA (m)
      case 0x49: break; // U -2.0 VA (m)
      case 0x4A: break;  // U -2.0 (m)   
	// Counters 
      case 0x0B: break; // L1 trigger CouNTer
      case 0x0C: break; // L2 trigger CouNTer
      case 0x0D: break; // Sampling CLK CouNTer
      case 0x0E: break; // DSTB CouNTer
	// Test mode 
      case 0x0F: break; // Test mode word
      case 0x10: break; // Undersampling ratio.
	// Configuration and status 
      case 0x11: break; // Config/Status Register 0
      case 0x12: break; // Config/Status Register 1
      case 0x13: break; // Config/Status Register 2
      case 0x14: break; // Config/Status Register 3
      case 0x15: break; // Free
	// Comands:
      case 0x16: break; // Latch L1, L2, SCLK Counters
      case 0x17: break; // Clear counters
      case 0x18: break; // Clear CSR1
      case 0x19: break; // rstb ALTROs
      case 0x1A: break; // rstb BC
      case 0x1B: break; // Start conversion
      case 0x1C: break; // Scan event length
      case 0x1D: break; // Read event length
      case 0x1E: break; // Start test mode
      case 0x1F: break; // Read acquisition memory
	// FMD
      case 0x20: break; // FMDD status
      case 0x21: break; // L0 counters
      case 0x22: break; // FMD: Wait to hold
      case 0x23: break; // FMD: L1 timeout
      case 0x24: break; // FMD: L2 timeout
      case 0x25: // FMD: Shift clk 
	shift_clk[board] = ((data >> 8 ) & 0xFF); 
	break; 
      case 0x26: // FMD: Strips 
	strip_low[board]  = ((data >> 0 ) & 0xFF); 
	strip_high[board] = ((data >> 8 ) & 0xFF);  
	break; 
      case 0x27: // FMD: Cal pulse 
	pulse_size[board]  =  ((data >> 8 ) & 0xFF);
	break; 
      case 0x28: break; // FMD: Shape bias
      case 0x29: break; // FMD: Shape ref
      case 0x2A: break; // FMD: Preamp ref
      case 0x2B: // FMD: Sample clk 
	sample_clk[board] = ((data >> 8 ) & 0xFF); 
	break; 
      case 0x2C: break; // FMD: Commands
      case 0x4B: // FMD: Cal events 
	pulse_length[board] = ((data >> 0 ) & 0xFF);
	break; 
      default: break;
	
      }
      AliFMDDebug(50, ("instruction 0x%x, dataword 0x%x",
		       instruction,dataWord));
    }
    
    UShort_t det,sector;
    Short_t strip;
    Char_t ring;
   
    const UInt_t boards[4] = {0,1,16,17};
    for(Int_t i=0;i<4;i++) {
      if(ddl==0 && (i==1 || i==3)) continue;

      UInt_t chip =0, channel=0;
      det = map->DDL2Detector(ddl);
      map->Channel2StripBase(boards[i], chip, channel, ring, sector, strip);
     
      UInt_t samplerate = 1;
      if(sample_clk[boards[i]] == 0) {
	if(ddl == 0) {
	  Int_t sample1 = sample_clk[boards[0]];
	  Int_t sample2 = sample_clk[boards[2]];	    
	  if(sample1) sample_clk[boards[i]] = sample1;
	  else sample_clk[boards[i]] = sample2;
	}
	else {
	  Int_t sample1 = sample_clk[boards[0]];
	  Int_t sample2 = sample_clk[boards[1]];
	  Int_t sample3 = sample_clk[boards[2]];
	  Int_t sample4 = sample_clk[boards[3]];
	  Int_t agreement = 0;
	  if(sample1 == sample2) agreement++;
	  if(sample1 == sample3) agreement++;
	  if(sample1 == sample4) agreement++;
	  if(sample2 == sample3) agreement++;
	  if(sample2 == sample4) agreement++;
	  if(sample3 == sample4) agreement++;
	    
	  Int_t idx = 0;
	  if(i<3) idx = i+1;
	  else  idx = i-1;
	  if(agreement == 3) {
	    sample_clk[boards[i]] = sample_clk[boards[idx]];
	    shift_clk[boards[i]] = shift_clk[boards[idx]];
	    strip_low[boards[i]] = strip_low[boards[idx]];
	    strip_high[boards[i]] = strip_high[boards[idx]];
	    pulse_length[boards[i]] = pulse_length[boards[idx]];
	    pulse_size[boards[i]] = pulse_size[boards[idx]];
	    AliFMDDebug(0, ("Vote taken for ddl %d, board 0x%x",
			    ddl,boards[i]));
	  }
	}
      } 
      
      if(sample_clk[boards[i]])
	samplerate = shift_clk[boards[i]]/sample_clk[boards[i]];
      sampleRate->Set(det,ring,sector,0,samplerate);
      stripRange->Set(det,ring,sector,0,
		      strip_low[boards[i]],strip_high[boards[i]]);
      
      AliFMDDebug(20, ("det %d, ring %c, ",det,ring));
      pulseLength.AddAt(pulse_length[boards[i]],
			GetHalfringIndex(det,ring,boards[i]/16));
      pulseSize.AddAt(pulse_size[boards[i]],
		      GetHalfringIndex(det,ring,boards[i]/16));
      
      
      
      AliFMDDebug(20, (": Board: 0x%02x\n"
		       "\tstrip_low  %3d, strip_high   %3d\n"
		       "\tshift_clk  %3d, sample_clk   %3d\n"
		       "\tpulse_size %3d, pulse_length %3d",
		       boards[i], 
		       strip_low[boards[i]], strip_high[boards[i]],
		       shift_clk[boards[i]], sample_clk[boards[i]],
		       pulse_size[boards[i]],pulse_length[boards[i]]));
    }
    
  }
  
  AliFMDParameters::Instance()->SetSampleRate(sampleRate);
  AliFMDParameters::Instance()->SetStripRange(stripRange);
  
  AliFMDDebug(0, ("End of SOD/EOD"));
  
  return kTRUE;
}
//____________________________________________________________________

UInt_t AliFMDRawReader::Get32bitWord(Int_t idx)
{
  // This method returns the 32 bit word at a given
  // position inside the raw data payload.
  // The 'index' points to the beginning of the next word.
  // The method is supposed to be endian (platform)
  // independent.
  if (!fData) {
    AliFatal("Raw data paylod buffer is not yet initialized !");
  }

  Int_t index = 4*idx;
  
  if (index < 4) {
    //  fRawReader->AddFatalErrorLog(k32bitWordReadErr,Form("pos = %d",index));
    //    PrintDebug();
    AliWarning(Form("Invalid raw data payload position (%d) !",index));
  }

  UInt_t word = 0;
   
  word  = fData[--index] << 24;
  word |= fData[--index] << 16;
  word |= fData[--index] << 8;
  word |= fData[--index] << 0 ;

  return word;
}
//_____________________________________________________________________ 
Int_t AliFMDRawReader::GetHalfringIndex(UShort_t det, Char_t ring, 
					UShort_t board) {

  UShort_t iring  =  (ring == 'I' ? 1 : 0);
  
  Int_t index = (((det-1) << 2) | (iring << 1) | (board << 0));
  
  return index-2;
  
}
//____________________________________________________________________
// 
// EOF
//
