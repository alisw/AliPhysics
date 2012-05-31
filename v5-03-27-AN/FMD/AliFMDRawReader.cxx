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
// #include "AliFMDRawStream.h"	// ALIFMDRAWSTREAM_H 
#include "AliRawReader.h"	// ALIRAWREADER_H 
#include "AliFMDRawReader.h"	// ALIFMDRAWREADER_H 
#include "AliFMDDebug.h"
#include "AliFMDCalibSampleRate.h"
#include "AliFMDCalibStripRange.h"
#include "AliFMDAltroMapping.h"
#include "AliFMDUShortMap.h"
// #include "AliFMDAltroIO.h"	// ALIFMDALTROIO_H 
#include "AliAltroRawStreamV3.h"
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
    fMinStrip(0),
    fMaxStrip(127), 
    fPreSamp(14+5),
    fSeen(0)
{
  // Default CTOR
  for (Int_t i = 0; i < 3; i++) { 
    fSampleRate[i]   = 0;
    fZeroSuppress[i] = kFALSE;
    fNoiseFactor[i]  = 1;
    fL1Phase[i]      = 0;
    fNErrors[i]      = 0;
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
  AliDebug(1,Form("Got a grand total of %d digits, wrote %d bytes to tree", 
		   array->GetEntriesFast(), nWrite));
  delete array;
}

//____________________________________________________________________
Int_t
AliFMDRawReader::NewDDL(AliAltroRawStreamV3& input, UShort_t& det)
{
  // Process a new DDL.  Sets the internal data members fZeroSuppress, 
  // fSampleRate, and fNoiseFactor based on information in the RCU trailer. 
  // 
  // Parameters:
  //   input Input stream
  //   det   On return, the detector number
  // 
  // Return value:
  //   negative value in case of problems, the DDL number otherwise

  // Get the DDL number
  UInt_t ddl = input.GetDDLNumber();
  AliDebug(2,Form("DDL number %d", ddl));

  // Note, previously, the ALTROCFG1 register was interpreted as 
  // 
  // Bits    Value    Description
  //   0- 3      0/1   1st Baseline filter, mode 
  //   4- 5   Over-1   2nd baseline filter, # of pre-samples
  //   6- 9   factor   2nd baseline filter, # of post-samples 
  //  10-          0   2nd baseline filter, enable
  //  11-12       00   Zero suppression, glitch filter mode
  //  13-15      001   Zero suppression, # of post samples
  //  16-17       01   Zero suppression, # of pre  samples
  //  18         0/1   Zero suppression, enable
  //
  // The interpretation used in AliAltroRawStreamerV3 - which
  // corresponds directly to ALTRO DPCFG register - is 
  //
  // Bits    Value  Description
  //   0- 3    0/1   1st Baseline filter, mode 
  //   4         0   Polarity (if '1', then "1's inverse")
  //   5- 6     01   Zero suppression, # of pre samples
  //   7-10   0001   Zero suppression, # of post samples
  //  11         0   2nd baseline filter, enable
  //  12-13     00   Zero suppression, glitch filter mode
  //  14-16 factor   2nd baseline filter, # of post-samples
  //  17-18     01   2nd baseline filter, # of pre-samples 
  //  19       0/1   Zero suppression, enable
  //
  //  Writing 'x' for variable values, that means we have the
  //  following patterns for the 2 cases 
  //
  //    bit #  20   16   12    8    4    0
  //     old    |0x01|0010|00xx|xxxx|xxxx|
  //     new    |x01x|xx00|0000|1010|xxxx|
  //
  //  That means that we can check if bits 10-13 are '1000' or
  //  '0000', which will tell us if the value was written with the
  //  new or the old interpretation.    That is, we can check that 
  //
  //    if (((altrocfg1 >> 10) & 0x8) == 0x8) { 
  //      // old interpretation 
  //    }
  //    else { 
  //      // New interpretation 
  //    }
  //
  // That means, that we should never 
  //
  //  - change the # of zero suppression post samples 
  //  - Turn on 2nd baseline correction 
  //  - Change the zero-suppression glitch filter mode
  //
  // This change as introduced in version 1.2 of Rcu++
  //
  UInt_t cfg1 = input.GetAltroCFG1();
  if (((cfg1 >> 10) & 0x8) == 0x8) {
    UInt_t cfg2 = input.GetAltroCFG2();
    AliDebug(3,Form("We have data from older MiniConf 0x%x cfg2=0x%08x", 
		    ((cfg1 >> 10) & 0x8), cfg2));
    fZeroSuppress[ddl] = (cfg1 >>  0) & 0x1;
    fNoiseFactor[ddl]  = (cfg1 >>  6) & 0xF;
    fSampleRate[ddl]   = (cfg2 >> 20) & 0xF;
  }
  else {
    AliDebug(3,Form("We have data from newer MiniConf 0x%x", 
		    ((cfg1 >> 10) & 0x8)));
    fZeroSuppress[ddl] = input.GetZeroSupp();
    // WARNING: We store the noise factor in the 2nd baseline
    // filters excluded post samples, since we'll never use that
    // mode. 
    fNoiseFactor[ddl]  = input.GetNPostsamples();
    // WARNING: We store the sample rate in the number of pre-trigger
    // samples, since we'll never use that mode.
    fSampleRate[ddl]     = input.GetNPretriggerSamples();
    // 
  }
  AliDebug(10,Form("Phase of DDL=%d is %g (%d)", ddl, input.GetL1Phase(),
		   input.GetAltroCFG2() & 0x1F));
  fL1Phase[ddl] = input.GetAltroCFG2() & 0x1F; // input.GetL1Phase();
  AliDebug(3,Form("RCU @ DDL %d zero suppression: %s", 
		   ddl, (fZeroSuppress[ddl] ? "yes" : "no")));
  AliDebug(3,Form("RCU @ DDL %d noise factor: %d", ddl,fNoiseFactor[ddl]));    
  AliDebug(3,Form("RCU @ DDL %d sample rate: %d", ddl,fSampleRate[ddl]));


  // Get Errors seen 
  Int_t nChAddrMismatch = input.GetNChAddrMismatch();
  Int_t nChLenMismatch  = input.GetNChLengthMismatch();
  if (nChAddrMismatch != 0) 
    AliWarning(Form("Got %d channels with address mis-matches for 0x%03x",
		    nChAddrMismatch, ddl));
  if (nChLenMismatch != 0) 
    AliWarning(Form("Got %d channels with length mis-matches for 0x%03x",
		    nChLenMismatch, ddl));

  // Map DDL number to the detector number 
  AliFMDParameters*    pars   = AliFMDParameters::Instance();
  AliFMDAltroMapping*  map    = pars->GetAltroMap();
  if (map->DDL2Detector(ddl) < 0) return -1;
  det = map->DDL2Detector(ddl);

  if (AliLog::GetDebugLevel("FMD", 0) > 5) 
    input.PrintRCUTrailer();
  return ddl;
}

//____________________________________________________________________
Int_t
AliFMDRawReader::NewChannel(const AliAltroRawStreamV3& input,  UShort_t det,
			    Char_t&  ring, UShort_t& sec, Short_t& strbase)
{
  // Processs a new channel.  Sets the internal data members
  // fMinStrip, fMaxStrip, and fPreSamp. 
  //
  // Parameter:
  //   input   Input stream
  //   ring    On return, the ring identifier 
  //   sec     On return, the sector number
  //   strbase On return, the strip base
  // 
  // Return value
  //   negative value in case of problems, hardware address otherwise

  // Get the hardware address, and map that to detector coordinates 
  UShort_t board, chip, channel;
  Int_t    ddl    = input.GetDDLNumber();
  Int_t    hwaddr = input.GetHWAddress();
  if (input.IsChannelBad()) { 
    const char* msg = Form("Ignoring channel %03d/0x%03x with errors", 
			   ddl, hwaddr); 
    AliWarning(msg); 
    if (AliDebugLevel() > 10) input.HexDumpChannel();
    fReader->AddMinorErrorLog(AliAltroRawStreamV3::kAltroPayloadErr,msg);
    fNErrors[ddl] += 1;
    return 0xFFFF;
  }
  
  AliFMDParameters*    pars   = AliFMDParameters::Instance();
  AliFMDAltroMapping*  map    = pars->GetAltroMap();
  // Map to hardware stuff 
  map->ChannelAddress(hwaddr, board, chip, channel);
  // Then try to map to detector address
  if (!map->Channel2StripBase(board, chip, channel, ring, sec, strbase)) {
    AliError(Form("Failed to get detector id from DDL %d, "
		  "hardware address 0x%03x", ddl, hwaddr));
    return -1;
  }
  AliDebug(4,Form("Board: 0x%02x, Altro: 0x%x, Channel: 0x%x", 
		  board, chip, channel));

  // Get the 'conditions'
  fMinStrip = pars->GetMinStrip(det, ring, sec, strbase);
  fMaxStrip = pars->GetMaxStrip(det, ring, sec, strbase);
  fPreSamp  = pars->GetPreSamples(det, ring, sec, strbase);
  if (fSampleRate[ddl] == 0) {
    AliDebug(3,Form("Get sample rate for RCU @ DDL %d from OCDB", ddl));
    fSampleRate[ddl] = pars->GetSampleRate(det, ring, sec, strbase);
  }
  AliDebug(3,Form("RCU @ DDL %d sample rate: %d", ddl,fSampleRate[ddl]));

  return hwaddr;
}

//____________________________________________________________________
Bool_t
AliFMDRawReader::NewBunch(const AliAltroRawStreamV3& input, 
			  UShort_t&  start, UShort_t& length)
{
  // 
  // Do some checks on the bunch data 
  // 
  Int_t    ddl      = input.GetDDLNumber();
  Int_t    hwaddr   = input.GetHWAddress();  
  UShort_t nSamples = input.GetNSamplesPerCh() + fPreSamp;
  UShort_t tstart   = input.GetStartTimeBin();
  length            = input.GetBunchLength();

  if (tstart >= nSamples) {
    const char* msg = Form("Bunch in %03d/0x%03x has an start time greater "
			   "than number of samples: 0x%x >= 0x%x", 
			   ddl, hwaddr, tstart, nSamples);
    AliWarning(msg);
    if (AliDebugLevel() > 10) input.HexDumpChannel();
    fReader->AddMinorErrorLog(AliAltroRawStreamV3::kAltroPayloadErr,msg);
    fNErrors[ddl]++;
    return false;
  }
  if ((int(tstart) - length + 1) < 0) { 
    const char* msg = Form("Bunch in %03d/0x%03x has an invalid length and "
			   "start time: 0x%x,0x%x (%d-%d+1=%d<0)", 
			   ddl, hwaddr, length, tstart, tstart, length, 
			   int(tstart)-length+1);
    AliWarning(msg);
    if (AliDebugLevel() > 10) input.HexDumpChannel();
    fReader->AddMinorErrorLog(AliAltroRawStreamV3::kAltroPayloadErr,msg);
    fNErrors[ddl]++;				
    return false;
  }
  if (tstart >= start) { 
    const char* msg = Form("Bunch in %03d/0x%03x has early start time: "
			   "0x%x >= 0x%x", ddl, hwaddr, tstart, start);
    AliWarning(msg);
    if (AliDebugLevel() > 10) input.HexDumpChannel();
    fReader->AddMinorErrorLog(AliAltroRawStreamV3::kAltroPayloadErr,msg);
    fNErrors[ddl]++;
    return false;
  }
  start = tstart;
  return true;
}

//____________________________________________________________________
Int_t
AliFMDRawReader::NewSample(const AliAltroRawStreamV3& input, 
			   Int_t i, UShort_t t, UShort_t sec,
			   UShort_t  strbase, Short_t&  str, UShort_t& samp)
{
  // Process a new timebin
  // 
  // Parameters:
  //   input   Input stream
  //   i       Index into bunch data
  //   t       Time
  //   strbase Base of strip numbers for this channel
  //   str     On return, the strip number
  //   samp    On return, the sample number
  // 
  // Return value
  //   negative value in case of problems, ADC value otherwise
  if (t < fPreSamp) return -1;

  Int_t           ddl  = input.GetDDLNumber();
  Int_t           hwa  = input.GetHWAddress();
  const UShort_t* data = input.GetSignals();
  Short_t         adc  = data[i];
  AliDebug(10,Form("0x%04x/0x%03x/%04d %4d", ddl, hwa, t, adc));

  AliFMDParameters*    pars   = AliFMDParameters::Instance();
  AliFMDAltroMapping*  map    = pars->GetAltroMap();

  samp            = 0;
  Short_t  stroff = 0;
  map->Timebin2Strip(sec, t, fPreSamp, fSampleRate[ddl], stroff, samp);
  str             = strbase + stroff;
      
  AliDebug(20,Form("0x%04x/0x%03x/%04d=%4d maps to strip %3d sample %d " 
		   "(pre: %d, min: %d, max: %d, rate: %d)",
		  ddl, hwa, t, adc, str, samp, fPreSamp, 
		  fMinStrip, fMaxStrip, fSampleRate[ddl]));
  if (str < 0) { 
    AliDebug(10,Form("Got presamples at timebin %d", i));
    return -1;
  }
	  
  // VA1 Local strip number 
  Short_t lstrip = (t - fPreSamp) / fSampleRate[ddl] + fMinStrip;
      
  AliDebug(15,Form("Checking if strip %d (%d) in range [%d,%d]", 
		   lstrip, str, fMinStrip, fMaxStrip));
  if (lstrip < fMinStrip || lstrip > fMaxStrip) {
    AliDebug(10,Form("Strip %03d-%d (%d,%d) from t=%d out of range (%3d->%3d)", 
		    str, samp, lstrip, stroff, t, fMinStrip, fMaxStrip));
    adc = -1;
  }
  // Possibly do pedestal subtraction of signal 
  if (adc > 1023) 
    AliWarning(Form("ADC value out of range: %4d", adc));
  return adc;
}

//____________________________________________________________________
Int_t
AliFMDRawReader::NextSample(UShort_t& det, Char_t&   rng, UShort_t& sec, 
			    UShort_t& str, UShort_t& sam, UShort_t& rat, 
			    Short_t&  adc, Bool_t&   zs,  UShort_t& fac)
{
  // Scan current event for next signal.   It returns kFALSE when
  // there's no more data in the event. 
  // 
  // Note, that this member function is in principle very fast, but
  // contains less error checking.  In particular, channels that have
  // bad bunches cannot be checked here.  Seeing a bad bunch will only
  // skip the remainder of the channel and not reset the already read
  // digits.   This is potentially dangerous. 
  //
  // Parameters: 
  //    det         On return, contain the detector number 
  //    rng         On return, contain the ring identifier 
  //    sec         On return, contain the sector number 
  //    str         On return, contain the strip number 
  //    sam         On return, contain the sample number 
  //    rat         On return, contain the sample rate 
  //    adc         On return, contain the ADC counts 
  //    zs          On return, contain the zero-supp. flag 
  //    fac         On return, contain the zero-supp. noise factor 
  // 
  // Return values: 
  //    0    No more data 
  //    -1   Read sample belongs to a bad bunch 
  //    >0   Good status - contains bit mask of values 
  //       Bit 1    New DDL
  //       Bit 2    New Channel
  //       Bit 3    New Bunch
  //       Bit 4    New Sample
  static AliAltroRawStreamV3 stream(fReader); //    = 0;
  static Int_t               ddl      = -1;
  static UShort_t            tdet     = 0;
  static Char_t              trng     = '\0';
  static UShort_t            tsec     = 0;
  static Short_t             tstr     = 0;   
  static Short_t             bstr     = -1;
  static UShort_t            tsam     = 0;   
  // static UInt_t           trate    = 0;
  static Int_t               hwaddr   = -1;
  static UShort_t            start    = 0;
  static UShort_t            length   = 0;
  static Short_t             t        = -1;
  static Int_t               i        = 0; 
  // First entry!
  if (stream.GetDDLNumber() < 0) { 
    fReader->Reset();
    fReader->Select("FMD");
    stream.Reset();
    stream.SelectRawData("FMD");
    stream.SetCheckAltroPayload(false);
    for (Int_t j = 0; j < kNDDL; j++) fNErrors[j] = 0;

    // Reset variables
    ddl    = -1;  
    // trate= 0;   
    tdet   = 0;   
    trng   = '\0';
    tsec   = 0;   
    tstr   = 0;  
    tsam   = -1;
    hwaddr = -1;
  }

  UShort_t ret = 0;
  do { 
    AliDebug(15,Form("t=%4d, start=%4d, length=%4d", t, start, length));
    if (t < start - length + 1) { 
      AliDebug(10,Form("Time t=%d < start-length+1=%d-%d+1 (%3d/0x%03x)", 
		       t, start, length, ddl, hwaddr));
      if (hwaddr > 0xFFF || 
	  hwaddr < 0 || 
	  !stream.NextBunch()) { 
	if (AliDebugLevel() >= 10 && hwaddr > 0xFFF) {
	  AliDebug(10,"Last channel read was marked bad");
	}
	if (AliDebugLevel() >= 10 && hwaddr < 0) {
	  AliDebug(10,"No more channels");
	}
	AliDebug(10,"No next bunch, or first entry");
	if (ddl < 0 || !stream.NextChannel()) { 
	  if (AliDebugLevel() >= 10 && ddl < 0) { 
	    AliDebug(10,"No DDL");
	  }
	  AliDebug(10,"No next channel, or first entry");
	  if (!stream.NextDDL()) {
	    AliDebug(10,"No more DDLs");
	    stream.Reset();
	    return 0;
	  }
	  ddl = NewDDL(stream, tdet);
	  AliDebug(5,Form("New DDL: %d (%d)", ddl, tdet));
	  ret |= 0x1;
	  continue;
	}
	hwaddr = NewChannel(stream, tdet, trng, tsec, bstr);
	if (hwaddr > 0xFFF) fNErrors[ddl] += 1;
	AliDebug(5,Form("New Channel: %3d/0x%03x", ddl, hwaddr));
	start  = 1024;
	ret |= 0x2;
	continue;
      }
      if (!NewBunch(stream, start, length)) { 
	// AliWarning(Form("Bad bunch in %3d/0x%03x read - "
	//                 "should progress to next channel "
	//                 "(t=%4d,start=%4d,length=%4d)", 
	//                 ddl, hwaddr, t,start, length));
	hwaddr = 0xFFFF; // Bad channel
	return -1;
      }
      AliDebug(5, Form("New bunch in  %3d/0x%03x: start=0x%03x, length=%4d", 
		       ddl, hwaddr, start, length));
      ret |= 0x4;
      t      = start;
      i      = 0;
      AliDebug(10,Form("Got new bunch FMD%d%c[%2d], bunch @ %d, length=%d", 
		       tdet, trng, tsec, start, length));
    }
    Int_t tadc = NewSample(stream, i, t, tsec, bstr, tstr, tsam);
    AliDebug(10,Form("New sample FMD%d%c[%2d,%3d]-%d = 0x%03x", 
		 tdet, trng, tsec, tstr, tsam, tadc));
    ret |= 0x8;
    if (tadc >= 0) { 
      det = tdet;
      rng = trng;
      sec = tsec;
      str = tstr;
      sam = tsam;
      adc = tadc;
      rat = fSampleRate[ddl];
      zs  = fZeroSuppress[ddl];
      fac = fNoiseFactor[ddl];
      t--;
      i++;
      AliDebug(10,Form("Returning FMD%d%c[%2d,%3d]-%d = 0x%03x (%d,%d,%d)",
		   det, rng, sec, str, sam, adc, rat, zs, fac));
      break;
    }
    t--;
    i++;
  } while (true);
  AliDebug(5,Form("Returning 0x%02x", ret));
  return ret;
}


//____________________________________________________________________
Int_t
AliFMDRawReader::NextSignal(UShort_t& det, Char_t&   rng, 
			    UShort_t& sec, UShort_t& str, 
			    Short_t&  adc, Bool_t&   zs, 
			    UShort_t& fac)
{
  // 
  // Get the next signal
  // 
  // Parameters:
  //    det  On return, the detector
  //    rng  On return, the ring
  //    sec  On return, the sector
  //    str  On return, the strip
  //    adc  On return, the ADC value
  //    zs   On return, whether zero-supp. is enabled
  //    fac  On return, the usd noise factor
  // 
  // Return:
  //    true if valid data is returned
  //
  Int_t ret = 0;
  do { 
    UShort_t samp, rate;
    if ((ret = NextSample(det, rng, sec, str, samp, rate, adc, zs, fac)) <= 0)
      return ret;

    Bool_t take = SelectSample(samp, rate);
    if (!take) continue;
    break;
  } while (true);
  return ret;
}

//____________________________________________________________________
Bool_t
AliFMDRawReader::SelectSample(UShort_t samp, UShort_t rate) 
{
  // Check if the passed sample is the one we need
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
  
//____________________________________________________________________
Bool_t
AliFMDRawReader::ReadAdcs(TClonesArray* array) 
{
  // Read ADC values from raw input into passed TClonesArray of AliFMDDigit
  // objects. 
  AliDebug(3,Form("Reading ADC values into a TClonesArray"));

  // Read raw data into the digits array, using AliFMDAltroReader. 
  if (!array) {
    AliError("No TClonesArray passed");
    return kFALSE;
  }
  const UShort_t kUShortMax = (1 << 16) - 1;
  fSeen.Reset(kUShortMax);
  for (Int_t ddl = 0; ddl < kNDDL; ddl++) fNErrors[ddl] = 0;

  AliAltroRawStreamV3  input(fReader);
  input.Reset();
  input.SetCheckAltroPayload(false);
  input.SelectRawData("FMD");
  
  // Loop over input RORCs
  while (input.NextDDL()) { 
    UShort_t det = 0;
    Int_t    ddl = NewDDL(input, det);
    if (ddl < 0) break;
    fNErrors[ddl] = 0;

    while (input.NextChannel()) { 
      // Get the hardware address, and map that to detector coordinates 
      Char_t   ring;
      UShort_t sec;
      Short_t  strbase;
      Int_t    hwaddr   = NewChannel(input, det, ring, sec, strbase);
      if (hwaddr < 0) break;
      if (hwaddr > 0xFFF) continue;  

      UShort_t start    = 0x3FF;
      Bool_t   errors   = false;
      Int_t    first    = -1;
      Int_t    last     = -1;
      // Loop over bunches 
      while (input.NextBunch()) { 
	// Get Lenght of bunch, and pointer to the data 
	const UShort_t* data   = input.GetSignals();
	UShort_t        length;
	if (!NewBunch(input, start, length)) {
	  errors = true;
	  break;
	}

      
	// Loop over the data and store it. 
	for (Int_t i = 0; i < length; i++) { 
	  // Time 
	  Short_t  str;
	  UShort_t samp;
	  Int_t    t    = start - i;
	  Int_t    adc  = NewSample(input, i, t, sec, strbase, str, samp);
	  if (adc < 0) continue;
	  UShort_t counts = adc;
      
	  AliDebug(10, Form("FMD%d%c[%02d,%03d]-%d: %4d", 
			    det, ring, sec, str, samp, counts));
	  // Check the cache of indicies
	  Int_t idx = fSeen(det, ring, sec, str);
	  AliFMDDigit* digit = 0;
	  if (idx == kUShortMax) { 
	    // We haven't seen this strip yet. 
	    fSeen(det, ring, sec, str) = idx = array->GetEntriesFast();
	    AliDebug(7,Form("making digit @ %5d for FMD%d%c[%2d,%3d]-%d "
			    "from %3d/0x%03x/%4d", 
			    idx, det, ring, sec, str, samp, ddl, hwaddr, t));
	    digit = new ((*array)[idx]) AliFMDDigit(det, ring, sec, str);
	    digit->SetDefaultCounts(fSampleRate[ddl]);
	  }
	  else {
	    digit = static_cast<AliFMDDigit*>(array->At(idx));
	  }
	  if (first < 0) first = idx;
	  last = idx;
	  AliDebug(10, Form("Setting FMD%d%c[%2d,%3d]-%d from timebin "
			    "%4d=%4d (%4d)", det, ring, sec, str, samp, t, 
			    counts, data[i]));
	  digit->SetCount(samp, counts);
	} // for (i)
      } // while (bunch)
      if (errors) { 
	AliWarning(Form("Channel %3d/0x%03x contain errors, "
			"resetting index %d to %d", ddl, hwaddr, first, last));
	if (first >= 0) {
	  for (Int_t i = first; i <= last; i++) { 
	    AliFMDDigit* digit = static_cast<AliFMDDigit*>(array->At(i));
	    for (Int_t j = 0; j < fSampleRate[ddl]; j++) {
	      AliDebug(10,Form("Resetting strip %s=%d",
			       digit->GetName(),digit->Counts()));
	      digit->SetCount(j, kBadSignal);
	    }
	  }
	}
      }
      // if (errors && (AliDebugLevel() > 0)) input.HexDumpChannel();
    } // while (channel)
  } // while (ddl)
  return kTRUE;
}
//____________________________________________________________________
Bool_t
AliFMDRawReader::ReadAdcs(AliFMDUShortMap& map) 
{
  // Read ADC values from raw input into passed TClonesArray of AliFMDDigit
  // objects. 
  AliDebug(3,Form("Reading ADC values into a map"));

  // const UShort_t kUShortMax = (1 << 16) - 1;
  for (Int_t ddl = 0; ddl < kNDDL; ddl++) fNErrors[ddl] = 0;

  AliAltroRawStreamV3  input(fReader);
  input.Reset();
  input.SetCheckAltroPayload(false);
  input.SelectRawData("FMD");
  
  // Loop over input RORCs
  while (input.NextDDL()) { 
    UShort_t det = 0;
    Int_t    ddl = NewDDL(input, det);
    if (ddl < 0) break;
    fNErrors[ddl] = 0;

    while (input.NextChannel()) { 
      // Get the hardware address, and map that to detector coordinates 
      Char_t   ring;
      UShort_t sec;
      Short_t  strbase;
      Int_t    hwaddr   = NewChannel(input, det, ring, sec, strbase);
      if (hwaddr < 0) break;
      if (hwaddr > 0xFFF) continue;  

      UShort_t start    = 0x3FF;
      Bool_t   errors   = false;
      Int_t    first    = -1;
      Int_t    last     = -1;
      // Loop over bunches 
      while (input.NextBunch()) { 
	// Get Lenght of bunch, and pointer to the data 
	// const UShort_t* data   = input.GetSignals();
	UShort_t        length;
	if (!NewBunch(input, start, length)) {
	  errors = true;
	  break;
	}

      
	// Loop over the data and store it. 
	for (Int_t i = 0; i < length; i++) { 
	  // Time 
	  Short_t  str;
	  UShort_t samp;
	  Int_t    t    = start - i;
	  Int_t    adc  = NewSample(input, i, t, sec, strbase, str, samp);
	  if (adc < 0) continue;
	  UShort_t counts = adc;
      
	  AliDebug(10, Form("FMD%d%c[%02d,%03d]-%d: %4d", 
			    det, ring, sec, str, samp, counts));
	  if (SelectSample(samp, fSampleRate[ddl]))
	    map(det,ring,sec,str) = counts; 
	  if (first < 0) first = str;
	  last = str;
	} // for (i)
      } // while (bunch)
      if (errors) { 
	AliWarning(Form("Channel %3d/0x%03x contain errors, "
			"resetting strips %d to %d", ddl, hwaddr, first, last));
	if (first >= 0) {
	  Int_t ds = first <= last ? 1 : -1;
	  for (Int_t i = first; i != last+ds; i += ds) { 
	    AliDebug(10, Form("Resetting strip FMD%d%c[%02d,%03d]=%d",
			      det,ring,sec,i,map(det,ring,sec,i)));
	    map(det,ring,sec,i) = kBadSignal;
	  }
	}
      }
    } // while (channel)
  } // while (ddl)
  return kTRUE;
}

//____________________________________________________________________
Bool_t AliFMDRawReader::ReadSODevent(AliFMDCalibSampleRate* sampleRate, 
				     AliFMDCalibStripRange* stripRange, 
				     TArrayS &pulseSize, 
				     TArrayS &pulseLength, 
				     Bool_t* detectors) 
{
  // 
  // Read SOD event into passed objects.
  // 
  // Parameters:
  //    samplerate   The sample rate object to fill
  //    striprange   The strip range object to fill
  //    pulseSize    The pulse size object to fill
  //    pulseLength  The pulse length (in events) object to fill
  // 
  // Return:
  //    @c true on success
  //  
  AliDebug(0,Form("Start of SOD/EOD"));
  
  UInt_t shift_clk[18];
  UInt_t sample_clk[18];
  UInt_t strip_low[18];
  UInt_t strip_high[18];
  UInt_t pulse_size[18];
  UInt_t pulse_length[18];  
  for (size_t i = 0; i < 18; i++) { 
    shift_clk[i]    = 0;
    sample_clk[i]   = 0;
    strip_low[i]    = 0;
    strip_high[i]   = 0;
    pulse_size[i]   = 0;
    pulse_length[i] = 0;
  }
  AliFMDParameters*   param = AliFMDParameters::Instance();
  AliFMDAltroMapping* map   = param->GetAltroMap();
  
  AliAltroRawStreamV3  streamer(fReader);
  streamer.Reset();
  streamer.SelectRawData("FMD");
  //fReader->GetDDLID();
  while (streamer.NextDDL()) {
    Int_t ddl   = streamer.GetDDLNumber();
    Int_t detID = fReader->GetDetectorID();
    if (detectors) detectors[map->DDL2Detector(ddl)-1] = kTRUE;
    AliDebug(0,Form(" From reader: DDL number is %d , det ID is %d",ddl,detID));
    
    ULong_t  nPayloadWords = streamer.GetRCUPayloadSizeInSOD();
    UChar_t* payloadData   = streamer.GetRCUPayloadInSOD();
    UInt_t*  payloadWords  = reinterpret_cast<UInt_t*>(payloadData);
    //UInt_t*   payloadWords  = streamer.GetRCUPayloadInSOD();

    //std::cout<<nPayloadWords<<"    "<<ddl<<std::endl;
    for (ULong_t i = 1; i <= nPayloadWords ; i++, payloadWords++) {
      UInt_t payloadWord = *payloadWords; // Get32bitWord(i);
    
      //std::cout<<i<<Form("  word: 0x%x",payloadWord)<<std::endl;
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
	payloadWords++;
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

      UInt_t dataWord      = *payloadWords;//Get32bitWord(i);
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
	AliDebug(30,Form("Read shift_clk=%d for board 0x%02x", 
			 shift_clk[board], board));
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
	AliDebug(30,Form("Read sample_clk=%d for board 0x%02x", 
			 sample_clk[board], board));
	break; 
      case 0x2C: break; // FMD: Commands
      case 0x4B: // FMD: Cal events 
	pulse_length[board] = ((data >> 0 ) & 0xFF);
	break; 
      default: break;
	
      }
      AliDebug(50,Form("instruction 0x%x, dataword 0x%x",
		       instruction,dataWord));
    } // End of loop over Result memory event
    
    UShort_t det    = 0;
    UShort_t sector = 0;
    Short_t  strip  = -1;
    Char_t   ring   = '\0';
   
    const UInt_t boards[4] = {0,1,16,17};
    for(Int_t i=0;i<4;i++) {
      if(ddl==0 && (i==1 || i==3)) continue;

      UInt_t chip =0, channel=0;
      det = map->DDL2Detector(ddl);
      map->Channel2StripBase(boards[i], chip, channel, ring, sector, strip);
     
      UInt_t samplerate = 0;
#if USE_VOTE
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
	    AliDebug(3,Form("Vote taken for ddl %d, board 0x%x",
			    ddl,boards[i]));
	  }
	}
      } 
#endif
      
      if(sample_clk[boards[i]])
	samplerate = shift_clk[boards[i]]/sample_clk[boards[i]];
      AliDebug(10,Form("Sample rate for board 0x%02x is %d", 
		      boards[i], samplerate));
      sampleRate->Set(det,ring,sector,0,samplerate);
      stripRange->Set(det,ring,sector,0,
		      strip_low[boards[i]],strip_high[boards[i]]);
      
      AliDebug(20,Form("det %d, ring %c, ",det,ring));
      pulseLength.AddAt(pulse_length[boards[i]],
			GetHalfringIndex(det,ring,boards[i]/16));
      pulseSize.AddAt(pulse_size[boards[i]],
		      GetHalfringIndex(det,ring,boards[i]/16));
      
      
      
      AliDebug(20,Form(": Board: 0x%02x\n"
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
  
  AliDebug(0,Form("End of SOD/EOD"));
  
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
					UShort_t board) const
{
  // 
  // Get short index for a given half-ring
  // 
  // Parameters:
  //    det   Detector number
  //    ring  Ring identifer
  //    board Board number
  // 
  // Return:
  //    
  //
  UShort_t iring  =  (ring == 'I' ? 1 : 0);
  
  Int_t index = (((det-1) << 2) | (iring << 1) | (board << 0));
  
  return index-2;
  
}

//____________________________________________________________________
// 
// EOF
//

