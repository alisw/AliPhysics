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
/** @file    AliFMDRawWriter.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:45:56 2006
    @brief   Class to write raw data 
*/
//____________________________________________________________________
//
// Class to write ADC values to a raw data file
//
// This class writes FMD Raw data to a file.   The sample rate (number
// of times the ALTRO ADC samples each pre-amp. channel - that is,
// data from a single strip), can be set via SetSampleRate. 
//
// Zero-suppression can be enabled by calling SetThreshold with a
// non-zero argument.   ADC values less than the value set will not be
// written to output.   Note, that if you use zero-suppression, you
// need to explicitly set the sample rate when reading back the data
// with AliFMDRawReader. 
// 
// This class uses the AliAltroBuffer class to write the data in the
// ALTRO format.  See the Exec member function for more information on
// that format.  
//
// #include <AliLog.h>		// ALILOG_H
#include "AliFMDDebug.h"        // Better debug macros
#include <AliLoader.h>		// ALILOADER_H
#include <AliAltroBufferV3.h>   // ALIALTROBUFFER_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDRawWriter.h"	// ALIFMDRAWREADER_H 
#include "AliFMDAltroMapping.h" // ALIFMDALTROMAPPING_H
// #include "AliFMDAltroIO.h"   // ALIFMDALTROWRITER_H
#include <TArrayI.h>		// ROOT_TArrayI
#include <TArrayF.h>		// ROOT_TArrayI
#include <TArrayC.h>		// ROOT_TArrayI
#include <TClonesArray.h>	// ROOT_TClonesArray
#include <TTree.h>
// #include <fstream>
#include "AliDAQ.h"

//____________________________________________________________________
ClassImp(AliFMDRawWriter)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDRawWriter::AliFMDRawWriter(AliFMD* fmd) 
  : TTask("FMDRawWriter", "Writer of Raw ADC values from the FMD"),
    fFMD(fmd),
    fSampleRate(0), 
    fChannelsPerAltro(0), 
    fThreshold(0)
{
  // CTOR 
  AliFMDDebug(5, ("Created AliFMDRawWriter object"));
}



//____________________________________________________________________
void
AliFMDRawWriter::Exec(Option_t*) 
{
  // Turn digits into raw data. 
  // 
  // Digits are read from the Digit branch, and processed to make
  // three DDL files, one for each of the sub-detectors FMD1, FMD2,
  // and FMD3. 
  //
  // The raw data files consists of a header, followed by ALTRO
  // formatted blocks.  
  // 
  //          +-------------+
  //          | Header      |
  //          +-------------+
  //          | ALTRO Block |
  //          | ...         |
  //          +-------------+
  //          DDL file 
  // 
  // An ALTRO formatted block, in the FMD context, consists of a
  // number of counts followed by a trailer. 
  // 
  //          +------------------+
  //          | Count            |
  //          | ...              |
  //          | possible fillers |
  //          +------------------+
  //          | Trailer          |
  //          +------------------+
  //          ALTRO block 
  // 
  // The counts are listed backwards, that is, starting with the
  // latest count, and ending in the first. 
  // 
  // Each count consist of 1 or more ADC samples of the VA1_ALICE
  // pre-amp. signal.  Just how many samples are used depends on
  // whether the ALTRO over samples the pre-amp.  Each sample is a
  // 10-bit word, and the samples are grouped into 40-bit blocks 
  //
  //          +------------------------------------+
  //          |  S(1)   | S(2)   | S(3)   | S(4)   |
  //          |  ...    | ...    | ...    | ...    |
  //          |  S(n)   | T(n)   | n+2    | 2AA    |
  //          +------------------------------------+
  //          Counts + possible filler 
  //
  // The trailer of the number of words of signales, the starting
  // strip number, the sector number, and the ring ID; each 10-bit
  // words,  packed into 40-bits. 
  // 
  //          +------------------------------------+
  //          |   2AAA   |  Len   |  A |  Address  |
  //          +------------------------------------+
  //          Trailer
  // 
  // Note, that this method assumes that the digits are ordered. 
  // 
  AliLoader* loader = fFMD->GetLoader();
  loader->LoadDigits("READ");
  TTree* digitTree = loader->TreeD();
  if (!digitTree) {
    AliError("no digit tree");
    return;
  }
  
  fFMD->SetTreeAddress();
  TClonesArray* digits = fFMD->Digits(); 
  // new TClonesArray("AliFMDDigit", 1000);
  // TBranch* digitBranch = digitTree->GetBranch(fFMD->GetName());
  // if (!digitBranch) {
  //   AliError(Form("no branch for %s", fFMD->GetName()));
  //   return;
  // }
  // digitBranch->SetAddress(&digits);
  
  Int_t nEvents = Int_t(digitTree->GetEntries());
  AliFMDDebug(5, ("Got a total of %5d events from tree", nEvents));
  for (Int_t event = 0; event < nEvents; event++) {
    fFMD->ResetDigits();
    digitTree->GetEvent(event);
    
    // Write out the digits
    WriteDigits(digits);
  }
  loader->UnloadDigits();
  //delete digits;
}

#if 1
//____________________________________________________________________
Long_t
AliFMDRawWriter::WriteDigits(TClonesArray* digits)
{
  // WRite an array of digits to disk file 
  Int_t nDigits = digits->GetEntries();
  if (nDigits < 1) return 0;
  AliFMDDebug(5, ("Got a total of %5d digits from tree", nDigits));

  AliFMDParameters* pars = AliFMDParameters::Instance();
  UShort_t threshold    = 0;
  UShort_t factor       = 0;
  UInt_t   prevddl      = 0xFFFF;
  UInt_t   prevaddr     = 0xFFF;
  // UShort_t prevStrip    = 0;
  
  // Which channel number in the ALTRO channel we're at 
  UShort_t nWords       = 0;
  UShort_t preSamples   = 0;
  UShort_t sampleRate   = 0;
  
  // A buffer to hold 1 ALTRO channel - Normally, one ALTRO channel
  // holds 128 VA1_ALICE channels, sampled at a rate of `sampleRate' 
  TArrayI data(pars->GetChannelsPerAltro() * 8);
  TArrayF peds(pars->GetChannelsPerAltro() * 8);
  TArrayF noise(pars->GetChannelsPerAltro() * 8);

  // The Altro buffer 
  AliAltroBufferV3* altro = 0;
  
  Int_t  totalWords = 0;
  Int_t  nCounts    = 0;
  Long_t nBits      = 0;

  // Loop over the digits in the event.  Note, that we assume the
  // the digits are in order in the branch.   If they were not, we'd
  // have to cache all channels before we could write the data to
  // the ALTRO buffer, or we'd have to set up a map of the digits. 
  UShort_t oldDet = 1000;
  for (Int_t i = 0; i < nDigits; i++) {
    // Get the digit
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));
    UShort_t det    = digit->Detector();
    Char_t   ring   = digit->Ring();
    UShort_t sector = digit->Sector();
    UShort_t strip  = digit->Strip();
    UShort_t ddl, addr, time;

    AliFMDDebug(15, ("Processing digit # %5d FMD%d%c[%2d,%3d]", 
		    i, det, ring, sector, strip));
    threshold  = pars->GetZeroSuppression(det, ring, sector, strip);
    sampleRate = pars->GetSampleRate(det, ring, sector, strip);
    preSamples = pars->GetPreSamples(det, ring, sector, strip);
    factor     = UShort_t(pars->GetPedestalFactor());

    if (det != oldDet) {
      AliFMDDebug(5, ("Got new detector: %d (was %d)", det, oldDet));
      oldDet = det;      
    }
    AliFMDDebug(15, ("Sample rate is %d", sampleRate));
    
    for (UShort_t j = 0; j < sampleRate; j++) { 
      if (!pars->Detector2Hardware(det,ring,sector,strip,j,ddl,addr,time)){
	AliError(Form("Failed to get hardware address for FMD%d%c[%2d,%3d]-%d", 
		      det, ring, sector, strip, j));
	continue;
      }
    
      AliFMDDebug(10, ("FMD%d%c[%2d,%3d]-%d-> 0x%x/0x%x/%04d", 
		       det, ring, sector, strip, j, ddl, addr, time));
      if (addr != prevaddr) {
	// Flush a channel to output 
	AliFMDDebug(5, ("Now hardware address 0x%x from FMD%d%c[%2d,%3d]-%d"
			 "(b: 0x%02x, a: 0x%01x, c: 0x%02x, t: %04d), "
			 "flushing old channel at 0x%x with %d words", 
			 addr, det, ring, sector, strip, j,
			 (addr >> 7), (addr >> 4) & 0x7, addr & 0xf, 
			 time, prevaddr, nWords));
	totalWords += nWords;
	ZeroSuppress(data.fArray, nWords, peds.fArray, noise.fArray, threshold);
	if (altro) 
	  /*nBits+=*/altro->WriteChannel(prevaddr,nWords,data.fArray,threshold);
	data.Reset(-1);
	peds.Reset(0);
	noise.Reset(0);
	nWords   = 0;
	prevaddr = addr;
      }
      if (ddl != prevddl) {
	AliFMDDebug(10, ("FMD: New DDL, was %d, now %d", prevddl, ddl));
	// If an altro exists, delete the object, flushing the data to
	// disk, and closing the file. 
	if (altro) { 
	  // When the first argument is false, we write the real
	  // header. 
	  AliFMDDebug(15, ("Closing output"));
	  WriteRCUTrailer(altro, prevddl, threshold > 0, factor, sampleRate); 
    
	  delete altro;
	  altro = 0;
	}
	prevddl = ddl;
	// Need to open a new DDL! 
	TString filename(AliDAQ::DdlFileName(fFMD ? fFMD->GetName() : "FMD",  ddl));
	AliFMDDebug(5, ("New altro buffer with DDL file %s", filename.Data()));
	// Create a new altro buffer - a `1' as the second argument
	// means `write mode' 
	altro = new AliAltroBufferV3(filename.Data());
	altro->SetMapping(pars->GetAltroMap());      
	// Write a dummy (first argument is true) header to the DDL
	// file - later on, when we close the file, we write the real
	// header
	altro->WriteDataHeader(kTRUE, kFALSE);
      }
    
      // Get the pedestal value 
      peds[time]  = pars->GetPedestal(det, ring, sector, strip);
      noise[time] = pars->GetPedestalWidth(det, ring, sector, strip);

      // Store the counts of the ADC in the channel buffer 
      AliFMDDebug(15, ("Storing FMD%d%c[%02d,%03d]-%d in timebin %d (%d)",
		      det, ring, sector, strip, j, time, preSamples));
      data[time] = digit->Count(j);
      nWords++;
      nCounts++;
      if (time == preSamples) {
	AliFMDDebug(15, ("Filling in %4d for %d presamples", 
			data[time], preSamples));
	for (int k = 0; k < preSamples; k++) { 
	  peds[k]  = peds[time];
	  noise[k] = noise[time];
	  data[k]  = data[time];
	}
	nWords += preSamples;
      }
    }
  }
  // Finally, we need to close the final ALTRO buffer if it wasn't
  // already 
  if (altro) {
    ZeroSuppress(data.fArray, nWords, peds.fArray, noise.fArray, threshold);
    if (nWords > 0) 
      /* nBits += */ altro->WriteChannel(prevaddr,nWords,data.fArray,threshold);
    WriteRCUTrailer(altro, prevddl, threshold > 0, factor, sampleRate); 
    delete altro;
  }
  AliFMDDebug(5, ("Wrote a total of %d words in %ld bytes for %d counts", 
		  nWords, nBits / 8, nCounts));
  return nBits;
}
//____________________________________________________________________
void
AliFMDRawWriter::WriteRCUTrailer(AliAltroBufferV3* altro,
				 UInt_t ddl,
				 Bool_t zs,
				 UShort_t factor,
				 UShort_t rate) const
{
  // Flush and write the data header
  altro->Flush();
  altro->WriteDataHeader(kFALSE, kFALSE);
    
  // Set parameters in RCU trailer. 
  // Zero-suppression flag
  altro->SetZeroSupp(zs); // bool
  // WARNING: We store the noise factor in the 2nd baseline
  // filters excluded post samples, since we'll never use that
  // mode. 
  altro->SetNPostsamples(factor); // 
  // WARNING: We store the sample rate in the number of pre-trigger
  // samples, since we'll never use that mode.
  altro->SetNPretriggerSamples(rate); // fSampleRate[ddl]
  // Active front-end cars 
  altro->SetActiveFECsA((ddl == 0 ? 0x1 : 0x3));
  altro->SetActiveFECsB((ddl == 0 ? 0x1 : 0x3));

  // Calculate number of samples 
  altro->SetNSamplesPerCh(rate * 128);
  AliDebug(5,Form("Writing RCU trailer @ DDL %d w/zs=%d, factor=%d, rate=%d",
		  ddl, zs > 0, factor, rate));
  altro->WriteRCUTrailer(ddl);
}

//____________________________________________________________________
void
AliFMDRawWriter::ZeroSuppress(Int_t*& data, Int_t nWords, 
			      const Float_t* peds, 
			      const Float_t* noise, UShort_t threshold) const
{
  // Simulate the ALTRO zero-suppression filter.  The data passed in
  // the first array is modified, such that all suppressed channels
  // are set to some value below threshold.  
  // 
  // If threshold == 0 zero suppression is considered disabled, and no
  // action is taken. 
  if (threshold <= 0) return;

  // const Short_t width  = 3;
  // If fPedSubtract is false, compare data-(ped+f*noise), if true
  // always modify data by -(ped+f*noise), and force negative values
  // to zero.
  Bool_t   pedSubtract = AliFMDParameters::Instance()->IsZSPedSubtract();
  UShort_t pre         = AliFMDParameters::Instance()->GetZSPreSamples();
  UShort_t post        = AliFMDParameters::Instance()->GetZSPostSamples();
  Float_t  factor      = AliFMDParameters::Instance()->GetPedestalFactor();
  
  TArrayC mask(nWords+1);
  for (Short_t i = 0; i < nWords; i++) { 
    Float_t            val     = data[i] - peds[i] - factor * noise[i];
    if (val < 0.5)     val     = 0;
    if (pedSubtract)   data[i] = Int_t(val) & 0x3FF;

    mask[i] = (val > threshold ? 1 : 0);
    AliFMDDebug(10, ("Comparing sample %d %d-%f-%f*%f=%f to %d -> %s", 
		     i, data[i], peds[i], factor, noise[i], val, threshold, 
		    (mask[i] ? "above" : "below")));
  }
  
  for (Short_t i = 0; i < nWords; i++) { 
    if (mask[i]) { // Signal above, so do nothing 
      AliFMDDebug(10, ("Sample %d, above", i));
      if (i < nWords-1 && !mask[i+1]) { 
	// After a valid sample.  Increase the pointer to the next
	// possible data, thereby skipping over the post-samples 
	AliFMDDebug(10, ("At sample %d, next is below, skipping %d to %d", 
			i, post, i+post));
	i += post;
      }
      continue;
    }
    
    Short_t lookahead = TMath::Min(Short_t(nWords), Short_t(i+pre));
    AliFMDDebug(10, ("Sample %d, below, look to %d", i, lookahead));
    if (mask[lookahead] && pre > 0) { 
      AliFMDDebug(10, ("Sample %d is a pre-sample to %d", i, lookahead));
      // We're in a presample, so don't modify the data, and increase
      // counter by the number of pre-samples 
      i += pre-1;
      continue;
    }
    
    // This sample must be surpressed 
    data[i] = threshold - 1;
  }
}

#else
//____________________________________________________________________
void
AliFMDRawWriter::WriteDigits(TClonesArray* digits)
{
  // 
  // Write digits to file 
  //
  Int_t nDigits = digits->GetEntries();
  if (nDigits < 1) return;

  AliFMDParameters*  pars    = AliFMDParameters::Instance();
  AliFMDAltroWriter* writer  = 0;
  Int_t          sampleRate  = -1;
  UShort_t       hwaddr      = 0;
  UShort_t       ddl         = 0;
  std::ofstream* file        = 0;
  Int_t          ret         = 0;
  for (Int_t i = 0; i < nDigits; i++) {
    // Get the digit
    AliFMDDigit* digit = static_cast<AliFMDDigit*>(digits->At(i));
    UInt_t   thisDDL, thisHwaddr;
    UShort_t det    = digit->Detector();
    Char_t   ring   = digit->Ring();
    UShort_t sector = digit->Sector();
    UShort_t strip  = digit->Strip();
    if (!pars->Detector2Hardware(det,ring,sector,strip,thisDDL,thisHwaddr)) {
      AliError(Form("Failed to get hardware address for FMD%d%c[%2d,%3d]",
		    det, ring, sector, strip));
      continue;
    }
    AliFMDDebug(40, ("Got DDL=%d and address=%d from FMD%d%c[%2d,%3d]", 
		     thisDDL, thisHwaddr, det, ring, sector, strip));
    // Check if we're still in the same channel
    if (thisHwaddr != hwaddr) {
      AliFMDDebug(30, ("Now hardware address 0x%x from FMD%d%c[%2d,%3d] "
		       "(board 0x%x, chip 0x%x, channel 0x%x)",
		       thisHwaddr, det, ring, sector, strip, 
		       (thisHwaddr >> 7), (thisHwaddr >> 4) & 0x7, 
		       thisHwaddr & 0xf));
      if (writer) writer->AddChannelTrailer(hwaddr);
      hwaddr = thisHwaddr;
    }
    // Check if we're still in the same detector (DDL)
    if (ddl != thisDDL) {
      if (writer) {
	AliFMDDebug(1, ("Closing altro writer %p", writer));
	if ((ret = writer->Close()) < 0) {
	  AliError(Form("Error: %s", writer->ErrorString(ret)));
	  return;
	}
	delete writer;
	writer = 0;
	file->close();
	delete file;
	file = 0;
      }
      ddl = thisDDL;
    }
    // If we haven't got a writer (either because none were made so
    // far, or because we've switch DDL), make one. 
    if (!writer) {
      AliFMDDebug(1, ("Opening new ALTRO writer w/file %s", 
		      AliDAQ::DdlFileName("FMD",ddl)));
      file   = new std::ofstream(AliDAQ::DdlFileName("FMD",ddl));
      if (!file || !*file) {
	AliFatal(Form("Failed to open file %s", 
		      AliDAQ::DdlFileName("FMD",ddl)));
	return;
      }
      writer  = new AliFMDAltroWriter(*file);
      writer->SetThreshold(pars->GetZeroSuppression(det, ring, sector, strip));
    }
    // Write out our signal
    sampleRate =  pars->GetSampleRate(det,ring,sector,strip);
    writer->AddSignal(digit->Count1());
    if (sampleRate >= 2) writer->AddSignal(digit->Count2());
    if (sampleRate >= 3) writer->AddSignal(digit->Count3());
  }
  if (writer) {
    writer->AddChannelTrailer(hwaddr);
    writer->Close();
    delete writer;
    file->close();
    delete file;
  }
}
#endif


  

//____________________________________________________________________
// 
// EOF
//

