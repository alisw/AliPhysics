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
#include "AliFMDDebug.h" // Better debug macros
#include <AliLoader.h>		// ALILOADER_H
#include <AliAltroBuffer.h>     // ALIALTROBUFFER_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDRawWriter.h"	// ALIFMDRAWREADER_H 
#include "AliFMDAltroMapping.h" // ALIFMDALTROMAPPING_H
// #include "AliFMDAltroIO.h"   // ALIFMDALTROWRITER_H
#include <TArrayI.h>		// ROOT_TArrayI
#include <TClonesArray.h>	// ROOT_TClonesArray
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
  
  TClonesArray* digits = new TClonesArray("AliFMDDigit", 1000);
  fFMD->SetTreeAddress();
  TBranch* digitBranch = digitTree->GetBranch(fFMD->GetName());
  if (!digitBranch) {
    AliError(Form("no branch for %s", fFMD->GetName()));
    return;
  }
  digitBranch->SetAddress(&digits);
  
  Int_t nEvents = Int_t(digitTree->GetEntries());
  AliFMDDebug(5, ("Got a total of %5d events from tree", nEvents));
  for (Int_t event = 0; event < nEvents; event++) {
    fFMD->ResetDigits();
    digitTree->GetEvent(event);
    
    // Write out the digits
    WriteDigits(digits);
  }
  loader->UnloadDigits();
}

#if 1
//____________________________________________________________________
void
AliFMDRawWriter::WriteDigits(TClonesArray* digits)
{
  // WRite an array of digits to disk file 
  Int_t nDigits = digits->GetEntries();
  if (nDigits < 1) return;
  AliFMDDebug(5, ("Got a total of %5d digits from tree", nDigits));

  AliFMDParameters* pars = AliFMDParameters::Instance();
  UShort_t threshold    = 0;
  UInt_t   prevddl      = 0xFFFF;
  UInt_t   prevaddr     = 0xFFF;
  // UShort_t prevStrip    = 0;
  
  // Which channel number in the ALTRO channel we're at 
  UShort_t nWords       = 0;
  UShort_t preSamples   = 0;
  
  // How many times the ALTRO Samples one VA1_ALICE channel 
  Int_t sampleRate      = 1;
  
  // A buffer to hold 1 ALTRO channel - Normally, one ALTRO channel
  // holds 128 VA1_ALICE channels, sampled at a rate of `sampleRate' 
  TArrayI data(pars->GetChannelsPerAltro() * 8);

  // The Altro buffer 
  AliAltroBuffer* altro = 0;
    
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
    UInt_t   ddl;
    UInt_t   addr;  
    if (det != oldDet) {
      AliFMDDebug(5, ("Got new detector: %d (was %d)", det, oldDet));
      oldDet = det;
    }
    AliFMDDebug(10, ("Processing digit # %5d FMD%d%c[%2d,%3d]", 
		    i, det, ring, sector, strip));
    threshold       = pars->GetZeroSuppression(det, ring, sector, strip);
    if (!pars->Detector2Hardware(det, ring, sector, strip, ddl, addr)) {
      AliError(Form("Failed to get hardware address for FMD%d%c[%2d,%3d]", 
		    det, ring, sector, strip));
      continue;
    }
    preSamples = pars->GetPreSamples(det, ring, sector, strip);
    
    AliFMDDebug(10, ("FMD%d%c[%2d,%3d]-> ddl: 0x%x addr: 0x%x", 
		    det, ring, sector, strip, ddl, addr));
    if (addr != prevaddr) {
      // Flush a channel to output 
      AliFMDDebug(15, ("Now hardware address 0x%x from FMD%d%c[%2d,%3d] "
		       "(board 0x%x, chip 0x%x, channel 0x%x), flushing old "
		       "channel at 0x%x with %d words", 
		       addr, det, ring, sector, strip, 
		       (addr >> 7), (addr >> 4) & 0x7, addr & 0xf, 
		       prevaddr, nWords));
      if (altro) altro->WriteChannel(prevaddr,nWords,data.fArray,threshold);
      nWords   = preSamples;
      prevaddr = addr;
      for (size_t j = 0; j < nWords; j++) data[j] = digit->Count(0);
    }
    if (ddl != prevddl) {
      AliFMDDebug(5, ("FMD: New DDL, was %d, now %d", prevddl, ddl));
      // If an altro exists, delete the object, flushing the data to
      // disk, and closing the file. 
      if (altro) { 
	// When the first argument is false, we write the real
	// header. 
	AliFMDDebug(15, ("Closing output"));
	altro->Flush();
	altro->WriteDataHeader(kFALSE, kFALSE);
	delete altro;
	altro = 0;
      }
      prevddl = ddl;
      // Need to open a new DDL! 
      TString filename(AliDAQ::DdlFileName(fFMD->GetName(),  ddl));
      AliFMDDebug(5, ("New altro buffer with DDL file %s", filename.Data()));
      // Create a new altro buffer - a `1' as the second argument
      // means `write mode' 
      altro = new AliAltroBuffer(filename.Data());
      altro->SetMapping(pars->GetAltroMap());      
      // Write a dummy (first argument is true) header to the DDL
      // file - later on, when we close the file, we write the real
      // header
      altro->WriteDataHeader(kTRUE, kFALSE);
    }
    
    // Store the counts of the ADC in the channel buffer 
    sampleRate = pars->GetSampleRate(det, ring, sector, strip);
    for (int s = 0; s < sampleRate; s++) {
      data[nWords] = digit->Count(s);
      nWords++;
    }
  }
  // Finally, we need to close the final ALTRO buffer if it wasn't
  // already 
  if (altro) {
    if (nWords > 0) altro->WriteChannel(prevaddr,nWords,data.fArray,threshold);
    altro->Flush();
    altro->WriteDataHeader(kFALSE, kFALSE);
    delete altro;
  }
}
#else
//____________________________________________________________________
void
AliFMDRawWriter::WriteDigits(TClonesArray* digits)
{
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
