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
/** @file    AliFMDAltroIO.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:27:06 2006
    @brief   Altro Input/Output
*/
//____________________________________________________________________
//                                                                          
// Mapping of ALTRO hardware channel to detector coordinates 
//
#include "AliFMDAltroIO.h"
#include <AliRawDataHeader.h>
#include <AliRawReader.h>
#include "AliLog.h"
#include <iostream>
#include <iomanip>
#define PRETTY_HEX(N,X) \
  "  0x" << std::setfill('0') << std::setw(N) << std::hex << X \
         << std::setfill(' ') << std::dec

//====================================================================
ClassImp(AliFMDAltroIO)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
const AliFMDAltroIO::W40_t AliFMDAltroIO::fgkTrailerMask = 
((AliFMDAltroIO::W40_t(0x2aaa) << 26) + (AliFMDAltroIO::W40_t(0xa) << 12));

//____________________________________________________________________
AliFMDAltroIO::AliFMDAltroIO() 
  : fBuffer(0), fIBuffer(0)
{}

//____________________________________________________________________
const char*
AliFMDAltroIO::ErrorString(Int_t err)  const
{
  switch (err) {
  case kNoError:    return "No error";                          break;
  case kBadFile:    return "Bad state after open/close file";   break;
  case kBadBits:    return "Bad bit offset specified";          break;
  case kBadRead:    return "Bad state after reading from file"; break;
  case kBadWrite:   return "Bad state after writing to file";   break;
  case kBadSeek:    return "Bad state after seeking in file";   break;
  case kBadTell:    return "Could not tell position in file";   break;
  case kBadTrailer: return "Bad trailer 40 bit word in file";   break;
  case kBadFill:    return "Bad fill word in file";             break;
  }
  return "Unknown";
}


//____________________________________________________________________
AliFMDAltroIO::W40_t
AliFMDAltroIO::ConcatW40(UShort_t n, const W10_t& w) const
{
  if (n > 3) return -kBadBits;
  return W40_t(w & 0x3ff) << (10 * n);
}

//____________________________________________________________________
AliFMDAltroIO::W10_t
AliFMDAltroIO::ExtractW10(UShort_t n, const W40_t w) const
{
  if (n > 3) return -kBadBits;
  return (w >> (10 * n)) & 0x3ff;
}

//====================================================================
ClassImp(AliFMDAltroReader)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDAltroReader::AliFMDAltroReader(std::istream& stream)
  : fInput(stream)
  // : fBuffer(buffer), fCurrent(n / 10 * sizeof(char))
{
  // fInput.open(filename);
  if (!fInput)      throw -kBadFile;
  fBegin   = fInput.tellg();
  if (fInput.bad()) throw -kBadTell;
  fInput.seekg(0, std::ios_base::end);
  if (fInput.bad()) throw -kBadSeek;
  fCurrent = fInput.tellg();
  if (fInput.bad()) throw -kBadTell;
#if 0
  fInput.seekg(fBegin);
  UShort_t i = 0;
  do {
    W40_t w = 0;
    fInput.read((char*)&w, 5);
    std::cout << std::setw(6) << i << ": " << PRETTY_HEX(10, w) << std::endl;
    i++;
  } while (!fInput.eof());
  fInput.seekg(fCurrent);
#endif
}

//____________________________________________________________________
Int_t
AliFMDAltroReader::ReadChannel(UShort_t& board, UShort_t& chip, 
			       UShort_t& channel, UShort_t& last, 
			       UShort_t* data) 
{
  UShort_t hwaddr;
  Int_t    ret = ReadChannel(hwaddr, last, data);
  board        = (hwaddr >>  7) & 0x1f;
  chip         = (hwaddr >>  4) & 0x3;
  channel      = hwaddr & 0xf;
  return ret;
}

//____________________________________________________________________
Int_t
AliFMDAltroReader::ReadChannel(UShort_t& hwaddr, UShort_t& last, 
			       UShort_t* data) 
{
  Int_t ret, tmp;
  AliDebug(15, Form("Reading a channel"));
  if ((ret = ExtractTrailer(hwaddr, last)) < 0) { 
    AliError(Form("Failed to read trailer: %s", ErrorString(-ret)));
    return ret;
  }
  AliDebug(15, Form("Now extracting bunches from %d 10 bit words", last));
  tmp     =  ExtractBunches(last, data); 
  if (tmp < 0) {
    AliError(Form("Failed to read bunches: %s", ErrorString(-tmp)));
    return tmp;
  }
  ret     += tmp;
  last    =  (last == 0 ? 0 : last - 2); 
  return ret;
}

//____________________________________________________________________
Int_t
AliFMDAltroReader::ExtractTrailer(UShort_t& hwaddr, UShort_t& last)
{
  AliDebug(15, "Extracting trailer");
  W40_t trailer = GetNextW40();
  if (trailer < 0) {
    AliError(Form("Trailer 0x%x is bad: %s", trailer, ErrorString(-trailer)));
    return trailer;
  }
  if (!IsTrailer(trailer)) { 
    AliError(Form("Bad trailer: 0x%08x", trailer));
    return -kBadTrailer;
  }
  last    = (trailer >> 16) & 0x3ff;
  hwaddr  = (trailer & 0xfff);
  return 4;
}

//____________________________________________________________________
Int_t
AliFMDAltroReader::ExtractBunches(UShort_t last, UShort_t* data) 
{
  Int_t ret;
  if ((ret = ExtractFillWords(last)) < 0) { 
    AliError(Form("Failed to read fill words: %s", ErrorString(-ret)));
    return ret;
  }
  while (last > 0) { 
    Int_t tmp = ExtractBunch(data);
    if (tmp <= 0) { 
      AliError(Form("Failed to extract bunch at %d: %s", 
		    last, ErrorString(-tmp)));
      return tmp;
    }
    ret  += tmp;
    last -= tmp;
  }
  return ret;
}

//____________________________________________________________________
Int_t
AliFMDAltroReader::ExtractFillWords(UShort_t last) 
{
  // Number of fill words 
  UShort_t nFill = (last % 4 == 0 ? 0 : 4 - last % 4);
  // Read the fill words 
  for (UShort_t i = 3; i >= 4 - nFill; i--) {
    W10_t f = GetNextW10();
    if (f != 0x2aa) return -kBadFill;
  }
  return nFill;
}

//____________________________________________________________________
Int_t
AliFMDAltroReader::ExtractBunch(UShort_t* data)
{
  Int_t ret =  0;
  W10_t l =  GetNextW10(); 
  if (l < 0) { 
    AliError(Form("Failed to read bunch length: %s", ErrorString(-l)));
    return l;
  }
  W10_t t =  GetNextW10(); 
  if (t < 0) { 
    AliError(Form("Failed to read bunch time: %s", ErrorString(-t)));
    return t;
  }
  ret     += 2;
  for (Int_t i = 2; i < l; i++) {
    W10_t s = GetNextW10();
    if (s < 0) { 
      AliError(Form("Failed to read bunch data: %s", ErrorString(-s)));
      return 2;
    }
    AliDebug(50,Form("Assigning to data[%d - (%d - 1)] = 0x%X", t, i, s));
    data[t - (i-1)] = s;
    ret++;
  }
  return ret;
}

//____________________________________________________________________
Bool_t
AliFMDAltroReader::IsTrailer(W40_t x) 
{
  return ((x & fgkTrailerMask) == fgkTrailerMask);
}

//____________________________________________________________________
Bool_t
AliFMDAltroReader::IsBof() 
{
  return fCurrent == fBegin;
}

//____________________________________________________________________
Int_t
AliFMDAltroReader::ReadW40() 
{
  fInput.seekg(fCurrent-std::istream::pos_type(5));
  if (fInput.bad()) return -kBadSeek;
  fCurrent = fInput.tellg();
  if (fInput.bad()) return -kBadTell;
  fInput.read((char*)&fBuffer, 5 * sizeof(char));
  if (fInput.bad()) return -kBadRead;
  fIBuffer = 4;
  AliDebug(15, Form("  0x%03x  0x%03x  0x%03x  0x%03x    0x%010x  %6d", 
		    ExtractW10(3, fBuffer), ExtractW10(2, fBuffer), 
		    ExtractW10(1, fBuffer), ExtractW10(0, fBuffer), 
		    fBuffer, fCurrent));
  return fCurrent;
}

//____________________________________________________________________
AliFMDAltroIO::W10_t
AliFMDAltroReader::GetNextW10()
{
  if (fIBuffer <= 0) {
    Int_t ret;
    if ((ret = ReadW40()) < 0) return ret;
  }
  fIBuffer--;
  W10_t w10 = ExtractW10(fIBuffer, fBuffer); 
  return w10;
}

//____________________________________________________________________
AliFMDAltroIO::W40_t
AliFMDAltroReader::GetNextW40() 
{
  W40_t w40 = 0;
  for (Int_t i = 3; i >= 0; i--) {
    W10_t tmp  =  GetNextW10();
    W40_t bits =  ConcatW40(i, tmp);
    if (bits < 0) return bits;
    w40        += bits;
  }
  return w40;
}

//====================================================================
ClassImp(AliFMDAltroWriter)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDAltroWriter::AliFMDAltroWriter(std::ostream& stream) 
  : fThreshold(0), fTotal(0), fOutput(stream)
{
  AliDebug(15, "New AliFMDAltroWriter object");
  fTime   = 0;
  fLength = 0;
  fLast   = 0;
  // Write a dummy header
  fHeader = fOutput.tellp();
  if (fOutput.bad()) throw -kBadTell;
  AliRawDataHeader header;
  fOutput.write((char*)(&header), sizeof(header));
  if (fOutput.bad()) throw -kBadWrite;
  fBegin = fOutput.tellp();
  if (fOutput.bad()) throw -kBadTell;
}

//____________________________________________________________________
Int_t
AliFMDAltroWriter::Flush() 
{
  if (fIBuffer == 0) return 0;
  fOutput.write((char*)&fBuffer, 5 * sizeof(char));
  if (fOutput.bad()) return -kBadWrite;
  // for (UShort_t i = 0; i < 4; i++) 
  //   std::cout << "\t" << PRETTY_HEX(3, ExtractW10(i, fBuffer));
  // std::cout << "\t" << PRETTY_HEX(10, fBuffer) << std::endl;
  fTotal   += 5;
  fIBuffer =  0;
  fBuffer  =  0;
  return 5;
}

//____________________________________________________________________
Int_t 
AliFMDAltroWriter::Close() 
{
  Flush();
  std::ostream::pos_type end = fOutput.tellp();
  if (fOutput.bad()) return -kBadTell;
  fOutput.seekp(fHeader, std::ios_base::beg);
  if (fOutput.bad()) return -kBadSeek;
  AliRawDataHeader header;
  header.fSize = (UShort_t(end) - fHeader);
  AliDebug(15, Form("Size set to %d (%d)", header.fSize, fTotal));
  header.SetAttribute(0);
  fOutput.write((char*)(&header), sizeof(header));
  if (fOutput.bad()) return -kBadWrite;
  fOutput.seekp(end);
  if (fOutput.bad()) return -kBadSeek;
  return sizeof(header);
}


//____________________________________________________________________
Int_t
AliFMDAltroWriter::AddSignal(UShort_t adc) 
{
  Int_t ret = 0;
  if (adc < fThreshold) 
    ret = AddBunchTrailer();
  else {
    ret = AddToBuffer(adc);
    fLength++;
  }
  fTime++;
  if (ret < 0) AliError(Form("Failed to add signal %x: %s", ErrorString(ret)));
  return ret;
}

//____________________________________________________________________
Int_t
AliFMDAltroWriter::AddChannelTrailer(UShort_t board, UShort_t chip, 
				     UShort_t channel)
{
  UInt_t hwaddr = (channel & 0xf)+((chip & 0x3) << 4)+((board & 0x1f) << 7);
  return AddChannelTrailer(hwaddr);
}

//____________________________________________________________________
Int_t
AliFMDAltroWriter::AddChannelTrailer(UInt_t hwaddr)
{
  Int_t ret =0, tmp;
  if ((tmp = AddBunchTrailer()) < 0) { 
    AliError(Form("Failed to bad bunch trailer: %s", ErrorString(tmp)));
    return tmp;
  }
  ret += tmp;
  if ((tmp = AddFillWords())    < 0) { 
    AliError(Form("Failed to bad fill words: %s", ErrorString(tmp)));
    return tmp;
  }
  ret += tmp;
  W40_t trailer = (fgkTrailerMask + hwaddr + ((fLast & 0x3ff) << 16));
  fBuffer = trailer;
  fIBuffer = 3;
  ret     += 4;
  if ((tmp = Flush()) < 0) {
    AliError(Form("Failed to flush: %s", ErrorString(tmp)));
    return tmp;
  }
  ret     += tmp;
  fTime   =  0;
  fLast   =  0;
  return ret;
}

//____________________________________________________________________
Int_t
AliFMDAltroWriter::AddToBuffer(UShort_t x) 
{
  W40_t tmp = ConcatW40(fIBuffer, x);
  if (tmp < 0) return tmp;
  fBuffer += tmp;
  fIBuffer++;
  fLast++;
  Int_t ret = 0;
  if (fIBuffer > 3 && (ret = Flush() < 0)) return ret;
  return 1;
}

//____________________________________________________________________
Int_t
AliFMDAltroWriter::AddBunchTrailer()
{
  if (fLength <= 0) return 0;    
  Int_t ret = 0, tmp;
  if ((tmp = AddToBuffer(fTime))     < 0) return tmp;
  ret += tmp;
  if ((tmp = AddToBuffer(fLength+2)) < 0) return tmp;
  ret += tmp;
  fLength = 0;
  return ret;
}

//____________________________________________________________________
Int_t
AliFMDAltroWriter::AddFillWords() 
{
  Int_t ret = 0, tmp;
  if (fIBuffer == 0) return ret;
  for (Int_t i = fIBuffer; i < 4; i++) { 
    if ((tmp = AddToBuffer(0x2aa)) < 0) return tmp;
    ret += tmp;
    fLast--; 
  }
  if ((tmp = Flush() < 0)) return tmp;
  return ret;
}

//_____________________________________________________________________________
//
// EOF
//
