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
/** @file    AliFMDAltroMapping.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:27:56 2006
    @brief   Map HW to detector 
*/
//____________________________________________________________________
//                                                                          
// Mapping of ALTRO hardware channel to detector coordinates 
//
// The hardware address consist of a DDL number and 12bits of ALTRO
// addresses.  The ALTRO address are formatted as follows. 
//
//    12              7         4            0
//    |---------------|---------|------------|
//    | Board #       | ALTRO # | Channel #  |
//    +---------------+---------+------------+
//
// The mapping is done purely by calculations.  In the future,
// however, we may need some hard-coded stuff, or an external file to
// read from.  
//
#include "AliFMDAltroMapping.h"		// ALIFMDALTROMAPPING_H
#include "AliFMDParameters.h"
#include "AliLog.h"
#include "AliFMDDebug.h"
#include <iostream>
#include <iomanip>

//____________________________________________________________________
ClassImp(AliFMDAltroMapping)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//_____________________________________________________________________________
AliFMDAltroMapping::AliFMDAltroMapping()
{
  // Constructor 
}


//_____________________________________________________________________________
Bool_t
AliFMDAltroMapping::ReadMapping()
{
  // Read map from file - not used
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t
AliFMDAltroMapping::CreateInvMapping()
{
  // Create inverse mapping - not used
  return kTRUE;
}


//____________________________________________________________________
Bool_t 
AliFMDAltroMapping::Channel2StripBase(UShort_t  board, UShort_t  altro, 
				      UShort_t  chan,  Char_t&   ring, 
				      UShort_t& sec,   Short_t&  str) const
{
  // Translate a hardware address to detector coordinates. 
  // The detector is simply 
  // 
  //    ddl + 1
  // 
  // The ring number, sector, and strip number is given by the addr
  // argument.  The address argument, has the following format 
  // 
  //   12            7          4          0
  //   +-------------+----------+----------+
  //   | Board       | ALTRO    | Channel  |
  //   +-------------+----------+----------+
  // 
  // The board number identifier among other things the ring.  There's
  // up to 4 boards per DDL, and the two first (0 and 16) corresponds
  // to the inner rings, while the two last (1 and 17) corresponds to
  // the outer rings. 
  // 
  // The board number and ALTRO number together identifies the sensor,
  // and hence.  The lower board number (0 or 16) are the first N / 2
  // sensors (where N is the number of sensors in the ring).  
  // 
  // There are 3 ALTRO's per card, and each ALTRO serves up to 4
  // sensors.  Which of sensor is determined by the channel number.
  // For the inner rings, the map is
  // 
  //    ALTRO 0, Channel  0 to  7   -> Sensor 0 or 5
  //    ALTRO 0, Channel  8 to 15   -> Sensor 1 or 6
  //    ALTRO 1, Channel  0 to  7   -> Sensor 2 or 7 
  //    ALTRO 2, Channel  0 to  7   -> Sensor 3 or 8
  //    ALTRO 2, Channel  8 to 15   -> Sensor 4 or 9
  // 
  // For the outer rings, the map is 
  // 
  //    ALTRO 0, Channel  0 to  3   -> Sensor 0 or 10
  //    ALTRO 0, Channel  4 to  7   -> Sensor 1 or 11
  //    ALTRO 0, Channel  8 to 11   -> Sensor 2 or 12
  //    ALTRO 0, Channel 12 to 15   -> Sensor 3 or 13
  //    ALTRO 1, Channel  0 to  3   -> Sensor 4 or 14
  //    ALTRO 1, Channel  4 to  7   -> Sensor 5 or 15
  //    ALTRO 2, Channel  0 to  3   -> Sensor 6 or 16
  //    ALTRO 2, Channel  4 to  7   -> Sensor 7 or 17
  //    ALTRO 2, Channel  8 to 11   -> Sensor 8 or 18
  //    ALTRO 2, Channel 12 to 15   -> Sensor 9 or 19
  // 
  // Which divison of the sensor we're in, depends on the channel
  // number only.  For the inner rings, the map is 
  // 
  //    Channel 0                   -> Sector 0, strips   0-127
  //    Channel 1                   -> Sector 1, strips 127-0
  //    Channel 3                   -> Sector 0, strips 128-255
  //    Channel 4                   -> Sector 1, strips 255-128
  //    Channel 5                   -> Sector 0, strips 256-383
  //    Channel 6                   -> Sector 1, strips 383-256
  //    Channel 7                   -> Sector 0, strips 384-511
  //    Channel 8                   -> Sector 1, strips 511-384
  // 
  // There are only half as many strips in the outer sensors, so there
  // only 4 channels are used for a full sensor.  The map is 
  // 
  //    Channel 0                   -> Sector 0, strips   0-127
  //    Channel 1                   -> Sector 1, strips 127-0
  //    Channel 3                   -> Sector 0, strips 128-255
  //    Channel 4                   -> Sector 1, strips 255-128
  // 
  // With this information, we can decode the hardware address to give
  // us detector coordinates, unique at least up a 128 strips.  We
  // return the first strip, as seen by the ALTRO channel, in the
  // given range.  
  //
  ring          =  Board2Ring(board);
  UShort_t fsec =  board < 16 ? 1 : 0;
  switch (ring) {
  case 'i':
  case 'I':
    sec = (fsec * 10 + (altro < 1 ? 0 : altro < 2 ? 4 : 6) 
	   + 2 * (chan / 8) + chan % 2);
    str = ((chan % 8) / 2) * 128;
    break;
  case 'o':
  case 'O': 
    sec = (fsec * 20 + (altro < 1 ? 0 : altro < 2 ? 8 : 12) 
	   + 2 * (chan / 4) + chan % 2);
    str = ((chan % 4) / 2) * 128;
    break;
  }
  if (sec % 2) str += 127;
  // AliFMDDebug(1, ("%02x/%x/%x Base strip = %d", board, altro, chan, str));
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDAltroMapping::Timebin2Strip(UShort_t  sec,         
				  UShort_t  timebin,
				  UShort_t  preSamples, 
				  UShort_t  sampleRate, 
				  Short_t&  stripOff,      
				  UShort_t& sample) const
{
  // Compute the strip off-set in the current channel from the sector
  // and timebin.   Also needed for this computation is the basic
  // offset in timebins, as well as the sample rat. 
  UShort_t t =  (timebin - preSamples);
  sample     =  (t % sampleRate);
  t          -= sample;
  stripOff   =  (sec % 2 ? -1 : 1) * t / sampleRate;
}

#if 0
//____________________________________________________________________
Bool_t 
AliFMDAltroMapping::Hardware2Detector(UInt_t    ddl,   UInt_t    board, 
				      UInt_t    altro, UInt_t    chan,
				      UShort_t& det,   Char_t&   ring, 
				      UShort_t& sec,   Short_t&  str) const
{
  // See also Hardware2Detector that requires 3 inputs 
  det = DDL2Detector(ddl);
  return Channel2StripBase(board, altro, chan, ring, sec, str);
}


//____________________________________________________________________
Bool_t 
AliFMDAltroMapping::Hardware2Detector(UInt_t    ddl, UInt_t    addr, 
				      UShort_t& det, Char_t&   ring, 
				      UShort_t& sec, Short_t&  str) const
{
  // Translate a hardware address to detector coordinates. 
  // 
  // See also Hardware2Detector that accepts 4 inputs 
  UShort_t board, altro, chan;
  ChannelAddress(addr, board, altro, chan);
  return Hardware2Detector(ddl,board, altro, chan, det,ring, sec, str);
}
#endif 

//____________________________________________________________________
Bool_t 
AliFMDAltroMapping::Hardware2Detector(UShort_t  ddl,     UShort_t    board, 
				      UShort_t  altro,   UShort_t    chan,
				      UShort_t  timebin, UShort_t    preSamples,
				      UShort_t  sampleRate,
				      UShort_t& det,     Char_t&   ring, 
				      UShort_t& sec,     Short_t&  str,
				      UShort_t& sam) const
{
  // Full conversion from hardware address, including timebin number,
  // to detector coordinates and sample number.  Note, that this
  // conversion depends on the oversampling rate and the number of
  // pre-samples 
  Short_t  baseStrip, stripOffset, tdet  = DDL2Detector(ddl);
  if (tdet < 0) return kFALSE;
  det = tdet;
  if (!Channel2StripBase(board, altro, chan, ring, sec, baseStrip)) 
    return kFALSE;
  Timebin2Strip(sec, timebin, preSamples, sampleRate, stripOffset, sam);
#if 0
  AliFMDDebug(1, ("0x%x/0x%02x/0x%x/0x%x/%04d -> FMD%d%c[%2d,%3d]-%d "
		  "(pre=%d,rate=%d,base=%3d,off=%3d)", 
		  ddl, 
		  board, 
		  altro, 
		  chan, 
		  timebin, 
		  det, 
		  ring, 
		  sec, 
		  str, 
		  sam,
		  preSamples, 
		  sampleRate,
		  baseStrip, 
		  stripOffset));
#endif
  str = baseStrip + stripOffset;
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDAltroMapping::Hardware2Detector(UShort_t  ddl,       UShort_t addr,
				      UShort_t  timebin,   UShort_t preSamples, 
				      UShort_t  sampleRate,
				      UShort_t& det,       Char_t&   ring, 
				      UShort_t& sec,       Short_t&  str,
				      UShort_t& sam) const
{
  // Translate a hardware address to detector coordinates. 
  // 
  // See also Hardware2Detector that accepts 4 inputs 
  UShort_t board, altro, chan;
  ChannelAddress(addr, board, altro, chan);
  return Hardware2Detector(ddl, board, altro, chan, 
			   timebin, preSamples, sampleRate, 
			   det, ring, sec, str, sam);
}


//____________________________________________________________________
Short_t
AliFMDAltroMapping::Sector2Board(Char_t ring, UShort_t sec) const
{
  switch (ring) { 
  case 'I': 
  case 'i':
    return (sec < 10 ? 16 : 0); // (sec / 10) * 16;
  case 'O': 
  case 'o': 
    return (sec < 20 ? 16 : 0) + 1; // (sec / 20) * 16 + 1;
  }
  return -1;
}

//_____________________________________________ _______________________
Bool_t
AliFMDAltroMapping::Strip2Channel(Char_t    ring,  UShort_t  sec,   
				  UShort_t  str,   UShort_t& board,
				  UShort_t& altro, UShort_t& chan) const
{
  // Translate detector coordinates to a hardware address.
  // The ddl is simply 
  // 
  //    (det - 1)
  // 
  // The ring number, sector, and strip number must be encoded into a
  // hardware address.  The address argument, will have the following
  // format on output
  // 
  //   12            7          4          0
  //   +-------------+----------+----------+
  //   | Board       | ALTRO    | Channel  |
  //   +-------------+----------+----------+
  // 
  // The board number is given by the ring and sector.  The inner
  // rings board 0 and 16, while the outer are 1 and 17.  Which of these
  // depends on the sector.  The map is 
  // 
  //    Ring I, sector  0- 9       ->   board 0 
  //    Ring I, sector 10-19       ->   board 16
  //    Ring O, sector  0-19       ->   board 1 
  //    Ring O, sector 20-39       ->   board 17
  // 
  // There are 3 ALTRO's per board.  The ALTRO number is given by the
  // sector number.  For the inner rings, these are given by
  // 
  //    Sector  0- 3  or  10-13    ->   ALTRO 0 
  //    Sector  4- 5  or  14-15    ->   ALTRO 1 
  //    Sector  6- 9  or  16-19    ->   ALTRO 2 
  // 
  // For the outers, it's given by 
  // 
  //    Sector  0- 7  or  20-27    ->   ALTRO 0 
  //    Sector  8-11  or  28-31    ->   ALTRO 1 
  //    Sector 12-19  or  32-39    ->   ALTRO 2 
  //
  // The channel number is given by the sector and strip number.  For
  // the inners, the map is 
  // 
  //    Sector  0, strips   0-127  ->   Channel  0
  //    Sector  0, strips 128-255  ->   Channel  2
  //    Sector  0, strips 256-383  ->   Channel  4
  //    Sector  0, strips 384-511  ->   Channel  6
  //    Sector  1, strips 127-  0  ->   Channel  1
  //    Sector  1, strips 255-128  ->   Channel  3
  //    Sector  1, strips 383-256  ->   Channel  5
  //    Sector  1, strips 511-384  ->   Channel  7
  //    Sector  2, strips   0-127  ->   Channel  8
  //    Sector  2, strips 128-255  ->   Channel 10
  //    Sector  2, strips 256-383  ->   Channel 12
  //    Sector  2, strips 384-511  ->   Channel 14
  //    Sector  3, strips 127-  0  ->   Channel  9
  //    Sector  3, strips 255-128  ->   Channel 11
  //    Sector  3, strips 383-256  ->   Channel 13
  //    Sector  3, strips 511-384  ->   Channel 15
  // 
  // and so on, up to sector 19.  For the outer, the map is 
  // 
  //    Sector  0, strips   0-127  ->   Channel  0
  //    Sector  0, strips 128-255  ->   Channel  2
  //    Sector  1, strips 127-  0  ->   Channel  1
  //    Sector  1, strips 255-128  ->   Channel  3
  //    Sector  2, strips   0-127  ->   Channel  4
  //    Sector  2, strips 128-255  ->   Channel  6
  //    Sector  3, strips 127-  0  ->   Channel  5
  //    Sector  3, strips 255-128  ->   Channel  7
  //    Sector  4, strips   0-127  ->   Channel  8
  //    Sector  4, strips 128-255  ->   Channel 10
  //    Sector  5, strips 127-  0  ->   Channel  9
  //    Sector  5, strips 255-128  ->   Channel 11
  //    Sector  6, strips   0-127  ->   Channel 12
  //    Sector  6, strips 128-255  ->   Channel 14
  //    Sector  7, strips 127-  0  ->   Channel 13
  //    Sector  7, strips 255-128  ->   Channel 15
  // 
  // and so on upto sector 40. 
  // 
  // With this information, we can decode the detector coordinates to
  // give us a unique hardware address 
  //
  UInt_t   tmp    = 0;
  UShort_t fboard = 0;
  switch (ring) {
  case 'I':
  case 'i':
    fboard =  sec < 10 ? 1 : 0;
    board  =  fboard * 16;
    altro  =  (sec % 10) < 4 ? 0 : (sec % 10) < 6 ? 1 : 2;
    tmp    =  (sec % 10) - (altro == 0 ? 0 : altro == 1 ? 4 : 6);
    chan   =  2  * (str / 128) + (sec % 2) + ((tmp / 2) % 2) * 8;
    break;
  case 'O':
  case 'o':
    fboard =  sec < 20 ? 1 : 0;
    board  =  fboard * 16 + 1;
    altro  =  (sec % 20) < 8 ? 0 : (sec % 20) < 12 ? 1 : 2;
    tmp    =  (sec % 20) - (altro == 0 ? 0 : altro == 1 ? 8 : 12);
    chan   =  2 * (str / 128) + (sec % 2) + ((tmp / 2) % 4) * 4;
    break;
  }
  return kTRUE;
}

//_____________________________________________ _______________________
UShort_t
AliFMDAltroMapping::Strip2Timebin(UShort_t sec, UShort_t strip, 
				  UShort_t sam, UShort_t preSamples, 
				  UShort_t sampleRate) const
{
  UShort_t timebin = preSamples;
  if (sec % 2)  timebin += (127 - (strip % 128)) * sampleRate;
  else          timebin += (strip % 128) * sampleRate;
  timebin += sam;
  return timebin;
}

#if 0
//_____________________________________________ _______________________
Bool_t 
AliFMDAltroMapping::Detector2Hardware(UShort_t  det,   Char_t    ring, 
				      UShort_t  sec,   UShort_t  str,
				      UShort_t& ddl,   UShort_t& board,
				      UShort_t& altro, UShort_t& chan) const
{
  ddl =  Detector2DDL(det);
  return Strip2Channel(ring, sec, str, board, altro, chan);
}


//____________________________________________________________________
Bool_t 
AliFMDAltroMapping::Detector2Hardware(UShort_t  det, Char_t    ring, 
				      UShort_t  sec, UShort_t  str,
				      UShort_t& ddl, UShort_t& addr) const
{
  // Translate detector coordinates to a hardware address.  
  // 
  // See also Detector2Hardware that returns 4 parameters.  
  UShort_t board = 0;
  UShort_t altro = 0;
  UShort_t chan  = 0;
  if (!Detector2Hardware(det,ring,sec,str,ddl,board,altro,chan)) return kFALSE;
  addr =  ChannelAddress(board, altro, chan);
  return kTRUE;
}
#endif

//____________________________________________________________________
Bool_t 
AliFMDAltroMapping::Detector2Hardware(UShort_t  det,        Char_t    ring, 
				      UShort_t  sec,        UShort_t  str,
				      UShort_t  sam, 
				      UShort_t  preSamples, 
				      UShort_t  sampleRate,
				      UShort_t& ddl,        UShort_t& board, 
				      UShort_t& altro,      UShort_t& channel, 
				      UShort_t& timebin) const
{
  ddl = Detector2DDL(det);
  if (!Strip2Channel(ring,sec,str,board,altro,channel)) return kFALSE;
  timebin = Strip2Timebin(sec, str, sam, preSamples, sampleRate);
  return kTRUE;
}


//____________________________________________________________________
Bool_t 
AliFMDAltroMapping::Detector2Hardware(UShort_t  det,        Char_t   ring, 
				      UShort_t  sec,        UShort_t str,
				      UShort_t  sam, 
				      UShort_t  preSamples, UShort_t sampleRate,
				      UShort_t& ddl,        UShort_t&  hwaddr, 
				      UShort_t& timebin) const
{
  UShort_t board = 0;
  UShort_t altro = 0;
  UShort_t chan  = 0;
  if (!Detector2Hardware(det, ring, sec, str, sam, 
			 preSamples, sampleRate, 
			 ddl, board, altro, chan, timebin)) return kFALSE;
  hwaddr =  ChannelAddress(board, altro, chan);
  return kTRUE;
}

      
//____________________________________________________________________
UInt_t 
AliFMDAltroMapping::ChannelAddress(UShort_t board, UShort_t altro, 
				   UShort_t channel) const
{
  return (((board & 0x1F) << 7) | ((altro & 0x7) << 4) | (channel & 0xF));
}

//____________________________________________________________________
void 
AliFMDAltroMapping::ChannelAddress(UShort_t hwaddr, UShort_t& board, 
				   UShort_t& altro, UShort_t& channel) const
{
  board   = ((hwaddr >> 7) & 0x1F);
  altro   = ((hwaddr >> 4) & 0x07);
  channel = ((hwaddr >> 0) & 0x0F);
}

//____________________________________________________________________
Int_t
AliFMDAltroMapping::GetHWAddress(Int_t sec, Int_t str, Int_t ring)
{
  // Return hardware address corresponding to sector sec, strip str,
  // and ring ring.   Mapping from TPC to FMD coordinates are 
  // 
  //   TPC     | FMD
  //   --------+------
  //   padrow  | sector
  //   pad     | strip
  //   sector  | ring 
  //
  Char_t r = Char_t(ring);
  UShort_t board, altro, channel;
  Strip2Channel(r, sec, str, board, altro, channel);
  return ChannelAddress(board, altro, channel);
}

//____________________________________________________________________
Int_t
AliFMDAltroMapping::GetPadRow(Int_t hwaddr) const
{
  // Return sector corresponding to hardware address hwaddr.  Mapping
  // from TPC to FMD coordinates are  
  // 
  //   TPC     | FMD
  //   --------+------
  //   padrow  | sector
  //   pad     | strip
  //   sector  | ring 
  //
  Char_t   ring;
  UShort_t board, altro, channel, sector;
  Short_t  baseStrip;
  ChannelAddress(hwaddr, board, altro, channel);
  if (!Channel2StripBase(board, altro, channel, ring, sector, baseStrip))
    return -1;
  return Int_t(sector);
}

//____________________________________________________________________
Int_t
AliFMDAltroMapping::GetPad(Int_t hwaddr) const
{
  // Return strip corresponding to hardware address hwaddr.  Mapping
  // from TPC to FMD coordinates are  
  // 
  //   TPC     | FMD
  //   --------+------
  //   padrow  | sector
  //   pad     | strip
  //   sector  | ring 
  //
  Char_t   ring;
  UShort_t board, altro, channel, sector;
  Short_t  baseStrip;
  ChannelAddress(hwaddr, board, altro, channel);
  if (!Channel2StripBase(board, altro, channel, ring, sector, baseStrip))
    return -1;
  return Int_t(baseStrip);
}
  
//____________________________________________________________________
Int_t
AliFMDAltroMapping::GetSector(Int_t hwaddr) const
{
  // Return ring corresponding to hardware address hwaddr.  Mapping
  // from TPC to FMD coordinates are  
  // 
  //   TPC     | FMD
  //   --------+------
  //   padrow  | sector
  //   pad     | strip
  //   sector  | ring 
  //
  Char_t   ring;
  UShort_t board, altro, channel, sector;
  Short_t  baseStrip;
  ChannelAddress(hwaddr, board, altro, channel);
  if (!Channel2StripBase(board, altro, channel, ring, sector, baseStrip))
    return -1;
  return Int_t(ring);
}

//____________________________________________________________________
void
AliFMDAltroMapping::Print(Option_t* option) const
{
  TString opt(option);
  opt.ToLower();
  UShort_t ddl, board, chip, chan, addr;
  UShort_t det, sec;
  Short_t  strBase;
  Char_t   rng;
  
  if (opt.Contains("hw") || opt.Contains("hardware")) { 
    std::cout << " DDL | Board | Chip | Chan | Address | Detector\n"
	      << "=====+=======+======+======+=========+===============" 
	      << std::endl;
    for (ddl = 0; ddl <= 2; ddl++) { 
      Int_t  boards[] = { 0, 16, (ddl == 0 ? 32 : 1), 17, 32};
      Int_t* ptr      = boards;
      det             = DDL2Detector(ddl);
      while ((board = *(ptr++)) < 32) { 
	for (chip = 0; chip <= 2; chip++) { 
	  UShort_t nchan = (chip == 1 ? 8 : 16);
	  for (chan = 0; chan < nchan; chan++) { 
	    Channel2StripBase(board, chip, chan, rng, sec, strBase);
	    addr = ChannelAddress(board, chip, chan);
	    std::cout << " "  
		      << std::setw(3) << ddl     << " | " 
		      << std::setfill('0')       << std::hex << " 0x"
		      << std::setw(2) << board   << " |  0x"
		      << std::setw(1) << chip    << " |  0x"
		      << std::setw(1) << chan    << " |   0x" 
		      << std::setw(3) << addr    << " | " 
		      << std::setfill(' ')       << std::dec << " FMD" 
		      << std::setw(1) << det     << rng << "[" 
		      << std::setw(2) << sec     << "," 
		      << std::setw(3) << strBase << "]" << std::endl;
	  } // for chan ...
	  if (chip == 2 && *ptr >= 32) continue;
	  std::cout << "     +       +      +      +         +              " 
		    << std::endl;
	} // for chip ... 
      } // while board 
      std::cout << "-----+-------+------+------+---------+---------------" 
		<< std::endl;
    } // for ddl ... 
  } // if hw 
  if (opt.Contains("det")) { 
    std::cout << " Detector      | DDL | Board | Chip | Chan | Address\n"
	      << "===============+=====+=======+======+======+========"
	      << std::endl;
    for (det = 1; det <= 3; det++) { 
      Char_t  rings[] = { 'I', (det == 1 ? '\0' : 'O'),'\0' };
      Char_t* ptr     = rings;
      ddl             = Detector2DDL(det);
      while ((rng = *(ptr++)) != '\0') { 
	UShort_t nsec = (rng == 'I' ?  20 :  40);
	UShort_t nstr = (rng == 'I' ? 512 : 256);
	for (sec = 0; sec < nsec; sec++) { 
	  for (strBase = 0; strBase < nstr; strBase += 128) {
	    Strip2Channel(rng, sec, strBase, board, chip, chan);
	    addr = ChannelAddress(board, chip, chan);
	    std::cout << std::setfill(' ')         << std::dec << " FMD" 
		      << std::setw(1) << det       << rng      << "[" 
		      << std::setw(2) << sec       << "," 
		      << std::setw(3) << strBase   << "] | " 
		      << std::setw(3) << ddl       << " |  0x"
		      << std::setfill('0')         << std::hex  
		      << std::setw(2) << board     << " |  0x"
		      << std::setw(1) << chip      << " |  0x"
		      << std::setw(1) << chan      << " |   0x" 
		      << std::setw(3) << addr      << std::endl;
	  } // for str ...
	} // for sec ... 
	if (*ptr == '\0') continue;
	std::cout << "               +     +       +      +      +        " 
		  << std::endl;
      } // while rng ... 
      std::cout << "---------------+-----+-------+------+------+--------" 
		<< std::endl;

    } // for det ... 
  } // if det 
}

//_____________________________________________________________________________
//
// EOF
//
