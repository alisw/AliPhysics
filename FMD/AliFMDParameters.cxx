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
// Forward Multiplicity Detector based on Silicon wafers. 
//
// This class is a singleton that handles various parameters of
// the FMD detectors.  
// Eventually, this class will use the Conditions DB to get the
// various parameters, which code can then request from here.
//                                                       
#include "AliLog.h"		   // ALILOG_H
#include "AliFMDParameters.h"	   // ALIFMDPARAMETERS_H
#include "AliFMDGeometry.h"	   // ALIFMDGEOMETRY_H
#include "AliFMDRing.h"	           // ALIFMDRING_H
#include "AliFMDCalibGain.h"       // ALIFMDCALIBGAIN_H
#include "AliFMDCalibPedestal.h"   // ALIFMDCALIBPEDESTAL_H
#include "AliFMDCalibSampleRate.h" // ALIFMDCALIBPEDESTAL_H
#include <Riostream.h>

//====================================================================
ClassImp(AliFMDParameters)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDParameters* AliFMDParameters::fgInstance = 0;

//____________________________________________________________________
AliFMDParameters* 
AliFMDParameters::Instance() 
{
  // Get static instance 
  if (!fgInstance) fgInstance = new AliFMDParameters;
  return fgInstance;
}

//____________________________________________________________________
AliFMDParameters::AliFMDParameters() 
  : fSiDeDxMip(1.664), 
    fFixedPulseGain(0), 
    fEdepMip(0),
    fZeroSuppression(0), 
    fSampleRate(0), 
    fPedestal(0), 
    fPulseGain(0), 
    fDeadMap(0)
{
  // Default constructor 
  SetVA1MipRange();
  SetAltroChannelSize();
  SetChannelsPerAltro();
  SetZeroSuppression();
  SetSampleRate();
  SetPedestal();
  SetPedestalWidth();
  SetPedestalFactor();
  SetThreshold();
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetThreshold() const
{
  if (!fPulseGain) return fFixedThreshold;
  return fPulseGain->Threshold();
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetPulseGain(UShort_t detector, Char_t ring, 
			       UShort_t sector, UShort_t strip) const
{
  // Returns the pulser calibrated gain for strip # strip in sector #
  // sector or ring id ring of detector # detector. 
  // 
  // For simulation, this is normally set to 
  // 
  //       VA1_MIP_Range 
  //    ------------------ * MIP_Energy_Loss
  //    ALTRO_channel_size
  // 
  if (!fPulseGain) { 
    if (fFixedPulseGain <= 0)
      fFixedPulseGain = fVA1MipRange * GetEdepMip() / fAltroChannelSize;
    return fFixedPulseGain;
  }  
  return fPulseGain->Value(detector, ring, sector, strip);
}

//__________________________________________________________________
Bool_t
AliFMDParameters::IsDead(UShort_t detector, Char_t ring, 
			 UShort_t sector, UShort_t strip) const
{
  if (!fDeadMap) return kFALSE;
  return fDeadMap->operator()(detector, ring, sector, strip);
}

//__________________________________________________________________
UShort_t
AliFMDParameters::GetZeroSuppression(UShort_t detector, Char_t ring, 
				     UShort_t sector, UShort_t strip) const
{
  if (!fZeroSuppression) return fFixedZeroSuppression;
  // Need to map strip to ALTRO chip. 
  return fZeroSuppression->operator()(detector, ring, sector, strip/128);
}

//__________________________________________________________________
UShort_t
AliFMDParameters::GetSampleRate(UShort_t ddl) const
{
  if (!fSampleRate) return fFixedSampleRate;
  // Need to map sector to digitizier card. 
  return fSampleRate->Rate(ddl);
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetPedestal(UShort_t detector, Char_t ring, 
			      UShort_t sector, UShort_t strip) const
{
  if (!fPedestal) return fFixedPedestal;
  return fPedestal->Value(detector, ring, sector, strip);
}

//__________________________________________________________________
Float_t
AliFMDParameters::GetPedestalWidth(UShort_t detector, Char_t ring, 
				   UShort_t sector, UShort_t strip) const
{
  if (!fPedestal) return fFixedPedestalWidth;
  return fPedestal->Width(detector, ring, sector, strip);
}
  
			      

//__________________________________________________________________
Float_t
AliFMDParameters::GetEdepMip() const 
{ 
  // Get energy deposited by a MIP in the silicon sensors
  if (fEdepMip <= 0){
    AliFMDGeometry* fmd = AliFMDGeometry::Instance();
    fEdepMip = (fSiDeDxMip 
		* fmd->GetRing('I')->GetSiThickness() 
		* fmd->GetSiDensity());
  }
  return fEdepMip;
}

//__________________________________________________________________
Bool_t
AliFMDParameters::Hardware2Detector(UInt_t ddl, UInt_t addr, UShort_t& det,
				    Char_t& ring, UShort_t& sec, 
				    UShort_t& str) const
{
  // Translate a hardware address to detector coordinates. 
  // The detector is simply 
  // 
  //    ddl - kBaseDDL + 1
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
  // up to 4 boards per DDL, and the two first (0 and 1) corresponds
  // to the inner rings, while the two last (2 and 3) corresponds to
  // the outer rings. 
  // 
  // The board number and ALTRO number together identifies the sensor,
  // and hence.  The lower board number (0 or 2) are the first N / 2
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
  //    Channel 1                   -> Sector 1, strips   0-127
  //    Channel 3                   -> Sector 0, strips 128-255
  //    Channel 4                   -> Sector 1, strips 128-255
  //    Channel 5                   -> Sector 0, strips 256-383
  //    Channel 6                   -> Sector 1, strips 256-383
  //    Channel 7                   -> Sector 0, strips 384-511
  //    Channel 8                   -> Sector 1, strips 384-511 
  // 
  // There are only half as many strips in the outer sensors, so there
  // only 4 channels are used for a full sensor.  The map is 
  // 
  //    Channel 0                   -> Sector 0, strips   0-127
  //    Channel 1                   -> Sector 1, strips   0-127
  //    Channel 3                   -> Sector 0, strips 128-255
  //    Channel 4                   -> Sector 1, strips 128-255
  // 
  // With this information, we can decode the hardware address to give
  // us detector coordinates, unique at least up a 128 strips.  We
  // return the first strip in the given range. 
  //
  det          =  (ddl - kBaseDDL) + 1;
  UInt_t board =  (addr >> 7) & 0x1F;
  UInt_t altro =  (addr >> 4) & 0x7;
  UInt_t chan  =  (addr & 0xf);
  if (board > 3) {
    AliError(Form("Invalid board address %d for the FMD", board));
    return kFALSE;
  }
  if (altro > 2) {
    AliError(Form("Invalid ALTRO address %d for the FMD digitizer %d", 
		  altro, board));
    return kFALSE;
  }
  ring         =  (board > 1 ? 'O' : 'I');
  UInt_t nsen  =  (ring == 'I' ? 10 : 20);
  UInt_t nsa   =  (ring == 'I' ? 2 : 4);   // Sensors per ALTRO
  UInt_t ncs   =  (ring == 'I' ? 8 : 4);   // Channels per sensor 
  UInt_t sen   =  (board % 2) * nsen / 2;  // Base for half-ring
  sen          += chan / ncs + (altro == 0 ? 0   : 
				altro == 1 ? nsa : UInt_t(1.5 * nsa));
  sec          =  2 * sen + (chan % 2);
  str          =  (chan % ncs) / 2 * 128;
  return kTRUE;
}

//__________________________________________________________________
Bool_t
AliFMDParameters::Detector2Hardware(UShort_t det, Char_t ring, UShort_t sec, 
				    UShort_t str, UInt_t& ddl, UInt_t& addr) 
  const
{
  // Translate detector coordinates to a hardware address.
  // The ddl is simply 
  // 
  //    kBaseDDL + (det - 1)
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
  // rings board 0 and 1, while the outer are 2 and 3.  Which of these
  // depends on the sector.  The map is 
  // 
  //    Ring I, sector  0- 9       ->   board 0 
  //    Ring I, sector 10-19       ->   board 1
  //    Ring O, sector  0-19       ->   board 2 
  //    Ring O, sector 20-39       ->   board 3
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
  //    Sector  1, strips   0-127  ->   Channel  1
  //    Sector  1, strips 128-255  ->   Channel  3
  //    Sector  1, strips 256-383  ->   Channel  5
  //    Sector  1, strips 384-511  ->   Channel  7
  //    Sector  2, strips   0-127  ->   Channel  8
  //    Sector  2, strips 128-255  ->   Channel 10
  //    Sector  2, strips 256-383  ->   Channel 12
  //    Sector  2, strips 384-511  ->   Channel 14
  //    Sector  3, strips   0-127  ->   Channel  9
  //    Sector  3, strips 128-255  ->   Channel 11
  //    Sector  3, strips 256-383  ->   Channel 13
  //    Sector  3, strips 384-511  ->   Channel 15
  // 
  // and so on, up to sector 19.  For the outer, the map is 
  // 
  //    Sector  0, strips   0-127  ->   Channel  0
  //    Sector  0, strips 128-255  ->   Channel  2
  //    Sector  1, strips   0-127  ->   Channel  1
  //    Sector  1, strips 128-255  ->   Channel  3
  //    Sector  2, strips   0-127  ->   Channel  4
  //    Sector  2, strips 128-255  ->   Channel  6
  //    Sector  3, strips   0-127  ->   Channel  5
  //    Sector  3, strips 128-255  ->   Channel  7
  //    Sector  4, strips   0-127  ->   Channel  8
  //    Sector  4, strips 128-255  ->   Channel 10
  //    Sector  5, strips   0-127  ->   Channel  9
  //    Sector  5, strips 128-255  ->   Channel 11
  //    Sector  6, strips   0-127  ->   Channel 12
  //    Sector  6, strips 128-255  ->   Channel 14
  //    Sector  7, strips   0-127  ->   Channel 13
  //    Sector  7, strips 128-255  ->   Channel 15
  // 
  // and so on upto sector 40. 
  // 
  // With this information, we can decode the detector coordinates to
  // give us a unique hardware address 
  //
  ddl          =  kBaseDDL + (det - 1);
  UInt_t nsen  =  (ring == 'I' ? 10 : 20);
  UInt_t nsa   =  (ring == 'I' ? 2 : 4);   // Sensors per ALTRO
  UInt_t ncs   =  (ring == 'I' ? 8 : 4);   // Channels per sensor 
  UInt_t bbase =  (ring == 'I' ? 0 : 2);
  UInt_t board =  bbase + sec / nsen;
  UInt_t lsen  =  (sec - (board - bbase) * nsen);
  UInt_t altro =  (lsen < 2 * nsa ? 0 : (lsen < 3 * nsa ? 1 : 2));
  UInt_t sbase =  (lsen < 2 * nsa ? 0 : (lsen < 3 * nsa ? 2*nsa : 3*nsa));
  UInt_t chan  =  (sec % 2) + (lsen-sbase) / 2 * ncs + 2 * str / 128;
  AliDebug(40, Form("\n"
		    "  chan = (%d %% 2) + (%d-%d) / %d * %d + 2 * %d / 128\n"
		    "       = %d + %d + %d = %d", 
		    sec, lsen, sbase, 2, ncs, str, 
		    (sec % 2), (lsen - sbase) / 2 * ncs, 
		    2 * str / 128, chan));
  addr         =  chan + (altro << 4) + (board << 7);
  
  return kTRUE;
}

  
  
  
//____________________________________________________________________
//
// EOF
//
