#ifndef ALIEMCALGEOPARAMS_H
#define ALIEMCALGEOPARAMS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEMCALGeoParams.h $ */

//////////////////////////////////////////////////////////
// class for holding various parameters; 
//////////////////////////////////////////////////////////

class AliEMCALGeoParams
{
public:

  // general geometry info
  static const int fgkEMCALModules = 12; // number of modules for EMCAL
  static const int fgkEMCALRows = 24; // number of rows per module for EMCAL
  static const int fgkEMCALCols = 48; // number of columns per module for EMCAL

  static const int fgkEMCALLEDRefs = 24; // number of LEDs (reference/monitors) per module for EMCAL; one per StripModule
  static const int fgkEMCALTempSensors = 8; // number Temperature sensors per module for EMCAL

  Int_t GetStripModule(Int_t iSM, Int_t iCol) const
    // Strip 0 is the one closest to the FEE crates; different for A (iColumn/2) and C sides
    { return ( (iSM%2==0) ? iCol/2 : AliEMCALGeoParams::fgkEMCALLEDRefs - 1 - iCol/2 ); }

  // also a few readout related variables:
  static const int fgkLastAltroDDL = 43; // 0..23 (i.e. 24) for EMCAL; 24..43 (i.e. 20) allocated for DCAL 
  static const int fgkSampleMax = 1023; // highest possible sample value (10-bit = 0x3ff)
  static const int fgkOverflowCut = 950; // saturation starts around here; also exist as private constant in AliEMCALRawUtils, should probably be replaced
  static const int fgkSampleMin = 0; // lowest possible sample value 

  // TRU numbers
  static const int fgkEMCALTRUsPerSM = 3; // number of TRU's in a SuperModule
  static const int fgkEMCAL2x2PerTRU = 96; // number of 2x2's in a TRU
  static const int fgkEMCALTRURows   = 4;
  static const int fgkEMCALTRUCols   = 24;
  
  // RAW/AliCaloAltroMapping provides the correspondence information between
  // an electronics HWAddress (Branch<<1 | FEC<<7 | ALTRO<<4 | Channel) 
  // for the RCUs and which tower (Column and Row) that corresponds to. 
  // For the cases when one doesn't have a Raw stream to decode the HW address
  // into the other FEE indices, we provide the needed simple methods here 
  // with arguments (within an RCU)
  Int_t GetHWAddress(Int_t iBranch, Int_t iFEC, Int_t iALTRO, Int_t iChannel) const
  { return ( (iBranch<<11) | (iFEC<<7) | (iALTRO<<4) | iChannel ); }; // 
  // and for converting back to the individual indices
  Int_t GetBranch(Int_t iHW) const { return ( (iHW>>11) & 0x1 ); }; // 
  Int_t GetFEC(Int_t iHW) const { return ( (iHW>>7) & 0xf ); }; // 
  Int_t GetAltro(Int_t iHW) const { return ( (iHW>>4) & 0x7 ); }; // 
  Int_t GetChannel(Int_t iHW) const { return ( iHW & 0xf ); }; // 

  // We can also encode a very similar CSP address
  Int_t GetCSPAddress(Int_t iBranch, Int_t iFEC, Int_t iCSP) const
  { return ( (iBranch<<11) | (iFEC<<7) | iCSP ); }; // 
  // and for converting back to the individual indices
  // Branch and FEC methods would just be the same as above
  Int_t GetCSPFromAddress(Int_t i) const { return ( i & 0x1f ); }; // 

  /* // Below is some placeholder info that can later be added
     // in AliEMCALGeometry, together with the Get methods just above 

  // But which CSP (0..31) corresponds to which ALTRO and Channel is not 
  // given anywhere (CSPs are used for APD biases etc).
  // So, we add a conversion method for that here also.
  // The order that the CSPs appear in the data is a bit funky so I include
  // a mapping array instead of some complicated function
  static const int fgkNCSP = 32;
  static const int fgkCspOrder[32] =
    { // just from ALTRO mapping of chips/channels to CSP
      11,  27,  10,  26,  24,   8,  25,   9, // ALTRO 0
      3,  19,   2,  18,  16,   0,  17,   1, // ALTRO 2
      4,  20,   5,  21,  23,   7,  22,   6, // ALTRO 3
      12,  28,  13,  29,  31,  15,  30,  14 // ALTRO 4
    };
  // This method is not used for reconstruction or so, but just for cross-
  // checks with the DCS for the APD biases. 
  int GetCSP(int iALTRO, int iChannel) const 
  { 
    int id = iChannel/2; // 2 channels per tower (low and high gain)
    int ichip = iALTRO;
    if (ichip>=2) { ichip--; } // there is no ALTRO 1; (0,2,3,4 -> 0,1,2,3)
    id += ichip*8; // 8 CSPs per ALTRO
    //return fgkCspOrder[id];
    return id;
  }

  */

};

#endif
