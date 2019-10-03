#ifndef ALIFMDSTRIPINDEX_H
#define ALIFMDSTRIPINDEX_H
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

// Struct to encode a strip address into one integer
// developed by Christian Holm Christensen (cholm@nbi.dk).
// 
// The functions are static to ensure applicability from
// anywhere. This is needed to smoothly store strip addresses in track
// references. 
//
// Added by Hans H. Dalsgaard (hans.dalsgaard@cern.ch) 


class AliFMDStripIndex
{
public:
  enum { 
    // Mask of ID
    kIdMask = 0x0007FFFF,
    // Mask of energy 
    kEMask  = 0xFFF80000,
    // Offset of energy 
    kEOffset = 19
  };
  /** 
   * Constructor
   * 
   */  
  AliFMDStripIndex() {}
  /** 
   * Destructor 
   * 
   */
  virtual ~AliFMDStripIndex() {}
  /** 
   * Pack an identifier from detector coordinates
   * 
   * @param det  Detector
   * @param rng  Ring
   * @param sec  Sector
   * @param str  Strip
   * 
   * @return Packed identifier 
   */
  static UInt_t Pack(UShort_t det, Char_t rng, UShort_t sec, UShort_t str) 
  {
    UInt_t irg  = (rng == 'I' || rng == 'i' ? 0 : 1);
    UInt_t id   = (((str & 0x1FF) <<  0) | 
		   ((sec & 0x03F) <<  9) | 
		   ((irg & 0x001) << 16) | 
		   ((det & 0x003) << 17));
    return (id & kIdMask);
  }
  /** 
   * Unpack an identifier to detector coordinates
   * 
   * @param id   Identifier to unpack
   * @param det  On return, the detector
   * @param rng  On return, the ring
   * @param sec  On return, the sector
   * @param str  On return, the strip
   */  
  static void Unpack(UInt_t id, 
		     UShort_t& det, Char_t& rng, UShort_t& sec, UShort_t& str)
  {
    UInt_t tmp = (kIdMask & id);
    str = ((tmp >>  0) & 0x1FF);
    sec = ((tmp >>  9) & 0x03F);
    rng = ((tmp >> 16) & 0x001) ? 'O' : 'I';
    det = ((tmp >> 17) & 0x003);
  }
  ClassDef(AliFMDStripIndex,1)
};
#endif
//
// Local Variables:
//   mode: C++
// End:
//
