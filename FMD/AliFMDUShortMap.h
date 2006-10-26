#ifndef ALIFMDUSHORTMAP_H
#define ALIFMDUSHORTMAP_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDUShortMap.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:48:18 2006
    @brief   Per strip of unisgned shorts (16 bit) data 
*/
// Map of an integer per strip
// This class stores one short unsigned integer (16 bits) per strip in
// the FMD detectors. 
#ifndef ALIFMDMAP_H
# include "AliFMDMap.h"
#endif 
//____________________________________________________________________
/** @class AliFMDUShortMap 
    @brief Map of an integer per strip
    @ingroup FMD_base
 */
class AliFMDUShortMap : public AliFMDMap
{
public:
  /** Copy constructor 
      @param other Object to copy from.  */
  AliFMDUShortMap(const AliFMDUShortMap& other);
  /** Constructor 
      @param maxDet  Number of detectors (3)
      @param maxRing Number of rings (2)
      @param maxSec  Number of sectors (40)
      @param maxStr  Number of strips (20) */
  AliFMDUShortMap(UShort_t maxDet = kMaxDetectors, 
		  UShort_t maxRing= kMaxRings, 
		  UShort_t maxSec = kMaxSectors, 
		  UShort_t maxStr = kMaxStrips);
  /** Destructor */
  virtual ~AliFMDUShortMap() { delete [] fData; }
  /** Assignment operator 
      @param other Object to assign from 
      @return reference to this object.  */
  AliFMDUShortMap& operator=(const AliFMDUShortMap& other);
  /** Reset to value 
      @param val Value to reset from */
  virtual void Reset(const UShort_t& val=UShort_t());
  /** Access operator 
      @param detector Detector 
      @param ring     Ring 
      @param sector   Sector  
      @param strip    Strip
      @return  reference value stored for the strip */
  virtual UShort_t& operator()(UShort_t detector, 
			       Char_t   ring, 
			       UShort_t sector, 
			       UShort_t strip);
  /** Access operator 
      @param detector Detector 
      @param ring     Ring 
      @param sector   Sector  
      @param strip    Strip
      @return  value stored for the strip */
  virtual const UShort_t& operator()(UShort_t detector, 
				     Char_t   ring, 
				     UShort_t sector, 
				     UShort_t strip) const;
 protected:
  UShort_t    fTotal; // Total number of entries 
  UShort_t* fData;  // [fTotal] The data 
  ClassDef(AliFMDUShortMap, 2) // Cache of edep,hit information per strip
};

#endif 
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//


