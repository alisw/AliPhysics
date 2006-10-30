#ifndef ALIFMDBOOLMAP_H
#define ALIFMDBOOLMAP_H
/* Copyright (c) 2004, ALICE Experiment @ CERN.
 * All rights reserved
 * See AliFMDBoolMap.cxx for full copyright notice
 * 
 * Created Mon Nov  8 12:51:51 2004 by Christian Holm Christensen
 */
/* $Id$ */
/** @file    AliFMDBoolMap.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:28:59 2006
    @brief   Per strip boolean map
*/
//__________________________________________________________
// 
// Map of Bool_t for each FMD strip
// Used in calibration and the like classes.
// Used amoung other things for dead-channel map
// 
#ifndef ALIFMDMAP_H
# include <AliFMDMap.h>
#endif

/** @class AliFMDBoolMap 
    @brief MAp of per strip boolean values. 
    @ingroup FMD_base
 */
class AliFMDBoolMap : public AliFMDMap
{
public:
  /** Copy constructor 
      @param other Object to copy from. */
  AliFMDBoolMap(const AliFMDBoolMap& other);
  /** Constructor 
      @param maxDet  Number of detectors (3)
      @param maxRing Number of rings (2)
      @param maxSec  Number of sectors (40)
      @param maxStr  Number of strips (20) */
  AliFMDBoolMap(UShort_t maxDet  = kMaxDetectors,
		UShort_t maxRing = kMaxRings,
		UShort_t maxSec  = kMaxSectors,
		UShort_t maxStr  = kMaxStrips);
  /** Destructor */
  virtual ~AliFMDBoolMap() { delete [] fData; }
  /** Assignment operator 
      @param other Object to assign from 
      @return reference to this object.  */
  AliFMDBoolMap& operator=(const AliFMDBoolMap& other);
  /** Reset to value 
      @param v Value to reset from */
  virtual void Reset(const Bool_t& v=Bool_t());
  /** Access operator 
      @param det   Detector 
      @param ring  Ring 
      @param sec   Sector  
      @param str   Strip
      @return  reference value stored for the strip */
  virtual Bool_t& operator()(UShort_t det,
			     Char_t   ring,
			     UShort_t sec,
			     UShort_t str);
  /** Access operator 
      @param det   Detector 
      @param ring  Ring 
      @param sec   Sector  
      @param str   Strip
      @return  value stored for the strip */
  virtual const Bool_t& operator()(UShort_t det,
				   Char_t   ring,
				   UShort_t sec,
				   UShort_t str) const;
protected:
  Int_t  fTotal; // Total number of entries 
  Bool_t* fData;  // [fTotal] The Data
  ClassDef(AliFMDBoolMap,3) // Map of Bool_t data per strip
};

#endif
//__________________________________________________________
// 
// Local Variables:
//   mode: C++
// End:
//
