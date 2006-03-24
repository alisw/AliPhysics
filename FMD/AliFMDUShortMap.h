#ifndef ALIFMDUSHORTMAP_H
#define ALIFMDUSHORTMAP_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
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
  AliFMDUShortMap(size_t maxDet = kMaxDetectors, 
		  size_t maxRing= kMaxRings, 
		  size_t maxSec = kMaxSectors, 
		  size_t maxStr = kMaxStrips);
  /** Destructor */
  virtual ~AliFMDUShortMap() { delete [] fData; }
  /** Assignment operator 
      @param other Object to assign from 
      @return reference to this object.  */
  AliFMDUShortMap& operator=(const AliFMDUShortMap& other);
  /** Reset to value 
      @param v Value to reset from */
  virtual void Reset(const UShort_t& val=UShort_t());
  /** Access operator 
      @param det   Detector 
      @param ring  Ring 
      @param sec   Sector  
      @param str   Strip
      @return  reference value stored for the strip */
  virtual UShort_t& operator()(UShort_t detector, 
			       Char_t   ring, 
			       UShort_t sector, 
			       UShort_t strip);
  /** Access operator 
      @param det   Detector 
      @param ring  Ring 
      @param sec   Sector  
      @param str   Strip
      @return  value stored for the strip */
  virtual const UShort_t& operator()(UShort_t detector, 
				     Char_t   ring, 
				     UShort_t sector, 
				     UShort_t strip) const;
 protected:
  size_t    fTotal; // Total number of entries 
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


