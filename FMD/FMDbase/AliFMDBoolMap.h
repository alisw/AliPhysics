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
  /** 
   * Constructor 
   */
  AliFMDBoolMap();
  /** Constructor 
      @param maxDet  Number of detectors (3)
      @param maxRing Number of rings (2)
      @param maxSec  Number of sectors (40)
      @param maxStr  Number of strips (20) */
  AliFMDBoolMap(UShort_t maxDet,
		UShort_t maxRing = 0,
		UShort_t maxSec  = 0,
		UShort_t maxStr  = 0);
  /** Destructor */
  virtual ~AliFMDBoolMap() { if (fData) delete [] fData; }
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
  Bool_t* Data() const { return fData; }
  Int_t   Total() const { return fTotal; }
  void*   Ptr() const { return reinterpret_cast<void*>(fData); }
protected:
  Int_t     MaxIndex()            const { return fTotal; }
  Bool_t    AtAsBool(Int_t idx)   const { return fData[idx]; }
  Bool_t&   AtAsBool(Int_t idx)         { return fData[idx]; }
  Bool_t    IsBool()              const { return kTRUE; }  
  Int_t     AtAsInt(Int_t idx)    const { return fData[idx] ? 1   : 0;   }
  Float_t   AtAsFloat(Int_t idx)  const { return fData[idx] ? 1.F : 0.F; }
  UShort_t  AtAsUShort(Int_t idx) const { return fData[idx] ? 1   : 0;   }
  Int_t&    AtAsInt(Int_t idx)          { return AliFMDMap::AtAsInt(idx);    }
  Float_t&  AtAsFloat(Int_t idx)        { return AliFMDMap::AtAsFloat(idx);  }
  UShort_t& AtAsUShort(Int_t idx)       { return AliFMDMap::AtAsUShort(idx); }
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
