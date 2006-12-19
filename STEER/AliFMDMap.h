#ifndef ALIFMDMAP_H
#define ALIFMDMAP_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ROOT_TObject
# include <TObject.h>
#endif 
class TFile;

//____________________________________________________________________
/** @class AliFMDMap
    @brief Base class for caches of per-strip information.
    @ingroup FMD_data
    This is used to index a strip. Data stored depends on derived
    class.  */
class AliFMDMap : public TObject 
{
public:
  enum { 
    /** Default maximum detector number */
    kMaxDetectors = 3, 
    /** Default maximum number of rings */
    kMaxRings     = 2, 
    /** Default maximum number of sectors */
    kMaxSectors   = 40, 
    /** Default maximum number of strips */
    kMaxStrips    = 512
  };
  /** Constructor 
      @param maxDet  Maximum allowed detector number
      @param maxRing Maximum number of rings
      @param maxSec  Maximum number of sectors
      @param maxStr  Maximum number of strips
      @return  */
  AliFMDMap(UShort_t maxDet = kMaxDetectors, 
	    UShort_t maxRing= kMaxRings, 
	    UShort_t maxSec = kMaxSectors, 
	    UShort_t maxStr = kMaxStrips);
  /** Destructor */
  virtual ~AliFMDMap() {}
  /** @return  Maximum detector number */
  UShort_t MaxDetectors() const { return fMaxDetectors; }
  /** @return  Maximum number of rings */
  UShort_t MaxRings()     const { return fMaxRings; }
  /** @return  Maximum number of sectors */
  UShort_t MaxSectors()   const { return fMaxSectors; }
  /** @return  Maximum number of strip */
  UShort_t MaxStrips()    const { return fMaxStrips; }
  /** Calculate, check, and return index for strip.  If the index is
      invalid, -1 is returned
      @param det   Detector number
      @param ring  Ring identifier
      @param sec   Sector number 
      @param str   Strip number
      @return  Unique index, or -1 in case of errors */
  Int_t  CheckIndex(UShort_t det, Char_t ring, 
		    UShort_t sec, UShort_t str) const;
  /** Check if we need UShort_t hack 
      @param file File this object was read from */
  void CheckNeedUShort(TFile* file);
  enum {
    /** In case of version 2 of this class, this bit should be set. */
    kNeedUShort = 14
  };
protected:
  /** Calculate index and return 
      @param det  Detector number
      @param ring Ring identifier 
      @param sec  Sector number 
      @param str  Strip number 
      @return  Index (not checked) */
  Int_t CalcIndex(UShort_t det, Char_t ring, 
		  UShort_t sec, UShort_t str) const;
  UShort_t fMaxDetectors;             // Maximum # of detectors
  UShort_t fMaxRings;                 // Maximum # of rings
  UShort_t fMaxSectors;               // Maximum # of sectors
  UShort_t fMaxStrips;                // Maximum # of strips

  ClassDef(AliFMDMap, 3) // Cache of per strip information
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


