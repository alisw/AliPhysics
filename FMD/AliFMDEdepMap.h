#ifndef ALIFMDEDEPMAP_H
#define ALIFMDEDEPMAP_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDEdepMap.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:39:50 2006
    @brief   Per strip map of energy deposited and number of hits 
*/
#ifndef ALIFMDMAP_H
# include "AliFMDMap.h"
#endif 
//____________________________________________________________________
/** @brief Cache of Energy deposited, hit information per strip.
    Contains a pair of energy deposited @c fEdep and 
    number of hits @c fN, @c fEdep is the summed energy deposition,
    and @c fN is the number of hits 
    @ingroup FMD_sim
*/
class AliFMDEdepHitPair 
{
public:
  Float_t  fEdep; // summed energy deposition
  UShort_t fN;    // Number of hits
  /** CTOR  */
  AliFMDEdepHitPair() : fEdep(0), fN(0) {}
  /** DTOR */
  virtual ~AliFMDEdepHitPair() {}
  /** Assignment operator 
      @param o Object to assign from 
      @return Reference to this object */
  AliFMDEdepHitPair& operator=(const AliFMDEdepHitPair& o) 
  { fEdep = o.fEdep; fN    = o.fN; return *this; }
  /** Copy CTOR 
      @param o Object to copy from */
  AliFMDEdepHitPair(const AliFMDEdepHitPair& o) : fEdep(o.fEdep), fN(o.fN) {}
  ClassDef(AliFMDEdepHitPair, 1)
};

//____________________________________________________________________
/** @brief Map of Energy deposited, hit information per strip.
    Contains a pair of energy deposited @c fEdep and 
    number of hits @c fN, @c fEdep is the summed energy deposition,
    and @c fN is the number of hits */
class AliFMDEdepMap : public AliFMDMap
{
public:
  /** Copy constructor 
      @param other Object to copy from. 
      @return  */
  AliFMDEdepMap(const AliFMDEdepMap& other);
  /** Constructor 
      @param maxDet  Number of detectors (3)
      @param maxRing Number of rings (2)
      @param maxSec  Number of sectors (40)
      @param maxStr  Number of strips (20) */
  AliFMDEdepMap(size_t maxDet = kMaxDetectors, 
		size_t maxRing= kMaxRings, 
		size_t maxSec = kMaxSectors, 
		size_t maxStr = kMaxStrips);
  /** DTOR */
  virtual ~AliFMDEdepMap() { delete [] fData; }
  AliFMDEdepMap& operator=(const AliFMDEdepMap& other);
  /** Reset to default */
  virtual void Reset();
  /** Reset to value 
      @param val Value to reset from */
  virtual void Reset(const AliFMDEdepHitPair& val);
  /** Access operator 
      @param detector Detector 
      @param ring     Ring 
      @param sector   Sector  
      @param strip    Strip
      @return  reference value stored for the strip */
  virtual AliFMDEdepHitPair& operator()(UShort_t detector, 
					Char_t   ring, 
					UShort_t sector, 
					UShort_t strip);
  /** Access operator 
      @param detector Detector 
      @param ring     Ring 
      @param sector   Sector  
      @param strip    Strip
      @return value stored for the strip */
  virtual const AliFMDEdepHitPair& operator()(UShort_t detector, 
					      Char_t   ring, 
					      UShort_t sector, 
					      UShort_t strip) const;
protected:
  size_t             fTotal; //  Total number of entries
  AliFMDEdepHitPair* fData;  //[fTotal] The data 
  ClassDef(AliFMDEdepMap, 2) // Cache of edep,hit information per strip
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


