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
    @ingroup FMD_sim
*/
//
// Contains a pair of energy deposited @c fEdep and number of hits @c
// fN, @c fEdep is the summed energy deposition, and @c fN is the
// number of hits  
#ifndef ALIFMDMAP_H
# include "AliFMDMap.h"
#endif 
#ifndef ALIFMDEDEPHITPAIR_H
# include <AliFMDEdepHitPair.h>
#endif


//____________________________________________________________________
/** @brief Map of Energy deposited, hit information per strip.
    Contains a pair of energy deposited @c fEdep and 
    number of hits @c fN, @c fEdep is the summed energy deposition,
    and @c fN is the number of hits 
    @ingroup FMD_sim
*/
class AliFMDEdepMap : public AliFMDMap
{
public:
  /** Copy constructor 
      @param other Object to copy from. 
      @return  */
  AliFMDEdepMap(const AliFMDEdepMap& other);
  /** Constructor */
  AliFMDEdepMap();
  /** Constructor 
      @param maxDet  Number of detectors (3)
      @param maxRing Number of rings (2)
      @param maxSec  Number of sectors (40)
      @param maxStr  Number of strips (20) */
  AliFMDEdepMap(UShort_t maxDet, 
		UShort_t maxRing= 0, 
		UShort_t maxSec = 0, 
		UShort_t maxStr = 0);
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
  void*     Ptr() const { return reinterpret_cast<void*>(fData); }
protected:
  Int_t     MaxIndex()            const { return fTotal; }
  Int_t     AtAsInt(Int_t idx)    const { return fData[idx].fN; }
  Int_t&    AtAsInt(Int_t idx)          { return AliFMDMap::AtAsInt(idx); }
  Bool_t    AtAsBool(Int_t idx)   const { return fData[idx].fEdep > 0; }
  Bool_t&   AtAsBool(Int_t idx)         { return AliFMDMap::AtAsBool(idx); }
  Float_t   AtAsFloat(Int_t idx)  const { return fData[idx].fEdep; }
  Float_t&  AtAsFloat(Int_t idx)        { return fData[idx].fEdep; }
  UShort_t  AtAsUShort(Int_t idx) const { return fData[idx].fN; }
  UShort_t& AtAsUShort(Int_t idx)       { return fData[idx].fN; }
  Bool_t    IsUShort()            const { return kTRUE; }  
  Bool_t    IsFloat()             const { return kTRUE; }  
  Int_t              fTotal; //  Total number of entries
  AliFMDEdepHitPair* fData;  //[fTotal] The data 
  ClassDef(AliFMDEdepMap, 3) // Cache of edep,hit information per strip
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


