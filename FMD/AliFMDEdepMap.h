#ifndef ALIFMDEDEPMAP_H
#define ALIFMDEDEPMAP_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIFMDMAP_H
# include <AliFMDMap.h>
#endif 
//____________________________________________________________________
//
// Cache of Energy deposited, hit information perr strip
//

//____________________________________________________________________
class AliFMDEdepHitPair 
{
 public:
  Float_t  fEdep;
  UShort_t fN;
  AliFMDEdepHitPair() : fEdep(0), fN(0) {}
};

//____________________________________________________________________
class AliFMDEdepMap : public AliFMDMap
{
public:
  AliFMDEdepMap(const AliFMDEdepMap& other);
  AliFMDEdepMap(size_t maxDet = kMaxDetectors, 
		size_t maxRing= kMaxRings, 
		size_t maxSec = kMaxSectors, 
		size_t maxStr = kMaxStrips);
  virtual ~AliFMDEdepMap() { delete [] fData; }
  AliFMDEdepMap& operator=(const AliFMDEdepMap& other);
  virtual void Clear(const AliFMDEdepHitPair& val=AliFMDEdepHitPair());
  virtual AliFMDEdepHitPair& operator()(UShort_t detector, 
				     Char_t   ring, 
				     UShort_t sector, 
				     UShort_t strip);
  virtual const AliFMDEdepHitPair& operator()(UShort_t detector, 
					   Char_t   ring, 
					   UShort_t sector, 
					   UShort_t strip) const;
protected:
  AliFMDEdepHitPair* fData;  // The data 
  ClassDef(AliFMDEdepMap, 1) // Cache of edep,hit information per strip
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


