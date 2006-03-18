#ifndef ALIFMDEDEPMAP_H
#define ALIFMDEDEPMAP_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIFMDMAP_H
# include "AliFMDMap.h"
#endif 
//____________________________________________________________________
//
// Cache of Energy deposited, hit information per strip.
// Contains a pair of fEdep and fN
// fEdep is the summed energy deposition, and fN is the number of hits
//

//____________________________________________________________________
class AliFMDEdepHitPair 
{
public:
  Float_t  fEdep; // summed energy deposition
  UShort_t fN;    // Number of hits
  AliFMDEdepHitPair() : fEdep(0), fN(0) {}
  AliFMDEdepHitPair& operator=(const AliFMDEdepHitPair& o) 
  { fEdep = o.fEdep; fN    = o.fN; return *this; }
  AliFMDEdepHitPair(const AliFMDEdepHitPair& o) : fEdep(o.fEdep), fN(o.fN) {}
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
  virtual void Reset();
  virtual void Reset(const AliFMDEdepHitPair& val);
  virtual AliFMDEdepHitPair& operator()(UShort_t detector, 
				     Char_t   ring, 
				     UShort_t sector, 
				     UShort_t strip);
  virtual const AliFMDEdepHitPair& operator()(UShort_t detector, 
					   Char_t   ring, 
					   UShort_t sector, 
					   UShort_t strip) const;
protected:
  size_t             fTotal; //  
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


