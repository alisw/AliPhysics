#ifndef ALIFMDFLOATMAP_H
#define ALIFMDFLOATMAP_H
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
// Array of floats indexed by strip identifier.
//
class AliFMDFloatMap : public AliFMDMap
{
public:
  AliFMDFloatMap(size_t  maxDet = kMaxDetectors, 
		 size_t  maxRing= kMaxRings, 
		 size_t  maxSec = kMaxSectors, 
		 size_t  maxStr = kMaxStrips);
  AliFMDFloatMap(const AliFMDFloatMap& o);
  virtual ~AliFMDFloatMap() { delete [] fData; }
  AliFMDFloatMap& operator=(const AliFMDFloatMap& o);
  virtual void Reset(const Float_t& v=Float_t());
  virtual Float_t& operator()(UShort_t det,
			      Char_t   ring,
			      UShort_t sec,
			      UShort_t str);
  virtual const Float_t& operator()(UShort_t det,
				    Char_t   ring,
				    UShort_t sec,
				    UShort_t str) const;
protected:
  Float_t* fData;
  ClassDef(AliFMDFloatMap,1) // Map of floats
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//

