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
class AliFMDUShortMap : public AliFMDMap
{
public:
  AliFMDUShortMap(const AliFMDUShortMap& other);
  AliFMDUShortMap(size_t maxDet = kMaxDetectors, 
		  size_t maxRing= kMaxRings, 
		  size_t maxSec = kMaxSectors, 
		  size_t maxStr = kMaxStrips);
  virtual ~AliFMDUShortMap() { delete [] fData; }
  AliFMDUShortMap& operator=(const AliFMDUShortMap& other);
  virtual void Clear(const UShort_t& val=UShort_t());
  virtual UShort_t& operator()(UShort_t detector, 
			       Char_t   ring, 
			       UShort_t sector, 
			       UShort_t strip);
  virtual const UShort_t& operator()(UShort_t detector, 
				     Char_t   ring, 
				     UShort_t sector, 
				     UShort_t strip) const;
 protected:
  UShort_t* fData;  // The data 
  ClassDef(AliFMDUShortMap, 1) // Cache of edep,hit information per strip
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


