#ifndef ALIFMDBOOLMAP_H
#define ALIFMDBOOLMAP_H
/* Copyright (c) 2004, ALICE Experiment @ CERN.
 * All rights reserved
 * See AliFMDBoolMap.cxx for full copyright notice
 * 
 * Created Mon Nov  8 12:51:51 2004 by Christian Holm Christensen
 */
/* $Id$ */
//__________________________________________________________
// 
// Map of Bool_t for each FMD strip
// Used in calibration and the like classes.
// Used amoung other things for dead-channel map
// 
#ifndef ALIFMDMAP_H
# include <AliFMDMap.h>
#endif

class AliFMDBoolMap : public AliFMDMap
{
public:
  AliFMDBoolMap(const AliFMDBoolMap& other);
  AliFMDBoolMap(size_t maxDet  = kMaxDetectors,
		size_t maxRing = kMaxRings,
		size_t maxSec  = kMaxSectors,
		size_t maxStr  = kMaxStrips);
  virtual ~AliFMDBoolMap() { delete [] fData; }
  AliFMDBoolMap& operator=(const AliFMDBoolMap& other);
  virtual void Reset(const Bool_t& v=Bool_t());
  virtual Bool_t& operator()(UShort_t det,
			     Char_t   ring,
			     UShort_t sec,
			     UShort_t str);
  virtual const Bool_t& operator()(UShort_t det,
				   Char_t   ring,
				   UShort_t sec,
				   UShort_t str) const;
protected:
  Bool_t* fData; // The Data
  ClassDef(AliFMDBoolMap,1) // Map of Bool_t data per strip
};

#endif
//__________________________________________________________
// 
// Local Variables:
//   mode: C++
// End:
//
