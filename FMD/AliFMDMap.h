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
//____________________________________________________________________
//
// Base class for caches of per-strip information.
//

class AliFMDMap : public TObject 
{
public:
  enum { 
    kMaxDetectors = 3, 
    kMaxRings     = 2, 
    kMaxSectors   = 40, 
    kMaxStrips    = 512
  };
  AliFMDMap(size_t maxDet = kMaxDetectors, 
	    size_t maxRing= kMaxRings, 
	    size_t maxSec = kMaxSectors, 
	    size_t maxStr = kMaxStrips);
  virtual ~AliFMDMap() {}
  Int_t CheckIndex(size_t det, Char_t ring, size_t sec, size_t str) const;
protected:
  size_t CalcIndex(size_t det, Char_t ring, size_t sec, size_t str) const;
  size_t fMaxDetectors;             // Maximum # of detectors
  size_t fMaxRings;                 // Maximum # of rings
  size_t fMaxSectors;               // Maximum # of sectors
  size_t fMaxStrips;                // Maximum # of strips
  ClassDef(AliFMDMap, 1) // Cache of per strip information
};

#ifdef  MAY_USE_TEMPLATES
#ifndef __VECTOR__
# include <vector>
#endif 
//____________________________________________________________________
//
// Class template for classes that cache per strip information.  
// Access to the data goes via 
//
//   Type& AliFMDMap<Type>::operator()(size_t detector,
//                                     Char_t ring, 
//                                     size_t sector,
//                                     size_t strip);
// 
// (as well as a const version of this member function). 
// The elements can be reset to the default value by calling 
// AliFMDMap<Type>::Clear().  This resets the values to `Type()'. 
//
template <typename Type> 
class AliFMDMap : public TObject 
{
public:
  AliFMDMap(size_t maxDet=3, size_t maxRing=2, size_t maxSec=40, 
	    size_t maxStr=512);
  virtual ~AliFMDMap() {}
  void Clear(const Type& val=Type());
  Type& operator()(size_t det, Char_t ring, size_t sec, size_t str);
  const Type& operator()(size_t det, Char_t ring, size_t sec, size_t str)const;
private:
  typedef std::vector<Type> ValueVector; // Type of container
  ValueVector fValues;                   // Contained values
  size_t      fMaxDetectors;             // Maximum # of detectors
  size_t      fMaxRings;                 // Maximum # of rings
  size_t      fMaxSectors;               // Maximum # of sectors
  size_t      fMaxStrips;                // Maximum # of strips
  
  size_t CalcIndex(size_t det, Char_t ring, size_t sec, size_t str) const;
  ClassDef(AliFMDMap, 0); // Map of FMD index's to values 
};


//____________________________________________________________________
template <typename Type>
inline 
AliFMDMap<Type>::AliFMDMap(size_t maxDet, 
			   size_t maxRing, 
			   size_t maxSec, 
			   size_t maxStr)
  : fValues(maxDet * maxRing * maxSec * maxStr), 
    fMaxDetectors(maxDet), 
    fMaxRings(maxRing), 
    fMaxSectors(maxSec), 
    fMaxStrips(maxStr)
{
  // Construct a map
  //
  // Parameters:
  //     maxDet       Maximum # of detectors
  //     maxRinf      Maximum # of rings
  //     maxSec       Maximum # of sectors
  //     maxStr       Maximum # of strips
}


//____________________________________________________________________
template <typename Type>
inline size_t 
AliFMDMap<Type>::CalcIndex(size_t det, Char_t ring, size_t sec, size_t str) const
{
  // Calculate index into storage from arguments. 
  // 
  // Parameters: 
  //     det       Detector #
  //     ring      Ring ID
  //     sec       Sector # 
  //     str       Strip # 
  //
  // Returns appropriate index into storage 
  //
  size_t ringi = (ring == 'I' ||  ring == 'i' ? 0 : 1);
  size_t idx = 
    (det + fMaxDetectors * (ringi + fMaxRings * (sec + fMaxSectors * str)));
  if (idx >= fMaxDetectors * fMaxRings * fMaxSectors * fMaxStrips) {
    Fatal("CalcIndex", "Index (%d,'%c',%d,%d) out of bounds, "
	  "in particular the %s index", 
	  det, ring, sec, str, 
	  (det >= fMaxDetectors ? "Detector" : 
	   (ringi >= fMaxRings ? "Ring" : 
	    (sec >= fMaxSectors ? "Sector" : "Strip"))));
    return 0;
  }
  return idx;
}

//____________________________________________________________________
template <typename Type>
inline void
AliFMDMap<Type>::Clear(const Type& val) 
{
  // Resets stored values to the default value for that type 
  for (size_t i = 0; i < fValues.size(); ++i) fValues[i] = val;
}

//____________________________________________________________________
template <typename Type>
inline Type& 
AliFMDMap<Type>::operator()(size_t det, Char_t ring, size_t sec, size_t str)
{
  // Parameters: 
  //     det       Detector #
  //     ring      Ring ID
  //     sec       Sector # 
  //     str       Strip # 
  //
  // Returns data[det][ring][sec][str]
  return fValues[CalcIndex(det, ring, sec, str)];
}

//____________________________________________________________________
template <typename Type>
inline const Type& 
AliFMDMap<Type>::operator()(size_t det, 
			    Char_t ring, 
			    size_t sec, 
			    size_t str) const
{
  // Parameters: 
  //     det       Detector #
  //     ring      Ring ID
  //     sec       Sector # 
  //     str       Strip # 
  //
  // Returns data[det][ring][sec][str]
  return fValues[CalcIndex(det, ring, sec, str)];
}


//____________________________________________________________________
// 
// Some specialisations 
//
typedef AliFMDMap<UShort_t> AliFMDAdcMap;
typedef AliFMDMap<std::pair<Float_t, UShort_t> > AliFMDEdepMap;

#endif
#endif 
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//


