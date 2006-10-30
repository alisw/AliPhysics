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
// This is used to index a strip. 
// Data stored depends on derived class. 
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
  AliFMDMap(UShort_t maxDet = kMaxDetectors, 
	    UShort_t maxRing= kMaxRings, 
	    UShort_t maxSec = kMaxSectors, 
	    UShort_t maxStr = kMaxStrips);
  virtual ~AliFMDMap() {}
  UShort_t MaxDetectors() const { return fMaxDetectors; }
  UShort_t MaxRings()     const { return fMaxRings; }
  UShort_t MaxSectors()   const { return fMaxSectors; }
  UShort_t MaxStrips()    const { return fMaxStrips; }
  Int_t  CheckIndex(UShort_t det, Char_t ring, 
		    UShort_t sec, UShort_t str) const;
protected:
  Int_t CalcIndex(UShort_t det, Char_t ring, 
		   UShort_t sec, UShort_t str) const;
  UShort_t fMaxDetectors;             // Maximum # of detectors
  UShort_t fMaxRings;                 // Maximum # of rings
  UShort_t fMaxSectors;               // Maximum # of sectors
  UShort_t fMaxStrips;                // Maximum # of strips
  ClassDef(AliFMDMap, 2) // Cache of per strip information
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
//   Type& AliFMDMap<Type>::operator()(UShort_t detector,
//                                     Char_t ring, 
//                                     UShort_t sector,
//                                     UShort_t strip);
// 
// (as well as a const version of this member function). 
// The elements can be reset to the default value by calling 
// AliFMDMap<Type>::Clear().  This resets the values to `Type()'. 
//
template <typename Type> 
class AliFMDMap : public TObject 
{
public:
  AliFMDMap(UShort_t maxDet=3, UShort_t maxRing=2, UShort_t maxSec=40, 
	    UShort_t maxStr=512);
  virtual ~AliFMDMap() {}
  void Clear(const Type& val=Type());
  Type& operator()(UShort_t det, Char_t ring, UShort_t sec, UShort_t str);
  const Type& operator()(UShort_t det, Char_t ring, 
			 UShort_t sec, UShort_t str)const;
private:
  typedef std::vector<Type> ValueVector; // Type of container
  ValueVector fValues;                   // Contained values
  UShort_t      fMaxDetectors;             // Maximum # of detectors
  UShort_t      fMaxRings;                 // Maximum # of rings
  UShort_t      fMaxSectors;               // Maximum # of sectors
  UShort_t      fMaxStrips;                // Maximum # of strips
  
  UShort_t CalcIndex(UShort_t det, Char_t ring, 
		     UShort_t sec, UShort_t str) const;
  ClassDef(AliFMDMap, 0); // Map of FMD index's to values 
};


//____________________________________________________________________
template <typename Type>
inline 
AliFMDMap<Type>::AliFMDMap(UShort_t maxDet, 
			   UShort_t maxRing, 
			   UShort_t maxSec, 
			   UShort_t maxStr)
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
inline UShort_t 
AliFMDMap<Type>::CalcIndex(UShort_t det, Char_t ring, 
			   UShort_t sec, UShort_t str) const
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
  UShort_t ringi = (ring == 'I' ||  ring == 'i' ? 0 : 1);
  UShort_t idx = 
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
  for (UShort_t i = 0; i < fValues.size(); ++i) fValues[i] = val;
}

//____________________________________________________________________
template <typename Type>
inline Type& 
AliFMDMap<Type>::operator()(UShort_t det, Char_t ring, 
			    UShort_t sec, UShort_t str)
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
AliFMDMap<Type>::operator()(UShort_t det, 
			    Char_t ring, 
			    UShort_t sec, 
			    UShort_t str) const
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


