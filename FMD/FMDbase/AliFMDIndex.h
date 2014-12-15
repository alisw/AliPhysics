#ifndef ALIFMDINDEX_H
#define ALIFMDINDEX_H
/** @file    AliFMDIndex.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:37:41 2006
    @brief   FMD detector coordinates
*/
//___________________________________________________________________
//
//  Class that holds an FMD index.  That is, it holds the detector
//  coordinates for a given strip:
//
//     Variable | Type     | Range   | Description
//     ---------+----------+---------+------------------
//     detector | UShort_t | 1-3     | Detector number 
//     ring     | Char_t   | 'I'/'O' | Ring identifier 
//     sector   | UShort_t | 0-39    | Sector number
//     strip    | UShort_t | 0-511   | Strip number
//
#ifndef ROOT_Rtypes
# include <Rtypes.h>
#endif
#ifndef ROOT_TObject
# include <TObject.h>
#endif
#ifndef ROOT_TString
# include <TString.h>
#endif
#include <iosfwd>

//____________________________________________________________________
/** @class AliFMDIndex AliFMDIndex.h <FMD/AliFMDIndex.h>
    @brief FMD detector coordinates 
    @ingroup FMD_base
 */
class AliFMDIndex 
{
public: 
  /** CTOR */
  AliFMDIndex();
  /** Copy CTOR 
      @param o Object to copy from */
  AliFMDIndex(const AliFMDIndex& o);
  /** Constrctor 
      @param detector Detector 
      @param ring     Ring
      @param sector   Sector
      @param strip    Strip */
  AliFMDIndex(UShort_t detector, 
	      Char_t   ring='\0', 
	      UShort_t sector=0, 
	      UShort_t strip=0);
  /** Assignment operator 
      @param o Object to assign from 
      @return Reference to this object  */
  AliFMDIndex& operator=(const AliFMDIndex& o);
  /** Comparison operator 
      @param o Object to compare to 
      @return @c true if these refer to the same index  */
  bool operator==(const AliFMDIndex& o) const;
  /** Comparison operator 
      @param o Object to compare to 
      @return @c true if this is smaller than @a o */
  bool operator<(const AliFMDIndex& o) const;
  /** DTOR */
  virtual ~AliFMDIndex() {}
  /** @return Detector # */
  UShort_t     Detector()	   const { return fDetector; }
  /** @return Ring ID */
  Char_t       Ring()	           const { return fRing;     }
  /** @return sector # */
  UShort_t     Sector()	           const { return fSector;   }
  /** @return strip # */
  UShort_t     Strip()	           const { return fStrip;    }
  /** @param x Detector # */
  void SetDetector(UShort_t x)	   { fHash = -1; fDetector = x; }
  /** @param x Ring ID */
  void SetRing(Char_t x)	   { fHash = -1; fRing = x; }
  /** @param x sector # */
  void SetSector(UShort_t x)	   { fHash = -1; fSector = x; }
  /** @param x strip # */
  void SetStrip(UShort_t x)	   { fHash = -1; fStrip = x; }
  /** Print information 
      @param opt Not used */
  virtual void Print(Option_t* opt="") const;
  /** @return Name */
  const char*  Name() const;
protected:
  Int_t Hash() const;
  UShort_t fDetector;      // (Sub) Detector # (1,2, or 3)
  Char_t   fRing;          // Ring ID ('I' or 'O')
  UShort_t fSector;        // Sector # (phi division)
  UShort_t fStrip;         // Strip # (radial division)
  mutable TString  fName;  //! Cached name
  mutable Int_t    fHash;  //! Cached hash value
  ClassDef(AliFMDIndex, 1) // Base class for FMD digits 
};

//____________________________________________________________________
class AliFMDObjIndex : public TObject, public AliFMDIndex
{
public:
  /** CTOR */
  AliFMDObjIndex() {}
  /** Copy CTOR 
      @param o Object to copy from */
  AliFMDObjIndex(const AliFMDObjIndex& o) : TObject(o), AliFMDIndex(o) {}
  /** Construct from a pure index
      @param o Object to copy from */
  explicit AliFMDObjIndex(const AliFMDIndex& o) : AliFMDIndex(o) {}
  /** Constrctor 
      @param detector Detector 
      @param ring     Ring
      @param sector   Sector
      @param strip    Strip */
  AliFMDObjIndex(UShort_t detector, 
		 Char_t   ring='\0', 
		 UShort_t sector=0, 
		 UShort_t strip=0) 
    : AliFMDIndex(detector, ring, sector, strip)
  {}
  /** DTOR */
  virtual ~AliFMDObjIndex() {}
  AliFMDObjIndex& operator=(const AliFMDObjIndex& o) 
  {
    if (&o == this) return *this;
    AliFMDIndex::operator=(o);
    return *this; 
  }
  /** @return name */
  virtual const char* GetName() const { return AliFMDIndex::Name(); }
  /** sort compare for TCollection's
      @param o Object to compare to
      @return  -1 if this is @e smaller than @a o, 0 if @e equal to 
      @a o, and 1 if this is @e larger than @a o */
  virtual Int_t Compare(const TObject* o) const;
  /** @return always true */
  Bool_t IsSortable() const { return kTRUE; }
  ClassDef(AliFMDObjIndex, 1) // Base class for FMD digits 
};
 
//____________________________________________________________________
inline
bool
AliFMDIndex::operator==(const AliFMDIndex& o) const
{
  return (o.Hash() == Hash());
}

//____________________________________________________________________
inline
bool
AliFMDIndex::operator<(const AliFMDIndex& rhs) const
{
  return (Hash() < rhs.Hash());
}

#if 0
//____________________________________________________________________
inline
bool
operator<(const AliFMDIndex& lhs, const AliFMDIndex& rhs)
{
  return (lhs.Detector() < rhs.Detector() ? true : 
	  (lhs.Ring()    < rhs.Ring()     ? true :
	   (lhs.Sector() < rhs.Sector()   ? true :
	    (lhs.Strip() < rhs.Strip()    ? true : false))));
}
#endif 
#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
//
// EOF
//
