#ifndef ALIFMDBASEDIGIT_H
#define ALIFMDBASEDIGIT_H
/** @file    AliFMDBaseDigit.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:37:41 2006
    @brief   Digits for the FMD 
    @ingroup FMD_base
*/
//___________________________________________________________________
//
//  Digits classes for the FMD
//  AliFMDBaseDigit - base class 
//  AliFMDDigit     - Normal (smeared) digit             
//  AliFMDSDigit    - Summable (non-smeared) digit             
//
#ifndef ROOT_TObject
# include <TObject.h>
#endif
#ifndef ROOT_TString
# include <TString.h>
#endif

//____________________________________________________________________
/** @class AliFMDBaseDigit AliFMDDigit.h <FMD/AliFMDDigit.h>
    @brief base class for digits 
    @ingroup FMD_base
 */
class AliFMDBaseDigit : public TObject 
{
public: 
  /** CTOR */
  AliFMDBaseDigit();
  /** Constrctor 
      @param detector Detector 
      @param ring     Ring
      @param sector   Sector
      @param strip    Strip */
  AliFMDBaseDigit(UShort_t detector, 
		  Char_t   ring='\0', 
		  UShort_t sector=0, 
		  UShort_t strip=0);
  /** DTOR */
  virtual ~AliFMDBaseDigit() {}
  /** @return Detector # */
  UShort_t     Detector()	   const { return fDetector; }
  /** @return Ring ID */
  Char_t       Ring()	           const { return fRing;     }
  /** @return sector # */
  UShort_t     Sector()	           const { return fSector;   }
  /** @return strip # */
  UShort_t     Strip()	           const { return fStrip;    }
  /** Print information 
      @param opt Not used */
  virtual void Print(Option_t* opt="") const;
  /** @return Name */
  const char*  GetName() const;
  /** @param rhs Other digit to compare to 
      @return -1 if this is less than  @a rhs, 0 if the refer to the
      same, and 1 if @a rhs is larger than this */
  Int_t Compare(const TObject* o) const;
  /** @return Always true */ 
  Bool_t IsSortable() const { return kTRUE; }
protected:
  ULong_t  Hash() const;
  UShort_t fDetector;  // (Sub) Detector # (1,2, or 3)
  Char_t   fRing;      // Ring ID ('I' or 'O')
  UShort_t fSector;    // Sector # (phi division)
  UShort_t fStrip;     // Strip # (radial division)
  mutable TString  fName;      //! Name (cached, but not stored) 
  ClassDef(AliFMDBaseDigit, 2) // Base class for FMD digits 
};

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
