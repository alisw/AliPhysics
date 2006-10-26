#ifndef ALIFMDDIGIT_H
#define ALIFMDDIGIT_H
/** @file    AliFMDDigit.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:37:41 2006
    @brief   Digits for the FMD 
*/
//___________________________________________________________________
//
//  Digits classes for the FMD
//  AliFMDBaseDigit - base class 
//  AliFMDDigit     - Normal (smeared) digit             
//  AliFMDSDigit    - Summable (non-smeared) digit             
//
#ifndef ALIFMDBASEDIGIT_H
# include <AliFMDBaseDigit.h>
#endif

//____________________________________________________________________
/** @class AliFMDDigit AliFMDDigit.h <FMD/AliFMDDigit.h>
    @brief class for digits 
    @ingroup FMD_base
 */
class AliFMDDigit : public AliFMDBaseDigit
{
public:
  /** CTOR */
  AliFMDDigit();
  /** Constrctor 
      @param detector Detector 
      @param ring     Ring
      @param sector   Sector
      @param strip    Strip 
      @param count    ADC (first sample)
      @param count2   ADC (second sample, or -1 if not used)
      @param count3   ADC (third sample, or -1 if not used) */
  AliFMDDigit(UShort_t detector, 
	      Char_t   ring='\0', 
	      UShort_t sector=0, 
	      UShort_t strip=0, 
	      UShort_t count=0, 
	      Short_t  count2=-1, 
	      Short_t  count3=-1);
  /** DTOR */
  virtual ~AliFMDDigit() {}
  /** @param i # of sample to get 
      @return sample # @a i */
  Int_t Count(UShort_t i=0) const;
  /** @return ADC count (first sample) */
  UShort_t Count1()                const { return fCount1;   }
  /** @return ADC count (second sample, or -1 if not used) */
  Short_t  Count2()                const { return fCount2;   }
  /** @return ADC count (third sample, or -1 if not used) */
  Short_t  Count3()                const { return fCount3;   }
  /** @return Canonical ADC counts */
  UShort_t Counts()                const;
  /** Print info 
      @param opt Not used */
  void     Print(Option_t* opt="") const;
  /** @return Title */
  const char* GetTitle() const;
protected:
  UShort_t fCount1;     // Digital signal 
  Short_t  fCount2;     // Digital signal (-1 if not used)
  Short_t  fCount3;     // Digital signal (-1 if not used)
  ClassDef(AliFMDDigit,1)     // Normal FMD digit
};

inline UShort_t 
AliFMDDigit::Counts() const 
{
  return fCount1 
    + (fCount2 >= 0 ? fCount2 : 0)
    + (fCount3 >= 0 ? fCount3 : 0);
}

inline Int_t
AliFMDDigit::Count(UShort_t i) const 
{
  switch (i) {
  case 0: return fCount1;
  case 1: return fCount2;
  case 2: return fCount3;
  }
  return -1;
}

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
