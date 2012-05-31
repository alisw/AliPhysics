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
#ifndef ROOT_TArrayI
# include <TArrayI.h>
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
  /** 
   * Constrctor 
   *
   * @param detector Detector 
   * @param ring     Ring
   * @param sector   Sector
   * @param strip    Strip 
   * @param count    ADC (first sample)
   * @param count2   ADC (second sample, or -1 if not used)
   * @param count3   ADC (third sample, or -1 if not used) 
   * @param refs     Track references
   */
  AliFMDDigit(UShort_t       detector, 
	      Char_t         ring='\0', 
	      UShort_t       sector=0, 
	      UShort_t       strip=0, 
	      UShort_t       count=0, 
	      Short_t        count2=-1, 
	      Short_t        count3=-1, 
	      Short_t        count4=-1, 
	      UShort_t       nrefs=0,
	      const Int_t*   refs=0);
  /** 
   * DTOR 
   */
  virtual ~AliFMDDigit() {}
  /** 
   * @param i # of sample to get 
   * 
   * @return sample # @a i 
   */
  Int_t Count(UShort_t i=0) const;
  /** 
   * 
   * @return ADC count (first sample) 
   */
  UShort_t Count1() const { return fCount1;   }
  /** 
   * 
   * @return ADC count (second sample, or -1 if not used) 
   */
  Short_t  Count2() const { return fCount2;   }
  /** 
   * 
   * @return ADC count (third sample, or -1 if not used) 
   */
  Short_t  Count3() const { return fCount3;   }
  /** 
   * 
   * @return ADC count (third sample, or -1 if not used) 
   */
  Short_t  Count4() const { return fCount4;   }
  /** 
   * 
   * @return Canonical ADC counts 
   */
  UShort_t Counts() const;
  /** 
   * Print info 
   * 
   * @param opt Not used 
   */
  void     Print(Option_t* opt="") const;
  /** 
   * 
   * @return Title 
   */
  const char* GetTitle() const;
  /** 
   * Set the count value 
   * 
   * @param s Sample number 
   * @param c Counts 
   */
  void SetCount(UShort_t s, Short_t c);
  /**
   * Initialize all counts to appropriate values for this oversampling
   * rate.  That is 
   *
   * @verbatim
   *     Rate | Sample 1 | Sample 2 | Sample 3 | Sample 4 
   *     -----+----------+----------+----------+----------
   *     1    | 0        | -1       | -1       | -1
   *     2    | 0        | 0        | -1       | -1
   *     3    | 0        | 0        | 0        | -1
   *     4    | 0        | 0        | 0        | 0
   * @endverbatim
   *
   * @param rate Oversampling rate 
   */
  void SetDefaultCounts(UShort_t rate);
protected:
  UShort_t fCount1;     // Digital signal 
  Short_t  fCount2;     // Digital signal (-1 if not used)
  Short_t  fCount3;     // Digital signal (-1 if not used)
  Short_t  fCount4;     // Digital signal (-1 if not used)
  ClassDef(AliFMDDigit,2)     // Normal FMD digit
};

inline void
AliFMDDigit::SetDefaultCounts(UShort_t rate)
{
  switch (rate) { 
  case 4: fCount4 = 0; // Fall through 
  case 3: fCount3 = 0; // Fall through 
  case 2: fCount2 = 0; // Fall through 
  case 1: fCount1 = 0;
    break;
  default: 
    fCount4 = fCount3 = fCount2 = fCount1 = 0;
    break;
  }
}
inline UShort_t 
AliFMDDigit::Counts() const 
{
  if (fCount4 >= 0) return fCount3;
  if (fCount3 >= 0) return fCount2;
  if (fCount2 >= 0) return fCount2;
  return fCount1;
}

inline Int_t
AliFMDDigit::Count(UShort_t i) const 
{
  switch (i) {
  case 0: return fCount1;
  case 1: return fCount2;
  case 2: return fCount3;
  case 3: return fCount4;
  }
  return -1;
}
inline void
AliFMDDigit::SetCount(UShort_t i, Short_t c)
{
  switch (i) {
  case 0: fCount1 = c; break;
  case 1: fCount2 = c; break;
  case 2: fCount3 = c; break;
  case 3: fCount4 = c; break;
  }
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
