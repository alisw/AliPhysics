#ifndef ALIFMDSDIGIT_H
#define ALIFMDSDIGIT_H
/** @file    AliFMDSDigit.h
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
// #ifndef ROOT_TArrayI
// # include <TArrayI.h>
// #endif

//____________________________________________________________________
/** @class AliFMDSDigit AliFMDDigit.h <FMD/AliFMDDigit.h>
    @brief class for summable digits 
    @ingroup FMD_base
 */
class AliFMDSDigit : public AliFMDBaseDigit
{
public:
  /** 
   * CTOR 
   */
  AliFMDSDigit();
  /** 
   * Constrctor 
   *
   * @param detector Detector 
   * @param ring     Ring
   * @param sector   Sector
   * @param strip    Strip 
   * @param edep     Energy deposited 
   * @param count    ADC (first sample)
   * @param count2   ADC (second sample, or -1 if not used)
   * @param count3   ADC (third sample, or -1 if not used) 
   */
  AliFMDSDigit(UShort_t       detector, 
	       Char_t         ring='\0', 
	       UShort_t       sector=0, 
	       UShort_t       strip=0, 
	       Float_t        edep=0,
	       UShort_t       count=0, 
	       Short_t        count2=-1, 
	       Short_t        count3=-1,
	       Short_t        count4=-1,
	       UShort_t       npart=0,
	       UShort_t       nprim=0, 
	       const Int_t*   refs=0);
  /** 
   * DTOR 
   */
  virtual ~AliFMDSDigit() {}
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
   * 
   * @return Energy deposited 
   */
  Float_t  Edep() const { return fEdep;     }
  /** 
   * 
   * @return Number of particles that hit this strip 
   */
  UShort_t NParticles() const { return fNParticles; }
  /** 
   * 
   * @return Number of primary particles that hit this strip 
   */
  UShort_t NPrimaries() const { return fNPrimaries; }
#if 0
  /** @return the track labels */
  const TArrayI& TrackLabels() const { return fLabels; }
#endif
  /** 
   * Print info 
   *
   * @param opt Not used 
   */
  void     Print(Option_t* opt="") const;
  /** 
   * Set the count value 
   * 
   * @param s Sample number 
   * @param c Counts 
   */
  void SetCount(UShort_t s, Short_t c);
protected:
  Float_t  fEdep;       // Energy deposited 
  UShort_t fCount1;     // Digital signal 
  Short_t  fCount2;     // Digital signal (-1 if not used)
  Short_t  fCount3;     // Digital signal (-1 if not used)
  Short_t  fCount4;     // Digital signal (-1 if not used)
  UShort_t fNParticles; // Total number of particles that hit this strip
  UShort_t fNPrimaries; // Number of primary particles that his this strip
#if 0
  TArrayI  fLabels;     // MC-truth track labels
#endif
  ClassDef(AliFMDSDigit,5)     // Summable FMD digit
};
  
inline UShort_t 
AliFMDSDigit::Counts() const 
{
  if (fCount4 >= 0) return fCount3;
  if (fCount3 >= 0) return fCount2;
  if (fCount2 >= 0) return fCount2;
  return fCount1;
}

inline void
AliFMDSDigit::SetCount(UShort_t i, Short_t c)
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
