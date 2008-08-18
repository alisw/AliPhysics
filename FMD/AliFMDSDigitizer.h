#ifndef ALIFMDSDIGITIZER_H
#define ALIFMDSDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
// Classses to make Hits into digits and summable digits. 
//    
//    Digits consists of
//    - Detector #
//    - Ring ID                                             
//    - Sector #     
//    - Strip #
//    - ADC count in this channel                                  
//
//    Summable digits consists of	
//    - Detector #
//    - Ring ID                                             
//    - Sector #     
//    - Strip #
//    - Total energy deposited in the strip
//    - ADC count in this channel                                  
//
/** @file    AliFMDSDigitizer.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:38:26 2006
    @brief   FMD Digitizers declaration
    @ingroup FMD_sim
*/
#ifndef ALIFMDBASEDIGITIZER_H
# include <AliFMDBaseDigitizer.h>
#endif

//====================================================================
/** @class AliFMDSDigitizer AliFMDDigitizer.h <FMD/AliFMDDigitizer.h>
    @brief Concrete implementation to make summable digits. 
    See also class documentation of AliFMDBaseDigitizer 
    @ingroup FMD_sim
 */
class AliFMDSDigitizer : public AliFMDBaseDigitizer 
{
public:
  /** CTOR */
  AliFMDSDigitizer() : AliFMDBaseDigitizer() {}
  /** CTOR 
      @param manager Manager of digitization */
  AliFMDSDigitizer(AliRunDigitizer * manager)
    : AliFMDBaseDigitizer(manager) 
  {}
  /** DTOR */
  virtual ~AliFMDSDigitizer() {}
protected:
  /** Output to disk 
      @param outFMD Loader
      @param fmd    AliFMD object */
  virtual void OutputTree(AliLoader* outFMD, AliFMD* fmd);
  /** Add a digit to output.
      @param fmd      Pointer to detector object
      @param detector Detector #
      @param ring     Ring ID
      @param sector   Sector number
      @param strip    Strip number
      @param edep     Energy deposited (not used)
      @param count1   ADC count 1
      @param count2   ADC count 2 (-1 if not used)
      @param count3   ADC count 3 (-1 if not used) */
  virtual void     AddDigit(AliFMD*  fmd,
			    UShort_t detector, 
			    Char_t   ring,
			    UShort_t sector, 
			    UShort_t strip, 
			    Float_t  edep, 
			    UShort_t count1, 
			    Short_t  count2, 
			    Short_t  count3,
			    Short_t  count4) const;
  ClassDef(AliFMDSDigitizer,0) // Make Summable Digits from Hits
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

