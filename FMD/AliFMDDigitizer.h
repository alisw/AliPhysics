#ifndef ALIFMDDIGITIZER_H
#define ALIFMDDIGITIZER_H
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
/** @file    AliFMDDigitizer.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:38:26 2006
    @brief   FMD Digitizers declaration
    @ingroup FMD_sim
*/
#ifndef ALIFMDBASEDIGITIZER_H
# include <AliFMDBaseDigitizer.h>
#endif

//====================================================================
class TClonesArray;
class AliFMD;
class AliLoader;
class AliRunLoader;
class AliFMDDigit;



//====================================================================
/** @class AliFMDDigitizer
    @brief Concrete digitizer to make digits from hits.  See also
    AliFMDBaseDigitizer documentation.  
    @ingroup FMD_sim
 */
class AliFMDDigitizer : public AliFMDBaseDigitizer 
{
public:
  /** CTOR */
  AliFMDDigitizer() : AliFMDBaseDigitizer() {}
  /** CTOR 
      @param manager Manager of digitization */
  AliFMDDigitizer(AliRunDigitizer * manager)
    : AliFMDBaseDigitizer(manager) {}
  /** DTOR */
  virtual ~AliFMDDigitizer() {}
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
      @param count3   ADC count 3 (-1 if not used) 
      @param count4   ADC count 4 (-1 if not used) */
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
  /** MAke a pedestal
      @param detector Detector #
      @param ring     Ring ID
      @param sector   Sector number
      @param strip    Strip number
      @return Random noise */
  virtual UShort_t MakePedestal(UShort_t  detector, 
				Char_t    ring, 
				UShort_t  sector, 
				UShort_t  strip) const;
  /** Check that digit data is consistent
      @param digit   Digit
      @param nhits   Number of hits
      @param counts  ADC counts */
  virtual void     CheckDigit(AliFMDDigit*    digit,
			      UShort_t        nhits,
			      const TArrayI&  counts);
  ClassDef(AliFMDDigitizer,1) // Make Digits from Hits
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

