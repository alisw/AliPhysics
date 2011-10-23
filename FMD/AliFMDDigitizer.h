#ifndef ALIFMDDIGITIZER_H
#define ALIFMDDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
// Classses to make SDigits into Digits 
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
/** 
 * @class AliFMDDigitizer
 * @brief Concrete digitizer to make digits from hits.  See also
 * AliFMDBaseDigitizer documentation.  
 *
 * @ingroup FMD_sim
 */
class AliFMDDigitizer : public AliFMDBaseDigitizer 
{
public:
  /** 
   * CTOR 
   */
  AliFMDDigitizer() : AliFMDBaseDigitizer() {}
  /** 
   * CTOR 
   *
   * @param input Input of digitization 
   */
  AliFMDDigitizer(AliDigitizationInput * digInput)
    : AliFMDBaseDigitizer(digInput) {}
  /** 
   * DTOR 
   */
  virtual ~AliFMDDigitizer() {}
  /** 
   * Initialise 
   */
  virtual Bool_t Init();
  /** 
   * Execute this digitizer.  
   * This member function will be called once per event by the passed
   * AliDigitizationInput* digInput object. 
   *
   * @param options Not used 
   */
  virtual void Digitize(Option_t* option="");
protected:
  /** 
   * Sum contributions from SDigits 
   *
   * @param sdigitsBranch Branch of SDigit data 
   */
  void SumContributions(TBranch* sdigitsBranch);
  
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

