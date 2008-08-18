#ifndef ALIFMDSSDIGITIZER_H
#define ALIFMDSSDIGITIZER_H
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
/** @file    AliFMDSSDigitizer.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:38:26 2006
    @brief   FMD Digitizers declaration
    @ingroup FMD_sim
*/
#ifndef ALIFMDDIGITIZER_H
# include <AliFMDDigitizer.h>
#endif

//====================================================================
class TClonesArray;
class AliFMD;
class AliLoader;
class AliRunLoader;
class AliFMDDigit;



//====================================================================
/** @class AliFMDSSDigitizer
    @brief Concrete digitizer to make digits from hits.  See also
    AliFMDBaseDigitizer documentation.  
    @ingroup FMD_sim
 */
class AliFMDSSDigitizer : public AliFMDDigitizer 
{
public:
  /** CTOR */
  AliFMDSSDigitizer() : AliFMDDigitizer() {}
  /** CTOR 
      @param manager Manager of digitization */
  AliFMDSSDigitizer(AliRunDigitizer * manager)
    : AliFMDDigitizer(manager) {}
  /** DTOR */
  virtual ~AliFMDSSDigitizer() {}
protected:
  /** Sum contributions from SDigits 
      @param fmd Pointr to loaded AliFMD  */
  void SumContributions(AliFMD* fmd);
  
  ClassDef(AliFMDSSDigitizer,1) // Make Digits from Hits
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

