#ifndef ALIFMDHITDIGITIZER_H
#define ALIFMDHITDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
// Classses to make Hits into digits and summable digits
//    
//    Digits consists of
//    - Detector #
//    - Ring ID                                             
//    - Sector #     
//    - Strip #
//    - ADC count in this channel
//
/** @file    AliFMDHitDigitizer.h
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
class AliStack;



//====================================================================
/** @class AliFMDHitDigitizer
    @brief Concrete digitizer to make digits from hits.  See also
    AliFMDBaseDigitizer documentation.  
    @ingroup FMD_sim
 */
class AliFMDHitDigitizer : public AliFMDBaseDigitizer
{
public:
  enum Output_t { 
    kDigits, 
    kSDigits
  };
    
  /** CTOR */
  AliFMDHitDigitizer() 
    : AliFMDBaseDigitizer(), 
      fOutput(kDigits), 
      fHoldTime(2e-6),
      fStack(0)
  {}
  /** CTOR 
      @param name Name */
  AliFMDHitDigitizer(AliFMD* fmd, Output_t  output);
  /** DTOR */
  virtual ~AliFMDHitDigitizer() {}
  /** Run over the input events (retrieved via run loader) */
  void Digitize(Option_t* option="");
  /** 
   * Set the end of integration
   * 
   * @param holdT Time when integration ends (nominally @f$
   *        2\mu{}s@f$) 
   */
  void SetHoldTime(Double_t holdT=2e-6) { fHoldTime = holdT; }
  /** 
   * Get the hold time 
   * 
   * @return Hold time in seconds
   */
  Double_t GetHoldTime() const { return fHoldTime; }
protected:
  /** Copy constructor 
      @param o Object to copy from */
  AliFMDHitDigitizer(const AliFMDHitDigitizer& o) 
    : AliFMDBaseDigitizer(o),
      fOutput(o.fOutput), 
      fHoldTime(2e-6),
      fStack(o.fStack)
  {}
  /** 
   * Assignment operator
   *
   * @param o Object to assign from 
   * @return Reference to this 
   */
  AliFMDHitDigitizer& operator=(const AliFMDHitDigitizer& o); 
  /** 
   * Make the output tree using the passed loader 
   *
   * @param loader 
   * @return The generated tree. 
   */
  TTree* MakeOutputTree(AliLoader* loader);
  /** Sum energy deposited contributions from each hit in a cache
      @param hitsBranch Branch in input tree */
  void SumContributions(TBranch* hitsBranch);
  /** Make a pedestal 
      @param detector Detector #
      @param ring     Ring ID
      @param sector   Sector #
      @param strip    Strip #
      @return Pedestal value */
  UShort_t MakePedestal(UShort_t  detector, 
			Char_t    ring, 
			UShort_t  sector, 
			UShort_t  strip) const;
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
  void AddDigit(UShort_t       detector, 
		Char_t         ring,
		UShort_t       sector, 
		UShort_t       strip, 
		Float_t        edep, 
		UShort_t       count1, 
		Short_t        count2, 
		Short_t        count3,
		Short_t        count4, 
		UShort_t       ntot, 
		UShort_t       nprim,
		const TArrayI& trackrefs) const;
  /** Check that digit data is consistent
      @param digit   Digit
      @param nhits   Number of hits
      @param counts  ADC counts */
  void CheckDigit(AliFMDDigit*    digit,
		  UShort_t        nhits,
		  const TArrayI&  counts);
  /** 
   * Store the data using the loader 
   *
   * @param loader The loader 
   */
  void StoreDigits(const AliLoader* loader);
  

  Output_t      fOutput;           // Output mode
  Double_t      fHoldTime;         // Stop of integration
  AliStack*     fStack;            // Kinematics

  ClassDef(AliFMDHitDigitizer,1) // Make Digits from Hits
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

