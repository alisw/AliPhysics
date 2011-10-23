#ifndef ALIFMDBASEDIGITIZER_H
#define ALIFMDBASEDIGITIZER_H
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
/** @file    AliFMDBaseDigitizer.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:38:26 2006
    @brief   FMD Digitizers declaration
    @ingroup FMD_sim
*/
#ifndef ALIDIGITIZER_H
# include <AliDigitizer.h>
#endif
#ifndef ALIRUNDIGITIZER_H
# include <AliDigitizationInput.h>
#endif
#ifndef ALIFMDEdepMAP_H
# include "AliFMDEdepMap.h"
#endif
//====================================================================
class TClonesArray;
class AliFMD;
class AliLoader;
class AliRunLoader;
class AliFMDDigit;


//====================================================================
/** @class AliFMDBaseDigitizer AliFMDDigitizer.h <FMD/AliFMDDigitizer>
    @brief Base class for digitizers.

    This class contains the procedures simulation ADC  signal for the
    Forward Multiplicity detector  : Hits->Digits and Hits->SDigits
    
    Digits consists of
    - Detector #
    - Ring ID                                             
    - Sector #     
    - Strip #
    - ADC count in this channel                                  

    Summable digits consists of	
    - Detector #
    - Ring ID                                             
    - Sector #     
    - Strip #
    - Total energy deposited in the strip
    - ADC count in this channel                                  

    As the Digits and SDigits have so much in common, the classes
    AliFMDDigitizer and AliFMDSDigitizer are implemented via a base
    class AliFMDBaseDigitizer.
    @verbatim
                    +---------------------+
                    | AliFMDBaseDigitizer |
                    +---------------------+
                              ^
                              |
                   +----------+---------+
                   |                    |
         +-----------------+     +------------------+
         | AliFMDDigitizer |	| AliFMDSDigitizer |
         +-----------------+	+------------------+
    @endverbatim
    These classes uses parameters fetched from the AliFMDParameters
    manager. 

    The shaping function of the VA1 is generally given by 
    @f[
    f(x) = A(1 - \exp(-Bx))
    @f]
    where A is the total charge collected in the pre-amp., and B is a
    paramter that depends on the shaping time of the VA1 circut.
    
    When simulating the shaping function of the VA1 pre-amp. chip, we
    have to take into account, that the shaping function depends on
    the previous value of read from the pre-amp.  

    That results in the following algorithm:
    @code
    last = 0;
    for (i=0; i < n_pre_amp_charge; i++) {
      charge = GetCharge(i);
      if (last < charge) 
        f(t) = (charge - last) * (1 - exp(-B * t)) + last
      else
        f(t) = (last - charge) * exp(-B * t) + charge)
      for (j=0; j < sample_rate; j++) 
        adc[j] = f(i / (# samples))
      last = charge
    }
    @endcode
    Here, the first loop is over all charges collected by the VA1
    chip, and the @c sample_rate is how many times the ALTRO ADC
    samples each of the 128  charges from the pre-amp. 

    The @c charge is the total charge @f$ Q@f$ collected by the VA1
    pre-amplifier for a strip.   @f$ Q@f$ is then given by 
    @f[
    Q = \frac{E}{e}\frac{S}{r}
    @f]
    where @f$ E@f$ is the total energy deposited in a silicon strip, 
    @f$ R@f$ is the dynamic range of the VA1 pre-amp, @f$ e@f$ is the 
    energy deposited by a single MIP, and @f$ S@f$ ALTRO channel size
    in each time step.

    The energy deposited per MIP is given by 
    @f$ 
    e = M \rho w 
    @f$
    where @f$ M@f$ is the universal number 
    @f$ 1.664 \mbox{keV}\mbox{cm}^{2}\mbox{g}^{-1}@f$, @f$ \rho@f$ is
    the density of silicon, and @f$ w@f$ is the depth of the silicon
    sensor.  

    The final ADC count is given by 
    @f[
    C' = C + P
    @f]
    where @f$ P@f$ is the (randomized) pedestal.

    This class uses the class template AliFMDEdepMap to make an
    internal cache of the energy deposted of the hits.  The class
    template is instantasized as 

    The first member of the values is the summed energy deposition in a
    given strip, while the second member of the values is the number of
    hits in a given strip.  Using the second member, it's possible to
    do some checks on just how many times a strip got hit, and what
    kind of error we get in our reconstructed hits.  Note, that this
    information is currently not written to the digits tree.  I think a
    QA (Quality Assurance) digit tree is better suited for that task.
    However, the information is there to be used in the future. 
    @ingroup FMD_sim
 */
class AliFMDBaseDigitizer : public AliDigitizer 
{
public:
  /** CTOR */
  AliFMDBaseDigitizer();
  /** Normal CTOR 
      @param manager Manager of digitization */
  AliFMDBaseDigitizer(AliDigitizationInput * digInp);
  /** Normal ctor 
      @param name Name 
      @param title Title */
  AliFMDBaseDigitizer(const Char_t* name, const Char_t* title);
  /** DTOR */
  virtual ~AliFMDBaseDigitizer();
   
  /** Initialize */
  virtual Bool_t Init();
  
  /** The response shape of the VA1 shaping circuit is approximently
      given by 
      @f[
      f(x) = A(1 - \exp(-Bx))
      @f]
      where @f$ A@f$ is the total charge collected by the pre-amp.,
      and @f$ B@f$ is parameter that depends on the shaping time of
      the @b VA1 pre-amp.  This member function sets the parameter @f$
      B@f$ 
      @param B */
  void     SetShapingTime(Float_t B=10) { fShapingTime = B;  }  
  /** @return Get the shaping time */
  Float_t  GetShapingTime()      const { return fShapingTime; }
  
  void SetStoreTrackRefs(Bool_t store=kTRUE) { fStoreTrackRefs = store; }
  Bool_t IsStoreTrackRefs() const { return fStoreTrackRefs; }
    
protected:
  /** For the stored energy contributions in the cache, convert the
      energy signal to ADC counts, and store the created digit in  
      the digits array
      @param fmd Pointer to detector */
  virtual void     DigitizeHits() const;
  /** Convert the total energy deposited to a (set of) ADC count(s).
      See also the class description for more details. 
      @param edep     Total energy deposited in detector
      @param last     Last charge collected in previous VA1 channnel
      @param detector Detector #
      @param ring     Ring ID
      @param sector   Sector #
      @param strip    Strip #
      @param counts   Array holding the counts on return */
  virtual void     ConvertToCount(Float_t   edep, 
				  Float_t   last,
				  UShort_t  detector, 
				  Char_t    ring, 
				  UShort_t  sector, 
				  UShort_t  strip,
				  TArrayI&  counts) const;
  /** Make a pedestal 
      @param detector Detector #
      @param ring     Ring ID
      @param sector   Sector #
      @param strip    Strip #
      @return Pedestal value */
  virtual UShort_t MakePedestal(UShort_t  detector, 
				Char_t    ring, 
				UShort_t  sector, 
				UShort_t  strip) const;
  /** Add noise to each sample */
  virtual void     AddNoise(TArrayI&) const {}

  /** Add edep contribution from (detector,ring,sector,strip) to cache */ 
  virtual void AddContribution(UShort_t detector, 
			       Char_t   ring, 
			       UShort_t sector, 
			       UShort_t strip, 
			       Float_t  edep, 
			       Bool_t   isPrimary,
			       Int_t    nTrackno,
			       Int_t*   tracknos);
  /** Add a digit to output */
  virtual void     AddDigit(UShort_t       detector, 
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
			    const TArrayI& refs) const;
  /** Make the output tree using the passed loader 
      @param loader 
      @return The generated tree. */
  virtual TTree* MakeOutputTree(AliLoader* loader);
  /** Store the data using the loader 
      @param loader The loader */
  virtual void StoreDigits(const AliLoader* loader);

  AliFMD*         fFMD;              // Detector object 
  AliRunLoader*   fRunLoader;	     //! Run loader
  AliFMDEdepMap   fEdep;             // Cache of Energy from hits 
  Float_t         fShapingTime;      // Shaping profile parameter
  Bool_t          fStoreTrackRefs;   // Wether to store track references
  mutable Int_t   fIgnoredLabels;    //! Number of labels not assigned 
  
  /** Copy CTOR 
      @param o object to copy from  */
  AliFMDBaseDigitizer(const AliFMDBaseDigitizer& o) 
    : AliDigitizer(o),
      fFMD(o.fFMD),
      fRunLoader(0),
      fEdep(o.fEdep),
      fShapingTime(o.fShapingTime),
      fStoreTrackRefs(o.fStoreTrackRefs), 
      fIgnoredLabels(o.fIgnoredLabels)
  {}
  /** 
   * Assignment operator
   * 
   * @return Reference to this object 
   */
  AliFMDBaseDigitizer& operator=(const AliFMDBaseDigitizer& o);

  ClassDef(AliFMDBaseDigitizer,5) // Base class for FMD digitizers
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

