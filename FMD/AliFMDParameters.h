#ifndef ALIFMDPARAMETERS_H
#define ALIFMDPARAMETERS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
//____________________________________________________________________
//
//  Singleton class to handle various parameters (not geometry) of the
//  FMD
//  Should get ata fromm Conditions DB.
//
/** @file    AliFMDParameters.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:44:43 2006
    @brief   Manager of FMD parameters
*/
#ifndef ROOT_TNamed
# include <TNamed.h>
#endif
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif
#ifndef ALIFMDUSHORTMAP_H
# include <AliFMDUShortMap.h>
#endif
#ifndef ALIFMDBOOLMAP_H
# include <AliFMDBoolMap.h>
#endif
typedef AliFMDUShortMap AliFMDCalibZeroSuppression;
typedef AliFMDBoolMap   AliFMDCalibDeadMap;
class AliFMDCalibPedestal;
class AliFMDCalibGain;
class AliFMDCalibSampleRate;
class AliFMDCalibStripRange;
class AliFMDAltroMapping;
//____________________________________________________________________
//
//  Singleton class to handle various parameters (not geometry) of the
//  FMD
//  Should get ata fromm Conditions DB.
//

/** @brief This class is a singleton that handles various parameters
    of the FMD detectors.  
    This class reads from the Conditions DB to get the various
    parameters, which code can then request from here. In that way,
    all code uses the same data, and the interface is consistent.
     
    Some of the parameter managed are 
    - @c fPedestal, @c fPedestalWidth
      Mean and width of the pedestal.  The pedestal is simulated
      by a Guassian, but derived classes my override MakePedestal
      to simulate it differently (or pick it up from a database).
    - @c fVA1MipRange
      The dymamic MIP range of the VA1_ALICE pre-amplifier chip 
    - @c fAltroChannelSize
      The largest number plus one that can be stored in one
      channel in one time step in the ALTRO ADC chip. 
    - @c fSampleRate
      How many times the ALTRO ADC chip samples the VA1_ALICE
      pre-amplifier signal.   The VA1_ALICE chip is read-out at
      10MHz, while it's possible to drive the ALTRO chip at
      25MHz.  That means, that the ALTRO chip can have time to
      sample each VA1_ALICE signal up to 2 times.  Although it's
      not certain this feature will be used in the production,
      we'd like have the option, and so it should be reflected in
      the code.

    @ingroup FMD_base
*/
class AliFMDParameters : public TNamed
{
public:
  /** Singleton access
      @return  single to */
  static AliFMDParameters* Instance();

  /** Initialize the manager.  This tries to read the parameters from
      CDB.  If that fails, the class uses the hard-coded parameters. 
   */
  void Init(Bool_t forceReInit=kFALSE);
  /** Print all parameters. 
      @param option Option string */
  void Print(Option_t* option="A") const;
  /** Draw parameters. 
      @param option What to draw. Should be one of 
      - dead	  Dead channels
      - threshold Threshold
      - gain	  Gain
      - pedestal  Pedestal
      - noise	  Noise (or pedestal width)
      - zero	  Zero suppression
      - rate	  Sampling rate (VA1 clock / ALTRO clock)
      - min	  Minimum strip read out
      - max 	  Maximum strip read out
      - map	  hardware address
  */
  void Draw(Option_t* option="pedestal");
  
  /** @{ */
  /** @name Set various `Fixed' parameters */
  /** @param r How many MIP signals we can fit in the VA1
      pre-amps. (default and design is 20) */
  void SetVA1MipRange(UShort_t r=20)          { fVA1MipRange = r; }
  /** @param s Maximum number of the ADC (ALTRO).  This is a 10 bit
      ADC so, the maximum number is 1024 */
  void SetAltroChannelSize(UShort_t s=1024)   { fAltroChannelSize = s;}
  /** @param size The number of strips multiplexed into one ALTRO
      channel. That is, how many strips is connected to one VA1
      pre-amp. */
  void SetChannelsPerAltro(UShort_t size=128) { fChannelsPerAltro = size; }
  /** @param f Factor to use for accepting a signal. */
  void SetPedestalFactor(Float_t f=3)         { fPedestalFactor = f; }
  /** @} */

  /** @{ */
  /** @name Set various variable parameter defaults */
  /** @param s Zero suppression threshold in ADC counts */
  void SetZeroSuppression(UShort_t s=0)       { fFixedZeroSuppression = s; }
  /** @param r How many times we oversample each strip. */
  void SetSampleRate(UShort_t r=1)            { fFixedSampleRate = (r>2?2:r);}
  /** @param p Pedestal value in ADC counts */
  void SetPedestal(Float_t p=10)              { fFixedPedestal = p; }
  /** @param w Pedestal width in ADC counts */
  void SetPedestalWidth(Float_t w=1)          { fFixedPedestalWidth = w; }
  /** @param t Threshold used for 1 MIP acceptance. */
  void SetThreshold(Float_t t=0)              { fFixedThreshold = t; }
  /** Range of strips read out 
      @param min Minimum strip number (0-127). 
      @param max Maximum strip number (0-127). */
  void SetStripRange(UShort_t min=0, UShort_t max=127);
  /** @} */

  /** @{ */
  /** @name Get `Fixed' various parameters */
  /** @return Number of MIP signals that fit inside a VA1 channel  */
  UShort_t GetVA1MipRange()          const { return fVA1MipRange; }
  /** @return The maximum count in the ADC */
  UShort_t GetAltroChannelSize()     const { return fAltroChannelSize; }
  /** @return Number of strips muliplexed into one ADC channel */
  UShort_t GetChannelsPerAltro()     const { return fChannelsPerAltro; }
  /** @return The average energy deposited by one MIP */
  Float_t  GetEdepMip()              const;
  /** @return The factor used of signal acceptance */
  Float_t  GetPedestalFactor()	     const { return fPedestalFactor; }
  /** @} */

  /** @{ */
  /** @name Get variable parameters */
  /** Whether the strip is considered dead
      @param detector Detector # (1-3)
      @param ring     Ring ID ('I' or 'O')
      @param sector   Sector number (0-39)
      @param strip    Strip number (0-511)
      @return @c true if the strip is considered dead, @c false if
      it's OK. */
  Bool_t   IsDead(UShort_t detector, 
		  Char_t ring, 
		  UShort_t sector, 
		  UShort_t strip) const;
  Float_t  GetThreshold() const;
  /** Gain of pre-amp. 
      @param detector Detector # (1-3)
      @param ring     Ring ID ('I' or 'O')
      @param sector   Sector number (0-39)
      @param strip    Strip number (0-511)
      @return Gain of pre-amp.  */
  Float_t  GetPulseGain(UShort_t detector, 
			Char_t ring, 
			UShort_t sector, 
			UShort_t strip) const;
  /** Get mean of pedestal
      @param detector Detector # (1-3)
      @param ring     Ring ID ('I' or 'O')
      @param sector   Sector number (0-39)
      @param strip    Strip number (0-511)
      @return Mean of pedestal */
  Float_t  GetPedestal(UShort_t detector, 
		       Char_t ring, 
		       UShort_t sector, 
		       UShort_t strip) const;
  /** Width of pedestal
      @param detector Detector # (1-3)
      @param ring     Ring ID ('I' or 'O')
      @param sector   Sector number (0-39)
      @param strip    Strip number (0-511)
      @return Width of pedestal */
  Float_t  GetPedestalWidth(UShort_t detector, 
			    Char_t ring, 
			    UShort_t sector, 
			    UShort_t strip) const;
  /** zero suppression threshold (in ADC counts)
      @param detector Detector # (1-3)
      @param ring     Ring ID ('I' or 'O')
      @param sector   Sector number (0-39)
      @param strip    Strip number (0-511)
      @return zero suppression threshold (in ADC counts) */
  UShort_t GetZeroSuppression(UShort_t detector, 
			      Char_t ring, 
			      UShort_t sector, 
			      UShort_t strip) const;
  /** Get the sampling rate
      @param detector Detector # (1-3)
      @param ring     Ring ID ('I' or 'O')
      @param sector   Sector number (0-39)
      @param strip    Strip number (0-511)
      @return The sampling rate */
  UShort_t GetSampleRate(UShort_t detector, 
			 Char_t ring, 
			 UShort_t sector, 
			 UShort_t strip) const;
  /** Get the minimum strip in the read-out range
      @param detector Detector # (1-3)
      @param ring     Ring ID ('I' or 'O')
      @param sector   Sector number (0-39)
      @param strip    Strip number (0-511)
      @return Minimum strip */
  UShort_t GetMinStrip(UShort_t detector, 
		       Char_t ring, 
		       UShort_t sector, 
		       UShort_t strip) const;
  /** Get the maximum strip in the read-out range
      @param detector Detector # (1-3)
      @param ring     Ring ID ('I' or 'O')
      @param sector   Sector number (0-39)
      @param strip    Strip number (0-511)
      @return Maximum strip */
  UShort_t GetMaxStrip(UShort_t detector, 
		       Char_t ring, 
		       UShort_t sector, 
		       UShort_t strip) const;
  /** Translate hardware address to detector coordinates 
      @param ddl      DDL number 
      @param addr     Hardware address
      @param det      On return, Detector # (1-3)
      @param ring     On return, Ring ID ('I' or 'O')
      @param sec      On return, Sector number (0-39)
      @param str      On return, Strip number (0-511)
      @return @c true on success. */
  Bool_t   Hardware2Detector(UInt_t ddl, UInt_t addr, UShort_t& det,
			     Char_t& ring, UShort_t& sec, UShort_t& str) const;
  /** Translate detector coordinates to hardware address 
      @param det      Detector # (1-3)
      @param ring     Ring ID ('I' or 'O')
      @param sec      Sector number (0-39)
      @param str      Strip number (0-511)
      @param ddl      On return, DDL number 
      @param addr     On return, Hardware address
      @return @c true on success. */
  Bool_t   Detector2Hardware(UShort_t det, Char_t ring, UShort_t sec, 
			     UShort_t str, UInt_t& ddl, UInt_t& addr) const;
  /** Get the map that translates hardware to detector coordinates 
      @return Get the map that translates hardware to detector
      coordinates */ 
  AliFMDAltroMapping* GetAltroMap() const;
  /** @} */

  static const char* PulseGainPath()       { return fgkPulseGain; }
  static const char* PedestalPath()        { return fgkPedestal; }
  static const char* DeadPath()            { return fgkDead; }
  static const char* SampleRatePath()      { return fgkSampleRate; }
  static const char* AltroMapPath()        { return fgkAltroMap; }
  static const char* ZeroSuppressionPath() { return fgkZeroSuppression; }
  static const char* StripRangePath()      { return fgkStripRange; }
protected:
  /** CTOR  */
  AliFMDParameters();
  /** CTOR  */
  AliFMDParameters(const AliFMDParameters& o) 
    : TNamed(o), 
      fIsInit(o.fIsInit),
      fkSiDeDxMip(o.fkSiDeDxMip),
      fVA1MipRange(o.fVA1MipRange),
      fAltroChannelSize(o.fAltroChannelSize),
      fChannelsPerAltro(o.fChannelsPerAltro),
      fPedestalFactor(o.fPedestalFactor),
      fFixedPedestal(o.fFixedPedestal),
      fFixedPedestalWidth(o.fFixedPedestalWidth),
      fFixedZeroSuppression(o.fFixedZeroSuppression),
      fFixedSampleRate(o.fFixedSampleRate),
      fFixedThreshold(o.fFixedThreshold),
      fFixedMinStrip(o.fFixedMinStrip),
      fFixedMaxStrip(o.fFixedMaxStrip),
      fFixedPulseGain(o.fFixedPulseGain),
      fEdepMip(o.fEdepMip),
      fZeroSuppression(o.fZeroSuppression),
      fSampleRate(o.fSampleRate),
      fPedestal(o.fPedestal),
      fPulseGain(o.fPulseGain),
      fDeadMap(o.fDeadMap),
      fAltroMap(o.fAltroMap),
      fStripRange(o.fStripRange)
  {}
  /** Assignement operator 
      @return Reference to this */
  AliFMDParameters& operator=(const AliFMDParameters&) { return *this; }
  /** DTOR */
  virtual ~AliFMDParameters() {}
  /** Singleton instance  */
  static AliFMDParameters* fgInstance;   // Static singleton instance
  /** Initialize gains.  Try to get them from CDB */
  void InitPulseGain();
  /** Initialize pedestals.  Try to get them from CDB */
  void InitPedestal();
  /** Initialize dead map.  Try to get it from CDB */
  void InitDeadMap();
  /** Initialize sample rates.  Try to get them from CDB */
  void InitSampleRate();
  /** Initialize zero suppression thresholds.  Try to get them from CDB */
  void InitZeroSuppression();
  /** Initialize hardware map.  Try to get it from CDB */
  void InitAltroMap();
  /** Initialize strip range.  Try to get it from CDB */
  void InitStripRange();

  Bool_t          fIsInit;               // Whether we've been initialised  

  static const char* fgkPulseGain;	 // Path to PulseGain calib object
  static const char* fgkPedestal;	 // Path to Pedestal calib object
  static const char* fgkDead;	         // Path to Dead calib object
  static const char* fgkSampleRate;	 // Path to SampleRate calib object
  static const char* fgkAltroMap;	 // Path to AltroMap calib object
  static const char* fgkZeroSuppression; // Path to ZeroSuppression cal object
  static const char* fgkStripRange;      // Path to strip range cal object
  const Float_t   fkSiDeDxMip;           // MIP dE/dx in Silicon
  UShort_t        fVA1MipRange;          // # MIPs the pre-amp can do    
  UShort_t        fAltroChannelSize;     // Largest # to store in 1 ADC ch.
  UShort_t        fChannelsPerAltro;     // Number of pre-amp. chan/adc chan.
  Float_t         fPedestalFactor;       // Number of pedestal widths

  Float_t         fFixedPedestal;        // Pedestal to subtract
  Float_t         fFixedPedestalWidth;   // Width of pedestal
  UShort_t        fFixedZeroSuppression; // Threshold for zero-suppression
  UShort_t        fFixedSampleRate;      // Times the ALTRO samples pre-amp.
  Float_t         fFixedThreshold;       // Threshold in ADC counts
  UShort_t        fFixedMinStrip;        // Minimum strip read-out
  UShort_t        fFixedMaxStrip;        // Maximum strip read-out 
  mutable Float_t fFixedPulseGain;       //! Gain (cached)
  mutable Float_t fEdepMip;              //! Cache of energy loss for a MIP
  
  AliFMDCalibZeroSuppression* fZeroSuppression; // Zero suppression from CDB
  AliFMDCalibSampleRate*      fSampleRate;      // Sample rate from CDB 
  AliFMDCalibPedestal*        fPedestal;        // Pedestals 
  AliFMDCalibGain*            fPulseGain;       // Pulser gain
  AliFMDCalibDeadMap*         fDeadMap;         // Pulser gain
  AliFMDAltroMapping*         fAltroMap;        // Map of hardware
  AliFMDCalibStripRange*      fStripRange;      // Strip range
  
  ClassDef(AliFMDParameters,5) // Manager of parameters
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//

