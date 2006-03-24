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
class AliFMDAltroMapping;

/** This class is a singleton that handles various parameters of the
    FMD detectors.  This class reads from the Conditions DB to get the
    various parameters, which code can then request from here. In that
    way, all code uses the same data, and the interface is consistent.
     
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
  void Init();
  
  /** @{ */
  /** @name Set various `Fixed' parameters */
  void SetVA1MipRange(UShort_t r=20)          { fVA1MipRange = r; }
  void SetAltroChannelSize(UShort_t s=1024)   { fAltroChannelSize = s;}
  void SetChannelsPerAltro(UShort_t size=128) { fChannelsPerAltro = size; }
  void SetPedestalFactor(Float_t f=3)         { fPedestalFactor = f; }
  /** @} */

  /** @{ */
  /** @name Set various variable parameter defaults */
  void SetZeroSuppression(UShort_t s=0)       { fFixedZeroSuppression = s; }
  void SetSampleRate(UShort_t r=1)            { fFixedSampleRate = (r>2?2:r);}
  void SetPedestal(Float_t p=10)              { fFixedPedestal = p; }
  void SetPedestalWidth(Float_t w=1)          { fFixedPedestalWidth = w; }
  void SetThreshold(Float_t t=0)              { fFixedThreshold = t; }
  /** @} */

  /** @{ */
  /** @name Get `Fixed' various parameters */
  UShort_t GetVA1MipRange()          const { return fVA1MipRange; }
  UShort_t GetAltroChannelSize()     const { return fAltroChannelSize; }
  UShort_t GetChannelsPerAltro()     const { return fChannelsPerAltro; }
  Float_t  GetEdepMip()              const;
  Float_t  GetPedestalFactor()	     const { return fPedestalFactor; }
  /** @} */

  /** @{ */
  /** @name Get variable parameters */
  Bool_t   IsDead(UShort_t detector, 
		  Char_t ring, 
		  UShort_t sector, 
		  UShort_t strip) const;
  Float_t  GetThreshold() const;
  Float_t  GetPulseGain(UShort_t detector, 
			Char_t ring, 
			UShort_t sector, 
			UShort_t strip) const;
  Float_t  GetPedestal(UShort_t detector, 
		       Char_t ring, 
		       UShort_t sector, 
		       UShort_t strip) const;
  Float_t  GetPedestalWidth(UShort_t detector, 
			    Char_t ring, 
			    UShort_t sector, 
			    UShort_t strip) const;
  UShort_t GetZeroSuppression(UShort_t detector, 
			      Char_t ring, 
			      UShort_t sector, 
			      UShort_t strip) const;
  UShort_t GetSampleRate(UShort_t ddl) const;

  Bool_t   Hardware2Detector(UInt_t ddl, UInt_t addr, UShort_t& det,
			     Char_t& ring, UShort_t& sec, UShort_t& str) const;
  Bool_t   Detector2Hardware(UShort_t det, Char_t ring, UShort_t sec, 
			     UShort_t str, UInt_t& ddl, UInt_t& addr) const;
  AliFMDAltroMapping* GetAltroMap() const;
  /** @} */

  enum { 
    kBaseDDL = 0x1000 // DDL offset for the FMD
  };
  static const char* fgkPulseGain;	 // Path to PulseGain calib object
  static const char* fgkPedestal;	 // Path to Pedestal calib object
  static const char* fgkDead;	         // Path to Dead calib object
  static const char* fgkSampleRate;	 // Path to SampleRate calib object
  static const char* fgkAltroMap;	 // Path to AltroMap calib object
  static const char* fgkZeroSuppression; // Path to ZeroSuppression cal object
protected:
  /** CTOR  */
  AliFMDParameters();
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

  Bool_t          fIsInit;               // Whether we've been initialised  

  const Float_t   fSiDeDxMip;            // MIP dE/dx in Silicon
  UShort_t        fVA1MipRange;          // # MIPs the pre-amp can do    
  UShort_t        fAltroChannelSize;     // Largest # to store in 1 ADC ch.
  UShort_t        fChannelsPerAltro;     // Number of pre-amp. chan/adc chan.
  Float_t         fPedestalFactor;       // Number of pedestal widths

  Float_t         fFixedPedestal;        // Pedestal to subtract
  Float_t         fFixedPedestalWidth;   // Width of pedestal
  UShort_t        fFixedZeroSuppression; // Threshold for zero-suppression
  UShort_t        fFixedSampleRate;      // Times the ALTRO samples pre-amp.
  Float_t         fFixedThreshold;       //
  mutable Float_t fFixedPulseGain;       //! Gain (cached)
  mutable Float_t fEdepMip;              //! Cache of energy loss for a MIP
  
  AliFMDCalibZeroSuppression* fZeroSuppression; // Zero suppression from CDB
  AliFMDCalibSampleRate*      fSampleRate;      // Sample rate from CDB 
  AliFMDCalibPedestal*        fPedestal;        // Pedestals 
  AliFMDCalibGain*            fPulseGain;       // Pulser gain
  AliFMDCalibDeadMap*         fDeadMap;         // Pulser gain
  AliFMDAltroMapping*         fAltroMap;        // Map of hardware
  
  ClassDef(AliFMDParameters,3)
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

