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

class AliFMDParameters : public TNamed
{
public:
  static AliFMDParameters* Instance();

  void Init();
  
  // Set various `Fixed' parameters 
  void SetVA1MipRange(UShort_t r=20)          { fVA1MipRange = r; }
  void SetAltroChannelSize(UShort_t s=1024)   { fAltroChannelSize = s;}
  void SetChannelsPerAltro(UShort_t size=128) { fChannelsPerAltro = size; }
  void SetPedestalFactor(Float_t f=3)         { fPedestalFactor = f; }

  // Set various variable parameter defaults
  void SetZeroSuppression(UShort_t s=0)       { fFixedZeroSuppression = s; }
  void SetSampleRate(UShort_t r=1)            { fFixedSampleRate = (r>2?2:r);}
  void SetPedestal(Float_t p=10)              { fFixedPedestal = p; }
  void SetPedestalWidth(Float_t w=1)          { fFixedPedestalWidth = w; }
  void SetThreshold(Float_t t=0)              { fFixedThreshold = t; }

  // Get `Fixed' various parameters
  UShort_t GetVA1MipRange()          const { return fVA1MipRange; }
  UShort_t GetAltroChannelSize()     const { return fAltroChannelSize; }
  UShort_t GetChannelsPerAltro()     const { return fChannelsPerAltro; }
  Float_t  GetEdepMip()              const;
  Float_t  GetPedestalFactor()	     const { return fPedestalFactor; }

  // Get variable parameters 
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
  AliFMDParameters();
  virtual ~AliFMDParameters() {}
  static AliFMDParameters* fgInstance;   // Static singleton instance

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

