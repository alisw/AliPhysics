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
//
#ifndef ROOT_TNamed
# include <TNamed.h>
#endif

class AliFMDParameters : public TNamed
{
public:
  static AliFMDParameters* Instance();
  
  // Set various parameters 
  void     SetVA1MipRange(UShort_t r=20)          { fVA1MipRange = r; }
  void     SetAltroChannelSize(UShort_t s=1024)   { fAltroChannelSize = s;}
  void     SetChannelsPerAltro(UShort_t size=128) { fChannelsPerAltro = size; }
  void     SetZeroSuppression(UShort_t s=0)       { fZeroSuppression = s; }
  void     SetSampleRate(UShort_t r=1)            { fSampleRate = (r>2?2:r);}
  void     SetPedestal(Float_t p=10)              { fPedestal = p; }
  void     SetPedestalWidth(Float_t w=1)          { fPedestalWidth = w; }
  void     SetPedestalFactor(Float_t f=3)         { fPedestalFactor = f; }

  // Get various parameters
  UShort_t GetVA1MipRange()          const { return fVA1MipRange; }
  UShort_t GetAltroChannelSize()     const { return fAltroChannelSize; }
  UShort_t GetChannelsPerAltro()     const { return fChannelsPerAltro; }
  UShort_t GetZeroSuppression()      const { return fZeroSuppression; }
  UShort_t GetSampleRate()           const { return fSampleRate; }
  Float_t  GetEdepMip()              const;
  Float_t  GetPedestal()	     const { return fPedestal; }
  Float_t  GetPedestalWidth()	     const { return fPedestalWidth; }
  Float_t  GetPedestalFactor()	     const { return fPedestalFactor; }

  enum { 
    kBaseDDL = 0x1000 // DDL offset for the FMD
  };
protected:
  AliFMDParameters();
  virtual ~AliFMDParameters() {}
  static AliFMDParameters* fgInstance;
  
  const Float_t fSiDeDxMip;        // MIP dE/dx in Silicon
  UShort_t      fVA1MipRange;      // # MIPs the pre-amp can do    
  UShort_t      fAltroChannelSize; // Largest # to store in 1 ADC ch.
  UShort_t      fChannelsPerAltro; // Number of pre-amp. channels/adc chan.
  UShort_t      fZeroSuppression;  // Threshold for zero-suppression
  UShort_t      fSampleRate;       // Times the ALTRO samples pre-amp.
  Float_t       fPedestal;         // Pedestal to subtract
  Float_t       fPedestalWidth;    // Width of pedestal
  Float_t       fPedestalFactor;   // Number of pedestal widths

  
  ClassDef(AliFMDParameters,1)
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

