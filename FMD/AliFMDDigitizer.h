#ifndef ALIFMDDIGITIZER_H
#define ALIFMDDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIDIGITIZER_H
# include <AliDigitizer.h>
#endif
#ifndef ALIRUNDIGITIZER_H
# include <AliRunDigitizer.h>
#endif
#ifndef ALIFMDMAP_H
# include <AliFMDMap.h>
#endif
#ifndef __UTILITY__
# include <utility>
#endif
// #ifndef ROOT_TArrayF
// # include <TArrayF.h>
// #endif

//====================================================================
class TClonesArray;
class AliFMD;
class AliLoader;
class AliRunLoader;

typedef AliFMDMap<std::pair<Float_t, UShort_t> > AliFMDEdepMap;

//====================================================================
class AliFMDBaseDigitizer : public AliDigitizer 
{
public:
  AliFMDBaseDigitizer();
  AliFMDBaseDigitizer(AliRunDigitizer * manager);
  AliFMDBaseDigitizer(const Char_t* name, const Char_t* title);
  virtual ~AliFMDBaseDigitizer();
   
  // Do the main work
  virtual Bool_t Init();

  // Extra member functions 
  void     SetVA1MipRange(UShort_t r=20) { fVA1MipRange = r; }
  void     SetAltroChannelSize(UShort_t s=1024) { fAltroChannelSize = s;}
  void     SetSampleRate(UShort_t r=1) { fSampleRate = (r>2 ? 2 : r); }
  void     SetShapingTime(Float_t t=10) { fShapingTime = t;  }
  
  UShort_t GetVA1MipRange()      const { return fVA1MipRange; }
  UShort_t GetAltroChannelSize() const { return fAltroChannelSize; }
  UShort_t GetSampleRate()       const { return fSampleRate; }
  Float_t  GetShapingTime()      const { return fShapingTime; }
protected:
  virtual void     SumContributions(AliFMD* fmd);
  virtual void     DigitizeHits(AliFMD* fmd) const;
  virtual void     ConvertToCount(Float_t   edep, 
				  Float_t   siThickness, 
				  Float_t   siDensity, 
				  TArrayI&  counts) const;
  virtual UShort_t MakePedestal() const { return 0; }
  virtual Float_t  ShapeIntegral(Float_t u, Float_t v) const;
  virtual void     AddNoise(TArrayI&) const {}
  virtual void     AddDigit(AliFMD*  /* fmd      */,
			    UShort_t /* detector */, 
			    Char_t   /* ring     */,
			    UShort_t /* sector   */, 
			    UShort_t /* strip    */, 
			    Float_t  /* edep     */, 
			    UShort_t /* count1   */, 
			    Short_t  /* count2   */, 
			    Short_t  /* count3   */) const {}

  AliRunLoader* fRunLoader;
  AliFMDEdepMap fEdep;             // Cache of Energy from hits 
  UShort_t      fVA1MipRange;      // How many MIPs the pre-amp can do    
  UShort_t      fAltroChannelSize; // Largest # to store in 1 ADC chan.
  UShort_t      fSampleRate;       // Times the ALTRO samples pre-amp.
  Float_t       fShapingTime;      // Shaping profile parameter
  
  enum { 
    kMaxDetectors = 3, 
    kMaxRings     = 2, 
    kMaxSectors   = 20, 
    kMaxStrips    = 512
  };
  ClassDef(AliFMDBaseDigitizer,0) // Base class for FMD digitizers
};

//====================================================================
class AliFMDDigitizer : public AliFMDBaseDigitizer 
{
public:
  AliFMDDigitizer();
  AliFMDDigitizer(AliRunDigitizer * manager);
  virtual ~AliFMDDigitizer() {}
  virtual void  Exec(Option_t* option=0);
  
   
  // Extra member functions 
  void     SetPedestal(Float_t mean=10, Float_t width=.5);
  void     GetPedestal(Float_t& mean,   Float_t& width) const;
protected:
  virtual void     AddDigit(AliFMD*  fmd,
			    UShort_t detector, 
			    Char_t   ring,
			    UShort_t sector, 
			    UShort_t strip, 
			    Float_t  edep, 
			    UShort_t count1, 
			    Short_t  count2, 
			    Short_t  count3) const;
  virtual UShort_t MakePedestal() const;
  virtual void     CheckDigit(Float_t         edep, 
			      UShort_t        nhits,
			      const TArrayI&  counts);
  Float_t       fPedestal;         // Mean of pedestal 
  Float_t       fPedestalWidth;    // Width of pedestal 
  ClassDef(AliFMDDigitizer,0) // Make Digits from Hits
};
//____________________________________________________________________
inline void 
AliFMDDigitizer::SetPedestal(Float_t mean, Float_t width) 
{
  fPedestal      = mean;
  fPedestalWidth = width;
}

//____________________________________________________________________
inline void 
AliFMDDigitizer::GetPedestal(Float_t& mean, Float_t& width)  const
{
  mean  = fPedestal;
  width = fPedestalWidth;
}


//====================================================================
class AliFMDSDigitizer : public AliFMDBaseDigitizer 
{
public:
  AliFMDSDigitizer();
  AliFMDSDigitizer(const Char_t* headerFile, const Char_t* sdigFile="");
  virtual ~AliFMDSDigitizer();
  virtual void  Exec(Option_t* option=0);
protected:
  virtual void     AddDigit(AliFMD*  fmd,
			    UShort_t detector, 
			    Char_t   ring,
			    UShort_t sector, 
			    UShort_t strip, 
			    Float_t  edep, 
			    UShort_t count1, 
			    Short_t  count2, 
			    Short_t  count3) const;
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

