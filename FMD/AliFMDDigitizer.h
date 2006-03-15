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
#ifndef ALIFMDEdepMAP_H
# include "AliFMDEdepMap.h"
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
class AliFMDDigit;


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
  void     SetShapingTime(Float_t t=10) { fShapingTime = t;  }  
  Float_t  GetShapingTime()      const { return fShapingTime; }
protected:
  virtual void     SumContributions(AliFMD* fmd);
  virtual void     DigitizeHits(AliFMD* fmd) const;
  virtual void     ConvertToCount(Float_t   edep, 
				  Float_t   last,
				  UShort_t  detector, 
				  Char_t    ring, 
				  UShort_t  sector, 
				  UShort_t  strip,
				  TArrayI&  counts) const;
  virtual UShort_t MakePedestal(UShort_t  detector, 
				Char_t    ring, 
				UShort_t  sector, 
				UShort_t  strip) const;
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

  AliRunLoader* fRunLoader;	   //! Run loader
  AliFMDEdepMap fEdep;             // Cache of Energy from hits 
  Float_t       fShapingTime;      // Shaping profile parameter
  
  AliFMDBaseDigitizer(const AliFMDBaseDigitizer& o) 
    : AliDigitizer(o) {}
  AliFMDBaseDigitizer& operator=(const AliFMDBaseDigitizer&) { return *this; }
  ClassDef(AliFMDBaseDigitizer,2) // Base class for FMD digitizers
};

//====================================================================
class AliFMDDigitizer : public AliFMDBaseDigitizer 
{
public:
  AliFMDDigitizer();
  AliFMDDigitizer(AliRunDigitizer * manager);
  virtual ~AliFMDDigitizer() {}
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
  virtual UShort_t MakePedestal(UShort_t  detector, 
				Char_t    ring, 
				UShort_t  sector, 
				UShort_t  strip) const;
  virtual void     CheckDigit(AliFMDDigit*    digit,
			      UShort_t        nhits,
			      const TArrayI&  counts);
  ClassDef(AliFMDDigitizer,1) // Make Digits from Hits
};

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

