// -*- mode: c++ -*- 
//
//  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
//  reserved. 
//
//  See cxx source for full Copyright notice                               
//
//  AliFMDReconstructor.h 
//  Task Class for making TreeR for FMD                        
//
//-- Authors: Evgeny Karpechev (INR) and Alla Maevskaia (INR)
//   Latest changes by Christian Holm Christensen <cholm@nbi.dk>
/*
    Reconstruct nember of particles in given group of pads for given
    FMDvolume determine by numberOfVolume ,
    numberOfMinSector,numberOfMaxSector, numberOfMinRing,
    numberOfMaxRing Reconstruction method choose dependence on number
    of empty pads
  */
/* $Id$ */
#ifndef ALIFMDRECONSTRUCTOR_H
#define ALIFMDRECONSTRUCTOR_H
#ifndef ALIRECONSTRUCTOR_H
# include <AliReconstructor.h>
#endif
#ifndef ALIFMDMAP_H
# include <AliFMDMap.h>
#endif

//____________________________________________________________________
class TClonesArray;
class AliFMD;
class AliLoader;
class AliRunLoader;
class AliFMDDigit;
class AliRawReader;
typedef AliFMDMap<UShort_t> AliFMDAdcMap;


//____________________________________________________________________
class AliFMDReconstructor: public AliReconstructor 
{
protected:
  mutable AliFMDAdcMap  fAdcs;
  mutable AliRunLoader* fRunLoader;
  mutable AliLoader*    fFMDLoader;
  mutable TClonesArray* fParticles;
  mutable AliFMD*       fFMD;
  
  Float_t               fDeltaEta;
  Float_t               fDeltaPhi;
  UShort_t              fThreshold;
  Float_t               fPedestal;
  Float_t               fPedestalWidth;
  mutable Int_t         fEmptyStrips;
  mutable Int_t         fTotalStrips;
  
  enum { 
    kMaxDetectors = 3, 
    kMaxRings     = 2, 
    kMaxSectors   = 20, 
    kMaxStrips    = 512
  };
  
public:
  AliFMDReconstructor();
  virtual ~AliFMDReconstructor() {}

  void         SetDeltaEta(Float_t deta=.1)  { fDeltaEta = deta;  }
  void         SetDeltaPhi(Float_t dphi=360) { fDeltaPhi = dphi;  } 
  void         SetThreshold(UShort_t t=6)    { fThreshold = t; }
  void         SetPedestal(Float_t mean=10, Float_t width=1);
     
  virtual void Reconstruct(AliRunLoader* runLoader) const;
  virtual void Reconstruct(AliRunLoader* runLoader,  
			   AliRawReader* rawReader) const;
  virtual void FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  
protected:
  virtual void     ProcessEvent(Int_t event, 
				AliRawReader* rawReader, 
				TClonesArray* digits) const;
  virtual Bool_t   ReadAdcs(TClonesArray* digits) const;
  virtual Bool_t   ReadAdcs(AliRawReader* rawReader) const;
  virtual void     ProcessDigit(AliFMDDigit* digit) const;
  virtual UShort_t SubtractPedestal(AliFMDDigit* digit) const;
  virtual void     ReconstructFromCache(Float_t zVertex) const;
  ClassDef(AliFMDReconstructor, 0)  // class for the FMD reconstruction
}; 
#endif
//____________________________________________________________________
//
// EOF
//
