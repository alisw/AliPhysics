#ifndef ALIFMDRECONSTRUCTOR_H
#define ALIFMDRECONSTRUCTOR_H
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

//____________________________________________________________________
//
// Class to do reconstruction of events based on the FMD data.  The
// class will do two kinds of reconstruction, one based on energy
// deposition, and one using hit patterns. 
//

// Header guards in the header files speeds up the compilation
// considerably.  Please leave them in. 
#ifndef ALIRECONSTRUCTOR_H
# include <AliReconstructor.h>
#endif
#ifndef ROOT_TObjArray
# include <TObjArray.h>
#endif
// #ifndef ALIFMDUSHORTMAP_H
// # include <AliFMDUShortMap.h>
// #endif

//____________________________________________________________________
class TClonesArray;
class AliFMD;
class AliLoader;
class AliRunLoader;
class AliFMDDigit;
class AliRawReader;

// typedef AliFMDUShortMap AliFMDAdcMap;


//____________________________________________________________________
class AliFMDReconstructor: public AliReconstructor 
{
public:
  AliFMDReconstructor();
  AliFMDReconstructor(const AliFMDReconstructor& other);
  virtual ~AliFMDReconstructor() {}
  AliFMDReconstructor& operator=(const AliFMDReconstructor& other);

  void         SetDeltaEta(Float_t deta=.1)  { fDeltaEta = deta;  }
  void         SetDeltaPhi(Float_t dphi=360) { fDeltaPhi = dphi;  } 
  void         SetThreshold(UShort_t t=6)    { fThreshold = t; }
  void         SetPedestal(Float_t mean=10, Float_t width=1, Float_t f=3);
     
  virtual void Reconstruct(AliRunLoader* runLoader) const;
  virtual void Reconstruct(AliRunLoader* runLoader,  
			   AliRawReader* rawReader) const;
  virtual void FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  
protected:
  virtual void     ProcessEvent(Int_t event, 
				AliRawReader* rawReader) const;
  virtual void     ProcessDigits(TClonesArray* digits) const;
  virtual UShort_t SubtractPedestal(AliFMDDigit* digit) const;
  virtual void     ReconstructFromCache() const;

  mutable AliRunLoader* fRunLoader;  // Run loader 
  mutable AliLoader*    fFMDLoader;  // FMD specific loader 
  mutable TClonesArray* fParticles;  // Array of particles 
  mutable AliFMD*       fFMD;        // Pointer to FMD manager 
  
  TObjArray	        fAlgorithms;    // Array of algorithms
  Float_t               fDeltaEta;      // Bin size in eta
  Float_t               fDeltaPhi;      // Bin size in phi
  UShort_t              fThreshold;     // Threshold for Poisson recon.
  Float_t               fPedestal;      // Pedestal to subtract
  Float_t               fPedestalWidth; // Width of pedestal
  Float_t               fPedestalFactor;// Number of pedestal widths 
  mutable Int_t         fEmptyStrips;   // Number of empty strips
  mutable Int_t         fTotalStrips;   // Total number of strips 
  mutable Float_t       fCurrentVertex; // Z-coordinate of primary vertex
  
  enum { 
    kMaxDetectors = 3,                  // Maximum number of sub-det.
    kMaxRings     = 2,                  // Maximum number of rings
    kMaxSectors   = 40,                 // Maximum number of sectors
    kMaxStrips    = 512                 // Maximum number of strips
  };
  
  ClassDef(AliFMDReconstructor, 0)  // class for the FMD reconstruction
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
