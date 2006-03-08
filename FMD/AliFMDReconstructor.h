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

//____________________________________________________________________
class TTree;
class TClonesArray;
class AliFMDDigit;
class AliRawReader;
class AliRunLoader;
class AliESD;
class AliESDFMD;

//____________________________________________________________________
class AliFMDReconstructor: public AliReconstructor 
{
public:
  AliFMDReconstructor();
  AliFMDReconstructor(const AliFMDReconstructor& other);
  virtual ~AliFMDReconstructor();
  AliFMDReconstructor& operator=(const AliFMDReconstructor& other);

  virtual void   Init(AliRunLoader* runLoader);
  virtual Bool_t HasDigitConversion() const { return kTRUE; }
  virtual void   ConvertDigits(AliRawReader* reader, TTree* digitsTree) const;
  virtual Bool_t HasLocalReconstruction() const { return kTRUE; }
  virtual void   Reconstruct(TTree* digitsTree, TTree* clusterTree) const;
  virtual void   FillESD(TTree* digitsTree, TTree* clusterTree, 
			 AliESD* esd) const;
  virtual void   SetESD(AliESD* esd) { fESD = esd; }
     
private:
  void Reconstruct(AliRawReader*, TTree*) const;
  void Reconstruct(AliRunLoader*) const;
  void Reconstruct(AliRunLoader*, AliRawReader*) const;
  void FillESD(AliRawReader*, TTree*, AliESD*) const;
  void FillESD(AliRunLoader*, AliESD*) const;
  void FillESD(AliRunLoader*, AliRawReader*, AliESD*) const;
  
protected:
  virtual void     ProcessDigits(TClonesArray* digits) const;
  virtual UShort_t SubtractPedestal(AliFMDDigit* digit) const;
  virtual Float_t  Adc2Energy(AliFMDDigit* digit, Float_t eta, 
			      UShort_t count) const;
  virtual Float_t  Energy2Multiplicity(AliFMDDigit* digit, Float_t edep) const;
  virtual void     PhysicalCoordinates(AliFMDDigit* digit, Float_t& eta, 
				       Float_t& phi) const;
  
  mutable TClonesArray* fMult;          // Cache of RecPoints
  mutable Int_t         fNMult;         // Number of entries in fMult 
  mutable TTree*        fTreeR;         // Output tree 
  mutable Float_t       fCurrentVertex; // Z-coordinate of primary vertex
  mutable AliESDFMD*    fESDObj;        // ESD output object
  AliESD*               fESD;
  
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
