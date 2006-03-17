#ifndef ALIFMDRECPOINT_H
#define ALIFMDRECPOINT_H

// Reconstracted Particles Class: has number of reconstructed
// particles in sectors from NumOfMinSector to NumberOfMaxSector()
// rings from NumOfMinRing to NumOfMaxRing for each FMDvolume 
//
#ifndef ROOT_TObject
# include <TObject.h>
#endif

class AliFMDRecPoint: public TObject
{
public:
  AliFMDRecPoint();
  AliFMDRecPoint(UShort_t detector,  Char_t   ring, 
	     UShort_t sector,    UShort_t strip, 
	     Float_t  eta,       Float_t  phi,
	     Float_t  edep,      Float_t  particles);
  virtual ~AliFMDRecPoint() {};

  UShort_t     Detector() const        { return fDetector; }
  Char_t       Ring() const            { return fRing; }
  UShort_t     Sector() const          { return fSector; }
  UShort_t     Strip() const           { return fStrip; }
  Float_t      Eta() const             { return fEta; }
  Float_t      Phi() const             { return fPhi; }
  Float_t      Edep() const            { return fEdep; }
  Float_t      Particles() const       { return fParticles; }
  virtual void Print(Option_t* opt="D") const;
  const char*  GetName()                const;
  const char*  GetTitle()               const;
protected:
  UShort_t fDetector;        // Detector #
  Char_t   fRing;            // Ring ID
  UShort_t fSector;          // Sector #
  UShort_t fStrip;           // Strip # 
  Float_t  fEta;             // Eta value 
  Float_t  fPhi;             // Phi value
  Float_t  fEdep;            // Energy deposited 
  Float_t  fParticles;       // Quasi-number of particles 

  ClassDef(AliFMDRecPoint,1)     // Base class for multiplicity data
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
