#ifndef ALIFMDMULTSTRIP_H
#define ALIFMDMULTSTRIP_H

// Reconstracted Particles Class: has number of reconstructed
// particles in sectors from NumOfMinSector to NumberOfMaxSector()
// rings from NumOfMinRing to NumOfMaxRing for each FMDvolume
//
#ifndef ALIFMDMULT_H
# include "AliFMDMult.h"
#endif

class AliFMDMultStrip: public AliFMDMult
{
public:
  AliFMDMultStrip();
  AliFMDMultStrip (UShort_t detector,  Char_t   ring, 
		   UShort_t sector,    UShort_t strip, 
		   Float_t  eta,       Float_t  phi,
		   Float_t  edep,      Float_t  particles, 
		   UShort_t method);
  virtual ~AliFMDMultStrip(){};

  UShort_t     Detector() const        { return fDetector; }
  Char_t       Ring() const            { return fRing; }
  UShort_t     Sector() const          { return fSector; }
  UShort_t     Strip() const           { return fStrip; }
  Float_t      Eta() const             { return fEta; }
  Float_t      Phi() const             { return fPhi; }
  Float_t      Edep() const            { return fEdep; }
  virtual void Print(Option_t* opt="D") const;
protected:
  UShort_t fDetector;        // Detector #
  Char_t   fRing;            // Ring ID
  UShort_t fSector;          // Sector #
  UShort_t fStrip;           // Strip # 
  Float_t  fEta;             // Eta value 
  Float_t  fPhi;             // Phi value
  Float_t  fEdep;            // Energy deposited 

  ClassDef(AliFMDMultStrip,1) // Rec. Multiplicity in a strip
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
