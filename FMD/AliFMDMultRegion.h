#ifndef ALIFMDMULTREGION_H
#define ALIFMDMULTREGION_H

/* Reconstracted Particles Class: has number of reconstructed
 * particles in sectors from NumOfMinSector to NumberOfMaxSector()
 * rings from NumOfMinRing to NumOfMaxRing for each FMDvolume
 */
#ifndef ALIFMDMULT_H
# include "AliFMDMult.h"
#endif

class AliFMDMultRegion: public AliFMDMult
{
public:
  AliFMDMultRegion();
  AliFMDMultRegion (UShort_t detector,  Char_t ring, 
		    UShort_t minSector, UShort_t maxSector, 
		    UShort_t minStrip,  UShort_t maxStrip, 
		    Float_t  minEta,    Float_t  maxEta, 
		    Float_t  minPhi,    Float_t  maxPhi,
		    Float_t  particles, UShort_t method);
  virtual ~AliFMDMultRegion(){};

  UShort_t     Detector() const        { return fDetector; }
  Char_t       Ring() const            { return fRing; }
  UShort_t     MinSector() const       { return fMinSector; }
  UShort_t     MaxSector() const       { return fMaxSector; }
  UShort_t     MinStrip() const        { return fMinStrip; }
  UShort_t     MaxStrip() const        { return fMaxStrip; }
  Float_t      MinEta() const          { return fMinEta; }
  Float_t      MaxEta() const          { return fMaxEta; }
  Float_t      MinPhi() const          { return fMinPhi; }
  Float_t      MaxPhi() const          { return fMaxPhi; }
  virtual void Print(Option_t* opt="") const;
protected:
  UShort_t fDetector;        // Detector #
  Char_t   fRing;            // Ring ID
  UShort_t fMinSector;       // First sector of this region
  UShort_t fMaxSector;       // Last sector of this region
  UShort_t fMinStrip;        // First strip of this region
  UShort_t fMaxStrip;        // Second strip of this region  
  Float_t  fMinEta;          // Least eta covered
  Float_t  fMaxEta;          // Largest eta covered
  Float_t  fMinPhi;          // Least phi covered
  Float_t  fMaxPhi;          // Largest phi covered

  ClassDef(AliFMDMultRegion,1) // Rec. Multiplicity in a eta,phi region
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
