//////////////////////////////////////////////////////////////////////
//
//  Forward Multiplicity Detector have to be reconstructed number of
//  particles in fixed pseudorapidity interval from fNumOfMinRing
//  to fNumOfMaxRing and phi interval from fNumOfMinSector to
//  fNumOfMaxSector
//
//////////////////////////////////////////////////////////////////////
#ifndef ALIFMDPARTICLES_H
# include "AliFMDParticles.h"
#endif
#ifndef __IOSTREAM__
# include <iostream>
#endif
#ifndef __IOMANIP__
# include <iomanip>
#endif

//____________________________________________________________________
ClassImp(AliFMDParticles)

//____________________________________________________________________
AliFMDParticles::AliFMDParticles()
  : fDetector(0),
    fRing('\0'),
    fMinSector(0),
    fMaxSector(0),
    fMinStrip(0),
    fMaxStrip(0),
    fMinEta(0),
    fMaxEta(0),
    fMinPhi(0),
    fMaxPhi(0),
    fParticles(0),
    fMethod(kNaive)
{}

//____________________________________________________________________
AliFMDParticles::AliFMDParticles(UShort_t detector,  Char_t ring, 
				 UShort_t minSector, UShort_t maxSector, 
				 UShort_t minStrip,  UShort_t maxStrip, 
				 Float_t  minEta,    Float_t  maxEta, 
				 Float_t  minPhi,    Float_t  maxPhi,
				 Float_t  particles, UShort_t method)
  : fDetector(detector),
    fRing(ring),
    fMinSector(minSector),
    fMaxSector(maxSector),
    fMinStrip(minStrip),
    fMaxStrip(maxStrip),
    fMinEta(minEta),
    fMaxEta(maxEta),
    fMinPhi(minPhi),
    fMaxPhi(maxPhi),
    fParticles(particles),
    fMethod(method)
{
  switch (fMethod) {
  case kPoission: 
  case kIterative: 
  case kNaive:    
    break;    
  default:        
    Warning("AliFMDParticles", "unknown method: %d", method);
    break;
  }
}


//____________________________________________________________________
void
AliFMDParticles::Print(Option_t* /* opt*/) const
{
  std::cout << "FMD Reconstructed particles: " << fParticles << "\n" 
	    << "  Detector:      FMD" << fDetector << fRing << "\n"
	    << "  Sector range:  [" << fMinSector << "," << fMaxSector << "\n"
	    << "  Strip range:   [" << fMinStrip << "," << fMaxStrip << "\n"
	    << "  Eta range:     [" << fMinEta << "," << fMaxEta << "\n"
	    << "  Phi range:     [" << fMinPhi << "," << fMaxPhi << "\n"
	    << "  Method:        " << std::flush;
  switch (fMethod) {
  case kPoission:  std::cout << "Poission"  << std::endl; break;
  case kIterative: std::cout << "Iterative" << std::endl; break;
  case kNaive:     std::cout << "Naive"     << std::endl; break; 
  default:         std::cout << "Unknown"   << std::endl; break;
  }
}

    
//____________________________________________________________________
//
// EOF
//
