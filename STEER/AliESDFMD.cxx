/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//____________________________________________________________________
//
//  ESD information from the FMD 
//  Contains information on:
//	Charged particle multiplicty per strip (rough estimate)
//	Psuedo-rapdity per strip
//  Latest changes by Christian Holm Christensen
//
#include "AliESDFMD.h"		// ALIFMDESD_H
#include "AliLog.h"		// ALILOG_H
#include "Riostream.h"		// ROOT_Riostream

//____________________________________________________________________
ClassImp(AliESDFMD)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif


//____________________________________________________________________
AliESDFMD::AliESDFMD()
  : fMultiplicity(AliFMDFloatMap::kMaxDetectors, 
		  AliFMDFloatMap::kMaxRings, 
		  AliFMDFloatMap::kMaxSectors, 
		  AliFMDFloatMap::kMaxStrips), 
    fEta(AliFMDFloatMap::kMaxDetectors, 
	 AliFMDFloatMap::kMaxRings, 
	 1,
	 AliFMDFloatMap::kMaxStrips)
{
  // Default CTOR
}
  
//____________________________________________________________________
AliESDFMD::AliESDFMD(const AliESDFMD& other)
  : TObject(other), 
    fMultiplicity(other.fMultiplicity),
    fEta(other.fEta)
{
  // Default CTOR
}

//____________________________________________________________________
AliESDFMD& 
AliESDFMD::operator=(const AliESDFMD& other)
{
  // Default CTOR
  fMultiplicity = other.fMultiplicity;
  fEta          = other.fEta;
  return *this;
}

//____________________________________________________________________
void
AliESDFMD::Clear(Option_t* )
{
  fMultiplicity.Reset(kInvalidMult);
  fEta.Reset(kInvalidEta);
}


//____________________________________________________________________
Float_t
AliESDFMD::Multiplicity(UShort_t detector, Char_t ring, UShort_t sector, 
			UShort_t strip) const
{
  // Return rough estimate of charged particle multiplicity in the
  // strip FMD<detector><ring>[<sector>,<strip>]. 
  // 
  // Note, that this should at most be interpreted as the sum
  // multiplicity of secondaries and primaries. 
  return fMultiplicity(detector, ring, sector, strip);
}

//____________________________________________________________________
Float_t
AliESDFMD::Eta(UShort_t detector, Char_t ring, UShort_t /* sector */, 
	       UShort_t strip) const
{
  // Return pseudo-rapidity of the strip
  // FMD<detector><ring>[<sector>,<strip>].  (actually, the sector
  // argument is ignored, as it is assumed that the primary vertex is
  // a (x,y) = (0,0), and that the modules are aligned with a
  // precision better than 2 degrees in the azimuthal angle). 
  // 
  return fEta(detector, ring, 0, strip);
}

//____________________________________________________________________
void
AliESDFMD::SetMultiplicity(UShort_t detector, Char_t ring, UShort_t sector, 
			   UShort_t strip, Float_t mult)
{
  // Return rough estimate of charged particle multiplicity in the
  // strip FMD<detector><ring>[<sector>,<strip>]. 
  // 
  // Note, that this should at most be interpreted as the sum
  // multiplicity of secondaries and primaries. 
  fMultiplicity(detector, ring, sector, strip) = mult;
}

//____________________________________________________________________
void
AliESDFMD::SetEta(UShort_t detector, Char_t ring, UShort_t /* sector */, 
		  UShort_t strip, Float_t eta)
{
  // Set pseudo-rapidity of the strip
  // FMD<detector><ring>[<sector>,<strip>].  (actually, the sector
  // argument is ignored, as it is assumed that the primary vertex is
  // a (x,y) = (0,0), and that the modules are aligned with a
  // precision better than 2 degrees in the azimuthal angle). 
  // 
  fEta(detector, ring, 0, strip) = eta;
}

//____________________________________________________________________
void
AliESDFMD::Print(Option_t* /* option*/) const
{
  // Print all information to standard output. 
  std::cout << "AliESDFMD:" << std::endl;
  for (UShort_t det = 1; det <= fMultiplicity.MaxDetectors(); det++) {
    for (UShort_t ir = 0; ir < fMultiplicity.MaxRings(); ir++) {
      Char_t ring = (ir == 0 ? 'I' : 'O');
      std::cout << "FMD" << det << ring << ":" << std::endl;
      for  (UShort_t sec = 0; sec < fMultiplicity.MaxSectors(); sec++) {
	std::cout << " Sector # " << sec << ":" << std::flush;
	for (UShort_t str = 0; str < fMultiplicity.MaxStrips(); str++) {
	  if (str % 6 == 0) std::cout << "\n  " << std::flush;
	  Float_t m = fMultiplicity(det, ring, sec, str);
	  Float_t e = fEta(det, ring, 0, str);
	  if (m == kInvalidMult && e == kInvalidEta) break;
	  if (m == kInvalidMult) std::cout << " ---- ";
	  else                   std::cout << Form("%6.3f", m);
	  std::cout << "/";
	  if (e == kInvalidEta)  std::cout << " ---- ";
	  else                   std::cout << Form("%-6.3f", e);
	  std::cout << " " << std::flush;
	}
	std::cout << std::endl;
      }
    }
  }
}


//____________________________________________________________________
//
// EOF
//
