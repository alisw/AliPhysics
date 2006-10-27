#ifndef ALIESDFMD_H
#define ALIESDFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//___________________________________________________________________
//
// AliESDFMD is the Event Summary Data entry for the FMD.  It contains
// a rough estimate of the charged particle multiplicity in each strip
// of the FMD.    It also contains the psuedo-rapidity of each strip.
// This is important, as it varies from event to event, due to a
// finite interaction point probability distribution. 
//
#ifndef ROOT_TObject
# include <TObject.h>
#endif
#ifndef ALIFMDFLOATMAP_H
# include "AliFMDFloatMap.h"
#endif

//___________________________________________________________________
class AliESDFMD : public TObject
{
public:
  AliESDFMD();
  AliESDFMD(const AliESDFMD& other);
  AliESDFMD& operator=(const AliESDFMD& other);
  virtual ~AliESDFMD() {}

  void Clear(Option_t *option="");
  Float_t Multiplicity(UShort_t detector, Char_t ring, 
		       UShort_t sector, UShort_t strip) const;
  Float_t Eta(UShort_t detector, Char_t ring, 
	      UShort_t sector, UShort_t strip) const;
  void SetMultiplicity(UShort_t detector, Char_t ring, 
		       UShort_t sector, UShort_t strip, 
		       Float_t mult);
  void SetEta(UShort_t detector, Char_t ring, 
	      UShort_t sector, UShort_t strip, 
	      Float_t mult);

  UShort_t MaxDetectors() const { return fMultiplicity.MaxDetectors(); }
  UShort_t MaxRings()     const { return fMultiplicity.MaxRings(); }
  UShort_t MaxSectors()   const { return fMultiplicity.MaxSectors(); }
  UShort_t MaxStrips()    const { return fMultiplicity.MaxStrips(); }
  void Print(Option_t* option="") const;
  enum {
    kInvalidMult = 1000
  };
  enum {
    kInvalidEta = 1000
  };
protected:
  AliFMDFloatMap fMultiplicity; // Psuedo multplicity per strip
  AliFMDFloatMap fEta;          // Psuedo-rapidity per strip
  
  ClassDef(AliESDFMD,2)  // ESD info from FMD
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
