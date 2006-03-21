#ifndef ALIFMDALTROMAPPING_H
#define ALIFMDALTROMAPPING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIALTROMAPPING_H
# include <AliAltroMapping.h>
#endif
#ifndef ALIFMDUSHORTMAP_H
# include "AliFMDUShortMap.h"
#endif
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif

class AliFMDAltroMapping : public AliAltroMapping
{
public:
  AliFMDAltroMapping();
  Bool_t Hardware2Detector(UInt_t    ddl, UInt_t    hwaddr, 
			   UShort_t& det, Char_t&   ring, 
			   UShort_t& sec, UShort_t& str) const;
  Bool_t Detector2Hardware(UShort_t  det, Char_t    ring, 
			   UShort_t  sec, UShort_t  str,
			   UInt_t&   ddl, UInt_t&   hwaddr) const;
  Int_t  GetHWAdress(Int_t sector, Int_t str, Int_t ring) const
  {
    return GetHWAddress(sector, str, ring);
  }
  Int_t  GetHWAddress(Int_t sector, Int_t str, Int_t ring) const;
  Int_t  GetPadRow(Int_t hwaddr) const;
  Int_t  GetPad(Int_t hwaddr) const;
  Int_t  GetSector(Int_t hwaddr) const;
protected:
  virtual Bool_t ReadMapping();
  virtual void   DeleteMappingArrays();
  
  ClassDef(AliFMDAltroMapping, 1) // Read raw FMD Altro data 
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
