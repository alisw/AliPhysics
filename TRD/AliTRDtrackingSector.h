#ifndef ALITRDTRACKINGSECTOR_H
#define ALITRDTRACKINGSECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingSector.h 22646 2007-11-29 18:13:40Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Data container for one TRD supermodule                                 // 
//                                                                        // 
// Authors:                                                               //
//                                                                        //
//    Marian Ivanov <M.Ivanov@gsi.de>                                     //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDGEOMETRY_H
#include "AliTRDgeometry.h"
#endif

class AliTRDCalDet;
class AliTRDReconstructor;
class AliTRDtrackingChamber;
class AliTRDtrackingSector 
{

public:
  AliTRDtrackingSector();
  AliTRDtrackingSector(AliTRDgeometry* geo, Int_t gs);
  virtual ~AliTRDtrackingSector(){;}
    
  void     Clear(const Option_t *opt = 0x0);
  Int_t    GetNChambers() const             { return fN; }
  Double_t GetX(Int_t pl) const                  { return pl >=0 && pl < AliTRDgeometry::kNlayer ? fX0[pl] : 0.; }
  AliTRDtrackingChamber* GetChamber(Int_t i) const  { return i>=0 && i < fN ? fChamber[i] : 0x0;  }
  AliTRDtrackingChamber* GetChamber(Int_t stack, Int_t plane, Bool_t build = kFALSE);
  AliTRDtrackingChamber** GetStack(Int_t stack);
  Int_t    GetSector() const {return fSector;}	

  void     Init(const AliTRDReconstructor *rec, const AliTRDCalDet *cal);
  void     Print(Option_t *opt = 0x0) const;
  
  void     SetGeometry(AliTRDgeometry *geo) {fGeom = geo;}
  
private:
  AliTRDtrackingSector(const AliTRDtrackingSector &/*t*/);
  AliTRDtrackingSector &operator=(const AliTRDtrackingSector &/*t*/);


private:
  Char_t         fSector;           // Sector# in AliTRDgeometry
  UChar_t        fN;                // Total number of chambers allocated
  Char_t         fIndex[AliTRDgeometry::kNdets]; // indexes of allocated chambers
  Float_t        fX0[AliTRDgeometry::kNlayer];      // average position of pad plane for each plane
  AliTRDgeometry *fGeom;            // Geometry
  AliTRDtrackingChamber *fChamber[AliTRDgeometry::kNdets];// chambers   
  AliTRDtrackingChamber *fStack[AliTRDgeometry::kNlayer]; //! temporary holding one stack

  ClassDef(AliTRDtrackingSector, 1) // TRD tracker container for one sector
};


#endif // ALITRDTRACKINGSECTOR_H

