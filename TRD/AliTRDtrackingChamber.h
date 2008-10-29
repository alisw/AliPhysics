#ifndef ALITRDTRACKINGCHAMBER_H
#define ALITRDTRACKINGCHAMBER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingChamber.h 22646 2007-11-29 18:13:40Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// Data container for one TRD chamber                                     // 
//                                                                        // 
// Authors:                                                               //
//                                                                        //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//                                                                        // 
////////////////////////////////////////////////////////////////////////////

#ifndef ALITRDCHAMBERTIMEBIN_H
#include "AliTRDchamberTimeBin.h"
#endif


#ifndef ALITRDSEED_H
#include "AliTRDseed.h"
#endif

class AliTRDCalDet;
class AliTRDgeometry;
class AliTRDReconstructor;
class AliTRDtrackingChamber : public TObject
{

public:
  enum{
    kNTimeBins = AliTRDseed::knTimebins
  };
  AliTRDtrackingChamber();
  virtual ~AliTRDtrackingChamber(){}
  
  void     Bootstrap(const AliTRDReconstructor *rec);
  Bool_t   Build(AliTRDgeometry *geo, const AliTRDCalDet *cal, Bool_t hlt = kFALSE);
  void     Clear(const Option_t *opt = 0x0);
  Int_t    GetDetector() const {return fDetector;}
  Int_t    GetNClusters() const;
  Double_t GetQuality();
  Bool_t   GetSeedingLayer(AliTRDchamberTimeBin *&layer, AliTRDgeometry *geo, const AliTRDReconstructor *rec);
  Float_t  GetX()        const {return fX0;}
  AliTRDchamberTimeBin* GetTB(int tb) {return tb >= 0 && tb < kNTimeBins ? &fTB[tb] : 0x0;}
  void     InsertCluster(AliTRDcluster *c, Int_t index);
  
  void     Print(Option_t *opt = 0x0) const;

  void     SetDetector(Int_t det) { fDetector = det;}
  void     SetOwner();

private:
  Int_t         fDetector;  // detector number
  Float_t       fX0;        // approximate position of the pad plane
  
  AliTRDchamberTimeBin fTB[kNTimeBins];    // time bins 
  
  
  ClassDef(AliTRDtrackingChamber, 1)  // TRD tracker container for one chamber
};

#endif  // ALITRDTRACKINGCHAMBER_H
