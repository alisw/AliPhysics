#ifndef ALIITSCORRECTSDDPOINTS_H
#define ALIITSCORRECTSDDPOINTS_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to apply SDD map corrections                            //
// for voltage divider shape and doping fluctuations             //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSsegmentationSDD.h"
#include "TObjArray.h"

class AliITSCorrectSDDPoints : public TObject {
 public:
  AliITSCorrectSDDPoints();
  AliITSCorrectSDDPoints(TString filname);
  ~AliITSCorrectSDDPoints();
  Float_t GetCorrection(Int_t modId, Float_t zloc, Float_t xloc) const;
  Float_t GetCorrectedXloc(Int_t modId, Float_t zloc, Float_t xloc) const{
    Float_t dx=GetCorrection(modId,zloc,xloc);
    return xloc-dx;
  }
 private:
  AliITSCorrectSDDPoints(const AliITSCorrectSDDPoints& csdd);
  AliITSCorrectSDDPoints& operator=(const AliITSCorrectSDDPoints& csdd);
 protected:
  TObjArray* fArrayOfMaps;                 // 520 AliITSCorrMapSDD objects
  AliITSsegmentationSDD* fSegmentationSDD; // SDD segmentation
  ClassDef(AliITSCorrectSDDPoints,0);
};
#endif
