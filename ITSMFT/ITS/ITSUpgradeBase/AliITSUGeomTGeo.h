#ifndef ALIITSUGEOMTGEO_H
#define ALIITSUGEOMTGEO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
//  AliITSUGeomTGeo is a simple interface class to TGeoManager         //
//  It is used in the simulation and reconstruction in order to        //
//  query the TGeo ITS geometry                                        //
//                                                                     //
//  author - cvetan.cheshkov@cern.ch                                   //
//  15/02/2007                                                         //
//  adapted to ITSupg 18/07/2012 - ruben.shahoyan@cern.ch              //
//  RS: in order to preserve the static character of the class but     //
//  make it dynamically access geometry, we need to check in every     //
//  method if the structures are initialized. To be converted to       //
//  singleton at later stage.                                          //
//                                                                     //
//  Note on the upgrade chip types:                                    //
//  The coarse type defines chips served by different classes,         //
//  like Pix. Each such a chip type can have kMaxSegmPerChipType       //
//  segmentations (pitch etc.) whose parameteres are stored in the     //
//  AliITSsegmentation derived class (like AliITSMFTSegmentationPix)     //
//  This allows to have in the setup chips served by the same          //
//  classes but with different segmentations.                          //
//  The full chip type is composed as:                                 //
//  CoarseType*kMaxSegmPerChipType + segmentationType                  //
//  The only requirement on the segmentationType that should be        //
//  < kMaxSegmPerChipType.                                             //
//  The methods like GetLayerChipTypeID return the full chip type      //
//                                                                     //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TGeoMatrix.h>
#include <TString.h>
#include <TObjArray.h>
#include "AliITSMFTAux.h"
#include "AliITSMFTSimuParam.h"
#include "AliITSMFTGeomTGeo.h"


class TGeoPNEntry;
class TDatime;
class AliITSMFTSegmentationPix;



class AliITSUGeomTGeo : public AliITSMFTGeomTGeo {
public:
  enum {kITSVNA, kITSVUpg}; // ITS version
  AliITSUGeomTGeo(Bool_t build = kFALSE, Bool_t loadSegmentations = kTRUE);
  AliITSUGeomTGeo(const AliITSUGeomTGeo &src);
  //virtual ~AliITSUGeomTGeo(); 
  
  Bool_t
  GetChipId(Int_t idx,Int_t &lay,Int_t &sta,Int_t &ssta,Int_t &mod,Int_t &chip)
  const;

  void Print(Option_t *opt="")  const;

private:
  void    BuildITS(Bool_t loadSegm);
  virtual TGeoHMatrix* ExtractMatrixSens(Int_t index) const;

  Int_t  fVersion;             // ITS Version 

  ClassDef(AliITSUGeomTGeo, 3) // ITS geometry based on TGeo
};

#endif
