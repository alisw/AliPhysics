#ifndef ALIITSMFTGEOMTGEO_H
#define ALIITSMFTGEOMTGEO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
//  AliITSMFTGeomTGeo is a simple interface class to TGeoManager       //
//  It is used in the simulation and reconstruction in order to        //
//  query the TGeo geometry                                            //
//                                                                     //
//  This is the base class for the AliITSUGeomTGeo and                 //
//  AliMFTGeomTGeo geometry manager classes                            //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>
#include "AliITSMFTAux.h"
#include <TGeoMatrix.h>
#include <TObjArray.h>

class TGeoPNEntry;
class TDatime;
class AliITSMFTSegmentationPix;



class AliITSMFTGeomTGeo : public TObject {

public:
  AliITSMFTGeomTGeo();
  AliITSMFTGeomTGeo(const AliITSMFTGeomTGeo &src);
  virtual ~AliITSMFTGeomTGeo(); 

  //
  Int_t  GetNChips() const {return fNChips;}  

  //
  // Attention: these are the matrices for the alignable volumes of the chips,
  // i.e. not necessarily the sensors
  TGeoHMatrix* GetMatrix(Int_t index)                                     const;
  Bool_t GetTranslation(Int_t index, Double_t t[3])                       const;
  Bool_t GetRotation(Int_t index, Double_t r[9])                          const;

  virtual Bool_t GetOrigMatrix(Int_t index, TGeoHMatrix &m) const = 0;
  
  Bool_t GetOrigTranslation(Int_t index, Double_t t[3])                   const;
  Bool_t GetOrigRotation(Int_t index, Double_t r[9])                      const;

  //
  const TGeoHMatrix* GetMatrixT2L(Int_t index);
  const TGeoHMatrix* GetMatrixSens(Int_t index);
  //
  Bool_t GetTrackingMatrix(Int_t index, TGeoHMatrix &m);

  //
  // Attention: these are transformations wrt sensitive volume!
  void   LocalToGlobal(Int_t index, const Double_t *loc, Double_t *glob);
  void   GlobalToLocal(Int_t index, const Double_t *glob, Double_t *loc);
  void   LocalToGlobalVect(Int_t index, const Double_t *loc, Double_t *glob);
  void   GlobalToLocalVect(Int_t index, const Double_t *glob, Double_t *loc);

  //
  const AliITSMFTSegmentationPix* GetSegmentationByID(Int_t id) const;
  TObjArray *GetSegmentations() const {return fSegm;}
  //

protected:
  AliITSMFTGeomTGeo& operator=(const AliITSMFTGeomTGeo &geom);

  virtual TGeoHMatrix* ExtractMatrixSens(Int_t index) const = 0;
  virtual TGeoPNEntry* GetPNEntry(Int_t index)        const = 0;

  void         FetchMatrices();
  void         CreateT2LMatrices();
  TGeoHMatrix* ExtractMatrixT2L(Int_t index)                      const;
  //
  //
protected:
  //
  Int_t  fNChips;              // The total number of chips
  TObjArray* fMatSens;         // Sensor's matrices pointers in the geometry
  TObjArray* fMatT2L;          // Tracking to Local matrices pointers in the geometry
  TObjArray* fSegm;            // segmentations
  //
  ClassDef(AliITSMFTGeomTGeo, 2) // ITS geometry based on TGeo
};

//_____________________________________________________________________________________________
inline const TGeoHMatrix* AliITSMFTGeomTGeo::GetMatrixSens(Int_t index)
{
  // access global to sensor matrix
  if (!fMatSens) FetchMatrices();
  return (TGeoHMatrix*)fMatSens->At(index);
}

//_____________________________________________________________________________________________
inline const TGeoHMatrix* AliITSMFTGeomTGeo::GetMatrixT2L(Int_t index)
{
  // access tracking to local matrix
  if (!fMatT2L) FetchMatrices();
  return (TGeoHMatrix*)fMatT2L->At(index);
}

//______________________________________________________________________
inline void AliITSMFTGeomTGeo::LocalToGlobal(Int_t index,const Double_t *loc, Double_t *glob)
{
  // sensor local to global 
  GetMatrixSens(index)->LocalToMaster(loc,glob);
}

//______________________________________________________________________
inline void AliITSMFTGeomTGeo::GlobalToLocal(Int_t index, const Double_t *glob, Double_t *loc)
{
  // global to sensor local 
  GetMatrixSens(index)->MasterToLocal(glob,loc);
}

//______________________________________________________________________
inline void AliITSMFTGeomTGeo::LocalToGlobalVect(Int_t index, const Double_t *loc, Double_t *glob)
{
  // sensor local to global 
  GetMatrixSens(index)->LocalToMasterVect(loc,glob);
}

//______________________________________________________________________
inline void AliITSMFTGeomTGeo::GlobalToLocalVect(Int_t index, const Double_t *glob, Double_t *loc)
{
  // global to sensor local
  GetMatrixSens(index)->MasterToLocalVect(glob,loc);
}

//_____________________________________________________________________________
inline const AliITSMFTSegmentationPix* AliITSMFTGeomTGeo::GetSegmentationByID(Int_t id) const 
{
  // get segmentation by ID
  return fSegm ? (AliITSMFTSegmentationPix*)fSegm->At(id) : 0;
}

#endif
