#ifndef ALIITSUGEOMTGEO_H
#define ALIITSUGEOMTGEO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
//  AliITSUGeomTGeo is a simple interface class to TGeoManager       //
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
//  Note on the upgrade detector types:                                //
//  The coarse type defines detectors served by different classes,     //
//  like Pix. Each such a detector type can have kMaxSegmPerDetType    //
//  segmentations (pitch etc.) whose parameteres are stored in the     //
//  AliITSsegmentation derived class (like AliITSUSegmentationPix)     //
//  This allows to have in the setup modules served by the same        //
//  classes but with different segmentations.                          //
//  The full detector type is composed as:                             //
//  CoarseType*kMaxSegmPerDetType + segmentationType                   //
//  The only requirement on the segmentationType that should be        //
//  < kMaxSegmPerDetType.                                              //
//  The methods like GetLayerDetTypeID return the full detector type   //
//                                                                     //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TGeoMatrix.h>
#include <TString.h>
#include <TObjArray.h>

class TGeoPNEntry;
class TDatime;
class AliITSsegmentation;

class AliITSUGeomTGeo : public TObject {

 public:
  enum {kITSVNA, kITSVUpg}; // ITS version
  enum {kDetTypePix=0, kNDetTypes, kMaxSegmPerDetType=10}; // defined detector types (each one can have different segmentations)
  //
  AliITSUGeomTGeo(Bool_t build = kFALSE, Bool_t loadSegmentations = kTRUE);
  virtual ~AliITSUGeomTGeo(); 
  AliITSUGeomTGeo(const AliITSUGeomTGeo &src);
  AliITSUGeomTGeo& operator=(const AliITSUGeomTGeo &geom);
  //
  Int_t  GetNModules()                                                    const {return fNModules;}
  Int_t  GetNDetectors(Int_t lay)                                         const {return fNDetectors[lay];}
  Int_t  GetNLadders(Int_t lay)                                           const {return fNLadders[lay];}
  Int_t  GetNLayers()                                                     const {return fNLayers;}
  
  Int_t  GetModuleIndex(Int_t lay,int detInLay)                           const {return GetFirstModIndex(lay)+detInLay;}
  Int_t  GetModuleIndex(Int_t lay,Int_t lad,Int_t det)                    const;
  Bool_t GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det)        const;
  Int_t  GetLayer(Int_t index)                                            const;
  Int_t  GetLadder(Int_t index)                                           const;
  Int_t  GetModIdInLayer(Int_t index)                                     const;
  Int_t  GetModIdInLadder(Int_t index)                                    const;
  //
  Int_t  GetLastModIndex(Int_t lay)                                       const {return fLastModIndex[lay];}
  Int_t  GetFirstModIndex(Int_t lay)                                      const {return (lay==0) ? 0:fLastModIndex[lay-1]+1;}
  //  
  const char *GetSymName(Int_t index)                                     const;
  const char *GetSymName(Int_t lay,Int_t lad,Int_t det)                   const;
  //
  // Attention: these are the matrices for the alignable volumes of the modules, i.e. not necessarily the sensors
  TGeoHMatrix* GetMatrix(Int_t index)                                     const;
  TGeoHMatrix* GetMatrix(Int_t lay,Int_t lad,Int_t det)                   const;
  Bool_t GetTranslation(Int_t index, Double_t t[3])                       const;
  Bool_t GetTranslation(Int_t lay,Int_t lad,Int_t det, Double_t t[3])     const;
  Bool_t GetRotation(Int_t index, Double_t r[9])                          const;
  Bool_t GetRotation(Int_t lay,Int_t lad,Int_t det, Double_t r[9])        const;
  Bool_t GetOrigMatrix(Int_t index, TGeoHMatrix &m)                       const;
  Bool_t GetOrigMatrix(Int_t lay,Int_t lad,Int_t det, TGeoHMatrix &m)     const;
  Bool_t GetOrigTranslation(Int_t index, Double_t t[3])                   const;
  Bool_t GetOrigTranslation(Int_t lay,Int_t lad,Int_t det, Double_t t[3]) const;
  Bool_t GetOrigRotation(Int_t index, Double_t r[9])                      const;
  Bool_t GetOrigRotation(Int_t lay,Int_t lad,Int_t det, Double_t r[9])    const;
  //
  const TGeoHMatrix* GetMatrixT2L(Int_t index);
  const TGeoHMatrix* GetMatrixT2L(Int_t lay,Int_t lad,Int_t det)  {return GetMatrixT2L( GetModuleIndex(lay,lad,det) );}
  const TGeoHMatrix* GetMatrixSens(Int_t index);
  const TGeoHMatrix* GetMatrixSens(Int_t lay,Int_t lad,Int_t det) {return GetMatrixSens( GetModuleIndex(lay,lad,det) );}
  //
  Bool_t GetTrackingMatrix(Int_t index, TGeoHMatrix &m);
  Bool_t GetTrackingMatrix(Int_t lay,Int_t lad,Int_t det, TGeoHMatrix &m);
  //
  // Attention: these are transformations wrt sensitive volume!
  void   LocalToGlobal(Int_t index, const Double_t *loc, Double_t *glob);
  void   LocalToGlobal(Int_t lay, Int_t lad, Int_t det,const Double_t *loc, Double_t *glob);
  //
  void   GlobalToLocal(Int_t index, const Double_t *glob, Double_t *loc);
  void   GlobalToLocal(Int_t lay, Int_t lad, Int_t det,const Double_t *glob, Double_t *loc);
  //
  void   LocalToGlobalVect(Int_t index, const Double_t *loc, Double_t *glob);
  void   GlobalToLocalVect(Int_t index, const Double_t *glob, Double_t *loc);
  Int_t  GetLayerDetTypeID(Int_t lr)                                         const;
  Int_t  GetModuleDetTypeID(Int_t id)                                        const;
  //
  const AliITSsegmentation* GetSegmentationByID(Int_t id)                    const;
  const AliITSsegmentation* GetSegmentation(Int_t lr)                        const;
  TObjArray*          GetSegmentations()                                     const {return (TObjArray*)fSegm;}
  virtual void Print(Option_t *opt="")  const;

  static const char* GetITSVolPattern()                                 {return fgkITSVolName;}
  static const char* GetITSLayerPattern()                               {return fgkITSLrName;}
  static const char* GetITSLadderPattern()                              {return fgkITSLadName;}
  static const char* GetITSModulePattern()                              {return fgkITSModName;}
  static const char* GetITSSensorPattern()                              {return fgkITSSensName;}
  static const char* GetITSsegmentationFileName()                       {return fgITSsegmFileName.Data();}
  static const char* GetDetTypeName(Int_t i);
  static void        SetITSsegmentationFileName(const char* nm)         {fgITSsegmFileName = nm;}
  static UInt_t      ComposeDetTypeID(UInt_t segmId);
  //
  static const char *ComposeSymNameITS();
  static const char *ComposeSymNameLayer(Int_t lr);
  static const char *ComposeSymNameLadder(Int_t lr, Int_t lad);
  static const char *ComposeSymNameModule(Int_t lr, Int_t lad, int det);
  //
  // hack to avoid using AliGeomManager
  Int_t              LayerToVolUID(Int_t lay,int detInLay)        const {return GetModuleIndex(lay,detInLay);}
  static Int_t       ModuleVolUID(Int_t mod)                            {return mod;}
  //
 protected:
  void         FetchMatrices();
  void         CreateT2LMatrices();
  TGeoHMatrix* ExtractMatrixT2L(Int_t index)                      const;
  TGeoHMatrix* ExtractMatrixSens(Int_t index)                     const;
  Bool_t       GetLayer(Int_t index,Int_t &lay,Int_t &index2)     const;
  TGeoPNEntry* GetPNEntry(Int_t index)                            const;
  Int_t        ExtractNumberOfDetectors(Int_t lay)                const;
  Int_t        ExtractNumberOfLadders(Int_t lay)                  const;
  Int_t        ExtractLayerDetType(Int_t lay)                     const;
  Int_t        ExtractNumberOfLayers()                            const;
  void         BuildITS(Bool_t loadSegm);
  //
 protected:
  //
  //
  Int_t  fVersion;             // ITS Version 
  Int_t  fNLayers;             // number of layers
  Int_t  fNModules;            // The total number of modules
  Int_t *fNLadders;            //[fNLayers] Array of the number of ladders/layer(layer)
  Int_t *fLrDetType;           //[fNLayers] Array of layer detector types
  Int_t *fNDetectors;          //[fNLayers] Array of the number of detector/ladder(layer)
  Int_t *fLastModIndex;        //[fNLayers] max ID of the detctor in the layer
  //
  TObjArray* fMatSens;         // Sensor's matrices pointers in the geometry
  TObjArray* fMatT2L;          // Tracking to Local matrices pointers in the geometry
  TObjArray* fSegm;            // segmentations
  //
  static const char*  fgkITSVolName;             // ITS mother volume name
  static const char*  fgkITSLrName;              // ITS Layer name
  static const char*  fgkITSLadName;             // ITS Ladder name 
  static const char*  fgkITSModName;             // ITS Module name 
  static const char*  fgkITSSensName;            // ITS Sensor name 
  static const char*  fgkITSDetTypeName[kNDetTypes]; // ITS upg detType Names
  //
  static TString      fgITSsegmFileName;         // file name for segmentations
  //
  ClassDef(AliITSUGeomTGeo, 1) // ITS geometry based on TGeo
};

//_____________________________________________________________________________________________
inline const char *AliITSUGeomTGeo::GetSymName(Int_t lay,Int_t lad,Int_t det) const    
{
  // sym name
  return GetSymName(GetModuleIndex(lay,lad,det));
}

//_____________________________________________________________________________________________
inline TGeoHMatrix* AliITSUGeomTGeo::GetMatrix(Int_t lay,Int_t lad,Int_t det) const    
{
  // module current matrix
  return GetMatrix(GetModuleIndex(lay,lad,det));
}

//_____________________________________________________________________________________________
inline Bool_t AliITSUGeomTGeo::GetTranslation(Int_t lay,Int_t lad,Int_t det, Double_t t[3]) const    
{
  // translation
  return GetTranslation(GetModuleIndex(lay,lad,det),t); 
}

//_____________________________________________________________________________________________
inline Bool_t AliITSUGeomTGeo::GetRotation(Int_t lay,Int_t lad,Int_t det, Double_t r[9]) const    
{
  // rot
  return GetRotation(GetModuleIndex(lay,lad,det),r); 
}

//_____________________________________________________________________________________________
inline Bool_t AliITSUGeomTGeo::GetOrigMatrix(Int_t lay,Int_t lad,Int_t det, TGeoHMatrix &m) const    
{
  // orig matrix
  return GetOrigMatrix(GetModuleIndex(lay,lad,det),m); 
}

//_____________________________________________________________________________________________
inline Bool_t AliITSUGeomTGeo::GetOrigTranslation(Int_t lay,Int_t lad,Int_t det, Double_t t[3]) const    
{
  // orig trans
  return GetOrigTranslation(GetModuleIndex(lay,lad,det),t); 
}

//_____________________________________________________________________________________________
inline Bool_t AliITSUGeomTGeo::GetOrigRotation(Int_t lay,Int_t lad,Int_t det, Double_t r[9]) const    
{
  // orig rot
  return GetOrigRotation(GetModuleIndex(lay,lad,det),r); 
}

//_____________________________________________________________________________________________
inline Bool_t AliITSUGeomTGeo::GetTrackingMatrix(Int_t lay,Int_t lad,Int_t det, TGeoHMatrix &m)
{
  // tracking mat
  return GetTrackingMatrix(GetModuleIndex(lay,lad,det),m); 
}

//_____________________________________________________________________________________________
inline Int_t  AliITSUGeomTGeo::GetLayerDetTypeID(Int_t lr) const  
{
  // detector type ID of layer
  return fLrDetType[lr];
}

//_____________________________________________________________________________________________
inline Int_t  AliITSUGeomTGeo::GetModuleDetTypeID(Int_t id) const  
{
  // detector type ID of module
  return GetLayerDetTypeID(GetLayer(id));
} 

//_____________________________________________________________________________________________
inline const TGeoHMatrix* AliITSUGeomTGeo::GetMatrixSens(Int_t index)
{
  // access global to sensor matrix
  if (!fMatSens) FetchMatrices();
  return (TGeoHMatrix*)fMatSens->At(index);
}

//_____________________________________________________________________________________________
inline const TGeoHMatrix* AliITSUGeomTGeo::GetMatrixT2L(Int_t index)
{
  // access tracking to local matrix
  if (!fMatT2L) FetchMatrices();
  return (TGeoHMatrix*)fMatT2L->At(index);
}

//______________________________________________________________________
inline void AliITSUGeomTGeo::LocalToGlobal(Int_t index,const Double_t *loc, Double_t *glob)
{
  // sensor local to global 
  GetMatrixSens(index)->LocalToMaster(loc,glob);
}

//______________________________________________________________________
inline void AliITSUGeomTGeo::GlobalToLocal(Int_t index, const Double_t *glob, Double_t *loc)
{
  // global to sensor local 
  GetMatrixSens(index)->MasterToLocal(glob,loc);
}

//______________________________________________________________________
inline void AliITSUGeomTGeo::LocalToGlobalVect(Int_t index, const Double_t *loc, Double_t *glob)
{
  // sensor local to global 
  GetMatrixSens(index)->LocalToMasterVect(loc,glob);
}

//______________________________________________________________________
inline void AliITSUGeomTGeo::GlobalToLocalVect(Int_t index, const Double_t *glob, Double_t *loc)
{
  // global to sensor local
  GetMatrixSens(index)->MasterToLocalVect(glob,loc);
}

//_____________________________________________________________________________________________
inline void AliITSUGeomTGeo::LocalToGlobal(Int_t lay, Int_t lad, Int_t det,const Double_t *loc, Double_t *glob)
{
  // Local2Master (sensor)
  LocalToGlobal(GetModuleIndex(lay,lad,det), loc, glob);
}

//_____________________________________________________________________________________________
inline void AliITSUGeomTGeo::GlobalToLocal(Int_t lay, Int_t lad, Int_t det,const Double_t *glob, Double_t *loc)
{
  // master2local (sensor)
  GlobalToLocal(GetModuleIndex(lay,lad,det), glob, loc);
}

//_____________________________________________________________________________________________
inline const char* AliITSUGeomTGeo::GetDetTypeName(Int_t i)
{
  if (i>=kNDetTypes) i/=kMaxSegmPerDetType; // full type is provided
  return fgkITSDetTypeName[i];
}

//_____________________________________________________________________________________________
inline const AliITSsegmentation* AliITSUGeomTGeo::GetSegmentationByID(Int_t id) const 
{
  // get segmentation by ID
  return fSegm ? (AliITSsegmentation*)fSegm->At(id) : 0;
}

//_____________________________________________________________________________________________
inline const AliITSsegmentation* AliITSUGeomTGeo::GetSegmentation(Int_t lr) const 
{
  // get segmentation of layer
  return fSegm ? (AliITSsegmentation*)fSegm->At( GetLayerDetTypeID(lr) ) : 0;
}

#endif
