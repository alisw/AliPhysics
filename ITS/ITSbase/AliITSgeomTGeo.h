#ifndef ALIITSGEOMTGEO_H
#define ALIITSGEOMTGEO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
//  AliITSgeomTGeo is a simple interface class to TGeoManager          //
//  It is used in the simulation and reconstruction in order to        //
//  query the TGeo ITS geometry                                        //
//                                                                     //
//  author - cvetan.cheshkov@cern.ch                                   //
//  15/02/2007                                                         //
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TGeoMatrix.h>

class TGeoPNEntry;

class AliITSgeomTGeo : public TObject {

 public:

  AliITSgeomTGeo() { } // Default constructor
  virtual ~AliITSgeomTGeo() { } // Destructor

  // This function returns the total number of ITS modules 
  static Int_t GetNModules() {return fgkNModules;}
  // This function returns the number of detectors/ladder for a given layer 
  static Int_t GetNDetectors(Int_t lay) {return fgkNDetectors[lay-1];}
  // This function returns the number of ladders for a given layer
  static Int_t GetNLadders(Int_t lay)   {return fgkNLadders[lay-1];}
  // This function returns the number of layers
  static Int_t GetNLayers()             {return kNLayers;}

  // Two methods to map module index to layer,ladder,detector indeces
  static Int_t  GetModuleIndex(Int_t lay,Int_t lad,Int_t det);
  static Bool_t GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det);

  static const char *GetSymName(Int_t index); // Get TGeoPNEntry symbolic name
  static const char *GetSymName(Int_t lay,Int_t lad,Int_t det)
    { return GetSymName(GetModuleIndex(lay,lad,det)); }
 
  // This function returns a pointer to the TGeoHMatrix (local->global)
  // of a given module index
  static TGeoHMatrix* GetMatrix(Int_t index);
  static TGeoHMatrix* GetMatrix(Int_t lay,Int_t lad,Int_t det)
    { return GetMatrix(GetModuleIndex(lay,lad,det)); }

  static Bool_t GetTranslation(Int_t index, Double_t t[3]);
  static Bool_t GetTranslation(Int_t lay,Int_t lad,Int_t det, Double_t t[3])
    { return GetTranslation(GetModuleIndex(lay,lad,det),t); }

  static Bool_t GetRotation(Int_t index, Double_t r[9]);
  static Bool_t GetRotation(Int_t lay,Int_t lad,Int_t det, Double_t r[9])
    { return GetRotation(GetModuleIndex(lay,lad,det),r); }

  // This function returns a pointer to the original TGeoHMatrix (local->global)
  // for a specific module index
  static Bool_t GetOrigMatrix(Int_t index, TGeoHMatrix &m);
  static Bool_t GetOrigMatrix(Int_t lay,Int_t lad,Int_t det, TGeoHMatrix &m)
    { return GetOrigMatrix(GetModuleIndex(lay,lad,det),m); }

  static Bool_t GetOrigTranslation(Int_t index, Double_t t[3]);
  static Bool_t GetOrigTranslation(Int_t lay,Int_t lad,Int_t det, Double_t t[3])
    { return GetOrigTranslation(GetModuleIndex(lay,lad,det),t); }

  static Bool_t GetOrigRotation(Int_t index, Double_t r[9]);
  static Bool_t GetOrigRotation(Int_t lay,Int_t lad,Int_t det, Double_t r[9])
    { return GetOrigRotation(GetModuleIndex(lay,lad,det),r); }

  static const TGeoHMatrix* GetTracking2LocalMatrix(Int_t index);
  static const TGeoHMatrix* GetTracking2LocalMatrix(Int_t lay,Int_t lad,Int_t det)
    { return GetTracking2LocalMatrix(GetModuleIndex(lay,lad,det)); }

  static Bool_t GetTrackingMatrix(Int_t index, TGeoHMatrix &m);
  static Bool_t GetTrackingMatrix(Int_t lay,Int_t lad,Int_t det, TGeoHMatrix &m)
    { return GetTrackingMatrix(GetModuleIndex(lay,lad,det),m); }

  static Bool_t LocalToGlobal(Int_t index, const Double_t *loc, Double_t *glob);
  static Bool_t LocalToGlobal(Int_t lay, Int_t lad, Int_t det,
			      const Double_t *loc, Double_t *glob)
    { return LocalToGlobal(GetModuleIndex(lay,lad,det), loc, glob);}

  static Bool_t GlobalToLocal(Int_t index, const Double_t *glob, Double_t *loc);
  static Bool_t GlobalToLocal(Int_t lay, Int_t lad, Int_t det,
			      const Double_t *glob, Double_t *loc)
    { return GlobalToLocal(GetModuleIndex(lay,lad,det), glob, loc);}

  static Bool_t LocalToGlobalVect(Int_t index, const Double_t *loc, Double_t *glob);
  static Bool_t GlobalToLocalVect(Int_t index, const Double_t *glob, Double_t *loc);

  enum {kNLayers = 6}; // The number of layers.

 private:

  static Bool_t       GetLayer(Int_t index,Int_t &lay,Int_t &index2);
  static TGeoPNEntry* GetPNEntry(Int_t index);

  AliITSgeomTGeo(const AliITSgeomTGeo &geom);     // Copy constructor
  AliITSgeomTGeo& operator=(const AliITSgeomTGeo &geom);// Assignment operator

  static const Int_t  fgkNModules;            // The total number of modules
  static const Int_t  fgkNLadders[kNLayers];  // Array of the number of ladders/layer(layer)
  static const Int_t  fgkNDetectors[kNLayers];// Array of the number of detector/ladder(layer)

  ClassDef(AliITSgeomTGeo, 0) // ITS geometry based on TGeo
};

#endif
