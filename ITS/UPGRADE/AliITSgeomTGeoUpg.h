#ifndef ALIITSGEOMTGEOUPG_H
#define ALIITSGEOMTGEOUPG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
//  AliITSgeomTGeoUpg is a simple interface class to TGeoManager       //
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
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TGeoMatrix.h>

class TGeoPNEntry;
class TDatime;

class AliITSgeomTGeoUpg : public TObject {

 public:
  enum {kNLayersOld = 6}; // The number of layers in OLD ITS.
  enum {kITSVNA, kITSVOld, kITSVUpg}; // ITS version

  AliITSgeomTGeoUpg() { } // Default constructor
  virtual ~AliITSgeomTGeoUpg() { }// Destructor

  // This function returns the total number of ITS modules 
  static Int_t GetNModules() {CheckInit(); return fgNModules;}
  // This function returns the number of detectors/ladder for a given layer 
  static Int_t GetNDetectors(Int_t lay) {CheckInit(); return (lay<1||lay>fgNLayers) ? 0:fgNDetectors[lay-1];}
  // This function returns the number of ladders for a given layer
  static Int_t GetNLadders(Int_t lay)   {CheckInit(); return (lay<1||lay>fgNLayers) ? 0:fgNLadders[lay-1];}
  // This function returns the number of layers
  static Int_t GetNLayers()             {CheckInit(); return fgNLayers;}

  // Two methods to map module index to layer,ladder,detector indeces
  static Int_t  GetModuleIndex(Int_t lay,Int_t lad,Int_t det);
  static Bool_t GetModuleId(Int_t index,Int_t &lay,Int_t &lad,Int_t &det);

  static const char *GetSymName(Int_t index); // Get TGeoPNEntry symbolic name
  static const char *GetSymName(Int_t lay,Int_t lad,Int_t det)
    { return GetSymName(GetModuleIndex(lay,lad,det)); }
 
  // This function returns a pointer to the TGeoHMatrix (local->global)
  // of a given module index
  static TGeoHMatrix* GetMatrix(Int_t index);
  static TGeoHMatrix* GetMatrix(Int_t lay,Int_t lad,Int_t det) { return GetMatrix(GetModuleIndex(lay,lad,det)); }

  static Bool_t GetTranslation(Int_t index, Double_t t[3]);
  static Bool_t GetTranslation(Int_t lay,Int_t lad,Int_t det, Double_t t[3]) { return GetTranslation(GetModuleIndex(lay,lad,det),t); }

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
  static void   CheckInit() { if (fgVersion==kITSVNA) BuildITS(); }
  static Int_t  GetLayerDetTypeID(Int_t lr) {CheckInit();  return (lr<1||lr>fgNLayers||!fgLrDetType) ? -1:fgLrDetType[lr-1];}

  static const char* GetITSVolPattern()        {return fgkITSVolName;}
  static const char* GetITSLayerPattern()      {return fgkITSLrName;}
  static const char* GetITSLadderPattern()     {return fgkITSLadName;}
  static const char* GetITSModulePattern()     {return fgkITSModName;}
  static const char* GetITSSensorPattern()     {return fgkITSSensName;}

 private:

  static Bool_t       GetLayer(Int_t index,Int_t &lay,Int_t &index2);
  static TGeoPNEntry* GetPNEntry(Int_t index);
  static Int_t        ExtractNumberOfDetectors(const Int_t lay);
  static Int_t        ExtractNumberOfLadders(const Int_t lay);
  static Int_t        ExtractLayerDetType(const Int_t lay);
  static Int_t        ExtractNumberOfLayers();
  static Bool_t       ReadVersionString(const Char_t *str,Int_t &maj,Int_t &min,TDatime &dt);
  static void         BuildITSUpg();
  static void         BuildITSOld();
  static void         BuildITS();
 //
  AliITSgeomTGeoUpg(const AliITSgeomTGeoUpg &geom);     // Copy constructor
  AliITSgeomTGeoUpg& operator=(const AliITSgeomTGeoUpg &geom);// Assignment operator

  static Int_t  fgVersion;             // ITS Version 
  static Int_t  fgNLayers;             // number of layers
  static Int_t  fgNModules;            // The total number of modules
  static Int_t *fgNLadders;            // Array of the number of ladders/layer(layer)
  static Int_t *fgLrDetType;           // Array of layer detector types
  static Int_t *fgNDetectors;          // Array of the number of detector/ladder(layer)


  // these are hardwired settings for old ITS
  static const Int_t  fgkNModulesOld;            // The total number of modules
  static const Int_t  fgkNLaddersOld[kNLayersOld];  // Array of the number of ladders/layer(layer)
  static const Int_t  fgkNDetectorsOld[kNLayersOld];// Array of the number of detector/ladder(layer)
  //
  static const char*  fgkITSVolName;             // ITS mother volume name
  static const char*  fgkITSLrName;              // ITS Layer name
  static const char*  fgkITSLadName;             // ITS Ladder name 
  static const char*  fgkITSModName;             // ITS Module name 
  static const char*  fgkITSSensName;            // ITS Sensor name 

  ClassDef(AliITSgeomTGeoUpg, 0) // ITS geometry based on TGeo
};

#endif
