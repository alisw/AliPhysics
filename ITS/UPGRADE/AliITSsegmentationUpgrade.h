#ifndef ALIITSSEGMENTATIONUPGRADE_H
#define ALIITSSEGMENTATIONUPGRADE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TArrayD.h>
#include <TObject.h>

//////////////////////////////////////////////////
//  Authors A.Mastroserio			//
//          C.Terrevoli				//
//          annalisa.mastroserio@cern.ch	//	
//	    cristina.terrevoli@ba.infn.it       //
// ITS Upgrade segmentation virtual base class  //
//                                              //
//////////////////////////////////////////////////


class AliITSsegmentationUpgrade : public TObject {
 public:

  AliITSsegmentationUpgrade();
  AliITSsegmentationUpgrade(TArrayD Radii, TArrayD widths, TArrayD Length);
  virtual ~AliITSsegmentationUpgrade(){}//dtor
    

  // Set Detector Segmentation Parameters
  virtual void SetSegmentation(Int_t ilayer, Double_t xsize, Double_t zsize); // x/z size in microns
  virtual void SetFullSegmentation(TArrayD xsize, TArrayD zsize); // x/z size in microns
  virtual void SetNSectors(Int_t nSect) {fNSectors=nSect;}
  virtual void GetNpad(Int_t ilayer, Int_t &nx, Int_t &nz);

  // Transformation from Geant cm detector center local coordinates
  // to detector segmentation/cell coordiantes starting from (0,0).
  Bool_t  GlobalToDet(Int_t ilayer, Double_t x,Double_t y,Double_t z,Double_t &xl,Double_t &zl) const;
  Bool_t  GlobalToDet(Int_t ilayer, Double_t x,Double_t y,Double_t z,Double_t &xl,Double_t &zl, Int_t &module) const;
  
  // Transformation from detector segmentation/cell coordiantes starting
  // from (0,0) to Geant cm detector center local coordinates.
  Bool_t  DetToGlobal(Int_t ilayer, Double_t xl,Double_t zl,Double_t &x,Double_t &y, Double_t &z) const;
  Bool_t  DetToGlobal(Int_t ilayer, Int_t module, Double_t xl,Double_t zl,Double_t &x,Double_t &y, Double_t &z) const;
  Bool_t  DetToPixID(Double_t xl, Double_t zl,Int_t layer, Int_t &nx, Int_t &nz) const;  
  Bool_t  DetToTrack(Int_t layer,Int_t module, Double_t xl, Double_t zl, Double_t &ytr, Double_t &ztr)const;
  Bool_t  DetToTrack2(Int_t layer,Int_t module, Double_t xl, Double_t zl, Double_t &ytr, Double_t &ztr)const;
   //
  // Get Detector Segmentation Parameters
  //
    
  Int_t GetIdIndex(Int_t layer, Int_t sector) const {return sector*100 + layer; }
  Int_t GetLayerFromIdIndex(Int_t id)const {return id%100; }
  Int_t GetModuleFromIdIndex(Int_t id)const {return id/100; }
  
  Double_t GetCellSizeX(Int_t ilayer){return fCellSizeX.At(ilayer);}
  Double_t GetCellSizeZ(Int_t ilayer){return fCellSizeZ.At(ilayer);}
  Double_t GetHalfLength(Int_t ilayer){return fHalfLength.At(ilayer);}
  Double_t GetRadius(Int_t ilayer) {return fMinRadius.At(ilayer);}
  Double_t GetAlpha(Int_t module) const; 
  Int_t GetModule(Double_t phi)const;  
  Int_t GetModule(Double_t x, Double_t y)const;  
  TArrayD GetFullCellSizeX() {return fCellSizeX;}
  TArrayD GetFullCellSizeZ() {return fCellSizeZ;}
  // Pixel size in x,z 
  virtual void GetSegmentation(Int_t ilayer,  Double_t &xsize, Double_t &zsize) const;
  Int_t GetNSectors() {return fNSectors;} 
  // layer thickness
  virtual Float_t GetThickness(Int_t ilayer) const {if(ilayer > fMinRadius.GetSize() || ilayer < 0) return -1; else return fMaxRadius.At(ilayer) - fMinRadius.At(ilayer);}
   
  static Int_t GetNLayers();  
 
 protected:

  TArrayD fCellSizeX;       //Size for each pixel in x -microns
  TArrayD fCellSizeZ;       //Size for each pixel in z -microns
  Int_t   fNSectors;
  TArrayD fMinRadius ;      // layer inner radius
  TArrayD fMaxRadius ;      // layer outer radius
  TArrayD fHalfLength   ;   // layer length
 private:
  AliITSsegmentationUpgrade(const AliITSsegmentationUpgrade &source);
  AliITSsegmentationUpgrade& operator=(const AliITSsegmentationUpgrade &source);

  ClassDef(AliITSsegmentationUpgrade,2) //Segmentation class for Upgrade 

    };

#endif

