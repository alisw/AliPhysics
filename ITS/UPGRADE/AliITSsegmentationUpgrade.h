#ifndef ALIITSSEGMENTATIONUPGRADE_H
#define ALIITSSEGMENTATIONUPGRADE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TArrayD.h>


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
  virtual void GetNpad(Int_t ilayer, Int_t &nx, Int_t &nz);

  // Transformation from Geant cm detector center local coordinates
  // to detector segmentation/cell coordiantes starting from (0,0).
  Bool_t  GlobalToDet(Int_t ilayer, Double_t x,Double_t y,Double_t z,Double_t &xl,Double_t &zl);
  // Transformation from detector segmentation/cell coordiantes starting
  // from (0,0) to Geant cm detector center local coordinates.
  Bool_t  DetToGlobal(Int_t ilayer, Double_t xl,Double_t zl,Double_t &x,Double_t &y, Double_t &z) const;
  
  //
  // Get Detector Segmentation Parameters
  //
    
  Double_t GetCellSizeX(Int_t ilayer){return fCellSizeX.At(ilayer);}
  Double_t GetCellSizeZ(Int_t ilayer){return fCellSizeZ.At(ilayer);}
  Double_t GetHalfLength(Int_t ilayer){return fHalfLength.At(ilayer);}
  Double_t GetRadius(Int_t ilayer) {return fMinRadius.At(ilayer);}
  
  TArrayD GetFullCellSizeX() {return fCellSizeX;}
  TArrayD GetFullCellSizeZ() {return fCellSizeZ;}
  // Pixel size in x,z 
  virtual void GetSegmentation(Int_t ilayer,  Double_t &xsize, Double_t &zsize) const;
   
  // layer thickness
  virtual Float_t GetThickness(Int_t ilayer) const {if(ilayer > fMinRadius.GetSize() || ilayer < 0) return -1; else return fMaxRadius.At(ilayer) - fMinRadius.At(ilayer);}
   
  static const Int_t GetNLayers();  
 
 protected:

  TArrayD fCellSizeX;       //Size for each pixel in x -microns
  TArrayD fCellSizeZ;       //Size for each pixel in z -microns
  TArrayD fMinRadius ;      // layer inner radius
  TArrayD fMaxRadius ;      // layer outer radius
  TArrayD fHalfLength   ;   // layer length
 private:
  AliITSsegmentationUpgrade(const AliITSsegmentationUpgrade &source);
  AliITSsegmentationUpgrade& operator=(const AliITSsegmentationUpgrade &source);

  ClassDef(AliITSsegmentationUpgrade,1) //Segmentation class for Upgrade 

    };

#endif

