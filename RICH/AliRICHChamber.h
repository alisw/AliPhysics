#ifndef ALIRICHCHAMBER_H
#define ALIRICHCHAMBER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObjArray.h>
#include <TRotMatrix.h>

#include "AliRICHTresholdMap.h"
#include "AliSegmentation.h"
#include "AliRICHGeometry.h"
#include "AliRICHResponse.h"

class AliRICHClusterFinder;

typedef enum {kMip, kCerenkov} ResponseType;

class AliRICHChamber : public TObject
{
 public:
    
  Int_t                fIndexMap[50];   //indeces of tresholds
  AliRICHTresholdMap*  fTresh;          //map of tresholds

 public:
    AliRICHChamber();
    AliRICHChamber(const AliRICHChamber & Chamber);
    ~AliRICHChamber(){}
//
// Set and get GEANT id  
    Int_t   GetGid()         {return fGid;}
    void    SetGid(Int_t id) {fGid=id;}
//  
// Initialisation and z-Position
    void    Init(Int_t id);
    // Set inner radius of sensitive volume 
    void SetRInner(Float_t rmin) {frMin=rmin;}
// Set outer radius of sensitive volum  
    void SetROuter(Float_t rmax) {frMax=rmax;}  
    
// Return inner radius of sensitive volume 
    Float_t RInner()            {return frMin;}
// Return outer radius of sensitive volum  
    Float_t ROuter()            {return frMax;}

    void    SetZPOS(Float_t p1) {fzPos=p1;}
    Float_t ZPosition()         {return fzPos;}

//
//Transformation from Global to local coordinates, chamber-dependant
    void LocaltoGlobal(Float_t pos[3],Float_t Localpos[3]);
    void GlobaltoLocal(Float_t pos[3],Float_t localpos[3]); 

//Generate pad dependent tresholds

    void GenerateTresholds();

    
//Setting chamber specific rotation matrices
    
    void SetChamberTransform(Float_t Trans1,Float_t Trans2,Float_t Trans3,TRotMatrix *Matrix)
	
	{
	    fChamberMatrix=Matrix;
	    fChamberTrans[0]=Trans1;
	    fChamberTrans[1]=Trans2;
	    fChamberTrans[2]=Trans3;
	}
    
    TRotMatrix * GetRotMatrix() {return fChamberMatrix;}
    
//Configure geometry model
    void    GeometryModel(AliRICHGeometry* thisGeometry){
      fGeometry=thisGeometry;
    }
    
    
// Configure response model
    void    ResponseModel(AliRICHResponse* thisResponse);
    
    //  
// Configure segmentation model
    void    SetSegmentationModel(AliSegmentation* thisSegmentation) {
	fSegmentation = thisSegmentation;
    }
    void    SetReconstructionModel(AliRICHClusterFinder *thisReconstruction) {
	fReconstruction = thisReconstruction;
    }

//  
//  Get reference to response model
    AliRICHResponse* GetResponseModel();
//  
//  Get reference to segmentation model
    AliSegmentation*  GetSegmentationModel() {
	return fSegmentation;
    }

//  Get reference to geometry model
    AliRICHGeometry*  GetGeometryModel() {
	return fGeometry;
    }
    

    AliSegmentation*  GetSegmentationModel(Int_t i) {
	return fSegmentation;
    }
    
    //
    AliRICHClusterFinder* &GetReconstructionModel() {return fReconstruction;}

// Member function forwarding to the segmentation and response models
//
// Calculate pulse height from energy loss  
    Float_t IntPH(Float_t eloss) {return fResponse->IntPH(eloss);}
    Float_t IntPH()              {return fResponse->IntPH();}
//  
// Ask segmentation if signal should be generated  
    Int_t   SigGenCond(Float_t x, Float_t y, Float_t z)
	{
	    return fSegmentation->SigGenCond(x, y, z);
	}

// Ask segmentation sector 
    Int_t   Sector(Float_t x, Float_t y)
	{
	    return fSegmentation->Sector(x, y);
	}   
    
//
// Initialisation of segmentation for hit  
    void   SigGenInit(Float_t x, Float_t y, Float_t z)
	{
	    fSegmentation->SigGenInit(x, y, z) ;
	}
// Configuration forwarding
//
    void   SetSigmaIntegration(Float_t p)
	{
	    fResponse->SetSigmaIntegration(p);
	}
    void   SetChargeSlope(Float_t p)
	{
	    fResponse->SetChargeSlope(p);
	}
    void   SetChargeSpread(Float_t p1, Float_t p2)
	{
	    fResponse->SetChargeSpread(p1,p2);
	}
    void   SetMaxAdc(Float_t p)
	{
	    fResponse->SetMaxAdc(p);
	}
    void   SetSqrtKx3(Float_t p)
	{
	    fResponse->SetSqrtKx3(p);
	}
    void   SetKx2(Float_t p)
	{
	    fResponse->SetKx2(p);
	}
    void   SetKx4(Float_t p)
	{
	    fResponse->SetKx4(p);
	}
    void   SetSqrtKy3(Float_t p)
	{
	    fResponse->SetSqrtKy3(p);
	}
    void   SetKy2(Float_t p)
	{
	    fResponse->SetKy2(p);
	}
    void   SetKy4(Float_t p)
	{
	    fResponse->SetKy4(p);
	}
    
    void   SetPitch(Float_t p)
	{
	    fResponse->SetPitch(p);
	}
    
    void   SetPadSize(Float_t p1, Float_t p2)
	{
	    fSegmentation->SetPadSize(p1,p2);
	}
    void   SetGapThickness(Float_t thickness)
      {
	fGeometry->SetGapThickness(thickness);
      } 
    void   SetProximityGapThickness(Float_t thickness)
      {
	fGeometry->SetProximityGapThickness(thickness);
      }
    void   SetQuartzLength(Float_t length)
      {
	fGeometry->SetQuartzLength(length);
      }
    void   SetQuartzWidth(Float_t width)
      {
	fGeometry->SetQuartzWidth(width);
      }
    void   SetQuartzThickness(Float_t thickness) 
      {
	fGeometry->SetQuartzThickness(thickness);
      }
    void   SetOuterFreonLength(Float_t length)
      {
	fGeometry->SetOuterFreonLength(length);
      }
    void   SetOuterFreonWidth(Float_t width)
      {
	fGeometry->SetOuterFreonWidth(width);
      }
    void   SetInnerFreonLength(Float_t length)
      {
	fGeometry->SetInnerFreonLength(length);
      }
    void   SetInnerFreonWidth(Float_t width) 
      {
	fGeometry->SetInnerFreonWidth(width);
      }
    void   SetFreonThickness(Float_t thickness)
      {
	fGeometry->SetFreonThickness(thickness);
      }

    AliRICHChamber& operator=(const AliRICHChamber& rhs);
    
//  
// Cluster formation method
    void   DisIntegration(Float_t eloss, Float_t xhit, Float_t yhit, Int_t&x, Float_t newclust[6][500], ResponseType res);
 private:
// GEANT volume if for sensitive volume of this
    Float_t frMin;                 // Minimum Chamber size
    Float_t frMax;                 // Maximum Chamber size 
    Int_t   fGid;                  // Id tag 
    Float_t fzPos;                 // z-position of this chamber

    TRotMatrix *fChamberMatrix;          //Rotation matrices for each chamber
    Float_t fChamberTrans[3];            //Translaction vectors for each chamber

    AliSegmentation               *fSegmentation;          //Segmentation model for each chamber
    AliRICHResponse               *fResponse;              //Response model for each chamber
    AliRICHGeometry               *fGeometry;              //Geometry model for each chamber
    AliRICHClusterFinder          *fReconstruction;        //Reconstruction model for each chamber
    ClassDef(AliRICHChamber,1)
};
#endif




