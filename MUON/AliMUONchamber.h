#ifndef MUONchamber_H
#define MUONchamber_H
#include "TObjArray.h"
#include "AliMUONSegRes.h"

class AliMUONClusterFinder;
class AliMUONresponse ;
class AliMUONsegmentation ;

class AliMUONchamber:
public TObject
{
 public:
    AliMUONchamber();
    ~AliMUONchamber(){}
//
// Set and get GEANT id  
  Int_t   GetGid()         {return fGid;}
  void    SetGid(Int_t id) {fGid=id;}
//  
// Initialisation and z-Position
  void    Init();
  void    SetZPOS(Float_t p1) {fzPos=p1;}
  Float_t ZPosition()         {return fzPos;}
// Set inner radius of sensitive volume 
  void SetRInner(Float_t rmin) {frMin=rmin;}
// Set outer radius of sensitive volum  
  void SetROuter(Float_t rmax) {frMax=rmax;}  

// Return inner radius of sensitive volume 
  Float_t RInner()            {return frMin;}
// Return outer radius of sensitive volum  
  Float_t ROuter()            {return frMax;}  
//  
// Configure response model
  void    ResponseModel(AliMUONresponse* thisResponse) {fResponse=thisResponse;}
//  
// Configure segmentation model
  void    SegmentationModel(Int_t i, AliMUONsegmentation* thisSegmentation) {
      (*fSegmentation)[i-1] = thisSegmentation;
  }
  void    ReconstructionModel(AliMUONClusterFinder *thisReconstruction) {
      fReconstruction = thisReconstruction;
  }
  
//  
//  Get reference to response model
  AliMUONresponse*          &GetResponseModel(){return fResponse;}
//
  AliMUONClusterFinder*     &GetReconstructionModel(){return fReconstruction;}
//  
//  Get reference to segmentation model
  AliMUONsegmentation*  GetSegmentationModel(Int_t isec) {
      return (AliMUONsegmentation *) (*fSegmentation)[isec-1];
  }
  TObjArray* GetChamberSegmentation(){return fSegmentation;}
  
  Int_t Nsec()              {return fnsec;}
  void  SetNsec(Int_t nsec) {fnsec=nsec;}
//
// Member function forwarding to the segmentation and response models
//
// Calculate pulse height from energy loss  
  Float_t IntPH(Float_t eloss) {return fResponse->IntPH(eloss);}
//  
// Ask segmentation if signal should be generated  
  Int_t   SigGenCond(Float_t x, Float_t y, Float_t z)
      {
	  if (fnsec==1) {
	      return ((AliMUONsegmentation*) (*fSegmentation)[0])
		  ->SigGenCond(x, y, z) ;
	  } else {
	      return (((AliMUONsegmentation*) (*fSegmentation)[0])
		      ->SigGenCond(x, y, z)) ||
		  (((AliMUONsegmentation*) (*fSegmentation)[1])
		   ->SigGenCond(x, y, z)) ;
	  }
  }
//
// Initialisation of segmentation for hit  
  void    SigGenInit(Float_t x, Float_t y, Float_t z)
      {
	  
	  if (fnsec==1) {
	      ((AliMUONsegmentation*) (*fSegmentation)[0])->SigGenInit(x, y, z) ;
	  } else {
	      ((AliMUONsegmentation*) (*fSegmentation)[0])->SigGenInit(x, y, z) ;
	      ((AliMUONsegmentation*) (*fSegmentation)[1])->SigGenInit(x, y, z) ;
	  }
      }

// Configuration forwarding
//
  void   SetSigmaIntegration(Float_t p1)         {fResponse->SetSigmaIntegration(p1);}
  void   SetChargeSlope(Float_t p1)              {fResponse->SetChargeSlope(p1);}
  void   SetChargeSpread(Float_t p1, Float_t p2) {fResponse->SetChargeSpread(p1,p2);}
  void   SetMaxAdc(Float_t p1)                   {fResponse->SetMaxAdc(p1);}

  void   SetPADSIZ(Int_t isec, Float_t p1, Float_t p2) {
      ((AliMUONsegmentation*) (*fSegmentation)[isec-1])->SetPADSIZ(p1,p2);
  }
//  
// Cluster formation method
  void   DisIntegration(Float_t, Float_t, Float_t, Int_t&x, Float_t newclust[6][500]);
    ClassDef(AliMUONchamber,1)
   void    InitGeo(Float_t z);

 private:
// GEANT volume if for sensitive volume of this chamber
  Int_t   fGid;
// z-position of this chamber
  Float_t fzPos; // z-position of chambers
  Int_t   fnsec; // number of segmentation zones
  Float_t frMin; // innermost sensitive radius
  Float_t frMax; // outermost sensitive radius
// The segmentation models for the cathode planes
// fnsec=1: one plane segmented, fnsec=2: both planes are segmented.

  TObjArray            *fSegmentation;
  AliMUONClusterFinder *fReconstruction;
  AliMUONresponse      *fResponse;
  
 public:
  Float_t fdGas; // half gaz gap
  Float_t fdAlu; // half Alu width
};

#endif
