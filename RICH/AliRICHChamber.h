#ifndef AliRICHChamber_h
#define AliRICHChamber_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <Riostream.h>

#include <TRotMatrix.h>
#include <TVector3.h>
#include <TMath.h>

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
// ctor & dtor      
   AliRICHChamber();                                  // default ctor
   AliRICHChamber(const AliRICHChamber & Chamber){}   // copy ctor 
   ~AliRICHChamber(){}                                // dtor
// The following staff is defined in AliRICHChamber.cxx:
   void LocaltoGlobal(Float_t pos[3],Float_t Localpos[3]);//Transformation from local to global coordinates, chamber-dependant
   void GlobaltoLocal(Float_t pos[3],Float_t localpos[3]);//Transformation from Global to local coordinates, chamber-dependant 
   void GenerateTresholds();                              //Generate pad dependent tresholds
   void DisIntegration(Float_t eloss, Float_t xhit, Float_t yhit, Int_t&x, Float_t newclust[6][500], ResponseType res);// Cluster formation method
// Inline methods:   
   void    Init(Int_t id)           {fSegmentation->Init(id);} // Recalculates all the values after some of them have been changed
   
   void    SetGid(Int_t id)         {fGid=id;}           // Set and get GEANT id  
   Int_t   GetGid()            const{return fGid;}       // Get GEANT id  

   void SetRInner(Float_t rmin)     {frMin=rmin;}        // Set inner radius of sensitive volume   
   Float_t RInner()            const{return frMin;}      // Return inner radius of sensitive volume 
   
   void SetROuter(Float_t rmax)     {frMax=rmax;}        // Set outer radius of sensitive volum  
   Float_t ROuter()            const{return frMax;}      // Return outer radius of sensitive volum  
 
   void    SetZPOS(Float_t p1)      {fzPos=p1;}
   Float_t ZPosition()         const{return fzPos;}
    
   void         SetChamberTransform(Float_t x,Float_t y,Float_t z,TRotMatrix *pRotMatrix) {fX=x; fY=y; fZ=z; fpRotMatrix=pRotMatrix;}
   TRotMatrix * GetRotMatrix()                                                    const   {return fpRotMatrix;}
   Float_t      GetX()                                                            const   {return fX;}
   Float_t      GetY()                                                            const   {return fY;}
   Float_t      GetZ()                                                            const   {return fZ;}    
   Float_t      GetOffset()                                                       const   {return TMath::Sqrt(fX*fX+fY*fY+fZ*fZ);}    
    
   void              SetGeometryModel(AliRICHGeometry* pRICHGeometry)           {fGeometry=pRICHGeometry;}        
   AliRICHGeometry*  GetGeometryModel()                                    const{return fGeometry;}
   
   void              SetResponseModel(AliRICHResponse* pRICHResponse)            {fResponse=pRICHResponse;}
   AliRICHResponse*  GetResponseModel()                                     const{return fResponse;}
   
   void              SetSegmentationModel(AliSegmentation* pRICHSegmentation)   {fSegmentation=pRICHSegmentation;}
   AliSegmentation*  GetSegmentationModel(Int_t i=0)                       const{return fSegmentation;}
   
   void                  SetReconstructionModel(AliRICHClusterFinder *pRICHReconstruction)    {fReconstruction=pRICHReconstruction;}
   AliRICHClusterFinder* &GetReconstructionModel()                                            {return fReconstruction;}

   void   SigGenInit(Float_t x, Float_t y, Float_t z)   {fSegmentation->SigGenInit(x, y, z) ;}
   Int_t  SigGenCond(Float_t x, Float_t y, Float_t z)	{return fSegmentation->SigGenCond(x, y, z);}
   Int_t  Sector(Float_t x, Float_t y)                  {return fSegmentation->Sector(x, y);} // Returns number of sector containing (x,y) position    
   void   SetPadSize(Float_t p1, Float_t p2)            {fSegmentation->SetPadSize(p1,p2);}
   
   Float_t IntPH(Float_t eloss, Float_t yhit)                        {return fResponse->IntPH(eloss,yhit);}
   Float_t IntPH(Float_t yhit)                                       {return fResponse->IntPH(yhit);}
   void   SetSigmaIntegration(Float_t p)                             {fResponse->SetSigmaIntegration(p);}
   void   SetChargeSlope(Float_t p)                                  {fResponse->SetChargeSlope(p);}
   void   SetChargeSpread(Float_t p1, Float_t p2)                    {fResponse->SetChargeSpread(p1,p2);}
   void   SetMaxAdc(Float_t p)                                       {fResponse->SetMaxAdc(p);}
   void   SetSqrtKx3(Float_t p)                                      {fResponse->SetSqrtKx3(p);}
   void   SetKx2(Float_t p)                                          {fResponse->SetKx2(p);}
   void   SetKx4(Float_t p)                                          {fResponse->SetKx4(p);}
   void   SetSqrtKy3(Float_t p)                                      {fResponse->SetSqrtKy3(p);}
   void   SetKy2(Float_t p)                                          {fResponse->SetKy2(p);}
   void   SetKy4(Float_t p)                                          {fResponse->SetKy4(p);}    
   void   SetPitch(Float_t p)                                        {fResponse->SetPitch(p);}
   void   SetWireSag(Int_t p)                                        {fResponse->SetWireSag(p);}
   void   SetVoltage(Int_t p)                                        {fResponse->SetVoltage(p);}    
   
   void   SetGapThickness(Float_t thickness)                         {fGeometry->SetGapThickness(thickness);} 
   void   SetProximityGapThickness(Float_t thickness)                {fGeometry->SetProximityGapThickness(thickness);}
   void   SetQuartzLength(Float_t length)                            {fGeometry->SetQuartzLength(length);}
   void   SetQuartzWidth(Float_t width)                              {fGeometry->SetQuartzWidth(width);}
   void   SetQuartzThickness(Float_t thickness)                      {fGeometry->SetQuartzThickness(thickness);}
   void   SetOuterFreonLength(Float_t length)                        {fGeometry->SetOuterFreonLength(length);}
   void   SetOuterFreonWidth(Float_t width)                          {fGeometry->SetOuterFreonWidth(width);}
   void   SetInnerFreonLength(Float_t length)                        {fGeometry->SetInnerFreonLength(length);}
   void   SetInnerFreonWidth(Float_t width)                          {fGeometry->SetInnerFreonWidth(width);}
   void   SetFreonThickness(Float_t thickness)                       {fGeometry->SetFreonThickness(thickness);}
   
   AliRICHChamber& operator=(const AliRICHChamber& rhs){return *this;}
   
   inline virtual void Print(Option_t *sOption)const;   

private:
   Float_t frMin;                                         // Minimum Chamber size
   Float_t frMax;                                         // Maximum Chamber size 
   Int_t   fGid;                                          // Id tag 
   Float_t fzPos;                                         // z-position of this chamber

   TRotMatrix *fpRotMatrix;                               // Rotation matrix of the chamber with respect to MRS 
   Float_t fX,fY,fZ;                                      // Position of the center of the chamber in MRS (cm)

   AliSegmentation               *fSegmentation;          // Segmentation model for each chamber
   AliRICHResponse               *fResponse;              // Response model for each chamber
   AliRICHGeometry               *fGeometry;              // Geometry model for each chamber
   AliRICHClusterFinder          *fReconstruction;        // Reconstruction model for each chamber
   ClassDef(AliRICHChamber,1)                             // A single RICH chamber desription
};
    
inline void AliRICHChamber::Print(Option_t *sOption)const
{
   TObject::Print(sOption);
   cout<<"X="<<fX<<endl;   
   cout<<"Y="<<fY<<endl;
   cout<<"Z="<<fZ<<endl;
   TVector3 vector3(fX,fY,fZ);
   cout<<"Offset="<<vector3.Mag()<<endl;
   cout<<"Polar angle="<<vector3.Theta()/TMath::Pi()*180<<endl;
   cout<<"Azimithal angle="<<vector3.Phi()/TMath::Pi()*180<<endl;
}// inline void AliRICHChamber::Print(Option_t *sOPtion)
     
#endif //AliRICHChamber_h
