#ifndef AliRICHChamber_h
#define AliRICHChamber_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TRotMatrix.h>
#include <TVector3.h>
#include <TMath.h>
#include <TRotation.h>
#include "AliRICHConst.h"

#include "AliRICHTresholdMap.h"
#include "AliSegmentation.h"
#include "AliRICHGeometry.h"
#include "AliRICHResponse.h"

class AliRICHClusterFinder;

typedef enum {kMip, kPhoton} ResponseType;
class AliRICHParam;

class AliRICHChamber : public TNamed
{
public:
    
   Int_t                fIndexMap[50];   //indeces of tresholds
   AliRICHTresholdMap*  fTresh;          //map of tresholds

public:
           AliRICHChamber();
           AliRICHChamber(Int_t iModuleN,AliRICHParam *pParam);
           AliRICHChamber(const AliRICHChamber &chamber):TNamed(chamber) {;}
  virtual ~AliRICHChamber()                                              {;}
           AliRICHChamber& operator=(const AliRICHChamber&)              {return *this;}
  
  TRotMatrix* RotMatrix()          const{return fpRotMatrix;}
  const char* RotMatrixName()      const{return "rot"+fName;}
  TRotation   Rot()                     {return fRot;}
  Double_t    Rho()                const{return fCenterV3.Mag();} 
  Double_t    ThetaD()             const{return fCenterV3.Theta()*kR2d;}    
  Double_t    PhiD()               const{return fCenterV3.Phi()*kR2d;}    
  Double_t    ThetaXd()            const{return fRot.ThetaX()*kR2d;}    
  Double_t    PhiXd()              const{return fRot.PhiX()*kR2d;}    
  Double_t    ThetaYd()            const{return fRot.ThetaY()*kR2d;}    
  Double_t    PhiYd()              const{return fRot.PhiY()*kR2d;}    
  Double_t    ThetaZd()            const{return fRot.ThetaZ()*kR2d;}    
  Double_t    PhiZd()              const{return fRot.PhiZ()*kR2d;}    
  void        RotateX(Double_t a)       {fRot.RotateX(a);fCenterV3.RotateX(a);}
  void        RotateY(Double_t a)       {fRot.RotateY(a);fCenterV3.RotateY(a);}
  void        RotateZ(Double_t a)       {fRot.RotateZ(a);fCenterV3.RotateZ(a);}
  Double_t    X()                  const{return fCenterV3.X();}  
  Double_t    Y()                  const{return fCenterV3.Y();}   
  Double_t    Z()                  const{return fCenterV3.Z();}
  TVector3    L2G(TVector3 v3)                      const{v3.Transform(fRot.Inverse());v3+=fCenterV3;return v3;}
  TVector3    L2G(Double_t x,Double_t y,Double_t z) const{TVector3 v3(x,y,z);return L2G(v3);}
  TVector3    G2L(TVector3 v3)     const{v3-=fCenterV3;  v3.Transform(fRot);return v3;}
  TVector3    G2L(Double_t x,Double_t y,Double_t z) const{TVector3 v3(x,y,z);return G2L(v3);}
//  TLorentzVector L2G(TLorentzVector v4)  const{v4.Transform(fRot.Inverse());v4+=fCenterV3;return v4;}???
  void        Print(Option_t *sOption)const;//virtual      
   
  void LocaltoGlobal(Float_t pos[3],Float_t Localpos[3]);//Transformation from local to global coordinates, chamber-dependant
  void GlobaltoLocal(Float_t pos[3],Float_t localpos[3]);//Transformation from Global to local coordinates, chamber-dependant 
  void GenerateTresholds();                              //Generate pad dependent tresholds
  void DisIntegration(Float_t eloss, Float_t xhit, Float_t yhit, Int_t&x, Float_t newclust[6][500], ResponseType res);// Cluster formation method
  void    Init(Int_t id)           {fSegmentation->Init(id);} // Recalculates all the values after some of them have been changed
  void              SetGeometryModel(AliRICHGeometry* pRICHGeometry)            {fGeometry=pRICHGeometry;}        
  AliRICHGeometry*  GetGeometryModel()                                     const{return fGeometry;}
  void              SetResponseModel(AliRICHResponse* pRICHResponse)            {fResponse=pRICHResponse;}
  AliRICHResponse*  GetResponseModel()                                     const{return fResponse;}
  void              SetSegmentationModel(AliSegmentation* pRICHSegmentation)    {fSegmentation=pRICHSegmentation;}
  AliSegmentation*  GetSegmentationModel()                                 const{return fSegmentation;}
  void                  SetReconstructionModel(AliRICHClusterFinder *pRICHReconstruction)    {fReconstruction=pRICHReconstruction;}
  AliRICHClusterFinder* &GetReconstructionModel()                                            {return fReconstruction;}
  void   SigGenInit(Float_t x, Float_t y, Float_t z)   {fSegmentation->SigGenInit(x, y, z) ;}
  Int_t  SigGenCond(Float_t x, Float_t y, Float_t z)   {return fSegmentation->SigGenCond(x, y, z);}
  Int_t  Sector(Float_t x, Float_t y)                  {return fSegmentation->Sector((Int_t)x, (Int_t)y);} // Returns number of sector containing (x,y) position    
  void   SetPadSize(Float_t p1, Float_t p2)            {fSegmentation->SetPadSize(p1,p2);}
  Double_t    GetX()               const{return fX;}
  Double_t    GetY()               const{return fY;}
  Double_t    GetZ()               const{return fZ;}    
  inline void SetCenter(Double_t x,Double_t y,Double_t z);
  TRotMatrix *GetRotMatrix()       const{return fpRotMatrix;}
  void        SetChamberTransform(Float_t x,Float_t y,Float_t z,TRotMatrix *pRotMatrix) {fX=x; fY=y; fZ=z; fpRotMatrix=pRotMatrix;}
  
protected:
  Float_t fX,fY,fZ;                                      // Position of the center of the chamber in MRS (cm)

  AliSegmentation               *fSegmentation;          //???Segmentation model for each chamber
  AliRICHResponse               *fResponse;              //???Response model for each chamber
  AliRICHGeometry               *fGeometry;              //???Geometry model for each chamber
  AliRICHClusterFinder          *fReconstruction;        //???Reconstruction model for each chamber
   
  TVector3      fCenterV3;        //chamber center position in MRS (cm)
  TRotation     fRot;             //chamber rotation in MRS
  TRotMatrix   *fpRotMatrix;      //rotation matrix of the chamber with respect to MRS 
  AliRICHParam *fpParam;          //main RICH parameters description  
  ClassDef(AliRICHChamber,2)      //single RICH chamber description
};//class AliRICHChamber

void AliRICHChamber::SetCenter(Double_t x,Double_t y,Double_t z)
{
  fCenterV3.SetXYZ(x,y,z);
  fX=x;fY=y;fZ=z;
}
  
#endif //AliRICHChamber_h
