#ifndef AliRICHChamber_h
#define AliRICHChamber_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TVector3.h>
#include <TMath.h>
#include <TRotation.h>
#include <TLorentzVector.h>
#include "AliRICHParam.h"
#include "AliSegmentation.h"
class AliRICHGeometry;
class AliRICHResponse;
class TRotMatrix;

typedef enum {kMip, kPhoton} ResponseType;
class AliRICHParam;

class AliRICHChamber : public TNamed
{
public:
           AliRICHChamber();
           AliRICHChamber(Int_t iModuleN,AliRICHParam *pParam);
           AliRICHChamber(const AliRICHChamber &chamber):TNamed(chamber) {;}
  virtual ~AliRICHChamber()                                              {;}
           AliRICHChamber& operator=(const AliRICHChamber&)              {return *this;}
  
  TRotMatrix* RotMatrix()          const{return fpRotMatrix;}
  TString RotMatrixName()          const{return "rot"+fName;}
  TRotation   Rot()                const{return fRot;}
  Double_t    Rho()                const{return fCenterV3.Mag();} 
  Double_t    ThetaD()             const{return fCenterV3.Theta()*TMath::RadToDeg();}    
  Double_t    PhiD()               const{return fCenterV3.Phi()*TMath::RadToDeg();}    
  Double_t    ThetaXd()            const{return fRot.ThetaX()*TMath::RadToDeg();}    
  Double_t    PhiXd()              const{return fRot.PhiX()*TMath::RadToDeg();}    
  Double_t    ThetaYd()            const{return fRot.ThetaY()*TMath::RadToDeg();}    
  Double_t    PhiYd()              const{return fRot.PhiY()*TMath::RadToDeg();}    
  Double_t    ThetaZd()            const{return fRot.ThetaZ()*TMath::RadToDeg();}    
  Double_t    PhiZd()              const{return fRot.PhiZ()*TMath::RadToDeg();}    
  void        RotateX(Double_t a)       {fRot.RotateX(a);fCenterV3.RotateX(a);fPcX3.RotateX(a);}
  void        RotateY(Double_t a)       {fRot.RotateY(a);fCenterV3.RotateY(a);fPcX3.RotateY(a);}
  void        RotateZ(Double_t a)       {fRot.RotateZ(a);fCenterV3.RotateZ(a);fPcX3.RotateZ(a);}
  Double_t    X()                  const{return fCenterV3.X();}  
  Double_t    Y()                  const{return fCenterV3.Y();}   
  Double_t    Z()                  const{return fCenterV3.Z();}
  TVector3    L2G(TVector3 x3)                       const{x3.Transform(fRot);x3+=fCenterV3;return x3;}
  TVector3    G2L(TVector3 x3)                       const{x3-=fCenterV3;x3.Transform(fRot.Inverse()); return x3;}
  inline TVector2  Glob2Loc(TVector3 x3, Bool_t isVector=kFALSE) const;
  TVector2    Glob2Loc(TLorentzVector x4,Bool_t isVector=kFALSE) const{return Glob2Loc(x4.Vect(),isVector);}
  TVector3    L2G(Double_t x,Double_t y,Double_t z)  const{return L2G(TVector3(x,y,z));}
  TVector3    G2L(TLorentzVector x4)                 const{return G2L(x4.Vect());}
  Float_t     G2Ly(TLorentzVector x4)                const{TVector3 x3=G2L(x4.Vect()); return x3.Z();}
  TVector3    G2L(Double_t x,Double_t y,Double_t z)  const{return G2L(TVector3(x,y,z));}
  Float_t     G2Lx(Double_t x,Double_t y,Double_t z) const{TVector3 x3=G2L(x,y,z); return x3.X();}
  Float_t     G2Ly(Double_t x,Double_t y,Double_t z) const{TVector3 x3=G2L(x,y,z); return x3.Z();}
  void        Print(Option_t *sOption)const;//virtual      
   
  void LocaltoGlobal(Float_t pos[3],Float_t localpos[3]) {
    //Transformation from local to global coordinates, chamber-dependant
    TVector3 buf = L2G(localpos[0],localpos[1],localpos[2]);
    pos[0]=buf.X();pos[1]=buf.Y();pos[2]=buf.Z();
  }
  void GlobaltoLocal(Float_t pos[3],Float_t localpos[3]) {
    //Transformation from Global to local coordinates, chamber-dependant 
    TVector3 buf = G2L(pos[0],pos[1],pos[2]);
    localpos[0]=buf.X();localpos[1]=buf.Y();localpos[2]=buf.Z();
  }
  inline void SetToZenith();
  TRotMatrix *GetRotMatrix()       const{return fpRotMatrix;}  
protected:
  TVector3      fCenterV3;        //chamber center position in MRS (cm)
  TVector3      fPcX3;            //PC center position in MRS (cm)
  TRotation     fRot;             //chamber rotation in MRS
  TRotMatrix   *fpRotMatrix;      //rotation matrix of the chamber with respect to MRS 
  AliRICHParam *fpParam;          //main RICH parameters description  
  ClassDef(AliRICHChamber,5)      //single RICH chamber description
};//class AliRICHChamber
//__________________________________________________________________________________________________
void AliRICHChamber::SetToZenith()
{
  fCenterV3.SetXYZ(0,AliRICHParam::Offset()-AliRICHParam::GapThickness()/2,0); 
  fPcX3.SetXYZ(0,AliRICHParam::Offset()-AliRICHParam::GapThickness()/2+5.276+0.25,0);   
}
//__________________________________________________________________________________________________
TVector2 AliRICHChamber::Glob2Loc(TVector3 x3,Bool_t isVector)const
{
  if(!isVector) x3-=fPcX3;
  x3.Transform(fRot.Inverse()); 
  return TVector2(x3.X(),x3.Z());//attention Y and Z are misplaced!
}
//__________________________________________________________________________________________________  
#endif //AliRICHChamber_h
