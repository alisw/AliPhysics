#ifndef AliRICHChamber_h
#define AliRICHChamber_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TVector3.h>
#include <TMath.h>
#include <TRotation.h>
#include <TLorentzVector.h>
#include "AliRICHParam.h"
class TRotMatrix;


class AliRICHChamber : public TNamed
{
public:
           AliRICHChamber():TNamed(),fpRotMatrix(0)                      {;}
           AliRICHChamber(Int_t iChamberN);
           AliRICHChamber(const AliRICHChamber &chamber):TNamed(chamber) {;}
  virtual ~AliRICHChamber()                                              {;}
           AliRICHChamber& operator=(const AliRICHChamber&)              {return *this;}
  
  TRotMatrix* RotMatrix()          const{return fpRotMatrix;}
  TString RotMatrixName()          const{return "rot"+fName;}
  TRotation   Rot()                const{return fRot;}
  Double_t    Rho()                const{return fCenterV3.Mag();}                                //gives  distance to chamber center in MRS
  Double_t    ThetaD()             const{return fCenterV3.Theta()*TMath::RadToDeg();}            //gives polar angle of chamber center in MRS
  Double_t    PhiD()               const{return fCenterV3.Phi()  *TMath::RadToDeg();}            //gives azimuthal angle of chamber center in MRS
  Double_t    ThetaXd()            const{return fRot.ThetaX()    *TMath::RadToDeg();}    
  Double_t    PhiXd()              const{return fRot.PhiX()      *TMath::RadToDeg();}    
  Double_t    ThetaYd()            const{return fRot.ThetaY()    *TMath::RadToDeg();}    
  Double_t    PhiYd()              const{return fRot.PhiY()      *TMath::RadToDeg();}    
  Double_t    ThetaZd()            const{return fRot.ThetaZ()    *TMath::RadToDeg();}    
  Double_t    PhiZd()              const{return fRot.PhiZ()      *TMath::RadToDeg();}    
  void        RotateX(Double_t a)       {fRot.RotateX(a);fCenterV3.RotateX(a);fPcX3.RotateX(a);} //rotate chamber around X by "a" degrees
  void        RotateY(Double_t a)       {fRot.RotateY(a);fCenterV3.RotateY(a);fPcX3.RotateY(a);} //rotate chamber around Y by "a" degrees
  void        RotateZ(Double_t a)       {fRot.RotateZ(a);fCenterV3.RotateZ(a);fPcX3.RotateZ(a);} //rotate chamber around Z by "a" degrees
  Double_t    X()                  const{return fCenterV3.X();}  
  Double_t    Y()                  const{return fCenterV3.Y();}   
  Double_t    Z()                  const{return fCenterV3.Z();}
  TVector2    Glob2Loc(TVector3 x3)const{x3-=fPcX3;x3.Transform(fRot.Inverse());return TVector2(x3.Z()+0.5*AliRICHParam::PcSizeX(),-x3.X()+0.5*AliRICHParam::PcSizeY());}//Y and Z are misplaced?????
  TVector3    Loc2Glob(TVector2 x2)const{TVector3 x3(-x2.Y()+0.5*AliRICHParam::PcSizeY(),0,x2.X()-0.5*AliRICHParam::PcSizeX());x3.Transform(fRot); x3+=fPcX3;return x3;}
  
  TVector2    Glob2Loc(TLorentzVector x4)            const{return Glob2Loc(x4.Vect());}
  
  void        Print(Option_t *sOption)const;//virtual      
   

  inline void SetToZenith();
  TRotMatrix *GetRotMatrix()       const{return fpRotMatrix;}  
protected:
  TVector3      fCenterV3;        //chamber center position in MRS (cm) 
  TVector3      fPcX3;            //PC center position in MRS (cm)
  TRotation     fRot;             //chamber rotation in MRS
  TRotMatrix   *fpRotMatrix;      //rotation matrix of the chamber with respect to MRS 
  ClassDef(AliRICHChamber,6)      //single RICH chamber description
};//class AliRICHChamber
//__________________________________________________________________________________________________
void AliRICHChamber::SetToZenith()
{
//Put the chamber to zenith. Position of PC is shifted in X-Z plane since the origin of chamber local system is in
//left hand down coner.     
  fCenterV3.SetXYZ(0,AliRICHParam::Offset()-AliRICHParam::GapThickness()/2           ,0); 
      fPcX3.SetXYZ(0,AliRICHParam::Offset()-AliRICHParam::GapThickness()/2+5.276+0.25,0);   
}
//__________________________________________________________________________________________________
#endif //AliRICHChamber_h
