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

  static Double_t AlphaFeedback(Int_t )      {return 0.020;}                              //determines number of feedback photons updated to 9/11/04 by Di Mauro
 
  TRotMatrix* RotMatrix()          const{return fpRotMatrix;}
  TString     RotMatrixName()      const{return "rot"+fName;}
  TRotation   Rot()                const{return fRot;}
  Double_t    Rho()                const{return fCenterX3.Mag();}                                //gives  distance to chamber center in MRS
  Double_t    ThetaD()             const{return fCenterX3.Theta()*TMath::RadToDeg();}            //gives polar angle of chamber center in MRS
  Double_t    PhiD()               const{return fCenterX3.Phi()  *TMath::RadToDeg();}            //gives azimuthal angle of chamber center in MRS
  Double_t    ThetaXd()            const{return fRot.ThetaX()    *TMath::RadToDeg();}    
  Double_t    PhiXd()              const{return fRot.PhiX()      *TMath::RadToDeg();}    
  Double_t    ThetaYd()            const{return fRot.ThetaY()    *TMath::RadToDeg();}    
  Double_t    PhiYd()              const{return fRot.PhiY()      *TMath::RadToDeg();}    
  Double_t    ThetaZd()            const{return fRot.ThetaZ()    *TMath::RadToDeg();}    
  Double_t    PhiZd()              const{return fRot.PhiZ()      *TMath::RadToDeg();}    
  void        RotX(Double_t a)       {a*=TMath::DegToRad();fRot.RotateX(a);fCenterX3.RotateX(a);fRadX3.RotateX(a);fPcX3.RotateX(a);}//degrees around X
  void        RotY(Double_t a)       {a*=TMath::DegToRad();fRot.RotateY(a);fCenterX3.RotateY(a);fRadX3.RotateY(a);fPcX3.RotateY(a);}//degrees around Y
  void        RotZ(Double_t a)       {a*=TMath::DegToRad();fRot.RotateZ(a);fCenterX3.RotateZ(a);fRadX3.RotateZ(a);fPcX3.RotateZ(a);}//degrees around Z
  TVector3    Rad()               const{return fRadX3;}         //provides center of radiator position in MRS, cm   
  TVector3    Pc()                const{return fPcX3;}          //provides center of photocathond position in MRS, cm
  TVector3    Center()            const{return fCenterX3;}      //provides center of chamber position in MRS, cm
  void        Print(Option_t *sOption)const;                    //virtual interface from TObject
//Transformations for photcathode plane  
  TVector2    Mrs2Pc(TVector3 x3)const{x3-=fPcX3;x3.Transform(fRot.Inverse());return TVector2(-x3.X()+0.5*AliRICHParam::PcSizeX(),x3.Y()+0.5*AliRICHParam::PcSizeY());}
  TVector3    Pc2Mrs(TVector2 x2)const{TVector3 x3(-x2.X()+0.5*AliRICHParam::PcSizeX(),x2.Y()-0.5*AliRICHParam::PcSizeY(),0);x3.Transform(fRot); x3+=fPcX3;return x3;}  
  TVector2    Mrs2Pc(TLorentzVector x4)            const{return Mrs2Pc(x4.Vect());}
//Transformations for radiator plane  
  TVector2    Mrs2Rad(TVector3 x3)const{x3-=fRadX3;x3.Transform(fRot.Inverse());return TVector2(-x3.X()+0.5*AliRICHParam::PcSizeX(),x3.Y()+0.5*AliRICHParam::PcSizeY());}
  TVector3    Rad2Mrs(TVector2 x2)const{TVector3 x3(-x2.X()+0.5*AliRICHParam::PcSizeX(),x2.Y()-0.5*AliRICHParam::PcSizeY(),0);x3.Transform(fRot); x3+=fRadX3;return x3;}  
  TVector3    PMrs2Loc(TVector3 p3)const{TVector3 ploc=Rot().Invert()*p3;ploc.SetXYZ(-ploc.Px(),ploc.Py(),ploc.Pz()); return ploc;}  
protected:
  TVector3      fCenterX3;        //chamber center position in MRS (cm) 
  TVector3      fRadX3;           //radiator entrance center position in MRS (cm)
  TVector3      fPcX3;            //PC center position in MRS (cm)
  TRotation     fRot;             //chamber rotation in MRS
  TRotMatrix   *fpRotMatrix;      //rotation matrix of the chamber with respect to MRS 
  ClassDef(AliRICHChamber,7)      //single RICH chamber description
};//class AliRICHChamber

#endif //AliRICHChamber_h
