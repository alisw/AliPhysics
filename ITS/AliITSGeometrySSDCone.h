#ifndef ALIITSGEOMETRYSSDCONE_H
#define ALIITSGEOMETRYSSDCONE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/*
  ITS SSD Cone Geometry. Version 11
*/
#include "AliITSBaseGeometry.h"
class TVector3;
class AliITS;

class  AliITSGeometrySSDCone : public AliITSBaseGeometry {
 public:
    AliITSGeometrySSDCone();
    AliITSGeometrySSDCone(AliITS *its,TVector3 &tran,const char *moth,
			  Int_t mat0);
    virtual ~AliITSGeometrySSDCone(){};
    void CreateG3Geometry(const char *moth,TVector3 &trans);
    void PositionG3Geometry(AliITSBaseVolParams &moth,Int_t cn,
			    TVector3 &t,Int_t irot);
    void CreateG3Materials();
    void BuildDisplayGeometry();
    void Print(ostream *os); // Prints out the contenes of this class
    void Read(istream *is); // Reads in the contenst of this class
    // mother volume.
 private: // functions Require at parts of Volume A to be already defined.
    // Retruns the value of Rmax corresponding to point z alone the line
    // defined by the two points p.Rmax(i1),p.ZAt(i1) and p.Rmax(i2),
    // p.ZAt(i2).
    Double_t RmaxFrom2Points(AliITSPConeData &p,Int_t i1,Int_t i2,Double_t z)
	{return p.Rmax(i2)+(p.Rmax(i1)-p.Rmax(i2))*(z-p.ZAt(i2))/
	     (p.ZAt(i1)-p.ZAt(i2));}
    // Retruns the value of Rmin corresponding to point z alone the line
    // defined by the two points p.Rmin(i1),p.ZAt(i1) and p.Rmin(i2),
    // p.ZAt(i2).
    Double_t RminFrom2Points(AliITSPConeData &p,Int_t i1,Int_t i2,Double_t z)
	{return p.Rmin(i2)+(p.Rmin(i1)-p.Rmin(i2))*(z-p.ZAt(i2))/

	     (p.ZAt(i1)-p.ZAt(i2));}
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p.Rmin(i1),p.ZAt(i1) and p.Rmin(i2),p.ZAt(i2)
    Double_t Zfrom2MinPoints(AliITSPConeData &p,Int_t i1,Int_t i2,Double_t r)
	{return p.ZAt(i2)+(p.ZAt(i1)-p.ZAt(i2))*(r-p.Rmin(i2))/
	     (p.Rmin(i1)-p.Rmin(i2));}
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p.Rmax(i1),p.ZAt(i1) and p.Rmax(i2),p.ZAt(i2)
    Double_t Zfrom2MaxPoints(AliITSPConeData &p,Int_t i1,Int_t i2,Double_t r)
	{return p.ZAt(i2)+(p.ZAt(i1)-p.ZAt(i2))*(r-p.Rmax(i2))/
	     (p.Rmax(i1)-p.Rmax(i2));}
    // General SSD Outer Cone surface equation Rmax.
    Double_t RmaxFromZpCone(AliITSPConeData &p,Double_t z,Double_t th=0.0)
	{return -fTantc*(z-p.ZAt(4))+p.Rmax(4)+th/fCostc;}
    // General SSD Inner Cone surface equation Rmin.
    Double_t RminFromZpCone(AliITSPConeData &p,Double_t z,Double_t th=0.0)
	{return -fTantc*(z-p.ZAt(3))+p.Rmin(3)+th/fCostc;}
    // Outer SSD cone surface
    Double_t RmaxFromZSSDcone(Double_t z,Double_t th=0.0) 
	{return RmaxFromZpCone(fA,z,th);}
    // Inner SSD cone surface
    Double_t RminFromZSSDcone(Double_t z,Double_t th=0.0) 
	{return RminFromZpCone(fA,z,th);}
    // General SSD Outer cone Surface equation for z.
    Double_t ZFromRmaxpCone(AliITSPConeData &p,Double_t r,Double_t th=0.0)
	{return p.ZAt(4)+(p.Rmax(4)+th/fCostc-r)/fTantc;}
    // General SSD Inner cone Surface equation for z.
    Double_t ZFromRminpCone(AliITSPConeData &p,Double_t r,Double_t th=0.0)
	{return p.ZAt(3)+(p.Rmin(3)+th/fCostc-r)/fTantc;}
    Double_t ZFromRmaxSSDcone(Double_t r,Double_t th=0.0)
	{return ZFromRmaxpCone(fA,r,th);}
    Double_t ZFromRminSSDcone(Double_t r,Double_t th=0.0)
	{return ZFromRminpCone(fA,r,th);}
    // Given a initial point z0,r0, the initial angle theta0, and the radius
    // of curvature, returns the point z1, r1 at the angle theta1. Theta
    // measured from the r axis in the clock wise direction [degrees].
    void RadiusOfCurvature(Double_t rc,Double_t theta0,Double_t z0,
			   Double_t r0,Double_t theta1,Double_t &z1,
			   Double_t &r1)
	{z1=rc*(Sind(theta1)-Sind(theta0))+z0;
	 r1=rc*(Cosd(theta1)-Cosd(theta0))+r0; return;}
 private:
    Double_t fThickness; //mm, Thickness of Rohacell+carbon fiber
    Double_t fCthick; //mm, Carbon finber thickness
    Double_t fRcurv; // mm, Radius of curvature.
    Double_t fTc; // angle of SSD cone [degrees].
    Double_t fSintc,fCostc,fTantc;
    Double_t fZ0,fZouterMilled,fZcylinder,fZposts;
    Double_t fRoutMax,fRoutHole,fRoutMin;
    Double_t fRholeMax,fRholeMin;
    Double_t fRpostMin,fdRpost,fZpostMax,fPhi0Post;
    Double_t fRinMax,fRinCylinder,fRinHole,fRinMin,fdZin;
    // Screws mounting SSD-SDD thermal/machanical cylinder
    Double_t fPhi0Screws,fRcylinderScrews,fDscrewHead;
    Double_t fDscrewShaft,fThScrewHeadHole;
    // The SDD mounting bracket SSD part of it.
    Int_t fNssdSupports;
    Double_t fPhi0SDDsupports;
    Double_t fRsddSupportPlate,fThSDDsupportPlate,fWsddSupportPlate;
    Int_t fNspoaks,fNinScrews,fNposts,fNmounts;
    Int_t fSSDcf; // SSD support cone Carbon Fiber materal number.
    Int_t fSSDfs; // SSD support cone inserto stesalite 4411w.
    Int_t fSSDfo; // SSD support cone foam, Rohacell 50A.
    Int_t fSSDsw; // SSD support cone screw material,Stainless steal
    Int_t fNcD; // number of screw ends (copy number) Volume D
    Int_t fNcE; // number of pin end (copy number) Volume E
    AliITSPConeData fA;//Poly-cone Volume A. Top part of SSD cone Carbon Fiber
    AliITSPConeData fB; // Poly-cone Volume B. Stesalite inside volume A.
    AliITSPConeData fC; // Poly-cone Volume C. Foam inside volume A.
    AliITSTubeData  fD; // Bolt holes willed with bolt
    AliITSTubeData  fE; // Pin holes willed with pin.
    AliITSPConeData fF;//Poly-cone Volume F. Foam in spoak reagion, inside A.
    // SSD support cone Spoaks.
    AliITSPConeData fG; // Poly-cone Volume G.
    AliITSPConeData fH; // Poly-cone Volume H.
    AliITSPConeData fI; // Poly-cone Volume I.
    AliITSPConeData fJ; // Poly-cone Volume J.
    AliITSPConeData fK; // Poly-cone Volume K.
    AliITSPConeData fL; // Poly-cone Volume L.
    AliITSPConeData fM; // Poly-cone Volume M.
    AliITSPConeData fN; // Poly-cone Volume N.
    AliITSPConeData fO; // Poly-cone Volume O.
    AliITSPConeData fP; // Poly-cone Volume P.
    AliITSPConeData fQ; // Poly-cone Volume Q.
    AliITSTubeData  fR; // Bolt holes willed with bolt
    AliITSTubeData  fS; // Pin holes willed with pin.
    AliITSPConeData fT; // Poly-cone Volume T.
    AliITSPConeData fU; // Poly-cone Volume U.

    ClassDef(AliITSGeometrySSDCone,1) // ITS SSD support cone geometry
	                              // version 11
};
// Input and Output functions for standard C++ input/output/
ostream &operator<<(ostream &os,AliITSGeometrySSDCone &source);
istream &operator>>(istream &is,AliITSGeometrySSDCone &source);

#endif

