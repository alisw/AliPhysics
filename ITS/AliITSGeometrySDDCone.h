#ifndef ALIITSGEOMETRYSDDCONE_H
#define ALIITSGEOMETRYSDDCONE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */

/*
  ITS SDD Cone Geometry. Version 11
*/
#include "AliITSBaseGeometry.h"
class TVector3;
class AliITS;
 
class  AliITSGeometrySDDCone : public AliITSBaseGeometry {
 public:
    AliITSGeometrySDDCone();
    AliITSGeometrySDDCone(AliITS *its,TVector3 *&tran,const char moth[3],Int_t mat0);
    virtual ~AliITSGeometrySDDCone(){};
    void CreateG3Geometry(const char moth[3],TVector3 &trans);
    void CreateG3Materials();
    void BuildDisplayGeometry();
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
    Int_t fSDDcf; // SSD support cone Carbon Fiber materal number.
    Int_t fSDDfs; // SSD support cone inserto stesalite 4411w.
    Int_t fSDDfo; // SSD support cone foam, Rohacell 50A.
    Int_t fSDDsw; // SSD support cone screw material,Stainless steal
    Int_t fNcse; // number of screw ends (copy number)
    Int_t fNcpe; // number of pin end (copy number)
    Int_t fNcst; // number of screw tops (copy number)
    AliITSPConeData fA;//Poly-cone Volume A.
    AliITSPConeData fB; // Poly-cone Volume B.
    AliITSPConeData fC; // Poly-cone Volume C.
    Double_t fDdZ,fDRmin,fDRmax; // Tube Volume D.
    Double_t fEdZ,fERmin,fERmax; // Tube Volume E.
    AliITSPConeData fF;//Poly-cone Volume F.
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
    Double_t fRdZ,fRRmin,fRRmax; // Tube Volume R
    Double_t fSdZ,fSRmin,fSRmax; // Tube Volume S
    AliITSPConeData fT; // Poly-cone Volume T.
    AliITSPConeData fU; // Poly-cone Volume U.
    

    ClassDef(AliITSGeometrySDDCone,1) // ITS SDD support cone geometry
	                              // version 1
};
 
#endif

