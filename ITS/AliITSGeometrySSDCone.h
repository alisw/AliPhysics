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
 
class  AliITSGeometrySSDCone : public AliITSBaseGeometry {
 public:
    AliITSGeometrySSDCone();
    AliITSGeometrySSDCone(TVector3 *&tran,const char moth[3],Int_t mat0);
    virtual ~AliITSGeometrySSDCone();
    void CreateG3Geometry(const char moth[3],TVector3 &trans);
    void CreateG3Materials();
    void BuildDisplayGeometry();
 private:
    Double_t th; //mm, Thickness of Rohacell+carbon fiber
    Double_t ct; //mm, Carbon finber thickness
    Double_t r; // mm, Radius of curvature.
    Double_t tc; // angle of SSD cone [degrees].
    Double_t sintc,costc,tantc;
    Double_t z0,zcylinder,zpost;
    Double_t Routmax,RoutHole,Routmin;
    Double_t Rholemax,Rholemin;
    Double_t RPostmin,dRPost,zpostmax,phi0post;
    Double_t Rinmax,Rincylinder,RinHole,Rinmin,dzin;
    Int_t nspoaks,ninscrews,npost,nmounts;
    Int_t SSDcf; // SSD support cone Carbon Fiber materal number.
    Int_t SSDfs; // SSD support cone inserto stesalite 4411w.
    Int_t SSDfo; // SSD support cone foam, Rohacell 50A.
    Int_t SSDsw; // SSD support cone screw material,Stainless steal
    Int_t ncse; // number of screw ends (copy number)
    Int_t ncpe; // number of pin end (copy number)
    Int_t ncst; // number of screw tops (copy number)
    Double_t dphia,phi0a;
    Int_t nza;
    Double_t za[7];
    Double_t rmina[7];
    Double_t rmaxa[7];
    Double_t dphib,phi0b;
    Int_t nzb;
    Double_t zb[6];
    Double_t rminb[6];
    Double_t rmaxb[6];
    Double_t dphic,phi0c;
    Int_t nzc;
    Double_t zc[4];
    Double_t rminc[4];
    Double_t rmaxc[4];
    Double_t dphid,phi0d;
    Int_t nzd;
    Double_t zd[4];
    Double_t rmind[4];
    Double_t rmaxd[4];
    Double_t dze;
    Double_t rmine;
    Double_t rmaxe;
    Double_t dze2;
    Double_t rmine2;
    Double_t rmaxe2;
    Double_t dphif,phi0f;
    Int_t nzf;
    Double_t zf[4];
    Double_t rminf[4];
    Double_t rmaxf[4];
    Double_t dphig,phi0g;
    Int_t nzg;
    Double_t zg[4];
    Double_t rming[4];
    Double_t rmaxg[4];
    Double_t dphih,phi0h;
    Int_t nzh;
    Double_t zh[4];
    Double_t rminh[4];
    Double_t rmaxh[4];
    Double_t dphii,phi0i;
    Int_t nzi;
    Double_t zi[8];
    Double_t rmini[8];
    Double_t rmaxi[8];
    Double_t dphij,phi0j;
    Int_t nzj;
    Double_t zj[4];
    Double_t rminj[4];
    Double_t rmaxj[4];
    Double_t dphik,phi0k;
    Int_t nzk;
    Double_t zk[7];
    Double_t rmink[7];
    Double_t rmaxk[7];
    Double_t dphil,phi0l;
    Int_t nzl;
    Double_t zl[4];
    Double_t rminl[4];
    Double_t rmaxl[4];
    Double_t dphim,phi0m;
    Int_t nzm;
    Double_t zm[4];
    Double_t rminm[4];
    Double_t rmaxm[4];
    Double_t dphin,phi0n;
    Int_t nzn;
    Double_t zn[4];
    Double_t rminn[4];
    Double_t rmaxn[4];
    Double_t dphio,phi0o;
    Int_t nzo;
    Double_t zo[3];
    Double_t rmino[3];
    Double_t rmaxo[3];
    Double_t dphip,phi0p;
    Int_t nzp;
    Double_t zp[3];
    Double_t rminp[3];
    Double_t rmaxp[3];

    ClassDef(AliITSGeometrySSDCone,1)// ITS SSD support cone geometry version 1
};
 
#endif
