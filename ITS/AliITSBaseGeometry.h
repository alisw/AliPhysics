#ifndef ALIITSMIXTURE_H
#define ALIITSMIXTURE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */
#include <TGeoMaterial.h>

class AliITSMixture : public TGeoMixture{
 public:
    AliITSMixture(){};
    AliITSMixture(const char *name,Int_t N,Double_t *w,TObjArray *m,
		  Double_t rho=-1.,Double_t radlen=0.,Double_t intleng=0.);
    virtual ~AliITSMixture(){};
 private:
    ClassDef(AliITSMixture,1) // Extension of TGeoMixture class
}
;
#endif

//#ifndef ALIITSARB8_H
//#define ALIITSARB8_H
///* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// * See cxx source for full Copyright notice                               */
//
///*
//  $Id$
// */
//
//#include <TGeoArb8.h>
//class AliITSArb8 : public TGeoArb8{
//    AliITSArb8(){};
//    virtual  ~AliITSArb8(){};
//    Double_t  GetVertice(Int_t i,Int_t ixy){return (GetVertices())[i][ixy];}const// Returns value of fXY[i][ixy]
//    Double_t  GetVerticeX(Int_t i){return GetVertice(i,0);}const// Returns value of fXY[i][0]
//    Double_t  GetVerticeY(Int_t i){return GetVertice(i,1);}const// Returns value of fXY[i][1]
//    Double_t& Dz(){return fDz;}// Returns address of fDz
//    Double_t& VerticeX(Int_t i){return fXY[i][0];}// Returns address of fXY[i][0]
//    Double_t& VerticeY(Int_t i){return fXY[i][1];}// Returns address of fXY[i][1]
// private:
//    ClassDef(AliITSArb8,1) // Extension of TGeoArb8 class
//}
//;
//#endif

//#ifndef ALIITSBBOX_H
//#define ALIITSBBOX_H
///* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// * See cxx source for full Copyright notice                               */
//
///*
//  $Id$
// */
//
//#include <TGeoBBox.h>
//class AliITSBBox : public TGeoBBox{
//    AliITSBBox(){};
//    virtual  ~AliITSBBox(){};
//    Double_t  GetOrigin(Int_t ixyz){return (GetOrigin())[ixyz];}const// Returns value of fOrigin[ixyz]
//    Double_t  GetX0(){return GetOrigin(0);}const// Returns value of fOrigin[0]
//    Double_t  GetY0(){return GetOrigin(1);}const// Returns value of fOrigin[1]
//   Double_t  GetY0(){return GetOrigin(2);}const// Returns value of fOrigin[2]
//    Double_t& Dx(){return fDX;}// Returns address of fDx
//    Double_t& Dy(){return fDY;}// Returns address of fDy
//    Double_t& Dz(){return fDZ;}// Returns address of fDz
//    Double_t& X0(){return fOrigin[0];}// Returns address of fOrigin[0]
//    Double_t& Y0(){return fOrigin[1];}// Returns address of fOrigin[1]
//    Double_t& Z0(){return fOrigin[2];}// Returns address of fOrigin[2]
// private:
//    ClassDef(AliITSBBox,1) // Extension of TGeoBBox class
//}
//;
//#endif

//#ifndef ALIITSCONE_H
//#define ALIITSCONE_H
///* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// * See cxx source for full Copyright notice                               */
//
///*
//  $Id$
// */
//
//#include <TGeoCone.h>
//class AliITSCone : public TGeoCone{
//    AliITSCone(){};
//    virtual  ~AliITSCone(){};
//    Double_t& Dz(){return fDz;}// Returns address of fDz
//    Double_t& Rmin1(){return fRmin1;}// Returns address of fRmin1
//    Double_t& Rmax1(){return fRmax1;}// Returns address of fRmax1
//    Double_t& Rmin2(){return fRmin2;}// Returns address of fRmin2
//    Double_t& Rmax2(){return fRmax2;}// Returns address of fRmax2
// private:
//    ClassDef(AliITSCone,1) // Extension of TGeoCone class
//}
//;
//#endif

//#ifndef ALIITSCONESEG_H
//#define ALIITSCONESEG_H
///* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// * See cxx source for full Copyright notice                               */
//
///*
//  $Id$
// */
//
//#include <TGeoCone.h>
//class AliITSConeSeg : public TGeoConeSeg{
//    AliITSConeSeg(){};
//    virtual  ~AliITSConeSeg(){};
//    Double_t& Phi1(){return fPhi1;}// Returns address of fPhi1
//    Double_t& Phi2(){return fPhi2;}// Returns address of fPhi2
// private:
//    ClassDef(AliITSConeSeg,1) // Extension of TGeoConeSeg class
//}
//;
//#endif

//#ifndef ALIITSCTUB_H
//#define ALIITSCTUB_H
///* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// * See cxx source for full Copyright notice                               */
//
///*
//  $Id$
// */
//
//#include <TGeoTube.h>
//
//class AliITSCtub : public TGeoCtub{
//    AliITSCtub(){};
//    virtual  ~AliITSCtub(){};
//    Double_t  GetNlow(Int_t ixyz){return (GetNlow())[ixyz]};const // Returns value GetNlow[ixyz]
//    Double_t  GetNhigh(Int_t ixyz){return (GetNhigh())[ixyz]};const // Returns value GetNhigh[ixyz]
//    Double_t& NlowX(){return fNlow[0];}// Returns address of fNlow[0]
//    Double_t& NlowY(){return fNlow[1];}// Returns address of fNlow[1]
//    Double_t& NlowZ(){return fNlow[2];}// Returns address of fNlow[2]
//    Double_t& NhighX(){return fNhigh[0];}// Returns address of fNhigh[0]
//    Double_t& NhighY(){return fNhigh[1];}// Returns address of fNhigh[1]
//    Double_t& NhighZ(){return fNhigh[2];}// Returns address of fNhigh[2]
// private:
//    ClassDef(AliITSCtub,1) // Extension of TGeoCtub class
//}
//;
//#endif

//#ifndef ALIITSGTRA_H
//#define ALIITSGTRA_H
///* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// * See cxx source for full Copyright notice                               */
//
///*
//  $Id$
// */
//
//#include <TGeoArb8.h>
//class AliITSGtra : public TGeoGtra{
//    AliITSGtra(){};
//    virtual  ~AliITSGtra(){};
//    Double_t& TwistAngle(){return fTwistAngle;}// Returns address of fTwistAngle
// private:
//    ClassDef(AliITSGtra,1) // Extension of TGeoGtra class
//}
//;
//#endif

// TGeoPara : TGeoBox

//#ifndef ALIITSPCON_H
//#define ALIITSPCON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */
/*
#include <TGeoPcon.h>
class AliITSpCon : public TGeoPcon{
 public:
    AliITSpCon(){};
    AliITSpCon(const char* name,Double_t phi,Double_t dphi,Int_t nz){
	TGeoPcon(name,phi,dphi,nz);};
    virtual  ~AliITSpCon(){;};
    // Returns value of fRmin[i]
    Double_t  GetRmin(Int_t i)const{return (((TGeoPcon*)this)->GetRmin())[i];}
    // Returns value of fRmax[i]
    Double_t  GetRmax(Int_t i)const{return (((TGeoPcon*)this)->GetRmax())[i];}
    // Returns value of fZ[i]
    Double_t  GetZ(Int_t i)const{return (((TGeoPcon*)this)->GetZ())[i];}
 private:
    ClassDef(AliITSpCon,1) // Extension of TGeoPcon class
}
;
#endif
*/
// TGeoPgon : TGeoPcon
// TGeoSphere : TGeoBBox
// TGeoTrap : TGeoArb8
// TGeoTrd1 : TGeoBBox
// TGeoTrd2 : TGeoBBox


//#ifndef ALIITSTUBE_H
//#define ALIITSTUBE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */
/*
#include <TGeoTube.h>
class AliITSTube : public TGeoTube{
 public:
    AliITSTube(){};
    AliITSTube(const char *name,double_t rmin,Double_t rmax,Double_t dz){
	TGeoTube(name,rmin,rmax,dz);};
    virtual  ~AliITSTube(){};
    Double_t& Rmin(){return fRmin;}// Returns address of fRmin
    Double_t& Rmax(){return fRmax;}// Returns address of fRmax
    Double_t& Dz(){return fDZ;}// Returns address of fDz
 private:
    ClassDef(AliITSTube,1) // Extension of TGeoTube class
}
;
#endif
*/
// TGeoTubeSeg : TGeoTube

//#ifndef ALIITSELTU_H
//#define ALIITSELTU_H
///* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// * See cxx source for full Copyright notice                               */
//
///*
//  $Id$
// */
//#include <TGeoTube.h>
//#include <TGeoEltu.h>
//class AliITSEltu : public TGeoEltu{
//    AliITSEltu(){};
//    virtual  ~AliITSEltu(){};
//    Double_t& A(){return TGeoTube::Rmin();}// Returns address of A
//    Double_t& B(){return TGeoTube::Rmax();}// Returns address of B
//    Double_t& Dz(){return TGeoTube::Dz();}// Returns address of Dz
// private:
//    ClassDef(AliITSEltu,1) // Extension of TGeoEltu class
//}
//;
//#endif
