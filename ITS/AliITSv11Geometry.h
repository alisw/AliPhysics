#ifndef ALIITSV11GEOMETRY_H
#define ALIITSV11GEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
 */
#include <TObject.h>
class TGeoArb8;
class TGeoPcon;
class TGeoTube;
class TGeoTubeSeg;
class TGeoConeSeg;
class TGeoBBox;

class AliITSv11Geometry : public TObject {
  public:
    AliITSv11Geometry(){fDebug=kTRUE;};
    AliITSv11Geometry(Bool_t debug){fDebug=debug;};
    virtual ~AliITSv11Geometry(){};
    //
    void SetDebug(){fDebug=kTRUE;}
    void SetNoDebug(){fDebug=kFALSE;}
    static Double_t RmaxFrom2Points(TGeoPcon *p,Int_t i1,Int_t i2,Double_t z);
    static Double_t RminFrom2Points(TGeoPcon *p,Int_t i1,Int_t i2,Double_t z);
    static Double_t RFrom2Points(Double_t *p,Double_t *Z,Int_t i1,Int_t i2,
                                 Double_t z);
    static Double_t Zfrom2MinPoints(TGeoPcon *p,Int_t i1,Int_t i2,Double_t r);
    static Double_t Zfrom2MaxPoints(TGeoPcon *p,Int_t i1,Int_t i2,Double_t r);
    static Double_t Zfrom2Points(Double_t *Z,Double_t *p,Int_t i1,Int_t i2,
                                 Double_t r);
    static Double_t RmaxFromZpCone(TGeoPcon *p,int ip,Double_t tc,Double_t z,
                                   Double_t th=0.0);
    static Double_t RmaxFromZpCone(TGeoPcon *p,Double_t tc,Double_t z,
                                   Double_t th=0.0){
        return RmaxFromZpCone(p,4,tc,z,th);};
    static Double_t RFromZpCone(Double_t *Rmax,Double_t *Z,int ip,Double_t tc,
                                Double_t z,Double_t th=0.0);
    static Double_t RmaxFromZpCone(Double_t *Rmax,Double_t *Z,Double_t tc,
                                   Double_t z,Double_t th=0.0){
        return RFromZpCone(Rmax,Z,4,tc,z,th);};
    static Double_t RminFromZpCone(TGeoPcon *p,Int_t ip,Double_t tc,Double_t z,
                                   Double_t th=0.0);
    static Double_t RminFromZpCone(TGeoPcon *p,Double_t tc,Double_t z,
                                   Double_t th=0.0){
        return RminFromZpCone(p,3,tc,z,th);};
    static Double_t RminFromZpCone(Double_t *Rmin,Double_t *Z,Double_t tc,
                                   Double_t z,Double_t th=0.0){
        return RFromZpCone(Rmin,Z,3,tc,z,th);};
    static Double_t ZFromRmaxpCone(TGeoPcon *p,int ip,Double_t tc,Double_t r,
                                   Double_t th=0.0);
    static Double_t ZFromRmaxpCone(TGeoPcon *p,Double_t tc,Double_t r,
                                   Double_t th=0.0)
        {return ZFromRmaxpCone(p,4,tc,r,th);};
    static Double_t ZFromRmaxpCone(Double_t *GetRmax,Double_t *GetZ,Int_t ip,
                                   Double_t tc,Double_t r,Double_t th=0.0);
    static Double_t ZFromRmaxpCone(Double_t *GetRmax,Double_t *GetZ,
                                   Double_t tc,Double_t r,Double_t th=0.0){
        return ZFromRmaxpCone(GetRmax,GetZ,4,tc,r,th);};
    static Double_t ZFromRminpCone(TGeoPcon *p,int ip,Double_t tc,Double_t r,
                                   Double_t th=0.0);
    static Double_t ZFromRminpCone(TGeoPcon *p,Double_t tc,Double_t r,
                                   Double_t th=0.0)
        {return ZFromRminpCone(p,3,tc,r,th);};
    static void InsidePoint(TGeoPcon *p,Int_t i1,Int_t i2,Int_t i3,
                            Double_t Cthick,TGeoPcon *q,Int_t j1,Bool_t max);
    static void InsidePoint(Double_t x0,Double_t y0,Double_t x1,Double_t y1,
                            Double_t x2,Double_t y2,Double_t c,
                            Double_t &x,Double_t &y);
    static void RadiusOfCurvature(Double_t rc,Double_t theta0,Double_t z0,
                                  Double_t r0,Double_t theta1,Double_t &z1,
                                  Double_t &r1);
    void printArb8(TGeoArb8 *A);
    void printPcon(TGeoPcon *A);
    void printTube(TGeoTube *A);
    void printTubeSeg(TGeoTubeSeg *A);
    void printConeSeg(TGeoConeSeg *A);
    void printBBox(TGeoBBox *A);
    Bool_t GetDebug(){return fDebug;}

  private:
    Bool_t fDebug; //! Debug flag
    ClassDef(AliITSv11Geometry,1) // Base class for ITS v11 geometry
};


// Units, Convert from k?? to cm,degree,GeV,seconds,
const Double_t kmm = 0.10; // Convert mm to TGeom's cm.
const Double_t kcm = 1.00; // Convert cv to TGeom's cm.
const Double_t kDegree = 1.0; // Convert degrees to TGeom's degrees
const Double_t kRadian = TMath::DegToRad(); // conver to Radians

#endif
