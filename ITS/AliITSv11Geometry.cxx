/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */
#include <stdio.h>
#include <stdlib.h>
// General Root includes
#include <Riostream.h>
#include <TMath.h>
#include <float.h>
#include <TFile.h>    // only required for Tracking function?
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TObjString.h>
// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include "AliITSBaseGeometry.h"
#include "AliITSv11Geometry.h"

ClassImp(AliITSv11Geometry)
//______________________________________________________________________
Double_t AliITSv11Geometry::RmaxFrom2Points(TGeoPcon *p,Int_t i1,Int_t i2,Double_t z){
    // functions Require at parts of Volume A to be already defined.
    // Retruns the value of Rmax corresponding to point z alone the line
    // defined by the two points p.Rmax(i1),p-GetZ(i1) and p->GetRmax(i2),
    // p->GetZ(i2).
    Double_t d0,d1,d2,r;

    d0 = p->GetRmax(i1)-p->GetRmax(i2);// cout <<"L263: d0="<<d0<<endl;
    d1 = z-p->GetZ(i2);// cout <<"L264: d1="<<d1<<endl;
    d2 = p->GetZ(i1)-p->GetZ(i2);// cout <<"L265: d2="<<d2<<endl;
    r  = p->GetRmax(i2) + d1*d0/d2;// cout <<"L266: r="<<r<<endl;
    return r;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RminFrom2Points(TGeoPcon *p,Int_t i1,Int_t i2,Double_t z){
    // Retruns the value of Rmin corresponding to point z alone the line
    // defined by the two points p->GetRmin(i1),p->GetZ(i1) and 
    // p->GetRmin(i2),  p->GetZ(i2).

    return p->GetRmin(i2)+(p->GetRmin(i1)-p->GetRmin(i2))*(z-p->GetZ(i2))/
     (p->GetZ(i1)-p->GetZ(i2));
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RFrom2Points(Double_t *p,Double_t *Z,Int_t i1,
                                 Int_t i2,Double_t z){
    // Retruns the value of Rmin corresponding to point z alone the line
    // defined by the two points p->GetRmin(i1),p->GetZ(i1) and 
    // p->GetRmin(i2), p->GetZ(i2).

    return p[i2]+(p[i1]-p[i2])*(z-Z[i2])/(Z[i1]-Z[i2]);
}
//______________________________________________________________________
Double_t AliITSv11Geometry::Zfrom2MinPoints(TGeoPcon *p,Int_t i1,Int_t i2,Double_t r){
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p->GetRmin(i1),p->GetZ(i1) and 
    // p->GetRmin(i2),p->GetZ(i2)

    return p->GetZ(i2)+(p->GetZ(i1)-p->GetZ(i2))*(r-p->GetRmin(i2))/
     (p->GetRmin(i1)-p->GetRmin(i2));
}
//______________________________________________________________________
Double_t AliITSv11Geometry::Zfrom2MaxPoints(TGeoPcon *p,Int_t i1,Int_t i2,Double_t r){
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p->GetRmax(i1),p->GetZ(i1) and 
    // p->GetRmax(i2),p->GetZ(i2)

    return p->GetZ(i2)+(p->GetZ(i1)-p->GetZ(i2))*(r-p->GetRmax(i2))/
     (p->GetRmax(i1)-p->GetRmax(i2));
}
//______________________________________________________________________
Double_t AliITSv11Geometry::Zfrom2Points(Double_t *Z,Double_t *p,Int_t i1,
                                 Int_t i2,Double_t r){
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p->GetRmax(i1),p->GetZ(i1) and 
    // p->GetRmax(i2),p->GetZ(i2)

    return Z[i2]+(Z[i1]-Z[i2])*(r-p[i2])/(p[i1]-p[i2]);
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RmaxFromZpCone(TGeoPcon *p,int ip,Double_t tc,Double_t z,
                                   Double_t th){
    // General SSD Outer Cone surface equation Rmax.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-p->GetZ(ip))+p->GetRmax(ip)+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RFromZpCone(Double_t *GetRmax,Double_t *GetZ,int ip,
                                   Double_t tc,Double_t z,Double_t th){
    // General SSD Outer Cone surface equation Rmax.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-GetZ[ip])+GetRmax[ip]+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RminFromZpCone(TGeoPcon *p,Int_t ip,Double_t tc,Double_t z,
                                   Double_t th){
    // General SSD Inner Cone surface equation Rmin.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-p->GetZ(ip))+p->GetRmin(ip)+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::ZFromRmaxpCone(TGeoPcon *p,int ip,Double_t tc,Double_t r,
                                   Double_t th){
    // General SSD Outer cone Surface equation for z.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return p->GetZ(ip)+(p->GetRmax(ip)+th/costc-r)/tantc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::ZFromRmaxpCone(Double_t *GetRmax,Double_t *GetZ,int ip,
                                   Double_t tc,Double_t r,Double_t th){
    // General SSD Outer cone Surface equation for z.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return GetZ[ip]+(GetRmax[ip]+th/costc-r)/tantc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::ZFromRminpCone(TGeoPcon *p,int ip,Double_t tc,Double_t r,
                                   Double_t th){
    // General SSD Inner cone Surface equation for z.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return p->GetZ(ip)+(p->GetRmin(ip)+th/costc-r)/tantc;
}
//______________________________________________________________________
void AliITSv11Geometry::RadiusOfCurvature(Double_t rc,Double_t theta0,Double_t z0,
                 Double_t r0,Double_t theta1,Double_t &z1,
                 Double_t &r1){
    // Given a initial point z0,r0, the initial angle theta0, and the radius
    // of curvature, returns the point z1, r1 at the angle theta1. Theta
    // measured from the r axis in the clock wise direction [degrees].
    Double_t sin0 = TMath::Sin(theta0*TMath::DegToRad());
    Double_t cos0 = TMath::Cos(theta0*TMath::DegToRad());
    Double_t sin1 = TMath::Sin(theta1*TMath::DegToRad());
    Double_t cos1 = TMath::Cos(theta1*TMath::DegToRad());

    z1 = rc*(sin1-sin0)+z0;
    r1 = rc*(cos1-cos0)+r0;
    return;
}
//______________________________________________________________________
void AliITSv11Geometry::InsidePoint(TGeoPcon *p,Int_t i1,Int_t i2,Int_t i3,
                            Double_t c,TGeoPcon *q,Int_t j1,Bool_t max){
    // Given two lines defined by the points i1, i2,i3 in the TGeoPcon 
    // class p that intersect at point p->GetZ(i2) return the point z,r 
    // that is Cthick away in the TGeoPcon class q. If points i1=i2
    // and max == kTRUE, then p->GetRmin(i1) and p->GetRmax(i2) are used.
    // if points i2=i3 and max=kTRUE then points p->GetRmax(i2) and
    // p->GetRmin(i3) are used. If i2=i3 and max=kFALSE, then p->GetRmin(i2)
    // and p->GetRmax(i3) are used.
    // Inputs:
    //    TGeoPcon  *p  Class where points i1, i2, and i3 are taken from
    //    Int_t     i1  First point in class p
    //    Int_t     i2  Second point in class p
    //    Int_t     i3  Third point in class p
    //    Double_t  c   Distance inside the outer surface/inner suface
    //                  that the point j1 is to be computed for.
    //    TGeoPcon  *q  Pointer to class for results to be put into.
    //    Int_t     j1  Point in class q where data is to be stored.
    //    Bool_t    max if kTRUE, then a Rmax value is computed,
    //                  else a Rmin valule is computed.
    // Output:
    //    TGeoPcon  *q  Pointer to class for results to be put into.
    // Return:
    //    none.
    Double_t x0,y0,x1,y1,x2,y2,x,y;

    if(max){
        c = -c; //cout <<"L394 c="<<c<<endl;
        y0 = p->GetRmax(i1);
        if(i1==i2) y0 = p->GetRmin(i1); //cout <<"L396 y0="<<y0<<endl;
        y1 = p->GetRmax(i2);  //cout <<"L397 y1="<<y1<<endl;
        y2 = p->GetRmax(i3); //cout <<"L398 y2="<<y2<<endl;
        if(i2==i3) y2 = p->GetRmin(i3); //cout <<"L399 y2="<<y2<<endl;
    }else{ // min
        y0 = p->GetRmin(i1); //cout <<"L401 y0="<<y0<<endl;
        y1 = p->GetRmin(i2); //cout <<"L402 y1="<<y1<<endl;
        y2 = p->GetRmin(i3);
        if(i2==i3) y2 = p->GetRmax(i3); //cout <<"L404 y2="<<y2<<endl;
    } // end if
    x0 = p->GetZ(i1); //cout <<"L406 x0="<<x0<<endl;
    x1 = p->GetZ(i2); //cout <<"L407 x1="<<x1<<endl;
    x2 = p->GetZ(i3); //cout <<"L408 x2="<<x2<<endl;
    //
    InsidePoint(x0,y0,x1,y1,x2,y2,c,x,y);
    q->Z(j1) = x;
    if(max) q->Rmax(j1) = y;
    else    q->Rmin(j1) = y;
    return;
}
//----------------------------------------------------------------------
void AliITSv11Geometry::InsidePoint(Double_t x0,Double_t y0,Double_t x1,Double_t y1,
                            Double_t x2,Double_t y2,Double_t c,
                            Double_t &x,Double_t &y){
    // Given two intersecting lines defined by the points (x0,y0), (x1,y1) and
    // (x1,y1), (x1,y2) {intersecting at (x1,y1)} the point (x,y) a distance
    // c away is returned such that two lines a distance c away from the
    // lines defined above intersect at (x,y).
    // Inputs:
    //    Double_t  x0 X point on the first intersecting sets of lines
    //    Double_t  y0 Y point on the first intersecting sets of lines
    //    Double_t  x1 X point on the first/second intersecting sets of lines
    //    Double_t  y1 Y point on the first/second intersecting sets of lines
    //    Double_t  x2 X point on the second intersecting sets of lines
    //    Double_t  y2 Y point on the second intersecting sets of lines
    //    Double_t  c  Distance the two sets of lines are from each other
    // Output:
    //    Double_t  x  X point for the intersecting sets of parellel lines
    //    Double_t  y  Y point for the intersecting sets of parellel lines
    // Return:
    //    none.
    Double_t dx01,dx12,dy01,dy12,R01,R12,m;
    dx01 = x0-x1; //cout <<"L410 dx01="<<dx01<<endl;
    dx12 = x1-x2; //cout <<"L411 dx12="<<dx12<<endl;
    dy01 = y0-y1; //cout <<"L412 dy01="<<dy01<<endl;
    dy12 = y1-y2; //cout <<"L413 dy12="<<dy12<<endl;
    R01  = TMath::Sqrt(dy01*dy01+dx01*dx01); //cout <<"L414 R01="<<R01<<endl;
    R12  = TMath::Sqrt(dy12*dy12+dx12*dx12); //cout <<"L415 R12="<<R12<<endl;
    m = dx12*dy01-dy12*dx01;
    if(m*m<DBL_EPSILON){ // m == n
        if(dy01==0.0){ // line are =
            x = x1+c; //cout <<"L419 x="<<x<<endl;
            y = y1; //cout <<"L420 y="<<y<<endl;
            return;
        }else if(dx01==0.0){
            x = x1;
            y = y1+c;
            return;
        }else{ // dx01!=0 and dy01 !=0.
            x = x1-0.5*c*R01/dy01; //cout <<"L434 x="<<x<<endl;
            y = y1+0.5*c*R01/dx01; //cout <<"L435 y="<<y<<endl;
        } // end if
        return;
    } //
    x = x1-c*(dx12*R01-dx01*R12)/m; //cout <<"L442 x="<<x<<endl;
    y = y1-c*(dy12*R01-dy01*R12)/m; //cout <<"L443 y="<<y<<endl;
    //cout <<"=============================================="<<endl;
    return;
}
//----------------------------------------------------------------------
void AliITSv11Geometry:: printArb8(TGeoArb8 *A){
    if(GetDebug()){
        cout << A->GetName() << ":";
        for(Int_t iii=0;iii<8;iii+=2){
            cout <<"("<<A->GetVertices()[iii]<<","
                 <<A->GetVertices()[iii+1]<<","<<-A->GetDz()<<")";
        } // end for iii
        for(Int_t iii=8;iii<16;iii+=2){
            cout <<"("<<A->GetVertices()[iii]<<","
                 <<A->GetVertices()[iii+1]<<","<<A->GetDz()<<")";
        } // end for iii
        cout << endl;
    } // end if
}
//----------------------------------------------------------------------
void AliITSv11Geometry:: printPcon(TGeoPcon *A){  
    if(GetDebug()) return;
    cout << A->GetName() << ": N=" << A->GetNz() << " Phi1=" << A->GetPhi1()
         << ", Dphi=" << A->GetDphi() << endl;
    cout << "i\t   Z   \t  Rmin \t  Rmax" << endl;
    for(Int_t iii=0;iii<A->GetNz();iii++){
        cout << iii << "\t" << A->GetZ(iii) << "\t" << A->GetRmin(iii)
             << "\t" << A->GetRmax(iii) << endl;
    } // end for iii
}
//----------------------------------------------------------------------
void AliITSv11Geometry::printTube(TGeoTube *A){
    if(GetDebug()) return;
    cout << A->GetName() <<": Rmin="<<A->GetRmin()
         <<" Rmax=" <<A->GetRmax()<<" Dz="<<A->GetDz()<<endl;
}
//----------------------------------------------------------------------
void AliITSv11Geometry::printTubeSeg(TGeoTubeSeg *A){
    if(GetDebug()) return;
    cout << A->GetName() <<": Phi1="<<A->GetPhi1()<<
        " Phi2="<<A->GetPhi2()<<" Rmin="<<A->GetRmin()
         <<" Rmax=" <<A->GetRmax()<<" Dz="<<A->GetDz()<<endl;
}
//----------------------------------------------------------------------
void AliITSv11Geometry::printConeSeg(TGeoConeSeg *A){
    if(GetDebug()) return;
    cout << A->GetName() <<": Phi1="<<A->GetPhi1()<<
        " Phi2="<<A->GetPhi2()<<" Rmin1="<<A->GetRmin1()
         <<" Rmax1=" <<A->GetRmax1()<<" Rmin2="<<A->GetRmin2()
         <<" Rmax2=" <<A->GetRmax2()<<" Dz="<<A->GetDz()<<endl;
}
//----------------------------------------------------------------------
void AliITSv11Geometry::printBBox(TGeoBBox *A){
    if(GetDebug()) return;
    cout << A->GetName() <<": Dx="<<A->GetDX()<<
        " Dy="<<A->GetDY()<<" Dz="<<A->GetDZ() <<endl;
}

