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

/*
 $Id$ 
*/


////////////////////////////////////////////////////////////////////////
//  This class is a base class for the ITS geometry version 11. It 
//  contains common/standard functions used in many places in defining 
//  the ITS geometry, version 11. Large posions of the ITS geometry, 
//  version 11, should be derived from this class so as to make maximum 
//  use of these common functions. This class also defines the proper 
//  conversion valuse such, to cm and degrees, such that the most usefull 
//  units, those used in the Engineering drawings, can be used.
////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <TMath.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TPolyMarker.h>
#include <TPolyLine.h>
#include "AliITSv11Geometry.h"

ClassImp(AliITSv11Geometry)
//______________________________________________________________________
Double_t AliITSv11Geometry::Yfrom2Points(Double_t x0,Double_t y0,
                                         Double_t x1,Double_t y1,
                                         Double_t x)const{
    // Given the two points (x0,y0) and (x1,y1) and the location x, returns
    // the value y corresponding to that point x on the line defined by the
    // two points.
    // Inputs:
    //    Double_t  x0  The first x value defining the line
    //    Double_t  y0  The first y value defining the line
    //    Double_t  x1  The second x value defining the line
    //    Double_t  y1  The second y value defining the line
    //    Double_t   x  The x value for which the y value is wanted.
    // Outputs:
    //    none.
    // Return:
    //    The value y corresponding to the point x on the line defined by
    //    the two points (x0,y0) and (x1,y1).

    if(x0==x1 && y0==y1) {
        printf("Error: AliITSv11Geometry::Yfrom2Ponts The two points are "
               "the same (%e,%e) and (%e,%e)",x0,y0,x1,y1);
        return 0.0;
    } // end if
    if(x0==x1){
        printf("Warning: AliITSv11Geometry::Yfrom2Points x0=%e == x1=%e. "
               "line vertical ""returning mean y",x0,x1);
        return 0.5*(y0+y1);
    }// end if x0==x1
    Double_t m = (y0-y1)/(x0-x1);
    return m*(x-x0)+y0;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::Xfrom2Points(Double_t x0,Double_t y0,
                                         Double_t x1,Double_t y1,
                                         Double_t y)const{
    // Given the two points (x0,y0) and (x1,y1) and the location y, returns
    // the value x corresponding to that point y on the line defined by the
    // two points.
    // Inputs:
    //    Double_t  x0  The first x value defining the line
    //    Double_t  y0  The first y value defining the line
    //    Double_t  x1  The second x value defining the line
    //    Double_t  y1  The second y value defining the line
    //    Double_t   y  The y value for which the x value is wanted.
    // Outputs:
    //    none.
    // Return:
    //    The value x corresponding to the point y on the line defined by
    //    the two points (x0,y0) and (x1,y1).

    if(x0==x1 && y0==y1) {
        printf("Error: AliITSv11Geometry::Yfrom2Ponts The two points are "
               "the same (%e,%e) and (%e,%e)",x0,y0,x1,y1);
        return 0.0;
    } // end if
    if(y0==y1){
        printf("Warrning: AliITSv11Geometry::Yfrom2Points y0=%e == y1=%e. "
               "line horizontal returning mean x",y0,y1);
        return 0.5*(x0+x1);
    }// end if y0==y1
    Double_t m = (x0-x1)/(y0-y1);
    return m*(y-y0)+x0;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RmaxFrom2Points(const TGeoPcon *p,Int_t i1,
                                            Int_t i2,Double_t z)const{
    // functions Require at parts of Volume A to be already defined.
    // Retruns the value of Rmax corresponding to point z alone the line
    // defined by the two points p.Rmax(i1),p-GetZ(i1) and p->GetRmax(i2),
    // p->GetZ(i2).
    // Inputs:
    //    TGeoPcon *p  The Polycone where the two points come from
    //    Int_t    i1  Point 1
    //    Int_t    i2  Point 2
    //    Double_t  z  The value of z for which Rmax is to be found
    // Outputs:
    //    none.
    // Return:
    //    Double_t Rmax the value corresponding to z
    Double_t d0,d1,d2,r;

    d0 = p->GetRmax(i1)-p->GetRmax(i2);// cout <<"L263: d0="<<d0<<endl;
    d1 = z-p->GetZ(i2);// cout <<"L264: d1="<<d1<<endl;
    d2 = p->GetZ(i1)-p->GetZ(i2);// cout <<"L265: d2="<<d2<<endl;
    r  = p->GetRmax(i2) + d1*d0/d2;// cout <<"L266: r="<<r<<endl;
    return r;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RminFrom2Points(const TGeoPcon *p,Int_t i1,
                                            Int_t i2,Double_t z)const{
    // Retruns the value of Rmin corresponding to point z alone the line
    // defined by the two points p->GetRmin(i1),p->GetZ(i1) and 
    // p->GetRmin(i2),  p->GetZ(i2).
    // Inputs:
    //    TGeoPcon *p  The Polycone where the two points come from
    //    Int_t    i1  Point 1
    //    Int_t    i2  Point 2
    //    Double_t  z  The value of z for which Rmax is to be found
    // Outputs:
    //    none.
    // Return:
    //    Double_t Rmax the value corresponding to z

    return p->GetRmin(i2)+(p->GetRmin(i1)-p->GetRmin(i2))*(z-p->GetZ(i2))/
     (p->GetZ(i1)-p->GetZ(i2));
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RFrom2Points(const Double_t *p,const Double_t *az,
                                         Int_t i1,Int_t i2,Double_t z)const{
    // Retruns the value of Rmin corresponding to point z alone the line
    // defined by the two points p->GetRmin(i1),p->GetZ(i1) and 
    // p->GetRmin(i2), p->GetZ(i2).
    // Inputs:
    //    Double_t az  Array of z values
    //    Double_t  r  Array of r values
    //    Int_t    i1  First Point in arrays
    //    Int_t    i2  Second Point in arrays
    //    Double_t z   Value z at which r is to be found
    // Outputs:
    //    none.
    // Return:
    //    The value r corresponding to z and the line defined by the two points

    return p[i2]+(p[i1]-p[i2])*(z-az[i2])/(az[i1]-az[i2]);
}
//______________________________________________________________________
Double_t AliITSv11Geometry::Zfrom2MinPoints(const TGeoPcon *p,Int_t i1,
                                            Int_t i2,Double_t r)const{
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p->GetRmin(i1),p->GetZ(i1) and 
    // p->GetRmin(i2),p->GetZ(i2)
    // Inputs:
    //    TGeoPcon *p  The Poly cone where the two points come from.
    //    Int_t    i1  First Point in arrays
    //    Int_t    i2  Second Point in arrays
    //    Double_t r   Value r min at which z is to be found
    // Outputs:
    //    none.
    // Return:
    //    The value z corresponding to r min and the line defined by 
    //    the two points

    return p->GetZ(i2)+(p->GetZ(i1)-p->GetZ(i2))*(r-p->GetRmin(i2))/
     (p->GetRmin(i1)-p->GetRmin(i2));
}
//______________________________________________________________________
Double_t AliITSv11Geometry::Zfrom2MaxPoints(const TGeoPcon *p,Int_t i1,
                                            Int_t i2,Double_t r)const{
    // Retruns the value of Z corresponding to point R alone the line
    // defined by the two points p->GetRmax(i1),p->GetZ(i1) and 
    // p->GetRmax(i2),p->GetZ(i2)
    // Inputs:
    //    TGeoPcon *p  The Poly cone where the two points come from.
    //    Int_t    i1  First Point in arrays
    //    Int_t    i2  Second Point in arrays
    //    Double_t r   Value r max at which z is to be found
    // Outputs:
    //    none.
    // Return:
    //    The value z corresponding to r max and the line defined by 
    //    the two points

    return p->GetZ(i2)+(p->GetZ(i1)-p->GetZ(i2))*(r-p->GetRmax(i2))/
     (p->GetRmax(i1)-p->GetRmax(i2));
}
//______________________________________________________________________
Double_t AliITSv11Geometry::Zfrom2Points(const Double_t *z,const Double_t *ar,
                                         Int_t i1,Int_t i2,Double_t r)const{
    // Retruns the value of z corresponding to point R alone the line
    // defined by the two points p->GetRmax(i1),p->GetZ(i1) and 
    // p->GetRmax(i2),p->GetZ(i2)
    // Inputs:
    //    Double_t  z  Array of z values
    //    Double_t ar  Array of r values
    //    Int_t    i1  First Point in arrays
    //    Int_t    i2  Second Point in arrays
    //    Double_t r   Value r at which z is to be found
    // Outputs:
    //    none.
    // Return:
    //    The value z corresponding to r and the line defined by the two points

    return z[i2]+(z[i1]-z[i2])*(r-ar[i2])/(ar[i1]-ar[i2]);
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RmaxFromZpCone(const TGeoPcon *p,int ip,
                                           Double_t tc,Double_t z,
                                           Double_t th)const{
    // General Outer Cone surface equation Rmax.
    // Intputs:
    //     TGeoPcon  *p   The poly cone where the initial point comes from
    //     Int_t     ip   The index in p to get the point location
    //     Double_t  tc   The angle of that part of the cone is at
    //     Double_t   z   The value of z to compute Rmax from
    //     Double_t  th   The perpendicular distance the parralell line is
    //                    from the point ip.
    // Outputs:
    //     none.
    // Return:
    //     The value Rmax correstponding to the line at angle th, offeset by
    //     th, and the point p->GetZ/Rmin[ip] at the location z.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-p->GetZ(ip))+p->GetRmax(ip)+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RFromZpCone(const Double_t *ar,
                                        const Double_t *az,int ip,
                                        Double_t tc,Double_t z,
                                        Double_t th)const{
    // General Cone surface equation R(z).
    // Intputs:
    //     Double_t  ar   The array of R values
    //     Double_t  az   The array of Z values
    //     Int_t     ip   The index in p to get the point location
    //     Double_t  tc   The angle of that part of the cone is at
    //     Double_t   z   The value of z to compute R from
    //     Double_t  th   The perpendicular distance the parralell line is
    //                    from the point ip.
    // Outputs:
    //     none.
    // Return:
    //     The value R correstponding to the line at angle th, offeset by
    //     th, and the point p->GetZ/Rmax[ip] at the locatin z.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-az[ip])+ar[ip]+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::RminFromZpCone(const TGeoPcon *p,Int_t ip,
                                           Double_t tc,Double_t z,
                                           Double_t th)const{
    // General Inner Cone surface equation Rmin.
    // Intputs:
    //     TGeoPcon  *p   The poly cone where the initial point comes from
    //     Int_t     ip   The index in p to get the point location
    //     Double_t  tc   The angle of that part of the cone is at
    //     Double_t   z   The value of z to compute Rmin from
    //     Double_t  th   The perpendicular distance the parralell line is
    //                    from the point ip.
    // Outputs:
    //     none.
    // Return:
    //     The value Rmin correstponding to the line at angle th, offeset by
    //     th, and the point p->GetZ/Rmin[ip] at the location z.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return -tantc*(z-p->GetZ(ip))+p->GetRmin(ip)+th/costc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::ZFromRmaxpCone(const TGeoPcon *p,int ip,
                                           Double_t tc,Double_t r,
                                           Double_t th)const{
    // General Outer cone Surface equation for z.
    // Intputs:
    //     TGeoPcon  *p   The poly cone where the initial point comes from
    //     Int_t     ip   The index in p to get the point location
    //     Double_t  tc   The angle of that part of the cone is at
    //     Double_t   r   The value of Rmax to compute z from
    //     Double_t  th   The perpendicular distance the parralell line is
    //                    from the point ip.
    // Outputs:
    //     none.
    // Return:
    //     The value Z correstponding to the line at angle th, offeset by
    //     th, and the point p->GetZ/Rmax[ip] at the location r.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return p->GetZ(ip)+(p->GetRmax(ip)+th/costc-r)/tantc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::ZFromRmaxpCone(const Double_t *ar,
                                           const Double_t *az,int ip,
                                           Double_t tc,Double_t r,
                                           Double_t th)const{
    // General Outer cone Surface equation for z.
    // Intputs:
    //     Double_t  ar   The array of R values
    //     Double_t  az   The array of Z values
    //     Int_t     ip   The index in p to get the point location
    //     Double_t  tc   The angle of that part of the cone is at
    //     Double_t   r   The value of Rmax to compute z from
    //     Double_t  th   The perpendicular distance the parralell line is
    //                    from the point ip.
    // Outputs:
    //     none.
    // Return:
    //     The value Z correstponding to the line at angle th, offeset by
    //     th, and the point p->GetZ/Rmax[ip] at the locatin r.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return az[ip]+(ar[ip]+th/costc-r)/tantc;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::ZFromRminpCone(const TGeoPcon *p,int ip,
                                           Double_t tc,Double_t r,
                                           Double_t th)const{
    // General Inner cone Surface equation for z.
    // Intputs:
    //     TGeoPcon  *p   The poly cone where the initial point comes from
    //     Int_t     ip   The index in p to get the point location
    //     Double_t  tc   The angle of that part of the cone is at
    //     Double_t   r   The value of Rmin to compute z from
    //     Double_t  th   The perpendicular distance the parralell line is
    //                    from the point ip.
    // Outputs:
    //     none.
    // Return:
    //     The value Z correstponding to the line at angle th, offeset by
    //     th, and the point p->GetZ/Rmin[ip] at the location r.
    Double_t tantc = TMath::Tan(tc*TMath::DegToRad());
    Double_t costc = TMath::Cos(tc*TMath::DegToRad());

    return p->GetZ(ip)+(p->GetRmin(ip)+th/costc-r)/tantc;
}
//______________________________________________________________________
void AliITSv11Geometry::RadiusOfCurvature(Double_t rc,Double_t theta0,
                                          Double_t z0,Double_t r0,
                                          Double_t theta1,Double_t &z1,
                                          Double_t &r1)const{
    // Given a initial point z0,r0, the initial angle theta0, and the radius
    // of curvature, returns the point z1, r1 at the angle theta1. Theta
    // measured from the r axis in the clock wise direction [degrees].
    // Inputs:
    //    Double_t rc     The radius of curvature
    //    Double_t theta0 The starting angle (degrees)
    //    Double_t z0     The value of z at theta0
    //    Double_t r0     The value of r at theta0
    //    Double_t theta1 The ending angle (degrees)
    // Outputs:
    //    Double_t &z1  The value of z at theta1
    //    Double_t &r1  The value of r at theta1
    // Return:
    //    none.

    z1 = rc*(TMath::Sin(theta1*TMath::DegToRad())-TMath::Sin(theta0*TMath::DegToRad()))+z0;
    r1 = rc*(TMath::Cos(theta1*TMath::DegToRad())-TMath::Cos(theta0*TMath::DegToRad()))+r0;
    return;
}
//______________________________________________________________________
void AliITSv11Geometry::InsidePoint(const TGeoPcon *p,Int_t i1,Int_t i2,
                                    Int_t i3,Double_t c,TGeoPcon *q,Int_t j1,
                                    Bool_t max)const{
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
void AliITSv11Geometry::InsidePoint(Double_t x0,Double_t y0,
                                    Double_t x1,Double_t y1,
                                    Double_t x2,Double_t y2,Double_t c,
                                    Double_t &x,Double_t &y)const{
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
    Double_t dx01,dx12,dy01,dy12,r01,r12,m;
    dx01 = x0-x1; //cout <<"L410 dx01="<<dx01<<endl;
    dx12 = x1-x2; //cout <<"L411 dx12="<<dx12<<endl;
    dy01 = y0-y1; //cout <<"L412 dy01="<<dy01<<endl;
    dy12 = y1-y2; //cout <<"L413 dy12="<<dy12<<endl;
    r01  = TMath::Sqrt(dy01*dy01+dx01*dx01); //cout <<"L414 r01="<<r01<<endl;
    r12  = TMath::Sqrt(dy12*dy12+dx12*dx12); //cout <<"L415 r12="<<r12<<endl;
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
            x = x1-0.5*c*r01/dy01; //cout <<"L434 x="<<x<<endl;
            y = y1+0.5*c*r01/dx01; //cout <<"L435 y="<<y<<endl;
        } // end if
        return;
    } //
    x = x1+c*(dx12*r01-dx01*r12)/m; //cout <<"L442 x="<<x<<endl;
    y = y1+c*(dy12*r01-dy01*r12)/m; //cout <<"L443 y="<<y<<endl;
    //cout <<"=============================================="<<endl;
    return;
}
//----------------------------------------------------------------------
void AliITSv11Geometry:: PrintArb8(const TGeoArb8 *a)const{
    // Prints out the content of the TGeoArb8. Usefull for debugging.
    // Inputs:
    //   TGeoArb8 *a
    // Outputs:
    //   none.
    // Return:
    //   none.

    if(!GetDebug()) return;
    printf("%s",a->GetName());
    a->InspectShape();
    return;
}
//----------------------------------------------------------------------
void AliITSv11Geometry:: PrintPcon(const TGeoPcon *a)const{
    // Prints out the content of the TGeoPcon. Usefull for debugging.
    // Inputs:
    //   TGeoPcon *a
    // Outputs:
    //   none.
    // Return:
    //   none.
  
    if(!GetDebug()) return;
    cout << a->GetName() << ": N=" << a->GetNz() << " Phi1=" << a->GetPhi1()
         << ", Dphi=" << a->GetDphi() << endl;
    cout << "i\t   Z   \t  Rmin \t  Rmax" << endl;
    for(Int_t iii=0;iii<a->GetNz();iii++){
        cout << iii << "\t" << a->GetZ(iii) << "\t" << a->GetRmin(iii)
             << "\t" << a->GetRmax(iii) << endl;
    } // end for iii
    return;
}
//----------------------------------------------------------------------
void AliITSv11Geometry::PrintTube(const TGeoTube *a)const{
    // Prints out the content of the TGeoTube. Usefull for debugging.
    // Inputs:
    //   TGeoTube *a
    // Outputs:
    //   none.
    // Return:
    //   none.

    if(!GetDebug()) return;
    cout << a->GetName() <<": Rmin="<<a->GetRmin()
         <<" Rmax=" <<a->GetRmax()<<" Dz="<<a->GetDz()<<endl;
    return;
}
//----------------------------------------------------------------------
void AliITSv11Geometry::PrintTubeSeg(const TGeoTubeSeg *a)const{
    // Prints out the content of the TGeoTubeSeg. Usefull for debugging.
    // Inputs:
    //   TGeoTubeSeg *a
    // Outputs:
    //   none.
    // Return:
    //   none.

    if(!GetDebug()) return;
    cout << a->GetName() <<": Phi1="<<a->GetPhi1()<<
        " Phi2="<<a->GetPhi2()<<" Rmin="<<a->GetRmin()
         <<" Rmax=" <<a->GetRmax()<<" Dz="<<a->GetDz()<<endl;
    return;
}
//----------------------------------------------------------------------
void AliITSv11Geometry::PrintConeSeg(const TGeoConeSeg *a)const{
    // Prints out the content of the TGeoConeSeg. Usefull for debugging.
    // Inputs:
    //   TGeoConeSeg *a
    // Outputs:
    //   none.
    // Return:
    //   none.

    if(!GetDebug()) return;
    cout << a->GetName() <<": Phi1="<<a->GetPhi1()<<
        " Phi2="<<a->GetPhi2()<<" Rmin1="<<a->GetRmin1()
         <<" Rmax1=" <<a->GetRmax1()<<" Rmin2="<<a->GetRmin2()
         <<" Rmax2=" <<a->GetRmax2()<<" Dz="<<a->GetDz()<<endl;
    return;
}
//----------------------------------------------------------------------
void AliITSv11Geometry::PrintBBox(const TGeoBBox *a)const{
    // Prints out the content of the TGeoBBox. Usefull for debugging.
    // Inputs:
    //   TGeoBBox *a
    // Outputs:
    //   none.
    // Return:
    //   none.

    if(!GetDebug()) return;
    cout << a->GetName() <<": Dx="<<a->GetDX()<<
        " Dy="<<a->GetDY()<<" Dz="<<a->GetDZ() <<endl;
    return;
}
//---------------------------------------------------------------------
void AliITSv11Geometry::DrawCrossSection(const TGeoPcon *p,
                            Int_t fillc,Int_t fills,
                            Int_t linec,Int_t lines,Int_t linew,
                            Int_t markc,Int_t marks,Float_t marksize)const{
    // Draws a cross sectional view of the TGeoPcon, Primarily for debugging.
    // A TCanvas should exist first.
    //  Inputs:
    //    TGeoPcon  *p  The TGeoPcon to be "drawn"
    //    Int_t  fillc  The fill color to be used
    //    Int_t  fills  The fill style to be used
    //    Int_t  linec  The line color to be used
    //    Int_t  lines  The line style to be used
    //    Int_t  linew  The line width to be used
    //    Int_t  markc  The markder color to be used
    //    Int_t  marks  The markder style to be used
    //    Float_t marksize The marker size
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t n=0,m=0,i=0;
    Double_t *z=0,*r=0;
    TPolyMarker *pts=0;
    TPolyLine   *line=0;

    n = p->GetNz();
    if(n<=0) return;
    m = 2*n+1;
    z = new Double_t[m];
    r = new Double_t[m];

    for(i=0;i<n;i++){
        z[i] = p->GetZ(i);
        r[i] = p->GetRmax(i);
        z[i+n] = p->GetZ(n-1-i);
        r[i+n] = p->GetRmin(n-1-i);
    } //  end for i
    z[n-1] = z[0];
    r[n-1] = r[0];

    line = new TPolyLine(n,z,r);
    pts  = new TPolyMarker(n,z,r);

    line->SetFillColor(fillc);
    line->SetFillStyle(fills);
    line->SetLineColor(linec);
    line->SetLineStyle(lines);
    line->SetLineWidth(linew);
    pts->SetMarkerColor(markc);
    pts->SetMarkerStyle(marks);
    pts->SetMarkerSize(marksize);

    line->Draw("f");
    line->Draw();
    pts->Draw();

    delete[] z;
    delete[] r;

    cout<<"Hit Return to continue"<<endl;
    cin >> n;
    delete line;
    delete pts;
    return;
}
//----------------------------------------------------------------------
