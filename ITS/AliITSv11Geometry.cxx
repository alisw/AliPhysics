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
#include <TArc.h>
#include <TLine.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TText.h>
#include <TObjArray.h>
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TArrayI.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoMaterial.h>
#include <TPolyMarker.h>
#include <TPolyLine.h>
#include "AliITSv11Geometry.h"

ClassImp(AliITSv11Geometry)

const Double_t AliITSv11Geometry::fgkmicron = 1.0E-4;
const Double_t AliITSv11Geometry::fgkmm = 0.10;
const Double_t AliITSv11Geometry::fgkcm = 1.00;
const Double_t AliITSv11Geometry::fgkDegree = 1.0;
const Double_t AliITSv11Geometry::fgkRadian = 180./3.14159265358979323846;

//----------------------------------------------------------------------
TGeoMixture * AliITSv11Geometry::CreateMixtureByVolume(const char* name,
                                 Int_t nel,const TArrayD *v,
                                 const TObjArray *mix,Double_t den){
    // Create a new TGeoMixture object based on a TObject array of
    // TGeoMixture/TGeoMaterials and their releative weight volume. For 
    // example, Consider TGeoMixture of 
    // Inputs:
    //    const char*     name Name of the new TGeoMixture
    //          Int_t     nel  The number of enteries
    //    const TArrayD   *v   Array of releative volume of each mixture
    //    const TObjArray *mix TObjArray holding the TGeoMixtures or 
    //                         TGeoMaterials which make up this new mixture.
    //          Double_t   den The density of this new mixture [g/cm^3].
    // Output:
    //    none.
    // Retrun:
    //    A pointer to a new instance of a TGeoMixture.
    Int_t i;
    TArrayD w(nel);

    for(i=0;i<nel;i++) if(v->At(i)>=0.0){
        w[i] = v->At(i) * ((TGeoMaterial*)(mix->At(i)))->GetDensity();
    }else{
        w[i] = 0.0;
    } // end if/for i
    return CreateMixtureByWeight(name,nel,&w,mix,den);
}
//______________________________________________________________________
TGeoMixture * AliITSv11Geometry::CreateMixtureByNumber(const char* name,
                                 Int_t nel,const TArrayI *w,
                                 const TObjArray *mix,Double_t den){
    // Create a new TGeoMixture object based on a TObject array of
    // TGeoMixture/TGeoMaterials and their releative number. For example,
    // Consider TGeoMixture of BGO (Bi_2O_3)_2(GeO_2)_3. Assume you have
    // defined Bismuth Oxide as Bi_2O_3 as one TGeoMixture ("Bi2O3") and 
    // Germainium Oxide as another ("GeO2"), then BGO is defined at
    // CreateMixtureByNumber("BGO",2,TArrayI(2,3),
    // TObjArray(TGeoMixture("Bi2O3"),TGeoMixture("GeO2")));
    // Inputs:
    //    const char* name     Name of the new TGeoMixture
    //          Int_t nel      The number of enteries
    //    const TArrayI *w     Array of releative number of each mixture
    //    const TObjArray *mix TObjArray holding the TGeoMixtures or 
    //                         TGeoMaterials which make up this new mixture.
    //          Double_t   den The density of this new mixture [g/cm^3].
    // Output:
    //    none.
    // Retrun:
    //    A pointer to a new instance of a TGeoMixture.
    Int_t i,j;
    Double_t a;
    TArrayD wa;
    TGeoMixture *mi;
    TGeoMaterial *ma;

    wa.Set(nel);
    wa.Reset();
    for(i=0;i<nel;i++) if(w->At(i)>0){
        mi = dynamic_cast<TGeoMixture*>(mix->At(i));
        if(mi!=0){ // Mixture
            a = 0.0;
            for(j=0;j<mi->GetNelements();j++) a+= mi->GetZmixt()[j];
            wa[i] =  ((Double_t)(w->At(i)))*a/((Double_t)(mi->GetNelements()));
        }else{ // Material
            ma = dynamic_cast<TGeoMaterial*>(mix->At(i));
            if(ma==0) continue; // don't know what this is.
            wa[i] = ((Double_t)(w->At(i)))*ma->GetA();
        } // end if
    } // end if/for i

    return CreateMixtureByWeight(name,nel,&wa,mix,den);
}
//----------------------------------------------------------------------
TGeoMixture * AliITSv11Geometry::CreateMixtureByWeight(const char* name,
                                 Int_t nel,const TArrayD *w,
                                 const TObjArray *mix,Double_t den){
    // Create a new TGeoMixture object based on a TObject array of
    // TGeoMixture/TGeoMaterials and their releative weight. For example,
    // Consider TGeoMixture of "standard shielding blocks" 52% O_2, 32.5% Si,
    // 6% Ca, 1.5% Na, 2% Fe, and 4% Al.
    // Inputs:
    //    const char* name     Name of the new TGeoMixture
    //          Int_t nel      The number of enteries
    //    const TArrayD *w     Array of releative Weights of each mixture
    //    const TObjArray *mix TObjArray holding the TGeoMixtures or 
    //                         TGeoMaterials which make up this new mixture.
    //          Double_t   den The density of this new mixture [g/cm^3].
    // Output:
    //    none.
    // Retrun:
    //    A pointer to a new instance of a TGeoMixture.
    Int_t i,j,k,n=0;
    Double_t s=0.0;
    TArrayD za,aa,wa;
    TGeoMixture *mi,*mixnew;
    TGeoMaterial *ma;
    Bool_t add=kTRUE;

    for(i=0;i<nel;i++) if(w->At(i)>0){
        mi = dynamic_cast<TGeoMixture*>(mix->At(i));
        if(mi!=0){
            for(j=0;j<mi->GetNelements();j++) n++;
        }else{
            ma = dynamic_cast<TGeoMaterial*>(mix->At(i));
            if(ma==0) continue; // don't know what this is.
            n++;
        } // end if
        s += w->At(i);
    } // end for i
    if(n<=0) return 0; // No elements found.
    za.Set(n); za.Reset();
    aa.Set(n); aa.Reset();
    wa.Set(n); wa.Reset();
    TMatrixD wb(n,nel); wb.Zero();
    //
    n = 0; // Now use for the number of enteries.
    for(i=0;i<nel;i++) {
        for(j=0;j<n;j++) wb(j,i) = 0.0; // zero out array
        if(w->At(i)>0){
            mi = dynamic_cast<TGeoMixture*>(mix->At(i));
            if(mi!=0){
                for(j=0;j<mi->GetNelements();j++){
                    for(k=0;k<n;k++){
                        if(za.At(k)==mi->GetZmixt()[j] && 
                           aa.At(k)==mi->GetAmixt()[j]){
                            add = kFALSE;
                            wb(k,i) = mi->GetWmixt()[j] * w->At(i)/s;
                            continue;
                        } // end if
                        if(add){
                            za[n]   = mi->GetZmixt()[j];
                            aa[n]   = mi->GetAmixt()[j];
                            wb(n,i) = mi->GetWmixt()[j] * w->At(i)/s;
                            n++;
                        } // end if
                        add = kTRUE;
                    } // end for k
                } // end for j
            }else{
                ma = dynamic_cast<TGeoMaterial*>(mix->At(i));
                if(ma==0) continue; // don't know what this is.
                for(k=0;k<n;k++){
                    if(za.At(k)==ma->GetZ() && aa.At(k)==ma->GetA()){
                        add = kFALSE;
                        wb(k,i) = w->At(i)/s;
                        continue;
                    } // end if
                    if(add){
                        za[n] = ma->GetZ();
                        aa[n] = ma->GetA();
                        wb(n,i) = w->At(i)/s;
                        n++;
                    } // end if
                    add = kTRUE;
                } // end for k
            } // end if
        } // end if
    } // end for i
    mixnew = new TGeoMixture(name,n,den);
    k = 0;
    for(i=0;i<n;i++) {
        for(j=0;j<nel;j++) wa.AddAt(wb(i,j),i);
        if(wa.At(i)<=0.0) continue;
        mixnew->DefineElement(k++,aa.At(i),za.At(i),wa.At(i));
    } // end for i
    //
    return mixnew;
}
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
//______________________________________________________________________
Bool_t AliITSv11Geometry::AngleOfIntersectionWithLine(Double_t x0,Double_t y0,
                                                      Double_t x1,Double_t y1,
                                                      Double_t xc,Double_t yc,
                                                      Double_t rc,Double_t &t0,
                                                      Double_t &t1)const{
    // Computes the angles, t0 and t1 corresponding to the intersection of
    // the line, defined by {x0,y0} {x1,y1}, and the circle, defined by
    // its center {xc,yc} and radius r. If the line does not intersect the
    // line, function returns kFALSE, otherwise it returns kTRUE. If the
    // line is tangent to the circle, the angles t0 and t1 will be the same.
    // Inputs:
    //   Double_t x0   X of first point defining the line
    //   Double_t y0   Y of first point defining the line
    //   Double_t x1   X of Second point defining the line
    //   Double_t y1   Y of Second point defining the line
    //   Double_t xc   X of Circle center point defining the line
    //   Double_t yc   Y of Circle center point defining the line
    //   Double_t r    radius of circle
    // Outputs:
    //   Double_t &t0  First angle where line intersects circle
    //   Double_t &t1  Second angle where line intersects circle
    // Return:
    //    kTRUE, line intersects circle, kFALSE line does not intersect circle
    //           or the line is not properly defined point {x0,y0} and {x1,y1}
    //           are the same point.
    Double_t dx,dy,cx,cy,s2,t[4];
    Double_t a0,b0,c0,a1,b1,c1,sinthp,sinthm,costhp,costhm;
    Int_t i,j;

    t0 = 400.0;
    t1 = 400.0;
    dx = x1-x0;
    dy = y1-y0;
    cx = xc-x0;
    cy = yc-y0;
    s2 = dx*dx+dy*dy;
    if(s2==0.0) return kFALSE;

    a0 = rc*rc*s2;
    if(a0==0.0) return kFALSE;
    b0 = 2.0*rc*dx*(dx*cy-cx*dy);
    c0 = dx*dx*cy*cy-2.0*dy*dx*cy*cx+cx*cx*dy*dy-rc*rc*dy*dy;
    c0 = 0.25*b0*b0/(a0*a0)-c0/a0;
    if(c0<0.0) return kFALSE;
    sinthp = -0.5*b0/a0+TMath::Sqrt(c0);
    sinthm = -0.5*b0/a0-TMath::Sqrt(c0);

    a1 = rc*rc*s2;
    if(a1==0.0) return kFALSE;
    b1 = 2.0*rc*dy*(dy*cx-dx*cy);
    c1 = dy*dy*cx*cx-2.0*dy*dx*cy*cx+dx*dx*cy*cy-rc*rc*dx*dx;
    c1 = 0.25*b1*b1/(a1*a1)-c1/a1;
    if(c1<0.0) return kFALSE;
    costhp = -0.5*b1/a1+TMath::Sqrt(c1);
    costhm = -0.5*b1/a1-TMath::Sqrt(c1);

    t[0] = t[1] = t[2] = t[3] = 400.;
    a0 = TMath::ATan2(sinthp,costhp); if(a0<0.0) a0 += 2.0*TMath::Pi();
    a1 = TMath::ATan2(sinthp,costhm); if(a1<0.0) a1 += 2.0*TMath::Pi();
    b0 = TMath::ATan2(sinthm,costhp); if(b0<0.0) b0 += 2.0*TMath::Pi();
    b1 = TMath::ATan2(sinthm,costhm); if(b1<0.0) b1 += 2.0*TMath::Pi();
    x1 = xc+rc*TMath::Cos(a0);
    y1 = yc+rc*TMath::Sin(a0);
    s2 = dx*(y1-y0)-dy*(x1-x0);
    if(s2*s2<DBL_EPSILON) t[0] = a0*TMath::RadToDeg();
    x1 = xc+rc*TMath::Cos(a1);
    y1 = yc+rc*TMath::Sin(a1);
    s2 = dx*(y1-y0)-dy*(x1-x0);
    if(s2*s2<DBL_EPSILON) t[1] = a1*TMath::RadToDeg();
    x1 = xc+rc*TMath::Cos(b0);
    y1 = yc+rc*TMath::Sin(b0);
    s2 = dx*(y1-y0)-dy*(x1-x0);
    if(s2*s2<DBL_EPSILON) t[2] = b0*TMath::RadToDeg();
    x1 = xc+rc*TMath::Cos(b1);
    y1 = yc+rc*TMath::Sin(b1);
    s2 = dx*(y1-y0)-dy*(x1-x0);
    if(s2*s2<DBL_EPSILON) t[3] = b1*TMath::RadToDeg();
    for(i=0;i<4;i++)for(j=i+1;j<4;j++){
        if(t[i]>t[j]) {t0 = t[i];t[i] = t[j];t[j] = t0;}
    } // end for i,j
    t0 = t[0];
    t1 = t[1];
    //
    return kTRUE;
}
//______________________________________________________________________
Double_t AliITSv11Geometry::AngleForRoundedCorners0(Double_t dx,Double_t dy,
                                                    Double_t sdr)const{
    // Basic function used to determine the ending angle and starting angles
    // for rounded corners given the relative distance between the centers
    // of the circles and the difference/sum of their radii. Case 0.
    // Inputs:
    //   Double_t dx    difference in x locations of the circle centers
    //   Double_t dy    difference in y locations of the circle centers
    //   Double_t sdr   difference or sum of the circle radii
    // Outputs:
    //   none.
    // Return:
    //   the angle in Degrees
    Double_t a,b;

    b = dy*dy+dx*dx-sdr*sdr;
    if(b<0.0) Error("AngleForRoundedCorners0",
                    "dx^2(%e)+dy^2(%e)-sdr^2(%e)=b=%e<0",dx,dy,sdr,b);
    b = TMath::Sqrt(b);
    a = -sdr*dy+dx*b;
    b = -sdr*dx-dy*b;
    return TMath::ATan2(a,b)*TMath::RadToDeg();
    
}
//______________________________________________________________________
Double_t AliITSv11Geometry::AngleForRoundedCorners1(Double_t dx,Double_t dy,
                                                    Double_t sdr)const{
    // Basic function used to determine the ending angle and starting angles
    // for rounded corners given the relative distance between the centers
    // of the circles and the difference/sum of their radii. Case 1.
    // Inputs:
    //   Double_t dx    difference in x locations of the circle centers
    //   Double_t dy    difference in y locations of the circle centers
    //   Double_t sdr   difference or sum of the circle radii
    // Outputs:
    //   none.
    // Return:
    //   the angle in Degrees
    Double_t a,b;

    b = dy*dy+dx*dx-sdr*sdr;
    if(b<0.0) Error("AngleForRoundedCorners1",
                    "dx^2(%e)+dy^2(%e)-sdr^2(%e)=b=%e<0",dx,dy,sdr,b);
    b = TMath::Sqrt(b);
    a = -sdr*dy-dx*b;
    b = -sdr*dx+dy*b;
    return TMath::ATan2(a,b)*TMath::RadToDeg();
    
}
//----------------------------------------------------------------------
void AliITSv11Geometry::AnglesForRoundedCorners(Double_t x0,Double_t y0,
                                                Double_t r0,Double_t x1,
                                                Double_t y1,Double_t r1,
                                                Double_t &t0,Double_t &t1)
    const{
    // Function to compute the ending angle, for arc 0, and starting angle,
    // for arc 1, such that a straight line will connect them with no
    // discontinuities.
    //Begin_Html
    /*
      <img src="picts/ITS/AliITSv11Geometry_AnglesForRoundedCorners.gif">
     */
    //End_Html
    // Inputs:
    //    Double_t x0  X Coordinate of arc 0 center.
    //    Double_t y0  Y Coordinate of arc 0 center.
    //    Double_t r0  Radius of curvature of arc 0. For signe see figure.
    //    Double_t x1  X Coordinate of arc 1 center.
    //    Double_t y1  Y Coordinate of arc 1 center.
    //    Double_t r1  Radius of curvature of arc 1. For signe see figure.
    // Outputs:
    //    Double_t t0  Ending angle of arch 0, with respect to x axis, Degrees.
    //    Double_t t1  Starting angle of arch 1, with respect to x axis, 
    //                 Degrees.
    // Return:
    //    none.
    Double_t t;

    if(r0>=0.0&&r1>=0.0) { // Inside to inside    ++
        t = AngleForRoundedCorners1(x1-x0,y1-y0,r1-r0);
        t0 = t1 = t;
        return;
    }else if(r0>=0.0&&r1<=0.0){ // Inside to Outside  +-
        r1 = -r1; // make positive
        t = AngleForRoundedCorners0(x1-x0,y1-y0,r1+r0);
        t0 = 180.0 + t;
        if(t0<0.0) t += 360.;
        if(t<0.0) t += 360.;
        t1 = t;
        return;
    }else if(r0<=0.0&&r1>=0.0){ // Outside to Inside  -+
        r0 = - r0; // make positive
        t = AngleForRoundedCorners1(x1-x0,y1-y0,r1+r0);
        t0 = 180.0 + t;
        if(t0>180.) t0 -= 360.;
        if(t >180.) t  -= 360.;
        t1 = t;
        return;
    }else if(r0<=0.0&&r1<=0.0) { // Outside to outside --
        r0 = -r0; // make positive
        r1 = -r1; // make positive
        t = AngleForRoundedCorners0(x1-x0,y1-y0,r1-r0);
        t0 = t1 = t;
        return;
    } // end if
    return;
}
//----------------------------------------------------------------------
void AliITSv11Geometry::MakeFigure1(Double_t x0,Double_t y0,Double_t r0,
                                    Double_t x1,Double_t y1,Double_t r1){
    // Function to create the figure discribing how the function 
    // AnglesForRoundedCorners works.
    //
    // Inputs:
    //    Double_t x0  X Coordinate of arc 0 center.
    //    Double_t y0  Y Coordinate of arc 0 center.
    //    Double_t r0  Radius of curvature of arc 0. For signe see figure.
    //    Double_t x1  X Coordinate of arc 1 center.
    //    Double_t y1  Y Coordinate of arc 1 center.
    //    Double_t r1  Radius of curvature of arc 1. For signe see figure.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Double_t t0[4],t1[4],xa0[4],ya0[4],xa1[4],ya1[4],ra0[4],ra1[4];
    Double_t xmin,ymin,xmax,ymax,h;
    Int_t j;

    for(j=0;j<4;j++) {
        ra0[j] = r0; if(j%2) ra0[j] = -r0;
        ra1[j] = r1; if(j>1) ra1[j] = -r1;
        AnglesForRoundedCorners(x0,y0,ra0[j],x1,y1,ra1[j],t0[j],t1[j]);
        xa0[j] = TMath::Abs(r0)*CosD(t0[j])+x0;
        ya0[j] = TMath::Abs(r0)*SinD(t0[j])+y0;
        xa1[j] = TMath::Abs(r1)*CosD(t1[j])+x1;
        ya1[j] = TMath::Abs(r1)*SinD(t1[j])+y1;
    } // end for j
    if(r0<0.0) r0 = -r0;
    if(r1<0.0) r1 = -r1;
    xmin = TMath::Min(x0 - r0,x1-r1);
    ymin = TMath::Min(y0 - r0,y1-r1);
    xmax = TMath::Max(x0 + r0,x1+r1);
    ymax = TMath::Max(y0 + r0,y1+r1);
    for(j=1;j<4;j++) {
        xmin = TMath::Min(xmin,xa0[j]);
        xmin = TMath::Min(xmin,xa1[j]);
        ymin = TMath::Min(ymin,ya0[j]);
        ymin = TMath::Min(ymin,ya1[j]);

        xmax = TMath::Max(xmax,xa0[j]);
        xmax = TMath::Max(xmax,xa1[j]);
        ymax = TMath::Max(ymax,ya0[j]);
        ymax = TMath::Max(ymax,ya1[j]);
    } // end for j
    if(xmin<0.0) xmin *= 1.1; else xmin *= 0.9;
    if(ymin<0.0) ymin *= 1.1; else ymin *= 0.9;
    if(xmax<0.0) xmax *= 0.9; else xmax *= 1.1;
    if(ymax<0.0) ymax *= 0.9; else ymax *= 1.1;
    j = (Int_t)(500.0*(ymax-ymin)/(xmax-xmin));
    TCanvas *can = new TCanvas("AliITSv11Geometry_AnglesForRoundedCorners",
                               "Figure for AliITSv11Geometry",500,j);
    h = ymax-ymin; if(h<0) h = -h;
    can->Range(xmin,ymin,xmax,ymax);
    TArc *c0 = new TArc(x0,y0,r0);
    TArc *c1 = new TArc(x1,y1,r1);
    TLine *line[4];
    TArrow *ar0[4];
    TArrow *ar1[4];
    for(j=0;j<4;j++){
        ar0[j] = new TArrow(x0,y0,xa0[j],ya0[j]);
        ar1[j] = new TArrow(x1,y1,xa1[j],ya1[j]);
        line[j] = new TLine(xa0[j],ya0[j],xa1[j],ya1[j]);
        ar0[j]->SetLineColor(j+1);
        ar0[j]->SetArrowSize(0.1*r0/h);
        ar1[j]->SetLineColor(j+1);
        ar1[j]->SetArrowSize(0.1*r1/h);
        line[j]->SetLineColor(j+1);
    } // end for j
    c0->Draw();
    c1->Draw();
    for(j=0;j<4;j++){
        ar0[j]->Draw();
        ar1[j]->Draw();
        line[j]->Draw();
    } // end for j
    TText *t = new TText();
    t->SetTextSize(0.02);
    Char_t txt[100];
    sprintf(txt,"(x0=%5.2f,y0=%5.2f)",x0,y0);
    t->DrawText(x0,y0,txt);
    sprintf(txt,"(x1=%5.2f,y1=%5.2f)",x1,y1);
    for(j=0;j<4;j++) {
        t->SetTextColor(j+1);
        t->DrawText(x1,y1,txt);
        sprintf(txt,"r0=%5.2f",ra0[j]);
        t->DrawText(0.5*(x0+xa0[j]),0.5*(y0+ya0[j]),txt);
        sprintf(txt,"r1=%5.2f",ra1[j]);
        t->DrawText(0.5*(x1+xa1[j]),0.5*(y1+ya1[j]),txt);
    } // end for j
}
