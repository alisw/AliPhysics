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

// This class Defines the Geometry for the ITS services and support cones
// outside of the ceneteral volume (except for the Ceneteral support 
// cylinders. Other classes define the rest of the ITS. Specificaly the ITS
// The SSD support cone,SSD Support centeral cylinder, SDD support cone,
// The SDD cupport centeral cylinder, the SPD Thermal Sheald, The supports
// and cable trays on both the RB26 (muon dump) and RB24 sides, and all of
// the cabling from the ladders/stave ends out past the TPC. 

/* $Id$ */
// General Root includes
#include <TMath.h>
// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include "AliITSv11GeometrySupport.h"

ClassImp(AliITSv11GeometrySupport)

#define SQ(A) (A)*(A)

//______________________________________________________________________
void AliITSv11GeometrySupport::SPDCone(TGeoVolume *moth){
    // Define the detail SPD support cone geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.

    SPDThermalSheald(moth);
}
//______________________________________________________________________
void AliITSv11GeometrySupport::SPDThermalSheald(TGeoVolume *moth){
    // Define the detail SPD Thermal Sheld geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
    // From ALICE-Thermal Screen (SPD) "Cylinder" file thermal-screen2_a3.ps
    // Volumes sA1,sA2,sA3,sAh1,sAh2,sAh3, and b1,b2,b3,bh1,bh2,bh3;
    // "CONE TRANSITION" file thermal-screen1_a3.ps Volumes sC1,sC2,sC3,
    // sCh1,sCh2, sCh3; "FLANGE" file thermal-screen4_a3.ps Volumes d,sDs,
    // sDw,sDws; and "HALF ASSEMBLY" file thermal-screen3_a3.ps. This object,
    // both halfs, are incased inside of a single minimum sized mother 
    // volume called M, which is a union of two parts sM1 and 4 copies of sM2.
    const Double_t ktscarbonFiberThA = 0.03*fgkmm; // 
    //const Double_t ktscarbonFiberThB = 0.10*fgkmm; //
    const Double_t ktscLengthB  = 50.0*fgkmm; //
    const Double_t ktscLengthA  = 900.0*fgkmm-2.0*ktscLengthB; //
    const Double_t ktscLengthC  = 290.0*fgkmm; //
    const Double_t ktscLengthD  = 15.0*fgkmm; //
    const Double_t ktscAngle    = 36.0*fgkDegree;//Rep. angle of cent. accordin
    const Double_t ktscRoutA    = 99.255*fgkmm; // Outer radii
    const Double_t ktscRinA     = 81.475*fgkmm; // Iner radii
    const Double_t ktscRoutB    = 99.955*fgkmm; // Outer radii
    const Double_t ktscRinB     = 80.775*fgkmm; // Iner radii
    const Double_t ktscRoutCp   = 390.0*fgkmm;  // Outer radii
    const Double_t ktscRinCp    = 373.0*fgkmm;  // Iner radii
    Double_t ktscRoutC,ktscRinC; // values need to be calculated
    const Double_t ktscRwingD   = 492.5*fgkmm;  // Outer radii
    const Double_t ktscRoutD    = 0.5*840.*fgkmm;// Outer radii
    const Double_t ktscRinD     = 373.0*fgkmm;  // Iner radii
    // angular wing
    const Double_t ktscAngleDD  = (60.*fgkmm/ktscRwingD)*fgkRadian;
                                                    // width of fill material
    const Double_t ktscAngleDDs = ((60.*fgkmm-2.*ktscarbonFiberThA)/
                                                  ktscRwingD)*fgkRadian;
    const Double_t ktscAngleD0  = 45.*fgkDegree;//Strting angle of wing
    const Double_t ktscoutSA    = 24.372*fgkmm; // The other one Calculated
    const Double_t ktscinLA     = 31.674*fgkmm; // The ohter one Calculated
    const Double_t ktscoutSB    = 24.596*fgkmm; // The other one Calculated
    const Double_t ktscinLB     = 31.453*fgkmm; // The ohter one Calculated
    const Double_t ktscoutSC    = 148.831*fgkmm;// The other one Calculated
    const Double_t ktscinLC     = 90.915*fgkmm; // The ohter one Calculated
    Int_t i,k;
    Double_t th;
    Double_t xo[7],yo[7],xi[7],yi[7];
    Double_t xbo[7],ybo[7],xbi[7],ybi[7];
    Double_t xco[7],yco[7],xci[7],yci[7];
    TGeoArb8 *sA1,*sA2,*sA3,*sAh1,*sAh2,*sAh3,*sB1,*sB2,*sB3,*sBh1,*sBh2,*sBh3;
    TGeoArb8 *sC1,*sC2,*sC3,*sCh1,*sCh2,*sCh3;
    TGeoPcon *sM1;
    TGeoTube  *sD,*sDs;
    TGeoTubeSeg *sDw,*sDws,*sM2;
    TGeoCompositeShape *sM;
    TGeoRotation *rot;
    TGeoTranslation *tranb,*tranbm,*tranc;
    TGeoTranslation *tranITSspdShealdVVt0;
    TGeoCombiTrans *rotITSspdShealdVVt1,*rotITSspdShealdVVt2;
    TGeoCombiTrans *rotITSspdShealdVVt3;
    TGeoMedium *medSPDcf  = 0; // SPD support cone Carbon Fiber materal number.
    TGeoMedium *medSPDfs  = 0; // SPD support cone inserto stesalite 4411w.
    TGeoMedium *medSPDfo  = 0; // SPD support cone foam, Rohacell 50A.
    TGeoMedium *medSPDss  = 0; // SPD support cone screw material,Stainless
    TGeoMedium *medSPDair = 0; // SPD support cone Air
    //TGeoMedium *medSPDal  = 0; // SPD support cone SDD mounting bracket Al

    ktscRoutC = TMath::Sqrt(ktscRoutCp*ktscRoutCp-0.25*ktscoutSC*ktscoutSC);
    ktscRinC  = TMath::Sqrt(ktscRinCp *ktscRinCp -0.25*ktscinLC *ktscinLC );
    sA1  = new TGeoArb8("ITS SPD Therm Screen Clyinder A1",0.5*ktscLengthA);
    sA2  = new TGeoArb8("ITS SPD Therm Screen Clyinder A2",0.5*ktscLengthA);
    sA3  = new TGeoArb8("ITS SPD Therm Screen Clyinder A3",0.5*ktscLengthA);
    sAh1 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ah1",0.5*ktscLengthA);
    sAh2 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ah2",0.5*ktscLengthA);
    sAh3 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ah3",0.5*ktscLengthA);
    sB1  = new TGeoArb8("ITS SPD Therm Screen Clyinder B1",0.5*ktscLengthB);
    sB2  = new TGeoArb8("ITS SPD Therm Screen Clyinder B2",0.5*ktscLengthB);
    sB3  = new TGeoArb8("ITS SPD Therm Screen Clyinder B3",0.5*ktscLengthB);
    sBh1 = new TGeoArb8("ITS SPD Therm Screen Cylinder Bh1",0.5*ktscLengthB);
    sBh2 = new TGeoArb8("ITS SPD Therm Screen Cylinder Bh2",0.5*ktscLengthB);
    sBh3 = new TGeoArb8("ITS SPD Therm Screen Cylinder Bh3",0.5*ktscLengthB);
    sC1  = new TGeoArb8("ITS SPD Therm Screen Clyinder C1",0.5*ktscLengthC);
    sC2  = new TGeoArb8("ITS SPD Therm Screen Clyinder C2",0.5*ktscLengthC);
    sC3  = new TGeoArb8("ITS SPD Therm Screen Clyinder C3",0.5*ktscLengthC);
    sCh1 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ch1",0.5*ktscLengthC);
    sCh2 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ch2",0.5*ktscLengthC);
    sCh3 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ch3",0.5*ktscLengthC);
    sD = new TGeoTube("ITS SPD Therm Screen Flange D",ktscRinD,ktscRoutD,
                    0.5*ktscLengthD);
    sDs = new TGeoTube("ITS SPD Therm Screen Flange fill Ds",
                      ktscRinD+ktscarbonFiberThA,ktscRoutD-ktscarbonFiberThA,
                      0.5*ktscLengthD);
    sDw = new TGeoTubeSeg("ITS SPD Therm Screen Flange Wing Dw",
                         ktscRoutD,ktscRwingD ,0.5*ktscLengthD,
                         ktscAngleD0-0.5*ktscAngleDD,
                         ktscAngleD0+0.5*ktscAngleDD);
    sDws = new TGeoTubeSeg("ITS SPD Therm Screen Flange Wing Fill Ds",
                          ktscRoutD,ktscRwingD-ktscarbonFiberThA,
                          0.5*ktscLengthD,ktscAngleD0-0.5*ktscAngleDDs,
                          ktscAngleD0+0.5*ktscAngleDDs);
    k = 0;
    for(i=-1;i<2;i++){
        th = ((Double_t)(i+1))*ktscAngle*fgkDegree;
        xo[k]  = ktscRoutA*SinD(th) - 0.5*ktscoutSA*CosD(th);
        yo[k]  = ktscRoutA*CosD(th) + 0.5*ktscoutSA*SinD(th);
        xi[k]  = ktscRinA *SinD(th) - 0.5*ktscinLA *CosD(th);
        yi[k]  = ktscRinA *CosD(th) + 0.5*ktscinLA *SinD(th);
        xbo[k] = ktscRoutB*SinD(th) - 0.5*ktscoutSB*CosD(th);
        ybo[k] = ktscRoutB*CosD(th) + 0.5*ktscoutSB*SinD(th);
        xbi[k] = ktscRinB *SinD(th) - 0.5*ktscinLB *CosD(th);
        ybi[k] = ktscRinB *CosD(th) + 0.5*ktscinLB *SinD(th);
        xco[k] = ktscRoutC*SinD(th) - 0.5*ktscoutSC*CosD(th);
        yco[k] = ktscRoutC*CosD(th) + 0.5*ktscoutSC*SinD(th);
        xci[k] = ktscRinC *SinD(th) - 0.5*ktscinLC *CosD(th);
        yci[k] = ktscRinC *CosD(th) + 0.5*ktscinLC *SinD(th);
        k++;
        xo[k]  = ktscRoutA*SinD(th) + 0.5*ktscoutSA*CosD(th);
        yo[k]  = ktscRoutA*CosD(th) - 0.5*ktscoutSA*SinD(th);
        xi[k]  = ktscRinA *SinD(th) + 0.5*ktscinLA *CosD(th);
        yi[k]  = ktscRinA *CosD(th) - 0.5*ktscinLA *SinD(th);
        xbo[k] = ktscRoutB*SinD(th) + 0.5*ktscoutSB*CosD(th);
        ybo[k] = ktscRoutB*CosD(th) - 0.5*ktscoutSB*SinD(th);
        xbi[k] = ktscRinB *SinD(th) + 0.5*ktscinLB *CosD(th);
        ybi[k] = ktscRinB *CosD(th) - 0.5*ktscinLB *SinD(th);
        xco[k] = ktscRoutC*SinD(th) + 0.5*ktscoutSC*CosD(th);
        yco[k] = ktscRoutC*CosD(th) - 0.5*ktscoutSC*SinD(th);
        xci[k] = ktscRinC *SinD(th) + 0.5*ktscinLC *CosD(th);
        yci[k] = ktscRinC *CosD(th) - 0.5*ktscinLC *SinD(th);
        k++;
    } // end for i
    xo[6] = xo[5];
    yo[6] = 0.0;
    xi[6] = xi[5];
    yi[6] = 0.0;
    xbo[6] = xbo[5];
    ybo[6] = 0.0;
    xbi[6] = xbi[5];
    ybi[6] = 0.0;
    xco[6] = xco[5];
    yco[6] = 0.0;
    xci[6] = xci[5];
    yci[6] = 0.0;
    if(GetDebug()){
        Info("SPDThermalSheald","i     \t  xo  yo    \t  xi yi     \t  xbo "
             "ybo   \t   xbi ybi  \t   xco yco   \t   xci yxi");
        for(i=0;i<7;i++){
            Info("SPDThermalSheald","%7d\t%7.4f,%7.4f\t%7.4f,%7.4f\t"
                 "%7.4f,%7.4f\t%7.4f,%7.4f\t%7.4f,%7.4f\t%7.4f,%7.4f",i,
                 xo[i],yo[i],xi[i],yi[i],
                 xbo[i],ybo[i],xbi[i],ybi[i],
                 xco[i],yco[i],xci[i],yci[i]);
        } // end for i
    } // end if GetDebug()
    //+++++++++++++++++++++++++
    sA1->SetVertex(0,xo[0],yo[0]);
    sA1->SetVertex(1,xo[1],yo[1]);
    sA1->SetVertex(2,xi[1],yi[1]);
    sA1->SetVertex(3,xi[0],yi[0]);
    //
    sA2->SetVertex(0,xo[1],yo[1]);
    sA2->SetVertex(1,xo[2],yo[2]);
    sA2->SetVertex(2,xi[2],yi[2]);
    sA2->SetVertex(3,xi[1],yi[1]);
    //
    sA3->SetVertex(0,xo[5],yo[5]);
    sA3->SetVertex(1,xo[6],yo[6]);
    sA3->SetVertex(2,xi[6],yi[6]);
    sA3->SetVertex(3,xi[5],yi[5]);
    //--------------------------
    sB1->SetVertex(0,xbo[0],ybo[0]);
    sB1->SetVertex(1,xbo[1],ybo[1]);
    sB1->SetVertex(2,xbi[1],ybi[1]);
    sB1->SetVertex(3,xbi[0],ybi[0]);
    //
    sB2->SetVertex(0,xbo[1],ybo[1]);
    sB2->SetVertex(1,xbo[2],ybo[2]);
    sB2->SetVertex(2,xbi[2],ybi[2]);
    sB2->SetVertex(3,xbi[1],ybi[1]);
    //
    sB3->SetVertex(0,xbo[5],ybo[5]);
    sB3->SetVertex(1,xbo[6],ybo[6]);
    sB3->SetVertex(2,xbi[6],ybi[6]);
    sB3->SetVertex(3,xbi[5],ybi[5]);
    //--------------------------
    sC1->SetVertex(0,xco[0],yco[0]);
    sC1->SetVertex(1,xco[1],yco[1]);
    sC1->SetVertex(2,xci[1],yci[1]);
    sC1->SetVertex(3,xci[0],yci[0]);
    //
    sC2->SetVertex(0,xco[1],yco[1]);
    sC2->SetVertex(1,xco[2],yco[2]);
    sC2->SetVertex(2,xci[2],yci[2]);
    sC2->SetVertex(3,xci[1],yci[1]);
    //
    sC3->SetVertex(0,xco[5],yco[5]);
    sC3->SetVertex(1,xco[6],yco[6]);
    sC3->SetVertex(2,xci[6],yci[6]);
    sC3->SetVertex(3,xci[5],yci[5]);
    // Defining the hole, filled with air
    Double_t lp1,lc1,x,y,x7[3],y7[3];
    lp1 = (xo[0]-xi[0])/(yo[0]-yi[0]);
    lc1 = xo[0]+0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(xo[0]-xi[0])+
                                            SQ(yo[0]-yi[0]))/(xo[0]-xi[0]);
    y = ktscRoutA-2.*ktscarbonFiberThA;
    x = lp1*(y-yo[0])+lc1;
    sAh1->SetVertex(0,x,y);
    sBh1->SetVertex(0,x,y);
    sCh1->SetVertex(4,x,y);
    y = ktscRinA+ktscarbonFiberThA;
    x = lp1*(y-yo[0])+lc1;
    sAh1->SetVertex(3,x,y);
    sBh1->SetVertex(3,x,y);
    x7[0] = x; y7[0] = y; // vortexing done after last point
    //sCh1->SetVertex(7,x,y);
    lp1 = (xo[1]-xi[1])/(yo[1]-yi[1]);
    lc1 = xo[1]-0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(xo[1]-xi[1])+
                                            SQ(yo[1]-yi[1]))/(xo[1]-xi[1]);
    y = ktscRoutA-2.*ktscarbonFiberThA;
    x = lp1*(y-yo[1])+lc1;
    sAh1->SetVertex(1,x,y);
    sBh1->SetVertex(1,x,y);
    sCh1->SetVertex(5,x,y);
    y = ktscRinA+ktscarbonFiberThA;
    x = lp1*(y-yo[1])+lc1;
    sAh1->SetVertex(2,x,y);
    sBh1->SetVertex(2,x,y);
    sCh1->SetVertex(6,x,y);
    //
    // The easist way to get the points for the hole in volume sA2 is to
    // rotate it to the Y axis where the y coordinates are easier to know
    // and then rotate it back.
    Double_t xp,yp,xa,ya,xb,yb;
    th = 0.5*ktscAngle;
    xa = CosD(th)*xo[1]-SinD(th)*yo[1];
    ya = SinD(th)*xo[1]+CosD(th)*yo[1];
    xb = CosD(th)*xi[1]-SinD(th)*yi[1];
    yb = SinD(th)*xi[1]+CosD(th)*yi[1];
    lp1 = (xa-xb)/(ya-yb);
    lc1 = xa+0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(xa-xb)+SQ(ya-yb))/(xa-xb);
    y = ya-ktscarbonFiberThA;
    x = lp1*(y-ya)+lc1;
    xp = CosD(-th)*x-SinD(-th)*y;
    yp = SinD(-th)*x+CosD(-th)*y;
    sAh2->SetVertex(0,xp,yp);
    sBh2->SetVertex(0,xp,yp);
    sCh2->SetVertex(4,xp,yp);
    y = yb+2.0*ktscarbonFiberThA;
    x = lp1*(y-ya)+lc1;
    xp = CosD(-th)*x-SinD(-th)*y;
    yp = SinD(-th)*x+CosD(-th)*y;
    sAh2->SetVertex(3,xp,yp);
    sBh2->SetVertex(3,xp,yp);
    x7[1] = x; y7[1] = y; // vortexing done after last point
    //sCh2->SetVertex(7,xp,yp);
    xa = CosD(th)*xo[2]-SinD(th)*yo[2];
    ya = SinD(th)*xo[2]+CosD(th)*yo[2];
    xb = CosD(th)*xi[2]-SinD(th)*yi[2];
    yb = SinD(th)*xi[2]+CosD(th)*yi[2];
    lp1 = (xa-xb)/(ya-yb);
    lc1 = xa-0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(xa-xb)+SQ(ya-yb))/(xa-xb);
    y = ya-ktscarbonFiberThA;
    x = lp1*(y-ya)+lc1;
    xp = CosD(-th)*x-SinD(-th)*y;
    yp = SinD(-th)*x+CosD(-th)*y;
    sAh2->SetVertex(1,xp,yp);
    sBh2->SetVertex(1,xp,yp);
    sCh2->SetVertex(5,xp,yp);
    y = yb+2.0*ktscarbonFiberThA;
    x = lp1*(y-ya)+lc1;
    xp = CosD(-th)*x-SinD(-th)*y;
    yp = SinD(-th)*x+CosD(-th)*y;
    sAh2->SetVertex(2,xp,yp);
    sBh2->SetVertex(2,xp,yp);
    sCh2->SetVertex(6,xp,yp);
    //
    lp1 = (yo[5]-yi[5])/(xo[5]-xi[5]);
    lc1 = yo[5]+0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(yo[5]-yi[5])+
                                            SQ(xo[5]-xi[5]))/(yo[5]-yi[5]);
    x = xo[5]-ktscarbonFiberThA;
    y = lp1*(x-xo[5])+lc1;
    sAh3->SetVertex(0,x,y);
    sBh3->SetVertex(0,x,y);
    sCh3->SetVertex(4,x,y);
    x = xi[5]+2.0*ktscarbonFiberThA;
    y = lp1*(x-xo[5])+lc1;
    sAh3->SetVertex(3,x,y);
    sBh3->SetVertex(3,x,y);
    x7[2] = x; y7[2] = y; // vortexing done after last point
    //sCh3->SetVertex(7,x,y);
    y = 2.0*ktscarbonFiberThA;
    x = xo[5]-ktscarbonFiberThA;
    sAh3->SetVertex(1,x,y);
    sBh3->SetVertex(1,x,y);
    sCh3->SetVertex(5,x,y);
    y = 2.0*ktscarbonFiberThA;
    x = xi[5]+2.0*ktscarbonFiberThA;
    sAh3->SetVertex(2,x,y);
    sBh3->SetVertex(2,x,y);
    sCh3->SetVertex(6,x,y);
    //
    for(i=0;i<4;i++){ // define points at +dz
     sA1->SetVertex(i+4,(sA1->GetVertices())[2*i],(sA1->GetVertices())[1+2*i]);
     sA2->SetVertex(i+4,(sA2->GetVertices())[2*i],(sA2->GetVertices())[1+2*i]);
     sA3->SetVertex(i+4,(sA3->GetVertices())[2*i],(sA3->GetVertices())[1+2*i]);
     //
     sB1->SetVertex(i+4,(sB1->GetVertices())[2*i],(sB1->GetVertices())[1+2*i]);
     sB2->SetVertex(i+4,(sB2->GetVertices())[2*i],(sB2->GetVertices())[1+2*i]);
     sB3->SetVertex(i+4,(sB3->GetVertices())[2*i],(sB3->GetVertices())[1+2*i]);
     // C's are a cone which must match up with B's.
     sC1->SetVertex(i+4,(sB1->GetVertices())[2*i],(sB1->GetVertices())[1+2*i]);
     sC2->SetVertex(i+4,(sB2->GetVertices())[2*i],(sB2->GetVertices())[1+2*i]);
     sC3->SetVertex(i+4,(sB3->GetVertices())[2*i],(sB3->GetVertices())[1+2*i]);
     //
     sAh1->SetVertex(i+4,(sAh1->GetVertices())[2*i],
                     (sAh1->GetVertices())[1+2*i]);
     sAh2->SetVertex(i+4,(sAh2->GetVertices())[2*i],
                     (sAh2->GetVertices())[1+2*i]);
     sAh3->SetVertex(i+4,(sAh3->GetVertices())[2*i],
                     (sAh3->GetVertices())[1+2*i]);
     //
     sBh1->SetVertex(i+4,(sBh1->GetVertices())[2*i],
                     (sBh1->GetVertices())[1+2*i]);
     sBh2->SetVertex(i+4,(sBh2->GetVertices())[2*i],
                     (sBh2->GetVertices())[1+2*i]);
     sBh3->SetVertex(i+4,(sBh3->GetVertices())[2*i],
                     (sBh3->GetVertices())[1+2*i]);
    } // end for
    //
    lp1 = (xco[0]-xci[0])/(yco[0]-yci[0]);
    lc1 = xco[0]+0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(xco[0]-xci[0])+
                                           SQ(yco[0]-yci[0]))/(xco[0]-xci[0]);
    y = ktscRoutC-2.*ktscarbonFiberThA;
    x = lp1*(y-yco[0])+lc1;
    sCh1->SetVertex(0,x,y);
    y = ktscRinC+ktscarbonFiberThA;
    x = lp1*(y-yci[0])+lc1;
    sCh1->SetVertex(2,x,y);
    lp1 = (xco[1]-xci[1])/(yco[1]-yci[1]);
    lc1 = xco[1]-0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(xco[1]-xci[1])+
                                           SQ(yco[1]-yci[1]))/(xco[1]-xci[1]);
    y = ktscRoutC-2.*ktscarbonFiberThA;
    x = lp1*(y-yco[1])+lc1;
    sCh1->SetVertex(1,x,y);
    y = ktscRinC+ktscarbonFiberThA;
    x = lp1*(y-yci[1])+lc1;
    sCh1->SetVertex(3,x,y);
    //
    th = 0.5*ktscAngle;
    xa = CosD(th)*xco[1]-SinD(th)*yco[1];
    ya = SinD(th)*xco[1]+CosD(th)*yco[1];
    xb = CosD(th)*xci[1]-SinD(th)*yci[1];
    yb = SinD(th)*xci[1]+CosD(th)*yci[1];
    lp1 = (xa-xb)/(ya-yb);
    lc1 = xa+0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(xa-xb)+SQ(ya-yb))/(xa-xb);
    y = ya-ktscarbonFiberThA;
    x = lp1*(y-ya)+lc1;
    xp = CosD(-th)*x-SinD(-th)*y;
    yp = SinD(-th)*x+CosD(-th)*y;
    yp = ya-ktscarbonFiberThA;
    xp = lp1*(y-ya)+lc1;
    sCh2->SetVertex(0,xp,yp);
    y = yb+2.0*ktscarbonFiberThA;
    x = lp1*(y-ya)+lc1;
    xp = CosD(-th)*x-SinD(-th)*y;
    yp = SinD(-th)*x+CosD(-th)*y;
    sCh2->SetVertex(2,xp,yp);
    xa = CosD(th)*xco[2]-SinD(th)*yco[2];
    ya = SinD(th)*xco[2]+CosD(th)*yco[2];
    xb = CosD(th)*xci[2]-SinD(th)*yci[2];
    yb = SinD(th)*xci[2]+CosD(th)*yci[2];
    lp1 = (xa-xb)/(ya-yb);
    lc1 = xa-0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(xa-xb)+SQ(ya-yb))/(xa-xb);
    y = ya-ktscarbonFiberThA;
    x = lp1*(y-ya)+lc1;
    xp = CosD(-th)*x-SinD(-th)*y;
    yp = SinD(-th)*x+CosD(-th)*y;
    sCh2->SetVertex(1,xp,yp);
    y = yb+2.0*ktscarbonFiberThA;
    x = lp1*(y-ya)+lc1;
    xp = CosD(-th)*x-SinD(-th)*y;
    yp = SinD(-th)*x+CosD(-th)*y;
    sCh2->SetVertex(3,xp,yp);
    //
    lp1 = (yco[5]-yci[5])/(xco[5]-xci[5]);
    lc1 = yco[5]+0.5*ktscarbonFiberThA*TMath::Sqrt(SQ(yco[5]-yci[5])+
                                          SQ(xco[5]-xci[5]))/(yco[5]-yci[5]);
    x = xco[5]-ktscarbonFiberThA;
    y = lp1*(x-xco[5])+lc1;
    sCh3->SetVertex(0,x,y);
    x = xci[5]+2.0*ktscarbonFiberThA;
    y = lp1*(x-xci[5])+lc1;
    sCh3->SetVertex(2,x,y);
    y = 2.0*ktscarbonFiberThA;
    x = xco[5]-ktscarbonFiberThA;
    sCh3->SetVertex(1,x,y);
    y = 2.0*ktscarbonFiberThA;
    x = xci[5]+2.0*ktscarbonFiberThA;
    sCh3->SetVertex(3,x,y);
    sCh1->SetVertex(7,x7[0],y7[0]); // 7th point most be done last ???
    sCh2->SetVertex(7,x7[1],y7[1]); // 7th point most be done last ???
    sCh3->SetVertex(7,x7[2],y7[2]); // 7th point most be done last ???
    //
    // Define Minimal volume to inclose this SPD Thermal Sheald.
    sM1 = new TGeoPcon("ITSspdShealdVV",0.0,360.0,9);
    sM1->Z(0)    = 0.5*ktscLengthA+ktscLengthB;
    sM1->Rmin(0) = ktscRinB;
    x = sB1->GetVertices()[0]; // [0][0]
    y = sB1->GetVertices()[1]; // [0][1]
    sM1->Rmax(0) = TMath::Sqrt(x*x+y*y);
    sM1->Z(1)    = sM1->GetZ(0)-ktscLengthB;
    sM1->Rmin(1) = sM1->GetRmin(0);
    sM1->Rmax(1) = sM1->GetRmax(0);
    sM1->Z(2)    = sM1->GetZ(1);
    sM1->Rmin(2) = ktscRinA;
    x = sA1->GetVertices()[0]; // [0]0]
    y = sA1->GetVertices()[1]; // [0][1]
    sM1->Rmax(2) = TMath::Sqrt(x*x+y*y);
    sM1->Z(3)    = -(sM1->GetZ(0)-ktscLengthB);
    sM1->Rmin(3) = sM1->GetRmin(2);
    sM1->Rmax(3) = sM1->GetRmax(2);
    sM1->Z(4)    = sM1->GetZ(3);
    sM1->Rmin(4) = sM1->GetRmin(1);
    sM1->Rmax(4) = sM1->GetRmax(1);
    sM1->Z(5)    = -(sM1->GetZ(0));
    sM1->Rmin(5) = sM1->GetRmin(0);
    sM1->Rmax(5) = sM1->GetRmax(0);
    sM1->Z(6)    = sM1->GetZ(5) - ktscLengthC;
    sM1->Rmin(6) = ktscRinC;
    x = sC1->GetVertices()[0]; // [0][0]
    y = sC1->GetVertices()[1]; // [0][1]
    sM1->Rmax(6) = TMath::Sqrt(x*x+y*y);
    sM1->Z(7)    = sM1->GetZ(6);
    sM1->Rmin(7) = sD->GetRmin();
    sM1->Rmax(7) = sD->GetRmax();
    sM1->Z(8)    = sM1->Z(7) - ktscLengthD;
    sM1->Rmin(8) = sM1->GetRmin(7);
    sM1->Rmax(8) = sM1->GetRmax(7);
    sM2 = new TGeoTubeSeg("ITSspdShealdWingVV",
                          sM1->GetRmax(8),sDw->GetRmax(),sDw->GetDz(),
                          sDw->GetPhi1(),sDw->GetPhi2());
    //
    x = 0.5*(sM1->GetZ(8) + sM1->GetZ(7));
    tranITSspdShealdVVt0 = new TGeoTranslation("ITSspdShealdVVt0",0.0,0.0,x);
    tranITSspdShealdVVt0->RegisterYourself();
    TGeoRotation rotz90("",0.0,0.0,90.0); // never registered.
    rotITSspdShealdVVt1 = new TGeoCombiTrans(*tranITSspdShealdVVt0,rotz90);
    rotITSspdShealdVVt1->SetName("ITSspdShealdVVt1");
    rotITSspdShealdVVt1->RegisterYourself();
    TGeoRotation rotz180("",0.0,0.0,180.0); // never registered
    rotITSspdShealdVVt2 = new TGeoCombiTrans(*tranITSspdShealdVVt0,rotz180);
    rotITSspdShealdVVt2->SetName("ITSspdShealdVVt2");
    rotITSspdShealdVVt2->RegisterYourself();
    TGeoRotation rotz270("",0.0,0.0,270.0); // never registered
    rotITSspdShealdVVt3 = new TGeoCombiTrans(*tranITSspdShealdVVt0,rotz270);
    rotITSspdShealdVVt3->SetName("ITSspdShealdVVt3");
    rotITSspdShealdVVt3->RegisterYourself();
    sM = new TGeoCompositeShape("ITS SPD Thermal sheald volume",
                                "(((ITSspdShealdVV+"
                                "ITSspdShealdWingVV:ITSspdShealdVVt0)+"
                                "ITSspdShealdWingVV:ITSspdShealdVVt1)+"
                                "ITSspdShealdWingVV:ITSspdShealdVVt2)+"
                                "ITSspdShealdWingVV:ITSspdShealdVVt3");
    //
    if(GetDebug()){
        tranITSspdShealdVVt0->Print();
        rotITSspdShealdVVt1->Print();
        rotITSspdShealdVVt2->Print();
        rotITSspdShealdVVt3->Print();
        sD->InspectShape();
        sDs->InspectShape();
        sDw->InspectShape();
        sDws->InspectShape();
        sA1->InspectShape();
        sAh1->InspectShape();
        sA2->InspectShape();
        sAh2->InspectShape();
        sA3->InspectShape();
        sAh3->InspectShape();
        sB1->InspectShape();
        sBh1->InspectShape();
        sB2->InspectShape();
        sBh2->InspectShape();
        sB3->InspectShape();
        sBh3->InspectShape();
        sC1->InspectShape();
        sCh1->InspectShape();
        sC2->InspectShape();
        sCh2->InspectShape();
        sC3->InspectShape();
        sCh3->InspectShape();
        sM1->InspectShape();
        sM2->InspectShape();
        sM->InspectShape();
    } // end if GetDebug
    //
    TGeoManager *mgr = gGeoManager;
    medSPDcf = mgr->GetMedium("ITSspdCarbonFiber");
    medSPDfs = mgr->GetMedium("ITSspdStaselite4411w");
    medSPDfo = mgr->GetMedium("ITSspdRohacell50A");
    medSPDss = mgr->GetMedium("ITSspdStainlessSteal");
    medSPDair= mgr->GetMedium("ITSspdAir");
    TGeoVolume *vA1,*vA2,*vA3,*vAh1,*vAh2,*vAh3;
    TGeoVolume *vB1,*vB2,*vB3,*vBh1,*vBh2,*vBh3;
    TGeoVolume *vC1,*vC2,*vC3,*vCh1,*vCh2,*vCh3;
    TGeoVolume *vD,*vDs,*vDw,*vDws,*vM;
    vM = new TGeoVolume("ITSspdThermalSheald",sM,medSPDair);
    vM->SetVisibility(kTRUE);
    vM->SetLineColor(7); // light Blue
    vM->SetLineWidth(1);
    vM->SetFillColor(vM->GetLineColor());
    vM->SetFillStyle(4090); // 90% transparent
    moth->AddNode(vM,1,0); ///////////////////// Virtual Volume ////////
    vA1 = new TGeoVolume("ITSspdCentCylA1CF",sA1,medSPDcf);
    vA1->SetVisibility(kTRUE);
    vA1->SetLineColor(4);
    vA1->SetLineWidth(1);
    vA2 = new TGeoVolume("ITSspdCentCylA2CF",sA2,medSPDcf);
    vA2->SetVisibility(kTRUE);
    vA2->SetLineColor(4);
    vA2->SetLineWidth(1);
    vA3 = new TGeoVolume("ITSspdCentCylA3CF",sA3,medSPDcf);
    vA3->SetVisibility(kTRUE);
    vA3->SetLineColor(4);
    vA3->SetLineWidth(1);
    vB1 = new TGeoVolume("ITSspdCentCylB1CF",sB1,medSPDcf);
    vB1->SetVisibility(kTRUE);
    vB1->SetLineColor(4);
    vB1->SetLineWidth(1);
    vB2 = new TGeoVolume("ITSspdCentCylB2CF",sB2,medSPDcf);
    vB2->SetVisibility(kTRUE);
    vB2->SetLineColor(4);
    vB2->SetLineWidth(1);
    vB3 = new TGeoVolume("ITSspdCentCylB3CF",sB3,medSPDcf);
    vB3->SetVisibility(kTRUE);
    vB3->SetLineColor(4);
    vB3->SetLineWidth(1);
    vC1 = new TGeoVolume("ITSspdCentCylC1CF",sC1,medSPDcf);
    vC1->SetVisibility(kTRUE);
    vC1->SetLineColor(4);
    vC1->SetLineWidth(1);
    vC2 = new TGeoVolume("ITSspdCentCylC2CF",sC2,medSPDcf);
    vC2->SetVisibility(kTRUE);
    vC2->SetLineColor(4);
    vC2->SetLineWidth(1);
    vC3 = new TGeoVolume("ITSspdCentCylC3CF",sC3,medSPDcf);
    vC3->SetVisibility(kTRUE);
    vC3->SetLineColor(4);
    vC3->SetLineWidth(1);
    vAh1 = new TGeoVolume("ITSspdCentCylA1AirA",sAh1,medSPDair);
    vAh1->SetVisibility(kTRUE);
    vAh1->SetLineColor(5); // Yellow
    vAh1->SetFillColor(vAh1->GetLineColor());
    vAh1->SetFillStyle(4090); // 90% transparent
    vAh2 = new TGeoVolume("ITSspdCentCylA2AirA",sAh2,medSPDair);
    vAh2->SetVisibility(kTRUE);
    vAh2->SetLineColor(5); // Yellow
    vAh2->SetFillColor(vAh2->GetLineColor());
    vAh2->SetFillStyle(4090); // 90% transparent
    vAh3 = new TGeoVolume("ITSspdCentCylA3AirA",sAh3,medSPDair);
    vAh3->SetVisibility(kTRUE);
    vAh3->SetLineColor(5); // Yellow
    vAh3->SetFillColor(vAh3->GetLineColor());
    vAh3->SetFillStyle(4090); // 90% transparent
    vBh1 = new TGeoVolume("ITSspdCentCylA1AirB",sBh1,medSPDair);
    vBh1->SetVisibility(kTRUE);
    vBh1->SetLineColor(5); // Yellow
    vBh1->SetFillColor(vBh1->GetLineColor());
    vBh1->SetFillStyle(4090); // 90% transparent
    vBh2 = new TGeoVolume("ITSspdCentCylA2AirB",sBh2,medSPDair);
    vBh2->SetVisibility(kTRUE);
    vBh2->SetLineColor(5); // Yellow
    vBh2->SetFillColor(vBh2->GetLineColor());
    vBh2->SetFillStyle(4090); // 90% transparent
    vBh3 = new TGeoVolume("ITSspdCentCylA3AirB",sBh3,medSPDair);
    vBh3->SetVisibility(kTRUE);
    vBh3->SetLineColor(5); // Yellow
    vBh3->SetFillColor(vBh3->GetLineColor());
    vBh3->SetFillStyle(4090); // 90% transparent
    vCh1 = new TGeoVolume("ITSspdCentCylA1AirC",sCh1,medSPDair);
    vCh1->SetVisibility(kTRUE);
    vCh1->SetLineColor(5); // Yellow
    vCh1->SetFillColor(vCh1->GetLineColor());
    vCh1->SetFillStyle(4090); // 90% transparent
    vCh2 = new TGeoVolume("ITSspdCentCylA2AirC",sCh2,medSPDair);
    vCh2->SetVisibility(kTRUE);
    vCh2->SetLineColor(5); // Yellow
    vCh2->SetFillColor(vCh2->GetLineColor());
    vCh2->SetFillStyle(4090); // 90% transparent
    vCh3 = new TGeoVolume("ITSspdCentCylA3AirC",sCh3,medSPDair);
    vCh3->SetVisibility(kTRUE);
    vCh3->SetLineColor(5); // Yellow
    vCh3->SetFillColor(vCh3->GetLineColor());
    vCh3->SetFillStyle(4090); // 90% transparent
    vD = new TGeoVolume("ITSspdCentCylA1CD",sD,medSPDcf);
    vD->SetVisibility(kTRUE);
    vD->SetLineColor(4);
    vD->SetLineWidth(1);
    vDw = new TGeoVolume("ITSspdCentCylA1CDw",sDw,medSPDcf);
    vDw->SetVisibility(kTRUE);
    vDw->SetLineColor(4);
    vDw->SetLineWidth(1);
    vDs = new TGeoVolume("ITSspdCentCylA1Dfill",sDs,medSPDfs);
    vDs->SetVisibility(kTRUE);
    vDs->SetLineColor(3); // Green
    vDs->SetFillColor(vDs->GetLineColor());
    vDs->SetFillStyle(4010); // 10% transparent
    vDws = new TGeoVolume("ITSspdCentCylA1DwingFill",sDws,medSPDfs);
    vDws->SetVisibility(kTRUE);
    vDws->SetLineColor(3); // Green
    vDws->SetFillColor(vDws->GetLineColor());
    vDws->SetFillStyle(4010); // 10% transparent
    //
    vA1->AddNode(vAh1,1,0);
    vA2->AddNode(vAh2,1,0);
    vA3->AddNode(vAh3,1,0);
    vB1->AddNode(vBh1,1,0);
    vB2->AddNode(vBh2,1,0);
    vB3->AddNode(vBh3,1,0);
    vC1->AddNode(vCh1,1,0);
    vC2->AddNode(vCh2,1,0);
    vC3->AddNode(vCh3,1,0);
    vD ->AddNode(vDs ,1,0);
    vDw->AddNode(vDws,1,0);
    //
    vM->AddNode(vA1,1,0);
    vM->AddNode(vA2,1,0);
    vM->AddNode(vA3,1,0);
    tranb  = new TGeoTranslation("",0.0,0.0,0.5*(ktscLengthA+ktscLengthB));
    tranbm = new TGeoTranslation("",0.0,0.0,0.5*(-ktscLengthA-ktscLengthB));
    vM->AddNode(vB1,1,tranb);
    vM->AddNode(vB2,1,tranb);
    vM->AddNode(vB3,1,tranb);
    vM->AddNode(vB1,2,tranbm);
    vM->AddNode(vB2,2,tranbm);
    vM->AddNode(vB3,2,tranbm);
    // Muon side (rsB26) is at -Z.
    tranc = new TGeoTranslation("",0.0,0.0,
                                0.5*(-ktscLengthA-ktscLengthB-ktscLengthC));
    vM->AddNode(vC1,1,tranc);
    vM->AddNode(vC2,1,tranc);
    vM->AddNode(vC3,1,tranc);
    vM->AddNode(vD,1,tranITSspdShealdVVt0);
    vM->AddNode(vDw,1,tranITSspdShealdVVt0);
    vM->AddNode(vDw,2,rotITSspdShealdVVt1);
    vM->AddNode(vDw,3,rotITSspdShealdVVt2);
    vM->AddNode(vDw,4,rotITSspdShealdVVt3);
    k=2;
    for(i=1;i<10;i++) {
        th = ((Double_t)i)*ktscAngle*fgkDegree;
        rot = new TGeoRotation("",0.0,0.0,th);
        vM->AddNode(vA1,i+1,rot);
        vM->AddNode(vB1,i+2,new TGeoCombiTrans(*tranb,*rot));
        vM->AddNode(vB1,i+12,new TGeoCombiTrans(*tranbm,*rot));
        vM->AddNode(vC1,i+1,new TGeoCombiTrans(*tranc,*rot));
        if(i!=0||i!=2||i!=7){
            vM->AddNode(vA2,k++,rot);
            vM->AddNode(vB2,k++,new TGeoCombiTrans(*tranb,*rot));
            vM->AddNode(vB2,k++,new TGeoCombiTrans(*tranbm,*rot));
            vM->AddNode(vC2,k++,new TGeoCombiTrans(*tranc,*rot));
        } // end if
        if(i==5) {
            vM->AddNode(vA3,2,rot);
            vM->AddNode(vB3,3,new TGeoCombiTrans(*tranb,*rot));
            vM->AddNode(vB3,4,new TGeoCombiTrans(*tranbm,*rot));
            vM->AddNode(vC3,2,new TGeoCombiTrans(*tranc,*rot));
        } // end if
    } // end for i
    rot = new TGeoRotation("",180.,0.0,0.0);
    vM->AddNode(vA3,3,rot);
    vM->AddNode(vB3,5,new TGeoCombiTrans(*tranb,*rot));
    vM->AddNode(vB3,6,new TGeoCombiTrans(*tranbm,*rot));
    vM->AddNode(vC3,3,new TGeoCombiTrans(*tranc,*rot));
    rot = new TGeoRotation("",180.,0.0,180.0);
    vM->AddNode(vA3,4,rot);
    vM->AddNode(vB3,7,new TGeoCombiTrans(*tranb,*rot));
    vM->AddNode(vB3,8,new TGeoCombiTrans(*tranbm,*rot));
    vM->AddNode(vC3,4,new TGeoCombiTrans(*tranc,*rot));
    if(GetDebug()){
        vA1->PrintNodes();
        vAh1->PrintNodes();
        vA2->PrintNodes();
        vAh2->PrintNodes();
        vA3->PrintNodes();
        vAh3->PrintNodes();
        vB1->PrintNodes();
        vBh1->PrintNodes();
        vB2->PrintNodes();
        vBh2->PrintNodes();
        vB3->PrintNodes();
        vBh3->PrintNodes();
        vC1->PrintNodes();
        vCh1->PrintNodes();
        vC2->PrintNodes();
        vCh2->PrintNodes();
        vC3->PrintNodes();
        vCh3->PrintNodes();
        vD->PrintNodes();
        vDs->PrintNodes();
        vDw->PrintNodes();
        vDws->PrintNodes();
        vM->PrintNodes();
    } // end if
}
//______________________________________________________________________
void AliITSv11GeometrySupport::SDDCone(TGeoVolume *moth){
    // Define the detail SDD support cone geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
    //
    // From Cilindro Centrale - Lavorazioni, ALR 0816/1 04/08/03 File
    // name SDD/Cilindro.hpgl
    const Double_t ktsLength       = 790.0*fgkmm; // Thermal Sheeld length
    const Double_t ktsInsertoLength= 15.0*fgkmm;    // ????
    const Double_t ktsOuterR       = 0.5*(220.+10.)*fgkmm; // ????
    const Double_t ktsInnerR       = 0.5*(220.-10.)*fgkmm; // ????
    const Double_t ktscarbonFiberth= 0.02*fgkmm;     // ????
    const Double_t ktsBoltDiameter = 6.0*fgkmm; // M6 screw
    const Double_t ktsBoltDepth    = 6.0*fgkmm; // in volume sC
    const Double_t ktsBoltRadius   = 0.5*220.*fgkmm; // Radius in volume sC
    const Double_t ktsBoltAngle0   = 0.0*fgkDegree; // Angle in volume sC
    const Double_t ktsBoltdAngle   = 30.0*fgkDegree; // Angle in Volume sC
    Double_t x,y,z,t,t0,rmin,rmax;
    Int_t i,n;
    TGeoTube *sA,*sB,*sC,*sD;
    TGeoTranslation *tran;
    TGeoRotation *rot;
    TGeoCombiTrans *rotran;
    TGeoMedium *medSDDcf,*medSDDfs,*medSDDfo,*medSDDss;

    sA = new TGeoTube("ITS SDD Central Cylinder",ktsInnerR,ktsOuterR,
                     0.5*ktsLength);
    sB = new TGeoTube("ITS SDD CC Foam",ktsInnerR+ktscarbonFiberth,
                    ktsOuterR-ktscarbonFiberth,
                    0.5*(ktsLength-2.0*ktsInsertoLength));
    sC = new TGeoTube("ITS SDD CC Inserto",ktsInnerR+ktscarbonFiberth,
                    ktsOuterR-ktscarbonFiberth,0.5*ktsLength);
    sD = new TGeoTube("ITS SDD CC M6 bolt end",0.0,0.5*ktsBoltDiameter,
                    0.5*ktsBoltDepth);
    if(GetDebug()){
        sA->InspectShape();
        sB->InspectShape();
        sC->InspectShape();
        sD->InspectShape();
    } // end if GetDebug
    //
    TGeoManager *mgr = gGeoManager;
    medSDDcf = mgr->GetMedium("ITSssdCarbonFiber");
    medSDDfs = mgr->GetMedium("ITSssdStaselite4411w");
    medSDDfo = mgr->GetMedium("ITSssdRohacell50A");
    medSDDss = mgr->GetMedium("ITSssdStainlessSteal");
    TGeoVolume *vA,*vB,*vC,*vD;
    vA = new TGeoVolume("ITSsddCentCylCF",sA,medSDDcf);
    vA->SetVisibility(kTRUE);
    vA->SetLineColor(4);
    vA->SetLineWidth(1);
    vA->SetFillColor(vA->GetLineColor());
    vA->SetFillStyle(4030); // 30% transparent
    vB = new TGeoVolume("ITSsddCentCylF",sB,medSDDfo);
    vB->SetVisibility(kTRUE);
    vB->SetLineColor(3);
    vB->SetLineWidth(1);
    vB->SetFillColor(vB->GetLineColor());
    vB->SetFillStyle(4050); // 50% transparent
    vC = new TGeoVolume("ITSsddCentCylSt",sC,medSDDfs);
    vC->SetVisibility(kTRUE);
    vC->SetLineColor(2);
    vC->SetLineWidth(1);
    vC->SetFillColor(vC->GetLineColor());
    vC->SetFillStyle(4050); // 50% transparent
    vD = new TGeoVolume("ITSsddCentCylSS",sD,medSDDss);
    vD->SetVisibility(kTRUE);
    vD->SetLineColor(1);
    vD->SetLineWidth(1);
    vD->SetFillColor(vD->GetLineColor());
    vD->SetFillStyle(4050); // 50% transparent
    //
    moth->AddNode(vA,1,0);
    vA->AddNode(vC,1,0);
    vC->AddNode(vB,1,0);
    n = (Int_t)((360.*fgkDegree)/ktsBoltdAngle);
    for(i=0;i<n;i++){
        t = ktsBoltAngle0+((Double_t)i)*ktsBoltdAngle;
        x = ktsBoltRadius*CosD(t);
        y = ktsBoltRadius*SinD(t);
        z = 0.5*(ktsLength-ktsBoltDepth);
        tran = new TGeoTranslation("",x,y,z);
        vC->AddNode(vD,i+1,tran);
        tran = new TGeoTranslation("",x,y,-z);
        vC->AddNode(vD,i+n+1,tran);
    } // end for i
    if(GetDebug()){
        vA->PrintNodes();
        vB->PrintNodes();
        vC->PrintNodes();
        vD->PrintNodes();
    } // end if
    // SDD Suport Cone
    //
    //
    const Double_t kconThickness    = 10.5*fgkmm;//Thickness Rohacell+car. fib.
    const Double_t kconCthick       = 1.5*fgkmm; // Carbon finber thickness
    const Double_t kconRcurv        = 15.0*fgkmm; // Radius of curvature.
    const Double_t kconTc           = 45.0; // angle of SDD cone [degrees].
    const Double_t kconZouterMilled = 23.0*fgkmm;
    const Double_t kconZcylinder    = 186.0*fgkmm;
    const Double_t kconZ0           = kconZcylinder + 0.5*ktsLength;
    //const Int_t kconNspoaks         = 12;
    //const Int_t kconNmounts         = 4;
    //const Double_t kconDmountAngle  = 9.0; // degrees
    const Double_t kconRoutMax      = 0.5*560.0*fgkmm;
    const Double_t kconRoutMin      = 0.5*539.0*fgkmm;
    // Holes in cone for cables
    const Double_t kconPhiHole1     = 0.0*fgkDegree;
    const Double_t kcondPhiHole1    = 25.0*fgkDegree;
    const Double_t kconRholeMax1    = 0.5*528.*fgkmm;
    const Double_t kconRholeMin1    = 0.5*464.*fgkmm;
    const Double_t kconPhiHole2     = 0.0*fgkDegree;
    const Double_t kcondPhiHole2    = 50.0*fgkDegree;
    const Double_t kconRholeMax2    = 0.5*375.*fgkmm;
    const Double_t kconRholeMin2    = 0.5*280.*fgkmm;
    //
    //const Int_t kconNpostsOut       = 6;
    //const Int_t kconNpostsIn        = 3;
    //const Double_t kconPhi0PostOut  = 0.0; // degree
    //const Double_t kconPhi0PostIn   = 0.0; // degree
    //const Double_t kcondRpostOut    = 16.0*fgkmm;
    //const Double_t kcondRpostIn     = 16.0*fgkmm;
    //const Double_t kconZpostMaxOut  = 116.0*fgkmm;
    //const Double_t kconZpostMaxIn   = 190.0*fgkmm;
    const Double_t kconRinMax       = 0.5*216*fgkmm;
    const Double_t kconRinCylinder  = 0.5*231.0*fgkmm;
    //const Double_t kconRinHole      = 0.5*220.0*fgkmm;
    const Double_t kconRinMin       = 0.5*210.0*fgkmm;
    const Double_t kcondZin         = 15.0*fgkmm; // ???
    const Double_t kSinkconTc       = SinD(kconTc);
    const Double_t kCoskconTc       = CosD(kconTc);
    const Double_t kTankconTc       = TanD(kconTc);
    //
    TGeoPcon *sE,*sF,*sG,*sH,*sI,*sJ,*sK;
    TGeoCompositeShape *sL,*sM,*sN;
    //
    Double_t dza = kconThickness/kSinkconTc-
        (kconRoutMax-kconRoutMin)/kTankconTc;
    if(dza<=0){ // The number or order of the points are in error for a proper
     // call to pcons!
     Error("SDDcone","The definition of the points for a call to PCONS is"
           " in error. abort.");
     return;
    } // end if
    sE = new TGeoPcon("ITSsddSuportConeCarbonFiberSurfaceE",0.0,360.0,12);
    sE->Z(0)    = 0.0;
    sE->Rmin(0) = kconRoutMin;
    sE->Rmax(0) = kconRoutMax;
    sE->Z(1)    = kconZouterMilled - dza;
    sE->Rmin(1) = sE->GetRmin(0);
    sE->Rmax(1) = sE->GetRmax(0);
    sE->Z(2)    = kconZouterMilled;
    sE->Rmax(2) = sE->GetRmax(0);
    RadiusOfCurvature(kconRcurv,0.,sE->GetZ(1),sE->GetRmin(1),kconTc,z,rmin);
    sE->Z(3)    = z;
    sE->Rmin(3) = rmin;
    sE->Rmin(2) = RminFrom2Points(sE,3,1,sE->GetZ(2));
    RadiusOfCurvature(kconRcurv,0.,sE->GetZ(2),sE->GetRmax(2),kconTc,z,rmax);
    sE->Z(4)    = z;
    sE->Rmax(4) = rmax;
    sE->Rmin(4) = RminFromZpCone(sE,3,kconTc,sE->GetZ(4),0.0);
    sE->Rmax(3) = RmaxFrom2Points(sE,4,2,sE->GetZ(3));
    sE->Rmin(7) = kconRinMin;
    sE->Rmin(8) = kconRinMin;
    RadiusOfCurvature(kconRcurv,90.0,0.0,kconRinMax,90.0-kconTc,z,rmax);
    sE->Rmax(8) = rmax;
    sE->Z(8)    = ZFromRmaxpCone(sE,4,kconTc,sE->GetRmax(8));
    sE->Z(9)    = kconZcylinder;
    sE->Rmin(9) = kconRinMin;
    sE->Z(10)    = sE->GetZ(9);
    sE->Rmin(10) = kconRinCylinder;
    sE->Rmin(11) = kconRinCylinder;
    sE->Rmax(11) = sE->GetRmin(11);
    rmin         = sE->GetRmin(8);
    RadiusOfCurvature(kconRcurv,90.0-kconTc,sE->GetZ(8),sE->GetRmax(8),90.0,
                      z,rmax);
    rmax = kconRinMax;
    sE->Z(11)    = z+(sE->GetZ(8)-z)*(sE->GetRmax(11)-rmax)/
                                           (sE->GetRmax(8)-rmax);
    sE->Rmax(9) = RmaxFrom2Points(sE,11,8,sE->GetZ(9));
    sE->Rmax(10) = sE->GetRmax(9);
    sE->Z(6)    = z-kcondZin;
    sE->Z(7)    = sE->GetZ(6);
    sE->Rmax(6) = RmaxFromZpCone(sE,4,kconTc,sE->GetZ(6));
    sE->Rmax(7) = sE->GetRmax(6);
    RadiusOfCurvature(kconRcurv,90.,sE->GetZ(6),0.0,90.0-kconTc,z,rmin);
    sE->Z(5)    = z;
    sE->Rmin(5) = RminFromZpCone(sE,3,kconTc,z);
    sE->Rmax(5) = RmaxFromZpCone(sE,4,kconTc,z);
    RadiusOfCurvature(kconRcurv,90.-kconTc,0.0,sE->Rmin(5),90.0,z,rmin);
    sE->Rmin(6) = rmin;
    // Inner Core, Inserto material
    sF = new TGeoPcon("ITSsddSuportConeInsertoStesaliteF",0.,360.0,9);
    sF->Z(0)    = sE->GetZ(0);
    sF->Rmin(0) = sE->GetRmin(0)+kconCthick;
    sF->Rmax(0) = sE->GetRmax(0)-kconCthick;
    sF->Z(1)    = sE->GetZ(1);
    sF->Rmin(1) = sF->GetRmin(0);
    sF->Rmax(1) = sF->GetRmax(0);
    sF->Z(2)    = sE->GetZ(2);
    sF->Rmax(2) = sF->GetRmax(1);
    RadiusOfCurvature(kconRcurv-kconCthick,0.,sF->GetZ(1),sF->GetRmax(1),
                      kconTc,z,rmin);
    sF->Z(3)    = z;
    sF->Rmin(3) = rmin;
    sF->Rmin(2) = RminFrom2Points(sF,3,1,sF->GetZ(2));
    RadiusOfCurvature(kconRcurv+kconCthick,0.,sF->GetZ(2),sF->GetRmax(2),
                      kconTc,z,rmax);
    sF->Z(4)    = z;
    sF->Rmax(4) = rmax;
    sF->Rmin(4) = RmaxFromZpCone(sE,2,kconTc,sF->GetZ(4),
                                                   -kconCthick);
    sF->Rmax(3) = RmaxFrom2Points(sF,4,2,sF->GetZ(3));
    sF->Rmin(7) = sE->GetRmin(7);
    sF->Rmin(8) = sE->GetRmin(8);
    sF->Z(6)    = sE->GetZ(6)+kconCthick;
    sF->Rmin(6) = sE->GetRmin(6);
    sF->Z(7)    = sF->GetZ(6);
    sF->Rmax(8) = sE->GetRmax(8)-kconCthick*kSinkconTc;
    RadiusOfCurvature(kconRcurv+kconCthick,90.0,sF->GetZ(6),sF->GetRmin(6),
                      90.0-kconTc,z,rmin);
    sF->Z(5)    = z;
    sF->Rmin(5) = rmin;
    sF->Rmax(5) = RmaxFromZpCone(sF,4,kconTc,z);
    sF->Rmax(6) = RmaxFromZpCone(sF,4,kconTc,sF->GetZ(6));
    sF->Rmax(7) = sF->GetRmax(6);
    sF->Z(8)    = ZFromRmaxpCone(sF,4,kconTc,sF->GetRmax(8),-kconCthick);
    // Inner Core, Inserto material
    sG = new TGeoPcon("ITSsddSuportConeFoamCoreG",0.0,360.0,4);
    RadiusOfCurvature(kconRcurv+kconCthick,0.0,sF->GetZ(1),sF->GetRmin(1),
                      kconTc,z,rmin);
    sG->Z(0)    = z;
    sG->Rmin(0) = rmin;
    sG->Rmax(0) = sG->GetRmin(0);
    sG->Z(1)    = sG->GetZ(0)+(kconThickness-2.0*kconCthick)/kSinkconTc;;
    sG->Rmin(1) = RminFromZpCone(sF,3,kconTc,sG->GetZ(1));
    sG->Rmax(1) = RmaxFromZpCone(sF,4,kconTc,sG->GetZ(1));
    sG->Z(2)    = sE->GetZ(5)-kconCthick;
    sG->Rmin(2) = RminFromZpCone(sF,3,kconTc,sG->GetZ(2));
    sG->Rmax(2) = RmaxFromZpCone(sF,4,kconTc,sG->GetZ(2));
    sG->Z(3)    = sF->GetZ(5)+(kconThickness-2.0*kconCthick)*kCoskconTc;
    sG->Rmax(3) = RmaxFromZpCone(sF,4,kconTc,sG->GetZ(3));
    sG->Rmin(3) = sG->GetRmax(3);
    //
    sH = new TGeoPcon("ITSsddSuportConeHoleH",kconPhiHole1,kcondPhiHole1,4);
    sH->Rmin(0) = kconRholeMax1;
    sH->Rmax(0) = sH->GetRmin(0);
    sH->Z(0)    = ZFromRminpCone(sE,3,kconTc,sH->GetRmin(0));
    sH->Rmax(1) = sH->GetRmax(0);
    sH->Z(1)    = ZFromRmaxpCone(sE,4,kconTc,sH->GetRmax(1));
    sH->Rmin(1) = RminFromZpCone(sE,3,kconTc,sH->GetZ(1));
    sH->Rmin(2) = kconRholeMin1;
    sH->Z(2)    = ZFromRminpCone(sE,3,kconTc,sH->GetRmin(2));
    sH->Rmax(2) = RmaxFromZpCone(sE,4,kconTc,sH->GetZ(2));
    sH->Rmin(3) = sH->GetRmin(2);
    sH->Rmax(3) = sH->GetRmin(3);
    sH->Z(3)    = ZFromRminpCone(sE,3,kconTc,sH->GetRmin(3));
    //
    x = kconCthick/(0.5*(kconRholeMax1+kconRholeMin1));
    t0 = kconPhiHole1 - x*fgkRadian;
    t  = kcondPhiHole1 + 2.0*x*fgkRadian;
    sI = new TGeoPcon("ITSsddSuportConeHoleI",t0,t,4);
    sI->Rmin(0) = kconRholeMax1+kconCthick;
    sI->Rmax(0) = sI->GetRmin(0);
    sI->Z(0)    = ZFromRminpCone(sF,3,kconTc,sI->GetRmin(0));
    sI->Rmax(1) = sI->GetRmax(0);
    sI->Z(1)    = ZFromRmaxpCone(sF,4,kconTc,sI->GetRmax(1));
    sI->Rmin(1) = RminFromZpCone(sF,3,kconTc,sI->GetZ(1));
    sI->Rmin(2) = kconRholeMin1-kconCthick;
    sI->Z(2)    = ZFromRminpCone(sF,3,kconTc,sI->GetRmin(2));
    sI->Rmax(2) = RmaxFromZpCone(sF,4,kconTc,sI->GetZ(2));
    sI->Rmin(3) = sI->GetRmin(2);
    sI->Rmax(3) = sI->GetRmin(3);
    sI->Z(3)    = ZFromRmaxpCone(sF,4,kconTc,sI->GetRmax(3));
    //
    sJ = new TGeoPcon("ITSsddSuportConeHoleJ",kconPhiHole2,
                                kcondPhiHole2,4);
    sJ->Rmin(0) = kconRholeMax2;
    sJ->Rmax(0) = sJ->GetRmin(0);
    sJ->Z(0)    = ZFromRminpCone(sE,3,kconTc,sJ->GetRmin(0));
    sJ->Rmax(1) = sJ->GetRmax(0);
    sJ->Z(1)    = ZFromRmaxpCone(sE,4,kconTc,sJ->GetRmax(1));
    sJ->Rmin(1) = RminFromZpCone(sE,3,kconTc,sJ->GetZ(1));
    sJ->Rmin(2) = kconRholeMin2;
    sJ->Z(2)    = ZFromRminpCone(sE,3,kconTc,sJ->GetRmin(2));
    sJ->Rmax(2) = RmaxFromZpCone(sE,4,kconTc,sJ->GetZ(2));
    sJ->Rmin(3) = sJ->GetRmin(2);
    sJ->Rmax(3) = sJ->GetRmin(3);
    sJ->Z(3)    = ZFromRmaxpCone(sE,4,kconTc,sJ->GetRmax(3));
    //
    x = kconCthick/(0.5*(kconRholeMax2+kconRholeMin2));
    t0 = kconPhiHole2 - x*fgkRadian;
    t  = kcondPhiHole2 + 2.0*x*fgkRadian;
    sK = new TGeoPcon("ITSsddSuportConeHoleK",t0,t,4);
    sK->Rmin(0) = kconRholeMax2+kconCthick;
    sK->Rmax(0) = sK->GetRmin(0);
    sK->Z(0)    = ZFromRminpCone(sF,3,kconTc,sK->GetRmin(0));
    sK->Rmax(1) = sK->GetRmax(0);
    sK->Z(1)    = ZFromRmaxpCone(sF,4,kconTc,sK->GetRmax(1));
    sK->Rmin(1) = RminFromZpCone(sF,3,kconTc,sK->GetZ(1));
    sK->Rmin(2) = kconRholeMin2-kconCthick;
    sK->Z(2)    = ZFromRminpCone(sF,3,kconTc,sK->GetRmin(2));
    sK->Rmax(2) = RmaxFromZpCone(sF,4,kconTc,sK->GetZ(2));
    sK->Rmin(3) = sK->GetRmin(2);
    sK->Rmax(3) = sK->GetRmin(3);
    sK->Z(3)    = ZFromRmaxpCone(sF,4,kconTc,sK->GetRmax(3));
    //
    rot = new TGeoRotation("ITSsddRotZ30",0.0,0.0,30.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ60",0.0,0.0,60.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ90",0.0,0.0,90.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ120",0.0,0.0,120.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ150",0.0,0.0,150.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ180",0.0,0.0,180.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ210",0.0,0.0,210.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ240",0.0,0.0,240.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ270",0.0,0.0,270.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ300",0.0,0.0,300.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    rot = new TGeoRotation("ITSsddRotZ330",0.0,0.0,330.0);
    rot->RegisterYourself();
    if(GetDebug()) rot->Print();
    sL = new TGeoCompositeShape("ITS SDD Suport Cone","((((((((((((((((("
                                "ITSsddSuportConeCarbonFiberSurfaceE -"
                                "ITSsddSuportConeHoleH)  -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ30) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ60) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ90) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ120) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ150) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ180) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ210) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ240) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ270) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ300) -"
                                "ITSsddSuportConeHoleH:ITSsddRotZ330) -"
                                "ITSsddSuportConeHoleJ)  -"
                                "ITSsddSuportConeHoleJ:ITSsddRotZ60) -"
                                "ITSsddSuportConeHoleJ:ITSsddRotZ120) -"
                                "ITSsddSuportConeHoleJ:ITSsddRotZ180) -"
                               "ITSsddSuportConeHoleJ:ITSsddRotZ240) -"
                                "ITSsddSuportConeHoleJ:ITSsddRotZ300");
    sM = new TGeoCompositeShape("ITS SDD Suport Cone Inserto Stesalite",
                                "((((((((((((((((("
                                "ITSsddSuportConeInsertoStesaliteF -"
                                "ITSsddSuportConeHoleI)  -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ30) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ60) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ90) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ120) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ150) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ180) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ210) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ240) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ270) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ300) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ330) -"
                                "ITSsddSuportConeHoleK)  -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ60) -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ120) -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ180) -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ240) -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ300");
    sN = new TGeoCompositeShape("ITS SDD Suport Cone Foam Core",
                                "((((((((((((((((("
                                "ITSsddSuportConeFoamCoreG -"
                                "ITSsddSuportConeHoleI)  -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ30) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ60) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ90) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ120) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ150) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ180) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ210) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ240) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ270) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ300) -"
                                "ITSsddSuportConeHoleI:ITSsddRotZ330) -"
                                "ITSsddSuportConeHoleK)  -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ60) -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ120) -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ180) -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ240) -"
                                "ITSsddSuportConeHoleK:ITSsddRotZ300");
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(GetDebug()){
        sE->InspectShape();
        sF->InspectShape();
        sG->InspectShape();
        sH->InspectShape();
        sI->InspectShape();
        sJ->InspectShape();
        sK->InspectShape();
        sL->InspectShape();
        sM->InspectShape();
        sN->InspectShape();
    } // end if GetDebug()
    //
    TGeoVolume *vL,*vM,*vN;
    vL = new TGeoVolume("ITSsddConeL",sL,medSDDcf);
    vL->SetVisibility(kTRUE);
    vL->SetLineColor(4);
    vL->SetLineWidth(1);
    vL->SetFillColor(vL->GetLineColor());
    vL->SetFillStyle(4000); // 0% transparent
    vM = new TGeoVolume("ITSsddConeM",sM,medSDDfs);
    vM->SetVisibility(kTRUE);
    vM->SetLineColor(2);
    vM->SetLineWidth(1);
    vM->SetFillColor(vM->GetLineColor());
    vM->SetFillStyle(4010); // 10% transparent
    vN = new TGeoVolume("ITSsddConeN",sN,medSDDfo);
    vN->SetVisibility(kTRUE);
    vN->SetLineColor(7);
    vN->SetLineWidth(1);
    vN->SetFillColor(vN->GetLineColor());
    vN->SetFillStyle(4050); // 50% transparent
    //
    vM->AddNode(vN,1,0);
    vL->AddNode(vM,1,0);
    tran = new TGeoTranslation("",0.0,0.0,-kconZ0);
    moth->AddNode(vL,1,tran);
    rot = new TGeoRotation("",0.0,180.0*fgkDegree,0.0);
    rotran = new TGeoCombiTrans("",0.0,0.0,kconZ0,rot);
    moth->AddNode(vL,2,rotran);
    if(GetDebug()){
        tran->Print();
        rot->Print();
        rotran->Print();
        vL->PrintNodes();
        vM->PrintNodes();
        vN->PrintNodes();
    } // end if
    delete rot;// rot not explicity used in AddNode functions.
}
//______________________________________________________________________
void AliITSv11GeometrySupport::SSDCone(TGeoVolume *moth){
    // Define the detail SSD support cone geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
    //
    Int_t i,j;
    Double_t t,t0,dt,x,y,z,vl[3],vg[3],x0,y0,rmin,rmax;
    TGeoMedium *medSSDcf  = 0; // SSD support cone Carbon Fiber materal number.
    TGeoMedium *medSSDfs  = 0; // SSD support cone inserto stesalite 4411w.
    TGeoMedium *medSSDfo  = 0; // SSD support cone foam, Rohacell 50A.
    TGeoMedium *medSSDss  = 0; // SSD support cone screw material,Stainless
    TGeoMedium *medSSDair = 0; // SSD support cone Air
    TGeoMedium *medSSDal  = 0; // SSD support cone SDD mounting bracket Al
    TGeoManager *mgr = gGeoManager;
    medSSDcf = mgr->GetMedium("ITSssdCarbonFiber");
    medSSDfs = mgr->GetMedium("ITSssdStaselite4411w");
    medSSDfo = mgr->GetMedium("ITSssdRohacell50A");
    medSSDss = mgr->GetMedium("ITSssdStainlessSteal");
    medSSDair= mgr->GetMedium("ITSssdAir");
    medSSDal = mgr->GetMedium("ITSssdAl");
    //
    // SSD Central cylinder/Thermal Sheald.
    const Double_t kcylZlength     = 1140.0*fgkmm; //
    const Double_t kcylZFoamlength = 1020.0*fgkmm; //
    const Double_t kcylROuter      = 0.5*595.0*fgkmm; //
    const Double_t kcylRInner      = 0.5*560.5*fgkmm; //
    const Double_t kcylCthick      = 0.64*fgkmm; //
    const Double_t kcylFoamThick   = 5.0*fgkmm; //
    const Double_t kcylRholes      = 0.5*575.0*fgkmm;
    const Double_t kcylZM6         = 6.0*fgkmm; //
    const Double_t kcylRM6         = 0.5*6.0*fgkmm;
    const Double_t kcylPhi0M6      = 0.0*fgkDegree;
    const Int_t    kcylNM6         = 40;
    const Double_t kcylZPin        = 10.0*fgkmm;
    const Double_t kcylRPin        = 0.5*4.0*fgkmm;
    const Double_t kcylPhi0Pin     = (90.0+4.5)*fgkDegree;
    const Int_t    kcylNPin        = 2;
    //
    TGeoPcon *sCA,*sCB;
    TGeoTube *sCC,*sCD,*sCE;
    //
    //Begin_Html
    /*
      <img src="picts/ITS/file_name.gif">
      <P>
      <FONT FACE'"TIMES">
      ITS SSD centreal support and thermal sheal cylinder.
      </FONT>
      </P>
     */
    //End_Html
    //
    sCC = new TGeoTube("ITS SSD Thermal Centeral Rohacell CylinderCC",
                       kcylROuter-kcylCthick-kcylFoamThick,
                       kcylROuter-kcylCthick,0.5*kcylZFoamlength);
    sCA = new TGeoPcon("ITS SSD Thermal Centeral Carbon Fiber CylinderCA",
                       0.0,360.0,6);
    sCB = new TGeoPcon("ITS SSD Thermal Centeral Stesalite CylinderCB",
                       0.0,360.0,6);
    sCA->Z(0)    = -0.5*kcylZlength;
    sCA->Rmin(0) = kcylRInner;
    sCA->Rmax(0) = kcylROuter;
    sCA->Z(1)    = sCA->GetZ(0) + kcylZM6;
    sCA->Rmin(1) = sCA->GetRmin(0);
    sCA->Rmax(1) = sCA->GetRmax(0);
    sCA->Z(2)    = -0.5*kcylZFoamlength;
    sCA->Rmin(2) = kcylROuter - 2.0*kcylCthick-kcylFoamThick;
    sCA->Rmax(2) = sCA->GetRmax(0);
    sCA->Z(3)    = -sCA->GetZ(2);
    sCA->Rmin(3) = sCA->GetRmin(2);
    sCA->Rmax(3) = sCA->GetRmax(2);
    sCA->Z(4)    = -sCA->GetZ(1);
    sCA->Rmin(4) = sCA->GetRmin(1);
    sCA->Rmax(4) = sCA->GetRmax(1);
    sCA->Z(5)    = -sCA->GetZ(0);
    sCA->Rmin(5) = sCA->GetRmin(0);
    sCA->Rmax(5) = sCA->GetRmax(0);
    //
    sCB->Z(0)    = sCA->GetZ(0);
    sCB->Rmin(0) = sCA->GetRmin(0) + kcylCthick;
    sCB->Rmax(0) = sCA->GetRmax(0) - kcylCthick;
    sCB->Z(1)    = sCA->GetZ(1);
    sCB->Rmin(1) = sCA->GetRmin(1) + kcylCthick;
    sCB->Rmax(1) = sCA->GetRmax(1) - kcylCthick;
    sCB->Z(2)    = sCA->GetZ(2);
    sCB->Rmin(2) = sCA->GetRmin(2) + kcylCthick;
    sCB->Rmax(2) = sCA->GetRmax(2) - kcylCthick;
    sCB->Z(3)    = sCA->GetZ(3);
    sCB->Rmin(3) = sCA->GetRmin(3) + kcylCthick;
    sCB->Rmax(3) = sCA->GetRmax(3) - kcylCthick;
    sCB->Z(4)    = sCA->GetZ(4);
    sCB->Rmin(4) = sCA->GetRmin(4) + kcylCthick;
    sCB->Rmax(4) = sCA->GetRmax(4) - kcylCthick;
    sCB->Z(5)    = sCA->GetZ(5);
    sCB->Rmin(5) = sCA->GetRmin(5) + kcylCthick;
    sCB->Rmax(5) = sCA->GetRmax(5) - kcylCthick;
    //
    sCD = new TGeoTube("ITS SSD Thermal Centeral Cylinder M6 screwCD",
                      0.0,kcylRM6,0.5*kcylZM6);
    sCE = new TGeoTube("ITS SSD Thermal Centeral Cylinder PinCE",
                      0.0,kcylRPin,0.5*kcylZPin);
    //
    if(GetDebug()){
        sCA->InspectShape();
        sCB->InspectShape();
        sCC->InspectShape();
        sCD->InspectShape();
        sCE->InspectShape();
    } // end if GetDegut()
    TGeoVolume *vCA,*vCB,*vCC,*vCD,*vCE;
    vCA = new TGeoVolume("ITSssdCentCylCA",sCA,medSSDcf);
    vCA->SetVisibility(kTRUE);
    vCA->SetLineColor(4); // blue
    vCA->SetLineWidth(1);
    vCA->SetFillColor(vCA->GetLineColor());
    vCA->SetFillStyle(4000); // 0% transparent
    vCB = new TGeoVolume("ITSssdCentCylCB",sCB,medSSDfs);
    vCB->SetVisibility(kTRUE);
    vCB->SetLineColor(2); // red
    vCB->SetLineWidth(1);
    vCB->SetFillColor(vCB->GetLineColor());
    vCB->SetFillStyle(4050); // 50% transparent
    vCC = new TGeoVolume("ITSssdCentCylCC",sCC,medSSDfo);
    vCC->SetVisibility(kTRUE);
    vCC->SetLineColor(3); // green
    vCC->SetLineWidth(1);
    vCC->SetFillColor(vCC->GetLineColor());
    vCC->SetFillStyle(4050); // 50% transparent
    vCD = new TGeoVolume("ITSssdCentCylCD",sCD,medSSDss);
    vCD->SetVisibility(kTRUE);
    vCD->SetLineColor(1); // black
    vCD->SetLineWidth(1);
    vCD->SetFillColor(vCD->GetLineColor());
    vCD->SetFillStyle(4000); // 0% transparent
    vCE = new TGeoVolume("ITSssdCentCylCE",sCE,medSSDss);
    vCE->SetVisibility(kTRUE);
    vCE->SetLineColor(1); // black
    vCE->SetLineWidth(1);
    vCE->SetFillColor(vCE->GetLineColor());
    vCE->SetFillStyle(4000); // 0% transparent
    // Insert Bolt and Pins in both the Cone and Cylinder at the same time.
    vCB->AddNode(vCC,1,0);
    vCA->AddNode(vCB,1,0);
    moth->AddNode(vCA,1,0);
    if(GetDebug()){
        vCA->PrintNodes();
        vCB->PrintNodes();
        vCC->PrintNodes();
        vCD->PrintNodes();
        vCE->PrintNodes();
    } // end if
    //
    // SSD Cone
    // Data from Drawings ALR 0743/2E "Supporto Globale Settore SSD" and 
    // ALR 0743/2A "Supporto Generale Settore SSD".
    //
    const Double_t kconThick            = 13.0*fgkmm; // Thickness of Cone.
    const Double_t kconCthick           = 0.75*fgkmm; // Car. finber thickness
    const Double_t kconRCurv0           = 10.0*fgkmm; // Radius of curvature.
    const Double_t kconRCurv1           = 25.0*fgkmm; // Radius of curvature.
    const Double_t kconT                = 39.0*fgkDegree; // angle of SSD cone.
    const Double_t kconZOuterRing       = 47.0*fgkmm;
    const Double_t kconZOuterRingMill   = kconZOuterRing-5.0*fgkmm;
    const Double_t kconZToCylinder      = 170.0*fgkmm;
    const Double_t kconZLengthMill      = 171.5*fgkmm;
    const Double_t kconZLength          = 176.5*fgkmm-
                                          (kconZOuterRing-kconZOuterRingMill);
    //const Double_t kconZInnerRing       = 161.5*fgkmm-
    //                                     (kconZOuterRing-kconZOuterRingMill);
    const Double_t kconZOuterRingInside = 30.25*fgkmm-
                                          (kconZOuterRing-kconZOuterRingMill);
    const Double_t kconZDisplacement    = kconZToCylinder + 0.5*kcylZlength;
    const Double_t kconROuterMax        = 0.5*985.0*fgkmm;
    const Double_t kconROuterMin        = 0.5*945.0*fgkmm;
    const Double_t kconRCylOuterMill    = 0.5*597.0*fgkmm;
    const Double_t kconRInnerMin        = 0.5*562.0*fgkmm;
    //const Double_t kconRCentCurv0       = 0.5*927.0*fgkmm;
    const Double_t kconRCentCurv1       = 0.5*593.0*fgkmm;
    const Double_t kconRCentCurv2       = 0.5*578.0*fgkmm;
    // Foam core.
    const Double_t kconRohacellL0       = 112.3*fgkmm;
    const Double_t kconRohacellL1       = 58.4*fgkmm;
    // Screws and pins in outer SSD cone ring
    const Double_t kconROutHoles        = 0.5*965.0*fgkmm;
    const Double_t kconRScrewM5by12     = 0.5*5.0*fgkmm;
    const Double_t kconLScrewM5by12     = 0.5*12.0*fgkmm;
    const Int_t    kconNScrewM5by12     = 2;
    const Double_t kconRPinO6           = 0.5*6.0*fgkmm;
    const Double_t kconLPinO6           = 0.5*10.0*fgkmm;
    const Int_t    kconNPinO6           = 3;
    const Int_t    kconNRailScrews      = 4;
    const Int_t    kconNRailPins        = 2;
    const Int_t    kconNmounts          = 4;
    const Double_t kconMountPhi0        = 9.0*fgkDegree; // degrees
    //
    const Double_t kconCableHoleROut    = 0.5*920.0*fgkmm;
    const Double_t kconCableHoleRinner  = 0.5*800.0*fgkmm;
    const Double_t kconCableHoleWidth   = 200.0*fgkmm;
    const Double_t kconCableHoleAngle   = 42.0*fgkDegree;
    //const Double_t kconCableHolePhi0    = 90.0/4.0*fgkDegree;
    //const Int_t    kconNCableHoles      = 8;
    const Double_t kconCoolHoleWidth    = 40.0*fgkmm;
    const Double_t kconCoolHoleHight    = 30.0*fgkmm;
    const Double_t kconCoolHoleRmin     = 350.0*fgkmm;
    //const Double_t kconCoolHolephi0     = 90.0/4.0*fgkDegree;
    //const Int_t    kconNCoolHoles       = 8;
    const Double_t kconMountHoleWidth   = 20.0*fgkmm;
    const Double_t kconMountHoleHight   = 20.0*fgkmm;
    const Double_t kconMountHoleRmin    = 317.5*fgkmm;
    //const Double_t kconMountHolephi0    = 0.0*fgkDegree;
    //const Int_t    kconNMountHoles      = 6;
    // SSD cone Wings with holes.
    const Double_t kconWingRmax         = 527.5*fgkmm;
    const Double_t kconWingWidth        = 70.0*fgkmm;
    const Double_t kconWingThick        = 10.0*fgkmm;
    const Double_t kconWingPhi0         = 45.0*fgkDegree;
    //const Int_t    kconNWings           = 4;
    // SSD-SDD Thermal/Mechanical cylinder mounts
    const Double_t kconRM6Head          = 8.0*fgkmm;
    const Double_t kconZM6Head          = 8.5*fgkmm;
    //
    // SSD-SDD Mounting bracket
    const Double_t ksupPRmin            = 0.5*539.0*fgkmm;// see SDD RoutMin
    const Double_t ksupPRmax            = 0.5*585.0*fgkmm;
    const Double_t ksupPZ               = 3.5*fgkmm;
    const Double_t ksupPPhi1            = (-0.5*70.*fgkmm/ksupPRmax)*fgkRadian;
    const Double_t ksupPPhi2            = -ksupPPhi1;
    //
    const Double_t kSinkconTc           = SinD(kconT);
    const Double_t kCoskconTc           = CosD(kconT);
    //
    TGeoPcon *sA0,*sB0,*sC0,*sF0,*sQ;
    TGeoConeSeg *sAh1,*sBh1;
    TGeoArb8 *sAh2,*sBh2;
    TGeoBBox *sAh3,*sBh3,*sAh4,*sBh4;
    TGeoConeSeg *sG,*sH;
    TGeoTubeSeg *sT;
    TGeoTube *sD,*sE,*sR,*sS;
    TGeoCompositeShape *sA,*sB,*sC,*sF;
    //
    // Lets start with the upper left outer carbon fiber surface.
    // Between za[2],rmaxa[2] and za[4],rmaxa[4] there is a curved section
    // given by rmaxa = rmaxa[2]-r*Sind(t) for 0<=t<=kconT and 
    // za = za[2] + r*Cosd(t) for 0<=t<=kconT. Simularly between za[1],rmina[1
    // and za[3],rmina[3] there is a curve section given by
    // rmina = rmina[1]-r*Sind(t) for 0<=t<=kconT and za = za[1]+r&Sind(t)
    // for t<=0<=kconT. These curves have been replaced by straight lines
    // between the equivelent points for simplicity.
    // Poly-cone Volume sA0. Top part of SSD cone Carbon Fiber.
    sA0 = new TGeoPcon("ITSssdSuportConeCarbonFiberSurfaceA0",0.0,360.0,15);
    sA0->Z(0)    = 0.0;
    sA0->Rmin(0) = kconROuterMin;
    sA0->Rmax(0) = kconROuterMax;
    sA0->Z(1)    = kconZOuterRingInside-kconRCurv0;
    sA0->Rmin(1) = sA0->GetRmin(0);
    sA0->Rmax(1) = sA0->GetRmax(0);
    sA0->Z(2)    = kconZOuterRingInside;
    sA0->Rmin(2) = sA0->GetRmin(1)-kconRCurv0;
    sA0->Rmax(2) = sA0->GetRmax(0);
    sA0->Z(3)    = sA0->GetZ(2);
    sA0->Rmin(3) = -1000; // See Below
    sA0->Rmax(3) = sA0->GetRmax(0);
    sA0->Z(4)    = kconZOuterRingMill-kconRCurv0;
    sA0->Rmin(4) = -1000; // See Below
    sA0->Rmax(4) = sA0->GetRmax(0);
    sA0->Z(5)    = kconZOuterRingMill;
    sA0->Rmin(5) = -1000; // See Below
    sA0->Rmax(5) = sA0->GetRmax(4) - kconRCurv0;
    sA0->Z(6)    = sA0->GetZ(5);
    sA0->Rmin(6) = -1000; // See Below
    sA0->Rmax(6) = -1000; // See Below
    sA0->Z(7)    = sA0->GetZ(6)+kconRCurv0*(1.-kCoskconTc);
    sA0->Rmin(7) = -1000; // See Below
    sA0->Rmax(7) = -1000; // See Below
    sA0->Z(8)    = -1000; // See Below
    sA0->Rmin(8) = kconRCentCurv2+kconRCurv1*kSinkconTc; // See Below
    sA0->Rmax(8) = -1000; // See Below
    sA0->Z(9)    = -1000; // See Below
    sA0->Rmin(9) = kconRCentCurv2;
    sA0->Rmax(9) = -1000; // See Below
    sA0->Z(10)   = -1000; // See Below
    sA0->Rmin(10)= kconRInnerMin;
    sA0->Rmax(10)= -1000; // See Below
    sA0->Z(11)   = kconZLengthMill-kconRCurv0*(1.0-kCoskconTc);
    sA0->Rmin(11)= sA0->GetRmin(10);
    sA0->Rmax(11)= kconRCentCurv1+kconRCurv0*kSinkconTc;
    sA0->Z(12)   = kconZToCylinder;
    sA0->Rmin(12)= sA0->GetRmin(10);
    sA0->Rmax(12)= -1000; // See Below
    sA0->Z(13)   = sA0->GetZ(12);
    sA0->Rmin(13)= kconRCylOuterMill;
    sA0->Rmax(13)= -1000; // See Below
    z            = kconZLengthMill;
    rmin         = kconRCentCurv1;
    rmax         = rmin;
    sA0->Z(14)   = -1000; // See Below
    sA0->Rmin(14)= sA0->GetRmin(13);
    sA0->Rmax(14)= sA0->GetRmin(14);
    // Compute values undefined above
    sA0->Z(14)   = Xfrom2Points(sA0->GetZ(11),sA0->GetRmax(11),z,rmax,
                               sA0->GetRmax(14));
    sA0->Z(8)    = ZFromRmaxpCone(sA0,11,90.-kconT,sA0->GetRmin(8),-kconThick);
    sA0->Rmax(8) = RmaxFromZpCone(sA0,11,90.-kconT,sA0->GetZ(8),0.0);
    sA0->Z(9)    = sA0->GetZ(8)+kconRCurv1*(1.-kCoskconTc);
    sA0->Z(10)   = sA0->GetZ(9);
    sA0->Rmin(3) = RminFromZpCone(sA0,8,90.-kconT,sA0->GetZ(3),0.0);
    sA0->Rmin(4) = RminFromZpCone(sA0,3,90.-kconT,sA0->GetZ(4),0.0);
    sA0->Rmin(5) = RminFromZpCone(sA0,3,90.-kconT,sA0->GetZ(5),0.0);
    sA0->Rmin(7) = RminFromZpCone(sA0,3,90.-kconT,sA0->GetZ(7),0.0);
    sA0->Rmax(7) = RmaxFromZpCone(sA0,11,90.-kconT,sA0->GetZ(7),0.0);
    sA0->Rmin(6) = sA0->GetRmin(5);
    sA0->Rmax(6) = RmaxFromZpCone(sA0,11,90.-kconT,sA0->GetZ(7),0.0);
    sA0->Rmax(9) = RmaxFromZpCone(sA0,11,90.-kconT,sA0->GetZ(9),0.0);
    sA0->Rmax(10)= sA0->GetRmax(9);
    t = TanD(270.+kconT);
    sA0->Rmax(12)= RmaxFrom2Points(sA0,11,14,sA0->GetZ(12));
    sA0->Rmax(13)= sA0->GetRmax(12);
    //
    // Poly-cone Volume B. Stesalite inside volume sA0.
    // Now lets define the Inserto Stesalite 4411w material volume.
    // Poly-cone Volume sA0. Top part of SSD cone Carbon Fiber.
    sB0 = new TGeoPcon("ITSssdSuportConeStaseliteB0",0.0,360.0,15);
    //
    sB0->Z(0)    = sA0->GetZ(0);
    sB0->Rmin(0) = sA0->GetRmin(0) + kconCthick;
    sB0->Rmax(0) = sA0->GetRmax(0) - kconCthick;
    InsidePoint(sA0,0,1,2,kconCthick,sB0,1,kFALSE); // Rmin
    sB0->Rmax(1) = sB0->Rmax(0);
    InsidePoint(sA0,1,2,3,kconCthick,sB0,2,kFALSE); // Rmin
    sB0->Rmax(2) = sB0->Rmax(0);
    InsidePoint(sA0,2,3,9,kconCthick,sB0,3,kFALSE);
    sB0->Rmax(3) = sB0->Rmax(0);
    InsidePoint(sA0,0,4,5,kconCthick,sB0,4,kTRUE); // Rmax
    sB0->Rmin(4) = -1000.; // see Bellow
    InsidePoint(sA0,4,5,6,kconCthick,sB0,5,kTRUE); // Rmax
    sB0->Rmin(5) = -1000.; // see Bellow
    InsidePoint(sA0,5,6,7,kconCthick,sB0,6,kTRUE); // Rmax
    sB0->Rmin(6) = -1000.; // see Bellow
    InsidePoint(sA0,6,7,11,kconCthick,sB0,7,kTRUE); // Rmax
    sB0->Rmin(7) = -1000.; // see Bellow
    InsidePoint(sA0,3,8,9,kconCthick,sB0,8,kFALSE); // Rmin
    sB0->Rmax(8) = -1000.; // see Bellow
    InsidePoint(sA0,8,9,10,kconCthick,sB0,9,kFALSE); // Rmin
    sB0->Rmax(9) = -1000.; // see Bellow
    sB0->Z(10)   = sA0->GetZ(10) + kconCthick;
    sB0->Rmin(10)= sA0->GetRmin(10);
    sB0->Rmax(10)= -1000.; // see Bellow
    InsidePoint(sA0,7,11,14,kconCthick,sB0,11,kTRUE); // Rmax
    sB0->Rmin(11)= sA0->GetRmin(10);
    sB0->Z(12)    = sA0->GetZ(12);
    sB0->Rmin(12)= sA0->GetRmin(12);
    sB0->Rmax(12)= -1000.; // see Bellow
    sB0->Z(13)   = sA0->GetZ(13);
    sB0->Rmin(13)= sA0->GetRmin(13);
    sB0->Rmax(13)= -1000.; // see Bellow
    sB0->Z(14)   = sA0->GetZ(14) - kconCthick;
    sB0->Rmin(14)= sA0->GetRmin(14);
    sB0->Rmax(14)= sB0->Rmin(14); // Close?
    sB0->Rmin(4) = RminFrom2Points(sB0,3,8,sB0->GetZ(4));
    sB0->Rmin(5) = RminFrom2Points(sB0,3,8,sB0->GetZ(5));
    sB0->Rmin(6) = sB0->GetRmin(5);
    sB0->Rmin(7) = RminFrom2Points(sB0,3,8,sB0->GetZ(7));
    sB0->Rmax(8) = RmaxFrom2Points(sB0,7,11,sB0->GetZ(8));
    sB0->Rmax(9) = RmaxFrom2Points(sB0,7,11,sB0->GetZ(9));
    sB0->Rmax(10)= sB0->GetRmax(9);
    sB0->Rmax(12)= RmaxFrom2Points(sB0,11,14,sB0->GetZ(12));
    sB0->Rmax(13)= RmaxFrom2Points(sB0,11,14,sB0->GetZ(13));
    //
    // Poly-cone Volume sC0. Foam inside volume sA0.
    // Now lets define the Rohacell foam material volume.
    sC0 = new TGeoPcon("ITSssdSuportConeRohacellC0",0.0,360.0,4);
    sC0->Z(1)    = sB0->GetZ(7);
    sC0->Rmax(1) = sB0->GetRmax(7);
    sC0->Rmin(1) = RminFrom2Points(sB0,3,8,sC0->GetZ(1));
    sC0->Rmin(0) = sC0->GetRmax(1);
    sC0->Rmax(0) = sC0->GetRmin(0);
    sC0->Z(0)    = Zfrom2MinPoints(sB0,3,8,sC0->Rmin(0));
    t = kconThick-2.0*kconCthick;
    sC0->Rmax(3) = sC0->GetRmax(0)-kCoskconTc*TMath::Sqrt(
                             kconRohacellL0*kconRohacellL0-t*t)+t*kSinkconTc;
    sC0->Rmin(3) = sC0->GetRmax(3);
    sC0->Z(3)    = ZFromRmaxpCone(sB0,11,90.-kconT,sC0->GetRmax(3),0.0);;
    sC0->Rmin(2) = sC0->GetRmin(3);
    sC0->Z(2)    = ZFromRminpCone(sB0,3,90.-kconT,sC0->GetRmin(2),0.0);
    sC0->Rmax(2) = RmaxFromZpCone(sB0,11,90.0-kconT,sC0->GetZ(2),0.0);
    //
    // Poly-cone Volume sF0.  Second Foam inside volume sA0.
    // Now lets define the Rohacell foam material volume.
    sF0 = new TGeoPcon("ITSssdSuportConeRohacellCF0",0.0,360.0,4);
    sF0->Z(2)    = sB0->GetZ(8);
    sF0->Rmin(2) = sB0->GetRmin(8);
    sF0->Rmax(2) = sB0->GetRmax(8);
    sF0->Z(0)    = sF0->GetZ(2)-kconRohacellL1*kSinkconTc;
    sF0->Rmin(0) = sF0->GetRmin(2)+kconRohacellL1*kCoskconTc;
    sF0->Rmax(0) = sF0->GetRmin(0);
    sF0->Z(1)    = ZFromRmaxpCone(sB0,11,90.-kconT,sF0->GetRmax(0),0.0);;
    sF0->Rmax(1) = sF0->GetRmax(0);
    sF0->Rmin(1) = RminFrom2Points(sB0,3,8,sF0->GetZ(1));
    sF0->Rmax(3) = sF0->GetRmin(2)+(kconThick-2.0*kconCthick)*kCoskconTc;
    sF0->Rmin(3) = sF0->GetRmax(3);
    sF0->Z(3)    = ZFromRmaxpCone(sB0,11,90.-kconT,sF0->GetRmax(3),0.0);
    // Holes for Cables to pass Through is created by the intersection
    // between a cone segment and an Arb8, One for the volume sA0 and a
    // larger one for the volumes sB0 and sC0, so that the surface is covered
    // in carbon figer (volume sA0).
    sAh1 = new TGeoConeSeg("ITSssdCableHoleAh1",
                           0.5*kconZLength,kconCableHoleRinner,
                           kconCableHoleROut,kconCableHoleRinner,
                           kconCableHoleROut,
                           90.-(0.5*kconCableHoleWidth/
                                kconCableHoleROut)*fgkRadian,
                           90.+(0.5*kconCableHoleWidth/
                                kconCableHoleROut)*fgkRadian);
    sBh1 = new TGeoConeSeg("ITSssdCableHoleBh1",0.5*kconZLength,
                           kconCableHoleRinner-kconCthick,
                           kconCableHoleROut+kconCthick,
                           kconCableHoleRinner-kconCthick,
                           kconCableHoleROut+kconCthick,
                           90.-(((0.5*kconCableHoleWidth+kconCthick)/
                                 (kconCableHoleROut+kconCthick)))*fgkRadian,
                           90.+(((0.5*kconCableHoleWidth+kconCthick)/
                                 (kconCableHoleROut+kconCthick)))*fgkRadian);
    x0 = sAh1->GetRmax1()*CosD(sAh1->GetPhi2());
    y0 = sAh1->GetRmax1()*SinD(sAh1->GetPhi2());
    sAh2 = new TGeoArb8("ITSssdCableHoleAh2",0.5*kconZLength);
    y  = sAh1->GetRmax1();
    x  = x0+(y-y0)/TanD(90.0+kconCableHoleAngle);
    sAh2->SetVertex(0,x,y);
    y  = sAh1->GetRmin1()*SinD(sAh1->GetPhi2());
    x  = x0+(y-y0)/TanD(90.0+kconCableHoleAngle);
    sAh2->SetVertex(3,x,y);
    x0 = sAh1->GetRmax1()*CosD(sAh1->GetPhi1());
    y0 = sAh1->GetRmax1()*SinD(sAh1->GetPhi1());
    y  = sAh1->GetRmax1();
    x  = x0+(y-y0)/TanD(90.0-kconCableHoleAngle);
    sAh2->SetVertex(1,x,y);
    y  = sAh1->GetRmin1()*SinD(sAh1->GetPhi1());
    x  = x0+(y-y0)/TanD(90.0-kconCableHoleAngle);
    sAh2->SetVertex(2,x,y);
    //
    x0 = sBh1->GetRmax1()*CosD(sBh1->GetPhi2());
    y0 = sBh1->GetRmax1()*SinD(sBh1->GetPhi2());
    sBh2 = new TGeoArb8("ITSssdCableHoleBh2",0.5*kconZLength);
    y  = sBh1->GetRmax1();
    x  = x0+(y-y0)/TanD(90.0+kconCableHoleAngle);
    sBh2->SetVertex(0,x,y);
    y  = sBh1->GetRmin1()*SinD(sBh1->GetPhi2());
    x  = x0+(y-y0)/TanD(90.0+kconCableHoleAngle);
    sBh2->SetVertex(3,x,y);
    x0 = sBh1->GetRmax1()*CosD(sBh1->GetPhi1());
    y0 = sBh1->GetRmax1()*SinD(sBh1->GetPhi1());
    y  = sBh1->GetRmax1();
    x  = x0+(y-y0)/TanD(90.0-kconCableHoleAngle);
    sBh2->SetVertex(1,x,y);
    y  = sBh1->GetRmin1()*SinD(sBh1->GetPhi1());
    x  = x0+(y-y0)/TanD(90.0-kconCableHoleAngle);
    sBh2->SetVertex(2,x,y);
    for(i=0;i<4;i++){ // define points at +dz
        sAh2->SetVertex(i+4,(sAh2->GetVertices())[2*i],
                           (sAh2->GetVertices())[1+2*i]);
        sBh2->SetVertex(i+4,(sBh2->GetVertices())[2*i],
                           (sBh2->GetVertices())[1+2*i]);
    } // end for i
    sAh3 = new TGeoBBox("ITSssdCoolingHoleAh3",0.5*kconCoolHoleWidth,
                        0.5*kconCoolHoleHight,kconZLength);
    sBh3 = new TGeoBBox("ITSssdCoolingHoleBh3",
                        0.5*kconCoolHoleWidth+kconCthick,
                        0.5*kconCoolHoleHight+kconCthick,kconZLength);
    sAh4 = new TGeoBBox("ITSssdMountingPostHoleAh4",0.5*kconMountHoleWidth,
                        0.5*kconMountHoleHight,0.5*kconZLength);
    z = sF0->GetZ(0)-sF0->GetZ(sF0->GetNz()-1);
    if(z<0.0) z = -z;
    sBh4 = new TGeoBBox("ITSssdMountingPostHoleBh4",
                        0.5*kconMountHoleWidth+kconCthick,
                        0.5*kconMountHoleHight+kconCthick,0.5*z);
    // SSD Cone Wings
    sG = new TGeoConeSeg("ITSssdWingCarbonFiberSurfaceG",
                         0.5*kconWingThick,kconROuterMax-kconCthick,
                         kconWingRmax,kconROuterMax-kconCthick,kconWingRmax,
                      kconWingPhi0-(0.5*kconWingWidth/kconWingRmax)*fgkRadian,
                      kconWingPhi0+(0.5*kconWingWidth/kconWingRmax)*fgkRadian);
    sH = new TGeoConeSeg("ITSssdWingStaseliteH",
                         0.5*kconWingThick-kconCthick,kconROuterMax-kconCthick,
                         kconWingRmax-kconCthick,
                         kconROuterMax-kconCthick,
                         kconWingRmax-kconCthick,
                         kconWingPhi0-((0.5*kconWingWidth-kconCthick)/
                                       (kconWingRmax-kconCthick))*fgkRadian,
                         kconWingPhi0+((0.5*kconWingWidth-kconCthick)/
                                       (kconWingRmax-kconCthick))*fgkRadian);
    // SDD support plate, SSD side.
    //Poly-cone Volume sT.
    sT = new TGeoTubeSeg("ITSssdsddMountingBracketT",ksupPRmin,ksupPRmax,
                         ksupPZ,ksupPPhi1,ksupPPhi2);
    //
    TGeoRotation *rotZ225 =new TGeoRotation("ITSssdConeZ225", 0.0,0.0, 22.5);
    rotZ225->RegisterYourself();
    TGeoRotation *rotZ675 =new TGeoRotation("ITSssdConeZ675", 0.0,0.0, 67.5);
    rotZ675->RegisterYourself();
    TGeoRotation *rotZ90  =new TGeoRotation("ITSssdConeZ90",  0.0,0.0, 90.0);
    rotZ90->RegisterYourself();
    TGeoRotation *rotZ1125=new TGeoRotation("ITSssdConeZ1125",0.0,0.0,112.5);
    rotZ1125->RegisterYourself();
    TGeoRotation *rotZ1575=new TGeoRotation("ITSssdConeZ1575",0.0,0.0,157.5);
    rotZ1575->RegisterYourself();
    TGeoRotation *rotZ180 =new TGeoRotation("ITSssdConeZ180", 0.0,0.0,180.0);
    rotZ180->RegisterYourself();
    TGeoRotation *rotZ2025=new TGeoRotation("ITSssdConeZ2025",0.0,0.0,202.5);
    rotZ2025->RegisterYourself();
    TGeoRotation *rotZ2475=new TGeoRotation("ITSssdConeZ2475",0.0,0.0,247.5);
    rotZ2475->RegisterYourself();
    TGeoRotation *rotZ270 =new TGeoRotation("ITSssdConeZ270", 0.0,0.0,270.0);
    rotZ270->RegisterYourself();
    TGeoRotation *rotZ2925=new TGeoRotation("ITSssdConeZ2925",0.0,0.0,292.5);
    rotZ2925->RegisterYourself();
    TGeoRotation *rotZ3375=new TGeoRotation("ITSssdConeZ3375",0.0,0.0,337.5);
    rotZ3375->RegisterYourself();
    //
    vl[0] = 0.0;vl[1] = kconCoolHoleRmin+0.5*kconCoolHoleHight;vl[2] = 0.0;
    rotZ225->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA225  = new TGeoCombiTrans("ITSssdConeTZ225",vg[0],
                                                     vg[1],vg[2],rotZ225);
    rotranA225->RegisterYourself();
    rotZ675->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA675  = new TGeoCombiTrans("ITSssdConeTZ675", vg[0],
                                                     vg[1],vg[2],rotZ675);
    rotranA675->RegisterYourself();
    rotZ1125->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA1125 = new TGeoCombiTrans("ITSssdConeTZ1125",vg[0],
                                                     vg[1],vg[2],rotZ1125);
    rotranA1125->RegisterYourself();
    rotZ1575->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA1575 = new TGeoCombiTrans("ITSssdConeTZ1575",vg[0],
                                                     vg[1],vg[2],rotZ1575);
    rotranA1575->RegisterYourself();
    rotZ2025->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA2025 = new TGeoCombiTrans("ITSssdConeTZ2025",vg[0],
                                                     vg[1],vg[2],rotZ2025);
    rotranA2025->RegisterYourself();
    rotZ2475->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA2475 = new TGeoCombiTrans("ITSssdConeTZ2475",vg[0],
                                                     vg[1],vg[2],rotZ2475);
    rotranA2475->RegisterYourself();
    rotZ2925->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA2925 = new TGeoCombiTrans("ITSssdConeTZ2925",vg[0],
                                                     vg[1],vg[2],rotZ2925);
    rotranA2925->RegisterYourself();
    rotZ3375->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA3375 = new TGeoCombiTrans("ITSssdConeTZ3375",vg[0],
                                                     vg[1],vg[2],rotZ3375);
    rotranA3375->RegisterYourself();
    TGeoRotation *rotZ30  = new TGeoRotation("ITSssdConeZ30", 0.0,0.0, 30.0);
    TGeoRotation *rotZ60  = new TGeoRotation("ITSssdConeZ60", 0.0,0.0, 60.0);
    //TGeoRotation *rotZ120 = new TGeoRotation("ITSssdConeZ120",0.0,0.0,120.0);
    TGeoRotation *rotZ150 = new TGeoRotation("ITSssdConeZ150",0.0,0.0,150.0);
    TGeoRotation *rotZ210 = new TGeoRotation("ITSssdConeZ210",0.0,0.0,210.0);
    //TGeoRotation *rotZ240 = new TGeoRotation("ITSssdConeZ240",0.0,0.0,240.0);
    TGeoRotation *rotZ300 = new TGeoRotation("ITSssdConeZ300",0.0,0.0,300.0);
    TGeoRotation *rotZ330 = new TGeoRotation("ITSssdConeZ330",0.0,0.0,330.0);
    vl[0] = kconMountHoleRmin+0.5*kconMountHoleHight; vl[1] = 0.0; vl[2] = 0.0;
    for(i=0;i<sF0->GetNz();i++) vl[2] += sF0->GetZ(i);
    vl[2] /= (Double_t)(sF0->GetNz());
    rotZ30->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA30 = new TGeoCombiTrans("ITSssdConeTZ30",vg[0],
                                                      vg[1],vg[2],rotZ30);
    rotranA30->RegisterYourself();
    rotZ90->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA90  = new TGeoCombiTrans("ITSssdConeTZ90", vg[0],
                                                     vg[1],vg[2],rotZ90);
    rotranA90->RegisterYourself();
    rotZ150->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA150 = new TGeoCombiTrans("ITSssdConeTZ150",vg[0],
                                                     vg[1],vg[2],rotZ150);
    rotranA150->RegisterYourself();
    rotZ210->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA210 = new TGeoCombiTrans("ITSssdConeTZ210",vg[0],
                                                     vg[1],vg[2],rotZ210);
    rotranA210->RegisterYourself();
    rotZ270->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA270 = new TGeoCombiTrans("ITSssdConeTZ270",vg[0],
                                                     vg[1],vg[2],rotZ270);
    rotranA270->RegisterYourself();
    rotZ330->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA330 = new TGeoCombiTrans("ITSssdConeTZ330",vg[0],
                                                     vg[1],vg[2],rotZ330);
    rotranA330->RegisterYourself();
    vl[0] = 0.0; vl[1] = 0.0; vl[2] = sA0->GetZ(10)+sT->GetDz();
    rotZ60->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranBrTZ60  = new TGeoCombiTrans("ITSssdConeBrTZ60",
                                                  vg[0],vg[1],vg[2],rotZ60);
    rotranBrTZ60->RegisterYourself();
    TGeoCombiTrans *rotranBrTZ180 = new TGeoCombiTrans("ITSssdConeBrTZ180",
                                                  vg[0],vg[1],vg[2],rotZ180);
    rotranBrTZ180->RegisterYourself();
    TGeoCombiTrans *rotranBrTZ300 = new TGeoCombiTrans("ITSssdConeBrTZ300",
                                                  vg[0],vg[1],vg[2],rotZ300);
    rotranBrTZ300->RegisterYourself();
    if(GetDebug()){
        rotZ225->Print();
        rotZ675->Print();
        rotZ90->Print();
        rotZ1125->Print();
        rotZ1575->Print();
        rotZ180->Print();
        rotZ2025->Print();
        rotZ2475->Print();
        rotZ270->Print();
        rotZ2925->Print();
        rotZ3375->Print();
        rotranA225->Print();
        rotranA675->Print();
        rotranA1125->Print();
        rotranA1575->Print();
        rotranA2025->Print();
        rotranA2475->Print();
        rotranA2925->Print();
        rotranA3375->Print();
        rotZ60->Print();
        rotZ300->Print();
        rotranA30->Print();
        rotranA90->Print();
        rotranA150->Print();
        rotranA210->Print();
        rotranA270->Print();
        rotranA330->Print();
        rotranBrTZ60->Print();
        rotranBrTZ180->Print();
        rotranBrTZ300->Print();
    } // end if GetDebug()
    sA = new TGeoCompositeShape("ITSssdSuportConeCarbonFiberSurfaceA",
        "(((((((((((((((((((((((((((("
        "ITSssdSuportConeCarbonFiberSurfaceA0 +"
        "ITSssdWingCarbonFiberSurfaceG) +"
        "ITSssdWingCarbonFiberSurfaceG:ITSssdConeZ90) +"
        "ITSssdWingCarbonFiberSurfaceG:ITSssdConeZ180) +"
        "ITSssdWingCarbonFiberSurfaceG:ITSssdConeZ270) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ225*ITSssdCableHoleAh2:ITSssdConeZ225)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ675*ITSssdCableHoleAh2:ITSssdConeZ675)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ1125*ITSssdCableHoleAh2:ITSssdConeZ1125)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ1575*ITSssdCableHoleAh2:ITSssdConeZ1575)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ2025*ITSssdCableHoleAh2:ITSssdConeZ2025)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ2475*ITSssdCableHoleAh2:ITSssdConeZ2475)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ2925*ITSssdCableHoleAh2:ITSssdConeZ2925)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ3375*ITSssdCableHoleAh2:ITSssdConeZ3375)) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ225) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ675) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ1125) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ1575) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ2025) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ2475) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ2925) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ3375) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ30) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ90) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ150) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ210) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ270) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ330) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ60) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ180) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ300"
        );
    sB = new TGeoCompositeShape("ITSssdSuportConeStaseliteB",
        "(((((((((((((((((((((((((((("
        "ITSssdSuportConeStaseliteB0 +"
        "ITSssdWingStaseliteH) +"
        "ITSssdWingStaseliteH:ITSssdConeZ90) +"
        "ITSssdWingStaseliteH:ITSssdConeZ180) +"
        "ITSssdWingStaseliteH:ITSssdConeZ270) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ225*ITSssdCableHoleBh2:ITSssdConeZ225)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ675*ITSssdCableHoleBh2:ITSssdConeZ675)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ1125*ITSssdCableHoleBh2:ITSssdConeZ1125)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ1575*ITSssdCableHoleBh2:ITSssdConeZ1575)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ2025*ITSssdCableHoleBh2:ITSssdConeZ2025)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ2475*ITSssdCableHoleBh2:ITSssdConeZ2475)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ2925*ITSssdCableHoleBh2:ITSssdConeZ2925)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ3375*ITSssdCableHoleBh2:ITSssdConeZ3375)) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ225) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ675) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ1125) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ1575) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ2025) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ2475) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ2925) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ3375) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ30) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ90) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ150) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ210) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ270) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ330) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ60) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ180) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ300"
        );
    sC = new TGeoCompositeShape("ITSssdSuportConeRohacellC",
      "((((((("
      "ITSssdSuportConeRohacellC0 -"
      "ITSssdCableHoleBh1:ITSssdConeZ225*ITSssdCableHoleBh2:ITSssdConeZ225) -"
      "ITSssdCableHoleBh1:ITSssdConeZ675*ITSssdCableHoleBh2:ITSssdConeZ675) -"
      "ITSssdCableHoleBh1:ITSssdConeZ1125*ITSssdCableHoleBh2:ITSssdConeZ1125) -"
      "ITSssdCableHoleBh1:ITSssdConeZ1575*ITSssdCableHoleBh2:ITSssdConeZ1575) -"
      "ITSssdCableHoleBh1:ITSssdConeZ2025*ITSssdCableHoleBh2:ITSssdConeZ2025) -"
      "ITSssdCableHoleBh1:ITSssdConeZ2475*ITSssdCableHoleBh2:ITSssdConeZ2475) -"
      "ITSssdCableHoleBh1:ITSssdConeZ2925*ITSssdCableHoleBh2:ITSssdConeZ2925) -"
      "ITSssdCableHoleBh1:ITSssdConeZ3375*ITSssdCableHoleBh2:ITSssdConeZ3375 "
        );
    sF = new TGeoCompositeShape("ITSssdSuportConeRohacellCF",
        "((((("
        "ITSssdSuportConeRohacellCF0 -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ30) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ90) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ150) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ210) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ270) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ330"
        );
    //
    // In volume SCB, th Inserto Stesalite 4411w material volume, there
    // are a number of Stainless steel screw and pin studs which will be
    // filled with screws/studs.
    sD = new TGeoTube("ITS Screw+stud used to mount things to the SSD "
                      "support cone",
                      0.0,kconRScrewM5by12,kconLScrewM5by12);
    sE = new TGeoTube("ITS pin used to mount things to the "
                      "SSD support cone",0.0,kconRPinO6,kconLPinO6);
    // Bolt heads holding the SSD-SDD tube to the SSD cone.
    // Bolt -- PolyCone
    //Poly-cone Volume sQ.
    sQ = new TGeoPcon("ITS SSD Thermal sheal M6 screw headQ",0.0,360.0,4);
    sQ->Z(0)    = sA0->GetZ(12);
    sQ->Rmin(0) = 0.0;
    sQ->Rmax(0) = kcylRM6;
    sQ->Z(1)    = sQ->GetZ(0) - kconZM6Head;
    sQ->Rmin(1) = 0.0;
    sQ->Rmax(1) = kcylRM6;
    sQ->Z(2)    = sQ->GetZ(1);
    sQ->Rmin(2) = 0.0;
    sQ->Rmax(2) = kconRM6Head;
    sQ->Z(3)    = sQ->GetZ(0)-ksupPZ;
    sQ->Rmin(3) = 0.0;
    sQ->Rmax(3) = 0.5*kconRM6Head;
    // air infront of bolt (stasolit Volume K) -- Tube
    sR = new TGeoTube("ITS Air in front of bolt (in stasolit)R",
                      sQ->GetRmin(3),sQ->GetRmax(3),0.5*(ksupPZ-kconCthick));
    // air infront of bolt (carbon fiber volume I) -- Tube
    sS = new TGeoTube("ITS Air in front of Stainless Steal Screw end, M6S",
                      sQ->GetRmin(3),sQ->GetRmax(3),0.5*kconCthick);
    //
    if(GetDebug()){
        sA0->InspectShape();
        sB0->InspectShape();
        sC0->InspectShape();
        sF0->InspectShape();
        sQ->InspectShape();
        sAh1->InspectShape();
        sBh1->InspectShape();
        sAh2->InspectShape();
        sBh2->InspectShape();
        sAh3->InspectShape();
        sBh3->InspectShape();
        sAh4->InspectShape();
        sBh4->InspectShape();
        sG->InspectShape();
        sH->InspectShape();
        sT->InspectShape();
        sD->InspectShape();
        sE->InspectShape();
        sR->InspectShape();
        sS->InspectShape();
        sA->InspectShape();
        sB->InspectShape();
        sC->InspectShape();
        sF->InspectShape();
    } // end if GetDebug()
    TGeoVolume *vA,*vB,*vC,*vD,*vE,*vF,*vQ,*vR,*vS,*vT;
    //
    vA = new TGeoVolume("ITSssdConeA",sA,medSSDcf); // Carbon Fiber
    vA->SetVisibility(kTRUE);
    vA->SetLineColor(4); // blue
    vA->SetLineWidth(1);
    vA->SetFillColor(vA->GetLineColor());
    vA->SetFillStyle(4050); // 50% transparent
    vB = new TGeoVolume("ITSssdConeB",sB,medSSDfs); // Staselite
    vB->SetVisibility(kTRUE);
    vB->SetLineColor(2); // red
    vB->SetLineWidth(1);
    vB->SetFillColor(vB->GetLineColor());
    vB->SetFillStyle(4050); // 50% transparent
    vC = new TGeoVolume("ITSssdConeC",sC,medSSDfo); // Rohacell
    vC->SetVisibility(kTRUE);
    vC->SetLineColor(3); // green
    vC->SetLineWidth(1);
    vC->SetFillColor(vC->GetLineColor());
    vC->SetFillStyle(4050); // 50% transparent
    vF = new TGeoVolume("ITSssdConeF",sF,medSSDfo); // Rohacell;
    vF->SetVisibility(kTRUE);
    vF->SetLineColor(3); // green
    vF->SetLineWidth(1);
    vF->SetFillColor(vF->GetLineColor());
    vF->SetFillStyle(4050); // 50% transparent
    vD = new TGeoVolume("ITSssdConeD",sD,medSSDss);
    vD->SetVisibility(kTRUE);
    vD->SetLineColor(1); // black
    vD->SetLineWidth(1);
    vD->SetFillColor(vD->GetLineColor());
    vD->SetFillStyle(4000); // 0% transparent
    vE = new TGeoVolume("ITSssdConeE",sE,medSSDss);
    vE->SetVisibility(kTRUE);
    vE->SetLineColor(1); // black
    vE->SetLineWidth(1);
    vE->SetFillColor(vE->GetLineColor());
    vE->SetFillStyle(4000); // 0% transparent
    vQ = new TGeoVolume("ITSssdConeQ",sQ,medSSDss);
    vQ->SetVisibility(kTRUE);
    vQ->SetLineColor(1); // black
    vQ->SetLineWidth(1);
    vQ->SetFillColor(vQ->GetLineColor());
    vQ->SetFillStyle(4000); // 0% transparent
    vR = new TGeoVolume("ITSssdConeR",sR,medSSDair);
    vR->SetVisibility(kTRUE);
    vR->SetLineColor(5); // yellow
    vR->SetLineWidth(1);
    vR->SetFillColor(vR->GetLineColor());
    vR->SetFillStyle(4090); // 90% transparent
    vS = new TGeoVolume("ITSssdConeS",sS,medSSDair);
    vS->SetVisibility(kTRUE);
    vS->SetLineColor(5); // yellow
    vS->SetLineWidth(1);
    vS->SetFillColor(vS->GetLineColor());
    vS->SetFillStyle(4090); // 90% transparent
    vT = new TGeoVolume("ITSssdsddMountingBracket",sT,medSSDal);
    vT->SetVisibility(kTRUE);
    vT->SetLineColor(5); // yellow
    vT->SetLineWidth(1);
    vT->SetFillColor(vT->GetLineColor());
    vT->SetFillStyle(4000); // 0% transparent
    //
    TGeoCombiTrans *rotran;
    TGeoTranslation *tran;
    tran = new TGeoTranslation("ITSssdConeTrans",0.0,0.0,-kconZDisplacement);
    TGeoRotation *rotY180 = new TGeoRotation("",0.0,180.0,0.0);
    TGeoCombiTrans *flip  = new TGeoCombiTrans("ITSssdConeFlip",
                                           0.0,0.0,kconZDisplacement,rotY180);
    delete rotY180;// rot not explicity used in AddNode functions.
    //
    //
    //
    //
    vA->AddNode(vB,1,0);
    vB->AddNode(vC,1,0);
    vB->AddNode(vF,1,0);
    moth->AddNode(vA,1,tran); // RB24 side
    moth->AddNode(vA,2,flip); // RB26 side (Absorber)
    //
    //
    //
    // Insert Bolt and Pins in both the Cone and Cylinder at the same time.
    Int_t nCopyCDv=0,nCopyCEv=0,nCopyQv=0,nCopyvR=0,nCopySv=0,nCopyTv=0;
    Int_t nCopyvD=0,nCopyvE=0;
    z = sCB->GetZ(0)-0.5*kcylZPin;
    dt = (360.0/((Double_t)kcylNPin));
    for(i=0;i<kcylNPin;i++){
        t = ((Double_t)i)*dt;
        x = kcylRholes*CosD(t+kcylPhi0Pin);
        y = kcylRholes*SinD(t+kcylPhi0Pin);
        tran = new TGeoTranslation("",x,y,z);
        vCB->AddNode(vCD,++nCopyCDv,tran);
        tran = new TGeoTranslation("",x,y,-z);
        vCB->AddNode(vCD,++nCopyCDv,tran);
    } // end for i
    dt = (360.0/((Double_t)kcylNM6));
    for(i=0;i<kcylNM6;i++){
        t = ((Double_t)i)*dt;
        x = kcylRholes*CosD(t+kcylPhi0M6);
        y = kcylRholes*SinD(t+kcylPhi0M6);
        z = sCB->GetZ(0)-0.5*kcylZM6;
        tran = new TGeoTranslation("",x,y,z);
        vCB->AddNode(vCE,++nCopyCEv,tran);
        tran = new TGeoTranslation("",x,y,-z);
        vCB->AddNode(vCE,++nCopyCEv,tran);
        tran = new TGeoTranslation("",x,y,0.0);
        vB->AddNode(vQ,++nCopyQv,tran);
        if(!((t<rotranBrTZ60->GetRotation()->GetPhiRotation()+sT->GetPhi2()&&
             t>rotranBrTZ60->GetRotation()->GetPhiRotation()-sT->GetPhi1())||
            (t<rotranBrTZ180->GetRotation()->GetPhiRotation()+sT->GetPhi2()&&
             t>rotranBrTZ180->GetRotation()->GetPhiRotation()-sT->GetPhi1())||
            (t<rotranBrTZ300->GetRotation()->GetPhiRotation()+sT->GetPhi2()&&
             t>rotranBrTZ300->GetRotation()->GetPhiRotation()-sT->GetPhi1()))){
            // If not at an angle where the bracket sT is located.
            tran = new TGeoTranslation("",x,y,sB0->GetZ(10)-sR->GetDz());
            vB->AddNode(vR,++nCopyvR,tran);
            tran = new TGeoTranslation("",x,y,sA0->GetZ(10)-sS->GetDz());
            vA->AddNode(vS,++nCopySv,tran);
        } // end if
    } // end for i
    // Add the mounting brackets to the RB24 side only.
    vl[0] = 0.0;
    vl[1] = 0.0;
    vl[2] = sA0->GetZ(10)+kconZDisplacement-sT->GetDz();
    rotZ60->LocalToMaster(vl,vg);
    rotran = new TGeoCombiTrans("",vg[0],vg[1],vg[2],rotZ60);
    moth->AddNode(vT,++nCopyTv,rotran);
    rotZ180->LocalToMaster(vl,vg);
    rotran = new TGeoCombiTrans("",vg[0],vg[1],vg[2],rotZ180);
    moth->AddNode(vT,++nCopyTv,rotran);
    rotZ300->LocalToMaster(vl,vg);
    rotran = new TGeoCombiTrans("",vg[0],vg[1],vg[2],rotZ300);
    moth->AddNode(vT,++nCopyTv,rotran);
    //
    Double_t da[] = {-3.5,-1.5,1.5,3.5};
    for(i=0;i<2;i++){ // Mounting for ITS-TPC bracket or ITS-Rails
        t0 = 180.*((Double_t)i);
        for(j=-kconNScrewM5by12/2;j<=kconNScrewM5by12/2;j++)if(j!=0){
                    //screws per ITS-TPC brkt
            t = t0 + 5.0*((Double_t)j);
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vD,++nCopyvD,tran);
        } // end or j
        for(j=-kconNPinO6/2;j<=kconNPinO6/2;j++){ // pins per ITS-TPC bracket
            t = t0 + 3.0*((Double_t)j);
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vE,++nCopyvE,tran);
        } // end or j
        t0 = (96.5+187.*((Double_t)i));
        for(j=0;j<kconNRailScrews;j++){ // screws per ITS-rail bracket
            t = t0+da[j];
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vD,++nCopyvD,tran);
        } // end or j
        t0 = (91.5+184.*((Double_t)i));
        for(j=-kconNRailPins/2;j<=kconNRailPins/2;j++)if(j!=0){ 
             // pins per ITS-rail bracket
            t = t0+(7.0*((Double_t)j));
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vE,++nCopyvE,tran);
        } // end or j
    } // end for i
    for(i=0;i<kconNmounts;i++){ 
                // mounting points for SPD-cone+Beam-pipe support
        t0 = (45.0+((Double_t)i)*360./((Double_t)kconNmounts));
        for(j=-1;j<=1;j++)if(j!=0){ // 2 screws per bracket
            t = t0+((Double_t)j)*0.5*kconMountPhi0;
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vD,++nCopyvD,tran);
        } // end for j
        for(j=0;j<1;j++){ // 1 pin per bracket
            t = t0;
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vE,++nCopyvE,tran);
        } // end for j
    } // end for i
    if(GetDebug()){
        vA->PrintNodes();
        vB->PrintNodes();
        vC->PrintNodes();
        vD->PrintNodes();
        vE->PrintNodes();
        vF->PrintNodes();
        vQ->PrintNodes();
        vR->PrintNodes();
        vS->PrintNodes();
        vT->PrintNodes();
    } // end if
}

//______________________________________________________________________
void AliITSv11GeometrySupport::ServicesCableSupport(TGeoVolume *moth){
    // Define the detail ITS cable support trays on both the RB24 and 
    // RB26 sides..
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
    // Based on the Drawings SSup_201A.jpg unless otherwise stated, 
    // Volumes A..., 
    TGeoMedium *medSUPcf    = 0; // SUP support cone Carbon Fiber materal nbr.
    TGeoMedium *medSUPfs    = 0; // SUP support cone inserto stesalite 4411w.
    TGeoMedium *medSUPfo    = 0; // SUP support cone foam, Rohacell 50A.
    TGeoMedium *medSUPss    = 0; // SUP support cone screw material,Stainless
    TGeoMedium *medSUPair   = 0; // SUP support cone Air
    TGeoMedium *medSUPal    = 0; // SUP support cone SDD mounting bracket Al
    TGeoMedium *medSUPwater = 0; // SUP support cone Water
    TGeoManager *mgr = gGeoManager;
    medSUPcf    = mgr->GetMedium("ITSssdCarbonFiber");
    medSUPfs    = mgr->GetMedium("ITSssdStaselite4411w");
    medSUPfo    = mgr->GetMedium("ITSssdRohacell50A");
    medSUPss    = mgr->GetMedium("ITSssdStainlessSteal");
    medSUPair   = mgr->GetMedium("ITSssdAir");
    medSUPal    = mgr->GetMedium("ITSssdAl");
    medSUPwater = mgr->GetMedium("ITSssdWater");
    //
    Int_t i,j;
    Double_t x,y,z,t,t0,dt,di,r;
    // RB 24 side
    const Double_t kfrm24Z0           = 900*fgkmm;//SSup_203A.jpg
    const Double_t kfrm24Thss         = 5.0*fgkmm;
    const Double_t kfrm24Rss          = 444.5*fgkmm-kfrm24Thss; //SSup_204A.jpg
    const Double_t kfrm24Width        = 10.0*fgkmm;
    const Double_t kfrm24Hight        = 10.0*fgkmm;
    const Double_t kfrm24Phi0         = 15.2*fgkDegree; // SSup_602A.jpg
    const Double_t kfrm24Phi1         = (90.0-7.6)*fgkDegree; // SSup_802A.jpg
    const Double_t kfrm24ZssSection   = (415.0-10.0)*fgkmm;
    const Int_t    kfrm24NZsections   = 4;
    const Int_t    kfrm24NPhiSections = 4;
    const Int_t    kfrm24NPhi         = 4;
    //
    TGeoTubeSeg *sM24,*sA24;
    TGeoBBox *sB24;
    sM24 = new TGeoTubeSeg("ITS sup Cable tray support frame mother volume "
                           "M24",kfrm24Rss,kfrm24Rss+kfrm24Thss,
                           0.5*(4.*kfrm24ZssSection+5*kfrm24Width),
                           kfrm24Phi0,kfrm24Phi1);
    sA24 = new TGeoTubeSeg("ITS sup Cable tray support frame radial section "
                           "A24",kfrm24Rss,kfrm24Rss+kfrm24Thss,
                           0.5*kfrm24Width,kfrm24Phi0,kfrm24Phi1);
    sB24 = new TGeoBBox("ITS sup Cable tray support frame Z section B24",
                        0.5*kfrm24Thss,0.5*kfrm24Hight,0.5*kfrm24ZssSection);
    if(GetDebug()){
        sM24->InspectShape();
        sA24->InspectShape();
        sB24->InspectShape();
    } // end if GetDebug()
    TGeoVolume *vA24,*vB24,*vM24;
    TGeoTranslation *tran;
    TGeoRotation    *rot;
    TGeoCombiTrans  *tranrot;
    //
    vA24 = new TGeoVolume("ITSsupFrameA24",sA24,medSUPss);
    vA24->SetVisibility(kTRUE);
    vA24->SetLineColor(1); // black
    vA24->SetLineWidth(1);
    vA24->SetFillColor(vA24->GetLineColor());
    vA24->SetFillStyle(4000); // 0% transparent
    vB24 = new TGeoVolume("ITSsupFrameB24",sB24,medSUPss);
    vB24->SetVisibility(kTRUE);
    vB24->SetLineColor(1); // black
    vB24->SetLineWidth(1);
    vB24->SetFillColor(vB24->GetLineColor());
    vB24->SetFillStyle(4000); // 0% transparent
    vM24 = new TGeoVolume("ITSsupFrameM24",sM24,medSUPair);
    vM24->SetVisibility(kTRUE);
    vM24->SetLineColor(7); // light blue
    vM24->SetLineWidth(1);
    vM24->SetFillColor(vM24->GetLineColor());
    vM24->SetFillStyle(4090); // 90% transparent
    //
    Int_t ncopyA24=1,ncopyB24=1;
    t0 = kfrm24Phi0;
    dt = (kfrm24Phi1-kfrm24Phi0)/((Double_t)kfrm24NPhiSections);
    for(i=0;i<=kfrm24NZsections;i++){
        di = (Double_t) i;
        z = -sM24->GetDz()+sA24->GetDz() + di*(kfrm24ZssSection+kfrm24Width);
        tran = new TGeoTranslation("",0.0,0.0,z);
        vM24->AddNode(vA24,ncopyA24++,tran);
        r = kfrm24Rss+sB24->GetDX();
        z = z + sA24->GetDz()+sB24->GetDZ();
       if(i<kfrm24NZsections) for(j=0;j<=kfrm24NPhiSections;j++){
            t = t0 + ((Double_t)j)*dt;
            rot = new TGeoRotation("",0.0,0.0,t);
            y = r*SinD(t);
            x = r*CosD(t);
            tranrot = new TGeoCombiTrans("",x,y,z,rot);
            delete rot;// rot not explicity used in AddNode functions.
            vM24->AddNode(vB24,ncopyB24++,tranrot);
        } // end for j
    } // end for i
    tran = new TGeoTranslation("",0.0,0.0,kfrm24Z0+sM24->GetDz());
    moth->AddNode(vM24,1,tran);
    for(i=1;i<kfrm24NPhi;i++){
        di = (Double_t) i;
        rot = new TGeoRotation("",0.0,0.0,90.0*di);
        tranrot = new TGeoCombiTrans("",0.0,0.0,kfrm24Z0+sM24->GetDz(),rot);
        delete rot;// rot not explicity used in AddNode functions.
        moth->AddNode(vM24,i+1,tranrot);
    } // end for i
    if(GetDebug()){
        vA24->PrintNodes();
        vB24->PrintNodes();
        vM24->PrintNodes();
    } // end if
    // Cable support tray 
    // Material is Aluminum
    //const Double_t kcsb24RSin       = TMath::Max(kfrm24Rss,444.5*fgkmm);
                                           // SSup_204A & SSup_206A
    //const Double_t kcb24RSAirout   = 459.5*fgkmm; // SSup_204A & SSup_206A
    //const Double_t kcb24RSout      = 494.5*fgkmm; // SSup_206A & SSup_204A
    //const Double_t kcb24RSPPout    = 550.0*fgkmm; // SSup_206A
    const Double_t kcb24LSPP       = 350.0*fgkmm; // SSup_202A
    const Double_t kcb24LS         = (2693.0-900.0)*fgkmm;//SSup_205A&SSup_207A
    const Double_t kcb24ThSwall    = 1.0*fgkmm; // SSup_209A & SSup_210A
    const Double_t kcb24WbS        = 42.0*fgkmm; // SSup_209A & SSup_210A
    //const Double_t kcb24WtS        = 46.9*fgkmm; // SSup_209A & SSup_210A
    const Double_t kcb24WcapS      = 50.0*fgkmm; // SSup_209A & SSup_210A
    //const Double_t kcb24WdS   = 41.0*fgkmm; //SSup_209A ? should be 41.469387
    const Double_t kcb24HS         = 50.0*fgkmm; // SSup_209A & SSup_210A
    const Double_t kcb24OutDcoolTub= 12.0*fgkmm; // SSup_209A
    const Double_t kcb24InDcoolTub = 10.0*fgkmm; // SSup_209A
    const Double_t kcbBlkNozInDS   = 6.0*fgkmm; // SSup_209A
    // The following are deduced or guessed at
    //const Double_t kcb24LtopLipS   = 6.0*fgkmm; // Guessed at.
    //const Double_t kcb24LdLipS     = 6.0*fgkmm; // Guessed at.
    //const Double_t kcb24HdS        = kcb24OutDcoolTub; //
    const Double_t kcb24BlkNozZS   = 6.0*fgkmm; // Guessed at.
    // Simplifided exterior shape. The side wall size is 2.5*thicker than
    // it should be (due to simplification).
    TGeoArb8 *sC24,*sD24,*sF24,*sH24;
    TGeoTube *sE24,*sG24;
    //
    sC24 = new TGeoArb8("ITS Sup Cable Tray Element C24",0.5*kcb24LS);
    sC24->SetVertex(0,-0.5*kcb24WcapS,kcb24HS+kcb24ThSwall);
    sC24->SetVertex(1,+0.5*kcb24WcapS,kcb24HS+kcb24ThSwall);
    sC24->SetVertex(2,+0.5*kcb24WbS,0.0);
    sC24->SetVertex(3,-0.5*kcb24WbS,0.0);
    sC24->SetVertex(4,-0.5*kcb24WcapS,kcb24HS+kcb24ThSwall);
    sC24->SetVertex(5,+0.5*kcb24WcapS,kcb24HS+kcb24ThSwall);
    sC24->SetVertex(6,+0.5*kcb24WbS,0.0);
    sC24->SetVertex(7,-0.5*kcb24WbS,0.0);
    sD24 = new TGeoArb8("ITS Sup Cable Tray lower Element D24",0.5*kcb24LS);
    // Because of question about the value of WdS24, compute what it
    // should be assuming cooling tube fixes hight of volume.
    x = kcb24OutDcoolTub*(0.5*kcb24WcapS-0.5*kcb24WbS-kcb24ThSwall)/
                                              (kcb24HS-kcb24ThSwall);
    sD24->SetVertex(0,-x,kcb24OutDcoolTub+kcb24ThSwall);
    sD24->SetVertex(1,+x,kcb24OutDcoolTub+kcb24ThSwall);
    sD24->SetVertex(2,+0.5*kcb24WbS-kcb24ThSwall,kcb24ThSwall);
    sD24->SetVertex(3,-0.5*kcb24WbS+kcb24ThSwall,kcb24ThSwall);
    sD24->SetVertex(4,-x,kcb24OutDcoolTub+kcb24ThSwall);
    sD24->SetVertex(5,+x,kcb24OutDcoolTub+kcb24ThSwall);
    sD24->SetVertex(6,+0.5*kcb24WbS-kcb24ThSwall,kcb24ThSwall);
    sD24->SetVertex(7,-0.5*kcb24WbS+kcb24ThSwall,kcb24ThSwall);
    sE24 = new TGeoTube("ITS Sup Cooling Tube E24",0.5*kcb24InDcoolTub,
                        0.5*kcb24OutDcoolTub,0.5*kcb24LS-kcb24BlkNozZS);
    sF24 = new TGeoArb8("ITS Sup Cable Tray lower Element block F24",
                        0.5*kcb24BlkNozZS);
    for(i=0;i<8;i++) sF24->SetVertex(i,sD24->GetVertices()[i*2+0],
                                      sD24->GetVertices()[i*2+1]); //
    sG24 = new TGeoTube("ITS Sup Cooling Tube hole in block G24",
                        0.0,0.5*kcbBlkNozInDS,0.5*kcb24BlkNozZS);
    sH24 = new TGeoArb8("ITS Sup Cable Tray upper Element H24",
                        0.5*(kcb24LS- kcb24LSPP));
    sH24->SetVertex(0,sC24->GetVertices()[0*2+0]+2.*kcb24ThSwall,
                     sC24->GetVertices()[0*2+1]-kcb24ThSwall);
    sH24->SetVertex(1,sC24->GetVertices()[1*2+0]-2.*kcb24ThSwall,
                     sC24->GetVertices()[1*2+1]-kcb24ThSwall);
    sH24->SetVertex(2,sD24->GetVertices()[1*2+0]-kcb24ThSwall,
                     sD24->GetVertices()[1*2+1]+kcb24ThSwall);
    sH24->SetVertex(3,sD24->GetVertices()[0*2+0]+kcb24ThSwall,
                     sD24->GetVertices()[0*2+1]+kcb24ThSwall);
    for(i=4;i<8;i++) sH24->SetVertex(i,sH24->GetVertices()[(i-4)*2+0],
                                      sH24->GetVertices()[(i-4)*2+1]); //
    if(GetDebug()){
        sC24->InspectShape();
        sD24->InspectShape();
        sF24->InspectShape();
        sH24->InspectShape();
        sE24->InspectShape();
        sG24->InspectShape();
    } // end if GetDebug()
    TGeoVolume *vC24,*vD24,*vE24,*vF24,*vGa24,*vGw24,*vH24;
    //
    vC24 = new TGeoVolume("ITSsupCableTrayC24",sC24,medSUPal);
    vC24->SetVisibility(kTRUE);
    vC24->SetLineColor(6); //
    vC24->SetLineWidth(1);
    vC24->SetFillColor(vC24->GetLineColor());
    vC24->SetFillStyle(4000); // 0% transparent
    vD24 = new TGeoVolume("ITSsupCableTrayLowerD24",sD24,medSUPair);
    vD24->SetVisibility(kTRUE);
    vD24->SetLineColor(6); //
    vD24->SetLineWidth(1);
    vD24->SetFillColor(vD24->GetLineColor());
    vD24->SetFillStyle(4000); // 0% transparent
    vE24 = new TGeoVolume("ITSsupCableTrayCoolTubeE24",sE24,medSUPss);
    vE24->SetVisibility(kTRUE);
    vE24->SetLineColor(6); //
    vE24->SetLineWidth(1);
    vE24->SetFillColor(vE24->GetLineColor());
    vE24->SetFillStyle(4000); // 0% transparent
    vF24 = new TGeoVolume("ITSsupCableTrayBlockF24",sF24,medSUPal);
    vF24->SetVisibility(kTRUE);
    vF24->SetLineColor(6); //
    vF24->SetLineWidth(1);
    vF24->SetFillColor(vF24->GetLineColor());
    vF24->SetFillStyle(4000); // 0% transparent
    vGw24 = new TGeoVolume("ITSsupCableTrayCoolantWaterG24",sG24,medSUPwater);
    vGw24->SetVisibility(kTRUE);
    vGw24->SetLineColor(6); //
    vGw24->SetLineWidth(1);
    vGw24->SetFillColor(vGw24->GetLineColor());
    vGw24->SetFillStyle(4000); // 0% transparent
    vGa24 = new TGeoVolume("ITSsupCableTrayCoolantAirG24",sG24,medSUPair);
    vGa24->SetVisibility(kTRUE);
    vGa24->SetLineColor(6); //
    vGa24->SetLineWidth(1);
    vGa24->SetFillColor(vGa24->GetLineColor());
    vGa24->SetFillStyle(4000); // 0% transparent
    vH24 = new TGeoVolume("ITSsupCableTrayUpperC24",sH24,medSUPair);
    vH24->SetVisibility(kTRUE);
    vH24->SetLineColor(6); //
    vH24->SetLineWidth(1);
    vH24->SetFillColor(vH24->GetLineColor());
    vH24->SetFillStyle(4000); // 0% transparent
    //
    tran = new TGeoTranslation("",-kcb24OutDcoolTub,
                               kcb24OutDcoolTub+kcb24ThSwall,0.0);
    vF24->AddNode(vGw24,1,tran);
    vD24->AddNode(vE24,1,tran);
    tran = new TGeoTranslation("",0.0,kcb24OutDcoolTub+kcb24ThSwall,0.0);
    vF24->AddNode(vGw24,2,tran);
    vD24->AddNode(vE24,2,tran);
    tran = new TGeoTranslation("",+kcb24OutDcoolTub,
                               kcb24OutDcoolTub+kcb24ThSwall,0.0);
    vF24->AddNode(vGw24,3,tran);
    vD24->AddNode(vE24,3,tran);
    tran = new TGeoTranslation("",0.0,0.0,0.5*kcb24LS-0.5*kcb24BlkNozZS);
    vD24->AddNode(vF24,1,tran);
    tran = new TGeoTranslation("",0.0,0.0,-(0.5*kcb24LS-0.5*kcb24BlkNozZS));
    vD24->AddNode(vF24,2,tran);
    vC24->AddNode(vD24,1,0);
    vC24->AddNode(vH24,1,0);
    if(GetDebug()){
        vC24->PrintNodes();
        vD24->PrintNodes();
        vE24->PrintNodes();
        vF24->PrintNodes();
        vGa24->PrintNodes();
        vGw24->PrintNodes();
        vH24->PrintNodes();
    } // end if GetDebug()
    //==================================================================
    //
    // RB 26 side
    const Double_t kfrm26Z0           = -900*fgkmm;//SSup_203A.jpg
    const Double_t kfrm26Thss         = 5.0*fgkmm;
    const Double_t kfrm26R0ss         = 444.5*fgkmm-kfrm26Thss; //SSup_204A.jpg
    const Double_t kfrm26R1ss         = 601.6*fgkmm-kfrm26Thss; //SSup_208A.jpg
    const Double_t kfrm26Width        = 10.0*fgkmm;
    //const Double_t kfrm26Hight       = 10.0*fgkmm;
    const Double_t kfrm26Phi0         = 15.2*fgkDegree; // SSup_602A.jpg
    const Double_t kfrm26Phi1         = (90.0-7.6)*fgkDegree; // SSup_802A.jpg
    const Double_t kfrm26ZssSection   = (415.0-10.0)*fgkmm;
    const Int_t    kfrm26NZsections   = 4;
    const Int_t    kfrm26NPhiSections = 4;
    const Int_t    kfrm26NPhi         = 4;
    TGeoConeSeg *sA26[kfrm26NZsections+1],*sM26;//Cylinderial support structure
    TGeoArb8     *sB26; // Cylinderial support structure
    Char_t name[100];
    Double_t r1,r2,m;

    sM26 = new TGeoConeSeg("ITS sup Cable tray support frame mother volume "
                          "M26",0.5*(4.*kfrm26ZssSection+5*kfrm26Width),
                          kfrm26R1ss,kfrm26R1ss+kfrm26Thss,
                          kfrm26R0ss,kfrm26R0ss+kfrm26Thss,
                          kfrm26Phi0,kfrm26Phi1);
    m = -((kfrm26R1ss-kfrm26R0ss)/
         (((Double_t)kfrm26NZsections)*(kfrm26ZssSection+kfrm26Width)));
    for(i=0;i<kfrm26NZsections+1;i++){
        di = ((Double_t) i)*(kfrm26ZssSection+kfrm26Width);
        sprintf(name,
                "ITS sup Cable tray support frame radial section A26[%d]",i);
        r1 = kfrm26R1ss+m*di;
        r2 = kfrm26R1ss+m*(di+kfrm26Width);
        sA26[i] = new TGeoConeSeg(name,0.5*kfrm26Width,r2,r2+kfrm26Thss,
                                 r1,r1+kfrm26Thss,kfrm26Phi0,kfrm26Phi1);
    } // end for i
    sB26 = new TGeoArb8("ITS sup Cable tray support frame Z section B26",
                       0.5*kfrm26ZssSection);
    r = 0.25*(sA26[0]->GetRmax1()+sA26[0]->GetRmin1()+
              sA26[1]->GetRmax2()+sA26[1]->GetRmin2());
    sB26->SetVertex(0,sA26[0]->GetRmax2()-r,+0.5*kfrm26Width);
    sB26->SetVertex(1,sA26[0]->GetRmax2()-r,-0.5*kfrm26Width);
    sB26->SetVertex(2,sA26[0]->GetRmin2()-r,-0.5*kfrm26Width);
    sB26->SetVertex(3,sA26[0]->GetRmin2()-r,+0.5*kfrm26Width);
    sB26->SetVertex(4,sA26[1]->GetRmax1()-r,+0.5*kfrm26Width);
    sB26->SetVertex(5,sA26[1]->GetRmax1()-r,-0.5*kfrm26Width);
    sB26->SetVertex(6,sA26[1]->GetRmin1()-r,-0.5*kfrm26Width);
    sB26->SetVertex(7,sA26[1]->GetRmin1()-r,+0.5*kfrm26Width);
    if(GetDebug()){
        for(i=0;i<kfrm26NZsections+1;i++) sA26[i]->InspectShape();
        sM26->InspectShape();
        sB26->InspectShape();
    } // end if GetDebug()
    //
    TGeoVolume *vA26[kfrm26NZsections+1],*vB26,*vM26;
    //
    for(i=0;i<kfrm26NZsections+1;i++){
        sprintf(name,"ITSsupFrameA26[%d]",i);
        vA26[i] = new TGeoVolume(name,sA26[i],medSUPss);
        vA26[i]->SetVisibility(kTRUE);
        vA26[i]->SetLineColor(1); // black
        vA26[i]->SetLineWidth(1);
        vA26[i]->SetFillColor(vA26[i]->GetLineColor());
        vA26[i]->SetFillStyle(4000); // 0% transparent
    } // end for i
    vB26 = new TGeoVolume("ITSsupFrameB26",sB26,medSUPss);
    vB26->SetVisibility(kTRUE);
    vB26->SetLineColor(1); // black
    vB26->SetLineWidth(1);
    vB26->SetFillColor(vB26->GetLineColor());
    vB26->SetFillStyle(4000); // 0% transparent
    vM26 = new TGeoVolume("ITSsupFrameM26",sM26,medSUPair);
    vM26->SetVisibility(kTRUE);
    vM26->SetLineColor(7); // light blue
    vM26->SetLineWidth(1);
    vM26->SetFillColor(vM26->GetLineColor());
    vM26->SetFillStyle(4090); // 90% transparent
    //
    Int_t ncopyB26=1;
    t0 = kfrm26Phi0;
    dt = (kfrm26Phi1-kfrm26Phi0)/((Double_t)kfrm26NPhiSections);
    for(i=0;i<=kfrm26NZsections;i++){
        di = ((Double_t) i)*(kfrm26ZssSection+kfrm26Width);
        z = -sM26->GetDz()+sA26[i]->GetDz() + di;
        tran = new TGeoTranslation("",0.0,0.0,z);
        vM26->AddNode(vA26[i],1,tran);
        z = z+sB26->GetDz();
        if(i<kfrm26NZsections)for(j=0;j<=kfrm26NPhiSections;j++){
            r = 0.25*(sA26[i]->GetRmax1()+sA26[i]->GetRmin1()+
                      sA26[i+1]->GetRmax2()+sA26[i+1]->GetRmin2());
            t = t0 + ((Double_t)j)*dt;
            rot = new TGeoRotation("",0.0,0.0,t);
            y = r*SinD(t);
            x = r*CosD(t);
            tranrot = new TGeoCombiTrans("",x,y,z,rot);
            delete rot; // rot not explicity used in AddNode functions.
            vM26->AddNode(vB26,ncopyB26++,tranrot);
        } // end for j
    } // end for i
    tran = new TGeoTranslation("",0.0,0.0,kfrm26Z0-sM26->GetDz());
    moth->AddNode(vM26,1,tran);
    for(i=1;i<kfrm26NPhi;i++){
        rot = new TGeoRotation("",0.0,0.0,90.0*((Double_t)i));
        tranrot = new TGeoCombiTrans(*tran,*rot);
        delete rot; // rot not explicity used in AddNode functions.
        moth->AddNode(vM26,i+1,tranrot);
    } // end for i
    if(GetDebug()){
        for(i=0;i<kfrm26NZsections+1;i++) vA26[i]->PrintNodes();
        vB26->PrintNodes();
        vM26->PrintNodes();
    } // end if
}
