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
#include "AliITSv11GeometrySupport.h"

ClassImp(AliITSv11GeometrySupport)

#define SQ(A) (A)*(A)

//______________________________________________________________________
void AliITSv11GeometrySupport::SPDCone(TGeoVolume *Moth){
    // Define the detail SPD support cone geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.

    SPDThermalSheald(Moth);
}
//______________________________________________________________________
void AliITSv11GeometrySupport::SPDThermalSheald(TGeoVolume *Moth){
    // Define the detail SPD Thermal Sheld geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
    // From ALICE-Thermal Screen (SPD) "Cylinder" file thermal-screen2_a3.ps
    // Volumes A1,A2,A2,Ah1,Ah2,Ah3, and B1,B2,B3,Bh1,Bh2,Bh3;
    // "CONE TRANSITION" file thermal-screen1_a3.ps Volumes C1,C2,C3,Ch1,Ch2,
    // Ch3; "FLANGE" file thermal-screen4_a3.ps Volumes D,Ds,Dw,Dws; and 
    // "HALF ASSEMBLY" file thermal-screen3_a3.ps. This object, both halfs,
    // are incased inside of a single minimum sized mother volume called M,
    // which is a union of two parts M1 and 4 copies of M2.
    const Double_t TSCarbonFiberThA = 0.03*kmm; // 
    //const Double_t TSCarbonFiberThB = 0.10*kmm; //
    const Double_t TSCLengthB  = 50.0*kmm; //
    const Double_t TSCLengthA  = 900.0*kmm-2.0*TSCLengthB; //
    const Double_t TSCLengthC  = 290.0*kmm; //
    const Double_t TSCLengthD  = 15.0*kmm; //
    const Double_t TSCAngle    = 36.0*kDegree;//Rep. angle of cent. accordin
    const Double_t TSCRoutA    = 99.255*kmm; // Outer radii
    const Double_t TSCRinA     = 81.475*kmm; // Iner radii
    const Double_t TSCRoutB    = 99.955*kmm; // Outer radii
    const Double_t TSCRinB     = 80.775*kmm; // Iner radii
    const Double_t TSCRoutCp   = 390.0*kmm;  // Outer radii
    const Double_t TSCRinCp    = 373.0*kmm;  // Iner radii
    Double_t TSCRoutC,TSCRinC; // values need to be calculated
    const Double_t TSCRwingD   = 492.5*kmm;  // Outer radii
    const Double_t TSCRoutD    = 0.5*840.*kmm;// Outer radii
    const Double_t TSCRinD     = 373.0*kmm;  // Iner radii
    const Double_t TSCAngleDD  = 60.*kmm/TSCRwingD/kRadian;// angular wing
                                                    // width of fill material
    const Double_t TSCAngleDDs = (60.*kmm-2.*TSCarbonFiberThA)/
                                                  TSCRwingD/kRadian;
    const Double_t TSCAngleD0  = 45.*kDegree;//Strting angle of wing
    const Double_t TSCoutSA    = 24.372*kmm; // The other one Calculated
    const Double_t TSCinLA     = 31.674*kmm; // The ohter one Calculated
    const Double_t TSCoutSB    = 24.596*kmm; // The other one Calculated
    const Double_t TSCinLB     = 31.453*kmm; // The ohter one Calculated
    const Double_t TSCoutSC    = 148.831*kmm;// The other one Calculated
    const Double_t TSCinLC     = 90.915*kmm; // The ohter one Calculated
    Int_t i,k;
    Double_t th;
    Double_t xo[7],yo[7],xi[7],yi[7];
    Double_t xbo[7],ybo[7],xbi[7],ybi[7];
    Double_t xco[7],yco[7],xci[7],yci[7];
    TGeoArb8 *A1,*A2,*A3,*Ah1,*Ah2,*Ah3,*B1,*B2,*B3,*Bh1,*Bh2,*Bh3;
    TGeoArb8 *C1,*C2,*C3,*Ch1,*Ch2,*Ch3;
    TGeoTube  *D,*Ds;
    TGeoTubeSeg *Dw,*Dws;
    TGeoCompositeShape *M;
    TGeoRotation *rot;
    TGeoTranslation *tranb,*tranbm,*tranc;
    TGeoTranslation *tranITSspdShealdVVt0;
    TGeoCombiTrans *rotITSspdShealdVVt1,*rotITSspdShealdVVt2;
    TGeoCombiTrans *rotITSspdShealdVVt3;
    TGeoMedium *SPDcf  = 0; // SPD support cone Carbon Fiber materal number.
    TGeoMedium *SPDfs  = 0; // SPD support cone inserto stesalite 4411w.
    TGeoMedium *SPDfo  = 0; // SPD support cone foam, Rohacell 50A.
    TGeoMedium *SPDss  = 0; // SPD support cone screw material,Stainless steal
    TGeoMedium *SPDair = 0; // SPD support cone Air
    //TGeoMedium *SPDal  = 0; // SPD support cone SDD mounting bracket Al

    TSCRoutC = TMath::Sqrt(TSCRoutCp*TSCRoutCp-0.25*TSCoutSC*TSCoutSC);
    TSCRinC  = TMath::Sqrt(TSCRinCp *TSCRinCp -0.25*TSCinLC *TSCinLC );
    A1  = new TGeoArb8("ITS SPD Therm Screen Clyinder A1",0.5*TSCLengthA);
    A2  = new TGeoArb8("ITS SPD Therm Screen Clyinder A2",0.5*TSCLengthA);
    A3  = new TGeoArb8("ITS SPD Therm Screen Clyinder A3",0.5*TSCLengthA);
    Ah1 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ah1",0.5*TSCLengthA);
    Ah2 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ah2",0.5*TSCLengthA);
    Ah3 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ah3",0.5*TSCLengthA);
    B1  = new TGeoArb8("ITS SPD Therm Screen Clyinder B1",0.5*TSCLengthB);
    B2  = new TGeoArb8("ITS SPD Therm Screen Clyinder B2",0.5*TSCLengthB);
    B3  = new TGeoArb8("ITS SPD Therm Screen Clyinder B3",0.5*TSCLengthB);
    Bh1 = new TGeoArb8("ITS SPD Therm Screen Cylinder Bh1",0.5*TSCLengthB);
    Bh2 = new TGeoArb8("ITS SPD Therm Screen Cylinder Bh2",0.5*TSCLengthB);
    Bh3 = new TGeoArb8("ITS SPD Therm Screen Cylinder Bh3",0.5*TSCLengthB);
    C1  = new TGeoArb8("ITS SPD Therm Screen Clyinder C1",0.5*TSCLengthC);
    C2  = new TGeoArb8("ITS SPD Therm Screen Clyinder C2",0.5*TSCLengthC);
    C3  = new TGeoArb8("ITS SPD Therm Screen Clyinder C3",0.5*TSCLengthC);
    Ch1 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ch1",0.5*TSCLengthC);
    Ch2 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ch2",0.5*TSCLengthC);
    Ch3 = new TGeoArb8("ITS SPD Therm Screen Cylinder Ch3",0.5*TSCLengthC);
    D = new TGeoTube("ITS SPD Therm Screen Flange D",TSCRinD,TSCRoutD,
                    0.5*TSCLengthD);
    Ds = new TGeoTube("ITS SPD Therm Screen Flange fill Ds",
                      TSCRinD+TSCarbonFiberThA,TSCRoutD-TSCarbonFiberThA,
                      0.5*TSCLengthD);
    printTube(D);
    printTube(Ds);
    Dw = new TGeoTubeSeg("ITS SPD Therm Screen Flange Wing Dw",
                         TSCRoutD,TSCRwingD ,0.5*TSCLengthD,
                         TSCAngleD0-0.5*TSCAngleDD,TSCAngleD0+0.5*TSCAngleDD);
    Dws = new TGeoTubeSeg("ITS SPD Therm Screen Flange Wing Fill Ds",
                          TSCRoutD,TSCRwingD-TSCarbonFiberThA,
                          0.5*TSCLengthD,TSCAngleD0-0.5*TSCAngleDDs,
                          TSCAngleD0+0.5*TSCAngleDDs);
    printTubeSeg(Dw);
    printTubeSeg(Dws);
    k = 0;
    for(i=-1;i<2;i++){
      th = ((Double_t)(i+1))*TSCAngle*kRadian;
      xo[k] = TSCRoutA*TMath::Sin(th) - 0.5*TSCoutSA*TMath::Cos(th);
      yo[k] = TSCRoutA*TMath::Cos(th) + 0.5*TSCoutSA*TMath::Sin(th);
      xi[k] = TSCRinA *TMath::Sin(th) - 0.5*TSCinLA *TMath::Cos(th);
      yi[k] = TSCRinA *TMath::Cos(th) + 0.5*TSCinLA *TMath::Sin(th);
      xbo[k] = TSCRoutB*TMath::Sin(th) - 0.5*TSCoutSB*TMath::Cos(th);
      ybo[k] = TSCRoutB*TMath::Cos(th) + 0.5*TSCoutSB*TMath::Sin(th);
      xbi[k] = TSCRinB *TMath::Sin(th) - 0.5*TSCinLB *TMath::Cos(th);
      ybi[k] = TSCRinB *TMath::Cos(th) + 0.5*TSCinLB *TMath::Sin(th);
      xco[k] = TSCRoutC*TMath::Sin(th) - 0.5*TSCoutSC*TMath::Cos(th);
      yco[k] = TSCRoutC*TMath::Cos(th) + 0.5*TSCoutSC*TMath::Sin(th);
      xci[k] = TSCRinC *TMath::Sin(th) - 0.5*TSCinLC *TMath::Cos(th);
      yci[k] = TSCRinC *TMath::Cos(th) + 0.5*TSCinLC *TMath::Sin(th);
      k++;
      xo[k] = TSCRoutA*TMath::Sin(th) + 0.5*TSCoutSA*TMath::Cos(th);
      yo[k] = TSCRoutA*TMath::Cos(th) - 0.5*TSCoutSA*TMath::Sin(th);
      xi[k] = TSCRinA *TMath::Sin(th) + 0.5*TSCinLA *TMath::Cos(th);
      yi[k] = TSCRinA *TMath::Cos(th) - 0.5*TSCinLA *TMath::Sin(th);
      xbo[k] = TSCRoutB*TMath::Sin(th) + 0.5*TSCoutSB*TMath::Cos(th);
      ybo[k] = TSCRoutB*TMath::Cos(th) - 0.5*TSCoutSB*TMath::Sin(th);
      xbi[k] = TSCRinB *TMath::Sin(th) + 0.5*TSCinLB *TMath::Cos(th);
      ybi[k] = TSCRinB *TMath::Cos(th) - 0.5*TSCinLB *TMath::Sin(th);
      xco[k] = TSCRoutC*TMath::Sin(th) + 0.5*TSCoutSC*TMath::Cos(th);
      yco[k] = TSCRoutC*TMath::Cos(th) - 0.5*TSCoutSC*TMath::Sin(th);
      xci[k] = TSCRinC *TMath::Sin(th) + 0.5*TSCinLC *TMath::Cos(th);
      yci[k] = TSCRinC *TMath::Cos(th) - 0.5*TSCinLC *TMath::Sin(th);
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
    cout.precision(4);
    cout.width(7);
    cout <<"i     \t  xo  yo    \t  xi yi     \t  xbo ybo   \t   xbi ybi  "
        "\t   xco yco   \t   xci yxi"<<endl;
    for(i=0;i<7;i++){
        cout << i <<"\t"<<xo[i]<<","<<yo[i];
        cout      <<"\t"<<xi[i]<<","<<yi[i];
        cout      <<"\t"<<xbo[i]<<","<<ybo[i];
        cout      <<"\t"<<xbi[i]<<","<<ybi[i];
        cout      <<"\t"<<xco[i]<<","<<yco[i];
        cout      <<"\t"<<xci[i]<<","<<yci[i];
        cout<<endl;}
    } // end if GetDebug()
    //+++++++++++++++++++++++++
    A1->SetVertex(0,xo[0],yo[0]);
    A1->SetVertex(1,xo[1],yo[1]);
    A1->SetVertex(2,xi[1],yi[1]);
    A1->SetVertex(3,xi[0],yi[0]);
    //
    A2->SetVertex(0,xo[1],yo[1]);
    A2->SetVertex(1,xo[2],yo[2]);
    A2->SetVertex(2,xi[2],yi[2]);
    A2->SetVertex(3,xi[1],yi[1]);
    //
    A3->SetVertex(0,xo[5],yo[5]);
    A3->SetVertex(1,xo[6],yo[6]);
    A3->SetVertex(2,xi[6],yi[6]);
    A3->SetVertex(3,xi[5],yi[5]);
    //--------------------------
    B1->SetVertex(0,xbo[0],ybo[0]);
    B1->SetVertex(1,xbo[1],ybo[1]);
    B1->SetVertex(2,xbi[1],ybi[1]);
    B1->SetVertex(3,xbi[0],ybi[0]);
    //
    B2->SetVertex(0,xbo[1],ybo[1]);
    B2->SetVertex(1,xbo[2],ybo[2]);
    B2->SetVertex(2,xbi[2],ybi[2]);
    B2->SetVertex(3,xbi[1],ybi[1]);
    //
    B3->SetVertex(0,xbo[5],ybo[5]);
    B3->SetVertex(1,xbo[6],ybo[6]);
    B3->SetVertex(2,xbi[6],ybi[6]);
    B3->SetVertex(3,xbi[5],ybi[5]);
    //--------------------------
    C1->SetVertex(0,xco[0],yco[0]);
    C1->SetVertex(1,xco[1],yco[1]);
    C1->SetVertex(2,xci[1],yci[1]);
    C1->SetVertex(3,xci[0],yci[0]);
    //
    C2->SetVertex(0,xco[1],yco[1]);
    C2->SetVertex(1,xco[2],yco[2]);
    C2->SetVertex(2,xci[2],yci[2]);
    C2->SetVertex(3,xci[1],yci[1]);
    //
    C3->SetVertex(0,xco[5],yco[5]);
    C3->SetVertex(1,xco[6],yco[6]);
    C3->SetVertex(2,xci[6],yci[6]);
    C3->SetVertex(3,xci[5],yci[5]);
    // Defining the hole, filled with air
    Double_t p1,c1,x,y,x7[3],y7[3];
    p1 = (xo[0]-xi[0])/(yo[0]-yi[0]);
    c1 = xo[0]+0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(xo[0]-xi[0])+
                                            SQ(yo[0]-yi[0]))/(xo[0]-xi[0]);
    y = TSCRoutA-2.*TSCarbonFiberThA;
    x = p1*(y-yo[0])+c1;
    Ah1->SetVertex(0,x,y);
    Bh1->SetVertex(0,x,y);
    Ch1->SetVertex(4,x,y);
    y = TSCRinA+TSCarbonFiberThA;
    x = p1*(y-yo[0])+c1;
    Ah1->SetVertex(3,x,y);
    Bh1->SetVertex(3,x,y);
    x7[0] = x; y7[0] = y; // vortexing done after last point
    //Ch1->SetVertex(7,x,y);
    p1 = (xo[1]-xi[1])/(yo[1]-yi[1]);
    c1 = xo[1]-0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(xo[1]-xi[1])+
                                            SQ(yo[1]-yi[1]))/(xo[1]-xi[1]);
    y = TSCRoutA-2.*TSCarbonFiberThA;
    x = p1*(y-yo[1])+c1;
    Ah1->SetVertex(1,x,y);
    Bh1->SetVertex(1,x,y);
    Ch1->SetVertex(5,x,y);
    y = TSCRinA+TSCarbonFiberThA;
    x = p1*(y-yo[1])+c1;
    Ah1->SetVertex(2,x,y);
    Bh1->SetVertex(2,x,y);
    Ch1->SetVertex(6,x,y);
    //
    // The easist way to get the points for the hole in volume A2 is to
    // rotate it to the Y axis where the y coordinates are easier to know
    // and then rotate it back.
    Double_t xp,yp,xa,ya,xb,yb;
    th = 0.5*TSCAngle*kRadian;
    xa = TMath::Cos(th)*xo[1]-TMath::Sin(th)*yo[1];
    ya = TMath::Sin(th)*xo[1]+TMath::Cos(th)*yo[1];
    xb = TMath::Cos(th)*xi[1]-TMath::Sin(th)*yi[1];
    yb = TMath::Sin(th)*xi[1]+TMath::Cos(th)*yi[1];
    p1 = (xa-xb)/(ya-yb);
    c1 = xa+0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(xa-xb)+SQ(ya-yb))/(xa-xb);
    y = ya-TSCarbonFiberThA;
    x = p1*(y-ya)+c1;
    xp = TMath::Cos(-th)*x-TMath::Sin(-th)*y;
    yp = TMath::Sin(-th)*x+TMath::Cos(-th)*y;
    Ah2->SetVertex(0,xp,yp);
    Bh2->SetVertex(0,xp,yp);
    Ch2->SetVertex(4,xp,yp);
    y = yb+2.0*TSCarbonFiberThA;
    x = p1*(y-ya)+c1;
    xp = TMath::Cos(-th)*x-TMath::Sin(-th)*y;
    yp = TMath::Sin(-th)*x+TMath::Cos(-th)*y;
    Ah2->SetVertex(3,xp,yp);
    Bh2->SetVertex(3,xp,yp);
    x7[1] = x; y7[1] = y; // vortexing done after last point
    //Ch2->SetVertex(7,xp,yp);
    xa = TMath::Cos(th)*xo[2]-TMath::Sin(th)*yo[2];
    ya = TMath::Sin(th)*xo[2]+TMath::Cos(th)*yo[2];
    xb = TMath::Cos(th)*xi[2]-TMath::Sin(th)*yi[2];
    yb = TMath::Sin(th)*xi[2]+TMath::Cos(th)*yi[2];
    p1 = (xa-xb)/(ya-yb);
    c1 = xa-0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(xa-xb)+SQ(ya-yb))/(xa-xb);
    y = ya-TSCarbonFiberThA;
    x = p1*(y-ya)+c1;
    xp = TMath::Cos(-th)*x-TMath::Sin(-th)*y;
    yp = TMath::Sin(-th)*x+TMath::Cos(-th)*y;
    Ah2->SetVertex(1,xp,yp);
    Bh2->SetVertex(1,xp,yp);
    Ch2->SetVertex(5,xp,yp);
    y = yb+2.0*TSCarbonFiberThA;
    x = p1*(y-ya)+c1;
    xp = TMath::Cos(-th)*x-TMath::Sin(-th)*y;
    yp = TMath::Sin(-th)*x+TMath::Cos(-th)*y;
    Ah2->SetVertex(2,xp,yp);
    Bh2->SetVertex(2,xp,yp);
    Ch2->SetVertex(6,xp,yp);
    //
    p1 = (yo[5]-yi[5])/(xo[5]-xi[5]);
    c1 = yo[5]+0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(yo[5]-yi[5])+
                                            SQ(xo[5]-xi[5]))/(yo[5]-yi[5]);
    x = xo[5]-TSCarbonFiberThA;
    y = p1*(x-xo[5])+c1;
    Ah3->SetVertex(0,x,y);
    Bh3->SetVertex(0,x,y);
    Ch3->SetVertex(4,x,y);
    x = xi[5]+2.0*TSCarbonFiberThA;
    y = p1*(x-xo[5])+c1;
    Ah3->SetVertex(3,x,y);
    Bh3->SetVertex(3,x,y);
    x7[2] = x; y7[2] = y; // vortexing done after last point
    //Ch3->SetVertex(7,x,y);
    y = 2.0*TSCarbonFiberThA;
    x = xo[5]-TSCarbonFiberThA;
    Ah3->SetVertex(1,x,y);
    Bh3->SetVertex(1,x,y);
    Ch3->SetVertex(5,x,y);
    y = 2.0*TSCarbonFiberThA;
    x = xi[5]+2.0*TSCarbonFiberThA;
    Ah3->SetVertex(2,x,y);
    Bh3->SetVertex(2,x,y);
    Ch3->SetVertex(6,x,y);
    //
    for(i=0;i<4;i++){ // define points at +dz
      A1->SetVertex(i+4,(A1->GetVertices())[2*i],(A1->GetVertices())[1+2*i]);
      A2->SetVertex(i+4,(A2->GetVertices())[2*i],(A2->GetVertices())[1+2*i]);
      A3->SetVertex(i+4,(A3->GetVertices())[2*i],(A3->GetVertices())[1+2*i]);
      //
      B1->SetVertex(i+4,(B1->GetVertices())[2*i],(B1->GetVertices())[1+2*i]);
      B2->SetVertex(i+4,(B2->GetVertices())[2*i],(B2->GetVertices())[1+2*i]);
      B3->SetVertex(i+4,(B3->GetVertices())[2*i],(B3->GetVertices())[1+2*i]);
      // C's are a cone which must match up with B's.
      C1->SetVertex(i+4,(B1->GetVertices())[2*i],(B1->GetVertices())[1+2*i]);
      C2->SetVertex(i+4,(B2->GetVertices())[2*i],(B2->GetVertices())[1+2*i]);
      C3->SetVertex(i+4,(B3->GetVertices())[2*i],(B3->GetVertices())[1+2*i]);
      //
      Ah1->SetVertex(i+4,(Ah1->GetVertices())[2*i],
                         (Ah1->GetVertices())[1+2*i]);
      Ah2->SetVertex(i+4,(Ah2->GetVertices())[2*i],
                         (Ah2->GetVertices())[1+2*i]);
      Ah3->SetVertex(i+4,(Ah3->GetVertices())[2*i],
                         (Ah3->GetVertices())[1+2*i]);
      //
      Bh1->SetVertex(i+4,(Bh1->GetVertices())[2*i],
                         (Bh1->GetVertices())[1+2*i]);
      Bh2->SetVertex(i+4,(Bh2->GetVertices())[2*i],
                         (Bh2->GetVertices())[1+2*i]);
      Bh3->SetVertex(i+4,(Bh3->GetVertices())[2*i],
                         (Bh3->GetVertices())[1+2*i]);
    } // end for
    //
    p1 = (xco[0]-xci[0])/(yco[0]-yci[0]);
    c1 = xco[0]+0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(xco[0]-xci[0])+
                                           SQ(yco[0]-yci[0]))/(xco[0]-xci[0]);
    y = TSCRoutC-2.*TSCarbonFiberThA;
    x = p1*(y-yco[0])+c1;
    Ch1->SetVertex(0,x,y);
    y = TSCRinC+TSCarbonFiberThA;
    x = p1*(y-yci[0])+c1;
    Ch1->SetVertex(2,x,y);
    p1 = (xco[1]-xci[1])/(yco[1]-yci[1]);
    c1 = xco[1]-0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(xco[1]-xci[1])+
                                           SQ(yco[1]-yci[1]))/(xco[1]-xci[1]);
    y = TSCRoutC-2.*TSCarbonFiberThA;
    x = p1*(y-yco[1])+c1;
    Ch1->SetVertex(1,x,y);
    y = TSCRinC+TSCarbonFiberThA;
    x = p1*(y-yci[1])+c1;
    Ch1->SetVertex(3,x,y);
    //
    th = 0.5*TSCAngle*kRadian;
    xa = TMath::Cos(th)*xco[1]-TMath::Sin(th)*yco[1];
    ya = TMath::Sin(th)*xco[1]+TMath::Cos(th)*yco[1];
    xb = TMath::Cos(th)*xci[1]-TMath::Sin(th)*yci[1];
    yb = TMath::Sin(th)*xci[1]+TMath::Cos(th)*yci[1];
    p1 = (xa-xb)/(ya-yb);
    c1 = xa+0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(xa-xb)+SQ(ya-yb))/(xa-xb);
    y = ya-TSCarbonFiberThA;
    x = p1*(y-ya)+c1;
    xp = TMath::Cos(-th)*x-TMath::Sin(-th)*y;
    yp = TMath::Sin(-th)*x+TMath::Cos(-th)*y;
    yp = ya-TSCarbonFiberThA;
    xp = p1*(y-ya)+c1;
    Ch2->SetVertex(0,xp,yp);
    y = yb+2.0*TSCarbonFiberThA;
    x = p1*(y-ya)+c1;
    xp = TMath::Cos(-th)*x-TMath::Sin(-th)*y;
    yp = TMath::Sin(-th)*x+TMath::Cos(-th)*y;
    Ch2->SetVertex(2,xp,yp);
    xa = TMath::Cos(th)*xco[2]-TMath::Sin(th)*yco[2];
    ya = TMath::Sin(th)*xco[2]+TMath::Cos(th)*yco[2];
    xb = TMath::Cos(th)*xci[2]-TMath::Sin(th)*yci[2];
    yb = TMath::Sin(th)*xci[2]+TMath::Cos(th)*yci[2];
    p1 = (xa-xb)/(ya-yb);
    c1 = xa-0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(xa-xb)+SQ(ya-yb))/(xa-xb);
    y = ya-TSCarbonFiberThA;
    x = p1*(y-ya)+c1;
    xp = TMath::Cos(-th)*x-TMath::Sin(-th)*y;
    yp = TMath::Sin(-th)*x+TMath::Cos(-th)*y;
    Ch2->SetVertex(1,xp,yp);
    y = yb+2.0*TSCarbonFiberThA;
    x = p1*(y-ya)+c1;
    xp = TMath::Cos(-th)*x-TMath::Sin(-th)*y;
    yp = TMath::Sin(-th)*x+TMath::Cos(-th)*y;
    Ch2->SetVertex(3,xp,yp);
    //
    p1 = (yco[5]-yci[5])/(xco[5]-xci[5]);
    c1 = yco[5]+0.5*TSCarbonFiberThA*TMath::Sqrt(SQ(yco[5]-yci[5])+
                                          SQ(xco[5]-xci[5]))/(yco[5]-yci[5]);
    x = xco[5]-TSCarbonFiberThA;
    y = p1*(x-xco[5])+c1;
    Ch3->SetVertex(0,x,y);
    x = xci[5]+2.0*TSCarbonFiberThA;
    y = p1*(x-xci[5])+c1;
    Ch3->SetVertex(2,x,y);
    y = 2.0*TSCarbonFiberThA;
    x = xco[5]-TSCarbonFiberThA;
    Ch3->SetVertex(1,x,y);
    y = 2.0*TSCarbonFiberThA;
    x = xci[5]+2.0*TSCarbonFiberThA;
    Ch3->SetVertex(3,x,y);
    Ch1->SetVertex(7,x7[0],y7[0]); // 7th point most be done last ???
    Ch2->SetVertex(7,x7[1],y7[1]); // 7th point most be done last ???
    Ch3->SetVertex(7,x7[2],y7[2]); // 7th point most be done last ???
    printArb8(A1);
    printArb8(Ah1);
    printArb8(A2);
    printArb8(Ah2);
    printArb8(A3);
    printArb8(Ah3);
    printArb8(B1);
    printArb8(Bh1);
    printArb8(B2);
    printArb8(Bh2);
    printArb8(B3);
    printArb8(Bh3);
    printArb8(C1);
    printArb8(Ch1);
    printArb8(C2);
    printArb8(Ch2);
    printArb8(C3);
    printArb8(Ch3);
    //
    // Define Minimal volume to inclose this SPD Thermal Sheald.
    TGeoPcon *M1 = new TGeoPcon("ITSspdShealdVV",0.0,360.0,9);
    M1->Z(0)    = 0.5*TSCLengthA+TSCLengthB;
    M1->Rmin(0) = TSCRinB;
    x = B1->GetVertices()[0]; // [0][0]
    y = B1->GetVertices()[1]; // [0][1]
    M1->Rmax(0) = TMath::Sqrt(x*x+y*y);
    M1->Z(1)    = M1->GetZ(0)-TSCLengthB;
    M1->Rmin(1) = M1->GetRmin(0);
    M1->Rmax(1) = M1->GetRmax(0);
    M1->Z(2)    = M1->GetZ(1);
    M1->Rmin(2) = TSCRinA;
    x = A1->GetVertices()[0]; // [0]0]
    y = A1->GetVertices()[1]; // [0][1]
    M1->Rmax(2) = TMath::Sqrt(x*x+y*y);
    M1->Z(3)    = -(M1->GetZ(0)-TSCLengthB);
    M1->Rmin(3) = M1->GetRmin(2);
    M1->Rmax(3) = M1->GetRmax(2);
    M1->Z(4)    = M1->GetZ(3);
    M1->Rmin(4) = M1->GetRmin(1);
    M1->Rmax(4) = M1->GetRmax(1);
    M1->Z(5)    = -(M1->GetZ(0));
    M1->Rmin(5) = M1->GetRmin(0);
    M1->Rmax(5) = M1->GetRmax(0);
    M1->Z(6)    = M1->GetZ(5) - TSCLengthC;
    M1->Rmin(6) = TSCRinC;
    x = C1->GetVertices()[0]; // [0][0]
    y = C1->GetVertices()[1]; // [0][1]
    M1->Rmax(6) = TMath::Sqrt(x*x+y*y);
    M1->Z(7)    = M1->GetZ(6);
    M1->Rmin(7) = D->GetRmin();
    M1->Rmax(7) = D->GetRmax();
    M1->Z(8)    = M1->Z(7) - TSCLengthD;
    M1->Rmin(8) = M1->GetRmin(7);
    M1->Rmax(8) = M1->GetRmax(7);
    TGeoTubeSeg *M2 = new TGeoTubeSeg("ITSspdShealdWingVV",
       M1->GetRmax(8),Dw->GetRmax(),Dw->GetDz(),Dw->GetPhi1(),Dw->GetPhi2());
    printTubeSeg(M2);
    //
    x = 0.5*(M1->GetZ(8) + M1->GetZ(7));
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
    M = new TGeoCompositeShape("ITS SPD Thermal sheald volume",
                              "(((ITSspdShealdVV+"
                              "ITSspdShealdWingVV:ITSspdShealdVVt0)+"
                              "ITSspdShealdWingVV:ITSspdShealdVVt1)+"
                              "ITSspdShealdWingVV:ITSspdShealdVVt2)+"
                              "ITSspdShealdWingVV:ITSspdShealdVVt3");
    //
    TGeoManager *mgr = gGeoManager;
    SPDcf = mgr->GetMedium("ITSspdCarbonFiber");
    SPDfs = mgr->GetMedium("ITSspdStaselite4411w");
    SPDfo = mgr->GetMedium("ITSspdRohacell50A");
    SPDss = mgr->GetMedium("ITSspdStainlessSteal");
    SPDair= mgr->GetMedium("ITSspdAir");
    TGeoVolume *A1v,*A2v,*A3v,*Ah1v,*Ah2v,*Ah3v;
    TGeoVolume *B1v,*B2v,*B3v,*Bh1v,*Bh2v,*Bh3v;
    TGeoVolume *C1v,*C2v,*C3v,*Ch1v,*Ch2v,*Ch3v;
    TGeoVolume *Dv,*Dsv,*Dwv,*Dwsv,*Mv;
    Mv = new TGeoVolume("ITSspdThermalSheald",M,SPDair);
    Mv->SetVisibility(kTRUE);
    Mv->SetLineColor(7); // light Blue
    Mv->SetLineWidth(1);
    Mv->SetFillColor(Mv->GetLineColor());
    Mv->SetFillStyle(4090); // 90% transparent
    Moth->AddNode(Mv,1,0); ///////////////////// Virtual Volume ////////
    A1v = new TGeoVolume("ITSspdCentCylA1CF",A1,SPDcf);
    A1v->SetVisibility(kTRUE);
    A1v->SetLineColor(4);
    A1v->SetLineWidth(1);
    A2v = new TGeoVolume("ITSspdCentCylA2CF",A2,SPDcf);
    A2v->SetVisibility(kTRUE);
    A2v->SetLineColor(4);
    A2v->SetLineWidth(1);
    A3v = new TGeoVolume("ITSspdCentCylA3CF",A3,SPDcf);
    A3v->SetVisibility(kTRUE);
    A3v->SetLineColor(4);
    A3v->SetLineWidth(1);
    B1v = new TGeoVolume("ITSspdCentCylB1CF",B1,SPDcf);
    B1v->SetVisibility(kTRUE);
    B1v->SetLineColor(4);
    B1v->SetLineWidth(1);
    B2v = new TGeoVolume("ITSspdCentCylB2CF",B2,SPDcf);
    B2v->SetVisibility(kTRUE);
    B2v->SetLineColor(4);
    B2v->SetLineWidth(1);
    B3v = new TGeoVolume("ITSspdCentCylB3CF",B3,SPDcf);
    B3v->SetVisibility(kTRUE);
    B3v->SetLineColor(4);
    B3v->SetLineWidth(1);
    C1v = new TGeoVolume("ITSspdCentCylC1CF",C1,SPDcf);
    C1v->SetVisibility(kTRUE);
    C1v->SetLineColor(4);
    C1v->SetLineWidth(1);
    C2v = new TGeoVolume("ITSspdCentCylC2CF",C2,SPDcf);
    C2v->SetVisibility(kTRUE);
    C2v->SetLineColor(4);
    C2v->SetLineWidth(1);
    C3v = new TGeoVolume("ITSspdCentCylC3CF",C3,SPDcf);
    C3v->SetVisibility(kTRUE);
    C3v->SetLineColor(4);
    C3v->SetLineWidth(1);
    Ah1v = new TGeoVolume("ITSspdCentCylA1AirA",Ah1,SPDair);
    Ah1v->SetVisibility(kTRUE);
    Ah1v->SetLineColor(5); // Yellow
    Ah1v->SetFillColor(Ah1v->GetLineColor());
    Ah1v->SetFillStyle(4090); // 90% transparent
    Ah2v = new TGeoVolume("ITSspdCentCylA2AirA",Ah2,SPDair);
    Ah2v->SetVisibility(kTRUE);
    Ah2v->SetLineColor(5); // Yellow
    Ah2v->SetFillColor(Ah2v->GetLineColor());
    Ah2v->SetFillStyle(4090); // 90% transparent
    Ah3v = new TGeoVolume("ITSspdCentCylA3AirA",Ah3,SPDair);
    Ah3v->SetVisibility(kTRUE);
    Ah3v->SetLineColor(5); // Yellow
    Ah3v->SetFillColor(Ah3v->GetLineColor());
    Ah3v->SetFillStyle(4090); // 90% transparent
    Bh1v = new TGeoVolume("ITSspdCentCylA1AirB",Bh1,SPDair);
    Bh1v->SetVisibility(kTRUE);
    Bh1v->SetLineColor(5); // Yellow
    Bh1v->SetFillColor(Bh1v->GetLineColor());
    Bh1v->SetFillStyle(4090); // 90% transparent
    Bh2v = new TGeoVolume("ITSspdCentCylA2AirB",Bh2,SPDair);
    Bh2v->SetVisibility(kTRUE);
    Bh2v->SetLineColor(5); // Yellow
    Bh2v->SetFillColor(Bh2v->GetLineColor());
    Bh2v->SetFillStyle(4090); // 90% transparent
    Bh3v = new TGeoVolume("ITSspdCentCylA3AirB",Bh3,SPDair);
    Bh3v->SetVisibility(kTRUE);
    Bh3v->SetLineColor(5); // Yellow
    Bh3v->SetFillColor(Bh3v->GetLineColor());
    Bh3v->SetFillStyle(4090); // 90% transparent
    Ch1v = new TGeoVolume("ITSspdCentCylA1AirC",Ch1,SPDair);
    Ch1v->SetVisibility(kTRUE);
    Ch1v->SetLineColor(5); // Yellow
    Ch1v->SetFillColor(Ch1v->GetLineColor());
    Ch1v->SetFillStyle(4090); // 90% transparent
    Ch2v = new TGeoVolume("ITSspdCentCylA2AirC",Ch2,SPDair);
    Ch2v->SetVisibility(kTRUE);
    Ch2v->SetLineColor(5); // Yellow
    Ch2v->SetFillColor(Ch2v->GetLineColor());
    Ch2v->SetFillStyle(4090); // 90% transparent
    Ch3v = new TGeoVolume("ITSspdCentCylA3AirC",Ch3,SPDair);
    Ch3v->SetVisibility(kTRUE);
    Ch3v->SetLineColor(5); // Yellow
    Ch3v->SetFillColor(Ch3v->GetLineColor());
    Ch3v->SetFillStyle(4090); // 90% transparent
    Dv = new TGeoVolume("ITSspdCentCylA1CD",D,SPDcf);
    Dv->SetVisibility(kTRUE);
    Dv->SetLineColor(4);
    Dv->SetLineWidth(1);
    Dwv = new TGeoVolume("ITSspdCentCylA1CDw",Dw,SPDcf);
    Dwv->SetVisibility(kTRUE);
    Dwv->SetLineColor(4);
    Dwv->SetLineWidth(1);
    Dsv = new TGeoVolume("ITSspdCentCylA1Dfill",Ds,SPDfs);
    Dsv->SetVisibility(kTRUE);
    Dsv->SetLineColor(3); // Green
    Dsv->SetFillColor(Dsv->GetLineColor());
    Dsv->SetFillStyle(4010); // 10% transparent
    Dwsv = new TGeoVolume("ITSspdCentCylA1DwingFill",Dws,SPDfs);
    Dwsv->SetVisibility(kTRUE);
    Dwsv->SetLineColor(3); // Green
    Dwsv->SetFillColor(Dwsv->GetLineColor());
    Dwsv->SetFillStyle(4010); // 10% transparent
    //
    A1v->AddNode(Ah1v,1,0);
    A2v->AddNode(Ah2v,1,0);
    A3v->AddNode(Ah3v,1,0);
    B1v->AddNode(Bh1v,1,0);
    B2v->AddNode(Bh2v,1,0);
    B3v->AddNode(Bh3v,1,0);
    C1v->AddNode(Ch1v,1,0);
    C2v->AddNode(Ch2v,1,0);
    C3v->AddNode(Ch3v,1,0);
    Dv ->AddNode(Dsv ,1,0);
    Dwv->AddNode(Dwsv,1,0);
    //
    Mv->AddNode(A1v,1,0);
    Mv->AddNode(A2v,1,0);
    Mv->AddNode(A3v,1,0);
    tranb  = new TGeoTranslation("",0.0,0.0,0.5*(TSCLengthA+TSCLengthB));
    tranbm = new TGeoTranslation("",0.0,0.0,0.5*(-TSCLengthA-TSCLengthB));
    Mv->AddNode(B1v,1,tranb);
    Mv->AddNode(B2v,1,tranb);
    Mv->AddNode(B3v,1,tranb);
    Mv->AddNode(B1v,2,tranbm);
    Mv->AddNode(B2v,2,tranbm);
    Mv->AddNode(B3v,2,tranbm);
    // Muon side (rb26) is at -Z.
    tranc = new TGeoTranslation("",0.0,0.0,
                                0.5*(-TSCLengthA-TSCLengthB-TSCLengthC));
    Mv->AddNode(C1v,1,tranc);
    Mv->AddNode(C2v,1,tranc);
    Mv->AddNode(C3v,1,tranc);
    Mv->AddNode(Dv,1,tranITSspdShealdVVt0);
    Mv->AddNode(Dwv,1,tranITSspdShealdVVt0);
    Mv->AddNode(Dwv,2,rotITSspdShealdVVt1);
    Mv->AddNode(Dwv,3,rotITSspdShealdVVt2);
    Mv->AddNode(Dwv,4,rotITSspdShealdVVt3);
    k=2;
    for(i=1;i<10;i++) {
        th = ((Double_t)i)*TSCAngle*kDegree;
        rot = new TGeoRotation("",0.0,0.0,th);
        Mv->AddNode(A1v,i+1,rot);
        Mv->AddNode(B1v,i+2,new TGeoCombiTrans(*tranb,*rot));
        Mv->AddNode(B1v,i+12,new TGeoCombiTrans(*tranbm,*rot));
        Mv->AddNode(C1v,i+1,new TGeoCombiTrans(*tranc,*rot));
        if(i!=0||i!=2||i!=7){
            Mv->AddNode(A2v,k++,rot);
            Mv->AddNode(B2v,k++,new TGeoCombiTrans(*tranb,*rot));
            Mv->AddNode(B2v,k++,new TGeoCombiTrans(*tranbm,*rot));
            Mv->AddNode(C2v,k++,new TGeoCombiTrans(*tranc,*rot));
        } // end if
        if(i==5) {
            Mv->AddNode(A3v,2,rot);
            Mv->AddNode(B3v,3,new TGeoCombiTrans(*tranb,*rot));
            Mv->AddNode(B3v,4,new TGeoCombiTrans(*tranbm,*rot));
            Mv->AddNode(C3v,2,new TGeoCombiTrans(*tranc,*rot));
        } // end if
    } // end for i
    rot = new TGeoRotation("",180.,0.0,0.0);
    Mv->AddNode(A3v,3,rot);
    Mv->AddNode(B3v,5,new TGeoCombiTrans(*tranb,*rot));
    Mv->AddNode(B3v,6,new TGeoCombiTrans(*tranbm,*rot));
    Mv->AddNode(C3v,3,new TGeoCombiTrans(*tranc,*rot));
    rot = new TGeoRotation("",180.,0.0,180.0);
    Mv->AddNode(A3v,4,rot);
    Mv->AddNode(B3v,7,new TGeoCombiTrans(*tranb,*rot));
    Mv->AddNode(B3v,8,new TGeoCombiTrans(*tranbm,*rot));
    Mv->AddNode(C3v,4,new TGeoCombiTrans(*tranc,*rot));
    if(GetDebug()){
        A1v->PrintNodes();
        Ah1v->PrintNodes();
        A2v->PrintNodes();
        Ah2v->PrintNodes();
        A3v->PrintNodes();
        Ah3v->PrintNodes();
        B1v->PrintNodes();
        Bh1v->PrintNodes();
        B2v->PrintNodes();
        Bh2v->PrintNodes();
        B3v->PrintNodes();
        Bh3v->PrintNodes();
        C1v->PrintNodes();
        Ch1v->PrintNodes();
        C2v->PrintNodes();
        Ch2v->PrintNodes();
        C3v->PrintNodes();
        Ch3v->PrintNodes();
        Dv->PrintNodes();
        Dsv->PrintNodes();
        Dwv->PrintNodes();
        Dwsv->PrintNodes();
        //Mv->PrintNodes();
    } // end if
}
//______________________________________________________________________
void AliITSv11GeometrySupport::SDDCone(TGeoVolume *Moth){
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
    const Double_t TSLength       = 790.0*kmm; // Thermal Sheeld length
    const Double_t TSInsertoLength= 15.0*kmm;    // ????
    const Double_t TSOuterR       = 0.5*(220.+10.)*kmm; // ????
    const Double_t TSInnerR       = 0.5*(220.-10.)*kmm; // ????
    const Double_t TSCarbonFiberth= 0.02*kmm;     // ????
    const Double_t TSBoltDiameter = 6.0*kmm; // M6 screw
    const Double_t TSBoltDepth    = 6.0*kmm; // in volume C
    const Double_t TSBoltRadius   = 0.5*220.*kmm; // Radius in volume C
    const Double_t TSBoltAngle0   = 0.0*kDegree; // Angle in volume C
    const Double_t TSBoltdAngle   = 30.0*kDegree; // Angle in Volume C
    Double_t x,y,z,t,t0;
    Int_t i,n;
    TGeoTube *A,*B,*C,*D;
    TGeoTranslation *tran;
    TGeoRotation *rot;
    TGeoCombiTrans *rotran;
    TGeoMedium *SDDcf,*SDDfs,*SDDfo,*SDDss;

    A = new TGeoTube("ITS SDD Central Cylinder",TSInnerR,TSOuterR,
                     0.5*TSLength);
    B = new TGeoTube("ITS SDD CC Foam",TSInnerR+TSCarbonFiberth,
                    TSOuterR-TSCarbonFiberth,
                    0.5*(TSLength-2.0*TSInsertoLength));
    C = new TGeoTube("ITS SDD CC Inserto",TSInnerR+TSCarbonFiberth,
                    TSOuterR-TSCarbonFiberth,0.5*TSLength);
    D = new TGeoTube("ITS SDD CC M6 bolt end",0.0,0.5*TSBoltDiameter,
                    0.5*TSBoltDepth);
    printTube(A);
    printTube(B);
    printTube(C);
    printTube(D);
    //
    TGeoManager *mgr = gGeoManager;
    SDDcf = mgr->GetMedium("ITSssdCarbonFiber");
    SDDfs = mgr->GetMedium("ITSssdStaselite4411w");
    SDDfo = mgr->GetMedium("ITSssdRohacell50A");
    SDDss = mgr->GetMedium("ITSssdStainlessSteal");
    TGeoVolume *Av,*Bv,*Cv,*Dv;
    Av = new TGeoVolume("ITSsddCentCylCF",A,SDDcf);
    Av->SetVisibility(kTRUE);
    Av->SetLineColor(4);
    Av->SetLineWidth(1);
    Av->SetFillColor(Av->GetLineColor());
    Av->SetFillStyle(4000); // 0% transparent
    Bv = new TGeoVolume("ITSsddCentCylF",B,SDDfo);
    Bv->SetVisibility(kTRUE);
    Bv->SetLineColor(3);
    Bv->SetLineWidth(1);
    Bv->SetFillColor(Bv->GetLineColor());
    Bv->SetFillStyle(4000); // 0% transparent
    Cv = new TGeoVolume("ITSsddCentCylSt",C,SDDfs);
    Cv->SetVisibility(kTRUE);
    Cv->SetLineColor(2);
    Cv->SetLineWidth(1);
    Cv->SetFillColor(Cv->GetLineColor());
    Cv->SetFillStyle(4000); // 0% transparent
    Dv = new TGeoVolume("ITSsddCentCylSS",D,SDDss);
    Dv->SetVisibility(kTRUE);
    Dv->SetLineColor(1);
    Dv->SetLineWidth(1);
    Dv->SetFillColor(Dv->GetLineColor());
    Dv->SetFillStyle(4000); // 0% transparent
    //
    Moth->AddNode(Av,1,0);
    Av->AddNode(Cv,1,0);
    Cv->AddNode(Bv,1,0);
    n = (Int_t)((360.*kDegree)/TSBoltdAngle);
    for(i=0;i<n;i++){
        t = TSBoltAngle0+((Double_t)i)*TSBoltdAngle;
        x = TSBoltRadius*TMath::Cos(t*kRadian);
        y = TSBoltRadius*TMath::Sin(t*kRadian);
        z = 0.5*(TSLength-TSBoltDepth);
        tran = new TGeoTranslation("",x,y,z);
        Cv->AddNode(Dv,i+1,tran);
        tran = new TGeoTranslation("",x,y,-z);
        Cv->AddNode(Dv,i+n+1,tran);
    } // end for i
    if(GetDebug()){
        Av->PrintNodes();
        Bv->PrintNodes();
        Cv->PrintNodes();
        Dv->PrintNodes();
    } // end if
    // SDD Suport Cone
    //
    //
    const Double_t Thickness = 10.5*kmm; // Thickness of Rohacell+carbon fiber
    const Double_t Cthick    = 1.5*kmm; // Carbon finber thickness
    const Double_t Rcurv     = 15.0*kmm; // Radius of curvature.
    const Double_t Tc        = 45.0; // angle of SDD cone [degrees].
    const Double_t Sintc = TMath::Sin(Tc*TMath::DegToRad());
    const Double_t Costc = TMath::Cos(Tc*TMath::DegToRad());
    const Double_t Tantc = TMath::Tan(Tc*TMath::DegToRad());
    const Double_t ZouterMilled = 23.0*kmm;
    const Double_t Zcylinder    = 186.0*kmm;
    const Double_t Z0           = Zcylinder + 0.5*TSLength;
    //const Int_t Nspoaks         = 12;
    //const Int_t Nmounts         = 4;
    //const Double_t DmountAngle  = 9.0; // degrees
    const Double_t RoutMax      = 0.5*560.0*kmm;
    const Double_t RoutMin      = 0.5*539.0*kmm;
    // Holes in cone for cables
    const Double_t PhiHole1     = 0.0*kDegree;
    const Double_t dPhiHole1    = 25.0*kDegree;
    const Double_t RholeMax1    = 0.5*528.*kmm;
    const Double_t RholeMin1    = 0.5*464.*kmm;
    const Double_t PhiHole2     = 0.0*kDegree;
    const Double_t dPhiHole2    = 50.0*kDegree;
    const Double_t RholeMax2    = 0.5*375.*kmm;
    const Double_t RholeMin2    = 0.5*280.*kmm;
    //
    //const Int_t NpostsOut       = 6;
    //const Int_t NpostsIn        = 3;
    //const Double_t Phi0PostOut  = 0.0; // degree
    //const Double_t Phi0PostIn   = 0.0; // degree
    //const Double_t dRpostOut    = 16.0*kmm;
    //const Double_t dRpostIn     = 16.0*kmm;
    //const Double_t ZpostMaxOut  = 116.0*kmm;
    //const Double_t ZpostMaxIn   = 190.0*kmm;
    const Double_t RinMax       = 0.5*216*kmm;
    const Double_t RinCylinder  = 0.5*231.0*kmm;
    //const Double_t RinHole      = 0.5*220.0*kmm;
    const Double_t RinMin       = 0.5*210.0*kmm;
    const Double_t dZin         = 15.0*kmm; // ???
    //
    Double_t dza = Thickness/Sintc-(RoutMax-RoutMin)/Tantc;
    Double_t Z,Rmin,Rmax; // Temp variables.
    if(dza<=0){ // The number or order of the points are in error for a proper
     // call to pcons!
     Error("SDDcone","The definition of the points for a call to PCONS is"
           " in error. abort.");
     return;
    } // end if
    TGeoPcon *E = new TGeoPcon("ITSsddSuportConeCarbonFiberSurfaceE",
                               0.0,360.0,12);
    E->Z(0)    = 0.0;
    E->Rmin(0) = RoutMin;
    E->Rmax(0) = RoutMax;
    E->Z(1)    = ZouterMilled - dza;
    E->Rmin(1) = E->GetRmin(0);
    E->Rmax(1) = E->GetRmax(0);
    E->Z(2)    = ZouterMilled;
    E->Rmax(2) = E->GetRmax(0);
    RadiusOfCurvature(Rcurv,0.,E->GetZ(1),E->GetRmin(1),Tc,Z,Rmin);
    E->Z(3)    = Z;
    E->Rmin(3) = Rmin;
    E->Rmin(2) = RminFrom2Points(E,3,1,E->GetZ(2));
    RadiusOfCurvature(Rcurv,0.,E->GetZ(2),E->GetRmax(2),Tc,Z,Rmax);
    E->Z(4)    = Z;
    E->Rmax(4) = Rmax;
    E->Rmin(4) = RminFromZpCone(E,Tc,E->GetZ(4),0.0);
    E->Rmax(3) = RmaxFrom2Points(E,4,2,E->GetZ(3));
    E->Rmin(7) = RinMin;
    E->Rmin(8) = RinMin;
    RadiusOfCurvature(Rcurv,90.0,0.0,RinMax,90.0-Tc,Z,Rmax);
    E->Rmax(8) = Rmax;
    E->Z(8)    = ZFromRmaxpCone(E,Tc,E->GetRmax(8));
    E->Z(9)    = Zcylinder;
    E->Rmin(9) = RinMin;
    E->Z(10)    = E->GetZ(9);
    E->Rmin(10) = RinCylinder;
    E->Rmin(11) = RinCylinder;
    E->Rmax(11) = E->GetRmin(11);
    Rmin        = E->GetRmin(8);
    RadiusOfCurvature(Rcurv,90.0-Tc,E->GetZ(8),E->GetRmax(8),90.0,Z,Rmax);
    Rmax = RinMax;
    E->Z(11)    = Z+(E->GetZ(8)-Z)*(E->GetRmax(11)-Rmax)/(E->GetRmax(8)-Rmax);
    E->Rmax(9) = RmaxFrom2Points(E,11,8,E->GetZ(9));
    E->Rmax(10) = E->GetRmax(9);
    E->Z(6)    = Z-dZin;
    E->Z(7)    = E->GetZ(6);
    E->Rmax(6) = RmaxFromZpCone(E,Tc,E->GetZ(6));
    E->Rmax(7) = E->GetRmax(6);
    RadiusOfCurvature(Rcurv,90.,E->GetZ(6),0.0,90.0-Tc,Z,Rmin);
    E->Z(5)    = Z;
    E->Rmin(5) = RminFromZpCone(E,Tc,Z);
    E->Rmax(5) = RmaxFromZpCone(E,Tc,Z);
    RadiusOfCurvature(Rcurv,90.-Tc,0.0,E->Rmin(5),90.0,Z,Rmin);
    E->Rmin(6) = Rmin;
    printPcon(E);
    // Inner Core, Inserto material
    TGeoPcon *F = new TGeoPcon("ITSsddSuportConeInsertoStesaliteF",
                               0.,360.0,9);
    F->Z(0)    = E->GetZ(0);
    F->Rmin(0) = E->GetRmin(0)+Cthick;
    F->Rmax(0) = E->GetRmax(0)-Cthick;
    F->Z(1)    = E->GetZ(1);
    F->Rmin(1) = F->GetRmin(0);
    F->Rmax(1) = F->GetRmax(0);
    F->Z(2)    = E->GetZ(2);
    F->Rmax(2) = F->GetRmax(1);
    RadiusOfCurvature(Rcurv-Cthick,0.,F->GetZ(1),F->GetRmax(1),Tc,Z,Rmin);
    F->Z(3)    = Z;
    F->Rmin(3) = Rmin;
    F->Rmin(2) = RminFrom2Points(F,3,1,F->GetZ(2));
    RadiusOfCurvature(Rcurv+Cthick,0.,F->GetZ(2),F->GetRmax(2),Tc,Z,Rmax);
    F->Z(4)    = Z;
    F->Rmax(4) = Rmax;
    F->Rmin(4) = RmaxFromZpCone(E,Tc,F->GetZ(4),-Cthick);
    F->Rmax(3) = RmaxFrom2Points(F,4,2,F->GetZ(3));
    F->Rmin(7) = E->GetRmin(7);
    F->Rmin(8) = E->GetRmin(8);
    F->Z(6)    = E->GetZ(6)+Cthick;
    F->Rmin(6) = E->GetRmin(6);
    F->Z(7)    = F->GetZ(6);
    F->Rmax(8) = E->GetRmax(8)-Cthick*Sintc;
    RadiusOfCurvature(Rcurv+Cthick,90.0,F->GetZ(6),F->GetRmin(6),90.0-Tc,
                      Z,Rmin);
    F->Z(5)    = Z;
    F->Rmin(5) = Rmin;
    F->Rmax(5) = RmaxFromZpCone(F,Tc,Z);
    F->Rmax(6) = RmaxFromZpCone(F,Tc,F->GetZ(6));
    F->Rmax(7) = F->GetRmax(6);
    F->Z(8)    = ZFromRmaxpCone(F,Tc,F->GetRmax(8),-Cthick);
    printPcon(F);
    // Inner Core, Inserto material
    TGeoPcon *G = new TGeoPcon("ITSsddSuportConeFoamCoreG",0.0,360.0,4);
    RadiusOfCurvature(Rcurv+Cthick,0.0,F->GetZ(1),F->GetRmin(1),Tc,Z,Rmin);
    G->Z(0)    = Z;
    G->Rmin(0) = Rmin;
    G->Rmax(0) = G->GetRmin(0);
    G->Z(1)    = G->GetZ(0)+(Thickness-2.0*Cthick)/Sintc;;
    G->Rmin(1) = RminFromZpCone(F,Tc,G->GetZ(1));
    G->Rmax(1) = RmaxFromZpCone(F,Tc,G->GetZ(1));
    G->Z(2)    = E->GetZ(5)-Cthick;
    G->Rmin(2) = RminFromZpCone(F,Tc,G->GetZ(2));
    G->Rmax(2) = RmaxFromZpCone(F,Tc,G->GetZ(2));
    G->Z(3)    = F->GetZ(5)+(Thickness-2.0*Cthick)*Costc;
    G->Rmax(3) = RmaxFromZpCone(F,Tc,G->GetZ(3));
    G->Rmin(3) = G->GetRmax(3);
    printPcon(G);
    //
    TGeoPcon *H = new TGeoPcon("ITSsddSuportConeHoleH",PhiHole1,dPhiHole1,4);
    H->Rmin(0) = RholeMax1;
    H->Rmax(0) = H->GetRmin(0);
    H->Z(0)    = ZFromRminpCone(E,Tc,H->GetRmin(0));
    H->Rmax(1) = H->GetRmax(0);
    H->Z(1)    = ZFromRmaxpCone(E,Tc,H->GetRmax(1));
    H->Rmin(1) = RminFromZpCone(E,Tc,H->GetZ(1));
    H->Rmin(2) = RholeMin1;
    H->Z(2)    = ZFromRminpCone(E,Tc,H->GetRmin(2));
    H->Rmax(2) = RmaxFromZpCone(E,Tc,H->GetZ(2));
    H->Rmin(3) = H->GetRmin(2);
    H->Rmax(3) = H->GetRmin(3);
    H->Z(3)    = ZFromRminpCone(E,Tc,H->GetRmin(3));
    printPcon(H);
    //
    x = Cthick/(0.5*(RholeMax1+RholeMin1));
    t0 = PhiHole1 - x/kRadian;
    t  = dPhiHole1 + 2.0*x/kRadian;
    TGeoPcon *I = new TGeoPcon("ITSsddSuportConeHoleI",t0,t,4);
    I->Rmin(0) = RholeMax1+Cthick;
    I->Rmax(0) = I->GetRmin(0);
    I->Z(0)    = ZFromRminpCone(F,Tc,I->GetRmin(0));
    I->Rmax(1) = I->GetRmax(0);
    I->Z(1)    = ZFromRmaxpCone(F,Tc,I->GetRmax(1));
    I->Rmin(1) = RminFromZpCone(F,Tc,I->GetZ(1));
    I->Rmin(2) = RholeMin1-Cthick;
    I->Z(2)    = ZFromRminpCone(F,Tc,I->GetRmin(2));
    I->Rmax(2) = RmaxFromZpCone(F,Tc,I->GetZ(2));
    I->Rmin(3) = I->GetRmin(2);
    I->Rmax(3) = I->GetRmin(3);
    I->Z(3)    = ZFromRmaxpCone(F,Tc,I->GetRmax(3));
    printPcon(I);
    //
    TGeoPcon *J = new TGeoPcon("ITSsddSuportConeHoleJ",PhiHole2,dPhiHole2,4);
    J->Rmin(0) = RholeMax2;
    J->Rmax(0) = J->GetRmin(0);
    J->Z(0)    = ZFromRminpCone(E,Tc,J->GetRmin(0));
    J->Rmax(1) = J->GetRmax(0);
    J->Z(1)    = ZFromRmaxpCone(E,Tc,J->GetRmax(1));
    J->Rmin(1) = RminFromZpCone(E,Tc,J->GetZ(1));
    J->Rmin(2) = RholeMin2;
    J->Z(2)    = ZFromRminpCone(E,Tc,J->GetRmin(2));
    J->Rmax(2) = RmaxFromZpCone(E,Tc,J->GetZ(2));
    J->Rmin(3) = J->GetRmin(2);
    J->Rmax(3) = J->GetRmin(3);
    J->Z(3)    = ZFromRmaxpCone(E,Tc,J->GetRmax(3));
    printPcon(J);
    //
    x = Cthick/(0.5*(RholeMax2+RholeMin2));
    t0 = PhiHole2 - x/kRadian;
    t  = dPhiHole2 + 2.0*x/kRadian;
    TGeoPcon *K = new TGeoPcon("ITSsddSuportConeHoleK",t0,t,4);
    K->Rmin(0) = RholeMax2+Cthick;
    K->Rmax(0) = K->GetRmin(0);
    K->Z(0)    = ZFromRminpCone(F,Tc,K->GetRmin(0));
    K->Rmax(1) = K->GetRmax(0);
    K->Z(1)    = ZFromRmaxpCone(F,Tc,K->GetRmax(1));
    K->Rmin(1) = RminFromZpCone(F,Tc,K->GetZ(1));
    K->Rmin(2) = RholeMin2-Cthick;
    K->Z(2)    = ZFromRminpCone(F,Tc,K->GetRmin(2));
    K->Rmax(2) = RmaxFromZpCone(F,Tc,K->GetZ(2));
    K->Rmin(3) = K->GetRmin(2);
    K->Rmax(3) = K->GetRmin(3);
    K->Z(3)    = ZFromRmaxpCone(F,Tc,K->GetRmax(3));
    printPcon(K);
    //
    TGeoCompositeShape *L,*M,*N;
    rot = new TGeoRotation("ITSsddRotZ30",0.0,0.0,30.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ60",0.0,0.0,60.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ90",0.0,0.0,90.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ120",0.0,0.0,120.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ150",0.0,0.0,150.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ180",0.0,0.0,180.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ210",0.0,0.0,210.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ240",0.0,0.0,240.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ270",0.0,0.0,270.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ300",0.0,0.0,300.0);
    rot->RegisterYourself();
    rot = new TGeoRotation("ITSsddRotZ330",0.0,0.0,330.0);
    rot->RegisterYourself();
    L = new TGeoCompositeShape("ITS SDD Suport Cone","((((((((((((((((("
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
    M = new TGeoCompositeShape("ITS SDD Suport Cone Inserto Stesalite",
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
    N = new TGeoCompositeShape("ITS SDD Suport Cone Foam Core",
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
    TGeoVolume *Lv,*Mv,*Nv;
    Lv = new TGeoVolume("ITSsddConeL",L,SDDcf);
    Lv->SetVisibility(kTRUE);
    Lv->SetLineColor(4);
    Lv->SetLineWidth(1);
    Lv->SetFillColor(Lv->GetLineColor());
    Lv->SetFillStyle(4000); // 0% transparent
    Mv = new TGeoVolume("ITSsddConeM",M,SDDfs);
    Mv->SetVisibility(kTRUE);
    Mv->SetLineColor(2);
    Mv->SetLineWidth(1);
    Mv->SetFillColor(Mv->GetLineColor());
    Mv->SetFillStyle(4010); // 10% transparent
    Nv = new TGeoVolume("ITSsddConeN",N,SDDfo);
    Nv->SetVisibility(kTRUE);
    Nv->SetLineColor(7);
    Nv->SetLineWidth(1);
    Nv->SetFillColor(Nv->GetLineColor());
    Nv->SetFillStyle(4050); // 50% transparent
    //
    Mv->AddNode(Nv,1,0);
    Lv->AddNode(Mv,1,0);
    tran = new TGeoTranslation("",0.0,0.0,-Z0);
    Moth->AddNode(Lv,1,tran);
    rot = new TGeoRotation("",0.0,180.0*kDegree,0.0);
    rotran = new TGeoCombiTrans("",0.0,0.0,Z0,rot);
    delete rot;// rot not explicity used in AddNode functions.
    Moth->AddNode(Lv,2,rotran);
    if(GetDebug()){
        Lv->PrintNodes();
        Mv->PrintNodes();
        Nv->PrintNodes();
    } // end if
}
//______________________________________________________________________
void AliITSv11GeometrySupport::SSDCone(TGeoVolume *Moth){
    // Define the detail SSD support cone geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
    //
    Int_t i,j;
    Double_t t,t0,dt,x,y,z,vl[3],vg[3],x0,y0;
    TGeoMedium *SSDcf  = 0; // SSD support cone Carbon Fiber materal number.
    TGeoMedium *SSDfs  = 0; // SSD support cone inserto stesalite 4411w.
    TGeoMedium *SSDfo  = 0; // SSD support cone foam, Rohacell 50A.
    TGeoMedium *SSDss  = 0; // SSD support cone screw material,Stainless steal
    TGeoMedium *SSDair = 0; // SSD support cone Air
    TGeoMedium *SSDal  = 0; // SSD support cone SDD mounting bracket Al
    TGeoManager *mgr = gGeoManager;
    SSDcf = mgr->GetMedium("ITSssdCarbonFiber");
    SSDfs = mgr->GetMedium("ITSssdStaselite4411w");
    SSDfo = mgr->GetMedium("ITSssdRohacell50A");
    SSDss = mgr->GetMedium("ITSssdStainlessSteal");
    SSDair= mgr->GetMedium("ITSssdAir");
    SSDal = mgr->GetMedium("ITSssdAl");
    //
    // SSD Central cylinder/Thermal Sheald.
    const Double_t CylZlength     = 1140.0*kmm; //
    const Double_t CylZFoamlength = 1020.0*kmm; //
    const Double_t CylROuter      = 0.5*595.0*kmm; //
    const Double_t CylRInner      = 0.5*560.5*kmm; //
    const Double_t CylCthick      = 0.64*kmm; //
    const Double_t CylFoamThick   = 5.0*kmm; //
    const Double_t CylRholes      = 0.5*575.0*kmm;
    const Double_t CylZM6         = 6.0*kmm; //
    const Double_t CylRM6         = 0.5*6.0*kmm;
    const Double_t CylPhi0M6      = 0.0*kDegree;
    const Int_t    CylNM6         = 40;
    const Double_t CylZPin        = 10.0*kmm;
    const Double_t CylRPin        = 0.5*4.0*kmm;
    const Double_t CylPhi0Pin     = (90.0+4.5)*kDegree;
    const Int_t    CylNPin        = 2;
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
    TGeoPcon *CA = new TGeoPcon("ITS SSD Thermal Centeral Carbon Fiber "
                                "CylinderCA",0.0,360.0,6);
    TGeoPcon *CB = new TGeoPcon("ITS SSD Thermal Centeral Stesalite "
                                "CylinderCB",0.0,360.0,6);
    TGeoTube *CC = new TGeoTube("ITS SSD Thermal Centeral Rohacell "
                                "CylinderCC",
                                CylROuter-CylCthick-CylFoamThick,
                                CylROuter-CylCthick,0.5*CylZFoamlength);
    CA->Z(0)    = -0.5*CylZlength;
    CA->Rmin(0) = CylRInner;
    CA->Rmax(0) = CylROuter;
    CA->Z(1)    = CA->GetZ(0) + CylZM6;
    CA->Rmin(1) = CA->GetRmin(0);
    CA->Rmax(1) = CA->GetRmax(0);
    CA->Z(2)    = -0.5*CylZFoamlength;
    CA->Rmin(2) = CylROuter - 2.0*CylCthick-CylFoamThick;
    CA->Rmax(2) = CA->GetRmax(0);
    CA->Z(3)    = -CA->GetZ(2);
    CA->Rmin(3) = CA->GetRmin(2);
    CA->Rmax(3) = CA->GetRmax(2);
    CA->Z(4)    = -CA->GetZ(1);
    CA->Rmin(4) = CA->GetRmin(1);
    CA->Rmax(4) = CA->GetRmax(1);
    CA->Z(5)    = -CA->GetZ(0);
    CA->Rmin(5) = CA->GetRmin(0);
    CA->Rmax(5) = CA->GetRmax(0);
    //
    CB->Z(0)    = CA->GetZ(0);
    CB->Rmin(0) = CA->GetRmin(0) + CylCthick;
    CB->Rmax(0) = CA->GetRmax(0) - CylCthick;
    CB->Z(1)    = CA->GetZ(1);
    CB->Rmin(1) = CA->GetRmin(1) + CylCthick;
    CB->Rmax(1) = CA->GetRmax(1) - CylCthick;
    CB->Z(2)    = CA->GetZ(2);
    CB->Rmin(2) = CA->GetRmin(2) + CylCthick;
    CB->Rmax(2) = CA->GetRmax(2) - CylCthick;
    CB->Z(3)    = CA->GetZ(3);
    CB->Rmin(3) = CA->GetRmin(3) + CylCthick;
    CB->Rmax(3) = CA->GetRmax(3) - CylCthick;
    CB->Z(4)    = CA->GetZ(4);
    CB->Rmin(4) = CA->GetRmin(4) + CylCthick;
    CB->Rmax(4) = CA->GetRmax(4) - CylCthick;
    CB->Z(5)    = CA->GetZ(5);
    CB->Rmin(5) = CA->GetRmin(5) + CylCthick;
    CB->Rmax(5) = CA->GetRmax(5) - CylCthick;
    //
    printPcon(CA);
    printPcon(CB);
    printTube(CC);
    //
    TGeoTube *CD = new TGeoTube("ITS SSD Thermal Centeral Cylinder M6 screwCD",
                                0.0,CylRM6,0.5*CylZM6);
    TGeoTube *CE = new TGeoTube("ITS SSD Thermal Centeral Cylinder PinCE",
                                0.0,CylRPin,0.5*CylZPin);
    //
    TGeoVolume *CAv,*CBv,*CCv,*CDv,*CEv;
    CAv = new TGeoVolume("ITSssdCentCylCA",CA,SSDcf);
    CAv->SetVisibility(kTRUE);
    CAv->SetLineColor(4); // blue
    CAv->SetLineWidth(1);
    CAv->SetFillColor(CAv->GetLineColor());
    CAv->SetFillStyle(4000); // 0% transparent
    CBv = new TGeoVolume("ITSssdCentCylCB",CB,SSDfs);
    CBv->SetVisibility(kTRUE);
    CBv->SetLineColor(2); // red
    CBv->SetLineWidth(1);
    CBv->SetFillColor(CBv->GetLineColor());
    CBv->SetFillStyle(4050); // 50% transparent
    CCv = new TGeoVolume("ITSssdCentCylCC",CC,SSDfo);
    CCv->SetVisibility(kTRUE);
    CCv->SetLineColor(3); // green
    CCv->SetLineWidth(1);
    CCv->SetFillColor(CCv->GetLineColor());
    CCv->SetFillStyle(4050); // 50% transparent
    CDv = new TGeoVolume("ITSssdCentCylCD",CD,SSDss);
    CDv->SetVisibility(kTRUE);
    CDv->SetLineColor(1); // black
    CDv->SetLineWidth(1);
    CDv->SetFillColor(CDv->GetLineColor());
    CDv->SetFillStyle(4000); // 0% transparent
    CEv = new TGeoVolume("ITSssdCentCylCE",CE,SSDss);
    CEv->SetVisibility(kTRUE);
    CEv->SetLineColor(1); // black
    CEv->SetLineWidth(1);
    CEv->SetFillColor(CEv->GetLineColor());
    CEv->SetFillStyle(4000); // 0% transparent
    // Insert Bolt and Pins in both the Cone and Cylinder at the same time.
    CBv->AddNode(CCv,1,0);
    CAv->AddNode(CBv,1,0);
    Moth->AddNode(CAv,1,0);
    if(GetDebug()){
        CAv->PrintNodes();
        CBv->PrintNodes();
        CCv->PrintNodes();
    } // end if
    //
    // SSD Cone
    // Data from Drawings ALR 0743/2E "Supporto Globale Settore SSD" and 
    // ALR 0743/2A "Supporto Generale Settore SSD".
    //
    const Double_t ConThick            = 13.0*kmm; // Thickness of Cone.
    const Double_t ConCthick           = 0.75*kmm; // Carbon finber thickness
    const Double_t ConRCurv            = 10.0*kmm; // Radius of curvature.
    const Double_t ConT                = 39.0*kDegree; // angle of SSD cone.
    const Double_t ConZOuterRing       = 47.0*kmm;
    const Double_t ConZOuterRingMill   = ConZOuterRing-5.0*kmm;
    const Double_t ConZToCylinder      = 170.0*kmm;
    const Double_t ConZLength          = 176.5*kmm-
                                          (ConZOuterRing-ConZOuterRingMill);
    const Double_t ConZInnerRing       = 161.5*kmm-
                                          (ConZOuterRing-ConZOuterRingMill);
    const Double_t ConZOuterRingInside = 30.25*kmm-
                                          (ConZOuterRing-ConZOuterRingMill);
    const Double_t ConZDisplacement    = ConZToCylinder + 0.5*CylZlength;
    const Double_t ConROuterMax        = 0.5*985.0*kmm;
    const Double_t ConROuterMin        = 0.5*945.0*kmm;
    const Double_t ConRCylOuterMill    = 0.5*597.0*kmm;
    const Double_t ConRInnerMin        = 0.5*564.0*kmm;
    const Double_t ConRCentCurv0       = 0.5*927.0*kmm;
    const Double_t ConRCentCurv1       = 0.5*593.0*kmm;
    //const Double_t ConRCentCurv2       = 0.5*578.0*kmm;
    // Foam core.
    const Double_t ConRohacellL0       = 112.3*kmm;
    const Double_t ConRohacellL1       = 58.4*kmm;
    // Screws and pins in outer SSD cone ring
    const Double_t ConROutHoles        = 0.5*965.0*kmm;
    const Double_t ConRScrewM5by12     = 0.5*5.0*kmm;
    const Double_t ConLScrewM5by12     = 0.5*12.0*kmm;
    const Int_t    ConNScrewM5by12     = 2;
    const Double_t ConRPinO6           = 0.5*6.0*kmm;
    const Double_t ConLPinO6           = 0.5*10.0*kmm;
    const Int_t    ConNPinO6           = 3;
    const Int_t    ConNRailScrews      = 4;
    const Int_t    ConNRailPins        = 2;
    const Int_t    ConNmounts          = 4;
    const Double_t ConMountPhi0        = 9.0*kDegree; // degrees
    // Holes in SSD cone, Ch* Cable Hole, Th* Tubing hole, and
    // Mh* mounting-post holes
    const Double_t ConCableHoleROut    = 0.5*920.0*kmm;
    const Double_t ConCableHoleRinner  = 0.5*800.0*kmm;
    const Double_t ConCableHoleWidth   = 200.0*kmm;
    const Double_t ConCableHoleAngle   = 42.0*kDegree;
    //const Double_t ConCableHolePhi0    = 90.0/4.0*kDegree;
    //const Int_t    ConNCableHoles      = 8;
    const Double_t ConCoolHoleWidth    = 40.0*kmm;
    const Double_t ConCoolHoleHight    = 30.0*kmm;
    const Double_t ConCoolHoleRmin     = 350.0*kmm;
    //const Double_t ConCoolHolephi0     = 90.0/4.0*kDegree;
    //const Int_t    ConNCoolHoles       = 8;
    const Double_t ConMountHoleWidth   = 20.0*kmm;
    const Double_t ConMountHoleHight   = 20.0*kmm;
    const Double_t ConMountHoleRmin    = 317.5*kmm;
    //const Double_t ConMountHolephi0    = 0.0*kDegree;
    //const Int_t    ConNMountHoles      = 6;
    // SSD cone Wings with holes.
    const Double_t ConWingRmax         = 527.5*kmm;
    const Double_t ConWingWidth        = 70.0*kmm;
    const Double_t ConWingThick        = 10.0*kmm;
    const Double_t ConWingPhi0         = 45.0*kDegree;
    //const Int_t    ConNWings           = 4;
    // SSD-SDD Thermal/Mechanical cylinder mounts
    const Double_t ConRM6Head          = 8.0*kmm;
    const Double_t ConZM6Head     = 8.5*kmm;
    //
    // SSD-SDD Mounting bracket
    const Double_t SupPRmin            = 0.5*539.0*kmm;// see SDD RoutMin
    const Double_t SupPRmax            = 0.5*585.0*kmm;
    const Double_t SupPZ               = 3.5*kmm;
    const Double_t SupPPhi1            = -0.5*70.0*kmm/SupPRmax*kRadian;
    const Double_t SupPPhi2            = -SupPPhi1;
    //
    const Double_t Sintc               = TMath::Sin(ConT*kRadian);
    const Double_t Costc               = TMath::Cos(ConT*kRadian);
    //
    // Lets start with the upper left outer carbon fiber surface.
    // Between za[2],rmaxa[2] and za[4],rmaxa[4] there is a curved section
    // given by rmaxa = rmaxa[2]-r*Sind(t) for 0<=t<=ConT and 
    // za = za[2] + r*Cosd(t) for 0<=t<=ConT. Simularly between za[1],rmina[1
    // and za[3],rmina[3] there is a curve section given by
    // rmina = rmina[1]-r*Sind(t) for 0<=t<=ConT and za = za[1]+r&Sind(t)
    // for t<=0<=ConT. These curves have been replaced by straight lines
    // between the equivelent points for simplicity.
    // Poly-cone Volume A0. Top part of SSD cone Carbon Fiber.
    TGeoPcon *A0 = new TGeoPcon("ITSssdSuportConeCarbonFiberSurfaceA0",
                               0.0,360.0,15);
    A0->Z(0)    = 0.0;
    A0->Rmin(0) = ConROuterMin;
    A0->Rmax(0) = ConROuterMax;
    A0->Z(1)    = ConZOuterRingInside-ConRCurv;
    A0->Rmin(1) = A0->GetRmin(0);
    A0->Rmax(1) = A0->GetRmax(0);
    A0->Z(2)    = ConZOuterRingInside;
    A0->Rmin(2) = ConROuterMin-ConRCurv;
    A0->Rmax(2) = A0->GetRmax(0);
    A0->Z(3)    = A0->GetZ(2);
    A0->Rmin(3) = -1000; // See Below
    A0->Rmax(3) = A0->GetRmax(0);
    A0->Z(4)    = ConZOuterRingMill-ConRCurv;
    A0->Rmin(4) = -1000; // See Below
    A0->Rmax(4) = A0->GetRmax(0);
    A0->Z(5)    = ConZOuterRingMill;
    A0->Rmin(5) = -1000; // See Below
    A0->Rmax(5) = A0->GetRmax(0) - ConRCurv;
    A0->Z(6)    = A0->GetZ(5);
    A0->Rmin(6) = -1000; // See Below
    A0->Rmax(6) = ConRCentCurv0;
    A0->Z(7)    = ConZOuterRingMill+ConRCurv*Sintc;
    A0->Rmin(7) = -1000; // See Below
    A0->Rmax(7) = ConRCentCurv0-ConRCurv*Costc;
    A0->Z(8)    = -1000; // See Below
    A0->Rmin(8) = ConRInnerMin;
    A0->Rmax(8) = -1000; // See Below
    A0->Z(9)    = ConZInnerRing;
    A0->Rmin(9) = -1000; // See Below
    A0->Rmax(9) = -1000; // See Below
    A0->Z(10)   = ConZInnerRing;
    A0->Rmin(10)= ConRInnerMin;
    A0->Rmax(10)= -1000; // See Below
    A0->Z(11)   = ConZLength-ConRCurv+ConRCurv*Costc;
    A0->Rmin(11)= ConRInnerMin;
    A0->Rmax(11)= ConRCentCurv1+ConRCurv*Sintc;
    A0->Z(12)   = ConZToCylinder;
    A0->Rmin(12)= ConRInnerMin;
    A0->Rmax(12)= -1000; // See Below
    A0->Z(13)   = ConZToCylinder;
    A0->Rmin(13)= ConRCylOuterMill;
    A0->Rmax(13)= -1000; // See Below
    A0->Z(14)   = -1000; // See Below
    A0->Rmin(14)= ConRCylOuterMill;
    A0->Rmax(14)= ConRCylOuterMill;
    // Compute values undefined above.
    RadiusOfCurvature(ConRCurv,0.0,A0->GetZ(9),A0->GetRmin(9),ConT,A0->Z(8),x);
    A0->Rmin(3) = RminFromZpCone(A0,8,90.-ConT,A0->GetZ(3),0.0);
    A0->Rmin(4) = RminFromZpCone(A0,3,90.-ConT,A0->GetZ(4),0.0);
    A0->Rmin(5) = RminFromZpCone(A0,3,90.-ConT,A0->GetZ(5),0.0);
    A0->Rmin(6) = A0->GetRmin(5);
    A0->Rmin(7) = RminFromZpCone(A0,3,90.-ConT,A0->GetZ(7),0.0);
    A0->Rmax(8) = RmaxFromZpCone(A0,4,90.-ConT,A0->GetZ(8),0.0);
    A0->Rmin(9) = RminFromZpCone(A0,3,90.-ConT,A0->GetZ(9),0.0);
    A0->Rmax(9) = RmaxFromZpCone(A0,4,90.-ConT,A0->GetZ(9),0.0);
    A0->Rmax(10)= RmaxFromZpCone(A0,4,90.-ConT,A0->GetZ(10),0.0);
    t = TMath::Tan((270.+ConT)*TMath::DegToRad());
    A0->Z(14)   = (ConRCylOuterMill-A0->GetRmax(4)+t*A0->GetZ(4))/t;
    A0->Rmax(12)= RmaxFrom2Points(A0,11,14,A0->GetZ(12));
    A0->Rmax(13)= RmaxFrom2Points(A0,11,14,A0->GetZ(13));
    printPcon(A0);
    //
    // Poly-cone Volume B. Stesalite inside volume A0.
    // Now lets define the Inserto Stesalite 4411w material volume.
    // Poly-cone Volume A0. Top part of SSD cone Carbon Fiber.
    TGeoPcon *B0 = new TGeoPcon("ITSssdSuportConeStaseliteB0",
                               0.0,360.0,15);
    //
    B0->Z(0)    = A0->GetZ(0);
    B0->Rmin(0) = A0->GetRmin(0) + ConCthick;
    B0->Rmax(0) = A0->GetRmax(0) - ConCthick;
    InsidePoint(A0,0,1,2,ConCthick,B0,1,kFALSE); // Rmin
    B0->Rmax(1) = B0->Rmax(0);
    InsidePoint(A0,1,2,3,ConCthick,B0,2,kFALSE); // Rmin
    B0->Rmax(2) = B0->Rmax(0);
    InsidePoint(A0,2,3,9,ConCthick,B0,3,kFALSE);
    B0->Rmax(3) = B0->Rmax(0);
    InsidePoint(A0,0,4,5,ConCthick,B0,4,kTRUE); // Rmax
    B0->Rmin(4) = -1000.; // see Bellow
    InsidePoint(A0,4,5,6,ConCthick,B0,5,kTRUE); // Rmax
    B0->Rmin(5) = -1000.; // see Bellow
    InsidePoint(A0,5,6,7,ConCthick,B0,6,kTRUE); // Rmax
    B0->Rmin(6) = -1000.; // see Bellow
    InsidePoint(A0,6,7,11,ConCthick,B0,7,kTRUE); // Rmax
    B0->Rmin(7) = -1000.; // see Bellow
    InsidePoint(A0,3,8,9,ConCthick,B0,8,kFALSE); // Rmin
    B0->Rmax(8) = -1000.; // see Bellow
    InsidePoint(A0,8,9,10,ConCthick,B0,9,kFALSE); // Rmin
    B0->Rmax(9) = -1000.; // see Bellow
    B0->Z(10)   = A0->GetZ(10) + ConCthick;
    B0->Rmin(10)= A0->GetRmin(10);
    B0->Rmax(10)= -1000.; // see Bellow
    InsidePoint(A0,7,11,14,ConCthick,B0,11,kTRUE); // Rmax
    B0->Rmin(11)= A0->GetRmin(10);
    B0->Z(2)    = A0->GetZ(12);
    B0->Rmin(12)= A0->GetRmin(12);
    B0->Rmax(12)= -1000.; // see Bellow
    B0->Z(13)   = A0->GetZ(13);
    B0->Rmin(13)= A0->GetRmin(13);
    B0->Rmax(13)= -1000.; // see Bellow
    B0->Z(14)   = A0->GetZ(14) - ConCthick;
    B0->Rmin(14)= A0->GetRmin(14);
    B0->Rmax(14)= B0->Rmin(14); // Close?
    B0->Rmin(4) = RminFrom2Points(B0,3,8,B0->GetZ(4));
    B0->Rmin(5) = RminFrom2Points(B0,3,8,B0->GetZ(5));
    B0->Rmin(6) = B0->GetRmin(5);
    B0->Rmin(7) = RminFrom2Points(B0,3,8,B0->GetZ(7));
    B0->Rmax(8) = RmaxFrom2Points(B0,7,11,B0->GetZ(8));
    B0->Rmax(9) = RmaxFrom2Points(B0,7,11,B0->GetZ(9));
    B0->Rmax(10)= B0->GetRmax(9);
    B0->Rmax(12)= RmaxFrom2Points(B0,11,14,B0->GetZ(12));
    B0->Rmax(13)= RmaxFrom2Points(B0,11,14,B0->GetZ(13));
    printPcon(B0);
    //
    // Poly-cone Volume C0. Foam inside volume A0.
    // Now lets define the Rohacell foam material volume.
    TGeoPcon *C0 = new TGeoPcon("ITSssdSuportConeRohacellC0",
                               0.0,360.0,4);
    C0->Z(1)    = B0->GetZ(7);
    C0->Rmax(1) = B0->GetRmax(7);
    C0->Rmin(1) = RminFrom2Points(B0,3,8,C0->GetZ(1));
    C0->Rmin(0) = C0->GetRmax(1);
    C0->Rmax(0) = C0->GetRmin(0);
    C0->Z(0)    = Zfrom2MinPoints(B0,3,8,C0->Rmin(0));
    C0->Z(3)    = C0->GetZ(0)+(ConThick-2.0*ConCthick+ConRohacellL0)*Costc;
    C0->Rmin(3) = C0->GetRmin(0)+(ConThick-2.0*ConCthick-ConRohacellL0)*Sintc;
    C0->Rmax(3) = C0->GetRmin(3);
    C0->Rmin(2) = C0->GetRmin(3);
    C0->Z(2)    = Zfrom2MinPoints(B0,3,8,C0->GetRmin(2));
    C0->Rmax(2) = RmaxFrom2Points(B0,4,11,C0->GetZ(2));
    printPcon(C0);
    //
    // Poly-cone Volume F.  Second Foam inside volume A0.
    // Now lets define the Rohacell foam material volume.
    TGeoPcon *F0 = new TGeoPcon("ITSssdSuportConeRohacellCF0",
                               0.0,360.0,4);
    F0->Z(2)    = B0->GetZ(8);
    F0->Rmin(2) = B0->GetRmin(8);
    F0->Rmax(2) = B0->GetRmax(8);
    F0->Z(0)    = F0->GetZ(2)-ConRohacellL1*Sintc;
    F0->Rmin(0) = F0->GetRmin(2)+ConRohacellL1*Costc;
    F0->Rmax(0) = F0->GetRmin(0);
    F0->Z(1)    = Zfrom2MaxPoints(B0,4,11,F0->GetRmax(0));
    F0->Rmax(1) = F0->GetRmax(0);
    F0->Rmin(1) = RminFrom2Points(B0,3,8,F0->GetZ(1));
    F0->Rmin(3) = F0->GetRmin(2)+(ConThick-2.0*ConCthick)*Costc;
    F0->Z(3)    = F0->GetZ(2)+(ConThick-2.0*ConCthick)*Sintc;
    F0->Rmax(3) = F0->GetRmin(3);
    printPcon(F0);
    // Holes for Cables to pass Through is created by the intersection
    // between a cone segment and an Arb8, One for the volume A0 and a
    // larger one for the volumes B0 and C0, so that the surface is covered
    // in carbon figer (volume A0).
    TGeoConeSeg *Ah1 = new TGeoConeSeg("ITSssdCableHoleAh1",
                                       0.5*ConZLength,ConCableHoleRinner,
                                       ConCableHoleROut,ConCableHoleRinner,
                                       ConCableHoleROut,
                                       90.-0.5*ConCableHoleWidth/
                                       ConCableHoleROut/kRadian,
                                       90.+0.5*ConCableHoleWidth/
                                       ConCableHoleROut/kRadian);
    TGeoConeSeg *Bh1 = new TGeoConeSeg("ITSssdCableHoleBh1",0.5*ConZLength,
                                       ConCableHoleRinner-ConCthick,
                                       ConCableHoleROut+ConCthick,
                                       ConCableHoleRinner-ConCthick,
                                       ConCableHoleROut+ConCthick,
                                       90.+((-0.5*ConCableHoleWidth-ConCthick)/
                                        (ConCableHoleROut-ConCthick))/kRadian,
                                       90.+((+0.5*ConCableHoleWidth-ConCthick)/
                                        (ConCableHoleROut-ConCthick))/kRadian);
    x0 = Ah1->GetRmax1()*TMath::Cos(Ah1->GetPhi2()*kRadian);
    y0 = Ah1->GetRmax1()*TMath::Sin(Ah1->GetPhi2()*kRadian);
    TGeoArb8 *Ah2 = new TGeoArb8("ITSssdCableHoleAh2",0.5*ConZLength);
    y  = Ah1->GetRmax1();
    x  = x0+(y-y0)/TMath::Tan((90.0+ConCableHoleAngle)*kRadian);
    Ah2->SetVertex(0,x,y);
    y  = Ah1->GetRmin1()*TMath::Sin(Ah1->GetPhi2()*kRadian);
    x  = x0+(y-y0)/TMath::Tan((90.0+ConCableHoleAngle)*kRadian);
    Ah2->SetVertex(3,x,y);
    x0 = Ah1->GetRmax1()*TMath::Cos(Ah1->GetPhi1()*kRadian);
    y0 = Ah1->GetRmax1()*TMath::Sin(Ah1->GetPhi1()*kRadian);
    y  = Ah1->GetRmax1();
    x  = x0+(y-y0)/TMath::Tan((90.0-ConCableHoleAngle)*kRadian);
    Ah2->SetVertex(1,x,y);
    y  = Ah1->GetRmin1()*TMath::Sin(Ah1->GetPhi1()*kRadian);
    x  = x0+(y-y0)/TMath::Tan((90.0-ConCableHoleAngle)*kRadian);
    Ah2->SetVertex(2,x,y);
    //
    x0 = Bh1->GetRmax1()*TMath::Cos(Bh1->GetPhi2()*kRadian);
    y0 = Bh1->GetRmax1()*TMath::Sin(Bh1->GetPhi2()*kRadian);
    TGeoArb8 *Bh2 = new TGeoArb8("ITSssdCableHoleBh2",0.5*ConZLength);
    y  = Bh1->GetRmax1();
    x  = x0+(y-y0)/TMath::Tan((90.0+ConCableHoleAngle)*kRadian);
    Bh2->SetVertex(0,x,y);
    y  = Bh1->GetRmin1()*TMath::Sin(Bh1->GetPhi2()*kRadian);
    x  = x0+(y-y0)/TMath::Tan((90.0+ConCableHoleAngle)*kRadian);
    Bh2->SetVertex(3,x,y);
    x0 = Bh1->GetRmax1()*TMath::Cos(Bh1->GetPhi1()*kRadian);
    y0 = Bh1->GetRmax1()*TMath::Sin(Bh1->GetPhi1()*kRadian);
    y  = Bh1->GetRmax1();
    x  = x0+(y-y0)/TMath::Tan((90.0-ConCableHoleAngle)*kRadian);
    Bh2->SetVertex(1,x,y);
    y  = Bh1->GetRmin1()*TMath::Sin(Bh1->GetPhi1()*kRadian);
    x  = x0+(y-y0)/TMath::Tan((90.0-ConCableHoleAngle)*kRadian);
    Bh2->SetVertex(2,x,y);
    for(i=0;i<4;i++){ // define points at +dz
        Ah2->SetVertex(i+4,(Ah2->GetVertices())[2*i],
                           (Ah2->GetVertices())[1+2*i]);
        Bh2->SetVertex(i+4,(Bh2->GetVertices())[2*i],
                           (Bh2->GetVertices())[1+2*i]);
    } // end for i
    TGeoBBox *Ah3 = new TGeoBBox("ITSssdCoolingHoleAh3",0.5*ConCoolHoleWidth,
                                 0.5*ConCoolHoleHight,0.5*ConZLength);
    TGeoBBox *Bh3 = new TGeoBBox("ITSssdCoolingHoleBh3",
                                 0.5*ConCoolHoleWidth+ConCthick,
                                 0.5*ConCoolHoleHight+ConCthick,
                                 0.5*ConZLength);
    TGeoBBox *Ah4 = new TGeoBBox("ITSssdMountingPostHoleAh4",
                                 0.5*ConMountHoleWidth,
                                 0.5*ConMountHoleHight,0.5*ConZLength);
    TGeoBBox *Bh4 = new TGeoBBox("ITSssdMountingPostHoleBh4",
                                 0.5*ConMountHoleWidth+ConCthick,
                                 0.5*ConMountHoleHight+ConCthick,
                                 0.5*ConZLength);
    printConeSeg(Ah1);
    printConeSeg(Bh1);
    printArb8(Ah2);
    printArb8(Bh2);
    printBBox(Ah3);
    printBBox(Bh3);
    printBBox(Ah4);
    printBBox(Bh4);
    // SSD Cone Wings
    TGeoConeSeg *G = new TGeoConeSeg("ITSssdWingCarbonFiberSurfaceG",
                                     0.5*ConWingThick,ConROuterMax-ConCthick,
                                     ConWingRmax,
                                     ConROuterMax-ConCthick,ConWingRmax,
                          ConWingPhi0-0.5*ConWingWidth/ConWingRmax*kRadian,
                          ConWingPhi0+0.5*ConWingWidth/ConWingRmax*kRadian);
    TGeoConeSeg *H = new TGeoConeSeg("ITSssdWingStaseliteH",
                 0.5*ConWingThick-ConCthick,ConROuterMax-ConCthick,
                                     ConWingRmax-ConCthick,
                                     ConROuterMax-ConCthick,
                                     ConWingRmax-ConCthick,
                                     ConWingPhi0-((0.5*ConWingWidth-ConCthick)/
                                             (ConWingRmax-ConCthick))*kRadian,
                                     ConWingPhi0+((0.5*ConWingWidth-ConCthick)/
                    (ConWingRmax-ConCthick))*kRadian);
    printConeSeg(G);
    printConeSeg(H);
    // SDD support plate, SSD side.
    //Poly-cone Volume T.
    TGeoTubeSeg *T = new TGeoTubeSeg("ITSssdsddMountingBracketT",
                                     SupPRmin,SupPRmax,
                                     SupPZ,SupPPhi1,
                                     SupPPhi2);
    printTubeSeg(T);
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
    vl[0] = 0.0;vl[1] = ConCoolHoleRmin+0.5*ConCoolHoleHight;vl[2] = 0.0;
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
    vl[0] = ConMountHoleRmin+0.5*ConMountHoleHight; vl[1] = 0.0; vl[2] = 0.0;
    rotZ30->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA30 = new TGeoCombiTrans("ITSssdConeTZ30",vl[0],
                                                      vl[1],vl[2],rotZ30);
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
    vl[0] = 0.0; vl[1] = 0.0; vl[2] = A0->GetZ(10)+T->GetDz();
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
    TGeoCompositeShape *A = new TGeoCompositeShape(
        "ITSssdSuportConeCarbonFiberSurfaceA","(((((((((((((((((((((((((((("
        "ITSssdSuportConeCarbonFiberSurfaceA0 +"
        "ITSssdWingCarbonFiberSurfaceG) +"
        "ITSssdWingCarbonFiberSurfaceG:ITSssdConeZ90) +"
        "ITSssdWingCarbonFiberSurfaceG:ITSssdConeZ180) +"
        "ITSssdWingCarbonFiberSurfaceG:ITSssdConeZ270) -"
        "(ITSssdCableHoleAh1*ITSssdCableHoleAh2):ITSssdConeZ225) -"
        "(ITSssdCableHoleAh1*ITSssdCableHoleAh2):ITSssdConeZ675) -"
        "(ITSssdCableHoleAh1*ITSssdCableHoleAh2):ITSssdConeZ1125) -"
        "(ITSssdCableHoleAh1*ITSssdCableHoleAh2):ITSssdConeZ1575) -"
        "(ITSssdCableHoleAh1*ITSssdCableHoleAh2):ITSssdConeZ2025) -"
        "(ITSssdCableHoleAh1*ITSssdCableHoleAh2):ITSssdConeZ2475) -"
        "(ITSssdCableHoleAh1*ITSssdCableHoleAh2):ITSssdConeZ2925) -"
        "(ITSssdCableHoleAh1*ITSssdCableHoleAh2):ITSssdConeZ3375) -"
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
    TGeoCompositeShape *B = new TGeoCompositeShape(
        "ITSssdSuportConeStaseliteB","(((((((((((((((((((((((((((("
        "ITSssdSuportConeStaseliteB0 +"
        "ITSssdWingStaseliteH) +"
        "ITSssdWingStaseliteH:ITSssdConeZ90) +"
        "ITSssdWingStaseliteH:ITSssdConeZ180) +"
        "ITSssdWingStaseliteH:ITSssdConeZ270) -"
        "(ITSssdCableHoleBh1*ITSssdCableHoleBh2):ITSssdConeZ225) -"
        "(ITSssdCableHoleBh1*ITSssdCableHoleBh2):ITSssdConeZ675) -"
        "(ITSssdCableHoleBh1*ITSssdCableHoleBh2):ITSssdConeZ1125) -"
        "(ITSssdCableHoleBh1*ITSssdCableHoleBh2):ITSssdConeZ1575) -"
        "(ITSssdCableHoleBh1*ITSssdCableHoleBh2):ITSssdConeZ2025) -"
        "(ITSssdCableHoleBh1*ITSssdCableHoleBh2):ITSssdConeZ2475) -"
        "(ITSssdCableHoleBh1*ITSssdCableHoleBh2):ITSssdConeZ2925) -"
        "(ITSssdCableHoleBh1*ITSssdCableHoleBh2):ITSssdConeZ3375) -"
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
    TGeoCompositeShape *C = new TGeoCompositeShape(
      "ITSssdSuportConeRohacellC","("
      "ITSssdSuportConeRohacellC0 -((((((("
      "ITSssdCableHoleBh1:ITSssdConeZ225*ITSssdCableHoleBh2:ITSssdConeZ225)-"
      "ITSssdCableHoleBh1:ITSssdConeZ675*ITSssdCableHoleBh2:ITSssdConeZ675)-"
      "ITSssdCableHoleBh1:ITSssdConeZ1125*ITSssdCableHoleBh2:ITSssdConeZ1125)-"
      "ITSssdCableHoleBh1:ITSssdConeZ1575*ITSssdCableHoleBh2:ITSssdConeZ1575)-"
      "ITSssdCableHoleBh1:ITSssdConeZ2025*ITSssdCableHoleBh2:ITSssdConeZ2025)-"
      "ITSssdCableHoleBh1:ITSssdConeZ2475*ITSssdCableHoleBh2:ITSssdConeZ2475)-"
      "ITSssdCableHoleBh1:ITSssdConeZ2925*ITSssdCableHoleBh2:ITSssdConeZ2925))"
        );
    TGeoCompositeShape *F = new TGeoCompositeShape(
        "ITSssdSuportConeRohacellCF","((((("
        "ITSssdSuportConeRohacellCF0 -("
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ30) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ90) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ150) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ210) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ270) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ330)"
        );
    //
    // In volume SCB, th Inserto Stesalite 4411w material volume, there
    // are a number of Stainless steel screw and pin studs which will be
    // filled with screws/studs.
    TGeoTube *D = new TGeoTube("ITS Screw+stud used to mount things to "
                       "the SSD support cone",
                               0.0,ConRScrewM5by12,ConLScrewM5by12);
    printTube(D);
    TGeoTube *E = new TGeoTube("ITS pin used to mount things to the "
                       "SSD support cone",0.0,ConRPinO6,ConLPinO6);
    printTube(E);
    // Bolt heads holding the SSD-SDD tube to the SSD cone.
    // Bolt -- PolyCone
    //Poly-cone Volume Q.
    TGeoPcon *Q = new TGeoPcon("ITS SSD Thermal sheal M6 screw headQ",
                               0.0,360.0,4);
    Q->Z(0)    = A0->GetZ(12);
    Q->Rmin(0) = 0.0;
    Q->Rmax(0) = CylRM6;
    Q->Z(1)    = Q->GetZ(0) + ConZM6Head;
    Q->Rmin(1) = 0.0;
    Q->Rmax(1) = CylRM6;
    Q->Z(2)    = Q->GetZ(1);
    Q->Rmin(2) = 0.0;
    Q->Rmax(2) = ConRM6Head;
    Q->Z(3)    = Q->GetZ(0)-SupPZ;
    Q->Rmin(3) = 0.0;
    Q->Rmax(3) = 0.5*ConRM6Head;
    printPcon(Q);
    // air infront of bolt (stasolit Volume K) -- Tube
    TGeoTube *R = new TGeoTube("ITS Air in front of bolt (in stasolit)R",
                               Q->GetRmin(3),Q->GetRmax(3),
                               0.5*(SupPZ-ConCthick));
    // air infront of bolt (carbon fiber volume I) -- Tube
    TGeoTube *S = new TGeoTube("ITS Air in front of Stainless Steal Screw "
                               "end, M6S",Q->GetRmin(3),Q->GetRmax(3),
                               0.5*ConCthick);
    printTube(S);
    //
    TGeoVolume *Av,*Bv,*Cv,*Dv,*Ev,*Fv,*Qv,*Rv,*Sv,*Tv;
    //
    Av = new TGeoVolume("ITSssdConeA",A,SSDcf); // Carbon Fiber
    Av->SetVisibility(kTRUE);
    Av->SetLineColor(4); // blue
    Av->SetLineWidth(1);
    Av->SetFillColor(Av->GetLineColor());
    Av->SetFillStyle(4000); // 0% transparent
    Bv = new TGeoVolume("ITSssdConeB",B,SSDfs); // Staselite
    Bv->SetVisibility(kTRUE);
    Bv->SetLineColor(2); // red
    Bv->SetLineWidth(1);
    Bv->SetFillColor(Bv->GetLineColor());
    Bv->SetFillStyle(4010); // 10% transparent
    Cv = new TGeoVolume("ITSssdConeC",C,SSDfo); // Rohacell
    Cv->SetVisibility(kTRUE);
    Cv->SetLineColor(3); // green
    Cv->SetLineWidth(1);
    Cv->SetFillColor(Cv->GetLineColor());
    Cv->SetFillStyle(4050); // 50% transparent
    Fv = new TGeoVolume("ITSssdConeF",F,SSDfo); // Rohacell;
    Fv->SetVisibility(kTRUE);
    Fv->SetLineColor(3); // green
    Fv->SetLineWidth(1);
    Fv->SetFillColor(Fv->GetLineColor());
    Fv->SetFillStyle(4050); // 50% transparent
    Dv = new TGeoVolume("ITSssdConeD",D,SSDss);
    Dv->SetVisibility(kTRUE);
    Dv->SetLineColor(1); // black
    Dv->SetLineWidth(1);
    Dv->SetFillColor(Dv->GetLineColor());
    Dv->SetFillStyle(4000); // 0% transparent
    Ev = new TGeoVolume("ITSssdConeE",E,SSDss);
    Ev->SetVisibility(kTRUE);
    Ev->SetLineColor(1); // black
    Ev->SetLineWidth(1);
    Ev->SetFillColor(Ev->GetLineColor());
    Ev->SetFillStyle(4000); // 0% transparent
    Qv = new TGeoVolume("ITSssdConeQ",Q,SSDss);
    Qv->SetVisibility(kTRUE);
    Qv->SetLineColor(1); // black
    Qv->SetLineWidth(1);
    Qv->SetFillColor(Qv->GetLineColor());
    Qv->SetFillStyle(4000); // 00% transparent
    Rv = new TGeoVolume("ITSssdConeR",R,SSDair);
    Rv->SetVisibility(kTRUE);
    Rv->SetLineColor(5); // yellow
    Rv->SetLineWidth(1);
    Rv->SetFillColor(Rv->GetLineColor());
    Rv->SetFillStyle(4090); // 90% transparent
    Sv = new TGeoVolume("ITSssdConeS",S,SSDair);
    Sv->SetVisibility(kTRUE);
    Sv->SetLineColor(5); // yellow
    Sv->SetLineWidth(1);
    Sv->SetFillColor(Sv->GetLineColor());
    Sv->SetFillStyle(4090); // 90% transparent
    Tv = new TGeoVolume("ITSssdsddMountingBracket",S,SSDal);
    Tv->SetVisibility(kTRUE);
    Tv->SetLineColor(5); // yellow
    Tv->SetLineWidth(1);
    Tv->SetFillColor(Tv->GetLineColor());
    Tv->SetFillStyle(4000); // 0% transparent
    //
    TGeoCombiTrans *rotran;
    TGeoTranslation *tran;
    tran = new TGeoTranslation("ITSssdConeTrans",0.0,0.0,-ConZDisplacement);
    TGeoRotation *rotY180 = new TGeoRotation("",0.0,180.0,0.0);
    TGeoCombiTrans *flip  = new TGeoCombiTrans("ITSssdConeFlip",
                                           0.0,0.0,ConZDisplacement,rotY180);
    delete rotY180;// rot not explicity used in AddNode functions.
    //
    //
    //
    //
    Av->AddNode(Bv,1,0);
    Bv->AddNode(Cv,1,0);
    Bv->AddNode(Fv,1,0);
    Moth->AddNode(Av,1,tran); // RB24 side
    Moth->AddNode(Av,2,flip); // RB26 side (Absorber)
    //
    //
    //
    // Insert Bolt and Pins in both the Cone and Cylinder at the same time.
    Int_t NcopyCDv=0,NcopyCEv=0,NcopyQv=0,NcopyRv=0,NcopySv=0,NcopyTv=0;
    Int_t NcopyDv=0,NcopyEv=0;
    z = CB->GetZ(0)-0.5*CylZPin;
    dt = (360.0/((Double_t)CylNPin));
    for(i=0;i<CylNPin;i++){
        t = ((Double_t)i)*dt;
        x = CylRholes*TMath::Cos((t+CylPhi0Pin)*kRadian);
        y = CylRholes*TMath::Sin((t+CylPhi0Pin)*kRadian);
        tran = new TGeoTranslation("",x,y,z);
        CBv->AddNode(CDv,++NcopyCDv,tran);
        tran = new TGeoTranslation("",x,y,-z);
        CBv->AddNode(CDv,++NcopyCDv,tran);
    } // end for i
    dt = (360.0/((Double_t)CylNM6));
    for(i=0;i<CylNM6;i++){
        t = ((Double_t)i)*dt;
        x = CylRholes*TMath::Cos((t+CylPhi0M6)*kRadian);
        y = CylRholes*TMath::Sin((t+CylPhi0M6)*kRadian);
        z = CB->GetZ(0)-0.5*CylZM6;
        tran = new TGeoTranslation("",x,y,z);
        CBv->AddNode(CEv,++NcopyCEv,tran);
        tran = new TGeoTranslation("",x,y,-z);
        CBv->AddNode(CEv,++NcopyCEv,tran);
        tran = new TGeoTranslation("",x,y,0.0);
        Bv->AddNode(Qv,++NcopyQv,tran);
        if(!((t<rotranBrTZ60->GetRotation()->GetPhiRotation()+T->GetPhi2()&&
              t>rotranBrTZ60->GetRotation()->GetPhiRotation()-T->GetPhi1())||
             (t<rotranBrTZ180->GetRotation()->GetPhiRotation()+T->GetPhi2()&&
              t>rotranBrTZ180->GetRotation()->GetPhiRotation()-T->GetPhi1())||
             (t<rotranBrTZ300->GetRotation()->GetPhiRotation()+T->GetPhi2()&&
              t>rotranBrTZ300->GetRotation()->GetPhiRotation()-T->GetPhi1()))){
            // If not at an angle where the bracket T is located.
            tran = new TGeoTranslation("",x,y,B0->GetZ(10)-R->GetDz());
            Bv->AddNode(Rv,++NcopyRv,tran);
            tran = new TGeoTranslation("",x,y,A0->GetZ(10)-S->GetDz());
            Av->AddNode(Sv,++NcopySv,tran);
        } // end if
    } // end for i
    // Add the mounting brackets to the RB24 side only.
    vl[0] = 0.0; vl[1] = 0.0, vl[2] = A0->GetZ(10)+ConZDisplacement-T->GetDz();
    rotZ60->LocalToMaster(vl,vg);
    rotran = new TGeoCombiTrans("",vg[0],vg[1],vg[2],rotZ60);
    Moth->AddNode(Tv,++NcopyTv,rotran);
    rotZ180->LocalToMaster(vl,vg);
    rotran = new TGeoCombiTrans("",vg[0],vg[1],vg[2],rotZ180);
    Moth->AddNode(Tv,++NcopyTv,rotran);
    rotZ300->LocalToMaster(vl,vg);
    rotran = new TGeoCombiTrans("",vg[0],vg[1],vg[2],rotZ300);
    Moth->AddNode(Tv,++NcopyTv,rotran);
    //
    Double_t da[] = {-3.5,-1.5,1.5,3.5};
    for(i=0;i<2;i++){ // Mounting for ITS-TPC bracket or ITS-Rails
        t0 = 180.*((Double_t)i)*kRadian;
        for(j=-ConNScrewM5by12/2;j<=ConNScrewM5by12/2;j++)if(j!=0){
                    //screws per ITS-TPC brkt
            t = t0 + 5.0*((Double_t)j)*kRadian;
            tran = new TGeoTranslation("",ConROutHoles*TMath::Cos(t),
                                          ConROutHoles*TMath::Sin(t),
                                          B0->GetZ(0)+D->GetDz());
            Bv->AddNode(Dv,++NcopyDv,tran);
        } // end or j
        for(j=-ConNPinO6/2;j<=ConNPinO6/2;j++){ // pins per ITS-TPC bracket
            t = t0 + 3.0*((Double_t)j)*kRadian;
            tran = new TGeoTranslation("",ConROutHoles*TMath::Cos(t),
                                          ConROutHoles*TMath::Sin(t),
                                          B0->GetZ(0)+D->GetDz());
            Bv->AddNode(Ev,++NcopyEv,tran);
        } // end or j
        t0 = (96.5+187.*((Double_t)i))*kRadian;
        for(j=0;j<ConNRailScrews;j++){ // screws per ITS-rail bracket
            t = t0+da[j]*kRadian;
            tran = new TGeoTranslation("",ConROutHoles*TMath::Cos(t),
                                          ConROutHoles*TMath::Sin(t),
                                          B0->GetZ(0)+D->GetDz());
            Bv->AddNode(Dv,++NcopyDv,tran);
        } // end or j
        t0 = (91.5+184.*((Double_t)i))*kRadian;
        for(j=-ConNRailPins/2;j<=ConNRailPins/2;j++)if(j!=0){ 
             // pins per ITS-rail bracket
            t = t0+(7.0*((Double_t)j))*kRadian;
            tran = new TGeoTranslation("",ConROutHoles*TMath::Cos(t),
                                          ConROutHoles*TMath::Sin(t),
                                          B0->GetZ(0)+D->GetDz());
            Bv->AddNode(Ev,++NcopyEv,tran);
        } // end or j
    } // end for i
    for(i=0;i<ConNmounts;i++){ 
                // mounting points for SPD-cone+Beam-pipe support
        t0 = (45.0+((Double_t)i)*360./((Double_t)ConNmounts))*kRadian;
        for(j=-1;j<=1;j++)if(j!=0){ // 2 screws per bracket
            t = t0+((Double_t)j)*0.5*ConMountPhi0*kRadian;
            tran = new TGeoTranslation("",ConROutHoles*TMath::Cos(t),
                                          ConROutHoles*TMath::Sin(t),
                                          B0->GetZ(0)+D->GetDz());
            Bv->AddNode(Dv,++NcopyDv,tran);
        } // end for j
        for(j=0;j<1;j++){ // 1 pin per bracket
            t = t0;
            tran = new TGeoTranslation("",ConROutHoles*TMath::Cos(t),
                                          ConROutHoles*TMath::Sin(t),
                                          B0->GetZ(0)+D->GetDz());
            Bv->AddNode(Ev,++NcopyEv,tran);
        } // end for j
    } // end for i
    if(GetDebug()){
        Av->PrintNodes();
        Bv->PrintNodes();
        Cv->PrintNodes();
        Dv->PrintNodes();
        Ev->PrintNodes();
        Fv->PrintNodes();
        Qv->PrintNodes();
        Rv->PrintNodes();
        Sv->PrintNodes();
        Tv->PrintNodes();
    } // end if
}

//______________________________________________________________________
void AliITSv11GeometrySupport::ServicesCableSupport(TGeoVolume *Moth){
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
    TGeoMedium *SUPcf    = 0; // SUP support cone Carbon Fiber materal number.
    TGeoMedium *SUPfs    = 0; // SUP support cone inserto stesalite 4411w.
    TGeoMedium *SUPfo    = 0; // SUP support cone foam, Rohacell 50A.
    TGeoMedium *SUPss    = 0; // SUP support cone screw material,Stainless
    TGeoMedium *SUPair   = 0; // SUP support cone Air
    TGeoMedium *SUPal    = 0; // SUP support cone SDD mounting bracket Al
    TGeoMedium *SUPwater = 0; // SUP support cone Water
    TGeoManager *mgr = gGeoManager;
    SUPcf    = mgr->GetMedium("ITSssdCarbonFiber");
    SUPfs    = mgr->GetMedium("ITSssdStaselite4411w");
    SUPfo    = mgr->GetMedium("ITSssdRohacell50A");
    SUPss    = mgr->GetMedium("ITSssdStainlessSteal");
    SUPair   = mgr->GetMedium("ITSssdAir");
    SUPal    = mgr->GetMedium("ITSssdAl");
    SUPwater = mgr->GetMedium("ITSssdWater");
    //
    Int_t i,j;
    Double_t x,y,z,t,t0,dt,di,r;

    // RB 24 side
    const Double_t Z024         = 900*kmm;//SSup_203A.jpg
    const Double_t ThssFrame24  = 5.0*kmm;
    const Double_t RssFrame24   = 444.5*kmm-ThssFrame24; // SSup_204A.jpg
    const Double_t WidthFrame24 = 10.0*kmm;
    const Double_t HightFrame24 = 10.0*kmm;
    const Double_t Phi0Frame24  = 15.2*kDegree; // SSup_602A.jpg
    const Double_t Phi1Frame24  = (90.0-7.6)*kDegree; // SSup_802A.jpg
    const Double_t ZssFrameSection24 = (415.0-10.0)*kmm;
    const Int_t    NZsections24      = 4;
    const Int_t    NPhiSections24    = 4;
    const Int_t    NFramesPhi24      = 4;
    //
    TGeoTubeSeg *M24 = new TGeoTubeSeg("ITS sup Cable tray support frame "
                                       "mother volume M24",
                                       RssFrame24,RssFrame24+ThssFrame24,
                                   0.5*(4.*ZssFrameSection24+5*WidthFrame24),
                                       Phi0Frame24,Phi1Frame24);
    TGeoTubeSeg *A24 = new TGeoTubeSeg("ITS sup Cable tray support frame "
                                       "radial section A24",
                          RssFrame24,RssFrame24+ThssFrame24,0.5*WidthFrame24,
                                       Phi0Frame24,Phi1Frame24);
    TGeoBBox *B24 = new TGeoBBox("ITS sup Cable tray support frame Z section "
                                 "B24",
                    0.5*ThssFrame24,0.5*HightFrame24,0.5*ZssFrameSection24);
    printTubeSeg(A24);
    printTubeSeg(M24);
    printBBox(B24);
    TGeoVolume *A24v,*B24v,*M24v;
    TGeoTranslation *tran;
    TGeoRotation    *rot;
    TGeoCombiTrans  *tranrot;
    //
    A24v = new TGeoVolume("ITSsupFrameA24",A24,SUPss);
    A24v->SetVisibility(kTRUE);
    A24v->SetLineColor(1); // black
    A24v->SetLineWidth(1);
    A24v->SetFillColor(A24v->GetLineColor());
    A24v->SetFillStyle(4000); // 0% transparent
    B24v = new TGeoVolume("ITSsupFrameB24",B24,SUPss);
    B24v->SetVisibility(kTRUE);
    B24v->SetLineColor(1); // black
    B24v->SetLineWidth(1);
    B24v->SetFillColor(B24v->GetLineColor());
    B24v->SetFillStyle(4000); // 0% transparent
    M24v = new TGeoVolume("ITSsupFrameM24",M24,SUPair);
    M24v->SetVisibility(kTRUE);
    M24v->SetLineColor(7); // light blue
    M24v->SetLineWidth(1);
    M24v->SetFillColor(M24v->GetLineColor());
    M24v->SetFillStyle(4090); // 90% transparent
    //
    Int_t NcA24=1,NcB24=1;
    t0 = Phi0Frame24;
    dt = (Phi1Frame24-Phi0Frame24)/((Double_t)NPhiSections24);
    for(i=0;i<=NZsections24;i++){
        di = (Double_t) i;
        z = -M24->GetDz()+A24->GetDz() + di*(ZssFrameSection24+WidthFrame24);
        tran = new TGeoTranslation("",0.0,0.0,z);
        M24v->AddNode(A24v,NcA24++,tran);
        r = RssFrame24+B24->GetDX();
        z = z + A24->GetDz()+B24->GetDZ();
       if(i<NZsections24) for(j=0;j<=NPhiSections24;j++){
            t = t0 + ((Double_t)j)*dt;
            rot = new TGeoRotation("",0.0,0.0,t);
            y = r*TMath::Sin(t*kRadian);
            x = r*TMath::Cos(t*kRadian);
            tranrot = new TGeoCombiTrans("",x,y,z,rot);
            delete rot;// rot not explicity used in AddNode functions.
            M24v->AddNode(B24v,NcB24++,tranrot);
        } // end for j
    } // end for i
    tran = new TGeoTranslation("",0.0,0.0,Z024+M24->GetDz());
    Moth->AddNode(M24v,1,tran);
    for(i=1;i<NFramesPhi24;i++){
        di = (Double_t) i;
        rot = new TGeoRotation("",0.0,0.0,90.0*di);
        tranrot = new TGeoCombiTrans("",0.0,0.0,Z024+M24->GetDz(),rot);
        delete rot;// rot not explicity used in AddNode functions.
        Moth->AddNode(M24v,i+1,tranrot);
    } // end for i
    if(GetDebug()){
        A24v->PrintNodes();
        B24v->PrintNodes();
        M24v->PrintNodes();
    } // end if
    // Cable support tray 
    // Material is Aluminum
    //const Double_t RS24in     = TMath::Max(RssFrame24,444.5*kmm);
                                           // SSup_204A & SSup_206A
    //const Double_t RS24Airout = 459.5*kmm; // SSup_204A & SSup_206A
    //const Double_t RS24out    = 494.5*kmm; // SSup_206A & SSup_204A
    //const Double_t RS24PPout  = 550.0*kmm; // SSup_206A
    const Double_t LS24PP     = 350.0*kmm; // SSup_202A
    const Double_t LS24       = (2693.0-900.0)*kmm; //SSup_205A & SSup_207A
    const Double_t ThS24wall  = 1.0*kmm; // SSup_209A & SSup_210A
    const Double_t WbS24      = 42.0*kmm; // SSup_209A & SSup_210A
    //const Double_t WtS24      = 46.9*kmm; // SSup_209A & SSup_210A
    const Double_t WcapS24    = 50.0*kmm; // SSup_209A & SSup_210A
    //const Double_t WdS24      = 41.0*kmm; //SSup_209A ? should be 41.46938776
    const Double_t HS24       = 50.0*kmm; // SSup_209A & SSup_210A
    const Double_t OutDcoolTub= 12.0*kmm; // SSup_209A
    const Double_t InDcoolTub = 10.0*kmm; // SSup_209A
    const Double_t BlkNozInDS24= 6.0*kmm; // SSup_209A
    // The following are deduced or guessed at
    //const Double_t LtopLipS24 = 6.0*kmm; // Guessed at.
    //const Double_t LdLipS24   = 6.0*kmm; // Guessed at.
    //const Double_t HdS24      = OutDcoolTub; //
    const Double_t BlkNozZS24 = 6.0*kmm; // Guessed at.
    // Simplifided exterior shape. The side wall size is 2.5*thicker than
    // it should be (due to simplification).
    TGeoArb8 *C24 = new TGeoArb8("ITS Sup Cable Tray Element C24",0.5*LS24);
    C24->SetVertex(0,-0.5*WcapS24,HS24+ThS24wall);
    C24->SetVertex(1,+0.5*WcapS24,HS24+ThS24wall);
    C24->SetVertex(2,+0.5*WbS24,0.0);
    C24->SetVertex(3,-0.5*WbS24,0.0);
    C24->SetVertex(4,-0.5*WcapS24,HS24+ThS24wall);
    C24->SetVertex(5,+0.5*WcapS24,HS24+ThS24wall);
    C24->SetVertex(6,+0.5*WbS24,0.0);
    C24->SetVertex(7,-0.5*WbS24,0.0);
    TGeoArb8 *D24 = new TGeoArb8("ITS Sup Cable Tray lower Element D24",
                                 0.5*LS24);
    // Because of question about the value of WdS24, compute what it
    // should be assuming cooling tube fixes hight of volume.
    x = OutDcoolTub*(0.5*WcapS24-0.5*WbS24-ThS24wall)/(HS24-ThS24wall);
    D24->SetVertex(0,-x,OutDcoolTub+ThS24wall);
    D24->SetVertex(1,+x,OutDcoolTub+ThS24wall);
    D24->SetVertex(2,+0.5*WbS24-ThS24wall,ThS24wall);
    D24->SetVertex(3,-0.5*WbS24+ThS24wall,ThS24wall);
    D24->SetVertex(4,-x,OutDcoolTub+ThS24wall);
    D24->SetVertex(5,+x,OutDcoolTub+ThS24wall);
    D24->SetVertex(6,+0.5*WbS24-ThS24wall,ThS24wall);
    D24->SetVertex(7,-0.5*WbS24+ThS24wall,ThS24wall);
    TGeoTube *E24 = new TGeoTube("ITS Sup Cooling Tube E24",0.5*InDcoolTub,
                                 0.5*OutDcoolTub,0.5*LS24-BlkNozZS24);
    TGeoArb8 *F24 = new TGeoArb8("ITS Sup Cable Tray lower Element block F24",
                                 0.5*BlkNozZS24);
    for(i=0;i<8;i++) F24->SetVertex(i,D24->GetVertices()[i*2+0],
                                      D24->GetVertices()[i*2+1]); //
    TGeoTube *G24 = new TGeoTube("ITS Sup Cooling Tube hole in block G24",
                                 0.0,0.5*BlkNozInDS24,0.5*BlkNozZS24);
    TGeoArb8 *H24 = new TGeoArb8("ITS Sup Cable Tray upper Element H24",
                                 0.5*(LS24- LS24PP));
    H24->SetVertex(0,C24->GetVertices()[0*2+0]+2.*ThS24wall,
                     C24->GetVertices()[0*2+1]-ThS24wall);
    H24->SetVertex(1,C24->GetVertices()[1*2+0]-2.*ThS24wall,
                     C24->GetVertices()[1*2+1]-ThS24wall);
    H24->SetVertex(2,D24->GetVertices()[1*2+0]-ThS24wall,
                     D24->GetVertices()[1*2+1]+ThS24wall);
    H24->SetVertex(3,D24->GetVertices()[0*2+0]+ThS24wall,
                     D24->GetVertices()[0*2+1]+ThS24wall);
    for(i=4;i<8;i++) H24->SetVertex(i,H24->GetVertices()[(i-4)*2+0],
                                      H24->GetVertices()[(i-4)*2+1]); //
    printArb8(C24);
    printArb8(D24);
    printTube(E24);
    printArb8(F24);
    printTube(G24);
    printArb8(H24);
    TGeoVolume *C24v,*D24v,*E24v,*F24v,*Ga24v,*Gw24v,*H24v;
    //
    C24v = new TGeoVolume("ITSsupCableTrayC24",C24,SUPal);
    C24v->SetVisibility(kTRUE);
    C24v->SetLineColor(6); //
    C24v->SetLineWidth(1);
    C24v->SetFillColor(C24v->GetLineColor());
    C24v->SetFillStyle(4000); // 0% transparent
    D24v = new TGeoVolume("ITSsupCableTrayLowerD24",D24,SUPair);
    D24v->SetVisibility(kTRUE);
    D24v->SetLineColor(6); //
    D24v->SetLineWidth(1);
    D24v->SetFillColor(D24v->GetLineColor());
    D24v->SetFillStyle(4000); // 0% transparent
    E24v = new TGeoVolume("ITSsupCableTrayCoolTubeE24",E24,SUPss);
    E24v->SetVisibility(kTRUE);
    E24v->SetLineColor(6); //
    E24v->SetLineWidth(1);
    E24v->SetFillColor(E24v->GetLineColor());
    E24v->SetFillStyle(4000); // 0% transparent
    F24v = new TGeoVolume("ITSsupCableTrayBlockF24",F24,SUPal);
    F24v->SetVisibility(kTRUE);
    F24v->SetLineColor(6); //
    F24v->SetLineWidth(1);
    F24v->SetFillColor(F24v->GetLineColor());
    F24v->SetFillStyle(4000); // 0% transparent
    Gw24v = new TGeoVolume("ITSsupCableTrayCoolantWaterG24",G24,SUPwater);
    Gw24v->SetVisibility(kTRUE);
    Gw24v->SetLineColor(6); //
    Gw24v->SetLineWidth(1);
    Gw24v->SetFillColor(Gw24v->GetLineColor());
    Gw24v->SetFillStyle(4000); // 0% transparent
    Ga24v = new TGeoVolume("ITSsupCableTrayCoolantAirG24",G24,SUPair);
    Ga24v->SetVisibility(kTRUE);
    Ga24v->SetLineColor(6); //
    Ga24v->SetLineWidth(1);
    Ga24v->SetFillColor(Ga24v->GetLineColor());
    Ga24v->SetFillStyle(4000); // 0% transparent
    H24v = new TGeoVolume("ITSsupCableTrayUpperC24",H24,SUPair);
    H24v->SetVisibility(kTRUE);
    H24v->SetLineColor(6); //
    H24v->SetLineWidth(1);
    H24v->SetFillColor(H24v->GetLineColor());
    H24v->SetFillStyle(4000); // 0% transparent
    //
    tran = new TGeoTranslation("",-OutDcoolTub,OutDcoolTub+ThS24wall,0.0);
    F24v->AddNode(Gw24v,1,tran);
    D24v->AddNode(E24v,1,tran);
    tran = new TGeoTranslation("",0.0,OutDcoolTub+ThS24wall,0.0);
    F24v->AddNode(Gw24v,2,tran);
    D24v->AddNode(E24v,2,tran);
    tran = new TGeoTranslation("",+OutDcoolTub,OutDcoolTub+ThS24wall,0.0);
    F24v->AddNode(Gw24v,3,tran);
    D24v->AddNode(E24v,3,tran);
    tran = new TGeoTranslation("",0.0,0.0,0.5*LS24-0.5*BlkNozZS24);
    D24v->AddNode(F24v,1,tran);
    tran = new TGeoTranslation("",0.0,0.0,-(0.5*LS24-0.5*BlkNozZS24));
    D24v->AddNode(F24v,2,tran);
    C24v->AddNode(D24v,1,0);
    C24v->AddNode(H24v,1,0);
    //==================================================================
    //
    // RB 26 side
    const Double_t Z026         = -900*kmm;//SSup_203A.jpg
    const Double_t ThssFrame26  = 5.0*kmm;
    const Double_t R0ssFrame26  = 444.5*kmm-ThssFrame26; // SSup_204A.jpg
    const Double_t R1ssFrame26  = 601.6*kmm-ThssFrame26; // SSup_208A.jpg
    const Double_t WidthFrame26 = 10.0*kmm;
    //const Double_t HightFrame26 = 10.0*kmm;
    const Double_t Phi0Frame26  = 15.2*kDegree; // SSup_602A.jpg
    const Double_t Phi1Frame26  = (90.0-7.6)*kDegree; // SSup_802A.jpg
    const Double_t ZssFrameSection26 = (415.0-10.0)*kmm;
    const Int_t    NZsections26      = 4;
    const Int_t    NPhiSections26    = 4;
    const Int_t    NFramesPhi26      = 4;
    TGeoConeSeg *A26[NZsections26+1],*M26; // Cylinderial support structure
    TGeoArb8     *B26; // Cylinderial support structure
    Char_t name[100];
    Double_t r1,r2,m;

    M26 = new TGeoConeSeg("ITS sup Cable tray support frame mother volume "
                          "M26",0.5*(4.*ZssFrameSection26+5*WidthFrame26),
                          R1ssFrame26,R1ssFrame26+ThssFrame26,
                          R0ssFrame26,R0ssFrame26+ThssFrame26,
                          Phi0Frame26,Phi1Frame26);
    m = -((R1ssFrame26-R0ssFrame26)/
         (((Double_t)NZsections26)*(ZssFrameSection26+WidthFrame26)));
    for(i=0;i<NZsections26+1;i++){
        di = ((Double_t) i)*(ZssFrameSection26+WidthFrame26);
        sprintf(name,
                "ITS sup Cable tray support frame radial section A26[%d]",i);
        r1 = R1ssFrame26+m*di;
        r2 = R1ssFrame26+m*(di+WidthFrame26);
        A26[i] = new TGeoConeSeg(name,0.5*WidthFrame26,r2,r2+ThssFrame26,
                                 r1,r1+ThssFrame26,Phi0Frame26,Phi1Frame26);
    } // end for i
    B26 = new TGeoArb8("ITS sup Cable tray support frame Z section B26",
                       0.5*ZssFrameSection26);
    r = 0.25*(A26[0]->GetRmax1()+A26[0]->GetRmin1()+
              A26[1]->GetRmax2()+A26[1]->GetRmin2());
    B26->SetVertex(0,A26[0]->GetRmax2()-r,+0.5*WidthFrame26);
    B26->SetVertex(1,A26[0]->GetRmax2()-r,-0.5*WidthFrame26);
    B26->SetVertex(2,A26[0]->GetRmin2()-r,-0.5*WidthFrame26);
    B26->SetVertex(3,A26[0]->GetRmin2()-r,+0.5*WidthFrame26);
    B26->SetVertex(4,A26[1]->GetRmax1()-r,+0.5*WidthFrame26);
    B26->SetVertex(5,A26[1]->GetRmax1()-r,-0.5*WidthFrame26);
    B26->SetVertex(6,A26[1]->GetRmin1()-r,-0.5*WidthFrame26);
    B26->SetVertex(7,A26[1]->GetRmin1()-r,+0.5*WidthFrame26);
    for(i=0;i<NZsections26+1;i++) printConeSeg(A26[i]);
    printConeSeg(M26);
    printArb8(B26);
    TGeoVolume *A26v[NZsections26+1],*B26v,*M26v;
    //
    for(i=0;i<NZsections26+1;i++){
        sprintf(name,"ITSsupFrameA26[%d]",i);
        A26v[i] = new TGeoVolume(name,A26[i],SUPss);
        A26v[i]->SetVisibility(kTRUE);
        A26v[i]->SetLineColor(1); // black
        A26v[i]->SetLineWidth(1);
        A26v[i]->SetFillColor(A26v[i]->GetLineColor());
        A26v[i]->SetFillStyle(4000); // 0% transparent
    } // end for i
    B26v = new TGeoVolume("ITSsupFrameB26",B26,SUPss);
    B26v->SetVisibility(kTRUE);
    B26v->SetLineColor(1); // black
    B26v->SetLineWidth(1);
    B26v->SetFillColor(B26v->GetLineColor());
    B26v->SetFillStyle(4000); // 0% transparent
    M26v = new TGeoVolume("ITSsupFrameM26",M26,SUPair);
    M26v->SetVisibility(kTRUE);
    M26v->SetLineColor(7); // light blue
    M26v->SetLineWidth(1);
    M26v->SetFillColor(M26v->GetLineColor());
    M26v->SetFillStyle(4090); // 90% transparent
    //
    Int_t NcB26=1;
    t0 = Phi0Frame26;
    dt = (Phi1Frame26-Phi0Frame26)/((Double_t)NPhiSections26);
    for(i=0;i<=NZsections26;i++){
        di = ((Double_t) i)*(ZssFrameSection26+WidthFrame26);
        z = -M26->GetDz()+A26[i]->GetDz() + di;
        tran = new TGeoTranslation("",0.0,0.0,z);
        M26v->AddNode(A26v[i],1,tran);
        z = z+B26->GetDz();
        if(i<NZsections26)for(j=0;j<=NPhiSections26;j++){
            r = 0.25*(A26[i]->GetRmax1()+A26[i]->GetRmin1()+
                      A26[i+1]->GetRmax2()+A26[i+1]->GetRmin2());
            t = t0 + ((Double_t)j)*dt;
            rot = new TGeoRotation("",0.0,0.0,t);
            y = r*TMath::Sin(t*kRadian);
            x = r*TMath::Cos(t*kRadian);
            tranrot = new TGeoCombiTrans("",x,y,z,rot);
            delete rot; // rot not explicity used in AddNode functions.
            M26v->AddNode(B26v,NcB26++,tranrot);
        } // end for j
    } // end for i
    tran = new TGeoTranslation("",0.0,0.0,Z026-M26->GetDz());
    Moth->AddNode(M26v,1,tran);
    for(i=1;i<NFramesPhi26;i++){
        rot = new TGeoRotation("",0.0,0.0,90.0*((Double_t)i));
        tranrot = new TGeoCombiTrans(*tran,*rot);
        delete rot; // rot not explicity used in AddNode functions.
        Moth->AddNode(M26v,i+1,tranrot);
    } // end for i
    if(GetDebug()){
        for(i=0;i<NZsections26+1;i++) A26v[i]->PrintNodes();
        B26v->PrintNodes();
        M26v->PrintNodes();
    } // end if
}
