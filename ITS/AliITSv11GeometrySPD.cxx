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
#include <Riostream.h>
#include <TMath.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TPolyLine.h>
// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoEltu.h>
#include <TGeoXtru.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
//#include <TGeoRotation.h>
//#include <TGeoCombiTrans.h>
//#include <TGeoTranslation.h>
#include "AliITSv11GeometrySPD.h"

ClassImp(AliITSv11GeometrySPD)

#define SQ(A) (A)*(A)

//______________________________________________________________________
void AliITSv11GeometrySPD::CarbonFiberSector(TGeoVolume *moth){
    // Define the detail SPD Carbon fiber support Sector geometry.
    // Based on the drawings ALICE-Pixel "Construzione Profilo Modulo"
    // March 25 2004 and ALICE-SUPPORTO "construzione Profilo Modulo"
    // Define Outside radii as negitive, Outside in the sence that the
    // center of the arc is outside of the object.
    // February 16 2004.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.
    TGeoManager *mgr = gGeoManager;
    TGeoMedium *medSPDcf  = 0; // SPD support cone Carbon Fiber materal number.
    //TGeoMedium *medSPDfs  = 0; // SPD support cone inserto stesalite 4411w.
    //TGeoMedium *medSPDfo  = 0; // SPD support cone foam, Rohacell 50A.
    TGeoMedium *medSPDss  = 0; // SPD support cone screw material,Stainless
    TGeoMedium *medSPDair = 0; // SPD support cone Air
    //TGeoMedium *medSPDal  = 0; // SPD support cone SDD mounting bracket Al
    TGeoMedium *medSPDcoolfl  = 0; // SPD cooling fluid, Freeon
    medSPDcf = mgr->GetMedium("ITSspdCarbonFiber");
    //medSPDfs = mgr->GetMedium("ITSspdStaselite4411w");
    //medSPDfo = mgr->GetMedium("ITSspdRohacell50A");
    medSPDss = mgr->GetMedium("ITSspdStainlessSteal");
    medSPDair= mgr->GetMedium("ITSspdAir");
    medSPDcoolfl= mgr->GetMedium("ITSspdCoolingFluid");
    //
    const Double_t ksecDz        = 0.5*500.0*fgkmm;
    const Double_t ksecLen       = 30.0*fgkmm;
    const Double_t ksecCthick    = 0.20*fgkmm;
    const Double_t ksecDipLength = 3.2*fgkmm;
    const Double_t ksecDipRadii  = 0.4*fgkmm;
    //const Double_t ksecCoolingTubeExtraDepth = 0.86*fgkmm;
    //
    const Double_t ksecX0   = -10.725*fgkmm;
    const Double_t ksecY0   = -14.853*fgkmm;
    const Double_t ksecR0   = -0.8*fgkmm; // Outside
    const Double_t ksecX1   = -13.187*fgkmm;
    const Double_t ksecY1   = -19.964*fgkmm;
    const Double_t ksecR1   = +0.6*fgkmm; // Inside
    const Double_t ksecDip0 = 5.9*fgkmm;
    //
    const Double_t ksecX2   = -3.883*fgkmm;
    const Double_t ksecY2   = -17.805*fgkmm;
    const Double_t ksecR2   = +0.80*fgkmm; // Inside Guess. 
    const Double_t ksecX3   = -3.123*fgkmm;
    const Double_t ksecY3   = -14.618*fgkmm;
    const Double_t ksecR3   = -0.6*fgkmm; // Outside
    const Double_t ksecDip1 = 8.035*fgkmm;
    //
    const Double_t ksecX4   = +11.280*fgkmm;
    const Double_t ksecY4   = -14.473*fgkmm;
    const Double_t ksecR4   = +0.8*fgkmm; // Inside
    const Double_t ksecX5   = +19.544*fgkmm;
    const Double_t ksecY5   = +10.961*fgkmm;
    const Double_t ksecR5   = +0.8*fgkmm; // Inside
    const Double_t ksecDip2 = 4.553*fgkmm;
    //
    const Double_t ksecX6   = +10.830*fgkmm;
    const Double_t ksecY6   = +16.858*fgkmm;
    const Double_t ksecR6   = +0.6*fgkmm; // Inside
    const Double_t ksecX7   = +11.581*fgkmm;
    const Double_t ksecY7   = +13.317*fgkmm;
    const Double_t ksecR7   = -0.6*fgkmm; // Outside
    const Double_t ksecDip3 = 6.978*fgkmm;
    //
    const Double_t ksecX8   = -0.733*fgkmm;
    const Double_t ksecY8   = +17.486*fgkmm;
    const Double_t ksecR8   = +0.6*fgkmm; // Inside
    const Double_t ksecX9   = +0.562*fgkmm;
    const Double_t ksecY9   = +14.486*fgkmm;
    const Double_t ksecR9   = -0.6*fgkmm; // Outside
    const Double_t ksecDip4 = 6.978*fgkmm;
    //
    const Double_t ksecX10  = -12.252*fgkmm;
    const Double_t ksecY10  = +16.298*fgkmm;
    const Double_t ksecR10  = +0.6*fgkmm; // Inside
    const Double_t ksecX11  = -10.445*fgkmm;
    const Double_t ksecY11  = +13.162*fgkmm;
    const Double_t ksecR11  = -0.6*fgkmm; // Outside
    const Double_t ksecDip5 = 6.978*fgkmm;
    //
    const Double_t ksecX12  = -22.276*fgkmm;
    const Double_t ksecY12  = +12.948*fgkmm;
    const Double_t ksecR12  = +0.85*fgkmm; // Inside
    //const Double_t ksecX13 = *fgkmm;
    //const Double_t ksecY13 = *fgkmm;
    const Double_t ksecR13  = -0.8*fgkmm; // Outside
    const Double_t ksecAngleSide13 = 36.0*fgkDegree;
    //
    const Int_t ksecNRadii = 14;
    const Int_t ksecNPointsPerRadii = 4;
    const Int_t ksecNCoolingTubeDips = 6;
    // Since the Rounded parts are aproximated by a regular polygon and
    // a cooling tube of the propper diameter must fit, a scaling factor
    // increases the size of the polygon for the tube to fit.
    //const Double_t ksecRCoolScale = 1./TMath::Cos(TMath::Pi()/
    //                                          (Double_t)ksecNPointsPerRadii);
    const Double_t ksecZEndLen  = 30.00*fgkmm;
    //const Double_t ksecZFlangLen= 45.00*fgkmm;
    const Double_t ksecTl       = 0.860*fgkmm;
    const Double_t ksecCthick2  = 0.600*fgkmm;
    const Double_t ksecCthick3  = 1.800*fgkmm;
    const Double_t ksecSidelen  = 22.00*fgkmm;
    const Double_t ksecSideD5   = 3.679*fgkmm;
    const Double_t ksecSideD12  = 7.066*fgkmm;
    const Double_t ksecRCoolOut = 2.400*fgkmm;
    const Double_t ksecRCoolIn  = 2.000*fgkmm;
    const Double_t ksecDl1      = 5.900*fgkmm;
    const Double_t ksecDl2      = 8.035*fgkmm;
    const Double_t ksecDl3      = 4.553*fgkmm;
    const Double_t ksecDl4      = 6.978*fgkmm;
    const Double_t ksecDl5      = 6.978*fgkmm;
    const Double_t ksecDl6      = 6.978*fgkmm;
    const Double_t ksecCoolTubeThick = 10.0*fgkmicron;
    //
    const Int_t ksecNPoints = (ksecNPointsPerRadii+1)*(ksecNRadii+
                               ksecNCoolingTubeDips) + 8;
    Double_t secX[ksecNRadii] = {ksecX0,ksecX1,ksecX2,ksecX3,ksecX4,ksecX5,
                                 ksecX6,ksecX7,ksecX8,ksecX9,ksecX10,ksecX11,
                                 ksecX12,-1000.0};
    Double_t secY[ksecNRadii] = {ksecY0,ksecY1,ksecY2,ksecY3,ksecY4,ksecY5,
                                 ksecY6,ksecY7,ksecY8,ksecY9,ksecY10,ksecY11,
                                 ksecY12,-1000.0};
    Double_t secR[ksecNRadii] = {ksecR0,ksecR1,ksecR2,ksecR3,ksecR4,ksecR5,
                                 ksecR6,ksecR7,ksecR8,ksecR9,ksecR10,ksecR11,
                                 ksecR12,ksecR13};
    Double_t secDip[ksecNRadii]={0.0,ksecDip0,0.0,ksecDip1,0.0,ksecDip2,0.0,
                                 ksecDip3,0.0,ksecDip4,0.0,ksecDip5,0.0,0.0};
    Double_t secX2[ksecNRadii+ksecNCoolingTubeDips] = {
        ksecX0,ksecX1,-1000.,ksecX2,ksecX3,-1000.,ksecX4,ksecX5,
        -1000.,ksecX6,ksecX7,-1000.,ksecX8,ksecX9,-1000.,
        ksecX10,ksecX11,-1000.,ksecX12,-1000.0};
    Double_t secY2[ksecNRadii+ksecNCoolingTubeDips] = {
        ksecY0,ksecY1,-1000.,ksecY2,ksecY3,-1000.,ksecY4,ksecY5,
        -1000.,ksecY6,ksecY7,-1000.,ksecY8,ksecY9,-1000.,
        ksecY10,ksecY11,-1000.,ksecY12,-1000.0};
    Double_t secR2[ksecNRadii+ksecNCoolingTubeDips] = {
        ksecR0,ksecR1,ksecRCoolOut,ksecR2,ksecR3,ksecRCoolOut,ksecR4,ksecR5,
        ksecRCoolOut,ksecR6,ksecR7,ksecRCoolOut,ksecR8,ksecR9,ksecRCoolOut,
        ksecR10,ksecR11,ksecRCoolOut,ksecR12,ksecR13};
    Double_t secDip2[ksecNRadii]={0.0,ksecDl1,0.0,ksecDl2,0.0,ksecDl3,0.0,
                                 ksecDl4,0.0,ksecDl5,0.0,ksecDl6,0.0,0.0};
    const Int_t ksecDipIndex[ksecNCoolingTubeDips] = {2,5,8,11,14,17};
    Double_t secAngleStart[ksecNRadii];
    Double_t secAngleEnd[ksecNRadii];
    Double_t secAngleStart2[ksecNRadii+ksecNCoolingTubeDips];
    Double_t secAngleEnd2[ksecNRadii+ksecNCoolingTubeDips];
    Double_t secAngleStart3[ksecNRadii+ksecNCoolingTubeDips];
    Double_t secAngleEnd3[ksecNRadii+ksecNCoolingTubeDips];
    Double_t xp[ksecNPoints],yp[ksecNPoints];
    TGeoXtru *sA0,*sA1,*sB0,*sB1;
    TGeoEltu *sTA0,*sTA1;
    TGeoTube *sTB0,*sTB1;
    TGeoRotation    *rot;
    TGeoTranslation *trans;
    TGeoCombiTrans  *rotrans;
    Double_t t,t0,t1,a,b,x0,y0,x1,y1;
    Int_t i,j,k;

    if(moth==0){
        Error("CarbonFiberSector","moth=%p",moth);
        return;
    } // end if moth==0
    SetDebug(3);

    secAngleStart[0] = 0.5*ksecAngleSide13;
    for(i=0;i<ksecNRadii-2;i++){
        AnglesForRoundedCorners(secX[i],secY[i],secR[i],
                                secX[i+1],secY[i+1],secR[i+1],t0,t1);
        secAngleEnd[i]     = t0;
        if(secR[i]>0.0&&secR[i+1]>0.0)if(secAngleStart[i]>secAngleEnd[i])
            secAngleEnd[i] += 360.0;
        secAngleStart[i+1] = t1;
    } // end for i
    secAngleEnd[ksecNRadii-2]   = secAngleStart[ksecNRadii-2] + 
                                     (secAngleEnd[10]-secAngleStart[10]);
    if(secAngleEnd[ksecNRadii-2]<0.0) secAngleEnd[ksecNRadii-2] += 360.0;
    secAngleStart[ksecNRadii-1] = secAngleEnd[ksecNRadii-2] - 180.0;
    secAngleEnd[ksecNRadii-1]   = secAngleStart[0];
    //
    i = 0;
    j = ksecNRadii-2;
    t0 = TanD(secAngleStart[i]-90.);
    t1 = TanD(secAngleEnd[j]-90.);
    t  = secY[i] - secY[j];
    // Note, secR[i=0] <0; secR[j=12]>0; and secR[j+1=13] <0
    t += (-secR[i]+secR[j+1])*SinD(secAngleStart[i]);
    t -= (secR[j]-secR[j+1])*SinD(secAngleEnd[j]);
    t += t1*secX[j] - t0*secX[i];
    t += t1*(secR[j]-secR[j+1])*CosD(secAngleEnd[j]);
    t -= t0*(-secR[i]+secR[j+1])*CosD(secAngleStart[i]);
    secX[ksecNRadii-1] = t/(t1-t0);
    secY[ksecNRadii-1] = TanD(90.+0.5*ksecAngleSide13)*
                          (secX[ksecNRadii-1]-secX[0]) + secY[0];
    //
    if(GetDebug(2)){
        cout <<"    X    \t  Y  \t  R  \t  S  \t  E"<<endl;
        for(i=0;i<ksecNRadii;i++){
            cout <<"{"<< secX[i] <<",";
            cout << secY[i] <<",";
            cout << secR[i] <<",";
            cout << secAngleStart[i] <<",";
            cout << secAngleEnd[i] <<"},"<< endl;
        } // end for i
    } // end if GetDebug
    //
    if(GetDebug(3)) cout <<"Double_t sA0[][";
    if(GetDebug(4)) cout <<"3]{";
    else if(GetDebug(3)) cout <<"2]{";
    j = -1;
    t0 = (Double_t)ksecNPointsPerRadii;
    for(i=0;i<ksecNRadii;i++){
        t1 = (secAngleEnd[i]-secAngleStart[i])/t0;
        if(GetDebug(5)) cout<<"t1="<< t1<<endl;
        for(k=0;k<=ksecNPointsPerRadii;k++){
            t=secAngleStart[i]+((Double_t)k)*t1;
            j++;
            xp[j] = TMath::Abs(secR[i])*CosD(t)+secX[i];
            yp[j] = TMath::Abs(secR[i])*SinD(t)+secY[i];
            if(GetDebug(3)){
                cout << "{"<<xp[j]<<","<<yp[j];
                if(GetDebug(4)) cout <<","<<t;
                cout <<"},";
            } // end if GetDebug
        } // end for k
        if(GetDebug(3)) cout << endl;
        t = secAngleEnd[i];
        a = ksecDipLength+2.0*(ksecDipRadii);
        b = secDip[i]-0.5*a;
        switch (i){
        case 1: case 5: // Dip0,2
            j++;
            xp[j] = xp[j-1]-b*CosD(t-90.);
            yp[j] = yp[j-1]-b*SinD(t-90.);
            j++;
            xp[j] = xp[j-1]-(2.0*ksecDipRadii)*CosD(t);
            yp[j] = yp[j-1]-(2.0*ksecDipRadii)*SinD(t);
            j++;
            xp[j] = xp[j-1]-a*CosD(t-90.);
            yp[j] = yp[j-1]-a*SinD(t-90.);
            j++;
            xp[j] = xp[j-1]+(2.0*ksecDipRadii)*CosD(t);
            yp[j] = yp[j-1]+(2.0*ksecDipRadii)*SinD(t);
            if(GetDebug(3))for(k=-3;k<=0;k++){
                cout << "{"<<xp[j+k]<<","<<yp[j+k];
                if(GetDebug(4)) cout <<","<<0.0;
                cout <<"},";
            } // end if GetDebug
            if(GetDebug(3)) cout << endl;
            break;
        case 3: case 7: case 9: case 11:// Dip 1,3,4,5
            j++;
            xp[j] = xp[j-1]+b*CosD(t-90.);
            yp[j] = yp[j-1]+b*SinD(t-90.);
            j++;
            xp[j] = xp[j-1]+(2.0*ksecDipRadii)*CosD(t);
            yp[j] = yp[j-1]+(2.0*ksecDipRadii)*SinD(t);
            j++;
            xp[j] = xp[j-1]+a*CosD(t-90.);
            yp[j] = yp[j-1]+a*SinD(t-90.);
            j++;
            xp[j] = xp[j-1]-(2.0*ksecDipRadii)*CosD(t);
            yp[j] = yp[j-1]-(2.0*ksecDipRadii)*SinD(t);
            if(GetDebug(3))for(k=-3;k<=0;k++){
                cout << "{"<<xp[j+k]<<","<<yp[j+k];
                if(GetDebug(4)) cout <<","<<0.0;
                cout <<"},";
            } // end if GetDebug
            if(GetDebug(3)) cout << endl;
            break;
        default:
            break;
        } // end switch
    } // end of i
    if(GetDebug(3)) cout<<"{"<<xp[0]<<","<<yp[0];
    if(GetDebug(4)) cout<<","<< secAngleStart[0];
    if(GetDebug(3)) cout<<"}} j="<<j<<endl;
    sA0 = new TGeoXtru(2);
    sA0->SetName("ITS SPD Carbon fiber support Sector A0");
    sA0->DefinePolygon(j+1,xp,yp);
    sA0->DefineSection(0,-ksecDz);
    sA0->DefineSection(1,ksecDz);
    //
    if(GetDebug(3)) cout <<"Double_t sA1[][{";
    if(GetDebug(4)) cout <<"3]{";
    else if(GetDebug(3)) cout <<"2]{";
    j = -1;
    t0 = (Double_t)ksecNPointsPerRadii;
    for(i=0;i<ksecNRadii;i++){
        t1 = (secAngleEnd[i]-secAngleStart[i])/t0;
        if(GetDebug(5)) cout<<"t1="<< t1<<endl;
        for(k=0;k<=ksecNPointsPerRadii;k++){
            t=secAngleStart[i]+((Double_t)k)*t1;
            j++;
            xp[j] = TMath::Abs(secR[i]-ksecCthick)*CosD(t)+secX[i];
            yp[j] = TMath::Abs(secR[i]-ksecCthick)*SinD(t)+secY[i];
            if(GetDebug(3)){
                cout << "{"<<xp[j]<<","<<yp[j];
                if(GetDebug(4)) cout <<","<<t;
                cout <<"},";
            } // end if GetDebug
        } // end for t
        if(GetDebug(3)) cout << endl;
        t = secAngleEnd[i];
        a = ksecDipLength+2.0*(ksecDipRadii+ksecCthick);
        b = secDip[i]-0.5*a;
        switch (i){
        case 1: case 5: // Dip0,2
            j++;
            xp[j] = xp[j-1]-b*CosD(t-90.);
            yp[j] = yp[j-1]-b*SinD(t-90.);
            j++;
            xp[j] = xp[j-1]-(2.0*ksecDipRadii)*CosD(t);
            yp[j] = yp[j-1]-(2.0*ksecDipRadii)*SinD(t);
            j++;
            xp[j] = xp[j-1]-a*CosD(t-90.);
            yp[j] = yp[j-1]-a*SinD(t-90.);
            j++;
            xp[j] = xp[j-1]+(2.0*ksecDipRadii)*CosD(t);
            yp[j] = yp[j-1]+(2.0*ksecDipRadii)*SinD(t);
            if(GetDebug(3))for(k=-3;k<=0;k++){
                cout << "{"<<xp[j+k]<<","<<yp[j+k];
                if(GetDebug(4)) cout <<",t="<<0.0;
                cout <<"},";
            } // end if GetDebug
            if(GetDebug(3)) cout << endl;
            break;
        case 3: case 7: case 9: case 11:// Dip 1,3,4,5
            j++;
            xp[j] = xp[j-1]+b*CosD(t-90.);
            yp[j] = yp[j-1]+b*SinD(t-90.);
            j++;
            xp[j] = xp[j-1]+(2.0*ksecDipRadii)*CosD(t);
            yp[j] = yp[j-1]+(2.0*ksecDipRadii)*SinD(t);
            j++;
            xp[j] = xp[j-1]+a*CosD(t-90.);
            yp[j] = yp[j-1]+a*SinD(t-90.);
            j++;
            xp[j] = xp[j-1]-(2.0*ksecDipRadii)*CosD(t);
            yp[j] = yp[j-1]-(2.0*ksecDipRadii)*SinD(t);
            if(GetDebug(3))for(k=-3;k<=0;k++){
                cout << "{"<<xp[j+k]<<","<<yp[j+k];
                if(GetDebug(4)) cout <<",t="<<0.0;
                cout <<"},";
            } // end if GetDebug
            if(GetDebug(3)) cout << endl;
            break;
        default:
            break;
        } // end switch
    } // end of i
    if(GetDebug(3)) cout<<"{"<<xp[0]<<","<<yp[0];
    if(GetDebug(4)) cout<<","<< secAngleStart[0];
    if(GetDebug(3)) cout<<"}} j="<<j<<endl;
    sA1 = new TGeoXtru(2);
    sA1->SetName("ITS SPD Carbon fiber support Sector Air A1");
    sA1->DefinePolygon(j+1,xp,yp);
    sA1->DefineSection(0,-ksecDz);
    sA1->DefineSection(1,ksecDz);
    //
    sTA0 = new TGeoEltu("ITS SPD Cooling Tube TA0",
                       ksecDipRadii,0.5*ksecDipLength,ksecDz);
    sTA1 = new TGeoEltu("ITS SPD Cooling Tube coolant TA1",
                        sTA0->GetA()-ksecCoolTubeThick,
                        sTA0->GetB()-ksecCoolTubeThick,ksecDz);
    //
    j = 0;
    for(i=0;i<ksecNRadii;i++){
        secAngleStart2[j] = secAngleStart[i];
        secAngleEnd2[j]   = secAngleEnd[i];
        secAngleStart3[j] = secAngleStart2[j];
        secAngleEnd3[j]   = secAngleEnd2[j];
        secX2[j]          = secX[i];
        secY2[j]          = secY[i];
        secR2[j]          = secR[i];
        j++;
        t = secAngleEnd[i];
        switch (i){
        case 1: case 5: // Tube 0,2
            x0 = secX[i] + TMath::Abs(secR[i])*CosD(t); // last point of turn
            y0 = secY[i] + TMath::Abs(secR[i])*SinD(t);
            x1 = x0-secDip2[i]*CosD(t-90.);  // center point of dip 
            y1 = y0-secDip2[i]*SinD(t-90.);  // along line
            secX2[j] = x1-ksecTl*CosD(t);  // location of circle center.
            secY2[j] = y1-ksecTl*SinD(t);
            x1 = secX[i+1]+TMath::Abs(secR[i+1])*CosD(secAngleStart2[i+1]);
            y1 = secY[i+1]+TMath::Abs(secR[i+1])*SinD(secAngleStart2[i+1]);
            //Find starting and ending angles, break if error.
            if(!AngleOfIntersectionWithLine(x0,y0,x1,y1,secX2[j],secY2[j],
            secR2[j],secAngleStart2[j],secAngleEnd2[j]))break;
            // thicknes of Carbon fiber over the cooling tubes.
            a  = ksecRCoolOut-ksecRCoolIn;
            // last point of turn
            x0 = secX[i]+(TMath::Abs(secR[i])-ksecCthick2)*CosD(t);
            y0 = secY[i]+(TMath::Abs(secR[i])-ksecCthick2)*SinD(t);
            x1 = secX[i+1]+(TMath::Abs(secR[i+1])-ksecCthick2)*CosD(secAngleStart3[i+1]);
            y1 = secY[i+1]+(TMath::Abs(secR[i+1])-ksecCthick2)*SinD(secAngleStart3[i+1]);
            //Find starting and ending angles, break if error.
            if(!AngleOfIntersectionWithLine(x0,y0,x1,y1,secX2[j],secY2[j],
            secR2[j]-a,secAngleStart3[j],secAngleEnd3[j]))break;
            if(i==1) { // Fix odd case.
                x0 = secAngleStart3[j];
                secAngleStart3[j] = secAngleEnd3[j] - 360.0;
                secAngleEnd3[j] = x0;
            } // end if i==1
            if(i==5) { // Fix odd case.
                x0 = secAngleStart2[j];
                secAngleStart2[j] = secAngleEnd2[j] - 360.0;
                secAngleEnd2[j] = x0;
                x0 = secAngleStart3[j];
                secAngleStart3[j] = secAngleEnd3[j] - 360.0;
                secAngleEnd3[j] = x0;
            } // end if i==5
            // Because a polygon is replacing the rounded surface, the
            // radius of the polygon must be larger to make room for the
            // cooling tube of the same size. The radio of the radii of
            // a circle fitting the inside/outside of a regualr polygon
            // is given by Cos(180/n) where n is the number of sides of the
            // regurla polygon. In this case, it is scalled for the partical
            // circles involved.
            secR2[j] = secR2[j]/CosD(0.5*(secAngleEnd[j]-secAngleStart[j])/
                                     ((Double_t)ksecNPointsPerRadii));
            j++;
            break;
        case 3: case 7: case 9: case 11:// Tube 1,2,4,5,6
            x0 = secX[i] + TMath::Abs(secR[i])*CosD(t); // last point of turn
            y0 = secY[i] + TMath::Abs(secR[i])*SinD(t);
            x1 = x0+secDip2[i]*CosD(t-90.);  // center point of dip 
            y1 = y0+secDip2[i]*SinD(t-90.);  // along line
            secX2[j] = x1+ksecTl*CosD(t);  // location of circle center.
            secY2[j] = y1+ksecTl*SinD(t);
            x1 = secX[i+1]+TMath::Abs(secR[i+1])*CosD(secAngleStart2[i+1]);
            y1 = secY[i+1]+TMath::Abs(secR[i+1])*SinD(secAngleStart2[i+1]);
            if(!AngleOfIntersectionWithLine(x0,y0,x1,y1,secX2[j],secY2[j],
            secR2[j],secAngleStart2[j],secAngleEnd2[j]))break;//don't intersect
            // thicknes of Carbon fiber over the cooling tubes.
            a  = ksecRCoolOut-ksecRCoolIn;
            // last point of turn
            x0 = secX[i] + (TMath::Abs(secR[i])+ksecCthick2)*CosD(t);
            y0 = secY[i] + (TMath::Abs(secR[i])+ksecCthick2)*SinD(t);
            x1 = secX[i+1]+(TMath::Abs(secR[i+1])-ksecCthick2)*CosD(secAngleStart3[i+1]);
            y1 = secY[i+1]+(TMath::Abs(secR[i+1])-ksecCthick2)*SinD(secAngleStart3[i+1]);
            //Find starting and ending angles, break if error.
            if(!AngleOfIntersectionWithLine(x0,y0,x1,y1,secX2[j],secY2[j],
            secR2[j]-a,secAngleStart3[j],secAngleEnd3[j]))break;
            if(i==7) { // Fix odd case
                x0 = secAngleStart2[j];
                secAngleStart2[j] = secAngleEnd2[j] - 360.0;
                secAngleEnd2[j] = x0;
                x0 = secAngleStart3[j];
                secAngleStart3[j] = secAngleEnd3[j] - 360.0;
                secAngleEnd3[j] = x0;
            } // end if i==7
            if(i==9) { // Fix odd case
                x0 = secAngleStart3[j];
                secAngleStart3[j] = secAngleEnd3[j] - 360.0;
                secAngleEnd3[j] = x0;
            } // end if i==7
            // Because a polygon is replacing the rounded surface, the
            // radius of the polygon must be larger to make room for the
            // cooling tube of the same size. The radio of the radii of
            // a circle fitting the inside/outside of a regualr polygon
            // is given by Cos(180/n) where n is the number of sides of the
            // regurla polygon. In this case, it is scalled for the partical
            // circles involved.
            secR2[j] = secR2[j]/CosD(0.5*(secAngleEnd[j]-secAngleStart[j])/
                                     ((Double_t)ksecNPointsPerRadii));
            j++;
            break;
        }// end switch
    } // end for i
    //
    if(GetDebug(2)){
        cout <<"    X2   \t  Y2 \t  R2 \t  S2 \t  E2"<<endl;
        for(i=0;i<j;i++){
            cout <<"{"<< secX2[i] <<",";
            cout << secY2[i] <<",";
            cout << secR2[i] <<",";
            cout << secAngleStart2[i] <<",";
            cout << secAngleEnd2[i] <<"}," <<  endl;
        } // end for i
    } // end if GetDebug
    if(GetDebug(2)){
        cout <<"    X2   \t  Y2 \t  R2 \t  S3 \t  E3"<<endl;
        for(i=0;i<j;i++){
            cout <<"{"<< secX2[i] <<",";
            cout << secY2[i] <<",";
            cout << secR2[i] <<",";
            cout << secAngleStart3[i] <<",";
            cout << secAngleEnd3[i] <<"}," <<  endl;
        } // end for i
    } // end if GetDebug
    if(GetDebug(3)) cout <<"Double_t sB0[][";
    if(GetDebug(4)) cout <<"3]{";
    else if(GetDebug(3)) cout <<"2]{";
    j = -1;
    t0 = (Double_t)ksecNPointsPerRadii;
    for(i=0;i<ksecNRadii+ksecNCoolingTubeDips;i++){
        t1 = (secAngleEnd2[i]-secAngleStart2[i])/t0;
        if(GetDebug(5)) cout<<"t1="<< t1<<endl;
        for(k=0;k<=ksecNPointsPerRadii;k++){
            t=secAngleStart2[i]+((Double_t)k)*t1;
            j++;
            xp[j] = TMath::Abs(secR2[i])*CosD(t)+secX2[i];
            yp[j] = TMath::Abs(secR2[i])*SinD(t)+secY2[i];
            if(GetDebug(3)){
                cout << "{"<<xp[j]<<","<<yp[j];
                if(GetDebug(4)) cout <<","<<t;
                cout <<"},";
            } // end if GetDebug
        } // end for k
        if(GetDebug(3)) cout << endl;
        if(i==6) { // add thicker side
            b = CosD(0.5*ksecAngleSide13);
            a = SinD(0.5*ksecAngleSide13);
            x1 = xp[j];
            y1 = yp[j];
            x0 = a*a*ksecX5 + b*b*x1 - (y1-ksecY5)*a*b;
            y0 = a*a*y1 + b*b*ksecY5 - (x1-ksecX5)*a*b;
            j++;
            xp[j+3] = x0 - ksecSideD5*a;
            yp[j+3] = y0 - ksecSideD5*b;
            xp[j+2] = xp[j+3] + (ksecCthick3-ksecCthick2)*b;
            yp[j+2] = yp[j+3] - (ksecCthick3-ksecCthick2)*a;
            xp[j+1] = xp[j+2] - ksecSidelen*a;
            yp[j+1] = yp[j+2] - ksecSidelen*b;
            xp[j]   = xp[j+1] - (ksecCthick3-ksecCthick2)*b;
            yp[j]   = yp[j+1] + (ksecCthick3-ksecCthick2)*a;
            j += 3;
            if(GetDebug(3))for(k=-3;k<=0;k++){
                cout << "{"<<xp[j+k]<<","<<yp[j+k];
                if(GetDebug(4)) cout <<",t="<<0.0;
                cout <<"},";
            } // end if GetDebug
            if(GetDebug(3)) cout << endl;
        } // end if i==6
        if(i==19) { // add thicker side
            // first propogate referece point 12 to -18 degree edge
            b = CosD(0.5*ksecAngleSide13);
            a = SinD(0.5*ksecAngleSide13);
            x1 = secX[0]+TMath::Abs(secR[0])*CosD(secAngleStart[0]);
            y1 = secY[0]+TMath::Abs(secR[0])*SinD(secAngleStart[0]);
            x0 = a*b*(y1-secY[ksecNRadii-2]) +
                     a*a*secX[ksecNRadii-2] + b*b*x1;
            y0 = a*b*(x1-secX[ksecNRadii-2]) +
                b*b*secY[ksecNRadii-2] + a*a*y1;
            j++;
            xp[j] = x0 + ksecSideD12*a;
            yp[j] = y0 - ksecSideD12*b;
            j++;
            xp[j] = xp[j-1] - (ksecCthick3-ksecCthick2)*b;
            yp[j] = yp[j-1] - (ksecCthick3-ksecCthick2)*a;
            j++;
            xp[j] = xp[j-1] + ksecSidelen*a;
            yp[j] = yp[j-1] - ksecSidelen*b;
            j++;
            xp[j] = xp[j-1] + (ksecCthick3-ksecCthick2)*b;
            yp[j] = yp[j-1] + (ksecCthick3-ksecCthick2)*a;
            if(GetDebug(3))for(k=-3;k<=0;k++){
                cout << "{"<<xp[j+k]<<","<<yp[j+k];
                if(GetDebug(4)) cout <<",t="<<0.0;
                cout <<"},";
            } // end if GetDebug
            if(GetDebug(3)) cout << endl;
        } // end if i==19
    } // end for i
    if(GetDebug(3)) cout<<"{"<<xp[0]<<","<<yp[0];
    if(GetDebug(4)) cout<<","<< secAngleStart2[0];
    if(GetDebug(3)) cout<<"}} j="<<j<<endl;
    sB0 = new TGeoXtru(2);
    sB0->SetName("ITS SPD Carbon fiber support Sector End B0");
    sB0->DefinePolygon(j+1,xp,yp);
    sB0->DefineSection(0,ksecDz);
    sB0->DefineSection(1,ksecDz+ksecZEndLen);
    //
    if(GetDebug(3)) cout <<"Double_t sB1[][{";
    if(GetDebug(4)) cout <<"3]{";
    else if(GetDebug(3)) cout <<"2]{";
    j = -1;
    t0 = (Double_t)ksecNPointsPerRadii;
    for(i=0;i<ksecNRadii+ksecNCoolingTubeDips;i++){
        t1 = (secAngleEnd3[i]-secAngleStart3[i])/t0;
        if(GetDebug(5)) cout<<"t1="<< t1<<endl;
        for(k=0;k<=ksecNPointsPerRadii;k++){
            t=secAngleStart3[i]+((Double_t)k)*t1;
            j++;
            x0 = TMath::Abs(secR2[i]-ksecCthick2);
            if(i==2||i==5||i==8||i==11||i==14||i==17){
                x0 = TMath::Abs(secR2[i]-ksecRCoolOut+ksecRCoolIn);/*
                if(k==0){// compute change in start and end angles to
                    // compensate for thickness of carbon fiber
                    y0 = (ksecCthick2-ksecRCoolOut*SinD(t))/
                        ksecRCoolIn; // sin th'
                    if(GetDebug(5))cout <<" k=0 t="<<t<<" y0="<<y0<<endl;
                    y1 = TMath::Sqrt(1-y0*y0);     // cos th'
                    t -= 180.*TMath::ASin(SinD(t)*y1-CosD(t)*y0)/TMath::Pi(); 
                }else if(k==ksecNPointsPerRadii) {
                    y0 = (ksecCthick2-ksecRCoolOut*SinD(t))/
                        ksecRCoolIn; // sin th'
                    if(GetDebug(5))cout <<" k="<<k<<" t="<<t<<" y0="<<y0<<endl;
                    y1 = TMath::Sqrt(1-y0*y0);     // cos th'
                    t += 180.*TMath::ASin(SinD(t)*y1-CosD(t)*y0)/TMath::Pi();
                } // end if
            */} // end if
            xp[j] = x0*CosD(t)+secX2[i];
            yp[j] = x0*SinD(t)+secY2[i];
            if(GetDebug(3)){
                cout << "{"<<xp[j]<<","<<yp[j];
                if(GetDebug(4)) cout <<","<<t;
                cout <<"},";
            } // end if GetDebug
        } // end for k
        if(GetDebug(3)) cout << endl;
    } // end for i
    if(GetDebug(3)) cout<<"{"<<xp[0]<<","<<yp[0];
    if(GetDebug(4)) cout<<","<< secAngleStart2[0];
    if(GetDebug(3)) cout<<"}} j="<<j<<endl;
    sB1 = new TGeoXtru(2);
    sB1->SetName("ITS SPD Carbon fiber support Sector Air End B1");
    sB1->DefinePolygon(j+1,xp,yp);
    sB1->DefineSection(0,ksecDz);
    sB1->DefineSection(1,ksecDz+ksecLen);
    sTB0 = new TGeoTube("ITS SPD Cooling Tube End TB0",0.0,
                       ksecRCoolIn,0.5*ksecLen);
    sTB1 = new TGeoTube("ITS SPD Cooling Tube End coolant TB0",0.0,
                       sTB0->GetRmax()-ksecCoolTubeThick,0.5*ksecLen);
    //
    if(GetDebug()){
        sA0->InspectShape();
        sA1->InspectShape();
        sB0->InspectShape();
        sB1->InspectShape();
    } // end if GetDebug
    //
    TGeoVolume *vM,*vA0,*vA1,*vTA0,*vTA1,*vB0,*vB1,*vTB0,*vTB1;
    vM = moth;
    vA0 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorA0",sA0,medSPDcf);
    vA0->SetVisibility(kTRUE);
    vA0->SetLineColor(4); // Blue
    vA0->SetLineWidth(1);
    vA0->SetFillColor(vA0->GetLineColor());
    vA0->SetFillStyle(4010); // 10% transparent
    vA1 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorAirA1",sA1,medSPDair);
    vA1->SetVisibility(kTRUE);
    vA1->SetLineColor(7); // light Blue
    vA1->SetLineWidth(1);
    vA1->SetFillColor(vA1->GetLineColor());
    vA1->SetFillStyle(4090); // 90% transparent
    vTA0 = new TGeoVolume("ITSSPDCoolingTubeTA0",sTA0,medSPDss);
    vTA0->SetVisibility(kTRUE);
    vTA0->SetLineColor(1); // Black
    vTA0->SetLineWidth(1);
    vTA0->SetFillColor(vTA0->GetLineColor());
    vTA0->SetFillStyle(4000); // 0% transparent
    vTA1 = new TGeoVolume("ITSSPDCoolingTubeFluidTA1",sTA1,medSPDcoolfl);
    vTA1->SetVisibility(kTRUE);
    vTA1->SetLineColor(6); // Purple
    vTA1->SetLineWidth(1);
    vTA1->SetFillColor(vTA1->GetLineColor());
    vTA1->SetFillStyle(4000); // 0% transparent
    vB0 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorEndB0",sB0,medSPDcf);
    vB0->SetVisibility(kTRUE);
    vB0->SetLineColor(4); // Blue
    vB0->SetLineWidth(1);
    vB0->SetFillColor(vB0->GetLineColor());
    vB0->SetFillStyle(4010); // 10% transparent
    vB1 = new TGeoVolume("ITSSPDCarbonFiberSupportSectorEndAirB1",sB1,medSPDair);
    vB1->SetVisibility(kTRUE);
    vB1->SetLineColor(7); // light Blue
    vB1->SetLineWidth(1);
    vB1->SetFillColor(vB1->GetLineColor());
    vB1->SetFillStyle(4090); // 90% transparent
    vTB0 = new TGeoVolume("ITSSPDCoolingTubeEndTB0",sTB0,medSPDss);
    vTB0->SetVisibility(kTRUE);
    vTB0->SetLineColor(1); // Black
    vTB0->SetLineWidth(1);
    vTB0->SetFillColor(vTB0->GetLineColor());
    vTB0->SetFillStyle(4000); // 0% transparent
    vTB1 = new TGeoVolume("ITSSPDCoolingTubeEndFluidTB1",sTB1,medSPDcoolfl);
    vTB1->SetVisibility(kTRUE);
    vTB1->SetLineColor(6); // Purple
    vTB1->SetLineWidth(1);
    vTB1->SetFillColor(vTB1->GetLineColor());
    vTB1->SetFillStyle(4000); // 0% transparent
    //
    vA0->AddNode(vA1,1,0); // Put air inside carbon fiber.
    vB0->AddNode(vB1,1,0); // Put air inside carbon fiber.
    vTA0->AddNode(vTA1,1,0); // Put air inside carbon fiber.
    vTB0->AddNode(vTB1,1,0); // Put air inside carbon fiber.
    for(i=0;i<ksecNCoolingTubeDips;i++){
        x0 = secX2[ksecDipIndex[i]];
        y0 = secY2[ksecDipIndex[i]];
        trans = new TGeoTranslation("",x0,y0,0.0);
        vB1->AddNode(vTB0,i+1,trans);
        rot = new TGeoRotation("",0.0,0.0,
                               TMath::RadToDeg()*TMath::ATan2(y0,x0));
        rotrans = new TGeoCombiTrans(*trans,*rot);
        vM->AddNode(vTA0,i+1,rotrans);
        delete rot;
    } // end for i
    vM->AddNode(vA0,1,0);
    vM->AddNode(vB0,1,0);
    vM->AddNode(vB0,2,new TGeoScale(1.0,1.0,-1.0)); // Reflection.
    if(GetDebug()){
        vA0->PrintNodes();
        vA1->PrintNodes();
        vB0->PrintNodes();
        vB1->PrintNodes();
    } // end if GetDebug
    //
}
//______________________________________________________________________
void AliITSv11GeometrySPD::HalfStave(TGeoVolume *moth){
    // Define the detail SPD Half Stave geometry.
    // Inputs:
    //   none.
    // Outputs:
    //  none.
    // Return:
    //  none.

    if(moth==0){
        Error("HalfStave","moth=%p",moth);
        return;
    } // end if moth==0
}
//----------------------------------------------------------------------
void AliITSv11GeometrySPD::CreateFigure0(const Char_t *filepath,
                                         const Char_t *type){
    // Creates Figure 0 for the documentation of this class. In this
    // specific case, it creates the X,Y cross section of the SPD suport
    // section, center and ends. The output is written to a standard
    // file name to the path specificed.
    // Inputs:
    //   const Char_t *filepath  Path where the figure is to be drawn
    //   const Char_t *type      The type of file, default is gif.
    // Output:
    //   none.
    // Return:
    //   none.
    TGeoXtru *sA0,*sA1,*sB0,*sB1;
    //TPolyMarker *pmA,*pmB;
    TPolyLine plA0,plA1,plB0,plB1;
    TCanvas *canvas;
    TLatex txt;
    Double_t x,y;

    sA0 = (TGeoXtru*) gGeomManager->GetVolume(
        "ITSSPDCarbonFiberSupportSectorA0_1")->GetShape();
    sA1 = (TGeoXtru*) gGeomManager->GetVolume(
        "ITSSPDCarbonFiberSupportSectorAirA1_1")->GetShape();
    sB0 = (TGeoXtru*) gGeomManager->GetVolume(
        "ITSSPDCarbonFiberSupportSectorEndB0_1")->GetShape();
    sB1 = (TGeoXtru*) gGeomManager->GetVolume(
        "ITSSPDCarbonFiberSupportSectorEndAirB1_1")->GetShape();
    //pmA = new TPolyMarker();
    //pmA.SetMarkerStyle(2); // +
    //pmA.SetMarkerColor(7); // light blue
    //pmB = new TPolyMarker();
    //pmB.SetMarkerStyle(5); // X
    //pmB.SetMarkerColor(6); // purple
    plA0.SetPolyline(sA0->GetNvert());
    plA0.SetLineColor(1); // black
    plA0.SetLineStyle(1);
    plA1.SetPolyline(sA1->GetNvert());
    plA1.SetLineColor(2); // red
    plA1.SetLineStyle(1);
    plB0.SetPolyLine(sB0.GetNvert());
    plB0.SetLineColor(3); // Green
    plB0.SetLineStyle(2);
    plB1.SetPolyLine(sB1->GetNvert());
    plB1.SetLineColor(4); // Blue
    plB1.SetLineStyle(2);
    //for(i=0;i<kNRadii;i++) pmA.SetPoint(i,xyB1p[i][0],xyB1p[i][1]);
    //for(i=0;i<kNRadii;i++) pmB.SetPoint(i,xyB1p[i][0],xyB1p[i][1]);
    for(i=0;i<sA0->GetNvert();i++) plA0.SetPoint(i,sA0->GetX(i),sA0->GetY(i));
    for(i=0;i<sA1->GetNvert();i++) plA1.SetPoint(i,sA1->GetX(i),sA1->GetY(i));
    for(i=0;i<sB0->GetNvert();i++) plB0.SetPoint(i,sB0->GetX(i),sB0->GetY(i));
    for(i=0;i<sB1->GetNvert();i++) plB1.SetPoint(i,sB1->GetX(i),sB1->GetY(i));
    canvas = new TCanvas("AliITSv11GeometrySPDFig0","",1000,1000);
    canvas.Range(-3.,-3.,3.,3.);
    txt.SetTextsize(0.05);
    txt.SetTextAlign(33);
    txt.SetTextColor(1);
    txt.Draw(2.9,2.9,"Section A-A outer Carbon Fiber surface");
    txt.SetTextColor(2);
    txt.Draw(2.9,2.5,"Section A-A Inner Carbon Fiber surface");
    txt.SetTextColor(3);
    txt.Draw(2.9,2.1,"Section E-E outer Carbon Fiber surface");
    txt.SetTextColor(4);
    txt.Draw(2.9,1.7,"Section E-E Inner Carbon Fiber surface");
    plA0.Draw();
    plA1.Draw();
    plB0.Draw();
    plB1.Draw();
    //pmA.Draw();
    //pmB.Draw();
    //
    Char_t chr[3];
    for(i=0;i<kNRadii;i++){
        sprintf(chr,"%2d",i);txt.DrawLatex(x-0.1,y,chr);
        sprintf(chr,"%8.4",);txt.DrawLatex(x,y,chr);
        sprintf(chr,"%8.4",);txt.DrawLatex(x+0.5,y,chr);
        sprintf(chr,"%8.4",);txt.DrawLatex(x+1.0,y,chr);
        sprintf(chr,"%8.4",);txt.DrawLatex(x+1.5,y,chr);
        sprintf(chr,"%8.4",);txt.DrawLatex(x+2.0,y,chr);
        if() txt.DrawLatex(x+2.5,y,"A-A/E-E");
        else txt.DrawLatex(x+2.5,y,"E-E");
    } // end for i
    txt.DrawLatex(x,y,"x_{c} mm");
    txt.DrawLatex(x+0.5,y,"y_{c} mm");
    txt.DrawLatex(x+1.0,y,"R mm");
    txt.DrawLatex(x+1.5,y,"#theta_{start}^{#circle}");
    txt.DrawLatex(x+2.0,y,"#theta_{end}^{#circle}");
    txt.DrawLatex(x+2.5,y,"Section");
    //
}
