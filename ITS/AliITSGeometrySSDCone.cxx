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
$Log$
$Id$
*/

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TTUBS.h>
#include <TPCON.h>
#include <TFile.h>    // only required for Tracking function?
#include <TCanvas.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TClonesArray.h>
#include <TBRIK.h>
#include <TSystem.h>
#include <TVector3.h>

#include "AliITSGeometrySSDCone.h"

ClassImp(AliITSGeometrySSDCone)

//______________________________________________________________________
AliITSGeometrySSDCone::AliITSGeometrySSDCone(){
    //Default Constructor for SSD Cone geometry

    SetScalemm();
}
//______________________________________________________________________
AliITSGeometrySSDCone::AliITSGeometrySSDCone(TVector3 *&tran,
					     const char moth[3],Int_t mat0){
    //Standard Constructor for SSD Cone geometry
    // Inputs:
    //   Double_t z0  Z-axis shift of this volume
    // Outputs:
    //   none.
    // Return:
    //   none.
    Double_t t; // some general angle and coordinates [degrees].

    th = 13.0; //mm, Thickness of Rohacell+carbon fiber
    ct=1.5; //mm, Carbon finber thickness
    r=15.0; // mm, Radius of curvature.
    tc=51.0; // angle of SSD cone [degrees].
    sintc=Sind(tc);costc=Cosd(tc);tantc=Tand(tc);
    z0=0.0;zcylinder=170.0;zpost=196.0;
    Routmax=0.5*985.0;RoutHole=0.5*965.0;Routmin=0.5*945.0;
    Rholemax=0.5*890.0;Rholemin=0.5*740.0;
    RPostmin=316.0;dRPost=23.0;zpostmax=196.0;phi0post=30.0;
    Rinmax=0.5*590.0;Rincylinder=0.5*597.0;RinHole=0.5*575.0;
    Rinmin=0.5*562.0;dzin=15.0;
    nspoaks=12;ninscrews=40;npost=4;nmounts=4;
    SSDcf=mat0+1; // SSD support cone Carbon Fiber materal number.
    SSDfs=mat0+2; // SSD support cone inserto stesalite 4411w.
    SSDfo=mat0+3; // SSD support cone foam, Rohacell 50A.
    SSDsw=mat0+4; // SSD support cone screw material,Stainless steal
    ncse=0; // number of screw ends (copy number)
    ncpe=0; // number of pin end (copy number)
    ncst=0; // number of screw tops (copy number)

    SetScalemm();
    // Lets start with the upper left outer carbon fiber surface.
    // Between za[2],rmaxa[2] and za[4],rmaxa[4] there is a curved section
    // given by rmaxa = rmaxa[2]-r*Sind(t) for 0<=t<=tc and 
    // za = za[2] + r*Cosd(t) for 0<=t<=tc. Simularly between za[1],rmina[1
    // and za[3],rmina[3] there is a curve section given by
    // rmina = rmina[1]-r*Sind(t) for 0<=t<=tc and za = za[1]+r&Sind(t)
    // for t<=0<=tc. These curves have been replaced by straight lines
    // between the equivelent points for simplicity.
    Double_t dza = th/sintc-(Routmax-Routmin)/tantc;
    if(dza<=0){ // The number or order of the points are in error for a proper
	// call to pcons!
	Error("SSDcone","The definition of the points for a call to PCONS is"
	      " in error. abort.");
	return;
    } // end if
    dphia=360.0;
    phi0a= 0.0;
    za[0]    = z0;
    rmina[0] = Routmin;
    rmaxa[0] = Routmax;
    za[1]    = za[0]+13.5-5.0 - dza; // za[2] - dza.
    rmina[1] = rmina[0];
    rmaxa[1] =rmaxa[0];
    za[2]    = za[0]+13.5-5.0; // From Drawing ALR-0767 and ALR-0767/3
    rmaxa[2] = rmaxa[0];
    za[3]    = za[1]+r*sintc;
    rmina[3] = rmina[1]-r*sintc;
    rmina[2] = rmina[1]+(rmina[3]-rmina[1])*(za[2]-za[1])/(za[3]-za[1]);
    za[4]    = za[2]+r*sintc;
    rmaxa[4] = rmaxa[2]-r*sintc;
    rmaxa[3] = rmaxa[2]+(rmaxa[4]-rmaxa[2])*(za[3]-za[2])/(za[4]-za[2]);
    rmina[5] = Rholemax;
    za[5]    = za[3]+(za[4]-za[3])*(rmina[5]-rmina[3])/(rmina[4]-rmina[3]);
    rmina[4] = rmina[3]+(rmina[5]-rmina[3])*(za[4]-za[3])/(za[5]-za[3]);
    za[6]    = th/sintc+za[5];
    rmina[6] = Rholemax;
    rmaxa[6] = rmina[6];
    rmaxa[5] = rmaxa[4]+(rmaxa[6]-rmaxa[4])*(za[5]-za[4])/(za[6]-za[4]);
    //
    // Now lets define the Inserto Stesalite 4411w material volume.
    dphib=360.0;
    phi0b= 0.0;
    zb[0] = z0;
    rminb[0] = rmina[0]+ct;
    rmaxb[0] = rmaxa[0]-ct;
    zb[1] = za[1];
    rminb[1] = rminb[0];
    rmaxb[1] = rmaxb[0];
    zb[2] = za[2];
    rmaxb[2] = rmaxb[1];
    zb[3] = za[4] - ct/sintc;
    rmaxb[3] = rmaxb[2] - (r-ct)*sintc;
    zb[4] = za[3]+ct/sintc;
    rminb[4] = rminb[1]-(r-ct)*sintc;
    rminb[2] = rminb[1]+(rminb[4]-rminb[1])*(zb[2]-zb[1])/(zb[4]-zb[1]);
    rminb[3] = rminb[1]+(rminb[4]-rminb[1])*(zb[3]-zb[1])/(zb[4]-zb[1]);
    zb[5] = zb[4]+(ct-2.*ct)/sintc;
    rminb[5] = rminb[4]+(ct-2.*ct)*tantc;
    rmaxb[5] = rminb[5];
    rmaxb[4] = rmaxb[3]+(rmaxb[5]-rmaxb[3])*(zb[4]-zb[3])/(zb[5]-zb[3]);
    //
    // Now lets define the Rohacell foam material volume.
    dphic=360.0;
    phi0c= 0.0;
    zc[0] = zb[4];
    rminc[0] = rminb[4];
    rmaxc[0] = rminc[0];
    zc[1] = zb[5];
    rmaxc[1] = rminb[5];
    zc[2] = za[5] + ct/sintc;
    rminc[2] = rmina[5]+ct; // leave space for carbon fiber covering hole.
    rminc[1] = rminc[0] +(rminc[2]-rminc[0])*(zc[1]-zc[0])/(zc[2]-zc[0]);
    zc[3] = za[6] - ct/sintc;
    rminc[3] = rmina[6]+ct;
    rmaxc[3] = rminc[3];
    rmaxc[2] = rmaxc[1]+(rmaxc[3]-rmaxc[1])*(zc[2]-zc[1])/(zc[3]-zc[1]);
    //
    // In volume SCB, th Inserto Stesalite 4411w material volume, there
    // are a number of Stainless steel screw and pin studs which will be
    // filled with screws/studs.
    rmine=0.0,rmaxe=6.0,dze=0.5*10.0; // mm
    rmine2=0.0;rmaxe2=6.0;dze2=0.5*12.0; // mm
    //
    // There is no carbon fiber between this upper left section and the
    // SSD spoaks. We remove it by replacing it with Rohacell foam.
    t = ct/(0.5*(Rholemax+Rholemin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    dphif = 5.0 - 2.0*t; // degrees
    phi0f = 12.5+t; // degrees see drawing ALR-0767.
    zf[0] = zc[2];
    rminf[0] = rminc[3];
    rmaxf[0] = rminf[0];
    rminf[1] = rmina[5];
    rmaxf[1] = rminf[0];
    zf[1] = zc[0]+(zc[2]-zc[0])*(rminf[1]-rminc[0])/(rminc[2]-rminc[0]);
    zf[2] = zc[3];
    rminf[2] = rminf[1];
    rmaxf[2] = rmaxf[1];
    zf[3] = zc[1]+(zc[3]-zc[1])*(rmaxf[3]-rmaxc[1])/(rmaxc[3]-rmaxc[1]);
    rminf[3] = rmina[5];
    rmaxf[3] = rminf[3];
    //=================================================================
    // Now for the spoak part of the SSD cone.
    // It is not posible to inclue the radius of curvature between
    // the spoak part and the upper left part of the SSD cone or lowwer right
    // part. This would be discribed by the following curves.
    // R = Rmax - (5mm)*Sin(t) phi = phi0+(5mm*180/(Pi*RoutHole))*Sin(t) 
    // where 0<=t<=90 For the inner curve a simular equiation holds.
    phi0g = 12.5; // degrees see drawing ALR-0767.
    dphig = 5.0; // degrees
    zg[0] = zb[5];
    rming[0] = rmina[5];
    rmaxg[0] = rming[0];
    zg[1] = za[6];
    rming[1] = -tantc*(zg[1]-za[3])+rmina[3];
    rmaxg[1] = rmaxg[0];
    rming[2] = Rholemin;
    zg[2] = za[3]-(rming[2]-rmina[3])/tantc;
    rmaxg[2] = -tantc*(zg[2]-za[4])+rmaxa[4];
    rming[3] = rming[2];
    rmaxg[3] = rming[3];
    zg[3] = za[4]-(rmaxg[3]-rmaxa[4])/tantc;
    // For the foam core.
    t = ct/(0.5*(Rholemax+Rholemin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    dphih = 5.0 - 2.0*t; // degrees
    phi0h = 12.5+t; // degrees see drawing ALR-0767.
    zh[0] = zf[2];
    rminh[0] = rming[0];
    rmaxh[0] = rmaxg[0];
    zh[1] = zf[3];
    rminh[1] = rming[1]-(ct/sintc-(zg[1]-zh[1]))*tantc;
    rmaxh[1] = rmaxh[0];
    zh[2] = zg[2]+ct/tantc;
    rminh[2] = rming[2];
    rmaxh[2] = rmaxg[2]-(ct/sintc-(zg[2]-zh[2]))*tantc;
    zh[3] = zg[3]-ct/sintc;
    rminh[3] = rminh[2];
    rmaxh[3] = rminh[3];
    //
    //==================================================================
    // Now for the Inner most part of the SSD cone.
    phi0i = 0.0;
    dphii = 360.0;
    Double_t za,rmina,rmaxa; // additional point not needed in call to pcons.
    zi[0] = zg[2];
    rmini[0] = rming[2];
    rmaxi[0] = rmini[0];
    zi[1] = zg[3];
    rmini[1] = -tantc*(zi[1]-za[3])+rmina[3];
    rmaxi[1] = rmaxi[0];
    rmini[5] = Rinmin;
    rmaxi[5] = Rinmax+r*sintc;
    zi[5] =za[4]+(rmaxa[4]-rmaxi[5])/tantc;
    za = zi[5]+r*costc;
    rmina = rmini[5];
    rmaxa = Rinmax;
    zi[3] = za-dzin;
    zi[2] = zi[3] -r*costc;
    rmini[2] = -tantc*(zi[2]-za[3])+rmina[3];
    rmaxi[2] = -tantc*(zi[2]-za[4])+rmaxa[4];
    rmini[3] = rmini[2] -r*costc;
    zi[4] = zi[3];
    rmini[4] = Rinmin;
    rmaxi[4] = -tantc*(zi[4]-za[4])+rmaxa[4];
    rmaxi[3] = rmaxi[4];
    zi[6] = zcylinder;
    rmini[6] = Rinmin;
    rmaxi[6] = rmaxi[5] - (zi[5]-zi[6])*(rmaxi[5]-rmaxa)/(zi[5]-za);
    zi[7] = zi[6];
    rmini[7] = Rincylinder;
    rmaxi[7] = rmaxi[6];
    rmini[8] = Rincylinder;
    rmaxi[8] = rmini[8];
    zi[8] = zi[5]+(rmaxi[8]-rmaxi[5])*(za-zi[5])/(rmaxa-rmaxi[5]);
    // Now for Inserto volume at the inner most radius.
    phi0k = 0.0;
    dphik = 360.0;
    zk[1] = zi[3]+ct;
    zk[0] = zk[1]-(r+ct)*costc;
    rmink[0] = rmini[3]+(r+ct)*sintc;
    rmaxk[0] = rmink[0];
    rmink[1] = rmini[3];
    zk[2] = zk[1];
    rmink[2] = rmini[6];
    rmaxk[2] = rmaxk[1];
    zk[3] = zk[0]+(th+2.0*ct)*costc;
    rmink[3] = rmini[6];
    rmaxk[3] = rmaxk[0]+(th+2.0*ct)*sintc;
    rmaxk[1] = rmaxk[0]+(rmaxk[3]-rmaxk[0])*(zk[1]-zk[0])/(zk[3]-zk[0]);
    rmink[4] = rmini[6];
    rmaxk[4] = rmaxi[5]-ct*sintc;
    zk[4] = zc[1]+(zc3[3]-zc[1])*(rmaxk[4]-rmaxc[1])/(rmaxc[3]-rmaxc[1]);
    zk[5] = zi[5]-r*costc-ct;
    rmink[5] = rmini[6];
    rmaxk[5] = rmini[8];
    zk[6] = zi[6];
    rmink[6] = rmini[6];
    rmaxk[6] = rmaxi[6];
    // Now for foam core at the inner most radius.
    phi0j = 0.0;
    dphij = 360.0;
    rminj[0] = rmini[0]-ct;
    zj[0] = zc[0]+(zc[2]-zc[0])*(rminj[0]-rminc[0])/(rminc[2]-rminc[0]);
    rmaxj[0] = rminj[0];
    rmaxj[1] = rmaxj[0];
    zj[1] = zc[1]+(zc[3]-zc[1])*(rmaxj[1]-rmaxc[1])/(rmaxc[3]-rmaxc[1]);
    rminj[1] = rminc[0]+(rminc[2]-rminc[0])*(zj[1]-zc[0])/(zc[2]-zc[0]);
    zj[2] = zk[0];
    rminj[2] = rmink[0];
    rmaxj[2] = rmaxc[1]+(rmaxc[3]-rmaxc[1])*(zj[2]-zc[1])/(zc[3]-zc[1]);
    zj[3] = zk[3];
    rminj[3] = rmaxk[3];
    rmaxj[3] = rminj[3];
    // Now for foam core at the top of the inner most radius where 
    // the spoaks are.
    t = ct/(0.5*(Rholemax+Rholemin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    dphil = 5.0 - 2.0*t; // degrees
    phi0l = 12.5+t; // degrees see drawing ALR-0767.
    zl[0] = zh[2];
    rminl[0] = rmini[0];
    rmaxl[0] = rminl[0];
    zl[1] = zj[0];
    rminl[1] = rminj[1];
    rmaxl[1] = rmaxi[0];
    zl[2] = zh[3];
    rminl[2] = rminl[1];
    rmaxl[2] = rmaxl[1];
    zl[3] = zj[1];
    rminl[3] = rminl[2];
    rmaxl[3] = rminl[3];
    // Now for the SSD mounting posts
    dphio = 180.0*dRPost/(RPostmin+0.5*dRPost)/TMath::Pi(); // degrees
    phi0o = phi0post; //
    rmino[0] = RPostmin+dRPost;
    rmaxo[0] = rmino[0];
    zo[0] = za[4]+(rmaxa[4]-rmaxo[0])/tantc;
    rmino[1] = RPostmin;
    zo[1] = za[4]+(rmaxa[4]-rmino[1])/tantc;
    rmaxo[1] = rmaxo[0];
    zo[2] = z0+zpostmax;
    rmino[2] = RPostmin;
    rmaxo[2] = rmino[2]+dRPost;
    // Now for the SSD mounting posts
    t = 180.0*ct/(RPostmin+0.5*dRPost)/TMath::Pi();
    dphip = dphio-2.0*t; // degrees
    phi0p = phio0+t; //
    rminp[0] = rmino[0]-ct;
    rmaxp[0] = rminp[0];
    zp[0] = za[4]+(rmaxa[4]-rmaxp[0])/tantc;
    rminp[1] = rmino[0]+ct;
    rmaxp[1] = rmino[0]-ct;
    zp[1] = za[4]+(rmaxa[4]-rminp[1])/tantc;
    rminp[2] = rminp[1];
    rmaxp[2] = rmaxp[1];
    zp[2] = z0-zpostmax;
    // This insrto continues into the SSD cone displacing the foam
    // and the carbon fiber surface at those points where the posts are.
    dphim=360.0;
    phi0m= 0.0;
    rminm[0] = RPostmin+dRPost-ct;
    rmaxm[0] = rminm[0];
    zm[0] = zj[0]+(zj[2]-zj[0])*(rminm[0]-rminj[0])/(rminj[2]-rminj[0]);
    rmaxm[1] = rmaxm[0];
    zm[1] = zj[1]+(zj[3]-zj[1])*(rmaxm[1]-rmaxj[1])/(rmaxj[3]-rmaxj[1]);
    rminm[2] = RPostmin+ct;
    zm[2] = zj[0]+(zj[2]-zj[0])*(rminm[2]-rminj[0])/(rminj[2]-rminm[0]);
    rmaxm[2] = rmaxj[1]+(rmaxj[3]-rmaxm[1])*(zm[2]-zj[1])/(zj[3]-zj[1]);
    rminm[3] = rminm[2];
    rmaxm[3] = rminm[3];
    dphin=360.0;
    phi0n= 0.0;
    zn[0] = zm[1];
    rminn[0] = rmaxm[1];
    rmaxn[0] = rminn[0];
    rmaxn[1] = rmaxn[0];
    zn[1] = za[4]+(rmaxa[4]-rmaxn[1])/tantc;
    rminn[1] = rmaxj[1]+(rmaxj[3]-rmaxj[1])*(zn[1]-zj[1])/(zj[3]-zj[1]);
    zn[2] = zm[3];
    rminn[2] = rminm[3];
    rmaxn[2] = -tantc*(zn[2]-za[4])+rmaxa[4];
    rminn[3] = rminn[2];
    rmaxn[3] = rminn[3];
    zn[3] = za[4]+(rmaxa[4]-rmaxn[3])/tantc;
}
//______________________________________________________________________
void AliITSGeometrySSDCone::CreateG3Geometry(const char moth[3],
					     TVector3 &trans){
    // Calls Geant 3 geometry inilization routines with the information
    // stored in this class.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.

    PolyCone("SCA","SSD Suport cone Carbon Fiber Surface outer left",
	     phi0a,dphia,nza,za,rmina,rmaxa,SSDcf);
    Pos("SCA",1,moth,trans.x(),trans.y(),trans.z(),0);
    XMatrix(1,180.0);
    Pos("SCA",2,moth,trans.x(),trans.y(),-trans.z(),1);
    PolyCone("SCB","SSD Suport cone Inserto Stesalite left edge",
	     phi0b,dphib,nzb,zb,rminb,rmaxb,SSDfs);
    Pos("SCB",1,"SCA",0.0,.0,0.0,0);
    PolyCone("SCC","SSD Suport cone Rohacell foam left edge",
	     phi0,dphi,nz,zc,rminc,rmaxc,SSDfo);
    Pos("SCC",1,"SCA",0.0,.0,0.0,0);
    Tube("SCD","Screw+stud used to mount things to the SSD support cone",
	 rmine,rmaxe,dze,SSDsw);
    Tube("SCE","pin used to mount things to the SSD support cone",
	 rmine2,rmaxe2,dze2,SSDsw);
    k=l=0;
    for(i=0;i<2;i++){ // position for ITS-TPC mounting brackets
	for(j=0;j<2;j++){ // 2 screws per bracket
	    ncse++;
	    t = -5.0+10.0*((Double_t)j)+180.*((Double_t)i);
	    x = RoutHole*Sind(t);
	    y = RoutHole*Cosd(t);
	    z = dz;
	    Pos("SCD",ncse,"SCB",x,y,z,0);
	} // end for j
	for(j=0;j<3;j++){ // 3 pins per bracket
	    ncpe++;
	    t = -3.0+3.0*((Double_t)j)+180.*((Double_t)i);
	    x = RoutHole*Sind(t);
	    y = RoutHole*Cosd(t);
	    z = dz;
	    Pos("SCE",ncpe,"SCB",x,y,z,0);
	} // end for j
    } // end for i
    for(i=0;i<2;i++){ // position for ITS-rail mounting brackets
	for(j=0;j<4;j++){ // 4 screws per bracket
	    a[4]={0.0,2.0,5.0,7.0}; // Relative angles.
	    ncse++;
	    t = 90.0-a[j]+187.*((Double_t)i);
	    x = RoutHole*Sind(t);
	    y = RoutHole*Cosd(t);
	    z = dz;
	    Pos("SCD",kncs,"SCB",x,y,z,0);
	} // end for j
	for(j=0;j<2;j++){ // 2 pins per bracket
	    ncpe++;
	    t = 88+7.0*((Double_t)j)+184.*((Double_t)i);
	    x = RoutHole*Sind(t);
	    y = RoutHole*Cosd(t);
	    z = dz;
	    Pos("SCE",ncse,"SCB",x,y,z,0);
	} // end for j
    } // end for i
    for(i=0;i<nmounts;i++){ // mounting holes/screws for beam pipe support
	// and SPD cone support (dump side,non-dump side has them to).
	for(j=0;j<2;j++){ // 2 screws per bracket
	    ncse++;
	    t = 180.*20./(RoutHole*TMath::Pi());
	    t = 45.0+((Doulbe_t)(j-1))*t+90.*((Double_t)i);
	    x = RoutHole*Sind(t);
	    y = RoutHole*Cosd(t);
	    z = dz;
	    Pos("SCD",ncse,"SCB",x,y,z,0);
	} // end for j
	for(j=0;j<1;j++){ // 1 pins per bracket
	    ncpe++;
	    t = 45.0+90.*((Double_t)i);
	    x = RoutHole*Sind(t);
	    y = RoutHole*Cosd(t);
	    z = dz;
	    Pos("SCE",ncpe,"SCB",x,y,z,0);
	} // end for j
    } // end for i
    PolyCone("SCF","SSD Suport cone Rohacell foam left edge",
	     phi0f,dphif,nzf,zf,rminf,rmaxf,SSDfo);
    Pos("SCF",1,"SCA",0.0,.0,0.0,0);
    for(i=1;i<nspoaks;i++){
	Zmatrix(irot+i,360./((Double_t)nspoaks));
	Pos("SCG",i+1,"SCA",0.0,.0,0.0,irot+i);
    } // end for i
    PolyCone("SCG","SSD spoak carbon fiber surfaces",
	     phi0g,dphig,nzg,zg,rming,rmaxc,SSDcf);
    Pos("SCG",i+1,"SCA",0.0,.0,0.0,0);
    for(i=1;i<nspoaks;i++){
	Pos("SCG",i+1,"SCA",0.0,.0,0.0,irot+i);
    } // end for i
    PolyCone("SCH","SSD spoak foam core",
	     phi0h,dphih,nzh,zh,rminh,rmaxh,SSDfo);
    Pos("SCH",1,"SCG",0.0,.0,0.0,0);
    PolyCone("SCI","SSD lower/inner right part of SSD cone",
	     phi0i,dphii,nzi,zci,rminci,rmaxci,SSDcf);
    Pos("SCI",1,moth,0.0,.0,0.0,0);
    PolyCone("SCK","SSD inner most inserto material",
	     phi0k,dphik,nzk,zk,rmink,rmaxk,SSDfs);
    Pos("SCK",1,"SCI",0.0,.0,0.0,0);
    PolyCone("SCJ","SSD inner most foam core",
	     phi0j,dphij,nzj,zj,rminj,rmaxj,SSDfo);
    Pos("SCJ",1,"SCI",0.0,.0,0.0,0);
    PolyCone("SCL","SSD inner most foam core",
	     phi0l,dphil,nzl,zl,rminl,rmaxl,SSDfo);
    Pos("SCL",1,"SCI",0.0,.0,0.0,0);
    for(i=1;i<nspoaks;i++){
	Pos("SCG",i+1,"SCA",0.0,.0,0.0,irot+i);
    } // end for i
    PolyCone("SCO","SSD mounting post, carbon fiber",
	     phi0,dphi,nz,zc,rminc,rmaxc,SSDcf);
    Pos("SCO",1,moth,0.0,.0,0.0,0);
    for(i=1;i<nposts;i++){
	Zmatrix(irot+i,360./((Double_t)nposts));
	Pos("SCO",i+1,moth,0.0,.0,0.0,irot+i);
    } // end for
    PolyCone("SCP","SSD mounting post, Inserto",
	     phi0p,dphip,nzp,zp,rminp,rmaxp,SSDfs);
    Pos("SCP",1,"SCO",0.0,.0,0.0,0);
    Pos("SCM",1,"SCJ",0.0,.0,0.0,0);
    Pos("SCN",1,"SCI",0.0,.0,0.0,0);
    for(i=1;i<nposts;i++){
	Pos("SCN",i+1,"SCJ",0.0,.0,0.0,irot+i);
	Pos("SCM",i+1,"SCI",0.0,.0,0.0,irot+i);
    } // end for i
    return;
}
//______________________________________________________________________
void CreateG3Materials(){
    // Fills the Geant 3 banks with Material and Medium definisions.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Returns:
    //   none.
    Double_t Z[5],W[5],dens;

    Z[0] = 1.; W[0] = 0.5; // Hydrogen Content
    Z[1] = 6.; W[1] = 0.5; // Carbon Content
    MixtureByWeight(SSDcf,"Carbon Fiber for SSD support cone",Z,W,dens,2);
    Z[0] = 1.; W[0] = 0.5; // Hydrogen Content
    Z[1] = 6.; W[1] = 0.5; // Carbon Content
    MixtureByWeight(SSDfs,"Inserto stealite 4411w for SSD support cone",
		    Z,W,dens,2);
    Z[0] = 1.; W[0] = 0.5; // Hydrogen Content
    Z[1] = 6.; W[1] = 0.5; // Carbon Content
    MixtureByWeight(SSDfo,"Foam core (Rohacell 50A) for SSD support cone",
		    Z,W,dens,2);
    Z[0] =  6.; W[0] = 0.5; // Carbon Content
    Z[1] = 25.; W[1] = 0.5; // Iron Content
    MixtureByWeight(SSDsw,"Stainless steal screw, pin, and stud material",
		    Z,W,dens,2);
}
//______________________________________________________________________
void AliITSGeometrySSDCone::BuildDisplayGeometry(){
    // Fill Root geometry banks for fast simple ITS simulation event
    // display. See Display.C, and related code, for more details.
    // Inputs:
    //    none.
    // Outputs:
    //   none.
    // Return:
    //  none.

    // No need to display ITS cones.
}
