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
Revision 1.1  2003/02/10 17:03:52  nilsen
New version and structure of ITS V11 geometry. Work still in progress.

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
#include <AliRun.h>
#include <AliITS.h>

#include "AliITSGeometrySSDCone.h"

ClassImp(AliITSGeometrySSDCone)

//______________________________________________________________________
AliITSGeometrySSDCone::AliITSGeometrySSDCone() : AliITSBaseGeometry(){
    //Default Constructor for SSD Cone geometry

    SetScalemm();
}
//______________________________________________________________________
AliITSGeometrySSDCone::AliITSGeometrySSDCone(AliITS *its,TVector3 &tran,
					     const char *moth,Int_t mat0):
    AliITSBaseGeometry(its,0){
    //Standard Constructor for SSD Cone geometry
    // Inputs:
    //   Double_t z0  Z-axis shift of this volume
    // Outputs:
    //   none.
    // Return:
    //   none.
    Double_t t; // some general angle and coordinates [degrees].
    Double_t Z,Rmin,Rmax; // additional point not needed in call to pcons.

    fThickness = 13.0; //mm, Thickness of Rohacell+carbon fiber
    fCthick=1.5; //mm, Carbon finber thickness
    fRcurv=15.0; // mm, Radius of curvature.
    fTc=51.0; // angle of SSD cone [degrees].
    fSintc=Sind(fTc);fCostc=Cosd(fTc);fTantc=Tand(fTc);
    fZ0=0.0;fZouterMilled=13.5-5.0;fZcylinder=170.0;fZposts=196.0;
    fNspoaks=12;fNposts=4;fNmounts=4;
    fRoutMax=0.5*985.0;fRoutHole=0.5*965.0;fRoutMin=0.5*945.0;
    fRholeMax=0.5*890.0;fRholeMin=0.5*740.0;
    fRpostMin=316.0;fdRpost=23.0;fZpostMax=196.0;fPhi0Post=30.0;
    fRinMax=0.5*590.0;fRinCylinder=0.5*597.0;fRinHole=0.5*575.0;
    fRinMin=0.5*562.0;fdZin=15.0;
    // SSD-SDD Thermal/Mechanical cylinder mounts
    fNinScrews=40;
    fPhi0Screws=0.5*360.0/((Double_t)fNinScrews);fRcylinderScrews= 0.5*570.0;
    fDscrewHead=8.0;fDscrewShaft=4.6;fThScrewHeadHole=8.5;
    // SDD mounting bracket, SSD part
    fNssdSupports=3;fPhi0SDDsupports=90.0;
    fRsddSupportPlate = 0.5*585.0;fThSDDsupportPlate=4.0;
    fWsddSupportPlate = 70.0;
    fSSDcf=26; // SSD support cone Carbon Fiber materal number.
    fSSDfs=25; // SSD support cone inserto stesalite 4411w.
    fSSDfo=68; // SSD support cone foam, Rohacell 50A.
    fSSDsw=14; // SSD support cone screw material,Stainless steal
    fNcD=0; // number of screw ends (copy number)
    fNcE=0; // number of pin end (copy number)

    SetScalemm();
    // Poly-cone Volume A. Top part of SSD cone Carbon Fiber.
    fA.Size(7,"SSD Suport cone Carbon Fiber Surface outer left");
     // Poly-cone Volume B. Stesalite inside volume A.
    fB.Size(6,"SSD Suport cone Inserto Stesalite left edge");
    // Poly-cone Volume C. Foam inside volume A.
    fC.Size(4,"SSD Suport cone Rohacell foam left edge");
    fD.SetName("Screw+stud used to mount things to the SSD support cone");
    fE.SetName("pin used to mount things to the SSD support cone");
    // Poly-cone Volume F. Foam in spoak reagion, inside volume A.
    fF.Size(4,"SSD Top Suport cone Rohacell foam Spoak");
    fG.Size(4,"SSD spoak carbon fiber surfaces"); // Poly-cone Volume G.
    fH.Size(4,"SSD support cone Rohacell foam Spoak"); // Poly-cone Volume H.
    fI.Size(9,"SSD lower/inner right part of SSD cone"); //Poly-cone Volume I.
    fJ.Size(4,"SSD inner most foam core"); // Poly-cone Volume J.
    fK.Size(7,"SSD inner most inserto material"); // Poly-cone Volume K.
    fL.Size(4,"SSD Bottom cone Rohacell foam Spoak"); // Poly-cone Volume L.
    fM.Size(4,"SSD mounting post foam substitute, Inserto");//Poly-cone Vol. M
    fN.Size(4,"SSD mounting post CF subsititute, Inserto");//Poly-cone Vol. N
    fO.Size(3,"SSD mounting post, carbon fiber"); // Poly-cone Volume O.
    fP.Size(3,"SSD mounting post, Inserto"); // Poly-cone Volume P.
    fQ.Size(4,"SSD Thermal sheal stainless steel bolts");//Poly-cone Volume Q.
    fR.SetName("Air in front of bolt (in stasolit)");
    fS.SetName("Air in front of Stainless Steal Screw end, N6");
    fT.Size(2,"SSD-SDD mounting bracket Inserto-> Al."); //Poly-cone Volume T.
    fU.Size(4,"SSD-SDD mounting bracket CF->Al."); // Poly-cone Volume U.
    // Lets start with the upper left outer carbon fiber surface.
    // Between za[2],rmaxa[2] and za[4],rmaxa[4] there is a curved section
    // given by rmaxa = rmaxa[2]-r*Sind(t) for 0<=t<=fTc and 
    // za = za[2] + r*Cosd(t) for 0<=t<=fTc. Simularly between za[1],rmina[1
    // and za[3],rmina[3] there is a curve section given by
    // rmina = rmina[1]-r*Sind(t) for 0<=t<=fTc and za = za[1]+r&Sind(t)
    // for t<=0<=fTc. These curves have been replaced by straight lines
    // between the equivelent points for simplicity.
    Double_t dza = fThickness/fSintc-(fRoutMax-fRoutMin)/fTantc;
    if(dza<=0){ // The number or order of the points are in error for a proper
	// call to pcons!
	Error("SSDcone","The definition of the points for a call to PCONS is"
	      " in error. abort.");
	return;
    } // end if
    fA.P0() = 0.0;
    fA.dP() = 360.0;
    fA.Z(0) = fZ0;
    fA.Rn(0) = fRoutMin;
    fA.Rx(0) = fRoutMax;
    fA.Z(1)  = fA.ZAt(0)+fZouterMilled - dza; // za[2] - dza.
    fA.Rn(1) = fA.Rmin(0);
    fA.Rx(1) = fA.Rmax(0);
    fA.Z(2)  = fA.ZAt(0)+fZouterMilled; //From Drawing ALR-0767 and ALR-0767/3
    fA.Rx(2) = fA.Rmax(0);
    RadiusOfCurvature(fRcurv,0.0,fA.ZAt(1),fA.Rmin(1),fTc,fA.Z(3),fA.Rn(3));
    fA.Rn(2) = RminFrom2Points(fA,3,1,fA.ZAt(2));
    RadiusOfCurvature(fRcurv,0.0,fA.ZAt(2),fA.Rmax(2),fTc,fA.Z(4),fA.Rx(4));
    fA.Rn(4) = RminFromZSSDcone(fA.ZAt(4));
    fA.Rx(3) = RmaxFrom2Points(fA,4,2,fA.ZAt(3));
    fA.Rn(5) = fRholeMax;
    fA.Z(5)  = Zfrom2MinPoints(fA,4,3,fA.Rmin(5));
    fA.Rx(5) = RmaxFromZSSDcone(fA.ZAt(5));
    fA.Rn(6) = fRholeMax;
    fA.Rx(6) = fA.Rmin(6);
    fA.Z(6)  = ZFromRmaxSSDcone(fA.Rmax(6));
    //
    // Now lets define the Inserto Stesalite 4411w material volume.
    fB.P0() = 0.0;
    fB.dP() = 360.0;
    fB.Z(0) = fA.ZAt(0);
    fB.Rn(0) = fA.Rmin(0)+fCthick;
    fB.Rx(0) = fA.Rmax(0)-fCthick;
    fB.Z(1)  = fA.ZAt(1);
    fB.Rn(1) = fB.Rmin(0);
    fB.Rx(1) = fB.Rmax(0);
    fB.Z(2)  = fA.ZAt(2);
    fB.Rx(2) = fB.Rmax(1);
    RadiusOfCurvature(fRcurv-fCthick,0.,fB.ZAt(2),fB.Rmax(2),
		                    fTc,fB.Z(3),fB.Rx(3));
    RadiusOfCurvature(fRcurv+fCthick,0.,fB.ZAt(1),fB.Rmin(1),
		                    fTc,fB.Z(4),fB.Rn(4));
    fB.Rn(2) = RminFrom2Points(fB,4,1,fB.ZAt(2));
    fB.Rn(3) = RminFrom2Points(fB,4,1,fB.ZAt(3));
    fB.Z(5)  = fB.ZAt(4)+(fThickness-2.0*fCthick)/fSintc;
    fB.Rn(5) = RmaxFromZSSDcone(fB.ZAt(5),-fCthick);
    fB.Rx(5) = fB.Rmin(5);
    fB.Rx(4) = RmaxFrom2Points(fB,5,3,fB.ZAt(4));
    //
    // Now lets define the Rohacell foam material volume.
    fC.P0() = 0.0;
    fC.dP() = 360.0;
    fC.Z(0) = fB.ZAt(4);
    fC.Rn(0) = fB.Rmin(4);
    fC.Rx(0) = fC.Rmin(0);
    fC.Z(1)  = fB.ZAt(5);
    fC.Rx(1) = fB.Rmin(5);
    fC.Rn(2) = fA.Rmin(5)+fCthick;//leave space for carbon fiber covering hole
    fC.Z(2)  = ZFromRminSSDcone(fC.Rn(2),+fCthick);
    fC.Rn(1) = RminFrom2Points(fC,2,0,fC.ZAt(1));
    fC.Rx(3) = fA.Rmin(6)+fCthick;
    fC.Rn(3) = fC.Rmax(3);
    fC.Z(3)  = ZFromRmaxSSDcone(fC.Rx(3),-fCthick);
    fC.Rx(2) = RmaxFrom2Points(fC,3,1,fC.ZAt(2));
    //
    // In volume SCB, th Inserto Stesalite 4411w material volume, there
    // are a number of Stainless steel screw and pin studs which will be
    // filled with screws/studs.
    fD.Rn()=0.0,fD.Rx()=6.0,fD.Z()=0.5*10.0; // mm
    fE.Rn()=0.0;fE.Rx()=6.0;fE.Z()=0.5*12.0; // mm
    //
    // There is no carbon fiber between this upper left section and the
    // SSD spoaks. We remove it by replacing it with Rohacell foam.
    t = fCthick/(0.5*(fRholeMax+fRholeMin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    fF.P0()  = 12.5+t; // degrees see drawing ALR-0767.
    fF.dP()  = 5.0 - 2.0*t; // degrees
    fF.Z(0)  = fC.ZAt(2);
    fF.Rn(0) = fC.Rmin(3);
    fF.Rx(0) = fF.Rmin(0);
    fF.Rn(1) = fA.Rmin(5);
    fF.Rx(1) = fF.Rmin(0);
    fF.Z(1)  = RminFromZSSDcone(fF.ZAt(1),+fCthick);
    fF.Z(2)  = fC.ZAt(3);
    fF.Rn(2) = fF.Rmin(1);
    fF.Rx(2) = fF.Rmax(1);
    fF.Rn(3) = fA.Rmin(6);
    fF.Rx(3) = fF.Rmin(3);
    fF.Z(3)  = ZFromRmaxSSDcone(fF.Rmax(3),-fCthick);
    //=================================================================
    // Now for the spoak part of the SSD cone.
    // It is not posible to inclue the radius of curvature between
    // the spoak part and the upper left part of the SSD cone or lowwer right
    // part. This would be discribed by the following curves.
    // R = Rmax - (5mm)*Sin(t) phi = phi0+(5mm*180/(Pi*fRoutHole))*Sin(t) 
    // where 0<=t<=90 For the inner curve a simular equiation holds.
    fG.P0()  = 12.5; // degrees see drawing ALR-0767.
    fG.dP()  = 5.0; // degrees
    fG.Z(0)  = fA.ZAt(5);
    fG.Rn(0) = fA.Rmin(5);
    fG.Rx(0) = fG.Rn(0);
    fG.Z(1)  = fA.ZAt(6);
    fG.Rn(1) = RminFromZSSDcone(fG.ZAt(1));
    fG.Rx(1) = fG.Rmax(0);
    fG.Rn(2) = fRholeMin;
    fG.Z(2)  = ZFromRminSSDcone(fG.Rmin(2));
    fG.Rx(2) = RmaxFromZSSDcone(fG.ZAt(2));
    fG.Rn(3) = fG.Rmin(2);
    fG.Rx(3) = fG.Rmin(3);
    fG.Z(3)  = ZFromRmaxSSDcone(fG.Rmax(3));
    // For the foam core.
    t = fCthick/(0.5*(fRholeMax+fRholeMin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    fH.P0()  = 5.0 - 2.0*t; // degrees
    fH.dP()  = 12.5+t; // degrees see drawing ALR-0767.
    fH.Z(0)  = fF.ZAt(1);
    fH.Rn(0) = fG.Rmin(0);
    fH.Rx(0) = fH.Rmin(0);
    fH.Z(1)  = fF.ZAt(3);
    fH.Rn(1) = RminFromZSSDcone(fH.Z(1),-fCthick);
    fH.Rx(1) = fH.Rmax(0);
    fH.Z(2)  = ZFromRminSSDcone(fG.Rmin(2),+fCthick);
    fH.Rn(2) = fG.Rmin(2);
    fH.Rx(2) = RmaxFromZSSDcone(fH.Z(2),-fCthick);
    fH.Z(3)  = ZFromRmaxSSDcone(fG.Rmin(3),-fCthick);
    fH.Rn(3) = fG.Rmin(3);
    fH.Rx(3) = fH.Rn(3);
    //
    //==================================================================
    // Now for the Inner most part of the SSD cone.
    fI.P0()  = 0.0;
    fI.dP()  = 360.0;
    fI.Z(0)  = fG.ZAt(2);
    fI.Rn(0) = fG.Rmin(2);
    fI.Rx(0) = fI.Rmin(0);
    fI.Z(1)  = fG.ZAt(3);
    fI.Rn(1) = RminFromZSSDcone(fI.ZAt(1));
    fI.Rx(1) = fI.Rmax(0);
    fI.Rn(4) = fRinMin;
    fI.Rn(5) = fRinMin;
    RadiusOfCurvature(fRcurv,90.0,0.0,fRinMin,90.0-fTc,Z,fI.Rx(5)); // z dummy
    fI.Z(5)  = ZFromRmaxSSDcone(fI.Rx(5));
    fI.Z(6)  = fZcylinder;
    fI.Rn(6) = fRinMin;
    fI.Z(7)  = fI.Z(6);
    fI.Rn(7) = fRinCylinder;
    fI.Rn(8) = fRinCylinder;
    fI.Rx(8) = fI.Rmin(8);
    Rmin = fI.Rmin(5);
    RadiusOfCurvature(fRcurv,90.0-fTc,fI.Z(5),fI.Rmax(5),90.0,Z,Rmax);
    Rmax = fRinMax;
    fI.Z(8)  = Z+(fI.ZAt(5)-Z)*(fI.Rmax(8)-Rmax)/(fI.Rmax(7)-Rmax);
    fI.Rx(6) = RmaxFrom2Points(fI,8,5,fI.ZAt(6));
    fI.Rx(7) = fI.Rmax(6);
    fI.Z(3)  = Z-fdZin;
    fI.Z(4)  = fI.ZAt(3);
    fI.Rx(3) = RmaxFromZSSDcone(fI.ZAt(3));
    fI.Rx(4) = fI.Rx(3);
    //rmin dummy
    RadiusOfCurvature(fRcurv,90.,fI.ZAt(3),0.,90.-fTc,fI.Z(2),Rmin);
    fI.Rn(2) = RminFromZSSDcone(fI.ZAt(2));
    fI.Rx(2) = RmaxFromZSSDcone(fI.ZAt(2));
    // z dummy
    RadiusOfCurvature(fRcurv,90.-fTc,0.0,fI.Rmin(2),90.0,Z,fI.Rn(3)); 
    // Now for Inserto volume at the inner most radius.
    fK.P0()  = 0.0;
    fK.dP()  = 360.0;
    fK.Z(1)  = fI.ZAt(3)+fCthick;
    fK.Rn(1) = fI.Rmin(3);
    fK.Z(2)  = fK.ZAt(1);
    fK.Rn(2) = fI.Rmin(4);
    fK.Rn(3) = fK.Rmin(2);
    fK.Rn(4) = fK.Rmin(2);
    fK.Rn(5) = fK.Rmin(2);
    fK.Rx(5) = fI.Rmin(8);
    fK.Z(6)  = fI.ZAt(6);
    fK.Rn(6) = fI.Rmin(6);
    fK.Rx(6) = fI.Rmin(7);
    RadiusOfCurvature(fRcurv+fCthick,90.0,fK.ZAt(1),fK.Rmin(1),
		                90.0-fTc,fK.Z(0),fK.Rn(0));
    fK.Rx(0) = fK.Rmin(0);
    fK.Z(3)  = fK.ZAt(0)+(fThickness+2.0*fCthick)*fCostc;;
    fK.Rx(3) = fK.Rmax(0)+(fThickness+2.0*fCthick)*fSintc;
    fK.Rx(1) = RmaxFrom2Points(fK,3,0,fK.ZAt(1));
    fK.Rx(2) = fK.Rmax(1);
    fK.Rx(4) = fI.Rmax(5)-fCthick*fSintc;
    fK.Z(4)  = ZFromRmaxSSDcone(fK.Rmax(4),-fCthick);
    fK.Z(5)  = fI.ZAt(5)-fRcurv*fCostc-fCthick;
    // Now for foam core at the inner most radius.
    fJ.P0() = 0.0;
    fJ.dP() = 360.0;
    fJ.Rn(0) = fI.Rmin(0)-fCthick;
    fJ.Z(0)  = ZFromRminSSDcone(fJ.Rmin(0),+fCthick);
    fJ.Rx(0) = fJ.Rmin(0);
    fJ.Rx(1) = fJ.Rmax(0);
    fJ.Z(1)  = ZFromRmaxSSDcone(fJ.Rmax(1),-fCthick);
    fJ.Rn(1) = RminFromZSSDcone(fJ.ZAt(1),-fCthick);
    fJ.Z(2)  = fK.ZAt(0);
    fJ.Rn(2) = fK.Rmin(0);
    fJ.Rx(2) = RmaxFromZSSDcone(fJ.ZAt(2),-fCthick);
    fJ.Z(3)  = fK.ZAt(3);
    fJ.Rn(3) = fK.Rmax(3);
    fJ.Rx(3) = fJ.Rmin(3);
    // Now for foam core at the top of the inner most radius where 
    // the spoaks are.
    t = fCthick/(0.5*(fRholeMax+fRholeMin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    fL.P0() = 5.0 - 2.0*t; // degrees
    fL.dP() = 12.5+t; // degrees see drawing ALR-0767.
    fL.Z(0) = fH.ZAt(2);
    fL.Rn(0) = fI.Rmin(0);
    fL.Rx(0) = fL.Rmin(0);
    fL.Z(1)  = fJ.ZAt(0);
    fL.Rn(1) = fJ.Rmin(1);
    fL.Rx(1) = fI.Rmax(1);
    fL.Z(2)  = fH.ZAt(3);
    fL.Rn(2) = fL.Rmin(1);
    fL.Rx(2) = fL.Rmax(1);
    fL.Z(3)  = fJ.ZAt(1);
    fL.Rn(3) = fL.Rmin(2);
    fL.Rx(3) = fL.Rmin(3);
    // Now for the SSD mounting posts
    fO.P0()  = 180.0*fdRpost/(fRpostMin+0.5*fdRpost)/TMath::Pi(); // degrees
    fO.dP()  = fPhi0Post; //
    fO.Rn(0) = fRpostMin+fdRpost;
    fO.Rx(0) = fO.Rmin(0);
    fO.Z(0)  = ZFromRmaxSSDcone(fO.Rmax(0));
    fO.Rn(1) = fRpostMin;
    fO.Z(1)  = ZFromRmaxSSDcone(fO.Rmin(1));
    fO.Rx(1) = fO.Rmax(0);
    fO.Z(2)  = fZ0+fZpostMax;
    fO.Rn(2) = fRpostMin;
    fO.Rx(2) = fO.Rmin(2)+fdRpost;
    // Now for the SSD mounting posts
    t = 180.0*fCthick/(fRpostMin+0.5*fdRpost)/TMath::Pi();
    fP.dP()  = fO.DPhi()-2.0*t; // degrees
    fP.P0()  = fO.Phi0()+t; //
    fP.Rn(0) = fO.Rmin(0)-fCthick;
    fP.Rx(0) = fP.Rmin(0);
    fP.Z(0)  = ZFromRmaxSSDcone(fP.Rmax(0));
    fP.Rn(1) = fO.Rmin(0)+fCthick;
    fP.Rx(1) = fO.Rmin(0)-fCthick;
    fP.Z(1)  = ZFromRminSSDcone(fP.Rmin(1));
    fP.Rn(2) = fP.Rmin(1);
    fP.Rx(2) = fP.Rmax(1);
    fP.Z(2)  = fZ0+fZpostMax;
    // This insrto continues into the SSD cone displacing the foam
    // and the carbon fiber surface at those points where the posts are.
    fM.P0()  = fP.Phi0();
    fM.dP()  = fP.DPhi();
    fM.Rn(0) = fRpostMin+fdRpost-fCthick;
    fM.Rx(0) = fM.Rmin(0);
    fM.Z(0)  = ZFromRminSSDcone(fM.Rmin(0),+fCthick);
    fM.Rx(1) = fM.Rmax(0);
    fM.Z(1)  = ZFromRmaxSSDcone(fM.Rmax(1),-fCthick);
    fM.Rn(1) = RminFromZSSDcone(fM.ZAt(1),+fCthick);
    fM.Rn(2) = fRpostMin+fCthick;
    fM.Z(2)  = ZFromRminSSDcone(fM.Rmin(2),+fCthick);
    fM.Rx(2) = RmaxFromZSSDcone(fM.ZAt(2),-fCthick);
    fM.Rn(3) = fM.Rmin(2);
    fM.Rx(3) = fM.Rmin(3);
    fM.Z(3)  = ZFromRmaxSSDcone(fM.Rmax(3),-fCthick);
    //
    fN.P0()  = fP.Phi0();
    fN.dP()  = fP.DPhi();
    fN.Z(0)  = fM.ZAt(1);
    fN.Rn(0) = fM.Rmax(1);
    fN.Rx(0) = fN.Rmin(0);
    fN.Rx(1) = fN.Rmax(0);
    fN.Z(1)  = ZFromRmaxSSDcone(fN.Rmax(1));
    fN.Rn(1) = RminFromZSSDcone(fN.ZAt(1),-fCthick);
    fN.Z(2)  = fM.ZAt(3);
    fN.Rn(2) = fM.Rmin(3);
    fN.Rx(2) = RmaxFromZSSDcone(fN.ZAt(2));
    fN.Rn(3) = fN.Rmin(2);
    fN.Rx(3) = fN.Rmin(3);
    fN.Z(3)  = ZFromRmaxSSDcone(fN.Rmax(3));
    // Bolt heads holding the SSD-SDD tube to the SSD cone.
    // Bolt -- PolyCone
    fQ.P0()  = 0.0;
    fQ.dP()  = 360.0;
    fQ.Z(0)  = fI.ZAt(4)-fThSDDsupportPlate;
    fQ.Rn(0) = 0.0;
    fQ.Rx(0) = 0.5*fDscrewHead;
    fQ.Z(1)  = fI.ZAt(4)-fThScrewHeadHole;
    fQ.Rn(1) = 0.0;
    fQ.Rx(1) = 0.5*fDscrewHead;
    fQ.Z(2)  = fQ.ZAt(1);
    fQ.Rn(2) = 0.0;
    fQ.Rx(2) = 0.5*fDscrewShaft;
    fQ.Z(3)  = fQ.ZAt(2);
    fQ.Rn(3) = 0.0;
    fQ.Rx(3) = fQ.Rmax(2);
    // air infront of bolt (stasolit Volume K) -- Tube
    fR.Z() = 0.5*(fThickness-fThScrewHeadHole);
    fR.Rn() = 0.0;
    fR.Rx() = 0.5*fDscrewHead;
    // air infront of bolt (carbon fiber volume I) -- Tube
    fS.Z() = 0.5*fThickness;
    fS.Rn() = 0.0;
    fS.Rx() = fR.Rmax();
    // SDD support plate, SSD side.
    fT.dP()  = 180.0*fWsddSupportPlate/(fRsddSupportPlate*TMath::Pi());
    fT.P0()  = fPhi0SDDsupports=90.0;
    fT.Z(0)  = fK.ZAt(2);
    fT.Rn(0) = fI.Rmin(4);
    fT.Rx(0) = fRsddSupportPlate;
    fT.Z(1)  = fI.ZAt(4) - fThSDDsupportPlate;
    fT.Rn(1) = fT.Rmin(0);
    fT.Rx(1) = fT.Rmax(0);
    //
    fU.dP() = fT.DPhi();
    fU.P0() = fT.Phi0();
    fU.Z(2) = fI.ZAt(4);
    fU.Rn(2) = fT.Rmin(0);
    fU.Rx(2) = fT.Rmax(0);
    fU.Z(3)  = fT.ZAt(0);
    fU.Rn(3) = fU.Rmin(2);
    fU.Rx(3) = fU.Rmax(2);
    fU.Z(1)  = fU.ZAt(2);
    fU.Rn(1) = fI.Rmin(3);
    fU.Rx(1) = fU.Rmax(3);
    fU.Rx(0) = fT.Rmax(0);
    fU.Rn(0) = fU.Rmax(0);
    fU.Z(0)  = Zfrom2MinPoints(fI,2,3,fU.Rmax(0));
    // Debuging
    Print(&cout);
}
//______________________________________________________________________
void AliITSGeometrySSDCone::CreateG3Geometry(const char *moth,
					     TVector3 &trans){
    // Calls Geant 3 geometry inilization routines with the information
    // stored in this class.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.

    PolyCone(fA,fSSDcf);
    PolyCone(fB,fSSDfs);
    PolyCone(fC,fSSDfo);
    Tube(fD,fSSDsw);
    Tube(fE,fSSDsw);
    PolyCone(fF,fSSDfo);
    PolyCone(fG,fSSDcf);
    PolyCone(fH,fSSDfo);
    PolyCone(fI,fSSDcf);
    PolyCone(fJ,fSSDfo);
    PolyCone(fK,fSSDfs);
    PolyCone(fL,fSSDfo);
    PolyCone(fM,fSSDfo);
    PolyCone(fN,fSSDfo);
    PolyCone(fO,fSSDcf);
    PolyCone(fP,fSSDfs);
    PolyCone(fQ,fSSDfo);
    Tube(fR,fSSDsw);
    Tube(fS,fSSDsw);
    PolyCone(fT,fSSDfo);
    PolyCone(fU,fSSDfo);
    return;
}
//______________________________________________________________________
void AliITSGeometrySSDCone::PositionG3Geometry(AliITSBaseVolParams &moth,
					       Int_t cn,TVector3 &trans,
					       Int_t irot){
    // Positions ths whole object at t with rotatin irot coply number cn
    // into volume moth.
    // Inputs:
    //   const AliITSBaseVolParams *moth   Mother volume where this object 
    //                                     is to be placed.
    //   Int_t      cn      Copy number.
    //   TVector3   &t      Translation vector for this whole volume
    //   Int_t      irot    rotation matrix number to be applyed to this
    //                      volume.
    // Output:
    //   none.
    // Return:
    //   none.
    Int_t i,j,k,l,irotSpoaks,irotPost;
    Double_t t;
    Bool_t init=kFALSE;
    TVector3 zero(0.0,0.0,0.0),v(0.0,0.0,0.0);

    if(cn<=0) return;
    if(cn==1) init=kTRUE;
    Pos(fA,cn,moth,trans,0);
    Pos(fI,cn,moth,trans,0);
    Pos(fG,fNspoaks*(cn-1)+1,fA,trans,0);
    irotSpoaks = irot;
    j = 0;
    for(i=fNspoaks*(cn-1)+2;i<fNspoaks*cn+1;i++){
	ZMatrix(++irot,((Double_t)j)*360./((Double_t)fNspoaks));
	Pos(fG,i,fA,trans,irot);
	j++;
    } // end for i
    Pos(fO,fNposts*(cn-1)+1,moth,trans,0);
    irotPost = irot;
    j = 0;
    for(i=fNposts*(cn-1)+2;i<fNposts*cn+1;i++){
	ZMatrix(++irot,((Double_t)j)*360./((Double_t)fNposts));
	Pos(fO,i,moth,trans,irot);
	j++;
    } // end for
    if(!init) return;
    // Inside volume A.
    Pos(fB,1,fA,zero,0);
    Pos(fC,1,fA,zero,0);
    // Inside Volume B
    k=l=0;
    for(i=0;i<2;i++){ // position for ITS-TPC mounting brackets
	for(j=0;j<2;j++){ // 2 screws per bracket
	    fNcD++;
	    t = -5.0+10.0*((Double_t)j)+180.*((Double_t)i);
	    v.SetX(fRoutHole*Sind(t));
	    v.SetY(fRoutHole*Cosd(t));
	    v.SetZ(fD.DzAt());
	    Pos(fD,fNcD,fB,v,0);
	} // end for j
	for(j=0;j<3;j++){ // 3 pins per bracket
	    fNcE++;
	    t = -3.0+3.0*((Double_t)j)+180.*((Double_t)i);
	    v.SetX(fRoutHole*Sind(t));
	    v.SetY(fRoutHole*Cosd(t));
	    v.SetZ(fE.DzAt());
	    Pos(fE,fNcE,fB,v,0);
	} // end for j
    } // end for i
    for(i=0;i<2;i++){ // position for ITS-rail mounting brackets
	for(j=0;j<4;j++){ // 4 screws per bracket
	    Double_t a[4]={0.0,2.0,5.0,7.0}; // Relative angles.
	    fNcD++;
	    t = 90.0-a[j]+187.*((Double_t)i);
	    v.SetX(fRoutHole*Sind(t));
	    v.SetY(fRoutHole*Cosd(t));
	    v.SetZ(fD.DzAt());
	    Pos(fD,fNcD,fB,v,0);
	} // end for j
	for(j=0;j<2;j++){ // 2 pins per bracket
	    fNcE++;
	    t = 88+7.0*((Double_t)j)+184.*((Double_t)i);
	    v.SetX(fRoutHole*Sind(t));
	    v.SetY(fRoutHole*Cosd(t));
	    v.SetZ(fE.DzAt());
	    Pos(fE,fNcE,fB,v,0);
	} // end for j
    } // end for i
    for(i=0;i<fNmounts;i++){ // mounting holes/screws for beam 
	// pipe support and SPD cone support (dump side,non-dump 
	// side has them to).
	for(j=0;j<2;j++){ // 2 screws per bracket
	    fNcD++;
	    t = 180.*20./(fRoutHole*TMath::Pi());
	    t = 45.0+((Double_t)(j-1))*t+90.*((Double_t)i);
	    v.SetX(fRoutHole*Sind(t));
	    v.SetY(fRoutHole*Cosd(t));
	    v.SetZ(fD.DzAt());
	    Pos(fD,fNcD,fB,v,0);
	} // end for j
	for(j=0;j<1;j++){ // 1 pins per bracket
	    fNcE++;
	    t = 45.0+90.*((Double_t)i);
	    v.SetX(fRoutHole*Sind(t));
	    v.SetY(fRoutHole*Cosd(t));
	    v.SetZ(fE.DzAt());
	    Pos(fE,fNcE,fB,v,0);
	} // end for j
    } // end for i
    Pos(fF,1,fA,zero,0);
    Pos(fL,1,fI,zero,0);
    for(i=1;i<fNspoaks;i++){
	Pos(fF,i+1,fA,zero,irotSpoaks+i);
	Pos(fL,i+1,fA,zero,irotSpoaks+i);
    } // end for i
    Pos(fH,1,fG,zero,0);
    Pos(fK,1,fI,zero,0);
    Pos(fJ,1,fI,zero,0);
    Pos(fP,1,fO,zero,0);
    Pos(fM,1,fJ,zero,0);
    Pos(fN,1,fI,zero,0);
    for(i=1;i<fNposts;i++){
	Pos(fN,i+1,fJ,zero,irotPost+i);
	Pos(fM,i+1,fI,zero,irotPost+i);
    } // end for i
    return;
}
//______________________________________________________________________
void AliITSGeometrySSDCone::CreateG3Materials(){
    // Fills the Geant 3 banks with Material and Medium definisions.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Returns:
    //   none.
    Int_t i;
    Int_t Z[5],N[5];
    Double_t W[5],dens;

    // epoxy
    dens = 10.*GetA(1)+13.*GetA(6)+3.*GetA(8);
    Z[0] = 1; W[0] = 10.*GetA(Z[0])/dens; // Hydrogen Content
    Z[1] = 6; W[1] = 13.*GetA(Z[1])/dens; // Carbon Content
    Z[2] = 8; W[2] =  3.*GetA(Z[2])/dens; // Oxegen
    // Carbon fiber is about 64% carbon fiber and 36% epoxy by volume.
    // Now need to add in the carbon fiber
    W[0] *= 0.36*GetA(Z[0]);
    W[1]  = 0.36*W[1]*dens*GetA(Z[1]) + 0.64*GetA(Z[1]);
    W[2] *= 0.36*GetA(Z[2]);
    // Renormilize the weights
    dens = 0.0;
    for(i=0;i<3;i++){dens += W[i];}
    for(i=0;i<3;i++){W[i] /= dens;}
    dens = 1.7; // grams/cm^3 taken as density of G10 PDG book.
    MixtureByWeight(fSSDcf,"Carbon Fiber for SSD support cone",Z,W,dens,3,0);
    // epoxy 
    dens = 10.*GetA(1)+13.*GetA(6)+3.*GetA(8);
    Z[0] = 1; W[0] = 10.*GetA(Z[0])/dens; // Hydrogen Content
    Z[1] = 6; W[1] = 13.*GetA(Z[1])/dens; // Carbon Content
    Z[2] = 8; W[2] =  3.*GetA(Z[2])/dens; // Oxegen
    Z[3] = 14;W[3] = 0.0; // no Silicon in epoxy.
    // glass fiber is about 64% carbon fiber and 36% epoxy by volume.
    // Now need to add in the glass fiber
    W[0] *= 0.36*GetA(Z[0]);
    W[1] *= 0.36*GetA(Z[1]);
    W[2]  = 0.36*W[2]*dens*GetA(Z[2]) + 0.64*2.0*GetA(Z[2]);
    W[3]  = 0.64*GetA(Z[3]); // Si
    // Renormilize the weights
    dens = 0.0;
    for(i=0;i<4;i++){dens += W[i];}
    for(i=0;i<4;i++){W[i] /= dens;}
    dens = 1.7; // grams/cm^3 taken as density of G10 PDG book.
    MixtureByWeight(fSSDfs,"Inserto stealite 4411w for SSD support cone",
		    Z,W,dens,4,0);
    // Rohacell 51 C14 H10 N2 O6 from Flavio Tosello
    // http://cesweb.grantadesign.com/demo/index.do
    Z[0] = 1; N[0] = 10; // Hydrogen Content
    Z[1] = 6; N[1] = 14; // Carbon Content
    Z[2] = 7; N[2] =  2; // Nitrogen Content
    Z[3] = 8; N[3] =  6; // Oxigen Content
    dens = 0.0513; // grams/cm^3 From Flavio Tosello 
    //  http://www.emkayplatics.co.uk/roh51.html
    MixtureByNumber(fSSDfo,"Foam core (Rohacell 51) for SSD support cone",
		    Z,N,dens,4,0);
    // Stainless steel. Temperary values.
    Z[0] =  6; W[0] = 0.5; // Carbon Content
    Z[1] = 25; W[1] = 0.5; // Iron Content
    dens = 7.87; // Grams/cm^3  density of iron used.
    MixtureByWeight(fSSDsw,"Stainless steal screw, pin, and stud material",
		    Z,W,dens,2,0);
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
//______________________________________________________________________
void AliITSGeometrySSDCone::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //   ostream *os pointer to the output stream
    // Outputs:
    //   none.
    // Return:
    //   none.

    *os << "Object AliITSGeometrySSDCone" << endl;
    *os << " Object fA" << endl << fA << endl;
    *os << " Object fB" << endl << fB << endl;
    *os << " Object fC" << endl << fC << endl;
    *os << " Object fD" << endl << fD << endl;
    *os << " Object fE" << endl << fE << endl;
    *os << " Object fF" << endl << fF << endl;
    *os << " Object fG" << endl << fG << endl;
    *os << " Object fH" << endl << fH << endl;
    *os << " Object fI" << endl << fI << endl;
    *os << " Object fJ" << endl << fJ << endl;
    *os << " Object fK" << endl << fK << endl;
    *os << " Object fL" << endl << fL << endl;
    *os << " Object fM" << endl << fM << endl;
    *os << " Object fN" << endl << fN << endl;
    *os << " Object fO" << endl << fO << endl;
    *os << " Object fP" << endl << fP << endl;
    *os << " Object fQ" << endl << fQ << endl;
    *os << " Object fR" << endl << fR << endl;
    *os << " Object fS" << endl << fS << endl;
    *os << " Object fT" << endl << fT << endl;
    *os << " Object fU" << endl << fU << endl;
    return;
}
//______________________________________________________________________
void AliITSGeometrySSDCone::Read(istream *is){
    // Read in data written with Print above.
    // Inputs:
    //   istream *is Input stream pointer
    // Output:
    //   none.
    // Return:
    //   none.
    char s[50];

    is->getline(s,49);
    is->getline(s,49);
    *is >> fA;
    is->getline(s,49);
    *is >> fB;
    is->getline(s,49);
    *is >> fC;
    is->getline(s,49);
    *is >> fD;
    is->getline(s,49);
    *is >> fE;
    is->getline(s,49);
    *is >> fF;
    is->getline(s,49);
    *is >> fG;
    is->getline(s,49);
    *is >> fH;
    is->getline(s,49);
    *is >> fI;
    is->getline(s,49);
    *is >> fJ;
    is->getline(s,49);
    *is >> fK;
    is->getline(s,49);
    *is >> fL;
    is->getline(s,49);
    *is >> fM;
    is->getline(s,49);
    *is >> fN;
    is->getline(s,49);
    *is >> fO;
    is->getline(s,49);
    *is >> fP;
    is->getline(s,49);
    *is >> fQ;
    is->getline(s,49);
    *is >> fR;
    is->getline(s,49);
    *is >> fS;
    is->getline(s,49);
    *is >> fT;
    is->getline(s,49);
    *is >> fU;
    return;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSGeometrySSDCone &s){
    // Operator << for C++ like output of AliITSGeometrySSDCone class.
    // Inputs:
    //   ostream &os              The output stream
    //   AliITSGeometrySSDCone &s The class to be outputed
    // Outputs:
    //   none.
    // Return:
    //  ostream &os  The address of the output stream

    s.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSGeometrySSDCone &s){
    // Operator >> for C++ like input of AliITSGeometrySSDCone class.
    // Inputs:
    //   istream &is              The input stream
    //   AliITSGeometrySSDCone &s The class to be inputed
    // Outputs:
    //   none.
    // Return:
    //  istream &is  The address of the input stream

    s.Read(&is);
    return is;
}
