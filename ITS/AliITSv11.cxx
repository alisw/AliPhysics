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

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  Inner Traking System version 11                                         //
//  This class contains the base procedures for the Inner Tracking System   //
//                                                                          //
// Authors: R. Barbera                                                      //
// version 6.                                                               //
// Created  2000.                                                           //
//                                                                          //
//  NOTE: THIS IS THE  SYMMETRIC PPR geometry of the ITS.                   //
// THIS WILL NOT WORK                                                       //
// with the geometry or module classes or any analysis classes. You are     //
// strongly encouraged to uses AliITSv5.                                    //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// See AliITSv11::StepManager().
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


#include "AliRun.h"
#include "AliMagF.h"
#include "AliConst.h"
#include "AliITSGeant3Geometry.h"
#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSv11.h"
#include "AliITSgeom.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "AliITSDetType.h"
#include "AliITSresponseSPD.h"
#include "AliITSresponseSDD.h"
#include "AliITSresponseSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSsimulationSPD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSsimulationSSD.h"
#include "AliITSClusterFinderSPD.h"
#include "AliITSClusterFinderSDD.h"
#include "AliITSClusterFinderSSD.h"


ClassImp(AliITSv11)

//______________________________________________________________________
AliITSv11::AliITSv11() : AliITS() {
////////////////////////////////////////////////////////////////////////
//    Standard default constructor for the ITS version 11.
////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
AliITSv11::AliITSv11(const char *title) : AliITS("ITS", title){
////////////////////////////////////////////////////////////////////////
//    Standard constructor for the ITS version 11.
////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
AliITSv11::~AliITSv11() {
////////////////////////////////////////////////////////////////////////
//    Standard destructor for the ITS version 11.
////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
void AliITSv11::Box(const char gnam[3],const TString &dis,
		    Double_t dx,Double_t dy,Double_t dz,Int_t med){
    // Interface to TMC->Gsvolu() for ITS bos geometries. Box with faces
    // perpendicular to the axes. It has 3 paramters. See SetScale() for
    // units. Default units are geant 3 [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t dx         half-length of box in x-axis
    //    Double_t dy         half-length of box in y-axis
    //    Double_t dz         half-length of box in z-axis
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[3];

    param[0] = fScale*dx;
    param[1] = fScale*dy;
    param[2] = fScale*dz;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"BOX ",fidmed[med],param,3);
}
//______________________________________________________________________
void AliITSv11::Trapezoid1(const char gnam[3],const TString &dis,
			   Double_t dxn,Double_t dxp,Double_t dy,Double_t dz,
			   Int_t med){
    // Interface to TMC->Gsvolu() for ITS TRD1 geometries. Trapezoid with the 
    // x dimension varing along z. It has 4 parameters. See SetScale() for
    // units. Default units are geant 3 [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t dxn        half-length along x at the z surface positioned 
    //                        at -DZ
    //    Double_t dxp        half-length along x at the z surface positioned 
    //                        at +DZ
    //    Double_t dy         half-length along the y-axis
    //    Double_t dz         half-length along the z-axis
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[4];

    param[0] = fScale*dxn;
    param[1] = fScale*dxp;
    param[2] = fScale*dy;
    param[3] = fScale*dz;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TRD1",fidmed[med],param,4);
}
//______________________________________________________________________
void AliITSv11::Trapezoid2(const char gnam[3],const TString &dis,Double_t dxn,
			   Double_t dxp,Double_t dyn,Double_t dyp,Double_t dz,
			   Int_t med){
    // Interface to TMC->Gsvolu() for ITS TRD2 geometries. Trapezoid with the 
    // x and y dimension varing along z. It has 5 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t dxn        half-length along x at the z surface positioned 
    //                        at -DZ
    //    Double_t dxp        half-length along x at the z surface positioned 
    //                        at +DZ
    //    Double_t dyn        half-length along x at the z surface positioned 
    //                        at -DZ
    //    Double_t dyp        half-length along x at the z surface positioned 
    //                        at +DZ
    //    Double_t dz         half-length along the z-axis
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[5];

    param[0] = fScale*dxn;
    param[1] = fScale*dxp;
    param[2] = fScale*dyn;
    param[3] = fScale*dyp;
    param[4] = fScale*dz;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TRD2",fidmed[med],param,5);
}
//______________________________________________________________________
void AliITSv11::Trapezoid(const char gnam[3],const TString &dis,Double_t dz,
			  Double_t thet,Double_t phi,Double_t h1,Double_t bl1,
			  Double_t tl1,Double_t alp1,Double_t h2,Double_t bl2,
			  Double_t tl2,Double_t alp2,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TRAP geometries. General Trapezoid, 
    // The faces perpendicular to z are trapezia and their centers are not 
    // necessarily on a line parallel to the z axis. This shape has 11 
    // parameters, but only cosidering that the faces should be planar, only 9 
    // are really independent. A check is performed on the user parameters and 
    // a message is printed in case of non-planar faces. Ignoring this warning 
    // may cause unpredictable effects at tracking time. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t dz         Half-length along the z-asix
    //    Double_t thet       Polar angle of the line joing the center of the 
    //                        face at -dz to the center of the one at dz 
    //                        [degree].
    //    Double_t phi        aximuthal angle of the line joing the center of 
    //                        the face at -dz to the center of the one at +dz 
    //                        [degree].
    //    Double_t h1         half-length along y of the face at -dz.
    //    Double_t bl1        half-length along x of the side at -h1 in y of 
    //                        the face at -dz in z.
    //    Double_t tl1        half-length along x of teh side at +h1 in y of 
    //                        the face at -dz in z.
    //    Double_t alp1       angle with respect to the y axis from the center 
    //                        of the side at -h1 in y to the cetner of the 
    //                        side at +h1 in y of the face at -dz in z 
    //                        [degree].
    //    Double_t h2         half-length along y of the face at +dz
    //    Double_t bl2        half-length along x of the side at -h2 in y of
    //                        the face at +dz in z.
    //    Double_t tl2        half-length along x of the side at _h2 in y of 
    //                        the face at +dz in z.
    //    Double_t alp2       angle with respect to the y axis from the center 
    //                        of the side at -h2 in y to the center of the 
    //                        side at +h2 in y of the face at +dz in z 
    //                        [degree].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[11];

    param[0] = fScale*dz;
    param[1] = thet;
    param[2] = phi;
    param[3] = fScale*h1;
    param[4] = fScale*bl1;
    param[5] = fScale*tl1;
    param[6] = alp1;
    param[7] = fScale*h2;
    param[8] = fScale*bl2;
    param[9] = fScale*tl2;
    param[10] = alp2;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TRAP",fidmed[med],param,11);
}
//______________________________________________________________________
void AliITSv11::Tube(const char gnam[3],const TString &dis,Double_t rmin,
		     Double_t rmax,Double_t dz,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TUBE geometries. Simple Tube. It has
    // 3 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t rmin       Inside Radius.
    //    Double_t rmax       Outside Radius.
    //    Double_t dz         half-length along the z-axis
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[3];

    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TUBE",fidmed[med],param,3);
}
//______________________________________________________________________
void AliITSv11::TubeSegment(const char gnam[3],const TString &dis,
			    Double_t rmin,Double_t rmax,Double_t dz,
			    Double_t phi1,Double_t phi2,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TUBE geometries. Phi segment of a 
    // tube. It has 5  parameters. Phi1 should be smaller than phi2. If this is
    // not the case, the system adds 360 degrees to phi2. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t rmin       Inside Radius.
    //    Double_t rmax       Outside Radius.
    //    Double_t dz         half-length along the z-axis
    //    Double_t phi1       Starting angle of the segment [degree].
    //    Double_t phi2       Ending angle of the segment [degree].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[5];

    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    param[3] = phi1;
    param[4] = phi2;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TUBS",fidmed[med],param,5);
}
//______________________________________________________________________
void AliITSv11::Cone(const char gnam[3],const TString &dis,Double_t dz,
		     Double_t rmin1,Double_t rmax1,Double_t rmin2,
		     Double_t rmax2,Int_t med){
    // Interface to TMC->Gsvolu() for ITS Cone geometries. Conical tube. It 
    // has 5 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t dz         half-length along the z-axis
    //    Double_t rmin1      Inside Radius at -dz.
    //    Double_t rmax1      Outside Radius at -dz.
    //    Double_t rmin2      inside radius at +dz.
    //    Double_t rmax2      outside radius at +dz.
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[5];

    param[0] = fScale*dz;
    param[1] = fScale*rmin1;
    param[2] = fScale*rmax1;
    param[3] = fScale*rmin2;
    param[4] = fScale*rmax2;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"CONS",fidmed[med],param,5);
}
//______________________________________________________________________
void AliITSv11::ConeSegment(const char gnam[3],const TString &dis,Double_t dz,
			    Double_t rmin1,Double_t rmax1,Double_t rmin2,
			    Double_t rmax2,Double_t phi1,Double_t phi2,
			    Int_t med){
    // Interface to TMC->Gsvolu() for ITS ConS geometries. One segment of a 
    // conical tube. It has 7 parameters. Phi1 should be smaller than phi2. If 
    // this is not the case, the system adds 360 degrees to phi2. See 
    // SetScale() for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t dz         half-length along the z-axis
    //    Double_t rmin1      Inside Radius at -dz.
    //    Double_t rmax1      Outside Radius at -dz.
    //    Double_t rmin2      inside radius at +dz.
    //    Double_t rmax2      outside radius at +dz.
    //    Double_t phi1       Starting angle of the segment [degree].
    //    Double_t phi2       Ending angle of the segment [degree].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[7];

    param[0] = fScale*dz;
    param[1] = fScale*rmin1;
    param[2] = fScale*rmax1;
    param[3] = fScale*rmin2;
    param[4] = fScale*rmax2;
    param[5] = phi1;
    param[6] = phi2;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"CONS",fidmed[med],param,7);
}
//______________________________________________________________________
void AliITSv11::Sphere(const char gnam[3],const TString &dis,Double_t rmin,
		       Double_t rmax,Double_t the1,Double_t the2,Double_t phi1,
		       Double_t phi2,Int_t med){
    // Interface to TMC->Gsvolu() for ITS SPHE geometries. Segment of a 
    // sphereical shell. It has 6 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t rmin       Inside Radius.
    //    Double_t rmax       Outside Radius.
    //    Double_t the1       staring polar angle of the shell [degree].
    //    Double_t the2       ending polar angle of the shell [degree].
    //    Double_t phui       staring asimuthal angle of the shell [degree].
    //    Double_t phi2       ending asimuthal angle of the shell [degree].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[6];

    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = the1;
    param[3] = the2;
    param[4] = phi1;
    param[5] = phi2;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"SPHE",fidmed[med],param,6);
}
//______________________________________________________________________
void AliITSv11::Parallelepiped(const char gnam[3],const TString &dis,
			       Double_t dx,Double_t dy,Double_t dz,
			       Double_t alph,Double_t thet,Double_t phi,
			       Int_t med){
    // Interface to TMC->Gsvolu() for ITS PARA geometries. Parallelepiped. It 
    // has 6 parameters. See SetScale() for units. Default units are geant 3 
    // [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t dx         half-length allong x-axis
    //    Double_t dy         half-length allong y-axis
    //    Double_t dz         half-length allong z-axis
    //    Double_t alpha      angle formed by the y axis and by the plane 
    //                        joining the center of teh faces parallel to the 
    //                        z-x plane at -dY and +dy [degree].
    //    Double_t thet       polar angle of the line joining the centers of 
    //                        the faces at -dz and +dz in z [degree].
    //    Double_t phi        azimuthal angle of teh line joing the centers of 
    //                        the faaces at -dz and +dz in z [degree].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[6];

    param[0] = fScale*dx;
    param[1] = fScale*dy;
    param[2] = fScale*dz;
    param[3] = alpha;
    param[4] = thet;
    param[5] = phi;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"PARA",fidmed[med],param,6);
}
//______________________________________________________________________
void AliITSv11::Polygon(const char gnam[3],const TString &dis,Double_t phi1,
			Double_t dphi,Int_t npdv,Int_t nz,Double_t *z,
			Double_t *rmin,Double_t *rmax,Double_t ,Int_t med){
    // Interface to TMC->Gsvolu() for ITS PGON geometry. Polygon It has 10 
    // parameters or more. See SetScale() for units. Default units are geant 3 
    // [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t phi1       the azimuthal angle at which the volume begins 
    //                        (angles are counted clouterclockwise) [degrees].
    //    Double_t dphi       opening angle of the volume, which extends from 
    //                        phi1 to phi1+dphi [degree].
    //    Int_t npdv          the number of sides of teh cross section between 
    //                        the given phi limits.
    //    Int_t nz            number of planes perpendicular to the z axis 
    //                        where the dimension of the section is given - 
    //                        this number should be at least 2 and NP triples 
    //                        of number must follow.
    //    Double_t *z         array [nz] of z coordiates of the sections..
    //    Double_t *rmin      array [nz] of radius of teh circle tangent to 
    //                        the sides of the inner polygon in teh 
    //                        cross-section.
    //    Double_t *rmax      array [nz] of radius of the circle tangent to 
    //                        the sides of the outer polygon in the 
    //                       cross-section.
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t *param;
    Int_t n,i;

    n = 4+3*nz;
    param = new Float_t[n]
    param[0] = phi1;
    param[1] = dphi;
    param[2] = (Float_t)npdv;
    param[3] = (Float_t)nz;
    for(i=0;i<nz;i++){
	param[4+3*i] = z[i];
	param[5+3*i] = rmin[i];
	param[6+3*i] = rmax[i];
    } // end for i
    name[0] = 'I';
    for(i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"PGON",fidmed[med],param,n);

    delete[] param;
}
//______________________________________________________________________
void AliITSv11::PolyCone(const char gnam[3],const TString &dis,Double_t phi1,
			 Double_t dphi,Int_t nz,Double_t *z,Double_t *rmin,
			 Double_t *rmax,Int_t med){
    // Interface to TMC->Gsvolu() for ITS PCON geometry. Poly-cone It has 9 
    // parameters or more. See SetScale() for units. Default units are geant 3 
    // [cm].
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t phi1       the azimuthal angle at which the volume begins 
    //                        (angles are counted clouterclockwise) [degrees].
    //    Double_t dphi       opening angle of the volume, which extends from 
    //                        phi1 to phi1+dphi [degree].
    //    Int_t nz            number of planes perpendicular to the z axis 
    //                        where the dimension of the section is given - 
    //                        this number should be at least 2 and NP triples 
    //                        of number must follow.
    //    Double_t *z         Array [nz] of z coordinate of the section.
    //    Double_t *rmin      Array [nz] of radius of teh inner circle in the 
    //                        cross-section.
    //    Double_t *rmax      Array [nz] of radius of the outer circle in the 
    //                        cross-section.
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t *param;
    Int_t n,i;

    n = 3+3*nz;
    param = new Float_t[n];
    param[0] = phi1;
    param[1] = dphi;
    param[2] = (Float_t) nz;
    for(i=0;i<nz;i++){
	param[3+3*i] = z[i];
	param[4+3*i] = rmin[i];
	param[5+3*i] = rmax[i];
    } // end for i
    name[0] = 'I';
    for(i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"PCON",fidmed[med],param,n);

    delete[] param;
}
//______________________________________________________________________
void AliITSv11::TubeElliptical(const char gnam[3],const TString &dis,
			       Double_t p1,Double_t p2,Double_t dz,Int_t med){
    // Interface to TMC->Gsvolu() for ITS ELTU geometries. Elliptical 
    // cross-section Tube. It has 3 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm]. The equation of the surface 
    // is x^2 * p1^-2 + y^2 * p2^-2 = 1.
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t p1         semi-axis of the elipse along x.
    //    Double_t p2         semi-axis of the elipse along y.
    //    Double_t dz         half-length along the z-axis
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[3];

    param[0] = fScale*p1;
    param[1] = fScale*p2;
    param[2] = fScale*dz;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"ELTU",fidmed[med],param,3);
}
//______________________________________________________________________
void AliITSv11::HyperbolicTube(const char gnam[3],const TString &dis,
			       Double_t rmin,Double_t rmax,Double_t dz,
			       Double_t thet,Int_t med){
    // Interface to TMC->Gsvolu() for ITS HYPE geometries. Hyperbolic tube. 
    // Fore example the inner and outer surfaces are hyperboloids, as would be 
    // foumed by a system of cylinderical wires which were then rotated 
    // tangentially about their centers. It has 4 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm]. The hyperbolic surfaces are 
    // given by r^2 = (ztan(thet)^2 + r(z=0)^2.
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t rmin       Inner radius at z=0 where tube is narrowest.
    //    Double_t rmax       Outer radius at z=0 where tube is narrowest.
    //    Double_t dz         half-length along the z-axis
    //    Double_t thet       stero angel of rotation of the two faces 
    //                       [degrees].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[4];

    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    param[3] = thet;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"HYPE",fidmed[med],param,4);
}
//______________________________________________________________________
void AliITSv11::TwistedTrapezoid(const char gnam[3],const TString &dis,
				 Double_t dz,Double_t thet,Double_t phi,
				 Double_t twist,Double_t h1,Double_t bl1,
				 Double_t tl1,Double_t apl1,Double_t h2,
				 Double_t bl2,Double_t tl2,Double_t apl2,
				 Int_t med){
    // Interface to TMC->Gsvolu() for ITS GTRA geometries. General twisted 
    // trapazoid. The faces perpendicular to z are trapazia and their centers 
    // are not necessarily on a line parallel to the z axis as the TRAP. 
    // Additionally, the faces may be twisted so that none of their edges are 
    // parallel. It is a TRAP shape, exept that it is twisted in the x-y plane 
    // as a function of z. The parallel sides perpendicular to the x axis are 
    // rotated with respect to the x axis by an angle TWIST, which is one of 
    // the parameters. The shape is defined by the eight corners and is assumed
    // to be constructed of straight lines joingin points on the boundry of the
    // trapezoidal face at Z=-dz to the coresponding points on the face at 
    // z=+dz. Divisions are not allowed. It has 12 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm]. Note: This shape suffers from
    // the same limitations than the TRAP. The tracking routines assume that 
    // the faces are planar, but htis constraint is not easily expressed in 
    // terms of the 12 parameters. Additionally, no check on th efaces is 
    // performed in this case. Users should avoid to use this shape as much as 
    // possible, and if they have to do so, they should make sure that the 
    // faces are really planes. If this is not the case, the result of the 
    // trasport is unpredictable. To accelerat ethe computations necessary for 
    // trasport, 18 additioanl parameters are calculated for this shape are
    // 1 DXODZ dx/dz of the line joing the centers of the faces at z=+_dz.
    // 2 DYODZ dy/dz of the line joing the centers of the faces at z=+_dz.
    // 3 XO1    x at z=0 for line joing the + on parallel side, perpendicular 
    //          corners at z=+_dz.
    // 4 YO1    y at z=0 for line joing the + on parallel side, + on 
    //          perpendicular corners at z=+-dz.
    // 5 DXDZ1  dx/dz for line joing the + on parallel side, + on 
    //          perpendicular corners at z=+-dz.
    // 6 DYDZ1  dy/dz for line joing the + on parallel side, + on 
    //          perpendicular corners at z=+-dz.
    // 7 X02    x at z=0 for line joing the - on parallel side, + on 
    //          perpendicular corners at z=+-dz.
    // 8 YO2    y at z=0 for line joing the - on parallel side, + on 
    //          perpendicular corners at z=+-dz.
    // 9 DXDZ2  dx/dz for line joing the - on parallel side, + on 
    //          perpendicular corners at z=+-dz.
    // 10 DYDZ2dy/dz for line joing the - on parallel side, + on 
    //          perpendicular corners at z=+-dz.
    // 11 XO3   x at z=0 for line joing the - on parallel side, - on 
    //          perpendicular corners at z=+-dz.
    // 12 YO3   y at z=0 for line joing the - on parallel side, - on 
    //          perpendicular corners at z=+-dz.
    // 13 DXDZ3 dx/dzfor line joing the - on parallel side, - on 
    //          perpendicular corners at z=+-dz.
    // 14 DYDZ3 dydz for line joing the - on parallel side, - on 
    //          perpendicular corners at z=+-dz.
    // 15 XO4   x at z=0 for line joing the + on parallel side, - on 
    //          perpendicular corners at z=+-dz.
    // 16 YO4   y at z=0 for line joing the + on parallel side, - on 
    //          perpendicular corners at z=+-dz.
    // 17 DXDZ4 dx/dz for line joing the + on parallel side, - on 
    //          perpendicular corners at z=+-dz.
    // 18 DYDZ4 dydz for line joing the + on parallel side, - on 
    //          perpendicular corners at z=+-dz.
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t dz         half-length along the z axis.
    //    Double_t thet       polar angle of the line joing the center of the 
    //                        face at -dz to the center of the one at +dz 
    //                        [degrees].
    //    Double_t phi        Azymuthal angle of teh line joing the centre of 
    //                        the face at -dz to the center of the one at +dz 
    //                        [degrees].
    //    Double_t twist      Twist angle of the faces parallel to the x-y 
    //                        plane at z=+-dz around an axis parallel to z 
    //                        passing through their centre [degrees].
    //    Double_t h1         Half-length along y of the face at -dz.
    //    Double_t bl1        half-length along x of the side -h1 in y of the 
    //                        face at -dz in z.
    //    Double_t tl1        half-length along x of the side at +h1 in y of 
    //                        the face at -dz in z.
    //    Double_t apl1       Angle with respect to the y ais from the center 
    //                        of the side at -h1 in y to the centere of the 
    //                        side at +h1 in y of the face at -dz in z 
    //                        [degrees].
    //    Double_t h2         half-length along the face at +dz.
    //    Double_t bl2        half-length along x of the side at -h2 in y of 
    //                        the face at -dz in z.
    //    Double_t tl2        half-length along x of the side at +h2 in y of 
    //                        the face at +dz in z.
    //    Double_t apl2       angle with respect to the y axis from the center 
    //                        of the side at -h2 in y to the center of the side
    //                        at +h2 in y of the face at +dz in z [degrees].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[12];

    param[0] = fScale*dz;
    param[1] = thet;
    param[2] = phi;
    param[3] = twist;
    param[4] = fScale*h1;
    param[5] = fScale*bl1;
    param[6] = fScale*tl1;
    param[7] = alp1;
    param[8] = fScale*h2;
    param[9] = fScale*bl2;
    param[10] = fScale*tl2;
    param[11] = alp2;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"GTRA",fidmed[med],param,12);
}
//______________________________________________________________________
void AliITSv11::CutTube(const char gnam[3],const TString &dis,Double_t rmin,
			Double_t rmax,Double_t dz,Double_t phi1,Double_t phi2,
			Double_t lx,Double_t ly,Double_t lz,Double_t hx,
			Double_t hy,Double_t hz,Int_t med){
    // Interface to TMC->Gsvolu() for ITS CTUB geometries. Cut tube. A tube cut
    // at the extremities with planes not necessarily perpendicular tot he z 
    // axis. It has 11 parameters. See SetScale() for units. Default units are 
    // geant 3 [cm]. phi1 should be smaller than phi2. If this is not the case,
    // the system adds 360 degrees to phi2.
    // Inputs:
    //    const char gnam[3]  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t rmin       Inner radius at z=0 where tube is narrowest.
    //    Double_t rmax       Outer radius at z=0 where tube is narrowest.
    //    Double_t dz         half-length along the z-axis
    //    Double_t dz         half-length along the z-axis
    //    Double_t phi1       Starting angle of the segment [degree].
    //    Double_t phi2       Ending angle of the segment [degree].
    //    Double_t lx         x component of a unit vector perpendicular to 
    //                        the face at -dz.
    //    Double_t ly         y component of a unit vector perpendicular to 
    //                        the face at -dz.
    //    Double_t lz         z component of a unit vector perpendicular to 
    //                        the face at -dz.
    //    Double_t hx         x component of a unit vector perpendicular to 
    //                        the face at +dz.
    //    Double_t hy         y component of a unit vector perpendicular to 
    //                        the face at +dz.
    //    Double_t hz         z component of a unit vector perpendicular to 
    //                        the face at +dz.
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[4];
    Float_t param[11];

    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    param[3] = phi1;
    param[4] = phi2;
    param[5] = lx;
    param[6] = ly;
    param[7] = lz;
    param[8] = hx;
    param[9] = hy;
    param[10] = hz;
    name[0] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"CTUB",fidmed[med],param,11);
}
//______________________________________________________________________
void AliITSv11::Pos(const char vol[3],Int_t cn,const char moth[3],Double_t x,
		    Double_t y,Double_t z,Int_t irot){
    // Place a copy of a volume previously defined by a call to GSVOLU inside 
    // its mother volulme moth.
    // Inputs:
    //   const char vol[3]  3 character geant volume name. The letter "I"
    //                      is appended to the front to indecate that this
    //                      is an ITS volume.
    //   const char moth[3] 3 character geant volume name of the mother volume 
    //                      in which vol will be placed. The letter "I" is 
    //                      appended to the front to indecate that this is an 
    //                      ITS volume.
    //   Double_t x         The x positon of the volume in the mother's 
    //                      reference system
    //   Double_t y         The y positon of the volume in the mother's 
    //                      reference system
    //   Double_t z         The z positon of the volume in the mother's 
    //                      reference system
    //   Int_t irot         the index for the rotation matrix to be used.
    //                      irot=-1 => unit rotation.
    // Outputs:
    //    none.
    // Return:
    //    none.
    char name[4],mother[4];
    Float_t param[3];
    Int_t r=0,i;

    param[0] = x;
    param[1] = y;
    param[2] = z;
    name[0] = 'I';
    for(i=0;i<3;i++) name[i+1] = vol[i];
    mother[0] = 'I';
    for(i=0;i<3;i++) mother[i+1] = moth[i];
    if(irot>=0) r=fidrot[irot];
    fMC->Gspos(name,mother,param[0],param[1],param[2],r,"ONLY");
}
//______________________________________________________________________
void AliITSv11::Matrix(Int_t irot,Double_t thet1,Double_t phi1,
		       Double_t thet2,Double_t phi2,
		       Double_t thet3,Double_t phi3){
    // Defines a Geant rotation matrix. checks to see if it is the unit
    // matrix. If so, then no additonal matrix is defined. Stores rotation 
    // matrix irot in the data structure JROTM. If the matrix is not 
    // orthonormal, it will be corrected by setting y' perpendicular to x' 
    // and z' = x' X y'. A warning message is printed in this case.
    // Inputs:
    //   Int_t irot     Intex specifing which rotation matrix.
    //   Double_t thet1 Polar angle for axisw x [degrees].
    //   Double_t phi1  azimuthal angle for axis x [degrees].
    //   Double_t thet12Polar angle for axisw y [degrees].
    //   Double_t phi2  azimuthal angle for axis y [degrees].
    //   Double_t thet3 Polar angle for axisw z [degrees].
    //   Double_t phi3  azimuthal angle for axis z [degrees].
    // Outputs:
    //    none.
    // Return:
    //    none.
    Float_t t1,p1,t2,p2,t3,p3;

    if(thet1==90.0&&phi1==0.0&&thet2==90.0&&phi2==90.0&&thet3==0.0&&phi3==0.0){
	fidrot[irot] = 0; // Unit matrix
    }else{
	t1 = thet1;
	p1 = phi1;
	t2 = thet2;
	p2 = phi2;
	t3 = thet3;
	p3 = phi3
	AliMatrix(fidrot[irot],t1,p1,t2,p2,t3,p3);
    } // end if
}
//______________________________________________________________________
void AliITSv11::Matrix(Int_t irot,Double_t rot[3][3]){
    // Defines a Geant rotation matrix. checks to see if it is the unit
    // matrix. If so, then no additonal matrix is defined. Stores rotation 
    // matrix irot in the data structure JROTM. If the matrix is not 
    // orthonormal, it will be corrected by setting y' perpendicular to x' 
    // and z' = x' X y'. A warning message is printed in this case.
    // Inputs:
    //   Int_t irot         Intex specifing which rotation matrix.
    //   Double_t rot[3][3] The 3 by 3 rotation matrix.
    // Outputs:
    //    none.
    // Return:
    //    none.

    if(rot[0][0]==1.0&&rot[1][1]==1.0&&rot[2][2]==1.0&&
       rot[0][1]==0.0&&rot[0][2]==0.0&&rot[1][0]==0.0&&
       rot[1][2]==0.0&&rot[2][0]==0.0&&rot[2][1]==0.0){
	fidrot[irot] = 0; // Unit matrix
    }else{
	Double_t si,c=180./TMath::Pi();
	Double_t ang[6];

	ang[1] = TMath::ATan2(rot[0][1],rot[0][0]);
	if(TMath::Cos(ang[1])!=0.0) si = rot[0][0]/TMath::Cos(ang[1]);
	else si = rot[0][1]/TMath::Sin(ang[1]);
	ang[0] = TMath::ATan2(si,rot[0][2]);

	ang[3] = TMath::ATan2(rot[1][1],rot[1][0]);
	if(TMath::Cos(ang[3])!=0.0) si = rot[1][0]/TMath::Cos(ang[3]);
	else si = rot[1][1]/TMath::Sin(ang[3]);
	ang[2] = TMath::ATan2(si,rot[1][2]);

	ang[5] = TMath::ATan2(rot[2][1],rot[2][0]);
	if(TMath::Cos(ang[5])!=0.0) si = rot[2][0]/TMath::Cos(ang[5]);
	else si = rot[2][1]/TMath::Sin(ang[5]);
	ang[4] = TMath::ATan2(si,rot[2][2]);

	for(Int_t i=0;i<6;i++) {ang[i] *= c; if(ang[i]<0.0) ang[i] += 360.;}
	AliMatrix(fidrot[irot],ang[0],ang[1],ang[2],ang[3],ang[4],ang[5]);
    } // end if
}
//______________________________________________________________________
void AliITSv11::Matrix(Int_t irot,Int_t axis,Double_t thet){
    // Defines a Geant rotation matrix. checks to see if it is the unit
    // matrix. If so, then no additonal matrix is defined. Stores rotation 
    // matrix irot in the data structure JROTM. If the matrix is not 
    // orthonormal, it will be corrected by setting y' perpendicular to x' 
    // and z' = x' X y'. A warning message is printed in this case.
    // Inputs:
    //   Int_t irot         Intex specifing which rotation matrix.
    //   Int_t axis         Axis about which rotation is to be done.
    //   Double_t thet      Angle to rotate by [degrees].
    // Outputs:
    //    none.
    // Return:
    //    none.

    if(thet==0.0){
	fidrot[irot] = 0; // Unit matrix
    }else{
	switch (irot) {
	case 0: //Rotate about x-axis, x-axis does not change.
	    AliMatrix(fidrot[irot],90.0,0.0,90.0+thet,90.0,thet,90.0);
	    break;
	case 1: //Rotate about y-axis, y-axis does not change.
	    AliMatrix(fidrot[irot],-90.0-thet,0.0,90.0,90.0,thet,90.0);
	    break;
	case 2: //Rotate about z-axis, z-axis does not change.
	    AliMatrix(fidrot[irot],90.0,thet,90.0,-thet-90.0,0.0,0.0);
	    break;
	default:
	    Error("Matrix","axis must be either 0, 1, or 2. for matrix=%d",
		  irot);
	    break;
    } // end if
}
//______________________________________________________________________
void AliITSv11::CreateGeometry(){
////////////////////////////////////////////////////////////////////////
//    This routine defines and Creates the geometry for version 9 of the ITS.
////////////////////////////////////////////////////////////////////////
}
