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
Revision 1.4  2003/01/29 16:01:14  nilsen
Update today's work.

Revision 1.3  2003/01/28 17:59:54  nilsen
Work continuing.

Revision 1.2  2003/01/26 14:35:15  nilsen
Some more geometry interface functions added and a start at the SSD support
cone geometry. Committed to allow easy updates of partical work between authors.

Revision 1.1  2003/01/20 23:32:49  nilsen
New ITS geometry. Only a Skeleton for now.

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
	} // end switch
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
Float_t AliITSv11::GetA(Int_t z){
    // Returns the isotopicaly averaged atomic number.
    // Inputs:
    //    Int_t z  Elemental number
    // Outputs:
    //    none.
    // Return:
    //    The atomic mass number.
    const Float_t A[]={ 1.00794 ,  4.0026902,  6.941   ,  9.012182 , 10.811   ,
                       12.01007 , 14.00674  , 15.9994  , 18.9984032, 20.1797  ,
                       22.98970 , 24.3050   , 26.981538, 28.0855   , 30.973761,
                       32.066   , 35.4527   , 39.948   , 39.0983   , 40.078   ,
                       44.95591 , 47.867    , 50.9415  , 51.9961   , 54.938049,
                       55.845   , 58.933200 , 58.6934  , 63.546    , 65.39    ,
                       69.723   , 72.61     , 74.92160 , 78.96     , 79.904   ,
                       83.80    , 85.4678   , 87.62    , 88.9085   , 91.224   ,
                       92.90638 , 95.94     , 97.907215, 101.07    ,102.90550 ,
                      106.42    ,107.8682   ,112.411   ,114.818    ,118.710   ,
                      121.760   ,127.60     ,126.90447 ,131.29     ,132.90545 ,
                      137.327   ,138.9055   ,140.116   ,140.90765  ,144.24    ,
                      144.912746,150.36     ,151.964   ,157.25     ,158.92534 ,
                      162.50     ,164.93032 ,167.26    ,168.93421  ,173.04    ,
                      174.967    ,178.49    ,180.9479 ,183.84      ,186.207   ,
                      190.23     ,192.217   ,195.078  ,196.96655   ,200.59    ,
                      204.3833   ,207.2     ,208.98038,208.982415  ,209.987131,
                      222.017570 ,223.019731,226.025402,227.027747 ,232.0381  ,
                      231.03588  238.0289};

    if(z<1||z>92){
	Error("GetA","z must be 0<z<93. z=%d",z);
	return 0.0;
    } // end if
    return A[z-1];
}
//______________________________________________________________________
Float_t AliITSv11::GetStandardMaxStepSize(Int_t istd){
    // Returns one of a set of standard Maximum Step Size values.
    // Inputs:
    //   Int_t istd  Index to indecate which standard.
    // Outputs:
    //    none.
    // Return:
    //    The appropreate standard Maximum Step Size value [cm].
    Float_t t[]={1.0, // default
	         0.0075, // Silicon detectors...
		 1.0, // Air in central detectors region
		 1.0  // Material in non-centeral region
    };
    return t[istd];
}
//______________________________________________________________________
Float_t AliITSv11::GetStandardThetaMax(Int_t istd){
    // Returns one of a set of standard Theata Max values.
    // Inputs:
    //   Int_t istd  Index to indecate which standard.
    // Outputs:
    //    none.
    // Return:
    //    The appropreate standard Theta max value [degrees].
    Float_t t[]={0.1, // default
	         0.1, // Silicon detectors...
		 0.1, // Air in central detectors region
		 1.0  // Material in non-centeral region
    };
    return t[istd];
}
//______________________________________________________________________
Float_t AliITSv11::GetStandardEfraction(Int_t istd){
    // Returns one of a set of standard E fraction values.
    // Inputs:
    //   Int_t istd  Index to indecate which standard.
    // Outputs:
    //    none.
    // Return:
    //    The appropreate standard E fraction value [#].
    Float_t t[]={0.1, // default
	         0.1, // Silicon detectors...
		 0.1, // Air in central detectors region
		 1.0  // Material in non-centeral region
    };
    return t[istd];
}
//______________________________________________________________________
void AliITSv11::Element(Int_t imat,const char* name,Int_t z,Double_t dens,
			Int_t istd){
    // Defines a Geant single element material and sets its Geant medium
    // proporties. The average atomic A is assumed to be given by their
    // natural abundances. Things like the radiation length are calculated
    // for you.
    // Inputs:
    //    Int_t imat       Material number.
    //    const char* name Material name. No need to add a $ at the end.
    //    Int_t z          The elemental number.
    //    Double_t dens    The density of the material [g/cm^3].
    //    Int_t istd       Defines which standard set of transport parameters
    //                     which should be used.
    // Output:
    //     none.
    // Return:
    //     none.
    Float_t rad,Z,A=GetA(z),tmax,stemax,deemax,epsilon;
    char *name2;
    Int_t len;

    len = strlng(name)+1;
    name2 = new char[len];
    strncpy(name2,name,len-1);
    name2[len-1] = '\0';
    name2[len-2] = '$';
    Z = (Float_t)z;
    rad = GetRadLength(z)/dens;
    AliMaterial(imat,name2,A,Z,dens,rad,0.0,0,0);
    tmax    = GetStandardTheataMax(istd);    // degree
    stemax  = GetStandardMaxStepSize(istd);  // cm
    deemax  = GetStandardEfraction(istd);     // #
    epsilon = GetStandardEpsilon(istd);
    AliMedium(imat,name2,imat,0,gAlice->Field()->Integ(),
	      gAlice->Field()->Max(),tmax,stemax,deemax,epsilon,0.0);
    delete[] name2;
}
//______________________________________________________________________
void AliITSv11::MixtureByWeight(Int_t imat,const char* name,Int_t *z,
				Double_t *w,Double_t dens,Int_t n,Int_t istd){
    // Defines a Geant material by a set of elements and weights, and sets 
    // its Geant medium proporties. The average atomic A is assumed to be 
    // given by their natural abundances. Things like the radiation length 
    // are calculated for you.
    // Inputs:
    //    Int_t imat       Material number.
    //    const char* name Material name. No need to add a $ at the end.
    //    Int_t *z         Array of The elemental numbers.
    //    Double_t *w      Array of relative weights.
    //    Double_t dens    The density of the material [g/cm^3].
    //    Int_t n          the number of elements making up the mixture.
    //    Int_t istd       Defines which standard set of transport parameters
    //                     which should be used.   
    // Output:
    //     none.
    // Return:
    //     none.
    Float_t rad,*Z,*A,tmax,stemax,deemax,epsilon;
    char *name2;
    Int_t len,i;
    Z = new Float_t[n];
    A = new Float_t[n];

    len = strlng(name)+1;
    name2 = new char[len];
    strncpy(name2,name,len-1);
    name2[len-1] = '\0';
    name2[len-2] = '$';
    for(i=0;i<n;i++){Z[i] = (Float_t)z[i];A[i] = (Float_t)GetA(z[i]);
                     W[i] = (Float_t)w[i]}
    rad = GetRadLength(z)/dens;
    AliMixture(imat,name2,A,Z,dens,n,W);
    tmax    = GetStandardTheataMax(istd);    // degree
    stemax  = GetStandardMaxStepSize(istd);  // cm
    deemax  = GetStandardEfraction(istd);     // #
    epsilon = GetStandardEpsilon(istd);
    AliMedium(imat,name2,imat,0,gAlice->Field()->Integ(),
	      gAlice->Field()->Max(),tmax,stemax,deemax,epsilon,0.0);
    delete[] name2;
}
//______________________________________________________________________
void AliITSv11::MixtureByNumber(Int_t imat,const char* name,Int_t *z,
				Int_t *w,Double_t dens,Int_t n,Int_t istd){
    // Defines a Geant material by a set of elements and number, and sets 
    // its Geant medium proporties. The average atomic A is assumed to be 
    // given by their natural abundances. Things like the radiation length 
    // are calculated for you.
    // Inputs:
    //    Int_t imat       Material number.
    //    const char* name Material name. No need to add a $ at the end.
    //    Int_t *z         Array of The elemental numbers.
    //    Int_t_t *w       Array of relative number.
    //    Double_t dens    The density of the material [g/cm^3].
    //    Int_t n          the number of elements making up the mixture.
    //    Int_t istd       Defines which standard set of transport parameters
    //                     which should be used.   
    // Output:
    //     none.
    // Return:
    //     none.
    Float_t rad,*Z,*A,tmax,stemax,deemax,epsilon;
    char *name2;
    Int_t len,i;
    Z = new Float_t[n];
    A = new Float_t[n];

    len = strlng(name)+1;
    name2 = new char[len];
    strncpy(name2,name,len-1);
    name2[len-1] = '\0';
    name2[len-2] = '$';
    for(i=0;i<n;i++){Z[i] = (Float_t)z[i];A[i] = (Float_t)GetA(z[i]);
                     W[i] = (Float_t)w[i]}
    rad = GetRadLength(z)/dens;
    AliMixture(imat,name2,A,Z,dens,-n,W);
    tmax    = GetStandardTheataMax(istd);    // degree
    stemax  = GetStandardMaxStepSize(istd);  // cm
    deemax  = GetStandardEfraction(istd);     // #
    epsilon = GetStandardEpsilon(istd);
    AliMedium(imat,name2,imat,0,gAlice->Field()->Integ(),
	      gAlice->Field()->Max(),tmax,stemax,deemax,epsilon,0.0);
    delete[] name2;
//______________________________________________________________________
void AliITSv11::SSDConeDetail(TVector3 &tran,const char moth[3],Int_t mat0){
    // Defines the volumes and materials for the ITS SSD Support cone.
    // Based on drawings ALR-0767 and ALR-0767/3. Units are in mm.
    // Inputs:
    //   Double_t zShift  The z shift to be applied to the final volume.
    // Outputs:
    //   none.
    // Return:
    //   none.
    Double_t th = 13.0; //mm, Thickness of Rohacell+carbon fiber
    Double_t ct=1.5; //mm, Carbon finber thickness
    Double_t r=15.0; // mm, Radius of curvature.
    Double_t tc=51.0; // angle of SSD cone [degrees].
    Double_t sintc=Sind(tc),costc=Cosd(tc),tantc=Tand(tc);
    Double_t z0=0.0,zcylinder=170.0,zpost=196.0;
    Double_t Routmax=0.5*985.0,RoutHole=0.5*965.0,Routmin=0.5*945.0;
    Double_t Rholemax=0.5*890.0,Rholemin=0.5*740.0;
    Double_t RPostmin=316.0,dRPost=23.0,zpostmax=196.0,phi0post=30.0;
    Double_t Rinmax=0.5*590.0,Rincylinder=0.5*597.0,RinHole=0.5*575.0,
	     Rinmin=0.5*562.0,dzin=15.0;
    Int_t nspoaks=12,ninscrews=40,npost=4,nmounts=4;
    Int_t SSDcf=man0+1; // SSD support cone Carbon Fiber materal number.
    Int_t SSDfs=mat0+2; // SSD support cone inserto stesalite 4411w.
    Int_t SSDfo=mat0+3; // SSD support cone foam, Rohacell 50A.
    Int_t SSDsw=mat0+4; // SSD support cone screw material,Stainless steal
    Int_t ncse=0; // number of screw ends (copy number)
    Int_t ncpe=0; // number of pin end (copy number)
    Int_t ncst=0; // number of screw tops (copy number)
    Double_t t; // some general angle [degrees].
    Double_t phi0=0.0,dphi=360.0,x,y,z;
    Int_t i,j,k,l,n,nz,nrad=0;

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
    nz = 7;
    Double_t *za    = new Double_t[nz];
    Double_t *rmina = new Double_t[nz];
    Double_t *rmaxa = new Double_t[nz];
    za[0]    = z0;
    rmina[0] = Routmin;
    rmaxa[0] = Routmax;
    za[1]    = za[0]+13.5-5.0 - dza; // za[2] - dza.
    rmina[1] = rmina[0];
    rmaxa[1] =rmaxa[0];
    za[2]    = za[0]+13.5-5.0; // From Drawing ALR-0767 and ALR-0767/3
    rmaxa[2] = rmaxa[0];
    za[3]    = za[1]+rc*sintc;
    rmina[3] = rmina[1]-rc*sintc;
    rmina[2] = rmina[1]+(rmina[3]-rmina[1])*(za[2]-za[1])/(za[3]-za[1]);
    za[4]    = za[2]+rc*sintc;
    rmaxa[4] = rmaxa[2]-rc*sintc;
    rmaxa[3] = rmaxa[2]+(rmaxa[4]-rmaxa[2])*(za[3]-za[2])/(za[4]-za[2]);
    rmina[5] = Rholemax;
    za[5]    = za[3]+(za[4]-za[3])*(rmina[5]-rmina[3])/(rmina[4]-rmina[3]);
    rmina[4] = rmina[3]+(rmina[5]-rmina[3])*(za[4]-za[3])/(za[5]-za[3]);
    za[6]    = th/sinth+za[5];
    rmina[6] = Rholemax;
    rmaxa[6] = rmina[6];
    rmaxa[5] = rmaxa[4]+(rmaxa[6]-rmaxa[4])*(za[5]-za[4])/(za[6]-za[4]);
    //
    PolyCone("SCA","SSD Suport cone Carbon Fiber Surface outer left",
	     phi0,dphi,nz,*z,*rmin,*rmax,SSDcf);
    Pos("SCA",1,moth,trans.x(),trans.y(),trans.z(),0);
    XMatrix(1,180.0);
    Pos("SCA",2,moth,trans.x(),trans.y(),-trans.z(),1);
    Za[0] = 1.; Wa[0] = ; // Hydrogen Content
    Za[1] = 6.; Wa[1] = ; // Carbon Content
    MixtureByWeight(SSDcf,"Carbon Fiber for SSD support cone",Z,W,dens,2);
    //
    // Now lets define the Inserto Stesalite 4411w material volume.
    nz = 6;
    Double_t *zb    = new Double_t[nz];
    Double_t *rminb = new Double_t[nz];
    Double_t *rmaxb = new Double_t[nz];
    zb[0] = z0;
    rminb[0] = rmina[0]+ct;
    rmaxb[0] = rmaxa[0]-ct;
    zb[1] = za[1];
    rminb[1] = rminb[0];
    rmaxb[1] = rmaxb[0];
    zb[2] = za[2];
    rmaxb[2] = rmaxb[1];
    zb[3] = za[4] - ct/sintc;
    rmaxb[3] = rmaxb[2] - (rc-ct)*sintc;
    zb[4] = za[3]+ct/sintc;
    rminb[4] = rminb[1]-(rc-ct)*sintc;
    rminb[2] = rminb[1]+(rminb[4]-rminb[1])*(zb[2]-zb[1])/(zb[4]-zb[1]);
    rminb[3] = rminb[1]+(rminb[4]-rminb[1])*(zb[3]-zb[1])/(zb[4]-zb[1]);
    zb[5] = zb[4]+(ct-2.*ct)/sintc;
    rminb[5] = rminb[4]+(ct-2.*ct)*tantc;
    rmaxb[5] = rminb[5];
    rmaxb[4] = rmaxb[3]+(rmaxb[5]-rmaxb[3])*(zb[4]-zb[3])/(zb[5]-zb[3]);
    PolyCone("SCB","SSD Suport cone Inserto Stesalite left edge",
	     phi0,dphi,nz,*zb,*rminb,*rmaxb,SSDfs);
    Pos("SCB",1,"SCA",0.0,.0,0.0,0);
    Za[0] = 1.; Wa[0] = ; // Hydrogen Content
    Za[1] = 6.; Wa[1] = ; // Carbon Content
    MixtureByWeight(SSDfs,"Inserto stealite 4411w for SSD support cone",
		    Z,W,dens,3);
    //
    // Now lets define the Rohacell foam material volume.
    nz = 4;
    Double_t *zc    = new Double_t[nz];
    Double_t *rminc = new Double_t[nz];
    Double_t *rmaxc = new Double_t[nz];
    zc[0] = zb[4];
    rminc[0] = rminb[4];
    rmaxc[0] = rmminc[0];
    zc[1] = zb[5];
    rmaxc[1] = rminb[5];
    zc[2] = za[5] + ct/sintc;
    rminc[2] = rmina[5]+ct; // leave space for carbon fiber covering hole.
    rminc[1] = rminc[0] +(rminc[2]-rminc[0])*(zc[1]-zc[0])/(zc[2]-zc[0]);
    zc[3] = za[6] - ct/sintc;
    rminc[3] = rmina[6]+ct;
    rmaxc[3] = rminc[3];
    rmaxc[2] = rmaxc[1]+(rmaxc[3]-rmaxc[1])*(zc[2]-zc[1])/(zc[3]-zc[1]);
    PolyCone("SCC","SSD Suport cone Rohacell foam left edge",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDfo);
    Pos("SCC",1,"SCA",0.0,.0,0.0,0);
    Za[0] = 1.; Wa[0] = ; // Hydrogen Content
    Za[1] = 6.; Wa[1] = ; // Carbon Content
    MixtureByWeight(SSDfo,"Foam core (Rohacell 50A) for SSD support cone",
		    Z,W,dens,3);
    //
    // In volume SCB, th Inserto Stesalite 4411w material volume, there
    // are a number of Stainless steel screw and pin studs which will be
    // filled with screws/studs.
    Double_t rmin=0.0,rmax=6.0,dz=0.5*10.0; // mm
    Tube("SCD","Screw+stud used to mount things to the SSD support cone",
	 rmin,rmax,dz,SSDsw);
    rmin=0.0;rmax=6.0;dz=0.5*12.0; // mm
    Tube("SCE","pin used to mount things to the SSD support cone",
	 rmin,rmax,dz,SSDsw);
    Za[0] =  6.; Wa[0] = ; // Carbon Content
    Za[1] = 25.; Wa[1] = ; // Iron Content
    MixtureByWeight(SSDsw,"Stainless steal screw, pin, and stud material",
		    Z,W,dens,3);
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
	    Double_t a[4]={0.0,2.0,5.0,7.0}; // Relative angles.
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
    //
    // There is no carbon fiber between this upper left section and the
    // SSD spoaks. We remove it by replacing it with Rohacell foam.
    t = ct/(0.5*(Rholemax+Rholemin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    dphi = 5.0 - 2.0*t; // degrees
    phi0 = 12.5+t; // degrees see drawing ALR-0767.
    nz = 4;
    Double_t *zf    = new Double_t[nz];
    Double_t *rminf = new Double_t[nz];
    Double_t *rmaxf = new Double_t[nz];
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
    PolyCone("SCF","SSD Suport cone Rohacell foam left edge",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDfo);
    Pos("SCF",1,"SCA",0.0,.0,0.0,0);
    for(i=1;i<nspoaks;i++){
	Zmatrix(irot+i,360./((Double_t)nspoaks));
	Pos("SCG",i+1,"SCA",0.0,.0,0.0,irot+i);
    } // end for i
    //=================================================================
    // Now for the spoak part of the SSD cone.
    // It is not posible to inclue the radius of curvature between
    // the spoak part and the upper left part of the SSD cone or lowwer right
    // part. This would be discribed by the following curves.
    // R = Rmax - (5mm)*Sin(t) phi = phi0+(5mm*180/(Pi*RoutHole))*Sin(t) 
    // where 0<=t<=90 For the inner curve a simular equiation holds.
    phi0 = 12.5; // degrees see drawing ALR-0767.
    dphi = 5.0; // degrees
    nz = 4;
    Double_t *zg    = new Double_t[nz];
    Double_t *rming = new Double_t[nz];
    Double_t *rmaxg = new Double_t[nz];
    zg[0] = zb[5];
    rming[0] = rmina[5];
    rmaxg[0] = rming[0];
    zg[1] = za[6];
    rming[1] = -thatc*(zg[1]-za[3])+rmina[3];
    rmaxg[1] = rmaxg[0];
    rming[2] = Rholemin;
    zg[2] = za[3]-(rming[2]-rmina[3])/tantc;
    rmaxg[2] = -thatc*(zg[2]-za[4])+rmaxa[4];
    rming[3] = rming[2];
    rmaxg[3] = rming[3];
    zg[3] = za[4]-(rmaxg[3]-rmaxa[4])/tantc;
    PolyCone("SCG","SSD spoak carbon fiber surfaces",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDcf);
    Pos("SCG",i+1,"SCA",0.0,.0,0.0,0);
    for(i=1;i<nspoaks;i++){
	Pos("SCG",i+1,"SCA",0.0,.0,0.0,irot+i);
    } // end for i
    // For the foam core.
    t = ct/(0.5*(Rholemax+Rholemin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    dphi = 5.0 - 2.0*t; // degrees
    phi0 = 12.5+t; // degrees see drawing ALR-0767.
    nz = 4;
    Double_t *zh    = new Double_t[nz];
    Double_t *rminh = new Double_t[nz];
    Double_t *rmaxh = new Double_t[nz];
    zh[0] = zf[2];
    rminh[0] = rming[0];
    rmaxh[0] = rmaxg[0];
    zh[1] = zf[3];
    rminh[1] = rming[1]-(ct/sintc-(zg[1]-zh[1]))*tantc;
    rmaxh[1] = rmaxh[0];
    zh[2] = zg[2]+ct/tanth;
    rminh[2] = rming[2];
    rmaxh[2] = rmaxg[2]-(ct/sintc-(zg[2]-zh[2]))*tantc;
    zh[3] = zg[3]-ct/sintc;
    rminh[3] = rminh[2];
    rmaxh[3] = rminh[3];
    PolyCone("SCH","SSD spoak foam core",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDfo);
    Pos("SCH",1,"SCG",0.0,.0,0.0,0);
    //
    //==================================================================
    // Now for the Inner most part of the SSD cone.
    phi0 = 0.0;
    dphi = 360.0;
    nz = 8;
    Double_t *zi = new Double_t[nz];
    Double_t *rmini = new Double_t[nz];
    Double_t *rmaxi = new Double_t[nz];
    Double_t za,rmina,rmaxa; // additional point not needed in call to pcons.
    zi[0] = zg[2];
    rmini[0] = rming[2];
    rmaxi[0] = rmini[0];
    zi[1] = zg[3];
    rmini[1] = -tantc*(zi[1]-za[3])+rmina[3];
    rmaxi[1] = rmaxi[0];
    rmini[5] = Rinmin;
    rmaxi[5] = Rinmax+rc*sintc;
    zi[5] =za[4]+(rmaxa[4]-rmaxi[5])/tantc;
    zia = zi[5]+rc*costc;
    rminia = rmini[5];
    rmaxia = Rinmax;
    zi[3] = zia-dzin;
    zi[2] = zi[3] -rc*costc;
    rmini[2] = -tantc*(zi[2]-za[3])+rmina[3];
    rmaxi[2] = -tantc*(zi[2]-za[4])+rmaxa[4];
    rmini[3] = rmini[2] -rc*costc;
    zi[4] = zi[3];
    rmini[4] = Rinmin;
    rmaxi[4] = -tantc*(zi[4]-za[4])+rmaxa[4];
    rmaxi[3] = rmaxi[4];
    zi[6] = zcylinder;
    rmini[6] = Rinmin;
    rmaxi[6] = rmaxi[5] - (zi[5]-zi[6])*(rmaxi[5]-rmaxa)/(zi[5]-zia);
    zi[7] = zi[6];
    rmini[7] = Rincylinder;
    rmaxi[7] = rmaxi[6];
    rmini[8] = Rincylinder;
    rmaxi[8] = Rini[8];
    zi[8] = zi[5]+(rmaxi[8]-rmaxi[5])*(zia-zi[5])/(rmaxa-rmaxi[5]);
    PolyCone("SCI","SSD lower/inner right part of SSD cone",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDcf);
    Pos("SCI",1,mother,0.0,.0,0.0,0);
    // Now for Inserto volume at the inner most radius.
    phi0 = 0.0;
    dphi = 360.0;
    nz = 7;
    Double_t *zk = new Double_t[nz];
    Double_t *rmink = new Double_t[nz];
    Double_t *rmaxk = new Double_t[nz];
    zk[1] = zi[3]+ct;
    zk[0] = zk[1]-(rc+ct)*costc;
    rmink[0] = rmini[3]+(rc+ct)*sintc;
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
    zk[5] = zi[5]-rc*costc-ct;
    rmink[5] = rmini[6];
    rmaxk[5] = rmini[8];
    zk[6] = zi[6];
    rmink[6] = rmini[6];
    rmaxk[6] = rmaxi[6];
    PolyCone("SCK","SSD inner most inserto material",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDfs);
    Pos("SCK",1,"SCI",0.0,.0,0.0,0);
    // Now for foam core at the inner most radius.
    phi0 = 0.0;
    dphi = 360.0;
    nz = 4;
    Double_t *zj = new Double_t[nz];
    Double_t *rminj = new Double_t[nz];
    Double_t *rmaxj = new Double_t[nz];
    rminj[0] = rmini[0]-ct;
    zj[0] = zc[0]+(zc[2]-zc[0])*(rminj[0]-rminc[0])/(rminc[2]-rminc[0]);
    rmaxj[0] = rminj[0];
    rmaxj[1] = rmaxj[0];
    zj[1] = zc[1]+(zc[3]-zc[1])*(rmaxj[1]-rmaxc[1])/(rmaxc[3]-rmaxc[1]);
    rminj[1] = rminc[0]+(rminc[2]-rminc[0])*(zj[1]-zc[0])/(zc[2]-zc[0]);
    zj[2] = zk[0];
    rminj[2] = rkmin[0];
    rmaxj[2] = rmaxc[1]+(rmaxc[3]-rmaxc[1])*(zj[2]-zc[1])/(zc[3]-zc[1]);
    zj[3] = zk[3];
    rminj[3] = rmaxk[3];
    rmaxj[3] = rminj[3];
    PolyCone("SCJ","SSD inner most foam core",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDfo);
    Pos("SCJ",1,"SCI",0.0,.0,0.0,0);
    // Now for foam core at the top of the inner most radius where 
    // the spoaks are.
    t = ct/(0.5*(Rholemax+Rholemin));// It is not posible to get the
    // carbon fiber thickness uniform in this phi direction. We can only
    // make it a fixed angular thickness.
    t *= 180.0/TMath::Pi();
    dphi = 5.0 - 2.0*t; // degrees
    phi0 = 12.5+t; // degrees see drawing ALR-0767.
    nz = 4;
    Double_t *zl = new Double_t[nz];
    Double_t *rminl = new Double_t[nz];
    Double_t *rmaxl = new Double_t[nz];
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
    PolyCone("SCL","SSD inner most foam core",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDfo);
    Pos("SCL",1,"SCI",0.0,.0,0.0,0);
    for(i=1;i<nspoaks;i++){
	Pos("SCG",i+1,"SCA",0.0,.0,0.0,irot+i);
    } // end for i
    // Now for the SSD mounting posts
    dphi = 180.0*dRPost/(RPostmin+0.5*dRPost)/TMath::Pi(); // degrees
    phi0 = phi0post; //
    nz = 3;
    Double_t *zo = new Double_t[nz];
    Double_t *rmino = new Double_t[nz];
    Double_t *rmaxo = new Double_t[nz];
    rmino[0] = RPostmin+dRPost;
    rmaxo[0] = rmino[0];
    zo[0] = za[4]+(rmaxa[4]-Rmaxo[0])/tantc;
    rmino[1] = RPostmin;
    zo[1] = za[4]+(rmaxa[4]-Rmino[1])/tantc;
    rmaxo[1] = rmaxo[0];
    zo[2] = z0+zpostmax;
    rmino[2] = RPostmin;
    rmaxo[2] = rmin0[2]+dRPost;
    PolyCone("SCO","SSD mounting post, carbon fiber",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDcf);
    Pos("SCO",1,mother,0.0,.0,0.0,0);
    for(i=1;i<nposts;i++){
	Zmatrix(irot+i,360./((Double_t)nposts));
	Pos("SCO",i+1,mother,0.0,.0,0.0,irot+i);
    } // end for i
    // Now for the SSD mounting posts
    t = 180.0*ct/(RPostmin+0.5*dRPost)/TMath::Pi();
    dphi = dphi-2.0*t; // degrees
    phi0 = phi0+t; //
    nz = 3;
    Double_t *zp = new Double_t[nz];
    Double_t *rminp = new Double_t[nz];
    Double_t *rmaxp = new Double_t[nz];
    rminp[0] = rmino[0]-ct;
    rmaxp[0] = rminp[0];
    zp[0] = za[4]+(rmaxa[4]-Rmaxp[0])/tantc;
    rminp[1] = rmino[0]+ct;
    rmaxp[1] = rmino[0]-ct;
    zp[1] = za[4]+(rmaxa[4]-Rminp[1])/tantc;
    rminp[2] = rminp[1];
    rmaxp[2] = rmaxp[1];
    zp[2] = z0-zpostmax;
    PolyCone("SCP","SSD mounting post, Inserto",
	     phi0,dphi,nz,*zc,*rminc,*rmaxc,SSDfs);
    Pos("SCP",1,"SCO",0.0,.0,0.0,0);
    // This insrto continues into the SSD cone displacing the foam
    // and the carbon fiber surface at those points where the posts are.
    nz = 4;
    Double_t *zm = new Double_t[nz];
    Double_t *rminm = new Double_t[nz];
    Double_t *rmaxm = new Double_t[nz];
    Double_t *zn = new Double_t[nz];
    Double_t *rminn = new Double_t[nz];
    Double_t *rmaxn = new Double_t[nz];
    rminm[0] = RPostmin+dRPost-ct;
    rmaxm[0] = rminm[0];
    zm[0] = zj[0]+(zj[2]-zj[0])*(rminm[0]-rminj[0])/(rminj[2]-rminj[0]);
    rmaxm[1] = rmaxm[0];
    zm[1] = zj[1]+(zj[3]-zj[1])*(rmaxm[1]-rmaxj[1])/(rmaxj[3]-rmaxj[1]);
    rminm[2] = RPostmin+ct;
    zm[2] = zj[0]+(zj[2]-zj[0])*(rminm[2]-rminj[0])/(rminj[2]-rminm[0]);
    rmaxm[2] = rjmax[1]+(rjmax[3]-rmnax[1])*(zm[2]-zj[1])/(zj[3]-zj[1]);
    rminm[3] = rminm[2];
    rmaxm[3] = rminm[3];
    zn[0] = zm[1];
    rminn[0] = rmaxm[1];
    rmaxn[0] = rminn[0];
    rmaxn[1] = rmaxn[0];
    zn[1] = za4+(rmaxa[4]-rmaxn[1])/tantc;
    rminn[1] = rmaxj[1]+(rmaxj[3]-rmaxj[1])*(zn[1]-zj[1])/(zj[3]-zj[1]);
    zn[2] = zm[3];
    rminn[2] = rminm[3];
    rmaxn[2] = -tantc*(zn[2]-za[4])+rmaxa[4];
    rminn[3] = rminn[2];
    rmaxn[3] = rminn[3];
    zn]3] = za[4]+(rmaxa[4]-rmaxn[3])/tantc;
    Pos("SCM",1,"SCJ",0.0,.0,0.0,0);
    Pos("SCN",1,"SCI",0.0,.0,0.0,0);
    for(i=1;i<nposts;i++){
	Pos("SCN",i+1,"SCJ",0.0,.0,0.0,irot+i);
	Pos("SCM",i+1,"SCI",0.0,.0,0.0,irot+i);
    } // end for i
    //
    //==============================================================
    delete[] za;delete[] rmina;delete[] rmaxa;
    delete[] zb;delete[] rminb;delete[] rmaxb;
    delete[] zc;delete[] rminc;delete[] rmaxc;
    delete[] zd;delete[] rmind;delete[] rmaxd;
    delete[] ze;delete[] rmine;delete[] rmaxe;
    delete[] zf;delete[] rminf;delete[] rmaxf;
    delete[] zg;delete[] rming;delete[] rmaxg;
    delete[] zh;delete[] rminh;delete[] rmaxh;
    delete[] zi;delete[] rmini;delete[] rmaxi;
    delete[] zj;delete[] rminj;delete[] rmaxj;
    // Set back to cm default scale before exiting.
    SetScalecm();
    return;
}
//______________________________________________________________________
void AliITSv11::CreateGeometry(){
    ////////////////////////////////////////////////////////////////////////
    // This routine defines and Creates the geometry for version 11 of the ITS.
    ////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
void AliITSv11::CreateMaterials(){
////////////////////////////////////////////////////////////////////////
  //
  // Create ITS materials
  //     This function defines the default materials used in the Geant
  // Monte Carlo simulations for the geometries AliITSv1, AliITSv3,
  // AliITSv11.
  // In general it is automatically replaced by
  // the CreatMaterials routine defined in AliITSv?. Should the function
  // CreateMaterials not exist for the geometry version you are using this
  // one is used. See the definition found in AliITSv5 or the other routine
  // for a complete definition.
  //
}
//______________________________________________________________________
void AliITSv11::InitAliITSgeom(){
    //     Based on the geometry tree defined in Geant 3.21, this
    // routine initilizes the Class AliITSgeom from the Geant 3.21 ITS geometry
    // sturture.
}

//______________________________________________________________________
void AliITSv11::Init(){
    ////////////////////////////////////////////////////////////////////////
    //     Initialise the ITS after it has been created.
    ////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
void AliITSv11::SetDefaults(){
    // sets the default segmentation, response, digit and raw cluster classes
}
//______________________________________________________________________
void AliITSv11::DrawModule(){
    ////////////////////////////////////////////////////////////////////////
    //     Draw a shaded view of the FMD version 11.
    ////////////////////////////////////////////////////////////////////////
}
//______________________________________________________________________
void AliITSv11::StepManager(){
    ////////////////////////////////////////////////////////////////////////
    //    Called for every step in the ITS, then calles the AliITShit class
    // creator with the information to be recoreded about that hit.
    //     The value of the macro ALIITSPRINTGEOM if set to 1 will allow the
    // printing of information to a file which can be used to create a .det
    // file read in by the routine CreateGeometry(). If set to 0 or any other
    // value except 1, the default behavior, then no such file is created nor
    // it the extra variables and the like used in the printing allocated.
    ////////////////////////////////////////////////////////////////////////
}
 
