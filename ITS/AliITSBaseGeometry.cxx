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

/*
  A base geometry class defining all of the ITS volumes that make up an ITS
geometry.
Auhors: B. S. Nilsen
Version 0
Created February 2003.
*/

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include <AliRun.h>
#include <AliMagF.h>
#include <AliConst.h>
#include "AliITSBaseGeometry.h"

ClassImp(AliITSBaseGeometry)
//______________________________________________________________________
AliITSBaseGeometry::AliITSBaseGeometry(){
    // Default construtor for the ITS Base Geometry class.
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.

    fScale = 1.0; // Default value.
    fits = 0; // zero pointers.
    fNCreates++; // incrament this creation counter.
}
//______________________________________________________________________
AliITSBaseGeometry::AliITSBaseGeometry(AliModule *its,Int_t iflag){
    // Standard construtor for the ITS Base Geometry class.
    // Inputs:
    //    Int_t iflag  flag to indecate specific swiches in the geometry
    // Outputs:
    //    none.
    // Return:
    //    none.

    fScale = 1.0; // Default value.
    fits = its; // get a copy of the pointer to the ITS.
    fNCreates++; // incrament this creation counter.
}
//______________________________________________________________________
AliITSBaseGeometry::~AliITSBaseGeometry(){
    // Standeard destructor for the ITS Base Geometry class.
    // Inputs:
    //    Int_t iflag  flag to indecate specific swiches in the geometry
    // Outputs:
    //    none.
    // Return:
    //    none.

    fits = 0; // This class does not own this class. It contaitns a pointer
    // to it for conveniance.
    fidmed = 0; // This class does not own this array of media indexs. It
    fNCreates--;
    if(fNCreates==0){ // Now delete the static members
	Int_t i;
	if(fVolName!=0){
	    for(i=0;i<fVolNameLast;i++) delete fVolName[i];
	    fVolNameSize = 0;
	    fVolNameLast = 0;
	    delete[] fVolName;
	}// end if
	delete[] fidrot;
	fidrotsize = fidrotlast = 0;
    }// end if
}
//______________________________________________________________________
Int_t AliITSBaseGeometry::AddVolName(const TString name){
    // Checks if the volume name already exist, if not it adds it to
    // the list of volume names and returns an index to that volume name.
    // it will create and expand the array of volume names as needed.
    // If the volume name already exists, it will give an error message and
    // return an index <0.
    // Inputs:
    //    const TString name  Volume name to be added to the list.
    // Outputs:
    //    none.
    // Return:
    //    The index where this volume name is stored.
    Int_t i;

    if(fVolName==0){ // must create array.
	fVolNameSize = 1000;
	fVolName = new TString[fVolNameSize];
	fVolNameLast = 0;
    } // end if
    for(i=0;i<fVolNameLast;i++) if(fVolName[i].CompareTo(name)==0){ // Error
	Error("AddVolName","Volume name already exists for volume %d",i);
	return -1;
    } // end for i
    if(fVolNameSize==fVolNameLast-1){ // Array is full must expand.
	Int_t size = fVolNameSize*2;
	TString *old = fVolName;
	fVolName = new TString[fVolNameSize];
	for(i=0;i<fVolNameLast;i++) fVolName[i] = old[i];
	delete[] old;
	fVolNameSize = size;
    } // end if
    if(strcmp(ITSIndexToITSG3name(fVolNameLast),"ITSV")==0){
	// Special Reserved Geant 3 volumen name. Skip it
	// fill it with explination for conveniance.
	fVolName[fVolNameLast] = "ITS Master Mother Volume";
	fVolNameLast++;
    } // end if
    fVolName[fVolNameLast] = name;
    fVolNameLast++;
    return fVolNameLast-1; // return the index
}
//______________________________________________________________________
char* AliITSBaseGeometry::ITSIndexToITSG3name(const Int_t i){
    // Given the ITS volume index i, it returns the Geant3 ITS volume
    // name. The valid characters must be in the range
    // '0' through 'Z'. This will include all upper case letter and the
    // numbers 0-9. In addition it does not will include the following simbols
    // ":;<=>?@"
    // Inputs:
    //    const Int_t i  the ITS volume index
    // Output:
    //    none.
    // Return:
    //    char[4] with the ITS volume name starting from "I000" to "IZZZ"
    const Int_t rangen=(Int_t)('9'-'0'+1); // range of numbers
    const Int_t rangel=(Int_t)('Z'-'A'+1); // range of letters
    const Int_t range = rangen+rangel; // the number of characters between 
                                       // 0-9 and A-Z.
    char a[4];
    Int_t j = i;

    a[0] = (char)('I');
    a[1] = (char)('0'+j/(range*range));
    if(a[1]>'9') a[1] += 'A'-'0'; // if it is a letter add in gap for simples.
    j -= range*range*(a[1]-'0');
    a[2] = (char)('0'+j/range);
    if(a[2]>'9') a[2] += 'A'-'0'; // if it is a letter add in gap for simples.
    j -= range*(a[2]-'0');
    a[3] = (char)('0'+j);
    if(a[3]>'9') a[3] += 'A'-'0'; // if it is a letter add in gap for simples.
    return a;
}
//______________________________________________________________________
Int_t AliITSBaseGeometry::ITSG3VnameToIndex(const char name[3])const{
    // Given the last three characters of the ITS Geant3 volume name,
    // this returns the index. The valid characters must be in the range
    // '0' through 'Z'. This will include all upper case letter and the
    // numbers 0-9. In addition it will include the following simbles
    // ":;<=>?@"
    // Inputs:
    //    const char name[3]  The last three characters of the ITS Geant3
    //                        volume name
    // Output:
    //    none.
    // Return:
    //    Int_t the index.
    const Int_t rangen=(Int_t)('9'-'0'+1); // range of numbers
    const Int_t rangel=(Int_t)('Z'-'A'+1); // range of letters
    const Int_t range = rangen+rangel; // the number of characters between 
                                       // 0-9 and A-Z.
    Int_t i,j;

    i = 0;
    for(j=3;j>-1;j--){
	if(isdigit(name[j])){ // number
	    i += (Int_t)(name[j]-'0')*TMath::Power(range,(Double_t)j);
	}else{ // Letter
	    i += (Int_t)(name[j]-'A'+rangen)*TMath::Power(range,(Double_t)j);
	} // end if
    } // end for j
    return i;
}
//______________________________________________________________________
TString AliITSBaseGeometry::GetVolName(const Int_t i)const{
    // Returns the volume name at a given index i. Index must be in
    // range and the array of volume names must exist. If there is an
    // error, a message is written and 0 is returned.
    // Inputs:
    //   const Int_t i Index
    // Output:
    //   none.
    // Return:
    //   A TString contianing the ITS volume name.

    if(i<0||i>=fVolNameLast){
	Error("GetVolName","Index=%d out of range but be witin 0<%d",i,
	      fVolName-1);
	return 0;
    } // end if Error
    return fVolName[i];
}
//______________________________________________________________________
Int_t AliITSBaseGeometry::GetVolumeIndex(const TString &a){
    // Return the index corresponding the the volume name a. If the
    // Volumen name is not found, return -1, and a warning message given.
    // Inputs:
    //   const TString &a  Name of volume for which index is wanted.
    // Output:
    //   none.
    // Return:
    //   Int_t Index corresponding the volume a. If not found -1 is returned.
    Int_t i;

    for(i=0;i<fVolNameLast;i++) if(fVolName[i].CompareTo(a)==0) return i;
    Info("GetVolumeIndex","Volume name %s not found",a.Data());
    return -1;
}
//______________________________________________________________________
void AliITSBaseGeometry::Box(const char gnam[3],const TString &dis,
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*dx;
    param[1] = fScale*dy;
    param[2] = fScale*dz;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"BOX ",fidmed[med],param,3);
}
//______________________________________________________________________
void AliITSBaseGeometry::Trapezoid1(const char gnam[3],const TString &dis,
				    Double_t dxn,Double_t dxp,Double_t dy,
				    Double_t dz,Int_t med){
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*dxn;
    param[1] = fScale*dxp;
    param[2] = fScale*dy;
    param[3] = fScale*dz;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TRD1",fidmed[med],param,4);
}
//______________________________________________________________________
void AliITSBaseGeometry::Trapezoid2(const char gnam[3],const TString &dis,
				    Double_t dxn,Double_t dxp,Double_t dyn,
				    Double_t dyp,Double_t dz,Int_t med){
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*dxn;
    param[1] = fScale*dxp;
    param[2] = fScale*dyn;
    param[3] = fScale*dyp;
    param[4] = fScale*dz;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TRD2",fidmed[med],param,5);
}
//______________________________________________________________________
void AliITSBaseGeometry::Trapezoid(const char gnam[3],const TString &dis,
				   Double_t dz,Double_t thet,Double_t phi,
				   Double_t h1,Double_t bl1,Double_t tl1,
				   Double_t alp1,Double_t h2,Double_t bl2,
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

    if(fidmed==0) SetMedArray();
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
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TRAP",fidmed[med],param,11);
}
//______________________________________________________________________
void AliITSBaseGeometry::Tube(const char gnam[3],const TString &dis,
			      Double_t rmin,Double_t rmax,Double_t dz,
			      Int_t med){
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TUBE",fidmed[med],param,3);
}
//______________________________________________________________________
void AliITSBaseGeometry::TubeSegment(const char gnam[3],const TString &dis,
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    param[3] = phi1;
    param[4] = phi2;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"TUBS",fidmed[med],param,5);
}
//______________________________________________________________________
void AliITSBaseGeometry::Cone(const char gnam[3],const TString &dis,
			      Double_t dz,Double_t rmin1,Double_t rmax1,
			      Double_t rmin2,Double_t rmax2,Int_t med){
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*dz;
    param[1] = fScale*rmin1;
    param[2] = fScale*rmax1;
    param[3] = fScale*rmin2;
    param[4] = fScale*rmax2;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"CONS",fidmed[med],param,5);
}
//______________________________________________________________________
void AliITSBaseGeometry::ConeSegment(const char gnam[3],const TString &dis,
				     Double_t dz,Double_t rmin1,Double_t rmax1,
				     Double_t rmin2,Double_t rmax2,
				     Double_t phi1,Double_t phi2,Int_t med){
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*dz;
    param[1] = fScale*rmin1;
    param[2] = fScale*rmax1;
    param[3] = fScale*rmin2;
    param[4] = fScale*rmax2;
    param[5] = phi1;
    param[6] = phi2;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"CONS",fidmed[med],param,7);
}
//______________________________________________________________________
void AliITSBaseGeometry::Sphere(const char gnam[3],const TString &dis,
				Double_t rmin,Double_t rmax,Double_t the1,
				Double_t the2,Double_t phi1,Double_t phi2,
				Int_t med){
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = the1;
    param[3] = the2;
    param[4] = phi1;
    param[5] = phi2;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"SPHE",fidmed[med],param,6);
}
//______________________________________________________________________
void AliITSBaseGeometry::Parallelepiped(const char gnam[3],const TString &dis,
					Double_t dx,Double_t dy,Double_t dz,
					Double_t alpha,Double_t thet,
					Double_t phi,Int_t med){
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*dx;
    param[1] = fScale*dy;
    param[2] = fScale*dz;
    param[3] = alpha;
    param[4] = thet;
    param[5] = phi;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"PARA",fidmed[med],param,6);
}
//______________________________________________________________________
void AliITSBaseGeometry::Polygon(const char gnam[3],const TString &dis,
				 Double_t phi1,Double_t dphi,Int_t npdv,
				 Int_t nz,Double_t *z,Double_t *rmin,
				 Double_t *rmax,Int_t med){
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

    if(fidmed==0) SetMedArray();
    n = 4+3*nz;
    param = new Float_t[n];
    param[0] = phi1;
    param[1] = dphi;
    param[2] = (Float_t)npdv;
    param[3] = (Float_t)nz;
    for(i=0;i<nz;i++){
	param[4+3*i] = z[i];
	param[5+3*i] = rmin[i];
	param[6+3*i] = rmax[i];
    } // end for i
    name[3] = 'I';
    for(i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"PGON",fidmed[med],param,n);

    delete[] param;
}
//______________________________________________________________________
void AliITSBaseGeometry::PolyCone(const char gnam[3],const TString &dis,
				  Double_t phi1,Double_t dphi,Int_t nz,
				  Double_t *z,Double_t *rmin,Double_t *rmax,
				  Int_t med){
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

    if(fidmed==0) SetMedArray();
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
    name[3] = 'I';
    for(i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"PCON",fidmed[med],param,n);

    delete[] param;
}
//______________________________________________________________________
void AliITSBaseGeometry::TubeElliptical(const char gnam[3],const TString &dis,
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*p1;
    param[1] = fScale*p2;
    param[2] = fScale*dz;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"ELTU",fidmed[med],param,3);
}
//______________________________________________________________________
void AliITSBaseGeometry::HyperbolicTube(const char gnam[3],const TString &dis,
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    param[3] = thet;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"HYPE",fidmed[med],param,4);
}
//______________________________________________________________________
void AliITSBaseGeometry::TwistedTrapezoid(const char gnam[3],
					  const TString &dis,
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

    if(fidmed==0) SetMedArray();
    param[0] = fScale*dz;
    param[1] = thet;
    param[2] = phi;
    param[3] = twist;
    param[4] = fScale*h1;
    param[5] = fScale*bl1;
    param[6] = fScale*tl1;
    param[7] = apl1;
    param[8] = fScale*h2;
    param[9] = fScale*bl2;
    param[10] = fScale*tl2;
    param[11] = apl2;
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"GTRA",fidmed[med],param,12);
}
//______________________________________________________________________
void AliITSBaseGeometry::CutTube(const char gnam[3],const TString &dis,
				 Double_t rmin,Double_t rmax,Double_t dz,
				 Double_t phi1,Double_t phi2,Double_t lx,
				 Double_t ly,Double_t lz,Double_t hx,
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

    if(fidmed==0) SetMedArray();
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
    name[3] = 'I';
    for(Int_t i=0;i<3;i++) name[i+1] = gnam[i];
    gMC->Gsvolu(name,"CTUB",fidmed[med],param,11);
}
//______________________________________________________________________
void AliITSBaseGeometry::Pos(const char vol[3],Int_t cn,const char moth[3],
			     Double_t x,Double_t y,Double_t z,Int_t irot){
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
    name[3] = 'I';
    for(i=0;i<3;i++) name[i+1] = vol[i];
    mother[3] = 'I';
    for(i=0;i<3;i++) mother[i+1] = moth[i];
    if(irot>=0) r=fidrot[irot];
    gMC->Gspos(name,1,mother,param[0],param[1],param[2],r,"ONLY");
}
//______________________________________________________________________
void AliITSBaseGeometry::Matrix(Int_t irot,Double_t thet1,Double_t phi1,
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
	p3 = phi3;
	fits->AliMatrix(fidrot[irot],t1,p1,t2,p2,t3,p3);
    } // end if
}
//______________________________________________________________________
void AliITSBaseGeometry::Matrix(Int_t irot,Int_t axis,Double_t thet){
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
	    fits->AliMatrix(fidrot[irot],90.0,0.0,90.0+thet,90.0,thet,90.0);
	    break;
	case 1: //Rotate about y-axis, y-axis does not change.
	    fits->AliMatrix(fidrot[irot],-90.0-thet,0.0,90.0,90.0,thet,90.0);
	    break;
	case 2: //Rotate about z-axis, z-axis does not change.
	    fits->AliMatrix(fidrot[irot],90.0,thet,90.0,-thet-90.0,0.0,0.0);
	    break;
	default:
	    Error("Matrix","axis must be either 0, 1, or 2. for matrix=%d",
		  irot);
	    break;
	} // end switch
    } // end if
}
//______________________________________________________________________
void AliITSBaseGeometry::Matrix(Int_t irot,Double_t rot[3][3]){
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
	fits->AliMatrix(fidrot[irot],ang[0],ang[1],ang[2],ang[3],
			ang[4],ang[5]);
    } // end if
}
//______________________________________________________________________
Float_t AliITSBaseGeometry::GetA(Int_t z){
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
                      162.50    ,164.93032  ,167.26    ,168.93421  ,173.04    ,
                      174.967   ,178.49     ,180.9479 ,183.84      ,186.207   ,
                      190.23    ,192.217    ,195.078  ,196.96655   ,200.59    ,
                      204.3833  ,207.2      ,208.98038,208.982415  ,209.987131,
                      222.017570,223.019731 ,226.025402,227.027747 ,232.0381  ,
                      231.03588 ,238.0289};

    if(z<1||z>92){
	Error("GetA","z must be 0<z<93. z=%d",z);
	return 0.0;
    } // end if
    return A[z-1];
}
//______________________________________________________________________
Float_t AliITSBaseGeometry::GetStandardMaxStepSize(Int_t istd){
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
Float_t AliITSBaseGeometry::GetStandardThetaMax(Int_t istd){
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
Float_t AliITSBaseGeometry::GetStandardEfraction(Int_t istd){
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
		 0.5  // Material in non-centeral region
    };
    return t[istd];
}
Float_t AliITSBaseGeometry::GetStandardEpsilon(Int_t istd){
    // Returns one of the standard Epsilon valuse
    // Inputs:
    //    Int_t istd  index of standard cuts to get
    // Output:
    //    none.
    // Return:
    //    Float_t the standard Epsilon cut value.
    Float_t t[]={1.0E-4, // default
		 1.0E-4, // Silicon detectors...
		 1.0E-4, // Air in central detector region
		 1.0E-3, // Material in non-cneteral regions
    };

    return t[istd];
}
//______________________________________________________________________
void AliITSBaseGeometry::Element(Int_t imat,const char* name,Int_t z,
				 Double_t dens,Int_t istd){
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

    len = strlen(name)+1;
    name2 = new char[len];
    strncpy(name2,name,len-1);
    name2[len-1] = '\0';
    name2[len-2] = '$';
    Z = (Float_t)z;
    rad = GetRadLength(z)/dens;
    fits->AliMaterial(imat,name2,A,Z,dens,rad,0.0,0,0);
    tmax    = GetStandardThetaMax(istd);    // degree
    stemax  = GetStandardMaxStepSize(istd);  // cm
    deemax  = GetStandardEfraction(istd);     // ratio
    epsilon = GetStandardEpsilon(istd);       //
    fits->AliMedium(imat,name2,imat,0,gAlice->Field()->Integ(),
		    gAlice->Field()->Max(),tmax,stemax,deemax,epsilon,0.0);
    delete[] name2;
}
//______________________________________________________________________
void AliITSBaseGeometry::MixtureByWeight(Int_t imat,const char* name,Int_t *z,
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
    Float_t *Z,*A,*W,tmax,stemax,deemax,epsilon;
    char *name2;
    Int_t len,i;
    Z = new Float_t[n];
    A = new Float_t[n];
    W = new Float_t[n];

    len = strlen(name)+1;
    name2 = new char[len];
    strncpy(name2,name,len-1);
    name2[len-1] = '\0';
    name2[len-2] = '$';
    for(i=0;i<n;i++){Z[i] = (Float_t)z[i];A[i] = (Float_t)GetA(z[i]);
                     W[i] = (Float_t)w[i];}
    fits->AliMixture(imat,name2,A,Z,dens,n,W);
    tmax    = GetStandardThetaMax(istd);    // degree
    stemax  = GetStandardMaxStepSize(istd);  // cm
    deemax  = GetStandardEfraction(istd);     // #
    epsilon = GetStandardEpsilon(istd);
    fits->AliMedium(imat,name2,imat,0,gAlice->Field()->Integ(),
	      gAlice->Field()->Max(),tmax,stemax,deemax,epsilon,0.0);
    delete[] name2;
    delete[] Z;
    delete[] A;
    delete[] W;
}
//______________________________________________________________________
void AliITSBaseGeometry::MixtureByNumber(Int_t imat,const char* name,Int_t *z,
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
    Float_t *Z,*A,*W,tmax,stemax,deemax,epsilon;
    char *name2;
    Int_t len,i;
    Z = new Float_t[n];
    A = new Float_t[n];
    W = new Float_t[n];

    len = strlen(name)+1;
    name2 = new char[len];
    strncpy(name2,name,len-1);
    name2[len-1] = '\0';
    name2[len-2] = '$';
    for(i=0;i<n;i++){Z[i] = (Float_t)z[i];A[i] = (Float_t)GetA(z[i]);
                     W[i] = (Float_t)w[i];}
    fits->AliMixture(imat,name2,A,Z,dens,-n,W);
    tmax    = GetStandardThetaMax(istd);    // degree
    stemax  = GetStandardMaxStepSize(istd);  // cm
    deemax  = GetStandardEfraction(istd);     // #
    epsilon = GetStandardEpsilon(istd);
    fits->AliMedium(imat,name2,imat,0,gAlice->Field()->Integ(),
		    gAlice->Field()->Max(),tmax,stemax,deemax,epsilon,0.0);
    delete[] name2;
    delete[] Z;
    delete[] A;
    delete[] W;
}
//______________________________________________________________________
Double_t AliITSBaseGeometry::RadLength(Int_t iz,Double_t a){
    // Computes the radiation length in accordance to the PDG 2000 Section
    // 23.4.1 p. 166. Transladed from the c code of Flavio Tosello.
    // Inputs:
    //    Int_t iz    The elemental number
    //    Dougle_t    The elemental average atomic mass number
    // Outputs:
    // Return:
    //    Double_t returns the radiation length of the element iz in
    //             [gm/cm^2].
    Double_t z = (Double_t)iz;
    Double_t alphaz = fAlpha*z;
    Double_t alphaz2 = alphaz*alphaz;
    Double_t c0 = +0.20206,c1 = -0.0369,c2 = +0.0083,c3 = -0.0020;
    Double_t z12,z23,l,lp,c;

    c = alphaz2*(1./(1.+alphaz2) + c0 + c1*alphaz2 + c2*alphaz2*alphaz2
		  +c3*alphaz2*alphaz2*alphaz2);
    z12 = TMath::Exp(TMath::Log(z)/3.0);
    z23 = z12*z12;
    switch (iz){
    case 1: //Hydrogen
	l  = 5.31;
	lp = 6.144;
	break;
    case 2: //Helium
	l  = 4.79;
	lp = 5,621;
	break;
    case 3: //Lithium
	l  = 4.74;
	lp = 5.805;
	break;
    case 4: //Berilium
	l  = 4.71;
	lp = 5.924;
	break;
    default: //Others
	l  = TMath::Log(184.15/z12);
	lp = TMath::Log(1194.0/z23);
	break;
    } // end switch
    Double_t re2,b,r,xz;

    re2 = fRe*fRe;
    b = 4.0*fAlpha*re2*fNa/a;
    r = b*z*(z*(l-c)+lp);
    xz = 1.0/r;
    return xz; // [gm/cm^2]
}
