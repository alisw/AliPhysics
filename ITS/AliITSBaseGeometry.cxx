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
#include <TVector3.h>
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

const Double_t AliITSBaseGeometry::fAlpha = 7.297352533e-3;
const Double_t AliITSBaseGeometry::fRe = 2.81794028e-13;
const Double_t AliITSBaseGeometry::fNa = 6.02214199e+23;
Int_t    AliITSBaseGeometry::fNCreates    = 0;
Int_t*   AliITSBaseGeometry::fidrot       = 0;
Int_t    AliITSBaseGeometry::fidrotsize   = 0;
Int_t    AliITSBaseGeometry::fidrotlast   = 0;
Int_t    AliITSBaseGeometry::fVolNameSize = 0;
Int_t    AliITSBaseGeometry::fVolNameLast = 0;
TString* AliITSBaseGeometry::fVolName     = 0;

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
    if(fNCreates==0){ // only for very first init
    } // end if
    fNCreates++; // incrament this creation counter.
}
//______________________________________________________________________
AliITSBaseGeometry::AliITSBaseGeometry(AliITS *its,Int_t iflag){
    // Standard construtor for the ITS Base Geometry class.
    // Inputs:
    //    Int_t iflag  flag to indecate specific swiches in the geometry
    // Outputs:
    //    none.
    // Return:
    //    none.

    fScale = iflag; // remove warning message for unused variable
    fScale = 1.0; // Default value.
    fits = its; // get a copy of the pointer to the ITS.
    if(fNCreates==0){ // only for very first init
	fidrotsize = ITSG3VnameToIndex("TSV")+1;
	fidrot = new Int_t[fidrotsize];
	fidrotlast = 0;
    } // end if
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
	fVolNameSize = 38624;
	fVolName = new TString[fVolNameSize];
	fVolNameLast = 0;
    } // end if
    for(i=0;i<fVolNameLast;i++) if(fVolName[i].CompareTo(name)==0){ // Error
	Error("AddVolName","Volume name already exists for volume %d name %s",
	      i,name.Data());
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
    i=ITSIndexToITSG3name(fVolNameLast);
    if(strcmp((char*)(&i),"ITSV")==0){
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
Int_t AliITSBaseGeometry::ITSIndexToITSG3name(const Int_t i){
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
    Int_t k;
    Byte_t *a = (Byte_t*) &k;
    Int_t j = i;

    k = 0;
    a[0] = (Byte_t)('I');
    a[1] = (Byte_t)('0'+j/(range*range));
    if(a[1]>'9') a[1] += 'A'-'9'-1;//if it is a letter add in gap for simples.
    j -= range*range*((Int_t)(j/(range*range)));
    a[2] = (Byte_t)('0'+j/range);
    if(a[2]>'9') a[2] += 'A'-'9'-1;//if it is a letter add in gap for simples.
    j -= range*((Int_t)(j/range));
    a[3] = (Byte_t)('0'+j);
    if(a[3]>'9') a[3] += 'A'-'9'-1;//if it is a letter add in gap for simples.
    return k;
}
//______________________________________________________________________
Int_t AliITSBaseGeometry::ITSG3VnameToIndex(const char *name){
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
    const Int_t rangen = (Int_t)('9'-'0'+1); // range of numbers
    const Int_t rangel = (Int_t)('Z'-'A'+1); // range of letters
    const Int_t range  = rangen+rangel; // the number of characters between
                                        // 0-9 + A-Z.
    Int_t i=0,j,k;

    k = strlen(name)-1;
    for(j=k;j>k-3;j--) if(isdigit(name[j])) // number
	i += (Int_t)((name[j]-'0')*TMath::Power((Double_t)range,
						(Double_t)(k-j)));
    else
	i += (Int_t)((name[j]-'A'+rangen)*TMath::Power((Double_t)range,
						       (Double_t)(k-j)));
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
void AliITSBaseGeometry::Box(const char *gnam,const TString &dis,
			     Double_t dx,Double_t dy,Double_t dz,Int_t med){
    // Interface to TMC->Gsvolu() for ITS box geometries. Box with faces
    // perpendicular to the axes. It has 3 paramters. See SetScale() for
    // units. Default units are geant 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t param[3];

    AddVolName(dis);
    param[0] = fScale*dx;
    param[1] = fScale*dy;
    param[2] = fScale*dz;
    G3name(gnam,name);
    gMC->Gsvolu(name,"BOX ",GetMed(med),param,3);
}
//______________________________________________________________________
void AliITSBaseGeometry::Box(AliITSBoxData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS box geometries. Box with faces
    // perpendicular to the axes. It has 3 paramters. See SetScale() for
    // units. Default units are geant 3 [cm].
    // Inputs:
    //    AliITSBoxData &d   Structure with the Box parameters defined.
    //    Int_t         med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[3];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.DxAt();
    param[1] = fScale*d.DyAt();
    param[2] = fScale*d.DzAt();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"BOX ",GetMed(med),param,3);
}
//______________________________________________________________________
void AliITSBaseGeometry::Trapezoid1(const char *gnam,const TString &dis,
				    Double_t dxn,Double_t dxp,Double_t dy,
				    Double_t dz,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TRD1 geometries. Trapezoid with the 
    // x dimension varing along z. It has 4 parameters. See SetScale() for
    // units. Default units are geant 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t param[4];

    AddVolName(dis);
    param[0] = fScale*dxn;
    param[1] = fScale*dxp;
    param[2] = fScale*dy;
    param[3] = fScale*dz;
    G3name(gnam,name);
    gMC->Gsvolu(name,"TRD1",GetMed(med),param,4);
}
//______________________________________________________________________
void AliITSBaseGeometry::Trapezoid1(AliITSTrapezoid1Data &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TRD1 geometries. Trapezoid with the 
    // x dimension varing along z. It has 4 parameters. See SetScale() for
    // units. Default units are geant 3 [cm].
    // Inputs:
    //    AliITSTrapezoid1Data &d   Structure with the Trapazoid data in it.
    //    Int_t                med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[4];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.DxAt(0);
    param[1] = fScale*d.DxAt(1);
    param[2] = fScale*d.DyAt();
    param[3] = fScale*d.DzAt();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"TRD1",GetMed(med),param,4);
}
//______________________________________________________________________
void AliITSBaseGeometry::Trapezoid2(const char *gnam,const TString &dis,
				    Double_t dxn,Double_t dxp,Double_t dyn,
				    Double_t dyp,Double_t dz,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TRD2 geometries. Trapezoid with the 
    // x and y dimension varing along z. It has 5 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t param[5];

    AddVolName(dis);
    param[0] = fScale*dxn;
    param[1] = fScale*dxp;
    param[2] = fScale*dyn;
    param[3] = fScale*dyp;
    param[4] = fScale*dz;
    G3name(gnam,name);
    gMC->Gsvolu(name,"TRD2",GetMed(med),param,5);
}
//______________________________________________________________________
void AliITSBaseGeometry::Trapezoid2(AliITSTrapezoid2Data &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TRD2 geometries. Trapezoid with the 
    // x and y dimension varing along z. It has 5 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    AliITSTrapezoid2Data &d   Structure with the Trapazoid data in it.
    //    Int_t                med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[5];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.DxAt(0);
    param[1] = fScale*d.DxAt(1);
    param[2] = fScale*d.DyAt(0);
    param[3] = fScale*d.DyAt(1);
    param[4] = fScale*d.DzAt();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"TRD2",GetMed(med),param,5);
}
//______________________________________________________________________
void AliITSBaseGeometry::Trapezoid(const char *gnam,const TString &dis,
				   Double_t dz,Double_t thet,Double_t phi,
				   Double_t h1,Double_t bl1,Double_t tl1,
				   Double_t alp1,Double_t h2,Double_t bl2,
				   Double_t tl2,Double_t alp2,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TRAP geometries. General Trapezoid, 
    // The faces perpendicular to z are trapezia and their centers are not 
    // necessarily on a line parallel to the z axis. This shape has 11 
    // parameters, but only cosidering that the faces should be planar, only
    // 9 are really independent. A check is performed on the user parameters 
    // and a message is printed in case of non-planar faces. Ignoring this
    // warning may cause unpredictable effects at tracking time. See 
    // SetScale() for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    //    Double_t alp1       angle with respect to the y axis from the 
    //                        center of the side at -h1 in y to the cetner 
    //                        of the side at +h1 in y of the face at -dz in z 
    //                        [degree].
    //    Double_t h2         half-length along y of the face at +dz
    //    Double_t bl2        half-length along x of the side at -h2 in y of
    //                        the face at +dz in z.
    //    Double_t tl2        half-length along x of the side at _h2 in y of 
    //                        the face at +dz in z.
    //    Double_t alp2       angle with respect to the y axis from the 
    //                        center of the side at -h2 in y to the center 
    //                        of the side at +h2 in y of the face at +dz in z 
    //                        [degree].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[11];

    AddVolName(dis);
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
    G3name(gnam,name);
    gMC->Gsvolu(name,"TRAP",GetMed(med),param,11);
}
//______________________________________________________________________
void AliITSBaseGeometry::Trapezoid(AliITSTrapezoidData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TRAP geometries. General Trapezoid, 
    // The faces perpendicular to z are trapezia and their centers are not 
    // necessarily on a line parallel to the z axis. This shape has 11 
    // parameters, but only cosidering that the faces should be planar, only
    // 9 are really independent. A check is performed on the user parameters 
    // and a message is printed in case of non-planar faces. Ignoring this
    // warning may cause unpredictable effects at tracking time. See 
    // SetScale() for units. Default units are geant 3 [cm].
    // Inputs:
    //    AliITSTrapezoidData &d   Structure with the Trapazoid data in it.
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[11];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.DzAt();
    param[1] = d.Theta();
    param[2] = d.Phi();
    param[3] = fScale*d.HAt(0);
    param[4] = fScale*d.Bl(0);
    param[5] = fScale*d.Tl(0);
    param[6] = d.Alpha(0);
    param[7] = fScale*d.HAt(1);
    param[8] = fScale*d.Bl(1);
    param[9] = fScale*d.Tl(1);
    param[10] = d.Alpha(1);
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"TRAP",GetMed(med),param,11);
}
//______________________________________________________________________
void AliITSBaseGeometry::TwistedTrapezoid(const char *gnam,
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
    // parallel. It is a TRAP shape, exept that it is twisted in the x-y 
    // plane as a function of z. The parallel sides perpendicular to the x 
    // axis are rotated with respect to the x axis by an angle TWIST, which 
    // is one of the parameters. The shape is defined by the eight corners 
    // and is assumed to be constructed of straight lines joingin points on 
    // the boundry of the trapezoidal face at Z=-dz to the coresponding 
    // points on the face at z=+dz. Divisions are not allowed. It has 12 
    // parameters. See SetScale() for units. Default units are geant 3 [cm].
    // Note: This shape suffers from the same limitations than the TRAP. The
    // tracking routines assume that the faces are planar, but htis
    // constraint is not easily expressed in terms of the 12 parameters. 
    // Additionally, no check on th efaces is performed in this case. Users 
    // should avoid to use this shape as much as possible, and if they have
    // to do so, they should make sure that the faces are really planes. 
    // If this is not the case, the result of the trasport is unpredictable. 
    // To accelerat ethe computations necessary for trasport, 18 additioanl 
    // parameters are calculated for this shape are 1 DXODZ dx/dz of the 
    // line joing the centers of the faces at z=+_dz. 2 DYODZ dy/dz of the 
    // line joing the centers of the faces at z=+_dz.
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
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    //    Double_t apl2       angle with respect to the y axis from the 
    //                        center of the side at -h2 in y to the center 
    //                        of the side at +h2 in y of the face at +dz in 
    //                        z [degrees].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[12];

    AddVolName(dis);
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
    G3name(gnam,name);
    gMC->Gsvolu(name,"GTRA",GetMed(med),param,12);
}
//______________________________________________________________________
void AliITSBaseGeometry::TwistedTrapezoid(AliITSTrapezoidTwistedData &d,
					  Int_t med){
    // Interface to TMC->Gsvolu() for ITS GTRA geometries. General twisted 
    // trapazoid. The faces perpendicular to z are trapazia and their centers 
    // are not necessarily on a line parallel to the z axis as the TRAP. 
    // Additionally, the faces may be twisted so that none of their edges are 
    // parallel. It is a TRAP shape, exept that it is twisted in the x-y 
    // plane as a function of z. The parallel sides perpendicular to the x 
    // axis are rotated with respect to the x axis by an angle TWIST, which 
    // is one of the parameters. The shape is defined by the eight corners 
    // and is assumed to be constructed of straight lines joingin points on 
    // the boundry of the trapezoidal face at Z=-dz to the coresponding 
    // points on the face at z=+dz. Divisions are not allowed. It has 12 
    // parameters. See SetScale() for units. Default units are geant 3 [cm].
    // Note: This shape suffers from the same limitations than the TRAP. The
    // tracking routines assume that the faces are planar, but htis
    // constraint is not easily expressed in terms of the 12 parameters. 
    // Additionally, no check on th efaces is performed in this case. Users 
    // should avoid to use this shape as much as possible, and if they have
    // to do so, they should make sure that the faces are really planes. 
    // If this is not the case, the result of the trasport is unpredictable. 
    // To accelerat ethe computations necessary for trasport, 18 additioanl 
    // parameters are calculated for this shape are 1 DXODZ dx/dz of the 
    // line joing the centers of the faces at z=+_dz. 2 DYODZ dy/dz of the 
    // line joing the centers of the faces at z=+_dz.
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
    //    AliITSTrapezoidTwistedData &d   Structure with the tube parameters
    //    Int_t                      med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[12];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.DzAt();
    param[1] = d.Theta();
    param[2] = d.Phi();
    param[3] = d.Twist();
    param[4] = fScale*d.HAt(0);
    param[5] = fScale*d.Bl(0);
    param[6] = fScale*d.Tl(0);
    param[7] = d.Alpha(0);
    param[8] = fScale*d.HAt(1);
    param[9] = fScale*d.Bl(1);
    param[10] = fScale*d.Tl(1);
    param[11] = d.Alpha(1);
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"GTRA",GetMed(med),param,12);
}
//______________________________________________________________________
void AliITSBaseGeometry::Tube(const char *gnam,const TString &dis,
			      Double_t rmin,Double_t rmax,Double_t dz,
			      Int_t med){
    // Interface to TMC->Gsvolu() for ITS TUBE geometries. Simple Tube. It has
    // 3 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t param[3];

    AddVolName(dis);
    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    G3name(gnam,name);
    gMC->Gsvolu(name,"TUBE",GetMed(med),param,3);
}
//______________________________________________________________________
void AliITSBaseGeometry::Tube(AliITSTubeData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TUBE geometries. Simple Tube. It has
    // 3 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    AliITSTubeData &d    Structure with the tube parameters
    //    Int_t          med   media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[3];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.Rmin();
    param[1] = fScale*d.Rmax();
    param[2] = fScale*d.DzAt();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"TUBE",GetMed(med),param,3);
}
//______________________________________________________________________
void AliITSBaseGeometry::TubeSegment(const char *gnam,const TString &dis,
				     Double_t rmin,Double_t rmax,Double_t dz,
				     Double_t phi1,Double_t phi2,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TUBE geometries. Phi segment of a 
    // tube. It has 5  parameters. Phi1 should be smaller than phi2. If this
    // is not the case, the system adds 360 degrees to phi2. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t param[5];

    AddVolName(dis);
    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    param[3] = phi1;
    param[4] = phi2;
    G3name(gnam,name);
    gMC->Gsvolu(name,"TUBS",GetMed(med),param,5);
}
//______________________________________________________________________
void AliITSBaseGeometry::TubeSegment(AliITSTubeSegData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS TUBE geometries. Phi segment of a 
    // tube. It has 5  parameters. Phi1 should be smaller than phi2. If this
    // is not the case, the system adds 360 degrees to phi2. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    AliITSTubeSegData &d   Structure with the tube parameters
    //    Int_t             med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[5];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.Rmin();
    param[1] = fScale*d.Rmax();
    param[2] = fScale*d.DzAt();
    param[3] = d.Phi0();
    param[4] = d.Phi1();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"TUBS",GetMed(med),param,5);
}
//______________________________________________________________________
void AliITSBaseGeometry::CutTube(const char *gnam,const TString &dis,
				 Double_t rmin,Double_t rmax,Double_t dz,
				 Double_t phi1,Double_t phi2,Double_t lx,
				 Double_t ly,Double_t lz,Double_t hx,
				 Double_t hy,Double_t hz,Int_t med){
    // Interface to TMC->Gsvolu() for ITS CTUB geometries. Cut tube. A tube 
    // cut at the extremities with planes not necessarily perpendicular to 
    // the z axis. It has 11 parameters. See SetScale() for units. Default 
    // units are geant 3 [cm]. phi1 should be smaller than phi2. If this is 
    // not the case, the system adds 360 degrees to phi2.
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t rmin       Inner radius at z=0 where tube is narrowest.
    //    Double_t rmax       Outer radius at z=0 where tube is narrowest.
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
    char name[5];
    Float_t param[11];

    AddVolName(dis);
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
    G3name(gnam,name);
    gMC->Gsvolu(name,"CTUB",GetMed(med),param,11);
}
//______________________________________________________________________
void AliITSBaseGeometry::CutTube(AliITSTubeCutData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS CTUB geometries. Cut tube. A tube 
    // cut at the extremities with planes not necessarily perpendicular to 
    // the z axis. It has 11 parameters. See SetScale() for units. Default 
    // units are geant 3 [cm]. phi1 should be smaller than phi2. If this is 
    // not the case, the system adds 360 degrees to phi2.
    // Inputs:
    //    AliITSTubeCutData &d    Structure with the tube parameters
    //    Int_t             med    media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[11];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.Rmin();
    param[1] = fScale*d.Rmax();
    param[2] = fScale*d.DzAt();
    param[3] = d.Phi0();
    param[4] = d.Phi1();
    param[5] = d.Normal(0,0);
    param[6] = d.Normal(0,1);
    param[7] = d.Normal(0,2);
    param[8] = d.Normal(1,0);
    param[9] = d.Normal(1,1);
    param[10] = d.Normal(1,2);
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"CTUB",GetMed(med),param,11);
}
//______________________________________________________________________
void AliITSBaseGeometry::TubeElliptical(const char *gnam,const TString &dis,
			       Double_t p1,Double_t p2,Double_t dz,Int_t med){
    // Interface to TMC->Gsvolu() for ITS ELTU geometries. Elliptical 
    // cross-section Tube. It has 3 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm]. The equation of the surface 
    // is x^2 * p1^-2 + y^2 * p2^-2 = 1.
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t param[3];

    AddVolName(dis);
    param[0] = fScale*p1;
    param[1] = fScale*p2;
    param[2] = fScale*dz;
    G3name(gnam,name);
    gMC->Gsvolu(name,"ELTU",GetMed(med),param,3);
}
//______________________________________________________________________
void AliITSBaseGeometry::TubeElliptical(AliITSTubeEllipticalData &d,
					Int_t med){
    // Interface to TMC->Gsvolu() for ITS ELTU geometries. Elliptical 
    // cross-section Tube. It has 3 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm]. The equation of the surface 
    // is x^2 * p1^-2 + y^2 * p2^-2 = 1.
    // Inputs:
    //    AliITSTubeElipticData &d  Structure with the tube parameters
    //    Int_t                med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[3];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.P0();
    param[1] = fScale*d.P1();
    param[2] = fScale*d.DzAt();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"ELTU",GetMed(med),param,3);
}
//______________________________________________________________________
void AliITSBaseGeometry::HyperbolicTube(const char *gnam,const TString &dis,
			       Double_t rmin,Double_t rmax,Double_t dz,
			       Double_t thet,Int_t med){
    // Interface to TMC->Gsvolu() for ITS HYPE geometries. Hyperbolic tube. 
    // Fore example the inner and outer surfaces are hyperboloids, as would 
    // be foumed by a system of cylinderical wires which were then rotated 
    // tangentially about their centers. It has 4 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm]. The hyperbolic surfaces are 
    // given by r^2 = (ztan(thet)^2 + r(z=0)^2.
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t param[4];

    AddVolName(dis);
    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = fScale*dz;
    param[3] = thet;
    G3name(gnam,name);
    gMC->Gsvolu(name,"HYPE",GetMed(med),param,4);
}
//______________________________________________________________________
void AliITSBaseGeometry::HyperbolicTube(AliITSTubeHyperbolicData &d,
					Int_t med){
    // Interface to TMC->Gsvolu() for ITS HYPE geometries. Hyperbolic tube. 
    // Fore example the inner and outer surfaces are hyperboloids, as would 
    // be foumed by a system of cylinderical wires which were then rotated 
    // tangentially about their centers. It has 4 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm]. The hyperbolic surfaces are 
    // given by r^2 = (ztan(thet)^2 + r(z=0)^2.
    // Inputs:
    //    AliITSTubeHyperbolicData &d  Structure with the tube parameters
    //    Int_t                    med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[4];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.Rmin();
    param[1] = fScale*d.Rmax();
    param[2] = fScale*d.DzAt();
    param[3] = d.Theta();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"HYPE",GetMed(med),param,4);
}
//______________________________________________________________________
void AliITSBaseGeometry::Cone(const char *gnam,const TString &dis,
			      Double_t dz,Double_t rmin1,Double_t rmax1,
			      Double_t rmin2,Double_t rmax2,Int_t med){
    // Interface to TMC->Gsvolu() for ITS Cone geometries. Conical tube. It 
    // has 5 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t param[5];

    AddVolName(dis);
    param[0] = fScale*dz;
    param[1] = fScale*rmin1;
    param[2] = fScale*rmax1;
    param[3] = fScale*rmin2;
    param[4] = fScale*rmax2;
    G3name(gnam,name);
    gMC->Gsvolu(name,"CONS",GetMed(med),param,5);
}
//______________________________________________________________________
void AliITSBaseGeometry::Cone(AliITSConeData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS Cone geometries. Conical tube. It 
    // has 5 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    AliITSConeData &d  Structure with the tube parameters
    //    Int_t         med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[5];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.DzAt();
    param[1] = fScale*d.Rmin0();
    param[2] = fScale*d.Rmax0();
    param[3] = fScale*d.Rmin1();
    param[4] = fScale*d.Rmax1();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"CONS",GetMed(med),param,5);
}
//______________________________________________________________________
void AliITSBaseGeometry::ConeSegment(const char *gnam,const TString &dis,
				     Double_t dz,Double_t rmin1,
				     Double_t rmax1,Double_t rmin2,
				     Double_t rmax2,Double_t phi1,
				     Double_t phi2,Int_t med){
    // Interface to TMC->Gsvolu() for ITS ConS geometries. One segment of a 
    // conical tube. It has 7 parameters. Phi1 should be smaller than phi2. 
    // If this is not the case, the system adds 360 degrees to phi2. See 
    // SetScale() for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that 
    //                        this is an ITS volume.
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
    char name[5];
    Float_t param[7];

    AddVolName(dis);
    param[0] = fScale*dz;
    param[1] = fScale*rmin1;
    param[2] = fScale*rmax1;
    param[3] = fScale*rmin2;
    param[4] = fScale*rmax2;
    param[5] = phi1;
    param[6] = phi2;
    G3name(gnam,name);
    gMC->Gsvolu(name,"CONS",GetMed(med),param,7);
}
//______________________________________________________________________
void AliITSBaseGeometry::ConeSegment(AliITSConeSegData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS ConS geometries. One segment of a 
    // conical tube. It has 7 parameters. Phi1 should be smaller than phi2. 
    // If this is not the case, the system adds 360 degrees to phi2. See 
    // SetScale() for units. Default units are geant 3 [cm].
    // Inputs:
    //    AliITSConeSegData &d   Structure with the tube parameters
    //    Int_t             med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[7];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.DzAt();
    param[1] = fScale*d.Rmin0();
    param[2] = fScale*d.Rmax0();
    param[3] = fScale*d.Rmin1();
    param[4] = fScale*d.Rmax1();
    param[5] = d.Phi0();
    param[6] = d.Phi1();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"CONS",GetMed(med),param,7);
}
//______________________________________________________________________
void AliITSBaseGeometry::PolyCone(const char *gnam,const TString &dis,
				  Double_t phi1,Double_t dphi,Int_t nz,
				  Double_t *z,Double_t *rmin,Double_t *rmax,
				  Int_t med){
    // Interface to TMC->Gsvolu() for ITS PCON geometry. Poly-cone It has 9 
    // parameters or more. See SetScale() for units. Default units are geant
    // 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t *param;
    Int_t n,i;

    AddVolName(dis);
    n = 3+3*nz;
    param = new Float_t[n];
    param[0] = phi1;
    param[1] = dphi;
    param[2] = (Float_t) nz;
    for(i=0;i<nz;i++){
	param[3+3*i] = fScale*z[i];
	param[4+3*i] = fScale*rmin[i];
	param[5+3*i] = fScale*rmax[i];
    } // end for i
    G3name(gnam,name);
    gMC->Gsvolu(name,"PCON",GetMed(med),param,n);

    delete[] param;
}
//______________________________________________________________________
void AliITSBaseGeometry::PolyCone(AliITSPConeData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS PCON geometry. Poly-cone It has 9 
    // parameters or more. See SetScale() for units. Default units are geant
    // 3 [cm].
    // Inputs:
    //    AliITSPConeData &d  Object with poly cone data stored in it.
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t *param;
    Int_t n,i,k;
    char *j = (char *) &k;

    n = 3+3*d.Nz();
    param = new Float_t[n];
    param[0] = d.Phi0();
    param[1] = d.DPhi();
    param[2] = (Float_t) d.Nz();
    for(i=0;i<d.Nz();i++){
	param[3+3*i] = fScale*d.ZAt(i);
	param[4+3*i] = fScale*d.Rmin(i);
	param[5+3*i] = fScale*d.Rmax(i);
    } // end for if
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"PCON",GetMed(med),param,n);

    delete[] param;
}
//______________________________________________________________________
void AliITSBaseGeometry::Sphere(const char *gnam,const TString &dis,
				Double_t rmin,Double_t rmax,Double_t the1,
				Double_t the2,Double_t phi1,Double_t phi2,
				Int_t med){
    // Interface to TMC->Gsvolu() for ITS SPHE geometries. Segment of a 
    // sphereical shell. It has 6 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    char name[5];
    Float_t param[6];

    AddVolName(dis);
    param[0] = fScale*rmin;
    param[1] = fScale*rmax;
    param[2] = the1;
    param[3] = the2;
    param[4] = phi1;
    param[5] = phi2;
    G3name(gnam,name);
    gMC->Gsvolu(name,"SPHE",GetMed(med),param,6);
}
//______________________________________________________________________
void AliITSBaseGeometry::Sphere(AliITSSphereData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS SPHE geometries. Segment of a 
    // sphereical shell. It has 6 parameters. See SetScale() 
    // for units. Default units are geant 3 [cm].
    // Inputs:
    //    AliITSSphereData &d   Structure with the tube parameters
    //    Int_t            med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[6];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.Rmin();
    param[1] = fScale*d.Rmax();
    param[2] = d.Theta0();
    param[3] = d.Theta1();
    param[4] = d.Phi0();
    param[5] = d.Phi1();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"SPHE",GetMed(med),param,6);
}
//______________________________________________________________________
void AliITSBaseGeometry::Parallelepiped(const char *gnam,const TString &dis,
					Double_t dx,Double_t dy,Double_t dz,
					Double_t alpha,Double_t thet,
					Double_t phi,Int_t med){
    // Interface to TMC->Gsvolu() for ITS PARA geometries. Parallelepiped. It 
    // has 6 parameters. See SetScale() for units. Default units are geant 3 
    // [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
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
    //    Double_t phi        azimuthal angle of teh line joing the centers 
    //                        of the faaces at -dz and +dz in z [degree].
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[6];

    AddVolName(dis);
    param[0] = fScale*dx;
    param[1] = fScale*dy;
    param[2] = fScale*dz;
    param[3] = alpha;
    param[4] = thet;
    param[5] = phi;
    G3name(gnam,name);
    gMC->Gsvolu(name,"PARA",GetMed(med),param,6);
}
//______________________________________________________________________
void AliITSBaseGeometry::Parallelepiped(AliITSParallelpipedData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS PARA geometries. Parallelepiped. It 
    // has 6 parameters. See SetScale() for units. Default units are geant 3 
    // [cm].
    // Inputs:
    //    AliITSParrellepipedData &d  Structre witht the volume data in it.
    //    Int_t                   med  media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t param[6];
    Int_t i,k;
    char *j = (char *) &k;

    param[0] = fScale*d.DxAt();
    param[1] = fScale*d.DyAt();
    param[2] = fScale*d.DzAt();
    param[3] = d.Alpha();
    param[4] = d.Theta();
    param[5] = d.Phi();
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"PARA",GetMed(med),param,6);
}
//______________________________________________________________________
void AliITSBaseGeometry::PolyGon(const char *gnam,const TString &dis,
				 Double_t phi1,Double_t dphi,Int_t npdv,
				 Int_t nz,Double_t *z,Double_t *rmin,
				 Double_t *rmax,Int_t med){
    // Interface to TMC->Gsvolu() for ITS PGON geometry. Polygon It has 10 
    // parameters or more. See SetScale() for units. Default units are geant
    // 3 [cm].
    // Inputs:
    //    const char *gnam  3 character geant volume name. The letter "I"
    //                        is appended to the front to indecate that this
    //                        is an ITS volume.
    //    TString &dis        String containging part discription.
    //    Double_t phi1       the azimuthal angle at which the volume begins 
    //                        (angles are counted clouterclockwise) [degrees].
    //    Double_t dphi       opening angle of the volume, which extends from 
    //                        phi1 to phi1+dphi [degree].
    //    Int_t npdv          the number of sides of teh cross section 
    //                        between the given phi limits.
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
    char name[5];
    Float_t *param;
    Int_t n,i;

    AddVolName(dis);
    n = 4+3*nz;
    param = new Float_t[n];
    param[0] = phi1;
    param[1] = dphi;
    param[2] = (Float_t)npdv;
    param[3] = (Float_t)nz;
    for(i=0;i<nz;i++){
	param[4+3*i] = fScale*z[i];
	param[5+3*i] = fScale*rmin[i];
	param[6+3*i] = fScale*rmax[i];
    } // end for i
    G3name(gnam,name);
    gMC->Gsvolu(name,"PGON",GetMed(med),param,n);

    delete[] param;
}
//______________________________________________________________________
void AliITSBaseGeometry::PolyGon(AliITSPGonData &d,Int_t med){
    // Interface to TMC->Gsvolu() for ITS PCON geometry. Poly-cone It has 9 
    // parameters or more. See SetScale() for units. Default units are geant
    // 3 [cm].
    // Inputs:
    //    AliITSPGonData &d  Object with poly cone data stored in it.
    //    Int_t    med        media index number.
    // Output:
    //    none.
    // Return.
    //    none.
    char name[5];
    Float_t *param;
    Int_t n,i,k;
    char *j = (char *) &k;

    n = 4+3*d.Nz();
    param = new Float_t[n];
    param[0] = d.Phi0();
    param[1] = d.DPhi();
    param[2] = (Float_t) d.NPhi();
    param[3] = (Float_t) d.Nz();
    for(i=0;i<d.Nz();i++){
	param[4+3*i] = fScale*d.ZAt(i);
	param[5+3*i] = fScale*d.Rmin(i);
	param[6+3*i] = fScale*d.Rmax(i);
    } // end for i
    d.SetVid(AddVolName((d.GetName())->Data()));
    k = ITSIndexToITSG3name(d.GetVid());
    for(i=0;i<4;i++) name[i] = j[i];
    name[4] = '\0';
    gMC->Gsvolu(name,"PGON",GetMed(med),param,n);

    delete[] param;
}
//______________________________________________________________________
void AliITSBaseGeometry::Pos(AliITSBaseVolParams &v,Int_t cn,
			     AliITSBaseVolParams &m,
			     TVector3 &t,Int_t irot){
    // Place a copy of a volume previously defined by a call to GSVOLU inside 
    // its mother volulme moth.
    // Inputs:
    //   const char vol[3]  3 character geant volume name. The letter "I"
    //                      is appended to the front to indecate that this
    //                      is an ITS volume.
    //   const char moth[3] 3 character geant volume name of the mother 
    //                      volume in which vol will be placed. The letter 
    //                      "I" is appended to the front to indecate that 
    //                      this is an ITS volume.
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
    char name[5],mother[5];
    Float_t param[3];
    Int_t r=0,i;
    char *n = (char*)&r;

    param[0] = fScale*t.X();
    param[1] = fScale*t.Y();
    param[2] = fScale*t.Z();
    r = ITSIndexToITSG3name(v.GetVid());
    for(i=0;i<4;i++) name[i] = n[i]; name[4] ='\0';
    r = ITSIndexToITSG3name(m.GetVid());
    for(i=0;i<4;i++) mother[i] = n[i]; mother[4] ='\0';
    if(irot>0) r = fidrot[irot]; else r=0;
    gMC->Gspos(name,cn,mother,param[0],param[1],param[2],r,"ONLY");
}
//______________________________________________________________________
void AliITSBaseGeometry::Pos(const char *vol,Int_t cn,const char *moth,
			     Double_t x,Double_t y,Double_t z,Int_t irot){
    // Place a copy of a volume previously defined by a call to GSVOLU inside 
    // its mother volulme moth.
    // Inputs:
    //   const char vol[3]  3 character geant volume name. The letter "I"
    //                      is appended to the front to indecate that this
    //                      is an ITS volume.
    //   const char moth[3] 3 character geant volume name of the mother 
    //                      volume in which vol will be placed. The letter 
    //                      "I" is appended to the front to indecate that 
    //                      this is an ITS volume.
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
    char name[5],mother[5];
    Float_t param[3];
    Int_t r=0;

    param[0] = fScale*x;
    param[1] = fScale*y;
    param[2] = fScale*z;
    G3name(vol,name);
    G3name(moth,mother);
    if(irot>0) r = fidrot[irot];
    gMC->Gspos(name,cn,mother,param[0],param[1],param[2],r,"ONLY");
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
    Float_t t1=0.0,p1=0.0,t2=0.0,p2=0.0,t3=0.0,p3=0.0;

    if(thet1==90.0&&phi1== 0.0&&
       thet2==90.0&&phi2==90.0&&
       thet3== 0.0&&phi3== 0.0){
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
    cout << "Matrix:: fidrot["<<irot<<"]="<<fidrot[irot];
    cout <<" angles="<<t1<<" "<<p1<<" "<<t2<<" "<<p2<<" "<<t3<< " "<<p3<<endl;
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
	switch (axis) {
	case 0: //Rotate about x-axis, x-axis does not change.
	    fits->AliMatrix(fidrot[irot],90.0,0.0,90.0+thet,90.0,thet,90.0);
	    /*
	    cout << "Matrix:: axis="<<axis<<" fidrot["<<irot<<"]=";
	    cout <<fidrot[irot];
	    cout <<" angles="<<90.0<<" "<<0.0<<" "<<90.0+thet<<" "<<90.0;
	    cout <<" "<<thet<< " "<<90.0<<endl;
	    */
	    break;
	case 1: //Rotate about y-axis, y-axis does not change.
	    fits->AliMatrix(fidrot[irot],90.0-thet,0.0,90.0,90.0,-thet,0.0);
	    /*
	    cout << "Matrix:: axis="<<axis<<" fidrot["<<irot<<"]=";
	    cout << fidrot[irot];
	    cout <<" angles="<<90.-thet<<" "<<0.0<<" "<<90.0<<" "<<90.0;
	    cout <<" "<<-thet<< " "<<0.0<<endl;
	    */
	    break;
	case 2: //Rotate about z-axis, z-axis does not change.
	    fits->AliMatrix(fidrot[irot],90.0,thet,90.0,90.+thet,0.0,0.0);
	    /*
	    cout << "Matrix:: axis="<<axis<<" fidrot["<<irot<<"]=";
	    cout <<fidrot[irot];
	    cout <<" angles="<<90.0<<" "<<thet<<" "<<90.0<<" "<<90.0+thet;
	    cout <<" "<<0.0<< " "<<0.0<<endl;
	    */
	    break;
	default:
	    Error("Matrix","axis must be either 0, 1, or 2. for matrix=%d",
		  irot);
	    /*
	    cout << "Matrix:: axis="<<axis<<" fidrot["<<irot<<"]=";
	    cout <<fidrot[irot];
	    cout <<" thet=" << thet<< endl;
	    */
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
	Double_t si,c=180./TMath::Pi();
	Double_t ang[6]={90.0,0.0,90.0,90.0,0.0,0.0};

    if(rot[0][0]==1.0&&rot[1][1]==1.0&&rot[2][2]==1.0&&
       rot[0][1]==0.0&&rot[0][2]==0.0&&rot[1][0]==0.0&&
       rot[1][2]==0.0&&rot[2][0]==0.0&&rot[2][1]==0.0){
	fidrot[irot] = 0; // Unit matrix
    }else{
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
    cout << "Matrix rot[3][3]:: fidrot["<<irot<<"]="<<fidrot[irot];
    cout <<" angles="<<ang[0]<<" "<<ang[1]<<" "<<ang[2]<<" "<<
	ang[3]<<" "<<ang[4]<< " "<<ang[5]<<endl;
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
    const Float_t A[]={
	  1.00794 ,  4.0026902,  6.941   ,  9.012182 , 10.811   , // H-B
         12.01007 , 14.00674  , 15.9994  , 18.9984032, 20.1797  , // C-Ne
         22.98970 , 24.3050   , 26.981538, 28.0855   , 30.973761, // Na-P
	 32.066   , 35.4527   , 39.948   , 39.0983   , 40.078   , // S-Ca
	 44.95591 , 47.867    , 50.9415  , 51.9961   , 54.938049, // Sc-Mn
	 55.845   , 58.933200 , 58.6934  , 63.546    , 65.39    , // Fe-Zn
	 69.723   , 72.61     , 74.92160 , 78.96     , 79.904   , // Ga-Br
	 83.80    , 85.4678   , 87.62    , 88.9085   , 91.224   , // Kr-Zr
	 92.90638 , 95.94     , 97.907215, 101.07    ,102.90550 , // Nb-Rh
	106.42    ,107.8682   ,112.411   ,114.818    ,118.710   , // Pd-Sn
	121.760   ,127.60     ,126.90447 ,131.29     ,132.90545 , // Sb-Cs
	137.327   ,138.9055   ,140.116   ,140.90765  ,144.24    , // La-Nd
	144.912746,150.36     ,151.964   ,157.25     ,158.92534 , // Pm-Tb
	162.50    ,164.93032  ,167.26    ,168.93421  ,173.04    , // Dy-Yb
	174.967   ,178.49     ,180.9479 ,183.84      ,186.207   , // Lu-Re
	190.23    ,192.217    ,195.078  ,196.96655   ,200.59    , // Os-Hg
	204.3833  ,207.2      ,208.98038,208.982415  ,209.987131, // Tl-At
	222.017570,223.019731 ,226.025402,227.027747 ,232.0381  , // Rn-Th
        231.03588 ,238.0289   }; // Pa,U

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
//______________________________________________________________________
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

    len = strlen(name)+2;
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
//======================================================================
ClassImp(AliITSBaseVolParams)
//______________________________________________________________________
void AliITSBaseVolParams::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

    *os<<"Volume Id="<<fVol<<" Copy="<<fCpn<<" Name: "<<fName<<endl;
}
//______________________________________________________________________
void AliITSBaseVolParams::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    is->get(s,10);
    *is >> fVol;
    is->get(s,6);
    *is >> fCpn;
    is->get(s,7);
    *is >> fName;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSBaseVolParams &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os            The output stream
    //    AliITSBaseVolParams &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSBaseVolParams &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is            The input stream
    //    AliITSBaseVolParams &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSBoxData)
//______________________________________________________________________
void AliITSBoxData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fDx=" << fDx << " fDy=" << fDy << " fDz=" << fDz << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSBoxData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);
    is->get(s,4);
    *is >> fDx;
    is->get(s,5);
    *is >> fDy;
    is->get(s,5);
    *is >> fDz;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSBoxData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os      The output stream
    //    AliITSBoxData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSBoxData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is      The input stream
    //    AliITSBoxData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSTrapezoid1Data)
//______________________________________________________________________
void AliITSTrapezoid1Data::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fDx[0]=" << fDx[0]<< " fDx[1]=" << fDx[1]  << " fDy=" << fDy;
    *os << " fDz=" << fDz << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSTrapezoid1Data::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);
    is->get(s,7);
    *is >> fDx[0];
    is->get(s,8);
    *is >> fDx[1];
    is->get(s,5);
    *is >> fDy;
    is->get(s,5);
    *is >> fDz;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTrapezoid1Data &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os      The output stream
    //    AliITSBoxData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTrapezoid1Data &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSPGonData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSTrapezoid2Data)
//______________________________________________________________________
void AliITSTrapezoid2Data::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os <<  "fDx[0]=" << fDx[0]<< " fDx[1]=" << fDx[1];
    *os << " fDy[0]=" << fDy[0] << " fDy[1]=" << fDy[1];
    *os << " fDz=" << fDz << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSTrapezoid2Data::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);
    is->get(s,7);
    *is >> fDx[0];
    is->get(s,8);
    *is >> fDx[1];
    is->get(s,8);
    *is >> fDy[0];
    is->get(s,8);
    *is >> fDy[1];
    is->get(s,5);
    *is >> fDz;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTrapezoid2Data &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os      The output stream
    //    AliITSBoxData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTrapezoid2Data &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSPGonData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSTrapezoidData)
//______________________________________________________________________
void AliITSTrapezoidData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fTheta=" << fTheta << " fPhi=" << fPhi << " fDz=" << fDz;
    *os << " fH[0]=" << fH[0]<< " fH[1]=" << fH[1];
    *os << " fBl[0]=" << fBl[0] << " fBl[1]=" << fBl[1];
    *os << " fTl[0]=" << fTl[0] << " fTl[1]=" << fTl[1];
    *os << " fAlp[0]=" << fAlp[0] << " fAlp[1]=" << fAlp[1];
    *os << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSTrapezoidData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);
    is->get(s,6);
    *is >> fTheta;
    is->get(s,6);
    *is >> fPhi;
    is->get(s,5);
    *is >> fDz;
    is->get(s,7);
    *is >> fH[0];
    is->get(s,7);
    *is >> fH[1];
    is->get(s,8);
    *is >> fBl[0];
    is->get(s,8);
    *is >> fBl[1];
    is->get(s,8);
    *is >> fTl[0];
    is->get(s,8);
    *is >> fTl[1];
    is->get(s,9);
    *is >> fAlp[0];
    is->get(s,9);
    *is >> fAlp[1];
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTrapezoidData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os      The output stream
    //    AliITSBoxData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTrapezoidData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSPGonData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSTrapezoidTwistedData)
//______________________________________________________________________
void AliITSTrapezoidTwistedData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fTheta=" << fTheta << " fPhi=" << fPhi << " fDz=" << fDz;
    *os << " fTwist=" << fTwist;
    *os << " fH[0]=" << fH[0]<< " fH[1]=" << fH[1];
    *os << " fBl[0]=" << fBl[0] << " fBl[1]=" << fBl[1];
    *os << " fTl[0]=" << fTl[0] << " fTl[1]=" << fTl[1];
    *os << " fAlp[0]=" << fAlp[0] << " fAlp[1]=" << fAlp[1];
    *os << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSTrapezoidTwistedData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);
    is->get(s,6);
    *is >> fTheta;
    is->get(s,6);
    *is >> fPhi;
    is->get(s,5);
    *is >> fDz;
    is->get(s,8);
    *is >> fTwist;
    is->get(s,7);
    *is >> fH[0];
    is->get(s,7);
    *is >> fH[1];
    is->get(s,8);
    *is >> fBl[0];
    is->get(s,8);
    *is >> fBl[1];
    is->get(s,8);
    *is >> fTl[0];
    is->get(s,8);
    *is >> fTl[1];
    is->get(s,9);
    *is >> fAlp[0];
    is->get(s,9);
    *is >> fAlp[1];
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTrapezoidTwistedData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os      The output stream
    //    AliITSBoxData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTrapezoidTwistedData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSPGonData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSTubeData)
//______________________________________________________________________
void AliITSTubeData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os <<"       Z        ,      Rmin      ,      Rmax      " << endl;
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << setprecision(16) << fDz <<"\t";
    *os << setprecision(16) << fRmin << "\t";
    *os << setprecision(16) << fRmax << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}

//______________________________________________________________________
void AliITSTubeData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);

    is->getline(s,49);
        *is >> fDz >> fRmin >> fRmax;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTubeData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os       The output stream
    //    AliITSTubeData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTubeData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSTubeData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSTubeSegData)
//______________________________________________________________________
void AliITSTubeSegData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fPhi0=" << fPhi0 << " fPhi1=" << fPhi1 << endl;
    *os <<"       Z        ,      Rmin      ,      Rmax      " << endl;
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << setprecision(16) << fDz <<"\t";
    *os << setprecision(16) << fRmin << "\t";
    *os << setprecision(16) << fRmax << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSTubeSegData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);

    is->get(s,6);
    *is >> fPhi0;
    is->get(s,7);
    *is >> fPhi1;
    is->getline(s,49);
	*is >> fDz >> fRmin >> fRmax;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTubeSegData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os       The output stream
    //    AliITSTubeData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTubeSegData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSTubeData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSTubeCutData)
//______________________________________________________________________
void AliITSTubeCutData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fPhi0=" << fPhi0 << " fPhi1=" << fPhi1;
    *os << " Norm0=("<<(fNorm[0])[0]<<" "<<(fNorm[0])[1]<<" "<<(fNorm[0])[2];
    *os << ") Norm1=("<<(fNorm[1])[0]<<" "<<(fNorm[1])[1]<<" "<<(fNorm[1])[2];
    *os << ")"<< endl;
    *os <<"       Z        ,      Rmin      ,      Rmax      " << endl;
    *os << setprecision(16) << fDz <<"\t";
    *os << setprecision(16) << fRmin << "\t";
    *os << setprecision(16) << fRmax << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSTubeCutData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);

    is->get(s,6);
    *is >> fPhi0;
    is->get(s,7);
    *is >> fPhi1;
    is->get(s,8);
    *is >> (fNorm[0])[0]>>(fNorm[0])[1]>>(fNorm[0])[2];
    is->get(s,9);
    *is >> (fNorm[1])[0]>>(fNorm[1])[1]>>(fNorm[1])[2];
    is->get(s,1);
    is->getline(s,49);
    *is >> fDz >> fRmin >> fRmax;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTubeCutData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os       The output stream
    //    AliITSTubeData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTubeCutData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSTubeData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}

//======================================================================
ClassImp(AliITSTubeEllipticalData)
//______________________________________________________________________
void AliITSTubeEllipticalData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os <<"       Z        ,  Semi-axis-x      ,  Semi-axis-y      " << endl;
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << setprecision(16) << fDz <<"\t";
    *os << setprecision(16) << fP0 << "\t";
    *os << setprecision(16) << fP1 << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}

//______________________________________________________________________
void AliITSTubeEllipticalData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);

    is->getline(s,49);
    *is >> fDz >> fP0 >> fP1;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTubeEllipticalData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os       The output stream
    //    AliITSTubeData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTubeEllipticalData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSTubeData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSTubeHyperbolicData)
//______________________________________________________________________
void AliITSTubeHyperbolicData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os <<"       Z               Rmin             Rmax         Theta"<<endl;
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << setprecision(16) << fDz <<"\t";
    *os << setprecision(16) << fRmin << "\t";
    *os << setprecision(16) << fRmax << "\t";
    *os << setprecision(16) << fTheta << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}

//______________________________________________________________________
void AliITSTubeHyperbolicData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);

    is->getline(s,49);
    *is >> fDz >> fRmin >> fRmax >> fTheta;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSTubeHyperbolicData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os       The output stream
    //    AliITSTubeData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSTubeHyperbolicData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSTubeData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSConeData)
//______________________________________________________________________
void AliITSConeData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os <<"       Z               Rmin             Rmax" << endl;
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << setprecision(16) << fDz <<"\t";
    *os << setprecision(16) << fRmin0 << "\t";
    *os << setprecision(16) << fRmax0 << endl;
    *os << setprecision(16) << fDz <<"\t";
    *os << setprecision(16) << fRmin1 << "\t";
    *os << setprecision(16) << fRmax1 << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSConeData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);

    is->getline(s,49);
    *is >> fDz >> fRmin0 >> fRmax0;
    *is >> fDz >> fRmin1 >> fRmax1;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSConeData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os       The output stream
    //    AliITSTubeData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSConeData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSTubeData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSConeSegData)
//______________________________________________________________________
void AliITSConeSegData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fPhi0=" << fPhi0 << " fPhi1=" << fPhi1 << endl;
    *os <<"       Z        ,      Rmin      ,      Rmax      " << endl;
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << setprecision(16) << fDz <<"\t";
    *os << setprecision(16) << fRmin0 << "\t";
    *os << setprecision(16) << fRmax0 << endl;
    *os << setprecision(16) << fDz <<"\t";
    *os << setprecision(16) << fRmin1 << "\t";
    *os << setprecision(16) << fRmax1 << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSConeSegData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);
    is->get(s,6);
    *is >> fPhi0;
    is->get(s,7);
    *is >> fPhi1;
    is->getline(s,49);
    *is >> fDz >> fRmin0 >> fRmax0;
    *is >> fDz >> fRmin1 >> fRmax1;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSConeSegData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os        The output stream
    //    AliITSConeSegData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSConeSegData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is        The input stream
    //    AliITSConeSegData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSPConeData)
//______________________________________________________________________
void AliITSPConeData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i;

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fNz=" << fNz << " fPhi0=" << fPhi0 << " fdPhi=" << fDphi << endl;
    *os <<"       Z        ,      Rmin      ,      Rmax      " << endl;
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    for(i=0;i<fNz;i++){
	*os << setprecision(16) << fZ[i] <<"\t";
	*os << setprecision(16) << fRmin[i] << "\t";
	*os << setprecision(16) << fRmax[i] << endl;
    } // end for i
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSPConeData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t i;
    char s[50];

    AliITSBaseVolParams::Read(is);
    is->get(s,4);
    *is >> fNz;
    is->get(s,6);
    *is >> fPhi0;
    is->get(s,6);
    *is >> fDphi;
    is->getline(s,49);
    Size(fNz);
    for(i=0;i<fNz;i++){
	*is >> fZ[i] >> fRmin[i] >> fRmax[i];
    } // end for i
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSPConeData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os        The output stream
    //    AliITSPConeData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSPConeData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is        The input stream
    //    AliITSPConeData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSSphereData)
//______________________________________________________________________
void AliITSSphereData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fTheta[0]=" << fTheta[0] << " fTheta[1]=" << fTheta[1] << endl;
    *os << "fPhi[0]=" << fPhi[0] << " fPhi[1]=" << fPhi[1] << endl;
    *os <<"      Rmin      ,      Rmax      " << endl;
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << setprecision(16) << fRmin << "\t";
    *os << setprecision(16) << fRmax << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSSphereData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);
    is->get(s,10);
    *is >> fTheta[0];
    is->get(s,11);
    *is >> fTheta[1];
    is->get(s,8);
    *is >> fPhi[0];
    is->get(s,9);
    *is >> fPhi[1];
    is->getline(s,49);
    *is >>fRmin >> fRmax;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSSphereData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os        The output stream
    //    AliITSPConeData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSSphereData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is        The input stream
    //    AliITSPConeData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSParallelpipedData)
//______________________________________________________________________
void AliITSParallelpipedData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fDx=" << fDx << " fDy=" << fDy << " fDz=" << fDz << endl;
    *os << "fAlpha=" << fAlpha << " fTheta=" << fTheta <<" fPhi="<<fPhi<<endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSParallelpipedData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    char s[50];

    AliITSBaseVolParams::Read(is);
    is->get(s,4);
    *is >> fDx;
    is->get(s,5);
    *is >> fDy;
    is->get(s,5);
    *is >> fDz;
    is->get(s,7);
    *is >> fAlpha;
    is->get(s,8);
    *is >> fTheta;
    is->get(s,6);
    *is >> fPhi;
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSParallelpipedData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os      The output stream
    //    AliITSBoxData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSParallelpipedData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is      The input stream
    //    AliITSBoxData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
//======================================================================
ClassImp(AliITSPGonData)
//______________________________________________________________________
void AliITSPGonData::Print(ostream *os){
    // Prints out the data kept in this class
    // Inputs:
    //    ostream *os The output stream pointer
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i;

#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    AliITSBaseVolParams::Print(os);
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "fNz=" << fNz << " fNphi=" << fNphi << " fPhi0=" << fPhi0;
    *os << " fdPhi=" << fDphi << endl;
    *os <<"       Z        ,      Rmin      ,      Rmax      " << endl;
    for(i=0;i<fNz;i++){
	*os << setprecision(16) << fZ[i] <<"\t";
	*os << setprecision(16) << fRmin[i] << "\t";
	*os << setprecision(16) << fRmax[i] << endl;
    } // end for i
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSPGonData::Read(istream *is){
    // Read in data kept in this class
    // Inputs:
    //   istream *is  the input stream
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t i;
    char s[50];

    AliITSBaseVolParams::Read(is);
    
    is->get(s,4);
    *is >> fNz;
    is->get(s,6);
    *is >> fNphi;
    is->get(s,6);
    *is >> fPhi0;
    is->get(s,6);
    *is >> fDphi;
    is->getline(s,49);

    Size(fNz);
    for(i=0;i<fNz;i++){
	*is >> fZ[i] >> fRmin[i] >> fRmax[i];
    } // end for i
}
//______________________________________________________________________
ostream &operator<<(ostream &os,AliITSPGonData &p){
    // Operator << for C++ like output
    // Inputs:
    //    ostream &os       The output stream
    //    AliITSPGonData &p The class to be outputed
    // Output:
    //    none.
    // Return:
    //    ostream &os        The output stream

    p.Print(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSPGonData &r){
    // Operator << for C++ like output
    // Inputs:
    //    istream &is       The input stream
    //    AliITSPGonData &r The class to be read in
    // Output:
    //    none.
    // Return:
    //    istream &is        The input stream

    r.Read(&is);
    return is;
}
