/* *************************************************************************
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
$Log:

*/

/* 
 *  Version: 0
 *  Written by Bjorn S. Nilsen with code stolen from Andreas Morsch's
 *  AliGeant3GeometryGUI class.
 */

#include <ctype.h>
#include "AliITSGeant3Geometry.h"

ClassImp(AliITSGeant3Geometry)

AliITSGeant3Geometry::AliITSGeant3Geometry(){
// Constructor

//  Store local copy of zebra bank entries
    TGeant3 *geant3 = (TGeant3*) gMC;
    if(geant3){
	fZlq    = geant3->Lq();
	fZq     = geant3->Q();
	fZiq    = geant3->Iq();
	fGclink = geant3->Gclink();
	fGcnum  = geant3->Gcnum();
	fGcvolu = geant3->Gcvolu();
    } // end if
}
//----------------------------------------------------------------------
Int_t AliITSGeant3Geometry::NChildren(Int_t idvol){
//
// Return number of children for volume idvol
    Int_t jvo = fZlq[fGclink->jvolum-idvol];
    Int_t nin = Int_t(fZq[jvo+3]);
    return nin;
}
//----------------------------------------------------------------------
Int_t AliITSGeant3Geometry::GetShape(Int_t idvol,Int_t &npar,Int_t &natt,
				     Float_t *par,Float_t *att){
// Returns the Geant shape number for a given volume
    Int_t ishape,i;

    Int_t jvo = fZlq[fGclink->jvolum-idvol];
//    Int_t nin = Int_t(fZq[jvo+3]);
    ishape = Int_t(fZq[jvo+2]);
    npar   = Int_t(fZq[jvo+5]);
    natt   = Int_t(fZq[jvo+6]);
    for(i=0;i<npar;i++) par[i] = fZq[jvo+7+i];
    for(i=0;i<natt;i++) att[i] = fZq[jvo+7+npar+i];
//    if(nin<0){ // devided volumes
//	printf("Divisions nin=%d\n",nin);
//    }else{
//	printf("SubVolumes nin=%d\n",nin);
//    } // end if
    return ishape;
}
//----------------------------------------------------------------------
Int_t AliITSGeant3Geometry::Child(Int_t idvol, Int_t idc){
//
// Return GEANT id of child number idc of volume idvol
    Int_t jvo = fZlq[fGclink->jvolum-idvol];
    Int_t nin=idc;
    Int_t jin = fZlq[jvo-nin];
    Int_t numb =  Int_t(fZq[jin +3]);
    if (numb > 1) {
	return -Int_t(fZq[jin+2]);
    } else {
	return Int_t(fZq[jin+2]);
    }
}
//----------------------------------------------------------------------
Int_t AliITSGeant3Geometry::Medium(Int_t idvol){
//
// Return medium number for volume idvol.
// If idvol is negative the volume results from a division.
    Int_t imed;
    if (idvol > 0) {
	Int_t jvo = fZlq[fGclink->jvolum-idvol];
	imed = Int_t(fZq[jvo+4]);
    } else {
	idvol=-idvol;
	Int_t jdiv = fZlq[fGclink->jvolum-idvol];
	Int_t ivin = Int_t ( fZq[jdiv+2]);
	Int_t jvin = fZlq[fGclink->jvolum-ivin];
	imed = Int_t (fZq[jvin+4]);
    }
    return imed;
}
//----------------------------------------------------------------------
Int_t AliITSGeant3Geometry::Material(Int_t idvol){
// Return material number for volume idvol.
// If idvol is negative the volume results from a division.

    Int_t imed=Medium(idvol);
    Int_t jtm  = fZlq[fGclink->jtmed-imed];
    return Int_t (fZq[jtm+6]);
}
//----------------------------------------------------------------------
void AliITSGeant3Geometry::GetGeometry(Int_t nlevels,Int_t *lnam,Int_t *lnum,
				       Double_t *xt,Double_t *r,Int_t &idshape,
				       Int_t &npar,Int_t &natt,Float_t *par,
				       Float_t *att,Int_t &imat,Int_t &imed){
//     Returns the MtoD/DtoM transformation for a given volume in the
//     volume tree LNAM for the given copy of those volumes in LNUM.
//     Input:
//           nlevels     the size of the arrays LNAM and LNUM used
//           lnam        the array of volume names through the tree of volumes
//           lnum        the copy number for the volumes above.
//     Output:
//           xt          Double precition array(3) for the coordiante
//                       translation for the DtoM/MtoD translation.
//           r           Double precition array(10) for the 3x3 rotation
//                       matrix for the DtoM/MtoD transformation. Element
//                       10 is a flag indicating if this matrix is the unit
//                       matrix
//           idshape     Geant3.21 volume shape index number
//           npar        Geant3.21 the number of shape parameters
//           natt        Geant3.21 the number of attributes
//           par         Geant3.21 an array of shape parameters
//           att         Geant3.21 an array of attributes
//           imat        Geant3.21 matrial index number for this volume
//           imed        Geant3.21 medium index number for this volume
//
//
    Int_t  ier,idvol,nlevel,i;
    TGeant3 *geant3 = (TGeant3*) gMC;
//
    fGcvolu->nlevel = 0;
    ier    = geant3->Glvolu(nlevels,lnam,lnum);
    nlevel = fGcvolu->nlevel;
    idvol  = fGcvolu->lvolum[nlevel-1];
    imed   = Medium(idvol);
    imat   = Material(idvol);
    idshape= GetShape(idvol,npar,natt,par,att);
    for(i=0;i<3;i++) xt[i] = fGcvolu->gtran[nlevel-1][i];
    for(i=0;i<10;i++) r[i] = fGcvolu->grmat[nlevel-1][i];
}
//----------------------------------------------------------------------
Int_t AliITSGeant3Geometry::StringToInt(char *name){
// converts up to a four letter char string to an equivelant int
    Int_t i,len = strlen(name);
    UInt_t iname=0;

    for(i=0;i<4;i++){
	if(i<len) iname += ( toupper(name[i]) << (8*(i)) );
	else iname += ( toupper(' ') << (8*(i)) );
    } // end for
    return (Int_t) iname;
}
