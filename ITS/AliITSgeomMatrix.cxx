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
Revision 1.7  2001/02/03 00:00:30  nilsen
New version of AliITSgeom and related files. Now uses automatic streamers,
set up for new formatted .det file which includes detector information.
Additional smaller modifications are still to come.

Revision 1.5  2000/10/02 16:32:35  barbera
Forward declaration added

Revision 1.1.2.6  2000/10/02 15:52:05  barbera
Forward declaration added

Revision 1.4  2000/09/07 17:30:45  nilsen
fixed a bug in SixAnglesFromMatrix.

Revision 1.3  2000/09/05 14:25:50  nilsen
Made fixes for HP compiler. All function parameter default values placed
in .h file. Fixed the usual problem with HP comilers and the "for(Int_t i..."
business. Replaced casting (Double_t [3][3]) to (Double_t (*)[3]) for HP.
Lastly removed all "const" before function parameters which were 2 dim. arrays,
because on HP root generates some strange code (?). Thanks Peter for the
changes.

Revision 1.2  2000/08/29 20:16:50  nilsen
New class for ITS coordiante transformations used by AliITSgeom nearly
exclusively.

Revision 1.1.2.1  2000/06/04 16:32:31  Nilsen
A new class to hold the matrix information needed by AliITSgeom.

*/
#include <iostream.h>
#include <iomanip.h>
#include <TMath.h>
#include <TBuffer.h>

#include "AliITSgeomMatrix.h"

ClassImp(AliITSgeomMatrix)
//----------------------------------------------------------------------
AliITSgeomMatrix::AliITSgeomMatrix(){
////////////////////////////////////////////////////////////////////////
// The Default constructor for the AliITSgeomMatrix class. By Default
// the angles of rotations are set to zero, meaning that the rotation
// matrix is the unit matrix. The translation vector is also set to zero
// as are the module id number. The detector type is set to -1 (an undefined
// value). The full rotation matrix is kept so that the evaluation 
// of a coordinate transformation can be done quickly and with a minimum
// of CPU overhead. The basic coordinate systems are the ALICE global
// coordinate system and the detector local coordinate system. In general
// this structure is not limited to just those two coordinate systems.
//Begin_Html
/*
<img src="picts/ITS/AliISgeomMatrix_L1.gif">
*/
//End_Html
////////////////////////////////////////////////////////////////////////
    Int_t i,j;

    fDetectorIndex = -1; // a value never defined.
    for(i=0;i<3;i++){
	fid[i] = 0;
	frot[i] = ftran[i] = 0.0;
	for(j=0;j<3;j++) fm[i][j] = 0.0;
    }// end for i
    fm[0][0] = fm[1][1] = fm[2][2] = 1.0;
}
//----------------------------------------------------------------------
AliITSgeomMatrix::AliITSgeomMatrix(const AliITSgeomMatrix &sourse){
////////////////////////////////////////////////////////////////////////
// The standard copy constructor. This make a full / proper copy of
// this class.
////////////////////////////////////////////////////////////////////////
	Int_t i,j;

	this->fDetectorIndex = sourse.fDetectorIndex;
	for(i=0;i<3;i++){
		this->fid[i]     = sourse.fid[i];
		this->frot[i]    = sourse.frot[i];
		this->ftran[i]   = sourse.ftran[i];
		for(j=0;j<3;j++) this->fm[i][j] = sourse.fm[i][j];
	}// end for i
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::operator=(const AliITSgeomMatrix &sourse){
////////////////////////////////////////////////////////////////////////
// The standard = operator. This make a full / proper copy of
// this class.
////////////////////////////////////////////////////////////////////////
	Int_t i,j;

	this->fDetectorIndex = sourse.fDetectorIndex;
	for(i=0;i<3;i++){
		this->fid[i]     = sourse.fid[i];
		this->frot[i]    = sourse.frot[i];
		this->ftran[i]   = sourse.ftran[i];
		for(j=0;j<3;j++) this->fm[i][j] = sourse.fm[i][j];
	}// end for i
}
//----------------------------------------------------------------------
AliITSgeomMatrix::AliITSgeomMatrix(const Int_t idt,const Int_t id[3],
		   const Double_t rot[3],const Double_t tran[3]){
////////////////////////////////////////////////////////////////////////
// This is a constructor for the AliITSgeomMatrix class. The matrix is
// defined by 3 standard rotation angles [radians], and the translation
// vector tran [cm]. In addition the layer, ladder, and detector number
// for this particular module and the type of module must be given.
// The full rotation matrix is kept so that the evaluation 
// of a coordinate transformation can be done quickly and with a minimum
// of CPU overhead. The basic coordinate systems are the ALICE global
// coordinate system and the detector local coordinate system. In general
// this structure is not limited to just those two coordinate systems.
//Begin_Html
/*
<img src="picts/ITS/AliISgeomMatrix_L1.gif">
*/
//End_Html
////////////////////////////////////////////////////////////////////////
	Int_t i;

	fDetectorIndex = idt; // a value never defined.
	for(i=0;i<3;i++){
		fid[i]   = id[i];
		frot[i]  = rot[i];
		ftran[i] = tran[i];
	}// end for i
	this->MatrixFromAngle();
}
//----------------------------------------------------------------------
AliITSgeomMatrix::AliITSgeomMatrix(const Int_t idt, const Int_t id[3],
                                   Double_t matrix[3][3],
                                   const Double_t tran[3]){
////////////////////////////////////////////////////////////////////////
// This is a constructor for the AliITSgeomMatrix class. The rotation matrix
// is given as one of the inputs, and the translation vector tran [cm]. In 
// addition the layer, ladder, and detector number for this particular
// module and the type of module must be given. The full rotation matrix
// is kept so that the evaluation of a coordinate transformation can be
// done quickly and with a minimum of CPU overhead. The basic coordinate
// systems are the ALICE global coordinate system and the detector local
// coordinate system. In general this structure is not limited to just
// those two coordinate systems.
//Begin_Html
/*
<img src="picts/ITS/AliISgeomMatrix_L1.gif">
*/
//End_Html
////////////////////////////////////////////////////////////////////////
	Int_t i,j;

	fDetectorIndex = idt; // a value never defined.
	for(i=0;i<3;i++){
		fid[i]   = id[i];
		ftran[i] = tran[i];
		for(j=0;j<3;j++) fm[i][j] = matrix[i][j];
	}// end for i
	this->AngleFromMatrix();
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::SixAnglesFromMatrix(Double_t *ang){
////////////////////////////////////////////////////////////////////////
// This function returns the 6 GEANT 3.21 rotation angles [degrees] in
// the array ang which must be at least [6] long.
////////////////////////////////////////////////////////////////////////
    Double_t si,c=180./TMath::Pi();

    ang[1] = TMath::ATan2(fm[0][1],fm[0][0]);
    if(TMath::Cos(ang[1])!=0.0) si = fm[0][0]/TMath::Cos(ang[1]);
    else si = fm[0][1]/TMath::Sin(ang[1]);
    ang[0] = TMath::ATan2(si,fm[0][2]);

    ang[3] = TMath::ATan2(fm[1][1],fm[1][0]);
    if(TMath::Cos(ang[3])!=0.0) si = fm[1][0]/TMath::Cos(ang[3]);
    else si = fm[1][1]/TMath::Sin(ang[3]);
    ang[2] = TMath::ATan2(si,fm[1][2]);

    ang[5] = TMath::ATan2(fm[2][1],fm[2][0]);
    if(TMath::Cos(ang[5])!=0.0) si = fm[2][0]/TMath::Cos(ang[5]);
    else si = fm[2][1]/TMath::Sin(ang[5]);
    ang[4] = TMath::ATan2(si,fm[2][2]);

    for(Int_t i=0;i<6;i++) {ang[i] *= c; if(ang[i]<0.0) ang[i] += 360.;}
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::MatrixFromSixAngles(const Double_t *ang){
////////////////////////////////////////////////////////////////////////
// Given the 6 GEANT 3.21 rotation angles [degree], this will compute and
// set the rotations matrix and 3 standard rotation angles [radians].
// These angles and rotation matrix are overwrite the existing values in
// this class.
////////////////////////////////////////////////////////////////////////
    Int_t    i,j;
    Double_t si,lr[9],c=TMath::Pi()/180.;

    si    = TMath::Sin(c*ang[0]);
    if(ang[0]== 90.0)                 si = +1.0;
    if(ang[0]==270.0)                 si = -1.0;
    if(ang[0]==  0.0||ang[0]==180.) si =  0.0;
    lr[0] = si * TMath::Cos(c*ang[1]);
    lr[1] = si * TMath::Sin(c*ang[1]);
    lr[2] = TMath::Cos(c*ang[0]);
    if(ang[0]== 90.0||ang[0]==270.) lr[2] =  0.0;
    if(ang[0]== 0.0)                  lr[2] = +1.0;
    if(ang[0]==180.0)                 lr[2] = -1.0;
//
    si    =  TMath::Sin(c*ang[2]);
    if(ang[2]== 90.0)                 si = +1.0; 
    if(ang[2]==270.0)                 si = -1.0;
    if(ang[2]==  0.0||ang[2]==180.) si =  0.0;
    lr[3] = si * TMath::Cos(c*ang[3]);
    lr[4] = si * TMath::Sin(c*ang[3]);
    lr[5] = TMath::Cos(c*ang[2]);
    if(ang[2]== 90.0||ang[2]==270.) lr[5] =  0.0;
    if(ang[2]==  0.0)                 lr[5] = +1.0;
    if(ang[2]==180.0)                 lr[5] = -1.0;
//
    si    = TMath::Sin(c*ang[4]);
    if(ang[4]== 90.0)                 si = +1.0;
    if(ang[4]==270.0)                 si = -1.0;
    if(ang[4]==  0.0||ang[4]==180.) si =  0.0;
    lr[6] = si * TMath::Cos(c*ang[5]);
    lr[7] = si * TMath::Sin(c*ang[5]);
    lr[8] = TMath::Cos(c*ang[4]);
    if(ang[4]== 90.0||ang[4]==270.0) lr[8] =  0.0;
    if(ang[4]==  0.0)                  lr[8] = +1.0;
    if(ang[4]==180.0)                  lr[8] = -1.0;
    // Normalize these elements and fill matrix fm.
    for(i=0;i<3;i++){// reuse si.
	si = 0.0;
	for(j=0;j<3;j++) si += lr[3*i+j]*lr[3*i+j];
	si = TMath::Sqrt(1./si);
	for(j=0;j<3;j++) fm[i][j] = si*lr[3*i+j];
    } // end for i
    this->AngleFromMatrix();
}
//----------------------------------------------------------------------
AliITSgeomMatrix::AliITSgeomMatrix(const Double_t rotd[6]/*degrees*/,
                                   const Int_t idt,const Int_t id[3],
		                   const Double_t tran[3]){
////////////////////////////////////////////////////////////////////////
// This is a constructor for the AliITSgeomMatrix class. The matrix is
// defined by the 6 GEANT 3.21 rotation angles [degrees], and the translation
// vector tran [cm]. In addition the layer, ladder, and detector number
// for this particular module and the type of module must be given.
// The full rotation matrix is kept so that the evaluation 
// of a coordinate transformation can be done quickly and with a minimum
// of CPU overhead. The basic coordinate systems are the ALICE global
// coordinate system and the detector local coordinate system. In general
// this structure is not limited to just those two coordinate systems.
//Begin_Html
/*
<img src="picts/ITS/AliISgeomMatrix_L1.gif">
*/
//End_Html
////////////////////////////////////////////////////////////////////////
    Int_t i;

    fDetectorIndex = idt; // a value never defined.
    for(i=0;i<3;i++){
	fid[i]   = id[i];
	ftran[i] = tran[i];
    }// end for i
    this->MatrixFromSixAngles(rotd);
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::AngleFromMatrix(){
////////////////////////////////////////////////////////////////////////
// Computes the angles from the rotation matrix up to a phase of 180 degrees.
////////////////////////////////////////////////////////////////////////
    Double_t rx,ry,rz;
    // get angles from matrix up to a phase of 180 degrees.

    rx = TMath::ATan2(fm[2][1],fm[2][2]);if(rx<0.0) rx += 2.0*TMath::Pi();
    ry = TMath::ASin(fm[0][2]);          if(ry<0.0) ry += 2.0*TMath::Pi();
    rz = TMath::ATan2(fm[1][1],fm[0][0]);if(rz<0.0) rz += 2.0*TMath::Pi();
    frot[0] = rx;
    frot[1] = ry;
    frot[2] = rz;
    return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::MatrixFromAngle(){
////////////////////////////////////////////////////////////////////////
// Computes the Rotation matrix from the angles [radians] kept in this
// class.
////////////////////////////////////////////////////////////////////////
   Double_t sx,sy,sz,cx,cy,cz;

   sx = TMath::Sin(frot[0]); cx = TMath::Cos(frot[0]);
   sy = TMath::Sin(frot[1]); cy = TMath::Cos(frot[1]);
   sz = TMath::Sin(frot[2]); cz = TMath::Cos(frot[2]);
   fm[0][0] =  cz*cy;             // fr[0]
   fm[0][1] = -cz*sy*sx - sz*cx;  // fr[1]
   fm[0][2] = -cz*sy*cx + sz*sx;  // fr[2]
   fm[1][0] =  sz*cy;             // fr[3]
   fm[1][1] = -sz*sy*sx + cz*cx;  // fr[4]
   fm[1][2] = -sz*sy*cx - cz*sx;  // fr[5]
   fm[2][0] =  sy;                // fr[6]
   fm[2][1] =  cy*sx;             // fr[7]
   fm[2][2] =  cy*cx;             // fr[8]

}
//----------------------------------------------------------------------
void AliITSgeomMatrix::GtoLPosition(const Double_t g0[3],Double_t l[3]){
////////////////////////////////////////////////////////////////////////
// Returns the local coordinates given the global coordinates [cm].
////////////////////////////////////////////////////////////////////////
	Int_t    i,j;
	Double_t g[3];

	for(i=0;i<3;i++) g[i] = g0[i] - ftran[i];
	for(i=0;i<3;i++){
		l[i] = 0.0;
		for(j=0;j<3;j++) l[i] += fm[i][j]*g[j];
		// g = R l + translation
	} // end for i
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGPosition(const Double_t l[3],Double_t g[3]){
////////////////////////////////////////////////////////////////////////
// Returns the global coordinates given the local coordinates [cm].
////////////////////////////////////////////////////////////////////////
	Int_t    i,j;

	for(i=0;i<3;i++){
		g[i] = 0.0;
		for(j=0;j<3;j++) g[i] += fm[j][i]*l[j];
		g[i] += ftran[i];
		// g = R^t l + translation
	} // end for i
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::GtoLMomentum(const Double_t g[3],Double_t l[3]){
////////////////////////////////////////////////////////////////////////
// Returns the local coordinates of the momentum given the global
// coordinates of the momentum. It transforms just like GtoLPosition
// except that the translation vector is zero.
////////////////////////////////////////////////////////////////////////
	Int_t    i,j;

	for(i=0;i<3;i++){
		l[i] = 0.0;
		for(j=0;j<3;j++) l[i] += fm[i][j]*g[j];
		// g = R l
	} // end for i
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGMomentum(const Double_t l[3],Double_t g[3]){
////////////////////////////////////////////////////////////////////////
// Returns the Global coordinates of the momentum given the local
// coordinates of the momentum. It transforms just like LtoGPosition
// except that the translation vector is zero.
////////////////////////////////////////////////////////////////////////
	Int_t    i,j;

	for(i=0;i<3;i++){
		g[i] = 0.0;
		for(j=0;j<3;j++) g[i] += fm[j][i]*l[j];
		// g = R^t l
	} // end for i
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::GtoLPositionError(Double_t g[3][3],
                                               Double_t l[3][3]){
////////////////////////////////////////////////////////////////////////
// Given an Uncertainty matrix in Global coordinates it is rotated so that 
// its representation in local coordinates can be returned. There is no
// effect due to the translation vector or its uncertainty.
////////////////////////////////////////////////////////////////////////
	Int_t    i,j,k,m;

	for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)for(m=0;m<3;m++)
		l[i][m] = fm[j][i]*g[j][k]*fm[k][m];
		// g = R^t l R
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGPositionError(Double_t l[3][3],
                                               Double_t g[3][3]){
////////////////////////////////////////////////////////////////////////
// Given an Uncertainty matrix in Local coordinates it is rotated so that 
// its representation in global coordinates can be returned. There is no
// effect due to the translation vector or its uncertainty.
////////////////////////////////////////////////////////////////////////
	Int_t    i,j,k,m;

	for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)for(m=0;m<3;m++)
		g[i][m] = fm[i][j]*l[j][k]*fm[m][k];
		// g = R l R^t
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::GtoLPositionTracking(const Double_t g0[3],
					    Double_t l[3]){
////////////////////////////////////////////////////////////////////////
// A slightly different coordinate system is used when tracking.
// This coordinate system is only relevant when the geometry represents
// the cylindrical ALICE ITS geometry. For tracking the Z axis is left
// alone but X -> -Y and Y -> X such that X always points out of the
// ITS Cylinder for every layer including layer 1 (where the detector 
// are mounted upside down).
//Begin_Html
/*
<img src="picts/ITS/AliITSgeomMatrix_T1.gif">
 */
//End_Html
////////////////////////////////////////////////////////////////////////
    Double_t l0[3];

    this->GtoLPosition(g0,l0);
    if(fid[0]==1){ // for layer 1 the detector are flipped upside down
	           // with respect to the others.
	l[0] = +l0[1];
	l[1] = -l0[0];
	l[2] = +l0[2];
    }else{
	l[0] = -l0[1];
	l[1] = +l0[0];
	l[2] = +l0[2];
    } // end if
    return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGPositionTracking(const Double_t l[3],
					    Double_t g[3]){
////////////////////////////////////////////////////////////////////////
// A slightly different coordinate system is used when tracking.
// This coordinate system is only relevant when the geometry represents
// the cylindrical ALICE ITS geometry. For tracking the Z axis is left
// alone but X -> -Y and Y -> X such that X always points out of the
// ITS Cylinder for every layer including layer 1 (where the detector 
// are mounted upside down).
//Begin_Html
/*
<img src="picts/ITS/AliITSgeomMatrix_T1.gif">
 */
//End_Html
////////////////////////////////////////////////////////////////////////
    Double_t l0[3];

    if(fid[0]==1){ // for layer 1 the detector are flipped upside down
	           // with respect to the others.
	l0[0] = -l[1];
	l0[1] = +l[0];
	l0[2] = +l[2];
    }else{
	l0[0] = +l[1];
	l0[1] = -l[0];
	l0[2] = +l[2];
    } // end if
    this->LtoGPosition(l0,g);
    return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::GtoLMomentumTracking(const Double_t g[3],
					    Double_t l[3]){
////////////////////////////////////////////////////////////////////////
// A slightly different coordinate system is used when tracking.
// This coordinate system is only relevant when the geometry represents
// the cylindrical ALICE ITS geometry. For tracking the Z axis is left
// alone but X -> -Y and Y -> X such that X always points out of the
// ITS Cylinder for every layer including layer 1 (where the detector 
// are mounted upside down).
//Begin_Html
/*
<img src="picts/ITS/AliITSgeomMatrix_T1.gif">
 */
//End_Html
////////////////////////////////////////////////////////////////////////
    Double_t l0[3];

    this->GtoLMomentum(g,l0);
    if(fid[0]==1){ // for layer 1 the detector are flipped upside down
	           // with respect to the others.
	l[0] = +l0[1];
	l[1] = -l0[0];
	l[2] = +l0[2];
    }else{
	l[0] = -l0[1];
	l[1] = +l0[0];
	l[2] = +l0[2];
    } // end if
    return;
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGMomentumTracking(const Double_t l[3],
					    Double_t g[3]){
////////////////////////////////////////////////////////////////////////
// A slightly different coordinate system is used when tracking.
// This coordinate system is only relevant when the geometry represents
// the cylindrical ALICE ITS geometry. For tracking the Z axis is left
// alone but X -> -Y and Y -> X such that X always points out of the
// ITS Cylinder for every layer including layer 1 (where the detector 
// are mounted upside down).
//Begin_Html
/*
<img src="picts/ITS/AliITSgeomMatrix_T1.gif">
 */
//End_Html
////////////////////////////////////////////////////////////////////////
    Double_t l0[3];

    if(fid[0]==1){ // for layer 1 the detector are flipped upside down
	           // with respect to the others.
	l0[0] = -l[1];
	l0[1] = +l[0];
	l0[2] = +l[2];
    }else{
	l0[0] = +l[1];
	l0[1] = -l[0];
	l0[2] = +l[2];
    } // end if
    this->LtoGMomentum(l0,g);
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::GtoLPositionErrorTracking(Double_t g[3][3],
						 Double_t l[3][3]){
////////////////////////////////////////////////////////////////////////
// A slightly different coordinate system is used when tracking.
// This coordinate system is only relevant when the geometry represents
// the cylindrical ALICE ITS geometry. For tracking the Z axis is left
// alone but X -> -Y and Y -> X such that X always points out of the
// ITS Cylinder for every layer including layer 1 (where the detector 
// are mounted upside down).
//Begin_Html
/*
<img src="picts/ITS/AliITSgeomMatrix_T1.gif">
 */
//End_Html
////////////////////////////////////////////////////////////////////////
	Int_t    i,j,k,m;
	Double_t Rt[3][3];
	Double_t A0[3][3] = {{0.,+1.,0.},{-1.,0.,0.},{0.,0.,+1.}};
	Double_t A1[3][3] = {{0.,-1.,0.},{+1.,0.,0.},{0.,0.,+1.}};

	if(fid[0]==1) for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)
	    Rt[i][k] = A0[i][j]*fm[j][k];
	else for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)
	    Rt[i][k] = A1[i][j]*fm[j][k];
	for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)for(m=0;m<3;m++)
		l[i][m] = Rt[j][i]*g[j][k]*Rt[k][m];
		// g = R^t l R
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGPositionErrorTracking(Double_t l[3][3],
						 Double_t g[3][3]){
////////////////////////////////////////////////////////////////////////
// A slightly different coordinate system is used when tracking.
// This coordinate system is only relevant when the geometry represents
// the cylindrical ALICE ITS geometry. For tracking the Z axis is left
// alone but X -> -Y and Y -> X such that X always points out of the
// ITS Cylinder for every layer including layer 1 (where the detector 
// are mounted upside down).
//Begin_Html
/*
<img src="picts/ITS/AliITSgeomMatrix_T1.gif">
 */
//End_Html
////////////////////////////////////////////////////////////////////////
	Int_t    i,j,k,m;
	Double_t Rt[3][3];
	Double_t A0[3][3] = {{0.,+1.,0.},{-1.,0.,0.},{0.,0.,+1.}};
	Double_t A1[3][3] = {{0.,-1.,0.},{+1.,0.,0.},{0.,0.,+1.}};

	if(fid[0]==1) for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)
	    Rt[i][k] = A0[i][j]*fm[j][k];
	else for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)
	    Rt[i][k] = A1[i][j]*fm[j][k];
	for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)for(m=0;m<3;m++)
		g[i][m] = Rt[i][j]*l[j][k]*Rt[m][k];
		// g = R l R^t
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::PrintTitles(ostream *os){
////////////////////////////////////////////////////////////////////////
// Standard output format for this class but it includes variable
// names and formatting that makes it easer to read.
////////////////////////////////////////////////////////////////////////
    Int_t i,j;

    *os << "fDetectorIndex=" << fDetectorIndex << " fid[3]={";
    for(i=0;i<3;i++) *os << fid[i]   << " ";
    *os << "} frot[3]={";
    for(i=0;i<3;i++) *os << frot[i]  << " ";
    *os << "} ftran[3]={";
    for(i=0;i<3;i++) *os << ftran[i] << " ";
    *os << "} fm[3][3]={";
    for(i=0;i<3;i++){for(j=0;j<3;j++){  *os << fm[i][j] << " ";} *os <<"}{";}
    *os << "}" << endl;
    return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::PrintComment(ostream *os){
////////////////////////////////////////////////////////////////////////
//  output format used by Print..
////////////////////////////////////////////////////////////////////////
    *os << "fDetectorIndex fid[0] fid[1] fid[2] ftran[0] ftran[1] ftran[2] ";
    *os << "fm[0][0]  fm[0][1]  fm[0][2]  fm[1][0]  fm[1][1]  fm[1][2]  ";
    *os << "fm[2][0]  fm[2][1]  fm[2][2] ";
    return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::Print(ostream *os){
////////////////////////////////////////////////////////////////////////
// Standard output format for this class.
////////////////////////////////////////////////////////////////////////
    Int_t i,j;
    Int_t fmt;

    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << fDetectorIndex << " ";
    for(i=0;i<3;i++) *os << fid[i]   << " ";
//    for(i=0;i<3;i++) *os << frot[i]  << " ";  // Redundant with fm[][].
    for(i=0;i<3;i++) *os << setprecision(16) << ftran[i] << " ";
    for(i=0;i<3;i++)for(j=0;j<3;j++)  *os << setprecision(16) << 
					  fm[i][j] << " ";
    *os << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::Read(istream *is){
////////////////////////////////////////////////////////////////////////
// Standard input format for this class.
////////////////////////////////////////////////////////////////////////
    Int_t i,j;

    *is >> fDetectorIndex;
    for(i=0;i<3;i++) *is >> fid[i];
//    for(i=0;i<3;i++) *is >> frot[i]; // Redundant with fm[][].
    for(i=0;i<3;i++) *is >> ftran[i];
    for(i=0;i<3;i++)for(j=0;j<3;j++)  *is >> fm[i][j];
    AngleFromMatrix(); // compute angles frot[].
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSgeomMatrix &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomMatrix &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//----------------------------------------------------------------------
