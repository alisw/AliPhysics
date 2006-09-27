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
$Id$ 
*/
////////////////////////////////////////////////////////////////////////
// This is the implementation file for AliITSgeomMatrix class. It 
// contains the routines to manipulate, setup, and queary the geometry 
// of a given ITS module. An ITS module may be one of at least three
// ITS detector technologies, Silicon Pixel, Drift, or Strip Detectors,
// and variations of these in size and/or layout. These routines let
// one go between ALICE global coordiantes (cm) to a given modules 
// specific local coordinates (cm).
////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>
#include <TBuffer.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyLine3D.h>
#include <TNode.h>
#include <TPCON.h>
#include <TBRIK.h>
#include <TXTRU.h>

#include "AliITSgeomMatrix.h"

ClassImp(AliITSgeomMatrix)
//----------------------------------------------------------------------
AliITSgeomMatrix::AliITSgeomMatrix():
TObject(),
fDetectorIndex(0), // Detector type index (like fShapeIndex was)
fid(),       // layer, ladder, detector numbers.
frot(),      //! vector of rotations about x,y,z [radians].
ftran(),     // Translation vector of module x,y,z.
fCylR(0.0),  //! R Translation in Cylinderical coordinates
fCylPhi(0.0),//! Phi Translation vector in Cylindrical coord.
fm(),        // Rotation matrix based on frot.
fPath(){     // Path in geometry to this module
    // The Default constructor for the AliITSgeomMatrix class. By Default
    // the angles of rotations are set to zero, meaning that the rotation
    // matrix is the unit matrix. The translation vector is also set to 
    // zero as are the module id number. The detector type is set to -1 
    // (an undefined value). The full rotation matrix is kept so that 
    // the evaluation  of a coordinate transformation can be done 
    // quickly and with a minimum of CPU overhead. The basic coordinate 
    // systems are the ALICE global coordinate system and the detector 
    // local coordinate system. In general this structure is not limited 
    // to just those two coordinate systems.
    //Begin_Html
    /*
      <img src="picts/ITS/AliITSgeomMatrix_L1.gif">
    */
    //End_Html
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A default constructes AliITSgeomMatrix class.
    Int_t i,j;

    fDetectorIndex = -1; // a value never defined.
    for(i=0;i<3;i++){
	fid[i] = 0;
	frot[i] = ftran[i] = 0.0;
	for(j=0;j<3;j++) fm[i][j] = 0.0;
	fCylR = fCylPhi = 0.0;
    }// end for i
    fm[0][0] = fm[1][1] = fm[2][2] = 1.0;
}
/*
//----------------------------------------------------------------------
AliITSgeomMatrix::AliITSgeomMatrix(const AliITSgeomMatrix &sourse) : 
    TObject(sourse){
    // The standard Copy constructor. This make a full / proper copy of
    // this class.
    // Inputs:
    //    AliITSgeomMatrix &source   The source of this copy
    // Outputs:
    //    none.
    // Return:
    //    A copy constructes AliITSgeomMatrix class.
	Int_t i,j;

	this->fDetectorIndex = sourse.fDetectorIndex;
	for(i=0;i<3;i++){
		this->fid[i]     = sourse.fid[i];
		this->frot[i]    = sourse.frot[i];
		this->ftran[i]   = sourse.ftran[i];
		this->fCylR      = sourse.fCylR;
		this->fCylPhi    = sourse.fCylPhi;
		for(j=0;j<3;j++) this->fm[i][j] = sourse.fm[i][j];
	}// end for i
     this->fPath   = sourse.fPath;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::operator=(const AliITSgeomMatrix &sourse){
    // The standard = operator. This make a full / proper copy of
    // this class.
    // The standard Copy constructor. This make a full / proper copy of
    // this class.
    // Inputs:
    //    AliITSgeomMatrix &source   The source of this copy
    // Outputs:
    //    none.
    // Return:
    //    A copy of the source AliITSgeomMatrix class.
	Int_t i,j;

	this->fDetectorIndex = sourse.fDetectorIndex;
	for(i=0;i<3;i++){
		this->fid[i]     = sourse.fid[i];
		this->frot[i]    = sourse.frot[i];
		this->ftran[i]   = sourse.ftran[i];
		this->fCylR      = sourse.fCylR;
		this->fCylPhi    = sourse.fCylPhi;
		for(j=0;j<3;j++) this->fm[i][j] = sourse.fm[i][j];
	}// end for i
     this->fPath   = sourse.fPath;
}
*/
//----------------------------------------------------------------------
AliITSgeomMatrix::AliITSgeomMatrix(Int_t idt,const Int_t id[3],
                        const Double_t rot[3],const Double_t tran[3]):
TObject(),
fDetectorIndex(idt), // Detector type index (like fShapeIndex was)
fid(),       // layer, ladder, detector numbers.
frot(),      //! vector of rotations about x,y,z [radians].
ftran(),     // Translation vector of module x,y,z.
fCylR(0.0),  //! R Translation in Cylinderical coordinates
fCylPhi(0.0),//! Phi Translation vector in Cylindrical coord.
fm(),        // Rotation matrix based on frot.
fPath(){     // Path in geometry to this moduel
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
      <img src="picts/ITS/AliITSgeomMatrix_L1.gif">
    */
    //End_Html
    // Inputs:
    //    Int_t idt        The detector index value
    //    Int_t id[3]      The layer, ladder, and detector numbers
    //    Double_t rot[3]  The 3 Cartician rotaion angles [radians]
    //    Double_t tran[3] The 3 Cartician translation distnaces
    // Outputs:
    //    none.
    // Return:
    //    A properly inilized AliITSgeomMatrix class.
    Int_t i;

    for(i=0;i<3;i++){
	fid[i]   = id[i];
	frot[i]  = rot[i];
	ftran[i] = tran[i];
    }// end for i
    fCylR   = TMath::Sqrt(ftran[0]*ftran[0]+ftran[1]*ftran[1]);
    fCylPhi = TMath::ATan2(ftran[1],ftran[0]);
    if(fCylPhi<0.0) fCylPhi += TMath::Pi();
    this->MatrixFromAngle();
}
//----------------------------------------------------------------------
AliITSgeomMatrix::AliITSgeomMatrix(Int_t idt, const Int_t id[3],
                                   Double_t matrix[3][3],
                                   const Double_t tran[3]):
TObject(),
fDetectorIndex(idt), // Detector type index (like fShapeIndex was)
fid(),       // layer, ladder, detector numbers.
frot(),      //! vector of rotations about x,y,z [radians].
ftran(),     // Translation vector of module x,y,z.
fCylR(0.0),  //! R Translation in Cylinderical coordinates
fCylPhi(0.0),//! Phi Translation vector in Cylindrical coord.
fm(),        // Rotation matrix based on frot.
fPath(){     // Path in geometry to this module
    // This is a constructor for the AliITSgeomMatrix class. The 
    // rotation matrix is given as one of the inputs, and the 
    // translation vector tran [cm]. In  addition the layer, ladder, 
    // and detector number for this particular module and the type of 
    // module must be given. The full rotation matrix is kept so that 
    // the evaluation of a coordinate transformation can be done quickly 
    // and with a minimum of CPU overhead. The basic coordinate systems 
    // are the ALICE global coordinate system and the detector local
    // coordinate system. In general this structure is not limited to just
    // those two coordinate systems.
    //Begin_Html
    /*
      <img src="picts/ITS/AliITSgeomMatrix_L1.gif">
    */
    //End_Html
    // Inputs:
    //    Int_t idt          The detector index value
    //    Int_t id[3]        The layer, ladder, and detector numbers
    //    Double_t rot[3][3] The 3x3 Cartician rotaion matrix
    //    Double_t tran[3]   The 3 Cartician translation distnaces
    // Outputs:
    //    none.
    // Return:
    //    A properly inilized AliITSgeomMatrix class.
    Int_t i,j;

    for(i=0;i<3;i++){
	fid[i]   = id[i];
	ftran[i] = tran[i];
	for(j=0;j<3;j++) fm[i][j] = matrix[i][j];
    }// end for i
    fCylR   = TMath::Sqrt(ftran[0]*ftran[0]+ftran[1]*ftran[1]);
    fCylPhi = TMath::ATan2(ftran[1],ftran[0]);
    if(fCylPhi<0.0) fCylPhi += TMath::Pi();
    this->AngleFromMatrix();
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::SixAnglesFromMatrix(Double_t *ang)const{
    // This function returns the 6 GEANT 3.21 rotation angles [degrees] in
    // the array ang which must be at least [6] long.
    // Inputs:
    //   none.
    // Outputs:
    //   Double_t ang[6]  The 6 Geant3.21 rotation angles. [degrees]
    // Return:
    //   noting
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
    // Given the 6 GEANT 3.21 rotation angles [degree], this will compute and
    // set the rotations matrix and 3 standard rotation angles [radians].
    // These angles and rotation matrix are overwrite the existing values in
    // this class.
    // Inputs:
    //   Double_t ang[6]  The 6 Geant3.21 rotation angles. [degrees]
    // Outputs:
    //   none.
    // Return:
    //   noting
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
                                   Int_t idt,const Int_t id[3],
                                   const Double_t tran[3]):
TObject(),
fDetectorIndex(idt),
fCylR(0.),
fCylPhi(0.),
fPath(){
    // This is a constructor for the AliITSgeomMatrix class. The matrix 
    // is defined by the 6 GEANT 3.21 rotation angles [degrees], and 
    // the translation vector tran [cm]. In addition the layer, ladder, 
    // and detector number for this particular module and the type of 
    // module must be given. The full rotation matrix is kept so that 
    // the evaluation  of a coordinate transformation can be done 
    // quickly and with a minimum of CPU overhead. The basic coordinate 
    // systems are the ALICE global coordinate system and the detector 
    // local coordinate system. In general this structure is not limited 
    // to just those two coordinate systems.
    //Begin_Html
    /*
      <img src="picts/ITS/AliITSgeomMatrix_L1.gif">
    */
    //End_Html
    // Inputs:
    //    Double_t rotd[6]  The 6 Geant 3.21 rotation angles [degrees]
    //    Int_t idt         The module Id number
    //    Int_t id[3]       The layer, ladder and detector number
    //    Double_t tran[3]  The translation vector
    Int_t i;

    for(i=0;i<3;i++){
	fid[i]   = id[i];
	ftran[i] = tran[i];
    }// end for i
    fCylR   = TMath::Sqrt(ftran[0]*ftran[0]+ftran[1]*ftran[1]);
    fCylPhi = TMath::ATan2(ftran[1],ftran[0]);
    if(fCylPhi<0.0) fCylPhi += TMath::Pi();
    this->MatrixFromSixAngles(rotd);
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::AngleFromMatrix(){
    // Computes the angles from the rotation matrix up to a phase of 
    // 180 degrees. The matrix used in AliITSgeomMatrix::MatrixFromAngle()
    // and  its inverse AliITSgeomMatrix::AngleFromMatrix() are defined in 
    // the following ways, R = Rz*Ry*Rx (M=R*L+T) where
    //     1   0   0       Cy  0 -Sy       Cz -Sz  0
    // Rx= 0   Cx -Sx  Ry=  0  1   0   Rz= Sz  Cz  0
    //     0   Sx  Cx      Sy  0  Cy        0   0  1
    // The choice of the since of S, comes from the choice between 
    // the rotation of the object or the coordinate system (view). I think
    // that this choice is the first, the rotation of the object.
    // Inputs:
    //   none
    // Outputs:
    //   none
    // Return:
    //   none
    Double_t rx,ry,rz;
    // get angles from matrix up to a phase of 180 degrees.

    rx = TMath::ATan2(fm[2][1],fm[2][2]);if(rx<0.0) rx += 2.0*TMath::Pi();
    ry = TMath::ASin(fm[0][2]);          if(ry<0.0) ry += 2.0*TMath::Pi();
    rz = TMath::ATan2(fm[1][0],fm[0][0]);if(rz<0.0) rz += 2.0*TMath::Pi();
    frot[0] = rx;
    frot[1] = ry;
    frot[2] = rz;
    return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::MatrixFromAngle(){
    // Computes the Rotation matrix from the angles [radians] kept in this
    // class. The matrix used in AliITSgeomMatrix::MatrixFromAngle() and 
    // its inverse AliITSgeomMatrix::AngleFromMatrix() are defined in 
    // the following ways, R = Rz*Ry*Rx (M=R*L+T) where
    //     1   0   0       Cy  0 -Sy       Cz -Sz  0
    // Rx= 0   Cx -Sx  Ry=  0  1   0   Rz= Sz  Cz  0
    //     0   Sx  Cx      Sy  0  Cy        0   0  1
    // The choice of the since of S, comes from the choice between 
    // the rotation of the object or the coordinate system (view). I think
    // that this choice is the first, the rotation of the object.
    // Inputs:
    //   none
    // Outputs:
    //   none
    // Return:
    //   none
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
void AliITSgeomMatrix::GtoLPosition(const Double_t g0[3],Double_t l[3]) const {
    // Returns the local coordinates given the global coordinates [cm].
    // Inputs:
    //   Double_t g[3]   The position represented in the ALICE 
    //                   global coordinate system
    // Outputs:
    //   Double_t l[3]  The poistion represented in the local
    //                  detector coordiante system
    // Return:
    //   none
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
void AliITSgeomMatrix::LtoGPosition(const Double_t l[3],Double_t g[3]) const {
    // Returns the global coordinates given the local coordinates [cm].
    // Inputs:
    //   Double_t l[3]   The poistion represented in the detector 
    //                   local coordinate system
    // Outputs:
    //   Double_t g[3]   The poistion represented in the ALICE
    //                   Global coordinate system
    // Return:
    //   none.
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
void AliITSgeomMatrix::GtoLMomentum(const Double_t g[3],Double_t l[3]) const{
    // Returns the local coordinates of the momentum given the global
    // coordinates of the momentum. It transforms just like GtoLPosition
    // except that the translation vector is zero.
    // Inputs:
    //   Double_t g[3] The momentum represented in the ALICE global 
    //                 coordinate system
    // Outputs:
    //   Double_t l[3] the momentum represented in the detector 
    //                 local coordinate system
    // Return:
    //   none.
	Int_t    i,j;

	for(i=0;i<3;i++){
		l[i] = 0.0;
		for(j=0;j<3;j++) l[i] += fm[i][j]*g[j];
		// g = R l
	} // end for i
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGMomentum(const Double_t l[3],Double_t g[3]) const {
    // Returns the Global coordinates of the momentum given the local
    // coordinates of the momentum. It transforms just like LtoGPosition
    // except that the translation vector is zero.
    // Inputs:
    //   Double_t l[3] the momentum represented in the detector 
    //                 local coordinate system
    // Outputs:
    //   Double_t g[3] The momentum represented in the ALICE global 
    //                 coordinate system
    // Return:
    //   none.
	Int_t    i,j;

	for(i=0;i<3;i++){
		g[i] = 0.0;
		for(j=0;j<3;j++) g[i] += fm[j][i]*l[j];
		// g = R^t l
	} // end for i
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::GtoLPositionError(const Double_t g[3][3],
                                         Double_t l[3][3]) const {
    // Given an Uncertainty matrix in Global coordinates it is 
    // rotated so that  its representation in local coordinates can 
    // be returned. There is no effect due to the translation vector 
    // or its uncertainty.
    // Inputs:
    //   Double_t g[3][3] The error matrix represented in the ALICE global 
    //                    coordinate system
    // Outputs:
    //   Double_t l[3][3] the error matrix represented in the detector 
    //                    local coordinate system
    // Return:
    //   none.
	Int_t    i,j,k,m;

	for(i=0;i<3;i++)for(m=0;m<3;m++){
	    l[i][m] = 0.0;
	    for(j=0;j<3;j++)for(k=0;k<3;k++)
		l[i][m] += fm[j][i]*g[j][k]*fm[k][m];
	} // end for i,m
	    // g = R^t l R
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGPositionError(const Double_t l[3][3],
                                               Double_t g[3][3]) const {
    // Given an Uncertainty matrix in Local coordinates it is rotated so that 
    // its representation in global coordinates can be returned. There is no
    // effect due to the translation vector or its uncertainty.
    // Inputs:
    //   Double_t l[3][3] the error matrix represented in the detector 
    //                    local coordinate system
    // Outputs:
    //   Double_t g[3][3] The error matrix represented in the ALICE global 
    //                    coordinate system
    // Return:
    //   none.
	Int_t    i,j,k,m;

	for(i=0;i<3;i++)for(m=0;m<3;m++){
	    g[i][m] = 0.0;
	    for(j=0;j<3;j++)for(k=0;k<3;k++)
		g[i][m] += fm[i][j]*l[j][k]*fm[m][k];
	} // end for i,m
	    // g = R l R^t
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::GtoLPositionTracking(const Double_t g[3],
                                            Double_t l[3]) const {
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
    // Inputs:
    //   Double_t g[3]   The position represented in the ALICE 
    //                   global coordinate system
    // Outputs:
    //   Double_t l[3]  The poistion represented in the local
    //                  detector coordiante system
    // Return:
    //   none
    Double_t l0[3];

    this->GtoLPosition(g,l0);
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
                                            Double_t g[3]) const {
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
    // Inputs:
    //   Double_t l[3]   The poistion represented in the detector 
    //                   local coordinate system
    // Outputs:
    //   Double_t g[3]   The poistion represented in the ALICE
    //                   Global coordinate system
    // Return:
    //   none.
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
                                            Double_t l[3]) const {
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
    // Inputs:
    //   Double_t g[3] The momentum represented in the ALICE global 
    //                 coordinate system
    // Outputs:
    //   Double_t l[3] the momentum represented in the detector 
    //                 local coordinate system
    // Return:
    //   none.
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
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGMomentumTracking(const Double_t l[3],
                                            Double_t g[3]) const {
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
    // Inputs:
    //   Double_t l[3] the momentum represented in the detector 
    //                 local coordinate system
    // Outputs:
    //   Double_t g[3] The momentum represented in the ALICE global 
    //                 coordinate system
    // Return:
    //   none.
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
void AliITSgeomMatrix::GtoLPositionErrorTracking(const Double_t g[3][3],
                                                 Double_t l[3][3]) const {
    // A slightly different coordinate system is used when tracking.
    // This coordinate system is only relevant when the geometry represents
    // the cylindrical ALICE ITS geometry. For tracking the Z axis is left
    // alone but X -> -Y and Y -> X such that X always points out of the
    // ITS Cylinder for every layer including layer 1 (where the detector 
    // are mounted upside down).
    //Begin_Html
    /*
      <img src="picts/ITS/AliITSgeomMatrix_TE1.gif">
    */
    //End_Html
    // Inputs:
    //   Double_t g[3][3] The error matrix represented in the ALICE global 
    //                    coordinate system
    // Outputs:
    //   Double_t l[3][3] the error matrix represented in the detector 
    //                    local coordinate system
    // Return:
	Int_t    i,j,k,m;
	Double_t rt[3][3];
	Double_t a0[3][3] = {{0.,+1.,0.},{-1.,0.,0.},{0.,0.,+1.}};
	Double_t a1[3][3] = {{0.,-1.,0.},{+1.,0.,0.},{0.,0.,+1.}};

	if(fid[0]==1) for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)
	    rt[i][k] = a0[i][j]*fm[j][k];
	else for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)
	    rt[i][k] = a1[i][j]*fm[j][k];
	for(i=0;i<3;i++)for(m=0;m<3;m++){
	    l[i][m] = 0.0;
	    for(j=0;j<3;j++)for(k=0;k<3;k++)
		l[i][m] += rt[j][i]*g[j][k]*rt[k][m];
	} // end for i,m
	    // g = R^t l R
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::LtoGPositionErrorTracking(const Double_t l[3][3],
                                                 Double_t g[3][3]) const {
    // A slightly different coordinate system is used when tracking.
    // This coordinate system is only relevant when the geometry represents
    // the cylindrical ALICE ITS geometry. For tracking the Z axis is left
    // alone but X -> -Y and Y -> X such that X always points out of the
    // ITS Cylinder for every layer including layer 1 (where the detector 
    // are mounted upside down).
    //Begin_Html
    /*
      <img src="picts/ITS/AliITSgeomMatrix_TE1.gif">
    */
    //End_Html
    // Inputs:
    //   Double_t l[3][3] the error matrix represented in the detector 
    //                    local coordinate system
    // Outputs:
    //   Double_t g[3][3] The error matrix represented in the ALICE global 
    //                    coordinate system
    // Return:
    //   none.
	Int_t    i,j,k,m;
	Double_t rt[3][3];
	Double_t a0[3][3] = {{0.,+1.,0.},{-1.,0.,0.},{0.,0.,+1.}};
	Double_t a1[3][3] = {{0.,-1.,0.},{+1.,0.,0.},{0.,0.,+1.}};

	if(fid[0]==1) for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)
	    rt[i][k] = a0[i][j]*fm[j][k];
	else for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)
	    rt[i][k] = a1[i][j]*fm[j][k];
	for(i=0;i<3;i++)for(m=0;m<3;m++){
	    g[i][m] = 0.0;
	    for(j=0;j<3;j++)for(k=0;k<3;k++)
		g[i][m] += rt[i][j]*l[j][k]*rt[m][k];
	} // end for i,m
	    // g = R l R^t
	return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::PrintTitles(ostream *os) const {
    // Standard output format for this class but it includes variable
    // names and formatting that makes it easer to read.
    // Inputs:
    //    ostream *os   The output stream to print the title on
    // Outputs:
    //    none.
    // Return:
    //    none.
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
void AliITSgeomMatrix::PrintComment(ostream *os) const {
    //  output format used by Print.
    // Inputs:
    //    ostream *os   The output stream to print the comments on
    // Outputs:
    //    none.
    // Return:
    //    none.
    *os << "fDetectorIndex fid[0] fid[1] fid[2] ftran[0] ftran[1] ftran[2] ";
    *os << "fm[0][0]  fm[0][1]  fm[0][2]  fm[1][0]  fm[1][1]  fm[1][2]  ";
    *os << "fm[2][0]  fm[2][1]  fm[2][2] ";
    return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::Print(ostream *os)const{
    // Standard output format for this class.
    // Inputs:
    //    ostream *os   The output stream to print the class data on
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,j;
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC || defined __xlC__
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << fDetectorIndex << " ";
    for(i=0;i<3;i++) *os << fid[i]   << " ";
//    for(i=0;i<3;i++) *os << frot[i]  << " ";  // Redundant with fm[][].
    for(i=0;i<3;i++) *os << setprecision(16) << ftran[i] << " ";
    for(i=0;i<3;i++)for(j=0;j<3;j++)  *os << setprecision(16) << 
					  fm[i][j] << " ";
    *os << fPath.Length()<< " ";
    for(i=0;i<fPath.Length();i++) *os << fPath[i];
    *os << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::Read(istream *is){
    // Standard input format for this class.
    // Inputs:
    //    istream *is   The input stream to read on
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i,j;

    *is >> fDetectorIndex;
    for(i=0;i<3;i++) *is >> fid[i];
//    for(i=0;i<3;i++) *is >> frot[i]; // Redundant with fm[][].
    for(i=0;i<3;i++) *is >> ftran[i];
    for(i=0;i<3;i++)for(j=0;j<3;j++)  *is >> fm[i][j];
    while(is->peek()==' ')is->get(); // skip white spaces
    if(isprint(is->peek())){ // old format did not have path.
	*is >> j; // string length
	fPath.Resize(j);
	for(i=0;i<j;i++) {*is >> fPath[i];}
    } // end if
    AngleFromMatrix(); // compute angles frot[].
    fCylR   = TMath::Sqrt(ftran[0]*ftran[0]+ftran[1]*ftran[1]);
    fCylPhi = TMath::ATan2(ftran[1],ftran[0]);
    if(fCylPhi<0.0) fCylPhi += TMath::Pi();
    return;
}
//______________________________________________________________________
void AliITSgeomMatrix::Streamer(TBuffer &R__b){
   // Stream an object of class AliITSgeomMatrix.
    // Inputs:
    //     TBuffer &R__b   The output buffer to stream data on.
    // Outputs:
    //    none.
    // Return:
    //    none.

    if (R__b.IsReading()) {
        AliITSgeomMatrix::Class()->ReadBuffer(R__b, this);
        fCylR   = TMath::Sqrt(ftran[0]*ftran[0]+ftran[1]*ftran[1]);
        fCylPhi = TMath::ATan2(ftran[1],ftran[0]);
        this->AngleFromMatrix();
        if(fCylPhi<0.0) fCylPhi += TMath::Pi();
    } else {
        AliITSgeomMatrix::Class()->WriteBuffer(R__b, this);
    } // end if
}
//______________________________________________________________________
void AliITSgeomMatrix::SetTranslation(const Double_t tran[3]){
    // Sets the translation vector and computes fCylR and fCylPhi.
    // Inputs:
    //   Double_t trans[3]   The translation vector to be used
    // Outputs:
    //   none.
    // Return:
    //   none.
    for(Int_t i=0;i<3;i++) ftran[i] = tran[i];
    fCylR   = TMath::Sqrt(ftran[0]*ftran[0]+ftran[1]*ftran[1]);
    fCylPhi = TMath::ATan2(ftran[1],ftran[0]);
    if(fCylPhi<0.0) fCylPhi += TMath::Pi();
}
//----------------------------------------------------------------------
TPolyLine3D* AliITSgeomMatrix::CreateLocalAxis() const {
    // This class is used as part of the documentation of this class
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   A pointer to a new TPolyLine3D object showing the 3 line
    //   segments that make up the this local axis in the global
    //   reference system.
    Float_t  gf[15];
    Double_t g[5][3];
    Double_t l[5][3]={{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,0.0},
                      {0.0,0.0,1.0}};
    Int_t i;

    for(i=0;i<5;i++) {
        LtoGPosition(l[i],g[i]);
        gf[3*i]=(Float_t)g[i][0];
        gf[3*i+1]=(Float_t)g[i][1];
        gf[3*i+2]=(Float_t)g[i][2];
    } // end for i
    return new TPolyLine3D(5,gf);
}
//----------------------------------------------------------------------
TPolyLine3D* AliITSgeomMatrix::CreateLocalAxisTracking() const {
    // This class is used as part of the documentation of this class
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   A pointer to a new TPolyLine3D object showing the 3 line
    //   segments that make up the this local axis in the global
    //   reference system.
    Float_t gf[15];
    Double_t g[5][3];
    Double_t l[5][3]={{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,0.0},
                      {0.0,0.0,1.0}};
    Int_t i;

    for(i=0;i<5;i++) {
        LtoGPositionTracking(l[i],g[i]);
        gf[3*i]=(Float_t)g[i][0];
        gf[3*i+1]=(Float_t)g[i][1];
        gf[3*i+2]=(Float_t)g[i][2];
    } // end for i
    return new TPolyLine3D(5,gf);
}
//----------------------------------------------------------------------
TNode* AliITSgeomMatrix::CreateNode(const Char_t *nodeName,
                                    const Char_t *nodeTitle,TNode *mother,
                                    TShape *shape,Bool_t axis) const {
    // Creates a node inside of the node mother out of the shape shape
    // in the position, with respect to mother, indecated by "this". If axis
    // is ture, it will insert an axis within this node/shape.
    // Inputs:
    //   Char_t *nodeName  This name of this node
    //   Char_t *nodeTitle This node title
    //   TNode  *mother    The node this node will be inside of/with respect to
    //   TShape *shape     The shape of this node
    //   Bool_t axis       If ture, a set of x,y,z axis will be included
    // Outputs:
    //   none.
    // Return:
    //   A pointer to "this" node.
    Double_t trans[3],matrix[3][3],*matr;
    TRotMatrix *rot = new TRotMatrix();
    TString name,title;

    matr = &(matrix[0][0]);
    this->GetTranslation(trans);
    this->GetMatrix(matrix);
    rot->SetMatrix(matr);
    //
    name = nodeName;
    title = nodeTitle;
    //
    mother->cd();
    TNode *node1 = new TNode(name.Data(),title.Data(),shape,trans[0],trans[1],trans[2],rot);
    if(axis){
        Int_t i,j;
        const Float_t kScale=0.5,kLw=0.2;
        Float_t xchar[13][2]={{0.5*kLw,1.},{0.,0.5*kLw},{0.5-0.5*kLw,0.5},
                              {0.,0.5*kLw},{0.5*kLw,0.},{0.5,0.5-0.5*kLw},
                              {1-0.5*kLw,0.},{1.,0.5*kLw},{0.5+0.5*kLw,0.5},
                              {1.,1.-0.5*kLw},{1.-0.5*kLw,1.},{0.5,0.5+0.5*kLw},
                              {0.5*kLw,1.}};
        Float_t ychar[10][2]={{.5-0.5*kLw,0.},{.5+0.5*kLw,0.},{.5+0.5*kLw,0.5-0.5*kLw},
                              {1.,1.-0.5*kLw},{1.-0.5*kLw,1.},{0.5+0.5*kLw,0.5},
                              {0.5*kLw,1.}   ,{0.,1-0.5*kLw} ,{0.5-0.5*kLw,0.5},
                              {.5-0.5*kLw,0.}};
        Float_t zchar[11][2]={{0.,1.},{0,1.-kLw},{1.-kLw,1.-kLw},{0.,kLw}   ,{0.,0.},
                              {1.,0.},{1.,kLw}  ,{kLw,kLw}      ,{1.,1.-kLw},{1.,1.},
                              {0.,1.}};
        for(i=0;i<13;i++)for(j=0;j<2;j++){
            if(i<13) xchar[i][j] = kScale*xchar[i][j];
            if(i<10) ychar[i][j] = kScale*ychar[i][j];
            if(i<11) zchar[i][j] = kScale*zchar[i][j];
        } // end for i,j
        TXTRU *axisxl = new TXTRU("x","x","text",12,2);
        for(i=0;i<12;i++) axisxl->DefineVertex(i,xchar[i][0],xchar[i][1]);
        axisxl->DefineSection(0,-0.5*kLw);axisxl->DefineSection(1,0.5*kLw);
        TXTRU *axisyl = new TXTRU("y","y","text",9,2);
        for(i=0;i<9;i++) axisyl->DefineVertex(i,ychar[i][0],ychar[i][1]);
        axisyl->DefineSection(0,-0.5*kLw);axisyl->DefineSection(1,0.5*kLw);
        TXTRU *axiszl = new TXTRU("z","z","text",10,2);
        for(i=0;i<10;i++) axiszl->DefineVertex(i,zchar[i][0],zchar[i][1]);
        axiszl->DefineSection(0,-0.5*kLw);axiszl->DefineSection(1,0.5*kLw);
        Float_t lxy[13][2]={{-0.5*kLw,-0.5*kLw},{0.8,-0.5*kLw},{0.8,-0.1},{1.0,0.0},
                            {0.8,0.1},{0.8,0.5*kLw},{0.5*kLw,0.5*kLw},{0.5*kLw,0.8},
                            {0.1,0.8},{0.0,1.0},{-0.1,0.8},{-0.5*kLw,0.8},
                            {-0.5*kLw,-0.5*kLw}};
        TXTRU *axisxy = new TXTRU("axisxy","axisxy","text",13,2);
        for(i=0;i<13;i++) axisxy->DefineVertex(i,lxy[i][0],lxy[i][1]);
        axisxy->DefineSection(0,-0.5*kLw);axisxy->DefineSection(1,0.5*kLw);
        Float_t lz[8][2]={{0.5*kLw,-0.5*kLw},{0.8,-0.5*kLw},{0.8,-0.1},{1.0,0.0},
                           {0.8,0.1},{0.8,0.5*kLw},{0.5*kLw,0.5*kLw},
                           {0.5*kLw,-0.5*kLw}};
        TXTRU *axisz = new TXTRU("axisz","axisz","text",8,2);
        for(i=0;i<8;i++) axisz->DefineVertex(i,lz[i][0],lz[i][1]);
        axisz->DefineSection(0,-0.5*kLw);axisz->DefineSection(1,0.5*kLw);
        //TRotMatrix *xaxis90= new TRotMatrix("xaixis90","",90.0, 0.0, 0.0);
        TRotMatrix *yaxis90= new TRotMatrix("yaixis90","", 0.0,90.0, 0.0);
        TRotMatrix *zaxis90= new TRotMatrix("zaixis90","", 0.0, 0.0,90.0);
        //
        node1->cd();
        title = name.Append("axisxy");
        TNode *nodeaxy = new TNode(title.Data(),title.Data(),axisxy);
        title = name.Append("axisz");
        TNode *nodeaz = new TNode(title.Data(),title.Data(),axisz,0.,0.,0.,yaxis90);
        TNode *textboxX0 = new TNode("textboxX0","textboxX0",axisxl,
                                    lxy[3][0],lxy[3][1],0.0);
        TNode *textboxX1 = new TNode("textboxX1","textboxX1",axisxl,
                                    lxy[3][0],lxy[3][1],0.0,yaxis90);
        TNode *textboxX2 = new TNode("textboxX2","textboxX2",axisxl,
                                    lxy[3][0],lxy[3][1],0.0,zaxis90);
        TNode *textboxY0 = new TNode("textboxY0","textboxY0",axisyl,
                                    lxy[9][0],lxy[9][1],0.0);
        TNode *textboxY1 = new TNode("textboxY1","textboxY1",axisyl,
                                    lxy[9][0],lxy[9][1],0.0,yaxis90);
        TNode *textboxY2 = new TNode("textboxY2","textboxY2",axisyl,
                                    lxy[9][0],lxy[9][1],0.0,zaxis90);
        TNode *textboxZ0 = new TNode("textboxZ0","textboxZ0",axiszl,
                                    0.0,0.0,lz[3][0]);
        TNode *textboxZ1 = new TNode("textboxZ1","textboxZ1",axiszl,
                                    0.0,0.0,lz[3][0],yaxis90);
        TNode *textboxZ2 = new TNode("textboxZ2","textboxZ2",axiszl,
                                    0.0,0.0,lz[3][0],zaxis90);
        nodeaxy->Draw();
        nodeaz->Draw();
        textboxX0->Draw();
        textboxX1->Draw();
        textboxX2->Draw();
        textboxY0->Draw();
        textboxY1->Draw();
        textboxY2->Draw();
        textboxZ0->Draw();
        textboxZ1->Draw();
        textboxZ2->Draw();
    } // end if
    mother->cd();
    return node1;
}
//----------------------------------------------------------------------
void AliITSgeomMatrix::MakeFigures() const {
    // make figures to help document this class
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    const Double_t kDx0=550.,kDy0=550.,kDz0=550.; // cm
    const Double_t kDx=1.0,kDy=0.300,kDz=3.0,kRmax=0.1; // cm
    Float_t l[5][3]={{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,0.0},
                      {0.0,0.0,1.0}};
    TCanvas *c = new TCanvas(kFALSE);// create a batch mode canvas.
    TView   *view = new TView(1); // Create Cartesian coordiante view
    TBRIK   *mother  = new TBRIK("Mother","Mother","void",kDx0,kDy0,kDz0);
    TBRIK   *det  = new TBRIK("Detector","","Si",kDx,kDy,kDz);
    TPolyLine3D *axis = new TPolyLine3D(5,&(l[0][0]));
    TPCON *arrow      = new TPCON("arrow","","air",0.0,360.,2);
    TRotMatrix *xarrow= new TRotMatrix("xarrow","",90.,0.0,0.0);
    TRotMatrix *yarrow= new TRotMatrix("yarrow","",0.0,90.,0.0);

    det->SetLineColor(0); // black
    det->SetLineStyle(1); // solid line
    det->SetLineWidth(2); // pixel units
    det->SetFillColor(1); // black
    det->SetFillStyle(4010); // window is 90% transparent
    arrow->SetLineColor(det->GetLineColor());
    arrow->SetLineWidth(det->GetLineWidth());
    arrow->SetLineStyle(det->GetLineStyle());
    arrow->SetFillColor(1); // black
    arrow->SetFillStyle(4100); // window is 100% opaque
    arrow->DefineSection(0,0.0,0.0,kRmax);
    arrow->DefineSection(1,2.*kRmax,0.0,0.0);
    view->SetRange(-kDx0,-kDy0,-kDz0,kDx0,kDy0,kDz0);
    //
    TNode *node0 = new TNode("NODE0","NODE0",mother);
    node0->cd();
    TNode *node1 = new TNode("NODE1","NODE1",det);
    node1->cd();
    TNode *nodex = new TNode("NODEx","NODEx",arrow,l[0][0],l[0][1],l[0][2],xarrow);
    TNode *nodey = new TNode("NODEy","NODEy",arrow,l[2][0],l[2][1],l[2][2],yarrow);
    TNode *nodez = new TNode("NODEz","NODEz",arrow,l[4][0],l[4][1],l[4][2]);
    //
    axis->Draw();
    nodex->Draw();
    nodey->Draw();
    nodez->Draw();
    
    //
    node0->cd();
    node0->Draw();
    c->Update();
    c->SaveAs("AliITSgeomMatrix_L1.gif");
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSgeomMatrix &p){
    // Standard output streaming function.
    // Inputs:
    //    ostream &os          The output stream to print the class data on
    //    AliITSgeomMatrix &p  This class
    // Outputs:
    //    none.
    // Return:
    //    none.

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomMatrix &r){
    // Standard input streaming function.
    // Inputs:
    //    ostream &os          The input stream to print the class data on
    //    AliITSgeomMatrix &p  This class
    // Outputs:
    //    none.
    // Return:
    //    none.

    r.Read(&is);
    return is;
}
//----------------------------------------------------------------------
