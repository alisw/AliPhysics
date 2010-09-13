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


////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliTPCCalibGlobalMisalignment class                                        //
// The class calculates the space point distortions due to simple         // 
// misalignments like shifts in caresian coordinates or a rotation        //
// of the TPC read out planes (A and C side)                              //
// Optionaly possible to use it for visualization of the alignemnt form the Alignment OCDB //
// fUseGeoManager has to be set to kTRUE to enable this option
//                                                                        //
// date: 06/05/2010                                                       //
// Authors: Stefan Rossegger, Jim Thomas, Magnus Mager                    //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalibGlobalMisalignment.h"
#include "TMath.h"
#include "TGeoMatrix.h"
#include "AliTPCROC.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include <TGeoPhysicalNode.h>


AliTPCCalibGlobalMisalignment::AliTPCCalibGlobalMisalignment()
  : AliTPCCorrection("mialign","Misalignment"),
    fXShift(0.),fYShift(0.),fZShift(0.),
    fRotPhiA(0.),fRotPhiC(0.),
    fdRPhiOffsetA(0.), fdRPhiOffsetC(0.), 
    fQuadrantDQ1(0), fQuadrantDQ2(0), fQuadrantQ2(0), fMatrixGlobal(0),
    fMatrixASide(0), fMatrixCSide(0),
    fUseGeomanager(kFALSE)
{
  //
  // default constructor
  //
}

AliTPCCalibGlobalMisalignment::~AliTPCCalibGlobalMisalignment() {
  //
  // default destructor
  //
  delete fQuadrantDQ1;   //OROC medium pads delta ly+ - ly-
  delete fQuadrantDQ2;   //OROC long   pads delta ly+ - ly-
  delete fQuadrantQ2;    //OROC long   pads - OROC medium pads

}
 
void AliTPCCalibGlobalMisalignment::SetQuadranAlign(const TVectorD *dq1, const TVectorD *dq2, const TVectorD *q2){
  //
  // Set quadrant alignment
  // 3 vectors for 36 (super) sectors
  //
  fQuadrantDQ1 = new TVectorD(*dq1);    //OROC medium pads delta ly+ - ly-
  fQuadrantDQ2 = new TVectorD(*dq2);;   //OROC long   pads delta ly+ - ly-
  fQuadrantQ2  = new TVectorD(*q2);;    //OROC long   pads - OROC medium pads
}

void AliTPCCalibGlobalMisalignment::SetGlobalAlign(const TGeoMatrix * matrixGlobal, const TGeoMatrix *matrixA, const TGeoMatrix *matrixC ){
  //
  // Set global misalignment as TGeoMatrix
  // 
  if (matrixGlobal) fMatrixGlobal = new TGeoHMatrix(*matrixGlobal);
  if (matrixA) fMatrixASide = new TGeoHMatrix(*matrixA);
  if (matrixC) fMatrixCSide = new TGeoHMatrix(*matrixC);
}


//void AliTPCCalibGlobalMisalignment::Init() {
//  //
// // Initialization funtion
//  //

//  // nothing to be initialized, results of this calibration class will go to the global aligment structure

//}

//void AliTPCCalibGlobalMisalignment::Update(const TTimeStamp &/*timeStamp*/) {
//  //
//  // Update function 
//  //
//
//  // nothing to be updated, results of this calibration class will go to the global aligment structure
//
//}



void AliTPCCalibGlobalMisalignment::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the simple correction due to a shift (in x,y,z) or an rotation of the TPC (around z)
  //  
  static AliTPCROC *tpcRoc =AliTPCROC::Instance();  
  Double_t r=0, phi=0;
  r   = TMath::Sqrt( x[0]*x[0] + x[1]*x[1] );
  phi = TMath::ATan2(x[1],x[0]);
  // Getsector number
  Double_t sec=TMath::Nint(-0.5+(phi*9./TMath::Pi()));
  if (sec<0) sec+=18;
  Int_t isec = TMath::Nint(sec);    
  if  (roc%36>=18) isec+=18;
  //
  // Get the point on the local coordiante frame
  //
  Double_t alpha=(sec+0.5)*TMath::Pi()/9;
  Double_t pos[3]={0,0,x[2]};
  pos[0]=  TMath::Cos(alpha)*x[0]+TMath::Sin(alpha)*x[1];
  pos[1]= -TMath::Sin(alpha)*x[0]+TMath::Cos(alpha)*x[1];
  if (pos[0]>tpcRoc->GetPadRowRadiiUp(0))  isec+=36;

  //
  // apply quadrant alignment if available - in local coordinate frame
  //
  Double_t posQG[3]={x[0],x[1],x[2]};
  if (fQuadrantDQ1){
    Double_t dly=0;
    Bool_t isQ1 = TMath::Abs(pos[0]-161)<28;
    Bool_t isQ2 = (pos[0]>189);
    if (isQ1){
      if (pos[1]>0.) dly+=(*fQuadrantDQ1)[isec];
      if (pos[1]<0.) dly-=(*fQuadrantDQ1)[isec];      
    }
    if (isQ2){
      dly+=(*fQuadrantQ2)[isec];
      if (pos[1]>0.) dly+=(*fQuadrantDQ2)[isec];
      if (pos[1]<0.) dly-=(*fQuadrantDQ2)[isec];
    }
    // Tranform the corrected point to the global frame
    posQG[0]=  TMath::Cos(alpha)*pos[0]-TMath::Sin(alpha)*(pos[1]+dly);
    posQG[1]=  TMath::Sin(alpha)*pos[0]+TMath::Cos(alpha)*(pos[1]+dly);
  }
  //
  // rotation of the read-out planes
  if  (roc%36<18) // A side
    phi += fRotPhiA;
  else         // C side
    phi += fRotPhiC;
  
  // Simply adding a constant dRPHi residual. PURELY FOR CALIBRATION PURPOSES
  if  (roc%36<18) // A side
    phi += fdRPhiOffsetA/r;
  else         // C side
    phi += fdRPhiOffsetC/r;
  
  dx[0] = r * TMath::Cos(phi) - x[0];
  dx[1] = r * TMath::Sin(phi) - x[1]; 
  dx[2] = 0.; 

  // Simple shifts
  dx[0] -= fXShift;
  dx[1] -= fYShift;
  dx[2] -= fZShift;
  // quadrant shifts
  dx[0] += (posQG[0]-x[0]);
  dx[1] += (posQG[1]-x[1]);
  //
  // alignment matrix in local frame
  //
  if (fUseGeomanager){ //loading from the OCDB
    Double_t posC[3] ={pos[0],pos[1],pos[2]};
    //
    //2. correct the point in the local frame
    AliTPCParam *param = AliTPCcalibDB::Instance()->GetParameters();
    if (!param){
      //AliFatal("OCDB not initialized");
    }
    TGeoHMatrix  *mat = param->GetClusterMatrix(isec);
    //
    if (mat) mat->LocalToMaster(pos,posC);
    Double_t posCG[3]={posC[0],posC[1],posC[2]};
    //3. tranform the corrected point to the global frame
    posCG[0]=  TMath::Cos(alpha)*posC[0]-TMath::Sin(alpha)*posC[1];
    posCG[1]=  TMath::Sin(alpha)*posC[0]+TMath::Cos(alpha)*posC[1];
    posCG[2]=  posC[2];
    //4. Add delta
    dx[0]+=posCG[0]-x[0];
    dx[1]+=posCG[1]-x[1];
    dx[2]+=posCG[2]-x[2];
  }
  if (fMatrixGlobal){
    // apply global alignment matrix
    Double_t ppos[3]={x[0],x[1],x[2]};
    Double_t pposC[3]={x[0],x[1],x[2]};
    fMatrixGlobal->LocalToMaster(ppos,pposC);
    dx[0]+=pposC[0]-ppos[0];
    dx[1]+=pposC[1]-ppos[1];
    dx[2]+=pposC[2]-ppos[2];
  }

  if (fMatrixASide && roc%36<18){
    // apply global alignment matrix
    Double_t ppos[3]={x[0],x[1],x[2]};
    Double_t pposC[3]={x[0],x[1],x[2]};
    fMatrixASide->LocalToMaster(ppos,pposC);
    dx[0]+=pposC[0]-ppos[0];
    dx[1]+=pposC[1]-ppos[1];
    dx[2]+=pposC[2]-ppos[2];
  }
  if (fMatrixCSide && roc%36>=18){
    // apply global alignment matrix
    Double_t ppos[3]={x[0],x[1],x[2]};
    Double_t pposC[3]={x[0],x[1],x[2]};
    fMatrixCSide->LocalToMaster(ppos,pposC);
    dx[0]+=pposC[0]-ppos[0];
    dx[1]+=pposC[1]-ppos[1];
    dx[2]+=pposC[2]-ppos[2];
  }
}

void AliTPCCalibGlobalMisalignment::Print(Option_t* /*option*/ ) const {
  //
  // Print function to check the settings 
  //
  printf("%s",GetTitle());  
  printf(" - Trivial Misalignments for calibration purposes: \n");
  printf(" - X-Shift: %1.3f cm, Y-Shift: %1.3f cm, Z-Shift: %1.3f cm \n",fXShift,fYShift,fZShift);
  printf(" - Phi-Rotations: A side: %1.5f rad, C side: %1.5f rad\n",fRotPhiA,fRotPhiC);
  printf(" - dRPhi offsets: A side: %1.5f cm, C side: %1.5f cm\n",fdRPhiOffsetA,fdRPhiOffsetC);
 
 
}
