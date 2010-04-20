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

//-----------------------------------------------------------------------------
/// \class AliMUONAlignment
/// Alignment class for the ALICE DiMuon spectrometer 
///
/// MUON specific alignment class which interface to AliMillepede.   
/// For each track ProcessTrack calculates the local and global derivatives
/// at each cluster and fill the corresponding local equations. Provide methods
/// for fixing or constraining detection elements for best results. 
///
/// \author Bruce Becker, Javier Castillo
//-----------------------------------------------------------------------------

#include "AliMUONAlignment.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVCluster.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMillepede.h"

#include "AliMpExMap.h"
#include "AliMpExMapIterator.h"

#include "AliAlignObjMatrix.h"
#include "AliLog.h"

#include "TMath.h"
#include "TMatrixDSym.h"

/// \cond CLASSIMP
ClassImp(AliMUONAlignment)
/// \endcond

  Int_t AliMUONAlignment::fgNDetElem = 4*2+4*2+18*2+26*2+26*2;
  Int_t AliMUONAlignment::fgNDetElemCh[10] = {4,4,4,4,18,18,26,26,26,26};
  Int_t AliMUONAlignment::fgSNDetElemCh[10] = {4,8,12,16,34,52,78,104,130,156};
  Int_t AliMUONAlignment::fgNParCh = 4;
  Int_t AliMUONAlignment::fgNTrkMod = 16;
  Int_t AliMUONAlignment::fgNCh = 10;
  Int_t AliMUONAlignment::fgNSt = 5;

AliMUONAlignment::AliMUONAlignment() 
  : TObject(),
    fBFieldOn(kTRUE),
    fStartFac(256.), 
    fResCutInitial(100.), 
    fResCut(100.),
    fMillepede(0),
    fTrackParamAtCluster(0),
    fTrack(0),
    fCluster(0),
    fTrackParam(0),
    fNGlobal(fgNDetElem*fgNParCh),
    fNLocal(4),
    fNStdDev(3),
    fDetElemId(0),
    fDetElemNumber(0),
    fPhi(0.),
    fCosPhi(1.),
    fSinPhi(0.),
    fTransform(0)
{
  /// Default constructor, setting define alignment parameters
  fSigma[0] = 1.0e-1;
  fSigma[1] = 1.0e-2;

  fDoF[0] = kTRUE;  fDoF[1] = kTRUE;  fDoF[2] = kTRUE;  fDoF[3] = kTRUE;
  fAllowVar[0] = 0.05;  fAllowVar[1] = 0.05;  fAllowVar[2] = 0.001;  fAllowVar[3] = 0.5;
  
  AliInfo(Form("fAllowVar[0]: %f\t fAllowVar[1]: %f\t fPhi: %f\t fgNDetElem: %i\t fNGlobal: %i\t fNLocal: %i",fAllowVar[0],fAllowVar[1],fPhi,fgNDetElem,fNGlobal,fNLocal));

  fMillepede = new AliMillepede();

  Init(fNGlobal, fNLocal, fNStdDev);

  ResetLocalEquation();
  AliInfo("Parameters initialized to zero");

}

AliMUONAlignment::~AliMUONAlignment() {
  /// Destructor
}

void AliMUONAlignment::Init(Int_t nGlobal,  /* number of global paramers */
			    Int_t nLocal,   /* number of local parameters */
			    Int_t nStdDev   /* std dev cut */ )
{
  /// Initialization of AliMillepede. Fix parameters, define constraints ...
  fMillepede->InitMille(nGlobal,nLocal,nStdDev,fResCut,fResCutInitial);

//  Bool_t bStOnOff[5] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
//  Bool_t bChOnOff[10] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
//  Bool_t bSpecLROnOff[2] = {kTRUE,kTRUE};

//   AllowVariations(bChOnOff);

  // Fix parameters or add constraints here
//   for (Int_t iSt=0; iSt<5; iSt++)
//     if (!bStOnOff[iSt]) FixStation(iSt+1);
//   for (Int_t iCh=0; iCh<10; iCh++)
//     if (!bChOnOff[iCh]) FixChamber(iCh+1);

//   FixHalfSpectrometer(bChOnOff,bSpecLROnOff);

  ResetConstraints();
  
  // Define global constrains to be applied
  // X, Y, P, XvsZ, YvsZ, PvsZ, XvsY, YvsY, PvsY
  Bool_t bVarXYT[9] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  Bool_t bDetTLBR[4] = {kFALSE,kTRUE,kFALSE,kTRUE};
  //  AddConstraints(bStOnOff,bVarXYT,bDetTLBR,bSpecLROnOff);

  // Other possible way to add constrains
  bVarXYT[0] = kFALSE; bVarXYT[1] = kFALSE; bVarXYT[2] = kTRUE;
  bDetTLBR[0] = kFALSE; bDetTLBR[1] = kTRUE; bDetTLBR[2] = kFALSE; bDetTLBR[3] = kFALSE;
//   AddConstraints(bStOnOff,bVarXYT,bDetTLBR);

  bVarXYT[0] = kTRUE; bVarXYT[1] = kTRUE; bVarXYT[2] = kFALSE;
  //  AddConstraints(bStOnOff,bVarXYT);
  
  // Set iterations
  if (fStartFac>1) fMillepede->SetIterations(fStartFac);          
}

void AliMUONAlignment::FixStation(Int_t iSt){
  /// Fix all detection elements of station iSt
  Int_t iDetElemFirst = (iSt>1) ? fgSNDetElemCh[2*(iSt-1)-1] : 0; 
  Int_t iDetElemLast = fgSNDetElemCh[2*(iSt)-1]; 
  for (Int_t i = iDetElemFirst; i < iDetElemLast; i++){    
    FixParameter(i*fgNParCh+0, 0.0);
    FixParameter(i*fgNParCh+1, 0.0);
    FixParameter(i*fgNParCh+2, 0.0);
    FixParameter(i*fgNParCh+3, 0.0);
  }
}

void AliMUONAlignment::FixChamber(Int_t iCh){
  /// Fix all detection elements of chamber iCh
  Int_t iDetElemFirst = (iCh>1) ? fgSNDetElemCh[iCh-2] : 0; 
  Int_t iDetElemLast = fgSNDetElemCh[iCh-1]; 
  for (Int_t i = iDetElemFirst; i < iDetElemLast; i++){    
    FixParameter(i*fgNParCh+0, 0.0);
    FixParameter(i*fgNParCh+1, 0.0);
    FixParameter(i*fgNParCh+2, 0.0);
    FixParameter(i*fgNParCh+3, 0.0);
  }
}

void AliMUONAlignment::FixDetElem(Int_t iDetElemId, TString sVarXYT){
  /// Fix a given detection element
  Int_t iDetElemNumber = iDetElemId%100;
  for (int iCh=0; iCh<iDetElemId/100-1; iCh++){
    iDetElemNumber += fgNDetElemCh[iCh];
  }
  if (sVarXYT.Contains("X")) { // X constraint
    FixParameter(iDetElemNumber*fgNParCh+0, 0.0);
  }
  if (sVarXYT.Contains("Y")) { // Y constraint
    FixParameter(iDetElemNumber*fgNParCh+1, 0.0);
  }
  if (sVarXYT.Contains("T")) { // T constraint
    FixParameter(iDetElemNumber*fgNParCh+2, 0.0);
  }
  if (sVarXYT.Contains("Z")) { // T constraint
    FixParameter(iDetElemNumber*fgNParCh+3, 0.0);
  }
}

void AliMUONAlignment::FixHalfSpectrometer(const Bool_t *lChOnOff, const Bool_t *lSpecLROnOff){
  /// Fix left or right detector
  for (Int_t i = 0; i < fgNDetElem; i++){    
    Int_t iCh=0;
    for (iCh=1; iCh<=fgNCh; iCh++){
      if (i<fgSNDetElemCh[iCh-1]) break;
    }
    if (lChOnOff[iCh-1]){
      Int_t lDetElemNumber = (iCh==1) ? i : i-fgSNDetElemCh[iCh-2];
      if (iCh>=1 && iCh<=4){
	if ((lDetElemNumber==1 || lDetElemNumber==2) && !lSpecLROnOff[0]){ // From track crossings
	  FixParameter(i*fgNParCh+0, 0.0);
	  FixParameter(i*fgNParCh+1, 0.0);
	  FixParameter(i*fgNParCh+2, 0.0);
	  FixParameter(i*fgNParCh+3, 0.0);
	}
	if ((lDetElemNumber==0 || lDetElemNumber==3) && !lSpecLROnOff[1]){ // From track crossings
	  FixParameter(i*fgNParCh+0, 0.0);
	  FixParameter(i*fgNParCh+1, 0.0);
	  FixParameter(i*fgNParCh+2, 0.0);
	  FixParameter(i*fgNParCh+3, 0.0);
	}
      }
      if (iCh>=5 && iCh<=6){
	if ((lDetElemNumber>=5&&lDetElemNumber<=13) && !lSpecLROnOff[0]){
	  FixParameter(i*fgNParCh+0, 0.0);
	  FixParameter(i*fgNParCh+1, 0.0);
	  FixParameter(i*fgNParCh+2, 0.0);
	  FixParameter(i*fgNParCh+3, 0.0);
	}
	if (((lDetElemNumber>=0&&lDetElemNumber<=4) || 
	     (lDetElemNumber>=14&&lDetElemNumber<=17)) && !lSpecLROnOff[1]){
	  FixParameter(i*fgNParCh+0, 0.0);
	  FixParameter(i*fgNParCh+1, 0.0);
	  FixParameter(i*fgNParCh+2, 0.0);
	  FixParameter(i*fgNParCh+3, 0.0);
	}
      }
      if (iCh>=7 && iCh<=10){
	if ((lDetElemNumber>=7&&lDetElemNumber<=19) && !lSpecLROnOff[0]){
	  FixParameter(i*fgNParCh+0, 0.0);
	  FixParameter(i*fgNParCh+1, 0.0);
	  FixParameter(i*fgNParCh+2, 0.0);
	  FixParameter(i*fgNParCh+3, 0.0);
	}
	if (((lDetElemNumber>=0&&lDetElemNumber<=6) || 
	     (lDetElemNumber>=20&&lDetElemNumber<=25)) && !lSpecLROnOff[1]){
	  FixParameter(i*fgNParCh+0, 0.0);
	  FixParameter(i*fgNParCh+1, 0.0);
	  FixParameter(i*fgNParCh+2, 0.0);
	  FixParameter(i*fgNParCh+3, 0.0);
	}
      }
    }
  }
}

void AliMUONAlignment::SetNonLinear(const Bool_t *lChOnOff, const Bool_t *lVarXYT){
  /// Set non linear parameter flag selected chambers and degrees of freedom
  for (Int_t i = 0; i < fgNDetElem; i++){    
    Int_t iCh=0;
    for (iCh=1; iCh<=fgNCh; iCh++){
      if (i<fgSNDetElemCh[iCh-1]) break;
    }
    if (lChOnOff[iCh-1]){
      if (lVarXYT[0]) { // X constraint
	SetNonLinear(i*fgNParCh+0);
      }
      if (lVarXYT[1]) { // Y constraint
	SetNonLinear(i*fgNParCh+1);
      }
      if (lVarXYT[2]) { // T constraint
	SetNonLinear(i*fgNParCh+2);
      }
      if (lVarXYT[3]) { // Z constraint
	SetNonLinear(i*fgNParCh+3);
      }
    }
  }
}

void AliMUONAlignment::AddConstraints(const Bool_t *lChOnOff, const Bool_t *lVarXYT){
  /// Add constraint equations for selected chambers and degrees of freedom 
  for (Int_t i = 0; i < fgNDetElem; i++){    
    Int_t iCh=0;
    for (iCh=1; iCh<=fgNCh; iCh++){
      if (i<fgSNDetElemCh[iCh-1]) break;
    }
    if (lChOnOff[iCh-1]){
      if (lVarXYT[0]) { // X constraint
	fConstraintX[i*fgNParCh+0]=1.0;
      }
      if (lVarXYT[1]) { // Y constraint
	fConstraintY[i*fgNParCh+1]=1.0;
      }
      if (lVarXYT[2]) { // T constraint
	fConstraintP[i*fgNParCh+2]=1.0;
      }
//       if (lVarXYT[3]) { // Z constraint
// 	fConstraintP[i*fgNParCh+3]=1.0;
//       }
    }
  }
  if (lVarXYT[0]) { // X constraint
    AddConstraint(fConstraintX,0.0);
  }
  if (lVarXYT[1]) { // Y constraint
    AddConstraint(fConstraintY,0.0);
  }
  if (lVarXYT[2]) { // T constraint
    AddConstraint(fConstraintP,0.0);
  }
//   if (lVarXYT[3]) { // Z constraint
//     AddConstraint(fConstraintP,0.0);
//   }
}

void AliMUONAlignment::AddConstraints(const Bool_t *lChOnOff, const Bool_t *lVarXYT, const Bool_t *lDetTLBR, const Bool_t *lSpecLROnOff){
  /// Add constraint equations for selected chambers, degrees of freedom and detector half 
  Double_t lDetElemLocX = 0.;
  Double_t lDetElemLocY = 0.;
  Double_t lDetElemLocZ = 0.;
  Double_t lDetElemGloX = 0.;
  Double_t lDetElemGloY = 0.;
  Double_t lDetElemGloZ = 0.;
  Double_t lMeanY = 0.;
  Double_t lSigmaY = 0.;
  Double_t lMeanZ = 0.;
  Double_t lSigmaZ = 0.;
  Int_t lNDetElem = 0;
  for (Int_t i = 0; i < fgNDetElem; i++){    
    Int_t iCh=0;
    for (iCh=1; iCh<=fgNCh; iCh++){
      if (i<fgSNDetElemCh[iCh-1]) break;
    }
    if (lChOnOff[iCh-1]){ 
      Int_t lDetElemNumber = (iCh==1) ? i : i-fgSNDetElemCh[iCh-2];
      Int_t lDetElemId = iCh*100+lDetElemNumber;
      fTransform->Local2Global(lDetElemId,lDetElemLocX,lDetElemLocY,lDetElemLocZ,
			       lDetElemGloX,lDetElemGloY,lDetElemGloZ);
      if (iCh>=1 && iCh<=4){
	if ((lDetElemNumber==1 || lDetElemNumber==2) && lSpecLROnOff[0]){ // From track crossings
	  lMeanY += lDetElemGloY;
	  lSigmaY += lDetElemGloY*lDetElemGloY;
	  lMeanZ += lDetElemGloZ;
	  lSigmaZ += lDetElemGloZ*lDetElemGloZ;
	  lNDetElem++;
	}
	if ((lDetElemNumber==0 || lDetElemNumber==3) && lSpecLROnOff[1]){ // From track crossings
	  lMeanY += lDetElemGloY;
	  lSigmaY += lDetElemGloY*lDetElemGloY;
	  lMeanZ += lDetElemGloZ;
	  lSigmaZ += lDetElemGloZ*lDetElemGloZ;
	  lNDetElem++;
	}
      }
      if (iCh>=5 && iCh<=6){
	if ((lDetElemNumber>=5&&lDetElemNumber<=13) && lSpecLROnOff[0]){
	  lMeanY += lDetElemGloY;
	  lSigmaY += lDetElemGloY*lDetElemGloY;
	  lMeanZ += lDetElemGloZ;
	  lSigmaZ += lDetElemGloZ*lDetElemGloZ;
	  lNDetElem++;
	}
	if (((lDetElemNumber>=0&&lDetElemNumber<=4) || 
	     (lDetElemNumber>=14&&lDetElemNumber<=17)) && lSpecLROnOff[1]){
	  lMeanY += lDetElemGloY;
	  lSigmaY += lDetElemGloY*lDetElemGloY;
	  lMeanZ += lDetElemGloZ;
	  lSigmaZ += lDetElemGloZ*lDetElemGloZ;
	  lNDetElem++;
	}
      }
      if (iCh>=7 && iCh<=10){
	if ((lDetElemNumber>=7&&lDetElemNumber<=19) && lSpecLROnOff[0]){
	  lMeanY += lDetElemGloY;
	  lSigmaY += lDetElemGloY*lDetElemGloY;
	  lMeanZ += lDetElemGloZ;
	  lSigmaZ += lDetElemGloZ*lDetElemGloZ;
	  lNDetElem++;
	}
	if (((lDetElemNumber>=0&&lDetElemNumber<=6) || 
	     (lDetElemNumber>=20&&lDetElemNumber<=25)) && lSpecLROnOff[1]){
	  lMeanY += lDetElemGloY;
	  lSigmaY += lDetElemGloY*lDetElemGloY;
	  lMeanZ += lDetElemGloZ;
	  lSigmaZ += lDetElemGloZ*lDetElemGloZ;
	  lNDetElem++;
	}
      }
    }
  }
  lMeanY /= lNDetElem;
  lSigmaY /= lNDetElem;
  lSigmaY = TMath::Sqrt(lSigmaY-lMeanY*lMeanY);
  lMeanZ /= lNDetElem;
  lSigmaZ /= lNDetElem;
  lSigmaZ = TMath::Sqrt(lSigmaZ-lMeanZ*lMeanZ);
  AliInfo(Form("Used %i DetElem, MeanZ= %f , SigmaZ= %f", lNDetElem,lMeanZ,lSigmaZ));  

  for (Int_t i = 0; i < fgNDetElem; i++){    
    Int_t iCh=0;
    for (iCh=1; iCh<=fgNCh; iCh++){
      if (i<fgSNDetElemCh[iCh-1]) break;
    }
    if (lChOnOff[iCh-1]){
      Int_t lDetElemNumber = (iCh==1) ? i : i-fgSNDetElemCh[iCh-2];
      Int_t lDetElemId = iCh*100+lDetElemNumber;
      fTransform->Local2Global(lDetElemId,lDetElemLocX,lDetElemLocY,lDetElemLocZ,
			       lDetElemGloX,lDetElemGloY,lDetElemGloZ);
      if (lVarXYT[0]) { // X constraint
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintXT,0); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintXL,0); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintXB,0); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintXR,0); // Right half
      }
      if (lVarXYT[1]) { // Y constraint
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintYT,1); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintYL,1); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintYB,1); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintYR,1); // Right half
      }
      if (lVarXYT[2]) { // P constraint
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintPT,2); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintPL,2); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintPB,2); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintPR,2); // Right half
      }
      if (lVarXYT[3]) { // X-Z shearing
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintXZT,0,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintXZL,0,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintXZB,0,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintXZR,0,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Right half
      }
      if (lVarXYT[4]) { // Y-Z shearing
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintYZT,1,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintYZL,1,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintYZB,1,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintYZR,1,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Right half
      }
      if (lVarXYT[5]) { // P-Z rotation
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintPZT,2,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintPZL,2,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintPZB,2,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintPZR,2,(lDetElemGloZ-lMeanZ)/lSigmaZ); // Right half
      }
      if (lVarXYT[6]) { // X-Y shearing
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintXYT,0,(lDetElemGloY-lMeanY)/lSigmaY); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintXYL,0,(lDetElemGloY-lMeanY)/lSigmaY); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintXYB,0,(lDetElemGloY-lMeanY)/lSigmaY); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintXYR,0,(lDetElemGloY-lMeanY)/lSigmaY); // Right half
      }
      if (lVarXYT[7]) { // Y-Y scaling
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintYYT,1,(lDetElemGloY-lMeanY)/lSigmaY); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintYYL,1,(lDetElemGloY-lMeanY)/lSigmaY); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintYYB,1,(lDetElemGloY-lMeanY)/lSigmaY); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintYYR,1,(lDetElemGloY-lMeanY)/lSigmaY); // Right half
      }
      if (lVarXYT[8]) { // P-Y rotation
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintPYT,2,(lDetElemGloY-lMeanY)/lSigmaY); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintPYL,2,(lDetElemGloY-lMeanY)/lSigmaY); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintPYB,2,(lDetElemGloY-lMeanY)/lSigmaY); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintPYR,2,(lDetElemGloY-lMeanY)/lSigmaY); // Right half
      }
    }
  }
  if (lVarXYT[0]) { // X constraint
    if (lDetTLBR[0]) AddConstraint(fConstraintXT,0.0); // Top half
    if (lDetTLBR[1]) AddConstraint(fConstraintXL,0.0); // Left half
    if (lDetTLBR[2]) AddConstraint(fConstraintXB,0.0); // Bottom half
    if (lDetTLBR[3]) AddConstraint(fConstraintXR,0.0); // Right half
  }
  if (lVarXYT[1]) { // Y constraint
    if (lDetTLBR[0]) AddConstraint(fConstraintYT,0.0); // Top half
    if (lDetTLBR[1]) AddConstraint(fConstraintYL,0.0); // Left half
    if (lDetTLBR[2]) AddConstraint(fConstraintYB,0.0); // Bottom half
    if (lDetTLBR[3]) AddConstraint(fConstraintYR,0.0); // Right half
  }
  if (lVarXYT[2]) { // T constraint
    if (lDetTLBR[0]) AddConstraint(fConstraintPT,0.0); // Top half
    if (lDetTLBR[1]) AddConstraint(fConstraintPL,0.0); // Left half
    if (lDetTLBR[2]) AddConstraint(fConstraintPB,0.0); // Bottom half
    if (lDetTLBR[3]) AddConstraint(fConstraintPR,0.0); // Right half
  }
  if (lVarXYT[3]) { // X-Z constraint
    if (lDetTLBR[0]) AddConstraint(fConstraintXZT,0.0); // Top half
    if (lDetTLBR[1]) AddConstraint(fConstraintXZL,0.0); // Left half
    if (lDetTLBR[2]) AddConstraint(fConstraintXZB,0.0); // Bottom half
    if (lDetTLBR[3]) AddConstraint(fConstraintXZR,0.0); // Right half
  }
  if (lVarXYT[4]) { // Y-Z constraint
    if (lDetTLBR[0]) AddConstraint(fConstraintYZT,0.0); // Top half
    if (lDetTLBR[1]) AddConstraint(fConstraintYZL,0.0); // Left half
    if (lDetTLBR[2]) AddConstraint(fConstraintYZB,0.0); // Bottom half
    if (lDetTLBR[3]) AddConstraint(fConstraintYZR,0.0); // Right half
  }
  if (lVarXYT[5]) { // P-Z constraint
    if (lDetTLBR[0]) AddConstraint(fConstraintPZT,0.0); // Top half
    if (lDetTLBR[1]) AddConstraint(fConstraintPZL,0.0); // Left half
    if (lDetTLBR[2]) AddConstraint(fConstraintPZB,0.0); // Bottom half
    if (lDetTLBR[3]) AddConstraint(fConstraintPZR,0.0); // Right half
  }
  if (lVarXYT[6]) { // X-Y constraint
    if (lDetTLBR[0]) AddConstraint(fConstraintXYT,0.0); // Top half
    if (lDetTLBR[1]) AddConstraint(fConstraintXYL,0.0); // Left half
    if (lDetTLBR[2]) AddConstraint(fConstraintXYB,0.0); // Bottom half
    if (lDetTLBR[3]) AddConstraint(fConstraintXYR,0.0); // Right half
  }
  if (lVarXYT[7]) { // Y-Y constraint
    if (lDetTLBR[0]) AddConstraint(fConstraintYYT,0.0); // Top half
    if (lDetTLBR[1]) AddConstraint(fConstraintYYL,0.0); // Left half
    if (lDetTLBR[2]) AddConstraint(fConstraintYYB,0.0); // Bottom half
    if (lDetTLBR[3]) AddConstraint(fConstraintYYR,0.0); // Right half
  }
  if (lVarXYT[8]) { // P-Y constraint
    if (lDetTLBR[0]) AddConstraint(fConstraintPYT,0.0); // Top half
    if (lDetTLBR[1]) AddConstraint(fConstraintPYL,0.0); // Left half
    if (lDetTLBR[2]) AddConstraint(fConstraintPYB,0.0); // Bottom half
    if (lDetTLBR[3]) AddConstraint(fConstraintPYR,0.0); // Right half
  }
}

void AliMUONAlignment::ConstrainT(Int_t lDetElem, Int_t lCh, Double_t *lConstraintT, Int_t iVar, Double_t /*lWeight*/) const{
  /// Set constrain equation for top half of spectrometer
  Int_t lDetElemNumber = (lCh==1) ? lDetElem : lDetElem-fgSNDetElemCh[lCh-2];
  if (lCh>=1 && lCh<=4){
    if (lDetElemNumber==0 || lDetElemNumber==1){ // From track crossings
      lConstraintT[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
  if (lCh>=5 && lCh<=6){
    if (lDetElemNumber>=0&&lDetElemNumber<=9){
      lConstraintT[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
  if (lCh>=7 && lCh<=10){
    if (lDetElemNumber>=0&&lDetElemNumber<=13){
      lConstraintT[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
}

void AliMUONAlignment::ConstrainL(Int_t lDetElem, Int_t lCh, Double_t *lConstraintL, Int_t iVar, Double_t lWeight) const{
  /// Set constrain equation for left half of spectrometer
  Int_t lDetElemNumber = (lCh==1) ? lDetElem : lDetElem-fgSNDetElemCh[lCh-2];
  if (lCh>=1 && lCh<=4){
    if (lDetElemNumber==1 || lDetElemNumber==2){ // From track crossings
      lConstraintL[lDetElem*fgNParCh+iVar]=lWeight;
    }
  }
  if (lCh>=5 && lCh<=6){
    if (lDetElemNumber>=5&&lDetElemNumber<=13){
      lConstraintL[lDetElem*fgNParCh+iVar]=lWeight;
    }
  }
  if (lCh>=7 && lCh<=10){
    if (lDetElemNumber>=7&&lDetElemNumber<=19){
      lConstraintL[lDetElem*fgNParCh+iVar]=lWeight;
    }
  }
}

void AliMUONAlignment::ConstrainB(Int_t lDetElem, Int_t lCh, Double_t *lConstraintB, Int_t iVar, Double_t /*lWeight*/) const{
  /// Set constrain equation for bottom half of spectrometer
  Int_t lDetElemNumber = (lCh==1) ? lDetElem : lDetElem-fgSNDetElemCh[lCh-2];
  if (lCh>=1 && lCh<=4){
    if (lDetElemNumber==2 || lDetElemNumber==3){ // From track crossings
      lConstraintB[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
  if (lCh>=5 && lCh<=6){
    if ((lDetElemNumber>=9&&lDetElemNumber<=17) || 
	(lDetElemNumber==0)){
      lConstraintB[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
  if (lCh>=7 && lCh<=10){
    if ((lDetElemNumber>=13&&lDetElemNumber<=25) || 
	(lDetElemNumber==0)){
      lConstraintB[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
}

void AliMUONAlignment::ConstrainR(Int_t lDetElem, Int_t lCh, Double_t *lConstraintR, Int_t iVar, Double_t lWeight) const{
  /// Set constrain equation for right half of spectrometer
  Int_t lDetElemNumber = (lCh==1) ? lDetElem : lDetElem-fgSNDetElemCh[lCh-2];
  if (lCh>=1 && lCh<=4){
    if (lDetElemNumber==0 || lDetElemNumber==3){ // From track crossings
      lConstraintR[lDetElem*fgNParCh+iVar]=lWeight;
    }
  }
  if (lCh>=5 && lCh<=6){
    if ((lDetElemNumber>=0&&lDetElemNumber<=4) || 
	(lDetElemNumber>=14&&lDetElemNumber<=17)){
      lConstraintR[lDetElem*fgNParCh+iVar]=lWeight;
    }
  }
  if (lCh>=7 && lCh<=10){
    if ((lDetElemNumber>=0&&lDetElemNumber<=6) || 
	(lDetElemNumber>=20&&lDetElemNumber<=25)){
      lConstraintR[lDetElem*fgNParCh+iVar]=lWeight;
    }
  }
}

void AliMUONAlignment::ResetConstraints(){
  /// Reset all constraint equations
  for (Int_t i = 0; i < fgNDetElem; i++){    
    fConstraintX[i*fgNParCh+0]=0.0;
    fConstraintX[i*fgNParCh+1]=0.0;
    fConstraintX[i*fgNParCh+2]=0.0;
    fConstraintY[i*fgNParCh+0]=0.0;
    fConstraintY[i*fgNParCh+1]=0.0;
    fConstraintY[i*fgNParCh+2]=0.0;
    fConstraintP[i*fgNParCh+0]=0.0;
    fConstraintP[i*fgNParCh+1]=0.0;
    fConstraintP[i*fgNParCh+2]=0.0;
    fConstraintXT[i*fgNParCh+0]=0.0;
    fConstraintXT[i*fgNParCh+1]=0.0;
    fConstraintXT[i*fgNParCh+2]=0.0;
    fConstraintYT[i*fgNParCh+0]=0.0;
    fConstraintYT[i*fgNParCh+1]=0.0;
    fConstraintYT[i*fgNParCh+2]=0.0;
    fConstraintPT[i*fgNParCh+0]=0.0;
    fConstraintPT[i*fgNParCh+1]=0.0;
    fConstraintPT[i*fgNParCh+2]=0.0;
    fConstraintXZT[i*fgNParCh+0]=0.0;
    fConstraintXZT[i*fgNParCh+1]=0.0;
    fConstraintXZT[i*fgNParCh+2]=0.0;
    fConstraintYZT[i*fgNParCh+0]=0.0;
    fConstraintYZT[i*fgNParCh+1]=0.0;
    fConstraintYZT[i*fgNParCh+2]=0.0;
    fConstraintPZT[i*fgNParCh+0]=0.0;
    fConstraintPZT[i*fgNParCh+1]=0.0;
    fConstraintPZT[i*fgNParCh+2]=0.0;
    fConstraintXYT[i*fgNParCh+0]=0.0;
    fConstraintXYT[i*fgNParCh+1]=0.0;
    fConstraintXYT[i*fgNParCh+2]=0.0;
    fConstraintYYT[i*fgNParCh+0]=0.0;
    fConstraintYYT[i*fgNParCh+1]=0.0;
    fConstraintYYT[i*fgNParCh+2]=0.0;
    fConstraintPYT[i*fgNParCh+0]=0.0;
    fConstraintPYT[i*fgNParCh+1]=0.0;
    fConstraintPYT[i*fgNParCh+2]=0.0;
    fConstraintXL[i*fgNParCh+0]=0.0;
    fConstraintXL[i*fgNParCh+1]=0.0;
    fConstraintXL[i*fgNParCh+2]=0.0;
    fConstraintYL[i*fgNParCh+0]=0.0;
    fConstraintYL[i*fgNParCh+1]=0.0;
    fConstraintYL[i*fgNParCh+2]=0.0;
    fConstraintPL[i*fgNParCh+0]=0.0;
    fConstraintPL[i*fgNParCh+1]=0.0;
    fConstraintPL[i*fgNParCh+2]=0.0;
    fConstraintXZL[i*fgNParCh+0]=0.0;
    fConstraintXZL[i*fgNParCh+1]=0.0;
    fConstraintXZL[i*fgNParCh+2]=0.0;
    fConstraintYZL[i*fgNParCh+0]=0.0;
    fConstraintYZL[i*fgNParCh+1]=0.0;
    fConstraintYZL[i*fgNParCh+2]=0.0;
    fConstraintPZL[i*fgNParCh+0]=0.0;
    fConstraintPZL[i*fgNParCh+1]=0.0;
    fConstraintPZL[i*fgNParCh+2]=0.0;
    fConstraintXYL[i*fgNParCh+0]=0.0;
    fConstraintXYL[i*fgNParCh+1]=0.0;
    fConstraintXYL[i*fgNParCh+2]=0.0;
    fConstraintYYL[i*fgNParCh+0]=0.0;
    fConstraintYYL[i*fgNParCh+1]=0.0;
    fConstraintYYL[i*fgNParCh+2]=0.0;
    fConstraintPYL[i*fgNParCh+0]=0.0;
    fConstraintPYL[i*fgNParCh+1]=0.0;
    fConstraintPYL[i*fgNParCh+2]=0.0;
    fConstraintXB[i*fgNParCh+0]=0.0;
    fConstraintXB[i*fgNParCh+1]=0.0;
    fConstraintXB[i*fgNParCh+2]=0.0;
    fConstraintYB[i*fgNParCh+0]=0.0;
    fConstraintYB[i*fgNParCh+1]=0.0;
    fConstraintYB[i*fgNParCh+2]=0.0;
    fConstraintPB[i*fgNParCh+0]=0.0;
    fConstraintPB[i*fgNParCh+1]=0.0;
    fConstraintPB[i*fgNParCh+2]=0.0;
    fConstraintXZB[i*fgNParCh+0]=0.0;
    fConstraintXZB[i*fgNParCh+1]=0.0;
    fConstraintXZB[i*fgNParCh+2]=0.0;
    fConstraintYZB[i*fgNParCh+0]=0.0;
    fConstraintYZB[i*fgNParCh+1]=0.0;
    fConstraintYZB[i*fgNParCh+2]=0.0;
    fConstraintPZB[i*fgNParCh+0]=0.0;
    fConstraintPZB[i*fgNParCh+1]=0.0;
    fConstraintPZB[i*fgNParCh+2]=0.0;
    fConstraintXYB[i*fgNParCh+0]=0.0;
    fConstraintXYB[i*fgNParCh+1]=0.0;
    fConstraintXYB[i*fgNParCh+2]=0.0;
    fConstraintYYB[i*fgNParCh+0]=0.0;
    fConstraintYYB[i*fgNParCh+1]=0.0;
    fConstraintYYB[i*fgNParCh+2]=0.0;
    fConstraintPYB[i*fgNParCh+0]=0.0;
    fConstraintPYB[i*fgNParCh+1]=0.0;
    fConstraintPYB[i*fgNParCh+2]=0.0;
    fConstraintXR[i*fgNParCh+0]=0.0;
    fConstraintXR[i*fgNParCh+1]=0.0;
    fConstraintXR[i*fgNParCh+2]=0.0;
    fConstraintYR[i*fgNParCh+0]=0.0;
    fConstraintYR[i*fgNParCh+1]=0.0;
    fConstraintYR[i*fgNParCh+2]=0.0;
    fConstraintPR[i*fgNParCh+0]=0.0;
    fConstraintPR[i*fgNParCh+1]=0.0;
    fConstraintPR[i*fgNParCh+2]=0.0;
    fConstraintXZR[i*fgNParCh+0]=0.0;
    fConstraintXZR[i*fgNParCh+1]=0.0;
    fConstraintXZR[i*fgNParCh+2]=0.0;
    fConstraintYZR[i*fgNParCh+0]=0.0;
    fConstraintYZR[i*fgNParCh+1]=0.0;
    fConstraintYZR[i*fgNParCh+2]=0.0;
    fConstraintPZR[i*fgNParCh+0]=0.0;
    fConstraintPZR[i*fgNParCh+1]=0.0;
    fConstraintPZR[i*fgNParCh+2]=0.0;
    fConstraintPZR[i*fgNParCh+0]=0.0;
    fConstraintPZR[i*fgNParCh+1]=0.0;
    fConstraintPZR[i*fgNParCh+2]=0.0;
    fConstraintXYR[i*fgNParCh+0]=0.0;
    fConstraintXYR[i*fgNParCh+1]=0.0;
    fConstraintXYR[i*fgNParCh+2]=0.0;
    fConstraintYYR[i*fgNParCh+0]=0.0;
    fConstraintYYR[i*fgNParCh+1]=0.0;
    fConstraintYYR[i*fgNParCh+2]=0.0;
    fConstraintPYR[i*fgNParCh+0]=0.0;
    fConstraintPYR[i*fgNParCh+1]=0.0;
    fConstraintPYR[i*fgNParCh+2]=0.0;
  }
}

void AliMUONAlignment::AddConstraint(Double_t *par, Double_t value) {
  /// Constrain equation defined by par to value
  fMillepede->SetGlobalConstraint(par, value);
  AliInfo("Adding constraint");
}

void AliMUONAlignment::InitGlobalParameters(Double_t *par) {
  /// Initialize global parameters with par array
  fMillepede->SetGlobalParameters(par);
  AliInfo("Init Global Parameters");
}
 
void AliMUONAlignment::FixParameter(Int_t iPar, Double_t value) {
  /// Parameter iPar is encourage to vary in [-value;value]. 
  /// If value == 0, parameter is fixed
  fMillepede->SetParSigma(iPar, value);
  if (TMath::Abs(value)<1e-4) AliInfo(Form("Parameter %i Fixed", iPar));
}

void AliMUONAlignment::ResetLocalEquation()
{
  /// Reset the derivative vectors
  for(int i=0; i<fNLocal; i++) {
    fLocalDerivatives[i] = 0.0;
  }
  for(int i=0; i<fNGlobal; i++) {
    fGlobalDerivatives[i] = 0.0;
  }
}

void AliMUONAlignment::AllowVariations(const Bool_t *bChOnOff) {
  /// Set allowed variation for selected chambers based on fDoF and fAllowVar
  for (Int_t iCh=1; iCh<=10; iCh++) {
    if (bChOnOff[iCh-1]) {
      Int_t iDetElemFirst = (iCh>1) ? fgSNDetElemCh[iCh-2] : 0; 
      Int_t iDetElemLast = fgSNDetElemCh[iCh-1]; 
      for (int i=0; i<fgNParCh; i++) {
	AliDebug(1,Form("fDoF[%d]= %d",i,fDoF[i]));    
	if (fDoF[i]) {
	  for (Int_t j=iDetElemFirst; j<iDetElemLast; j++){    
	    FixParameter(j*fgNParCh+i, fAllowVar[i]);
	  }
	}
      }
    }
  }
}

void AliMUONAlignment::SetNonLinear(Int_t iPar  /* set non linear flag */ ) {
  /// Set nonlinear flag for parameter iPar
  fMillepede->SetNonLinear(iPar);
  AliInfo(Form("Parameter %i set to non linear", iPar));
}

void AliMUONAlignment::LocalEquationX() {
  /// Define local equation for current track and cluster in x coor. measurement
  // set local derivatives
  SetLocalDerivative(0, fCosPhi);
  SetLocalDerivative(1, fCosPhi * (fTrackPos[2] - fTrackPos0[2]));
  SetLocalDerivative(2, fSinPhi);
  SetLocalDerivative(3, fSinPhi * (fTrackPos[2] - fTrackPos0[2]));

  // set global derivatives
  SetGlobalDerivative(fDetElemNumber*fgNParCh+0, -fCosPhi);
  SetGlobalDerivative(fDetElemNumber*fgNParCh+1, -fSinPhi);
  if (fBFieldOn){
    SetGlobalDerivative(fDetElemNumber*fgNParCh+2,
			-fSinPhi*(fTrackPos[0]-fDetElemPos[0]) 
			+fCosPhi*(fTrackPos[1]-fDetElemPos[1]));
  }
  else {
    SetGlobalDerivative(fDetElemNumber*fgNParCh+2,
			-fSinPhi*(fTrackPos0[0]+fTrackSlope0[0]*
				  (fTrackPos[2]-fTrackPos0[2])-fDetElemPos[0]) 
			+fCosPhi*(fTrackPos0[1]+fTrackSlope0[1]*
				  (fTrackPos[2]-fTrackPos0[2])-fDetElemPos[1]));
  }
  SetGlobalDerivative(fDetElemNumber*fgNParCh+3, 
		      fCosPhi*fTrackSlope0[0]+fSinPhi*fTrackSlope0[1]);

  fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, fMeas[0], fSigma[0]);
}

void AliMUONAlignment::LocalEquationY() {
  /// Define local equation for current track and cluster in y coor. measurement
  // set local derivatives
  SetLocalDerivative(0,-fSinPhi);
  SetLocalDerivative(1,-fSinPhi * (fTrackPos[2] - fTrackPos0[2]));
  SetLocalDerivative(2, fCosPhi);
  SetLocalDerivative(3, fCosPhi * (fTrackPos[2] - fTrackPos0[2]));

  // set global derivatives
  SetGlobalDerivative(fDetElemNumber*fgNParCh+0,  fSinPhi);
  SetGlobalDerivative(fDetElemNumber*fgNParCh+1, -fCosPhi);
  if (fBFieldOn){
  SetGlobalDerivative(fDetElemNumber*fgNParCh+2,
 		      -fCosPhi*(fTrackPos[0]-fDetElemPos[0])
 		      -fSinPhi*(fTrackPos[1]-fDetElemPos[1]));
  }
  else {
    SetGlobalDerivative(fDetElemNumber*fgNParCh+2,
			-fCosPhi*(fTrackPos0[0]+fTrackSlope0[0]*
				  (fTrackPos[2]-fTrackPos0[2])-fDetElemPos[0])
			-fSinPhi*(fTrackPos0[1]+fTrackSlope0[1]*
				  (fTrackPos[2]-fTrackPos0[2])-fDetElemPos[1]));
  }
  SetGlobalDerivative(fDetElemNumber*fgNParCh+3,
		      -fSinPhi*fTrackSlope0[0]+fCosPhi*fTrackSlope0[1]);

  fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, fMeas[1], fSigma[1]);
}

void AliMUONAlignment::FillRecPointData() {
  /// Get information of current cluster
  fClustPos[0] = fCluster->GetX();
  fClustPos[1] = fCluster->GetY();
  fClustPos[2] = fCluster->GetZ();
  fTransform->Global2Local(fDetElemId,fClustPos[0],fClustPos[1],fClustPos[2],
			    fClustPosLoc[0],fClustPosLoc[1],fClustPosLoc[2]);
}

void AliMUONAlignment::FillTrackParamData() {
  /// Get information of current track at current cluster
  fTrackPos[0] = fTrackParam->GetNonBendingCoor();
  fTrackPos[1] = fTrackParam->GetBendingCoor();
  fTrackPos[2] = fTrackParam->GetZ();
  fTrackSlope[0] = fTrackParam->GetNonBendingSlope();
  fTrackSlope[1] = fTrackParam->GetBendingSlope();
  fTransform->Global2Local(fDetElemId,fTrackPos[0],fTrackPos[1],fTrackPos[2],
			    fTrackPosLoc[0],fTrackPosLoc[1],fTrackPosLoc[2]);
}

void AliMUONAlignment::FillDetElemData() {
  /// Get information of current detection element
  Double_t lDetElemLocX = 0.;
  Double_t lDetElemLocY = 0.;
  Double_t lDetElemLocZ = 0.;
  fDetElemId = fCluster->GetDetElemId();
  fDetElemNumber = fDetElemId%100;
  for (int iCh=0; iCh<fDetElemId/100-1; iCh++){
    fDetElemNumber += fgNDetElemCh[iCh];
  }
  fTransform->Local2Global(fDetElemId,lDetElemLocX,lDetElemLocY,lDetElemLocZ,
			   fDetElemPos[0],fDetElemPos[1],fDetElemPos[2]);
}

void AliMUONAlignment::ProcessTrack(AliMUONTrack * track) {
  /// Process track; Loop over clusters and set local equations
  fTrack = track;
  // get tclones arrays.
  fTrackParamAtCluster = fTrack->GetTrackParamAtCluster();
  
  // get size of arrays
  Int_t nTrackParam = fTrackParamAtCluster->GetEntries();
  AliDebug(1,Form("Number of track param entries : %i ", nTrackParam));

  for(Int_t iCluster=0; iCluster<nTrackParam; iCluster++) {
    fTrackParam = (AliMUONTrackParam *) fTrack->GetTrackParamAtCluster()->At(iCluster);
    fCluster = fTrackParam->GetClusterPtr();
    if (!fCluster || !fTrackParam) continue;
    // fill local variables for this position --> one measurement
    FillDetElemData();
    FillRecPointData();
    FillTrackParamData();         
//     if (fDetElemId<500) continue;
    fTrackPos0[0]      = fTrackPos[0]; 	  
    fTrackPos0[1]      = fTrackPos[1]; 	  
    fTrackPos0[2]      = fTrackPos[2]; 	  
    fTrackSlope0[0] = fTrackSlope[0]; 
    fTrackSlope0[1] = fTrackSlope[1];   
    break;
  }

  for(Int_t iCluster=0; iCluster<nTrackParam; iCluster++) {
    // and get new pointers
    fTrackParam = (AliMUONTrackParam *) fTrack->GetTrackParamAtCluster()->At(iCluster);
    fCluster = fTrackParam->GetClusterPtr();
    if (!fCluster || !fTrackParam) continue;
    // fill local variables for this position --> one measurement
    FillDetElemData();        
    FillRecPointData();
    FillTrackParamData();
//     if (fDetElemId<500) continue;
    AliDebug(1,Form("cluster: %i", iCluster));
    AliDebug(1,Form("x: %f\t y: %f\t z: %f\t DetElemID: %i\t ", fClustPos[0], fClustPos[1], fClustPos[2], fDetElemId));
    AliDebug(1,Form("fDetElemPos[0]: %f\t fDetElemPos[1]: %f\t fDetElemPos[2]: %f\t DetElemID: %i\t ", fDetElemPos[0],fDetElemPos[1],fDetElemPos[2], fDetElemId));

    AliDebug(1,Form("Track Parameter: %i", iCluster));
    AliDebug(1,Form("x: %f\t y: %f\t z: %f\t slopex: %f\t slopey: %f", fTrackPos[0], fTrackPos[1], fTrackPos[2], fTrackSlope[0], fTrackSlope[1]));
    AliDebug(1,Form("x0: %f\t y0: %f\t z0: %f\t slopex0: %f\t slopey0: %f", fTrackPos0[0], fTrackPos0[1], fTrackPos0[2], fTrackSlope0[0], fTrackSlope0[1]));
    
    fCosPhi = TMath::Cos(fPhi);
    fSinPhi = TMath::Sin(fPhi);
    if (fBFieldOn){
      fMeas[0] = fTrackPos[0] - fClustPos[0];
      fMeas[1] = fTrackPos[1] - fClustPos[1];
    }
    else {
      fMeas[0] = - fClustPos[0];
      fMeas[1] = - fClustPos[1];
    }
    AliDebug(1,Form("fMeas[0]: %f\t fMeas[1]: %f\t fSigma[0]: %f\t fSigma[1]: %f", fMeas[0], fMeas[1], fSigma[0], fSigma[1]));    
    // Set local equations
    LocalEquationX();
    LocalEquationY();
  }
}

void AliMUONAlignment::LocalFit(Int_t iTrack, Double_t *lTrackParam, Int_t lSingleFit) {
  /// Call local fit for this tracks
  Int_t iRes = fMillepede->LocalFit(iTrack,lTrackParam,lSingleFit);
  if (iRes && !lSingleFit) {
    fMillepede->SetNLocalEquations(fMillepede->GetNLocalEquations()+1);
  }
}

void AliMUONAlignment::GlobalFit(Double_t *parameters,Double_t *errors,Double_t *pulls) {
  /// Call global fit; Global parameters are stored in parameters
  fMillepede->GlobalFit(parameters,errors,pulls);

  AliInfo("Done fitting global parameters!");
  for (int iGlob=0; iGlob<fgNDetElem; iGlob++){
    printf("%d\t %f\t %f\t %f \n",iGlob,parameters[iGlob*fgNParCh+0],parameters[iGlob*fgNParCh+1],parameters[iGlob*fgNParCh+2]);
  }
}

Double_t AliMUONAlignment::GetParError(Int_t iPar) {
  /// Get error of parameter iPar
  Double_t lErr = fMillepede->GetParError(iPar);
  return lErr;
}

void AliMUONAlignment::PrintGlobalParameters() {
  /// Print global parameters
  fMillepede->PrintGlobalParameters();
}

//_________________________________________________________________________
TGeoCombiTrans AliMUONAlignment::ReAlign(const TGeoCombiTrans & transform, const double *lMisAlignment) const
{
  /// Realign given transformation by given misalignment and return the misaligned transformation
  
  Double_t cartMisAlig[3] = {0,0,0};
  Double_t angMisAlig[3] = {0,0,0};
//   const Double_t *trans = transform.GetTranslation();
//   TGeoRotation *rot;
//   // check if the rotation we obtain is not NULL
//   if (transform.GetRotation()) {    
//     rot = transform.GetRotation();
//   }
//   else {    
//     rot = new TGeoRotation("rot");
//   }			// default constructor.

  cartMisAlig[0] = -TMath::Sign(1.0,transform.GetRotationMatrix()[0])*lMisAlignment[0];
  cartMisAlig[1] = -TMath::Sign(1.0,transform.GetRotationMatrix()[4])*lMisAlignment[1];
  cartMisAlig[2] = -TMath::Sign(1.0,transform.GetRotationMatrix()[8])*lMisAlignment[3];
  angMisAlig[2] = -TMath::Sign(1.0,transform.GetRotationMatrix()[0]*transform.GetRotationMatrix()[4])*lMisAlignment[2]*180./TMath::Pi();

  TGeoTranslation deltaTrans(cartMisAlig[0], cartMisAlig[1], cartMisAlig[2]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(angMisAlig[0]);
  deltaRot.RotateY(angMisAlig[1]);
  deltaRot.RotateZ(angMisAlig[2]);

  TGeoCombiTrans deltaTransf(deltaTrans,deltaRot);
  TGeoHMatrix newTransfMat = transform * deltaTransf;

  return TGeoCombiTrans(newTransfMat);
}

//______________________________________________________________________
AliMUONGeometryTransformer *
AliMUONAlignment::ReAlign(const AliMUONGeometryTransformer * transformer,
			    const double *misAlignments, Bool_t verbose)
			    
{
  /// Returns a new AliMUONGeometryTransformer with the found misalignments
  /// applied. 

  // Takes the internal geometry module transformers, copies them
  // and gets the Detection Elements from them.
  // Takes misalignment parameters and applies these
  // to the local transform of the Detection Element
  // Obtains the global transform by multiplying the module transformer
  // transformation with the local transformation 
  // Applies the global transform to a new detection element
  // Adds the new detection element to a new module transformer
  // Adds the new module transformer to a new geometry transformer
  // Returns the new geometry transformer

  Double_t lModuleMisAlignment[3] = {0.,0.,0.};
  Double_t lDetElemMisAlignment[4] = {0.,0.,0.,0.};
  Int_t iDetElemId = 0;
  Int_t iDetElemNumber = 0;

  AliMUONGeometryTransformer *newGeometryTransformer =
    new AliMUONGeometryTransformer();
  for (Int_t iMt = 0; iMt < transformer->GetNofModuleTransformers(); iMt++) {
    // module transformers    
    const AliMUONGeometryModuleTransformer *kModuleTransformer =
      transformer->GetModuleTransformer(iMt, true);
      
    AliMUONGeometryModuleTransformer *newModuleTransformer =
      new AliMUONGeometryModuleTransformer(iMt);
    newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
    
    TGeoCombiTrans moduleTransform =
      TGeoCombiTrans(*kModuleTransformer->GetTransformation());
    // New module transformation
    TGeoCombiTrans newModuleTransform = ReAlign(moduleTransform,lModuleMisAlignment);
    newModuleTransformer->SetTransformation(newModuleTransform);
    
    // Get delta transformation: 
    // Tdelta = Tnew * Told.inverse
    TGeoHMatrix deltaModuleTransform = 
      AliMUONGeometryBuilder::Multiply(newModuleTransform, 
				       kModuleTransformer->GetTransformation()->Inverse());    
    // Create module mis alignment matrix
    newGeometryTransformer
      ->AddMisAlignModule(kModuleTransformer->GetModuleId(), deltaModuleTransform);
    
    AliMpExMap *detElements = kModuleTransformer->GetDetElementStore();
    
    if (verbose)
      AliInfo(Form("%i DEs in old GeometryStore  %i",detElements->GetSize(), iMt));

    TIter next(detElements->CreateIterator());
    AliMUONGeometryDetElement* detElement;
    Int_t iDe(-1);
    while ( ( detElement = static_cast<AliMUONGeometryDetElement*>(next()) ) )
    {
      ++iDe;
      // make a new detection element
      AliMUONGeometryDetElement *newDetElement =
	new AliMUONGeometryDetElement(detElement->GetId(),
				      detElement->GetVolumePath());
      TString lDetElemName(detElement->GetDEName());
      lDetElemName.ReplaceAll("DE","");
      iDetElemId = lDetElemName.Atoi();
      iDetElemNumber = iDetElemId%100;
      for (int iCh=0; iCh<iDetElemId/100-1; iCh++){
	iDetElemNumber += fgNDetElemCh[iCh];
      }
      for (int i=0; i<fgNParCh; i++) {
	lDetElemMisAlignment[i] = 0.0;
	if (iMt<fgNTrkMod) {
	  AliInfo(Form("iMt %i, iCh %i, iDe %i, iDeId %i, iDeNb %i, iPar %i",iMt, iDetElemId/100, iDe, iDetElemId, iDetElemNumber, iDetElemNumber*fgNParCh+i));
	  lDetElemMisAlignment[i] =  misAlignments[iDetElemNumber*fgNParCh+i];
	}
      }
      // local transformation of this detection element.
      TGeoCombiTrans localTransform
	= TGeoCombiTrans(*detElement->GetLocalTransformation());
      TGeoCombiTrans newLocalTransform = ReAlign(localTransform,lDetElemMisAlignment);
      newDetElement->SetLocalTransformation(newLocalTransform);	  

      // global transformation
      TGeoHMatrix newGlobalTransform =
	AliMUONGeometryBuilder::Multiply(newModuleTransform,
					 newLocalTransform);
      newDetElement->SetGlobalTransformation(newGlobalTransform);
	  
      // add this det element to module
      newModuleTransformer->GetDetElementStore()->Add(newDetElement->GetId(),
						      newDetElement);

      // In the Alice Alignment Framework misalignment objects store
      // global delta transformation
      // Get detection "intermediate" global transformation
      TGeoHMatrix newOldGlobalTransform = newModuleTransform * localTransform;
      // Get detection element global delta transformation: 
      // Tdelta = Tnew * Told.inverse
      TGeoHMatrix  deltaGlobalTransform
	= AliMUONGeometryBuilder::Multiply(newGlobalTransform, 
					   newOldGlobalTransform.Inverse());
	  
      // Create mis alignment matrix
      newGeometryTransformer
	->AddMisAlignDetElement(detElement->GetId(), deltaGlobalTransform);
    }
      
    if (verbose)
      AliInfo(Form("Added module transformer %i to the transformer", iMt));
    newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
  }
  return newGeometryTransformer;
}

//______________________________________________________________________
void AliMUONAlignment::SetAlignmentResolution(const TClonesArray* misAlignArray, Int_t rChId, Double_t rChResX, Double_t rChResY, Double_t rDeResX, Double_t rDeResY){
  /// Set alignment resolution to misalign objects to be stored in CDB
  Int_t chIdMin = (rChId<0)? 0 : rChId;
  Int_t chIdMax = (rChId<0)? 9 : rChId;
  Double_t chResX = rChResX;
  Double_t chResY = rChResY;
  Double_t deResX = rDeResX;
  Double_t deResY = rDeResY;

  TMatrixDSym mChCorrMatrix(6);
  mChCorrMatrix[0][0]=chResX*chResX;
  mChCorrMatrix[1][1]=chResY*chResY;
  //  mChCorrMatrix.Print();

  TMatrixDSym mDECorrMatrix(6);
  mDECorrMatrix[0][0]=deResX*deResX;
  mDECorrMatrix[1][1]=deResY*deResY;
  //  mDECorrMatrix.Print();

  AliAlignObjMatrix *alignMat = 0x0;

  for(Int_t chId=chIdMin; chId<=chIdMax; chId++) {
    TString chName1;
    TString chName2;
    if (chId<4){
      chName1 = Form("GM%d",chId);
      chName2 = Form("GM%d",chId);
    } else {
      chName1 = Form("GM%d",4+(chId-4)*2);
      chName2 = Form("GM%d",4+(chId-4)*2+1);
    }
    
    for (int i=0; i<misAlignArray->GetEntries(); i++) {
      alignMat = (AliAlignObjMatrix*)misAlignArray->At(i);
      TString volName(alignMat->GetSymName());
      if((volName.Contains(chName1)&&
	  ((volName.Last('/')==volName.Index(chName1)+chName1.Length())||
	   (volName.Length()==volName.Index(chName1)+chName1.Length())))||
	 (volName.Contains(chName2)&&
	  ((volName.Last('/')==volName.Index(chName2)+chName2.Length())||
	   (volName.Length()==volName.Index(chName2)+chName2.Length())))){
	volName.Remove(0,volName.Last('/')+1);
	if (volName.Contains("GM")) {
	  //	alignMat->Print("NULL");
	  alignMat->SetCorrMatrix(mChCorrMatrix);
	} else if (volName.Contains("DE")) {
	  //	alignMat->Print("NULL");
	  alignMat->SetCorrMatrix(mDECorrMatrix);
	}
      }
    }
  }
}
