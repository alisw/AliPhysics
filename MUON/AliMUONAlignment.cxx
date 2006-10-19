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
/// \class AliMUONAlignmentV5
/// Alignment class fro the ALICE DiMuon spectrometer 
///
/// MUON specific alignment class which interface to AliMillepede.   
/// For each track ProcessTrack calculates the local and global derivatives
/// at each hit and fill the corresponding local equations. Provide methods for
/// fixing or constraining detection elements for best results. 
///
/// \author Bruce Becker, Javier Castillo
//-----------------------------------------------------------------------------

#include "AliMUONAlignment.h"
#include "AliMUONTrack.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrackParam.h"
#include "AliMUONHitForRec.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMUONConstants.h"
#include "AliMillepede.h"

#include "AliMpExMap.h"

#include "AliLog.h"

#include "TSystem.h"

ClassImp(AliMUONAlignment)
  Int_t AliMUONAlignment::fgNDetElem = 4*2+4*2+18*2+26*2+26*2;
  Int_t AliMUONAlignment::fgNDetElemCh[10] = {4,4,4,4,18,18,26,26,26,26};
  Int_t AliMUONAlignment::fgSNDetElemCh[10] = {4,8,12,16,34,52,78,104,130,156};
  Int_t AliMUONAlignment::fgNParCh = 3;
  Int_t AliMUONAlignment::fgNCh = 10;
  Int_t AliMUONAlignment::fgNSt = 5;

AliMUONAlignment::AliMUONAlignment() 
  : TObject(),
    fBFieldOn(kTRUE),
    fStartFac(16.), 
    fResCutInitial(100.), 
    fResCut(100.),
    fMillepede(0),
    fTrackParamAtHit(0),
    fHitForRecAtHit(0),
    fTrack(0),
    fRecHit(0),
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

  fDoF[0] = kTRUE;  fDoF[1] = kTRUE;  fDoF[2] = kTRUE;
  fAllowVar[0] = 0.05;  fAllowVar[1] = 0.05;  fAllowVar[2] = 0.001;
  
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

  Bool_t bStOnOff[5] = {kFALSE,kFALSE,kTRUE,kTRUE,kTRUE};

  AllowVariations(bStOnOff);

  // Fix parameters or add constraints here
  for (Int_t iSt=0; iSt<5; iSt++)
    if (!bStOnOff[iSt]) FixStation(iSt+1);

  Bool_t bVarXYT[3] = {kFALSE,kTRUE,kFALSE};
  Bool_t bDetTLBR[4] = {kFALSE,kTRUE,kFALSE,kTRUE};
  ResetConstraints();
  AddConstraints(bStOnOff,bVarXYT,bDetTLBR);
  bVarXYT[0] = kFALSE; bVarXYT[1] = kFALSE; bVarXYT[2] = kTRUE;
  bDetTLBR[0] = kFALSE; bDetTLBR[1] = kTRUE; bDetTLBR[2] = kFALSE; bDetTLBR[3] = kFALSE;
  AddConstraints(bStOnOff,bVarXYT,bDetTLBR);
  bVarXYT[0] = kFALSE; bVarXYT[1] = kFALSE; bVarXYT[2] = kTRUE;
  AddConstraints(bStOnOff,bVarXYT);
  
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
  }
}

void AliMUONAlignment::SetNonLinear(Bool_t *lStOnOff,Bool_t *lVarXYT){
  /// Set non linear parameter flag selected stations and degrees of freedom
  for (Int_t i = 0; i < fgNDetElem; i++){    
    Int_t iCh=0;
    for (iCh=1; iCh<=fgNCh; iCh++){
      if (i<fgSNDetElemCh[iCh-1]) break;
    }
    Int_t iSt = lStOnOff[(iCh-1)/2] ? (iCh+1)/2 : 0; 
    if (iSt){
      if (lVarXYT[0]) { // X constraint
	SetNonLinear(i*fgNParCh+0);
      }
      if (lVarXYT[1]) { // Y constraint
	SetNonLinear(i*fgNParCh+1);
      }
      if (lVarXYT[2]) { // T constraint
	SetNonLinear(i*fgNParCh+2);
      }
    }
  }
}

void AliMUONAlignment::AddConstraints(Bool_t *lStOnOff,Bool_t *lVarXYT){
  /// Add constraint equations for selected stations and degrees of freedom 
  for (Int_t i = 0; i < fgNDetElem; i++){    
    Int_t iCh=0;
    for (iCh=1; iCh<=fgNCh; iCh++){
      if (i<fgSNDetElemCh[iCh-1]) break;
    }
    Int_t iSt = lStOnOff[(iCh-1)/2] ? (iCh+1)/2 : 0; 
    if (iSt){
      if (lVarXYT[0]) { // X constraint
	fConstraintX[i*fgNParCh+0]=1.0;
      }
      if (lVarXYT[1]) { // Y constraint
	fConstraintY[i*fgNParCh+1]=1.0;
      }
      if (lVarXYT[2]) { // T constraint
	fConstraintP[i*fgNParCh+2]=1.0;
      }
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
}

void AliMUONAlignment::AddConstraints(Bool_t *lStOnOff,Bool_t *lVarXYT, Bool_t *lDetTLBR){
  /// Add constraint equations for selected stations, degrees of freedom detector half 
  for (Int_t i = 0; i < fgNDetElem; i++){    
    Int_t iCh=0;
    for (iCh=1; iCh<=fgNCh; iCh++){
      if (i<fgSNDetElemCh[iCh-1]) break;
    }
    Int_t iSt = lStOnOff[(iCh-1)/2] ? (iCh+1)/2 : 0; 
    if (iSt){
      if (lVarXYT[0]) { // X constraint
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintXT,0); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintXL,0); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintXB,0); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintXR,0); // Right half
      }
      if (lVarXYT[1]) { // X constraint
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintYT,1); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintYL,1); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintYB,1); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintYR,1); // Right half
      }
      if (lVarXYT[2]) { // X constraint
	if (lDetTLBR[0]) ConstrainT(i,iCh,fConstraintPT,2); // Top half
	if (lDetTLBR[1]) ConstrainL(i,iCh,fConstraintPL,2); // Left half
	if (lDetTLBR[2]) ConstrainB(i,iCh,fConstraintPB,2); // Bottom half
	if (lDetTLBR[3]) ConstrainR(i,iCh,fConstraintPR,2); // Right half
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
}

void AliMUONAlignment::ConstrainT(Int_t lDetElem, Int_t lCh, Double_t *lConstraintT, Int_t iVar){
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

void AliMUONAlignment::ConstrainL(Int_t lDetElem, Int_t lCh, Double_t *lConstraintL, Int_t iVar){
  /// Set constrain equation for left half of spectrometer
  Int_t lDetElemNumber = (lCh==1) ? lDetElem : lDetElem-fgSNDetElemCh[lCh-2];
  if (lCh>=1 && lCh<=4){
    if (lDetElemNumber==1 || lDetElemNumber==2){ // From track crossings
      lConstraintL[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
  if (lCh>=5 && lCh<=6){
    if (lDetElemNumber>=5&&lDetElemNumber<=13){
      lConstraintL[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
  if (lCh>=7 && lCh<=10){
    if (lDetElemNumber>=7&&lDetElemNumber<=19){
      lConstraintL[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
}

void AliMUONAlignment::ConstrainB(Int_t lDetElem, Int_t lCh, Double_t *lConstraintB, Int_t iVar){
  /// Set constrain equation for bottom half of spectrometer
  Int_t lDetElemNumber = (lCh==1) ? lDetElem : lDetElem-fgSNDetElemCh[lCh-2];
  if (lCh>=1 && lCh<=4){
    if (lDetElemNumber==2 && lDetElemNumber==3){ // From track crossings
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

void AliMUONAlignment::ConstrainR(Int_t lDetElem, Int_t lCh, Double_t *lConstraintR, Int_t iVar){
  /// Set constrain equation for right half of spectrometer
  Int_t lDetElemNumber = (lCh==1) ? lDetElem : lDetElem-fgSNDetElemCh[lCh-2];
  if (lCh>=1 && lCh<=4){
    if (lDetElemNumber==0 && lDetElemNumber==3){ // From track crossings
      lConstraintR[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
  if (lCh>=5 && lCh<=6){
    if ((lDetElemNumber>=0&&lDetElemNumber<=4) || 
	(lDetElemNumber>=14&&lDetElemNumber<=17)){
      lConstraintR[lDetElem*fgNParCh+iVar]=1.0;
    }
  }
  if (lCh>=7 && lCh<=10){
    if ((lDetElemNumber>=0&&lDetElemNumber<=6) || 
	(lDetElemNumber>=20&&lDetElemNumber<=25)){
      lConstraintR[lDetElem*fgNParCh+iVar]=1.0;
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
      fConstraintXL[i*fgNParCh+0]=0.0;
      fConstraintXL[i*fgNParCh+1]=0.0;
      fConstraintXL[i*fgNParCh+2]=0.0;
      fConstraintYL[i*fgNParCh+0]=0.0;
      fConstraintYL[i*fgNParCh+1]=0.0;
      fConstraintYL[i*fgNParCh+2]=0.0;
      fConstraintPL[i*fgNParCh+0]=0.0;
      fConstraintPL[i*fgNParCh+1]=0.0;
      fConstraintPL[i*fgNParCh+2]=0.0;
      fConstraintXB[i*fgNParCh+0]=0.0;
      fConstraintXB[i*fgNParCh+1]=0.0;
      fConstraintXB[i*fgNParCh+2]=0.0;
      fConstraintYB[i*fgNParCh+0]=0.0;
      fConstraintYB[i*fgNParCh+1]=0.0;
      fConstraintYB[i*fgNParCh+2]=0.0;
      fConstraintPB[i*fgNParCh+0]=0.0;
      fConstraintPB[i*fgNParCh+1]=0.0;
      fConstraintPB[i*fgNParCh+2]=0.0;
      fConstraintXR[i*fgNParCh+0]=0.0;
      fConstraintXR[i*fgNParCh+1]=0.0;
      fConstraintXR[i*fgNParCh+2]=0.0;
      fConstraintYR[i*fgNParCh+0]=0.0;
      fConstraintYR[i*fgNParCh+1]=0.0;
      fConstraintYR[i*fgNParCh+2]=0.0;
      fConstraintPR[i*fgNParCh+0]=0.0;
      fConstraintPR[i*fgNParCh+1]=0.0;
      fConstraintPR[i*fgNParCh+2]=0.0;
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
  if (value==0) AliInfo(Form("Parameter %i Fixed", iPar));
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

void AliMUONAlignment::AllowVariations(Bool_t *bStOnOff) {
  /// Set allowed variation for selected stations based on fDoF and fAllowVar
  for (Int_t iSt=1; iSt<=5; iSt++) {
    if (bStOnOff[iSt-1]) {
      Int_t iDetElemFirst = (iSt>1) ? fgSNDetElemCh[2*(iSt-1)-1] : 0; 
      Int_t iDetElemLast = fgSNDetElemCh[2*(iSt)-1]; 
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
  /// Define local equation for current track and hit in x coor. measurement
  // set local derivatives
  SetLocalDerivative(0, fCosPhi);
  SetLocalDerivative(1, fCosPhi * (fTrackPos[2] - fTrackPos0[2]));
  SetLocalDerivative(2, fSinPhi);
  SetLocalDerivative(3, fSinPhi * (fTrackPos[2] - fTrackPos0[2]));

  // set global derivatives
  SetGlobalDerivative(fDetElemNumber*fgNParCh+0, -1.);
  SetGlobalDerivative(fDetElemNumber*fgNParCh+1,  0.);
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

  fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, fMeas[0], fSigma[0]);
}

void AliMUONAlignment::LocalEquationY() {
  /// Define local equation for current track and hit in y coor. measurement
  // set local derivatives
  SetLocalDerivative(0,-fSinPhi);
  SetLocalDerivative(1,-fSinPhi * (fTrackPos[2] - fTrackPos0[2]));
  SetLocalDerivative(2, fCosPhi);
  SetLocalDerivative(3, fCosPhi * (fTrackPos[2] - fTrackPos0[2]));

  // set global derivatives
  SetGlobalDerivative(fDetElemNumber*fgNParCh+0,  0.);
  SetGlobalDerivative(fDetElemNumber*fgNParCh+1, -1.);
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

  fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, fMeas[1], fSigma[1]);
}

void AliMUONAlignment::FillRecPointData() {
  /// Get information of current hit
  fClustPos[0] = fRecHit->GetNonBendingCoor();
  fClustPos[1] = fRecHit->GetBendingCoor();
  fClustPos[2] = fRecHit->GetZ();
  fTransform->Global2Local(fDetElemId,fClustPos[0],fClustPos[1],fClustPos[2],
			    fClustPosLoc[0],fClustPosLoc[1],fClustPosLoc[2]);
}

void AliMUONAlignment::FillTrackParamData() {
  /// Get information of current track at current hit
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
  fDetElemId = fRecHit->GetDetElemId();
  fDetElemNumber = fDetElemId%100;
  for (int iCh=0; iCh<fDetElemId/100-1; iCh++){
    fDetElemNumber += fgNDetElemCh[iCh];
  }
  fTransform->Local2Global(fDetElemId,lDetElemLocX,lDetElemLocY,lDetElemLocZ,
			   fDetElemPos[0],fDetElemPos[1],fDetElemPos[2]);
  if (fDetElemId/100 < 5){    
    fSigma[0] = 3.0e-2;
    fSigma[1] = 3.0e-2;
  }
  else {
    fSigma[0] = 1.0e-1;
    fSigma[1] = 1.0e-2;    
  }  
}

void AliMUONAlignment::ProcessTrack(AliMUONTrack * track) {
  /// Process track; Loop over hits and set local equations
  fTrack = track;
  // get tclones arrays.
  fTrackParamAtHit = fTrack->GetTrackParamAtHit();
  fHitForRecAtHit = fTrack->GetHitForRecAtHit();
  
  // get size of arrays
  Int_t nTrackParam = fTrackParamAtHit->GetEntries();
  Int_t nHitForRec = fHitForRecAtHit->GetEntries();
  AliInfo(Form("Number of track param entries : %i ", nTrackParam));
  AliInfo(Form("Number of hit for rec entries : %i ", nHitForRec));

  for(Int_t iHit=0; iHit<nHitForRec; iHit++) {
    fRecHit = (AliMUONHitForRec *) fHitForRecAtHit->At(iHit);
    fTrackParam = (AliMUONTrackParam *) fTrack->GetTrackParamAtHit()->At(iHit);
    if (!fRecHit || !fTrackParam) continue;
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

  for(Int_t iHit=0; iHit<nHitForRec; iHit++) {
    // and get new pointers
    fRecHit = (AliMUONHitForRec *) fHitForRecAtHit->At(iHit);
    fTrackParam = (AliMUONTrackParam *) fTrack->GetTrackParamAtHit()->At(iHit);
    if (!fRecHit || !fTrackParam) continue;
    // fill local variables for this position --> one measurement
    FillDetElemData();        
    FillRecPointData();
    FillTrackParamData();
//     if (fDetElemId<500) continue;
    AliDebug(1,Form("cluster: %i", iHit));
    AliDebug(1,Form("x: %f\t y: %f\t z: %f\t DetElemID: %i\t ", fClustPos[0], fClustPos[1], fClustPos[2], fDetElemId));
    AliDebug(1,Form("fDetElemPos[0]: %f\t fDetElemPos[1]: %f\t fDetElemPos[2]: %f\t DetElemID: %i\t ", fDetElemPos[0],fDetElemPos[1],fDetElemPos[2], fDetElemId));

    AliDebug(1,Form("Track Parameter: %i", iHit));
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
TGeoCombiTrans AliMUONAlignment::ReAlign(const TGeoCombiTrans & transform, double *detElemMisAlignment) const
{
  /// Realign given transformation by given misalignment and return the misaligned transformation
  
  Double_t cartMisAlig[3] = {0,0,0};
  Double_t angMisAlig[3] = {0,0,0};
  const Double_t *trans = transform.GetTranslation();
  TGeoRotation *rot;
  // check if the rotation we obtain is not NULL
  if (transform.GetRotation()) {    
    rot = transform.GetRotation();
  }
  else {    
    rot = new TGeoRotation("rot");
  }			// default constructor.

  cartMisAlig[0] = -detElemMisAlignment[0];
  cartMisAlig[1] = -detElemMisAlignment[1];
  angMisAlig[2] = -detElemMisAlignment[2]*180./TMath::Pi();

  TGeoTranslation newTrans(cartMisAlig[0] + trans[0], cartMisAlig[1] + trans[1], cartMisAlig[2] + trans[2]);
  
  rot->RotateX(angMisAlig[0]);
  rot->RotateY(angMisAlig[1]);
  rot->RotateZ(angMisAlig[2]);

  return TGeoCombiTrans(newTrans, *rot);
}

//______________________________________________________________________
AliMUONGeometryTransformer *
AliMUONAlignment::ReAlign(const AliMUONGeometryTransformer * transformer,
			    double *misAlignments, Bool_t verbose)
			    
{
  /////////////////////////////////////////////////////////////////////
  //   Takes the internal geometry module transformers, copies them
  // and gets the Detection Elements from them.
  // Takes misalignment parameters and applies these
  // to the local transform of the Detection Element
  // Obtains the global transform by multiplying the module transformer
  // transformation with the local transformation 
  // Applies the global transform to a new detection element
  // Adds the new detection element to a new module transformer
  // Adds the new module transformer to a new geometry transformer
  // Returns the new geometry transformer

  Double_t lDetElemMisAlignment[3] = {0.,0.,0.};

  AliMUONGeometryTransformer *newGeometryTransformer =
    new AliMUONGeometryTransformer(kTRUE);
  for (Int_t iMt = 0; iMt < transformer->GetNofModuleTransformers(); iMt++)
    {				// module transformers
      
      const AliMUONGeometryModuleTransformer *kModuleTransformer =
	transformer->GetModuleTransformer(iMt, true);
      
      AliMUONGeometryModuleTransformer *newModuleTransformer =
	new AliMUONGeometryModuleTransformer(iMt);
      newGeometryTransformer->AddModuleTransformer(newModuleTransformer);

      TGeoCombiTrans moduleTransform =
	TGeoCombiTrans(*kModuleTransformer->GetTransformation());
      TGeoCombiTrans *newModuleTransform = new TGeoCombiTrans(moduleTransform);	
              // same module transform as the previous one 
	      // no mis align object created
      newModuleTransformer->SetTransformation(moduleTransform);

      AliMpExMap *detElements = kModuleTransformer->GetDetElementStore();

      if (verbose)
	AliInfo(Form
		("%i DEs in old GeometryStore  %i",
		 detElements->GetSize(), iMt));
      Int_t iBase = !iMt ? 0 : fgSNDetElemCh[iMt-1];
      for (Int_t iDe = 0; iDe < detElements->GetSize(); iDe++)
	{			// detection elements.
	  AliMUONGeometryDetElement *detElement =
	    (AliMUONGeometryDetElement *) detElements->GetObject(iDe);
	  if (!detElement)
	    AliFatal("Detection element not found.");

	  /// make a new detection element
	  AliMUONGeometryDetElement *newDetElement =
	    new AliMUONGeometryDetElement(detElement->GetId(),
					  detElement->GetVolumePath());
	  for (int i=0; i<fgNParCh; i++) {
	    AliInfo(Form("iMt %i, iBase %i, iDe %i, iPar %i",iMt, iBase, iDe, (iBase+iDe)*fgNParCh+i));
	    lDetElemMisAlignment[i] = 
	      (iMt<fgNCh) ? misAlignments[(iBase+iDe)*fgNParCh+i] : 0.;  	      
	  }
	  // local transformation of this detection element.
          TGeoCombiTrans localTransform
	    = TGeoCombiTrans(*detElement->GetLocalTransformation());
	  TGeoCombiTrans newLocalTransform = ReAlign(localTransform,lDetElemMisAlignment);
          newDetElement->SetLocalTransformation(newLocalTransform);	  

	  // global transformation
	  TGeoHMatrix newGlobalTransform =
	    AliMUONGeometryBuilder::Multiply(*newModuleTransform,
	                                      newLocalTransform);
	  newDetElement->SetGlobalTransformation(newGlobalTransform);
	  
	  // add this det element to module
	  newModuleTransformer->GetDetElementStore()->Add(newDetElement->GetId(),
							  newDetElement);
          // Get delta transformation: 
	  // Tdelta = Tnew * Told.inverse
	  TGeoHMatrix  deltaTransform
	    = AliMUONGeometryBuilder::Multiply(
	        newGlobalTransform, 
		detElement->GetGlobalTransformation()->Inverse());
	  
	  // Create mis alignment matrix
	  newGeometryTransformer
	    ->AddMisAlignDetElement(detElement->GetId(), deltaTransform);
	}
      if (verbose)
	AliInfo(Form("Added module transformer %i to the transformer", iMt));
      newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
    }
  return newGeometryTransformer;
}

