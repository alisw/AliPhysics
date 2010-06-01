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
// AliTPCExBEffective class                                                   //
// Correct for the rest of ExB effect which are not covered by physical models
//
// Motivation:
//   ExB correction: 
//      dr    =  c0* integral(Er/Ez) + c1* integral(Erphi/Ez)
//      drphi = -c1* integral(Er/Ez) + c0* integral(Erphi/Ez)
//   Where:      
//   wt = Bz*(k*vdrift/E)           ~ 0.3 at B=0.5 T 
//   c0 = 1/(1+T2*T2*wt*wt) 
//   c1 = T1*wt/(1+T1*T1*wt*wt)
//   
// Residual integral(Er/Ez,Erphi/Ez) obtained comparing the B field 0 and B field +-0.5 T setting
// minimizing track matching residuals 
// delta(Er/Ez) ~ sum[ poln(r) * polm(z) * cos(n,phi)] 
//  
////////////////////////////////////////////////////////////////////////////
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"

#include "TMath.h"
#include "AliTPCROC.h"
#include "AliTPCExBEffective.h"
ClassImp(AliTPCExBEffective)

AliTPCExBEffective::AliTPCExBEffective()
  : AliTPCCorrection("ExB_effective","ExB effective"),
    fC0(1.),fC1(0.), 
    fPolynomA(0),
    fPolynomC(0),
    fPolynomValA(0),
    fPolynomValC(0)
{
  //
  // default constructor
  //
}

AliTPCExBEffective::~AliTPCExBEffective() {
  //
  // default destructor
  //
}



void AliTPCExBEffective::Init() {
  //
  // Initialization funtion
  //
  
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);


}

void AliTPCExBEffective::Update(const TTimeStamp &/*timeStamp*/) {
  //
  // Update function 
  //
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);


}



void AliTPCExBEffective::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction due conical shape
  //   
  if (!fPolynomA) return;
  AliTPCROC * calROC = AliTPCROC::Instance();
  const Double_t kRTPC0  =calROC->GetPadRowRadii(0,0);
  const Double_t kRTPC1  =calROC->GetPadRowRadii(36,calROC->GetNRows(36)-1);
  Float_t rmiddle=(kRTPC0+kRTPC1)/2.;
  //
  Double_t phi      = TMath::ATan2(x[1],x[0]);
  Double_t r        = TMath::Sqrt(x[1]*x[1]+x[0]*x[0]);
  Double_t driftN   = 1.-TMath::Abs(x[2])/calROC->GetZLength(0);  // drift from 0 to 1
  Double_t localxN  = 2*(r-rmiddle)/(kRTPC1-kRTPC0);         // normalize local x position
  //
  Double_t erez = 0;
  Double_t erphiez = 0;
  if (roc%36<18)  erez= GetSum(*fPolynomA, *fPolynomValA, localxN, driftN, phi,0);
  if (roc%36>=18) erez= GetSum(*fPolynomC, *fPolynomValC, localxN, driftN, phi,0);
  if (roc%36<18)  erphiez= GetSum(*fPolynomA, *fPolynomValA, localxN, driftN, phi,1);
  if (roc%36>=18) erphiez= GetSum(*fPolynomC, *fPolynomValC, localxN, driftN, phi,1);

  Double_t dr    =   fC0 * erez + fC1 * erphiez;
  Double_t drphi =  -fC1 * erez + fC0 * erphiez;
  //
  dx[0]= TMath::Cos(phi)*dr-TMath::Sin(phi)*drphi;
  dx[1]= TMath::Sin(phi)*dr+TMath::Cos(phi)*drphi;
  dx[2]= 0;

}



Double_t AliTPCExBEffective::GetSum(const TMatrixD& mpol, const TMatrixD&mcoef, Double_t r, Double_t drift, Double_t phi, Int_t coord) const {
  //
  //
  //
  Int_t npols=mpol.GetNrows();
  Double_t sum=0;
  for (Int_t ipol=0;ipol<npols; ipol++){
    Double_t pR = 1, pD=1, pPhi=1;
    Int_t icoord   = TMath::Nint(mpol(ipol,0));
    if (icoord!=coord) continue;
    Int_t npolR    = TMath::Nint(mpol(ipol,1));
    Int_t npolD    = TMath::Nint(mpol(ipol,2));
    Int_t npolPhi  = TMath::Nint(mpol(ipol,3));
    Double_t coef=mcoef(ipol,0);
    //
    for (Int_t ipolR=1; ipolR<=npolR; ipolR++) pR*=r;           // use simple polynoms
    for (Int_t ipolD=1; ipolD<=npolD; ipolD++) pD*=drift;       // use simple polynoms
    pPhi=TMath::Cos(npolPhi*phi);
    sum+= pR*pD*pPhi*coef;
  }
  return sum;
}
 

void AliTPCExBEffective::SetPolynoms(const TMatrixD *polA,const TMatrixD *polC){
  //
  // Set correction polynom - coefficients
  //
  fPolynomA = new TMatrixD(*polA);
  fPolynomC = new TMatrixD(*polC);
}

void AliTPCExBEffective::SetCoeficients(const TMatrixD *valA,const TMatrixD *valC){
  //
  // Set correction polynom - coefficients
  //
  fPolynomValA = new TMatrixD(*valA);
  fPolynomValC = new TMatrixD(*valC);
}

