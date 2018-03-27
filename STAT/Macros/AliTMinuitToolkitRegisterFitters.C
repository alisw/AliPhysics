/*
  .L $AliRoot_SRC/STAT/Macros/AliTMinuitToolkitRegisterFitters.C+
RegisterGaussFitters();

 */

#include "AliTMinuitToolkit.h"
#include "TF1.h"
#include "TMath.h"
#include "TMatrixD.h"

Double_t GaussPol0(Double_t* x, Double_t *p){
  Double_t gauss=p[0]*TMath::Gaus(x[0],p[1],p[2],kFALSE);
  return p[3]+gauss;
}
Double_t GaussPol1(Double_t* x, Double_t *p){
  Double_t gauss=p[0]*TMath::Gaus(x[0],p[1],p[2],kFALSE);
  return p[3]+p[4]*(x[0]-p[1]) + gauss;
}


void RegisterGaussFitters(){
  TF1 *likeGausCachy = new TF1("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike,-10,10,2);
   likeGausCachy->SetParameters(0.8,1);

  TF1 *fGaus0 = new TF1("fgaus0",GaussPol0,-10,10,4);
  TF1 *fGaus1 = new TF1("fgaus1",GaussPol1,-10,10,5);
  TMatrixD initPar0(4,4);
  initPar0(0,0)=1; initPar0(0,1)=1;
  initPar0(1,0)=0; initPar0(1,1)=1;
  initPar0(2,0)=1; initPar0(2,1)=1;
  initPar0(3,0)=1; initPar0(3,1)=1;
  TMatrixD initPar1(5,4);
  initPar1(0,0)=1; initPar1(0,1)=1;
  initPar1(1,0)=0; initPar1(1,1)=1;
  initPar1(2,0)=1; initPar1(2,1)=1;
  initPar1(3,0)=1; initPar1(3,1)=1;
  initPar1(4,0)=1; initPar1(4,1)=1;

  AliTMinuitToolkit * fitterGausPol0 = new AliTMinuitToolkit();
  fitterGausPol0->SetVerbose(0x1); fitterGausPol0->SetFitFunction(fGaus0,kTRUE);
  fitterGausPol0->SetLogLikelihoodFunction(likeGausCachy);
  fitterGausPol0->SetInitialParam(&initPar0);
  AliTMinuitToolkit::SetPredefinedFitter("gausPol0",fitterGausPol0);
  AliTMinuitToolkit * fitterGausPol1 = new AliTMinuitToolkit();
  fitterGausPol1->SetVerbose(0x1); fitterGausPol1->SetFitFunction(fGaus1,kTRUE);
  fitterGausPol1->SetLogLikelihoodFunction(likeGausCachy);
  fitterGausPol1->SetInitialParam(&initPar1);
  AliTMinuitToolkit::SetPredefinedFitter("gausPol1",fitterGausPol1);


}