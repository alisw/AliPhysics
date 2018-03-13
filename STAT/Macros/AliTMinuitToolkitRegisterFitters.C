/*
  .L $AliRoot_SRC/STAT/Macros/AliTMinuitToolkitRegisterFitters.C+

 */

#include "AliTMinuitToolkit.h"
#include "AliTMinuitToolkit.h"

Double_t GaussPol0(Double_t* x, Double_t *p){
  Double_t gauss=p[0]*TMath::Gaus(x[0],p[1],p[2],kFALSE);
  return p[3]+gaus;
}
Double_t GaussPol1(Double_t* x, Double_t *p){
  Double_t gauss=p[0]*TMath::Gaus(x[0],p[1],p[2],kFALSE);
  return p[3]+p[4]*(x[0]-p[1]) + gaus;
}


void RegisterGaussFitteres(){
  TF1 *fGaus = new TF1("fgaus","gaus");
  TMatrixD initPar(3,4);
  initPar(0,0)=0; initPar(0,1)=1;
  initPar(1,0)=0; initPar(1,1)=1;
  initPar(2,0)=1; initPar(2,1)=1;
  AliTMinuitToolkit * fitterGR = new AliTMinuitToolkit();
  fitterGR->SetVerbose(0x1); fitterGR->SetFitFunction(fGaus,kTRUE);
  fitterGR->SetLogLikelihoodFunction(likeGausCachy);
  fitterGR->SetInitialParam(&initPar);
  AliTMinuitToolkit::SetPredefinedFitter("gausR",fitterGR);

}