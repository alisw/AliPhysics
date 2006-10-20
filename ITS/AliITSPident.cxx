//////////////////////////////////////////////////////////////////////////
// USER Class for PID in the ITS                                          //
//The PID is based on the likelihood of all the four ITS' layers,         //
//without using the truncated mean for the dE/dx. The response            //
//functions for each layer are convoluted Landau-Gaussian functions.      // 
// Origin: Elena Bruna bruna@to.infn.it,, Massimo Masera masera@to.infn.it//
//////////////////////////////////////////////////////////////////////////
#include "AliITSPident.h"

ClassImp(AliITSPident)
  //_______________________________________________________________________
AliITSPident::AliITSPident():
fMom(0),
fdEdx(0),
fPBayesp(0),
fPBayesk(0),
fPBayespi(0),
fPPriorip(0),
fPPriorik(0),
fPPrioripi(0),
fPPriorie(0),
fInvPt(0)
{
  // default constructor
  for (Int_t i=0;i<4;i++){
    fCondFunProLay[i]=0;
    fCondFunKLay[i]=0;
    fCondFunPiLay[i]=0;
  }
}
//_______________________________________________________________________
AliITSPident::~AliITSPident(){
  // destructor
}
//______________________________________________________________________
AliITSPident::AliITSPident(const AliITSPident &ob) :TObject(ob),
fMom(ob.fMom),
fdEdx(ob.fdEdx),
fPBayesp(ob.fPBayesp),
fPBayesk(ob.fPBayesk),
fPBayespi(ob.fPBayespi),
fPPriorip(ob.fPPriorip),
fPPriorik(ob.fPPriorik),
fPPrioripi(ob.fPPrioripi),
fPPriorie(ob.fPPriorie),
fInvPt(ob.fInvPt)  
{
  // Copy constructor
}

//______________________________________________________________________
AliITSPident& AliITSPident::operator=(const AliITSPident&  ob){
  // Assignment operator
  this->~AliITSPident();
  new(this) AliITSPident(ob);
  return *this;
}


//_______________________________________________________________________
AliITSPident::AliITSPident(Double_t mom,Double_t invPt,Double_t dEdx,AliITSSteerPid *sp,Float_t *Qlay,Float_t priorip,Float_t priorik,Float_t prioripi,Float_t priorie):
fMom(mom),
fdEdx(dEdx),
fPBayesp(0),
fPBayesk(0),
fPBayespi(0),
fPPriorip(priorip),
fPPriorik(priorik),
fPPrioripi(prioripi),
fPPriorie(priorie),
fInvPt(invPt){

  //test

  for(Int_t la=0;la<4;la++){//loop on layers
    Double_t parp[3];Double_t park[3];Double_t parpi[3];
    sp->GetParFitLayer(la,fMom,parp,park,parpi);

    Double_t range[6];
    range[0]=0.3*parp[1];
    range[1]=2.*parp[1];

    range[2]=0.3*park[1];
    range[3]=2.*park[1];
    
    range[4]=0.3*parpi[1];
    range[5]=2.*parpi[1];
      CookFunItsLay(la,0,parp,Qlay[la],fMom,range[0],range[1],"fPro");
      CookFunItsLay(la,1,park,Qlay[la],fMom,range[2],range[3],"fKao");
      CookFunItsLay(la,2,parpi,Qlay[la],fMom,range[4],range[5],"fPi");
    
  }

  Float_t prior[4];Double_t condFun[4][3];

  prior[0]=fPPriorip;
  prior[1]=fPPriorik;
  prior[2]=fPPrioripi;
  prior[3]=fPPriorie;
  for(Int_t la=0;la<4;la++){
    condFun[la][0]= fCondFunProLay[la];
    condFun[la][1]= fCondFunKLay[la]; 
    condFun[la][2]= fCondFunPiLay[la];
    
  }

  fPBayesp=CookCombinedBayes(condFun,prior,0);
  fPBayesk=CookCombinedBayes(condFun,prior,1); 
  fPBayespi=CookCombinedBayes(condFun,prior,2); 
 }


//_______________________________________________________________________
AliITSPident::AliITSPident(AliITStrackV2 *trackITS,AliITSSteerPid *sp,Float_t *Qlay,Float_t priorip,Float_t priorik,Float_t prioripi,Float_t priorie):
fMom(0),
fdEdx(0),
fPBayesp(0),
fPBayesk(0),
fPBayespi(0),
fPPriorip(priorip),
fPPriorik(priorik),
fPPrioripi(prioripi),
fPPriorie(priorie),
fInvPt(0)
{
  //
  Double_t xr;
  Double_t par[5];
  trackITS->GetExternalParameters(xr,par);
  if (par[4]!=0) {
    Float_t lamb=TMath::ATan(par[3]);
    fMom=1/(TMath::Abs(par[4])*TMath::Cos(lamb));
  }
 
  fInvPt=par[4];
 
  for(Int_t la=0;la<4;la++){//loop on layers
    Double_t parp[3];Double_t park[3];Double_t parpi[3];
    sp->GetParFitLayer(la,fMom,parp,park,parpi);
    Double_t range[8];
    range[0]=0.3*parp[1];
    range[1]=2.*parp[1];

    range[2]=0.3*park[1];
    range[3]=2.*park[1];
    
    range[4]=0.3*parpi[1];
    range[5]=2.*parpi[1];
      CookFunItsLay(la,0,parp,Qlay[la],fMom,range[0],range[1],"fPro");
      CookFunItsLay(la,1,park,Qlay[la],fMom,range[2],range[3],"fKao");
      CookFunItsLay(la,2,parpi,Qlay[la],fMom,range[4],range[5],"fPi");
  }

  Float_t prior[4];Double_t condFun[4][3];

  prior[0]=fPPriorip;
  prior[1]=fPPriorik;
  prior[2]=fPPrioripi;
  prior[3]=fPPriorie;
  for(Int_t la=0;la<4;la++){
    condFun[la][0]= fCondFunProLay[la];
    condFun[la][1]= fCondFunKLay[la]; 
    condFun[la][2]= fCondFunPiLay[la];
    
  }

  fPBayesp=CookCombinedBayes(condFun,prior,0);
  fPBayesk=CookCombinedBayes(condFun,prior,1); 
  fPBayespi=CookCombinedBayes(condFun,prior,2); 
  fdEdx=trackITS->GetdEdx();

}
//_______________________________________________________________________
void AliITSPident::PrintParameters() const{
 //print parameters
  cout<<"___________________________\n";
  cout<<"Track Local Momentum = "<<"  "<<fMom<<endl;
  cout<<"Track dEdx = "<<"  "<<fdEdx<<endl;
  cout<<"Inv Pt = "<<fInvPt<<endl;
}

//_______________________________________________________________________
Double_t AliITSPident::Langaufun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

//_______________________________________________________________________
Double_t AliITSPident::Langaufun2(Double_t *x, Double_t *par){
  //normalized langaufun
  return 1/par[4]*Langaufun(x,par);
}
//_______________________________________________________________________
Double_t AliITSPident::Langaufunnorm(Double_t *x, Double_t *par){
   //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  
  //par[2]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[2];
  xupp = x[0] + sc * par[2];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[2]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[2]);
  }

  return (step * sum * invsq2pi / par[2]);
}
//_______________________________________________________________________
Double_t AliITSPident::Gaus2(Double_t *x, Double_t *par){
  //normalized gaussian function
  return 1/(sqrt(2*TMath::Pi())*par[1])*TMath::Gaus(x[0],par[0],par[1]);
}
//_______________________________________________________________________
void AliITSPident::CookFunItsLay(Int_t lay,Int_t opt,Double_t *par,Double_t dedx,Double_t mom,Double_t rangei,Double_t rangef,TString comment){
  //it gives the response functions
 TF1 funLay(comment,Langaufunnorm,rangei,rangef,3);
 funLay.SetParameters(par);
 Double_t condFun=funLay.Eval(dedx);
  if(opt==0){
    fCondFunProLay[lay]=condFun;
    if(mom<0.4 && dedx<100)fCondFunProLay[lay]=0.00001;
    if(mom<0.4 &&dedx<50)fCondFunProLay[lay]=0.0000001;
  }
  if(opt==1){
    fCondFunKLay[lay]=condFun; 
    if(mom<0.25 && dedx<100)fCondFunKLay[lay]=0.00001;
    if(mom<0.4 &&dedx<30)fCondFunKLay[lay]=0.0000001;   
  }
  if(opt==2){
    fCondFunPiLay[lay]=condFun;
    if(mom<0.6 &&dedx<20)fCondFunPiLay[lay]=0.001;
  }

}
//_______________________________________________________________________
Float_t AliITSPident::CookCombinedBayes(Double_t condfun[][3],Float_t *prior,Int_t part)const {
  //Bayesian combined PID in the ITS
  Int_t test=0; 
  Float_t bayes;
  Float_t pprior[4]={0,0,0,0};
  for(Int_t j=0;j<4;j++)pprior[j]=prior[j];
  pprior[2]+=pprior[3];//prior for electrons summed to priors for pions
  for(Int_t i=0;i<4;i++){//layer
    if (condfun[i][0]>0 || condfun[i][1]>0 ||condfun[i][2]>0) test++; 
  }
  if(test==4){
    if ((pprior[0]!=0 || pprior[1]!=0 ||pprior[2]!=0)&&CookSum(condfun,pprior)!=0){

      bayes=pprior[part]*CookProd(condfun,part)*1/CookSum(condfun,pprior);
  
    }
    else bayes=-100;
  }
  else bayes=-100;
  return bayes;
}
//_______________________________________________________________________
Float_t AliITSPident::CookProd(Double_t condfun[][3],Int_t part)const{
  //
  Float_t p=1;
  for(Int_t lay=0;lay<4;lay++){
  
    p=p*condfun[lay][part];
  }

  return p;
}
//_______________________________________________________________________
Float_t AliITSPident::CookSum(Double_t condfun[][3],Float_t *prior)const{
  //
  Float_t sum=0;
  Float_t pprior[4]={0,0,0,0};
  for(Int_t j=0;j<4;j++)pprior[j]=prior[j];
  pprior[2]+=pprior[3];//prior for electrons summed to priors for pions
  for(Int_t i=0;i<3;i++){//sum over the particles--electrons excluded
    sum+=pprior[i]*CookProd(condfun,i);
  }
  return sum;

}
