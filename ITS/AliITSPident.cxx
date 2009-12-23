//////////////////////////////////////////////////////////////////////////
// USER Class for PID in the ITS                                          //
//The PID is based on the likelihood of all the four ITS' layers,         //
//without using the truncated mean for the dE/dx. The response            //
//functions for each layer are convoluted Landau-Gaussian functions.      // 
// Origin: Elena Bruna bruna@to.infn.it,, Massimo Masera masera@to.infn.it//
//////////////////////////////////////////////////////////////////////////
#include "AliESDtrack.h"
#include "AliITSPidParams.h"
#include "AliPID.h"
#include "AliITSPident.h"

ClassImp(AliITSPident)
  //_______________________________________________________________________
AliITSPident::AliITSPident():
fPBayesp(0),
fPBayesk(0),
fPBayespi(0),
fPPriorip(0),
fPPriorik(0),
fPPrioripi(0),
fPPriorie(0)
{
  // default constructor
  for (Int_t i=0;i<8;i++){
    fCondFunProLay[i]=0;
    fCondFunKLay[i]=0;
    fCondFunPiLay[i]=0;
  }
}
//_______________________________________________________________________
AliITSPident::~AliITSPident(){
  // destructor
}
//_______________________________________________________________________
AliITSPident::AliITSPident(Double_t mom,AliITSPidParams *pars, Double_t *Qlay,Double_t priorip,Double_t priorik,Double_t prioripi,Double_t priorie):
fPBayesp(0),
fPBayesk(0),
fPBayespi(0),
fPPriorip(priorip),
fPPriorik(priorik),
fPPrioripi(prioripi),
fPPriorie(priorie)
{
  //
  CalculateResponses(mom,pars,Qlay);
}
//__________________________________________________________________________________________
AliITSPident::AliITSPident(AliESDtrack *track,AliITSPidParams *pars,Double_t priorip,Double_t priorik,Double_t prioripi,Double_t priorie):
fPBayesp(0),
fPBayesk(0),
fPBayespi(0),
fPPriorip(priorip),
fPPriorik(priorik),
fPPrioripi(prioripi),
fPPriorie(priorie)
{
  //
  Double_t mom=track->GetP();
  Double_t Qlay[4]={0.,0.,0.,0.};
  track->GetITSdEdxSamples(Qlay);
  CalculateResponses(mom,pars,Qlay);
} 
//_______________________________________________________________________
void AliITSPident::CalculateResponses(Double_t mom,AliITSPidParams *pars, Double_t *Qlay){
  // calculates conditional probabilities

  for (Int_t i=0;i<8;i++){
    fCondFunProLay[i]=-1;
    fCondFunKLay[i]=-1;
    fCondFunPiLay[i]=-1;
  }

  for(Int_t iLay=0; iLay<4; iLay++){//loop on layers (3=SDD inner 6=SSD outer)
    Double_t dedx=Qlay[iLay];
    if(dedx>0){
      fCondFunProLay[iLay]=pars->GetLandauGausNorm(dedx,AliPID::kProton,mom,iLay+3);
      if(mom<0.4 && dedx<100)fCondFunProLay[iLay]=0.00001;
      if(mom<0.4 &&dedx<50)fCondFunProLay[iLay]=0.0000001;
      fCondFunKLay[iLay]=pars->GetLandauGausNorm(dedx,AliPID::kKaon,mom,iLay+3);
      if(mom<0.25 && dedx<100)fCondFunKLay[iLay]=0.00001;
      if(mom<0.4 &&dedx<30)fCondFunKLay[iLay]=0.0000001;   
      fCondFunPiLay[iLay]=pars->GetLandauGausNorm(dedx,AliPID::kPion,mom,iLay+3);
      if(mom<0.6 &&dedx<20)fCondFunPiLay[iLay]=0.001;
    }
  }
  Double_t prior[4];
  Double_t condFun[8][3];

  prior[0]=fPPriorip;
  prior[1]=fPPriorik;
  prior[2]=fPPrioripi;
  prior[3]=fPPriorie;
  for(Int_t iLay=0;iLay<8;iLay++){
    condFun[iLay][0]= fCondFunProLay[iLay];
    condFun[iLay][1]= fCondFunKLay[iLay]; 
    condFun[iLay][2]= fCondFunPiLay[iLay];
  }

  fPBayesp=CookCombinedBayes(condFun,prior,0);
  fPBayesk=CookCombinedBayes(condFun,prior,1); 
  fPBayespi=CookCombinedBayes(condFun,prior,2); 
}
//_______________________________________________________________________
Double_t AliITSPident::GetProdCondFunPro() const {
  //Product of conditional probability functions for protons
    Double_t rv=1.;
    for(Int_t i=0;i<8;i++){
      Double_t fun=GetCondFunPro(i);
      if(fun>=0)rv*=fun;
    }
    return rv;
}//_______________________________________________________________________
Double_t AliITSPident::GetProdCondFunK() const {
  //Product of conditional probability functions for kaons
    Double_t rv=1.; 
    for(Int_t i=0;i<8;i++){
      Double_t fun=GetCondFunK(i);
      if(fun>=0)rv*=fun;
    }
    return rv;
}
//_______________________________________________________________________
Double_t AliITSPident::GetProdCondFunPi() const {
  //Product of conditional probability functions for pions
  Double_t rv=1.; 
  for(Int_t i=0;i<8;i++){
    Double_t fun=GetCondFunPi(i);
    if(fun>=0)rv*=fun;
  }
  return rv;
}
//_______________________________________________________________________
Double_t AliITSPident::CookCombinedBayes(Double_t condfun[][3],Double_t *prior,Int_t part)const {
  //Bayesian combined PID in the ITS
  Int_t test=0; 
  Double_t bayes;
  Double_t pprior[4]={0,0,0,0};
  for(Int_t j=0;j<4;j++)pprior[j]=prior[j];
  pprior[2]+=pprior[3];//prior for electrons summed to priors for pions
  for(Int_t i=0;i<8;i++){//layer
    if (condfun[i][0]>0 || condfun[i][1]>0 ||condfun[i][2]>0) test++; 
  }

  if(test>0){
    Double_t sum=CookSum(condfun,pprior);
    if ((pprior[0]!=0 || pprior[1]!=0 ||pprior[2]!=0)&&sum!=0.){

      bayes=pprior[part]*CookProd(condfun,part)*1/sum;
  
    }
    else bayes=-100;
  }
  else bayes=-100;
  return bayes;
}
//_______________________________________________________________________
Double_t AliITSPident::CookProd(Double_t condfun[][3],Int_t part)const{
  //
  Double_t p=1;
  for(Int_t lay=0;lay<8;lay++){
    if(condfun[lay][part]>0.) p=p*condfun[lay][part];
  }
  return p;
}
//_______________________________________________________________________
Double_t AliITSPident::CookSum(Double_t condfun[][3],Double_t *prior)const{
  //
  Double_t sum=0.;
  Double_t pprior[4]={0,0,0,0};
  for(Int_t j=0;j<4;j++) pprior[j]=prior[j];
  pprior[2]+=pprior[3];//prior for electrons summed to priors for pions
  for(Int_t i=0;i<3;i++){//sum over the particles--electrons excluded
    sum+=pprior[i]*CookProd(condfun,i);
  }
  return sum;

}
