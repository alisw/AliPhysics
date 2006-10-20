#ifndef ALIITSPIDENT_H
#define ALIITSPIDENT_H

/////////////////////////////////////////////////////////////////////////
//Class for PID in the ITS                                               //
//The PID is based on the likelihood of all the four ITS' layers,        //
//without using the truncated mean for the dE/dx. The response           //
//functions for each layer are convoluted Landau-Gaussian functions.     // 
// Origin: Elena Bruna bruna@to.infn.it, Massimo Masera masera@to.infn.it//
////////////////////////////////////////////////////////////////////////
#include "AliITStrackV2.h"
#include <TObject.h>
#include <Riostream.h>
#include <TF1.h>
#include <TTree.h>
#include "AliITSPidParItem.h"
#include "AliITSSteerPid.h"
class AliITSPident : public TObject{

 public:
  AliITSPident();
  AliITSPident(Double_t mom,Double_t invPt,Double_t dEdx,AliITSSteerPid *sp,Float_t *Qlay,Float_t priorip=0.066,Float_t priorik=0.103,Float_t prioripi=0.83,Float_t priorie=0.001);
  
  AliITSPident(AliITStrackV2 *trackITS,AliITSSteerPid *sp,Float_t *Qlay,Float_t priorip=0.066,Float_t priorik=0.103,Float_t prioripi=0.83,Float_t priorie=0.001);

  virtual ~AliITSPident();
  Float_t GetP() const {return fMom;}//local momentum (GeV/c)
  Double_t GetdEdx() const {return fdEdx;} //dEdx (ADC units)
  Double_t GetCondFunPro(Int_t lay) const {
    return fCondFunProLay[lay];
  }
  Double_t GetProdCondFunPro() const {
    Double_t rv=1.; for(Int_t i=0;i<4;i++)rv*=GetCondFunPro(i);
    return rv;
  }
  Double_t GetCondFunK(Int_t lay) const {
    return fCondFunKLay[lay];
  }
  Double_t GetProdCondFunK() const {
    Double_t rv=1.; for(Int_t i=0;i<4;i++)rv*=GetCondFunK(i);
    return rv;
  }
  Double_t GetCondFunPi(Int_t lay) const {
    return fCondFunPiLay[lay];
  }
  Double_t GetProdCondFunPi() const {
    Double_t rv=1.; for(Int_t i=0;i<4;i++)rv*=GetCondFunPi(i);
    return rv;
  }
  Float_t GetInvPt() const {return fInvPt;}
  void PrintParameters() const;
  Float_t GetPBayesp()const {return fPBayesp;}
  Float_t GetPBayesk()const {return fPBayesk;}
  Float_t GetPBayespi()const {return fPBayespi;}
  Float_t GetPPriorip() const {return fPPriorip;}
  Float_t GetPPriorik() const {return fPPriorik;}
  Float_t GetPPrioripi() const {return fPPrioripi;}
  Float_t GetPPriorie() const {return fPPriorie;}
  
  static Double_t Langaufun(Double_t *x, Double_t *par);
  static Double_t Langaufun2(Double_t *x, Double_t *par);
  static Double_t Langaufunnorm(Double_t *x, Double_t *par);
  static Double_t Gaus2(Double_t *x, Double_t *par);
  
 private:

  void CookFunItsLay(Int_t lay,Int_t opt,Double_t *par,Double_t dedx,Double_t mom,Double_t rangei,Double_t rangef,TString comment);
  void CookBayes(Double_t *condfun,Float_t *prior);
  Float_t CookCombinedBayes(Double_t condfun[][3],Float_t *prior,Int_t part)const;
  Float_t CookProd(Double_t condfun[][3],Int_t part) const;
  Float_t CookSum(Double_t condfun[][3],Float_t *prior) const;
  AliITSPident(const AliITSPident &ob); // copy constructor
  AliITSPident& operator=(const AliITSPident & ob); // ass. op.

  Float_t fMom;                   // Particle momentum
  Double_t fdEdx;                 // Particle dE/dx
  Double_t fCondFunProLay[4];     // one for each silicon layer
  Double_t fCondFunKLay[4];       // cond. prob. function kaons per layer
  Double_t fCondFunPiLay[4];      // cond. prob. function pions per layer
  Float_t fPBayesp;               // Bayes prob. 
  Float_t fPBayesk;               // Bayes prob. for kaons
  Float_t fPBayespi;              // Bayes prob. for pions 
  Float_t fPPriorip;              // Priori prob. 
  Float_t fPPriorik;              // Priori prob. for kaons
  Float_t fPPrioripi;             // Priori prob. for pions
  Float_t fPPriorie;              // Priori prob. for electrons
  Float_t fInvPt;                 // Particle Inverse Transverse momentum

  ClassDef(AliITSPident,1);
};
#endif
