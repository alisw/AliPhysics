#ifndef ALIITSPIDENT_H
#define ALIITSPIDENT_H

/////////////////////////////////////////////////////////////////////////
//Class for PID in the ITS                                               //
//The PID is based on the likelihood of all the four ITS' layers,        //
//without using the truncated mean for the dE/dx. The response           //
//functions for each layer are convoluted Landau-Gaussian functions.     // 
// Origin: Elena Bruna bruna@to.infn.it, Massimo Masera masera@to.infn.it//
////////////////////////////////////////////////////////////////////////

#include <TObject.h>
class AliESDtrack;
class TF1;
class AliITSPidParams;

class AliITSPident : public TObject{

 public:
  AliITSPident();
  AliITSPident(Double_t mom,AliITSPidParams *pars,Double_t *Qlay,Double_t priorip=0.066,Double_t priorik=0.103,Double_t prioripi=0.83,Double_t priorie=0.001);
  
  AliITSPident(AliESDtrack *track,AliITSPidParams *pars,Double_t priorip=0.066,Double_t priorik=0.103,Double_t prioripi=0.83,Double_t priorie=0.001);

  virtual ~AliITSPident();
  Double_t GetCondFunPro(Int_t lay) const {
    return fCondFunProLay[lay];
  }
  Double_t GetProdCondFunPro() const;
  Double_t GetCondFunK(Int_t lay) const {
    return fCondFunKLay[lay];
  }
  Double_t GetProdCondFunK() const;
  Double_t GetCondFunPi(Int_t lay) const {
    return fCondFunPiLay[lay];
  }
  Double_t GetProdCondFunPi() const;

  Double_t GetPBayesp()const {return fPBayesp;}
  Double_t GetPBayesk()const {return fPBayesk;}
  Double_t GetPBayespi()const {return fPBayespi;}
  Double_t GetPPriorip() const {return fPPriorip;}
  Double_t GetPPriorik() const {return fPPriorik;}
  Double_t GetPPrioripi() const {return fPPrioripi;}
  Double_t GetPPriorie() const {return fPPriorie;}
  
 private:

  void CalculateResponses(Double_t mom,AliITSPidParams *pars, Double_t *Qlay);
  Double_t CookCombinedBayes(Double_t condfun[][3],Double_t *prior,Int_t part)const;
  Double_t CookProd(Double_t condfun[][3],Int_t part) const;
  Double_t CookSum(Double_t condfun[][3],Double_t *prior) const;

  AliITSPident(const AliITSPident &ob); // copy constructor
  AliITSPident& operator=(const AliITSPident & ob); // ass. op.

  Double_t fCondFunProLay[8];     // one for each silicon layer
  Double_t fCondFunKLay[8];       // cond. prob. function kaons per layer
  Double_t fCondFunPiLay[8];      // cond. prob. function pions per layer
  Double_t fPBayesp;               // Bayes prob. 
  Double_t fPBayesk;               // Bayes prob. for kaons
  Double_t fPBayespi;              // Bayes prob. for pions 
  Double_t fPPriorip;              // Priori prob. 
  Double_t fPPriorik;              // Priori prob. for kaons
  Double_t fPPrioripi;             // Priori prob. for pions
  Double_t fPPriorie;              // Priori prob. for electrons

  ClassDef(AliITSPident,4);
};
#endif
