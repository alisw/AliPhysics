#ifndef ALIHBTLLWEIGHTTHEORFCTN_H
#define ALIHBTLLWEIGHTTHEORFCTN_H
/* $Id$ */

//This function allows to obtain Q_inv correlation function with weights
//calculated by Lednicky's alghorithm.
//Numerator is filled with weighted events. Weights are attributed to simulated particles.
//Weights are calculated with corresponding simulated particles momenta.
//Denominator is filled with mixing unweighted simulated particles.
//One needs only simulated pairs, so 
//this function is of class AliHBTOnePairFctn1D.
//Author Ludmila Malinina JINR (malinina@sunhe.jinr.ru)

#include "AliHBTFunction.h"

class AliHBTLLWeights;

class AliHBTLLWeightTheorQInvFctn: public AliHBTOnePairFctn1D
{

  public:
  AliHBTLLWeightTheorQInvFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual  ~AliHBTLLWeightTheorQInvFctn(){}

  TH1* GetResult(); 
  void   ProcessSameEventParticles(AliHBTPair* partpair);
      
  Double_t GetValue(AliHBTPair* partpair)
    { return partpair->GetQInv();} //isn't used

  ClassDef(AliHBTLLWeightTheorQInvFctn,1)
};

class AliHBTLLWeightTheorQOutFctn: public AliHBTOnePairFctn1D
{

  public:
  AliHBTLLWeightTheorQOutFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual  ~AliHBTLLWeightTheorQOutFctn(){}

  TH1* GetResult(); 
  void   ProcessSameEventParticles(AliHBTPair* partpair);
      
  Double_t GetValue(AliHBTPair* partpair)
    { return partpair->GetQOutCMSLC();} //isn't used

  ClassDef(AliHBTLLWeightTheorQOutFctn,1)
};

class AliHBTLLWeightTheorQSideFctn: public AliHBTOnePairFctn1D
{

  public:
  AliHBTLLWeightTheorQSideFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual  ~AliHBTLLWeightTheorQSideFctn(){}

  TH1* GetResult(); 
  void   ProcessSameEventParticles(AliHBTPair* partpair);
      
  Double_t GetValue(AliHBTPair* partpair)
    { return partpair->GetQSideCMSLC();} //isn't used

  ClassDef(AliHBTLLWeightTheorQSideFctn,1)
};

class AliHBTLLWeightTheorQLongFctn: public AliHBTOnePairFctn1D
{

  public:
  AliHBTLLWeightTheorQLongFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual  ~AliHBTLLWeightTheorQLongFctn(){}

  TH1* GetResult(); 
  void   ProcessSameEventParticles(AliHBTPair* partpair);
      
  Double_t GetValue(AliHBTPair* partpair)
    { return partpair->GetQLongCMSLC();} //isn't used

  ClassDef(AliHBTLLWeightTheorQLongFctn,1)
};
    
#endif
