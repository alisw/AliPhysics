#ifndef ALIHBTLLWEIGHTQINVFCTN_H
#define ALIHBTLLWEIGHTQINVFCTN_H

/* $Id$ */


#include "AliHBTFunction.h"


class AliHBTLLWeights;

class AliHBTLLWeightQInvFctn: public AliHBTTwoPairFctn1D
{
 public:
  AliHBTLLWeightQInvFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual  ~AliHBTLLWeightQInvFctn(){};
  TH1* GetResult(); 
  
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) const
    { return trackpair->GetQInv()-partpair->GetQInv();} //isn't use                                                                    
  ClassDef(AliHBTLLWeightQInvFctn,1)
};
/*************************************************************************************/ 

class AliHBTLLWeightQOutFctn: public AliHBTTwoPairFctn1D
{

  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightQOutFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0):
    AliHBTTwoPairFctn1D(nbins,maxXval,minXval){}
  virtual ~AliHBTLLWeightQOutFctn(){};
  TH1* GetResult();
 protected:
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)         
    { return trackpair->GetQOutCMSLC()-partpair->GetQOutCMSLC();} //isn't use                                                                    
  ClassDef(AliHBTLLWeightQOutFctn,1)
 
};
/*************************************************************************************/ 
class AliHBTLLWeightQLongFctn: public AliHBTTwoPairFctn1D
{
  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightQLongFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0):
    AliHBTTwoPairFctn1D(nbins,maxXval,minXval){}
  virtual ~AliHBTLLWeightQLongFctn(){};
  TH1* GetResult();
 protected:
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)
    { return trackpair->GetQLongCMSLC()-partpair->GetQLongCMSLC();} //isn't used

  ClassDef(AliHBTLLWeightQLongFctn,1)
 
};
/*************************************************************************************/ 
class AliHBTLLWeightQSideFctn: public AliHBTTwoPairFctn1D
{
  // friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightQSideFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0):
    AliHBTTwoPairFctn1D(nbins,maxXval,minXval){}
  virtual ~AliHBTLLWeightQSideFctn(){};
  TH1* GetResult();
 protected:
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)         
    { return trackpair->GetQLongCMSLC()-partpair->GetQLongCMSLC();} //isn't used

  ClassDef(AliHBTLLWeightQSideFctn,1) 
};
/*************************************************************************************/ 
class AliHBTLLWeightTwoKStarFctn: public AliHBTTwoPairFctn1D
{
  // friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightTwoKStarFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0):
    AliHBTTwoPairFctn1D(nbins,maxXval,minXval){}
  virtual ~AliHBTLLWeightTwoKStarFctn(){};
  TH1* GetResult();
 protected:
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)         
    { return trackpair->GetKStar()-partpair->GetKStar();} //isn't used
  ClassDef(AliHBTLLWeightTwoKStarFctn,1) 

};
#endif
