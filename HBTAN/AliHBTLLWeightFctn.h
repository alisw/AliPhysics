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
  
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)
    { return trackpair->GetQInv()-partpair->GetQInv();} //isn't use                                                                    
  ClassDef(AliHBTLLWeightQInvFctn,1)
};
/*************************************************************************************/ 

class AliHBTLLWeightQOutFctn: public AliHBTTwoPairFctn1D
{

  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightQOutFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual ~AliHBTLLWeightQOutFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)         
    { return trackpair->GetQOutCMSLC()-partpair->GetQOutCMSLC();} //isn't use                                                                    
  ClassDef(AliHBTLLWeightQOutFctn,1)
 
};
/*************************************************************************************/ 
  
class AliHBTLLWeightQLongFctn: public AliHBTTwoPairFctn1D
{
  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightQLongFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual ~AliHBTLLWeightQLongFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)
    { return trackpair->GetQLongCMSLC()-partpair->GetQLongCMSLC();} //isn't used

  ClassDef(AliHBTLLWeightQLongFctn,1)
 
};
/*************************************************************************************/ 
  
class AliHBTLLWeightQSideFctn: public AliHBTTwoPairFctn1D
{
  // friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightQSideFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual ~AliHBTLLWeightQSideFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)         
    { return trackpair->GetQLongCMSLC()-partpair->GetQLongCMSLC();} //isn't used

  ClassDef(AliHBTLLWeightQSideFctn,1) 
};
/*************************************************************************************/ 
  
class AliHBTLLWeightTwoKStarFctn: public AliHBTTwoPairFctn1D
{
  // friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightTwoKStarFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual ~AliHBTLLWeightTwoKStarFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)         
    { return trackpair->GetKStar()-partpair->GetKStar();} //isn't used
  ClassDef(AliHBTLLWeightTwoKStarFctn,1) 

};
/*************************************************************************************/ 

class AliHBTLLWeightQOutQSideFctn: public AliHBTTwoPairFctn2D
{

  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightQOutQSideFctn(Int_t nxbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                              Int_t nybins = 100, Double_t maxYval = 0.15, Double_t minYval = 0.0);
  virtual ~AliHBTLLWeightQOutQSideFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  void GetValues(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/, Double_t& /*x*/, Double_t& /*y*/){}
  ClassDef(AliHBTLLWeightQOutQSideFctn,1)
 
};
/*************************************************************************************/ 

class AliHBTLLWeightQOutQLongFctn: public AliHBTTwoPairFctn2D
{

  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightQOutQLongFctn(Int_t nxbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                              Int_t nybins = 100, Double_t maxYval = 0.15, Double_t minYval = 0.0);
  virtual ~AliHBTLLWeightQOutQLongFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  void GetValues(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/, Double_t& /*x*/, Double_t& /*y*/){}
  ClassDef(AliHBTLLWeightQOutQLongFctn,1)
 
};

/*************************************************************************************/ 

class AliHBTLLWeightQSideQLongFctn: public AliHBTTwoPairFctn2D
{

  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTLLWeightQSideQLongFctn(Int_t nxbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                              Int_t nybins = 100, Double_t maxYval = 0.15, Double_t minYval = 0.0);
  virtual ~AliHBTLLWeightQSideQLongFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  void GetValues(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/, Double_t& /*x*/, Double_t& /*y*/){}
  ClassDef(AliHBTLLWeightQSideQLongFctn,1)
 
};

#endif
