#ifndef ALIHBTWeightQINVFCTN_H
#define ALIHBTWeightQINVFCTN_H

/* $Id$ */
//_________________________________________________________________________
//
// class AliHBTWeightQInvFctn
// class AliHBTWeightQOutFctn
// class AliHBTWeightQSideFctn
// class AliHBTWeightQLongFctn
// 
// This class allows to obtain Q_inv correlation function with weights
// calculated by Lednicky's alghorithm.
// Numerator is filled with weighted events. Weights are attributed to reconstructed tracks.
// Weights are calculated with corresponding simulated particles momenta.
// Denominator is filled with mixing unweighted reconstructed tracks.
// One needs both pairs 
// (simulated and recontructed), thus function is of class AliHBTTwoPairFctn1D.
// Author: Ludmila Malinina, JINR (malinina@sunhe.jinr.ru)
//
////////////////////////////////////////////////////////////////////////////////

#include "AliHBTFunction.h"


class AliHBTWeights;

class AliHBTWeightQInvFctn: public AliHBTTwoPairFctn1D, public AliHBTCorrelFunction
{
 public:
  AliHBTWeightQInvFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual  ~AliHBTWeightQInvFctn(){};
  TH1* GetResult(); 
  
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) const
    { return trackpair->GetQInv()-partpair->GetQInv();} //isn't use                                                                    
  ClassDef(AliHBTWeightQInvFctn,3)
};
/*************************************************************************************/ 

class AliHBTWeightQOutFctn: public AliHBTTwoPairFctn1D, public AliHBTCorrelFunction
{

  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTWeightQOutFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual ~AliHBTWeightQOutFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) const
    { return trackpair->GetQOutLCMS()-partpair->GetQOutLCMS();} //isn't use                                                                    
  ClassDef(AliHBTWeightQOutFctn,3)
 
};
/*************************************************************************************/ 
  
class AliHBTWeightQLongFctn: public AliHBTTwoPairFctn1D, public AliHBTCorrelFunction
{
  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTWeightQLongFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual ~AliHBTWeightQLongFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) const
    { return trackpair->GetQLongLCMS()-partpair->GetQLongLCMS();} //isn't used

  ClassDef(AliHBTWeightQLongFctn,3)
 
};
/*************************************************************************************/ 
  
class AliHBTWeightQSideFctn: public AliHBTTwoPairFctn1D, public AliHBTCorrelFunction
{
  // friend class AliHBTOnePairFctn1D;
 public:
  AliHBTWeightQSideFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual ~AliHBTWeightQSideFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) const
    { return trackpair->GetQLongLCMS()-partpair->GetQLongLCMS();} //isn't used

  ClassDef(AliHBTWeightQSideFctn,3) 
};
/*************************************************************************************/ 
  
class AliHBTWeightTwoKStarFctn: public AliHBTTwoPairFctn1D, public AliHBTCorrelFunction
{
  // friend class AliHBTOnePairFctn1D;
 public:
  AliHBTWeightTwoKStarFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
  virtual ~AliHBTWeightTwoKStarFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) const
    { return trackpair->GetKStar()-partpair->GetKStar();} //isn't used
  ClassDef(AliHBTWeightTwoKStarFctn,3) 

};
/*************************************************************************************/ 

class AliHBTWeightQOutQSideFctn: public AliHBTTwoPairFctn2D, public AliHBTCorrelFunction
{

  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTWeightQOutQSideFctn(Int_t nxbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                              Int_t nybins = 100, Double_t maxYval = 0.15, Double_t minYval = 0.0);
  virtual ~AliHBTWeightQOutQSideFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  void GetValues(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/, Double_t& /*x*/, Double_t& /*y*/) const{}
  ClassDef(AliHBTWeightQOutQSideFctn,3)
 
};
/*************************************************************************************/ 

class AliHBTWeightQOutQLongFctn: public AliHBTTwoPairFctn2D, public AliHBTCorrelFunction
{

  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTWeightQOutQLongFctn(Int_t nxbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                              Int_t nybins = 100, Double_t maxYval = 0.15, Double_t minYval = 0.0);
  virtual ~AliHBTWeightQOutQLongFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  void GetValues(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/, Double_t& /*x*/, Double_t& /*y*/) const{}
  ClassDef(AliHBTWeightQOutQLongFctn,3)
 
};
/*************************************************************************************/ 

class AliHBTWeightQSideQLongFctn: public AliHBTTwoPairFctn2D, public AliHBTCorrelFunction
{

  //  friend class AliHBTOnePairFctn1D;
 public:
  AliHBTWeightQSideQLongFctn(Int_t nxbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0,
                             Int_t nybins = 100, Double_t maxYval = 0.15, Double_t minYval = 0.0);
  virtual ~AliHBTWeightQSideQLongFctn(){};
  TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
 protected:
  void GetValues(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/, Double_t& /*x*/, Double_t& /*y*/) const{}
  ClassDef(AliHBTWeightQSideQLongFctn,3)
 
};
/*************************************************************************************/ 

class AliHBTWeightQOutSQideQLongFctn: public AliHBTTwoPairFctn3D, public AliHBTCorrelFunction
{

  public:
    AliHBTWeightQOutSQideQLongFctn(Int_t nXbins = 100, Double_t maxXval = 0.3, Double_t minXval = 0.0,
                                Int_t nYbins = 100, Double_t maxYval = 0.3, Double_t minYval = 0.0,
                                Int_t nZbins = 100, Double_t maxZval = 0.3, Double_t minZval = 0.0);
    virtual  ~AliHBTWeightQOutSQideQLongFctn(){}

    TH1* GetResult();
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  
  protected:
    void GetValues(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/, 
                   Double_t& /*x*/, Double_t& /*y*/, Double_t& /*z*/)  const{ } 

    ClassDef(AliHBTWeightQOutSQideQLongFctn,1)
};


#endif
