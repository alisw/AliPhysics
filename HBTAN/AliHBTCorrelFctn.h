#ifndef ALIHBTCORRELFUNCTION_H
#define ALIHBTCORRELFUNCTION_H

#include "AliHBTFunction.h"
#include "AliHBTParticle.h"
#include <Riostream.h>
//Set of functions:
//   Q Invaraint Correlation Function
//   Invariant Mass Function
//
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//Piotr.Skowronski@cern.ch

/*************************************************************************************/ 
class AliHBTQInvCorrelFctn: public AliHBTOnePairFctn1D
{
//Q Invaraint Correlation Function
//1D two particle function 
 public:
   AliHBTQInvCorrelFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTQInvCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair){return pair->GetQInv();}
  public:
    ClassDef(AliHBTQInvCorrelFctn,1)
 
};
/*************************************************************************************/ 

class AliHBTQOutCMSLCCorrelFctn: public AliHBTOnePairFctn1D
{
//Q OutCMSLCaraint Correlation Function
//1D two particle function 
 public:
   AliHBTQOutCMSLCCorrelFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTQOutCMSLCCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair){return TMath::Abs(pair->GetQOutCMSLC());}
  public:
    ClassDef(AliHBTQOutCMSLCCorrelFctn,1)
 
};
/*************************************************************************************/ 

class AliHBTQLongCMSLCCorrelFctn: public AliHBTOnePairFctn1D
{
//Q LongCMSLCaraint Correlation Function
//1D two particle function 
 public:
   AliHBTQLongCMSLCCorrelFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTQLongCMSLCCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair){return TMath::Abs(pair->GetQLongCMSLC());}
  public:
    ClassDef(AliHBTQLongCMSLCCorrelFctn,1)
 
};
/*************************************************************************************/ 

class AliHBTQSideCMSLCCorrelFctn: public AliHBTOnePairFctn1D
{
//Q SideCMSLCaraint Correlation Function
//1D two particle function 
 public:
   AliHBTQSideCMSLCCorrelFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTQSideCMSLCCorrelFctn(){}
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair){return TMath::Abs(pair->GetQSideCMSLC());}
  public:
    ClassDef(AliHBTQSideCMSLCCorrelFctn,1)
 
};
/*************************************************************************************/ 

class AliHBTInvMassCorrelFctn: public AliHBTOnePairFctn1D
{
//   Invariant Mass Function 
 public:
   AliHBTInvMassCorrelFctn(Int_t nbins = 2000, Double_t maxXval = 2., Double_t minXval = 0.0);
   virtual ~AliHBTInvMassCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair) { return pair->GetInvMass();}
  public:
    ClassDef(AliHBTInvMassCorrelFctn,1)
 
};

/*************************************************************************************/ 

class AliHBTTwoKStarCorrelFctn: public AliHBTOnePairFctn1D
{
//   Correlation Function of 2*KStar
 public:
   AliHBTTwoKStarCorrelFctn(Int_t nbins = 200, Double_t maxXval = 0.15, Double_t minXval = 0.0);
   virtual ~AliHBTTwoKStarCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair) { return 2.0*pair->GetKStar();}
  public:
    ClassDef(AliHBTTwoKStarCorrelFctn,1)
 
};


#endif
