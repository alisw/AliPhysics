#ifndef ALIHBTCORRELFUNCTION_H
#define ALIHBTCORRELFUNCTION_H

#include "AliHBTFunction.h"
#include "AliHBTParticle.h"
//Set of functions:
//   Q Invaraint Correlation Function
//   Invariant Mass Function
//
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//Piotr.Skowronski@cern.ch


class AliHBTQInvCorrelFctn: public AliHBTTwoPartFctn1D
{
//Q Invaraint Correlation Function
//1D two particle function 
 public:
   AliHBTQInvCorrelFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0):
                        AliHBTTwoPartFctn1D(nbins,maxXval,minXval){}
   virtual ~AliHBTQInvCorrelFctn(){};
   TH1* GetResult();
 protected:
   Double_t GetValue(AliHBTPair * pair){return pair->GetQInv();}
  public:
    ClassDef(AliHBTQInvCorrelFctn,1)
 
};


class AliHBTInvMassCorrelFctn: public AliHBTTwoPartFctn1D
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



#endif
