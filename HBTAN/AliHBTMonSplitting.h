#ifndef ALIHBTMONSPLITTING_H
#define ALIHBTMONSPLITTING_H

#include "AliHBTFunction.h"

class AliHBTMonSplittingQosl: public AliHBTTwoPairFctn3D, public AliHBTCorrelFunction
{
//histograms qout-qside-qlong of splitted tracks - found more than ones

  public:
   AliHBTMonSplittingQosl(Int_t nXbins = 50, Double_t maxXval = 0.05, Double_t minXval = -0.05,
                          Int_t nYbins = 50, Double_t maxYval = 0.05, Double_t minYval = -0.05,
                          Int_t nZbins = 50, Double_t maxZval = 0.05, Double_t minZval = -0.05 );
  virtual  ~AliHBTMonSplittingQosl(){}

  TH1*   GetResult() {return fNumerator;}
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/) {}
  
  protected:
    void GetValues(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/, 
                   Double_t& /*x*/, Double_t& /*y*/, Double_t& /*z*/)  const{ } 

    ClassDef(AliHBTMonSplittingQosl,1)
};


class AliHBTMonSplittingDptDthetaDphi: public AliHBTTwoPairFctn3D, public AliHBTCorrelFunction
{
//histograms qout-qside-qlong of splitted tracks - found more than ones

  public:
   AliHBTMonSplittingDptDthetaDphi(Int_t nXbins = 50, Double_t maxXval = 0.05, Double_t minXval = -0.05,
                          Int_t nYbins = 50, Double_t maxYval = 0.05, Double_t minYval = -0.05,
                          Int_t nZbins = 50, Double_t maxZval = 0.05, Double_t minZval = -0.05 );
  virtual  ~AliHBTMonSplittingDptDthetaDphi(){}

  TH1*   GetResult() {return fNumerator;}
  void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
  void   ProcessDiffEventParticles(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/) {}
  
  protected:
    void GetValues(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/, 
                   Double_t& /*x*/, Double_t& /*y*/, Double_t& /*z*/)  const{ } 

    ClassDef(AliHBTMonSplittingDptDthetaDphi,1)
};



#endif
