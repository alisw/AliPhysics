#ifndef ALIHBTLLWEIGHTQINVFCTN_H
#define ALIHBTLLWEIGHTQINVFCTN_H

/* $Id$ */


#include "AliHBTFunction.h"


class AliHBTLLWeights;

class AliHBTLLWeightQInvFctn: public AliHBTTwoPairFctn1D
{
  friend class AliHBTOnePairFctn1D;

  public:
      AliHBTLLWeightQInvFctn(Int_t nbins = 100, Double_t maxXval = 0.15, Double_t minXval = 0.0);
      virtual  ~AliHBTLLWeightQInvFctn(){};
      TH1* GetResult(); 

      void   ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      void   ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
      
      Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair) const
	{ return trackpair->GetQInv()-partpair->GetQInv();} //isn't use                                                                    
  private:

     ClassDef(AliHBTLLWeightQInvFctn,1)
};
  
#endif
