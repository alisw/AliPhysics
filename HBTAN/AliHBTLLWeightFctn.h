//This function allows to obtain Q_inv correlation function with weights
//calculated by Lednicky's alghorithm.
//Numerator is filled with weighted events. Weights are attributed to reconstructed tracks.
//Weights are calculated with corresponding simulated particles momenta.
//Denominator is filled with mixing unweighted reconstructed tracks.
//One needs both pairs 
//(simulated and recontructed), thus function is of class AliHBTTwoPairFctn1D.

#ifndef ALIHBTLLWEIGHTFCTN_H
#define ALIHBTLLWEIGHTFCTN_H
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
      
      Double_t GetValue(AliHBTPair* trackpair, AliHBTPair* partpair)         
           { return trackpair->GetQInv()-partpair->GetQInv();} //isn't use                                                                    
	

  protected:

  private:
  public:
     ClassDef(AliHBTLLWeightQInvFctn,1)
};
  
#endif
