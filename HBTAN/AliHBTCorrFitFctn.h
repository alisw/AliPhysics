#ifndef ALIHBTCORRFITFCTN_H
#define ALIHBTCORRFITFCTN_H
//____________________________________________________________
///////////////////////////////////////////////////////////////
//                                                           //
// class AliHBTCorrFitFctn                                   //
//                                                           //
//                                                           //
///////////////////////////////////////////////////////////////

#include "AliHBTFunction.h"
#include <TNtuple.h>

class AliHBTCorrFitFctn: public AliHBTTwoPairFctn1D
{
//Q Invaraint Correlation Function
//It writes Ntuple that is input for CorrFit
 public:
   AliHBTCorrFitFctn(Int_t nbins = 300, Double_t maxXval = 0.3, Double_t minXval = 0.0);
   virtual ~AliHBTCorrFitFctn(){delete fNtuple;}
   void ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
   void ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair);
   
   TH1* GetResult();
   void WriteFunction();
 protected:
   Double_t GetValue(AliHBTPair* /*trackpair*/, AliHBTPair* /*partpair*/) const {return 0.0;}//not usable
   
   TNtuple* fNtuple;//ntuple for storig pairs
   Int_t    fNPairsFitArea;//number of pairs in fitting area
   Int_t    fNPairsNormArea;//number of pairs in normalization area
  public:
   ClassDef(AliHBTCorrFitFctn,1)
};

#endif
