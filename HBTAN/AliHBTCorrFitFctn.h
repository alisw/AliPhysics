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

class AliHBTCorrFitFctn: public AliHBTOnePairFctn1D
{
//Q Invaraint Correlation Function
//It writes Ntuple that is input for CorrFit
 public:
   AliHBTCorrFitFctn();
   AliHBTCorrFitFctn(Int_t  fit, Int_t  norm);
   virtual ~AliHBTCorrFitFctn(){delete fNtuple;}
   void ProcessSameEventParticles(AliHBTPair* trackpair);
   void ProcessDiffEventParticles(AliHBTPair* trackpair);
   
   TH1* GetResult();
   void WriteFunction();

   void SetMaxNumberOfPairs(Int_t  fit, Int_t  norm){fNMaxPairsFitArea = fit; fNMaxPairsNormArea = norm;}
   void SetFitRange(Float_t max) {fFitRangeMax = max;}
   void SetNormalizationRange(Float_t min, Float_t max) { fNormRangeMin = min; fNormRangeMax= max;}
 protected:
   Double_t GetValue(AliHBTPair* /*trackpair*/) const {return 0.0;}//not usable
   
   TNtuple* fNtuple;//ntuple for storig pairs

   Int_t    fNPairsFitArea;//current number of pairs in fitting area
   Int_t    fNMaxPairsFitArea;//current number of pairs in fitting area
   Float_t  fFitRangeMax;
   
   Int_t    fNPairsNormArea;//number of pairs in normalization area
   Int_t    fNMaxPairsNormArea;//number of pairs in normalization area
   Float_t  fNormRangeMin;
   Float_t  fNormRangeMax;
  public:
   ClassDef(AliHBTCorrFitFctn,1)
};

#endif
