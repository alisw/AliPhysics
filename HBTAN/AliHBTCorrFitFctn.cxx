#include "AliHBTCorrFitFctn.h"
//____________________________________________________________
///////////////////////////////////////////////////////////////
//                                                           //
// class AliHBTCorrFitFctn                                   //
//                                                           //
//                                                           //
///////////////////////////////////////////////////////////////

ClassImp(AliHBTCorrFitFctn)

/*****************************************************************/

AliHBTCorrFitFctn::AliHBTCorrFitFctn(Int_t nbins, Double_t maxXval, Double_t minXval):
 AliHBTTwoPairFctn1D(nbins,maxXval,minXval),
 fNtuple(new TNtuple("pair", "pair", "px1:py1:pz1:e1:px2:py2:pz2:e2")),
 fNPairsFitArea(0),
 fNPairsNormArea(0)
{
//ctor
 fWriteNumAndDen = kTRUE;//change default behaviour
 Rename("wqinvcfCorrFit","Lednicky Weught Theoretical Q_{inv} Correlation Function");
} 
/*****************************************************************/

void AliHBTCorrFitFctn::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
 //Fills the numerator using pair from the same event
   partpair = CheckPair(partpair);
   if(partpair == 0x0) return; 
   trackpair = CheckPair(trackpair);
   if(trackpair == 0x0) return; 
   
   Double_t q = trackpair->GetQInv();
   Bool_t fill = kFALSE;
   
   Double_t weight = partpair->GetWeight();
   fNumerator->Fill(q,weight);
   
   if ( (q < 0.15) && (fNPairsFitArea < 2.e+5))
     {
       fNPairsFitArea++;
       fill = kTRUE;
     }

   if ( (q > 0.15) && (q < 0.3) && (fNPairsFitArea < 1.e+5))
     {
       fNPairsNormArea++;
       fill = kTRUE;
     }
   
   if (fill)
    {  
      const AliVAODParticle& p1 = *(trackpair->Particle1());
      const AliVAODParticle& p2 = *(trackpair->Particle2());
      fNtuple->Fill(p1.Px(),p1.Py(),p1.Pz(),p1.E(),
                    p2.Px(),p2.Py(),p2.Pz(),p2.E());
    }
}
/****************************************************************/

void  AliHBTCorrFitFctn::ProcessDiffEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  // Fills the denominator using mixed pairs
  trackpair = CheckPair(trackpair);
  partpair  = CheckPair(partpair);
  if ( trackpair && partpair)
  {
     fDenominator->Fill(trackpair->GetQInv());
  }
}
/*****************************************************************/

TH1* AliHBTCorrFitFctn::GetResult()
{
//returns ratio of numerator and denominator
 return GetRatio(Scale());
}
/**************************************************************/

void AliHBTCorrFitFctn::WriteFunction()
{
  //writes a function 
  AliHBTFunction::WriteFunction();
  fNtuple->Write(0,TObject::kOverwrite);
}
