#include "AliHBTPair.h"
#include "AliHBTParticle.h"

ClassImp(AliHBTPair)

/************************************************************************/
AliHBTPair::AliHBTPair(Bool_t rev):
 fPart1(0x0),
 fPart2(0x0)
 {
//value of rev defines if it is Swaped
//if you pass kTRUE swpaped pair will NOT be created
//though you wont be able to get the swaped pair from this pair

  if(!rev) fSwapedPair = new AliHBTPair(kTRUE); //if false create swaped pair
  else fSwapedPair = 0x0; //if true set the pointer to NULL
  
 }
/************************************************************************/

AliHBTPair::AliHBTPair(AliHBTParticle* part1, AliHBTParticle* part2, Bool_t rev):
 fPart1(part1),
 fPart2(part2)
 {
//value of rev defines if it is Swaped
//if you pass kTRUE swpaped pair will NOT be created
//though you wont be able to get the swaped pair from this pair

  if(!rev) fSwapedPair = new AliHBTPair(part2,part1,kTRUE); //if false create swaped pair
  else fSwapedPair = 0x0; //if true set the pointer to NULL
  
 }
/************************************************************************/

Double_t AliHBTPair::GetInvMass()
{
//Returns qinv value for a pair
  if(fInvMassNotCalc)
   {
     CalculateInvMassSqr(); //method is inline so we not waste th time for jumping into method 
     
     if(fInvMassSqr<0)  fInvMass = TMath::Sqrt(-fInvMassSqr);
     else fInvMass = TMath::Sqrt(fInvMassSqr); 
     
     fInvMassNotCalc = kFALSE;
   }
  return fInvMass;
}
/************************************************************************/
Double_t AliHBTPair::GetQSideCMSLC()
{
  //return Q Side in Central Of Mass System in Longitudialy Comoving Frame
 
  if (fQSideCMSLCNotCalc)
   {
    fQSideCMSLC = (fPart1->Px()*fPart2->Py()-fPart2->Px()*fPart1->Py())/GetKt();
    fQSideCMSLCNotCalc = kFALSE;
   }
  return fQSideCMSLC;
}
/************************************************************************/
Double_t AliHBTPair::GetQOutCMSLC()
{
 if(fQOutCMSLCNotCalc)
  {
   CalculateSums();
   CalculateDiffs();
   Double_t k2 = fPxSum*fPxDiff+fPySum*fPyDiff;
   fQOutCMSLC = 0.5*k2/GetKt();
   fQOutCMSLCNotCalc = kFALSE;
  }
 return fQOutCMSLC;
}
/************************************************************************/
Double_t AliHBTPair::GetQLongCMSLC()
{
 if (fQLongCMSLCNotCalc)
  {
    CalculateSums();
    CalculateDiffs();
    Double_t beta = fPzSum/fESum;
    Double_t gamma = 1.0/TMath::Sqrt(1.0 - beta*beta);
    fQLongCMSLC = gamma*( fPzDiff - beta*fEDiff );
    fQLongCMSLCNotCalc = kFALSE;
  }
 return fQLongCMSLC; 
}
/************************************************************************/
Double_t AliHBTPair::GetKt()
{
  if(fKtNotCalc)
   { 
     CalculateSums();
     fKt =  0.5*TMath::Hypot(fPxSum,fPySum);
     fKtNotCalc = kFALSE;
   }
  return fKt;
}
/************************************************************************/

Double_t AliHBTPair::GetKStar()
{
  if (fKStarNotCalc)
   { 
    
    CalculateSums();

    Double_t Ptrans = fPxSum*fPxSum + fPySum*fPySum;
    Double_t Mtrans = fESum*fESum - fPzSum*fPzSum;
    Double_t Pinv =   TMath::Sqrt(Mtrans - Ptrans);

    Double_t Q = ( fPart1->GetMass()*fPart1->GetMass() - fPart2->GetMass()*fPart2->GetMass())/Pinv;
    
    CalculateQInvL();
    
    Q = TMath::Sqrt( Q*Q - fQInvL);
    fKStar = Q/2.;
    fKStarNotCalc = kFALSE;
   }
  return fKStar;
}
/************************************************************************/

Double_t AliHBTPair::GetQInv()
{
  if(fQInvNotCalc)
   {
    CalculateQInvL();
    fQInv = TMath::Sqrt(TMath::Abs(fQInvL));
    fQInvNotCalc = kFALSE;
   }
  return fQInv;
}

/************************************************************************/
/************************************************************************/

/************************************************************************/







