#include "AliHBTPair.h"
#include "AliHBTParticle.h"

ClassImp(AliHBTPair)

/************************************************************************/
AliHBTPair::AliHBTPair(Bool_t rev):
 fPart1(0x0),
 fPart2(0x0),
 fSwapedPair(0x0),
 fQSideCMSLC(0.0),
 fQSideCMSLCNotCalc(kTRUE),
 fQOutCMSLC(0.0),
 fQOutCMSLCNotCalc(kTRUE),
 fQLongCMSLC(0.0),
 fQLongCMSLCNotCalc(kTRUE),
 fQInv(0.0),
 fQInvNotCalc(kTRUE),
 fInvMass(0.0),
 fInvMassNotCalc(kTRUE),
 fKt(0.0),
 fKtNotCalc(kTRUE),
 fKStar(0.0),
 fKStarNotCalc(kTRUE),
 fPInv(0.0),
 fQSide(0.0),
 fOut(0.0),
 fQLong(0.0),
 fMt(0.0),
 fMtNotCalc(kTRUE),
 fInvMassSqr(0.0),
 fMassSqrNotCalc(kTRUE),
 fQInvL(0.0),
 fQInvLNotCalc(kTRUE),
 fPxSum(0.0),
 fPySum(0.0),
 fPzSum(0.0),
 fESum(0.0),
 fSumsNotCalc(kTRUE),
 fPxDiff(0.0),
 fPyDiff(0.0),
 fPzDiff(0.0),
 fEDiff(0.0),
 fDiffsNotCalc(kTRUE),
 fGammaCMSLC(0.0),
 fGammaCMSLCNotCalc(kTRUE),
 fChanged(kTRUE)
 {
//value of rev defines if it is Swaped
//if you pass kTRUE swpaped pair will NOT be created
//though you wont be able to get the swaped pair from this pair

  if(!rev) fSwapedPair = new AliHBTPair(kTRUE); //if false create swaped pair
  
 }
/************************************************************************/

AliHBTPair::AliHBTPair(AliHBTParticle* part1, AliHBTParticle* part2, Bool_t rev):
 fPart1(part1),
 fPart2(part2),
 fSwapedPair(0x0),
 fQSideCMSLC(0.0),
 fQSideCMSLCNotCalc(kTRUE),
 fQOutCMSLC(0.0),
 fQOutCMSLCNotCalc(kTRUE),
 fQLongCMSLC(0.0),
 fQLongCMSLCNotCalc(kTRUE),
 fQInv(0.0),
 fQInvNotCalc(kTRUE),
 fInvMass(0.0),
 fInvMassNotCalc(kTRUE),
 fKt(0.0),
 fKtNotCalc(kTRUE),
 fKStar(0.0),
 fKStarNotCalc(kTRUE),
 fPInv(0.0),
 fQSide(0.0),
 fOut(0.0),
 fQLong(0.0),
 fMt(0.0),
 fMtNotCalc(kTRUE),
 fInvMassSqr(0.0),
 fMassSqrNotCalc(kTRUE),
 fQInvL(0.0),
 fQInvLNotCalc(kTRUE),
 fPxSum(0.0),
 fPySum(0.0),
 fPzSum(0.0),
 fESum(0.0),
 fSumsNotCalc(kTRUE),
 fPxDiff(0.0),
 fPyDiff(0.0),
 fPzDiff(0.0),
 fEDiff(0.0),
 fDiffsNotCalc(kTRUE),
 fGammaCMSLC(0.0),
 fGammaCMSLCNotCalc(kTRUE),
 fChanged(kTRUE)
 {
//value of rev defines if it is Swaped
//if you pass kTRUE swpaped pair will NOT be created
//though you wont be able to get the swaped pair from this pair

  if(!rev) fSwapedPair = new AliHBTPair(part2,part1,kTRUE); //if false create swaped pair
  
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

   if (fPart1->GetMass() != fPart2->GetMass())
    {
/*    
      //STAR algorithm
      Double_t beta  = fPzSum/fESum;
      Double_t gamma = GetGammaToCMSLC();
      Double_t el = gamma * (fPart1->Energy() - beta * fPart1->Pz());
      Double_t x  = ( fPart1->Px()*fPxSum + fPart1->Py()*fPySum) / ( 2.0*GetKt() );
      beta  = 2.0*GetKt()/GetMt();
      gamma = GetMt()/GetQInv();
      fQOutCMSLC = gamma * (x - beta * el);
*/

      //beta=fPzSum/fESum;    // Longit. V == beta
      Double_t beta=fPzSum/fESum;
      Double_t gamma = GetGammaToCMSLC();
      
      Double_t cosphi=fPxSum/(2.0*GetKt());  // cos(phi)
      Double_t sinphi=fPySum/(2.0*GetKt()); // sin(phi)
      
//      ROTATE(part1Px,part1Py,SPHI,CPHI,part1Px,part1Py);//ROT8
//      ROTATE(part2Px,part2Py,SPHI,CPHI,part2Px,part2Py);//ROT8
      Double_t tmp;
      tmp = fPart1->Px()*cosphi + fPart1->Py()*sinphi;
      Double_t part1Py = fPart1->Py()*cosphi - fPart1->Px()*sinphi;
      Double_t part1Px = tmp;

      tmp = fPart2->Px()*cosphi + fPart2->Py()*sinphi;
      Double_t part2Py = fPart2->Py()*cosphi - fPart2->Px()*sinphi;
      Double_t part2Px = tmp;
      
      
//      LTR(part1Pz,E1,beta,GetGammaToCMSLC(),part1Pz,E1a);
//      LTR(part2Pz,E2,beta,GetGammaToCMSLC(),part2Pz,E2a);
      Double_t part1Pz=gamma*(fPart1->Pz()-beta*fPart1->Energy());
      Double_t part2Pz=gamma*(fPart2->Pz()-beta*fPart2->Energy());

      Double_t part1P2=part1Px*part1Px+part1Py*part1Py+part1Pz*part1Pz;
      Double_t part2P2=part2Px*part2Px+part2Py*part2Py+part2Pz*part2Pz;
      Double_t part1E=TMath::Sqrt(fPart1->GetMass()*fPart1->GetMass()+part1P2);
      Double_t part2E=TMath::Sqrt(fPart2->GetMass()*fPart2->GetMass()+part2P2);
      Double_t sumE=part1E+part2E;
      Double_t sumPx=part1Px+part2Px;
      Double_t sumPy=part1Py+part2Py;
      Double_t sumPZ=part1Pz+part2Pz;
      Double_t sumP2=sumPx*sumPx+sumPy*sumPy+sumPZ*sumPZ;

      Double_t relmass=TMath::Sqrt(sumE*sumE-sumP2);
      Double_t hf = (fPart1->GetMass()*fPart1->GetMass() - fPart2->GetMass()*fPart2->GetMass())/(relmass*relmass);
      fQOutCMSLC=(part1Px-part2Px);//== id
      fQOutCMSLC=fQOutCMSLC-sumPx*hf; //sumPx == fPxSum ale po rotacji i transf
    }
   else
    {
      Double_t k2 = fPxSum*fPxDiff+fPySum*fPyDiff;
      fQOutCMSLC = 0.5*k2/GetKt();
   // if (non-id)  fQOutCMSLC=fQOutCMSLC - sumPx*HF;
    }

    
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
    fQLongCMSLC = GetGammaToCMSLC() * ( fPzDiff - beta*fEDiff );
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

    Double_t Q = (fPart1->GetMass()*fPart1->GetMass() - fPart2->GetMass()*fPart2->GetMass())/Pinv;
    
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
//returns Qinv 
//warning for non-id particles you want to use 2*KStar
  if(fQInvNotCalc)
   {
    CalculateQInvL();
    fQInv = TMath::Sqrt(TMath::Abs(fQInvL));
    fQInvNotCalc = kFALSE;
   }
  return fQInv;
}
/************************************************************************/

Double_t AliHBTPair::GetGammaToCMSLC()
{
  if(fGammaCMSLCNotCalc)
   {
     CalculateSums();
     Double_t beta = fPzSum/fESum;
     fGammaCMSLC = 1.0/TMath::Sqrt(1.0 - beta*beta);
     fGammaCMSLCNotCalc = kFALSE;
   }
  return fGammaCMSLC;
}
/************************************************************************/

Double_t AliHBTPair::GetMt()
{
  if (fMtNotCalc)
   {
     CalculateSums();
     fMt = TMath::Sqrt(fESum*fESum - fPzSum*fPzSum);
     fMtNotCalc = kFALSE;
   }
   return fMt;
}
/************************************************************************/







