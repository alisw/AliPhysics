/* $Id$ */

//------------------------------------------------------------------
// This class introduces the weight's calculation 
// according to the Lednicky's algorithm.
// The detailed description of the algorithm can be found 
// in comments to fortran code:
// fsiw.f, fsiini.f
// Author:
//------------------------------------------------------------------

#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom.h>

#include "AliHBTLLWeights.h"
#include "AliHBTPair.h"
#include "AliHBTParticle.h"
#include "WLedCOMMONS.h"

/*******************************************************************/
/******      ROUTINES    USED    FOR     COMMUNUCATION      ********/
/********************     WITH      FORTRAN     ********************/
/*******************************************************************/
#ifndef WIN32
# define led_bldata led_bldata_
# define fsiini fsiini_
# define ltran12 ltran12_
# define fsiw fsiw_
# define type_of_call
#else
# define led_bldata LED_BLDATA
# define fsiini FSIINI
# define ltran12 LTRAN12
# define fsiw FSIW
# define type_of_call _stdcall
#endif
/****************************************************************/
extern "C" void type_of_call led_bldata(); 
extern "C" void type_of_call fsiini();
extern "C" void type_of_call ltran12();
extern "C" void type_of_call fsiw();
/**************************************************************/

ClassImp(AliHBTLLWeights)  
 
AliHBTLLWeights* AliHBTLLWeights::fgLLWeights=NULL; 

AliHBTLLWeights::AliHBTLLWeights()
{
  // Default Constructor 
  fPID1 = 0;
  fPID2 = 0;
  SetRandomPosition();
  SetColWithResidNuclSwitch();
  SetStrongInterSwitch();
  SetQuantumStatistics();
  SetColoumb();
  SetTest();
  
}


AliHBTLLWeights* AliHBTLLWeights::Instance()
{
  // Instantiates new object or returns a pointer to already exitsing one
  if (fgLLWeights) {
    return fgLLWeights;
  } else {
    fgLLWeights = new AliHBTLLWeights();
    return fgLLWeights;
  }
}
                                            

Double_t AliHBTLLWeights::GetWeight(const AliHBTPair* partpair)
{
  // Returns the weignt of the pair "partpair"
  AliHBTParticle *part1 = partpair->Particle1();
  AliHBTParticle *part2 = partpair->Particle2();

  if ( (part1 == 0x0) || (part2 == 0x0))
    {
      Error("GetWeight","Null particle pointer");
      return 0.0;
    }
  
   
  Double_t part1Momentum[]={part1->Px(),part1->Py(),part1->Pz()};
  Double_t part2Momentum[]={part2->Px(),part2->Py(),part2->Pz()};
  
  if ( (part1->Px() == part2->Px()) && 
       (part1->Py() == part2->Py()) && 
       (part1->Pz() == part2->Pz()) )
    {
      return 0.0;
    }
  
             
  if ((!fRandomPosition) && 
      (part1->Vx()  == part2->Vx()) && (part1->Vy()  == part2->Vy())
      && (part1->Vz()  == part2->Vz()) )
    {        
      return 0.0;
    }
  
  
  
  FSI_MOM.P1X=part1Momentum[0];
  FSI_MOM.P1Y=part1Momentum[1];
  FSI_MOM.P1Z=part1Momentum[2];
  
  FSI_MOM.P2X=part2Momentum[0];
  FSI_MOM.P2Y=part2Momentum[1];
  FSI_MOM.P2Z=part2Momentum[2];
  
  if (fRandomPosition){
    
    Double_t rxcm = fsigma*gRandom->Gaus();
    Double_t rycm = fsigma*gRandom->Gaus();
    Double_t rzcm = fsigma*gRandom->Gaus();   
    
    FSI_PRF.X=rxcm*fwcons;
    FSI_PRF.Y=rycm*fwcons;
    FSI_PRF.Z=rzcm*fwcons;
    FSI_PRF.T=0.;
    
    Double_t rps=rxcm*rxcm+rycm*rycm+rzcm*rzcm;
    Double_t rp=TMath::Sqrt(rps);
    FSI_PRF.RP=rp;
    FSI_PRF.RPS=rps;
    
  }        
  
  ltran12();
  fsiw();

  if(flambda<1){
    if(gRandom->Rndm()<(1-flambda))LEDWEIGHT.WEIN=1.;}
  
  return LEDWEIGHT.WEIN;
}

/************************************************************/
void AliHBTLLWeights::Init()
{
  //---------------------------------------------------------------------  
  //initial parameters of model
  
  FSI_NS.NS = fApproximationModel;      
  
  if(!ftest){LEDWEIGHT.ITEST=0;}
  
  if(ftest){
    LEDWEIGHT.ITEST=1;
    if(fColoumbSwitch){FSI_NS.ICH =1;}
    else{FSI_NS.ICH=0;}
    if(fStrongInterSwitch){FSI_NS.ISI=1;}
    else{FSI_NS.ISI=0;}
    if(fQuantStatSwitch){FSI_NS.IQS=1;}
    else{FSI_NS.IQS=0;}
    if(fColWithResidNuclSwitch){FSI_NS.I3C=1;}
    else{FSI_NS.I3C=0;}
  }
  
  if(fRandomPosition){LEDWEIGHT.IRANPOS=1;}
  else{LEDWEIGHT.IRANPOS=0;}
  
  
  if ( (fPID1 == 0) || (fPID2 == 0) )
    {
      Fatal("Init","Particles types are not set");
      return;//pro forma
    }
  
  FSI_NS.LL = GetPairCode(fPID1,fPID2);
  
  if (FSI_NS.LL == 0) 
    {
      Fatal("Init","Particles types are not supported");
      return;//pro forma
    }
  
  
  TParticlePDG* tpart1 = TDatabasePDG::Instance()->GetParticle(fPID1);
  if (tpart1 == 0x0)
    {
      Fatal("init","We can not find particle with ID=%d in our DataBase",fPID1);
      return;
    }
  
  FSI_POC.AM1=tpart1->Mass();
  FSI_POC.C1=tpart1->Charge(); 
  
  TParticlePDG* tpart2 = TDatabasePDG::Instance()->GetParticle(fPID2);
  //mlv
  
  
  
  if (tpart2 == 0x0)
    {
      Fatal("init","We can not find particle with ID=%d in our DataBase",fPID2);
      return;
    }
  
  FSI_POC.AM2=tpart2->Mass();
  FSI_POC.C1=tpart2->Charge();
  
  led_bldata();
  fsiini();
  
  
  //constants for radii simulation 
  
  if(fRandomPosition){
    fsigma =TMath::Sqrt(2.)*fRadius;     
    fwcons =FSI_CONS.W;
  } 
} 

Int_t AliHBTLLWeights::GetPairCode(const AliHBTPair* partpair)
{
  // Return the code of the pair "partpair"
  return GetPairCode(partpair->Particle1()->GetPdgCode(),partpair->Particle2()->GetPdgCode());
}

Int_t AliHBTLLWeights::GetPairCode(Int_t pid1,Int_t pid2)
{
  //   pairCode   1  2  3  4   5    6   7  8  9 10  11  12  13  14 15 16 17 18  19  20   21   22  23 24 25 26    27     28
  //   hpid:      n  p  n alfa pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-  d  d    t   t    K0  K0  d p  p      p      n
  //   lpid:      n  p  p alfa pi- pi0 pi+ d  d  K-  K+  p   p  K- K+ p  p   d alfa  t  alfa  K0  K0b t t alfa lambda lambda
  //   NS=1 y/n:  +  +  +  +   +    -   -  -  -  -   -   -   -  -  -  -  -   -  -    -    -    -  -   - -  -      -      -
  
  //alphas, deuterons and tyts are NOT supported here
  
  Int_t chargefactor = 1;
  Int_t hpid; //pid in higher row
  Int_t lpid; //pid in lower row
  Int_t code; //pairCode
  
  Bool_t swap;
  
  //determine the order of selcetion in switch  
  if (TMath::Abs(pid1) < TMath::Abs(pid2) ) 
    {
      if (pid1<0) chargefactor=-1;
      hpid=pid2*chargefactor;
      lpid=pid1*chargefactor;
      swap = kFALSE;
    } 
  else 
    {
      if (pid2<0) chargefactor=-1;
      hpid=pid1*chargefactor;
      lpid=pid2*chargefactor;
      swap = kTRUE;
    }
  
  //mlv
  hpid=pid1;
  lpid=pid2;
  
  
  //Determine the pair code
  switch (hpid) //switch on first  particle id
    {
    case kNeutron:
      switch (lpid)
	{
	case kNeutron: 
	  code = 1;  //neutron neutron
	  break;
	  
	case kProton: 
	  code = 3;  //neutron proton
	  break;
	  
	case kLambda0: 
	  code = 28;  //neutron lambda
	  break;
	  
	default: 
	  return 0; //given pair not supported
	  break;
	}
      break;
      
    case kProton:
      switch (lpid)
	{
	case kProton:
	  code = 2; //proton proton
	  break;
	  
	case kLambda0: 
	  code = 27;//proton lambda
	  break;
	  
	default: 
	  return 0; //given pair not supported
	  break;
	  
	}
      break;
      
    case kPiPlus:
      
      switch (lpid)
	{
	case kPiPlus:
	  code = 7; //piplus piplus
	  break;
	  
	case kPiMinus:
	  code = 5; //piplus piminus
	  break;
	  
	case kKMinus:
	  code = 10; //piplus Kminus
	  break;
	  
	case kKPlus:
	  code = 11; //piplus Kplus
	  break;
	  
	case kProton:
	  code = 12; //piplus proton
	  chargefactor*=-1;
	  break;
	  
	default: 
	  return 0; //given pair not supported
	  break;
	}
      break;
    case kPi0:
      switch (lpid)
	{
	case kPi0:
	  code = 6;
	  break;
	  
	default: 
	  return 0; //given pair not supported
	  break;
	}
      break;
      
    case kKPlus:
      switch (lpid)
	{
	case kKMinus:
	  code = 14; //Kplus Kminus
	  break;
	  
	case kKPlus:
	  code = 15; //Kplus Kplus
	  break;
	  
	case kProton:
	  code = 16; //Kplus proton
	  break;
	  
	default: 
	  return 0; //given pair not supported
	  break;
	}
      break;
      
    case kKMinus:
      switch (lpid)
	{
	case kProton:
	  code = 17; //Kminus proton
	  chargefactor*=1;
	  break;
	  
	default: 
	  return 0; //given pair not supported
	  break;
	}
      break;
      
    case kK0:
      switch (lpid)
	{
	case kK0:
	  code = 2; //Kzero Kzero
	  break;
	  
	case kK0Bar:
	  code = 17; //Kzero KzeroBar
	  break;
	  
	default: 
	  return 0; //given pair not supported
	  break;
	}
      break;
      
    default: return 0;
    }
  return code;
}

