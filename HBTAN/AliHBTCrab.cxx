#include "AliHBTCrab.h"
//______________________________________________________________________
/////////////////////////////////////////////////////////////////////////
//                                                                     //
// AliRoot wrapper to CRAB                                             //
// taken from http://www.nscl.msu.edu/~pratt/freecodes/crab/home.html  //
// written by Scott Pratt                                              //
//                                                                     //
//                                                                     //
/////////////////////////////////////////////////////////////////////////
 
#include "AliHBTPair.h"
#include "AliVAODParticle.h"
#include "TDatabasePDG.h"
#include <TMath.h>

#include "volya_complex.h"

AliHBTCrab* AliHBTCrab::fgCrab = 0x0;

const Double_t AliHBTCrab::fgkWcons = 1./0.1973;
const Double_t AliHBTCrab::fgkROOT2=1.41421356237309504880;
#ifdef __DECCXX
const complex AliHBTCrab::fgkCI(0.0,1.0);
#else
const doublecomplex AliHBTCrab::fgkCI(0.0,1.0);
#endif

/************************************************************/

AliHBTCrab* AliHBTCrab::Instance()
{
// returns instance of class
 if (fgCrab == 0x0)
  {
    fgCrab = new AliHBTCrab();
  }
 return fgCrab;
}
//===================================================================

void AliHBTCrab::Set()
{
//sets this as weighitng class
 Info("Set","Setting CRAB as Weighing Class");
 
 if ( fgWeights == 0x0 )  
  {
    fgWeights = AliHBTCrab::Instance();
    return;
  }  
 if ( fgWeights == AliHBTCrab::Instance() ) return;
 delete fgWeights;
 fgWeights = AliHBTCrab::Instance();
}
//===================================================================

AliHBTCrab::AliHBTCrab():
fBreitWigner(kFALSE),
fReducedMom(kTRUE),
fMaxMomentum(100.0)
{
  //ctor
  if(fgCrab)
   {
     Fatal("AliHBTCrab","Do not use constructor directly. Use Instance() instead.");
   }
}
//===================================================================

AliHBTCrab::AliHBTCrab(const AliHBTCrab &/*source*/):
AliHBTWeights(),
fBreitWigner(kFALSE),
fReducedMom(kTRUE),
fMaxMomentum(100.0)
{
  //ctor
}
//===================================================================
AliHBTCrab & AliHBTCrab::operator=(const AliHBTCrab& /*source*/)
{
//cpy constructor
 return *AliHBTCrab::Instance();
}

void AliHBTCrab::Init(Int_t pid1,Int_t pid2)
{
//Initialization method
  fMass1 = TDatabasePDG::Instance()->GetParticle(pid1)->Mass();
  fMass2 = TDatabasePDG::Instance()->GetParticle(pid2)->Mass();
  fInteractionWsym = 1.0;
  fInteractionWanti = 0.0;
  fInteractionWnosym = 0.0;
  fInteractionDelk  = 1.0;
  fInteractionNkmax = 100;
     
  fPid1 = pid1;   
  fPid2 = pid2;
}
//===================================================================

Bool_t AliHBTCrab::SetConfig(const AliHBTPair* pair)
{
//returns the SetConfig
  
  Int_t pdg1 = pair->Particle1()->GetPdgCode();
  Int_t pdg2 = pair->Particle2()->GetPdgCode();
  
  if ( ( pdg1 == fPid1) && ( pdg2  == fPid2) ) return kFALSE;
  else Init (pdg1,pdg2);
  
  return kTRUE;
}
//===================================================================

Double_t AliHBTCrab::GetWeight(AliHBTPair* partpair)
{
//returns the weight
  Double_t qred, r, qdotr, mom;
  Int_t test;
  
  SetConfig(partpair);
  
  GetComQuantities(partpair, &qred, &r, &qdotr, &mom, &test);
  if(test==0) 
   {
     Info("GetWeight","Test is 0");
   }
  Double_t corr = CorrCalc(qred,qdotr,r);
  
  return corr;
}
//===================================================================

void AliHBTCrab::GetComQuantities(const AliHBTPair* pair, 
       double *qred,double *r,double *qdotr,double *mom, int *test)
{
//************************************
//  ALICE //

 double p1[4]; 
 double p2[4]; 
 double r1[4]; 
 double r2[4]; 
 
 static const Double_t kCmToFm = 1.e13;
// static const Double_t cmtoOneOverGeV = kCmToFm*fgkWcons;
 
 AliVAODParticle *part1 = pair->Particle1();
 AliVAODParticle *part2 = pair->Particle2();

 p1[0] = part1->E()*1000.0;
 p1[1] = part1->Px()*1000.0;
 p1[2] = part1->Py()*1000.0;
 p1[3] = part1->Pz()*1000.0;
 
 p2[0] = part2->E()*1000.0;
 p2[1] = part2->Px()*1000.0;
 p2[2] = part2->Py()*1000.0;
 p2[3] = part2->Pz()*1000.0;

 r1[0] = part1->T();
 r1[1] = part1->Vx()*kCmToFm;
 r1[2] = part1->Vy()*kCmToFm;
 r1[3] = part1->Vz()*kCmToFm;

 r2[0] = part2->T();
 r2[1] = part2->Vx()*kCmToFm;
 r2[2] = part2->Vy()*kCmToFm;
 r2[3] = part2->Vz()*kCmToFm;

//  END OF ALICE STUFF

// This code is written by Scott Pratt
// taken from http://www.nscl.msu.edu/~pratt/freecodes/crab/home.html
  int alpha;
  double kdotr;
  double momtest;
  if (fReducedMom)
   {
     momtest=4.0*fMaxMomentum*fMaxMomentum;
   }  
  else
   {
     momtest=fMaxMomentum*fMaxMomentum;
   }  
   
  double ptot2,pdotr,pp,rr;
  
  if ( part1->GetPdgCode() == part2->GetPdgCode() )
   {
    *test=1;
    *mom=-(p2[0]-p1[0])*(p2[0]-p1[0]);
    for(alpha=1;alpha<4;alpha++){
      *mom=*mom+(p2[alpha]-p1[alpha])*(p2[alpha]-p1[alpha]);
    }
   //#if ! defined MIXED_PAIRS_FOR_DENOM
   //  if(*mom>momtest){
   //    *test=0;
   //    return;
   //  }
   //#endif
    pp=(p1[0]+p2[0]);
    rr=(r2[0]-r1[0]);
    pdotr=pp*rr;
    kdotr=(p2[0]-p1[0])*rr;
    ptot2=pp*pp;
    *r=-rr*rr;
    for(alpha=1;alpha<4;alpha++){
      pp=(p1[alpha]+p2[alpha]);
      rr=(r2[alpha]-r1[alpha]);
      pdotr=pdotr-pp*rr;
      kdotr=kdotr-(p2[alpha]-p1[alpha])*rr;
      ptot2=ptot2-pp*pp;
      *r=*r+rr*rr;
    }
    *mom=sqrt(*mom);
    *qred=0.5**mom;

    if (fReducedMom)
     {
      *mom=*qred;
     } 

    *qdotr=0.5*kdotr;
    *r=sqrt(*r+pdotr*pdotr/ptot2);
   }  
  else //identical
  {
  //  const double  kdotp=fMass2*fMass2-fMass1*fMass1;
   const double  kdotp = part2->Mass()*part2->Mass()- part1->Mass()*part1->Mass();
   *test=1;
   *mom=-(p2[0]-p1[0])*(p2[0]-p1[0]);
   ptot2=(p1[0]+p2[0])*(p1[0]+p2[0]);
   for(alpha=1;alpha<4;alpha++){
     *mom=*mom+(p2[alpha]-p1[alpha])*(p2[alpha]-p1[alpha]);
     ptot2=ptot2-(p1[alpha]+p2[alpha])*(p1[alpha]+p2[alpha]);
   }
   *mom=*mom+kdotp*kdotp/ptot2;
  //#if ! defined MIXED_PAIRS_FOR_DENOM
  //  if(*mom>momtest){
  //    *test=0;
  //    return;
  //  }
  //#endif
   pp=(p1[0]+p2[0]);
   rr=(r2[0]-r1[0]);
   pdotr=pp*rr;
   kdotr=(p2[0]-p1[0])*rr;
   *r=-rr*rr;
   for(alpha=1;alpha<4;alpha++){
     pp=(p1[alpha]+p2[alpha]);
     rr=(r2[alpha]-r1[alpha]);
     pdotr=pdotr-pp*rr;
     kdotr=kdotr-(p2[alpha]-p1[alpha])*rr;
     *r=*r+rr*rr;
   }
   kdotr=(-kdotr+kdotp*pdotr/ptot2);
   *mom=sqrt(*mom);
   *qred=0.5**mom;

   if (fReducedMom)
    {
     *mom=*qred;
    } 

   *qdotr=0.5*kdotr;
   *r=sqrt(*r+pdotr*pdotr/ptot2);
  }//not identical  

  return;
}


//===================================================================

double  AliHBTCrab::CorrCalc(double trueqred,double trueqdotr,double truer)
{
//#define REDUCED_MOM
// This code is written by Scott Pratt
// taken from http://www.nscl.msu.edu/~pratt/freecodes/crab/home.html
  double eta,arg,corr0;
//  double xx,xxprime,xxjj,p1,zk;
//  int jj,kk,ipart,ipartcount,ispin;
  int kk;
  double wsymleftover,wantileftover,wnosymleftover;
  double qred,qdotr,r;
//  const double rmass=fMass1*fMass2/(fMass1+fMass2);
#ifdef __DECCXX
  complex cphi1,cphi2,cphis,cphia;
#else
  doublecomplex cphi1,cphi2,cphis,cphia;
#endif

  arg=trueqdotr/197.323-2.0*TMath::Pi()*TMath::Floor(trueqdotr/(197.323*2.0*TMath::Pi()));
  cphi1=exp(fgkCI*arg);
  cphis=fgkROOT2*real(cphi1);
  cphia=fgkCI*fgkROOT2*imag(cphi1);
  corr0=real(fInteractionWsym*cphis*conj(cphis)
	     +fInteractionWanti*cphia*conj(cphia)
	     +fInteractionWnosym*cphi1*conj(cphi1));
  goto OUTSIDE_INTERACTION_RANGE;

#ifdef REDUCED_MOM
  kk=(int)TMath::Floor(trueqred/fInteractionDelk);
  qred=(0.5+kk)*fInteractionDelk;
#else
  kk=(int)TMath::Floor(2.0*trueqred/fInteractionDelk);
  qred=(0.5+kk)*fInteractionDelk/2.0;
#endif
  qdotr=trueqdotr*qred/trueqred;
  if(kk>=fInteractionNkmax){
    corr0=1.0;
    goto OUTSIDE_INTERACTION_RANGE;
  }
  r=truer;

  eta=0.0;
  arg=qdotr/197.323-2.0*TMath::Pi()*TMath::Floor(qdotr/(197.323*2.0*TMath::Pi()));
  cphi1=exp(fgkCI*arg);
  cphi2=conj(cphi1);

  cphis=(cphi1+cphi2)/fgkROOT2;
  cphia=(cphi1-cphi2)/fgkROOT2;
  corr0=0.0;
  /* If there are corrections for strong interactions, add the
     change for each partial wave.  If npartial = 0 then there
     are no strong int. corrections. */
  wsymleftover=fInteractionWsym;
  wantileftover=fInteractionWanti;
  wnosymleftover=fInteractionWnosym;

  corr0=corr0+real(wsymleftover*cphis*conj(cphis)
		   +wantileftover*cphia*conj(cphia)
		   +wnosymleftover*cphi1*conj(cphi1));
OUTSIDE_INTERACTION_RANGE:

#ifdef BREIT_WIGNER
  corr0=corr0+bwcalc(trueqred,truer);
#endif

  return corr0;
}

#ifdef __DECCXX
complex AliHBTCrab::CGamma(complex c){
/* This calc.s gamma functions which are in the form gamma(n+i*y)
   where n is an int and y is real. */
// This code is written by Scott Pratt
// taken from http://www.nscl.msu.edu/~pratt/freecodes/crab/home.html
#else
doublecomplex AliHBTCrab::CGamma(doublecomplex c){
/* This calc.s gamma functions which are in the form gamma(n+i*y)
   where n is an int and y is real. */
// This code is written by Scott Pratt
// taken from http://www.nscl.msu.edu/~pratt/freecodes/crab/home.html
#endif
#ifdef __DECCXX
  complex cg,cphase;
#else
  doublecomplex cg,cphase;
#endif
  int mm,j;
  double x,y,phase,delp,cgmag;
  x=real(c);
  y=imag(c);
  phase=-TMath::E()*y;
  for(j=1;j<=100000;j++){
    delp=(y/(double)j)-atan(y/(double)j);
    phase=phase+delp;
    if(TMath::Abs(delp)<1E-10) goto CGamma_ESCAPE;
  }
  printf("oops not accurate enough, increase jmax\n");
CGamma_ESCAPE:
  phase=phase-2.0*TMath::Pi()*TMath::Floor(phase/(2.0*TMath::Pi()));
  cphase=exp(fgkCI*phase);
  cgmag=sqrt(TMath::Pi()*y/sinh(TMath::Pi()*y));
  mm=(int)TMath::Floor(x+0.5);
  cg=cgmag*cphase;
  if(mm<1){
    for(j=1;j<=-mm+1;j++){
      cg=cg/(1.0+(double)(-j)+fgkCI*y);
    }
  }
  if(mm>1) {
    for(j=1;j<=mm-1;j++){
      cg=cg*((double)(j)+fgkCI*y);
    }
  }
  return cg;
}
