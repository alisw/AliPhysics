#include "AliHBTLLWeights.h"
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//_________________________________________________________________________
///////////////////////////////////////////////////////////////////////////
//
//  class AliHBTLLWeights
//
//  This class introduces the weight's calculation 
//  according to the Lednicky's algorithm.
//  
//  
//  fsiw.f, fsiini.f  
//
//  Description from fortran code by author R. Lednicky
//
//  Calculates final state interaction (FSI) weights                     
//  WEIF = weight due to particle - (effective) nucleus FSI (p-N)
//  WEI  = weight due to p-p-N FSI
//  WEIN = weight due to p-p FSI; note that WEIN=WEI if I3C=0;
//                                note that if I3C=1 the calculation of
//                                WEIN can be skipped by putting J=0
//.......................................................................
//  Correlation Functions:
//  CF(p-p-N)   = sum(WEI)/sum(WEIF)
//  CF(p-p)     = sum(WEIN)/sum(1); here the nucleus is completely
//                                  inactive
//  CF(p-p-"N") = sum(WEIN*WEIF')/sum(WEIF'), where WEIN and WEIF'
//                are not correlated (calculated at different emission
//                points, e.g., for different events);
//                thus here the nucleus affects one-particle
//                spectra but not the correlation
//.......................................................................
//  User must supply data file <fn> on unit NUNIT (e.g. =11) specifying
//  LL   : particle pair
//  NS   : approximation used to calculate Bethe-Salpeter amplitude
//  ITEST: test switch
//         If ITEST=1 then also following parameters are required
//  ICH  : 1(0) Coulomb interaction between the two particles ON (OFF)
//  IQS  : 1(0) quantum statistics for the two particles ON (OFF)
//  ISI  : 1(0) strong interaction between the two particles ON (OFF)
//  I3C  : 1(0) Coulomb interaction with residual nucleus ON (OFF)
//  This data file can contain other information useful for the user.
//  It is read by subroutines READINT4 and READREA8(4) (or READ_FILE).
//  -------------------------------------------------------------------
//-   LL       1  2  3  4   5    6   7  8  9 10  11  12  13  14 15 16 17
//-   part. 1: n  p  n alfa pi+ pi0 pi+ n  p pi+ pi+ pi+ pi- K+ K+ K+ K-
//-   part. 2: n  p  p alfa pi- pi0 pi+ d  d  K-  K+  p   p  K- K+ p  p
//   NS=1 y/n: +  +  +  +   +    -   -  -  -  -   -   -   -  -  -  -  -
//  -------------------------------------------------------------------
//-   LL       18  19 20  21  22 23 24 25 26    27     28
//-   part. 1:  d  d   t  t   K0 K0  d p  p      p      n
//-   part. 2:  d alfa t alfa K0 K0b t t alfa lambda lambda
//   NS=1 y/n:  -  -   -  -   -  -   - -  -      +      +
//  -------------------------------------------------------------------
//   NS=1  Square well potential,
//   NS=3  not used
//   NS=4  scattered wave approximated by the spherical wave,
//   NS=2  same as NS=4 but the approx. of equal emission times in PRF
//       not required (t=0 approx. used in all other cases).
//   Note: if NS=2,4, the B-S amplitude diverges at zero distance r* in
//       the two-particle c.m.s.; user can specify a cutoff AA in
//       SUBROUTINE FSIINI, for example:
//       IF(NS.EQ.2.OR.NS.EQ.4)AA=5.D0 !! in 1/GeV --> AA=1. fm
//  ------------------------------------------------------------------
//  ITEST=1 any values of parameters ICH, IQS, ISI, I3C are allowed
//          and should be given in data file <fn>
//  ITEST=0 physical values of these parameters are put automatically
//          in FSIINI (their values are not required in data file)
//=====================================================================
//  At the beginning of calculation user should call FSIINI,
//  which reads LL, NS, ITEST (and eventually ICH, IQS, ISI, I3C)
//  and ializes various parameters.
//  In particular the constants in
//    COMMON/FSI_CONS/PI,PI2,SPI,DR,W
//  may be useful for the user:
//   W=1/.1973D0    ! from fm to 1/GeV
//   PI=4*DATAN(1.D0)
//   PI2=2*PI
//   SPI=TMath::Sqrt(PI)
//   DR=180.D0/PI   ! from radian to degree
//    _______________________________________________________
//  !! |Important note: all real quantities are assumed REAL*8 | !!
//    -------------------------------------------------------
//  For each event user should fill in the following information
//  in COMMONs (all COMMONs in FSI calculation start with FSI_):
//  ...................................................................
//   COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
//  Only
//       AMN  = mass of the effective nucleus   [GeV/c**2]
//       CN   = charge of the effective nucleus [elem. charge units]
//  are required
//  ...................................................................
//   COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1, !part. momenta in the rest frame 
//  1               P2X,P2Y,P2Z,E2,P2  !of effective nucleus (NRF)
//  Only the components
//                      PiX,PiY,PiZ  [GeV/c]
//  in NRF are required.
//  To make the corresponding Lorentz transformation user can use the
//  subroutines LTRAN and LTRANB
//  ...................................................................
//  COMMON/FSI_COOR/X1,Y1,Z1,T1,R1,     ! 4-coord. of emission
//  1               X2,Y2,Z2,T2,R2      ! points in NRF
//  The componets
//                     Xi,Yi,Zi  [fm]
//  and emission times
//                        Ti   [fm/c]
//  should be given in NRF with the origin assumed at the center
//  of the effective nucleus. If the effect of residual nucleus is
//  not calculated within FSIW, the NRF can be any fixed frame.
//  --------------------------------------------------------------------
//  Before calling FSIW the user must call
//   CALL LTRAN12
//  Besides Lorentz transformation to pair rest frame:
//  (p1-p2)/2 --> k* it also transforms 4-coordinates of
//  emission points from fm to 1/GeV and calculates Ei,Pi and Ri.
//  Note that |k*|=AK in COMMON/FSI_PRF/
//  --------------------------------------------------------------------
//  After making some additional filtering using k* (say k* < k*max)
//  or direction of vector k*,
//  user can finally call FSIW to calculate the FSI weights
//  to be used to construct the correlation function
//======================================================================


/*******************************************************************/
/******      ROUTINES    USED    FOR     COMMUNUCATION      ********/
/********************     WITH      FORTRAN     ********************/
/*******************************************************************/
#ifndef WIN32
# define led_bldata led_bldata_
# define fsiini fsiini_
# define ltran12 ltran12_
# define boosttoprf boosttoprf_
# define fsiw fsiw_
# define setpdist setpdist_
# define type_of_call
#else
# define led_bldata LED_BLDATA
# define fsiini FSIINI
# define ltran12 LTRAN12
# define boosttoprf BOOSTTOPRF
# define fsiw FSIW
# define setpdist SETPDIST
# define type_of_call _stdcall
#endif
/****************************************************************/
extern "C" void type_of_call led_bldata(); 
extern "C" void type_of_call fsiini();
extern "C" void type_of_call ltran12();
extern "C" void type_of_call boosttoprf();
extern "C" void type_of_call fsiw();
extern "C" void type_of_call setpdist(Double_t& r);
/**************************************************************/

#include "AliHBTPair.h"
#include "AliVAODParticle.h"
#include "WLedCOMMONS.h"
#include <TRandom.h>   
#include <TMath.h>     
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>


ClassImp(AliHBTLLWeights)  
 
AliHBTLLWeights* AliHBTLLWeights::fgLLWeights = 0x0; 
const Double_t AliHBTLLWeights::fgkWcons = 1./0.1973;

AliHBTLLWeights::AliHBTLLWeights():
 fTest(kTRUE),
 fColoumbSwitch(kTRUE),
 fQuantStatSwitch(kTRUE),
 fStrongInterSwitch(kTRUE),
 fColWithResidNuclSwitch(kFALSE),
 fNuclMass(0.0),
 fNuclCharge(0.0),
 fRandomPosition(kNone),
 fRadius(0.0),
 fOneMinusLambda(0.0),
 fPID1(0),
 fPID2(0),
 fSigma(0.0)
{
// Default Constructor 
  if (fgLLWeights)
   Fatal("AliHBTLLWeights","LLWeights already instatiated. Use AliHBTLLWeights::Instance()");
}
/**************************************************************/

AliHBTLLWeights::AliHBTLLWeights(const AliHBTLLWeights &/*source*/):
 AliHBTWeights(),
 fTest(kTRUE),
 fColoumbSwitch(kTRUE),
 fQuantStatSwitch(kTRUE),
 fStrongInterSwitch(kTRUE),
 fColWithResidNuclSwitch(kFALSE),
 fNuclMass(0.0),
 fNuclCharge(0.0),
 fRandomPosition(kNone),
 fRadius(0.0),
 fOneMinusLambda(0.0),
 fPID1(0),
 fPID2(0),
 fSigma(0.0)
{
  //Copy ctor needed by the coding conventions but not used
  Fatal("AliHBTLLWeights","copy ctor not implemented");
}
/************************************************************/

AliHBTLLWeights& AliHBTLLWeights::operator=(const AliHBTLLWeights& /*source*/)
{
  //Assignment operator needed by the coding conventions but not used
  Fatal("AliHBTLLWeights","assignment operator not implemented");
  return * this;
}
/************************************************************/

AliHBTLLWeights* AliHBTLLWeights::Instance()
{     
// returns instance of class 
 if (fgLLWeights) 
  {
    return fgLLWeights;
  } 
 else 
  {
   fgLLWeights = new AliHBTLLWeights();            
   return fgLLWeights; 
  } 
}     
/************************************************************/

void AliHBTLLWeights::Set()
{
 //sets this as weighitng class
 Info("Set","Setting Lednicky-Lyuboshitz as Weighing Class");
 
 if ( fgWeights == 0x0 )  
  {
    fgWeights = AliHBTLLWeights::Instance();
    return;
  }  
 if ( fgWeights == AliHBTLLWeights::Instance() ) return;
 delete fgWeights;
 fgWeights = AliHBTLLWeights::Instance();
}
/************************************************************/

Double_t AliHBTLLWeights::GetWeight(AliHBTPair* partpair)
{
// calculates weight for a pair
  static const Double_t kcmtofm = 1.e13;
  static const Double_t kcmtoOneOverGeV = kcmtofm*fgkWcons;  
  
  AliVAODParticle *part1 = partpair->Particle1();
  AliVAODParticle *part2 = partpair->Particle2();

  if ( (part1 == 0x0) || (part2 == 0x0))
   {
     Error("GetWeight","Null particle pointer");
     return 0.0;
   }

  if ( fPID1 != part1->GetPdgCode() ) return 1.0; 
  if ( fPID2 != part2->GetPdgCode() ) return 1.0; 

//takes a lot of time
  if ( (part1->Px() == part2->Px()) && 
       (part1->Py() == part2->Py()) && 
       (part1->Pz() == part2->Pz()) )
   {
     return 0.0;
   }

  if ((!fRandomPosition) && 
      (part1->Vx()  == part2->Vx()) && 
      (part1->Vy()  == part2->Vy()) && 
      (part1->Vz()  == part2->Vz()) )
    {        
      return 0.0;
    }

  if(fOneMinusLambda)//implemetation of non-zero intetcept parameter 
   {
     if( gRandom->Rndm() < fOneMinusLambda ) return 1.0;
   }
    

  //this must  be after ltran12 because it would overwrite what we set below
  switch (fRandomPosition)
   {
    case kGaussianQInv:
     {
       // Set Momenta in the common block
       FSI_MOM.P1X = part1->Px();
       FSI_MOM.P1Y = part1->Py();
       FSI_MOM.P1Z = part1->Pz();

       FSI_MOM.P2X = part2->Px();
       FSI_MOM.P2Y = part2->Py();
       FSI_MOM.P2Z = part2->Pz();
      
       //boost it to PRF
       ltran12();
       boosttoprf();
       
       //Set particle positions now so they are Gaussian in PRF
       RandomPairDistances();
       
       FSI_PRF.RPS= FSI_COOR.X2*FSI_COOR.X2 + FSI_COOR.Y2*FSI_COOR.Y2 +FSI_COOR.Z2*FSI_COOR.Z2;
       FSI_PRF.RP=TMath::Sqrt(FSI_PRF.RPS);
       
       break;
      } 
    case kGaussianOSL:
      { 
       //boost pair to LCMS and set such a momenta in the common block
       Int_t retv = SetMomentaInLCMS(part1,part2);
       if (retv) return 1;
       //random particle positions/distance so they are Gaussian in LCMS
       RandomPairDistances();
       
       //Boost to PRF
       ltran12();
       boosttoprf();
    
       break;
      }
    case kNone:
    default:
       //set momenta and paricle positions as they are
       FSI_MOM.P1X = part1->Px();
       FSI_MOM.P1Y = part1->Py();
       FSI_MOM.P1Z = part1->Pz();

       FSI_MOM.P2X = part2->Px();
       FSI_MOM.P2Y = part2->Py();
       FSI_MOM.P2Z = part2->Pz();

       FSI_COOR.X1 = part1->Vx()*kcmtoOneOverGeV;
       FSI_COOR.Y1 = part1->Vy()*kcmtoOneOverGeV;
       FSI_COOR.Z1 = part1->Vz()*kcmtoOneOverGeV;
       FSI_COOR.T1 = part1->T()*kcmtoOneOverGeV;

       FSI_COOR.X2 = part2->Vx()*kcmtoOneOverGeV;
       FSI_COOR.Y2 = part2->Vy()*kcmtoOneOverGeV;
       FSI_COOR.Z2 = part2->Vz()*kcmtoOneOverGeV;
       FSI_COOR.T2 = part2->T()*kcmtoOneOverGeV;

       ltran12();
       boosttoprf();
       
       break;
   }

  fsiw();
  return LEDWEIGHT.WEIN;
}
/************************************************************/

void AliHBTLLWeights::Init()
{
//initial parameters of model

  FSI_NS.NS = fApproximationModel;      
  
  LEDWEIGHT.ITEST = fTest;  
  if(fTest)
   {
     FSI_NS.ICH = fColoumbSwitch;
     FSI_NS.ISI = fStrongInterSwitch;
     FSI_NS.IQS = fQuantStatSwitch;
     FSI_NS.I3C = fColWithResidNuclSwitch;
     LEDWEIGHT.IRANPOS = fRandomPosition;
   }
 
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

  Info("Init","Setting PIDs %d %d. LL Code is %d",fPID1,fPID2,FSI_NS.LL);


  TParticlePDG* tpart1 = TDatabasePDG::Instance()->GetParticle(fPID1);
  if (tpart1 == 0x0)
   {
     Fatal("init","We can not find particle with ID=%d in PDG DataBase",fPID1);
     return;
   }
      
  FSI_POC.AM1=tpart1->Mass();
  FSI_POC.C1=tpart1->Charge(); 

  TParticlePDG* tpart2 = TDatabasePDG::Instance()->GetParticle(fPID2);
//lv
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

  if(fRandomPosition)
   {
     fSigma =TMath::Sqrt(2.)*fRadius;     
   } 
} 
/************************************************************/

Int_t AliHBTLLWeights::GetPairCode(const AliHBTPair* partpair)
{
//returns Code corresponding to that pair
 return GetPairCode(partpair->Particle1()->GetPdgCode(),partpair->Particle2()->GetPdgCode());
}
/************************************************************/

Int_t AliHBTLLWeights::GetPairCode(Int_t pid1,Int_t pid2)
{
// returns code corresponding to the pair of PIDs
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
/************************************************************/

void AliHBTLLWeights::SetTest(Bool_t rtest)
{
  //Sets fTest member
  fTest = rtest;
} 
/************************************************************/

void AliHBTLLWeights::SetColoumb(Bool_t col)
{
  // (ICH in fortran code) Coulomb interaction between the two particles ON (OFF)
  fColoumbSwitch = col;
}
/************************************************************/

void AliHBTLLWeights::SetQuantumStatistics(Bool_t qss)
{
  //IQS: quantum statistics for the two particles ON (OFF) 
  //if non-identical particles automatically off
  fQuantStatSwitch = qss;
}
/************************************************************/

void AliHBTLLWeights::SetStrongInterSwitch(Bool_t sis)
{
  //ISI: strong interaction between the two particles ON (OFF)
  fStrongInterSwitch = sis;
}
/************************************************************/

void AliHBTLLWeights::SetColWithResidNuclSwitch(Bool_t crn)
{
  //I3C: Coulomb interaction with residual nucleus ON (OFF)  
  fColWithResidNuclSwitch = crn;
}
/************************************************************/

void AliHBTLLWeights::SetApproxModel(Int_t ap)
{
  //sets  Model of Approximation (NS in Fortran code)
  fApproximationModel=ap;
}
/************************************************************/
     
void AliHBTLLWeights::SetRandomPosition(ERandomizationWay rw)
{ 
// rw can be: kGaussianQInv - so 1D Qinv correlation function has a Gaussian shape 
//                            (and also Qlong and Qside, but not Qout)
//            kGaussianOSL - so 3D Qout-Qside-Qlong correlation function has a Gaussian shape
//            kNone  - no randomization performed (DEFAULT)
 fRandomPosition = rw;
}
/************************************************************/

void AliHBTLLWeights::SetR1dw(Double_t R)
{
  //spherical source model radii
  fRadius=R;
}
/************************************************************/

void AliHBTLLWeights::SetParticlesTypes(Int_t pid1, Int_t pid2)
{
  //set AliRoot particles types   
  fPID1 = pid1; 
  fPID2 = pid2;
}
/************************************************************/
    
void AliHBTLLWeights::SetNucleusCharge(Double_t ch)
{
  // not used now  (see comments in fortran code)
  fNuclCharge=ch;
}
/************************************************************/

void AliHBTLLWeights::SetNucleusMass(Double_t mass)
{
  // (see comments in fortran code)
  fNuclMass=mass;
}
/************************************************************/

Int_t AliHBTLLWeights::SetMomentaInLCMS(AliVAODParticle* part1, AliVAODParticle* part2)
{
//sets paricle momenta in the common block in LCMS frame
//---> Particle energies ---------

  Double_t p1x = part1->Px();
  Double_t p1y = part1->Py();
  Double_t p1z = part1->Pz();
      
  Double_t p2x = part2->Px();
  Double_t p2y = part2->Py();
  Double_t p2z = part2->Pz();

  Double_t am1 = part1->Mass();
  Double_t am2 = part2->Mass();
  
  Double_t p1s=p1x*p1x+p1y*p1y+p1z*p1z;
  Double_t p2s=p2x*p2x+p2y*p2y+p2z*p2z;
  Double_t e1=TMath::Sqrt(am1*am1+p1s);
  Double_t e2=TMath::Sqrt(am2*am2+p2s);
//---> pair parameters -----------
  Double_t e12=e1+e2;       // energy
  Double_t p12x=p1x+p2x;    // px
  Double_t p12y=p1y+p2y;    // py
  Double_t p12z=p1z+p2z;    // pz
  Double_t p12s=p12x*p12x+p12y*p12y+p12z*p12z;
  Double_t p12 =TMath::Sqrt(p12s);// momentum
  Double_t v12 =p12/e12;    // velocity

  Double_t cth =p12z/p12;   // cos(theta)
//  Double_t sth =TMath::Sqrt(1. - cth*cth); //sin
  Double_t v12z=v12*cth;    // longit. v
  Double_t gamz=1./TMath::Sqrt(1.-v12z*v12z);
//--      v12t=v12*sth    // transv. v in cms (not needed)


  Double_t p12ts=p12x*p12x+p12y*p12y;
  Double_t p12t=TMath::Sqrt(p12ts); //pt
//===> azimuthal rotation (pt||x) ============
   if(p12t != 0.0) 
    {
     Double_t cphi=p12x/p12t;  // cos(phi)
     Double_t sphi=p12y/p12t;  // sin(phi)
     Rotate(p1x,p1y,sphi,cphi,p1x,p1y);       
     Rotate(p2x,p2y,sphi,cphi,p2x,p2y);
//     Rotate(x1,y1,sphi,cphi,x1,y1);
//     Rotate(x2,y2,sphi,cphi,x2,y2);
     }    
   else // rotation impossible 
    {
      return 1;
    }

//===> co-moving ref. frame       ============
   
   Double_t nothing;
   Boost(p1z,e1,v12z,gamz,p1z,nothing);
   Boost(p2z,e2,v12z,gamz,p2z,nothing);

//   p1s=p1x*p1x+p1y*p1y+p1z*p1z;
//   p2s=p2x*p2x+p2y*p2y+p2z*p2z;
//   e1=TMath::Sqrt(am1*am1+p1s);
//   e2=TMath::Sqrt(am2*am2+p2s);
//   Boost(z1,t1,v12z,gamz,z1,t1);
//   Boost(z2,t2,v12z,gamz,z2,t2);


     FSI_MOM.P1X = p1x;
     FSI_MOM.P1Y = p1y;
     FSI_MOM.P1Z = p1z;

     FSI_MOM.P2X = p2x;
     FSI_MOM.P2Y = p2y;
     FSI_MOM.P2Z = p2z;
     
//     Info("Simulate"," %f %f %f ",p1x + p2x, p1y + p2y, p1z + p2z);
     return 0;
}

/**********************************************************/  

void AliHBTLLWeights::RandomPairDistances()
{
 //randomizes pair distances
  Double_t rxcm = fSigma*gRandom->Gaus();
  Double_t rycm = fSigma*gRandom->Gaus();
  Double_t rzcm = fSigma*gRandom->Gaus();


  FSI_COOR.X1 = 0;
  FSI_COOR.Y1 = 0;
  FSI_COOR.Z1 = 0;
  FSI_COOR.T1 = 0;

  FSI_COOR.X2 = rxcm*fgkWcons;
  FSI_COOR.Y2 = rycm*fgkWcons;
  FSI_COOR.Z2 = rzcm*fgkWcons;
  FSI_COOR.T2 = 0;

}


void AliHBTLLWeights::Rotate(Double_t x, Double_t y, Double_t sphi, Double_t cphi, Double_t& xr, Double_t& yr)
{
//rotates 
 yr=y*cphi-x*sphi;
 xr=x*cphi+y*sphi;
}

void AliHBTLLWeights::Boost(Double_t z, Double_t t, Double_t beta, Double_t gamma, Double_t& zt, Double_t& yt)
{
//boosts
  zt=gamma*(z-beta*t);
  yt=gamma*(t-beta*z);
}
