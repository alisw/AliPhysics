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

/* $Id$ */

//////////////////////////////////////////////////////////////////////
//                                                                  //         //
//      Generator of spectator nucleons (either protons or neutrons)//
//        computes beam crossing and divergence and Fermi momentum  //
//                                                                  //
/////////////////////////////////////////////////////////////////////

#include <assert.h>

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TMCProcess.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TVector3.h>

#include "AliConst.h"
#include "AliGenZDC.h"
#include "AliRun.h"
#include "AliMC.h"
 
ClassImp(AliGenZDC)
 
//_____________________________________________________________________________
AliGenZDC::AliGenZDC()
   :AliGenerator(),
  fIpart(0),
  fCosx(0),  	
  fCosy(0),  	
  fCosz(0),  	
  fPseudoRapidity(0),		
  fFermiflag(0),	
  fBeamDiv(0),	
  fBeamCrossAngle(0),
  fBeamCrossPlane(0),
  fDebugOpt(0)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliGenZDC::AliGenZDC(Int_t npart)
   :AliGenerator(npart),
  fIpart(kNeutron),
  fCosx(0.),  	
  fCosy(0.),  	
  fCosz(1.),  	
  fPseudoRapidity(0.),		
  fFermiflag(1),	
  fBeamDiv(0.000032),	
  fBeamCrossAngle(0.0001),
  fBeamCrossPlane(2),
  fDebugOpt(0)
{
  //
  // Standard constructor
  //
  fName = "AliGenZDC";
  fTitle = "Generation of Test Particles for ZDCs";
  
  for(Int_t i=0; i<201; i++){
     fProbintp[i] = 0;
     fProbintn[i] = 0;
     fPp[i] = 0;
  }
}

//_____________________________________________________________________________
void AliGenZDC::Init()
{
  //Initialize Fermi momentum distributions for Pb-Pb
  //
  printf("\n\n		AliGenZDC initialization:\n");
  printf("   Particle: %d, Track cosines: x = %f, y = %f, z = %f \n", 
  	 fIpart,fCosx,fCosy,fCosz);
  printf("   Fermi flag = %d, Beam divergence = %f, Crossing angle "
         "= %f, Crossing plane = %d\n\n", fFermiflag, fBeamDiv, fBeamCrossAngle,
	 fBeamCrossPlane);

  FermiTwoGaussian(208.);
}  
  
//_____________________________________________________________________________
void AliGenZDC::Generate()
{
  //
  // Generate one trigger (n or p)
  //
  Int_t i;

  Double_t mass, pLab[3], fP0, fP[3], fBoostP[3], ddp[3], dddp0, dddp[3]; 
  Float_t  fPTrack[3], ptot = fPMin;
  Int_t nt;
  
  if(fPseudoRapidity==0.){ 
    pLab[0] = ptot*fCosx;
    pLab[1] = ptot*fCosy;
    pLab[2] = ptot*fCosz;
  }
  else{
    Float_t scang = 2*TMath::ATan(TMath::Exp(-(fPseudoRapidity)));
    pLab[0] = -ptot*TMath::Sin(scang);
    pLab[1] = 0.;
    pLab[2] = ptot*TMath::Cos(scang);
  }
  for(i=0; i<=2; i++) fP[i] = pLab[i];  
  if(fDebugOpt == 1){
    printf("\n\n		Particle momentum before divergence and crossing\n");
    for(i=0; i<=2; i++)printf(" 	pLab[%d] = %f\n",i,pLab[i]);
  }
  
  // Beam divergence and crossing angle
  if(fBeamCrossAngle!=0.) {
    BeamDivCross(1, pLab);
    for(i=0; i<=2; i++) fP[i] = pLab[i];
  }
  if(fBeamDiv!=0.) {
    BeamDivCross(0, pLab);
    for(i=0; i<=2; i++) fP[i] = pLab[i];
  }

  // If required apply the Fermi momentum
  if(fFermiflag==1){
    if((fIpart==kProton) || (fIpart==kNeutron))
      ExtractFermi(fIpart, ddp);
    mass=TDatabasePDG::Instance()->GetParticle(fIpart)->Mass();
    fP0 = TMath::Sqrt(fP[0]*fP[0]+fP[1]*fP[1]+fP[2]*fP[2]+mass*mass);
    for(i=0; i<=2; i++) dddp[i] = ddp[i];
    dddp0 = TMath::Sqrt(dddp[0]*dddp[0]+dddp[1]*dddp[1]+dddp[2]*dddp[2]+mass*mass);
    
    TVector3 b(fP[0]/fP0, fP[1]/fP0, fP[2]/fP0);
    TLorentzVector pFermi(dddp[0], dddp[1], dddp[2], dddp0);

    pFermi.Boost(b);
    for(i=0; i<=2; i++){
       fBoostP[i] = pFermi[i];
       fP[i] = pFermi[i];
    }

  }
  
  for(i=0; i<=2; i++) fPTrack[i] = fP[i];
      
  Float_t polar[3] = {0,0,0};
  gAlice->GetMCApp()->PushTrack(fTrackIt,-1,fIpart,fPTrack,fOrigin.GetArray(),polar,0,
  		   kPPrimary,nt);
  // -----------------------------------------------------------------------
  if(fDebugOpt == 1){
    printf("\n\n		Track momentum:\n");
    printf("\n	 fPTrack = %f, %f, %f \n",fPTrack[0],fPTrack[1],fPTrack[2]);
  }
  else if(fDebugOpt == 2){
    FILE *file;
    if((file = fopen("SpectMomentum.dat","a")) == NULL){
      printf("Cannot open file  SpectMomentum.dat\n");
      return;
    }
    fprintf(file," %f \t %f \t %f \n",fPTrack[0],fPTrack[1],fPTrack[2]);
    fclose(file);
  }
    
}

//_____________________________________________________________________________
void AliGenZDC::FermiTwoGaussian(Float_t A)
{
//
// Momenta distributions according to the "double-gaussian"
// distribution (Ilinov) - equal for protons and neutrons
//

   Double_t sig1 = 0.113;
   Double_t sig2 = 0.250;
   Double_t alfa = 0.18*(TMath::Power((A/12.),(Float_t)1/3));
   Double_t xk = (2*k2PI)/((1.+alfa)*(TMath::Power(k2PI,1.5)));
   
   for(Int_t i=1; i<=200; i++){
      Double_t p = i*0.005;
      fPp[i] = p;
      Double_t e1 = (p*p)/(2.*sig1*sig1);
      Double_t e2 = (p*p)/(2.*sig2*sig2);
      Double_t f1 = TMath::Exp(-(e1));
      Double_t f2 = TMath::Exp(-(e2));
      Double_t probp = xk*p*p*(f1/(TMath::Power(sig1,3.))+
                      alfa*f2/(TMath::Power(sig2,3.)))*0.005;
      fProbintp[i] = fProbintp[i-1] + probp;
      fProbintn[i] = fProbintp[i];
   }
   if(fDebugOpt == 1){
     printf("\n\n		Initialization of Fermi momenta distribution \n");
     //for(Int_t i=0; i<=200; i++)
     //   printf(" fProbintp[%d] = %f, fProbintn[%d] = %f\n",i,fProbintp[i],i,fProbintn[i]);
   }
} 
//_____________________________________________________________________________
void AliGenZDC::ExtractFermi(Int_t id, Double_t *ddp)
{
//
// Compute Fermi momentum for spectator nucleons
//
  
  Int_t i;
  Float_t xx = gRandom->Rndm();
  assert ( id==kProton || id==kNeutron );
  if(id==kProton){
    for(i=1; i<=200; i++){
       if((xx>=fProbintp[i-1]) && (xx<fProbintp[i])) break;
       }
  }
  else {
    for(i=0; i<=200; i++){
       if((xx>=fProbintn[i-1]) && (xx<fProbintn[i])) break;
       }
   }
         Float_t pext = fPp[i]+0.001;
	 Float_t phi = k2PI*(gRandom->Rndm());
	 Float_t cost = (1.-2.*(gRandom->Rndm()));
	 Float_t tet = TMath::ACos(cost);
	 ddp[0] = pext*TMath::Sin(tet)*TMath::Cos(phi);
	 ddp[1] = pext*TMath::Sin(tet)*TMath::Sin(phi);
	 ddp[2] = pext*cost;

  if(fDebugOpt == 1){
    printf("\n\n		Extraction of Fermi momentum\n");
    printf("\n 	pxFermi = %f  pyFermi = %f  pzFermi = %f \n",ddp[0],ddp[1],ddp[2]); 
  }
}

//_____________________________________________________________________________
void AliGenZDC::BeamDivCross(Int_t icross, Double_t *pLab)
{
  // Applying beam divergence and crossing angle
  //
  Double_t tetpart, fipart, tetdiv=0, fidiv=0, angleSum[2], tetsum, fisum;
  Double_t rvec;

  Double_t pmq = 0.;
  Int_t i;
  for(i=0; i<=2; i++) pmq = pmq+pLab[i]*pLab[i];
  Double_t pmod = TMath::Sqrt(pmq);

  if(icross==0){      // ##### Beam divergence
    rvec = gRandom->Gaus(0.0,1.0);
    tetdiv = fBeamDiv * TMath::Abs(rvec);
    fidiv = (gRandom->Rndm())*k2PI;
  }
  else if(icross==1){ // ##### Crossing angle
    if(fBeamCrossPlane==0){
      tetdiv = 0.;
      fidiv = 0.;
    }
    else if(fBeamCrossPlane==1){     // Horizontal crossing plane
      tetdiv = fBeamCrossAngle;
      fidiv = 0.;
    }
    else if(fBeamCrossPlane==2){     // Vertical crossing plane
      tetdiv = fBeamCrossAngle;
      fidiv = k2PI/4.;
    }
  }

  tetpart = TMath::ATan2(TMath::Sqrt(pLab[0]*pLab[0]+pLab[1]*pLab[1]),pLab[2]);
  if(pLab[1]!=0. || pLab[0]!=0.) fipart = TMath::ATan2(pLab[1],pLab[0]);
  else fipart = 0.;
  if(fipart<0.) {fipart = fipart+k2PI;}
  tetdiv = tetdiv*kRaddeg;
  fidiv = fidiv*kRaddeg;
  tetpart = tetpart*kRaddeg;
  fipart = fipart*kRaddeg;
  AddAngle(tetpart,fipart,tetdiv,fidiv,angleSum);
  tetsum = angleSum[0];
  fisum  = angleSum[1];
  tetsum = tetsum*kDegrad;
  fisum = fisum*kDegrad;
  pLab[0] = pmod*TMath::Sin(tetsum)*TMath::Cos(fisum);
  pLab[1] = pmod*TMath::Sin(tetsum)*TMath::Sin(fisum);
  pLab[2] = pmod*TMath::Cos(tetsum);
  if(fDebugOpt == 1){
    if(icross==0) printf("\n\n		Beam divergence \n");
    else  	  printf("\n\n		Beam crossing \n");
    for(i=0; i<=2; i++)printf(" 	pLab[%d] = %f\n",i,pLab[i]);
  }
}
  
//_____________________________________________________________________________
void  AliGenZDC::AddAngle(Double_t theta1, Double_t phi1, Double_t theta2,
  	       Double_t phi2, Double_t *angleSum)
{ 
  // Calculating the sum of 2 angles
  Double_t temp, conv, cx, cy, cz, ct1, st1, ct2, st2, cp1, sp1, cp2, sp2;
  Double_t rtetsum, tetsum, fisum;
  
  temp = -1.;
  conv = 180./TMath::ACos(temp);
  
  ct1 = TMath::Cos(theta1/conv);
  st1 = TMath::Sin(theta1/conv);
  cp1 = TMath::Cos(phi1/conv);
  sp1 = TMath::Sin(phi1/conv);
  ct2 = TMath::Cos(theta2/conv);
  st2 = TMath::Sin(theta2/conv);
  cp2 = TMath::Cos(phi2/conv);
  sp2 = TMath::Sin(phi2/conv);
  cx = ct1*cp1*st2*cp2+st1*cp1*ct2-sp1*st2*sp2;
  cy = ct1*sp1*st2*cp2+st1*sp1*ct2+cp1*st2*sp2;
  cz = ct1*ct2-st1*st2*cp2;
  
  rtetsum = TMath::ACos(cz);
  tetsum = conv*rtetsum;
  if(tetsum==0. || tetsum==180.){
    fisum = 0.;
    return;
  }
  temp = cx/TMath::Sin(rtetsum);
  if(temp>1.) temp=1.;
  if(temp<-1.) temp=-1.;
  fisum = conv*TMath::ACos(temp);
  if(cy<0) {fisum = 360.-fisum;}
  angleSum[0] = tetsum;
  angleSum[1] = fisum;
}  

