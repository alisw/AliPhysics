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

/*
$Log$
Revision 1.1  2000/07/10 13:58:01  fca
New version of ZDC from E.Scomparin & C.Oppedisano

Revision 1.7  2000/01/19 17:17:40  fca

Revision 1.6  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "AliGenZDC.h"
#include "AliConst.h"
#include "AliPDG.h"
#include "AliRun.h"
 
ClassImp(AliGenZDC)
 
//_____________________________________________________________________________
AliGenZDC::AliGenZDC()
   :AliGenerator()
{
  //
  // Default constructor
  //
  fIpart = 0;
}

//_____________________________________________________________________________
AliGenZDC::AliGenZDC(Int_t npart)
   :AliGenerator(npart)
{
  //
  // Standard constructor
  //
  fName = "AliGenZDC";
  fTitle = "Generation of Test Particles for ZDCs";
  fIpart = kNeutron;
  fCosx  = 0.;
  fCosy  = 0.;
  fCosz  = 1.;
  fPseudoRapidity = 0.;
  fFermiflag = 1;
  // LHC values for beam divergence and crossing angle
  fBeamDiv = 0.000032;
  fBeamCrossAngle = 0.0001;
  fBeamCrossPlane = 2;
}

//_____________________________________________________________________________
void AliGenZDC::Init()
{
  printf("		Initializing AliGenZDC\n");
  printf("	Fermi flag = %d, Beam Divergence = %f, Crossing Angle "
         "= %f, Crossing Plane = %d\n\n", fFermiflag, fBeamDiv, fBeamCrossAngle,
	 fBeamCrossPlane);
  //Initialize Fermi momentum distributions for Pb-Pb
  FermiTwoGaussian(207.,82.,fPp,fProbintp,fProbintn);
}  
  
//_____________________________________________________________________________
void AliGenZDC::Generate()
{
  //
  // Generate one trigger (n or p)
  //
  Int_t i;

  Double_t mass, pLab[3], balp0, balp[3], ddp[3], dddp0, dddp[3];
  Float_t ptot = fPMin;
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
  for(i=0; i<=2; i++){
     fP[i] = pLab[i];
  }
  
  // Beam divergence and crossing angle
  if(fBeamDiv!=0.) {BeamDivCross(0,fBeamDiv,fBeamCrossAngle,fBeamCrossPlane,pLab);}
  if(fBeamCrossAngle!=0.) {BeamDivCross(1,fBeamDiv,fBeamCrossAngle,fBeamCrossPlane,pLab);}
  
  // If required apply the Fermi momentum
  if(fFermiflag==1){
    if((fIpart==kProton) || (fIpart==kNeutron)){
      ExtractFermi(fIpart,fPp,fProbintp,fProbintn,ddp);
    }
    if(fIpart==kProton) {mass = 0.93956563;}
    if(fIpart==kNeutron) {mass = 0.93827231;}
//  printf(" pLABx = %f  pLABy = %f  pLABz = %f \n",pLab[0],pLab[1],pLab[2]); 
    for(i=0; i<=2; i++){
       balp[i] = -pLab[i];
    }
    balp0 = TMath::Sqrt(pLab[0]*pLab[0]+pLab[1]*pLab[1]+pLab[2]*pLab[2]+mass*mass);
    for(i=0; i<=2; i++){
       dddp[i] = ddp[i];
    }
    dddp0 = TMath::Sqrt(dddp[0]*dddp[0]+dddp[1]*dddp[1]+dddp[2]*dddp[2]+mass*mass);
    
    TVector3 b(balp[0]/balp0, balp[1]/balp0, balp[2]/balp0);
    TLorentzVector pFermi(dddp[0], dddp[1], dddp[2], dddp0);

//    printf(" pmu -> pLABx = %f  pLABy = %f  pLABz = %f  E = %f\n",
//           balp[0],balp[1],balp[2],balp0); 
//    printf(" Beta -> bx = %f  by = %f  bz = %f\n", b[0], b[1], b[2]);  
//    printf(" pFermi -> px = %f, py = %f, pz = %f\n", pFermi[0], pFermi[1], pFermi[2]);
    
    pFermi.Boost(b);

//    printf(" Boosted momentum -> px = %f, py = %f, pz = %f\n",
//	     pFermi[0], pFermi[1], pFermi[2]);
    for(i=0; i<=2; i++){
       fBoostP[i] = pFermi[i];
    }

  }
    
  Float_t polar[3] = {0,0,0};
  gAlice->SetTrack(fTrackIt,-1,fIpart,fBoostP,fOrigin.GetArray(),polar,0,
  		   "Primary",nt);
}

//_____________________________________________________________________________
void AliGenZDC::FermiTwoGaussian(Double_t A, Float_t Z, Double_t* fPp, Double_t*
  		fProbintp, Double_t* fProbintn)
{
//
// Momenta distributions according to the "double-gaussian"
// distribution (Ilinov) - equal for protons and neutrons
//
//   printf("		Initialization of Fermi momenta distribution\n");
   fProbintp[0] = 0;
   fProbintn[0] = 0;
   Double_t sig1 = 0.113;
   Double_t sig2 = 0.250;
   Double_t alfa = 0.18*(TMath::Power((A/12.),(Float_t)1/3));
   Double_t xk = (2*k2PI)/((1.+alfa)*(TMath::Power(k2PI,1.5)));
   
   for(Int_t i=1; i<=200; i++){
      Double_t p = i*0.005;
      fPp[i] = p;
//      printf(" fPp[%d] = %f\n",i,fPp[i]);
      Double_t e1 = (p*p)/(2.*sig1*sig1);
      Double_t e2 = (p*p)/(2.*sig2*sig2);
      Double_t f1 = TMath::Exp(-(e1));
      Double_t f2 = TMath::Exp(-(e2));
      Double_t probp = xk*p*p*(f1/(TMath::Power(sig1,3.))+
                      alfa*f2/(TMath::Power(sig2,3.)))*0.005;
//      printf(" 	probp = %f\n",probp);
      fProbintp[i] = fProbintp[i-1] + probp;
      fProbintn[i] = fProbintp[i];
//      printf(" fProbintp[%d] = %f, fProbintp[%d] = %f\n",i,fProbintp[i],i,fProbintn[i]);
   }
} 
//_____________________________________________________________________________
void AliGenZDC::ExtractFermi(Int_t id, Double_t* fPp, Double_t* fProbintp,
                Double_t* fProbintn, Double_t* ddp)
{
//
// Compute Fermi momentum for spectator nucleons
//
  Int_t i;
  Float_t xx = gRandom->Rndm();
  if(id==kProton){
    for(i=0; i<=200; i++){
       if((xx>=fProbintp[i-1]) && (xx<fProbintp[i])) break;
       }
  }
  else if(id==kNeutron){
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
}

//_____________________________________________________________________________
void AliGenZDC::BeamDivCross(Int_t icross, Float_t fBeamDiv, Float_t fBeamCrossAngle, 
                Int_t fBeamCrossPlane, Double_t* pLab)
{
  Double_t tetpart, fipart, tetdiv, fidiv, angleSum[2], tetsum, fisum, dplab[3];
  Double_t rvec;

  Int_t i;
  
  Double_t pmq = 0.;
  for(i=0; i<=2; i++){
     dplab[i] = pLab[i];
     pmq = pmq+pLab[i]*pLab[i];
  }
  Double_t pmod = TMath::Sqrt(pmq);
//  printf("	pmod = %f\n",pmod);

//  printf("	icross = %d, fBeamDiv = %f\n",icross,fBeamDiv);
  if(icross==0){
    rvec = gRandom->Gaus(0.0,1.0);
    tetdiv = fBeamDiv * TMath::Abs(rvec);
    fidiv = (gRandom->Rndm())*k2PI;
  }
  else if(icross==1){
    if(fBeamCrossPlane==0.){
      tetdiv = 0.;
      fidiv = 0.;
    }
    else if(fBeamCrossPlane==1.){
      tetdiv = fBeamCrossAngle;
      fidiv = 0.;
    }
    else if(fBeamCrossPlane==2.){
      tetdiv = fBeamCrossAngle;
      fidiv = k2PI/4.;
    }
  }
//  printf("	tetdiv = %f, fidiv = %f\n",tetdiv,fidiv);
  tetpart = TMath::ATan(TMath::Sqrt(dplab[0]*dplab[0]+dplab[1]*dplab[1])/dplab[2]);
  if(dplab[1]!=0. || dplab[0]!=0.){
    fipart = TMath::ATan2(dplab[1],dplab[0]);
  }
  else{
    fipart = 0.;
  }
  if(fipart<0.) {fipart = fipart+k2PI;}
//  printf("	tetpart = %f, fipart = %f\n",tetpart,fipart);
  tetdiv = tetdiv*kRaddeg;
  fidiv = fidiv*kRaddeg;
  tetpart = tetpart*kRaddeg;
  fipart = fipart*kRaddeg;
  AddAngle(tetpart,fipart,tetdiv,fidiv,angleSum);
  tetsum = angleSum[0];
  fisum  = angleSum[1];
//  printf("	tetsum = %f, fisum = %f\n",tetsum,fisum);
  tetsum = tetsum*kDegrad;
  fisum = fisum*kDegrad;
  pLab[0] = pmod*TMath::Sin(tetsum)*TMath::Cos(fisum);
  pLab[1] = pmod*TMath::Sin(tetsum)*TMath::Sin(fisum);
  pLab[2] = pmod*TMath::Cos(tetsum);
//  printf("	pLab[0] = %f pLab[1] = %f pLab[2] = %f \n\n",
//         pLab[0],pLab[1],pLab[2]);
  for(i=0; i<=2; i++){
     fDivP[i] = pLab[i];
  }
}
  
//_____________________________________________________________________________
void  AliGenZDC::AddAngle(Double_t theta1, Double_t phi1, Double_t theta2,
  	       Double_t phi2, Double_t* angleSum)
{
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
//  printf("	AddAngle -> tetsum = %f, fisum = %f\n",tetsum, fisum); 
  angleSum[0] = tetsum;
  angleSum[1] = fisum;
}  

