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

// --- Standard libraries
#include <stdlib.h>

// --- ROOT system
#include <TRandom.h>
#include <TF1.h>

// --- AliRoot classes
#include "AliZDCFragment.h"
 
ClassImp(AliZDCFragment)
   
int comp(const void *i,const void *j) {return *(int *)i - *(int *)j;}


//_____________________________________________________________________________
AliZDCFragment::AliZDCFragment()
{
  //
  // Default constructor
  //
  fB = 0;
}

//_____________________________________________________________________________
AliZDCFragment::AliZDCFragment(Float_t b)
     : TNamed(" "," ")
{
  //
  // Standard constructor
  //
  fB = b;
  fZbAverage = 0;
  fNimf = 0;
  fZmax = 0;
  fTau = 0;
  for(Int_t i=0; i<=99; i++){
     fZZ[i] = 0;
     fNN[i] = 0;
  }
  fNalpha = 0;
  fZtot = 0;
  fNtot = 0;
  
}

//_____________________________________________________________________________
void AliZDCFragment::GenerateIMF(Int_t* fZZ, Int_t &fNalpha)
{
   // Coefficients of polynomial for average number of IMF
   const Float_t  ParamNimf[5]={0.011236,1.8364,56.572,-116.24,58.289}; 
   // Coefficients of polynomial for fluctuations on average number of IMF
   const Float_t  ParamFluctNimf[4]={-0.13176,2.9392,-5.2147,2.3092}; 
   // Coefficients of polynomial for average maximum Z of fragments
   const Float_t  ParamZmax[4]={0.16899,14.203,-2.8284,65.036}; 
   // Coefficients of polynomial for fluctuations on maximum Z of fragments
   const Float_t  ParamFluctZmax[5]={0.013782,-0.17282,1.5065,1.0654,-2.4317}; 
   // Coefficients of polynomial for exponent tau of fragments Z distribution
   const Float_t  ParamTau[3]={6.7233,-15.85,13.047};  
   //Coefficients of polynomial for average number of alphas
   const Float_t  ParamNalpha[4]={-0.68554,39.605,-68.311,30.165}; 
   // Coefficients of polynomial for fluctuations on average number of alphas
   const Float_t  ParamFluctNalpha[5]={0.283,6.2141,-17.113,17.394,-6.6084}; 
   // Coefficients of function for Pb nucleus skin
   const Float_t  ParamSkinPb[2]={0.93,11.05};
   
   // Thickness of nuclear surface
   const Float_t  NuclearThick = 0.52;
   // Maximum impact parameter for U [r0*A**(1/3)]
   const Float_t  bMaxU = 14.87;
   // Maximum impact parameter for Pb [r0*A**(1/3)]
   const Float_t  bMaxPb = 14.22;
   // Z of the projectile
   const Float_t  ZProj = 82.;
   
   // From b(Pb) to b(U)
   Float_t  bU = fB*bMaxU/bMaxPb;
    
   // From b(U) to Zbound(U) 
   // --- A.Schuttauf et al, Nuc.Phys. A607 (1996) 457 ---------------
   // From geometrical consideration and from dsigma/dZbound for U+U,
   // which is approx. constant, the constant value is found  
   // integrating the nucleus cross surface from 0 to bmax=R1+R2 where 
   // R = 1.2*A**(1/3). This value has been measured in Aladin (U+U).
   Float_t  ZbU = bU*bU*TMath::Pi()/7.48;
   
   //  Rescale Zbound for Pb
   fZbAverage = ZProj/92.*ZbU;
   
   // Zbound is proportional to b**2 up to b < bMaxPb-2*NuclearThick
   // and then it is an increasing exponential, imposing that at 
   // b=bMaxPb-2NuclearThick the two functions have the same derivative
   Float_t bCore = bMaxPb-2*NuclearThick;
   if(fB>bCore){
     fZbAverage=ZProj*(1.-TMath::Exp(-ParamSkinPb[0]*(fB-ParamSkinPb[1])));
   }
   if(fZbAverage>ZProj) fZbAverage = ZProj;
   Float_t ZbNorm = fZbAverage/ZProj;
   Float_t bNorm = fB/bMaxPb;
   
   // From Zbound to <Nimf>,<Zmax>,tau
   // Polinomial fits to Aladin distribution
   // --- A.Schuttauf et al, Nuc.Phys. A607 (1996) 457.
   Float_t AverageNimf = ParamNimf[0]+ParamNimf[1]*ZbNorm+ParamNimf[2]*
           TMath::Power(ZbNorm,2)+ParamNimf[3]*TMath::Power(ZbNorm,3)+
	   ParamNimf[4]*TMath::Power(ZbNorm,4);
   
   // Add fluctuation: from Singh et al. 
   Float_t FluctNimf = ParamFluctNimf[0]+ParamFluctNimf[1]*ZbNorm+
           ParamFluctNimf[2]*TMath::Power(ZbNorm,2)+ParamFluctNimf[3]
	   *TMath::Power(ZbNorm,3);
   Float_t xx = gRandom->Gaus(0.0,1.0);
   FluctNimf = FluctNimf*xx;
   fNimf = Int_t(AverageNimf+FluctNimf);
   Float_t y = gRandom->Rndm();
   if(y < ((AverageNimf+FluctNimf)-fNimf)) fNimf += 1;
   if(fNimf ==0 && ZbNorm>0.75) fNimf = 1;
   
   Float_t AverageZmax = ParamZmax[0]+ParamZmax[1]*ZbNorm+ParamZmax[2]*
           TMath::Power(ZbNorm,2)+ParamZmax[3]*TMath::Power(ZbNorm,3);
   fTau = ParamTau[0]+ParamTau[1]*ZbNorm+ParamTau[2]*TMath::Power(ZbNorm,2);
   
   // Add fluctuation to mean value of Zmax (see Hubele)
   Float_t FluctZmax = ParamFluctZmax[0]+ParamFluctZmax[1]*ZbNorm+
           ParamFluctZmax[2]*TMath::Power(ZbNorm,2)+ParamFluctZmax[3]*
	   TMath::Power(ZbNorm,3)+ParamFluctZmax[4]*TMath::Power(ZbNorm,4);
   FluctZmax = FluctZmax*ZProj/6.;
   Float_t xg = gRandom->Gaus(0.0,1.0);
   FluctZmax = FluctZmax*xg;
   fZmax = AverageZmax+FluctZmax;
   if(fZmax>ZProj) fZmax = ZProj;
   
//   printf("\n\n ------------------------------------------------------------");   
//   printf("\n Generation of nuclear fragments\n");   
//   printf("\n fNimf = %d\n", fNimf);   
//   printf("\n fZmax = %f\n", fZmax); 

   // Find the number of alpha particles 
   // from Singh et al. : Pb+emulsion
   Float_t AverageAlpha = ParamNalpha[0]+ParamNalpha[1]*ZbNorm+
           ParamNalpha[2]*TMath::Power(ZbNorm,2)+ParamNalpha[3]*
	   TMath::Power(ZbNorm,3);
   Float_t FluctAlpha = ParamFluctNalpha[0]+ParamFluctNalpha[1]*
           ZbNorm+ParamFluctNalpha[2]*TMath::Power(ZbNorm,2)+
	   ParamFluctNalpha[3]*TMath::Power(ZbNorm,3)+
	   ParamFluctNalpha[4]*TMath::Power(ZbNorm,4);
   Float_t xxx = gRandom->Gaus(0.0,1.0);
   FluctAlpha = FluctAlpha*xxx;
   fNalpha = Int_t(AverageAlpha+FluctAlpha);
   Float_t yy = gRandom->Rndm();
   if(yy < ((AverageAlpha+FluctAlpha)-fNalpha)) fNalpha += 1;

   // 2 possibilities:
   // 1) for bNorm < 0.9 ==> first remove alphas, then fragments
   // 2) for bNorm > 0.9 ==> first remove fragments, then alphas

   Int_t Choice = 0;
   Float_t ZbFrag = 0, SumZ = 0.;

   if(bNorm<=0.9) {
   // remove alpha from zbound to find zbound for fragments  (Z>=3)
     ZbFrag = fZbAverage-fNalpha*2;
     Choice = 1;
   }
   else {
     ZbFrag = fZbAverage;
     Choice = 0;
   }
//   printf("\n Choice = %d, fZbAverage = %f, ZbFrag = %f \n", Choice, fZbAverage, ZbFrag);
   
   
   // Check if ZbFrag < fZmax
   if(ZbFrag<=fZmax) {
     if(fNimf>0 && ZbFrag>=2){
       fNimf = 1;
       fZZ[0] = Int_t(ZbFrag);
       SumZ = ZbFrag;
     }
     else {
       fNimf = 0;
     }
     return;
   }
   
   // Prepare the exponential charge distribution dN/dZ
   if(fZmax <= 0.01) {
     fNimf = 0;
     return;
   }
   if(fNimf == 0) {
     fNimf = 0;
     return;
   }
   
   TF1 *funTau = new TF1("funTau","1./(x**[0])",0.01,fZmax);
   funTau->SetParameter(0,fTau);

   // Extract randomly the charge of the fragments from the distribution
 
   Float_t zz[fNimf];
   for(Int_t j=0; j<fNimf; j++){
      zz[j] =0;
   }
   for(Int_t i=0; i<fNimf; i++){
      zz[i] = Float_t(funTau->GetRandom());
//      printf("\n	zz[%d] = %f \n",i,zz[i]);
   }
   delete funTau;
   
   // Sorting vector in ascending order with C function QSORT 
   qsort((void*)zz,fNimf,sizeof(float),comp);

   
//   for(Int_t i=0; i<fNimf; i++){
//      printf("\n After sorting -> zz[%d] = %f \n",i,zz[i]);
//   }
   
   // Rescale the maximum charge to fZmax
   for(Int_t j=0; j<fNimf; j++){
     fZZ[j] = Int_t (zz[j]*fZmax/zz[fNimf-1]);
     if(fZZ[j]<3) fZZ[j] = 3;
//     printf("\n 	fZZ[%d] = %d \n",j,fZZ[j]);
   }
   
   // Check that the sum of the bound charges is not > than Zbound-Zalfa
   
   for(Int_t ii=0; ii<fNimf; ii++){
     SumZ += fZZ[ii];
   }
   
   Int_t k = 0;
   if(SumZ>ZbFrag){
     for(Int_t i=0; i< fNimf; i++){
       k += 1;
       SumZ -= fZZ[i];
       if(SumZ<=ZbFrag){
         fNimf -= (i+1);
         break;
       }
     }
   }
   else {
     if(Choice == 1) return;
     Int_t iDiff = Int_t((ZbFrag-SumZ)/2);
     if(iDiff<fNalpha){
       fNalpha=iDiff;
       return;
     }
     else{
       return;
     }
   }

   fNimf += k;
   for(Int_t i=0; i<fNimf; i++){
     fZZ[i] = fZZ[i+k];
   }
   fNimf -= k;
   
   SumZ=0;
   for(Int_t i=0; i<fNimf; i++){
     SumZ += fZZ[i];
   }
   
}

//_____________________________________________________________________________
void AliZDCFragment::AttachNeutrons(Int_t *fZZ, Int_t *fNN, Int_t &fZtot,Int_t &fNtot)
{
   const Float_t AIon[68]={1.87612,2.80943,3.7284,5.60305,6.53536,
     		     6.53622,8.39479,9.32699,10.2551,11.17793,
     		     13.04378,14.89917,17.6969,18.62284,21.41483,
     		     22.34193,25.13314,26.06034,28.85188,29.7818,
     		     32.57328,33.50356,36.29447,37.22492,41.87617,
     		     44.66324,47.45401,48.38228,51.17447,52.10307,
     		     54.89593,53.96644,58.61856,59.54963,68.85715,
     		     74.44178,78.16309,81.88358,83.74571,91.19832,
     		     98.64997,106.10997,111.68821,122.86796,
     		     128.45793,
     		     130.32111,141.51236,
     		     141.55,146.477,148.033,152.699,153.631,
     		     155.802,157.357,162.022,162.984,166.2624,
     		     168.554,171.349,173.4536,177.198,179.0518,
     		     180.675,183.473,188.1345,190.77,193.729,
     		     221.74295};
   const Int_t ZIon[68]={1,1,2,3,3,
     		     4,4,5,5,6,
     		     7,8,9,10,11,
     		     12,13,14,15,16,
     		     17,18,19,20,21,
     		     22,23,24,25,26,
     		     27,28,29,30,32,
     		     34,36,38,40,42,
     		     46,48,50,54,56,
     		     58,62,
     		     63,64,65,66,67,
     		     68,69,70,71,72,
     		     73,74,75,76,77,
     		     78,79,80,81,82,
      		     92};
    
   Int_t iZ, iA;  
//   printf("\n fNimf=%d\n",fNimf);  

   for(Int_t i=0; i<fNimf; i++) {
      for(Int_t j=0; j<68; j++) {
        iZ = ZIon[j];
	if((fZZ[i]-iZ) == 0){
	  iA = Int_t(AIon[j]/0.93149432+0.5);
	  fNN[i] = iA - iZ;
          break;
	}
	else if((fZZ[i]-iZ) < 0){
	  fZZ[i] = ZIon[j-1];
	  iA = Int_t (AIon[j-1]/0.93149432+0.5);
	  fNN[i] = iA - ZIon[j-1];
          break;
	}
      }
      fZtot += fZZ[i];
      fNtot += fNN[i];
   }		     
   

}
