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


// ******************************************************************
//
//  	Class for nuclear fragments formation
//
// ******************************************************************

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
AliZDCFragment::AliZDCFragment():
  fB(0),
  fZbAverage(0),
  fNimf(0),
  fZmax(0),
  fTau(0),
  fNalpha(0),
  fZtot(0),
  fNtot(0)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliZDCFragment::AliZDCFragment(Float_t b): 
  TNamed(" "," "),
  fB(b),
  fZbAverage(0),
  fNimf(0),
  fZmax(0),
  fTau(0),
  fNalpha(0),
  fZtot(0),
  fNtot(0)
{
  //
  // Standard constructor
  //
  for(Int_t i=0; i<=99; i++){
     fZZ[i] = 0;
     fNN[i] = 0;
  }
  
}

//_____________________________________________________________________________
void AliZDCFragment::GenerateIMF()
{

   // Loop variables
  Int_t i,j;

   // Coefficients of polynomial for average number of IMF
   const Float_t  kParamNimf[5]={0.011236,1.8364,56.572,-116.24,58.289}; 
   // Coefficients of polynomial for fluctuations on average number of IMF
   const Float_t  kParamFluctNimf[4]={-0.13176,2.9392,-5.2147,2.3092}; 
   // Coefficients of polynomial for average maximum Z of fragments
   //const Float_t  kParamZmax[4]={0.16899,14.203,-2.8284,65.036}; 
   const Float_t  kParamZmax[4]={0.16899,14.203,-2.8284,70.5}; 
   // Coefficients of polynomial for fluctuations on maximum Z of fragments
   const Float_t  kParamFluctZmax[5]={0.013782,-0.17282,1.5065,1.0654,-2.4317}; 
   // Coefficients of polynomial for exponent tau of fragments Z distribution
   const Float_t  kParamTau[3]={6.7233,-15.85,13.047};  
   //Coefficients of polynomial for average number of alphas
   const Float_t  kParamNalpha[4]={-0.68554,39.605,-68.311,30.165}; 
   // Coefficients of polynomial for fluctuations on average number of alphas
   const Float_t  kParamFluctNalpha[5]={0.283,6.2141,-17.113,17.394,-6.6084}; 
   // Coefficients of function for Pb nucleus skin
   const Float_t  kParamSkinPb[2]={0.762408, 20.};
   
   // Thickness of nuclear surface
   const Float_t  kNuclearThick = 0.52;
   // Maximum impact parameter for U [r0*A**(1/3)]
   const Float_t  kbMaxU = 14.87;
   // Maximum impact parameter for Pb [r0*A**(1/3)]
   //const Float_t  kbMaxPb = 14.22+4*kNuclearThick;
   const Float_t  kbMaxPb = 14.22;
   // Z of the projectile
   const Float_t  kZProj = 82.;
   
   // From b(Pb) to b(U)
   if(fB>kbMaxPb) fB = 2*kbMaxPb-fB;
   
   Float_t  bU = fB*kbMaxU/kbMaxPb;
    
   // From b(U) to Zbound(U) 
   // --- A.Schuttauf et al, Nuc.Phys. A607 (1996) 457 ---------------
   // From geometrical consideration and from dsigma/dZbound for U+U,
   // which is approx. constant, the constant value is found  
   // integrating the nucleus cross surface from 0 to bmax=R1+R2 where 
   // R = 1.2*A**(1/3). This value has been measured in Aladin (U+U).
   Float_t  zbU = bU*bU*TMath::Pi()/7.48;
   
   //  Rescale Zbound for Pb
   fZbAverage = kZProj/92.*zbU;
   
   // Zbound is proportional to b**2 up to b < kbMaxPb-2*kNuclearThick
   // and then it is an increasing exponential, imposing that at 
   // b=kbMaxPb-2kNuclearThick the two functions have the same derivative
   //Float_t bCore = kbMaxPb-2*kNuclearThick;
   if(fB>kbMaxPb){
     fZbAverage = TMath::Exp(-kParamSkinPb[0]*(fB-kParamSkinPb[1]));
     printf(" b = %1.2f fm   Z_bound %1.2f\n", fB, fZbAverage);
   }
   if(fZbAverage>kZProj) fZbAverage = kZProj;
   Float_t zbNorm = fZbAverage/kZProj;
   Float_t bNorm = fB/kbMaxPb;
   
   // From Zbound to <Nimf>,<Zmax>,tau
   // Polinomial fits to Aladin distribution
   // --- A.Schuttauf et al, Nuc.Phys. A607 (1996) 457.
   Float_t averageNimf = kParamNimf[0]+kParamNimf[1]*zbNorm+kParamNimf[2]*
           TMath::Power(zbNorm,2)+kParamNimf[3]*TMath::Power(zbNorm,3)+
	   kParamNimf[4]*TMath::Power(zbNorm,4);
   
   // Add fluctuation: from Singh et al. 
   Float_t fluctNimf = kParamFluctNimf[0]+kParamFluctNimf[1]*zbNorm+
           kParamFluctNimf[2]*TMath::Power(zbNorm,2)+kParamFluctNimf[3]
	   *TMath::Power(zbNorm,3);
   Float_t xx = gRandom->Gaus(0.0,1.0);
   fluctNimf = fluctNimf*xx;
   fNimf = Int_t(averageNimf+fluctNimf);
   Float_t y = gRandom->Rndm();
   if(y < ((averageNimf+fluctNimf)-fNimf)) fNimf += 1;
   if(fNimf ==0 && zbNorm>0.75) fNimf = 1;
   
   Float_t averageZmax = kParamZmax[0]+kParamZmax[1]*zbNorm+kParamZmax[2]*
           TMath::Power(zbNorm,2)+kParamZmax[3]*TMath::Power(zbNorm,3);
   fTau = kParamTau[0]+kParamTau[1]*zbNorm+kParamTau[2]*TMath::Power(zbNorm,2);
   
   // Add fluctuation to mean value of Zmax (see Hubele)
   Float_t fluctZmax = kParamFluctZmax[0]+kParamFluctZmax[1]*zbNorm+
           kParamFluctZmax[2]*TMath::Power(zbNorm,2)+kParamFluctZmax[3]*
	   TMath::Power(zbNorm,3)+kParamFluctZmax[4]*TMath::Power(zbNorm,4);
   fluctZmax = fluctZmax*kZProj/6.;
   Float_t xg = gRandom->Gaus(0.0,1.0);
   fluctZmax = fluctZmax*xg;
   fZmax = (averageZmax+fluctZmax);
   if(fZmax>kZProj) fZmax = kZProj;
   
//   printf("\n\n ------------------------------------------------------------");   
//   printf("\n Generation of nuclear fragments\n");   
//   printf("\n fNimf = %d\n", fNimf);   
//   printf("\n fZmax = %f\n", fZmax); 

   // Find the number of alpha particles 
   // from Singh et al. : Pb+emulsion
   Float_t averageAlpha = kParamNalpha[0]+kParamNalpha[1]*zbNorm+
           kParamNalpha[2]*TMath::Power(zbNorm,2)+kParamNalpha[3]*
	   TMath::Power(zbNorm,3);
   Float_t fluctAlpha = kParamFluctNalpha[0]+kParamFluctNalpha[1]*
           zbNorm+kParamFluctNalpha[2]*TMath::Power(zbNorm,2)+
	   kParamFluctNalpha[3]*TMath::Power(zbNorm,3)+
	   kParamFluctNalpha[4]*TMath::Power(zbNorm,4);
   Float_t xxx = gRandom->Gaus(0.0,1.0);
   fluctAlpha = fluctAlpha*xxx;
   fNalpha = Int_t(averageAlpha+fluctAlpha);
   Float_t yy = gRandom->Rndm();
   if(yy < ((averageAlpha+fluctAlpha)-fNalpha)) fNalpha += 1;

   // 2 possibilities:
   // 1) for bNorm < 0.9 ==> first remove alphas, then fragments
   // 2) for bNorm > 0.9 ==> first remove fragments, then alphas

   Int_t choice = 0;
   Float_t zbFrag = 0, sumZ = 0.;

   if(bNorm<=0.9) {
   // remove alpha from zbound to find zbound for fragments  (Z>=3)
     zbFrag = fZbAverage-fNalpha*2;
     choice = 1;
   }
   else {
     zbFrag = fZbAverage;
     choice = 0;
   }
//   printf("\n choice = %d, fZbAverage = %f, zbFrag = %f \n", choice, fZbAverage, zbFrag);
   
   
   // Check if zbFrag < fZmax
   if(zbFrag<=fZmax) {
     if(fNimf>0 && zbFrag>=2){
       fNimf = 1;
       fZZ[0] = Int_t(zbFrag);
       sumZ = zbFrag;
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
 
   Float_t * zz = new Float_t[fNimf];
   for(j=0; j<fNimf; j++){
      zz[j] =0;
   }
   for(i=0; i<fNimf; i++){
      zz[i] = Float_t(funTau->GetRandom());
//      printf("\n	zz[%d] = %f \n",i,zz[i]);
   }
   delete funTau;
   
   // Sorting vector in ascending order with C function QSORT 
   qsort((void*)zz,fNimf,sizeof(Float_t),comp);

   
//   for(Int_t i=0; i<fNimf; i++){
//      printf("\n After sorting -> zz[%d] = %f \n",i,zz[i]);
//   }
   
   // Rescale the maximum charge to fZmax
   for(j=0; j<fNimf; j++){
     fZZ[j] = Int_t (zz[j]*fZmax/zz[fNimf-1]);
     if(fZZ[j]<3) fZZ[j] = 3;
//     printf("\n 	fZZ[%d] = %d \n",j,fZZ[j]);
   }

   delete[] zz;
   
   // Check that the sum of the bound charges is not > than Zbound-Zalfa
   
   for(Int_t ii=0; ii<fNimf; ii++){
     sumZ += fZZ[ii];
   }
   
   Int_t k = 0;
   if(sumZ>zbFrag){
     for(i=0; i< fNimf; i++){
       k += 1;
       sumZ -= fZZ[i];
       if(sumZ<=zbFrag){
         fNimf -= (i+1);
         break;
       }
     }
   }
   else {
     if(choice == 1) return;
     Int_t iDiff = Int_t((zbFrag-sumZ)/2);
     if(iDiff<fNalpha){
       fNalpha=iDiff;
       return;
     }
     else{
       return;
     }
   }

   fNimf += k;
   for(i=0; i<fNimf; i++){
     fZZ[i] = fZZ[i+k];
   }
   fNimf -= k;
   
   sumZ=0;
   for(i=0; i<fNimf; i++){
     sumZ += fZZ[i];
   }
   
}

//_____________________________________________________________________________
void AliZDCFragment::AttachNeutrons()
{
//
// Prepare nuclear fragment by attaching a suitable number of neutrons
//
   const Float_t kAIon[68]={1.87612,2.80943,3.7284,5.60305,6.53536,
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
   const Int_t kZIon[68]={1,1,2,3,3,
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
        iZ = kZIon[j];
	if((fZZ[i]-iZ) == 0){
	  iA = Int_t(kAIon[j]/0.93149432+0.5);
	  fNN[i] = iA - iZ;
          break;
	}
	else if((fZZ[i]-iZ) < 0){
	  fZZ[i] = kZIon[j-1];
	  iA = Int_t (kAIon[j-1]/0.93149432+0.5);
	  fNN[i] = iA - kZIon[j-1];
          break;
	}
      }
      fZtot += fZZ[i];
      fNtot += fNN[i];
   }		     
   

}

//_____________________________________________________________________________
Float_t AliZDCFragment::DeuteronNumber()
{
    // Calculates the fraction of deuterum nucleus produced
    //
    Float_t deuteronProdPar[2] = {-0.068,0.0385};
    Float_t deutNum = deuteronProdPar[0] + deuteronProdPar[1]*fB;
    if(deutNum<0.) deutNum = 0.;
    return deutNum;
}
