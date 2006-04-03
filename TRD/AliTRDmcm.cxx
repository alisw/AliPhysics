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
Revision 1.1.1.1  2004/08/19 14:58:11  vulpescu
CVS head

Revision 1.1.1.1  2004/08/18 07:47:17  vulpescu
test

*/

#include <TMath.h>
#include <TH2F.h>

#include "AliTRDmcm.h"
#include "AliTRDtrigParam.h"

ClassImp(AliTRDmcm)

//_____________________________________________________________________________
AliTRDmcm::AliTRDmcm() 
{

  //
  // AliTRDmcm default constructor
  //

  fTrigParam = 0;

  fNtrk = 0;
  for (Int_t i = 0; i < kMaxTrackletsPerMCM; i++) {
    fTrkIndex[i] = 0;
  }
  fRobId = 0;
  fChaId = 0;
  fRow = 0;
  fColFirst = 0;
  fColLast = 0;
  for (Int_t i = 0; i < kMcmCol; i++) {
    fPadHits[i] = 0;
    for (Int_t j = 0; j < kMcmTBmax; j++) {
      fADC[i][j] = 0.0;
      fIsClus[i][j] = kFALSE;
    }
  }
  fPadThr = 0;
  fClusThr = 0;
  fTime1 = 0;
  fTime2 = 0;
  fNtrkSeeds = 0;
  for (Int_t i = 0; i < kMaxTrackletsPerMCM; i++) {
    fSeedCol[i] = -1;
  }

  fR1 = 0.0;
  fR2 = 0.0;
  fC1 = 0.0;
  fC2 = 0.0;
  fPedestal = 0.0;

  fId = 0;

}

//_____________________________________________________________________________
AliTRDmcm::AliTRDmcm(AliTRDtrigParam *trigp, const Int_t id) 
{
  //
  // AliTRDmcm constructor
  //

  fTrigParam = trigp;

  fNtrk = 0;
  for (Int_t i = 0; i < kMaxTrackletsPerMCM; i++) {
    fTrkIndex[i] = 0;
  }
  fRobId = 0;
  fChaId = 0;
  fRow = 0;
  fColFirst = 0;
  fColLast = 0;
  for (Int_t i = 0; i < kMcmCol; i++) {
    fPadHits[i] = 0;
    for (Int_t j = 0; j < kMcmTBmax; j++) {
      fADC[i][j] = 0.0;
      fIsClus[i][j] = kFALSE;
    }
  }
  fPadThr  = fTrigParam->GetPadThr();
  fClusThr = fTrigParam->GetClusThr();
  fTime1 = fTrigParam->GetTime1();
  fTime2 = fTrigParam->GetTime2();
  fNtrkSeeds = 0;
  for (Int_t i = 0; i < kMaxTrackletsPerMCM; i++) {
    fSeedCol[i] = -1;
  }
  
  fR1 = 0.0;
  fR2 = 0.0;
  fC1 = 0.0;
  fC2 = 0.0;
  fPedestal = 0.0;

  fTrigParam->GetFilterParam(fR1,fR2,fC1,fC2,fPedestal);

  fId = id;

}

//_____________________________________________________________________________
AliTRDmcm::~AliTRDmcm() 
{

  //
  // AliTRDmcm destructor
  //

}

//_____________________________________________________________________________
void AliTRDmcm::AddTrk(const Int_t id) 
{
  //
  // Add a tracklet index
  //

  fTrkIndex[fNtrk] = id;
  fNtrk++;

  return;

}

//_____________________________________________________________________________
void AliTRDmcm::Reset()
{
  //
  // Reset MCM data
  //

  for (Int_t i = 0; i < kMcmCol; i++) {
    fPadHits[i] = 0;
    for (Int_t j = 0; j < kMcmTBmax; j++) {
      fADC[i][j] = 0.0;
      fIsClus[i][j] = kFALSE;
    }
  }
  for (Int_t i = 0; i < kMaxTrackletsPerMCM; i++) {
    fSeedCol[i] = -1;
  }
  
}

//_____________________________________________________________________________
Bool_t AliTRDmcm::Run()
{
  //
  // Run MCM
  //

  if ( fTrigParam->GetDebugLevel() > 1 ) printf("AliTRDmcm::Run MCM %d\n",Id());

  Float_t Amp[3] = {0.0, 0.0, 0.0};
  Int_t   nClus;
  Int_t   ClusCol[kMcmCol/2];
  Float_t ClusAmp[kMcmCol/2];
  Float_t VeryLarge;
  Int_t   ClusMin = -1;
  
  for (Int_t iTime = fTime1; iTime <= fTime2; iTime++) {  // main TB loop

    // find clusters...
    nClus = 0;
    for (Int_t iCol = 1; iCol < (kMcmCol-1); iCol++) {
      Amp[0] = fADC[iCol-1][iTime];
      Amp[1] = fADC[iCol  ][iTime];
      Amp[2] = fADC[iCol+1][iTime];
      if (IsCluster(Amp)) {
	fIsClus[iCol][iTime] = kTRUE;
	ClusCol[nClus] = iCol;
	ClusAmp[nClus] = Amp[0]+Amp[1]+Amp[2];
	nClus++;
	if (nClus == kMcmCol/2) {
	  printf("Too many clusters in time bin %2d MCM %d...\n",iTime,Id());
	  //return kFALSE;
	  break;
	}
      }
    }

    // ...but no more than six...
    if (nClus > (Int_t)kSelClus) {
      for (Int_t j = kSelClus/2; j < nClus-kSelClus/2; j++) {
	fIsClus[ClusCol[j]][iTime] = kFALSE;
      } 
    }

    // ...and take the largest four.

    Int_t nClusPlus = nClus - kMaxClus;
    for (Int_t iPlus = 0; iPlus < nClusPlus; iPlus++ ) {
      VeryLarge = 1.E+10;
      for (Int_t i = 0; i < nClus; i++) {
	if (fIsClus[ClusCol[i]][iTime]) {
	  if (ClusAmp[i] <= VeryLarge) {
	    VeryLarge = ClusAmp[i];
	    ClusMin = i;
	  }
	}
      }
      fIsClus[ClusCol[ClusMin]][iTime] = kFALSE;
    }

    AddTimeBin(iTime);

  }  // end main TB loop
  
  if (fTrigParam->GetDebugLevel() > 1) {
    for (Int_t i = fTime1; i <= fTime2; i++) {
      printf("%2d: ",i);
      for (Int_t j = 0; j < kMcmCol; j++) {
	printf("%1d ",fIsClus[j][i]);
      }
      printf("\n");
    }
    printf("PadHits:     ");
    for (Int_t iPad = 0; iPad < kMcmCol; iPad++) {
      printf("%2d ",fPadHits[iPad]);
    }
    printf("\n");
  }
  
  if ((fNtrkSeeds = CreateSeeds())) {

    return kTRUE;

  }

  return kFALSE;

}
//_____________________________________________________________________________
Int_t AliTRDmcm::CreateSeeds()
{
  //
  // Make column seeds (from Falk Lesser, ex KIP)
  //

  if ( fTrigParam->GetDebugLevel() > 1 ) printf("AliTRDmcm::CreateSeeds MCM %d \n",Id());

  Int_t nSeeds = 0;

  // working array for hit sums
  Int_t fHit2padSum[2][kMcmCol];            

  // initialize the array
  for( Int_t i = 0; i < 2; i++ ) {
    for( Int_t j = 0; j < kMcmCol; j++ ) {
      if( i == 0 ) {
	fHit2padSum[i][j] = j; 
      } else {
	fHit2padSum[i][j] = -1; 
      }
    }
  }

  Int_t Sum10 = fTrigParam->GetSum10();
  Int_t Sum12 = fTrigParam->GetSum12();

  // build the 2padSum
  Int_t Nsum2seed = 0;
  for( Int_t i = 0; i < kMcmCol; i++ ) {
    if( i < (kMcmCol-1) ) {
      if( (fPadHits[i] >= Sum10) && ((fPadHits[i] + fPadHits[i+1]) >= Sum12) ) {
	fHit2padSum[1][i] = fPadHits[i] + fPadHits[i+1]; 
      } else {
	fHit2padSum[1][i] = -1;
      }
    } else {
      if ( fPadHits[i] >= Sum12 ) {
	fHit2padSum[1][i] = fPadHits[i]; 
      } else {
	fHit2padSum[1][i] = -1;
      }
    }
    if (fHit2padSum[1][i] > 0) Nsum2seed++;
  }

  if (fTrigParam->GetDebugLevel() > 1) {
    printf("fHit2padSum: ");
    for( Int_t i = 0; i < kMcmCol; i++ ) {
      printf("%2d ",fHit2padSum[0][i]);
    }
    printf("\n");
    printf("             ");
    for( Int_t i = 0; i < kMcmCol; i++ ) {
      printf("%2d ",fHit2padSum[1][i]);
    }
    printf("\n");
  }

  // sort the sums in decreasing order of the amplitude	
  Sort(kMcmCol,&fHit2padSum[0][0],&fHit2padSum[1][0],1);

  if (fTrigParam->GetDebugLevel() > 1) {
    printf("fHit2padSum: ");
    for( Int_t i = 0; i < kMcmCol; i++ ) {
      printf("%2d ",fHit2padSum[0][i]);
    }
    printf("\n");
    printf("             ");
    for( Int_t i = 0; i < kMcmCol; i++ ) {
      printf("%2d ",fHit2padSum[1][i]);
    }
    printf("\n");
  }

  // arrange (maximum number of) candidates in increasing order of the column number
  nSeeds = TMath::Min(Nsum2seed,kMaxTrackletsPerMCM);
  Sort(nSeeds,&fHit2padSum[1][0],&fHit2padSum[0][0],0);

  for (Int_t i = 0; i < nSeeds; i++) {
    fSeedCol[i] = fHit2padSum[0][i];
  }

  if (fTrigParam->GetDebugLevel() > 1) {
    printf("Found %d seeds before multiple rejection. \n",nSeeds);
    printf("fHit2padSum: ");
    for( Int_t i = 0; i < kMcmCol; i++ ) {
      printf("%2d ",fHit2padSum[0][i]);
    }
    printf("\n");
    printf("             ");
    for( Int_t i = 0; i < kMcmCol; i++ ) {
      printf("%2d ",fHit2padSum[1][i]);
    }
    printf("\n");
  }

  // reject multiple found tracklets
  Int_t imax = nSeeds-1;
  for (Int_t i = 0; i < imax; i++) {

    if ((fHit2padSum[0][i]+1) == fHit2padSum[0][i+1]) {
      nSeeds--;
      if (fHit2padSum[1][i] >= fHit2padSum[1][i+1]) {
	if (fTrigParam->GetDebugLevel() > 1) 
	  printf("Reject seed %1d in col %02d. \n",i,fHit2padSum[0][i+1]);
	fSeedCol[i+1] = -1;
      } else {
	if (fTrigParam->GetDebugLevel() > 1) 
	  printf("Reject seed %1d in col %02d. \n",i,fHit2padSum[0][i]);
	fSeedCol[i] = -1;
      }
    }

  }

  if ( fTrigParam->GetDebugLevel() > 1 ) {
    printf("Found %d seeds in MCM %d ",nSeeds,Id());
    for (Int_t i = 0; i < (imax+1); i++) {
      if (fSeedCol[i] >= 0) printf(", %02d ",fSeedCol[i]);
    }
    printf("\n");
  }

  return nSeeds;

}

//_____________________________________________________________________________
void AliTRDmcm::Sort(Int_t nel, Int_t *x1, Int_t *x2, Int_t dir)
{

  // Sort two parallel vectors (x1[nel], x2[nel]) after the second one (x2)
  // in the direction: dir = 0 ascending order
  //                   dir = 1 descending order

  Bool_t sort;
  Int_t tmp1, tmp2;

  if ( dir == 0 ) {

    do { 
      sort = kTRUE;
      for ( Int_t i = 0; i < (nel-1); i++ )
	if ( x2[i+1] < x2[i] ) {
	  tmp2 = x2[i]; 
	  x2[i] = x2[i+1]; 
	  x2[i+1] = tmp2;
	  tmp1 = x1[i]; 
	  x1[i] = x1[i+1]; 
	  x1[i+1] = tmp1;
	  sort = kFALSE;
	}
    } while ( sort == kFALSE );

  }

  if ( dir == 1 ) {

    do { 
      sort = kTRUE;
      for ( Int_t i = 0; i < (nel-1); i++ )
	if ( x2[i+1] > x2[i] ) {
	  tmp2 = x2[i]; 
	  x2[i] = x2[i+1]; 
	  x2[i+1] = tmp2;
	  tmp1 = x1[i]; 
	  x1[i] = x1[i+1]; 
	  x1[i+1] = tmp1;
	  sort = kFALSE;
	}
    } while ( sort == kFALSE );

  }

}

//_____________________________________________________________________________
void AliTRDmcm::AddTimeBin(const Int_t iTime)
{
  //
  // Build column seeds
  //

  for (Int_t iPad = 1; iPad < (kMcmCol-1); iPad++) {
    if (fIsClus[iPad][iTime]) {
      fPadHits[iPad]++;
    }
  }

}

//_____________________________________________________________________________
Bool_t AliTRDmcm::IsCluster(Float_t amp[3])
{
  //
  // Find if the amplitudes amp[0], amp[1], amp[2] are a cluster
  //

  // -> shape
  if (amp[0] > amp[1] || amp[2] > amp[1]) return kFALSE;

  // -> cluster amplitude
  if ((amp[0]+amp[1]+amp[2]) < fClusThr) return kFALSE;

  // -> pad amplitude
  if (amp[0] < fPadThr && amp[2] < fPadThr) return kFALSE;

  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDmcm::Filter(Int_t nexp, Int_t ftype)
{
  //
  // exponential filter
  //

  Double_t sour[kMcmTBmax];
  Double_t Dtarg[kMcmTBmax];
  Int_t    Itarg[kMcmTBmax];

  switch(ftype) {

  case 0:

    for (Int_t iCol = 0; iCol < kMcmCol; iCol++) {
      for (Int_t iTime = 0; iTime < kMcmTBmax; iTime++) {
	sour[iTime] = fADC[iCol][iTime];
      }
      DeConvExpA(sour,Dtarg,kMcmTBmax,nexp);
      for (Int_t iTime = 0; iTime < kMcmTBmax; iTime++) {
	//fADC[iCol][iTime] = (Int_t)TMath::Max(0.0,Dtarg[iTime]);
	fADC[iCol][iTime] = TMath::Max(0.0,Dtarg[iTime]);
      }
    }
    break;

  case 1:

    for (Int_t iCol = 0; iCol < kMcmCol; iCol++) {
      for (Int_t iTime = 0; iTime < kMcmTBmax; iTime++) {
	sour[iTime] = fADC[iCol][iTime];
      }
      DeConvExpD(sour,Itarg,kMcmTBmax,nexp);
      for (Int_t iTime = 0; iTime < kMcmTBmax; iTime++) {
	fADC[iCol][iTime] = Itarg[iTime];
      }
    }
    break;

  case 2:

    for (Int_t iCol = 0; iCol < kMcmCol; iCol++) {
      for (Int_t iTime = 0; iTime < kMcmTBmax; iTime++) {
	sour[iTime] = fADC[iCol][iTime];
      }
      DeConvExpMI(sour,Dtarg,kMcmTBmax);
      for (Int_t iTime = 0; iTime < kMcmTBmax; iTime++) {
	//fADC[iCol][iTime] = (Int_t)TMath::Max(0.0,Dtarg[iTime]);
	fADC[iCol][iTime] = TMath::Max(0.0,Dtarg[iTime]);
      }
    }
    break;

  default:

    printf("Invalid filter type %d ! \n",ftype);
    return;

  }

}

//_____________________________________________________________________________
void AliTRDmcm::DeConvExpA(Double_t *source, Double_t *target, Int_t n, Int_t nexp) 
{

  Double_t rates[2];
  Double_t coefficients[2];

  // initialize (coefficient = alpha, rates = lambda)
  
  // FilterOpt.C (aliroot@pel:/homel/aliroot/root/work/beamt/CERN02)
  Double_t R1, R2, C1, C2;
  R1 = (Double_t)fR1;
  R2 = (Double_t)fR2;
  C1 = (Double_t)fC1;
  C2 = (Double_t)fC2;
  
  coefficients[0] = C1;
  coefficients[1] = C2;

  Double_t Dt = 0.100;
  rates[0] = TMath::Exp(-Dt/(R1));
  rates[1] = TMath::Exp(-Dt/(R2));
  
  Int_t i, k;
  Double_t reminder[2];
  Double_t correction, result;

  /* attention: computation order is important */
  correction=0.0;
  for ( k=0; k<nexp; k++ ) reminder[k]=0.0;
    
  for ( i=0; i<n; i++ ) {
    result = ( source[i] - correction );    // no rescaling
    target[i] = result;
    
    for ( k=0; k<nexp; k++ ) reminder[k] = rates[k] * ( reminder[k] + coefficients[k] * result);
      
    correction=0.0;
    for ( k=0; k<nexp; k++ ) correction += reminder[k];
  }
  
}

//_____________________________________________________________________________
void AliTRDmcm::DeConvExpD(Double_t *source, Int_t *target, Int_t n, Int_t nexp) {

  Int_t fAlpha_l, fAlpha_s, fLambda_l, fLambda_s, fTailPed;
  Int_t iAlpha_l, iAlpha_s, iLambda_l, iLambda_s;

  // FilterOpt.C (aliroot@pel:/homel/aliroot/root/work/beamt/CERN02)
  // initialize (coefficient = alpha, rates = lambda)

  fLambda_l = 0; fAlpha_l = 0; fLambda_s = 0; fAlpha_s = 0;
  iLambda_l = 0; iAlpha_l = 0; iLambda_s = 0; iAlpha_s = 0;

  Double_t Dt = 0.100;

  Double_t R1, R2, C1, C2;
  R1 = (Double_t)fR1;
  R2 = (Double_t)fR2;
  C1 = (Double_t)fC1;
  C2 = (Double_t)fC2;

  fLambda_l = (Int_t)((TMath::Exp(-Dt/R1)-0.75)*2048.0);
  fLambda_s = (Int_t)((TMath::Exp(-Dt/R2)-0.25)*2048.0);
  iLambda_l = fLambda_l & 0x01FF; iLambda_l |= 0x0600; 	//  9 bit paramter + fixed bits
  iLambda_s = fLambda_s & 0x01FF; iLambda_s |= 0x0200; 	//  9 bit paramter + fixed bits

  if (nexp == 1) {
    fAlpha_l = (Int_t)(C1*2048.0);
    iAlpha_l = fAlpha_l & 0x03FF;				// 10 bit paramter
  }
  if (nexp == 2) {
    fAlpha_l = (Int_t)(C1*2048.0);
    fAlpha_s = (Int_t)((C2-0.5)*2048.0);
    iAlpha_l = fAlpha_l & 0x03FF;				// 10 bit paramter
    iAlpha_s = fAlpha_s & 0x03FF; iAlpha_s |= 0x0400;	        // 10 bit paramter + fixed bits
  }
  
  Double_t iAl = iAlpha_l / 2048.0;	       // alpha L: correspondence to floating point numbers
  Double_t iAs = iAlpha_s / 2048.0;	       // alpha S: correspondence to floating point numbers
  Double_t iLl = iLambda_l / 2048.0;	       // lambda L: correspondence to floating point numbers
  Double_t iLs = iLambda_s / 2048.0;	       // lambda S: correspondence to floating point numbers

  Int_t h1,h2;
  Int_t rem1, rem2;
  Int_t correction, result;
  Int_t iFactor = ((Int_t)fPedestal)<<2;

  Double_t xi = 1-(iLl*iAs + iLs*iAl);		    // calculation of equilibrium values of the
  rem1 = (Int_t)((iFactor/xi) * ((1-iLs)*iLl*iAl)); // internal registers to prevent switch on effects.
  rem2 = (Int_t)((iFactor/xi) * ((1-iLl)*iLs*iAs));
  
  // further initialization
  if ( (rem1 + rem2) > 0x0FFF ) {
    correction = 0x0FFF;
  } else {
    correction = (rem1 + rem2) & 0x0FFF;
  }

  fTailPed = iFactor - correction;

  for (Int_t i=0; i<n; i++) {

    result = ( (Int_t)source[i] - correction );
    if ( result<0 ) {			        
      result = 0;
    }

    target[i] = result;
                                                        
    h1 = ( rem1 + ((iAlpha_l * result) >> 11) );
    if ( h1 > 0x0FFF ) {
      h1 = 0x0FFF;
    } else {
      h1 &= 0x0FFF;
    }

    h2 = ( rem2 + ((iAlpha_s * result) >> 11));
    if ( h2 > 0x0FFF ) {
      h2 = 0x0FFF;
    } else {
      h2 &= 0x0FFF;
    }
  
    rem1 = (iLambda_l * h1 ) >> 11;
    rem2 = (iLambda_s * h2 ) >> 11;
    
    if ( (rem1 + rem2) > 0x0FFF ) {
      correction = 0x0FFF;
    } else {
      correction = (rem1 + rem2) & 0x0FFF;
    }
  }

}

//_____________________________________________________________________________
void AliTRDmcm::DeConvExpMI(Double_t *source, Double_t *target, Int_t n) {

   Double_t Sig1[100], Sig2[100], Sig3[100];//, Sig4[100];
   for (Int_t i = 0; i < n; i++) Sig1[i] = source[i];

   Float_t Dt = 0.100;

   //Float_t lambda0 = 9.8016*Dt;  // short
   //Float_t lambda1 = 1.0778*Dt;  // long

   Float_t lambda0 = (1.0/fR2)*Dt;
   Float_t lambda1 = (1.0/fR1)*Dt;

   TailMakerSpline(Sig1,Sig2,lambda0,n);
   TailCancelationMI(Sig2,Sig3,0.7,lambda1,n);

   for (Int_t i = 0; i < n; i++) target[i] = Sig3[i];

}

//______________________________________________________________________________
void AliTRDmcm::TailMakerSpline(Double_t *ampin, Double_t *ampout, Double_t lambda, Int_t n) {

   Double_t l = TMath::Exp(-lambda*0.5);
   //
   //
   Double_t in[1000];
   Double_t out[1000];
   // initialize in[] and out[] goes 0 ... 2*n+19
   for (Int_t i=0; i<n*2+20; i++){
     in[i] = 0;
     out[i]= 0;
   }

   // in[] goes 0, 1
   in[0] = ampin[0];
   in[1] = (ampin[0]+ampin[1])*0.5;
   
   // add charge to the end
   for (Int_t i=0; i<22; i++){
     // in[] goes 2*n-2, 2*n-1, ... , 2*n+19 
     in[2*(n-1)+i] = ampin[n-1];
   }

   // use arithm mean
   for (Int_t i=1; i<n-1; i++){
     // in[] goes 2, 3, ... , 2*n-4, 2*n-3
     in[2*i]    = ampin[i];
     in[2*i+1]  = ((ampin[i]+ampin[i+1]))/2.;
   }
   /*
   // add charge to the end
   for (Int_t i=-2; i<22; i++){
     // in[] goes 2*n-4, 2*n-3, ... , 2*n+19 
     in[2*(n-1)+i] = ampin[n-1];
   }

   // use spline mean
   for (Int_t i=1; i<n-2; i++){
     // in[] goes 2, 3, ... , 2*n-6, 2*n-5
     in[2*i]    = ampin[i];
     in[2*i+1]  = (9.*(ampin[i]+ampin[i+1])-(ampin[i-1]+ampin[i+2]))/16;
   }
   */
   //
   Double_t temp;
   out[2*n]      = in[2*n];
   temp          = 0;
   for (int i=2*n; i>=0; i--){
     out[i]      = in[i]   + temp;
     temp        = l*(temp+in[i]);
   }

   //
   for (int i=0;i<n;i++){
     //ampout[i]      = out[2*i+1];  // org
     ampout[i]      = out[2*i];
   }

}

//______________________________________________________________________________
void AliTRDmcm::TailCancelationMI(Double_t *ampin, Double_t *ampout, Double_t norm, Double_t lambda, Int_t n) {

   Double_t L = TMath::Exp(-lambda*0.5);
   Double_t K = L*(1.-norm*lambda*0.5);
   //
   //
   Double_t in[1000];
   Double_t out[1000];
   // initialize in[] and out[] goes 0 ... 2*n+19
   for (Int_t i=0; i<n*2+20; i++){
     in[i] = 0;
     out[i]= 0;
   }

   // in[] goes 0, 1
   in[0] = ampin[0];
   in[1] = (ampin[0]+ampin[1])*0.5;

   // add charge to the end
   for (Int_t i=-2; i<22; i++){
     // in[] goes 2*n-4, 2*n-3, ... , 2*n+19 
     in[2*(n-1)+i] = ampin[n-1];
   }

   //
   for (Int_t i=1; i<n-2; i++){
     // in[] goes 2, 3, ... , 2*n-6, 2*n-5
     in[2*i]    = ampin[i];
     in[2*i+1]  = (9.*(ampin[i]+ampin[i+1])-(ampin[i-1]+ampin[i+2]))/16.;
     //in[2*i+1]  = ((ampin[i]+ampin[i+1]))/2.;
   }

   //
   Double_t temp;
   out[0]     = in[0];
   temp       = in[0];
   for (int i=1; i<=2*n; i++){
     out[i]      = in[i]   + (K-L)*temp;
     temp        = in[i]   +  K*temp;
   }
   //
   //
   for (int i=0; i<n; i++){
     //ampout[i]      = out[2*i+1];  // org
     //ampout[i]      = TMath::Max(out[2*i+1], 0.0);  // org
     ampout[i]      = TMath::Max(out[2*i], 0.0);
   }

}

