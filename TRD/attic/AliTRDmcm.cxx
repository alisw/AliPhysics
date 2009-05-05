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

/* $Id: AliTRDmcm.cxx 29514 2008-10-26 10:24:38Z hristov $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//                                                                        //
//  Multi Chip Module exponential filter and tracklet finder              //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TMath.h>

#include "AliLog.h"

#include "AliTRDmcm.h"
#include "AliTRDtrigParam.h"

ClassImp(AliTRDmcm)

//_____________________________________________________________________________
AliTRDmcm::AliTRDmcm() 
  :TObject()
  ,fNtrk(0)
  ,fRobId(0)
  ,fChaId(0)
  ,fRow(0)
  ,fColFirst(0)
  ,fColLast(0)
  ,fTime1(0)
  ,fTime2(0)
  ,fClusThr(0)
  ,fPadThr(0)
  ,fNtrkSeeds(0)
  ,fR1(0)
  ,fR2(0)
  ,fC1(0)
  ,fC2(0)
  ,fPedestal(0)
  ,fId(0)
{
  //
  // AliTRDmcm default constructor
  //

  Int_t i = 0;
  Int_t j = 0;

  for (i = 0; i < kMaxTrackletsPerMCM; i++) {
    fTrkIndex[i] = 0;
    fSeedCol[i]  = -1;
  }
  for (i = 0; i < kMcmCol; i++) {
    fPadHits[i]  = 0;
    for (j = 0; j < kMcmTBmax; j++) {
      fADC[i][j]    = 0.0;
      fIsClus[i][j] = kFALSE;
    }
  }

}

//_____________________________________________________________________________
AliTRDmcm::AliTRDmcm(Int_t id) 
  :TObject()
  ,fNtrk(0)
  ,fRobId(0)
  ,fChaId(0)
  ,fRow(0)
  ,fColFirst(0)
  ,fColLast(0)
  ,fTime1(0)
  ,fTime2(0)
  ,fClusThr(0)
  ,fPadThr(0)
  ,fNtrkSeeds(0)
  ,fR1(0)
  ,fR2(0)
  ,fC1(0)
  ,fC2(0)
  ,fPedestal(0)
  ,fId(id)
{
  //
  // AliTRDmcm constructor
  //

  Int_t i = 0;
  Int_t j = 0;

  for (i = 0; i < kMaxTrackletsPerMCM; i++) {
    fTrkIndex[i] = 0;
    fSeedCol[i]  = -1;
  }
  for (i = 0; i < kMcmCol; i++) {
    fPadHits[i]  = 0;
    for (j = 0; j < kMcmTBmax; j++) {
      fADC[i][j]    = 0.0;
      fIsClus[i][j] = kFALSE;
    }
  }

  fTime1   = AliTRDtrigParam::Instance()->GetTime1();
  fTime2   = AliTRDtrigParam::Instance()->GetTime2();
  fClusThr = AliTRDtrigParam::Instance()->GetClusThr();
  fPadThr  = AliTRDtrigParam::Instance()->GetPadThr();
  
  AliTRDtrigParam::Instance()->GetFilterParam(fR1,fR2,fC1,fC2,fPedestal);

}

//_____________________________________________________________________________
AliTRDmcm::AliTRDmcm(const AliTRDmcm &m) 
  :TObject(m)
  ,fNtrk(m.fNtrk)
  ,fRobId(m.fRobId)
  ,fChaId(m.fChaId)
  ,fRow(m.fRow)
  ,fColFirst(m.fColFirst)
  ,fColLast(m.fColLast)
  ,fTime1(m.fTime1)
  ,fTime2(m.fTime2)
  ,fClusThr(m.fClusThr)
  ,fPadThr(m.fPadThr)
  ,fNtrkSeeds(m.fNtrkSeeds)
  ,fR1(m.fR1)
  ,fR2(m.fR2)
  ,fC1(m.fC1)
  ,fC2(m.fC2)
  ,fPedestal(m.fPedestal)
  ,fId(m.fId)
{
  //
  // AliTRDmcm copy constructor
  //

  Int_t i = 0;
  Int_t j = 0;

  for (i = 0; i < kMaxTrackletsPerMCM; i++) {
    ((AliTRDmcm &) m).fTrkIndex[i] = 0;
    ((AliTRDmcm &) m).fSeedCol[i]  = -1;
  }
  for (i = 0; i < kMcmCol; i++) {
    ((AliTRDmcm &) m).fPadHits[i]  = 0;
    for (j = 0; j < kMcmTBmax; j++) {
      ((AliTRDmcm &) m).fADC[i][j]    = 0.0;
      ((AliTRDmcm &) m).fIsClus[i][j] = kFALSE;
    }
  }

}

//_____________________________________________________________________________
AliTRDmcm::~AliTRDmcm() 
{
  //
  // AliTRDmcm destructor
  //

}

//_____________________________________________________________________________
AliTRDmcm &AliTRDmcm::operator=(const AliTRDmcm &m)
{
  //
  // Assignment operator
  //

  if (this != &m) ((AliTRDmcm &) m).Copy(*this); 

  return *this;

}

//_____________________________________________________________________________
void AliTRDmcm::Copy(TObject &m) const
{
  //
  // Copy function
  //

  Int_t i = 0;
  Int_t j = 0;

  ((AliTRDmcm &) m).fNtrk       = fNtrk;
  ((AliTRDmcm &) m).fRobId      = fRobId;
  ((AliTRDmcm &) m).fChaId      = fChaId;
  ((AliTRDmcm &) m).fRow        = fRow;
  ((AliTRDmcm &) m).fColFirst   = fColFirst;
  ((AliTRDmcm &) m).fColLast    = fColLast;
  ((AliTRDmcm &) m).fTime1      = fTime1;
  ((AliTRDmcm &) m).fTime2      = fTime2;
  ((AliTRDmcm &) m).fClusThr    = fClusThr;
  ((AliTRDmcm &) m).fPadThr     = fPadThr;
  ((AliTRDmcm &) m).fNtrkSeeds  = fNtrkSeeds;
  ((AliTRDmcm &) m).fR1         = fR1;
  ((AliTRDmcm &) m).fR2         = fR2;
  ((AliTRDmcm &) m).fC1         = fC1;
  ((AliTRDmcm &) m).fC2         = fC2;
  ((AliTRDmcm &) m).fPedestal   = fPedestal;
  ((AliTRDmcm &) m).fId         = fId;

  for (i = 0; i < kMaxTrackletsPerMCM; i++) {
    ((AliTRDmcm &) m).fTrkIndex[i] = 0;
    ((AliTRDmcm &) m).fSeedCol[i]  = -1;
  }
  for (i = 0; i < kMcmCol; i++) {
    ((AliTRDmcm &) m).fPadHits[i]  = 0;
    for (j = 0; j < kMcmTBmax; j++) {
      ((AliTRDmcm &) m).fADC[i][j]    = 0.0;
      ((AliTRDmcm &) m).fIsClus[i][j] = kFALSE;
    }
  }

}

//_____________________________________________________________________________
void AliTRDmcm::AddTrk(Int_t id) 
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

  Int_t i = 0;
  Int_t j = 0;

  for (i = 0; i < kMcmCol; i++) {
    fPadHits[i] = 0;
    for (j = 0; j < kMcmTBmax; j++) {
      fADC[i][j]    = 0.0;
      fIsClus[i][j] = kFALSE;
    }
  }
  for (i = 0; i < kMaxTrackletsPerMCM; i++) {
    fSeedCol[i] = -1;
  }
  
}

//_____________________________________________________________________________
Bool_t AliTRDmcm::Run()
{
  //
  // Run MCM
  //

  AliDebug(2,(Form("Run MCM %d\n",Id())));

  Int_t   iTime = 0;
  Int_t   iCol  = 0;
  Int_t   iPlus = 0;
  Int_t   i     = 0;
  Int_t   j     = 0;

  Float_t amp[3] = { 0.0, 0.0, 0.0 };
  Int_t   nClus;
  Int_t   clusCol[kMcmCol/2];
  Float_t clusAmp[kMcmCol/2];
  Float_t veryLarge;
  Int_t   clusMin = -1;
  
  // Main TB loop
  for (iTime = fTime1; iTime <= fTime2; iTime++) {  

    // Find clusters...
    nClus = 0;
    for (iCol = 1; iCol < (kMcmCol-1); iCol++) {
      amp[0] = fADC[iCol-1][iTime];
      amp[1] = fADC[iCol  ][iTime];
      amp[2] = fADC[iCol+1][iTime];
      if (IsCluster(amp)) {
	fIsClus[iCol][iTime] = kTRUE;
	clusCol[nClus]       = iCol;
	clusAmp[nClus]       = amp[0]+amp[1]+amp[2];
	nClus++;
	if (nClus == kMcmCol/2) {
	  AliWarning(Form("Too many clusters in time bin %2d MCM %d...\n",iTime,Id()));
	  break;
	}
      }
    }

    // ...but no more than six...
    if (nClus > (Int_t) kSelClus) {
      for (j = kSelClus/2; j < nClus-kSelClus/2; j++) {
	fIsClus[clusCol[j]][iTime] = kFALSE;
      } 
    }

    // ...and take the largest four.
    Int_t nClusPlus = nClus - kMaxClus;
    for (iPlus = 0; iPlus < nClusPlus; iPlus++ ) {
      veryLarge = 1.E+10;
      for (i = 0; i < nClus; i++) {
	if (fIsClus[clusCol[i]][iTime]) {
	  if (clusAmp[i] <= veryLarge) {
	    veryLarge = clusAmp[i];
	    clusMin = i;
	  }
	}
      }
      fIsClus[clusCol[clusMin]][iTime] = kFALSE;
    }

    AddTimeBin(iTime);

  }  // end main TB loop
    
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

  Int_t i      = 0;
  Int_t j      = 0;
  Int_t nSeeds = 0;

  AliDebug(2,Form("AliTRDmcm::CreateSeeds MCM %d \n",Id()));

  // Working array for hit sums
  Int_t fHit2padSum[2][kMcmCol];            

  // Initialize the array
  for (i = 0; i < 2; i++) {
    for (j = 0; j < kMcmCol; j++ ) {
      if (i == 0) {
	fHit2padSum[i][j] = j; 
      } else {
	fHit2padSum[i][j] = -1; 
      }
    }
  }

  Int_t sum10 = AliTRDtrigParam::Instance()->GetSum10();
  Int_t sum12 = AliTRDtrigParam::Instance()->GetSum12();

  // Build the 2padSum
  Int_t nsum2seed = 0;
  for (i = 0; i < kMcmCol; i++) {
    if (i < (kMcmCol-1)) {
      if ((fPadHits[i] >= sum10) && ((fPadHits[i] + fPadHits[i+1]) >= sum12)) {
	fHit2padSum[1][i] = fPadHits[i] + fPadHits[i+1]; 
      } 
      else {
	fHit2padSum[1][i] = -1;
      }
    } 
    else {
      if (fPadHits[i] >= sum12) {
	fHit2padSum[1][i] = fPadHits[i]; 
      } 
      else {
	fHit2padSum[1][i] = -1;
      }
    }
    if (fHit2padSum[1][i] > 0) {
      nsum2seed++;
    }
  }

  // sort the sums in decreasing order of the amplitude	
  Sort(kMcmCol,&fHit2padSum[0][0],&fHit2padSum[1][0],1);

  // arrange (maximum number of) candidates in increasing order of the column number
  nSeeds = TMath::Min(nsum2seed,kMaxTrackletsPerMCM);
  Sort(nSeeds,&fHit2padSum[1][0],&fHit2padSum[0][0],0);

  for (i = 0; i < nSeeds; i++) {
    fSeedCol[i] = fHit2padSum[0][i];
  }

  // reject multiple found tracklets
  Int_t imax = nSeeds - 1;
  for (i = 0; i < imax; i++) {

    if ((fHit2padSum[0][i]+1) == fHit2padSum[0][i+1]) {
      nSeeds--;
      if (fHit2padSum[1][i] >= fHit2padSum[1][i+1]) {
	AliDebug(2,Form("Reject seed %1d in col %02d. \n",i,fHit2padSum[0][i+1]));
	fSeedCol[i+1] = -1;
      } 
      else {
	AliDebug(2,Form("Reject seed %1d in col %02d. \n",i,fHit2padSum[0][i]));
	fSeedCol[i]   = -1;
      }
    }

  }

  return nSeeds;

}

//_____________________________________________________________________________
void AliTRDmcm::Sort(Int_t nel, Int_t *x1, Int_t *x2, Int_t dir) const
{
  //
  // Sort two parallel vectors (x1[nel], x2[nel]) after the second one (x2)
  // in the direction: dir = 0 ascending order
  //                   dir = 1 descending order
  //

  Int_t  i = 0;
  Bool_t sort;
  Int_t  tmp1;
  Int_t  tmp2;

  if (dir == 0) {

    do { 
      sort = kTRUE;
      for (i = 0; i < (nel-1); i++) {
	if (x2[i+1] < x2[i]) {
	  tmp2    = x2[i]; 
	  x2[i]   = x2[i+1]; 
	  x2[i+1] = tmp2;
	  tmp1    = x1[i]; 
	  x1[i]   = x1[i+1]; 
	  x1[i+1] = tmp1;
	  sort    = kFALSE;
	}
      }
    } while (sort == kFALSE);

  }

  if (dir == 1) {

    do { 
      sort = kTRUE;
      for (i = 0; i < (nel-1); i++) {
	if (x2[i+1] > x2[i]) {
	  tmp2    = x2[i]; 
	  x2[i]   = x2[i+1]; 
	  x2[i+1] = tmp2;
	  tmp1    = x1[i]; 
	  x1[i]   = x1[i+1]; 
	  x1[i+1] = tmp1;
	  sort    = kFALSE;
	}
      }
    } while (sort == kFALSE);

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
Bool_t AliTRDmcm::IsCluster(Float_t amp[3]) const
{
  //
  // Find if the amplitudes amp[0], amp[1], amp[2] are a cluster
  //

  // -> shape
  if (amp[0] > amp[1] || amp[2] > amp[1]) {
    return kFALSE;
  }

  // -> cluster amplitude
  if ((amp[0]+amp[1]+amp[2]) < fClusThr) {
    return kFALSE;
  }

  // -> pad amplitude
  if (amp[0] < fPadThr && amp[2] < fPadThr) {
    return kFALSE;
  }

  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDmcm::Filter(Int_t nexp, Int_t ftype)
{
  //
  // Exponential filter
  //

  Int_t iCol  = 0;
  Int_t iTime = 0;

  Double_t sour[kMcmTBmax];
  Double_t dtarg[kMcmTBmax];
  Int_t    itarg[kMcmTBmax];

  switch(ftype) {

   case 0:

    for (iCol = 0; iCol < kMcmCol; iCol++) {
      for (iTime = 0; iTime < kMcmTBmax; iTime++) {
	sour[iTime] = fADC[iCol][iTime];
      }
      DeConvExpA(sour,dtarg,kMcmTBmax,nexp);
      for (iTime = 0; iTime < kMcmTBmax; iTime++) {
	fADC[iCol][iTime] = TMath::Max(0.0,dtarg[iTime]);
      }
    }
    break;

  case 1:

    for (iCol = 0; iCol < kMcmCol; iCol++) {
      for (iTime = 0; iTime < kMcmTBmax; iTime++) {
	sour[iTime] = fADC[iCol][iTime];
      }
      DeConvExpD(sour,itarg,kMcmTBmax,nexp);
      for (iTime = 0; iTime < kMcmTBmax; iTime++) {
	fADC[iCol][iTime] = itarg[iTime];
      }
    }
    break;

  case 2:

    for (iCol = 0; iCol < kMcmCol; iCol++) {
      for (iTime = 0; iTime < kMcmTBmax; iTime++) {
	sour[iTime] = fADC[iCol][iTime];
      }
      DeConvExpMI(sour,dtarg,kMcmTBmax);
      for (iTime = 0; iTime < kMcmTBmax; iTime++) {
	fADC[iCol][iTime] = TMath::Max(0.0,dtarg[iTime]);
      }
    }
    break;

  default:

    AliError(Form("Invalid filter type %d ! \n",ftype));
    return;

  }

}

//_____________________________________________________________________________
void AliTRDmcm::DeConvExpA(Double_t *source, Double_t *target, Int_t n, Int_t nexp) 
{
  //
  // Exponential filter "analog"
  //

  Int_t    i = 0;
  Int_t    k = 0;
  Double_t reminder[2];
  Double_t correction;
  Double_t result;
  Double_t rates[2];
  Double_t coefficients[2];

  // Initialize (coefficient = alpha, rates = lambda)
  
  // FilterOpt.C (aliroot@pel:/homel/aliroot/root/work/beamt/CERN02)
  Double_t r1, r2, c1, c2;
  r1 = (Double_t)fR1;
  r2 = (Double_t)fR2;
  c1 = (Double_t)fC1;
  c2 = (Double_t)fC2;
  
  coefficients[0] = c1;
  coefficients[1] = c2;

  Double_t dt = 0.1;
  rates[0] = TMath::Exp(-dt/(r1));
  rates[1] = TMath::Exp(-dt/(r2));

  // Attention: computation order is important
  correction = 0.0;
  for (k = 0; k < nexp; k++) {
    reminder[k] = 0.0;
  }
    
  for (i = 0; i < n; i++) {

    result    = (source[i] - correction);    // no rescaling
    target[i] = result;
    
    for (k = 0; k < nexp; k++) {
      reminder[k] = rates[k] * (reminder[k] + coefficients[k] * result);
    }
      
    correction = 0.0;
    for (k = 0; k < nexp; k++) {
      correction += reminder[k];
    }

  }
  
}

//_____________________________________________________________________________
void AliTRDmcm::DeConvExpD(Double_t *source, Int_t *target, Int_t n, Int_t nexp) 
{
  //
  // Exponential filter "digital"
  //

  Int_t i = 0;

  Int_t fAlphaL;
  Int_t fAlphaS;
  Int_t fLambdaL;
  Int_t fLambdaS;
  Int_t fTailPed;

  Int_t iAlphaL;
  Int_t iAlphaS;
  Int_t iLambdaL;
  Int_t iLambdaS;

  // FilterOpt.C (aliroot@pel:/homel/aliroot/root/work/beamt/CERN02)
  // initialize (coefficient = alpha, rates = lambda)

  fLambdaL = 0; 
  fAlphaL  = 0; 
  fLambdaS = 0; 
  fAlphaS  = 0;
  iLambdaL = 0; 
  iAlphaL  = 0; 
  iLambdaS = 0; 
  iAlphaS  = 0;

  Double_t dt = 0.1;

  Double_t r1;
  Double_t r2;
  Double_t c1;
  Double_t c2;
  r1 = (Double_t) fR1;
  r2 = (Double_t) fR2;
  c1 = (Double_t) fC1;
  c2 = (Double_t) fC2;

  fLambdaL = (Int_t)((TMath::Exp(-dt/r1) - 0.75) * 2048.0);
  fLambdaS = (Int_t)((TMath::Exp(-dt/r2) - 0.25) * 2048.0);
  iLambdaL = fLambdaL & 0x01FF; iLambdaL |= 0x0600; 	//  9 bit paramter + fixed bits
  iLambdaS = fLambdaS & 0x01FF; iLambdaS |= 0x0200; 	//  9 bit paramter + fixed bits

  if (nexp == 1) {
    fAlphaL = (Int_t) (c1 * 2048.0);
    iAlphaL = fAlphaL & 0x03FF;				// 10 bit paramter
  }
  if (nexp == 2) {
    fAlphaL = (Int_t) (c1 * 2048.0);
    fAlphaS = (Int_t) ((c2 - 0.5) * 2048.0);
    iAlphaL = fAlphaL & 0x03FF;				// 10 bit paramter
    iAlphaS = fAlphaS & 0x03FF; iAlphaS |= 0x0400;	        // 10 bit paramter + fixed bits
  }
  
  Double_t iAl = iAlphaL  / 2048.0;	       // alpha L: correspondence to floating point numbers
  Double_t iAs = iAlphaS  / 2048.0;	       // alpha S: correspondence to floating point numbers
  Double_t iLl = iLambdaL / 2048.0;	       // lambda L: correspondence to floating point numbers
  Double_t iLs = iLambdaS / 2048.0;	       // lambda S: correspondence to floating point numbers

  Int_t h1;
  Int_t h2;
  Int_t rem1;
  Int_t rem2;
  Int_t correction;
  Int_t result;
  Int_t iFactor = ((Int_t) fPedestal) << 2;

  Double_t xi = 1 - (iLl*iAs + iLs*iAl);	     // Calculation of equilibrium values of the
  rem1 = (Int_t) ((iFactor/xi) * ((1-iLs)*iLl*iAl)); // Internal registers to prevent switch on effects.
  rem2 = (Int_t) ((iFactor/xi) * ((1-iLl)*iLs*iAs));
  
  // further initialization
  if ((rem1 + rem2) > 0x0FFF) {
    correction = 0x0FFF;
  } 
  else {
    correction = (rem1 + rem2) & 0x0FFF;
  }

  fTailPed = iFactor - correction;

  for (i = 0; i < n; i++) {

    result = ((Int_t)source[i] - correction);
    if (result < 0) {			        
      result = 0;
    }

    target[i] = result;
                                                        
    h1 = (rem1 + ((iAlphaL * result) >> 11));
    if (h1 > 0x0FFF) {
      h1 = 0x0FFF;
    } 
    else {
      h1 &= 0x0FFF;
    }

    h2 = (rem2 + ((iAlphaS * result) >> 11));
    if (h2 > 0x0FFF) {
      h2 = 0x0FFF;
    } 
    else {
      h2 &= 0x0FFF;
    }
  
    rem1 = (iLambdaL * h1 ) >> 11;
    rem2 = (iLambdaS * h2 ) >> 11;
    
    if ((rem1 + rem2) > 0x0FFF) {
      correction = 0x0FFF;
    } 
    else {
      correction = (rem1 + rem2) & 0x0FFF;
    }

  }

}

//_____________________________________________________________________________
void AliTRDmcm::DeConvExpMI(Double_t *source, Double_t *target, Int_t n) 
{
  //
  // Exponential filter (M. Ivanov)
  //

  Int_t i = 0;

  Double_t sig1[100];
  Double_t sig2[100];
  Double_t sig3[100];

  for (i = 0; i < n; i++) {
    sig1[i] = source[i];
  }

  Float_t dt = 0.1;

  Float_t lambda0 = (1.0 / fR2) * dt;
  Float_t lambda1 = (1.0 / fR1) * dt;

  TailMakerSpline(sig1,sig2,lambda0,n);
  TailCancelationMI(sig2,sig3,0.7,lambda1,n);

  for (i = 0; i < n; i++) {
    target[i] = sig3[i];
  }

}

//______________________________________________________________________________
void AliTRDmcm::TailMakerSpline(Double_t *ampin, Double_t *ampout, Double_t lambda, Int_t n) 
{
  //
  // Special filter (M. Ivanov)
  //

  Int_t    i = 0;

  Double_t l = TMath::Exp(-lambda*0.5);
  Double_t in[1000];
  Double_t out[1000];

  // Initialize in[] and out[] goes 0 ... 2*n+19
  for (i = 0; i < n*2+20; i++) {
    in[i]  = 0;
    out[i] = 0;
  }

  // in[] goes 0, 1
  in[0] = ampin[0];
  in[1] = (ampin[0] + ampin[1]) * 0.5;
   
  // Add charge to the end
  for (i = 0; i < 22; i++) {
    // in[] goes 2*n-2, 2*n-1, ... , 2*n+19 
    in[2*(n-1)+i] = ampin[n-1];
  }

  // Use arithmetic mean
  for (i = 1; i < n-1; i++) {
    // in[] goes 2, 3, ... , 2*n-4, 2*n-3
    in[2*i]   = ampin[i];
    in[2*i+1] = ((ampin[i]+ampin[i+1]))/2.;
  }

  Double_t temp;
  out[2*n]    = in[2*n];
  temp        = 0;
  for (i = 2*n; i >= 0; i--) {
    out[i]    = in[i] + temp;
    temp      = l*(temp+in[i]);
  }

  for (i = 0; i < n; i++){
    //ampout[i] = out[2*i+1];  // org
    ampout[i] = out[2*i];
  }

}

//______________________________________________________________________________
void AliTRDmcm::TailCancelationMI(Double_t *ampin, Double_t *ampout
                                , Double_t norm, Double_t lambda, Int_t n) 
{
  //
  // Special filter (M. Ivanov)
  //

  Int_t    i = 0;

  Double_t l = TMath::Exp(-lambda*0.5);
  Double_t k = l*(1.0 - norm*lambda*0.5);
  Double_t in[1000];
  Double_t out[1000];

  // Initialize in[] and out[] goes 0 ... 2*n+19
  for (i = 0; i < n*2+20; i++) {
    in[i]  = 0;
    out[i] = 0;
  }

  // in[] goes 0, 1
  in[0] = ampin[0];
  in[1] = (ampin[0]+ampin[1])*0.5;

  // Add charge to the end
  for (i =-2; i < 22; i++) {
    // in[] goes 2*n-4, 2*n-3, ... , 2*n+19 
    in[2*(n-1)+i] = ampin[n-1];
  }

  for (i = 1; i < n-2; i++) {
    // in[] goes 2, 3, ... , 2*n-6, 2*n-5
    in[2*i]    = ampin[i];
    in[2*i+1]  = (9.0 * (ampin[i]+ampin[i+1]) - (ampin[i-1]+ampin[i+2])) / 16.0;
    //in[2*i+1]  = ((ampin[i]+ampin[i+1]))/2.0;
  }

  Double_t temp;
  out[0] = in[0];
  temp   = in[0];
  for (i = 1; i <= 2*n; i++) {
    out[i] = in[i] + (k-l)*temp;
    temp   = in[i] +  k   *temp;
  }

  for (i = 0; i < n; i++) {
    //ampout[i] = out[2*i+1];  // org
    //ampout[i] = TMath::Max(out[2*i+1],0.0);  // org
    ampout[i] = TMath::Max(out[2*i],0.0);
  }

}

