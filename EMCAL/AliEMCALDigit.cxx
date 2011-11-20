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

//_________________________________________________________________________
//  EMCAL digit: 
//      A Digit is the sum of the energy lost in an EMCAL Tower
//      It also stores information on Primary, and enterring particle
//      tracknumbers Digits are created using AliEMCALSDigitizer, followed
//      by AliEMCALDigitizer 
//
//*-- Author: Sahal Yacoob (LBL)
// based on : AliPHOSDigit
//__________________________________________________________________________

// --- ROOT system ---
#include <Riostream.h>
#include <TMath.h>

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"

ClassImp(AliEMCALDigit)

//____________________________________________________________________________
AliEMCALDigit::AliEMCALDigit() :  
AliDigitNew(), 
  fAmpFloat(0.),
  fNSamples(0),
  fSamples(0x0),
  fNSamplesHG(0),
  fSamplesHG(0x0),
  fNprimary(0),
  fNMaxPrimary(5),
  fPrimary(0x0),
  fDEPrimary(0x0),
  fNiparent(0),
  fNMaxiparent(5), 
  fIparent(0x0),
  fDEParent(0x0),
  fMaxIter(0),
  fTime(0.), 
  fTimeR(0.),
  fChi2(0.),
  fNDF(0),
  fDigitType(kUnknown),
  fAmpCalib(-1)
{
  // default ctor 

  // Need to initialise for reading old files
  fPrimary   = new Int_t[fNMaxPrimary] ;
  fDEPrimary = new Float_t[fNMaxPrimary] ;
  fIparent   = new Int_t[fNMaxiparent] ; 
  fDEParent  = new Float_t[fNMaxiparent] ; 
  for ( Int_t i = 0; i < fNMaxPrimary ; i++) {
    fPrimary[i]    = -1 ;
    fDEPrimary[i]  = 0 ;
  } 

  for ( Int_t i = 0; i < fNMaxiparent ; i++) {
    fIparent[i]  = -1 ;
    fDEParent[i] = 0 ;
  }
}

//____________________________________________________________________________
AliEMCALDigit::AliEMCALDigit(Int_t primary, Int_t iparent, Int_t id, Float_t digEnergy, Float_t time, Int_t type, Int_t index, Float_t chi2, Int_t ndf, Float_t dE) 
  : AliDigitNew(),
    fAmpFloat(digEnergy),
    fNSamples(0),
    fSamples(0x0),
    fNSamplesHG(0),
    fSamplesHG(0x0),
    fNprimary(0),
    fNMaxPrimary(25),
    fPrimary(0x0),
    fDEPrimary(0x0),
    fNiparent(0),
    fNMaxiparent(150),
    fIparent(0x0),
    fDEParent(0x0),
    fMaxIter(5),
    fTime(time),
    fTimeR(time),
    fChi2(chi2),
    fNDF(ndf),
    fDigitType(type),
    fAmpCalib(-1)
{  
  // ctor with all data 

  // data memebrs of the base class (AliNewDigit)
  fAmp         = 0;
  fId          = id ;
  fIndexInList = index ; 

  // data member  
  fPrimary   = new Int_t[fNMaxPrimary] ;
  fDEPrimary = new Float_t[fNMaxPrimary] ;
  fIparent   = new Int_t[fNMaxiparent] ; 
  fDEParent  = new Float_t[fNMaxiparent] ; 
  if( primary != -1){
    fNprimary    = 1 ; 
    fPrimary[0]  = primary ;  
    fDEPrimary[0] = dE ;  
    fNiparent   = 1 ;
    fIparent[0] = iparent ;    
    fDEParent[0] = dE ;  
  }
  else{  //If the contribution of this primary smaller than fDigitThreshold (AliEMCALv1)
    fNprimary    = 0 ; 
    fPrimary[0]  = -1 ;
    fDEPrimary[0]  = 0 ;
    fNiparent   = 0 ;
    fIparent[0] = -1 ;  
    fDEParent[0]  = 0 ;  
  }
  Int_t i ;
  for ( i = 1; i < fNMaxPrimary ; i++) {
    fPrimary[i]  = -1 ; 
    fDEPrimary[i]  = 0 ;
  } 

  for ( i = 1; i< fNMaxiparent ; i++) {
    fIparent[i] = -1 ;  
    fDEParent[i] = 0 ;
  }  
}

//____________________________________________________________________________
AliEMCALDigit::AliEMCALDigit(const AliEMCALDigit & digit) 
  : AliDigitNew(digit),
    fAmpFloat(digit.fAmpFloat),	
    fNSamples(digit.fNSamples),
    fSamples(),
    fNSamplesHG(digit.fNSamplesHG),
    fSamplesHG(),
    fNprimary(digit.fNprimary),
    fNMaxPrimary(digit.fNMaxPrimary),
    fPrimary(),
    fDEPrimary(),
    fNiparent(digit.fNiparent),
    fNMaxiparent(digit.fNMaxiparent),
    fIparent(),
    fDEParent(),
    fMaxIter(digit.fMaxIter),
    fTime(digit.fTime),
    fTimeR(digit.fTimeR), 
    fChi2(digit.fChi2), 
    fNDF(digit.fNDF),
    fDigitType(digit.fDigitType),
    fAmpCalib(digit.fAmpCalib)
{
  // copy ctor
  // data memebrs of the base class (AliNewDigit)
  fAmp         = digit.fAmp ;
  fId          = digit.fId;
  fIndexInList = digit.fIndexInList ; 

  // data members  
  if (fSamples )   delete [] fSamples ;   fSamples   = NULL ;
  if (fSamplesHG ) delete [] fSamplesHG ; fSamplesHG = NULL ;
  if (fPrimary  )  delete [] fPrimary ;   fPrimary   = NULL ;
  if (fDEPrimary)  delete [] fDEPrimary ; fDEPrimary = NULL ;
  if (fIparent )   delete [] fIparent ;   fIparent   = NULL ;
  if (fDEParent)   delete [] fDEParent ;  fDEParent  = NULL ;

  if (fNSamples){
    fSamples     = new Int_t[fNSamples]; 
    for (Int_t i=0; i < digit.fNSamples; i++) fSamples[i] = digit.fSamples[i];
  }
  
  if (fNSamplesHG){
    fSamplesHG   = new Int_t[fNSamplesHG]; 
    for (Int_t i=0; i < digit.fNSamplesHG; i++) fSamplesHG[i] = digit.fSamplesHG[i];
  }
  
 
  if (fNMaxPrimary){
    fPrimary   = new Int_t  [fNMaxPrimary] ;  
    fDEPrimary = new Float_t[fNMaxPrimary] ;
    for ( Int_t i = 0; i < fNMaxPrimary ; i++) {
      fPrimary[i]   = digit.fPrimary[i] ;
      fDEPrimary[i] = digit.fDEPrimary[i] ;
    }
  }
  
  if (fNMaxiparent){
    fIparent  = new Int_t  [fNMaxiparent] ;
    fDEParent = new Float_t[fNMaxiparent] ;
    for (Int_t j = 0; j< fNMaxiparent ; j++) {
    fIparent[j]  = digit.fIparent[j] ;
    fDEParent[j] = digit.fDEParent[j] ;
    }
  }  
 
}

//____________________________________________________________________________
AliEMCALDigit::~AliEMCALDigit() 
{
  // Delete array of primaries if any

  if (fSamples )   delete [] fSamples ;   fSamples   = NULL ;
  if (fSamplesHG ) delete [] fSamplesHG ; fSamplesHG = NULL ;
  if (fPrimary  )  delete [] fPrimary ;   fPrimary   = NULL ;
  if (fDEPrimary)  delete [] fDEPrimary ; fDEPrimary = NULL ;
  if (fIparent )   delete [] fIparent ;   fIparent   = NULL ;
  if (fDEParent)   delete [] fDEParent ;  fDEParent  = NULL ;

}


//____________________________________________________________________________
void AliEMCALDigit::Clear(const Option_t*) 
{
  // Delete array of primaries if any
  if (fSamples )   delete [] fSamples ;   fSamples   = NULL ;
  if (fSamplesHG ) delete [] fSamplesHG ; fSamplesHG = NULL ;
  if (fPrimary  )  delete [] fPrimary ;   fPrimary   = NULL ;
  if (fDEPrimary)  delete [] fDEPrimary ; fDEPrimary = NULL ;
  if (fIparent )   delete [] fIparent ;   fIparent   = NULL ;
  if (fDEParent)   delete [] fDEParent ;  fDEParent  = NULL ;

}

//____________________________________________________________________________
Int_t AliEMCALDigit::Compare(const TObject * obj) const
{
  // Compares two digits with respect to its Id
  // to sort according increasing Id

  Int_t rv = 2 ;

  AliEMCALDigit * digit = (AliEMCALDigit *)obj ; 

  Int_t iddiff = fId - digit->GetId() ; 

  if ( iddiff > 0 ) 
    rv = 1 ;
  else if ( iddiff < 0 )
    rv = -1 ; 
  else
    rv = 0 ;
  
  return rv ; 

}

//____________________________________________________________________________
Float_t AliEMCALDigit::GetEta() const
{ 
  //return pseudorapidity for this digit
  // should be change in EMCALGeometry - 19-nov-04
  Float_t eta=-10., phi=-10.;
  Int_t id = GetId();
  const AliEMCALGeometry *g = AliEMCALGeometry::GetInstance();
  g->EtaPhiFromIndex(id,eta,phi);
  return eta ;
}

//____________________________________________________________________________
Float_t AliEMCALDigit::GetPhi() const
{ 
  //return phi coordinate of digit
  // should be change in EMCALGeometry - 19-nov-04
  Float_t eta=-10., phi=-10.;
  Int_t id = GetId();
  const AliEMCALGeometry *g = AliEMCALGeometry::GetInstance();
  g->EtaPhiFromIndex(id,eta,phi);
  return phi ;
}

//____________________________________________________________________________
Bool_t AliEMCALDigit::GetFALTROSample(const Int_t iSample, Int_t& timeBin, Int_t& amp) const
{
	//Get FALTRO sample in time bin iSample
	if (iSample >= fNSamples || iSample < 0 || fDigitType==kTrigger) return kFALSE;
	
	amp     =  fSamples[iSample] & 0xFFF;
	timeBin = (fSamples[iSample] >> 12) & 0xFF;
	
	return kTRUE;
}
//____________________________________________________________________________
Bool_t AliEMCALDigit::GetALTROSampleLG(const Int_t iSample, Int_t& timeBin, Int_t& amp) const
{
	//Get Low Gain ALTRO sample in time bin iSample
	if (iSample >= fNSamples || iSample < 0 || fDigitType==kLG) return kFALSE;
	
	amp     =  fSamples[iSample] & 0xFFF;
	timeBin = (fSamples[iSample] >> 12) & 0xFF;
	
	return kTRUE;
}

//____________________________________________________________________________
void AliEMCALDigit::SetALTROSamplesLG(const Int_t nSamples, Int_t *samples)
{
    //Set array of ALTRO samples, Low Gain or FALTRO 
	fNSamples = nSamples;
	fSamples  = new Int_t[fNSamples]; 
	for (Int_t i=0; i < fNSamples; i++) fSamples[i] = samples[i];
}

//____________________________________________________________________________
void AliEMCALDigit::SetALTROSamplesHG(const Int_t nSamples, Int_t *samples)
{
	//Set array of ALTRO samples, High Gain.
	fNSamplesHG = nSamples;
	fSamplesHG  = new Int_t[fNSamplesHG]; 
	for (Int_t i=0; i < fNSamplesHG; i++) fSamplesHG[i] = samples[i];
}


//____________________________________________________________________________
Bool_t AliEMCALDigit::GetALTROSampleHG(const Int_t iSample, Int_t& timeBin, Int_t& amp) const
{
	//Get High Gain ALTRO sample in time bin iSample
	if (iSample >= fNSamplesHG || iSample < 0 || fDigitType==kHG) return kFALSE;
	
	amp     =  fSamplesHG[iSample] & 0xFFF;
	timeBin = (fSamplesHG[iSample] >> 12) & 0xFF;
	
	return kTRUE;
}


//____________________________________________________________________________
Int_t AliEMCALDigit::GetPrimary(Int_t index) const
{
  // retrieves the primary particle number given its index in the list 
  if ( (index <= fNprimary) && (index > 0)){
    return fPrimary[index-1] ;
  } 

  return -1 ; 
}

//____________________________________________________________________________
Float_t AliEMCALDigit::GetDEPrimary(Int_t index) const
{
  // retrieves the primary particle energy contribution 
  // given its index in the list 
  if ( (index <= fNprimary) && (index > 0)){
    return fDEPrimary[index-1] ;
  } 

  return 0 ; 
  
}

//____________________________________________________________________________
Int_t AliEMCALDigit::GetIparent(Int_t index) const
{
  // retrieves the primary particle number given its index in the list 
  if ( index <= fNiparent && index > 0){
    return fIparent[index-1] ;
  } 

  return -1 ; 
  
}

//____________________________________________________________________________
Float_t AliEMCALDigit::GetDEParent(Int_t index) const
{
  // retrieves the parent particle energy contribution 
  // given its index in the list 
  if ( (index <= fNiparent) && (index > 0)){
    return fDEParent[index-1] ;
  } 

  return 0; 
}

//____________________________________________________________________________
void AliEMCALDigit::ShiftPrimary(Int_t shift){
  //shifts primary number to BIG offset, to separate primary in different TreeK
  Int_t index  ;
  for(index = 0; index <fNprimary; index++ ){
    fPrimary[index] = fPrimary[index]+ shift * 10000000   ;}
  for(index =0; index <fNiparent; index++){
    fIparent[index] = fIparent[index] + shift * 10000000 ;}
}


//____________________________________________________________________________
AliEMCALDigit& AliEMCALDigit::operator= (const AliEMCALDigit &digit)
{
  // assignment operator
  
  if(&digit == this) return *this;
  
  fAmpFloat    = digit.fAmpFloat;	
  fNSamples    = digit.fNSamples;
  fNSamplesHG  = digit.fNSamplesHG;
  fNprimary    = digit.fNprimary;
  fNMaxPrimary = digit.fNMaxPrimary;
  fNiparent    = digit.fNiparent;
  fNMaxiparent = digit.fNMaxiparent;
  fMaxIter     = digit.fMaxIter;
  fTime        = digit.fTime;
  fTimeR       = digit.fTimeR; 
  fChi2        = digit.fChi2; 
  fNDF         = digit.fNDF;
  fDigitType   = digit.fDigitType;
  fAmpCalib    = digit.fAmpCalib;
  fAmp         = digit.fAmp ;
  fId          = digit.fId;
  fIndexInList = digit.fIndexInList ; 
  
  // data members  
  if (fSamples )   delete [] fSamples ;   fSamples   = NULL ;
  if (fSamplesHG ) delete [] fSamplesHG ; fSamplesHG = NULL ;
  if (fPrimary  )  delete [] fPrimary ;   fPrimary   = NULL ;
  if (fDEPrimary)  delete [] fDEPrimary ; fDEPrimary = NULL ;
  if (fIparent )   delete [] fIparent ;   fIparent   = NULL ;
  if (fDEParent)   delete [] fDEParent ;  fDEParent  = NULL ;
  
  if (fNSamples){
    fSamples     = new Int_t[fNSamples]; 
    for (Int_t i=0; i < digit.fNSamples; i++) fSamples[i] = digit.fSamples[i];
  }
  
  if (fNSamplesHG){
    fSamplesHG   = new Int_t[fNSamplesHG]; 
    for (Int_t i=0; i < digit.fNSamplesHG; i++) fSamplesHG[i] = digit.fSamplesHG[i];
  }
  
  
  if (fNMaxPrimary){
    fPrimary   = new Int_t  [fNMaxPrimary] ;  
    fDEPrimary = new Float_t[fNMaxPrimary] ;
    for ( Int_t i = 0; i < fNMaxPrimary ; i++) {
      fPrimary[i]   = digit.fPrimary[i] ;
      fDEPrimary[i] = digit.fDEPrimary[i] ;
    }
  }
  
  if (fNMaxiparent){
    fIparent  = new Int_t  [fNMaxiparent] ;
    fDEParent = new Float_t[fNMaxiparent] ;
    for (Int_t j = 0; j< fNMaxiparent ; j++) {
      fIparent[j]  = digit.fIparent[j] ;
      fDEParent[j] = digit.fDEParent[j] ;
    }
  }  
  
  return *this;
  
}  

//____________________________________________________________________________
Bool_t AliEMCALDigit::operator==(AliEMCALDigit const & digit) const 
{
  // Two digits are equal if they have the same Id
  
  if ( fId == digit.fId ) 
    return kTRUE ;
  else 
    return kFALSE ;
}
 
//____________________________________________________________________________
AliEMCALDigit AliEMCALDigit::operator+(const AliEMCALDigit &digit) 
{
  // Adds the amplitude of digits and completes the list of primary particles
  // if amplitude is larger than 
  
  fAmpFloat += digit.fAmpFloat ;
  for (Int_t i=0; i < fNSamples  ; i++) fSamples[i]   += digit.fSamples[i];
  for (Int_t i=0; i < fNSamplesHG; i++) fSamplesHG[i] += digit.fSamplesHG[i];
  
  fAmp    += digit.fAmp ;
  if(fTime > digit.fTime)
    fTime = digit.fTime ;
  if (digit.fTimeR < fTimeR)
    fTimeR = digit.fTimeR ; 
  
  Int_t max1 = fNprimary ; 
  Int_t max2 = fNiparent ;  
  Int_t index ; 
  for (index = 0 ; index < digit.fNprimary  ; index++){
    Bool_t newPrim = kTRUE ;
    Int_t old ;
    for ( old = 0 ; (old < max1) && newPrim; old++) { //already have this primary?
      if(fPrimary[old] == digit.fPrimary[index]) {
        newPrim = kFALSE;
        fDEPrimary[old] += digit.fDEPrimary[index];
      }
    }
    if (newPrim) {
      if(max1<fNMaxPrimary){ 
        fPrimary[max1] = digit.fPrimary[index] ; 
        fDEPrimary[max1] = digit.fDEPrimary[index] ; 
        fNprimary++ ;
        max1++;
      }
      if(fNprimary==fNMaxPrimary) {
        
        TString mess = " NMaxPrimary  =  " ; 
        mess += fNMaxPrimary ; 
        mess += " is too small" ; 
        AliFatal(mess.Data()) ; 
        
      }
    }
  }
  
  for (index = 0 ; index < digit.fNiparent ; index++){
    Bool_t newParent = kTRUE ;
    Int_t old ;
    for ( old = 0 ; (old < max2) && newParent; old++) { //already have this primary?
      if(fIparent[old] == digit.fIparent[index]) {
        newParent = kFALSE;
        fDEParent[old] += digit.fDEParent[index];
      }
    }
    if(newParent){
      if(max2<fNMaxiparent) { 
        fIparent[max2] = digit.fIparent[index] ; 
        fDEParent[max2] = digit.fDEParent[index] ; 
        fNiparent++ ;
        max2++;
      }
      if(fNiparent==fNMaxiparent) {
        
        TString mess = " NMaxiparent  =  " ; 
        mess += fNMaxiparent ; 
        mess += " is too small" ; 
        AliFatal(mess.Data()) ; 
        
      }
    }
  }
  
  return *this ;
}

//____________________________________________________________________________
AliEMCALDigit AliEMCALDigit::operator*(Float_t factor) 
{
  // Multiplies the amplitude by a factor
  
  //Float_t tempo = static_cast<Float_t>(fAmp) ; 
  //tempo *= factor ; 
  //fAmp = static_cast<Int_t>(TMath::Floor(tempo)) ; 
	
  fAmpFloat *= factor;
  for (Int_t i=0; i < fNSamples  ; i++) fSamples[i]   = Int_t(factor*fSamples[i]);
  for (Int_t i=0; i < fNSamplesHG; i++) fSamplesHG[i] = Int_t(factor*fSamplesHG[i]);

  for(Int_t i=0; i < fNprimary; i++) 
    fDEPrimary[i] *= factor;
  for(Int_t i=0; i < fNiparent; i++) 
    fDEParent[i] *= factor;

  return *this ;
}

//____________________________________________________________________________
ostream& operator << ( ostream& out , const AliEMCALDigit & digit)
{
  // Prints the data of the digit
  
  out << "ID " << digit.fId << " Energy = " << digit.fAmp <<  " Time = " << digit.fTime << endl ; 
  for(Int_t i=0;i<digit.fNprimary;i++) 
    out << "Primary " << i+1 << " = " << digit.fPrimary[i] 
	<< " : DE " << digit.fDEPrimary[i] << endl ;
   
  for(Int_t j=0;j<digit.fNiparent;j++)
    out << "Iparent " << j+1 << " = " << digit.fIparent[j] 
	<< " : DE " << digit.fDEParent[j] << endl ;
  out << "Position in list = " << digit.fIndexInList << endl ; 
  return out ;
}

//____________________________________________________________________________
void AliEMCALDigit::Print(const Option_t* /*opt*/) const
{
    //Print
	printf("===\nDigit id: %4d / Energy %2.3f ; Time %e ; Time samples %d ; Chi2 %2.3f, NDF %d, Type? %d \n",
		   fId, fAmpFloat,fTime, fNSamples, fChi2, fNDF, fDigitType);
	if(fDigitType==kTrigger){
		printf("FALTRO: ");
		for (Int_t i=0; i < GetNFALTROSamples(); i++) 
		{
			Int_t timeBin, amp;
			GetFALTROSample(i, timeBin, amp);
			printf(" (%d,%d) ",timeBin,amp);
		}
		printf("\n");
	}//trigger
	else{
		printf("ALTRO, Low Gain: ");
		for (Int_t i=0; i < GetNALTROSamplesLG(); i++) 
		{
			Int_t timeBin, amp;
			GetALTROSampleLG(i, timeBin, amp);
			printf(" (%d,%d) ",timeBin,amp);
		}
		printf("\n");
		printf("ALTRO, High Gain: ");
		for (Int_t i=0; i < GetNALTROSamplesHG(); i++) 
		{
			Int_t timeBin, amp;
			GetALTROSampleHG(i, timeBin, amp);
			printf(" (%d,%d) ",timeBin,amp);
		}
		printf("\n");
	}//trigger
	
}



