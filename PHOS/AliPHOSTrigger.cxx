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
//
//  Class for trigger analysis.
//  Digits are grouped in TRU's (16x28 crystals). The algorithm searches all  
//  possible 4x4 crystal combinations and per each TRU, adding the digits 
//  amplitude and finding the maximum. Maximums are compared to triggers threshold 
//  and they are set.
//  FIRST ATTEMPT TO MAKE A TRIGGER CLASS. IT WILL CHANGE WHEN CENTRAL TRIGGER CLASS FIXES 
//  THINGS
// 
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---


// --- ALIROOT system ---
#include "AliPHOSTrigger.h" 
#include "AliPHOSGeometry.h"
#include "AliPHOSGetter.h" 


ClassImp(AliPHOSTrigger)

//____________________________________________________________________________
  AliPHOSTrigger::AliPHOSTrigger() : TObject(), 
				     fL0(kFALSE), fL1Low(kFALSE), fL1Medium(kFALSE), 
				     fL1High(kFALSE), fL0Threshold(50), fL1LowThreshold(1200),
				     fL1MediumThreshold(12000), fL1HighThreshold(30000)  
{
  // default ctor
  Print("") ; 
}

//____________________________________________________________________________
AliPHOSTrigger::AliPHOSTrigger(const AliPHOSTrigger & trig) : TObject(trig)
{

  // cpy ctor

  fL0                = trig.fL0 ;
  fL1Low             = trig.fL1Low ;
  fL1Medium          = trig.fL1Medium ;
  fL1High            = trig.fL1High ;
  fL0Threshold       = trig.fL0Threshold ; 
  fL1LowThreshold    = trig.fL1LowThreshold ;
  fL1MediumThreshold = trig.fL1MediumThreshold ;
  fL1HighThreshold   = trig.fL1HighThreshold ;

}

//____________________________________________________________________________
TClonesArray *  AliPHOSTrigger::FillTRU(const TClonesArray * digits){

  //Orders digits ampitudes list in 8 TRUs (16x28 crystals) per module. 
  //Each TRU is a TMatrixD, and they are kept in TClonesArrays.

  //Initilize variables

  const AliPHOSGeometry * geom = AliPHOSGetter::Instance()->PHOSGeometry() ;

  TClonesArray * matrix = new TClonesArray("TMatrixD",1000);
  for(Int_t k = 0; k < 40 ; k++){
    TMatrixD  * trus = new TMatrixD(16,28) ;
    for(Int_t i = 0; i < 16; i++)
      for(Int_t j = 0; j < 28; j++)
	(*trus)(i,j) = 0.0;
    new((*matrix)[k]) TMatrixD(*trus) ;
  }
  
  AliPHOSDigit * dig ;
  
  //Declare different variables
  Int_t relid[4] ; 
  Float_t amp = 0;
  
  for(Int_t idig = 0 ; idig < digits->GetEntriesFast() ; idig++){
    
    dig = static_cast<AliPHOSDigit *>(digits->At(idig)) ;
    amp = dig->GetAmp() ; //Energy of the digit (arbitrary units)	    
    geom->AbsToRelNumbering(dig->GetId(), relid) ;//Transform digit number into 4 numbers
    //relid[0] = module
    //relid[1] = EMC (0) or CPV (-1)
    //relid[2] = row <= 64 (fNPhi)
    //relid[3] = column <= 56 (fNZ)
    
    if(relid[1] == 0){
      
      Int_t ntru = 1;
      Int_t row  = 1;
      Int_t col  = 1;
      //Check which TRU in the module. It is divided in a (4,2) matrix.
      //Fill the TRU matrix (16,28)
      if(relid[3] > 28)
	col = 2;
      
      if(relid[2] > 16 && relid[2] <= 32)
	row = 2;
      
      if(relid[2] > 32 && relid[2] <= 48)
	row = 3;
      
      if(relid[2] > 48)
	row = 4;
      
      ntru     = col*row + (relid[0]-1)*8 - 1;
      TMatrixD * trus = dynamic_cast<TMatrixD *>(matrix->At(ntru)) ;
      
      Int_t nrow = (relid[2]-1) - (row-1) * 16 ;	
      Int_t ncol = (relid[3]-1) - (col-1) * 28 ;
     
      (*trus)(nrow,ncol) = amp ;
    }
  }
  return matrix;
}

//____________________________________________________________________________
void AliPHOSTrigger::InitTriggers()  
{
  //Initialize Boolean variables per each event

  fL0                = kFALSE ;
  fL1Low             = kFALSE ;
  fL1Medium          = kFALSE ;
  fL1High            = kFALSE ;

}
//____________________________________________________________________________
void AliPHOSTrigger::MakeSlidingCell(const TClonesArray * trus, const Int_t mod, 
				     Float_t *ampmax){

  //Sums energy of all possible 4x4 crystals per each TRU. Fast signal 
  //in the experiment is given by 2x2 crystals, for this reason we loop 
  //inside the TRU crystals by 2. 
 
  Float_t amp = 0 ;

  for(Int_t i = 0 ; i < 8 ; i++)   
    ampmax[i] = 0 ;

  //Loop over all TRUS in a module
  for(Int_t itru = 0 + (mod - 1) * 8 ; itru < mod*8 ; itru++){
    TMatrixD * tru = dynamic_cast<TMatrixD *>(trus->At(itru)) ;
    //Sliding 2x2 cell
    //ampmax[itru-(mod-1)*8]=0.0;       
    for(Int_t irow = 0 ; irow < 16 ; irow += 2){ 
      for(Int_t icol = 0 ; icol < 28 ; icol += 2){
	amp = 0;
	if( (irow+4) < 16 && (icol+4) < 28){//Avoid exit the TRU
	  for(Int_t i = irow; i < irow + 4 ; i++){
	    for(Int_t j = icol; j < icol + 4 ; j++){
	      amp += (*tru)(i,j) ;
	    }
	  }
	}
	if(amp > ampmax[itru-(mod-1)*8])
	  ampmax[itru-(mod-1)*8] = amp ;
      }
    }
  }
}

//____________________________________________________________________________
void AliPHOSTrigger::MakeTrigger() 
{

  //Main Method to select triggers.

  InitTriggers() ; //Initialize triggers to kFALSE

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  
  //Take the digits list and declare digits pointers
  TClonesArray * digits = gime->Digits() ;
  //Fill TRU Matrix  
  TClonesArray * trus = FillTRU(digits);

  //Do Cell Sliding and select Trigger
  Float_t max [8] ;
  for(Int_t imod = 1 ; imod <= 5 ; imod++) {
    MakeSlidingCell(trus, imod, max);
    SetTriggers(max) ;
  }
  
}

//____________________________________________________________________________
void AliPHOSTrigger::Print(const Option_t * opt) const 
{

  //Prints main parameters
 
  if(! opt)
    return;

    AliInfo("PHOS trigger information:") ; 
    printf( "                         Threshold for LO %d\n", fL0Threshold) ;  
    printf( "                     Low Threshold for L1 %d\n", fL1LowThreshold) ;  
    printf( "                  Medium Threshold for L1 %d\n", fL1MediumThreshold) ;  
    printf( "                    High Threshold for L1 %d\n", fL1HighThreshold) ;  
    if ( IsL0Set() ) 
      printf("                         LO is set\n") ; 
    else 
      printf("                         LO is not set\n") ;  
    if ( IsL1LowSet() ) 
      printf("                         L1 Low is set\n") ; 
    else 
      printf("                         L1 Low is not set\n") ;  
    if ( IsL1MediumSet() ) 
      printf("                         L1 Medium is set\n") ; 
    else 
      printf("                         L1 Medium is not set\n") ;  
    if ( IsL1HighSet() ) 
      printf("                         L1 High is set\n") ; 
    else 
      printf("                         L1 High is not set\n") ; 

}

//____________________________________________________________________________
void AliPHOSTrigger::SetTriggers(const Float_t * amp)  
{

  //Checks the maximum amplitude per each TRU and compares with the 
  //different triggers thresholds

  Float_t max = 0;
  for(Int_t i = 0 ; i < 8 ; i++){
    if(max < amp[i] )
      max = amp[i] ;
  }

  if(max >= fL0Threshold){
    SetL0();
    if(max >= fL1LowThreshold){
      SetL1Low(); 
      if(max >= fL1MediumThreshold){
	SetL1Medium(); 
	if(max >= fL1HighThreshold){
	  SetL1High();
	}
      }
    }
  }
}
