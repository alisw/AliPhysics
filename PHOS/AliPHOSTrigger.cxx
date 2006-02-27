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
//  Digits are grouped in TRU's (16x28 crystals ordered fNTRUPhi x fNTRUZ). 
//  The algorithm searches all possible 4x4 cell combinations per each TRU, 
//  adding the digits amplitude and finding the maximum. Maximums are compared 
//  to triggers threshold and they are set. Thresholds need to be fixed. 
//  Usage:
//
//  //Inside the event loop
//  AliEMCALTrigger *tr = new AliEMCALTrigger();//Init Trigger
//  tr->SetL0MBPbPbThreshold(500);
//  tr->SetL0MBppThreshold(100);
//  tr->SetL1JetLowPtThreshold(1000);
//  tr->SetL1JetMediumPtThreshold(10000);
//  tr->SetL1JetHighPtThreshold(20000);
//  tr->Trigger(); //Execute Trigger
//  tr->Print(""); //Print result, with "deb" option all data members 
//  //are printed
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TMatrixD.h"

// --- ALIROOT system ---
#include "AliPHOSTrigger.h" 
#include "AliPHOSGeometry.h"
#include "AliPHOSGetter.h" 
#include "AliTriggerInput.h"


ClassImp(AliPHOSTrigger)

//______________________________________________________________________
AliPHOSTrigger::AliPHOSTrigger()
  : AliTriggerDetector(), fNTRU(8), fNTRUZ(2), fNTRUPhi(4), 
    fL0MBPbPbThreshold(500), fL0MBppThreshold(50), 
    fL1JetLowPtThreshold(1200), fL1JetHighPtThreshold(30000) 
{
  //ctor

   SetName("PHOS");
   CreateInputs();
   
   //Print("all") ; 
}



//____________________________________________________________________________
AliPHOSTrigger::AliPHOSTrigger(const AliPHOSTrigger & trig) 
  : AliTriggerDetector(trig) 
{

  // cpy ctor
  fNTRU                 = trig.fNTRU ; 
  fNTRUZ                = trig.fNTRUZ ; 
  fNTRUPhi              = trig.fNTRUPhi ; 
  fL0MBPbPbThreshold    = trig.fL0MBPbPbThreshold ; 
  fL0MBppThreshold      = trig.fL0MBppThreshold ; 
  fL1JetLowPtThreshold  = trig.fL1JetLowPtThreshold ;
  fL1JetHighPtThreshold = trig.fL1JetHighPtThreshold ;

}

//----------------------------------------------------------------------
void AliPHOSTrigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "PHOS_MB_PbPb_L0", "PHOS PbPb Minimum Bias L0",  0x01 ) );
   fInputs.AddLast( new AliTriggerInput( "PHOS_MB_pp_L0", "PHOS pp Minimum Bias L0",  0x02 ) );
   fInputs.AddLast( new AliTriggerInput( "PHOS_PbPb_JetHPt_L1", "PHOS PbPb Jet High Pt L1", 0x04 ) );
   fInputs.AddLast( new AliTriggerInput( "PHOS_PbPb_JetLPt_L1", "PHOS PbPb Jet Low Pt L1", 0x08 ) );
 
}

//____________________________________________________________________________
TClonesArray *  AliPHOSTrigger::FillTRU(const TClonesArray * digits,
					const AliPHOSGeometry * geom,
					const Int_t nModules, 
					const Int_t nCrystalsPhi, 
					const Int_t nCrystalsZ) const {

  //Orders digits ampitudes list in fNTRU TRUs (28x16 crystals) per module. 
  //Each TRU is a TMatrixD, and they are kept in TClonesArrays. The number of 
  //TRU in phi is fNTRUPhi, and the number of TRU in eta is fNTRUZ. 

  //Check data members
  
  if(fNTRUZ*fNTRUPhi != fNTRU)
    Error("FillTRU"," Wrong number of TRUS per Z or Phi");

  //Initilize variables
  //List of TRU matrices initialized to 0.

  TClonesArray * matrix = new TClonesArray("TMatrixD",1000);
  for(Int_t k = 0; k < fNTRU*nModules ; k++){
    TMatrixD  * trus = new TMatrixD(nCrystalsPhi,nCrystalsZ) ;
    for(Int_t i = 0; i < nCrystalsPhi; i++)
      for(Int_t j = 0; j < nCrystalsZ; j++)
	(*trus)(i,j) = 0.0;
    new((*matrix)[k]) TMatrixD(*trus) ;
  }
  
  AliPHOSDigit * dig ;
  
  //Declare variables
  Int_t relid[4] ; 
  Float_t amp = 0;

  //Digits loop to fill TRU matrices with amplitudes.

  for(Int_t idig = 0 ; idig < digits->GetEntriesFast() ; idig++){
    
    dig = static_cast<AliPHOSDigit *>(digits->At(idig)) ;
    amp = dig->GetAmp() ; //Energy of the digit (arbitrary units)	    
    geom->AbsToRelNumbering(dig->GetId(), relid) ;
    //Transform digit number into 4 numbers
    //relid[0] = module
    //relid[1] = EMC (0) or CPV (-1)
    //relid[2] = row <= 64 (fNPhi)
    //relid[3] = column <= 56 (fNZ)
    
    if(relid[1] == 0){//Not CPV, Only EMC digits
     
      //Check to which TRU in the supermodule belongs the cell. 
      //Supermodules are divided in a TRU matrix of dimension 
      //(fNTRUPhi,fNTRUZ).
      //Each TRU is a cell matrix of dimension (nCrystalsPhi,nCrystalsZ)
      
      //First calculate the row and column in the supermodule 
      //of the TRU to which the cell belongs.

      Int_t col   = (relid[3]-1)/nCrystalsZ+1; 
      Int_t row   = (relid[2]-1)/nCrystalsPhi+1; 
      Int_t itru  = col*row + (relid[0]-1)*fNTRU - 1;

      //Fill TRU matrix with cell values

      TMatrixD * trus = dynamic_cast<TMatrixD *>(matrix->At(itru)) ;

    //Calculate row and column of the cell inside the TRU with number itru

      Int_t nrow = (relid[2]-1) - (row-1) *  nCrystalsPhi;	
      Int_t ncol = (relid[3]-1) - (col-1) *  nCrystalsZ;
     
      (*trus)(nrow,ncol) = amp ;
    }
  }
  return matrix;
}

//____________________________________________________________________________
void AliPHOSTrigger::MakeSlidingCell(const TClonesArray * trus, 
				     const Int_t mod,  
				     const Int_t nCrystalsPhi,
				     const Int_t nCrystalsZ,
				     Float_t *ampmax){

  //Sums energy of all possible 4x4 crystals per each TRU. Fast signal 
  //in the experiment is given by 2x2 crystals, for this reason we loop 
  //inside the TRU crystals by 2. 
 
  Float_t amp = 0 ;

  for(Int_t i = 0 ; i < fNTRU ; i++)   
    ampmax[i] = 0 ;

  //Loop over all TRUS in a module
  for(Int_t itru = 0 + (mod - 1) * fNTRU ; itru < mod*fNTRU ; itru++){
    TMatrixD * tru = dynamic_cast<TMatrixD *>(trus->At(itru)) ;
    //Sliding 2x2 cell
    //ampmax[itru-(mod-1)*8]=0.0;       
    for(Int_t irow = 0 ; irow <  nCrystalsPhi; irow += 2){ 
      for(Int_t icol = 0 ; icol < nCrystalsZ ; icol += 2){
	amp = 0;
	if( (irow+3) < nCrystalsPhi && (icol+3) < nCrystalsZ){//Avoid exit the TRU
	  for(Int_t i = irow; i < irow + 4 ; i++){
	    for(Int_t j = icol; j < icol + 4 ; j++){
	      amp += (*tru)(i,j) ;
	    }
	  }
	}
	if(amp > ampmax[itru-(mod-1)*fNTRU])
	  ampmax[itru-(mod-1)*fNTRU] = amp ;
      }
    }
  }
}

//____________________________________________________________________________
void AliPHOSTrigger::Trigger() 
{

  //Main Method to select triggers.
  //Getter
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
 
  //Get Geometry
  const AliPHOSGeometry * geom = AliPHOSGetter::Instance()->PHOSGeometry() ;

  //Define some useful parameters
  Int_t nModules     = geom->GetNModules();
  Int_t nCrystalsPhi = geom->GetNPhi()/fNTRUPhi ;// 64/4=16
  Int_t nCrystalsZ   = geom->GetNZ()/fNTRUZ ;// 56/2=28

  //Take the digits list and declare digits pointers
  TClonesArray * digits = gime->Digits() ;

  //Fill TRU Matrix  
  TClonesArray * trus = FillTRU(digits,geom, nModules, nCrystalsPhi, 
				nCrystalsZ) ;

  //Do Cell Sliding and select Trigger
  Float_t max [8] ;
  for(Int_t imod = 1 ; imod <= nModules ; imod++) {
    MakeSlidingCell(trus, imod, nCrystalsPhi, nCrystalsZ, max);
//     cout<<" Max Amplitude: mod " << imod <<" TRU0 " << max[0] 
// 	<<" TRU1 " << max[1] <<" TRU2 " << max[2] <<" TRU3 " << max[3]
//   	<<" TRU4 " << max[4] <<" TRU5 " << max[5] <<" TRU6 " << max[6]
// 	<<" TRU7 " << max[7] <<endl;
  
    SetTriggers(max) ;
  }
  
}

//____________________________________________________________________________
void AliPHOSTrigger::Print(const Option_t * opt) const 
{

  //Prints main parameters
 
  if(! opt)
    return;
  AliTriggerInput* in = 0x0 ;

  AliInfo("PHOS trigger information:") ; 
  printf( "             Threshold for pp MB LO %d\n", fL0MBppThreshold) ;  
  in = (AliTriggerInput*)fInputs.FindObject( "PHOS_MB_pp_L0" );
  if(in->GetValue())
    printf( "             *** PHOS MB pp LO is set ***\n") ; 

  printf( "             Threshold for PbPb MB LO %d\n", fL0MBPbPbThreshold) ;  
  in = (AliTriggerInput*)fInputs.FindObject( "PHOS_MB_PbPb_L0" );
  if(in->GetValue())
    printf( "             *** PHOS MB PbPb LO is set ***\n") ; 
  
  printf( "             Jet Low Pt Threshold for PbPb L1 %d\n", fL1JetLowPtThreshold) ;
  in = (AliTriggerInput*)fInputs.FindObject( "PHOS_PbPb_JetLPt_L1" );
  if(in->GetValue())
    printf( "             *** PHOS Jet Low Pt for PbPb L1 is set ***\n") ;

  printf( "             Jet High Pt Threshold for L1 %d\n", fL1JetHighPtThreshold) ;  
  in = (AliTriggerInput*) fInputs.FindObject( "PHOS_PbPb_JetHPt_L1" );
  if(in->GetValue())
    printf( "             *** PHOS Jet High Pt for PbPb L1 is set ***\n") ;

  if(strstr(opt,"all")){
    printf( "                         Number of TRUs %d\n", fNTRU) ;
    printf( "                         Number of crystals in Z in TRUs %d\n", 
	    fNTRUZ) ;
    printf( "                         Number of crystals in Phi in TRUs %d\n", 
	    fNTRUPhi) ;
  }   
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
 
  if(max >= fL0MBppThreshold)
    SetInput("PHOS_MB_pp_L0");
  if(max >= fL0MBPbPbThreshold)
    SetInput("PHOS_MB_PbPb_L0");
  if(max >= fL1JetLowPtThreshold)
    SetInput("PHOS_PbPb_JetLPt_L1"); 
  if(max >= fL1JetHighPtThreshold)
    SetInput("PHOS_PbPb_JetHPt_L1");
  
}
