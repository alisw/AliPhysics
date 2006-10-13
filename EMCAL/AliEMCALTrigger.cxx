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
/* $Log $ */

//_________________________________________________________________________  
//
//  Class for trigger analysis.
//  Digits are grouped in TRU's (384 cells? ordered fNTRUPhi x fNTRUEta). 
//  The algorithm searches all possible 4x4 cell combinations per each TRU, 
//  adding the digits amplitude and finding the maximum. Maximums are compared 
//  to triggers threshold and they are set. Thresholds need to be fixed. 
//  Last 2 modules are half size in Phi, I considered that the number of TRU
//  is maintained for the last modules but decision not taken. If different, 
//  then this must be changed. 
//  Usage:
//
//  //Inside the event loop
//  AliEMCALTrigger *tr = new AliEMCALTrigger();//Init Trigger
//  tr->SetL0Threshold(100);
//  tr->SetL1JetLowPtThreshold(1000);
//  tr->SetL1JetMediumPtThreshold(10000);
//  tr->SetL1JetHighPtThreshold(20000);
//  tr->Trigger(); //Execute Trigger
//  tr->Print(""); //Print results
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "TMatrixD.h"

// --- ALIROOT system ---

#include "AliRun.h" 
#include "AliRunLoader.h" 
#include "AliTriggerInput.h"
#include "AliEMCAL.h" 
#include "AliEMCALLoader.h" 
#include "AliEMCALDigit.h"
#include "AliEMCALTrigger.h" 
#include "AliEMCALGeometry.h"

ClassImp(AliEMCALTrigger)

//______________________________________________________________________
AliEMCALTrigger::AliEMCALTrigger()
  : AliTriggerDetector(),
    f2x2MaxAmp(-1), f2x2CellPhi(-1),  f2x2CellEta(-1),
    f2x2SM(0),
    f4x4MaxAmp(-1), f4x4CellPhi(-1),  f4x4CellEta(-1),
    f4x4SM(0),
    fADCValuesHigh4x4(0x0),fADCValuesLow4x4(0x0),
    fADCValuesHigh2x2(0x0),fADCValuesLow2x2(0x0),
    fDigitsList(0x0),
    fL0Threshold(100),fL1JetLowPtThreshold(200),
    fL1JetMediumPtThreshold(500), fL1JetHighPtThreshold(1000),
    fSimulation(kTRUE)
{
  //ctor 

  SetName("EMCAL");
  CreateInputs();
   
   //Print("") ; 
}



//____________________________________________________________________________
AliEMCALTrigger::AliEMCALTrigger(const AliEMCALTrigger & trig) 
  : AliTriggerDetector(trig),
    f2x2MaxAmp(trig.f2x2MaxAmp), 
    f2x2CellPhi(trig.f2x2CellPhi),  
    f2x2CellEta(trig.f2x2CellEta),
    f2x2SM(trig.f2x2SM),
    f4x4MaxAmp(trig.f4x4MaxAmp), 
    f4x4CellPhi(trig.f4x4CellPhi),  
    f4x4CellEta(trig.f4x4CellEta),
    f4x4SM(trig.f4x4SM),
    fADCValuesHigh4x4(trig.fADCValuesHigh4x4),
    fADCValuesLow4x4(trig.fADCValuesLow4x4),
    fADCValuesHigh2x2(trig.fADCValuesHigh2x2),
    fADCValuesLow2x2(trig.fADCValuesLow2x2),
    fDigitsList(trig.fDigitsList),
    fL0Threshold(trig.fL0Threshold),
    fL1JetLowPtThreshold(trig.fL1JetLowPtThreshold),
    fL1JetMediumPtThreshold(trig.fL1JetMediumPtThreshold), 
    fL1JetHighPtThreshold(trig.fL1JetHighPtThreshold),
    fSimulation(trig.fSimulation)
{
  // cpy ctor
}

//----------------------------------------------------------------------
void AliEMCALTrigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "EMCAL_L0",       "EMCAL L0", 0x02 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_JetHPt_L1","EMCAL Jet High Pt L1",   0x04 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_JetMPt_L1","EMCAL Jet Medium Pt L1", 0x08 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_JetLPt_L1","EMCAL Jet Low Pt L1",    0x016 ) );
 
}

//____________________________________________________________________________
void AliEMCALTrigger::MakeSlidingCell(const TClonesArray * amptrus, const TClonesArray * timeRtrus, const Int_t isupermod,TMatrixD *ampmax2, TMatrixD *ampmax4, AliEMCALGeometry *geom){
  
  //Sums energy of all possible 2x2 (L0) and 4x4 (L1) cells per each TRU. 
  //Fast signal in the experiment is given by 2x2 cells, 
  //for this reason we loop inside the TRU cells by 2. 

  //Declare and initialize variables
  Int_t nCellsPhi  = geom->GetNPhi()*2/geom->GetNTRUPhi() ;
  if(isupermod > 9)
    nCellsPhi =  nCellsPhi / 2 ; //Half size SM. Not Final.
  // 12(tow)*2(cell)/1 TRU, cells in Phi in one TRU
  Int_t nCellsEta  = geom->GetNEta()*2/geom->GetNTRUEta() ;
  // 24(mod)*2(tower)/3 TRU, cells in Eta in one TRU
  Int_t nTRU          = geom->GetNTRU();//3 TRU per super module

  Float_t amp2 = 0 ;
  Float_t amp4 = 0 ; 
  for(Int_t i = 0; i < 4; i++){
    for(Int_t j = 0; j < nTRU; j++){
      (*ampmax2)(i,j) = -1;
      (*ampmax4)(i,j) = -1;
    }
  }

  //Create matrix that will contain 2x2 amplitude sums
  //used to calculate the 4x4 sums
  TMatrixD  * tru2x2 = new TMatrixD(nCellsPhi/2,nCellsEta/2) ;
  for(Int_t i = 0; i < nCellsPhi/2; i++)
    for(Int_t j = 0; j < nCellsEta/2; j++)
      (*tru2x2)(i,j) = -1;
  
  //Loop over all TRUS in a supermodule
  for(Int_t itru = 0 + isupermod * nTRU ; itru < (isupermod+1)*nTRU ; itru++) {
    TMatrixD * amptru   = dynamic_cast<TMatrixD *>(amptrus->At(itru)) ;
    TMatrixD * timeRtru = dynamic_cast<TMatrixD *>(timeRtrus->At(itru)) ;
    Int_t mtru = itru-isupermod*nTRU ; //Number of TRU in Supermodule
   
    //Sliding 2x2, add 2x2 amplitudes (NOT OVERLAP)
    for(Int_t irow = 0 ; irow <  nCellsPhi; irow += 2){ 
      for(Int_t icol = 0 ; icol < nCellsEta ; icol += 2){
	amp2 = (*amptru)(irow,icol)+(*amptru)(irow+1,icol)+
	  (*amptru)(irow,icol+1)+(*amptru)(irow+1,icol+1);

	//Fill matrix with added 2x2 crystals for use in 4x4 sums
	(*tru2x2)(irow/2,icol/2) = amp2 ;
	//Select 2x2 maximum sums to select L0 
	if(amp2 > (*ampmax2)(0,mtru)){
	  (*ampmax2)(0,mtru) = amp2 ; 
	  (*ampmax2)(1,mtru) = irow;
	  (*ampmax2)(2,mtru) = icol;
	}
      }
    }
    
    //Find most recent time in the selected 2x2 cell
    (*ampmax2)(3,mtru) = 1 ;
    Int_t row2 =  static_cast <Int_t> ((*ampmax2)(1,mtru));
    Int_t col2 =  static_cast <Int_t> ((*ampmax2)(2,mtru));
    for(Int_t i = 0; i<2; i++){
      for(Int_t j = 0; j<2; j++){
	if((*amptru)(row2+i,col2+j) > 0 &&  (*timeRtru)(row2+i,col2+j)> 0){	  
	  if((*timeRtru)(row2+i,col2+j) <  (*ampmax2)(3,mtru)  )
	    (*ampmax2)(3,mtru) =  (*timeRtru)(row2+i,col2+j);
	}
      }
    }
    
    //Sliding 4x4, add 4x4 amplitudes  (OVERLAP)
    for(Int_t irow = 0 ; irow <  nCellsPhi/2; irow++){ 
      for(Int_t icol = 0 ; icol < nCellsEta/2 ; icol++){
	if( (irow+1) < nCellsPhi/2 && (icol+1) < nCellsEta/2){//Avoid exit the TRU
	  amp4 = (*tru2x2)(irow,icol)+(*tru2x2)(irow+1,icol)+
	    (*tru2x2)(irow,icol+1)+(*tru2x2)(irow+1,icol+1);
	  //Select 4x4 maximum sums to select L1 
	  if(amp4 > (*ampmax4)(0,mtru)){
	    (*ampmax4)(0,mtru) = amp4 ; 
	    (*ampmax4)(1,mtru) = irow*2;
	    (*ampmax4)(2,mtru) = icol*2;
	  }
	}
      }
    }
    
    //Find most recent time in selected 4x4 cell
    (*ampmax4)(3,mtru) = 1 ;
    Int_t row4 =  static_cast <Int_t> ((*ampmax4)(1,mtru));
    Int_t col4 =  static_cast <Int_t> ((*ampmax4)(2,mtru));
    for(Int_t i = 0; i<4; i++){
      for(Int_t j = 0; j<4; j++){
	if((*amptru)(row4+i,col4+j) > 0 &&  (*timeRtru)(row4+i,col4+j)> 0){
	  if((*timeRtru)(row4+i,col4+j) <  (*ampmax4)(3,mtru)  )
	    (*ampmax4)(3,mtru) =  (*timeRtru)(row4+i,col4+j);
	}
      }
    }
  }
}

//____________________________________________________________________________
void AliEMCALTrigger::Print(const Option_t * opt) const 
{
  
  //Prints main parameters
  
  if(! opt)
    return;
  AliTriggerInput* in = 0x0 ;

  printf( "             Maximum Amplitude after Sliding Cell, \n") ; 
  printf( "               -2x2 cells sum (not overlapped): %10.2f, in Super Module %d\n",
	  f2x2MaxAmp,f2x2SM) ; 
   printf( "               -2x2 from row %d to row %d and from column %d to column %d\n", f2x2CellPhi, f2x2CellPhi+2, f2x2CellEta, f2x2CellEta+2) ; 
  printf( "               -4x4 cells sum (overlapped)    : %10.2f, in Super Module %d\n",
	  f4x4MaxAmp,f4x4SM) ; 
  printf( "               -4x4 from row %d to row %d and from column %d to column %d\n", f4x4CellPhi, f4x4CellPhi+4, f4x4CellEta, f4x4CellEta+4) ; 
  printf( "             Threshold for LO %10.2f\n", 
	  fL0Threshold) ;  
  in = (AliTriggerInput*)fInputs.FindObject( "EMCAL_L0" );
  if(in->GetValue())
    printf( "             *** EMCAL LO is set ***\n") ; 
  
  printf( "             Jet Low Pt Threshold for L1 %10.2f\n", 
	  fL1JetLowPtThreshold) ;
  in = (AliTriggerInput*)fInputs.FindObject( "EMCAL_JetLPt_L1" );
  if(in->GetValue())
    printf( "             *** EMCAL Jet Low Pt for L1 is set ***\n") ;

  printf( "             Jet Medium Pt Threshold for L1 %10.2f\n", 
	  fL1JetMediumPtThreshold) ;  
  in = (AliTriggerInput*) fInputs.FindObject( "EMCAL_JetMPt_L1" );
  if(in->GetValue())
    printf( "             *** EMCAL Jet Medium Pt for L1 is set ***\n") ;

  printf( "             Jet High Pt Threshold for L1 %10.2f\n", 
	  fL1JetHighPtThreshold) ;  
  in = (AliTriggerInput*) fInputs.FindObject( "EMCAL_JetHPt_L1" );
  if(in->GetValue())
    printf( "             *** EMCAL Jet High Pt for L1 is set ***\n") ;

}

//____________________________________________________________________________
void AliEMCALTrigger::SetTriggers(const Int_t iSM, const TMatrixD *ampmax2, 
				  const TMatrixD *ampmax4, AliEMCALGeometry *geom)  
{

  //Checks the 2x2 and 4x4 maximum amplitude per each TRU and 
  //compares with the different L0 and L1 triggers thresholds
  Float_t max2[] = {-1,-1,-1,-1} ;
  Float_t max4[] = {-1,-1,-1,-1} ;
  Int_t   itru2  = -1 ;
  Int_t   itru4  = -1 ;

  //Find maximum summed amplitude of all the TRU 
  //in a Super Module
    for(Int_t i = 0 ; i < geom->GetNTRU() ; i++){
      if(max2[0] < (*ampmax2)(0,i) ){
	max2[0] =  (*ampmax2)(0,i) ; // 2x2 summed max amplitude
	max2[1] =  (*ampmax2)(1,i) ; // corresponding phi position in TRU
	max2[2] =  (*ampmax2)(2,i) ; // corresponding eta position in TRU
	max2[3] =  (*ampmax2)(3,i) ; // corresponding most recent time
	itru2   = i ;
      }
      if(max4[0] < (*ampmax4)(0,i) ){
	max4[0] =  (*ampmax4)(0,i) ; // 4x4 summed max amplitude
	max4[1] =  (*ampmax4)(1,i) ; // corresponding phi position in TRU
	max4[2] =  (*ampmax4)(2,i) ; // corresponding eta position in TRU
	max4[3] =  (*ampmax4)(3,i) ; // corresponding most recent time
	itru4   = i ;
      }
    }

  //--------Set max amplitude if larger than in other Super Modules------------
  Float_t maxtimeR2 = -1 ;
  Float_t maxtimeR4 = -1 ;
  AliRunLoader *rl  = AliRunLoader::GetRunLoader();
  AliRun * gAlice   = rl->GetAliRun(); 
  AliEMCAL * emcal  = (AliEMCAL*)gAlice->GetDetector("EMCAL");
  Int_t nTimeBins = emcal->GetRawFormatTimeBins() ;

  //Set max of 2x2 amplitudes and select L0 trigger
  if(max2[0] > f2x2MaxAmp ){
    f2x2MaxAmp  = max2[0] ;
    f2x2SM      = iSM ;
    maxtimeR2   = max2[3] ;
    geom->GetCellPhiEtaIndexInSModuleFromTRUIndex(itru2, 
						  static_cast<Int_t>(max2[1]),
						  static_cast<Int_t>(max2[2]),
						  f2x2CellPhi,f2x2CellEta) ;
    
    //Transform digit amplitude in Raw Samples
    fADCValuesLow2x2  = new Int_t[nTimeBins];
    fADCValuesHigh2x2 = new Int_t[nTimeBins];
    emcal->RawSampledResponse(maxtimeR2, f2x2MaxAmp, fADCValuesHigh2x2, fADCValuesLow2x2) ; 
    
    //Set Trigger Inputs, compare ADC time bins until threshold is attained
    //Set L0
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHigh2x2[i] >= fL0Threshold          || fADCValuesLow2x2[i] >= fL0Threshold){
	SetInput("EMCAL_L0") ;
	break;
      }
    }
  }
  
  //------------Set max of 4x4 amplitudes and select L1 trigger---------
  if(max4[0] > f4x4MaxAmp ){
    f4x4MaxAmp  = max4[0] ;
    f4x4SM      = iSM ;
    maxtimeR4   = max4[3] ;
    geom->GetCellPhiEtaIndexInSModuleFromTRUIndex(itru4, 
						  static_cast<Int_t>(max4[1]),
						  static_cast<Int_t>(max4[2]),
						  f4x4CellPhi,f4x4CellEta) ; 
    //Transform digit amplitude in Raw Samples
    fADCValuesHigh4x4 = new Int_t[nTimeBins];
    fADCValuesLow4x4  = new Int_t[nTimeBins];
    emcal->RawSampledResponse(maxtimeR4, f4x4MaxAmp, fADCValuesHigh4x4, fADCValuesLow4x4) ;
    
    //Set Trigger Inputs, compare ADC time bins until threshold is attained
    //SetL1 Low
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHigh4x4[i] >= fL1JetLowPtThreshold  || fADCValuesLow4x4[i] >= fL1JetLowPtThreshold){
	SetInput("EMCAL_JetLPt_L1") ;
	break; 
      }
    }
    
    //SetL1 Medium
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHigh4x4[i] >= fL1JetMediumPtThreshold || fADCValuesLow4x4[i] >= fL1JetMediumPtThreshold){
	SetInput("EMCAL_JetMPt_L1") ;
	break;
      }
    }
    
    //SetL1 High
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHigh4x4[i] >= fL1JetHighPtThreshold || fADCValuesLow4x4[i] >= fL1JetHighPtThreshold){
	SetInput("EMCAL_JetHPt_L1") ;
	break;
      }
    }
  } 
}

//____________________________________________________________________________
void AliEMCALTrigger::Trigger() 
{
  //Main Method to select triggers.
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));

  //Load EMCAL Geometry
  rl->LoadgAlice(); 
  AliRun * gAlice = rl->GetAliRun(); 
  AliEMCAL * emcal  = (AliEMCAL*)gAlice->GetDetector("EMCAL");
  AliEMCALGeometry * geom = emcal->GetGeometry();
  
  if (geom==0)
    AliFatal("Did not get geometry from EMCALLoader");


  //Define parameters
  Int_t nSuperModules = geom->GetNumberOfSuperModules() ; //12 SM in EMCAL
  Int_t nTRU          = geom->GetNTRU();//3 TRU per super module

  //Intialize data members each time the trigger is called in event loop
  f2x2MaxAmp = -1; f2x2CellPhi = -1;  f2x2CellEta = -1;
  f4x4MaxAmp = -1; f4x4CellPhi = -1;  f4x4CellEta = -1;

  //Take the digits list if simulation
  if(fSimulation){
    rl->LoadDigits("EMCAL");
    fDigitsList = emcalLoader->Digits() ;
  }
  if(!fDigitsList)
    AliFatal("Digits not found !") ;
  
  //Take the digits list 
  
  //Fill TRU Matrix  
  TClonesArray * amptrus   = new TClonesArray("TMatrixD",1000);
  TClonesArray * timeRtrus = new TClonesArray("TMatrixD",1000);

  geom->FillTRU(fDigitsList, amptrus, timeRtrus) ;

  //Do Cell Sliding and select Trigger
  //Initialize varible that will contain maximum amplitudes and 
  //its corresponding cell position in eta and phi, and time.
  TMatrixD  * ampmax2 = new TMatrixD(4,nTRU) ;
  TMatrixD  * ampmax4 = new TMatrixD(4,nTRU) ;

  for(Int_t iSM = 0 ; iSM < nSuperModules ; iSM++) {
    //Do 2x2 and 4x4 sums, select maximums. 
    MakeSlidingCell(amptrus, timeRtrus, iSM, ampmax2, ampmax4, geom);
  
    //Set the trigger
    SetTriggers(iSM, ampmax2, ampmax4, geom) ;
  }
}
