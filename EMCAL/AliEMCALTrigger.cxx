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
//  Digits are grouped in TRU's (384 cells ordered fNTRUPhi x fNTRUEta). 
//  The algorithm searches all possible 2x2 and nxn (n is a multiple of 4) cell 
//  combinations per each TRU,  adding the digits amplitude and finding the 
//  maximum. Maxima are compared to triggers threshold and they are set. 
//  Thresholds need to be fixed. Last 2 modules are half size in Phi, I considered 
//  that the number of TRU is maintained for the last modules but decision not taken. 
//  If different, then this must be changed. 
//  Usage:
//
//  //Inside the event loop
//  AliEMCALTrigger *tr = new AliEMCALTrigger();//Init Trigger
//  tr->SetL0Threshold(100); //Arbitrary threshold values
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
    fnxnMaxAmp(-1), fnxnCellPhi(-1),  fnxnCellEta(-1),
    fnxnSM(0),
    fADCValuesHighnxn(0x0),fADCValuesLownxn(0x0),
    fADCValuesHigh2x2(0x0),fADCValuesLow2x2(0x0),
    fDigitsList(0x0),
    fL0Threshold(100),fL1JetLowPtThreshold(200),
    fL1JetMediumPtThreshold(500), fL1JetHighPtThreshold(1000),
    fPatchSize(1), fSimulation(kTRUE)
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
    fnxnMaxAmp(trig.fnxnMaxAmp), 
    fnxnCellPhi(trig.fnxnCellPhi),  
    fnxnCellEta(trig.fnxnCellEta),
    fnxnSM(trig.fnxnSM),
    fADCValuesHighnxn(trig.fADCValuesHighnxn),
    fADCValuesLownxn(trig.fADCValuesLownxn),
    fADCValuesHigh2x2(trig.fADCValuesHigh2x2),
    fADCValuesLow2x2(trig.fADCValuesLow2x2),
    fDigitsList(trig.fDigitsList),
    fL0Threshold(trig.fL0Threshold),
    fL1JetLowPtThreshold(trig.fL1JetLowPtThreshold),
    fL1JetMediumPtThreshold(trig.fL1JetMediumPtThreshold), 
    fL1JetHighPtThreshold(trig.fL1JetHighPtThreshold),
    fPatchSize(trig.fPatchSize), fSimulation(trig.fSimulation)
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
void AliEMCALTrigger::MakeSlidingCell(const TClonesArray * amptrus, const TClonesArray * timeRtrus, const Int_t isupermod,TMatrixD *ampmax2, TMatrixD *ampmaxn, AliEMCALGeometry *geom){
  
  //Sums energy of all possible 2x2 (L0) and nxn (L1) cells per each TRU. 
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
  Float_t ampn = 0 ; 
  for(Int_t i = 0; i < 4; i++){
    for(Int_t j = 0; j < nTRU; j++){
      (*ampmax2)(i,j) = -1;
      (*ampmaxn)(i,j) = -1;
    }
  }

  //Create matrix that will contain 2x2 amplitude sums
  //used to calculate the nxn sums
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

	//Fill matrix with added 2x2 crystals for use in nxn sums
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

    //Sliding nxn, add nxn amplitudes  (OVERLAP)
    if(fPatchSize > 0){
      for(Int_t irow = 0 ; irow <  nCellsPhi/2; irow++){ 
	for(Int_t icol = 0 ; icol < nCellsEta/2 ; icol++){
	  ampn = 0;
	  if( (irow+fPatchSize) < nCellsPhi/2 && (icol+fPatchSize) < nCellsEta/2){//Avoid exit the TRU
	    for(Int_t i = 0 ; i <= fPatchSize ; i++)
	      for(Int_t j = 0 ; j <= fPatchSize ; j++)
		ampn += (*tru2x2)(irow+i,icol+j);
	    //Select nxn maximum sums to select L1 
	    if(ampn > (*ampmaxn)(0,mtru)){
	      (*ampmaxn)(0,mtru) = ampn ; 
	      (*ampmaxn)(1,mtru) = irow*2;
	      (*ampmaxn)(2,mtru) = icol*2;
	    }
	  }
	}
      }
      
      //Find most recent time in selected nxn cell
      (*ampmaxn)(3,mtru) = 1 ;
      Int_t rown =  static_cast <Int_t> ((*ampmaxn)(1,mtru));
      Int_t coln =  static_cast <Int_t> ((*ampmaxn)(2,mtru));
      for(Int_t i = 0; i<4*fPatchSize; i++){
	for(Int_t j = 0; j<4*fPatchSize; j++){
	  if( (rown+i) < nCellsPhi && (coln+j) < nCellsEta/2){//Avoid exit the TRU
	    if((*amptru)(rown+i,coln+j) > 0 &&  (*timeRtru)(rown+i,coln+j)> 0){
	      if((*timeRtru)(rown+i,coln+j) <  (*ampmaxn)(3,mtru)  )
		(*ampmaxn)(3,mtru) =  (*timeRtru)(rown+i,coln+j);
	    }
	  }
	}
      }
    }
    else {  
	(*ampmaxn)(0,mtru) =  (*ampmax2)(0,mtru); 
	(*ampmaxn)(1,mtru) =  (*ampmax2)(1,mtru);
	(*ampmaxn)(2,mtru) =  (*ampmax2)(2,mtru);
	(*ampmaxn)(3,mtru) =  (*ampmax2)(3,mtru);
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
 
  if(fPatchSize > 0){
    printf( "             Patch Size, n x n: %d x %d cells\n",4*fPatchSize, 4*fPatchSize);
    printf( "               -nxn cells sum (overlapped)    : %10.2f, in Super Module %d\n",
	    fnxnMaxAmp,fnxnSM) ; 
    printf( "               -nxn from row %d to row %d and from column %d to column %d\n", fnxnCellPhi, fnxnCellPhi+4*fPatchSize, fnxnCellEta, fnxnCellEta+4*fPatchSize) ; 
  }
 
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
				  const TMatrixD *ampmaxn, AliEMCALGeometry *geom)  
{

  //Checks the 2x2 and nxn maximum amplitude per each TRU and 
  //compares with the different L0 and L1 triggers thresholds
  Float_t max2[] = {-1,-1,-1,-1} ;
  Float_t maxn[] = {-1,-1,-1,-1} ;
  Int_t   itru2  = -1 ;
  Int_t   itrun  = -1 ;

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
      if(maxn[0] < (*ampmaxn)(0,i) ){
	maxn[0] =  (*ampmaxn)(0,i) ; // nxn summed max amplitude
	maxn[1] =  (*ampmaxn)(1,i) ; // corresponding phi position in TRU
	maxn[2] =  (*ampmaxn)(2,i) ; // corresponding eta position in TRU
	maxn[3] =  (*ampmaxn)(3,i) ; // corresponding most recent time
	itrun   = i ;
      }
    }

  //--------Set max amplitude if larger than in other Super Modules------------
  Float_t maxtimeR2 = -1 ;
  Float_t maxtimeRn = -1 ;
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
  
  //------------Set max of nxn amplitudes and select L1 trigger---------
  if(maxn[0] > fnxnMaxAmp ){
    fnxnMaxAmp  = maxn[0] ;
    fnxnSM      = iSM ;
    maxtimeRn   = maxn[3] ;
    geom->GetCellPhiEtaIndexInSModuleFromTRUIndex(itrun, 
						  static_cast<Int_t>(maxn[1]),
						  static_cast<Int_t>(maxn[2]),
						  fnxnCellPhi,fnxnCellEta) ; 
    //Transform digit amplitude in Raw Samples
    fADCValuesHighnxn = new Int_t[nTimeBins];
    fADCValuesLownxn  = new Int_t[nTimeBins];
    emcal->RawSampledResponse(maxtimeRn, fnxnMaxAmp, fADCValuesHighnxn, fADCValuesLownxn) ;
    
    //Set Trigger Inputs, compare ADC time bins until threshold is attained
    //SetL1 Low
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHighnxn[i] >= fL1JetLowPtThreshold  || fADCValuesLownxn[i] >= fL1JetLowPtThreshold){
	SetInput("EMCAL_JetLPt_L1") ;
	break; 
      }
    }
    
    //SetL1 Medium
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHighnxn[i] >= fL1JetMediumPtThreshold || fADCValuesLownxn[i] >= fL1JetMediumPtThreshold){
	SetInput("EMCAL_JetMPt_L1") ;
	break;
      }
    }
    
    //SetL1 High
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHighnxn[i] >= fL1JetHighPtThreshold || fADCValuesLownxn[i] >= fL1JetHighPtThreshold){
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
  fnxnMaxAmp = -1; fnxnCellPhi = -1;  fnxnCellEta = -1;

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
  TMatrixD  * ampmaxn = new TMatrixD(4,nTRU) ;

  for(Int_t iSM = 0 ; iSM < nSuperModules ; iSM++) {
    //Do 2x2 and nxn sums, select maximums. 
    MakeSlidingCell(amptrus, timeRtrus, iSM, ampmax2, ampmaxn, geom);
  
    //Set the trigger
    SetTriggers(iSM, ampmax2, ampmaxn, geom) ;
  }
}
