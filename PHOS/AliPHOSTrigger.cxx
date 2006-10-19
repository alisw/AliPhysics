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
//  Class for trigger analysis.
//  Digits are grouped in TRU's (Trigger Units). A TRU consist of 16x28 
//  crystals ordered fNTRUPhi x fNTRUZ. The algorithm searches all possible 
//  2x2 and nxn (n multiple of 4) crystal combinations per each TRU, adding the 
//  digits amplitude and  finding the maximum. Maxima are transformed in ADC 
//  time samples. Each time bin is compared to the trigger threshold until it is larger 
//  and then, triggers are set. Thresholds need to be fixed. 
//  Usage:
//
//  //Inside the event loop
//  AliEMCALTrigger *tr = new AliEMCALTrigger();//Init Trigger
//  tr->SetL0Threshold(100);
//  tr->SetL1JetLowPtThreshold(1000);
//  tr->SetL1JetHighPtThreshold(20000);
//  tr->Trigger(); //Execute Trigger
//  tr->Print(""); //Print result, with "deb" option all data members 
//  //are printed
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
//#include "TMatrixD.h"

// --- ALIROOT system ---
#include "AliPHOS.h"
#include "AliPHOSTrigger.h" 
#include "AliPHOSGeometry.h"
#include "AliPHOSGetter.h" 
#include "AliTriggerInput.h"
//#include "AliLog.h"

ClassImp(AliPHOSTrigger)

//______________________________________________________________________
AliPHOSTrigger::AliPHOSTrigger()
  : AliTriggerDetector(),
    f2x2MaxAmp(-1), f2x2CrystalPhi(-1),  f2x2CrystalEta(-1), f2x2SM(0),
    fnxnMaxAmp(-1), fnxnCrystalPhi(-1),  fnxnCrystalEta(-1), fnxnSM(0),
    fADCValuesHighnxn(0), fADCValuesLownxn(0),
    fADCValuesHigh2x2(0), fADCValuesLow2x2(0), fDigitsList(0),
    fL0Threshold(50), fL1JetLowPtThreshold(200), fL1JetHighPtThreshold(500),
    fNTRU(8), fNTRUZ(2), fNTRUPhi(4),  fPatchSize(1), fSimulation(kTRUE)
{
  //ctor
  fADCValuesHighnxn = 0x0; //new Int_t[fTimeBins];
  fADCValuesLownxn  = 0x0; //new Int_t[fTimeBins];
  fADCValuesHigh2x2 = 0x0; //new Int_t[fTimeBins];
  fADCValuesLow2x2  = 0x0; //new Int_t[fTimeBins];

  SetName("PHOS");
  CreateInputs();
  
  //Print("") ; 
}

//____________________________________________________________________________
AliPHOSTrigger::AliPHOSTrigger(const AliPHOSTrigger & trig) : 
  AliTriggerDetector(trig),
  f2x2MaxAmp(trig.f2x2MaxAmp),
  f2x2CrystalPhi(trig.f2x2CrystalPhi),
  f2x2CrystalEta(trig.f2x2CrystalEta),
  f2x2SM(trig.f2x2SM),
  fnxnMaxAmp(trig.fnxnMaxAmp),
  fnxnCrystalPhi(trig.fnxnCrystalPhi),
  fnxnCrystalEta(trig.fnxnCrystalEta),
  fnxnSM(trig.fnxnSM),
  fADCValuesHighnxn(trig.fADCValuesHighnxn),
  fADCValuesLownxn(trig.fADCValuesLownxn),
  fADCValuesHigh2x2(trig.fADCValuesHigh2x2),
  fADCValuesLow2x2(trig.fADCValuesLow2x2),
  fDigitsList(trig.fDigitsList),
  fL0Threshold(trig.fL0Threshold),
  fL1JetLowPtThreshold(trig.fL1JetLowPtThreshold),
  fL1JetHighPtThreshold(trig.fL1JetHighPtThreshold),
  fNTRU(trig.fNTRU),
  fNTRUZ(trig.fNTRUZ),
  fNTRUPhi(trig.fNTRUPhi),
  fPatchSize(trig.fPatchSize), fSimulation(trig.fSimulation)
{
  // cpy ctor
}

//_________________________________________________________________________
AliPHOSTrigger & AliPHOSTrigger::operator = (const AliPHOSTrigger &)
{
  Fatal("operator =", "no implemented");
  return *this;
}

void AliPHOSTrigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "PHOS_L0",       "PHOS L0", 0x02 ) );
   fInputs.AddLast( new AliTriggerInput( "PHOS_JetHPt_L1","PHOS Jet High Pt L1", 0x04 ) );
   fInputs.AddLast( new AliTriggerInput( "PHOS_JetLPt_L1","PHOS Jet Low Pt L1", 0x08 ) );
 
}

//____________________________________________________________________________
void AliPHOSTrigger::FillTRU(const TClonesArray * digits, const AliPHOSGeometry * geom, TClonesArray * ampmatrix, TClonesArray * timeRmatrix) const {

  //Orders digits ampitudes list and times in fNTRU TRUs (28x16 crystals) 
  //per module. Each TRU is a TMatrixD, and they are kept in TClonesArrays. 
  //In a module, the number of TRU in phi is fNTRUPhi, and the number of 
  //TRU in eta is fNTRUZ. 

  //Check data members
  
  if(fNTRUZ*fNTRUPhi != fNTRU)
    Error("FillTRU"," Wrong number of TRUS per Z or Phi");

  //Initilize and declare variables
  Int_t nModules     = geom->GetNModules();
  Int_t nCrystalsPhi = geom->GetNPhi()/fNTRUPhi ;// 64/4=16
  Int_t nCrystalsZ   = geom->GetNZ()/fNTRUZ ;// 56/2=28
  Int_t relid[4] ; 
  Float_t amp   = -1;
  Float_t timeR = -1;
  Int_t id      = -1;

  //List of TRU matrices initialized to 0.
  for(Int_t k = 0; k < fNTRU*nModules ; k++){
    TMatrixD  * amptrus   = new TMatrixD(nCrystalsPhi,nCrystalsZ) ;
    TMatrixD  * timeRtrus = new TMatrixD(nCrystalsPhi,nCrystalsZ) ;
    for(Int_t i = 0; i < nCrystalsPhi; i++){
      for(Int_t j = 0; j < nCrystalsZ; j++){
	(*amptrus)(i,j)   = 0.0;
	(*timeRtrus)(i,j) = 0.0;
      }
    }
    new((*ampmatrix)[k])   TMatrixD(*amptrus) ;
    new((*timeRmatrix)[k]) TMatrixD(*timeRtrus) ; 
  }
  
  AliPHOSDigit * dig ;
 
  //Digits loop to fill TRU matrices with amplitudes.

  for(Int_t idig = 0 ; idig < digits->GetEntriesFast() ; idig++){
    
    dig    = static_cast<AliPHOSDigit *>(digits->At(idig)) ;
    amp    = dig->GetEnergy() ;   // Energy of the digit 
    id     = dig->GetId() ;    // Id label of the cell
    timeR  = dig->GetTimeR() ; // Earliest time of the digit
    geom->AbsToRelNumbering(id, relid) ;
    //Transform digit number into 4 numbers
    //relid[0] = module
    //relid[1] = EMC (0) or CPV (-1)
    //relid[2] = row <= 64 (fNPhi)
    //relid[3] = column <= 56 (fNZ)
    
    if(relid[1] == 0){//Not CPV, Only EMC digits
     
      //Check to which TRU in the supermodule belongs the crystal. 
      //Supermodules are divided in a TRU matrix of dimension 
      //(fNTRUPhi,fNTRUZ).
      //Each TRU is a crystal matrix of dimension (nCrystalsPhi,nCrystalsZ)
      
      //First calculate the row and column in the supermodule 
      //of the TRU to which the crystal belongs.
      Int_t col   = (relid[3]-1)/nCrystalsZ+1; 
      Int_t row   = (relid[2]-1)/nCrystalsPhi+1;
 
      //Calculate label number of the TRU  
      Int_t itru  = (row-1) + (col-1)*fNTRUPhi + (relid[0]-1)*fNTRU ;

      //Fill TRU matrix with crystal values
      TMatrixD * amptrus   = dynamic_cast<TMatrixD *>(ampmatrix->At(itru)) ;
      TMatrixD * timeRtrus = dynamic_cast<TMatrixD *>(timeRmatrix->At(itru)) ;

      //Calculate row and column of the crystal inside the TRU with number itru
      Int_t irow = (relid[2]-1) - (row-1) *  nCrystalsPhi;	
      Int_t icol = (relid[3]-1) - (col-1) *  nCrystalsZ;
      
      (*amptrus)(irow,icol)   = amp ;
      (*timeRtrus)(irow,icol) = timeR ;
    }
  }
}

//______________________________________________________________________
void AliPHOSTrigger::GetCrystalPhiEtaIndexInModuleFromTRUIndex(const Int_t itru,const Int_t iphitru,const Int_t ietatru,Int_t &iphiMod,Int_t &ietaMod,const AliPHOSGeometry* geom) const 
{
  // This method transforms the (eta,phi) index of a crystals in a 
  // TRU matrix into Super Module (eta,phi) index.
  
  // Calculate in which row and column in which the TRU are 
  // ordered in the SM
  Int_t col = itru/ fNTRUPhi + 1;
  Int_t row = itru - (col-1)*fNTRUPhi + 1;
  
  //Calculate the (eta,phi) index in SM
  Int_t nCrystalsPhi = geom->GetNPhi()/fNTRUPhi;
  Int_t nCrystalsZ   = geom->GetNZ()/fNTRUZ;
  
  iphiMod = nCrystalsPhi*(row-1) + iphitru + 1 ;
  ietaMod = nCrystalsZ*(col-1) + ietatru + 1 ;

}
//____________________________________________________________________________
void AliPHOSTrigger::MakeSlidingCell(const TClonesArray * amptrus, const TClonesArray * timeRtrus, const Int_t imod, TMatrixD *ampmax2, TMatrixD *ampmaxn, const AliPHOSGeometry *geom){
  //Sums energy of all possible 2x2 (L0) and nxn (L1) crystals per each TRU. 
  //Fast signal in the experiment is given by 2x2 crystals, 
  //for this reason we loop inside the TRU crystals by 2. 
 
  //Declare and initialize varibles
  Int_t nCrystalsPhi = geom->GetNPhi()/fNTRUPhi ;// 64/4=16
  Int_t nCrystalsZ   = geom->GetNZ()/fNTRUZ ;// 56/2=28
  Float_t amp2 = 0 ;
  Float_t ampn = 0 ; 
  for(Int_t i = 0; i < 4; i++){
    for(Int_t j = 0; j < fNTRU; j++){
      (*ampmax2)(i,j) = -1;
      (*ampmaxn)(i,j) = -1;
    }
  }

  //Create matrix that will contain 2x2 amplitude sums
  //used to calculate the nxn sums
  TMatrixD  * tru2x2 = new TMatrixD(nCrystalsPhi/2,nCrystalsZ/2) ;
  for(Int_t i = 0; i < nCrystalsPhi/2; i++)
    for(Int_t j = 0; j < nCrystalsZ/2; j++)
      (*tru2x2)(i,j) = 0.0;
    
  //Loop over all TRUS in a module
  for(Int_t itru = 0 + (imod - 1) * fNTRU ; itru < imod*fNTRU ; itru++){
    TMatrixD * amptru   = dynamic_cast<TMatrixD *>(amptrus->At(itru)) ;
    TMatrixD * timeRtru = dynamic_cast<TMatrixD *>(timeRtrus->At(itru)) ;
    Int_t mtru = itru-(imod-1)*fNTRU ; //Number of TRU in Module
    
    //Sliding 2x2, add 2x2 amplitudes (NOT OVERLAP)
    for(Int_t irow = 0 ; irow <  nCrystalsPhi; irow += 2){ 
      for(Int_t icol = 0 ; icol < nCrystalsZ ; icol += 2){
	amp2 = (*amptru)(irow,icol)+(*amptru)(irow+1,icol)+
	  (*amptru)(irow,icol+1)+(*amptru)(irow+1,icol+1);
	
	//Fill new matrix with added 2x2 crystals for use in nxn sums
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

    //Sliding nxn, add nxn amplitudes (OVERLAP)
    if(fPatchSize > 0){
      for(Int_t irow = 0 ; irow <  nCrystalsPhi/2; irow++){ 
	for(Int_t icol = 0 ; icol < nCrystalsZ/2 ; icol++){
	  ampn = 0;
	  if( (irow+fPatchSize) < nCrystalsPhi/2 && (icol+fPatchSize) < nCrystalsZ/2){//Avoid exit the TRU
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
	  if( (rown+i) < nCrystalsPhi && (coln+j) < nCrystalsZ/2){//Avoid exit the TRU
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
void AliPHOSTrigger::Print(const Option_t * opt) const 
{

  //Prints main parameters
 
  if(! opt)
    return;
  AliTriggerInput* in = 0x0 ;

  printf( "             Maximum Amplitude after Sliding Crystal, \n") ; 
  printf( "               -2x2 crystals sum (not overlapped): %10.2f, in Super Module %d\n",
	  f2x2MaxAmp,f2x2SM) ; 
  printf( "               -2x2 from row %d to row %d and from column %d to column %d\n", f2x2CrystalPhi, f2x2CrystalPhi+2, f2x2CrystalEta, f2x2CrystalEta+2) ; 

  if(fPatchSize > 0){
    printf( "             Patch Size, n x n: %d x %d cells\n",4*fPatchSize, 4*fPatchSize);
    printf( "               -nxn crystals sum (overlapped)    : %10.2f, in Super Module %d\n",
	    fnxnMaxAmp,fnxnSM) ; 
    printf( "               -nxn from row %d to row %d and from column %d to column %d\n", fnxnCrystalPhi, fnxnCrystalPhi+4, fnxnCrystalEta, fnxnCrystalEta+4) ; 
  }
  printf( "             Threshold for LO %10.1f\n", 
	  fL0Threshold) ;  
  
  printf( "             Threshold for LO %10.2f\n", fL0Threshold) ;  
  in = (AliTriggerInput*)fInputs.FindObject( "PHOS_L0" );
  if(in->GetValue())
    printf( "             *** PHOS LO is set ***\n") ; 
  
  printf( "             Jet Low Pt Threshold for L1 %10.2f\n", fL1JetLowPtThreshold) ;
  in = (AliTriggerInput*)fInputs.FindObject( "PHOS_JetLPt_L1" );
  if(in->GetValue())
    printf( "             *** PHOS Jet Low Pt for L1 is set ***\n") ;
  
  printf( "             Jet High Pt Threshold for L1 %10.2f\n", fL1JetHighPtThreshold) ;  
  in = (AliTriggerInput*) fInputs.FindObject( "PHOS_JetHPt_L1" );
  if(in->GetValue())
    printf( "              *** PHOS Jet High Pt for L1 is set ***\n") ;
 
}

//____________________________________________________________________________
void AliPHOSTrigger::SetTriggers(const Int_t iMod, const TMatrixD * ampmax2, const TMatrixD * ampmaxn, const AliPHOSGeometry *geom)  
{
  //Checks the 2x2 and nxn maximum amplitude per each TRU and compares 
  //with the different L0 and L1 triggers thresholds

  //Initialize variables
  Float_t max2[] = {-1,-1,-1,-1} ;
  Float_t maxn[] = {-1,-1,-1,-1} ;
  Int_t   itru2  = -1 ;
  Int_t   itrun  = -1 ;


  //Find maximum summed amplitude of all the TRU 
  //in a Module
  for(Int_t i = 0 ; i < fNTRU ; i++){
    if(max2[0] < (*ampmax2)(0,i) ){
      max2[0] =  (*ampmax2)(0,i) ; // 2x2 summed max amplitude
      max2[1] =  (*ampmax2)(1,i) ; // corresponding phi position in TRU
      max2[2] =  (*ampmax2)(2,i) ; // corresponding eta position in TRU
      max2[3] =  (*ampmax2)(3,i) ; // corresponding most recent time
      itru2   = i ; // TRU number
    }
    if(maxn[0] < (*ampmaxn)(0,i) ){
      maxn[0] =  (*ampmaxn)(0,i) ; // nxn summed max amplitude
      maxn[1] =  (*ampmaxn)(1,i) ; // corresponding phi position in TRU
      maxn[2] =  (*ampmaxn)(2,i) ; // corresponding eta position in TRU
      maxn[3] =  (*ampmaxn)(3,i) ; // corresponding most recent time
      itrun   = i ; // TRU number
    }
  }
  
  //Set max amplitude if larger than in other Modules
  Float_t maxtimeR2 = -1 ;
  Float_t maxtimeRn = -1 ;
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ;
  AliPHOS * phos  = gime->PHOS();
  Int_t nTimeBins = phos->GetRawFormatTimeBins() ;
 
  //Set max 2x2 amplitude and select L0 trigger
  if(max2[0] > f2x2MaxAmp ){
    f2x2MaxAmp  = max2[0] ;
    f2x2SM      = iMod ;
    maxtimeR2   = max2[3] ;
    GetCrystalPhiEtaIndexInModuleFromTRUIndex(itru2,static_cast<Int_t>(max2[1]),static_cast<Int_t>(max2[2]),f2x2CrystalPhi,f2x2CrystalEta,geom) ;
    
    //Transform digit amplitude in Raw Samples
    fADCValuesLow2x2  = new Int_t[nTimeBins];
    fADCValuesHigh2x2 = new Int_t[nTimeBins];
    
    phos->RawSampledResponse(maxtimeR2, f2x2MaxAmp, fADCValuesHigh2x2, fADCValuesLow2x2) ; 
    
    //Set Trigger Inputs, compare ADC time bins until threshold is attained
    //Set L0
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHigh2x2[i] >= fL0Threshold          || fADCValuesLow2x2[i] >= fL0Threshold){
	SetInput("PHOS_L0") ;
	break;
      }
    }
//     for(Int_t i = 0 ; i < 256 ; i++)
//       if(fADCValuesLow2x2[i]!=0||fADCValuesHigh2x2[i]!=0)
// 	cout<< "2x2 Time Bin "<<i
// 	    <<"; 2x2 Low Gain  "<<fADCValuesLow2x2[i]
// 	    <<"; 2x2 High Gain "<<fADCValuesHigh2x2[i]<<endl;
  }

  //Set max nxn amplitude and select L1 triggers
  if(maxn[0] > fnxnMaxAmp ){
    fnxnMaxAmp  = maxn[0] ;
    fnxnSM      = iMod ;
    maxtimeRn   = maxn[3] ;
    GetCrystalPhiEtaIndexInModuleFromTRUIndex(itrun,static_cast<Int_t>(maxn[1]),static_cast<Int_t>(maxn[2]),fnxnCrystalPhi,fnxnCrystalEta,geom) ; 
    
    //Transform digit amplitude in Raw Samples
    fADCValuesHighnxn = new Int_t[nTimeBins];
    fADCValuesLownxn  = new Int_t[nTimeBins];
    phos->RawSampledResponse(maxtimeRn, fnxnMaxAmp, fADCValuesHighnxn, fADCValuesLownxn) ;
    
    //Set Trigger Inputs, compare ADC time bins until threshold is attained
    //SetL1 Low
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHighnxn[i] >= fL1JetLowPtThreshold  || fADCValuesLownxn[i] >= fL1JetLowPtThreshold){
	SetInput("PHOS_JetLPt_L1") ;
	break; 
      }
    }
    //SetL1 High
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHighnxn[i] >= fL1JetHighPtThreshold || fADCValuesLownxn[i] >= fL1JetHighPtThreshold){
	SetInput("PHOS_JetHPt_L1") ;
	break;
      }
    }
//     for(Int_t i = 0 ; i < 256 ; i++)
//       if(fADCValuesLownxn[i]!=0||fADCValuesHighnxn[i]!=0)
// 	cout<< "nxn Time Bin "<<i
// 	    <<"; nxn Low Gain  "<<fADCValuesLownxn[i]
// 	    <<"; nxn High Gain "<<fADCValuesHighnxn[i]<<endl;
  }
}

//____________________________________________________________________________
void AliPHOSTrigger::Trigger() 
{

  //Main Method to select triggers.
  AliRunLoader *rl = gAlice->GetRunLoader();
  //Getter
  AliPHOSGetter * gime = AliPHOSGetter::Instance( rl->GetFileName() ) ;
  //AliPHOSGetter * gime = AliPHOSGetter::Instance() ;

  //Get Geometry
  const AliPHOSGeometry * geom = AliPHOSGetter::Instance()->PHOSGeometry() ;
   
  //Define parameters
  Int_t nModules     = geom->GetNModules();

  //Intialize data members each time the trigger is called in event loop
  f2x2MaxAmp = -1; f2x2CrystalPhi = -1;  f2x2CrystalEta = -1;
  fnxnMaxAmp = -1; fnxnCrystalPhi = -1;  fnxnCrystalEta = -1;

  //Take the digits list if simulation
  if(fSimulation)
    fDigitsList = gime->Digits() ;
  
  if(!fDigitsList)
    AliFatal("Digits not found !") ;
  
  //Fill TRU Matrix  
  TClonesArray * amptrus   = new TClonesArray("TMatrixD",1000);
  TClonesArray * timeRtrus = new TClonesArray("TMatrixD",1000);
  FillTRU(fDigitsList,geom,amptrus, timeRtrus) ;

  //Do Crystal Sliding and select Trigger
  //Initialize varible that will contain maximum amplitudes and 
  //its corresponding cell position in eta and phi, and time.
  TMatrixD  * ampmax2 = new TMatrixD(4,fNTRU) ;
  TMatrixD  * ampmaxn = new TMatrixD(4,fNTRU) ;

  for(Int_t imod = 1 ; imod <= nModules ; imod++) {
    //Do 2x2 and nxn sums, select maximums. 
    MakeSlidingCell(amptrus, timeRtrus, imod, ampmax2, ampmaxn, geom);
    //Set the trigger
    SetTriggers(imod,ampmax2,ampmaxn, geom) ;
  }
}
