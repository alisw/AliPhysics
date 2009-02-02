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
//  Class for trigger analysis.
//  Digits are grouped in TRU's (Trigger Units). A TRU consist of 16x28 
//  crystals ordered fNTRUPhi x fNTRUZ. The algorithm searches all possible 
//  2x2 and nxn (n multiple of 2) crystal combinations per each TRU, adding the 
//  digits amplitude and  finding the maximum. If found, look if it is isolated.
//  Maxima are transformed in ADC time samples. Each time bin is compared to the trigger 
//  threshold until it is larger and then, triggers are set. Thresholds need to be fixed. 
//  Usage:
//
//  //Inside the event loop
//  AliPHOSTrigger *tr = new AliPHOSTrigger();//Init Trigger
//  tr->SetL0Threshold(100);
//  tr->SetL1JetLowPtThreshold(1000);
//  tr->SetL1JetMediumPtThreshold(10000);
//  tr->SetL1JetHighPtThreshold(20000);
//  ....
//  tr->Trigger(); //Execute Trigger
//  tr->Print(""); //Print data members after calculation.
//  
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TMath.h"

// --- ALIROOT system ---
#include "AliConfig.h"
#include "AliPHOS.h"
#include "AliPHOSTrigger.h" 
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h" 
#include "AliPHOSLoader.h" 
#include "AliPHOSPulseGenerator.h" 
#include "AliTriggerInput.h"


ClassImp(AliPHOSTrigger)

//______________________________________________________________________
AliPHOSTrigger::AliPHOSTrigger()
  : AliTriggerDetector(),
    f2x2MaxAmp(-1), f2x2CrystalPhi(-1),  f2x2CrystalEta(-1), f2x2SM(0),
    fnxnMaxAmp(-1), fnxnCrystalPhi(-1),  fnxnCrystalEta(-1), fnxnSM(0),
    fADCValuesHighnxn(0), fADCValuesLownxn(0),
    fADCValuesHigh2x2(0), fADCValuesLow2x2(0), fDigitsList(0),
    fAmptrus(0), fAmpmods(0), fTimeRtrus(0),
    fL0Threshold(50), fL1JetLowPtThreshold(200),   fL1JetMediumPtThreshold(500),  
    fL1JetHighPtThreshold(1000),
    fNTRU(8), fNTRUZ(2), fNTRUPhi(4), 
    fNCrystalsPhi(16),
    fNCrystalsZ(28),
    fPatchSize(1), fIsolPatchSize(1), 
    f2x2AmpOutOfPatch(-1), fnxnAmpOutOfPatch(-1), 
    f2x2AmpOutOfPatchThres(2),  fnxnAmpOutOfPatchThres(2), //2 GeV out of patch 
    fIs2x2Isol(kFALSE), fIsnxnIsol(kFALSE),  
    fSimulation(kTRUE), fIsolateInModule(kTRUE)
{
  //ctor
  fADCValuesHighnxn = 0x0; //new Int_t[fTimeBins];
  fADCValuesLownxn  = 0x0; //new Int_t[fTimeBins];
  fADCValuesHigh2x2 = 0x0; //new Int_t[fTimeBins];
  fADCValuesLow2x2  = 0x0; //new Int_t[fTimeBins];

  SetName("PHOS");
  CreateInputs();
  
  fAmptrus   = new TClonesArray("TMatrixD",1000);
  fAmpmods   = new TClonesArray("TMatrixD",1000);
  fTimeRtrus = new TClonesArray("TMatrixD",1000);
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
  fAmptrus(trig.fAmptrus), fAmpmods(trig.fAmpmods), fTimeRtrus(trig.fTimeRtrus),
  fL0Threshold(trig.fL0Threshold),
  fL1JetLowPtThreshold(trig.fL1JetLowPtThreshold),
  fL1JetMediumPtThreshold(trig.fL1JetMediumPtThreshold), 
  fL1JetHighPtThreshold(trig.fL1JetHighPtThreshold),
  fNTRU(trig.fNTRU),
  fNTRUZ(trig.fNTRUZ),
  fNTRUPhi(trig.fNTRUPhi),
  fNCrystalsPhi(trig.fNCrystalsPhi),
  fNCrystalsZ(trig. fNCrystalsZ),
  fPatchSize(trig.fPatchSize),
  fIsolPatchSize(trig.fIsolPatchSize), 
  f2x2AmpOutOfPatch(trig.f2x2AmpOutOfPatch), 
  fnxnAmpOutOfPatch(trig.fnxnAmpOutOfPatch), 
  f2x2AmpOutOfPatchThres(trig.f2x2AmpOutOfPatchThres),  
  fnxnAmpOutOfPatchThres(trig.fnxnAmpOutOfPatchThres), 
  fIs2x2Isol(trig.fIs2x2Isol),
  fIsnxnIsol(trig.fIsnxnIsol),  
  fSimulation(trig.fSimulation), 
  fIsolateInModule(trig.fIsolateInModule)
{
  // cpy ctor
}

//_________________________________________________________________________
AliPHOSTrigger::~AliPHOSTrigger() 
{
  // dtor
  
  if(fADCValuesHighnxn)delete [] fADCValuesHighnxn;
  if(fADCValuesLownxn)delete [] fADCValuesLownxn;
  if(fADCValuesHigh2x2)delete []  fADCValuesHigh2x2;
  if(fADCValuesLow2x2)delete [] fADCValuesLow2x2;
  // fDigitsList is now ours...
  if(fAmptrus)   { fAmptrus->Delete()  ; delete fAmptrus  ; }
  if(fAmpmods)   { fAmpmods->Delete()  ; delete fAmpmods  ; }
  if(fTimeRtrus) { fTimeRtrus->Delete(); delete fTimeRtrus; }
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

   TString name = GetName();
   
   fInputs.AddLast( new AliTriggerInput( "PHOS_L0",       name, 0 ) );
   fInputs.AddLast( new AliTriggerInput( "PHOS_JetHPt_L1",name, 1 ) );
   fInputs.AddLast( new AliTriggerInput( "PHOS_JetMPt_L1",name, 1 ) );
   fInputs.AddLast( new AliTriggerInput( "PHOS_JetLPt_L1",name, 1 ) );
 
}

//____________________________________________________________________________
void AliPHOSTrigger::FillTRU(const TClonesArray * digits, const AliPHOSGeometry * geom) const {

  //Orders digits ampitudes list and times in fNTRU TRUs (28x16 crystals) 
  //per module. Each TRU is a TMatrixD, and they are kept in TClonesArrays. 
  //In a module, the number of TRU in phi is fNTRUPhi, and the number of 
  //TRU in eta is fNTRUZ. Also fill a matrix with all amplitudes in module for isolation studies. 

  //Check data members
  
  if(fNTRUZ*fNTRUPhi != fNTRU)
    Error("FillTRU"," Wrong number of TRUS per Z or Phi");

  //Initilize and declare variables
  Int_t nModules     = geom->GetNModules();
  Int_t relid[4] ; 
  Float_t amp   = -1;
  Float_t timeR = -1;
  Int_t id      = -1;

  //List of TRU matrices initialized to 0.
  for(Int_t k = 0; k < fNTRU*nModules ; k++){
    TMatrixD   amptrus(fNCrystalsPhi,fNCrystalsZ) ;
    TMatrixD   timeRtrus(fNCrystalsPhi,fNCrystalsZ) ;
    for(Int_t i = 0; i < fNCrystalsPhi; i++){
      for(Int_t j = 0; j < fNCrystalsZ; j++){
	amptrus(i,j)   = 0.0;
	timeRtrus(i,j) = 0.0;
      }
    }
    new((*fAmptrus)[k])   TMatrixD(amptrus) ;
    new((*fTimeRtrus)[k]) TMatrixD(timeRtrus) ; 
  }

  //List of Modules matrices initialized to 0.
  Int_t nmodphi = geom->GetNPhi();
  Int_t nmodz = geom->GetNZ();
  
  for(Int_t k = 0; k < nModules ; k++){
    TMatrixD  ampmods(nmodphi,nmodz) ;
    for(Int_t i = 0; i < nmodphi; i++){
      for(Int_t j = 0; j < nmodz; j++){
	ampmods(i,j)   = 0.0;
      }
    }
    new((*fAmpmods)[k])   TMatrixD(ampmods) ;
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
      //############# TRU ###################
      //Check to which TRU in the supermodule belongs the crystal. 
      //Supermodules are divided in a TRU matrix of dimension 
      //(fNTRUPhi,fNTRUZ).
      //Each TRU is a crystal matrix of dimension (fNCrystalsPhi,fNCrystalsZ)
      
      //First calculate the row and column in the supermodule 
      //of the TRU to which the crystal belongs.
      Int_t col   = (relid[3]-1)/fNCrystalsZ+1; 
      Int_t row   = (relid[2]-1)/fNCrystalsPhi+1;
 
      //Calculate label number of the TRU  
      Int_t itru  = (row-1) + (col-1)*fNTRUPhi + (relid[0]-1)*fNTRU ;

      //Fill TRU matrix with crystal values
      TMatrixD * amptrus   = dynamic_cast<TMatrixD *>(fAmptrus  ->At(itru)) ;
      TMatrixD * timeRtrus = dynamic_cast<TMatrixD *>(fTimeRtrus->At(itru)) ;

      //Calculate row and column of the crystal inside the TRU with number itru
      Int_t irow = (relid[2]-1) - (row-1) *  fNCrystalsPhi;	
      Int_t icol = (relid[3]-1) - (col-1) *  fNCrystalsZ;
      
      (*amptrus)(irow,icol)   = amp ;
      (*timeRtrus)(irow,icol) = timeR ;

      //####################MODULE MATRIX ##################
      TMatrixD * ampmods   = dynamic_cast<TMatrixD *>(fAmpmods->At(relid[0]-1)) ;
      (*ampmods)(relid[2]-1,relid[3]-1)   = amp ;
    }
  }
}

//______________________________________________________________________
void AliPHOSTrigger::GetCrystalPhiEtaIndexInModuleFromTRUIndex(Int_t itru, Int_t iphitru, Int_t ietatru,Int_t &iphiMod,Int_t &ietaMod) const 
{
  // This method transforms the (eta,phi) index of a crystals in a 
  // TRU matrix into Super Module (eta,phi) index.
  
  // Calculate in which row and column in which the TRU are 
  // ordered in the SM
  Int_t col = itru/ fNTRUPhi + 1;
  Int_t row = itru - (col-1)*fNTRUPhi + 1;
  
  //Calculate the (eta,phi) index in SM
  
  iphiMod = fNCrystalsPhi*(row-1) + iphitru + 1 ;
  ietaMod = fNCrystalsZ*(col-1)   + ietatru + 1 ;

}

//____________________________________________________________________________
Bool_t AliPHOSTrigger::IsPatchIsolated(Int_t iPatchType, const Int_t imod, const Int_t mtru, const Float_t maxamp, const Int_t maxphi, const Int_t maxeta) {

  //Calculate if the maximum patch found is isolated, find amplitude around maximum (2x2 or nxn) patch, 
  //inside isolation patch . iPatchType = 0 means calculation for 2x2 patch, 
  //iPatchType = 1 means calculation for nxn patch.
  //In the next table there is an example of the different options of patch size and isolation patch size:
  //                                                                                 Patch Size (fPatchSize)
  //                                                             0                          1                                  2
  //          fIsolPatchSize                 2x2 (not overlap)   4x4 (overlapped)        6x6(overlapped) ...
  //                   1                                       4x4                      8x8                              10x10
  //                   2                                       6x6                     12x12                           14x14    
  //                   3                                       8x8                     16x16                           18x18
                          
  Bool_t b = kFALSE;
  Float_t amp = 0;
 
 //Get matrix of TRU or Module with maximum amplitude patch.
  Int_t itru = mtru+imod*fNTRU ; //number of tru, min 0 max 8*5.
  TMatrixD * ampmatrix   = 0x0;
  Int_t colborder = 0;
  Int_t rowborder = 0;

  if(fIsolateInModule){
    ampmatrix = dynamic_cast<TMatrixD *>(fAmpmods->At(imod)) ;
    rowborder = fNCrystalsPhi*fNTRUPhi;
    colborder = fNCrystalsZ*fNTRUZ;
    AliDebug(2,"Isolate trigger in Module");
  }
  else{
    ampmatrix = dynamic_cast<TMatrixD *>(fAmptrus->At(itru)) ;
    rowborder = fNCrystalsPhi;
    colborder = fNCrystalsZ;
    AliDebug(2,"Isolate trigger in TRU");
  }

  //Define patch cells
  Int_t isolcells = fIsolPatchSize*(1+iPatchType);
  Int_t ipatchcells = 2*(1+fPatchSize*iPatchType);
  Int_t minrow =  maxphi - isolcells;
  Int_t mincol =  maxeta - isolcells;
  Int_t maxrow =  maxphi + isolcells + ipatchcells;
  Int_t maxcol = maxeta +  isolcells + ipatchcells;

  AliDebug(2,Form("Number of added Isol Cells %d, Patch Size %d",isolcells, ipatchcells));
  AliDebug(2,Form("Patch: minrow %d, maxrow %d, mincol %d, maxcol %d",minrow,maxrow,mincol,maxcol));
  
  if(minrow < 0 || mincol < 0 || maxrow > rowborder || maxcol > colborder){
    AliDebug(1,Form("Out of Module/TRU range, cannot isolate patch"));
    return kFALSE;
  }

  //Add amplitudes in all isolation patch
  for(Int_t irow = minrow ; irow <  maxrow; irow ++)
    for(Int_t icol = mincol ; icol < maxcol ; icol ++)
      amp += (*ampmatrix)(irow,icol);

  AliDebug(2,Form("Type %d, Maximum amplitude %f, patch+isol square %f",iPatchType, maxamp, amp));

  if(TMath::Nint(amp*1E5) < TMath::Nint(maxamp*1E5)){
    AliError(Form("Bad sum: Type %d, Maximum amplitude %f, patch+isol square %f",iPatchType, maxamp, amp));
    return kFALSE;
  }
  else
    amp-=maxamp; //Calculate energy in isolation patch that do not comes from maximum patch.
  
  AliDebug(2, Form("Maximum amplitude %f, Out of patch %f",maxamp, amp));

  //Fill isolation amplitude data member and say if patch is isolated.
  if(iPatchType == 0){ //2x2 case
    f2x2AmpOutOfPatch = amp;   
    if(amp < f2x2AmpOutOfPatchThres)
      b=kTRUE;
  }
  else  if(iPatchType == 1){ //nxn case
    fnxnAmpOutOfPatch = amp;   
    if(amp < fnxnAmpOutOfPatchThres)
      b=kTRUE;
  }

  return b;

}


//____________________________________________________________________________
void AliPHOSTrigger::MakeSlidingCell(const Int_t imod, TMatrixD &ampmax2, TMatrixD &ampmaxn)
{
  //Sums energy of all possible 2x2 (L0) and nxn (L1) crystals per each TRU. 
  //Fast signal in the experiment is given by 2x2 crystals, 
  //for this reason we loop inside the TRU crystals by 2. 
 
  //Declare and initialize varibles
  Float_t amp2 = 0 ;
  Float_t ampn = 0 ; 
  for(Int_t i = 0; i < 4; i++){
    for(Int_t j = 0; j < fNTRU; j++){
      ampmax2(i,j) = -1;
      ampmaxn(i,j) = -1;
    }
  }

  //Create matrix that will contain 2x2 amplitude sums
  //used to calculate the nxn sums
  TMatrixD  tru2x2(fNCrystalsPhi/2,fNCrystalsZ/2) ;
  for(Int_t i = 0; i < fNCrystalsPhi/2; i++)
    for(Int_t j = 0; j < fNCrystalsZ/2; j++)
      tru2x2(i,j) = -1.;
    
  //Loop over all TRUS in a module
  for(Int_t itru = 0 + imod  * fNTRU ; itru < (imod+1)*fNTRU ; itru++){
    TMatrixD * amptru   = dynamic_cast<TMatrixD *>(fAmptrus  ->At(itru)) ;
    TMatrixD * timeRtru = dynamic_cast<TMatrixD *>(fTimeRtrus->At(itru)) ;
    Int_t mtru = itru-imod*fNTRU ; //Number of TRU in Module
    
    //Sliding 2x2, add 2x2 amplitudes (NOT OVERLAP)
    for(Int_t irow = 0 ; irow <  fNCrystalsPhi; irow += 2){ 
      for(Int_t icol = 0 ; icol < fNCrystalsZ ; icol += 2){
	amp2 = (*amptru)(irow,icol)+(*amptru)(irow+1,icol)+
	  (*amptru)(irow,icol+1)+(*amptru)(irow+1,icol+1);
	//Fill new matrix with added 2x2 crystals for use in nxn sums
	tru2x2(irow/2,icol/2) = amp2 ;
	//Select 2x2 maximum sums to select L0 
	if(amp2 > ampmax2(0,mtru)){
	  ampmax2(0,mtru) = amp2 ; 
	  ampmax2(1,mtru) = irow;
	  ampmax2(2,mtru) = icol;
	}
      }
    }

    //Find most recent time in the selected 2x2 cell
    ampmax2(3,mtru) = 1 ;
    Int_t row2 =  static_cast <Int_t> (ampmax2(1,mtru));
    Int_t col2 =  static_cast <Int_t> (ampmax2(2,mtru));
    for(Int_t i = 0; i<2; i++){
      for(Int_t j = 0; j<2; j++){
	if((*amptru)(row2+i,col2+j) > 0 &&  (*timeRtru)(row2+i,col2+j)> 0){	  
	  if((*timeRtru)(row2+i,col2+j) <  ampmax2(3,mtru)  )
	    ampmax2(3,mtru) =  (*timeRtru)(row2+i,col2+j);
	}
      }
    }

    //Sliding nxn, add nxn amplitudes (OVERLAP)
    if(fPatchSize > 0){
      for(Int_t irow = 0 ; irow <  fNCrystalsPhi/2; irow++){ 
	for(Int_t icol = 0 ; icol < fNCrystalsZ/2 ; icol++){
	  ampn = 0;
	  if( (irow+fPatchSize) < fNCrystalsPhi/2 && (icol+fPatchSize) < fNCrystalsZ/2){//Avoid exit the TRU
	    for(Int_t i = 0 ; i <= fPatchSize ; i++)
	      for(Int_t j = 0 ; j <= fPatchSize ; j++)
		ampn += tru2x2(irow+i,icol+j);
	    //Select nxn maximum sums to select L1 
	    if(ampn > ampmaxn(0,mtru)){
	      ampmaxn(0,mtru) = ampn ; 
	      ampmaxn(1,mtru) = irow*2;
	      ampmaxn(2,mtru) = icol*2;
	    }
	  }
	}
      }
      
      //Find most recent time in selected nxn cell
      ampmaxn(3,mtru) = 1 ;
      Int_t rown =  static_cast <Int_t> (ampmaxn(1,mtru));
      Int_t coln =  static_cast <Int_t> (ampmaxn(2,mtru));
      for(Int_t i = 0; i<4*fPatchSize; i++){
	for(Int_t j = 0; j<4*fPatchSize; j++){
	  if( (rown+i) < fNCrystalsPhi && (coln+j) < fNCrystalsZ/2){//Avoid exit the TRU
	    if((*amptru)(rown+i,coln+j) > 0 &&  (*timeRtru)(rown+i,coln+j)> 0){
	      if((*timeRtru)(rown+i,coln+j) <  ampmaxn(3,mtru)  )
		ampmaxn(3,mtru) =  (*timeRtru)(rown+i,coln+j);
	    }
	  }
	}
      }
    }
    else {  
	ampmaxn(0,mtru) =  ampmax2(0,mtru); 
	ampmaxn(1,mtru) =  ampmax2(1,mtru);
	ampmaxn(2,mtru) =  ampmax2(2,mtru);
	ampmaxn(3,mtru) =  ampmax2(3,mtru);
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
  printf( "               -2x2 Isolation Patch %d x %d, Amplitude out of 2x2 patch is %f, threshold %f, Isolated? %d \n", 
  	  2*fIsolPatchSize+2, 2*fIsolPatchSize+2,  f2x2AmpOutOfPatch,  f2x2AmpOutOfPatchThres,static_cast<Int_t> (fIs2x2Isol)) ; 
  if(fPatchSize > 0){
    printf( "             Patch Size, n x n: %d x %d cells\n",2*(fPatchSize+1), 2*(fPatchSize+1));
    printf( "               -nxn crystals sum (overlapped)    : %10.2f, in Super Module %d\n",
	    fnxnMaxAmp,fnxnSM) ; 
    printf( "               -nxn from row %d to row %d and from column %d to column %d\n", fnxnCrystalPhi, fnxnCrystalPhi+4*fPatchSize, fnxnCrystalEta, fnxnCrystalEta+4*fPatchSize) ; 
    printf( "               -nxn Isolation Patch %d x %d, Amplitude out of nxn patch is %f, threshold %f, Isolated? %d \n", 
	    4*fIsolPatchSize+2*(fPatchSize+1),4*fIsolPatchSize+2*(fPatchSize+1) ,  fnxnAmpOutOfPatch,  fnxnAmpOutOfPatchThres,static_cast<Int_t> (fIsnxnIsol) ) ; 
  }

  printf( "             Isolate in Module? %d\n",  
          fIsolateInModule) ;  

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
  
  printf( "             Jet Medium Pt Threshold for L1 %10.2f\n", fL1JetMediumPtThreshold) ;
  in = (AliTriggerInput*)fInputs.FindObject( "PHOS_JetMPt_L1" );
  if(in->GetValue())
    printf( "             *** PHOS Jet Medium Pt for L1 is set ***\n") ;
  
  printf( "             Jet High Pt Threshold for L1 %10.2f\n", fL1JetHighPtThreshold) ;  
  in = (AliTriggerInput*) fInputs.FindObject( "PHOS_JetHPt_L1" );
  if(in->GetValue())
    printf( "              *** PHOS Jet High Pt for L1 is set ***\n") ;
  
}

//____________________________________________________________________________
void AliPHOSTrigger::SetTriggers(const Int_t iMod, const TMatrixD & ampmax2, const TMatrixD & ampmaxn)  
{
  //Checks the 2x2 and nxn maximum amplitude per each TRU and compares 
  //with the different L0 and L1 triggers thresholds. It finds if maximum amplitudes are isolated.

  //Initialize variables
  Float_t max2[] = {-1,-1,-1,-1} ;
  Float_t maxn[] = {-1,-1,-1,-1} ;
  Int_t   mtru2  = -1 ;
  Int_t   mtrun  = -1 ;


  //Find maximum summed amplitude of all the TRU 
  //in a Module
  for(Int_t i = 0 ; i < fNTRU ; i++){
    if(max2[0] < ampmax2(0,i) ){
      max2[0] =  ampmax2(0,i) ; // 2x2 summed max amplitude
      max2[1] =  ampmax2(1,i) ; // corresponding phi position in TRU
      max2[2] =  ampmax2(2,i) ; // corresponding eta position in TRU
      max2[3] =  ampmax2(3,i) ; // corresponding most recent time
      mtru2   = i ; // TRU number in module
    }
    if(maxn[0] < ampmaxn(0,i) ){
      maxn[0] =  ampmaxn(0,i) ; // nxn summed max amplitude
      maxn[1] =  ampmaxn(1,i) ; // corresponding phi position in TRU
      maxn[2] =  ampmaxn(2,i) ; // corresponding eta position in TRU
      maxn[3] =  ampmaxn(3,i) ; // corresponding most recent time
      mtrun   = i ; // TRU number in module
    }
  }
  
  //Set max amplitude if larger than in other Modules
  Float_t maxtimeR2 = -1 ;
  Float_t maxtimeRn = -1 ;
  // Create a shaper pulse object
  AliPHOSPulseGenerator pulse ;
  Int_t nTimeBins = pulse.GetRawFormatTimeBins() ;
 
  //Set max 2x2 amplitude and select L0 trigger
  if(max2[0] > f2x2MaxAmp ){
    f2x2MaxAmp  = max2[0] ;
    f2x2SM      = iMod ;
    maxtimeR2   = max2[3] ;
    GetCrystalPhiEtaIndexInModuleFromTRUIndex(mtru2,
					      static_cast<Int_t>(max2[1]),
					      static_cast<Int_t>(max2[2]),
					      f2x2CrystalPhi,f2x2CrystalEta) ;
    
    //Isolated patch?
    if(fIsolateInModule)
      fIs2x2Isol =  IsPatchIsolated(0, iMod, mtru2,  f2x2MaxAmp, f2x2CrystalPhi,f2x2CrystalEta) ;
    else
      fIs2x2Isol =  IsPatchIsolated(0, iMod, mtru2,  f2x2MaxAmp,  static_cast<Int_t>(max2[1]), static_cast<Int_t>(max2[2])) ;

    //Transform digit amplitude in Raw Samples
    if (fADCValuesLow2x2 == 0) {
      fADCValuesLow2x2  = new Int_t[nTimeBins];
    }
    if(!fADCValuesHigh2x2) fADCValuesHigh2x2 = new Int_t[nTimeBins];

    
    pulse.SetAmplitude(f2x2MaxAmp);
    pulse.SetTZero(maxtimeR2);
    pulse.MakeSamples();
    pulse.GetSamples(fADCValuesHigh2x2, fADCValuesLow2x2) ; 
    
    //Set Trigger Inputs, compare ADC time bins until threshold is attained
    //Set L0
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHigh2x2[i] >= fL0Threshold || fADCValuesLow2x2[i] >= fL0Threshold) {
	SetInput("PHOS_L0") ;
	break;
      }
    }
  }

  //Set max nxn amplitude and select L1 triggers
  if(maxn[0] > fnxnMaxAmp  && fPatchSize > 0){
    fnxnMaxAmp  = maxn[0] ;
    fnxnSM      = iMod ;
    maxtimeRn   = maxn[3] ;
    GetCrystalPhiEtaIndexInModuleFromTRUIndex(mtrun,
					      static_cast<Int_t>(maxn[1]),
					      static_cast<Int_t>(maxn[2]),
					      fnxnCrystalPhi,fnxnCrystalEta) ; 
    
    //Isolated patch?
    if(fIsolateInModule)
      fIsnxnIsol =  IsPatchIsolated(1, iMod, mtrun,  fnxnMaxAmp, fnxnCrystalPhi, fnxnCrystalEta) ;
    else
      fIsnxnIsol =  IsPatchIsolated(1, iMod, mtrun,  fnxnMaxAmp,  static_cast<Int_t>(maxn[1]), static_cast<Int_t>(maxn[2])) ;

    //Transform digit amplitude in Raw Samples
    if (fADCValuesHighnxn == 0) {
      fADCValuesHighnxn = new Int_t[nTimeBins];
      fADCValuesLownxn  = new Int_t[nTimeBins];
    }

    pulse.SetAmplitude(fnxnMaxAmp);
    pulse.SetTZero(maxtimeRn);
    pulse.MakeSamples();
    pulse.GetSamples(fADCValuesHighnxn, fADCValuesLownxn) ;
    
    //Set Trigger Inputs, compare ADC time bins until threshold is attained
    //SetL1 Low
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHighnxn[i] >= fL1JetLowPtThreshold  || fADCValuesLownxn[i] >= fL1JetLowPtThreshold){
	SetInput("PHOS_JetLPt_L1") ;
	break; 
      }
    }
    //SetL1 Medium
    for(Int_t i = 0 ; i < nTimeBins ; i++){
      if(fADCValuesHighnxn[i] >= fL1JetMediumPtThreshold  || fADCValuesLownxn[i] >= fL1JetMediumPtThreshold){
	SetInput("PHOS_JetMPt_L1") ;
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
  }
}

//____________________________________________________________________________
void AliPHOSTrigger::Trigger(TClonesArray *digits) 
{
  //Main Method to select triggers.

  fDigitsList = digits;
  DoIt() ; 
}

//____________________________________________________________________________
void AliPHOSTrigger::DoIt()
{
  // does the trigger job

  AliRunLoader* rl = AliRunLoader::Instance() ;
  AliPHOSLoader * phosLoader = dynamic_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));
  
  // Get PHOS Geometry object
  AliPHOSGeometry *geom;
  if (!(geom = AliPHOSGeometry::GetInstance())) 
        geom = AliPHOSGeometry::GetInstance("IHEP","");
   
  //Define parameters
  Int_t nModules     = geom->GetNModules();
  fNCrystalsPhi = geom->GetNPhi()/fNTRUPhi ;// 64/4=16
  fNCrystalsZ   = geom->GetNZ()/fNTRUZ ;// 56/2=28

  //Intialize data members each time the trigger is called in event loop
  f2x2MaxAmp = -1; f2x2CrystalPhi = -1;  f2x2CrystalEta = -1;
  fnxnMaxAmp = -1; fnxnCrystalPhi = -1;  fnxnCrystalEta = -1;

  //Take the digits list if simulation
  if(fSimulation)
    fDigitsList = phosLoader->Digits() ;
  
  if(!fDigitsList)
    AliFatal("Digits not found !") ;
  
  //Fill TRU Matrix  
//   TClonesArray * amptrus   = new TClonesArray("TMatrixD",1000);
//   TClonesArray * ampmods   = new TClonesArray("TMatrixD",1000);
//   TClonesArray * timeRtrus = new TClonesArray("TMatrixD",1000);
  FillTRU(fDigitsList,geom) ;

  //Do Crystal Sliding and select Trigger
  //Initialize varible that will contain maximum amplitudes and 
  //its corresponding cell position in eta and phi, and time.
  TMatrixD   ampmax2(4,fNTRU) ;
  TMatrixD   ampmaxn(4,fNTRU) ;

  for(Int_t imod = 0 ; imod < nModules ; imod++) {

    //Do 2x2 and nxn sums, select maximums. 
    MakeSlidingCell(imod, ampmax2, ampmaxn);
    //Set the trigger
    SetTriggers(imod,ampmax2,ampmaxn) ;
  }

  fAmptrus->Delete();
//   delete amptrus; amptrus=0;
  fAmpmods->Delete();
//   delete ampmods; ampmods=0;
  fTimeRtrus->Delete();
//   delete timeRtrus; timeRtrus=0;
  //Print();

}
