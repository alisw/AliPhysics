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
//  Digits are grouped in TRU's  (Trigger Units). A TRU consists of 384 
//  cells ordered fNTRUPhi x fNTRUEta. The algorithm searches all possible 2x2 
//  and nxn (n is a multiple of 2) cell combinations per each TRU,  adding the 
//  digits amplitude and finding the maximum. If found, look if it is isolated.
//  Maxima are transformed in ADC time samples. Each time bin is compared to the trigger 
//  threshold until it is larger and then, triggers are set. Thresholds need to be fixed.  
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
//  ...
//  tr->Trigger(); //Execute Trigger
//  tr->Print(""); //Print results
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

// --- ALIROOT system ---
#include "AliRun.h" 
#include "AliRunLoader.h" 
#include "AliTriggerInput.h"
#include "AliEMCAL.h" 
#include "AliEMCALLoader.h" 
#include "AliEMCALDigit.h"
#include "AliEMCALTrigger.h" 
#include "AliEMCALGeometry.h"
#include "AliEMCALRawUtils.h"

ClassImp(AliEMCALTrigger)

//______________________________________________________________________
AliEMCALTrigger::AliEMCALTrigger()
  : AliTriggerDetector(), fGeom(0),
    f2x2MaxAmp(-1), f2x2CellPhi(-1),  f2x2CellEta(-1),
    f2x2SM(0),
    fnxnMaxAmp(-1), fnxnCellPhi(-1),  fnxnCellEta(-1),
    fnxnSM(0),
    fADCValuesHighnxn(0),fADCValuesLownxn(0),
    fADCValuesHigh2x2(0),fADCValuesLow2x2(0),
    fDigitsList(0),
    fL0Threshold(100),fL1JetLowPtThreshold(200),
    fL1JetMediumPtThreshold(500), fL1JetHighPtThreshold(1000),
    fPatchSize(1),  fIsolPatchSize(1), 
    f2x2AmpOutOfPatch(-1), fnxnAmpOutOfPatch(-1), 
    f2x2AmpOutOfPatchThres(100000),  fnxnAmpOutOfPatchThres(100000), 
    fIs2x2Isol(kFALSE), fIsnxnIsol(kFALSE),  
    fSimulation(kTRUE), fIsolateInSuperModule(kTRUE)
{
  //ctor 
  fADCValuesHighnxn = 0x0; //new Int_t[fTimeBins];
  fADCValuesLownxn  = 0x0; //new Int_t[fTimeBins];
  fADCValuesHigh2x2 = 0x0; //new Int_t[fTimeBins];
  fADCValuesLow2x2  = 0x0; //new Int_t[fTimeBins];

  SetName("EMCAL");
  CreateInputs();
   
   //Print("") ; 
}



//____________________________________________________________________________
AliEMCALTrigger::AliEMCALTrigger(const AliEMCALTrigger & trig) 
  : AliTriggerDetector(trig),
    fGeom(trig.fGeom),
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
    fPatchSize(trig.fPatchSize),
    fIsolPatchSize(trig.fIsolPatchSize), 
    f2x2AmpOutOfPatch(trig.f2x2AmpOutOfPatch), 
    fnxnAmpOutOfPatch(trig.fnxnAmpOutOfPatch), 
    f2x2AmpOutOfPatchThres(trig.f2x2AmpOutOfPatchThres),  
    fnxnAmpOutOfPatchThres(trig.fnxnAmpOutOfPatchThres), 
    fIs2x2Isol(trig.fIs2x2Isol),
    fIsnxnIsol(trig.fIsnxnIsol),  
    fSimulation(trig.fSimulation),
    fIsolateInSuperModule(trig.fIsolateInSuperModule)
{
  // cpy ctor
}

AliEMCALTrigger::~AliEMCALTrigger() {
  delete [] fADCValuesHighnxn; 
  delete [] fADCValuesLownxn;
  delete [] fADCValuesHigh2x2;
  delete [] fADCValuesLow2x2;
}

//----------------------------------------------------------------------
void AliEMCALTrigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "EMCAL_L0",       "EMCAL", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_JetHPt_L1","EMCAL", 1 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_JetMPt_L1","EMCAL", 1 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_JetLPt_L1","EMCAL", 1 ) );
 
}

//____________________________________________________________________________
Bool_t AliEMCALTrigger::IsPatchIsolated(Int_t iPatchType, const TClonesArray * ampmatrixes, const Int_t iSM, const Int_t mtru, const Float_t maxamp, const Int_t maxphi, const Int_t maxeta) {

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
  Int_t itru = mtru + iSM * fGeom->GetNTRU(); //number of tru, min 0 max 8*5.
  TMatrixD * ampmatrix   = 0x0;
  Int_t colborder = 0;
  Int_t rowborder = 0;
  
  if(fIsolateInSuperModule){
    ampmatrix = dynamic_cast<TMatrixD *>(ampmatrixes->At(iSM)) ;
    rowborder = fGeom->GetNPhi()*2; 
    colborder = fGeom->GetNZ()*2;
    AliDebug(2,"Isolate trigger in Module");
  }
  else{
    ampmatrix = dynamic_cast<TMatrixD *>(ampmatrixes->At(itru)) ;
    rowborder = fGeom->GetNCellsInTRUPhi();
    colborder = fGeom->GetNCellsInTRUEta();
    AliDebug(2,"Isolate trigger in TRU");
  }
  
  //Define patch cells
  Int_t isolcells   = fIsolPatchSize*(1+iPatchType);
  Int_t ipatchcells = 2*(1+fPatchSize*iPatchType);
  Int_t minrow      = maxphi - isolcells;
  Int_t mincol      = maxeta - isolcells;
  Int_t maxrow      = maxphi + isolcells + ipatchcells;
  Int_t maxcol      = maxeta + isolcells + ipatchcells;

  if (minrow < 0)
    minrow = 0;
  if (mincol < 0)
    mincol = 0;
  if (maxrow > rowborder)
    maxrow = rowborder;
  if (maxcol > colborder)
    maxcol = colborder;
  
  AliDebug(2,Form("Number of added Isol Cells %d, Patch Size %d",isolcells, ipatchcells));
  AliDebug(2,Form("Patch: minrow %d, maxrow %d, mincol %d, maxcol %d",minrow,maxrow,mincol,maxcol));
  
  //Add amplitudes in all isolation patch
  for(Int_t irow = minrow ; irow <  maxrow; irow ++)
    for(Int_t icol = mincol ; icol < maxcol ; icol ++)
      amp += (*ampmatrix)(irow,icol);
  
  AliDebug(2,Form("Type %d, Maximum amplitude %f, patch+isol square %f",iPatchType, maxamp, amp));

  if(amp < maxamp){
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
void AliEMCALTrigger::MakeSlidingCell(const TClonesArray * amptrus, const TClonesArray * timeRtrus, const Int_t isupermod,TMatrixD &ampmax2, TMatrixD &ampmaxn){
  
  //Sums energy of all possible 2x2 (L0) and nxn (L1) cells per each TRU. 
  //Fast signal in the experiment is given by 2x2 cells, 
  //for this reason we loop inside the TRU cells by 2. 

  //Declare and initialize variables
  Int_t nCellsPhi  = fGeom->GetNCellsInTRUPhi();
  if(isupermod > 9)
    nCellsPhi =  nCellsPhi / 2 ; //Half size SM. Not Final.
  // 12(tow)*2(cell)/1 TRU, cells in Phi in one TRU
  Int_t nCellsEta  = fGeom->GetNCellsInTRUEta();
  Int_t nTRU  = fGeom->GetNTRU();
  // 24(mod)*2(tower)/3 TRU, cells in Eta in one TRU
  //Int_t nTRU          = geom->GeNTRU();//3 TRU per super module

  Float_t amp2 = 0 ;
  Float_t ampn = 0 ; 
  for(Int_t i = 0; i < 4; i++){
    for(Int_t j = 0; j < nTRU; j++){
      ampmax2(i,j) = -1;
      ampmaxn(i,j) = -1;
    }
  }

  //Create matrix that will contain 2x2 amplitude sums
  //used to calculate the nxn sums
  TMatrixD tru2x2(nCellsPhi/2,nCellsEta/2) ;
  for(Int_t i = 0; i < nCellsPhi/2; i++)
    for(Int_t j = 0; j < nCellsEta/2; j++)
      tru2x2(i,j) = -1;
  
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

	//Fill matrix with added 2x2 cells for use in nxn sums
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

    //Sliding nxn, add nxn amplitudes  (OVERLAP)
    if(fPatchSize > 0){
      for(Int_t irow = 0 ; irow <  nCellsPhi/2; irow++){ 
	for(Int_t icol = 0 ; icol < nCellsEta/2 ; icol++){
	  ampn = 0;
	  if( (irow+fPatchSize) < nCellsPhi/2 && (icol+fPatchSize) < nCellsEta/2){//Avoid exit the TRU
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
	  if( (rown+i) < nCellsPhi && (coln+j) < nCellsEta){//Avoid exit the TRU
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
  printf( "               -2x2 Isolation Patch %d x %d, Amplitude out of 2x2 patch is %f, threshold %f, Isolated? %d \n", 
  	  2*fIsolPatchSize+2, 2*fIsolPatchSize+2,  f2x2AmpOutOfPatch,  f2x2AmpOutOfPatchThres,static_cast<Int_t> (fIs2x2Isol)) ; 
  if(fPatchSize > 0){
    printf( "             Patch Size, n x n: %d x %d cells\n",2*(fPatchSize+1), 2*(fPatchSize+1));
    printf( "               -nxn cells sum (overlapped)    : %10.2f, in Super Module %d\n",
	    fnxnMaxAmp,fnxnSM) ; 
    printf( "               -nxn from row %d to row %d and from column %d to column %d\n", fnxnCellPhi, fnxnCellPhi+4*fPatchSize, fnxnCellEta, fnxnCellEta+4*fPatchSize) ; 
    printf( "               -nxn Isolation Patch %d x %d, Amplitude out of nxn patch is %f, threshold %f, Isolated? %d \n", 
	    4*fIsolPatchSize+2*(fPatchSize+1),4*fIsolPatchSize+2*(fPatchSize+1) ,  fnxnAmpOutOfPatch,  fnxnAmpOutOfPatchThres,static_cast<Int_t> (fIsnxnIsol) ) ; 
  }
  
  printf( "             Isolate in SuperModule? %d\n",  
          fIsolateInSuperModule) ;  

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
void AliEMCALTrigger::SetTriggers(const TClonesArray * ampmatrix,const Int_t iSM, 
				  const TMatrixD &ampmax2, 
				  const TMatrixD &ampmaxn)  
{

  //Checks the 2x2 and nxn maximum amplitude per each TRU and 
  //compares with the different L0 and L1 triggers thresholds
  Float_t max2[] = {-1,-1,-1,-1} ;
  Float_t maxn[] = {-1,-1,-1,-1} ;
  Int_t   mtru2  = -1 ;
  Int_t   mtrun  = -1 ;

  Int_t nTRU = fGeom->GetNTRU();

  //Find maximum summed amplitude of all the TRU 
  //in a Super Module
    for(Int_t i = 0 ; i < nTRU ; i++){
      if(max2[0] < ampmax2(0,i) ){
	max2[0] =  ampmax2(0,i) ; // 2x2 summed max amplitude
	max2[1] =  ampmax2(1,i) ; // corresponding phi position in TRU
	max2[2] =  ampmax2(2,i) ; // corresponding eta position in TRU
	max2[3] =  ampmax2(3,i) ; // corresponding most recent time
	mtru2   = i ;
      }
      if(maxn[0] < ampmaxn(0,i) ){
	maxn[0] =  ampmaxn(0,i) ; // nxn summed max amplitude
	maxn[1] =  ampmaxn(1,i) ; // corresponding phi position in TRU
	maxn[2] =  ampmaxn(2,i) ; // corresponding eta position in TRU
	maxn[3] =  ampmaxn(3,i) ; // corresponding most recent time
	mtrun   = i ;
      }
    }

  //--------Set max amplitude if larger than in other Super Modules------------
  Float_t maxtimeR2 = -1 ;
  Float_t maxtimeRn = -1 ;
  static AliEMCALRawUtils rawUtil;
  Int_t nTimeBins = rawUtil.GetRawFormatTimeBins() ;

  //Set max of 2x2 amplitudes and select L0 trigger
  if(max2[0] > f2x2MaxAmp ){
    f2x2MaxAmp  = max2[0] ;
    f2x2SM      = iSM ;
    maxtimeR2   = max2[3] ;
    fGeom->GetCellPhiEtaIndexInSModuleFromTRUIndex(mtru2, 
						  static_cast<Int_t>(max2[1]),
						  static_cast<Int_t>(max2[2]),
						  f2x2CellPhi,f2x2CellEta) ;
    
    //Isolated patch?
    if(fIsolateInSuperModule)
      fIs2x2Isol =  IsPatchIsolated(0, ampmatrix, iSM, mtru2,  f2x2MaxAmp, f2x2CellPhi,f2x2CellEta) ;
    else
      fIs2x2Isol =  IsPatchIsolated(0, ampmatrix, iSM, mtru2,  f2x2MaxAmp,  static_cast<Int_t>(max2[1]), static_cast<Int_t>(max2[2])) ;

    //Transform digit amplitude in Raw Samples
    if (fADCValuesLow2x2 == 0) {
      fADCValuesLow2x2  = new Int_t[nTimeBins];
      fADCValuesHigh2x2 = new Int_t[nTimeBins];
    }
    rawUtil.RawSampledResponse(maxtimeR2, f2x2MaxAmp, fADCValuesHigh2x2, fADCValuesLow2x2) ; 
    
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
    fGeom->GetCellPhiEtaIndexInSModuleFromTRUIndex(mtrun, 
						  static_cast<Int_t>(maxn[1]),
						  static_cast<Int_t>(maxn[2]),
						  fnxnCellPhi,fnxnCellEta) ; 
    
    //Isolated patch?
    if(fIsolateInSuperModule)
      fIsnxnIsol =  IsPatchIsolated(1, ampmatrix, iSM, mtrun,  fnxnMaxAmp, fnxnCellPhi, fnxnCellEta) ;
    else
      fIsnxnIsol =  IsPatchIsolated(1, ampmatrix, iSM, mtrun,  fnxnMaxAmp,  static_cast<Int_t>(maxn[1]), static_cast<Int_t>(maxn[2])) ;
    
    //Transform digit amplitude in Raw Samples
    if (fADCValuesLownxn == 0) {
      fADCValuesHighnxn = new Int_t[nTimeBins];
      fADCValuesLownxn  = new Int_t[nTimeBins];
    }
    rawUtil.RawSampledResponse(maxtimeRn, fnxnMaxAmp, fADCValuesHighnxn, fADCValuesLownxn) ;
    
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
void AliEMCALTrigger::FillTRU(const TClonesArray * digits, TClonesArray * ampmatrix, TClonesArray * ampmatrixsmod, TClonesArray * timeRmatrix) {

//  Orders digits ampitudes list in fNTRU TRUs (384 cells) per supermodule. 
//  Each TRU is a TMatrixD, and they are kept in TClonesArrays. The number of 
//  TRU in phi is fNTRUPhi, and the number of TRU in eta is fNTRUEta.
//  Last 2 modules are half size in Phi, I considered that the number of TRU
//  is maintained for the last modules but decision not taken. If different, 
//  then this must be changed. Also fill a matrix with all amplitudes in supermodule for isolation studies. 
 
  //Initilize and declare variables
  //List of TRU matrices initialized to 0.
  Int_t nPhi = fGeom->GetNPhi();
  Int_t nZ = fGeom->GetNZ();
  Int_t nTRU = fGeom->GetNTRU();
  Int_t nTRUPhi = fGeom->GetNTRUPhi();
  Int_t nCellsPhi  = fGeom->GetNCellsInTRUPhi();
  Int_t nCellsPhi2 = fGeom->GetNCellsInTRUPhi();
  Int_t nCellsEta  = fGeom->GetNCellsInTRUEta();

  Int_t id       = -1; 
  Float_t amp    = -1;
  Float_t timeR  = -1;
  Int_t iSupMod  = -1;
  Int_t nModule  = -1;
  Int_t nIphi    = -1;
  Int_t nIeta    = -1;
  Int_t iphi     = -1;
  Int_t ieta     = -1;

  //List of TRU matrices initialized to 0.
  Int_t nSup = fGeom->GetNumberOfSuperModules();
  for(Int_t k = 0; k < nTRU*nSup; k++){
    TMatrixD amptrus(nCellsPhi,nCellsEta) ;
    TMatrixD timeRtrus(nCellsPhi,nCellsEta) ;
    // Do we need to initialise? I think TMatrixD does it by itself...
    for(Int_t i = 0; i < nCellsPhi; i++){
      for(Int_t j = 0; j < nCellsEta; j++){
	amptrus(i,j) = 0.0;
	timeRtrus(i,j) = 0.0;
      }
    }
    new((*ampmatrix)[k])   TMatrixD(amptrus) ;
    new((*timeRmatrix)[k]) TMatrixD(timeRtrus) ; 
  }
  
  //List of Modules matrices initialized to 0.
  for(Int_t k = 0; k < nSup ; k++){
    TMatrixD  ampsmods( nPhi*2, nZ*2) ;
    for(Int_t i = 0; i <  nPhi*2; i++){
      for(Int_t j = 0; j < nZ*2; j++){
	ampsmods(i,j)   = 0.0;
      }
    }
    new((*ampmatrixsmod)[k])   TMatrixD(ampsmods) ;
  }

  AliEMCALDigit * dig ;
  
  //Digits loop to fill TRU matrices with amplitudes.
  for(Int_t idig = 0 ; idig < digits->GetEntriesFast() ; idig++){
    
    dig = dynamic_cast<AliEMCALDigit *>(digits->At(idig)) ;
    amp    = dig->GetAmp() ;   // Energy of the digit (arbitrary units)
    id     = dig->GetId() ;    // Id label of the cell
    timeR  = dig->GetTimeR() ; // Earliest time of the digit
   
    //Get eta and phi cell position in supermodule
    Bool_t bCell = fGeom->GetCellIndex(id, iSupMod, nModule, nIphi, nIeta) ;
    if(!bCell)
      Error("FillTRU","Wrong cell id number") ;
    
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);

    //Check to which TRU in the supermodule belongs the cell. 
    //Supermodules are divided in a TRU matrix of dimension 
    //(fNTRUPhi,fNTRUEta).
    //Each TRU is a cell matrix of dimension (nCellsPhi,nCellsEta)

    //First calculate the row and column in the supermodule 
    //of the TRU to which the cell belongs.
    Int_t col   = ieta/nCellsEta; 
    Int_t row   = iphi/nCellsPhi; 
    if(iSupMod > 9)
      row   = iphi/nCellsPhi2; 
    //Calculate label number of the TRU
    Int_t itru  = row + col*nTRUPhi + iSupMod*nTRU ;  
 
    //Fill TRU matrix with cell values
    TMatrixD * amptrus   = dynamic_cast<TMatrixD *>(ampmatrix->At(itru)) ;
    TMatrixD * timeRtrus = dynamic_cast<TMatrixD *>(timeRmatrix->At(itru)) ;

    //Calculate row and column of the cell inside the TRU with number itru
    Int_t irow = iphi - row *  nCellsPhi;
    if(iSupMod > 9)
      irow = iphi - row *  nCellsPhi2;
    Int_t icol = ieta - col *  nCellsEta;
    
    (*amptrus)(irow,icol) = amp ;
    (*timeRtrus)(irow,icol) = timeR ;
    
    //####################SUPERMODULE MATRIX ##################
    TMatrixD * ampsmods   = dynamic_cast<TMatrixD *>(ampmatrixsmod->At(iSupMod)) ;
    (*ampsmods)(iphi,ieta)   = amp ;
    
  }
}
//____________________________________________________________________________
void AliEMCALTrigger::Trigger() 
{
  //Main Method to select triggers.
  AliRunLoader *runLoader = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (runLoader->GetDetectorLoader("EMCAL"));
 
  //Load EMCAL Geometry
  if (runLoader->GetAliRun() && runLoader->GetAliRun()->GetDetector("EMCAL"))
    fGeom = dynamic_cast<AliEMCAL*>(runLoader->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  if (fGeom == 0)
    fGeom = AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaulGeometryName());

  if (fGeom==0)
    AliFatal("Did not get geometry from EMCALLoader");
  
  //Define parameters
  Int_t nSuperModules = fGeom->GetNumberOfSuperModules() ; //12 SM in EMCAL
  Int_t nTRU       = fGeom->GetNTRU();    //3 TRU per super module

  //Intialize data members each time the trigger is called in event loop
  f2x2MaxAmp = -1; f2x2CellPhi = -1;  f2x2CellEta = -1;
  fnxnMaxAmp = -1; fnxnCellPhi = -1;  fnxnCellEta = -1;

  //Take the digits list if simulation
  if(fSimulation){
    runLoader->LoadDigits("EMCAL");
    fDigitsList = emcalLoader->Digits() ;
  }
  if(!fDigitsList)
    AliFatal("Digits not found !") ;
  
  //Take the digits list 
  
  //Fill TRU Matrix  
  TClonesArray * amptrus   = new TClonesArray("TMatrixD",1000);
  TClonesArray * ampsmods  = new TClonesArray("TMatrixD",1000);
  TClonesArray * timeRtrus = new TClonesArray("TMatrixD",1000);
  
  FillTRU(fDigitsList, amptrus, ampsmods, timeRtrus) ;

  //Do Cell Sliding and select Trigger
  //Initialize varible that will contain maximum amplitudes and 
  //its corresponding cell position in eta and phi, and time.
  TMatrixD ampmax2(4,nTRU) ;
  TMatrixD ampmaxn(4,nTRU) ;
  
  for(Int_t iSM = 0 ; iSM < nSuperModules ; iSM++) {
    //Do 2x2 and nxn sums, select maximums. 
    MakeSlidingCell(amptrus, timeRtrus, iSM, ampmax2, ampmaxn);
    
    //Set the trigger
    if(fIsolateInSuperModule)
      SetTriggers(ampsmods,iSM,ampmax2,ampmaxn) ;
    if(!fIsolateInSuperModule)
      SetTriggers(amptrus,iSM,ampmax2,ampmaxn) ;
  }
  
  amptrus->Delete();
  delete amptrus; amptrus = 0;
  ampsmods->Delete();
  delete ampsmods; ampsmods = 0;
  timeRtrus->Delete();
  delete timeRtrus; timeRtrus = 0;
  //Print();

}
