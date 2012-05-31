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
//  modules ordered fNTRUPhi x fNTRUEta. The algorithm searches all possible 2x2 
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
//  tr->SetL1GammaLowPtThreshold(1000);
//  tr->SetL1GammaMediumPtThreshold(10000);
//  tr->SetL1GammaHighPtThreshold(20000);
//  ...
//  tr->Trigger(); //Execute Trigger
//  tr->Print(""); //Print results
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////

#include <cassert>

// --- ROOT system ---
#include <TTree.h>
#include <TBranch.h>
#include <TBrowser.h>
#include <TH2F.h>

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
#include "AliLog.h"
#include "AliCaloConstants.h"
#include "AliEMCALRawResponse.h"

using namespace CALO;

ClassImp(AliEMCALTrigger)

TString AliEMCALTrigger::fgNameOfJetTriggers("EMCALJetTriggerL1");

//______________________________________________________________________
AliEMCALTrigger::AliEMCALTrigger()
  : AliTriggerDetector(), fGeom(0),
    f2x2MaxAmp(-1), f2x2ModulePhi(-1),  f2x2ModuleEta(-1),
    f2x2SM(0),
    fnxnMaxAmp(-1), fnxnModulePhi(-1),  fnxnModuleEta(-1),
    fnxnSM(0),
    fADCValuesHighnxn(0),fADCValuesLownxn(0),
    fADCValuesHigh2x2(0),fADCValuesLow2x2(0),
    fDigitsList(0),
    fL0Threshold(100),fL1GammaLowPtThreshold(200),
    fL1GammaMediumPtThreshold(500), fL1GammaHighPtThreshold(1000),
    fPatchSize(1),  fIsolPatchSize(1), 
    f2x2AmpOutOfPatch(-1), fnxnAmpOutOfPatch(-1), 
    f2x2AmpOutOfPatchThres(100000),  fnxnAmpOutOfPatchThres(100000), 
    fIs2x2Isol(kFALSE), fIsnxnIsol(kFALSE),  
    fSimulation(kTRUE), fIsolateInSuperModule(kTRUE), fTimeKey(kFALSE),
    fAmpTrus(0),fTimeRtrus(0),fAmpSMods(0),
    fTriggerPosition(6), fTriggerAmplitudes(4), 
    fNJetPatchPhi(3), fNJetPatchEta(3), fNJetThreshold(3), fL1JetThreshold(0), fJetMaxAmp(0),
    fAmpJetMatrix(0), fJetMatrixE(0), fAmpJetMax(6,1), fVZER0Mult(0.)
{
  //ctor 
  fADCValuesHighnxn = 0x0; //new Int_t[fTimeBins];
  fADCValuesLownxn  = 0x0; //new Int_t[fTimeBins];
  fADCValuesHigh2x2 = 0x0; //new Int_t[fTimeBins];
  fADCValuesLow2x2  = 0x0; //new Int_t[fTimeBins];

  SetName("EMCAL");
  // Define jet threshold - can not change from outside now
  // fNJetThreshold  = 7; // For MB Pythia suppression
  fNJetThreshold  = 10;   // Hijing  
  fL1JetThreshold = new Double_t[fNJetThreshold];
  if(fNJetThreshold == 7) {
    fL1JetThreshold[0] =  5./0.0153;
    fL1JetThreshold[1] =  8./0.0153;
    fL1JetThreshold[2] = 10./0.0153;
    fL1JetThreshold[3] = 12./0.0153;
    fL1JetThreshold[4] = 13./0.0153;
    fL1JetThreshold[5] = 14./0.0153;
    fL1JetThreshold[6] = 15./0.0153;
  } else if(fNJetThreshold == 10) {
    Double_t thGev[10]={5.,8.,10., 12., 13.,14.,15., 17., 20., 25.};
    for(Int_t i=0; i<fNJetThreshold; i++) fL1JetThreshold[i] =  thGev[i]/0.0153;
  } else {
    fL1JetThreshold[0] =  5./0.0153;
    fL1JetThreshold[1] = 10./0.0153;
    fL1JetThreshold[2] = 15./0.0153;
    fL1JetThreshold[3] = 20./0.0153;
    fL1JetThreshold[4] = 25./0.0153;
  }
  //
  CreateInputs();

  fInputs.SetName("TriggersInputs");   
   //Print("") ; 
}

//____________________________________________________________________________
AliEMCALTrigger::~AliEMCALTrigger() {
	
  //dtor
	
  if(GetTimeKey()) {
    delete [] fADCValuesHighnxn; 
    delete [] fADCValuesLownxn;
    delete [] fADCValuesHigh2x2;
    delete [] fADCValuesLow2x2;
  }
  if(fAmpTrus)      {fAmpTrus->Delete();   delete fAmpTrus;}
  if(fTimeRtrus)    {fTimeRtrus->Delete(); delete fTimeRtrus;}
  if(fAmpSMods)     {fAmpSMods->Delete();  delete fAmpSMods;}
  if(fAmpJetMatrix) delete fAmpJetMatrix;
  if(fJetMatrixE)   delete fJetMatrixE;
  if(fL1JetThreshold) delete [] fL1JetThreshold;
}

//----------------------------------------------------------------------
void AliEMCALTrigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;

   // Second parameter should be detector name = "EMCAL"
   TString det("EMCAL"); // Apr 29, 2008
   fInputs.AddLast( new AliTriggerInput( det+"_L0",          det, 0x02) );
   fInputs.AddLast( new AliTriggerInput( det+"_GammaHPt_L1", det, 0x04 ) );
   fInputs.AddLast( new AliTriggerInput( det+"_GammaMPt_L1", det, 0x08 ) );
   fInputs.AddLast( new AliTriggerInput( det+"_GammaLPt_L1", det, 0x016 ) );
   fInputs.AddLast( new AliTriggerInput( det+"_JetHPt_L1", det, 0x032 ) );
   fInputs.AddLast( new AliTriggerInput( det+"_JetMPt_L1", det, 0x048 ) );
   fInputs.AddLast( new AliTriggerInput( det+"_JetLPt_L1", det, 0x064 ) );

   if(fNJetThreshold<=0) return;
   // Jet Trigger(s)
   UInt_t level = 0x032;
   for(Int_t i=0; i<fNJetThreshold; i++ ) {
     TString name(GetNameOfJetTrigger(i));
     TString title("EMCAL Jet triger L1 :"); // unused now
     // 0.0153 - hard coded now
     title += Form("Th %i(%5.1f GeV) :", (Int_t)fL1JetThreshold[i], fL1JetThreshold[i] * 0.0153); 
     title += Form("patch %ix%i~(%3.2f(#phi)x%3.2f(#eta)) ", 
		 fNJetPatchPhi, fNJetPatchEta, 0.11*(fNJetPatchPhi), 0.11*(fNJetPatchEta)); 
     fInputs.AddLast( new AliTriggerInput(name, det, UChar_t(level)) );
     level *= 2;
   }
 
}

//____________________________________________________________________________
Bool_t AliEMCALTrigger::IsPatchIsolated(Int_t iPatchType, const TClonesArray * ampmatrixes, const Int_t iSM, const Int_t mtru, const Float_t maxamp, const Int_t maxphi, const Int_t maxeta) {

  // Nov 8, 2007 
  // EMCAL RTU size is 4modules(phi) x 24modules (eta)
  // So maximum size of patch is 4modules x 4modules (EMCAL L0 trigger). 
  // Calculate if the maximum patch found is isolated, find amplitude around maximum (2x2 or nxn) patch, 
  // inside isolation patch . iPatchType = 0 means calculation for 2x2 patch, 
  // iPatchType = 1 means calculation for nxn patch.
  // In the next table there is an example of the different options of patch size and isolation patch size:
  //                                                           Patch Size (fPatchSize)
  //                                           0                          1              
  //          fIsolPatchSize    0             2x2 (not overlap)   4x4 (overlapped)       
  //                            1             4x4                      8x8               
                          
  Bool_t b = kFALSE;
  if(!ampmatrixes) return kFALSE;
  
  // Get matrix of TRU or Module with maximum amplitude patch.
  Int_t itru = mtru + iSM * fGeom->GetNTRU(); //number of tru, min 0 max 3*12=36.
  TMatrixD * ampmatrix   = 0x0;
  Int_t colborder = 0;
  Int_t rowborder = 0;
  static int keyPrint = 0;
  if(keyPrint) AliDebug(2,Form(" IsPatchIsolated : iSM %i mtru %i itru %i maxphi %i maxeta %i \n", iSM, mtru, itru, maxphi, maxeta));
  
  if(fIsolateInSuperModule){ // ?
    ampmatrix = dynamic_cast<TMatrixD *>(ampmatrixes->At(iSM)) ;
    rowborder = fGeom->GetNPhi();
    colborder = fGeom->GetNZ();
    AliDebug(2,"Isolate trigger in Module");
  } else{
    ampmatrix = dynamic_cast<TMatrixD *>(ampmatrixes->At(itru)) ;
    rowborder = fGeom->GetNModulesInTRUPhi();
    colborder = fGeom->GetNModulesInTRUEta();
    AliDebug(2,"Isolate trigger in TRU");
  }
  if(iSM>9) rowborder /= 2; // half size in phi
  
  if(!ampmatrixes || !ampmatrix){
    AliError("Could not recover the matrix with the amplitudes");
    return kFALSE;
  }
  
  //Define patch modules - what is this ??
  Int_t isolmodules   = fIsolPatchSize*(1+iPatchType);
  Int_t ipatchmodules = 2*(1+fPatchSize*iPatchType);
  Int_t minrow      = maxphi - isolmodules;
  Int_t mincol      = maxeta - isolmodules;
  Int_t maxrow      = maxphi + isolmodules + ipatchmodules;
  Int_t maxcol      = maxeta + isolmodules + ipatchmodules;

  minrow =  minrow<0?0 :minrow;
  mincol =  mincol<0?0 :mincol;

  maxrow =  maxrow>rowborder?rowborder :maxrow;
  maxcol =  maxcol>colborder?colborder :maxcol;
  
  //printf("%s\n",Form("Number of added Isol Modules %d, Patch Size %d",isolmodules, ipatchmodules));
  //printf("%s\n",Form("Patch: minrow %d, maxrow %d, mincol %d, maxcol %d",minrow,maxrow,mincol,maxcol));
  //  AliDebug(2,Form("Number of added Isol Modules %d, Patch Size %d",isolmodules, ipatchmodules));
  //AliDebug(2,Form("Patch: minrow %d, maxrow %d, mincol %d, maxcol %d",minrow,maxrow,mincol,maxcol));
  
  //Add amplitudes in all isolation patch
  Float_t amp = 0.;
  for(Int_t irow = minrow ; irow <  maxrow; irow ++)
    for(Int_t icol = mincol ; icol < maxcol ; icol ++)
      amp += (*ampmatrix)(irow,icol);
  
  AliDebug(2,Form("Type %d, Maximum amplitude %f, patch+isol square %f",iPatchType, maxamp, amp));

  if(amp < maxamp){
    //    AliError(Form("Bad sum: Type %d, Maximum amplitude %f, patch+isol square %f",iPatchType, maxamp, amp));
    //    ampmatrix->Print();
    return kFALSE;
  } else {
    amp-=maxamp; //Calculate energy in isolation patch that do not comes from maximum patch.
  }
  
  AliDebug(2, Form("Maximum amplitude %f, Out of patch %f",maxamp, amp));

  //Fill isolation amplitude data member and say if patch is isolated.
  if(iPatchType == 0){ //2x2 case
    f2x2AmpOutOfPatch = amp;   
    if(amp < f2x2AmpOutOfPatchThres) b=kTRUE;
  } else  if(iPatchType == 1){ //nxn case
    fnxnAmpOutOfPatch = amp;   
    if(amp < fnxnAmpOutOfPatchThres) b=kTRUE;
  }

  if(keyPrint) AliDebug(2,Form(" IsPatchIsolated - OUT \n"));

  return b;

}

/*
//____________________________________________________________________________
void AliEMCALTrigger::MakeSlidingCell(const TClonesArray * amptrus, const TClonesArray * timeRtrus, const Int_t isupermod,TMatrixD &ampmax2, TMatrixD &ampmaxn){
  
  //Sums energy of all possible 2x2 (L0) and nxn (L1) modules per each TRU. 
  //Fast signal in the experiment is given by 2x2 modules, 
  //for this reason we loop inside the TRU modules by 2. 

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
*/
//____________________________________________________________________________
void AliEMCALTrigger::MakeSlidingTowers(const TClonesArray * amptrus, const TClonesArray * timeRtrus, 
                                        const Int_t isupermod,TMatrixD &ampmax2, TMatrixD &ampmaxn){
  
  // Output from module (2x2 cells from one module)
  Int_t nModulesPhi  = fGeom->GetNModulesInTRUPhi(); // now 4 modules (3 div in phi)
  if(isupermod > 9)
    nModulesPhi =  nModulesPhi / 2 ; // Half size SM. Not Final.
  // 
  Int_t nModulesEta  = fGeom->GetNModulesInTRUEta(); // now 24 modules (no division in eta)
  Int_t nTRU         = fGeom->GetNTRU();
  static int keyPrint = 0;
  if(keyPrint) AliDebug(2,Form("MakeSlidingTowers : nTRU %i nModulesPhi %i nModulesEta %i ", 
                               nTRU, nModulesPhi, nModulesEta ));
  
  Float_t amp2 = 0 ;
  Float_t ampn = 0 ; 
  for(Int_t i = 0; i < 4; i++){
    for(Int_t j = 0; j < nTRU; j++){
      ampmax2(i,j) = ampmaxn(i,j) = -1;
    }
  }
  
  // Create matrix that will contain 2x2 amplitude sums
  // used to calculate the nxn sums
  TMatrixD tru2x2(nModulesPhi/2,nModulesEta/2);
  
  // Loop over all TRUS in a supermodule
  for(Int_t itru = 0 + isupermod * nTRU ; itru < (isupermod+1)*nTRU ; itru++) {
    TMatrixD * amptru   = dynamic_cast<TMatrixD *>(amptrus->At(itru)) ;
    TMatrixD * timeRtru = dynamic_cast<TMatrixD *>(timeRtrus->At(itru)) ;
    Int_t mtru = itru - isupermod*nTRU ; // Number of TRU in Supermodule !!
    
    if(!amptru || !timeRtru){
      AliError("Amplitude or Time TRU matrix not available");
      return;
    }
    
    // Sliding 2x2, add 2x2 amplitudes (NOT OVERLAP)
    for(Int_t irow = 0 ; irow <  nModulesPhi; irow +=2){ 
      for(Int_t icol = 0 ; icol < nModulesEta ; icol +=2){
        amp2 = (*amptru)(irow,icol) +(*amptru)(irow+1,icol)+
        (*amptru)(irow,icol+1)+(*amptru)(irow+1,icol+1);
        
        //Fill matrix with added 2x2 towers for use in nxn sums
        tru2x2(irow/2,icol/2) = amp2 ;
        //Select 2x2 maximum sums to select L0 
        if(amp2 > ampmax2(0,mtru)){
          ampmax2(0,mtru) = amp2 ; 
          ampmax2(1,mtru) = irow;
          ampmax2(2,mtru) = icol;
        }
      }
    }
    
    ampmax2(3,mtru) = 0.;
    if(GetTimeKey()) {
      // Find most recent time in the selected 2x2 towers
      Int_t row2 =  static_cast <Int_t> (ampmax2(1,mtru));
      Int_t col2 =  static_cast <Int_t> (ampmax2(2,mtru));
      for(Int_t i = 0; i<2; i++){
        for(Int_t j = 0; j<2; j++){
          if((*amptru)(row2+i,col2+j) > 0 &&  (*timeRtru)(row2+i,col2+j)> 0){	  
            if((*timeRtru)(row2+i,col2+j) >  ampmax2(3,mtru)  )
              ampmax2(3,mtru) =  (*timeRtru)(row2+i,col2+j); // max time
          }
        }
      }
    }
    
    //Sliding nxn, add nxn amplitudes  (OVERLAP)
    if(fPatchSize > 0){
      for(Int_t irow = 0 ; irow <  nModulesPhi/2; irow++){ 
        for(Int_t icol = 0 ; icol < nModulesEta/2; icol++){
          ampn = 0;
          if( (irow+fPatchSize) < nModulesPhi/2 && (icol+fPatchSize) < nModulesEta/2){ //Avoid exit the TRU
            for(Int_t i = 0 ; i <= fPatchSize ; i++)
              for(Int_t j = 0 ; j <= fPatchSize ; j++)
                ampn += tru2x2(irow+i,icol+j);
            //Select nxn maximum sums to select L1 
            if(ampn > ampmaxn(0,mtru)){
              ampmaxn(0,mtru) = ampn ; 
              ampmaxn(1,mtru) = irow;
              ampmaxn(2,mtru) = icol;
            }
          }
        }
      }
      
      ampmaxn(3,mtru) = 0.; // Was 1 , I don't know why
      if(GetTimeKey()) {
        //Find most recent time in selected nxn cell
        Int_t rown =  static_cast <Int_t> (ampmaxn(1,mtru));
        Int_t coln =  static_cast <Int_t> (ampmaxn(2,mtru));
        for(Int_t i = 0; i<4*fPatchSize; i++){
          for(Int_t j = 0; j<4*fPatchSize; j++){
            if( (rown+i) < nModulesPhi && (coln+j) < nModulesEta){//Avoid exit the TRU
              if((*amptru)(rown+i,coln+j) > 0 &&  (*timeRtru)(rown+i,coln+j)> 0){
                if((*timeRtru)(rown+i,coln+j) >  ampmaxn(3,mtru)  )
                  ampmaxn(3,mtru) =  (*timeRtru)(rown+i,coln+j); // max time
              }
            }
          }
        }
      }
    } else { // copy 2x2 to nxn  
      ampmaxn(0,mtru) =  ampmax2(0,mtru); 
      ampmaxn(1,mtru) =  ampmax2(1,mtru);
      ampmaxn(2,mtru) =  ampmax2(2,mtru);
      ampmaxn(3,mtru) =  ampmax2(3,mtru);
    }
  }
  if(keyPrint) AliDebug(2,Form(" : MakeSlidingTowers -OUt \n"));
}

//____________________________________________________________________________
void AliEMCALTrigger::Print(const Option_t * opt) const 
{
  
  //Prints main parameters
  
  if(! opt)
    return;
  AliTriggerInput* in = 0x0 ;
  AliInfo(Form(" fSimulation %i (input option) : #digits %i\n", fSimulation, fDigitsList->GetEntries()));
  AliInfo(Form(" fTimeKey    %i  \n ", fTimeKey));

  AliInfo(Form("\t Maximum Amplitude after Sliding Cell, \n")) ; 
  AliInfo(Form("\t -2x2 cells sum (not overlapped): %10.2f, in Super Module %d\n",
	  f2x2MaxAmp,f2x2SM)) ; 
  AliInfo(Form("\t -2x2 from row %d to row %d and from column %d to column %d\n", f2x2ModulePhi, f2x2ModulePhi+2, f2x2ModuleEta, f2x2ModuleEta+2)); 
  AliInfo(Form("\t -2x2 Isolation Patch %d x %d, Amplitude out of 2x2 patch is %f, threshold %f, Isolated? %d \n", 2*fIsolPatchSize+2, 2*fIsolPatchSize+2,  f2x2AmpOutOfPatch,  f2x2AmpOutOfPatchThres,static_cast<Int_t> (fIs2x2Isol))); 
  if(fPatchSize > 0){
    AliInfo(Form("\t Patch Size, n x n: %d x %d cells\n",2*(fPatchSize+1), 2*(fPatchSize+1)));
    AliInfo(Form("\t -nxn cells sum (overlapped)    : %10.2f, in Super Module %d\n", fnxnMaxAmp,fnxnSM)); 
    AliInfo(Form("\t -nxn from row %d to row %d and from column %d to column %d\n", fnxnModulePhi, fnxnModulePhi+4*fPatchSize, fnxnModuleEta, fnxnModuleEta+4*fPatchSize)) ; 
    AliInfo(Form("\t -nxn Isolation Patch %d x %d, Amplitude out of nxn patch is %f, threshold %f, Isolated? %d \n", 4*fIsolPatchSize+2*(fPatchSize+1),4*fIsolPatchSize+2*(fPatchSize+1) ,  fnxnAmpOutOfPatch,  fnxnAmpOutOfPatchThres,static_cast<Int_t> (fIsnxnIsol) )); 
  }
  
  AliInfo(Form("\t Isolate in SuperModule? %d\n", fIsolateInSuperModule)) ;  
  AliInfo(Form("\t Threshold for LO %10.2f\n", fL0Threshold));  

  in = (AliTriggerInput*)fInputs.FindObject( "EMCAL_L0" );
  if(in->GetValue())
    AliInfo(Form("\t *** EMCAL LO is set ***\n")); 
  
  AliInfo(Form("\t Gamma Low Pt Threshold for L1 %10.2f\n", fL1GammaLowPtThreshold));
  in = (AliTriggerInput*)fInputs.FindObject( "EMCAL_GammaLPt_L1" );
  if(in->GetValue())
    AliInfo(Form("\t *** EMCAL Gamma Low Pt for L1 is set ***\n"));

  AliInfo(Form("\t Gamma Medium Pt Threshold for L1 %10.2f\n", fL1GammaMediumPtThreshold));  
  in = (AliTriggerInput*) fInputs.FindObject( "EMCAL_GammaMPt_L1" );
  if(in->GetValue())
    AliInfo(Form("\t *** EMCAL Gamma Medium Pt for L1 is set ***\n"));

  AliInfo(Form("\t Gamma High Pt Threshold for L1 %10.2f\n", fL1GammaHighPtThreshold));  
  in = (AliTriggerInput*) fInputs.FindObject( "EMCAL_GammaHPt_L1" );
  if(in->GetValue())
    AliInfo(Form("\t *** EMCAL Gamma High Pt for L1 is set ***\n")) ;

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
  // Int_t nTimeBins = rawUtil.GetRawFormatTimeBins() ;
  Int_t nTimeBins = TIMEBINS;  //changed by PTH

  //Set max of 2x2 amplitudes and select L0 trigger
  if(max2[0] > f2x2MaxAmp ){
    //    if(max2[0] > 5) printf(" L0 : iSM %i: max2[0] %5.0f :  max2[3] %5.0f  (maxtimeR2) \n", 
    //			   iSM, max2[0], max2[3]);
    f2x2MaxAmp  = max2[0] ;
    f2x2SM      = iSM ;
    maxtimeR2   = max2[3] ;
    fGeom->GetModulePhiEtaIndexInSModuleFromTRUIndex(mtru2, 
						  static_cast<Int_t>(max2[1]),
						  static_cast<Int_t>(max2[2]),
						   f2x2ModulePhi,f2x2ModuleEta);
    
    //Isolated patch?
    if(fIsolateInSuperModule)
      fIs2x2Isol =  IsPatchIsolated(0, ampmatrix, iSM, mtru2,  f2x2MaxAmp, f2x2ModulePhi,f2x2ModuleEta) ;
    else
      fIs2x2Isol =  IsPatchIsolated(0, ampmatrix, iSM, mtru2,  f2x2MaxAmp,  static_cast<Int_t>(max2[1]), static_cast<Int_t>(max2[2])) ;

    if(GetTimeKey()) {
    //Transform digit amplitude in Raw Samples
      if (fADCValuesLow2x2 == 0) {
        fADCValuesLow2x2  = new Int_t[nTimeBins];
        fADCValuesHigh2x2 = new Int_t[nTimeBins];
      }
      //printf(" maxtimeR2 %12.5e (1)\n", maxtimeR2);
      //  rawUtil.RawSampledResponse(maxtimeR2 * AliEMCALRawUtils::GetRawFormatTimeBin(), 
      //				 f2x2MaxAmp, fADCValuesHigh2x2, fADCValuesLow2x2) ; 
       
      //  rawUtil.RawSampledResponse(maxtimeR2*TIMEBINMAX/TIMEBINS, 
      //				 f2x2MaxAmp, fADCValuesHigh2x2, fADCValuesLow2x2) ; 
      
      AliEMCALRawResponse::RawSampledResponse( maxtimeR2*TIMEBINMAX/TIMEBINS, 
					       f2x2MaxAmp, fADCValuesHigh2x2, fADCValuesLow2x2) ; 

    // Set Trigger Inputs, compare ADC time bins until threshold is attained
    // Set L0
      for(Int_t i = 0 ; i < nTimeBins ; i++){
      //      printf(" fADCValuesHigh2x2[%i] %i : %i \n", i, fADCValuesHigh2x2[i], fADCValuesLow2x2[i]); 
        if(fADCValuesHigh2x2[i] >= fL0Threshold          || fADCValuesLow2x2[i] >= fL0Threshold){
	  SetInput("EMCAL_L0") ;
	  break;
        }
      }
    } else {
      // Nov 5 - no analysis of time information
      if(f2x2MaxAmp >= fL0Threshold) { // should add the low amp too
        SetInput("EMCAL_L0");
      }
    }
  }
  
  //------------Set max of nxn amplitudes and select L1 trigger---------
  if(maxn[0] > fnxnMaxAmp ){
    fnxnMaxAmp  = maxn[0] ;
    fnxnSM      = iSM ;
    maxtimeRn   = maxn[3] ;
    fGeom->GetModulePhiEtaIndexInSModuleFromTRUIndex(mtrun, 
						  static_cast<Int_t>(maxn[1]),
						  static_cast<Int_t>(maxn[2]),
						  fnxnModulePhi,fnxnModuleEta) ; 
    
    //Isolated patch?
    if(fIsolateInSuperModule)
      fIsnxnIsol =  IsPatchIsolated(1, ampmatrix, iSM, mtrun,  fnxnMaxAmp, fnxnModulePhi, fnxnModuleEta) ;
    else
      fIsnxnIsol =  IsPatchIsolated(1, ampmatrix, iSM, mtrun,  fnxnMaxAmp,  static_cast<Int_t>(maxn[1]), static_cast<Int_t>(maxn[2])) ;
    
    if(GetTimeKey()) {
    //Transform digit amplitude in Raw Samples
      if (fADCValuesLownxn == 0) {
        fADCValuesHighnxn = new Int_t[nTimeBins];
        fADCValuesLownxn  = new Int_t[nTimeBins];
      }
      //  rawUtil.RawSampledResponse(maxtimeRn * AliEMCALRawUtils::GetRawFormatTimeBin(), 
      //   fnxnMaxAmp, fADCValuesHighnxn, fADCValuesLownxn) ;

      //rawUtil.RawSampledResponse(maxtimeRn*TIMEBINMAX/TIMEBINS, 
      //				 fnxnMaxAmp, fADCValuesHighnxn, fADCValuesLownxn) ;

      AliEMCALRawResponse::RawSampledResponse (maxtimeRn*TIMEBINMAX/TIMEBINS, 
					       fnxnMaxAmp, fADCValuesHighnxn, fADCValuesLownxn) ;
      
    //Set Trigger Inputs, compare ADC time bins until threshold is attained
    //SetL1 Low
      for(Int_t i = 0 ; i < nTimeBins ; i++){
        if(fADCValuesHighnxn[i] >= fL1GammaLowPtThreshold  || fADCValuesLownxn[i] >= fL1GammaLowPtThreshold){
	  SetInput("EMCAL_GammaLPt_L1") ;
	  break; 
        }
      }
    
    //SetL1 Medium
      for(Int_t i = 0 ; i < nTimeBins ; i++){
        if(fADCValuesHighnxn[i] >= fL1GammaMediumPtThreshold || fADCValuesLownxn[i] >= fL1GammaMediumPtThreshold){
	  SetInput("EMCAL_GammaMPt_L1") ;
	  break;
        }
      }
    
    //SetL1 High
      for(Int_t i = 0 ; i < nTimeBins ; i++){
        if(fADCValuesHighnxn[i] >= fL1GammaHighPtThreshold || fADCValuesLownxn[i] >= fL1GammaHighPtThreshold){
	  SetInput("EMCAL_GammaHPt_L1") ;
	  break;
        }
      }
    } else {
      // Nov 5 - no analysis of time information
      if(fnxnMaxAmp >= fL1GammaLowPtThreshold) { // should add the low amp too
	SetInput("EMCAL_GammaLPt_L1") ; //SetL1 Low
      }
      if(fnxnMaxAmp >= fL1GammaMediumPtThreshold) { // should add the low amp too
	SetInput("EMCAL_GammaMPt_L1") ; //SetL1 Medium
      }
      if(fnxnMaxAmp >= fL1GammaHighPtThreshold) { // should add the low amp too
	SetInput("EMCAL_GammaHPt_L1") ; //SetL1 High
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
 
//  Initilize and declare variables
//  List of TRU matrices initialized to 0.
//  printf("<I> AliEMCALTrigger::FillTRU() started : # digits %i\n", digits->GetEntriesFast());

// Nov 2, 2007.
// One input per EMCAL module so size of matrix is reduced by 4 (2x2 division case) 

  Int_t nPhi        = fGeom->GetNPhi();
  Int_t nZ          = fGeom->GetNZ();
  Int_t nTRU        = fGeom->GetNTRU();
  //  Int_t nTRUPhi     = fGeom->GetNTRUPhi();
  Int_t nModulesPhi  = fGeom->GetNModulesInTRUPhi();
  Int_t nModulesPhi2 = fGeom->GetNModulesInTRUPhi();
  Int_t nModulesEta  = fGeom->GetNModulesInTRUEta();
  //  printf("<I> AliEMCALTrigger::FillTRU() nTRU %i  nTRUPhi %i : nModulesPhi %i nModulesEta %i \n", 
  //	 nTRU, nTRUPhi, nModulesPhi, nModulesEta);

  Int_t id       = -1; 
  Float_t amp    = -1;
  Float_t timeR  = -1;
  Int_t iSupMod  = -1;
  Int_t nModule  = -1;
  Int_t nIphi    = -1;
  Int_t nIeta    = -1;
  Int_t iphi     = -1;
  Int_t ieta     = -1;
  // iphim, ietam - module indexes in SM 
  Int_t iphim    = -1;
  Int_t ietam    = -1;

  //List of TRU matrices initialized to 0.
  Int_t nSup = fGeom->GetNumberOfSuperModules();
  for(Int_t k = 0; k < nTRU*nSup; k++){
    TMatrixD amptrus(nModulesPhi,nModulesEta) ;
    TMatrixD timeRtrus(nModulesPhi,nModulesEta) ;
    // Do we need to initialise? I think TMatrixD does it by itself...
    for(Int_t i = 0; i < nModulesPhi; i++){
      for(Int_t j = 0; j < nModulesEta; j++){
	amptrus(i,j) = 0.0;
	timeRtrus(i,j) = 0.0;
      }
    }
    new((*ampmatrix)[k])   TMatrixD(amptrus) ;
    new((*timeRmatrix)[k]) TMatrixD(timeRtrus) ; 
  }
  
  // List of Modules matrices initialized to 0.
  for(Int_t k = 0; k < nSup ; k++){
    int mphi = nPhi;
    //    if(nSup>9) mphi = nPhi/2; // the same size
    TMatrixD  ampsmods( mphi, nZ);
    for(Int_t i = 0; i <  mphi; i++){
      for(Int_t j = 0; j < nZ; j++){
	ampsmods(i,j)   = 0.0;
      }
    }
    new((*ampmatrixsmod)[k]) TMatrixD(ampsmods) ;
  }

  AliEMCALDigit * dig ;
  
  //Digits loop to fill TRU matrices with amplitudes.
  for(Int_t idig = 0 ; idig < digits->GetEntriesFast() ; idig++){
    
    dig = dynamic_cast<AliEMCALDigit *>(digits->At(idig)) ;
    if(dig){
      amp    = Float_t(dig->GetAmplitude()); // Energy of the digit (arbitrary units)
      id     = dig->GetId() ;          // Id label of the cell
      timeR  = dig->GetTimeR() ;       // Earliest time of the digit
      if(amp<=0.0) AliDebug(1,Form(" id %i amp %f \n", id, amp));
      // printf(" FILLTRU : timeR %10.5e time %10.5e : amp %10.5e \n", timeR, dig->GetTime(), amp);
      // Get eta and phi cell position in supermodule
      Bool_t bCell = fGeom->GetCellIndex(id, iSupMod, nModule, nIphi, nIeta) ;
      if(!bCell)
        AliError(Form("%i Wrong cell id number %i ", idig, id)) ;
      
      fGeom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);
      // iphim, ietam - module indexes in SM
      fGeom->GetModuleIndexesFromCellIndexesInSModule(iSupMod,iphi,ieta, iphim, ietam, nModule); 
      //if(iSupMod >9) 
      //printf("iSupMod %i nModule %i iphi %i  ieta %i  iphim %i  ietam %i \n",
      //iSupMod,nModule, iphi, ieta, iphim, ietam); 
      
      // Check to which TRU in the supermodule belongs the cell. 
      // Supermodules are divided in a TRU matrix of dimension 
      // (fNTRUPhi,fNTRUEta).
      // Each TRU is a cell matrix of dimension (nModulesPhi,nModulesEta)
      
      // First calculate the row and column in the supermodule 
      // of the TRU to which the cell belongs.
      Int_t row   = iphim / nModulesPhi;
      Int_t col   = ietam / nModulesEta;
      //Calculate label number of the TRU
      Int_t itru  = fGeom->GetAbsTRUNumberFromNumberInSm(row, col, iSupMod);
      
      //Fill TRU matrix with cell values
      TMatrixD * amptrus   = dynamic_cast<TMatrixD *>(ampmatrix->At(itru)) ;
      TMatrixD * timeRtrus = dynamic_cast<TMatrixD *>(timeRmatrix->At(itru)) ;
      
      if(!amptrus || !timeRtrus){
        AliError("Could not recover the TRU matrix with amplitudes or times");
      }
      else{
        //Calculate row and column of the module inside the TRU with number itru
        Int_t irow = iphim - row * nModulesPhi;
        if(iSupMod > 9)
          irow = iphim - row *  nModulesPhi2; // size of matrix the same
        Int_t icol = ietam - col * nModulesEta;
      
        (*amptrus)(irow,icol)  += amp ;
        if((*timeRtrus)(irow,icol) <0.0 || (*timeRtrus)(irow,icol) <= timeR){ // ??
          (*timeRtrus)(irow,icol) = timeR ;
        }
      }
      //printf(" ieta %i iphi %i iSM %i || col %i row %i : itru %i -> amp %f\n", 
      //	   ieta, iphi, iSupMod, col, row, itru, amp);     
      //####################SUPERMODULE MATRIX ##################
      TMatrixD * ampsmods   = dynamic_cast<TMatrixD *>(ampmatrixsmod->At(iSupMod)) ;
      if(!ampsmods){
        AliError("Could not recover the matrix per SM");
        continue;
      }
      (*ampsmods)(iphim,ietam)  += amp ;
      //    printf(" id %i iphim %i ietam %i SM %i : irow %i icol %i itru %i : amp %6.0f\n", 
      //id, iphim, ietam, iSupMod, irow, icol, itru, amp); 
    }
    else AliError("Could not recover the digit");
  }
  //assert(0);
  //printf("<I> AliEMCALTrigger::FillTRU() is ended \n");
}

//____________________________________________________________________________
void AliEMCALTrigger::Trigger() 
{
  //Main Method to select triggers.
  TH1::AddDirectory(0);

  AliRunLoader *runLoader = AliRunLoader::Instance();
  AliEMCALLoader *emcalLoader = 0;
  if(runLoader) {
    emcalLoader = dynamic_cast<AliEMCALLoader*>(runLoader->GetDetectorLoader("EMCAL"));
  }
 
  //Load EMCAL Geometry
  if (runLoader && runLoader->GetAliRun()){
    AliEMCAL* emcal = dynamic_cast<AliEMCAL*>(runLoader->GetAliRun()->GetDetector("EMCAL"));
    if(emcal)fGeom = emcal->GetGeometry();
  }
  
  if (!fGeom) 	 
    fGeom = AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName()); 	 

  if (!fGeom)
    AliFatal("Did not get geometry from EMCALLoader");
  
  //Define parameters
  Int_t nSuperModules = fGeom->GetNumberOfSuperModules() ; //12 SM in EMCAL
  Int_t nTRU       = fGeom->GetNTRU();    // 3 TRU per super module

  //Intialize data members each time the trigger is called in event loop
  f2x2MaxAmp = -1; f2x2ModulePhi = -1;  f2x2ModuleEta = -1;
  fnxnMaxAmp = -1; fnxnModulePhi = -1;  fnxnModuleEta = -1;

  // Take the digits list if simulation
  if(fSimulation && runLoader && emcalLoader){ // works than run seperate macros
    runLoader->LoadDigits("EMCAL");
    fDigitsList = emcalLoader->Digits() ;
    runLoader->LoadSDigits("EMCAL");
  }
  // Digits list should be set by method SetDigitsList(TClonesArray * digits)
  if(!fDigitsList)
    AliFatal("Digits not found !") ;
  
  //Take the digits list 
  
  // Delete old if unzero
  if(fAmpTrus)     {fAmpTrus->Delete();   delete fAmpTrus;}
  if(fTimeRtrus)   {fTimeRtrus->Delete(); delete fTimeRtrus;}
  if(fAmpSMods)    {fAmpSMods->Delete();  delete fAmpSMods;}
  // Fill TRU and SM matrix    
  fAmpTrus   = new TClonesArray("TMatrixD",nTRU);
  fAmpTrus->SetName("AmpTrus");
  fTimeRtrus = new TClonesArray("TMatrixD",nTRU);
  fTimeRtrus->SetName("TimeRtrus");
  fAmpSMods  = new TClonesArray("TMatrixD",nSuperModules);
  fAmpSMods->SetName("AmpSMods");
  
  FillTRU(fDigitsList, fAmpTrus, fAmpSMods, fTimeRtrus);

  // Jet stuff - only one case, no freedom here
  if(fGeom->GetNEtaSubOfTRU() == 6) {
    if(fAmpJetMatrix) {delete fAmpJetMatrix; fAmpJetMatrix=0;}
    if(fJetMatrixE)   {delete fJetMatrixE; fJetMatrixE=0;}

    fAmpJetMatrix = new TMatrixD(17,12); // 17-phi(row), 12-eta(col)
    fJetMatrixE = new TH2F("fJetMatrixE"," E of max patch in (#phi,#eta)", 
    17, 80.*TMath::DegToRad(), (180.+20.*2/3.)*TMath::DegToRad(), 12, -0.7, 0.7);
    for(Int_t row=0; row<fAmpJetMatrix->GetNrows(); row++) {
      for(Int_t col=0; col<fAmpJetMatrix->GetNcols(); col++) {
        (*fAmpJetMatrix)(row,col) = 0.;
      }
    }
    FillJetMatrixFromSMs(fAmpSMods, fAmpJetMatrix, fGeom);
  }
  if(!CheckConsistentOfMatrixes()) assert(0);

  // Do Tower Sliding and select Trigger
  // Initialize varible that will contain maximum amplitudes and 
  // its corresponding tower position in eta and phi, and time.
  TMatrixD ampmax2(4,nTRU) ; // 0-max amp, 1-irow, 2-icol, 3-timeR
  TMatrixD ampmaxn(4,nTRU) ;
  
  for(Int_t iSM = 0 ; iSM < nSuperModules ; iSM++) {
    //Do 2x2 and nxn sums, select maximums. 

    MakeSlidingTowers(fAmpTrus, fTimeRtrus, iSM, ampmax2, ampmaxn);
    
    // Set the trigger
    if(fIsolateInSuperModule) // here some discripency between tru and SM
      SetTriggers(fAmpSMods,iSM,ampmax2,ampmaxn) ;
    if(!fIsolateInSuperModule)
      SetTriggers(fAmpTrus,iSM,ampmax2,ampmaxn) ;
  }
  
  // Do patch sliding and select Jet Trigger
  // 0-max amp-meanFromVZERO(if), 1-irow, 2-icol, 3-timeR, 
  // 4-max amp , 5-meanFromVZERO (Nov 25, 2007)
  // fAmpJetMax(6,1)
  MakeSlidingPatch((*fAmpJetMatrix), fNJetPatchPhi, fAmpJetMax); // no timing information here

  //Print();
  // fDigitsList = 0;
}

//____________________________________________________________________________
void AliEMCALTrigger::GetTriggerInfo(TArrayF &triggerPosition, TArrayF &triggerAmplitudes) const
{
  // Template - should be defined; Nov 5, 2007
  triggerPosition[0]   = 0.; 
  triggerAmplitudes[0] = 0.;
}

//____________________________________________________________________________
void AliEMCALTrigger::FillJetMatrixFromSMs(TClonesArray *ampmatrixsmod, TMatrixD* jetMat, AliEMCALGeometry *g)
{
  // Nov 5, 2007
  // Fill matrix for jet trigger from SM matrixes of modules
  //
  static int keyPrint = 0;

  if(ampmatrixsmod==0 || jetMat==0 || g==0) return;
  Double_t amp = 0.0, ampSum=0.0;

  Int_t nEtaModSum = g->GetNZ()   / g->GetNEtaSubOfTRU(); // should be 4 
  Int_t nPhiModSum = g->GetNPhi() / g->GetNTRUPhi();      // should be 4 

  if(keyPrint) AliDebug(2,Form("%s",Form(" AliEMCALTrigger::FillJetMatrixFromSMs | nEtaModSum %i : nPhiModSum %i \n", nEtaModSum, nPhiModSum)));
  Int_t jrow=0, jcol=0;  // indexes of jet matrix
  Int_t nEtaSM=0, nPhiSM=0;
  for(Int_t iSM=0; iSM<ampmatrixsmod->GetEntries(); iSM++) {
    TMatrixD * ampsmods   = dynamic_cast<TMatrixD *>(ampmatrixsmod->At(iSM));
    
    if(!ampsmods) return;
    
    Int_t nrow = ampsmods->GetNrows();
    Int_t ncol = ampsmods->GetNcols();
    //printf("%s",Form(" ######## SM %i : nrow %i : ncol %i ##### \n", iSM, nrow, ncol));
    for(Int_t row=0; row<nrow; row++) {
      for(Int_t col=0; col<ncol; col++) {
        amp  = (*ampsmods)(row,col);
        nPhiSM = iSM / 2; 
        nEtaSM = iSM % 2;
        if       (amp>0.0) {
           if(keyPrint) AliDebug(2,Form("%s",Form(" ** nPhiSm %i : nEtaSM %i : row %2.2i : col %2.2i -> ", nPhiSM, nEtaSM, row, col))); 
          if(nEtaSM == 0) { // positive Z
            jrow = 3*nPhiSM + row/nPhiModSum;
            jcol = 6 + col / nEtaModSum;
	  } else {         // negative Z
            if(iSM<=9) jrow = 3*nPhiSM + 2 - row/nPhiModSum;
            else       jrow = 3*nPhiSM + 1 - row/nPhiModSum; // half size 
            jcol = 5 - col / nEtaModSum;
	  }
	  if(keyPrint) AliDebug(2,Form("%s",Form(" jrow %2.2i : jcol %2.2i : amp %f (jetMat) \n", jrow, jcol, amp))); 

          (*jetMat)(jrow,jcol) += amp;
          ampSum += amp; // For controling
        } else if(amp<0.0) {
          AliDebug(1,Form(" jrow %2.2i : jcol %2.2i : amp %f (jetMat: amp<0) \n", jrow, jcol, amp)); 
	  assert(0);
        }
      }
    }
  } // cycle on SM
  if(ampSum <= 0.0) AliDebug(1,Form("ampSum %f (<=0.0) ", ampSum));
}

//____________________________________________________________________________
void AliEMCALTrigger::MakeSlidingPatch(const TMatrixD &jm, const Int_t nPatchSize, TMatrixD &ampJetMax)
{
  // Sliding patch : nPatchSize x nPatchSize (OVERLAP)
  static int keyPrint = 0;
  if(keyPrint) AliDebug(2,Form(" AliEMCALTrigger::MakeSlidingPatch() was started \n"));
  Double_t ampCur = 0.0, e=0.0;
  ampJetMax(0,0)  = 0.0;
  ampJetMax(3,0)  = 0.0; // unused now
  ampJetMax(4,0)  = ampJetMax(5,0)  = 0.0;
  for(Int_t row=0; row<fAmpJetMatrix->GetNrows(); row ++) {
    for(Int_t col=0; col<fAmpJetMatrix->GetNcols(); col++) {
      ampCur = 0.;
      // check on patch size
      if( (row+nPatchSize-1) <  fAmpJetMatrix->GetNrows() && (col+nPatchSize-1) < fAmpJetMatrix->GetNcols()){
        for(Int_t i = 0 ; i < nPatchSize ; i++) {
	  for(Int_t j = 0 ; j < nPatchSize ; j++) {
            ampCur += jm(row+i, col+j);
          }
        } // end cycle on patch
        if(ampCur > ampJetMax(0,0)){
	  ampJetMax(0,0) = ampCur; 
	  ampJetMax(1,0) = row;
	  ampJetMax(2,0) = col;
        }
      } // check on patch size
    }
  }
  if(keyPrint) AliDebug(2,Form(" ampJetMax %i row %2i->%2i col %2i->%2i \n", Int_t(ampJetMax(0,0)), Int_t(ampJetMax(1,0)), Int_t(ampJetMax(1,0))+nPatchSize-1, Int_t(ampJetMax(2,0)), Int_t(ampJetMax(2,0))+nPatchSize-1));

  Double_t eCorrJetMatrix=0.0;
  if(fVZER0Mult > 0.0) {
  // Correct patch energy (adc) and jet patch matrix energy
    Double_t meanAmpBG = GetMeanEmcalPatchEnergy(Int_t(fVZER0Mult), nPatchSize)/0.0153;
    ampJetMax(4,0) = ampJetMax(0,0);
    ampJetMax(5,0) = meanAmpBG;

    Double_t eCorr     = ampJetMax(0,0) - meanAmpBG; 
    AliDebug(2,Form(" ampJetMax(0,0) %f meanAmpBG %f  eCorr %f : ampJetMax(4,0) %f \n",
	   ampJetMax(0,0), meanAmpBG, eCorr, ampJetMax(5,0)));
    ampJetMax(0,0)     = eCorr;
    // --
    eCorrJetMatrix = GetMeanEmcalEnergy(Int_t(fVZER0Mult)) / 208.;
  }
  // Fill patch energy matrix
  for(int row=Int_t(ampJetMax(1,0)); row<Int_t(ampJetMax(1,0))+nPatchSize; row++) {
    for(int col=Int_t(ampJetMax(2,0)); col<Int_t(ampJetMax(2,0))+nPatchSize; col++) {
      e = Double_t(jm(row,col)*0.0153);  // 0.0153 - hard coded now
      if(eCorrJetMatrix > 0.0) { // BG subtraction case
        e -= eCorrJetMatrix;
        fJetMatrixE->SetBinContent(row+1, col+1, e);
      } else if(e > 0.0) {
        fJetMatrixE->SetBinContent(row+1, col+1, e);
      }
    }
  }
  // PrintJetMatrix();
  // Set the jet trigger(s), multiple threshold now, Nov 19,2007
  for(Int_t i=0; i<fNJetThreshold; i++ ) {
    if(ampJetMax(0,0) >= fL1JetThreshold[i]) {
      SetInput(GetNameOfJetTrigger(i)); 
    }
  }
}

//____________________________________________________________________________
Double_t AliEMCALTrigger::GetEmcalSumAmp() const 
{ 
  // Return sum of amplidutes from EMCal
  // Used calibration coefficeint for transition to energy
  return fAmpJetMatrix >0 ?fAmpJetMatrix->Sum() :0.0;
}

//____________________________________________________________________________
void AliEMCALTrigger::PrintJetMatrix() const
{
  //  fAmpJetMatrix : (17,12); // 17-phi(row), 12-eta(col)
  if(fAmpJetMatrix == 0) return;

  AliInfo(Form("\n ####  jetMatrix : (%i,%i) ##### \n ", 
	       fAmpJetMatrix->GetNrows(), fAmpJetMatrix->GetNcols()));
  PrintMatrix(*fAmpJetMatrix);
}

//____________________________________________________________________________
void AliEMCALTrigger::PrintAmpTruMatrix(Int_t ind) const
{
 // Print matrix with TRU patches
  TMatrixD * tru = dynamic_cast<TMatrixD *>(fAmpTrus->At(ind));
  if(tru == 0) return;
  AliInfo(Form("\n ####  Amp TRU matrix(%i) : (%i,%i) ##### \n ", 
	 ind, tru->GetNrows(), tru->GetNcols()));
  PrintMatrix(*tru);
}

//____________________________________________________________________________
void AliEMCALTrigger::PrintAmpSmMatrix(Int_t ind) const
{
	// Print matrix with SM amplitudes
	TMatrixD * sm = dynamic_cast<TMatrixD *>(fAmpSMods->At(ind));
  if(sm == 0) return;
  AliInfo(Form("\n ####  Amp SM matrix(%i) : (%i,%i) ##### \n ", 
	 ind, sm->GetNrows(), sm->GetNcols()));
  PrintMatrix(*sm);
}

//____________________________________________________________________________
void AliEMCALTrigger::PrintMatrix(const TMatrixD &mat) const
{
  //Print matrix object
  for(Int_t col=0; col<mat.GetNcols(); col++) AliInfo(Form(" %3i ", col));
  AliInfo(Form("\n -- \n"));
  for(Int_t row=0; row<mat.GetNrows(); row++) {
    AliInfo(Form(" row:%2i ", row));
    for(Int_t col=0; col<mat.GetNcols(); col++) {
      AliInfo(Form(" %4i", (Int_t)mat(row,col)));
    }
    AliInfo("\n");
  }
}

//____________________________________________________________________________
Bool_t AliEMCALTrigger::CheckConsistentOfMatrixes(const Int_t pri)
{
  // Check consitency of matrices
  Double_t sumSM = 0.0, smCur=0.0;
  Double_t sumTru=0.0,  sumTruInSM = 0.0, truSum=0.0;
  //  Bool_t key = kTRUE;
 
  for(Int_t i=0; i<fAmpSMods->GetEntries(); i++) {
    TMatrixD * sm = dynamic_cast<TMatrixD *>(fAmpSMods->At(i));
    if(sm) {
      smCur  = sm->Sum();
      sumSM += smCur;

      sumTruInSM = 0.0;
      for(Int_t itru=0; itru<3; itru++) { // Cycle on tru inside SM
        Int_t ind = 3*i + itru;
        TMatrixD *tru = dynamic_cast<TMatrixD *>(fAmpTrus->At(ind)); 
        if(tru) {
          truSum = tru->Sum();
	  sumTruInSM += truSum;
        }
      }
      sumTru += sumTruInSM;

      if(sumTruInSM != smCur) {
        AliDebug(1,Form(" sm %i : smCur %f -> sumTruInSM %f \n", i, smCur, sumTruInSM));
        return kFALSE;
      }
    }
  }
  Double_t sumJetMat = fAmpJetMatrix->Sum();
	if(pri || TMath::Abs(sumSM-sumTru)>0.0001 || TMath::Abs(sumSM-sumJetMat) > 0.0001) 
   AliDebug(1,Form(" sumSM %f : sumTru %f : sumJetMat %f \n", sumSM, sumTru, sumJetMat)); 
	if(TMath::Abs(sumSM - sumTru)>0.0001 || TMath::Abs(sumSM-sumJetMat) > 0.0001) return kFALSE; 
  else                                       return kTRUE; 
}

//____________________________________________________________________________
void AliEMCALTrigger::Browse(TBrowser* b)
{
  //Browse.
  if(&fInputs)      b->Add(&fInputs);
  if(fAmpTrus)      b->Add(fAmpTrus);
  if(fTimeRtrus)    b->Add(fTimeRtrus);
  if(fAmpSMods)     b->Add(fAmpSMods);
  if(fAmpJetMatrix) b->Add(fAmpJetMatrix);
  if(fJetMatrixE)   b->Add(fJetMatrixE);
  //  if(c) b->Add(c);
}
