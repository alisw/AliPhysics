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
//  Digits are grouped in TRU's (384 cells? ordered fNTRUPhi x fNTRUEta). 
//  The algorithm searches all possible 4x4 cell combinations per each TRU, 
//  adding the digits amplitude and finding the maximum. Maximums are compared 
//  to triggers threshold and they are set. Thresholds need to be fixed. 
//  Last 2 modules are half size but they are treated as fullsize, then their 
//  TRU should be smaller. When this is fixed, class will be updated. 
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
//  tr->Print(""); //Print results
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TMatrixD.h"

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
    fL0MBPbPbThreshold(500), fL0MBppThreshold(50),fL1JetLowPtThreshold(1000), 
    fL1JetMediumPtThreshold(10000), fL1JetHighPtThreshold(20000) 
{
  //ctor

   SetName("EMCAL");
   CreateInputs();
   
   //Print("all") ; 
}



//____________________________________________________________________________
AliEMCALTrigger::AliEMCALTrigger(const AliEMCALTrigger & trig) 
  : AliTriggerDetector(trig) 
{

  // cpy ctor
  
  fL0MBPbPbThreshold       = trig.fL0MBPbPbThreshold ; 
  fL0MBppThreshold         = trig.fL0MBppThreshold ; 
  fL1JetLowPtThreshold     = trig.fL1JetLowPtThreshold ;
  fL1JetMediumPtThreshold  = trig.fL1JetMediumPtThreshold ;
  fL1JetHighPtThreshold    = trig.fL1JetHighPtThreshold ;

}

//----------------------------------------------------------------------
void AliEMCALTrigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "EMCAL_MB_PbPb_L0",     "EMCAL PbPb Minimum Bias L0",  0x01 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_MB_pp_L0",       "EMCAL pp Minimum Bias L0",    0x02 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_PbPb_JetHPt_L1", "EMCAL PbPb Jet High Pt L1",   0x04 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_PbPb_JetMPt_L1", "EMCAL PbPb Jet Medium Pt L1", 0x08 ) );
   fInputs.AddLast( new AliTriggerInput( "EMCAL_PbPb_JetLPt_L1", "EMCAL PbPb Jet Low Pt L1",    0x016 ) );
 
}

//____________________________________________________________________________
void AliEMCALTrigger::MakeSlidingCell(const TClonesArray * trus , 
				      const Int_t isupermod,
				      const Int_t nTRU,
				      const Int_t nCellsPhi,
				      const Int_t nCellsEta,
				      Float_t *ampmax){
  
  //Sums energy of all possible 4x4 cells per each TRU. Fast signal 
  //in the experiment is given by 2x2 cells, for this reason we loop 
  //inside the TRU cells by 2. 
  
  Float_t amp = 0 ;
  
  for(Int_t i = 0 ; i < nTRU ; i++)   
    ampmax[i] = 0 ;
  
  //Loop over all TRUS in the seleted supermodule
  for(Int_t itru = 0 + (isupermod - 1) * nTRU ; itru < isupermod*nTRU ; itru++)
    {
      TMatrixD * tru = dynamic_cast<TMatrixD *>(trus->At(itru)) ;
      
      //Sliding 2x2 cell       
      for(Int_t irow = 0 ; irow <  nCellsPhi; irow += 2){ 
	for(Int_t icol = 0 ; icol < nCellsEta ; icol += 2){
	  amp = 0;
	  if( (irow+3) < nCellsPhi && (icol+3) < nCellsEta){//Avoid exit the TRU
	    for(Int_t i = irow; i < irow + 4 ; i++){
	      for(Int_t j = icol; j < icol + 4 ; j++){
		amp += (*tru)(i,j) ;
	      }
	    }
	  }
	  if(amp > ampmax[itru-(isupermod-1)*nTRU])
	    ampmax[itru-(isupermod-1)*nTRU] = amp ;
	  
	}
      }
    }
}

//____________________________________________________________________________
void AliEMCALTrigger::Trigger() 
{
  
  //Main Method to select triggers.
  //Loader
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));
  
  //Load EMCAL Geometry
  rl->LoadgAlice(); 
  AliRun * gAlice = rl->GetAliRun(); 
  AliEMCAL * emcal  = (AliEMCAL*)gAlice->GetDetector("EMCAL");
  AliEMCALGeometry * geom = emcal->GetGeometry();
  //AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance(); 

  if (geom==0)
    AliFatal("Did not get geometry from EMCALLoader");
 
  //Define some useful parameters
 
  Int_t nSuperModules = geom->GetNumberOfSuperModules() ; //12 SM in EMCAL
  Int_t nTRU          = geom->GetNTRU();//3 TRU per super module
  Int_t nCellsPhi  = geom->GetNPhi()*2/geom->GetNTRUPhi() ;
  // 12(tow)*2(cell)/1 TRU, cells in Phi in one TRU
  Int_t nCellsEta  = geom->GetNEta()*2/geom->GetNTRUEta() ;
  // 24(mod)*2(tower)/3 TRU, cells in Eta in one TRU

  //Info("Trigger","nSuperModules %d, nCellsPhi %d, nCellsEta %d",
  //     nSuperModules, nCellsPhi,nCellsEta);
  
  //Take the digits list and declare digits pointers
  TClonesArray * digits = emcalLoader->Digits(); //gime->Digits() ;
  
  TClonesArray * trus = geom->FillTRU(digits) ;

  //Do Cell Sliding and select Trigger
  Float_t max [10] ;
  for(Int_t iSM = 1 ; iSM <= nSuperModules ; iSM++) {
    
    MakeSlidingCell(trus, iSM, nTRU, nCellsPhi, nCellsEta, max);

    //cout<<"Max Amplitude in SuperMod "<<iSM<<" TRU1 "<<max[0]
    //<<" TRU2 "<<max[1]<<" TRU3 "<<max[2]<<endl;
    SetTriggers(max, nTRU) ;
  }
  
}

//____________________________________________________________________________
void AliEMCALTrigger::Print(const Option_t * opt) const 
{
  
  //Prints main parameters
  
  if(! opt)
    return;
  AliTriggerInput* in = 0x0 ;
  
  
  AliInfo("EMCAL trigger information:") ; 
  printf( "             Threshold for pp MB LO %d\n", fL0MBppThreshold) ;  
  in = (AliTriggerInput*)fInputs.FindObject( "EMCAL_MB_pp_L0" );
  if(in->GetValue())
    printf( "             *** EMCAL MB pp LO is set ***\n") ; 
  
  printf( "             Threshold for PbPb MB LO %d\n", fL0MBPbPbThreshold) ;  
  in = (AliTriggerInput*)fInputs.FindObject( "EMCAL_MB_PbPb_L0" );
  if(in->GetValue())
    printf( "             *** EMCAL MB PbPb LO is set ***\n") ; 
  
  printf( "             Jet Low Pt Threshold for PbPb L1 %d\n", fL1JetLowPtThreshold) ;
  in = (AliTriggerInput*)fInputs.FindObject( "EMCAL_PbPb_JetLPt_L1" );
  if(in->GetValue())
    printf( "             *** EMCAL Jet Low Pt for PbPb L1 is set ***\n") ;

  printf( "             Jet Medium Pt Threshold for L1 %d\n", fL1JetMediumPtThreshold) ;  
  in = (AliTriggerInput*) fInputs.FindObject( "EMCAL_PbPb_JetMPt_L1" );
  if(in->GetValue())
    printf( "             *** EMCAL Jet Medium Pt for PbPb L1 is set ***\n") ;

  printf( "             Jet High Pt Threshold for L1 %d\n", fL1JetHighPtThreshold) ;  
  in = (AliTriggerInput*) fInputs.FindObject( "EMCAL_PbPb_JetHPt_L1" );
  if(in->GetValue())
    printf( "             *** EMCAL Jet High Pt for PbPb L1 is set ***\n") ;

}

//____________________________________________________________________________
void AliEMCALTrigger::SetTriggers(const Float_t * amp, const Int_t nTRU)  
{

  //Checks the maximum amplitude per each TRU and compares with the 
  //different triggers thresholds

  Float_t max = 0;
  for(Int_t i = 0 ; i < nTRU ; i++){
    if(max < amp[i] )
      max = amp[i] ;
  }
 
  if(max >= fL0MBppThreshold)
    SetInput("EMCAL_MB_pp_L0");
  if(max >= fL0MBPbPbThreshold)
    SetInput("EMCAL_MB_PbPb_L0");
  if(max >= fL1JetLowPtThreshold)
    SetInput("EMCAL_PbPb_JetLPt_L1"); 
  if(max >= fL1JetMediumPtThreshold)
    SetInput("EMCAL_PbPb_JetMPt_L1"); 
  if(max >= fL1JetHighPtThreshold)
    SetInput("EMCAL_PbPb_JetHPt_L1");

}
