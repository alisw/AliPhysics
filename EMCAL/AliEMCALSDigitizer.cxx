/*************************************************************************
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
// This is a TTask that makes SDigits out of Hits
// A Summable Digits is the sum of all hits originating 
// from one in one tower of the EMCAL 
// A threshold for assignment of the primary to SDigit is applied 
// SDigits are written to TreeS, branch "EMCAL"
// AliEMCALSDigitizer with all current parameters is written 
// to TreeS branch "AliEMCALSDigitizer".
// Both branches have the same title. If necessary one can produce 
// another set of SDigits with different parameters. Two versions
// can be distunguished using titles of the branches.
// User case:
// root [0] AliEMCALSDigitizer * s = new AliEMCALSDigitizer("galice.root")
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
// root [1] s->ExecuteTask()
//             // Makes SDigitis for all events stored in galice.root
// root [2] s->SetPedestalParameter(0.001)
//             // One can change parameters of digitization
// root [3] s->SetSDigitsBranch("Redestal 0.001")
//             // and write them into the new branch
// root [4] s->ExeciteTask("deb all tim")
//             // available parameters:
//             deb - print # of produced SDigitis
//             deb all  - print # and list of produced SDigits
//             tim - print benchmarking information
//
//*-- Author : Sahal Yacoob (LBL)
// based on  : AliPHOSSDigitzer 
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TBenchmark.h"
	       //#include "TObjectTable.h"

// --- Standard library ---
#include "stdlib.h"

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGetter.h"
#include "AliEMCALHit.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"
	       //#include "AliMemoryWatcher.h"

ClassImp(AliEMCALSDigitizer)
           
//____________________________________________________________________________ 
  AliEMCALSDigitizer::AliEMCALSDigitizer():TTask("","") 
{
  // ctor
  fFirstEvent = fLastEvent  = 0 ;  
  fDefaultInit = kTRUE ; 
}

//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer(const char * alirunFileName, 
				       const char * eventFolderName):
  TTask("EMCAL"+AliConfig::Instance()->GetSDigitizerTaskName(), alirunFileName),
  fEventFolderName(eventFolderName)
{
  // ctor
  fFirstEvent = fLastEvent  = 0 ; // runs one event by defaut  
  Init();
  InitParameters() ; 
  fDefaultInit = kFALSE ; 
}


//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer(const AliEMCALSDigitizer & sd) : TTask(sd) {
  //cpy ctor 

  fFirstEvent    = sd.fFirstEvent ; 
  fLastEvent     = sd.fLastEvent ;
  fA             = sd.fA ;
  fB             = sd.fB ;
  fECPrimThreshold  = sd.fECPrimThreshold ;
  fSDigitsInRun  = sd.fSDigitsInRun ;
  SetName(sd.GetName()) ; 
  SetTitle(sd.GetTitle()) ; 
  fEventFolderName = sd.fEventFolderName;
}


//____________________________________________________________________________ 
AliEMCALSDigitizer::~AliEMCALSDigitizer() {
  // dtor
  AliEMCALGetter * gime = 
    //   AliEMCALGetter::Instance(GetTitle(), fEventFolderName.Data());  
    AliEMCALGetter::Instance();  
  gime->EmcalLoader()->CleanSDigitizer();
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::Init(){
  // Initialization: open root-file, allocate arrays for hits and sdigits,
  // attach task SDigitizer to the list of EMCAL tasks
  // 
  // Initialization can not be done in the default constructor
  //============================================================= YS
  //  The initialisation is now done by the getter

  fInit = kTRUE ; 
   
  AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName.Data());  
  if ( gime == 0 ) {
    Fatal("Init", "Could not obtain the Getter objectfor file %s and event %s !", GetTitle(), fEventFolderName.Data()) ;  
    return ;
  } 
  
  TString opt("SDigits") ; 
  if(gime->VersionExists(opt) ) { 
    Error( "Init", "Give a version name different from %s", fEventFolderName.Data() ) ;
    fInit = kFALSE ; 
  }
  
  gime->PostSDigitizer(this);
  gime->EmcalLoader()->GetSDigitsDataLoader()->GetBaseTaskLoader()->SetDoNotReload(kTRUE);
  
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::InitParameters()
{
  // Initializes parameters
  fA                      = 0;
  fB                      = 10000000.;

  AliEMCALGetter * gime = AliEMCALGetter::Instance() ;
  const AliEMCALGeometry * geom = gime->EMCALGeometry() ; 
  if (geom->GetSampling() == 0.) {
    Fatal("InitParameters", "Sampling factor not set !") ; 
  }
//   else
//     Info("InitParameters", "Sampling factor set to %f", geom->GetSampling()) ; 
  
  // this threshold corresponds approximately to 100 MeV
  fECPrimThreshold     = 100E-3;
}

//____________________________________________________________________________
void AliEMCALSDigitizer::Exec(Option_t *option) 
{ 
  // Collects all hit of the same tower into digits
  if (strstr(option, "print") ) {
    Print() ; 
    return ; 
  }
  
  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALSDigitizer");

  //AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle()) ;
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ;

  //switch off reloading of this task while getting event
  if (!fInit) { // to prevent overwrite existing file
    AliError( Form("Give a version name different from %s", fEventFolderName.Data()) ) ;
    return ;
    }

  if (fLastEvent == -1) 
    fLastEvent = gime->MaxEvent() - 1 ;
  else 
    fLastEvent = TMath::Min(fFirstEvent, gime->MaxEvent()); // only ine event at the time
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;

  Int_t ievent ;   

  //AliMemoryWatcher memwatcher;

  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    gime->Event(ievent,"H") ;  
    TTree * treeS = gime->TreeS(); 
    TClonesArray * hits = gime->Hits() ; 
    TClonesArray * sdigits = gime->SDigits() ;
    sdigits->Clear();
    Int_t nSdigits = 0 ;
    //Now make SDigits from hits, for EMCAL it is the same, so just copy    
    Int_t nPrim =  static_cast<Int_t>((gime->TreeH())->GetEntries()) ; 
    // Attention nPrim is the number of primaries tracked by Geant 
    // and this number could be different to the number of Primaries in TreeK;
    Int_t iprim ;
    for ( iprim = 0 ; iprim < nPrim ; iprim++ ) { 
      //=========== Get the EMCAL branch from Hits Tree for the Primary iprim
      gime->Track(iprim) ;
      Int_t i;
      for ( i = 0 ; i < hits->GetEntries() ; i++ ) {
	AliEMCALHit * hit = dynamic_cast<AliEMCALHit*>(hits->At(i)) ;
	AliEMCALDigit * curSDigit = 0 ;
	AliEMCALDigit * sdigit = 0 ;
	Bool_t newsdigit = kTRUE; 

	// Assign primary number only if deposited energy is significant
	AliEMCALGeometry * geom = gime->EMCALGeometry() ; 
	if( geom->IsInECA(hit->GetId()) ){
	  if( hit->GetEnergy() >  fECPrimThreshold )
	    curSDigit =  new AliEMCALDigit( hit->GetPrimary(),
					    hit->GetIparent(), hit->GetId(), 
					    Digitize(hit->GetEnergy()), hit->GetTime() ) ;
	  else
	    curSDigit =  new AliEMCALDigit( -1               , 
					    -1               ,
					    hit->GetId(), 
					    Digitize(hit->GetEnergy()), hit->GetTime() ) ;
	}
	else{
	  newsdigit = kFALSE;
	}

	Int_t check = 0 ;

	for(check= 0; check < nSdigits ; check++) {
	  sdigit = dynamic_cast<AliEMCALDigit *>(sdigits->At(check)) ;
          if(curSDigit != 0){
	    if( sdigit->GetId() == curSDigit->GetId()) { // Are we in the same ECAL tower ?              
	    *sdigit = *sdigit + *curSDigit;
	    newsdigit = kFALSE;
	    }
	  }
	}
	if (newsdigit) {
	  new((*sdigits)[nSdigits])  AliEMCALDigit(*curSDigit);
	  nSdigits++ ;  
	}
	delete curSDigit ; 
      }  // loop over all hits (hit = deposited energy/entering particle)
    } // loop over iprim
    sdigits->Sort() ;
    
    nSdigits = sdigits->GetEntriesFast() ;
    fSDigitsInRun += nSdigits ;  
    sdigits->Expand(nSdigits) ;

    Int_t nPrimarymax = -1 ; 
    Int_t i ;
    for (i = 0 ; i < sdigits->GetEntriesFast() ; i++) { 
      AliEMCALDigit * sdigit = dynamic_cast<AliEMCALDigit *>(sdigits->At(i)) ;
      sdigit->SetIndexInList(i) ;
    }
    
    for (i = 0 ; i < sdigits->GetEntriesFast() ; i++) {   
      if (((dynamic_cast<AliEMCALDigit *>(sdigits->At(i)))->GetNprimary()) > nPrimarymax)
	nPrimarymax = ((dynamic_cast<AliEMCALDigit *>(sdigits->At(i)))->GetNprimary()) ;
    }
    
    // Now write SDigits    
    //First list of sdigits

    Int_t bufferSize = 32000 ;    
    TBranch * sdigitsBranch = treeS->Branch("EMCAL",&sdigits,bufferSize);
 
    sdigitsBranch->Fill() ;

    gime->WriteSDigits("OVERWRITE");
    
    //NEXT - SDigitizer

    gime->WriteSDigitizer("OVERWRITE");
    
    if(strstr(option,"deb"))
      PrintSDigits(option) ;  
    
    //gObjectTable->Print() ; 
    //memwatcher.Watch(ievent); 
  }// event loop

  Unload();
  
  gime->EmcalLoader()->GetSDigitsDataLoader()->GetBaseTaskLoader()->SetDoNotReload(kTRUE);
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALSDigitizer"); 
    printf("Exec: took %f seconds for SDigitizing %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALSDigitizer"), gBenchmark->GetCpuTime("EMCALSDigitizer")/nEvents ) ; 
  }

  //TFile f("out.root","RECREATE");
  //memwatcher.WriteToFile(); 
  //f.Close();
}


//__________________________________________________________________
void AliEMCALSDigitizer::Print()const
{ 
  // Prints parameters of SDigitizer
  printf("Print: \n------------------- %s -------------", GetName() ) ; 
  printf("   Writing SDigits to branch with title  %s\n", fEventFolderName.Data()) ;
  printf("   with digitization parameters  A = %f\n", fA) ; 
  printf("                                 B = %f\n", fB) ;
  printf("   Threshold for EC Primary assignment= %f\n", fECPrimThreshold)  ; 
  printf("---------------------------------------------------\n") ;

}

//__________________________________________________________________
Bool_t AliEMCALSDigitizer::operator==( AliEMCALSDigitizer const &sd )const
{
  // Equal operator.
  // SDititizers are equal if their pedestal, slope and threshold are equal
  if( (fA==sd.fA)&&(fB==sd.fB)&&
      (fECPrimThreshold==sd.fECPrimThreshold))
    return kTRUE ;
  else
    return kFALSE ;
}

//__________________________________________________________________
void AliEMCALSDigitizer::PrintSDigits(Option_t * option)
{
  //Prints list of digits produced at the current pass of AliEMCALDigitizer
  
  
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  const TClonesArray * sdigits = gime->SDigits() ; 
  
  printf("\n") ;  
  printf("event %i", gAlice->GetEvNumber()) ;
  printf("\n      Number of entries in SDigits list %i", sdigits->GetEntriesFast()); 
  if(strstr(option,"all")||strstr(option,"EMC")){
    
    //loop over digits
    AliEMCALDigit * digit;
    printf("\n   Id  Amplitude    Time          Index Nprim: Primaries list \n") ;    
    Int_t index ;
    char * tempo = new char[8192]; 
    for (index = 0 ; index < sdigits->GetEntries() ; index++) {
      digit = dynamic_cast<AliEMCALDigit *>( sdigits->At(index) ) ;
      sprintf(tempo, "\n%6d  %8d    %6.5e %4d      %2d :",
	      digit->GetId(), digit->GetAmp(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
      printf(tempo); 
      
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
	sprintf(tempo, "%d ",digit->GetPrimary(iprimary+1) ) ; 
	printf(tempo); 
      }  	 
    }
    delete tempo ;
  }
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::Unload() const
{
  // Unload Hits and SDigits from the folder
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  AliEMCALLoader * loader = gime->EmcalLoader() ; 
  loader->UnloadHits() ; 
  loader->UnloadSDigits() ; 
}
