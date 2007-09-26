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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.51  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.50  2006/08/28 10:01:56  kharlov
 * Effective C++ warnings fixed (Timur Pocheptsov)
 *
 * Revision 1.49  2006/05/10 06:42:53  kharlov
 * Remove redundant loop over primaries
 *
 * Revision 1.48  2006/04/22 10:30:17  hristov
 * Add fEnergy to AliPHOSDigit and operate with EMC amplitude in energy units (Yu.Kharlov)
 *
 * Revision 1.47  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// This is a TTask that makes SDigits out of Hits
// The name of the TTask is also the title of the branch that will contain 
// the created SDigits
// The title of the TTAsk is the name of the file that contains the hits from
// which the SDigits are created
// A Summable Digits is the sum of all hits originating 
// from one primary in one active cell
// A threshold for assignment of the primary to SDigit is applied 
// SDigits are written to TreeS, branch "PHOS"
// AliPHOSSDigitizer with all current parameters is written 
// to TreeS branch "AliPHOSSDigitizer".
// Both branches have the same title. If necessary one can produce 
// another set of SDigits with different parameters. Two versions
// can be distunguished using titles of the branches.
// User case:
//  root [0] AliPHOSSDigitizer * s = new AliPHOSSDigitizer("galice.root")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//  root [1] s->ExecuteTask()
//             // Makes SDigitis for all events stored in galice.root
//  root [2] s->SetPedestalParameter(0.001)
//             // One can change parameters of digitization
// root [3] s->SetSDigitsBranch("Pedestal 0.001")
//             // and write them into the new branch
// root [4] s->ExecuteTask("deb all tim")
//             // available parameters:
//             deb - print # of produced SDigitis
//             deb all  - print # and list of produced SDigits
//             tim - print benchmarking information
//
//*-- Author :  Dmitri Peressounko (SUBATECH & KI) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TBenchmark.h"
		//#include "TObjectTable.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSGeometry.h" 
#include "AliPHOSDigit.h"
#include "AliPHOSGetter.h"
#include "AliPHOSHit.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSQualAssDataMaker.h" 

ClassImp(AliPHOSSDigitizer)

           
//____________________________________________________________________________ 
AliPHOSSDigitizer::AliPHOSSDigitizer() : 
  TTask("",""),
  fA(0.f), fB(0.f),
  fPrimThreshold(0.f),
  fDefaultInit(kTRUE),
  fEventFolderName(""),
  fInit(kFALSE),
  fSDigitsInRun(0),
  fFirstEvent(0),
  fLastEvent(0), 
  fQADM (0x0)
{
  // ctor
  // Intialize the quality assurance data maker 	
}

//____________________________________________________________________________ 
AliPHOSSDigitizer::AliPHOSSDigitizer(const char * alirunFileName, 
				     const char * eventFolderName):
  TTask("PHOS"+AliConfig::Instance()->GetSDigitizerTaskName(), alirunFileName),
  fA(0.f), fB(0.f),
  fPrimThreshold(0.f),
  fDefaultInit(kFALSE),
  fEventFolderName(eventFolderName),
  fInit(kFALSE),
  fSDigitsInRun(0),
  fFirstEvent(0),
  fLastEvent(0), 
  fQADM (0x0)
{
  // ctor
  InitParameters() ; 
  Init();
  fDefaultInit = kFALSE ; 
  // Intialize the quality assurance data maker 	
  GetQualAssDataMaker()->Init(AliQualAss::kHITS) ;
  GetQualAssDataMaker()->Init(AliQualAss::kSDIGITS) ; 
}

//____________________________________________________________________________
AliPHOSSDigitizer::AliPHOSSDigitizer(const AliPHOSSDigitizer& sd) :
  TTask(sd.GetName(), sd.GetTitle()),
  fA(sd.fA), fB(sd.fB),
  fPrimThreshold(sd.fPrimThreshold),
  fDefaultInit(kFALSE),
  fEventFolderName(sd.fEventFolderName),
  fInit(kFALSE),
  fSDigitsInRun(sd.fSDigitsInRun),
  fFirstEvent(sd.fFirstEvent),
  fLastEvent(sd.fLastEvent), 
  fQADM (sd.fQADM)
{ 
  // cpy ctor
  // Intialize the quality assurance data maker 	
  GetQualAssDataMaker()->Init(AliQualAss::kHITS) ;
  GetQualAssDataMaker()->Init(AliQualAss::kSDIGITS) ; 
}


//_____________________________________________________________________________
AliPHOSSDigitizer& AliPHOSSDigitizer::operator = (const AliPHOSSDigitizer& qa)
{
// assignment operator

  this->~AliPHOSSDigitizer();
  new(this) AliPHOSSDigitizer(qa);
  return *this;
}

//____________________________________________________________________________ 
AliPHOSSDigitizer::~AliPHOSSDigitizer() {
  //dtor
  //  AliPHOSGetter * gime =
  //  AliPHOSGetter::Instance(GetTitle(),fEventFolderName.Data());  
  AliPHOSGetter * gime =
    AliPHOSGetter::Instance();  
  gime->PhosLoader()->CleanSDigitizer();
  delete fQADM ; 
}

//____________________________________________________________________________ 
void AliPHOSSDigitizer::Init()
{
  // Uses the getter to access the required files
  
  fInit = kTRUE ; 
  
  AliPHOSGetter * gime = AliPHOSGetter::Instance(GetTitle(), fEventFolderName.Data());  
  if ( gime == 0 ) {
    Fatal("Init" ,"Could not obtain the Getter object for file %s and event %s !", GetTitle(), fEventFolderName.Data()) ;  
    return ;
  } 
  
  TString opt("SDigits") ; 
  if(gime->VersionExists(opt) ) { 
    Error( "Init", "Give a version name different from %s", fEventFolderName.Data() ) ;
    fInit = kFALSE ; 
  }

  gime->PostSDigitizer(this);
  gime->PhosLoader()->GetSDigitsDataLoader()->GetBaseTaskLoader()->SetDoNotReload(kTRUE);
 
  fQADM = new AliPHOSQualAssDataMaker() ;  

}

//____________________________________________________________________________ 
void AliPHOSSDigitizer::InitParameters()
{ 
  // initializes the parameters for digitization
  fA             = 0;
  fB             = 10000000.;
  fPrimThreshold = 0.01 ;
  fSDigitsInRun  = 0 ;
}

//____________________________________________________________________________
void AliPHOSSDigitizer::Exec(Option_t *option) 
{ 
  // Steering method to produce summable digits for events
  // in the range from fFirstEvent to fLastEvent.
  // This range is optionally set by SetEventRange().
  // if fLastEvent=-1 (by default), then process events until the end.
  //
  // Summable digit is a sum of all hits in the same active
  // volume into digit
  
  if (strstr(option, "print") ) {
    Print() ; 
    return ; 
  }

  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSSDigitizer");
  
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ;

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
  
  Int_t ievent, i;

  //AliMemoryWatcher memwatcher;

  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {    
    gime->Event(ievent,"H") ;   
    TTree * treeS = gime->TreeS(); 
    TClonesArray * hits = gime->Hits() ;
    TClonesArray * sdigits = gime->SDigits() ;
    sdigits->Clear();
    Int_t nSdigits = 0 ;
    //Now make SDigits from hits, for PHOS it is the same, so just copy    
    for ( i = 0 ; i < hits->GetEntries() ; i++ ) {
      AliPHOSHit * hit = dynamic_cast<AliPHOSHit *>(hits->At(i)) ;
      // Assign primary number only if contribution is significant
      
      if( hit->GetEnergy() > fPrimThreshold)
	new((*sdigits)[nSdigits]) AliPHOSDigit(hit->GetPrimary(),hit->GetId(),
					       hit->GetEnergy() ,hit->GetTime()) ;
      else
	new((*sdigits)[nSdigits]) AliPHOSDigit(-1               ,hit->GetId(), 
					       hit->GetEnergy() ,hit->GetTime()) ;
      nSdigits++ ;	
      
    }
 
    sdigits->Sort() ;

    nSdigits = sdigits->GetEntriesFast() ;

    fSDigitsInRun += nSdigits ;  
    sdigits->Expand(nSdigits) ;

    for (i = 0 ; i < nSdigits ; i++) { 
      AliPHOSDigit * digit = dynamic_cast<AliPHOSDigit *>(sdigits->At(i)) ; 
      digit->SetIndexInList(i) ;     
    }

    // make Quality Assurance data

    GetQualAssDataMaker()->Exec(AliQualAss::kHITS, hits) ; 
    GetQualAssDataMaker()->Exec(AliQualAss::kSDIGITS, sdigits) ; 


    //Now write SDigits

    
    //First list of sdigits

    Int_t bufferSize = 32000 ;
    TBranch * sdigitsBranch = treeS->Branch("PHOS",&sdigits,bufferSize);

    sdigitsBranch->Fill() ;

    gime->WriteSDigits("OVERWRITE");

    //Next - SDigitizer

    gime->WriteSDigitizer("OVERWRITE");
    //gObjectTable->Print() ; 
  
    if(strstr(option,"deb"))
      PrintSDigits(option) ;
      
    //memwatcher.Watch(ievent); 
  }// event loop
  
  //Write the quality assurance data 
  GetQualAssDataMaker()->Finish(AliQualAss::kHITS) ;
  GetQualAssDataMaker()->Finish(AliQualAss::kSDIGITS) ;

  Unload();

  //  gime->PhosLoader()->GetSDigitsDataLoader()->GetBaseTaskLoader()->SetDoNotReload(kTRUE);
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSSDigitizer");
    Info("Exec","   took %f seconds for SDigitizing  %f seconds per event",
	 gBenchmark->GetCpuTime("PHOSSDigitizer"), gBenchmark->GetCpuTime("PHOSSDigitizer")/nEvents) ;
  }
  
  //TFile f("out.root","RECREATE");
  //memwatcher.WriteToFile(); 
  //f.Close();
}

//__________________________________________________________________
void AliPHOSSDigitizer::Print(const Option_t *)const
{
  // Prints parameters of SDigitizer
  Info("Print", "\n------------------- %s -------------", GetName() ) ; 
  printf("   Writing SDigits to branch with title  %s\n", fEventFolderName.Data()) ;
  printf("   with digitization parameters  A = %f\n", fA) ; 
  printf("                                 B = %f\n", fB) ;
  printf("   Threshold for Primary assignment= %f\n", fPrimThreshold)  ; 
  printf("---------------------------------------------------\n") ;
  
}

//__________________________________________________________________
Bool_t AliPHOSSDigitizer::operator==( AliPHOSSDigitizer const &sd )const
{
  // Equal operator.
  // SDititizers are equal if their pedestal, slope and threshold are equal

  if( (fA==sd.fA)&&(fB==sd.fB)&&(fPrimThreshold==sd.fPrimThreshold))
    return kTRUE ;
  else
    return kFALSE ;
}

//__________________________________________________________________
void AliPHOSSDigitizer::PrintSDigits(Option_t * option)
{
  // Prints list of digits produced in the current pass of AliPHOSDigitizer


  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  const TClonesArray * sdigits = gime->SDigits() ;
  
  Info( "\nPrintSDigits", "event # %d %d sdigits", gAlice->GetEvNumber(), sdigits->GetEntriesFast() ) ; 

  if(strstr(option,"all")||strstr(option,"EMC")){
    
    //loop over digits
    AliPHOSDigit * digit;
    printf("\nEMC sdigits\n") ; 
    Int_t maxEmc = gime->PHOSGeometry()->GetNModules()*gime->PHOSGeometry()->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; (index < sdigits->GetEntriesFast()) && 
	 ((dynamic_cast<AliPHOSDigit *> (sdigits->At(index)))->GetId() <= maxEmc) ; index++) {
      digit = dynamic_cast<AliPHOSDigit *>( sdigits->At(index) ) ;
      //  if(digit->GetNprimary() == 0) 
      // 	continue;
//       printf("%6d  %8d    %6.5e %4d      %2d :\n", // YVK
      printf("%6d  %.4f    %6.5e %4d      %2d :\n",
  	      digit->GetId(), digit->GetEnergy(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
	printf("%d ",digit->GetPrimary(iprimary+1) ) ; 
      }  
    }    
  }

  if(strstr(option,"all")||strstr(option,"CPV")){
    
    //loop over CPV digits
    AliPHOSDigit * digit;
    printf("\nCPV sdigits\n") ; 
    Int_t maxEmc = gime->PHOSGeometry()->GetNModules()*gime->PHOSGeometry()->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; index < sdigits->GetEntriesFast(); index++) {
      digit = dynamic_cast<AliPHOSDigit *>( sdigits->At(index) ) ;
      if(digit->GetId() > maxEmc){
	printf("\n%6d  %8d    %4d      %2d :",
		digit->GetId(), digit->GetAmp(), digit->GetIndexInList(), digit->GetNprimary()) ;  
	Int_t iprimary;
	for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
	  printf("%d ",digit->GetPrimary(iprimary+1) ) ; 
	}
      }    
    }
  }
}

//____________________________________________________________________________ 
void AliPHOSSDigitizer::Unload() const
{
  // Unloads the objects from the folder
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  AliPHOSLoader * loader = gime->PhosLoader() ; 
  loader->UnloadHits() ; 
  loader->UnloadSDigits() ; 
}
