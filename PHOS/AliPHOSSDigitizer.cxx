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
 * Revision 1.56  2007/10/19 18:04:29  schutz
 * The standalone QA data maker is called from AliSimulation and AliReconstruction outside the event loop; i.e. re-reading the data. The QA data making in the event loop has been commented out.
 *
 * Revision 1.55  2007/10/14 21:08:10  schutz
 * Introduced the checking of QA results from previous step before entering the event loop
 *
 * Revision 1.54  2007/10/10 09:05:10  schutz
 * Changing name QualAss to QA
 *
 * Revision 1.53  2007/09/30 17:08:20  schutz
 * Introducing the notion of QA data acquisition cycle (needed by online)
 *
 * Revision 1.52  2007/09/26 14:22:18  cvetan
 * Important changes to the reconstructor classes. Complete elimination of the run-loaders, which are now steered only from AliReconstruction. Removal of the corresponding Reconstruct() and FillESD() methods.
 *
 * Revision 1.51  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
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
//  root [1] s->Digitize()
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
//-- Author :  Dmitri Peressounko (SUBATECH & KI) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TBenchmark.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSGeometry.h" 
#include "AliPHOSDigit.h"
#include "AliPHOSLoader.h"
#include "AliPHOSHit.h"
#include "AliPHOSSDigitizer.h"

ClassImp(AliPHOSSDigitizer)

           
//____________________________________________________________________________ 
AliPHOSSDigitizer::AliPHOSSDigitizer() : 
  TNamed("",""),
  fPrimThreshold(0.f),
  fDefaultInit(kTRUE),
  fEventFolderName(""),
  fInit(kFALSE),
  fSDigitsInRun(0),
  fFirstEvent(0),
  fLastEvent(0)
{
  // ctor
  // Intialize the quality assurance data maker 	
}

//____________________________________________________________________________ 
AliPHOSSDigitizer::AliPHOSSDigitizer(const char * alirunFileName, 
				     const char * eventFolderName):
  TNamed("PHOSSDigitizer", alirunFileName),
  fPrimThreshold(0.f),
  fDefaultInit(kFALSE),
  fEventFolderName(eventFolderName),
  fInit(kFALSE),
  fSDigitsInRun(0),
  fFirstEvent(0),
  fLastEvent(0)
{
  // ctor
  InitParameters() ; 
  Init();
  fDefaultInit = kFALSE ; 
}

//____________________________________________________________________________
AliPHOSSDigitizer::AliPHOSSDigitizer(const AliPHOSSDigitizer& sd) :
  TNamed(sd.GetName(), sd.GetTitle()),
  fPrimThreshold(sd.fPrimThreshold),
  fDefaultInit(kFALSE),
  fEventFolderName(sd.fEventFolderName),
  fInit(kFALSE),
  fSDigitsInRun(sd.fSDigitsInRun),
  fFirstEvent(sd.fFirstEvent),
  fLastEvent(sd.fLastEvent)
{ 
  // cpy ctor
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
AliPHOSSDigitizer::~AliPHOSSDigitizer() 
{
  //dtor
}

//____________________________________________________________________________ 
void AliPHOSSDigitizer::Init()
{
  // Uses the Loader to access the required files
  
  fInit = kTRUE ; 
  
  //to prevent cleaning of this object while GetEvent is called
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  if (!rl) rl = AliRunLoader::Open(GetTitle(), fEventFolderName) ; 
}

//____________________________________________________________________________ 
void AliPHOSSDigitizer::InitParameters()
{ 
  // initializes the parameters for digitization
  fPrimThreshold = 0.01 ;
  fSDigitsInRun  = 0 ;
}

//____________________________________________________________________________
void AliPHOSSDigitizer::Digitize(Option_t *option) 
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
  
/*
	// check the QA result for RAWS
  AliQAv1 * qa = AliQAv1::Instance(AliQAv1::kPHOS) ; 
  if ( qa->IsSet(AliQAv1::kPHOS, AliQAv1::kRAW, AliQAv1::kFATAL)) {
	AliFatal("QA status in RAW was Fatal") ;
  } else if ( qa->IsSet(AliQAv1::kPHOS, AliQAv1::kRAW, AliQAv1::kERROR)) {
	AliError("QA status in RAW was Error") ;
  } else if ( qa->IsSet(AliQAv1::kPHOS, AliQAv1::kRAW, AliQAv1::kWARNING) ) {
	AliWarning("QA status in RAW was Warning") ;
  } else if ( qa->IsSet(AliQAv1::kPHOS, AliQAv1::kRAW, AliQAv1::kINFO) ) {
	AliInfo("QA status in RAW was Info") ;
  }
*/
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  AliPHOSLoader * phosLoader = static_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));
  
  //switch off reloading of this task while getting event
  if (!fInit) { // to prevent overwrite existing file
    AliError( Form("Give a version name different from %s", fEventFolderName.Data()) ) ;
    return ;
  }
  
  if (fLastEvent == -1) 
    fLastEvent = rl->GetNumberOfEvents() - 1 ;
  else 
    fLastEvent = TMath::Min(fFirstEvent, rl->GetNumberOfEvents()); // only one event at the time 
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;
  
  Int_t ievent, i;

  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {    
    rl->GetEvent(ievent) ;
    TTree * treeS = phosLoader->TreeS(); 
    if(!treeS){
      phosLoader->MakeTree("S");
      treeS = phosLoader->TreeS();
    }

    phosLoader->CleanHits() ; 
    phosLoader->LoadHits("READ") ;

    TClonesArray * hits    = phosLoader->Hits() ;
    TClonesArray * sdigits = phosLoader->SDigits() ;
    if( !sdigits ) {
      phosLoader->MakeSDigitsArray() ; 
      sdigits = phosLoader->SDigits() ;
    }
    sdigits->Clear();
    Int_t nSdigits = 0 ;

    //Now make SDigits from hits, for PHOS it is the same, so just copy    
    for ( i = 0 ; i < hits->GetEntries() ; i++ ) {
      
      AliPHOSHit * hit = static_cast<AliPHOSHit *>(hits->At(i)) ;
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
      AliPHOSDigit * digit = static_cast<AliPHOSDigit *>(sdigits->At(i)) ;
      digit->SetIndexInList(i) ;     
    }
    
//    // make Quality Assurance data
//
//    if (GetQADataMaker()->IsCycleDone() ) {
//      GetQADataMaker()->EndOfCycle(AliQAv1::kHITS) ; 
//	  GetQADataMaker()->EndOfCycle(AliQAv1::kSDIGITS) ; 
//      GetQADataMaker()->StartOfCycle(AliQAv1::kHITS) ; 
//	  GetQADataMaker()->StartOfCycle(AliQAv1::kSDIGITS, kTRUE) ; 
//   }
//    GetQADataMaker()->Exec(AliQAv1::kHITS, hits) ; 
//    GetQADataMaker()->Exec(AliQAv1::kSDIGITS, sdigits) ; 
//    GetQADataMaker()->Increment() ;
	
    //Now write SDigits

    
    //First list of sdigits

    Int_t bufferSize = 32000 ;
    TBranch * sdigitsBranch = treeS->Branch("PHOS",&sdigits,bufferSize);
    sdigitsBranch->Fill() ;

    phosLoader->WriteSDigits("OVERWRITE");

    if(strstr(option,"deb"))
      PrintSDigits(option) ;
      
  }// event loop
  
//  //Write the quality assurance data 
//  GetQADataMaker()->EndOfCycle(AliQAv1::kHITS) ;    
//  GetQADataMaker()->EndOfCycle(AliQAv1::kSDIGITS) ;    
//  GetQADataMaker()->Finish() ;

  Unload();

  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSSDigitizer");
    Info("Exec","   took %f seconds for SDigitizing  %f seconds per event",
	 gBenchmark->GetCpuTime("PHOSSDigitizer"), 
	 gBenchmark->GetCpuTime("PHOSSDigitizer")/nEvents) ;
  }
  
}

//__________________________________________________________________
void AliPHOSSDigitizer::Print(const Option_t *)const
{
  // Prints parameters of SDigitizer
  Info("Print", "\n------------------- %s -------------", GetName() ) ; 
  printf("   Writing SDigits to branch with title  %s\n", fEventFolderName.Data()) ;
  printf("   Threshold for Primary assignment= %f\n", fPrimThreshold)  ; 
  printf("---------------------------------------------------\n") ;
  
}

//__________________________________________________________________
Bool_t AliPHOSSDigitizer::operator==( AliPHOSSDigitizer const &sd )const
{
  // Equal operator.
  // SDititizers are equal if their pedestal, slope and threshold are equal

  if(fPrimThreshold==sd.fPrimThreshold)
    return kTRUE ;
  else
    return kFALSE ;
}

//__________________________________________________________________
void AliPHOSSDigitizer::PrintSDigits(Option_t * option)
{
  // Prints list of digits produced in the current pass of AliPHOSDigitizer

  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  AliPHOSLoader * phosLoader = static_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));

  // Get PHOS Geometry object
  AliPHOSGeometry *geom;
  if (!(geom = AliPHOSGeometry::GetInstance())) 
        geom = AliPHOSGeometry::GetInstance("IHEP","");

  const TClonesArray * sdigits = phosLoader->SDigits() ;
  
  Info( "\nPrintSDigits", "event # %d %d sdigits", 
	gAlice->GetEvNumber(), sdigits->GetEntriesFast() ) ; 

  if(strstr(option,"all")||strstr(option,"EMC")){
    
    //loop over digits
    AliPHOSDigit * digit;
    printf("\nEMC sdigits\n") ; 
    Int_t maxEmc = geom->GetNModules() * geom->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; (index < sdigits->GetEntriesFast()) && 
	 ((static_cast<AliPHOSDigit *> (sdigits->At(index)))->GetId() <= maxEmc) ; index++) {
      digit = static_cast<AliPHOSDigit *>( sdigits->At(index) ) ;
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
    Int_t maxEmc = geom->GetNModules() * geom->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; index < sdigits->GetEntriesFast(); index++) {
      digit = static_cast<AliPHOSDigit *>( sdigits->At(index) ) ;
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
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEventFolderName) ;
  AliPHOSLoader * phosLoader = static_cast<AliPHOSLoader*>(rl->GetLoader("PHOSLoader"));
  
  phosLoader->UnloadHits() ; 
  phosLoader->UnloadSDigits() ; 
}
