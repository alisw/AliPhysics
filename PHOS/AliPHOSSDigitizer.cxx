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
#include "TFile.h"
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TBenchmark.h"
#include "TGeometry.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliHeader.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSGetter.h"
#include "AliPHOSHit.h"
#include "AliPHOSSDigitizer.h"

ClassImp(AliPHOSSDigitizer)

           
//____________________________________________________________________________ 
  AliPHOSSDigitizer::AliPHOSSDigitizer():TTask("","") 
{
  // ctor
  InitParameters() ;
  fDefaultInit = kTRUE ; 
}

//____________________________________________________________________________ 
AliPHOSSDigitizer::AliPHOSSDigitizer(const char * headerFile, const char * sDigitsTitle, const Bool_t toSplit):
TTask(sDigitsTitle, headerFile)
{
  // ctor
  InitParameters() ; 
  fToSplit = toSplit ;
  Init();
  fDefaultInit = kFALSE ; 
}

//____________________________________________________________________________ 
AliPHOSSDigitizer::~AliPHOSSDigitizer()
{
  // dtor
  
  fSplitFile = 0 ; 
}

//____________________________________________________________________________ 
void AliPHOSSDigitizer::Init()
{
  // Initialization: open root-file, allocate arrays for hits and sdigits,
  // attach task SDigitizer to the list of PHOS tasks
  // 
  // Initialization can not be done in the default constructor
  //============================================================= YS
  //  The initialisation is now done by AliPHOSGetter
  
  if( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(GetTitle(), GetName(),fToSplit) ;  
  if ( gime == 0 ) {
    Error("Init" ,"Could not obtain the Getter object !") ;  
    return ;
  } 
  
  gime->PostSDigits( GetName(), GetTitle() ) ; 

  fSplitFile = 0 ;
  if(fToSplit){
    // construct the name of the file as /path/PHOS.SDigits.root
    // First - extract full path if necessary
    TString sDigitsFileName(GetTitle()) ;
    Ssiz_t islash = sDigitsFileName.Last('/') ;
    if(islash<sDigitsFileName.Length())
      sDigitsFileName.Remove(islash+1,sDigitsFileName.Length()) ;
    else
      sDigitsFileName="" ;
    // Next - append the file name 
    sDigitsFileName+="PHOS.SDigits." ;
    if((strcmp(GetName(),"Default")!=0)&&(strcmp(GetName(),"")!=0)){
      sDigitsFileName+=GetName() ;
      sDigitsFileName+="." ;
    }
    sDigitsFileName+="root" ;
    // Finally - check if the file already opened or open the file
    fSplitFile = static_cast<TFile*>(gROOT->GetFile(sDigitsFileName.Data()));   
    if(!fSplitFile)
      fSplitFile =  TFile::Open(sDigitsFileName.Data(),"update") ;
  }

  TString sdname(GetName() );
  sdname.Append(":") ;
  sdname.Append(GetTitle() ) ;
  SetName(sdname) ;
  gime->PostSDigitizer(this) ;
}

//____________________________________________________________________________ 
void AliPHOSSDigitizer::InitParameters()
{ 
  fA             = 0;
  fB             = 10000000.;
  fPrimThreshold = 0.01 ;
  fSDigitsInRun  = 0 ;
  fSplitFile     = 0 ; 
  fToSplit       = kFALSE ;
}

//____________________________________________________________________________
void AliPHOSSDigitizer::Exec(Option_t *option) 
{ 
  // Collects all hits in the same active volume into digit

  if( strcmp(GetName(), "") == 0 )
    Init() ;
  
  if (strstr(option, "print") ) {
    Print("") ; 
    return ; 
  }

  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSSDigitizer");
  
  //Check, if this branch already exits
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  if(gime->BranchExists("SDigits") ) 
    return;   

  TString sdname(GetName()) ;
  sdname.Remove(sdname.Index(GetTitle())-1) ;
   
  Int_t nevents = gime->MaxEvent() ; 
  Int_t ievent ;
  for(ievent = 0; ievent < nevents; ievent++){
    gime->Event(ievent,"H") ;
    const TClonesArray * hits = gime->Hits() ;
    TClonesArray * sdigits = gime->SDigits(sdname.Data()) ;
    sdigits->Clear();
    Int_t nSdigits = 0 ;
    
    //Now make SDigits from hits, for PHOS it is the same, so just copy    
    Int_t nPrim =  static_cast<Int_t>((gAlice->TreeH())->GetEntries()) ; 
    // Attention nPrim is the number of primaries tracked by Geant 
    // and this number could be different to the number of Primaries in TreeK;
    Int_t iprim ;
    for (iprim = 0 ; iprim < nPrim ; iprim ++) { 
      //=========== Get the PHOS branch from Hits Tree for the Primary iprim
      gime->Track(iprim) ;
      Int_t i;
      for ( i = 0 ; i < hits->GetEntries() ; i++ ) {
	AliPHOSHit * hit = dynamic_cast<AliPHOSHit *>(hits->At(i)) ;
	// Assign primary number only if contribution is significant
	
	if( hit->GetEnergy() > fPrimThreshold)
	  new((*sdigits)[nSdigits]) AliPHOSDigit(hit->GetPrimary(),hit->GetId(),
						 Digitize(hit->GetEnergy()), hit->GetTime()) ;
	else
	  new((*sdigits)[nSdigits]) AliPHOSDigit( -1              , hit->GetId(), 
						  Digitize(hit->GetEnergy()), hit->GetTime()) ;
	nSdigits++ ;	
	
      }
      
    } // loop over iprim
    
    sdigits->Sort() ;
    
    nSdigits = sdigits->GetEntriesFast() ;
    fSDigitsInRun += nSdigits ;  
    sdigits->Expand(nSdigits) ;
    
    Int_t i ;
    for (i = 0 ; i < nSdigits ; i++) { 
      AliPHOSDigit * digit = dynamic_cast<AliPHOSDigit *>(sdigits->At(i)) ; 
      digit->SetIndexInList(i) ;     
    }

    //Now write SDigits
    
    if((gAlice->TreeS() == 0)|| (fSplitFile))  
      gAlice->MakeTree("S", fSplitFile);
    
    if(fSplitFile)
      fSplitFile->cd() ;

    //First list of sdigits
    Int_t bufferSize = 32000 ;
    TBranch * sdigitsBranch = gAlice->TreeS()->Branch("PHOS",&sdigits,bufferSize);
    sdigitsBranch->SetTitle(sdname);
    
    //Next - SDigitizer
    Int_t splitlevel = 0 ;
    AliPHOSSDigitizer * sd = this ;
    TBranch * sdigitizerBranch = gAlice->TreeS()->Branch("AliPHOSSDigitizer","AliPHOSSDigitizer",
							 &sd,bufferSize,splitlevel);
    sdigitizerBranch->SetTitle(sdname);
    
    sdigitsBranch->Fill() ;
    sdigitizerBranch->Fill() ;

    gAlice->TreeS()->AutoSave() ;
        
    if(strstr(option,"deb"))
      PrintSDigits(option) ;
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSSDigitizer");
    Info("Exec","   took %f seconds for SDigitizing  %f seconds per event",
	 gBenchmark->GetCpuTime("PHOSSDigitizer"), gBenchmark->GetCpuTime("PHOSSDigitizer")/nevents) ;
  }
}

//__________________________________________________________________
void AliPHOSSDigitizer::SetSDigitsBranch(const char * title )
{
  // Setting title to branch SDigits 

  TString stitle(title) ;

  // check if branch with title already exists
  TBranch * sdigitsBranch    = 
    static_cast<TBranch*>(gAlice->TreeS()->GetListOfBranches()->FindObject("PHOS")) ; 
  TBranch * sdigitizerBranch =  
    static_cast<TBranch*>(gAlice->TreeS()->GetListOfBranches()->FindObject("AliPHOSSDigitizer")) ;
  const char * sdigitsTitle    = sdigitsBranch ->GetTitle() ;  
  const char * sdigitizerTitle = sdigitizerBranch ->GetTitle() ;
  if ( stitle.CompareTo(sdigitsTitle)==0 || stitle.CompareTo(sdigitizerTitle)==0 ){
    Error("SetSDigitsBranch", "Cannot overwrite existing branch with title %s", title) ;
    return ;
  }
  
  Info("SetSDigitsBranch", "-> Changing SDigits file from %s to %s", GetName(), title) ;

  SetName(title) ; 
    
  // Post to the WhiteBoard
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  gime->PostSDigits( title, GetTitle()) ; 
}


//__________________________________________________________________
void AliPHOSSDigitizer::Print(Option_t* option)const
{
  // Prints parameters of SDigitizer
  TString message ; 
  message  = "\n------------------- %s -------------\n" ;  
  message += "   Writing SDigits to branch with title  %s\n" ;
  message += "   with digitization parameters  A = %f\n" ; 
  message += "                                 B = %f\n" ;
  message += "   Threshold for Primary assignment= %f\n" ; 
  message += "---------------------------------------------------\n" ;
  Info("Print", message.Data(),  GetName(),  GetName(), fA, fB, fPrimThreshold ) ;
  
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


  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  TString sdname(GetName()) ;
  sdname.Remove(sdname.Index(GetTitle())-1) ;
  const TClonesArray * sdigits = gime->SDigits(sdname.Data()) ;

  TString message ; 
  message  = "\nAliPHOSSDigitiser: event " ;
  message += gAlice->GetEvNumber(); 
  message += "\n      Number of entries in SDigits list " ;  
  message += sdigits->GetEntriesFast() ; 
  
  if(strstr(option,"all")||strstr(option,"EMC")){
    
    //loop over digits
    AliPHOSDigit * digit;
    message += "\nEMC sdigits\n" ;
    message += "Digit Id    Amplitude     Index     Nprim  Primaries list\n" ;       
    Int_t maxEmc = gime->PHOSGeometry()->GetNModules()*gime->PHOSGeometry()->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; (index < sdigits->GetEntriesFast()) && 
	 ((dynamic_cast<AliPHOSDigit *> (sdigits->At(index)))->GetId() <= maxEmc) ; index++) {
      digit = dynamic_cast<AliPHOSDigit *>( sdigits->At(index) ) ;
      if(digit->GetNprimary() == 0) 
	continue;
      message += digit->GetId() ; 
      message += digit->GetAmp() ;
      message += digit->GetIndexInList() ;
      message += digit->GetNprimary() ;
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++)
	message += digit->GetPrimary(iprimary+1) ;
    }    
  }

  if(strstr(option,"all")||strstr(option,"CPV")){
    
    //loop over CPV digits
    AliPHOSDigit * digit;
    
    message += "CPV sdigits\n" ;
    message += "Digit Id Amplitude Index Nprim  Primaries list\n" ;
    Int_t maxEmc = gime->PHOSGeometry()->GetNModules()*gime->PHOSGeometry()->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; index < sdigits->GetEntriesFast(); index++) {
      digit = dynamic_cast<AliPHOSDigit *>( sdigits->At(index) ) ;
      if(digit->GetId() > maxEmc){
	message += "\n" ; 
	message += digit->GetId() ; 
	message += "   " ; 
	message += digit->GetAmp() ;
	message += "   " ; 
	message += digit->GetIndexInList() ;
	message += "   " ; 
	message += digit->GetNprimary() ;
	message += "   " ; 
	Int_t iprimary;
	for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
	  message += digit->GetPrimary(iprimary+1) ;
	  message += "," ; 
	}
      }    
    }
  }
  Info("PrintSDigits", message.Data() ) ;
}

//____________________________________________________________________________ 
void AliPHOSSDigitizer::UseHitsFrom(const char * filename)
{
  SetTitle(filename) ; 
  Init() ; 
}
