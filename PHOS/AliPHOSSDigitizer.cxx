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
// --- Standard library ---
#include <iomanip.h>

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGetter.h"
#include "AliPHOSHit.h"
#include "AliPHOSSDigitizer.h"


ClassImp(AliPHOSSDigitizer)

           
//____________________________________________________________________________ 
  AliPHOSSDigitizer::AliPHOSSDigitizer():TTask("","") 
{
  // ctor
  fA             = 0;
  fB             = 10000000.;
  fPrimThreshold = 0.01 ;
  fSDigitsInRun  = 0 ; 
}

//____________________________________________________________________________ 
AliPHOSSDigitizer::AliPHOSSDigitizer(const char * headerFile, const char * sDigitsTitle):TTask(sDigitsTitle, headerFile)
{
  // ctor
  fA             = 0;
  fB             = 10000000.;
  fPrimThreshold = 0.01 ;
  fSDigitsInRun  = 0 ; 
  Init();
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
  
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(GetTitle(), GetName()) ;     
  if ( gime == 0 ) {
    cerr << "ERROR: AliPHOSSDigitizer::Init -> Could not obtain the Getter object !" << endl ; 
    return ;
  } 
  
  gime->PostSDigits( GetName(), GetTitle() ) ; 

  TString sdname(GetName() );
  sdname.Append(":") ;
  sdname.Append(GetTitle() ) ;
  SetName(sdname) ;
  gime->PostSDigitizer(this) ;
 
  // create a folder on the white board //YSAlice/WhiteBoard/SDigits/PHOS/headerFile/sdigitsTitle
    
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
  gAlice->GetEvent(0) ;
  if(gAlice->TreeS() ) {
    TObjArray * lob = static_cast<TObjArray*>(gAlice->TreeS()->GetListOfBranches()) ;
    TIter next(lob) ; 
    TBranch * branch = 0 ;  
    Bool_t phosfound = kFALSE, sdigitizerfound = kFALSE ; 
    
    while ( (branch = (static_cast<TBranch*>(next()))) && (!phosfound || !sdigitizerfound) ) {
      if ( (strcmp(branch->GetName(), "PHOS")==0) && (strcmp(branch->GetTitle(), GetName())==0) ) 
	phosfound = kTRUE ;
      
      else if ( (strcmp(branch->GetName(), "AliPHOSSDigitizer")==0) && (strcmp(branch->GetTitle(), GetName())==0) ) 
	sdigitizerfound = kTRUE ; 
    }
    
    if ( phosfound || sdigitizerfound ) {
      cerr << "WARNING: AliPHOSSDigitizer::Exec -> SDigits and/or SDigitizer branch with name " << GetName() 
	   << " already exits" << endl ;
      return ; 
    }   
  }  

  TString sdname(GetName()) ;
  sdname.Remove(sdname.Index(GetTitle())-1) ;
    
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 

  Int_t nevents = (Int_t) gAlice->TreeE()->GetEntries() ; 
  
  Int_t ievent ;
  for(ievent = 0; ievent < nevents; ievent++){
    gime->Event(ievent,"H") ;

    TClonesArray * sdigits = gime->SDigits(sdname.Data()) ;
    sdigits->Clear();
    Int_t nSdigits = 0 ;
    
    
    //Now make SDigits from hits, for PHOS it is the same, so just copy    
    Int_t itrack ;
    for (itrack=0; itrack < gAlice->GetNtrack(); itrack++){
  
      //=========== Get the PHOS branch from Hits Tree for the Primary track itrack
      gime->Track(itrack) ;
      TClonesArray * hits = gime->Hits() ;
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
    } // loop over tracks
    
    sdigits->Sort() ;
    
    nSdigits = sdigits->GetEntriesFast() ;
    
    sdigits->Expand(nSdigits) ;
    Int_t i ;
    for (i = 0 ; i < nSdigits ; i++) { 
      AliPHOSDigit * digit = dynamic_cast<AliPHOSDigit *>(sdigits->At(i)) ; 
      digit->SetIndexInList(i) ;     
    }

    if(gAlice->TreeS() == 0)
      gAlice->MakeTree("S") ;
    
    //Make (if necessary) branches    
    char * file =0;
    if(gSystem->Getenv("CONFIG_SPLIT_FILE")){ //generating file name
      file = new char[strlen(gAlice->GetBaseFile())+20] ;
      sprintf(file,"%s/PHOS.SDigits.root",gAlice->GetBaseFile()) ;
    }
    
    TDirectory *cwd = gDirectory;
    
    //First list of sdigits
    Int_t bufferSize = 32000 ;    
    TBranch * sdigitsBranch = gAlice->TreeS()->Branch("PHOS",&sdigits,bufferSize);
    sdigitsBranch->SetTitle(sdname);
    if (file) {
      sdigitsBranch->SetFile(file);
      TIter next( sdigitsBranch->GetListOfBranches());
      TBranch * subbr;
      while ((subbr=static_cast<TBranch*>(next()))) {
	subbr->SetFile(file);
      }   
      cwd->cd();
    } 
      
    //Next - SDigitizer
    Int_t splitlevel = 0 ;
    AliPHOSSDigitizer * sd = this ;
    TBranch * sdigitizerBranch = gAlice->TreeS()->Branch("AliPHOSSDigitizer","AliPHOSSDigitizer",
					       &sd,bufferSize,splitlevel); 
    sdigitizerBranch->SetTitle(sdname);
    if (file) {
      sdigitizerBranch->SetFile(file);
      TIter next( sdigitizerBranch->GetListOfBranches());
      TBranch * subbr ;
      while ((subbr=static_cast<TBranch*>(next()))) {
	subbr->SetFile(file);
      }   
      cwd->cd();
      delete file;
    }

    sdigitsBranch->Fill() ;
    sdigitizerBranch->Fill() ;
    gAlice->TreeS()->Write(0,TObject::kOverwrite) ;
    
    if(strstr(option,"deb"))
      PrintSDigits(option) ;
    
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSSDigitizer");
    cout << "AliPHOSSDigitizer:" << endl ;
    cout << "   took " << gBenchmark->GetCpuTime("PHOSSDigitizer") << " seconds for SDigitizing " 
	 <<  gBenchmark->GetCpuTime("PHOSSDigitizer")/nevents << " seconds per event " << endl ;
    cout << endl ;
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
    cerr << "ERROR: AliPHOSSdigitizer::SetSDigitsBranch -> Cannot overwrite existing branch with title " << title << endl ;
    return ;
  }
  
  cout << "AliPHOSSdigitizer::SetSDigitsBranch -> Changing SDigits file from " << GetName() << " to " << title << endl ;

  SetName(title) ; 
    
  // Post to the WhiteBoard
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  gime->PostSDigits( title, GetTitle()) ; 
}

//__________________________________________________________________
void AliPHOSSDigitizer::Print(Option_t* option)const
{
  // Prints parameters of SDigitizer
  cout << "------------------- "<< GetName() << " -------------" << endl ;
  cout << "   Writing SDigits to branch with title  " << GetName() << endl ;
  cout << "   with digitization parameters  A = " << fA << endl ;
  cout << "                                 B = " << fB << endl ;
  cout << "   Threshold for Primary assignment= " << fPrimThreshold << endl ; 
  cout << "---------------------------------------------------"<<endl ;
  
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
  TClonesArray * sdigits = gime->SDigits(sdname.Data()) ; 

  cout << "AliPHOSSDigitizer: event "<<gAlice->GetEvNumber() << endl ;
  cout << "       Number of entries in SDigits list  " << sdigits->GetEntriesFast() << endl ;
  cout << endl ;

  fSDigitsInRun +=  sdigits->GetEntriesFast() ; 
  
  if(strstr(option,"all")){// print all digits
    
    //loop over digits
    AliPHOSDigit * digit;
    cout << "SDigit Id " << " Amplitude " <<  " Index "  <<  " Nprim " << " Primaries list " <<  endl;    
    Int_t index ;
    for (index = 0 ; index < sdigits->GetEntries() ; index++) {
      digit = (AliPHOSDigit * )  sdigits->At(index) ;
      cout << setw(8)  <<  digit->GetId() << " "  << 	setw(3)  <<  digit->GetAmp() <<   "  "  
	   << setw(6)  <<  digit->GetIndexInList() << "  "   
	   << setw(5)  <<  digit->GetNprimary() <<"  ";
      
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++)
	cout << setw(5)  <<  digit->GetPrimary(iprimary+1) << "  ";
      cout << endl;  	 
    }
    
  }
}

//____________________________________________________________________________ 
void AliPHOSSDigitizer::UseHitsFrom(const char * filename)
{
  SetTitle(filename) ; 
  Init() ; 
}
