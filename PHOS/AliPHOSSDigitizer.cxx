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
// A Summable Digits is the sum of all hits originating 
// from one primary in one active cell
// A threshold for assignment of the primary to SDigit is applied 
// SDigits are written to TreeS, branch "PHOS"
// AliPHOSSDigitizer with all current parameters is written 
// to TreeS branch "AliPHOSSDigitizer".
// Both branches, "PHOS" and "AliPHOSSDigitizer", are written to the same
// file, and therefore, changing branch file name one can produce several
// versions of SDigitization from the same hits.
// 
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
#include "AliPHOSHit.h"
#include "AliPHOSSDigitizer.h"


ClassImp(AliPHOSSDigitizer)

           
//____________________________________________________________________________ 
  AliPHOSSDigitizer::AliPHOSSDigitizer():TTask("AliPHOSSDigitizer","") 
{
  // ctor
  fA = 0;
  fB = 10000000. ;
  fPrimThreshold = 0.01 ;
  fNevents = 0 ;     
  fSDigits = 0 ;
  fHits = 0 ;
  fIsInitialized = kFALSE ;

}

//____________________________________________________________________________ 
AliPHOSSDigitizer::AliPHOSSDigitizer(const char* headerFile, const char *sDigitsTitle):TTask("AliPHOSSDigitizer","")
{
  // ctor
  fA = 0;
  fB = 10000000.;
  fPrimThreshold = 0.01 ;
  fNevents = 0 ;      
  fSDigitsTitle = sDigitsTitle ;
  fHeadersFile = headerFile ;
  fSDigits = new TClonesArray("AliPHOSDigit",1000);
  fHits    = new TClonesArray("AliPHOSHit",1000);

  TFile * file = (TFile*) gROOT->GetFile(fHeadersFile.Data() ) ;
  
  //File was not opened yet
  if(file == 0){
    file = new TFile(fHeadersFile.Data(),"update") ;
    gAlice = (AliRun *) file->Get("gAlice") ;
  }
  
  //add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
  
  fIsInitialized = kTRUE ;
}

//____________________________________________________________________________ 
AliPHOSSDigitizer::~AliPHOSSDigitizer()
{
  // dtor
  if(fSDigits)
    delete fSDigits ;
  if(fHits)
    delete fHits ;
}
//____________________________________________________________________________ 
void AliPHOSSDigitizer::Init(){
  //Initialization can not be done in the default constructor

  if(!fIsInitialized){

    if(fHeadersFile.IsNull())
      fHeadersFile="galice.root" ;

    TFile * file = (TFile*) gROOT->GetFile(fHeadersFile.Data() ) ;
    
    //if file was not opened yet, read gAlice
    if(file == 0){
      file = new TFile(fHeadersFile.Data(),"update") ;
      gAlice = (AliRun *) file->Get("gAlice") ;
    }
    
    fHits    = new TClonesArray("AliPHOSHit",1000);
    fSDigits = new TClonesArray("AliPHOSDigit",1000);
    
    // add Task to //root/Tasks folder
    TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
    roottasks->Add(this) ; 
    
    fIsInitialized = kTRUE ;
  }
}
//____________________________________________________________________________
void AliPHOSSDigitizer::Exec(Option_t *option) { 
  //Collects all hits in the same active volume into digit
  
  if(!fIsInitialized)
    Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSSDigitizer");
  
  fNevents = (Int_t) gAlice->TreeE()->GetEntries() ; 
  
  Int_t ievent ;
  for(ievent = 0; ievent < fNevents; ievent++){
    gAlice->GetEvent(ievent) ;
    gAlice->SetEvent(ievent) ;
    
    if(gAlice->TreeH()==0){
      cout << "AliPHOSSDigitizer: There is no Hit Tree" << endl;
      return ;
    }
    
    //set address of the hits 
    TBranch * branch = gAlice->TreeH()->GetBranch("PHOS");
    if (branch) 
      branch->SetAddress(&fHits);
    else{
      cout << "ERROR in AliPHOSSDigitizer: "<< endl ;
      cout << "      no branch PHOS in TreeH"<< endl ;
      cout << "      do nothing " << endl ;
      return ;
    }
    
    fSDigits->Clear();
    Int_t nSdigits = 0 ;
    
    
    //Now made SDigits from hits, for PHOS it is the same, so just copy    
    Int_t itrack ;
    for (itrack=0; itrack<gAlice->GetNtrack(); itrack++){
      
      //=========== Get the Hits Tree for the Primary track itrack
      gAlice->ResetHits();    
      gAlice->TreeH()->GetEvent(itrack);
      
      Int_t i;
      for ( i = 0 ; i < fHits->GetEntries() ; i++ ) {
	AliPHOSHit * hit = (AliPHOSHit*)fHits->At(i) ;
	AliPHOSDigit * newdigit ;

	// Assign primary number only if contribution is significant
	if( hit->GetEnergy() > fPrimThreshold)
	  newdigit = new AliPHOSDigit( hit->GetPrimary(), hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
	else
	  newdigit = new AliPHOSDigit( -1               , hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
	
	new((*fSDigits)[nSdigits]) AliPHOSDigit(* newdigit) ;
	nSdigits++ ;  
	
	delete newdigit ;    
      } 
      
    } // loop over tracks
    
    fSDigits->Sort() ;
    
    nSdigits = fSDigits->GetEntriesFast() ;
    fSDigits->Expand(nSdigits) ;
    
    Int_t i ;
    for (i = 0 ; i < nSdigits ; i++) { 
      AliPHOSDigit * digit = (AliPHOSDigit *) fSDigits->At(i) ; 
      digit->SetIndexInList(i) ;     
    }

    if(gAlice->TreeS() == 0)
      gAlice->MakeTree("S") ;
    
    //check, if this branch already exits?
    TBranch * sdigitsBranch = 0;
    TBranch * sdigitizerBranch = 0;
    
    TObjArray * branches = gAlice->TreeS()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t phosNotFound = kTRUE ;
    Bool_t sdigitizerNotFound = kTRUE ;
    
    for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
      
      if(phosNotFound){
	sdigitsBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp("PHOS",sdigitsBranch->GetName())==0 ) &&
	    (fSDigitsTitle.CompareTo(sdigitsBranch->GetTitle()) == 0) )
	  phosNotFound = kFALSE ;
      }
      if(sdigitizerNotFound){
	sdigitizerBranch = (TBranch *) branches->At(ibranch) ;
	if( (strcmp(sdigitizerBranch->GetName(),"AliPHOSSDigitizer") == 0)&&
	    (fSDigitsTitle.CompareTo(sdigitizerBranch->GetTitle()) == 0) )
	  sdigitizerNotFound = kFALSE ;
      }
    }

    if(!(sdigitizerNotFound && phosNotFound)){
      cout << "AliPHOSSdigitizer error:" << endl ;
      cout << "Can not overwrite existing branches: do not write" << endl ;
      return ;
    }
    
    //Make (if necessary) branches    
    char * file =0;
    if(gSystem->Getenv("CONFIG_SPLIT_FILE")){ //generating file name
      file = new char[strlen(gAlice->GetBaseFile())+20] ;
      sprintf(file,"%s/PHOS.SDigits.root",gAlice->GetBaseFile()) ;
    }
    
    TDirectory *cwd = gDirectory;
    
    //First list of sdigits
    Int_t bufferSize = 32000 ;    
    sdigitsBranch = gAlice->TreeS()->Branch("PHOS",&fSDigits,bufferSize);
    sdigitsBranch->SetTitle(fSDigitsTitle.Data());
    if (file) {
      sdigitsBranch->SetFile(file);
      TIter next( sdigitsBranch->GetListOfBranches());
      while ((sdigitsBranch=(TBranch*)next())) {
	sdigitsBranch->SetFile(file);
      }   
      cwd->cd();
    } 
      
    //second - SDigitizer
    Int_t splitlevel = 0 ;
    AliPHOSSDigitizer * sd = this ;
    sdigitizerBranch = gAlice->TreeS()->Branch("AliPHOSSDigitizer","AliPHOSSDigitizer",
					       &sd,bufferSize,splitlevel); 
    sdigitizerBranch->SetTitle(fSDigitsTitle.Data());
    if (file) {
      sdigitizerBranch->SetFile(file);
      TIter next( sdigitizerBranch->GetListOfBranches());
      while ((sdigitizerBranch=(TBranch*)next())) {
	sdigitizerBranch->SetFile(file);
      }   
      cwd->cd();
      delete file;
    }

    gAlice->TreeS()->Fill() ;
    gAlice->TreeS()->Write(0,TObject::kOverwrite) ;
    
    if(strstr(option,"deb"))
      PrintSDigits(option) ;
    
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSSDigitizer");
    cout << "AliPHOSSDigitizer:" << endl ;
    cout << "   took " << gBenchmark->GetCpuTime("PHOSSDigitizer") << " seconds for SDigitizing " 
	 <<  gBenchmark->GetCpuTime("PHOSSDigitizer")/fNevents << " seconds per event " << endl ;
    cout << endl ;
  }
  
  
}
//__________________________________________________________________
void AliPHOSSDigitizer::SetSDigitsBranch(const char * title ){
  //Seting title to branch SDigits 
  if(!fSDigitsTitle.IsNull())
    cout << "AliPHOSSdigitizer: changing SDigits file from " <<fSDigitsTitle.Data() << " to " << title << endl ;
  fSDigitsTitle=title ;
}
//__________________________________________________________________
void AliPHOSSDigitizer::Print(Option_t* option)const{
  cout << "------------------- "<< GetName() << " -------------" << endl ;
  cout << "   Writing SDigitis to branch with title  " << fSDigitsTitle.Data() << endl ;
  cout << "   with digitization parameters  A = " << fA << endl ;
  cout << "                                 B = " << fB << endl ;
  cout << "   Threshold for Primary assignment= " << fPrimThreshold << endl ; 
  cout << "---------------------------------------------------"<<endl ;
  
}
//__________________________________________________________________
Bool_t AliPHOSSDigitizer::operator==( AliPHOSSDigitizer const &sd )const{
  if( (fA==sd.fA)&&(fB==sd.fB)&&(fPrimThreshold==sd.fPrimThreshold))
    return kTRUE ;
  else
    return kFALSE ;
}
//__________________________________________________________________
void AliPHOSSDigitizer::PrintSDigits(Option_t * option){
  //Prints list of digits produced at the current pass of AliPHOSDigitizer
  
  cout << "AliPHOSSDigitizer: " << endl ;
  cout << "       Number of entries in SDigits list  " << fSDigits->GetEntriesFast() << endl ;
  cout << endl ;
  
  if(strstr(option,"all")){// print all digits
    
    //loop over digits
    AliPHOSDigit * digit;
    cout << "SDigit Id " << " Amplitude " <<  " Index "  <<  " Nprim " << " Primaries list " <<  endl;    
    Int_t index ;
    for (index = 0 ; index < fSDigits->GetEntries() ; index++) {
      digit = (AliPHOSDigit * )  fSDigits->At(index) ;
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
