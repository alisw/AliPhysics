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
#include "AliEMCALDigit.h"
#include "AliEMCALHit.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALv1.h"

ClassImp(AliEMCALSDigitizer)

           
//____________________________________________________________________________ 
  AliEMCALSDigitizer::AliEMCALSDigitizer():TTask("AliEMCALSDigitizer","") 
{
  // ctor
  fA = 0;
  fB = 10000000. ;
  fPrimThreshold = 0.001 ;
  fNevents = 0 ;     
  fSDigits = 0 ;
  fHits = 0 ;
  fIsInitialized = kFALSE ;

}

//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer(const char* headerFile, const char *sDigitsTitle):TTask("AliEMCALSDigitizer","")
{
  // ctor
  fA = 0;
  fB = 10000000.;
  fPrimThreshold = 0.001 ;
  fNevents = 0 ;      
  fSDigitsTitle = sDigitsTitle ;
  fHeadersFile = headerFile ;
  fSDigits = new TClonesArray("AliEMCALDigit",1000);
  fHits    = new TClonesArray("AliEMCALHit",1000);

  TFile * file = (TFile*) gROOT->GetFile(fHeadersFile.Data() ) ;
  
  //File was not opened yet
  if(file == 0){
    if(fHeadersFile.Contains("rfio"))
      file =	TFile::Open(fHeadersFile,"update") ;
    else
      file = new TFile(fHeadersFile.Data(),"update") ;
    gAlice = (AliRun *) file->Get("gAlice") ;
  }
  
  //add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
  
  fIsInitialized = kTRUE ;
}

//____________________________________________________________________________ 
AliEMCALSDigitizer::~AliEMCALSDigitizer()
{
  // dtor
  if(fSDigits)
    delete fSDigits ;
  if(fHits)
    delete fHits ;
}
//____________________________________________________________________________ 
void AliEMCALSDigitizer::Init(){
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
    
    fHits    = new TClonesArray("AliEMCALHit",1000);
    fSDigits = new TClonesArray("AliEMCALDigit",1000);
    
    // add Task to //root/Tasks folder
    TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
    roottasks->Add(this) ; 
    
    fIsInitialized = kTRUE ;
  }
}
//____________________________________________________________________________
void AliEMCALSDigitizer::Exec(Option_t *option) { 
  //Collects all hits in the same active volume into digit
  
  if(!fIsInitialized)
    Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALSDigitizer");
  
  fNevents = (Int_t) gAlice->TreeE()->GetEntries() ; 
  
  Int_t ievent ;
  for(ievent = 0; ievent < fNevents; ievent++){
    gAlice->GetEvent(ievent) ;
    gAlice->SetEvent(ievent) ;
    
    if(gAlice->TreeH()==0){
      cout << "AliEMCALSDigitizer: There is no Hit Tree" << endl;
      return ;
    }
    
    //set address of the hits 
    TBranch * branch = gAlice->TreeH()->GetBranch("EMCAL");
    if (branch) 
      branch->SetAddress(&fHits);
    else{
      cout << "ERROR in AliEMCALSDigitizer: "<< endl ;
      cout << "      no branch EMCAL in TreeH"<< endl ;
      cout << "      do nothing " << endl ;
      return ;
    }
    
    fSDigits->Clear();
    Int_t nSdigits = 0 ;
    
    
    //Now made SDigits from hits, for EMCAL it is the same, so just copy    
    Int_t itrack ;
    for (itrack=0; itrack < gAlice->GetNtrack(); itrack++){
      
      //=========== Get the EMCAL branch from Hits Tree for the Primary track itrack
      branch->GetEntry(itrack,0);
      AliEMCAL * EMCAL = (AliEMCAL *) gAlice->GetDetector("EMCAL") ;
      AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance( EMCAL->GetGeometry()->GetName(), EMCAL->GetGeometry()->GetTitle() );
 
      Int_t i;
      for ( i = 0 ; i < fHits->GetEntries() ; i++ ) {
	AliEMCALHit * hit = (AliEMCALHit*)fHits->At(i) ;
        AliEMCALDigit curSDigit = AliEMCALDigit(1,1,1,1);
        AliEMCALDigit *sdigit ;
 
	// Assign primary number only if contribution is significant
     
        if( hit->GetEnergy() > fPrimThreshold)
          curSDigit =  AliEMCALDigit( hit->GetPrimary(), hit->GetIparent(), (((hit->GetId()/geom->GetNPhi())%geom->GetNZ()+1) * (hit->GetId()%(geom->GetNPhi()+1))), Digitize( hit->GetEnergy() ) ) ;
        else
          curSDigit =  AliEMCALDigit( -1               , -1               ,(((hit->GetId()/geom->GetNPhi())%geom->GetNZ()+1) * (hit->GetId()%(geom->GetNPhi()+1))), Digitize( hit->GetEnergy() ) ) ;
	
        for(Int_t check= 0; check < nSdigits; check++) {
          sdigit = (AliEMCALDigit *)fSDigits->At(check);
          if( (((hit->GetId()/geom->GetNPhi())%geom->GetNZ()) == ((sdigit->GetId()/geom->GetNPhi())%geom->GetNZ()))  && ((hit->GetId()%geom->GetNPhi()) == (sdigit->GetId()%geom->GetNPhi()))) 
           { 
             *sdigit = *sdigit + curSDigit ;
           }
     else 
         { new((*fSDigits)[nSdigits])  AliEMCALDigit(curSDigit);
	  nSdigits++ ; } 
	}
  

        if( hit->GetEnergy() > fPrimThreshold)
          curSDigit =  AliEMCALDigit( hit->GetPrimary(), hit->GetIparent(), ((geom->GetNZ() * geom->GetNPhi()) + ((hit->GetId()/geom->GetNPhi())%geom->GetNZ()+1) * (hit->GetId()%(geom->GetNPhi()+1))), Digitize( hit->GetEnergy() ) ) ;
        else
          curSDigit =  AliEMCALDigit( -1               , -1               ,((geom->GetNZ() * geom->GetNPhi()) + ((hit->GetId()/geom->GetNPhi())%geom->GetNZ()+1) * (hit->GetId()%(geom->GetNPhi()+1))), Digitize( hit->GetEnergy() ) ) ;
 
      if((hit->GetId()/geom->GetNPhi()) < (2*geom->GetNZ())) 
       {
        for(Int_t check= 0; check < nSdigits; check++) {
          sdigit = (AliEMCALDigit *)fSDigits->At(check);
          if( (((hit->GetId()/geom->GetNPhi())%geom->GetNZ()) == ((sdigit->GetId()/geom->GetNPhi())%geom->GetNZ()))  && ((hit->GetId()%geom->GetNPhi()) == (sdigit->GetId()%geom->GetNPhi()))) 
           { 
             *sdigit = *sdigit + curSDigit ;
           }
     else 
         { new((*fSDigits)[nSdigits])  AliEMCALDigit(curSDigit);
	  nSdigits++ ; } 
	}
      } 
     } 
    } // loop over tracks
    
    fSDigits->Sort() ;
    
    nSdigits = fSDigits->GetEntriesFast() ;
    fSDigits->Expand(nSdigits) ;
    
    Int_t i ;
    for (i = 0 ; i < nSdigits ; i++) { 
      AliEMCALDigit * digit = (AliEMCALDigit *) fSDigits->At(i) ; 
      digit->SetIndexInList(i) ;     
    }

    if(gAlice->TreeS() == 0)
      gAlice->MakeTree("S") ;
    
    //check, if this branch already exits?
    TBranch * sdigitsBranch = 0;
    TBranch * sdigitizerBranch = 0;
    
    TObjArray * branches = gAlice->TreeS()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t emcalNotFound = kTRUE ;
    Bool_t sdigitizerNotFound = kTRUE ;
    
    for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
      
      if(emcalNotFound){
	sdigitsBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp("EMCAL",sdigitsBranch->GetName())==0 ) &&
	    (fSDigitsTitle.CompareTo(sdigitsBranch->GetTitle()) == 0) )
	  emcalNotFound = kFALSE ;
      }
      if(sdigitizerNotFound){
	sdigitizerBranch = (TBranch *) branches->At(ibranch) ;
	if( (strcmp(sdigitizerBranch->GetName(),"AliEMCALSDigitizer") == 0)&&
	    (fSDigitsTitle.CompareTo(sdigitizerBranch->GetTitle()) == 0) )
	  sdigitizerNotFound = kFALSE ;
      }
    }

    if(!(sdigitizerNotFound && emcalNotFound)){
      cout << "AliEMCALSdigitizer error:" << endl ;
      cout << "Can not overwrite existing branches: do not write" << endl ;
      return ;
    }
    
    //Make (if necessary) branches    
    char * file =0;
    if(gSystem->Getenv("CONFIG_SPLIT_FILE")){ //generating file name
      file = new char[strlen(gAlice->GetBaseFile())+20] ;
      sprintf(file,"%s/EMCAL.SDigits.root",gAlice->GetBaseFile()) ;
    }
    
    TDirectory *cwd = gDirectory;
    
    //First list of sdigits
    Int_t bufferSize = 32000 ;    
    sdigitsBranch = gAlice->TreeS()->Branch("EMCAL",&fSDigits,bufferSize);
    sdigitsBranch->SetTitle(fSDigitsTitle.Data());
    if (file) {
      sdigitsBranch->SetFile(file);
      TIter next( sdigitsBranch->GetListOfBranches());
      TBranch * subbr;
      while ((subbr=(TBranch*)next())) {
	subbr->SetFile(file);
      }   
      cwd->cd();
    } 
      
    //second - SDigitizer
    Int_t splitlevel = 0 ;
    AliEMCALSDigitizer * sd = this ;
    sdigitizerBranch = gAlice->TreeS()->Branch("AliEMCALSDigitizer","AliEMCALSDigitizer",
					       &sd,bufferSize,splitlevel); 
    sdigitizerBranch->SetTitle(fSDigitsTitle.Data());
    if (file) {
      sdigitizerBranch->SetFile(file);
      TIter next( sdigitizerBranch->GetListOfBranches());
      TBranch * subbr ;
      while ((subbr=(TBranch*)next())) {
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
    gBenchmark->Stop("EMCALSDigitizer");
    cout << "AliEMCALSDigitizer:" << endl ;
    cout << "   took " << gBenchmark->GetCpuTime("EMCALSDigitizer") << " seconds for SDigitizing " 
	 <<  gBenchmark->GetCpuTime("EMCALSDigitizer")/fNevents << " seconds per event " << endl ;
    cout << endl ;
  }
  
  
}
//__________________________________________________________________
void AliEMCALSDigitizer::SetSDigitsBranch(const char * title ){
  //Seting title to branch SDigits 
  if(!fSDigitsTitle.IsNull())
    cout << "AliEMCALSdigitizer: changing SDigits file from " <<fSDigitsTitle.Data() << " to " << title << endl ;
  fSDigitsTitle=title ;
}
//__________________________________________________________________
void AliEMCALSDigitizer::Print(Option_t* option)const{
  cout << "------------------- "<< GetName() << " -------------" << endl ;
  cout << "   Writing SDigitis to branch with title  " << fSDigitsTitle.Data() << endl ;
  cout << "   with digitization parameters  A = " << fA << endl ;
  cout << "                                 B = " << fB << endl ;
  cout << "   Threshold for Primary assignment= " << fPrimThreshold << endl ; 
  cout << "---------------------------------------------------"<<endl ;
  
}
//__________________________________________________________________
Bool_t AliEMCALSDigitizer::operator==( AliEMCALSDigitizer const &sd )const{
  if( (fA==sd.fA)&&(fB==sd.fB)&&(fPrimThreshold==sd.fPrimThreshold))
    return kTRUE ;
  else
    return kFALSE ;
}
//__________________________________________________________________
void AliEMCALSDigitizer::PrintSDigits(Option_t * option){
  //Prints list of digits produced at the current pass of AliEMCALDigitizer
  
  cout << "AliEMCALSDigitizer: " << endl ;
  cout << "       Number of entries in SDigits list  " << fSDigits->GetEntriesFast() << endl ;
  cout << endl ;
  
  if(strstr(option,"all")){// print all digits
    
    //loop over digits
    AliEMCALDigit * digit;
    cout << "SDigit Id " << " Amplitude " <<  " Index "  <<  " Nprim " << " Primaries list " <<  endl;    
    Int_t index ;
    for (index = 0 ; index < fSDigits->GetEntries() ; index++) {
      digit = (AliEMCALDigit * )  fSDigits->At(index) ;
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
