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

/* $Id:  */

/* $Log:
*/

//_________________________________________________________________________
//  A singleton. This class should be used in the analysis stage to get 
//  reconstructed objects: Digits, RecPoints, TrackSegments and RecParticles,
//  instead of directly reading them from galice.root file. This container 
//  ensures, that one reads Digits, made of these particular digits, RecPoints, 
//  made of these particular RecPoints, TrackSegments and RecParticles. 
//  This becomes non trivial if there are several identical branches, produced with
//  different set of parameters. Currently This class only Retrieves Hits, Digits, and SDigits. 
//
//  An example of how to use (see also class AliEMCALAnalyser):
//  AliEMCALGetter * gime = AliEMCALGetter::GetInstance("galice.root","test") ;
//     ................
//  please->GetEvent(event) ;    // reads new event from galice.root
//                  
//*-- Author: Sahal Yacoob (LBL) 
// based on : AliPHOSGetter
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TObjString.h"
#include "TFolder.h"

// --- Standard library ---
#include <iostream.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliConfig.h"
#include "AliEMCALGetter.h"
#include "AliEMCALv1.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"

ClassImp(AliEMCALGetter)
  
  AliEMCALGetter * AliEMCALGetter::fgObjGetter = 0 ; 

//____________________________________________________________________________ 
AliEMCALGetter::AliEMCALGetter(const char* headerFile, const char* branchTitle )
{
  //Initialize  all lists

  fHeaderFile         = headerFile ; 
  fBranchTitle        = branchTitle ;
  fSDigitsTitle       = branchTitle ; 
  fDigitsTitle        = branchTitle ; 

  fPrimaries = new TObjArray(1) ;
  fModuleFolder  = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Configuration/Modules"));
  fHitsFolder    = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/Hits")); 
  fSDigitsFolder = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/SDigits")); 
  fDigitsFolder  = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Event/Data")); 
  fTasksFolder   = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Tasks")) ; 

  if ( fHeaderFile != "aliroot"  ) { // to call the getter without a file

    //open headers file
    TFile * file = static_cast<TFile*>(gROOT->GetFile(fHeaderFile.Data() ) ) ;
    
    if(file == 0){    //if file was not opened yet, read gAlice
      if(fHeaderFile.Contains("rfio")) // if we read file using HPSS
	file =	TFile::Open(fHeaderFile.Data(),"update") ;
      else
	file = new TFile(fHeaderFile.Data(),"update") ;
      
      if (!file->IsOpen()) {
	cerr << "ERROR : AliEMCALGetter::AliEMCALGetter -> Cannot open " << fHeaderFile.Data() << endl ; 
	abort() ; 
      }
      
      gAlice = static_cast<AliRun *>(file->Get("gAlice")) ;
      
      if (!gAlice) {
	cerr << "ERROR : AliEMCALGetter::AliEMCALGetter -> Cannot find gAlice in " << fHeaderFile.Data() << endl ; 
	abort() ; 
      }
      if (!EMCAL()) {
	  cout << "INFO: AliEMCALGetter -> Posting EMCAL to Folders" << endl ; 
	AliConfig * conf = AliConfig::Instance() ; 
	conf->Add(static_cast<AliDetector*>(gAlice->GetDetector("EMCAL"))) ; 
 	conf->Add(static_cast<AliModule*>(gAlice->GetDetector("EMCAL"))) ; 
     }
    }

  }
}
//____________________________________________________________________________ 
AliEMCALGetter::~AliEMCALGetter(){

}

//____________________________________________________________________________ 
void AliEMCALGetter::CreateWhiteBoard() const
{

}

//____________________________________________________________________________ 
AliEMCALGetter * AliEMCALGetter::GetInstance()
{
  // Returns the pointer of the unique instance already defined
  
  AliEMCALGetter * rv = 0 ;
  if ( fgObjGetter )
    rv = fgObjGetter ;
  else
    cout << "AliEMCALGetter::GetInstance ERROR: not yet initialized" << endl ;

  return rv ;
}

//____________________________________________________________________________ 
AliEMCALGetter * AliEMCALGetter::GetInstance(const char* headerFile,
					   const char* branchTitle)
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed 

  if ( fgObjGetter )    
    if((fgObjGetter->fBranchTitle.CompareTo(branchTitle) == 0) && 
       (fgObjGetter->fHeaderFile.CompareTo(headerFile)==0))
      return fgObjGetter ;
    else
      fgObjGetter->~AliEMCALGetter() ;  // delete it if already exists another version
  
  fgObjGetter = new AliEMCALGetter(headerFile,branchTitle) ; 
  
  // Posts a few item to the white board (folders)
  // fgObjGetter->CreateWhiteBoard() ;
    
  return fgObjGetter ; 
  
}

//____________________________________________________________________________ 
const AliEMCALv0 * AliEMCALGetter::EMCAL() 
{
  // returns the EMCAL object 
  //AliEMCALv0 * emcal = dynamic_cast<AliEMCALv0 *>(gAlice->GetDetector("EMCAL")) ; 
  AliEMCALv0 * emcal = dynamic_cast<AliEMCALv1 *>(fModuleFolder->FindObject("EMCAL")) ;  
  if (!emcal) 
      cout << "WARNING: AliEMCALGetter::EMCAL -> EMCAL module not found in Folders" << endl ; 
  return emcal ; 
}  

//____________________________________________________________________________ 
const AliEMCALGeometry * AliEMCALGetter::EMCALGeometry() 
{
  AliEMCALGeometry * rv = 0 ; 
  if (EMCAL() )
   rv =  EMCAL()->GetGeometry() ;
  return rv ; 
} 

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostHits(void) const 
{  //------- Hits ----------------------

  // the hierarchy is //Folders/RunMC/Event/Data/EMCAL/Hits
  
  TFolder * emcalFolder = dynamic_cast<TFolder*>(fHitsFolder->FindObject("EMCAL")) ; 
  if ( !emcalFolder ) {
      cout << "WARNING: AliEMCALGetter::Post H -> Folder //" << fHitsFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post H -> Adding Folder //" << fHitsFolder << "/EMCAL/"  << endl;
    emcalFolder = fHitsFolder->AddFolder("EMCAL", "Hits from EMCAL") ; 
  }    
  TClonesArray *hits=  new TClonesArray("AliEMCALHit",1000) ;
  hits->SetName("Hits") ;
  emcalFolder->Add(hits) ; 
  
  return kTRUE;
} 

//____________________________________________________________________________ 
void * AliEMCALGetter::HitsRef(void) const 
{  //------- Hits ----------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/EMCAL/Hits
  if ( !fHitsFolder ) {
    cerr << "ERROR: AliEMCALGetter::Post H -> Folder //" << fHitsFolder << " not found!" << endl;
    return 0;
  }    
 
  TFolder * emcalFolder = dynamic_cast<TFolder *>(fHitsFolder->FindObject("EMCAL")) ;
  if ( !emcalFolder ) {
    cerr << "ERROR: AliEMCALGetter::Post HRef -> Folder //" << fHitsFolder << "/EMCAL/ not found!" << endl;  
    return 0;
  }
 
  TObject * h = emcalFolder->FindObject("Hits") ;
  if(!h) {
    cerr << "ERROR: AliEMCALGetter::HRef -> " << emcalFolder->GetName() << "/Hits not found !" << endl ; 
    return 0 ;
  }
  else
    return static_cast<void *>(emcalFolder->GetListOfFolders()->GetObjectRef(h)) ;
}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostSDigits(const char * name, const char * headerFile) const 
{  //---------- SDigits -------------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/EMCAL/SDigits/headerFile/sdigitsname
  // because you can have sdigits from several hit files for mixing
  
  TFolder * emcalFolder = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("EMCAL")) ;
  if ( !emcalFolder ) {
      cout << "WARNING: AliEMCALGetter::Post S -> Folder //" << fSDigitsFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post S -> Adding Folder //" << fHitsFolder << "/EMCAL/" << endl;
    emcalFolder = fSDigitsFolder->AddFolder("EMCAL", "SDigits from EMCAL") ; 
  }    
  TString subdir(headerFile) ;
  TFolder * emcalSubFolder = dynamic_cast<TFolder*>(emcalFolder->FindObject(subdir)) ; 
  if ( !emcalSubFolder ) 
    emcalSubFolder = emcalFolder->AddFolder(subdir, ""); 
  
  TObject * sd  = emcalSubFolder->FindObject(name); 
  if ( sd ) {
      cerr <<"INFO: AliEMCALGetter::Post S -> Folder " << subdir 
	   << " already exists!" << endl ;  
  }else{
    TClonesArray * sdigits = new TClonesArray("AliEMCALDigit",1000) ;
    sdigits->SetName(name) ;
    emcalSubFolder->Add(sdigits) ;
  }
  
  return kTRUE;
} 
//____________________________________________________________________________ 
void * AliEMCALGetter::SDigitsRef(const char * name, const char * file) const 
{  //------- SDigits ----------------------
  
  // the hierarchy is //Folders/RunMC/Event/Data/EMCAL/SDigits/filename/SDigits

  if ( !fSDigitsFolder ) {
    cerr << "ERROR: AliEMCALGetter::Post SRef -> Folder //" << fSDigitsFolder << " not found!" << endl;
    return 0;
  }    
 
  TFolder * emcalFolder = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject("EMCAL")) ;
  if ( !emcalFolder ) {
    cerr << "ERROR: AliEMCALGetter::Post SRef -> Folder //" << fSDigitsFolder << "/EMCAL/ not found!" << endl;
    return 0;
  }

  TFolder * emcalSubFolder = 0 ;
  if(file)
    emcalSubFolder = dynamic_cast<TFolder *>(emcalFolder->FindObject(file)) ;
  else
    emcalSubFolder = dynamic_cast<TFolder *>(emcalFolder->FindObject(fHeaderFile)) ;
  
  if(!emcalSubFolder) {
    cerr << "ERROR: AliEMCALGetter::Post SRef -> Folder //Folders/RunMC/Event/Data/EMCAL/" << file << "not found!" << endl;
    return 0;
  }

  TObject * dis = emcalSubFolder->FindObject(name) ;
  if(!dis)
    return 0 ;
  else
    return static_cast<void *>(emcalSubFolder->GetListOfFolders()->GetObjectRef(dis)) ;

}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostSDigitizer(AliEMCALSDigitizer * sdigitizer) const 
{  //---------- SDigitizer -------------------------
    
  // the hierarchy is //Folders/Tasks/SDigitizer/EMCAL/sdigitsname


  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ; 

  if ( !sd ) {
    cerr << "ERROR: AliEMCALGetter::Post Ser -> Task //" << fTasksFolder << "/SDigitizer not found!" << endl;
    return kFALSE ;
  }        
  TTask * emcal = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
      cout <<"WARNING: AliEMCALGetter::Post Ser ->//" << fTasksFolder << "/SDigitizer/EMCAL/ not found!" << endl;  
      cout <<"INFO: AliEMCALGetter::Post Ser -> Adding //" << fTasksFolder << "/SDigitizer/EMCAL/" << endl;
    emcal = new TTask("EMCAL", "") ; 
    sd->Add(emcal) ; 
  } 
  AliEMCALSDigitizer * emcalsd  = dynamic_cast<AliEMCALSDigitizer *>(emcal->GetListOfTasks()->FindObject( sdigitizer->GetName() )); 
  if (emcalsd) { 
      cout << "INFO: AliEMCALGetter::Post Ser -> Task " << sdigitizer->GetName() << " already exists" << endl ; 
    emcal->GetListOfTasks()->Remove(emcalsd) ;
  }
  emcal->Add(sdigitizer) ;	
  return kTRUE; 
  
}

//____________________________________________________________________________ 
void * AliEMCALGetter::SDigitizerRef(const char * name) const 
{  

  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ; 
  if ( !sd ) {
    cerr << "ERROR: AliEMCALGetter::Post SerRef -> Task //" << fTasksFolder << "/SDigitizer not found!" << endl;
    abort();
  }        

  TTask * emcal = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    cerr <<"ERROR: AliEMCALGetter::Post SerRef ->  //" << fTasksFolder << "/SDigitizer/EMCAL not found!" << endl;
    abort();
  }        

  TTask * task = dynamic_cast<TTask*>(emcal->GetListOfTasks()->FindObject(name)) ; 

  return static_cast<void *>(emcal->GetListOfTasks()->GetObjectRef(task)) ;

}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostSDigitizer(const char * name, const char * file) const 
{  //---------- SDigitizer -------------------------
  
 // the hierarchy is //Folders/Tasks/SDigitizer/EMCAL/sdigitsname


  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ; 
  if ( !sd ) {
    cerr << "ERROR: AliEMCALGetter::Post Ser -> Task //" << fTasksFolder << "/SDigitizer not found!" << endl;
    return kFALSE ;
  }        

  TTask * emcal = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
      cout <<"WARNING: AliEMCALGetter::Post Ser ->  //" << fTasksFolder << "/SDigitizer/EMCAL/ not found!" << endl;
      cout <<"INFO: AliEMCALGetter::Post Ser -> Adding  //" << fTasksFolder << "/SDigitizer/EMCAL" << endl;
    emcal = new TTask("EMCAL", "") ; 
    sd->Add(emcal) ; 
  } 

  TString sdname(name) ;
  sdname.Append(":") ;
  sdname.Append(file);
  AliEMCALSDigitizer * emcalsd  = dynamic_cast<AliEMCALSDigitizer *>(emcal->GetListOfTasks()->FindObject( sdname )); 
  if (!emcalsd) {
    emcalsd = new AliEMCALSDigitizer() ;  
    //Note, we can not call constructor with parameters: it will call Getter and scrud up everething
    emcalsd->SetName(sdname) ;
    emcalsd->SetTitle(file) ;
    emcal->Add(emcalsd) ;	
  }
  return kTRUE; 
  
}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostDigits(const char * name) const 
{  //---------- Digits -------------------------

  // the hierarchy is //Folders/Run/Event/Data/EMCAL/SDigits/name

  TFolder * emcalFolder  = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("EMCAL")) ;

  if ( !emcalFolder ) {
      cout << "WARNING: AliEMCALGetter::Post D -> Folder //" << fDigitsFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post D -> Adding Folder //" << fDigitsFolder << "/EMCAL/" << endl;
    emcalFolder = fDigitsFolder->AddFolder("EMCAL", "Digits from EMCAL") ;  
  }    
 
  TObject*  dig = emcalFolder->FindObject( name ) ;
  if ( !dig ) {
    TClonesArray * digits = new TClonesArray("AliEMCALDigit",1000) ;
    digits->SetName(name) ;
    emcalFolder->Add(digits) ;  
  }
  return kTRUE; 
}

//____________________________________________________________________________ 
void * AliEMCALGetter::DigitsRef(const char * name) const 
{ //------- Digits ----------------------
  
  // the hierarchy is //Folders/Run/Event/Data/EMCAL/Digits/name

  if ( !fDigitsFolder ) {
    cerr << "ERROR: AliEMCALGetter::Post DRef -> Folder //" << fDigitsFolder << " not found!" << endl;
    return 0;
  }    
  
  TFolder * emcalFolder  = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("EMCAL")) ; 
  if ( !emcalFolder ) {
    cerr << "ERROR: AliEMCALGetter::DRef -> Folder //" << fDigitsFolder << "/EMCAL/ not found!" << endl;
    return 0;
  }    

  TObject * d = emcalFolder->FindObject(name) ;
  if(!d)
    return 0 ;
  else
    return static_cast<void *>(emcalFolder->GetListOfFolders()->GetObjectRef(d)) ;

}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostDigitizer(AliEMCALDigitizer * digitizer) const 
{  //---------- Digitizer -------------------------
  
  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 

  if ( !sd ) {
    cerr << "ERROR: AliEMCALGetter::Post Der -> Task //" << fTasksFolder << "/Digitizer not found!" << endl;
    return kFALSE ;
  }        
  TTask * emcal = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
      cout <<"WARNING: AliEMCALGetter::Post Der ->  //" << fTasksFolder << "/Digitizer/EMCAL not found!" << endl;
      cout <<"INFO: AliEMCALGetter::Post Der -> Adding //" << fTasksFolder << "/Digitizer/EMCAL" << endl; 
    emcal = new TTask("EMCAL", "") ; 
    sd->Add(emcal) ; 
  } 

    AliEMCALDigitizer * emcald = dynamic_cast<AliEMCALDigitizer*>(emcal->GetListOfTasks()->FindObject(digitizer->GetName())) ; 
    if (emcald) { 
      emcald->Delete() ;
      emcal->GetListOfTasks()->Remove(emcald) ;
    }
    emcal->Add(digitizer) ; 
    return kTRUE; 
}  

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostDigitizer(const char * name) const 
{  //---------- Digitizer -------------------------
  
 // the hierarchy is //Folders/Tasks/SDigitizer/EMCAL/sdigitsname

  TTask * d  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 
  if ( !d ) {
    cerr << "ERROR: AliEMCALGetter::Post Der -> Task //" << fTasksFolder << "/Digitizer not found!" << endl;
    return kFALSE ;
  }        

  TTask * emcal = dynamic_cast<TTask*>(d->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
      cout <<"WARNING: AliEMCALGetter::Post Der -> //" << fTasksFolder << "/Digitizer/EMCAL not found!" << endl; 
      cout <<"INFO: AliEMCALGetter::Post Der -> Adding //" << fTasksFolder << "/Digitizer/EMCAL" << endl;
    emcal = new TTask("EMCAL", "") ; 
    d->Add(emcal) ; 
} 

  AliEMCALDigitizer * emcald = dynamic_cast<AliEMCALDigitizer*>(emcal->GetListOfTasks()->FindObject(name)) ; 
  if (!emcald) { 
    emcald = new AliEMCALDigitizer() ;
    emcald->SetName(fDigitsTitle) ;
    emcald->SetTitle(fHeaderFile) ;
    emcal->Add(emcald) ;
  }
  return kTRUE;  
}

//____________________________________________________________________________ 
void * AliEMCALGetter::DigitizerRef(const char * name) const 
{  
  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 
  if ( !sd ) {
    cerr << "ERROR: AliEMCALGetter::Post DerRef -> Task //" << fTasksFolder << "/Digitizer not found!" << endl;
    abort();
  }        

  TTask * emcal = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    cerr <<"ERROR: AliEMCALGetter::Post DerRef ->  //" << fTasksFolder << "/Digitizer/EMCAL" << endl;
    abort();
  }        

  TTask * task = dynamic_cast<TTask*>(emcal->GetListOfTasks()->FindObject(name)) ; 

  return static_cast<void *>(emcal->GetListOfTasks()->GetObjectRef(task)) ;

}
 
//____________________________________________________________________________ 
const TParticle * AliEMCALGetter::Primary(Int_t index) const
{
  // Return primary particle numbered by <index>

  if(index < 0) 
    return 0 ;
  
  Int_t primaryIndex = index % 10000000 ; 
  Int_t primaryList = (Int_t ) ((index-primaryIndex)/10000000.)  ;
  
  if ( primaryList > 0  ) {
      cout << " Getter does not support currently Mixing of primary " << endl ;
      cout << "   can not return primary: " << index<< " (list "<< primaryList<< " primary # " << primaryIndex << " )"<<endl ;
    return 0;
  }
  
  return gAlice->Particle(primaryIndex) ;
  
}

//____________________________________________________________________________ 
void AliEMCALGetter::ReadTreeD()
{
  // Read the digit tree gAlice->TreeD()  
  if(gAlice->TreeD()== 0){
    cerr <<   "ERROR: AliEMCALGetter::ReadTreeD: can not read TreeD " << endl ;
  return ;
  }
 cout << "hello" << endl;  
  TObjArray * lob = static_cast<TObjArray*>(gAlice->TreeD()->GetListOfBranches()) ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * digitsbranch = 0 ; 
  TBranch * digitizerbranch = 0 ; 
  Bool_t emcalfound = kFALSE, digitizerfound = kFALSE ; 
  
  while ( (branch = static_cast<TBranch*>(next())) && (!emcalfound || !digitizerfound) ) {
    if ( (strcmp(branch->GetName(), "EMCAL")==0) && (strcmp(branch->GetTitle(), fDigitsTitle)==0) ) {
      digitsbranch = branch ; 
      emcalfound = kTRUE ;
    }
    else if ( (strcmp(branch->GetName(), "AliEMCALDigitizer")==0) && (strcmp(branch->GetTitle(), fDigitsTitle)==0) ) {
      digitizerbranch = branch ; 
      digitizerfound = kTRUE ; 
    }
  }

  if ( !emcalfound || !digitizerfound ) {
      cout << "WARNING: AliEMCALGetter::ReadTreeD -> Cannot find Digits and/or Digitizer with name " 
	   << fDigitsTitle << endl ;
    return ; 
  }   
 
  //read digits
  if(!Digits(fDigitsTitle) ) 
    PostDigits(fDigitsTitle);
  digitsbranch->SetAddress(DigitsRef(fDigitsTitle)) ;
  digitsbranch->GetEntry(0) ;
  
  
  // read  the Digitizer
  if(!Digitizer(fDigitsTitle))
    PostDigitizer(fDigitsTitle) ;
  digitizerbranch->SetAddress(DigitizerRef(fDigitsTitle)) ;
  digitizerbranch->GetEntry(0) ;
 
  
}

//____________________________________________________________________________ 
void AliEMCALGetter::ReadTreeH()
{
  // Read the first entry of EMCAL branch in hit tree gAlice->TreeH()

  if(gAlice->TreeH()== 0){
    cerr <<   "ERROR: AliEMCALGetter::ReadTreeH: -> Cannot read TreeH " << endl ;
    return ;
  }
  
  TBranch * hitsbranch = static_cast<TBranch*>(gAlice->TreeH()->GetBranch("EMCAL")) ;
  if ( !hitsbranch ) {
      cout << "WARNING:  AliEMCALGetter::ReadTreeH -> Cannot find branch EMCAL" << endl ; 
    return ;
  }
  if(!Hits())
    PostHits() ;

  hitsbranch->SetAddress(HitsRef()) ;

  hitsbranch->GetEntry(0) ;

}

//____________________________________________________________________________ 
void AliEMCALGetter::Track(Int_t itrack)
{
  // Read the first entry of EMCAL branch in hit tree gAlice->TreeH()

  if(gAlice->TreeH()== 0){
    cerr <<   "ERROR: AliEMCALGetter::ReadTreeH: -> Cannot read TreeH " << endl ;
    return ;
  }
  
  TBranch * hitsbranch = dynamic_cast<TBranch*>(gAlice->TreeH()->GetListOfBranches()->FindObject("EMCAL")) ;
  if ( !hitsbranch ) {
      cout << "WARNING:  AliEMCALGetter::ReadTreeH -> Cannot find branch EMCAL" << endl ; 
    return ;
  }  
  if(!Hits())
    PostHits() ;
  hitsbranch->SetAddress(HitsRef()) ;
  hitsbranch->GetEntry(itrack) ;


}
//____________________________________________________________________________ 
void AliEMCALGetter::ReadTreeS(Int_t event)
{
  // Read the summable digits tree gAlice->TreeS()  
  
  // loop over all opened files and read their SDigits to the White Board
  TFolder * emcalF = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject("EMCAL")) ;
  if (!emcalF) 
    emcalF = fSDigitsFolder->AddFolder("EMCAL", "SDigits from EMCAL") ; 
  TCollection * folderslist = emcalF->GetListOfFolders() ; 
  
  //Add current file to list if it is not there yet
  if ( (fHeaderFile != "aliroot") && ( !folderslist->Contains(fHeaderFile) ) ){
    emcalF->AddFolder(fHeaderFile, ""); 
  }
    
  TIter next(folderslist) ; 
  TFolder * folder = 0 ; 
  TFile * file; 
  TTree * treeS = 0;
  while ( (folder = static_cast<TFolder*>(next())) ) {
    if(fHeaderFile.CompareTo(folder->GetName()) == 0 ) 
      {treeS=gAlice->TreeS() ;
        cout << "ReadTreeS  "<<  gAlice->TreeS()  <<endl ;}
    else{
     cout << " AliEMCALGetter::ReadTreeS 2 " <<  folder->GetName() << endl ; 
     file = static_cast<TFile*>(gROOT->GetFile(folder->GetName())); 
      file->cd() ;
      
      // Get SDigits Tree header from file
      TString treeName("TreeS") ;
      treeName += event ; 
      treeS = dynamic_cast<TTree*>(gDirectory->Get(treeName.Data()));
    }
    if(treeS==0){
      cerr << "ERROR: AliEMCALGetter::ReadTreeS There is no SDigit Tree" << endl;
      return ;
    }
    
    //set address of the SDigits and SDigitizer
    TBranch   * sdigitsBranch    = 0;
    TBranch   * sdigitizerBranch = 0;
    TBranch   * branch           = 0 ;  
    TObjArray * lob = static_cast<TObjArray*>(treeS->GetListOfBranches()) ;
    TIter next(lob) ; 
    Bool_t emcalfound = kFALSE, sdigitizerfound = kFALSE ; 
    
    while ( (branch = static_cast<TBranch*>(next())) && (!emcalfound || !sdigitizerfound) ) {
 
     if ( (strcmp(branch->GetName(), "EMCAL")==0) && (strcmp(branch->GetTitle(), fSDigitsTitle)==0) ) {
	emcalfound = kTRUE ;
	sdigitsBranch = branch ; 
       cout << "sdigitsbranch found = " << branch << endl ; 
       }
      
      else if ( (strcmp(branch->GetName(), "AliEMCALSDigitizer")==0) && (strcmp(branch->GetTitle(), fSDigitsTitle)==0) ) {
	sdigitizerfound = kTRUE ; 
	sdigitizerBranch = branch ;
       cout << "sdigitizerbranch found = " << branch << endl ; 
      }
    }
    if ( !emcalfound || !sdigitizerfound ) {
	cout << "WARNING: AliEMCALDigitizer::ReadSDigits -> Digits and/or Digitizer branch with name " << GetName() 
	     << " not found" << endl ;
      return ; 
    }   
    
    if ( !folder->FindObject(fSDigitsTitle) )  
     { PostSDigits(fSDigitsTitle,folder->GetName()) ;
       cout << "Posting SDigits " << endl << endl ;}  
    sdigitsBranch->SetAddress(SDigitsRef(fSDigitsTitle,folder->GetName())) ;
    
    sdigitsBranch->GetEntry(0) ;
    
    TString sdname(fSDigitsTitle) ;
    cout << sdname << endl ;
    sdname+=":" ;
    sdname+=folder->GetName() ;
    if(!SDigitizer(sdname) ) 
      PostSDigitizer(fSDigitsTitle,folder->GetName()) ;
    sdigitizerBranch->SetAddress(SDigitizerRef(sdname)) ;
    sdigitizerBranch->GetEntry(0) ;
    
  }    
  
  // After SDigits have been read from all files, return to the first one
  
  next.Reset();
  folder = static_cast<TFolder*>(next());
  if(folder){
    file   = static_cast<TFile*>(gROOT->GetFile(folder->GetName())); 
    file   ->cd() ;
  }
  
}
//____________________________________________________________________________ 
void AliEMCALGetter::ReadTreeS(TTree * treeS, Int_t input)
{  // Read the summable digits fron treeS()  


  TString filename("mergefile") ;
  filename+= input ;

  TFolder * emcalFolder = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("EMCAL")) ; 
  if ( !emcalFolder ) { 
   emcalFolder = fSDigitsFolder->AddFolder("EMCAL", "SDigits from EMCAL") ; 
  } 
  TFolder * folder=(TFolder*)emcalFolder->FindObject(filename) ;
  //set address of the SDigits and SDigitizer
  TBranch   * sdigitsBranch    = 0;
  TBranch   * sdigitizerBranch = 0;
  TBranch   * branch           = 0 ;  
  TObjArray * lob = (TObjArray*)treeS->GetListOfBranches() ;
  TIter next(lob) ; 
  Bool_t emcalfound = kFALSE, sdigitizerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!emcalfound || !sdigitizerfound) ) {
    if ( strcmp(branch->GetName(), "EMCAL")==0) {
      emcalfound = kTRUE ;
      sdigitsBranch = branch ; 
    }
    
    else if ( strcmp(branch->GetName(), "AliEMCALSDigitizer")==0) {
      sdigitizerfound = kTRUE ; 
      sdigitizerBranch = branch ;
    }
  }
  if ( !emcalfound || !sdigitizerfound ) {
      cout << "WARNING: AliEMCALGetter::ReadTreeS -> Digits and/or Digitizer branch not found" << endl ;
    return ; 
  }   
  
  if (!folder || !(folder->FindObject(sdigitsBranch->GetTitle()) ) )
    PostSDigits(sdigitsBranch->GetTitle(),filename) ;

  sdigitsBranch->SetAddress(SDigitsRef(sdigitsBranch->GetTitle(),filename)) ;
  
  TString sdname(sdigitsBranch->GetTitle()) ;
  sdname+=":" ;
  sdname+=filename ;
  if(!SDigitizer(sdigitsBranch->GetTitle()) )
    PostSDigitizer(sdigitsBranch->GetTitle(),filename) ;
  sdigitizerBranch->SetAddress(SDigitizerRef(sdname)) ;
  
  sdigitsBranch->GetEntry(0) ;
  sdigitizerBranch->GetEntry(0) ;
  
}    


//____________________________________________________________________________ 
void AliEMCALGetter::ReadPrimaries()
{
  // Reads specific branches of primaries
  
  fNPrimaries = gAlice->GetNtrack();
  
  //   //Check, is it necessary to open new files
  //   TArrayI* events = fDigitizer->GetCurrentEvents() ; 
  //   TClonesArray * filenames = fDigitizer->GetHeadersFiles() ;
//   Int_t input ;
//   for(input = 0; input < filenames->GetEntriesFast(); input++){

//     TObjString * filename = (TObjString *) filenames->At(input) ;

//     //Test, if this file already open
//     TFile *file = (TFile*) gROOT->GetFile( filename->GetString() ) ;
//     if(file == 0)
//       file = new TFile( filename->GetString()) ;
//     file->cd() ;
    
//     // Get Kine Tree from file
// //     char treeName[20];
// //     sprintf(treeName,"TreeK%d",events->At(input));
// //     TTree * treeK = (TTree*)gDirectory->Get(treeName);
// //     if (treeK) 
// //       treeK->SetBranchAddress("Particles", &fParticleBuffer);
// //     else    
// //       cout << "AliEMCALGetter: cannot find Kine Tree for event:" << events->At(input) << endl;

// //     // Create the particle stack
// //     if(!fParticles) fParticles = new TClonesArray("TParticle",1000);
// //     // Build the pointer list
// //     if(fParticleMap) {     <----
// //       fParticleMap->Clear();
// //       fParticleMap->Expand(treeK->GetEntries());
// //     } else
// //       fParticleMap = new TObjArray(treeK->GetEntries());
    
//     // From gAlice->Particle(i) 


// //   if(!(*fParticleMap)[i]) {
// //     Int_t nentries = fParticles->GetEntries();
    
// //     // algorithmic way of getting entry index
// //     // (primary particles are filled after secondaries)
// //     Int_t entry;
// //     if (i<fHeader.GetNprimary())
// //       entry = i+fHeader.GetNsecondary();
// //     else 
// //       entry = i-fHeader.GetNprimary();
      
// //     // only check the algorithmic way and give
// //     // the fatal error if it is wrong
// //     if (entry != fParticleFileMap[i]) {
// //       Fatal("Particle",
// //         "!!!! The algorithmic way is WRONG: !!!\n entry: %d map: %d",
// // 	entry, fParticleFileMap[i]); 
// //     }  
      
// //     fTreeK->GetEntry(fParticleFileMap[i]);
// //     new ((*fParticles)[nentries]) TParticle(*fParticleBuffer);
// //     fParticleMap->AddAt((*fParticles)[nentries],i);
// //   }
// //   return (TParticle *) (*fParticleMap)[i];

   
    
//   }


//   //scan over opened files and read corresponding TreeK##

  return ;
}
//____________________________________________________________________________ 
void AliEMCALGetter::Event(const Int_t event, const char* opt)
{
  // Reads the content of all Tree's S, D and R
  
  if (event >= gAlice->TreeE()->GetEntries() ) {
    cerr << "ERROR: AliEMCALGetter::Event -> " << event << " not found in TreeE!" << endl ; 
    return ; 
  }
  gAlice->GetEvent(event) ;

  if(strstr(opt,"H") )
   {cout<<"Reading TreeH" << endl ; 
   ReadTreeH() ;}
  
  if(strstr(opt,"S") )
    { cout << "Reading TreeS" << endl ;
    ReadTreeS(event) ;}

  if( strstr(opt,"D") )
    ReadTreeD() ;

  if( strstr(opt,"R") )
//    ReadTreeR() ;

  if( strstr(opt,"Q") )
//    ReadTreeQA() ;

  if( strstr(opt,"P") )
    ReadPrimaries() ;

}

//____________________________________________________________________________ 
const TObject * AliEMCALGetter::ReturnO(TString what, TString name, TString file) const 
{
  // get the object named "what" from the folder
  // folders are named like //Folders

  if ( file.IsNull() ) 
    file = fHeaderFile ; 

  TFolder * folder = 0 ;
  TObject * emcalO  = 0 ; 

  //  if ( name.IsNull() ) {
  if ( what.CompareTo("Hits") == 0 ) {
    folder = dynamic_cast<TFolder *>(fHitsFolder->FindObject("EMCAL")) ; 
    if (folder) 
      emcalO  = dynamic_cast<TObject *>(folder->FindObject("Hits")) ;  
  }
  else if ( what.CompareTo("SDigits") == 0 ) { 
    TString path = "EMCAL/" + file  ; 
    folder = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject(path.Data())) ; 
    if (folder) { 
    cout << "folder found" << endl ;
    if (name.IsNull())
	name = fSDigitsTitle ; 
      emcalO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    }
  }
  else if ( what.CompareTo("Digits") == 0 ){
    folder = dynamic_cast<TFolder *>(fDigitsFolder->FindObject("EMCAL")) ; 
    if (folder) { 
      if (name.IsNull())
	name = fDigitsTitle ; 
      emcalO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    } 
  }
  if (!emcalO) {
      cerr << "ERROR : AliEMCALGetter::ReturnO -> Object " << what << " not found in " << folder->GetName() << endl ; 
    return 0 ;
  }
  return emcalO ;
}
  
//____________________________________________________________________________ 
const TTask * AliEMCALGetter::ReturnT(TString what, TString name) const 
{
  // get the TTask named "what" from the folder
  // folders are named like //Folders/Tasks/what/EMCAL/name

  TString search(what) ; 

  TTask * tasks = dynamic_cast<TTask*>(fTasksFolder->FindObject(search)) ; 

  if (!tasks) {
    cerr << "ERROR: AliEMCALGetter::ReturnT -> Task " << what << " not found!" << endl ;  
    return 0 ; 
  }

  TTask * emcalT = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if (!emcalT) { 
    cerr << "ERROR: AliEMCALGetter::ReturnT -> Task " << what << "/EMCAL not found!" << endl ;  
    return 0 ; 
  }
  
  TList * list = emcalT->GetListOfTasks() ; 
 
  if (what.CompareTo("SDigitizer") == 0) {  
    if ( name.IsNull() )
      name =  fSDigitsTitle ; 
  } else  if (what.CompareTo("Digitizer") == 0){ 
    if ( name.IsNull() )
      name =  fDigitsTitle ;
  }
  
  TIter it(list) ;
  TTask * task = 0 ; 
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(name))
      return task ;
  }
  
    cout << "WARNING: AliEMCALGetter::ReturnT -> Task " << search << "/" << name << " not found!" << endl ; 
  return 0 ;
}
