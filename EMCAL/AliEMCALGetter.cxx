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
   29.05.2001 Yuri Kharlov:
              Everywhere reading the treese TTree->GetEvent(i)
              is replaced by reading the branches TBranch->GetEntry(0)
*/

//_________________________________________________________________________
//  A singleton. This class should be used in the analysis stage to get 
//  reconstructed objects: Digits, RecPoints, TrackSegments and RecParticles,
//  instead of directly reading them from galice.root file. This container 
//  ensures, that one reads Digits, made of these particular digits, RecPoints, 
//  made of these particular RecPoints, TrackSegments and RecParticles. 
//  This becomes non trivial if there are several identical branches, produced with
//  different set of parameters. 
//
//  An example of how to use (see also class AliEMCALAnalyser):
//  AliEMCALGetter * gime = AliEMCALGetter::GetInstance("galice.root","test") ;
//  for(Int_t irecp = 0; irecp < gime->NRecParticles() ; irecp++)
//     AliEMCALRecParticle * part = gime->RecParticle(1) ;
//     ................
//  please->GetEvent(event) ;    // reads new event from galice.root
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//*--         Completely redesigned by Dmitri Peressounko March 2001  
//
//*-- YS June 2001 : renamed the original AliEMCALIndexToObject and make
//*--         systematic usage of TFolders without changing the interface        
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
#include "AliEMCALHit.h"
#include "AliEMCALv1.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALClusterizer.h"
#include "AliEMCALClusterizerv1.h"
//#include "AliEMCALTrackSegmentMaker.h"
//#include "AliEMCALTrackSegmentMakerv1.h"
//#include "AliEMCALTrackSegment.h"
//#include "AliEMCALPID.h" 
//#include "AliEMCALPIDv1.h" 
#include "AliEMCALGeometry.h"

ClassImp(AliEMCALGetter)
  
  AliEMCALGetter * AliEMCALGetter::fgObjGetter = 0 ; 
  TFile * AliEMCALGetter::fFile = 0 ; 

//____________________________________________________________________________ 
AliEMCALGetter::AliEMCALGetter(const char* headerFile, const char* branchTitle, const Option_t * rw)
{
  //Initialize  all lists

  fDebug = 0 ; 
  
  fAlice = 0 ; 
 
  fHeaderFile         = headerFile ; 
  fBranchTitle        = branchTitle ;
  fSDigitsTitle       = branchTitle ; 
  fDigitsTitle        = branchTitle ; 
  fRecPointsTitle     = branchTitle ; 
  //fRecParticlesTitle  = branchTitle ; 
  //fTrackSegmentsTitle = branchTitle ; 

  fPrimaries = new TObjArray(1) ;

  fModuleFolder    = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Configuration/Modules")); 
  fPrimariesFolder = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data")); 
  fHitsFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/Hits")); 
  fSDigitsFolder   = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/SDigits")); 
  fDigitsFolder    = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Event/Data")); 
  fRecoFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Event/RecData")); 
  //fQAFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Conditions/QA")); 
  fTasksFolder     = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Tasks")) ; 

  fFailed = kFALSE ; 
 		   
  if ( fHeaderFile != "aliroot"  ) { // to call the getter without a file

    //open headers file
    fFile = static_cast<TFile*>(gROOT->GetFile(fHeaderFile.Data() ) ) ;
    if(!fFile){    //if file was not opened yet, read gAlice
      if ( fHeaderFile.Contains("_") ) {
	cerr << "AliPHOSGetter::AliPHOSGetter -> Invalid file name (_ not allowed) " << fHeaderFile.Data() << endl ;
	abort() ; 
      }
      fFile =	TFile::Open(fHeaderFile.Data(),rw) ;   
      
      if (!fFile->IsOpen()) {
	cerr << "ERROR : AliEMCALGetter::AliEMCALGetter -> Cannot open " << fHeaderFile.Data() << endl ; 
       	fFailed = kTRUE ;
        return ;  
      }
      gAlice = static_cast<AliRun *>(fFile->Get("gAlice")) ;
    }
  }
  

  if (!gAlice) {
    cerr << "ERROR : AliEMCALGetter::AliEMCALGetter -> Cannot find gAlice in " << fHeaderFile.Data() << endl ; 
    fFailed = kTRUE ;
    return ; 
  }
  if (!EMCAL()) {
    if (fDebug)
      cout << "INFO: AliEMCALGetter -> Posting EMCAL to Folders" << endl ; 
    if (gAlice->GetDetector("EMCAL")) {
      AliConfig * conf = AliConfig::Instance() ;
      conf->Add(static_cast<AliDetector*>(gAlice->GetDetector("EMCAL"))) ; 
      conf->Add(static_cast<AliModule*>(gAlice->GetDetector("EMCAL"))) ;
    }
    else 
      cerr << "ERROR: AliEMCALGetter -> detector EMCAL not found" << endl ;
  }
  
  fDebug=0;
}
//____________________________________________________________________________ 
AliEMCALGetter::~AliEMCALGetter()
{
  if (fPrimaries) {
    fPrimaries->Delete() ; 
    delete fPrimaries ; 
  }

  TFolder * emcalF = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject("EMCAL")) ;
  TCollection * folderslist = emcalF->GetListOfFolders() ; 
  TIter next(folderslist) ; 
  TFolder * folder = 0 ; 
  while ( (folder = static_cast<TFolder*>(next())) ) 
    emcalF->Remove(folder) ; 

  delete fFile ; 
  fFile = 0 ;

  fgObjGetter = 0 ; 

}

//____________________________________________________________________________ 
void AliEMCALGetter::CreateWhiteBoard() const
{

}

//____________________________________________________________________________ 
void AliEMCALGetter::CloseFile()
{
  delete gAlice ; 
  gAlice = 0 ; 
  delete fAlice ; 
  fAlice = 0 ; 
}

//____________________________________________________________________________ 
const TFolder * AliEMCALGetter::Folder(const TString what) const {

  // returns the EMCAL folder required by what
  // what = hits, sdigits, digits

  if ( what == "hits" ) 
    return dynamic_cast<const TFolder *>(fHitsFolder->FindObject("EMCAL")) ; 
  else if ( what == "sdigits" ) 
    return  dynamic_cast<const TFolder *>(fSDigitsFolder->FindObject("EMCAL")) ; 
  else if ( what == "digits" ) 
    return  dynamic_cast<const TFolder *>(fDigitsFolder->FindObject("EMCAL")) ; 
  else {
    cerr << "ERROR: AliEMCALGetter::GetFolder -> " << what.Data() << " illegal option (hits, sdigits, digits) " << endl ; 
    return 0 ; 
  }
}
//____________________________________________________________________________ 
AliEMCALGetter * AliEMCALGetter::GetInstance()
{
  // Returns the pointer of the unique instance already defined
  
  if ( fgObjGetter ) {
    fFile->cd() ; 
    return fgObjGetter ;
  }
  else {
    // cout << "AliEMCALGetter::GetInstance ERROR: not yet initialized" << endl ;
    return 0 ;
  }
}

//____________________________________________________________________________ 
AliEMCALGetter * AliEMCALGetter::GetInstance(const char* headerFile,
					   const char* branchTitle, const Option_t * rw)
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed 

  if ( fgObjGetter && fFile->IsOpen()) // an instance exists and the file is still open        
    if((fgObjGetter->fBranchTitle.CompareTo(branchTitle) == 0) && 
       (fgObjGetter->fHeaderFile.CompareTo(headerFile)==0)) {
      fFile->cd() ; 
      return fgObjGetter ;
    }
    else // another file than the existing one is required, scratch the getter
      fgObjGetter->~AliEMCALGetter() ;  // delete it already exists another version

  fgObjGetter = new AliEMCALGetter(headerFile,branchTitle, rw) ; 

  if (fgObjGetter->HasFailed() ) 
    fgObjGetter = 0 ; 
    
  fFile->cd() ; 
  return fgObjGetter ; 
  
}

//____________________________________________________________________________ 
const AliEMCALv1 * AliEMCALGetter::EMCAL() 
{
  // returns the EMCAL object 
  AliEMCALv1 * emcal = dynamic_cast<AliEMCALv1*>(fModuleFolder->FindObject("EMCAL")) ;  
  if (!emcal) 
    if (fDebug)
      cout << "WARNING: AliEMCALGetter::EMCAL -> EMCAL module not found in Folders" << endl ; 
  return emcal ; 
}  

//____________________________________________________________________________ 
AliEMCALGeometry * AliEMCALGetter::EMCALGeometry() 
{
  AliEMCALGeometry * rv = 0 ; 
  if (EMCAL() )
    rv =  EMCAL()->GetGeometry() ;
  return rv ; 
} 

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostPrimaries(void) const 
{  //------- Primaries ----------------------

  // the hierarchy is //Folders/RunMC/Event/Data/Primaries
  
  TFolder * primariesFolder = dynamic_cast<TFolder*>(fPrimariesFolder->FindObject("Primaries")) ; 
  if ( !primariesFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post Primaries -> Folder //" << fPrimariesFolder->GetName() << "/Primaries/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post Primaries -> Adding Folder //" << fPrimariesFolder->GetName() << "/Primaries/"  << endl;
    }
    primariesFolder = fPrimariesFolder->AddFolder("Primaries", "Primaries particles from TreeK") ; 
  }    
  TClonesArray *primaries=  new TClonesArray("TParticle",1000) ;
  primaries->SetName("Primaries") ;
  primariesFolder->Add(primaries) ; 
  
  return kTRUE;
} 

//____________________________________________________________________________ 
TObject** AliEMCALGetter::PrimariesRef(void) const 
{  //------- Primaries ----------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/Primaries
  if ( !fPrimariesFolder ) {
    cerr << "ERROR: AliEMCALGetter::PrimariesRef -> Folder //" << fPrimariesFolder << " not found!" << endl;
    abort() ;
  }    
 
  TFolder * primariesFolder = dynamic_cast<TFolder *>(fPrimariesFolder->FindObject("Primaries")) ;
  if ( !primariesFolder ) {
    cerr << "ERROR: AliEMCALGetter::PrimariesRef -> Folder //" << fPrimariesFolder << "/Primaries/ not found!" << endl;  
    abort() ;
  }
 
  TObject * p = primariesFolder->FindObject("Primaries") ;
  if(!p) {
    cerr << "ERROR: AliEMCALGetter::PrimariesRef -> " << primariesFolder->GetName() << "/Primaries not found !" << endl ; 
    abort() ;
  }
  else
    return primariesFolder->GetListOfFolders()->GetObjectRef(p) ;
}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostHits(void) const 
{  //------- Hits ----------------------

  // the hierarchy is //Folders/RunMC/Event/Data/EMCAL/Hits
  
  TFolder * emcalFolder = dynamic_cast<TFolder*>(fHitsFolder->FindObject("EMCAL")) ; 
  if ( !emcalFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post H -> Folder //" << fHitsFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post H -> Adding Folder //" << fHitsFolder << "/EMCAL/"  << endl;
    }
    emcalFolder = fHitsFolder->AddFolder("EMCAL", "Hits from EMCAL") ; 
  }    
  TClonesArray *hits=  new TClonesArray("AliEMCALHit",1000) ;
  hits->SetName("Hits") ;
  emcalFolder->Add(hits) ; 
  
  return kTRUE;
} 

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::HitsRef(void) const 
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
    return emcalFolder->GetListOfFolders()->GetObjectRef(h) ;
}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostSDigits(const char * name, const char * headerFile) const 
{  //---------- SDigits -------------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/EMCAL/SDigits/headerFile/sdigitsname
  // because you can have sdigits from several hit files for mixing
  
  TFolder * emcalFolder = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("EMCAL")) ;
  if ( !emcalFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post S -> Folder //" << fSDigitsFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post S -> Adding Folder //" << fHitsFolder << "/EMCAL/" << endl;
    }
    emcalFolder = fSDigitsFolder->AddFolder("EMCAL", "SDigits from EMCAL") ; 
  }    
  TString subdir(headerFile) ;
  subdir.ReplaceAll("/", "_") ; 
  TFolder * emcalSubFolder = dynamic_cast<TFolder*>(emcalFolder->FindObject(subdir)) ; 
  if ( !emcalSubFolder ) 
    emcalSubFolder = emcalFolder->AddFolder(subdir, ""); 
  
  TObject * sd  = emcalSubFolder->FindObject(name); 
  if ( sd ) {
    if (fDebug)
      cerr <<"INFO: AliEMCALGetter::Post S -> Folder " << subdir 
	   << " already exists!" << endl ;  
  }else{
    TClonesArray * sdigits = new TClonesArray("AliEMCALDigit",1) ;
    sdigits->SetName(name) ;
    emcalSubFolder->Add(sdigits) ;
  }
  
  return kTRUE;
} 
//____________________________________________________________________________ 
TObject ** AliEMCALGetter::SDigitsRef(const char * name, const char * file) const 
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
    return emcalSubFolder->GetListOfFolders()->GetObjectRef(dis) ;

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
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Ser ->//" << fTasksFolder << "/SDigitizer/EMCAL/ not found!" << endl;  
      cout <<"INFO: AliEMCALGetter::Post Ser -> Adding //" << fTasksFolder << "/SDigitizer/EMCAL/" << endl;
    }
    emcal = new TTask("EMCAL", "") ; 
    sd->Add(emcal) ; 
  } 
  AliEMCALSDigitizer * emcalsd  = dynamic_cast<AliEMCALSDigitizer *>(emcal->GetListOfTasks()->FindObject( sdigitizer->GetName() )); 
  if (emcalsd) { 
    if (fDebug)
      cout << "INFO: AliEMCALGetter::Post Ser -> Task " << sdigitizer->GetName() << " already exists" << endl ; 
    emcal->GetListOfTasks()->Remove(emcalsd) ;
  }
  emcal->Add(sdigitizer) ;	
  return kTRUE; 
  
}

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::SDigitizerRef(const char * name) const 
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

  return emcal->GetListOfTasks()->GetObjectRef(task) ;

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
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Ser ->  //" << fTasksFolder << "/SDigitizer/EMCAL/ not found!" << endl;
      cout <<"INFO: AliEMCALGetter::Post Ser -> Adding  //" << fTasksFolder << "/SDigitizer/EMCAL" << endl;
    }
    emcal = new TTask("EMCAL", "") ; 
    sd->Add(emcal) ; 
  } 

  TString sdname(name) ;
  sdname.Append(":") ;
  sdname.Append(file);
  sdname.ReplaceAll("/","_") ; 
  AliEMCALSDigitizer * emcalsd  = dynamic_cast<AliEMCALSDigitizer *>(emcal->GetListOfTasks()->FindObject( sdname )); 
  if (!emcalsd) {
    emcalsd = new AliEMCALSDigitizer() ;  
    //Note, we can not call constructor with parameters: it will call Getter and screw up everething
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
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post D -> Folder //" << fDigitsFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post D -> Adding Folder //" << fDigitsFolder << "/EMCAL/" << endl;
    }
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
TObject ** AliEMCALGetter::DigitsRef(const char * name) const 
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
    return emcalFolder->GetListOfFolders()->GetObjectRef(d) ;

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
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Der ->  //" << fTasksFolder << "/Digitizer/EMCAL not found!" << endl;
      cout <<"INFO: AliEMCALGetter::Post Der -> Adding //" << fTasksFolder << "/Digitizer/EMCAL" << endl; 
    }
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
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Der -> //" << fTasksFolder << "/Digitizer/EMCAL not found!" << endl; 
      cout <<"INFO: AliEMCALGetter::Post Der -> Adding //" << fTasksFolder << "/Digitizer/EMCAL" << endl;
    }
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
TObject ** AliEMCALGetter::DigitizerRef(const char * name) const 
{  
  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 
  if ( !sd ) {
    cerr << "ERROR: AliEMCALGetter::Post DerRef -> Task //" << fTasksFolder->GetName() << "/Digitizer not found!" << endl;
    abort();
  }        

  TTask * emcal = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    cerr <<"ERROR: AliEMCALGetter::Post DerRef ->  //" << fTasksFolder->GetName() << "/Digitizer/EMCAL" << endl;
    abort();
  }        

  TTask * task = dynamic_cast<TTask*>(emcal->GetListOfTasks()->FindObject(name)) ; 

  return emcal->GetListOfTasks()->GetObjectRef(task) ;

}
 
//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostRecPoints(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/TowerRecPoints/name
  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/PreShoRecPoints/name

  TFolder * emcalFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL")) ; 
  
  if ( !emcalFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post RPo -> Folder //" << fRecoFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post Rpo -> Adding Folder //" << fRecoFolder << "/EMCAL/" << endl;
    }
    emcalFolder = fRecoFolder->AddFolder("EMCAL", "Reconstructed data from EMCAL") ;  
  }    
  
  // Tower RecPoints 
  TFolder * emcalRPoTowerFolder  = dynamic_cast<TFolder*>(emcalFolder->FindObject("TowerRecPoints")) ;
  if ( !emcalRPoTowerFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post RPo -> Folder //" << fRecoFolder << "/EMCAL/TowerRecPoints/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post Rpo -> Adding Folder //" << fRecoFolder << "/EMCAL/TowerRecPoints not found!" << endl;
    }
    emcalRPoTowerFolder = emcalFolder->AddFolder("TowerRecPoints", "Tower RecPoints from EMCAL") ;  
  }    
  
  TObject * erp = emcalFolder->FindObject( name ) ;
  if ( !erp )   {
    TObjArray * towerrp = new TObjArray(100) ;
    towerrp->SetName(name) ;
    emcalRPoTowerFolder->Add(towerrp) ;  
  }

  // Pre Shower RecPoints 
  TFolder * emcalRPoPreShoFolder  = dynamic_cast<TFolder*>(emcalFolder->FindObject("PreShoRecPoints")) ;
  if ( !emcalRPoPreShoFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post RPo -> Folder //" << fRecoFolder << "/EMCAL/PreShoRecPoints/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post Rpo -> Adding Folder //" << fRecoFolder << "/EMCAL/PreShoRecPoints/" << endl;
    }
    emcalRPoPreShoFolder = emcalFolder->AddFolder("PreShoRecPoints", "PreSho RecPoints from EMCAL") ;  
  }    
  
  TObject * crp =  emcalRPoPreShoFolder->FindObject( name ) ;
  if ( !crp )   {
    TObjArray * preshorp = new TObjArray(100) ;
    preshorp->SetName(name) ;
    emcalRPoPreShoFolder->Add(preshorp) ;  
  }
  return kTRUE; 
}

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::TowerRecPointsRef(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/TowerRecPoints/name
   
  if ( !fRecoFolder ) {
    cerr << "ERROR: AliEMCALGetter::TowerRecPointsRef -> Folder //" << fRecoFolder << " not found!" << endl;
    abort() ; 
  }    

  TFolder * towerFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/TowerRecPoints")) ; 
  if ( !towerFolder ) {
    cerr << "ERROR: AliEMCALGetter::TowerRecPointsRef -> Folder //" << fRecoFolder << "/EMCAL/TowerRecPoints/ not found!" << endl;
    abort() ;
  }    


  TObject * trp = towerFolder->FindObject(name ) ;
  if ( !trp )   {
    cerr << "ERROR: AliEMCALGetter::TowerRecPointsRef -> Object " << name << " not found!" << endl  ;
    abort() ; 
  }
  return towerFolder->GetListOfFolders()->GetObjectRef(trp) ;

} 

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::PreShoRecPointsRef(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/PreShoRecPoints/name
   
  if ( !fRecoFolder ) {
    cerr << "ERROR: AliEMCALGetter::PreShoRecPointsRef -> Folder //" << fRecoFolder << " not found!" << endl;
    abort() ; 
  }    

  TFolder * preshoFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/PreShoRecPoints")) ; 
  if ( !preshoFolder ) {
    cerr << "ERROR: AliEMCALGetter::PreShoRecPointsRef -> Folder //" << fRecoFolder << "/EMCAL/PreShoRecPoints/" << endl;
    abort() ;
  }    

  TObject * prp = preshoFolder->FindObject(name ) ;
  if ( !prp )   {
    cerr << "ERROR: AliEMCALGetter::PreShoRecPointsRef -> Object " << name << " not found! " << endl ; 
    abort() ;
  }
  return preshoFolder->GetListOfFolders()->GetObjectRef(prp) ;

} 

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostClusterizer(AliEMCALClusterizerv1 * clu) const 
{ // ------------------ AliEMCALClusterizer ------------------------
  
  // the hierarchy is //Folders/Tasks/Reconstructioner/EMCAL/sdigitsname

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliEMCALGetter::Post Rer -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
        
  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Rer -> //" << fTasksFolder << "/ReconstructionerEMCAL not found!" << endl; 
      cout <<"INFO: AliEMCALGetter::Post Rer -> Adding //" << fTasksFolder << "/Reconstructioner/EMCAL" << endl; 
    }
    emcal = new TTask("EMCAL", "") ; 
    tasks->Add(emcal) ; 
  } 

  AliEMCALClusterizerv1 * emcalcl = dynamic_cast<AliEMCALClusterizerv1*>(emcal->GetListOfTasks()->FindObject(clu->GetName())) ; 
  if (emcalcl) { 
    if (fDebug)
      cout << "INFO: AliEMCALGetter::Post Rer -> Task " << clu->GetName() << " already exists" << endl ; 
    emcalcl->Delete() ; 
    emcal->GetListOfTasks()->Remove(emcalcl) ;
  }
  emcal->Add(clu) ;      
  return kTRUE; 
} 

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::ClusterizerRef(const char * name) const 
{ // ------------------ AliEMCALClusterizer ------------------------
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliEMCALGetter::Post RerRef -> Task //" << fTasksFolder->GetName() << "/Reconstructioner not found!" << endl;
    abort() ;
  }        
        
  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    cerr <<"WARNING: AliEMCALGetter::Post RerRef -> //" << fTasksFolder->GetName() << "/Reconstructioner/EMCAL" << endl; 
    abort() ; 
  }   

  TList * l = emcal->GetListOfTasks() ; 
  TIter it(l) ;
  TTask * task ;
  TTask * clu = 0 ;
  TString cluname(name) ;
  cluname+=":clu-" ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(cluname)){
      clu = task ;
      break ;
    }
  }

  if(clu) 
    return l->GetObjectRef(clu) ;
  else {
    cerr << "ERROR: AliEMCALGetter::Post RerRef -> task " << task->GetName() << " not found! " << endl ; 
    abort() ;
  }
}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostClusterizer(const char * name) const 
{ // ------------------ AliEMCALClusterizer ------------------------

  // the hierarchy is //Folders/Tasks/Reconstructioner/EMCAL/sdigitsname
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliEMCALGetter::Post Rer -> Task//" << fTasksFolder->GetName() << "/Reconstructioner not found!" << endl; 
    return kFALSE ;
  }        
  
  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Rer -> //" << fTasksFolder << "/Reconstructioner/EMCAL not found!" << endl;
      cout <<"INFO: AliEMCALGetter::Post Rer -> Adding //" << fTasksFolder << "/Reconstructioner/EMCAL" << endl;
    }
    emcal = new TTask("EMCAL", "") ; 
    tasks->Add(emcal) ; 
  } 

  TList * l = emcal->GetListOfTasks() ;   
  TIter it(l) ;
  TString clun(name) ;
  clun+=":clu" ; 
  TTask * task ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(clun))
      return kTRUE ;
  }

  AliEMCALClusterizerv1 * emcalcl = new AliEMCALClusterizerv1() ;
  emcal->Add(emcalcl) ;
  return kTRUE; 
  
}

//____________________________________________________________________________ 
/*Bool_t AliEMCALGetter::PostTrackSegments(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/TrackSegments/name

  TFolder * emcalFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL")) ; 
  
  if ( !emcalFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post TS -> Folder //" << fRecoFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post TS -> Adding Folder //" << fRecoFolder << "/EMCAL" << endl;
    }
    emcalFolder = fRecoFolder->AddFolder("EMCAL", "Reconstructed data from EMCAL") ;  
  }    

  TFolder * emcalTSFolder  = dynamic_cast<TFolder*>(emcalFolder->FindObject("TrackSegments")) ;
  if ( !emcalTSFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post TS -> Folder//" << fRecoFolder << "/EMCAL/TrackSegments/ not found!" << endl; 
      cout << "INFO:    AliEMCALGetter::Post TS -> Adding Folder //" << fRecoFolder << "/EMCAL/TrackSegments/" << endl; 
    }
    emcalTSFolder = emcalFolder->AddFolder("TrackSegments", "TrackSegments from EMCAL") ;  
  }    
  
  TObject * tss =  emcalTSFolder->FindObject( name ) ;
  if (!tss) {
    TClonesArray * ts = new TClonesArray("AliEMCALTrackSegment",100) ;
    ts->SetName(name) ;
    emcalTSFolder->Add(ts) ;  
  }
  return kTRUE; 
} 

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::TrackSegmentsRef(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/TrackSegments/name

 if ( !fRecoFolder ) {
    cerr << "ERROR: AliEMCALGetter::TrackSegmentsRef -> Folder //" << fRecoFolder << "not found!" << endl;
    abort() ; 
  }    

  TFolder * emcalFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/TrackSegments")) ; 
  if ( !emcalFolder ) {
    cerr << "ERROR: AliEMCALGetter::TrackSegmentsRef -> Folder //" << fRecoFolder << "/EMCAL/TrackSegments/ not found!" << endl;
    abort();
  }    
  
  TObject * tss =  emcalFolder->FindObject(name) ;
  if (!tss) {
    cerr << "ERROR: AliEMCALGetter::TrackSegmentsRef -> object " << name << " not found! " << endl ;  
    abort() ;  
  }
  return emcalFolder->GetListOfFolders()->GetObjectRef(tss) ;
} 

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostTrackSegmentMaker(AliEMCALTrackSegmentMaker * tsmaker) const 
{ //------------Track Segment Maker ------------------------------
  
  // the hierarchy is //Folders/Tasks/Reconstructioner/EMCAL/sdigitsname

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliEMCALGetter::Post Ter -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
        
  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Rer -> //" << fTasksFolder << "/Reconstructioner/EMCAL not found!" << endl; 
      cout <<"INFO: AliEMCALGetter::Post Rer -> Adding //" << fTasksFolder << "/Reconstructioner/EMCAL" << endl;
    }
    emcal = new TTask("EMCAL", "") ; 
    tasks->Add(emcal) ; 
  } 

  AliEMCALTrackSegmentMaker * emcalts = 
    dynamic_cast<AliEMCALTrackSegmentMaker*>(emcal->GetListOfTasks()->FindObject(tsmaker->GetName())) ; 
  if (emcalts) { 
    emcalts->Delete() ;
    emcal->GetListOfTasks()->Remove(emcalts) ;
  }
  emcal->Add(tsmaker) ;      
  return kTRUE; 
  
} 
//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostTrackSegmentMaker(const char * name) const 
{ //------------Track Segment Maker ------------------------------
  
  // the hierarchy is //Folders/Tasks/Reconstructioner/EMCAL/sdigitsname
  
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 
  
  if ( !tasks ) {
    cerr << "ERROR: AliEMCALGetter::Post Ter -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
  
  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Rer -> //" << fTasksFolder << "/Reconstructioner/EMCAL not found!" << endl; 
      cout <<"INFO: AliEMCALGetter::Post Rer -> Adding //" << fTasksFolder->GetName() << "/Reconstructioner/EMCAL" << endl;
    }
    emcal = new TTask("EMCAL", "") ; 
    tasks->Add(emcal) ; 
  } 

  AliEMCALTrackSegmentMaker * emcalts = 
    dynamic_cast<AliEMCALTrackSegmentMaker*>(emcal->GetListOfTasks()->FindObject(name)) ; 
  if (!emcalts) { 
    emcalts = new AliEMCALTrackSegmentMakerv1() ;
    TString tsn(name);
    tsn+=":tsm-v1" ;
    emcalts->SetName(tsn) ;
    emcal->Add(emcalts) ;      
  }
  return kTRUE; 
  
} 

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::TSMakerRef(const char * name) const 
{ //------------Track Segment Maker ------------------------------
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliEMCALGetter::TSLakerRef TerRef -> Task //" << fTasksFolder->GetName() << "/Reconstructioner not found!" << endl;
    abort() ;
  }        
        
  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    cerr <<"WARNING: AliEMCALGetter::TSMakerRef TerRef -> //" << fTasksFolder->GetName() << "/Reconstructioner/EMCAL not found!" << endl; 
    abort() ; 
  }   

  TList * l = emcal->GetListOfTasks() ; 
  TIter it(l) ;
  TTask * task ;
  TTask * tsm = 0 ;
  TString tsmname(name) ;
  tsmname+=":tsm-" ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(tsmname)){
      tsm = task ;
      break ;
    }
  }
  
  if(tsm) 
    return l->GetObjectRef(tsm) ;
  else {
    cerr << "ERROR: AliEMCALGetter::TSLakerRef -> task " << task->GetName() << " not found! " << endl ; 
    abort() ;
  }
} 

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostRecParticles(const char * name) const 
{  // -------------------- RecParticles ------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/TrackSegments/name

  TFolder * emcalFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL")) ; 
  
  if ( !emcalFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post RPa -> Folder //" << fRecoFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post Rpa -> Adding Folder //" << fRecoFolder << "/EMCAL/" << endl;
    }
    emcalFolder = fRecoFolder->AddFolder("EMCAL", "Reconstructed data from EMCAL") ;  
  }    

 TFolder * emcalRPaFolder  = dynamic_cast<TFolder*>(emcalFolder->FindObject("RecParticles")) ;
  if ( !emcalRPaFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post RPa -> Folder //" << fRecoFolder << "/EMCAL/RecParticles/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post RPa -> Adding Folder //" << fRecoFolder << "/EMCAL/RecParticles/" << endl;
    }
    emcalRPaFolder = emcalFolder->AddFolder("RecParticles", "RecParticles from EMCAL") ;  
  } 

  TObject * rps = emcalRPaFolder->FindObject( name )  ;
  if ( !rps ) {
    TClonesArray * rp = new TClonesArray("AliEMCALRecParticle",100) ;
    rp->SetName(name) ;    
    emcalRPaFolder->Add(rp) ;  
  }
  return kTRUE; 
} 

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::RecParticlesRef(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/TrackSegments/name

 if ( !fRecoFolder ) {
    cerr << "ERROR: AliEMCALGetter::RecParticlesRef -> Folder//" << fRecoFolder << " not found!" << endl; 
    abort() ; 
  }    

  TFolder * emcalFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/RecParticles")) ; 
  if ( !emcalFolder ) {
    cerr << "ERROR: AliEMCALGetter::RecParticlesRef -> Folder //" << fRecoFolder << "/EMCAL/RecParticles/ not found!" << endl;
    abort() ;
  }    

  TObject * tss =  emcalFolder->FindObject(name) ;
  if (!tss) {
    cerr << "ERROR: AliEMCALGetter::RecParticlesRef -> object " << name << " not found! " << endl ; 
    abort() ;  
  }
  return emcalFolder->GetListOfFolders()->GetObjectRef(tss) ;
}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostPID(AliEMCALPID * pid) const 
{      // ------------AliEMCAL PID -----------------------------

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliEMCALGetter::Post Per -> Task //" << fTasksFolder->GetName() << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
  
  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Per -> //" << fTasksFolder->GetName() << "/Reconstructioner/EMCAL not found!" << endl; 
      cout <<"INFO: AliEMCALGetter::Post Per -> Adding //" << fTasksFolder->GetName() << "/Reconstructioner/EMCAL" << endl;
    }
    emcal = new TTask("EMCAL", "") ; 
    tasks->Add(emcal) ; 
  } 

  AliEMCALPID * emcalpid = dynamic_cast<AliEMCALPID*>(emcal->GetListOfTasks()->FindObject(pid->GetName())) ; 
  if (emcalpid) { 
    if (fDebug)
      cout << "INFO: AliEMCALGetter::Post Per -> Task " << pid->GetName()
	   << " already exists" << endl ; 
    emcal->GetListOfTasks()->Remove(emcalpid) ;
  }
  
  emcal->Add(pid) ;      
  return kTRUE; 
} 

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostPID(const char * name) const 
{     
  // the hierarchy is //Folders/Tasks/Reconstructioner/EMCAL/sdigitsname
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliEMCALGetter::Post Per -> Task //" << fTasksFolder->GetName() << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
  
  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    if (fDebug) {
      cout <<"WARNING: AliEMCALGetter::Post Per -> //" << fTasksFolder->GetName() << "/Reconstructioner/EMCAL not found!" << endl; 
      cout <<"INFO: AliEMCALGetter::Post Per -> Adding //" << fTasksFolder->GetName() << "/Reconstructioner/EMCAL" << endl;
    }
    emcal = new TTask("EMCAL", "") ; 
    tasks->Add(emcal) ; 
  } 

  TList * l = emcal->GetListOfTasks() ;   
  TIter it(l) ;
  TString pidname(name) ;
  pidname+=":pid" ; 
  TTask * task ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(pidname))
      return kTRUE ;
  }
 
  AliEMCALPIDv1 * emcalpid = new AliEMCALPIDv1() ;
  pidname+="-v1" ;
  emcalpid->SetName(pidname) ;
  emcal->Add(emcalpid) ;      
  
  return kTRUE; 
} 

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::PIDRef(const char * name) const 
{ //------------PID ------------------------------

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliEMCALGetter::PIDRef PerRef -> Task //" << fTasksFolder->GetName() << "/Reconstructioner not found!" << endl;
    abort() ;
  }        
        
  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 
  if ( !emcal )  {
    cerr <<"WARNING: AliEMCALGetter::PIDRef PerRef -> //" << fTasksFolder->GetName() << "/ReconstructionerEMCAL not found!" << endl; 
    abort() ; 
  }   
  
  TList * l = emcal->GetListOfTasks() ; 
  TIter it(l) ;
  TTask * task ;
  TTask * pid = 0 ;
  TString pidname(name) ;
  pidname+=":pid-" ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(pidname)){
      pid = task ;
      break ;
    }
  }
  
  if(pid) 
    return l->GetObjectRef(pid) ;
  else {
    cerr << "ERROR: AliEMCALGetter::PIDRef -> task " << task->GetName() << " not found! " << endl ;  
    abort() ;
  }
} 

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::PostQA(void) const 
{ // ------------------ QA ---------------------------------

  // the hierarchy is //Folders/Run/Conditions/QA/EMCAL/alarmsName

  TFolder * emcalFolder = dynamic_cast<TFolder*>(fQAFolder->FindObject("EMCAL")) ; 
  if ( !emcalFolder ) {
    if (fDebug) {
      cout << "WARNING: AliEMCALGetter::Post Q -> Folder //" << fQAFolder << "/EMCAL/ not found!" << endl;
      cout << "INFO:    AliEMCALGetter::Post Q -> Adding Folder //" << fQAFolder << "/EMCAL/" << endl;
    }
    emcalFolder = fQAFolder->AddFolder("EMCAL", "QA from EMCAL") ; 
  }      

  return kTRUE;
}

//____________________________________________________________________________ 
TObject ** AliEMCALGetter::AlarmsRef(void) const 
{  //------- Alarms ----------------------

  
  // the hierarchy is //Folders/Run/Conditions/QA/EMCAL
  if ( !fQAFolder ) {
    cerr << "ERROR: AliEMCALGetter::AlarmsRef QRef -> Folder //" << fQAFolder << " not found!" << endl;
    abort() ;
  }    
 
  TFolder * emcalFolder = dynamic_cast<TFolder *>(fQAFolder->FindObject("EMCAL")) ;
  if ( !emcalFolder ) {
    cerr << "ERROR: AliEMCALGetter::AlarmsRef QRef -> Folder //" << fQAFolder << "/EMCAL/ not found!" << endl;
    abort() ;
  }
   
  return fQAFolder->GetListOfFolders()->GetObjectRef(emcalFolder) ;
}
*/
//____________________________________________________________________________ 
const TParticle * AliEMCALGetter::Primary(Int_t index) const
{
  // Return primary particle numbered by <index>
  
  if(index < 0) 
    return 0 ;
  TParticle *  p = 0 ;
  if (fAlice) 
    p = fAlice->Particle(index) ; 
  else 
    p = gAlice->Particle(index) ; 
  //   if (p->GetFirstMother() != -1 ) {
  //     cout << "AliEMCALGetter::Primary : Not a primary " << endl ; 
  //   }
  
  return p ; 
      
}

//____________________________________________________________________________ 
void AliEMCALGetter::ReadTreeD()
{
  // Read the digit tree gAlice->TreeD()  

  TTree * treeD = gAlice->TreeD() ;
  
  if(!treeD) { // TreeD not found in header file
    
    if (fDebug) 
      cout <<   "WARNING: AliEMCALGetter::ReadTreeD -> Cannot find TreeD in " << fHeaderFile << endl ;
    
    TString searchFileName("") ; 
    
    if (Digitizer())  // Digitizer found in header file
      searchFileName = Digitizer()->GetTitle() ; 
    
    else if (Clusterizer())  // Clusterizer found in header file
      searchFileName = Clusterizer()->GetDigitsFileName() ; 
    
    if (treeD = TreeD(searchFileName)) { //found TreeD in the file which contains the hits
      if (fDebug) 
	cout << "INFO: AliEMCALGetter::ReadTreeD -> TreeD found in " << searchFileName.Data() << endl ; 
      
    } else {
      cerr << "ERROR: AliEMCALGetter::ReadTreeD -> TreeD not found " << endl ; 
      return ;
    }   
  }
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
    if (fDebug)
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

  TTree * treeH = gAlice->TreeH() ;
  
  if(!treeH) {// TreeH not found in header file
    
    if (fDebug) 
      cout <<   "WARNING: AliEMCALGetter::ReadTreeH -> Cannot find TreeH in " << fHeaderFile << endl ;
    
    TString searchFileName("") ; 
    
    if (SDigitizer())  // SDigitizer found in header file
      searchFileName = SDigitizer()->GetTitle() ;
    
    else if (Digitizer())  // Digitizer found in header file
      searchFileName = Digitizer()->GetHitsFileName() ; 
    
    else if (Clusterizer())  // Clusterizer found in header file
      searchFileName = Clusterizer()->GetHitsFileName() ; 
    
    if (treeH = TreeH(searchFileName)) { //found TreeH in the file which contains the hits
      if (fDebug) 
	cout << "INFO: AliEMCALGetter::ReadTreeH -> TreeH found in " << searchFileName.Data() << endl ; 
      
    } else {
      cerr << "ERROR: AliEMCALGetter::ReadTreeH -> TreeH not found " << endl ; 
      return ;
    }  
  }  
  
  TBranch * hitsbranch = static_cast<TBranch*>(gAlice->TreeH()->GetBranch("EMCAL")) ;
  if ( !hitsbranch ) {
    if (fDebug)
      cout << "WARNING:  AliEMCALGetter::ReadTreeH -> Cannot find branch EMCAL" << endl ; 
    return ;
  }
  if(!Hits())
    PostHits() ;
 
  if (hitsbranch->GetEntries() > 1 ) {
    TClonesArray * tempo =  new TClonesArray("AliEMCALHit",1000) ;
    TClonesArray * hits = dynamic_cast<TClonesArray*>(*HitsRef()) ; 
    hitsbranch->SetAddress(&tempo) ;
    Int_t index = 0 ; 
    Int_t i = 0 ;
    for (i = 0 ; i < hitsbranch->GetEntries() ; i++) {
      hitsbranch->GetEntry(i) ;
      Int_t j = 0 ; 
      for ( j = 0 ; j < tempo->GetEntries() ; j++) { 
	const AliEMCALHit * hit = static_cast<const AliEMCALHit *>(tempo->At(j)) ; 
	new((*hits)[index]) AliEMCALHit( *hit ) ;
	index++ ; 
      }
    }
    delete tempo ; 
  }
  else {
  hitsbranch->SetAddress(HitsRef()) ;
  hitsbranch->GetEntry(0) ;
  }
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
    if (fDebug)
      cout << "WARNING:  AliEMCALGetter::ReadTreeH -> Cannot find branch EMCAL" << endl ; 
    return ;
  }  
  if(!Hits())
    PostHits() ;
  hitsbranch->SetAddress(HitsRef()) ;
  hitsbranch->GetEntry(itrack) ;

  
}


//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeK(TString filename)  
{

  // returns TreeK from file filename
  // usefull in case of split file

  if ( filename.IsNull() ) 
    filename = fHeaderFile ; 

  TFile * file = 0 ; 
  // file = static_cast<TFile*>(gROOT->GetFile(filename.Data() ) ) ;
  if (!file) {  // file not open yet
    //   file->Close() ; 
    file = TFile::Open(filename.Data(), "read") ; 
    delete fAlice ; 
    fAlice = static_cast<AliRun *>(file->Get("gAlice")) ; 
  }

  TString treeName("TreeK") ; 
  treeName += EventNumber()  ; 
  TTree * tree = static_cast<TTree *>(file->Get(treeName.Data())) ;
  if (!tree && fDebug)  
    cout << "WARNING: AliEMCALGetter::TreeK -> " << treeName.Data() << " not found in " << filename.Data() << endl ; 
  
  return tree ; 		      
}

//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeH(TString filename)  
{

  // returns TreeH from file filename
  // usefull in case of split file

  if ( filename.IsNull() ) 
    filename = fHeaderFile ; 

  TFile * file = 0 ; 
  file = static_cast<TFile*>(gROOT->GetFile(filename.Data() ) ) ;
  if (!file) { // file not open yet
    file = TFile::Open(filename.Data(), "read") ; 
  }
  TString treeName("TreeH") ; 
  treeName += EventNumber()  ; 
  TTree * tree = static_cast<TTree *>(file->Get(treeName.Data())) ;
  if (!tree && fDebug)  
    cout << "WARNING: AliEMCALGetter::TreeH -> " << treeName.Data() << " not found in " << filename.Data() << endl ; 
  
  return tree ; 		      
}

//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeS(TString filename)  
{

  // returns TreeS from file filename
  // usefull in case of split file

  if ( filename.IsNull() ) 
    filename = fHeaderFile ; 

  TFile * file = 0 ; 
  file = static_cast<TFile*>(gROOT->GetFile(filename.Data() ) ) ;
  if (!file) { // file not open yet
    file = TFile::Open(filename.Data(), "read") ; 
  }
  TString treeName("TreeS") ; 
  treeName += EventNumber()  ; 
  TTree * tree = static_cast<TTree *>(file->Get(treeName.Data())) ;
  if (!tree && fDebug)  
    cout << "WARNING: AliEMCALGetter::TreeS -> " << treeName.Data() << " not found in " << filename.Data() << endl ; 
  
  return tree ; 		      
}

//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeD(TString filename)  
{

  // returns TreeD from file filename
  // usefull in case of split file

  if ( filename.IsNull() ) 
    filename = fHeaderFile ; 

  TFile * file = 0 ; 
  file = static_cast<TFile*>(gROOT->GetFile(filename.Data() ) ) ;
  if (!file) { // file not open yet
    file = TFile::Open(filename.Data(), "read") ; 
  }
  TString treeName("TreeD") ; 
  treeName += EventNumber()  ; 
  TTree * tree = static_cast<TTree *>(file->Get(treeName.Data())) ;
  if (!tree && fDebug)  
    cout << "WARNING: AliEMCALGetter::TreeD -> " << treeName.Data() << " not found in " << filename.Data() << endl ; 
  
  return tree ; 		      
}

//____________________________________________________________________________ 
// void AliEMCALGetter::ReadTreeQA()
//{
  // Read the digit tree gAlice->TreeQA()
  // so far only EMCAL knows about this Tree  

//   if(EMCAL()->TreeQA()== 0){
//     cerr <<   "ERROR: AliEMCALGetter::ReadTreeQA: can not read TreeQA " << endl ;
//     return ;
//   }
  
//   TBranch * qabranch = EMCAL()->TreeQA()->GetBranch("EMCAL") ; 
//   if (!qabranch) { 
//     if (fDebug)
//       cout << "WARNING: AliEMCALGetter::ReadTreeQA -> Cannot find QA Alarms for EMCAL" << endl ;
//     return ; 
//   }   
  
//   if(!Alarms())
//     PostQA() ; 

//   qabranch->SetAddress(AlarmsRef()) ;

//   qabranch->GetEntry(0) ;
 
//   PostQA("EMCAL") ; 
//   TFolder * alarmsF = Alarms() ; 
//   alarmsF->Clear() ; 
//   qabranch->SetAddress(&alarmsF) ;
//   qabranch->GetEntry(0) ;
  
//}
  
//____________________________________________________________________________ 
void AliEMCALGetter::ReadTreeR()
{
      // Read the reconstrunction tree gAlice->TreeR()

  if(gAlice->TreeR()== 0){
    if (fDebug) 
      cout <<   "WARNING: AliEMCALGetter::ReadTreeR: can not read TreeR " << endl ;
    return ;
  }
  
  // RecPoints 
  TObjArray * lob = static_cast<TObjArray*>(gAlice->TreeR()->GetListOfBranches()) ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * towerbranch = 0 ; 
  TBranch * preshobranch = 0 ; 
  TBranch * clusterizerbranch = 0 ; 
  Bool_t emcaltowerrpfound = kFALSE, emcalpreshorpfound = kFALSE, clusterizerfound = kFALSE ; 
  
  while ( (branch = static_cast<TBranch*>(next())) && (!emcaltowerrpfound || !emcalpreshorpfound || !clusterizerfound) ) 
    if(strcmp(branch->GetTitle(), fRecPointsTitle)==0) {
      if ( strcmp(branch->GetName(), "EMCALTowerRP")==0) {
	towerbranch = branch ; 
	emcaltowerrpfound = kTRUE ;
      }
      else if ( strcmp(branch->GetName(), "EMCALPreShoRP")==0) {
	preshobranch = branch ; 
	emcalpreshorpfound = kTRUE ;
      }
      else if(strcmp(branch->GetName(), "AliEMCALClusterizer")==0){
	clusterizerbranch = branch ; 
	clusterizerfound = kTRUE ; 
      }
    }

  if ( !emcaltowerrpfound ) {
    if (fDebug)
      cout << "WARNING: AliEMCALGetter::ReadTreeR -> Cannot find TowerRecPoints with title " 
	   << fRecPointsTitle << endl ;
    return ; 
  }   
  if ( !emcalpreshorpfound ) {
    if (fDebug)
      cout << "WARNING: AliEMCALGetter::ReadTreeR -> Cannot find PreShoRecPoints with title " 
	   << fRecPointsTitle << endl ;
    return ; 
  }   
  if ( !clusterizerfound ) {
    if (fDebug)
      cout << "WARNING: AliEMCALGetter::ReadTreeR -> Can not find Clusterizer with title " 
	   << fRecPointsTitle << endl ;
    return ; 
  }   
  
  // Read and Post the RecPoints
  if(!TowerRecPoints(fRecPointsTitle) )
    PostRecPoints(fRecPointsTitle) ;
  towerbranch->SetAddress(TowerRecPointsRef(fRecPointsTitle)) ;
  towerbranch->GetEntry(0) ;

  preshobranch->SetAddress(PreShoRecPointsRef(fRecPointsTitle)) ;
  preshobranch->GetEntry(0) ;
  
  if(!Clusterizer(fRecPointsTitle) )
    PostClusterizer(fRecPointsTitle) ;
  
  clusterizerbranch->SetAddress(ClusterizerRef(fRecPointsTitle)) ;
  clusterizerbranch->GetEntry(0) ;
 
  
  //------------------- TrackSegments ---------------------
//   next.Reset() ; 
//   TBranch * tsbranch = 0 ; 
//   TBranch * tsmakerbranch = 0 ; 
//   Bool_t emcaltsfound = kFALSE, tsmakerfound = kFALSE ; 
    
//   while ( (branch = static_cast<TBranch*>(next())) && (!emcaltsfound || !tsmakerfound) ) 
//     if(strcmp(branch->GetTitle(), fTrackSegmentsTitle)==0)  {
//       if ( strcmp(branch->GetName(), "EMCALTS")==0){
// 	tsbranch = branch ; 
// 	emcaltsfound = kTRUE ;
//       }
//       else if(strcmp(branch->GetName(), "AliEMCALTrackSegmentMaker")==0) {
// 	tsmakerbranch = branch ; 
// 	tsmakerfound  = kTRUE ; 
//       }
//     }
  
//   if ( !emcaltsfound || !tsmakerfound ) {
//     if (fDebug)
//       cout << "WARNING: AliEMCALGetter::ReadTreeR -> Cannot find TrackSegments and/or TrackSegmentMaker with name "
// 	   << fTrackSegmentsTitle << endl ;
//     return ; 
//   } 
  
//   // Read and Post the TrackSegments
//   if(!TrackSegments(fTrackSegmentsTitle))
//     PostTrackSegments(fTrackSegmentsTitle) ;
//   tsbranch->SetAddress(TrackSegmentsRef(fTrackSegmentsTitle)) ;
//   tsbranch->GetEntry(0) ;
  
//   // Read and Post the TrackSegment Maker
//   if(!TrackSegmentMaker(fTrackSegmentsTitle))
//     PostTrackSegmentMaker(fTrackSegmentsTitle) ;
//   tsmakerbranch->SetAddress(TSMakerRef(fTrackSegmentsTitle)) ;
//   tsmakerbranch->GetEntry(0) ;
  
  
//   //------------ RecParticles ----------------------------
//   next.Reset() ; 
//   TBranch * rpabranch = 0 ; 
//   TBranch * pidbranch = 0 ; 
//   Bool_t emcalrpafound = kFALSE, pidfound = kFALSE ; 
  
//   while ( (branch = static_cast<TBranch*>(next())) && (!emcalrpafound || !pidfound) ) 
//     if(strcmp(branch->GetTitle(), fRecParticlesTitle)==0) {   
//       if ( strcmp(branch->GetName(), "EMCALRP")==0) {   
// 	rpabranch = branch ; 
// 	emcalrpafound = kTRUE ;
//       }
//       else if (strcmp(branch->GetName(), "AliEMCALPID")==0) {
// 	pidbranch = branch ; 
// 	pidfound  = kTRUE ; 
//       }
//     }
  
//   if ( !emcalrpafound || !pidfound ) {
//     if (fDebug)
//       cout << "WARNING: AliEMCALGetter::ReadTreeR -> Cannot find RecParticles and/or PID with name " 
// 	   << fRecParticlesTitle << endl ;
//     return ; 
//   } 
  
//   // Read and Post the RecParticles
//   if(!RecParticles(fRecParticlesTitle))
//     PostRecParticles(fRecParticlesTitle) ;
//   rpabranch->SetAddress(RecParticlesRef(fRecParticlesTitle)) ;
//   rpabranch->GetEntry(0) ;
  
//   // Read and Post the PID
//   if(!PID(fRecParticlesTitle))
//     PostPID(fRecParticlesTitle) ;
//   pidbranch->SetAddress(PIDRef(fRecParticlesTitle)) ;
//   pidbranch->GetEntry(0) ;
  
  
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

  TString subdir(fHeaderFile) ;
  subdir.ReplaceAll("/","_") ; 

  if ( (subdir != "aliroot") && ( !folderslist->Contains(subdir) ) ){
    emcalF->AddFolder(subdir, ""); 
  }

  TIter next(folderslist) ; 
  TFolder * folder = 0 ; 
  TFile * file; 
  TTree * treeS = 0;
  while ( (folder = static_cast<TFolder*>(next())) ) {
   TString fileName(folder->GetName()) ; 
   fileName.ReplaceAll("_","/") ; 
   if(fHeaderFile.CompareTo(fileName) == 0 ) 
      treeS=gAlice->TreeS() ;
    else{
      file = static_cast<TFile*>(gROOT->GetFile(fileName)); 
      file->cd() ;
      
      // Get SDigits Tree header from file
      TString treeName("TreeS") ;
      treeName += event ; 
      treeS = dynamic_cast<TTree*>(gDirectory->Get(treeName.Data()));
    }
   if(!treeS){ // TreeS not found in header file
     
     if (fDebug)
       cout << "WARNING: AliEMCALGetter::ReadTreeS -> Cannot find TreeS in " << fHeaderFile << endl;
     
     TString searchFileName("") ; 
     
     if (SDigitizer())  // SDigitizer found in header file
       searchFileName = SDigitizer()->GetTitle() ;
     
     else if (Digitizer())  // Digitizer found in header file
       searchFileName = Digitizer()->GetSDigitsFileName() ; 
     
     else if (Clusterizer())  // Clusterizer found in header file
       searchFileName = Clusterizer()->GetSDigitsFileName() ; 
     
     if (treeS = TreeS(searchFileName)) { //found TreeS in the file which contains the hits
       if (fDebug) 
	 cout << "INFO: AliEMCALGetter::ReadTreeS -> TreeS found in " << searchFileName.Data() << endl ; 
       
     } else {
       cerr << "ERROR: AliEMCALGetter::ReadTreeS -> TreeS not found " << endl ; 
       return ;
     }
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
     }
     
     else if ( (strcmp(branch->GetName(), "AliEMCALSDigitizer")==0) && (strcmp(branch->GetTitle(), fSDigitsTitle)==0) ) {
       sdigitizerfound = kTRUE ; 
       sdigitizerBranch = branch ;
     }
   }
   if ( !emcalfound || !sdigitizerfound ) {
     if (fDebug)
       cout << "WARNING: AliEMCALGetter::ReadSDigits -> Digits and/or Digitizer branch with name " <<  fSDigitsTitle
	    << " not found" << endl ;
     return ; 
   }   
   
   if ( !folder->FindObject(fSDigitsTitle) )  
     PostSDigits(fSDigitsTitle,folder->GetName()) ;
   sdigitsBranch->SetAddress(SDigitsRef(fSDigitsTitle,folder->GetName())) ;
   sdigitsBranch->GetEntry(0) ;
   
   TString sdname(fSDigitsTitle) ;
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
    TString fileName(folder->GetName()) ; 
    fileName.ReplaceAll("_","/") ; 
    file   = static_cast<TFile*>(gROOT->GetFile(fileName)); 
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
    if (fDebug)
      cout << "WARNING: AliEMCALGetter::ReadTreeS -> Digits and/or Digitizer branch not found" << endl ;
    return ; 
  }   
  
  if (!folder || !(folder->FindObject(sdigitsBranch->GetTitle()) ) )
    PostSDigits(sdigitsBranch->GetTitle(),filename) ;

  sdigitsBranch->SetAddress(SDigitsRef(sdigitsBranch->GetTitle(),filename)) ;
  sdigitsBranch->GetEntry(0) ;
  
  TString sdname(sdigitsBranch->GetTitle()) ;
  sdname+=":" ;
  sdname+=filename ;
  if(!SDigitizer(sdigitsBranch->GetTitle()) )
    PostSDigitizer(sdigitsBranch->GetTitle(),filename) ;

  sdigitizerBranch->SetAddress(SDigitizerRef(sdname)) ;
  sdigitizerBranch->GetEntry(0) ;
  
}    


//____________________________________________________________________________ 
void AliEMCALGetter::ReadPrimaries()
{
  // Reads specific branches of primaries
  
  TClonesArray * ar = 0  ; 
  if(! (ar = Primaries()) ) { 
    PostPrimaries() ;
    ar = Primaries() ; 
  }
  ar->Delete() ; 
  
  if (TreeK(fHeaderFile)) { // treeK found in header file
    if (fDebug) 
      cout << "INFO: AliEMCALGetter::ReadPrimaries -> TreeK found in " << fHeaderFile.Data() << endl ; 
    fNPrimaries = gAlice->GetNtrack() ; 
  
  } else { // treeK not found in header file

    TString searchFileName("") ; 

    if (SDigitizer())  // SDigitizer found in header file
      searchFileName = SDigitizer()->GetTitle() ;

    else if (Digitizer())  // Digitizer found in header file
      searchFileName = Digitizer()->GetHitsFileName() ; 

    else if (Clusterizer())  // Clusterizer found in header file
      searchFileName = Clusterizer()->GetHitsFileName() ; 
    
    if (TreeK(searchFileName)) { //found TreeK in the file which contains the hits
      if (fDebug) 
	cout << "INFO: AliEMCALGetter::ReadPrimaries -> TreeK found in " << searchFileName.Data() << endl ; 
      fAlice->GetEvent(EventNumber()) ; 
      fNPrimaries = fAlice->GetNtrack() ; 
      
    } else {
      cerr << "ERROR: AliEMCALGetter::ReadPrimaries -> TreeK not  found " << endl ; 
      return ;
    }
    
  }
  Int_t index = 0 ; 
  for (index = 0 ; index < fNPrimaries; index++) { 
    new ((*ar)[index]) TParticle(*(Primary(index)));
  }
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

  // Carefull the order of reading is important 
  
  if( strstr(opt,"R") )
    ReadTreeR() ;

  if( strstr(opt,"D") )
    ReadTreeD() ;

  if(strstr(opt,"S") )
    ReadTreeS(event) ;

  if(strstr(opt,"H") )
    ReadTreeH() ;

 // if( strstr(opt,"Q") )
  //  ReadTreeQA() ;

  if( strstr(opt,"P") || (strcmp(opt,"")==0) )
    ReadPrimaries() ;
  
}

//____________________________________________________________________________ 
TObject * AliEMCALGetter::ReturnO(TString what, TString name, TString file) const 
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
    file.ReplaceAll("/","_") ; 
    TString path = "EMCAL/" + file  ; 
    folder = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject(path.Data())) ; 
    if (folder) { 
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
  else if ( what.CompareTo("TowerRecPoints") == 0 ) {
    folder = dynamic_cast<TFolder *>(fRecoFolder->FindObject("EMCAL/TowerRecPoints")) ; 
    if (folder) { 
      if (name.IsNull())
	name = fRecPointsTitle ; 
      emcalO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    } 
  }
  else if ( what.CompareTo("PreShoRecPoints") == 0 ) {
    folder = dynamic_cast<TFolder *>(fRecoFolder->FindObject("EMCAL/PreShoRecPoints")) ; 
    if (folder) { 
      if (name.IsNull())
	name = fRecPointsTitle ; 
      emcalO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    }   
  }
  /*
  else if ( what.CompareTo("TrackSegments") == 0 ) {
    folder = dynamic_cast<TFolder *>(fRecoFolder->FindObject("EMCAL/TrackSegments")) ; 
    if (folder) { 
      if (name.IsNull())
	name = fTrackSegmentsTitle ; 
      emcalO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    }   
  }
  else if ( what.CompareTo("RecParticles") == 0 ) {
    folder = dynamic_cast<TFolder *>(fRecoFolder->FindObject("EMCAL/RecParticles")) ; 
   if (folder) { 
      if (name.IsNull())
	name = fRecParticlesTitle ; 
      emcalO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    }   
 }
  else if ( what.CompareTo("Alarms") == 0 ){ 
    if (name.IsNull() ) 
      emcalO = dynamic_cast<TObject *>(fQAFolder->FindObject("EMCAL")) ;  
    else {
      folder = dynamic_cast<TFolder *>(fQAFolder->FindObject("EMCAL")) ; 
      if (!folder) 
	emcalO = 0 ; 
      else 
	emcalO = dynamic_cast<TObject *>(folder->FindObject(name)) ;  
    }
  }
*/
  if (!emcalO) {
    if(fDebug)
      cout << "WARNING : AliEMCALGetter::ReturnO -> Object " << what << " not found in " << folder->GetName() << endl ; 
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
  if ( what.CompareTo("Clusterizer") == 0 ) 
    search = "Reconstructioner" ; 
  else if ( what.CompareTo("TrackSegmentMaker") == 0 ) 
    search = "Reconstructioner" ; 
  else if ( what.CompareTo("PID") == 0 ) 
    search = "Reconstructioner" ; 
  else if ( what.CompareTo("QATasks") == 0 ) 
    search = "QA" ; 
  
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
  } else  if (what.CompareTo("Clusterizer") == 0){ 
    if ( name.IsNull() )
      name =  fRecPointsTitle ;
    name.Append(":clu") ;
  }
 //  else  if (what.CompareTo("TrackSegmentMaker") == 0){ 
//     if ( name.IsNull() )
//       name =  fTrackSegmentsTitle ;
//     name.Append(":tsm") ;
//   }
//   else  if (what.CompareTo("PID") == 0){ 
//     if ( name.IsNull() )
//       name =  fRecParticlesTitle ;
//     name.Append(":pid") ;
//   }
//   else  if (what.CompareTo("QATasks") == 0){ 
//     if ( name.IsNull() )
//       return emcalT ;
//   }
  
  TIter it(list) ;
  TTask * task = 0 ; 
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(name)){
    return task ;}
  }
  
  if(fDebug)
    cout << "WARNING: AliEMCALGetter::ReturnT -> Task " << search << "/" << name << " not found!" << endl ; 
  return 0 ;
}


//____________________________________________________________________________ 
void AliEMCALGetter::RemoveTask(TString opt, TString name) const 
{
  // remove a task from the folder
  // path is fTasksFolder/SDigitizer/EMCAL/name

  TTask * task  = 0 ; 
  TTask * emcal = 0 ; 
  TList * lofTasks = 0 ; 

  if (opt == "S") { // SDigitizer
    task = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ;
    if (!task) 
      return ; 
  }
  
  else if (opt == "D") { // Digitizer
    task = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ;
    if (!task) 
      return ; 
  }
  else if (opt == "C") { // Clusterizer
    task = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ;
    if (!task) 
      return ; 
  }
  else {
    cerr << "WARNING: AliEMCALGetter::RemoveTask -> Unknown option " << opt.Data() << endl ; 
    return ; 
  }
    
  emcal =  dynamic_cast<TTask*>(task->GetListOfTasks()->FindObject("EMCAL")) ;
  if (!emcal)
    return ; 
  lofTasks = emcal->GetListOfTasks() ;
  if (!lofTasks) 
    return ; 
  TObject * obj = lofTasks->FindObject(name) ; 
  if (obj) 
    lofTasks->Remove(obj) ;
}
  
//____________________________________________________________________________ 
void AliEMCALGetter::RemoveObjects(TString opt, TString name) const 
{
  // remove SDigits from the folder
  // path is fSDigitsFolder/fHeaderFileName/name

  TFolder * emcal     = 0 ; 
  TFolder * emcalmain = 0 ; 

  if (opt == "H") { // Hits
    emcal = dynamic_cast<TFolder*>(fHitsFolder->FindObject("EMCAL")) ;
    if (!emcal) 
      return ;
    name = "Hits" ; 
  }
  
  else if ( opt == "S") { // SDigits
    emcalmain = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("EMCAL")) ;
    if (!emcalmain) 
      return ;
    emcal = dynamic_cast<TFolder*>(emcalmain->FindObject(fHeaderFile)) ;
    if (!emcal) 
      return ;
  }

  else if (opt == "D") { // Digits
    emcal = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("EMCAL")) ;
    if (!emcal) 
      return ;
  }

  else if (opt == "RT") { // Tower RecPoints
    emcal = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/TowerRecPoints")) ;
    if (!emcal) 
      return ;
  }

  else if (opt == "RP") { // Preshower RecPoints
    emcal = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/PreShoRecPoints")) ;
    if (!emcal) 
      return ;
  }

  else {
    cerr << "WARNING: AliEMCALGetter::RemoveObjects -> Unknown option " << opt.Data() << endl ; 
    return ; 
  }
  
  TObjArray * ar  = dynamic_cast<TObjArray*>(emcal->FindObject(name)) ; 
  if (ar) { 
    emcal->Remove(ar) ;
    ar->Delete() ; 
    delete ar ; 
  }

  if (opt == "S") 
    emcalmain->Remove(emcal) ; 
  
}

//____________________________________________________________________________ 
void AliEMCALGetter::RemoveSDigits() const 
{
  TFolder * emcal= dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("EMCAL")) ;
  if (!emcal) 
    return ;
  
  emcal->SetOwner() ; 
  emcal->Clear() ; 
}
