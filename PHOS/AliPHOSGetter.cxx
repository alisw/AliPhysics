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
//  An example of how to use (see also class AliPHOSAnalyser):
//  AliPHOSGetter * gime = AliPHOSGetter::GetInstance("galice.root","test") ;
//  for(Int_t irecp = 0; irecp < gime->NRecParticles() ; irecp++)
//     AliPHOSRecParticle * part = gime->RecParticle(1) ;
//     ................
//  please->GetEvent(event) ;    // reads new event from galice.root
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//*--         Completely redesigned by Dmitri Peressounko March 2001  
//
//*-- YS June 2001 : renamed the original AliPHOSIndexToObject and make
//*--         systematic usage of TFolders without changing the interface        
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TObjString.h"
#include "TFolder.h"
#include "TParticle.h"

// --- Standard library ---
#include <iostream.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliConfig.h"
#include "AliPHOSGetter.h"
#include "AliPHOS.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSClusterizer.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSPID.h" 
#include "AliPHOSPIDv1.h" 
#include "AliPHOSGeometry.h"

ClassImp(AliPHOSGetter)
  
  AliPHOSGetter * AliPHOSGetter::fgObjGetter = 0 ; 

//____________________________________________________________________________ 
AliPHOSGetter::AliPHOSGetter(const char* headerFile, const char* branchTitle )
{
  //Initialize  all lists

  fHeaderFile         = headerFile ; 
  fBranchTitle        = branchTitle ;
  fSDigitsTitle       = branchTitle ; 
  fDigitsTitle        = branchTitle ; 
  fRecPointsTitle     = branchTitle ; 
  fRecParticlesTitle  = branchTitle ; 
  fTrackSegmentsTitle = branchTitle ; 

  fPrimaries = new TObjArray(1) ;

  fModuleFolder    = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Configuration/Modules")); 
  fPrimariesFolder = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data")); 
  fHitsFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/Hits")); 
  fSDigitsFolder   = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/SDigits")); 
  fDigitsFolder    = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Event/Data")); 
  fRecoFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Event/RecData")); 
  fQAFolder        = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Conditions/QA")); 
  fTasksFolder     = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Tasks")) ; 

  if ( fHeaderFile != "aliroot"  ) { // to call the getter without a file

    //open headers file
    TFile * file = static_cast<TFile*>(gROOT->GetFile(fHeaderFile.Data() ) ) ;
    
    if(file == 0){    //if file was not opened yet, read gAlice
      if(fHeaderFile.Contains("rfio")) // if we read file using HPSS
	file =	TFile::Open(fHeaderFile.Data(),"update") ;
      else
	file = new TFile(fHeaderFile.Data(),"update") ;
      
      if (!file->IsOpen()) {
	cerr << "ERROR : AliPHOSGetter::AliPHOSGetter -> Cannot open " << fHeaderFile.Data() << endl ; 
	abort() ; 
      }
      
      gAlice = static_cast<AliRun *>(file->Get("gAlice")) ;
    }
  }

  if (!gAlice) {
    cerr << "ERROR : AliPHOSGetter::AliPHOSGetter -> Cannot find gAlice in " << fHeaderFile.Data() << endl ; 
    abort() ; 
  }
  if (!PHOS()) {
    if (fDebug)
      cout << "INFO: AliPHOSGetter -> Posting PHOS to Folders" << endl ; 
    AliConfig * conf = AliConfig::Instance() ; 
    conf->Add(static_cast<AliDetector*>(gAlice->GetDetector("PHOS"))) ; 
    conf->Add(static_cast<AliModule*>(gAlice->GetDetector("PHOS"))) ; 
  }
  
  fDebug=0;
}
//____________________________________________________________________________ 
AliPHOSGetter::~AliPHOSGetter(){

}

//____________________________________________________________________________ 
void AliPHOSGetter::CreateWhiteBoard() const
{

}

//____________________________________________________________________________ 
AliPHOSGetter * AliPHOSGetter::GetInstance()
{
  // Returns the pointer of the unique instance already defined
  
  AliPHOSGetter * rv = 0 ;
  if ( fgObjGetter )
    rv = fgObjGetter ;
  else
    cout << "AliPHOSGetter::GetInstance ERROR: not yet initialized" << endl ;

  return rv ;
}

//____________________________________________________________________________ 
AliPHOSGetter * AliPHOSGetter::GetInstance(const char* headerFile,
					   const char* branchTitle)
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed 

  if ( fgObjGetter )    
    if((fgObjGetter->fBranchTitle.CompareTo(branchTitle) == 0) && 
       (fgObjGetter->fHeaderFile.CompareTo(headerFile)==0))
      return fgObjGetter ;
    else
      fgObjGetter->~AliPHOSGetter() ;  // delete it if already exists another version
  
  fgObjGetter = new AliPHOSGetter(headerFile,branchTitle) ; 
  
  // Posts a few item to the white board (folders)
  // fgObjGetter->CreateWhiteBoard() ;
    
  return fgObjGetter ; 
  
}

//____________________________________________________________________________ 
const AliPHOS * AliPHOSGetter::PHOS() 
{
  // returns the PHOS object 
  AliPHOS * phos = dynamic_cast<AliPHOS*>(fModuleFolder->FindObject("PHOS")) ;  
  if (!phos) 
    if (fDebug)
      cout << "WARNING: AliPHOSGetter::PHOS -> PHOS module not found in Folders" << endl ; 
  return phos ; 
}  

//____________________________________________________________________________ 
const AliPHOSGeometry * AliPHOSGetter::PHOSGeometry() 
{
  AliPHOSGeometry * rv = 0 ; 
  if (PHOS() )
    rv =  PHOS()->GetGeometry() ;
  return rv ; 
} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostPrimaries(void) const 
{  //------- Primaries ----------------------

  // the hierarchy is //Folders/RunMC/Event/Data/Primaries
  
  TFolder * primariesFolder = dynamic_cast<TFolder*>(fPrimariesFolder->FindObject("Primaries")) ; 
  if ( !primariesFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post Primaries -> Folder //" << fPrimariesFolder << "/Primaries/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post Primaries -> Adding Folder //" << fPrimariesFolder << "/Primaries/"  << endl;
    }
    primariesFolder = fPrimariesFolder->AddFolder("Primaries", "Primaries particles from TreeK") ; 
  }    
  TClonesArray *primaries=  new TClonesArray("TParticle",1000) ;
  primaries->SetName("Primaries") ;
  primariesFolder->Add(primaries) ; 
  
  return kTRUE;
} 

//____________________________________________________________________________ 
void * AliPHOSGetter::PrimariesRef(void) const 
{  //------- Primaries ----------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/Primaries
  if ( !fPrimariesFolder ) {
    cerr << "ERROR: AliPHOSGetter::Post PrimariesRef -> Folder //" << fPrimariesFolder << " not found!" << endl;
    return 0;
  }    
 
  TFolder * primariesFolder = dynamic_cast<TFolder *>(fPrimariesFolder->FindObject("Primaries")) ;
  if ( !primariesFolder ) {
    cerr << "ERROR: AliPHOSGetter::Post PrimariesRef -> Folder //" << fPrimariesFolder << "/Primaries/ not found!" << endl;  
    return 0;
  }
 
  TObject * p = primariesFolder->FindObject("Primaries") ;
  if(!p) {
    cerr << "ERROR: AliPHOSGetter::PrimariesRef -> " << primariesFolder->GetName() << "/Primaries not found !" << endl ; 
    return 0 ;
  }
  else
    return static_cast<void *>(primariesFolder->GetListOfFolders()->GetObjectRef(p)) ;
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostHits(void) const 
{  //------- Hits ----------------------

  // the hierarchy is //Folders/RunMC/Event/Data/PHOS/Hits
  
  TFolder * phosFolder = dynamic_cast<TFolder*>(fHitsFolder->FindObject("PHOS")) ; 
  if ( !phosFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post H -> Folder //" << fHitsFolder << "/PHOS/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post H -> Adding Folder //" << fHitsFolder << "/PHOS/"  << endl;
    }
    phosFolder = fHitsFolder->AddFolder("PHOS", "Hits from PHOS") ; 
  }    
  TClonesArray *hits=  new TClonesArray("AliPHOSHit",1000) ;
  hits->SetName("Hits") ;
  phosFolder->Add(hits) ; 
  
  return kTRUE;
} 

//____________________________________________________________________________ 
void * AliPHOSGetter::HitsRef(void) const 
{  //------- Hits ----------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/PHOS/Hits
  if ( !fHitsFolder ) {
    cerr << "ERROR: AliPHOSGetter::Post H -> Folder //" << fHitsFolder << " not found!" << endl;
    return 0;
  }    
 
  TFolder * phosFolder = dynamic_cast<TFolder *>(fHitsFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) {
    cerr << "ERROR: AliPHOSGetter::Post HRef -> Folder //" << fHitsFolder << "/PHOS/ not found!" << endl;  
    return 0;
  }
 
  TObject * h = phosFolder->FindObject("Hits") ;
  if(!h) {
    cerr << "ERROR: AliPHOSGetter::HRef -> " << phosFolder->GetName() << "/Hits not found !" << endl ; 
    return 0 ;
  }
  else
    return static_cast<void *>(phosFolder->GetListOfFolders()->GetObjectRef(h)) ;
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostSDigits(const char * name, const char * headerFile) const 
{  //---------- SDigits -------------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/PHOS/SDigits/headerFile/sdigitsname
  // because you can have sdigits from several hit files for mixing
  
  TFolder * phosFolder = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post S -> Folder //" << fSDigitsFolder << "/PHOS/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post S -> Adding Folder //" << fHitsFolder << "/PHOS/" << endl;
    }
    phosFolder = fSDigitsFolder->AddFolder("PHOS", "SDigits from PHOS") ; 
  }    
  TString subdir(headerFile) ;
  TFolder * phosSubFolder = dynamic_cast<TFolder*>(phosFolder->FindObject(subdir)) ; 
  if ( !phosSubFolder ) 
    phosSubFolder = phosFolder->AddFolder(subdir, ""); 
  
  TObject * sd  = phosSubFolder->FindObject(name); 
  if ( sd ) {
    if (fDebug)
      cerr <<"INFO: AliPHOSGetter::Post S -> Folder " << subdir 
	   << " already exists!" << endl ;  
  }else{
    TClonesArray * sdigits = new TClonesArray("AliPHOSDigit",1) ;
    sdigits->SetName(name) ;
    phosSubFolder->Add(sdigits) ;
  }
  
  return kTRUE;
} 
//____________________________________________________________________________ 
void * AliPHOSGetter::SDigitsRef(const char * name, const char * file) const 
{  //------- SDigits ----------------------
  
  // the hierarchy is //Folders/RunMC/Event/Data/PHOS/SDigits/filename/SDigits

  if ( !fSDigitsFolder ) {
    cerr << "ERROR: AliPHOSGetter::Post SRef -> Folder //" << fSDigitsFolder << " not found!" << endl;
    return 0;
  }    
 
  TFolder * phosFolder = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) {
    cerr << "ERROR: AliPHOSGetter::Post SRef -> Folder //" << fSDigitsFolder << "/PHOS/ not found!" << endl;
    return 0;
  }

  TFolder * phosSubFolder = 0 ;
  if(file)
    phosSubFolder = dynamic_cast<TFolder *>(phosFolder->FindObject(file)) ;
  else
    phosSubFolder = dynamic_cast<TFolder *>(phosFolder->FindObject(fHeaderFile)) ;
  
  if(!phosSubFolder) {
    cerr << "ERROR: AliPHOSGetter::Post SRef -> Folder //Folders/RunMC/Event/Data/PHOS/" << file << "not found!" << endl;
    return 0;
  }

  TObject * dis = phosSubFolder->FindObject(name) ;
  if(!dis)
    return 0 ;
  else
    return static_cast<void *>(phosSubFolder->GetListOfFolders()->GetObjectRef(dis)) ;

}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostSDigitizer(AliPHOSSDigitizer * sdigitizer) const 
{  //---------- SDigitizer -------------------------
    
  // the hierarchy is //Folders/Tasks/SDigitizer/PHOS/sdigitsname


  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ; 

  if ( !sd ) {
    cerr << "ERROR: AliPHOSGetter::Post Ser -> Task //" << fTasksFolder << "/SDigitizer not found!" << endl;
    return kFALSE ;
  }        
  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Ser ->//" << fTasksFolder << "/SDigitizer/PHOS/ not found!" << endl;  
      cout <<"INFO: AliPHOSGetter::Post Ser -> Adding //" << fTasksFolder << "/SDigitizer/PHOS/" << endl;
    }
    phos = new TTask("PHOS", "") ; 
    sd->Add(phos) ; 
  } 
  AliPHOSSDigitizer * phossd  = dynamic_cast<AliPHOSSDigitizer *>(phos->GetListOfTasks()->FindObject( sdigitizer->GetName() )); 
  if (phossd) { 
    if (fDebug)
      cout << "INFO: AliPHOSGetter::Post Ser -> Task " << sdigitizer->GetName() << " already exists" << endl ; 
    phos->GetListOfTasks()->Remove(phossd) ;
  }
  phos->Add(sdigitizer) ;	
  return kTRUE; 
  
}

//____________________________________________________________________________ 
void * AliPHOSGetter::SDigitizerRef(const char * name) const 
{  

  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ; 
  if ( !sd ) {
    cerr << "ERROR: AliPHOSGetter::Post SerRef -> Task //" << fTasksFolder << "/SDigitizer not found!" << endl;
    abort();
  }        

  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post SerRef ->  //" << fTasksFolder << "/SDigitizer/PHOS not found!" << endl;
    abort();
  }        

  TTask * task = dynamic_cast<TTask*>(phos->GetListOfTasks()->FindObject(name)) ; 

  return static_cast<void *>(phos->GetListOfTasks()->GetObjectRef(task)) ;

}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostSDigitizer(const char * name, const char * file) const 
{  //---------- SDigitizer -------------------------
  
 // the hierarchy is //Folders/Tasks/SDigitizer/PHOS/sdigitsname


  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ; 
  if ( !sd ) {
    cerr << "ERROR: AliPHOSGetter::Post Ser -> Task //" << fTasksFolder << "/SDigitizer not found!" << endl;
    return kFALSE ;
  }        

  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Ser ->  //" << fTasksFolder << "/SDigitizer/PHOS/ not found!" << endl;
      cout <<"INFO: AliPHOSGetter::Post Ser -> Adding  //" << fTasksFolder << "/SDigitizer/PHOS" << endl;
    }
    phos = new TTask("PHOS", "") ; 
    sd->Add(phos) ; 
  } 

  TString sdname(name) ;
  sdname.Append(":") ;
  sdname.Append(file);
  AliPHOSSDigitizer * phossd  = dynamic_cast<AliPHOSSDigitizer *>(phos->GetListOfTasks()->FindObject( sdname )); 
  if (!phossd) {
    phossd = new AliPHOSSDigitizer() ;  
    //Note, we can not call constructor with parameters: it will call Getter and scrud up everething
    phossd->SetName(sdname) ;
    phossd->SetTitle(file) ;
    phos->Add(phossd) ;	
  }
  return kTRUE; 
  
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostDigits(const char * name) const 
{  //---------- Digits -------------------------

  // the hierarchy is //Folders/Run/Event/Data/PHOS/SDigits/name

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("PHOS")) ;

  if ( !phosFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post D -> Folder //" << fDigitsFolder << "/PHOS/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post D -> Adding Folder //" << fDigitsFolder << "/PHOS/" << endl;
    }
    phosFolder = fDigitsFolder->AddFolder("PHOS", "Digits from PHOS") ;  
  }    
 
  TObject*  dig = phosFolder->FindObject( name ) ;
  if ( !dig ) {
    TClonesArray * digits = new TClonesArray("AliPHOSDigit",1000) ;
    digits->SetName(name) ;
    phosFolder->Add(digits) ;  
  }
  return kTRUE; 
}

//____________________________________________________________________________ 
void * AliPHOSGetter::DigitsRef(const char * name) const 
{ //------- Digits ----------------------
  
  // the hierarchy is //Folders/Run/Event/Data/PHOS/Digits/name

  if ( !fDigitsFolder ) {
    cerr << "ERROR: AliPHOSGetter::Post DRef -> Folder //" << fDigitsFolder << " not found!" << endl;
    return 0;
  }    
  
  TFolder * phosFolder  = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("PHOS")) ; 
  if ( !phosFolder ) {
    cerr << "ERROR: AliPHOSGetter::DRef -> Folder //" << fDigitsFolder << "/PHOS/ not found!" << endl;
    return 0;
  }    

  TObject * d = phosFolder->FindObject(name) ;
  if(!d)
    return 0 ;
  else
    return static_cast<void *>(phosFolder->GetListOfFolders()->GetObjectRef(d)) ;

}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostDigitizer(AliPHOSDigitizer * digitizer) const 
{  //---------- Digitizer -------------------------
  
  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 

  if ( !sd ) {
    cerr << "ERROR: AliPHOSGetter::Post Der -> Task //" << fTasksFolder << "/Digitizer not found!" << endl;
    return kFALSE ;
  }        
  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Der ->  //" << fTasksFolder << "/Digitizer/PHOS not found!" << endl;
      cout <<"INFO: AliPHOSGetter::Post Der -> Adding //" << fTasksFolder << "/Digitizer/PHOS" << endl; 
    }
    phos = new TTask("PHOS", "") ; 
    sd->Add(phos) ; 
  } 

    AliPHOSDigitizer * phosd = dynamic_cast<AliPHOSDigitizer*>(phos->GetListOfTasks()->FindObject(digitizer->GetName())) ; 
    if (phosd) { 
      phosd->Delete() ;
      phos->GetListOfTasks()->Remove(phosd) ;
    }
    phos->Add(digitizer) ; 
    return kTRUE; 
}  

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostDigitizer(const char * name) const 
{  //---------- Digitizer -------------------------
  
 // the hierarchy is //Folders/Tasks/SDigitizer/PHOS/sdigitsname

  TTask * d  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 
  if ( !d ) {
    cerr << "ERROR: AliPHOSGetter::Post Der -> Task //" << fTasksFolder << "/Digitizer not found!" << endl;
    return kFALSE ;
  }        

  TTask * phos = dynamic_cast<TTask*>(d->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Der -> //" << fTasksFolder << "/Digitizer/PHOS not found!" << endl; 
      cout <<"INFO: AliPHOSGetter::Post Der -> Adding //" << fTasksFolder << "/Digitizer/PHOS" << endl;
    }
    phos = new TTask("PHOS", "") ; 
    d->Add(phos) ; 
} 

  AliPHOSDigitizer * phosd = dynamic_cast<AliPHOSDigitizer*>(phos->GetListOfTasks()->FindObject(name)) ; 
  if (!phosd) { 
    phosd = new AliPHOSDigitizer() ;
    phosd->SetName(fDigitsTitle) ;
    phosd->SetTitle(fHeaderFile) ;
    phos->Add(phosd) ;
  }
  return kTRUE;  
}

//____________________________________________________________________________ 
void * AliPHOSGetter::DigitizerRef(const char * name) const 
{  
  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 
  if ( !sd ) {
    cerr << "ERROR: AliPHOSGetter::Post DerRef -> Task //" << fTasksFolder << "/Digitizer not found!" << endl;
    abort();
  }        

  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post DerRef ->  //" << fTasksFolder << "/Digitizer/PHOS" << endl;
    abort();
  }        

  TTask * task = dynamic_cast<TTask*>(phos->GetListOfTasks()->FindObject(name)) ; 

  return static_cast<void *>(phos->GetListOfTasks()->GetObjectRef(task)) ;

}
 
//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostRecPoints(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/EMCARecPoints/name
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/CPVRecPoints/name

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS")) ; 
  
  if ( !phosFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post RPo -> Folder //" << fRecoFolder << "/PHOS/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post Rpo -> Adding Folder //" << fRecoFolder << "/PHOS/" << endl;
    }
    phosFolder = fRecoFolder->AddFolder("PHOS", "Reconstructed data from PHOS") ;  
  }    
  
  // EMCA RecPoints 
  TFolder * phosRPoEMCAFolder  = dynamic_cast<TFolder*>(phosFolder->FindObject("EMCARecPoints")) ;
  if ( !phosRPoEMCAFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post RPo -> Folder //" << fRecoFolder << "/PHOS/EMCARecPoints/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post Rpo -> Adding Folder //" << fRecoFolder << "/PHOS/EMCARecPoints not found!" << endl;
    }
    phosRPoEMCAFolder = phosFolder->AddFolder("EMCARecPoints", "EMCA RecPoints from PHOS") ;  
  }    
  
  TObject * erp = phosFolder->FindObject( name ) ;
  if ( !erp )   {
    TObjArray * emcrp = new TObjArray(100) ;
    emcrp->SetName(name) ;
    phosRPoEMCAFolder->Add(emcrp) ;  
  }

  // CPV RecPoints 
  TFolder * phosRPoCPVFolder  = dynamic_cast<TFolder*>(phosFolder->FindObject("CPVRecPoints")) ;
  if ( !phosRPoCPVFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post RPo -> Folder //" << fRecoFolder << "/PHOS/CPVRecPoints/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post Rpo -> Adding Folder //" << fRecoFolder << "/PHOS/CPVRecPoints/" << endl;
    }
    phosRPoCPVFolder = phosFolder->AddFolder("CPVRecPoints", "CPV RecPoints from PHOS") ;  
  }    
  
  TObject * crp =  phosRPoCPVFolder->FindObject( name ) ;
  if ( !crp )   {
    TObjArray * cpvrp = new TObjArray(100) ;
    cpvrp->SetName(name) ;
    phosRPoCPVFolder->Add(cpvrp) ;  
  }
  return kTRUE; 
}

//____________________________________________________________________________ 
void * AliPHOSGetter::EmcRecPointsRef(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/EMCARecPoints/name
   
  if ( !fRecoFolder ) {
    cerr << "ERROR: AliPHOSGetter::EmcRecPointsRef -> Folder //" << fRecoFolder << " not found!" << endl;
    return 0 ; 
  }    

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/EMCARecPoints")) ; 
  if ( !phosFolder ) {
    cerr << "ERROR: AliPHOSGetter::EmcRecPointsRef -> Folder //" << fRecoFolder << "/PHOS/EMCARecPoints/ not found!" << endl;
    return 0;
  }    


  TObject * erp = phosFolder->FindObject(name ) ;
  if ( !erp )   {
    return 0 ;
  }
  return static_cast<void *>(phosFolder->GetListOfFolders()->GetObjectRef(erp)) ;

} 

//____________________________________________________________________________ 
void * AliPHOSGetter::CpvRecPointsRef(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/CPVRecPoints/name
   
  if ( !fRecoFolder ) {
    cerr << "ERROR: AliPHOSGetter::EmcRecPointsRef -> Folder //" << fRecoFolder << " not found!" << endl;
    return 0 ; 
  }    

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/CPVRecPoints")) ; 
  if ( !phosFolder ) {
    cerr << "ERROR: AliPHOSGetter::CpvRecPointsRef -> Folder //" << fRecoFolder << "/PHOS/CPVRecPoints/" << endl;
    return 0;
  }    

  TObject * crp = phosFolder->FindObject(name ) ;
  if ( !crp )   {
    return 0 ;
  }
  return static_cast<void *>(phosFolder->GetListOfFolders()->GetObjectRef(crp)) ;

} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostClusterizer(AliPHOSClusterizer * clu) const 
{ // ------------------ AliPHOSClusterizer ------------------------
  
  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliPHOSGetter::Post Rer -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Rer -> //" << fTasksFolder << "/ReconstructionerPHOS not found!" << endl; 
      cout <<"INFO: AliPHOSGetter::Post Rer -> Adding //" << fTasksFolder << "/Reconstructioner/PHOS" << endl; 
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  AliPHOSClusterizer * phoscl = dynamic_cast<AliPHOSClusterizer*>(phos->GetListOfTasks()->FindObject(clu->GetName())) ; 
  if (phoscl) { 
    if (fDebug)
      cout << "INFO: AliPHOSGetter::Post Rer -> Task " << clu->GetName() << " already exists" << endl ; 
    phoscl->Delete() ; 
    phos->GetListOfTasks()->Remove(phoscl) ;
  }
  phos->Add(clu) ;      
  return kTRUE; 
} 

//____________________________________________________________________________ 
void * AliPHOSGetter::ClusterizerRef(const char * name) const 
{ // ------------------ AliPHOSClusterizer ------------------------
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliPHOSGetter::Post RerRef -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    cerr <<"WARNING: AliPHOSGetter::Post RerRef -> //" << fTasksFolder << "/Reconstructioner/PHOS" << endl; 
    return 0 ; 
  }   

  TList * l = phos->GetListOfTasks() ; 
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
    return static_cast<void *>(l->GetObjectRef(clu)) ;
  else
    return 0 ;
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostClusterizer(const char * name) const 
{ // ------------------ AliPHOSClusterizer ------------------------

  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliPHOSGetter::Post Rer -> Task//" << fTasksFolder << "/Reconstructioner not found!" << endl; 
    return kFALSE ;
  }        
  
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Rer -> //" << fTasksFolder << "/Reconstructioner/PHOS not found!" << endl;
      cout <<"INFO: AliPHOSGetter::Post Rer -> Adding //" << fTasksFolder << "/Reconstructioner/PHOS" << endl;
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  AliPHOSClusterizer * phoscl = new AliPHOSClusterizerv1() ;
  TString clun(name) ;
  clun+=":clu-v1" ;
  phoscl->SetName(clun) ;
  phos->Add(phoscl) ;
  return kTRUE; 
  
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostTrackSegments(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/TrackSegments/name

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS")) ; 
  
  if ( !phosFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post TS -> Folder //" << fRecoFolder << "/PHOS/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post TS -> Adding Folder //" << fRecoFolder << "/PHOS" << endl;
    }
    phosFolder = fRecoFolder->AddFolder("PHOS", "Reconstructed data from PHOS") ;  
  }    

  TFolder * phosTSFolder  = dynamic_cast<TFolder*>(phosFolder->FindObject("TrackSegments")) ;
  if ( !phosTSFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post TS -> Folder//" << fRecoFolder << "/PHOS/TrackSegments/ not found!" << endl; 
      cout << "INFO:    AliPHOSGetter::Post TS -> Adding Folder //" << fRecoFolder << "/PHOS/TrackSegments/" << endl; 
    }
    phosTSFolder = phosFolder->AddFolder("TrackSegments", "TrackSegments from PHOS") ;  
  }    
  
  TObject * tss =  phosTSFolder->FindObject( name ) ;
  if (!tss) {
    TClonesArray * ts = new TClonesArray("AliPHOSTrackSegment",100) ;
    ts->SetName(name) ;
    phosTSFolder->Add(ts) ;  
  }
  return kTRUE; 
} 

//____________________________________________________________________________ 
void * AliPHOSGetter::TrackSegmentsRef(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/TrackSegments/name

 if ( !fRecoFolder ) {
    cerr << "ERROR: AliPHOSGetter::TrackSegmentsRef -> Folder //" << fRecoFolder << "not found!" << endl;
    return 0 ; 
  }    

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/TrackSegments")) ; 
  if ( !phosFolder ) {
    cerr << "ERROR: AliPHOSGetter::TrackSegmentsRef -> Folder //" << fRecoFolder << "/PHOS/TrackSegments/ not found!" << endl;
    return 0;
  }    
  
  TObject * tss =  phosFolder->FindObject(name) ;
  if (!tss) {
    return 0 ;  
  }
  return static_cast<void *>(phosFolder->GetListOfFolders()->GetObjectRef(tss)) ;
} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostTrackSegmentMaker(AliPHOSTrackSegmentMaker * tsmaker) const 
{ //------------Track Segment Maker ------------------------------
  
  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliPHOSGetter::Post Ter -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Rer -> //" << fTasksFolder << "/Reconstructioner/PHOS not found!" << endl; 
      cout <<"INFO: AliPHOSGetter::Post Rer -> Adding //" << fTasksFolder << "/Reconstructioner/PHOS" << endl;
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  AliPHOSTrackSegmentMaker * phosts = 
    dynamic_cast<AliPHOSTrackSegmentMaker*>(phos->GetListOfTasks()->FindObject(tsmaker->GetName())) ; 
  if (phosts) { 
    phosts->Delete() ;
    phos->GetListOfTasks()->Remove(phosts) ;
  }
  phos->Add(tsmaker) ;      
  return kTRUE; 
  
} 
//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostTrackSegmentMaker(const char * name) const 
{ //------------Track Segment Maker ------------------------------
  
  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname
  
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 
  
  if ( !tasks ) {
    cerr << "ERROR: AliPHOSGetter::Post Ter -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
  
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Rer -> //" << fTasksFolder << "/Reconstructioner/PHOS not found!" << endl; 
      cout <<"INFO: AliPHOSGetter::Post Rer -> Adding //" << fTasksFolder << "/Reconstructioner/PHOS" << endl;
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  AliPHOSTrackSegmentMaker * phosts = 
    dynamic_cast<AliPHOSTrackSegmentMaker*>(phos->GetListOfTasks()->FindObject(name)) ; 
  if (!phosts) { 
    phosts = new AliPHOSTrackSegmentMakerv1() ;
    TString tsn(name);
    tsn+=":tsm-v1" ;
    phosts->SetName(tsn) ;
    phos->Add(phosts) ;      
  }
  return kTRUE; 
  
} 

//____________________________________________________________________________ 
void * AliPHOSGetter::TSMakerRef(const char * name) const 
{ //------------Track Segment Maker ------------------------------
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliPHOSGetter::Post TerRef -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    cerr <<"WARNING: AliPHOSGetter::Post TerRef -> //" << fTasksFolder << "/Reconstructioner/PHOS not found!" << endl; 
    return 0 ; 
  }   

  TList * l = phos->GetListOfTasks() ; 
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
    return static_cast<void *>(l->GetObjectRef(tsm)) ;
  else
    return 0 ;
  
} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostRecParticles(const char * name) const 
{  // -------------------- RecParticles ------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/TrackSegments/name

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS")) ; 
  
  if ( !phosFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post RPa -> Folder //" << fRecoFolder << "/PHOS/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post Rpa -> Adding Folder //" << fRecoFolder << "/PHOS/" << endl;
    }
    phosFolder = fRecoFolder->AddFolder("PHOS", "Reconstructed data from PHOS") ;  
  }    

 TFolder * phosRPaFolder  = dynamic_cast<TFolder*>(phosFolder->FindObject("RecParticles")) ;
  if ( !phosRPaFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post RPa -> Folder //" << fRecoFolder << "/PHOS/RecParticles/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post RPa -> Adding Folder //" << fRecoFolder << "/PHOS/RecParticles/" << endl;
    }
    phosRPaFolder = phosFolder->AddFolder("RecParticles", "RecParticles from PHOS") ;  
  } 

  TObject * rps = phosRPaFolder->FindObject( name )  ;
  if ( !rps ) {
    TClonesArray * rp = new TClonesArray("AliPHOSRecParticle",100) ;
    rp->SetName(name) ;    
    phosRPaFolder->Add(rp) ;  
  }
  return kTRUE; 
} 

//____________________________________________________________________________ 
void * AliPHOSGetter::RecParticlesRef(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/TrackSegments/name

 if ( !fRecoFolder ) {
    cerr << "ERROR: AliPHOSGetter::RecParticlesRef -> Folder//" << fRecoFolder << " not found!" << endl; 
    return 0 ; 
  }    

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/RecParticles")) ; 
  if ( !phosFolder ) {
    cerr << "ERROR: AliPHOSGetter::RecParticlesRef -> Folder //" << fRecoFolder << "/PHOS/RecParticles/ not found!" << endl;
    return 0;
  }    

  TObject * tss =  phosFolder->FindObject(name  ) ;
  if (!tss) {
    return 0 ;  
  }
  return static_cast<void *>(phosFolder->GetListOfFolders()->GetObjectRef(tss)) ;
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostPID(AliPHOSPID * pid) const 
{      // ------------AliPHOS PID -----------------------------

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliPHOSGetter::Post Per -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
  
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Per -> //" << fTasksFolder << "/Reconstructioner/PHOS not found!" << endl; 
      cout <<"INFO: AliPHOSGetter::Post Per -> Adding //" << fTasksFolder << "/Reconstructioner/PHOS" << endl;
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  AliPHOSPID * phospid = dynamic_cast<AliPHOSPID*>(phos->GetListOfTasks()->FindObject(pid->GetName())) ; 
  if (phospid) { 
    if (fDebug)
      cout << "INFO: AliPHOSGetter::Post Per -> Task " << pid->GetName()
	   << " already exists" << endl ; 
    phos->GetListOfTasks()->Remove(phospid) ;
  }
  
  phos->Add(pid) ;      
  return kTRUE; 
} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostPID(const char * name) const 
{     
  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliPHOSGetter::Post Per -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
  
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      cout <<"WARNING: AliPHOSGetter::Post Per -> //" << fTasksFolder << "/Reconstructioner/PHOS not found!" << endl; 
      cout <<"INFO: AliPHOSGetter::Post Per -> Adding //" << fTasksFolder << "/Reconstructioner/PHOS" << endl;
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  TList * l = phos->GetListOfTasks() ;   
  TIter it(l) ;
  TString pidname(name) ;
  pidname+=":pid" ; 
  TTask * task ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(pidname))
      return kTRUE ;
  }
 
  AliPHOSPIDv1 * phospid = new AliPHOSPIDv1() ;
  pidname+="-v1" ;
  phospid->SetName(pidname) ;
  phos->Add(phospid) ;      
  
  return kTRUE; 
} 

//____________________________________________________________________________ 
void * AliPHOSGetter::PIDRef(const char * name) const 
{ //------------PID ------------------------------

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    cerr << "ERROR: AliPHOSGetter::Post PerRef -> Task //" << fTasksFolder << "/Reconstructioner not found!" << endl;
    return kFALSE ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    cerr <<"WARNING: AliPHOSGetter::Post PerRef -> //" << fTasksFolder << "/ReconstructionerPHOS not found!" << endl; 
    return 0 ; 
  }   
  
  TList * l = phos->GetListOfTasks() ; 
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
    return static_cast<void *>(l->GetObjectRef(pid)) ;
  else
    return 0 ;
  
} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostQA(void) const 
{ // ------------------ QA ---------------------------------

  // the hierarchy is //Folders/Run/Conditions/QA/PHOS/alarmsName

  TFolder * phosFolder = dynamic_cast<TFolder*>(fQAFolder->FindObject("PHOS")) ; 
  if ( !phosFolder ) {
    if (fDebug) {
      cout << "WARNING: AliPHOSGetter::Post Q -> Folder //" << fQAFolder << "/PHOS/ not found!" << endl;
      cout << "INFO:    AliPHOSGetter::Post Q -> Adding Folder //" << fQAFolder << "/PHOS/" << endl;
    }
    phosFolder = fQAFolder->AddFolder("PHOS", "QA from PHOS") ; 
  }      

  return kTRUE;
}

//____________________________________________________________________________ 
void * AliPHOSGetter::AlarmsRef(void) const 
{  //------- Alarms ----------------------

  
  // the hierarchy is //Folders/Run/Conditions/QA/PHOS
  if ( !fQAFolder ) {
    cerr << "ERROR: AliPHOSGetter::Post QRef -> Folder //" << fQAFolder << " not found!" << endl;
    return 0;
  }    
 
  TFolder * phosFolder = dynamic_cast<TFolder *>(fQAFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) {
    cerr << "ERROR: AliPHOSGetter::Post QRef -> Folder //" << fQAFolder << "/PHOS/ not found!" << endl;
    return 0;
  }
   
  return static_cast<void *>(fQAFolder->GetListOfFolders()->GetObjectRef(phosFolder)) ;
}

//____________________________________________________________________________ 
const TParticle * AliPHOSGetter::Primary(Int_t index) const
{
  // Return primary particle numbered by <index>

  if(index < 0) 
    return 0 ;
  TParticle *  p = gAlice->Particle(index) ; 
  if (p->GetFirstMother() != -1 ) {
    cout << "AliPHOSGetter::Primary : Not a primary " << endl ; 
  }
  return p ; 
  
  
//   Int_t primaryIndex = index % 10000000 ; 
//   Int_t primaryList = (Int_t ) ((index-primaryIndex)/10000000.)  ;
  
//   if ( primaryList > 0  ) {
//     if (fDebug) {
//       cout << " Getter does not support currently Mixing of primary " << endl ;
//       cout << "   can not return primary: " << index<< " (list "<< primaryList<< " primary # " << primaryIndex << " )"<<endl ;
//     }
//     return 0;
//   }
  
//   return gAlice->Particle(primaryIndex) ;
  
}

//____________________________________________________________________________ 
const TParticle * AliPHOSGetter::Secondary(TParticle* p, Int_t index) const
{
  // Return first (index=1) or second (index=2) secondary particle of primary particle p 

  if(index <= 0) 
    return 0 ;
  if(index > 2)
    return 0 ;

  if(p) {
  Int_t daughterIndex = p->GetDaughter(index-1) ; 
  return  gAlice->Particle(daughterIndex) ; 
  }
  else
    return 0 ;
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeD()
{
  // Read the digit tree gAlice->TreeD()  
  if(gAlice->TreeD()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeD: can not read TreeD " << endl ;
  return ;
  }
  
  TObjArray * lob = static_cast<TObjArray*>(gAlice->TreeD()->GetListOfBranches()) ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * digitsbranch = 0 ; 
  TBranch * digitizerbranch = 0 ; 
  Bool_t phosfound = kFALSE, digitizerfound = kFALSE ; 
  
  while ( (branch = static_cast<TBranch*>(next())) && (!phosfound || !digitizerfound) ) {
    if ( (strcmp(branch->GetName(), "PHOS")==0) && (strcmp(branch->GetTitle(), fDigitsTitle)==0) ) {
      digitsbranch = branch ; 
      phosfound = kTRUE ;
    }
    else if ( (strcmp(branch->GetName(), "AliPHOSDigitizer")==0) && (strcmp(branch->GetTitle(), fDigitsTitle)==0) ) {
      digitizerbranch = branch ; 
      digitizerfound = kTRUE ; 
    }
  }

  if ( !phosfound || !digitizerfound ) {
    if (fDebug)
      cout << "WARNING: AliPHOSGetter::ReadTreeD -> Cannot find Digits and/or Digitizer with name " 
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
void AliPHOSGetter::ReadTreeH()
{
  // Read the first entry of PHOS branch in hit tree gAlice->TreeH()

  if(gAlice->TreeH()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeH: -> Cannot read TreeH " << endl ;
    return ;
  }
  
  TBranch * hitsbranch = static_cast<TBranch*>(gAlice->TreeH()->GetBranch("PHOS")) ;
  if ( !hitsbranch ) {
    if (fDebug)
      cout << "WARNING:  AliPHOSGetter::ReadTreeH -> Cannot find branch PHOS" << endl ; 
    return ;
  }
  if(!Hits())
    PostHits() ;

  hitsbranch->SetAddress(HitsRef()) ;

  hitsbranch->GetEntry(0) ;

}

//____________________________________________________________________________ 
void AliPHOSGetter::Track(Int_t itrack)
{
  // Read the first entry of PHOS branch in hit tree gAlice->TreeH()

  if(gAlice->TreeH()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeH: -> Cannot read TreeH " << endl ;
    return ;
  }
  
  TBranch * hitsbranch = dynamic_cast<TBranch*>(gAlice->TreeH()->GetListOfBranches()->FindObject("PHOS")) ;
  if ( !hitsbranch ) {
    if (fDebug)
      cout << "WARNING:  AliPHOSGetter::ReadTreeH -> Cannot find branch PHOS" << endl ; 
    return ;
  }  
  if(!Hits())
    PostHits() ;
  hitsbranch->SetAddress(HitsRef()) ;
  hitsbranch->GetEntry(itrack) ;


}
//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeQA()
{
  // Read the digit tree gAlice->TreeQA()
  // so far only PHOS knows about this Tree  

  if(PHOS()->TreeQA()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeQA: can not read TreeQA " << endl ;
    return ;
  }
  
  TBranch * qabranch = PHOS()->TreeQA()->GetBranch("PHOS") ; 
  if (!qabranch) { 
    if (fDebug)
      cout << "WARNING: AliPHOSGetter::ReadTreeQA -> Cannot find QA Alarms for PHOS" << endl ;
    return ; 
  }   
  
  if(!Alarms())
    PostQA() ; 

  qabranch->SetAddress(AlarmsRef()) ;

  qabranch->GetEntry(0) ;
 
//   PostQA("PHOS") ; 
//   TFolder * alarmsF = Alarms() ; 
//   alarmsF->Clear() ; 
//   qabranch->SetAddress(&alarmsF) ;
//   qabranch->GetEntry(0) ;
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeR()
{
  // Read the reconstrunction tree gAlice->TreeR()

  if(gAlice->TreeR()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeR: can not read TreeR " << endl ;
    return ;
  }
  
  // RecPoints 
  TObjArray * lob = static_cast<TObjArray*>(gAlice->TreeR()->GetListOfBranches()) ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * emcbranch = 0 ; 
  TBranch * cpvbranch = 0 ; 
  TBranch * clusterizerbranch = 0 ; 
  Bool_t phosemcrpfound = kFALSE, phoscpvrpfound = kFALSE, clusterizerfound = kFALSE ; 
  
  while ( (branch = static_cast<TBranch*>(next())) && (!phosemcrpfound || !phoscpvrpfound || !clusterizerfound) ) 
    if(strcmp(branch->GetTitle(), fRecPointsTitle)==0) {
      if ( strcmp(branch->GetName(), "PHOSEmcRP")==0) {
	emcbranch = branch ; 
	phosemcrpfound = kTRUE ;
      }
      else if ( strcmp(branch->GetName(), "PHOSCpvRP")==0) {
	cpvbranch = branch ; 
	phoscpvrpfound = kTRUE ;
      }
      else if(strcmp(branch->GetName(), "AliPHOSClusterizer")==0){
	clusterizerbranch = branch ; 
	clusterizerfound = kTRUE ; 
      }
    }

  if ( !phosemcrpfound ) {
    if (fDebug)
      cout << "WARNING: AliPHOSGetter::ReadTreeR -> Cannot find EmcRecPoints with title " 
	   << fRecPointsTitle << endl ;
    return ; 
  }   
  if ( !phoscpvrpfound ) {
    if (fDebug)
      cout << "WARNING: AliPHOSGetter::ReadTreeR -> Cannot find CpvRecPoints with title " 
	   << fRecPointsTitle << endl ;
    return ; 
  }   
  if ( !clusterizerfound ) {
    if (fDebug)
      cout << "WARNING: AliPHOSGetter::ReadTreeR -> Can not find Clusterizer with title " 
	   << fRecPointsTitle << endl ;
    return ; 
  }   
  
  // Read and Post the RecPoints
  if(!EmcRecPoints(fRecPointsTitle) )
    PostRecPoints(fRecPointsTitle) ;
  emcbranch->SetAddress(EmcRecPointsRef(fRecPointsTitle)) ;
  emcbranch->GetEntry(0) ;

  cpvbranch->SetAddress(CpvRecPointsRef(fRecPointsTitle)) ;
  cpvbranch->GetEntry(0) ;
  
  if(!Clusterizer(fRecPointsTitle) )
    PostClusterizer(fRecPointsTitle) ;
  clusterizerbranch->SetAddress(ClusterizerRef(fRecPointsTitle)) ;
  clusterizerbranch->GetEntry(0) ;
 
  
  //------------------- TrackSegments ---------------------
  next.Reset() ; 
  TBranch * tsbranch = 0 ; 
  TBranch * tsmakerbranch = 0 ; 
  Bool_t phostsfound = kFALSE, tsmakerfound = kFALSE ; 
    
  while ( (branch = static_cast<TBranch*>(next())) && (!phostsfound || !tsmakerfound) ) 
    if(strcmp(branch->GetTitle(), fTrackSegmentsTitle)==0)  {
      if ( strcmp(branch->GetName(), "PHOSTS")==0){
	tsbranch = branch ; 
	phostsfound = kTRUE ;
      }
      else if(strcmp(branch->GetName(), "AliPHOSTrackSegmentMaker")==0) {
	tsmakerbranch = branch ; 
	tsmakerfound  = kTRUE ; 
      }
    }
  
  if ( !phostsfound || !tsmakerfound ) {
    if (fDebug)
      cout << "WARNING: AliPHOSGetter::ReadTreeR -> Cannot find TrackSegments and/or TrackSegmentMaker with name "
	   << fTrackSegmentsTitle << endl ;
    return ; 
  } 
  
  // Read and Post the TrackSegments
  if(!TrackSegments(fTrackSegmentsTitle))
    PostTrackSegments(fTrackSegmentsTitle) ;
  tsbranch->SetAddress(TrackSegmentsRef(fTrackSegmentsTitle)) ;
  tsbranch->GetEntry(0) ;
  
  // Read and Post the TrackSegment Maker
  if(!TrackSegmentMaker(fTrackSegmentsTitle))
    PostTrackSegmentMaker(fTrackSegmentsTitle) ;
  tsmakerbranch->SetAddress(TSMakerRef(fTrackSegmentsTitle)) ;
  tsmakerbranch->GetEntry(0) ;
  
  
  //------------ RecParticles ----------------------------
  next.Reset() ; 
  TBranch * rpabranch = 0 ; 
  TBranch * pidbranch = 0 ; 
  Bool_t phosrpafound = kFALSE, pidfound = kFALSE ; 
  
  while ( (branch = static_cast<TBranch*>(next())) && (!phosrpafound || !pidfound) ) 
    if(strcmp(branch->GetTitle(), fRecParticlesTitle)==0) {   
      if ( strcmp(branch->GetName(), "PHOSRP")==0) {   
	rpabranch = branch ; 
	phosrpafound = kTRUE ;
      }
      else if (strcmp(branch->GetName(), "AliPHOSPID")==0) {
	pidbranch = branch ; 
	pidfound  = kTRUE ; 
      }
    }
  
  if ( !phosrpafound || !pidfound ) {
    if (fDebug)
      cout << "WARNING: AliPHOSGetter::ReadTreeR -> Cannot find RecParticles and/or PID with name " 
	   << fRecParticlesTitle << endl ;
    return ; 
  } 
  
  // Read and Post the RecParticles
  if(!RecParticles(fRecParticlesTitle))
    PostRecParticles(fRecParticlesTitle) ;
  rpabranch->SetAddress(RecParticlesRef(fRecParticlesTitle)) ;
  rpabranch->GetEntry(0) ;
  
  // Read and Post the PID
  if(!PID(fRecParticlesTitle))
    PostPID(fRecParticlesTitle) ;
  pidbranch->SetAddress(PIDRef(fRecParticlesTitle)) ;
  pidbranch->GetEntry(0) ;
  
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeS(Int_t event)
{
  // Read the summable digits tree gAlice->TreeS()  
  
  // loop over all opened files and read their SDigits to the White Board
  TFolder * phosF = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject("PHOS")) ;
  if (!phosF) 
    phosF = fSDigitsFolder->AddFolder("PHOS", "SDigits from PHOS") ; 
  TCollection * folderslist = phosF->GetListOfFolders() ; 
  
  //Add current file to list if it is not there yet
  if ( (fHeaderFile != "aliroot") && ( !folderslist->Contains(fHeaderFile) ) ){
    phosF->AddFolder(fHeaderFile, ""); 
  }
    
  TIter next(folderslist) ; 
  TFolder * folder = 0 ; 
  TFile * file; 
  TTree * treeS = 0;
  while ( (folder = static_cast<TFolder*>(next())) ) {
    if(fHeaderFile.CompareTo(folder->GetName()) == 0 ) 
      treeS=gAlice->TreeS() ;
    else{
      file = static_cast<TFile*>(gROOT->GetFile(folder->GetName())); 
      file->cd() ;
      
      // Get SDigits Tree header from file
      TString treeName("TreeS") ;
      treeName += event ; 
      treeS = dynamic_cast<TTree*>(gDirectory->Get(treeName.Data()));
    }
    if(treeS==0){
      cerr << "ERROR: AliPHOSGetter::ReadTreeS There is no SDigit Tree" << endl;
      return ;
    }
    
    //set address of the SDigits and SDigitizer
    TBranch   * sdigitsBranch    = 0;
    TBranch   * sdigitizerBranch = 0;
    TBranch   * branch           = 0 ;  
    TObjArray * lob = static_cast<TObjArray*>(treeS->GetListOfBranches()) ;
    TIter next(lob) ; 
    Bool_t phosfound = kFALSE, sdigitizerfound = kFALSE ; 
    
    while ( (branch = static_cast<TBranch*>(next())) && (!phosfound || !sdigitizerfound) ) {
      if ( (strcmp(branch->GetName(), "PHOS")==0) && (strcmp(branch->GetTitle(), fSDigitsTitle)==0) ) {
	phosfound = kTRUE ;
	sdigitsBranch = branch ; 
      }
      
      else if ( (strcmp(branch->GetName(), "AliPHOSSDigitizer")==0) && (strcmp(branch->GetTitle(), fSDigitsTitle)==0) ) {
	sdigitizerfound = kTRUE ; 
	sdigitizerBranch = branch ;
      }
    }
    if ( !phosfound || !sdigitizerfound ) {
      if (fDebug)
	cout << "WARNING: AliPHOSDigitizer::ReadSDigits -> Digits and/or Digitizer branch with name " << GetName() 
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
    file   = static_cast<TFile*>(gROOT->GetFile(folder->GetName())); 
    file   ->cd() ;
  }
  
}
//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeS(TTree * treeS, Int_t input)
{  // Read the summable digits fron treeS()  


  TString filename("mergefile") ;
  filename+= input ;

  TFolder * phosFolder = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("PHOS")) ; 
  if ( !phosFolder ) { 
   phosFolder = fSDigitsFolder->AddFolder("PHOS", "SDigits from PHOS") ; 
  } 
  TFolder * folder=(TFolder*)phosFolder->FindObject(filename) ;
  //set address of the SDigits and SDigitizer
  TBranch   * sdigitsBranch    = 0;
  TBranch   * sdigitizerBranch = 0;
  TBranch   * branch           = 0 ;  
  TObjArray * lob = (TObjArray*)treeS->GetListOfBranches() ;
  TIter next(lob) ; 
  Bool_t phosfound = kFALSE, sdigitizerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phosfound || !sdigitizerfound) ) {
    if ( strcmp(branch->GetName(), "PHOS")==0) {
      phosfound = kTRUE ;
      sdigitsBranch = branch ; 
    }
    
    else if ( strcmp(branch->GetName(), "AliPHOSSDigitizer")==0) {
      sdigitizerfound = kTRUE ; 
      sdigitizerBranch = branch ;
    }
  }
  if ( !phosfound || !sdigitizerfound ) {
    if (fDebug)
      cout << "WARNING: AliPHOSGetter::ReadTreeS -> Digits and/or Digitizer branch not found" << endl ;
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
void AliPHOSGetter::ReadPrimaries()
{

  
  // Read the first entry of PHOS branch in hit tree gAlice->TreeH()
  TTree * particleTree = gAlice->TreeK()  ; 
  if( ! particleTree ){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreePrimaries: -> Cannot read TreeK " << endl ;
    return ;
  }
  
  TBranch * particlesbranch = static_cast<TBranch*>(particleTree->GetBranch("Particles")) ;
  if ( !particlesbranch ) {
    if (fDebug)
      cout << "WARNING:  AliPHOSGetter::ReadTreePrimaries -> Cannot find branch Particles" << endl ; 
    return ;
  }

  TParticle * particle = 0; 
  TClonesArray * ar ; 
  if(! (ar = Primaries()) ) { 
    PostPrimaries() ;
    ar = Primaries() ; 
  }
  ar->Delete() ; 
  particlesbranch->SetAddress(&particle) ;
  
  Int_t index = 0 ; 
  for (index = 0 ; index < particleTree->GetEntries(); index++) { 
    particlesbranch->GetEntry(index) ;
    new ((*ar)[index]) TParticle(*particle);
  }
  fNPrimaries= ar->GetEntries() ; 
  
  // Reads specific branches of primaries
  // fNPrimaries = gAlice->GetNtrack();
  
  //Check, is it necessary to open new files
 //  TArrayI* events = fDigitizer->GetCurrentEvents() ; 
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
//     char treeName[20];
//     sprintf(treeName,"TreeK%d",events->At(input));
//     TTree * treeK = (TTree*)gDirectory->Get(treeName);
//     if (treeK) 
//       treeK->SetBranchAddress("Particles", &fParticleBuffer);
//     else    
//       cout << "AliPHOSGetter: cannot find Kine Tree for event:" << events->At(input) << endl;

//     // Create the particle stack
//     if(!fParticles) fParticles = new TClonesArray("TParticle",1000);
//     // Build the pointer list
//     if(fParticleMap) {     <----
//       fParticleMap->Clear();
//       fParticleMap->Expand(treeK->GetEntries());
//     } else
//       fParticleMap = new TObjArray(treeK->GetEntries());
    
//     // From gAlice->Particle(i) 


//   if(!(*fParticleMap)[i]) {
//     Int_t nentries = fParticles->GetEntries();
    
//     // algorithmic way of getting entry index
//     // (primary particles are filled after secondaries)
//     Int_t entry;
//     if (i<fHeader.GetNprimary())
//       entry = i+fHeader.GetNsecondary();
//     else 
//       entry = i-fHeader.GetNprimary();
      
//     // only check the algorithmic way and give
//     // the fatal error if it is wrong
//     if (entry != fParticleFileMap[i]) {
//       Fatal("Particle",
//         "!!!! The algorithmic way is WRONG: !!!\n entry: %d map: %d",
// 	entry, fParticleFileMap[i]); 
//     }  
      
//     fTreeK->GetEntry(fParticleFileMap[i]);
//     new ((*fParticles)[nentries]) TParticle(*fParticleBuffer);
//     fParticleMap->AddAt((*fParticles)[nentries],i);
//   }
//   return (TParticle *) (*fParticleMap)[i];

   
    
// //   }


// //   //scan over opened files and read corresponding TreeK##

//   return ;
}
//____________________________________________________________________________ 
void AliPHOSGetter::Event(const Int_t event, const char* opt)
{
  // Reads the content of all Tree's S, D and R
  
  if (event >= gAlice->TreeE()->GetEntries() ) {
    cerr << "ERROR: AliPHOSGetter::Event -> " << event << " not found in TreeE!" << endl ; 
    return ; 
  }
  gAlice->GetEvent(event) ; 
 
  if(strstr(opt,"H") )
    ReadTreeH() ;
  
  if(strstr(opt,"S") )
    ReadTreeS(event) ;

  if( strstr(opt,"D") )
    ReadTreeD() ;

  if( strstr(opt,"R") )
    ReadTreeR() ;

  if( strstr(opt,"Q") )
    ReadTreeQA() ;

  //  if( strstr(opt,"P") || (strcmp(opt,"")==0) )
  //  ReadPrimaries() ;

}

//____________________________________________________________________________ 
TObject * AliPHOSGetter::ReturnO(TString what, TString name, TString file) const 
{
  // get the object named "what" from the folder
  // folders are named like //Folders

  if ( file.IsNull() ) 
    file = fHeaderFile ; 

  TFolder * folder = 0 ;
  TObject * phosO  = 0 ; 

  //  if ( name.IsNull() ) {
  if ( what.CompareTo("Primaries") == 0 ) {
    folder = dynamic_cast<TFolder *>(fPrimariesFolder->FindObject("Primaries")) ; 
    if (folder) 
      phosO  = dynamic_cast<TObject *>(folder->FindObject("Primaries")) ;  
  }
  else if ( what.CompareTo("Hits") == 0 ) {
    folder = dynamic_cast<TFolder *>(fHitsFolder->FindObject("PHOS")) ; 
    if (folder) 
      phosO  = dynamic_cast<TObject *>(folder->FindObject("Hits")) ;  
  }
  else if ( what.CompareTo("SDigits") == 0 ) { 
    TString path = "PHOS/" + file  ; 
    folder = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject(path.Data())) ; 
    if (folder) { 
      if (name.IsNull())
	name = fSDigitsTitle ; 
      phosO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    }
  }
  else if ( what.CompareTo("Digits") == 0 ){
    folder = dynamic_cast<TFolder *>(fDigitsFolder->FindObject("PHOS")) ; 
    if (folder) { 
      if (name.IsNull())
	name = fDigitsTitle ; 
      phosO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    } 
  }
  else if ( what.CompareTo("EmcRecPoints") == 0 ) {
    folder = dynamic_cast<TFolder *>(fRecoFolder->FindObject("PHOS/EMCARecPoints")) ; 
    if (folder) { 
      if (name.IsNull())
	name = fRecPointsTitle ; 
      phosO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    } 
  }
  else if ( what.CompareTo("CpvRecPoints") == 0 ) {
    folder = dynamic_cast<TFolder *>(fRecoFolder->FindObject("PHOS/CPVRecPoints")) ; 
    if (folder) { 
      if (name.IsNull())
	name = fRecPointsTitle ; 
      phosO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    }   
  }
  else if ( what.CompareTo("TrackSegments") == 0 ) {
    folder = dynamic_cast<TFolder *>(fRecoFolder->FindObject("PHOS/TrackSegments")) ; 
    if (folder) { 
      if (name.IsNull())
	name = fTrackSegmentsTitle ; 
      phosO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    }   
  }
  else if ( what.CompareTo("RecParticles") == 0 ) {
    folder = dynamic_cast<TFolder *>(fRecoFolder->FindObject("PHOS/RecParticles")) ; 
   if (folder) { 
      if (name.IsNull())
	name = fRecParticlesTitle ; 
      phosO  = dynamic_cast<TObject *>(folder->FindObject(name)) ; 
    }   
 }
  else if ( what.CompareTo("Alarms") == 0 ){ 
    if (name.IsNull() ) 
      phosO = dynamic_cast<TObject *>(fQAFolder->FindObject("PHOS")) ;  
    else {
      folder = dynamic_cast<TFolder *>(fQAFolder->FindObject("PHOS")) ; 
      if (!folder) 
	phosO = 0 ; 
      else 
	phosO = dynamic_cast<TObject *>(folder->FindObject(name)) ;  
    }
  }
  if (!phosO) {
    if(fDebug)
      cerr << "ERROR : AliPHOSGetter::ReturnO -> Object " << what << " not found in " << folder->GetName() << endl ; 
    return 0 ;
  }

  return phosO ;
}
  
//____________________________________________________________________________ 
const TTask * AliPHOSGetter::ReturnT(TString what, TString name) const 
{
  // get the TTask named "what" from the folder
  // folders are named like //Folders/Tasks/what/PHOS/name

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
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << what << " not found!" << endl ;  
    return 0 ; 
  }

  TTask * phosT = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if (!phosT) { 
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << what << "/PHOS not found!" << endl ;  
    return 0 ; 
  }
  
  TList * list = phosT->GetListOfTasks() ; 
 
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
  else  if (what.CompareTo("TrackSegmentMaker") == 0){ 
    if ( name.IsNull() )
      name =  fTrackSegmentsTitle ;
    name.Append(":tsm") ;
  }
  else  if (what.CompareTo("PID") == 0){ 
    if ( name.IsNull() )
      name =  fRecParticlesTitle ;
    name.Append(":pid") ;
  }
  else  if (what.CompareTo("QATasks") == 0){ 
    if ( name.IsNull() )
      return phosT ;
  }
  
  TIter it(list) ;
  TTask * task = 0 ; 
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(name))
      return task ;
  }
  
  if(fDebug)
    cout << "WARNING: AliPHOSGetter::ReturnT -> Task " << search << "/" << name << " not found!" << endl ; 
  return 0 ;
}
