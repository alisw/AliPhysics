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
/* $Log:
   08.2002 Dmitri Peressounko:

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
//  gime->Event(event) ;    // reads new event from galice.root
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

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliConfig.h"
#include "AliPHOSGetter.h"
#include "AliPHOS.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSPIDv1.h" 
#include "AliPHOSGeometry.h"
#include "AliPHOSRaw2Digits.h"
//#include "AliPHOSCalibrationDB.h"
#include "AliPHOSBeamTestEvent.h"
ClassImp(AliPHOSGetter)
  
  AliPHOSGetter * AliPHOSGetter::fgObjGetter = 0 ; 
  TFile * AliPHOSGetter::fFile = 0 ; 

//____________________________________________________________________________ 
AliPHOSGetter::AliPHOSGetter(const char* headerFile, const char* branchTitle, const Bool_t toSplit )
{
  // This is the ctor called by GetInstance and the only one that can be used 

  if( fHeaderFile.Contains("_") ) {
    Fatal("AliPHOSGetter", "Invalid file name (_ not allowed) %s", fHeaderFile.Data() ) ;
  }

  //Initialize  all data

  fFailed = kFALSE ;   
  fDebug  = 0 ; 
  fAlice  = 0 ; 
  fBTE    = 0 ;

  fToSplit    = toSplit ;
  fHeaderFile = headerFile ; 

  fPrimaries = new TObjArray(1) ;

  fModuleFolder    = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Configuration/Modules")); 
  fPrimariesFolder = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data")); 
  fHitsFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/Hits")); 
  fSDigitsFolder   = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/SDigits")); 
  fDigitsFolder    = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Event/Data")); 
  fRecoFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Event/RecData")); 
  fQAFolder        = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Conditions/QA")); 
  fTasksFolder     = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Tasks")) ; 

  //Set titles to branches and create PHOS specific folders
  SetTitle(branchTitle) ;

  if ( fHeaderFile != "aliroot"  ) { // to call the getter without a file
    //open headers file
    fFile = static_cast<TFile*>(gROOT->GetFile(fHeaderFile.Data() ) ) ;

    if(!fFile) {    //if file was not opened yet, read gAlice
      fFile = TFile::Open(fHeaderFile.Data(),"update") ;
      if (!fFile->IsOpen()) {
	Error("AliPHOSGetter", "Cannot open %s", fHeaderFile.Data() ) ; 
	fFailed = kTRUE ;
        return ;  
      }
    }
    gAlice = dynamic_cast<AliRun *>(fFile->Get("gAlice")) ;
  }
  
  if (!gAlice) {
    Error("AliPHOSGetter", "Cannot find gAlice in %s", fHeaderFile.Data() ) ; 
    fFailed = kTRUE ;
    return ; 
  }
  if (!PHOS()) {
    if (fDebug)
      Info("AliPHOSGetter", "-> Posting PHOS to Folders") ; 
    if (gAlice->GetDetector("PHOS")) {
      AliConfig * conf = AliConfig::Instance() ; 
      conf->Add(static_cast<AliDetector*>(gAlice->GetDetector("PHOS"))) ; 
      conf->Add(static_cast<AliModule*>(gAlice->GetDetector("PHOS"))) ; 
    }
    else 
      Error("AliPHOSGetter", "detector PHOS not found") ;    
  }

  fDebug=0;
}

//____________________________________________________________________________ 
AliPHOSGetter::~AliPHOSGetter(){

  if (fPrimaries) {
    fPrimaries->Delete() ; 
    delete fPrimaries ; 
  }

  TFolder * phosF = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject("PHOS")) ;
  TCollection * folderslist = phosF->GetListOfFolders() ; 
  TIter next(folderslist) ; 
  TFolder * folder = 0 ; 
  while ( (folder = static_cast<TFolder*>(next())) ) 
    phosF->Remove(folder) ; 
  
  if (fFile) { 
    fFile->Close() ;  
    delete fFile ;
    fFile = 0 ;
  }
  fgObjGetter = 0 ; 
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::CloseFile()
{
  if(gAlice)
    delete gAlice ;  
  gAlice = 0 ; 
  if(fAlice)
    delete fAlice ; 
  fAlice = 0 ; 
}

//____________________________________________________________________________ 
const TFolder * AliPHOSGetter::Folder(const TString what) const {

  // returns the PHOS folder required by what
  // what = hits, sdigits, digits

  if ( what == "hits" ) 
    return dynamic_cast<const TFolder *>(fHitsFolder->FindObject("PHOS")) ; 
  else if ( what == "sdigits" ) 
    return  dynamic_cast<const TFolder *>(fSDigitsFolder->FindObject("PHOS")) ; 
  else if ( what == "digits" ) 
    return  dynamic_cast<const TFolder *>(fDigitsFolder->FindObject("PHOS")) ; 
  else {
    Error("GetFolder", " %s illegal option (hits, sdigits, digits)", what.Data() ) ; 
    return 0 ; 
  }
}

//____________________________________________________________________________ 
AliPHOSGetter * AliPHOSGetter::GetInstance()
{
  // Returns the pointer of the unique instance already defined
  
  if ( fgObjGetter ) {
    return fgObjGetter ;
  }
  else {
    //Warning("GetInstance", "not yet initialized") ;
    return 0 ; 
  }
}

//____________________________________________________________________________ 
AliPHOSGetter * AliPHOSGetter::GetInstance(const char* headerFile,
					   const char* branchTitle,
                                           const Bool_t toSplit)
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed

  if(!fgObjGetter){
    fgObjGetter = new AliPHOSGetter(headerFile,branchTitle,toSplit) ;
    if(fgObjGetter->fFailed)
      return 0;
    else
      return fgObjGetter ;
  }

  //First checks, if header file already opened
  if(!fgObjGetter->fFile){
     fgObjGetter = new AliPHOSGetter(headerFile,branchTitle,toSplit) ;
    if(fgObjGetter->fFailed)
      return 0;
    else
      return fgObjGetter ;
  }

  if(fgObjGetter->fHeaderFile.CompareTo(headerFile)==0){ //Opened the same header file
    if((fgObjGetter->fBranchTitle.CompareTo(branchTitle) == 0)&&   //Open the same branch title
       (toSplit==fgObjGetter->fToSplit)){                          //Nothing should be cleaned
    }
    else{ //Clean all data and AliPHOS...zers
      if(fgObjGetter->fToSplit)
	fgObjGetter->CloseSplitFiles() ;
      fgObjGetter->CleanWhiteBoard() ;
      fgObjGetter->fToSplit = toSplit ;
      fgObjGetter->SetTitle(branchTitle) ;
    }
  }
  else{  //Close already opened files, clean memory and open new header file
    if(gAlice){ //should first delete gAlice, then close file
      //Should be in dtor of PHOS, but if one changes path ...
      fgObjGetter->fModuleFolder->Remove(fgObjGetter->fModuleFolder->FindObject("PHOS")) ; 
      delete gAlice ;      
    }
    if(fgObjGetter->fFile){
      fgObjGetter->fFile->Close() ;
      fgObjGetter->fFile=0;
    }
    if(fgObjGetter->fToSplit)
      fgObjGetter->CloseSplitFiles() ;
    fgObjGetter->CleanWhiteBoard() ;    
    fgObjGetter = new AliPHOSGetter(headerFile,branchTitle,toSplit) ;

  }
  return fgObjGetter ; 
  
}

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::BranchExists(const TString recName) const
{
  //Looks in the tree Tree"name" if branch with current name olready exists

  TString filename("") ;
  TString name, dataname, zername ;
  if(recName == "SDigits"){
    filename=fSDigitsFileName ;
    name = "TreeS0" ;
    dataname = "PHOS" ;
    zername = "AliPHOSSDigitizer" ;
  }
  else
    if(recName == "Digits"){
      filename=fDigitsFileName ;
      name = "TreeD0" ;
      dataname = "PHOS" ;
      zername = "AliPHOSDigitizer" ;
    }
    else
      if(recName == "RecPoints"){
	filename=fRecPointsFileName ;
	name = "TreeR0" ;
	dataname = "PHOSEmcRP" ;
	zername = "AliPHOSClusterizer" ;
      }
      else
	if(recName == "TrackSegments"){
	  filename=fTrackSegmentsFileName ;
	  name = "TreeR0" ;
	  dataname = "PHOSTS" ;
	  zername = "AliPHOSTrackSegmentMaker" ;
	}	 
	else
	  if(recName == "RecParticles"){
	    filename= fRecParticlesFileName ;
	    name = "TreeR0" ;
	    dataname = "PHOSRP" ;
	    zername = "AliPHOSPID" ;
	  }
	  else
	    return kFALSE ;

  TFile * file ;
  TTree * tree ;
  if(fToSplit){
    file = static_cast<TFile*>(gROOT->GetFile(filename.Data() ) ) ;
    if(!file)
      file = TFile::Open(fSDigitsFileName.Data(),"update");
  }
  else
    file = fFile ;

  tree = (TTree *)file->Get(name.Data()) ;
  
  if(!tree ) 
    return kFALSE ;

  TObjArray * lob = static_cast<TObjArray*>(tree->GetListOfBranches()) ;
  TIter next(lob) ; 
  TBranch * branch = 0 ;  
  TString titleName(fBranchTitle);
  titleName+=":";
  while ((branch = (static_cast<TBranch*>(next())))) {
    TString branchName(branch->GetName() ) ; 
    TString branchTitle(branch->GetTitle() ) ;  
    if ( branchName.BeginsWith(dataname) && branchTitle.BeginsWith(fBranchTitle) ){  
      Warning("BranchExists", "branch %s with title %s already exists in %s", dataname.Data(), fBranchTitle.Data(), name.Data() ) ;
      return kTRUE ;
    }
    if ( branchName.BeginsWith(zername) &&  branchTitle.BeginsWith(titleName) ){
      Warning("BranchExists", "branch AliPHOS... with title %s already exists in %s", branch->GetTitle(), name.Data() ) ;     
      return kTRUE ; 
    }
  }
  //We can't delete three if gAlice points to it... To be redisigned somehow???!!!
  if(!fToSplit){
    if(name.Contains("TreeS"))
      if(tree!=gAlice->TreeS())
	tree->Delete();
    if(name.Contains("TreeD"))
      if(tree!=gAlice->TreeD())
	tree->Delete();
    if(name.Contains("TreeR"))
      if(tree!=gAlice->TreeR())
	tree->Delete();    
  }
  return kFALSE ;
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::ListBranches(Int_t event) const  
{
  
  TBranch * branch = 0 ; 
  if (gAlice->GetEvent(event) == -1)
    return ; 
  
  TTree * t =  gAlice->TreeH() ; 
  if(t){
    Info("ListBranches", "-> ****** Hits    : ") ; 
    TObjArray * lob = t->GetListOfBranches() ;
    TIter next(lob) ; 
    while ( (branch = static_cast<TBranch*>(next())) )
      Info("ListBranches", "             %s", branch->GetName()) ; 
  } else
    Warning("ListBranches", "TreeH not found for event %d", event ) ;  
  
  t = gAlice->TreeS() ;
  if(t){
    Info("ListBranches", "-> ****** SDigits : ") ; 
    TObjArray * lob = t->GetListOfBranches() ;
    TIter next(lob) ; 
    while ( (branch = static_cast<TBranch*>(next())) )
      Info("ListBranches", "             %s %s", branch->GetName(), branch->GetTitle()) ; 
  } else 
    Warning("ListBranches", "TreeS not found for event %d", event)  ;  
  
  
  t = gAlice->TreeD() ;
  if(t){
    Info("ListBranches", "-> ****** Digits  : ") ; 
    TObjArray * lob = t->GetListOfBranches() ;
    TIter next(lob) ; 
    while ( (branch = static_cast<TBranch*>(next())) )
      Info("ListBranches", "             %s %s", branch->GetName(), branch->GetTitle()) ; 
  } else 
    Warning("ListBranches", "TreeD not found for event %d", event) ;  
  

  t = gAlice->TreeR() ;
  if(t){
    Info("ListBranches", "-> ****** Recon   : ") ; 
    TObjArray * lob = t->GetListOfBranches() ;
    TIter next(lob) ; 
    while ( (branch = static_cast<TBranch*>(next())) )
     Info("ListBranches", "             %s %s", branch->GetName(), branch->GetTitle()) ; 
  } else 
    Warning("ListBranches", "TreeR not found for event %d", event) ;  
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::NewBranch(TString name, Int_t event)  
{
  fBranchTitle = fSDigitsTitle = fDigitsTitle = fRecPointsTitle = fTrackSegmentsTitle = fRecParticlesTitle =  name ; 
  Event(event) ; 
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::NewFile(TString name)  
{
  fHeaderFile = name ; 
  fFile->Close() ; 
  fFailed = kFALSE; 

  fFile = static_cast<TFile*>(gROOT->GetFile(fHeaderFile.Data() ) ) ;
  if(!fFile) {    //if file was not opened yet, read gAlice
    fFile = TFile::Open(fHeaderFile.Data(),"update") ;
    if (!fFile->IsOpen()) {
      Error("NewFile", "Cannot open %s", fHeaderFile.Data() ) ; 
      fFailed = kTRUE ;
      return fFailed ;  
    }
    gAlice = static_cast<AliRun *>(fFile->Get("gAlice")) ;
  } 
  
  if (!gAlice) {
    Error("AliPHOSGetter", "Cannot find gAlice in %s", fHeaderFile.Data() ) ; 
    fFailed = kTRUE ;
    return fFailed ; 
  }
  return fFailed ; 
}

//____________________________________________________________________________ 
const AliPHOS * AliPHOSGetter::PHOS() 
{
  // returns the PHOS object 
  AliPHOS * phos = dynamic_cast<AliPHOS*>(fModuleFolder->FindObject("PHOS")) ;  
  if (!phos) 
    if (fDebug)
      Warning("PHOS", "-> PHOS module not found in Folders") ; 
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
const Bool_t AliPHOSGetter::PostPrimaries(void) const 
{  //------- Primaries ----------------------

  // the hierarchy is //Folders/RunMC/Event/Data/Primaries
  
  TFolder * primariesFolder = dynamic_cast<TFolder*>(fPrimariesFolder->FindObject("Primaries")) ; 
  if ( !primariesFolder ) {
    if (fDebug) {
      Warning("PostPrimaries", "-> Folder //%s/Primaries/ not found!", fPrimariesFolder->GetName()) ;
      Info("PostPrimaries", "-> Adding Folder //%s/Primaries", fPrimariesFolder->GetName()) ;
    }
    primariesFolder = fPrimariesFolder->AddFolder("Primaries", "Primaries particles from TreeK") ; 
  }    
  TClonesArray *primaries=  new TClonesArray("TParticle",1000) ;
  primaries->SetName("Primaries") ;
  primariesFolder->Add(primaries) ; 
  
  return kTRUE;
} 

//____________________________________________________________________________ 
TObject** AliPHOSGetter::PrimariesRef(void) const 
{  //------- Primaries ----------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/Primaries
  if ( !fPrimariesFolder ) {
    Fatal("PrimariesRef", "Folder //%s not found", fPrimariesFolder) ;
  }    
 
  TFolder * primariesFolder = dynamic_cast<TFolder *>(fPrimariesFolder->FindObject("Primaries")) ;
  if ( !primariesFolder ) {
    Fatal("PrimariesRef", "Folder //%s/Primaries/ not found", fPrimariesFolder) ;  
  }
 
  TObject * p = primariesFolder->FindObject("Primaries") ;
  if(!p) {
    Fatal("PrimariesRef","%s /Primaries not found !", primariesFolder->GetName() ) ; 
  }

  return primariesFolder->GetListOfFolders()->GetObjectRef(p) ;
}

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostHits(void) const 
{  //------- Hits ----------------------

  // the hierarchy is //Folders/RunMC/Event/Data/PHOS/Hits
  
  TFolder * phosFolder = dynamic_cast<TFolder*>(fHitsFolder->FindObject("PHOS")) ; 
  if ( !phosFolder ) {
    if (fDebug) {
      Warning("PostHits", "-> Folder //%s/PHOS/ not found!", fHitsFolder) ;
      Info("PostHits", "-> Adding Folder //%s/PHOS/", fHitsFolder) ;
    }
    phosFolder = fHitsFolder->AddFolder("PHOS", "Hits from PHOS") ; 
  }    
  TClonesArray *hits=  new TClonesArray("AliPHOSHit",1000) ;
  hits->SetName("Hits") ;
  phosFolder->Add(hits) ; 
  
  return kTRUE;
} 

//____________________________________________________________________________ 
TObject** AliPHOSGetter::HitsRef(void) const 
{  //------- Hits ----------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/PHOS/Hits
  if ( !fHitsFolder ) {
    Fatal("HitsRef", "Folder //%s not found !", fHitsFolder) ;
  }    
 
  TFolder * phosFolder = dynamic_cast<TFolder *>(fHitsFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) {
    Fatal("HitsRef", "Folder //%s/PHOS/ not found !", fHitsFolder) ;  
  }
 
  TObject * h = phosFolder->FindObject("Hits") ;
  if(!h) {
    Fatal("HitsRef", "%s/Hits not fount !", phosFolder->GetName() ) ; 
  }

  return phosFolder->GetListOfFolders()->GetObjectRef(h) ;
}

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostSDigits(const char * name, const char * headerFile) const 
{  //---------- SDigits -------------------------

  
  // the hierarchy is //Folders/RunMC/Event/Data/PHOS/SDigits/headerFile/sdigitsname
  // because you can have sdigits from several hit files for mixing
  
  TFolder * phosFolder = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) {
    if (fDebug) {
      Warning("PostSDigits", "-> Folder //%s/PHOS/ not found!", fSDigitsFolder) ;
      Info("PostSDigits", "-> Adding Folder //%s/PHOS/", fHitsFolder) ;
    }
    phosFolder = fSDigitsFolder->AddFolder("PHOS", "SDigits from PHOS") ; 
  }    

  
  TString subdir(headerFile) ;
  if(fToSplit){
    subdir.Remove(subdir.Last('/')+1,subdir.Length()) ;
    subdir.ReplaceAll("/","_") ; 
    subdir+="PHOS.SDigits." ;
    if(name && (strcmp(name,"Default")!=0)){
      subdir+=name ;
      subdir+="." ;
    }
    subdir+="root" ;
  }
    
  TFolder * phosSubFolder = dynamic_cast<TFolder*>(phosFolder->FindObject(subdir)) ; 
  if ( !phosSubFolder ) 
    phosSubFolder = phosFolder->AddFolder(subdir, ""); 
  

  TObject * sd  = phosSubFolder->FindObject(name); 
  if ( !sd ) {
    TClonesArray * sdigits = new TClonesArray("AliPHOSDigit",1) ;
    sdigits->SetName(name) ;
    phosSubFolder->Add(sdigits) ;
  }
  
  return kTRUE;
} 
//____________________________________________________________________________ 
TObject** AliPHOSGetter::SDigitsRef(const char * name, const char * foldername) const 
{  //------- SDigits ----------------------
  
  // the hierarchy is //Folders/RunMC/Event/Data/PHOS/SDigits/filename/SDigits

  if ( !fSDigitsFolder ) {
    Fatal("SDigitsRef", "Folder //%s not found !", fSDigitsFolder) ;
  }    
 
  TFolder * phosFolder = static_cast<TFolder *>(fSDigitsFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) {
    Fatal("SDigitsRef", "Folder //%s/PHOS not found !", fSDigitsFolder) ;
  }

  TFolder * phosSubFolder = 0 ;

  phosSubFolder = dynamic_cast<TFolder *>(phosFolder->FindObject(foldername)) ;
  
  if(!phosSubFolder) {
    Fatal("SDigitsRef", "Folder //Folders/RunMC/Event/Data/PHOS/%s not found !", foldername) ;
  }

  TObject * dis = phosSubFolder->FindObject(name) ;
  if(!dis){
    Fatal("SDigitsRef", "object %s not found !", name) ; 
  }

  return phosSubFolder->GetListOfFolders()->GetObjectRef(dis) ;

}

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostSDigitizer(AliPHOSSDigitizer * sdigitizer) const 
{  //---------- SDigitizer -------------------------
    
  // the hierarchy is //Folders/Tasks/SDigitizer/PHOS/sdigitsname


  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ; 

  if ( !sd ) {
    Error("PostDigitizer", "Task //%s/SDigitizer not found !", fTasksFolder) ;
    return kFALSE ;
  }        
  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Warning("PostSDigitizer", "->//%s/SDigitizer/PHOS/ not found!", fTasksFolder) ;  
      Info("PostSDigitizer", "-> Adding //%s/SDigitizer/PHOS/", fTasksFolder) ;
    }
    phos = new TTask("PHOS", "") ; 
    sd->Add(phos) ; 
  } 
  AliPHOSSDigitizer * phossd  = dynamic_cast<AliPHOSSDigitizer *>(phos->GetListOfTasks()->FindObject( sdigitizer->GetName() )); 
  if (phossd) { 
    if (fDebug)
      Info("PostSDigitizer", "-> Task %s already exists", sdigitizer->GetName()) ; 
    phos->GetListOfTasks()->Remove(phossd) ;
  }
  phos->Add(sdigitizer) ;	
  return kTRUE; 
  
}

//____________________________________________________________________________ 
TObject** AliPHOSGetter::SDigitizerRef(const char * name) const 
{  

  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ; 
  if ( !sd ) {
    Fatal("SDigitizerRef", "Task //%s/SDigitizer not found !", fTasksFolder) ;
  }        

  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    Fatal("SDigitizerRef", "//%s/SDigitizer/PHOS not found !", fTasksFolder) ;
  }        

  TTask * task = dynamic_cast<TTask*>(phos->GetListOfTasks()->FindObject(name)) ; 

  return phos->GetListOfTasks()->GetObjectRef(task) ;

}

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostSDigitizer(const char * name, const char * file) const 
{  //---------- SDigitizer -------------------------
  
 // the hierarchy is //Folders/Tasks/SDigitizer/PHOS/sdigitsname

  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ; 
  if ( !sd ) {
    Error("PostSDigitizer", "Task //%s/SDigitizer not found !", fTasksFolder) ;
    return kFALSE ;
  }        

  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Error("PostSDigitizer", "->  //%s/SDigitizer/PHOS/ not found!", fTasksFolder) ;
      Info("PostSDigitizer", "-> Adding  //%s/SDigitizer/PHOS", fTasksFolder) ;
    }
    phos = new TTask("PHOS", "") ; 
    sd->Add(phos) ; 
  } 

  TString sdname(name) ;
  sdname.Append(":") ;
  sdname.Append(file);
  sdname.ReplaceAll("/","_") ; 
  AliPHOSSDigitizer * phossd  = dynamic_cast<AliPHOSSDigitizer *>(phos->GetListOfTasks()->FindObject( sdname )); 
  if (!phossd) {
    phossd = new AliPHOSSDigitizer() ;  
    //Note, we can not call constructor with parameters: it will call Getter and scew up everething
    phossd->SetName(sdname) ;
    phossd->SetTitle(file) ;
    phos->Add(phossd) ;	
  }
  return kTRUE; 
  
}

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostDigits(const char * name) const 
{  //---------- Digits -------------------------

  // the hierarchy is //Folders/Run/Event/Data/PHOS/SDigits/name

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("PHOS")) ;

  if ( !phosFolder ) {
    if (fDebug) {
      Warning("PostDigitizer", "-> Folder //%s/PHOS/ not found!", fDigitsFolder) ;
      Info("PostDigitizer", "-> Adding Folder //%s/PHOS/", fDigitsFolder) ;
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
TObject** AliPHOSGetter::DigitsRef(const char * name) const 
{ //------- Digits ----------------------
  
  // the hierarchy is //Folders/Run/Event/Data/PHOS/Digits/name

  if ( !fDigitsFolder ) {
    Fatal("DigitsRef", "Folder //%s not found !", fDigitsFolder) ;
  }    
  
  TFolder * phosFolder  = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("PHOS")) ; 
  if ( !phosFolder ) {
    Fatal("DigitsRef", "Folder //%s/PHOS/ not found !", fDigitsFolder) ;
  }    

  TObject * d = phosFolder->FindObject(name) ;
  if(!d) {
    Fatal("DigitsRef", "object %s not found !", name) ; 
  }

  return phosFolder->GetListOfFolders()->GetObjectRef(d) ;

}

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostDigitizer(AliPHOSDigitizer * digitizer) const 
{  //---------- Digitizer -------------------------
  
  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 

  if ( !sd ) {
    Error("PostDigitizer", "Task //%s/Digitizer not found !", fTasksFolder) ;
    return kFALSE ;
  }        
  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Error("PostDigitizer", "//%s/Digitizer/PHOS not found!", fTasksFolder) ;
      Info("PostDigitizer", "Adding //%s/Digitizer/PHOS", fTasksFolder) ; 
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
const Bool_t AliPHOSGetter::PostDigitizer(const char * name) const 
{  //---------- Digitizer -------------------------
  
 // the hierarchy is //Folders/Tasks/SDigitizer/PHOS/sdigitsname

  TTask * d  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 
  if ( !d ) {
    Error("PostDigitizer", "Task //%s/Digitizer not found !", fTasksFolder) ;
    return kFALSE ;
  }        

  TTask * phos = dynamic_cast<TTask*>(d->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Warning("PostDigitizer", "//%s/Digitizer/PHOS not found!", fTasksFolder) ; 
      Info("PostDigitizer", "Adding //%s/Digitizer/PHOS", fTasksFolder) ;
    }
    phos = new TTask("PHOS", "") ; 
    d->Add(phos) ; 
} 

  TTask * phosd = dynamic_cast<TTask*>(phos->GetListOfTasks()->FindObject(fDigitsTitle)) ; 
  if (!phosd) { 
    if(strcmp(name, "Digitizer")==0){
      phosd = new AliPHOSDigitizer() ;
      phosd->SetName(fDigitsTitle) ;
      phosd->SetTitle(fHeaderFile) ;
      phos->Add(phosd) ;
    } 
    else{
      phosd = new AliPHOSRaw2Digits() ;
      phosd->SetName(fDigitsTitle) ;
      phosd->SetTitle(fHeaderFile) ;
      phos->Add(phosd) ;
    }      
  }
  return kTRUE;  
}

//____________________________________________________________________________ 
TObject** AliPHOSGetter::DigitizerRef(const char * name) const 
{  
  TTask * sd  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ; 
  if ( !sd ) {
    Fatal("DigitizerRef", "Task //%s/Digitizer not found !", fTasksFolder) ;
  }        

  TTask * phos = dynamic_cast<TTask*>(sd->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    Fatal("DigitizerRef", "//%s/Digitizer/PHOS", fTasksFolder) ;
  }        

  TTask * task = dynamic_cast<TTask*>(phos->GetListOfTasks()->FindObject(name)) ; 

  return phos->GetListOfTasks()->GetObjectRef(task) ;

}
 
//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostRecPoints(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/EMCARecPoints/name
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/CPVRecPoints/name

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS")) ; 
  
  if ( !phosFolder ) {
    if (fDebug) {
      Warning("PostRecPoints", "-> Folder //%s/PHOS/ not found!", fRecoFolder->GetName()) ;
      Info("PostRecPoints", "-> Adding Folder //%s/PHOS/", fRecoFolder->GetName()) ;
    }
    phosFolder = fRecoFolder->AddFolder("PHOS", "Reconstructed data from PHOS") ;  
  }    
  
  // EMCA RecPoints 
  TFolder * phosRPoEMCAFolder  = dynamic_cast<TFolder*>(phosFolder->FindObject("EMCARecPoints")) ;
  if ( !phosRPoEMCAFolder ) {
    if (fDebug) {
      Warning("PostRecPoints", "-> Folder //%s/PHOS/EMCARecPoints/ not found!", fRecoFolder->GetName()) ;
      Info("PostRecPoints", "-> Adding Folder //%s/PHOS/EMCARecPoints", fRecoFolder->GetName()) ;
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
      Warning("PostRecPoints", "-> Folder //%s/PHOS/CPVRecPoints/ not found!", fRecoFolder->GetName()) ;
      Info("PostRecPoints", "-> Adding Folder //%s/PHOS/CPVRecPoints/", fRecoFolder->GetName()) ;
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
TObject** AliPHOSGetter::EmcRecPointsRef(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/EMCARecPoints/name
   
  if ( !fRecoFolder ) {
    Fatal("EmcRecPointsRef", "Folder //%s not found !", fRecoFolder->GetName() ) ;
  }    

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/EMCARecPoints")) ; 
  if ( !phosFolder ) {
     Fatal("EmcRecPointsRef", "Folder //%s/PHOS/EMCARecPoints/ not found !", fRecoFolder->GetName() ) ;
  }    


  TObject * erp = phosFolder->FindObject(name ) ;
  if ( !erp )   {
    Fatal("EmcRecPointsRef", "object %s not found !", name) ; 
  }
  return phosFolder->GetListOfFolders()->GetObjectRef(erp) ;
  
} 

//____________________________________________________________________________ 
TObject** AliPHOSGetter::CpvRecPointsRef(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/CPVRecPoints/name
   
  if ( !fRecoFolder ) {
    Fatal("CpvRecPointsRef", "Folder //%s not found !", fRecoFolder->GetName() ) ;
  }    

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/CPVRecPoints")) ; 
  if ( !phosFolder ) {
    Fatal("CpvRecPointsRef", "Folder //%s/PHOS/CPVRecPoints/ not found !", fRecoFolder->GetName() ) ;
  }    

  TObject * crp = phosFolder->FindObject(name ) ;
  if ( !crp )   {
    Fatal("CpvRecPointsRef", "object %s nott found", name ) ; 
  }
  return phosFolder->GetListOfFolders()->GetObjectRef(crp) ;
  
} 

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostClusterizer(AliPHOSClusterizer * clu) const 
{ // ------------------ AliPHOSClusterizer ------------------------
  
  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    Error("PostClusterizer", "Task //%s/Reconstructioner not found !", fTasksFolder) ;
    return kFALSE ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Warning("PostClusterizer", "//%s/Reconstructioner/PHOS not found!", fTasksFolder) ; 
      Info("PostClusterizer", "Adding //%s/Reconstructioner/PHOS", fTasksFolder) ; 
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  AliPHOSClusterizer * phoscl = dynamic_cast<AliPHOSClusterizer*>(phos->GetListOfTasks()->FindObject(clu->GetName())) ; 
  if (phoscl) { 
    if (fDebug)
      Info("PostClusterizer", "Task %s already exists", clu->GetName()) ; 
    phoscl->Delete() ; 
    phos->GetListOfTasks()->Remove(phoscl) ;
  }
  phos->Add(clu) ;      
  return kTRUE; 
} 

//____________________________________________________________________________ 
TObject** AliPHOSGetter::ClusterizerRef(const char * name) const 
{ // ------------------ AliPHOSClusterizer ------------------------
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    Fatal("ClusterizerRef", "Task //%s/Reconstructioner not found !", fTasksFolder->GetName() ) ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    Fatal("ClusterizerRef", " //%s/Reconstructioner/PHOS not founf !", fTasksFolder->GetName() ) ; 
  }   

  TList * l = phos->GetListOfTasks() ; 
  TIter it(l) ;
  TTask * task ;
  TTask * clu = 0 ;
  TString cluname(name) ;
  cluname+=":clu" ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(cluname)){
      clu = task ;
      break ;
    }
  }

  if(!clu) {
    Fatal("ClusterizerRef", "Task //%s/Reconstructioner/clusterizer/%s not found", fTasksFolder->GetName(), name) ;
  }

  return l->GetObjectRef(clu) ;

}

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostClusterizer(const char * name) const 
{ // ------------------ AliPHOSClusterizer ------------------------

  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    Error("PostClusterizer", "Task//%s/Reconstructioner not found !", fTasksFolder) ; 
    return kFALSE ;
  }        
  
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Warning("PostClusterizer", "//%s/Reconstructioner/PHOS not found!", fTasksFolder) ;
      Info("PostClusterizer", "Adding //%s/Reconstructioner/PHOS", fTasksFolder) ;
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  TList * l = phos->GetListOfTasks() ;   
  TIter it(l) ;
  TString clun(name) ;
  clun+=":clu" ; 
  TTask * task ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(clun))
      return kTRUE ;
  }

  AliPHOSClusterizerv1 * phoscl = new AliPHOSClusterizerv1() ;
  clun+="-v1" ; 
  phoscl->SetName(clun) ;
  phoscl->SetTitle(fHeaderFile) ;
  phos->Add(phoscl) ;
  return kTRUE; 
  
}

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostTrackSegments(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/TrackSegments/name

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS")) ; 
  
  if ( !phosFolder ) {
    if (fDebug) {
      Warning("PostTrackSegments", "-> Folder //%s/PHOS/ not found", fRecoFolder->GetName()) ;
      Info("PostTrackSegments", "-> Adding Folder //%s/PHOS", fRecoFolder->GetName()) ;
    }
    phosFolder = fRecoFolder->AddFolder("PHOS", "Reconstructed data from PHOS") ;  
  }    

  TFolder * phosTSFolder  = dynamic_cast<TFolder*>(phosFolder->FindObject("TrackSegments")) ;
  if ( !phosTSFolder ) {
    if (fDebug) {
      Warning("PostTrackSegments", "-> Folder //%s/PHOS/TrackSegments/ not found!", fRecoFolder->GetName() ) ; 
      Info("PostTrackSegments", "-> Adding Folder //%s/PHOS/TrackSegments/", fRecoFolder->GetName()) ; 
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
TObject** AliPHOSGetter::TrackSegmentsRef(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/TrackSegments/name

 if ( !fRecoFolder ) {
    Fatal("TrackSegmentsRef", "Folder //%s not found !", fRecoFolder->GetName() ) ;
  }    

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/TrackSegments")) ; 
  if ( !phosFolder ) {
    Fatal("TrackSegmentsRef", "Folder //%s/PHOS/TrackSegments/ not found !", fRecoFolder->GetName() ) ;
  }    
  
  TObject * tss =  phosFolder->FindObject(name) ;
  if (!tss) {
    Fatal("TrackSegmentsRef", "object %s not found !", name) ;  
  }
  return phosFolder->GetListOfFolders()->GetObjectRef(tss) ;
} 

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostTrackSegmentMaker(AliPHOSTrackSegmentMaker * tsmaker) const 
{ //------------Track Segment Maker ------------------------------
  
  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    Error("PostTrackSegmentMaker", "Task //%s/Reconstructioner not found !", fTasksFolder) ;
    return kFALSE ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Warning("PostTrackSegmentMaker", "//%s/Reconstructioner/PHOS not found!", fTasksFolder) ; 
      Info("PostTrackSegmentMaker", "Adding //%s/Reconstructioner/PHOS", fTasksFolder) ;
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
const Bool_t AliPHOSGetter::PostTrackSegmentMaker(const char * name) const 
{ //------------Track Segment Maker ------------------------------
  
  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 
  
  if ( !tasks ) {
    Error("PostTrackSegmentMaker", "Task //%s/Reconstructioner not found !", fTasksFolder->GetName() ) ;
    return kFALSE ;
  }        
  
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Warning("PostTrackSegmentMaker", "//%s/Reconstructioner/PHOS not found!", fTasksFolder->GetName() ) ; 
      Info("PostTrackSegmentMaker", "Adding //%s/Reconstructioner/PHOS", fTasksFolder->GetName()) ;
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  TList * l = phos->GetListOfTasks() ;   
  TIter it(l) ;
  TString tsn(name);
  tsn+=":tsm" ; 
  TTask * task ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(tsn))
      return kTRUE ;
  }
  
  AliPHOSTrackSegmentMakerv1 * phosts = new AliPHOSTrackSegmentMakerv1() ;
  tsn+="-v1" ;
  phosts->SetName(tsn) ;
  phosts->SetTitle(fHeaderFile) ;
  phos->Add(phosts) ;      
  return kTRUE; 
  
} 

//____________________________________________________________________________ 
TObject** AliPHOSGetter::TSMakerRef(const char * name) const 
{ //------------Track Segment Maker ------------------------------
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    Fatal("TSMakerRef", "Task //%s/Reconstructioner not found !", fTasksFolder->GetName() ) ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    Fatal("TSMakerRef", "//%s/Reconstructioner/PHOS not found !", fTasksFolder->GetName() ) ; 
  }   

  TList * l = phos->GetListOfTasks() ; 
  TIter it(l) ;
  TTask * task ;
  TTask * tsm = 0 ;
  TString tsmname(name) ;
  tsmname+=":tsm" ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(tsmname)){
      tsm = task ;
      break ;
    }
  }
  
  if(!tsm) {
   Fatal("TSMakerRef", "Task //%s/Reconstructioner/PHOS/TrackSegmentMarker/%s not found !", fTasksFolder->GetName(),  name) ;
  }
 
  return l->GetObjectRef(tsm) ;

} 

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostRecParticles(const char * name) const 
{  // -------------------- RecParticles ------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/RecParticles/name

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS")) ; 
  
  if ( !phosFolder ) {
    if (fDebug) {
      Warning("PostRecParticles", "-> Folder //%s/PHOS/ not found!", fRecoFolder->GetName()) ;
      Info("PostRecParticles", "-> Adding Folder //%s/PHOS/", fRecoFolder->GetName()) ;
    }
    phosFolder = fRecoFolder->AddFolder("PHOS", "Reconstructed data from PHOS") ;  
  }    

 TFolder * phosRPaFolder  = dynamic_cast<TFolder*>(phosFolder->FindObject("RecParticles")) ;
  if ( !phosRPaFolder ) {
    if (fDebug) {
      Warning("PostRecParticles", "-> Folder //%s/PHOS/RecParticles/ not found!", fRecoFolder->GetName()) ;
      Info("PostRecParticles", "-> Adding Folder //%s/PHOS/RecParticles/", fRecoFolder->GetName()) ;
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
TObject** AliPHOSGetter::RecParticlesRef(const char * name) const 
{ // ---------------RecParticles -----------------------------------
  
  // the hierarchy is //Folders/Run/Event/RecData/PHOS/TrackSegments/name

 if ( !fRecoFolder ) {
    Fatal("RecParticlesRef", "Folder//%s not found !", fRecoFolder->GetName() ) ; 
  }    

  TFolder * phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/RecParticles")) ; 
  if ( !phosFolder ) {
    Fatal("RecParticlesRef", "Folder //%s/PHOS/RecParticles/ not found !", fRecoFolder->GetName() ) ;
  }    

  TObject * tss =  phosFolder->FindObject(name  ) ;
  if (!tss) {
    Fatal("RecParticlesRef", "object %s not found !", name) ; 
  }
  return phosFolder->GetListOfFolders()->GetObjectRef(tss) ;
}
//____________________________________________________________________________ 
const UShort_t AliPHOSGetter::EventPattern(void){
  if(fBTE)
    return fBTE->GetPattern() ;
  else
    return 0 ;
}
//____________________________________________________________________________ 
const Float_t AliPHOSGetter::BeamEnergy(void){
  if(fBTE)
    return fBTE->GetBeamEnergy() ;
  else
    return 0 ;
}
//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostPID(AliPHOSPID * pid) const 
{      // ------------AliPHOS PID -----------------------------

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    Error("PostPID", "Task //%s/Reconstructioner not found !", fTasksFolder) ;
    return kFALSE ;
  }        
  
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Warning("PostPID", "//%s/Reconstructioner/PHOS not found!", fTasksFolder) ; 
      Info("PostPID", "Adding //%s/Reconstructioner/PHOS", fTasksFolder) ;
    }
    phos = new TTask("PHOS", "") ; 
    tasks->Add(phos) ; 
  } 

  AliPHOSPID * phospid = dynamic_cast<AliPHOSPID*>(phos->GetListOfTasks()->FindObject(pid->GetName())) ; 
  if (phospid) { 
    if (fDebug)
      Info("PostPID", "-> Task %s qlready exists", pid->GetName()) ; 
    phos->GetListOfTasks()->Remove(phospid) ;
  }
  
  phos->Add(pid) ;      
  return kTRUE; 
} 

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostPID(const char * name) const 
{     
  // the hierarchy is //Folders/Tasks/Reconstructioner/PHOS/sdigitsname
  
  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    Error("PostPID", "Task //%s/Reconstructioner not found !", fTasksFolder->GetName() ) ;
    return kFALSE ;
  }        
  
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    if (fDebug) {
      Warning("PostPID", "//%s/Reconstructioner/PHOS not found!", fTasksFolder->GetName()) ; 
      Info("PostPID", "Adding //%s/Reconstructioner/PHOS", fTasksFolder->GetName()) ;
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
  phospid->SetTitle(fHeaderFile) ;
  phos->Add(phospid) ;      
  
  return kTRUE; 
} 

//____________________________________________________________________________ 
TObject** AliPHOSGetter::PIDRef(const char * name) const 
{ //------------PID ------------------------------

  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  if ( !tasks ) {
    Fatal("PIDRef", "Task //%s/Reconstructioner not found !", fTasksFolder->GetName() ) ;
  }        
        
  TTask * phos = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if ( !phos )  {
    Fatal("PIDRef", "//%s/Reconstructioner/PHOS not found !", fTasksFolder->GetName() ) ; 
  }   
  
  TList * l = phos->GetListOfTasks() ; 
  TIter it(l) ;
  TTask * task ;
  TTask * pid = 0 ;
  TString pidname(name) ;
  pidname+=":pid" ;
  while((task = static_cast<TTask *>(it.Next()) )){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(pidname)){
      pid = task ;
      break ;
    }
  }
  
  if(!pid) {
    Fatal("PIDRef", "Task //%s/Reconstructioner/PHOS/PID/%s not found !", fTasksFolder->GetName(), name) ;
  }
  
    return l->GetObjectRef(pid) ;
} 

//____________________________________________________________________________ 
const Bool_t AliPHOSGetter::PostQA(void) const 
{ // ------------------ QA ---------------------------------

  // the hierarchy is //Folders/Run/Conditions/QA/PHOS/alarmsName

  TFolder * phosFolder = dynamic_cast<TFolder*>(fQAFolder->FindObject("PHOS")) ; 
  if ( !phosFolder ) {
    if (fDebug) {
      Warning("PostQA", "-> Folder //%s/PHOS/ not found!", fQAFolder) ;
      Info("PostQA", "-> Adding Folder //%s/PHOS", fQAFolder) ;
    }
    phosFolder = fQAFolder->AddFolder("PHOS", "QA from PHOS") ; 
  }      

  return kTRUE;
}

//____________________________________________________________________________ 
TObject** AliPHOSGetter::AlarmsRef(void) const 
{  //------- Alarms ----------------------

  
  // the hierarchy is //Folders/Run/Conditions/QA/PHOS
  if ( !fQAFolder ) {
    Fatal("AlarmsRef", "Folder //%s not found !", fQAFolder) ;
  }    
 
  TFolder * phosFolder = dynamic_cast<TFolder *>(fQAFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) {
    Fatal("AlarmsRef", "Folder //%s/PHOS/ not found !", fQAFolder) ;
  }
   
  return fQAFolder->GetListOfFolders()->GetObjectRef(phosFolder) ;
}


//____________________________________________________________________________ 
TTree * AliPHOSGetter::TreeK(TString filename)  
{

  // returns TreeK from file filename
  // usefull in case of split file

  if ( filename.IsNull() ) 
    filename = fHeaderFile ; 

  TFile * file = 0 ; 
  file = static_cast<TFile*>(gROOT->GetFile(filename.Data() ) ) ;
//   if (file && (filename != fHeaderFile) ) {  // file already open 
//     file->Close() ; 
//     //delete fAlice ; 
//   }
  if(!file || !file->IsOpen())    
    file = TFile::Open(filename.Data(), "read") ;
  if(filename != fHeaderFile ){
    fAlice = dynamic_cast<AliRun *>(file->Get("gAlice")) ;
  } 
  TString treeName("TreeK") ; 
  treeName += EventNumber()  ; 
  TTree * tree = dynamic_cast<TTree *>(file->Get(treeName.Data())) ;
  if (!tree && fDebug)  
    Warning("TreeK", "-> %s not found in %s", treeName.Data(), filename.Data()) ; 
  
  return tree ; 		      
}

//____________________________________________________________________________ 
TTree * AliPHOSGetter::TreeH(TString filename)  
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
    Warning("TreeH", "-> %s not found in %s", treeName.Data(), filename.Data()) ; 
  
  return tree ; 		      
}

//____________________________________________________________________________ 
TTree * AliPHOSGetter::TreeS(TString filename)  
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
    Warning("TreeS", "-> %s not found in %s", treeName.Data(), filename.Data() ); 
  
  return tree ; 		      
}

//____________________________________________________________________________ 
TTree * AliPHOSGetter::TreeD(TString filename)  
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
    Warning("TreeD", "-> %s not found in %s", treeName.Data(), filename.Data()) ; 
  
  return tree ; 		      
}

//____________________________________________________________________________ 
const TParticle * AliPHOSGetter::Primary(Int_t index) const 
{
  // Return primary particle numbered by <index>

  if(index < 0) 
    return 0 ;
  TParticle *  p = 0 ;
  if (fAlice) 
    p = fAlice->Particle(index) ; 
  else 
    p = gAlice->Particle(index) ; 
  
  return p ; 
    
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
Int_t AliPHOSGetter::ReadTreeD(const Int_t event)
{
  // Read the digit tree gAlice->TreeD()  
  
  TTree * treeD ;
  if(fToSplit){
    TFile * file = static_cast<TFile*>(gROOT->GetFile(fDigitsFileName)); 
    if(!file) 
      file = TFile::Open(fDigitsFileName) ;      
    // Get Digits Tree header from file
    TString treeName("TreeD") ;
    treeName += event ; 
    treeD = dynamic_cast<TTree*>(file->Get(treeName.Data()));
    if(!treeD){ // TreeD not found in header file
      if (fDebug)
	Warning("ReadTreeD", "-> Cannot find TreeD in %s", fDigitsFileName.Data()) ;
      return 1;
    }
  }
  else
    treeD = gAlice->TreeD() ;
  
  TObjArray * lob = static_cast<TObjArray*>(treeD->GetListOfBranches()) ;
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
    else if ( ((strcmp(branch->GetName(), "AliPHOSDigitizer")==0)||
	       (strcmp(branch->GetName(), "AliPHOSRaw2Digits")==0)) &&
	      (strcmp(branch->GetTitle(), fDigitsTitle)==0) ) {
      digitizerbranch = branch ; 
      digitizerfound = kTRUE ; 
    }
  }
  
  if ( !phosfound || !digitizerfound ) {
    if (fDebug)
      Warning("ReadTreeD", "-> Cannot find Digits and/or Digitizer with name %s", fDigitsTitle.Data()) ;
    return 2; 
  }   
  
  //read digits
  if(!Digits(fDigitsTitle) ) 
    PostDigits(fDigitsTitle);
  else
    Digits(fDigitsTitle)->Clear() ;
  digitsbranch->SetAddress(DigitsRef(fDigitsTitle)) ;
  digitsbranch->GetEntry(0) ;
  
  // read  the Digitizer
  if(Digitizer()){
    if(strcmp(Digitizer()->IsA()->GetName(),digitizerbranch->GetName())!=0){
      RemoveTask("D", fDigitsTitle) ;
      if(strcmp(digitizerbranch->GetName(), "AliPHOSDigitizer")==0)
	PostDigitizer("Digitizer") ;
      else
	PostDigitizer("Raw2Digits") ;
    }
  }
  else{
    if(strcmp(digitizerbranch->GetName(), "AliPHOSDigitizer")==0)
      PostDigitizer("Digitizer") ;
    else
      PostDigitizer("Raw2Digits") ;
  }
    

  digitizerbranch->SetAddress(DigitizerRef(fDigitsTitle)) ;
  digitizerbranch->GetEntry(0) ;
  

//   if((!fcdb)&&(strcmp(digitizerbranch->GetName(), "AliPHOSRaw2Digits")==0))
//     ReadCalibrationDB("Primordial","beamtest.root") ;


  if(gAlice->TreeD()!=treeD)
    treeD->Delete();

  return 0 ; 
}
//____________________________________________________________________________ 
//void AliPHOSGetter::ReadCalibrationDB(const char * database,const char * filename){
//
//  if(fcdb && (strcmp(database,fcdb->GetTitle())==0))
//    return ;
//
//  TFile * file = gROOT->GetFile(filename) ;
//  if(!file)
//    file = TFile::Open(filename);
//  if(!file){
//    Error ("ReadCalibrationDB", "Cannot open file %s", filename) ;
//    return ;
//  }
//  if(fcdb)
//    fcdb->Delete() ;
//  fcdb = dynamic_cast<AliPHOSCalibrationDB *>(file->Get("AliPHOSCalibrationDB")) ;
//  if(!fcdb)
//    Error ("ReadCalibrationDB", "No database %s in file %s", database, filename) ;
//}

//____________________________________________________________________________ 
Int_t AliPHOSGetter::ReadTreeH()
{
  // Read the first entry of PHOS branch in hit tree gAlice->TreeH()
  
  TTree * treeH = gAlice->TreeH() ;

  if(!treeH) {// TreeH not found in header file
 
    if (fDebug) 
      Warning("ReadTreeH", "-> Cannot find TreeH in %s", fHeaderFile.Data() ) ;
    
    TString searchFileName("PHOS.Hits") ; 
    if((strcmp(fBranchTitle.Data(),"Default")!=0)&&(strcmp(fBranchTitle.Data(),"")!=0)){
      searchFileName+="." ;
      searchFileName += fBranchTitle ;
    }
    searchFileName+=".root" ;
    
    if ( (treeH = TreeH(searchFileName)) ) { //found TreeH in the file which contains the hits
      if (fDebug) 
	Info("ReadTreeH", "-> TreeH found in %s", searchFileName.Data()) ; 
      
    } else {
      Error("ReadTreeH", "TreeH not found") ; 
      return 1;
    }  
  }
  
  TBranch * hitsbranch = static_cast<TBranch*>(treeH->GetBranch("PHOS")) ;
  if ( !hitsbranch ) {
    if (fDebug)
      Warning("ReadTreeH", "-> Cannot find branch PHOS") ; 
    return 2;
  }
  if(!Hits())
    PostHits() ;

  if (hitsbranch->GetEntries() > 1 ) {
    (dynamic_cast<TClonesArray*> (*HitsRef()))->Clear() ;
    TClonesArray * tempo =  new TClonesArray("AliPHOSHit",1000) ;
    TClonesArray * hits = dynamic_cast<TClonesArray*>(*HitsRef()) ; 
    hitsbranch->SetAddress(&tempo) ;
    Int_t index = 0 ; 
    Int_t i = 0 ;
    for (i = 0 ; i < hitsbranch->GetEntries() ; i++) {
      hitsbranch->GetEntry(i) ;
      Int_t j = 0 ; 
      for ( j = 0 ; j < tempo->GetEntries() ; j++) { 
	const AliPHOSHit * hit = static_cast<const AliPHOSHit *>(tempo->At(j)) ; 
	new((*hits)[index]) AliPHOSHit( *hit ) ;
	index++ ; 
      }
    }
    delete tempo ; 
  }
  else {
    (dynamic_cast<TClonesArray*> (*HitsRef()))->Clear() ;
    hitsbranch->SetAddress(HitsRef()) ;
    hitsbranch->GetEntry(0) ;
  }
  return 0 ; 
}

//____________________________________________________________________________ 
void AliPHOSGetter::Track(const Int_t itrack) 
{
  // Read the first entry of PHOS branch in hit tree gAlice->TreeH()

  if(gAlice->TreeH()== 0){
    Error("Track", "Cannot read TreeH") ;
    return ;
  }
  
  TBranch * hitsbranch = dynamic_cast<TBranch*>(gAlice->TreeH()->GetListOfBranches()->FindObject("PHOS")) ;
  if ( !hitsbranch ) {
    if (fDebug)
      Warning("Track", "Cannot find branch PHOS") ; 
    return ;
  }  
  if(!Hits())
    PostHits() ;

  (dynamic_cast<TClonesArray*> (*HitsRef()))->Clear() ;
  hitsbranch->SetAddress(HitsRef()) ;
  hitsbranch->GetEntry(itrack) ;

}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeQA()
{
  // Read the digit tree gAlice->TreeQA()
  // so far only PHOS knows about this Tree  

  if(PHOS()->TreeQA()== 0){
    Error("ReadTreeQA", "Cannot read TreeQA") ;
    return ;
  }
  
  TBranch * qabranch = PHOS()->TreeQA()->GetBranch("PHOS") ; 
  if (!qabranch) { 
    if (fDebug)
      Warning("ReadTreeQA", "Cannot find QA Alarms for PHOS");
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
Int_t AliPHOSGetter::ReadTreeR(const Int_t event)
{
  // Read the reconstrunction tree gAlice->TreeR()
  // A particularity has been introduced here :
  //  if gime->Event(ievent,"R") is called branches with the current title are read, the current title
  //   being for example give in AliPHOSPID(fileName, title)
  //  if gime(Event(ievent, "RA") is called the title of the branches is not checked anymore, "A" stands for any
  // This is a feature needed by PID to be able to reconstruct several times particles (each time a ther title is given)
  // from a given set of TrackSegments (with a given name)
  // This is why any is NOT used to read the branch of RecParticles
  // any migh have become obsolete : to be checked
  // See AliPHOSPIDv1    

  //first - clean if necessary
  if(EmcRecPoints(fRecPointsTitle)){
    EmcRecPoints(fRecPointsTitle)->Delete() ;
    CpvRecPoints(fRecPointsTitle)->Delete() ;
  }
  //clear TrackSegments
  if(TrackSegments(fTrackSegmentsTitle))
    TrackSegments(fTrackSegmentsTitle)->Clear() ;
  if(RecParticles(fRecParticlesTitle))
    RecParticles(fRecParticlesTitle)->Clear() ;
	
  TTree * treeR ;
  if(fToSplit){
    TFile * file = static_cast<TFile*>(gROOT->GetFile(fRecPointsFileName)); 
    if(!file) 
      file = TFile::Open(fRecPointsFileName) ;      
    // Get Digits Tree header from file
    TString treeName("TreeR") ;
    treeName += event ; 
    treeR = dynamic_cast<TTree*>(file->Get(treeName.Data()));
    if(!treeR){ // TreeR not found in header file
      if (fDebug)
	Warning("ReadTreeD", "-> Cannot find TreeR in %s", fRecPointsFileName.Data()) ;
      return 1;
    }
  }
  else
    treeR = gAlice->TreeR() ;
  
  // RecPoints 
  TObjArray * lob = static_cast<TObjArray*>(treeR->GetListOfBranches()) ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * emcbranch = 0 ; 
  TBranch * cpvbranch = 0 ; 
  TBranch * clusterizerbranch = 0 ; 
  Bool_t phosemcrpfound = kFALSE, phoscpvrpfound = kFALSE, clusterizerfound = kFALSE ; 

  
  while ( (branch = static_cast<TBranch*>(next())) && (!phosemcrpfound || !phoscpvrpfound || !clusterizerfound) ) {
    if(strcmp(branch->GetTitle(), fRecPointsTitle)==0 ) {
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
  }

  if ( !phosemcrpfound || !phoscpvrpfound || !clusterizerfound) {
    if (fDebug)
      Warning("ReadTreeR", "-> Cannot find RecPoints and/or Clusterizer with name %s", fRecPointsTitle.Data() ) ;
 
  } else { 
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
  }
  
  //------------------- TrackSegments ---------------------
  next.Reset() ; 
  TBranch * tsbranch = 0 ; 
  TBranch * tsmakerbranch = 0 ; 
  Bool_t phostsfound = kFALSE, tsmakerfound = kFALSE ; 
  while ( (branch = static_cast<TBranch*>(next())) && (!phostsfound || !tsmakerfound) ) {
    if(strcmp(branch->GetTitle(), fTrackSegmentsTitle)==0 )  {
      if ( strcmp(branch->GetName(), "PHOSTS")==0){
	tsbranch = branch ; 
	phostsfound = kTRUE ;
      }
      else if(strcmp(branch->GetName(), "AliPHOSTrackSegmentMaker")==0) {
	tsmakerbranch = branch ; 
	tsmakerfound  = kTRUE ; 
      }
    }
  }

  if ( !phostsfound || !tsmakerfound ) {
    if (fDebug)
      Warning("ReadTreeR", "-> Cannot find TrackSegments and/or TrackSegmentMaker with name %s", fTrackSegmentsTitle.Data() ) ;
  } else { 
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
 }
  
  
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
      Warning("ReadTreeR", "-> Cannot find RecParticles and/or PID with name %s", fRecParticlesTitle.Data() ) ; 
  } else { 
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

  if(gAlice->TreeR()!=treeR)
    treeR->Delete();
  return 0 ; 
}

//____________________________________________________________________________ 
Int_t AliPHOSGetter::ReadTreeS(const Int_t event)
{
  // Reads the SDigits treeS from all files  
  // Files, which should be opened are listed in phosF
  // So, first get list of files
  TFolder * phosF = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject("PHOS")) ;
  if (!phosF) 
    phosF = fSDigitsFolder->AddFolder("PHOS", "SDigits from PHOS") ; 
  TCollection * folderslist = phosF->GetListOfFolders() ; 
  
  // Now iterate over the list of files and read TreeS into Whiteboard
  TIter next(folderslist) ; 
  TFolder * folder = 0 ; 
  TFile * file; 
  TTree * treeS = 0;
  while ( (folder = static_cast<TFolder*>(next())) ) {
    TString fileName("") ;
    fileName = folder->GetName() ; 
    fileName.ReplaceAll("_","/") ; 
    file = static_cast<TFile*>(gROOT->GetFile(fileName)); 
    if(!file) 
      file = TFile::Open(fileName) ;      
    // Get SDigits Tree header from file
    TString treeName("TreeS") ;
    treeName += event ; 
    treeS = dynamic_cast<TTree*>(file->Get(treeName.Data()));

    if(!treeS){ // TreeS not found in header file
      if (fDebug)
	Warning("ReadTreeS", "-> Cannot find TreeS in %s", fileName.Data()) ;
      return 1;
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
      
      else if ( (strcmp(branch->GetName(), "AliPHOSSDigitizer")==0) && 
		(strcmp(branch->GetTitle(), fSDigitsTitle)==0) ) {
	sdigitizerfound = kTRUE ; 
	sdigitizerBranch = branch ;
      }
    }
    if ( !phosfound || !sdigitizerfound ) {
      if (fDebug)
	Warning("ReadSDigits", "-> Digits and/or Digitizer branch with name %s not found", GetName()) ;
      return 2; 
    }   
    
    if ( !folder->FindObject(fSDigitsTitle) ){  
      TClonesArray * sdigits = new TClonesArray("AliPHOSDigit",1) ;
      sdigits->SetName(fSDigitsTitle) ;
      folder->Add(sdigits) ;
    }

    ((TClonesArray*) (*SDigitsRef(fSDigitsTitle,folder->GetName())))->Clear() ;
    sdigitsBranch->SetAddress(SDigitsRef(fSDigitsTitle,folder->GetName())) ;
    sdigitsBranch->GetEntry(0) ;
    
    TString sdname(fSDigitsTitle) ;
    sdname+=":" ;
    sdname+=folder->GetName() ;
    if(!SDigitizer(sdname) ) 
      PostSDigitizer(fSDigitsTitle,folder->GetName()) ;
    sdigitizerBranch->SetAddress(SDigitizerRef(sdname)) ;
    sdigitizerBranch->GetEntry(0) ; 
    if(gAlice->TreeS()!=treeS)
      treeS->Delete();
  }    
  return 0 ; 
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
      Warning("ReadTreeS", "-> Digits and/or Digitizer branch not found") ;
    return ; 
  }   
  
  if (!folder || !(folder->FindObject(sdigitsBranch->GetTitle()) ) )
    PostSDigits(sdigitsBranch->GetTitle(),filename) ;

  sdigitsBranch->SetAddress(SDigitsRef(sdigitsBranch->GetTitle(),folder->GetName())) ;
  sdigitsBranch->GetEntry(0) ;
  
  TString sdname(sdigitsBranch->GetTitle()) ;
  sdname+=":" ;
  sdname+=filename ;
  
  if(!SDigitizer(sdigitsBranch->GetTitle()) )
    PostSDigitizer(sdigitsBranch->GetTitle(),filename) ;
  sdigitizerBranch->SetAddress(SDigitizerRef(sdname)) ;
  sdigitizerBranch->GetEntry(0) ;
  if(gAlice->TreeS()!=treeS)
    treeS->Delete();
}    


//____________________________________________________________________________ 
void AliPHOSGetter::ReadPrimaries()
{
  // a lot simplified.... if 2 files are opened then we have a problem

  TClonesArray * ar = 0  ; 
  if(! (ar = Primaries()) ) { 
    PostPrimaries() ;
    ar = Primaries() ; 
  }
  ar->Delete() ; 
  
  if (TreeK(fHeaderFile)) { // treeK found in header file
    if (fDebug) 
      Info("ReadPrimaries", "-> TreeK found in %s", fHeaderFile.Data() ); 
    fNPrimaries = gAlice->GetNtrack() ; 
    fAlice = 0 ; 
  
  } else { // treeK not found in header file
    Error("ReadPrimaries", "TreeK not found") ; 
    return ;  
  }

  Int_t index = 0 ; 
  for (index = 0 ; index < fNPrimaries; index++) { 
    new ((*ar)[index]) TParticle(*(Primary(index)));
  }
}

//____________________________________________________________________________ 
void AliPHOSGetter::Event(const Int_t event, const char* opt)  
{
  // Reads the content of all Tree's S, D and R

  if (event >= gAlice->TreeE()->GetEntries() ) {
    Error("Event", "%d not found in TreeE !", event) ; 
    return ; 
  }

  TBranch * btb = gAlice->TreeE()->GetBranch("AliPHOSBeamTestEvent") ;
  if(btb){
    if(!fBTE)
      fBTE = new AliPHOSBeamTestEvent() ;
    btb->SetAddress(&fBTE) ;
    btb->GetEntry(event) ;
  }
  else{
    if(fBTE){
      delete fBTE ;
      fBTE = 0 ;
    }
  }

  Bool_t any = kFALSE ; 
  if (strstr(opt,"A") ) // do not check the title of the branches
    any = kTRUE; 

  gAlice->GetEvent(event) ; 

  if( strstr(opt,"R") )
    ReadTreeR(event) ;

  if( strstr(opt,"D") )
    ReadTreeD(event) ;

  if(strstr(opt,"S") )
    ReadTreeS(event) ;

  if(strstr(opt,"H") )
    ReadTreeH() ;
   
  if( strstr(opt,"Q") )
    ReadTreeQA() ;

  if( strstr(opt,"P") || (strcmp(opt,"")==0) )
    ReadPrimaries() ;

 
}

//____________________________________________________________________________ 
TObject * AliPHOSGetter::ReturnO(TString what, TString name, TString file) const 
{
  // get the object named "what" from the folder
  // folders are named like //Folders

  if ( file.IsNull() ) 
    file = fHeaderFile ; 
  if( name.IsNull() )
    name = fBranchTitle ;

  TFolder * folder = 0 ;
  TObject * phosO  = 0 ; 

  if ( what.CompareTo("Primaries") == 0 ) {
    folder = dynamic_cast<TFolder *>(fPrimariesFolder->FindObject("Primaries")) ; 
    if (folder) 
      phosO  = dynamic_cast<TObject *>(folder->FindObject("Primaries")) ;  
    else 
      return 0 ; 
  }
  else if ( what.CompareTo("Hits") == 0 ) {
    folder = dynamic_cast<TFolder *>(fHitsFolder->FindObject("PHOS")) ; 
    if (folder) 
      phosO  = dynamic_cast<TObject *>(folder->FindObject("Hits")) ;  
  }
  else if ( what.CompareTo("SDigits") == 0 ) {
    if(fToSplit){
      file.Remove(file.Last('/')+1,file.Length()-file.Last('/')-1) ;
      file.ReplaceAll("/","_") ; 
      file+="PHOS.SDigits." ;
      if(name && (strcmp(name,"Default")!=0)){
	file+=name ;
	file+="." ;
      }
      file+="root" ;
    }
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
      Warning("ReturnO", "Object %s not found in PHOS", what.Data() ) ; 
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
    Error("ReturnT", "Task %s not found !", what.Data() ) ;  
    return 0 ; 
  }

  TTask * phosT = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("PHOS")) ; 
  if (!phosT) { 
     Error("ReturnT", "Task %s/PHOS not found !", what.Data() ) ;  
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
    Warning("ReturnT", "-> Task %s/PHOS/%s not found", search.Data(), name.Data() ) ; 
  return 0 ;
}

//____________________________________________________________________________ 
void AliPHOSGetter::RemoveTask(TString opt, TString name) const 
{
  // remove a task from the folder
  // path is fTasksFolder/SDigitizer/PHOS/name
  
  TTask * task = 0 ; 
  TTask * phos = 0 ; 
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
  else if (opt == "C" || opt == "T" || opt == "P"  ) { // Clusterizer, TrackSegmentMaker, PID
    task = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ;
    if (!task) 
      return ; 
  }
  else {
    Warning("RemoveTask", "Unknown option %s", opt.Data() ); 
    return ; 
  }
  phos =  dynamic_cast<TTask*>(task->GetListOfTasks()->FindObject("PHOS")) ;
  if (!phos)
    return ; 
  lofTasks = phos->GetListOfTasks() ;
  if (!lofTasks) 
    return ; 
  TObject * obj = lofTasks->FindObject(name) ; 
  if (obj) 
    lofTasks->Remove(obj) ;
   
}

//____________________________________________________________________________ 
void AliPHOSGetter::RemoveObjects(TString opt, TString name) const 
{
  // remove SDigits from the folder
  // path is fSDigitsFolder/fHeaderFileName/name

  TFolder * phos     = 0 ; 
  TFolder * phosmain = 0 ; 

  if (opt == "H") { // Hits
    phos = dynamic_cast<TFolder*>(fHitsFolder->FindObject("PHOS")) ;
    if (!phos) 
      return ;
    name = "Hits" ; 
  }

  else if ( opt == "S") { // SDigits
    phosmain = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("PHOS")) ;
    if (!phosmain) 
      return ;
    phos = dynamic_cast<TFolder*>(phosmain->FindObject(fHeaderFile)) ;
    if (!phos) 
      return ;
  }
  
  else if (opt == "D") { // Digits
    phos = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("PHOS")) ;
    if (!phos) 
      return ;
  }

  else if (opt == "RE") { // EMCARecPoints
    phos = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/EMCARecPoints")) ;
    if (!phos) 
      return ;
  }

  else if (opt == "RC") { // CPVRecPoints
    phos = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/CPVRecPoints")) ;
    if (!phos) 
      return ;
  }  

  else if (opt == "T") { // TrackSegments
    phos = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/TrackSegments")) ;
    if (!phos) 
      return ;
  }

  else if (opt == "P") { // RecParticles
    phos = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/RecParticles")) ;
    if (!phos) 
      return ;
  }
  
  else {
    Warning("RemoveObjects", "Unknown option %s", opt.Data() ) ; 
    return ; 
  }
  
  TObjArray * ar  = dynamic_cast<TObjArray*>(phos->FindObject(name)) ; 
  if (ar) { 
    phos->Remove(ar) ;
    ar->Delete() ; 
    delete ar ; 
  }

  if (opt == "S") 
    phosmain->Remove(phos) ; 
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::RemoveSDigits() const 
{
  TFolder * phos= dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("PHOS")) ;
  if (!phos) 
    return ;
  
  phos->SetOwner() ; 
  phos->Clear() ; 
}

//____________________________________________________________________________ 
void AliPHOSGetter::CleanWhiteBoard(void){

  TFolder * phosmain = 0 ; 
  TFolder * phos ;
  TObjArray * ar ;
  TList * lofTasks = 0 ; 
  TTask * task = 0 ; 
  TTask * phost = 0 ; 
  
  // Hits  
  phos = dynamic_cast<TFolder*>(fHitsFolder->FindObject("PHOS")) ;
  if (phos){  
    TObjArray * ar  = dynamic_cast<TObjArray*>(phos->FindObject("Hits")) ; 
    if (ar) { 
      phos->Remove(ar) ;
      ar->Delete() ; 
      delete ar ; 
    }
  }
  
  // SDigits
  phosmain = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("PHOS")) ;
  if (phosmain){ 
    phos = dynamic_cast<TFolder*>(phosmain->FindObject(fHeaderFile)) ;
    if (phos) {
      ar  = dynamic_cast<TObjArray*>(phos->FindObject(fSDigitsTitle)) ; 
      if (ar) { 
	phos->Remove(ar) ;
	ar->Delete() ; 
	delete ar ; 
      }
    }
    phosmain->Remove(phos) ; 
  }

  
  // Digits
  phos = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("PHOS")) ;
  if (phos){ 
    ar  = dynamic_cast<TObjArray*>(phos->FindObject(fDigitsTitle)) ; 
    if (ar) { 
      phos->Remove(ar) ;
      ar->Delete() ; 
      delete ar ; 
    }
  }


  // EMCARecPoints
  phos = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/EMCARecPoints")) ;
  if (phos){ 
    ar  = dynamic_cast<TObjArray*>(phos->FindObject(fRecPointsTitle)) ; 
    if (ar) { 
      phos->Remove(ar) ;
      ar->Delete() ; 
      delete ar ; 
    }
  }

  
  // CPVRecPoints
  phos = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/CPVRecPoints")) ;
  if (phos){ 
    ar  = dynamic_cast<TObjArray*>(phos->FindObject(fRecPointsTitle)) ; 
    if (ar) { 
      phos->Remove(ar) ;
      ar->Delete() ; 
      delete ar ; 
    }
  }  

  
  // TrackSegments
  phos = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/TrackSegments")) ;
  if (phos) { 
    ar  = dynamic_cast<TObjArray*>(phos->FindObject(fTrackSegmentsTitle)) ; 
    if (ar) { 
      phos->Remove(ar) ;
      ar->Delete() ; 
      delete ar ; 
    }
  }
  


  // RecParticles
  phos = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS/RecParticles")) ;
  if (phos){ 
    ar  = dynamic_cast<TObjArray*>(phos->FindObject(fRecParticlesTitle)) ; 
    if (ar) { 
      phos->Remove(ar) ;
      ar->Delete() ; 
      delete ar ; 
    }
  }


  //---- Now Tasks ----------- 

  TObject * obj ;
  TString sdname(fSDigitsTitle);
  
  // Digitizer
  task = dynamic_cast<TTask*>(fTasksFolder->FindObject("Digitizer")) ;
  if (task){ 
    phost =  dynamic_cast<TTask*>(task->GetListOfTasks()->FindObject("PHOS")) ;
    if (phost){
      lofTasks = phost->GetListOfTasks() ;
      if (lofTasks){ 
	obj = lofTasks->FindObject(sdname.Data()) ; 
	if (obj) 
	  lofTasks->Remove(obj) ;
      }
    }      
  }
  

  sdname.Append(":") ;
  // Clusterizer, TrackSegmentMaker, PID
  task = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ;
  if (task){ 
    phost =  dynamic_cast<TTask*>(task->GetListOfTasks()->FindObject("PHOS")) ;
    if (phost){
      lofTasks = phost->GetListOfTasks() ;
      TIter next(lofTasks);
      while((obj=next())){ 
	TString oname(obj->GetName()) ;
	if (oname.BeginsWith(sdname)){ 
	  lofTasks->Remove(obj) ;
	}
      }
    }  
  }


  // SDigitizer
  sdname.Append(fHeaderFile) ;
  task = dynamic_cast<TTask*>(fTasksFolder->FindObject("SDigitizer")) ;
  if (task) {
    phost =  dynamic_cast<TTask*>(task->GetListOfTasks()->FindObject("PHOS")) ;
    if (phost){
      lofTasks = phost->GetListOfTasks() ;
      if (lofTasks){ 
	obj = lofTasks->FindObject(sdname.Data()) ; 
	if (obj) 
	  lofTasks->Remove(obj) ;
      }
    }
  }  

}
//____________________________________________________________________________ 
void AliPHOSGetter::SetTitle(const char * branchTitle ) 
{
  fBranchTitle        = branchTitle ;
  fSDigitsTitle       = branchTitle ; 
  fDigitsTitle        = branchTitle ; 
  fRecPointsTitle     = branchTitle ; 
  fRecParticlesTitle  = branchTitle ; 
  fTrackSegmentsTitle = branchTitle ; 
  if(fToSplit){
    //First - extract full path if necessary
    TString sFileName(fHeaderFile) ;
    Ssiz_t islash = sFileName.Last('/') ;
    if(islash<sFileName.Length())
      sFileName.Remove(islash+1,sFileName.Length()) ;
    else
      sFileName="" ;
    //Now construct file names
    fSDigitsFileName       = sFileName ;
    fDigitsFileName        = sFileName ; 
    fRecPointsFileName     = sFileName ; 
    fRecParticlesFileName  = sFileName ; 
    fTrackSegmentsFileName = sFileName ; 
    fSDigitsFileName      += "PHOS.SDigits." ;
    fDigitsFileName       += "PHOS.Digits." ; 
    fRecPointsFileName    += "PHOS.RecData." ; 
    fTrackSegmentsFileName+= "PHOS.RecData." ; 
    fRecParticlesFileName += "PHOS.RecData." ; 
    if((strcmp(fBranchTitle.Data(),"Default")!=0)&&(strcmp(fBranchTitle.Data(),"")!=0)){
      fSDigitsFileName      += fBranchTitle ;
      fSDigitsFileName      += "." ;
      fDigitsFileName       += fBranchTitle ; 
      fDigitsFileName       += "." ; 
      fRecPointsFileName    += fBranchTitle ; 
      fRecPointsFileName    += "." ; 
      fRecParticlesFileName += fBranchTitle ; 
      fRecParticlesFileName += "." ; 
      fTrackSegmentsFileName+= fBranchTitle ; 
      fTrackSegmentsFileName+= "." ; 
    }
    fSDigitsFileName      += "root" ;
    fDigitsFileName       += "root" ; 
    fRecPointsFileName    += "root" ; 
    fRecParticlesFileName += "root" ; 
    fTrackSegmentsFileName+= "root" ; 
  }else{
    fSDigitsFileName       = fHeaderFile ;

    fDigitsFileName        = "" ; 
    fRecPointsFileName     = "" ; 
    fRecParticlesFileName  = "" ; 
    fTrackSegmentsFileName = "" ; 
  }
  TFolder * phosFolder ; 
  phosFolder = dynamic_cast<TFolder*>(fHitsFolder->FindObject("PHOS")) ; 
  if ( !phosFolder ) 
    phosFolder = fHitsFolder->AddFolder("PHOS", "Hits from PHOS") ; 

  phosFolder = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) 
    phosFolder = fSDigitsFolder->AddFolder("PHOS", "SDigits from PHOS") ; 

  //Make folder for SDigits
  fSDigitsFileName.ReplaceAll("/","_") ;
  phosFolder->AddFolder(fSDigitsFileName.Data(),"");

  phosFolder  = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("PHOS")) ;
  if ( !phosFolder ) 
    phosFolder = fDigitsFolder->AddFolder("PHOS", "Digits from PHOS") ;  

  phosFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("PHOS")) ; 
  if ( !phosFolder )
    phosFolder = fRecoFolder->AddFolder("PHOS", "Reconstructed data from PHOS") ;  
  
}
//____________________________________________________________________________ 
void AliPHOSGetter::CloseSplitFiles(void){
  TFile * file ;
  file = static_cast<TFile*>(gROOT->GetFile(fSDigitsFileName.Data() ) ) ;
  if(file)
    file->Close() ;
  file = static_cast<TFile*>(gROOT->GetFile(fDigitsFileName.Data() ) ) ;
  if(file)
    file->Close() ;
  file = static_cast<TFile*>(gROOT->GetFile(fRecPointsFileName.Data() ) ) ;
  if(file)
    file->Close() ;
  file = static_cast<TFile*>(gROOT->GetFile(fTrackSegmentsFileName.Data() ) ) ;
  if(file)
    file->Close() ;
  file = static_cast<TFile*>(gROOT->GetFile(fRecParticlesFileName.Data() ) ) ;
  if(file)
    file->Close() ;

}
