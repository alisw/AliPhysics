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

//*-- YS August 2002 : clone PHOS as closely as possible and intoduction

//                     of new  IO (à la PHOS)

//////////////////////////////////////////////////////////////////////////////



// --- ROOT system ---

#include "TFile.h"

#include "TTree.h"

#include "TROOT.h"

#include "TObjString.h"

#include "TFolder.h"

#include "TParticle.h"



// --- Standard library ---

#include <Riostream.h>



// --- AliRoot header files ---

#include "AliRun.h"

#include "AliConfig.h"

#include "AliEMCALGetter.h"

#include "AliEMCAL.h"

#include "AliEMCALDigitizer.h"

#include "AliEMCALSDigitizer.h"

#include "AliEMCALClusterizerv1.h"

//#include "AliEMCALTrackSegmentMakerv1.h"

//#include "AliEMCALTrackSegment.h"

//#include "AliEMCALPIDv1.h" 

#include "AliEMCALGeometry.h"



ClassImp(AliEMCALGetter)

  

  AliEMCALGetter * AliEMCALGetter::fgObjGetter = 0 ; 

  TFile * AliEMCALGetter::fFile = 0 ; 



//____________________________________________________________________________ 

AliEMCALGetter::AliEMCALGetter(const char* headerFile, const char* branchTitle, const Bool_t toSplit)

{

  // This is the ctor called by GetInstance and the only one that can be used 

  

  if ( fHeaderFile.Contains("_") ) {

    cerr << "AliEMCALGetter::AliEMCALGetter -> Invalid file name (_ not allowed) " << fHeaderFile.Data() << endl ;

    abort() ; 

  }

  

  //Initialize  all data



  fFailed = kFALSE ;   

  fDebug  = 0 ; 

  fToSplit    = toSplit ;

  fHeaderFile = headerFile ; 



  fPrimaries = new TObjArray(1) ;



  fModuleFolder    = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Configuration/Modules")); 

  fPrimariesFolder = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data")); 

  fHitsFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/Hits")); 

  fSDigitsFolder   = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/RunMC/Event/Data/SDigits")); 

  fDigitsFolder    = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Event/Data")); 

  fRecoFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Event/RecData")); 

  //fQAFolder      = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Run/Conditions/QA")); 

  fTasksFolder     = dynamic_cast<TFolder*>(gROOT->FindObjectAny("Folders/Tasks")) ; 



  //Set titles to branches and create PHOS specific folders

  SetTitle(branchTitle) ;

  

  if ( fHeaderFile != "aliroot"  ) { // to call the getter without a file

    //open headers file

    fFile = static_cast<TFile*>(gROOT->GetFile(fHeaderFile.Data() ) ) ;

  

    if(!fFile){    //if file was not opened yet, read gAlice

      fFile = TFile::Open(fHeaderFile.Data(), "update") ;   

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



  if (fFile) { 

    fFile->Close() ;  

    delete fFile ; 

    fFile = 0 ;

  }

  fgObjGetter = 0 ; 



}



//____________________________________________________________________________ 

void AliEMCALGetter::CloseFile()

{

  delete gAlice ;  

  gAlice = 0 ; 

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

    return fgObjGetter ;

  }

  else {

    // cout << "AliEMCALGetter::GetInstance ERROR: not yet initialized" << endl ;

    return 0 ;

  }

}



//____________________________________________________________________________ 

AliEMCALGetter * AliEMCALGetter::GetInstance(const char* headerFile,

					     const char* branchTitle,

					     const Bool_t toSplit)

{

  // Creates and returns the pointer of the unique instance

  // Must be called only when the environment has changed 



  if(!fgObjGetter){

    fgObjGetter = new AliEMCALGetter(headerFile,branchTitle,toSplit) ;

    if(fgObjGetter->fFailed)

      return 0;

    else

      return fgObjGetter ;

  }



  //First checks, if header file already opened

  if(!fgObjGetter->fFile){

     fgObjGetter = new AliEMCALGetter(headerFile,branchTitle,toSplit) ;

    if(fgObjGetter->fFailed)

      return 0;

    else

      return fgObjGetter ;

  }

 

  if(fgObjGetter->fHeaderFile.CompareTo(headerFile)==0){ //Opened the same header file   

    if((fgObjGetter->fBranchTitle.CompareTo(branchTitle) == 0)&&   //Open the same branch title

       (toSplit==fgObjGetter->fToSplit)){                          //Nothing should be cleaned

      return fgObjGetter ;

    }

    else{ //Clean all data and AliEMCAL...zers

      if(fgObjGetter->fToSplit)

	fgObjGetter->CloseSplitFiles() ;	  

      fgObjGetter->CleanWhiteBoard() ;

      fgObjGetter->fToSplit = toSplit ;

      fgObjGetter->SetTitle(branchTitle) ;

      return fgObjGetter ;

    }

  }

  else{  //Close already opened files, clean memory and open new header file

    if(gAlice)

      delete gAlice ;

    gAlice= 0;

    if(fgObjGetter->fFile){

      fgObjGetter->fFile->Close() ;

      fgObjGetter->fFile=0;

    }

    if(fgObjGetter->fToSplit)

      fgObjGetter->CloseSplitFiles() ;

    fgObjGetter->CleanWhiteBoard() ;    

    fgObjGetter = new AliEMCALGetter(headerFile,branchTitle,toSplit) ;

    return fgObjGetter ; 

  }

  return fgObjGetter ; 

  

}



//____________________________________________________________________________ 

const Bool_t AliEMCALGetter::BranchExists(const TString recName) const

{

  //Looks in the tree Tree"name" if branch with current name olready exists

  

  TString filename("") ;

  TString name, dataname, zername;

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

      if(recName =="RecPoints"){

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

      cerr << "WARNING: AliEMCALGetter::BranchExists -> branch " << dataname.Data() << " with title " << fBranchTitle << " already exits in "

	   << name.Data() << endl;

      return kTRUE ;

    }

    if ( branchName.BeginsWith(zername) &&  branchTitle.BeginsWith(titleName) ){

      cerr << "WARNING: AliEMCALGetter::BranchExists -> branch AliEMCAL... with title " << branch->GetTitle() << " already exits in "

	   << name.Data() << endl;     

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

void AliEMCALGetter::ListBranches(Int_t event) const  

{

  

  TBranch * branch = 0 ; 

  if (gAlice->GetEvent(event) == -1)

    return ; 



  TTree * t =  gAlice->TreeH() ; 

  if(t){

    cout << "INFO: AliEMCALGetter::ListBranches -> ****** Hits    : " << endl ; 

    TObjArray * lob = t->GetListOfBranches() ;

    TIter next(lob) ; 

    while ( (branch = static_cast<TBranch*>(next())) )

      cout << "             " << branch->GetName() << endl ; 

  } else 

    cerr << "WARNING::AliEMCALGetter::ListBranches -> TreeH not found for event " << event << endl ;  

  

  

  t = gAlice->TreeS() ;

  if(t){

    cout << "INFO: AliEMCALGetter::ListBranches -> ****** SDigits : " << endl ; 

    TObjArray * lob = t->GetListOfBranches() ;

    TIter next(lob) ; 

    while ( (branch = static_cast<TBranch*>(next())) )

      cout << "             " << branch->GetName() << " " << branch->GetTitle() << endl ; 

  } else 

    cerr << "WARNING::AliEMCALGetter::ListBranches -> TreeS not found for event " << event << endl ;  

    

  t = gAlice->TreeD() ;

  if(t){

    cout << "INFO: AliEMCALGetter::ListBranches -> ****** Digits  : " << endl ; 

    TObjArray * lob = t->GetListOfBranches() ;

    TIter next(lob) ; 

    while ( (branch = static_cast<TBranch*>(next())) )

      cout << "             " << branch->GetName() << " " << branch->GetTitle() << endl ; 

  } else 

    cerr << "WARNING::AliEMCALGetter::ListBranches -> TreeD not found for event " << event << endl ;  

  



  t = gAlice->TreeR() ;

  if(t){

    cout << "INFO: AliEMCALGetter::ListBranches -> ****** Recon   : " << endl ; 

    TObjArray * lob = t->GetListOfBranches() ;

    TIter next(lob) ; 

    while ( (branch = static_cast<TBranch*>(next())) )

      cout << "             " << branch->GetName() << " " << branch->GetTitle() << endl ; 

  } else 

    cerr << "WARNING::AliEMCALGetter::ListBranches -> TreeR not found for event " << event << endl ;  

  

}



//____________________________________________________________________________ 

void AliEMCALGetter::NewBranch(TString name, Int_t event)  

{

  fBranchTitle = fSDigitsTitle = fDigitsTitle = fRecPointsTitle = fTrackSegmentsTitle = fRecParticlesTitle =  name ; 

  Event(event) ; 

}



//____________________________________________________________________________ 

Bool_t AliEMCALGetter::NewFile(TString name)  

{

  fHeaderFile = name ; 

  fFile->Close() ; 

  fFailed = kFALSE; 



  fFile = static_cast<TFile*>(gROOT->GetFile(fHeaderFile.Data() ) ) ;

  if(!fFile) {    //if file was not opened yet, read gAlice

    fFile = TFile::Open(fHeaderFile.Data(),"update") ;

    if (!fFile->IsOpen()) {

      cerr << "ERROR : AliEMCALGetter::NewFile -> Cannot open " << fHeaderFile.Data() << endl ; 

      fFailed = kTRUE ;

      return fFailed ;  

    }

    gAlice = static_cast<AliRun *>(fFile->Get("gAlice")) ;

  } 

  

  if (!gAlice) {

    cerr << "ERROR : AliEMCALGetter::AliEMCALGetter -> Cannot find gAlice in " << fHeaderFile.Data() << endl ; 

    fFailed = kTRUE ;

    return fFailed ; 

  }

  return fFailed ; 

}



//____________________________________________________________________________ 

const AliEMCAL * AliEMCALGetter::EMCAL() 

{

  // returns the EMCAL object 

  AliEMCAL * emcal = dynamic_cast<AliEMCAL*>(fModuleFolder->FindObject("EMCAL")) ;  

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

const Bool_t AliEMCALGetter::PostPrimaries(void) const 

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

const Bool_t AliEMCALGetter::PostHits(void) const 

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

const Bool_t AliEMCALGetter::PostSDigits(const char * name, const char * headerFile) const 

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

  if ( !sd ) {

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

    abort() ;

  }    

 

  TFolder * emcalFolder = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject("EMCAL")) ;

  if ( !emcalFolder ) {

    cerr << "ERROR: AliEMCALGetter::Post SRef -> Folder //" << fSDigitsFolder << "/EMCAL/ not found!" << endl;

    abort() ;

  }



  TFolder * emcalSubFolder = 0 ;

  if(file)

    emcalSubFolder = dynamic_cast<TFolder *>(emcalFolder->FindObject(file)) ;

  else

    emcalSubFolder = dynamic_cast<TFolder *>(emcalFolder->FindObject(fHeaderFile)) ;

  

  if(!emcalSubFolder) {

    cerr << "ERROR: AliEMCALGetter::Post SRef -> Folder //Folders/RunMC/Event/Data/EMCAL/" << file << "not found!" << endl;

    abort() ;

  }



  TObject * dis = emcalSubFolder->FindObject(name) ;

  if(!dis) {

    cerr << "ERROR: AliEMCALGetter::DigitesSRef -> object " << name << " not found! " << endl ; 

    abort() ;

  }

  else

    return emcalSubFolder->GetListOfFolders()->GetObjectRef(dis) ;



}



//____________________________________________________________________________ 

const Bool_t AliEMCALGetter::PostSDigitizer(AliEMCALSDigitizer * sdigitizer) const 

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

const Bool_t AliEMCALGetter::PostSDigitizer(const char * name, const char * file) const 

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

const Bool_t AliEMCALGetter::PostDigits(const char * name) const 

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

    abort() ;

  }    

  

  TFolder * emcalFolder  = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("EMCAL")) ; 

  if ( !emcalFolder ) {

    cerr << "ERROR: AliEMCALGetter::DRef -> Folder //" << fDigitsFolder << "/EMCAL/ not found!" << endl;

    abort() ;

  }    



  TObject * d = emcalFolder->FindObject(name) ;

  if(!d) {

    cerr << "ERROR: AliEMCALGetter::DigitsRef -> object " << name << " not found! " << endl ; 

    abort() ;

  }   

  else

    return emcalFolder->GetListOfFolders()->GetObjectRef(d) ;



}



//____________________________________________________________________________ 

const Bool_t AliEMCALGetter::PostDigitizer(AliEMCALDigitizer * digitizer) const 

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

const Bool_t AliEMCALGetter::PostDigitizer(const char * name) const 

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

const Bool_t AliEMCALGetter::PostRecPoints(const char * name) const 

{ // -------------- RecPoints -------------------------------------------

  

  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/TowerRecPoints/name

  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/PreShowerRecPoints/name



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

  TFolder * emcalRPoPreShoFolder  = dynamic_cast<TFolder*>(emcalFolder->FindObject("PreShowerRecPoints")) ;

  if ( !emcalRPoPreShoFolder ) {

    if (fDebug) {

      cout << "WARNING: AliEMCALGetter::Post RPo -> Folder //" << fRecoFolder << "/EMCAL/PreShowerRecPoints/ not found!" << endl;

      cout << "INFO:    AliEMCALGetter::Post Rpo -> Adding Folder //" << fRecoFolder << "/EMCAL/PreShowerRecPoints/" << endl;

    }

    emcalRPoPreShoFolder = emcalFolder->AddFolder("PreShowerRecPoints", "PreSho RecPoints from EMCAL") ;  

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

TObject ** AliEMCALGetter::PreShowerRecPointsRef(const char * name) const 

{ // -------------- RecPoints -------------------------------------------

  

  // the hierarchy is //Folders/Run/Event/RecData/EMCAL/PreShowerRecPoints/name

   

  if ( !fRecoFolder ) {

    cerr << "ERROR: AliEMCALGetter::PreShowerRecPointsRef -> Folder //" << fRecoFolder << " not found!" << endl;

    abort() ; 

  }    



  TFolder * preshoFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/PreShowerRecPoints")) ; 

  if ( !preshoFolder ) {

    cerr << "ERROR: AliEMCALGetter::PreShowerRecPointsRef -> Folder //" << fRecoFolder << "/EMCAL/PreShowerRecPoints/" << endl;

    abort() ;

  }    



  TObject * prp = preshoFolder->FindObject(name ) ;

  if ( !prp )   {

    cerr << "ERROR: AliEMCALGetter::PreShowerRecPointsRef -> Object " << name << " not found! " << endl ; 

    abort() ;

  }

  return preshoFolder->GetListOfFolders()->GetObjectRef(prp) ;



} 



//____________________________________________________________________________ 

const Bool_t AliEMCALGetter::PostClusterizer(AliEMCALClusterizer * clu) const 

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

const Bool_t AliEMCALGetter::PostClusterizer(const char * name) const 

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

  clun+="-v1" ; 

  emcalcl->SetName(clun) ;

  emcalcl->SetTitle(fHeaderFile) ; 

  emcal->Add(emcalcl) ;

  return kTRUE; 

  

}



//____________________________________________________________________________ 

/*const Bool_t AliEMCALGetter::PostTrackSegments(const char * name) const 

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

const Bool_t AliEMCALGetter::PostTrackSegmentMaker(AliEMCALTrackSegmentMaker * tsmaker) const 

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

const Bool_t AliEMCALGetter::PostTrackSegmentMaker(const char * name) const 

{ //------------Track Segment Maker ------------------------------

  

  // the hierarchy is //Folders/Tasks/Reconstructioner/EMCAL/sdigitsname



  TTask * tasks  = dynamic_cast<TTask*>(fTasksFolder->FindObject("Reconstructioner")) ; 

  

  if ( !tasks ) {

    cerr << "ERROR: AliEMCALGetter::Post Ter -> Task //" << fTasksFolder->GetName() << "/Reconstructioner not found!" << endl;

    return kFALSE ;

  }        

  

  TTask * emcal = dynamic_cast<TTask*>(tasks->GetListOfTasks()->FindObject("EMCAL")) ; 

  if ( !emcal )  {

    if (fDebug) {

      cout <<"WARNING: AliEMCALGetter::Post Ter -> //" << fTasksFolder->GetName() << "/Reconstructioner/EMCAL not found!" << endl; 

      cout <<"INFO: AliEMCALGetter::Post Ter -> Adding //" << fTasksFolder->GetName() << "/Reconstructioner/EMCAL" << endl;

    }

    emcal = new TTask("EMCAL", "") ; 

    tasks->Add(emcal) ; 

  } 



  TList * l = emcal->GetListOfTasks() ;   

  TIter it(l) ;

  TString tsn(name);

  tsn+=":tsm" ; 

  TTask * task ;

  while((task = static_cast<TTask *>(it.Next()) )){

    TString taskname(task->GetName()) ;

    if(taskname.BeginsWith(tsn))

      return kTRUE ;

  }

  

  AliEMCALTrackSegmentMakerv1 * emcalts = new AliEMCALTrackSegmentMakerv1() ;

  tsn+="-v1" ;

  emcalts->SetName(tsn) ;

  emcalts->SetTitle(fHeaderFile) ;

  emcal->Add(emcalts) ;      

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

const Bool_t AliEMCALGetter::PostRecParticles(const char * name) const 

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

{ // ---------------RecParticles -----------------------------------

  

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

const Bool_t AliEMCALGetter::PostPID(AliEMCALPID * pid) const 

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

const Bool_t AliEMCALGetter::PostPID(const char * name) const 

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

  emcalpid->SetTitle(fHeaderFile) ;

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

const Bool_t AliEMCALGetter::PostQA(void) const 

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

TTree * AliEMCALGetter::TreeK(TString filename)  

{



  // returns TreeK from file filename

  // usefull in case of split file



  if ( filename.IsNull() ) 

    filename = fHeaderFile ; 



  TFile * file = 0 ; 

  file = static_cast<TFile*>(gROOT->GetFile(filename.Data() ) ) ;

  if (!file)  {  // file not yet open 

    file = TFile::Open(filename.Data(), "read") ; 

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

const TParticle * AliEMCALGetter::Primary(Int_t index) const

{

  // Return primary particle numbered by <index>

  

  if(index < 0) 

    return 0 ;

  TParticle *  p = 0 ;

  p = gAlice->Particle(index) ; 
  

  return p ; 

      

}



//____________________________________________________________________________ 

const TParticle * AliEMCALGetter::Secondary(TParticle* p, Int_t index) const

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

Int_t AliEMCALGetter::ReadTreeD(const Int_t event)

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

	cout << "WARNING: AliEMCALGetter::ReadTreeD -> Cannot find TreeD in " << fDigitsFileName.Data() << endl;

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

    return 2; 

  }   

 

  //read digits

  if(!Digits(fDigitsTitle) ) 

    PostDigits(fDigitsTitle);

  digitsbranch->SetAddress(DigitsRef(fDigitsTitle)) ;

  digitsbranch->GetEntry(0) ;

  

  

  // read  the Digitizer

  RemoveTask("D", fDigitsTitle) ; // I do not understand why I need that 

  if(!Digitizer(fDigitsTitle))

    PostDigitizer(fDigitsTitle) ;

  digitizerbranch->SetAddress(DigitizerRef(fDigitsTitle)) ;

  digitizerbranch->GetEntry(0) ;

 

  lob  ->Delete();

  if(gAlice->TreeD()!=treeD)

    treeD->Delete();

  return 0 ; 

}



//____________________________________________________________________________ 

Int_t AliEMCALGetter::ReadTreeH()

{

  // Read the first entry of EMCAL branch in hit tree gAlice->TreeH()



  TTree * treeH = gAlice->TreeH() ;

  

  if(!treeH) {// TreeH not found in header file

    

    if (fDebug) 

      cout <<   "WARNING: AliEMCALGetter::ReadTreeH -> Cannot find TreeH in " << fHeaderFile << endl ;

    

    TString searchFileName("EMCAL.HITS") ; 

    if((strcmp(fBranchTitle.Data(),"Default")!=0)&&(strcmp(fBranchTitle.Data(),"")!=0)){

      searchFileName+="." ;

      searchFileName += fBranchTitle ;

    }

    searchFileName+=".root" ;

    

    if ( (treeH = TreeH(searchFileName)) ) { //found TreeH in the file which contains the hits

      if (fDebug) 

	cout << "INFO: AliEMCALGetter::ReadTreeH -> TreeH found in " << searchFileName.Data() << endl ; 

      

    } else {

      cerr << "ERROR: AliEMCALGetter::ReadTreeH -> TreeH not found " << endl ; 

      return 1;

    }  

  }

  

  TBranch * hitsbranch = static_cast<TBranch*>(gAlice->TreeH()->GetBranch("EMCAL")) ;

  if ( !hitsbranch ) {

    if (fDebug)

      cout << "WARNING:  AliEMCALGetter::ReadTreeH -> Cannot find branch EMCAL" << endl ; 

    return 2;

  }

  if(!Hits())

    PostHits() ;

 

  if (hitsbranch->GetEntries() > 1 ) {

    (dynamic_cast<TClonesArray*> (*HitsRef()))->Clear() ;

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

  (dynamic_cast<TClonesArray*> (*HitsRef()))->Clear() ;

  hitsbranch->SetAddress(HitsRef()) ;

  hitsbranch->GetEntry(0) ;

  }

  return 0 ; 

}



//____________________________________________________________________________ 

void AliEMCALGetter::Track(const Int_t itrack)

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



  (dynamic_cast<TClonesArray*> (*HitsRef()))->Clear() ;

  hitsbranch->SetAddress(HitsRef()) ;

  hitsbranch->GetEntry(itrack) ;

  

}



//____________________________________________________________________________ 

void AliEMCALGetter::ReadTreeQA()

{

  if (fDebug)

    cout << "WARNING : " << ClassName() << "::ReadTreeQA -> not implemented" << endl ; 

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

  

}

  

//____________________________________________________________________________ 

Int_t AliEMCALGetter::ReadTreeR(const Int_t event)

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

  // See AliEMCALPIDv1    



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

	cout << "WARNING: AliEMCALGetter::ReadTreeD -> Cannot find TreeR in " << fRecPointsFileName.Data() << endl;

      return 1;

    }

  }

  else

    treeR = gAlice->TreeR() ;

  

  // RecPoints 

  TObjArray * lob = static_cast<TObjArray*>(treeR->GetListOfBranches()) ;

  TIter next(lob) ; 

  TBranch * branch = 0 ; 

  TBranch * towerbranch = 0 ; 

  TBranch * preshowerbranch = 0 ; 

  TBranch * clusterizerbranch = 0 ; 

  Bool_t emcalTowerRPfound = kFALSE, emcalPreShoRPfound = kFALSE, clusterizerfound = kFALSE ; 



  

  while ( (branch = static_cast<TBranch*>(next())) && (!emcalTowerRPfound || !emcalPreShoRPfound || !clusterizerfound) ) {

    if(strcmp(branch->GetTitle(), fRecPointsTitle)==0 ) {

      if ( strcmp(branch->GetName(), "EMCALTowerRP")==0) {

	towerbranch = branch ; 

	emcalTowerRPfound = kTRUE ;

      }

      else if ( strcmp(branch->GetName(), "EMCALPreShoRP")==0) {

	preshowerbranch = branch ; 

	emcalPreShoRPfound = kTRUE ;

      }

      else if(strcmp(branch->GetName(), "AliEMCALClusterizer")==0){

	clusterizerbranch = branch ; 

	clusterizerfound = kTRUE ; 

      }

    }

  }



  if ( !emcalTowerRPfound || !emcalPreShoRPfound || !clusterizerfound) {

    if (fDebug)

      cout << "WARNING: AliEMCALGetter::ReadTreeR -> Cannot find RecPoints and/or Clusterizer with name " 

	   << fRecPointsTitle << endl ;

 

  } else { 

    if(!TowerRecPoints(fRecPointsTitle) ) 

      PostRecPoints(fRecPointsTitle) ;

    

    towerbranch->SetAddress(TowerRecPointsRef(fRecPointsTitle)) ;

    towerbranch->GetEntry(0) ;



    preshowerbranch->SetAddress(PreShowerRecPointsRef(fRecPointsTitle)) ; 

    preshowerbranch->GetEntry(0) ;  

    

    if(!Clusterizer(fRecPointsTitle) )

      PostClusterizer(fRecPointsTitle) ;

    

    clusterizerbranch->SetAddress(ClusterizerRef(fRecPointsTitle)) ;

    clusterizerbranch->GetEntry(0) ;

  }

  

//   //------------------- TrackSegments ---------------------

//   next.Reset() ; 

//   TBranch * tsbranch = 0 ; 

//   TBranch * tsmakerbranch = 0 ; 

//   Bool_t emcaltsfound = kFALSE, tsmakerfound = kFALSE ; 

//   while ( (branch = static_cast<TBranch*>(next())) && (!emcaltsfound || !tsmakerfound) ) {

//     if(strcmp(branch->GetTitle(), fTrackSegmentsTitle)==0 )  {

//       if ( strcmp(branch->GetName(), "EMCALTS")==0){

// 	tsbranch = branch ; 

// 	emcaltsfound = kTRUE ;

//       }

//       else if(strcmp(branch->GetName(), "AliEMCALTrackSegmentMaker")==0) {

// 	tsmakerbranch = branch ; 

// 	tsmakerfound  = kTRUE ; 

//       }

//     }

//   }



//   if ( !emcaltsfound || !tsmakerfound ) {

//     if (fDebug)

//       cout << "WARNING: AliEMCALGetter::ReadTreeR -> Cannot find TrackSegments and/or TrackSegmentMaker with name "

// 	   << fTrackSegmentsTitle << endl ;

//   } else { 

//     // Read and Post the TrackSegments

//     if(!TrackSegments(fTrackSegmentsTitle))

//       PostTrackSegments(fTrackSegmentsTitle) ;

//     tsbranch->SetAddress(TrackSegmentsRef(fTrackSegmentsTitle)) ;

//     tsbranch->GetEntry(0) ;



//     // Read and Post the TrackSegment Maker

//     if(!TrackSegmentMaker(fTrackSegmentsTitle))

//       PostTrackSegmentMaker(fTrackSegmentsTitle) ;

//     tsmakerbranch->SetAddress(TSMakerRef(fTrackSegmentsTitle)) ;

//     tsmakerbranch->GetEntry(0) ;

//  }

  

  

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

//   } else { 

//     // Read and Post the RecParticles

//     if(!RecParticles(fRecParticlesTitle)) 

//       PostRecParticles(fRecParticlesTitle) ;

//     rpabranch->SetAddress(RecParticlesRef(fRecParticlesTitle)) ;

//     rpabranch->GetEntry(0) ;

//     // Read and Post the PID

//     if(!PID(fRecParticlesTitle))

//       PostPID(fRecParticlesTitle) ;

//     pidbranch->SetAddress(PIDRef(fRecParticlesTitle)) ;

//     pidbranch->GetEntry(0) ;

//   }



  if(gAlice->TreeR()!=treeR)

    treeR->Delete();

  return 0 ; 

}



//____________________________________________________________________________ 

Int_t AliEMCALGetter::ReadTreeS(Int_t event)

{

  // Reads the SDigits treeS from all files  

  // Files, which should be opened are listed in emcalF

  // So, first get list of files

  TFolder * emcalF = dynamic_cast<TFolder *>(fSDigitsFolder->FindObject("EMCAL")) ;

  if (!emcalF) 

    emcalF = fSDigitsFolder->AddFolder("EMCAL", "SDigits from EMCAL") ; 

  TCollection * folderslist = emcalF->GetListOfFolders() ; 

  

  // Now iterate over the list of files and read TreeS into Whiteboard

  TIter next(folderslist) ; 

  TFolder * folder = 0 ; 

  TFile * file; 

  TTree * treeS = 0;

  while ( (folder = static_cast<TFolder*>(next())) ) {

    TString fileName("") ;

    if(fToSplit)

      fileName = folder->GetTitle() ;

    else

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

	cout << "WARNING: AliEMCALGetter::ReadTreeS -> Cannot find TreeS in " << fileName.Data() << endl;

      return 1;

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

      

      else if ( (strcmp(branch->GetName(), "AliEMCALSDigitizer")==0) && 

		(strcmp(branch->GetTitle(), fSDigitsTitle)==0) ) {

	sdigitizerfound = kTRUE ; 

	sdigitizerBranch = branch ;

      }

    }

    if ( !emcalfound || !sdigitizerfound ) {

      if (fDebug)

	cout << "WARNING: AliEMCALDigitizer::ReadSDigits -> Digits and/or Digitizer branch with name " << GetName() 

	     << " not found" << endl ;

      return 2; 

    }   

    

    if ( !folder->FindObject(fSDigitsTitle) )  

      PostSDigits(fSDigitsTitle,folder->GetName()) ;



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

  if(gAlice->TreeS()!=treeS)

    treeS->Delete();

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



    cerr << "ERROR: AliEMCALGetter::ReadPrimaries -> TreeK not  found " << endl ; 

    return ;

    

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

 

  if( strstr(opt,"P") )

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



  if ( what.CompareTo("Primaries") == 0 ) {

    folder = dynamic_cast<TFolder *>(fPrimariesFolder->FindObject("Primaries")) ; 

    if (folder) 

      emcalO  = dynamic_cast<TObject *>(folder->FindObject("Primaries")) ;  

    else 

      return 0 ; 

  }

  else if ( what.CompareTo("Hits") == 0 ) {

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

  else if ( what.CompareTo("PreShowerRecPoints") == 0 ) {

    folder = dynamic_cast<TFolder *>(fRecoFolder->FindObject("EMCAL/PreShowerRecPoints")) ; 

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

    emcal = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/PreShowerRecPoints")) ;

    if (!emcal) 

      return ;

  }



  else if (opt == "T") { // TrackSegments

    emcal = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/TrackSegments")) ;

    if (!emcal) 

      return ;

  }



  else if (opt == "P") { // RecParticles

    emcal = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/RecParticles")) ;

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



//____________________________________________________________________________ 

void AliEMCALGetter::CleanWhiteBoard(void){



  TFolder * emcalmain = 0 ; 

  TFolder * emcal ;

  TObjArray * ar ;

  TList * lofTasks = 0 ; 

  TTask * task = 0 ; 

  TTask * emcalt = 0 ; 

  

  // Hits  

  emcal = dynamic_cast<TFolder*>(fHitsFolder->FindObject("EMCAL")) ;

  if (emcal){  

    TObjArray * ar  = dynamic_cast<TObjArray*>(emcal->FindObject("Hits")) ; 

    if (ar) { 

      emcal->Remove(ar) ;

      ar->Delete() ; 

      delete ar ; 

    }

  }

  

  // SDigits

  emcalmain = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("EMCAL")) ;

  if (emcalmain){ 

    emcal = dynamic_cast<TFolder*>(emcalmain->FindObject(fHeaderFile)) ;

    if (emcal) {

      ar  = dynamic_cast<TObjArray*>(emcal->FindObject(fSDigitsTitle)) ; 

      if (ar) { 

	emcal->Remove(ar) ;

	ar->Delete() ; 

	delete ar ; 

      }

    }

    emcalmain->Remove(emcal) ; 

  }



  

  // Digits

  emcal = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("EMCAL")) ;

  if (emcal){ 

    ar  = dynamic_cast<TObjArray*>(emcal->FindObject(fDigitsTitle)) ; 

    if (ar) { 

      emcal->Remove(ar) ;

      ar->Delete() ; 

      delete ar ; 

    }

  }





  // TowerRecPoints

  emcal = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/TowerRecPoints")) ;

  if (emcal){ 

    ar  = dynamic_cast<TObjArray*>(emcal->FindObject(fRecPointsTitle)) ; 

    if (ar) { 

      emcal->Remove(ar) ;

      ar->Delete() ; 

      delete ar ; 

    }

  }



  

  // PreShowerRecPoints

  emcal = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/PreShowerRecPoints")) ;

  if (emcal){ 

    ar  = dynamic_cast<TObjArray*>(emcal->FindObject(fRecPointsTitle)) ; 

    if (ar) { 

      emcal->Remove(ar) ;

      ar->Delete() ; 

      delete ar ; 

    }

  }  



  

  // TrackSegments

  emcal = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/TrackSegments")) ;

  if (emcal) { 

    ar  = dynamic_cast<TObjArray*>(emcal->FindObject(fTrackSegmentsTitle)) ; 

    if (ar) { 

      emcal->Remove(ar) ;

      ar->Delete() ; 

      delete ar ; 

    }

  }

  





  // RecParticles

  emcal = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL/RecParticles")) ;

  if (emcal){ 

    ar  = dynamic_cast<TObjArray*>(emcal->FindObject(fRecParticlesTitle)) ; 

    if (ar) { 

      emcal->Remove(ar) ;

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

    emcalt =  dynamic_cast<TTask*>(task->GetListOfTasks()->FindObject("EMCAL")) ;

    if (emcalt){

      lofTasks = emcalt->GetListOfTasks() ;

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

    emcalt =  dynamic_cast<TTask*>(task->GetListOfTasks()->FindObject("EMCAL")) ;

    if (emcalt){

      lofTasks = emcalt->GetListOfTasks() ;

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

    emcalt =  dynamic_cast<TTask*>(task->GetListOfTasks()->FindObject("EMCAL")) ;

    if (emcalt){

      lofTasks = emcalt->GetListOfTasks() ;

      if (lofTasks){ 

	obj = lofTasks->FindObject(sdname.Data()) ; 

	if (obj) 

	  lofTasks->Remove(obj) ;

      }

    }

  }  



}

//____________________________________________________________________________ 

void AliEMCALGetter::SetTitle(const char * branchTitle ) 

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

    fSDigitsFileName      += "EMCAL.SDigits." ;

    fDigitsFileName       += "EMCAL.Digits." ; 

    fRecPointsFileName    += "EMCAL.RecData." ; 

    fTrackSegmentsFileName+= "EMCAL.RecData." ; 

    fRecParticlesFileName += "EMCAL.RecData." ; 

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

    fSDigitsFileName       = "" ; 

    fDigitsFileName        = "" ; 

    fRecPointsFileName     = "" ; 

    fRecParticlesFileName  = "" ; 

    fTrackSegmentsFileName = "" ; 

  }

  TFolder * emcalFolder ; 

  emcalFolder = dynamic_cast<TFolder*>(fHitsFolder->FindObject("EMCAL")) ; 

  if ( !emcalFolder ) 

    emcalFolder = fHitsFolder->AddFolder("EMCAL", "Hits from EMCAL") ; 



  emcalFolder = dynamic_cast<TFolder*>(fSDigitsFolder->FindObject("EMCAL")) ;

  if ( !emcalFolder ) 

    emcalFolder = fSDigitsFolder->AddFolder("EMCAL", "SDigits from EMCAL") ; 

  

  //Make folder for SDigits

  TString subdir(fHeaderFile) ;

  subdir.ReplaceAll("/","_") ;

  emcalFolder->AddFolder(subdir, fSDigitsFileName.Data());





  emcalFolder  = dynamic_cast<TFolder*>(fDigitsFolder->FindObject("EMCAL")) ;

  if ( !emcalFolder ) 

    emcalFolder = fDigitsFolder->AddFolder("EMCAL", "Digits from EMCAL") ;  



  emcalFolder  = dynamic_cast<TFolder*>(fRecoFolder->FindObject("EMCAL")) ; 

  if ( !emcalFolder )

    emcalFolder = fRecoFolder->AddFolder("EMCAL", "Reconstructed data from EMCAL") ;  

  



}

//____________________________________________________________________________ 

void AliEMCALGetter::CloseSplitFiles(void){

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

