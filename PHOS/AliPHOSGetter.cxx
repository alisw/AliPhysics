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

// --- Standard library ---
#include <iostream.h>

// --- AliRoot header files ---

#include "AliRun.h"
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

  if ( fHeaderFile != "aliroot"  ) { // to call the getter without a file

    //open headers file
    TFile * file = (TFile*) gROOT->GetFile(fHeaderFile.Data() ) ;
    
    if(file == 0){    //if file was not opened yet, read gAlice
      if(fHeaderFile.Contains("rfio")) // if we read file using HPSS
	file =	TFile::Open(fHeaderFile.Data(),"update") ;
      else
	file = new TFile(fHeaderFile.Data(),"update") ;
      
      if (!file->IsOpen()) {
	cerr << "ERROR : AliPHOSGetter::AliPHOSGetter -> Cannot open " << fHeaderFile.Data() << endl ; 
	abort() ; 
      }
      
      gAlice = (AliRun *) file->Get("gAlice") ;
      
      if (!gAlice) {
	cerr << "ERROR : AliPHOSGetter::AliPHOSGetter -> Cannot find gAlice in " << fHeaderFile.Data() << endl ; 
	abort() ; 
      }
    }
  }
  fDebug=0;
}
//____________________________________________________________________________ 
AliPHOSGetter::~AliPHOSGetter(){
  //Here we remove all TFolders and TTasks with current title and file name
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  if(aliceF){
    //Hits:  the hierarchy is //YSALICE/WhiteBoard/Hits/PHOS/...
    TFolder * hitsF  = (TFolder*)aliceF->FindObject("WhiteBoard/Hits/PHOS") ; 
    if(hitsF){
      TObject * h = hitsF->FindObject("hits") ;
      hitsF->Remove(h) ;
    }
    
    //SDigits:  the hierarchy is //YSALICE/WhiteBoard/SDigits/PHOS/...
    TFolder * sdigitsF  = (TFolder*)aliceF->FindObject("WhiteBoard/SDigits/PHOS") ; 
    if(sdigitsF){
      TCollection* l = sdigitsF->GetListOfFolders() ;
      TIter it(l) ;
      TObject * sf ;
      while((sf = it.Next()) )
	sdigitsF->RecursiveRemove(sf) ;
    }
  }
}
//____________________________________________________________________________ 
void AliPHOSGetter::CreateWhiteBoard() const
{
  // Posts a few item to the white board (folders)
  
  // -- the geometry
  if(!PostGeometry() ) abort() ; 
  
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
  fgObjGetter->CreateWhiteBoard() ;
    
  return fgObjGetter ; 
  
}

//____________________________________________________________________________ 
 const  AliPHOS * AliPHOSGetter::PHOS() const 
{
  // returns the PHOS object 
  return ( (AliPHOS*)gAlice->GetDetector("PHOS") ); 
}  

//____________________________________________________________________________ 
  const AliPHOSGeometry *  AliPHOSGetter::PHOSGeometry() const 
{
  // retrieves the geometr from the folder

  TString path("YSAlice/WhiteBoard/Geometry/PHOS/") ;
  path += PHOS()->GetTitle() ;
  return (AliPHOSGeometry*)gROOT->FindObjectAny(path.Data()) ; 
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostGeometry() const 
{  //--------Geometry --------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/WhiteBoard/Geometry/PHOS/
  TFolder * geomF  = (TFolder*)aliceF->FindObject("WhiteBoard/Geometry/PHOS") ; 
  if ( !geomF ) {
    cerr << "ERROR: AliPHOSGetter::Post G -> Folder WhiteBoard/Geometry/PHOS/" << " not found!" << endl;
    return kFALSE ;
  }    
  else {
    AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance(PHOS()->GetTitle(),"") ;  
    geomF->Add((TObject*)geom) ; 
  }
  return kTRUE;   
}
//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostHits(void) const 
{  //------- Hits ----------------------

  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/WhiteBoard/SDigits/PHOS/Hits
  TFolder * hitsF  = (TFolder*)aliceF->FindObject("WhiteBoard/Hits/PHOS") ; 
  if ( !hitsF ) {
    cerr << "ERROR: AliPHOSGetter::Post H -> Folder WhiteBoard/Hits/PHOS/" << " not found!" << endl;
    return kFALSE ;
  }    
 
  TObject * h = hitsF->FindObject("Hits") ;

  if(!h){
    TClonesArray *hits=  new TClonesArray("AliPHOSHit",1000) ;
    hits->SetName("Hits") ;
    hitsF->Add(hits) ; 
  }
 
  return kTRUE;
} 
//____________________________________________________________________________ 
TClonesArray ** AliPHOSGetter::HitsRef(void) const 
{  //------- Hits ----------------------

  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/WhiteBoard/SDigits/PHOS/Hits
  TFolder * hitsF  = (TFolder*)aliceF->FindObject("WhiteBoard/Hits/PHOS") ; 
  if ( !hitsF ) {
    cerr << "ERROR: AliPHOSGetter::Post H -> Folder WhiteBoard/Hits/PHOS/" << " not found!" << endl;
    return 0;
  }    
 
  TObject * h = hitsF->FindObject("Hits") ;
  if(!h)
    return 0 ;
  else
    return (TClonesArray **) hitsF->GetListOfFolders()->GetObjectRef(h) ;

}
//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostSDigits(const char * name, const char * headerFile) const 
{  //---------- SDigits -------------------------

  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/WhiteBoard/SDigits/PHOS/headerFile/sdigitsname
  // because you can have sdigits from several hit files for mixing

  TFolder * sdigitsF  = (TFolder*)aliceF->FindObject("WhiteBoard/SDigits/PHOS") ;
  TString subdir(headerFile) ;
  TFolder * sdigitsF2 = (sdigitsF2=(TFolder*)sdigitsF->FindObject(subdir)) ; 
  if ( !sdigitsF2 ) 
    sdigitsF2 = sdigitsF->AddFolder(subdir, ""); 

  TObject * sd  = sdigitsF2->FindObject(name ); 
  if ( sd ) {
    if (fDebug)
      cerr <<"INFO: AliPHOSGetter::Post S -> Folder " << subdir 
	   << " already exists!" << endl ;  
  }else{
    TClonesArray * sdigits = new TClonesArray("AliPHOSDigit",1000) ;
    sdigits->SetName(name) ;
    sdigitsF2->Add(sdigits) ;
  }
  
  return kTRUE;
} 
//____________________________________________________________________________ 
TClonesArray ** AliPHOSGetter::SDigitsRef(const char * name, const char * file) const 
{  //------- Hits ----------------------

  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/WhiteBoard/SDigits/PHOS/Hits
  TFolder * sdisF  = (TFolder*)aliceF->FindObject("WhiteBoard/SDigits/PHOS") ; 
  if ( !sdisF ) {
    cerr << "ERROR: AliPHOSGetter::SDRef -> Folder WhiteBoard/SDigits/PHOS/" << " not found!" << endl;
    return 0;
  }    
  TFolder * fileF ;
  if(file)
    fileF = (TFolder *) sdisF->FindObject(file) ;
  else
    fileF = (TFolder *) sdisF->FindObject(fHeaderFile) ;

  if(!fileF)
    abort() ;

  TObject * dis = fileF->FindObject(name) ;
  if(!dis)
    return 0 ;
  else
    return (TClonesArray **) fileF->GetListOfFolders()->GetObjectRef(dis) ;

}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostSDigitizer(AliPHOSSDigitizer * sdigitizer) const 
{  //---------- SDigitizer -------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/tasks/Digitizer/PHOS/sdigitsname
  TTask * sdigitsF  = (TTask*)aliceF->FindObject("tasks/SDigitizer") ; 
  if ( !sdigitsF ) {
    cerr << "ERROR: AliPHOSGetter::Post Ser -> Task tasks/SDigitizer" << " not found!" << endl;
    return kFALSE ;
  }        
  TTask * phos = (TTask*)sdigitsF->GetListOfTasks()->FindObject("PHOS") ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Ser -> tasks/SDigitizer/PHOS" << " not found!" << endl; 
    return kFALSE ; 
  } 
  AliPHOSSDigitizer * phossd  = (AliPHOSSDigitizer *) phos->GetListOfTasks()->FindObject( sdigitizer->GetName() ); 
  if (phossd) { 
    if (fDebug)
      cout << "INFO: AliPHOSGetter::Post Ser -> Task " << sdigitizer->GetName() << " already exists" << endl ; 
    phos->GetListOfTasks()->Remove(phossd) ;
  }
  phos->Add(sdigitizer) ;	
  return kTRUE; 
  
}
//____________________________________________________________________________ 
AliPHOSSDigitizer ** AliPHOSGetter::SDigitizerRef(const char * name) const 
{  

  TString path("tasks/SDigitizer") ;
 
  TFolder * aliceF  = (TFolder*)gROOT ->FindObjectAny("YSAlice") ; 
  TTask   * aliceT  = (TTask*)  aliceF->FindObject(path) ; 

  if (!aliceT) {
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << path << " not found!" << endl ;  
    abort() ;
  }

  TTask   * phosT   = (TTask*)  aliceT->GetListOfTasks()->FindObject("PHOS") ; 
  if (!phosT) { 
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << path << "/PHOS not found!" << endl ;  
    abort() ;
  }
  TList * l = phosT->GetListOfTasks() ; 

  TTask * task = (TTask*)l->FindObject(name) ; 

  return (AliPHOSSDigitizer **) l->GetObjectRef(task) ;

}


//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostSDigitizer(const char * name, const char * file) const 
{  //---------- SDigitizer -------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/tasks/Digitizer/PHOS/sdigitsname
  TTask * sdigitsF  = (TTask*)aliceF->FindObject("tasks/SDigitizer") ; 
  if ( !sdigitsF ) {
    cerr << "ERROR: AliPHOSGetter::Post Ser -> Task tasks/SDigitizer" << " not found!" << endl;
    return kFALSE ;
  }        
  TTask * phos = (TTask*)sdigitsF->GetListOfTasks()->FindObject("PHOS") ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Ser -> tasks/SDigitizer/PHOS" << " not found!" << endl; 
    return kFALSE ; 
  } 
  TString sdname(name) ;
  sdname.Append(":") ;
  sdname.Append(file);
  AliPHOSSDigitizer * phossd  = (AliPHOSSDigitizer *) phos->GetListOfTasks()->FindObject( sdname ); 
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

  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/WhiteBoard/Digits/PHOS/digitsname
  TFolder * digitsF  = (TFolder*)aliceF->FindObject("WhiteBoard/Digits/PHOS") ;
  if ( !digitsF ) {
    cerr << "ERROR: AliPHOSGetter::Post D -> Folder WhiteBoard/Digits/PHOS/" << " not found!" << endl;
    return kFALSE ; 
  }    
  digitsF->SetTitle("") ; 
  TObject*  dig = digitsF->FindObject( name ) ;
  if ( !dig ) {
    TClonesArray * digits = new TClonesArray("AliPHOSDigit",1000) ;
    digits->SetName(name) ;
    digitsF->Add(digits) ;  
  }
  return kTRUE; 
}
//____________________________________________________________________________ 
TClonesArray ** AliPHOSGetter::DigitsRef(const char * name) const 
{  

  TFolder * digitsF  = (TFolder*)gROOT->FindObjectAny("YSAlice/WhiteBoard/Digits/PHOS") ; 
  
   if ( !digitsF ) {
    cerr << "ERROR: AliPHOSGetter::DRef -> Folder WhiteBoard/Digits/PHOS/" << " not found!" << endl;
    return 0;
  }    

  TObject * d = digitsF->FindObject(name) ;
  if(!d)
    return 0 ;
  else
    return (TClonesArray **) digitsF->GetListOfFolders()->GetObjectRef(d) ;

}

 
//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostDigitizer(AliPHOSDigitizer * digitizer) const 
{  //---------- Digitizer -------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/tasks/Digitizer/PHOS/digitsname
  TTask * digitsF  = (TTask*)aliceF->FindObject("tasks/Digitizer") ; 
  if ( !digitsF ) {
    cerr << "ERROR: AliPHOSGetter::Post Der -> Task tasks/Digitizer" << " not found!" << endl;
    return kFALSE ; 
  }        
  TTask * phos = (TTask*)digitsF->GetListOfTasks()->FindObject("PHOS") ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Der -> tasks/Digitizer/PHOS" << " not found!" << endl; 
    return kFALSE ; 
  } else {
    AliPHOSDigitizer * phosd = (AliPHOSDigitizer*)phos->GetListOfTasks()->FindObject(digitizer->GetName()) ; 
    if (phosd) { 
      phosd->Delete() ;
      phos->GetListOfTasks()->Remove(phosd) ;
    }
    phos->Add((TTask*)digitizer) ; 
    return kTRUE; 
  } 
} 
//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostDigitizer(const char * name) const 
{  //---------- Digitizer -------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/tasks/Digitizer/PHOS/digitsname
  TTask * digitsF  = (TTask*)aliceF->FindObject("tasks/Digitizer") ; 
  if ( !digitsF ) {
    cerr << "ERROR: AliPHOSGetter::Post Der -> Task tasks/Digitizer" << " not found!" << endl;
    return kFALSE ; 
  }        
  TTask * phos = (TTask*)digitsF->GetListOfTasks()->FindObject("PHOS") ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Der -> tasks/Digitizer/PHOS" << " not found!" << endl; 
    return kFALSE ; 
  }

  AliPHOSDigitizer * phosd = (AliPHOSDigitizer*)phos->GetListOfTasks()->FindObject(name) ; 
  if (!phosd) { 
    phosd = new AliPHOSDigitizer() ;
    phosd->SetName(fDigitsTitle) ;
    phosd->SetTitle(fHeaderFile) ;
    phos->Add(phosd) ;
  }
  return kTRUE;  
}

//____________________________________________________________________________ 
AliPHOSDigitizer ** AliPHOSGetter::DigitizerRef(const char * name) const 
{  

  TString path("tasks/Digitizer") ;
 
  TFolder * aliceF  = (TFolder*)gROOT ->FindObjectAny("YSAlice") ; 
  TTask   * aliceT  = (TTask*)  aliceF->FindObject(path) ; 

  if (!aliceT) {
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << path << " not found!" << endl ;  
    abort() ;
  }

  TTask   * phosT   = (TTask*)  aliceT->GetListOfTasks()->FindObject("PHOS") ; 
  if (!phosT) { 
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << path << "/PHOS not found!" << endl ;  
    abort() ;
  }
  TList * l = phosT->GetListOfTasks() ; 

  TTask * task = (TTask*)l->FindObject(name) ; 

  return (AliPHOSDigitizer **) l->GetObjectRef(task) ;

}
 
//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostRecPoints(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/WhiteBoard/RecPoints/PHOS/emc/recpointsname  
  TFolder * emcrpF  = (TFolder*)aliceF->FindObject("WhiteBoard/RecPoints/PHOS/emc") ; 
  
  if ( !emcrpF ) {
    cerr << "ERROR: AliPHOSGetter::Post R -> Folder WhiteBoard/RecPoints/PHOS/emc" 
	 << " not found!" << endl;
    return kFALSE ; 
  }    
  emcrpF->SetTitle("") ;
  TObject * erp = emcrpF->FindObject(name ) ;
  if ( !erp )   {
    TObjArray * emcrp = new TObjArray(100) ;
    emcrp->SetName(name) ;
    emcrpF->Add(emcrp) ;  
  }

  TFolder * cpvrpF  = (TFolder*)aliceF->FindObject("WhiteBoard/RecPoints/PHOS/cpv") ; 
  
  if ( !cpvrpF ) {
    cerr << "ERROR: AliPHOSGetter::Post R -> Folder WhiteBoard/RecPoints/PHOS/cpv" 
	 << " not found!" << endl;
    return kFALSE ; 
  }    
  cpvrpF->SetTitle("") ; 
  TObject * crp =  cpvrpF->FindObject( name ) ;
  if ( !crp )   {
    TObjArray * cpvrp = new TObjArray(100) ;
    cpvrp->SetName(name) ;
    cpvrpF->Add(cpvrp) ;  
  }
  return kTRUE; 
}

//____________________________________________________________________________ 
TObjArray ** AliPHOSGetter::EmcRecPointsRef(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  TFolder * emcrpF  = (TFolder*)gROOT->FindObjectAny("YSAlice/WhiteBoard/RecPoints/PHOS/emc") ; 
   
  if ( !emcrpF ) {
    cerr << "ERROR: AliPHOSGetter::EmcRecPointsRef -> Folder WhiteBoard/RecPoints/PHOS/emc" 
	 << " not found!" << endl;
    return 0 ; 
  }    

  TObject * erp = emcrpF->FindObject(name ) ;
  if ( !erp )   {
    return 0 ;
  }
  return (TObjArray **) emcrpF->GetListOfFolders()->GetObjectRef(erp) ;

} 
//____________________________________________________________________________ 
TObjArray ** AliPHOSGetter::CpvRecPointsRef(const char * name) const 
{ // -------------- RecPoints -------------------------------------------
  
  TFolder * cpvrpF  = (TFolder*)gROOT->FindObjectAny("YSAlice/WhiteBoard/RecPoints/PHOS/cpv") ; 
   
  if ( !cpvrpF ) {
    cerr << "ERROR: AliPHOSGetter::CpvRecPointsRef -> Folder WhiteBoard/RecPoints/PHOS/cpv" 
	 << " not found!" << endl;
    return 0 ; 
  }    
  TObject * crp = cpvrpF->FindObject(name ) ;
  if ( !crp )   {
    return 0 ;
  }
  return (TObjArray **) cpvrpF->GetListOfFolders()->GetObjectRef(crp) ;

} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostClusterizer(AliPHOSClusterizer * clu) const 
{ // ------------------ AliPHOSClusterizer ------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/recpointsname
  TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
  if ( !reF ) {
    cerr << "ERROR: AliPHOSGetter::Post Rer -> Task tasks/Reconstructioner" << " not found!" << endl;
    return kFALSE ; 
  } 
       
  TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
  
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Rer -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
    return kFALSE ; 
  } else {
    AliPHOSClusterizer * phoscl = (AliPHOSClusterizer*)phos->GetListOfTasks()->FindObject(clu->GetName()) ; 
    if (phoscl) { 
      if (fDebug)
	cout << "INFO: AliPHOSGetter::Post Rer -> Task " << clu->GetName() << " already exists" << endl ; 
      phos->GetListOfTasks()->Remove(phoscl) ;
    }
    phos->Add(clu) ;      
    return kTRUE; 
  } 
}

//____________________________________________________________________________ 
AliPHOSClusterizer ** AliPHOSGetter::ClusterizerRef(const char * name) const 
{ // ------------------ AliPHOSClusterizer ------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/recpointsname
  TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
  if ( !reF ) {
    cerr << "ERROR: AliPHOSGetter::Post Rer -> Task tasks/Reconstructioner" << " not found!" << endl;
    return kFALSE ; 
  } 
       
  TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
  
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Rer -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
    return 0 ; 
  }
  TList * l = phos->GetListOfTasks() ; 
  TIter it(l) ;
  TTask * task ;
  TTask * clu = 0 ;
  TString cluname(name) ;
  cluname+=":clu-" ;
  while((task = (TTask *)it.Next()) ){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(cluname)){
      clu = task ;
      break ;
    }
  }

  if(clu) 
    return (AliPHOSClusterizer **) l->GetObjectRef(clu) ;
  else
    return 0 ;
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostClusterizer(const char * name) const 
{ // ------------------ AliPHOSClusterizer ------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/recpointsname
  TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
  if ( !reF ) {
    cerr << "ERROR: AliPHOSGetter::Post Rer -> Task tasks/Reconstructioner" << " not found!" << endl;
    return kFALSE ; 
  } 
       
  TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
  
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Rer -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
    return kFALSE ; 
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
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  // the hierarchy is //YSALICE/WhiteBoard/TrackSegments/PHOS/tracksegmentsname  
  TFolder * tracksegmentsF  = (TFolder*)aliceF->FindObject("WhiteBoard/TrackSegments/PHOS") ; 
  tracksegmentsF->SetTitle("") ;  
  if ( !tracksegmentsF) {
    cerr << "ERROR: AliPHOSGetter::Post T -> Folder WhiteBoard/TrackSegments/PHOS" << " not found!" << endl;
    return kFALSE ; 
  }    
  TObject * tss =  tracksegmentsF->FindObject(name  ) ;
  if (!tss) {
    TClonesArray * ts = new TClonesArray("AliPHOSTrackSegment",100) ;
    ts->SetName(name) ;
    tracksegmentsF->Add(ts) ;  
  }
  return kTRUE; 
} 

//____________________________________________________________________________ 
TClonesArray ** AliPHOSGetter::TrackSegmentsRef(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  TFolder * phosF  = (TFolder*)gROOT->FindObjectAny("YSAlice/WhiteBoard/TrackSegments/PHOS") ; 
  if ( !phosF) {
    cerr << "ERROR: AliPHOSGetter::TrackSegmentsRef -> Folder WhiteBoard/TrackSegments/PHOS" << " not found!" << endl;
    return 0 ; 
  }    
  
  TObject * tss =  phosF->FindObject(name) ;
  if (!tss) {
    return 0 ;  
  }
  return (TClonesArray **) phosF->GetListOfFolders()->GetObjectRef(tss) ;
} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostTrackSegmentMaker(AliPHOSTrackSegmentMaker * tsmaker) const 
{ //------------Track Segment Maker ------------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/tracksegmentsname
  TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
  if ( !reF ) {
    cerr << "ERROR: AliPHOSGetter::Post Ter -> Task tasks/Reconstructioner" << " not found!" << endl;
    return kFALSE ; 
  }        
  TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Ter -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
    return kFALSE ; 
  }
  
  AliPHOSTrackSegmentMaker * phosts = 
    (AliPHOSTrackSegmentMaker*)phos->GetListOfTasks()->FindObject(tsmaker->GetName()) ; 
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
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/tracksegmentsname
  TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
  if ( !reF ) {
    cerr << "ERROR: AliPHOSGetter::Post Ter -> Task tasks/Reconstructioner" << " not found!" << endl;
    return kFALSE ; 
  }        
  TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Ter -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
    return kFALSE ; 
  } 
  AliPHOSTrackSegmentMaker * phosts = 
    (AliPHOSTrackSegmentMaker*)phos->GetListOfTasks()->FindObject(name) ; 
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
AliPHOSTrackSegmentMaker ** AliPHOSGetter::TSMakerRef(const char * name) const 
{ //------------Track Segment Maker ------------------------------
  
  TTask * reF  = (TTask*)gROOT->FindObjectAny("YSAlice/tasks/Reconstructioner") ; 
  if ( !reF ) {
    cerr << "ERROR: AliPHOSGetter::TrackSegmentMakerRef -> Task tasks/Reconstructioner" 
	 << " not found!" << endl;
    return 0 ; 
  }        
  
  TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ;   
  if ( !phos ) {
    cerr << "ERROR: AliPHOSGetter::TrackSegmentMakerRef -> Task tasks/Reconstructioner/PHOS" 
	 << " not found!" << endl;
    return 0 ; 
  }        
  TList * l = phos->GetListOfTasks() ; 
  TIter it(l) ;
  TTask * task ;
  TTask * tsm = 0 ;
  TString tsmname(name) ;
  tsmname+=":tsm-" ;
  while((task = (TTask *)it.Next()) ){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(tsmname)){
      tsm = task ;
      break ;
    }
  }
  
  if(tsm) 
    return (AliPHOSTrackSegmentMaker **) l->GetObjectRef(tsm) ;
  else
    return 0 ;
  
} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostRecParticles(const char * name) const 
{  // -------------------- RecParticles ------------------------
  
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // the hierarchy is //YSALICE/WhiteBoard/RecParticles/PHOS/recparticlesname  
  TFolder * recparticlesF  = (TFolder*)aliceF->FindObject("WhiteBoard/RecParticles/PHOS") ; 
  recparticlesF->SetTitle("") ;  
  if ( !recparticlesF) {
    cerr << "ERROR: AliPHOSGetter::Post P -> Folder WhiteBoard/RecParticles/PHOS" << " not found!" << endl;
    return kFALSE ; 
  }    
  TObject * rps = recparticlesF->FindObject( name )  ;
  if ( !rps ) {
    TClonesArray * rp = new TClonesArray("AliPHOSRecParticle",100) ;
    rp->SetName(name) ;    
    recparticlesF->Add(rp) ;  
  }
  return kTRUE; 
} 
//____________________________________________________________________________ 
TClonesArray ** AliPHOSGetter::RecParticlesRef(const char * name) const 
{ // ---------------TrackSegments -----------------------------------
  
  TFolder * tsF  = (TFolder*)gROOT->FindObjectAny("YSAlice/WhiteBoard/RecParticles/PHOS") ; 
  if ( !tsF) {
    cerr << "ERROR: AliPHOSGetter::RecParticlesRef -> Folder WhiteBoard/RecParticles/PHOS" << " not found!" << endl;
    return 0 ; 
  }    
  TObject * tss =  tsF->FindObject(name  ) ;
  if (!tss) {
    return 0 ;  
  }
  return (TClonesArray **) tsF->GetListOfFolders()->GetObjectRef(tss) ;
}

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostPID(AliPHOSPID * pid) const 
{     
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // ------------AliPHOS PID -----------------------------
  // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/recparticlesname
  TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
  if ( !reF ) {
    cerr << "ERROR: AliPHOSGetter::Post Per -> Task tasks/Reconstructioner" << " not found!" << endl;
    return kFALSE ; 
  }        
  TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Per -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
    return kFALSE ; 
  } 
  AliPHOSPID * phospid = (AliPHOSPID*)phos->GetListOfTasks()->FindObject(pid->GetName()) ; 
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
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // ------------AliPHOS PID -----------------------------
  // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/recparticlesname
  TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
  if ( !reF ) {
    cerr << "ERROR: AliPHOSGetter::Post Per -> Task tasks/Reconstructioner" << " not found!" << endl;
    return kFALSE ; 
  }        
  TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
  if ( !phos )  {
    cerr <<"ERROR: AliPHOSGetter::Post Per -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
    return kFALSE ; 
  }

  TList * l = phos->GetListOfTasks() ;   
  TIter it(l) ;
  TString pidname(name) ;
  pidname+=":pid" ; 
  TTask * task ;
  while((task = (TTask *)it.Next()) ){
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
AliPHOSPID ** AliPHOSGetter::PIDRef(const char * name) const 
{ //------------PID ------------------------------

  TTask * reF  = (TTask*)gROOT->FindObjectAny("YSAlice/tasks/Reconstructioner") ; 
  if ( !reF ) {
    cerr << "ERROR: AliPHOSGetter::PIDRef -> Task tasks/Reconstructioner" 
	 << " not found!" << endl;
    return 0 ; 
  }        
  
  TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ;   
  
  if ( !phos ) {
    cerr << "ERROR: AliPHOSGetter::PIDRef -> Task tasks/Reconstructioner/PHOS" 
	 << " not found!" << endl;
    return 0 ; 
  }        
  TList * l = phos->GetListOfTasks() ; 
  TIter it(l) ;
  TTask * task ;
  TTask * pid = 0 ;
  TString pidname(name) ;
  pidname+=":pid-" ;
  while((task = (TTask *)it.Next()) ){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(pidname)){
      pid = task ;
      break ;
    }
  }
  
  if(pid) 
    return (AliPHOSPID **) l->GetObjectRef(pid) ;
  else
    return 0 ;
  
} 

//____________________________________________________________________________ 
Bool_t AliPHOSGetter::PostQA( const char * name) const 
{     
  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  
  // ------------------ QA ---------------------------------
  // the hierarchy is //YSALICE/WhiteBoard/Alarms/PHOS/
  
  TFolder * alarmsF  = new TFolder() ;
  TString alarmsName ; 
  if (name) 
    alarmsName = name ;
  else       
    alarmsName = "Alarm with no name" ;	
  alarmsF->SetName( alarmsName.Data() ) ; 
  alarmsF->SetTitle("") ;  
  TFolder * qaaF  = (TFolder*)aliceF->FindObject("WhiteBoard/QAAlarms") ; 
  if ( !qaaF) {
    cerr << "ERROR: AliPHOSGetter::Post QA -> Folder WhiteBoard/QAAlarms/" << " not found!" << endl;
    return kFALSE; 
  }    
  if ( qaaF->FindObject( alarmsName.Data() ) ) 
    qaaF->RecursiveRemove(  qaaF->FindObject( alarmsName.Data() ) ) ; 
  
  qaaF->Add(alarmsF) ;  
  
  return kTRUE;
}

//____________________________________________________________________________ 
const TParticle * AliPHOSGetter::Primary(Int_t index) const
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
void AliPHOSGetter::ReadTreeD()
{
  // Read the digit tree gAlice->TreeD()  
  if(gAlice->TreeD()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeD: can not read TreeD " << endl ;
  return ;
  }
  
  TObjArray * lob = (TObjArray*)gAlice->TreeD()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * digitsbranch = 0 ; 
  TBranch * digitizerbranch = 0 ; 
  Bool_t phosfound = kFALSE, digitizerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phosfound || !digitizerfound) ) {
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
  
  TBranch * hitsbranch = (TBranch*)gAlice->TreeH()->GetBranch("PHOS") ;
  if ( !hitsbranch ) {
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
  
  TBranch * hitsbranch = (TBranch*)gAlice->TreeH()->GetListOfBranches()->FindObject("PHOS") ;
  if ( !hitsbranch ) {
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
    cout << "WARNING: AliPHOSGetter::ReadTreeQA -> Cannot find QA Alarms for PHOS" << endl ;
    return ; 
  }   
  
  // Post the QA Alarms
  PostQA("PHOS") ; 
  TFolder * alarmsF = Alarms() ; 
  alarmsF->Clear() ; 
  qabranch->SetAddress(&alarmsF) ;
  qabranch->GetEntry(0) ;
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeR()
{
  // Read the reconstrunction tree gAlice->TreeR()

  if(gAlice->TreeR()== 0){
    cout <<   "ERROR: AliPHOSGetter::ReadTreeR: can not read TreeR " << endl ;
    return ;
  }
  
  // RecPoints 
  TObjArray * lob = (TObjArray*)gAlice->TreeR()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * emcbranch = 0 ; 
  TBranch * cpvbranch = 0 ; 
  TBranch * clusterizerbranch = 0 ; 
  Bool_t phosemcrpfound = kFALSE, phoscpvrpfound = kFALSE, clusterizerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phosemcrpfound || !phoscpvrpfound || !clusterizerfound) ) 
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
    cout << "WARNING: AliPHOSGetter::ReadTreeR -> Cannot find EmcRecPoints with title " 
	 << fRecPointsTitle << endl ;
    return ; 
  }   
  if ( !phoscpvrpfound ) {
    cout << "WARNING: AliPHOSGetter::ReadTreeR -> Cannot find CpvRecPoints with title " 
	 << fRecPointsTitle << endl ;
    return ; 
  }   
  if ( !clusterizerfound ) {
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
    
  while ( (branch = (TBranch*)next()) && (!phostsfound || !tsmakerfound) ) 
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
  
  while ( (branch = (TBranch*)next()) && (!phosrpafound || !pidfound) ) 
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
  TFolder * phosF = (TFolder *)gROOT->FindObjectAny("YSAlice/WhiteBoard/SDigits/PHOS") ;
  TCollection * folderslist = phosF->GetListOfFolders() ; 
  
  //Add current file to list if it is not there yet
  if ( (fHeaderFile != "aliroot") && ( !folderslist->Contains(fHeaderFile) ) ){
    phosF->AddFolder(fHeaderFile, ""); 
    folderslist = phosF->GetListOfFolders() ;
  }
  
  TIter next(folderslist) ; 
  TFolder * folder = 0 ; 
  TFile * file; 
  TTree * treeS = 0;
  while ( (folder = (TFolder*)next()) ) {
    if(fHeaderFile.CompareTo(folder->GetName()) == 0 ) 
      treeS=gAlice->TreeS() ;
    else{
      file = (TFile*)gROOT->GetFile(folder->GetName()); 
      file->cd() ;
      
      // Get SDigits Tree header from file
      TString treeName("TreeS") ;
      treeName += event ; 
      treeS = (TTree*)gDirectory->Get(treeName.Data());
    }
    if(treeS==0){
      cerr << "ERROR: AliPHOSGetter::ReadTreeS There is no SDigit Tree" << endl;
      return ;
    }
    
    //set address of the SDigits and SDigitizer
    TBranch   * sdigitsBranch    = 0;
    TBranch   * sdigitizerBranch = 0;
    TBranch   * branch           = 0 ;  
    TObjArray * lob = (TObjArray*)treeS->GetListOfBranches() ;
    TIter next(lob) ; 
    Bool_t phosfound = kFALSE, sdigitizerfound = kFALSE ; 
    
    while ( (branch = (TBranch*)next()) && (!phosfound || !sdigitizerfound) ) {
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
  folder = (TFolder*)next();
  if(folder){
    file   = (TFile*)gROOT->GetFile(folder->GetName()); 
    file   ->cd() ;
  }
  
}
//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeS(TTree * treeS, Int_t input)
{  // Read the summable digits fron treeS()  

  TString filename("mergefile") ;
  filename+= input ;
  TFolder * phosF =(TFolder*)gROOT->FindObjectAny("YSAlice/WhiteBoard/SDigits/PHOS")  ;   
  TFolder * folder=(TFolder*)phosF->FindObject(filename) ;
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
    cout << "WARNING: AliPHOSGetter::ReadTreeS -> Digits and/or Digitizer branch not found" << endl ;
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
void AliPHOSGetter::ReadPrimaries()
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
// //       cout << "AliPHOSGetter: cannot find Kine Tree for event:" << events->At(input) << endl;

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
void AliPHOSGetter::Event(Int_t event, const char* opt)
{
  // Reads the content of all Tree's S, D and R
  
  if ( event > gAlice->TreeE()->GetEntries() ) {
    cerr << "ERROR: AliPHOSGetter::Event -> There are only " 
	 << gAlice->TreeE()->GetEntries() << " events in this file" << endl ; 
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

  if( strstr(opt,"P") )
    ReadPrimaries() ;

}

//____________________________________________________________________________ 
const TObject * AliPHOSGetter::ReturnO(TString what, TString name, TString file) const 
{
  // get the object named "what" from the folder
  // folders are named like //YSAlice/WhiteBoard/what/PHOS/name

  if ( file.IsNull() ) 
    file = fHeaderFile ; 
  TString path("WhiteBoard/") ;
  if ( name.IsNull() ) {
    if ( what.CompareTo("Hits") == 0 ) {
      path += what ; 
      path += "/PHOS/"; 
      path += what ; 
    }
    else if ( what.CompareTo("SDigits") == 0 ) { 
      path += what ; 
      path += "/PHOS/"; 
      path += file ; 
      path += "/" ; 
      path += fSDigitsTitle ; 
    }
    else if ( what.CompareTo("Digits") == 0 ){
      path += what ; 
      path += "/PHOS/"; 
      path += fDigitsTitle ;
    } 
    else if ( what.CompareTo("EmcRecPoints") == 0 ) {
      path += "RecPoints/PHOS/";  
      path += "emc/" ; 
      path += fRecPointsTitle ; 
    }
    else if ( what.CompareTo("CpvRecPoints") == 0 ) {
      path += "RecPoints/PHOS/";  
      path += "cpv/" ; 
      path += fRecPointsTitle ; 
    }
    else if ( what.CompareTo("TrackSegments") == 0 ) {
      path += "TrackSegments/PHOS/";  
      path += fTrackSegmentsTitle ; 
    }  
    else if ( what.CompareTo("RecParticles") == 0 ) {
      path += "RecParticles/PHOS/";  
      path += fRecParticlesTitle ; 
    }  
     else if ( what.CompareTo("Alarms") == 0 ) {
      path += "QAAlarms/PHOS";   
    }  
  } 
  else {
    if ( what.CompareTo("SDigits") == 0 ) { 
      path += what ; 
      path += "/PHOS/"; 
      path += file ; 
      path += "/" ; 
      path += name ;
    } 
    else if ( what.CompareTo("Digits") == 0 ) {
      path += what ; 
      path += "/PHOS/"; 
      path += name ; 
    }
    else if ( what.CompareTo("EmcRecPoints") == 0 ) {
      path += "RecPoints/PHOS/";  
      path += "emc/" ; 
      path += name ; 
    }
    else if ( what.CompareTo("CpvRecPoints") == 0 ) {
      path += "RecPoints/PHOS/";  
      path += "cpv/" ; 
      path += name ; 
    }  
    else if ( what.CompareTo("TrackSegments") == 0 ) {
      path += "TrackSegments/PHOS/";  
      path += name ; 
    } 
    else if ( what.CompareTo("RecParticles") == 0 ) {
      path += "RecParticles/PHOS/";  
      path += name ; 
    } 
    else if ( what.CompareTo("Alarms") == 0 ) {
      path += "QAAlarms/PHOS/";  
      path += name ; 
    } 
  }
  path.Prepend("YSAlice/") ;
  TObject * phosO  = (TObject*)gROOT->FindObjectAny(path) ; 
  if (!phosO) {
    if(fDebug)
      cerr << "ERROR : AliPHOSGetter::ReturnO -> Object " << path << " not found!" << endl ; 
    return 0 ;
  }
  return phosO ;
}
  
//____________________________________________________________________________ 
const TTask * AliPHOSGetter::ReturnT(TString what, TString name) const 
{
  // get the TTask named "what" from the folder
  // folders are named like //YSAlice/Tasks/what/PHOS/name

  TString path("tasks") ;
 
  if ( what.CompareTo("SDigitizer") == 0 ) 
    path += "/SDigitizer" ;
  else if ( what.CompareTo("Digitizer") == 0 ) 
    path += "/Digitizer" ; 
  else if ( what.CompareTo("Clusterizer") == 0 ) 
    path += "/Reconstructioner" ; 
  else if ( what.CompareTo("TrackSegmentMaker") == 0 ) 
    path += "/Reconstructioner" ; 
  else if ( what.CompareTo("PID") == 0 ) 
    path += "/Reconstructioner" ; 
  else if ( what.CompareTo("QATasks") == 0 ) 
    path += "/QA" ; 

  TFolder * aliceF  = (TFolder*)gROOT ->FindObjectAny("YSAlice") ; 
  TTask   * aliceT  = (TTask*)  aliceF->FindObject(path) ; 

  if (!aliceT) {
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << path << " not found!" << endl ;  
    abort() ;
  }

  TTask   * phosT   = (TTask*)  aliceT->GetListOfTasks()->FindObject("PHOS") ; 
  if (!phosT) { 
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << path << "/PHOS not found!" << endl ;  
    abort() ;
  }
  TList * l = phosT->GetListOfTasks() ; 
 
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
  
  TIter it(l) ;
  TTask * task = 0 ; 
  while((task = (TTask *)it.Next()) ){
    TString taskname(task->GetName()) ;
    if(taskname.BeginsWith(name))
      return task ;
  }
  
  if(fDebug)
    cout << "WARNING: AliPHOSGetter::ReturnT -> Task " << path << "/" << name << " not found!" << endl ; 
  return 0 ;
}
