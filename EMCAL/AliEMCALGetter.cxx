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
//  gime->Event(event) ;    // reads new event from galice.root
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//*--         Completely redesigned by Dmitri Peressounko March 2001  
//
//*-- YS June 2001 : renamed the original AliEMCALIndexToObject and make
//*--         systematic usage of TFolders without changing the interface        
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>


// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCAL.h"
#include "AliEMCALGetter.h"
#include "AliEMCALLoader.h"
#include "AliHeader.h"  
#include "AliMC.h"
#include "AliRunLoader.h"
#include "AliStack.h"  

ClassImp(AliEMCALGetter)
  
AliEMCALGetter * AliEMCALGetter::fgObjGetter = 0 ; 
AliEMCALLoader * AliEMCALGetter::fgEmcalLoader = 0;
Int_t AliEMCALGetter::fgDebug = 0;
TString AliEMCALGetter::fVersion = "";

//  TFile * AliEMCALGetter::fgFile = 0 ; 

//____________________________________________________________________________ 
AliEMCALGetter::AliEMCALGetter(const char* headerFile, const char* version, Option_t * openingOption) 
{
  // ctor only called by Instance()

  // initialize data members
  SetDebug(0) ; 
  //fBTE = 0 ; 
  
  fLoadingStatus = "" ; 

  fgObjGetter=this;
  
  OpenFile(headerFile,version,openingOption);
}


//____________________________________________________________________________ 
AliEMCALGetter::~AliEMCALGetter()
{
  AliRunLoader * rl = AliRunLoader::GetRunLoader(fgEmcalLoader->GetTitle());
  delete rl;
  fgEmcalLoader = 0 ;
  fgObjGetter = 0; 
  fVersion = "";
}

//____________________________________________________________________________ 
void AliEMCALGetter::OpenFile(const char* headerFile, const char* version, Option_t * openingOption) {
  fVersion = version;
  AliRunLoader* rl = AliRunLoader::GetRunLoader(version) ; 
  if (!rl) {
    rl = AliRunLoader::Open(headerFile, version, openingOption);
    if (!rl) {
      Fatal("AliEMCALGetter", "Could not find the Run Loader for %s - %s",headerFile, version) ; 
      return ;
    } 
    if (rl->GetAliRun() == 0x0) {
      rl->LoadgAlice();
      gAlice = rl->GetAliRun(); // should be removed
    }
  } 
  fgEmcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetLoader("EMCALLoader"));
  if ( !fgEmcalLoader ) 
    Error("AliEMCALGetter", "Could not find EMCALLoader") ; 
  else 
    fgEmcalLoader->SetTitle(version);
}

//____________________________________________________________________________ 
void AliEMCALGetter::Reset()
{
  // resets things in case the getter is called consecutively with different files
  // the EMCAL Loader is already deleted by the Run Loader

}

//____________________________________________________________________________ 
AliEMCALClusterizer * AliEMCALGetter::Clusterizer()
{ 
  // return pointer to Clusterizer Tree
  AliEMCALClusterizer * rv ; 
  rv =  dynamic_cast<AliEMCALClusterizer *>(EmcalLoader()->Reconstructioner()) ;
  if (!rv) {
    Event(0, "R") ; 
    rv =  dynamic_cast<AliEMCALClusterizer*>(EmcalLoader()->Reconstructioner()) ;
  }
  return rv ; 
}


//____________________________________________________________________________ 
TClonesArray * AliEMCALGetter::Digits() const  
{
  // asks the Loader to return the Digits container 

  TClonesArray * rv = 0 ; 
  rv = EmcalLoader()->Digits() ; 

  if( !rv ) {
    EmcalLoader()->MakeDigitsArray() ; 
    rv = EmcalLoader()->Digits() ;
  }
  return rv ; 
}

//____________________________________________________________________________ 
AliEMCALDigitizer * AliEMCALGetter::Digitizer() 
{ 
  // return pointer to Digitizer Tree
  AliEMCALDigitizer * rv ; 
  rv =  dynamic_cast<AliEMCALDigitizer *>(EmcalLoader()->Digitizer()) ;
  if (!rv) {
    Event(0, "D") ; 
    rv =  dynamic_cast<AliEMCALDigitizer *>(EmcalLoader()->Digitizer()) ;
  }
  return rv ; 
}

//____________________________________________________________________________ 
TObjArray * AliEMCALGetter::ECARecPoints() const 
{
  // asks the Loader to return the EMC RecPoints container 

  TObjArray * rv = 0 ; 
  
  rv = EmcalLoader()->ECARecPoints() ; 
  if (!rv) {
    EmcalLoader()->MakeRecPointsArray() ;
    rv = EmcalLoader()->ECARecPoints() ; 
  }
  return rv ; 
}

//____________________________________________________________________________ 
TClonesArray * AliEMCALGetter::TrackSegments() const 
{
  // asks the Loader to return the TrackSegments container 

  TClonesArray * rv = 0 ; 
  
  rv = EmcalLoader()->TrackSegments() ; 
  if (!rv) {
    EmcalLoader()->MakeTrackSegmentsArray() ;
    rv = EmcalLoader()->TrackSegments() ; 
  }
  return rv ; 
}

//____________________________________________________________________________ 
AliEMCALTrackSegmentMaker * AliEMCALGetter::TrackSegmentMaker() 
{ 
  // return pointer to TrackSegmentMaker Tree
  AliEMCALTrackSegmentMaker * rv ; 
  rv =  dynamic_cast<AliEMCALTrackSegmentMaker *>(EmcalLoader()->TrackSegmentMaker()) ;
  if (!rv) {
    Event(0, "T") ; 
    rv =  dynamic_cast<AliEMCALTrackSegmentMaker *>(EmcalLoader()->TrackSegmentMaker()) ;
  }
  return rv ; 
}

//____________________________________________________________________________ 
TClonesArray * AliEMCALGetter::RecParticles() const  
{
  // asks the Loader to return the TrackSegments container 

  TClonesArray * rv = 0 ; 
  
  rv = EmcalLoader()->RecParticles() ; 
  if (!rv) {
    EmcalLoader()->MakeRecParticlesArray() ;
    rv = EmcalLoader()->RecParticles() ; 
  }
  return rv ; 
}
//____________________________________________________________________________ 
void AliEMCALGetter::Event(Int_t event, const char* opt) 
{
  // Reads the content of all Tree's S, D and R

  if ( event >= MaxEvent() ) {
    Error("Event", "%d not found in TreeE !", event) ; 
    return ; 
  }

  AliRunLoader * rl = AliRunLoader::GetRunLoader(EmcalLoader()->GetTitle());
  // checks if we are dealing with test-beam data
//   TBranch * btb = rl->TreeE()->GetBranch("AliEMCALBeamTestEvent") ;
//   if(btb){
//     if(!fBTE)
//       fBTE = new AliEMCALBeamTestEvent() ;
//     btb->SetAddress(&fBTE) ;
//     btb->GetEntry(event) ;
//   }
//   else{
//     if(fBTE){
//       delete fBTE ;
//       fBTE = 0 ;
//     }
//   }

  // Loads the type of object(s) requested
  
  rl->GetEvent(event) ;

  if( strstr(opt,"X") || (strcmp(opt,"")==0) )
    ReadPrimaries() ;

  if(strstr(opt,"H") )
    ReadTreeH();

  if(strstr(opt,"S") )
    ReadTreeS() ;

  if( strstr(opt,"D") )
    ReadTreeD() ;

  if( strstr(opt,"R") )
    ReadTreeR() ;

  if( strstr(opt,"T") )
    ReadTreeT() ;

  if( strstr(opt,"P") )
    ReadTreeP() ;
 
}


//____________________________________________________________________________ 
Int_t AliEMCALGetter::EventNumber() const
  {
  // return the current event number
  AliRunLoader * rl = AliRunLoader::GetRunLoader(EmcalLoader()->GetTitle());
  return static_cast<Int_t>(rl->GetEventNumber()) ;   
}

//____________________________________________________________________________ 
  TClonesArray * AliEMCALGetter::Hits() const   
{
  // asks the loader to return  the Hits container 
  
  TClonesArray * rv = 0 ; 
  
  rv = EmcalLoader()->Hits() ; 
  if ( !rv ) {
    EmcalLoader()->LoadHits("read"); 
    rv = EmcalLoader()->Hits() ; 
  }
  return rv ; 
}

//____________________________________________________________________________ 
AliEMCALGetter * AliEMCALGetter::Instance(const char* alirunFileName, const char* version, Option_t * openingOption) 
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed

  if(!fgObjGetter){ // first time the getter is called 
    fgObjGetter = new AliEMCALGetter(alirunFileName, version, openingOption) ;
  }
  else { // the getter has been called previously
    AliRunLoader * rl = AliRunLoader::GetRunLoader(fVersion);
    if (rl == 0)  fgObjGetter->OpenFile(alirunFileName, version, openingOption) ; 
    else if ( rl->GetFileName() == alirunFileName ) {// the alirunFile has the same name
      // check if the file is already open
      TFile * galiceFile = dynamic_cast<TFile *>(gROOT->FindObject(rl->GetFileName()) ) ; 
      
      if ( !galiceFile ) 
	fgObjGetter->OpenFile(alirunFileName, version, openingOption);
      
      else {  // the file is already open check the version name
	TString currentVersionName = rl->GetEventFolder()->GetName() ; 
	TString newVersionName(version) ; 
	if (currentVersionName == newVersionName) 
	  if(fgDebug)
	    ::Warning( "Instance", "Files with version %s already open", currentVersionName.Data() ) ;  
	else {    
	  fgEmcalLoader->SetTitle(version); fVersion = version;
	}
      }
    }
    else { 
      AliRunLoader * rl = AliRunLoader::GetRunLoader(fVersion);
      if ( strstr(version, AliConfig::GetDefaultEventFolderName()) ) // false in case of merging
	delete rl ; 
      fgObjGetter->OpenFile(alirunFileName, version, openingOption) ;      
    }
  }
  if (!fgObjGetter) 
    ::Error("AliEMCALGetter::Instance", "Failed to create the EMCAL Getter object") ;
  else 
    if (fgDebug)
      Print() ;

  return fgObjGetter ;
}

//____________________________________________________________________________ 
AliEMCALGetter *  AliEMCALGetter::Instance()
{
  // Returns the pointer of the unique instance already defined
  
  if(!fgObjGetter && fgDebug)
     ::Warning("AliEMCALGetter::Instance", "Getter not initialized") ;

   return fgObjGetter ;
           
}

//____________________________________________________________________________ 
Int_t AliEMCALGetter::MaxEvent() const 
{
  // returns the number of events in the run (from TE)

  AliRunLoader * rl = AliRunLoader::GetRunLoader(EmcalLoader()->GetTitle());
  return static_cast<Int_t>(rl->GetNumberOfEvents()) ; 
}

//____________________________________________________________________________ 
TParticle * AliEMCALGetter::Primary(Int_t index) const
{
  AliRunLoader * rl = AliRunLoader::GetRunLoader(EmcalLoader()->GetTitle());
  return rl->Stack()->Particle(index) ; 
} 

//____________________________________________________________________________ 
Int_t AliEMCALGetter::NPrimaries() const
{
  AliRunLoader * rl = AliRunLoader::GetRunLoader(EmcalLoader()->GetTitle());
  return (rl->GetHeader())->GetNtrack(); 
} 

//____________________________________________________________________________ 
AliEMCAL * AliEMCALGetter:: EMCAL() const  
{
  // returns the EMCAL object 
  AliEMCAL * emcal = dynamic_cast<AliEMCAL*>(EmcalLoader()->GetModulesFolder()->FindObject("EMCAL")) ;  
  if (!emcal) 
    if (fgDebug)
      Warning("EMCAL", "EMCAL module not found in module folders: %s", EmcalLoader()->GetModulesFolder()->GetName() ) ; 
  return emcal ; 
}  



//____________________________________________________________________________ 
AliEMCALPID * AliEMCALGetter::PID() 
{ 
  // return pointer to PID Tree
  AliEMCALPID * rv ; 
  rv =  dynamic_cast<AliEMCALPID *>(EmcalLoader()->PIDTask()) ;
  if (!rv) {
    Event(0, "P") ; 
    rv =  dynamic_cast<AliEMCALPID *>(EmcalLoader()->PIDTask()) ;
  }
  return rv ; 
}

//____________________________________________________________________________ 
AliEMCALGeometry * AliEMCALGetter::EMCALGeometry() const 
{
  // Returns EMCAL geometry

  AliEMCALGeometry * rv = 0 ; 
  if (EMCAL() )
    rv =  EMCAL()->GetGeometry() ;
  return rv ; 
} 

//____________________________________________________________________________ 
void  AliEMCALGetter::Print() 
{
  // Print usefull information about the getter
    
  AliRunLoader * rl = AliRunLoader::GetRunLoader(fgEmcalLoader->GetTitle());
  ::Info("Print", "gAlice file is %s -- version name is %s", (rl->GetFileName()).Data(), rl->GetEventFolder()->GetName() ) ; 
}

//____________________________________________________________________________ 
void AliEMCALGetter::ReadPrimaries()  
{
  // Read Primaries from Kinematics.root
  
  AliRunLoader * rl = AliRunLoader::GetRunLoader(EmcalLoader()->GetTitle());
  
  // gets kine tree from the root file (Kinematics.root)
  if ( ! rl->TreeK() )  // load treeK the first time
    rl->LoadKinematics() ;
  
  if (fgDebug) 
    Info("ReadTreeK", "Found %d particles in event # %d", NPrimaries(), EventNumber() ) ; 
}

//____________________________________________________________________________ 
Int_t AliEMCALGetter::ReadTreeD()
{
  // Read the Digits
  
   EmcalLoader()->CleanDigits() ; 
  // gets TreeD from the root file (EMCAL.Digits.root)
  //if ( !IsLoaded("D") ) {
    EmcalLoader()->LoadDigits("UPDATE") ;
    EmcalLoader()->LoadDigitizer("UPDATE") ;
    //  SetLoaded("D") ; 
    //} 
  return Digits()->GetEntries() ; 
}

//____________________________________________________________________________ 
Int_t AliEMCALGetter::ReadTreeH()
{
  // Read the Hits
  EmcalLoader()->CleanHits() ; 
  // gets TreeH from the root file (EMCAL.Hit.root)
  //if ( !IsLoaded("H") ) {
    EmcalLoader()->LoadHits("READ") ;
    SetLoaded("H") ; 
    //}  
  return Hits()->GetEntries() ; 
}

//____________________________________________________________________________ 
Int_t AliEMCALGetter::ReadTreeR()
{
  // Read the RecPoints

   EmcalLoader()->CleanRecPoints() ; 
  // gets TreeR from the root file (EMCAL.RecPoints.root)
  //if ( !IsLoaded("R") ) {
    EmcalLoader()->LoadRecPoints("UPDATE") ;
    EmcalLoader()->LoadClusterizer("UPDATE") ;
    //  SetLoaded("R") ; 
    //}

  return ECARecPoints()->GetEntries() ; 
}

//____________________________________________________________________________ 
Int_t AliEMCALGetter::ReadTreeT()
{
  // Read the TrackSegments
  
  EmcalLoader()->CleanTracks() ; 
  // gets TreeT from the root file (EMCAL.TrackSegments.root)
  //if ( !IsLoaded("T") ) {
    EmcalLoader()->LoadTracks("UPDATE") ;
    EmcalLoader()->LoadTrackSegmentMaker("UPDATE") ;
    //  SetLoaded("T") ; 
    //}

  return TrackSegments()->GetEntries() ; 
}
//____________________________________________________________________________ 
Int_t AliEMCALGetter::ReadTreeP()
{
  // Read the RecParticles
  
  EmcalLoader()->CleanRecParticles() ; 
  // gets TreeP from the root file (EMCAL.RecParticles.root)
  //  if ( !IsLoaded("P") ) {
    EmcalLoader()->LoadRecParticles("UPDATE") ;
    EmcalLoader()->LoadPID("UPDATE") ;
    //  SetLoaded("P") ; 
    // }

  return RecParticles()->GetEntries() ; 
}
//____________________________________________________________________________ 
Int_t AliEMCALGetter::ReadTreeS()
{
  // Read the SDigits
  
  //  EmcalLoader()->CleanSDigits() ; 
  // gets TreeS from the root file (EMCAL.SDigits.root)
  // if ( !IsLoaded("S") ) {
    EmcalLoader()->LoadSDigits("READ") ;
    EmcalLoader()->LoadSDigitizer("READ") ;
    //  SetLoaded("S") ; 
    //}

  return SDigits()->GetEntries() ; 
}

//____________________________________________________________________________ 
TClonesArray * AliEMCALGetter::SDigits() const  
{
  // asks the Loader to return the Digits container 

  TClonesArray * rv = 0 ; 
  
  rv = EmcalLoader()->SDigits() ; 
  if (!rv) {
    EmcalLoader()->MakeSDigitsArray() ;
    rv = EmcalLoader()->SDigits() ; 
  }
  return rv ; 
}

//____________________________________________________________________________ 
AliEMCALSDigitizer * AliEMCALGetter::SDigitizer() 
{ 
  // Return pointer to SDigitizer task
  AliEMCALSDigitizer * rv ; 
  rv =  dynamic_cast<AliEMCALSDigitizer *>(EmcalLoader()->SDigitizer()) ;
  if (!rv) {
    Event(0, "S") ; 
    rv =  dynamic_cast<AliEMCALSDigitizer *>(EmcalLoader()->SDigitizer()) ;
  }
  return rv ; 
}

//____________________________________________________________________________ 
TParticle * AliEMCALGetter::Secondary(const TParticle* p, Int_t index) const
{
  // Return first (index=1) or second (index=2) secondary particle of primary particle p 

  if(index <= 0) 
    return 0 ;
  if(index > 2)
    return 0 ;

  if(p) {
  Int_t daughterIndex = p->GetDaughter(index-1) ; 
  AliRunLoader * rl = AliRunLoader::GetRunLoader(EmcalLoader()->GetTitle());
  return  rl->GetAliRun()->GetMCApp()->Particle(daughterIndex) ; 
  }
  else
    return 0 ;
}

//____________________________________________________________________________ 
void AliEMCALGetter::Track(Int_t itrack) 
{
  // Read the first entry of EMCAL branch in hit tree gAlice->TreeH()
 
 AliRunLoader * rl = AliRunLoader::GetRunLoader(EmcalLoader()->GetTitle());

  if( !TreeH() ) // load treeH the first time
    rl->LoadHits() ;

  // first time create the container
  TClonesArray * hits = Hits() ; 
  if ( hits ) 
    hits->Clear() ; 

  TBranch * emcalbranch = dynamic_cast<TBranch*>(TreeH()->GetBranch("EMCAL")) ; 
  emcalbranch->SetAddress(&hits) ;
  emcalbranch->GetEntry(itrack) ;
}

//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeD() const 
{
  // return pointer to Digits Tree
  TTree * rv = 0 ; 
  rv = EmcalLoader()->TreeD() ; 
  if ( !rv ) {
    EmcalLoader()->MakeTree("D");
    rv = EmcalLoader()->TreeD() ;
  } 
  
  return rv ; 
}

//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeH() const 
{
  // return pointer to Hits Tree
  TTree * rv = 0 ; 
  rv = EmcalLoader()->TreeH() ; 
  if ( !rv ) {
    EmcalLoader()->MakeTree("H");
    rv = EmcalLoader()->TreeH() ;
  } 
  
  return rv ; 
}

//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeR() const 
{
  // return pointer to RecPoints Tree

  TTree * rv = 0 ; 
  rv = EmcalLoader()->TreeR() ; 
  if ( !rv ) {
    EmcalLoader()->MakeTree("R");
    rv = EmcalLoader()->TreeR() ;
  } 
  
  return rv ; 
}

//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeT() const 
{  
  // return pointer to TrackSegments Tree
  TTree * rv = 0 ; 
  rv = EmcalLoader()->TreeT() ; 
  if ( !rv ) {
    EmcalLoader()->MakeTree("T");
    rv = EmcalLoader()->TreeT() ;
  } 
  
  return rv ; 
}
//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeP() const 
{
  // return pointer to RecParticles Tree
  TTree * rv = 0 ; 
  rv = EmcalLoader()->TreeP() ; 
  if ( !rv ) {
    EmcalLoader()->MakeTree("P");
    rv = EmcalLoader()->TreeP() ;
  } 
  
  return rv ; 
}

//____________________________________________________________________________ 
TTree * AliEMCALGetter::TreeS() const 
{
  // return pointer to SDigits Tree
  TTree * rv = 0 ; 
  rv = EmcalLoader()->TreeS() ; 
  if ( !rv ) {
    EmcalLoader()->MakeTree("S");
    rv = EmcalLoader()->TreeS() ;
  } 
  
  return rv ; 
}

//____________________________________________________________________________ 
Bool_t AliEMCALGetter::VersionExists(TString & opt) const
{
  // checks if the version with the present name already exists in the same directory

  Bool_t rv = kFALSE ;
 
  AliRunLoader * rl = AliRunLoader::GetRunLoader(EmcalLoader()->GetTitle());
  TString version( rl->GetEventFolder()->GetName() ) ; 

  opt.ToLower() ; 
  
  if ( opt == "sdigits") {
    // add the version name to the root file name
    TString fileName( EmcalLoader()->GetSDigitsFileName() ) ; 
    if (version != AliConfig::GetDefaultEventFolderName()) // only if not the default folder name 
      fileName = fileName.ReplaceAll(".root", "") + "_" + version + ".root" ;
    if ( !(gSystem->AccessPathName(fileName)) ) { 
      Warning("VersionExists", "The file %s already exists", fileName.Data()) ;
      rv = kTRUE ; 
    }
    EmcalLoader()->SetSDigitsFileName(fileName) ;
  }

  if ( opt == "digits") {
    // add the version name to the root file name
    TString fileName( EmcalLoader()->GetDigitsFileName() ) ; 
    if (version != AliConfig::GetDefaultEventFolderName()) // only if not the default folder name 
      fileName = fileName.ReplaceAll(".root", "") + "_" + version + ".root" ;
    if ( !(gSystem->AccessPathName(fileName)) ) {
      Warning("VersionExists", "The file %s already exists", fileName.Data()) ;  
      rv = kTRUE ; 
    }
  }

  return rv ;

}

//____________________________________________________________________________ 
UShort_t AliEMCALGetter::EventPattern(void) const
{
  // Return the pattern (trigger bit register) of the beam-test event
//   if(fBTE)
//     return fBTE->GetPattern() ;
//   else
    return 0 ;
}
//____________________________________________________________________________ 
Float_t AliEMCALGetter::BeamEnergy(void) const
{
  // Return the beam energy of the beam-test event
//   if(fBTE)
//     return fBTE->GetBeamEnergy() ;
//   else
    return 0 ;
}
