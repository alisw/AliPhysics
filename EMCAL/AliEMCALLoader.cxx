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

#include "TTree.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALLoader.h"
#include "AliEMCAL.h"
#include "AliEMCALHit.h"
#include "AliEMCALGetter.h"

ClassImp(AliEMCALLoader)
  
  
const TString AliEMCALLoader::fgkHitsName("HITS");//Name for TClonesArray with hits from one event
const TString AliEMCALLoader::fgkSDigitsName("SDIGITS");//Name for TClonesArray 
const TString AliEMCALLoader::fgkDigitsName("DIGITS");//Name for TClonesArray 
const TString AliEMCALLoader::fgkECARecPointsName("ECARECPOINTS");//Name for TClonesArray 
const TString AliEMCALLoader::fgkTracksName("TRACKS");//Name for TClonesArray 
const TString AliEMCALLoader::fgkRecParticlesName("RECPARTICLES");//Name for TClonesArray

const TString AliEMCALLoader::fgkECARecPointsBranchName("EMCALECARP");//Name for branch with ECA Reconstructed Points
const TString AliEMCALLoader::fgkTrackSegmentsBranchName("EMCALTS");//Name for branch with TrackSegments
const TString AliEMCALLoader::fgkRecParticlesBranchName("EMCALRP");//Name for branch with Reconstructed Particles

//____________________________________________________________________________ 
AliEMCALLoader::AliEMCALLoader()
{
  fDebug = 0;
  fRecParticlesLoaded = kFALSE;
}

//____________________________________________________________________________ 
AliEMCALLoader::AliEMCALLoader(const Char_t *detname,const Char_t *eventfoldername):
  AliLoader(detname,eventfoldername)
{
  fDebug=0;
  fRecParticlesLoaded = kFALSE;
}

//____________________________________________________________________________ 
AliEMCALLoader::~AliEMCALLoader()
{
  //remove and delete arrays
  Clean(fgkHitsName);
  Clean(fgkSDigitsName);
  Clean(fgkDigitsName);
  Clean(fgkECARecPointsName);
  Clean(fgkTracksName);
  Clean(fgkRecParticlesName);
  CleanFolders() ; 
 // set to 0x0 the objgetter in AliGetter ... weird isn it !
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  if (gime) 
    gime->Reset() ;
}

//____________________________________________________________________________ 
void AliEMCALLoader::CleanFolders()
{
  CleanRecParticles();
  AliLoader::CleanFolders();
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::SetEvent()
{
  //Cleans loaded stuff and and sets Files and Directories
  // do not post any data to folder/tasks

  
  Int_t retval = AliLoader::SetEvent();
  if (retval)
   {
     Error("SetEvent","AliLoader::SetEvent returned error");
     return retval;
   }


  if (Hits()) Hits()->Clear();
  if (SDigits()) SDigits()->Clear();
  if (Digits()) Digits()->Clear();
  if (ECARecPoints()) ECARecPoints()->Clear();
  if (TrackSegments()) TrackSegments()->Clear();
  if (RecParticles()) RecParticles()->Clear();
   
  return 0;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::GetEvent()
{
  //Overloads GetEvent method called by AliRunLoader::GetEvent(Int_t) method
  //to add Rec Particles specific for EMCAL
  
  //First call the original method to get whatever from std. setup is needed
  Int_t retval;
  
  retval = AliLoader::GetEvent();
  if (retval)
    {
      Error("GetEvent","AliLoader::GetEvent returned error");
      return retval;
    }
  
  if (GetHitsDataLoader()->GetBaseDataLoader()->IsLoaded()) ReadHits();
  if (GetSDigitsDataLoader()->GetBaseDataLoader()->IsLoaded()) ReadSDigits();
  if (GetDigitsDataLoader()->GetBaseDataLoader()->IsLoaded()) ReadDigits();
  if (GetRecPointsDataLoader()->GetBaseDataLoader()->IsLoaded()) ReadRecPoints();
  if (GetTracksDataLoader()->GetBaseDataLoader()->IsLoaded()) ReadTracks();
  if (GetRecParticlesDataLoader()->GetBaseDataLoader()->IsLoaded()) ReadRecParticles();


  //Now, check if RecPart were loaded  
  return 0;
}

//____________________________________________________________________________ 
const AliEMCAL * AliEMCALLoader::EMCAL() 
{
  // returns the EMCAL object 
  AliEMCAL * emcal = dynamic_cast<AliEMCAL*>(GetModulesFolder()->FindObject(fDetectorName));
  if ( emcal == 0x0) 
    if (fDebug)
      cout << "WARNING: AliEMCALLoader::EMCAL -> EMCAL module not found in Folders" << endl ; 
  return emcal ; 
}  

//____________________________________________________________________________ 
const AliEMCALGeometry * AliEMCALLoader::EMCALGeometry() 
{
  // Gets the EMCAL Geometry object 
  AliEMCALGeometry * rv = 0 ; 
  if (EMCAL() )
    rv =  EMCAL()->GetGeometry();
  return rv ; 
} 

//____________________________________________________________________________ 
Int_t AliEMCALLoader::LoadHits(Option_t* opt)
{  
  //------- Hits ----------------------
  //Overload (extends) LoadHits implemented in AliLoader
  //
  Int_t res;
  
  //First call the AliLoader's method to send the TreeH to folder
  res = AliLoader::LoadHits(opt);
  
  if (res)
   {//oops, error
     Error("LoadHits","AliLoader::LoadHits returned error");
     return res;
   }

  //read the data from tree in folder and send it to folder
  res = ReadHits();
  return 0;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::LoadSDigits(Option_t* opt)
{  
  // Loads the SDigits array in the folder structure
  Int_t res;
  //First call the AliLoader's method to send the TreeS to folder
  res = AliLoader::LoadSDigits(opt);
  if (res)
    {//oops, error
      Error("PostSDigits","AliLoader::LoadSDigits returned error");
      return res;
    }
  return ReadSDigits();
  
} 

//____________________________________________________________________________ 
Int_t AliEMCALLoader::LoadDigits(Option_t* opt)
{ 
  // Loads the Digits array in the folder structure
  
  Int_t res;
  //First call the AliLoader's method to send the TreeS to folder
  res = AliLoader::LoadDigits(opt);
  if (res)
    {//oops, error
      Error("LoadDigits","AliLoader::LoadDigits returned error");
      return res;
    }
  return ReadDigits();
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::LoadRecPoints(Option_t* opt) 
{ 
  // Loads the RecPoints array in the folder structure
  Int_t res;
  //First call the AliLoader's method to send the TreeR to folder
  res = AliLoader::LoadRecPoints(opt);
  if (res)
    {//oops, error
      Error("LoadRecPoints","AliLoader::LoadRecPoints returned error");
      return res;
    }
  
  TFolder * emcalFolder = GetDetectorDataFolder();
  if ( emcalFolder  == 0x0 ) 
    {
      Error("LoadRecPoints","Can not get detector data folder");
      return 1;
    }
  return ReadRecPoints();
}

//____________________________________________________________________________ 
Int_t  AliEMCALLoader::LoadTracks(Option_t* opt)
{
  //Loads Tracks: Open File, Reads Tree and posts, Read Data and Posts
  if (GetDebug()) 
    printf("LoadTracks: opt = %s",opt);
  if (fTracksLoaded)
    {
      Warning("LoadTracks","Tracks are already loaded");
      return 0;
    }
  Int_t res;
  //First call the AliLoader's method to send the TreeS to folder
  if (GetTracksDataLoader()->GetBaseLoader(0)->IsLoaded() == kFALSE) 
    {//tracks can be loaded by LoadRecPoints
      res = AliLoader::LoadTracks(opt);
      if (res)
	{//oops, error
	  Error("LoadTracks","AliLoader::LoadTracks returned error");
	  return res;
	}
    }
  res = ReadTracks();
  if (res)
    {
      Error("LoadTracks","Error occured while reading Tracks");
      return res;
    }
  
  fTracksLoaded = kTRUE;
  return 0;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::LoadRecParticles(Option_t* opt) 
{ 
  // Loads the RecParticles array in the folder structure
  Int_t res;
  //First call the AliLoader's method to send the TreeS to folder
  res = AliLoader::LoadRecParticles(opt);
  if (res)
    {//oops, error
      Error("LoadRecParticles","AliLoader::LoadRecParticles returned error");
      return res;
    }
  
  TFolder * emcalFolder = GetDetectorDataFolder();
  if ( emcalFolder  == 0x0 ) 
    {
      Error("PostDigits","Can not get detector data folder");
      return 1;
    }
  return ReadRecParticles();
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::PostHits()
{
  // Post Hits
  Int_t reval = AliLoader::PostHits();
  if (reval)
    {
     Error("PostHits","AliLoader::  returned error");
     return reval;
    }
  return ReadHits();
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::PostSDigits()
{
  // Posts the SDigits array to the folder structure
  Int_t reval = AliLoader::PostSDigits();
  if (reval)
   {
     Error("PostSDigits","AliLoader::PostSDigits  returned error");
     return reval;
   }
  return ReadSDigits();
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::PostDigits()
{
  // Posts the Digits array to the folder structure
  Int_t reval = AliLoader::PostDigits();
  if (reval)
    {
      Error("PostDigits","AliLoader::PostDigits  returned error");
      return reval;
    }
  return ReadDigits();
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::PostRecPoints()
{
  // Posts the RecPoints array to the folder structure
  Int_t reval = AliLoader::PostRecPoints();
  if (reval)
   {
     Error("PostRecPoints","AliLoader::PostRecPoints  returned error");
     return reval;
   }
  return ReadRecPoints();
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::PostRecParticles()
{
  // Posts the RecParticles array to the folder structure
  
  Int_t reval = AliLoader::PostRecParticles();
  if (reval)
    {
      Error("PostRecParticles","AliLoader::PostRecParticles  returned error");
      return reval;
    }
  return ReadRecParticles();
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::PostTracks()
{
  // Posts the Tracks array to the folder structure
  Int_t reval = AliLoader::PostTracks();
  if (reval)
    {
      Error("PostTracks","AliLoader::PostTracks  returned error");
      return reval;
    }
  return ReadTracks();
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::ReadHits()
{
  // If there is no Clones Array in folder creates it and sends to folder
  // then tries to read
  // Reads the first entry of EMCAL branch in hit tree TreeH()
  // Reads data from TreeH and stores it in TClonesArray that sits in DetectorDataFolder
  //
  TObject** hitref = HitsRef();
  if(hitref == 0x0)
    {
      MakeHitsArray();
     hitref = HitsRef();
    }
  
  TClonesArray* hits = dynamic_cast<TClonesArray*>(*hitref);
  
  TTree* treeh = TreeH();
  
  if(treeh == 0)
    {
      Error("ReadHits"," Cannot read TreeH from folder");
      return 1;
    }
  
  TBranch * hitsbranch = treeh->GetBranch(fDetectorName);
  if (hitsbranch == 0) 
    {
      Error("ReadHits"," Cannot find branch EMCAL"); 
      return 1;
    }
  
  if (GetDebug()) 
    printf("ReadHits: Reading Hits");
  
  if (hitsbranch->GetEntries() > 1)
    {
      TClonesArray * tempo =  new TClonesArray("AliEMCALHit",1000);
      
      hitsbranch->SetAddress(&tempo);
      Int_t index = 0 ; 
      Int_t i = 0 ;
      for (i = 0 ; i < hitsbranch->GetEntries(); i++) 
	{
	  hitsbranch->GetEntry(i) ;
	  Int_t j = 0 ;
	  for ( j = 0 ; j < tempo->GetEntries() ; j++) 
	    {
	      AliEMCALHit* hit = (AliEMCALHit*)tempo->At(j); 
	      new((*hits)[index]) AliEMCALHit( *hit ) ;
	      index++ ; 
	    }
	}
      tempo->Delete() ; 
      delete tempo;
    }
  else 
    {
      hitsbranch->SetAddress(hitref);
      hitsbranch->GetEntry(0) ;
    }
  
  return 0;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::ReadSDigits()
{
  // Read the summable digits tree TreeS():
  // Check if TClones is in folder
  // if not create and add to folder
  // connect to tree if available
  // Read the data
  
  TObject** sdref = SDigitsRef();
  if(sdref == 0x0)
    {
      MakeSDigitsArray();
      sdref = SDigitsRef();
    }
  
  TTree * treeS = TreeS();
  if(treeS==0)
    {
      //May happen if file is truncated or new in LoadSDigits
      //Error("ReadSDigits","There is no SDigit Tree");
      return 0;
    }
  
  TBranch * branch = treeS->GetBranch(fDetectorName);
  if (branch == 0) 
    {//easy, maybe just a new tree
      //Error("ReadSDigits"," Cannot find branch EMCAL"); 
      return 0;
    }
  
  branch->SetAddress(SDigitsRef());
  branch->GetEntry(0);
  return 0;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::ReadDigits()
{
  // Read the summable digits tree TreeS():
  // Check if TClones is in folder
  // if not create and add to folder
  // connect to tree if available
  // Read the data
  
  TObject** dref = DigitsRef();
  if(dref == 0x0)
    {//if there is not array in folder, create it and put it there
      MakeDigitsArray();
      dref = DigitsRef();
    }
  
  TTree * treeD = TreeD();
  if(treeD==0)
    {
      //May happen if file is truncated or new in LoadSDigits
      //Error("ReadDigits","There is no Digit Tree");
      return 0;
    }
  
  TBranch * branch = treeD->GetBranch(fDetectorName);
  if (branch == 0) 
    {//easy, maybe just a new tree
      //Error("ReadDigits"," Cannot find branch ",fDetectorName.Data()); 
      return 0;
    }
  
  branch->SetAddress(dref);//connect branch to buffer sitting in folder
  branch->GetEntry(0);//get first event 
  
  return 0;  
}

//____________________________________________________________________________ 
void AliEMCALLoader::UnloadRecParticles()
{
  // Unloads the RecParticles array fromthe folder structure
  fRecParticlesLoaded = kFALSE;
  CleanRecParticles();
  if (fTracksLoaded == kFALSE) UnloadTracks();
}

//____________________________________________________________________________ 
void AliEMCALLoader::UnloadTracks()
{
  // Unloads the Tracks array fromthe folder structure
  CleanTracks();//free the memory
  //in case RecPart are loaded we can not onload tree and close the file
  if (fRecParticlesLoaded == kFALSE) AliLoader::UnloadTracks();
  fTracksLoaded = kFALSE;//mark that nobody needs them
}

//____________________________________________________________________________ 
void AliEMCALLoader::Track(Int_t itrack)
{
  // Read the first entry of EMCAL branch in hit tree gAlice->TreeH()
  if(TreeH()== 0)
    {
      if (LoadHits())
	{
	  Error("Track","Can not load hits.");
	  return;
	} 
    }
  
  TBranch * hitsbranch = dynamic_cast<TBranch*>(TreeH()->GetListOfBranches()->FindObject("EMCAL")) ;
  if ( !hitsbranch ) {
    if (fDebug)
      cout << "WARNING:  AliEMCALLoader::ReadTreeH -> Cannot find branch EMCAL" << endl ; 
    return ;
  }  
  if(!Hits()) PostHits();
  
  hitsbranch->SetAddress(HitsRef());
  hitsbranch->GetEntry(itrack);
  
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::ReadRecPoints()
{
  //Creates and posts to folder an array container, 
  //connects branch in tree (if exists), and reads data to array
  
  MakeRecPointsArray();
   
  TObjArray * eca = 0x0 ;  

  TTree * treeR = TreeR();
  
  if(treeR==0)
    {
      //May happen if file is truncated or new in LoadSDigits
      return 0;
    }
  
  Int_t retval = 0;

  TBranch * ecabranch = treeR->GetBranch(fgkECARecPointsBranchName);
  if (ecabranch == 0x0)
    {
      Error("ReadRecPoints","Can not get branch with ECA Rec. Points named %s",fgkECARecPointsBranchName.Data());
      retval = 2;
    }
  else
   {
     ecabranch->SetAddress(&eca);
     ecabranch->GetEntry(0) ;
   }


  Int_t ii ; 

  Int_t maxeca = eca->GetEntries() ; 
  for ( ii= 0 ; ii < maxeca ; ii++ ) 
    ECARecPoints()->Add(eca->At(ii)) ;

  return retval;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::ReadTracks()
{
  //Creates and posts to folder an array container, 
  //connects branch in tree (if exists), and reads data to arry
  
  TObject** trkref = TracksRef();
  if ( trkref == 0x0 )   
    {//Create and post array
      MakeTrackSegmentsArray();
      trkref = TracksRef();
    }
  
  TTree * treeT = TreeT();
  if(treeT==0)
    {
      //May happen if file is truncated or new in LoadSDigits, or the file is in update mode, 
      //but tracking was not performed yet for a current event
      //Error("ReadTracks","There is no Tree with Tracks");
      return 0;
    }
  
  TBranch * branch = treeT->GetBranch(fgkTrackSegmentsBranchName);
  if (branch == 0) 
    {//easy, maybe just a new tree
      Error("ReadTracks"," Cannot find branch named %s",fgkTrackSegmentsBranchName.Data());
      return 0;
    }
  
  branch->SetAddress(trkref);//connect branch to buffer sitting in folder
  branch->GetEntry(0);//get first event 
  
  return 0;
}

//____________________________________________________________________________ 
Int_t AliEMCALLoader::ReadRecParticles()
{
  //Reads Reconstructed  Particles from file
  //Creates and posts to folder an array container, 
  //connects branch in tree (if exists), and reads data to arry
  
  TObject** recpartref = RecParticlesRef();
  
  if ( recpartref == 0x0 )   
    {//Create and post array
      MakeRecParticlesArray();
      recpartref = RecParticlesRef();
    }
  
  TTree * treeP = TreeP();
  if(treeP==0)
    {
      //May happen if file is truncated or new in LoadSDigits, 
      //or the file is in update mode, 
      //but tracking was not performed yet for a current event
      //     Error("ReadRecParticles","There is no Tree with Tracks and Reconstructed Particles");
      return 0;
    }
  
  TBranch * branch = treeP->GetBranch(fgkRecParticlesBranchName);
  if (branch == 0) 
    {//easy, maybe just a new tree
      Error("ReadRecParticles"," Cannot find branch %s",fgkRecParticlesBranchName.Data()); 
      return 0;
    }
  
  branch->SetAddress(recpartref);//connect branch to buffer sitting in folder
  branch->GetEntry(0);//get first event 
  
  return 0;
}

//____________________________________________________________________________ 
AliEMCALGeometry* AliEMCALLoader::GetEMCALGeometry()
{
  //returns EMCAL geometry from gAlice 
  //static Method used by some classes where it is not convienient to pass eventfoldername
  if (gAlice == 0x0)
    return 0x0;
  AliEMCAL* emcal=dynamic_cast<AliEMCAL*>(gAlice->GetDetector("EMCAL"));
  if (emcal == 0x0)
    return 0x0;
  return emcal->GetGeometry();
}

//____________________________________________________________________________ 
AliEMCALLoader* AliEMCALLoader::GetEMCALLoader(const  char* eventfoldername)
{
  // Get an instance of the EMCALLoader object
  AliRunLoader* rn  = AliRunLoader::GetRunLoader(eventfoldername);
  if (rn == 0x0)
    {
      cerr<<"Error: <AliEMCALLoader::GetEMCALLoader>: "
	  << "Can not find Run Loader in folder "<<eventfoldername<<endl;
      return 0x0;
    }
  return dynamic_cast<AliEMCALLoader*>(rn->GetLoader("EMCALLoader"));
}

//____________________________________________________________________________ 
Bool_t AliEMCALLoader::BranchExists(const TString& recName)
{ 
  // Check is branch exists
  if (fBranchTitle.IsNull()) return kFALSE;
  TString dataname, zername ;
  TTree* tree;
  if(recName == "SDigits") {
    tree = TreeS();
    dataname = GetDetectorName();
    zername = "AliEMCALSDigitizer" ;
  }
  else if(recName == "Digits"){
    tree = TreeD();
    dataname = GetDetectorName();
    zername = "AliEMCALDigitizer" ;
  }
  else if(recName == "ECARecPoints"){
    tree = TreeR();
    dataname = fgkECARecPointsBranchName;
    zername = "AliEMCALClusterizer" ;
  }
  else if(recName == "TrackSegments"){
    tree = TreeT();
    dataname = fgkTrackSegmentsBranchName;
    zername = "AliEMCALTrackSegmentMaker";
  }        
  else if(recName == "RecParticles"){
    tree = TreeP();
    dataname = fgkRecParticlesBranchName;
    zername = "AliEMCALPID";
  }
  else
    return kFALSE ;
  
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
      Warning("BranchExists","branch %s  with title  %s ",dataname.Data(),fBranchTitle.Data());
      return kTRUE ;
    }
    if ( branchName.BeginsWith(zername) &&  branchTitle.BeginsWith(titleName) ){
      Warning("BranchExists","branch AliEMCAL... with title  %s ",branch->GetTitle());
      return kTRUE ; 
    }
  }
  return kFALSE ;
  
}

//____________________________________________________________________________ 
void AliEMCALLoader::SetBranchTitle(const TString& btitle)
{
  // Gives a name to a branch in the folder structure
  if (btitle.CompareTo(fBranchTitle) == 0) return;
  fBranchTitle = btitle;
  ReloadAll();
}

//____________________________________________________________________________ 
void AliEMCALLoader::CleanHits()
{
  // Clean hits	
  AliLoader::CleanHits();
  //Clear an array 
  TClonesArray* hits = Hits();
  if (hits) 
    hits->Clear();
}

//____________________________________________________________________________ 
void AliEMCALLoader::CleanSDigits()
{
  // Cleans the SDigits array in the folder structure
  AliLoader::CleanSDigits();
  TClonesArray* sdigits = SDigits();
  if (sdigits) sdigits->Clear();
  
}
//____________________________________________________________________________ 

void AliEMCALLoader::CleanDigits()
{
  // Cleans the Digits array in the folder structure
  AliLoader::CleanDigits();
  TClonesArray* digits = Digits();
  if (digits) digits->Clear();
}

//____________________________________________________________________________ 
void AliEMCALLoader::CleanRecPoints()
{
  // Cleans the RecPoints array in the folder structure
  AliLoader::CleanRecPoints();
  TObjArray* recpoints = ECARecPoints();
  if (recpoints) recpoints->Clear();

}

//____________________________________________________________________________ 
void AliEMCALLoader::CleanTracks()
{
  // Cleans the Tracks array in the folder structure
  AliLoader::CleanTracks();//tree
  //and clear the array
  TClonesArray* tracks = TrackSegments();
  if (tracks) 
    tracks->Clear();  
}

//____________________________________________________________________________ 
void AliEMCALLoader::CleanRecParticles()
{
  // Cleans the RecParticles array in the folder structure
  TClonesArray *recpar = RecParticles();
  if (recpar) 
    recpar->Clear();
}

//____________________________________________________________________________ 
// void AliEMCALLoader::ReadCalibrationDB(const char * database,const char * filename)
// {

//   if(fcdb && (strcmp(database,fcdb->GetTitle())==0))
//     return ;

//   TFile * file = gROOT->GetFile(filename) ;
//   if(!file)
//     file = TFile::Open(filename);
//   if(!file){
//     Error ("ReadCalibrationDB", "Cannot open file %s", filename) ;
//     return ;
//   }
//   if(fcdb)
//     fcdb->Delete() ;
//   fcdb = dynamic_cast<AliEMCALCalibrationDB *>(file->Get("AliEMCALCalibrationDB")) ;
//   if(!fcdb)
//     Error ("ReadCalibrationDB", "No database %s in file %s", database, filename) ;
// }
//____________________________________________________________________________ 

// AliEMCALSDigitizer*  AliEMCALLoader::EMCALSDigitizer() 
// { 
// //return EMCAL SDigitizer
//  return  dynamic_cast<AliEMCALSDigitizer*>(SDigitizer()) ;
// }

//____________________________________________________________________________ 
void AliEMCALLoader::MakeHitsArray()
{
  // Create the array for Hits
  if (Hits()) return;
  TClonesArray* hits = new TClonesArray("AliEMCALHit",1000);
  hits->SetName(fgkHitsName);
  GetDetectorDataFolder()->Add(hits);
}

//____________________________________________________________________________ 
void AliEMCALLoader::MakeSDigitsArray()
{
  // Create the array for SDigits
  if ( SDigits()) return;
  TClonesArray* sdigits = new TClonesArray("AliEMCALDigit",1);
  sdigits->SetName(fgkSDigitsName);
  GetDetectorDataFolder()->Add(sdigits);
}

//____________________________________________________________________________ 
void AliEMCALLoader::MakeDigitsArray()
{
  // Create the array for Digits
  if ( Digits()) return;
  TClonesArray* digits = new TClonesArray("AliEMCALDigit",1);
  digits->SetName(fgkDigitsName);
  GetDetectorDataFolder()->Add(digits);
  
}

//____________________________________________________________________________ 
void AliEMCALLoader::MakeRecPointsArray()
{
  // Make recpoints array
  if ( ECARecPoints() == 0x0) {
    if (GetDebug()>9) 
      printf("MakeRecPointsArray: Making array for ECA");
    TObjArray* eca = new TObjArray(100) ;
    eca->SetName(fgkECARecPointsName) ;
    GetDetectorDataFolder()->Add(eca);
   }
}

//____________________________________________________________________________ 
void AliEMCALLoader::MakeTrackSegmentsArray()
{
  // Create the array for TrackSegments
  if ( TrackSegments()) 
    return;
  TClonesArray * ts = new TClonesArray("AliEMCALTrackSegment",100) ;
  ts->SetName(fgkTracksName);
  GetDetectorDataFolder()->Add(ts);
}

//____________________________________________________________________________ 
void AliEMCALLoader::MakeRecParticlesArray()
{  
  // Create the array for RecParticles
  if ( RecParticles()) return;
  TClonesArray * rp = new TClonesArray("AliEMCALRecParticle",100) ;
  rp->SetName(fgkRecParticlesName);
  GetDetectorDataFolder()->Add(rp);
}
