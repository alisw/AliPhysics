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

/// \class AliMUONSimData
///
/// Class containing MUON data: hits, digits, rawclusters, globaltrigger, localtrigger, etc ..
/// The classe makes the lik between the MUON data lists and the event trees from loaders
///
/// \author Gines Martinez, Subatech,  September 2003
///

#include "AliMUONSimData.h"
#include "AliMUONDataIterator.h"
#include "AliMUONConstants.h"
#include "AliMUONHit.h"
#include "AliMUONDigit.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONRawCluster.h"

// This is from rec, classes in base should not depend on rec !!!
//#include "AliMUONTrack.h"
//#include "AliMUONTriggerTrack.h"

#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliLog.h"

#include <TString.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <Riostream.h>
#include <TFile.h>

/// \cond CLASSIMP
ClassImp(AliMUONSimData)
/// \endcond
 
//_____________________________________________________________________________
AliMUONSimData::AliMUONSimData()
  : AliMUONData(),
    fHits(0x0),
    //fSDigits(0x0),
    fNhits(0)
    //fNSdigits(0x0)
{
/// Default constructor
}
//_____________________________________________________________________________
AliMUONSimData::AliMUONSimData(AliLoader * loader, const char* name, const char* title)
  : AliMUONData(loader, name, title),
    fHits(0x0),
    //fSDigits(0x0),
    fNhits(0)
    //fNSdigits(0x0)
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMUONSimData::AliMUONSimData(const char* galiceFile)
  : AliMUONData(galiceFile, "MUONFolderSim"),
    fHits(0x0),
    //fSDigits(0x0),
    fNhits(0)
    //fNSdigits(0x0)
{
/// Constructor for loading data from gAlice file
}

//_____________________________________________________________________________
AliMUONSimData::~AliMUONSimData()
{
/// Destructor for AliMUONSimData
  if (fHits) {
    fHits->Delete();
    delete fHits;
  }
/*  
  if (fSDigits) {
    fSDigits->Delete();
    delete fSDigits;
  }
*/  
}
//____________________________________________________________________________
void AliMUONSimData::AddHit(Int_t fIshunt, Int_t track, Int_t detElemId, 
			 Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
			 Float_t tof, Float_t momentum, Float_t theta, 
			 Float_t phi, Float_t length, Float_t destep,
			 Float_t Xref,Float_t Yref,Float_t Zref)
{
 /// Add new hit to the hit list

  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliMUONHit(fIshunt, track, detElemId, 
				  idpart, X, Y, Z, 
				  tof, momentum, theta, 
				  phi, length, destep,
				  Xref,Yref,Zref);
}
/*
//_____________________________________________________________________________
void AliMUONSimData::AddSDigit(Int_t id, const AliMUONDigit& Sdigit)
{
/// Add a MUON Sdigit to the list of SDigits of the detection plane id

  TClonesArray &lSdigits = * SDigits(id) ; 
  new(lSdigits[fNSdigits[id]++]) AliMUONDigit(Sdigit);
}

//____________________________________________________________________________
TClonesArray*  AliMUONSimData::SDigits(Int_t DetectionPlane) const
{
/// Getting List of SDigits

  if (fSDigits)
    return ( (TClonesArray*) fSDigits->At(DetectionPlane) );
  else
    return NULL;
}
*/
//____________________________________________________________________________
void AliMUONSimData::FillOwn(Option_t* option)
{
/// Method to fill the trees

  const char *cH   = strstr(option,"H");
  //const char *cS   = strstr(option,"S");   // SDigits branches in TreeS

  //const char *cRP  = strstr(option,"RP");  // Reconstructed Particle in TreeP
  
  // Filling TreeH
  if ( TreeH() && cH ) 
  {
    TreeH()->Fill();
  }  
/* 
  // Filling TreeS
  if ( TreeS() && cS) 
  {
    TreeS()->Fill();
  }
*/  
}

//_____________________________________________________________________________
void AliMUONSimData::MakeOwnBranch(Option_t* option)
{
/// Create Tree branches for the MUON.

  const Int_t kBufferSize = 4000;
  char branchname[30];
  
  //Setting Data Container
  SetDataContainer(option);  

  const char *cH   = strstr(option,"H");
  //const char *cS   = strstr(option,"S");   // Digits branches in TreeS
  
  TBranch * branch = 0x0;
  
  // Creating Branches for Hits
  if (TreeH() && cH) {
    sprintf(branchname,"%sHits",GetName());  
    branch = TreeH()->GetBranch(branchname);
    if (branch) {  
      AliInfo(Form("MakeBranch","Branch %s is already in tree.",branchname));
      return ;
    }
    branch = TreeH()->Branch(branchname,&fHits,kBufferSize);
    //Info("MakeBranch","Making Branch %s for hits \n",branchname);
  }  
/*  
  //Creating Branches for SDigits
  if (TreeS() && cS ) {
    // one branch for Sdigits per chamber
    for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
      sprintf(branchname,"%sSDigits%d",GetName(),iDetectionPlane+1);
      branch = 0x0;
      branch = TreeS()->GetBranch(branchname);
      if (branch) {  
        AliInfo(Form("Branch %s is already in tree.",branchname));
        return;
      }
      TClonesArray * sdigits = SDigits(iDetectionPlane); 
      branch = TreeS()->Branch(branchname, &sdigits, kBufferSize,1);
      //Info("MakeBranch","Making Branch %s for sdigits in detection plane %d\n",branchname,iDetectionPlane+1);
    }
  }
*/  
}

//____________________________________________________________________________
void AliMUONSimData::SetOwnDataContainer(Option_t* option)
{
/// Setting data containers of muon data

  const char *cH   = strstr(option,"H");
  //const char *cS   = strstr(option,"S");   // SDigits
                                           //const char *cRP  = strstr(option,"RP");  // Reconstructed Particles  
  AliDebug(1,Form("option=%s",option));
  //
  // Clones array for hits
  if ( cH ) {
    if (fHits == 0x0) {
      fHits     = new TClonesArray("AliMUONHit",1000);
    }
    ResetHits();
  }
/*  
  //
  // Container for Sdigits
  if (cS) {
    if (fSDigits == 0x0) { 
      AliDebug(1,"Creating fSDigits TObjArray");
      fSDigits = new TObjArray(AliMUONConstants::NCh());
      fNSdigits= new Int_t[AliMUONConstants::NCh()];
      for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
	TClonesArray* a = new TClonesArray("AliMUONDigit",10000);
	a->SetOwner();
	fSDigits->AddAt(a,i);
	AliDebug(1,Form("fSDigits[%d]=%p",i,a));
        fNSdigits[i]=0;
      }
    }
    else {
      AliDebug(1,Form("fSDigits already there = %p",fSDigits));
    }
    ResetSDigits();
  }  
*/
}

//____________________________________________________________________________
void AliMUONSimData::SetOwnTreeAddress(Option_t* option)
{
  // Setting Data containers
  SetOwnDataContainer(option);

/// Setting Addresses to the events trees

  const char *cH   = strstr(option,"H");
  //const char *cS   = strstr(option,"S");   // SDigits branches in TreeS
  
  // Set branch address for the Hits, Digits, RawClusters, GlobalTrigger and LocalTrigger Tree.
  char branchname[30];
  TBranch * branch = 0x0;
  
  AliDebug(1,Form("option=%s",option));
  //
  // Branch address for hit tree
  if (TreeH() && fHits && cH) {
    sprintf(branchname,"%sHits",GetName());  
    branch = TreeH()->GetBranch(branchname);
    if (branch) {
      //      Info("SetTreeAddress","(%s) Setting for Hits",GetName());
      branch->SetAddress(&fHits);
    }
    else { //can be invoked before branch creation
      //AliWarning(Form("(%s) Failed for Hits. Can not find branch in tree.",GetName()));
    }
  }
/*    
  //
  // Branch address for Sdigit tree
  if (TreeS() && fSDigits && cS) {
    AliDebug(1,"Setting branch addresses");
    for (int i=0; i<AliMUONConstants::NCh(); i++) {
      sprintf(branchname,"%sSDigits%d",GetName(),i+1);
      if (fSDigits) {
        AliDebug(1,Form("TreeS=%p for ich=%d branchname=%s",
                        TreeS(),i,branchname));
        branch = TreeS()->GetBranch(branchname);
        TClonesArray * sdigits = SDigits(i);
        if (branch) branch->SetAddress( &sdigits );
        else AliWarning(Form("(%s) Failed for SDigits Detection plane %d. Can not find branch in tree.",GetName(),i));
      }
    }
  }
*/  
}

//____________________________________________________________________________
void AliMUONSimData::Fill(Option_t* option)
{
/// Method to fill the trees

  AliMUONData::Fill(option);
  FillOwn(option);
}

//_____________________________________________________________________________
void AliMUONSimData::MakeBranch(Option_t* option)
{
/// Create Tree branches for the MUON.

  AliMUONData::MakeBranch(option);
  MakeOwnBranch(option);
}

//____________________________________________________________________________
void AliMUONSimData::SetDataContainer(Option_t* option)
{
/// Setting data containers of muon data

  AliMUONData::SetDataContainer(option);
  SetOwnDataContainer(option);
}

//____________________________________________________________________________
void AliMUONSimData::SetTreeAddress(Option_t* option)
{
  // Setting Data containers
  SetDataContainer(option);
  
  AliMUONData::SetTreeAddress(option);
  SetOwnTreeAddress(option);
}

//____________________________________________________________________________
Int_t          
AliMUONSimData::GetNtracks() const      
{
/// Get number of entries in hits three

  Int_t ntrk = 0;
  if (fLoader && fLoader->TreeH())
    ntrk = (Int_t) fLoader->TreeH()->GetEntries();
  return ntrk;
}
/*
//____________________________________________________________________________
void AliMUONSimData::ResetSDigits()
{
/// Reset number of Sdigits and the Sdigits array for this detector

    if (fSDigits == 0x0) return;
    for ( int i=0;i<AliMUONConstants::NCh();i++ ) {
      if ((*fSDigits)[i])    ((TClonesArray*)fSDigits->At(i))->Clear();
      if (fNSdigits)  fNSdigits[i]=0;
    }
}
*/
//______________________________________________________________________________
void AliMUONSimData::ResetHits()
{
/// Reset number of clusters and the cluster array for this detector

  fNhits   = 0;
  if (fHits) fHits->Clear();
}

//_____________________________________________________________________________
void 
AliMUONSimData::DumpKine(Int_t event2Check)
{
/// Dump kinematics

  fRunLoader->LoadKinematics("READ");

  Int_t nevents = fRunLoader->GetNumberOfEvents();
  for (Int_t ievent=0; ievent<nevents; ievent++) {  // Event loop
    if ( event2Check != 0 ) ievent=event2Check;

    // Getting event ievent
    fRunLoader->GetEvent(ievent); 

    // Stack of particle for this event
    AliStack* stack = fRunLoader->Stack();

    Int_t nparticles = (Int_t) fRunLoader->Stack()->GetNtrack();
    printf(">>> Event %d, Number of particles is %d \n", ievent, nparticles);

    for (Int_t iparticle=0; iparticle<nparticles; iparticle++) {
      stack->Particle(iparticle)->Print("");  
    }
    if (event2Check!=0) ievent=nevents;
  }
  fRunLoader->UnloadKinematics();
}


//_____________________________________________________________________________
void 
AliMUONSimData::DumpHits(Int_t event2Check, Option_t* opt)
{
/// Dump hits

  fLoader->LoadHits("READ");

  // Event loop
  Int_t nevents = fRunLoader->GetNumberOfEvents();
  for (Int_t ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);

    // Getting event ievent
    fRunLoader->GetEvent(ievent); 
    SetTreeAddress("H");

    // Track loop
    Int_t ntracks = (Int_t) GetNtracks();
    for (Int_t itrack=0; itrack<ntracks; itrack++) {
      //Getting List of Hits of Track itrack
      GetTrack(itrack);

      Int_t nhits = (Int_t) Hits()->GetEntriesFast();
      printf(">>> Track %d, Number of hits %d \n",itrack,nhits);
      for (Int_t ihit=0; ihit<nhits; ihit++) {
	AliMUONHit* mHit = static_cast<AliMUONHit*>(Hits()->At(ihit));
	mHit->Print(opt);
      }
      ResetHits();
    }
    if (event2Check!=0) ievent=nevents;
  }
  fLoader->UnloadHits();
}
/*
//_____________________________________________________________________________
void 
AliMUONSimData::DumpSDigits(Int_t event2Check, Option_t* opt)
{
/// Dump SDigits

  fLoader->LoadSDigits("READ");
  
  // Event loop
  Int_t nevents = fRunLoader->GetNumberOfEvents();
  for (Int_t ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);

    // Getting event ievent
    fRunLoader->GetEvent(ievent);
    SetTreeAddress("S");
    GetSDigits();

    // Loop on chambers
    Int_t nchambers = AliMUONConstants::NCh(); ;
    for (Int_t ichamber=0; ichamber<nchambers; ichamber++) {
      TClonesArray* digits = SDigits(ichamber);

      // Loop on Sdigits
      Int_t ndigits = (Int_t)digits->GetEntriesFast();
      for (Int_t idigit=0; idigit<ndigits; idigit++) {
        AliMUONDigit* mDigit = static_cast<AliMUONDigit*>(digits->At(idigit));
        mDigit->Print(opt);
      }
    }
    ResetSDigits();
    if (event2Check!=0) ievent=nevents;
  }
  fLoader->UnloadSDigits();
}
*/
