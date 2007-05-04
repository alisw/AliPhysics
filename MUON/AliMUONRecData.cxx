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

/// \class AliMUONRecData
///
/// Class containing MUON data: hits, digits, rawclusters, globaltrigger, localtrigger, etc ..
/// The classe makes the lik between the MUON data lists and the event trees from loaders
///
/// \author Gines Martinez, Subatech,  September 2003
///

#include "AliMUONRecData.h"
#include "AliMUONDataIterator.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONTriggerTrack.h"

#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliLog.h"

#include <TString.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <Riostream.h>
#include <TFile.h>

/// \cond CLASSIMP
ClassImp(AliMUONRecData)
/// \endcond
 
//_____________________________________________________________________________
AliMUONRecData::AliMUONRecData()
  : AliMUONData(),
    fRawClusters(0x0),
    fRecTracks(0x0),
    fRecTriggerTracks(0x0),
    fNrawclusters(0x0),
    fNrectracks(0),
    fNrectriggertracks(0)
{
/// Default constructor
}
//_____________________________________________________________________________
AliMUONRecData::AliMUONRecData(AliLoader * loader, const char* name, const char* title)
  : AliMUONData(loader, name, title),
    fRawClusters(0x0),
    fRecTracks(0x0),
    fRecTriggerTracks(0x0),
    fNrawclusters(0x0),
    fNrectracks(0),
    fNrectriggertracks(0)
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMUONRecData::AliMUONRecData(const char* galiceFile)
  : AliMUONData(galiceFile, "MUONFolder"),
    fRawClusters(0x0),
    fRecTracks(0x0),
    fRecTriggerTracks(0x0),
    fNrawclusters(0x0),
    fNrectracks(0),
    fNrectriggertracks(0)
{
/// Constructor for loading data from gAlice file
}

//_____________________________________________________________________________
AliMUONRecData::~AliMUONRecData()
{
/// Destructor for AliMUONRecData

  if (fRawClusters) {
    fRawClusters->Delete();
    delete fRawClusters;
  }
  if (fRecTracks){
    fRecTracks->Delete();
    delete fRecTracks;
  }
  if (fRecTriggerTracks){
    fRecTriggerTracks->Delete();
    delete fRecTriggerTracks;
  }
}
//_____________________________________________________________________________
void AliMUONRecData::AddRawCluster(Int_t id, const AliMUONRawCluster& c)
{
/// Add a MUON rawcluster to the list in the detection plane id

  TClonesArray &lrawcl = *((TClonesArray*) fRawClusters->At(id));
  new(lrawcl[fNrawclusters[id]++]) AliMUONRawCluster(c);
}
//_____________________________________________________________________________
void AliMUONRecData::AddRecTrack(const AliMUONTrack& track)
{
/// Add a MUON rectrack

  TClonesArray &lrectracks = *fRecTracks;
  new(lrectracks[fNrectracks++]) AliMUONTrack(track);
}
//_____________________________________________________________________________
void AliMUONRecData::AddRecTriggerTrack(const AliMUONTriggerTrack& triggertrack)
{
/// Add a MUON triggerrectrack

  TClonesArray &lrectriggertracks = *fRecTriggerTracks;  
  new(lrectriggertracks[fNrectriggertracks++]) AliMUONTriggerTrack(triggertrack);
}
//____________________________________________________________________________
Bool_t   AliMUONRecData::IsRawClusterBranchesInTree()
{
/// Checking if there are RawCluster Branches In TreeR

  if (TreeR()==0x0) {
    AliError("No treeR in memory");
    return kFALSE;
  }
  else {
     char branchname[30];
     sprintf(branchname,"%sRawClusters1",GetName());
     TBranch * branch = 0x0;
     branch = TreeR()->GetBranch(branchname);
     if (branch)  return kTRUE;
     else return kFALSE;    
  }
}
//____________________________________________________________________________
Bool_t   AliMUONRecData::IsTrackBranchesInTree()
{
/// Checking if there are Track Branches In TreeT
  if (TreeT()==0x0) {
    AliError("No treeT in memory");
    return kFALSE;
  }
  else {
     char branchname[30];
     sprintf(branchname,"%sTrack",GetName());
     TBranch * branch = 0x0;
     branch = TreeT()->GetBranch(branchname);
     if (branch)  return kTRUE;
     else return kFALSE;    
  }
}
//____________________________________________________________________________
Bool_t   AliMUONRecData::IsTriggerBranchesInTree()
{
/// Checking if there are Trigger Branches In TreeR
 if (TreeR()==0x0) {
    AliError("No treeR in memory");
    return kFALSE;
  }
  else {
     char branchname[30];
     sprintf(branchname,"%sLocalTrigger",GetName());
     TBranch * branch = 0x0;
     branch = TreeR()->GetBranch(branchname);
     if (branch)  return kTRUE;
     else return kFALSE;    
  }
}
//____________________________________________________________________________
Bool_t   AliMUONRecData::IsTriggerTrackBranchesInTree()
{
/// Checking if there are TriggerTrack Branches In TreeT
  if (TreeT()==0x0) {
    AliError("No treeT in memory");
    return kFALSE;
  }
  else {
     char branchname[30];
     sprintf(branchname,"%sTriggerTrack",GetName());
     TBranch * branch = 0x0;
     branch = TreeT()->GetBranch(branchname);
     if (branch)  return kTRUE;
     else return kFALSE;    
  }
}
//____________________________________________________________________________
void AliMUONRecData::FillOwn(Option_t* option)
{
/// Method to fill the trees

  const char *cRC  = strstr(option,"RC");  // RawCluster branches in TreeR
  const char *cRT  = strstr(option,"RT");  // Reconstructed Track in TreeT
  const char *cRL = strstr(option,"RL");   // Reconstructed Trigger Track in TreeT
  const char *cTC = strstr(option,"TC");   // global and local Trigger branches Copy in TreeR
  
  char branchname[30];
  TBranch * branch = 0x0;

  // Filling TreeR
  
  if ( TreeR() && cRC && cTC )
  {
    TreeR()->Fill();
  }
  else
  {  
    if ( TreeR()  && cRC ) 
    {
      if ( IsTriggerBranchesInTree() ) 
      {
      // Branch per branch filling
        for (int i=0; i<AliMUONConstants::NTrackingCh(); i++) 
        {
          sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
          branch = TreeR()->GetBranch(branchname);
          branch->Fill();
        }
      }
      else  
      {
        TreeR()->Fill();
      }
    }
    
    if ( TreeR()  && cTC) 
    {
      if (IsRawClusterBranchesInTree()) 
      {
        // Branch per branch filling
        sprintf(branchname,"%sLocalTrigger",GetName());
        branch = TreeR()->GetBranch(branchname); 
        branch->Fill();
	sprintf(branchname,"%sRegionalTrigger",GetName());
        branch = TreeR()->GetBranch(branchname); 
        branch->Fill();
        sprintf(branchname,"%sGlobalTrigger",GetName());
        branch = TreeR()->GetBranch(branchname);
        branch->Fill();
      }
      else
      {
        TreeR()->Fill();
      }
    }
  }

  // Filling TreeT
  
  if ( TreeT() && cRT && cRL )
  {
    TreeT()->Fill();
  }
  else
  {
    if ( TreeT() && cRT ) 
    {
      if (IsTriggerTrackBranchesInTree()) 
      {
        sprintf(branchname,"%sTrack",GetName());  
        branch = TreeT()->GetBranch(branchname);
        branch->Fill();
      }
      else 
      {
        TreeT()->Fill();
      }
    }

    if ( TreeT() && cRL ) 
    {
      if (IsTrackBranchesInTree()) 
      {
        sprintf(branchname,"%sTriggerTrack",GetName());  
        branch = TreeT()->GetBranch(branchname);
        branch->Fill();
      }    
      else 
      {
        TreeT()->Fill();
      }
    }
  }
}

//_____________________________________________________________________________
void AliMUONRecData::MakeOwnBranch(Option_t* option)
{
/// Create Tree branches for the MUON.

  const Int_t kBufferSize = 4000;
  char branchname[30];
  
  //Setting Data Container
  SetDataContainer(option);  

  const char *cRC  = strstr(option,"RC");  // RawCluster branches in TreeR
  const char *cRT  = strstr(option,"RT");  // Reconstructed Track in TreeT
  const char *cRL  = strstr(option,"RL");  // Reconstructed Trigger Track in TreeT
                                           //const char *cRP  = strstr(option,"RP");  // Reconstructed Particle in TreeP
  const char *cTC = strstr(option,"TC");   // global and local Trigger branches Copy in TreeR
  
  TBranch * branch = 0x0;
  
  if (TreeR() && cRC ) {
    //  one branch for raw clusters per tracking detection plane
    //        
    Int_t i; 
    for (i=0; i<AliMUONConstants::NTrackingCh() ;i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);	
      branch = 0x0;
      branch = TreeR()->GetBranch(branchname);
      if (branch) {  
        AliInfo(Form("Branch %s is already in tree.",branchname));
        return;
      }
      branch = TreeR()->Branch(branchname, &((*fRawClusters)[i]),kBufferSize);
      //Info("MakeBranch","Making Branch %s for rawcluster in detection plane %d\n",branchname,i+1);
    }
  }
  
  if (TreeR() && cTC ) {
    //
    // one branch for global trigger
    //
    sprintf(branchname,"%sGlobalTrigger",GetName());
    branch = 0x0;
    branch = TreeR()->GetBranch(branchname);
    if (branch) {  
      AliInfo(Form("Branch GlobalTrigger is already in treeR."));
      return ;
    }
    branch = TreeR()->Branch(branchname, &fGlobalTrigger, kBufferSize);
    //Info("MakeBranch", "Making Branch %s for Global Trigger\n",branchname);

  //
    // one branch for regional trigger
    //  
    sprintf(branchname,"%sRegionalTrigger",GetName());
    branch = 0x0;
    branch = TreeR()->GetBranch(branchname);
    if (branch) {  
      AliInfo(Form("Branch RegionalTrigger is already in treeR."));
      return;
    }
    branch = TreeR()->Branch(branchname, &fRegionalTrigger, kBufferSize);
     
    //
    // one branch for local trigger
    //  
    sprintf(branchname,"%sLocalTrigger",GetName());
    branch = 0x0;
    branch = TreeR()->GetBranch(branchname);
    if (branch) {  
      AliInfo(Form("Branch LocalTrigger is already in treeR."));
      return;
    }
    branch = TreeR()->Branch(branchname, &fLocalTrigger, kBufferSize);
    //Info("MakeBranch", "Making Branch %s for Global Trigger\n",branchname);  
  }
  
  if (TreeT() && cRT ) {
    sprintf(branchname,"%sTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) {  
      AliInfo(Form("Branch %s is already in tree.",GetName()));
      return ;
    }
    branch = TreeT()->Branch(branchname,&fRecTracks,kBufferSize);
    //Info("MakeBranch","Making Branch %s for tracks \n",branchname);
  }  
  // trigger tracks
  if (TreeT() && cRL ) {
    sprintf(branchname,"%sTriggerTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) {  
      AliInfo(Form("Branch %s is already in tree.",GetName()));
      return ;
    }
    branch = TreeT()->Branch(branchname,&fRecTriggerTracks,kBufferSize);
    //Info("MakeBranch","Making Branch %s for trigger tracks \n",branchname);
  }  
}
//____________________________________________________________________________
void AliMUONRecData::SetOwnDataContainer(Option_t* option)
{
/// Setting data containers of muon data

  const char *cRC  = strstr(option,"RC");  // RawCluster
  const char *cRT  = strstr(option,"RT");  // Reconstructed Tracks
  const char *cRL  = strstr(option,"RL");  // Reconstructed Trigger Tracks
                                           //const char *cRP  = strstr(option,"RP");  // Reconstructed Particles  
  const char *cTC = strstr(option,"TC");   // global and local Trigger branches Copy in TreeR

  AliDebug(1,Form("option=%s",option));

  //
  // Containers for rawclusters, globaltrigger and local trigger tree
  if (cRC ) {
    if (fRawClusters == 0x0) {
      fRawClusters = new TObjArray(AliMUONConstants::NTrackingCh());
      fNrawclusters= new Int_t[AliMUONConstants::NTrackingCh()];
      for (Int_t i=0; i<AliMUONConstants::NTrackingCh();i++) {
	TClonesArray* tca = new TClonesArray("AliMUONRawCluster",10000);
	tca->SetOwner();
        fRawClusters->AddAt(tca,i); 
        fNrawclusters[i]=0;
      }
    }
    // ResetRawClusters(); 
    // It breaks the correct functioning of the combined reconstruction (AZ)
    
  }
  if (cTC ) {
    if (fLocalTrigger == 0x0) {
      fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);
    }
   if (fRegionalTrigger == 0x0) {
      fRegionalTrigger  = new TClonesArray("AliMUONRegionalTrigger",16);
    }
    if (fGlobalTrigger== 0x0) {
      fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1); 
    }
    // ResetTrigger(); 
    // This is not necessary here since trigger info ins copied from digits info on flight to RecPoint output
  }

  //
  // Containers for rectracks and rectrigger tracks
  if ( cRT ) {
    if (fRecTracks == 0x0)  {
      fRecTracks  = new TClonesArray("AliMUONTrack",100);
    }
    ResetRecTracks();
  }
  if (cRL) {
    if (fRecTriggerTracks == 0x0 && cRL)  {
      fRecTriggerTracks  = new TClonesArray("AliMUONTriggerTrack",100);
    }
    ResetRecTriggerTracks();
  }  
}

//____________________________________________________________________________
void AliMUONRecData::SetOwnTreeAddress(Option_t* option)
{
  // Setting Data containers
  SetOwnDataContainer(option);

/// Setting Addresses to the events trees

  const char *cRC  = strstr(option,"RC");  // RawCluster branches in TreeR
  const char *cRT  = strstr(option,"RT");  // Reconstructed Track in TreeT
  const char *cRL  = strstr(option,"RL");  // Reconstructed Trigger Track in TreeT
                                           //const char *cRP  = strstr(option,"RP");  // Reconstructed Particle in TreeP
  const char *cTC = strstr(option,"TC");   // global and local Trigger branches Copy in TreeR
  
  // Set branch address for the Hits, Digits, RawClusters, GlobalTrigger and LocalTrigger Tree.
  char branchname[30];
  TBranch * branch = 0x0;
  
  AliDebug(1,Form("option=%s",option));

  //
  // Branch address for rawclusters, globaltrigger and local trigger tree
  if ( TreeR()  && fRawClusters && cRC && !strstr(cRC,"RCC")) {
    for (int i=0; i<AliMUONConstants::NTrackingCh(); i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
      if (fRawClusters) {
        branch = TreeR()->GetBranch(branchname);
        if (branch) branch->SetAddress( &((*fRawClusters)[i]) );
        else AliWarning(Form("(%s) Failed for RawClusters Detection plane %d. Can not find branch in tree.",GetName(),i));
      }
    }
  }
  if ( TreeR()  && fLocalTrigger && cTC) {
    sprintf(branchname,"%sLocalTrigger",GetName());
    branch = TreeR()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fLocalTrigger);
    else AliWarning(Form("(%s) Failed for LocalTrigger. Can not find branch in treeR.",GetName()));
  }
  if ( TreeR()  && fRegionalTrigger && cTC) {
    sprintf(branchname,"%sRegionalTrigger",GetName());
    branch = TreeR()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fRegionalTrigger);
    else AliWarning(Form("(%s) Failed for RegionalTrigger. Can not find branch in treeR.",GetName()));
  }
  if ( TreeR() && fGlobalTrigger && cTC) {
    sprintf(branchname,"%sGlobalTrigger",GetName());
    branch = TreeR()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fGlobalTrigger);
    else AliWarning(Form("(%s) Failed for GlobalTrigger. Can not find branch in treeR.",GetName()));
  }

  // Rec Trakcs
  if ( TreeT() && fRecTracks && cRT ) {
    sprintf(branchname,"%sTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fRecTracks);
    else AliWarning(Form("(%s) Failed for Tracks. Can not find branch in tree.",GetName()));
  }
  // Trigger tracks
  if ( TreeT() && fRecTriggerTracks && cRL ) {
    sprintf(branchname,"%sTriggerTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fRecTriggerTracks);
    else AliWarning(Form("(%s) Failed for Trigger Tracks. Can not find branch in tree.",GetName()));
  }
}

//____________________________________________________________________________
void AliMUONRecData::Fill(Option_t* option)
{
/// Method to fill the trees

  AliMUONData::Fill(option);
  FillOwn(option);
}  

//_____________________________________________________________________________
void AliMUONRecData::MakeBranch(Option_t* option)
{
/// Create Tree branches for the MUON.

  AliMUONData::MakeBranch(option);
  MakeOwnBranch(option);
}  

//____________________________________________________________________________
void AliMUONRecData::SetDataContainer(Option_t* option)
{
/// Setting data containers of muon data
  AliMUONData::SetDataContainer(option);
  SetOwnDataContainer(option);
}

//____________________________________________________________________________
void AliMUONRecData::SetTreeAddress(Option_t* option)
{
  AliMUONData::SetTreeAddress(option);
  SetOwnTreeAddress(option);
}  

//____________________________________________________________________________
TClonesArray*  AliMUONRecData::RawClusters(Int_t DetectionPlane)
{
/// Getting Raw Clusters

  if (fRawClusters) 
    return ( (TClonesArray*) fRawClusters->At(DetectionPlane) );
  else
    return NULL;
}

//_______________________________________________________________________________
void AliMUONRecData::ResetRawClusters()
{
/// Reset number of raw clusters and the raw clust array for this detector

  for ( int i=0;i<AliMUONConstants::NTrackingCh();i++ ) {
    if ((*fRawClusters)[i])    ((TClonesArray*)fRawClusters->At(i))->Clear();
    if (fNrawclusters)  fNrawclusters[i]=0;
  }
}

//____________________________________________________________________________
void AliMUONRecData::ResetRecTracks()
{
/// Reset tracks information

  fNrectracks = 0;
  if (fRecTracks) fRecTracks->Delete(); // necessary to delete in case of memory allocation
}

//____________________________________________________________________________
void AliMUONRecData::ResetRecTriggerTracks()
{
/// Reset tracks information

  fNrectriggertracks = 0;
  if (fRecTriggerTracks) fRecTriggerTracks->Delete(); // necessary to delete in case of memory allocation
}

//_____________________________________________________________________________
void 
AliMUONRecData::DumpRecPoints(Int_t event2Check, Option_t* opt) 
{
/// Dump rec points

  fLoader->LoadRecPoints("READ");

  // Event loop
  Int_t nevents = fRunLoader->GetNumberOfEvents();
  for (Int_t ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);

    // Getting event ievent
    fRunLoader->GetEvent(ievent);
    Int_t nchambers = AliMUONConstants::NTrackingCh();
    SetTreeAddress("RC,TC"); 
    GetRawClusters();

    // Loop on chambers
    for (Int_t ichamber=0; ichamber<nchambers; ichamber++) {
      char branchname[30];    
      sprintf(branchname,"MUONRawClusters%d",ichamber+1);
      //printf(">>>  branchname %s\n",branchname);

      // Loop on rec points
      Int_t nrecpoints = (Int_t) RawClusters(ichamber)->GetEntriesFast();
      // printf(">>> Chamber %2d, Number of recpoints = %6d \n",ichamber+1, nrecpoints);
      for (Int_t irecpoint=0; irecpoint<nrecpoints; irecpoint++) {
	AliMUONRawCluster* mRecPoint = static_cast<AliMUONRawCluster*>(RawClusters(ichamber)->At(irecpoint));
	mRecPoint->Print(opt);
      }
    }
    ResetRawClusters();
    if (event2Check!=0) ievent=nevents;
  }
  fLoader->UnloadRecPoints();
}

//_____________________________________________________________________________
void 
AliMUONRecData::DumpTracks(Int_t event2Check, Option_t* opt) 
{
  /// Dump tracks
  
  fLoader->LoadTracks("READ");
  
  // Event loop
  Int_t nevents = fRunLoader->GetNumberOfEvents();
  
  for (Int_t ievent=0; ievent<nevents; ievent++) 
  {
    if (event2Check!=0) ievent=event2Check;
  
    printf(">>> Event %d \n",ievent);
    
    // Getting event ievent
    fRunLoader->GetEvent(ievent);
    
    Int_t nchambers = AliMUONConstants::NTrackingCh();
    SetTreeAddress("RT"); 
    GetRecTracks();
    
    // Loop on chambers
    for (Int_t ichamber=0; ichamber<nchambers; ++ichamber) 
    {
      // Loop on tracks
      TClonesArray* tracks = RecTracks();
      Int_t ntracks = (Int_t) tracks->GetEntriesFast();

      for (Int_t i=0; i<ntracks; ++i) 
      {
        AliMUONTrack* mTrack = static_cast<AliMUONTrack*>(tracks->At(i));
        mTrack->Print(opt);
      }
    }
    ResetRecTracks();
    if (event2Check!=0) ievent=nevents;
  }
  fLoader->UnloadTracks();
}

//_____________________________________________________________________________
void 
AliMUONRecData::DumpRecTrigger(Int_t event2Check, 
                            Int_t write, Bool_t readFromRP)
{
/// Reads and dumps trigger objects from MUON.RecPoints.root

  TClonesArray * globalTrigger;
  TClonesArray * localTrigger;
  
  // Do NOT print out all the info if the loop runs over all events 
  Int_t printout = (event2Check == 0 ) ? 0 : 1 ;  

  // Book a ntuple for more detailled studies
  TNtuple *tupleGlo = new TNtuple("TgtupleGlo","Global Trigger Ntuple","ev:global:slpt:shpt:uplpt:uphpt:lplpt:lplpt");
  TNtuple *tupleLoc = new TNtuple("TgtupleLoc","Local Trigger Ntuple","ev:LoCircuit:LoStripX:LoDev:StripY:LoLpt:LoHpt:y11:y21:x11");

  // counters
  Int_t sLowpt=0,sHighpt=0;
  Int_t uSLowpt=0,uSHighpt=0;
  Int_t lSLowpt=0,lSHighpt=0;

  AliMUONTriggerCrateStore* crateManager = new AliMUONTriggerCrateStore();   
  crateManager->ReadFromFile();

  AliMUONGeometryTransformer* transformer = new AliMUONGeometryTransformer(kFALSE);
  transformer->ReadGeometryData("volpath.dat", "geometry.root");

  TClonesArray*  triggerCircuit = new TClonesArray("AliMUONTriggerCircuit", 234);

  for (Int_t i = 0; i < AliMUONConstants::NTriggerCircuit(); i++)  {
      AliMUONTriggerCircuit* c = new AliMUONTriggerCircuit();
      c->SetTransformer(transformer);
      c->Init(i,*crateManager);
      TClonesArray& circuit = *triggerCircuit;
      new(circuit[circuit.GetEntriesFast()])AliMUONTriggerCircuit(*c);
      delete c;
  }
  
  Char_t fileName[30];
  if (!readFromRP) {
      AliInfoStream() << " reading from digits \n";
      fLoader->LoadDigits("READ");
      sprintf(fileName,"TriggerCheckFromDigits.root");
  } else {
      AliInfoStream() << " reading from RecPoints \n";
      fLoader->LoadRecPoints("READ");
      sprintf(fileName,"TriggerCheckFromRP.root");
  }

  
  AliMUONGlobalTrigger *gloTrg(0x0);
  AliMUONLocalTrigger *locTrg(0x0);

  Int_t nevents = fRunLoader->GetNumberOfEvents();
  for (Int_t ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    if (ievent%100==0 || event2Check) 
      AliInfoStream() << "Processing event " << ievent << endl;
    fRunLoader->GetEvent(ievent);
    
    if (!readFromRP) {
	SetTreeAddress("D,GLT"); 
	GetTriggerD();
    } else {    
	SetTreeAddress("RC,TC"); 
	GetTrigger();
    }

    globalTrigger = GlobalTrigger();
    localTrigger = LocalTrigger();
    
    Int_t nglobals = (Int_t) globalTrigger->GetEntriesFast(); // should be 1
    Int_t nlocals  = (Int_t) localTrigger->GetEntriesFast(); // up to 234
    if (printout) printf("###################################################\n");
    if (printout) printf("event %d nglobal %d nlocal %d \n",ievent,nglobals,nlocals);

    for (Int_t iglobal=0; iglobal<nglobals; iglobal++) { // Global Trigger
      gloTrg = static_cast<AliMUONGlobalTrigger*>(globalTrigger->At(iglobal));
      
      sLowpt+=gloTrg->SingleLpt() ;
      sHighpt+=gloTrg->SingleHpt() ;
      uSLowpt+=gloTrg->PairUnlikeLpt(); 
      uSHighpt+=gloTrg->PairUnlikeHpt();
      lSLowpt+=gloTrg->PairLikeLpt(); 
      lSHighpt+=gloTrg->PairLikeHpt();
      
      if (printout) gloTrg->Print("full");

    } // end of loop on Global Trigger

    for (Int_t ilocal=0; ilocal<nlocals; ilocal++) { // Local Trigger
      locTrg = static_cast<AliMUONLocalTrigger*>(localTrigger->At(ilocal));

      Bool_t xTrig=kFALSE;
      Bool_t yTrig=kFALSE;

      if ( locTrg->LoSdev()==1 && locTrg->LoDev()==0 && 
	   locTrg->LoStripX()==0) xTrig=kFALSE; // no trigger in X
      else xTrig=kTRUE;                         // trigger in X
      if (locTrg->LoTrigY()==1 && 
	  locTrg->LoStripY()==15 ) yTrig = kFALSE; // no trigger in Y
      else yTrig = kTRUE;                          // trigger in Y

      if (xTrig && yTrig) { // make Trigger Track if trigger in X and Y
	  
	  if (printout) locTrg->Print("full");
	  
	  AliMUONTriggerCircuit* circuit = (AliMUONTriggerCircuit*)triggerCircuit->At(locTrg->LoCircuit()-1); 
	  
	  tupleLoc->Fill(ievent,locTrg->LoCircuit(),locTrg->LoStripX(),locTrg->LoDev(),locTrg->LoStripY(),locTrg->LoLpt(),locTrg->LoHpt(),circuit->GetY11Pos(locTrg->LoStripX()),circuit->GetY21Pos(locTrg->LoStripX()+locTrg->LoDev()+1),circuit->GetX11Pos(locTrg->LoStripY()));
      }
      
    } // end of loop on Local Trigger
    
    // fill ntuple
    tupleGlo->Fill(ievent,nglobals,gloTrg->SingleLpt(),gloTrg->SingleHpt(),gloTrg->PairUnlikeLpt(),gloTrg->PairUnlikeHpt(),gloTrg->PairLikeLpt(),gloTrg->PairLikeHpt());
    
    ResetTrigger();
    if (event2Check!=0) ievent=nevents;
  } // end loop on event  
  
  // Print out summary if loop ran over all event
  if (!event2Check){

    printf("\n");
    printf("=============================================\n");
    printf("================  SUMMARY  ==================\n");
    printf("\n");
    printf("Total number of events processed %d \n", (event2Check==0) ? nevents : 1);
    printf("\n");
    printf(" Global Trigger output       Low pt  High pt\n");
    printf(" number of Single           :\t");
    printf("%i\t%i\t",sLowpt,sHighpt);
    printf("\n");
    printf(" number of UnlikeSign pair  :\t"); 
    printf("%i\t%i\t",uSLowpt,uSHighpt);
    printf("\n");
    printf(" number of LikeSign pair    :\t");  
    printf("%i\t%i\t",lSLowpt,lSHighpt);
    printf("\n");
    printf("=============================================\n");
    fflush(stdout);
  }
  
  if (write){
      TFile *myFile = new TFile(fileName, "RECREATE");
      tupleGlo->Write();
      tupleLoc->Write();
      myFile->Close();
  }

  fLoader->UnloadRecPoints();

  delete crateManager;
  delete transformer;
  delete triggerCircuit;
  
}
