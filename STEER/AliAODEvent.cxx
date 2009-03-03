/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     AOD base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TROOT.h>
#include <TTree.h>
#include <TFolder.h>
#include <TFriendElement.h>
#include <TProcessID.h>
#include <TCollection.h>
#include "Riostream.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODTrack.h"

ClassImp(AliAODEvent)

// definition of std AOD member names
  const char* AliAODEvent::fAODListName[kAODListN] = {"header",
						      "tracks",
						      "vertices",
						      "v0s",
						      "cascades",
						      "tracklets",
						      "jets",
						      "emcalCells",
						      "phosCells",
						      "caloClusters",
						      "fmdClusters",
						      "pmdClusters"
						      
};
//______________________________________________________________________________
AliAODEvent::AliAODEvent() :
  AliVEvent(),
  fAODObjects(new TList()),
  fAODFolder(0),
  fConnected(kFALSE),
  fHeader(0),
  fTracks(0),
  fVertices(0),
  fV0s(0),
  fCascades(0),
  fTracklets(0),
  fJets(0),
  fEmcalCells(0),
  fPhosCells(0),
  fCaloClusters(0),
  fFmdClusters(0),
  fPmdClusters(0)
{
  // default constructor
}

//______________________________________________________________________________
AliAODEvent::AliAODEvent(const AliAODEvent& aod):
  AliVEvent(aod),
  fAODObjects(new TList()),
  fAODFolder(new TFolder()),
  fConnected(kFALSE),
  fHeader(new AliAODHeader(*aod.fHeader)),
  fTracks(new TClonesArray(*aod.fTracks)),
  fVertices(new TClonesArray(*aod.fVertices)),
  fV0s(new TClonesArray(*aod.fV0s)),
  fCascades(new TClonesArray(*aod.fCascades)),
  fTracklets(new AliAODTracklets(*aod.fTracklets)),
  fJets(new TClonesArray(*aod.fJets)),
  fEmcalCells(new AliAODCaloCells(*aod.fEmcalCells)),
  fPhosCells(new AliAODCaloCells(*aod.fPhosCells)),
  fCaloClusters(new TClonesArray(*aod.fCaloClusters)),
  fFmdClusters(new TClonesArray(*aod.fFmdClusters)),
  fPmdClusters(new TClonesArray(*aod.fPmdClusters))
{
  // Copy constructor
  AddObject(fHeader);
  AddObject(fTracks);
  AddObject(fVertices);
  AddObject(fV0s);
  AddObject(fCascades);
  AddObject(fTracklets);
  AddObject(fJets);
  AddObject(fEmcalCells);
  AddObject(fPhosCells);
  AddObject(fCaloClusters);
  AddObject(fFmdClusters);
  AddObject(fPmdClusters);
  fConnected = aod.fConnected;
  GetStdContent();
}

//______________________________________________________________________________
AliAODEvent & AliAODEvent::operator=(const AliAODEvent& aod) {

    // Assignment operator

    if(&aod == this) return *this;
    AliVEvent::operator=(aod);

    fAODObjects      = new TList();
    fAODFolder       = new TFolder();
    fConnected       = aod.fConnected;
    fHeader          = new AliAODHeader(*aod.fHeader);
    fTracks          = new TClonesArray(*aod.fTracks);
    fVertices        = new TClonesArray(*aod.fVertices);
    fV0s             = new TClonesArray(*aod.fV0s);
    fCascades        = new TClonesArray(*aod.fCascades);
    fTracklets       = new AliAODTracklets(*aod.fTracklets);
    fJets            = new TClonesArray(*aod.fJets);
    fEmcalCells      = new AliAODCaloCells(*aod.fEmcalCells);
    fPhosCells       = new AliAODCaloCells(*aod.fPhosCells);
    fCaloClusters    = new TClonesArray(*aod.fCaloClusters);
    fFmdClusters     = new TClonesArray(*aod.fFmdClusters);
    fPmdClusters     = new TClonesArray(*aod.fPmdClusters);
    
    fAODObjects = new TList();
    
    AddObject(fHeader);
    AddObject(fTracks);
    AddObject(fVertices);
    AddObject(fV0s);
    AddObject(fCascades);
    AddObject(fTracklets);
    AddObject(fJets);
    AddObject(fEmcalCells);
    AddObject(fPhosCells);
    AddObject(fCaloClusters);
    AddObject(fFmdClusters);
    AddObject(fPmdClusters);
    GetStdContent();
    return *this;
}


//______________________________________________________________________________
AliAODEvent::~AliAODEvent() 
{
// destructor
    if(fAODObjects&&!fConnected)
    {
	delete fAODObjects;
	fAODObjects = 0;
    }

    delete fAODFolder;
}

//______________________________________________________________________________
void AliAODEvent::AddObject(TObject* obj) 
{
  // Add an object to the list of objects.
  // Please be aware that in order to increase performance you should
  // refrain from using TObjArrays (if possible). Use TClonesArrays, instead.
  
  fAODObjects->AddLast(obj);
}

//______________________________________________________________________________
void AliAODEvent::RemoveObject(TObject* obj) 
{
  // Removes an object from the list of objects.
  
  fAODObjects->Remove(obj);
}

//______________________________________________________________________________
TObject *AliAODEvent::FindListObject(const char *objName)
{
  // Return the pointer to the object with the given name.

  return fAODObjects->FindObject(objName);
}

//______________________________________________________________________________
void AliAODEvent::CreateStdContent() 
{
  // create the standard AOD content and set pointers

  // create standard objects and add them to the TList of objects
  AddObject(new AliAODHeader());
  AddObject(new TClonesArray("AliAODTrack", 0));
  AddObject(new TClonesArray("AliAODVertex", 0));
  AddObject(new TClonesArray("AliAODv0", 0));
  AddObject(new TClonesArray("AliAODcascade", 0));
  AddObject(new AliAODTracklets());
  AddObject(new TClonesArray("AliAODJet", 0));
  AddObject(new AliAODCaloCells());
  AddObject(new AliAODCaloCells());
  AddObject(new TClonesArray("AliAODCaloCluster", 0));
  AddObject(new TClonesArray("AliAODFmdCluster", 0));
  AddObject(new TClonesArray("AliAODPmdCluster", 0));
  // set names
  SetStdNames();

  // read back pointers
  GetStdContent();
  CreateStdFolders();
  return;
}

void  AliAODEvent::MakeEntriesReferencable()
{
    // Make all entries referencable in a subsequent process
    //
    TIter next(fAODObjects);
    TObject* obj;
    while ((obj = next()))
    {
	if(obj->InheritsFrom("TCollection"))
	    {
		AssignIDtoCollection((TCollection*)obj);
	    }
    }
}

//______________________________________________________________________________
void AliAODEvent::SetStdNames()
{
  // introduce the standard naming

  if(fAODObjects->GetEntries()==kAODListN){
    for(int i = 0;i < fAODObjects->GetEntries();i++){
      TObject *fObj = fAODObjects->At(i);
      if(fObj->InheritsFrom("TNamed")){
	((TNamed*)fObj)->SetName(fAODListName[i]);
      }
      else if(fObj->InheritsFrom("TClonesArray")){
	((TClonesArray*)fObj)->SetName(fAODListName[i]);
      }
    }
  }
  else{
    printf("%s:%d SetStdNames() Wrong number of Std Entries \n",(char*)__FILE__,__LINE__);
  }
} 

void AliAODEvent::CreateStdFolders()
{
    // Create the standard folder structure
    fAODFolder = gROOT->GetRootFolder()->AddFolder("AOD", "AOD");
    if(fAODObjects->GetEntries()==kAODListN){
	for(int i = 0;i < fAODObjects->GetEntries();i++){
	    TObject *fObj = fAODObjects->At(i);
	    if(fObj->InheritsFrom("TClonesArray")){
		fAODFolder->AddFolder(fAODListName[i], fAODListName[i], (TCollection*) fObj);
	    } else {
		fAODFolder->AddFolder(fAODListName[i], fAODListName[i], 0);
	    }
	}
    }
    else{
	printf("%s:%d CreateStdFolders() Wrong number of Std Entries \n",(char*)__FILE__,__LINE__);
    }
} 

//______________________________________________________________________________
void AliAODEvent::GetStdContent()
{
  // set pointers for standard content

  fHeader        = (AliAODHeader*)fAODObjects->FindObject("header");
  fTracks        = (TClonesArray*)fAODObjects->FindObject("tracks");
  fVertices      = (TClonesArray*)fAODObjects->FindObject("vertices");
  fV0s           = (TClonesArray*)fAODObjects->FindObject("v0s");
  fCascades      = (TClonesArray*)fAODObjects->FindObject("cascades");
  fTracklets     = (AliAODTracklets*)fAODObjects->FindObject("tracklets");
  fJets          = (TClonesArray*)fAODObjects->FindObject("jets");
  fEmcalCells    = (AliAODCaloCells*)fAODObjects->FindObject("emcalCells");
  fPhosCells     = (AliAODCaloCells*)fAODObjects->FindObject("phosCells");
  fCaloClusters  = (TClonesArray*)fAODObjects->FindObject("caloClusters");
  fFmdClusters   = (TClonesArray*)fAODObjects->FindObject("fmdClusters");
  fPmdClusters   = (TClonesArray*)fAODObjects->FindObject("pmdClusters");
}

//______________________________________________________________________________
void AliAODEvent::ResetStd(Int_t trkArrSize, 
			   Int_t vtxArrSize, 
			   Int_t v0ArrSize,
			   Int_t cascadeArrSize,
			   Int_t jetSize, 
			   Int_t caloClusSize, 
			   Int_t fmdClusSize, 
			   Int_t pmdClusSize
			   )
{
  // deletes content of standard arrays and resets size 
  
  fTracks->Delete();
  if (trkArrSize > fTracks->GetSize()) 
    fTracks->Expand(trkArrSize);

  fVertices->Delete();
  if (vtxArrSize > fVertices->GetSize()) 
    fVertices->Expand(vtxArrSize);
  	 
  fV0s->Delete();
  if (v0ArrSize > fV0s->GetSize()) 
    fV0s->Expand(v0ArrSize);
  
  fCascades->Delete();
  if (cascadeArrSize > fCascades->GetSize()) 
    fCascades->Expand(cascadeArrSize);
  
  fJets->Delete();
  if (jetSize > fJets->GetSize())
    fJets->Expand(jetSize);

  fCaloClusters->Delete();
  if (caloClusSize > fCaloClusters->GetSize()) 
    fCaloClusters->Expand(caloClusSize);

  fFmdClusters->Delete();
  if (fmdClusSize > fFmdClusters->GetSize()) 
    fFmdClusters->Expand(fmdClusSize);

  fPmdClusters->Delete();
  if (pmdClusSize > fPmdClusters->GetSize()) 
    fPmdClusters->Expand(pmdClusSize);

  // Reset the tracklets
  fTracklets->DeleteContainer();
  fPhosCells->DeleteContainer();  
  fEmcalCells->DeleteContainer();

}

void AliAODEvent::ClearStd()
{
  // clears the standard arrays
  fHeader        ->RemoveQTheta();
  fTracks        ->Delete();
  fVertices      ->Delete();
  fV0s           ->Delete();
  fCascades      ->Delete();
  fTracklets     ->DeleteContainer();
  fJets          ->Delete();
  fEmcalCells    ->DeleteContainer();
  fPhosCells     ->DeleteContainer();
  fCaloClusters  ->Delete();
  fFmdClusters   ->Clear();
  fPmdClusters   ->Clear();
}

//_________________________________________________________________
Int_t AliAODEvent::GetPHOSClusters(TRefArray *clusters) const
{
  // fills the provided TRefArray with all found phos clusters
  
  clusters->Clear();
  
  AliAODCaloCluster *cl = 0;
  Bool_t first = kTRUE;
  for (Int_t i = 0; i < GetNCaloClusters() ; i++) {
    if ( (cl = GetCaloCluster(i)) ) {
      if (cl->IsPHOSCluster()){
	if(first) {
	  new (clusters) TRefArray(TProcessID::GetProcessWithUID(cl)); 
	  first=kFALSE;
	}
	clusters->Add(cl);
	//printf("IsPHOS cluster %d, E %2.3f Size: %d \n",i,cl->E(),clusters->GetEntriesFast());
      }
    }
  }
  return clusters->GetEntriesFast();
}

//_________________________________________________________________
Int_t AliAODEvent::GetEMCALClusters(TRefArray *clusters) const
{
  // fills the provided TRefArray with all found emcal clusters

  clusters->Clear();
  cout<<"AOD event 1: nclus "<<GetNCaloClusters()<<endl;
  AliAODCaloCluster *cl = 0;
  Bool_t first = kTRUE;
  for (Int_t i = 0; i < GetNCaloClusters(); i++) {
    if ( (cl = GetCaloCluster(i)) ) {
      if (cl->IsEMCALCluster()){
	if(first) {
	  new (clusters) TRefArray(TProcessID::GetProcessWithUID(cl)); 
	  first=kFALSE;
	}
	clusters->Add(cl);
	//printf("IsEMCal cluster %d, E %2.3f Size: %d \n",i,cl->E(),clusters->GetEntriesFast());
      }
    }
  }
  return clusters->GetEntriesFast();
}


//______________________________________________________________________________
Int_t AliAODEvent::GetMuonTracks(TRefArray *muonTracks) const
{
  // fills the provided TRefArray with all found muon tracks

  muonTracks->Clear();

  AliAODTrack *track = 0;
  for (Int_t iTrack = 0; iTrack < GetNTracks(); iTrack++) {
    if ((track = GetTrack(iTrack))->IsMuonTrack()) {
      muonTracks->Add(track);
    }
  }
  
  return muonTracks->GetEntriesFast();
}


void AliAODEvent::ReadFromTree(TTree *tree, Option_t* opt /*= ""*/)
{
  // Connects aod event to tree
  
  if(!tree){
    Printf("%s %d AliAODEvent::ReadFromTree() Zero Pointer to Tree \n",(char*)__FILE__,__LINE__);
    return;
  }
  // load the TTree
  if(!tree->GetTree())tree->LoadTree(0);

  // Try to find AliAODEvent
  AliAODEvent *aodEvent = 0;
  aodEvent = (AliAODEvent*)tree->GetTree()->GetUserInfo()->FindObject("AliAODEvent");
  if(aodEvent){
    // Check if already connected to tree
    TList* connectedList = (TList*) (tree->GetUserInfo()->FindObject("AODObjectsConnectedToTree"));
    if (connectedList && (strcmp(opt, "reconnect"))) {
	// If connected use the connected list of objects
	printf("Delete and reconnect \n");
	
	fAODObjects->Delete();
	fAODObjects = connectedList;
	GetStdContent(); 
	fConnected = kTRUE;
	return;
    } 
    // Connect to tree
    // prevent a memory leak when reading back the TList
    if (!(strcmp(opt, "reconnect"))) fAODObjects->Delete();
    delete fAODObjects;
    fAODObjects = 0;
    // create a new TList from the UserInfo TList... 
    // copy constructor does not work...
    fAODObjects = (TList*)(aodEvent->GetList()->Clone());
    fAODObjects->SetOwner(kFALSE);
    if(fAODObjects->GetEntries()<kAODListN){
      printf("%s %d AliAODEvent::ReadFromTree() TList contains less than the standard contents %d < %d \n",
	     (char*)__FILE__,__LINE__,fAODObjects->GetEntries(),kAODListN);
    }
    //
    // Let's find out whether we have friends
    TList* friendL = tree->GetTree()->GetListOfFriends();
    if (friendL) 
    {
	TIter next(friendL);
	TFriendElement* fe;
	while ((fe = (TFriendElement*)next())){
	    aodEvent = (AliAODEvent*)(fe->GetTree()->GetUserInfo()->FindObject("AliAODEvent"));
	    if (!aodEvent) {
		printf("No UserInfo on tree \n");
	    } else {

		TList* objL = (TList*)(aodEvent->GetList()->Clone());
		printf("Get list of object from tree %d !!\n", objL->GetEntries());
		TIter nextobject(objL);
		TObject* obj =  0;
		while((obj = nextobject()))
		{
		    printf("Adding object from friend %s !\n", obj->GetName());
		    fAODObjects->Add(obj);
		} // object "branch" loop
	    } // has userinfo  
	} // friend loop
    } // has friends	
	    

// set the branch addresses
    TIter next(fAODObjects);
    TNamed *el;
    while((el=(TNamed*)next())){
      TString bname(el->GetName());
      // check if branch exists under this Name
      TBranch *br = tree->GetTree()->GetBranch(bname.Data());
      if(br){
	tree->SetBranchAddress(bname.Data(),fAODObjects->GetObjectRef(el));
      } else {
	  br = tree->GetBranch(Form("%s.",bname.Data()));
	  if(br){
	      tree->SetBranchAddress(Form("%s.",bname.Data()),fAODObjects->GetObjectRef(el));
	  }
	  else{
	      printf("%s %d AliAODEvent::ReadFromTree() No Branch found with Name %s. \n",
		     (char*)__FILE__,__LINE__,bname.Data());
	  }	
      }
    }
    GetStdContent();
    // when reading back we are not owner of the list 
    // must not delete it
    fAODObjects->SetOwner(kFALSE);
    fAODObjects->SetName("AODObjectsConnectedToTree");
    // we are not owner of the list objects 
    // must not delete it
    tree->GetUserInfo()->Add(fAODObjects);
    fConnected = kTRUE;
  }// no aodEvent
  else {
    // we can't get the list from the user data, create standard content
    // and set it by hand
    CreateStdContent();
    TIter next(fAODObjects);
    TNamed *el;
    while((el=(TNamed*)next())){
      TString bname(el->GetName());    
      tree->SetBranchAddress(bname.Data(),fAODObjects->GetObjectRef(el));
    }
    GetStdContent();
    // when reading back we are not owner of the list 
    // must not delete it
    fAODObjects->SetOwner(kFALSE);
  }
}

//______________________________________________________________________________
void AliAODEvent::Print(Option_t *) const
{
  // Something meaningful should be implemented here.
  
  return;
}

void AliAODEvent::AssignIDtoCollection(TCollection* col)
{
    // Static method which assigns a ID to each object in a collection
    // In this way the objects are marked as referenced and written with 
    // an ID. This has the advantage that TRefs to this objects can be 
    // written by a subsequent process.
    TIter next(col);
    TObject* obj;
    while ((obj = next()))
	TProcessID::AssignID(obj);
}
