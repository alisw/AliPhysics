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

//-----------------------------------------------------------------
//           Implementation of the ESD class
//   This is the class to deal with during the phisical analysis of data
//   This class is generated directly by the reconstruction methods
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "TList.h"
#include <TNamed.h>

#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDVZERO.h"
#include "AliESDHLTtrack.h"
#include "AliESDFMD.h"
#include "AliESD.h"


ClassImp(AliESDEvent)

//______________________________________________________________________________
AliESDEvent::AliESDEvent():
  fESDObjects(new TList()),
  fESDRun(0),
  fHeader(0),
  fESDZDC(0),
  fESDFMD(0),
  fESDVZERO(0),
  fESDTZERO(0),
  fSPDVertex(0),
  fPrimaryVertex(0),
  fSPDMult(0),
  fPHOSTrigger(0),
  fEMCALTrigger(0),
  fTracks(0),
  fMuonTracks(0),
  fPmdTracks(0),
  fTrdTracks(0),
  fV0s(0),  
  fCascades(0),
  fKinks(0),
  fCaloClusters(0),
  fErrorLogs(0),
  fESDOld(0),
  fEMCALClusters(0), 
  fFirstEMCALCluster(-1),
  fPHOSClusters(0), 
  fFirstPHOSCluster(-1)
{
}
//______________________________________________________________________________
AliESDEvent::AliESDEvent(const AliESDEvent& esd):
  TObject(esd),
  fESDObjects(new TList()),
  fESDRun(new AliESDRun(*esd.fESDRun)),
  fHeader(new AliESDHeader(*esd.fHeader)),
  fESDZDC(new AliESDZDC(*esd.fESDZDC)),
  fESDFMD(new AliESDFMD(*esd.fESDFMD)),
  fESDVZERO(new AliESDVZERO(*esd.fESDVZERO)),
  fESDTZERO(new AliESDTZERO(*esd.fESDTZERO)),
  fSPDVertex(new AliESDVertex(*esd.fSPDVertex)),
  fPrimaryVertex(new AliESDVertex(*esd.fPrimaryVertex)),
  fSPDMult(new AliMultiplicity(*esd.fSPDMult)),
  fPHOSTrigger(new AliESDCaloTrigger(*esd.fPHOSTrigger)),
  fEMCALTrigger(new AliESDCaloTrigger(*esd.fEMCALTrigger)),
  fTracks(new TClonesArray(*esd.fTracks)),
  fMuonTracks(new TClonesArray(*esd.fMuonTracks)),
  fPmdTracks(new TClonesArray(*esd.fPmdTracks)),
  fTrdTracks(new TClonesArray(*esd.fTrdTracks)),
  fV0s(new TClonesArray(*esd.fV0s)),  
  fCascades(new TClonesArray(*esd.fCascades)),
  fKinks(new TClonesArray(*esd.fKinks)),
  fCaloClusters(new TClonesArray(*esd.fCaloClusters)),
  fErrorLogs(new TClonesArray(*esd.fErrorLogs)),
  fESDOld(new AliESD(*esd.fESDOld)),
  fEMCALClusters(esd.fEMCALClusters), 
  fFirstEMCALCluster(esd.fFirstEMCALCluster),
  fPHOSClusters(esd.fPHOSClusters), 
  fFirstPHOSCluster(esd.fFirstPHOSCluster)

{
  // CKB init in the constructor list and only add here ...
  AddObject(fESDRun);
  AddObject(fHeader);
  AddObject(fESDZDC);
  AddObject(fESDFMD);
  AddObject(fESDVZERO);
  AddObject(fESDTZERO);
  AddObject(fSPDVertex);
  AddObject(fPrimaryVertex);
  AddObject(fSPDMult);
  AddObject(fPHOSTrigger);
  AddObject(fEMCALTrigger);
  AddObject(fTracks);
  AddObject(fMuonTracks);
  AddObject(fPmdTracks);
  AddObject(fTrdTracks);
  AddObject(fV0s);
  AddObject(fCascades);
  AddObject(fKinks);
  AddObject(fCaloClusters);
  AddObject(fErrorLogs);

  GetStdContent();

}

//______________________________________________________________________________
AliESDEvent & AliESDEvent::operator=(const AliESDEvent& source) {

  // Assignment operator

  if(&source == this) return *this;
  TObject::operator=(source);

  fESDRun = new AliESDRun(*source.fESDRun);
  fHeader = new AliESDHeader(*source.fHeader);
  fESDZDC = new AliESDZDC(*source.fESDZDC);
  fESDFMD = new AliESDFMD(*source.fESDFMD);
  fESDVZERO = new AliESDVZERO(*source.fESDVZERO);
  fESDTZERO = new AliESDTZERO(*source.fESDTZERO);
  fSPDVertex = new AliESDVertex(*source.fSPDVertex);
  fPrimaryVertex = new AliESDVertex(*source.fPrimaryVertex);
  fSPDMult = new AliMultiplicity(*source.fSPDMult);
  fPHOSTrigger = new AliESDCaloTrigger(*source.fPHOSTrigger);
  fEMCALTrigger = new AliESDCaloTrigger(*source.fEMCALTrigger);
  fTracks = new TClonesArray(*source.fTracks);
  fMuonTracks = new TClonesArray(*source.fMuonTracks);
  fPmdTracks = new TClonesArray(*source.fPmdTracks);
  fTrdTracks = new TClonesArray(*source.fTrdTracks);
  fV0s = new TClonesArray(*source.fV0s);
  fCascades = new TClonesArray(*source.fCascades);
  fKinks = new TClonesArray(*source.fKinks);
  fCaloClusters = new TClonesArray(*source.fCaloClusters);
  fErrorLogs = new TClonesArray(*source.fErrorLogs);
  fESDOld = new AliESD(*source.fESDOld);
  // CKB this way?? or 
  // or AddObject(  fESDZDC = new AliESDZDC(*source.fESDZDC));

  fESDObjects = new TList();
  AddObject(fESDRun);
  AddObject(fHeader);
  AddObject(fESDZDC);
  AddObject(fESDFMD);
  AddObject(fESDVZERO);
  AddObject(fESDTZERO);
  AddObject(fSPDVertex);
  AddObject(fPrimaryVertex);
  AddObject(fSPDMult);
  AddObject(fPHOSTrigger);
  AddObject(fEMCALTrigger);
  AddObject(fTracks);
  AddObject(fMuonTracks);
  AddObject(fPmdTracks);
  AddObject(fTrdTracks);
  AddObject(fV0s);
  AddObject(fCascades);
  AddObject(fKinks);
  AddObject(fCaloClusters);
  AddObject(fErrorLogs);


  fEMCALClusters = source.fEMCALClusters;
  fFirstEMCALCluster = source.fFirstEMCALCluster;
  fPHOSClusters = source.fPHOSClusters;
  fFirstPHOSCluster = source.fFirstPHOSCluster;



  return *this;

}


//______________________________________________________________________________
AliESDEvent::~AliESDEvent()
{
  //
  // Standard destructor
  //

  // everthing on the list gets deleted automatically
  delete fESDObjects;
  fESDObjects = 0;
}

//______________________________________________________________________________
void AliESDEvent::Reset()
{

  
  // Reset the standard contents
  ResetStdContent(); 

  if(fESDOld)fESDOld->Reset();
  // call reset for user supplied data?
}

void AliESDEvent::ResetStdContent()
{
  // Reset the standard contents
  if(fESDRun) fESDRun->Reset();
  if(fHeader) fHeader->Reset();
  if(fESDZDC) fESDZDC->Reset();
  if(fESDFMD) fESDFMD->Clear(); // why clear.... need consistend names
  // if(fESDVZERO) fESDVZERO->; // NOT IMPLEMENTED 
  //  if(fESDVZERO) new (fESDVZERO) AliESDVZERO();
  if(fESDTZERO) fESDTZERO->Reset(); 
  // CKB no clear/reset implemented
  if(fSPDVertex){
    new (fSPDVertex) AliESDVertex();
    fSPDVertex->SetName("SPDVertex");
  }
  if(fPrimaryVertex){
    new (fPrimaryVertex) AliESDVertex();
    fPrimaryVertex->SetName("PrimaryVertex");
  }
  if(fSPDMult)new (fSPDMult) AliMultiplicity();
  if(fPHOSTrigger)fPHOSTrigger->Reset(); 
  if(fEMCALTrigger)fEMCALTrigger->Reset(); 
  if(fTracks)fTracks->Clear();
  if(fMuonTracks)fMuonTracks->Clear();
  if(fPmdTracks)fPmdTracks->Clear();
  if(fTrdTracks)fTrdTracks->Clear();
  if(fV0s)fV0s->Clear();
  if(fCascades)fCascades->Clear();
  if(fKinks)fKinks->Clear();
  if(fCaloClusters)fCaloClusters->Clear();
  if(fErrorLogs) fErrorLogs->Clear();

  fEMCALClusters=0; 
  fFirstEMCALCluster=-1; 
  fPHOSClusters=0; 
  fFirstPHOSCluster=-1; 
}


Int_t AliESDEvent::AddV0(const AliESDv0 *v) {
  //
  // Add V0
  //
  TClonesArray &fv = *fV0s;
  Int_t idx=fV0s->GetEntriesFast();
  new(fv[idx]) AliESDv0(*v);
  return idx;
}  

//______________________________________________________________________________
void AliESDEvent::Print(Option_t *) const 
{
  //
  // Print header information of the event
  //
  printf("ESD run information\n");
  printf("Event # in file %d Bunch crossing # %d Orbit # %d Period # %d Run # %d Trigger %lld Magnetic field %f \n",
	 GetEventNumberInFile(),
	 GetBunchCrossNumber(),
	 GetOrbitNumber(),
	 GetPeriodNumber(),
	 GetRunNumber(),
	 GetTriggerMask(),
	 GetMagneticField() );
  printf("Vertex: (%.4f +- %.4f, %.4f +- %.4f, %.4f +- %.4f) cm\n",
	   fPrimaryVertex->GetXv(), fPrimaryVertex->GetXRes(),
	   fPrimaryVertex->GetYv(), fPrimaryVertex->GetYRes(),
	   fPrimaryVertex->GetZv(), fPrimaryVertex->GetZRes());
    printf("Mean vertex in RUN: X=%.4f Y=%.4f cm\n",
	   GetDiamondX(),GetDiamondY());
    printf("SPD Multiplicity. Number of tracklets %d \n",
           fSPDMult->GetNumberOfTracklets());
  printf("Number of tracks: \n");
  printf("                 charged   %d\n", GetNumberOfTracks());
  printf("                 muon      %d\n", GetNumberOfMuonTracks());
  printf("                 pmd       %d\n", GetNumberOfPmdTracks());
  printf("                 trd       %d\n", GetNumberOfTrdTracks());
  printf("                 v0        %d\n", GetNumberOfV0s());
  printf("                 cascades  %d\n", GetNumberOfCascades());
  printf("                 kinks     %d\n", GetNumberOfKinks());
  printf("                 CaloClusters %d\n", GetNumberOfCaloClusters());
  printf("                 phos      %d\n", GetNumberOfPHOSClusters());
  printf("                 emcal     %d\n", GetNumberOfEMCALClusters());
  printf("                 FMD       %s\n", (fESDFMD ? "yes" : "no"));
  printf("                 VZERO     %s\n", (fESDVZERO ? "yes" : "no"));
}

void AliESDEvent::SetESDfriend(const AliESDfriend *ev) {
  //
  // Attaches the complementary info to the ESD
  //
  if (!ev) return;

  // to be sure that we set the tracks also
  // in case of old esds 
  if(fESDOld)CopyFromOldESD();

  Int_t ntrk=ev->GetNumberOfTracks();
 
  for (Int_t i=0; i<ntrk; i++) {
    const AliESDfriendTrack *f=ev->GetTrack(i);
    GetTrack(i)->SetFriendTrack(f);
  }
}

Int_t  AliESDEvent::AddTrack(const AliESDtrack *t) {
    // Add track
    TClonesArray &ftr = *fTracks;
    AliESDtrack * track = new(ftr[fTracks->GetEntriesFast()])AliESDtrack(*t);
    track->SetID(fTracks->GetEntriesFast()-1);
    return  track->GetID();    
}

Int_t AliESDEvent::AddKink(const AliESDkink *c) {
    // Add kink
    TClonesArray &fk = *fKinks;
    AliESDkink * kink = new(fk[fKinks->GetEntriesFast()]) AliESDkink(*c);
    kink->SetID(fKinks->GetEntriesFast()); // CKB different from the other imps..
    return fKinks->GetEntriesFast()-1;
}

Int_t AliESDEvent::AddCaloCluster(const AliESDCaloCluster *c) {
    // Add calocluster
    TClonesArray &fc = *fCaloClusters;
    AliESDCaloCluster *clus = new(fc[fCaloClusters->GetEntriesFast()]) AliESDCaloCluster(*c);
    clus->SetID(fCaloClusters->GetEntriesFast()-1);
    return fCaloClusters->GetEntriesFast()-1;
  }


void AliESDEvent::SetFMDData(AliESDFMD * obj) { 
  // use already allocated space
  if(fESDFMD){
    new(fESDFMD) AliESDFMD(*obj); 
  }
}

void AliESDEvent::SetVZEROData(AliESDVZERO * obj){ 
  // use already allocated space
  if(fESDVZERO)
    new(fESDVZERO) AliESDVZERO(*obj);
}

void AliESDEvent::GetESDfriend(AliESDfriend *ev) const {
  //
  // Extracts the complementary info from the ESD
  //
  if (!ev) return;

  Int_t ntrk=GetNumberOfTracks();

  for (Int_t i=0; i<ntrk; i++) {
    const AliESDtrack *t=GetTrack(i);
    const AliESDfriendTrack *f=t->GetFriendTrack();
    ev->AddTrack(f);
  }
}


void AliESDEvent::AddObject(TObject* obj) 
{
  // Add an object to the list of object.
  // Please be aware that in order to increase performance you should
  // refrain from using TObjArrays (if possible). Use TClonesArrays, instead.
  fESDObjects->SetOwner(kTRUE);
  fESDObjects->AddLast(obj);
}


void AliESDEvent::GetStdContent() 
{
  // set pointers for standard content
  // eventually get by name?
  fESDRun = (AliESDRun*)fESDObjects->At(kESDRun);
  fHeader = (AliESDHeader*)fESDObjects->At(kHeader);
  fESDZDC = (AliESDZDC*)fESDObjects->At(kESDZDC);
  fESDFMD = (AliESDFMD*)fESDObjects->At(kESDFMD);
  fESDVZERO = (AliESDVZERO*)fESDObjects->At(kESDVZERO);
  fESDTZERO = (AliESDTZERO*)fESDObjects->At(kESDTZERO);
  fSPDVertex = (AliESDVertex*)fESDObjects->At(kSPDVertex);
  fPrimaryVertex = (AliESDVertex*)fESDObjects->At(kPrimaryVertex);
  fSPDMult =       (AliMultiplicity*)fESDObjects->At(kSPDMult);
  fPHOSTrigger = (AliESDCaloTrigger*)fESDObjects->At(kPHOSTrigger);
  fEMCALTrigger = (AliESDCaloTrigger*)fESDObjects->At(kEMCALTrigger);
  fTracks = (TClonesArray*)fESDObjects->At(kTracks);
  fMuonTracks = (TClonesArray*)fESDObjects->At(kMuonTracks);
  fPmdTracks = (TClonesArray*)fESDObjects->At(kPmdTracks);
  fTrdTracks = (TClonesArray*)fESDObjects->At(kTrdTracks);
  fV0s = (TClonesArray*)fESDObjects->At(kV0s);
  fCascades = (TClonesArray*)fESDObjects->At(kCascades);
  fKinks = (TClonesArray*)fESDObjects->At(kKinks);
  fCaloClusters = (TClonesArray*)fESDObjects->At(kCaloClusters);
  fErrorLogs = (TClonesArray*)fESDObjects->At(kErrorLogs);

}

void AliESDEvent::SetStdNames(){
  // Set the names of the standard contents
  fSPDVertex->SetName("SPDVertex");
  fPrimaryVertex->SetName("PrimaryVertex");
  fPHOSTrigger->SetName("PHOSTrigger");
  fEMCALTrigger->SetName("EMCALTrigger");
  fTracks->SetName("Tracks");
  fMuonTracks->SetName("MuonTracks");
  fPmdTracks->SetName("PmdTracks");
  fTrdTracks->SetName("TrdTracks");
  fV0s->SetName("V0s");
  fCascades->SetName("Cascades");
  fKinks->SetName("Kinks");
  fCaloClusters->SetName("CaloClusters");

} 

void AliESDEvent::CreateStdContent() 
{
  // create the standard AOD content and set pointers

  // create standard objects and add them to the TList of objects
  AddObject(new AliESDRun());
  AddObject(new AliESDHeader());
  AddObject(new AliESDZDC());
  AddObject(new AliESDFMD());
  AddObject(new AliESDVZERO());
  AddObject(new AliESDTZERO());
  AddObject(new AliESDVertex());
  AddObject(new AliESDVertex());
  AddObject(new AliMultiplicity());
  AddObject(new AliESDCaloTrigger());
  AddObject(new AliESDCaloTrigger());
  AddObject(new TClonesArray("AliESDtrack",0));
  AddObject(new TClonesArray("AliESDMuonTrack",0));
  AddObject(new TClonesArray("AliESDPmdTrack",0));
  AddObject(new TClonesArray("AliESDTrdTrack",0));
  AddObject(new TClonesArray("AliESDv0",0));
  AddObject(new TClonesArray("AliESDcascade",0));
  AddObject(new TClonesArray("AliESDkink",0));
  AddObject(new TClonesArray("AliESDCaloCluster",0));
  AddObject(new TClonesArray("AliRawDataErrorLog",0));

  // check the order of the indices against enum...

  // read back pointers
  GetStdContent();
  // set names
  SetStdNames();

}

TObject* AliESDEvent::FindListObject(const char *name){
  if(fESDObjects)return fESDObjects->FindObject(name);
  return 0;
} 

void AliESDEvent::ReadFromTree(TTree *tree){
  
 
  // if we find the "ESD" branch on the tree we do have the old structure
  if(tree->GetBranch("ESD")){
      char ** address = (char **)(tree->GetBranch("ESD")->GetAddress());
      if (!address) {
	  printf("%s %d AliESDEvent::ReadFromTree() Reading old Tree \n",(char*)__FILE__,__LINE__);
	  tree->SetBranchAddress("ESD",&fESDOld);
      } else {
	  printf("%s %d AliESDEvent::ReadFromTree() Reading old Tree \n",(char*)__FILE__,__LINE__);
	  printf("%s %d Branch already connected. Using existing branch address. \n",(char*)__FILE__,__LINE__);
	  fESDOld = (AliESD*) (*address);
      }
      
      
      CreateStdContent(); // create for copy
      // when reading back we are not owner of the list 
      // must not delete it
      fESDObjects->SetOwner(kFALSE);
      return;
  }

  fESDOld = 0;


  // Try to find AliESDEvent
  AliESDEvent *esdEvent = 0;
  esdEvent = (AliESDEvent*)tree->GetTree()->GetUserInfo()->FindObject("AliESDEvent");
  //esdEvent = (AliESDEvent*)tree->GetUserInfo()->FindObject("AliESDEvent");

  if(esdEvent){   
      // Check if already connected to tree
//      TList* connectedList = (TList*) (tree->GetTree()->GetUserInfo()->FindObject("ESDObjectsConnectedToTree"));
      TList* connectedList = (TList*) (tree->GetUserInfo()->FindObject("ESDObjectsConnectedToTree"));
      if (connectedList) {
	  // If connected use the connected list if objects
	  fESDObjects->Delete();
	  fESDObjects = connectedList;
	  GetStdContent();
	  return;
      }
      // Connect to tree
    if(fESDObjects->GetEntries()!=0){
      // this should not happen here put a warning?
    }
    // prevent a memory leak when reading back the TList
    delete fESDObjects;
    fESDObjects = 0;
    // create a new TList from the UserInfo TList... 
    // copy constructor does not work...
    fESDObjects = (TList*)(esdEvent->GetList()->Clone());
    if(fESDObjects->GetEntries()<kESDListN){
      printf("%s %d AliESDEvent::ReadFromTree() TList contains less than the standard contents %d < %d \n",
	     (char*)__FILE__,__LINE__,fESDObjects->GetEntries(),kESDListN);
    }
    // set the branch addresses
    TIter next(fESDObjects);
    TNamed *el;
    while((el=(TNamed*)next())){
      TString bname(el->GetName());
      
      if(bname.CompareTo("AliESDfriend")==0)
	{
	  // AliESDfriend does not have a name ...
	  tree->SetBranchAddress("ESDfriend.",fESDObjects->GetObjectRef(el));
	}
      else{
	tree->SetBranchAddress(bname.Data(),fESDObjects->GetObjectRef(el));
      }
    }
    GetStdContent();
    // when reading back we are not owner of the list 
    // must not delete it
    fESDObjects->SetOwner(kFALSE);
    // Add list to user info
    fESDObjects->SetName("ESDObjectsConnectedToTree");
    tree->GetUserInfo()->Add(fESDObjects);
  }// no esdEvent
  else {
    // we can't get the list from the user data, create standard content
    // and set it by hand (no ESDfriend at the moment
    CreateStdContent();
    TIter next(fESDObjects);
    TNamed *el;
    while((el=(TNamed*)next())){
      TString bname(el->GetName());    
      tree->SetBranchAddress(bname.Data(),fESDObjects->GetObjectRef(el));
    }
    GetStdContent();
    // when reading back we are not owner of the list 
    // must not delete it
    fESDObjects->SetOwner(kFALSE);
  }



}


void AliESDEvent::CopyFromOldESD()
{
  // Method which copies over everthing from the old esd structure to the 
  // new  

  if(fESDOld){
    ResetStdContent();
     // Run
    SetRunNumber(fESDOld->GetRunNumber());
    SetPeriodNumber(fESDOld->GetPeriodNumber());
    SetMagneticField(fESDRun->GetMagneticField());
  
    // leave out diamond ...
    // SetDiamond(const AliESDVertex *vertex) { fESDRun->SetDiamond(vertex);}

    // header
    SetTriggerMask(fESDOld->GetTriggerMask());
    SetOrbitNumber(fESDOld->GetOrbitNumber());
    SetTimeStamp(fESDOld->GetTimeStamp());
    SetEventType(fESDOld->GetEventType());
    SetEventNumberInFile(fESDOld->GetEventNumberInFile());
    SetBunchCrossNumber(fESDOld->GetBunchCrossNumber());
    SetTriggerCluster(fESDOld->GetTriggerCluster());

    // ZDC

    SetZDC(fESDOld->GetZDCN1Energy(),
	   fESDOld->GetZDCP1Energy(),
	   fESDOld->GetZDCEMEnergy(),
	   fESDOld->GetZDCN2Energy(),
	   fESDOld->GetZDCP2Energy(),
	   fESDOld->GetZDCParticipants());

    // FMD
    
    SetFMDData(fESDOld->GetFMDData());

    // T0

    SetT0zVertex(fESDOld->GetT0zVertex());
    SetT0(fESDOld->GetT0());
    //  leave amps out

    // VZERO
    SetVZEROData(fESDOld->GetVZEROData());

    SetVertex(fESDOld->GetVertex());

    SetPrimaryVertex(fESDOld->GetPrimaryVertex());

    SetMultiplicity(fESDOld->GetMultiplicity());
    
    for(int i = 0;i<fESDOld->GetNumberOfTracks();i++){
      AddTrack(fESDOld->GetTrack(i));
    }

    for(int i = 0;i<fESDOld->GetNumberOfMuonTracks();i++){
      AddMuonTrack(fESDOld->GetMuonTrack(i));
    }

    for(int i = 0;i<fESDOld->GetNumberOfPmdTracks();i++){
      AddPmdTrack(fESDOld->GetPmdTrack(i));
    }

    for(int i = 0;i<fESDOld->GetNumberOfTrdTracks();i++){
      AddTrdTrack(fESDOld->GetTrdTrack(i));
    }

    for(int i = 0;i<fESDOld->GetNumberOfV0s();i++){
      AddV0(fESDOld->GetV0(i));
    }

    for(int i = 0;i<fESDOld->GetNumberOfCascades();i++){
      AddCascade(fESDOld->GetCascade(i));
    }

    for(int i = 0;i<fESDOld->GetNumberOfKinks();i++){
      AddKink(fESDOld->GetKink(i));
    }


    for(int i = 0;i<fESDOld->GetNumberOfCaloClusters();i++){
      AddCaloCluster(fESDOld->GetCaloCluster(i));
    }
  }// if fesdold
}



