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
// AliMUONData classes
// Class containing MUON data: hits, digits, rawclusters, globaltrigger, localtrigger, etc ..
// The classe makes the lik between the MUON data lists and the event trees from loaders
// Gines Martinez, Subatech,  September 2003
//

//Root includes
#include "TNamed.h"
//AliRoot include
#include "AliRun.h"
#include "AliMC.h" 
#include "AliLoader.h" 
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONHit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONTriggerTrack.h"
ClassImp(AliMUONData)
 
//_____________________________________________________________________________
  AliMUONData::AliMUONData():
    TNamed(),
    fLoader(0x0),
    fHits(0x0),
    fDigits(0x0),
    fSDigits(0x0),
    fRawClusters(0x0),
    fGlobalTrigger(0x0),
    fLocalTrigger(0x0),
    fRecTracks(0x0),
    fRecTriggerTracks(0x0),
    fNhits(0),
    fNdigits(0x0),
    fNSdigits(0x0),
    fNrawclusters(0x0),
    fNglobaltrigger(0),
    fNlocaltrigger(0),
    fNrectracks(0),
    fNrectriggertracks(0),
    fSplitLevel(0)
{
  // Default constructor
}
//_____________________________________________________________________________
AliMUONData::AliMUONData(AliLoader * loader, const char* name, const char* title):
  TNamed(name,title),
    fLoader(loader),
    fHits(0x0),
    fDigits(0x0),
    fSDigits(0x0),
    fRawClusters(0x0),
    fGlobalTrigger(0x0),
    fLocalTrigger(0x0),
    fRecTracks(0x0),
    fRecTriggerTracks(0x0),
    fNhits(0),
    fNdigits(0x0),
    fNSdigits(0x0),
    fNrawclusters(0x0),
    fNglobaltrigger(0),
    fNlocaltrigger(0),
    fNrectracks(0),
    fNrectriggertracks(0),
    fSplitLevel(0)
{
  // Constructor for AliMUONData

//   fHits          = new TClonesArray("AliMUONHit",1000);
//   fNhits         = 0;
//   fDigits        = new TObjArray(AliMUONConstants::NCh());
//   fNdigits       = new Int_t[AliMUONConstants::NCh()];
//   for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
//     fDigits->AddAt(new TClonesArray("AliMUONDigit",10000),iDetectionPlane); 
//     fNdigits[iDetectionPlane]=0;
//   }
//   fRawClusters   = new TObjArray(AliMUONConstants::NTrackingCh());
//   fNrawclusters  = new Int_t[AliMUONConstants::NTrackingCh()];
//   for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NTrackingCh();iDetectionPlane++) {
//     fRawClusters->AddAt(new TClonesArray("AliMUONRawCluster",10000),iDetectionPlane); 
//     fNrawclusters[iDetectionPlane]=0;
//   }
//   fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1);    
//   fNglobaltrigger =0;
//   fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);   
//   fNlocaltrigger = 0;
//   fRecTracks     = new TClonesArray("AliMUONTrack", 100);
//   fNrectracks    = 0; // really needed or GetEntriesFast sufficient ????


}

//_____________________________________________________________________________
AliMUONData::AliMUONData(const AliMUONData& rMUONData):TNamed(rMUONData)
{
// Protected copy constructor

  Fatal("AliMUONData", "Not implemented.");
}

//_____________________________________________________________________________
AliMUONData::~AliMUONData()
{
  // Destructor for AliMUONData
  if (fHits) {
    fHits->Delete();
    delete fHits;
  }
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
  }
  if (fSDigits) {
    fSDigits->Delete();
    delete fSDigits;
  }
  if (fRawClusters) {
    fRawClusters->Delete();
    delete fRawClusters;
  }
  if (fGlobalTrigger){
    fGlobalTrigger->Delete();
    delete fGlobalTrigger;
  }  
  if (fLocalTrigger){
    fLocalTrigger->Delete();
    delete fLocalTrigger;
  }
  if (fRecTracks){
    fRecTracks->Delete();
    delete fRecTracks;
  }
  if (fRecTriggerTracks){
    fRecTriggerTracks->Delete();
    delete fRecTriggerTracks;
  }
  //detructor 
}

//_____________________________________________________________________________
AliMUONData& AliMUONData::operator=(const AliMUONData& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  Fatal("operator=", "Not implemented.");
    
  return *this;  
}    
          

//_____________________________________________________________________________
void AliMUONData::AddDigit(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits)
{
  //
  // Add a MUON digit to the list of Digits of the detection plane id
  //
  TClonesArray &ldigits = * Digits(id) ; 
  new(ldigits[fNdigits[id]++]) AliMUONDigit(tracks,charges,digits);
}
//_____________________________________________________________________________
void AliMUONData::AddDigit(Int_t id, const AliMUONDigit& digit)
{
  //
  // Add a MUON digit to the list of Digits of the detection plane id
  //
  TClonesArray &ldigits = * Digits(id) ; 
  new(ldigits[fNdigits[id]++]) AliMUONDigit(digit);
}
//_____________________________________________________________________________
void AliMUONData::AddSDigit(Int_t id, Int_t *tracks, Int_t *charges, Int_t *Sdigits)
{
  //
  // Add a MUON Sdigit to the list of SDigits of the detection plane id
  //
  TClonesArray &lSdigits = * SDigits(id) ; 
  new(lSdigits[fNSdigits[id]++]) AliMUONDigit(tracks,charges,Sdigits);
}
//_____________________________________________________________________________
void AliMUONData::AddSDigit(Int_t id, const AliMUONDigit& Sdigit)
{
  //
  // Add a MUON Sdigit to the list of SDigits of the detection plane id
  //
  TClonesArray &lSdigits = * SDigits(id) ; 
  new(lSdigits[fNSdigits[id]++]) AliMUONDigit(Sdigit);
}
//_____________________________________________________________________________
void AliMUONData::AddGlobalTrigger(Int_t *singlePlus, Int_t *singleMinus,
				   Int_t *singleUndef,
				   Int_t *pairUnlike, Int_t *pairLike)
{
  // add a MUON Global Trigger to the list (only one GlobalTrigger per event !)
  TClonesArray &globalTrigger = *fGlobalTrigger;
  new(globalTrigger[fNglobaltrigger++]) 
    AliMUONGlobalTrigger(singlePlus, singleMinus,  singleUndef, pairUnlike, pairLike);
}
//_____________________________________________________________________________
void AliMUONData::AddGlobalTrigger(const AliMUONGlobalTrigger& trigger )
{
  // add a MUON Global Trigger to the list (only one GlobalTrigger per event !)
  TClonesArray &globalTrigger = *fGlobalTrigger;
  new(globalTrigger[fNglobaltrigger++]) AliMUONGlobalTrigger(trigger);
}
//_____________________________________________________________________________
void AliMUONData::AddHit(Int_t fIshunt, Int_t track, Int_t iChamber, 
			 Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
			 Float_t tof, Float_t momentum, Float_t theta, 
			 Float_t phi, Float_t length, Float_t destep)
{
  // Add new hit to the hit list
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliMUONHit(fIshunt, track, iChamber, 
				  idpart, X, Y, Z, 
				  tof, momentum, theta, 
				  phi, length, destep);
}
//____________________________________________________________________________
void AliMUONData::AddHit(Int_t fIshunt, Int_t track, Int_t iChamber, 
			 Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
			 Float_t tof, Float_t momentum, Float_t theta, 
			 Float_t phi, Float_t length, Float_t destep,
			 Float_t Xref,Float_t Yref,Float_t Zref)
{
 // Add new hit to the hit list
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliMUONHit(fIshunt, track, iChamber, 
				  idpart, X, Y, Z, 
				  tof, momentum, theta, 
				  phi, length, destep,
				  Xref,Yref,Zref);
}
//____________________________________________________________________________
void AliMUONData::AddHit(const AliMUONHit& hit)
{
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliMUONHit(hit);
}
//____________________________________________________________________________
void AliMUONData::AddLocalTrigger(Int_t *localtr)
{
  // add a MUON Local Trigger to the list
  TClonesArray &localTrigger = *fLocalTrigger;
  new(localTrigger[fNlocaltrigger++]) AliMUONLocalTrigger(localtr);
}
//____________________________________________________________________________
void AliMUONData::AddLocalTrigger(const  AliMUONLocalTrigger& trigger)
{
  // add a MUON Local Trigger to the list
  TClonesArray &localTrigger = *fLocalTrigger;
  new(localTrigger[fNlocaltrigger++]) AliMUONLocalTrigger(trigger);
}
//_____________________________________________________________________________
void AliMUONData::AddRawCluster(Int_t id, const AliMUONRawCluster& c)
{
  //
  // Add a MUON rawcluster to the list in the detection plane id
  //
  TClonesArray &lrawcl = *((TClonesArray*) fRawClusters->At(id));
  new(lrawcl[fNrawclusters[id]++]) AliMUONRawCluster(c);
}
//_____________________________________________________________________________
void AliMUONData::AddRecTrack(const AliMUONTrack& track)
{
  //
  // Add a MUON rectrack
  //
  TClonesArray &lrectracks = *fRecTracks;
  new(lrectracks[fNrectracks++]) AliMUONTrack(track);
  //  printf("TTTTTT %d ,\n",((AliMUONTrack*)fRecTracks->At(fNrectracks-1))->GetNTrackHits());
}
//_____________________________________________________________________________
void AliMUONData::AddRecTriggerTrack(const AliMUONTriggerTrack& triggertrack)
{
  //
  // Add a MUON triggerrectrack
  //
  TClonesArray &lrectriggertracks = *fRecTriggerTracks;  
  new(lrectriggertracks[fNrectriggertracks++]) AliMUONTriggerTrack(triggertrack);
  //  printf("TTTTTT %d ,\n",((AliMUONTrack*)fRecTracks->At(fNrectracks-1))->GetNTrackHits());
}

//____________________________________________________________________________
TClonesArray*  AliMUONData::Digits(Int_t DetectionPlane) 
{
  //Getting List of Digits
  if (fDigits)
    return ( (TClonesArray*) fDigits->At(DetectionPlane) );
  else
    return NULL;
}
//____________________________________________________________________________
TClonesArray*  AliMUONData::SDigits(Int_t DetectionPlane) 
{
  //Getting List of SDigits
  if (fSDigits)
    return ( (TClonesArray*) fSDigits->At(DetectionPlane) );
  else
    return NULL;
}
//____________________________________________________________________________
Bool_t   AliMUONData::IsRawClusterBranchesInTree()
{
  // Checking if there are RawCluster Branches In TreeR
  if (TreeR()==0x0) {
    Error("TreeR","No treeR in memory");
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
Bool_t   AliMUONData::IsTriggerBranchesInTree()
{
  // Checking if there are Trigger Branches In TreeR
 if (TreeR()==0x0) {
    Error("TreeR","No treeR in memory");
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
Bool_t   AliMUONData::IsTrackBranchesInTree()
{
  // Checking if there are Track Branches In TreeT
  if (TreeT()==0x0) {
    Error("TreeT","No treeT in memory");
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
Bool_t   AliMUONData::IsTriggerTrackBranchesInTree()
{
  // Checking if there are TriggerTrack Branches In TreeT
  if (TreeT()==0x0) {
    Error("TreeT","No treeT in memory");
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
void AliMUONData::Fill(Option_t* option)
{
  // Method to fill the trees
  const char *cH   = strstr(option,"H");
  const char *cD   = strstr(option,"D");   // Digits branches in TreeD
  const char *cS   = strstr(option,"S");   // SDigits branches in TreeS
  const char *cRC  = strstr(option,"RC");  // RawCluster branches in TreeR
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger branches in TreeR
  const char *cRT  = strstr(option,"RT");  // Reconstructed Track in TreeT
  const char *cRL = strstr(option,"RL");  // Reconstructed Trigger Track in TreeT

  //const char *cRP  = strstr(option,"RP");  // Reconstructed Particle in TreeP
  
  char branchname[30];
  TBranch * branch = 0x0;

  //
  // Filling TreeH
  if ( TreeH() && cH ) {
    TreeH()->Fill();
  }  
  //
  // Filling TreeD
  if ( TreeD() && cD) {
    TreeD()->Fill();
  }
  // Filling TreeS
  if ( TreeS() && cS) {
    TreeS()->Fill();
  }

  //
  // filling rawclusters
  if ( TreeR()  && cRC ) {
    if ( IsTriggerBranchesInTree() ) {
      // Branch per branch filling
      for (int i=0; i<AliMUONConstants::NTrackingCh(); i++) {
	sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	branch = TreeR()->GetBranch(branchname);
	branch->Fill();
      }
    }
    else  TreeR()->Fill();
  }
  
 //
  // filling trigger 
  if ( TreeR()  && cGLT) {
    if (IsRawClusterBranchesInTree()) {
      // Branch per branch filling
      sprintf(branchname,"%sLocalTrigger",GetName());
      branch = TreeR()->GetBranch(branchname); 
      branch->Fill();
      sprintf(branchname,"%sGlobalTrigger",GetName());
      branch = TreeR()->GetBranch(branchname);
      branch->Fill();
    }
    else  TreeR()->Fill();
  }
  //
  // filling tracks
  if ( TreeT() && cRT ) {
    if (IsTriggerTrackBranchesInTree()) {
	sprintf(branchname,"%sTrack",GetName());  
	branch = TreeT()->GetBranch(branchname);
	branch->Fill();
    }
    else  TreeT()->Fill();
  }
  // filling trigger tracks
  if ( TreeT() && cRL ) {
    if (IsTrackBranchesInTree()) {
	sprintf(branchname,"%sTriggerTrack",GetName());  
	branch = TreeT()->GetBranch(branchname);
	branch->Fill();
    }    
    else TreeT()->Fill();
  }
//   if ( TreeT() && cRL ) {
//     sprintf(branchname,"%sTrackTrig",GetName());  
//     TreeT()->Fill();
//   }
}
//_____________________________________________________________________________
void AliMUONData::MakeBranch(Option_t* option)
{
  //
  // Create Tree branches for the MUON.
  //
  const Int_t kBufferSize = 4000;
  char branchname[30];
  

  const char *cH   = strstr(option,"H");
  const char *cD   = strstr(option,"D");   // Digits branches in TreeD
  const char *cS   = strstr(option,"S");   // Digits branches in TreeS
  const char *cRC  = strstr(option,"RC");  // RawCluster branches in TreeR
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger branches in TreeR
  const char *cRT  = strstr(option,"RT");  // Reconstructed Track in TreeT
  const char *cRL  = strstr(option,"RL");  // Reconstructed Trigger Track in TreeT
  //const char *cRP  = strstr(option,"RP");  // Reconstructed Particle in TreeP

  TBranch * branch = 0x0;
  
  // Creating Branches for Hits
  if (TreeH() && cH) {

    if (fHits == 0x0)  {
	fHits = new TClonesArray("AliMUONHit",1000);
// 	if (gAlice->GetMCApp())
// 	  gAlice->GetMCApp()->AddHitList (fHits);
    }
	    
    fNhits = 0;
    sprintf(branchname,"%sHits",GetName());  
    branch = TreeH()->GetBranch(branchname);
    if (branch) {  
      Info("MakeBranch","Branch %s is already in tree.",GetName());
      return ;
    }
    branch = TreeH()->Branch(branchname,&fHits,kBufferSize);
    //Info("MakeBranch","Making Branch %s for hits \n",branchname);
  }  
  
  //Creating Branches for Digits
  if (TreeD() && cD ) {
    // one branch for digits per chamber
    if (fDigits  == 0x0) {
      fDigits  = new TObjArray(AliMUONConstants::NCh());
      for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
	TClonesArray * tca = new TClonesArray("AliMUONDigit",10000);
	tca->SetOwner();
	fDigits->AddAt(tca,iDetectionPlane); 
      }
    }
    if (fNdigits == 0x0) {
      fNdigits = new Int_t[AliMUONConstants::NCh()];
      for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
	fNdigits[iDetectionPlane]=0;
      }
    }
    for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
      sprintf(branchname,"%sDigits%d",GetName(),iDetectionPlane+1);
      branch = 0x0;
      branch = TreeD()->GetBranch(branchname);
      if (branch) {  
	Info("MakeBranch","Branch %s is already in tree.",GetName());
	return;
      }
      TClonesArray * digits = Digits(iDetectionPlane); 
      branch = TreeD()->Branch(branchname, &digits, kBufferSize,1);
      //Info("MakeBranch","Making Branch %s for digits in detection plane %d\n",branchname,iDetectionPlane+1);
      }
  }
  
  //Creating Branches for SDigits
  if (TreeS() && cS ) {
    // one branch for Sdigits per chamber
    if (fSDigits  == 0x0) {
      fSDigits  = new TObjArray(AliMUONConstants::NCh());
      for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
	TClonesArray * tca = new TClonesArray("AliMUONDigit",10000);
	tca->SetOwner();
	fSDigits->AddAt(tca,iDetectionPlane); 
      }
    }
    if (fNSdigits == 0x0) {
      fNSdigits = new Int_t[AliMUONConstants::NCh()];
      for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
	fNSdigits[iDetectionPlane]=0;
      }
    }
    for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
      sprintf(branchname,"%sSDigits%d",GetName(),iDetectionPlane+1);
      branch = 0x0;
      branch = TreeS()->GetBranch(branchname);
      if (branch) {  
	Info("MakeBranch","Branch %s is already in tree.",GetName());
	return;
      }
      TClonesArray * Sdigits = SDigits(iDetectionPlane); 
      branch = TreeS()->Branch(branchname, &Sdigits, kBufferSize,1);
      //Info("MakeBranch","Making Branch %s for Sdigits in detection plane %d\n",branchname,iDetectionPlane+1);
      }
  }

  if (TreeR() && cRC ) {
    //  one branch for raw clusters per tracking detection plane
    //        
    Int_t i;
    if (fRawClusters == 0x0) {
      fRawClusters = new TObjArray(AliMUONConstants::NTrackingCh());
      for (Int_t i=0; i<AliMUONConstants::NTrackingCh();i++) {
	TClonesArray * tca = new TClonesArray("AliMUONRawCluster",1000);
	tca->SetOwner();
	fRawClusters->AddAt(tca,i); 
      }
    }

    if (fNrawclusters == 0x0) {
      fNrawclusters= new Int_t[AliMUONConstants::NTrackingCh()];
      for (Int_t i=0; i<AliMUONConstants::NTrackingCh();i++) {
	fNrawclusters[i]=0;
      }
    }
    
    for (i=0; i<AliMUONConstants::NTrackingCh() ;i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);	
      branch = 0x0;
      branch = TreeR()->GetBranch(branchname);
      if (branch) {  
	Info("MakeBranch","Branch %s is already in tree.",GetName());
	return;
      }
      branch = TreeR()->Branch(branchname, &((*fRawClusters)[i]),kBufferSize);
      //Info("MakeBranch","Making Branch %s for rawcluster in detection plane %d\n",branchname,i+1);
    }
  }

  if (TreeR() && cGLT ) {
    //
    // one branch for global trigger
    //
    sprintf(branchname,"%sGlobalTrigger",GetName());
    branch = 0x0;
    
    if (fGlobalTrigger == 0x0) {
      fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger"); 
      fNglobaltrigger = 0;
    }
    branch = TreeR()->GetBranch(branchname);
    if (branch) {  
      Info("MakeBranch","Branch %s is already in tree.",GetName());
      return ;
    }
    branch = TreeR()->Branch(branchname, &fGlobalTrigger, kBufferSize);
    //Info("MakeBranch", "Making Branch %s for Global Trigger\n",branchname);
    
    //
    // one branch for local trigger
    //  
    sprintf(branchname,"%sLocalTrigger",GetName());
    branch = 0x0;
    
    if (fLocalTrigger == 0x0) {
      fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);
      fNlocaltrigger = 0;
    }
    branch = TreeR()->GetBranch(branchname);
    if (branch) {  
      Info("MakeBranch","Branch %s is already in tree.",GetName());
      return;
    }
    branch = TreeR()->Branch(branchname, &fLocalTrigger, kBufferSize);
    //Info("MakeBranch", "Making Branch %s for Global Trigger\n",branchname);  
  }
  
  if (TreeT() && cRT ) {
    if (fRecTracks == 0x0)  fRecTracks = new TClonesArray("AliMUONTrack",100);
    fNrectracks = 0;
    sprintf(branchname,"%sTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) {  
      Info("MakeBranch","Branch %s is already in tree.",GetName());
      return ;
    }
    branch = TreeT()->Branch(branchname,&fRecTracks,kBufferSize);
    //Info("MakeBranch","Making Branch %s for tracks \n",branchname);
  }  
// trigger tracks
  if (TreeT() && cRL ) {
    if (fRecTriggerTracks == 0x0)  fRecTriggerTracks = new TClonesArray("AliMUONTriggerTrack",100);
    fNrectriggertracks = 0;
    sprintf(branchname,"%sTriggerTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) {  
      Info("MakeBranch","Branch %s is already in tree.",GetName());
      return ;
    }
    branch = TreeT()->Branch(branchname,&fRecTriggerTracks,kBufferSize);
    //Info("MakeBranch","Making Branch %s for trigger tracks \n",branchname);
  }  
}
//____________________________________________________________________________
TClonesArray*  AliMUONData::RawClusters(Int_t DetectionPlane)
{
  // Getting Raw Clusters
  if (fRawClusters) 
    return ( (TClonesArray*) fRawClusters->At(DetectionPlane) );
  else
    return NULL;
}
//____________________________________________________________________________
TClonesArray*  AliMUONData::LocalTrigger()
{
  // Getting Local Trigger
  if (fLocalTrigger) 
    return ( (TClonesArray*) fLocalTrigger );
  else
    return NULL;
}
//____________________________________________________________________________
TClonesArray*  AliMUONData::GlobalTrigger()
{
  // Getting Global Trigger
  if (fGlobalTrigger) 
    return ( (TClonesArray*) fGlobalTrigger );
  else
    return NULL;
}
//____________________________________________________________________________
void AliMUONData::ResetDigits()
{
    //
    // Reset number of digits and the digits array for this detector
    //
    if (fDigits == 0x0) return;
    for ( int i=0;i<AliMUONConstants::NCh();i++ ) {
      if ((*fDigits)[i])    ((TClonesArray*)fDigits->At(i))->Clear();
      if (fNdigits)  fNdigits[i]=0;
    }
}
//____________________________________________________________________________
void AliMUONData::ResetSDigits()
{
    //
    // Reset number of Sdigits and the Sdigits array for this detector
    //
    if (fSDigits == 0x0) return;
    for ( int i=0;i<AliMUONConstants::NCh();i++ ) {
      if ((*fSDigits)[i])    ((TClonesArray*)fSDigits->At(i))->Clear();
      if (fNSdigits)  fNSdigits[i]=0;
    }
}
//______________________________________________________________________________
void AliMUONData::ResetHits()
{
  // Reset number of clusters and the cluster array for this detector
  fNhits   = 0;
  if (fHits) fHits->Clear();
}
//_______________________________________________________________________________
void AliMUONData::ResetRawClusters()
{
    // Reset number of raw clusters and the raw clust array for this detector
    //
  for ( int i=0;i<AliMUONConstants::NTrackingCh();i++ ) {
    if ((*fRawClusters)[i])    ((TClonesArray*)fRawClusters->At(i))->Clear();
    if (fNrawclusters)  fNrawclusters[i]=0;
  }
}
//_______________________________________________________________________________
void AliMUONData::ResetTrigger()
{
  //  Reset Local and Global Trigger 
  fNglobaltrigger = 0;
  if (fGlobalTrigger) fGlobalTrigger->Clear();
  fNlocaltrigger = 0;
  if (fLocalTrigger) fLocalTrigger->Clear();
}
//____________________________________________________________________________
void AliMUONData::ResetRecTracks()
{
  // Reset tracks information
  fNrectracks = 0;
  if (fRecTracks) fRecTracks->Clear();
}
//____________________________________________________________________________
void AliMUONData::ResetRecTriggerTracks()
{
  // Reset tracks information
  fNrectriggertracks = 0;
  if (fRecTriggerTracks) fRecTriggerTracks->Clear();
}
//_____________________________________________________________________________
void AliMUONData::SetTreeAddress(Option_t* option)
{
  //Setting Addresses to the events trees
  const char *cH   = strstr(option,"H");
  const char *cD   = strstr(option,"D");   // Digits branches in TreeD
  const char *cS   = strstr(option,"S");   // SDigits branches in TreeS
  const char *cRC  = strstr(option,"RC");  // RawCluster branches in TreeR
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger branches in TreeR
  const char *cRT  = strstr(option,"RT");  // Reconstructed Track in TreeT
  const char *cRL  = strstr(option,"RL");  // Reconstructed Trigger Track in TreeT
  //const char *cRP  = strstr(option,"RP");  // Reconstructed Particle in TreeP
  
  // Set branch address for the Hits, Digits, RawClusters, GlobalTrigger and LocalTrigger Tree.
  char branchname[30];
  TBranch * branch = 0x0;

  //
  // Branch address for hit tree
  if ( TreeH() && cH ) {
      if (fHits == 0x0) {
	fHits     = new TClonesArray("AliMUONHit",1000);
	//	if (gAlice->GetMCApp())
	//  gAlice->GetMCApp()->AddHitList (fHits);  Moved to AliMUON
    }
    fNhits =0;
  } 
  if (TreeH() && fHits && cH) {
    sprintf(branchname,"%sHits",GetName());  
    branch = TreeH()->GetBranch(branchname);
    if (branch) {
      //      Info("SetTreeAddress","(%s) Setting for Hits",GetName());
      branch->SetAddress(&fHits);
    }
    else { //can be invoked before branch creation
      Warning("SetTreeAddress","(%s) Failed for Hits. Can not find branch in tree.",GetName());
    }
  }
  
  //
  // Branch address for digit tree
  if ( TreeD() && cD) {
    if (fDigits == 0x0) { 
      fDigits = new TObjArray(AliMUONConstants::NCh());
      fNdigits= new Int_t[AliMUONConstants::NCh()];
      for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
	fDigits->AddAt(new TClonesArray("AliMUONDigit",10000),i); 
	fNdigits[i]=0;
      }
    }
  }

  if (TreeD() && fDigits && cD) {
    for (int i=0; i<AliMUONConstants::NCh(); i++) {
      sprintf(branchname,"%sDigits%d",GetName(),i+1);
      if (fDigits) {
	branch = TreeD()->GetBranch(branchname);
	TClonesArray * digits = Digits(i);
	if (branch) {
	  branch->SetAddress( &digits );
	}
	else Warning("SetTreeAddress","(%s) Failed for Digits Detection plane %d. Can not find branch in tree.",GetName(),i);
      }
    }
  }
  //
  // Branch address for Sdigit tree
  if ( TreeS() && cS) {
    if (fSDigits == 0x0) { 
      fSDigits = new TObjArray(AliMUONConstants::NCh());
      fNSdigits= new Int_t[AliMUONConstants::NCh()];
      for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
	fSDigits->AddAt(new TClonesArray("AliMUONDigit",10000),i); 
	fNSdigits[i]=0;
      }
    }
  }

  if (TreeS() && fSDigits && cS) {
    for (int i=0; i<AliMUONConstants::NCh(); i++) {
      sprintf(branchname,"%sSDigits%d",GetName(),i+1);
      if (fSDigits) {
	branch = TreeS()->GetBranch(branchname);
	TClonesArray * Sdigits = SDigits(i);
	if (branch) branch->SetAddress( &Sdigits );
	else Warning("SetTreeAddress","(%s) Failed for SDigits Detection plane %d. Can not find branch in tree.",GetName(),i);
      }
    }
  }
  
  //
  // Branch address for rawclusters, globaltrigger and local trigger tree
  if (TreeR() ) {
    if (fRawClusters == 0x0 && cRC) {
      fRawClusters = new TObjArray(AliMUONConstants::NTrackingCh());
      fNrawclusters= new Int_t[AliMUONConstants::NTrackingCh()];
      for (Int_t i=0; i<AliMUONConstants::NTrackingCh();i++) {
	fRawClusters->AddAt(new TClonesArray("AliMUONRawCluster",10000),i); 
	fNrawclusters[i]=0;
      }
    }
    if (fLocalTrigger == 0x0 && cGLT) {
      fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);
    }
    if (fGlobalTrigger== 0x0 && cGLT) {
        fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1); 
    }

  }
  if ( TreeR()  && fRawClusters && cRC) {
    for (int i=0; i<AliMUONConstants::NTrackingCh(); i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
      if (fRawClusters) {
	branch = TreeR()->GetBranch(branchname);
	if (branch) branch->SetAddress( &((*fRawClusters)[i]) );
	else Warning("SetTreeAddress","(%s) Failed for RawClusters Detection plane %d. Can not find branch in tree.",GetName(),i);
      }
    }
  }
  if ( TreeR()  && fLocalTrigger && cGLT) {
    sprintf(branchname,"%sLocalTrigger",GetName());
    branch = TreeR()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fLocalTrigger);
    else Warning("SetTreeAddress","(%s) Failed for LocalTrigger. Can not find branch in tree.",GetName());
  }
  if ( TreeR() && fGlobalTrigger && cGLT) {
    sprintf(branchname,"%sGlobalTrigger",GetName());
    branch = TreeR()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fGlobalTrigger);
    else Warning("SetTreeAddress","(%s) Failed for LocalTrigger. Can not find branch in tree.",GetName());
  }

  if ( TreeT() ) {
    if (fRecTracks == 0x0 && cRT)  {
      fRecTracks  = new TClonesArray("AliMUONTrack",100);
    }

  }
  if ( TreeT() && fRecTracks && cRT ) {
    sprintf(branchname,"%sTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fRecTracks);
    else Warning("SetTreeAddress","(%s) Failed for Tracks. Can not find branch in tree.",GetName());
  }
// trigger tracks
  if ( TreeT() ) {
    if (fRecTriggerTracks == 0x0 && cRL)  {
      fRecTriggerTracks  = new TClonesArray("AliMUONTriggerTrack",100);
    }

  }
  if ( TreeT() && fRecTriggerTracks && cRL ) {
    sprintf(branchname,"%sTriggerTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fRecTriggerTracks);
    else Warning("SetTreeAddress","(%s) Failed for Trigger Tracks. Can not find branch in tree.",GetName());
  }


}
//_____________________________________________________________________________
